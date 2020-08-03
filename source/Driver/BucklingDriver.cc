#include "BucklingDriver.hh"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <boost/filesystem.hpp>

#include "Driver/SolveSteps/StdSolveStep.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/SimState.hh"
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/Domain.hh"
#include "MatVec/SBM_Matrix.hh"
#include "OLAS/solver/BaseEigenSolver.hh"
#include "OLAS/algsys/AlgebraicSys.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Optimization.hh"
#include "PDE/StdPDE.hh"

#ifdef USE_ARPACK
#include "OLAS/external/arpack/ArpackEigenSolver.hh"
#endif

using std::cout;
using std::setw;
using std::string;

/**
 * This class solves a buckling problem, which is a general eigenvalue problem:
 * \f$(K - \lambda G) * v = 0\f$
 * with
 * K              stiffness matrix (positive definite)
 * G              geometric stiffness matrix (indefinite)
 * \f$\lambda\f$  eigenvalue
 * v              eigenvector
 * For this problem we have \f$v^T * d * G * v = 1\f$.
 * However, one often uses a reformulated problem, which is obtained by
 * multiplying the equation with \f$-1/\lambda\f$:
 * \f$(G - 1/\lambda * K) * v = 0\f$
 * Now for the eigenvectors \f$v^T * K * v = 1\f$.
 * Unfortunately, in the reformulated problem load factors will be sorted
 * in decreasing order, such that we have to search for the biggest
 * eigenvalue in optimization, which is hard for our eigenvalue problem solvers.
 */

namespace CoupledField {

DEFINE_LOG(buckD, "bucklingDriver")

// forward declaration
class BaseVector;

// ***************
//   Constructor
// ***************
BucklingDriver::BucklingDriver(UInt sequenceStep,
                              bool isPartOfSequence,
                              shared_ptr<SimState> state,
                              Domain *domain,
                              PtrParamNode paramNode,
                              PtrParamNode infoNode)
  : SingleDriver( sequenceStep, isPartOfSequence, state, domain, paramNode, infoNode )
{
  // set analysistype
  analysis_ = BasePDE::BUCKLING;
  // input parameter
  inputMethod_ = NONE;
  numMode_ = 0;
  valueShift_ = 0.0;
  calcModes_ = false;
  minVal_ = 0.0;
  maxVal_ = 0.0;
  modeNormalization_ = BaseEigenSolver::NONE;
  isStoredSymmetric_ = false;

  numEigenValues_ = -1;
  eigenValues = new Vector<Complex>();
  errBounds = new Vector<Complex>();

  loadFactors_ = new Vector<Complex>();

  //specifying parameter node
  param_ = param_->Get("buckling");
  info_ = info_->Get("buckling");

  // some variables needed later
  writeAllSteps_ = false;
}

BucklingDriver::~BucklingDriver() {

}

void BucklingDriver::Init(bool restart) {
  // read required parameters from parameter node
  param_->GetValue("inverse", isInverseProblem_, ParamNode::INSERT);
  param_->GetValue("numModes", numMode_, ParamNode::INSERT);
  param_->GetValue("valueShift", valueShift_, ParamNode::INSERT);
  param_->GetValue("minVal", minVal_, ParamNode::INSERT);
  param_->GetValue("maxVal", maxVal_, ParamNode::INSERT);
  param_->GetValue("allowPostProc", writeAllSteps_, ParamNode::PASS);
  calcModes_  = param_->Has("calcModes");

  if (calcModes_) {
    // determine type of mode normalization, and transform into ENUM
    std::string normString = "solver";
    param_->Get("calcModes")->GetValue("normalization", normString, ParamNode::PASS);
    if (normString == "solver") {
      modeNormalization_ = BaseEigenSolver::NONE;
    }
    else if (normString == "max") {
      modeNormalization_ = BaseEigenSolver::MAX;
    }
    else if (normString == "norm") {
      modeNormalization_ = BaseEigenSolver::NORM;
    }
    else {
      EXCEPTION("Specified mode normalization '" + normString + "' not implemented");
    }
  }
  // set definition of PDE relevant matrices
  InitializePDEs();

}

void BucklingDriver::SolveProblem() {
  // define result handler
  ResultHandler *resHandler = domain_->GetResultHandler();

  // we write the results only if we are not optimization. Optimization writes the results by itself via calling StoreResults().
  // But we call PrintResult() to do the output to the info.xml even in case of optimization.

  UInt numSteps = 1;
  // see comments in StaticDriver::SolveProblem() for the interplay with optimization
  if (!domain->GetOptimization())
    resHandler->BeginMultiSequenceStep(sequenceStep_, analysis_, numSteps); // optimization does it by itself

  if (writeAllSteps_ || isPartOfSequence_)
    simState_->BeginMultiSequenceStep(sequenceStep_, analysis_);

  // trigger calculation
  ptPDE_->WriteGeneralPDEdefines();
  StdSolveStep *sstep = dynamic_cast<StdSolveStep*>(ptPDE_->GetSolveStep());

  // set the mode normalization
  sstep->GetAlgSys()->GetEigenSolver()->SetModeNormalization(modeNormalization_);

  // initialize algebraic system
  sstep->GetAlgSys()->InitSol();
  sstep->GetAlgSys()->InitMatrix();
  sstep->GetAssemble()->AssembleMatrices();
  sstep->GetAlgSys()->ExportLinSys(true, false, false); // if asked, export matrices

  // define input method, which is needed in CalcValues()
  // inputMethod_ is 1 for minVal & maxVal, 2 for numModes & shiftMode
  if (minVal_ < maxVal_) {
    inputMethod_ = 1;
  }
  else if (minVal_ > maxVal_) {
    EXCEPTION("minVal is bigger than maxVal. Please correct XML input file.");
  }
  else if (numMode_ > 0) {
    inputMethod_ = 2;
  }
  else {
    // this case should not be possible in the XML schema
    EXCEPTION("Input method cannot be determined.")
  }

  // actually solve problem
  CalcValues();
  SortModes(false);
  if(calcModes_)
    StoreModes();

  PrintResult();

  // in optimization we write the results via StoreResults() because
  // we don't necessarily write every forward step.
  if(!domain->GetOptimization())
  {
    StoreResults(1, -1.0);
    handler_->FinishMultiSequenceStep();

    if (!isPartOfSequence_)
      handler_->Finalize(); // to be called only once in a HDF5 lifetime!
  }

  if (writeAllSteps_ || isPartOfSequence_)
    simState_->FinishMultiSequenceStep(true);
}

void BucklingDriver::CalcValues() {
  BaseSolveStep *step = ptPDE_->GetSolveStep();
  StdSolveStep *sstep = dynamic_cast<StdSolveStep*>(step);

  // get active solver
  BaseEigenSolver *solver = sstep->GetAlgSys()->GetEigenSolver();

  // get type of solver
  solverType_ = solver->GetEigenSolverName();

  // get needed matrices
  SBM_Matrix *stiffMat = sstep->GetAlgSys()->GetMatrix(STIFFNESS);
  SBM_Matrix *geoStiffMat = sstep->GetAlgSys()->GetMatrix(GEOMETRIC_STIFFNESS);

  // check matrices dimensions
  if (geoStiffMat->GetNumCols() > 1) {
    EXCEPTION("only implemented for SBM matrices with a single block")
  }
  assert(stiffMat->GetNumCols() == stiffMat->GetNumRows());
  assert(geoStiffMat->GetNumCols() == geoStiffMat->GetNumRows());

  //check storage type
  bool isReal = false;
  solver->CheckMatrix(isReal, isStoredSymmetric_, *(stiffMat->GetPointer(0, 0)));

  // unfortunately the implemented eigenSolvers have serious usage restrictions.
  // currently only FEAST with non-symmetric matrices gives relabel solutions
  if (solverType_ == BaseEigenSolver::FEAST && isStoredSymmetric_ == true) {
    EXCEPTION("Wrong matrix storage type. Please define an non-symmetric storage type.");
  }
  else if (solverType_ == BaseEigenSolver::FEAST && inputMethod_ != 1) {
    EXCEPTION("FEAST eigenSolver does currently not support input as numModes & shiftMode.");
  }
  else if (solverType_ == BaseEigenSolver::ARPACK && inputMethod_ != 2) {
    EXCEPTION("ARPACK eigenSolver does currently not support input as minValue & maxValue.");
  }
  else if (solverType_ == BaseEigenSolver::PHIST) {
    EXCEPTION("PHIST is currently not supported.");
  }

/* TODO FEAST
 * - add inputMethod_ 2, as numModes / shiftMode
 */
/* TODO ARPACK
 * - add inputMethod_ 1, as min/maxValue
 * - add CalcEigenValues()
 */
/* TODO PHIST
 * - add non-symmetric Mode (?)
 */

  if (isStoredSymmetric_ || solverType_ == BaseEigenSolver::ARPACK) {
    eigenValues = new Vector<Double>();
    errBounds = new Vector<Double>();
    loadFactors_ = new Vector<Double>();
  } // else vectors are complex, see constructor

  // solve \f$(G - 1/\lambda * K) * v = 0\f$
  // unfortunately the implemented eigenSolvers do not conform to a single scheme
  if (inputMethod_ == 1) { // inputMethod_ = minVal & maxVal
    assert(solverType_ ==  BaseEigenSolver::FEAST);
    if (isInverseProblem_)
      solver->Setup(*(geoStiffMat->GetPointer(0, 0)), *(stiffMat->GetPointer(0, 0)), false);
    else
      solver->Setup(*(stiffMat->GetPointer(0, 0)), *(geoStiffMat->GetPointer(0, 0)), false);
    solver->CalcEigenValues(*eigenValues, *errBounds, minVal_, maxVal_);
  }
  else if (inputMethod_ == 2) {// inputMethod_ = numMode & valueShift
    assert(solverType_ == BaseEigenSolver::ARPACK);
    isStoredSymmetric_ = true;
    // this is necessary for the original problem, but ill lead to wrong results
    // for the reformulated problem due to wrong ordering of the eigenvalues
    if (valueShift_ == 0.0) {
      valueShift_ = 0.1;
      info_->Get(ParamNode::HEADER)->SetWarning("valueShift = 0 should not be used for buckling. Changed to 0.1.");
    }
    if (isInverseProblem_)
      solver->Setup(*(geoStiffMat->GetPointer(0, 0)), *(stiffMat->GetPointer(0, 0)), numMode_, 1.0/valueShift_, true, false);
    else
    {
#ifdef USE_ARPACK
      dynamic_cast<ArpackEigenSolver*>(solver)->SetComputeMode(ArpackMatInterface::BUCKLING);
#endif
      solver->Setup(*(stiffMat->GetPointer(0, 0)), *(geoStiffMat->GetPointer(0, 0)), numMode_, valueShift_, true, false);
    }
    solver->CalcEigenValues(*eigenValues, *errBounds, numMode_, valueShift_);
  }
  else {
    EXCEPTION("input method not known")
  }

  numEigenValues_ = eigenValues->GetSize();

  // If no eigenvalue at all converged, just leave
  if(numEigenValues_ == 0) {
    EXCEPTION( "Did not find any eigenvalue!" );
  }

  LOG_DBG(buckD) << "CV: eigenvalues = " << eigenValues->ToString();

  loadFactors_->Resize(numEigenValues_);
  for (UInt i = 0; i < numEigenValues_; i++)
  {
    if (isStoredSymmetric_ || solverType_ == BaseEigenSolver::ARPACK) {
      if (isInverseProblem_)
        loadFactors_->SetEntry(i, 1.0/(eigenValues->GetDoubleEntry(i)));
      else
        loadFactors_->SetEntry(i, eigenValues->GetDoubleEntry(i));
    }
    else
    {
      if (isInverseProblem_)
        loadFactors_->SetEntry(i, 1.0/(eigenValues->GetComplexEntry(i)));
      else
        loadFactors_->SetEntry(i, eigenValues->GetComplexEntry(i));
    }
  }

  LOG_TRACE(buckD) << "CV: loadFactors = " << loadFactors_->ToString();

  if(!domain->GetOptimization())
    cout << "\n++ Finished solving eigenvalue Problem." << std::endl;
}

void BucklingDriver::StoreModes() {
  Vector<Double> loadFactorsRealPart = GetRealPartOfLoadFactors();

  StdSolveStep *sstep = dynamic_cast<StdSolveStep*>(ptPDE_->GetSolveStep());

  for(unsigned int ev = 0; ev < numEigenValues_; ev++)
  {
    Double currentLoadFactor = loadFactorsRealPart[modeOrder_[ev]];

    // store the eigen mode in algsys sol_
    sstep->SetActStep(ev);
    sstep->SetActFreq(currentLoadFactor);
    if (isStoredSymmetric_ || solverType_ == BaseEigenSolver::ARPACK)
    {
      sstep->SetActFreq(loadFactors_->GetDoubleEntry(modeOrder_[ev]));
      ptPDE_->GetDomain()->GetMathParser()->SetValue(MathParser::GLOB_HANDLER, "f", loadFactors_->GetDoubleEntry(modeOrder_[ev]));
    }
    else
    {
      sstep->SetActFreq(loadFactors_->GetComplexEntry(modeOrder_[ev]).real());
      ptPDE_->GetDomain()->GetMathParser()->SetValue(MathParser::GLOB_HANDLER, "f", loadFactors_->GetComplexEntry(modeOrder_[ev]).real());
    }
    sstep->GetEigenMode(modeOrder_[ev]); // this stores the eigen mode result in algsys sol_

    LOG_DBG(buckD) << "stored mode " << ev;
  }
}

void BucklingDriver::PrintResult() {
  unsigned int numConverged = eigenValues->GetSize();

  // Issue warning, if number of converged eigenvalues differs from
  // number of requested ones
  if( numMode_ > 0 && numConverged != numMode_ ) {
    WARN( "Only " << numConverged << " eigenvalues of "
        << numMode_ << " converged. To improve convergence, either "
        << "reduce the number of eigenvalues or the tolerance." );
  }

  int fieldwidth = 23; // field width

  // console output (reduced form bloch)
  if(!domain->GetOptimization() || domain->GetOptimization()->GetOptimizerType() == Optimization::EVALUATE_INITIAL_DESIGN)
  {
    cout << "\n";
    cout << " Mode | ";
    cout << setw(fieldwidth) << "Load Factor" << " | ";
    cout << setw(fieldwidth) << "Error Bounds" << " | ";
    if (!isStoredSymmetric_) {
      cout << setw(fieldwidth) << "Imaginary Part" << " | ";
    }
    cout << "\n";
  }

  // also log via info node.
  // new result only on detailed output and for bloch case only for step=0
  PtrParamNode res;
  bool append = progOpts->DoDetailedInfo();
  if(append)
    res = info_->Get("result", ParamNode::APPEND);
  else  {
    // create or take the last
    ParamNodeList l = info_->GetList("result");
    res = l.IsEmpty() ? info_->Get("result") : l.Last();
  }
  if(domain->GetOptimization())
    res->Get("iteration")->SetValue(domain->GetOptimization()->GetCurrentIteration());

  // the problem might be that we have no detailed output, hence the res might already exist.
  // then, numConverged might be smaller than an older optimization iteration and the old stuff > numConverged remains
  res->GetChildren().Resize(0);

  // is not set for pure simulation
  if(analysis_id_.ToString() != "")
    res->Get("id")->SetValue(analysis_id_.ToString());

  assert(modeOrder_.GetSize() > 0); // set by SortModes()

  for(unsigned int i=0; i < numConverged; i++)
  {
    // command line output only when not optimizing
    if(!domain->GetOptimization() || domain->GetOptimization()->GetOptimizerType() == Optimization::EVALUATE_INITIAL_DESIGN)
    {
      cout << setw(5) << i + 1 << " | ";
      if (isStoredSymmetric_ || solverType_ == BaseEigenSolver::ARPACK) {
        cout << setw(fieldwidth) << loadFactors_->GetDoubleEntry(modeOrder_[i]) << " | ";
        cout << setw(fieldwidth) << errBounds->GetDoubleEntry(modeOrder_[i]) << " | " << "\n";
      }
      else {
        Complex ev = loadFactors_->GetComplexEntry(modeOrder_[i]);
        cout << setw(fieldwidth) << ev.real() << " | ";
        cout << setw(fieldwidth) << errBounds->GetComplexEntry(modeOrder_[i]).real() << " | ";
        cout << setw(fieldwidth) << ev.imag() << " | " << "\n";
      }
    }

    // res got cleared above, so we don't overlap with previous in case numConverged is too small
    PtrParamNode mode = res->Get("mode", ParamNode::APPEND);

    mode->Get("nr")->SetValue(i+1); // not the mode but eigenvalue in list
    if (isStoredSymmetric_) {
      mode->Get("loadfactor")->SetValue(loadFactors_->GetDoubleEntry(i),15);
      mode->Get("errorbound")->SetValue(errBounds->GetDoubleEntry(i));
    }
    else {
      mode->Get("loadfactor")->SetValue(loadFactors_->GetComplexEntry(i).real(),15);
      mode->Get("errorbound")->SetValue(errBounds->GetComplexEntry(i).real());
    }
  }
}

unsigned int BucklingDriver::StoreResults(unsigned int stepNum, double step_val)
{
  LOG_DBG2(buckD) << "SR: step=" << stepNum << " val=" << step_val;

  if(calcModes_)
  {
    Vector<Double> loadFactorsRealPart = GetRealPartOfLoadFactors();

    // check if negative Eigenvalues are present
    for (unsigned int ev = 0; ev < numEigenValues_; ev++) {
      if (loadFactorsRealPart[modeOrder_[ev]] < 0) {
        cout << "\n++ WARNING negative proportionality factor will be displayed as positiv in Paraview: " << loadFactorsRealPart[modeOrder_[ev]] << "\n";
      }
    }

    int digs =  boost::lexical_cast<string>((int)numEigenValues_).size() + 2;
    double sig = std::pow((float) 10.0, -digs); // 1e-2 -> 10 ^ -2 -> a compiler complained with simply 10

    for(unsigned int ev = 0; ev < numEigenValues_; ev++)
    {
      // save_value is the "time" value displayed in paraview
      double save_value = -1.0;
      if(domain->GetOptimization() && (domain->GetOptimization()->me->DoHomogenization()
          || domain->GetOptimization()->GetOptimizerType() != Optimization::EVALUATE_INITIAL_DESIGN) )
      {
        // time is step.step_val nr
        save_value = step_val + (ev + 1) * sig; // +1 for one based
        LOG_DBG3(buckD) << "SR: total=" << numEigenValues_ << " digs=" << digs << " sig=" << sig << " count=" << (ev + 1);
      }
      else {
        Vector<Double> loadFactorsRealPart = GetRealPartOfLoadFactors();
        save_value = std::abs(loadFactorsRealPart[modeOrder_[ev]]);
      }

      LOG_DBG2(buckD) << "SR: ev=" << ev << " stepNum=" << stepNum << " save_value=" << save_value;

      // store mode in pde
      if(domain->GetOptimization())
      {
        ErsatzMaterial* em = dynamic_cast<ErsatzMaterial*>(domain->GetOptimization());
        Context* context = domain->GetOptimization()->context;
        Excitation* excite = context->GetExcitation();
        StateSolution* ss = em->GetForwardStates().Get(excite, NULL, ev);
        ss->Write(context->pde); // forward is function NULL
      }

      // write modes to file
      handler_->BeginStep(stepNum, save_value);
      ptPDE_->WriteResultsInFile(stepNum, save_value);
      handler_->FinishStep();

      if(writeAllSteps_ || isPartOfSequence_)
        simState_->WriteStep(stepNum, save_value);

      if (!GetResultHandler()->streamOnly)
        stepNum++;
    }
    if(!domain->GetOptimization())
      cout << "\n++ Finished calculating and storing buckling modes." << std::endl;
  }

  if (!GetResultHandler()->streamOnly)
    return stepNum-1;
  else
    return stepNum;
}

void BucklingDriver::ExportModes() {
  StdSolveStep *sstep = dynamic_cast<StdSolveStep*>(ptPDE_->GetSolveStep());

  // generates a index-array modeOrder_ containing the mode indices sorted by ascending eigenvalue
  SortModes(true);

  // export solution, i.e. modes
  PtrParamNode els = sstep->GetAlgSys()->GetExportLinSysParam();
  if (els && els->Get("solution")->As<bool>()) {
    BaseMatrix::OutputFormat vec_format = BaseMatrix::outputFormat.Parse(els->Get("vecFormat")->As<std::string>());
    std::string base = els->Has("baseName") ? els->Get("baseName")->As<std::string>() : progOpts->GetSimName();
    if(domain->GetDriver()->GetAnalysisId().ToString(true) != "") {
      base += "_" + domain->GetDriver()->GetAnalysisId().ToString(true);
    }
    for (UInt i=0; i< eigenValues->GetSize();i++) {
      Vector<Complex> mode;
      sstep->GetAlgSys()->GetEigenSolver()->GetNormalizedEigenMode(modeOrder_[i], mode);
      mode.Export( base + "_mode_" + lexical_cast<std::string>(modeOrder_[i]+1), vec_format);
//      sstep->GetAlgSys()->GetEigenSolver()->GetNormalizedEigenMode(modeOrder_[i], mode, false);
//      mode.Export( base + "_mode-left_" + lexical_cast<std::string>(modeOrder_[i]+1), vec_format);
    }
  }
}

void BucklingDriver::SortModes(bool inAbs) {
  Vector<Double> loadFactorsRealPart = GetRealPartOfLoadFactors();

  if (inAbs) {
    for (unsigned int i = 0; i < loadFactorsRealPart.GetSize(); i++) {
      loadFactorsRealPart[i] = std::abs(loadFactorsRealPart[i]);
    }
  }

  modeOrder_.Resize(loadFactorsRealPart.GetSize());
  std::size_t n(0);
  // allocate modeOrder_
  std::generate(std::begin(modeOrder_), std::end(modeOrder_), [&] {return n++;});
  // sort it by value
  std::sort(std::begin(modeOrder_), std::end(modeOrder_), [&](int i1, int i2) {
    return loadFactorsRealPart[i1] < loadFactorsRealPart[i2];
  });
}

Vector<Double> BucklingDriver::GetRealPartOfLoadFactors() {
  Vector<Double> loadFactorsRealPart;

  loadFactorsRealPart.Resize(numEigenValues_);
  if (isStoredSymmetric_ || solverType_ == BaseEigenSolver::ARPACK) {
    loadFactorsRealPart = dynamic_cast<Vector<Double>&>(*loadFactors_);
  }
  else
  {
    for (unsigned int mI = 0; mI < numEigenValues_; mI++)
      loadFactorsRealPart[mI] = loadFactors_->GetComplexEntry(mI).real();
  }
  return loadFactorsRealPart;
}

} // end of namespace
