#include "BucklingDriver.hh"

#include <def_use_arpack.hh>

#include <boost/filesystem.hpp>
#include <cmath>
#include <iomanip>
#include <iostream>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/SimState.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Domain/Domain.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "MatVec/SBM_Matrix.hh"
#include "OLAS/algsys/AlgebraicSys.hh"
#include "OLAS/solver/BaseEigenSolver.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Excitation.hh"
#include "Optimization/Optimization.hh"
#include "PDE/StdPDE.hh"

#ifdef USE_ARPACK
#include "OLAS/external/arpack/ArpackEigenSolver.hh"
#endif

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

  numEV_ = -1;
  eigenValues = new Vector<Double>();
  errors = new Vector<Double>();
  loadFactors_ = new Vector<Double>();

  //specifying parameter node
  param_ = param_->Get("buckling");
  info_ = info_->Get("buckling");

  // some variables needed later
  writeAllSteps_ = false;
  isInverseProblem_ = false;

  solver = NULL;
  solverType_ = BaseEigenSolver::NO_EIGENSOLVER;
}

BucklingDriver::~BucklingDriver() {
  delete eigenValues;
  delete errors;
  delete loadFactors_;
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

  if (minVal_ < 0)
    EXCEPTION("minVal has to be non-negative");
  if (maxVal_ < 0 || maxVal_ < minVal_)
    EXCEPTION("maxVal has to be non-negative and larger than minVal");
  if (valueShift_ < 0)
    EXCEPTION("valueShift has to be non-negative");

  if (numMode_ > 0 && valueShift_ == 0.0) {
    valueShift_ = 0.1;
    info_->Get(ParamNode::HEADER)->SetWarning("valueShift = 0 should not be used for buckling. Changed to 0.1.");
  }

  // has to be treated here, else multiple calls to SolveProblem (e.g. during gradient check)
  // would invert valueShift_ all the time
  if (isInverseProblem_) {
    valueShift_ = 1.0/valueShift_;
    // we cannot divide by zero, so slightly shift the value
    if (std::abs(minVal_) < 1e-10) minVal_ = -1e-10;
    if (std::abs(maxVal_) < 1e-10) maxVal_ =  1e-10;
    double tmp = minVal_;
    minVal_ = 1.0/maxVal_;
    maxVal_ = 1.0/tmp;
  }
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
  if (numMode_ > 0) {
    inputMethod_ = 2;
  }
  else if (minVal_ != maxVal_) {
    inputMethod_ = 1;
  } else {
    // this case should not be possible in the XML schema
    EXCEPTION("Input method cannot be determined.")
  }

  // actually solve problem
  SetupSolver();
  CalcValues();

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

void BucklingDriver::SetToStepValue(UInt stepNum, Double stepVal) {
  // ensure that this method is only called if simState has input
  if (! simState_->HasInput()) {
    EXCEPTION("Can only set external time step, if simulation state " << "is read from external file");
  }
  // Set current eigenvalue in the mathParser
  domain_->GetMathParser()->SetValue(MathParser_GLOB_HANDLER, "f", stepVal);
  domain_->GetMathParser()->SetValue(MathParser_GLOB_HANDLER, "step", stepNum);
}


void BucklingDriver::SetupSolver() {
  BaseSolveStep *step = ptPDE_->GetSolveStep();
  StdSolveStep *sstep = dynamic_cast<StdSolveStep*>(step);

  // get active solver
  solver = sstep->GetAlgSys()->GetEigenSolver();

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
  if (solverType_ == BaseEigenSolver::FEAST && inputMethod_ != 1) {
    EXCEPTION("FEAST eigenSolver does currently not support input as numModes & shiftMode.");
  }
  else if (solverType_ == BaseEigenSolver::ARPACK && inputMethod_ != 2) {
    EXCEPTION("ARPACK eigenSolver does currently not support input as minValue & maxValue.");
  }
  else if (solverType_ == BaseEigenSolver::PHIST) {
    EXCEPTION("PHIST is currently not supported.");
  }

  // feast handles non symmetric problems as complex problems
  if (!isStoredSymmetric_ && solverType_ == BaseEigenSolver::FEAST) {
    delete eigenValues;
    delete errors;
    delete loadFactors_;
    eigenValues = new Vector<Complex>();
    errors = new Vector<Complex>();
    loadFactors_ = new Vector<Complex>();
  }

  solver->GetSetupTimer()->Start();

  // setup \f$(G - 1/\lambda * K) * v = 0\f$
  if (isInverseProblem_)
    solver->Setup(*(geoStiffMat->GetPointer(0, 0)), *(stiffMat->GetPointer(0, 0)), false);
  else
#ifdef USE_ARPACK
    // arpack has a special buckling mode
    if (!isInverseProblem_ && solver->GetEigenSolverName() == BaseEigenSolver::ARPACK)
      // this is necessary for the original problem, but might lead to wrong results
      // for the reformulated problem due to wrong ordering of the eigenvalues
      dynamic_cast<ArpackEigenSolver*>(solver)->Setup(*(stiffMat->GetPointer(0, 0)), *(geoStiffMat->GetPointer(0, 0)), ArpackMatInterface::ComputeMode::BUCKLING, false);
    else
#endif
      solver->Setup(*(stiffMat->GetPointer(0, 0)), *(geoStiffMat->GetPointer(0, 0)), false);

  solver->GetSetupTimer()->Stop();
}

void BucklingDriver::CalcValues(unsigned int recursionCount) {
  solver->GetSolveTimer()->Start();
  if (inputMethod_ == 1) { // inputMethod_ = minVal & maxVal
    assert(solverType_ ==  BaseEigenSolver::FEAST);
    solver->CalcEigenValues(*eigenValues, *errors, minVal_, maxVal_);
  }
  else if (inputMethod_ == 2) {// inputMethod_ = numMode & valueShift
    assert(solverType_ == BaseEigenSolver::ARPACK);
    solver->CalcEigenValues(*eigenValues, *errors, numMode_, valueShift_);
  }
  solver->GetSolveTimer()->Stop();

  numEV_ = eigenValues->GetSize();

  // If no eigenvalue at all converged, just leave
  if(numEV_ == 0) {
    if(recursionCount < 3) {
      if (inputMethod_ == 1) { // inputMethod_ = minVal & maxVal
        minVal_ = minVal_ / 2.0;
        maxVal_ = maxVal_ * 2.0;
      }
      else if (inputMethod_ == 2) {// inputMethod_ = numMode & valueShift
        numMode_ = numMode_ * 2.0;
      }
      CalcValues(++recursionCount);
    }
    else {
      string str = "Did not find any eigenvalue!";
      if(solverType_ == BaseEigenSolver::FEAST)
        str.append("\nTry a smaller or larger interval.");
      EXCEPTION(str);
    }
  }

  LOG_DBG3(buckD) << "CV: eigenvalues = " << eigenValues->ToString();

  Vector<Double> eigenValuesRealPart = GetRealPartOfVector(eigenValues);

  SortModes(false);

  // treat negative eigenvalues
  // this will never happen for feast, which has minVal as input
  if (eigenValuesRealPart.Min() < 0 && recursionCount < 3)
  {
    LOG_DBG3(buckD) << "CV: Found negative eigenvalue";

    assert(eigenValues->IsComplex() == false); // this part is written only for real eigenValues

    // remove negative eigenvalues and corresponding modes
// DO WE REALLY WANT TO DO THIS? negative ev -> force had to be reversed for buckling
//    Vector<Double> tmp1;
//    StdVector<unsigned int> tmp2;
//    for (unsigned int ev = 0; ev < numEV_; ev++) {
//      if (eigenValuesRealPart[ev] >= 0) {
//        tmp1.Push_back(eigenValuesRealPart[ev]);
//        tmp2.push_back(modeOrder_[ev]);
//      }
//    }
//    (*eigenValues) = tmp1; //copy
//    modeOrder_ = tmp2;
//    numEV_ = eigenValues->GetSize();
//
//    LOG_DBG3(buckD) << "CV: eigenValues = " << eigenValues->ToString();
//    LOG_DBG3(buckD) << "CV: modeOrder_ = " << modeOrder_.ToString();

    // if all eigenvalues were negative, search for more
    if (numEV_ == 0 || eigenValuesRealPart.Max() < 0)
    {
      if (inputMethod_ == 1) { // inputMethod_ = minVal & maxVal
        minVal_ = minVal_ / 2.0;
        maxVal_ = maxVal_ * 2.0;
      }
      else if (inputMethod_ == 2) {// inputMethod_ = numMode & valueShift
        numMode_ = numMode_ * 2.0;
      }

      CalcValues(++recursionCount);
    }
  }

  domain->GetOptimization()->context->num_eigenmodes = numEV_;

  // convert eigenvalues to load factors (inverse problem: lf=1/ev, else: lf=ev)
  loadFactors_->Resize(numEV_);
  for (UInt i = 0; i < numEV_; i++)
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

  LOG_DBG(buckD) << "CV: loadFactors = " << loadFactors_->ToString();

  if(!domain->GetOptimization())
    std::cout << "\n++ Finished solving eigenvalue Problem." << std::endl;
}

void BucklingDriver::StoreMode(unsigned int index) {
  Vector<Double> loadFactorsRealPart = GetRealPartOfVector(loadFactors_);
  Double currentLoadFactor = loadFactorsRealPart[modeOrder_[index]];

  // store the eigen mode in algsys sol_
  StdSolveStep *sstep = dynamic_cast<StdSolveStep*>(ptPDE_->GetSolveStep());
  sstep->SetActStep(index);
  sstep->SetActFreq(currentLoadFactor);
  ptPDE_->GetDomain()->GetMathParser()->SetValue(MathParser::GLOB_HANDLER, "f", currentLoadFactor);
  sstep->GetEigenMode(modeOrder_[index]); // this stores the eigen mode in algsys sol_

  LOG_DBG(buckD) << "stored mode " << index;
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

  // console output
  if(!domain->GetOptimization() || domain->GetOptimization()->GetOptimizerType() == Optimization::EVALUATE_INITIAL_DESIGN)
  {
    std::cout << "\n";
    std::cout << " Mode | ";
    std::cout << setw(fieldwidth) << "Load Factor" << " | ";
    std::cout << setw(fieldwidth) << "Error" << " | ";
    if (!isStoredSymmetric_ && solverType_ == BaseEigenSolver::FEAST) {
      std::cout << setw(fieldwidth) << "Imaginary Part" << " | ";
    }
    std::cout << "\n";
  }

  // also log via info node.
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
      std::cout << setw(5) << i + 1 << " | ";
      if (isStoredSymmetric_ || solverType_ == BaseEigenSolver::ARPACK) {
        std::cout << setw(fieldwidth) << std::setprecision(10) << loadFactors_->GetDoubleEntry(modeOrder_[i]) << " | ";
        std::cout << setw(fieldwidth) << errors->GetDoubleEntry(modeOrder_[i]) << " | " << "\n";
      }
      else {
        Complex ev = loadFactors_->GetComplexEntry(modeOrder_[i]);
        std::cout << setw(fieldwidth) << ev.real() << " | ";
        std::cout << setw(fieldwidth) << errors->GetComplexEntry(modeOrder_[i]).real() << " | ";
        std::cout << setw(fieldwidth) << ev.imag() << " | " << "\n";
      }
    }

    // res got cleared above, so we don't overlap with previous in case numConverged is too small
    PtrParamNode mode = res->Get("mode", ParamNode::APPEND);

    mode->Get("nr")->SetValue(i+1); // not the mode but eigenvalue in list
    if (isStoredSymmetric_ || solverType_ == BaseEigenSolver::ARPACK) {
      mode->Get("loadfactor")->SetValue(loadFactors_->GetDoubleEntry(i),15);
      mode->Get("error")->SetValue(errors->GetDoubleEntry(i));
    }
    else {
      mode->Get("loadfactor")->SetValue(loadFactors_->GetComplexEntry(i).real(),15);
      mode->Get("error")->SetValue(errors->GetComplexEntry(i).real());
    }
  }
}

unsigned int BucklingDriver::StoreResults(unsigned int stepNum, double step_val)
{
  LOG_DBG2(buckD) << "SR: step=" << stepNum << " val=" << step_val;

  if(calcModes_)
  {
    Vector<Double> loadFactorsRealPart = GetRealPartOfVector(loadFactors_);

    // check if negative Eigenvalues are present
    for (unsigned int ev = 0; ev < numEV_; ev++) {
      if (loadFactorsRealPart[modeOrder_[ev]] < 0) {
        std::cout << "\n++ WARNING negative proportionality factor will be displayed as positiv in Paraview: " << loadFactorsRealPart[modeOrder_[ev]] << "\n";
      }
    }

    int digs =  boost::lexical_cast<string>((int)numEV_).size() + 2;
    double sig = std::pow((float) 10.0, -digs); // 1e-2 -> 10 ^ -2 -> a compiler complained with simply 10

    for(unsigned int ev = 0; ev < numEV_; ev++)
    {
      // save_value is the "time" value displayed in paraview
      double save_value = -1.0;
      if(domain->GetOptimization() && (domain->GetOptimization()->me->DoHomogenization()
          || domain->GetOptimization()->GetOptimizerType() != Optimization::EVALUATE_INITIAL_DESIGN) )
      {
        // time is step.step_val nr
        save_value = step_val + (ev + 1) * sig; // +1 for one based
        LOG_DBG3(buckD) << "SR: total=" << numEV_ << " digs=" << digs << " sig=" << sig << " count=" << (ev + 1);
      }
      else
        save_value = std::abs(loadFactorsRealPart[modeOrder_[ev]]);

      LOG_DBG2(buckD) << "SR: ev=" << ev << " stepNum=" << stepNum << " save_value=" << save_value;

      // store mode in pde
      if(domain->GetOptimization())
      {
        ErsatzMaterial* em = dynamic_cast<ErsatzMaterial*>(domain->GetOptimization());
        Context* context = domain->GetOptimization()->context;
        Excitation* excite = context->GetExcitation();
        StateSolution* ss = em->forward.Get(excite, NULL, modeOrder_[ev]);
        ss->Write(context->pde); // forward is function NULL
      }
      else
        if(calcModes_)
          StoreMode(ev);

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
      std::cout << "\n++ Finished calculating and storing buckling modes." << std::endl;
  }

  if (!GetResultHandler()->streamOnly)
    return stepNum-1;
  else
    return stepNum;
}

void BucklingDriver::SortModes(bool inAbs) {
  Vector<Double> realPart = GetRealPartOfVector(eigenValues);

  if (inAbs) {
    for (unsigned int i = 0; i < realPart.GetSize(); i++) {
      realPart[i] = std::abs(realPart[i]);
    }
  }

  modeOrder_.Resize(numEV_);
  std::size_t n(0);
  LOG_DBG3(buckD) << "numEV_ = " << numEV_;

  // allocate modeOrder_
  std::generate(std::begin(modeOrder_), std::end(modeOrder_), [&] {return n++;});
  LOG_DBG3(buckD) << "SM: modeOrder_ before sorting = " << modeOrder_.ToString();

  // sort it by value
  std::sort(std::begin(modeOrder_), std::end(modeOrder_), [&](int i1, int i2) {
    if (isInverseProblem_)
      return realPart[i1] > realPart[i2];
    else
      return realPart[i1] < realPart[i2];
  });

  LOG_DBG3(buckD) << "SM: modeOrder_  after sorting = " << modeOrder_.ToString();
}

Vector<Double> BucklingDriver::GetRealPartOfVector(SingleVector* vec) {
  Vector<Double> realPart;

  LOG_DBG3(buckD) << "GRPOV: vec = " << vec->ToString();

  unsigned int sz = vec->GetSize();
  realPart.Resize(sz);

  if (vec->IsComplex())
  {
    for (unsigned int mI = 0; mI < sz; mI++)
      realPart[mI] = vec->GetComplexEntry(mI).real();
  }
  else
    realPart = dynamic_cast<Vector<Double>&>(*vec);

  LOG_DBG3(buckD) << "GRPOV: realPart = " << realPart.ToString();
  return realPart;
}

} // end of namespace
