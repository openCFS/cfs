#include "BucklingDriver.hh"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <boost/filesystem.hpp>

#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/SingleDriver.hh"
#include "Domain/Domain.hh"
#include "Optimization/Optimization.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/SimState.hh"
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "OLAS/solver/BaseEigenSolver.hh"
#include "OLAS/algsys/AlgebraicSys.hh"
#include "MatVec/SBM_Matrix.hh"

#include "PDE/StdPDE.hh"

using std::cout;
using std::setw;
using std::string;

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
                              PtrParamNode infoNode) :
    SingleDriver(sequenceStep, isPartOfSequence, state, domain, paramNode, infoNode) {

  // set analysistype
  analysis_ = BasePDE::BUCKLING;
  // input parameter
  numMode_ = 0;
  valueShift_ = 0.0;
  calcModes_ = false;
  minVal_ = 0.0;
  maxVal_ = 0.0;
  modeNormalization_ = BaseEigenSolver::NONE;
  isSymmetrical_ = false;

  //specifying parameter node
  param_ = param_->Get("buckling");
  info_ = info_->Get("buckling");

  // some variables needed later
  save_step_ = 1; //see calcModes()
  writeAllSteps_ = false;
}

BucklingDriver::~BucklingDriver() {

}

void BucklingDriver::Init(bool restart) {
  // read required parameters from parameter node
  param_->GetValue("numModes", numMode_, ParamNode::INSERT);
  param_->GetValue("valueShift", valueShift_, ParamNode::INSERT);
  param_->GetValue("minVal", minVal_, ParamNode::INSERT);
  param_->GetValue("maxVal", maxVal_, ParamNode::INSERT);
  param_->GetValue("allowPostProc", writeAllSteps_, ParamNode::PASS);
  calcModes_  = param_->Has("calcModes");

  if(calcModes_){
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

double BucklingDriver::GetPropFactor(unsigned int idx) const {  // maybe this function could be removed
  return eigenValues_[idx];;
}

void BucklingDriver::SolveProblem() {

  // define Result Handler
  ResultHandler *resHandler = domain_->GetResultHandler();

  //TODO integrate optimization
  // see comments in StaticDriver::SolveProblem() for the interplay with optimization
  UInt numSteps = 1;
  if (!domain->GetOptimization())
    resHandler->BeginMultiSequenceStep(sequenceStep_, analysis_, numSteps); // optimization does it by itself

  if (writeAllSteps_ || isPartOfSequence_)
    simState_->BeginMultiSequenceStep(sequenceStep_, analysis_);

  // Trigger calculation
  ptPDE_->WriteGeneralPDEdefines();
  BaseSolveStep *step = ptPDE_->GetSolveStep();

  // set the mode normalization
  dynamic_cast<StdSolveStep*>(step)->GetAlgSys()->GetEigenSolver()->SetModeNormalization(modeNormalization_);

  // initialize Algebraic System
  StdSolveStep *sstep = dynamic_cast<StdSolveStep*>(step);
  sstep->GetAlgSys()->InitSol();
  sstep->GetAlgSys()->InitMatrix();
  sstep->GetAssemble()->AssembleMatrices();
  sstep->GetAlgSys()->ExportLinSys(true, false, false); // if asked, export matrices

  // define input method, which is needed in CalcValues()
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
    EXCEPTION("this case should not be possible in the XML schema")
  }

  // actually solve problem
  CalcValues();
  SortModes(false);
  PrintResult();
  if (calcModes_) {
    CalcMode();
  }

  // finish step
  handler_->FinishMultiSequenceStep();

  if (!isPartOfSequence_)
    handler_->Finalize(); // to be called only once in a HDF5 lifetime!

  if (writeAllSteps_ || isPartOfSequence_)
    simState_->FinishMultiSequenceStep(true);
}

void BucklingDriver::CalcValues() {

  BaseSolveStep *step = ptPDE_->GetSolveStep();
  StdSolveStep *sstep = dynamic_cast<StdSolveStep*>(step);
  //StdSolveStep *sstep = ptPDE_->GetSolveStep();

  // get active solver
  BaseEigenSolver *solver = sstep->GetAlgSys()->GetEigenSolver();

  // get type of solver
  BaseEigenSolver::EigenSolverType solverType;
  solverType = solver->GetEigenSolverName();

  // get needed matrices
  SBM_Matrix *tangStiffMat = sstep->GetAlgSys()->GetMatrix(TANGENTIAL_STIFFNESS);
  SBM_Matrix *stiffMat = sstep->GetAlgSys()->GetMatrix(STIFFNESS);

  // check matrices dimensions
  UInt i = tangStiffMat->GetNumCols();
  if (i > 1) {
    EXCEPTION("only implemented for SBM matrices with a single block")
  }
  assert(tangStiffMat->GetNumCols() == tangStiffMat->GetNumRows());
  assert(stiffMat->GetNumCols() == stiffMat->GetNumRows());

  //check storage type
  bool isReal = false;
  solver->CheckMatrix(isReal, isSymmetrical_, *(stiffMat->GetPointer(0, 0)));


  // unfortunately the implemented eigenSolvers have serious usage restrictions.
  // currently only FEAST with non-symmetric matrices gives relabel solutions
  if(solverType == BaseEigenSolver::FEAST && isSymmetrical_ == true){
    EXCEPTION("Wrong Matrix storage type. Please define an Non-symmetrical storage type.");
  }
  else if(solverType == BaseEigenSolver::FEAST && inputMethod_ != 1){
    EXCEPTION("FEAST eigenSolver does currently not support input as numModes & shiftMode.");
  }
  else if(solverType == BaseEigenSolver::ARPACK){
    EXCEPTION("The needed dsdrv5 driver for ARPACK solver has not been implemented, yet.");
  }
  else if(solverType == BaseEigenSolver::PHIST){
    EXCEPTION("PHIST is currently not supported.");
  }

/* TODO FEAST
 * - add inputMethod_ 2, as numModes / shiftMode
 */
/* TODO ARPACK
 * - add inputMethod_ 1, as min/maxValue
 * - add CalcEigenValues()
 * - add dsdrv5 driver, and create a "setProblemType()"--like function to change iparam[6] = 4 and use sad driver
 */
/* TODO PHIST
 * - add non-symmetric Mode (?)
 */



  // unfortunately the implemented eigenSolvers do not conform to a single scheme
  if (inputMethod_ == 1) { // inputMethod_ = minVal & maxVal
    solver->Setup(*(stiffMat->GetPointer(0, 0)), *(tangStiffMat->GetPointer(0, 0)), false);

    if (isSymmetrical_) {
      solver->CalcEigenValues(eigenValues_, errBounds_, minVal_, maxVal_);
    }
    else {
      solver->CalcEigenValues(eigenValuesComplex_, errBoundsComplex_, minVal_, maxVal_);
    }

  }
  else if (inputMethod_ == 2) {// inputMethod_ = numMode & valueShift
    solver->Setup(*(stiffMat->GetPointer(0, 0)), *(tangStiffMat->GetPointer(0, 0)), numMode_, valueShift_, true, false);

    if (isSymmetrical_) {
      solver->CalcEigenValues(eigenValues_, errBounds_, numMode_, valueShift_);
    }
    else {
      solver->CalcEigenValues(eigenValuesComplex_, errBoundsComplex_, numMode_, valueShift_);
    }

  }
  else {
    EXCEPTION("input method not known")
  }
  std::cout << "\n" << "++ " << "Finished solving eigenvalue Problem." << "\n";
}

void BucklingDriver::PrintResult() {
  // print results to console

  int n = 23; // field width
  cout << "\n";
  cout << " Mode | ";
  cout << setw(n) << "Proportionality Factor" << " | ";
  cout << setw(n) << "Error Bounds" << " | ";
  if (isSymmetrical_ == false){
    cout << setw(n) << "Imaginary Part" << " | ";
  }
  cout << "\n";
  for (unsigned int i = 0; i < modeOrder_.GetSize(); i++) {
    cout << setw(5) << i + 1 << " | ";
    if (isSymmetrical_) {
      cout << setw(n) << eigenValues_[modeOrder_[i]] << " | ";
      cout << setw(n) << errBounds_[modeOrder_[i]] << " | " << "\n";
    }
    else {
      cout << setw(n) << eigenValuesComplex_[modeOrder_[i]].real() << " | ";
      cout << setw(n) << errBoundsComplex_[modeOrder_[i]].real() << " | ";
      cout << setw(n) << eigenValuesComplex_[modeOrder_[i]].imag() << " | " << "\n";
    }
  }
}

void BucklingDriver::CalcMode() {
  // calculate buckling modes

  BaseSolveStep *step = ptPDE_->GetSolveStep();
  StdSolveStep *sstep = dynamic_cast<StdSolveStep*>(step);

  Vector<Double> eigenValuesForModes;
  unsigned int numEigenValues = 0;

  if (isSymmetrical_) {
    numEigenValues = eigenValues_.GetSize();
    eigenValuesForModes = eigenValues_;
  }
  else{
    numEigenValues = eigenValuesComplex_.GetSize();
    eigenValuesForModes.Resize(numEigenValues);
    for (unsigned int mI = 0; mI < numEigenValues; mI++) {
      eigenValuesForModes[mI] = eigenValuesComplex_[mI].real();
    }
  }

  //Sort Eigenvalues with absolute values. Paraview can not display negative Eigenvalues.
  SortModes(true);

  //Check if negative Eigenvalues are present
  for (unsigned int mI = 0; mI < numEigenValues; mI++) {
    if(eigenValuesForModes[modeOrder_[mI]] < 0){
      std::cout << "\n" << "++ " << "WARNING negative proportionality factor will be displayed as positiv in Paraview: " << eigenValuesForModes[mI] << "\n";
    }
  }

  for (unsigned int mI = 0; mI < numEigenValues; mI++) {

    Double currentEigenValue = eigenValuesForModes[modeOrder_[mI]];

    sstep->SetActStep(mI);
    sstep->SetActFreq(currentEigenValue);
    ptPDE_->GetDomain()->GetMathParser()->SetValue(MathParser::GLOB_HANDLER, "eigenValue", currentEigenValue);
    sstep->GetEigenMode(modeOrder_[mI]); // this stores the eigen mode result in AlgSys's sol_


    double save_value = std::abs(currentEigenValue); // Paraview needs an increasing series of save_value :(

    LOG_DBG(buckD)<< " mI=" << mI << " save_step_=" << save_step_ << " save_value=" << save_value;

    handler_->BeginStep(save_step_, save_value);
    ptPDE_->WriteResultsInFile(save_step_, save_value);
    handler_->FinishStep();

    if (writeAllSteps_ || isPartOfSequence_)
      simState_->WriteStep(save_step_, save_value);

    save_step_++;

  }
  //Reverse order of modeOrder_
  SortModes(false);
  std::cout << "\n" << "++ " << "Finished calculating and storing buckling modes." << "\n";
}

} // end of namespace
