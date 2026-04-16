#include "EigenValueDriver.hh"

#include <iostream>
#include <iomanip>
#include <cmath>
#include <filesystem>

#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/SingleDriver.hh"
#include "Domain/Domain.hh"
#include "Optimization/Optimization.hh"
#include "Utils/mathParser/mathParser.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/SimState.hh"
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "OLAS/solver/BaseEigenSolver.hh"
#include "OLAS/external/arpack/ArpackEigenSolver.hh"
#include "OLAS/algsys/AlgebraicSys.hh"
#include "MatVec/SBM_Matrix.hh"

#include "PDE/StdPDE.hh"

using std::cout;
using std::setw;
using std::string;

//Definition of EVP Types
#define STANDARD 1
#define GENERALIZED 2
#define QUADRATIC 3

namespace CoupledField {

DEFINE_LOG(evd, "eigenValueDriver")

// forward declaration
class BaseVector;

// ***************
//   Constructor
// ***************
EigenValueDriver::EigenValueDriver(UInt sequenceStep,
                              bool isPartOfSequence,
                              shared_ptr<SimState> state,
                              Domain *domain,
                              PtrParamNode paramNode,
                              PtrParamNode infoNode) :
    SingleDriver(sequenceStep, isPartOfSequence, state, domain, paramNode, infoNode) {

  // set analysistype
  analysis_ = BasePDE::EIGENVALUE;
  // input parameter
  inputMethod_ = NONE;
  numValue_ = 0;
  shiftPoint_ = 0.0;
  shiftPoint_Real_ = 0;
  shiftPoint_Imag_ = 0;
  calcModes_ = false;
  minVal_ = 0.0;
  maxVal_ = 0.0;
  modeNormalization_ = BaseEigenSolver::NONE;
  isStoredSymmetric_ = false;
  rightEigenvectors_ = true;

  //specifying parameter node
  param_ = param_->Get("eigenValue");
  info_ = info_->Get("eigenValue");

  // some variables needed later
  save_step_ = 1; //see calcModes()
  writeAllSteps_ = false;
  evpType_ = NONE;
  matrixC_ = NOTYPE;
  matrixB_ = NOTYPE;
  matrixA_ = NOTYPE;
}

EigenValueDriver::~EigenValueDriver() {

}
// ***************
//   Initiation of EigenValueDriver
//	 Reading parameters of XML-file
// ***************
void EigenValueDriver::Init(bool restart) {
  // read required parameters from parameter node
	if(param_->Has("valuesAround"))
	{
      param_->Get("valuesAround")->GetValue("number", numValue_, ParamNode::INSERT);
      param_->Get("valuesAround")->Get("shiftPoint")->GetValue("Real", shiftPoint_Real_, ParamNode::INSERT);
      param_->Get("valuesAround")->Get("shiftPoint")->GetValue("Imag", shiftPoint_Imag_, ParamNode::PASS);
      shiftPoint_ = Complex(shiftPoint_Real_,shiftPoint_Imag_);
	}
  if(param_->Has("inInterval"))
    {
      param_->Get("inInterval")->GetValue("min", minVal_, ParamNode::INSERT);
      param_->Get("inInterval")->GetValue("max", maxVal_, ParamNode::INSERT);
    }
  solverDefined_ = param_->Has("solverDefined");
  param_->GetValue("allowPostProc", writeAllSteps_, ParamNode::PASS);
  calcModes_  = param_->Has("eigenVectors");

  if (calcModes_) {
    // determine type of mode normalization, and transform into ENUM
    std::string normString = "solver";
    param_->Get("eigenVectors")->GetValue("normalization", normString, ParamNode::PASS);
    if (normString == "none") {
      modeNormalization_ = BaseEigenSolver::NONE;
    }
    else if (normString == "unit") {
      modeNormalization_ = BaseEigenSolver::MAX;
    }
    else if (normString == "norm") {
      modeNormalization_ = BaseEigenSolver::NORM;
    }
    else {
      EXCEPTION("Specified mode normalization '" + normString + "' not implemented");
    }
    // check if left or right eigenvectors should be written in hdf5
    std::string sideString = "";
    param_->Get("eigenVectors")->GetValue("side", sideString, ParamNode::PASS);
    if (sideString == "right") {
      rightEigenvectors_ = true;
    }
    else if (sideString == "left") {
      rightEigenvectors_ = false;
    }
    else {
      EXCEPTION("Select left or right eigenvectors!")
    }
  }

  //determine the problemType and allocate matrices
  if(param_->Has("problemType"))
  {
	  if(param_->Get("problemType")->Has("Standard"))
	  {
		  evpType_ = STANDARD;
		  std::string strMatrixA = "";
		  param_->Get("problemType")->Get("Standard")->GetValue("Matrix", strMatrixA, ParamNode::INSERT);

		  if(strMatrixA == "mass"){
		  		  matrixA_ = MASS;
		  }
		  if(strMatrixA == "damping"){
		  matrixA_ = DAMPING;
		  }
		  if(strMatrixA == "stiffness"){
		  matrixA_ = STIFFNESS;
		  }
	  }

	  else if(param_->Get("problemType")->Has("Generalized"))
	  {
		  evpType_ = GENERALIZED;
		  std::string strMatrixA = "";
		  std::string strMatrixB = "";
		  param_->Get("problemType")->Get("Generalized")->GetValue("aMatrix", strMatrixA, ParamNode::INSERT);
		  param_->Get("problemType")->Get("Generalized")->GetValue("bMatrix", strMatrixB, ParamNode::INSERT);

		  if(strMatrixA == "mass"){
		  		  matrixA_ = MASS;
		  }
		  if(strMatrixA == "damping"){
		  matrixA_ = DAMPING;
		  }
		  if(strMatrixA == "stiffness"){
		  matrixA_ = STIFFNESS;
		  }

		  if(strMatrixB == "mass"){
				  matrixB_ = MASS;
		  }
		  if(strMatrixB == "damping"){
		  matrixB_ = DAMPING;
		  }
		  if(strMatrixB == "stiffness"){
		  matrixB_ = STIFFNESS;
		  }
	  }

	  else if(param_->Get("problemType")->Has("Quadratic"))
	  {
		  evpType_ = QUADRATIC;
		  std::string strMatrixA = "";
		  std::string strMatrixB = "";
		  std::string strMatrixC = "";

		  param_->Get("problemType")->Get("Quadratic")->GetValue("quadratic", strMatrixA, ParamNode::INSERT);
		  param_->Get("problemType")->Get("Quadratic")->GetValue("linear", strMatrixB, ParamNode::INSERT);
		  param_->Get("problemType")->Get("Quadratic")->GetValue("constant", strMatrixC, ParamNode::INSERT);

		  if(strMatrixA == "mass"){
		  	matrixA_ = MASS;
		  }
		  if(strMatrixA == "damping"){
		    matrixA_ = DAMPING;
		  }
		  if(strMatrixA == "stiffness"){
		    matrixA_ = STIFFNESS;
		  }

		  if(strMatrixB == "mass"){
				matrixB_ = MASS;
		  }
		  if(strMatrixB == "damping"){
		    matrixB_ = DAMPING;
		  }
		  if(strMatrixB == "stiffness"){
		    matrixB_ = STIFFNESS;
		  }

		  if(strMatrixC == "mass"){
				 matrixC_ = MASS;
		  }
		  if(strMatrixC == "damping"){
		    matrixC_ = DAMPING;
		  }
		  if(strMatrixC == "stiffness"){
		    matrixC_ = STIFFNESS;
		  }
	  }
  }
  else
  {
      EXCEPTION("problemType not defined!");
  }

  // set definition of PDE relevant matrices
  InitializePDEs();

  // read sorting information
  Enum<SortingMethodType> sortingMethod;
  sortingMethod.SetName("EigenValueDriver::SortingMethodType");
  sortingMethod.Add(ABS,"abs");
  sortingMethod.Add(REAL,"real");
  sortingMethod.Add(IMAG,"imag");
  sortingMethod_ = sortingMethod.Parse(param_->Get("sort")->As<std::string>());

}
// ***************
//   Solve problem step
// ***************
void EigenValueDriver::SolveProblem() {
  // define result handler
  ResultHandler *resHandler = domain_->GetResultHandler();

  // see comments in StaticDriver::SolveProblem() for the interplay with optimization
  numEV_ = 1;
  if (!domain->GetOptimization())
    resHandler->BeginMultiSequenceStep(sequenceStep_, analysis_, numEV_); // optimization does it by itself

  if (writeAllSteps_ || isPartOfSequence_)
    simState_->BeginMultiSequenceStep(sequenceStep_, analysis_);

  // trigger calculation
  ptPDE_->WriteGeneralPDEdefines();
  BaseSolveStep *step = ptPDE_->GetSolveStep();

  // set the mode normalization
  dynamic_cast<StdSolveStep*>(step)->GetAlgSys()->GetEigenSolver()->SetModeNormalization(modeNormalization_);

  // initialize algebraic system
  StdSolveStep *sstep = dynamic_cast<StdSolveStep*>(step);
  sstep->GetAlgSys()->InitSol();
  sstep->GetAlgSys()->InitMatrix();
  sstep->GetAssemble()->AssembleMatrices();
  sstep->GetAlgSys()->ExportLinSys(true, false, false); // if asked, export matrices

  // define input method, which is needed in CalcValues()
  // inputMethod_ is 1 for minVal & maxVal, 2 for numModes & shiftMode, 3 for solver defined (currently only FEAST custom contour)
  if (minVal_ < maxVal_) {
    inputMethod_ = 1;
  }
  else if (minVal_ > maxVal_) {
    EXCEPTION("minVal is bigger than maxVal. Please correct XML input file.");
  }
  else if (numValue_ > 0) {
    inputMethod_ = 2;
  }
  else if ( solverDefined_ ) {
    inputMethod_ = 3;
  }
  else {
    // this case should not be possible in the XML schema
    EXCEPTION("Input method cannot be determined.")
  }

  // actually solve problem
  CalcEigenValues();
  // populate errBounds
  if( errBounds_.GetSize() == 0 && errBoundsComplex_.GetSize() > 0 ) { // error bounds are complex, set the real ones to abs
    errBounds_.Resize(errBoundsComplex_.GetSize());
    for (unsigned int i = 0; i < errBoundsComplex_.GetSize(); i++) {
      errBounds_[i] = std::abs(errBoundsComplex_[i]);
    }
  } else if ( errBounds_.GetSize() > 0 && errBoundsComplex_.GetSize() == 0 ) { // err bounds are real, set complex accordingly
    errBoundsComplex_.Resize(errBounds_.GetSize());
    for (unsigned int i = 0; i < errBounds_.GetSize(); i++) {
      errBoundsComplex_[i] = Complex(errBounds_[i],0.0);
    }
  } else {
    EXCEPTION("Both error bounds are empty? No Eigenvalues found?")
  }

  numEV_ = eigenValuesComplex_.GetSize();
  SortModes();
  PrintResult();
  ToInfo();
  // store to hdf5
  if (calcModes_) {
	  StoreResults(1, -1.0);
  }
  // store to textfile
  PtrParamNode els = sstep->GetAlgSys()->GetExportLinSysParam();
  if (els) {
    if(els->Get("solution")->As<bool>()) {
      BaseMatrix::OutputFormat vec_format = MatrixOutputFormatEnum.Parse(els->Get("vecFormat")->As<std::string>());
      std::string base = els->Has("baseName") ? els->Get("baseName")->As<std::string>() : progOpts->GetSimName();
      if(domain->GetDriver()->GetAnalysisId().ToString(true) != ""){
        base += "_" + domain->GetDriver()->GetAnalysisId().ToString(true);
      }
      for (UInt i=0; i<numEV_; i++) {
        Vector<Complex> mode;
        // TU Wien Variant with normalized eigenmodes
        sstep->GetAlgSys()->GetEigenSolver()->GetNormalizedEigenMode(modeOrder_[i],mode);
        mode.Export( base + "_mode_" + std::to_string(i+1), vec_format);
        sstep->GetAlgSys()->GetEigenSolver()->GetNormalizedEigenMode(modeOrder_[i],mode,false);
        mode.Export( base + "_mode-left_" + std::to_string(i+1), vec_format);
      }
    }
  }

  // finish step
  handler_->FinishMultiSequenceStep();

  if (!isPartOfSequence_)
    handler_->Finalize(); // to be called only once in a HDF5 lifetime!

  if (writeAllSteps_ || isPartOfSequence_)
    simState_->FinishMultiSequenceStep(true);
}
// ***************
//   Calculation of the eigenvalues
// ***************
void EigenValueDriver::CalcEigenValues() {
  BaseSolveStep *step = ptPDE_->GetSolveStep();
  StdSolveStep *sstep = dynamic_cast<StdSolveStep*>(step);

  // get active solver
  BaseEigenSolver *solver = sstep->GetAlgSys()->GetEigenSolver();

  SBM_Matrix *matrixA;
  SBM_Matrix *matrixB;
  SBM_Matrix *matrixC;
  // switch matrix setup based on EVP type
  switch(evpType_) {
    // Standard EVP
    case STANDARD:
    {
      matrixA = sstep->GetAlgSys()->GetMatrix(matrixA_);
      if (matrixA->GetNumCols()>1) {
          EXCEPTION("only implemented for SBM matrices with a single block")
      }
      assert(matrixA->GetNumCols() == matrixA->GetNumRows());
      solver->Setup(*(matrixA->GetPointer(0, 0)), false);
      break;
    }
    //Generalized EVP
    case GENERALIZED:
    {
      matrixA = sstep->GetAlgSys()->GetMatrix(matrixA_);

      matrixB = sstep->GetAlgSys()->GetMatrix(matrixB_);
      if (matrixA->GetNumCols()>1) {
          EXCEPTION("only implemented for SBM matrices with a single block")
      }
      assert(matrixA->GetNumCols() == matrixA->GetNumRows());
      assert(matrixB->GetNumCols() == matrixB->GetNumRows());
      solver->Setup(*(matrixA->GetPointer(0, 0)), *(matrixB->GetPointer(0, 0)), false);
      bool isReal;
      solver->CheckMatrix(isReal, isStoredSymmetric_, *(matrixB->GetPointer(0, 0)));
      break;
    }

    //Quadratic EVP
    case QUADRATIC:
    {
      matrixA = sstep->GetAlgSys()->GetMatrix(matrixA_);
      matrixB = sstep->GetAlgSys()->GetMatrix(matrixB_);
      matrixC = sstep->GetAlgSys()->GetMatrix(matrixC_);
      if (matrixA->GetNumCols()>1) {
          EXCEPTION("only implemented for SBM matrices with a single block")
      }
      //Check Matrix dimensions
      assert(matrixA->GetNumCols() == matrixA->GetNumRows());
      assert(matrixB->GetNumCols() == matrixB->GetNumRows());
      assert(matrixC->GetNumCols() == matrixC->GetNumRows());

      solver->Setup(*(matrixC->GetPointer(0, 0)), *(matrixB->GetPointer(0, 0)), *(matrixA->GetPointer(0, 0)));
      break;
    }
  }

  // check if complex eigenvalues
  bool complexEV = solver->HasComplexEigenvalues();

  // switch CalcEigenValues method based on input
  switch(inputMethod_){
    // minVal + maxVal
    case 1: {
      if (complexEV) {
        sstep->GetAlgSys()->GetEigenSolver()->CalcEigenValues(eigenValuesComplex_,errBoundsComplex_,minVal_,maxVal_);
      }
      else {
        sstep->GetAlgSys()->GetEigenSolver()->CalcEigenValues(eigenValues_,errBounds_,minVal_,maxVal_);
      }
      break;
    }
    // numMode + freqShift
    case 2: {
      solver->CalcEigenValues(eigenValuesComplex_, errBounds_, numValue_, shiftPoint_);
      break;
    }
    // solver defined
    case 3: {
      if (complexEV) {
        sstep->GetAlgSys()->GetEigenSolver()->CalcEigenValues(eigenValuesComplex_,errBoundsComplex_);
      }
      else {
        sstep->GetAlgSys()->GetEigenSolver()->CalcEigenValues(eigenValues_,errBounds_);
      }
      break;
    }
  }
}
// ***************
//   Printing the eigenvalues
// ***************
void EigenValueDriver::PrintResult() {
  int n = 23; // field width
  cout << "\n";
  cout << " Mode | ";
  cout << setw(n) << "Real Part" << " | ";
  cout << setw(n) << "Error Bounds" << " | ";
  if (eigenValuesComplex_.GetSize() > 0 ) {
    cout << setw(n) << "Imaginary Part" << " | ";
  }
  cout << "\n";
  for (unsigned int i = 0; i < modeOrder_.GetSize(); i++) {
    cout << setw(5) << i + 1 << " | ";
    if (eigenValues_.GetSize() > 0) {
      cout << setw(n) << eigenValues_[modeOrder_[i]] << " | ";
      cout << setw(n) << errBounds_[modeOrder_[i]] << " | " << "\n";
    }
    else {
      cout << setw(n) << eigenValuesComplex_[modeOrder_[i]].real() << " | ";
      cout << setw(n) << errBounds_[modeOrder_[i]] << " | ";
      cout << setw(n) << eigenValuesComplex_[modeOrder_[i]].imag() << " | " << "\n";
    }
  }
}
// ***************
//   Saving the eigenvalues to info.xml
// ***************
void EigenValueDriver::ToInfo() {
  // saving of eigenvalues to .info.xml
  // this is especially necessary for complex eigenvalues, since only their magnitude is saved in the hdf5 file
  PtrParamNode res = info_->Get("result", ParamNode::APPEND);
  PtrParamNode evals = res->Get("eigenvalues", ParamNode::APPEND);
  for (unsigned int i = 0; i < modeOrder_.GetSize(); i++) {
    if (eigenValues_.GetSize() > 0) {
      // save in .info.xml
      PtrParamNode mode = evals->Get("mode", ParamNode::APPEND);
      mode->Get("nr")->SetValue(i+1);
      mode->Get("real")->SetValue(eigenValues_[modeOrder_[i]]);
    }
    else {
      // save in .info.xml
      PtrParamNode mode = evals->Get("mode", ParamNode::APPEND);
      mode->Get("nr")->SetValue(i+1);
      mode->Get("real")->SetValue(eigenValuesComplex_[modeOrder_[i]].real());
      mode->Get("imag")->SetValue(eigenValuesComplex_[modeOrder_[i]].imag());
    }
  }
}
unsigned int EigenValueDriver::StoreResults(unsigned int stepNum, double step_val)
{
  // stepNum and step_val are ignored
  LOG_DBG(evd) << "SR step=" << stepNum << " val=" << step_val;

  // generates a index-array modeOrder_ containing the mode indices sorted by ascending Frequency value
  SortModes();
  for(unsigned int fi=0; fi < eigenValuesComplex_.GetSize(); fi++)
  {
    // Phase 2: calculate eigenmodes
    if(save_step_)
    {
      ptPDE_->GetSolveStep()->SetActStep(fi);
      ptPDE_->GetSolveStep()->SetActFreq(eigenValuesComplex_[modeOrder_[fi]].imag());
      ptPDE_->GetDomain()->GetMathParser()->SetValue(MathParser::GLOB_HANDLER, "f", eigenValuesComplex_[modeOrder_[fi]].imag());
      if (rightEigenvectors_) {
        ptPDE_->GetSolveStep()->GetEigenMode(modeOrder_[fi]); // store right eigenmode
      }
      else{
        ptPDE_->GetSolveStep()->GetEigenMode(modeOrder_[fi], false); // store left eigenmode
      }
      

      // stupid paraview needs an increasing series of save_value :(

      double save_value = -1.0;
      save_value = GetSortValue(eigenValuesComplex_[modeOrder_[fi]]);

      LOG_DBG(evd) <<  "fi=" << fi << " save_step_=" << save_step_ << " save_value=" << save_value;

      handler_->BeginStep(save_step_, save_value);
      ptPDE_->WriteResultsInFile(save_step_, save_value);
      handler_->FinishStep();

      if(writeAllSteps_ || isPartOfSequence_)
        simState_->WriteStep(save_step_, save_value);

      save_step_++;
    }
  }
  for(unsigned int fi=0; fi < eigenValues_.GetSize(); fi++)
  {
    // Phase 2: calculate eigenmodes
    if(save_step_)
    {
      ptPDE_->GetSolveStep()->SetActStep(fi);
      ptPDE_->GetSolveStep()->SetActFreq(eigenValues_[modeOrder_[fi]]);
      ptPDE_->GetDomain()->GetMathParser()->SetValue(MathParser::GLOB_HANDLER, "f", eigenValues_[modeOrder_[fi]]);
      ptPDE_->GetSolveStep()->GetEigenMode(modeOrder_[fi]); // this stores the eigen mode result in AlgSys's sol_

      // stupid paraview needs an increasing series of save_value :(

      double save_value = -1.0;
      save_value =std::abs(eigenValues_[modeOrder_[fi]]);

      LOG_DBG(evd) <<  "fi=" << fi << " save_step_=" << save_step_ << " save_value=" << save_value;

      handler_->BeginStep(save_step_, save_value);
      ptPDE_->WriteResultsInFile(save_step_, save_value);
      handler_->FinishStep();

      if(writeAllSteps_ || isPartOfSequence_)
        simState_->WriteStep(save_step_, save_value);

      save_step_++;
    }
  }

  if (!GetResultHandler()->streamOnly)
    return stepNum-1;
  else
    return stepNum;
}
// ***************
//   Sort modes by imaginary part
// ***************
void EigenValueDriver::SortModes() {
  Vector<Double> RealeigenValues;
  Vector<Complex> ComplexeigenValues;
  if (eigenValues_.GetSize()>0) {
    RealeigenValues = eigenValues_;
    modeOrder_.Resize(RealeigenValues.GetSize());
    std::size_t n(0);
    // allocate modeOrder_
    std::generate(std::begin(modeOrder_), std::end(modeOrder_), [&] {return n++;});
    // sort it by value
    std::sort(std::begin(modeOrder_), std::end(modeOrder_), [&](int i1, int i2) {
      return RealeigenValues[i1] < RealeigenValues[i2];});
  }
  else {
	  ComplexeigenValues.Resize(eigenValuesComplex_.GetSize());
	  ComplexeigenValues = eigenValuesComplex_;
	  modeOrder_.Resize(ComplexeigenValues.GetSize());
	  std::size_t n(0);
	  // allocate modeOrder_
	  std::generate(std::begin(modeOrder_), std::end(modeOrder_), [&] {return n++;});
	  // sort it by value
	  std::sort(std::begin(modeOrder_), std::end(modeOrder_), [&](int i1, int i2) {
	  return GetSortValue(ComplexeigenValues[i1]) < GetSortValue(ComplexeigenValues[i2]);});
  }
  }

void EigenValueDriver::SetToStepValue(UInt stepNum, Double stepVal) {
  // ensure that this method is only called if simState has input
  if (! simState_->HasInput()) {
    EXCEPTION("Can only set external time step, if simulation state " << "is read from external file");
  }
  // Set current eigenvalue in the mathParser
  domain_->GetMathParser()->SetValue(MathParser::GLOB_HANDLER, "f", stepVal);
  domain_->GetMathParser()->SetValue(MathParser::GLOB_HANDLER, "step", stepNum);

}

double EigenValueDriver::GetSortValue(Complex entry)
{
  switch (sortingMethod_) {
    case REAL:
      return entry.real();
    case IMAG:
      return entry.imag();
    case ABS:
      return std::abs(entry);
    default:
      return std::abs(entry);
  }
}

}




