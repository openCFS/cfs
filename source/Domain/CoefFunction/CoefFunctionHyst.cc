#include "CoefFunctionHyst.hh"

// classes for function / spline approximation
#include "Materials/Models/Preisach.hh"
#include "Materials/Models/VectorPreisachv10.hh"
#include "FeBasis/FeFunctions.hh"
#include "FeBasis/FeSpace.hh"
#include "Forms/Operators/BaseBOperator.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include <cmath>

namespace CoupledField {

  DECLARE_LOG(coeffcthyst)
  DEFINE_LOG(coeffcthyst, "coeffcthyst")

CoefFunctionHyst::CoefFunctionHyst( BaseMaterial* const material,
                    shared_ptr<ElemList> actSDList,
                    PtrCoefFct dependency,
					SubTensorType tensorType, MaterialType matType) : CoefFunction() {

  // this type of coefficient is nonlinear (i.e. solution dependent)
  dimType_ = VECTOR;
  dependType_ = SOLUTION;
  isAnalytic_ = false;
  isComplex_  = false;
  dependCoef_ = dependency;
  tensorType_ = tensorType;
  matType_    = matType;
  material_   = material;
  tol_ = 1e-14;

  /*
   * returnFreeFieldTensor_ -> if true, getTensor returns either eps0/nu0
   * returnInitialTensor_ -> if true, getTensor returns starting value for eps/nu
   *
   * if both are false -> getTensor returns eps0/nu0 + deltaMat
   */
  returnFreeFieldTensor_ = false;
  returnInitialTensor_ = true;

  /*
   * if returnZeroValues_ = true, getVector and getScalar will return 0
   */
  returnZeroValues_ = false;

  /*
   * if useDeltaY = true, deltaMat = (Y_new - Y_old)/denominator
   * else, deltaMat = Y_new/denominator
   *
   * if useDeltaX = true, deltaMat = nominator/(X_new - X_old)
   * else, DeltaMat = nominator/X_new
   *
   * if useNextToLastTS_ = true, Y_old = YnextToLastTS_/YnextToLastTSVEC_
   *                         and X_old = XnextToLastTS_/XnextToLastTSVEC_
   * else, Y_old = YpreviousIt_/YpreviousItVEC_
   *   and X_old = XpreviousIt_/XpreviousItVEC_
   */
  useDeltaY_ = true;
  useDeltaX_ = true;
  useNextToLastTS_ = true;

  /*
   * if hystMemoryLocked_ = true, changes to hysteresis operator will be per-
   *                          formed on a temporal storage only; no permanent
   *                          change to hysteresis memory
   * if hystDirectionLocked_ = true, the rotation direction of vector Preisach
   *                          models will not get updated; the old rotation
   *                          state is kept; the switching state is changed as
   *                          usual
   */
  hystMemoryLocked_ = false;
  hystDirectionLocked_ = false;

  allowBMP_ = false;
  Init(material, actSDList);
}

CoefFunctionHyst::~CoefFunctionHyst(){
  ;
}

void CoefFunctionHyst::Init( BaseMaterial* const material,
                             shared_ptr<ElemList> actSDList) {

	// set parameters
	Double Xsat, Ysat;
	material->GetScalar(Xsat, X_SATURATION, Global::REAL);
	material->GetScalar(Ysat, Y_SATURATION, Global::REAL);
	Matrix<Double> weights;
	material->GetTensor(weights,  PREISACH_WEIGHTS, Global::REAL);

	xSat_ = Xsat;
	ySat_ = Ysat;
	numRows_ = weights.GetNumRows();

	bool isVirgin = true;
	UInt numElemSD = actSDList->GetSize();

	std::string dimTypeStr;
	material->GetScalar(dimTypeStr, PREISACH_DIM);

	// use variable methodType_ to distinguish between scalar and vector modell
	// do not confuse this with dimType_ !
	if(dimTypeStr == "SCALAR"){
	  methodType_ = SCALAR;
	} else if(dimTypeStr == "VECTOR"){
	  methodType_ = VECTOR;
	}

	if(methodType_ == SCALAR){
    //get direction
    std::string str;
    material->GetScalar(str, P_DIRECTION);
    Directions dir;
    String2Enum(str,dir);
    dirP_ = dir;

    hyst_ = new Preisach(numElemSD, Xsat, Ysat, weights, isVirgin);

    dim_ = dependCoef_->GetVecSize();
    //dim_ = 1;
    evalVersion_ = 0;

	} else if(methodType_ == VECTOR){

	  dim_ = dependCoef_->GetVecSize();

	  Double angDistance;
	  Matrix<Double> easyAxis_Matrix;
	  Vector<Double> easyAxis = Vector<Double>(dim_);
    int printOut;
    int bmpResolution;

	  material->GetScalar(rotRes_, ROT_RESISTANCE, Global::REAL);
	  material->GetScalar(angDistance, ANG_DISTANCE, Global::REAL);

	  material->GetScalar(printOut,PRINT_PREISACH);
	  material->GetScalar(bmpResolution,PRINT_PREISACH_RESOLUTION);

	  /*
	   * if printOut > 0 -> activate output of overlaid switching and rotation state; output every printOut timestep
	   * bmpResolution -> number of pixels (std = 1000)
	   */
	  printOut_ = (UInt)printOut;
	  bmpResolution_ = (UInt)bmpResolution;

    int eval;
    material->GetScalar(eval, EVAL_VERSION);

    evalVersion_ = (UInt) eval;
    bool classical;

    if(evalVersion_ == 1){
      classical = true; // original vector preisach model -> sutor2012

      hyst_ = new VectorPreisachv10_ListApproach(numElemSD, Xsat, Ysat,
                                                 weights, rotRes_, dim_, isVirgin,
                                                 classical, angDistance);
    } else if(evalVersion_ == 2){
      classical = false; // revised vector preisach model -> sutor2015

      hyst_ = new VectorPreisachv10_ListApproach(numElemSD, Xsat, Ysat,
                                                 weights, rotRes_, dim_, isVirgin,
                                                 classical, angDistance);
    } else if(evalVersion_ == 10){
      classical = true; // original vector preisach model -> sutor2015; matrix based implementation

      hyst_ = new VectorPreisachv10_MatrixApproach(numElemSD, Xsat, Ysat,
                                                 weights, rotRes_, dim_, isVirgin,
                                                 classical, angDistance);
    } else if(evalVersion_ == 20){
      classical = false; // revised vector preisach model -> sutor2015; matrix based implementation

      hyst_ = new VectorPreisachv10_MatrixApproach(numElemSD, Xsat, Ysat,
                                                 weights, rotRes_, dim_, isVirgin,
                                                 classical, angDistance);
    } else {
      EXCEPTION("evalVersion has to be one of the following: \n "
          "1: classical vector model (sutor2012) \n"
          "2: revised vector model (sutor2015) [DEFAULT] \n"
          "10: classical vector model (sutor2012) - Matrix implementation, only for reference \n"
          "20: revised vector model (sutor2015) - Matrix implementation, only for reference \n")
    }
	}

	// set map: global to local element number
	EntityIterator it = actSDList->GetIterator();
	UInt iel = 0;
	UInt globalElNr;
	for ( it.Begin(); !it.IsEnd(); it++, iel++) {
		globalElNr = it.GetElem()->elemNum;
		globalElem2Local_[globalElNr] = iel;
	}

	//store subdomain list of elements
	SDList_ = actSDList;

	//allocate memory for previous results, needed for the
	//effective material parameter formulation
	if(methodType_ == SCALAR){
    XpreviousEval_.Resize(numElemSD);
    YpreviousEval_.Resize(numElemSD);
    XpreviousEval_.Init();
    YpreviousEval_.Init();

    /*
     * initialize previous x values to a normally unreachable value, such that we can trigger
     * a first calcuation in each case
     */
    //XpreviousIt_.Init(-inf);
    /*
     * In new version, we use XpreviousIt_ directly during the first computation
     * so it is not allowed to be inf at the point; use 0.0 instead;
     */
    XpreviousIt_.Resize(numElemSD);
    YpreviousIt_.Resize(numElemSD);
    XpreviousIt_.Init();
    YpreviousIt_.Init();

    XnextToLastTS_.Resize(numElemSD);
    YnextToLastTS_.Resize(numElemSD);
    XnextToLastTS_.Init();
    YnextToLastTS_.Init();

	} else if(methodType_ == VECTOR){
	  XpreviousEvalVEC_ = new Vector<Double>[numElemSD];
	  YpreviousEvalVEC_ = new Vector<Double>[numElemSD];

    XpreviousItVEC_ = new Vector<Double>[numElemSD];
    YpreviousItVEC_ = new Vector<Double>[numElemSD];

    XnextToLastTSVEC_ = new Vector<Double>[numElemSD];
    YnextToLastTSVEC_ = new Vector<Double>[numElemSD];

    dXpreviousItVEC_ = new Vector<Double>[numElemSD];
    dYpreviousItVEC_ = new Vector<Double>[numElemSD];

    XcurrentItVEC_ = new Vector<Double>[numElemSD];
    YcurrentItVEC_ = new Vector<Double>[numElemSD];

	  for(UInt i = 0; i < numElemSD; i++){
	    XpreviousEvalVEC_[i].Resize(dim_);
	    XpreviousEvalVEC_[i].Init();

	    YpreviousEvalVEC_[i].Resize(dim_);
	    YpreviousEvalVEC_[i].Init();

      XpreviousItVEC_[i].Resize(dim_);
      XpreviousItVEC_[i].Init();

      YpreviousItVEC_[i].Resize(dim_);
      YpreviousItVEC_[i].Init();

      XnextToLastTSVEC_[i].Resize(dim_);
      XnextToLastTSVEC_[i].Init();

      YnextToLastTSVEC_[i].Resize(dim_);
      YnextToLastTSVEC_[i].Init();

      dXpreviousItVEC_[i].Resize(dim_);
      dXpreviousItVEC_[i].Init();

      dYpreviousItVEC_[i].Resize(dim_);
      dYpreviousItVEC_[i].Init();

      XcurrentItVEC_[i].Resize(dim_);
      XcurrentItVEC_[i].Init();

      YcurrentItVEC_[i].Resize(dim_);
      YcurrentItVEC_[i].Init();
	  }
	}

	//get initial permittivity tensor
	//we just take the first element
	it.Begin();
  const Elem * el = it.GetElem();
	LocPoint lp = Elem::shapes[el->type].midPointCoord;
	LocPointMapped lpm;
	shared_ptr<ElemShapeMap> esm = it.GetGrid()->GetElemShapeMap(el, true);
	lpm.Set(lp, esm, 0.0);

	// to calculate differential materialproperties, we need to know e0 / nu0
  if(material_->GetMaterialDatabaseName() == "Electrostatic"){
    rev_mat_fac_ = 8.854187817e-12; //eps0
    PDEName_ = "Electrostatic";
  } else if(material_->GetMaterialDatabaseName() == "Electromagnetics"){
    rev_mat_fac_ = 795774.7155; //nu0
    PDEName_ = "Electromagnetics";
  } else {
    EXCEPTION("Currently only Electrostatics and Electromagnetics are supported");
  }


  PtrCoefFct matCoef = material_->GetTensorCoefFnc(matType_, tensorType_,
                                               Global::REAL, false);
  /*
   * matInitialTensor_ is the tensor to be returned, when
   * returnInitialTensor_ is true
   */
  matCoef->GetTensor(matInitialTensor_, lpm);

  /*
   * matFreeFieldTensor_ is the tensor to be returned, when
   * returnFreeFieldTensor_ is true
   */
  matFreeFieldTensor_ = Matrix<Double>(dim_,dim_);
  matFreeFieldTensor_.Init();
  for(UInt i = 0; i < dim_; i++){
    matFreeFieldTensor_[i][i] = rev_mat_fac_;
  }

  matDeltaTensor_ = new Matrix<Double>[numElemSD];
	for(UInt k = 0; k < numElemSD; k++){
	  matDeltaTensor_[k] = matFreeFieldTensor_;
	}
}


std::string CoefFunctionHyst::ToString() const {
  return "CoefFunctionHyst";
}


void CoefFunctionHyst::SetPreviousHystVals(bool setNextToLastTS_) {
/*
 * This function evaluates the hysteresis operator with the CURRENT solution
 * as input (obtained directly from equation system); the CURRENT solution X
 * and the resulting output Y of the hysteresis operator are stored as
 *      XpreviousIt_[idx] / XpreviousItVEC_[idx]
 *      YpreviousIt_[idx] / YpreviousItVEC_[idx]
 * if setNextToLastTS_ = false
 * or as
 *      XnextToLastTS_[idx] / XnextToLastTSVEC_[idx]
 *      YnextToLastTS_[idx] / YnextToLastTSVEC_[idx]
 * if setNextToLastTS_ = true
 *
 * Usage:
 * 1. this function shall be called AFTER setting up the equation system
 *        (as this might require the old values which are stored in
 *        XnextToLastTS_, YnextToLastTS_) but BEFORE the new solution of the
 *        equation system is computed, i.e. before calling algsys->solve
 * 2. during the FIRST iteration of a timestep, this function shall be called
 *        twice, once with setNextToLastTS_ = true and with the flag set to
 *        false;
 * 3. during later iteration of a timestep, this function shall only be called
 *        with setNextToLastTS_ = false
 */

  if(hystMemoryLocked_){
    WARN("SetPreviousHystVals was called while Hysteresis operator is locked. This should not be the case.")
    return;
  }

  if (returnZeroValues_){
    WARN("SetPreviousHystVals was called with returnZeroValues_ set to true! This should not be the case.")
    return;
  }

  EntityIterator it = SDList_->GetIterator();

  //Vector<Double> elemSol;
  for(it.Begin(); !it.IsEnd(); it++) {
    const Elem * el = it.GetElem();
    LocPoint lp = Elem::shapes[el->type].midPointCoord;
    LocPointMapped lpm;
    shared_ptr<ElemShapeMap> esm = it.GetGrid()->GetElemShapeMap(el, true);
    lpm.Set(lp, esm, 0.0);
    UInt idx = globalElem2Local_[el->elemNum];

    if(methodType_ == SCALAR){
      Double X, Y;
      ComputeXY(lpm, X, Y);
      if(setNextToLastTS_){
        XnextToLastTS_[idx] = X;
        YnextToLastTS_[idx] = Y;
      } else {
        XpreviousIt_[idx] = X;
        YpreviousIt_[idx] = Y;
      }

    } else if(methodType_ == VECTOR){
      Vector<Double> X,Y;
      ComputeXY_vec(lpm, X, Y);

      if(setNextToLastTS_){
        //std::cout << "Set previous TS values" << std::endl;
        XnextToLastTSVEC_[idx] = X;
        YnextToLastTSVEC_[idx] = Y;
      } else {
        //std::cout << "Set previous IT values" << std::endl;
        XpreviousItVEC_[idx] = X;
        YpreviousItVEC_[idx] = Y;
      }
    }
  }
}

void CoefFunctionHyst::SetFlag(std::string flagName,bool newState){

  if(flagName == "returnFreeFieldTensor"){
    /*
     * flagName == "returnFreeFieldTensor"
     * on: getTensor will return rev_mat_fac (eps0,nu0);
     *      will set flag returnInitialTensor_ to Off/false
     * off: getTensor will return rev_mat_fac (eps0,nu0) + deltaMat
     */
    returnFreeFieldTensor_ = newState;
    if(newState == true){
      returnInitialTensor_ = false;
    }
  } else if (flagName == "returnInitialTensor"){
    /*
     * flagName == "returnInitialTensor"
     * on: getTensor will return a predefined initial value for eps/nu,
     *      e.g. eps at working point E = 0 on the major loop or 10*eps0
     *        -> an appropriate value has to be found out
     *     will set flag returnFreeFieldTensor_ to Off/false
     * off: getTensor will return rev_mat_fac (eps0,nu0) + deltaMat
     */
    returnInitialTensor_ = newState;
    if(newState == true){
      returnFreeFieldTensor_ = false;
    }
  } else if (flagName == "addFreeFieldTensorToDeltaMat"){
    /*
     * flagName == "returnFreeFieldTensorToDeltaMat"
     * on: getTensor will return rev_mat_fac (eps0,nu0) + deltaMat
     * off: getTensor will return delatMat + initial value for eps/nu
     */
    addFreeFieldTensorToDeltaMat_ = newState;
    addInitialTensorToDeltaMat_ = !newState;
  } else if (flagName == "addInitialTensorToDeltaMat"){
    /*
     * flagName == "addInitalValuesToDeltaMat"
     * on: getTensor will return delatMat + initial value for eps/nu
     * off: getTensor will return rev_mat_fac (eps0,nu0) + deltaMat
     */
    addInitialTensorToDeltaMat_ = newState;
    addFreeFieldTensorToDeltaMat_ = !newState;

  } else if (flagName == "useNextToLastTS"){
    /*
     * flagName == "useNextToLastTS"
     * on: old values Y_old and X_old used in deltaComputation
     *      (see flags useDeltaY and useDeltaX) will be the
     *      ones of the last timestep rather than from the last
     *      iteration. I.e.
     *       Y_old = YnextToLastTS_ , YnextToLastTSVEC_
     *       X_old = XnextToLastTS_ , XnextToLastTSVEC_
     * off:
     *     old values of Y_old and X_old will be taken from
     *     last iteration
     *       Y_old = YpreviousIt_ , YpreviousItVEC_
     *       X_old = XpreviousIt_ , XpreviousItVEC_
     */
    useNextToLastTS_ = newState;
  } else if (flagName == "useDeltaYInNominator"){
    /*
     * flagName == "useDeltaYInNominator"
     * on: deltaMat will be computed with (Y_new - Y_old)/denominator, where
     *      Y_new is the current output of the hysteresis operator and
     *      Y_old is the value of the last iteraton
     *      special case: first iteration of a timestep:
     *      Y_new = result of last timestep
     *      Y_old = result of next to last timestep (or something else?!)
     * off: deltaMat will be computed using ONLY Y_new in the nominator
     */
    useDeltaY_ = newState;
  } else if (flagName == "useDeltaXInDenominator"){
    /*
     * flagName == "useDeltaXInDenominator"
     * on: deltaMat will be computed with nominator/(X_new - X_old), where
     *      X_new is the current input to the hysteresis operator and
     *      X_old is the value of the last iteraton
     *      special case: first iteration of a timestep:
     *      X_new = input of last timestep = current solution before system
     *              gets solved anew
     *      X_old = input of next to last timestep (or something else?!)
     * off: deltaMat will be computed using ONLY X_new in the denominator
     */
    useDeltaX_ = newState;
  } else if (flagName == "lockHysteresisMemory"){
    /*
     * flagName == "lockHysteresisMemory"
     * on: memory of hysteresis operator will be locked, i.e. all changes due to
     *      input will be performed on a temporal storage (needed e.g. during
     *      linesearch where different steppings shall be tested)
     * off: input will lead to permanent changes of hysteresis operator
     */
    hystMemoryLocked_ = newState;
  } else if (flagName == "lockHysteresisDirection"){
    /*
     * flagName == "lockHysteresisDirection"
     * on: rotation state of VECTOR hysteresis is kept fixed
     * off: rotation state of VECTOR hysteresis operator is allowed to follow
     *      current input
     */
    hystDirectionLocked_ = newState;
  } else if (flagName == "returnZeroValues"){
    /*
     * flagName == "returnZeroValues"
     * on: getScalar and getTensor will only return 0; needed if hysteresis
     *      operator was put on rhs of system but a call to AssembleNonLinRhs
     *      should not return a hysteresis value; hysteresis operator will not
     *      be called at all!
     * off: getScalar and getTensor will trigger an evaluation of hysteresis
     *      operator and return the result
     */
    returnZeroValues_ = newState;
  } else if (flagName == "allowBMP"){
    /*
     * flagName == "allowBMP_"
     * on: computeXY_vec outputs the overlaid switching and rotation state of
     *      the vector Preisach model as a BMP image, but only if printOut was
     *      set to true, too (printOut is to be set via material file)
     * off: disables output of bmp, even if printOut was set to true
     */
    allowBMP_ = newState;
  }

  else {
    WARN("flagName = "<<flagName<<" not known");
    return;
  }
}


void CoefFunctionHyst::GetScalar( Double& valY,const LocPointMapped& lpm ) {

  if(methodType_ != SCALAR){
    EXCEPTION("Only implemented for scalar model");
  }

  if(returnZeroValues_){
    valY = 0;
    return;
  }

  Double valX;
  /*
   * ComputeXY checks if values were already used; in that case, no recomputation is done
   */
  ComputeXY(lpm, valX, valY);
}

void CoefFunctionHyst::GetVector( Vector<Double>& valY,const LocPointMapped& lpm) {

  if(returnZeroValues_){
    valY.Resize(dependCoef_->GetVecSize());
    valY.Init();
    return;
  }

  if(methodType_ == SCALAR){
    Double valX_scal;
    Double valY_scal;
    /*
     * ComputeXY checks if values were already used; in that case, no recomputation is done
     */
    ComputeXY(lpm, valX_scal, valY_scal);

    valY.Resize(dependCoef_->GetVecSize());
    valY.Init();
    valY[dirP_] = valY_scal;
  }
  else if(methodType_ == VECTOR){
    Vector<Double> valX;
    /*
     * ComputeXY checks if values were already used; in that case, no recomputation is done
     */
    ComputeXY_vec(lpm, valX, valY);

    std::cout << "GetVector: " << std::endl;
    std::cout << std::setprecision(16) << std::scientific << "E_in:  " << valX.ToString() << std::endl;
    std::cout << std::setprecision(16) << std::scientific << "P_out:  " << valY.ToString() << std::endl;

    for(UInt i = 0; i < dim_; i++){
      valX[i] += 4.0;
      valY[i] += 4.0;
    }

    for(UInt i = 0; i < dim_; i++){
      valX[i] -= 4.0;
      valY[i] -= 4.0;
    }

    std::cout << "GetVector (mod): " << std::endl;
    std::cout << std::setprecision(16) << std::scientific << "E_in:  " << valX.ToString() << std::endl;
    std::cout << std::setprecision(16) << std::scientific << "P_out:  " << valY.ToString() << std::endl;

  } else {
    EXCEPTION("methodType_ " << methodType_ << " unknown")
  }
}


void CoefFunctionHyst::GetTensor( Matrix<Double>& tensor,
                                  const LocPointMapped& lpm ) {

  if(returnFreeFieldTensor_){
    tensor = matFreeFieldTensor_;
    return;
  }

  if(returnInitialTensor_){
    tensor = matInitialTensor_;
    return;
  }

  /*
   * calculate deltaMatrix
   */
  const Elem * el = lpm.ptEl;

  UInt idx = globalElem2Local_[el->elemNum];

  if(methodType_ == SCALAR){

    Double Xcurrent, Ycurrent;
    Double Xold, Yold;
    Double dX,dY;

    ComputeXY(lpm, Xcurrent, Ycurrent);

    if(useNextToLastTS_){
      Xold = XnextToLastTS_[idx];
      Yold = YnextToLastTS_[idx];
    } else {
      Xold = XpreviousIt_[idx];
      Yold = YpreviousIt_[idx];
    }

    if(useDeltaX_){
      dX = Xcurrent - Xold;
    } else {
      dX = Xcurrent;
    }

    if(useDeltaY_){
      dY = Ycurrent - Yold;
    } else {
      dY = Ycurrent;
    }

   // std::cout << "dX, dY " << dX << " " << dY << std::endl;

    Vector<Double> dX_vec = Vector<Double>(dim_);
    dX_vec.Init();
    dX_vec[dirP_] = dX;

    Vector<Double> dY_vec = Vector<Double>(dim_);
    dY_vec.Init();
    dY_vec[dirP_] = dY;

    CreateDeltaMatrix(dX_vec,dY_vec,matDeltaTensor_[idx]);

    tensor = matDeltaTensor_[idx];

  } else if(methodType_ == VECTOR ){

    /*
     * New approach started Jan. 25th 2017
     *
     * Idea: Calculate deltaMaterial using a not yet fixed method (e.g. solving a least square relation between
     *        vector deltaY and deltaX) but use this deltaMaterial ONLY for the system matrix
     *       The right hand side will be updated using NOT
     *          deltaMat * X_new
     *       but instead
     *          rev_mat_factor * X_new + Y_new
     *       In the later case, this function only returns a
     *          rev_mat_factor * eyes(dim_)
     *       -> Y_new has to be put onto the rhs in the corresponding PDE!
     * Idea behind idea:
     *       On rhs we are actually not allowed to use deltaMat as material relation.
     *       We could use Y_new/X_new instead, but that would fail at remanence.
     *       -> for detailed description of problem please see CoefFunctionHyst.hh
     */

    Vector<Double> Xcurrent, Ycurrent;
    Vector<Double> Xold, Yold;
    Vector<Double> dX_vec,dY_vec;

    ComputeXY_vec(lpm, Xcurrent, Ycurrent);

    if(useNextToLastTS_){
      Xold = XnextToLastTSVEC_[idx];
      Yold = YnextToLastTSVEC_[idx];
    } else {
      Xold = XpreviousItVEC_[idx];
      Yold = YpreviousItVEC_[idx];
    }

    dX_vec = Xcurrent;
    if(useDeltaX_){
      dX_vec.Add(-1.0,Xold);
    }

    dY_vec = Ycurrent;
    if(useDeltaY_){
      dY_vec.Add(-1.0,Yold);
    }

    CreateDeltaMatrix(dX_vec,dY_vec,matDeltaTensor_[idx]);

    tensor = matDeltaTensor_[idx];

//     /*
//      * check if the change changed
//      */
//     Vector<Double> tmpX = dX;
//     Vector<Double> tmpY = dY;
//     tmpX -= dXpreviousItVEC_[idx];
//     tmpY -= dYpreviousItVEC_[idx];
//
//     if((dX.NormL2()/Xcurrent.NormL2() < -1e-17) || dY.NormL2()/Ycurrent.NormL2() < -1e-17){
//       //std::cout << "Reuse old value as rel norm is too small" << std::endl;
//       /*
//        * reuse old matrix
//        */
//       tensor = matDeltaTensor_[idx];
//     } else {
//       if((tmpX.NormL2() < -tol_) && (tmpY.NormL2() < -tol_)){
//         //std::cout << "Reuse old value as change to last dX/dY is too small" << std::endl;
//         /*
//          * reuse old matrix
//          */
//         tensor = matDeltaTensor_[idx];
//       } else {
//         //std::cout << "Compute new value" << std::endl;
//         /*
//          * calculate new matrix
//          */
//         CreateDeltaMatrix(dX,dY,matDeltaTensor_[idx],idx);
//
//         /*
//          * store last values
//          */
//         //dXpreviousItVEC_[idx] = dX;
//         //dYpreviousItVEC_[idx] = dY;
//         tensor = matDeltaTensor_[idx];
//       }
//     }
   }

  //std::cout << "Returning tensor: " << std::endl;
  //std::cout << tensor.ToString() << std::endl;
}


void CoefFunctionHyst::ComputeXY( const LocPointMapped& lpm, Double& X, Double& Y) {

  // evaluate vector of dependency
  Vector<Double> elemSol;

  //we just need the value in the element mid point!
  const Elem * el = lpm.ptEl;
  LocPoint lp = Elem::shapes[el->type].midPointCoord;
  LocPointMapped copylpm;
  shared_ptr<ElemShapeMap> esm = lpm.shapeMap;
  copylpm.Set(lp, esm, 0.0);

  dependCoef_->GetVector(elemSol, copylpm);

  //std::cout << "ElemSol: " << elemSol.ToString() << std::endl;

  UInt idx = globalElem2Local_[el->elemNum];
  X = elemSol[dirP_];

  /*
   * Compare to input of last EVALUATION (not iteration!)
   */
  //std::cout << "X-XpreviousEval_[idx]: " << X-XpreviousEval_[idx] << std::endl;

  if( abs(X-XpreviousEval_[idx]) < -tol_) {
    /*
     * dangerous!
     * this seems to save time, but in fact it may lead to strange errors!
     * description:
     *  two input may be very close together (diff < tol) so that old value is
     *  reused; during the next iteration diff might only be slightly larger than
     *  tol so that a new value gets computed which should not lead to issues as
     *  the the new value is so close to the old one; however, the new output
     *  might have a strong jump in it
     * reason:
     *  although two input values might have roughly the same value, they might
     *  be different in terms of their maximum-minimum type
     *  if we do not update the list for the everett function in a proper way,
     *  a new input (even though it only changed slightly) might erase some
     *  values of the list, leading to jumps in the output
     * solution:
     *  do not check here, but let Preisach operator do the work!
     */
    Y = YpreviousEval_[idx];
  } else {
    /*
     * check different material cases
     * 1. electrostatics:
     *      elemSol = E -> Elem solution can be used directly as input to
     *                      hysteresis operator
     * 2. magnetics:
     *      elemSol = B -> estimate corresponding H by H = nu0*(B - J_last)
     *                      with J_last =
     *                      if(useNextToLastTS_):  YnextToLastTSVEC_[idx];
     *                      else: YpreviousItVEC_[idx];
     *
     * Note: for magnetics, we cannot use H as elemSol directly!
     * Reason: magneticPDE tries to compute H by nu*B; however, to compute
     *          nu, getTensor is called which itself calls ComputeXY; this
     *          ends in an infinite recursion
     */
    if(PDEName_ == "Electromagnetics"){
      if(useNextToLastTS_){
        X -= YnextToLastTS_[idx]; //B-J_last_iteration
      } else {
        X -= YpreviousIt_[idx]; //B-J_last_iteration
      }

      X = X*rev_mat_fac_; //(B-J_last_iteration)*nu0
    }

    /*
     * In older versions, we stored the normalized/cutted version of X to
     * avoid recalculations once X stays above saturation.
     * Normally, that does not lead to large computations as the Hysteresis operator
     * does the see the cutting and compares this value with the already entered history.
     *
     * For the computation of the deltaMaterial tensor, a cutted X leads to dX = 0 in saturation,
     * so that dY/dX becomes a problem. Of course, dY should be 0, too, and we can simply check for
     * dX = 0. However, to be consistent with the vector case, we store the uncutted input.
     */
    XpreviousEval_[idx] = X;

    /*
     * NOTE:
     * for magnetics: X_in = H, X = B, Y = J = mu0 * M
     * for electrostatics: X_in = E, X = E, Y = P
     */
    Y = hyst_->computeValueAndUpdate(X, idx, !hystMemoryLocked_);

    YpreviousEval_[idx] = Y;
  }

  std::cout << " <- in " << X << std::endl;
  std::cout << " -> out " << Y << std::endl;
}

void CoefFunctionHyst::ComputeXY_vec( const LocPointMapped& lpm, Vector<Double>& X, Vector<Double>& Y) {

  // evaluate vector of dependency
  Vector<Double> elemSol;

  //we just need the value in the element mid point!
  const Elem * el = lpm.ptEl;
  LocPoint lp = Elem::shapes[el->type].midPointCoord;
  LocPointMapped copylpm;
  shared_ptr<ElemShapeMap> esm = lpm.shapeMap;
  copylpm.Set(lp, esm, 0.0);
  dependCoef_->GetVector(elemSol, copylpm);

  UInt idx = globalElem2Local_[el->elemNum];

  X = elemSol;

  Vector<Double> tmp = X;
  /*
   * Compare to last value that was evaluated at this point
   */
  tmp -= XpreviousEvalVEC_[idx];
  if(tmp.NormL2() < -tol_){
    /*
     * no change in input -> (hopefully) no change in output; in case of hysteresis it is possible
     *                          that the same input leads to a different output due to the internal state
     *                          but as the operator was evaluated with the same value the last time it was called,
     *                          it cannot have changed its memory
     */
    Y = YpreviousEvalVEC_[idx];
  } else {

    /*
     * Do not normalize here
     * -> the Hysteresis operator can handle input > xSat_; it will simply drive the system to its bounds
     *    (i.e. completely fills the Preisach plane)
     *    -> Normalizing not needed
     */

//    if(X.NormL2() > xSat_){
//      //std::cout << "Normalized X to xSat" << std::endl;
//      X = X*xSat_/X.NormL2();
//    }

    /*
     * check different material cases
     * 1. electrostatics:
     *      elemSol = E -> Elem solution can be used directly as input to
     *                      hysteresis operator
     * 2. magnetics:
     *      elemSol = B -> estimate corresponding H by H = nu0*(B - J_last)
     *                      with J_last =
     *                      if(useNextToLastTS_):  YnextToLastTSVEC_[idx];
     *                      else: YpreviousItVEC_[idx];
     *
     * Note: for magnetics, we cannot use H as elemSol directly!
     * Reason: magneticPDE tries to compute H by nu*B; however, to compute
     *          nu, getTensor is called which itself calls ComputeXY; this
     *          ends in an infinite recursion
     */
    if(PDEName_ == "Electromagnetics"){
      if(useNextToLastTS_){
        X -= YnextToLastTSVEC_[idx]; //B-J_last_iteration
      } else {
        X -= YpreviousItVEC_[idx]; //B-J_last_iteration
      }

      X = X*rev_mat_fac_; //(B-J_last_iteration)*nu0
    }

    /*
     * NOTE:
     * for magnetics: X_in = H, X = B, Y = J = mu0 * M
     * for electrostatics: X_in = E, X = E, Y = P
     */
    Y = hyst_->computeValue_vec(X, idx, !hystMemoryLocked_,!hystDirectionLocked_);

    /*
     * store the value from last evaluation
     * Note: there is a difference between
     * XpreviousEvalVEC_[idx] and XpreviousItVEC_[idx]
     * The later one is only set via SetPreviousHystValues
     */
    XpreviousEvalVEC_[idx] = X;
    YpreviousEvalVEC_[idx] = Y;

    /*!
     * Print out of Hysteresis state
     * ONLY FOR DEBUGGING REASONS!
     * USE CAREFULLY! MASSIVE IMAGE OUTPUT!
     */
    static UInt cnt = 0;
    static UInt firstIdx = 1;

    if((printOut_ > 0) && (allowBMP_ == true) ){
      if((cnt%printOut_==0)&&(idx==firstIdx)){
        std::cout << "Outputting bmp" << std::endl;
        std::stringstream filenamebuf;
        filenamebuf << "Switch_Elem"<<firstIdx<<"_Step" << std::setfill('0') << std::setw(5)<<cnt<<"_v"<<evalVersion_<<"_numRows"<<numRows_<<".bmp";
        hyst_->switchingStateToBmp(bmpResolution_,filenamebuf.str(),idx,true);
      }

      std::cout << "idx: " << idx << std::endl;
      std::cout << "firstIdx: " << firstIdx << std::endl;

      if(idx == firstIdx){
        cnt++;
        allowBMP_ = false;
      }
      /*
       * disable output until reset
       * -> otherwise we would get two images for each timestep if P and D are computed
       */
    }
  }
}

void CoefFunctionHyst::CreateDeltaMatrix(Vector<Double>& dX,Vector<Double>& dY, Matrix<Double>& deltaMat){

  /*
   * Function for calculation of deltaMatrix
   *  options (to be triggered via SetFlag):
   *    useDeltaY_, useDeltaX_, useNextToLastTS_
   * Note: these parameter will determine the values of dX and dY!
   *    useDeltaX_ = true -> dX = X_current - X_old
   *               = false > dX = X_current
   *    useDeltaY_ = true -> dY = Y_current - Y_old
   *               = false > dY = Y_current
   *    useNextToLastTS_ = true -> X_old/Y_old = XnextToLastTS_/YnextToLastTS_
   *                     = false > X_old/Y_old = XpreviousIt_/YpreviousIt_
   *
   */

  std::string version_name = "DirectDivisionAbs"; // "L1_Regularized_LeastSquares";// "Unregularized_LeastSquares";
  Double scaling = 0.0;
  if(PDEName_ == "Electromagnetics"){
    /*
     * Electromagnetics (J = magnetic polarization = mu_0 * M, rev_mat_fac_ = nu0)
     *
     * useDeltaY_ = false:
     * H_n+1 = rev_mat_fac_ * (B_n+1 - J_n+1)
     *       = rev_mat_fac_ * (B_n + deltaB - J_n+1)
     *       = rev_mat_fac_ * B_n + (rev_mat_fac_ - rev_mat_fac_*J_n+1/deltaB)*deltaB
     *
     * -> deltaMatrix = rev_mat_fac_ - rev_mat_fac_*J_n+1/deltaB
     *
     * useDeltaY_ = true:
     * H_n+1 = rev_mat_fac_ * (B_n+1 - J_n+1)
     *       = rev_mat_fac_ * (B_n + deltaB - J_n - deltaJ)
     *       = rev_mat_fac_ * (B_n - J_n) + (rev_mat_fac_ - rev_mat_fac_*deltaJ/deltaB)*deltaB
     *
     * -> deltaMatrix = rev_mat_fac_ - rev_mat_fac_*deltaJ/deltaB
     *
     *
     *  -> scaling = -1 * rev_mat_fac
     */
    scaling = -rev_mat_fac_;

  } else if(PDEName_ == "Electrostatic"){
    /*
     * Electrostatic (P = electric polarization, rev_mat_fac_ = eps0)
     *
     * useDeltaY_ = false:
     * D_n+1 = rev_mat_fac_ * E_n+1 + P_n+1
     *       = rev_mat_fac_ * (E_n + deltaE) + P_n+1
     *       = rev_mat_fac_ * E_n + (rev_mat_fac_ + P_n+1/deltaE)*deltaE
     *
     * -> deltaMatrix = rev_mat_fac_ + P_n+1/deltaE
     *
     * useDeltaY_ = true:
     * HDn+1 = rev_mat_fac_ * E_n+1 + P_n+1
     *       = rev_mat_fac_ * (E_n + deltaE) + P_n + deltaP
     *       = rev_mat_fac_ * E_n + P_n + (rev_mat_fac_ + deltaP/deltaE)*deltaE
     *
     * -> deltaMatrix = rev_mat_fac_ + deltaP/deltaE
     *
     *
     *  -> scaling = +1
     */
    scaling = 1.0;

  }

  /*
   * reset deltaMatrix
   */
  deltaMat.Resize(dim_,dim_);
  deltaMat.Init();

  if(version_name == "DirectDivision"){
    /*
     * Take the most straight forward approach
     * eps_ii = dY_used_i/dX_used_i
     * eps_ij = 0 (i != j)
     *
     * Open question:
     *  what happens if dX_used_i = 0?
     *    dY_used = 0 -> 0?
     *    dY_used != 0 -> 0?
     */
    Double addValue;
    for(UInt i = 0; i < dim_; i++){
      if(addFreeFieldTensorToDeltaMat_){
        addValue = rev_mat_fac_;
      } else if(addInitialTensorToDeltaMat_){
        addValue = matInitialTensor_[i][i];
      } else {
        EXCEPTION("Either addFreeFieldTensorToDeltaMat_ or addInitialTensorToDeltaMat_ should be true");
      }

      if(dX[i] != 0){
        // TODO: should we use addValue for scaling in magnetics, too?
        //        or should we always take nu0?
        deltaMat[i][i] = addValue + scaling*(dY[i]/dX[i]);
      } else {
        //take initial values? -> open question
        deltaMat[i][i] = addValue;
      }
    }
  } else if(version_name == "DiretDivisionAbs"){
    /*
     * same as DirectDivision except for using abs(dY/dX)
     */
    Double addValue;
    for(UInt i = 0; i < dim_; i++){
      if(addFreeFieldTensorToDeltaMat_){
        addValue = rev_mat_fac_;
      } else if(addInitialTensorToDeltaMat_){
        addValue = matInitialTensor_[i][i];
      } else {
        EXCEPTION("Either addFreeFieldTensorToDeltaMat_ or addInitialTensorToDeltaMat_ should be true");
      }

      if(dX[i] != 0){
        // TODO: should we use addValue for scaling in magnetics, too?
        //        or should we always take nu0?
        deltaMat[i][i] = addValue + scaling*abs(dY[i]/dX[i]);
      } else {
        //take initial values? -> open question
        deltaMat[i][i] = addValue;
      }
    }
  } else {
    EXCEPTION("Computation of deltaMatrix via "<<version_name<<" not implemented");
  }

// OLD VERSIONS, additional computation methods, most do not work
//  Double rev_mat_fac_used;
//  Double scaling;
//  Vector<Double> dX_used;
//  Vector<Double> dY_used;
//  bool add_rev = true;
//
//  if(compute_inverse_){
//    /*
//     * Electromagnetics (J = magnetic polarization = mu_0 * M)
//     * H_n+1 = 1/rev_mat_fac_ * (B_n+1 - J_n+1)
//     *       = 1/rev_mat_fac_ * (B_n + deltaB - J_n+1)
//     *       = 1/rev_mat_fac_ * B_n + (1/rev_mat_fac_ - 1/rev_mat_fac_*J_n+1/deltaB)*deltaB
//     *
//     * -> deltaMatrix = 1/rev_mat_fac_ - 1/rev_mat_fac_*J_n+1/deltaB
//     */
//
//    rev_mat_fac_used = 1/rev_mat_fac_;
//    scaling = -rev_mat_fac_used;
//    dX_used = dX; // make sure that deltaB is passed to this function rather than deltaH
//    dY_used = dY; // -> apply scaling later on as we might use DirectDivisionAbs as method
//                  // which takes abs(dY/dX) -> here we would loose the -1 factor
//  } else {
//    /*
//     * Electrostatics:
//     * D_n+1 = rev_mat_fac_ * E_n+1 + P_n+1
//     *       = rev_mat_fac_ * E_n + rev_mat_fac_ * deltaE + P_n+1
//     *       = rev_mat_fac_ * E_n + (rev_mat_fac_ + P_n+1/deltaE)*deltaE
//     *
//     * -> deltaMatrix = rev_mat_fac + P_n+1/deltaE
//     */
//    scaling = 1;
//    rev_mat_fac_used = rev_mat_fac_;
//    dX_used = dX;
//    dY_used = dY;
//  }
//  Matrix<Double> deltaMatBak = Matrix<Double>(deltaMat);
//  /*
//   * reset deltaMatrix
//   */
//  deltaMat.Resize(dim_,dim_);
//  deltaMat.Init();
//
//  std::string version_name = "DirectDivisionAbs"; // "L1_Regularized_LeastSquares";// "Unregularized_LeastSquares";
//
//  if(version_name == "DirectDivision_OnlyMaxDirectionMod"){
//    /*
//     * Take the most straight forward approach
//     * eps_ii = dY_used_i/dX_used_i
//     * eps_ij = 0 (i != j)
//     * but apply it only to the maximal component
//     */
//    Double dy_max = 0.0;
//    Double tmp = 0.0;
//    UInt idx = 0;
//    for(UInt i = 0; i < dim_; i++){
//      if(dX_used[i] != 0){
//        tmp = abs(dY_used[i]);
//      }
//      if(tmp > dy_max){
//        dy_max = tmp;
//        idx = i;l
//      }
//    }
//
//    for(UInt i = 0; i < dim_; i++){
//      if( (dX_used[i] != 0) && (abs(dY_used[i]/dy_max)>1e-2) ){
//        deltaMat[i][i] = dY_used[i]/dX_used[i];
//      }
//    }
//  }
//
//  else if(version_name == "DirectDivision_OnlyMaxDirection"){
//    /*
//     * Take the most straight forward approach
//     * eps_ii = dY_used_i/dX_used_i
//     * eps_ij = 0 (i != j)
//     * but apply it only to the maximal component
//     */
//    Double dy_max = 0.0;
//    Double tmp = 0.0;
//    UInt idx = 0;
//    for(UInt i = 0; i < dim_; i++){
//      if(dX_used[i] != 0){
//        tmp = abs(dY_used[i]);
//      }
//      if(tmp > dy_max){
//        dy_max = tmp;
//        idx = i;
//      }
//    }
//
//    if(dX_used[idx] != 0){
//      deltaMat[idx][idx] = dY_used[idx]/dX_used[idx];
//    }
//
//  }
//  else if(version_name == "DirectDivision"){
//        add_rev = false;
//        /*
//         * Take the most straight forward approach
//         * eps_ii = dY_used_i/dX_used_i
//         * eps_ij = 0 (i != j)
//         *
//         * Open question:
//         *  what happens if dX_used_i = 0?
//         *    dY_used = 0 -> 0?
//         *    dY_used != 0 -> 0?
//         */
//        for(UInt i = 0; i < dim_; i++){
//          if(dX_used[i] != 0){
//            std::cout << "dX_used[i] = " << dX_used[i] << std::endl;
//            std::cout << "dY_used[i] = " << dY_used[i] << std::endl;
//            deltaMat[i][i] = rev_mat_fac_used + scaling*(dY_used[i]/dX_used[i]);
//            //deltaMat[i][i] = scaling*abs(dY_used[i]/dX_used[i]);
//          } else {
//            std::cout << "dX_used[i] = 0" << std::endl;
//            std::cout << "dY_used[i] = " << dY_used[i] << std::endl;
//            std::cout << "Reuse old value as dX == 0" << std::endl;
//            deltaMat[i][i] = deltaMatBak[i][i];
//          }
//        }
//      }
//  else if(version_name == "DirectDivisionAbs"){
//      add_rev = false;
//      /*
//       * Take the most straight forward approach
//       * eps_ii = dY_used_i/dX_used_i
//       * eps_ij = 0 (i != j)
//       *
//       * Open question:
//       *  what happens if dX_used_i = 0?
//       *    dY_used = 0 -> 0?
//       *    dY_used != 0 -> 0?
//       */
//      for(UInt i = 0; i < dim_; i++){
//        if(dX_used[i] != 0){
//          std::cout << "dX_used[i] = " << dX_used[i] << std::endl;
//          std::cout << "dY_used[i] = " << dY_used[i] << std::endl;
//          deltaMat[i][i] = rev_mat_fac_used + scaling*abs(dY_used[i]/dX_used[i]);
//          //deltaMat[i][i] = scaling*abs(dY_used[i]/dX_used[i]);
//        } else {
//          std::cout << "dX_used[i] = 0" << std::endl;
//          std::cout << "dY_used[i] = " << dY_used[i] << std::endl;
//          std::cout << "Reuse old value as dX == 0" << std::endl;
//          deltaMat[i][i] = deltaMatBak[i][i];
//        }
//      }
//    }
//  else if(version_name == "DirectDivisionMod"){
//    /*
//     * Take the most straight forward approach
//     * eps_ii = dY_used_i/dX_used_i
//     * eps_ij = 0 (i != j)
//     *
//     * Open question:
//     *  what happens if dX_used_i = 0?
//     *    dY_used = 0 -> 0?
//     *    dY_used != 0 -> 0?
//     */
//    for(UInt i = 0; i < dim_; i++){
//      if( (dX_used[i] != 0) && (abs(dY_used[i]) > 1e-15) ) {
//        std::cout << "dX_used[i] = " << dX_used[i] << std::endl;
//        std::cout << "dY_used[i] = " << dY_used[i] << std::endl;
//      //if( (dX_used[i] != 0) && (abs(dY_used[i]) > 1e-16) ) {
//        deltaMat[i][i] = dY_used[i]/dX_used[i];
//      } else {
//
//        std::cout << "dX_used[i] = " << dX_used[i] << std::endl;
//        std::cout << "dY_used[i] = " << dY_used[i] << std::endl;
//      }
//    }
//  }
//  else if(version_name == "DirectDivisionMod2"){
//    /*
//     * Take the most straight forward approach
//     * eps_ii = dY_used_i/dX_used_i
//     * eps_ij = 0 (i != j)
//     *
//     * Open question:
//     *  what happens if dX_used_i = 0?
//     *    dY_used = 0 -> 0?
//     *    dY_used != 0 -> 0?
//     */
//    for(UInt i = 0; i < dim_; i++){
//      /*
//       * set tolerance for dX such large, that in the worst case of only rounding
//       * errors in dY, i.e. dY = tol_y, the resulting eps will not be orders of
//       * magnitude larger than eps
//       * -> this happens especially in the beginning of very simple setups
//       * (only 1 fix direction)
//       * -> due to rounding errors in dY, we suddenly get 1e-20/1e-12=1e-8 >>> eps0 = 1e-12
//       * -> even worse, this is different for each element (some might habve 1e-20, others
//       * 2e-19 or 3e-18
//       * -> huge differences in material tensor
//       * -> upcoming iterations will magnify this problem as they use this warped
//       * material tensor as basis
//       */
//      Double tol_y = 1e-15;
//      if( (abs(dX_used[i]) > tol_y/rev_mat_fac_used) && (abs(dY_used[i]) > tol_y) ) {
//      //if( (dX_used[i] != 0) && (abs(dY_used[i]) > 1e-16) ) {
//        deltaMat[i][i] = dY_used[i]/dX_used[i];
//      } else {
//        std::cout << "dX_used[i] = 0" << std::endl;
//        std::cout << "dY_used[i] = " << dY_used[i] << std::endl;
//      }
//    }
//  }
//  else if(version_name == "L1_Regularized_LeastSquares"){
//    /*
//     * Perform the shooting algorithm to solve least squares problem
//     *  y = Ax
//     * with a L1 regularization, i.e. norm_L1(x) < lambda
//     * with norm_L1(x) = sum_i abs(x_i)
//     *
//     * Definition:
//     * 2D case:
//     *
//     *  deltaMat = [ [d11 d12]; [d21 d22] ] = [ [d11 d12]; [d12 d22] ]
//     *  dX = [dX1 dX2]
//     *  dY = [dY1 dY2]
//     *
//     *  x = [d11 d22 d12] (Voigt)
//     *  A = [[dX1 0 dX2]
//     *       [0 dX2 dX1]]
//     *  AT = [[dX1   0]
//     *        [  0 dX2]
//     *        [dX2 dX1]]
//     *  y = [dY1 dY2]
//     *
//     * 3D case:
//     *
//     *  deltaMat = [ [d11 d12 d13]; [d21 d22 d23]; [d31 d32 d33] ] = [ [d11 d12 d13]; [d12 d22 d23]; [d13 d23 d33] ]
//     *  dX = [dX1 dX2 dX3]
//     *  dY = [dY1 dY2 dY3]
//     *
//     * x = [d11 d22 d33 d23 d13 d12] (Voigt)
//     * A = [[dX1   0   0   0 dX3 dX2]
//     *      [  0 dX2   0 dX3   0 dX1]
//     *      [  0   0 dX3 dX2 dX1   0]]
//     * AT = [[dX1   0   0]
//     *       [  0 dX2   0]
//     *       [  0   0 dX3]
//     *       [  0 dX3 dX2]
//     *       [dX3   0 dX1]
//     *       [dX2 dX1   0]]
//     * y = [dY1 dY2 dY3]
//     *
//     *
//     * Shooting algorithm:
//     *
//     *  start with initial guess for x
//     *  precompute ATA = A^T*A
//     *  set lambda
//     *  S = Vector<Double>(length(x))
//     *
//     *  while( norm_L2(y-Ax) > tol )
//     *    y = Ax
//     *    S = -A^T*y + ATA*x
//     *    for(j = 0; j < length(x); j++)
//     *      S(j) -= ATA(j,j)*x(j)
//     *
//     *    for(j = 0; j < length(x); j++)
//     *      if S(J) > lambda
//     *        x(j) = (lambda - S(J))/ATA(j,j)
//     *      else if S(J) < -lambda
//     *        x(j) = (-lambda - S(J))/ATA(j,j)
//     *      else
//     *        x(j) = 0.0
//     */
//
//
//    Matrix<Double> AT;
//    UInt dim2 = 0;
//
//    if(dim_ == 2){
//      // number of unknowns of the material tensor
//      dim2 = 3;
//      /*
//       * start with transposed matrix (such that we can use MultT later on
//       *
//       *  x = [d11 d22 d21] (Voigt)
//       *  A = [[dX1 0   dX2]
//       *       [  0 dX2 dX1]]
//       *  AT = [[dX1   0]
//       *        [  0 dX2]
//       *        [dX2 dX1]]
//       *  y = [dY1 dY2]
//       */
//      AT = Matrix<Double>(dim2,dim_);
//      AT.Init();
//      AT[0][0] = dX_used[0];
//      AT[1][1] = dX_used[1];
//      AT[2][0] = dX_used[1];
//      AT[2][1] = dX_used[0];
//
//    } else if (dim_ == 3) {
//      // number of unknowns of the material tensor
//      dim2 = 6;
//      /*
//       * x = [d11 d22 d33 d32 d31 d21] (Voigt)
//       * A = [[dX1   0   0   0 dX3 dX2]
//       *      [  0 dX2   0 dX3   0 dX1]
//       *      [  0   0 dX3 dX2 dX1   0]]
//       * AT = [[dX1   0   0]
//       *       [  0 dX2   0]
//       *       [  0   0 dX3]
//       *       [  0 dX3 dX2]
//       *       [dX3   0 dX1]
//       *       [dX2 dX1   0]]
//       * y = [dY1 dY2 dY3]
//       */
//      AT = Matrix<Double>(dim2,dim_);
//      AT.Init();
//      AT[0][0] = dX_used[0];
//      AT[1][1] = dX_used[1];
//      AT[2][2] = dX_used[2];
//      AT[3][1] = dX_used[2];
//      AT[3][2] = dX_used[1];
//      AT[4][0] = dX_used[2];
//      AT[4][2] = dX_used[0];
//      AT[5][0] = dX_used[1];
//      AT[5][1] = dX_used[0];
//    }
//
//    Matrix<Double> ATA = Matrix<Double>(dim2,dim2);
//    Vector<Double> x = Vector<Double>(dim2);
//    Vector<Double> S = Vector<Double>(dim2);
//    Vector<Double> ATy = Vector<Double>(dim2);
//    Vector<Double> Ax = Vector<Double>(dim_);
//    Vector<Double> res = Vector<Double>(dim2);
//    Matrix<Double> A;
//    AT.Transpose(A);
//
//    AT.Mult(A,ATA);
//
//    for(UInt j = 0; j < dim_; j++){
//      /*
//       * take rev_mat_fac_used*Identtiy as initial guess for delta matrix
//       * i.e. set the first dim_ entries of x to rev_mat_fac_used
//       */
//      x[j] = rev_mat_fac_used;
//    }
//
//    Double lambda = 1e-12;
//    Double tol = 1e-14;
//    UInt maxIters = 2;
//
//    /*
//     * compute initial residual
//     */
//    Double x_l1Norm = 0.0;
//    Double res_Norm = 0.0;
//    for(UInt i = 0; i < x.GetSize(); i++){
//      x_l1Norm += abs(x[i]);
//    }
//    res = -dY_used;
//    A.Mult(x,Ax);
//    res += Ax;
//    res_Norm = res.NormL2() + lambda*x_l1Norm;
//
//    /*
//     * precompute ATy
//     */
//    AT.Mult(dY_used,ATy);
//
//    UInt cnt = 0;
//    while(res_Norm > tol){
//      ATA.Mult(x,S);
//      S -= ATy;
//      for(UInt j = 0; j < dim2; j++){
//        S[j] -= ATA[j][j]*x[j];
//
//        Double tmp;
//        if(S[j] > lambda ){
//          tmp = ( lambda - S[j] )/ATA[j][j];
//        } else if(S[j] < -lambda ){
//          tmp = ( -lambda - S[j] )/ATA[j][j];
//        } else {
//          tmp = 0.0;
//        }
//        x[j] = tmp; //0.1*x[j]+0.9*tmp;
//      }
//
//      cnt++;
//      std::cout << "after iteration " << cnt << std::endl;
//
//      /*
//       * calc residual
//       */
//      x_l1Norm = 0.0;
//      for(UInt i = 0; i < x.GetSize(); i++){
//        x_l1Norm += abs(x[i]);
//      }
//      std::cout << "x: " << x.ToString() << std::endl;
//
//      res = -dY_used;
//      A.Mult(x,Ax);
//      res += Ax;
//
//      std::cout << "y: " << dY_used.ToString() << std::endl;
//      std::cout << "A: " << A.ToString() << std::endl;
//      std::cout << "Ax: " << Ax.ToString() << std::endl;
//      std::cout << "y-Ax: " << res.ToString() << std::endl;
//      res_Norm = res.NormL2() + lambda*x_l1Norm;
//
//      std::cout << "|y-Ax|: " << res.NormL2() << std::endl;
//      std::cout << "x_l1Norm: " << x_l1Norm << std::endl;
//      std::cout << "lambda*x_l1Norm: " << x_l1Norm*lambda << std::endl;
//      std::cout << "Res norm: " << res_Norm << std::endl;
//      if(cnt > maxIters){
//        break;
//      }
//    }
//
//    // sort vector into matrix
//    for(UInt i = 0; i < dim_; i++){
//      for(UInt j = 0; j <= i; j++){
//        if(i == j){
//          deltaMat[i][i] = x[i];
//        } else {
//          deltaMat[i][j] = x[dim2-i-j];
//          deltaMat[j][i] = deltaMat[i][j];
//        }
//      }
//    }
//
//  }
//
//  else if(version_name == "Unregularized_LeastSquares"){
//
//    /*
//     * Get vector norms
//     */
//    Double norm_dX = dX_used.NormL2();
//
//    /*
//     * If norm_dX < tol_, it is a bad idea to compute dY/dX
//     * -> in that case, norm_dY should be quite small also
//     * (can be checked); however, return simply matrix
//     * with rev_mat_fac_ on diagonal
//     */
//    if(norm_dX > tol_){
//      //std::cout << "NormDX = " << norm_dX << " is larger than tol" << std::endl;
//      std::cout << "dX: " << dX_used << std::endl;
//      std::cout << "dY: " << dY_used << std::endl;
//      /*
//       * Solve least square problem to determine delta matrix [deltaMat], such that
//       * dY = [deltaMat]*dX
//       *
//       * Assumptions: deltaMat = symmetric
//       * Procedure (for 2D, uncoupled case):l
//       *  1. definition
//       *  deltaMat = [ [d11 d12]; [d21 d22] ] = [ [d11 d12]; [d12 d22] ]
//       *  dX = [dX1 dX2]
//       *  dY = [dY1 dY2]
//       *
//       *  x = [d11 d22 d12] (Voigt)
//       *  A = [[dX1 0 dX2]
//       *       [0 dX2 dX1]]
//       *  y = [dY1 dY2]
//       *
//       *  2. solve least squares problem
//       *    Ax = y
//       *  by solving
//       *    (A*A^T)z = y
//       *  and setting
//       *    x = (A^T)z
//       *
//       *  3. create deltaMat from x
//       *
//       *  Note: solution may lead to very small diagonal entries for deltamat
//       *        in that case consider a regularization, i.e. solve
//       *          (omega*I + A*A^T)z = y
//       *        with omega being a regularization parameter
//       *
//       */
//
//      Matrix<Double> AT;
//      UInt dim2 = 0;
//
//      if(dim_ == 2){
//        // number of unknowns of the material tensor
//        dim2 = 3;
//        /*
//         * start with transposed matrix (such that we can use MultT later on
//         *
//         *  x = [d11 d22 d21] (Voigt)
//         *  A = [[dX1 0   dX2]
//         *       [  0 dX2 dX1]]
//         *  AT = [[dX1   0]
//         *        [  0 dX2]
//         *        [dX2 dX1]]
//         *  y = [dY1 dY2]
//         */
//        AT = Matrix<Double>(dim2,dim_);
//        AT.Init();
//        AT[0][0] = dX_used[0];
//        AT[1][1] = dX_used[1];
//        AT[2][0] = dX_used[1];
//        AT[2][1] = dX_used[0];
//
//      } else if (dim_ == 3) {
//        // number of unknowns of the material tensor
//        dim2 = 6;
//        /*
//         * x = [d11 d22 d33 d32 d31 d21] (Voigt)
//         * A = [[dX1   0   0   0 dX3 dX2]
//         *      [  0 dX2   0 dX3   0 dX1]
//         *      [  0   0 dX3 dX2 dX1   0]]
//         * AT = [[dX1   0   0]
//         *       [  0 dX2   0]
//         *       [  0   0 dX3]
//         *       [  0 dX3 dX2]
//         *       [dX3   0 dX1]
//         *       [dX2 dX1   0]]
//         * y = [dY1 dY2 dY3]
//         */
//        AT = Matrix<Double>(dim2,dim_);
//        AT.Init();
//        AT[0][0] = dX_used[0];
//        AT[1][1] = dX_used[1];
//        AT[2][2] = dX_used[2];
//        AT[3][1] = dX_used[2];
//        AT[3][2] = dX_used[1];
//        AT[4][0] = dX_used[2];
//        AT[4][2] = dX_used[0];
//        AT[5][0] = dX_used[1];
//        AT[5][1] = dX_used[0];
//      }
//
//      Matrix<Double> AAT = Matrix<Double>(dim_,dim_);
//      /*
//       * perform AAT = (AT)^T * AT
//       */
//      AT.MultT(AT,AAT);
//
//      Matrix<Double> AATI = Matrix<Double>(dim_,dim_);
//
//      // add regualrization
//  //      double omega = 0.0;
//  //      for(UInt i = 0; i < dim_; i++){
//  //        AAT[i][i] += omega;
//  //      }
//
//      AAT.Invert(AATI);
//
//      Vector<Double> z = Vector<Double>(dim_);
//      Vector<Double> x = Vector<Double>(dim2);
//
//      // solve least squares problem
//      z = AATI*dY_used;
//
//      // recover vector of unknowns
//      x = AT*z;
//
//      // sort vector into matrix
//      for(UInt i = 0; i < dim_; i++){
//        for(UInt j = 0; j <= i; j++){
//          if(i == j){
//            deltaMat[i][i] = x[i];
//          } else {
//            deltaMat[i][j] = x[dim2-i-j];
//            deltaMat[j][i] = deltaMat[i][j];
//          }
//        }
//      }
//    } else {
//      std::cout << "NormDX = " << norm_dX << " is smaller than tol" << std::endl;
//      //norm_dX < tol_ -> only return scaled unit matrix
//      // -> in the inverse case, where dX_used represents the difference of the Hysteresis output, norm_dX becomes
//      // 0 once the material is in saturation
//    }
//  }
//
////  std::cout << "deltaMat before adding std term: " << std::endl;
////  std::cout << deltaMat.ToString() << std::endl;
//
//  // add rev_mat_fac_used after computation of matDeltaTensorl
//  //deltaMat += matFreeFieldTensor_;
//
////  /*
////   * apply scaling
////   */
////  deltaMat = deltaMat*scaling;
////
//  /*
//   * add rev_mat_fac
//   */
////  for(UInt i = 0; i < dim_; i++){
////    deltaMat[i][i] += rev_mat_fac_used;
////  }
//
//  if(add_rev){
//  for(UInt i = 0; i < dim_; i++){
//    if(deltaMat[i][i] < 0){
//      std::cout << "Ever needed?" << std::endl;
//      deltaMat[i][i] -= rev_mat_fac_used;
//    } else {
//      deltaMat[i][i] += rev_mat_fac_used;
//    }
//  }
//  }
}

}// namespace
