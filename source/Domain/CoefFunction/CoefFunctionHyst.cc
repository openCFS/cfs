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

  // if set to false; changes to hysteresis operator will not be stored
  overwrite_ = true;
  overwriteDirection_ = true;
  deltaComputation_ = true;
  useNextToLastTS_ = true;

  Init(material, actSDList);
}

CoefFunctionHyst::~CoefFunctionHyst(){
  ;
}

void CoefFunctionHyst::Init( BaseMaterial* const material,
                             shared_ptr<ElemList> actSDList) {

	// set the parameters
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

  //Double inf = std::numeric_limits<double>::infinity();

	//allocate memory for previous results, needed for the
	//effective material parameter formulation
	if(methodType_ == SCALAR){
    XpreviousEval_.Resize(numElemSD);
    YpreviousEval_.Resize(numElemSD);
    XpreviousEval_.Init();
    YpreviousEval_.Init();

    XpreviousIt_.Resize(numElemSD);
    YpreviousIt_.Resize(numElemSD);
    /*
     * initialize previous x values to a normally unreachable value, such that we can trigger
     * a first calcuation in each case
     */
    //XpreviousIt_.Init(-inf);
    /*
     * In new version, we use XpreviousIt_ directly during the first computation
     * so it is not allowed to be inf at the point; use 0.0 instead;
     */
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

    dXpreviousItVEC_ = new Vector<Double>[numElemSD];
    dYpreviousItVEC_ = new Vector<Double>[numElemSD];

    XcurrentItVEC_ = new Vector<Double>[numElemSD];
    YcurrentItVEC_ = new Vector<Double>[numElemSD];

    XnextToLastTSVEC_ = new Vector<Double>[numElemSD];
    YnextToLastTSVEC_ = new Vector<Double>[numElemSD];

	  for(UInt i = 0; i < numElemSD; i++){
	    XpreviousEvalVEC_[i].Resize(dim_);
	    XpreviousEvalVEC_[i].Init();

	    YpreviousEvalVEC_[i].Resize(dim_);
	    YpreviousEvalVEC_[i].Init();

      XpreviousItVEC_[i].Resize(dim_);
      XpreviousItVEC_[i].Init();

      YpreviousItVEC_[i].Resize(dim_);
      YpreviousItVEC_[i].Init();

      dXpreviousItVEC_[i].Resize(dim_);
      dXpreviousItVEC_[i].Init();

      dYpreviousItVEC_[i].Resize(dim_);
      dYpreviousItVEC_[i].Init();

      XcurrentItVEC_[i].Resize(dim_);
      XcurrentItVEC_[i].Init();

      YcurrentItVEC_[i].Resize(dim_);
      YcurrentItVEC_[i].Init();

      XnextToLastTSVEC_[i].Resize(dim_);
      XnextToLastTSVEC_[i].Init();

      YnextToLastTSVEC_[i].Resize(dim_);
      YnextToLastTSVEC_[i].Init();
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

	// to calculate differential materialproperties, we need to know e0 / mu0
  if(material_->GetMaterialDatabaseName() == "Electrostatic"){
    rev_mat_fac_ = 8.854187817e-12;
    compute_inverse_ = false;
  } else if(material_->GetName() == "Magnetic"){
    rev_mat_fac_ = 1.2566370614e-6;
    compute_inverse_ = true;
  } else {
    EXCEPTION("Currently only Electrostatics and Magnetics are supported");
  }

  /*
   * Initialize non-delta matrix
   * -> scaled unit matrix
   */
  matTensor_ = Matrix<Double>(dim_,dim_);
  matTensor_.Init();
  for(UInt i = 0; i < dim_; i++){
    for(UInt j = 0; j < dim_; j++){
      if(i == j){
        if(compute_inverse_){
          matTensor_[i][j] = 1/rev_mat_fac_;
        } else {
          matTensor_[i][j] = rev_mat_fac_;
        }
      }
    }
  }

	PtrCoefFct matCoef = material_->GetTensorCoefFnc(matType_, tensorType_,
			                                         Global::REAL, false);

	matCoef->GetTensor(matTensorStart_, lpm);

	matDeltaTensor_ = new Matrix<Double>[numElemSD];

	for(UInt k = 0; k < numElemSD; k++){
	  //std::cout << "Initialize matDeltaTensor with: " << std::endl;
	  //std::cout << matTensorStart_.ToString() << std::endl;
	  matDeltaTensor_[k] = matTensor_;

	  /*
	   * former approach
	   */
//	  matDeltaTensor_[k] = matTensorStart_;
//
//	  for(UInt j = 0; j < matTensorStart_.GetNumRows(); j++){
//	    for(UInt i = 0; i < matTensorStart_.GetNumCols(); i++){
//	      if(i != j){
//	          matDeltaTensor_[k][j][i] = 0.0;
//	      }
//	    }
//	  }
	}
}

void CoefFunctionHyst::GetTensor( Matrix<Double>& tensor,
                                  const LocPointMapped& lpm ) {

  //std::cout << "Get Back Material Tensor" << std::endl;

  const Elem * el = lpm.ptEl;

  UInt idx = globalElem2Local_[el->elemNum];

  bool fullYInsteadOfDY = true;

  /*
   * set flag to false to avoid recalculationl
   */
  if(methodType_ == SCALAR){
    /*
     * Check if this approach is sufficient or if we have to adapt it
     * so that it returns only rev_mat_fac_ for terms on rhs
     * -> see CoefFunctionHyst.hh for detailed consideration
     */

    Double Xval, Ycurrent;

    ComputeXY(lpm, Xval, Ycurrent, overwrite_);
    //UInt version = 1;

    Double dX,dY;

    if(useNextToLastTS_){
      dX = Xval - XnextToLastTS_[idx];
      dY = Ycurrent - YnextToLastTS_[idx];
    } else {
      dX = Xval - XpreviousIt_[idx];
      dY = Ycurrent - YpreviousIt_[idx];
    }

    if(fullYInsteadOfDY){
      dY = Ycurrent;
    }

    Vector<Double> dX_vec = Vector<Double>(dim_);
    dX_vec.Init();
    dX_vec[dirP_] = dX;

    Vector<Double> dY_vec = Vector<Double>(dim_);
    dY_vec.Init();
    dY_vec[dirP_] = dY;

    if(deltaComputation_ == false){
      //std::cout << "deltaComputation_ = false" << std::endl;
      /*
       * here we return only the scaled unit matrix
       */
      tensor = matTensor_;
    } else {
      CreateDeltaMatrix(dX_vec,dY_vec,matDeltaTensor_[idx],idx);

      /*
       * store last values
       */
      //dXpreviousItVEC_[idx] = dX;
      //dYpreviousItVEC_[idx] = dY;
      tensor = matDeltaTensor_[idx];
    }

// geht nicht:
//
//    if (abs(dX_used) < tol_) {
//
//      if(version == 0){
//      //to be discussed!!
//        tensor = matTensorStart_;
//      } else {
//        for ( UInt i=0; i<matDeltaTensor_[idx].GetNumRows(); i++){
//          matDeltaTensor_[idx][i][i] = rev_mat_fac_used;
//        }
//        tensor = matDeltaTensor_[idx];
//      }
//    } else {
//      /*
//       * seems to work without abs -> why does it not work in vector case, too?
//       */
//      matDiff = dY_used / dX_used;
//      for ( UInt i=0; i<matDeltaTensor_[idx].GetNumRows(); i++){
//        matDeltaTensor_[idx][i][i] = rev_mat_fac_used;
//      }
//      matDeltaTensor_[idx][dirP_][dirP_] += matDiff;
//      tensor = matDeltaTensor_[idx];
//    }

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

     if(deltaComputation_ == false){
       //std::cout << "deltaComputation_ = false" << std::endl;
       /*
        * here we return only the scaled unit matrix
        */
       tensor = matTensor_;
     } else {
       //std::cout << "deltaComputation_ = true" << std::endl;
       /*
        * Compute that actutal deltaMaterial matrix
        *
        *  deltaMat = deltaY/deltaX + rev_mat_fac
        *
        */
       /*
        * Get vectors dX and dY
        */
       Vector<Double> Xcurrent, Ycurrent;

       /*
        * Note that ComputeXY_vec will check if Xcurrent was already used by comparing it to
        * XcurrentIt_
        * -> avoid reevaluation
        */

       std::cout << "GetTensor -> ComputeXY_vec with overwrite = " << overwrite_ << std::endl;

       //ComputeXY_vec(lpm, Xcurrent, Ycurrent, true);
       ComputeXY_vec(lpm, Xcurrent, Ycurrent, overwrite_);

       Vector<Double> dX = Xcurrent;
       Vector<Double> dY = Ycurrent;
       /*
        * check change to old iteration (not to old time step!)
        */

       if(useNextToLastTS_){
         /*
          * use difference to next to last timestep
          */
         dX.Add(-1.0,XnextToLastTSVEC_[idx]);
         dY.Add(-1.0,YnextToLastTSVEC_[idx]);
       } else {
         /*
          * use difference to lastIteration
          */
         dX.Add(-1.0,XpreviousItVEC_[idx]);
         dY.Add(-1.0,YpreviousItVEC_[idx]);
       }

       if(fullYInsteadOfDY){
         dY = Ycurrent;
       }

//       std::cout << "X: " << Xcurrent.ToString() << std::endl;
//       std::cout << "Y: " << Ycurrent.ToString() << std::endl;
//
//       std::cout << "old_X: " << XpreviousItVEC_[idx].ToString() << std::endl;
//       std::cout << "old_Y: " << YpreviousItVEC_[idx].ToString() << std::endl;

       /*
        * check if the change changed
        */
       Vector<Double> tmpX = dX;
       Vector<Double> tmpY = dY;
       tmpX -= dXpreviousItVEC_[idx];
       tmpY -= dYpreviousItVEC_[idx];
//
     //  std::cout << "dX: " << dX.ToString() << std::endl;
     //  std::cout << "dY: " << dY.ToString() << std::endl;

//       std::cout << "dX_norm_rel: " << dX.NormL2()/Xcurrent.NormL2() << std::endl;
//       std::cout << "dY_norm_rel: " << dY.NormL2()/Ycurrent.NormL2() << std::endl;
//
//       std::cout << "old_dX: " << dXpreviousItVEC_[idx].ToString() << std::endl;
//       std::cout << "old_dY: " << dYpreviousItVEC_[idx].ToString() << std::endl;

       if((dX.NormL2()/Xcurrent.NormL2() < -1e-17) || dY.NormL2()/Ycurrent.NormL2() < -1e-17){
         //std::cout << "Reuse old value as rel norm is too small" << std::endl;
         /*
          * reuse old matrix
          */
         tensor = matDeltaTensor_[idx];
       } else {
         if((tmpX.NormL2() < -tol_) && (tmpY.NormL2() < -tol_)){
           //std::cout << "Reuse old value as change to last dX/dY is too small" << std::endl;
           /*
            * reuse old matrix
            */
           tensor = matDeltaTensor_[idx];
         } else {
           //std::cout << "Compute new value" << std::endl;
           /*
            * calculate new matrix
            */
           CreateDeltaMatrix(dX,dY,matDeltaTensor_[idx],idx);

           /*
            * store last values
            */
           //dXpreviousItVEC_[idx] = dX;
           //dYpreviousItVEC_[idx] = dY;
           tensor = matDeltaTensor_[idx];
         }
       }
     }
   }
  std::cout << "Returning tensor: " << std::endl;
  std::cout << tensor.ToString() << std::endl;
}

void CoefFunctionHyst::GetScalar( Double& valY,const LocPointMapped& lpm ) {

  if(methodType_ != SCALAR){
    EXCEPTION("Only implemented for scalar model");
  }
  Double valX;
  /*
   * ComputeXY checks if values were already used; in that case, no recomputation is done
   */
  ComputeXY(lpm, valX, valY, overwrite_);
}

void CoefFunctionHyst::GetVector( Vector<Double>& valY,const LocPointMapped& lpm) {

  if(methodType_ == SCALAR){
    Double valX_scal;
    Double valY_scal;
    /*
     * ComputeXY checks if values were already used; in that case, no recomputation is done
     */
    ComputeXY(lpm, valX_scal, valY_scal, overwrite_);

    valY.Resize(dependCoef_->GetVecSize());
    valY.Init();
    valY[dirP_] = valY_scal;
  }
  else if(methodType_ == VECTOR){
    Vector<Double> valX;
    /*
     * ComputeXY checks if values were already used; in that case, no recomputation is done
     */
    ComputeXY_vec(lpm, valX, valY, overwrite_);

    std::cout << "GetVector: " << std::endl;
    std::cout << "E_in:  " << valX.ToString() << std::endl;
    std::cout << "P_out:  " << valY.ToString() << std::endl;

  } else {
    EXCEPTION("methodType_ " << methodType_ << " unknown")
  }
}

void CoefFunctionHyst::SetPreviousHystVals(bool setNextToLastTS_too) {
/*
 * This functions set the values of the last ITERATION
 * (for transient simulations, the value of the last time step is simply the last
 * iteration value, too)
 */

  if(overwrite_ == false){
    WARN("SetPreviousHystVals was called while Hysteresis operator is locked. This should not be the case.")
    return;

  } else {
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
        ComputeXY(lpm, X, Y, overwrite_);
        XpreviousIt_[idx] = X;
        YpreviousIt_[idx] = Y;

        if(setNextToLastTS_too){
          XnextToLastTS_[idx] = X;
          YnextToLastTS_[idx] = Y;
        }

      } else if(methodType_ == VECTOR){
        Vector<Double> X,Y;
        ComputeXY_vec(lpm, X, Y, overwrite_);
        /*
         * set value of last iteration
         */
       // std::cout << "Set previous hyst values" << std::endl;
        XpreviousItVEC_[idx] = X;
        YpreviousItVEC_[idx] = Y;

        if(setNextToLastTS_too){
          XnextToLastTSVEC_[idx] = X;
          YnextToLastTSVEC_[idx] = Y;
        }

       // std::cout << "Prev X: "<< XpreviousItVEC_[idx].ToString() << std::endl;
       // std::cout << "Prev Y: "<< YpreviousItVEC_[idx].ToString() << std::endl;
      }
    }
  }
}

void CoefFunctionHyst::ComputeXY_vec( const LocPointMapped& lpm, Vector<Double>& X, Vector<Double>& Y, bool overwrite) {

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

  bool truncate = !true;

  //std::cout << "Input X before truncation: " << X.ToString() << std::endl;
  if(truncate){
    Double tolInv = 1e14;
    Double tol = 1e-14;
    Double X_absmax = 0.0;
    int idx = -1;
    for(int i = 0; i < (int) dim_; i++){
      if(abs(X[i])>X_absmax){
        X_absmax = abs(X[i]);
        idx = i;
      }
    }

    if(idx != -1){
      /*
       * input not zero
       * -> truncate smaller entries
       */
      for(int i = 0; i < (int) dim_; i++){
        if(i == idx){
          continue;
        }
        Double X_mod = X[i]*tolInv;;
       //std::cout << "1. X_mod: " << X_mod << std::endl;
        X_mod = X_mod/X_absmax;
        //std::cout << "2. X_mod: " << X_mod << std::endl;
        X_mod = trunc(X_mod);
        //std::cout << "3. X_mod: " << X_mod << std::endl;
        X_mod *= X_absmax;
        //std::cout << "4. X_mod: " << X_mod << std::endl;
        X_mod *= tol;
        //std::cout << "5. X_mod: " << X_mod << std::endl;
        X[i] = X_mod;
      }
    }
  }
  //std::cout << "Input X after truncation: " << X.ToString() << std::endl;

  std::cout << "ComputeXY_vec: " << std::endl;
  std::cout << "overwrite_ " << overwrite_ << std::endl;
  std::cout << "overwrite " << overwrite << std::endl;

  Vector<Double> tmp = X;
  /*
   * Compare to last value that was evaluated at this point
   */
  tmp -= XpreviousEvalVEC_[idx];
  if(tmp.NormL2() < -tol_){
    std::cout << "Reuse last input! No reevaluation needed. " << std::endl;
    /*
     * no change in input -> (hopefully) no change in output; in case of hysteresis it is possible
     *                          that the same input leads to a different output due to the internal state
     *                          but as the operator was evaluated with the same value the last time it was called,
     *                          it cannot have changed its memory
     */
    Y = YpreviousEvalVEC_[idx];
  } else {
    if(overwrite != overwrite_){
      std::cout << "overwrite != overwrite_ and reevaluation is needed!" << std::endl;
    }
    /*
     * !!!!! Do not normalize here !!!!!
     * 1. the Hysteresis operator can handle input > xSat_; it will simply drive the system to its bounds
     *    (i.e. completely fills the Preisach plane)
     *    -> Normalizing not needed
     * 2. by restricting the input X to the saturation value xSat_, we will derive different
     *    slopes by computing dY/dX below
     *    -> wrong results!
     */

//    if(X.NormL2() > xSat_){
//      //std::cout << "Normalized X to xSat" << std::endl;
//      X = X*xSat_/X.NormL2();
//    }

    Y = hyst_->computeValue_vec(X, idx, overwrite,overwriteDirection_);

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
    static UInt firstIdx = idx;

    if(printOut_ > 0){
      if((cnt%printOut_==0)&&(idx==firstIdx)){
        std::cout << "Outputting bmp" << std::endl;
        std::stringstream filenamebuf;
        filenamebuf << "Switch_Elem"<<firstIdx<<"_Step" << std::setfill('0') << std::setw(5)<<cnt<<"_v"<<evalVersion_<<"_numRows"<<numRows_<<".bmp";
        hyst_->switchingStateToBmp(bmpResolution_,filenamebuf.str(),idx,true);
      }
    }
    if(idx == firstIdx){
      cnt++;
    }
  }
}

void CoefFunctionHyst::ComputeXY( const LocPointMapped& lpm, Double& X, Double& Y, bool overwrite) {

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
	X = elemSol[dirP_];

  std::cout << "GetVector (SCALAR): in " << X << std::endl;

	/*
	 * Compare to input of last EVALUATION (not iteration!)
	 */
  std::cout << "X-XpreviousEval_[idx]: " << X-XpreviousEval_[idx] << std::endl;

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
	   *  do not check here, but let Preisach operator do the work
	   */
	  std::cout << "Input did not change -> return old value" << std::endl;
	  Y = YpreviousEval_[idx];
	} else {
	  /*
	   * if input is too large (e.g. at electrode corners) deltaX may be huge
	   * whereas deltaY might still be small; the resulting deltaY/deltaX seems to cause trouble
	   * with adapted values, scalar model using delta material converges (at least better)
	   */

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

	  if(X > xSat_){
	    std::cout << "X > xSat_" << std::endl;
	    X = xSat_;
	  } else if(X < -xSat_){
	    std::cout << "X < xSat_" << std::endl;
	    X = -xSat_;
	  }

	  Y = hyst_->computeValueAndUpdate(X, idx, overwrite);

	  YpreviousEval_[idx] = Y;
	}

	std::cout << " -> out " << Y << std::endl;
}

std::string CoefFunctionHyst::ToString() const {
  return "CoefFunctionHyst";
}

void CoefFunctionHyst::CreateDeltaMatrix(Vector<Double>& dX,Vector<Double>& dY, Matrix<Double>& deltaMat, UInt idxElem){

  /* TODO:
   * For magnetics, we have to return the inverse as we need \nu instead of \mu.
   * -> Idea: compute dX/dY then add unit matrix with 1/\mu_0 on diagonal
   */
  if(material_->GetMaterialDatabaseName() == "Magnetic"){
    EXCEPTION("Not implemented yet for magnetics! Fix problems with electrostatics first!")
  }

  Double rev_mat_fac_used;
  Vector<Double> dX_used;
  Vector<Double> dY_used;

  if(compute_inverse_){
    rev_mat_fac_used = 1/rev_mat_fac_;
    dX_used = dY;
    dY_used = dX;
  } else {
    rev_mat_fac_used = rev_mat_fac_;
    dX_used = dX;
    dY_used = dY;
  }

  Matrix<Double> deltaMatBak = Matrix<Double>(deltaMat);
  /*
   * reset deltaMatrix
   */
  deltaMat.Resize(dim_,dim_);
  deltaMat.Init();

  std::string version_name = "DirectDivisionAbs"; // "L1_Regularized_LeastSquares";// "Unregularized_LeastSquares";

  if(version_name == "DirectDivision_OnlyMaxDirectionMod"){
    /*
     * Take the most straight forward approach
     * eps_ii = dY_used_i/dX_used_i
     * eps_ij = 0 (i != j)
     * but apply it only to the maximal component
     */
    Double dy_max = 0.0;
    Double tmp = 0.0;
    UInt idx = 0;
    for(UInt i = 0; i < dim_; i++){
      if(dX_used[i] != 0){
        tmp = abs(dY_used[i]);
      }
      if(tmp > dy_max){
        dy_max = tmp;
        idx = i;
      }
    }

    for(UInt i = 0; i < dim_; i++){
      if( (dX_used[i] != 0) && (abs(dY_used[i]/dy_max)>1e-2) ){
        deltaMat[i][i] = dY_used[i]/dX_used[i];
      }
    }
  }

  else if(version_name == "DirectDivision_OnlyMaxDirection"){
    /*
     * Take the most straight forward approach
     * eps_ii = dY_used_i/dX_used_i
     * eps_ij = 0 (i != j)
     * but apply it only to the maximal component
     */
    Double dy_max = 0.0;
    Double tmp = 0.0;
    UInt idx = 0;
    for(UInt i = 0; i < dim_; i++){
      if(dX_used[i] != 0){
        tmp = abs(dY_used[i]);
      }
      if(tmp > dy_max){
        dy_max = tmp;
        idx = i;
      }
    }

    if(dX_used[idx] != 0){
      deltaMat[idx][idx] = dY_used[idx]/dX_used[idx];
    }

  }
  else if(version_name == "DirectDivision"){
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
    for(UInt i = 0; i < dim_; i++){
      if(dX_used[i] != 0){
        std::cout << "dX_used[i] = " << dX_used[i] << std::endl;
        std::cout << "dY_used[i] = " << dY_used[i] << std::endl;
        deltaMat[i][i] = dY_used[i]/dX_used[i];
      } else {
        std::cout << "dX_used[i] = " << dX_used[i] << std::endl;
        std::cout << "dY_used[i] = " << dY_used[i] << std::endl;
        std::cout << "Reuse old value as dX == 0" << std::endl;
        deltaMat[i][i] = deltaMatBak[i][i];
      }
    }
  }
  else if(version_name == "DirectDivisionAbs"){
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
      for(UInt i = 0; i < dim_; i++){
        if(dX_used[i] != 0){
          std::cout << "dX_used[i] = " << dX_used[i] << std::endl;
          std::cout << "dY_used[i] = " << dY_used[i] << std::endl;
          deltaMat[i][i] = abs(dY_used[i]/dX_used[i]);
        } else {
          std::cout << "dX_used[i] = 0" << std::endl;
          std::cout << "dY_used[i] = " << dY_used[i] << std::endl;
          std::cout << "Reuse old value as dX == 0" << std::endl;
          deltaMat[i][i] = deltaMatBak[i][i];
        }
      }
    }
  else if(version_name == "DirectDivisionMod"){
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
    for(UInt i = 0; i < dim_; i++){
      if( (dX_used[i] != 0) && (abs(dY_used[i]) > 1e-15) ) {
        std::cout << "dX_used[i] = " << dX_used[i] << std::endl;
        std::cout << "dY_used[i] = " << dY_used[i] << std::endl;
      //if( (dX_used[i] != 0) && (abs(dY_used[i]) > 1e-16) ) {
        deltaMat[i][i] = dY_used[i]/dX_used[i];
      } else {

        std::cout << "dX_used[i] = " << dX_used[i] << std::endl;
        std::cout << "dY_used[i] = " << dY_used[i] << std::endl;
      }
    }
  }
  else if(version_name == "DirectDivisionMod2"){
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
    for(UInt i = 0; i < dim_; i++){
      /*
       * set tolerance for dX such large, that in the worst case of only rounding
       * errors in dY, i.e. dY = tol_y, the resulting eps will not be orders of
       * magnitude larger than eps
       * -> this happens especially in the beginning of very simple setups
       * (only 1 fix direction)
       * -> due to rounding errors in dY, we suddenly get 1e-20/1e-12=1e-8 >>> eps0 = 1e-12
       * -> even worse, this is different for each element (some might habve 1e-20, others
       * 2e-19 or 3e-18
       * -> huge differences in material tensor
       * -> upcoming iterations will magnify this problem as they use this warped
       * material tensor as basis
       */
      Double tol_y = 1e-15;
      if( (abs(dX_used[i]) > tol_y/rev_mat_fac_used) && (abs(dY_used[i]) > tol_y) ) {
      //if( (dX_used[i] != 0) && (abs(dY_used[i]) > 1e-16) ) {
        deltaMat[i][i] = dY_used[i]/dX_used[i];
      } else {
        std::cout << "dX_used[i] = 0" << std::endl;
        std::cout << "dY_used[i] = " << dY_used[i] << std::endl;
      }
    }
  }
  else if(version_name == "L1_Regularized_LeastSquares"){
    /*
     * Perform the shooting algorithm to solve least squares problem
     *  y = Ax
     * with a L1 regularization, i.e. norm_L1(x) < lambda
     * with norm_L1(x) = sum_i abs(x_i)
     *
     * Definition:
     * 2D case:
     *
     *  deltaMat = [ [d11 d12]; [d21 d22] ] = [ [d11 d12]; [d12 d22] ]
     *  dX = [dX1 dX2]
     *  dY = [dY1 dY2]
     *
     *  x = [d11 d22 d12] (Voigt)
     *  A = [[dX1 0 dX2]
     *       [0 dX2 dX1]]
     *  AT = [[dX1   0]
     *        [  0 dX2]
     *        [dX2 dX1]]
     *  y = [dY1 dY2]
     *
     * 3D case:
     *
     *  deltaMat = [ [d11 d12 d13]; [d21 d22 d23]; [d31 d32 d33] ] = [ [d11 d12 d13]; [d12 d22 d23]; [d13 d23 d33] ]
     *  dX = [dX1 dX2 dX3]
     *  dY = [dY1 dY2 dY3]
     *
     * x = [d11 d22 d33 d23 d13 d12] (Voigt)
     * A = [[dX1   0   0   0 dX3 dX2]
     *      [  0 dX2   0 dX3   0 dX1]
     *      [  0   0 dX3 dX2 dX1   0]]
     * AT = [[dX1   0   0]
     *       [  0 dX2   0]
     *       [  0   0 dX3]
     *       [  0 dX3 dX2]
     *       [dX3   0 dX1]
     *       [dX2 dX1   0]]
     * y = [dY1 dY2 dY3]
     *
     *
     * Shooting algorithm:
     *
     *  start with initial guess for x
     *  precompute ATA = A^T*A
     *  set lambda
     *  S = Vector<Double>(length(x))
     *
     *  while( norm_L2(y-Ax) > tol )
     *    y = Ax
     *    S = -A^T*y + ATA*x
     *    for(j = 0; j < length(x); j++)
     *      S(j) -= ATA(j,j)*x(j)
     *
     *    for(j = 0; j < length(x); j++)
     *      if S(J) > lambda
     *        x(j) = (lambda - S(J))/ATA(j,j)
     *      else if S(J) < -lambda
     *        x(j) = (-lambda - S(J))/ATA(j,j)
     *      else
     *        x(j) = 0.0
     */


    Matrix<Double> AT;
    UInt dim2 = 0;

    if(dim_ == 2){
      // number of unknowns of the material tensor
      dim2 = 3;
      /*
       * start with transposed matrix (such that we can use MultT later on
       *
       *  x = [d11 d22 d21] (Voigt)
       *  A = [[dX1 0   dX2]
       *       [  0 dX2 dX1]]
       *  AT = [[dX1   0]
       *        [  0 dX2]
       *        [dX2 dX1]]
       *  y = [dY1 dY2]
       */
      AT = Matrix<Double>(dim2,dim_);
      AT.Init();
      AT[0][0] = dX_used[0];
      AT[1][1] = dX_used[1];
      AT[2][0] = dX_used[1];
      AT[2][1] = dX_used[0];

    } else if (dim_ == 3) {
      // number of unknowns of the material tensor
      dim2 = 6;
      /*
       * x = [d11 d22 d33 d32 d31 d21] (Voigt)
       * A = [[dX1   0   0   0 dX3 dX2]
       *      [  0 dX2   0 dX3   0 dX1]
       *      [  0   0 dX3 dX2 dX1   0]]
       * AT = [[dX1   0   0]
       *       [  0 dX2   0]
       *       [  0   0 dX3]
       *       [  0 dX3 dX2]
       *       [dX3   0 dX1]
       *       [dX2 dX1   0]]
       * y = [dY1 dY2 dY3]
       */
      AT = Matrix<Double>(dim2,dim_);
      AT.Init();
      AT[0][0] = dX_used[0];
      AT[1][1] = dX_used[1];
      AT[2][2] = dX_used[2];
      AT[3][1] = dX_used[2];
      AT[3][2] = dX_used[1];
      AT[4][0] = dX_used[2];
      AT[4][2] = dX_used[0];
      AT[5][0] = dX_used[1];
      AT[5][1] = dX_used[0];
    }

    Matrix<Double> ATA = Matrix<Double>(dim2,dim2);
    Vector<Double> x = Vector<Double>(dim2);
    Vector<Double> S = Vector<Double>(dim2);
    Vector<Double> ATy = Vector<Double>(dim2);
    Vector<Double> Ax = Vector<Double>(dim_);
    Vector<Double> res = Vector<Double>(dim2);
    Matrix<Double> A;
    AT.Transpose(A);

    AT.Mult(A,ATA);

    for(UInt j = 0; j < dim_; j++){
      /*
       * take rev_mat_fac_used*Identtiy as initial guess for delta matrix
       * i.e. set the first dim_ entries of x to rev_mat_fac_used
       */
      x[j] = rev_mat_fac_used;
    }

    Double lambda = 1e-12;
    Double tol = 1e-14;
    UInt maxIters = 2;

    /*
     * compute initial residual
     */
    Double x_l1Norm = 0.0;
    Double res_Norm = 0.0;
    for(UInt i = 0; i < x.GetSize(); i++){
      x_l1Norm += abs(x[i]);
    }
    res = -dY_used;
    A.Mult(x,Ax);
    res += Ax;
    res_Norm = res.NormL2() + lambda*x_l1Norm;

    /*
     * precompute ATy
     */
    AT.Mult(dY_used,ATy);

    UInt cnt = 0;
    while(res_Norm > tol){
      ATA.Mult(x,S);
      S -= ATy;
      for(UInt j = 0; j < dim2; j++){
        S[j] -= ATA[j][j]*x[j];

        Double tmp;
        if(S[j] > lambda ){
          tmp = ( lambda - S[j] )/ATA[j][j];
        } else if(S[j] < -lambda ){
          tmp = ( -lambda - S[j] )/ATA[j][j];
        } else {
          tmp = 0.0;
        }
        x[j] = tmp; //0.1*x[j]+0.9*tmp;
      }

      cnt++;
      std::cout << "after iteration " << cnt << std::endl;

      /*
       * calc residual
       */
      x_l1Norm = 0.0;
      for(UInt i = 0; i < x.GetSize(); i++){
        x_l1Norm += abs(x[i]);
      }
      std::cout << "x: " << x.ToString() << std::endl;

      res = -dY_used;
      A.Mult(x,Ax);
      res += Ax;

      std::cout << "y: " << dY_used.ToString() << std::endl;
      std::cout << "A: " << A.ToString() << std::endl;
      std::cout << "Ax: " << Ax.ToString() << std::endl;
      std::cout << "y-Ax: " << res.ToString() << std::endl;
      res_Norm = res.NormL2() + lambda*x_l1Norm;

      std::cout << "|y-Ax|: " << res.NormL2() << std::endl;
      std::cout << "x_l1Norm: " << x_l1Norm << std::endl;
      std::cout << "lambda*x_l1Norm: " << x_l1Norm*lambda << std::endl;
      std::cout << "Res norm: " << res_Norm << std::endl;
      if(cnt > maxIters){
        break;
      }
    }

    // sort vector into matrix
    for(UInt i = 0; i < dim_; i++){
      for(UInt j = 0; j <= i; j++){
        if(i == j){
          deltaMat[i][i] = x[i];
        } else {
          deltaMat[i][j] = x[dim2-i-j];
          deltaMat[j][i] = deltaMat[i][j];
        }
      }
    }

  }

  else if(version_name == "Unregularized_LeastSquares"){

    /*
     * Get vector norms
     */
    Double norm_dX = dX_used.NormL2();

    /*
     * If norm_dX < tol_, it is a bad idea to compute dY/dX
     * -> in that case, norm_dY should be quite small also
     * (can be checked); however, return simply matrix
     * with rev_mat_fac_ on diagonal
     */
    if(norm_dX > tol_){
      //std::cout << "NormDX = " << norm_dX << " is larger than tol" << std::endl;
      std::cout << "dX: " << dX_used << std::endl;
      std::cout << "dY: " << dY_used << std::endl;
      /*
       * Solve least square problem to determine delta matrix [deltaMat], such that
       * dY = [deltaMat]*dX
       *
       * Assumptions: deltaMat = symmetric
       * Procedure (for 2D, uncoupled case):l
       *  1. definition
       *  deltaMat = [ [d11 d12]; [d21 d22] ] = [ [d11 d12]; [d12 d22] ]
       *  dX = [dX1 dX2]
       *  dY = [dY1 dY2]
       *
       *  x = [d11 d22 d12] (Voigt)
       *  A = [[dX1 0 dX2]
       *       [0 dX2 dX1]]
       *  y = [dY1 dY2]
       *
       *  2. solve least squares problem
       *    Ax = y
       *  by solving
       *    (A*A^T)z = y
       *  and setting
       *    x = (A^T)z
       *
       *  3. create deltaMat from x
       *
       *  Note: solution may lead to very small diagonal entries for deltamat
       *        in that case consider a regularization, i.e. solve
       *          (omega*I + A*A^T)z = y
       *        with omega being a regularization parameter
       *
       * TODO: extend for coupled case
       */

      Matrix<Double> AT;
      UInt dim2 = 0;

      if(dim_ == 2){
        // number of unknowns of the material tensor
        dim2 = 3;
        /*
         * start with transposed matrix (such that we can use MultT later on
         *
         *  x = [d11 d22 d21] (Voigt)
         *  A = [[dX1 0   dX2]
         *       [  0 dX2 dX1]]
         *  AT = [[dX1   0]
         *        [  0 dX2]
         *        [dX2 dX1]]
         *  y = [dY1 dY2]
         */
        AT = Matrix<Double>(dim2,dim_);
        AT.Init();
        AT[0][0] = dX_used[0];
        AT[1][1] = dX_used[1];
        AT[2][0] = dX_used[1];
        AT[2][1] = dX_used[0];

      } else if (dim_ == 3) {
        // number of unknowns of the material tensor
        dim2 = 6;
        /*
         * x = [d11 d22 d33 d32 d31 d21] (Voigt)
         * A = [[dX1   0   0   0 dX3 dX2]
         *      [  0 dX2   0 dX3   0 dX1]
         *      [  0   0 dX3 dX2 dX1   0]]
         * AT = [[dX1   0   0]
         *       [  0 dX2   0]
         *       [  0   0 dX3]
         *       [  0 dX3 dX2]
         *       [dX3   0 dX1]
         *       [dX2 dX1   0]]
         * y = [dY1 dY2 dY3]
         */
        AT = Matrix<Double>(dim2,dim_);
        AT.Init();
        AT[0][0] = dX_used[0];
        AT[1][1] = dX_used[1];
        AT[2][2] = dX_used[2];
        AT[3][1] = dX_used[2];
        AT[3][2] = dX_used[1];
        AT[4][0] = dX_used[2];
        AT[4][2] = dX_used[0];
        AT[5][0] = dX_used[1];
        AT[5][1] = dX_used[0];
      }

      Matrix<Double> AAT = Matrix<Double>(dim_,dim_);
      /*
       * perform AAT = (AT)^T * AT
       */
      AT.MultT(AT,AAT);

      Matrix<Double> AATI = Matrix<Double>(dim_,dim_);

      // add regualrization
  //      double omega = 0.0;
  //      for(UInt i = 0; i < dim_; i++){
  //        AAT[i][i] += omega;
  //      }

      AAT.Invert(AATI);

      Vector<Double> z = Vector<Double>(dim_);
      Vector<Double> x = Vector<Double>(dim2);

      // solve least squares problem
      z = AATI*dY_used;

      // recover vector of unknowns
      x = AT*z;

      // sort vector into matrix
      for(UInt i = 0; i < dim_; i++){
        for(UInt j = 0; j <= i; j++){
          if(i == j){
            deltaMat[i][i] = x[i];
          } else {
            deltaMat[i][j] = x[dim2-i-j];
            deltaMat[j][i] = deltaMat[i][j];
          }
        }
      }
    } else {
      std::cout << "NormDX = " << norm_dX << " is smaller than tol" << std::endl;
      //norm_dX < tol_ -> only return scaled unit matrix
      // -> in the inverse case, where dX_used represents the difference of the Hysteresis output, norm_dX becomes
      // 0 once the material is in saturation
    }
  }

//  std::cout << "deltaMat before adding std term: " << std::endl;
//  std::cout << deltaMat.ToString() << std::endl;

  // add rev_mat_fac_used after computation of matDeltaTensorl
  //deltaMat += matTensor_;

  for(UInt i = 0; i < dim_; i++){
    if(deltaMat[i][i] < 0){
      deltaMat[i][i] -= rev_mat_fac_used;
    } else {
      deltaMat[i][i] += rev_mat_fac_used;
    }
  }

}

}// namespace
