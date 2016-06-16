#include "CoefFunctionHyst.hh"

// classes for function / spline approximation
#include "Materials/Models/Preisach.hh"
#include "Materials/Models/VectorPreisach.hh"
#include "Materials/Models/VectorPreisachv7.hh"
#include "FeBasis/FeFunctions.hh"
#include "FeBasis/FeSpace.hh"
#include "Forms/Operators/BaseBOperator.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "DataInOut/Logging/LogConfigurator.hh"


namespace CoupledField {

  DECLARE_LOG(coeffcthyst)
  DEFINE_LOG(coeffcthyst, "coeffcthyst")

CoefFunctionHyst::CoefFunctionHyst( BaseMaterial* const material,
                    shared_ptr<ElemList> actSDList,
                    PtrCoefFct dependency,
					SubTensorType tensorType,
					MaterialType matType) : CoefFunction() {

  // this type of coefficient is nonlinear (i.e. solution dependent)
  dimType_ = VECTOR;
  dependType_ = SOLUTION;
  isAnalytic_ = false;
  isComplex_  = false;
  dependCoef_ = dependency;
  tensorType_ = tensorType;
  matType_    = matType;
  material_   = material;
  tol_ = 1e-15;

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

	 // std::cout << "Using Scalar Hysteresis" << std::endl;

    //get direction
    std::string str;
    material->GetScalar(str, P_DIRECTION);
    Directions dir;
    String2Enum(str,dir);
    dirP_ = dir;

    hyst_ = new Preisach(numElemSD, Xsat, Ysat, weights, isVirgin);

    dim_ = 1;
    evalVersion_ = 0;
	} else if(methodType_ == VECTOR){

	  //std::cout << "Using Vector Hysteresis" << std::endl;

	  material->GetScalar(rotRes_, ROT_RESISTANCE, Global::REAL);

    int eval;
    material->GetScalar(eval, EVAL_VERSION);

    evalVersion_ = (UInt) eval;
    int isTesting;
    material->GetScalar(isTesting, IS_TESTING);

	  dim_ = dependCoef_->GetVecSize();
	  if(evalVersion_ == 7){
	    //std::cout << "Using Vector Hystersis v7" << std::endl;
	    hyst_ = new VectorPreisachv7(numElemSD, Xsat, Ysat, weights,rotRes_,dim_, isVirgin, isTesting!=0, (UInt) evalVersion_);
	  } else {
	    //std::cout << "Using Vector Hystersis v" << evalVersion_ << std::endl;
	    hyst_ = new VectorPreisach(numElemSD, Xsat, Ysat, weights,rotRes_,dim_, isVirgin, isTesting!=0, (UInt) evalVersion_);
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

  Double inf = std::numeric_limits<double>::infinity();

	//allocate memory for previous results, needed for the
	//effective material parameter formulation
	if(methodType_ == SCALAR){
    Xprevious_.Resize(numElemSD);
    Yprevious_.Resize(numElemSD);
    Xprevious_.Init();
    Yprevious_.Init();

    XpreviousIt_.Resize(numElemSD);
    YpreviousIt_.Resize(numElemSD);
    /*
     * initialize previous x values to a normally unreachable value, such that we can trigger
     * a first calcuation in each case
     */
    XpreviousIt_.Init(-inf);
    YpreviousIt_.Init();

	} else if(methodType_ == VECTOR){
	  XpreviousVEC_ = new Vector<Double>[numElemSD];
	  YpreviousVEC_ = new Vector<Double>[numElemSD];

    XpreviousItVEC_ = new Vector<Double>[numElemSD];
    YpreviousItVEC_ = new Vector<Double>[numElemSD];

	  for(UInt i = 0; i < numElemSD; i++){
	    XpreviousVEC_[i].Resize(dim_);
	    XpreviousVEC_[i].Init();

      YpreviousVEC_[i].Resize(dim_);
      YpreviousVEC_[i].Init();

      XpreviousItVEC_[i].Resize(dim_);
      XpreviousItVEC_[i].Init(-inf);

      YpreviousItVEC_[i].Resize(dim_);
      YpreviousItVEC_[i].Init();
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

	PtrCoefFct matCoef = material_->GetTensorCoefFnc(matType_, tensorType_,
			                                         Global::REAL, false);

	matCoef->GetTensor(matTensorStart_, lpm);
	
	matDeltaTensor_ = new Matrix<Double>[numElemSD];

	for(UInt k = 0; k < numElemSD; k++){
	  //std::cout << "Initialize matDeltaTensor with: " << std::endl;
	  //std::cout << matTensorStart_.ToString() << std::endl;
	  matDeltaTensor_[k] = matTensorStart_;
	  for(UInt j = 0; j < matTensorStart_.GetNumRows(); j++){
	    for(UInt i = 0; i < matTensorStart_.GetNumCols(); i++){
	      if(i != j){
	          matDeltaTensor_[k][j][i] = 0.0;
	      }
	    }
	  }
	}

	  // to calculate differential materialproperties, we need to know e0/ mu0
  if(material_->GetMaterialDatabaseName() == "Electrostatic"){
    rev_mat_fac_ = 8.854187817e-12;
  } else if(material_->GetName() == "Magnetic"){
    rev_mat_fac_ = 1.2566370614e-6;
  } else {
  }
  //std::cout << "Init CoefFncHyst passed" << std::endl;

}

void CoefFunctionHyst::GetTensor( Matrix<Double>& tensor,
                                  const LocPointMapped& lpm ) {


 // std::cout << "CoefFunctionHyst::GetTensor called" << std::endl;
  const Elem * el = lpm.ptEl;

  UInt idx = globalElem2Local_[el->elemNum];

    /*
     * set flag to false to avoid recalculation
     */
    if(methodType_ == SCALAR){

      Double Xval, Ycurrent, matDiff;
      ComputeXY(lpm, Xval, Ycurrent);
      UInt version = 1;

      //compute differential material parameter
      Double dX = Xval - Xprevious_[idx];

      // D = eps0 + P
      // Ycurrent-Yprevious = deltaP!
      Double dY = Ycurrent -Yprevious_[idx] + rev_mat_fac_*dX;

      //if ( (abs(dY) < rev_mat_fac_) || (abs(dX) < 1e-10) ) {
      if ( (abs(dY) < rev_mat_fac_) || (abs(dX) < 1e-10) ) {

        if(version == 0){
        //to be discussed!!
          tensor = matTensorStart_;
        } else {
          for ( UInt i=0; i<matDeltaTensor_[idx].GetNumRows(); i++){
            matDeltaTensor_[idx][i][i] = rev_mat_fac_;
          }
          tensor = matDeltaTensor_[idx];
        }
      } else {
        /*
         * seems to work without abs -> why does it not work in vector case, too?
         */
        matDiff = dY / dX;
        for ( UInt i=0; i<matDeltaTensor_[idx].GetNumRows(); i++){
          matDeltaTensor_[idx][i][i] = matDiff;
        }
        //matDeltaTensor_[dirP_][dirP_] = matDiff;
        tensor = matDeltaTensor_[idx];
      }

    } else if(methodType_ == VECTOR ){

      /*
       * NEW 22.2.16
       * create matrix R, such that dY = R*dX
       *
       * Algorithm:
       *  1. find vector v such that dY x v and dX x v both != 0; extend vector to 3D!
       *  2. create m1 = dX x v; m2 = v x m1; m3 = m1 x m2
       *            n1 = dY x v; n2 = v x n1; n3 = n1 x n2
       *     note: vectors should be normalized to length 1
       *  3. create matrices M = [m1; m2; m3] and N = [n1; n2; n3]
       *  4. R = N*inv(M)*norm(dY)/norm(dX)
       * (5. check if dY = R*dX
       *
       */

      /*
       * Get vectors dX and dY
       */
      Vector<Double> Xval, Ycurrent;

      /*
       * TODO:
       * BIG PROBLEM:
       *  the hysteresis operator has a memory of the past
       *  this is of course what we want but it means that each input-nonsense
       *  will be stored too
       *  this is currently happening in the deltaMaterial method if the non-diagonal
       *  elements are non-zero. In that case the whole methods works up to a certain point.
       *  At that point, the solution value explodes because of some unknown reason. Although
       *  the values will be normal later on, huge parts of the Preisach plane flipped. From now
       *  on the polarization/magnetization will have some kind of offset which prohibits the
       *  method to produce useful results.
       *  Maybe one has to "reset" the Preisach operator to the state before the iterations occured and
       *  then set it to the resulting value. This will however not prohibit that the solution process
       *  will converge to wrong results.
       *  Better solution: Prohibit explosion of input values!
       *
       */

      /*
       * DEBUGGING
       */
      // get region id
      UInt regionId = el->regionId;
      /*
       * for the test example we know that regionId=6 is the observer at the singularity and regionId=5 is the observer at the line of sym
       */

      if(regionId == 5){
         std::cout << "Element at middleline" << std::endl;
       } else if(regionId == 6){
         std::cout << "Element at singularity" << std::endl;
       } //else {

//      if((regionId == 5) || (regionId == 6)){
//        std::cout << "XpreIt: " << XpreviousItVEC_[idx].ToString() << std::endl;
//        std::cout << "YpreIt: " << YpreviousItVEC_[idx].ToString() << std::endl;
//        std::cout << "XpreTS: " << XpreviousVEC_[idx].ToString() << std::endl;
//        std::cout << "YpreTS: " << YpreviousVEC_[idx].ToString() << std::endl;
//      }

      ComputeXY_vec(lpm, Xval, Ycurrent, true);

      Vector<Double> dX = Xval;
      dX.Add(-1.0,XpreviousVEC_[idx]);

      Vector<Double> dY = Ycurrent;
      dY.Add(-1.0,YpreviousVEC_[idx]);

      dY.Add(rev_mat_fac_,dX);

      //matDeltaTensor_[idx].Init();
      CreateDeltaMatrix(dX,dY,matDeltaTensor_[idx],4,idx);

      if((regionId == 5) || (regionId == 6)){
        std::cout << "Xin: " << Xval.ToString() << std::endl;
        std::cout << "Yout: " << Ycurrent.ToString() << std::endl;
//        std::cout << "XpreIt (after call): " << XpreviousItVEC_[idx].ToString() << std::endl;
//        std::cout << "YpreIt (after call): " << YpreviousItVEC_[idx].ToString() << std::endl;
        std::cout << "dX: " << dX.ToString() << std::endl;
        std::cout << "dY: " << dY.ToString() << std::endl;
        std::cout << matDeltaTensor_[idx].ToString() << std::endl;
      }


  //    std::cout << "### v0 ###" << std::endl;
   //   std::cout << "Resulting matrix after v0: " << std::endl;
   //   std::cout << matDeltaTensor_.ToString() << std::endl;
  //    std::cout << "dY - matDeltaTensor *dX:" << std::endl;
  //    std::cout << dY - matDeltaTensor_*dX << std::endl;
  //    std::cout << std::endl;

  //    matDeltaTensor_.Init();
  //    CreateDeltaMatrix(dX,dY,matDeltaTensor_,4);
  //
  //    if(idx == 30){
  //
  ////    std::cout << "### v4 ###" << std::endl;
  //   std::cout << "Resulting matrix after v4: " << std::endl;
  //   std::cout << matDeltaTensor_.ToString() << std::endl;
  //    }
  //    std::cout << "dY - matDeltaTensor *dX:" << std::endl;
  //    std::cout << dY - matDeltaTensor_*dX << std::endl;
  //    std::cout << std::endl;
    //  matDeltaTensor_.Init();
    //  CreateDeltaMatrix(dX,dY,matDeltaTensor_,2);

  //    std::cout << "Resulting matrix after v2: " << std::endl;
  //    std::cout << matDeltaTensor_.ToString() << std::endl;

    //  matDeltaTensor_.Init();
    //  CreateDeltaMatrix(dX,dY,matDeltaTensor_,4);

  //
  //
  //    //Vector<Double> Xval, Ycurrent;
  //    //Double matDiff;
  //    ComputeXY_vec(lpm, Xval, Ycurrent);
  //
  //    Vector<Double> dXvec = Xval;
  //    dXvec.Add(-1.0,XpreviousVEC_[idx]);
  //    Vector<Double> dYvec = Ycurrent;
  //    dYvec.Add(-1.0,YpreviousVEC_[idx]);
  //
  //    Double dXNorm = dXvec.NormL2();
  //    Double dYNorm = dYvec.NormL2();
  //
  //    //std::cout << "dXvec = " << dXvec.ToString() << std::endl;
  //    //std::cout << "dXvec.NornL2() = " << dXNorm << std::endl;
  //    //std::cout << "dYvec = " << dYvec.ToString() << std::endl;
  //    //std::cout << "dYvec.NornL2() = " << dYNorm << std::endl;
  //
  //    Double dXmax = 0.0;
  //    UInt idxMax = 0;
  //    for(UInt i = 0; i < dim_; i++ ){
  //      if(abs(dXvec[i]) > dXmax){
  //        dXmax = abs(dXvec[i]);
  //        idxMax = i;
  //      }
  //    }
  //
  //    // assume that a material cannot have more than mat_factor = rev_rel_max*rev_mat_fac_
  //    Double rev_rel_max = 1e10;
  //    Double eps = 1e-10; // tolerance
  //    // reset matrix first;
  //    for(UInt i = 0; i < dim_; i++ ){
  //      for(UInt j = 0; j < dim_; j++){
  //        if(i != j){
  //          // Initialize nondiagonal entries
  //          matDeltaTensor_ = 0.0;
  //        }
  //      }
  //    }
  //
  //
  //    for(UInt i = 0; i < dim_; i++ ){
  //
  //      if(abs(dYvec[i]) < rev_mat_fac_){
  //        //std::cout << "Case 1: Will use eps0 = " << rev_mat_fac_ << std::endl;
  //        //std::cout << "dY = " << dYvec[i] << " / dX = " << dXvec[i] << " / dY/dX = " << dYvec[i]/dXvec[i] << std::endl;
  //
  //        matDeltaTensor_[i][i] = rev_mat_fac_;
  //      } else if(abs(dXvec[i]) < eps){
  //        //std::cout << "Case 2: Will use old value = " << matDeltaTensor_[i][i] << std::endl;
  //        //std::cout << "dY = " << dYvec[i] << " / dX = " << dXvec[i] << " / dY/dX = " << dYvec[i]/dXvec[i] << std::endl;
  //      } else if (abs(dYvec[i])/abs(dXvec[i]) > rev_rel_max*rev_mat_fac_ ){
  //     //    std::cout << "abs(dXvec[i]) < eps and abs(dYvec[i]) > rev_rel_max*eps*rev_mat_fac_!" << std::endl;
  //     //    std::cout << "dYvec[i] can thus not depend soleley on dXvec[i]" << std::endl;
  //     //    std::cout << "Therefore will add matrix entry to the offdiagonal relating dYvec[i] to max of dXvec[i]" << std::endl;
  //        matDeltaTensor_[i][i] = (dYvec[i]+rev_mat_fac_*dXmax)/dXmax;
  //        //std::cout << "Case 3: Will use scaling to dXmax and setting entry " << i << ", " << idxMax << " to " << matDeltaTensor_[i][idxMax] << std::endl;
  //        //std::cout << "dY = " << dYvec[i] << " / dX = " << dXvec[i] << " / dY/dX = " << dYvec[i]/dXvec[i] << std::endl;
  //      } else {
  //        matDeltaTensor_[i][i] = (dYvec[i]+rev_mat_fac_*dXvec[i])/dXvec[i];
  //        //std::cout << "Case 4: Will use std procedure and set " << matDeltaTensor_[i][idxMax] << std::endl;
  //        //std::cout << "dY = " << dYvec[i] << " / dX = " << dXvec[i] << " / dY/dX = " << dYvec[i]/dXvec[i] << std::endl;
  //      }
  //
  //    }
  //
  //    std::cout << "Resulting matrix old: " << std::endl;
  //    std::cout << matDeltaTensor_.ToString() << std::endl;

  //
  //    for(UInt i = 0; i < dim_; i++ ){
  //      Double dX = dXvec[i];
  //      Double dY = dYvec[i] + rev_mat_fac_*dX;
  //      Double tol = 1e-8; // tolerance
  //
  //      std::cout << "Component = " << i << std::endl;
  //      if ( (abs(dY) < rev_mat_fac_)  ) {
  //        std::cout << "Case 1: Will use eps0 = " << rev_mat_fac_ << std::endl;
  //        std::cout << "dY = " << dY << " / dX = " << dX << " / dY/dX = " << dY/dX << std::endl;
  //      //if ( (abs(dX) < 1e-10) ) {
  //        //to be discussed!!
  //        /*
  //         * dont do the tensor = matTensorStart_! it might override the previous iteration which already set useful values!
  //         */
  //        //tensor = matTensorStart_;
  //        matDeltaTensor_[i][i] = rev_mat_fac_;
  //      } else if ((abs(dX) < tol)){
  //        /* the following case might occur (and occurs here and there)
  //         *
  //         * dX = 0 but dY != 0
  //         * this can be the case if e.g.
  //         * X1 = 0.5,1 and X2 = 0.5,2
  //         * as long as abs(X) is not larger than saturation, X1 and X2 will be unchanged during calculation
  //         * we see that X1[0]-X2[0] = 0
  //         * however, the direction of X has changed and such Y might have changed, too!
  //         * this is why we get dY != 0 but dX = 0
  //         *
  //         * -> this kind of normalization is problematic as it would allow eps values of inf size
  //         *
  //         */
  //        std::cout << "Case 2: Will use old value = " << matDeltaTensor_[i][i] << std::endl;
  //        std::cout << "dY = " << dY << " / dX = " << dX << " / dY/dX = " << dY/dX << std::endl;
  //      } else {
  //        matDiff = abs(dY/dX);
  //        //std::cout << "Mat diff / component " << matDiff << " / " << i << std::endl;
  //        matDeltaTensor_[i][i] = matDiff;
  //        std::cout << "Case 3: Will use new diff = " << matDeltaTensor_[i][i] << std::endl;
  //        std::cout << "dY = " << dY << " / dX = " << dX << " / dY/dX = " << dY/dX << std::endl;
  //      }
  //    }
      tensor = matDeltaTensor_[idx];
    }
}

void CoefFunctionHyst::GetScalar( Double& valY,
                                  const LocPointMapped& lpm ) {

  if(methodType_ != SCALAR){
    EXCEPTION("Only implemented for scalar model");
  }
	Double valX;
	ComputeXY(lpm, valX, valY);
}

void CoefFunctionHyst::GetVector( Vector<Double>& valY,const LocPointMapped& lpm ) {

 // std::cout << "CoefFunctionHyst::GetVector called" << std::endl;
  if(methodType_ == VECTOR){
    Vector<Double> valX;
    ComputeXY_vec(lpm, valX, valY, true);
  } else if(methodType_ == SCALAR){
    Double valX_scal;
    Double valY_scal;
    ComputeXY(lpm, valX_scal, valY_scal);
    valY.Resize(dependCoef_->GetVecSize());
    valY.Init();
    valY[dirP_] = valY_scal;
  }
}


void CoefFunctionHyst::SetPreviousHystVals() {

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
      Xprevious_[idx] = X;
      Yprevious_[idx] = Y;
		} else if(methodType_ == VECTOR){
		  Vector<Double> X,Y;
		  ComputeXY_vec(lpm, X, Y, true);
		  //std::cout << "SetPreviousHystVals X / Y " << X.ToString() << " / " << Y.ToString() << std::endl;
      XpreviousVEC_[idx] = X;
      YpreviousVEC_[idx] = Y;
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

  Vector<Double> tmp = X - XpreviousItVEC_[idx];
  //std::cout << "Normdiff to old input: " << tmp.NormL2() << std::endl;
  if(tmp.NormL2() < tol_){
    //std::cout << "tmp.NormL2() < tol_ -> reuse last state!" << std::endl;
      //std::cout << "Input did not change since last time -> use old value" << std::endl;
      Y = YpreviousItVEC_[idx];
  } else {



  /*
   * does not work so good as in scalar case
   * Idea: X itself can be larger than xSat insider Preisach Operator
   * the only problem is the calculation on deltaMat
   * -> normalize after Y was computed
   */
  /*
    Double Xmax = 0.0;
    std::cout << "X pre norm: " << X << std::endl;
    for(UInt i = 0; i < X.GetSize(); i++){
      if(abs(X[i]) > Xmax){
        Xmax = abs(X[i]);
      }
    }
    if(Xmax > xSat_){
      X = X * xSat_/Xmax;
    }
    std::cout << "X post norm: " << X << std::endl;
  */
  /*
   * Not so good as we change the direction of the vector

    for(UInt i = 0; i < X.GetSize(); i++){
      if(X[i] > xSat_){
        X[i] = xSat_;
      } else if(X[i] < -xSat_){
        X[i] = -xSat_;
      }
    }
  */
    /*
    if(X.NormL2() > xSat_){
      X = X * xSat_/X.NormL2();
    }*/

    /*
     * We should not normalize before calling the vector preisach as large values above sat will help to
     * align the rotation operators!
     */

    /*
     * It should make no difference if X or X_norm is provided to the function as long as eX stays the same
     * Reason:
     *  X will be normalized only, if abs(X)>xSat
     *  in that case, xthres = (x/xSat)^rot will be >= 1 for all rot > 0
     *  the algorithm switches for xthres > alpha and -xthres < beta
     *  as alphamax = 1 and betamax = 1, all domains will rotate in the case  that abs(x)>xsat and in the case that abs(x) was normalized to xsat
     *  the switching will rotate all rotationoperatores to x/abs(x); here it will again make no difference if x was scaled previously
     *  finally the switching operators are updated by the value of x in direction of the rotationoperartors
     *  as all rotationoperators will be in direction of x, the value xpar = x dot rotationoperators = x*x/abs(x) = abs(x)
     *  again abs(x) = 1 will be enough to switch all operators, i.e. a previous normalization of x wont change the behaviour either
     */

  /* Test: call computeValue_vec once with X and once with X normalized to xsat/abs(x)
   * compare y afterwards
   */




  //  if(X.NormL2() >= xSat_){
  //    ////std::cout << "compute Y with X = " << X << std::endl;
  //    Y = hyst_->computeValue_vec(X, idx);
  //    Vector<Double> Xtemp = Vector<Double>(dim_);
  //    Xtemp.Add(xSat_/X.NormL2(),X);
  //    //std::cout << "compute Ytemp with Xnorm = " << Xtemp << std::endl;
  //    Vector<Double>  Ytemp = hyst_->computeValue_vec(Xtemp, idx);
  //    //std::cout << "Y = " << Y << std::endl;
  //    //std::cout << "Ytemp = " << Ytemp << std::endl;
  //
  //    /*
  //     * it really makes no difference
  //     * -> use it
  //     */
  //  } else {
  //
  //    Y = hyst_->computeValue_vec(X, idx);
  //  }

  //  for(UInt i = 0; i < X.GetSize(); i++){
  //    if(X[i] > xSat_){
  //      X[i] = xSat_;
  //    } else if(X[i] < -xSat_){
  //      X[i] = -xSat_;
  //    }
  //  }

    if(X.NormL2() > xSat_){
      //std::cout << "Normalized X to xSat" << std::endl;
      X = X*xSat_/X.NormL2();
    }

  //  std::cout << "Try to compute value" << std::endl;
    Y = hyst_->computeValue_vec(X, idx);

  //
  //  /*
  //   * the following normalization does not work
  //*/
  //
  //  Double Xmax = 0.0;
  //  for(UInt i = 0; i < X.GetSize(); i++){
  //    if(abs(X[i]) > Xmax){
  //      Xmax = abs(X[i]);
  //    }
  //  }
  ///*
  //  if(Xmax > xSat_){
  //    X = X * xSat_/Xmax;
  //  }
  //*/
  ///*
  // * it converges also with this method
  // * still result look bad
  //*/
  //
  //  // get value of X in direction of Y
  //  Vector<Double> eY = Vector<Double>(dim_);
  //  eY.Init();
  //  eY.Add(1.0/Y.NormL2(),Y);
  //  Double valX = X * eY;
  //
  //  //std::cout << "Xmax: " << Xmax << std::endl;
  //  //std::cout << "abs(X): " << X.NormL2() << std::endl;
  //  //std::cout << "X cdot eY: " << valX << std::endl;
  //
  //  /* THIS DOES NOT WORK YET!!!! */
  //  /* TODO: FIX! */
  //
  //
  //   if(abs(valX) > xSat_){
  //     X = X * xSat_/abs(valX);
  //   }
  //  /*
  //   if(X.NormL2() > xSat_){
  //     X = X * xSat_/X.NormL2();
  //   }
  //*/
  //  /*
  // * it converges with this normaiization below; however results look not so nice
  //  */
  ///*
  //  for(UInt i = 0; i < X.GetSize(); i++){
  //     if(X[i] > xSat_){
  //       X[i] = xSat_;
  //     } else if(X[i] < -xSat_){
  //       X[i] = -xSat_;
  //     }
  //   }
  //*/

    //std::cout << "X: " << X << std::endl;
    //std::cout << "Y: " << Y << std::endl;


    /*
     * store X and Y for future iterations
     * NOTE: we explicitly store the cutted version of X as this is the quantity that matters during the
     * evaluation of the hysteresis operator!
     */
    XpreviousItVEC_[idx] = X;
    YpreviousItVEC_[idx] = Y;
//    std::cout << "X: " << X << std::endl;
//    std::cout << "Y: " << Y << std::endl;


    static UInt cnt = 0;
    static UInt firstIdx = idx;
//    std::cout << "cnt: " << cnt << std::endl;
//    std::cout << "firstIdx: " << firstIdx << std::endl;
//    std::cout << "idx: " << idx << std::endl;

    //std::cout << "CNT: " << cnt << std::endl;
    if((cnt==100)&&(idx==firstIdx)){
      if((evalVersion_ == 0)||(evalVersion_ == 6)||(evalVersion_ == 7)){
        std::cout << "Outputting bmp" << std::endl;
        std::stringstream filenamebuf;
        filenamebuf << "Switch_Elem"<<firstIdx<<"_Step" << std::setfill('0') << std::setw(5)<<cnt<<"_v"<<evalVersion_<<"_numRows"<<numRows_<<".bmp";
        hyst_->switchingStateToBmp(10000,filenamebuf.str(),idx,true);
      } else {
        std::cout << "No output of bmp as these versions do only partially update the switching matrix." << std::endl;
      }

//
//      if(evalVersion_ == 7){
//        std::stringstream filenamebuf;
//        filenamebuf << "Rot_Elem"<<firstIdx<<"_Step" << std::setfill('0') << std::setw(5)<<cnt<<"_v"<<evalVersion_<<"_numRows"<<numRows_<<".bmp";
//        hyst_->rotationStateToBmp(1000,filenamebuf.str(),idx);
//      }

    }

    if(idx == firstIdx){
      //std::cout << "cnt++" << std::endl;
      cnt++;
    }
  }
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

	UInt idx = globalElem2Local_[el->elemNum];
	X = elemSol[dirP_];

	if(X == XpreviousIt_[idx]){
	  //std::cout << "Input did not change since last time -> use old value" << std::endl;
	  Y = YpreviousIt_[idx];
	} else {
	  /*
	   * if input is too large (e.g. at electrode corners) deltaX may be huge
	   * whereas deltaY might still be small; the resulting deltaY/deltaX seems to cause trouble
	   * with adapted values, scalar model using delta material converges (at least better)
	   */
	  if(X > xSat_){
	    X = xSat_;
	  } else if(X < -xSat_){
	    X = -xSat_;
	  }

	  Y = hyst_->computeValueAndUpdate(X, idx);

	  /*
	   * store X and Y for future iterations
	   * NOTE: we explicitly store the cutted version of X as this is the quantity that matters during the
	   * evaluation of the hysteresis operator!
	   */
	  XpreviousIt_[idx] = X;
	  YpreviousIt_[idx] = Y;
	  //std::cout << "X: " << X << std::endl;
	  //std::cout << "Y: " << Y << std::endl;

	}
}

std::string CoefFunctionHyst::ToString() const {
  return "CoefFunctionHyst";
}

void CoefFunctionHyst::CreateDeltaMatrix(Vector<Double>& dX,Vector<Double>& dY, Matrix<Double>& deltaMat, UInt version, UInt idxElem){

  //std::cout << "Using version: " << version << std::endl;

    /*
     * Check vector norms first
     */
    Double norm_dX = dX.NormL2();
    Double norm_dY = dY.NormL2();

//    if((dX.NormL2() == 0) || (dY.NormL2() == 0)){
//      /*
//       * return rev_mat_fac_ on diagonal
//       */
//      for(UInt i = 0; i < dim_; i++){
//        deltaMat[i][i] = rev_mat_fac_;
//      }
//    } else
    if(version == 0){
      /*
       * std procedure, which is also used in scalar case
       *
       */
      for ( UInt i=0; i<deltaMat.GetNumRows(); i++){

          if ( (abs(dY[i]) < 1e-18) || (abs(dX[i]) < 1e-18) ) {
                //to be discussed!!
            //deltaMat = matTensorStart_;
                //tensor = matDeltaTensor_;
                deltaMat[i][i] = matDeltaTensor_[idxElem][i][i]; // MUCH BETTER THAN rev_mat_fac_;
          }
          else {
            //compute new one
            //currently we only consider detlaP / deltaE instead of deltaD / deltaE
            /*
             * works only if abs(dY/dX) is used
             * WHY???
             */
            deltaMat[i][i] =  abs( dY[i] / dX[i] );
          }
      }
    } else if(version == 1){

      if((dX.NormL2()  == 0) || (dY.NormL2() < rev_mat_fac_)){
       // std::cout << "Case1" << std::endl;
            //to be discussed!!
        //deltaMat = matTensorStart_;
            //tensor = matDeltaTensor_;
              for(UInt i = 0; i < dim_; i++){
                deltaMat[i][i] = matDeltaTensor_[idxElem][i][i]; //rev_mat_fac_;
              }
      } else {
       /*
        * solve least min problem to solve underdetermined system
        *  Ax = b
        *
        * with A = mat(dX) being:
        * 2d: [[dx1   0 dx2]
        *      [  0 dx2 dx1]]
        * 3d: [[dx1   0   0   0  dx3 dx2]
        *      [  0 dx2   0  dx3   0 dx1]
        *      [  0   0 dx3  dx2 dx1   0]]
        *
        * x = vec(deltaMat) - using Voigts notation:
        * 2d: [dm11 dm22 dm21]
        * 3d: [dm11 dm22 dm33 dm32 dm31 dm21]
        *
        * b = vec(dY):
        * 2d: [dy1 dy2]
        * 3d: [dy1 dy2 dy3]
        *
        * there are two ways to solve the system (which should be equal from mathematical point of view
        * but due to calculations of inverse matrices and so on we get quite different results)
        *
        * 1. solve
        *  A^T A x = A^T b
        *        x = (A^T A)^-1 A^T b
        * 2. solve
        *  A A^T z = b
        *        x = A^T z = A^T (A A^T)^-1 b
        *
        * (equivalence: (A^T A)^-1 A^T b ?= A^T (A A^T)^-1 b
        *   1. multiply both sides with (A^T A)
        *     (A^T A) (A^T A)^-1 A^T b = (A^T A) A^T (A A^T)^-1 b
        *     A^T b = (A^T A) A^T (A A^T)^-1 b
        *   2. reoder brackets on rhs
        *     A^T b = A^T (A A^T) (A A^T)^-1 b
        *     A^T b = A^T b
        *
        * comparison of 1. and 2.:
        *  i) if all dxi != 0, the resulting matrix from 1. is close to the one obtained by the std procedure, i.e.
        *     dmii = dyi/dxi; dmij = 0 for i!=j
        *  ii) if one dxi == 0, the matrix A^T A is singular, i.e. 1. does not work
        *  iii) 1. can be regularized by adding omega*I to A^T A before inverting; the resulting solution will be close
        *     to the one of the second method, does not solve the original system anymore but does work even if one dx1 == 0
        *  iv) 2. is solvable even if one dx1 == 0 as A A^T does not become singular
        *  v) 2. results in matrices which are quite small on the diagonal , thus leading to problems during the solution
        *     of the FE system (at least it didn't work so far)
        *  vi) 1. is harder to solve as A^T A is larger than A A^T
        *
        * suggested procedure:
        *   use regularized version of 1. i.e. solve
        *   (A^T A + omega I) x = A^T b
        *
        *   with
        *     omega = 0, if all dxi != 0
        *     omega > 0  if one dxi = 0
        *
        *   (for all dxi = 0, set dmii = rev_mat_fac; dmij = 0 for i!=j)
        *
        * -> does not work; A^T A has extrem bad condition
        * -> inversion/solution fails often
        * -> even matrix vector multiplication leads to extremly unstable results (the xth entry behind the . can
        *     change the resulting vector from e^-6 to e^-2 etc)
        *
        * new procedure:
        *   use std procedure, i.e. dmii = dyi/dxi as long as dxi > threshold
        *   if all dxi < threshold -> use matrix with rev_mat_fac as diagonal entries
        *   otherwise use 2.
        *
        *
        */

        /*
         * I create matrix A from input dX
         */
        Matrix<Double> A;
        UInt dim2 = 0;
        if(dim_ == 2){
          // number of unknowns of the material tensor
          dim2 = 3;
          /*
           *  x = [d11 d22 d21] (Voigt)
           *  A = [[dX1 0   dX2]
           *       [  0 dX2 dX1]]
           *  y = [dY1 dY2]
           */
          A = Matrix<Double>(dim_,dim2);
          A.Init();
          A[0][0] = dX[0];
          A[1][1] = dX[1];
          A[0][2] = dX[1];
          A[1][2] = dX[0];

        } else if (dim_ == 3) {
          // number of unknowns of the material tensor
          dim2 = 6;
          /*
           * x = [d11 d22 d33 d32 d31 d21] (Voigt)
           * A = [[dX1   0   0   0 dX3 dX2]
           *      [  0 dX2   0 dX3   0 dX1]
           *      [  0   0 dX3 dX2 dX1   0]]
           * y = [dY1 dY2 dY3]
           */
          A = Matrix<Double>(dim2,dim_);
          A.Init();
          A[0][0] = dX[0];
          A[0][4] = dX[2];
          A[0][5] = dX[1];
          A[1][1] = dX[1];
          A[1][3] = dX[2];
          A[1][5] = dX[0];
          A[2][2] = dX[2];
          A[2][3] = dX[1];
          A[2][4] = dX[0];
        }

        /*
         * Get matrix A^T A and rhs A^T b
         */

        Matrix<Double> ATA = Matrix<Double>(dim2,dim2);
        ATA.Init();
        Vector<Double> ATb = Vector<Double>(dim2);
        ATb.Init();

        A.MultT(A,ATA);
        A.MultT(dY,ATb);

        /*
         * Check if regularization is needed
         */
        bool needed = false;
        Double omega = rev_mat_fac_;
        for(UInt i = 0; i < dim_; i++){
          if(dX[i] == 0){
            needed = true;
          }
        }
        if(needed){
          std::cout << "Adding regularization" << std::endl;
          for(UInt i = 0; i < dim2; i++){
            ATA[i][i] += omega;
          }
        }

        std::cout << "ATA: " << ATA.ToString() << std::endl;
        std::cout << "ATb: " << ATb.ToString() << std::endl;



        /*
         * Solve system via Gaussian eleminiation -> implementation without pivoting, thus not stable
         * try implementation of invert instead -> even worse
         */

        Vector<Double> x = Vector<Double>(dim2);
        Matrix<Double> ATAI = Matrix<Double>(dim2,dim2);
        ATAI.Init();
        ATA.Invert(ATAI);

        std::cout << "ATAI: " << ATAI.ToString() << std::endl;

        x.Init();
        x = ATAI*ATb;
        //ATA.DirectSolve(x,ATb);

        std::cout << "Vector x: " << x << std::endl;

        x.Init();
        ATAI.Mult(ATb,x);

        std::cout << "Vector x (after MULT): " << x << std::endl;

        Double scal = 0.0;
        dX.Inner(dY,scal);
        std::cout << "Result of scalar product: " << scal << std::endl;


        /*
         * Distribute entries of x into deltaMat
         */
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

//
//    //    std::cout << "Case2" << std::endl;
//      /*
//       * Solve least square problem to determine delta matrix [deltaMat], such that
//       * dY = [deltaMat]*dX
//       *
//       * Assumptions: deltaMat = symmetric
//       * Procedure (for 2D, uncoupled case):
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
//       *
//       * TODO: extend for coupled case
//       */
//
//      Matrix<Double> AT;
//      UInt dim2 = 0;
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
//        AT[0][0] = dX[0];
//        AT[1][1] = dX[1];
//        AT[2][0] = dX[1];
//        AT[2][1] = dX[0];
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
//        AT[0][0] = dX[0];
//        AT[1][1] = dX[1];
//        AT[2][2] = dX[2];
//        AT[3][1] = dX[2];
//        AT[3][2] = dX[1];
//        AT[4][0] = dX[2];
//        AT[4][2] = dX[0];
//        AT[5][0] = dX[1];
//        AT[5][1] = dX[0];
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
//      //std::cout << "AAT (pre reg): " << AAT << std::endl;
//
//      // add regualrization
//      double omega = 0.0;
//      for(UInt i = 0; i < dim_; i++){
//        AAT[i][i] += omega;
//      }
//
//      //std::cout << "AAT (post reg): " << AAT << std::endl;
//
//      AAT.Invert(AATI);
//
//      Vector<Double> z = Vector<Double>(dim_);
//      Vector<Double> x = Vector<Double>(dim2);
//
//      // solve least squares problem
//      z = AATI*dY;
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
//
//      // check result
//      Vector<Double> tmp;
//      tmp = deltaMat*dX;
//
////      for(UInt i = 0; i < dim_; i++){
////        //for(UInt j = 0; j < dim_; j++){
////            deltaMat[i][i] = abs(deltaMat[i][i]);
////        //}
////      }
//
//      for(UInt i = 0; i < dim_; i++){
//        for(UInt j = 0; j < dim_; j++){
//          if(i != j){
//            /*
//             * it converges but only if offdiagonal entries are 0
//             * what is the problem with the offdiagonal entries????
//             */
//           // deltaMat[i][j] = 0; //abs(deltaMat[i][j]);
//
//            //deltaMat[i][j] = -deltaMat[i][j];
//          }
//        }
//      }

//      std::cout << "dX: " << dX << std::endl;
//      std::cout << "dY: " << dY << std::endl;
//      std::cout << "deltaMat: " << deltaMat << std::endl;
//      std::cout << "tmp - dY: " << tmp-dY << std::endl;
      //std::cout << "deltaMat: " << deltaMat << std::endl;
      }
    } else if(version == 5){
      /*
       * check if std dx/dy gives useful results
       * the definition of useful is not clear by now
       * it appears, that the method does converge if we use abs(dy/dx) but notFSetPr
       *  if dy/dx is used
       * however it may very likely be the case that abs(dy/dx) != dy/dx as the polarization may force dy into
       * a different direction than dx.
       * in that case dy/dx becomes negative and method converges to strange values
       * this may not be because of the negative sign but maybe due to the large value which results from a
       * relative large dy divided by a small dx (although the magnitude is no problem if abs() is used)
       *
       * Idea for version 4:
       *  used dy/dx as long as the quotient can be calculated (dx!=0) and the value is positive
       *  otherwise use the least squares approach
       */
      for ( UInt i=0; i<deltaMat.GetNumRows(); i++){

        if ((abs(dX[i]) < 1e-10)||(abs(dY[i])<rev_mat_fac_) ) {
          deltaMat[i][i] = rev_mat_fac_;
        } else {
          deltaMat[i][i] = abs(dY[i]/dX[i]);
        }
      }
//
//      bool working = true;
//      bool zeroInput = true;
//      //deltaMat.Init();
//      for ( UInt i=0; i<deltaMat.GetNumRows(); i++){
//
//              if ((abs(dX[i]) == 0) ) {
//                working = false;CreateDeltaMatrix
//              } else {
//                zeroInput = false;
////
////                if (dY[i]/dX[i] < 0){
////                  working = false;
////                }
//              }
//      }
//      if(zeroInput){
//        /*
//         * use rev_mat_factor
//         */
////        for ( UInt i=0; i<deltaMat.GetNumRows(); i++){
////          //std::cout << "Using delta0" << std::endl;
////          deltaMat[i][i] = rev_mat_fac_;
////        }
//        deltaMat = matDeltaTensor_[idxElem];
//      } else {
//
//        if(working){
//          //std::cout << "Using std" << std::endl;
//          for ( UInt i=0; i<deltaMat.GetNumRows(); i++){
//
//            deltaMat[i][i] = dY[i]/dX[i];
//          }
//        } else {
//          /*
//           * take old values
//           */
//          //std::cout << "Take old value" << std::endl;
//          //std::cout << matDeltaTensor_[idxElem] << std::endl;
//          deltaMat = matDeltaTensor_[idxElem];
//        }
//      }

    } else if(version == 4){
      /*
       * check if std dx/dy gives useful results
       * the definition of useful is not clear by now
       * it appears, that the method does converge if we use abs(dy/dx) but not if dy/dx is used
       * however it may very likely be the case that abs(dy/dx) != dy/dx as the polarization may force dy into
       * a different direction than dx.
       * in that case dy/dx becomes negative and method converges to strange values
       * this may not be because of the negative sign but maybe due to the large value which results from a
       * relative large dy divided by a small dx (although the magnitude is no problem if abs() is used)
       *
       * Idea for version 4:
       *  used dy/dx as long as the quotient can be calculated (dx!=0) and the value is positive
       *  otherwise use the least squares approach
       */
      bool working = true;
      bool zeroInput = true;

      for ( UInt i=0; i<deltaMat.GetNumRows(); i++){

              if ((abs(dX[i]) == 0) ) {
                working = false;
              } else {
                zeroInput = false;

                if (dY[i]/dX[i] < 0){
                  working = false;
                }
              }
      }
      if(zeroInput){
        /*
         * use rev_mat_factor
         */
        for ( UInt i=0; i<deltaMat.GetNumRows(); i++){
          //std::cout << "Using delta0" << std::endl;
          deltaMat[i][i] = matDeltaTensor_[idxElem][i][i]; //rev_mat_fac_;
        }
      } else {

      if(working){
        //std::cout << "Using std" << std::endl;
        for ( UInt i=0; i<deltaMat.GetNumRows(); i++){

          deltaMat[i][i] = dY[i]/dX[i];
        }
      } else {
//        std::cout << "Using rotation" << std::endl;
//
//        Vector<Double> dX_ext;
//        Vector<Double> dY_ext;
//        if(dim_ == 2){
//          /*
//           * Algorithm only works for 3D Vectors, so in 2D add 0 at the end
//           */
//          dX_ext = Vector<Double>(3);
//          dX_ext[0] = dX[0];
//          dX_ext[1] = dX[1];
//          dY_ext = Vector<Double>(3);
//          dY_ext[0] = dY[0];
//          dY_ext[1] = dY[1];
//        } else {
//          dX_ext = dX;
//          dY_ext = dY;
//        }
//
//        Vector<Double> tmp;
//        dX_ext.CrossProduct(dY_ext,tmp);
//
//        dX_ext = dX_ext/norm_dX;
//        dY_ext = dY_ext/norm_dY;
//
//        /*
//         * Find vector
//         *  here we have to check three vectors in the worst case
//         *  note that the three vectors are arbitrary as long as they are independent of each other
//         */
//
//        Vector<Double> e1 = Vector<Double>(3);
//        Vector<Double> e2 = Vector<Double>(3);
//        Vector<Double> e3 = Vector<Double>(3);
//
//        /*
//         * test the three unit vectors [1,0,0], [0,1,0] and [0,0,1]
//         */
//        e1[0] = 1.0;
//        //e1[1] = 1.0;
//        //e1[2] = 1.0;
//        //e1 = e1/e1.NormL2();
//        e2[1] = 1.0;
//        e3[2] = 1.0;
//
//        bool found = false;
//        Vector<Double> v;
//        Vector<Double> vRotdX;
//        Vector<Double> vRotdY;
//
//        for(UInt i = 0; i < 3; i++){
//
//          if(i == 0){
//            v = e1;
//          } else if(i == 1){
//            v = e2;
//          } else {
//            v = e3;
//          }
//
//          v.CrossProduct(dX_ext,vRotdX);
//    //      std::cout << "v: " << v.ToString() << std::endl;
//    //      std::cout << "dX_ext: " << dX_ext.ToString() << std::endl;
//    //      std::cout << "vRotdX: " << vRotdX.ToString() << std::endl;
//          v.CrossProduct(dY_ext,vRotdY);
//
//          if(vRotdX.NormL2() != 0 && vRotdY.NormL2() != 0){
//            vRotdX = vRotdX/vRotdX.NormL2();
//            vRotdY = vRotdY/vRotdY.NormL2();
//            found = true;
//          //  std::cout << "Using: " << v << std::endl;
//            break;
//          }
//        }
//
//        if(found == false){
//          EXCEPTION("No suitable vector found. Check input!");
//          /*
//           * should only happen if one of the input vectors in 0,0,0
//           */
//        }
//
//        /*
//         * Create Matrices
//         */
//        Vector<Double> m1 = vRotdX;
//        Vector<Double> m2;
//        dX_ext.CrossProduct(vRotdX,m2);
//        Vector<Double> m3;
//        m1.CrossProduct(m2,m3);
//
//        Vector<Double> n1 = vRotdY;
//        Vector<Double> n2;
//        dY_ext.CrossProduct(vRotdY,n2);
//        Vector<Double> n3;
//        n1.CrossProduct(n2,n3);
//
//    //    std::cout << "m1: " << m1.ToString() << std::endl;
//    //    std::cout << "m1.NormL2()" << m1.NormL2() << std::endl;
//    //    std::cout << "m2: " << m2.ToString() << std::endl;
//    //    std::cout << "m2.NormL2()" << m2.NormL2() << std::endl;
//    //    std::cout << "m3: " << m3.ToString() << std::endl;
//    //    std::cout << "m3.NormL2()" << m3.NormL2() << std::endl;
//
//
//        Vector<Double> tmp1, tmp2;
//        Matrix<Double> M = Matrix<Double>(3,3);
//        Matrix<Double> N = Matrix<Double>(3,3);
//        for(UInt i = 0; i < 3; i++){
//          for(UInt j = 0; j < 3; j++){
//            if(i == 0){
//              tmp1 = m1;
//              tmp2 = n1;
//            } else if(i == 1){
//              tmp1 = m2;
//              tmp2 = n2;
//            } else {
//              tmp1 = m3;
//              tmp2 = n3;
//            }
//            /*
//             * mi are the columns of M
//             * ni are the columns of N
//             */
//            M[j][i] = tmp1[j];
//            N[j][i] = tmp2[j];
//          }
//        }
//
//        //std::cout << "M: " << M.ToString() << std::endl;
//
//        /*
//         * Create rotation matrix
//         */
//        Matrix<Double> M_inv;
//        M.Invert(M_inv);
//
//        Matrix<Double> R = Matrix<Double>(3,3);
//        N.Mult(M_inv,R);
//
//        if(dim_ == 2){
//          /*
//           * take only upper 2 x 2 submatrix
//           */
//          Matrix<Double> R_tmp = R;
//          R.Resize(2,2);
//          R[0][0] = R_tmp[0][0];
//          R[0][1] = R_tmp[0][1];
//          R[1][0] = R_tmp[1][0];
//          R[1][1] = R_tmp[1][1];
//
//    //      std::cout << "R_tmp: " << R_tmp.ToString() << std::endl;
//    //      std::cout << "R: " << R.ToString() << std::endl;
//        }
//
//        /*
//         * Scale Matrix to correct amplitudes
//         */
//        Double fact = norm_dY/norm_dX;
//        R = R*fact;
//
//  //      /*
//  //       * Check for correctness
//  //       */
//        Vector<Double> test = Vector<Double>(dim_);
//        R.Mult(dX,test);
//  //      std::cout << "R: " << R.ToString() << std::endl;
//  //      std::cout << "dX: " << dX.ToString() << std::endl;
//  //      std::cout << "dY: " << dY.ToString() << std::endl;
//  //      std::cout << "test: " << test.ToString() << std::endl;
//
//        deltaMat = R;
//
//        /*
//         * it converges but only if offdiagonal entries are 0
//         * what is the problem with the offdiagonal entries????
//         */
//        for(UInt i = 0; i < dim_; i++){
//          for(UInt j = 0; j < dim_; j++){
//            if(i != j){
//              /*
//               * it converges but only if offdiagonal entries are 0
//               * what is the problem with the offdiagonal entries????
//               */
//            //  deltaMat[i][j] = 0; //abs(deltaMat[i][j]);
//
//            //  deltaMat[i][j] = 0; //-deltaMat[i][j];
//            }
//          }
//        }


        //std::cout << "Using least squares" << std::endl;
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
          AT[0][0] = dX[0];
          AT[1][1] = dX[1];
          AT[2][0] = dX[1];
          AT[2][1] = dX[0];

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
          AT[0][0] = dX[0];
          AT[1][1] = dX[1];
          AT[2][2] = dX[2];
          AT[3][1] = dX[2];
          AT[3][2] = dX[1];
          AT[4][0] = dX[2];
          AT[4][2] = dX[0];
          AT[5][0] = dX[1];
          AT[5][1] = dX[0];
        }

        Matrix<Double> AAT = Matrix<Double>(dim_,dim_);
        /*
         * perform AAT = (AT)^T * AT
         */
        AT.MultT(AT,AAT);

        Matrix<Double> AATI = Matrix<Double>(dim_,dim_);

        AAT.Invert(AATI);

        Vector<Double> z = Vector<Double>(dim_);
        Vector<Double> x = Vector<Double>(dim2);

        // solve least squares problem
        z = AATI*dY;

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
     }
     }



    } else if(version == 2){
      if((dX.NormL2()  < 1e-10) || (dY.NormL2() < rev_mat_fac_)){
                 //to be discussed!!
             //deltaMat = matTensorStart_;
                 //tensor = matDeltaTensor_;
                   for(UInt i = 0; i < dim_; i++){
                     deltaMat[i][i] = rev_mat_fac_;
                   }
      } else {
      /*
       * create matrix R, such that dY = R*dX
       *
       * Algorithm:
       *  1. find vector v such that dY x v and dX x v both != 0; extend vector to 3D!
       *  2. create m1 = dX x v; m2 = v x m1; m3 = m1 x m2
       *            n1 = dY x v; n2 = v x n1; n3 = n1 x n2
       *     note: vectors should be normalized to length 1
       *  3. create matrices M = [m1; m2; m3] and N = [n1; n2; n3]
       *  4. R = N*inv(M)*norm(dY)/norm(dX)
       * (5. check if dY = R*dX
       *
       * -> matrix is not symmetric!
       */
      /*
       * extend vectors to 3D such that we can evaluate the cross product
       */

      Vector<Double> dX_ext;
      Vector<Double> dY_ext;
      if(dim_ == 2){
        /*
         * Algorithm only works for 3D Vectors, so in 2D add 0 at the end
         */
        dX_ext = Vector<Double>(3);
        dX_ext[0] = dX[0];
        dX_ext[1] = dX[1];
        dY_ext = Vector<Double>(3);
        dY_ext[0] = dY[0];
        dY_ext[1] = dY[1];
      } else {
        dX_ext = dX;
        dY_ext = dY;
      }

      Vector<Double> tmp;
      dX_ext.CrossProduct(dY_ext,tmp);

      //std::cout << "tmp.NormL2: " << tmp.NormL2() << std::endl;

      bool use_rot = true;
      for(UInt i = 0; i < dim_; i++){
        if((abs(dX[i]) == 0 ) && (dY[i] != 0)){
          use_rot = true;
        }
      }

      if(use_rot == false){

        for(UInt i = 0; i < dim_; i++){
          if(abs(dX[i]) > 1e-10){
            deltaMat[i][i] = abs(dY[i]/dX[i]);
          } else {
//            if(i == 0){
//              deltaMat[i][i] = abs(dY[1]/dX[1]);
//            } else {
//              deltaMat[i][i] = abs(dY[0]/dX[0]);
//            }
         //   std::cout << "Taking old value" << std::endl;
         //   std::cout << "dX = " << dX.ToString() << std::endl;
            //deltaMat[i][i] = rev_mat_fac_;
          }

        }
      } else {

       // std::cout << "Taking full rot" << std::endl;

        /*
         * normalize vectors
         */
        dX_ext = dX_ext/norm_dX;
        dY_ext = dY_ext/norm_dY;

        /*
         * Find vector
         *  here we have to check three vectors in the worst case
         *  note that the three vectors are arbitrary as long as they are independent of each other
         */

        Vector<Double> e1 = Vector<Double>(3);
        Vector<Double> e2 = Vector<Double>(3);
        Vector<Double> e3 = Vector<Double>(3);

        /*
         * test the three unit vectors [1,0,0], [0,1,0] and [0,0,1]
         */
        e1[0] = 1.0;
        //e1[1] = 1.0;
        //e1[2] = 1.0;
        //e1 = e1/e1.NormL2();
        e2[1] = 1.0;
        e3[2] = 1.0;

        bool found = false;
        Vector<Double> v;
        Vector<Double> vRotdX;
        Vector<Double> vRotdY;

        for(UInt i = 0; i < 3; i++){

          if(i == 0){
            v = e1;
          } else if(i == 1){
            v = e2;
          } else {
            v = e3;
          }

          v.CrossProduct(dX_ext,vRotdX);
    //      std::cout << "v: " << v.ToString() << std::endl;
    //      std::cout << "dX_ext: " << dX_ext.ToString() << std::endl;
    //      std::cout << "vRotdX: " << vRotdX.ToString() << std::endl;
          v.CrossProduct(dY_ext,vRotdY);

          if(vRotdX.NormL2() != 0 && vRotdY.NormL2() != 0){
            vRotdX = vRotdX/vRotdX.NormL2();
            vRotdY = vRotdY/vRotdY.NormL2();
            found = true;
          //  std::cout << "Using: " << v << std::endl;
            break;
          }
        }

        if(found == false){
          EXCEPTION("No suitable vector found. Check input!");
          /*
           * should only happen if one of the input vectors in 0,0,0
           */
        }

        /*
         * Create Matrices
         */
        Vector<Double> m1 = vRotdX;
        Vector<Double> m2;
        dX_ext.CrossProduct(vRotdX,m2);
        Vector<Double> m3;
        m1.CrossProduct(m2,m3);

        Vector<Double> n1 = vRotdY;
        Vector<Double> n2;
        dY_ext.CrossProduct(vRotdY,n2);
        Vector<Double> n3;
        n1.CrossProduct(n2,n3);

    //    std::cout << "m1: " << m1.ToString() << std::endl;
    //    std::cout << "m1.NormL2()" << m1.NormL2() << std::endl;
    //    std::cout << "m2: " << m2.ToString() << std::endl;
    //    std::cout << "m2.NormL2()" << m2.NormL2() << std::endl;
    //    std::cout << "m3: " << m3.ToString() << std::endl;
    //    std::cout << "m3.NormL2()" << m3.NormL2() << std::endl;


        Vector<Double> tmp1, tmp2;
        Matrix<Double> M = Matrix<Double>(3,3);
        Matrix<Double> N = Matrix<Double>(3,3);
        for(UInt i = 0; i < 3; i++){
          for(UInt j = 0; j < 3; j++){
            if(i == 0){
              tmp1 = m1;
              tmp2 = n1;
            } else if(i == 1){
              tmp1 = m2;
              tmp2 = n2;
            } else {
              tmp1 = m3;
              tmp2 = n3;
            }
            /*
             * mi are the columns of M
             * ni are the columns of N
             */
            M[j][i] = tmp1[j];
            N[j][i] = tmp2[j];
          }
        }

        //std::cout << "M: " << M.ToString() << std::endl;

        /*
         * Create rotation matrix
         */
        Matrix<Double> M_inv;
        M.Invert(M_inv);

        Matrix<Double> R = Matrix<Double>(3,3);
        N.Mult(M_inv,R);

        if(dim_ == 2){
          /*
           * take only upper 2 x 2 submatrix
           */
          Matrix<Double> R_tmp = R;
          R.Resize(2,2);
          R[0][0] = R_tmp[0][0];
          R[0][1] = R_tmp[0][1];
          R[1][0] = R_tmp[1][0];
          R[1][1] = R_tmp[1][1];

    //      std::cout << "R_tmp: " << R_tmp.ToString() << std::endl;
    //      std::cout << "R: " << R.ToString() << std::endl;
        }

        /*
         * Scale Matrix to correct amplitudes
         */
        Double fact = norm_dY/norm_dX;
        R = R*fact;

  //      /*
  //       * Check for correctness
  //       */
        Vector<Double> test = Vector<Double>(dim_);
        R.Mult(dX,test);
  //      std::cout << "R: " << R.ToString() << std::endl;
  //      std::cout << "dX: " << dX.ToString() << std::endl;
  //      std::cout << "dY: " << dY.ToString() << std::endl;
  //      std::cout << "test: " << test.ToString() << std::endl;

        deltaMat = R;

        /*
         * it converges but only if offdiagonal entries are 0
         * what is the problem with the offdiagonal entries????
         */
        for(UInt i = 0; i < dim_; i++){
          for(UInt j = 0; j < dim_; j++){
            if(i != j){
              /*
               * it converges but only if offdiagonal entries are 0
               * what is the problem with the offdiagonal entries????
               */
              deltaMat[i][j] = 0; //abs(deltaMat[i][j]);

              //deltaMat[i][j] = -deltaMat[i][j];
            }
          }
        }

        //std::cout << "deltaMat: " << deltaMat << std::endl;
       }
     }
    }
   }

}// namespace
