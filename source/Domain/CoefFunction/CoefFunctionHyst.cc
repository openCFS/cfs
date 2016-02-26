#include "CoefFunctionHyst.hh"

// classes for function / spline approximation
#include "Materials/Models/Preisach.hh"
#include "Materials/Models/VectorPreisach.hh"
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

	  std::cout << "Using Scalar Hystersis" << std::endl;

    //get direction
    std::string str;
    material->GetScalar(str, P_DIRECTION);
    Directions dir;
    String2Enum(str,dir);
    dirP_ = dir;

    hyst_ = new Preisach(numElemSD, Xsat, Ysat, weights, isVirgin);

    dim_ = 1;
	} else if(methodType_ == VECTOR){

	  std::cout << "Using Vector Hystersis" << std::endl;

	  material->GetScalar(rotRes_, ROT_RESISTANCE, Global::REAL);


	  dim_ = dependCoef_->GetVecSize();
	  hyst_ = new VectorPreisach(numElemSD, Xsat, Ysat, weights,rotRes_,dim_, isVirgin);
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
    Xprevious_.Resize(numElemSD);
    Yprevious_.Resize(numElemSD);
    Xprevious_.Init();
    Yprevious_.Init();
	} else if(methodType_ == VECTOR){
	  XpreviousVEC_ = new Vector<Double>[numElemSD];
	  YpreviousVEC_ = new Vector<Double>[numElemSD];

	  for(UInt i = 0; i < numElemSD; i++){
	    XpreviousVEC_[i].Resize(dim_);
	    XpreviousVEC_[i].Init();

      YpreviousVEC_[i].Resize(dim_);
      YpreviousVEC_[i].Init();
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
	matDeltaTensor_ = matTensorStart_;
	//std::cout << "Tensor:\n " << matTensorStart_ << std::endl;
	
	  // to calculate differential materialproperties, we need to know e0/ mu0
  if(material_->GetMaterialDatabaseName() == "Electrostatic"){
    rev_mat_fac_ = 8.854187817e-12;
  } else if(material_->GetName() == "Magnetic"){
    rev_mat_fac_ = 1.2566370614e-6;
  } else {
  }

}


void CoefFunctionHyst::GetTensor( Matrix<Double>& tensor,
                                  const LocPointMapped& lpm ) {


 // std::cout << "CoefFunctionHyst::GetTensor called" << std::endl;
  const Elem * el = lpm.ptEl;
  UInt idx = globalElem2Local_[el->elemNum];


  if(methodType_ == SCALAR){

    Double Xval, Ycurrent, matDiff;
    ComputeXY(lpm, Xval, Ycurrent);

    //std::cout << "E=" << Xval << "P=" << Ycurrent  << std::endl;

    //compute differential material parameter
    Double dX = Xval - Xprevious_[idx];
    // D = eps0 + P
    // Ycurrent-Yprevious = deltaP!
    Double dY = Ycurrent -Yprevious_[idx] + rev_mat_fac_*dX;
//    std::cout << "dY: " << dY << " / dX: " << dX << std::endl;
//    std::cout << "(dY/dX)/eps0 = " << dY/dX/rev_mat_fac_ << std::endl;
    if ( (abs(dY) < rev_mat_fac_) || (abs(dX) < 1e-10) ) {
      //to be discussed!!
      tensor = matTensorStart_;
      //tensor = matDeltaTensor_;
    }
    else {
      //compute new one
      //currently we only consider detlaP / deltaE instead of deltaD / deltaE
      matDiff = abs(dY) / abs(dX);
      for ( UInt i=0; i<matDeltaTensor_.GetNumRows(); i++){
        matDeltaTensor_[i][i] = matDiff;
      }
      //matDeltaTensor_[dirP_][dirP_] = matDiff;
      tensor = matDeltaTensor_;
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
    ComputeXY_vec(lpm, Xval, Ycurrent);

    Vector<Double> dX = Xval;
    dX.Add(-1.0,XpreviousVEC_[idx]);
    Vector<Double> dY = Ycurrent;
    dY.Add(-1.0,YpreviousVEC_[idx]);
    // add rev part
    dY.Add(rev_mat_fac_,dX);

    //matDeltaTensor_.Init();
    //CreateDeltaMatrix(dX,dY,matDeltaTensor_);

    for ( UInt i=0; i<matDeltaTensor_.GetNumRows(); i++){

        if ( (abs(dY[i]) < rev_mat_fac_) || (abs(dX[i]) < 1e-10) ) {
              //to be discussed!!
              tensor = matTensorStart_;
              //tensor = matDeltaTensor_;
        }
        else {
          //compute new one
          //currently we only consider detlaP / deltaE instead of deltaD / deltaE
          matDeltaTensor_[i][i] = abs(dY[i]) / abs(dX[i]);

          //matDeltaTensor_[dirP_][dirP_] = matDiff;
          tensor = matDeltaTensor_;
        }
    }
//    std::cout << "Resulting matrix new: " << std::endl;
//    std::cout << matDeltaTensor_.ToString() << std::endl;

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
    tensor = matDeltaTensor_;
  }

  //std::cout << "CoefFncHyst was called; returning: " << tensor << std::endl;

}//ELEC_FIELD_INTENS


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
    ComputeXY_vec(lpm, valX, valY);
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
		  ComputeXY_vec(lpm, X, Y);
		  //std::cout << "SetPreviousHystVals X / Y " << X.ToString() << " / " << Y.ToString() << std::endl;
      XpreviousVEC_[idx] = X;
      YpreviousVEC_[idx] = Y;
		}
	}
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
    X = X*xSat_/X.NormL2();
  }

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

  //std::cout << "X: " << X << std::endl;
  //std::cout << "Y: " << Y << std::endl;

}



std::string CoefFunctionHyst::ToString() const {
  return "CoefFunctionHyst";
}

void CoefFunctionHyst::CreateDeltaMatrix(Vector<Double>& dX,Vector<Double>& dY, Matrix<Double>& deltaMat){
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
     * Check vector norms
     */
    Double norm_dX = dX.NormL2();
    Double norm_dY = dY.NormL2();
    //deltaMat.Resize(dim_,dim_);

    if((dX.NormL2() == 0) || (dY.NormL2() == 0)){
      /*
       * return rev_mat_fac_ on diagonal
       */
      for(UInt i = 0; i < dim_; i++){
        deltaMat[i][i] = rev_mat_fac_;
      }
    } else {
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

      bool use_rot = false;
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
            std::cout << "Taking old value" << std::endl;
            std::cout << "dX = " << dX.ToString() << std::endl;
            //deltaMat[i][i] = rev_mat_fac_;
          }

        }

      } else {

        std::cout << "Taking full rot" << std::endl;

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
            std::cout << "Using: " << v << std::endl;
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
       }
     }
   }

}// namespace
