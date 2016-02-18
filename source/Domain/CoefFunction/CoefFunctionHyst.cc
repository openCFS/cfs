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
  dimType_ = TENSOR;
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

	bool isVirgin = true;
	UInt numElemSD = actSDList->GetSize();

	std::string dimTypeStr;
	material->GetScalar(dimTypeStr, PREISACH_DIM);

	if(dimTypeStr == "SCALAR"){
	  dimType_ = SCALAR;
	} else if(dimTypeStr == "VECTOR"){
	  dimType_ = VECTOR;
	}

	if(dimType_ == SCALAR){

	  std::cout << "Using Scalar Hystersis" << std::endl;

    //get direction
    std::string str;
    material->GetScalar(str, P_DIRECTION);
    Directions dir;
    String2Enum(str,dir);
    dirP_ = dir;

    hyst_ = new Preisach(numElemSD, Xsat, Ysat, weights, isVirgin);

    dim_ = 1;
	} else if(dimType_ == VECTOR){

	  std::cout << "Using Vector Hystersis" << std::endl;

	  Double rotationalResistance;
	  material->GetScalar(rotationalResistance, ROT_RESISTANCE, Global::REAL);


	  dim_ = dependCoef_->GetVecSize();
	  hyst_ = new VectorPreisach(numElemSD, Xsat, Ysat, weights,rotationalResistance,dim_, isVirgin);
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
	if(dimType_ == SCALAR){
    Xprevious_.Resize(numElemSD);
    Yprevious_.Resize(numElemSD);
    Xprevious_.Init();
    Yprevious_.Init();
	} else if(dimType_ == VECTOR){
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
}


void CoefFunctionHyst::GetTensor( Matrix<Double>& tensor,
                                  const LocPointMapped& lpm ) {


 // std::cout << "CoefFunctionHyst::GetTensor called" << std::endl;
  const Elem * el = lpm.ptEl;
  UInt idx = globalElem2Local_[el->elemNum];

  double rev_mat_fac = 0.0;
  // to calculate differential materialproperties, we need to know e0/ mu0
  if(material_->GetMaterialDatabaseName() == "Electrostatic"){
    rev_mat_fac = 8.854187817e-12;
  } else if(material_->GetName() == "Magnetic"){
    rev_mat_fac = 1.2566370614e-6;
  } else {
  }

  if(dimType_ == SCALAR){

    Double Xval, Ycurrent, matDiff;
    ComputeXY(lpm, Xval, Ycurrent);

    //std::cout << "E=" << Xval << "P=" << Ycurrent  << std::endl;

    //compute differential material parameter
    Double dX = Xval - Xprevious_[idx];
    // D = eps0 + P
    // Ycurrent-Yprevious = deltaP!
    Double dY = Ycurrent -Yprevious_[idx] + rev_mat_fac*dX;

    if ( (abs(dY) < 1e-12) || (abs(dX) < 1e-10) ) {
      //to be discussed!!
      //tensor = matTensorStart_;
      tensor = matDeltaTensor_;
    }
    else {
      //compute new one
      //currently we only consider detlaP / deltaE instead of deltaD / deltaE
      matDiff = dY / dX;
      for ( UInt i=0; i<matDeltaTensor_.GetNumRows(); i++)
        matDeltaTensor_[i][i] = matDiff;
      tensor = matDeltaTensor_;
    }

  } else if(dimType_ == VECTOR ){

    Vector<Double> Xval, Ycurrent;
    Double matDiff;
    ComputeXY_vec(lpm, Xval, Ycurrent);

    for(UInt i = 0; i < dim_; i++ ){

      Double dX = Xval[i] - XpreviousVEC_[idx][i];
      Double dY = Ycurrent[i] - YpreviousVEC_[idx][i] + rev_mat_fac*dX;

      if ( (abs(dY) < rev_mat_fac) || (abs(dX) < 1e-10) ) {
      //if ( (abs(dX) < 1e-10) ) {
        //to be discussed!!
        //matDeltaTensor_[i][i] = rev_mat_fac;
      } else {
        matDiff = dY/dX;
        //std::cout << "Mat diff / component " << matDiff << " / " << i << std::endl;
        matDeltaTensor_[i][i] = matDiff;
      }
    }
    tensor = matDeltaTensor_;
  }

}//ELEC_FIELD_INTENS


void CoefFunctionHyst::GetScalar( Double& valY,
                                  const LocPointMapped& lpm ) {

  if(dimType_ != SCALAR){
    EXCEPTION("Only implemented for scalar model");
  }
	Double valX;
	ComputeXY(lpm, valX, valY);
}

void CoefFunctionHyst::GetVector( Vector<Double>& valY,const LocPointMapped& lpm ) {

 // std::cout << "CoefFunctionHyst::GetVector called" << std::endl;
  if(dimType_ == VECTOR){
    Vector<Double> valX;
    ComputeXY_vec(lpm, valX, valY);
  } else if(dimType_ == SCALAR){
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

		if(dimType_ == SCALAR){
      Double X, Y;
      ComputeXY(lpm, X, Y);
      Xprevious_[idx] = X;
      Yprevious_[idx] = Y;
		} else if(dimType_ == VECTOR){
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
  Y = hyst_->computeValue_vec(X, idx);

 // std::cout << "X: " << X << std::endl;
 // std::cout << "Y: " << Y << std::endl;

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
	Y = hyst_->computeValueAndUpdate(X, idx);

  //std::cout << "X: " << X << std::endl;
  //std::cout << "Y: " << Y << std::endl;

}



std::string CoefFunctionHyst::ToString() const {
  return "CoefFunctionHyst";
}

} // namespace
