#include "CoefFunctionHyst.hh"

// classes for function / spline approximation
#include "Materials/Models/Preisach.hh"
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

	//get direction
	std::string str;
	material->GetScalar(str, P_DIRECTION);
	Directions dir;
	String2Enum(str,dir);
	dirP_ = dir;


	bool isVirgin = true;
	UInt numElemSD = actSDList->GetSize();

	hyst_ = new Preisach(numElemSD, Xsat, Ysat, weights, isVirgin);

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
	Xprevious_.Resize(numElemSD);
	Yprevious_.Resize(numElemSD);
	Xprevious_.Init();
	Yprevious_.Init();

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

	Double Xval, Ycurrent, matDiff;;
	ComputeXY(lpm, Xval, Ycurrent);
	const Elem * el = lpm.ptEl;
	UInt idx = globalElem2Local_[el->elemNum];

	//std::cout << "E=" << Xval << "P=" << Ycurrent  << std::endl;

	//compute differential material parameter
	Double dX = Xval - Xprevious_[idx];
	Double dY = Ycurrent -Yprevious_[idx];

	if ( (abs(dY) < 1e-12) || (abs(dX) < 1e-10) ) {
		//to be discussed!!
		//tensor = matTensorStart_;
		tensor = matDeltaTensor_;
	}
	else {
		//compute new one
		matDiff = dY / dX;
		for ( UInt i=0; i<matDeltaTensor_.GetNumRows(); i++)
			matDeltaTensor_[i][i] = matDiff;
		tensor = matDeltaTensor_;
	}
}


void CoefFunctionHyst::GetScalar( Double& valY,
                                  const LocPointMapped& lpm ) {

	Double valX;
	ComputeXY(lpm, valX, valY);
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

		Double X, Y;
		ComputeXY(lpm, X, Y);
		UInt idx = globalElem2Local_[el->elemNum];
		Xprevious_[idx] = X;
		Yprevious_[idx] = Y;
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
	Y = hyst_->computeValueAndUpdate(X, idx);
}



std::string CoefFunctionHyst::ToString() const {
  return "CoefFunctionHyst";
}

} // namespace
