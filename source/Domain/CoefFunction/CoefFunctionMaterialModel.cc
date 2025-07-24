#include "DataInOut/Logging/LogConfigurator.hh"
#include "General/Environment.hh"
#include <cmath>

#include <boost/bind/bind.hpp>
#include <boost/lexical_cast.hpp>
#include <Domain/CoefFunction/CoefFunctionMaterialModel.hh>

//Include the new Models here
#include "Materials/Models/Jiles.hh"
#include "Materials/Models/Model.hh"

#include <list>

#include "MatVec/Vector.hh"
#include "CoefFunctionMaterialModel.hh"

namespace CoupledField {

DEFINE_LOG(cfjc, "CoefFunctionMaterialModel")

// ==========================================================================
//  COEFFICIENT FUNCTION CACHE
// ==========================================================================
// Example for using Coefficient Function Cache:
// See MechPDE::DefinePostProcResults() MechStress, MechPrincipalStress
template<class TYPE> CoefFunctionMaterialModel<TYPE>::CoefFunctionMaterialModel() :
    CoefFunction() {
  dependType_ = SOLUTION;

  isAnalytic_ = false;
  isComplex_ = false;
  //this should change, depending on the model
  dimType_ = SCALAR;

}

template<class TYPE> void CoefFunctionMaterialModel<TYPE>::Init( PtrCoefFct depCoef, std::string modelName, UInt dim){

  LOG_DBG(cfjc)
  << "Init CoefFunctionMaterialModel" << std::endl;

  depCoef_ = depCoef;

  modelName_ = modelName;

  spaceDim_ = dim;

  if (modelName_ == "JilesAthertonModel") {

    dimType_ = SCALAR;

    static Jiles JilesModel;
    matModel_ = &JilesModel;
  } else if(modelName_ == "EBHysteresisModel"){
    dimType_ = TENSOR;

    static EBHysteresis EBHysteresisModel;
    matModel_ = &EBHysteresisModel;
  } else if(modelName_ == "invEBHysteresisModel"){
    dimType_ = TENSOR;

    static invEBHysteresis invEBHysteresisModel;
    matModel_ = &invEBHysteresisModel;
  } else {

    EXCEPTION("Model not implemented! ("<< modelName_<<")")
  }

}

template<class TYPE> CoefFunctionMaterialModel<TYPE>::~CoefFunctionMaterialModel() {
}

template<class T> void CoefFunctionMaterialModel<T>::InitModel(
    std::map<std::string, double> ParameterMap, UInt numElems) {

  matModel_->Init(ParameterMap, numElems, spaceDim_);

}

template<class T> void CoefFunctionMaterialModel<T>::InitModel(
    std::map<std::string, double> ParameterMap, shared_ptr<ElemList> entityList) {

  matModel_->Init(ParameterMap, entityList, spaceDim_);

}

template<class T> void CoefFunctionMaterialModel<T>::InitModel(
  std::map<std::string, double> ParameterMap, std::map<std::string, string> StringParameterMap, shared_ptr<ElemList> entityList) {

matModel_->Init(ParameterMap, StringParameterMap, entityList, spaceDim_);

}


template <class TYPE>
void CoefFunctionMaterialModel<TYPE>::UpdateHistoryValues()
{
  matModel_->UpdateStates();
}

template <class TYPE>
void CoefFunctionMaterialModel<TYPE>::AllowUpdates(bool allow)
{
  matModel_->AllowUpdates(allow);
}


template <class T>
void CoefFunctionMaterialModel<T>::GetScalar(
    Double &coefScalar, const LocPointMapped &lpm)
{
    Vector<Complex> DependentVec;

    depCoef_->GetVector(DependentVec, lpm);

    // Can i do this in less codelines? Complex vector to Real Vector??
    Vector<Double> RealDependentVec;

    RealDependentVec.Resize(3);
    RealDependentVec.Init(0);

    RealDependentVec[0] = std::real(DependentVec[0]);
    RealDependentVec[1] = std::real(DependentVec[1]);
    RealDependentVec[2] = std::real(DependentVec[2]);

    coefScalar = matModel_->ComputeMaterialParameter(RealDependentVec, lpm.ptEl->elemNum);

    LOG_DBG(cfjc)
        << "NrElem = :" << lpm.ptEl->elemNum << std::endl;
    LOG_DBG(cfjc)
        << "E = :[" << RealDependentVec.ToString() << "]" << std::endl;
    LOG_DBG(cfjc)
        << "Epsilon = :" << coefScalar << std::endl;
}

template<class T> void CoefFunctionMaterialModel<T>::GetTensor(
    Matrix<Double> &coefTensor, const LocPointMapped &lpm) {
  Vector<Complex> DependentVec;

  depCoef_->GetVector(DependentVec, lpm);

  Vector<Double> RealDependentVec;

  RealDependentVec.Resize(spaceDim_);
  RealDependentVec.Init(0);

  for(UInt i = 0 ; i < spaceDim_; ++i){
    RealDependentVec[i] = std::real(DependentVec[i]); 
  }

  coefTensor = matModel_->ComputeTensorialMaterialParameter(RealDependentVec, lpm.ptEl->elemNum);

  LOG_DBG(cfjc)
  << "NrElem = :" << lpm.ptEl->elemNum << std::endl;
  LOG_DBG(cfjc)
  << "E = :[" << RealDependentVec.ToString() << "]" << std::endl;
  LOG_DBG(cfjc)
  << "Epsilon = :" << coefTensor << std::endl;
}

template<class T> void CoefFunctionMaterialModel<T>::GetVector(
    Vector<Double> &coefVector, const LocPointMapped &lpm) {
  Vector<Complex> DependentVec;
  depCoef_->GetVector(DependentVec, lpm);

  Vector<Double> RealDependentVec;

  RealDependentVec.Resize(spaceDim_);
  RealDependentVec.Init(0);

  for(UInt i = 0 ; i < spaceDim_; ++i){
    RealDependentVec[i] = std::real(DependentVec[i]); 
  }
  if (modelName_ == "EBHysteresisModel") {
    if(stressCoef_){
      // multiscale version, which requires mechanical stress input
      coefVector = matModel_->GetFluxDensity(RealDependentVec, lpm.ptEl->elemNum, lpm, stressCoef_);
    }else{
      coefVector = matModel_->GetFluxDensity(RealDependentVec, lpm.ptEl->elemNum);
    }
  } else if (modelName_ == "invEBHysteresisModel")  {
    coefVector = matModel_->GetFieldIntensity(RealDependentVec, lpm.ptEl->elemNum);
  } else {
    coefVector = matModel_->GetFluxDensity(RealDependentVec, lpm.ptEl->elemNum);
  }


  
  LOG_DBG(cfjc)
  << "NrElem = :" << lpm.ptEl->elemNum << std::endl;
  LOG_DBG(cfjc)
  << "Flux Density = :" << coefVector << std::endl;
}


template class CoefFunctionMaterialModel<Double> ;
template class CoefFunctionMaterialModel<Complex> ;

} // end of namespace
