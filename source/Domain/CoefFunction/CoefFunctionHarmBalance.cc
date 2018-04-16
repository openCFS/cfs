// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <def_expl_templ_inst.hh>
#include <limits>

#include "CoefFunctionHarmBalance.hh"
#include <boost/bind.hpp>
#include <boost/lexical_cast.hpp>
#include "Domain/Mesh/Grid.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"

#include "Domain/CoefFunction/CoefXpr.hh"


namespace CoupledField
{
  // ==========================================================================
  //  COEFFICIENT FUNCTION HARMONIC BALANCING
  // ==========================================================================
  template<class T>
  CoefFunctionHarmBalance<T>::CoefFunctionHarmBalance(shared_ptr<BaseFeFunction> feFct,
                                                      const StdVector<RegionIdType>& regions,
                                                      const std::map<RegionIdType, BaseMaterial*>& materials,
                                                       Grid* ptGrid,
                                                      PtrCoefFct magFluxCoef):CoefFunction() {

    feFct_ = dynamic_pointer_cast<FeFunction<T> >(feFct);
    regions_ = regions;
    materials_ = materials;
    ptGrid_ = ptGrid;

    magFluxCoef_ = magFluxCoef;

    maxInt_ = std::numeric_limits<unsigned int>::max();

    // Reset the frequency results
    freqResult_.Resize(0,0);


    // For the callback mechanism
    mp_ = domain->GetMathParser();
    harmHandle_ = mp_->GetNewHandle(true);
    // Bind the handle to the correct expression
    mp_->SetExpr(harmHandle_,"mhFlag");
    // Flag for first time calculation
    mp_->SetValue(MathParser::GLOB_HANDLER, "mhFlag", std::numeric_limits<unsigned int>::max());
    // register callback mechanism if expression changes
    mp_->AddExpChangeCallBack( boost::bind(&CoefFunctionHarmBalance<T>::UpdateHarm, this ), harmHandle_ );


    // For the callback mechanism
    solHandle_ = mp_->GetNewHandle(true);
    // Bind the handle to the correct expression
    mp_->SetExpr(solHandle_,"updateTrigger");
    // register callback mechanism if expression changes
    mp_->AddExpChangeCallBack( boost::bind(&CoefFunctionHarmBalance<T>::UpdateSolution, this ), solHandle_ );

  }

  template<class T>
  CoefFunctionHarmBalance<T>::
  ~CoefFunctionHarmBalance() {}

  template<class T>
  void CoefFunctionHarmBalance<T>::
  GetScalar(T& coefScal, const LocPointMapped& lpm ){
    EXCEPTION("CoefFunctionHarmBalance::GetScalar NOT IMPLEMENTED YET")
  }

  template<class T>
  void CoefFunctionHarmBalance<T>::
  GetVector(Vector<T>& coefVec, const LocPointMapped& lpm){
    EXCEPTION("CoefFunctionHarmBalance::GetVector NOT IMPLEMENTED YET")
  }

  template<class T>
  void CoefFunctionHarmBalance<T>::
  GetTensor(Matrix<T>& coefMat, const LocPointMapped& lpm){
    EXCEPTION("CoefFunctionHarmBalance::GetTensor NOT IMPLEMENTED YET")
  }

  template<class T>
  UInt CoefFunctionHarmBalance<T>::
  GetVecSize() const{
    EXCEPTION("CoefFunctionHarmBalance::GetVecSize NOT IMPLEMENTED YET")
  }

  template<class T>
  std::string CoefFunctionHarmBalance<T>::
  ToString() const {
    return "CoefFunctionHarmBalance";
  }


  template<class T>
  PtrCoefFct CoefFunctionHarmBalance<T>::
  GenerateMatCoefFnc(const UInt& iRegion,
                     const std::string& name){
    if( name != "Reluctivity"){
      EXCEPTION("CoefFunctionHarmBalance::GenerateMatCoefFnc unrecognized"
                "name of expected quantity:"<<name)
    }

    RegionIdType actRegion;
    BaseMaterial * actMat = NULL;

    actRegion = regions_[iRegion];
    actMat    = materials_[actRegion];

    PtrCoefFct ret = NULL;

    // Check if this is the initial computation, if so, then we act like the material
    // was linear
    if(mp_->Eval(harmHandle_) == maxInt_ ){
      // return the linear reluctivity
      ret = actMat->GetScalCoefFnc(MAG_RELUCTIVITY,Global::REAL );
    }else{
      EXCEPTION("CoefFunctionHarmBalance<T>::GenerateMatCoefFnc Nonlinear part not yet implemented");
      //ret = actMat->GetScalCoefFncNonLin( MAG_RELUCTIVITY, Global::REAL, magFluxCoef_);
    }


    return ret;
  }

  template<class T>
  void CoefFunctionHarmBalance<T>::
  UpdateHarm(){
    std::cout<< this->mp_->GetExpr(harmHandle_  )<<std::endl;
    std::cout<< this->mp_->Eval(harmHandle_) <<std::endl;
  }

  template<class T>
  void CoefFunctionHarmBalance<T>::
  UpdateSolution(){
    std::cout<< this->mp_->GetExpr(solHandle_  )<<std::endl;
    std::cout<< this->mp_->Eval(solHandle_) <<std::endl;

    //here read and cache the solution vector from the PDE



  }


  // Explicit template instantiation
  #ifdef EXPLICIT_TEMPLATE_INSTANTIATION
    template class CoefFunctionHarmBalance<Double>;
    template class CoefFunctionHarmBalance<Complex>;
  #endif

} // end of namespace


