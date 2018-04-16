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

    regionRegistration_ = false;

    maxInt_ = std::numeric_limits<unsigned int>::max();

    updateIter_ = 0;

    // Reset the frequency results
    freqResult_.Resize(0,0);

    numRegions_ = regions.GetSize();

    elemListPerRegion_.Resize(0);

    numElems_ = 0;


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
    mp_->SetExpr(solHandle_,"cacheResult");
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
  RegisterElemsInRegion(shared_ptr<ElemList> actSDList){
    // Check if the number of allowed regions is exceeded
    if(regionRegistration_ == true){
      EXCEPTION("CoefFunctionHarmBalance::RegisterElemsInRegion Error when trying \n"
          "to set a region, which was probably already registered!");
    }

    elemListPerRegion_.Push_back(actSDList);
    numElems_ += actSDList->GetSize();

    // Check and set the flag if all regions are set
    if(elemListPerRegion_.GetSize() == numRegions_) {
      regionRegistration_ = true;
    }
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

    /*
     * 1. Check if the update callback was called more often than 2*N_+1
     *    if so, throw an exception
     *    2.true We can perform the inverse FFT of the nu(t) signal
     *           This should be somehow done in the MHTimeFreqResult class
     *           But how can we give the solution to this class?
     *           Maybe create an instance right here in this class
     *           and us it as the nu-object in accordance to the python tool
     *    2.false Cache the solution vectors until 1. is true...I think that
     *            is more performant because afterwards we can loop over all
     *            vectors in parallel and evaluate B where we need
     *            it IN PARALLEL.
     *            But therefore we need to really copy the vectors...
     *            we'll see how the memory is affected
     */


    //=================================================================================
    //=============================== Just for testing purposes =======================
    if( this->mp_->Eval(solHandle_)== 0 ){
      // testwise dummy evaluation of magnetic flux density at some element point

      EntityIterator it = elemListPerRegion_[0]->GetIterator();
      for(  it.Begin(); !it.IsEnd(); it++  ){
        const Elem * el = it.GetElem();

        LocPoint lp = Elem::shapes[el->type].midPointCoord;
        std::cout<<lp.coord.ToString()<<std::endl;
        LocPointMapped lpm;
        shared_ptr<ElemShapeMap> esm = it.GetGrid()->GetElemShapeMap(el, true);
        lpm.Set(lp, esm, 0.0);

        Vector<Complex> dummyFluxVec;
        magFluxCoef_->GetVector(dummyFluxVec, lpm);
        std::cout<<dummyFluxVec.ToString()<<std::endl;
      }
      //==============================================================================



    }else{
      EXCEPTION("This case is not yet handled");
    }


  }


  // Explicit template instantiation
  #ifdef EXPLICIT_TEMPLATE_INSTANTIATION
    template class CoefFunctionHarmBalance<Double>;
    template class CoefFunctionHarmBalance<Complex>;
  #endif

} // end of namespace


