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

#include "FeBasis/FeSpace.hh"

namespace CoupledField
{
  // ==========================================================================
  //  COEFFICIENT FUNCTION HARMONIC BALANCING
  // ==========================================================================
  template<class T>
  CoefFunctionHarmBalance<T>::CoefFunctionHarmBalance(shared_ptr<BaseFeFunction> feFct,
                                                      shared_ptr<FeSpace> feSpc,
                                                      const StdVector<RegionIdType>& regions,
                                                      const std::map<RegionIdType, BaseMaterial*>& materials,
                                                       Grid* ptGrid,
                                                      PtrCoefFct magFluxCoef):CoefFunction() {

    feFct_ = dynamic_pointer_cast<FeFunction<T> >(feFct);

    fluxFESpace_ = dynamic_pointer_cast<FeSpace >(feSpc);

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
    regionList_.Resize(0);

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
    // Flag for first time calculation
    mp_->SetValue(MathParser::GLOB_HANDLER, "cacheResult", -1);
    // register callback mechanism if expression changes
    mp_->AddExpChangeCallBack( boost::bind(&CoefFunctionHarmBalance<T>::UpdateSolution, this ), solHandle_ );



    isAnalytic_ = false;
    isComplex_ = false;

    // take care, it's changed throughout the class!!
    dimType_ = CoefFunction::SCALAR;

  }

  template<class T>
  CoefFunctionHarmBalance<T>::
  ~CoefFunctionHarmBalance() {}


  template<class T>
  PtrCoefFct CoefFunctionHarmBalance<T>::
  GenerateMatCoefFnc(const UInt& iRegion,
                     const std::string& name,
                     const bool nonLin){
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

      // Also set the correct nonlinear reluctivity evaluation CoefFunctions here
      // ONLY FOR THE REGIONS WHICH ARE REALLY NONLINEAR!!!
      if( nonLin ){
        // we need the real part of this CoefFunction GetVector because the
        // GetScalCoefFncNonLin method can only handle real values

        dimType_ = VECTOR;
        PtrCoefFct realThis = CoefFunction::Generate( mp_, Global::REAL, CoefXprUnaryOp( mp_, (PtrCoefFct)this, CoefXpr::OP_RE));
        //PtrCoefFct realThis = CoefFunction::Generate( mp_, Global::REAL, CoefXprUnaryOp( mp_, magFluxCoef_, CoefXpr::OP_RE));
        JUSTATEST_ = realThis;
        nonLinNuCoefMap_[actRegion] = actMat->GetScalCoefFncNonLin( MAG_RELUCTIVITY, Global::REAL, realThis);

      }else{
        nonLinNuCoefMap_[actRegion] = ret;
      }

    }else{
      EXCEPTION("CoefFunctionHarmBalance<T>::GenerateMatCoefFnc Nonlinear part not yet implemented");
    }
    return ret;
  }


  template<class T>
  void CoefFunctionHarmBalance<T>::
  RegisterElemsInRegion(shared_ptr<ElemList> actSDList,
                        const UInt& iRegion){
    // Check if the number of allowed regions is exceeded
    if(regionRegistration_ == true){
      EXCEPTION("CoefFunctionHarmBalance::RegisterElemsInRegion Error when trying \n"
          "to set a region, which was probably already registered!");
    }

    elemListPerRegion_.Push_back(actSDList);
    regionList_.Push_back(iRegion);
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
     *
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
     *
     * UPDATE:  It does not mind if we evaluate the magnetic flux density
     *          only on the element midpoint because if we calculate it
     *          on integrations points, we get the same value for every
     *          integration point of one element, since the specific
     *          coef function, which is necessary to compute the B field
     *          calls FeFunctions::GetElemSolution ... which is one B-value
     *          for one element
     */

    // Loop over all regions and over every integration point in this region
    // and cache the magnetic flux density result





    if( this->mp_->Eval(solHandle_) >= 0 ){

      // These integration parameters can vary between regions
      IntegOrder order;
      IntScheme::IntegMethod method;
      StdVector<LocPoint> intPoints;
      StdVector<Double> weights;
      LocPointMapped lpm;
      Vector<Complex> dummyFluxVec;

      // Loop over every region
      for(UInt i = 0; i < elemListPerRegion_.GetSize(); ++i){

        EntityIterator it = elemListPerRegion_[i]->GetIterator();
        RegionIdType regId = regions_[ regionList_[i] ];
        fluxFESpace_->GetFe(it, method, order);
        // Loop over every element in that region
        for(it.Begin(); !it.IsEnd(); it++){
          const Elem * el = it.GetElem();
          shared_ptr<ElemShapeMap> esm = it.GetGrid()->GetElemShapeMap(el, true);
          // Where do we evaluate the magnetic flux density?
          // Element or integration point -> see UPDATE from above
          LocPoint lp = Elem::shapes[el->type].midPointCoord;
          lpm.Set(lp, esm, 0.0);
          magFluxCoef_->GetVector(dummyFluxVec, lpm);
//          Evaluate the nu right here!!!!!!!!


std::cout<<"dummyFluxVec.ToString()"<<dummyFluxVec.ToString()<<std::endl;

          Vector<Double> test;
          JUSTATEST_->GetVector(test, lpm);
std::cout<<"test.ToString()"<<test.ToString()<<std::endl;


          Double dummy;
          nonLinNuCoefMap_[regId]->GetScalar(dummy, lpm);
std::cout<<"dummy.ToString()"<<dummy<<std::endl;


        }
      }

EXCEPTION("WE STILL NEED TO CACHE IT!");

    }
  }



  template<class T>
  void CoefFunctionHarmBalance<T>::
  GetScalar(T& coefScal, const LocPointMapped& lpm ){
    EXCEPTION("CoefFunctionHarmBalance::GetScalar NOT IMPLEMENTED YET")
  }

  template<class T>
  void CoefFunctionHarmBalance<T>::
  GetVector(Vector<T>& coefVec, const LocPointMapped& lpm){

    // This is just for testing, don't take anything for granted from here !!!!!!!!!
    Vector<Complex> tmp;
    magFluxCoef_->GetVector(tmp, lpm);
    coefVec = tmp.GetPart(Global::REAL);
    //std::cout<<"=========================================================\n"
    //           "========================================================="<<std::endl;
    //EXCEPTION("CoefFunctionHarmBalance::GetVector NOT IMPLEMENTED YET")
  }

  template<class T>
  void CoefFunctionHarmBalance<T>::
  GetTensor(Matrix<T>& coefMat, const LocPointMapped& lpm){
    EXCEPTION("CoefFunctionHarmBalance::GetTensor NOT IMPLEMENTED YET")
  }

  template<class T>
  UInt CoefFunctionHarmBalance<T>::
  GetVecSize() const{
    // spatial dimension (dimension of B vector)
    return ptGrid_->GetDim();
  }

  template<class T>
  std::string CoefFunctionHarmBalance<T>::
  ToString() const {
    return "CoefFunctionHarmBalance";
  }



  // Explicit template instantiation
  #ifdef EXPLICIT_TEMPLATE_INSTANTIATION
    template class CoefFunctionHarmBalance<Double>;
    template class CoefFunctionHarmBalance<Complex>;
  #endif

} // end of namespace


