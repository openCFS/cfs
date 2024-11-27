// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#include <limits>

#include "CoefFunctionHarmBalance.hh"
#include <boost/lexical_cast.hpp>
#include "Utils/mathParser/mathParser.hh"
#include "Domain/Mesh/Grid.hh"
#include "Domain/CoordinateSystems/CoordSystem.hh"

#include "Domain/CoefFunction/CoefXpr.hh"

#include "FeBasis/FeSpace.hh"

#include "DataInOut/Logging/LogConfigurator.hh"
#include <Domain/CoefFunction/CoefFunctionMaterialModel.hh>

namespace CoupledField
{
DEFINE_LOG(coeffctharmbalance, "coeffctharmbalance")

  // ==========================================================================
  //  COEFFICIENT FUNCTION HARMONIC BALANCING
  // ==========================================================================
template<class T>
  CoefFunctionHarmBalance<T>::CoefFunctionHarmBalance():CoefFunction(){
  isMH_ = false;
}


  template<class T>
  void CoefFunctionHarmBalance<T>::Init(shared_ptr<BaseFeFunction> feFct,
                                                      shared_ptr<FeSpace> feSpc,
                                                      const StdVector<RegionIdType>& regions,
                                                      std::map<RegionIdType, BaseMaterial*>& materials,
                                                      Grid* ptGrid,
                                                      PtrCoefFct magFluxCoef,
                                                      const UInt& N,
                                                      const UInt& M,
                                                      const Double& baseFreq,
                                                      const UInt& nFFT,
                                                      std::string modelName,
                                                      PtrCoefFct matModelCoef){

    // Initialize the hbRegion_ vector. Therefore loop over every region
    // and assign material...and roughly init the remaining in the struct
    hbRegion_.Resize(regions.GetSize());

    // Inform the coeffunction that this is indeed a multiharmonic analysis
    isMH_ = true;

    for(UInt iRegion = 0; iRegion < regions.GetSize(); ++iRegion){
      HBRegionHelper& st = hbRegion_[iRegion];
      st.elemListPerRegion = NULL;
      RegionIdType actRegion = regions[iRegion];
      st.material = materials[actRegion];
      st.nonLinNuCoefMap = NULL;
      st.linNuCoefMap = NULL;
      st.region = regions[iRegion];
    }

    //! For the initialization from the PDE-class,
    //! Calls the Init method in MHTimeFreqResult
    freqTimeRes_.Init(N, M, baseFreq, nFFT, domain);

    feFct_ = dynamic_pointer_cast<FeFunction<T> >(feFct);

    fluxFESpace_ = dynamic_pointer_cast<FeSpace >(feSpc);

    ptGrid_ = ptGrid;

    magFluxCoef_ = magFluxCoef;

    regionRegistration_ = false;

    maxInt_ = std::numeric_limits<int>::max();

    N_ = N;

    M_ = M;

    updateIter_ = 0;

    numRegions_ = regions.GetSize();

    numRegisteredRegions_ = 0;

    numElems_ = 0;

    nuFreqTmp_.Resize(0);

    // For the callback mechanism
    mp_ = domain->GetMathParser();
    cashHandle_ = mp_->GetNewHandle(true);
    // Bind the handle to the correct expression
    mp_->SetExpr(cashHandle_,"finishCash");
    // Flag for first time calculation
    mp_->SetValue(MathParser::GLOB_HANDLER, "finishCash", std::numeric_limits<unsigned int>::max());
    //mp_->SetValue(MathParser::GLOB_HANDLER, "finishCash", 0);
    // register callback mechanism if expression changes
    mp_->AddExpChangeCallBack( boost::bind(&CoefFunctionHarmBalance<T>::FinishCash, this ), cashHandle_ );


    // For the callback mechanism
    solHandle_ = mp_->GetNewHandle(true);
    // Bind the handle to the correct expression
    mp_->SetExpr(solHandle_,"cacheResult");
    // Flag for first time calculation
    mp_->SetValue(MathParser::GLOB_HANDLER, "cacheResult", -1);
    // register callback mechanism if expression changes
    mp_->AddExpChangeCallBack( boost::bind(&CoefFunctionHarmBalance<T>::CashTimeResult, this ), solHandle_ );

    // For the correct choice of nu's in GetScalar method
    harmonicHandle_ = mp_->GetNewHandle(true);
    // Bind the handle to the correct expression
    mp_->SetExpr(harmonicHandle_,"harmonicHandle");
    // Flag for first time calculation
    mp_->SetValue(MathParser::GLOB_HANDLER, "harmonicHandle", maxInt_);


    // Variables from CoefFunction base-class
    isAnalytic_ = false;
    isComplex_ = true;
    // take care, it's changed throughout the class!!
    dimType_ = CoefFunction::SCALAR;

    modelName_=modelName;
    matModelCoef_= matModelCoef;
  }

  template<class T>
  CoefFunctionHarmBalance<T>::
  ~CoefFunctionHarmBalance() {
    if( isMH_ ){
      mp_->ReleaseHandle(cashHandle_);
      mp_->ReleaseHandle(solHandle_);
      mp_->ReleaseHandle(harmonicHandle_);
    }

  }


  template<class T>
  PtrCoefFct CoefFunctionHarmBalance<T>::
  GenerateMatCoefFnc(const UInt& iRegion,
                     const std::string& name,
                     const bool nonLin,
                     shared_ptr<ElemList> actSDList){
    if( name != "Reluctivity" && name!="Permittivity"){
      EXCEPTION("CoefFunctionHarmBalance::GenerateMatCoefFnc unrecognized"
                "name of expected quantity:"<<name)
    }

    if(name == "Permittivity"){
      RegionIdType actRegion;
          BaseMaterial * actMat = NULL;

          HBRegionHelper & regStruc = hbRegion_[iRegion];

          actRegion = regStruc.region;
          actMat    = regStruc.material;

          // Register the elements
          this->RegisterElemsInRegion(actSDList, iRegion, regStruc);

          PtrCoefFct ret = NULL;


            // Also set the correct nonlinear materialparameter evaluation CoefFunctions here
            // Only for the regions which are really nonlinear!
            if( nonLin && (modelName_ == "nonlinearCurve" )){
              LOG_DBG(coeffctharmbalance) << "Generating nonlinear multiharmonic "
                  "material coefficient function for region"<< actRegion <<" with material "<< actMat;
              regStruc.isNonLin = true;

              // we need the real part of this CoefFunction GetVector because the
              // GetScalCoefFncNonLin method can only handle real values
              // We dont need "Bfield" in electrostatics because we already give the elecFieldCoef
              // shared_ptr<CoefFunction> BField = NULL;
              // BField.reset(new CoefFunctionHarmBalanceEval<Double>(magFluxCoef_, ptGrid_));
              // regStruc.nonLinNuCoefMap = actMat->GetScalCoefFncNonLin( ELEC_PERMITTIVITY_SCALAR, Global::REAL, BField);
              // it seems like i get the same result...

              regStruc.nonLinNuCoefMap = actMat->GetScalCoefFncNonLin( ELEC_PERMITTIVITY_SCALAR, Global::REAL, magFluxCoef_);
              // for the initial solution we also need the linear nu
              regStruc.linNuCoefMap = actMat->GetScalCoefFnc( ELEC_PERMITTIVITY_SCALAR, Global::REAL);

              dimType_ = SCALAR;
              ret = (PtrCoefFct)this;
            } else if(nonLin && (modelName_ != "nonlinearCurve" )){
              LOG_DBG(coeffctharmbalance) << "Generating hysteresis multiharmonic "
                  "material coefficient function for region"<< actRegion <<" with material "<< actMat;
              regStruc.isNonLin = true;

              if(modelName_ == "JilesAthertonModel" ){
                std::map<std::string, double> ParameterMap;

                actMat->GetScalar(ParameterMap["Ps"], ELEC_PS_JILES, Global::REAL );
                actMat->GetScalar(ParameterMap["alpha"], ELEC_ALPHA_JILES, Global::REAL );
                actMat->GetScalar(ParameterMap["a"], ELEC_A_JILES, Global::REAL );
                actMat->GetScalar(ParameterMap["k"], ELEC_K_JILES, Global::REAL );
                actMat->GetScalar(ParameterMap["c"], ELEC_C_JILES, Global::REAL );
                ParameterMap["isMH"] = 1.0;
                matModelCoef_->InitModel( ParameterMap, actSDList->GetSize());
              }

              regStruc.nonLinNuCoefMap = matModelCoef_; //actMat->GetScalCoefFncModel(matModelCoef_);
              // for the initial solution we also need the linear nu
              regStruc.linNuCoefMap = actMat->GetScalCoefFnc( ELEC_PERMITTIVITY_SCALAR, Global::REAL);

              dimType_ = SCALAR;
              ret = (PtrCoefFct)this;

            }else{
              regStruc.isNonLin = false;
              regStruc.nonLinNuCoefMap = actMat->GetScalCoefFnc(ELEC_PERMITTIVITY_SCALAR,Global::REAL );
              //TODO Clean this up, it's the same as nonLinNuCoefMap for the linear material
              regStruc.linNuCoefMap = actMat->GetScalCoefFnc(ELEC_PERMITTIVITY_SCALAR,Global::REAL );
              ret = (PtrCoefFct)this;
            }
            return ret;
    }else{
    RegionIdType actRegion;
    BaseMaterial * actMat = NULL;


    HBRegionHelper & regStruc = hbRegion_[iRegion];


    actRegion = regStruc.region;
    actMat    = regStruc.material;

    // Register the elements
    this->RegisterElemsInRegion(actSDList, iRegion, regStruc);

    PtrCoefFct ret = NULL;


      // Also set the correct nonlinear reluctivity evaluation CoefFunctions here
      // Only for the regions which are really nonlinear!
      if( nonLin ){
        LOG_DBG(coeffctharmbalance) << "Generating nonlinear multiharmonic "
            "material coefficient function for region"<< actRegion <<" with material "<< actMat;
        regStruc.isNonLin = true;

        // we need the real part of this CoefFunction GetVector because the
        // GetScalCoefFncNonLin method can only handle real values
        shared_ptr<CoefFunction> BField = NULL;
        BField.reset(new CoefFunctionHarmBalanceEval<Double>(magFluxCoef_, ptGrid_));

        regStruc.nonLinNuCoefMap = actMat->GetScalCoefFncNonLin( MAG_RELUCTIVITY_SCALAR, Global::REAL, BField);
        // for the initial solution we also need the linear nu
        regStruc.linNuCoefMap = actMat->GetScalCoefFnc( MAG_RELUCTIVITY_SCALAR, Global::REAL);

        dimType_ = SCALAR;
        ret = (PtrCoefFct)this;
      }else{
        regStruc.isNonLin = false;
        regStruc.nonLinNuCoefMap = actMat->GetScalCoefFnc(MAG_RELUCTIVITY_SCALAR,Global::REAL );
        //TODO Clean this up, it's the same as nonLinNuCoefMap for the linear material
        regStruc.linNuCoefMap = actMat->GetScalCoefFnc(MAG_RELUCTIVITY_SCALAR,Global::REAL );
        ret = (PtrCoefFct)this;
      }
    return ret;
    }
  }


  template<class T>
  void CoefFunctionHarmBalance<T>::
  RegisterElemsInRegion(shared_ptr<ElemList> actSDList,
                        const UInt& iRegion,
                        HBRegionHelper& regStruc){
    // Check if the number of allowed regions is exceeded
    if(regionRegistration_ == true){
      EXCEPTION("CoefFunctionHarmBalance::RegisterElemsInRegion Error when trying \n"
          "to set a region, which was probably already registered!");
    }

    LOG_DBG(coeffctharmbalance) << "Registering region "<<regStruc.region<<" with "<< actSDList->GetSize()<<" elements";

    regStruc.elemListPerRegion = actSDList;

    numElems_ += actSDList->GetSize();

    ++numRegisteredRegions_;

    LOG_DBG(coeffctharmbalance) << "Until now, a total of "<< numElems_ <<" elements was registered";

    // Check and set the flag if all regions are set
    if(numRegisteredRegions_ == numRegions_) {
      regionRegistration_ = true;
      LOG_DBG(coeffctharmbalance) << "\t All regions were registered, total number of elements = "<< numElems_;
      LOG_DBG(coeffctharmbalance) << "Setting CoefFuncionHarmBalance-intern solution vector for nu(timestep)";

      // Now set the vector containing nu of elements
      nuFreqTmp_.Resize(numElems_, 0.0);
      UInt elemIterator = 0;
      EntityIterator it;
      for(UInt i = 0; i < hbRegion_.GetSize(); ++i){
        HBRegionHelper& regStruc = hbRegion_[i];
        // obtain the iterator to loop over the elements of the region
        it = regStruc.elemListPerRegion->GetIterator();
        // Loop over every element in that region
        for(it.Begin(); !it.IsEnd(); it++){
          positionOfElem_[it.GetElem()->elemNum] = elemIterator;
          ++elemIterator;
        }
      }

      // Sanity check
      if(elemIterator != numElems_){
        EXCEPTION("CoefFunctionHarmBalance<T>::RegisterElemsInRegion This should really not happen!!");
      }

      // Initialize the time result
      freqTimeRes_.InitTimeResult(elemIterator);
      LOG_DBG(coeffctharmbalance) << "\t Finished setting CoefFuncionHarmBalance-intern solution vector for nu(timestep)";
    }
  }


  template<class T>
  void CoefFunctionHarmBalance<T>::
  FinishCash(){
    LOG_DBG(coeffctharmbalance) << "\t UpdateHarm()";

    // Just for the logging output give the index of solution vector to print
    UInt index = 9;
    if (IS_LOG_ENABLED(coeffctharmbalance, dbg3)) {
      freqTimeRes_.PrintTimeResults(index);
    }

    UInt c = this->mp_->Eval(cashHandle_);
    if( c != std::numeric_limits<unsigned int>::max() ){
      // Perform FFT of time-signal
      freqTimeRes_.TimeToFourier();
    }

    if (IS_LOG_ENABLED(coeffctharmbalance, dbg3)) {
      freqTimeRes_.PrintFreqResults(index);
    }
  }

  template<class T>
  void CoefFunctionHarmBalance<T>::
  CashTimeResult(){
    LOG_DBG(coeffctharmbalance) << "\t Starting UpdateSolution() in CoefFunctionHarmBalance";

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
    // and cache the magnetic reluctivity result


    // These integration parameters can vary between regions
    LocPoint lp;
    EntityIterator it;
    LocPointMapped lpm;
    Double nuOfB = 0.0;
    shared_ptr<ElemShapeMap> esm;


    UInt elemIterator = 0;

    // This is the first call of this method, therefore we
    // store some stuff, we need frequently in the following iterations
    if( this->mp_->Eval(solHandle_) >= 0 ){
          // Loop over every region
          for(UInt i = 0; i < hbRegion_.GetSize(); ++i){
            HBRegionHelper& regStruc = hbRegion_[i];

            // obtain the iterator to loop over the elements of the region
            it = regStruc.elemListPerRegion->GetIterator();

            // Loop over every element in that region
            for(it.Begin(); !it.IsEnd(); it++){
              const Elem * el = it.GetElem();
              // esm = it.GetGrid()->GetElemShapeMap(el, true);
              // esm.reset(el->ptrShapeMap);
              // esm = el->ptrShapeMap;
              shared_ptr<ElemShapeMap> esm = (el)->GetElemShapeMap(it.GetGrid(), true);

              // Where do we evaluate the magnetic flux density?
              // Element or integration point -> see UPDATE from above
              lp = Elem::shapes[el->type].midPointCoord;
              lpm.Set(lp, esm, 0.0);

              // Evaluate the nu(B)
              regStruc.nonLinNuCoefMap->GetScalar(nuOfB, lpm);

              // now cache it for the element
              nuFreqTmp_[ elemIterator ] = nuOfB;
//              std::cout << "----CashResult----" << std::endl;
//              std::cout << "ElemIter = " << elemIterator << std::endl;
//              std::cout << "Epsilon = " << nuOfB << std::endl;
              ++elemIterator;
            }
            LOG_DBG(coeffctharmbalance) << "nu of virtual timestep "
                <<this->mp_->Eval(solHandle_)<<" : "<<nuFreqTmp_.ToString();
          } // loop over every region


      // Now fill the appropriate time result in MHTimeFreqResult
      freqTimeRes_.SetTimeResult(nuFreqTmp_, this->mp_->Eval(solHandle_) );

    }


  }

  template<>
  void CoefFunctionHarmBalance<Double>::
  GetScalar(Double& coefScal, const LocPointMapped& lpm ){
    EXCEPTION("CoefFunctionHarmBalance<T>::GetScalar called in Double version, this is not allowed!")
  }

  template<>
  void CoefFunctionHarmBalance<Complex>::
  GetScalar(Complex& coefScal, const LocPointMapped& lpm ){
    Double coefScalReal = 0.0;

    RegionIdType elemReg = lpm.ptEl->regionId;
    int harmonic = this->mp_->Eval(harmonicHandle_);

    if( harmonic == maxInt_ ){
      // For the initial multiharmonic iteration
      // loop over regions and get the region of the lpm
      for(auto reg : hbRegion_){
        if ( lpm.isSurface ){
          const MortarNcSurfElem* e = dynamic_cast<const MortarNcSurfElem*>(lpm.ptEl);
          if( reg.region == e->ptSecondary->ptVolElems[0]->regionId ){
            reg.linNuCoefMap->GetScalar(coefScalReal, lpm);
            coefScal = (Complex)coefScalReal;
            break;
          }
        }else{
          if( reg.region == elemReg){
            reg.linNuCoefMap->GetScalar(coefScalReal, lpm);
            coefScal = (Complex)coefScalReal;
            break;
          }
        }
      }
    }else{
      const Vector<Complex>& fR = freqTimeRes_.GetFreqResult( N_ + harmonic);
      if ( lpm.isSurface ){
        const MortarNcSurfElem* e = dynamic_cast<const MortarNcSurfElem*>(lpm.ptEl);
        coefScal = fR[ positionOfElem_[e->ptSecondary->ptVolElems[0]->elemNum] ];
        elemReg = e->ptSecondary->ptVolElems[0]->regionId;
        LOG_DBG(coeffctharmbalance) <<"nu for Nitsche interface with volume region " << ptGrid_->GetRegionName( elemReg )<<
                                    " in harmonic "<< harmonic <<" = " <<coefScal;
      }else{
        coefScal = fR[ positionOfElem_[lpm.ptEl->elemNum] ];
        LOG_DBG(coeffctharmbalance) <<"nu for region " << ptGrid_->GetRegionName(elemReg)<<" in harmonic "<< harmonic <<" = " <<coefScal;
      }
    }


  }




  template<class T>
  void CoefFunctionHarmBalance<T>::
  GetVector(Vector<T>& coefVec, const LocPointMapped& lpm){
    EXCEPTION("CoefFunctionHarmBalance::GetVector NOT IMPLEMENTED HERE");
  }


  template<class T>
  void CoefFunctionHarmBalance<T>::
  GetTensor(Matrix<T>& coefMat, const LocPointMapped& lpm){
    EXCEPTION("CoefFunctionHarmBalance::GetTensor NOT IMPLEMENTED YET");
  }

  template<class T>
  UInt CoefFunctionHarmBalance<T>::
  GetVecSize() const{
    // spatial dimension (dimension of B vector)
    return ptGrid_->GetDim();
  }

  // ==========================================================================
  //  COEFFICIENT FUNCTION HARMONIC BALANCING EVALUATE B FIELD
  // ==========================================================================
  template<class T>
  CoefFunctionHarmBalanceEval<T>::CoefFunctionHarmBalanceEval(PtrCoefFct magFluxCoef,
                                                              Grid* ptGrid):CoefFunction() {

    mp_ = domain->GetMathParser();

    ptGrid_ = ptGrid;


    //PtrCoefFct realBField = CoefFunction::Generate( mp_, Global::REAL, CoefXprUnaryOp( mp_, magFluxCoef, CoefXpr::OP_RE));

    magFluxCoef_ = magFluxCoef;

    // Variables from CoefFunction base-class
    isAnalytic_ = false;
    isComplex_ = false;
    // take care, it's changed throughout the class!!
    dimType_ = CoefFunction::VECTOR;
  }

  template<class T>
  CoefFunctionHarmBalanceEval<T>::
  ~CoefFunctionHarmBalanceEval() {}

  template<class T>
  void CoefFunctionHarmBalanceEval<T>::
  GetVector(Vector<T>& coefVec, const LocPointMapped& lpm){
      // This is just for testing, don't take anything for granted from here !!!!!!!!!
      Vector<Complex> tmp;
      magFluxCoef_->GetVector(tmp, lpm);
      coefVec = tmp.GetPart(Global::REAL);
  }

  template<class T>
  void CoefFunctionHarmBalanceEval<T>::
  GetTensor(Matrix<T>& coefMat, const LocPointMapped& lpm){
    EXCEPTION("CoefFunctionHarmBalance::GetTensor NOT IMPLEMENTED YET")
  }

  template<class T>
  void CoefFunctionHarmBalanceEval<T>::
  GetScalar(T& coefMat, const LocPointMapped& lpm){
    EXCEPTION("CoefFunctionHarmBalance::GetScalar NOT IMPLEMENTED YET")
  }

  template<class T>
  UInt CoefFunctionHarmBalanceEval<T>::
  GetVecSize() const{
    // spatial dimension (dimension of B vector)
    return ptGrid_->GetDim();
  }


  // Explicit template instantiation
  template class CoefFunctionHarmBalance<Double>;
  template class CoefFunctionHarmBalance<Complex>;
  template class CoefFunctionHarmBalanceEval<Double>;
  template class CoefFunctionHarmBalanceEval<Complex>;


} // end of namespace
