#include <fstream>

#include "MagEdgePDE.hh"
#include "MagEdgeSpecialAVPDE.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Utils/Coil.hh"
#include "Utils/SmoothSpline.hh"
#include "Utils/LinInterpolate.hh"

#include "Driver/Assemble.hh"
#include "FeBasis/HCurl/FeSpaceHCurlHi.hh"
#include "FeBasis/FeSpaceConst.hh"

#include "FeBasis/HCurl/HCurlElems.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "Domain/CoefFunction/CoefFunctionExpression.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefXpr.hh"



// forms
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/BiLinWrappedLinForm.hh"
#include "Forms/BiLinForms/SingleEntryBiLinInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Forms/Operators/IdentityOperator.hh"

//time stepping
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"

#include "Driver/MultiHarmonicDriver.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
namespace CoupledField {

// declare class specific logging stream
DEFINE_LOG(magEdgeSpecialAVPde, "magEdgeSpecialAVPde")


  // **************
  //  Constructor
  // **************
  MagEdgeSpecialAVPDE::MagEdgeSpecialAVPDE( Grid * aptgrid, PtrParamNode paramNode,
                          PtrParamNode infoNode,
                          shared_ptr<SimState> simState, Domain* domain )
    :MagEdgePDE( aptgrid, paramNode, infoNode, simState, domain ) {

    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "magneticEdgeSpecialAV";
    pdematerialclass_ = ELECTROMAGNETIC;

    if(analysistype_ == BasePDE::AnalysisType::STATIC){
      EXCEPTION("You are using the special A-V formulation for a static analysis...do you really know what you're doing?")
    }
//    if(analysistype_ == BasePDE::AnalysisType::TRANSIENT){
//      EXCEPTION("Special A-V formulation is currently not implemented...but it's no big deal ;)")
//    }

    useModifiedAVVoltageFormulation_ = false;
    useModifiedAVCurrentFormulation_ = false;
    PtrParamNode coilNode = myParam_->Get( "coilList", ParamNode::PASS );
    if ( coilNode ){
      ParamNodeList coilNodes = coilNode->GetChildren();
      for( UInt k = 0; k < coilNodes.GetSize(); k++ ){
        if( coilNodes[k]->Has("source") ){
          std::string exType = coilNodes[k]->Get("source")->Get("type")->As<std::string>();
          if( exType == "voltage" ){
            EXCEPTION("MagEdgeSpecialAVPDE can only have specialvoltage or specialcurrent excitation!");
          }else if( exType == "specialvoltage"){
            useModifiedAVVoltageFormulation_ = true;
          }else if ( exType == "specialcurrent"){
            useModifiedAVCurrentFormulation_ = true;
          }
        }
      }
    }

  }


  // *************
  //  Destructor
  // *************
  MagEdgeSpecialAVPDE::~MagEdgeSpecialAVPDE() {
  }


  void MagEdgeSpecialAVPDE::ReadSpecialBCs() {

    this->ReadCoils();

  }


  // *****************************
  //  Definition of Integrators
  // *****************************
  void MagEdgeSpecialAVPDE::DefineIntegrators() {

    DefineStandardIntegrators();

    this->DefineCoilIntegrators();
  }// end DefineIntegrators


  void MagEdgeSpecialAVPDE::DefineCoilIntegrators(){
    // ============================
    // COIL INTEGRATORS
    // ============================
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;

    shared_ptr<BaseFeFunction> feFunc = feFunctions_[MAG_POTENTIAL];
    shared_ptr<FeSpace> feSpace = feFunc->GetFeSpace();



    std::map<Coil::IdType, shared_ptr<Coil> >::iterator coilIt;
    coilIt = coils_.begin();
    for( ; coilIt != coils_.end(); coilIt++ ) {
      Coil& actCoil = *(coilIt->second);
      // run over all parts
      std::map<RegionIdType,shared_ptr<Coil::Part> >::iterator partIt;
      partIt = actCoil.parts_.begin();


      if( actCoil.sourceType_ == Coil::SPECIALVOLTAGE ||
          actCoil.sourceType_ == Coil::SPECIALVOLTAGE_MULTHARM) {
        /*  ==================================================
         *            SPECIALVOLTAGE FORMULATION
            ================================================== */
        for( partIt = actCoil.parts_.begin();
            partIt != actCoil.parts_.end();
            partIt++ ) {
          Coil::Part & actPart = *(partIt->second);
          RegionIdType actRegion = partIt->first;
          shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
          actSDList->SetRegion( actRegion );


          /*
            for the non-multiharmonic case it is a simple PtrCoefFct type sufficient
            but due to the flexibility to use multiharmonic excitation,
            we need a map of PtrCoefFct's, where the key is the harmonic number.
            So in the non-multiharmonic case, we simply have one key "0"
          */
          std::map<UInt, PtrCoefFct> jFct;
          if( actCoil.sourceType_ == Coil::SPECIALVOLTAGE ){
            CoefXprVecScalOp iVec = CoefXprVecScalOp(mp_, actPart.jUnitVec, actCoil.srcVal_, CoefXpr::OP_MULT);
            PtrCoefFct iFct = CoefFunction::Generate(mp_, part, iVec);
            CoefXprVecScalOp jVec = CoefXprVecScalOp(mp_, iFct, boost::lexical_cast<std::string>(actPart.wireCrossSect), CoefXpr::OP_DIV);
            jFct[0] = CoefFunction::Generate(mp_, part, jVec);
          } else if( actCoil.sourceType_ == Coil::SPECIALVOLTAGE_MULTHARM ){
            if( (actCoil.srcValMH_.find(0) != actCoil.srcValMH_.end()) && (!dynamic_cast<MultiHarmonicDriver*>(domain_->GetSingleDriver())->fullSystem_) ){
              EXCEPTION("You specified an excitation in harmonic 0 but didn't set the"
                        "<fullSystem> tag to true in the analysis section!");
            }

            // loop over all inserted harmonics of the xml file; type is std::map<Integer, PtrCoefFct>
            for(auto& h : actCoil.srcValMH_){
              CoefXprVecScalOp iVec = CoefXprVecScalOp(mp_, actPart.jUnitVec, h.second  , CoefXpr::OP_MULT);
              PtrCoefFct iFct = CoefFunction::Generate(mp_, part, iVec);
              CoefXprVecScalOp jVec = CoefXprVecScalOp(mp_, iFct, boost::lexical_cast<std::string>(actPart.wireCrossSect), CoefXpr::OP_DIV);
              jFct[h.first] = CoefFunction::Generate(mp_, part, jVec);
            }
          }

          Double conductivity = 0.0;
          materials_[actRegion]->GetScalar(conductivity,MAG_CONDUCTIVITY_SCALAR,Global::REAL);
          PtrCoefFct conduccoef = CoefFunction::Generate(mp_, Global::REAL, lexical_cast<std::string>(conductivity));

          // This switch is necessary because we need to select which harmonic component
          // we want to assemble in the rhs
          if( actCoil.sourceType_ == Coil::SPECIALVOLTAGE_MULTHARM ){
            // Multiharmonic Case
            LinearForm* curInt;
            for(auto& h : jFct){
              CoefXprVecScalOp uVec = CoefXprVecScalOp(mp_, h.second, conduccoef, CoefXpr::OP_MULT);
              PtrCoefFct sigma_gradV = CoefFunction::Generate(mp_, part, uVec);
              coilCurrentDens_[actRegion] = h.second;
              curInt = new BUIntegrator<Complex>( new IdentityOperator<FeHCurl,3,1,Complex>(),
                                              -1.0, sigma_gradV, updatedGeo_);

              curInt->SetName("CoilIntegrator");
              curInt->SetHarm(h.first);
              LinearFormContext * coilContext = new LinearFormContext( curInt );
              coilContext->SetEntities( actSDList );
              coilContext->SetFeFunction( feFunc );
              assemble_->AddLinearForm( coilContext );
            }
          }else{
            // Classic Case
            CoefXprVecScalOp uVec = CoefXprVecScalOp(mp_, jFct[0], conduccoef, CoefXpr::OP_MULT);
            PtrCoefFct sigma_gradV = CoefFunction::Generate(mp_, part, uVec);
            coilCurrentDens_[actRegion] = jFct[0];

            LinearForm* curInt;
            if( isComplex_ ){
               curInt = new BUIntegrator<Complex>( new IdentityOperator<FeHCurl,3,1,Complex>(),
                   -1.0, sigma_gradV, updatedGeo_);
            }else{
              curInt = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,3,1,Double>(),
                  -1.0, sigma_gradV, updatedGeo_);
            }

            curInt->SetName("CoilIntegrator");
            LinearFormContext * coilContext = new LinearFormContext( curInt );
            coilContext->SetEntities( actSDList );
            coilContext->SetFeFunction( feFunc );
            assemble_->AddLinearForm( coilContext );
          }
        } // loop: parts

      }else if( actCoil.sourceType_ == Coil::SPECIALCURRENT ){
        /*  ==================================================
         *            SPECIALCURRENT FORMULATION
            ================================================== */

        shared_ptr<CoilList> singleCoilList( new CoilList( ptGrid_ ) );
        singleCoilList->AddCoil( coilIt->second );
        feFunctions_[COIL_VOLTAGE_INTEGRAL]->AddEntityList( singleCoilList );

        for( partIt = actCoil.parts_.begin();
             partIt != actCoil.parts_.end();
             partIt++ ) {

          Coil::Part & actPart = *(partIt->second);
          CoefXprVecScalOp eJscaledOp = CoefXprVecScalOp( mp_, actPart.jUnitVec, boost::lexical_cast<std::string>(actPart.wireCrossSect), CoefXpr::OP_DIV );
          PtrCoefFct eJscaled = CoefFunction::Generate( mp_, part, eJscaledOp );

          shared_ptr<ElemList> actSDList( new ElemList( ptGrid_ ) );
          RegionIdType actRegion = partIt->first;
          actSDList->SetRegion( actRegion );


          shared_ptr<FeSpace> fsV = feFunctions_[COIL_VOLTAGE_INTEGRAL]->GetFeSpace();
          fsV->InsertElemsToCoilList(actSDList, singleCoilList);

          coilCurrentDens_[actRegion] = eJscaled;


          Double conductivity = 0.0;
          materials_[actRegion]->GetScalar(conductivity,MAG_CONDUCTIVITY_SCALAR,Global::REAL);
          PtrCoefFct conduccoef = CoefFunction::Generate(mp_, Global::REAL, lexical_cast<std::string>(conductivity));
          CoefXprVecScalOp uVec = CoefXprVecScalOp(mp_, eJscaled, conduccoef, CoefXpr::OP_MULT);
          PtrCoefFct sigma_gradV = CoefFunction::Generate(mp_, part, uVec);

          // === UPPER RIGHT PART (coupling) ===
          LinearForm* upperInt;
          if( isComplex_ ) {
            upperInt = new BUIntegrator<Complex>( new IdentityOperator<FeHCurl,3,1,Complex>(),
                1.0, sigma_gradV, updatedGeo_);
          } else {
            upperInt = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,3,1,Double>(),
                1.0, sigma_gradV, updatedGeo_);
          }
          upperInt->SetName("CoilVoltCouplInt");

          bool assembleTransposed = false;
          BiLinearForm* pseudoBiLin = new BiLinWrappedLinForm( upperInt, assembleTransposed );
          BiLinFormContext* currCoilContext = new BiLinFormContext( pseudoBiLin, DAMPING );
          currCoilContext->SetEntities( actSDList, singleCoilList );
          currCoilContext->SetFeFunctions( feFunc, feFunctions_[COIL_VOLTAGE_INTEGRAL] );
          currCoilContext->SetCounterPart(true);
          assemble_->AddBiLinearForm( currCoilContext );

          // === Lower Diagonal PART (single line per coil) ===
          // compute coefficient: gamma*Gard(V).Grad(V') where Gard(V).Grad(V') is read from the first step
          PtrCoefFct totR = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_, lexical_cast<std::string>(gradVsource_[partIt->second]), conduccoef, CoefXpr::OP_MULT) );
          SingleEntryBiLinInt * totRBiLin = new SingleEntryBiLinInt(1, totR, mp_);
          totRBiLin->SetName( "LowerDiagIntegrator" );
          BiLinFormContext* totRcontext = new BiLinFormContext( totRBiLin, DAMPING );
          totRcontext->SetEntities( singleCoilList, singleCoilList );
          totRcontext->SetFeFunctions( feFunctions_[COIL_VOLTAGE_INTEGRAL], feFunctions_[COIL_VOLTAGE_INTEGRAL] );
          assemble_->AddBiLinearForm( totRcontext );

        } // loop: parts

        // set I on the rhs
        LinearForm* currInt = new SingleEntryInt( actCoil.srcVal_ );
        currInt->SetName( "CoilCurrentLoadInt" );
        LinearFormContext* currContext = new LinearFormContext( currInt );
        currContext->SetEntities( singleCoilList );
        currContext->SetFeFunction( feFunctions_[COIL_VOLTAGE_INTEGRAL] );
        assemble_->AddLinearForm( currContext );
      } // if: current / voltage driven
    } // loop: coils
  } // end DefineCoilIntegrators

  


  // ======================================================
  // TIME-STEPPING SECTION
  // ======================================================

  void MagEdgeSpecialAVPDE::InitTimeStepping() {
	// Use complete implicit scheme
    Double gamma = 1.0;
    GLMScheme * scheme = new Trapezoidal(gamma);
    TimeSchemeGLM::NonLinType nlType = (nonLin_)? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType) );

    feFunctions_[MAG_POTENTIAL]->SetTimeScheme(myScheme);

    if( useModifiedAVCurrentFormulation_ ){
      // Important: Create a new time scheme just for the current unknowns, as otherwise the
      // size of the vectors does not match!
      GLMScheme * scheme2 = new Trapezoidal(gamma);
      shared_ptr<BaseTimeScheme> myScheme2(new TimeSchemeGLM(scheme2, 0, nlType) );

      feFunctions_[COIL_VOLTAGE_INTEGRAL]->SetTimeScheme(myScheme2);
    }

  }

  // ******************************************************
  //   Query parameter object for information about coils
  // ******************************************************
  void MagEdgeSpecialAVPDE::ReadCoils() {
    // Check if the element "coils" is present at all.
    // Otherwise leave
    PtrParamNode coilNode = myParam_->Get( "coilList", ParamNode::PASS );
    PtrParamNode coilInfoNode = myInfo_->Get( "coilList", ParamNode::PASS );
    if ( !coilNode )
      return;

    // Get single coil nodes
    ParamNodeList coilNodes = coilNode->GetChildren();

    // Trigger reading in of definitions
    Global::ComplexPart cplx = isComplex_ ? Global::COMPLEX : Global::REAL;
    if( coilNodes.GetSize() > 0 ) {
      for( UInt i = 0; i < coilNodes.GetSize(); i++ ) {

        // get coil and id
        std::string coilId = coilNodes[i]->Get("id")->As<std::string>();

        // Check if coil with same ID already exists
        if( coils_.find(coilId) != coils_.end() ) {
          EXCEPTION("A coil with ID '" << coilId << "' was already defined.")
        }

        // Create new coil
        shared_ptr<Coil> actCoil( new Coil( coilNodes[i], coilInfoNode,
            ptGrid_, mp_, cplx ) );
        coils_[coilId] = actCoil;

        // Associate mapping of coil parts with regions
        std::map<RegionIdType, shared_ptr<Coil::Part> >::const_iterator it;
        for( it = actCoil->parts_.begin(); it != actCoil->parts_.end(); it++ ) {
          coilRegions_[it->first] = actCoil;
        }
      }

      // Insert the current densities which are defined externally (simulation or sequence step).
      // This is done here because it is impossible for the coil to use a PDE pointer.
      // We have to distinguish between external current density direction and external source.
      // External source includes the direction, but not vice versa. Therefore, the external source
      // must be stored per part anyway, although it counts for the whole coil. Additionally, the
      // parts need the regions and coef functions.
      std::map<Coil::IdType, shared_ptr<Coil> >::iterator coilIt;
      for( coilIt = coils_.begin(); coilIt != coils_.end(); ++coilIt ){
        std::map<shared_ptr<Coil::Part>, PtrParamNode >::iterator extPartIt;
        for( extPartIt = coilIt->second->partsExtJDir_.begin();
            extPartIt != coilIt->second->partsExtJDir_.end(); ++extPartIt ){
          PtrParamNode extNode = extPartIt->second;
          // determine if normalise is set
          bool normalise = true;
          if ( extNode->Has("normalise") ) {
            if ( extNode->Get("normalise")->As<std::string>() == "no" ) {
              normalise = false;
            }
          }
          if( useModifiedAVVoltageFormulation_ || useModifiedAVCurrentFormulation_){
            if(extPartIt->first->wireCrossSect != 1.0){
              EXCEPTION("In the SPECIALVOLTAGE or SPECIALCURRENT version, the wireCrossSect of the coil needs to be 1.0!!");
            }
            if(!extNode->Has("normalise") || extNode->Get("normalise")->As<std::string>() == "yes"){
              EXCEPTION("In the SPECIALVOLTAGE or SPECIALCURRENT version, the normalization of current needs to be switched OFF!!");
            }
          }


          // If we use the special current A-V formulation, we also need the history
          // postprocessing result of the elecCurrentPDE
          std::map<shared_ptr<Coil::Part>, PtrParamNode >::iterator extgVIt;
          for( extgVIt = coilIt->second->partsExtIntgVgV_.begin();
              extgVIt != coilIt->second->partsExtIntgVgV_.end(); ++extgVIt ){
            PtrParamNode extNode = extgVIt->second;

            shared_ptr<CoefFunctionMulti> currDens(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, isComplex_));
            for( UInt k_reg = 0; k_reg < extgVIt->first->regions.GetSize(); ++k_reg ){
              Vector<Double> histVal;
              std::string regionName = ptGrid_->GetRegion().ToString(extgVIt->first->regions[k_reg]);
              ReadUserHistValues(extNode, ResultInfo::SCALAR, histVal, regionName);
              gradVsource_[extgVIt->first] = histVal[0];
            }
          }

          shared_ptr<CoefFunctionMulti> unitCurrDens(new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1,
              isComplex_));
          shared_ptr<CoefFunctionMulti> currDens(new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1,
              isComplex_));
          for( UInt k_reg = 0; k_reg < extPartIt->first->regions.GetSize(); ++k_reg ){
            std::string regName = ptGrid_->regionData[extPartIt->first->regions[k_reg]].name;
            shared_ptr<EntityList> elems;
            elems = ptGrid_->GetEntityList( EntityList::ELEM_LIST, regName );
            PtrCoefFct regCurrDens; // ReadUserFieldValues assigns a value to this
            StdVector<std::string> vecComponents;
            vecComponents = "x", "y", "z";
            std::set<UInt> definedDofs; // ReadUserFieldValues assigns a value to this
            bool updateGeo; // ReadUserFieldValues assigns a value to this
            ReadUserFieldValues(elems,extNode,vecComponents,
                ResultInfo::VECTOR,isComplex_,regCurrDens,
                definedDofs,updateGeo);
            // take the read values and normalise to a length of 1
            PtrCoefFct unitDir;
            if ( normalise ) {
              cplx = Global::REAL;
              CoefXprUnaryOp dirAbsOp = CoefXprUnaryOp( mp_, regCurrDens, CoefXpr::OP_NORM );
              PtrCoefFct dirAbs = CoefFunction::Generate( mp_, cplx, dirAbsOp );
              CoefXprVecScalOp unitOp = CoefXprVecScalOp( mp_, regCurrDens, dirAbs, CoefXpr::OP_DIV );
              unitDir = CoefFunction::Generate( mp_, cplx, unitOp );
            }
            else {
              unitDir = regCurrDens;
            }
            unitCurrDens->AddRegion(extPartIt->first->regions[k_reg],unitDir);
            if( coilIt->second->sourceType_ == Coil::EXTERNAL ){
              currDens->AddRegion(extPartIt->first->regions[k_reg],regCurrDens);
            }
          }
          extPartIt->first->jUnitVec = unitCurrDens;
          if( coilIt->second->sourceType_ == Coil::EXTERNAL ){
            coilPartsExtJ_[extPartIt->first] = currDens;
          }
        }

      }
    }
  }

  void MagEdgeSpecialAVPDE::DefinePrimaryResults() {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    // === MAGNETIC VECTOR POTENTIAL ===
    shared_ptr<ResultInfo> potInfo(new ResultInfo);
    potInfo->resultType = MAG_POTENTIAL;
    potInfo->dofNames = vecComponents;
    potInfo->unit = "Vs/m";
    potInfo->definedOn = ResultInfo::ELEMENT;
    potInfo->entryType = ResultInfo::VECTOR;
    potInfo->SetFeFunction(feFunctions_[MAG_POTENTIAL]);

    feFunctions_[MAG_POTENTIAL]->SetResultInfo(potInfo);
    DefineFieldResult( feFunctions_[MAG_POTENTIAL], potInfo );

    // === COIL VOLTAGE ===
    if( useModifiedAVCurrentFormulation_ ){
      // primary result is the time integral of the coil voltage
      shared_ptr<ResultInfo> voltIntInfo(new ResultInfo);
      voltIntInfo->resultType = COIL_VOLTAGE_INTEGRAL;
      voltIntInfo->dofNames = "";
      voltIntInfo->unit = MapSolTypeToUnit(COIL_VOLTAGE_INTEGRAL);
      voltIntInfo->definedOn = ResultInfo::COIL;
      voltIntInfo->entryType = ResultInfo::SCALAR;
      voltIntInfo->SetFeFunction(feFunctions_[COIL_VOLTAGE_INTEGRAL]);

      feFunctions_[COIL_VOLTAGE_INTEGRAL]->SetResultInfo(voltIntInfo);
      DefineFieldResult( feFunctions_[COIL_VOLTAGE_INTEGRAL], voltIntInfo );

      //  coil voltage is the time derivative
      shared_ptr<ResultInfo> voltInfo(new ResultInfo);
      voltInfo->resultType = COIL_VOLTAGE;
      voltInfo->dofNames = "";
      voltInfo->unit = MapSolTypeToUnit(COIL_VOLTAGE);
      voltInfo->definedOn = ResultInfo::COIL;
      voltInfo->entryType = ResultInfo::SCALAR;
      availResults_.insert( voltInfo );
      DefineTimeDerivResult( voltInfo->resultType, 1, COIL_VOLTAGE_INTEGRAL );
    }

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[MAG_POTENTIAL] = "fluxParallel";
    idbcSolNameMap_[MAG_POTENTIAL] = "potential";

    // === PERMEABILITY ===
    shared_ptr<ResultInfo> permeability ( new ResultInfo );
    permeability->resultType = MAG_ELEM_PERMEABILITY;
    permeability->dofNames = "";
    permeability->unit = "Vs/Am";
    permeability->definedOn = ResultInfo::ELEMENT;
    permeability->entryType = ResultInfo::SCALAR;
    // In multiharmonic analysis we have complex permeability
    if(analysistype_ == MULTIHARMONIC){
      shared_ptr<CoefFunctionMulti> permFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, true));
      matCoefs_[MAG_ELEM_PERMEABILITY] = permFct;
      DefineFieldResult(permFct, permeability);
    }else{
      shared_ptr<CoefFunctionMulti> permFct(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, false));
      matCoefs_[MAG_ELEM_PERMEABILITY] = permFct;
      DefineFieldResult(permFct, permeability);
    }

  }

  void MagEdgeSpecialAVPDE::DefinePostProcResults() {

    StdVector<std::string> vecComponents;
    vecComponents = "x", "y", "z";

    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL];

    // === TIME DERIVATIVES OF PRIMARY RESULTS ===
    if( analysistype_ == TRANSIENT || analysistype_ == HARMONIC  || analysistype_ == MULTIHARMONIC) {
      // === MAGNETIC VECTOR POTENTIAL - 1ST DERIVATIVE ===
      shared_ptr<ResultInfo> aDot(new ResultInfo);
      aDot->resultType = MAG_POTENTIAL_DERIV1;
      aDot->dofNames = vecComponents;
      aDot->unit = "V/m";
      aDot->definedOn = ResultInfo::ELEMENT;
      aDot->entryType = ResultInfo::VECTOR;
      aDot->SetFeFunction(feFunctions_[MAG_POTENTIAL]);
      availResults_.insert( aDot );
      DefineTimeDerivResult( MAG_POTENTIAL_DERIV1, 1, MAG_POTENTIAL );

    }

    /*
     * in Init_NonLin() (see MagEdgePDE) we call Init_HystCoefs which
     * has to define the magnetic flux density already; defining it
     * a second time is not allowed, so check first
     */
    if(!fluxDensityDefined_){
      DefineMagFluxDensity();
    }

    PtrCoefFct bFunc = this->GetCoefFct(MAG_FLUX_DENSITY); // we need it for Lorenz force and Maxwell force

    // === MAGNETIC ENERGY ===
    shared_ptr<ResultInfo> energy(new ResultInfo);
    energy->resultType = MAG_ENERGY;
    energy->dofNames = "";
    energy->unit = "Ws";
    energy->definedOn = ResultInfo::REGION;
    energy->entryType = ResultInfo::SCALAR;
    availResults_.insert( energy );
    shared_ptr<ResultFunctor> energyFunc;
    if( isComplex_ ) {
      energyFunc.reset(new EnergyResultFunctor<Complex>(feFct, energy, 0.5));
    } else {
      energyFunc.reset(new EnergyResultFunctor<Double>(feFct, energy, 0.5));
    }
    resultFunctors_[MAG_ENERGY] = energyFunc;
    stiffFormFunctors_.insert(energyFunc);



    // === RESULTS RELATED TO TIME DERIVATIVES ===
    shared_ptr<CoefFunctionFormBased> jFunc;
    shared_ptr<CoefFunction> jPowerDensFunc;
    if( analysistype_ != STATIC ) {

      // === ELECTRIC FIELD INTENSITY ===
      shared_ptr<ResultInfo> elecIntens(new ResultInfo);
      elecIntens->resultType = ELEC_FIELD_INTENSITY;
      elecIntens->SetVectorDOFs(dim_, isaxi_);
      elecIntens->dofNames = vecComponents;
      elecIntens->unit = "V/m";
      elecIntens->definedOn = ResultInfo::ELEMENT;
      elecIntens->entryType = ResultInfo::VECTOR;
      shared_ptr<CoefFunctionMulti> elecIntensFunc( new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1, isComplex_));
      DefineFieldResult( elecIntensFunc, elecIntens );


      // === EDDY CURRENT DENSITY ===
      shared_ptr<ResultInfo> eddy(new ResultInfo);
      eddy->resultType = MAG_EDDY_CURRENT_DENSITY;
      eddy->dofNames = vecComponents;
      eddy->unit = "A/m^2";
      eddy->definedOn = ResultInfo::ELEMENT;
      eddy->entryType = ResultInfo::VECTOR;
      shared_ptr<CoefFunctionMulti> eddyFunc( new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1, isComplex_));
      DefineFieldResult( eddyFunc, eddy );

      GenerateLorentzForceResults(vecComponents, eddyFunc, bFunc, part, feFct);

      // === EDDY CURRENT (SURFACE RESULT) ===
      shared_ptr<ResultInfo> ec(new ResultInfo());
      ec->resultType = MAG_EDDY_CURRENT;
      ec->dofNames = "";
      ec->unit = "A";
      ec->definedOn = ResultInfo::SURF_REGION;
      ec->entryType = ResultInfo::SCALAR;
      availResults_.insert( ec );
      // first, create normal mapping
      shared_ptr<CoefFunctionSurf> ncd(new CoefFunctionSurf(true, 1.0, ec));
      surfCoefFcts_[ncd] = eddyFunc;
      // then, integrate values
      shared_ptr<ResultFunctor> eddyCurrentFunc;
      if( isComplex_ ) {
        eddyCurrentFunc.reset(new ResultFunctorIntegrate<Complex>(ncd, feFct, ec ) );
      } else {
        eddyCurrentFunc.reset(new ResultFunctorIntegrate<Double>(ncd, feFct, ec ) );
      }
      resultFunctors_[MAG_EDDY_CURRENT] = eddyCurrentFunc;



      // === JOULE LOSS Power DENSITY INTEGRATED OVER PERIOD  (in the harmonic case)===
      shared_ptr<ResultInfo> jld(new ResultInfo);
      jld->resultType = MAG_JOULE_LOSS_POWER_DENSITY;
      jld->dofNames = "";
      jld->unit = "W/m^3";
      jld->definedOn = ResultInfo::ELEMENT;
      jld->entryType = ResultInfo::SCALAR;
      shared_ptr<CoefFunctionMulti> jldCoef(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, isComplex_));
      DefineFieldResult( jldCoef, jld );


      // === JOULE LOSS POWER INTEGRATED OVER PERIOD (in the harmonic case) ===
      shared_ptr<ResultInfo> jldRes(new ResultInfo());
      jldRes->resultType = MAG_JOULE_LOSS_POWER;
      jldRes->dofNames = "";
      jldRes->unit = "W";
      jldRes->definedOn = ResultInfo::REGION;
      jldRes->entryType = ResultInfo::SCALAR;
      availResults_.insert( jldRes );
      shared_ptr<ResultFunctor> jldFunc;
      shared_ptr<CoefFunctionMulti> coreLossCoef(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));
      if( isComplex_ ){
        jldFunc.reset( new ResultFunctorIntegrate<Complex>(jldCoef, feFct, jldRes) );
      } else {
        jldFunc.reset( new ResultFunctorIntegrate<Double>(jldCoef, feFct, jldRes) );
      }
      resultFunctors_[MAG_JOULE_LOSS_POWER] = jldFunc;
      // it is an integrated result but we need to save the coef function
      // somewhere for the finalization
      fieldCoefs_[MAG_CORE_LOSS] = coreLossCoef;


    }// end if !static

    GenerateMaxwellForce(vecComponents, bFunc, feFct);

    GenerateVWPForce(vecComponents, bFunc, feFct);

    // === PERMEABILITY  ===
    shared_ptr<ResultInfo> perm(new ResultInfo);
    perm->resultType = MAG_ELEM_PERMEABILITY;
    perm->dofNames = "";
    perm->unit = "Vs/Am";
    perm->definedOn = ResultInfo::ELEMENT;
    perm->entryType = ResultInfo::SCALAR;

    PtrCoefFct muFunc = CoefFunction::Generate( mp_, part,
        CoefXprUnaryOp( mp_, reluc_, CoefXpr::OP_INV ) );
    DefineFieldResult( muFunc, perm );

    // === RELUCTIVITY  ===
    shared_ptr<ResultInfo> reluc(new ResultInfo);
    reluc->resultType = MAG_ELEM_RELUCTIVITY;
    reluc->dofNames = "";
    reluc->unit = "Am/Vs";
    reluc->definedOn = ResultInfo::ELEMENT;
    reluc->entryType = ResultInfo::SCALAR;
    DefineFieldResult( reluc_, reluc );


  }

  void MagEdgeSpecialAVPDE::FinalizePostProcResults() {
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;

    // Initialize standard postprocessing results
    SinglePDE::FinalizePostProcResults();


    // ============ ELECTRIC FIELD INTENSITY ============
    // Assemble coefficient function for
    // E = -\frac{\partial A}{\partial t} - \grad V
    // whereas \grad V = U * \grad V_0  and V_0 is the unit source
    // electric scalar potential (in our case it's the solution
    // of the electrokinetic problem
    // Differentiate two cases:
    //      1) No conductor/inductor: e.g. air domain, there V=0
    //         and E = -\frac{\partial A}{\partial t}
    //      2) Inductor region: E = -\frac{\partial A}{\partial t} - U * grad V_0
    //         where U is either the prescribed voltage (specialvoltage excitation)
    //         or the unknown solution quantity COIL_VOLTAGE, which is the case
    //         when we use a prescribed total current (specialcurrent excitation)

    shared_ptr<CoefFunctionMulti> elecIntensCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[ELEC_FIELD_INTENSITY]);

    StdVector<RegionIdType>::iterator regIt = regions_.Begin();
    PtrCoefFct constMinusOne = CoefFunction::Generate( mp_, Global::REAL, "-1.0");
    regIt = regions_.Begin();
    for( ; regIt != regions_.End(); ++regIt ){
      // Check if region is an inductor
      CoilRegionMap::iterator cIt = coilRegions_.find(*regIt);
      if(cIt != coilRegions_.end()){
        // It's a coil/inductor/conductor, we need to consider the electric scalar potential V

        // if specialvoltage excitation, we already scaled the grad(V_0) with U in coilCurrentDens_
        // if specialcurrent is used, we have to multiply the solution U by our own
        PtrCoefFct h1 = NULL;
        if( useModifiedAVCurrentFormulation_ ){
          h1 = CoefFunction::Generate( mp_, part,
              CoefXprVecScalOp( mp_, coilCurrentDens_[*regIt], GetCoefFct(COIL_VOLTAGE), CoefXpr::OP_MULT ) );
        }else{
          // U from voltage excitation is already included here
          if( analysistype_ == MULTIHARMONIC){
            PtrCoefFct tmp = CoefFunction::Generate( mp_, part,
                CoefXprVecScalOp(mp_, coilCurrentDens_[*regIt],
                    lexical_cast<std::string>(0.5), CoefXpr::OP_MULT));
            h1 = tmp;
          }else{
            h1 = coilCurrentDens_[*regIt];
          }
        }

        PtrCoefFct h = CoefFunction::Generate( mp_, part,
                       CoefXprBinOp( mp_, h1, GetCoefFct( MAG_POTENTIAL_DERIV1 ), CoefXpr::OP_ADD ) );
        PtrCoefFct h2 = CoefFunction::Generate( mp_, part, CoefXprVecScalOp(mp_, h, constMinusOne, CoefXpr::OP_MULT));
        elecIntensCoef->AddRegion(*regIt, h2);
      }else{
        // It's a standard region, we can skip the electric scalar potential V
        PtrCoefFct h2 = CoefFunction::Generate( mp_, part, CoefXprVecScalOp(mp_,
                                GetCoefFct( MAG_POTENTIAL_DERIV1), constMinusOne, CoefXpr::OP_MULT));
        elecIntensCoef->AddRegion(*regIt, h2);
      }// end if coil region or not
    }// loop over regions




    // ============ EDDY CURRENT DENSITY ============
    // J = conductivity * E
    shared_ptr<CoefFunctionMulti> eddyCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_EDDY_CURRENT_DENSITY]);
    regIt = regions_.Begin();
    Double conductivity = 0.0;
    for( ; regIt != regions_.End(); ++regIt ){
      materials_[*regIt]->GetScalar(conductivity,MAG_CONDUCTIVITY_SCALAR,Global::REAL);
      PtrCoefFct conduccoef = CoefFunction::Generate(mp_, Global::REAL, lexical_cast<std::string>(conductivity));
      PtrCoefFct eddyC = CoefFunction::Generate( mp_, part,
                         CoefXprVecScalOp(mp_, GetCoefFct( ELEC_FIELD_INTENSITY), conduccoef, CoefXpr::OP_MULT));
      eddyCoef->AddRegion(*regIt, eddyC);
    }// loop over regions






    // === EDDY CURRENT (JOULE) LOSS DENSITY INTEGRATED===
    /*  The Joule loss power averaved over
     *  one period T of the time history
     *  the following notation is used A(t) = Re(A exp(jwt)) , A = Ar + j Ai , and ' stands for complex conjugate
     *    P_mean = 1/T \int_0^T E(t)*J(t) dt = 1/2 Re(J E') = 1/2 Re(E J') = 1/4(J*E'+J'*E) = 1/2(Jr*Er + Ji*Ei)
     *  with the electric field E(t) and the total current density J(t).
     */
    if( analysistype_ == HARMONIC || analysistype_ == MULTIHARMONIC){
      shared_ptr<CoefFunctionMulti> eddyLossCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_JOULE_LOSS_POWER_DENSITY]);

      regIt = regions_.Begin();
      // for the sake of simplicity we should real with the total current density
      for( ; regIt != regions_.End(); ++regIt ) {
        RegionIdType actRegion = *regIt;
        // OP_MULT_CONJ computes A in B'
        PtrCoefFct conjEinJ = CoefFunction::Generate( mp_, part,
            CoefXprBinOp( mp_, GetCoefFct(MAG_EDDY_CURRENT_DENSITY),
                               GetCoefFct(ELEC_FIELD_INTENSITY), CoefXpr::OP_MULT_CONJ) );
        PtrCoefFct conjJinE = CoefFunction::Generate( mp_, part,
                    CoefXprBinOp( mp_, GetCoefFct(ELEC_FIELD_INTENSITY),
                                       GetCoefFct(MAG_EDDY_CURRENT_DENSITY), CoefXpr::OP_MULT_CONJ) );

        // to make the power consistent between harmonic and multiharmonic
        PtrCoefFct halfCoef;
        if ( analysistype_ == MULTIHARMONIC){
            halfCoef = CoefFunction::Generate( mp_, part, "0.5");
        } else {
            halfCoef = CoefFunction::Generate( mp_, part, "0.25");
        }
        PtrCoefFct tmp = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_, conjEinJ, conjJinE, CoefXpr::OP_ADD) );

        eddyLossCoef->AddRegion(actRegion, CoefFunction::Generate( mp_, part,  CoefXprBinOp(mp_, tmp, halfCoef, CoefXpr::OP_MULT) ));
      }
    }else if( analysistype_ == TRANSIENT){
      shared_ptr<CoefFunctionMulti> eddyLossCoef = dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_JOULE_LOSS_POWER_DENSITY]);

      regIt = regions_.Begin();
      // for the sake of simplicity we should real with the total current density
      for( ; regIt != regions_.End(); ++regIt ) {
        RegionIdType actRegion = *regIt;
        // OP_MULT_CONJ computes A in B'
        PtrCoefFct conjJinE = CoefFunction::Generate( mp_, part,
                    CoefXprBinOp( mp_, GetCoefFct(ELEC_FIELD_INTENSITY),
                                       GetCoefFct(MAG_EDDY_CURRENT_DENSITY), CoefXpr::OP_MULT) );
        PtrCoefFct tmp = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_, conjJinE, conjJinE, CoefXpr::OP_ADD) );

        eddyLossCoef->AddRegion(actRegion, tmp);
      }
    }

  }

  std::map<SolutionType, shared_ptr<FeSpace> >
  MagEdgeSpecialAVPDE::CreateFeSpaces(const std::string& formulation,
                             PtrParamNode infoNode ) {
    //ok default case so we create grid based approximation H1 elements
    //and standard Gauss integration
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H_CURL"){
      PtrParamNode potSpaceNode = infoNode->Get("magPotential");
      crSpaces[MAG_POTENTIAL] =
          FeSpace::CreateInstance(myParam_, potSpaceNode, FeSpace::HCURL, ptGrid_ );
      crSpaces[MAG_POTENTIAL]->Init(solStrat_);
    }else{
      EXCEPTION("The formulation " << formulation
                << "of magnetic edge PDE is not known!");
    }


    if( useModifiedAVCurrentFormulation_ ){
      PtrParamNode voltSpaceNode = infoNode->Get("coilVoltage");
      crSpaces[COIL_VOLTAGE_INTEGRAL] = FeSpace::CreateInstance(myParam_, voltSpaceNode, FeSpace::CONSTANT, ptGrid_, true);
      crSpaces[COIL_VOLTAGE_INTEGRAL]->Init(solStrat_);
    }
    return crSpaces;
  }

} // end of namespace


