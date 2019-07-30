// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MagBasePDE.hh"
#include "Utils/Coil.hh"
#include "Domain/CoefFunction/CoefFunctionHyst.hh"
#include "Driver/SolveSteps/SolveStepHyst.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Domain/CoefFunction/CoefFunctionExpression.hh"
#include "Domain/CoefFunction/CoefXpr.hh"

#include "Domain/Results/ResultFunctor.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"

#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/LinForms/BDUInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Forms/BiLinForms/BiLinWrappedLinForm.hh"
#include "Driver/MultiHarmonicDriver.hh"

namespace CoupledField
{

  MagBasePDE::MagBasePDE( Grid *aGrid, PtrParamNode paramNode,
          PtrParamNode infoNode,
          shared_ptr<SimState> simState, Domain* domain )
  : SinglePDE( aGrid, paramNode, infoNode, simState, domain) {

    isMixed_ = false;
    regionApproxSet_ = false;
    anyRegionHasConductivity_ = false;
    formulation_ = MagBasePDE::BASE;
    fluxDensityDefined_ = false;

    is3d_ = domain_->GetParamRoot()->Get("domain")->Get("geometryType")->As<std::string>() == "3d";

    // initialize material coef functions covering all regions
    reluc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, dim_, dim_, isComplex_));
    relucTensor_.reset(new CoefFunctionMulti(CoefFunction::TENSOR, dim_, dim_, false));
    conduc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1, isComplex_));

    // determine if there are coils excited by voltage
    hasVoltCoils_ = false;
    PtrParamNode coilNode = myParam_->Get( "coilList", ParamNode::PASS );
    if ( coilNode ){
      ParamNodeList coilNodes = coilNode->GetChildren();
      for( UInt k = 0; k < coilNodes.GetSize(); k++ ){
        if( coilNodes[k]->Has("source") ){
          std::string exType = coilNodes[k]->Get("source")->Get("type")->As<std::string>();
          if( exType == "voltage" ){
            hasVoltCoils_ = true;
            break;
          }
        }
      }
    }

    hysteresisCoefs_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
    // Magnetization, Polarization and Field Intensity
    // this is the point where the actual storage is created
    hysteresisPolarization_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
    hysteresisMagnetization_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
    hysteresisFieldIntensity_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));

  }

  MagBasePDE::~MagBasePDE() {
  }

  /*
   * Common functions for NODAL and EDGE formulation
   */
  void MagBasePDE::DefineSolveStep() {
    // isHysteresis_ is set in SinglePDE.cc during InitNonLin
		if(isHysteresis_){
			solveStep_ = new SolveStepHyst(*this);
		} else {
			solveStep_ = new StdSolveStep(*this);
		}
  }

  void MagBasePDE::InitNonLin() {
    SinglePDE::InitNonLin();

    InitHystCoefs();
  }

  void MagBasePDE::InitTimeStepping() {
    // Use complete implicit scheme
    Double gamma = 1.0;
    GLMScheme * scheme = new Trapezoidal(gamma);

    TimeSchemeGLM::NonLinType nlType = (nonLin_ || isHysteresis_)? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType) );
    feFunctions_[MAG_POTENTIAL]->SetTimeScheme(myScheme);

    // Important: Create a new time scheme for each additional feFunction
    // NEW: from NACS - copy stepping scheme from mag potential
    shared_ptr<TimeSchemeGLM> mainScheme = dynamic_pointer_cast<TimeSchemeGLM>(
            feFunctions_[MAG_POTENTIAL]->GetTimeScheme());
    assert(mainScheme);

    if( hasVoltCoils_ ){
      shared_ptr<BaseTimeScheme> myScheme2(new TimeSchemeGLM(*mainScheme));
      feFunctions_[COIL_CURRENT]->SetTimeScheme(myScheme2);
    }

    if( isMixed_ ){
      shared_ptr<BaseTimeScheme> myScheme3(new TimeSchemeGLM(*mainScheme));
      feFunctions_[ELEC_POTENTIAL]->SetTimeScheme(myScheme3);
    }
  }

  void MagBasePDE::ReadSpecialBCs() {
    // --------------------------------------------------------------------
    //   Get information about coils and open files for measurement coils
    // --------------------------------------------------------------------
    ReadCoils();
  }

  /*
   * Coil specific functions
   */
  shared_ptr<Coil> MagBasePDE::GetCoilById(const Coil::IdType& id) {
    return coils_.at(id);
  }

  void MagBasePDE::ReadCoils() {

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

      // Adjust printing of coil information to info node
      // WARN("Adapt printing of coils to InfoNode");
    }
  }

  /*
   * Copy-paste from DefinePostProcessingResults
   * > setup as stand-alone so that it can be called in InitHystCoefs which requires the FluxDensity to be defined
   */
  void MagBasePDE::DefineMagFluxDensity(){
    StdVector<std::string> vecComponents, aVecComponents;
    if( dim_ == 3 ) {
      vecComponents = "x", "y", "z";
      aVecComponents = "x", "y", "z";
    }
    else if( isaxi_ ) {
      vecComponents = "r", "z";
      aVecComponents = "phi";
    }
    else {
      vecComponents = "x", "y";
      aVecComponents = "z";
    }

    shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL];

    // === MAGNETIC FLUX DENSITY ===
    shared_ptr<ResultInfo> fluxDens(new ResultInfo);
    fluxDens->resultType = MAG_FLUX_DENSITY;
    fluxDens->dofNames = vecComponents;
    fluxDens->unit = "Vs/m^2";
    fluxDens->definedOn = ResultInfo::ELEMENT;
    fluxDens->entryType = ResultInfo::VECTOR;
    availResults_.insert( fluxDens );
    shared_ptr<CoefFunctionFormBased> bFunc;
    if( isComplex_ ) {
      bFunc.reset(new CoefFunctionBOp<Complex>(feFct, fluxDens));
    } else {
      bFunc.reset(new CoefFunctionBOp<Double>(feFct, fluxDens));
    }
    DefineFieldResult( bFunc, fluxDens );
    stiffFormCoefs_.insert(bFunc);

    fluxDensityDefined_ = true;
  }


  /*
   * Hysteresis specific functions
   */
  void MagBasePDE::InitHystCoefs() {
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;

    SubTensorType tensorType;
    if ( dim_ == 3  || subType_ == "2.5d" ) {
      tensorType = FULL;
    } else {
      if ( isaxi_ == true ) {
        tensorType = AXI;
      } else {
        // 2d: plane case
        tensorType = PLANE_STRAIN;
      }
    }

    std::map<RegionIdType, BaseMaterial*>::iterator it;
    shared_ptr<FeSpace> mySpace = feFunctions_[MAG_POTENTIAL]->GetFeSpace();
    // we set the dependency of the hyst operator to the actual
    // quantity that is computed (B) although the standard Preisach
    // as well as the VectorPreisach model need H as input
    if(!fluxDensityDefined_){
      DefineMagFluxDensity();
    }
    PtrCoefFct magFieldCoef = this->GetCoefFct(MAG_FLUX_DENSITY);

    for ( it = materials_.begin(); it != materials_.end(); it++ ) {

      // Set current region and material
      actRegion = it->first;
      actSDMat = it->second;

      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );

      // ==========================================================
      // New implementation
      // ==========================================================
      // --- Set the approximation for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId, integId);

      SDLists_[actRegion] = actSDList;

      StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[actRegion];
      if ( nonLinTypes.Find(HYSTERESIS) != -1 ){
//        shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
//        actSDList->SetRegion( actRegion );

        // check if hysteresis is defined in material file, too!
        std::string hystType;
        actSDMat->GetScalar(hystType, HYST_MODEL);

        if(hystType == "none"){
          std::string warnmsg = "Hysteresis set on region " + regionName + " but no hysteresis model was defined in mat file. Skip.";
          regionNonLinTypes_[actRegion] = NO_NONLINEARITY;
          WARN(warnmsg);
        } else {

          // please note:
          //  in magnetics, the hysteresis model is supposed to return the
          //  magnetic polarizaiton J_P = mu*M
          //  (in older versions, the magnetization M was returned!)
          PtrCoefFct hystPol(new CoefFunctionHyst( actSDMat, actSDList,
                  magFieldCoef,tensorType,MAG_RELUCTIVITY,mySpace));
//          std::cout << "hysteresisCoefs_->AddRegion( actRegion, hystPol);" << std::endl;
          hysteresisCoefs_->AddRegion( actRegion, hystPol);

          PtrCoefFct hystOutput = hystPol->GenerateOutputCoefFnc("MagPolarization");
          hysteresisPolarization_->AddRegion( actRegion, hystOutput);

          PtrCoefFct hystOutput2 = hystPol->GenerateOutputCoefFnc("MagMagnetization");
          hysteresisMagnetization_->AddRegion( actRegion, hystOutput2);

          PtrCoefFct hystOutput3 = hystPol->GenerateOutputCoefFnc("MagFieldIntensityHyst");
          hysteresisFieldIntensity_->AddRegion( actRegion, hystOutput3);

        }
      }
    }
    regionApproxSet_ = true;
  }

  void MagBasePDE::DefineCoilIntegrators(Double scaling) {


    // ============================
    // COIL INTEGRATORS >  based on EDGE form
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

      if( actCoil.sourceType_ == Coil::SPECIALCURRENT ||
          actCoil.sourceType_ == Coil::SPECIALVOLTAGE){
        EXCEPTION("For specialvoltage or specialcurrent excitation, please use specialA-V formulation!!");
      }

      if(( actCoil.sourceType_ == Coil::CURRENT )||
         ( actCoil.sourceType_ == Coil::CURRENT_MULTHARM )||
         ( actCoil.sourceType_ == Coil::EXTERNAL )) {

        if( ( actCoil.sourceType_ == Coil::CURRENT_MULTHARM ) && (formulation_ != MagBasePDE::EDGE )){
          EXCEPTION("Multiharmonic coils only tested/used for EDGE formulation so far");
        }

        /*
        =====================================================
         1) CURRENT driven coils OR EXTERNAL current density

         Ref: M. Kaltenbacher, Numer. Sim. of. Mech.
              Sens. and Act., 2nd edition, p. 131ff
        =====================================================
        */

        for( partIt = actCoil.parts_.begin();
            partIt != actCoil.parts_.end();
            partIt++ ) {
          Coil::Part & actPart = *(partIt->second);
          RegionIdType actRegion = partIt->first;
          shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
          actSDList->SetRegion( actRegion );

          /*
            generate source current vector
            for the non-multiharmonic case it is a simple PtrCoefFct type
            but due to the flexibility to use multiharmonic excitation,
            we need a map of PtrCoefFct's, where the key is the harmonic number.
            So in the non-multiharmonic case, we simply have one key "0"
          */
          std::map<UInt, PtrCoefFct> jFct;
          if( actCoil.sourceType_ == Coil::CURRENT ){
            CoefXprVecScalOp iVec = CoefXprVecScalOp(mp_, actPart.jUnitVec, actCoil.srcVal_,
                                                   CoefXpr::OP_MULT);
            PtrCoefFct iFct = CoefFunction::Generate(mp_, part, iVec);

            CoefXprVecScalOp jVec = CoefXprVecScalOp(mp_, iFct, boost::lexical_cast<std::string>(actPart.wireCrossSect),
                                                   CoefXpr::OP_DIV);
            jFct[0] = CoefFunction::Generate(mp_, part, jVec);
          } else if( actCoil.sourceType_ == Coil::CURRENT_MULTHARM ){

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
          } else {
            jFct[0] = coilPartsExtJ_[partIt->second];
          }


          // This switch is necessary because we need to select which harmonic component
          // we want to assemble in the rhs
          if( actCoil.sourceType_ == Coil::CURRENT_MULTHARM ){
            // Multiharmonic Case
            LinearForm* curInt;
            for(auto& h : jFct){
              coilCurrentDens_[actRegion] = h.second;
              curInt = GetCurrentDensityInt( scaling, h.second );
//              if( isComplex_ ) {
//                curInt = new BUIntegrator<Complex>( new IdentityOperator<FeHCurl,3,1,Complex>(),
//                                                    1.0, h.second, updatedGeo_);
//              }
//              else {
//                curInt = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,3,1,Double>(),
//                                                   1.0, h.second, updatedGeo_);
//              }
              curInt->SetName("CoilIntegrator");
              curInt->SetHarm(h.first);
              LinearFormContext * coilContext = new LinearFormContext( curInt );
              coilContext->SetEntities( actSDList );
              coilContext->SetFeFunction( feFunc );
              assemble_->AddLinearForm( coilContext );
            }
          }else{
            // Classic Case
            LinearForm* curInt;
            coilCurrentDens_[actRegion] = jFct[0];
            curInt = GetCurrentDensityInt( scaling, jFct[0] );
//            if( isComplex_ ) {
//              curInt = new BUIntegrator<Complex>( new IdentityOperator<FeHCurl,3,1,Complex>(),
//                                                  1.0, jFct[0], updatedGeo_);
//            }
//            else {
//              curInt = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,3,1,Double>(),
//                                                 1.0, jFct[0], updatedGeo_);
//            }
            curInt->SetName("CoilIntegrator");
            LinearFormContext * coilContext = new LinearFormContext( curInt );
            coilContext->SetEntities( actSDList );
            coilContext->SetFeFunction( feFunc );
            assemble_->AddLinearForm( coilContext );

          }
        } // loop: parts

      }else{

        /*
        ============================================
         2) VOLTAGE driven coils

         Ref: M. Kaltenbacher, Numer. Sim. of. Mech.
              Sens. and Act., 2nd edition, p. 211ff
        ============================================
         The coupled equation system in this case looks like
           ( M_A     0 ) ( A_dot ) + ( K_A -f_A ) ( A ) = ( 0 )
           ( (f_A)^T 0 ) ( i_dot )   ( 0     R  ) ( i )   ( u )
        */

        std::string totRstr = "";
        shared_ptr<CoilList> singleCoilList( new CoilList( ptGrid_ ) );
        singleCoilList->AddCoil( coilIt->second );
        feFunctions_[COIL_CURRENT]->AddEntityList( singleCoilList );

        for( partIt = actCoil.parts_.begin();
             partIt != actCoil.parts_.end();
             partIt++ ) {

          Coil::Part & actPart = *(partIt->second);

          if( totRstr.empty() ){
            totRstr = actPart.resistance;
          } else {
            totRstr += " + " + actPart.resistance;
          }

          CoefXprVecScalOp eJscaledOp = CoefXprVecScalOp( mp_, actPart.jUnitVec,
              boost::lexical_cast<std::string>(actPart.wireCrossSect), CoefXpr::OP_DIV );
          PtrCoefFct eJscaled = CoefFunction::Generate( mp_, part, eJscaledOp );

          shared_ptr<ElemList> actSDList( new ElemList( ptGrid_ ) );
          RegionIdType actRegion = partIt->first;
          actSDList->SetRegion( actRegion );

          // implementation of coil current density is difficult because of FeSpaceConst;
          // it looks simple: J = I/Gamma_c, where Gamma_c is the coil cross section;
          // 1) but the FeSpaceConst does not have elements and the CoefFunction asks
          //    for elements in order to evaluate its expression (FeFunction::GetScalar);
          //    although the origin of this problem does not seem to be the FeFunction;
          //    the coil result coil current works properly
          // 2) the automatic calculation of the cross section of the coil
          //    must be implemented because we need the number of turns or the coil cross
          //    section (only the winding cross section is not enough!)
          //    or
          //    number of turns or coil cross section could be additionally specified in xml;
          // However, the effects are not severe: Current density results are not available
          // for coils with voltage sources, which is generally not interesting anyway.
          // With the 2 points resolved, the code could look like:
          /*CoefXprVecScalOp testOp = CoefXprVecScalOp( mp_, actPart.jUnitVec,
            GetCoefFct( COIL_CURRENT ), CoefXpr::OP_MULT );
          PtrCoefFct test = CoefFunction::Generate( mp_, part, testOp );*/
          // now the division by the cross section would be necessary
          // the unit vector of the current density is added as dummy so that the other
          // current density results do not get mixed up
          coilCurrentDens_[actRegion] = actPart.jUnitVec;

          // === -f_A ===
          LinearForm* psiDotInt = GetCurrentDensityInt(-1.0*scaling, eJscaled );
//          LinearForm* psiDotInt;
//          if( isComplex_ ) {
//            psiDotInt = new BUIntegrator<Complex>( new IdentityOperator<FeHCurl,3,1,Complex>(),
//                -1.0, eJscaled, updatedGeo_);
//          } else {
//            psiDotInt = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,3,1,Double>(),
//                -1.0, eJscaled, updatedGeo_);
//          }
          psiDotInt->SetName("CoilVoltCouplInt");

          bool assembleTransposed = false;
          BiLinearForm* pseudoBiLin = new BiLinWrappedLinForm( psiDotInt, assembleTransposed );
          BiLinFormContext* voltCoilContext = new BiLinFormContext( pseudoBiLin, STIFFNESS );
          voltCoilContext->SetEntities( actSDList, singleCoilList );
          voltCoilContext->SetFeFunctions( feFunc, feFunctions_[COIL_CURRENT] );
          voltCoilContext->SetCounterPart(false);
          assemble_->AddBiLinearForm( voltCoilContext );

          // === (f_A)^T ===
          if( analysistype_ != STATIC ){
            LinearForm* psiDotIntT = GetCurrentDensityInt(scaling, eJscaled);
//            LinearForm* psiDotIntT;
//            if( isComplex_ ) {
//              psiDotIntT = new BUIntegrator<Complex>( new IdentityOperator<FeHCurl,3,1,Complex>(),
//                  1.0, eJscaled, updatedGeo_);
//            } else {
//              psiDotIntT = new BUIntegrator<Double>( new IdentityOperator<FeHCurl,3,1,Double>(),
//                  1.0, eJscaled, updatedGeo_);
//            }
            psiDotIntT->SetName("CoilVoltCouplIntTransposed");

            assembleTransposed = true;
            BiLinearForm* pseudoBiLinT = new BiLinWrappedLinForm( psiDotIntT, assembleTransposed );
            BiLinFormContext* voltCoilContextT = new BiLinFormContext( pseudoBiLinT, DAMPING );
            voltCoilContextT->SetEntities( singleCoilList, actSDList );
            voltCoilContextT->SetFeFunctions( feFunctions_[COIL_CURRENT], feFunc );
            voltCoilContextT->SetCounterPart(false);
            assemble_->AddBiLinearForm( voltCoilContextT );
          }

        } // loop: parts

        // === R ===
        PtrCoefFct totR = CoefFunction::Generate( mp_, part, totRstr, "0.0" );
        LinearForm* totRint = new SingleEntryInt( totR );
        totRint->SetName( "CoilResistanceInt" );
        BiLinearForm* totRBiLin = new BiLinWrappedLinForm( totRint, false );
        BiLinFormContext* totRcontext = new BiLinFormContext( totRBiLin, STIFFNESS );
        totRcontext->SetEntities( singleCoilList, singleCoilList );
        totRcontext->SetFeFunctions( feFunctions_[COIL_CURRENT], feFunctions_[COIL_CURRENT] );
        totRcontext->SetCounterPart(false);
        assemble_->AddBiLinearForm( totRcontext );

        // === u ===
        LinearForm* voltInt = new SingleEntryInt( actCoil.srcVal_ );
        voltInt->SetName( "CoilVoltageLoadInt" );
        LinearFormContext* voltContext = new LinearFormContext( voltInt );
        voltContext->SetEntities( singleCoilList );
        voltContext->SetFeFunction( feFunctions_[COIL_CURRENT] );
        assemble_->AddLinearForm( voltContext );

      } // if: current / voltage driven
    } // loop: coils
  }

} // end of namespace

//
//
//#include <cmath>
//#include <boost/tr1/type_traits.hpp>
//
//#include "Driver/Assemble.hh"
//
//// include fespaces
//#include "FeBasis/H1/H1Elems.hh"
//#include "FeBasis/HCurl/HCurlElems.hh"
//
//// new integrator concept
//#include "Forms/Operators/IdentityOperator.hh"
//#include "Forms/LinForms/BDUInt.hh"
//#include "Forms/LinForms/BUInt.hh"
//#include "Forms/LinForms/SingleEntryInt.hh"
//#include "Forms/BiLinForms/BiLinWrappedLinForm.hh"
//
//#include "Domain/ElemMapping/EntityLists.hh"
//#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
//
//// new postprocessing concept
////#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
////#include "Domain/CoefFunction/CoefFunctionExpression.hh"
//#include "Domain/CoefFunction/CoefXpr.hh"
//#include <Domain/CoefFunction/FeFunctionBiotSavart.hh>
//
//
//namespace NACS {
//
//MagBasePDE::MagBasePDE( Grid *aGrid, PtrParamNode paramNode,
//                        PtrParamNode infoNode,
//                        shared_ptr<SimState> simState, Domain* domain )
// : SinglePDE( aGrid, paramNode, infoNode, simState, domain)
//{
//
//  //! Always use updated Lagrangian formulation
//  updatedGeo_        = true;
//  pdematerialclass_  = BM::ELECTROMAGNETIC;
//
//  hasVoltCoils_ = false;
//
//  // check if we have a 3d setup
//  is3d_ = domain_->GetParamRoot()->Get("domain")->Get("geometryType")->As<std::string>() == "3d";
//
//  reluc_.reset(new CoefFunctionMulti(CoefFunction::TENSOR, dim_, dim_, isComplex_));
//  conduc_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1,1, isComplex_));
//  magnetization_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_, 1, isComplex_, true));
//
//  // Read the coils already here, as we need the information about
//  // voltage-loaded coils already during creation of FeSpaces and integrators
//  coilInfo_.reset(new ParamNode(ParamNode::APPEND, ParamNode::ELEMENT, false));
//  coilInfo_->SetName("coilList");
//  ReadCoils();
//}
//
//MagBasePDE::~MagBasePDE() {
//
//}
//
//shared_ptr<Coil> MagBasePDE::GetCoilById(const Coil::IdType& id) {
//  return coils_.at(id);
//}
//
//void MagBasePDE::InitNonLin() {
//  SinglePDE::InitNonLin();
//}
//
//void MagBasePDE::DefineCoilIntegrators() {
//
//  // obtain space / feFunction for magnetic vector potential
//  shared_ptr<BaseFeFunction> myFct = feFunctions_[MAG_POTENTIAL];
//
//  // ============================
//  // COIL INTEGRATORS
//  // ============================
//  Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
//  std::map<Coil::IdType, shared_ptr<Coil> >::iterator coilIt;
//  coilIt = coils_.begin();
//  for( ; coilIt != coils_.end(); coilIt++ ) {
//    Coil& actCoil = *(coilIt->second);
//    // run over all parts
//    std::map<RegionIdType,shared_ptr<Coil::Part> >::iterator partIt;
//    partIt = actCoil.parts_.begin();
//    if( actCoil.sourceType_ == Coil::CURRENT ) {
//      /*
//         ============================================
//          1) CURRENT driven coils
//
//          Ref: M. Kaltenbacher, Numer. Sim. of. Mech.
//               Sens. and Act., 2nd edition, p. 131ff
//         ============================================
//       */
//      for( partIt = actCoil.parts_.begin();
//          partIt != actCoil.parts_.end();
//          partIt++ ) {
//        Coil::Part & actPart = *(partIt->second);
//        RegionIdType actRegion = partIt->first;
//        shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
//        actSDList->SetRegion( actRegion );
//
//
//
//        // generate source current vector by multiplying J with I
//        CoefXprVecScalOp jVec = CoefXprVecScalOp(mp_, actPart.jUnitVec, actCoil.srcVal_,
//                                                 CoefXpr::OP_MULT);
//        PtrCoefFct jFct = CoefFunction::Generate(mp_, part, jVec);
//
//        LinearForm* curInt = GetCurrentDensityInt( 1.0, jFct );
//
//        // store coil current density in map
//        coilCurrentDens_[actRegion] = jFct;
//
//        curInt->SetName("CoilIntegrator");
//        LinearFormContext * coilContext =
//            new LinearFormContext( curInt );
//        coilContext->SetEntities( actSDList );
//        coilContext->SetFeFunction( myFct );
//        assemble_->AddLinearForm( coilContext );
//
//        if (computeMode_ != COMPUTE_FULL && !isComplex_) {
//          ScalarInput scin;
//          scin.name = "Current_Coil_" + actCoil.coilId_ + "_" +
//                      ptGrid_->GetRegion().ToString(actRegion);
//          scin.type = COIL_CURRENT;
//          scin.integrator = GetCurrentDensityInt(1.0, actPart.jUnitVec);
//          scin.context.reset(new LinearFormContext(scin.integrator));
//          scin.context->SetEntities(actSDList);
//          scin.context->SetFeFunction(myFct);
//          scalarInputs_.push_back(scin);
//        }
//      } // loop: parts
//    } else {
//
//      /*
//         ============================================
//          2) VOLTAGE driven coils
//
//          Ref: M. Kaltenbacher, Numer. Sim. of. Mech.
//               Sens. and Act., 2nd edition, p. 211ff
//         ============================================
//          The coupled equation system in this case looks like
//            ( M_A     0 ) ( A_dot ) + ( K_A -f_A  ) ( A ) = ( 0 )
//            ( (f_A)^T 0 ) ( i_dot )   ( 0     R   ) ( i )   ( u )
//       */
//
//      shared_ptr<CoilList> singleCoilList( new CoilList( ptGrid_ ) );
//      singleCoilList->AddCoil( coilIt->second );
//      feFunctions_[COIL_CURRENT]->AddEntityList( singleCoilList );
//      std::string totRstr = coilIt->second->resistance;
//      for( partIt = actCoil.parts_.begin();
//          partIt != actCoil.parts_.end();
//          partIt++ ) {
//
//        Coil::Part & actPart = *(partIt->second);
//        shared_ptr<ElemList> actSDList( new ElemList( ptGrid_ ) );
//        RegionIdType actRegion = partIt->first;
//        actSDList->SetRegion( actRegion );
//
//        // === -f_A ===
//        LinearForm* psiDotInt = GetCurrentDensityInt(-1.0, actPart.jUnitVec );
//        psiDotInt->SetName("CoilVoltCouplInt");
//
//        bool assembleTransposed = false;
//        BiLinearForm* pseudoBiLin = new BiLinWrappedLinForm( psiDotInt, assembleTransposed );
//        BiLinFormContext* voltCoilContext = new BiLinFormContext( pseudoBiLin, STIFFNESS );
//        voltCoilContext->SetEntities( actSDList, singleCoilList );
//        voltCoilContext->SetFeFunctions( myFct, feFunctions_[COIL_CURRENT] );
//        voltCoilContext->SetCounterPart(false);
//        assemble_->AddBiLinearForm( voltCoilContext );
//
//        // === (f_A)^T ===
//        if( analysistype_ != STATIC ){
//          LinearForm* psiDotIntT = GetCurrentDensityInt(1.0, actPart.jUnitVec);
//          psiDotIntT->SetName("CoilVoltCouplIntTransposed");
//
//          assembleTransposed = true;
//          BiLinearForm* pseudoBiLinT = new BiLinWrappedLinForm( psiDotIntT, assembleTransposed );
//          BiLinFormContext* voltCoilContextT = new BiLinFormContext( pseudoBiLinT, DAMPING );
//          voltCoilContextT->SetEntities( singleCoilList, actSDList );
//          voltCoilContextT->SetFeFunctions( feFunctions_[COIL_CURRENT], myFct );
//          voltCoilContextT->SetCounterPart(false);
//          assemble_->AddBiLinearForm( voltCoilContextT );
//        }
//
//        // CoefFunction for current density of this coil
//        PtrCoefFct jFct;
//        if( isComplex_ ) {
//          jFct.reset( new CoefFunctionJVoltCoil<Complex>(ptGrid_, coilIt->second,
//                                                         feFunctions_[COIL_CURRENT]));
//        } else {
//          jFct.reset( new CoefFunctionJVoltCoil<Double>( ptGrid_, coilIt->second,
//                                                        feFunctions_[COIL_CURRENT]));
//        }
//        coilCurrentDens_[actRegion] = jFct;
//
//
//      } // loop: parts
//
//      // === R ===
//      PtrCoefFct totR = CoefFunction::Generate( mp_, part, totRstr, "0.0" );
//      LinearForm* totRint = new SingleEntryInt( totR );
//      totRint->SetName( "CoilResistanceInt" );
//      BiLinearForm* totRBiLin = new BiLinWrappedLinForm( totRint, false );
//      BiLinFormContext* totRcontext = new BiLinFormContext( totRBiLin, STIFFNESS );
//      totRcontext->SetEntities( singleCoilList, singleCoilList );
//      totRcontext->SetFeFunctions( feFunctions_[COIL_CURRENT], feFunctions_[COIL_CURRENT] );
//      totRcontext->SetCounterPart(false);
//      assemble_->AddBiLinearForm( totRcontext );
//
//      // === u ===
//      LinearForm* voltInt = new SingleEntryInt( actCoil.srcVal_ );
//      voltInt->SetName( "CoilVoltageLoadInt" );
//      LinearFormContext* voltContext = new LinearFormContext( voltInt );
//      voltContext->SetEntities( singleCoilList );
//      voltContext->SetFeFunction( feFunctions_[COIL_CURRENT] );
//      assemble_->AddLinearForm( voltContext );
//
//    } //if: current / voltage driven
//  } // loop: coils
//
//  // append existing coil info node to the PDE's info node
//  infoNode_->AddChildNode(coilInfo_);
//}
//
//
//
//void MagBasePDE::DefineCoilFeSpace(std::map<SolutionType, shared_ptr<FeSpace> >& crSpaces,
//                                  PtrParamNode infoNode) {
//  if( hasVoltCoils_ ){
//    PtrParamNode voltSpaceNode = infoNode->Get("coilCurrent");
//    crSpaces[COIL_CURRENT] =
//        FeSpace::CreateInstance(myParam_, voltSpaceNode, FeSpace::CONSTANT, ptGrid_);
//    crSpaces[COIL_CURRENT]->Init(solStrat_);
//  }
//}
//
//
//void MagBasePDE::DefinePrimaryCoilResults() {
//
//
////  if( hasVoltCoils_ ){
////    // === COIL CURRENT ===
////    shared_ptr<ResultInfo> currentInfo(new ResultInfo);
////    currentInfo->resultType = COIL_CURRENT;
////    currentInfo->dofNames = "";
////    currentInfo->unit = "A";
////    currentInfo->definedOn = ResultInfo::COIL;
////    currentInfo->entryType = ResultInfo::SCALAR;
////
////
////    feFunctions_[COIL_CURRENT]->SetResultInfo(currentInfo);
////    DefineFieldResult( feFunctions_[COIL_CURRENT], currentInfo );
////  }
//}
//
//void MagBasePDE::DefinePostProcCoilResults() {
//
//
//  StdVector<std::string> vecComponents;
//  if( dim_ == 3 ) {
//    vecComponents = "x", "y", "z";
//  }
//  else if( isaxi_ ) {
//    vecComponents = "phi";
//  }
//  else {
//    vecComponents = "z";
//  }
//
//  // === COIL CURRENT DENSITY ===
//  shared_ptr<ResultInfo> ccd(new ResultInfo);
//  ccd->resultType = MAG_COIL_CURRENT_DENSITY;
//  ccd->dofNames = vecComponents;
//  ccd->unit = "A/m^2";
//  ccd->definedOn = ResultInfo::ELEMENT;
//  ccd->entryType = ResultInfo::VECTOR;
//  availResults_.insert( ccd );
//  shared_ptr<CoefFunctionMulti> ccdCoef(new CoefFunctionMulti(CoefFunction::VECTOR, vecComponents.GetSize(),1,
//                                                              isComplex_));
////  // loop over all coil coefficients and add contribution to coef
////  std::map<RegionIdType, PtrCoefFct>::iterator coilIt = coilCurrentDens_.begin();
////  for( ; coilIt != coilCurrentDens_.end(); ++coilIt ) {
////    ccdCoef->AddRegion( coilIt->first, coilIt->second);
////  }
//  DefineFieldResult( ccdCoef, ccd );
//
//  // Obtain CoefFunctionf for potential and time derivative
//  PtrCoefFct coefA = GetCoefFct(MAG_POTENTIAL);
//  PtrCoefFct coefADot;
//  if( analysistype_ != STATIC ) {
//    coefADot = GetCoefFct(MAG_POTENTIAL_DERIV1);
//  }
//  shared_ptr<BaseFeFunction> myFct = feFunctions_[MAG_POTENTIAL];
//  shared_ptr<BaseFeFunction> currentFunc;
//  if( hasVoltCoils_) {
//    currentFunc = feFunctions_[COIL_CURRENT];
//  }
//
//  // === COIL CURRENT ===
//  shared_ptr<ResultInfo> cc (new ResultInfo);
//  cc->resultType = COIL_CURRENT;
//  cc->dofNames = "";
//  cc->unit = "A";
//  cc->definedOn = ResultInfo::COIL;
//  cc->entryType = ResultInfo::SCALAR;
//  availResults_.insert( cc );
//  shared_ptr<ResultFunctor> ccFunctor;
//  if( isComplex_ ) {
//    ccFunctor.reset(new CoilResultFunctor<Complex>(mp_, ptGrid_, myFct,cc, coils_, coefA,
//                                                   coefADot, currentFunc, conduc_ ));
//
//  } else {
//    ccFunctor.reset(new CoilResultFunctor<Double>(mp_, ptGrid_,myFct,cc, coils_, coefA,
//                                                  coefADot, currentFunc, conduc_ ));
//  }
//  resultFunctors_[COIL_CURRENT] = ccFunctor;
//
//  if(hasVoltCoils_) {
//    feFunctions_[COIL_CURRENT]->SetResultInfo(cc);
//  }
//
//
//  // === COIL LINKED FLUX ===
//  shared_ptr<ResultInfo> lFlux(new ResultInfo);
//  lFlux->resultType = COIL_LINKED_FLUX;
//  lFlux->dofNames = "";
//  lFlux->unit = "Vs";
//  lFlux->definedOn = ResultInfo::COIL;
//  lFlux->entryType = ResultInfo::SCALAR;
//  availResults_.insert( lFlux );
//  shared_ptr<ResultFunctor> lfFunctor;
//  if( isComplex_ ) {
//    lfFunctor.reset(new CoilResultFunctor<Complex>(mp_, ptGrid_,myFct, lFlux, coils_, coefA, coefADot,
//                                                   currentFunc, conduc_ ));
//
//  } else {
//    lfFunctor.reset(new CoilResultFunctor<Double>(mp_,ptGrid_, myFct, lFlux, coils_, coefA, coefADot,
//                                                  currentFunc, conduc_ ));
//  }
//  resultFunctors_[COIL_LINKED_FLUX] = lfFunctor;
//
//
//  // === COIL INDUCTANCE ===
//  shared_ptr<ResultInfo> induct(new ResultInfo);
//  induct->resultType = COIL_INDUCTANCE;
//  induct->dofNames = "";
//  induct->unit = "Vs/A";
//  induct->definedOn = ResultInfo::COIL;
//  induct->entryType = ResultInfo::SCALAR;
//  availResults_.insert( induct );
//  shared_ptr<ResultFunctor> ciFunctor;
//  if( isComplex_ ) {
//    ciFunctor.reset(new CoilResultFunctor<Complex>(mp_, ptGrid_,myFct, induct, coils_, coefA, coefADot,
//                                                   currentFunc, conduc_ ));
//
//  } else {
//    ciFunctor.reset(new CoilResultFunctor<Double>(mp_,ptGrid_, myFct, induct, coils_, coefA, coefADot,
//                                                  currentFunc, conduc_ ));
//  }
//  resultFunctors_[COIL_INDUCTANCE] = ciFunctor;
//
//  // === COIL INDUCED VOLTAGE ===
//  shared_ptr<ResultInfo> volt(new ResultInfo);
//  volt->resultType = COIL_INDUCED_VOLTAGE;
//  volt->dofNames = "";
//  volt->unit = "V";
//  volt->definedOn = ResultInfo::COIL;
//  volt->entryType = ResultInfo::SCALAR;
//  availResults_.insert( volt );
//  shared_ptr<ResultFunctor> ivFunctor;
//  if( isComplex_ ) {
//    ivFunctor.reset(new CoilResultFunctor<Complex>(mp_, ptGrid_,myFct, volt, coils_, coefA, coefADot,
//                                                   currentFunc, conduc_ ));
//
//  } else {
//    ivFunctor.reset(new CoilResultFunctor<Double>(mp_,ptGrid_, myFct, volt, coils_, coefA, coefADot,
//                                                  currentFunc, conduc_ ));
//  }
//  resultFunctors_[COIL_INDUCED_VOLTAGE] = ivFunctor;
//}
//
//void MagBasePDE::FinalizePostProcCoilResults() {
//  // === COIL CURRENT DENSITY ===
//  shared_ptr<CoefFunctionMulti> ccdCoef =
//      dynamic_pointer_cast<CoefFunctionMulti>(fieldCoefs_[MAG_COIL_CURRENT_DENSITY]);
//  // loop over all coil coefficients and add contribution to coef
//  std::map<RegionIdType, PtrCoefFct>::iterator coilIt = coilCurrentDens_.begin();
//  for( ; coilIt != coilCurrentDens_.end(); ++coilIt ) {
//    ccdCoef->AddRegion( coilIt->first, coilIt->second);
//  }
//}
//
//void MagBasePDE::ReadCoils() {
//  // Check if the element "coils" is present at all.
//  // Otherwise leave
//  PtrParamNode coilNode = myParam_->Get( "coilList", ParamNode::PASS );
//  if ( !coilNode )
//    return;
//
//  // Get ParamNodes of standard coils
//  ParamNodeList coilNodes = coilNode->GetList("coil");
//
//  // Trigger reading in of definitions
//  if( coilNodes.GetSize() > 0 ) {
//    for( UInt i = 0; i < coilNodes.GetSize(); i++ ) {
//      Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
//
//      // get coil and id
//      std::string coilId = coilNodes[i]->Get("id")->As<std::string>();
//
//      // Check if coil with same ID already exists
//      if( coils_.find(coilId) != coils_.end() ) {
//        EXCEPTION("A coil with ID '" << coilId << "' was already defined.")
//      }
//
//      // Create new coil
//      shared_ptr<Coil> actCoil( new Coil( coilNodes[i], coilInfo_,
//                                          ptGrid_, mp_, part) );
//      coils_[coilId] = actCoil;
//
//      // Associate mapping of coil parts with regions
//      std::map<RegionIdType, shared_ptr<Coil::Part> >::const_iterator it;
//      for( it = actCoil->parts_.begin(); it != actCoil->parts_.end(); it++ ) {
//        coilRegions_[it->first] = actCoil;
//      }
//
//      // check if coil is voltage-driven
//      if( actCoil->sourceType_ == Coil::VOLTAGE ) {
//        hasVoltCoils_ = true;
//      }
//    }
//  }
//
//  if (coilNode->Has("biotSavartCoil")) {
//    // initialize Biot-Savart function
//    if (IsComplex()) {
//      FeFunctionBiotSavart<Complex> *bs =
//          new FeFunctionBiotSavart<Complex>(coilNode, coilInfo_, ptGrid_, mp_);
//      bs->SetFormulation(FeFunctionBiotSavart<Complex>::VEC_POT);
//      biotSavart_.reset(bs);
//    }
//    else {
//      FeFunctionBiotSavart<Double> *bs =
//          new FeFunctionBiotSavart<Double>(coilNode, coilInfo_, ptGrid_, mp_);
//      bs->SetFormulation(FeFunctionBiotSavart<Double>::VEC_POT);
//      biotSavart_.reset(bs);
//    }
//    biotSavart_->SetPDE(this);
//  }
//}
//
//void MagBasePDE::InitCoilTimeStepping() {
//  if( hasVoltCoils_ ){
//    // Important: Create a new time scheme just for the current unknowns, as otherwise the
//    // size of the vectors does not match!
//    shared_ptr<TimeSchemeGLM> mainScheme = dynamic_pointer_cast<TimeSchemeGLM>(
//        feFunctions_[MAG_POTENTIAL]->GetTimeScheme());
//    assert(mainScheme);
//
//    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(*mainScheme));
//
//    feFunctions_[COIL_CURRENT]->SetTimeScheme(myScheme);
//  }
//}
//
//
//// ========================================================================
////  Class for calculating coil results
//// ========================================================================
//
//template<class TYPE>
//CoilResultFunctor<TYPE>::
//CoilResultFunctor(shared_ptr<MathParser> mp,
//                  Grid * ptGrid,
//                  shared_ptr<BaseFeFunction> feFct,
//                  shared_ptr<ResultInfo> info,
//                  const std::map<Coil::IdType, shared_ptr<Coil> >& coilMap,
//                  PtrCoefFct potential,
//                  PtrCoefFct potentialD1,
//                  shared_ptr<BaseFeFunction> currentFct,
//                  shared_ptr<CoefFunctionMulti> conductivity )
//                  : ResultFunctor(info) {
//  mp_ = mp;
//  curTime_ = 0.0;
//  prevTime_ = 0.0;
//  feFct_ = feFct;
//  ptGrid_ = ptGrid;
//  coils_ = coilMap;
//  potential_ = potential;
//  potentialD1_ = potentialD1;
//  currentFunc_ = currentFct;
//  conduc_ = conductivity;
//  Global::ComplexPart type;
//  bool isComplex = std::tr1::is_same<TYPE,Complex>::value;
//  type = isComplex ? Global::COMPLEX : Global::REAL;
//
//  aTimesJ_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1,
//                                       isComplex ));
//  aD1TimesJ_.reset(new CoefFunctionMulti(CoefFunction::SCALAR, 1, 1,
//                                         isComplex ));
//  mphT_ = mp_->GetNewHandle();
//  mp_->SetExpr(mphT_, "t");
//
//  // Loop over all coils
//  std::map<Coil::IdType, shared_ptr<Coil> >::const_iterator coilIt;
//  coilIt = coils_.begin();
//  for( ; coilIt != coils_.end(); coilIt++ ) {
//
//    Coil::IdType id = coilIt->first;
//    shared_ptr<Coil> actCoil = coilIt->second;
//    const std::map<RegionIdType, shared_ptr<Coil::Part> >&
//    parts = actCoil->parts_;
//    std::map<RegionIdType, shared_ptr<Coil::Part> >::const_iterator partIt;
//    partIt = parts.begin();
//
//    // Loop over all part regions
//    for( ; partIt != parts.end(); partIt++ ) {
//      RegionIdType actRegion = partIt->first;
//      const Coil::Part & actPart = *(partIt->second);
//      PtrCoefFct jUnit = actPart.jUnitVec;
//      // Create CoefFunction (A * e_J)
//      CoefXprBinOp aExpr(mp_, potential_, jUnit, CoefXpr::OP_MULT);
//      aTimesJ_->AddRegion(actRegion,
//                          CoefFunction::Generate( mp, type, aExpr ));
//      if( potentialD1_) {
//        // Create CoefFunction (dA/dt * e_J)
//        CoefXprBinOp aDotExpr(mp_, potentialD1_, jUnit, CoefXpr::OP_MULT);
//        aD1TimesJ_->AddRegion(actRegion,
//                              CoefFunction::Generate( mp, type, aDotExpr ));
//      }
//    } // loop: parts
//
//    // initialize coil current with zero
//    prevCoilCurrent_[actCoil] = 0.0;
//  } // loop: coils
//
//
//  // Create result functors
//  fluxInt_.reset(new ResultFunctorIntegrate<TYPE>(aTimesJ_, feFct_,
//                                                  resultInfo_));
//
//  fluxD1Int_.reset(new ResultFunctorIntegrate<TYPE>(aD1TimesJ_,
//                                                    feFct_,
//                                                    resultInfo_));
//}
//
//template<class TYPE>
//CoilResultFunctor<TYPE>::~CoilResultFunctor() {
//}
//
//template<class TYPE>
//void CoilResultFunctor<TYPE>::EvalResult(shared_ptr<BaseResult> res ) {
//  assert( res->GetResultInfo()->resultType == resultInfo_->resultType);
//
//  Result<TYPE>& actSol = static_cast<Result<TYPE>& >(*res);
//  EntityIterator coilIt = actSol.GetEntityList()->GetIterator();
//  Vector<TYPE>& vec = actSol.GetVector();
//  vec.Resize( coilIt.GetSize() );
//
//  // Loop over coils
//  for( coilIt.Begin(); !coilIt.IsEnd(); coilIt++ ) {
//    TYPE coilRes = TYPE(0.0);
//
//    // switch depending on result type
//    switch( resultInfo_->resultType ) {
//      case COIL_INDUCTANCE:
//        coilRes = CalcSingleCoilInductance(coilIt);
//        break;
//
//      case COIL_INDUCED_VOLTAGE:
//        coilRes = CalcSingleCoilVoltage(coilIt);
//        break;
//
//      case COIL_CURRENT:
//        coilRes = CalcSingleCoilCurrent(coilIt);
//        break;
//
//      case COIL_LINKED_FLUX:
//        coilRes = CalcSingleCoilFlux(coilIt, false);
//        break;
//
//      default:
//        EXCEPTION("Result "
//            << SolutionTypeEnum.ToString(resultInfo_->resultType)
//            << " can not be computed for coils!");
//        break;
//    }
//
//    vec[coilIt.GetPos()]  = coilRes;
//  } // loop: coils
//}
//
//
//template<class TYPE>
//TYPE CoilResultFunctor<TYPE>::CalcSingleCoilCurrent( EntityIterator& coilIt) {
//  // If there is a new time step we move the stored coil currents to the
//  // previous step.
//  Double t = 0.0;
//  try {
//    t = mp_->Eval(mphT_);
//  } catch (...) {} // if this is a harmonic analysis, just use t=0
//  if (t > curTime_) {
//    prevCoilCurrent_ = coilCurrent_;
//    coilCurrent_.clear();
//    prevTime_ = curTime_;
//    curTime_ = t;
//  }
//
//  shared_ptr<Coil> actCoil = coilIt.GetCoil();
//  // switch depending on source type of coil
//  TYPE actCurrent = 0.0;
//  LocPointMapped lpm;
//  if( actCoil->sourceType_ == Coil::CURRENT  ) {
//    actCoil->srcVal_->GetScalar(actCurrent, lpm);
//  } else  if( actCoil->sourceType_ == Coil::VOLTAGE ) {
//    // get hold of voltage FeSpace and get voltage of given coil
//    Vector<TYPE> sol;
//    currentFunc_->GetEntitySolution(sol, coilIt);
//    actCurrent = sol[0];
//  }
//
//  coilCurrent_[actCoil] = actCurrent;
//  return actCurrent;
//}
//
//template<class TYPE>
//TYPE CoilResultFunctor<TYPE>::CalcSingleCoilVoltage( EntityIterator& coilIt ) {
// return CalcSingleCoilFlux(coilIt, true);
//}
//
//template<class TYPE>
//TYPE CoilResultFunctor<TYPE>::CalcSingleCoilInductance( EntityIterator& coilIt ) {
//  // The inductance of a coil is given as L = Psi / I
//  TYPE current = CalcSingleCoilCurrent(coilIt);
//
//  if (curTime_ > 0.0)  { // transient case
//    shared_ptr<Coil> actCoil = coilIt.GetCoil();
//    TYPE totalFluxD1 = CalcSingleCoilFlux(coilIt, true);
//    TYPE currentD1 = (current - prevCoilCurrent_[actCoil]) / (curTime_ - prevTime_);
//    if ( std::abs(currentD1) > EPS) {
//      return std::abs(totalFluxD1 / currentD1);
//    }
//    else {
//      return TYPE(0.0);
//    }
//  }
//  else { // static or harmonic case
//    TYPE totalFlux = CalcSingleCoilFlux(coilIt, false);
//    if( std::abs(current) > EPS ) {
//      return totalFlux / current;
//    } else {
//      return TYPE(0.0);
//    }
//  }
//}
//
//template<class TYPE>
//TYPE CoilResultFunctor<TYPE>::CalcSingleCoilFlux( EntityIterator& coilIt,
//                                                  bool timeDeriv ) {
//  shared_ptr<Coil> actCoil = coilIt.GetCoil();
//
//  // Assemble an entity list with all regions
//  shared_ptr<NameList> nameList(new NameList(ptGrid_));
//
//  // Loop over all parts
//  StdVector<std::string> regions;
//  std::map<RegionIdType, shared_ptr<Coil::Part> >::const_iterator partIt;
//  partIt = actCoil->parts_.begin();
//  for( ; partIt != actCoil->parts_.end(); partIt++ ) {
//    regions.Push_back(ptGrid_->GetRegion().ToString(partIt->first));
//  }
//  nameList->SetNames(regions);
//  shared_ptr<Result<TYPE> > tmp(new Result<TYPE>());
//  tmp->SetResultInfo(this->resultInfo_);
//  tmp->SetEntityList(nameList);
//  Vector<TYPE>& vec = tmp->GetVector();
//  vec.Resize( regions.GetSize() );
//  vec.Init();
//
//  TYPE sum = 0.0;
//  if( timeDeriv ) {
//    fluxD1Int_->EvalResult(tmp);
//    for( UInt i = 0; i < vec.GetSize(); ++i ) {
//      // U_ind = - dPsi/dt
//      sum -= vec[i];
//    }
//  } else {
//    fluxInt_->EvalResult(tmp);
//
//    for( UInt i = 0; i < vec.GetSize(); ++i ) {
//      sum += vec[i];
//    }
//  }
//  return sum;
//}
//
//// ========================================================================
////  CoefFunction for current density of voltage driven coils
//// ========================================================================
//
////! Coefficient function for current density of voltage driven coil
//template<class TYPE>
//CoefFunctionJVoltCoil<TYPE>::
//CoefFunctionJVoltCoil(Grid* ptGrid, shared_ptr<Coil> coil,
//                      shared_ptr<BaseFeFunction> currentFnc)
//{
//  // this type of coefficient is always constant
//  dependType_ = GENERAL;
//  isAnalytic_ = false;
//  isComplex_ = std::tr1::is_same<TYPE,Complex>::value;
//  dimType_ = VECTOR;
//  if( ptGrid->GetDim() == 2) {
//    dim_ = 1;
//  } else {
//    dim_ = 3;
//  }
//
//  coil_ = coil;
//  coilList_.reset(new CoilList(ptGrid));
//  coilList_->AddCoil(coil);
//  currentFunc_ = currentFnc;
//}
//
//template<class TYPE>
//CoefFunctionJVoltCoil<TYPE>::~CoefFunctionJVoltCoil() {}
//
//template<class TYPE>
//std::string  CoefFunctionJVoltCoil<TYPE>::
//ToString() const {
//  return "J-of-Coil "+this->coil_->coilId_;
//}
//
//template<class TYPE>
//void  CoefFunctionJVoltCoil<TYPE>::
//GetVector(Vector<TYPE>& vec, const LocPointMapped& lpm ) {
//  Vector<TYPE> temp;
//  EntityIterator it = coilList_->GetIterator();
//  it.Begin();
//
//  // obtain scaled unit direction
//  vec.Init();
//  RegionIdType actRegion = lpm.ptEl->regionId;
//  coil_->parts_[actRegion]->jUnitVec->GetVector(vec, lpm);
//
//  // multiply by absolute current
//  currentFunc_->GetEntitySolution(temp, it);
//  vec *= temp[0];
//}
//
//// ========================================================================
////  CoefFunction for magnetic field intensity
//// ========================================================================
//
//template<class TYPE>
//CoefFunctionMagFieldIntensity<TYPE>::
//CoefFunctionMagFieldIntensity(PtrCoefFct hField,
//                              PtrCoefFct magnetization,
//                              PtrCoefFct reluctivity)
//{
//  dimType_ = VECTOR;
//  dependType_ = SOLUTION;
//  isAnalytic_ = false;
//  isComplex_ = std::tr1::is_same<TYPE,Complex>::value;
//
//  hField_ = hField;
//  magnetization_ = magnetization;
//  reluc_ = reluctivity;
//}
//
//template<class TYPE>
//void CoefFunctionMagFieldIntensity<TYPE>::GetVector(Vector<TYPE>& vec,
//                                                    const LocPointMapped& lpm)
//{
//  hField_->GetVector(vec, lpm);
//
//  Vector<TYPE> remanence;
//  magnetization_->GetVector(remanence, lpm);
//
//  Matrix<TYPE> nu;
//  reluc_->GetTensor(nu, lpm);
//
//  vec -= nu * remanence;
//}
//
//template<class TYPE>
//UInt CoefFunctionMagFieldIntensity<TYPE>::GetVecSize() const {
//  return hField_->GetVecSize();
//}
//
//
//// ========================================================================
////  Maxwell stress tensor
//// ========================================================================
//template<class FE>
//CoefFunctionMaxwell<FE>::CoefFunctionMaxwell(PtrCoefFct Bfield,
//    std::map<SolutionType, shared_ptr<CoefFunctionMulti> > matCoefs,
//    Grid* ptGrid, MagBasePDE::MagFormulationType formulation,
//    shared_ptr<FeSpace> feSpace)
//{
//  bField_ = Bfield;
//  matCoef_ = matCoefs;
//  ptGrid_ = ptGrid;
//  feSpace_ = feSpace;
//
//  dimType_ = VECTOR;
//  dependType_ = SOLUTION;
//}
//
//template<class FE>
//void CoefFunctionMaxwell<FE>::GetVector(Vector<Double>& coefVec,
//                                    const LocPointMapped& lpm)
//{
//  assert(this->dimType_ == VECTOR);
//
//  Double permeability = 0.0;
//
//  matCoef_[MAG_ELEM_PERMEABILITY]->GetScalar(permeability, lpm);
//
//  //UInt dim = bField_->GetVecSize();
//  UInt dim = ptGrid_->GetDim();
//  UInt dimBmat;
//  UInt numNodes = Elem::shapes[lpm.shapeMap->GetElem()->type].numNodes;
//  Vector<Double> Bvec; //(numNodes * dim), B2(numNodes);
//  Double B2;
//  Vector<Double> maxwell;
//  if (dim == 3) {
//    dimBmat = 6;
//  }
//  /*else if (ptGrid_->IsAxi()) {
//    dimBmat = 4;
//  }*/
//  else {
//    dimBmat = 3;
//  }
//  maxwell.Resize(dimBmat * numNodes);
//
//  // Get magnetic flux density in all nodes of this element
//  // WARNING: This will only work for elements with local support
//  Vector<Double> globCoord;
//  //Vector<Double> coefAtNode;
//  LocPoint lpNode;
//  LocPointMapped lpmNode;
//  for (UInt node = 0; node < numNodes; ++node) {
//    lpNode.coord = Elem::shapes[lpm.ptEl->type].nodeCoords[node];
//    lpmNode.Set(lpNode, lpm.shapeMap, 1.0);
//    bField_->GetVector(Bvec, lpmNode);
//    /*for (UInt dof = 0; dof < dim; ++dof) {
//      Bvec[node*dim+dof] = coefAtNode[dof];
//    }
//    B2[node] = coefAtNode.Inner();*/
//    B2 = Bvec.Inner();
//
//    if (dim == 3) {
//      maxwell[6*node+0] = Bvec[0]*Bvec[0] - 0.5*B2;
//      maxwell[6*node+1] = Bvec[1]*Bvec[1] - 0.5*B2;
//      maxwell[6*node+2] = Bvec[2]*Bvec[2] - 0.5*B2;
//      maxwell[6*node+3] = Bvec[1]*Bvec[2];
//      maxwell[6*node+4] = Bvec[0]*Bvec[2];
//      maxwell[6*node+5] = Bvec[0]*Bvec[1];
//    }
//    /*else if (ptGrid_->IsAxi()) {
//      maxwell[4*node+0] = Bvec[0]*Bvec[0] - 0.5*B2;
//      maxwell[4*node+1] = Bvec[1]*Bvec[1] - 0.5*B2;
//      maxwell[4*node+2] = Bvec[0]*Bvec[1];
//      maxwell[4*node+3] = Bvec[2]*Bvec[2] - 0.5*B2;
//    }*/
//    else{
//      maxwell[3*node+0] = Bvec[0]*Bvec[0] - 0.5*B2;
//      maxwell[3*node+1] = Bvec[1]*Bvec[1] - 0.5*B2;
//      maxwell[3*node+2] = Bvec[0]*Bvec[1];
//    }
//  }
//
//  maxwell *= -1.0 / permeability;
//
//  FE *fe = static_cast<FE*>(feSpace_->GetFe(lpm.ptEl->elemNum));
//
//  Matrix<Double> xiDx, bMat;
//  fe->GetGlobDerivShFnc( xiDx, lpm, lpm.shapeMap->GetElem() );
//
//  const UInt numFncs = xiDx.GetNumRows();
//
//  // Set correct size of matrix B and initialise with zeros
//  bMat.Resize(dim, numFncs * dimBmat);
//  bMat.Init();
//
//  UInt iFunc, d;
//  for (iFunc = 0; iFunc < numFncs; ++iFunc) {
//    for (d = 0; d < dim; ++d) {
//      bMat[d][dimBmat*iFunc+d] = xiDx[iFunc][d];
//    }
//  }
//
//  if (dim == 3) {
//    for (iFunc = 0; iFunc < numFncs; ++iFunc) {
//      bMat[1][dimBmat*iFunc+dim  ] = xiDx[iFunc][2];
//      bMat[2][dimBmat*iFunc+dim  ] = xiDx[iFunc][1];
//
//      bMat[0][dimBmat*iFunc+dim+1] = xiDx[iFunc][2];
//      bMat[2][dimBmat*iFunc+dim+1] = xiDx[iFunc][0];
//
//      bMat[0][dimBmat*iFunc+dim+2] = xiDx[iFunc][1];
//      bMat[1][dimBmat*iFunc+dim+2] = xiDx[iFunc][0];
//    }
//  }
//  else {
//    for (iFunc = 0; iFunc < numFncs; ++iFunc) {
//      bMat[0][dimBmat*iFunc+dim] = xiDx[iFunc][1];
//      bMat[1][dimBmat*iFunc+dim] = xiDx[iFunc][0];
//    }
//    /*if (ptGrid_->IsAxi()) {
//      Matrix<Double> shape;
//      fe->GetShFnc( shape, lpm.lp, lpm.shapeMap->GetElem() );
//      const Double oneOverR = 1.0 / lpm.globPoint[0];
//    }*/
//  }
//
//  coefVec = bMat * maxwell;
//}
//
//template<class FE>
//void CoefFunctionMaxwell<FE>::GetVector(Vector<Complex>& coefMat,
//                                    const LocPointMapped& lpm ) {
//  assert(this->dimType_ == VECTOR);
//
//  EXCEPTION("CoefFunctionMaxwell for Harmonic Analysis not implemented");
//
//}
//
//
//// explicit template instantiation
//template class CoefFunctionMagFieldIntensity<Double>;
//template class CoefFunctionMagFieldIntensity<Complex>;
//template class CoefFunctionMaxwell<FeH1>;
//template class CoefFunctionMaxwell<FeHCurl>;
//
//} // end of namespace
