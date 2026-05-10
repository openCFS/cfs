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
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefFunctionConst.hh"

#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/LinForms/BDUInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Forms/BiLinForms/BiLinWrappedLinForm.hh"
#include "Driver/MultiHarmonicDriver.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"

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

    magnetizationSet_ = false;
    coilOptimization_ = false;
	
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
    // > moved to InitMagnetization
    // > Reason: during constructor, isComplex_ is not set yet
//    polarization_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
//    magnetization_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
//    fieldIntensity_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));

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

    InitMagnetization();
  }

  void MagBasePDE::InitTimeStepping()
  {
    // Check if time integration is defined in XML input
    PtrParamNode transientNode = myParam_->GetParent()->GetParent()->Get("analysis")->Get("transient", ParamNode::PASS);
    PtrParamNode integrationScheme = transientNode->Get("integrationScheme", ParamNode::PASS);

    // Use complete implicit scheme by default
    GLMScheme* scheme = nullptr;
    if (integrationScheme)
    {
      scheme = GetXmlDefinedScheme(integrationScheme);
    }
    else
    {
      scheme = new Trapezoidal(1.0);
    }

    TimeSchemeGLM::NonLinType nlType = (nonLin_ || isHysteresis_) ? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType));
    feFunctions_[MAG_POTENTIAL]->SetTimeScheme(myScheme);

    // Important: Create a new time scheme for each additional feFunction
    // NEW: from NACS - copy stepping scheme from mag potential
    shared_ptr<TimeSchemeGLM> mainScheme = dynamic_pointer_cast<TimeSchemeGLM>(
        feFunctions_[MAG_POTENTIAL]->GetTimeScheme());
    assert(mainScheme);

    if (hasVoltCoils_)
    {
      shared_ptr<BaseTimeScheme> myScheme2(new TimeSchemeGLM(*mainScheme));
      feFunctions_[COIL_CURRENT]->SetTimeScheme(myScheme2);
    }

    if (isMixed_)
    {
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
    /*
     * Moved reset of CoefFunctionMulti here as isComplex_ is not known in constructor; thus
     * the inserted CoefFunctions do not fit; this will throw EXCEPTION( "All coefficient functions must have the same complexType");
     */
    polarization_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
    magnetization_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
    fieldIntensity_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
    
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
    fluxDens->SetFeFunction(feFunctions_[MAG_POTENTIAL]);
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
   * Inits coefficient functions for H (B/mu_0 - M), M and P (mu_0*M)
   * > M can be hysteretic (using e.g. Preisach), constant (defined in mat.xml) or 0 (default)
   * > H and P will be defined based on M
   */
  void MagBasePDE::InitMagnetization() {
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
        actSDMat->GetString(hystType, HYST_MODEL);

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
                  magFieldCoef, tensorType, MAG_RELUCTIVITY_TENSOR, mySpace ));
//          std::cout << "hysteresisCoefs_->AddRegion( actRegion, hystPol);" << std::endl;
          hysteresisCoefs_->AddRegion( actRegion, hystPol);

          PtrCoefFct hystOutput = hystPol->GenerateOutputCoefFnc("MagPolarization");
          polarization_->AddRegion( actRegion, hystOutput);

          PtrCoefFct hystOutput2 = hystPol->GenerateOutputCoefFnc("MagMagnetization");
          magnetization_->AddRegion( actRegion, hystOutput2);

          PtrCoefFct hystOutput3 = hystPol->GenerateOutputCoefFnc("MagFieldIntensityHyst");
          fieldIntensity_->AddRegion( actRegion, hystOutput3);
//          std::string warnmsg = "Hysteretic M set on region " + regionName;
//          std::cout << warnmsg << std::endl;
        }
      } 
      else {
        // check material for fixed magnetization
        int hasFixedMagnetization = 0;
        materials_[actRegion]->GetScalar(hasFixedMagnetization,PRESCRIBED_MAGNETIZATION);
        CoefFunctionConst<Double>* magFnc = new CoefFunctionConst<Double>();
        CoefFunctionConst<Double>* polFnc = new CoefFunctionConst<Double>();
        
        CoefFunctionConst<Complex>* magFncC = new CoefFunctionConst<Complex>();
        CoefFunctionConst<Complex>* polFncC = new CoefFunctionConst<Complex>();
        
        Vector<Double> magVec = Vector<Double>(dim_);
        magVec.Init();
        Vector<Double> polVec = Vector<Double>(dim_);
        polVec.Init();
        
        Vector<Complex> magVecC = Vector<Complex>(dim_);
        magVecC.Init();
        Vector<Complex> polVecC = Vector<Complex>(dim_);
        polVecC.Init();
        
        if(hasFixedMagnetization == 1){
//          std::string warnmsg = "Fixed M set on region " + regionName;
//          std::cout << warnmsg << std::endl;
          Double magX,magY,magZ;
          
          materials_[actRegion]->GetScalar(magX,PRESCRIBED_MAGNETIZATION_X,Global::REAL);
          materials_[actRegion]->GetScalar(magY,PRESCRIBED_MAGNETIZATION_Y,Global::REAL);
          materials_[actRegion]->GetScalar(magZ,PRESCRIBED_MAGNETIZATION_Z,Global::REAL);
          
          magVec[0] = magX;
          magVec[1] = magY;
          if(dim_ == 3){
            magVec[2] = magZ;
          }
          
          magVecC[0] = Complex(magX);
          magVecC[1] = Complex(magY);
          if(dim_ == 3){
            magVecC[2] = Complex(magZ);
          }

          /*
           * according to definition P = mu_0 * M
           * however, one could also think about using mu as given in material file
           * TODO: decide and implement!
           */
          Double mu_0 = 4*M_PI*1e-7;
          polVec.Add(mu_0,magVec);
          polVecC.Add(Complex(mu_0),magVecC);
        } 
//        else {
////          std::string warnmsg = "No M set on region " + regionName;
////          std::cout << warnmsg << std::endl;
//          // just set 0-vector
//        }
        magFnc->SetVector(magVec);
        polFnc->SetVector(polVec);
        
        magFncC->SetVector(magVecC);
        polFncC->SetVector(polVecC);
        
        shared_ptr<CoefFunctionConst<Double> > magFncSharedPtr(magFnc);
        shared_ptr<CoefFunctionConst<Double> > polFncSharedPtr(polFnc);
        
        shared_ptr<CoefFunctionConst<Complex> > magFncSharedPtrC(magFncC);
        shared_ptr<CoefFunctionConst<Complex> > polFncSharedPtrC(polFncC);
         
        if(hasFixedMagnetization == 1){
          // mRHSRegions_ is for excitation on rhs
          if(isComplex_){
            mRHSRegions_[actRegion] = magFncSharedPtrC;
          } else {
            mRHSRegions_[actRegion] = magFncSharedPtr;
                      
//            std::cout << "Prescribed magnetization vector: " << magFncSharedPtr->ToString() << std::endl;
//          
          }
        }
        
//        std::cout << "Part 0: setup of rhs magnetization" << std::endl;
//        
        // magnetization_ and polarization are for output
        if(isComplex_){
          magnetization_->AddRegion( actRegion, magFncSharedPtrC);
          polarization_->AddRegion( actRegion, polFncSharedPtrC);
        } else {
          magnetization_->AddRegion( actRegion, magFncSharedPtr);
          polarization_->AddRegion( actRegion, polFncSharedPtr);
        }

//        std::cout << "Part 1: setup of magnetization and polarization done" << std::endl;
//        
        /*
         * FIELD INTENSITY
         */
        // define field intensity on non-hysteretic region, too!
        shared_ptr<BaseFeFunction> feFct = feFunctions_[MAG_POTENTIAL];
        
        StdVector<std::string> vecComponents;
        if( dim_ == 3 ) {
          vecComponents = "x", "y", "z";
        }
        else if( isaxi_ ) {
          vecComponents = "r", "z";
        }
        else {
          vecComponents = "x", "y";
        }
        
        shared_ptr<ResultInfo> magIntens ( new ResultInfo );
        magIntens->resultType = MAG_FIELD_INTENSITY;
        magIntens->dofNames = vecComponents;
        //magIntens->SetVectorDOFs(dim_, isaxi_);
        magIntens->unit = "A/m";
        magIntens->definedOn = ResultInfo::ELEMENT;
        magIntens->entryType = ResultInfo::VECTOR;
        
        shared_ptr<CoefFunctionFormBased> magIntensFunc;
        if( isComplex_ ) {
          magIntensFunc.reset(new CoefFunctionFlux<Complex>(feFct, magIntens));
        } else {
          magIntensFunc.reset(new CoefFunctionFlux<Double>(feFct, magIntens));
        }
//        DefineFieldResult( magIntensFunc, magIntens );
        
        // subtract fixed magnetization for output!
        if(hasFixedMagnetization == 1){
          Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
          CoefXprBinOp sub(mp_, magIntensFunc, mRHSRegions_[actRegion], CoefXpr::OP_SUB);
          PtrCoefFct magIntensFuncCorrected(CoefFunction::Generate(mp_, part, sub));
          fieldIntensity_->AddRegion( actRegion, magIntensFuncCorrected);
//          stiffFormCoefs_.insert(magIntensFuncCorrected);
        } else {
          fieldIntensity_->AddRegion( actRegion, magIntensFunc);
//          stiffFormCoefs_.insert(magIntensFunc);
        }
        stiffFormCoefs_.insert(magIntensFunc);
//        std::cout << "Part 2: setup of field intensity done" << std::endl;
      }
    }
//    std::cout << "Magnetization Initiatlized!" << std::endl;
    magnetizationSet_ = true;
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
          LinearForm* curInt = NULL;

          /*
            generate source current vector
            for the non-multiharmonic case it is a simple PtrCoefFct type
            but due to the flexibility to use multiharmonic excitation,
            we need a map of PtrCoefFct's, where the key is the harmonic number.
            So in the non-multiharmonic case, we simply have one key "0"
          */
          CoefFunctionOpt* cfoc = NULL; // we might do optimization and then we have such a thing
          std::map<UInt, PtrCoefFct> jFct;
          if( actCoil.sourceType_ == Coil::CURRENT ){
            CoefXprVecScalOp iVec = CoefXprVecScalOp(mp_, actPart.jUnitVec, actCoil.srcVal_, CoefXpr::OP_MULT);
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
            for(auto& h : jFct){
              coilCurrentDens_[actRegion] = h.second;
              curInt = GetCurrentDensityInt( scaling, h.second );
              curInt->SetName("CoilIntegrator");
              curInt->SetHarm(h.first);
              LinearFormContext * coilContext = new LinearFormContext( curInt );
              coilContext->SetEntities( actSDList );
              coilContext->SetFeFunction( feFunc );
              assemble_->AddLinearForm( coilContext );
            }
          }else{
            // Classic Case

            if(actCoil.coilOptimization_ == true)
            {
              coilOptimization_ = true;
              if(domain->HasDesign())
              {
                cfoc = new CoefFunctionOpt(domain->GetDesign(), jFct[0], NO_MATERIAL, this); // no explicit material
                jFct[0].reset(cfoc);
              }
            }

            coilCurrentDens_[actRegion] = jFct[0];
            curInt = GetCurrentDensityInt( scaling, jFct[0] );
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

  void MagBasePDE::GenerateLorentzForceResults(CoupledField::StdVector<std::string> &vecComponents, shared_ptr<CoupledField::CoefFunctionMulti> &tcdCoef,
    CoupledField::PtrCoefFct &bFunc, CoupledField::Global::ComplexPart &part, shared_ptr<CoupledField::BaseFeFunction> &feFct) {

    if (analysistype_ == HARMONIC || analysistype_ == MULTIHARMONIC)
    {
      CoefXpr::OpType op_cr = isaxi_ ? CoefXpr::OP_CROSS_AXI : CoefXpr::OP_CROSS;

      // === LORENTZ FORCE DENSITY HARMONIC CONSTANT ===
      shared_ptr<ResultInfo> lfd_static(new ResultInfo);
      lfd_static->resultType = MAG_FORCE_LORENTZ_DENSITY_STATIC;
      lfd_static->dofNames = vecComponents;
      lfd_static->unit = "N/m^3";
      lfd_static->definedOn = ResultInfo::ELEMENT;
      lfd_static->entryType = ResultInfo::VECTOR;

      // assemble coefficient function F_L = 0.5 * (Re(J) X Re(B) + Im(J) X Im(B))  
      PtrCoefFct tcdRe = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprUnaryOp( mp_, tcdCoef, CoefXpr::OP_RE ) );
      PtrCoefFct bFuncRe = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprUnaryOp( mp_, bFunc, CoefXpr::OP_RE ) );
      PtrCoefFct tcdIm = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprUnaryOp( mp_, tcdCoef, CoefXpr::OP_IM ) );
      PtrCoefFct bFuncIm = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprUnaryOp( mp_, bFunc, CoefXpr::OP_IM ) );
      
      PtrCoefFct realCrossProd = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp( mp_, tcdRe, bFuncRe, op_cr ) );
      PtrCoefFct imagCrossProd = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp( mp_, tcdIm, bFuncIm, op_cr ) );
      
      PtrCoefFct sumCrossProds = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp( mp_, realCrossProd, imagCrossProd, CoefXpr::OP_ADD ) );
      PtrCoefFct lfdFunc_static = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp(mp_, "0.5", sumCrossProds, CoefXpr::OP_MULT ) ); 

      // scale solution for MULTIHARMONIC case
      // consider that the quantities in multiharmonic are 1/2 of the harmonic ones,
      // the total force is obtained by summing all the harmonics
      if (analysistype_ == MULTIHARMONIC) 
      {
        lfdFunc_static = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp(mp_, "2", lfdFunc_static, CoefXpr::OP_MULT ));
      }
      
      DefineFieldResult( lfdFunc_static, lfd_static);

      // === LORENTZ FORCE DENSITY HARMONIC TIME ===
      shared_ptr<ResultInfo> lfd_harmonic(new ResultInfo);
      lfd_harmonic->resultType = MAG_FORCE_LORENTZ_DENSITY_HARMONIC;
      lfd_harmonic->dofNames = vecComponents;
      lfd_harmonic->unit = "N/m^3";
      lfd_harmonic->definedOn = ResultInfo::ELEMENT;
      lfd_harmonic->entryType = ResultInfo::VECTOR;

      // assemble coefficient function F_L = 0.5 * (Re(J) X Re(B) - Im(J) X Im(B)) +
      //                                     0.5j * (Im(J) X Re(B) + Re(J) X Im(B))
      
      PtrCoefFct ImJReBCrossProd = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp( mp_, tcdIm, bFuncRe, op_cr ) );
      PtrCoefFct ReJImBCrossProd = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp( mp_, tcdRe, bFuncIm, op_cr ) );
      
      // 0.5 * (Re(J) X Re(B) - Im(J) X Im(B))
      PtrCoefFct subReAmp = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp ( mp_, realCrossProd, imagCrossProd, CoefXpr::OP_SUB ) );
      subReAmp = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp(mp_, "0.5", subReAmp, CoefXpr::OP_MULT ));
      
      // 0.5j * (Im(J) X Re(B) + Re(J) X Im(B))
      PtrCoefFct halfJ = CoefFunction::Generate( mp_, Global::COMPLEX, "0", "0.5");
      PtrCoefFct sumImAmp = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp( mp_, ImJReBCrossProd, ReJImBCrossProd, CoefXpr::OP_ADD ) );
      sumImAmp = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp(mp_, halfJ, sumImAmp, CoefXpr::OP_MULT ));

      PtrCoefFct lfdFunc_harmonic = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp( mp_, subReAmp, sumImAmp, CoefXpr::OP_ADD));

      // scale solution for MULTIHARMONIC case
      // consider that the quantities in multiharmonic are 1/2 of the harmonic ones,
      // the total force is obtained by summing all the harmonics
      if (analysistype_ == MULTIHARMONIC) 
      {
        lfdFunc_harmonic = CoefFunction::Generate( mp_, Global::COMPLEX, CoefXprBinOp(mp_, "2", lfdFunc_harmonic, CoefXpr::OP_MULT ));
      }
      
      DefineFieldResult( lfdFunc_harmonic, lfd_harmonic);
    }
    else {
      // === LORENTZ FORCE DENSITY - STATIC AND TRANSIENT ===
      shared_ptr<ResultInfo> lfd(new ResultInfo);
      lfd->resultType = MAG_FORCE_LORENTZ_DENSITY;
      lfd->dofNames = vecComponents;
      lfd->unit = "N/m^3";
      lfd->definedOn = ResultInfo::ELEMENT;
      lfd->entryType = ResultInfo::VECTOR;
      availResults_.insert( lfd );

      // assemble coefficient function F_L = J X B

      // switch type of cross-product depending on dimensionality
      CoefXpr::OpType op = isaxi_ ? CoefXpr::OP_CROSS_AXI : CoefXpr::OP_CROSS;
      PtrCoefFct lfdFunc = CoefFunction::Generate( mp_, part,
              CoefXprBinOp(mp_,  tcdCoef, bFunc, op) );
      DefineFieldResult( lfdFunc, lfd);
    
      // === LORENTZ FORCE (TOTAL) ===
      shared_ptr<ResultInfo> lf(new ResultInfo);
      lf->resultType = MAG_FORCE_LORENTZ;
      lf->dofNames = vecComponents;
      lf->unit = "N";
      lf->definedOn = ResultInfo::REGION;
      lf->entryType = ResultInfo::VECTOR;
      availResults_.insert( lf );

      // build result functor for integration
      shared_ptr<ResultFunctor> lfFunc;
      if( isComplex_ ) {
        lfFunc.reset(new ResultFunctorIntegrate<Complex>(lfdFunc, feFct, lf ) );
      } else {
        lfFunc.reset(new ResultFunctorIntegrate<Double>(lfdFunc, feFct, lf ) );
      }
      resultFunctors_[MAG_FORCE_LORENTZ] = lfFunc;
    }
    
  }

  void MagBasePDE::GenerateMaxwellForce(CoupledField::StdVector<std::string> &vecComponents,
    CoupledField::PtrCoefFct &bFunc, shared_ptr<CoupledField::BaseFeFunction> &feFct) {

    if (analysistype_ == HARMONIC || analysistype_ == MULTIHARMONIC)
    {
      // Note: The positive normal direction is defined as the inward facing one (factor = -1.0).
      // For MULTIHARMONIC, field amplitudes are 1/2 of the harmonic ones;
      // the Maxwell stress is quadratic in B, so we apply the same factor-of-2 scaling
      // convention used for the Lorentz force to recover the correct total force
      // by summing contributions from all harmonics.
      if (analysistype_ == MULTIHARMONIC) {
        WARN("For multiharmonic analysis, Maxwell force density is not properly validated.");
      }
      Double normalFactor = (analysistype_ == MULTIHARMONIC) ? -2.0 : -1.0;

      // === MAXWELL FORCE DENSITY HARMONIC STATIC (DC, time-averaged) ===
      // T_DC = (1/2μ) * {(B_re·n)B_re + (B_im·n)B_im - (1/2)|B̂|²n}
      shared_ptr<ResultInfo> mfd_static(new ResultInfo);
      mfd_static->resultType = MAG_FORCE_MAXWELL_DENSITY_STATIC;
      mfd_static->dofNames = vecComponents;
      mfd_static->unit = "N/m^2";
      mfd_static->definedOn = ResultInfo::SURF_ELEM;
      mfd_static->entryType = ResultInfo::VECTOR;
      availResults_.insert( mfd_static );

      shared_ptr<CoefFunctionSurfMaxwell> maxForceDens_static( new CoefFunctionSurfMaxwell(false, matCoefs_, ptGrid_, normalFactor, mfd_static) );
      DefineFieldResult( maxForceDens_static, mfd_static);
      surfCoefFcts_[maxForceDens_static] = bFunc;

      // === MAXWELL FORCE DENSITY HARMONIC TIME (2ω complex amplitude) ===
      // T_2ω = (1/2μ) * {(B̂·n)B̂ - (1/2)(B̂·B̂)n}  (non-Hermitian products)
      shared_ptr<ResultInfo> mfd_harmonic(new ResultInfo);
      mfd_harmonic->resultType = MAG_FORCE_MAXWELL_DENSITY_HARMONIC;
      mfd_harmonic->dofNames = vecComponents;
      mfd_harmonic->unit = "N/m^2";
      mfd_harmonic->definedOn = ResultInfo::SURF_ELEM;
      mfd_harmonic->entryType = ResultInfo::VECTOR;
      availResults_.insert( mfd_harmonic );

      shared_ptr<CoefFunctionSurfMaxwell> maxForceDens_harmonic(
        new CoefFunctionSurfMaxwell(false, matCoefs_, ptGrid_, normalFactor, mfd_harmonic));
      DefineFieldResult( maxForceDens_harmonic, mfd_harmonic);
      surfCoefFcts_[maxForceDens_harmonic] = bFunc;
      WARN( "WARNING: For multiharmonic analysis, Maxwell force density is not properly validated." ); 
    
      // === MAXWELL FORCE STATIC (TOTAL) ===
      shared_ptr<ResultInfo> mf_static(new ResultInfo);
      mf_static->resultType = MAG_FORCE_MAXWELL_STATIC;
      mf_static->dofNames = vecComponents;
      mf_static->unit = "N";
      mf_static->definedOn = ResultInfo::SURF_REGION;
      mf_static->entryType = ResultInfo::VECTOR;
      availResults_.insert( mf_static );

      shared_ptr<ResultFunctor> mfFunc_static(
        new ResultFunctorIntegrate<Complex>(maxForceDens_static, feFct, mf_static));
      resultFunctors_[MAG_FORCE_MAXWELL_STATIC] = mfFunc_static;

      // === MAXWELL FORCE HARMONIC (TOTAL) ===
      shared_ptr<ResultInfo> mf_harmonic(new ResultInfo);
      mf_harmonic->resultType = MAG_FORCE_MAXWELL_HARMONIC;
      mf_harmonic->dofNames = vecComponents;
      mf_harmonic->unit = "N";
      mf_harmonic->definedOn = ResultInfo::SURF_REGION;
      mf_harmonic->entryType = ResultInfo::VECTOR;
      availResults_.insert( mf_harmonic );

      shared_ptr<ResultFunctor> mfFunc_harmonic(
        new ResultFunctorIntegrate<Complex>(maxForceDens_harmonic, feFct, mf_harmonic));
      resultFunctors_[MAG_FORCE_MAXWELL_HARMONIC] = mfFunc_harmonic;

    }
    else {
      // === MAXWELL FORCE DENSITY (STATIC / TRANSIENT) ===
      shared_ptr<ResultInfo> mfd(new ResultInfo);
      mfd->resultType = MAG_FORCE_MAXWELL_DENSITY;
      mfd->dofNames = vecComponents;
      mfd->unit = "N/m^3";
      mfd->definedOn = ResultInfo::SURF_ELEM;
      mfd->entryType = ResultInfo::VECTOR;
      availResults_.insert( mfd );
      // Note: The positive normal direction in this case is defined as the
      //       inward facing one.
      shared_ptr<CoefFunctionSurfMaxwell> maxForceDens(new CoefFunctionSurfMaxwell(false, matCoefs_, ptGrid_, -1.0, mfd));
      DefineFieldResult( maxForceDens, mfd);
      surfCoefFcts_[maxForceDens] = bFunc;

      // === MAXWELL FORCE (TOTAL) ===
      shared_ptr<ResultInfo> mf(new ResultInfo);
      mf->resultType = MAG_FORCE_MAXWELL;
      mf->dofNames = vecComponents;
      mf->unit = "N";
      mf->definedOn = ResultInfo::SURF_REGION;
      mf->entryType = ResultInfo::VECTOR;
      availResults_.insert( mf );

      // build result functor for integration
      shared_ptr<ResultFunctor> mfFunc;
      if( isComplex_ ) {
        mfFunc.reset(new ResultFunctorIntegrate<Complex>(maxForceDens, feFct, mf ) );
      } else {
        mfFunc.reset(new ResultFunctorIntegrate<Double>(maxForceDens, feFct, mf ) );
      }
      resultFunctors_[MAG_FORCE_MAXWELL] = mfFunc;
    }
  }


  void MagBasePDE::GenerateVWPForce(CoupledField::StdVector<std::string> &vecComponents,
    CoupledField::PtrCoefFct &bFunc, shared_ptr<CoupledField::BaseFeFunction> &feFct) {
    if( (analysistype_ != HARMONIC) && (analysistype_ != MULTIHARMONIC) ) {
      // === VIRTUAL WORK PRINCIPLE FORCE (TOTAL) ===
      shared_ptr<ResultInfo> vwp(new ResultInfo);
      vwp->resultType = MAG_FORCE_VWP;
      vwp->dofNames = vecComponents;
      vwp->unit = "N";
      vwp->definedOn = ResultInfo::SURF_REGION;
      vwp->entryType = ResultInfo::VECTOR;
      availResults_.insert( vwp );
      
      // define and save coefFunction
      shared_ptr<CoefFunctionSurfVWP> vwpForce(new CoefFunctionSurfVWP(false, matCoefs_, 1.0, vwp));
      surfCoefFcts_[vwpForce] = bFunc;
      
      // build result functor for integration
      shared_ptr<ResultFunctor> vwpFunc;
      if( isComplex_ ) {
        vwpFunc.reset(new ResultFunctorVWP<Complex>(vwpForce, feFct, vwp, ptGrid_ ) );
      } else {
        vwpFunc.reset(new ResultFunctorVWP<Double>(vwpForce, feFct, vwp, ptGrid_ ) );
      }
      resultFunctors_[MAG_FORCE_VWP] = vwpFunc;
    }
  }

  void MagBasePDE::ReadRegionVelocityField(std::string velocityId, shared_ptr<ElemList> actSDList, RegionIdType actRegion, bool& coefUpdateGeo){
    // Get result info object for flow
    shared_ptr<ResultInfo> velInfo = GetResultInfo(MEAN_FLUIDMECH_VELOCITY);

    // Add the region information
    PtrParamNode velNode = myParam_->Get("velocityList")->GetByVal("velocity","name",velocityId.c_str());

    // Read velocity coefficient function for this region and add it to velocity functor
    PtrCoefFct regionMoving;
    std::set<UInt> definedDofs;

    //we assume that velocity is real
    ReadUserFieldValues( actSDList, velNode, velInfo->dofNames, velInfo->entryType, isComplex_, regionMoving, definedDofs, coefUpdateGeo );
    VelocityCoef_->AddRegion( actRegion, regionMoving );
  }

} // end of namespace

