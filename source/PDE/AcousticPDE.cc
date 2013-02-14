#include "AcousticPDE.hh"

#include "General/defs.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"


//new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/BDUInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"
#include "Forms/Operators/ConvectiveOperator.hh"
#include "Forms/Operators/ConvectivePierceOperator.hh"

#include "FeBasis/FeFunctions.hh"
#include "Utils/StdVector.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"

#include "Domain/Results/ResultFunctor.hh"

#include "Driver/Assemble.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionPML.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"

#include <boost/lexical_cast.hpp>
#include <cmath>

#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Materials/AcousticMaterial.hh"

namespace CoupledField{

  DECLARE_LOG(acousticpde)
   DEFINE_LOG(acousticpde, "pde.acoustic")


  AcousticPDE::AcousticPDE( Grid* aGrid, PtrParamNode paramNode )
              : SinglePDE( aGrid, paramNode ){

    pdename_           = "acoustic";
    pdematerialclass_  = FLUID;
    maxTimeDerivOrder_ = 2;
    nonLin_            = false;

    isMechCoupled_ = false;

    //check for pressure or potential formulation
    std::string pdeFormulation = myParam_->Get("formulation")->As<std::string>();
    if(pdeFormulation == "default"){
      formulation_ = ACOU_PRESSURE;
    }else{
      formulation_ = SolutionTypeEnum.Parse(pdeFormulation);
    }
  }

  std::map<SolutionType, shared_ptr<FeSpace> >
  AcousticPDE::CreateFeSpaces( const std::string&  formulation,
                  PtrParamNode infoNode ){

    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H1"){
      std::string form = SolutionTypeEnum.ToString(formulation_);
      PtrParamNode potSpaceNode = infoNode->Get(form);
      crSpaces[formulation_] =
        FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[formulation_]->Init(solStrat_);
    }else{
      EXCEPTION("The formulation " << formulation << "of acoustic PDE is not known!");
    }
    return crSpaces;
  }

  void AcousticPDE::DefineIntegrators(){

    RegionIdType actRegion;
    // BaseMaterial * actSDMat = NULL;

    //type of geometry
    std::string geometryType;
    param->Get("domain")->GetValue("geometryType", geometryType );

    // convert to tensor type
    // SubTensorType tensorType = FULL;
    if (geometryType == "plane") {
      // tensorType = PLANE_STRAIN;
    } else if (geometryType == "axi") {
      // tensorType = AXI;
      isaxi_ = true;
    }

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    shared_ptr<FeSpace> mySpace = feFunctions_[formulation_]->GetFeSpace();

    //flag indicating frequency PML formulation
    bool harmonicPML = false;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      // Set current region and material
      actRegion = it->first;
      // actSDMat = it->second;

      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );

      // --- Set the FE ansatz for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId,integId);
      
      //=======================================================================
      // Generate coefficient functions
      //=======================================================================

      PtrCoefFct dens = materials_[actRegion]->GetScalCoefFnc( DENSITY, Global::REAL );
      PtrCoefFct blk = materials_[actRegion]->GetScalCoefFnc( ACOU_BULK_MODULUS, Global::REAL );
      
      // c0 = sqrt(bulk_modulus / density)
      PtrCoefFct c0 = 
          CoefFunction::Generate( Global::REAL,
                                  CoefXprUnaryOp( CoefXprBinOp(blk, dens, CoefXpr::OP_DIV),
                                  CoefXpr::OP_SQRT) );
      
      // store coefficient functions
      matCoefs_[ELEM_DENSITY]->AddRegion(actRegion, dens);
      matCoefs_[ACOU_ELEM_SPEED_OF_SOUND]->AddRegion( actRegion, c0);

      // if pde couples with mechanic, we have to multiply the density by -1
      PtrCoefFct factor;
      if ( isMechCoupled_ == true && formulation_ != ACOU_PRESSURE ) {
        // Important: In case of a general / quadratic EV problem, we must
        // ensure to have a "positive definite" matrix, i.e. we are not allowed
        // to multiply all matrices by -1!
        std::string stringFac = (analysistype_ != EIGENFREQUENCY) ? "-1.0" : "1.0";

        factor = CoefFunction::Generate(Global::REAL,
                                        CoefXprBinOp(dens, stringFac, CoefXpr::OP_MULT ) );
      } else {
        factor = CoefFunction::Generate(Global::REAL, "1.0");
      }
          
      // build coefficient for mass matrix as (factor / (c0*c0))
      PtrCoefFct coeffM =
          CoefFunction::Generate(Global::REAL,
                                 CoefXprBinOp( factor,
                                               CoefXprBinOp(c0,c0,CoefXpr::OP_MULT),
                                               CoefXpr::OP_DIV ) );

      // ====================================================================
      // Take account for pml (frequency domain only)
      // ====================================================================
      shared_ptr<CoefFunction> coeffPMLScal, coeffPMLVec;
      //just coeffunctions for the transformation of jacobians
      shared_ptr<CoefFunction> coeffPMLfactor;
      shared_ptr<CoefFunction> coeffPMLMass;
      std::string dampId;
      curRegNode->GetValue("dampingId",dampId);

      if(dampId != ""){
        if(analysistype_ == HARMONIC){
          PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal("pml","id",dampId.c_str());
          coeffPMLVec.reset(new CoefFunctionPML<Complex>(pmlNode,c0,actSDList,regions_,true));
          coeffPMLScal.reset(new CoefFunctionPML<Complex>(pmlNode,c0,actSDList,regions_,false));
          // store pml factor
          matCoefs_[PML_DAMP_FACTOR]->AddRegion(actRegion, coeffPMLVec);
          coeffPMLfactor = CoefFunction::Generate(Global::COMPLEX,
                                            CoefXprBinOp(factor,coeffPMLScal,CoefXpr::OP_MULT));
          coeffPMLMass = CoefFunction::Generate(Global::COMPLEX,
                                            CoefXprBinOp(coeffM,coeffPMLScal,CoefXpr::OP_MULT));
          harmonicPML = true;
        }else{
          EXCEPTION("PML is only available in the frequency domain right now");
        }
      }else{
        harmonicPML = false;
      }


      // ====================================================================
      // standard stiffness integrator
      // ====================================================================
      BaseBDBInt * stiffInt = NULL;
      if( dim_ == 2 ) {
        if(harmonicPML){
          stiffInt = new BBInt<Complex>(new ScaledGradientOperator<FeH1,2,Complex>(), coeffPMLfactor, 1.0 );
          stiffInt->SetBCoefFunctionOpB(coeffPMLVec);
        }else{
          stiffInt = new BBInt<Double>(new GradientOperator<FeH1,2>(), factor, 1.0 );
        }
      }
      else{
        if(harmonicPML){
          stiffInt = new BBInt<Complex>(new ScaledGradientOperator<FeH1,3,Complex>(),  coeffPMLfactor, 1.0 );
          stiffInt->SetBCoefFunctionOpB(coeffPMLVec);
        }else{
          stiffInt = new BBInt<Double>(new GradientOperator<FeH1,3>(), factor, 1.0 );
        }
      }
      stiffInt->SetName("LaplaceIntegrator");

      BiLinFormContext * stiffIntDescr =
        new BiLinFormContext(stiffInt, STIFFNESS );

      feFunctions_[formulation_]->AddEntityList( actSDList );

      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions(feFunctions_[formulation_],feFunctions_[formulation_]);
      stiffInt->SetFeSpace( feFunctions_[formulation_]->GetFeSpace());

      assemble_->AddBiLinearForm( stiffIntDescr );
      // Important: Add bdb-integrator to global list, as we need them later
      // for calculation of postprocessing results
      bdbInts_[actRegion] = stiffInt;

      // ====================================================================
      // standard mass integrator
      // ====================================================================

      BaseBDBInt *massInt = NULL;
      
      if(dim_==2){
        if(harmonicPML)
          massInt = new BBInt<Complex>(new IdentityOperator<FeH1,2,1,Complex>(), coeffPMLMass, 1.0 );
        else
          massInt = new BBInt<Double>(new IdentityOperator<FeH1,2,1,Double>,coeffM, 1.0 );
      }else{
        if(harmonicPML)
          massInt = new BBInt<Complex>(new IdentityOperator<FeH1,3,1,Complex>(), coeffPMLMass, 1.0 );
        else
          massInt = new BBInt<Double>(new IdentityOperator<FeH1,3,1,Double>, coeffM, 1.0 );
      }

      massInt->SetName("MassIntegrator");
      massInt->SetFeSpace( feFunctions_[formulation_]->GetFeSpace() );

      BiLinFormContext *massContext =  new BiLinFormContext(massInt, MASS );

      massContext->SetEntities( actSDList, actSDList );
      massContext->SetFeFunctions( feFunctions_[formulation_],feFunctions_[formulation_]);
      assemble_->AddBiLinearForm( massContext );
      // Important: Add mass-integrator to global list, as we need them later
      // for calculation of postprocessing results
      massInts_[actRegion] = massInt;

      // ====================================================================
      // check for flow (Pierce equation)
      // ====================================================================
      std::string flowId = curRegNode->Get("flowId")->As<std::string>();
      if(flowId != "") {
        if(dampId != ""){
          EXCEPTION("PML not available for flow domains!");
        }
        if ( formulation_ != ACOU_POTENTIAL )
          EXCEPTION("Pierce-Equation just possible in velocity potential formulation" );

        // Get result info object for flow
        shared_ptr<ResultInfo> flowInfo = GetResultInfo(MEAN_FLUIDMECH_VELOCITY); 
        
        //Add the region information
        PtrParamNode flowNode = myParam_->Get("flowList")->GetByVal("flow","name",flowId.c_str());

        // Read coefficient flow coefficient function for this region and add it to flow functor
        PtrCoefFct regionFlow;
        std::set<UInt> definedDofs;
        ReadUserFieldValues( regionName, flowNode, flowInfo->dofNames, flowInfo->entryType, 
                             isComplex_, regionFlow, definedDofs );
        meanFlowCoef_->AddRegion( actRegion, regionFlow );

        //now create the integrators
        BiLinearForm *convectiveStiff = NULL;
        BiLinearForm *convectiveDamp = NULL;
        if( dim_ == 2 ) {
          convectiveDamp  = new ABInt<>(new IdentityOperator<FeH1,2,1>(),
                                        new ConvectiveOperator<FeH1,2,1>(), coeffM, 2.0);
          convectiveStiff = new ABInt<>(new ConvectivePierceOperator<FeH1,2,1>(),
                                        new ConvectiveOperator<FeH1,2,1>(),coeffM, -1.0);
        } else {
          convectiveDamp  = new ABInt<>(new IdentityOperator<FeH1,3,1>(),
                                        new ConvectiveOperator<FeH1,3,1>(), coeffM, 2.0);
          convectiveStiff = new ABInt<>(new ConvectivePierceOperator<FeH1,3,1>(),
                                        new ConvectiveOperator<FeH1,3,1>(), coeffM, -1.0);
        }
        convectiveDamp->SetBCoefFunctionOpB(meanFlowCoef_);
        convectiveDamp->SetName("convectiveDampPierce");
        convectiveStiff->SetBCoefFunctionOpB(meanFlowCoef_);
        convectiveStiff->SetName("convectiveStiffPierce");

        BiLinFormContext *convectiveContextDamp  =  new BiLinFormContext(convectiveDamp, DAMPING );
        BiLinFormContext *convectiveContextStiff =  new BiLinFormContext(convectiveDamp, STIFFNESS );

        convectiveContextDamp->SetEntities( actSDList, actSDList );
        convectiveContextDamp->SetFeFunctions( feFunctions_[formulation_],feFunctions_[formulation_]);
        convectiveContextStiff->SetEntities( actSDList, actSDList );
        convectiveContextStiff->SetFeFunctions( feFunctions_[formulation_],feFunctions_[formulation_]);

        assemble_->AddBiLinearForm( convectiveContextDamp );
        assemble_->AddBiLinearForm( convectiveContextStiff );
      }
    }
  }

  void AcousticPDE::DefineSurfaceIntegrators( ){
    //========================================================================================
    // ABC boundaries
    //========================================================================================
    PtrParamNode bcNode = myParam_->Get( "bcsAndLoads", ParamNode::PASS );
    if( bcNode ) {
      ParamNodeList abcNodes = bcNode->GetList( "absorbingBCs" );

      for( UInt i = 0; i < abcNodes.GetSize(); i++ ) {
        std::string regionName = abcNodes[i]->Get("name")->As<std::string>();
        shared_ptr<EntityList> actSDList =  ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );
        std::string volRegName = abcNodes[i]->Get("volumeRegion")->As<std::string>();
        RegionIdType sRegId = ptGrid_->GetRegion().Parse(regionName);

        // --- Set the FE ansatz for the current region ---
        PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",volRegName.c_str());
        std::string polyId = curRegNode->Get("polyId")->As<std::string>();
        std::string integId = curRegNode->Get("integId")->As<std::string>();
        feFunctions_[formulation_]->GetFeSpace()->SetRegionApproximation(sRegId, polyId,integId);

        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);
        // c0 = sqrt(bulk_modulus / density)
        PtrCoefFct dens = materials_[aRegion]->GetScalCoefFnc( DENSITY, Global::REAL );
        PtrCoefFct blk = materials_[aRegion]->GetScalCoefFnc( ACOU_BULK_MODULUS, Global::REAL );
        PtrCoefFct c0 = 
            CoefFunction::Generate( Global::REAL,
                                    CoefXprUnaryOp( CoefXprBinOp(blk, dens, 
                                                                 CoefXpr::OP_DIV),
                                                    CoefXpr::OP_SQRT) );
        // In the case of acou-mech coupling we have to multiply the
        // abc-Integrator matrix with -1
        std::string factor = "1.0";
        if ( isMechCoupled_ == true && formulation_ !=  ACOU_PRESSURE ) {
          factor = "-1.0";
        }
        
        // factor for damping matrix: factor / c0
        PtrCoefFct coeffDamp = 
        CoefFunction::Generate(Global::REAL,
                               CoefXprBinOp(factor, c0, CoefXpr::OP_DIV ) );
        BiLinearForm * abcInt = NULL;
        if( dim_ == 2 ) {
          abcInt = new BBInt<>(new IdentityOperator<FeH1,2,1>(), coeffDamp, 1.0 );
        } else {
          abcInt = new BBInt<>(new IdentityOperator<FeH1,3,1>(), coeffDamp, 1.0 );
        }

        abcInt->SetName("abcIntegrator");
        BiLinFormContext *abcContext = new BiLinFormContext(abcInt, DAMPING );

        abcContext->SetEntities( actSDList, actSDList );
        abcContext->SetFeFunctions( feFunctions_[formulation_] , feFunctions_[formulation_]);
        feFunctions_[formulation_]->AddEntityList( actSDList );
        assemble_->AddBiLinearForm( abcContext );
      }
    }


  }

  void AcousticPDE::DefineRhsLoadIntegrators() {
    LOG_TRACE(acousticpde) << "Defining rhs load integrators for acoustic PDE";
    
    // Get FESpace and FeFunction of mechanical displacement
    shared_ptr<BaseFeFunction> myFct = feFunctions_[formulation_];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    StdVector<std::string> empty;
    LinearForm * lin = NULL;

    // obtain density
    shared_ptr<CoefFunctionMulti> densFct  = matCoefs_[ELEM_DENSITY];
    shared_ptr<CoefFunctionSurf> surfDens(new CoefFunctionSurf(false));
    surfDens->SetVolumeCoefs( densFct->GetRegionCoefs() );
    
           
    // In the case of acou-mech coupling we have to multiply the
    // integrators by -densiy
    Double scalFactor = 1.0;
    std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );
    
    // ===========================
    //  NORMAL VELOCITY (surface)
    // ===========================
    LOG_DBG(acousticpde) << "Reading normal velocity"; 
    
    ReadRhsExcitation( "normalVelocity", empty, ResultInfo::SCALAR, isComplex_,
                          ent, coef );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      
      // ensure that list contains only surface elements
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if( elemDim != (dim_-1) ) {
        EXCEPTION("Normal velocity can only be defined on surface elements");
      }
      PtrCoefFct exValue;
      if ( isMechCoupled_ == true && formulation_ !=  ACOU_PRESSURE ) {
        scalFactor = -1.0;
        exValue = 
            CoefFunction::Generate(part,
                                   CoefXprBinOp(coef[i],surfDens, CoefXpr::OP_MULT) );
      } else {
        exValue = coef[i];
      }
      
      if( formulation_ == ACOU_POTENTIAL ) {
        if( dim_ == 2) {
          if(isComplex_) {
            lin = new BUIntegrator<IdentityOperator<FeH1,2,1>,
                Complex,true>(scalFactor, exValue, volRegions);
          } else {
            lin = new BUIntegrator<IdentityOperator<FeH1,2,1>, 
                Double,true>(scalFactor, exValue, volRegions);
          }
        } else  {
          if(isComplex_) {
            lin = new BUIntegrator<IdentityOperator<FeH1,3,1>, 
                Complex,true>(1.0, exValue, volRegions);
          } else {
            lin = new BUIntegrator<IdentityOperator<FeH1,3,1>, 
                Double,true>(1.0, exValue , volRegions);
          }
        }
      } else if( formulation_ == ACOU_PRESSURE && isComplex_) {

        // in this case the pressure can be related to the 
        // normal velocity as 
        // p_n = - j*Omega*v_n*rho
        PtrCoefFct tmp = CoefFunction::Generate(part,"0.0",
                                                "-2*pi*f");
        PtrCoefFct tmp2 = 
            CoefFunction::Generate(part, 
                                   CoefXprBinOp(tmp, exValue, CoefXpr::OP_MULT) );
        exValue = CoefFunction::Generate(part, 
                                         CoefXprBinOp(tmp2,surfDens, 
                                                      CoefXpr::OP_MULT) );
        if( dim_ == 2) {
          lin = new BUIntegrator<IdentityOperator<FeH1,2,1>,
              Complex,true>(1.0, exValue, volRegions);
        } else  {
          lin = new BUIntegrator<IdentityOperator<FeH1,3,1>, 
              Complex,true>(1.0, exValue, volRegions);
        }

      } else {
        EXCEPTION( "Normal velocity can only be prescribed for potential "
            << "formulation or for harmonic pressure formulation" )
      }
      lin->SetName("NormalVelocityIntegrator");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    }

    
    // ===========================
    //  TOTAL VELOCITY (surface)
    // ===========================
    LOG_DBG(acousticpde) << "Reading total velocity";
    StdVector<std::string> vecDofNames;
    if(dim_ == 3)
      vecDofNames = "x", "y", "z";
    if(dim_ == 2 && !isaxi_)
      vecDofNames = "x", "y";
    if(dim_ == 2 && isaxi_)
      vecDofNames = "r", "z";
    ReadRhsExcitation( "velocity", vecDofNames, ResultInfo::VECTOR, isComplex_,
                       ent, coef );
    
    for( UInt i = 0; i < ent.GetSize(); ++i ) {

      // ensure that list contains only surface elements
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if( elemDim != (dim_-1) ) {
        EXCEPTION("Normal velocity can only be defined on surface elements");
      }

      PtrCoefFct exValue;
      if ( isMechCoupled_ == true && formulation_ !=  ACOU_PRESSURE ) {
        scalFactor = -1.0;
        exValue = 
            CoefFunction::Generate(part,
                                   CoefXprBinOp(coef[i],surfDens, CoefXpr::OP_MULT) );
      } else {
        exValue = coef[i];
      }
      if( formulation_ == ACOU_POTENTIAL ) {
        if( dim_ == 2) {
          if(isComplex_) {
            lin = new BUIntegrator<IdentityOperatorNormal<FeH1,2>,
                Complex,true>(scalFactor, exValue, volRegions);
          } else {
            lin = new BUIntegrator<IdentityOperatorNormal<FeH1,2>, 
                Double,true>(scalFactor, exValue, volRegions);
          }
        } else  {
          if(isComplex_) {
            lin = new BUIntegrator<IdentityOperatorNormal<FeH1,3>, 
                Complex,true>(scalFactor, exValue, volRegions);
          } else {
            lin = new BUIntegrator<IdentityOperatorNormal<FeH1,3>, 
                Double,true>(scalFactor, exValue, volRegions);
          }
        }
      } else if( formulation_ == ACOU_PRESSURE && isComplex_) {

        // in this case the pressure can be related to the 
        // normal velocity as 
        // p_n = - j*Omega*v_n*rho
        PtrCoefFct tmp = CoefFunction::Generate(part,"0.0",
                                                "-2*pi*f");
        PtrCoefFct tmp2 = 
            CoefFunction::Generate(part, 
                                   CoefXprBinOp(tmp, exValue, CoefXpr::OP_MULT) );
        exValue = CoefFunction::Generate(part, 
                                         CoefXprBinOp(surfDens, tmp2, 
                                                      CoefXpr::OP_MULT) );
        if( dim_ == 2) {
          lin = new BUIntegrator<IdentityOperatorNormal<FeH1,2>,
              Complex,true>(scalFactor, exValue, volRegions);
        } else  {
          lin = new BUIntegrator<IdentityOperatorNormal<FeH1,3>, 
              Complex,true>(scalFactor, exValue, volRegions);
        }
      } else {
        EXCEPTION( "Normal velocity can only be prescribed for potential "
            << "formulation or for harmonic pressure formulation" )
      }
      lin->SetName("VelocityIntegrator");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    }
    
    // ===========================
    //  PRESSURE (surface) 
    // ===========================
    // Only possible for harmonic potential formulation
    ReadRhsExcitation( "pressure", empty, ResultInfo::SCALAR, isComplex_,
                       ent, coef );
    if( formulation_ == ACOU_POTENTIAL && isComplex_ ) { 
      for( UInt i = 0; i < ent.GetSize(); ++i ) {

        // ensure that list contains only surface elements
        EntityIterator it = ent[i]->GetIterator();
        UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
        if( elemDim != (dim_-1) ) {
          EXCEPTION("Pressure can only be defined on surface elements");
        }
        PtrCoefFct exValue;
        if ( isMechCoupled_ == true && formulation_ !=  ACOU_PRESSURE ) {
          scalFactor = -1.0;
          exValue = 
              CoefFunction::Generate(part,
                                     CoefXprBinOp(coef[i],surfDens, CoefXpr::OP_MULT) );
        } else {
          exValue = coef[i];
        }
        // psi_n = j * 1 / (omega*rho) * p
        PtrCoefFct tmp = 
             CoefFunction::Generate(Global::COMPLEX, "0","1/(2*pi*f)");
         PtrCoefFct tmp2 = 
             CoefFunction::Generate(Global::COMPLEX,
                                    CoefXprBinOp(tmp, surfDens, CoefXpr::OP_DIV ) );
         exValue = CoefFunction::Generate(part, 
                                          CoefXprBinOp(exValue, tmp2, 
                                                       CoefXpr::OP_MULT) );
         shared_ptr<InhomDirichletBc> actBc ( new InhomDirichletBc );
         actBc->entities = ent[i];
         actBc->result = myFct->GetResultInfo();
         actBc->value = exValue;
         actBc->dofs.insert(0);
         myFct->AddInhomDirichletBc(actBc);
      } // loop entities
    } // if


  }

  void AcousticPDE::DefineSolveStep(){
    solveStep_ = new StdSolveStep(*this);
  }

  void AcousticPDE::DefinePrimaryResults(){

    // === Primary result according to definition ===
    shared_ptr<ResultInfo> res1( new ResultInfo);
    if ( formulation_ ==  ACOU_PRESSURE) {
      res1->resultType = ACOU_PRESSURE;
      res1->dofNames = "";
      res1->unit = "Pa";
    } else {
      res1->resultType = ACOU_POTENTIAL;
      res1->dofNames = "";
      res1->unit = "m^2/s";
    }
    res1->definedOn = ResultInfo::NODE;
    res1->entryType = ResultInfo::SCALAR;
    feFunctions_[formulation_]->SetResultInfo(res1);
    results_.Push_back( res1 );
    res1->SetFeFunction(feFunctions_[formulation_]);
    DefineFieldResult( feFunctions_[formulation_], res1 );
    
    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    if ( formulation_ ==  ACOU_PRESSURE) {
      hdbcSolNameMap_[ACOU_PRESSURE] = "soundSoft";
      idbcSolNameMap_[ACOU_PRESSURE] = "pressure";
    } else {
      hdbcSolNameMap_[ACOU_POTENTIAL] = "soundSoft";
      idbcSolNameMap_[ACOU_POTENTIAL] = "potential";
    }
    
    // === ACOUSTIC RHS ===
    shared_ptr<ResultInfo> rhs ( new ResultInfo );
    rhs->resultType = ACOU_RHS_LOAD;
    rhs->dofNames = "";
    rhs->unit = "?";
    rhs->definedOn = results_[0]->definedOn;
    rhs->entryType = ResultInfo::SCALAR;
    availResults_.insert( rhs );

    //creates the mean flow
    StdVector<std::string> vecDofNames;
    std::string geometryType;
    param->Get("domain")->GetValue("geometryType", geometryType );

    if( geometryType == "3d" ) {
      vecDofNames = "x", "y", "z";
    } else if( geometryType == "plane" ) {
      vecDofNames = "x", "y";
    } else if( geometryType == "axi" ) {
      vecDofNames = "r", "z";
    }

    CreateMeanFlowFunction(vecDofNames);
    
    // ==============================================
    //  Define CoefFunctions for material parameters
    // ==============================================
    
    // === DENSITY ===
    shared_ptr<ResultInfo> density ( new ResultInfo );
    density->resultType = ELEM_DENSITY;
    density->dofNames = "";
    density->unit = "kg/m^3";
    density->definedOn = ResultInfo::ELEMENT;
    density->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> densFct(new CoefFunctionMulti());
    matCoefs_[ELEM_DENSITY] = densFct;
    DefineFieldResult(densFct, density);
    
    // === SPEED OF SOUND ===
    shared_ptr<ResultInfo> sos ( new ResultInfo );
    sos->resultType = ACOU_ELEM_SPEED_OF_SOUND;
    sos->dofNames = "";
    sos->unit = "m/s";
    sos->definedOn = ResultInfo::ELEMENT;
    sos->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionMulti> sosFct(new CoefFunctionMulti());
    matCoefs_[ACOU_ELEM_SPEED_OF_SOUND] = sosFct;
    DefineFieldResult(sosFct, sos);
    
    // === PML DAMPING FACTORS ===
    //if( matCoefs_.find(PML_DAMP_FACTOR) != matCoefs_.end() ) {
    shared_ptr<ResultInfo> pml ( new ResultInfo );
    pml->resultType = PML_DAMP_FACTOR;
    pml->dofNames = vecDofNames;
    //pml->dofNames = "";
    pml->unit = "";
    pml->definedOn = ResultInfo::ELEMENT;
    pml->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionMulti> pmlFct(new CoefFunctionMulti());
    matCoefs_[PML_DAMP_FACTOR] = pmlFct;
    DefineFieldResult(pmlFct, pml);
    //}
  }
  
  void AcousticPDE::DefinePostProcResults(){
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    shared_ptr<BaseFeFunction> feFct = feFunctions_[formulation_];
    shared_ptr<ResultInfo> res1 = feFct->GetResultInfo();
    
    std::string geometryType;
    param->Get("domain")->GetValue("geometryType", geometryType );
    StdVector<std::string> vecDofNames;
    if( geometryType == "3d" ) {
      vecDofNames = "x", "y", "z";
    } else if( geometryType == "plane" ) {
      vecDofNames = "x", "y";
    } else if( geometryType == "axi" ) {
      vecDofNames = "r", "z";
    }
    
    // === PRESSURE / POTENTIAL - 1.DERIVATIVE ===
    shared_ptr<ResultInfo> deriv1(new ResultInfo);
    if( formulation_ == ACOU_POTENTIAL ) {
      deriv1->resultType = ACOU_POTENTIAL_DERIV_1;
      deriv1->dofNames = "";
      deriv1->unit = "m^2/s^2";
    } else {
      deriv1->resultType = ACOU_PRESSURE_DERIV_1;
      deriv1->dofNames = "";
      deriv1->unit = "Pa/s";
    }
    deriv1->entryType = res1->entryType;
    deriv1->definedOn = res1->definedOn;
    availResults_.insert( deriv1 );
    DefineTimeDerivResult( deriv1->resultType, 1, formulation_ );

    // === PRESSURE / POTENTIAL - 2.DERIVATIVE ===
    shared_ptr<ResultInfo> deriv2(new ResultInfo);
    if( formulation_ == ACOU_POTENTIAL ) {
      deriv2->resultType = ACOU_POTENTIAL_DERIV_2;
      deriv2->dofNames = "";
      deriv2->unit = "m^2/s^3";
    } else {
      deriv2->resultType = ACOU_PRESSURE_DERIV_2;
      deriv2->dofNames = "";
      deriv2->unit = "Pa/s^2";
    }
    deriv2->entryType = res1->entryType;
    deriv2->definedOn = res1->definedOn;
    availResults_.insert( deriv2 );
    DefineTimeDerivResult( deriv2->resultType, 2, formulation_ );

    
    // define acoustic pressure, depending on formulation
    PtrCoefFct presFct;
    if( formulation_ == ACOU_POTENTIAL ) {
      shared_ptr<ResultInfo> pres(new ResultInfo);
      pres->resultType = ACOU_PRESSURE;
      pres->dofNames = "";
      pres->unit = "Pa";
      pres->entryType = ResultInfo::SCALAR;
      pres->definedOn = ResultInfo::ELEMENT;
      // Define pressure as p = rho * dPsi/dt
      PtrCoefFct potD1Fct= this->GetCoefFct( ACOU_POTENTIAL_DERIV_1 );
      PtrCoefFct densFct = this->GetCoefFct( ELEM_DENSITY);
      presFct = 
          CoefFunction::Generate(part,
                                 CoefXprBinOp(potD1Fct, densFct, CoefXpr::OP_MULT ) );
      DefineFieldResult(presFct, pres); 
    } else{
      presFct = this->GetCoefFct(ACOU_PRESSURE);
    }
    
    
    shared_ptr<ResultInfo> vel, velNormal, intensity, power, pres;
    PtrCoefFct intensFct, velFct;
    shared_ptr<CoefFunctionSurf> sNormIntens, velFctNormal;
    shared_ptr<CoefFunctionFormBased>  presGradFct, velFctPot;
    shared_ptr<ResultFunctor> powerFct;
    
    // some results are only available in potential and / or
    // harmonic simulation
    if( formulation_ == ACOU_POTENTIAL ||
        isComplex_ ){
      // === ACOU_PRESSURE ===
      pres.reset(new ResultInfo);
      pres->resultType = ACOU_PRESSURE;
      pres->dofNames = "";
      pres->unit = "Pa";
      pres->entryType = ResultInfo::SCALAR;
      pres->definedOn = ResultInfo::ELEMENT;
      
      // === ACOU_VELOCITY ===
      vel.reset(new ResultInfo);
      vel->resultType = ACOU_VELOCITY;
      vel->dofNames = vecDofNames;
      vel->unit = "m/s";
      vel->entryType = ResultInfo::VECTOR;
      vel->definedOn = ResultInfo::ELEMENT;
      
      // === ACOU_NORMAL_VELOCITY ===
      velNormal.reset(new ResultInfo);
      velNormal->resultType = ACOU_NORMAL_VELOCITY;
      velNormal->dofNames = "";
      velNormal->unit = "m/s";
      velNormal->entryType = ResultInfo::SCALAR;
      velNormal->definedOn = ResultInfo::SURF_ELEM;
      
      // === ACOU_INTENSITY ===
      intensity.reset(new ResultInfo);
      intensity->resultType = ACOU_INTENSITY;
      intensity->dofNames = vecDofNames;
      intensity->unit = "W/m^2";
      intensity->entryType = ResultInfo::VECTOR;
      intensity->definedOn = ResultInfo::ELEMENT;
      
      // === ACOU_POWER ===
      power.reset(new ResultInfo);
      power->resultType = ACOU_POWER;
      power->dofNames = "";
      power->unit = "W";
      power->entryType = ResultInfo::SCALAR;
      power->definedOn = ResultInfo::SURF_REGION;
    }
    
    if( formulation_ == ACOU_POTENTIAL ) {
      // --------------------------
      //  POTENTENTIAL FORMULATION
      // --------------------------
      // === ACOU_VELOCITY ===
      // Velocity v = - grad Psi
      if( isComplex_ ) {
        velFctPot.reset(new CoefFunctionBOp<Complex>(feFct, vel, -1.0));
      } else {
        velFctPot.reset(new CoefFunctionBOp<Double>(feFct, vel, -1.0));
      }
      velFct = velFctPot;
      DefineFieldResult( velFct, vel );
      
      // === ACOU_NORMAL_VELOCITY ===
      velFctNormal.reset(new CoefFunctionSurf(true));
      DefineFieldResult(velFctNormal, velNormal);
      
      // === ACOU_INTENSITY ===
      // Intensity I = p * conj(v)
      intensFct = 
          CoefFunction::Generate(part,
                                 CoefXprBinOp(presFct, velFct, CoefXpr::OP_MULT_CONJ ) );
      DefineFieldResult(intensFct, intensity);
      
      // === ACOU_POWER ===
      // Power p = \int_Gamma I *n dGamma
      // a) first, define surface normal intensity
      sNormIntens.reset(new CoefFunctionSurf(true));
      
      // b) then, integrate values
      if( isComplex_ ) {
        powerFct.reset(new ResultFunctorIntegrate<Complex>(sNormIntens, 
                                                            feFct, power ) );
      } else {
        powerFct.reset(new ResultFunctorIntegrate<Double>(sNormIntens, 
                                                           feFct, power ) );
      }
      resultFunctors_[ACOU_POWER] = powerFct;
      availResults_.insert(power);
    }
    
    if( formulation_ == ACOU_PRESSURE && isComplex_ ) {
      // --------------------------------
      //  COMPLEX & PRESSURE FORMULATION
      // --------------------------------
      
      // === ACOU_VELOCITY ===
      
      // Velocity v = j* 1/(omega*rho) * grad(p)
      // a) define gradient of pressure
      if( isComplex_ ) {
        presGradFct.reset(new CoefFunctionBOp<Complex>(feFct, vel, 1.0));  
      } else {
        presGradFct.reset(new CoefFunctionBOp<Double>(feFct, vel, 1.0));
      }
      
      // b) multiply by factor
      PtrCoefFct densFct = this->GetCoefFct( ELEM_DENSITY);
      PtrCoefFct factor = 
          CoefFunction::Generate(Global::COMPLEX, "0","1/(2*pi*f)");
      PtrCoefFct factor2 = 
          CoefFunction::Generate(Global::COMPLEX,
                                 CoefXprBinOp(factor, densFct, CoefXpr::OP_DIV ) );
      velFct = 
          CoefFunction::Generate( part,
                                  CoefXprBinOp( factor2, presGradFct, CoefXpr::OP_MULT ) );
      DefineFieldResult( velFct, vel );

      // === ACOU_NORMAL_VELOCITY ===
      velFctNormal.reset(new CoefFunctionSurf(true));
      DefineFieldResult(velFctNormal, velNormal);
      
      // === ACOU_INTENSITY ===
      // Intensity I = p * v
      intensFct = 
          CoefFunction::Generate(part,
                                 CoefXprBinOp(presFct, velFct, CoefXpr::OP_MULT_CONJ ) );
      DefineFieldResult(intensFct, intensity);

      // === ACOU_POWER ===
      // Power p = \int_Gamma I *n dGamma
      // a) first, define surface normal intensity
      sNormIntens.reset(new CoefFunctionSurf(true));

      // b) then, integrate values
      powerFct.reset(new ResultFunctorIntegrate<Complex>(sNormIntens, 
                                                         feFct, power ) );
      resultFunctors_[ACOU_POWER] = powerFct;
      availResults_.insert(power);
    }
    
    // ============================
    // Initialize result functors:
    // ============================
    // collect the flux volume flux coefficient function for each region
    std::map<RegionIdType, PtrCoefFct> intensCoefs, velCoefs;

    // 1) Loop over all BDB-integrators
    std::map<RegionIdType, BaseBDBInt*>::iterator it = bdbInts_.begin();
    for( ; it != bdbInts_.end(); ++it ) {
      RegionIdType region = it->first;
      BaseBDBInt* bdb = it->second;
      if( velFctPot)
        velFctPot->AddIntegrator(bdb, region);
      if( presGradFct)
        presGradFct->AddIntegrator(bdb, region);
      intensCoefs[region] = intensFct;
      velCoefs[region] = velFct;
    }
    
    if( intensFct) {
      sNormIntens->SetVolumeCoefs(intensCoefs);
    }
    if( velFct) {
      velFctNormal->SetVolumeCoefs(velCoefs);
    }
  }

   void AcousticPDE::CreateMeanFlowFunction(StdVector<std::string> dofNames){
     //// === MEAN FLUIDMECH VELOCITY ===
     shared_ptr<ResultInfo> flowvelocity( new ResultInfo);
     flowvelocity->resultType = MEAN_FLUIDMECH_VELOCITY;
     flowvelocity->dofNames = dofNames;
     flowvelocity->unit = "m/s";

     flowvelocity->definedOn = ResultInfo::NODE;
     flowvelocity->entryType = ResultInfo::VECTOR;
     
     meanFlowCoef_.reset(new CoefFunctionMulti());
     DefineFieldResult( meanFlowCoef_, flowvelocity );

     results_.Push_back( flowvelocity );
     availResults_.insert( flowvelocity );

   }

  //! Init the time stepping
  void AcousticPDE::InitTimeStepping(){
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(TimeSchemeGLM::NEWMARK, 0) );

    feFunctions_[formulation_]->SetTimeScheme(myScheme);

  }
}
