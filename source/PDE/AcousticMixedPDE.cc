// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     AcousticMixedPDE.cc
 *       \brief    <Description>
 *
 *       \date     Feb 2, 2012
 *       \author   ahueppe
 */
//================================================================================================


#include "AcousticMixedPDE.hh"

#include <boost/lexical_cast.hpp>
#include <cmath>

#include "General/defs.hh"

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/ConvectiveOperator.hh"
#include "Forms/Operators/NormalFluxOperators.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"

#include "FeBasis/FeFunctions.hh"

#include "FeBasis/H1/FeSpaceH1.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"
#include "FeBasis/L2/FeSpaceL2.hh"
#include "FeBasis/L2/FeSpaceL2Nodal.hh"


#include "Domain/Results/ResultFunctor.hh"
#include "Domain/Results/ExternalFieldFunctors.hh"

#include "Driver/Assemble.hh"

#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Materials/AcousticMaterial.hh"

namespace CoupledField{

  DECLARE_LOG(acousticmixedpde)
   DEFINE_LOG(acousticmixedpde, "pde.acousticmixed")

   AcousticMixedPDE::AcousticMixedPDE( Grid* aGrid, PtrParamNode paramNode )
               : SinglePDE( aGrid, paramNode ){

     pdename_           = "acousticMixed";
     pdematerialclass_  = FLUID;
     maxTimeDerivOrder_ = 1;
     nonLin_            = false;
     InitialCondition_  = 0.0;
     usePiola_          = false;
     penalized_         = true;
     doFluxTerm_        = true;

   }

  std::map<SolutionType, shared_ptr<FeSpace> >
  AcousticMixedPDE::CreateFeSpaces( const std::string&  formulation,
                    PtrParamNode infoNode ){

    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "MixedPiola" || formulation == "Mixed"){
      std::string form = SolutionTypeEnum.ToString(ACOU_PRESSURE);
      PtrParamNode potSpaceNode = infoNode->Get(form);
      crSpaces[ACOU_PRESSURE] =
              FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1);
      crSpaces[ACOU_VELOCITY] =
          FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::L2);

      crSpaces[ACOU_PRESSURE]->Init(solStrat_);
      crSpaces[ACOU_VELOCITY]->Init(solStrat_);
    }else{
      EXCEPTION("The formulation " << formulation << "of acousticMixed PDE is not known!");
    }
    usePiola_ = myParam_->Get("usePiolaTransform")->As<bool>();
    penalized_ = myParam_->Get("penalized")->As<bool>();
    doFluxTerm_ = myParam_->Get("fluxTerm")->As<bool>();
    return crSpaces;
  }

  void AcousticMixedPDE::DefineIntegrators(){
    if(dim_==2){
      if(isComplex_)
        DefineIntegratorsTempl<Complex,2>();
      else
        DefineIntegratorsTempl<Double,2>();
    }else if(dim_==3){
      if(isComplex_)
        DefineIntegratorsTempl<Complex,2>();
      else
        DefineIntegratorsTempl<Double,3>();
    }
  }

  template<class DATA_TYPE, UInt DIM>
  void AcousticMixedPDE::DefineIntegratorsTempl(){

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
    shared_ptr<FeSpace> spaceP = feFunctions_[ACOU_PRESSURE]->GetFeSpace();
    shared_ptr<FeSpace> spaceV = feFunctions_[ACOU_VELOCITY]->GetFeSpace();

    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      // Set current region and material
      actRegion = it->first;
      // Get current region name
      std::string regionName = ptgrid_->GetRegion().ToString(actRegion);

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );

      // --- Set the FE ansatz for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      spaceP->SetRegionApproximation(actRegion, polyId,integId);
      spaceV->SetRegionApproximation(actRegion, polyId,integId);

      Double density=0.0;
      Double compressibility=0.0;

      materials_[actRegion]->GetScalar(density,DENSITY,Global::REAL);
      materials_[actRegion]->GetScalar(compressibility,ACOU_BULK_MODULUS,Global::REAL);

    // Not needed at the moment. Commented out due to gcc 4.6.
#if 0
      Double c0=0.0;
      c0 = std::sqrt(compressibility/density);
#endif

      feFunctions_[ACOU_PRESSURE]->AddEntityList( actSDList );
      feFunctions_[ACOU_VELOCITY]->AddEntityList( actSDList );


      // ====================================================================
      // stiffness integrators
      // ====================================================================
      // --------------------------------------------------------------------
      //  VERSION 1: K_PV Integrator (upper off-diagonal integrator)
      // --------------------------------------------------------------------

      shared_ptr<CoefFunction> coeffKPV
                = CoefFunction::Generate(Global::REAL, "1.0");
      BiLinearForm * stiffIntPV = NULL;

      if(usePiola_)
        stiffIntPV = new ABInt<GradientOperator<FeH1,DIM,DATA_TYPE>,IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE> , DATA_TYPE >(coeffKPV,1.0 );
      else
        stiffIntPV = new ABInt< GradientOperator<FeH1,DIM,DATA_TYPE> , IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>,  DATA_TYPE >(coeffKPV,1.0 );

      stiffIntPV->SetName("MixedStiffIntPV");
      //
      BiLinFormContext *stiffContPV = new BiLinFormContext(stiffIntPV, STIFFNESS );
      
      stiffContPV->SetEntities( actSDList, actSDList );
      stiffContPV->SetFeFunctions( feFunctions_[ACOU_PRESSURE],feFunctions_[ACOU_VELOCITY]);

      //stiffContPV->SetCounterPart(true);
      assemble_->AddBiLinearForm( stiffContPV );

      
      // --------------------------------------------------------------------
      //  VERSION 2: K_VP = K_PV^T Integrator (lower off-diagonal integrator)
      // --------------------------------------------------------------------
      shared_ptr<CoefFunction> coeffKVP
                = CoefFunction::Generate(Global::REAL, "1.0");
      BiLinearForm * stiffIntVP = NULL;

      if(usePiola_)
        stiffIntVP = new ABInt< IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE> , //
                                GradientOperator<FeH1,DIM,DATA_TYPE>, //
                                DATA_TYPE >(coeffKVP,-1.0 );
      else
        stiffIntVP = new ABInt< IdentityOperator<FeH1,DIM,DIM,DATA_TYPE> , //
                                GradientOperator<FeH1,DIM,DATA_TYPE>,  //
                                DATA_TYPE  >(coeffKVP,-1.0 );

      stiffIntVP->SetName("MixedStiffIntVP");
      BiLinFormContext *stiffContVP = new BiLinFormContext(stiffIntVP, STIFFNESS );

      stiffContVP->SetEntities( actSDList, actSDList );
      stiffContVP->SetFeFunctions( feFunctions_[ACOU_VELOCITY],feFunctions_[ACOU_PRESSURE]);
      assemble_->AddBiLinearForm( stiffContVP );

      // ====================================================================
      // mass integrators
      // ====================================================================
      shared_ptr<CoefFunction> coeffMPP
                = CoefFunction::Generate(Global::REAL, lexical_cast<std::string>(1.0));

      BiLinearForm *massIntPP = NULL;

      massIntPP = new BBInt<IdentityOperator<FeH1,DIM,1,DATA_TYPE> , DATA_TYPE, DATA_TYPE>(coeffMPP, 1.0 / (compressibility) );

      massIntPP->SetName("PressureTimeInt");
      //massIntPP->SetFeSpace( feFunctions_[ACOU_PRESSURE]->GetFeSpace() );

      BiLinFormContext *massContextPP =  new BiLinFormContext(massIntPP, DAMPING );

      massContextPP->SetEntities( actSDList, actSDList );
      massContextPP->SetFeFunctions( feFunctions_[ACOU_PRESSURE],feFunctions_[ACOU_PRESSURE]);
      assemble_->AddBiLinearForm( massContextPP );

      shared_ptr<CoefFunction> coeffMVV
                = CoefFunction::Generate(Global::REAL, lexical_cast<std::string>(1.0));

      BiLinearForm *massIntVV = NULL;
      if( usePiola_)
        massIntVV = new BBInt<IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>, DATA_TYPE, DATA_TYPE >(coeffMVV, density);
      else
        massIntVV = new BBInt<IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>, DATA_TYPE, DATA_TYPE >(coeffMVV, density);

      massIntVV->SetName("VelocityTimeInt");
      //massIntVV->SetFeSpace( feFunctions_[ACOU_VELOCITY]->GetFeSpace() );

      BiLinFormContext *massContextVV =  new BiLinFormContext(massIntVV, DAMPING );

      massContextVV->SetEntities( actSDList, actSDList );
      massContextVV->SetFeFunctions( feFunctions_[ACOU_VELOCITY],feFunctions_[ACOU_VELOCITY]);
      assemble_->AddBiLinearForm( massContextVV );

      //======================================================================
      // CHECK FOR FLOW
      //=====================================================================
      std::string flowId = curRegNode->Get("flowId")->As<std::string>();
      if(flowId != ""){
        //Add the region information
        PtrParamNode flowNode = myParam_->Get("flowList")->GetByVal("flow","name",flowId.c_str());
        if(isComplex_){
          shared_ptr<FieldFunctor<Complex> > fct = dynamic_pointer_cast<FieldFunctor<Complex> >(meanFlowFunctor_);
          fct->AddRegion(actRegion,flowNode);
        }else{
          shared_ptr<FieldFunctor<Double> > fct = dynamic_pointer_cast<FieldFunctor<Double> >(meanFlowFunctor_);
          fct->AddRegion(actRegion,flowNode);
        }

        //now create the integrators
        BiLinearForm *convectiveVV = NULL;
        BiLinearForm *convectivePP = NULL;
        Double convFactor = density;
        if(doFluxTerm_)
          convFactor *= 0.5;

        if(usePiola_){
          convectiveVV = new ABInt<IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>, ConvectiveOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>, DATA_TYPE >(coeffMVV, convFactor);
        } else {
          convectiveVV = new ABInt<IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>, ConvectiveOperator<FeH1,DIM,DIM,DATA_TYPE>, DATA_TYPE >(coeffMVV, convFactor);
        }
        convectivePP = new ABInt<IdentityOperator<FeH1,DIM,1,DATA_TYPE>, ConvectiveOperator<FeH1,DIM,1,DATA_TYPE>, DATA_TYPE >(coeffMPP, 1.0 / (compressibility));

        convectiveVV->SetBCoefFunctionOpB(meanFlowFunctor_);
        convectivePP->SetBCoefFunctionOpB(meanFlowFunctor_);

        convectiveVV->SetName("convectiveVV");
        convectivePP->SetName("convectivePP");

        BiLinFormContext *convectiveContextVV =  new BiLinFormContext(convectiveVV, STIFFNESS );
        BiLinFormContext *convectiveContextPP =  new BiLinFormContext(convectivePP, STIFFNESS );

        convectiveContextVV->SetEntities( actSDList, actSDList );
        convectiveContextVV->SetFeFunctions( feFunctions_[ACOU_VELOCITY],feFunctions_[ACOU_VELOCITY]);
        convectiveContextPP->SetEntities( actSDList, actSDList );
        convectiveContextPP->SetFeFunctions( feFunctions_[ACOU_PRESSURE],feFunctions_[ACOU_PRESSURE]);
        assemble_->AddBiLinearForm( convectiveContextVV );
        assemble_->AddBiLinearForm( convectiveContextPP );

        //============================================================================================
        // ADD THE PENALIZATION TERM
        //============================================================================================
        if(penalized_){
          //two forms one with opposing elements the other just like mass...
          Double formFactor = density * curRegNode->Get("penalizationFactor")->As<Double>();
          BiLinearForm *convectiveVOpp = NULL;
          BiLinearForm *convectiveV = NULL;
          BiLinearForm *exteriorVV = NULL;
          std::set<RegionIdType> volRegion;
          volRegion.insert(actRegion);
          if( usePiola_ ) {
            convectiveVOpp = new SurfaceBBInt<IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>,DATA_TYPE, DATA_TYPE >(coeffMVV, -0.5*formFactor,volRegion);
            convectiveV    = new SurfaceBBInt<IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>,DATA_TYPE, DATA_TYPE >(coeffMVV, 0.5*formFactor,volRegion);
            exteriorVV    = new SurfaceBBInt<IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>,DATA_TYPE, DATA_TYPE >(coeffMVV, -0.5*formFactor,volRegion);
          } else {
            convectiveVOpp = new BBInt<IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>,DATA_TYPE, DATA_TYPE >(coeffMVV, -0.5 * formFactor);
            convectiveV    = new BBInt<IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>,DATA_TYPE, DATA_TYPE >(coeffMVV, 0.5*formFactor);
            exteriorVV    = new BBInt<IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>,DATA_TYPE, DATA_TYPE >(coeffMVV, -0.5*formFactor);
          }
          convectiveV->SetName("penaltyMass");
          convectiveVOpp->SetName("penalityOpposite");
          exteriorVV->SetName("exteriorRegion");
          //now we obtain the entity lists
          shared_ptr<EntityList> list,oppositList,extList;
          feFunctions_[ACOU_VELOCITY]->GetFeSpace()->GetInteriorSurfaceElems(actRegion,list,oppositList);
          feFunctions_[ACOU_VELOCITY]->GetFeSpace()->GetExteriorSurfaceElems(actRegion,extList);

          BiLinFormContext * penaltyContextVOpp =  new BiLinFormContext(convectiveVOpp, STIFFNESS );
          BiLinFormContext * penaltyContextV =  new BiLinFormContext(convectiveV, STIFFNESS );
          BiLinFormContext * penaltyContextExt =  new BiLinFormContext(exteriorVV, STIFFNESS );

          penaltyContextVOpp->SetEntities( list, oppositList );
          penaltyContextV->SetEntities( list, list );
          penaltyContextExt->SetEntities(extList,extList);

          penaltyContextVOpp->SetFeFunctions( feFunctions_[ACOU_VELOCITY],feFunctions_[ACOU_VELOCITY]);
          penaltyContextV->SetFeFunctions( feFunctions_[ACOU_VELOCITY],feFunctions_[ACOU_VELOCITY]);
          penaltyContextExt->SetFeFunctions( feFunctions_[ACOU_VELOCITY],feFunctions_[ACOU_VELOCITY]);
          assemble_->AddBiLinearForm( penaltyContextVOpp );
          assemble_->AddBiLinearForm( penaltyContextV );
          assemble_->AddBiLinearForm( penaltyContextExt );

        }

        if(doFluxTerm_){
          BiLinearForm *fluxTerm = NULL;
          BiLinearForm *convectiveVVTrans = NULL;
          std::set<RegionIdType> volRegion;
          volRegion.insert(actRegion);

          if( usePiola_ ) {
            fluxTerm = new SurfaceABInt< IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>,//
                                         NormalFluxOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>,//
                                         DATA_TYPE >(coeffMVV, -0.5 * density,volRegion);

            convectiveVVTrans = new ABInt< ConvectiveOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>,//
                                           IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>,//
                                           DATA_TYPE >(coeffMVV, -0.5 * density);
          } else {
            fluxTerm = new SurfaceABInt<IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>,//
                                        NormalFluxOperator<FeH1,DIM,DIM,DATA_TYPE>,//
                                        DATA_TYPE >(coeffMVV, -0.5 * density,volRegion);

            convectiveVVTrans = new ABInt< ConvectiveOperator<FeH1,DIM,DIM,DATA_TYPE>, //
                                           IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>, //
                                           DATA_TYPE >(coeffMVV, -0.5 * density);
          }

          convectiveVVTrans->SetBCoefFunctionOpA(meanFlowFunctor_);
          fluxTerm->SetBCoefFunctionOpB(meanFlowFunctor_);

          fluxTerm->SetName("fluxTerm");
          convectiveVVTrans->SetName("convectiveVVtrans");

          //now we obtain the entity lists
          shared_ptr<EntityList> list,oppositList,extList;
          feFunctions_[ACOU_VELOCITY]->GetFeSpace()->GetInteriorSurfaceElems(actRegion,list,oppositList);
          feFunctions_[ACOU_VELOCITY]->GetFeSpace()->GetExteriorSurfaceElems(actRegion,extList);

          BiLinFormContext * fluxContext =  new BiLinFormContext(fluxTerm, STIFFNESS );
          BiLinFormContext * convectiveTransContext =  new BiLinFormContext(convectiveVVTrans, STIFFNESS );

          fluxContext->SetEntities( list, oppositList );
          convectiveTransContext->SetEntities( actSDList, actSDList );

          fluxContext->SetFeFunctions( feFunctions_[ACOU_VELOCITY],feFunctions_[ACOU_VELOCITY]);
          convectiveTransContext->SetFeFunctions( feFunctions_[ACOU_VELOCITY],feFunctions_[ACOU_VELOCITY]);
          assemble_->AddBiLinearForm( fluxContext );
          assemble_->AddBiLinearForm( convectiveTransContext );

        }
      }
    }
  }

   void AcousticMixedPDE::DefineRhsLoadIntegrators(){
     if(dim_==2){
       if(isComplex_)
         DefineRhsLoadIntegratorsTempl<Complex,2>();
       else
         DefineRhsLoadIntegratorsTempl<Double,2>();
     }else if(dim_==3){
       if(isComplex_)
         DefineRhsLoadIntegratorsTempl<Complex,3>();
       else
         DefineRhsLoadIntegratorsTempl<Double,3>();
     }
   }

  template<class DATA_TYPE,UInt DIM>
  void AcousticMixedPDE::DefineRhsLoadIntegratorsTempl() {
    LOG_TRACE(acousticmixedpde) << "Defining rhs load integrators for acoustic mixed PDE";

    // Get FESpace and FeFunctions
    shared_ptr<BaseFeFunction> pressureFct = feFunctions_[ACOU_PRESSURE];
    shared_ptr<BaseFeFunction> velocityFct = feFunctions_[ACOU_VELOCITY];
    shared_ptr<FeSpace> pSpace = pressureFct->GetFeSpace();
    shared_ptr<FeSpace> vSpace = velocityFct->GetFeSpace();

    StdVector<shared_ptr<EntityList> > ent;
    StdVector<shared_ptr<CoefFunction> > coef;
    LinearForm * lin = NULL;
    StdVector<std::string> vDofNames = velocityFct->GetResultInfo()->dofNames;


    LOG_DBG(acousticmixedpde) << "Reading loads for mass conservation equation";
    ReadRhsExcitation( "massEquationLoad", pressureFct->GetResultInfo()->dofNames, ResultInfo::SCALAR, isComplex_, ent, coef );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {

      lin = new BUIntegrator<IdentityOperator<FeH1,DIM,1,DATA_TYPE>, DATA_TYPE >(1.0,coef[i]);

      lin->SetName("massEquationInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(pressureFct);
      assemble_->AddLinearForm(ctx);
    }

    LOG_DBG(acousticmixedpde) << "Reading loads for momentum conservation equation";
    ReadRhsExcitation( "momentumEquationLoad", vDofNames, ResultInfo::VECTOR, isComplex_, ent, coef );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      if(usePiola_){
        lin = new BUIntegrator<IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>, DATA_TYPE >(1.0,coef[i]);
      }else{
        lin = new BUIntegrator<IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>, DATA_TYPE >(1.0,coef[i]);
      }

      lin->SetName("momentumEquationInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(velocityFct);
      assemble_->AddLinearForm(ctx);
    }
  }
  void AcousticMixedPDE::DefineSurfaceIntegrators( ){
    //========================================================================================
    // ABC boundaries
    //========================================================================================
    PtrParamNode bcNode = myParam_->Get( "bcsAndLoads", ParamNode::PASS );
    if( bcNode ) {
      ParamNodeList abcNodes = bcNode->GetList( "absorbingBCs" );

      for( UInt i = 0; i < abcNodes.GetSize(); i++ ) {
        std::string regionName = abcNodes[i]->Get("name")->As<std::string>();
        shared_ptr<EntityList> actSDList =  ptgrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );
        std::string volRegName = abcNodes[i]->Get("volumeRegion")->As<std::string>();
        RegionIdType sRegId = ptgrid_->GetRegion().Parse(regionName);

        // --- Set the FE ansatz for the current region ---
        PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",volRegName.c_str());
        std::string polyId = curRegNode->Get("polyId")->As<std::string>();
        std::string integId = curRegNode->Get("integId")->As<std::string>();
        feFunctions_[ACOU_PRESSURE]->GetFeSpace()->SetRegionApproximation(sRegId, polyId,integId);

        RegionIdType aRegion = ptgrid_->GetRegion().Parse(volRegName);
        //
        Double density=0.0,compressibility=0.0,c0;

        //
        materials_[aRegion]->GetScalar(density,DENSITY,Global::REAL);
        materials_[aRegion]->GetScalar(compressibility,ACOU_BULK_MODULUS,Global::REAL);
        c0 = std::sqrt(compressibility/density);

        shared_ptr<CoefFunction> coeffKPV
                  = CoefFunction::Generate(Global::REAL, "1.0");

        BiLinearForm * abcInt = NULL;

        if( dim_ == 2 ) {
          if(isComplex_)
            abcInt = new BBInt<IdentityOperator<FeH1,2,1,Complex>, Complex,Complex >(coeffKPV,1.0/(density*c0) );
          else
            abcInt = new BBInt<IdentityOperator<FeH1,2,1,Double>, Double,Double >(coeffKPV,1.0/(density*c0) );
        } else {
          if(isComplex_)
            abcInt = new BBInt<IdentityOperator<FeH1,3,1,Complex>, Complex,Complex >(coeffKPV,1.0/(density*c0) );
          else
            abcInt = new BBInt<IdentityOperator<FeH1,3,1,Double>, Double,Double >(coeffKPV,1.0/(density*c0) );
        }

        abcInt->SetName("abcIntegrator");
        BiLinFormContext *abcContext = new BiLinFormContext(abcInt, STIFFNESS );

        abcContext->SetEntities( actSDList, actSDList );
        abcContext->SetFeFunctions( feFunctions_[ACOU_PRESSURE] , feFunctions_[ACOU_PRESSURE]);
        feFunctions_[ACOU_PRESSURE]->AddEntityList( actSDList );
        assemble_->AddBiLinearForm( abcContext );
      }
    }
  }

  void AcousticMixedPDE::DefineSolveStep(){
      solveStep_ = new StdSolveStep(*this);
    }

    void AcousticMixedPDE::SetInitialCondition(){

    }

    void AcousticMixedPDE::CalcResults( shared_ptr<BaseResult> result ){
      switch (result->GetResultInfo()->resultType ) {
      case ACOU_PRESSURE:
        feFunctions_[ACOU_PRESSURE]->ExtractResult( result );
        break;
      case ACOU_VELOCITY:
        feFunctions_[ACOU_PRESSURE]->ExtractResult( result );
        break;
      case MEAN_FLUIDMECH_VELOCITY:
        meanFlowFunctor_->EvalResult( result );
        break;
      default:
        WARN( "Resulttype not computable by acoustic PDE" );
        break;
      }
    }

    //!  Define available postprocessing results
    void AcousticMixedPDE::DefinePrimaryResults(){

      // Check for subType
      StdVector<std::string> velDofNames;

      std::string geometryType;
      param->Get("domain")->GetValue("geometryType", geometryType );

      if( geometryType == "3d" ) {
        velDofNames = "x", "y", "z";
      } else if( geometryType == "plane" ) {
        velDofNames = "x", "y";
      } else if( geometryType == "axi" ) {
        velDofNames = "r", "z";
      }

      // === Primary result according to definition ===
      //PRESSURE
      shared_ptr<ResultInfo> pressure( new ResultInfo);
      pressure->resultType = ACOU_PRESSURE;
      pressure->dofNames = "";
      pressure->unit = "Pa";

      pressure->definedOn = ResultInfo::NODE;
      pressure->entryType = ResultInfo::SCALAR;
      feFunctions_[ACOU_PRESSURE]->SetResultInfo(pressure);
      results_.Push_back( pressure );
      availResults_.insert( pressure );

      pressure->SetFeFunction(feFunctions_[ACOU_PRESSURE]);

      //VELOCITY
      shared_ptr<ResultInfo> velocity( new ResultInfo);
      velocity->resultType = ACOU_VELOCITY;
      velocity->dofNames = velDofNames;
      velocity->unit = "m/s";

      velocity->definedOn = ResultInfo::NODE;
      velocity->entryType = ResultInfo::VECTOR;
      feFunctions_[ACOU_VELOCITY]->SetResultInfo(velocity);
      results_.Push_back( velocity );
      availResults_.insert( velocity );

      velocity->SetFeFunction(feFunctions_[ACOU_VELOCITY]);

      //// === ACOUSTIC RHS ===
      //shared_ptr<ResultInfo> rhs ( new ResultInfo );
      //rhs->resultType = ACOU_RHS_LOAD;
      //rhs->dofNames = "";
      //rhs->unit = "?";
      //rhs->definedOn = results_[0]->definedOn;
      //rhs->entryType = ResultInfo::SCALAR;
      //availResults_.insert( rhs );
      //postProcResults_[ACOU_RHS_LOAD] = ACOU_RHS_LOAD;


      CreateMeanFlowFunction(velDofNames);
   }

   void AcousticMixedPDE::CreateMeanFlowFunction(StdVector<std::string> dofNames){
     //// === MEAN FLUIDMECH VELOCITY ===
     shared_ptr<ResultInfo> flowvelocity( new ResultInfo);
     flowvelocity->resultType = MEAN_FLUIDMECH_VELOCITY;
     flowvelocity->dofNames = dofNames;
     flowvelocity->unit = "m/s";

     flowvelocity->definedOn = ResultInfo::NODE;
     flowvelocity->entryType = ResultInfo::VECTOR;



     shared_ptr<BaseFeFunction> meanFunction;
     std::string form = SolutionTypeEnum.ToString(MEAN_FLUIDMECH_VELOCITY);
     PtrParamNode feSpaceNode = infoNode_->Get("feSpaces");
     PtrParamNode potSpaceNode = feSpaceNode->Get(form);
     shared_ptr<FeSpace> tmpSpace = FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1);

     if(isComplex_){
       meanFunction.reset(new FeFunction<Complex>());
       meanFunction->SetFeSpace(tmpSpace);
       meanFunction->SetResultInfo(flowvelocity);
       meanFunction->SetGrid(ptgrid_);
       meanFunction->SetPDE(this);
       flowvelocity->SetFeFunction(meanFunction);
       if(dim_==2)
         meanFlowFunctor_.reset(new ExternalFieldFunctor<IdentityOperator<FeH1,2,2,Complex>,Complex >(meanFunction,flowvelocity));
       else
         meanFlowFunctor_.reset(new ExternalFieldFunctor<IdentityOperator<FeH1,3,3,Complex>,Complex >(meanFunction,flowvelocity));
     }else{
       meanFunction.reset(new FeFunction<Double>());
       meanFunction->SetFeSpace(tmpSpace);
       meanFunction->SetResultInfo(flowvelocity);
       meanFunction->SetGrid(ptgrid_);
       meanFunction->SetPDE(this);
       flowvelocity->SetFeFunction(meanFunction);
       if(dim_==2)
         meanFlowFunctor_.reset(new ExternalFieldFunctor<IdentityOperator<FeH1,2,2,Double>,Double >(meanFunction,flowvelocity));
       else
         meanFlowFunctor_.reset(new ExternalFieldFunctor<IdentityOperator<FeH1,3,3,Double>,Double >(meanFunction,flowvelocity));
     }

     results_.Push_back( flowvelocity );
     availResults_.insert( flowvelocity );

   }

   void AcousticMixedPDE::InitTimeStepping(){
     shared_ptr<BaseTimeScheme> mySchemeV(new TimeSchemeGLM(TimeSchemeGLM::TRAPEZOIDAL, 0) );
     shared_ptr<BaseTimeScheme> mySchemeP(new TimeSchemeGLM(TimeSchemeGLM::TRAPEZOIDAL, 0) );

     feFunctions_[ACOU_PRESSURE]->SetTimeScheme(mySchemeP);
     feFunctions_[ACOU_VELOCITY]->SetTimeScheme(mySchemeV);

   }
}

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template void AcousticMixedPDE::DefineRhsLoadIntegratorsTempl<Double,2>();
  template void AcousticMixedPDE::DefineRhsLoadIntegratorsTempl<Complex,2>();
  template void AcousticMixedPDE::DefineRhsLoadIntegratorsTempl<Double,3>();
  template void AcousticMixedPDE::DefineRhsLoadIntegratorsTempl<Complex,3>();

  template void AcousticMixedPDE::DefineIntegratorsTempl<Complex,2>();
  template void AcousticMixedPDE::DefineIntegratorsTempl<Double,2>();
  template void AcousticMixedPDE::DefineIntegratorsTempl<Complex,3>();
  template void AcousticMixedPDE::DefineIntegratorsTempl<Double,3>();
#endif
