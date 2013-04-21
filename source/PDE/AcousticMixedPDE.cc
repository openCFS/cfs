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

#include "FeBasis/H1/FeSpaceH1Nodal.hh"
#include "FeBasis/L2/FeSpaceL2.hh"
#include "FeBasis/L2/FeSpaceL2Nodal.hh"


#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/Results/ResultFunctor.hh"

#include "Driver/Assemble.hh"

#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Materials/AcousticMaterial.hh"

namespace CoupledField{

  DECLARE_LOG(acousticmixedpde)
   DEFINE_LOG(acousticmixedpde, "pde.acousticmixed")

   AcousticMixedPDE::AcousticMixedPDE( Grid* aGrid, PtrParamNode paramNode,
                                       PtrParamNode infoNode,
                                       shared_ptr<SimState> simState, 
                                       Domain* domain)
               : SinglePDE( aGrid, paramNode, infoNode, simState, domain ){

     pdename_           = "acousticMixed";
     pdematerialclass_  = FLUID;
     nonLin_            = false;
     usePiola_          = false;
     penalized_         = true;
     doFluxTerm_        = true;
     
     //! Always use total Lagrangian formulation 
     updatedGeo_        = false;

   }

  std::map<SolutionType, shared_ptr<FeSpace> >
  AcousticMixedPDE::CreateFeSpaces( const std::string&  formulation,
                    PtrParamNode infoNode ){

    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "MixedPiola" || formulation == "Mixed"){
      std::string form = SolutionTypeEnum.ToString(ACOU_PRESSURE);
      PtrParamNode potSpaceNode = infoNode->Get(form);
      crSpaces[ACOU_PRESSURE] =
              FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[ACOU_VELOCITY] =
          FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::L2, ptGrid_);

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
    domain_->GetParamRoot()->Get("domain")->GetValue("geometryType", geometryType );

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
    BaseMaterial * actMat = NULL;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      // Set current region and material
      actRegion = it->first;
      actMat = it->second;
      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
      
      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );

      // --- Set the FE ansatz for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId;
      curRegNode->GetValue("polyId", polyId);
      std::string integId;
      curRegNode->GetValue("integId", integId);
      spaceP->SetRegionApproximation(actRegion, polyId,integId);
      spaceV->SetRegionApproximation(actRegion, polyId,integId);


      // Obtain density and compressibility as coefficient functions
      PtrCoefFct density = actMat->GetScalCoefFnc( DENSITY, Global::REAL );
      PtrCoefFct compressibility = actMat->GetScalCoefFnc( ACOU_BULK_MODULUS, Global::REAL );

      feFunctions_[ACOU_PRESSURE]->AddEntityList( actSDList );
      feFunctions_[ACOU_VELOCITY]->AddEntityList( actSDList );


      // ====================================================================
      // stiffness integrators
      // ====================================================================
      // --------------------------------------------------------------------
      //  VERSION 1: K_PV Integrator (upper off-diagonal integrator)
      // --------------------------------------------------------------------

      PtrCoefFct coeffKPV = CoefFunction::Generate(mp_, Global::REAL, "1.0");
      BiLinearForm * stiffIntPV = NULL;

      if(usePiola_)
        stiffIntPV = new ABInt<DATA_TYPE>(new GradientOperator<FeH1,DIM,DATA_TYPE>(),
                                          new IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>() , 
                                          coeffKPV,1.0, updatedGeo_ );
      else
        stiffIntPV = new ABInt<DATA_TYPE>(new  GradientOperator<FeH1,DIM,DATA_TYPE>() ,
                                          new IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>(),  
                                          coeffKPV,1.0, updatedGeo_ );

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
      PtrCoefFct coeffKVP = CoefFunction::Generate(mp_, Global::REAL, "1.0");
      BiLinearForm * stiffIntVP = NULL;

      if(usePiola_)
        stiffIntVP = 
            new ABInt<DATA_TYPE >(new IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>() , //
                                  new GradientOperator<FeH1,DIM,DATA_TYPE>(), //
                                  coeffKVP,-1.0, updatedGeo_ );
      else
        stiffIntVP = 
            new ABInt<DATA_TYPE>(new  IdentityOperator<FeH1,DIM,DIM,DATA_TYPE> , //
                                 new GradientOperator<FeH1,DIM,DATA_TYPE>(),  //
                                 coeffKVP,-1.0, updatedGeo_ );

      stiffIntVP->SetName("MixedStiffIntVP");
      BiLinFormContext *stiffContVP = new BiLinFormContext(stiffIntVP, STIFFNESS );

      stiffContVP->SetEntities( actSDList, actSDList );
      stiffContVP->SetFeFunctions( feFunctions_[ACOU_VELOCITY],feFunctions_[ACOU_PRESSURE]);
      assemble_->AddBiLinearForm( stiffContVP );

      // ====================================================================
      // mass integrators
      // ====================================================================
      // coefficent for mass integrator: 1.0 / compressibility
      PtrCoefFct coeffMPP
          = CoefFunction::Generate(mp_, Global::REAL, 
                                   CoefXprBinOp(mp_, "1.0", compressibility, CoefXpr::OP_DIV) );

      BiLinearForm *massIntPP = NULL;

      massIntPP = new BBInt<DATA_TYPE>(new IdentityOperator<FeH1,DIM,1,DATA_TYPE>() , 
                                       coeffMPP, 1.0, updatedGeo_);

      massIntPP->SetName("PressureTimeInt");
      //massIntPP->SetFeSpace( feFunctions_[ACOU_PRESSURE]->GetFeSpace() );

      BiLinFormContext *massContextPP =  new BiLinFormContext(massIntPP, DAMPING );

      massContextPP->SetEntities( actSDList, actSDList );
      massContextPP->SetFeFunctions( feFunctions_[ACOU_PRESSURE],feFunctions_[ACOU_PRESSURE]);
      assemble_->AddBiLinearForm( massContextPP );

      BiLinearForm *massIntVV = NULL;
      if( usePiola_)
        massIntVV = new BBInt<DATA_TYPE>(new IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>(), 
                                         density, 1.0, updatedGeo_);
      else
        massIntVV = new BBInt<DATA_TYPE>(new IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>(),
                                         density, 1.0, updatedGeo_);

      massIntVV->SetName("VelocityTimeInt");
      //massIntVV->SetFeSpace( feFunctions_[ACOU_VELOCITY]->GetFeSpace() );

      BiLinFormContext *massContextVV =  new BiLinFormContext(massIntVV, DAMPING );

      massContextVV->SetEntities( actSDList, actSDList );
      massContextVV->SetFeFunctions( feFunctions_[ACOU_VELOCITY],feFunctions_[ACOU_VELOCITY]);
      assemble_->AddBiLinearForm( massContextVV );

      //======================================================================
      // CHECK FOR FLOW
      //=====================================================================
      std::string flowId;
      curRegNode->GetValue("flowId", flowId);
      if(flowId != ""){
        
        // Get result info object for flow
        shared_ptr<ResultInfo> flowInfo = GetResultInfo(MEAN_FLUIDMECH_VELOCITY);
        
        
        //Add the region information
        PtrParamNode flowNode = myParam_->Get("flowList")->GetByVal("flow","name",flowId.c_str());
        
        // Read coefficient flow coefficient function for this region and add it to flow functor
        PtrCoefFct regionFlow;
        std::set<UInt> definedDofs;
        bool coefUpdateGeo;
        ReadUserFieldValues( actSDList, flowNode, flowInfo->dofNames, flowInfo->entryType, 
                             isComplex_, regionFlow, definedDofs, coefUpdateGeo );
        meanFlowCoef_->AddRegion( actRegion, regionFlow );
        regionFlow->AddEntityList(actSDList);

        //now create the integrators
        BiLinearForm *convectiveVV = NULL;
        BiLinearForm *convectivePP = NULL;
        PtrCoefFct convFactor;
        
        if(doFluxTerm_) {
          convFactor = 
              CoefFunction::Generate(mp_, Global::REAL, 
                                     CoefXprBinOp(mp_, "0.5", density, CoefXpr::OP_MULT ) );
        } else {
          convFactor = density;
        }

        if(usePiola_){
          convectiveVV = 
              new ABInt<DATA_TYPE>( new IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>(), 
                                    new ConvectiveOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>(), 
                                    convFactor, 1.0, coefUpdateGeo);
        } else {
          convectiveVV = 
              new ABInt<DATA_TYPE>(new IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>(), 
                                   new ConvectiveOperator<FeH1,DIM,DIM,DATA_TYPE>(), 
                                   convFactor, 1.0, coefUpdateGeo);
        }
        
        // Coefficient for mass integrator: 1.0 / compressibility
        PtrCoefFct coeffPP
        = CoefFunction::Generate(mp_, Global::REAL, 
                                 CoefXprBinOp(mp_, "1.0", compressibility, CoefXpr::OP_DIV) );
        
        convectivePP = 
            new ABInt<DATA_TYPE>(new IdentityOperator<FeH1,DIM,1,DATA_TYPE>(), 
                                 new ConvectiveOperator<FeH1,DIM,1,DATA_TYPE>(),
                                 coeffPP, 1.0 , coefUpdateGeo);

        convectiveVV->SetBCoefFunctionOpB(meanFlowCoef_);
        convectivePP->SetBCoefFunctionOpB(meanFlowCoef_);

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
          Double penaltyFactor = 0.0;
          curRegNode->GetValue("penalizationFactor", penaltyFactor);
          
          // formFactor = density * penaltyFactor
          PtrCoefFct formFactor = 
              CoefFunction::Generate(mp_, Global::REAL,
                                     CoefXprBinOp(mp_, lexical_cast<std::string>(penaltyFactor),
                                                  density, CoefXpr::OP_MULT ) );
          BiLinearForm *convectiveVOpp = NULL;
          BiLinearForm *convectiveV = NULL;
          BiLinearForm *exteriorVV = NULL;
          std::set<RegionIdType> volRegion;
          volRegion.insert(actRegion);
          if( usePiola_ ) {
            convectiveVOpp = new SurfaceBBInt<DATA_TYPE>(new IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>(),
                                                         formFactor, -0.5,volRegion, updatedGeo_);
            convectiveV    = new SurfaceBBInt<DATA_TYPE>(new IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>(),
                                                         formFactor, 0.5, volRegion, updatedGeo_);
            exteriorVV     = new SurfaceBBInt<DATA_TYPE>(new IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>(),
                                                         formFactor, -0.5,volRegion, updatedGeo_);
          } else {
            convectiveVOpp = new BBInt<DATA_TYPE>(new IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>(),
                                                  formFactor, -0.5, updatedGeo_ );
            convectiveV    = new BBInt<DATA_TYPE>(new IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>(), 
                                                  formFactor, 0.5, updatedGeo_);
            exteriorVV     = new BBInt<DATA_TYPE>(new IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>(), 
                                                  formFactor, -0.5, updatedGeo_);
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
            fluxTerm = new SurfaceABInt<DATA_TYPE>(new IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>(),//
                                                   new NormalFluxOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>(),//
                                                   density, -0.5 ,volRegion, updatedGeo_);

            convectiveVVTrans = new ABInt<DATA_TYPE>(new  ConvectiveOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>(),//
                                                     new IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>(),//
                                                     density, -0.5, updatedGeo_);
          } else {
            fluxTerm = new SurfaceABInt<DATA_TYPE>(new IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>(),//
                                                   new NormalFluxOperator<FeH1,DIM,DIM,DATA_TYPE>(),//
                                                   density, -0.5,volRegion, updatedGeo_);

            convectiveVVTrans = new ABInt<DATA_TYPE>(new  ConvectiveOperator<FeH1,DIM,DIM,DATA_TYPE>(), //
                                                     new IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>(), //
                                                     density, -0.5, updatedGeo_);
          }

          convectiveVVTrans->SetBCoefFunctionOpA(meanFlowCoef_);
          fluxTerm->SetBCoefFunctionOpB(meanFlowCoef_);

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
    StdVector<PtrCoefFct > coef;
    LinearForm * lin = NULL;
    StdVector<std::string> vDofNames = velocityFct->GetResultInfo()->dofNames;

    bool coefUpdateGeo;
    LOG_DBG(acousticmixedpde) << "Reading loads for mass conservation equation";
    ReadRhsExcitation( "massEquationLoad", pressureFct->GetResultInfo()->dofNames, 
                       ResultInfo::SCALAR, isComplex_, ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {

      lin = new BUIntegrator<IdentityOperator<FeH1,DIM,1,DATA_TYPE>, DATA_TYPE >(1.0,coef[i],coefUpdateGeo);

      lin->SetName("massEquationInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(pressureFct);
      assemble_->AddLinearForm(ctx);
    }

    LOG_DBG(acousticmixedpde) << "Reading loads for momentum conservation equation";
    ReadRhsExcitation( "momentumEquationLoad", vDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      if(usePiola_){
        lin = new BUIntegrator<IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>, DATA_TYPE >(1.0,coef[i],coefUpdateGeo);
      }else{
        lin = new BUIntegrator<IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>, DATA_TYPE >(1.0,coef[i],coefUpdateGeo);
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
        shared_ptr<EntityList> actSDList =  ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );
        std::string volRegName = abcNodes[i]->Get("volumeRegion")->As<std::string>();
        RegionIdType sRegId = ptGrid_->GetRegion().Parse(regionName);

        // --- Set the FE ansatz for the current region ---
        PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",volRegName.c_str());
        std::string polyId = curRegNode->Get("polyId")->As<std::string>();
        std::string integId = curRegNode->Get("integId")->As<std::string>();
        feFunctions_[ACOU_PRESSURE]->GetFeSpace()->SetRegionApproximation(sRegId, polyId,integId);

        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);
        //
        PtrCoefFct density = materials_[aRegion]->GetScalCoefFnc( DENSITY, Global::REAL );
        PtrCoefFct compressibility = materials_[aRegion]->GetScalCoefFnc( ACOU_BULK_MODULUS, Global::REAL );
        
        // c0 = sqrt(bulk_modulus / density)
        PtrCoefFct c0 = 
            CoefFunction::Generate(mp_,  Global::REAL,
                                    CoefXprUnaryOp( mp_, CoefXprBinOp(mp_, compressibility, density, CoefXpr::OP_DIV),
                                    CoefXpr::OP_SQRT) );

        PtrCoefFct coeffKPV
        = CoefFunction::Generate(mp_, Global::REAL,
                                 CoefXprBinOp(mp_, "1.0", 
                                              CoefXprBinOp(mp_, density, c0, CoefXpr::OP_MULT ),
                                              CoefXpr::OP_DIV) ) ;
              
        BiLinearForm * abcInt = NULL;

        if( dim_ == 2 ) {
          if(isComplex_)
            abcInt = new BBInt<Complex>(new IdentityOperator<FeH1,2,1,Complex>(), coeffKPV,1.0, updatedGeo_ );
          else
            abcInt = new BBInt<Double>(new IdentityOperator<FeH1,2,1,Double>(), coeffKPV,1.0, updatedGeo_ );
        } else {
          if(isComplex_)
            abcInt = new BBInt<Complex>(new IdentityOperator<FeH1,3,1,Complex>(), coeffKPV,1.0, updatedGeo_ );
          else
            abcInt = new BBInt<Double>(new IdentityOperator<FeH1,3,1,Double>(), coeffKPV,1.0, updatedGeo_ );
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

    //!  Define available postprocessing results
    void AcousticMixedPDE::DefinePrimaryResults(){

      // Check for subType
      StdVector<std::string> velDofNames;

      std::string geometryType;
      domain_->GetParamRoot()->Get("domain")->GetValue("geometryType", geometryType );

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
      DefineFieldResult( feFunctions_[ACOU_PRESSURE], pressure );

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
      DefineFieldResult( feFunctions_[ACOU_VELOCITY], velocity );

      //  Define xml-names of Dirichlet BCs
      // -----------------------------------
      hdbcSolNameMap_[ACOU_PRESSURE] = "soundSoft";
      hdbcSolNameMap_[ACOU_VELOCITY] = "soundHard";
      idbcSolNameMap_[ACOU_PRESSURE] = "pressure";
      idbcSolNameMap_[ACOU_VELOCITY] = "velocity";
      
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
     results_.Push_back( flowvelocity );
     availResults_.insert( flowvelocity );
     meanFlowCoef_.reset(new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1, 
                                               isComplex_));
     
     DefineFieldResult( meanFlowCoef_, flowvelocity );
   }

   void AcousticMixedPDE::InitTimeStepping(){
     shared_ptr<BaseTimeScheme> mySchemeV(new TimeSchemeGLM(GLMScheme::TRAPEZOIDAL, 0) );
     shared_ptr<BaseTimeScheme> mySchemeP(new TimeSchemeGLM(GLMScheme::TRAPEZOIDAL, 0) );

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
