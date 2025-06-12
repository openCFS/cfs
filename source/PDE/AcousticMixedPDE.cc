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
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/BiLinForms/ADBInt.hh"
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
#include "Domain/CoefFunction/CoefFunctionPML.hh"
#include "Domain/CoefFunction/CoefFunctionCompound.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"

#include "Domain/Results/ResultFunctor.hh"

#include "Driver/Assemble.hh"

#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Materials/AcousticMaterial.hh"

namespace CoupledField{

   DEFINE_LOG(acousticmixedpde, "pde.acousticmixed")

   AcousticMixedPDE::AcousticMixedPDE( Grid* aGrid, PtrParamNode paramNode,
                                       PtrParamNode infoNode,
                                       shared_ptr<SimState> simState, 
                                       Domain* domain)
               : SinglePDE( aGrid, paramNode, infoNode, simState, domain ){

     pdename_           = "acousticMixed";
     pdematerialclass_  = ACOUSTIC;
     nonLin_            = false;
     usePiola_          = false;
     penalized_         = true;
     doFluxTerm_        = true;
     
     //! Always use total Lagrangian formulation 
     updatedGeo_        = false;
     isTimeDomPML_      = false;
     isTaylorHood_      = false;
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
    }else if(formulation == "TaylorHood"){
      std::string form = SolutionTypeEnum.ToString(ACOU_PRESSURE);
      PtrParamNode potSpaceNode = infoNode->Get(form);
      crSpaces[ACOU_PRESSURE] =
              FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[ACOU_VELOCITY] =
          FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);

      crSpaces[ACOU_PRESSURE]->Init(solStrat_);
      crSpaces[ACOU_VELOCITY]->Init(solStrat_);
      isTaylorHood_ = true;
    }else{
      EXCEPTION("The formulation " << formulation << "of acousticMixed PDE is not known!");
    }
    usePiola_ = myParam_->Get("usePiolaTransform")->As<bool>();
    penalized_ = myParam_->Get("penalized")->As<bool>();
    doFluxTerm_ = myParam_->Get("fluxTerm")->As<bool>();

    if(isTimeDomPML_){
      PtrParamNode vectorPML = infoNode->Get("TransientPMLVectorAuxVar");
      crSpaces[ACOU_PMLAUXVEC] =
          FeSpace::CreateInstance(myParam_,vectorPML,FeSpace::H1, ptGrid_);
      crSpaces[ACOU_PMLAUXVEC]->Init(solStrat_);
    }
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
    isaxi_ = ptGrid_->IsAxi();

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    shared_ptr<FeSpace> spaceP = feFunctions_[ACOU_PRESSURE]->GetFeSpace();
    shared_ptr<FeSpace> spaceV = feFunctions_[ACOU_VELOCITY]->GetFeSpace();
    BaseMaterial * actMat = NULL;
    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion++){
      actRegion = regions_[iRegion];
      actMat    = materials_[actRegion];
      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
      
      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );

      // --- Set the FE ansatz for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());

      std::string polyIdV;
      std::string polyIdP;
      std::string integIdV;
      std::string integIdP;

      if(isTaylorHood_){
        polyIdP = "presP";
        polyIdV = "velP";

        integIdP = "presI";
        integIdV = "velI";
      }else{
        curRegNode->GetValue("polyId", polyIdP);
        curRegNode->GetValue("integId", integIdP);
        integIdV = integIdP;
        polyIdV = polyIdP;
      }
      spaceP->SetRegionApproximation(actRegion, polyIdP,integIdP);
      spaceV->SetRegionApproximation(actRegion, polyIdV,integIdV);


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
        stiffIntPV = new ABInt<DATA_TYPE>(new GradientOperator<FeH1,DIM,1,DATA_TYPE>(),
                                          new IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>() , 
                                          coeffKPV,-1.0, updatedGeo_ );
      else
        stiffIntPV = new ABInt<DATA_TYPE>(new  GradientOperator<FeH1,DIM,1,DATA_TYPE>() ,
                                          new IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>(),  
                                          coeffKPV,-1.0, updatedGeo_ );

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
                                  new GradientOperator<FeH1,DIM,1,DATA_TYPE>(), //
                                  coeffKVP,1.0, updatedGeo_ );
      else
        stiffIntVP = 
            new ABInt<DATA_TYPE>(new  IdentityOperator<FeH1,DIM,DIM,DATA_TYPE> , //
                                 new GradientOperator<FeH1,DIM,1,DATA_TYPE>(),  //
                                 coeffKVP,1.0, updatedGeo_ );

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


      //=====================================================================
      //check for PML
      //=====================================================================
      if( dampingList_[actRegion] == PML ) {
        std::string dampId;
        curRegNode->GetValue("dampingId",dampId);
        if(analysistype_ == HARMONIC){
          EXCEPTION("Harmonic PML not implemented yet")
        }else{
          if(dim_==2)
            DefineTransientPMLInts<2>(actSDList,dampId);
          else
            DefineTransientPMLInts<3>(actSDList,dampId);
        }
      }

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

        //now create the integrators
        BiLinearForm *convectiveVV = NULL;
        BiLinearForm *convectivePP = NULL;

        PtrCoefFct convFactor;
        Double factorPPT = 1.0;
        if(doFluxTerm_) {
          convFactor = 
              CoefFunction::Generate(mp_, Global::REAL, 
                                     CoefXprBinOp(mp_, "0.5", density, CoefXpr::OP_MULT ) );
          factorPPT = 0.5;
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
                                 coeffPP, factorPPT , coefUpdateGeo);



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

          shared_ptr<ResultInfo> flowvelocityNormal( new ResultInfo);
          flowvelocityNormal->resultType = MEAN_FLUIDMECH_VELOCITY_NORMAL;
          flowvelocityNormal->dofNames = "";
          flowvelocityNormal->unit = "m/s";

          flowvelocityNormal->definedOn = ResultInfo::NODE;
          flowvelocityNormal->entryType = ResultInfo::SCALAR;


          shared_ptr<CoefFunctionSurf> normVel = 
              shared_ptr<CoefFunctionSurf>(new CoefFunctionSurf(true, 1.0, 
                                                                flowvelocityNormal));

          normVel->AddVolumeCoef(actRegion,meanFlowCoef_);

          PtrCoefFct formFactor2 =
              CoefFunction::Generate(mp_, Global::REAL,
                                     CoefXprBinOp(mp_,normVel, formFactor,
                                          CoefXpr::OP_MULT ) );
          PtrCoefFct formFactor3 =
              CoefFunction::Generate(mp_, Global::REAL,
                  CoefXprUnaryOp(mp_,formFactor2,CoefXpr::OP_NORM ) );

          BiLinearForm *convectiveVOpp = NULL;
          BiLinearForm *convectiveV = NULL;
          BiLinearForm *exteriorVV = NULL;
          std::set<RegionIdType> volRegion;
          volRegion.insert(actRegion);
          if( usePiola_ ) {
            convectiveVOpp = new SurfaceBBInt<DATA_TYPE>(new IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>(),
                formFactor3, -0.5,volRegion, updatedGeo_);
            convectiveV   = new SurfaceBBInt<DATA_TYPE>(new IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>(),
                formFactor3, 0.5, volRegion, updatedGeo_);
            exteriorVV     = new SurfaceBBInt<DATA_TYPE>(new IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>(),
                formFactor3, 0.5,volRegion, updatedGeo_);
          } else {
            convectiveVOpp = new BBInt<DATA_TYPE>(new IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>(),
                formFactor3, -0.5, updatedGeo_ );
            convectiveV    = new BBInt<DATA_TYPE>(new IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>(),
                formFactor3, 0.5, updatedGeo_);
            exteriorVV    = new BBInt<DATA_TYPE>(new IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>(),
                formFactor3, 0.5, updatedGeo_);
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

          penaltyContextVOpp->SetEntities(  oppositList,list );
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
          BiLinearForm *convectivePPT = NULL;

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

          convectivePPT = new ABInt<DATA_TYPE>(new ConvectiveOperator<FeH1,DIM,1,DATA_TYPE>(),
                                                 new IdentityOperator<FeH1,DIM,1,DATA_TYPE>(),
                                                 coeffPP, -0.5 , coefUpdateGeo);

          convectiveVVTrans->SetBCoefFunctionOpA(meanFlowCoef_);
          fluxTerm->SetBCoefFunctionOpB(meanFlowCoef_);

          fluxTerm->SetName("fluxTerm");
          convectiveVVTrans->SetName("convectiveVVtrans");

          convectivePPT->SetBCoefFunctionOpA(meanFlowCoef_);
          convectivePPT->SetName("convectivePPTransposed");

          //now we obtain the entity lists
          shared_ptr<EntityList> list,oppositList,extList;
          feFunctions_[ACOU_VELOCITY]->GetFeSpace()->GetInteriorSurfaceElems(actRegion,list,oppositList);
          feFunctions_[ACOU_VELOCITY]->GetFeSpace()->GetExteriorSurfaceElems(actRegion,extList);

          BiLinFormContext * fluxContext =  new BiLinFormContext(fluxTerm, STIFFNESS );
          BiLinFormContext * convectiveTransContext =  new BiLinFormContext(convectiveVVTrans, STIFFNESS );
          BiLinFormContext * convectiveContextPPT =  new BiLinFormContext(convectivePPT, STIFFNESS );

          fluxContext->SetEntities( oppositList,list  );
          convectiveTransContext->SetEntities( actSDList, actSDList );

          fluxContext->SetFeFunctions( feFunctions_[ACOU_VELOCITY],feFunctions_[ACOU_VELOCITY]);
          convectiveTransContext->SetFeFunctions( feFunctions_[ACOU_VELOCITY],feFunctions_[ACOU_VELOCITY]);
          convectiveContextPPT->SetEntities( actSDList, actSDList );
          convectiveContextPPT->SetFeFunctions( feFunctions_[ACOU_PRESSURE],feFunctions_[ACOU_PRESSURE]);
          assemble_->AddBiLinearForm( fluxContext );
          assemble_->AddBiLinearForm( convectiveTransContext );
          assemble_->AddBiLinearForm( convectiveContextPPT );
        }
      }
    }
  }

  template<UInt DIM>
  void AcousticMixedPDE::DefineTransientPMLInts(shared_ptr<ElemList> eList, std::string id){

    //define some material coeffunction as above...
    PtrCoefFct factor = CoefFunction::Generate( mp_, Global::REAL, "1.0");
    PtrCoefFct dens = materials_[eList->GetRegion()]->GetScalCoefFnc( DENSITY, Global::REAL );
    PtrCoefFct blk = materials_[eList->GetRegion()]->GetScalCoefFnc( ACOU_BULK_MODULUS, Global::REAL );
    // c0 = sqrt(bulk_modulus / density)
    PtrCoefFct c0 =
        CoefFunction::Generate( mp_,  Global::REAL,
                                CoefXprUnaryOp( mp_, CoefXprBinOp(mp_, blk, dens, CoefXpr::OP_DIV),
                                CoefXpr::OP_SQRT) );


    PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal("pml","id",id.c_str());
    shared_ptr<CoefFunction> coeffPMLVec;
    coeffPMLVec.reset( new CoefFunctionPML<Double>(pmlNode,c0,eList,regions_,true) );

    matCoefs_[PML_DAMP_FACTOR]->AddRegion(eList->GetRegion(), coeffPMLVec);

    //the tensorial PML
    shared_ptr<CoefFunctionCompound<Double> > coefA(new CoefFunctionCompound<Double>(mp_));
    //vector PML
    shared_ptr<CoefFunctionCompound<Double> > coefB(new CoefFunctionCompound<Double>(mp_));
    shared_ptr<CoefFunctionCompound<Double> > coefC(new CoefFunctionCompound<Double>(mp_));

    // --- Set the FE ansatz for the current region ---
    PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",eList->GetName().c_str());
    shared_ptr<FeSpace> vecSpace = feFunctions_[ACOU_PMLAUXVEC]->GetFeSpace();
    if(isTaylorHood_){
      vecSpace->SetRegionApproximation(eList->GetRegion(), "velP" ,"velI");
    }else{
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      vecSpace->SetRegionApproximation(eList->GetRegion(), polyId ,integId);
    }



    std::map<std::string, PtrCoefFct> vars;
    std::map<std::string, PtrCoefFct> var;
    var["a"]  = coeffPMLVec;
    vars["a"]  = coeffPMLVec;
    vars["b"] = dens;

    StdVector<std::string> matAReal;
    StdVector<std::string> matBReal;
    StdVector<std::string> matCReal;

    if(DIM == 3 ){
      const std::string Amat[] =  { "a_0_R", "0.0" , "0.0" , "0.0" , "a_1_R" , "0.0" , "0.0" , "0.0" , "a_2_R" };
      const std::string Bmat[] =  { "a_0_R", "a_1_R" , "a_2_R" };
      const std::string Cmat[] =  { "a_0_R * b_R", "0.0" , "0.0" , "0.0" , "a_1_R * b_R" , "0.0" , "0.0" , "0.0" , "a_2_R * b_R" };
      matAReal.Import(Amat,9);
      matBReal.Import(Bmat,3);
      matCReal.Import(Cmat,9);
    }else{
      const std::string Amat[] =  { "a_0_R", "0.0" , "0.0" , "a_1_R" };
      const std::string Bmat[] = {"a_0_R", "a_1_R"};
      const std::string Cmat[] =  { "a_0_R * b_R", "0.0" , "0.0" , "a_1_R * b_R" };
      matAReal.Import(Amat,4);
      matBReal.Import(Bmat,2);
      matCReal.Import(Cmat,4);
    }
    coefA->SetTensor(matAReal,DIM,DIM,var);
    coefB->SetTensor(matBReal,1,DIM,var);
    coefC->SetTensor(matCReal,DIM,DIM,vars);

    BaseBDBInt *  dampQP = new ADBInt<>(new IdentityOperator<FeH1,DIM,1>(), new IdentityOperator<FeH1,DIM,DIM>(), coefB, 1.0, updatedGeo_ );
    BaseBDBInt *  dampVV = new BDBInt<>(new IdentityOperator<FeH1,DIM,DIM>(), coefC, 1.0, updatedGeo_ );
    BaseBDBInt *  timeQQ = new BBInt<>(new IdentityOperator<FeH1,DIM,DIM>(), factor , 1.0 , updatedGeo_ );
    BaseBDBInt *  dampQQ = new BDBInt<>(new IdentityOperator<FeH1,DIM,DIM>(),coefA, 1.0 , updatedGeo_ );
    BaseBDBInt *  dampQV = new ABInt<>(new IdentityOperator<FeH1,DIM,DIM>(), new GradientOperator<FeH1,DIM,DIM>(), factor , 1.0, updatedGeo_ );

    dampQP->SetName("dampQP");
    dampVV->SetName("dampVV");
    timeQQ->SetName("timeQQ");
    dampQQ->SetName("dampQQ");
    dampQV->SetName("dampQV");

    BiLinFormContext *contextDampQP = new BiLinFormContext(dampQP,STIFFNESS);
    BiLinFormContext *contextDampVV = new BiLinFormContext(dampVV,STIFFNESS);
    BiLinFormContext *contextTimeQQ = new BiLinFormContext(timeQQ,DAMPING);
    BiLinFormContext *contextDampQQ = new BiLinFormContext(dampQQ,STIFFNESS);
    BiLinFormContext *contextDampQV = new BiLinFormContext(dampQV,STIFFNESS);

    contextDampQP->SetEntities( eList, eList );
    contextDampVV->SetEntities( eList, eList );
    contextTimeQQ->SetEntities( eList, eList );
    contextDampQQ->SetEntities( eList, eList );
    contextDampQV->SetEntities( eList, eList );

    contextDampQP->SetFeFunctions(feFunctions_[ACOU_PRESSURE],feFunctions_[ACOU_PMLAUXVEC]);
    contextDampVV->SetFeFunctions(feFunctions_[ACOU_VELOCITY],feFunctions_[ACOU_VELOCITY]);
    contextTimeQQ->SetFeFunctions(feFunctions_[ACOU_PMLAUXVEC],feFunctions_[ACOU_PMLAUXVEC]);
    contextDampQQ->SetFeFunctions(feFunctions_[ACOU_PMLAUXVEC],feFunctions_[ACOU_PMLAUXVEC]);
    contextDampQV->SetFeFunctions(feFunctions_[ACOU_PMLAUXVEC],feFunctions_[ACOU_VELOCITY]);

    assemble_->AddBiLinearForm(contextDampQP);
    assemble_->AddBiLinearForm(contextDampVV);
    assemble_->AddBiLinearForm(contextDampQQ);
    assemble_->AddBiLinearForm(contextTimeQQ);
    assemble_->AddBiLinearForm(contextDampQV);

    feFunctions_[ACOU_PMLAUXVEC]->AddEntityList( eList );

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

      if(coef[i]->IsConservative()){
        this->rhsFeFunctions_[ACOU_PRESSURE]->AddLoadCoefFunction(coef[i], ent[i]);
      }else{
        lin = new BUIntegrator<DATA_TYPE>( new IdentityOperator<FeH1,DIM,1,DATA_TYPE>(),
                                           1.0,coef[i],coefUpdateGeo);

        lin->SetName("massEquationInt");
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( ent[i] );
        ctx->SetFeFunction(pressureFct);
        assemble_->AddLinearForm(ctx);
      }
    }

    LOG_DBG(acousticmixedpde) << "Reading loads for momentum conservation equation";
    ReadRhsExcitation( "momentumEquationLoad", vDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      if(coef[i]->IsConservative()){

         //check for PIOLA crap... we do not support this here right now
         if(usePiola_)
           EXCEPTION("Piola mapping not supported for conservative interpolation...");

         this->rhsFeFunctions_[ACOU_VELOCITY]->AddLoadCoefFunction(coef[i], ent[i]);
      }else{
        if(usePiola_){
          lin = new BUIntegrator<DATA_TYPE>( new IdentityOperatorPiola<FeH1,DIM,DIM,DATA_TYPE>(),
                                             1.0,coef[i],coefUpdateGeo);
        }else{
          lin = new BUIntegrator<DATA_TYPE>( new IdentityOperator<FeH1,DIM,DIM,DATA_TYPE>(),
                                             1.0,coef[i],coefUpdateGeo);
        }

        lin->SetName("momentumEquationInt");
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( ent[i] );
        ctx->SetFeFunction(velocityFct);
        assemble_->AddLinearForm(ctx);
      }
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

      if( ptGrid_->GetDim() == 3 ) {
        velDofNames = "x", "y", "z";
      } else {
        if( ptGrid_->IsAxi() ) {
          velDofNames = "r", "z";
        } else {
          velDofNames = "x", "y";
        }
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
      
      // === ACOUSTIC MASS RHS ===
      shared_ptr<ResultInfo> rhsP ( new ResultInfo );
      rhsP->resultType = ACOU_MIXED_MASS_LOAD;
      rhsP->dofNames = "";
      rhsP->unit = "?";
      rhsP->definedOn = ResultInfo::NODE;
      rhsP->entryType = ResultInfo::SCALAR;
      rhsP->SetFeFunction(feFunctions_[ACOU_PRESSURE]);
      DefineFieldResult( this->rhsFeFunctions_[ACOU_PRESSURE], rhsP );
      results_.Push_back( rhsP );
      availResults_.insert( rhsP );

      // === ACOUSTIC MOMENTUM RHS ===
      shared_ptr<ResultInfo> rhsV ( new ResultInfo );
      rhsV->resultType = ACOU_MIXED_MOMENTUM_LOAD;
      rhsV->dofNames = velDofNames;
      rhsV->unit = "?";
      rhsV->definedOn = ResultInfo::NODE;
      rhsV->entryType = ResultInfo::VECTOR;;
      rhsV->SetFeFunction(feFunctions_[ACOU_VELOCITY]);
      DefineFieldResult( this->rhsFeFunctions_[ACOU_VELOCITY], rhsV );
      results_.Push_back( rhsV );
      availResults_.insert( rhsV );


      // === PML DAMPING FACTORS ===
      //if( matCoefs_.find(PML_DAMP_FACTOR) != matCoefs_.end() ) {
      shared_ptr<ResultInfo> pml ( new ResultInfo );
      pml->resultType = PML_DAMP_FACTOR;
      pml->dofNames = velDofNames;
      //pml->dofNames = "";
      pml->unit = "";
      pml->definedOn = ResultInfo::ELEMENT;
      pml->entryType = ResultInfo::VECTOR;
      pml->SetFeFunction(feFunctions_[ACOU_VELOCITY]);
      shared_ptr<CoefFunctionMulti> pmlFct(new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1,
                                                                 isComplex_));
      matCoefs_[PML_DAMP_FACTOR] = pmlFct;
      DefineFieldResult(pmlFct, pml);

      // === PML AUX Variables ===
      if(this->isTimeDomPML_){
        shared_ptr<ResultInfo> pmlVec ( new ResultInfo );
        pmlVec->resultType = ACOU_PMLAUXVEC;
        pmlVec->dofNames = velDofNames;
        pmlVec->unit = "-";
        pmlVec->definedOn = ResultInfo::NODE;
        pmlVec->entryType = ResultInfo::VECTOR;
        feFunctions_[ACOU_PMLAUXVEC]->SetResultInfo(pmlVec);
        results_.Push_back( pmlVec );
        pmlVec->SetFeFunction(feFunctions_[ACOU_PMLAUXVEC]);
        DefineFieldResult( feFunctions_[ACOU_PMLAUXVEC], pmlVec );
      }

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
     flowvelocity->SetFeFunction(feFunctions_[ACOU_VELOCITY]);
     results_.Push_back( flowvelocity );
     availResults_.insert( flowvelocity );
     meanFlowCoef_.reset(new CoefFunctionMulti(CoefFunction::VECTOR,dim_,1, 
                                               isComplex_));
     
     DefineFieldResult( meanFlowCoef_, flowvelocity );
   }

   void AcousticMixedPDE::InitTimeStepping(){
     shared_ptr<BaseTimeScheme> mySchemeV(new TimeSchemeGLM(GLMScheme::BDF2, 0) );
     shared_ptr<BaseTimeScheme> mySchemeP(new TimeSchemeGLM(GLMScheme::BDF2, 0) );

     feFunctions_[ACOU_PRESSURE]->SetTimeScheme(mySchemeP);
     feFunctions_[ACOU_VELOCITY]->SetTimeScheme(mySchemeV);

     if(this->isTimeDomPML_){
       shared_ptr<BaseTimeScheme> mySchemeQ(new TimeSchemeGLM(GLMScheme::BDF2, 0) );
       feFunctions_[ACOU_PMLAUXVEC]->SetTimeScheme(mySchemeQ);
     }
   }

   void AcousticMixedPDE::ReadDampingInformation() {
     std::map<std::string, DampingType> idDampType;

     // try to get dampingList
     PtrParamNode dampListNode = myParam_->Get( "dampingList", ParamNode::PASS );
     if( dampListNode ) {

       // get specific damping nodes
       ParamNodeList dampNodes = dampListNode->GetChildren();

       for( UInt i = 0; i < dampNodes.GetSize(); i++ ) {

         std::string dampString = dampNodes[i]->GetName();
         std::string actId = dampNodes[i]->Get("id")->As<std::string>();

         // determine type of damping
         DampingType actType;
         String2Enum( dampString, actType );

         // store damping type string
         idDampType[actId] = actType;
       }
     }

     // Run over all region and set entry in "regionNonLinId"
     ParamNodeList regionNodes =
         myParam_->Get("regionList")->GetChildren();

     RegionIdType actRegionId;
     std::string actRegionName, actDampingId;

     for (UInt k = 0; k < regionNodes.GetSize(); k++) {
       regionNodes[k]->GetValue( "name", actRegionName );
       regionNodes[k]->GetValue( "dampingId", actDampingId );
       if( actDampingId == "" )
         continue;

       actRegionId = ptGrid_->GetRegion().Parse( actRegionName );

       // Check actDampingId was already registerd
       if( idDampType.count( actDampingId ) == 0 ) {
         EXCEPTION( "Damping with id '" << actDampingId
                    << "' was not defined in 'dampingList'" );
       }

       dampingList_[actRegionId] = idDampType[actDampingId];
       if(dampingList_[actRegionId] == PML &&
          analysistype_ == BasePDE::TRANSIENT ) {
         isTimeDomPML_ = true;
       }
     }
   }
}

template void AcousticMixedPDE::DefineRhsLoadIntegratorsTempl<Double,2>();
template void AcousticMixedPDE::DefineRhsLoadIntegratorsTempl<Complex,2>();
template void AcousticMixedPDE::DefineRhsLoadIntegratorsTempl<Double,3>();
template void AcousticMixedPDE::DefineRhsLoadIntegratorsTempl<Complex,3>();

template void AcousticMixedPDE::DefineIntegratorsTempl<Complex,2>();
template void AcousticMixedPDE::DefineIntegratorsTempl<Double,2>();
template void AcousticMixedPDE::DefineIntegratorsTempl<Complex,3>();
template void AcousticMixedPDE::DefineIntegratorsTempl<Double,3>();

template void AcousticMixedPDE::DefineTransientPMLInts<2>(shared_ptr<ElemList>, std::string);
template void AcousticMixedPDE::DefineTransientPMLInts<3>(shared_ptr<ElemList>, std::string);
