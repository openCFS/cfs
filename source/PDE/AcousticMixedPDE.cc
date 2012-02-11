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
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"

#include "FeBasis/FeFunctions.hh"

#include "FeBasis/H1/FeSpaceH1.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"
#include "FeBasis/L2/FeSpaceL2.hh"
#include "FeBasis/L2/FeSpaceL2Nodal.hh"


#include "Domain/Results/ResultFunctor.hh"

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

   }

  std::map<SolutionType, shared_ptr<FeSpace> >
  AcousticMixedPDE::CreateFeSpaces( const std::string&  formulation,
                    PtrParamNode infoNode ){

    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "MixedPiola"){
      usePiola_ = true;
      std::string form = SolutionTypeEnum.ToString(ACOU_PRESSURE);
      PtrParamNode potSpaceNode = infoNode->Get(form);
      crSpaces[ACOU_PRESSURE] =
              FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1);
      crSpaces[ACOU_VELOCITY] =
          FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::L2);

      crSpaces[ACOU_PRESSURE]->Init();
      crSpaces[ACOU_VELOCITY]->Init();
    }else{
      EXCEPTION("The formulation " << formulation << "of acousticMixed PDE is not known!");
    }
    return crSpaces;
  }

  void AcousticMixedPDE::DefineIntegrators(){

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

      Double density=0.0,compressibility=0.0,c0=0.0;

      materials_[actRegion]->GetScalar(density,DENSITY,Global::REAL);
      materials_[actRegion]->GetScalar(compressibility,ACOU_BULK_MODULUS,Global::REAL);
      c0 = std::sqrt(compressibility/density);

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
      if( dim_ == 2 ) {
        stiffIntPV = new ABInt<GradientOperator<FeH1,2> , IdentityOperatorPiola<FeH1,2,2> >(coeffKPV,1.0 );
      } else {
        stiffIntPV = new ABInt<GradientOperator<FeH1,3> , IdentityOperatorPiola<FeH1,3,3> >(coeffKPV,1.0 );
      }
      stiffIntPV->SetName("MixedStiffIntPV");
      //stiffIntPV->SetFeSpace( feFunctions_[ACOU_PRESSURE]->GetFeSpace(), feFunctions_[ACOU_VELOCITY]->GetFeSpace() );
      BiLinFormContext *stiffContPV = new BiLinFormContext(stiffIntPV, STIFFNESS );
      
      stiffContPV->SetEntities( actSDList, actSDList );
      stiffContPV->SetFeFunctions( feFunctions_[ACOU_PRESSURE],feFunctions_[ACOU_VELOCITY]);
      stiffContPV->SetCounterPart(true);
      assemble_->AddBiLinearForm( stiffContPV );

      
      // --------------------------------------------------------------------
      //  VERSION 2: K_VP = K_PV^T Integrator (lower off-diagonal integrator)
      // --------------------------------------------------------------------
//      shared_ptr<CoefFunction> coeffKVP
//                = CoefFunction::Generate(Global::REAL, "1.0");
//      BiLinearForm * stiffIntVP = NULL;
//      if( dim_ == 2 ) {
//        stiffIntVP = new ABInt< IdentityOperatorPiola<FeH1,2,2> , GradientOperator<FeH1,2> >(coeffKVP,1.0 );
//      } else {
//        stiffIntVP = new ABInt< IdentityOperatorPiola<FeH1,3,3> , GradientOperator<FeH1,3> >(coeffKVP,1.0 );
//      }
//      stiffIntVP->SetName("MixedStiffIntVP");
//      //stiffIntVP->SetFeSpace( feFunctions_[ACOU_PRESSURE]->GetFeSpace(), feFunctions_[ACOU_VELOCITY]->GetFeSpace() );
//      BiLinFormContext *stiffContVP = new BiLinFormContext(stiffIntVP, STIFFNESS );
//
//      stiffContVP->SetEntities( actSDList, actSDList );
//      stiffContVP->SetFeFunctions( feFunctions_[ACOU_VELOCITY],feFunctions_[ACOU_PRESSURE]);
//      stiffContVP->SetCounterPart(true);
//      assemble_->AddBiLinearForm( stiffContVP );

      // ====================================================================
      // mass integrators
      // ====================================================================
      shared_ptr<CoefFunction> coeffMPP
                = CoefFunction::Generate(Global::REAL, lexical_cast<std::string>(1.0));

      BiLinearForm *massIntPP = NULL;
      if( dim_ == 2 ) {
        massIntPP = new BBInt<IdentityOperator<FeH1,2,1> >(coeffMPP, -1.0 / (density * c0*c0) );
      } else {
        massIntPP = new BBInt<IdentityOperator<FeH1,3,1> >(coeffMPP, -1.0 / (density * c0*c0) );
      }
      massIntPP->SetName("PressureTimeInt");
      //massIntPP->SetFeSpace( feFunctions_[ACOU_PRESSURE]->GetFeSpace() );

      BiLinFormContext *massContextPP =  new BiLinFormContext(massIntPP, DAMPING );

      massContextPP->SetEntities( actSDList, actSDList );
      massContextPP->SetFeFunctions( feFunctions_[ACOU_PRESSURE],feFunctions_[ACOU_PRESSURE]);
      assemble_->AddBiLinearForm( massContextPP );

      shared_ptr<CoefFunction> coeffMVV
                = CoefFunction::Generate(Global::REAL, lexical_cast<std::string>(1.0));

      BiLinearForm *massIntVV = NULL;
      if( dim_ == 2 ) {
        massIntVV = new BBInt<IdentityOperatorPiola<FeH1,2,2> >(coeffMVV, density);
      } else {
        massIntVV = new BBInt<IdentityOperatorPiola<FeH1,3,3> >(coeffMVV, density);
      }
      massIntVV->SetName("VelocityTimeInt");
      //massIntVV->SetFeSpace( feFunctions_[ACOU_VELOCITY]->GetFeSpace() );

      BiLinFormContext *massContextVV =  new BiLinFormContext(massIntVV, DAMPING );

      massContextVV->SetEntities( actSDList, actSDList );
      massContextVV->SetFeFunctions( feFunctions_[ACOU_VELOCITY],feFunctions_[ACOU_VELOCITY]);
      assemble_->AddBiLinearForm( massContextVV );
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
          abcInt = new BBInt<IdentityOperator<FeH1,2,1> >(coeffKPV,-1.0/(density*c0) );
        } else {
          abcInt = new BBInt<IdentityOperator<FeH1,3,1> >(coeffKPV,-1.0/(density*c0) );
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
      velocity->unit = "m^2/s";

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

   }

   void AcousticMixedPDE::InitTimeStepping(){
     shared_ptr<BaseTimeScheme> mySchemeV(new TimeSchemeGLM(TimeSchemeGLM::TRAPEZOIDAL, 0) );
     shared_ptr<BaseTimeScheme> mySchemeP(new TimeSchemeGLM(TimeSchemeGLM::TRAPEZOIDAL, 0) );

     feFunctions_[ACOU_PRESSURE]->SetTimeScheme(mySchemeP);
     feFunctions_[ACOU_VELOCITY]->SetTimeScheme(mySchemeV);

   }
}
