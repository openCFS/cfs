// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <stdlib.h>
#include <complex>
#include <memory>
#include <sstream>
#include <utility>

#include "CoupledPDE/pdecoupling.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/Logging/log.hpp"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Scripting/scriptable.hh"
#include "DataInOut/WriteInfo.hh"
#include "Domain/ansatzFct.hh"
#include "Domain/domain.hh"
#include "Domain/elem.hh"
#include "Domain/entityList.hh"
#include "Domain/grid.hh"
#include "Domain/resultInfo.hh"
#include "Domain/surfElem.hh"
#include "Driver/assemble.hh"
#include "Driver/formsContext.hh"
#include "Driver/solveStepFluidMech.hh"
#include "Elements/basefe.hh"
#include "Forms/baseForm.hh"
#include "Forms/fluidMechMassInt.hh"
#include "Forms/fluidMechShearStress.hh"
#include "Forms/fluidMechStiffInt.hh"
#include "Forms/linElastInt.hh"
#include "Forms/linPressureInt.hh"
#include "Forms/linSurfForm.hh"
#include "Forms/linearForm.hh"
#include "Forms/massInt.hh"
#include "Forms/mechStressStrain.hh"
#include "General/Enum.hh"
#include "General/exception.hh"
#include "MatVec/SingleVector.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/exprt/xpr2.hh"
#include "Materials/baseMaterial.hh"
#include "OLAS/algsys/basesystem.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/StdPDE.hh"
#include "PDE/basePDE.hh"
#include "PDE/eqnMap.hh"
#include "PDE/timestepping.hh"
#include "Utils/basenodestoresol.hh"
#include "Utils/result.hh"
#include "Utils/tools.hh"
#include "bdf2.hh"
#include "fluidMechPDE.hh"

namespace CoupledField
{

DECLARE_LOG(fluidmechpde)
DEFINE_LOG(fluidmechpde, "fluidmechpde")


FluidMechPDE::FluidMechPDE(Grid * aptgrid, PtrParamNode paramNode )
:SinglePDE( aptgrid, paramNode ) {

  pdename_          = "fluidMech";
  pdematerialclass_ = FLOW;
  maxTimeDerivOrder_ = 2;

  nonLin_      = false;

  needSolPrev_ = true;

  //set coordinate Update true to obtain the updated Lagrange formulation
  coordUp_ = true;

  gridSol_ = new NodeStoreSol<Double>;
  // ****************************
  // DETERMINE GEOMETRY
  // ****************************

  // Get problem geometry and PDE subtype
  myParam_->GetValue("subType", subType_ );

  std::string probGeo;
  param->Get("domain")->GetValue("geometryType", probGeo );

  // Set number of degrees of freedom and
  // ensure that subtype fits to problem geometry
  if ( subType_ == "3d" && probGeo == "3d" ) {
    EXCEPTION("3d not implement yet.");
    stressDim_ = 6;
    Info->PrintF("", "=== 3D PROBLEM\n");
  } else if ( subType_ == "axi" && probGeo == "axi" ) {
    isaxi_ = true;
    stressDim_ = 4;
    Info->PrintF("", "=== AXISYSMMETRIC PROBLEM\n");
  } else if ( subType_ == "plane" && probGeo == "plane" ) {
    stressDim_ = 3;
    Info->PrintF("", "=== PLANE PROBLEM\n");
  } else {
    EXCEPTION( "Subtype '" <<  subType_ << "' of PDE '"
        <<  pdename_ <<  "' does not fit to problem  geometry '"
        << probGeo << "'"; );
  }

  approxType_ = myParam_->Get("type")->As<std::string>();
  stabilizationType_ = myParam_->Get("stabilizationType")->As<std::string>();
  movingMesh_ = myParam_->Get("movingMesh")->As<bool>();
  saveAcouSrc_= myParam_->Get("aeroAcoustics")->As<bool>();

  PtrParamNode nonLinNode;
  nonLinNode = myParam_->Get("nonLinear", ParamNode::PASS );
  if ( !nonLinNode ) {
    EXCEPTION("Navier-Stokes in non linear and needs a non linear solver. Please provide one.");
  }
  // how to solve non linearity
  nonLinMethod_ = nonLinNode->Get("method")->As<std::string>();
  if (nonLinMethod_ != "newton")
  {
    EXCEPTION("non linear solver \"" << nonLinMethod_ << "\" not implemented");
  }
  if (approxType_ != "lagrange")
  {
    EXCEPTION("approximation type \"" << approxType_ << "\" not implemented");
  }

  forceFac_ = 1.0;

  // Register functions for scripting
  //RegisterFunctions();
}

FluidMechPDE::~FluidMechPDE()
{

}

void FluidMechPDE::ReadSpecialBCs() {

  // read pressure information
  ReadPressureLoads();

  // ***************************************************************
  //   ABC for flow
  // ***************************************************************
  PtrParamNode bcNode = myParam_->Get( "bcsAndLoads", ParamNode::PASS );
  
  if( bcNode ) {
    ParamNodeList abcNodes = bcNode->GetList( "absorbingBCs" );
    
    for ( UInt i=0, n=abcNodes.GetSize(); i<n; ++i )
      {
        std::string regionName = abcNodes[i]->Get("name")->As<std::string>();
        absBCs_.Push_back( ptgrid_
                           ->GetEntityList( EntityList::SURF_ELEM_LIST,
                                            regionName, EntityList::REGION ) );
        Info->PrintF( pdename_,
                      "Apply Absorbing Boundary Conditions on surfaceRegion '%s'\n",
                      regionName.c_str() );
      }
  }

}

void FluidMechPDE::InitNonLin()
{
  // ==============================================================
  // NOTE: We assume that for a fluid mechanical PDE all regions
  //       either are non-linear due to the convective term!
  // ==============================================================

  nonLin_ = true;

}


void FluidMechPDE::DefineIntegrators()
{
  Double density, dynamicViscosity, kinematicViscosity;

  //for viscous damping
  Double startPos, endPos, dampFactor;   //for viscous damping
  std::string dampingType;

  // help variables for parameter checking
  RegionIdType actRegion;
  BaseMaterial * actSDMat = NULL;

  //volume integrators

  std::map<RegionIdType, BaseMaterial*>::iterator it;
  for ( it = materials_.begin(); it != materials_.end(); it++ )
  {
    // Set current region and material
    actRegion = it->first;
    actSDMat = it->second;

    actSDMat->GetScalar(density,DENSITY,Global::REAL);
    actSDMat->GetScalar(dynamicViscosity,DYNAMIC_VISCOSITY,Global::REAL);
    actSDMat->GetScalar(kinematicViscosity,KINEMATIC_VISCOSITY,Global::REAL);

    // create new entity list
    shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
    actSDList->SetRegion( actRegion );

    // get current region name and get grip on paramNode
    std::string actRegionName;
    actRegionName = ptgrid_->GetRegion().ToString(actRegion);
    const UInt spaceDim = ptgrid_->GetDim();
    if (spaceDim != 2)
    {
      EXCEPTION("Wrong space dim!");
    }

    // Check for viscous damping; currently defined by PML structure
    bool isViscousDamp = false;
    startPos = 0.0;
    endPos   = 0.0;
    if ( dampingList_[actRegion] == PML ) {
      isViscousDamp = true;
      PtrParamNode actRegionNode =
        myParam_->Get("regionList")->GetByVal( "region", "name", actRegionName );

      //read data for PML layer
      Matrix<Double> inner;
      Matrix<Double> outer;
      std::string coordSysId;
      
      std::string id = actRegionNode->Get("dampingId")->As<std::string>();
      PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal("pml", "id", id);
      ReadDataPML(dampingType, inner, dampFactor, coordSysId, pmlNode );

      //extract data, we need for viscous damping
      startPos = inner[0][0];
      endPos   = inner[1][0];
      //std::cout << "start: " << startPos << "  end: " << endPos << std::endl;
    }


    // ==============  add stiffness matrix ===============

    //****************************************************************
    //******** stabilized ********************************************
    //****************************************************************
    BaseForm * bilinearStiff_UV = GetStiffIntegrator_UV(density,kinematicViscosity);
    if ( isViscousDamp )
      bilinearStiff_UV->InitSpongeLayer(startPos, endPos, dampFactor, dampingType);

    BiLinFormContext * actStiffContext_UV = new BiLinFormContext(bilinearStiff_UV, STIFFNESS );
    actStiffContext_UV->SetCounterPart( false );
    actStiffContext_UV->SetPtPdes(this, this);
    actStiffContext_UV->SetResults( results_[0],results_[0],actSDList, actSDList );
    assemble_->AddBiLinearForm( actStiffContext_UV );
    // use SUPG ********************************************************************************************************
    //******************************************************************************************************************
    if( stabilizationType_ != "none"){	
      BaseForm * bilinearStiff_PQ = GetStiffIntegrator_PQ(density,kinematicViscosity);
      if ( isViscousDamp )
        bilinearStiff_PQ->InitSpongeLayer(startPos, endPos, dampFactor, dampingType);

      BiLinFormContext * actStiffContext_PQ = new BiLinFormContext(bilinearStiff_PQ, STIFFNESS );
      actStiffContext_PQ->SetCounterPart( false );
      actStiffContext_PQ->SetPtPdes(this, this);
      actStiffContext_PQ->SetResults( results_[1],results_[1],actSDList, actSDList );
      assemble_->AddBiLinearForm( actStiffContext_PQ );
    }

    //******************************************************************************************************************
    if (approxType_!= "taylorHood"){
      BaseForm * bilinearStiff_UQ = GetStiffIntegrator_UQ(density,kinematicViscosity);
      if ( isViscousDamp )
        bilinearStiff_UQ->InitSpongeLayer(startPos, endPos, dampFactor, dampingType);

      BiLinFormContext * actStiffContext_UQ = new BiLinFormContext(bilinearStiff_UQ, STIFFNESS );
      actStiffContext_UQ->SetCounterPart( false );
      actStiffContext_UQ->SetPtPdes(this, this);
      actStiffContext_UQ->SetResults( results_[1],results_[0],actSDList, actSDList );
      assemble_->AddBiLinearForm( actStiffContext_UQ );
    }
    //******************************************************************************************************************
    BaseForm * bilinearStiff_PV = GetStiffIntegrator_PV(density,kinematicViscosity);
    if ( isViscousDamp )
      bilinearStiff_PV->InitSpongeLayer(startPos, endPos, dampFactor, dampingType);

    BiLinFormContext * actStiffContext_PV = new BiLinFormContext(bilinearStiff_PV, STIFFNESS );
    // We also need to set the transposed of the coupling
    // matrix to the lower diagonal side
    actStiffContext_PV->SetCounterPart( false );
    if (approxType_== "taylorHood")
      actStiffContext_PV->SetCounterPart( true );
    actStiffContext_PV->SetPtPdes(this, this);
    actStiffContext_PV->SetResults( results_[0],results_[1],actSDList, actSDList );
    assemble_->AddBiLinearForm( actStiffContext_PV );

    // ==============  add mass matrix ====================
    if (analysistype_ == TRANSIENT)
    {
      isInstationary_=true;
      //****************************************************************
      //******** stabilized ********************************************
      //****************************************************************
      //****************************************************************
      BaseForm * bilinearMass_UV = GetMassIntegrator_UV(density,kinematicViscosity);
      if ( isViscousDamp )
        bilinearMass_UV->InitSpongeLayer(startPos, endPos, dampFactor, dampingType);

      BiLinFormContext * actMassContext_UV = new BiLinFormContext(bilinearMass_UV, MASS );
      actMassContext_UV->SetCounterPart( false );
      actMassContext_UV->SetPtPdes(this, this);
      actMassContext_UV->SetResults(results_[0],results_[0],actSDList, actSDList );
      assemble_->AddBiLinearForm( actMassContext_UV );

      if(stabilizationType_ != "none")
      {	
        //****************************************************************
        BaseForm * bilinearMass_UQ = GetMassIntegrator_UQ(density,kinematicViscosity);
        if ( isViscousDamp )
          bilinearMass_UQ->InitSpongeLayer(startPos, endPos, dampFactor, dampingType);

        BiLinFormContext * actMassContext_UQ = new BiLinFormContext(bilinearMass_UQ, MASS );
        actMassContext_UQ->SetCounterPart( false );
        actMassContext_UQ->SetPtPdes(this, this);
        actMassContext_UQ->SetResults(results_[1],results_[0],actSDList,actSDList);
        assemble_->AddBiLinearForm( actMassContext_UQ );
        //****************************************************************
      }
      //****************************************************************
      // This is for updating the Newton RHS
      BaseForm * bilinearAux_UV = GetAuxIntegrator_UV(density, kinematicViscosity);
      if ( isViscousDamp )
        bilinearAux_UV->InitSpongeLayer(startPos, endPos, dampFactor, dampingType);

      BiLinFormContext * actAuxContext_UV = new BiLinFormContext(bilinearAux_UV, AUXILIARY );
      actAuxContext_UV->SetCounterPart( false );
      actAuxContext_UV->SetPtPdes(this, this);
      actAuxContext_UV->SetResults(results_[0], results_[0], actSDList, actSDList);
      assemble_->AddBiLinearForm( actAuxContext_UV );
      //***********V*****************************************************
      //****************************************************************
      // This is for updating the Newton RHS
      BaseForm * bilinearAux_UQ = GetAuxIntegrator_UQ(density, kinematicViscosity);
      if ( isViscousDamp )
        bilinearAux_UQ->InitSpongeLayer(startPos, endPos, dampFactor, dampingType);

      BiLinFormContext * actAuxContext_UQ = new BiLinFormContext(bilinearAux_UQ, AUXILIARY );
      actAuxContext_UQ->SetCounterPart( false );
      actAuxContext_UQ->SetPtPdes(this, this);
      actAuxContext_UQ->SetResults(results_[1], results_[0], actSDList, actSDList);
      assemble_->AddBiLinearForm( actAuxContext_UQ );
      //***********V*****************************************************
    }
    eqnMap_->AddResult( *results_[0], actSDList );
    eqnMap_->AddResult( *results_[1], actSDList );
  } // loop over materials

  //surface integrators
  //RHS-part
  for (UInt actSF = 0; actSF < pressSurf_.GetSize(); actSF++)
  {
    LinearSurfForm * rhsSrcSurf =
      new PressureLinForm(pressVals_[actSF], pressPhase_[actSF],
          subType_, isaxi_ );
    rhsSrcSurf->SetVoluInfo( materials_ );

    LinearFormContext * pressRhs =
      new LinearFormContext( rhsSrcSurf );
    pressRhs->SetPtPde( this );
    pressRhs->SetResult( results_[0], pressSurf_[actSF] );
    assemble_->AddLinearForm( pressRhs );

    // Give entities and result to equation numbering class
    // and solution class
    eqnMap_->AddResult( *results_[0], pressSurf_[actSF] );
  }

  // **********************************************************************
  //   surface-integration: Absorbing boundaries
  // **********************************************************************
  if ( absBCs_.GetSize() > 0 ) {
    for (UInt actSD=0, numSD=absBCs_.GetSize(); actSD<numSD; ++actSD) {
      std::cout << "Do ANC" << std::endl;
       FluidMechAbsorbingFlow *linearAbsorb =
        new FluidMechAbsorbingFlow(kinematicViscosity);
       linearAbsorb->SetVoluInfo( materials_ );
       linearAbsorb->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
      
      LinearFormContext * abcContext =
        new LinearFormContext( linearAbsorb );
      abcContext->SetPtPde(this);
      abcContext->SetResult( results_[0], absBCs_[actSD] );
      assemble_->AddLinearForm( abcContext );
      
      // Give result to equation numbering class
      eqnMap_->AddResult( *results_[0], absBCs_[actSD] );
    }
  }

}


BaseForm * FluidMechPDE::GetStiffIntegrator(BaseMaterial* actSDMat,
    RegionIdType regionId,
    bool reducedInt) {
  // Get region name
  std::string regionName = ptgrid_->GetRegion().ToString(regionId);

  BaseForm * bilinearStiff = NULL;

  //transform the type
  SubTensorType tensorType;
  String2Enum(subType_,tensorType);

  bilinearStiff = new linElastInt(actSDMat, tensorType);

  return bilinearStiff;
}


BaseForm * FluidMechPDE::GetStiffIntegrator_UV(Double density,
    Double kinematicViscosity)
{
  BaseForm * bilinearStiff = NULL;
  if (subType_ == "plane")
  {
    if( stabilizationType_ != "none" && approxType_ == "lagrange" )
    {
      if (nonLinMethod_ == "newton")
      {
        bilinearStiff = new FluidMechPlaneStiffNewtonInt_UV(density, kinematicViscosity,
            stabilParams_, movingMesh_, stabilizationType_);
        bilinearStiff->setTimeStepping(TS_alg_);
      } else
        EXCEPTION("non linear solver \"" << nonLinMethod_ << "\" not implemented");
    } else
      EXCEPTION("Unknown discretizationType in fluidMech PDE! ");
  } else
    EXCEPTION("Unknown subtype in fluidMech PDE! ");

  bilinearStiff->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
  bilinearStiff->SetGridSolution( dynamic_cast<NodeStoreSol<Double>&>(*gridSol_ ));
  return bilinearStiff;
}

BaseForm * FluidMechPDE::GetStiffIntegrator_PQ(Double density,Double kinematicViscosity) {
  BaseForm * bilinearStiff = NULL;
  if (subType_ == "plane")
  {
    if (nonLinMethod_ == "newton")
    {
      bilinearStiff = new FluidMechPlaneStiffNewtonInt_PQ(density, kinematicViscosity, \
          stabilParams_, movingMesh_, stabilizationType_);
      bilinearStiff->setTimeStepping(TS_alg_);
    } else
      EXCEPTION("non linear solver \"" << nonLinMethod_ << "\" not implemented");
  } else
    EXCEPTION("Unknown subtype in fluidMech PDE! ");
  bilinearStiff->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
  bilinearStiff->SetGridSolution( dynamic_cast<NodeStoreSol<Double>&>(*gridSol_ ));
  return bilinearStiff;
}

BaseForm * FluidMechPDE::GetStiffIntegrator_UQ(Double density,Double kinematicViscosity) {
  BaseForm * bilinearStiff = NULL;
  if (subType_ == "plane")
  {
    if( stabilizationType_ != "none" && approxType_ == "lagrange" )
    {
      if (nonLinMethod_ == "newton")
      {
        bilinearStiff = new FluidMechPlaneStiffNewtonInt_UQ(density, kinematicViscosity, \
            stabilParams_, movingMesh_, stabilizationType_);
        bilinearStiff->setTimeStepping(TS_alg_);
      } else
        EXCEPTION("non linearity not supported: " << nonLinMethod_);
    } else {
      EXCEPTION("Unknown discretizationType in fluidMech PDE! ");
    }
  } else {
    EXCEPTION("Unknown subtype in fluidMech PDE! ");
  }
  bilinearStiff->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
  bilinearStiff->SetGridSolution( dynamic_cast<NodeStoreSol<Double>&>(*gridSol_ ));
  return bilinearStiff;
}

BaseForm * FluidMechPDE::GetStiffIntegrator_PV(Double density,Double kinematicViscosity)
{
  BaseForm * bilinearStiff = NULL;
  if (subType_ == "plane")
  {
    if( stabilizationType_ != "none" && approxType_ == "lagrange" )
    {
      if (nonLinMethod_ == "newton")
      {
        bilinearStiff = new FluidMechPlaneStiffNewtonInt_PV(density, kinematicViscosity, \
            stabilParams_, movingMesh_, stabilizationType_);
        bilinearStiff->setTimeStepping(TS_alg_);
      } else
        EXCEPTION("non linear solver \"" << nonLinMethod_ << "\" not implemented");
    } else {
      EXCEPTION("Unknown discretizationType in fluidMech PDE! ");
    }
  } else {
    EXCEPTION("Unknown subtype in fluidMech PDE! ");
  }

  bilinearStiff->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
  bilinearStiff->SetGridSolution( dynamic_cast<NodeStoreSol<Double>&>(*gridSol_ ));
  return bilinearStiff;
}

BaseForm * FluidMechPDE::GetMassIntegrator_UV(Double density,Double kinematicViscosity)
{
  BaseForm * bilinearMass = NULL;

  if (subType_ == "plane")
  {
    if( stabilizationType_ != "none" && approxType_ == "lagrange" )
    {
      if (nonLinMethod_ == "newton")
      {
        bilinearMass = new FluidMechPlaneMassNewtonInt_UV(density, kinematicViscosity, \
            stabilParams_, movingMesh_, stabilizationType_);
        bilinearMass->setTimeStepping(TS_alg_);
      } else
        EXCEPTION("non linear solver \"" << nonLinMethod_ << "\" not implemented");
    } else
      EXCEPTION("Unknown discretizationType in fluidMech PDE! ");
  } else
    EXCEPTION("Unknown subtype in fluidMech PDE! ");

  bilinearMass->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
  bilinearMass->SetGridSolution( dynamic_cast<NodeStoreSol<Double>&>(*gridSol_ ));
  return bilinearMass;
}

BaseForm * FluidMechPDE::GetMassIntegrator_UQ(Double density,Double kinematicViscosity)
{
  BaseForm * bilinearMass = NULL;

  if (subType_ == "plane")
  {
    if (nonLinMethod_ == "newton")
    {
      bilinearMass = new FluidMechPlaneMassNewtonInt_UQ(density, kinematicViscosity, \
          stabilParams_, movingMesh_, stabilizationType_);
      bilinearMass->setTimeStepping(TS_alg_);
    } else
      EXCEPTION("non linear solver \"" << nonLinMethod_ << "\" not implemented");
  } else
    EXCEPTION("Unknown subtype in fluidMech PDE! ");

  bilinearMass->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_));
  bilinearMass->SetGridSolution( dynamic_cast<NodeStoreSol<Double>&>(*gridSol_));
  return bilinearMass;
}

BaseForm * FluidMechPDE::GetAuxIntegrator_UV(Double density,Double kinematicViscosity)
{
  BaseForm * bilinearAux = NULL;

  if (subType_ == "plane")
  {
    if (nonLinMethod_ == "newton")
    {
      bilinearAux = new FluidMechPlaneIntNewton_RhsUV(density, kinematicViscosity, \
          stabilParams_, movingMesh_, stabilizationType_);
      bilinearAux->setTimeStepping(TS_alg_);
    } else
      EXCEPTION("Non linear method: \"" << nonLinMethod_<< "\" not implemented.");
  } else {
      EXCEPTION("Not implemented for subtype: \"" << subType_);
  }

  bilinearAux->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_));
  bilinearAux->SetGridSolution( dynamic_cast<NodeStoreSol<Double>&>(*gridSol_));
  return bilinearAux;
}

BaseForm * FluidMechPDE::GetAuxIntegrator_UQ(Double density,Double kinematicViscosity)
{
  BaseForm * bilinearAux = NULL;

  if (subType_ == "plane")
  {
    if (nonLinMethod_ == "newton")
    {
      bilinearAux = new FluidMechPlaneIntNewton_RhsUQ(density, kinematicViscosity, \
          stabilParams_, movingMesh_, stabilizationType_);
      bilinearAux->setTimeStepping(TS_alg_);
    } else
      EXCEPTION("Non linear method: \"" << nonLinMethod_<< "\" not implemented.");
  } else {
      EXCEPTION("Not implemented for subtype: \"" << subType_);
  }

  bilinearAux->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_));
  bilinearAux->SetGridSolution( dynamic_cast<NodeStoreSol<Double>&>(*gridSol_));
  return bilinearAux;
}

void FluidMechPDE::DefineSolveStep()
{
  solveStep_ = new SolveStepFluidMech(*this);
}




// ======================================================
// COUPLING SECTION
// ======================================================


void FluidMechPDE::InitCoupling(PDECoupling * Coupling)
{
  isIterCoupled_ = true;
  ptCoupling_   = Coupling;

  const UInt numOutputCouplings = ptCoupling_->GetNumOutputCouplings();

  StdVector<StdVector<UInt> > elemNodeToCouplingNode_tmp;
  elemNodeToCouplingNode_.Resize(numOutputCouplings);
  elemNodeToCouplingNode_.Init();

  // input couplings
  for (UInt i=0; i<ptCoupling_->GetNumInputCouplings(); i++) {
    //          // check for input mechanic velocity
    //          if (ptCoupling_->GetInputQuantity(i) == MECH_VELOCITY) {
    //            numCouplingBcs_ += (dim_ * ptCoupling_->GetInputNumNodes(i));
    //          }
    // check for input smooth velocity
    if (ptCoupling_->GetInputQuantity(i) == SMOOTH_VELOCITY) {
      numCouplingBcs_ += (dim_ * ptCoupling_->GetInputNumNodes(i));
    }
  }

  for (UInt i = 0; i < numOutputCouplings; i++) {
    if (ptCoupling_->GetOutputQuantity(i) == FLUIDMECH_FORCE) {
      // Intialize the memory of the coupling values
      ptCoupling_->CreateCouplingVector(i,isComplex_);

      // now since we need a incremental formulation,
      //  initialize some necessary vectors
      isIncrFormulation_ = true;
    } else if (ptCoupling_->GetOutputQuantity(i) == ACOU_RHSVAL) {
      // Intialize the memory of the coupling values
      ptCoupling_->CreateCouplingVector(i,isComplex_);

      //now since we need a incremental formulation, initialize some necessary vectors
      isIncrFormulation_ = true;
    } else if (ptCoupling_->GetOutputQuantity(i) == ACOUSURF_RHSVAL) {
      // Intialize the memory of the coupling values
      ptCoupling_->CreateCouplingVector(i,isComplex_);

      //now since we need a incremental formulation, initialize some necessary vectors
      isIncrFormulation_ = true;
    } else {
      EXCEPTION("unknown coupling quantitiy");
    }

  }
}


void FluidMechPDE::CalcInputCoupling() {

  std::string errMsg;
  StdVector<UInt> * nodes;
  SingleVector * val;
  Integer eqnNr;
  UInt couplingDof;
  // Determine maximal allowed equation number for algebraic system
  Integer maxAllowedEqn = (Integer)algsys_->GetDimension();

  // at first, check if this PDE is iterative coupled
  if (isIterCoupled_ == false)
    return;

  // Outer loop over all INPUT coupling terms
  for (UInt i=0; i<ptCoupling_->GetNumInputCouplings(); i++) {

    //    ptCoupling_ = &ptCoupling_[i];
    ptCoupling_->GetInputValues(i, val);
    couplingDof = ptCoupling_->GetInputDof(i);

    // Up to now, Coupling is only possible with
    // Real valued solutions
    Vector<Double> const & help = dynamic_cast<Vector<Double>&>(*val);

    switch(ptCoupling_->GetInputType(i)) {

      // -------------------
      // COORDINATE COUPLING
      // -------------------
      case COORD:

        ptCoupling_->GetInputNodes(i, nodes);
        ptgrid_->SetNodeOffset( *nodes , help );

        break;

        // -------------------
        // RHS COUPLING
        // -------------------
      case RHS:
        ptCoupling_->GetInputNodes(i, nodes);

        try {
          for ( UInt dof = 0; dof < couplingDof; dof++ ) {
            for ( UInt j = 0; j < nodes->GetSize(); j++ ) {
              eqnNr = eqnMap_->GetNodeEqn((*nodes)[j],dof+1);
              if ( eqnNr != 0 && eqnNr <= maxAllowedEqn ) {
                algsys_->SetNodeRHS( help[ dof + couplingDof * j ], pdeId_,
                    eqnNr);
              }

              LOG_DBG(fluidmechpde) << "CalcInputCoupling: " << "(" << pdename_ << ") "
                << "Refused to pass " << "eqnNr=" << eqnNr << " to SetNodeRHS()."
                << " numLastFreeDof=" << maxAllowedEqn;
            }
          }
        } catch (Exception& ex ) {
          RETHROW_EXCEPTION(ex, pdename_ << ": " << std:: endl 
              << "Could not get equation number.");
        }
        break;

        // -----------------------
        // InhomDirichlet COUPLING
        // -----------------------
      case ID_BC:
        // How do we want to treat inhomogeneous Dirichlet boundary
        // conditions?
        {
          //check whether penalty formulation is used
          // with elimination this does not work currentlly
          bool usePenalty = true;
          std::string aux;
          StdVector<std::string> keyVec;
          StdVector<std::string> attrVec;
          StdVector<std::string> valVec;
          //           keyVec  = "linearSystems", "system", "setup", "idbcHandling";
          //           attrVec = "", "name", "";
          //           valVec  = "", pdename_, "";
          //           params->Get( keyVec, attrVec, valVec, aux );
          //           usePenalty = aux == "penalty" ? true : false;
          usePenalty = true;
          Info->PrintF( pdename_, "Treating IDBCs using '%s' approach\n",
              aux.c_str() );
          if ( usePenalty == false ) {
            EXCEPTION("Cannot use inhom. Dirichlet coupling together with "
                << "IDBC elimination, since the equation numbering "
                << "object does currently not have the information "
                << "required to put those values at the end of the "
                << "equation number interval! Please set idbcHandling="
                << '"' << "penalty" << '"' << " in your xml-file");
          }
        }

        // Set flag that the boundary conditions have to be incorporated
        updateCouplingBCs_ = true;

        ptCoupling_->GetInputNodes(i, nodes);

        //         LOG_DBG(fluidMechPDE) << "Input FLUIDMECH_VELOCITY\n"
        //                               << "couplingnodes:\n" << (*nodes)
        //                               << "values:\n" << help << std::endl;

        try {
          for ( UInt dof = 0; dof < dim_; dof++ )
          {
            for ( UInt j = 0; j < nodes->GetSize(); j++)
            {
              //eqnNr = eqnMap_->GetNodeEqn((*nodes)[j],dof+1 );
              eqnNr = eqnMap_->GetNodeEqn(*results_[0],(*nodes)[j],dof+1);

              if (eqnNr==0)
              {
                EXCEPTION( "The specified coupling node has no equation number");
              }

              algsys_->SetDirichlet( pdeId_, eqnNr, \
                                     help[dof+j*dim_] );
            }
          }
        } catch (Exception& ex ) {
          RETHROW_EXCEPTION(ex, pdename_ << ": " << std:: endl 
              << "Could not get equation number.");
        }
        break;

      case MAT:
        EXCEPTION( "Not implemented yet" );
        break;

      case GRID_VEL:
        ptCoupling_->GetInputNodes(i, nodes);

        //         LOG_DBG(fluidMechPDE) << "Input GRID_VELOCITY\n"
        //                               << "couplingnodes:\n" << (*nodes)
        //                               << "values:\n" << help << std::endl;

        //gridSol_->CouplingToNodeSolution(help,(*nodes));
        try {
          for ( UInt dof = 0; dof < dim_; dof++ ) {
            for ( UInt j = 0; j < nodes->GetSize();j++) {
              eqnNr = eqnMap_->GetNodeEqn(*results_[0],(*nodes)[j],dof+1);
              //gridSol_->data_[abs(eqnNr-1)]=help[dof+j*dim_];
              gridSol_->Set(eqnNr,help[dof+j*dim_]);
            }
          }
        } catch (Exception& ex ) {
          RETHROW_EXCEPTION(ex, pdename_ << ": " << std:: endl 
              << "Could not get equation number.");
        }
        break;
    }
    //       for (UInt ii=0; ii<help.GetSize(); ii++)
    //         Info->PrintF( pdename_, "help[%d]=%e\n", ii, help[ii]);

  }
}

void FluidMechPDE::CalcOutputCoupling()
{
  UInt dof;
  SolutionType quantity;
  StdVector<Elem*> * couplingElems = NULL;
  StdVector<UInt> * couplingNodes = NULL;
  SingleVector * temp_values = NULL;
  // TODO: Check if this is still needed
  //SingleVector * temp_oldValues = NULL;
  //UInt regionCount = 0;

  //        Info->PrintF( "FACTORIZATION", "The fluid mechanical forces will be factorized by forceFactor = %e\n",forceFac_);

  // at first, check if this PDE is iterative coupled
  if (isIterCoupled_ == false)
    return;

  // loop over all output coupling quantities
  for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++) {

    //Welche Größe ist die Koppelgröße??
    //gibt einen Solution Type zurück
    quantity = ptCoupling_->GetOutputQuantity(i);
    ptCoupling_->GetOutputValues(i, temp_values);
    // TODO: Check if this is still needed
    // ptCoupling_->GetOutputOldValues(i, temp_oldValues);

    // hard coded cast, since coupling is only possible with
    // real valued entries
    Vector<Double> * values = dynamic_cast<Vector<Double>*>(temp_values);
    // TODO: Check if this is still needed
    // Vector<Double> * oldValues = dynamic_cast<Vector<Double>*>(temp_oldValues);
    Vector<Double> auxValues;
    auxValues.Resize(values->GetSize());
    //Info->PrintF( pdename_, "Values_Size=%d\n", values->GetSize());
    //for(UInt i=0; i< values->GetSize(); i++)
    //Info->PrintF( pdename_, "Values[%d]=%d\n", i, (*values)[i]);


    switch(ptCoupling_->GetOutputType(i))
    {
      case NODE:
        ptCoupling_->GetOutputNodes(i, couplingNodes);

        if (quantity == FLUIDMECH_FORCE)
        {
          ptCoupling_->GetOutputElements(i, couplingElems);
          dof = ptCoupling_->GetOutputDof(i);
          CalcMechCouplingRHS(couplingElems, *couplingNodes, auxValues, dof);
          //(*values)=(auxValues*forceFac_) + ((*oldValues)*(1.0-forceFac_));
          (*values)=(auxValues*forceFac_);
          Info->PrintF( pdename_, "Coupling Output Type=Force\n");
        }
        if (quantity == ACOU_RHSVAL) {
          if(saveAcouSrc_ == true)
          {
            acou_src_.NodeSolutionToCoupling(*values, *couplingNodes, 1);
            Info->PrintF( pdename_, "Coupling Output Type=acouRHSval\n");
          } else {
              EXCEPTION("Sorry, but acouRHSval is not defined as a "
                  << "NodeResult in your XML, please do so");
          }

        }
        if (quantity == ACOUSURF_RHSVAL)
        {
          ptCoupling_->GetOutputElements(i, couplingElems);
          dof = ptCoupling_->GetOutputDof(i);
          CalcAcouSurfSourceCouplingRHS(couplingElems, \
                        *couplingNodes, (*values), dof);
          Info->PrintF( pdename_, "Coupling Output Type=acouSurfRHSVal\n");
        }
        break;

      case ELEM:
        EXCEPTION("No Element coupling output");
    }
  }
}

void FluidMechPDE::GetPresSolVecOfElement( Vector<Double>& elemSol,
    StdVector<UInt>& connecth ) {

  // fluidMech pressure of element nodes
  elemSol.Resize(1 * connecth.GetSize());
  elemSol.Init();
  Integer eqnNr;
  //UInt presDof;

  NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
  Vector<Double> & sol = solhelp->GetAlgSysVector();

  try {
    for(UInt actNode=0; actNode<connecth.GetSize(); actNode++) {
      //eqnNr = eqnMap_->GetNodeEqn(connecth[actNode],presDof);
      eqnNr = eqnMap_->GetNodeEqn(*results_[1],connecth[actNode],1);
      if (eqnNr!= 0) {
        elemSol[actNode] = sol[(abs(eqnNr-1))];
        //elemSol[actNode] = presSol[(abs(eqnNr-1))];
      } else {
        elemSol[actNode] = 0.0;
      }
    }
  } catch (Exception& ex ) {
    RETHROW_EXCEPTION(ex, pdename_ << ": " << std:: endl 
        << "Could not get equation number.");
  }
}

template <class TYPE>
void FluidMechPDE::GetVeloSolMatOfElement( Matrix<TYPE>& elemSol,
    StdVector<UInt>& connecth ) {

  // fluidMech velocities of element nodes
  elemSol.Resize(dim_,connecth.GetSize());
  elemSol.Init();
  Integer eqnNr;
  //UInt veloEndDof;

  // Determine maximal allowed equation number for algebraic system
  Integer maxAllowedEqn = (Integer)algsys_->GetDimension();

  NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
  Vector<Double> & sol = solhelp->GetAlgSysVector();

  try {
    for(UInt actNode=0; actNode<connecth.GetSize(); actNode++)
    {
      for(UInt actDof=0; actDof<dim_; actDof++)
      {
        eqnNr = eqnMap_->GetNodeEqn(*results_[0],connecth[actNode],actDof+1);
        if ( eqnNr != 0 )
        {
          elemSol[actDof][actNode] = sol[eqnNr-1];
        } else if ( eqnNr <= maxAllowedEqn ) {
          elemSol[actDof][actNode] = 0.0;
        } else {
          EXCEPTION( "wrong equation Number");
        }
      }
    }
  } catch (Exception& ex ) {
    RETHROW_EXCEPTION(ex, pdename_ << ": " << std:: endl 
        << "Could not get equation number.");
  }
}

void FluidMechPDE::GetVeloDeriv1SolVecOfElement( Vector<Double>& sol,
    const EntityIterator& it)
{
  if (analysistype_ == TRANSIENT) {
    StdVector<Integer> eqns;
    eqnMap_->GetEqns( eqns, *results_[0], it );
    sol.Resize( eqns.GetSize() );
    sol.Init();

    const Vector<Double> & sol_der1 = getDeriv(FIRST_DERIV);

    for ( UInt i = 0; i < eqns.GetSize(); i++ ) {
      if ( eqns[i] != 0 ) {
        sol[i] = sol_der1[abs(eqns[i])-1];
      } else {
        sol[i] = 0.0;
      }
    }
  }
}

void FluidMechPDE::AcouSourceCalc()
{
  /******************************************************************************************
   * Hier die Akustik-Flow Kopplung
   ******************************************************************************************/
  //Hier die Funktion CalcElemVector4Quad aufrufen

  if(saveAcouSrc_==true){

    UInt ii, elsize = 0;
    Integer eqnNr;

    StdVector<Elem*> elemssd;
    StdVector<UInt> connecth;
    StdVector<Integer> connect_PDE;
    Vector<Double>  Result;

    Matrix<Double> ptCoordNodes;
    Matrix<Double> flowdata;

    BaseFE * ptEl;

    //Reset the values in nodestresol, due to the usage of the function "Add"
    acou_src_.Init();

    // loop over all subdomains
    for (UInt isd=0; isd<subdoms_.GetSize(); isd++)
    {
      //Get the density of the current subdomain
      Double density;
      BaseMaterial * actMat = materials_[subdoms_[isd]];
      actMat->GetScalar(density,DENSITY,Global::REAL);

      ElemList actSDList(ptgrid_ );
      actSDList.SetRegion( subdoms_[isd] );
      EntityIterator it = actSDList.GetIterator();
      //           if(ptgrid_->RegionIdToName(subdoms_[isd]) == "fluid"){
      // loop over all elements of subdomain
      for ( it.Begin(); !it.IsEnd(); it++ ) {
        //Pointer to the reference element
        const Elem * actCoupleElem = it.GetElem();
        ptEl = actCoupleElem->ptElem;

        //generate the object LinearFlowNoiseInt which
        //is defined in linearForm.cc
        LinearFlowNoiseInt * linear_load = new LinearFlowNoiseInt(ptEl);

        //Number of Nodes of the Element
        elsize=(actCoupleElem->connect).GetSize();
        //the length of the connecth is set to the number of nodes
        //because it stores the nodesnumbers of the element
        connecth.Resize(elsize);
        for (ii=0; ii<elsize; ii++)
          connecth[ii]=(actCoupleElem->connect)[ii];


        //Schreibt die Koordinaten der Knoten des Elements in ptCoordNodes
        Matrix<Double> ptCoordNodes;
        ptgrid_->GetElemNodesCoord(ptCoordNodes, connecth);

        //as sad as it is, but we have to bring the fluid calculation results
        //in Matrix notation here :-(
        GetVeloSolMatOfElement(flowdata, connecth);
        //in case, that we really use the sorting funtionality, for the results, in CalcElemVector4Quad
        //we have to add one row for the pressure!
        //flowdata = elem_data;

        //Now we have all the stuff needed to calculate the quadrupole source term,
        //CalcElemVector4Quad, gets the coordinates of the Nodes, the Node Numbers, the flowdata
        //and returns Result
        //              LOG_DBG(fluidMechPDE) << "The current Element is:\n" << actCoupleElem->elemNum << std::endl;
        linear_load->CalcElemVec_withCFSVel(ptCoordNodes, flowdata, Result, density);
        // Now the results are written out as a gid result
        try {
          for(ii = 0; ii< ptEl->GetNumNodes(); ii++)
          {
            eqnNr = eqnMap_->GetNodeEqn(*results_[1],connecth[ii],1);
            acou_src_.Add(eqnNr,Result[ii]);
          }
        } catch (Exception& ex ) {
          RETHROW_EXCEPTION(ex, pdename_ << ": " << std:: endl 
              << "Could not get equation number.");
        }
        delete linear_load;
      }
      //          }
    }
  }
}


void FluidMechPDE::CalcMechCouplingRHS( StdVector<Elem*> * couplingElems,
    StdVector<UInt> & couplingNodes,
    Vector<Double>& elemCouplingSols,
    UInt couplingdof )
{
  Double density = 0.0;
  Double dynamicViscosity = 0.0;
  Double kinematicViscosity = 0.0;
  Double sign = 0.0;
  // TODO: Check if this is still needed
  // Double auxPres;
  Integer matIndex = -1;
  Elem * ptVolElem = NULL;
  Matrix<Double> ptCoord, ptCoordVol;
  Vector<Double> presSol, normal;
  Matrix<Double> veloSol;

  elemCouplingSols.Init();
  for (UInt actElem=0; actElem<couplingElems->GetSize(); actElem++)
  {
    // Perform cast from volume element to surface element, since
    // fluidMech-Mech coupling makes only sense on surface elements
    SurfElem * actCoupleElem =
      dynamic_cast<SurfElem*> ((*couplingElems)[actElem]);

    if (actCoupleElem == NULL)
    {
      EXCEPTION( "No elements found for coupling!");
    }

    StdVector<UInt> & connecth = actCoupleElem->connect;
    ptgrid_->GetElemNodesCoord( ptCoord, connecth, true );

    // Try to find according region for first neighbouring volume
    // element of the surface element
    matIndex = subdoms_.Find(actCoupleElem->ptVolElem1->regionId);

    // If first volume element does not belong to fluidMech PDE, try the
    // second one
    if ( matIndex == -1 )
    {
      matIndex = subdoms_.Find(actCoupleElem->ptVolElem2->regionId);
      ptVolElem = actCoupleElem->ptVolElem2;
      sign = actCoupleElem->normalSign;

    } else {
      ptVolElem = actCoupleElem->ptVolElem1;
      sign = -1.0 * actCoupleElem->normalSign;
    }

    if ( matIndex == -1)
    {
      EXCEPTION("FluidMechPDE::CalcMechCouplingRHS: The two volume element neighbours of surface element Nr.actCoupleElem->elemNum do not belong to my regions!");
    }

    // Assign correct density, dynamic viscosity and kinematic viscosity
    materials_[subdoms_[matIndex]]->GetScalar(density, DENSITY, Global::REAL);
    materials_[subdoms_[matIndex]]->GetScalar(dynamicViscosity, DYNAMIC_VISCOSITY, Global::REAL);
    materials_[subdoms_[matIndex]]->GetScalar(kinematicViscosity, KINEMATIC_VISCOSITY, Global::REAL);

    // the normal vector points outwards of the MECHANICAL domain
    // (see. Kaltenbacher, "Num. Sim. of Mechatr. Act. & Sens." chapter 8.2)
    ptgrid_->CalcSurfNormal(normal, *actCoupleElem);
    normal *= sign;

    //transform the type
    SubTensorType type;
    String2Enum(subType_,type);
    if (type == PLANE)
      type = PLANE_STRAIN;

    std::auto_ptr<FluidMechShearStress<Double> > shearStress =
      std::auto_ptr<FluidMechShearStress<Double> >(new FluidMechShearStress<Double>(materials_[subdoms_[matIndex]],type) );

    ElemList temp1( ptgrid_ );
    temp1.SetElement( actCoupleElem );
    EntityIterator surfIter = temp1.GetIterator();

    ElemList temp2( ptgrid_ );
    temp2.SetElement( ptVolElem );
    EntityIterator volIter = temp2.GetIterator();

    // force has to be added on RHS with negative sign
    Matrix<Double> pressureForce;
    Matrix<Double> viscStress;

    //fluidMech-pressure
    StdVector<UInt> & connecthVol = ptVolElem->connect;
    GetPresSolVecOfElement(presSol, connecthVol);
    shearStress->SetActElemPresSol(presSol);
    ptgrid_->GetElemNodesCoord( ptCoordVol, connecthVol, true );
    GetVeloSolMatOfElement(veloSol, connecthVol);
    shearStress->SetActElemSol(veloSol);

    shearStress->CalcPressureForce(pressureForce, volIter, surfIter, normal );
    shearStress->CalcShearStressVec(viscStress, volIter, surfIter, normal );

    //***********************************************************************************************//

    for (UInt actNode=0; actNode < connecth.GetSize(); actNode++)
    {
      UInt nodePos = 0;
      while (connecth[actNode] != couplingNodes[nodePos] \
                     && nodePos < couplingNodes.GetSize())
      {
        nodePos++;
      }
      const UInt nodePos_couplingdof = nodePos * couplingdof;
      for (UInt actDof = 0; actDof < couplingdof; ++actDof)
      {
        elemCouplingSols[nodePos_couplingdof +actDof] += \
          pressureForce[actNode][actDof] + viscStress[actNode][actDof];
      }
    }
  }//end loop over elem
}

void FluidMechPDE::CalcAcouSurfSourceCouplingRHS( StdVector<Elem*> * couplingElems,
    StdVector<UInt> & couplingNodes,
    Vector<Double>& elemCouplingSols,
    UInt couplingdof ) {
  Double density = 0.0;
  Double sign = 0.0;
  Integer matIndex = -1;
  Matrix<Double> ptCoord, elemMat;
  Vector<Double> normal, actElemAcouRHS;
  Vector<Double> veloDerivSol, nveloDerivSol;

  elemCouplingSols.Init();
  for (UInt actElem=0; actElem<couplingElems->GetSize(); actElem++) {

    // Perform cast from volume element to surface element, since
    // acouSurfSource coupling makes only sense on surface elements
    SurfElem * actCoupleElem =
      dynamic_cast<SurfElem*> ((*couplingElems)[actElem]);

    if (actCoupleElem == NULL) {
      EXCEPTION( "No elements found for coupling!");
    }

    // TODO: Check if this is still needed
    // BaseFE * ptElem = actCoupleElem->ptElem;

    StdVector<UInt> & connecth = actCoupleElem->connect;
    ptgrid_->GetElemNodesCoord( ptCoord, connecth, true );

    //Double surfArea;
    //surfArea=ptElem->CalcVolume(ptCoord, isaxi_);

    // Try to find according region for first neighbouring volume
    // element of the surface element
    matIndex = subdoms_.Find(actCoupleElem->ptVolElem1->regionId);

    // If first volume element does not belong to fluidMech PDE, try the
    // second one
    if ( matIndex == -1 ) {
      matIndex = subdoms_.Find(actCoupleElem->ptVolElem2->regionId);
      sign = -1.0*actCoupleElem->normalSign;

    } else {
      sign = actCoupleElem->normalSign;
    }

    if ( matIndex == -1) {
      EXCEPTION("FluidMechPDE::CalcAcouSurfSourceCouplingRHS: The two volume element neighbours of surface element do not belong to my regions!");
    }

    materials_[subdoms_[matIndex]]->GetScalar(density,DENSITY,Global::REAL);

    std::auto_ptr<BaseForm> bilinear_mass = std::auto_ptr<BaseForm>(new MassInt(1.0,isaxi_,coordUp_ ));

    ElemList temp( ptgrid_ );
    temp.SetElement( actCoupleElem );
    EntityIterator it = temp.GetIterator();
    bilinear_mass->CalcElementMatrix(elemMat, it, it);

    GetVeloDeriv1SolVecOfElement(veloDerivSol, it);
    nveloDerivSol.Resize(connecth.GetSize());

    // the normal vector points outwards of the MECHANICAL domain
    // (see. Kaltenbacher, "Num. Sim. of Mechatr. Act. & Sens." chapter 8.2)
    ptgrid_->CalcSurfNormal(normal, *actCoupleElem);
    normal *= sign;
    for (UInt actNode=0; actNode < connecth.GetSize(); actNode++)
      for (UInt actDof=0; actDof<dim_; actDof++)
        nveloDerivSol[actNode] += veloDerivSol[actDof + actNode*dim_] * normal[actDof];


    actElemAcouRHS = elemMat * nveloDerivSol;

    // force has to be added on RHS with negative sign
    actElemAcouRHS *= (-1.0 * density);



    for (UInt actNode=0; actNode<connecth.GetSize(); actNode++) {
      UInt nodePos = 0;

      while(connecth[actNode] != couplingNodes[nodePos] &&
          nodePos < couplingNodes.GetSize()) {
        nodePos++;
      }
      elemCouplingSols[nodePos] += actElemAcouRHS[actNode];
    }
  }//end loop over elem
}

bool FluidMechPDE::HasOutput(SolutionType output)
{
  if ( output == FLUIDMECH_VELOCITY || output == FLUIDMECH_PRESSURE \
      || output == FLUIDMECH_FORCE \
      || output == ACOU_RHSVAL || output == ACOUSURF_RHSVAL )
  {
    return true;
  }
  return false;
}



// ======================================================
// TIME STEPPING SECTION
// ======================================================


void FluidMechPDE::InitTimeStepping()
{
  // timestepping formulation
  TS_alg_ = new Bdf2( algsys_);
}

// ***********************************************************************
//   Obtain information on desired output quantities from parameter file
// ***********************************************************************
void FluidMechPDE::ReadSpecialResults() {

  //gridSol_ = NULL;
  gridSol_->SetNumSolutions(1);
  gridSol_->SetNumNodes(numPDENodes_);
  gridSol_->SetSolutionTypeName(FLUIDMECH_VELOCITY,GRID_VELOCITY);
  gridSol_->SetNumDofs(dim_);
  gridSol_->SetResults( results_ );
  gridSol_->SetPtrEQNData( eqnMap_.get(), ptgrid_);
  gridSol_->SetRegions( subdoms_ );
  gridSol_->Init();

  if(saveAcouSrc_ == true){
    acou_src_.SetNumSolutions(1);
    acou_src_.SetSolutionTypeName(FLUIDMECH_PRESSURE,ACOU_RHSVAL);
    acou_src_.SetNumNodes(numPDENodes_);
    acou_src_.SetNumDofs(1);   //dim_);
    acou_src_.SetResults(results_);
    acou_src_.SetPtrEQNData( eqnMap_.get(), ptgrid_);
    acou_src_.SetRegions( subdoms_ );
    acou_src_.Init();
  }

}


// *********************************************
//   Check what type of damping should be used
// *********************************************
void FluidMechPDE::ReadDampingInformation( ) {

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

    if( regionNodes.GetSize() > 0 ) {
      Info->PrintF( pdename_, "Damping in following region(s)\n" );
    }

    for (UInt k = 0; k < regionNodes.GetSize(); k++) {
      regionNodes[k]->GetValue( "name", actRegionName );
      regionNodes[k]->GetValue( "dampingId", actDampingId );
      if( actDampingId == "" )
        continue;

      actRegionId = ptgrid_->GetRegion().Parse(actRegionName);

      // Check actDampingId was already registerd
      if( idDampType.count( actDampingId ) == 0 ) {
        EXCEPTION( "Damping with id '" << actDampingId
                   << "' was not defined in 'dampingList'" );
      }

      dampingList_[actRegionId] = idDampType[actDampingId];

      // Log to info file
      std::string dampString;
      Enum2String( dampingList_[actRegionId], dampString );
      Info->PrintF( pdename_, " %s: %s\n", actRegionName.c_str(),
                    dampString.c_str() );
    }
    Info->PrintF( pdename_, "\n" );
  }


void FluidMechPDE::DefineAvailResults() {

  //the size of the matrix with the stabilization parameters is currently
  //Num2DElements/Num3DElements long.
  //Bec currently i have no idea how to access local elem numbers in the integrator
  stabilParams_.Resize(ptgrid_->GetNumElems()+1,4);
  stabilParams_.InitValue(-1.0);    	    	

  // Check for subType
  StdVector<std::string> dispDofNames, stressDofNames;

  if( subType_ == "3d" ) {
    dispDofNames = "x", "y", "z";
    stressDofNames = "xx", "yy", "zz", "yz", "xz", "xy";

  } else if( subType_ == "plane" ) {
    dispDofNames = "x", "y";
    stressDofNames = "xx", "yy", "xy";

  } else if( subType_ == "axi" ) {
    dispDofNames = "r", "z";
    stressDofNames = "rr", "zz", "rz", "phiphi";

  }

  // === FLUIDMECHANIC VELOCITY ===
  shared_ptr<ResultInfo> vel(new ResultInfo);
  vel->resultType = FLUIDMECH_VELOCITY;
  vel->dofNames = dispDofNames;
  vel->unit = "m/s";
  vel->entryType = ResultInfo::VECTOR;

  // check if problem is lagrange, legendre or taylorHood
  if (approxType_ == "lagrange")
  {	
    shared_ptr<AnsatzFct> fct(new LagrangeFct);
    vel->fctType = fct;
    vel->definedOn = ResultInfo::NODE;
  } else if ( approxType_ == "legendre" ) {
    shared_ptr<LegendreFct> fct(new LegendreFct);
    UInt order = myParam_->Get("order")->As<UInt>();
    fct->SetIsoOrder( order );
    vel->fctType = fct;
    vel->definedOn = ResultInfo::PFEM;
  } else
    EXCEPTION( "Wrong approximation type" );	

  results_.Push_back( vel );
  availResults_.insert( vel);

  // === FLUIDMECHANIC PRESSURE ===
  shared_ptr<ResultInfo> pres(new ResultInfo);
  pres->resultType = FLUIDMECH_PRESSURE;
  pres->dofNames = "";
  pres->unit = "N/m^2";
  pres->entryType = ResultInfo::SCALAR;

  // check if problem is lagrange or legendre
  if (approxType_ == "lagrange")
  {
    shared_ptr<AnsatzFct> fct(new LagrangeFct);
    pres->fctType = fct;
    pres->definedOn = ResultInfo::NODE;
  } else if (approxType_ == "legendre") {
    // define Legendre type
    shared_ptr<LegendreFct> fct(new LegendreFct);
    UInt order = myParam_->Get("order")->As<UInt>();
    fct->SetIsoOrder( order );
    pres->fctType = fct;
    pres->definedOn = ResultInfo::PFEM;
  } else {
    EXCEPTION( "Wrong approximation type" );	
  }

  results_.Push_back( pres );
  availResults_.insert( pres);

  //      // === FLUIDMECHANIC STRESS ===
  //      shared_ptr<ResultInfo> stress(new ResultInfo);
  //      stress->resultType = FLUIDMECH_STRESS;
  //      stress->dofNames = stressDofNames;
  //      stress->unit =  "N/m^2";
  //      stress->entryType = ResultInfo::TENSOR;
  //      stress->definedOn = ResultInfo::ELEMENT;
  //      stress->fctType = shared_ptr<ConstFct>(new ConstFct() );
  //      availResults_.insert( stress );
  //
  //      // === FLUIDMECHANIC STRAIN RATE===
  //      shared_ptr<ResultInfo> strainRate(new ResultInfo);
  //      strainRate->resultType = FLUIDMECH_STRAINRATE;
  //      strainRate->dofNames = stressDofNames;
  //      strainRate->unit =  "";
  //      strainRate->entryType = ResultInfo::TENSOR;
  //      strainRate->definedOn = ResultInfo::ELEMENT;
  //      strainRate->fctType = shared_ptr<ConstFct>(new ConstFct() );
  //      availResults_.insert( strainRate );
  //
  //      // === FLUIDMECHANIC ENERGY ===
  //      shared_ptr<ResultInfo> energy(new ResultInfo);
  //      energy->resultType = FLUIDMECH_ENERGY;
  //      energy->dofNames = "";
  //      energy->unit = "Ws";
  //      energy->entryType = ResultInfo::SCALAR;
  //      energy->definedOn = ResultInfo::REGION;
  //      energy->fctType = shared_ptr<ConstFct>(new ConstFct() );
  //      availResults_.insert( energy );
  //
  //      // === ACOUSTIC SOURCES (Lighthill) ===
  //      shared_ptr<ResultInfo> acouSrc(new ResultInfo);
  //      acouSrc->resultType = ACOU_RHS_LOAD;
  //      acouSrc->dofNames = "";
  //      acouSrc->unit = "";
  //      acouSrc->entryType = ResultInfo::SCALAR;
  //      acouSrc->definedOn = ResultInfo::NODE;
  //      acouSrc->fctType = shared_ptr<ConstFct>(new ConstFct() );
  //      availResults_.insert( acouSrc );
  //
  //
  //      // === STABILIZATION PARAMETERS ===
  //      shared_ptr<ResultInfo> stabilParam(new ResultInfo);
  //      stabilParam->resultType = FLUIDMECH_STABILPARAM;
  //      stabilParam->dofNames = "";
  //      stabilParam->unit = "";
  //      stabilParam->entryType = ResultInfo::VECTOR;
  //      stabilParam->definedOn = ResultInfo::ELEMENT;
  //      stabilParam->fctType = shared_ptr<ConstFct>(new ConstFct() );
  //      availResults_.insert( stabilParam );

  // === OPT_RESULT_1/2/3 ===
  // this is added via the optimization stuff in DesignSpace.
}

void FluidMechPDE::CalcResults( shared_ptr<BaseResult> result ) {

  switch (result->GetResultInfo()->resultType ) {

    case FLUIDMECH_VELOCITY:
      if( isComplex_ ) {
        ExtractResult<Complex>( result, sol_ );
      } else {
        ExtractResult<Double>( result, sol_ );
      }
      break;

    case FLUIDMECH_PRESSURE:
      if( isComplex_ ) {
        ExtractResult<Complex>( result, sol_ );
      } else {
        ExtractResult<Double>( result, sol_ );
      }
      break;

    case ACOU_RHS_LOAD:
      if( isComplex_ ) {
        ExtractRhsResult<Complex>( result, results_[0] );
      } else {
        ExtractRhsResult<Double>( result, results_[0] );
      }
      break;

    case FLUIDMECH_STRESS:
      if( isComplex_ ) {
        CalcStresses<Complex>( result );
      } else {
        CalcStresses<Double>( result );
      }
      break;

    case FLUIDMECH_STRAINRATE:
      if( isComplex_ ) {
        CalcStrainRates<Complex>( result );
      } else {
        CalcStrainRates<Double>( result );
      }
      break;

    case FLUIDMECH_ENERGY:
      if( isComplex_ ) {
        CalcEnergy<Complex>( result );
      } else {
        CalcEnergy<Double>( result );
      }
      break;

      // the actual case is given in the result info in result
    case OPT_RESULT_1:
    case OPT_RESULT_2:
    case OPT_RESULT_3:
      // design should work, this is checked in AvailabeResults()
      domain->GetErsatzMaterial()->ExtractResults(result, isComplex_);
      break;


    default:
      WARN( "Resulttype not computable by mechanic PDE" );
  }
}


template <class TYPE>
void FluidMechPDE::CalcStresses( shared_ptr<BaseResult> res ) {

  // TODO: Why the hell is this copy and paste from MechPDE??

  //get the correct bilinear form
  Vector<Double> intPoint;

  //transform the type
  SubTensorType type;
  String2Enum(subType_,type);

  Vector<TYPE> elemStress;
  elemStress.Resize(stressDim_);
  elemStress.Init();

  Result<TYPE> &  actRes =
    dynamic_cast<Result<TYPE>&>(*res);
  EntityIterator it = actRes.GetEntityList()->GetIterator();

  Vector<TYPE> & actVal = actRes.GetVector();
  actVal.Resize( actRes.GetEntityList()->GetSize() * stressDim_ );

  // Fetch material: As we assume, that all elements belong to
  // one and the same region, we simply take the subdomain of the first
  // element
  it.Begin();
  BaseMaterial* actSDMat = materials_[it.GetElem()->regionId];
  MechStressStrain<TYPE> *stress = new MechStressStrain<TYPE>(actSDMat,type);

  // loop over elements
  for ( it.Begin(); !it.IsEnd(); it++ ) {
    it.GetElem()->ptElem->GetCoordMidPoint(intPoint);

    //set element solution
    Matrix<TYPE> elSol;
    sol_->GetElemSolutionAsMatrix(elSol, it);
    stress->SetActElemSol(elSol);

    //calculates the stress
    stress->SetIntPoint(intPoint);
    stress->CalcStressVec(elemStress,1,it);
    stress->UnsetIntPoint();
    for( UInt iDof = 0; iDof < stressDim_; iDof++ ) {
      actVal[it.GetPos()*stressDim_ + iDof] = elemStress[iDof];
    }
  }

  delete stress;
}

template <class TYPE>
void FluidMechPDE::CalcStrainRates( shared_ptr<BaseResult> res ) {
  //transform the subType of the pde
  SubTensorType type;
  String2Enum(subType_,type);

  Vector<Double> intPoint;
  Vector<TYPE> elemStrainRate;
  elemStrainRate.Resize(stressDim_);
  elemStrainRate.Init();

  Result<TYPE> &  actRes =
    dynamic_cast<Result<TYPE>&>(*res);
  EntityIterator it = actRes.GetEntityList()->GetIterator();

  Vector<TYPE> & actVal = actRes.GetVector();
  actVal.Resize( actRes.GetEntityList()->GetSize() * stressDim_ );

  // Fetch material: As we assume, that all elements belong to
  // one and the same region, we simply take the subdomain of the first
  // element
  it.Begin();
  BaseMaterial* actSDMat = materials_[it.GetElem()->regionId];
  MechStressStrain<TYPE> *strainRate = new MechStressStrain<TYPE>(actSDMat,type);

  // loop over elements
  for ( it.Begin(); !it.IsEnd(); it++ ) {
    it.GetElem()->ptElem->GetCoordMidPoint(intPoint);

    StdVector<UInt> connecth = it.GetElem()->connect;

    //set element solution
    Matrix<TYPE> elSol;
    GetVeloSolMatOfElement(elSol, connecth);

    strainRate->SetActElemSol(elSol);

    //calculates the element strainRate
    strainRate->SetIntPoint(intPoint);
    strainRate->CalcStrainVec(elemStrainRate,1,it);
    strainRate->UnsetIntPoint();
    for( UInt iDof = 0; iDof < stressDim_; iDof++ ) {
      actVal[it.GetPos()*stressDim_ + iDof] = elemStrainRate[iDof];
    }
  }

  delete strainRate;
}




// ======================================================
// POSTPROCESSING
// ======================================================
template<class TYPE>
void FluidMechPDE::ComputeVolDefSurf( shared_ptr<BaseResult> vals ) {


  // convert result and get region iterator
  Result<TYPE> &  actSol =
    dynamic_cast<Result<TYPE>&>(*vals);
  EntityIterator regionIt = actSol.GetEntityList()->GetIterator();

  // resize vector
  Vector<TYPE> & actVal = actSol.GetVector();
  actVal.Resize( actSol.GetEntityList()->GetSize() );

  // get dof for direction
  std::string dir = volAboveDefSurfDir_[actSol.GetEntityList()];
  Integer dof = results_[0]->dofNames.Find( dir );
  if( dof < 0 ) {
    EXCEPTION( "ComputeVolDefSurf: dof = '" << dir << "' is not "
        << "allowed! Must be one of 'x', 'y' or 'z'" );
  }

  Vector<TYPE> elemDisp, elemDispDof;
  TYPE actVolume;

  // Loop over regions
  for( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ ) {
    SurfElemList actSDList(ptgrid_ );
    actSDList.SetRegion( regionIt.GetRegion() );
    EntityIterator elemIt = actSDList.GetIterator();

    actVolume = 0.0;
    // Loop over elements
    for( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {

      BaseFE * ptSurfEl = elemIt.GetSurfElem()->ptElem;
      StdVector<UInt> connecth = elemIt.GetSurfElem()->connect;

      Matrix<Double> ptSurfCoord;
      ptgrid_->GetElemNodesCoord( ptSurfCoord, connecth, false );

      GetSolVecOfElement( elemDisp, elemIt, results_[0] );
      elemDispDof.Resize(connecth.GetSize() );
      UInt actEntry = dof;
      for( UInt i = 0; i < connecth.GetSize(); i++ ) {
        elemDispDof[i] = elemDisp[actEntry];
        actEntry += dim_;
      }
      actVolume += ComputeVolElem<TYPE>( ptSurfEl,
          ptSurfCoord, elemDispDof );

    }
    actVal[regionIt.GetPos()] = actVolume;
  }

}

template<class TYPE>
TYPE FluidMechPDE::ComputeVolElem(BaseFE * surfEl, Matrix<Double>& surfCoord,
    Vector<TYPE> disp) {


  TYPE elemVol, averageDis;
  UInt nrSurfNodes = surfEl->GetNumNodes();


  //compute average displacedment
  averageDis = 0;
  for (UInt i=0; i<nrSurfNodes; i++) {
    averageDis += disp[i];
  }
  averageDis /= (Double)nrSurfNodes;

  //compute the deformed volume
  elemVol = averageDis * surfEl->CalcVolume(surfCoord,isaxi_);

  return elemVol;
}



void FluidMechPDE::ReadPressureLoads() {

  // try to get bcsAndLoads node
  PtrParamNode bcNode = myParam_->Get("bcsAndLoads", ParamNode::PASS);
  if( !bcNode )
    return;
  ParamNodeList presNodes = bcNode->GetList("pressure");

  // iterate over all pressure definitions
  std::string name, value, phase;
  for( UInt i = 0; i < presNodes.GetSize(); i++ ) {

    presNodes[i]->GetValue("name", name );
    presNodes[i]->GetValue("value", value );
    presNodes[i]->GetValue("phase", phase );

    pressSurf_.Push_back(ptgrid_->GetEntityList( EntityList::SURF_ELEM_LIST,
          name, EntityList::REGION ) );
    pressVals_.Push_back( value );
    pressPhase_.Push_back( phase );

  }
}


  template <class TYPE>
void FluidMechPDE::CalcEnergy( shared_ptr<BaseResult> vals )
{

  Matrix<Double> elemmat;
  Vector<TYPE> help, elvel;
  TYPE energy ;

  Result<TYPE> &  actSol =
    dynamic_cast<Result<TYPE>&>(*vals);

  // resize vector
  Vector<TYPE> & actVal = actSol.GetVector();
  actVal.Resize( actSol.GetEntityList()->GetSize() );

  // loop over regions
  EntityIterator regionIt = actSol.GetEntityList()->GetIterator();
  for( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ ) {

    //get material
    BaseMaterial* actSDMat = materials_[regionIt.GetRegion()];

    // get bilinear stiffness integrator
    BaseForm * bilinear_stiff = GetStiffIntegrator(actSDMat,
        regionIt.GetRegion() );
    ElemList actSDList(ptgrid_ );
    actSDList.SetRegion( regionIt.GetRegion() );
    EntityIterator it = actSDList.GetIterator();
    energy = 0.0;

    // Loop over elements
    for ( it.Begin(); !it.IsEnd(); it++ ) {
      bilinear_stiff->SetAnsatzFct( results_[0]->fctType );
      bilinear_stiff->CalcElementMatrix(elemmat,it,it);
      sol_->GetElemSolution(elvel, it);
      help = elemmat * elvel;
      energy += ( help * elvel) * 0.5;
    }
    actVal[regionIt.GetPos()] = energy;
    delete bilinear_stiff;
  }
}

void FluidMechPDE::InitStabParams(){

  if( stabilizationType_ != "none"){
    //the size of the matrix with the stabilization parameters is currently
    //Num2DElements/Num3DElements long.
    //Bec currently i have no idea how to access local elem numbers in the integrator
    //std::cerr << "init\n";
    for(UInt i=0; i<ptgrid_->GetNumElems(); i++)
      stabilParams_[i][3]=(-1.0);
  }
}

void FluidMechPDE::PrintStabParams(){

  if( stabilizationType_ != "none"){
    LOG_DBG(fluidmechpde) << "Stabilization parameters\n"
      << stabilParams_ << std::endl;
  }
}


void FluidMechPDE::RegisterFunctions() {

  typedef FctPointer<FluidMechPDE> FCPT;
  StdVector<ArgList> a;
  StdVector<FCPT*> pt;
  StdVector<std::string> name;

  // --- ReadRegionLoad ---
  a.Push_back();
  a.Last().RegisterParam( "name", ArgList::STRING );
  a.Last().RegisterParam( "value", ArgList::STRING );
  a.Last().RegisterParam( "dof", ArgList::STRING );
  a.Last().RegisterParam( "coordSysId", ArgList::STRING );
  a.Last().RegisterParam( "type", ArgList::STRING );
  pt.Push_back( new FCPT ( this, &FluidMechPDE::ReadRegionLoads ) );
  name.Push_back( "setRegionLoad" );


  // Now register all functions with scripting
  for (UInt i = 0; i < pt.GetSize(); i++ ) {
    Script_RegisterFct(name[i], pt[i], a[i] );
  }

}


} // end namespace CoupledField
