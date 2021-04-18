// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

//#include <boost/smart_ptr/shared_ptr.hpp>
#include "LinFlowMechCoupling.hh"

#include "PDE/SinglePDE.hh"
#include "PDE/PerturbedFlowPDE.hh"
#include "PDE/MechPDE.hh"
#include "CoupledPDE/BasePairCoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "Materials/BaseMaterial.hh"
#include "Driver/FormsContexts.hh"
#include "Driver/Assemble.hh"

#include "Domain/Mesh/NcInterfaces/MortarInterface.hh"

// include fespaces
#include "FeBasis/H1/H1Elems.hh"

// new integrator concept
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/BiLinForms/BiLinearForm.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"
#include "Forms/Operators/SurfaceOperators.hh"
#include "Forms/Operators/SurfaceNormalStressOperator.hh"
#include "Domain/CoefFunction/CoefXpr.hh"

namespace CoupledField {



// ***************
//   Constructor
// ***************
LinFlowMechCoupling::LinFlowMechCoupling( SinglePDE *pde1, SinglePDE *pde2,
    PtrParamNode paramNode,
    PtrParamNode infoNode,
    shared_ptr<SimState> simState,
    Domain* domain)
: BasePairCoupling( pde1, pde2, paramNode, infoNode, simState, domain ),
  lmOrderSameAsVel_(true), IsLagrangeMultiplierMethod_(false){

  couplingName_ = "linFlowMechDirect";
  materialClass_ = FLOW;

  formulation_ = NO_SOLUTION_TYPE;

  // determine subtype from mechanic pde
  pde2_->GetParamNode()->GetValue( "subType", subType_ );

  nonLin_ = false;

  // Initialize nonlinearities
  InitNonLin();
  if( paramNode->Has("IsLagrangeMultiplierMethod") ) {
    IsLagrangeMultiplierMethod_ =  paramNode->Get("IsLagrangeMultiplierMethod")->As<bool>();
    if( paramNode->Has("lmOrderSameAsVel") ) {
      lmOrderSameAsVel_ =  paramNode->Get("lmOrderSameAsVel")->As<bool>();
    }
  }
}


// **************
//   Destructor
// **************
LinFlowMechCoupling::~LinFlowMechCoupling() {
}


// *********************
//   DefineIntegrators
// *********************
void LinFlowMechCoupling::DefineIntegrators() {

  MathParser * mp = domain_->GetMathParser();

  // fixed in Domain::CreateDirectCoupledPDEs: pde1 is fluidMechLin and
  // pde2 is mechanic
  //
  shared_ptr<BaseFeFunction> velFct  = pde1_->GetFeFunction(FLUIDMECH_VELOCITY);

  shared_ptr<BaseFeFunction> dispFct = pde2_->GetFeFunction(MECH_DISPLACEMENT);

  // Create coefficient functions
  std::map< RegionIdType, PtrCoefFct > oneFuncs;

  shared_ptr<FeSpace> dispSpace = dispFct->GetFeSpace();
  shared_ptr<FeSpace> velSpace = velFct->GetFeSpace();

  std::set< RegionIdType > flowRegions;
  std::map<RegionIdType, BaseMaterial*> flowMaterials;
  flowMaterials = pde1_->GetMaterialData();
  std::map<RegionIdType, BaseMaterial*>::iterator it, end;
  it = flowMaterials.begin();
  end = flowMaterials.end();

  for ( ; it != end; it++ ) {
    RegionIdType volRegId = it->first;
    flowRegions.insert(volRegId);

    oneFuncs[volRegId] = CoefFunction::Generate(mp, Global::REAL,lexical_cast<std::string>(1.0));
    shared_ptr<ElemList> actSDList( new ElemList( ptGrid_ ) );
    actSDList->SetRegion( volRegId );
  }

  if (IsLagrangeMultiplierMethod_){
    shared_ptr<BaseFeFunction> presFct = pde1_->GetFeFunction(FLUIDMECH_PRESSURE);
    shared_ptr<BaseFeFunction> lagrangeMultFct = feFunctions_[LAGRANGE_MULT];
    shared_ptr<FeSpace> presSpace = presFct->GetFeSpace();
    shared_ptr<FeSpace> lagrangeMultSpace = lagrangeMultFct->GetFeSpace();
    //      go over all coupled surfaces
    for ( UInt actSD = 0, n = entityLists_.GetSize(); actSD < n; actSD++ ) {

      shared_ptr<SurfElemList> actSDList =dynamic_pointer_cast<SurfElemList>(entityLists_[actSD]);
      RegionIdType region = actSDList->GetRegion();

      velFct->AddEntityList(actSDList);
      presFct->AddEntityList(actSDList);
      dispFct->AddEntityList(actSDList);
      lagrangeMultFct->AddEntityList(actSDList);

      velSpace->SetRegionApproximation(region, "velPolyId", "velIntegId");
      presSpace->SetRegionApproximation(region, "presPolyId", "presIntegId");
      dispSpace->SetRegionApproximation(region, "default", "default");
      if(lmOrderSameAsVel_) {
        lagrangeMultSpace->SetRegionApproximation(region, "velPolyId", "velIntegId");
      }
      else {
        lagrangeMultSpace->SetRegionApproximation(region, "presPolyId", "velIntegId");
      }

      // Integrator being assembled into damping (first time deriv.) matrix; first part
      // of additional equation guaranteeing continuity of velocities
      DefineDampingIntegrators("LinFlowMechDampingLMVelCouplingInt", dispFct, lagrangeMultFct,
          actSDList, oneFuncs, flowRegions);

      // These integrators gets assembled into the stiffness matrix of the mechanic PDE,
      // LinFlowPDE and the additional equation guaranteeing continuity of velocities
      // equation for continuity of velocities)
      DefineStiffnessIntegrators("LinFlowMechStiff", dispFct, velFct, lagrangeMultFct,
          actSDList, oneFuncs, oneFuncs, oneFuncs, flowRegions);
    }
  }
  else{ //Non conforming integrations for Nitsche coupling
    ParamNodeList ncList = myParam_->Get("ncInterfaceList", ParamNode::PASS)->GetList("ncInterface");
    for(UInt i = 0; i < ncList.GetSize(); i++){
      //      std::cout<< ncList[i]->Get("ncInterface/name")->As<std::string>() <<std::endl;
      std::string ncName = ncList[i]->Get("name")->As<std::string>();
      shared_ptr<BaseNcInterface> ncIf = ptGrid_->GetNcInterface(ptGrid_->GetNcInterfaceId(ncName));
      if (!ncIf)
      {
        EXCEPTION("No interface with the name '" << ncName << "' found!");
      }
      std::string formulation = ncList[i]->Get("formulation")->As<std::string>();
      if (formulation == "Nitsche"){
        MortarInterface * nitscheIf = dynamic_cast<MortarInterface*>(ncIf.get());
        assert(nitscheIf);
        //in case of Nitsche coupling edge/face information is required
        this->ptGrid_->MapEdges();
        this->ptGrid_->MapFaces();
        RegionIdType volMasterId = nitscheIf->GetMasterVolRegion();
        RegionIdType volSlaveId = nitscheIf->GetSlaveVolRegion();

        //we set here the penalty factor
        Double beta = ncList[i]->Get("nitscheFactor")->As<Double>();
        // create new entity list
        shared_ptr<ElemList> actNCSDList = ncIf->GetElemList();
        // in case of mechanical PDE, we need the material tensor
        shared_ptr<CoefFunction > coefMech;
        SubTensorType tensorType = NO_TENSOR;

        if ( subType_ == "3d" )
          tensorType = FULL;
        else if ( subType_ == "axi" )
          tensorType = AXI;
        else if ( subType_ == "planeStrain" )
          tensorType = PLANE_STRAIN;
        else if ( subType_ == "planeStress" )
          tensorType = PLANE_STRESS;

        //get correct scaling of penalty term
        StdVector<Vector<Double> > points(1);
        Vector<Double> p1(dim_);
        p1.Init();
        points[0]= p1;
        bool isMaterialComplex_( false );

        // extract stiffness tensor from mech pde for stress tensor
        std::map<RegionIdType, BaseMaterial*>  mechMat;
        mechMat = pde2_->GetMaterialData();

        // Here the mechanic PDE must be considered as MASTER
        if ( isMaterialComplex_ ) {
          coefMech = mechMat[volMasterId]->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, tensorType, Global::COMPLEX);
          StdVector< Matrix<Complex> > mat;
          coefMech->GetTensorValuesAtCoords(points, mat, this->ptGrid_);
        }
        else {
          coefMech = mechMat[volMasterId]->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, tensorType, Global::REAL);
          StdVector< Matrix<Double> > mat;
          coefMech->GetTensorValuesAtCoords(points, mat, this->ptGrid_);
        }

        //get correct scaling of penalty term
        Double Beta = 0.0;

        // extract material from LinFlow pde
        shared_ptr<CoefFunction > shearViscosity;
        shared_ptr<CoefFunction > bulkViscosity;
        shearViscosity = flowMaterials[volSlaveId]->GetScalCoefFnc(FLUID_DYNAMIC_VISCOSITY, Global::REAL);
        bulkViscosity = flowMaterials[volSlaveId]->GetScalCoefFnc(FLUID_BULK_VISCOSITY, Global::REAL);
        PtrCoefFct density_f = flowMaterials[volSlaveId]->GetScalCoefFnc(DENSITY, Global::REAL);
        PtrCoefFct density_s = mechMat[volMasterId]->GetScalCoefFnc(DENSITY, Global::REAL);

        PtrCoefFct shearPlusBulkCoeff = CoefFunction::Generate(mp,Global::REAL,CoefXprBinOp(mp,shearViscosity,bulkViscosity, CoefXpr::OP_ADD)); //sum of the viscosities
        LocPointMapped map;
        Double shearPlusBulk = 0;
        shearPlusBulkCoeff->GetScalar(shearPlusBulk,map);
        if (shearPlusBulk == 0) {    //make sure viscosities are not zero
          shearPlusBulkCoeff = CoefFunction::Generate( mp, Global::REAL, "1");
        }

        // densratiosqr = sqrt(density_f/density_s)
        PtrCoefFct densratiosqr = CoefFunction::Generate(mp, Global::REAL,
            CoefXprUnaryOp(mp,
                CoefXprBinOp(mp,density_f ,density_s ,CoefXpr::OP_DIV)
                ,CoefXpr::OP_SQRT));

        //  betaCoeff = densratiosqr*(\lambda+\mu)
        PtrCoefFct betaCoeff = CoefFunction::Generate(mp,Global::REAL,CoefXprBinOp(mp,shearPlusBulkCoeff,densratiosqr, CoefXpr::OP_MULT));
        betaCoeff->GetScalar(Beta,map);

        // Create final Nitsche's factor
        beta *= Beta;

        PtrCoefFct coefOne = CoefFunction::Generate( mp, Global::REAL, "1");
        PtrCoefFct coefMinusOne = CoefFunction::Generate( mp, Global::REAL, "-1.0");
        // Note:
        // Define the bilinear forms for the Nitsche coupling consists of six terms
        // The terms are divided into four categories:
        // 1) The bilinear forms coupling mechanical displacement with LinFlow velocity
        // 2) The Penalty integrals coupling mechanical displacement with LinFlow velocity
        // 3) The bilinear forms coupling LinFlow velocity with mechanical displacement
        // 4) The Penalty integrals coupling LinFlow velocity with mechanical displacement
        // assume phi is displacement test function and u is displacement
        // psi is velocity test function and v LinFlow velocity

        BiLinearForm* int_PhiM_duM = NULL;   // Term 1: u'M_duM : phi * sigma[u].n
        BiLinearForm* penalty_PhiM_uM = NULL;// term 2  u'M_uM : beta*phi*du/dt
        BiLinearForm* penalty_PhiM_vS = NULL;// term 3 u'M_vS :-beta*phi*v
        BiLinearForm* int_PsiS_duM = NULL;   // Term 4: v'S_duM : Psi*sigma[u].n
        BiLinearForm* penalty_PsiS_vS = NULL;// Term 5: v'S_vS : beta*Psi.v
        BiLinearForm* penalty_PsiS_uM = NULL;// Term 6: v'S_uM : -beta*Psi.du/dt

        BaseBOperator * sNSOp1 = NULL;
        BaseBOperator * sNSOp2 = NULL;
        // ====================================================================
        //  PART ONE
        // ====================================================================
        // --------------------------
        // Term 1 : phi * sigma[u].n
        // --------------------------
        if (dim_ == 3)
        {
          sNSOp1 = new SurfaceNormalStressOperator<FeH1, 3, 3,Complex>(subType_ ,false);
          sNSOp1->SetCoefFunction(coefMech);
          int_PhiM_duM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 3, 3>(),
              sNSOp1, coefMinusOne, 1.0, BiLinearForm::MASTER_MASTER, false, true, false);
        }
        else if (dim_ == 2 && subType_ == "2.5d")
        {
          sNSOp1 = new SurfaceNormalStressOperator<FeH1, 2, 3,Complex>(subType_ ,false);
          sNSOp1->SetCoefFunction(coefMech);
          int_PhiM_duM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 3>(),
              sNSOp1, coefMinusOne, 1.0, BiLinearForm::MASTER_MASTER, false, true, false);
        }
        else
        {
          sNSOp1 = new SurfaceNormalStressOperator<FeH1, 2, 2,Complex>(subType_ ,false);
          sNSOp1->SetCoefFunction(coefMech);
          int_PhiM_duM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 2>(),
              sNSOp1, coefMinusOne, 1.0, BiLinearForm::MASTER_MASTER, false, true, false);
        }
        // ====================================================================
        //  PART TWO
        // ====================================================================
        // --------------------------
        // Term 2 : beta*phi*du/dt
        // --------------------------
        if (dim_ == 3)
        {
          penalty_PhiM_uM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 3, 3>(),
              new SurfaceIdentityOperator<FeH1, 3, 3>(),
              coefOne, beta, BiLinearForm::MASTER_MASTER , false, true, true);
        }
        else if (dim_ == 2 && subType_ == "2.5d")
        {
          penalty_PhiM_uM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 3>(),
              new SurfaceIdentityOperator<FeH1, 2, 3>(),
              coefOne, beta, BiLinearForm::MASTER_MASTER, false, true, true);
        }
        else
        {

          penalty_PhiM_uM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 2>(),
              new SurfaceIdentityOperator<FeH1, 2, 2>(),
              coefOne, beta, BiLinearForm::MASTER_MASTER, false, true, true);
        }
        // --------------------------
        // Term 3 : -beta*phi*v
        // --------------------------
        if (dim_ == 3)
        {
          penalty_PhiM_vS = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 3, 3>(),
              new SurfaceIdentityOperator<FeH1, 3, 3>(),
              coefMinusOne, beta,BiLinearForm::MASTER_SLAVE , false, true, true);
        }
        else if (dim_ == 2 && subType_ == "2.5d")
        {
          penalty_PhiM_vS = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 3>(),
              new SurfaceIdentityOperator<FeH1, 2, 3>(),
              coefMinusOne, beta, BiLinearForm::MASTER_SLAVE , false, true, true);
        }
        else
        {
          penalty_PhiM_vS = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 2>(),
              new SurfaceIdentityOperator<FeH1, 2, 2>(),
              coefMinusOne, beta, BiLinearForm::MASTER_SLAVE , false, true, true);
        }
        // ====================================================================
        //  PART THREE
        // ====================================================================
        // --------------------------
        // Term 4 : Psi*sigma[u].n
        // --------------------------
        if (dim_ == 3)
        {
          sNSOp2 = new SurfaceNormalStressOperator<FeH1, 3, 3,Complex>(subType_ ,false);
          sNSOp2->SetCoefFunction(coefMech);
          int_PsiS_duM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 3, 3>(),sNSOp2,
              coefOne, 1.0, BiLinearForm::SLAVE_MASTER, false, true, false);
        }
        else if (dim_ == 2 && subType_ == "2.5d")
        {
          sNSOp2 = new SurfaceNormalStressOperator<FeH1, 2, 3,Complex>(subType_ ,false);
          sNSOp2->SetCoefFunction(coefMech);
          int_PsiS_duM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 3>(),sNSOp2,
              coefOne, 1.0, BiLinearForm::SLAVE_MASTER, false, true, false);
        }
        else
        {
          sNSOp2 = new SurfaceNormalStressOperator<FeH1, 2, 2,Complex>(subType_ ,false);
          sNSOp2->SetCoefFunction(coefMech);
          int_PsiS_duM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 2>(),sNSOp2,
              coefOne, 1.0, BiLinearForm::SLAVE_MASTER, false, false, false);
        }
        // ====================================================================
        //  PART FOUR
        // ====================================================================
        // --------------------------
        // Term 5 : beta*Psi.v
        // --------------------------
        if (dim_ == 3)
        {
          penalty_PsiS_vS = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 3, 3>(),
              new SurfaceIdentityOperator<FeH1, 3, 3>(),
              coefOne, beta, BiLinearForm::SLAVE_SLAVE, false, true, true);
        }
        else if (dim_ == 2 && subType_ == "2.5d")
        {
          penalty_PsiS_vS = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 3>(),
              new SurfaceIdentityOperator<FeH1, 2, 3>(),
              coefOne, beta, BiLinearForm::SLAVE_SLAVE, false, true, true);
        }
        else
        {
          penalty_PsiS_vS = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 2>(),
              new SurfaceIdentityOperator<FeH1, 2, 2>(),
              coefOne, beta, BiLinearForm::SLAVE_SLAVE, false, false, true);
        }
        // --------------------------
        // Term 6 : -beta*Psi.du/dt
        // --------------------------

        if (dim_ == 3)
        {
          penalty_PsiS_uM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 3, 3>(),
              new SurfaceIdentityOperator<FeH1, 3, 3>(),
              coefMinusOne, beta, BiLinearForm::SLAVE_MASTER, false, true, true);
        }
        else if (dim_ == 2 && subType_ == "2.5d")
        {
          penalty_PsiS_uM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 3>(),
              new SurfaceIdentityOperator<FeH1, 2, 3>(),
              coefMinusOne, beta, BiLinearForm::SLAVE_MASTER, false, true, true);
        }
        else
        {
          penalty_PsiS_uM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 2>(),
              new SurfaceIdentityOperator<FeH1, 2, 2>(),
              coefMinusOne, beta, BiLinearForm::SLAVE_MASTER, false, false, true);
        }

        int_PhiM_duM->SetName("int_PhiM_duM");
        penalty_PhiM_uM->SetName("penalty_PhiM_uM");
        penalty_PhiM_vS->SetName("penalty_PhiM_vS");
        int_PsiS_duM->SetName("int_PsiS_duM");
        penalty_PsiS_vS->SetName("penalty_PsiS_vS");
        penalty_PsiS_uM->SetName("penalty_PsiS_uM");

        // define contexts for bilinear forms
        SurfaceBiLinFormContext* descr_PhiM_duM = new SurfaceBiLinFormContext(int_PhiM_duM, STIFFNESS, BiLinearForm::MASTER_MASTER);
        SurfaceBiLinFormContext* descr_PhiM_uM = new SurfaceBiLinFormContext(penalty_PhiM_uM, DAMPING , BiLinearForm::MASTER_MASTER);
        SurfaceBiLinFormContext* descr_PhiM_vS = new SurfaceBiLinFormContext(penalty_PhiM_vS, STIFFNESS, BiLinearForm::MASTER_SLAVE);
        SurfaceBiLinFormContext* descr_PsiS_duM = new SurfaceBiLinFormContext(int_PsiS_duM, STIFFNESS, BiLinearForm::SLAVE_MASTER);
        SurfaceBiLinFormContext* descr_PsiS_vS = new SurfaceBiLinFormContext(penalty_PsiS_vS, STIFFNESS, BiLinearForm::SLAVE_SLAVE);
        SurfaceBiLinFormContext* descr_PsiS_uM = new SurfaceBiLinFormContext(penalty_PsiS_uM, DAMPING, BiLinearForm::SLAVE_MASTER);

        descr_PhiM_duM->SetEntities(actNCSDList, actNCSDList);
         descr_PhiM_uM->SetEntities(actNCSDList, actNCSDList);
        descr_PhiM_vS->SetEntities(actNCSDList, actNCSDList);
        descr_PsiS_duM->SetEntities(actNCSDList, actNCSDList);
        descr_PsiS_vS->SetEntities(actNCSDList, actNCSDList);
        descr_PsiS_uM->SetEntities(actNCSDList, actNCSDList);

        descr_PhiM_duM->SetFeFunctions(dispFct, velFct);
        descr_PhiM_uM->SetFeFunctions(dispFct, velFct);
        descr_PhiM_vS->SetFeFunctions(dispFct, velFct);
        descr_PsiS_duM->SetFeFunctions(dispFct, velFct);
        descr_PsiS_vS->SetFeFunctions(dispFct, velFct);
        descr_PsiS_uM->SetFeFunctions(dispFct, velFct);

        assemble_->AddBiLinearForm(descr_PhiM_duM);
        assemble_->AddBiLinearForm(descr_PhiM_uM);
        assemble_->AddBiLinearForm(descr_PhiM_vS);
        assemble_->AddBiLinearForm(descr_PsiS_duM);
        assemble_->AddBiLinearForm(descr_PsiS_vS);
        assemble_->AddBiLinearForm(descr_PsiS_uM);
      }
      else if (formulation == "Mortar")
      {
        // we do nothing here
      } // end mortar
      else
      {
        EXCEPTION("Unknown formulation: '" << formulation << "'!");
      }
    }
  }
}

void LinFlowMechCoupling::DefineDampingIntegrators(const std::string& name,
    shared_ptr<BaseFeFunction>& dispFct,
    shared_ptr<BaseFeFunction>& lmFct,
    shared_ptr<SurfElemList>& actSDList,
    const std::map< RegionIdType, PtrCoefFct >& oneFuncs,
    const std::set< RegionIdType >& flowRegions){
  BiLinearForm * dampInt = NULL;

  if( subType_ == "axi" ) {
    dampInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
        new IdentityOperator<FeH1,2,2>(),
        oneFuncs, 1.0, flowRegions);
  } else if( subType_ == "planeStrain" ) {
    dampInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
        new IdentityOperator<FeH1,2,2>(),
        oneFuncs, 1.0, flowRegions);
  } else if( subType_ == "planeStress" ) {
    dampInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
        new IdentityOperator<FeH1,2,2>(),
        oneFuncs, 1.0, flowRegions);
  } else if( subType_ == "3d") {
    dampInt = new SurfaceABInt<>(new IdentityOperator<FeH1,3,3>(),
        new IdentityOperator<FeH1,3,3>(),
        oneFuncs, 1.0, flowRegions);
  } else {
    EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
  }

  dampInt->SetName(name);
  BiLinFormContext * context =
      new BiLinFormContext(dampInt, DAMPING );

  context->SetEntities( actSDList, actSDList );
  context->SetFeFunctions( lmFct, dispFct );
  context->SetCounterPart(false);

  assemble_->AddBiLinearForm( context );
}

void LinFlowMechCoupling::DefineStiffnessIntegrators(const std::string& name,
    shared_ptr<BaseFeFunction>& dispFct,
    shared_ptr<BaseFeFunction>& velFct,
    shared_ptr<BaseFeFunction>& lmFct,
    shared_ptr<SurfElemList>& actSDList,
    const std::map< RegionIdType, PtrCoefFct >& densityFuncs,
    const std::map< RegionIdType, PtrCoefFct >& muFuncs,
    const std::map< RegionIdType, PtrCoefFct >& oneFuncs,
    const std::set< RegionIdType >& flowRegions){
  BiLinearForm * stiffInt = NULL;

  // LM-velocity integrator in row of LM and column of lin. flow velocity
  std::string intName  = name + "LMVelCouplingInt";
  if( subType_ == "axi" ) {
    stiffInt = new SurfaceABInt<>( new IdentityOperator<FeH1,2,2>(),
        new IdentityOperator<FeH1,2,2>(),
        oneFuncs, -1.0, flowRegions);
  } else if( subType_ == "planeStrain" ) {
    stiffInt = new SurfaceABInt<>(new  IdentityOperator<FeH1,2,2>(),
        new IdentityOperator<FeH1,2,2>(),
        oneFuncs, -1.0, flowRegions);
  } else if( subType_ == "planeStress" ) {
    stiffInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
        new IdentityOperator<FeH1,2,2>(),
        oneFuncs, -1.0, flowRegions);
  } else if( subType_ == "3d") {
    stiffInt = new SurfaceABInt<>( new IdentityOperator<FeH1,3,3>(),
        new IdentityOperator<FeH1,3,3>(),
        oneFuncs, -1.0, flowRegions);
  } else {
    EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
  }

  stiffInt->SetName(intName);
  BiLinFormContext * context =
      new BiLinFormContext( stiffInt, STIFFNESS );

  context->SetEntities( actSDList, actSDList );
  context->SetFeFunctions( lmFct, velFct );
  context->SetCounterPart(false);

  assemble_->AddBiLinearForm( context );

  // Displacement-LM integrator in row of displacement and column of LM
  intName  = name + "DispLMCouplingInt";
  if( subType_ == "axi" ) {
    stiffInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
        new IdentityOperator<FeH1,2,2>(),
        oneFuncs, 1.0, flowRegions);
  } else if( subType_ == "planeStrain" ) {
    stiffInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
        new IdentityOperator<FeH1,2,2>(),
        oneFuncs, 1.0, flowRegions);
  } else if( subType_ == "planeStress" ) {
    stiffInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
        new IdentityOperator<FeH1,2,2>(),
        oneFuncs, 1.0, flowRegions);
  } else if( subType_ == "3d") {
    stiffInt = new SurfaceABInt<>(new IdentityOperator<FeH1,3,3>(),
        new IdentityOperator<FeH1,3,3>(),
        oneFuncs, 1.0, flowRegions);
  } else {
    EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
  }

  stiffInt->SetName(intName);
  context = new BiLinFormContext( stiffInt, STIFFNESS );

  context->SetEntities( actSDList, actSDList );
  context->SetFeFunctions( dispFct, lmFct );
  context->SetCounterPart(false);

  assemble_->AddBiLinearForm( context );

  // Velocity-LM integrator in row of flow velocity and column of LM
  intName  = name + "VelLMCouplingInt";
  if( subType_ == "axi" ) {
    stiffInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
        new IdentityOperator<FeH1,2,2>(),
        oneFuncs, -1.0, flowRegions);
    //          (densityFuncs, 1.0, flowRegions);
  } else if( subType_ == "planeStrain" ) {
    stiffInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
        new IdentityOperator<FeH1,2,2>(),
        oneFuncs, -1.0, flowRegions);
    //        (densityFuncs, 1.0, flowRegions);
  } else if( subType_ == "planeStress" ) {
    stiffInt = new SurfaceABInt<>(new IdentityOperator<FeH1,2,2>(),
        new IdentityOperator<FeH1,2,2>(),
        oneFuncs, -1.0, flowRegions );
    //        (densityFuncs, 1.0, flowRegions);
  } else if( subType_ == "3d") {
    stiffInt = new SurfaceABInt<>(new IdentityOperator<FeH1,3,3>(),
        new IdentityOperator<FeH1,3,3>(),
        oneFuncs, -1.0, flowRegions );
    //        (densityFuncs, 1.0, flowRegions);
  } else {
    EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
  }

  stiffInt->SetName(intName);
  context = new BiLinFormContext( stiffInt, STIFFNESS );

  context->SetEntities( actSDList, actSDList );
  context->SetFeFunctions( velFct, lmFct );
  context->SetCounterPart(false);

  assemble_->AddBiLinearForm( context );

}

void LinFlowMechCoupling::DefineAvailResults() {
  REFACTOR
}

void LinFlowMechCoupling::DefinePrimaryResults() {
  // Check for subType
  if (IsLagrangeMultiplierMethod_){// === LAGRANGE MULTIPLIER ===
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
    shared_ptr<ResultInfo> res1( new ResultInfo);
    res1->resultType = LAGRANGE_MULT;

    res1->dofNames = velDofNames;
    res1->unit = "Pa";
    res1->definedOn = ResultInfo::NODE;
    res1->entryType = ResultInfo::VECTOR;
    feFunctions_[LAGRANGE_MULT]->SetResultInfo(res1);
    results_.Push_back( res1 );
    availResults_.insert( res1 );

    res1->SetFeFunction(feFunctions_[LAGRANGE_MULT]);

    DefineFieldResult( feFunctions_[LAGRANGE_MULT], res1 );
  }
}

void LinFlowMechCoupling::CreateFeSpaces( const std::string&  type,
    PtrParamNode infoNode,
    std::map<SolutionType, shared_ptr<FeSpace> >& crSpaces) {

  if (IsLagrangeMultiplierMethod_){// if the coupling is conforming via lagrange multiplier
    //we need a lagrange multiplier
    formulation_ = LAGRANGE_MULT;
    PtrParamNode spaceNode;

    if(lmOrderSameAsVel_) {
      spaceNode = infoNode->Get(SolutionTypeEnum.ToString(FLUIDMECH_VELOCITY));
    }
    else {
      spaceNode = infoNode->Get(SolutionTypeEnum.ToString(FLUIDMECH_PRESSURE));
    }
    crSpaces[formulation_] =
        FeSpace::CreateInstance(myParam_, spaceNode, FeSpace::H1, ptGrid_);

    crSpaces[formulation_]->SetLagrSurfSpace();

    crSpaces[formulation_]->Init(solStrat_);
  }
}
}
