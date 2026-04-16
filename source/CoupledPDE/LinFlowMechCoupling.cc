// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

//#include <boost/smart_ptr/shared_ptr.hpp>
#include "LinFlowMechCoupling.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/PerturbedFlowPDE.hh"
#include "PDE/LinFlowPDE.hh"
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
  lmOrderSameAsVel_(true), IsLagrangeMultiplierMethod_(false), hasMortarIface_(false){

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

  PtrParamNode ncIfListNode = myParam_->Get("ncInterfaceList", ParamNode::PASS);
  if(ncIfListNode) {
    ParamNodeList ncListMortar =
        ncIfListNode->GetListByVal("ncInterface", "formulation", "Mortar");
    hasMortarIface_ = (ncListMortar.GetSize() > 0);
  }

  // LinFlow FE basis functions
  presPolyId_ = pde1_->GetParamNode()->Get("presPolyId")->As<std::string>();
  velPolyId_  = pde1_->GetParamNode()->Get("velPolyId")->As<std::string>();
  // LinFlow integration order
  presIntegId_ = pde1_->GetParamNode()->Get("presIntegId")->As<std::string>();
  velIntegId_  = pde1_->GetParamNode()->Get("velIntegId")->As<std::string>();

  if( lmOrderSameAsVel_ ) {
	lagrangeMultPolyId_ = velPolyId_;
	lagrangeMultIntegId_ = velIntegId_;
  }
  else {
	lagrangeMultPolyId_ = presPolyId_;
	lagrangeMultIntegId_ = presIntegId_;
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

  // Get LinFlowPDE
  LinFlowPDE* linflowPDE = dynamic_cast<LinFlowPDE*>(pde1_);
  // Get balance of momentum sign of LinFlowPDE
  double linFlowBalanceOfMomentumSign = linflowPDE->GetBalanceOfMomentumSign();

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

    oneFuncs[volRegId] = CoefFunction::Generate(mp, Global::REAL, "1.0");
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

      velSpace->SetRegionApproximation(region, velPolyId_, velIntegId_);
      presSpace->SetRegionApproximation(region, presPolyId_, presIntegId_);
      dispSpace->SetRegionApproximation(region, "default", "default");
      lagrangeMultSpace->SetRegionApproximation(region, lagrangeMultPolyId_, lagrangeMultIntegId_);

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
  else{ //Non conforming integrations for Nitsche or Mortar coupling
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
        IncludeSymmetrizationTerm_ = ncList[i]->Get("IncludeSymmetrizationTerm")->As<bool>();
        MortarInterface * nitscheIf = dynamic_cast<MortarInterface*>(ncIf.get());
        assert(nitscheIf);
        //in case of Nitsche coupling edge/face information is required
        this->ptGrid_->MapEdges();
        this->ptGrid_->MapFaces();
        RegionIdType volMasterId = nitscheIf->GetPrimaryVolRegion();
        RegionIdType volSlaveId = nitscheIf->GetSecondaryVolRegion();

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
        //get correct scaling of penalty term
        Double matVal = 0.0;
        if ( isMaterialComplex_ ) {
          coefMech = mechMat[volMasterId]->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, tensorType, Global::COMPLEX);
          StdVector< Matrix<Complex> > mat;
          coefMech->GetTensorValuesAtCoords(points, mat, this->ptGrid_);
          for (UInt i = 0, numRows = mat[0].GetNumRows(); i < numRows; ++i) {
            matVal += mat[0][i][i].real();
          }

        matVal /= (Double) mat[0].GetNumRows();
        }
        else {
          coefMech = mechMat[volMasterId]->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, tensorType, Global::REAL);
          StdVector< Matrix<Double> > mat;
          coefMech->GetTensorValuesAtCoords(points, mat, this->ptGrid_);
          for (UInt i = 0, numRows = mat[0].GetNumRows(); i < numRows; ++i) {
            matVal += mat[0][i][i];
          }

        matVal /= (Double) mat[0].GetNumRows();
        }
        //estimate of solid elasticity for calculating penalty factor
        std::string strmatVal = boost::lexical_cast<std::string>(matVal);
        PtrCoefFct young = CoefFunction::Generate(mp, Global::REAL, strmatVal);
        //get shear viscosity for calculating penalty factor
        PtrCoefFct shearViscosity = flowMaterials[volSlaveId]->GetScalCoefFnc(FLUID_DYNAMIC_VISCOSITY, Global::REAL);
        shared_ptr<CoefFunction > omega;
        if(analysisType_ == BasePDE::HARMONIC){
          omega = CoefFunction::Generate(mp, Global::REAL, "2*pi*f"); // Angular velocity
        } else {
          EXCEPTION("Nitsche method is only implemented for harmonic simulation"); 
        }
        RegionIdType surfMasterId = nitscheIf->GetPrimarySurfRegion();
        RegionIdType surfSlaveId = nitscheIf->GetSecondarySurfRegion();

        shared_ptr<SurfElemList> elMaster(new SurfElemList(ptGrid_)),
                                 elSlave(new SurfElemList(ptGrid_));
        elMaster->SetRegion(surfMasterId);
        elSlave->SetRegion(surfSlaveId);

        WARN("Developer Info: For LinFlowMech coupling it is recommended to use hdf5 input" );
        // GetElem gives an error (in some cases e.g, PlaneWave2D) while using CDB input
        // Using volMasterId would solve this problem, but produces the same error using hdf5 input
        const Elem* ptElem1 = elMaster->GetElem(surfMasterId);
        const Elem* ptElem2 = elSlave->GetElem(surfSlaveId);

        shared_ptr<ElemShapeMap> esm1T = elMaster->GetGrid()->GetElemShapeMap(ptElem1,true);
        shared_ptr<ElemShapeMap> esm2T = elSlave->GetGrid()->GetElemShapeMap(ptElem2,true);

        //Get minimum element size from solid
        Double min1,max1; 
        esm1T->GetMaxMinEdgeLength(max1, min1);
        std::string str1 = boost::lexical_cast<std::string>(min1);
        PtrCoefFct h_s = CoefFunction::Generate( mp, Global::REAL, str1);
        //Get minimum element size from LinFlow 
        Double min2,max2;
        esm2T->GetMaxMinEdgeLength(max2, min2);
        std::string str2 = boost::lexical_cast<std::string>(min2);
        PtrCoefFct h_f = CoefFunction::Generate(mp, Global::REAL, str2);

        //  betaFluid = shearViscosity/h_f
        PtrCoefFct betaFluid = CoefFunction::Generate(mp,Global::REAL,
        CoefXprBinOp(mp,shearViscosity,h_f, CoefXpr::OP_DIV));

        //  betaSolid= YoungModulus/h_s
        PtrCoefFct betaSolid = CoefFunction::Generate(mp,Global::REAL,
        CoefXprBinOp(mp,young,h_s, CoefXpr::OP_DIV));

        //  betaSolidFreq = YoungModulus/(h_s*omega)
        PtrCoefFct betaSolidFreq = CoefFunction::Generate(mp,Global::REAL,
        CoefXprBinOp(mp,young,CoefXprBinOp(mp,
                h_s,
                omega,
                CoefXpr::OP_MULT
                ), CoefXpr::OP_DIV));
        //create scaling of penalty term
        shared_ptr<CoefFunction> beta_scalingCoef;
        if (IncludeSymmetrizationTerm_ == true){
          beta_scalingCoef = CoefFunction::Generate(mp,Global::REAL,
          CoefXprBinOp(mp,betaFluid,CoefXprBinOp(mp,
                betaSolid,
                betaSolidFreq,
                CoefXpr::OP_ADD
                ),CoefXpr::OP_ADD));
        }else{
          beta_scalingCoef = CoefFunction::Generate(mp,Global::REAL,
          CoefXprBinOp(mp,betaFluid,betaSolidFreq,CoefXpr::OP_ADD));
        }

        PtrCoefFct coefOne = CoefFunction::Generate(mp, Global::REAL, "1");
        PtrCoefFct coefMinusOne = CoefFunction::Generate(mp, Global::REAL, "-1.0");
        
        // Note:
        // Define the bilinear forms for the Nitsche coupling consists of six terms
        // The terms are divided into four categories:
        // 1) The bilinear forms coupling mechanical displacement with LinFlow velocity
        // 2) The Penalty integrals coupling mechanical displacement with LinFlow velocity
        // 3) The bilinear forms coupling LinFlow velocity with mechanical displacement
        // 4) The Penalty integrals coupling LinFlow velocity with mechanical displacement
        // assume phi is displacement test function and u is displacement
        // psi is velocity test function and v LinFlow velocity
        // 5) Symmetrization terms are seperately added to formulation

        BiLinearForm* int_PhiM_duM = NULL;   // Term 1: u'M_duM : phi * sigma[u].n
        BiLinearForm* penalty_PhiM_uM = NULL;// term 2  u'M_uM : beta*phi*du/dt
        BiLinearForm* penalty_PhiM_vS = NULL;// term 3 u'M_vS :-beta*phi*v
        BiLinearForm* int_PsiS_duM = NULL;   // Term 4: v'S_duM : Psi*sigma[u].n
        BiLinearForm* penalty_PsiS_vS = NULL;// Term 5: v'S_vS : beta*Psi.v
        BiLinearForm* penalty_PsiS_uM = NULL;// Term 6: v'S_uM : -beta*Psi.du/dt

        BaseBOperator * sNSOp1 = NULL;
        BaseBOperator * sNSOp2 = NULL;
        BaseBOperator * sNSOp3 = NULL;
        BaseBOperator * sNSOp4 = NULL;
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
              sNSOp1, coefMinusOne, 1.0, BiLinearForm::PRIM_PRIM, false, true, false);
        }
        else if (dim_ == 2 && subType_ == "2.5d")
        {
          sNSOp1 = new SurfaceNormalStressOperator<FeH1, 2, 3,Complex>(subType_ ,false);
          sNSOp1->SetCoefFunction(coefMech);
          int_PhiM_duM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 3>(),
              sNSOp1, coefMinusOne, 1.0, BiLinearForm::PRIM_PRIM, false, true, false);
        }
        else
        {
          sNSOp1 = new SurfaceNormalStressOperator<FeH1, 2, 2,Complex>(subType_ ,false);
          sNSOp1->SetCoefFunction(coefMech);
          int_PhiM_duM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 2>(),
              sNSOp1, coefMinusOne, 1.0, BiLinearForm::PRIM_PRIM, false, true, false);
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
              beta_scalingCoef, beta, BiLinearForm::PRIM_PRIM , false, true, true);
        }
        else if (dim_ == 2 && subType_ == "2.5d")
        {
          penalty_PhiM_uM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 3>(),
              new SurfaceIdentityOperator<FeH1, 2, 3>(),
              beta_scalingCoef, beta, BiLinearForm::PRIM_PRIM, false, true, true);
        }
        else
        {

          penalty_PhiM_uM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 2>(),
              new SurfaceIdentityOperator<FeH1, 2, 2>(),
              beta_scalingCoef, beta, BiLinearForm::PRIM_PRIM, false, true, true);
        }
        // --------------------------
        // Term 3 : -beta*phi*v
        // --------------------------
        if (dim_ == 3)
        {
          penalty_PhiM_vS = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 3, 3>(),
              new SurfaceIdentityOperator<FeH1, 3, 3>(),
              beta_scalingCoef, -beta,BiLinearForm::PRIM_SEC , false, true, true);
        }
        else if (dim_ == 2 && subType_ == "2.5d")
        {
          penalty_PhiM_vS = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 3>(),
              new SurfaceIdentityOperator<FeH1, 2, 3>(),
              beta_scalingCoef, -beta, BiLinearForm::PRIM_SEC , false, true, true);
        }
        else
        {
          penalty_PhiM_vS = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 2>(),
              new SurfaceIdentityOperator<FeH1, 2, 2>(),
              beta_scalingCoef, -beta, BiLinearForm::PRIM_SEC , false, true, true);
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
              coefOne, linFlowBalanceOfMomentumSign, BiLinearForm::SEC_PRIM, false, true, false);
        }
        else if (dim_ == 2 && subType_ == "2.5d")
        {
          sNSOp2 = new SurfaceNormalStressOperator<FeH1, 2, 3,Complex>(subType_ ,false);
          sNSOp2->SetCoefFunction(coefMech);
          int_PsiS_duM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 3>(),sNSOp2,
              coefOne, linFlowBalanceOfMomentumSign, BiLinearForm::SEC_PRIM, false, true, false);
        }
        else
        {
          sNSOp2 = new SurfaceNormalStressOperator<FeH1, 2, 2,Complex>(subType_ ,false);
          sNSOp2->SetCoefFunction(coefMech);
          int_PsiS_duM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 2>(),sNSOp2,
              coefOne, linFlowBalanceOfMomentumSign, BiLinearForm::SEC_PRIM, false, false, false);
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
              beta_scalingCoef, linFlowBalanceOfMomentumSign*beta, BiLinearForm::SEC_SEC, false, true, true);
        }
        else if (dim_ == 2 && subType_ == "2.5d")
        {
          penalty_PsiS_vS = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 3>(),
              new SurfaceIdentityOperator<FeH1, 2, 3>(),
              beta_scalingCoef, linFlowBalanceOfMomentumSign*beta, BiLinearForm::SEC_SEC, false, true, true);
        }
        else
        {
          penalty_PsiS_vS = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 2>(),
              new SurfaceIdentityOperator<FeH1, 2, 2>(),
              beta_scalingCoef, linFlowBalanceOfMomentumSign*beta, BiLinearForm::SEC_SEC, false, false, true);
        }
        // --------------------------
        // Term 6 : -beta*Psi.du/dt
        // --------------------------

        if (dim_ == 3)
        {
          penalty_PsiS_uM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 3, 3>(),
              new SurfaceIdentityOperator<FeH1, 3, 3>(),
              beta_scalingCoef, -linFlowBalanceOfMomentumSign*beta, BiLinearForm::SEC_PRIM, false, true, true);
        }
        else if (dim_ == 2 && subType_ == "2.5d")
        {
          penalty_PsiS_uM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 3>(),
              new SurfaceIdentityOperator<FeH1, 2, 3>(),
              beta_scalingCoef, -linFlowBalanceOfMomentumSign*beta, BiLinearForm::SEC_PRIM, false, true, true);
        }
        else
        {
          penalty_PsiS_uM = new SurfaceNitscheABInt<Double,Double>
          (new SurfaceIdentityOperator<FeH1, 2, 2>(),
              new SurfaceIdentityOperator<FeH1, 2, 2>(),
              beta_scalingCoef, -linFlowBalanceOfMomentumSign*beta, BiLinearForm::SEC_PRIM, false, false, true);
        }

        int_PhiM_duM->SetName("int_PhiM_duM");
        penalty_PhiM_uM->SetName("penalty_PhiM_uM");
        penalty_PhiM_vS->SetName("penalty_PhiM_vS");
        int_PsiS_duM->SetName("int_PsiS_duM");
        penalty_PsiS_vS->SetName("penalty_PsiS_vS");
        penalty_PsiS_uM->SetName("penalty_PsiS_uM");


        // define contexts for bilinear forms
        SurfaceBiLinFormContext* descr_PhiM_duM = new SurfaceBiLinFormContext(int_PhiM_duM, STIFFNESS, BiLinearForm::PRIM_PRIM);
        SurfaceBiLinFormContext* descr_PhiM_uM = new SurfaceBiLinFormContext(penalty_PhiM_uM, DAMPING , BiLinearForm::PRIM_PRIM);
        SurfaceBiLinFormContext* descr_PhiM_vS = new SurfaceBiLinFormContext(penalty_PhiM_vS, STIFFNESS, BiLinearForm::PRIM_SEC);
        SurfaceBiLinFormContext* descr_PsiS_duM = new SurfaceBiLinFormContext(int_PsiS_duM, STIFFNESS, BiLinearForm::SEC_PRIM);
        SurfaceBiLinFormContext* descr_PsiS_vS = new SurfaceBiLinFormContext(penalty_PsiS_vS, STIFFNESS, BiLinearForm::SEC_SEC);
        SurfaceBiLinFormContext* descr_PsiS_uM = new SurfaceBiLinFormContext(penalty_PsiS_uM, DAMPING, BiLinearForm::SEC_PRIM);

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

        if (IncludeSymmetrizationTerm_ == true){//Symmetrization terms is added
          BiLinearForm* int_dPhiM_vS = NULL;    // Term 7: du'M_vS : sigma[u'].n*v
          BiLinearForm* int_dPhiM_uM = NULL;    // Term 8: du'M_Um : -sigma[u'].n*du/dt

          // --------------------------
          // Term 7 : du'M_vS : sigma[u'].n*v
          // --------------------------
          if (dim_ == 3)
          {
            sNSOp3 = new SurfaceNormalStressOperator<FeH1, 3, 3,Complex>(subType_ ,false);
            sNSOp3->SetCoefFunction(coefMech);
            int_dPhiM_vS = new SurfaceNitscheABInt<Double,Double>
            (sNSOp3, new SurfaceIdentityOperator<FeH1, 3, 3>(),
                coefOne, 1.0, BiLinearForm::PRIM_SEC, false, true, false);
          }
          else if (dim_ == 2 && subType_ == "2.5d")
          {
            sNSOp3 = new SurfaceNormalStressOperator<FeH1, 2, 3,Complex>(subType_ ,false);
            sNSOp3->SetCoefFunction(coefMech);
            int_dPhiM_vS = new SurfaceNitscheABInt<Double,Double>
            (sNSOp3, new SurfaceIdentityOperator<FeH1, 2, 3>(),
                coefOne, 1.0, BiLinearForm::PRIM_SEC, false, true, false);
          }
          else
          {
            sNSOp3 = new SurfaceNormalStressOperator<FeH1, 2, 2,Complex>(subType_ ,false);
            sNSOp3->SetCoefFunction(coefMech);
            int_dPhiM_vS = new SurfaceNitscheABInt<Double,Double>
            (sNSOp3, new SurfaceIdentityOperator<FeH1, 2, 2>(),
                coefOne, 1.0, BiLinearForm::PRIM_SEC, false, true, false);
          }

          // --------------------------
          // Term 8 : du'M_Um : -sigma[u'].n*du/dt
          // --------------------------
          if (dim_ == 3)
          {
            sNSOp4 = new SurfaceNormalStressOperator<FeH1, 3, 3,Complex>(subType_ ,false);
            sNSOp4->SetCoefFunction(coefMech);
            int_dPhiM_uM = new SurfaceNitscheABInt<Double,Double>
            (sNSOp4, new SurfaceIdentityOperator<FeH1, 3, 3>(),
                coefMinusOne, 1.0, BiLinearForm::PRIM_PRIM, false, true, false);
          }
          else if (dim_ == 2 && subType_ == "2.5d")
          {
            sNSOp4 = new SurfaceNormalStressOperator<FeH1, 2, 3,Complex>(subType_ ,false);
            sNSOp4->SetCoefFunction(coefMech);
            int_dPhiM_uM = new SurfaceNitscheABInt<Double,Double>
            (sNSOp4, new SurfaceIdentityOperator<FeH1, 2, 3>(),
                coefMinusOne, 1.0, BiLinearForm::PRIM_PRIM, false, true, false);
          }
          else
          {
            sNSOp4 = new SurfaceNormalStressOperator<FeH1, 2, 2,Complex>(subType_ ,false);
            sNSOp4->SetCoefFunction(coefMech);
            int_dPhiM_uM = new SurfaceNitscheABInt<Double,Double>
            (sNSOp4, new SurfaceIdentityOperator<FeH1, 2, 2>(),
                coefMinusOne, 1.0, BiLinearForm::PRIM_PRIM, false, true, false);
          }


          int_dPhiM_vS->SetName("int_dPhiM_vS");
          int_dPhiM_uM->SetName("int_dPhiM_uM");

          SurfaceBiLinFormContext* descr_dPhiM_vS = new SurfaceBiLinFormContext(int_dPhiM_vS, STIFFNESS, BiLinearForm::PRIM_SEC);
          SurfaceBiLinFormContext* descr_dPhiM_uM = new SurfaceBiLinFormContext(int_dPhiM_uM, DAMPING, BiLinearForm::PRIM_PRIM);

          descr_dPhiM_vS->SetEntities(actNCSDList, actNCSDList);
          descr_dPhiM_uM->SetEntities(actNCSDList, actNCSDList);

          descr_dPhiM_vS->SetFeFunctions(dispFct, velFct);
          descr_dPhiM_uM->SetFeFunctions(dispFct, velFct);

          assemble_->AddBiLinearForm(descr_dPhiM_vS);
          assemble_->AddBiLinearForm(descr_dPhiM_uM);
        }
      }
      else if (formulation == "Mortar")
      {
        shared_ptr<BaseFeFunction> presFct = pde1_->GetFeFunction(FLUIDMECH_PRESSURE);
        shared_ptr<BaseFeFunction> lagrangeMultFct = feFunctions_[LAGRANGE_MULT];
        shared_ptr<FeSpace> presSpace = presFct->GetFeSpace();
        shared_ptr<FeSpace> lagrangeMultSpace = lagrangeMultFct->GetFeSpace();

        MortarInterface * mortarIf = dynamic_cast<MortarInterface*>(ncIf.get());
        assert(mortarIf);

        RegionIdType surfMasterId = mortarIf->GetPrimarySurfRegion();
        RegionIdType surfSlaveId = mortarIf->GetSecondarySurfRegion();

        // create ElemLists for master & slave surfaces
        shared_ptr<SurfElemList> elMaster(new SurfElemList(ptGrid_)),
                                 elSlave(new SurfElemList(ptGrid_));
        elMaster->SetRegion(surfMasterId);
        elSlave->SetRegion(surfSlaveId);

        velFct->AddEntityList(elSlave);
        presFct->AddEntityList(elSlave);
        dispFct->AddEntityList(elMaster);
        lagrangeMultFct->AddEntityList(elSlave);

        velSpace->SetRegionApproximation(surfSlaveId, velPolyId_, velIntegId_);
        presSpace->SetRegionApproximation(surfSlaveId, presPolyId_, presIntegId_);
        dispSpace->SetRegionApproximation(surfMasterId, "default", "default");
        lagrangeMultSpace->SetRegionApproximation(surfSlaveId, lagrangeMultPolyId_, lagrangeMultIntegId_);


        PtrCoefFct coefOne = CoefFunction::Generate( mp, Global::REAL, "1");

        // Coupling across NC interface -> SurfaceMortarABInt
        DefineMortarIntNC("LinFlowMechDampingLMDispCouplingIntNC", lagrangeMultFct, dispFct,
                          mortarIf, 1.0, coefOne, DAMPING, BiLinearForm::SEC_PRIM);
        DefineMortarIntNC("LinFlowMechStiffDispLMCouplingIntNC", dispFct, lagrangeMultFct,
                          mortarIf, 1.0, coefOne, STIFFNESS, BiLinearForm::PRIM_SEC);
        // Coupling on secondary (slave) side -> BBInt
        DefineMortarIntNCSecondary("LinFlowMechStiffLMVelCouplingIntNC", lagrangeMultFct, velFct,
                          elSlave, -1.0, coefOne, STIFFNESS);
        DefineMortarIntNCSecondary("LinFlowMechStiffVelLMCouplingIntNC", velFct, lagrangeMultFct,
                          elSlave, -linFlowBalanceOfMomentumSign, coefOne, STIFFNESS);
      }
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

void LinFlowMechCoupling::DefineMortarIntNC(const std::string& name,
    shared_ptr<BaseFeFunction>& fct1,
    shared_ptr<BaseFeFunction>& fct2,
    MortarInterface* mortarIf,
    Double factor,
    PtrCoefFct scalCoef,
    FEMatrixType matType,
    BiLinearForm::CouplingDirection cplDir){

  shared_ptr<ElemList> actSDList = mortarIf->GetElemList();

  RegionIdType volMasterId = mortarIf->GetPrimaryVolRegion();
  RegionIdType volSlaveId = mortarIf->GetSecondaryVolRegion();
  bool coplanar = mortarIf->IsCoplanar();

  BiLinearForm * cplInt = NULL;

  if( subType_ == "axi" || subType_ == "planeStrain" || subType_ == "planeStress" ) {
    cplInt = new SurfaceMortarABInt<>(new IdentityOperator<FeH1,2,2>(),
        new IdentityOperator<FeH1,2,2>(),
        scalCoef, factor, volMasterId, volSlaveId, coplanar, false, cplDir);
  } else if( subType_ == "3d") {
    cplInt = new SurfaceMortarABInt<>(new IdentityOperator<FeH1,3,3>(),
        new IdentityOperator<FeH1,3,3>(),
        scalCoef, factor, volMasterId, volSlaveId, coplanar, false, cplDir);
  } else {
    EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
  }

  cplInt->SetName(name);
  NcBiLinFormContext * context =
      new NcBiLinFormContext( cplInt, matType );

  context->SetEntities( actSDList, actSDList );
  context->SetFeFunctions( fct1, fct2 );
  context->SetCounterPart(false);

  assemble_->AddBiLinearForm( context );
  mortarIf->RegisterIntegrator( context );
}

void LinFlowMechCoupling::DefineMortarIntNCSecondary(const std::string& name,
    shared_ptr<BaseFeFunction>& fct1,
    shared_ptr<BaseFeFunction>& fct2,
    shared_ptr<SurfElemList>& actSDList,
    Double factor,
    PtrCoefFct scalCoef,
    FEMatrixType matType){

  BiLinearForm * cplInt = NULL;

  if( subType_ == "axi" || subType_ == "planeStrain" || subType_ == "planeStress" ) {
    cplInt = new BBInt<>(new IdentityOperator<FeH1,2,2>(),
        scalCoef, factor, false);
  } else if( subType_ == "3d") {
    cplInt = new BBInt<>(new IdentityOperator<FeH1,3,3>(),
        scalCoef, factor, false);
  } else {
    EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
  }

  cplInt->SetName(name);
  BiLinFormContext * context =
      new BiLinFormContext( cplInt, matType );

  context->SetEntities( actSDList, actSDList );
  context->SetFeFunctions( fct1, fct2 );
  context->SetCounterPart(false);

  assemble_->AddBiLinearForm( context );
}

void LinFlowMechCoupling::DefineAvailResults() {
  REFACTOR
}

void LinFlowMechCoupling::DefinePrimaryResults() {
  // Check for subType
  if(IsLagrangeMultiplierMethod_ || hasMortarIface_) {// === LAGRANGE MULTIPLIER ===
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

  // coupling via Lagrange multiplier (conforming or non-conforming Mortar)
  if (IsLagrangeMultiplierMethod_ || hasMortarIface_) {
    formulation_ = LAGRANGE_MULT;
    PtrParamNode spaceNode;
    // for the Mortar method, the LM space must be on the "secondary" side (LinFlow)
    PtrParamNode ParamNodeLM = IsLagrangeMultiplierMethod_ ? myParam_ : pde1_->GetParamNode();

    if(lmOrderSameAsVel_) {
      spaceNode = infoNode->Get(SolutionTypeEnum.ToString(FLUIDMECH_VELOCITY));
    }
    else {
      spaceNode = infoNode->Get(SolutionTypeEnum.ToString(FLUIDMECH_PRESSURE));
    }
    
    crSpaces[formulation_] =
        FeSpace::CreateInstance(ParamNodeLM, spaceNode, FeSpace::H1, ptGrid_);

    crSpaces[formulation_]->SetLagrSurfSpace();

    crSpaces[formulation_]->Init(solStrat_);
  }
}

}
