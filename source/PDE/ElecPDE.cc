#include <fstream>
#include <iostream>
#include <sstream>
#include <math.h>
#include <string>
#include <set>

#include "ElecPDE.hh"

#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefFunctionPML.hh"
#include "Domain/Mesh/NcInterfaces/MortarInterface.hh"
#include "Domain/CoefFunction/CoefFunctionHyst.hh"
#include "Utils/StdVector.hh"
#include "Driver/SolveSteps/SolveStepElec.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Driver/Assemble.hh"
#include "Utils/ApproxData.hh"
#include "Utils/SmoothSpline.hh"
#include "Materials/Models/Hysteresis.hh"
#include "Materials/Models/Preisach.hh"
#include "Materials/Models/VectorPreisach.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"
#include "FeBasis/FeFunctions.hh"

// new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/DivOperator.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdentityOperatorNormal.hh"
#include "Forms/Operators/SurfaceOperators.hh"
#include "Forms/Operators/SurfaceNormalFluxDensityOperator.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"


namespace CoupledField {

  DECLARE_LOG(elecpde)
  DEFINE_LOG(elecpde, "pde.electrostatic")

  // ***************
  //   Constructor
  // ***************
  ElecPDE::ElecPDE( Grid* aptgrid, PtrParamNode paramNode,
                    PtrParamNode infoNode,
                    shared_ptr<SimState> simState, Domain* domain )
    :SinglePDE( aptgrid, paramNode, infoNode, simState, domain ) {


    // =====================================================================
    // set solution information
    // =====================================================================
    pdename_          = "electrostatic";
    pdematerialclass_ = ELECTROSTATIC;
 
    nonLin_         = false;
    nonLinMaterial_ = false;
//    isAlwaysStatic_ = true;
    isPiezoCoupled_ = false;
    
    //! Always use updated Lagrangian formulation 
    updatedGeo_     = true;
    
    
    // Check the subtype of the problem
    paramNode->GetValue("subType", subType_);
  }
  
  void ElecPDE::InitNonLin() {

    SinglePDE::InitNonLin();
  }

  void ElecPDE::ReadDampingInformation()
  {

  }

  void ElecPDE::DefineIntegrators() {

    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;
    
    // bool upLagrangeForm = true;

    //transform the type
    SubTensorType tensorType;

    if ( dim_ == 3  || subType_ == "2.5d" ) {
      tensorType = FULL;
    }
    else {
      if ( isaxi_ == true ) {
        tensorType = AXI;
      }
      else {
        // 2d: plane case
        tensorType = PLANE_STRAIN;
      }
    }


    // if the pde is piezo-coupled, the electrostatic entries
    // have to be multiplied with -1
    std::string factor = "1.0";
    if ( isPiezoCoupled_ == true )
      factor = "-1.0";  

    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    shared_ptr<FeSpace> mySpace = feFunctions_[ELEC_POTENTIAL]->GetFeSpace();

    //flag indicating frequency PML formulation
    bool harmonicPML = false;

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

      // Take account for pml in frequency domain
      // 'coeffPMLScal' is the function used to scale the material tensor. If no PML is defined, it's unity
      PtrCoefFct coefPMLScal;
      PtrCoefFct coefPMLVec, speedOfSnd;
      PtrParamNode pmlNode;
      std::string pmlFormul;
      //
      
      if (dampingList_[actRegion] == PML)
      {
        if (analysistype_ == HARMONIC)
        {
          harmonicPML = true;
          std::string dampId;
          curRegNode->GetValue("dampingId", dampId);
          pmlNode = myParam_->Get("dampingList")->GetByVal("pml", "id", dampId.c_str());
          pmlFormul = pmlNode->Get("formulation")->As<std::string>();

          // speed of sound is set to equal '1.0'
          speedOfSnd = CoefFunction::Generate(mp_, Global::REAL, "1.0");
          if (pmlFormul == "classic")
          {
            coefPMLScal.reset(new CoefFunctionPML<Complex>(pmlNode, speedOfSnd,
                              ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, false));
            coefPMLVec.reset(new CoefFunctionPML<Complex>(pmlNode, speedOfSnd,
                             ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, true));
          }
          else if (pmlFormul == "shifted")
          {
            coefPMLScal.reset(new CoefFunctionShiftedPML<Complex>(pmlNode, speedOfSnd,
                              ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, false));
            coefPMLVec.reset(new CoefFunctionShiftedPML<Complex>(pmlNode, speedOfSnd,
                             ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, true));
          }
          else
          {
            EXCEPTION("Unknown PML-formulation '" << pmlFormul << "'");
          }
        }
        else
          EXCEPTION("Not implemented yet");
      }
      else
      {
        harmonicPML = false;
      }
      
      // ----- standard real-valued stiffness integrator
      BaseBDBInt* stiffInt = NULL;
      if (harmonicPML)
      {
        stiffInt = GetStiffIntegrator(actSDMat, tensorType, actRegion, coefPMLScal);
        stiffInt->GetBOp()->SetCoefFunction(coefPMLVec);
      }
      else
        stiffInt = GetStiffIntegrator(actSDMat, tensorType, actRegion);

      stiffInt->SetName("LinElecIntegrator");
      BiLinFormContext * stiffIntDescr =
          new BiLinFormContext(stiffInt, STIFFNESS );
      feFunctions_[ELEC_POTENTIAL]->AddEntityList( actSDList );

      //stiffIntDescr->SetPtPdes(this, this);
      stiffIntDescr->SetEntities( actSDList, actSDList );
      stiffIntDescr->SetFeFunctions( feFunctions_[ELEC_POTENTIAL],feFunctions_[ELEC_POTENTIAL]);
      stiffInt->SetFeSpace( feFunctions_[ELEC_POTENTIAL]->GetFeSpace());

      assemble_->AddBiLinearForm( stiffIntDescr );
      // Important: Add bdb-integrator to global list, as we need them later
      // for calculation of postprocessing results
      bdbInts_[actRegion] = stiffInt;
      
    }
    
    // Define integrators for composite materials
    // (only for subType "flatShell")
    std::map<RegionIdType, Composite>::iterator compIt;
    for( compIt=compositeMaterials_.begin(); compIt!=compositeMaterials_.end();
         compIt++ ) {
      
      REFACTOR;
    }

    // define integrators for electric impedances
    DefineImpedanceIntegrators();
  }
  
  void ElecPDE::DefineSurfaceIntegrators()
  {
    PtrParamNode bcNode = myParam_->Get("bcsAndLoads", ParamNode::PASS);
    this->ptGrid_->MapEdges();
    this->ptGrid_->MapFaces();
    if (!bcNode)
      return;

    ParamNodeList blochNodesList = bcNode->GetList("blochPeriodic");
    for (UInt i = 0; i < blochNodesList.GetSize(); i++)
    {
      std::string str_value = blochNodesList[i]->Get("factor_value")->As<std::string>();
      std::string str_phase = blochNodesList[i]->Get("factor_phase")->As<std::string>();
      std::string formulation = blochNodesList[i]->Get("formulation")->As<std::string>();

      // propagation factor \gamma from xml-file
      std::string str_real, str_imag;
      str_real = AmplPhaseToReal(str_value, str_phase, true);
      str_imag = AmplPhaseToImag(str_value, str_phase, true);

      PtrCoefFct factor = CoefFunction::Generate(mp_, Global::COMPLEX, str_real, str_imag);
      PtrCoefFct one = CoefFunction::Generate(mp_, Global::COMPLEX, "1.0", "0.0");

      ParamNodeList regionsList = blochNodesList[i]->GetList("region");
      for (UInt j = 0; j < regionsList.GetSize(); j++)
      {
        std::string ncRegionName = regionsList[j]->Get("name")->As<std::string>();
        shared_ptr<BaseNcInterface> ncIf = ptGrid_->GetNcInterface(ptGrid_->GetNcInterfaceId(ncRegionName));
        if (!ncIf)
        {
          EXCEPTION("No interface with the name '" << ncRegionName << "' found!");
        }
        shared_ptr<MortarInterface> mortarIf = boost::dynamic_pointer_cast<MortarInterface>(ncIf);
        assert(mortarIf);

        // TODO: WHAT THE HELL? IT DOESN'T WORK WITH pz = -1...
        // To be checked later during the work with piezo-coupling
        Double pz = 1.0;
        if (isPiezoCoupled_ == true)
          pz = 1.0;

        if (formulation == "Nitsche")
        {
          PtrCoefFct matDataTensorMas, matDataTensorSla, matData;
          RegionIdType volMasterId = mortarIf->GetMasterVolRegion();
          RegionIdType volSlaveId = mortarIf->GetSlaveVolRegion();

          matDataTensorMas = regionPermittivity_[volMasterId];
          matDataTensorSla = regionPermittivity_[volSlaveId];
          assert(matDataTensorMas);
          assert(matDataTensorSla);

          std::string nitFac = blochNodesList[i]->Get("nitscheFactor")->As<std::string>();
          Double nitscheFactor = lexical_cast<Double>(nitFac);
          // master & slave penalty integrals
          BiLinearForm *pnlt_PhiM_PsiM = NULL;
          BiLinearForm *pnlt_PhiM_PsiS = NULL;
          BiLinearForm *pnlt_PhiS_PsiM = NULL;
          BiLinearForm *pnlt_PhiS_PsiS = NULL;
          // master & slave integrals with normal derivatives
          BiLinearForm *flux_DPhiM_PsiM = NULL;
          BiLinearForm *flux_DPhiM_PsiS = NULL;
          BiLinearForm *flux_PhiM_DPsiM = NULL;
          BiLinearForm *flux_PhiS_DPsiM = NULL;

          shared_ptr<ElemList> actSDList = ncIf->GetElemList();
          Double beta;
          BiLinearForm::CouplingDirection cplDir;
          PtrCoefFct factorSqr = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, factor, factor, CoefXpr::OP_MULT));

          // obtain a proper scaling for the penalty terms
          StdVector<Vector<Double> > points(1);
          Vector<Double> p1(dim_);
          p1.Init();
          points[0]= p1;

          if (matDataTensorMas->IsComplex())
          {
            Complex tmp(0.0, 0.0);
            StdVector<Matrix<Complex> > valsM, valsS;
            matDataTensorMas->GetTensorValuesAtCoords(points, valsM, this->ptGrid_);
            matDataTensorSla->GetTensorValuesAtCoords(points, valsS, this->ptGrid_);
            for (UInt k = 0; k < valsM[0].GetNumRows(); k++)
            {
              tmp += valsM[0][k][k]*conj(valsM[0][k][k]);
              tmp += valsS[0][k][k]*conj(valsS[0][k][k]);
            }
            beta = sqrt(0.5*tmp.real())*nitscheFactor;
          }
          else
          {
            Double tmp(0.0);
            StdVector<Matrix<Double> > valsM, valsS;
            matDataTensorMas->GetTensorValuesAtCoords(points, valsM, this->ptGrid_);
            matDataTensorSla->GetTensorValuesAtCoords(points, valsS, this->ptGrid_);
            for (UInt k = 0; k < valsM[0].GetNumRows(); k++)
            {
              tmp += valsM[0][k][k]*valsM[0][k][k];
              tmp += valsS[0][k][k]*valsS[0][k][k];
            }
            beta = sqrt(0.5*tmp)*nitscheFactor;
          }

          PtrCoefFct coefFuncPMLVec, coefFuncPMLScl;
          if (dampingList_[volMasterId] == PML)
          {
            if (analysistype_ == HARMONIC)
            {
              std::string regionName = ptGrid_->GetRegion().ToString(volMasterId);
              std::string dampId, pmlFormul;

              PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region", "name", regionName.c_str());
              curRegNode->GetValue("dampingId", dampId);
              PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal("pml", "id", dampId.c_str());

              // speed of sound is set to equal '1.0'
              PtrCoefFct speedOfSnd = CoefFunction::Generate(mp_, Global::REAL, "1.0");
              pmlFormul = pmlNode->Get("formulation")->As<std::string>();

              if (pmlFormul == "classic")
              {
                coefFuncPMLVec.reset(new CoefFunctionPML<Complex>(pmlNode, speedOfSnd,
                                     ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, true));
                coefFuncPMLScl.reset(new CoefFunctionPML<Complex>(pmlNode, speedOfSnd,
                                     ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, false));
              }
              else if (pmlFormul == "shifted")
              {
                coefFuncPMLVec.reset(new CoefFunctionShiftedPML<Complex>(pmlNode, speedOfSnd,
                                     ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, true));
                coefFuncPMLScl.reset(new CoefFunctionShiftedPML<Complex>(pmlNode, speedOfSnd,
                                     ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, false));
              }
              else
              {
                EXCEPTION("Unknown PML-formulation '" << pmlFormul << "'");
              }

              matData = CoefFunction::Generate(mp_, Global::COMPLEX,
                                               CoefXprTensScalOp(mp_, matDataTensorMas, coefFuncPMLScl, CoefXpr::OP_MULT));
            }
            else
              EXCEPTION("The analysis type '" << this->analysistype_ << "' is not supported!");
          }
          else
            matData = matDataTensorMas;

          // define bilinear forms for Nitsche coupling
          cplDir = BiLinearForm::MASTER_MASTER;
          pnlt_PhiM_PsiM = GetPenaltyIntegrator(factor, beta*pz, cplDir);
          flux_DPhiM_PsiM = GetFluxIntegrator(one, coefFuncPMLVec, -1.0*pz, cplDir, true);
          flux_PhiM_DPsiM = GetFluxIntegrator(factor, coefFuncPMLVec, -1.0*pz, cplDir, false);
          flux_DPhiM_PsiM->SetBCoefFunctionOpA(matData);
          flux_PhiM_DPsiM->SetBCoefFunctionOpB(matData);

          cplDir = BiLinearForm::MASTER_SLAVE;
          pnlt_PhiM_PsiS = GetPenaltyIntegrator(factorSqr, -beta*pz, cplDir);
          flux_DPhiM_PsiS = GetFluxIntegrator(factor, coefFuncPMLVec, 1.0*pz, cplDir, true);
          flux_DPhiM_PsiS->SetBCoefFunctionOpA(matData);

          cplDir = BiLinearForm::SLAVE_MASTER;
          pnlt_PhiS_PsiM = GetPenaltyIntegrator(one, -beta*pz, cplDir);
          flux_PhiS_DPsiM = GetFluxIntegrator(one, coefFuncPMLVec, 1.0*pz, cplDir, false);
          flux_PhiS_DPsiM->SetBCoefFunctionOpB(matData);

          cplDir = BiLinearForm::SLAVE_SLAVE;
          pnlt_PhiS_PsiS = GetPenaltyIntegrator(factor, beta, cplDir);

          // master-master
          pnlt_PhiM_PsiM->SetName("pnlt_PhiM_PsiM");
          flux_DPhiM_PsiM->SetName("flux_DPhiM_PsiM");
          flux_PhiM_DPsiM->SetName("flux_PhiM_DPsiM");
          //master-slave
          pnlt_PhiM_PsiS->SetName("pnlt_PhiM_PsiS");
          flux_DPhiM_PsiS->SetName("flux_DPhiM_PsiS");
          // slave-master
          pnlt_PhiS_PsiM->SetName("pnlt_PhiS_PsiM");
          flux_PhiS_DPsiM->SetName("flux_PhiS_DPsiM");
          //slave-slave
          pnlt_PhiS_PsiS->SetName("pnlt_PhiS_PsiS");

          cplDir = BiLinearForm::MASTER_MASTER;
          SurfaceBiLinFormContext *pnlt_PhiM_PsiM_cont = new SurfaceBiLinFormContext(pnlt_PhiM_PsiM, STIFFNESS, cplDir);
          SurfaceBiLinFormContext *flux_DPhiM_PsiM_cont = new SurfaceBiLinFormContext(flux_DPhiM_PsiM, STIFFNESS, cplDir);
          SurfaceBiLinFormContext *flux_PhiM_DPsiM_cont = new SurfaceBiLinFormContext(flux_PhiM_DPsiM, STIFFNESS, cplDir);
          cplDir = BiLinearForm::MASTER_SLAVE;
          SurfaceBiLinFormContext *pnlt_PhiM_PsiS_cont = new SurfaceBiLinFormContext(pnlt_PhiM_PsiS, STIFFNESS, cplDir);
          SurfaceBiLinFormContext *flux_DPhiM_PsiS_cont = new SurfaceBiLinFormContext(flux_DPhiM_PsiS, STIFFNESS, cplDir);
          cplDir = BiLinearForm::SLAVE_MASTER;
          SurfaceBiLinFormContext *pnlt_PhiS_PsiM_cont = new SurfaceBiLinFormContext(pnlt_PhiS_PsiM, STIFFNESS, cplDir);
          SurfaceBiLinFormContext *flux_PhiS_DPsiM_cont = new SurfaceBiLinFormContext(flux_PhiS_DPsiM, STIFFNESS, cplDir);
          cplDir = BiLinearForm::SLAVE_SLAVE;
          SurfaceBiLinFormContext *pnlt_PhiS_PsiS_cont = new SurfaceBiLinFormContext(pnlt_PhiS_PsiS, STIFFNESS, cplDir);

          pnlt_PhiM_PsiM_cont->SetEntities(actSDList, actSDList);
          flux_DPhiM_PsiM_cont->SetEntities(actSDList, actSDList);
          flux_PhiM_DPsiM_cont->SetEntities(actSDList, actSDList);
          pnlt_PhiM_PsiS_cont->SetEntities(actSDList, actSDList);
          flux_DPhiM_PsiS_cont->SetEntities(actSDList, actSDList);
          pnlt_PhiS_PsiM_cont->SetEntities(actSDList, actSDList);
          flux_PhiS_DPsiM_cont->SetEntities(actSDList, actSDList);
          pnlt_PhiS_PsiS_cont->SetEntities(actSDList, actSDList);

          pnlt_PhiM_PsiM_cont->SetFeFunctions(feFunctions_[ELEC_POTENTIAL], feFunctions_[ELEC_POTENTIAL]);
          flux_DPhiM_PsiM_cont->SetFeFunctions(feFunctions_[ELEC_POTENTIAL], feFunctions_[ELEC_POTENTIAL]);
          flux_PhiM_DPsiM_cont->SetFeFunctions(feFunctions_[ELEC_POTENTIAL], feFunctions_[ELEC_POTENTIAL]);
          pnlt_PhiM_PsiS_cont->SetFeFunctions(feFunctions_[ELEC_POTENTIAL], feFunctions_[ELEC_POTENTIAL]);
          flux_DPhiM_PsiS_cont->SetFeFunctions(feFunctions_[ELEC_POTENTIAL], feFunctions_[ELEC_POTENTIAL]);
          pnlt_PhiS_PsiM_cont->SetFeFunctions(feFunctions_[ELEC_POTENTIAL], feFunctions_[ELEC_POTENTIAL]);
          flux_PhiS_DPsiM_cont->SetFeFunctions(feFunctions_[ELEC_POTENTIAL], feFunctions_[ELEC_POTENTIAL]);
          pnlt_PhiS_PsiS_cont->SetFeFunctions(feFunctions_[ELEC_POTENTIAL], feFunctions_[ELEC_POTENTIAL]);

          assemble_->AddBiLinearForm(pnlt_PhiM_PsiM_cont);
          assemble_->AddBiLinearForm(flux_DPhiM_PsiM_cont);
          assemble_->AddBiLinearForm(flux_PhiM_DPsiM_cont);
          assemble_->AddBiLinearForm(pnlt_PhiM_PsiS_cont);
          assemble_->AddBiLinearForm(flux_DPhiM_PsiS_cont);
          assemble_->AddBiLinearForm(pnlt_PhiS_PsiM_cont);
          assemble_->AddBiLinearForm(flux_PhiS_DPsiM_cont);
          assemble_->AddBiLinearForm(pnlt_PhiS_PsiS_cont);
        } // end nitsche
        else if (formulation == "Mortar")
        {
          shared_ptr<SurfElemList> surfMasterGrid(new SurfElemList(ptGrid_)), surfSlaveGrid(new SurfElemList(ptGrid_));
          surfMasterGrid->SetRegion(mortarIf->GetMasterSurfRegion());
          surfSlaveGrid->SetRegion(mortarIf->GetSlaveSurfRegion());

          // --- Set the approximation for Lagrange Multipliers for the current region ---
          RegionIdType regId = surfSlaveGrid->GetRegion();
          std::string polyId = regionsList[j]->Get("polyId")->As<std::string>();
          std::string integId = regionsList[j]->Get("integId")->As<std::string>();
          feFunctions_[LAGRANGE_MULT]->GetFeSpace()->SetRegionApproximation(regId, polyId, integId);

          BiLinearForm *intOne1 = 0;
          BiLinearForm *intFactor1 = 0;
          BiLinearForm *intOne2 = 0;
          BiLinearForm *intFactor2 = 0;

          if (dim_ == 2)
          {
            intOne1 = new SurfaceMortarABInt<Complex, Complex>(new IdentityOperator<FeH1, 2, 1, Complex>(),
                                                               new IdentityOperator<FeH1, 2, 1, Complex>(),
                                                               one, pz,
                                                               mortarIf->GetSlaveVolRegion(), mortarIf->GetMasterVolRegion(),
                                                               mortarIf->IsPlanar(), updatedGeo_, BiLinearForm::SLAVE_MASTER);
            intFactor1 = new BBInt<Complex, Complex>(new IdentityOperator<FeH1, 2, 1, Complex>(), factor, -pz, updatedGeo_);
            intOne2 = new BBInt<Complex, Complex>(new IdentityOperator<FeH1, 2, 1, Complex>(), one, pz, updatedGeo_);
            intFactor2 = new SurfaceMortarABInt<Complex, Complex>(new IdentityOperator<FeH1, 2, 1, Complex>(),
                                                                  new IdentityOperator<FeH1, 2, 1, Complex>(),
                                                                  factor, -pz,
                                                                  mortarIf->GetMasterVolRegion(), mortarIf->GetSlaveVolRegion(),
                                                                  mortarIf->IsPlanar(), updatedGeo_, BiLinearForm::MASTER_SLAVE);
          }
          else
          {
            intOne1 = new SurfaceMortarABInt<Complex, Complex>(new IdentityOperator<FeH1, 3, 1, Complex>(),
                                                               new IdentityOperator<FeH1, 3, 1, Complex>(),
                                                               one, pz,
                                                               mortarIf->GetSlaveVolRegion(), mortarIf->GetMasterVolRegion(),
                                                               mortarIf->IsPlanar(), updatedGeo_, BiLinearForm::SLAVE_MASTER);
            intFactor1 = new BBInt<Complex, Complex>(new IdentityOperator<FeH1, 3, 1, Complex>(), factor, -pz, updatedGeo_);
            intOne2 = new BBInt<Complex, Complex>(new IdentityOperator<FeH1, 3, 1, Complex>(), one, pz, updatedGeo_);
            intFactor2 = new SurfaceMortarABInt<Complex, Complex>(new IdentityOperator<FeH1, 3, 1, Complex>(),
                                                                  new IdentityOperator<FeH1, 3, 1, Complex>(),
                                                                  factor, -pz,
                                                                  mortarIf->GetMasterVolRegion(), mortarIf->GetSlaveVolRegion(),
                                                                  mortarIf->IsPlanar(), updatedGeo_, BiLinearForm::MASTER_SLAVE);
          }

          intOne1->SetName("master1Elec");
          intFactor1->SetName("slave1Elec");
          intOne2->SetName("master2Elec");
          intFactor2->SetName("slave2Elec");

          // (1) weak form of the periodic boundary conditions
          BiLinFormContext *lagPotContMaster = new BiLinFormContext(intFactor1, STIFFNESS);
          NcBiLinFormContext *lagPotContSlave = new NcBiLinFormContext(intOne1, STIFFNESS);

          lagPotContMaster->SetEntities(surfSlaveGrid, surfSlaveGrid);
          lagPotContSlave->SetEntities(mortarIf->GetElemList(), mortarIf->GetElemList());
          lagPotContMaster->SetFeFunctions(feFunctions_[LAGRANGE_MULT], feFunctions_[ELEC_POTENTIAL]);
          lagPotContSlave->SetFeFunctions(feFunctions_[LAGRANGE_MULT], feFunctions_[ELEC_POTENTIAL]);

          assemble_->AddBiLinearForm(lagPotContMaster);
          assemble_->AddBiLinearForm(lagPotContSlave);

          // (2) weak form of the boundary integrals of the PDE
          BiLinFormContext *potLagContMaster = new BiLinFormContext(intOne2, STIFFNESS);
          NcBiLinFormContext *potLagContSlave = new NcBiLinFormContext(intFactor2, STIFFNESS);

          potLagContMaster->SetEntities(surfSlaveGrid, surfSlaveGrid);
          potLagContSlave->SetEntities(mortarIf->GetElemList(), mortarIf->GetElemList());
          potLagContMaster->SetFeFunctions(feFunctions_[ELEC_POTENTIAL], feFunctions_[LAGRANGE_MULT]);
          potLagContSlave->SetFeFunctions(feFunctions_[ELEC_POTENTIAL], feFunctions_[LAGRANGE_MULT]);

          assemble_->AddBiLinearForm(potLagContMaster);
          assemble_->AddBiLinearForm(potLagContSlave);

          feFunctions_[LAGRANGE_MULT]->AddEntityList(surfSlaveGrid);
          feFunctions_[ELEC_POTENTIAL]->AddEntityList(surfMasterGrid);
          feFunctions_[ELEC_POTENTIAL]->AddEntityList(surfSlaveGrid);
        } // end mortar
        else
        {
          EXCEPTION("Unknown formulation: '" << formulation << "'!");
        }
      }
    }
  }

  void ElecPDE::DefineRhsLoadIntegrators() {
    LOG_TRACE(elecpde) << "Defining rhs load integrators for electrostatic PDE";

    // Get FESpace and FeFunction of electric potential
    shared_ptr<BaseFeFunction> myFct = feFunctions_[ELEC_POTENTIAL];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();

    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    LinearForm * lin = NULL;
    StdVector<std::string> dofNames;

    // Note; in the piezoelectric case we have to multiply by -1
     Double factor = 1.0;
     if ( isPiezoCoupled_ )
       factor = -1.0;

    
    // Flag, if coefficient function lives on updated geoemtry
    bool coefUpdateGeo = true;
    
    // =========================
    //  Charges (volume, nodal)
    // =========================
    LOG_DBG(elecpde) << "Reading charges";
    ReadRhsExcitation( "charge", dofNames, ResultInfo::VECTOR, 
                       isComplex_, ent, coef,coefUpdateGeo );

    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {

        // ---------------
        //  Nodal Charges 
        // ---------------
        // Nodal charge must be constant
        if( coef[i]->GetDependency() == CoefFunction::GENERAL ) {
          EXCEPTION("Nodal charges must not be spatial dependent");
        }

        UInt numNodes = ent[i]->GetSize();
        Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;  
        coef[i] = CoefFunction::Generate(mp_, part, 
                                         CoefXprBinOp(mp_, coef[i], 
                                         boost::lexical_cast<std::string>(numNodes*factor), CoefXpr::OP_DIV) );

        lin = new SingleEntryInt(coef[i]);
        lin->SetName("NodalChargeInt");
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( ent[i] );
        ctx->SetFeFunction(myFct);
        assemble_->AddLinearForm(ctx);
        myFct->AddEntityList(ent[i]);
      } else {
        // --------------------------
        //  Surface / Volume Charges 
        // --------------------------
        if( coef[i]->GetDependency() == CoefFunction::GENERAL ) {
          EXCEPTION("The total prescribed charge must no be spatial dependent");
        }
        // "Divide" the total charge by the volume / surface of the current entity list
        Double volume = ptGrid_->CalcVolumeOfEntityList( ent[i], false );
        Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;  
        coef[i] = CoefFunction::Generate(mp_, part, 
                   CoefXprVecScalOp(mp_, coef[i], 
                   boost::lexical_cast<std::string>(volume), CoefXpr::OP_DIV) );
        if(isComplex_) {
          lin = new BUIntegrator<Complex>( new IdentityOperator<FeH1>(),
                                           Complex(factor), coef[i],coefUpdateGeo);
        } else  {
          lin = new BUIntegrator<Double>( new IdentityOperator<FeH1>(),
                                          factor, coef[i],coefUpdateGeo);
        }
        lin->SetName("ChargeDensityInt");
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( ent[i] );
        ctx->SetFeFunction(myFct);
        assemble_->AddLinearForm(ctx);
        myFct->AddEntityList(ent[i]);
      }
    } // for

    // ================
    //  CHARGE DENSITY 
    // ================
    LOG_DBG(elecpde) << "Reading charge densities";
    ReadRhsExcitation( "chargeDensity", dofNames, 
                       ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Charge density must be defined on elements")
      }
      if(isComplex_) {
        lin = new BUIntegrator<Complex>( new IdentityOperator<FeH1>(),
                                         Complex(factor), coef[i], coefUpdateGeo);
      } else  {
        lin = new BUIntegrator<Double>( new IdentityOperator<FeH1>(),
                                        factor, coef[i], coefUpdateGeo);
      }
      lin->SetName("ChargeDensityInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
    } // for
    
    // =================================
    //  Polarisation -> from hysteresis (VOLUME)
    // =================================

    //check for hysteresis
    if ( isHysteresis_ && isHysteresisFixPoint_ == true ) {
      LOG_DBG(elecpde) << "Putting polarisation to rhs";

      std::map<RegionIdType,PtrCoefFct > regionCoefs = hysteresisCoefs_->GetRegionCoefs();
      std::map<RegionIdType, shared_ptr<CoefFunction> > ::iterator it;
      for( it = regionCoefs.begin(); it != regionCoefs.end(); it++) {

        // get regionIdType
        RegionIdType curReg = it->first;

        // get SDList
        shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
        actSDList->SetRegion( curReg );

        // currently we evaluate P only at the midpoint -> fullevaluation inside BUIntegrator has to be false
        bool fullevaluation = false;

        //NOTE: 
        //If this is red correctly, we assume an integrator of the form - \int_\Omega \varphi \nabla \cdot \vec{P} d\Omega
        //This cannot be acchieved with the code below as it really defines 
        // -\int_\Omega \nabla \cdot \varphi \cdot \vec{P} d\Omega which does not make sense.
        //The operator in BUInt is supposed to work only on the test function.
        //Two possible ways to come around this:
        // 1. Integrate by parts and replace DivOperator by GradientOperator for the scalar test
        //    function, keeping in mind the boundary term.
        // 2. Modify the coefFunction (Material?) to return the divergence directly and
        //    set an identity operator on the scalar test function.

        /* for possibility two, a test with the following can be a solution 
        it->second->SetDerivativeOperation(CoefFunction::VECTOR_DIVERGENCE);
        if(isComplex_) {
          if( dim_ == 2 ) {
          // we need -factor as we put +divP to the rhs
            lin = new BUIntegrator<Complex>( new IdentityOperator<FeH1>(),
                                           Complex(-1*factor),it->second,  coefUpdateGeo, fullevaluation);
         ...... and so on

         But make sure, that the result is as expected. The SetDerivativeOperation is only
         implemented in some CoefFunctions...
         */


        if(isComplex_) {
          if( dim_ == 2 ) {
          // we need -factor as we put +divP to the rhs
            lin = new BUIntegrator<Complex>( new DivOperator<FeH1,2,Complex>(),
                                           Complex(-1*factor),it->second,  coefUpdateGeo, fullevaluation);
          } else {
            lin = new BUIntegrator<Complex>( new DivOperator<FeH1,3,Complex>(),
                                            Complex(-1*factor),it->second,  coefUpdateGeo, fullevaluation);
          }
        } else  {
          if( dim_ == 2 ) {
            // we need -factor as we put +divP to the rhs
            lin = new BUIntegrator<Double>( new DivOperator<FeH1,2,Double>(),
                                            (-1*factor),it->second,  coefUpdateGeo, fullevaluation);
          } else {
            lin = new BUIntegrator<Double>( new DivOperator<FeH1,3,Double>(),
                                             (-1*factor),it->second,  coefUpdateGeo, fullevaluation);
          }
        }

        lin->SetName("rhs_polarization");
        lin->SetSolDependent();
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( actSDList );
        ctx->SetFeFunction(myFct);
        assemble_->AddLinearForm(ctx);
        // Add entity list will add nothing, if entities were already assigned
        //myFct->AddEntityList(actSDList);
      }
    }

    // =================================
    //  FLUX DENSITY (VOLUME / SURFACE) 
    // =================================
    LOG_DBG(elecpde) << "Reading prescribed flux density";
    StdVector<std::string> vecDofNames;
    if(dim_ == 3 || subType_== "2.5d")
    {
      vecDofNames = "x", "y", "z";
    }
    else
    {
      if(dim_ == 2 && !isaxi_)
        vecDofNames = "x", "y";
      if(dim_ == 2 && isaxi_)
        vecDofNames = "r", "z";
    }
    ReadRhsExcitation( "fluxDensity", vecDofNames, 
                       ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo );
    std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Charge density must be defined on elements")
      }

      // determine dimension
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if( elemDim == (dim_-1) ) {
        // === SURFACE ===
        if( dim_ == 2) {
          if(isComplex_) {
            lin = new BUIntegrator<Complex, true>( new IdentityOperatorNormal<FeH1,2>(),
                                                   Complex(factor), coef[i], volRegions, coefUpdateGeo);
          } else  {
            lin = new BUIntegrator<Double, true>( new IdentityOperatorNormal<FeH1,2>(),
                                                  factor, coef[i], volRegions, coefUpdateGeo);
          } 
        }else {
          if(isComplex_) {
            lin = new BUIntegrator<Complex, true>( new IdentityOperatorNormal<FeH1,3>(),
                                                   Complex(factor), coef[i], volRegions, coefUpdateGeo);
          } else  {
            lin = new BUIntegrator<Double, true>( new IdentityOperatorNormal<FeH1,3>(),
                                                  factor, coef[i], volRegions, coefUpdateGeo);
          } 
        }
        lin->SetName("SurfaceFluxDensityInt");
      } else {
        // === VOLUME ===
        EXCEPTION("Elec Flux density in the volume not yet implemented");
      }

      
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    } // for
  }

  BaseBDBInt * ElecPDE::GetStiffIntegrator( BaseMaterial* actSDMat,
                                            SubTensorType tensorType,
                                            RegionIdType regionId ) {

    BaseBDBInt * integ = NULL;
    bool isComplex = complexMatData_[regionId];

    shared_ptr<CoefFunction > curCoef;
    //get possible nonlinearities defined in this region
    StdVector<NonLinType> nonLinTypes = regionNonLinTypes_[regionId];
    if ( nonLinTypes.Find(HYSTERESIS) != -1 || nonLinTypes.Find(HYSTERESIS_FIXPOINT) != -1 ){
      /* for both the delta material method as well as the std fixpoint method we have to know
       * which regions are affected by hystersis
       */

      shared_ptr<CoefFunction > curCoef_tmp;
    	// create new entity list
    	shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
    	actSDList->SetRegion( regionId );
    	PtrCoefFct elecFieldCoef = this->GetCoefFct(ELEC_FIELD_INTENSITY);

    	curCoef_tmp.reset(new CoefFunctionHyst( actSDMat, actSDList,
						                    elecFieldCoef,tensorType,ELEC_PERMITTIVITY));

    	hysteresisCoefs_->AddRegion( regionId, curCoef_tmp);

      if ( nonLinTypes.Find(HYSTERESIS_FIXPOINT) != -1 ){
        /*
         * here we treat the case: D = eps0*E + P
         * P will be put on the rhs, for stiffness integrator we need just eps0, so we reset curCoef to eps0
         */
        std::string eps0 = "8.854187817e-12";
        StdVector<std::string> realVal = StdVector<std::string>(dim_*dim_);
        realVal.Init("0.0");
        realVal[0] = eps0;
        if(dim_ == 2){
          realVal[3] = eps0;
        } else if(dim_ == 3){
          realVal[4] = eps0;
          realVal[8] = eps0;
        }
        StdVector<std::string> imagVal = StdVector<std::string>(dim_*dim_);
        imagVal.Init("0.0");

        curCoef = CoefFunction::Generate(mp_, Global::REAL, dim_, dim_, realVal, imagVal);

        std::cout << "Using FixPoint Hystersis" << std::endl;
        std::cout << "Attention: FixPoint Hysteresis just applies Preisach to given field. Hysteresis does not influence the result! " << std::endl;

        isHysteresisFixPoint_ = true;
      } else {

        std::cout << "Using DeltaMaterial Hystersis" << std::endl;

        curCoef = curCoef_tmp;
        isHysteresisFixPoint_ = false;
      }

    }
    else {

    	if ( isComplex ) {
    		curCoef = actSDMat->GetTensorCoefFnc( ELEC_PERMITTIVITY,tensorType,
                                                  Global::COMPLEX);
    	}
    	else {
    		curCoef = actSDMat->GetTensorCoefFnc( ELEC_PERMITTIVITY,tensorType,
                                                  Global::REAL);
    	}
    }
    
    // store coefficient function for later use (e.g. in boundary integrators)
    regionPermittivity_[regionId] = curCoef;

    // Note; in the piezoelectric case we have to multiply by -1
    Double factor = 1.0;
    if ( isPiezoCoupled_ )
    	factor = -1.0;


    if ( isComplex && !isHysteresis_) {
    	if( dim_ == 2 ) {
    		integ = new BDBInt<Complex,Complex >(new GradientOperator<FeH1,2,1,Complex>(),
    				curCoef, factor, updatedGeo_ );
    	} else {
    		integ = new BDBInt<Complex,Complex >(new GradientOperator<FeH1,3,1,Complex>(),
    				curCoef, factor, updatedGeo_ );
    	}
    }
    else {
    	if( dim_ == 2 ) {
    		integ = new BDBInt<>(new GradientOperator<FeH1,2> (),
    				curCoef, factor, updatedGeo_ );
    	} else {
    		integ = new BDBInt<>(new GradientOperator<FeH1,3>(),
    				curCoef, factor, updatedGeo_ );
    	}
    }

    return integ;
  }
  
  BaseBDBInt* ElecPDE::GetStiffIntegrator(BaseMaterial* actSDMat, SubTensorType tensorType,
                                          RegionIdType regionId, PtrCoefFct scalingFactor)
  {
    BaseBDBInt* integ = NULL;

    shared_ptr<CoefFunction > curCoef;
    curCoef = actSDMat->GetTensorCoefFnc(ELEC_PERMITTIVITY, tensorType, Global::COMPLEX);
    
    // store coefficient function for later use (e.g. in boundary integrators)
    regionPermittivity_[regionId] = curCoef;
    PtrCoefFct curCoefScl = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprTensScalOp(mp_, curCoef, scalingFactor, CoefXpr::OP_MULT_TENSOR));

    // Note; in the piezoelectric case we have to multiply by -1
    Double factor = 1.0;
    if (isPiezoCoupled_)
      factor = -1.0;

    BaseBOperator* bOp = NULL;
    if (dim_ == 2)
    {
      if (subType_ == "2.5d")
        bOp = new ScaledGradientOperator2p5D<FeH1, 2, Complex>();
      else
        bOp = new ScaledGradientOperator<FeH1, 2, Complex>();
    }
    else
      bOp = new ScaledGradientOperator<FeH1, 3, Complex>();

    integ = new BDBInt<Complex, Complex>(bOp, curCoefScl, factor, updatedGeo_);

    return integ;
  }
  
  BiLinearForm* ElecPDE::GetFluxIntegrator(PtrCoefFct scalCoefFunc, PtrCoefFct coefFnc, Complex factor,
                                           BiLinearForm::CouplingDirection cplDir, bool fluxOpA)
  {
    BiLinearForm* integ = NULL;
    BaseBOperator *fluxOp = NULL, *idOp = NULL;

    // The flux operator will implement either a scaled gradient operator if a non-zero 'coefFnc' has bee passed,
    // or a gradient operator otherwise
    if (dim_ == 2)
    {
      idOp = new SurfaceIdentityOperator<FeH1, 2, 1>();
      if (coefFnc)
        fluxOp = new SurfaceNormalFluxDensityOperator<FeH1, 2, 1, Complex>(subType_, coefFnc);
      else
        fluxOp = new SurfaceNormalFluxDensityOperator<FeH1, 2, 1, Double>(subType_);
    }
    else
    {
      idOp = new SurfaceIdentityOperator<FeH1, 3, 1>();
      if (coefFnc)
        fluxOp = new SurfaceNormalFluxDensityOperator<FeH1, 3, 1, Complex>(subType_, coefFnc);
      else
        fluxOp = new SurfaceNormalFluxDensityOperator<FeH1, 3, 1, Double>(subType_);
    }

    // Check whether we have a du/dn*v or a u*dv/dn bilinear form
    if (fluxOpA)
      integ = new SurfaceNitscheABInt<Complex, Complex>(fluxOp, idOp, scalCoefFunc, factor, cplDir, updatedGeo_);
    else
      integ = new SurfaceNitscheABInt<Complex, Complex>(idOp, fluxOp, scalCoefFunc, factor, cplDir, updatedGeo_);

    return integ;
  }

  BiLinearForm* ElecPDE::GetPenaltyIntegrator(PtrCoefFct scalCoefFunc, Complex factor, BiLinearForm::CouplingDirection cplDir)
  {
    BiLinearForm* integ = NULL;

    if (dim_ == 2)
    {
      integ = new SurfaceNitscheABInt<Complex, Complex>(new SurfaceIdentityOperator<FeH1, 2, 1>(),
                                                        new SurfaceIdentityOperator<FeH1, 2, 1>(),
                                                        scalCoefFunc, factor, cplDir, updatedGeo_, false, true);
    }
    else
    {
      integ = new SurfaceNitscheABInt<Complex, Complex>(new SurfaceIdentityOperator<FeH1, 3, 1>(),
                                                        new SurfaceIdentityOperator<FeH1, 3, 1>(),
                                                        scalCoefFunc, factor, cplDir, updatedGeo_, false, true);
    }

    return integ;
  }

  void ElecPDE::DefineImpedanceIntegrators() {
    
    // The definition of the impedances is taken from:

    // Jian S. Wang and Dale F. Ostergaard:
    // A Finite Element-Electric Circuit Coupled Simulation Method for
    // Piezoelectric Transducer
    // IEEE Ultrasonics Symposium, 1999

    // =======================================================================
    // Integrators for electric impedances
    // =======================================================================
    
    for( UInt i = 0; i < impedances_.GetSize(); i++ ) {
      
      REFACTOR;
    }
  }

  void ElecPDE::DefineSolveStep() {
    solveStep_ = new StdSolveStep(*this);
  }

  void ElecPDE::InitTimeStepping() {
	  Double gamma = 1.0;
	  GLMScheme * scheme = new Trapezoidal(gamma);

	  //TimeSchemeGLM::NonLinType nlType = (nonLin_ || isHysteresisFixPoint_)? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
	  shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0) );
	  feFunctions_[ELEC_POTENTIAL]->SetTimeScheme(myScheme);
  }

  void ElecPDE::FinilizeAfterTimeStep() {

	  //check for hysteresis
	  if ( isHysteresis_ && isHysteresisFixPoint_ == false ) {
		  //set current values to previous values for hysteresis operator
		  //needed for the next time step
		  std::map<RegionIdType,PtrCoefFct > regionCoefs = hysteresisCoefs_->GetRegionCoefs();
		  std::map<RegionIdType, shared_ptr<CoefFunction> > ::iterator it;
		  for( it = regionCoefs.begin(); it != regionCoefs.end(); it++) {
			  it->second->SetPreviousHystVals();
		  }
	  }
  }

  void ElecPDE::ReadSpecialBCs( ) 
  {
     //ReadImpedances();
  }
  
  void ElecPDE::ReadImpedances( ) 
  {
    REFACTOR;
  }
  
  void ElecPDE::CalcPolarizationField( shared_ptr<BaseResult> res )
  {
    REFACTOR;
  }



  template <class TYPE>
  void ElecPDE::CalcCharges( shared_ptr<BaseResult> res ) {
    REFACTOR;
    // This will be also moved to the resultFunctor class,
    // as soon as the surface B-operators are defined
  }



  void ElecPDE::SetPiezoCoupling()
  {
  
    isPiezoCoupled_ = true;

  }
  
  

  void ElecPDE::DefinePrimaryResults() {
    
    shared_ptr<BaseFeFunction> feFct = feFunctions_[ELEC_POTENTIAL];
    
    // Electric Potential
    shared_ptr<ResultInfo> res1( new ResultInfo);
    res1->resultType = ELEC_POTENTIAL;
    
    // check for special subtype 
    if( subType_  != "flatShell" ) {
      res1->dofNames = "";
      res1->unit = "V";
      res1->definedOn = ResultInfo::NODE;
      res1->entryType = ResultInfo::SCALAR;
    } else {

      // check number of composite materials:
      // Currently, we can handle just one electrostatic composite
      // material, as we would have a different number of electric
      // unknowns for different regions
      if( compositeMaterials_.size() > 1 ) {
        WARN("Currently the electrostatic flatShell PDE can only "
            "handle ONE type of composite material!");
      }

      Composite & actComp = compositeMaterials_.begin()->second;
      UInt numLaminas = actComp.thickness.GetSize();
      for( UInt i=0; i<numLaminas; i++ ) {
        std::string dofName = "ep";
        dofName += lexical_cast<std::string>(i+1);
        res1->dofNames.Push_back( dofName );
      }
      res1->unit = "V";
      res1->definedOn = ResultInfo::ELEMENT;
      res1->entryType = ResultInfo::SCALAR;
    }
    feFunctions_[ELEC_POTENTIAL]->SetResultInfo(res1);

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[ELEC_POTENTIAL] = "ground";
    idbcSolNameMap_[ELEC_POTENTIAL] = "potential";
    
    res1->SetFeFunction(feFunctions_[ELEC_POTENTIAL]);
    results_.Push_back( res1 );
    DefineFieldResult( feFunctions_[ELEC_POTENTIAL], res1 );
    
    // Electric RHS Load
    // create new resultDof object
    shared_ptr<ResultInfo> rhs ( new ResultInfo );
    rhs->resultType = ELEC_RHS_LOAD;
    rhs->dofNames = "";
    rhs->unit = "C";
    rhs->definedOn = results_[0]->definedOn;
    rhs->entryType = ResultInfo::SCALAR;
    rhsFeFunctions_[ELEC_POTENTIAL]->SetResultInfo(rhs);
    availResults_.insert( rhs );
    DefineFieldResult( rhsFeFunctions_[ELEC_POTENTIAL], rhs );

    // Check if an FE space for Lagrange multiplier has been created.
    // If any, define results for Lagrange multiplier in case of p.b.c.
    PtrParamNode lagSpaceNode = infoNode_->Get("feSpaces", ParamNode::PASS)->Get("lagrangeMultiplier", ParamNode::PASS);
    if (lagSpaceNode)
    {
      // <!> This is a hack, because there is no cross points handling
      hdbcSolNameMap_[LAGRANGE_MULT] = "ground";
      idbcSolNameMap_[LAGRANGE_MULT] = "potential";
      //
      shared_ptr<ResultInfo> lagMultElec(new ResultInfo);
      lagMultElec->resultType = LAGRANGE_MULT;
      lagMultElec->dofNames = "";
      lagMultElec->unit = "C/m^2";
      lagMultElec->entryType = ResultInfo::SCALAR;
      lagMultElec->SetFeFunction(feFunctions_[LAGRANGE_MULT]);
      lagMultElec->definedOn = ResultInfo::NODE;
      feFunctions_[LAGRANGE_MULT]->SetResultInfo(lagMultElec);

      results_.Push_back(lagMultElec);
      DefineFieldResult(feFunctions_[LAGRANGE_MULT], lagMultElec);
    }
  }
  
  void ElecPDE::DefineNcIntegrators() {
    StdVector< NcInterfaceInfo >::iterator ncIt = ncInterfaces_.Begin(),
                                           endIt = ncInterfaces_.End();
    for ( ; ncIt != endIt; ++ncIt ) {
      switch (ncIt->type) {
      case NC_MORTAR:
        if (dim_ == 2)
          DefineMortarCoupling<2,1>(ELEC_POTENTIAL, *ncIt);
        else
          DefineMortarCoupling<3,1>(ELEC_POTENTIAL, *ncIt);
        break;
      case NC_NITSCHE:
        EXCEPTION("ncInterface of Nitsche type is not implemented for ElecPDE");
        break;
      default:
        EXCEPTION("Unknown type of ncInterface");
        break;
      }
    }
  }
  
  void ElecPDE::DefinePostProcResults() {

    shared_ptr<BaseFeFunction> feFct = feFunctions_[ELEC_POTENTIAL];
    bool is2p5 = (subType_ == "2.5d");
    
    // === ELECTRIC FIELD INTENSITY ===
    shared_ptr<ResultInfo> ef ( new ResultInfo );
    ef->resultType = ELEC_FIELD_INTENSITY;
    ef->SetVectorDOFs(dim_, isaxi_, is2p5);
    ef->unit = "V/m";
    ef->definedOn = ResultInfo::ELEMENT;
    ef->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionFormBased> eFunc;
    if( isComplex_ ) {
      eFunc.reset(new CoefFunctionBOp<Complex>(feFct, ef, -1.0));
    } else {
      eFunc.reset(new CoefFunctionBOp<Double>(feFct, ef, -1.0));
    }
    DefineFieldResult( eFunc, ef );
    stiffFormCoefs_.insert(eFunc);
    
    // === ELECTRIC FLUX DENSITY ===
    shared_ptr<ResultInfo> flux ( new ResultInfo );
    flux->resultType = ELEC_FLUX_DENSITY;
    flux->SetVectorDOFs(dim_, isaxi_, is2p5);
    flux->unit = "C/m^2";
    flux->definedOn = ResultInfo::ELEMENT;
    flux->entryType = ResultInfo::VECTOR;
    shared_ptr<CoefFunctionFormBased> fluxFunc;
    if( isComplex_ ) {
      fluxFunc.reset(new CoefFunctionFlux<Complex>(feFct, flux, -1.0));
    } else {
      fluxFunc.reset(new CoefFunctionFlux<Double>(feFct, flux, -1.0));
    }
    DefineFieldResult( fluxFunc, flux );
    stiffFormCoefs_.insert(fluxFunc);

    if ( isHysteresis_){
	   hysteresisCoefs_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
	   shared_ptr<ResultInfo> elecP ( new ResultInfo );
	   elecP->resultType = ELEC_POLARIZATION;
	   elecP->SetVectorDOFs(dim_, isaxi_, is2p5);
	   elecP->unit = "C/m^2";
	   elecP->definedOn = ResultInfo::ELEMENT;
	   elecP->entryType = ResultInfo::VECTOR;
	   DefineFieldResult( hysteresisCoefs_, elecP );
	   availResults_.insert( elecP );
    }

    // === ELECTRIC SURFACE CHARGE DENSITY ===
    shared_ptr<ResultInfo> chargeD(new ResultInfo);
    chargeD->resultType = ELEC_CHARGE_DENSITY;
    chargeD->dofNames = "";
    chargeD->unit = "C/m^2";
    chargeD->definedOn = ResultInfo::SURF_ELEM;
    chargeD->entryType = ResultInfo::SCALAR;
    availResults_.insert( chargeD );
    // the coefficient function is defined later
    // Note: The positive normal direction in this case is defined as the
    //       inward facing one. 
    shared_ptr<CoefFunctionSurf> sChargeDens(new CoefFunctionSurf(true, -1.0, chargeD));
    DefineFieldResult( sChargeDens, chargeD);
    surfCoefFcts_[sChargeDens] = fluxFunc;
    
    // === TOTAL ELECTRIC CHARGE ===
    shared_ptr<ResultInfo> charge(new ResultInfo);
    charge->resultType = ELEC_CHARGE;
    charge->dofNames = "";
    charge->unit = "C";
    charge->definedOn = ResultInfo::SURF_REGION;
    charge->entryType = ResultInfo::SCALAR;
    availResults_.insert( charge );
    // build result functor for integration
    shared_ptr<ResultFunctor> chargeFunc;
    if( isComplex_ ) {
      chargeFunc.reset(new ResultFunctorIntegrate<Complex>(sChargeDens, feFct, charge ) );
    } else {
      chargeFunc.reset(new ResultFunctorIntegrate<Double>(sChargeDens, feFct, charge ) );
    }
    resultFunctors_[ELEC_CHARGE] = chargeFunc;

    // === ELECTRIC ENERGY DENSITY ===
    shared_ptr<ResultInfo> ed ( new ResultInfo );
    ed->resultType = ELEC_ENERGY_DENSITY;
    ed->dofNames = "";
    ed->unit = "Ws/m^3";
    ed->definedOn = ResultInfo::ELEMENT;
    ed->entryType = ResultInfo::SCALAR;
    shared_ptr<CoefFunctionFormBased> edFunc;
    
    // for both BdBKernel and EnergyResultFunctor, we need to apply the -1 factor
    // to get right sign in the results (even though the energy results are not really usable in the coupled case as they neglect the influnce of the coupled pde)
    Double factor = 1.0;
    if ( isPiezoCoupled_ ){
      factor = -1.0;
    }
    
    if( isComplex_ ) {
      edFunc.reset(new CoefFunctionBdBKernel<Complex>(feFct, factor*0.5));
    } else {
      edFunc.reset(new CoefFunctionBdBKernel<Double>(feFct, factor*0.5));
    }
    
    DefineFieldResult( edFunc, ed );
    stiffFormCoefs_.insert(edFunc);
    
    // Electric energy
    shared_ptr<ResultInfo> energy( new ResultInfo );
    energy->resultType = ELEC_ENERGY;
    energy->definedOn = ResultInfo::REGION;
    energy->entryType = ResultInfo::SCALAR;
    energy->dofNames = "";
    energy->unit = "Ws";
    availResults_.insert( energy );
    shared_ptr<ResultFunctor> energyFunc;
    if( isComplex_ ) {
      energyFunc.reset(new EnergyResultFunctor<Complex>(feFct, energy,factor*0.5));
    } else {
      energyFunc.reset(new EnergyResultFunctor<Double>(feFct, energy,factor*0.5));
    }
    resultFunctors_[ELEC_ENERGY] = energyFunc;
    stiffFormFunctors_.insert(energyFunc);
    
  }

  
  std::map<SolutionType, shared_ptr<FeSpace> >
  ElecPDE::CreateFeSpaces(const std::string& formulation, PtrParamNode infoNode) {
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    if(formulation == "default" || formulation == "H1"){
      PtrParamNode potSpaceNode = infoNode->Get("elecPotential");
      crSpaces[ELEC_POTENTIAL] =
        FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[ELEC_POTENTIAL]->Init(solStrat_);

      CoupledField::ParamNodeList blochListMortar = myParam_->Get("bcsAndLoads")->GetListByVal("blochPeriodic", "formulation", "Mortar");
      if (blochListMortar.GetSize())
      {
        // Create FE-Space for Lagrange multiplier
        PtrParamNode lagSpaceNode = infoNode->Get("lagrangeMultiplier");
        crSpaces[LAGRANGE_MULT] = FeSpace::CreateInstance(myParam_, lagSpaceNode, FeSpace::H1, ptGrid_);
        crSpaces[LAGRANGE_MULT]->SetLagrSurfSpace();
        crSpaces[LAGRANGE_MULT]->Init(solStrat_);
      }
    }else{
      EXCEPTION("The formulation " << formulation << "of electric PDE is not known!");
    }
    return crSpaces;
  }
}
