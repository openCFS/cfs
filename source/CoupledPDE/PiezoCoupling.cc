// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "PiezoCoupling.hh"

#include "PDE/SinglePDE.hh"
#include "CoupledPDE/BasePairCoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Enum.hh"
#include "Materials/BaseMaterial.hh"
#include "Driver/FormsContexts.hh"
#include "Driver/Assemble.hh"

#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "Domain/CoefFunction/CoefFunctionPML.hh"
#include "Domain/Mesh/NcInterfaces/MortarInterface.hh"

//transient simulations
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Driver/TransientDriver.hh"

// include fespaces
#include "FeBasis/H1/H1Elems.hh"

// new integrator concept
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/BiLinForms/ADBInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/Operators/GradientOperator.hh"
#include "Forms/Operators/StrainOperator.hh"
#include "Forms/Operators/SurfaceOperators.hh"
#include "Forms/Operators/SurfaceNormalStressOperator.hh"
#include "Forms/Operators/SurfaceNormalFluxDensityOperator.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"

namespace CoupledField {
  // ***************
  //   Constructor
  // ***************
  PiezoCoupling::PiezoCoupling( SinglePDE *pde1, SinglePDE *pde2,
          PtrParamNode paramNode,
          PtrParamNode infoNode,
          shared_ptr<SimState> simState,
          Domain* domain)
  : BasePairCoupling( pde1, pde2, paramNode, infoNode, simState, domain )
  {
    //		std::cout << "PiezoCoupling - Constructor" << std::endl;
		
    couplingName_ = "piezoDirect";
    materialClass_ = PIEZO;

    // determine subtype from mechanic pde
    pde1_->GetParamNode()->GetValue( "subType", subType_ );

    // at this point, the single pdes are not set up yet
    // > we cannot find out if something is nonlinear or not
    // > postpone to definePostProcResults
    nonLin_ = false;
    isHyst_ = false;
    irrStrainsSet_ = false;
    
    // Initialize nonlinearities
    InitNonLin();
  }
   // **************
  //   Destructor
  // **************
  PiezoCoupling::~PiezoCoupling() {
  }
  // *********************
  //   DefineIntegrators
  // *********************
  void PiezoCoupling::DefineIntegrators() {

    // at this point, the single PDEs have been set up and we can no check
    // if we are nonlinear or hysteretic
    nonLin_ = pde1_->IsNonLin() || pde2_->IsNonLin();
    // check if elec pde is hysteretic
    isHyst_ = false;
    if(pde1_->IsHysteresis()){
      EXCEPTION("Currently only the elec PDE may be hysteretic");
    } else if(pde2_->IsHysteresis()){
      isHyst_ = true;
    }
    // get hold of both feFunctions
    shared_ptr<BaseFeFunction> dispFct = pde1_->GetFeFunction(MECH_DISPLACEMENT);
    shared_ptr<BaseFeFunction> elecFct = pde2_->GetFeFunction(ELEC_POTENTIAL);
    shared_ptr<FeSpace> dispSpace = dispFct->GetFeSpace();
    shared_ptr<FeSpace> elecSpace = elecFct->GetFeSpace();
    // Global::ComplexPart matType = Global::REAL;
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;
    
    // Define integrators for "standard" materials
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    //flag indicating frequency PML formulation
    bool harmonicPML = false;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      
      // Set current region and material
      actRegion = it->first;
      actSDMat = it->second;
      // matType = Global::REAL;

      //transform the type
      SubTensorType tensorType;
      String2Enum(subType_,tensorType);
           // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );
      
      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region", "name", regionName.c_str());
      // Take account for pml in frequency domain
      // 'coeffPMLScal' is the function, the material tensor is to be scaled by. if PML isn't defined, it's unity
      PtrCoefFct coefPMLScal, coefPMLVec;
      PtrCoefFct speedOfSnd;
      PtrParamNode pmlNode;
      //
      if ((pde1_->GetDamping(actRegion) == PML) && (pde1_->GetDamping(actRegion) == PML))
      {
        if ((pde1_->GetAnalysisType() == BasePDE::HARMONIC) && (pde1_->GetAnalysisType() == BasePDE::HARMONIC))
        {
          MathParser* mp = pde1_->GetDomain()->GetMathParser();
          std::string pmlFormul;
          harmonicPML = true;
          std::string dampId;
          curRegNode->GetValue("dampingId", dampId);
          pmlNode = pde1_->GetParamNode()->Get("dampingList")->GetByVal("pml", "id", dampId.c_str());
          pmlFormul = pmlNode->Get("formulation")->As<std::string>();
          // speed of sound is set to equal '1.0'
          speedOfSnd = CoefFunction::Generate(mp, Global::REAL, "1.0");
          if (pmlFormul == "classic")
          {
            coefPMLScal.reset(new CoefFunctionPML<Complex>(pmlNode, speedOfSnd, actSDList, pde1_->GetRegions(), false));
            coefPMLVec.reset(new CoefFunctionPML<Complex>(pmlNode, speedOfSnd, actSDList, pde1_->GetRegions(), true));
          }
          else if (pmlFormul == "shifted")
          {
            coefPMLScal.reset(new CoefFunctionShiftedPML<Complex>(pmlNode, speedOfSnd, actSDList, pde1_->GetRegions(), false));
            coefPMLVec.reset(new CoefFunctionShiftedPML<Complex>(pmlNode, speedOfSnd, actSDList, pde1_->GetRegions(), true));
          }
          else{
            EXCEPTION("Unknown PML-formulation '" << pmlFormul << "'")
					}
        }
        else{
          EXCEPTION("Not implemented yet");
				}
      }
      else
      {
        harmonicPML = false;
      }

			if(isHyst_){
        //				std::cout << "Hyst case -> check if region is hyst" << std::endl;

				BaseBDBInt* mechToElecInt = NULL;
				BaseBDBInt* elecToMechInt = NULL;
				
				bool regionIsHyst = GetStiffIntegratorHyst( actSDMat,actRegion,complexMatData_[actRegion],&elecToMechInt, &mechToElecInt);
				
				if(regionIsHyst){
          //					std::cout << "Region is hyst" << std::endl;
					
          //					if(mechToElecInt == NULL){
          //						std::cout << "mechToElecInt == NULL!!" << std::endl;
          //					}
          //					if(elecToMechInt == NULL){
          //						std::cout << "elecToMechInt == NULL!!" << std::endl;
          //					}
					
					elecToMechInt->SetName("PiezoCouplingIntElecToMech");
					BiLinFormContext * stiffIntDescrElecMech =
                  new BiLinFormContext(elecToMechInt, STIFFNESS );

					stiffIntDescrElecMech->SetEntities( actSDList, actSDList );
					stiffIntDescrElecMech->SetFeFunctions( dispFct, elecFct );
					stiffIntDescrElecMech->SetCounterPart(false);

					assemble_->AddBiLinearForm( stiffIntDescrElecMech );
          //					std::cout << "first part put into place" << std::endl;
					// remember own bilinearform 
					bdbInts_[actRegion] = elecToMechInt;
					
					mechToElecInt->SetName("PiezoCouplingIntMechToElec");
					BiLinFormContext * stiffIntDescrMechElec =
                  new BiLinFormContext(mechToElecInt, STIFFNESS );

					stiffIntDescrMechElec->SetEntities( actSDList, actSDList );
					stiffIntDescrMechElec->SetFeFunctions( elecFct, dispFct );
					stiffIntDescrMechElec->SetCounterPart(false);

					assemble_->AddBiLinearForm( stiffIntDescrMechElec );

					// remember own bilinearform 
					bdbIntsCounterpart_[actRegion] = mechToElecInt;
					considerCounterpart_ = true;
          //					std::cout << "everything put into place" << std::endl;
				} else {
					// std linear case
					elecToMechInt->SetName("PiezoCouplingInt");
					BiLinFormContext * stiffIntDescr =
                  new BiLinFormContext(elecToMechInt, STIFFNESS );

					stiffIntDescr->SetEntities( actSDList, actSDList );
					stiffIntDescr->SetFeFunctions( dispFct, elecFct );
					stiffIntDescr->SetCounterPart(true);

					assemble_->AddBiLinearForm( stiffIntDescr );
						// remember own bilinearform
					bdbInts_[actRegion] = elecToMechInt;
				}
				
			} else {
				// ==================================================================
				//  STANDARD COUPLING INTEGRATOR
				// ==================================================================
				BaseBDBInt * stiffInt = NULL;
				if (harmonicPML)
				{
					stiffInt = GetStiffIntegrator(actSDMat, actRegion, true, coefPMLScal);
					// we need to set the coefficient function 'coefPMLVec' to both A-operator = ScaledStrainOp
					// and B-operator = ScaledGradOp of the bilinear form 'stiffInt'
					stiffInt->SetBCoefFunctionOpA(coefPMLVec);
					stiffInt->GetBOp()->SetCoefFunction(coefPMLVec);
				}
				else
				{
					stiffInt = GetStiffIntegrator(actSDMat, actRegion, complexMatData_[actRegion]);
				}
				stiffInt->SetName("PiezoCouplingInt");
				BiLinFormContext * stiffIntDescr =
                new BiLinFormContext(stiffInt, STIFFNESS );

				stiffIntDescr->SetEntities( actSDList, actSDList );
				stiffIntDescr->SetFeFunctions( dispFct, elecFct );
				stiffIntDescr->SetCounterPart(true);

				assemble_->AddBiLinearForm( stiffIntDescr );

				// remember own bilinearform 
				bdbInts_[actRegion] = stiffInt;

			}      
    }

    // handling periodic boundary conditions (if any)
    DefinePBCIntegrators(dispFct, elecFct);
		
		DefineRhsLoadIntegrators();
  }
  void PiezoCoupling::DefinePBCIntegrators(shared_ptr<BaseFeFunction>& fe1, shared_ptr<BaseFeFunction>& fe2)
  {
    PtrParamNode bcNode = myParam_->Get("bcsAndLoads", ParamNode::PASS);
    this->ptGrid_->MapEdges();
    this->ptGrid_->MapFaces();
    if (!bcNode)
      return;

    MathParser* mp = pde1_->GetDomain()->GetMathParser();
       ParamNodeList blochNodesList = bcNode->GetList("blochPeriodic");
    for (UInt i = 0; i < blochNodesList.GetSize(); i++)
    {
      std::string str_value = blochNodesList[i]->Get("factor_value")->As<std::string>();
      std::string str_phase = blochNodesList[i]->Get("factor_phase")->As<std::string>();
      // propagation factor \gamma from xml-file
      std::string str_real, str_imag;
      str_real = AmplPhaseToReal(str_value, str_phase, true);
      str_imag = AmplPhaseToImag(str_value, str_phase, true);

      PtrCoefFct factor = CoefFunction::Generate(mp, Global::COMPLEX, str_real, str_imag);
      PtrCoefFct one = CoefFunction::Generate(mp, Global::REAL, "1.0", "0.0");
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

        PtrCoefFct matData, tmpData;
        RegionIdType volMasterId = mortarIf->GetPrimaryVolRegion();

        SubTensorType tensorType = NO_TENSOR;
        if( subType_ == "planeStrain" || subType_ == "planeStress" ) {
          tensorType = PLANE_STRAIN;
        } else if( subType_ == "axi" ){
          tensorType = AXI;
        } else if (subType_ == "3d" || subType_ == "2.5d" ){
          tensorType = FULL;
        } else {
          EXCEPTION( "Unknown subtype '" << subType_ << "'" );
        }

        std::string formulation = regionsList[j]->Get("formulation")->As<std::string>();
        if (formulation == "Nitsche")
        {
          std::string nitFac = regionsList[j]->Get("nitscheFactor")->As<std::string>();
          shared_ptr<ElemList> actSDList = ncIf->GetElemList();

          PtrCoefFct coefFuncPMLVec, coefFuncPMLScl;
          if ((pde1_->GetDamping(volMasterId) == PML) && (pde2_->GetDamping(volMasterId) == PML))
          {
            if ((pde1_->GetAnalysisType() == BasePDE::HARMONIC) && (pde1_->GetAnalysisType() == BasePDE::HARMONIC))
            {
              std::string regionName = ptGrid_->GetRegion().ToString(volMasterId);
              std::string dampId, pmlFormul;

              PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region", "name", regionName.c_str());
              curRegNode->GetValue("dampingId", dampId);
              PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal("pml", "id", dampId.c_str());

              // speed of sound is set to equal '1.0'
              PtrCoefFct speedOfSnd = CoefFunction::Generate(mp, Global::REAL, "1.0");
              pmlFormul = pmlNode->Get("formulation")->As<std::string>();

              if (pmlFormul == "classic")
              {
                coefFuncPMLVec.reset(new CoefFunctionPML<Complex>(pmlNode, speedOfSnd,
                        ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), pde1_->GetRegions(), true));
                coefFuncPMLScl.reset(new CoefFunctionPML<Complex>(pmlNode, speedOfSnd,
                        ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), pde1_->GetRegions(), false));
              }
              else if (pmlFormul == "shifted")
              {
                coefFuncPMLVec.reset(new CoefFunctionShiftedPML<Complex>(pmlNode, speedOfSnd,
                        ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), pde1_->GetRegions(), true));
                coefFuncPMLScl.reset(new CoefFunctionShiftedPML<Complex>(pmlNode, speedOfSnd,
                        ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), pde1_->GetRegions(), false));
              }
              else
              {
                EXCEPTION("Unknown PML-formulation '" << pmlFormul << "'");
              }

              tmpData = materials_[volMasterId]->GetTensorCoefFnc(PIEZO_TENSOR, tensorType, Global::COMPLEX, true);
              matData = CoefFunction::Generate(mp, Global::COMPLEX, CoefXprTensScalOp(mp, tmpData, coefFuncPMLScl, CoefXpr::OP_MULT));
            }
            else
              matData = materials_[volMasterId]->GetTensorCoefFnc(PIEZO_TENSOR, tensorType, Global::REAL, true);
          }

          // mechanical displacement (fe1) as an unknown and electric potential (fe2) as a test function
          BiLinearForm* int_PhiM_DUM = NULL; // PM' * n([e]BmUM)
          BiLinearForm* int_PhiS_DUM = NULL; // PS' * n([e]BmUM)
          BiLinearForm* int_DPhiM_UM = NULL; // n([e`]BePM') * UM
          BiLinearForm* int_DPhiM_US = NULL; // n([e`]BePM') * US

          // define bilinear forms coupling mechanical displacement with electric potential
          if (matData->IsComplex())
          {
            int_PhiM_DUM = GetNormalPiezoFluxIntegrator<Complex>(one, coefFuncPMLVec, -1.0, BiLinearForm::PRIM_PRIM, false);
            int_PhiS_DUM = GetNormalPiezoFluxIntegrator<Complex>(factor, coefFuncPMLVec, 1.0, BiLinearForm::SEC_PRIM, false);
            int_DPhiM_UM = GetNormalPiezoStrainIntegrator<Complex>(factor, coefFuncPMLVec, -1.0, BiLinearForm::PRIM_PRIM, true);
            int_DPhiM_US = GetNormalPiezoStrainIntegrator<Complex>(one, coefFuncPMLVec, 1.0, BiLinearForm::PRIM_SEC, true);
          }
          else
          {
            int_PhiM_DUM = GetNormalPiezoFluxIntegrator<Double>(one, coefFuncPMLVec, -1.0, BiLinearForm::PRIM_PRIM, false);
            int_PhiS_DUM = GetNormalPiezoFluxIntegrator<Complex>(factor, coefFuncPMLVec, 1.0, BiLinearForm::SEC_PRIM, false);
            int_DPhiM_UM = GetNormalPiezoStrainIntegrator<Complex>(factor, coefFuncPMLVec, -1.0, BiLinearForm::PRIM_PRIM, true);
            int_DPhiM_US = GetNormalPiezoStrainIntegrator<Double>(one, coefFuncPMLVec, 1.0, BiLinearForm::PRIM_SEC, true);
          }

          int_PhiM_DUM->SetBCoefFunctionOpB(matData);
          int_PhiM_DUM->SetName("int_PhiM_DUM");
          int_PhiS_DUM->SetBCoefFunctionOpB(matData);
          int_PhiS_DUM->SetName("int_PhiS_DUM");
          int_DPhiM_UM->SetBCoefFunctionOpA(matData);
          int_DPhiM_UM->SetName("int_DPhiM_UM");
          int_DPhiM_US->SetBCoefFunctionOpA(matData);
          int_DPhiM_US->SetName("int_DPhiM_US");

          // define contexts for bilinear forms
          SurfaceBiLinFormContext* descr_PhiM_DUM = new SurfaceBiLinFormContext(int_PhiM_DUM, STIFFNESS, BiLinearForm::PRIM_PRIM);
          SurfaceBiLinFormContext* descr_PhiS_DUM = new SurfaceBiLinFormContext(int_PhiS_DUM, STIFFNESS, BiLinearForm::SEC_PRIM);
          SurfaceBiLinFormContext* descr_DPhiM_UM = new SurfaceBiLinFormContext(int_DPhiM_UM, STIFFNESS, BiLinearForm::PRIM_PRIM);
          SurfaceBiLinFormContext* descr_DPhiM_US = new SurfaceBiLinFormContext(int_DPhiM_US, STIFFNESS, BiLinearForm::PRIM_SEC);

          descr_PhiM_DUM->SetEntities(actSDList, actSDList);
          descr_PhiS_DUM->SetEntities(actSDList, actSDList);
          descr_DPhiM_UM->SetEntities(actSDList, actSDList);
          descr_DPhiM_US->SetEntities(actSDList, actSDList);
          descr_PhiM_DUM->SetFeFunctions(fe2, fe1);
          descr_PhiS_DUM->SetFeFunctions(fe2, fe1);
          descr_DPhiM_UM->SetFeFunctions(fe2, fe1);
          descr_DPhiM_US->SetFeFunctions(fe2, fe1);
          assemble_->AddBiLinearForm(descr_PhiM_DUM);
          assemble_->AddBiLinearForm(descr_PhiS_DUM);
          assemble_->AddBiLinearForm(descr_DPhiM_UM);
          assemble_->AddBiLinearForm(descr_DPhiM_US);

          // electric potential (fe2) as an unknown and mechanical displacement (fe1) as a test function
          BiLinearForm* int_UM_DPhiM = NULL; // UM' * n([e`]BePM)
          BiLinearForm* int_US_DPhiM = NULL; // US' * n([e`]BePM)
          BiLinearForm* int_DUM_PhiM = NULL; // n([e]BmUM') * PM
          BiLinearForm* int_DUM_PhiS = NULL; // n([e]BmUM') * PS

          // define bilinear forms coupling electric potential with mechanical displacement
          if (matData->IsComplex())
          {
            int_UM_DPhiM = GetNormalPiezoStrainIntegrator<Complex>(one, coefFuncPMLVec, -1.0, BiLinearForm::PRIM_PRIM, false);
            int_US_DPhiM = GetNormalPiezoStrainIntegrator<Complex>(factor, coefFuncPMLVec, 1.0, BiLinearForm::SEC_PRIM, false);
            int_DUM_PhiM = GetNormalPiezoFluxIntegrator<Complex>(factor, coefFuncPMLVec, -1.0, BiLinearForm::PRIM_PRIM, true);
            int_DUM_PhiS = GetNormalPiezoFluxIntegrator<Complex>(one, coefFuncPMLVec, 1.0, BiLinearForm::PRIM_SEC, true);
          }
          else
          {
            int_UM_DPhiM = GetNormalPiezoStrainIntegrator<Double>(one, coefFuncPMLVec, -1.0, BiLinearForm::PRIM_PRIM, false);
            int_US_DPhiM = GetNormalPiezoStrainIntegrator<Complex>(factor, coefFuncPMLVec, 1.0, BiLinearForm::SEC_PRIM, false);
            int_DUM_PhiM = GetNormalPiezoFluxIntegrator<Complex>(factor, coefFuncPMLVec, -1.0, BiLinearForm::PRIM_PRIM, true);
            int_DUM_PhiS = GetNormalPiezoFluxIntegrator<Double>(one, coefFuncPMLVec, 1.0, BiLinearForm::PRIM_SEC, true);
          }


          int_UM_DPhiM->SetBCoefFunctionOpB(matData);
          int_UM_DPhiM->SetName("int_UM_DPhiM");
          int_US_DPhiM->SetBCoefFunctionOpB(matData);
          int_US_DPhiM->SetName("int_US_DPhiM");
          int_DUM_PhiM->SetBCoefFunctionOpA(matData);
          int_DUM_PhiM->SetName("int_DUM_PhiM");
          int_DUM_PhiS->SetBCoefFunctionOpA(matData);
          int_DUM_PhiS->SetName("int_DUM_PhiS");

          // define contexts for bilinear forms
          SurfaceBiLinFormContext* descr_UM_DPhiM = new SurfaceBiLinFormContext(int_UM_DPhiM, STIFFNESS, BiLinearForm::PRIM_PRIM);
          SurfaceBiLinFormContext* descr_US_DPhiM = new SurfaceBiLinFormContext(int_US_DPhiM, STIFFNESS, BiLinearForm::SEC_PRIM);
          SurfaceBiLinFormContext* descr_DUM_PhiM = new SurfaceBiLinFormContext(int_DUM_PhiM, STIFFNESS, BiLinearForm::PRIM_PRIM);
          SurfaceBiLinFormContext* descr_DUM_PhiS = new SurfaceBiLinFormContext(int_DUM_PhiS, STIFFNESS, BiLinearForm::PRIM_SEC);

          descr_UM_DPhiM->SetEntities(actSDList, actSDList);
          descr_US_DPhiM->SetEntities(actSDList, actSDList);
          descr_DUM_PhiM->SetEntities(actSDList, actSDList);
          descr_DUM_PhiS->SetEntities(actSDList, actSDList);

          descr_UM_DPhiM->SetFeFunctions(fe1, fe2);
          descr_US_DPhiM->SetFeFunctions(fe1, fe2);
          descr_DUM_PhiM->SetFeFunctions(fe1, fe2);
          descr_DUM_PhiS->SetFeFunctions(fe1, fe2);

          assemble_->AddBiLinearForm(descr_UM_DPhiM);
          assemble_->AddBiLinearForm(descr_US_DPhiM);
          assemble_->AddBiLinearForm(descr_DUM_PhiM);
          assemble_->AddBiLinearForm(descr_DUM_PhiS);
        } // end nitsche
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
  
  void PiezoCoupling::DefinePostProcResults() {
    
//		std::cout << "PiezoCoupling - DefinePostProcResults" << std::endl;
		
    shared_ptr<BaseFeFunction> dispFct = pde1_->GetFeFunction(MECH_DISPLACEMENT);
    shared_ptr<BaseFeFunction> elecFct = pde2_->GetFeFunction(ELEC_POTENTIAL);
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
    MathParser * mp = domain_->GetMathParser();
    
    StdVector<std::string> stressComponents;
    if( subType_ == "3d" || subType_ == "2.5d") {
      stressComponents = "xx", "yy", "zz", "yz", "xz", "xy";
    } else if( subType_ == "planeStrain" ) {
      stressComponents = "xx", "yy", "xy";
    } else if( subType_ == "planeStress" ) {
      stressComponents = "xx", "yy", "xy";
    } else if( subType_ == "axi" ) {
      stressComponents = "rr", "zz", "rz", "phiphi";
    }
    
    // === Electric Flux Density ===
    shared_ptr<ResultInfo> flux ( new ResultInfo );
    flux->resultType = ELEC_FLUX_DENSITY;
    flux->SetVectorDOFs(dim_, isaxi_, subType_ == "2.5d");
    flux->unit = "C/m^2";
    flux->definedOn = ResultInfo::ELEMENT;
    flux->entryType = ResultInfo::VECTOR;
    availResults_.insert( flux );
    
    // The electric flux density in the coupled case calculates as 
    // D = [e]*S + [eps]*E
    // a) electric part -> take from electrostatic PDE (should contain contributions from hysteresis, too)
    PtrCoefFct coefElecD = pde2_->GetCoefFct( ELEC_FLUX_DENSITY);
    
    // b) coupling part -> use own ADB-form
    Double cplFactor = 1.0;
    shared_ptr<CoefFunctionFormBased> cplFunc;
		// if counterpart was stored separately, it is already transposed, i.e. we
		// do not need have to tell cplFunc that we want to transpose
		bool transposed = !considerCounterpart_;
		if(transposed){
			if( isComplex_ ) {
				cplFunc.reset(new CoefFunctionFlux<Complex,true>(dispFct, flux, cplFactor));
			} else {
				cplFunc.reset(new CoefFunctionFlux<Double,true>(dispFct, flux, cplFactor));
			}
		} else {
			if( isComplex_ ) {
				cplFunc.reset(new CoefFunctionFlux<Complex>(dispFct, flux, cplFactor));
			} else {
				cplFunc.reset(new CoefFunctionFlux<Double>(dispFct, flux, cplFactor));
			}
		}

    // Build compound coefficient function for flux density
    PtrCoefFct coefFlux = CoefFunction::Generate(mp,part,
            CoefXprBinOp(mp,coefElecD, cplFunc,
            CoefXpr::OP_ADD) );
    // D = epsE + P  + eSr - eSi
    // in hyst case, we have to subtract [e]S_irr!
    if( !isHyst_){
      DefineFieldResult( coefFlux, flux );
    } else {
      if(irrStrainsSet_ == false){
        EXCEPTION("Hyst coef functions not initialized yet");
      }
      PtrCoefFct coefFlux_full = CoefFunction::Generate(mp,part,
              CoefXprBinOp(mp,coefFlux, coupledIrrStrains_,
              CoefXpr::OP_SUB) );
      DefineFieldResult( coefFlux_full, flux );
    }

    // === irrevesible strains and stresses ===
    // check if at least one region has hysteretic behavior
    if ( isHyst_){
      // Note: in elecPDE we have to setup storage for hyst output in definePostProcResults
      // as this function is called prior to DefineIntegrators
      // in pizeoCase, the postprocresults will be defined after the integrators, so
      // here we have to setup the storage earlier in GetStiffIntegratorHyst
      //      irrStrains_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, numStrainEntries,1,isComplex_));
      //      irrStresses_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, numStrainEntries,1,isComplex_));
      shared_ptr<ResultInfo> irrStrainInfo(new ResultInfo);
      irrStrainInfo->resultType = MECH_IRR_STRAIN;
      irrStrainInfo->dofNames = stressComponents;
      irrStrainInfo->unit =  "";
      irrStrainInfo->entryType = ResultInfo::TENSOR;
      // strains retruns them in voigt vector
      irrStrainInfo->definedOn = ResultInfo::ELEMENT;
      DefineFieldResult( irrStrains_, irrStrainInfo );
      availResults_.insert( irrStrainInfo );
      shared_ptr<ResultInfo> irrStressInfo(new ResultInfo);
      irrStressInfo->resultType = MECH_IRR_STRESS;
      irrStressInfo->dofNames = stressComponents;
      irrStressInfo->unit =  "N/m^2";
      irrStressInfo->entryType = ResultInfo::TENSOR;
      // strains retruns them in voigt vector
      irrStressInfo->definedOn = ResultInfo::ELEMENT;
      DefineFieldResult( irrStresses_, irrStressInfo );
      availResults_.insert( irrStressInfo );
    }
       // ==== Mechanic Stress ===
    shared_ptr<ResultInfo> stress(new ResultInfo);
    stress->resultType = MECH_STRESS;
    stress->dofNames = stressComponents;
    stress->unit =  "N/m^2";
    stress->entryType = ResultInfo::TENSOR;
    stress->definedOn = ResultInfo::ELEMENT;
    availResults_.insert( stress );
    
    // The mechanic stress calculates as 
    // [sigma] = [c]*S - [e]*E
    // a) mechanic  -> take from mechanic PDE
    PtrCoefFct coefMechSigma = pde1_->GetCoefFct( MECH_STRESS );
        
    // b) coupling part [e]*E = - [e] grad(V_p) -> use own ADB-form
    shared_ptr<CoefFunctionFormBased> stressCplFunc;
    // Note: As the "B"-operator of the ADB operator only calculates 
    //       the gradient and we need the electric field (E = - grad(V_p)),
    //       we have to multiply the entries by -1.0. However, due to the 
    //       constitutive law, an additional -1 is needed, so we just
    //       need scaling factor 1.0.
    Double stressCplFactor = 1.0;
    if( isComplex_ ) {
      stressCplFunc.reset(new CoefFunctionFlux<Complex>(elecFct, stress, stressCplFactor));
    } else {
      stressCplFunc.reset(new CoefFunctionFlux<Double>(elecFct, stress, stressCplFactor));
    }
    // in case of hysteresis, stresses are computed via
    // [sigma] = [c]*S_revesible - [e]*E
    // > from mech PDE we get [sigma] = [c]*Bu = [c]*S_total = [c]*(S_reversible + S_irreversible)
    // > subtract [c]*S_irreversible
    PtrCoefFct coefStress;
    if ( isHyst_ ){
      PtrCoefFct coefStress_full = CoefFunction::Generate(mp,part,
              CoefXprBinOp(mp,coefMechSigma, stressCplFunc,
              CoefXpr::OP_ADD) );
      coefStress = CoefFunction::Generate(mp,part,
              CoefXprBinOp(mp,coefStress_full, irrStressesVector_,
              CoefXpr::OP_SUB) );
    } else {
      coefStress = CoefFunction::Generate(mp,part,
              CoefXprBinOp(mp,coefMechSigma, stressCplFunc,
              CoefXpr::OP_ADD) );
    }

    DefineFieldResult(coefStress, stress);

    // === Electric Charge Density (surface) ===
    shared_ptr<ResultInfo> chargeD(new ResultInfo);
    chargeD->resultType = ELEC_SURFACE_CHARGE_DENSITY;
    chargeD->dofNames = "";
    chargeD->unit = MapSolTypeToUnit(ELEC_SURFACE_CHARGE_DENSITY);
    chargeD->definedOn = ResultInfo::SURF_ELEM;
    chargeD->entryType = ResultInfo::SCALAR;
    availResults_.insert( chargeD );
    // Note: The positive normal direction in this case is defined as the
    //       inward facing one.
    shared_ptr<CoefFunctionSurf> sChargeDens(new CoefFunctionSurf(true, -1.0, chargeD));
    DefineFieldResult( sChargeDens, chargeD);
    
    // === Electric Charge (integrated) ===
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
      chargeFunc.reset(new ResultFunctorIntegrate<Complex>(sChargeDens, dispFct, charge ) );
    } else {
      chargeFunc.reset(new ResultFunctorIntegrate<Double>(sChargeDens, dispFct, charge ) );
    }
    resultFunctors_[ELEC_CHARGE] = chargeFunc;
    
    
    // ============================
    // Initialize result functors:
    // ============================
    // 1) Loop over all BDB-integrators
    std::map<RegionIdType, BaseBDBInt*>::iterator it = bdbInts_.begin();
    for( ; it != bdbInts_.end(); ++it ) {
      RegionIdType region = it->first;
			
			if(considerCounterpart_){
				// treat elecMech and mechElec separately
				// bdbInts_ stores bdbIntegrator for mechanical pde (elecToMechInt)
				BaseBDBInt* bdb = it->second;
				stressCplFunc->AddIntegrator(bdb, region);
        //				std::cout << "StressCpl success" << std::endl;
				
				// bdbIntsCounterpart_ stores bdeIntegrator for electric pde (mechToElecInt)
				BaseBDBInt* bdbCounterpart = bdbIntsCounterpart_[region];
				cplFunc->AddIntegrator(bdbCounterpart, region);
        //				std::cout << "cpl success" << std::endl;
				sChargeDens->AddVolumeCoef(region, coefFlux);
			} else {
				BaseBDBInt* bdb = it->second;
				// 2) pass integrators to functors
				cplFunc->AddIntegrator(bdb, region);
				stressCplFunc->AddIntegrator(bdb, region);
				sChargeDens->AddVolumeCoef(region, coefFlux);
			}
    }
  }
  
  void PiezoCoupling::DefineAvailResults() {

    //
    //    // === ELECTRIC CHARGE ===
    //    shared_ptr<ResultInfo> charge(new ResultInfo);
    //    charge->resultType = ELEC_CHARGE;
    //    charge->dofNames = "";
    //    charge->unit = "C";
    //    charge->definedOn = ResultInfo::SURF_ELEM;
    //    charge->entryType = ResultInfo::SCALAR;
    //    charge->fctType = shared_ptr<ConstFct>(new ConstFct() );
    //    availResults_.insert( charge );
  }

	// for hysteresis we need rhs loads that have to be defined in coupling
	// as it needs informattion about both pdes
	void PiezoCoupling::DefineRhsLoadIntegrators() {

    //		std::cout << "PiezoCoupling - DefineLOADIntegrators" << std::endl;
		
		shared_ptr<BaseFeFunction> dispFct = pde1_->GetFeFunction(MECH_DISPLACEMENT);
    shared_ptr<BaseFeFunction> elecFct = pde2_->GetFeFunction(ELEC_POTENTIAL);
      LinearForm * lin = NULL;
		LinearForm * linElec = NULL;
    // Flag, if coefficient function lives on updated geoemtry
    bool coefUpdateGeo = true;
		
		SubTensorType tensorType = NO_TENSOR;
    if( subType_ == "planeStrain" || subType_ == "planeStress" ) {
      tensorType = PLANE_STRAIN;
    } else if( subType_ == "axi" ){
      tensorType = AXI;
    } else if (subType_ == "3d" || subType_ == "2.5d"){
      tensorType = FULL;
    } else {
      EXCEPTION( "Unknown subtype '" << subType_ << "'" );
    }

		if(isHyst_){
			// get all regions with hysteresis information from elecPDE
			shared_ptr<CoefFunctionMulti> hysteresisCoefs = pde2_->GetHystCoefs();
			BaseMaterial * actSDMat = NULL;
			
      std::map<RegionIdType,PtrCoefFct > regionCoefs = hysteresisCoefs->GetRegionCoefs();
      std::map<RegionIdType, shared_ptr<CoefFunction> > ::iterator it;
      for( it = regionCoefs.begin(); it != regionCoefs.end(); it++) {

        // get regionIdType
        RegionIdType curReg = it->first;
				PtrCoefFct regionHystOperator = it->second;

				actSDMat = materials_[curReg];
				
        // get SDList
        shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
        actSDList->SetRegion( curReg );
				// per default we evaluate at each integration point
				// this behavior can be changed by setting a flag in coefFncHyst
        bool fullevaluation = true;
				Double factor = 1.0;
				
				// ------------------------
				//  Obtain linear material
				// ------------------------
				shared_ptr<CoefFunction > couplingCoef;
				if( isComplex_ ) {
					couplingCoef = actSDMat->GetTensorCoefFnc(PIEZO_TENSOR, tensorType, 
                  Global::COMPLEX, false  );
				} else {
					couplingCoef = actSDMat->GetTensorCoefFnc(PIEZO_TENSOR, tensorType, 
                  Global::REAL, false );
				}	
				
        //				std::cout << "couplCoef" << std::endl;
        //				std::cout << couplingCoef->ToString() << std::endl;
				
				// extract stiffness tensor from mech pde
				std::map<RegionIdType, BaseMaterial*>  mechMat;
				mechMat = pde1_->GetMaterialData();

				shared_ptr<CoefFunction > stiffCoef;
				stiffCoef = mechMat[curReg]->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, tensorType, 
                Global::REAL, true );
				
        // important: pass information about stiffness and coupling tensor to hystOperator
        // (if not already done)
        regionHystOperator->SetElastAndCouplTensor(stiffCoef, couplingCoef);
				
				// mech RHS
				shared_ptr<CoefFunction> mechRHS = regionHystOperator->GenerateRHSCoefFnc("IrrStressForMechPDE");

				if ( isComplex_ ) {
					if( subType_ == "axi" ) {
						lin = new BUIntegrator<Complex>( new StrainOperatorAxi<FeH1,Complex>(),
                    Complex(factor),mechRHS,  coefUpdateGeo, fullevaluation);

					} else if( subType_ == "planeStrain" ) {
						lin = new BUIntegrator<Complex>( new StrainOperator2D<FeH1,Complex>(),
                    Complex(factor),mechRHS,  coefUpdateGeo, fullevaluation);

					} else if( subType_ == "planeStress" ) {
						lin = new BUIntegrator<Complex>( new StrainOperator2D<FeH1,Complex>(),
                    Complex(factor),mechRHS,  coefUpdateGeo, fullevaluation);

					} else if( subType_ == "3d") {
						lin = new BUIntegrator<Complex>( new StrainOperator3D<FeH1,Complex>(),
                    Complex(factor),mechRHS,  coefUpdateGeo, fullevaluation);

					} else if( subType_ == "2.5d") {
						lin = new BUIntegrator<Complex>( new StrainOperator2p5D<FeH1,Complex>(),
                    Complex(factor),mechRHS,  coefUpdateGeo, fullevaluation);
					} else {
						EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
					}
				} else {
					if( subType_ == "axi" ) {
						lin = new BUIntegrator<Double>( new StrainOperatorAxi<FeH1>(),
                    (factor),mechRHS, coefUpdateGeo, fullevaluation);

					} else if( subType_ == "planeStrain" ) {
						lin = new BUIntegrator<Double>( new StrainOperator2D<FeH1>(),
                    (factor),mechRHS, coefUpdateGeo, fullevaluation);

					} else if( subType_ == "planeStress" ) {
						lin = new BUIntegrator<Double>( new StrainOperator2D<FeH1>(),
                    (factor),mechRHS, coefUpdateGeo, fullevaluation);
					} else if( subType_ == "3d") {
						lin = new BUIntegrator<Double>( new StrainOperator3D<FeH1>(),
                    (factor),mechRHS, coefUpdateGeo, fullevaluation);
					} else if( subType_ == "2.5d") {
						lin = new BUIntegrator<Double>( new StrainOperator2p5D<FeH1>(),
                    (factor),mechRHS, coefUpdateGeo, fullevaluation);

					} else {
						EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
					}
				}

        lin->SetName("IrrStressForMechPDE");
        lin->SetSolDependent();
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( actSDList );
        ctx->SetFeFunction(dispFct);
        
				assemble_->AddLinearForm(ctx);
        // Add entity list will add nothing, if entities were already assigned
        dispFct->AddEntityList(actSDList);
				
        // ++++++++
				// elec RHS
        // > note: elecPDE will get multiplied by -1 in case of piezo coupling
        // > rhs term will also switch sign
        factor = -1.0;
				shared_ptr<CoefFunction> elecRHS = regionHystOperator->GenerateRHSCoefFnc("CoupledIrrStrainForElecPDE");
				
        if(isComplex_) {
          if( dim_ == 2 ) {
            //            lin = new BUIntegrator<Complex>( new GradientOperator<FeH1,2,1,Complex>(),
            //                                           Complex(factor),it->second,  coefUpdateGeo, fullevaluation);
						linElec = new BUIntegrator<Complex>( new GradientOperator<FeH1,2,1,Complex>(),
                    Complex(factor),elecRHS, coefUpdateGeo, fullevaluation);
          } else {
            linElec = new BUIntegrator<Complex>( new GradientOperator<FeH1,3,1,Complex>(),
                    Complex(factor),elecRHS, coefUpdateGeo, fullevaluation);
          }
        } else  {
          if( dim_ == 2 ) {
            linElec = new BUIntegrator<Double>( new GradientOperator<FeH1,2> (),
                    (factor),elecRHS, coefUpdateGeo, fullevaluation);
          } else {
            linElec = new BUIntegrator<Double>( new GradientOperator<FeH1,3> (),
                    (factor),elecRHS, coefUpdateGeo, fullevaluation);
          }
        }

        linElec->SetName("CoupledIrrStrainForElecPDE");
        linElec->SetSolDependent();
        LinearFormContext *ctxElec = new LinearFormContext( linElec );
        ctxElec->SetEntities( actSDList );
        ctxElec->SetFeFunction(elecFct);
        
				assemble_->AddLinearForm(ctxElec);
        // Add entity list will add nothing, if entities were already assigned
        elecFct->AddEntityList(actSDList);

      }
    }
  }

  bool PiezoCoupling::GetStiffIntegratorHyst( BaseMaterial* actSDMat,
          RegionIdType regionId, bool isComplex,
          BaseBDBInt** elecToMechInt, BaseBDBInt** mechToElecInt) {
    
    SubTensorType tensorType = NO_TENSOR;
    UInt numStrainEntries = 3;
    if( subType_ == "planeStrain" || subType_ == "planeStress" ) {
      tensorType = PLANE_STRAIN;
    } else if( subType_ == "axi" ){
      tensorType = AXI;
      numStrainEntries = 4;
    } else if (subType_ == "3d" || subType_ == "2.5d"){
      tensorType = FULL;
      numStrainEntries = 6;
    } else {
      EXCEPTION( "Unknown subtype '" << subType_ << "'" );
    }

    // setup arrays for output > has to be done here as we insert further below
    // in elecPDE this step is done in definePostProcResults, but in case of piezoCoupling
    // GetStiffIntegratorHyst is called prior to definePostProcResults
    //    if(!irrStrainsSet_){
    //      irrStrains_.reset(new CoefFunctionMulti(CoefFunction::TENSOR, numStrainEntries,1,isComplex_));
    //      irrStresses_.reset(new CoefFunctionMulti(CoefFunction::TENSOR, numStrainEntries,1,isComplex_));
    //      irrStrainsSet_ = true;
    //    }
    // setup arrays for output > has to be done here as we insert further below
    // in elecPDE this step is done in definePostProcResults, but in case of piezoCoupling
    // GetStiffIntegratorHyst is called prior to definePostProcResults
    if(!irrStrainsSet_){
      irrStrains_.reset(new CoefFunctionMulti(CoefFunction::TENSOR, dim_,dim_,isComplex_));
      irrStresses_.reset(new CoefFunctionMulti(CoefFunction::TENSOR, dim_,dim_,isComplex_));
      irrStressesVector_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, numStrainEntries,1,isComplex_));
      coupledIrrStrains_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, dim_,1,isComplex_));
      irrStrainsSet_ = true;
    }
       // ------------------------
    //  Obtain linear material
    // ------------------------		
		shared_ptr<CoefFunction > couplingCoef;
    if( isComplex ) {
      couplingCoef = actSDMat->GetTensorCoefFnc(PIEZO_TENSOR, tensorType, 
              Global::COMPLEX, false  );
    } else {
      couplingCoef = actSDMat->GetTensorCoefFnc(PIEZO_TENSOR, tensorType, 
              Global::REAL, false );
    }	
		
    //		std::cout << "couplingCoefTensor: " << couplingCoef->ToString() << std::endl;
    //
		// check if current region is hysteretic
		std::map<RegionIdType, StdVector<NonLinType> > nonLinMap;
		nonLinMap = pde2_->GetNonLinRegionTypes();

		StdVector<NonLinType> nonLinTypes = nonLinMap[regionId];
    if ( nonLinTypes.Find(HYSTERESIS) != -1 ){

			// extract stiffness tensor from mech pde
			std::map<RegionIdType, BaseMaterial*>  mechMat;
			mechMat = pde1_->GetMaterialData();

			shared_ptr<CoefFunction > stiffCoef;
			stiffCoef = mechMat[regionId]->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, tensorType, 
              Global::REAL, true );

			// get hyst operator for current region from elec pde
			PtrCoefFct hystOperator;
			shared_ptr<CoefFunctionMulti> hystCoefs;
			hystCoefs = pde2_->GetHystCoefs();
			hystOperator =hystCoefs->GetRegionCoef(regionId);
				
      // IMPORTANT: pass information about stiffCoef and couplCoef to hystOperrator
      // there it will not only be used for the Following MatCoefFnc but also used
      // in the single pdes to get correct matrix entries!
      hystOperator->SetElastAndCouplTensor(stiffCoef, couplingCoef);
      
			// create hyst material
			// note: curCoef holds the coupling tensor
			// hyst case might be unsymmetric
			// > in case of deltaFormulation, elecToMech holds an addition deltaS/deltaE
			PtrCoefFct mechToElec = hystOperator->GenerateMatCoefFnc("CouplingMechToElec");
			PtrCoefFct elecToMech = hystOperator->GenerateMatCoefFnc("CouplingElecToMech");
			
      //			std::cout << "check size of matcoeffnc" << std::endl;
//			UInt numRows, numCols;
//			mechToElec->GetTensorSize(numRows,numCols);
      //			std::cout << "mechToElec: " << numRows << " " << numCols << std::endl;
			
//			elecToMech->GetTensorSize(numRows,numCols);
      //			std::cout << "elecToMech: " << numRows << " " << numCols << std::endl;
			
			if ( isComplex ) {
				if( subType_ == "axi" ) {
					*elecToMechInt = new ADBInt<Complex>(new StrainOperatorAxi<FeH1,Complex>(),
                  new GradientOperator<FeH1,2,1,Complex>(),
                  elecToMech, 1.0, true );
					*mechToElecInt = new ADBInt<Complex>(new GradientOperator<FeH1,2,1,Complex>(),
                  new StrainOperatorAxi<FeH1,Complex>(),
                  mechToElec, 1.0, true );
				} else if( subType_ == "planeStrain" ) {
					*elecToMechInt = new ADBInt<Complex>(new StrainOperator2D<FeH1,Complex>(),
                  new GradientOperator<FeH1,2,1,Complex>(),
                  elecToMech, 1.0, true);
					*mechToElecInt = new ADBInt<Complex>(new GradientOperator<FeH1,2,1,Complex>(),
                  new StrainOperator2D<FeH1,Complex>(),
                  mechToElec, 1.0, true);
				} else if( subType_ == "planeStress" ) {
					*elecToMechInt = new ADBInt<Complex>(new StrainOperator2D<FeH1,Complex>(),
                  new GradientOperator<FeH1,2,1,Complex>(),
                  elecToMech, 1.0, true);
					*mechToElecInt = new ADBInt<Complex>(new GradientOperator<FeH1,2,1,Complex>(),
                  new StrainOperator2D<FeH1,Complex>(),
                  mechToElec, 1.0, true);
				} else if( subType_ == "3d") {
					*elecToMechInt = new ADBInt<Complex>(new StrainOperator3D<FeH1,Complex>(),
                  new GradientOperator<FeH1,3,1,Complex>(),
                  elecToMech, 1.0, true);
					*mechToElecInt = new ADBInt<Complex>(new GradientOperator<FeH1,3,1,Complex>(),
                  new StrainOperator3D<FeH1,Complex>(),
                  mechToElec, 1.0, true);
				} else if( subType_ == "2.5d") {
					*elecToMechInt = new ADBInt<Complex>(new StrainOperator2p5D<FeH1, Complex>(),
                  new GradientOperator2p5D<FeH1, 2, 1, Complex>(),
                  elecToMech, 1.0, true);
					*mechToElecInt = new ADBInt<Complex>(new GradientOperator2p5D<FeH1, 2, 1, Complex>(),
                  new StrainOperator2p5D<FeH1, Complex>(),
                  mechToElec, 1.0, true);
				} else {
					EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
				}
			}
			else {
				if( subType_ == "axi" ) {
					*elecToMechInt = new ADBInt<Double>(new StrainOperatorAxi<FeH1>(),
                  new GradientOperator<FeH1,2>(),
                  elecToMech, 1.0, true );
					*mechToElecInt = new ADBInt<Double>(new GradientOperator<FeH1,2>(),
                  new StrainOperatorAxi<FeH1>(),
                  mechToElec, 1.0, true );
				} else if( subType_ == "planeStrain" ) {
					*elecToMechInt = new ADBInt<Double>(new StrainOperator2D<FeH1>(),
                  new GradientOperator<FeH1,2>(),
                  elecToMech, 1.0, true);
					*mechToElecInt = new ADBInt<Double>(new GradientOperator<FeH1,2>(),
                  new StrainOperator2D<FeH1>(),
                  mechToElec, 1.0, true);
				} else if( subType_ == "planeStress" ) {
					*elecToMechInt = new ADBInt<Double>(new StrainOperator2D<FeH1>(),
                  new GradientOperator<FeH1,2>(),
                  elecToMech, 1.0, true);
					*mechToElecInt = new ADBInt<Double>(new GradientOperator<FeH1,2>(), 
                  new StrainOperator2D<FeH1>(),
                  mechToElec, 1.0, true);
				} else if( subType_ == "3d") {
					*elecToMechInt = new ADBInt<Double>(new StrainOperator3D<FeH1>(),
                  new GradientOperator<FeH1,3>(),
                  elecToMech, 1.0, true);
          *mechToElecInt = new ADBInt<Double>(new GradientOperator<FeH1,3>(),
                  new StrainOperator3D<FeH1>(),
                  mechToElec, 1.0, true);
				} else if( subType_ == "2.5d") {
					*elecToMechInt = new ADBInt<Double>(new StrainOperator2p5D<FeH1, Double>(),
                  new GradientOperator2p5D<FeH1, 2, 1, Double>(),
                  elecToMech, 1.0, true);
					*mechToElecInt = new ADBInt<Double>(new GradientOperator2p5D<FeH1, 2, 1, Double>(),
                  new StrainOperator2p5D<FeH1, Double>(),
                  mechToElec, 1.0, true);
				} else {
					EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
				}
			}

      // 5.6.2018: create output coef function at this point, too > compare to elecPDE
      //      PtrCoefFct irrStrains = hystOperator->GenerateOutputCoefFnc("IrrStrainsPiezo");
      //      irrStrains_->AddRegion( regionId, irrStrains);
      //      PtrCoefFct irrStresses = hystOperator->GenerateOutputCoefFnc("IrrStressesPiezo");
      //      irrStresses_->AddRegion( regionId, irrStresses);
      PtrCoefFct irrStrains = hystOperator->GenerateOutputCoefFnc("IrrStrainsPiezo_TensorForm");
      irrStrains_->AddRegion( regionId, irrStrains);
      PtrCoefFct irrStressesVec = hystOperator->GenerateOutputCoefFnc("IrrStressesPiezo_VectorForm");
      irrStressesVector_->AddRegion( regionId, irrStressesVec);
      PtrCoefFct irrStresses = hystOperator->GenerateOutputCoefFnc("IrrStressesPiezo_TensorForm");
      irrStresses_->AddRegion( regionId, irrStresses);
      PtrCoefFct coupledIrrStrains = hystOperator->GenerateOutputCoefFnc("CoupledIrrStrainsPiezo");
      coupledIrrStrains_->AddRegion( regionId, coupledIrrStrains);

			return true;
		} else {
			// NON-HYST CASE 
			*elecToMechInt = GetStiffIntegrator( actSDMat, regionId, isComplex );
			return false;
		}
  }
  
  BaseBDBInt *
  PiezoCoupling::GetStiffIntegrator( BaseMaterial* actSDMat,
          RegionIdType regionId,
          bool isComplex ) {
    
    // Get region name
    std::string regionName = ptGrid_->GetRegion().ToString( regionId );
       SubTensorType tensorType = NO_TENSOR;
    if( subType_ == "planeStrain" || subType_ == "planeStress" ) {
      tensorType = PLANE_STRAIN;
    } else if( subType_ == "axi" ){
      tensorType = AXI;
    } else if (subType_ == "3d" || subType_ == "2.5d"){
      tensorType = FULL;
    } else {
      EXCEPTION( "Unknown subtype '" << subType_ << "'" );
    }
    // ------------------------
    //  Obtain linear material
    // ------------------------
    shared_ptr<CoefFunction > curCoef;
    if( isComplex ) {
      curCoef = actSDMat->GetTensorCoefFnc(PIEZO_TENSOR, tensorType, 
              Global::COMPLEX, true  );
    } else {
      curCoef = actSDMat->GetTensorCoefFnc(PIEZO_TENSOR, tensorType, 
              Global::REAL, true );
    }	
		BaseBDBInt * integ = NULL;
		
		// NON-HYST CASE -> Symmetric
		// ----------------------------------------
		//  Determine correct stiffness integrator 
		// ----------------------------------------

		// NOTE: here we have to couple +Bu with +GradV as in the constitutive equations Bu and -E are coupled!
		// -> no factor -1 here
		if ( isComplex ) {
			if( subType_ == "axi" ) {
				integ = new ADBInt<Complex>(new StrainOperatorAxi<FeH1,Complex>(),
                new GradientOperator<FeH1,2,1,Complex>(),
                curCoef, 1.0, true );
			} else if( subType_ == "planeStrain" ) {
				integ = new ADBInt<Complex>(new StrainOperator2D<FeH1,Complex>(),
                new GradientOperator<FeH1,2,1,Complex>(),
                curCoef, 1.0, true);
			} else if( subType_ == "planeStress" ) {
				integ = new ADBInt<Complex>(new StrainOperator2D<FeH1,Complex>(),
                new GradientOperator<FeH1,2,1,Complex>(),
                curCoef, 1.0, true);
			} else if( subType_ == "3d") {
				integ = new ADBInt<Complex>(new StrainOperator3D<FeH1,Complex>(),
                new GradientOperator<FeH1,3,1,Complex>(),
                curCoef, 1.0, true);
			} else if( subType_ == "2.5d") {
				integ = new ADBInt<Complex>(new StrainOperator2p5D<FeH1, Complex>(),
                new GradientOperator2p5D<FeH1, 2, 1, Complex>(),
                curCoef, 1.0, true);
			} else {
				EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
			}
		}
		else {
			if( subType_ == "axi" ) {
				integ = new ADBInt<Double>(new StrainOperatorAxi<FeH1>(),
                new GradientOperator<FeH1,2>(),
                curCoef, 1.0, true );
			} else if( subType_ == "planeStrain" ) {
				integ = new ADBInt<Double>(new StrainOperator2D<FeH1>(),
                new GradientOperator<FeH1,2>(),
                curCoef, 1.0, true);
			} else if( subType_ == "planeStress" ) {
				integ = new ADBInt<Double>(new StrainOperator2D<FeH1>(),
                new GradientOperator<FeH1,2>(),
                curCoef, 1.0, true);
			} else if( subType_ == "3d") {
				integ = new ADBInt<Double>(new StrainOperator3D<FeH1>(),
                new GradientOperator<FeH1,3>(),
                curCoef, 1.0, true);
			} else if( subType_ == "2.5d") {
				integ = new ADBInt<Double>(new StrainOperator2p5D<FeH1, Double>(),
                new GradientOperator2p5D<FeH1, 2, 1, Double>(),
                curCoef, 1.0, true);
			} else {
				EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
			}
		}

    return integ;

  }
  BaseBDBInt* PiezoCoupling::GetStiffIntegrator(BaseMaterial* actSDMat, RegionIdType regionId, bool isComplex, PtrCoefFct scalingFactor)
  {
    SubTensorType tensorType = NO_TENSOR;
    if (subType_ == "planeStrain" || subType_ == "planeStress")
      tensorType = PLANE_STRAIN;
    else if (subType_ == "axi")
      tensorType = AXI;
    else if (subType_ == "3d" || subType_ == "2.5d")
      tensorType = FULL;
    else
      EXCEPTION( "Unknown subtype '" << subType_ << "'" );
    // ------------------------
    //  Obtain linear material
    // ------------------------
    MathParser* mp = pde1_->GetDomain()->GetMathParser();
    shared_ptr<CoefFunction > curCoef, curCoefScl;
    curCoef = actSDMat->GetTensorCoefFnc(PIEZO_TENSOR, tensorType, Global::COMPLEX, true);
    curCoefScl = CoefFunction::Generate(mp, Global::COMPLEX, CoefXprTensScalOp(mp, curCoef, scalingFactor, CoefXpr::OP_MULT));

    // ----------------------------------------
    //  Determine correct stiffness integrator
    // ----------------------------------------
       // NOTE: here we have to couple +Bu with +GradV as in the constitutive equations Bu and -E are coupled!
    // -> no factor -1 here
    BaseBDBInt * integ = NULL;
    if (subType_ == "planeStrain" || subType_ == "planeStress")
      integ = new ADBInt<Complex>(new ScaledStrainOperator2D<FeH1, Complex>(), new ScaledGradientOperator<FeH1, 2, Complex>(),
              curCoefScl, 1.0, true);
    else if (subType_ == "3d")
      integ = new ADBInt<Complex>(new ScaledStrainOperator3D<FeH1, Complex>(), new ScaledGradientOperator<FeH1, 3, Complex>(),
              curCoefScl, 1.0, true);
    else if (subType_ == "2.5d")
      integ = new ADBInt<Complex>(new ScaledStrainOperator2p5D<FeH1, Complex>(), new ScaledGradientOperator2p5D<FeH1, 2, Complex>(),
              curCoefScl, 1.0, true);
    else
      EXCEPTION("Subtype '" << subType_ << "' unknown for mechanic physic");

    return integ;
  }

  template<typename DATA_TYPE>
  BiLinearForm* PiezoCoupling::GetNormalPiezoFluxIntegrator(PtrCoefFct scalCoefFucn, PtrCoefFct coefFuncPMLVec,
          Double factor, BiLinearForm::CouplingDirection cplDir, bool fluxOpA)
  {
    BiLinearForm* integ = NULL;
    BaseBOperator *fluxOp = NULL, *idOp = NULL;

    if (dim_ == 3)
    {
      if (coefFuncPMLVec)
        fluxOp = new SurfaceNormalPiezoFluxOperator<FeH1, 3, 3, Complex>(subType_, coefFuncPMLVec);
      else
        fluxOp = new SurfaceNormalPiezoFluxOperator<FeH1, 3, 3, Double>(subType_);
      idOp = new SurfaceMultiIdOp<FeH1, 3, 1>();
    }
    else if (dim_ == 2 && subType_ == "2.5d")
    {
      if (coefFuncPMLVec)
        fluxOp = new SurfaceNormalPiezoFluxOperator<FeH1, 2, 3, Complex>(subType_, coefFuncPMLVec);
      else
        fluxOp = new SurfaceNormalPiezoFluxOperator<FeH1, 2, 3, Double>(subType_);
      idOp = new SurfaceMultiIdOp<FeH1, 3, 1>();
    }
    else
    {
      if (coefFuncPMLVec)
        fluxOp = new SurfaceNormalPiezoFluxOperator<FeH1, 2, 2, Complex>(subType_, coefFuncPMLVec);
      else
        fluxOp = new SurfaceNormalPiezoFluxOperator<FeH1, 2, 2, Double>(subType_);
      idOp = new SurfaceMultiIdOp<FeH1, 2, 1>();
    }
    if (fluxOpA)
      integ = new SurfaceNitscheABInt<DATA_TYPE, DATA_TYPE>(fluxOp, idOp, scalCoefFucn, factor, cplDir, true);
    else
      integ = new SurfaceNitscheABInt<DATA_TYPE, DATA_TYPE>(idOp, fluxOp, scalCoefFucn, factor, cplDir, true);

    return integ;
  }

  template<typename DATA_TYPE>
  BiLinearForm* PiezoCoupling::GetNormalPiezoStrainIntegrator(PtrCoefFct scalCoefFucn, PtrCoefFct coefFuncPMLVec,
          Double factor, BiLinearForm::CouplingDirection cplDir, bool fluxOpA)
  {
    BiLinearForm* integ = NULL;
    BaseBOperator *fluxOp = NULL, *idOp = NULL;
       if (dim_ == 3)
    {
      if (coefFuncPMLVec)
        fluxOp = new SurfaceNormalPiezoStrainOperator<FeH1, 3, 1, Complex>(subType_, coefFuncPMLVec);
      else
        fluxOp = new SurfaceNormalPiezoStrainOperator<FeH1, 3, 1, Double>(subType_);
      idOp = new SurfaceIdentityOperator<FeH1, 3, 3>();
    }
    else if (dim_ == 2 && subType_ == "2.5d")
    {
      if (coefFuncPMLVec)
        fluxOp = new SurfaceNormalPiezoStrainOperator<FeH1, 2, 1, Complex>(subType_, coefFuncPMLVec);
      else
        fluxOp = new SurfaceNormalPiezoStrainOperator<FeH1, 2, 1, Double>(subType_);
      idOp = new SurfaceIdentityOperator<FeH1, 2, 3>();
    }
    else
    {
      if (coefFuncPMLVec)
        fluxOp = new SurfaceNormalPiezoStrainOperator<FeH1, 2, 1, Complex>(subType_, coefFuncPMLVec);
      else
        fluxOp = new SurfaceNormalPiezoStrainOperator<FeH1, 2, 1, Double>(subType_);
      idOp = new SurfaceIdentityOperator<FeH1, 2, 2>();
    }

    if (fluxOpA)
      integ = new SurfaceNitscheABInt<DATA_TYPE, DATA_TYPE>(fluxOp, idOp, scalCoefFucn, factor, cplDir, true);
    else
      integ = new SurfaceNitscheABInt<DATA_TYPE, DATA_TYPE>(idOp, fluxOp, scalCoefFucn, factor, cplDir, true);

    return integ;
  }

  void PiezoCoupling::InitTimeStepping(){
    if ( analysisType_ == BasePDE::TRANSIENT ) {
      Double dt;
      dt = dynamic_cast<TransientDriver*>(domain_->GetSingleDriver())
              ->GetDeltaT();

      //in this case we additionally need to define
      //a timestepping for the elecPDE
      shared_ptr<BaseFeFunction> elecFct = pde2_->GetFeFunction(ELEC_POTENTIAL);

      TimeSchemeGLM::NonLinType nlType = (nonLin_ || isHyst_)? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;

      //get possible alpha parameter of HHT-time stepping from mechanical PDE
      PtrParamNode mechParam = pde1_->GetMyParam();
      Double alpha = mechParam->Get("timeStepAlpha")->As<Double>();

      GLMScheme * scheme1 = new Newmark(0.5,0.25,alpha); 
      shared_ptr<BaseTimeScheme> elecScheme(new TimeSchemeGLM(scheme1, 0, nlType) );

      elecFct->SetTimeScheme(elecScheme);
      elecFct->GetTimeScheme()->Init(elecFct->GetSingleVector(),dt);
    }
  }

  template BiLinearForm* PiezoCoupling::GetNormalPiezoFluxIntegrator<Double>(PtrCoefFct, PtrCoefFct, Double, BiLinearForm::CouplingDirection, bool);
  template BiLinearForm* PiezoCoupling::GetNormalPiezoStrainIntegrator<Double>(PtrCoefFct, PtrCoefFct, Double, BiLinearForm::CouplingDirection, bool);
  template BiLinearForm* PiezoCoupling::GetNormalPiezoFluxIntegrator<Complex>(PtrCoefFct, PtrCoefFct, Double, BiLinearForm::CouplingDirection, bool);
  template BiLinearForm* PiezoCoupling::GetNormalPiezoStrainIntegrator<Complex>(PtrCoefFct, PtrCoefFct, Double, BiLinearForm::CouplingDirection, bool);

}

