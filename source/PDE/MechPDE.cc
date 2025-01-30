// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "MechPDE.hh"

#include <sstream>
#include <iomanip>
#include <set>
#include <boost/tokenizer.hpp>

#include "../Forms/Operators/PiolaStressOperator.hh"
#include "General/defs.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Driver/Assemble.hh"

#include "DataInOut/SimState.hh"
#include "DataInOut/SimInOut/hdf5/SimInputHDF5.hh"
#include "DataInOut/SimInOut/hdf5/SimOutputHDF5.hh"

// include elements
#include "FeBasis/FeFunctions.hh"
#include "FeBasis/H1/H1Elems.hh"

// new integrator concept
#include "Forms/BiLinForms/BDBInt.hh"
#include "Forms/BiLinForms/BBInt.hh"
#include "Forms/BiLinForms/ABInt.hh"
#include "Forms/BiLinForms/ICModesInt.hh"
#include "Forms/LinForms/SingleEntryInt.hh"
#include "Forms/LinForms/BUInt.hh"
#include "Forms/LinForms/BDUInt.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/IdOpNormalStrain.hh"
#include "Forms/Operators/IdentityOperatorNormalTrans.hh"
#include "Forms/Operators/StrainOperator.hh"
#include "Forms/Operators/SurfaceNormalStressOperator.hh"
#include "Forms/Operators/SurfaceOperators.hh"
#include "Forms/BiLinForms/SingleEntryBiLinInt.hh"

// new postprocessing concept
#include "Domain/Results/ResultFunctor.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"

#include "Domain/CoefFunction/CoefXpr.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefFunctionCache.hh"
#include "Domain/CoefFunction/CoefFunctionConst.hh"
#include "Domain/CoefFunction/CoefFunctionPML.hh"
#include "Domain/CoefFunction/CoefFunctionMapping.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include "Domain/CoefFunction/CoefFunctionCompound.hh"
#include "Domain/CoefFunction/CoefFunctionSurf.hh"
#include "Domain/CoefFunction/CoefFunctionOpt.hh"
#include "Driver/SolveSteps/StdSolveStep.hh"
#include "Driver/TimeSchemes/TimeSchemeGLM.hh"
#include "Driver/EigenFrequencyDriver.hh"
#include "Driver/BucklingDriver.hh"

#include "Domain/Mesh/NcInterfaces/MortarInterface.hh"

#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Context.hh"
#include "Optimization/Excitation.hh"

#include "Domain/Results/ResultInfo.hh"

namespace CoupledField {

  DEFINE_LOG(mechpde, "mechpde")

          /** the static test strain enum mapping */
  Enum<MechPDE::TestStrain> MechPDE::testStrain;
  
  MechPDE::MechPDE(Grid * aptgrid, PtrParamNode paramNode,PtrParamNode infoNode,
          shared_ptr<SimState> simState, Domain* domain )
  :SinglePDE( aptgrid, paramNode, infoNode, simState, domain ) {
    pdename_          = "mechanic";
    pdematerialclass_ = MECHANIC;
    
    nonLin_        = false;
    nonLinMaterial_= false;
    
    //! Always use total Lagrangian formulation 
    updatedGeo_        = false;
    
    // ****************************
    // DETERMINE GEOMETRY
    // ****************************
    
    // Get problem geometry and PDE subtype
    myParam_->GetValue("subType", subType_ );
    
    std::string probGeo = domain_->GetParamRoot()->Get("domain")->Get("geometryType")->As<std::string>();
    
    // Set number of degrees of freedom and
    // ensure that subtype fits to problem geometry
    if ( subType_ == "3d" && probGeo == "3d" ) {
      stressDim_ = 6;
      tensorType_ = FULL;
      dofNames_ = "x", "y", "z";
    }
    else if (subType_ == "2.5d" && probGeo == "plane")
    {
      stressDim_ = 6;
      tensorType_ = FULL;
      dofNames_ = "x", "y", "z";
    }
    else if ( subType_ == "axi" && probGeo == "axi" ) {
      isaxi_ = true;
      stressDim_ = 4;
      tensorType_ = AXI;
      dofNames_ = "r", "z";
    }
    else if ( subType_ == "planeStrain" && probGeo == "plane" ) {
      stressDim_ = 3;
      tensorType_ = PLANE_STRAIN;
      dofNames_ = "x", "y";
    }
    
    else if ( subType_ == "planeStress" && probGeo == "plane" ) {
      stressDim_ = 3;
      tensorType_ = PLANE_STRESS;
      dofNames_ = "x", "y";
    }
    
    else if ( subType_ == "flatShell" ) {
      stressDim_ = 3;
      dofNames_ = "x", "y";
    }
    else
      EXCEPTION("Subtype '" <<  subType_ << "' of PDE '" <<  pdename_ <<  "' does not fit to problem  geometry '" << probGeo << "'");
    
    // Sanity check: 3D can only be computed if 3D elements are present/
    if(subType_ == "3d" && ptGrid_->GetNumElemOfDim(3) == 0)
      EXCEPTION("Can not calculate 3D mechanics without 3D elements in the grid!");
    
    // thermal stress coefFunction
    thermalStress_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, stressDim_, 1, isComplex_, true));
    
    // thermal strain coefFunction
    thermalStrain_.reset(new CoefFunctionMulti(CoefFunction::VECTOR, stressDim_, 1, isComplex_, true));
  }
  
  
  MechPDE::~MechPDE()
  {
    
  }

  void MechPDE::ReadDampingInformation()
  {
    // get regions of current PDE
    ParamNodeList regionParamNodes = myParam_->Get("regionList")->GetChildren();
    // corresponding region id
    RegionIdType actRegionId;
    // corresponding region name and damping id
    std::string actRegionName, actDampingId;
    // try to get the dampingList and return if it is not specified
    PtrParamNode dampListNode = myParam_->Get("dampingList", ParamNode::PASS);
    if (dampListNode) {
      // get the individual damping nodes
      ParamNodeList dampNodes = dampListNode->GetChildren();
      // map of damping ids from the xml and corresponding damping types
      std::map<std::string, DampingType> idDampType;

      // Run over all region param nodes and assign the required damping information
      for (UInt iRegion = 0; iRegion < regionParamNodes.GetSize(); ++iRegion) {
        regionParamNodes[iRegion]->GetValue("name", actRegionName);
        regionParamNodes[iRegion]->GetValue("dampingId", actDampingId);

        // pass if no damping is defined
        if (actDampingId == "")
          continue;

        // parse region id from region name
        actRegionId = ptGrid_->GetRegion().Parse(actRegionName);

        // now, read the damping information from the dampingList
        for (UInt iChild = 0; iChild < dampNodes.GetSize(); ++iChild) {
          std::string dampString = dampNodes[iChild]->GetName();
          std::string actId = dampNodes[iChild]->Get("id")->As<std::string>();
          // only consider damping information for the current damping id
          if (actId != actDampingId)
            continue;

          // determine type of damping
          DampingType actType;
          String2Enum(dampString, actType);

          // store damping type string
          idDampType[actId] = actType;
          // break after the information is set as only one damping ID per region is possible
          break;
        }

        // check actDampingId was indeed registerd above
        if (idDampType.count(actDampingId) == 0)
          EXCEPTION("Damping with id '" << actDampingId << "' of region '" << actRegionName << "' was not found. Is it defined in the 'dampingList'?");

        // assign damping type to the region
        dampingList_[actRegionId] = idDampType[actDampingId];

        // if Rayleigh damping is specified, parse and store the additional damping information
        if (dampingList_[actRegionId] == RAYLEIGH) {
          RaylDampingData actRaylCoeffs;
          materials_[actRegionId]->GetRayleighCoeffStrings(actRaylCoeffs.alpha, actRaylCoeffs.beta);
          regionRaylDamping_[actRegionId] = actRaylCoeffs;
          PtrParamNode in = infoNode_->Get(ParamNode::HEADER)->GetByVal("region", "name", ptGrid_->GetRegion().ToString(actRegionId));
          in->Get("alpha_M")->SetValue(actRaylCoeffs.alpha);
          in->Get("alpha_K")->SetValue(actRaylCoeffs.beta);
        }
        else if (dampingList_[actRegionId] == ADAPTED_LOSS_TANGENS_DELTA) {
          if (!(analysistype_ == BasePDE::HARMONIC))
            EXCEPTION("Adapted loss tangent delta damping is only allowed for harmonic analysis.");
          RaylDampingData actRaylCoeffs;
          materials_[actRegionId]->GetFreqAdaptedRayleighCoeffStrings(actRaylCoeffs.alpha, actRaylCoeffs.beta);
          regionRaylDamping_[actRegionId] = actRaylCoeffs;
          PtrParamNode in = infoNode_->Get(ParamNode::HEADER)->GetByVal("region", "name", ptGrid_->GetRegion().ToString(actRegionId));
          in->Get("alpha_M")->SetValue(actRaylCoeffs.alpha);
          in->Get("alpha_K")->SetValue(actRaylCoeffs.beta);
        }
        else if (dampingList_[actRegionId] == GLOBAL_RAYLEIGH) {
          EXCEPTION("Global Rayleigh damping is not yet implemented.");
          if (dampNodes.GetSize() != 1)
            EXCEPTION("Global Rayleigh damping does not allow for additional damping nodes defined.");
          RaylDampingData actRaylCoeffs;
          actRaylCoeffs.alpha = dampNodes[0]->Get("alpha")->As<std::string>();
          actRaylCoeffs.beta = dampNodes[0]->Get("beta")->As<std::string>();
          regionRaylDamping_[actRegionId] = actRaylCoeffs;
          PtrParamNode in = infoNode_->Get(ParamNode::HEADER)->GetByVal("region", "name", ptGrid_->GetRegion().ToString(actRegionId));
          in->Get("alpha_M")->SetValue(actRaylCoeffs.alpha);
          in->Get("alpha_K")->SetValue(actRaylCoeffs.beta);
        }
        // if Kelvin-Voigt damping is specified, read the viscous tensor and assign its coeffunction to the region
        else if (dampingList_[actRegionId] == KELVIN_VOIGT) {
          regionKelvinVoigtDamping_[actRegionId] = materials_[actRegionId]->GetTensorCoefFnc(MECH_KV_VISCOUS_TENSOR, tensorType_, Global::REAL);
        }
      }
    }
  }

  void MechPDE::ReadSoftening() {
    
    // Check if softeningList node is present
    std::string type, id;
    std::map<std::string, std::string> idSoftTypeMap;
    PtrParamNode softListNode = myParam_->Get("softeningList", ParamNode::PASS );
    PtrParamNode softInfo;
    if( softListNode ) {
      
      // Get child elements and read data
      ParamNodeList softNodes = softListNode->GetChildren();
      for( UInt i = 0; i < softNodes.GetSize(); i++ ) {
        type = softNodes[i]->GetName();
        softNodes[i]->GetValue( "id", id );
        idSoftTypeMap[id] = type;
      }
      
      if( softNodes.GetSize() ) {
        softInfo = infoNode_->Get("softeningList");
      }
    }
    
    // Now iterate over all regions and check for softening type
    ParamNodeList regionNodes =
            myParam_->Get("regionList")->GetList("region");
    for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
      
      // get region Name
      std::string regionName = regionNodes[i]->Get("name")->As<std::string>();
      RegionIdType regionId = ptGrid_->GetRegion().Parse( regionName );
      
      // get softeningId of current region
      id = "";
      regionNodes[i]->GetValue( "softeningId", id, ParamNode::PASS );
      if (id == "") continue;
      
      // try to find related softening definition
      if( idSoftTypeMap.find( id) == idSoftTypeMap.end() ) {
        EXCEPTION( "Softening with id '" << id << "' for region '"
                << regionName << "' could not be found in "
                << "softeningList " );
      }
      // assign correct softening type to current region and map to log
      regionSoftening_[regionId] = idSoftTypeMap[id];
      PtrParamNode regionNode = softInfo->Get("region",ParamNode::APPEND);
      regionNode->Get("name")->SetValue(regionName);
      regionNode->Get("type")->SetValue(idSoftTypeMap[id]);
    }
  }
  
  void MechPDE::InitNonLin()
  {
    SinglePDE::InitNonLin();
  }
  
  
  void MechPDE::DefineIntegrators() {
    
    // =======================================
    //  Get information about softening types
    // =======================================
    ReadSoftening();
    
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;
    
    // Get FESpace and FeFunction of mechanical displacement
    shared_ptr<BaseFeFunction> myFct = feFunctions_[MECH_DISPLACEMENT];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
    
    //flag indicating frequency PML formulation
    bool harmonicPML = false;
    // determine if we do bloch eigenfrequency analysis
    bool do_bloch = domain->GetDriver()->DoBlochModeEigenfrequency();
    
    unsigned int rows = 0;
    unsigned int cols = 0;
    materials_.begin()->second->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, tensorType_, Global::REAL)->GetTensorSize(rows, cols);
    assert(rows > 0 && cols > 0);
    
    //  Loop over all regions
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for(UInt iRegion = 0; iRegion < regions_.GetSize() ; iRegion++){
      actRegion = regions_[iRegion];
      actSDMat    = materials_[actRegion];
      
      // in case we are not 3D, we add the plane stress/strain tensor to .info.xml
      // in BaseMaterial::ToInfo() only the full tensor is stored 
      if(tensorType_ != FULL)
      {
        // append new when a region has another material
        PtrParamNode in = infoNode_->GetByVal("material", "name", actSDMat->GetName()); 
        PtrCoefFct subTensor = actSDMat->GetSubTensorCoefFnc(MECH_STIFFNESS_TENSOR, tensorType_, Global::COMPLEX);
        actSDMat->StoreTensor(in->Get("subtensor/tensor"), subTensor); // seems to create a string, hence we cannot add type to subsensor directly
        in->Get("subtensor/type")->SetValue(subType_);
      }
      
      // Get current region name
      std::string regionName = ptGrid_->GetRegion().ToString(actRegion);
      
      // complex material or bloch mode with complex B-matrices
      bool isComplex = do_bloch | complexMatData_[actRegion];
      
      // get list of nonlinearities
      StdVector<NonLinType> & nonLinTypes = regionNonLinTypes_[actRegion];
      
      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptGrid_ ) );
      actSDList->SetRegion( actRegion );
      
      // --------------------------
      //  Set region approximation
      // --------------------------
      // --- Set the approximation for the current region ---
      PtrParamNode curRegNode = myParam_->Get("regionList")->GetByVal("region","name",regionName.c_str());
      std::string polyId = curRegNode->Get("polyId")->As<std::string>();
      std::string integId = curRegNode->Get("integId")->As<std::string>();
      mySpace->SetRegionApproximation(actRegion, polyId,integId);
      
      // Take account of pml in frequency domain
      // 'coeffPMLScal' is the function, the material tensor is to be scaled by. if PML isn't defined, it's unity
      PtrCoefFct coefPMLScalPos, coefPMLScalNeg, coeffMass; // "ngeative" version for stiffness integrator
      PtrCoefFct coefPMLVec;
      PtrCoefFct posOne = CoefFunction::Generate(mp_, Global::REAL, "1.0");
      PtrCoefFct negOne = CoefFunction::Generate(mp_, Global::REAL, "-1.0");
      PtrParamNode pmlNode;
      std::string pmlFormul;
      
      // ====================================================================
      //  Standard Linear Stiffness
      // ====================================================================
      if( !nonLin_ ) {
        if (dampingList_[actRegion] == PML) {
          if (analysistype_ == HARMONIC) {
            harmonicPML = true;
            std::string dampId;
            curRegNode->GetValue("dampingId", dampId);
            pmlNode = myParam_->Get("dampingList")->GetByVal("pml", "id", dampId.c_str());
            pmlFormul = pmlNode->Get("formulation")->As<std::string>();
            
            // speed of sound is set to equal '-1.0' or '1.0' depending on stiffness or mass integrator
            if (pmlFormul == "classic")
            {
              coefPMLScalPos.reset(new CoefFunctionPML<Complex>(pmlNode, posOne,
                      ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, false));
              coefPMLScalNeg.reset(new CoefFunctionPML<Complex>(pmlNode, negOne,
                      ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, false));
              coefPMLVec.reset(new CoefFunctionPML<Complex>(pmlNode, posOne,
                      ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, true));
            }
            else if (pmlFormul == "shifted")
            {
              coefPMLScalPos.reset(new CoefFunctionShiftedPML<Complex>(pmlNode, posOne,
                      ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, false));
              coefPMLScalNeg.reset(new CoefFunctionShiftedPML<Complex>(pmlNode, negOne, // this used to be posOne, no idea what is correct, like this it works in the uniaxial case ... //ftoth
                      ptGrid_->GetEntityList(EntityList::ELEM_LIST, regionName), regions_, false));
              coefPMLVec.reset(new CoefFunctionShiftedPML<Complex>(pmlNode, posOne,
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
        
        // ====================================================================
        // Take account for mapping
        // ====================================================================
        PtrCoefFct factor = CoefFunction::Generate( mp_, Global::REAL, "1.0");
        shared_ptr<CoefFunction> coeffMAPScal, coeffMAPVec;
        bool isMapping = false;
        if( dampingList_[actRegion] == MAPPING ) {
          std::string dampId;
          curRegNode->GetValue("dampingId",dampId);
          if(analysistype_ == HARMONIC){
            EXCEPTION("Harmonic analysis not allowed!");
          }else{
            PtrParamNode mapNode = myParam_->Get("dampingList")->GetByVal("mapping","id",dampId.c_str());
            coeffMAPVec.reset(new CoefFunctionMapping<Double>(mapNode,factor,actSDList,regions_,true));
            coeffMAPScal.reset(new CoefFunctionMapping<Double>(mapNode,factor,actSDList,regions_,false));
            isMapping = true;
          }
        }
        
        BaseBDBInt* stiffInt;
        
        // We use a stiffness integrator that implements the scaled strain operator if a PML or Mapping domain is considered.
        // Otherwise, we proceed with the normal stiffness integrator.
        if(isMapping){
          stiffInt =  GetStiffIntegrator(actSDMat, actRegion, false, coeffMAPScal);
          stiffInt->SetBCoefFunctionOpA(coeffMAPVec);
        }
        else if (harmonicPML)
        {
          stiffInt =  GetStiffIntegrator(actSDMat, actRegion, harmonicPML, coefPMLScalNeg);
          stiffInt->SetBCoefFunctionOpA(coefPMLVec);
        }
        else
        {
          // in the optimization case the coef function will be CoefFunctionOpt
          stiffInt =  GetStiffIntegrator(actSDMat, actRegion, isComplex);
        }
        
        stiffInt->SetName("LinElastInt");
        stiffInt->SetFeSpace( mySpace);
        
        BiLinFormContext * stiffIntDescr = new BiLinFormContext(stiffInt, STIFFNESS );
        stiffIntDescr->SetEntities( actSDList, actSDList );
        stiffIntDescr->SetFeFunctions( myFct, myFct );
        assemble_->AddBiLinearForm(stiffIntDescr);
        bdbInts_.insert(std::pair<RegionIdType, BaseBDBInt *>(actRegion, stiffInt));
        LOG_DBG(mechpde) << "Add Lin BDB" << std::endl;

        // check for Rayleigh damping (stiffness part)
        if (dampingList_[actRegion] == RAYLEIGH || dampingList_[actRegion] == ADAPTED_LOSS_TANGENS_DELTA || dampingList_[actRegion] == GLOBAL_RAYLEIGH) {
          RaylDampingData &actDamp = (regionRaylDamping_[actRegion]);
          stiffIntDescr->SetSecDestMat(DAMPING, actDamp.beta);
        }
        // check for Kelvin-Voigt damping
        else if (dampingList_[actRegion] == KELVIN_VOIGT) {
          if (harmonicPML) {
            EXCEPTION("Kelvin-Voigt damping and harmonic PML is not implemented/tested.");
          }
          if (isMapping) {
            EXCEPTION("Kelvin-Voigt damping and mapping is not implemented/tested.");
          }
          // add damping integrators the same way as stiffness integrators are defined
          BaseBDBInt *dampInt;
          dampInt = GetKelvinVoigtDampingIntegrator(actSDMat, actRegion, isComplex);
          dampInt->SetName("LinDampInt");
          dampInt->SetFeSpace(mySpace);

          BiLinFormContext *dampIntDescr = new BiLinFormContext(dampInt, DAMPING);
          dampIntDescr->SetEntities(actSDList, actSDList);
          dampIntDescr->SetFeFunctions(myFct, myFct);
          assemble_->AddBiLinearForm(dampIntDescr);
          LOG_TRACE(mechpde) << "Add Kelvin Voigt Damp BDB" << std::endl;
        }
      }

      // ====================================================================
      //  Geometric Nonlinear Stiffness
      // ====================================================================
      else if ( nonLinTypes.Find(GEOMETRIC) != -1 ) {
    	  //nonlinear (overall) stiffness matrix
    	  BaseBDBInt *nlBInt = NULL;
    	  PtrCoefFct stiffCoeff;
    	  if( isComplex ) {
    		  stiffCoeff = actSDMat->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, tensorType_, Global::COMPLEX);
    	  }
    	  else {
    		  stiffCoeff = actSDMat->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, tensorType_, Global::REAL);
    	  }
    	  regionStiffness_[actRegion] = stiffCoeff;
        
    	  if (subType_ == "axi") {
    		  nlBInt = new BDBInt<Double>(new NonLinStrainOperatorAxi<FeH1, Double>(myFct), stiffCoeff, 1.0, false);
    	  }
    	  else if (subType_ == "planeStrain" || subType_ == "planeStress") {
     	    LOG_TRACE(mechpde) << "Add NonLinStrainOperator2D" << std::endl;
    		  nlBInt = new BDBInt<Double>(new NonLinStrainOperator2D<FeH1, Double>(myFct), stiffCoeff, 1.0, false);
    	  }
    	  else if (subType_ =="3d") {
    		  nlBInt = new BDBInt<Double>(new NonLinStrainOperator3D<FeH1, Double>(myFct), stiffCoeff, 1.0, false);
    	  }
    	  else {
    		  assert(false);
    	  }
        
    	  nlBInt->SetSolDependent(true);
        nlBInt->SetFeSpace(mySpace);
        nlBInt->SetName("NonLinearStrainInt");
        
        BiLinFormContext *nlContext = new BiLinFormContext(nlBInt, STIFFNESS);
        nlContext->SetEntities(actSDList, actSDList);
        nlContext->SetFeFunctions(myFct, myFct);
        assemble_->AddBiLinearForm(nlContext);

        // check for damping
        if (dampingList_[actRegion] == RAYLEIGH || dampingList_[actRegion] == ADAPTED_LOSS_TANGENS_DELTA || dampingList_[actRegion] == GLOBAL_RAYLEIGH) {
          RaylDampingData &actDamp = (regionRaylDamping_[actRegion]);
          nlContext->SetSecDestMat(DAMPING, actDamp.beta);
        }
        else if (dampingList_[actRegion] == KELVIN_VOIGT) {
          EXCEPTION("Kelvin-Voigt damping and non-linear elasticity is not implemented/tested.");
        }

        // Important: Add bdb-integrator to global list, as we need them later
        // for calculation of postprocessing results.
        bdbInts_.insert( std::pair<RegionIdType, BaseBDBInt*>(actRegion,nlBInt) );
        
        if( nonLinMethod_ == NEWTON ) {
          //here we define the tangent matrix
          BaseBDBInt *piolaInt = NULL;
          
          PtrCoefFct piolaTensor(new CoefFunction2ndPiolaTensor(tensorType_,
                  stiffCoeff, myFct));
          
          if (subType_ == "axi") {
            piolaInt = new BDBInt<Double>(new PiolaStressOperatorAxi<FeH1>(false),
                    piolaTensor, 1.0, false);
          }
          else if (subType_ == "planeStrain" || subType_ == "planeStress") {
            piolaInt = new BDBInt<Double>(new PiolaStressOperator<FeH1, 2>(false),
                    piolaTensor, 1.0, false);
          }
          else if (subType_ =="3d") {
            piolaInt = new BDBInt<Double>(new PiolaStressOperator<FeH1, 3>(false),
                    piolaTensor, 1.0, false);
          }
          else {
            assert(false);
          }
          
          //very important: ist the tangent stiffness matrix
          piolaInt->SetNewtonBiLinearForm();
          
          piolaInt->SetSolDependent(true);
          piolaInt->SetFeSpace(mySpace);
          piolaInt->SetName("NonLinearPiolaInt");
          
          BiLinFormContext *piolaContext = new BiLinFormContext(piolaInt, STIFFNESS);
          piolaContext->SetEntities(actSDList, actSDList);
          piolaContext->SetFeFunctions(myFct, myFct);
          assemble_->AddBiLinearForm(piolaContext);
        }
      }
      
      
      //prestressing
      PtrParamNode preStressNode;
      PtrParamNode bcNode = this->myParam_->Get("bcsAndLoads");
      if(bcNode){
        preStressNode = bcNode->GetByVal("preStress","region",regionName.c_str(),ParamNode::PASS);
      }
      
      if( preStressNode ){
        ParamNodeList pbcList = bcNode->GetList("blochPeriodic");
        // complex prestressing can be due to: complex material; bloch mode with complex B-matrices; PML; periodic BC
        bool complexPre = do_bloch || harmonicPML || (pbcList.GetSize() > 0);
        
        PtrCoefFct preStressFct = CreatePreStressFct(complexPre, preStressNode);
        regionPreStress_[actRegion] = preStressFct;
        
        if(do_bloch || subType_ == "axi"){
          EXCEPTION("Prestressing is not available for Block periodic or axi-symmetric computations");
        }
        
        BaseBDBInt *preStressInt = NULL;
        if (harmonicPML)
        {
          preStressInt = GetPreStressIntegrator(preStressFct, actRegion, harmonicPML, coefPMLScalNeg);
          preStressInt->SetBCoefFunctionOpA(coefPMLVec);
        }
        else
          preStressInt = GetPreStressIntegrator(preStressFct, actRegion, complexPre);
        
        preStressInt->SetName("PreStressInt");
        preStressInt->SetFeSpace( mySpace );
        
        BiLinFormContext *preStressContext =  new BiLinFormContext( preStressInt, STIFFNESS );
        preStressContext->SetEntities( actSDList, actSDList );
        preStressContext->SetFeFunctions( myFct, myFct );
        
        assemble_->AddBiLinearForm( preStressContext );

        if (dampingList_[actRegion] == KELVIN_VOIGT) {
          EXCEPTION("Kelvin-Voigt damping and prestress is not implemented/tested.");
        }
      }
      
      // ====================================================================
      //  Standard Mass Integrator
      // ====================================================================
      PtrCoefFct densCoeff = actSDMat->GetScalCoefFnc( DENSITY,Global::REAL );
      PtrCoefFct densCoeffScaled;
      
      // when we do optimization we wrap the original CoefFunction. Don't check for region to handle dim-1 pressure on dim elements
      if(domain->HasDesign())
      {
        CoefFunctionOpt* tmpFnc = new CoefFunctionOpt(domain->GetDesign(), densCoeff, DENSITY, this);
        densCoeff.reset(tmpFnc);
      }
      
      BaseBDBInt *massInt = NULL;
      
      if ( harmonicPML ) {
        // mass integrator in a PML region: density coefficient function is scaled by sx*sy*sz
        densCoeffScaled = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprBinOp(mp_, densCoeff, coefPMLScalPos, CoefXpr::OP_MULT));
        if (dim_ == 2 && subType_ != "2.5d")
          massInt = new BBInt<Complex, Complex>(new IdentityOperator<FeH1, 2, 2>(), densCoeffScaled, 1.0, updatedGeo_);
        else
          massInt = new BBInt<Complex, Complex>(new IdentityOperator<FeH1, 3, 3>(), densCoeffScaled, 1.0, updatedGeo_);
      }
      else {
        // complex mass integrator due to complex bloch stiffness matrix
        if ( do_bloch )  {
          if (dim_ == 2 && subType_ != "2.5d")
            massInt = new BBInt<Complex, Complex>(new IdentityOperator<FeH1, 2, 2>(), densCoeff, 1.0);
          else
            massInt = new BBInt<Complex, Complex>(new IdentityOperator<FeH1, 3, 3>(), densCoeff, 1.0);
        }
        else {
          if (dim_ == 2 && subType_ != "2.5d")
            massInt = new BBInt<>(new IdentityOperator<FeH1, 2, 2>(), densCoeff, 1.0);
          else
            massInt = new BBInt<>(new IdentityOperator<FeH1, 3, 3>(), densCoeff, 1.0);
        }
      }
      
      massInt->SetName("MassInt");
      massInt->SetFeSpace( mySpace );
      
      // the integrator has a coef function but for the optimization case the opt coef needs to know also the integrator
      if(domain->HasDesign())
        dynamic_pointer_cast<CoefFunctionOpt>(densCoeff)->SetForm(massInt);
      
      BiLinFormContext *massContext =  new BiLinFormContext( massInt, MASS );
      massContext->SetEntities( actSDList, actSDList );
      massContext->SetFeFunctions( myFct, myFct );

      // Check for Rayleigh damping (mass part)
      if (dampingList_[actRegion] == RAYLEIGH || dampingList_[actRegion] == ADAPTED_LOSS_TANGENS_DELTA || dampingList_[actRegion] == GLOBAL_RAYLEIGH) {
        RaylDampingData &actDamp = regionRaylDamping_[actRegion];
        massContext->SetSecDestMat(DAMPING, actDamp.alpha);
      }

      // Important: Add mass-integrator to global list, as we need them later
      // for calculation of postprocessing results
      massInts_[actRegion] = massInt;
      assemble_->AddBiLinearForm( massContext );
      
      
      // in the end, at the region to the feFunction      // to be implemented
      myFct->AddEntityList( actSDList );


      // ====================================================================
      //  Geometric / Stress Stiffness Matrix for Buckling Analysis
      // ====================================================================
      if (analysistype_ == BUCKLING) {
        bool refLoadFound = false;
        PtrParamNode referenceLoadNode;
        PtrParamNode loadNode = this->myParam_->Get("bcsAndLoads");

        // check if BC are set at all. If not set, the eigenSolver may not be able to solve problem.
        if(!loadNode->Has("fix")){
          std::cout << "\n" << "++ " << "WARNING no fixed BC found. This can lead to non positive definite stiffness matrix" << "\n";
        }

        PtrParamNode variant_referenceLoadNode;

        if (loadNode)
          referenceLoadNode = loadNode->GetByVal("referenceLoad", "region", regionName.c_str(), ParamNode::PASS);

        if (referenceLoadNode)
          variant_referenceLoadNode = referenceLoadNode->Get("referenceDisplacement", ParamNode::PASS);
        if (variant_referenceLoadNode) {
          // generate displacement dependent stiffness matrix

          // TODO currently only conventional buckling analysis is implemented
          // see nonLinearStressOperator for implementation
          // add displacement dependent matrix to GEOMETRIC_STIFFNESS Matrix
          // (!! set bilinear factor to -1 !!)
          EXCEPTION("referenceLoad with displacement is not yet implemented.")
          refLoadFound = true;
        }
        else {
          // generate stress dependent stiffness matrix
          PtrCoefFct preStressFct = NULL;

          if (domain->GetOptimization()) {
            preStressFct = CreatePreStressFct(false, NULL);
            refLoadFound = true;
          }
          else {
            // check if current region has referenceLoad
            if (referenceLoadNode)
              variant_referenceLoadNode = referenceLoadNode->Get("referenceStress", ParamNode::PASS);
            if (variant_referenceLoadNode) {
              preStressFct = CreatePreStressFct(false, referenceLoadNode);
              refLoadFound = true;
            }
          }
          if (preStressFct != NULL) {
            // factor = -1 to account for different sign in the general eigenvalue formulation [A+w(-B)]x = 0 => Ax = wBx
            BaseBDBInt* preStressInt = GetPreStressIntegrator(preStressFct, actRegion, false, -1.0);

            preStressInt->SetName("PreStressInt");
            preStressInt->SetFeSpace(mySpace);

            // add stress dependent matrix to GEOMETRIC_STIFFNESS matrix
            BiLinFormContext *stressContext = new BiLinFormContext(preStressInt, GEOMETRIC_STIFFNESS);
            stressContext->SetEntities(actSDList, actSDList);
            stressContext->SetFeFunctions(myFct, myFct);

            assemble_->AddBiLinearForm(stressContext);
          }
        }

        if (refLoadFound == false)
          EXCEPTION("No reference load defined.")
      }// end Geometric Stiffness Matrix

    }
    
    // ====================================================================
    //  Concentrated Mechanical Network Elements
    // ====================================================================
    DefineConcentratedElems();
  }

  
  void MechPDE::DefineSurfaceIntegrators( ){

    PtrParamNode bcNode = myParam_->Get( "bcsAndLoads", ParamNode::PASS );
    if( bcNode ) {
      //========================================================================================
      // ABC boundaries
      //========================================================================================
      ParamNodeList abcNodes = bcNode->GetList( "absorbingBCs" );
      for( UInt i = 0; i < abcNodes.GetSize(); i++ ) {
        std::string regionName = abcNodes[i]->Get("name")->As<std::string>();
        shared_ptr<EntityList> actSDList =  ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST,regionName );
        std::string volRegName = abcNodes[i]->Get("volumeRegion")->As<std::string>();
        
        RegionIdType aRegion = ptGrid_->GetRegion().Parse(volRegName);
        
        std::string factor = "1.0";
        std::string factor2 = "2.0";
        //c_long = sqrt( (youngModulus*(1-poissonRatio)) / (density*(1-poissonRatio-2*poissonRatio^2)) )
        PtrCoefFct dens = materials_[aRegion]->GetScalCoefFnc(DENSITY,Global::REAL);
        PtrCoefFct poisson = materials_[aRegion]->GetScalCoefFnc(MECH_POISSON,Global::REAL);
        PtrCoefFct young = materials_[aRegion]->GetScalCoefFnc(MECH_EMODULUS,Global::REAL);
        
        // create nominator first
        PtrCoefFct nominator1 = CoefFunction::Generate(mp_, Global::REAL,
                CoefXprBinOp(mp_, //nominator
                young,
                CoefXprBinOp(mp_,
                factor,
                poisson,
                CoefXpr::OP_SUB
                ),
                CoefXpr::OP_MULT
                ) );
        // now create denominator
        PtrCoefFct denomiator1 = CoefFunction::Generate(mp_, Global::REAL,
                CoefXprBinOp(mp_, //denominator)
                dens,
                CoefXprBinOp(mp_,
                factor, // we calc 1 - (poisson + 2* poisson^2))
                CoefXprBinOp(mp_,
                poisson,
                CoefXprBinOp(mp_,
                factor2,
                CoefXprBinOp(mp_,
                poisson,
                poisson,
                CoefXpr::OP_MULT
                ),
                CoefXpr::OP_MULT
                ),
                CoefXpr::OP_ADD
                ),
                CoefXpr::OP_SUB
                ),
                CoefXpr::OP_MULT
                ) );
        
        
        PtrCoefFct cl = CoefFunction::Generate(mp_, Global::REAL,
                CoefXprUnaryOp(mp_,  //sqrt
                CoefXprBinOp(mp_, //div 
                nominator1,
                denomiator1,
                CoefXpr::OP_DIV
                ),
                CoefXpr::OP_SQRT) );
        
        //c_trans = sqrt( (youngModulus) / (2*density*(1+poissonRatio)) )
        PtrCoefFct ct = CoefFunction::Generate(mp_,Global::REAL,
                CoefXprUnaryOp(mp_, //sqrt
                CoefXprBinOp(mp_,              
                young,
                CoefXprBinOp(mp_,
                factor2,            
                CoefXprBinOp(mp_,
                dens,
                CoefXprBinOp(mp_,
                factor,
                poisson,
                CoefXpr::OP_ADD
                ),
                CoefXpr::OP_MULT
                ),
                CoefXpr::OP_MULT
                ),
                CoefXpr::OP_DIV
                ),
                CoefXpr::OP_SQRT) );
        
        // scaling factor for normal component: density * cl
        // scaling factor for tangential component: density * ct
        // the density however is a scalar so we can simply pass it as scalar coefficient function
        // this reduces the normal and tangential damping to cl and ct
        // to apply these values to the correct component it is useful to apply the scaling direcly to the ScaledIdNormalStrainOperator.
        // this operator will be used to create the B and B_transposed matrix so we have to pass only the ROOT of the scaling factors
        // to apply correct damping in normal and tangential directon we need to scale the normal vector
        // this is done by using IdOpNormalStrainScaled
        // however the resulting operator will be used for B and B_transped inside SurfaceBBInt
        // therefore we have to scale the normal just with the root the scaling factors
        
        PtrCoefFct coeffDampNormal = CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_,cl,CoefXpr::OP_SQRT) );
        
        PtrCoefFct coeffDampTangential = CoefFunction::Generate(mp_, Global::REAL, CoefXprUnaryOp(mp_,ct,CoefXpr::OP_SQRT) );                 
        
        BiLinearForm * abcInt = NULL;
        std::set<RegionIdType> volRegions;
        volRegions.insert(aRegion);
        
        if( dim_ == 2 ) {
          abcInt = new SurfaceBBInt<>(new ScaledIdNormalStrainOperator<FeH1,2,Double>(coeffDampNormal,coeffDampTangential), dens, 1.0, volRegions, updatedGeo_ );
        } else {
          abcInt = new SurfaceBBInt<>(new ScaledIdNormalStrainOperator<FeH1,3,Double>(coeffDampNormal,coeffDampTangential), dens, 1.0, volRegions, updatedGeo_ );
        }
        
        abcInt->SetName("abcIntegrator");
        BiLinFormContext *abcContext = new BiLinFormContext(abcInt, DAMPING );
        
        abcContext->SetEntities( actSDList, actSDList );
        abcContext->SetFeFunctions( feFunctions_[MECH_DISPLACEMENT] , feFunctions_[MECH_DISPLACEMENT]);
        feFunctions_[MECH_DISPLACEMENT]->AddEntityList( actSDList );
        assemble_->AddBiLinearForm( abcContext );
      }

      //========================================================================================
      // Normal X : assumes a normal traction proportional to the normal Y
      // X/Y = Stiffness/Displacement or Damping/Velocity or Mass/Acceleration
      // Essentially the mechanic boundary traction arising in the weak form u'*t_n
      // is replaced by u'*k u_n = u'*(k u*n n) = k u'*n u*n
      // where k is the parameter read from the input file, and u represents the unknown Y
      //========================================================================================
      typedef std::pair<std::string, FEMatrixType> normalBCtype;
      std::vector<normalBCtype> normalBCs;
      normalBCs.push_back(std::make_pair("normalStiffness", STIFFNESS));
      normalBCs.push_back(std::make_pair("normalDamping", DAMPING));
      normalBCs.push_back(std::make_pair("normalMass", MASS));
      for( std::vector<normalBCtype>::iterator it = normalBCs.begin() ; it != normalBCs.end(); ++it ){
        std::string xmlName = it->first;
        FEMatrixType feMat = it->second;
        LOG_DBG(mechpde) << "Reading '" << xmlName << "' definition";
        StdVector<shared_ptr<EntityList> > ent;
        StdVector<PtrCoefFct > kCoef;
        StdVector<std::string> volumeRegions;
        ReadRhsExcitation( xmlName , feFunctions_[MECH_DISPLACEMENT]->GetResultInfo()->dofNames, ResultInfo::SCALAR, ent, kCoef, updatedGeo_, volumeRegions);
        for( UInt i = 0; i < ent.GetSize(); ++i ) {
          // get the volume region for defining the correct normal direction
          RegionIdType aRegion = ptGrid_->GetRegion().Parse(volumeRegions[i]);
          std::set<RegionIdType> volRegion;
          volRegion.insert(aRegion);
          // check type of entitylist
          if (ent[i]->GetType() == EntityList::NODE_LIST) {
            EXCEPTION( xmlName << " must be defined on (surface) elements")
          }
          if ( kCoef[i]->IsComplex() && !(isComplex_) ) {
            EXCEPTION( xmlName << " is defied as complex but PDE is not")
          }
          // setup the integrator for: u'*t_n = u'*k u_n = u'*(k u*n n) = k u'*n u*n
          BiLinearForm * tangInt = NULL;
          if( kCoef[i]->IsComplex() ) {
            if (dim_ == 2){
              tangInt = new SurfaceBBInt<Complex,Double>(new IdentityOperatorNormalTrans<FeH1,2,2>(), kCoef[i], Complex(1.0,0), volRegion, updatedGeo_ );
            } else {
              tangInt = new SurfaceBBInt<Complex,Double>(new IdentityOperatorNormalTrans<FeH1,3,3>(), kCoef[i], Complex(1.0,0), volRegion, updatedGeo_ );
            }
          } else {
            if (dim_ == 2){
              tangInt = new SurfaceBBInt<>(new IdentityOperatorNormalTrans<FeH1,2,2>(), kCoef[i], 1.0, volRegion, updatedGeo_ );
            } else {
              tangInt = new SurfaceBBInt<>(new IdentityOperatorNormalTrans<FeH1,3,3>(), kCoef[i], 1.0, volRegion, updatedGeo_ );
            }
          }
          tangInt->SetName(xmlName + "Integrator");
          BiLinFormContext *tangContext = new BiLinFormContext(tangInt, feMat );
          tangContext->SetEntities( ent[i], ent[i]);
          tangContext->SetFeFunctions( feFunctions_[MECH_DISPLACEMENT], feFunctions_[MECH_DISPLACEMENT]);
          feFunctions_[MECH_DISPLACEMENT]->AddEntityList( ent[i] );
          assemble_->AddBiLinearForm( tangContext );
        }
      }

      // Bloch-periodic boundary conditions
      this->ptGrid_->MapEdges();
      this->ptGrid_->MapFaces();
      
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
        PtrCoefFct one = CoefFunction::Generate(mp_, Global::REAL, "1.0", "0.0");
        
        ParamNodeList regionsList = blochNodesList[i]->GetList("region");
        std::string volMasterName, volSlaveName;
        for (UInt j = 0; j < regionsList.GetSize(); j++)
        {
          std::string ncRegionName = regionsList[j]->Get("name")->As<std::string>();
          shared_ptr<BaseNcInterface> ncIf = ptGrid_->GetNcInterface(ptGrid_->GetNcInterfaceId(ncRegionName));
          if (!ncIf)
          {
            EXCEPTION("No interface with the name '" << ncRegionName << "' found!");
          }
          shared_ptr<MortarInterface> mortarIf = dynamic_pointer_cast<MortarInterface>(ncIf);
          assert(mortarIf);
          
          PtrCoefFct matDataTensorMas, matDataTensorSla, matData;
          RegionIdType volMasterId = mortarIf->GetPrimaryVolRegion();
          RegionIdType volSlaveId = mortarIf->GetSecondaryVolRegion();
          
          matDataTensorMas = regionStiffness_[volMasterId];
          matDataTensorSla = regionStiffness_[volSlaveId];
          assert(matDataTensorMas);
          assert(matDataTensorSla);
          
          if (formulation == "Nitsche")
          {
            std::string nitFac = blochNodesList[i]->Get("nitscheFactor")->As<std::string>();
            Double nitscheFactor = std::stod(nitFac);
            // master & slave penalty integrals
            BiLinearForm *pnlt_uM_vM = NULL;
            BiLinearForm *pnlt_uM_vS = NULL;
            BiLinearForm *pnlt_uS_vM = NULL;
            BiLinearForm *pnlt_uS_vS = NULL;
            // master & slave integrals with normal derivatives
            BiLinearForm *flux_DuM_vM = NULL;
            BiLinearForm *flux_DuM_vS = NULL;
            BiLinearForm *flux_uM_DvM = NULL;
            BiLinearForm *flux_uS_DvM = NULL;
            
            shared_ptr<ElemList> actSDList = ncIf->GetElemList();
            Double beta;
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
            
            //check for softening
            bool icModes = false;
            if (regionSoftening_[volMasterId] == "icModesTW" || regionSoftening_[volSlaveId]  == "icModesTW")
              icModes = true;
            
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
            // penalty integrators
            pnlt_uM_vM = GetPenaltyIntegrator<Complex>(factor, beta, BiLinearForm::PRIM_PRIM);
            pnlt_uM_vS = GetPenaltyIntegrator<Complex>(factorSqr, -beta, BiLinearForm::PRIM_SEC);
            pnlt_uS_vM = GetPenaltyIntegrator<Double>(one, -beta, BiLinearForm::SEC_PRIM);
            pnlt_uS_vS = GetPenaltyIntegrator<Complex>(factor, beta, BiLinearForm::SEC_SEC);
            // flux integrators
            if (matData->IsComplex())
            {
              flux_DuM_vM = GetFluxIntegrator<Complex>(one, coefFuncPMLVec, -1.0, BiLinearForm::PRIM_PRIM, true, icModes);
              flux_uM_DvM = GetFluxIntegrator<Complex>(factor, coefFuncPMLVec, -1.0, BiLinearForm::PRIM_PRIM, false, icModes);
              flux_DuM_vS = GetFluxIntegrator<Complex>(factor, coefFuncPMLVec, 1.0, BiLinearForm::PRIM_SEC, true, icModes);
              flux_uS_DvM = GetFluxIntegrator<Complex>(one, coefFuncPMLVec, 1.0, BiLinearForm::SEC_PRIM, false, icModes);
            }
            else
            {
              flux_DuM_vM = GetFluxIntegrator<Double>(one, coefFuncPMLVec, -1.0, BiLinearForm::PRIM_PRIM, true, icModes);
              flux_uM_DvM = GetFluxIntegrator<Complex>(factor, coefFuncPMLVec, -1.0, BiLinearForm::PRIM_PRIM, false, icModes);
              flux_DuM_vS = GetFluxIntegrator<Complex>(factor, coefFuncPMLVec, 1.0, BiLinearForm::PRIM_SEC, true, icModes);
              flux_uS_DvM = GetFluxIntegrator<Double>(one, coefFuncPMLVec, 1.0, BiLinearForm::SEC_PRIM, false, icModes);
            }
            
            // pass material data to the flux operators
            flux_DuM_vM->SetBCoefFunctionOpA(matData);
            flux_uM_DvM->SetBCoefFunctionOpB(matData);
            flux_DuM_vS->SetBCoefFunctionOpA(matData);
            flux_uS_DvM->SetBCoefFunctionOpB(matData);
            
            // master-master
            pnlt_uM_vM->SetName("pnlt_uM_vM");
            flux_DuM_vM->SetName("flux_DuM_vM");
            flux_uM_DvM->SetName("flux_uM_DvM");
            //master-slave
            pnlt_uM_vS->SetName("pnlt_uM_vS");
            flux_DuM_vS->SetName("flux_DuM_vS");
            // slave-master
            pnlt_uS_vM->SetName("pnlt_uS_vM");
            flux_uS_DvM->SetName("flux_uS_DvM");
            //slave-slave
            pnlt_uS_vS->SetName("pnlt_uS_vS");
            
            // BiLinearForm::PRIM_PRIM;
            SurfaceBiLinFormContext *pnlt_uM_vM_cont = new SurfaceBiLinFormContext(pnlt_uM_vM, STIFFNESS, BiLinearForm::PRIM_PRIM);
            SurfaceBiLinFormContext *flux_DuM_vM_cont = new SurfaceBiLinFormContext(flux_DuM_vM, STIFFNESS, BiLinearForm::PRIM_PRIM);
            SurfaceBiLinFormContext *flux_uM_DvM_cont = new SurfaceBiLinFormContext(flux_uM_DvM, STIFFNESS, BiLinearForm::PRIM_PRIM);
            // BiLinearForm::PRIM_SEC;
            SurfaceBiLinFormContext *pnlt_uM_vS_cont = new SurfaceBiLinFormContext(pnlt_uM_vS, STIFFNESS, BiLinearForm::PRIM_SEC);
            SurfaceBiLinFormContext *flux_DuM_vS_cont = new SurfaceBiLinFormContext(flux_DuM_vS, STIFFNESS, BiLinearForm::PRIM_SEC);
            // BiLinearForm::SEC_PRIM;
            SurfaceBiLinFormContext *pnlt_uS_vM_cont = new SurfaceBiLinFormContext(pnlt_uS_vM, STIFFNESS, BiLinearForm::SEC_PRIM);
            SurfaceBiLinFormContext *flux_uS_DvM_cont = new SurfaceBiLinFormContext(flux_uS_DvM, STIFFNESS, BiLinearForm::SEC_PRIM);
            // BiLinearForm::SEC_SEC;
            SurfaceBiLinFormContext *pnlt_uS_vS_cont = new SurfaceBiLinFormContext(pnlt_uS_vS, STIFFNESS, BiLinearForm::SEC_SEC);
            
            pnlt_uM_vM_cont->SetEntities(actSDList, actSDList);
            flux_DuM_vM_cont->SetEntities(actSDList, actSDList);
            flux_uM_DvM_cont->SetEntities(actSDList, actSDList);
            pnlt_uM_vS_cont->SetEntities(actSDList, actSDList);
            flux_DuM_vS_cont->SetEntities(actSDList, actSDList);
            pnlt_uS_vM_cont->SetEntities(actSDList, actSDList);
            flux_uS_DvM_cont->SetEntities(actSDList, actSDList);
            pnlt_uS_vS_cont->SetEntities(actSDList, actSDList);
            
            pnlt_uM_vM_cont->SetFeFunctions(feFunctions_[MECH_DISPLACEMENT], feFunctions_[MECH_DISPLACEMENT]);
            flux_DuM_vM_cont->SetFeFunctions(feFunctions_[MECH_DISPLACEMENT], feFunctions_[MECH_DISPLACEMENT]);
            flux_uM_DvM_cont->SetFeFunctions(feFunctions_[MECH_DISPLACEMENT], feFunctions_[MECH_DISPLACEMENT]);
            pnlt_uM_vS_cont->SetFeFunctions(feFunctions_[MECH_DISPLACEMENT], feFunctions_[MECH_DISPLACEMENT]);
            flux_DuM_vS_cont->SetFeFunctions(feFunctions_[MECH_DISPLACEMENT], feFunctions_[MECH_DISPLACEMENT]);
            pnlt_uS_vM_cont->SetFeFunctions(feFunctions_[MECH_DISPLACEMENT], feFunctions_[MECH_DISPLACEMENT]);
            flux_uS_DvM_cont->SetFeFunctions(feFunctions_[MECH_DISPLACEMENT], feFunctions_[MECH_DISPLACEMENT]);
            pnlt_uS_vS_cont->SetFeFunctions(feFunctions_[MECH_DISPLACEMENT], feFunctions_[MECH_DISPLACEMENT]);
            
            assemble_->AddBiLinearForm(pnlt_uM_vM_cont);
            assemble_->AddBiLinearForm(flux_DuM_vM_cont);
            assemble_->AddBiLinearForm(flux_uM_DvM_cont);
            assemble_->AddBiLinearForm(pnlt_uM_vS_cont);
            assemble_->AddBiLinearForm(flux_DuM_vS_cont);
            assemble_->AddBiLinearForm(pnlt_uS_vM_cont);
            assemble_->AddBiLinearForm(flux_uS_DvM_cont);
            assemble_->AddBiLinearForm(pnlt_uS_vS_cont);
            
            // check for prestressing
            // we have to do it once again in order to create a complex-valued coefficient function;
            // the use of a real-valued function goes well for an analytical function, but not for a compound one
            PtrParamNode preStressNode = bcNode->GetByVal("preStress", "region",
                    ptGrid_->GetRegion().ToString(volMasterId), ParamNode::PASS);
            if (preStressNode)
            {
              PtrCoefFct preStressFct = regionPreStress_[volMasterId];
              if (coefFuncPMLScl)
              {
                preStressFct = CoefFunction::Generate(mp_, Global::COMPLEX,
                        CoefXprTensScalOp(mp_, regionPreStress_[volMasterId], coefFuncPMLScl, CoefXpr::OP_MULT));
              }
              
              // master & slave integrals with normal derivatives
              BiLinearForm *preStrFlux_DuM_vM = NULL;
              BiLinearForm *preStrFlux_DuM_vS = NULL;
              BiLinearForm *preStrFlux_uM_DvM = NULL;
              BiLinearForm *preStrFlux_uS_DvM = NULL;
              
              // define bilinear forms for Nitsche coupling
              if (preStressFct->IsComplex())
              {
                preStrFlux_DuM_vM = GetFluxIntegrator<Complex>(one, coefFuncPMLVec, -1.0, BiLinearForm::PRIM_PRIM, true, icModes, true);
                preStrFlux_uM_DvM = GetFluxIntegrator<Complex>(factor, coefFuncPMLVec, -1.0, BiLinearForm::PRIM_PRIM, false, icModes, true);
                preStrFlux_DuM_vS = GetFluxIntegrator<Complex>(factor, coefFuncPMLVec, 1.0, BiLinearForm::PRIM_SEC, true, icModes, true);
                preStrFlux_uS_DvM = GetFluxIntegrator<Complex>(one, coefFuncPMLVec, 1.0, BiLinearForm::SEC_PRIM, false, icModes, true);
              }
              else
              {
                preStrFlux_DuM_vM = GetFluxIntegrator<Double>(one, coefFuncPMLVec, -1.0, BiLinearForm::PRIM_PRIM, true, icModes, true);
                preStrFlux_uM_DvM = GetFluxIntegrator<Complex>(factor, coefFuncPMLVec, -1.0, BiLinearForm::PRIM_PRIM, false, icModes, true);
                preStrFlux_DuM_vS = GetFluxIntegrator<Complex>(factor, coefFuncPMLVec, 1.0, BiLinearForm::PRIM_SEC, true, icModes, true);
                preStrFlux_uS_DvM = GetFluxIntegrator<Double>(one, coefFuncPMLVec, 1.0, BiLinearForm::SEC_PRIM, false, icModes, true);
              }
              
              preStrFlux_DuM_vM->SetBCoefFunctionOpA(preStressFct);
              preStrFlux_uM_DvM->SetBCoefFunctionOpB(preStressFct);
              preStrFlux_DuM_vS->SetBCoefFunctionOpA(preStressFct);
              preStrFlux_uS_DvM->SetBCoefFunctionOpB(preStressFct);
              
              // master-master
              preStrFlux_DuM_vM->SetName("preStrFlux_DuM_vM");
              preStrFlux_uM_DvM->SetName("preStrFlux_uM_DvM");
              //master-slave
              preStrFlux_DuM_vS->SetName("preStrFlux_DuM_vS");
              // slave-master
              preStrFlux_uS_DvM->SetName("preStrFlux_uS_DvM");
              
              SurfaceBiLinFormContext *preStrFlux_DuM_vM_cont = new SurfaceBiLinFormContext(preStrFlux_DuM_vM, STIFFNESS, BiLinearForm::PRIM_PRIM);
              SurfaceBiLinFormContext *preStrFlux_uM_DvM_cont = new SurfaceBiLinFormContext(preStrFlux_uM_DvM, STIFFNESS, BiLinearForm::PRIM_PRIM);
              SurfaceBiLinFormContext *preStrFlux_DuM_vS_cont = new SurfaceBiLinFormContext(preStrFlux_DuM_vS, STIFFNESS, BiLinearForm::PRIM_SEC);
              SurfaceBiLinFormContext *preStrFlux_uS_DvM_cont = new SurfaceBiLinFormContext(preStrFlux_uS_DvM, STIFFNESS, BiLinearForm::SEC_PRIM);
              
              preStrFlux_DuM_vM_cont->SetEntities(actSDList, actSDList);
              preStrFlux_uM_DvM_cont->SetEntities(actSDList, actSDList);
              preStrFlux_DuM_vS_cont->SetEntities(actSDList, actSDList);
              preStrFlux_uS_DvM_cont->SetEntities(actSDList, actSDList);
              
              preStrFlux_DuM_vM_cont->SetFeFunctions(feFunctions_[MECH_DISPLACEMENT], feFunctions_[MECH_DISPLACEMENT]);
              preStrFlux_uM_DvM_cont->SetFeFunctions(feFunctions_[MECH_DISPLACEMENT], feFunctions_[MECH_DISPLACEMENT]);
              preStrFlux_DuM_vS_cont->SetFeFunctions(feFunctions_[MECH_DISPLACEMENT], feFunctions_[MECH_DISPLACEMENT]);
              preStrFlux_uS_DvM_cont->SetFeFunctions(feFunctions_[MECH_DISPLACEMENT], feFunctions_[MECH_DISPLACEMENT]);
              
              assemble_->AddBiLinearForm(preStrFlux_DuM_vM_cont);
              assemble_->AddBiLinearForm(preStrFlux_uM_DvM_cont);
              assemble_->AddBiLinearForm(preStrFlux_DuM_vS_cont);
              assemble_->AddBiLinearForm(preStrFlux_uS_DvM_cont);
            }
          } // end nitsche
          else if (formulation == "Mortar")
          {
            EXCEPTION("Mortar coupling is not implemented for the PDE '" << this->GetName() << "'. Use Nitsche's one!");
          } // end mortar
          else
          {
            EXCEPTION("Unknown formulation: '" << formulation << "'!");
          }
        }
      }
    }
  }
  
  
  
  
  void MechPDE::DefineNcIntegrators() {
    StdVector< NcInterfaceInfo >::iterator ncIt = ncInterfaces_.Begin(),
            endIt = ncInterfaces_.End();
    for ( ; ncIt != endIt; ++ncIt ) {
      switch (ncIt->type) {
        case NC_MORTAR:
          if (dim_ == 2)
            DefineMortarCoupling<2,2>(MECH_DISPLACEMENT, *ncIt);
          else
            DefineMortarCoupling<3,3>(MECH_DISPLACEMENT, *ncIt);
          break;
        case NC_NITSCHE:
          if(dim_ == 2)
            DefineNitscheCoupling<2,2>(MECH_DISPLACEMENT, *ncIt);
          else
            DefineNitscheCoupling<3,3>(MECH_DISPLACEMENT, *ncIt);
          break;
        default:
          EXCEPTION("Unknown type of ncInterface");
          break;
      }
    }
  }
  
  void MechPDE::DefineConcentratedElems() {
    
    // Get FESpace and FeFunction of mechanical displacement
    shared_ptr<BaseFeFunction> myFct = feFunctions_[MECH_DISPLACEMENT];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
    
    // try to get bcsAndLoads node
    PtrParamNode bcNode = myParam_->Get("bcsAndLoads", ParamNode::PASS);
    if( !bcNode )
      return;
    
    // fetch parameter node specifying spring
    ParamNodeList springNodes = bcNode->GetList("concentratedElem");
    
    // Iterate over all springs
    std::string name, dofName;
    std::string massVal, dampVal, stiffVal;
    for( UInt i = 0; i < springNodes.GetSize(); i++ ) {
      
      // get data from node
      springNodes[i]->GetValue( "name", name );
      springNodes[i]->GetValue( "dof", dofName );
      springNodes[i]->GetValue( "massValue", massVal );
      springNodes[i]->GetValue( "dampingValue", dampVal );
      springNodes[i]->GetValue( "stiffnessValue", stiffVal );
      
      UInt dof = results_[0]->GetDofIndex( dofName );
      
      shared_ptr<EntityList> nodes = ptGrid_->GetEntityList(EntityList::NODE_LIST, name);
      UInt numNodes = nodes->GetSize();
      
      // Ensure, that only lists with 1 or 2 nodes are present
      if( numNodes > 2 ) {
        if(myFct->HasConstraint(name,dof)) {
          numNodes = 1; // is technically reduced to one node
        } else {
          WARN( "Concentrated mechanical element on '"  << name << "' is omitted, as it consists of more than " << "2 nodes!"; );
          continue;
        }
      }
      
      StdVector<FEMatrixType> matrices;
      matrices = STIFFNESS, MASS, DAMPING;
      StdVector<std::string> vals;
      vals = stiffVal, massVal, dampVal;
      
      if( numNodes == 1 ) {
        // ============================
        //  POINT CONCENTRATED ELEMENT
        // ============================
        
        for( UInt i = 0; i < matrices.GetSize(); ++i ) {
          // if value is zero, just continue
          if( vals[i] == "0" || vals[i] == "0.0" ) {
            continue;
          }
          
          SingleEntryBiLinInt * myInt = new SingleEntryBiLinInt(dim_, vals[i], dof, mp_);
          BiLinFormContext * intCtx = new BiLinFormContext(myInt, matrices[i]);
          intCtx->SetEntities( nodes, nodes );
          intCtx->SetFeFunctions( myFct, myFct );
          
          assemble_->AddBiLinearForm( intCtx );
        } // loop over mass/stiffness/damp
        
      } else {
        // =================================
        //  PAIR-WISE CONCENTRATED ELEMENTS
        // =================================
        
        // extract both nodes and put them into a new node list
        shared_ptr<NodeList> node1(new NodeList(ptGrid_));
        shared_ptr<NodeList> node2(new NodeList(ptGrid_));
        StdVector<UInt> tmp(1);
        EntityIterator it = nodes->GetIterator();
        tmp[0] = it.GetNode();
        node1->SetNodes( tmp );
        it++;
        tmp[0] = it.GetNode();
        node2->SetNodes( tmp );
        
        // loop over stiffness, mass and dampingvalue
        std::string diagVal, offDiagVal;
        for( UInt i = 0; i < matrices.GetSize(); ++ i ) {
          // if value is zero, just continue
          if( vals[i] == "0" || vals[i] == "0.0" ) {
            continue;
          }
          if( matrices[i] == STIFFNESS ) {
            diagVal = "1.0 * (" + vals[i] + ")";
            offDiagVal = "-1.0 * (" + vals[i] + ")";
          } else {
            diagVal = "(" + vals[i] + ") * 1.0 / 3.0";
            offDiagVal = "(" + vals[i] + ") * 1.0 / 6.0";
          }
          
          // a) diagonal entries
          SingleEntryBiLinInt * diagInt1 = new SingleEntryBiLinInt( dim_, diagVal, dof, mp_);
          BiLinFormContext * diagCtx1 =
                  new BiLinFormContext( diagInt1, matrices[i] );
          diagCtx1->SetEntities( node1, node1 );
          diagCtx1->SetFeFunctions( myFct, myFct );
          assemble_->AddBiLinearForm( diagCtx1 );
          
          SingleEntryBiLinInt * diagInt2 = new SingleEntryBiLinInt( dim_, diagVal, dof, mp_);
          BiLinFormContext * diagCtx2 =
                  new BiLinFormContext( diagInt2, matrices[i] );
          diagCtx2->SetEntities( node2, node2 );
          diagCtx2 ->SetFeFunctions( myFct, myFct );
          assemble_->AddBiLinearForm( diagCtx2 );
          
          // b) off-diagonal entries
          SingleEntryBiLinInt * offDiagInt = 
                  new SingleEntryBiLinInt( dim_, offDiagVal, dof, mp_);
          BiLinFormContext * offDiagCtx = 
                  new BiLinFormContext( offDiagInt, matrices[i] );
          offDiagCtx->SetEntities( node1, node2 );
          offDiagCtx->SetFeFunctions( myFct, myFct );
          offDiagCtx->SetCounterPart( true );
          assemble_->AddBiLinearForm( offDiagCtx );
        } // loop matrix types
      } // if: 1 / 2 nodes
      
    } // loop concentrated elements
    
  }
  
  void MechPDE::DefineRhsLoadIntegrators(PtrParamNode input)
  {
    LOG_TRACE(mechpde) << "Defining rhs load integrators for mechanic PDE";
    
    // Get FESpace and FeFunction of mechanical displacement
    shared_ptr<BaseFeFunction> myFct = feFunctions_[MECH_DISPLACEMENT];
    shared_ptr<FeSpace> mySpace = myFct->GetFeSpace();
    
    StdVector<shared_ptr<EntityList> > ent;
    StdVector<PtrCoefFct > coef;
    LinearForm * lin = NULL;
    StdVector<std::string> dispDofNames = myFct->GetResultInfo()->dofNames;
    
    /* Flag, if coefficient function lives on updated geometry (updated 
     * Lagrangian formulation). For analytically prescribed values
     * (pressure, force, stress), this is in general not the case. 
     * In the coupled case, the magnetic force density however could live
     * on an updated geometry. In this case, the ReadRhsExcitation()
     * method will set the coefUpdateGeo flag to true, so the RHS-
     * integrator also has to use this value.
     * 
     */
    bool coefUpdateGeo = false;
    
    // ========================
    //  FORCES (volume, nodal)
    // ========================
    ReadRhsExcitation("force", dispDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo, input);
    
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      LOG_DBG(mechpde) << "DRLI: reading force " << i << " coef=" << coef[i]->GetName() << " val=" << coef[i]->ToString() << " dep=" << CoefFunction::coefDependType.ToString(coef[i]->GetDependency());

      // check type of entity list
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        // --------------
        //  Nodal Forces 
        // --------------
        UInt numNodes = ent[i]->GetSize();
        // If more than one node is defined, we divide the total force by the number
        // of nodes to ensure that the total force is applied, independent of the 
        // number of nodes.
        // Note that is also happens when we have an expression as function.
        // Run with -d and check the .info.xml
        if(numNodes > 1 && coef[i]->DoNormalize()) {
          Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;  
          coef[i] = CoefFunction::Generate(mp_, part, CoefXprVecScalOp(mp_, coef[i],
                  boost::lexical_cast<std::string>(numNodes), CoefXpr::OP_DIV) );
        }
        
        lin = new SingleEntryInt(coef[i]);
        lin->SetName("NodalForceInt");
        LinearFormContext *ctx = new LinearFormContext( lin );
        ctx->SetEntities( ent[i] );
        ctx->SetFeFunction(myFct);
        assemble_->AddLinearForm(ctx);
        
      } else {
        
        // -------------------------
        //  Surface / Volume Forces 
        // -------------------------
        EXCEPTION("Not yet implemented: Surface / Volume Forces")
                
                // Same issue here as above: We need to "divide" the total force by the
                // area / volume to get the force density.
      }
    } // for
    
    // ===============
    //  FORCE DENSITY 
    // ===============
    LOG_DBG(mechpde) << "Reading force densities";
    
    ReadRhsExcitation("forceDensity", dispDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo, input);
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Force density must be defined on elements")
      }
      
      if( dim_ == 2) {
        if(isComplex_) {
          lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,2,2>(),
                  Complex(1.0), coef[i], coefUpdateGeo);
        } else {
          if(coef[i]->IsComplex()){
            // rhs excitation comes from a harmonic simulation (actually a real result with zero imaginary part)
            lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,2,2>(), Complex(1.0), coef[i], coefUpdateGeo, true, coef[i]->IsComplex() );
            WARN("<forceDensity> has complex-valued input. "
              "Only the real part is used. This is typically intended for mean values "
              "from harmonic analyses. Please verify if this behavior is appropriate for your use case.");
          }else{
            lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,2,2>(), 1.0, coef[i], coefUpdateGeo, true, coef[i]->IsComplex() );
          }
        }
      } else  {
        if(isComplex_) {
          lin = new BUIntegrator<Complex>(new IdentityOperator<FeH1,3,3>(),
                  Complex(1.0), coef[i], coefUpdateGeo);
        } else {
          if(coef[i]->IsComplex()){
            // rhs excitation comes from a harmonic simulation (actually a real result with zero imaginary part)
            lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,3,3>(), Complex(1.0), coef[i], coefUpdateGeo, true, coef[i]->IsComplex() );
            WARN("<forceDensity> has complex-valued input. "
              "Only the real part is used. This is typically intended for mean values "
              "from harmonic analyses. Please verify if this behavior is appropriate for your use case.");
          }else{
            lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,3,3>(), 1.0, coef[i], coefUpdateGeo, true, coef[i]->IsComplex() );
          }
        }
      }
      lin->SetName("ForceDensityInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    } // for

    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;

    // ===============
    //  VOLUME ACCELERATION
    // ===============
    LOG_DBG(mechpde) << "Reading volume acceleration";

    ReadRhsExcitation("volumeAcceleration", dispDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo, input);
    for (UInt i = 0; i < ent.GetSize(); ++i) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Volume acceleration must be defined on elements")
      }
      // get the material density
      RegionIdType curRegionId = ent[i]->GetRegion();
      BaseMaterial *curMaterial = materials_[curRegionId];
      PtrCoefFct density;
      if (this->IsMaterialComplex()) {
        // complex material definition
        EXCEPTION("Volume acceleration is not implemented for complex-valued materials");
      }
      else {
        density = curMaterial->GetScalCoefFnc(DENSITY, part);
      }
      // multiply by the acceleration to get the force density
      PtrCoefFct forceDens = CoefFunction::Generate(mp_, part, CoefXprVecScalOp(mp_, coef[i], density, CoefXpr::OP_MULT));

      if (isComplex_) {
        // frequency domain or eigenvalue analysis
        EXCEPTION("This feature is untested. Consider adding a test case before activating it.");
        if (dim_ == 2) {
          lin = new BUIntegrator<Complex>(new IdentityOperator<FeH1, 2, 2>(), Complex(1.0), forceDens, coefUpdateGeo);
        }
        else {
          lin = new BUIntegrator<Complex>(new IdentityOperator<FeH1, 3, 3>(), Complex(1.0), forceDens, coefUpdateGeo);
        }
      }
      else {
        if (coef[i]->IsComplex()) {
          // rhs excitation comes from a harmonic simulation (actually a real result with zero imaginary part)
          EXCEPTION("This feature is untested. Consider adding a test case before activating it.");
          WARN("<volumeAcceleration> has complex-valued input. "
               "Only the real part is used. This is typically intended for mean values "
               "from harmonic analyses. Please verify if this behavior is appropriate for your use case.");
          if (dim_ == 2) {
            lin = new BUIntegrator<Complex>(new IdentityOperator<FeH1, 2, 2>(), Complex(1.0), forceDens, coefUpdateGeo);
          }
          else {
            lin = new BUIntegrator<Complex>(new IdentityOperator<FeH1, 3, 3>(), Complex(1.0), forceDens, coefUpdateGeo, true, forceDens->IsComplex());
          }
        }
        else {
          if (dim_ == 2) {
            lin = new BUIntegrator<Double>(new IdentityOperator<FeH1, 2, 2>(), 1.0, forceDens, coefUpdateGeo, true, forceDens->IsComplex());
          }
          else {
            lin = new BUIntegrator<Double>(new IdentityOperator<FeH1, 3, 3>(), 1.0, forceDens, coefUpdateGeo, true, forceDens->IsComplex());
          }
        }
      }
      lin->SetName("VolumeAccelerationInt");
      LinearFormContext *ctx = new LinearFormContext(lin);
      ctx->SetEntities(ent[i]);
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    } // for

    // ===============
    //  magneticFluxDensity (couples in case of magnetostrictive material)
    // ===============
    LOG_DBG(mechpde) << "Reading magnetic flux density";
    
    ReadRhsExcitation( "magFluxDensity", dispDofNames, ResultInfo::VECTOR, isComplex_, 
            ent, coef, coefUpdateGeo );
    
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      
      // get region and material from entity
      RegionIdType curRegionId = ent[i]->GetRegion();
      BaseMaterial* curMaterial = NULL;
      bool complexMat = complexMatData_[curRegionId];
      curMaterial = materials_[curRegionId];
      
      //Set tensor type
      SubTensorType tensorType = GetSubTensorType();
      
      shared_ptr<CoefFunction > curCoef;
      
      
      if( complexMat ) {
        curCoef = curMaterial->GetTensorCoefFnc(MAGNETOSTRICTION_TENSOR_h_mech, tensorType, 
                Global::COMPLEX, true  );
      } else {
        curCoef = curMaterial->GetTensorCoefFnc(MAGNETOSTRICTION_TENSOR_h_mech, tensorType, 
                Global::REAL, true );
      }
      
      // std::cout << "coefUpdateGeoMech " << coefUpdateGeo << std::endl;
      
      if ( complexMat) {
        // explanation for factor 1.0:
        // the used coupling equations for magnotostriction are:
        // H = -hS + nuB
        // sigma = cS - hB
        // Bringing the term -hB to the right hand side would lead to +hB BUT
        // due to the prior weak formulation we lost one -1 factor so that we have -hB on the RHS.
        // However, we have -div(sigma) as in the mechanical case such that we get an additional -1 again leading to +hB
        //Complex factor = Complex(1.0);
        Complex factor = Complex(1.0);
        if( subType_ == "axi" ) {

          lin = new BDUIntegrator<StrainOperatorAxi<FeH1,Complex>, Complex> (factor,coef[i],curCoef,coefUpdateGeo);
          /*  		
           lin = new BUIntegrator<Complex> (new StrainOperatorAxi<FeH1,Complex>(),
           Complex(-1.0), CoefFunction::Generate( mp_, part, 
           CoefXprBinOp(mp_,  curCoef, coef[i], CoefXpr::OP_MULT) ), coefUpdateGeo);*/
          
          /*	
           integ = new BUIntegrator<Complex><Complex>(new StrainOperatorAxi<FeH1,Complex>(),
           new CurlOperatorAxi<Complex>(),
           //new GradientOperator<FeH1,2,1,Complex>(),
           curCoef, factor, true );*/
        } else if( subType_ == "planeStrain" ) {
          lin = new BDUIntegrator<StrainOperator2D<FeH1,Complex>, Complex> (factor,coef[i],curCoef,coefUpdateGeo);
          /*
           lin = new BUIntegrator<Complex> (new StrainOperator2D<FeH1,Complex>(),
           Complex(-1.0), CoefFunction::Generate( mp_, part, 
           CoefXprBinOp(mp_,  curCoef, coef[i], CoefXpr::OP_MULT) ), coefUpdateGeo);*/
          
        } else if( subType_ == "planeStress" ) {
          lin = new BDUIntegrator<StrainOperator2D<FeH1,Complex>, Complex> (factor,coef[i],curCoef,coefUpdateGeo);
          /*
           lin = new BUIntegrator<Complex> (new StrainOperator2D<FeH1,Complex>(),
           Complex(-1.0), CoefFunction::Generate( mp_, part, 
           CoefXprBinOp(mp_,  curCoef, coef[i], CoefXpr::OP_MULT) ), coefUpdateGeo);*/
          
        } else if( subType_ == "3d") {
          lin = new BDUIntegrator<StrainOperator3D<FeH1,Complex>, Complex> (factor,coef[i],curCoef,coefUpdateGeo);
          /*lin = new BUIntegrator<Complex> (new StrainOperator3D<FeH1,Complex>(),
           Complex(-1.0), CoefFunction::Generate( mp_, part,
           CoefXprBinOp(mp_,  curCoef, coef[i], CoefXpr::OP_MULT) ), coefUpdateGeo); */
        } else {
          EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
        }
      }
      else {
        //Double factor = 1.0;
        Double factor = 1.0;
        if( subType_ == "axi" ) {		
          
          lin = new BDUIntegrator<StrainOperatorAxi<FeH1,Double>, Double> (factor,coef[i],curCoef,coefUpdateGeo);
          /*  		
           lin = new BUIntegrator<Double> (new StrainOperatorAxi<FeH1,Double>(),
           Double(-1.0), CoefFunction::Generate( mp_, part, 
           CoefXprBinOp(mp_,  curCoef, coef[i], CoefXpr::OP_MULT) ), coefUpdateGeo);*/
          
          /*	
           integ = new BUIntegrator<Double><Double>(new StrainOperatorAxi<FeH1,Double>(),
           new CurlOperatorAxi<Double>(),
           //new GradientOperator<FeH1,2,1,Double>(),
           curCoef, factor, true );*/
        } else if( subType_ == "planeStrain" ) {
          lin = new BDUIntegrator<StrainOperator2D<FeH1,Double>, Double> (factor,coef[i],curCoef,coefUpdateGeo);
          /*
           lin = new BUIntegrator<Double> (new StrainOperator2D<FeH1,Double>(),
           Double(-1.0), CoefFunction::Generate( mp_, part, 
           CoefXprBinOp(mp_,  curCoef, coef[i], CoefXpr::OP_MULT) ), coefUpdateGeo);*/
          
        } else if( subType_ == "planeStress" ) {
          lin = new BDUIntegrator<StrainOperator2D<FeH1,Double>, Double> (factor,coef[i],curCoef,coefUpdateGeo);
          /*
           lin = new BUIntegrator<Double> (new StrainOperator2D<FeH1,Double>(),
           Double(-1.0), CoefFunction::Generate( mp_, part, 
           CoefXprBinOp(mp_,  curCoef, coef[i], CoefXpr::OP_MULT) ), coefUpdateGeo);*/
          
        } else if( subType_ == "3d") {
          lin = new BDUIntegrator<StrainOperator3D<FeH1,Double>, Double> (factor,coef[i],curCoef,coefUpdateGeo);
          /*lin = new BUIntegrator<Double> (new StrainOperator3D<FeH1,Double>(),
           Double(-1.0), CoefFunction::Generate( mp_, part,
           CoefXprBinOp(mp_,  curCoef, coef[i], CoefXpr::OP_MULT) ), coefUpdateGeo); */
        } else {
          EXCEPTION( "Subtype '" << subType_ << "' unknown for mechanic physic" );
        }
      }
      
      lin->SetName("magneticFluxDensityInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
      
    } // for
    
    // ==================
    //  THERMAL STRAIN
    // ==================
    LOG_DBG(mechpde) << "Reading thermal strain definition";
    StdVector<PtrCoefFct > tCoef;
    ReadRhsExcitation("thermalStrain", dispDofNames, ResultInfo::SCALAR, isComplex_, ent, tCoef, coefUpdateGeo, input);
    
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Thermal strain must be defined on elements")
      }
      
      RegionIdType myRegionId = ent[i]->GetRegion();
      // set sub tensor type
      SubTensorType subType = GetSubTensorType();
      
      // get stiffness tensor and thermal expansion coefficient (reduced to problem dim)
      PtrCoefFct cCoef, aCoef;
      if( isComplex_ ) {
        cCoef = materials_[myRegionId]->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, subType, Global::COMPLEX);
        aCoef = materials_[myRegionId]->GetSubVectorCoefFnc(MECH_THERMAL_EXPANSION_TENSOR, subType, Global::COMPLEX);
      }
      else {
        cCoef = materials_[myRegionId]->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, subType, Global::REAL);
        aCoef = materials_[myRegionId]->GetSubVectorCoefFnc(MECH_THERMAL_EXPANSION_TENSOR, subType, Global::REAL);
      }
      
      LOG_DBG3(mechpde)<< "Generating dT" << std::endl;
      // get reference temperature and compute dT = T - T_ref
      PtrCoefFct refTemp = materials_[myRegionId]->GetScalCoefFnc(MECH_TE_REFTEMPERATURE, part);
      PtrCoefFct dT = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_,tCoef[i],refTemp,CoefXpr::OP_SUB));
      
      LOG_DBG3(mechpde)<< "Generating thermalStrainCoef" << std::endl;
      // compute the thermal strain eps_th = alpha*dT

      // Using CoefXprBinOp instead of CoefXprVecScalOp
      PtrCoefFct thermalStrainCoef = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_,dT,aCoef,CoefXpr::OP_MULT));
      // store the coefFunction for postprocessing
      thermalStrain_->AddRegion( myRegionId, thermalStrainCoef);
      
      // Compute thermal stress (C*alpha)*dT
      // Note: the order of combining coef functions is important: C*alpha can usually be evaluated analytically.
      // Therefore, we combine it first. Then it is multiplied with dT, which might be a CoefFuctionGrid. 

      PtrCoefFct Calpha = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_,cCoef,aCoef,CoefXpr::OP_MULT));

      LOG_DBG3(mechpde)<< "Generating thermalStressCoef" << std::endl;
      // Using CoefXprBinOp instead of CoefXprVecScalOp
      PtrCoefFct thermalStressCoef = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_,dT,Calpha, CoefXpr::OP_MULT));

      // store the coefFunction for postprocessing
      thermalStress_->AddRegion( myRegionId, thermalStressCoef);
      
      BaseBOperator * bOp = GetStrainOperator( isComplex_, false );
      if(isComplex_) {
        lin = new BUIntegrator<Complex>(bOp, 1.0, thermalStressCoef, coefUpdateGeo);
      }
      else {
        lin = new BUIntegrator<Double>(bOp, 1.0, thermalStressCoef, coefUpdateGeo);
      }
      lin->SetName("ThermalStrainInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    } // for
    
    //================
    //ELECTRIC STRAIN (for forwardcoupling the piezo problem)
    //================
    LOG_DBG(mechpde)<< "Reading electric strain definition";
    //excitation
    StdVector<PtrCoefFct > exCoef;
    //direction
    StdVector<PtrCoefFct > dirCoef;
    //dependent
    StdVector<PtrCoefFct > depCoef;
    //formulation
    std::string formulation {"none"};
    /*
     * !!! Can currently only defined once in XML-schema!
     * Can be extended by ParamList->Get("bcsAndLoads")->GetChildren() and loop over children and
     * define rhs for each elecStrain entry
     *
     * Here we want to implement a pretty complex rhs, which is defined by:
     *
     * * excitation quantity: The main excitation. In our case the electric field
     * * direction: Just needed if you need to rotate the material tensor. Idea is that here you could
     *           specify for example the polarization, and the material tensor gets rotated in direction
     *           of the polarization
     * * scaling: Depends out of two parts:
     *            -Dependent scaling field quantity
     *            -formulation
     *
     * So the RHS looks something like this (without any testfunctions and diff operators):
     *       f(P) * Q(D) e Q(D) * E
     *
     *       scaling:
     *       f(P) = factor * |P| -> have to define P in dependentCoef (for example polarization
     *
     *       direction D: -> for example D = P/|P|
     *       Q(D) is the rotation matrix witch rotates e in D
     *
     *       E is the excitation vector.
    */

    //setting default PtrParamNode
    PtrParamNode elecStrainNode = myParam_->Get("bcsAndLoads");

    //be aware that node1->Get("ThisDoNotExcists,PASS)->Get("This->ThrowsError!")

    if(myParam_->Get("bcsAndLoads")->Has("elecStrain")){
      elecStrainNode = myParam_->Get("bcsAndLoads")->Get("elecStrain");

      //read in all the CoefFunctions and check if they are defined on the same entity
      //Excitation is checked as last one, because data is overwritten
      //and this is the only one which is not optional

      std::string entName {};

      ReadRhsExcitation("direction", dispDofNames, ResultInfo::VECTOR, isComplex_, ent,
          dirCoef, coefUpdateGeo,  elecStrainNode);

      //check if its defined not more than 1
      //if so, save Name
      if( ent.size() > 1 ){
        EXCEPTION("Electric Strain can (currently) only be assigned once in bcsAndLoads!")
      } else if (ent.size() == 1){
        entName= ent[0]->GetName();
      }

      ReadRhsExcitation("dependency", dispDofNames, ResultInfo::VECTOR, isComplex_, ent,
          depCoef, coefUpdateGeo, elecStrainNode->Get("scaling", ParamNode::PASS));

      //check if its defined not more than 1
      if( ent.size() > 1 ){
        EXCEPTION("Electric Strain can (currently) only be assigned once in bcsAndLoads!")
      } else if (ent.size() == 1){
        if(entName.empty()){
          entName= ent[0]->GetName();
        }
        else if(entName != ent[0]->GetName()){
          EXCEPTION("All electric strain sub-quantities must be defined on same region!")
        }
      }

      ReadRhsExcitation("excitation", dispDofNames, ResultInfo::VECTOR, isComplex_, ent,
          exCoef, coefUpdateGeo, elecStrainNode);

      if( ent.size() > 1 ){
        EXCEPTION("Electric Strain can (currently) only be assigned once in bcsAndLoads!")
      } else {
        if(entName.empty()){
          entName= ent[0]->GetName();
        }
        else if(ent[0]->GetName() != entName){
          EXCEPTION("All electric strain sub-quantities must be defined on same region!")
        }
      }

      //Get specific paramNode
      if(elecStrainNode->Has("scaling") && elecStrainNode->Get("scaling")->Has("formulation")){
        formulation = elecStrainNode->Get("scaling")->Get("formulation")->GetChild()->GetName();
      }

      LOG_DBG3(mechpde) << "Formulation: " << formulation << std::endl;

      // check type of entitylist
      if (ent[0]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("electric strain must be defined on elements") // Needed?
      }

      PtrCoefFct tmpCoef = exCoef[0];

      //Default value
      std::string factorString = "1";

      //if formulation exists, there must be also a depCoef
      if (formulation == "linear"){
        //linear scaling
        PtrCoefFct depNormCoef = CoefFunction::Generate( mp_, part, CoefXprUnaryOp(mp_,depCoef[0],CoefXpr::OP_NORM));
        //resulting coefFunction
        tmpCoef = CoefFunction::Generate( mp_, part, CoefXprBinOp(mp_,depNormCoef,exCoef[0],CoefXpr::OP_MULT));
        //can i read it directly as coeffunction?
        factorString = elecStrainNode->Get("scaling")->Get("formulation")->GetChild()->Get("factor")->As<std::string>();
      }

      PtrCoefFct factor = CoefFunction::Generate(mp_, part,factorString);

      PtrCoefFct resCoef = CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_,factor,tmpCoef,CoefXpr::OP_MULT));

      RegionIdType myRegionId = ent[0]->GetRegion();

      // set sub tensor type
      SubTensorType subType = GetSubTensorType();

      // Getting the piezoelectric couple tensor
      PtrCoefFct eCoef;
      if( isComplex_ ) {
        eCoef = materials_[myRegionId]->GetTensorCoefFnc(PIEZO_TENSOR, subType, Global::COMPLEX, true);
      }
      else {
        eCoef = materials_[myRegionId]->GetTensorCoefFnc(PIEZO_TENSOR, subType, Global::REAL, true);
      }

      //TODO: rotate e depending on dirCoef
      //===================================
      //Insert here a cool thingy which does that
      //neatRotationFunction(eCoef,dirCoef)

      if (dim_ == 3){
        if(isComplex_) {
          lin = new BDUIntegrator<StrainOperator3D<FeH1,Complex>, Complex> (1 ,resCoef,eCoef,coefUpdateGeo);
        } else {
          lin = new BDUIntegrator<StrainOperator3D<FeH1,Complex>, Double> (1 ,resCoef,eCoef,coefUpdateGeo);
        }
      } else {
        EXCEPTION("Currently only 3d case implemented! Dont worry its just copy pase.")
      }

      lin->SetName("ElecStrainInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[0] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[0]);
    }

    // ===============
    //  PRESSURE 
    // ===============
    LOG_DBG(mechpde) << "Reading mechanical pressure";
    StdVector<std::string> empty;
    ReadRhsExcitation("pressure", empty, ResultInfo::SCALAR, isComplex_, ent, coef, coefUpdateGeo, input);
    std::set<RegionIdType> volRegions (regions_.Begin(), regions_.End() );
    
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Mechanical pressure must be defined on elements")
      }
      
      // Factor for the pressure:
      // The pressure is by definition in the opposite direction as the 
      // normal stress, i.e. a positive pressure means a compressive stress
      // (<0). Thus we have to take the minus sign into account
      const Double presFac = -1.0;
      
      if( dim_ == 2) {
        if(isComplex_) {
          lin = new BUIntegrator<Complex, true> ( new IdentityOperatorNormalTrans<FeH1,2,1>(),
                  Complex(presFac), coef[i],
                  volRegions, coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double,true> ( new IdentityOperatorNormalTrans<FeH1,2,1>(),
                  presFac, coef[i], volRegions,
                  coefUpdateGeo);
        }
      } else  {
        if(isComplex_) {
          lin = new BUIntegrator<Complex, true> ( new IdentityOperatorNormalTrans<FeH1,3,1>(),
                  Complex(presFac), coef[i],
                  volRegions, coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double, true> ( new IdentityOperatorNormalTrans<FeH1,3,1>(),
                  presFac, coef[i],
                  volRegions, coefUpdateGeo);
        }
      }
      lin->SetName("PressureInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    } // for



    // ===============
    //  MECH_NORMAL_STRESS
    // ===============
    LOG_DBG(mechpde) << "Reading mechanical normal stress";
    ReadRhsExcitation("mechNormalStress", dispDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo, input);

    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Mechanical normal stress must be defined on elements")
      }
      // ensure that list contains only surface elements
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if( elemDim != (dim_-1) ) {
        EXCEPTION("Mechanical normal stress can only be defined on surface elements");
      }

      if( dim_ == 2) {
        if(isComplex_) {
          lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,2,2>(),
                  Complex(1.0), coef[i], coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,2,2>(),
                  1.0, coef[i],coefUpdateGeo);
        }
      } else  {
        if(isComplex_) {
          lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,3,3>(),
                  Complex(1.0), coef[i], coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,3,3>(),
                  1.0, coef[i], coefUpdateGeo);
        }
      }
      lin->SetName("normalStressInt");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    } // for
    
    
    
    // ==================
    //  SURFACE TRACTION  
    // ==================
    LOG_DBG(mechpde) << "Reading surface tractions";
    
    ReadRhsExcitation("traction", dispDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo, input);
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      // check type of entitylist
      if (ent[i]->GetType() == EntityList::NODE_LIST) {
        EXCEPTION("Surface traction must be defined on elements")
      }
      // ensure that list contains only surface elements
      EntityIterator it = ent[i]->GetIterator();
      UInt elemDim = Elem::shapes[it.GetElem()->type].dim;
      if( elemDim != (dim_-1) ) {
        EXCEPTION("Surface traction can only be defined on surface elements");
      }

      // check if we are iteratively coupled to the LinFlow-PDE and adapt sign if necessary
      std::string couplName;
      // go through the paramNodes to get the name of the coupled quantity (if there is any)
      ParamNodeList tracNodeList = myParam_->Get("bcsAndLoads")->GetList("traction");
      // ensure we have the correct node - preceeding checks are not necessary because this has already been handled by SinglePDE::ReadEntities
      PtrParamNode tracNode = tracNodeList[i];
      if( tracNode ) {
        PtrParamNode couplNode = tracNode->Get("coupling",ParamNode::PASS);
        if( couplNode ) {
          PtrParamNode quantNode = couplNode->Get("quantity",ParamNode::PASS);
          if( quantNode ) {
            quantNode->GetValue("name", couplName );
          }
        }
        
      }
      // for the fluidMechSurfaceTraction we have to multiply by -1 to be consistent with the normal vectors
      Double tracFac;
      if ( couplName == "fluidMechSurfaceTraction") {
        tracFac = -1.0;
      } else {
        tracFac = 1.0;
      }



      if( dim_ == 2) {
        if(isComplex_) {
          lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,2,2>(),
                  Complex(tracFac), coef[i], coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,2,2>(),
                  tracFac, coef[i],coefUpdateGeo);
        }
      } else  {
        if(isComplex_) {
          lin = new BUIntegrator<Complex> ( new IdentityOperator<FeH1,3,3>(),
                  Complex(tracFac), coef[i], coefUpdateGeo);
        } else {
          lin = new BUIntegrator<Double> ( new IdentityOperator<FeH1,3,3>(),
                  tracFac, coef[i], coefUpdateGeo);
        }
      }
      lin->SetName("TractionIntegrator");
      LinearFormContext *ctx = new LinearFormContext( lin );
      ctx->SetEntities( ent[i] );
      ctx->SetFeFunction(myFct);
      assemble_->AddLinearForm(ctx);
      myFct->AddEntityList(ent[i]);
    } // for
    
    // ==================
    //  RHS NODAL VALUES
    // ==================
    LOG_DBG(mechpde) << "Reading direct right hand side values";
    
    ReadRhsExcitation("rhsValues", dispDofNames, ResultInfo::VECTOR, isComplex_, ent, coef, coefUpdateGeo, input);
    
    for( UInt i = 0; i < ent.GetSize(); ++i ) {
      //for non-linear simulations we might need a conservative interpolation in each timestep...
      coef[i]->SetConservative(true);
      this->rhsFeFunctions_[MECH_DISPLACEMENT]->AddLoadCoefFunction(coef[i], ent[i]);
    }
    
    // =============
    //  TEST STRAIN
    // =============
    // the test strains are for homogenization, they make only sense combined with optimization but they can
    // also be set for simulation such that the impact can be studied
    if(myParam_->Has("bcsAndLoads/testStrain"))
    {
      ParamNodeList tsl = myParam_->Get("bcsAndLoads")->GetList("testStrain");
      for(unsigned int i = 0; i < tsl.GetSize(); i++)
        DefineTestStrainIntegrator(testStrain.Parse(tsl[i]->Get("strain")->As<std::string>()));
    }
  }
  
  void MechPDE::DefineTestStrainIntegrator(const TestStrain test, StdVector<LinearFormContext*>* linForms)
  {
    LOG_DBG(mechpde) << "DTSI: test=" << test << " lf=" << (linForms ? linForms->ToString() : "null");
    
    shared_ptr<BaseFeFunction> myFct = feFunctions_[MECH_DISPLACEMENT];
    
    // to generate the vector values coef function for the test strain we need scalar const coef functions for zero and one
    PtrCoefFct one  = CoefFunction::Generate(mp_, Global::REAL, "1.0");
    PtrCoefFct zero = CoefFunction::Generate(mp_, Global::REAL, "0.0");
    
    StdVector<PtrCoefFct> strain(dim_ == 2 ? 3 : 6);
    strain.Init(zero);
    strain[dim_ == 2 && test == MechPDE::XY ? MechPDE::Z : test] = one; // xy goes to the third element (z) for 2D
    LOG_DBG(mechpde) << "DTSI: idx=" << (dim_ == 2 && test == MechPDE::XY ? MechPDE::Z : test) << " -> one";
    
    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for(it = materials_.begin(); it != materials_.end(); it++)
    {
      // Set current region and material
      RegionIdType actRegion = it->first;
      
      shared_ptr<ElemList> actSDList(new ElemList(domain->GetGrid()));
      actSDList->SetRegion( actRegion );
      
      PtrCoefFct ts = CoefFunction::Generate(mp_, Global::REAL, strain);
      assert(regionStiffness_[actRegion]->GetDimType() == CoefFunction::TENSOR);
      assert(ts->GetDimType() == CoefFunction::VECTOR);
      LinearForm* lin = NULL;
      if(dim_ == 3)
        lin = new BDUIntegrator<StrainOperator3D<FeH1,double>, double>(1.0, ts, regionStiffness_[actRegion], false); // no updateGeo
      else
        lin = new BDUIntegrator<StrainOperator2D<FeH1,double>, double>(1.0, ts, regionStiffness_[actRegion], false); // no updateGeo
      LinearFormContext* ctx = new LinearFormContext(lin);
      ctx->SetEntities(actSDList);
      ctx->SetFeFunction(myFct);
      
      if(linForms != NULL)
        linForms->Push_back(ctx);
      else
        assemble_->AddLinearForm(ctx);
    }
  }
  
  BaseBDBInt* MechPDE::GetStiffIntegrator(BaseMaterial* actSDMat, RegionIdType regionId, bool isComplex)
  {
    // Get region name
    std::string regionName = ptGrid_->GetRegion().ToString( regionId );
    
    // ------------------------
    //  Obtain linear material
    // ------------------------
    shared_ptr<CoefFunction > curCoef;
    if( isComplex )
      curCoef = actSDMat->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, tensorType_, Global::COMPLEX);
    else
      curCoef = actSDMat->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, tensorType_, Global::REAL);
    
    // when we do optimization we wrap the original CoefFunction. Don't check for region to handle dim-1 pressure on dim elements
    if(domain->HasDesign())
    {
      CoefFunctionOpt* tmpFnc = new CoefFunctionOpt(domain->GetDesign(), curCoef, MECH_STIFFNESS_TENSOR, this); // takes double and complex
      curCoef.reset(tmpFnc);
    }
    
    // store coefficient function for later use (e.g. in boundary integrators)
    regionStiffness_[regionId] = curCoef;
    
    // ----------------------------------------
    //  Determine correct stiffness integrator 
    // ----------------------------------------
    
    BaseBDBInt* integ = NULL;
    BaseBOperator* bOp = GetStrainOperator(isComplex, false);
    
    if( regionSoftening_[regionId] == "icModesTW") {
      // ===================
      //  ICModes Softening 
      // ===================
      BaseBOperator * gOp = GetStrainOperator( isComplex, true );
      if (isComplex ) 
        integ = new ICModesInt<Complex>(bOp, gOp, curCoef, 1.0);
      else
        integ = new ICModesInt<Double>(bOp, gOp, curCoef, 1.0);
      
    } else {
      // ====================
      //  Standard Stiffness 
      // ====================
      if (isComplex ) 
        integ = new BDBInt<Complex>(bOp, curCoef, 1.0, updatedGeo_);
      else
        integ = new BDBInt<Double>(bOp, curCoef, 1.0, updatedGeo_);
    }
    
    // the integrator has a coef function but for the optimization case the opt coef needs to know also the integrator
    if(domain->HasDesign())
      dynamic_pointer_cast<CoefFunctionOpt>(curCoef)->SetForm(integ);
    
    return integ;
  }
  
  BaseBDBInt* MechPDE::GetStiffIntegrator(BaseMaterial* actSDMat, RegionIdType regionId, bool isComplex, PtrCoefFct scalingFactor)
  {
    BaseBDBInt* integ = NULL;
    BaseBOperator* bOp = NULL;
    
    // ------------------------
    //  Obtain linear material
    // ------------------------
    shared_ptr<CoefFunction > curCoef;
    
    // ----------------------------------------
    //  Determine correct stiffness integrator 
    // ----------------------------------------
    if (isComplex){
      curCoef = actSDMat->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, tensorType_, Global::COMPLEX);
      
      // store coefficient function for later use (e.g. in boundary integrators)
      regionStiffness_[regionId] = curCoef;
      PtrCoefFct curCoefScl = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprTensScalOp(mp_, curCoef, scalingFactor, CoefXpr::OP_MULT_TENSOR));
      
      if ((subType_ == "planeStrain") || (subType_ == "planeStress"))
        bOp = new ScaledStrainOperator2D<FeH1, Complex>();
      else if (subType_ == "2.5d")
        bOp = new ScaledStrainOperator2p5D<FeH1, Complex>();
      else if (subType_ == "3d")
        bOp = new ScaledStrainOperator3D<FeH1, Complex>();
      else
        EXCEPTION("Scaled strain operator is not implemented for the subtype: " << subType_);
      
      if (regionSoftening_[regionId] == "icModesTW")
      {
        BaseBOperator* gOp = bOp->Clone();
        integ = new ICModesInt<Complex>(bOp, gOp, curCoefScl, 1.0);
      }
      else
        integ = new BDBInt<Complex, Complex>(bOp, curCoefScl, 1.0, updatedGeo_);
      
    } else{
      curCoef = actSDMat->GetTensorCoefFnc(MECH_STIFFNESS_TENSOR, tensorType_, Global::REAL);
      
      // store coefficient function for later use (e.g. in boundary integrators)
      regionStiffness_[regionId] = curCoef;
      PtrCoefFct curCoefScl = CoefFunction::Generate(mp_, Global::REAL, CoefXprTensScalOp(mp_, curCoef, scalingFactor, CoefXpr::OP_MULT_TENSOR));
      
      if ((subType_ == "planeStrain") || (subType_ == "planeStress"))
        bOp = new ScaledStrainOperator2D<FeH1, Double>();
      else if (subType_ == "2.5d")
        bOp = new ScaledStrainOperator2p5D<FeH1, Double>();
      else if (subType_ == "3d")
        bOp = new ScaledStrainOperator3D<FeH1, Double>();
      else
        EXCEPTION("Scaled strain operator is not implemented for the subtype: " << subType_);
      
      if (regionSoftening_[regionId] == "icModesTW")
      {
        BaseBOperator* gOp = bOp->Clone();
        integ = new ICModesInt<Double>(bOp, gOp, curCoefScl, 1.0);
      }
      else
        integ = new BDBInt<Double, Double>(bOp, curCoefScl, 1.0, updatedGeo_);
    }
    

    return integ;
  }

  BaseBDBInt *MechPDE::GetKelvinVoigtDampingIntegrator(BaseMaterial *actSDMat, RegionIdType regionId, bool isComplex)
  {
    // Obtain linear viscous tensor (has been read in ReadDampingInformation)
    PtrCoefFct curCoef = regionKelvinVoigtDamping_[regionId];

    // Determine correct damping integrator (more or less copy-paste from GetStiffIntegrator)
    BaseBDBInt *integ = nullptr;
    BaseBOperator *bOp = GetStrainOperator(isComplex, false);

    if (regionSoftening_[regionId] == "icModesTW") {
      // ICModes Softening
      EXCEPTION("Kelvin-Voigt damping and icModes is not implemented/tested.");
      BaseBOperator *gOp = GetStrainOperator(isComplex, true);
      if (isComplex)
        integ = new ICModesInt<Complex>(bOp, gOp, curCoef, 1.0);
      else
        integ = new ICModesInt<Double>(bOp, gOp, curCoef, 1.0);
    }
    else {
      // Standard Kelvin Voigt-Damping (like stiffness)
      if (isComplex) {
        integ = new BDBInt<Complex>(bOp, curCoef, 1.0, updatedGeo_);
      }
      else {
        integ = new BDBInt<Double>(bOp, curCoef, 1.0, updatedGeo_);
      }
    }
    return integ;
  }

  BaseBDBInt* MechPDE::GetPreStressIntegrator(PtrCoefFct preStressFct, RegionIdType regionId, bool isComplex, Double factor)
  {
    BaseBDBInt *preStressInt = NULL;
    if(dim_==2 && subType_ != "2.5d"){
      if( regionSoftening_[regionId] == "icModesTW") {
        if(isComplex){
          preStressInt = new ICModesInt<Complex>(new PiolaStressOperator<FeH1,2,Complex>(true),
                  new PiolaStressOperator<FeH1,2,Complex>(true),preStressFct,factor);
        }else{
          preStressInt = new ICModesInt<Double>(new PiolaStressOperator<FeH1,2,Double>(true),
                  new PiolaStressOperator<FeH1,2,Double>(true),preStressFct,factor);
        }
      }else{
        if(isComplex){
          preStressInt = new BDBInt<Complex, Complex>(new PiolaStressOperator<FeH1,2,Complex>(false),preStressFct,factor);
        }else{
          preStressInt = new BDBInt<Double, Double>(new PiolaStressOperator<FeH1,2,Double>(false),preStressFct,factor);
        }
      }
    }else if (dim_==2 && subType_ == "2.5d"){
      if( regionSoftening_[regionId] == "icModesTW") {
        if(isComplex){
          preStressInt = new ICModesInt<Complex>(new PreStressOperator2p5D<FeH1,2,3,Complex>(true),
                  new PreStressOperator2p5D<FeH1,2,3,Complex>(true),preStressFct,factor);
        }else{
          preStressInt = new ICModesInt<Double>(new PreStressOperator2p5D<FeH1,2,3,Double>(true),
                  new PreStressOperator2p5D<FeH1,2,3,Double>(true),preStressFct,factor);
        }
      }else{
        if(isComplex){
          preStressInt = new BDBInt<Complex, Complex>(new PreStressOperator2p5D<FeH1,2,3,Complex>(false),preStressFct,factor);
        }else{
          preStressInt = new BDBInt<Double, Double>(new PreStressOperator2p5D<FeH1,2,3,Double>(false),preStressFct,factor);
        }
      }
    }else{
      if( regionSoftening_[regionId] == "icModesTW") {
        if(isComplex){
          preStressInt = new ICModesInt<Complex>(new PiolaStressOperator<FeH1,3,Complex>(true),
                  new PiolaStressOperator<FeH1,3,Complex>(true),preStressFct,factor);
        }else{
          preStressInt = new ICModesInt<Double>(new PiolaStressOperator<FeH1,3,Double>(true),
                  new PiolaStressOperator<FeH1,3,Double>(true),preStressFct,factor);
        }
      }else{
        if(isComplex){
          preStressInt = new BDBInt<Complex, Complex>(new PiolaStressOperator<FeH1,3,Complex>(false),preStressFct,factor);
        }else{
          preStressInt = new BDBInt<Double, Double>(new PiolaStressOperator<FeH1,3,Double>(false),preStressFct,factor);
        }
      }
    }

    // the integrator has a coef function but for the optimization case the opt coef needs to know also the integrator
    if (domain->HasDesign())
      dynamic_pointer_cast<CoefFunctionOpt>(preStressFct)->SetForm(preStressInt);

    return preStressInt;
  }
  
  BaseBDBInt* MechPDE::GetPreStressIntegrator(PtrCoefFct preStressFct, RegionIdType regionId, bool isComplex, PtrCoefFct scalingFactor)
  {
    BaseBDBInt* integ = NULL;
    BaseBOperator* bOp = NULL;
    
    PtrCoefFct curCoefScl = CoefFunction::Generate(mp_, Global::COMPLEX, CoefXprTensScalOp(mp_, preStressFct,
            scalingFactor, CoefXpr::OP_MULT_TENSOR));
    
    // ----------------------------------------
    //  Determine correct prestress integrator
    // ----------------------------------------
    if ((subType_ == "planeStrain") || (subType_ == "planeStress"))
      bOp = new ScaledPreStressOperator<FeH1, 2, Complex>();
    else if (subType_ == "2.5d")
      bOp = new ScaledPreStressOperator2p5D<FeH1, 2, 3, Complex>();
    else if (subType_ == "3d")
      bOp = new ScaledPreStressOperator<FeH1, 3, Complex>();
    else
      EXCEPTION("Scaled strain operator is not implemented for the subtype: " << subType_);
    
    if (regionSoftening_[regionId] == "icModesTW")
    {
      BaseBOperator* gOp = bOp->Clone();
      integ = new ICModesInt<Complex>(bOp, gOp, curCoefScl, 1.0);
    }
    else
      integ = new BDBInt<Complex, Complex>(bOp, curCoefScl, 1.0, updatedGeo_);
    
    return integ;
  }
  
  BaseBOperator* MechPDE::GetStrainOperator(bool isComplex, bool icModes)
  {
    BaseBOperator* bOp = NULL;
    // determine if we do bloch eigenfrequency analysis
    bool do_bloch = domain->GetDriver()->DoBlochModeEigenfrequency();
    assert(!(do_bloch && !isComplex));
    
    
    if(isComplex)
    {
      if( subType_ == "planeStrain" )
      {
        if(do_bloch)
        {
          bOp = new StrainOperatorBloch2D<FeH1,Complex>(icModes);
          EigenFrequencyDriver* efd = dynamic_cast<EigenFrequencyDriver*>(domain->GetSingleDriver());
          dynamic_cast<StrainOperatorBloch2D<FeH1,Complex>* >(bOp)->SetWaveVector(efd->GetCurrentWaveVector());
        }
        else
          bOp = new StrainOperator2D<FeH1,Complex>(icModes);
      }
      if( subType_ == "axi" )
        bOp = new StrainOperatorAxi<FeH1,Complex>(icModes);
      if(subType_ == "planeStress")
        bOp = new StrainOperator2D<FeH1,Complex>(icModes);
      if(subType_ == "3d")
      {
        if(do_bloch)
        {
          bOp = new StrainOperatorBloch3D<FeH1,Complex>(icModes);
          EigenFrequencyDriver* efd = dynamic_cast<EigenFrequencyDriver*>(domain->GetSingleDriver());
          dynamic_cast<StrainOperatorBloch3D<FeH1,Complex>* >(bOp)->SetWaveVector(efd->GetCurrentWaveVector());
        }
        else
          bOp = new StrainOperator3D<FeH1,Complex>(icModes);
      }
      if(subType_ == "2.5d")
        bOp = new StrainOperator2p5D<FeH1,Complex>(icModes);
    }
    else // not complex
    {
      if(subType_ == "axi")
        bOp = new StrainOperatorAxi<FeH1,Double>(icModes);
      if(subType_ == "planeStrain")
        bOp = new StrainOperator2D<FeH1,Double>(icModes);
      if(subType_ == "planeStress")
        bOp = new StrainOperator2D<FeH1,Double>(icModes);
      if(subType_ == "3d")
        bOp = new StrainOperator3D<FeH1,Double>(icModes);
      if(subType_ == "2.5d")
        bOp = new StrainOperator2p5D<FeH1,Double>(icModes);
    }
    
    if(bOp == NULL)
      EXCEPTION("strain operator not implemented for analysis type");
    
    return bOp;
  }
  
  template<typename DATA_TYPE>
  BiLinearForm* MechPDE::GetPenaltyIntegrator(PtrCoefFct scalCoefFunc, Double factor, BiLinearForm::CouplingDirection cplDir)
  {
    BiLinearForm* integ = NULL;
    
    if (dim_ == 2)
    {
      if (subType_ == "2.5d")
        integ = new SurfaceNitscheABInt<DATA_TYPE, DATA_TYPE>(new SurfaceIdentityOperator<FeH1, 2, 3>(),
                new SurfaceIdentityOperator<FeH1, 2, 3>(),
                scalCoefFunc, factor, cplDir, updatedGeo_, false, true);
      else
        integ = new SurfaceNitscheABInt<DATA_TYPE, DATA_TYPE>(new SurfaceIdentityOperator<FeH1, 2, 2>(),
                new SurfaceIdentityOperator<FeH1, 2, 2>(),
                scalCoefFunc, factor, cplDir, updatedGeo_, false, true);
    }
    else
    {
      integ = new SurfaceNitscheABInt<DATA_TYPE, DATA_TYPE>(new SurfaceIdentityOperator<FeH1, 3, 3>(),
              new SurfaceIdentityOperator<FeH1, 3, 3>(),
              scalCoefFunc, factor, cplDir, updatedGeo_, false, true);
    }
    
    return integ;
  }
  
  
  template<typename DATA_TYPE>
  BiLinearForm* MechPDE::GetFluxIntegrator(PtrCoefFct scalCoefFucn, PtrCoefFct coefFuncPMLVec, Double factor,
          BiLinearForm::CouplingDirection cplDir, bool fluxOpA, bool icModes, bool preStress)
  {
    BiLinearForm* integ = NULL;
    BaseBOperator *fluxOp = NULL, *idOp = NULL;
    
    // The flux operator will implement either a scaled differential operator if a non-zero 'coefFnc' has bee passed,
    // or a normal differential operator otherwise. The differential operator can be either 'SurfaceNormalStress' or
    // 'SurfaceNormalPreStress', depending on the flag 'preStress'
    if (dim_ == 3)
    {
      idOp = new SurfaceIdentityOperator<FeH1, 3, 3>();
      if (coefFuncPMLVec)
      {
        if (preStress)
          fluxOp = new SurfaceNormalPreStressOperator<FeH1, 3, 3, Complex>(subType_, coefFuncPMLVec, icModes);
        else
          fluxOp = new SurfaceNormalStressOperator<FeH1, 3, 3, Complex>(subType_, coefFuncPMLVec, icModes);
      }
      else
      {
        if (preStress)
          fluxOp = new SurfaceNormalPreStressOperator<FeH1, 3, 3, Double>(subType_, icModes);
        else
          fluxOp = new SurfaceNormalStressOperator<FeH1, 3, 3, Double>(subType_, icModes);
      }
    }
    else if (dim_ == 2 && subType_ == "2.5d")
    {
      idOp = new SurfaceIdentityOperator<FeH1, 2, 3>();
      if (coefFuncPMLVec)
      {
        if (preStress)
          fluxOp = new SurfaceNormalPreStressOperator<FeH1, 2, 3, Complex>(subType_, coefFuncPMLVec, icModes);
        else
          fluxOp = new SurfaceNormalStressOperator<FeH1, 2, 3, Complex>(subType_, coefFuncPMLVec, icModes);
      }
      else
      {
        if (preStress)
          fluxOp = new SurfaceNormalPreStressOperator<FeH1, 2, 3, Double>(subType_, icModes);
        else
          fluxOp = new SurfaceNormalStressOperator<FeH1, 2, 3, Double>(subType_, icModes);
      }
    }
    else
    {
      idOp = new SurfaceIdentityOperator<FeH1, 2, 2>();
      if (coefFuncPMLVec)
      {
        if (preStress)
          fluxOp = new SurfaceNormalPreStressOperator<FeH1, 2, 2, Complex>(subType_, coefFuncPMLVec, icModes);
        else
          fluxOp = new SurfaceNormalStressOperator<FeH1, 2, 2, Complex>(subType_, coefFuncPMLVec, icModes);
      }
      else
      {
        if (preStress)
          fluxOp = new SurfaceNormalPreStressOperator<FeH1, 2, 2, Double>(subType_, icModes);
        else
          fluxOp = new SurfaceNormalStressOperator<FeH1, 2, 2, Double>(subType_, icModes);
      }
    }
    
    // Check whether we have a du/dn*v or a u*dv/dn bilinear form
    if (fluxOpA)
      integ = new SurfaceNitscheABInt<DATA_TYPE, DATA_TYPE>(fluxOp, idOp, scalCoefFucn, factor, cplDir, updatedGeo_);
    else
      integ = new SurfaceNitscheABInt<DATA_TYPE, DATA_TYPE>(idOp, fluxOp, scalCoefFucn, factor, cplDir, updatedGeo_);
    
    return integ;
  }
  
  void MechPDE::DefineSolveStep()
  {
    solveStep_ = new StdSolveStep(*this);
  }

  const Matrix<double>& MechPDE::GetVonMisesMatrix(int dim)
  {
    Matrix<double>& m = dim == 2 ? vonMisesMatrix_2d_ : vonMisesMatrix_3d_;
    if(m.GetNumRows() == 0)
    {
      // Kocvara and Stingl; 2007 -> von Mises Stress = stress^T * M * stress
      if(dim == 2)
      {
        m.Resize(3,3);
        m.Init();
        m[0][0] = 1.0;
        m[1][1] = 1.0;
        m[2][2] = 3.0;
        m[0][1] = -0.5;
        m[1][0] = -0.5;

      }
      else
      {
        m.Resize(6,6);
        m.Init();
        m[0][0] = 2.0;
        m[1][1] = 2.0;
        m[2][2] = 2.0;
        m[3][3] = 6.0;
        m[4][4] = 6.0;
        m[5][5] = 6.0;
        m[0][1] = -1.0;
        m[0][2] = -1.0;
        m[1][2] = -1.0;
        m[1][0] = -1.0;
        m[2][0] = -1.0;
        m[2][1] = -1.0;
      }
    }
    return m;
  }

  const Matrix<double>& MechPDE::GetHillMandelMatrix(int dim)
  {
    Matrix<double>& m = dim == 2 ? hillMandelMatrix_2d_ : hillMandelMatrix_3d_;
    if(m.GetNumRows() == 0)
    {
      // Hill Mandel norm = stress^T * M * stress
      if(dim == 2)
      {
        m.Resize(3,3);
        m.Init();
        m[0][0] = 1.0;
        m[1][1] = 1.0;
        m[2][2] = 2.0;
      }
      else
      {
        m.Resize(6,6);
        m.Init();
        m[0][0] = 1.0;
        m[1][1] = 1.0;
        m[2][2] = 1.0;
        m[3][3] = 2.0;
        m[4][4] = 2.0;
        m[5][5] = 2.0;
      }
    }
    return m;
  }

  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================
  void MechPDE::InitTimeStepping()
  {
    // Check if time integration is defined in XML input
    PtrParamNode transientNode = myParam_->GetParent()->GetParent()->Get("analysis")->Get("transient", ParamNode::PASS);
    PtrParamNode integrationScheme = transientNode->Get("integrationScheme", ParamNode::PASS);

    PtrParamNode timeStepAlphaNode = this->myParam_->Get("timeStepAlpha", ParamNode::PASS);
    if (integrationScheme && timeStepAlphaNode)
      throw Exception("Both 'integrationScheme' and 'timeStepAlpha' are specified for the mechanical PDE. "
                      "Please use 'integrationScheme' only, as it provides more flexibility and "
                      "supersedes the legacy 'timeStepAlpha' parameter.");

    GLMScheme* scheme = nullptr;
    if (integrationScheme)
    {
      scheme = GetXmlDefinedScheme(integrationScheme);
    }
    else
    {
      Double alpha = this->myParam_->Get("timeStepAlpha")->As<Double>();
      scheme = new Newmark(0.5, 0.25, alpha);
    }

    TimeSchemeGLM::NonLinType nlType = (nonLin_) ? TimeSchemeGLM::INCREMENTAL : TimeSchemeGLM::NONE;
    shared_ptr<BaseTimeScheme> myScheme(new TimeSchemeGLM(scheme, 0, nlType));
    feFunctions_[MECH_DISPLACEMENT]->SetTimeScheme(myScheme);
  }
  
  void MechPDE::DefinePrimaryResults()
  {
    // Check for subType
    StdVector<std::string> stressDofNames;
    
    if(subType_ == "3d" || subType_ == "2.5d")
      stressDofNames = "xx", "yy", "zz", "yz", "xz", "xy";
    else if(subType_ == "planeStrain")
      stressDofNames = "xx", "yy", "xy";
    else if(subType_ == "planeStress")
      stressDofNames = "xx", "yy", "xy";
    else if(subType_ == "axi")
      stressDofNames = "rr", "zz", "rz", "phiphi";
    else if(subType_ == "2.5d")
      stressDofNames = "xx", "yy", "zz", "yz", "xz", "xy";
    
    // === MECHANIC DISPLACEMENT ===
    shared_ptr<ResultInfo> disp(new ResultInfo);
    disp->resultType = MECH_DISPLACEMENT;
    disp->dofNames = dofNames_;
    disp->unit = "m";
    disp->entryType = ResultInfo::VECTOR;
    disp->SetFeFunction(feFunctions_[MECH_DISPLACEMENT]);
    disp->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_DISPLACEMENT);
    feFunctions_[MECH_DISPLACEMENT]->SetResultInfo(disp);

    // -----------------------------------
    //  Define xml-names of Dirichlet BCs
    // -----------------------------------
    hdbcSolNameMap_[MECH_DISPLACEMENT] = "fix";
    idbcSolNameMap_[MECH_DISPLACEMENT] = "displacement";
    idbcSolNameMapD1_[MECH_DISPLACEMENT] = "velocity";
    idbcSolNameMapD2_[MECH_DISPLACEMENT] = "acceleration";
    
    // this defines the primary unknown
    results_.Push_back( disp );
    
    // define functor for interpolation of the field
    shared_ptr<BaseFeFunction> feFct = feFunctions_[MECH_DISPLACEMENT];
    DefineFieldResult( feFct, disp );
  }
  
  void MechPDE::DefinePostProcResults() {
    Global::ComplexPart part = isComplex_ ? Global::COMPLEX : Global::REAL;
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
    StdVector<std::string > dispDofNames;
    dispDofNames = feFunctions_[MECH_DISPLACEMENT]->GetResultInfo()->dofNames;
    shared_ptr<BaseFeFunction> feFct = feFunctions_[MECH_DISPLACEMENT];
    shared_ptr<BaseFeFunction> vFct;
    shared_ptr<BaseFeFunction> aFct;
    if ( analysistype_ != STATIC && analysistype_ != BUCKLING) {
      // === MECHANIC VELOCITY ===
      // Velocity \bm{v} = \frac{\partial \bm{u}} {\partial t}
      shared_ptr<ResultInfo> vel(new ResultInfo);
      vel->resultType = MECH_VELOCITY;
      vel->dofNames = dispDofNames;
      vel->unit = "m/s";
      vel->entryType = ResultInfo::VECTOR;
      vel->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_VELOCITY);
      availResults_.insert( vel );
      DefineTimeDerivResult( MECH_VELOCITY, 1, MECH_DISPLACEMENT );
      vFct = timeDerivFeFunctions_[MECH_VELOCITY];
      //feFunctions_[MECH_VELOCITY] = vFct;

      // === MECHANIC VELOCITY ELEMRES===
      shared_ptr<ResultInfo> velElem(new ResultInfo);
      velElem->resultType = MECH_VELOCITY_ELEM;
      velElem->dofNames = dispDofNames;
      velElem->unit = "m/s";
      velElem->entryType = ResultInfo::VECTOR;
      velElem->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_VELOCITY_ELEM);
      availResults_.insert( velElem );
      PtrCoefFct velFct= this->GetCoefFct( MECH_VELOCITY );
      PtrCoefFct velFctCoef;
      PtrCoefFct constOne = CoefFunction::Generate( mp_, Global::REAL, "0.0");
      velFctCoef  = CoefFunction::Generate( mp_, Global::REAL, CoefXprBinOp(mp_, constOne, velFct , CoefXpr::OP_ADD ) );
      DefineFieldResult(velFctCoef, velElem);

      // === MECHANIC ACCELERATION ===
      // Acceleration \bm{a} = \frac{\partial^{2} \bm{u}} {\partial t^{2}}
      shared_ptr<ResultInfo> acc(new ResultInfo);
      acc->resultType = MECH_ACCELERATION;
      acc->dofNames = dispDofNames;
      acc->unit = "m/s^2";
      acc->entryType = ResultInfo::VECTOR;
      acc->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_ACCELERATION);
      availResults_.insert( acc );
      DefineTimeDerivResult( MECH_ACCELERATION, 2, MECH_DISPLACEMENT );
      aFct = timeDerivFeFunctions_[MECH_ACCELERATION];
      //feFunctions_[MECH_ACCELERATION] = aFct;
    } // end guard for not static and not buckling
    
    // === MECHANIC RHS ===
    shared_ptr<ResultInfo> rhs(new ResultInfo);
    rhs->resultType = MECH_RHS_LOAD;
    rhs->dofNames = dispDofNames;
    rhs->unit = "N";
    rhs->entryType = ResultInfo::VECTOR;
    rhs->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_RHS_LOAD);
    rhsFeFunctions_[MECH_DISPLACEMENT]->SetResultInfo(rhs);
    DefineFieldResult( rhsFeFunctions_[MECH_DISPLACEMENT], rhs );
    
    // === MECHANIC STRESS ===
    // Stress \bm{\sigma} = \left[ \bm{C} \right] : {\bm{s}}
    shared_ptr<ResultInfo> stress(new ResultInfo);
    stress->resultType = MECH_STRESS;
    stress->dofNames = stressComponents;
    stress->unit =  "N/m^2";
    stress->entryType = ResultInfo::TENSOR;
    stress->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_STRESS);
    PtrCoefFct stressCoef;
    
    shared_ptr<CoefFunctionFormBased> sigmaFunc;
    if( isComplex_ ) {
      sigmaFunc.reset(new CoefFunctionFlux<Complex>(feFct, stress));
    } else {
      sigmaFunc.reset(new CoefFunctionFlux<Double>(feFct, stress));
    }
    // check for thermal strains
    bool isThermalStrain = false;
    PtrParamNode bcsNode = myParam_->Get("bcsAndLoads", ParamNode::PASS );
    if( bcsNode ) {
      if( bcsNode->Has("thermalStrain") ) {
        isThermalStrain = true;
      }
    }
    if ( isThermalStrain )  {
      //! Cauchy stress is [c] ( Bu - alpha DeltaT )
      stressCoef =
              CoefFunction::Generate( mp_, part,
              CoefXprBinOp(mp_,sigmaFunc,thermalStress_,CoefXpr::OP_SUB));
    }
    else {
      stressCoef = sigmaFunc;
    }
    DefineFieldResult( stressCoef, stress );
    stiffFormCoefs_.insert(sigmaFunc);
    
    //If stresses should be cached - take the following lines and go with stressCoefCache
    //shared_ptr<CoefFunctionCache<Double> > stressCoefCache(new CoefFunctionCache<Double>(feFct, stress, stressCoef));
    //DefineFieldResult(stressCoefCache, stress);
    //stiffFormCoefs_.insert(sigmaFunc);
    
    // === MECHANIC PRINCIPAL STRESS ===
    //This result type should not be called in the xml-file. Call e.g. Mech_Principal_Stress_Max instead.
    
    shared_ptr<ResultInfo> principalstress(new ResultInfo);
    
    principalstress->resultType = MECH_PRINCIPAL_STRESS;
    if( subType_ == "3d" || subType_ == "2.5d") {
      principalstress->dofNames = "x", "y", "z", "x", "y", "z", "x", "y", "z", "3", "2", "1";
    } else if( subType_ == "planeStress" || subType_ == "planeStrain" ) {
      principalstress->dofNames = "x", "y", "x", "y", "2", "1";
    }
    principalstress->unit =  "N/m^2";
    principalstress->entryType = ResultInfo::VECTOR;
    principalstress->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_PRINCIPAL_STRESS);
    
    
    //Variante principal stress CACHE, stress NO CACHE
    shared_ptr<CoefFunctionEigen> prinStressCoef(new CoefFunctionEigen(feFct, principalstress, stressCoef));
    shared_ptr<CoefFunctionCache<Double> > prinStressCoefCache(new CoefFunctionCache<Double>(feFct, principalstress, prinStressCoef));
    DefineFieldResult( prinStressCoefCache, principalstress);
    
    // === MECHANIC MINIMUM PRINCIPAL STRESS - VECTOR ===
    shared_ptr<ResultInfo> principalstress_min(new ResultInfo);
    principalstress_min->resultType = MECH_PRINCIPAL_STRESS_MIN;
    principalstress_min->dofNames = dispDofNames;
    principalstress_min->unit = "N/m^2";
    principalstress_min->entryType = ResultInfo::VECTOR;
    principalstress_min->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_PRINCIPAL_STRESS_MIN);
    
    // === MECHANIC MINIMUM PRINCIPAL STRESS - SCALAR ===
    shared_ptr<ResultInfo> principalstress_min_scal(new ResultInfo);
    principalstress_min_scal->resultType = MECH_PRINCIPAL_STRESS_MIN_SCAL;
    principalstress_min_scal->dofNames = "min";
    principalstress_min_scal->unit = "N/m^2";
    principalstress_min_scal->entryType = ResultInfo::SCALAR;
    principalstress_min_scal->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_PRINCIPAL_STRESS_MIN_SCAL);
    
    // === MECHANIC MAXIMUM PRINCIPAL STRESS - VECTOR ===
    shared_ptr<ResultInfo> principalstress_max(new ResultInfo);
    principalstress_max->resultType = MECH_PRINCIPAL_STRESS_MAX;
    principalstress_max->dofNames = dispDofNames;
    principalstress_max->unit = "N/m^2";
    principalstress_max->entryType = ResultInfo::VECTOR;
    principalstress_max->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_PRINCIPAL_STRESS_MAX);
    
    // === MECHANIC MAXIMUM PRINCIPAL STRESS - SCALAR ===
    shared_ptr<ResultInfo> principalstress_max_scal(new ResultInfo);
    principalstress_max_scal->resultType = MECH_PRINCIPAL_STRESS_MAX_SCAL;
    principalstress_max_scal->dofNames = "max";
    principalstress_max_scal->unit = "N/m^2";
    principalstress_max_scal->entryType = ResultInfo::SCALAR;
    principalstress_max_scal->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_PRINCIPAL_STRESS_MAX_SCAL);
    
    // === MECHANIC MEDIUM PRINCIPAL STRESS - VECTOR ===
    shared_ptr<ResultInfo> principalstress_med(new ResultInfo);
    principalstress_med->resultType = MECH_PRINCIPAL_STRESS_MED;
    principalstress_med->dofNames = dispDofNames;
    principalstress_med->unit = "N/m^2";
    principalstress_med->entryType = ResultInfo::VECTOR;
    principalstress_med->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_PRINCIPAL_STRESS_MED);
    
    // === MECHANIC MEDIUM PRINCIPAL STRESS - SCALAR ===
    shared_ptr<ResultInfo> principalstress_med_scal(new ResultInfo);
    principalstress_med_scal->resultType = MECH_PRINCIPAL_STRESS_MED_SCAL;
    principalstress_med_scal->dofNames = "med";
    principalstress_med_scal->unit = "N/m^2";
    principalstress_med_scal->entryType = ResultInfo::SCALAR;
    principalstress_med_scal->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_PRINCIPAL_STRESS_MED_SCAL);
    
    shared_ptr<CoefFunctionCompound<Double> > prinStressCoefMin(new CoefFunctionCompound<Double>(mp_));
    shared_ptr<CoefFunctionCompound<Double> > prinStressCoefMax(new CoefFunctionCompound<Double>(mp_));
    shared_ptr<CoefFunctionCompound<Double> > prinStressCoefMed(new CoefFunctionCompound<Double>(mp_));
    shared_ptr<CoefFunctionCompound<Double> > prinStressCoefMinScal(new CoefFunctionCompound<Double>(mp_));
    shared_ptr<CoefFunctionCompound<Double> > prinStressCoefMaxScal(new CoefFunctionCompound<Double>(mp_));
    shared_ptr<CoefFunctionCompound<Double> > prinStressCoefMedScal(new CoefFunctionCompound<Double>(mp_));
    
    std::map<std::string, PtrCoefFct> var2;
    //Take Coefficients of CoefFctCache
    var2["d"] = prinStressCoefCache;
    
    if ( !isComplex_ ) {
      
      StdVector<std::string> coefMin;
      StdVector<std::string> coefMax;
      StdVector<std::string> coefMed;
      std::string CoefMinScal, CoefMaxScal, CoefMedScal;
      
      if ( stressDim_==3 ) {
        const std::string prinStressCoefMinStr[] = {"d_0_R", "d_1_R"};
        const std::string prinStressCoefMedStr[] = {"0", "0"};
        const std::string prinStressCoefMaxStr[] = {"d_2_R", "d_3_R"};
        CoefMinScal = "d_4_R";
        CoefMedScal = "0";
        CoefMaxScal = "d_5_R";
        coefMin.Import(prinStressCoefMinStr, 2);
        coefMax.Import(prinStressCoefMaxStr, 2);
        coefMed.Import(prinStressCoefMedStr, 2);
      }
      
      else if ( stressDim_==6 ) {
        const std::string prinStressCoefMinStr[] = {"d_0_R", "d_1_R", "d_2_R"};
        const std::string prinStressCoefMedStr[] = {"d_3_R", "d_4_R", "d_5_R"};
        const std::string prinStressCoefMaxStr[] = {"d_6_R", "d_7_R", "d_8_R"};
        CoefMinScal = "d_9_R";
        CoefMedScal = "d_10_R";
        CoefMaxScal = "d_11_R";
        coefMin.Import(prinStressCoefMinStr, 3);
        coefMax.Import(prinStressCoefMaxStr, 3);
        coefMed.Import(prinStressCoefMedStr, 3);
      }
      else if ( stressDim_==4 ) {
        const std::string prinStressCoefMinStr[] = {"0", "0"};
        const std::string prinStressCoefMedStr[] = {"0", "0"};
        const std::string prinStressCoefMaxStr[] = {"0", "0"};
        CoefMinScal = "0";
        CoefMedScal = "0";
        CoefMaxScal = "0";
        coefMin.Import(prinStressCoefMinStr, 2);
        coefMax.Import(prinStressCoefMaxStr, 2);
        coefMed.Import(prinStressCoefMedStr, 2);
        WARN("No implementation for principal stress in axi-formulation yet.");
      }
      
      else {
        EXCEPTION( "Wrong dimension for principal stress: in DefinePostprocResults");
      }
      prinStressCoefMin->SetVector(coefMin, var2);
      prinStressCoefMax->SetVector(coefMax, var2);
      prinStressCoefMed->SetVector(coefMed, var2);
      prinStressCoefMinScal->SetScalar(CoefMinScal, var2);
      prinStressCoefMaxScal->SetScalar(CoefMaxScal, var2);
      prinStressCoefMedScal->SetScalar(CoefMedScal, var2);
      DefineFieldResult( prinStressCoefMin, principalstress_min );
      DefineFieldResult( prinStressCoefMax, principalstress_max );
      DefineFieldResult( prinStressCoefMed, principalstress_med );
      DefineFieldResult( prinStressCoefMinScal, principalstress_min_scal );
      DefineFieldResult( prinStressCoefMaxScal, principalstress_max_scal );
      DefineFieldResult( prinStressCoefMedScal, principalstress_med_scal );
    }
    
    // === THERMOMECHANICAL STRESS ===
    if ( isThermalStrain )  {
      shared_ptr<ResultInfo> stress(new ResultInfo);
      stress->resultType = MECH_THERMAL_STRESS;
      stress->dofNames = stressComponents;
      stress->unit =  "N/m^2";
      stress->entryType = ResultInfo::TENSOR;
      stress->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_THERMAL_STRESS);
      DefineFieldResult( thermalStress_, stress );
    }
    
    // === VON_MISES_STRESS  ===
    // Von Mises Stress \sigma_{\mathrm{v}} = {\sqrt{{\frac{1} {2}}((\sigma_{11} - \sigma_{22})^{2} +
    //                      (\sigma_{22} - \sigma_{33})^{2} + (\sigma_{33} - \sigma_{11})^{2} +
    //                      3(\sigma_{12}^{2}  + \sigma_{23}^{2} + \sigma_{31}^{2})}}    
    shared_ptr<ResultInfo> misesStress(new ResultInfo);
    misesStress->resultType = VON_MISES_STRESS;
    misesStress->dofNames = "";
    misesStress->unit = "N/m^2";
    misesStress->entryType = ResultInfo::SCALAR;
    misesStress->definedOn = ResultInfo::MapSolTypeToDefinedOn(VON_MISES_STRESS);
    if ( !isComplex_ ) {
      shared_ptr<CoefFunctionCompound<Double> > misesCoef(new CoefFunctionCompound<Double>(mp_));
      std::map<std::string, PtrCoefFct> var;
      std::string misesStr;
      
      var["s"] = stressCoef;
      if ( stressDim_==3 ) {
        // 2D plane strain or stress case
        misesStr = "sqrt(   s_0_R^2 + s_1_R^2 \
                          - s_0_R*s_1_R \
                          + 3.0*s_2_R^2 )";
      }
      else if ( stressDim_==4 ) {
        // 2D axi case
        misesStr = "sqrt(   s_0_R^2 + s_1_R^2 + s_3_R^2\
                             - s_0_R*s_1_R - s_0_R*s_3_R - s_1_R*s_3_R  \
                             + 3.0*s_2_R^2 )";
      }
      else if ( stressDim_==6 ) {
        // 3D case
        misesStr = "sqrt(   s_0_R^2 + s_1_R^2 + s_2_R^2\
                           - s_0_R*s_1_R - s_0_R*s_2_R - s_1_R*s_2_R  \
                           + 3.0*(s_3_R^2 + s_4_R^2 + s_5_R^2) )";
      }
      else
        EXCEPTION( "Wrong dimesnion for stress: in DefinePostprocResults");
      
      misesCoef->SetScalar(misesStr,var);
      DefineFieldResult( misesCoef, misesStress );
    }
    
    // === MECHANIC STRAIN ===
    // Strain \bm{s}  = \frac{1} {2}(\nabla \bm{u} + \nabla \bm{u}^{{\mathrm{T}}})
    shared_ptr<ResultInfo> strain(new ResultInfo);
    strain->resultType = MECH_STRAIN;
    strain->dofNames = stressComponents;
    strain->unit =  "";
    strain->entryType = ResultInfo::TENSOR;
    strain->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_STRAIN);
    shared_ptr<CoefFunctionFormBased> strainFunc;
    if( isComplex_ ) {
      strainFunc.reset(new CoefFunctionBOp<Complex>(feFct, strain));
    } else {
      strainFunc.reset(new CoefFunctionBOp<Double>(feFct, strain));
    }
    DefineFieldResult( strainFunc, strain );
    stiffFormCoefs_.insert(strainFunc);
    
    // === MECHANIC PRINCIPAL STRAIN ===
    
    shared_ptr<ResultInfo> principalstrain(new ResultInfo);
    
    principalstrain->resultType = MECH_PRINCIPAL_STRAIN;
    if( subType_ == "3d" || subType_ == "2.5d") {
      principalstrain->dofNames = "x", "y", "z", "x", "y", "z", "x", "y", "z", "3", "2", "1";
    } else if( subType_ == "planeStress" || subType_ == "planeStrain" ) {
      principalstrain->dofNames = "x", "y", "x", "y", "2", "1";
    }
    principalstrain->unit =  "";
    principalstrain->entryType = ResultInfo::VECTOR;
    principalstrain->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_PRINCIPAL_STRAIN);
    
    // === MECHANIC MINIMUM PRINCIPAL STRAIN - VECTOR ===
    shared_ptr<ResultInfo> principalstrain_min(new ResultInfo);
    principalstrain_min->resultType = MECH_PRINCIPAL_STRAIN_MIN;
    principalstrain_min->dofNames = dispDofNames;
    principalstrain_min->unit = "";
    principalstrain_min->entryType = ResultInfo::VECTOR;
    principalstrain_min->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_PRINCIPAL_STRAIN_MIN);
    
    // === MECHANIC MINIMUM PRINCIPAL STRAIN - SCALAR ===
    shared_ptr<ResultInfo> principalstrain_min_scal(new ResultInfo);
    principalstrain_min_scal->resultType = MECH_PRINCIPAL_STRAIN_MIN_SCAL;
    principalstrain_min_scal->dofNames = "";
    principalstrain_min_scal->unit = "";
    principalstrain_min_scal->entryType = ResultInfo::SCALAR;
    principalstrain_min_scal->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_PRINCIPAL_STRAIN_MIN_SCAL);
    
    // === MECHANIC MAXIMUM PRINCIPAL STRAIN - VECTOR ===
    shared_ptr<ResultInfo> principalstrain_max(new ResultInfo);
    principalstrain_max->resultType = MECH_PRINCIPAL_STRAIN_MAX;
    principalstrain_max->dofNames = dispDofNames;
    principalstrain_max->unit = "";
    principalstrain_max->entryType = ResultInfo::VECTOR;
    principalstrain_max->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_PRINCIPAL_STRAIN_MAX);
    
    // === MECHANIC MAXIMUM PRINCIPAL STRAIN - SCALAR ===
    shared_ptr<ResultInfo> principalstrain_max_scal(new ResultInfo);
    principalstrain_max_scal->resultType = MECH_PRINCIPAL_STRAIN_MAX_SCAL;
    principalstrain_max_scal->dofNames = "";
    principalstrain_max_scal->unit = "";
    principalstrain_max_scal->entryType = ResultInfo::SCALAR;
    principalstrain_max_scal->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_PRINCIPAL_STRAIN_MAX_SCAL);
    
    // === MECHANIC MEDIUM PRINCIPAL STRAIN - VECTOR ===
    shared_ptr<ResultInfo> principalstrain_med(new ResultInfo);
    principalstrain_med->resultType = MECH_PRINCIPAL_STRAIN_MED;
    principalstrain_med->dofNames = dispDofNames;
    principalstrain_med->unit = "";
    principalstrain_med->entryType = ResultInfo::VECTOR;
    principalstrain_med->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_PRINCIPAL_STRAIN_MED);
    
    // === MECHANIC MEDIUM PRINCIPAL STRAIN - SCALAR ===
    shared_ptr<ResultInfo> principalstrain_med_scal(new ResultInfo);
    principalstrain_med_scal->resultType = MECH_PRINCIPAL_STRAIN_MED_SCAL;
    principalstrain_med_scal->dofNames = "";
    principalstrain_med_scal->unit = "";
    principalstrain_med_scal->entryType = ResultInfo::SCALAR;
    principalstrain_med_scal->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_PRINCIPAL_STRAIN_MED_SCAL);
    
    shared_ptr<CoefFunctionEigen> prinStrainCoef(new CoefFunctionEigen(feFct, principalstrain, strainFunc));
    shared_ptr<CoefFunctionCache<Double> > prinStrainCoefCache(new CoefFunctionCache<Double>(feFct, principalstrain, prinStrainCoef));
    
    DefineFieldResult( prinStrainCoefCache, principalstrain);
    
    shared_ptr<CoefFunctionCompound<Double> > prinStrainCoefMin(new CoefFunctionCompound<Double>(mp_));
    shared_ptr<CoefFunctionCompound<Double> > prinStrainCoefMax(new CoefFunctionCompound<Double>(mp_));
    shared_ptr<CoefFunctionCompound<Double> > prinStrainCoefMed(new CoefFunctionCompound<Double>(mp_));
    shared_ptr<CoefFunctionCompound<Double> > prinStrainCoefMinScal(new CoefFunctionCompound<Double>(mp_));
    shared_ptr<CoefFunctionCompound<Double> > prinStrainCoefMaxScal(new CoefFunctionCompound<Double>(mp_));
    shared_ptr<CoefFunctionCompound<Double> > prinStrainCoefMedScal(new CoefFunctionCompound<Double>(mp_));
    
    std::map<std::string, PtrCoefFct> var3;
    var3["e"] = prinStrainCoefCache;
    
    if ( !isComplex_ ) {
      
      StdVector<std::string> coefMinReal;
      StdVector<std::string> coefMaxReal;
      StdVector<std::string> coefMedReal;
      std::string CoefMinScal, CoefMaxScal, CoefMedScal;
      
      if ( stressDim_==3 ) {
        const std::string prinStrainCoefMinStr[] = {"e_0_R", "e_1_R"};
        const std::string prinStrainCoefMaxStr[] = {"e_2_R", "e_3_R"};
        const std::string prinStrainCoefMedStr[] = {"0", "0"};
        CoefMinScal = "e_4_R";
        CoefMedScal = "0";
        CoefMaxScal = "e_5_R";
        coefMinReal.Import(prinStrainCoefMinStr, 2);
        coefMaxReal.Import(prinStrainCoefMaxStr, 2);
        coefMedReal.Import(prinStrainCoefMedStr, 2);
      }
      
      else if ( stressDim_==6 ) {
        const std::string prinStrainCoefMinStr[] = {"e_0_R", "e_1_R", "e_2_R"};
        const std::string prinStrainCoefMedStr[] = {"e_3_R", "e_4_R", "e_5_R"};
        const std::string prinStrainCoefMaxStr[] = {"e_6_R", "e_7_R", "e_8_R"};
        CoefMinScal = "e_9_R";
        CoefMedScal = "e_10_R";
        CoefMaxScal = "e_11_R";
        coefMinReal.Import(prinStrainCoefMinStr, 3);
        coefMaxReal.Import(prinStrainCoefMaxStr, 3);
        coefMedReal.Import(prinStrainCoefMedStr, 3);
      }
      
      else if ( stressDim_==4 ) {
        const std::string prinStrainCoefMinStr[] = {"0", "0"};
        const std::string prinStrainCoefMedStr[] = {"0", "0"};
        const std::string prinStrainCoefMaxStr[] = {"0", "0"};
        CoefMinScal = "0";
        CoefMedScal = "0";
        CoefMaxScal = "0";
        coefMinReal.Import(prinStrainCoefMinStr, 2);
        coefMaxReal.Import(prinStrainCoefMaxStr, 2);
        coefMedReal.Import(prinStrainCoefMedStr, 2);
        WARN("No implementation for principal strain in axi-formulation yet.");
      }
      
      else {
        EXCEPTION( "Wrong dimension for principal strain: in DefinePostprocResults");
      }
      prinStrainCoefMin->SetVector(coefMinReal, var3);
      prinStrainCoefMax->SetVector(coefMaxReal, var3);
      prinStrainCoefMed->SetVector(coefMedReal, var3);
      prinStrainCoefMinScal->SetScalar(CoefMinScal, var3);
      prinStrainCoefMaxScal->SetScalar(CoefMaxScal, var3);
      prinStrainCoefMedScal->SetScalar(CoefMedScal, var3);
      DefineFieldResult( prinStrainCoefMin, principalstrain_min );
      DefineFieldResult( prinStrainCoefMax, principalstrain_max );
      DefineFieldResult( prinStrainCoefMed, principalstrain_med );
      DefineFieldResult( prinStrainCoefMinScal, principalstrain_min_scal );
      DefineFieldResult( prinStrainCoefMaxScal, principalstrain_max_scal );
      DefineFieldResult( prinStrainCoefMedScal, principalstrain_med_scal );
    }
    
    if ( isThermalStrain )  {
      shared_ptr<ResultInfo> thermalStrain(new ResultInfo);
      thermalStrain->resultType = MECH_THERMAL_STRAIN;
      thermalStrain->dofNames = stressComponents;
      thermalStrain->unit =  "";
      thermalStrain->entryType = ResultInfo::TENSOR;
      thermalStrain->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_THERMAL_STRAIN);
      DefineFieldResult( thermalStrain_, thermalStrain );
    }
    
    PtrCoefFct intensFct;
    shared_ptr<CoefFunctionFormBased> kedFunc;
    shared_ptr<ResultFunctor> keFunc;
    shared_ptr<CoefFunctionSurf> sNormStructIntens;
    shared_ptr<ResultFunctor> powerFunc;
    if ( analysistype_ != STATIC && analysistype_ != BUCKLING) {
      // === MECHANIC STRUCTURAL INTENSTIY ===
      shared_ptr<ResultInfo> intens(new ResultInfo);
      intens->resultType = MECH_STRUCT_INTENSTIY;
      intens->dofNames = dispDofNames ;
      intens->unit =  "N/ms";
      intens->entryType = ResultInfo::VECTOR;
      intens->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_STRUCT_INTENSTIY);
      
      // The mechanic structural intensity calculates as
      // I = -[stress] * conj(v)
      PtrCoefFct velFnc = this->GetCoefFct( MECH_VELOCITY );
      // define temporary function, without the -1 sign
      PtrCoefFct intensTmp = CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_, sigmaFunc, velFnc, CoefXpr::OP_MULT_VOIGT_TENSOR_VEC_CONJ));
      intensFct = CoefFunction::Generate(mp_, part, CoefXprBinOp(mp_,  "-1.0", intensTmp , CoefXpr::OP_MULT));
      DefineFieldResult( intensFct, intens );
      
      
      // === MECHANIC KINETIC ENERGY DENSITY ===
      // kinetic Energy Density e_{\mathrm{kin}} = {\frac{1} {2}} \rho \bm{v} \cdot \bm{v}
      shared_ptr<ResultInfo> kinEnergyDens(new ResultInfo);
      kinEnergyDens->resultType = MECH_KIN_ENERGY_DENS;
      kinEnergyDens->dofNames = "";
      kinEnergyDens->unit = "Ws/m^3";
      kinEnergyDens->entryType = ResultInfo::SCALAR;
      kinEnergyDens->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_KIN_ENERGY_DENS);
      if( isComplex_ ) {
        kedFunc.reset(new CoefFunctionBdBKernel<Complex>(vFct, 0.5));
      } else {
        kedFunc.reset(new CoefFunctionBdBKernel<Double>(vFct, 0.5));
      }
      DefineFieldResult( kedFunc, kinEnergyDens );
      massFormCoefs_.insert(kedFunc);
      
      // === MECHANIC KINETIC ENERGY ===
      // kinetic Energy E_{\mathrm{kin}} = \int_{\Omega} e_{\mathrm{kin}} \mathrm{d} \Omega
      shared_ptr<ResultInfo> kinEnergy(new ResultInfo);
      kinEnergy->resultType = MECH_KIN_ENERGY;
      kinEnergy->dofNames = "";
      kinEnergy->unit = "Ws";
      kinEnergy->entryType = ResultInfo::SCALAR;
      kinEnergy->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_KIN_ENERGY);
      availResults_.insert ( kinEnergy );
      if( isComplex_ ) {
        keFunc.reset(new EnergyResultFunctor<Complex>(vFct, kinEnergy, 0.5));
      } else {
        keFunc.reset(new EnergyResultFunctor<Double>(vFct, kinEnergy, 0.5));
      }
      resultFunctors_[MECH_KIN_ENERGY] = keFunc;
      massFormFunctors_.insert(keFunc);
      
      // === MECHANIC NORMAL STRUCTURAL INTENSTIY ===
      shared_ptr<ResultInfo> intensNormal(new ResultInfo);
      intensNormal->resultType = MECH_NORMAL_STRUCT_INTENSITY;
      intensNormal->dofNames = "" ;
      intensNormal->unit =  "N/ms";
      intensNormal->entryType = ResultInfo::SCALAR;
      intensNormal->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_NORMAL_STRUCT_INTENSITY);
      sNormStructIntens.reset(new CoefFunctionSurf(true, 1.0, intensNormal));
      DefineFieldResult( sNormStructIntens, intensNormal );
      surfCoefFcts_[sNormStructIntens] = intensFct;
      
      // === MECHANIC POWER ===
      shared_ptr<ResultInfo> power(new ResultInfo);
      power->resultType = MECH_POWER;
      power->dofNames = "";
      power->unit = "W";
      power->entryType = ResultInfo::SCALAR;
      power->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_POWER);
      availResults_.insert( power );
      // then, integrate values
      if(isComplex_)
        powerFunc.reset(new ResultFunctorIntegrate<Complex>(sNormStructIntens, feFct, power));
      else
        powerFunc.reset(new ResultFunctorIntegrate<Double>(sNormStructIntens, feFct, power));
      
      resultFunctors_[MECH_POWER] = powerFunc;
    }
    
    // === scalar product of displacement or quadratic displacement norm. For topology gradient for Bloch mode analysis (Nazarov) ===
    shared_ptr<ResultInfo> quadDisp(new ResultInfo);
    quadDisp->resultType = MECH_QUAD_DISP;
    quadDisp->dofNames = "";
    quadDisp->unit = "m^2";
    quadDisp->entryType = ResultInfo::SCALAR;
    quadDisp->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_QUAD_DISP);
    shared_ptr<CoefFunctionFormBased> quad;
    if(isComplex_)
      quad.reset(new CoefFunctionQuadSol<Complex>(feFct));
    else
      quad.reset(new CoefFunctionQuadSol<double>(feFct));
    DefineFieldResult(quad, quadDisp);
    
    // === integrated quadratic displacement ===
    shared_ptr<ResultInfo> quadDispSum(new ResultInfo);
    quadDispSum->resultType = MECH_QUAD_DISP_SUM;
    quadDispSum->dofNames = "";
    quadDispSum->unit = "m^2";
    quadDispSum->entryType = ResultInfo::SCALAR;
    quadDispSum->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_QUAD_DISP_SUM);
    availResults_.insert( quadDispSum );
    shared_ptr<ResultFunctor> qdsFunc;
    if(isComplex_)
      qdsFunc.reset(new ResultFunctorIntegrate<Complex>(quad, feFct, quadDispSum));
    else
      qdsFunc.reset(new ResultFunctorIntegrate<Double>(quad, feFct, quadDispSum));
    resultFunctors_[MECH_QUAD_DISP_SUM] = qdsFunc;
    
    // === Hermitian dyadic strain prodcut for topology gradient for bloch mode analysis (Nazarov) ===
    shared_ptr<ResultInfo> dyadicStrain(new ResultInfo);
    dyadicStrain->resultType = MECH_DYADIC_STRAIN;
    dyadicStrain->dofNames = "e11", "e12", "e13", "e21", "e22", "e23", "e31", "e32", "e33";
    dyadicStrain->unit = "m/m";
    dyadicStrain->entryType = ResultInfo::TENSOR;
    dyadicStrain->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_DYADIC_STRAIN);
    shared_ptr<CoefFunctionFormBased> dyadic;
    if(isComplex_)
      dyadic.reset(new CoefFunctionDyadicStrain<Complex>(feFct));
    else
      dyadic.reset(new CoefFunctionDyadicStrain<double>(feFct));
    DefineFieldResult(dyadic, dyadicStrain);
    stiffFormCoefs_.insert(dyadic);
    
    // === integrated MECH_DYADIC_STRAIN ===
    shared_ptr<ResultInfo> dyadicStrainSum(new ResultInfo);
    dyadicStrainSum->resultType = MECH_DYADIC_STRAIN_SUM;
    dyadicStrainSum->dofNames = "e11", "e12", "e13", "e21", "e22", "e23", "e31", "e32", "e33";
    dyadicStrainSum->unit = "m/m";
    dyadicStrainSum->entryType = ResultInfo::TENSOR;
    dyadicStrainSum->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_DYADIC_STRAIN_SUM);
    availResults_.insert( dyadicStrainSum );
    shared_ptr<ResultFunctor> dssFunc;
    if(isComplex_)
      dssFunc.reset(new ResultFunctorIntegrate<Complex>(dyadic, feFct, dyadicStrainSum));
    else
      dssFunc.reset(new ResultFunctorIntegrate<Double>(dyadic, feFct, dyadicStrainSum));
    resultFunctors_[MECH_DYADIC_STRAIN_SUM] = dssFunc;

    // === MECHANIC DEFORMATION ENERGY DENSITY ===
    shared_ptr<ResultInfo> defEnergyDens(new ResultInfo);
    defEnergyDens->resultType = MECH_DEFORM_ENERGY_DENS;
    defEnergyDens->dofNames = "";
    defEnergyDens->unit = "Ws/m^3";
    defEnergyDens->entryType = ResultInfo::SCALAR;
    defEnergyDens->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_DEFORM_ENERGY_DENS);
    shared_ptr<CoefFunctionFormBased> dedFunc;
    if( isComplex_ ) {
      dedFunc.reset(new CoefFunctionBdBKernel<Complex>(feFct, 0.5));
    } else {
      dedFunc.reset(new CoefFunctionBdBKernel<Double>(feFct, 0.5));
    }
    DefineFieldResult( dedFunc, defEnergyDens );
    stiffFormCoefs_.insert(dedFunc);
    
    // === MECHANIC DEFORMATION ENERGY ===
    // deformation Energy E_{\mathrm{def}} = {\frac{1} {2}} <\bm{u},\bm{K u}>
    shared_ptr<ResultInfo> defEnergy(new ResultInfo);
    defEnergy->resultType = MECH_DEFORM_ENERGY;
    defEnergy->dofNames = "";
    defEnergy->unit = "Ws";
    defEnergy->entryType = ResultInfo::SCALAR;
    defEnergy->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_DEFORM_ENERGY);
    availResults_.insert( defEnergy );
    shared_ptr<ResultFunctor> deFunc;
    if(isComplex_)
      deFunc.reset(new EnergyResultFunctor<Complex>(feFct, defEnergy, 0.5));
    else
      deFunc.reset(new EnergyResultFunctor<Double>(feFct, defEnergy, 0.5));
    
    resultFunctors_[MECH_DEFORM_ENERGY] = deFunc;
    stiffFormFunctors_.insert(deFunc);
    
    // === MECHANIC TOTAL ENERGY DENSITY ===
    shared_ptr<ResultInfo> totEnergyDens(new ResultInfo);
    totEnergyDens->resultType = MECH_TOTAL_ENERGY_DENS;
    totEnergyDens->dofNames = "";
    totEnergyDens->unit = "Ws";
    totEnergyDens->entryType = ResultInfo::SCALAR;
    totEnergyDens->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_TOTAL_ENERGY_DENS);
    shared_ptr<CoefFunction> tedFunc;
    // in static and buckling analysis, the total energy density equals the deformation one
    if (analysistype_ == STATIC  || analysistype_ == BUCKLING)
      tedFunc = dedFunc;
    else
      tedFunc = CoefFunction::Generate(mp_, part, CoefXprBinOp( mp_, dedFunc, kedFunc, CoefXpr::OP_ADD) );
    DefineFieldResult(tedFunc, totEnergyDens);
    
    // === MECHANIC TOTALENERGY ===
    shared_ptr<ResultInfo> tEnergy(new ResultInfo);
    tEnergy->resultType = MECH_TOTAL_ENERGY;
    tEnergy->dofNames = "";
    tEnergy->unit = "Ws";
    tEnergy->entryType = ResultInfo::SCALAR;
    tEnergy->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_TOTAL_ENERGY);
    availResults_.insert( tEnergy );
    shared_ptr<ResultFunctor> teFunc;
    if( isComplex_ ) {
      teFunc.reset(new ResultFunctorIntegrate<Complex>(tedFunc, feFct, tEnergy ) );
    } else {
      teFunc.reset(new ResultFunctorIntegrate<Double>(tedFunc, feFct, tEnergy ) );
    }
    resultFunctors_[MECH_TOTAL_ENERGY] = teFunc;
    
    
    // === MECHANIC DISPLACED SURFACE VOLUME ===
    shared_ptr<ResultInfo> dispNormal, dispVol;
    shared_ptr<CoefFunctionSurf> dispFctNormal;
    
    //normal mechanical displacement
    dispNormal.reset(new ResultInfo);
    dispNormal->resultType = MECH_NORMAL_DISPLACEMENT;
    dispNormal->dofNames = "";
    dispNormal->unit = "m";
    dispNormal->entryType = ResultInfo::SCALAR;
    dispNormal->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_NORMAL_DISPLACEMENT);
    
    dispFctNormal.reset(new CoefFunctionSurf(true, 1.0, dispNormal));
    DefineFieldResult(dispFctNormal, dispNormal);
    surfCoefFcts_[dispFctNormal] = feFct;
    
    dispVol.reset(new ResultInfo);
    dispVol->resultType = MECH_DEF_SURF_VOLUME;
    dispVol->dofNames = "";
    dispVol->unit = "m^3";
    dispVol->entryType = ResultInfo::SCALAR;
    dispVol->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_DEF_SURF_VOLUME);
    // Integrate normal displacement
    shared_ptr<ResultFunctor> dispVolFct;
    if(isComplex_)
      dispVolFct.reset(new ResultFunctorIntegrate<Complex>(dispFctNormal, feFct, dispVol));
    else
      dispVolFct.reset(new ResultFunctorIntegrate<Double>(dispFctNormal, feFct, dispVol));
    resultFunctors_[MECH_DEF_SURF_VOLUME] = dispVolFct;
    availResults_.insert(dispVol);
    
    // === MECHANIC_NORMAL_STRESS ===
    shared_ptr<ResultInfo> normalStressInfo;
    shared_ptr<CoefFunctionSurf> normalStressFct;
    normalStressInfo.reset(new ResultInfo);
    normalStressInfo->resultType = MECH_NORMAL_STRESS;
    normalStressInfo->dofNames = dispDofNames;
    normalStressInfo->unit = "Pa";
    normalStressInfo->entryType = ResultInfo::VECTOR;
    normalStressInfo->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_NORMAL_STRESS);
    
    normalStressFct.reset(new CoefFunctionSurf(true, 1.0, normalStressInfo));
    DefineFieldResult(normalStressFct, normalStressInfo);
    surfCoefFcts_[normalStressFct] = sigmaFunc;
    
    // === MECHANIC REACTION FORCE (= integral of surface traction, i.e. normal stress from above, over the surface region ) ===
    shared_ptr<ResultInfo> reactionForceInfo;
    reactionForceInfo.reset(new ResultInfo);
    reactionForceInfo->resultType = MECH_FORCE;
    reactionForceInfo->dofNames = dispDofNames;
    reactionForceInfo->unit = "N";
    reactionForceInfo->entryType = ResultInfo::VECTOR;
    reactionForceInfo->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_FORCE);
    // Integrate surface traction
    shared_ptr<ResultFunctor> reactionForceFct;
    if(isComplex_)
        reactionForceFct.reset(new ResultFunctorIntegrate<Complex>(normalStressFct, feFct, reactionForceInfo));
    else
        reactionForceFct.reset(new ResultFunctorIntegrate<Double>(normalStressFct, feFct, reactionForceInfo));
    resultFunctors_[MECH_FORCE] = reactionForceFct;
    availResults_.insert(reactionForceInfo);

    // === MECH_TENSOR_HILL_MANDEL converts to HillMandel notation
    shared_ptr<ResultInfo> mech_tensor_hm(new ResultInfo);
    mech_tensor_hm->resultType = MECH_TENSOR_HILL_MANDEL;
    if(dim_ == 2)
      mech_tensor_hm->dofNames = "e11", "e22", "e33", "e23", "e13", "e12";
    else
      mech_tensor_hm->dofNames = "e11", "e22", "e33", "e44", "e55", "e66", "e56", "e46", "e36", "e26", "e16", "e45", "e35", "e25", "e15", "e34", "e24", "e14", "e23", "e13", "e12";

    mech_tensor_hm->unit = "Pa";
    mech_tensor_hm->entryType = ResultInfo::TENSOR;
    mech_tensor_hm->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_TENSOR_HILL_MANDEL);
    shared_ptr<CoefFunctionFormBased> stiff_coef_hm;
    if(isComplex_) // does not really handle the case where only some regions have complex material
      stiff_coef_hm.reset(new CoefFunctionHomogenization<Complex, App::MECH>(feFct, HILL_MANDEL));
    else
      stiff_coef_hm.reset(new CoefFunctionHomogenization<double, App::MECH>(feFct, HILL_MANDEL));
    DefineFieldResult(stiff_coef_hm, mech_tensor_hm);
    stiffFormCoefs_.insert(stiff_coef_hm); // will define the forms
    
    // === MECH_TENSOR for free and parameterized material optimization but generally it simply returns the tensor
    shared_ptr<ResultInfo> mech_tensor(new ResultInfo);
    mech_tensor->resultType = MECH_TENSOR;
    if(dim_ == 2)
      mech_tensor->dofNames = "e11", "e22", "e33", "e23", "e13", "e12";
    else
      mech_tensor->dofNames = "e11", "e22", "e33", "e44", "e55", "e66", "e56", "e46", "e36", "e26", "e16", "e45", "e35", "e25", "e15", "e34", "e24", "e14", "e23", "e13", "e12";

    mech_tensor->unit = "Pa";
    mech_tensor->entryType = ResultInfo::TENSOR;
    mech_tensor->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_TENSOR);
    shared_ptr<CoefFunctionFormBased> stiff_coef;
    if(isComplex_) // does not really handle the case where only some regions have complex material
      stiff_coef.reset(new CoefFunctionHomogenization<Complex, App::MECH>(feFct, VOIGT));
    else
      stiff_coef.reset(new CoefFunctionHomogenization<double, App::MECH>(feFct, VOIGT));
    DefineFieldResult(stiff_coef, mech_tensor);
    stiffFormCoefs_.insert(stiff_coef); // will define the forms
    
    // optimization results are provided in DesignSpace::ExtractResults()
    
    // === MECH_PSEUDO_DENISTY ===
    shared_ptr<ResultInfo> mpd(new ResultInfo);
    mpd->resultType = MECH_PSEUDO_DENSITY;
    mpd->entryType = ResultInfo::SCALAR;
    mpd->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_PSEUDO_DENSITY);
    mpd->dofNames = "";
    mpd->fromOptimization = true;
    DefineFieldResult(shared_ptr<FeFunction<double> >(new FeFunction<double>(NULL)), mpd); // the fe-function is only a dummy
    
    // === PHYSICAL_PSEUDO_DENISTY ===
    shared_ptr<ResultInfo> ppd(new ResultInfo);
    ppd->resultType = PHYSICAL_PSEUDO_DENSITY;
    ppd->entryType = ResultInfo::SCALAR;
    ppd->definedOn = ResultInfo::MapSolTypeToDefinedOn(PHYSICAL_PSEUDO_DENSITY);
    ppd->dofNames = "";
    ppd->fromOptimization = true;
    DefineFieldResult(shared_ptr<FeFunction<double> >(new FeFunction<double>(NULL)), ppd);
    
    // === MECH_TENSOR_TRACE for free and parameterized material optimization
    shared_ptr<ResultInfo> mtt(new ResultInfo);
    mtt->resultType = MECH_TENSOR_TRACE;
    mtt->dofNames = "Tr(E)";
    mtt->unit = "Pa";
    mtt->entryType = ResultInfo::SCALAR;
    mtt->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_TENSOR_TRACE);
    mtt->fromOptimization = true;
    DefineFieldResult(shared_ptr<FeFunction<double> >(new FeFunction<double>(NULL)), mtt);
    
    // === MECH_SHAPE for ms optimization ===
    shared_ptr<ResultInfo> ms(new ResultInfo);
    ms->resultType = MECH_SHAPE;
    ms->dofNames = dispDofNames;
    ms->unit = "m";
    ms->entryType = ResultInfo::VECTOR;
    ms->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_SHAPE);
    ms->fromOptimization = true;
    DefineFieldResult(shared_ptr<FeFunction<double> >(new FeFunction<double>(NULL)), ms);
    
    // === MECH_ELEM_VOL for free and parameterized material optimization
    shared_ptr<ResultInfo> mev(new ResultInfo);
    mev->resultType = MECH_ELEM_VOL;
    mev->dofNames = "elemVol";
    mev->unit = "";
    mev->entryType = ResultInfo::SCALAR;
    mev->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_ELEM_VOL);
    mev->fromOptimization = true;
    DefineFieldResult(shared_ptr<FeFunction<double> >(new FeFunction<double>(NULL)), mev);

    // === MECH_ELEM_POROSITY for parameterized material optimization (currently only necessary for two-scale optimization)
    shared_ptr<ResultInfo> mep(new ResultInfo);
    mep->resultType = MECH_ELEM_POROSITY;
    mep->dofNames = "elemPorosity";
    mep->unit = "";
    mep->entryType = ResultInfo::SCALAR;
    mep->definedOn = ResultInfo::MapSolTypeToDefinedOn(MECH_ELEM_POROSITY);
    mep->fromOptimization = true;
    DefineFieldResult(shared_ptr<FeFunction<double> >(new FeFunction<double>(NULL)), mep);

    // the OPT_RESULT_* are added via the optimization stuff in DesignSpace.
  }
  
  std::map<SolutionType, shared_ptr<FeSpace> > MechPDE::CreateFeSpaces(const std::string& formulation, PtrParamNode infoNode)
  {
    std::map<SolutionType, shared_ptr<FeSpace> > crSpaces;
    
    if( formulation == "default" || formulation == "H1" )
    {
      PtrParamNode potSpaceNode = infoNode->Get("mechDisplacement");
      crSpaces[MECH_DISPLACEMENT] = FeSpace::CreateInstance(myParam_,potSpaceNode,FeSpace::H1, ptGrid_);
      crSpaces[MECH_DISPLACEMENT]->Init(solStrat_);
    }
    else
      EXCEPTION( "The formulation " << formulation << "of the mechanic PDE is not known!" );
    return crSpaces;
  }
  
  //+++++++++++++++++++++++++++++++++++++++++++++++++
  // PreStressing CoefFunction Creation
  
  //Helper
  void MakeBigPreStressVector(StdVector<std::string>& bigVec, StdVector<std::string> smallVec, UInt dim) {
    // this is in accordance with Voigt notation
    if(dim ==2){
      bigVec.Resize(16);
      bigVec.Init("0.0");
      bigVec[0] = smallVec[0];
      bigVec[1] = smallVec[2];
      bigVec[4] = smallVec[2];
      bigVec[5] = smallVec[1];
      bigVec[10] = smallVec[0];
      bigVec[11] = smallVec[2];
      bigVec[14] = smallVec[2];
      bigVec[15] = smallVec[1];
    }else if (dim ==3){
      bigVec.Resize(81);
      bigVec.Init("0.0");
      for(UInt i=0;i<3;i++){
        bigVec[i*30+0] = smallVec[0];
        bigVec[i*30+1] = smallVec[5];
        bigVec[i*30+2] = smallVec[4];
        
        bigVec[i*30+9] = smallVec[5];
        bigVec[i*30+10] = smallVec[1];
        bigVec[i*30+11] = smallVec[3];
        
        bigVec[i*30+18] = smallVec[4];
        bigVec[i*30+19] = smallVec[3];
        bigVec[i*30+20] = smallVec[2];
      }
    }
  }

  PtrCoefFct MechPDE::CreatePreStressFct( bool isComplex, PtrParamNode stressNode){

    PtrCoefFct   coef;
    PtrParamNode inputNode;
    
    UInt dimPre = dim_;
    // in 2.5D case the dimension of the prestress vector is the same as that in 3D, i.e. = 3
    if (subType_ == "2.5d")
      dimPre = 3;
    
    if(stressNode != NULL && stressNode->Has("prescribedLHS")){

      inputNode = stressNode->Get("prescribedLHS",ParamNode::PASS);
      //TODO: This does not support coordinate systems. If this is needed,
      // one possibility would be to create a stress tensor coeffunction first, apply coordinate systems and then blow it up
      // according to the space dimension
      
      //execute strTok and cast
      typedef boost::tokenizer<boost::char_separator<char> >
      tokenizer;
      boost::char_separator<char> sep(" ");

      StdVector<std::string> preVecR;
      StdVector<std::string> preVecI;
      
      std::string valueRStr =  inputNode->Get("value")->As<std::string>();
      tokenizer tokensR(valueRStr, sep);
      for (tokenizer::iterator tok_iter = tokensR.begin();
              tok_iter != tokensR.end(); ++tok_iter){
        preVecR.Push_back(*tok_iter);
      }
      
      //some consistency checks
      if(dimPre == 2 && preVecR.GetSize() != 3){
        Exception("For a 2D simulation, we expect 3 values for real preStress tensor in Voigt notation.");
      }
      if(dimPre == 3 && preVecR.GetSize() != 6){
        Exception("For a 2.5D or a 3D simulation, we expect 6 values for real preStress tensor in Voigt notation.");
      }
      
      //first we create the big tensor for real values
      StdVector<std::string> bigVecR;
      StdVector<std::string> bigVecI;
      MakeBigPreStressVector(bigVecR,preVecR,dimPre);
      
      if(isComplex){
        //create just an empty tensor
        preVecI.Resize(preVecR.GetSize());
        preVecI.Init("0.0");
        MakeBigPreStressVector(bigVecI,preVecI,dimPre);
      }
      
      if(isComplex){
        coef =  CoefFunction::Generate(mp_,Global::COMPLEX,dimPre*dimPre,dimPre*dimPre,bigVecR,bigVecI);
      }else{
        coef =  CoefFunction::Generate(mp_,Global::REAL,dimPre*dimPre,dimPre*dimPre,bigVecR);
      }
      return coef;

    }
    else if(stressNode == NULL || stressNode->Has("computeLHS") || stressNode->Has("referenceStress")){

      // only real valued coefFunctions supported!
      PtrCoefFct stressVec;

      if(stressNode == NULL) {
        // dummy for PostInit
        // we update this with the actual values in each iteration by Excitation::SetStressCoefFct
        stressVec.reset(new CoefFunctionConst<Double>());
        Vector<Double> vec(dimPre == 2 ? 3 : 6);
        vec.Init(0.0);
        dynamic_cast< CoefFunctionConst<Double>* > (stressVec.get())->SetVector(vec);
      }
      else {
        if (stressNode->Has("computeLHS"))
          inputNode = stressNode->Get("computeLHS",ParamNode::PASS);
        else if(stressNode->Has("referenceStress"))
          inputNode = stressNode->Get("referenceStress",ParamNode::PASS);

        UInt aSStep = 0;

        //redefine if user passes the argument
        if( inputNode->Get("sequenceStep",ParamNode::PASS) )
          aSStep = inputNode->Get("sequenceStep",ParamNode::PASS)->As<UInt>();

        if(aSStep < 1) {
          // GetPreceeding sequence step
          aSStep = domain_->GetDriver()->GetActSequenceStep();
          aSStep--;
        }

        stressVec = GetStressCoefFromSeqStep(aSStep);
      }

      std::map<std::string, PtrCoefFct> var;
      var["a"]  = stressVec;

      StdVector<std::string> preVecR;
      StdVector<std::string> preVecI;
      StdVector<std::string> bigVecR;
      StdVector<std::string> bigVecI;

      if(dimPre == 2){
        const std::string vecR[] = { "a_0_R" , "a_1_R" , "a_2_R" };
        preVecR.Import(vecR,3);
        // convert vector preVecR to 4th-order tensor with 2^4 entries
        MakeBigPreStressVector(bigVecR,preVecR,dimPre);
        if(isComplex){
          const std::string vecI[] = { "0.0" , "0.0" , "0.0" };
          preVecI.Import(vecI,3);
          MakeBigPreStressVector(bigVecI,preVecI,dimPre);
        }
      } else {
        const std::string vecR[] = { "a_0_R" , "a_1_R" , "a_2_R" , "a_3_R" , "a_4_R" , "a_5_R"};
        preVecR.Import(vecR,6);
        // convert vector preVecR to 4th-order tensor with 3^4 entries
        MakeBigPreStressVector(bigVecR,preVecR,dimPre);
        if(isComplex){
          const std::string vecI[] = { "0.0" , "0.0" , "0.0" , "0.0" , "0.0" , "0.0" };
          preVecI.Import(vecI,6);
          MakeBigPreStressVector(bigVecI,preVecI,dimPre);
        }
      }

      //create the coefFunction object
      if(isComplex){
        coef.reset(new CoefFunctionCompound<Complex>(mp_));
        CoefFunctionCompound<Complex>*  stressTens = dynamic_cast< CoefFunctionCompound<Complex>* > (coef.get());
        stressTens->SetTensor(bigVecR,bigVecI,dimPre*dimPre,dimPre*dimPre,var);
      } else {
        coef.reset(new CoefFunctionCompound<Double>(mp_));
        CoefFunctionCompound<Double>*  stressTens = dynamic_cast< CoefFunctionCompound<Double>* > (coef.get());
        stressTens->SetTensor(bigVecR,dimPre*dimPre,dimPre*dimPre,var);
      }

      if(domain->HasDesign())
      {
        CoefFunctionOpt* tmpFnc = new CoefFunctionOpt(domain->GetDesign(), coef, NO_MATERIAL, this); // pre stress is not really a material type
        coef.reset(tmpFnc);
      }

      return coef;
    } else {
      EXCEPTION("Cannot read definition of prestressing!");
      return coef;
    }
  }

  PtrCoefFct MechPDE::GetStressCoefFromSeqStep(UInt seqStep){
    //This function uses mostly the simState algorithms from SinglePDE
    //TODO: This is the third(?) time this code is used (see ReadUserFieldValues and ReadInitialConditions).
    //Define some function to handle feFunction extraction from other sequence steps or external simulations
    
    Domain * inDomain = NULL;
    PtrCoefFct stressVec;
    //Get Stress CoefFunction from previous state
    shared_ptr<SimState> inState(new SimState(true, domain_));
    
    PtrParamNode icInfo = infoNode_->Get("Prestressing");
    PtrParamNode isInfo = icInfo->Get("ComputedLHS");
    try{
      std::string fileName = simState_->GetOutputWriter()->GetFileName().string();
      PtrParamNode node(new ParamNode());
      PtrParamNode infoNode = ParamNode::GenerateWriteNode("", "", ParamNode::APPEND); // empty filename means we don't write and ignore ParamNode::ToFile()
      shared_ptr<SimInputHDF5> in;
      in.reset(new SimInputHDF5(fileName, node, infoNode));
      inState->SetInputHdf5Reader(in);
      SimState::GridMap gridMap = domain_->GetGridMap();


      inDomain = inState->GetDomain(seqStep, gridMap);
      Double stepVal = 0.0;
      UInt lastStepNum = 0;
      inState->GetLastStepNum(seqStep, lastStepNum, stepVal);
      // log to info node
      isInfo->Get("inputSequenceStep")->SetValue(seqStep);
      isInfo->Get("inputStepNumber")->SetValue(lastStepNum);
      // update to last step number
      inState->SetInterpolation(SimState::CONSTANT, mp_, analysistype_, 0);
      inState->UpdateToStep(seqStep, lastStepNum);
      SinglePDE * inPDE = inDomain->GetSinglePDE(pdename_);
      if( inPDE->GetAnalysisType() != STATIC && inPDE->GetAnalysisType() != TRANSIENT){
        EXCEPTION("Prestressing is only supported for a preceding transient or static analysis");
      }
      
      // Directly acquire the mechanical stress from the previous sequence step
      stressVec = inPDE->GetCoefFct(MECH_STRESS);
      
      // Store the data input for later. It will be destroyed in the destructor of the SinglePDE
      inputs_[inState] = inDomain;
    } catch (Exception& e) {
      if( inState ) {
        inState->Finalize();
        inState.reset();
      }
      if(inDomain)
        delete inDomain;
      
      RETHROW_EXCEPTION(e, "Cannot obtain mechanic Stress coefficient function from last sequence step for prestressing."
              << "' from sequenceStep " << seqStep );
    }
    return stressVec;
  }
  
  MechPDE::CoefFunction2ndPiolaTensor::
  CoefFunction2ndPiolaTensor(SubTensorType &subType,
          PtrCoefFct stiffness,
          shared_ptr<BaseFeFunction> displ)
  : CoefFunction(),
          tensorType_(subType),
          stiffCoef_(stiffness),
          linOp_(NULL),
          nonLinOp_(NULL)
  {
    switch (tensorType_) {
      case FULL:
        linOp_ = new StrainOperator3D<FeH1,Double>(false);
        nonLinOp_ = new NonLinStrainOperator3D<FeH1,Double>(displ);
        break;
      case AXI:
        linOp_ = new StrainOperatorAxi<FeH1,Double>(false);
        nonLinOp_ = new NonLinStrainOperatorAxi<FeH1,Double>(displ);
        break;
      case PLANE_STRAIN:
      case PLANE_STRESS:
        linOp_ = new StrainOperator2D<FeH1,Double>(false);
        nonLinOp_ = new NonLinStrainOperator2D<FeH1,Double>(displ);
        break;
      default:
        EXCEPTION("Unknown mechanical subtype: " << tensorType_ << std::endl);
    }
    
    if (displ->IsComplex()) {
      dispCoefComplex_ = dynamic_pointer_cast< FeFunction<Complex> >(displ);
    }
    else {
      dispCoefReal_ = dynamic_pointer_cast< FeFunction<Double> >(displ);
    }
    if (!dispCoefComplex_ && !dispCoefReal_) {
      EXCEPTION("Could not cast BaseFeFunction to FeFunction!");
    }
    
    dimType_ = TENSOR;
    dependType_ = SOLUTION;
    isAnalytic_ = false;
    isComplex_ = false;
  }
  
  MechPDE::CoefFunction2ndPiolaTensor::~CoefFunction2ndPiolaTensor() {
    if (linOp_)
      delete linOp_;
    if (nonLinOp_)
      delete nonLinOp_;
  }
  
  void MechPDE::CoefFunction2ndPiolaTensor::
  GetTensorSize( UInt& numRows, UInt& numCols ) const {
    switch (tensorType_) {
      case FULL:
        numRows = StrainOperator3D<FeH1,Double>::DIM_D_MAT;
        numCols = StrainOperator3D<FeH1,Double>::DIM_D_MAT;
        break;
      case AXI:
        numRows = StrainOperatorAxi<FeH1,Double>::DIM_D_MAT;
        numCols = StrainOperatorAxi<FeH1,Double>::DIM_D_MAT;
        break;
      case PLANE_STRAIN:
      case PLANE_STRESS:
        numRows = StrainOperator2D<FeH1,Double>::DIM_D_MAT;
        numCols = StrainOperator2D<FeH1,Double>::DIM_D_MAT;
        break;
      default:
        EXCEPTION("Unknown mechanical subtype: " << tensorType_ << std::endl);
    }
  }
  
  void MechPDE::CoefFunction2ndPiolaTensor::
  GetVector(Vector<Double>& vec, const LocPointMapped& lpm ) {
    Vector<Double> disp, linStrain, nlStrain, totalStrain;
    Matrix<Double> stiff;
    BaseFE *ptFe = dispCoefReal_->GetFeSpace()->GetFe(lpm.ptEl->elemNum);
    
    dispCoefReal_->GetElemSolution(disp, lpm.ptEl);
    stiffCoef_->GetTensor(stiff, lpm);
    
    linOp_->ApplyOp(linStrain, lpm, ptFe, disp);
    
    nonLinOp_->ApplyOp(nlStrain, lpm, ptFe, disp);
    nlStrain *= 0.5;
    
    totalStrain = linStrain + nlStrain;
    
    vec = stiff * totalStrain;
  }
  
  void MechPDE::CoefFunction2ndPiolaTensor::
  GetTensor(Matrix<Double> &tensor, const LocPointMapped &lpm) {
    Vector<Double> voigtVec;
    GetVector(voigtVec, lpm);
    
    switch (tensorType_) {
      case PLANE_STRAIN:
      case PLANE_STRESS:
        tensor.Resize(4, 4);
        tensor.Init();
        
        tensor[0][0] = voigtVec[0];
        tensor[0][1] = voigtVec[2];
        tensor[1][0] = voigtVec[2];
        tensor[1][1] = voigtVec[1];
        tensor[2][2] = voigtVec[0];
        tensor[2][3] = voigtVec[2];
        tensor[3][2] = voigtVec[2];
        tensor[3][3] = voigtVec[1];
        break;
        
      case AXI:
        tensor.Resize(5, 5);
        tensor.Init();
        
        tensor[0][0] = voigtVec[0];
        tensor[0][1] = voigtVec[2];
        tensor[1][0] = voigtVec[2];
        tensor[1][1] = voigtVec[1];
        tensor[2][2] = voigtVec[0];
        tensor[2][3] = voigtVec[2];
        tensor[3][2] = voigtVec[2];
        tensor[3][3] = voigtVec[1];
        tensor[4][4] = voigtVec[3];
        break;
        
      case FULL:
        tensor.Resize(9, 9);
        tensor.Init();
        
        for (UInt i=0; i<3 ; ++i) {
          tensor[i*3+0][i*3+0] = voigtVec[0];
          tensor[i*3+0][i*3+1] = voigtVec[3];
          tensor[i*3+0][i*3+2] = voigtVec[5];
          
          tensor[i*3+1][i*3+0] = voigtVec[3];
          tensor[i*3+1][i*3+1] = voigtVec[1];
          tensor[i*3+1][i*3+2] = voigtVec[4];
          
          tensor[i*3+2][i*3+0] = voigtVec[5];
          tensor[i*3+2][i*3+1] = voigtVec[4];
          tensor[i*3+2][i*3+2] = voigtVec[2];
        }
        break;
        
      default:
        EXCEPTION("Unknown mechanical subtype: " << tensorType_ << std::endl);
    }
  }
  
  SubTensorType MechPDE::GetSubTensorType(){
    if(subType_ =="axi"){return AXI;}
    if(subType_ =="planeStrain"){return PLANE_STRAIN;}
    if(subType_ =="planeStress"){return PLANE_STRAIN;}
    if(subType_ =="3d"){return FULL;}
    else{EXCEPTION( "Unknown subtype '" << subType_ << "'" );}
  }

  template BiLinearForm* MechPDE::GetPenaltyIntegrator<Double>(PtrCoefFct, Double, BiLinearForm::CouplingDirection);
  template BiLinearForm* MechPDE::GetPenaltyIntegrator<Complex>(PtrCoefFct, Double, BiLinearForm::CouplingDirection);
  template BiLinearForm* MechPDE::GetFluxIntegrator<Double>(PtrCoefFct, PtrCoefFct, Double, BiLinearForm::CouplingDirection, bool, bool, bool);
  template BiLinearForm* MechPDE::GetFluxIntegrator<Complex>(PtrCoefFct, PtrCoefFct, Double, BiLinearForm::CouplingDirection, bool, bool, bool);
} // end namespace CoupledField
