// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "mechPDE.hh"

#include <sstream>
#include <iomanip>
#include <set>

#include "Forms/mechStressStrain.hh"
#include "Forms/PMLInt.hh"
#include "Forms/linViscoElastInt.hh"
#include "Forms/linElastInt.hh"
#include "Forms/FlatShellStiffInt.hh"
#include "Forms/FlatShellMassInt.hh"
#include "Forms/massInt.hh"
#include "Forms/nonConformingInt.hh"
#include "Forms/linPressureInt.hh"
#include "Forms/singleEntryInt.hh"
#include "Forms/linSurfStressInt.hh"
#include "Forms/piezoPolarizationMatrixRHSInt.hh"
#include "Driver/assemble.hh"
#include "newmark.hh"
#include "newmarkFracDampMech.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/resultHandler.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"
#include "Domain/ansatzFct.hh"
#include "Driver/stdSolveStep.hh"
#include "Driver/solveStepMech.hh"
#include "Utils/SmoothSpline.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/SIMP.hh"

#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh"
#endif

namespace CoupledField {

DECLARE_LOG(mechpde)
DEFINE_LOG(mechpde, "mechpde")

/** the static test strain enum mapping */
Enum<MechPDE::TestStrain> MechPDE::testStrain;

MechPDE::MechPDE(Grid * aptgrid, PtrParamNode paramNode )
    :SinglePDE( aptgrid, paramNode ) {

    pdename_          = "mechanic";
    pdematerialclass_ = MECHANIC;
    maxTimeDerivOrder_ = 2;

    fracDamping_   = false;
    effectiveMass_ = false;
    nonLin_        = false;
    nonLinMaterial_= false;
    isHeatCoupled_ = false;

    needSolPrev_ = true;

    firstTime_ = true;
    useAitken_ = myParam_->Get("useAitken")->As<bool>();
    displFac_  = myParam_->Get("fsiRelaxationParam")->As<Double>();
    aitkenOmega_ = displFac_;
    aitkenOmegaPrevIter_ = displFac_;
    FSI_ = myParam_->Get("fsi")->As<bool>();

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
      stressDim_ = 6;
    }
    else if ( subType_ == "axi" && probGeo == "axi" ) {
      isaxi_ = true;
      stressDim_ = 4;
    }
    else if ( subType_ == "planeStrain" && probGeo == "plane" ) {
      stressDim_ = 3;
    }

    else if ( subType_ == "planeStress" && probGeo == "plane" ) {
      stressDim_ = 3;
    }

    else if ( subType_ == "flatShell" ) {
      stressDim_ = 3;
    }
    else {
      EXCEPTION( "Subtype '" <<  subType_ << "' of PDE '"
                 <<  pdename_ <<  "' does not fit to problem  geometry '"
                 << probGeo << "'"; );
    }
    
    // Sanity check: 3D can only be computed if 3D elements are present
    if( subType_ == "3d" && ptgrid_->GetNumElemOfDim(3) == 0 ) {
      EXCEPTION("Can not calculate 3D mechanics without 3D elements in the grid!");
    }

    //check for prestressing
    //    ReadPreStressing();

    // Register functions for scripting
    RegisterFunctions();
  }

  MechPDE::~MechPDE()
  {

  }

  void MechPDE::ReadDampingInformation( )
  {

    fracMemory_ = 0;
    bool identical = true; // i.e. same type of damping for all regions
    Integer firstFrac=-1;

    std::map<std::string, DampingType> idDampType;
    std::map<std::string, shared_ptr<RaylDampingData> > idRaylData;

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

        // make special things for fractional damping
        if( actType == FRACTIONAL ) {

          fracDamping_ = true;
          // Find first region containing fractional damping
          if ( firstFrac < 0 )
            firstFrac = i;

          // Gather additional information for fractional damping model
          std::string fracAlg = "gl";
          dampNodes[i]->GetValue( "fracAlg", fracAlg, ParamNode::PASS );

          std::string interpol = "no";
          dampNodes[i]->GetValue( "interpolation", interpol, ParamNode::PASS );
          UInt fracMem = 1;
          dampNodes[i]->GetValue( "fracMemory", fracMem, ParamNode::PASS );

          if  ( fracAlg == "gl" ) {

            // Include fracAlg and interpolation info in dampingList
            if (interpol == "no" )
              actType = FRACTIONAL_GL;
            else {
              EXCEPTION("Till now no interpolation is allowed in "
                        << "mechanics fractional damnping!" );
            }
          }

          // up to now take maximum of fracMemory
          if ( fracMem > fracMemory_ )
            fracMemory_ = fracMem;
        }

        else if( actType == RAYLEIGH ) {
          // set data for Rayleigh damping
          shared_ptr<RaylDampingData> actRaylDamp(new RaylDampingData());
          actRaylDamp->alpha = "0.0";
          actRaylDamp->beta = "0.0";
          actRaylDamp->adjustDamping = true;
          actRaylDamp->ratioDeltaF = 0.01;
          actRaylDamp->freq = 0.0;

          dampNodes[i]->GetValue( "freq", actRaylDamp->freq, ParamNode::PASS);
          dampNodes[i]->GetValue( "ratioDeltaF", actRaylDamp->ratioDeltaF, ParamNode::PASS );
          dampNodes[i]->GetValue( "adjustDamping", actRaylDamp->adjustDamping, ParamNode::PASS );
          idRaylData[actId] = actRaylDamp;
        }

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

      actRegionId = ptgrid_->GetRegion().Parse( actRegionName );

      // Check actDampingId was already registerd
      if( idDampType.count( actDampingId ) == 0 ) {
        EXCEPTION( "Damping with id '" << actDampingId
                   << "' was not defined in 'dampingList'" );
      }

      dampingList_[actRegionId] = idDampType[actDampingId];
      if ( dampingList_[actRegionId] == RAYLEIGH ){
        RaylDampingData actRayl = *(idRaylData[actDampingId]); 
        Double dampFreq;
        
        if( actRayl.freq == 0.0 ) {
        materials_[actRegionId]->GetScalar(dampFreq,RAYLEIGH_FREQUENCY,Global::REAL);
        } else { 
          dampFreq = actRayl.freq;
        }
        
        // Compute Rayleigh damping parameters
        materials_[actRegionId]->
         ComputeRayleighDamping( actRayl.alpha, actRayl.beta,
                                 dampFreq, actRayl.ratioDeltaF, 
                                 actRayl.adjustDamping, isComplex_ );
        regionRaylDamping_[actRegionId] = actRayl;

        PtrParamNode in = infoNode_->Get(ParamNode::HEADER)->GetByVal("region", "name", domain->GetGrid()->GetRegion().ToString(actRegionId));
        in->Get("alpha_M")->SetValue(actRayl.alpha);
        in->Get("alpha_K")->SetValue(actRayl.beta);
      }
    }

    // Check, if all entries are identical
    for ( UInt i = 1; i < dampingList_.size(); i++ ) {
      if ( dampingList_[subdoms_[i-1]] != dampingList_[subdoms_[i]] ) {
        identical = false;
        break;
      }
    }

    // Fractional damping can only be enabled, if all regions are damped
    // this way. Otherwise an error is thrown.
    if ( fracDamping_ == true ) {
      if ( identical == true ) {
        fracDamping_ = true;
        Info->PrintF(pdename_, "Memory size for fractional damping  is: %d\n",
                     fracMemory_ );
      } else {
        EXCEPTION( "Fractional damping can only be used if it is applied for "
                   << "ALL regions of the mechanical domain!\n"
                   << "Please check your parameter file!"; );
      }
    }
    Info->PrintF( pdename_, "\n" );
  }

  void MechPDE::ReadSpecialBCs() {

    // read volume force definition
    ReadRegionLoads();

    // check for prestressing
    ReadPreStressing();
    
    // read surface stress information
    ReadSurfStress();

    // read pressure information
    ReadPressureLoads();

  }


  void MechPDE::DefineSprings() {

    // try to get bcsAndLoads node
    PtrParamNode bcNode = myParam_->Get("bcsAndLoads", ParamNode::PASS);
    if( !bcNode )
      return;

    // fetch parameter node specifying spring
    ParamNodeList springNodes =
      bcNode->GetList("spring");

    // Iterate over all springs
    std::string name, dofName;
    Double massVal, dampVal, stiffVal;
    bool relStiff; // relative stiffness
    for( UInt i = 0; i < springNodes.GetSize(); i++ ) {

      // get data from node
      springNodes[i]->GetValue( "name", name );
      springNodes[i]->GetValue( "dof", dofName );
      springNodes[i]->GetValue( "massValue", massVal );
      springNodes[i]->GetValue( "dampingValue", dampVal );
      springNodes[i]->GetValue( "stiffnessValue", stiffVal );
      relStiff = springNodes[i]->Get("relStiffness")->As<bool>(); // other Get is ambiguous

      UInt dof = results_[0]->GetDofIndex( dofName );

      LOG_DBG2(mechpde) << "Spring: name=" << name << " dof=" << dofName << " massVal="
                    << massVal << "dampingValue=" << dampVal << "stiffnessValue="
                    << stiffVal << " relStiffness=" << relStiff;
      shared_ptr<NodeList> spNode (new NodeList(ptgrid_) );
      spNode->SetNamedNodes( name );

      // stiffness value
      if( stiffVal > EPS ) {
        // check for relative stiffness which is used for mechanism optimization
        if(relStiff) {
          double matStiffness = 0.0;

          if(materials_.size() != 1)
            EXCEPTION("relative springstiffness only for one region implemented");
          materials_.begin()->second->GetScalar(matStiffness, MECH_EMODULUS, Global::REAL );
          LOG_DBG2(mechpde) << "Set relative spring stiffness  to " << stiffVal << "*" << matStiffness;
          stiffVal *= matStiffness;
        }
        SingleEntryInt * stiffInt =
          new SingleEntryInt( lexical_cast<std::string>(stiffVal),  dof, dim_ );
        BiLinFormContext * stiffIntContext =
          new BiLinFormContext( stiffInt, STIFFNESS );
        stiffIntContext->SetPtPdes(this, this);
        stiffIntContext->SetResults( results_[0], results_[0],
                                     spNode, spNode );
        assemble_->AddBiLinearForm( stiffIntContext );
      }

      // mass Value
      if ( massVal > EPS ) {
        SingleEntryInt * massInt =
          new SingleEntryInt( lexical_cast<std::string>(massVal), dof, dim_ );
        BiLinFormContext * massIntContext =
          new BiLinFormContext( massInt, MASS );
        massIntContext->SetPtPdes(this, this);
        massIntContext->SetResults( results_[0], results_[0],
                                    spNode, spNode );
        assemble_->AddBiLinearForm( massIntContext );
      }

      // damping value
      if( dampVal > EPS ) {
        SingleEntryInt * dampInt =
          new SingleEntryInt( lexical_cast<std::string>(dampVal), dof, dim_ );
        BiLinFormContext * dampIntContext =
          new BiLinFormContext( dampInt, DAMPING );
        dampIntContext->SetPtPdes(this, this);
        dampIntContext->SetResults( results_[0], results_[0],
                                    spNode, spNode );
        assemble_->AddBiLinearForm( dampIntContext );
      }
    }
  }

  void MechPDE::ReadSoftening() {

    // Check if softeningList node is present
    std::string type, id;
    std::map<std::string, std::string> idSoftTypeMap;
    PtrParamNode softListNode = myParam_->Get("softeningList", ParamNode::PASS );
    if( softListNode ) {

      // Get child elements and read data
      ParamNodeList softNodes = softListNode->GetChildren();
      for( UInt i = 0; i < softNodes.GetSize(); i++ ) {
        type = softNodes[i]->GetName();
        softNodes[i]->GetValue( "id", id );
        idSoftTypeMap[id] = type;
      }

      if( softNodes.GetSize() ) {
        Info->PrintF( pdename_, "Applying softening for regions:\n" );
      }
    }

    // Now iterate over all regions and check for softening type
    ParamNodeList regionNodes =
      myParam_->Get("regionList")->GetList("region");
    for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {

      // get region Name
      std::string regionName = regionNodes[i]->Get("name")->As<std::string>();
      RegionIdType regionId = ptgrid_->GetRegion().Parse( regionName );

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
      Info->PrintF( pdename_, " %s: %s\n", regionName.c_str(),
                    (idSoftTypeMap[id]).c_str() );
    }
  }

  void MechPDE::InitNonLin()
  {

    nonLin_ = false;

    // ==============================================================
    // NOTE: Currently we can only treat geometric non-linearity and
    //       we assume that for a mechanic PDE all regions either
    //       are linear or non-linear!
    // ==============================================================

    // Check, if "nonLinList" is present
    PtrParamNode nonLinListNode = myParam_->Get("nonLinList", ParamNode::PASS );
    if( nonLinListNode) {

      // Get nonlinear types
      ParamNodeList nonLinNodes = nonLinListNode->GetChildren();
      for( UInt i = 0; i < nonLinNodes.GetSize(); i++ ) {

        std::string actTypeString = nonLinNodes[i]->GetName();
        std::string actId = nonLinNodes[i]->Get("id")->As<std::string>();

        NonLinType actType;
        String2Enum( actTypeString, actType );
        nonLinIdType_[actId] = actType;
      }
    }

    // Run over all regions and get entry in "regionNonLinId"
    ParamNodeList regionNodes =
      myParam_->Get("regionList")->GetChildren();

    RegionIdType actRegionId;
    std::string actRegionName, actNonLinId;

    if( regionNodes.GetSize() > 0 ) {
      Info->PrintF( pdename_, "Non-linearity in following region(s)\n" );
    }
    for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {

      // get data
      regionNodes[i]->GetValue( "name", actRegionName );
      regionNodes[i]->GetValue( "nonLinId", actNonLinId );

      if( actNonLinId == "" )
        continue;

      actRegionId = ptgrid_->GetRegion().Parse( actRegionName );

      // Check nonLinId was already registerd
      if( nonLinIdType_.find( actNonLinId) == nonLinIdType_.end() ) {
        EXCEPTION( "NonLinearity with id '" << actNonLinId
                   << "' was not defined in 'nonLinList'" );
      }

      regionNonLinId_[actRegionId] = actNonLinId;

      // ger related type of nonlinearity
      NonLinType actType = nonLinIdType_[actNonLinId];
      regionNonLinType_[actRegionId] = actType;

      // check type
      if( actType == GEOMETRIC ) {
        nonLin_ = true;
      }

      if( actType == MATERIAL ) {
        nonLin_ = true;
        nonLinMaterial_ = true;
      }

      // Log to info file
      std::string nonLinString;
      Enum2String( nonLinIdType_[actNonLinId], nonLinString );
      Info->PrintF( pdename_, " %s: %s\n", actRegionName.c_str(),
                    nonLinString.c_str() );

    }

    // Check, if nonlinear type is the same for all regions
    if( regionNonLinType_.size() > 1 ) {
      std::map<RegionIdType, NonLinType>::iterator it;
      it = regionNonLinType_.begin();
      NonLinType firstType = it->second;
      for( ; it != regionNonLinType_.end(); it++ ) {
        if( it->second != firstType ) {
          EXCEPTION( "Non-linearity should be the same for all regions!" );
        }
      }
    }

    // ------------------------------------------
    //   Get information on reduced integration
    // ------------------------------------------
    if ( nonLin_ == true ) {
      for ( UInt i = 0; i < regionNodes.GetSize(); i++ ) {
        if ( regionNodes[i]->Get("reducedInt")->As<std::string>() == "yes" ) {
          EXCEPTION( "Currently we do not support non-linearity with "
                     << "reduced integration!" );
        }
      }
    }
    Info->PrintF( pdename_, "\n" );
  }


  void MechPDE::DefineIntegrators() {

    bool isRegionPML;

    // =======================================
    //  Get information about softening types
    // =======================================
    ReadSoftening();

    // help variables for parameter checking
    RegionIdType actRegion;
    BaseMaterial* actSDMat = NULL;

    // volume integrators

    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {

      isRegionPML = false;

      // Set current region and material
      actRegion = it->first;
      actSDMat = it->second;

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );

      // get current region name and get grip of paramNode
      std::string actRegionName;
      actRegionName = ptgrid_->GetRegion().ToString( actRegion );
      // isMicroPiezo = true means, that the bilinear-forms 
      // will be defined in PiezoCupling.cc!!
      bool isMicroPiezo = IsRegionMicroPiezo( actRegionName );


      // Safety check: In general it does not make sense to allow for 
      // the same region RAYLEIGH damping AND complex-valued material
      // parameters, as both of them specify a damping behaviour at
      // the smae time
      if( dampingList_[actRegion] == RAYLEIGH &&
          complexMatData_[actRegion] ) {
        WARN("You have defined for region '" << actRegionName
            << "' both RAYLEIGH damping, "
            << "as well as complex-valued material parameters, which is "
            << "in general not recommended. Note, that the stiffness "
            << "proportional damping part is just computed from the real"
            << "part of the stiffness." );
      }
      
      //================= Check for Perfectly matched layers ====================//
      if ( dampingList_[actRegion] == PML &&  !isMicroPiezo ) {
        if ( analysistype_ != HARMONIC ) {
          EXCEPTION( "PML just supported for Harmonic-Analysis" );
        }

        isRegionPML = true;

        //read data for PML layer

        //type of PML damping
        std::string dampingTypePML;

        // inner / outer region
        Matrix<Double> inner;
        Matrix<Double> outer;
        std::string coordSysId;
        
        //damping factor
        Double dampPML;

        PtrParamNode actRegionNode =
          myParam_->Get("regionList")->GetByVal( "region", "name", actRegionName );

        std::string id = actRegionNode->Get("dampingId")->As<std::string>();
        PtrParamNode pmlNode = myParam_->Get("dampingList")->GetByVal("pml", "id", id);
        ReadDataPML(dampingTypePML, inner, dampPML, coordSysId, pmlNode );

        GetPMLLayerData(inner, outer, actRegion, coordSysId );

        //====================================================================
        //	 stiffness integrator for PML
        //====================================================================

        std::string formsType = "laplaceInt";
        SubTensorType tensorType;
        String2Enum(subType_,tensorType);

        BaseForm * bilinearStiffPML =
          new MechPMLInt(formsType, actSDMat, dampingTypePML, dampPML, tensorType);

        bilinearStiffPML->SetPosPML(inner,outer, coordSysId);

        BiLinFormContext * stiffContextPML =
          new BiLinFormContext( bilinearStiffPML, STIFFNESS );

        stiffContextPML->SetPtPdes(this, this);
        stiffContextPML->SetResults( results_[0], results_[0],
                                     actSDList, actSDList );
        // stiffContextReal->SetEntryType(matType);
        assemble_->AddBiLinearForm( stiffContextPML);


        //====================================================================
        //	 mass integrator for PML
        //====================================================================
        formsType = "massInt";

        double density;
        actSDMat->GetScalar(density,DENSITY,Global::REAL);

        //set real part
        PMLInt * bilinearMassPML =
          new PMLInt( formsType, density, dampingTypePML, dampPML, isaxi_ );

        bilinearMassPML->SetNrDofs( dim_ );
        bilinearMassPML->SetPosPML(inner,outer, coordSysId);

        BiLinFormContext * massContextPML =
          new BiLinFormContext( bilinearMassPML, MASS);

        massContextPML->SetPtPdes(this, this);
        massContextPML->SetResults( results_[0], results_[0],
                                    actSDList, actSDList );
        // massContextReal->SetEntryType( matType );
        assemble_->AddBiLinearForm( massContextPML );

      } // end of pml part

      else {
        // ==============  add "standard" linear stiffness ===========================     
        if ( regionNonLinType_[actRegion] != MATERIAL &&  !isMicroPiezo ) { 
          BaseForm * bilinearStiff = GetStiffIntegrator(actSDMat, actRegion);

          //check  for softening!
          if ( regionSoftening_.count(actRegion) ) {
            bilinearStiff->SetSofteningModel( regionSoftening_[actRegion] );
          }

          BiLinFormContext * actIntDescrStiff =
            new BiLinFormContext(bilinearStiff, STIFFNESS );

          actIntDescrStiff->SetPtPdes(this, this);
          actIntDescrStiff->SetResults( results_[0], results_[0],
                                        actSDList, actSDList );

          //check for damping
          if ( dampingList_[actRegion] == RAYLEIGH ) {
            RaylDampingData & actDamp = (regionRaylDamping_[actRegion]);
            actIntDescrStiff->SetSecDestMat(DAMPING, actDamp.beta );
          }

          assemble_->AddBiLinearForm( actIntDescrStiff );

          // check for complex valued material parameter
          if( complexMatData_[actRegion] ) {
            BaseForm * bilinearStiffImag = GetStiffIntegrator(actSDMat,
                                                              actRegion);

            Global::ComplexPart matType = Global::IMAG;
            bilinearStiffImag->SetMatDataType(matType);

            //check  for softening!
            if ( regionSoftening_.count(actRegion) ) {
              bilinearStiffImag->SetSofteningModel( regionSoftening_[actRegion] );
            }

            BiLinFormContext * actIntDescrStiffImag =
              new BiLinFormContext(bilinearStiffImag, STIFFNESS );

            actIntDescrStiffImag->SetPtPdes(this, this);
            actIntDescrStiffImag->SetResults( results_[0], results_[0],
                                              actSDList, actSDList );
            actIntDescrStiffImag->SetEntryType(matType);

            assemble_->AddBiLinearForm(actIntDescrStiffImag  );
          }
        }
      } // end linear stiffness
      // ==============  add nonlinear stiffness ============================
      if ( (nonLin_ && regionNonLinType_[actRegion] == GEOMETRIC) &&  !isMicroPiezo ) {

        BaseForm *nLinPart1 = NULL;
        BaseForm *nLinPart2 = NULL;

        if (subType_ == "3d")
          {
            nLinPart1 = new nLinMech3dInt_BNonLin(actSDMat);
            nLinPart2 = new nLinMech3dInt_PiolaStress(actSDMat);
          }

        else if (subType_ == "planeStrain")
          {
            nLinPart1 = new nLinMechPlaneStrainInt_BNonLin(actSDMat);
            nLinPart2 = new nLinMechPlaneStrainInt_PiolaStress(actSDMat);

          }
        else if (subType_ == "axi")
          {
            nLinPart1 = new nLinMechAxiInt_BNonLin(actSDMat);
            nLinPart2 = new nLinMechAxiInt_PiolaStress(actSDMat);

          }
        // pass solution pointer to nonlinear forms
        nLinPart1->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
        nLinPart2->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));

        BiLinFormContext * stiffNL1Descr =
          new BiLinFormContext(nLinPart1, STIFFNESS );

        stiffNL1Descr->SetPtPdes(this, this);
        stiffNL1Descr->SetResults( results_[0], results_[0],
                                   actSDList, actSDList );
        assemble_->AddBiLinearForm(stiffNL1Descr);

        BiLinFormContext * stiffNL2Descr =
          new BiLinFormContext(nLinPart2, STIFFNESS );

        stiffNL2Descr->SetPtPdes(this, this);
        stiffNL2Descr->SetResults( results_[0], results_[0],
                                   actSDList, actSDList );
        assemble_->AddBiLinearForm(stiffNL2Descr);

        // assemble prestress, if in config-file given
        //    if (preStressVal_)
        //      AssemblePreStressMat(ptEl, connect_PDE, ptCoord,
        //      actMatData, elDisp);
      }

      if ( (nonLin_ && regionNonLinType_[actRegion] == MATERIAL ) &&  !isMicroPiezo) {

        BaseForm *nLinMaterial = NULL;

        std::string nlfnc;
        materials_[actRegion]->GetScalar(nlfnc,NONLIN_DATA_NAME);

        ApproxData *nlinFnc = new SmoothSpline(nlfnc);
        //            ApproxData *nlinFnc = new LinInterpolate(nlfnc);
        nlinFnc->CalcBestParameter();
        nlinFnc->CalcApproximation();

        if (subType_ == "3d")
          {
            nLinMaterial = new nLinMech3dInt_Material(nlinFnc, actSDMat);
          }

        nLinMaterial->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));

        BiLinFormContext * stiffNLMaterialDescr =
          new BiLinFormContext(nLinMaterial, STIFFNESS );

        stiffNLMaterialDescr->SetPtPdes(this, this);
        stiffNLMaterialDescr->SetResults( results_[0], results_[0],
                                          actSDList, actSDList );
        assemble_->AddBiLinearForm(stiffNLMaterialDescr);

      }

      // ==============  add mass ===========================================

      if ( !isRegionPML ) {
        BiLinFormContext * actIntDescr = NULL;
        if ( subType_ != "flatShell" )
        {
          // we let the integrator obtain the density by itself, this way we can easily ensure that
          // everything goes right when we do bimaterial topology optimization
          BaseForm::MaterialDescriptor md(BaseForm::MaterialDescriptor::SCALAR, MECHANIC, DENSITY,Global::REAL);
          MassInt* bilinearMass = new MassInt(actSDMat, md, dim_, isaxi_);
          if ( diagMass_ ) {
            // diagonal mass matrix
            bilinearMass->SetDiagMass();
          }

          actIntDescr =  new BiLinFormContext(bilinearMass, MASS );

        } else {

          FlatShellMassInt * bilinearMass = new FlatShellMassInt(actSDMat, false);

          // Obtain thickness and penalty dof
          PtrParamNode actRegionNode =
            myParam_->GetByVal( "region", "name", actRegionName );

          Double thickness = actRegionNode->Get("thickness")->As<Double>();
          bilinearMass->SetThickness( thickness );

          // Get penalty value for drilling dof of region
          Double penaltyDof = actRegionNode->Get("penaltyDof")->As<Double>();
          bilinearMass->SetPenaltyDof( penaltyDof );

          actIntDescr = new BiLinFormContext(bilinearMass, MASS);
        }

        actIntDescr->SetPtPdes(this, this);
        actIntDescr->SetResults( results_[0], results_[0],
                                 actSDList, actSDList );

        // Check for damping (mass part)
        if ( dampingList_[actRegion] == RAYLEIGH ) {
          RaylDampingData & actDamp = regionRaylDamping_[actRegion];
          actIntDescr->SetSecDestMat( DAMPING, actDamp.alpha );
        }

        assemble_->AddBiLinearForm(actIntDescr);
      }



      // ==================== RHS ===========================================
      if (nonLin_ && regionNonLinType_[actRegion] != MATERIAL ) {

        LinearForm * rhsSource = new nLinMech_linFormInt(actSDMat, isaxi_);
        rhsSource->SetSolution( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
        LinearFormContext * nLinRhs =
          new LinearFormContext( rhsSource );
        nLinRhs->SetPtPde( this );
        nLinRhs->SetResult( results_[0], actSDList );
        assemble_->AddLinearForm( nLinRhs );
      }

      if (  preStressList_[actRegion] == "RHS" ) {
        Vector<Double> preStress = preStressVal_[actRegion];

        //transform the type
        SubTensorType tensorType;
        String2Enum(subType_,tensorType);
        LinearForm * rhsStress = new AddStressRHSInt(actSDMat, preStress, tensorType);

        LinearFormContext * linRhs =
          new LinearFormContext( rhsStress );
        linRhs->SetPtPde( this );
        linRhs->SetResult( results_[0], actSDList );
        assemble_->AddLinearForm( linRhs );
      }
      
      //check for prestressing
      if (  preStressList_[actRegion] == "RHS" ) {
        EXCEPTION("RHS prestressing is not implemented ");
      }
      else if ( preStressList_[actRegion] == "prescribedLHS"
          || preStressList_[actRegion] == "computeLHS" ) {
        //
        // given prestress state assembled to system matrix
        //
        Vector<Double> preStress = preStressVal_[actRegion];
        BaseForm * bilinearPreStress;
        if (subType_ == "planeStrain")
          bilinearPreStress = new PreStressIntPlaneStrain( actSDMat );
        else if (subType_ == "axi")
          bilinearPreStress = new PreStressIntAxi( actSDMat );
        else if (subType_ == "3d")
          bilinearPreStress = new PreStressInt3D( actSDMat );
        else 
          EXCEPTION("Unknown subtype in mech PDE!");               

        // check if prestress is given; if yes, then set it to the bilinear form        
        if (  preStressList_[actRegion] == "prescribedLHS" ) {
          dynamic_cast<PreStressInt*>(bilinearPreStress)->SetPreStress( preStress ); 
        }
        else if (  preStressList_[actRegion] == "computeLHS" ) {
          if ( analysistype_ == HARMONIC )
            dynamic_cast<PreStressInt*>(bilinearPreStress)->
            SetMechDisp( dynamic_cast<NodeStoreSol<Complex>&>(*sol_ ));
          else 
            dynamic_cast<PreStressInt*>(bilinearPreStress)->
            SetMechDisp( dynamic_cast<NodeStoreSol<Double>&>(*sol_ ));
        }

        BiLinFormContext * actIntDescrPre =
            new BiLinFormContext(bilinearPreStress, STIFFNESS );
        actIntDescrPre->SetPtPdes(this, this);
        actIntDescrPre->SetResults( results_[0], results_[0],
                                    actSDList, actSDList );

        assemble_->AddBiLinearForm(actIntDescrPre);
      }

      // Give entities and result to equation numbering class
      // and solution class
      eqnMap_->AddResult( *results_[0], actSDList );
    }


    // define Integrators for composite materials
    std::map<RegionIdType, Composite>::iterator compIt;
    for( compIt=compositeMaterials_.begin(); compIt!=compositeMaterials_.end();
         compIt++ ) {

      // Get current subdomain and composite material
      RegionIdType actRegion = compIt->first;
      Composite & composite = compIt->second;

      // get current region name and get grip of paramNode
      std::string actRegionName;
      actRegionName = ptgrid_->GetRegion().ToString( actRegion );
      PtrParamNode actRegionNode =
        myParam_->Get("regionList")->GetByVal( "region", "name", actRegionName );


      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );


      // Get penalty value for drilling dof of region
      Double penaltyDof = actRegionNode->Get("penaltyDof")->As<Double>();


      // ==============  add stiffness ===========================================

      // Now define integrator for this region
      FlatShellStiffInt * myInt = new FlatShellStiffInt(&composite, false);
      myInt->SetPenaltyDof( penaltyDof );
      BiLinFormContext * actIntDescrStiff =
        new BiLinFormContext( myInt, STIFFNESS);

      //check for damping
      if ( dampingList_[actRegion] == RAYLEIGH && 
          complexMatData_[actRegion] == false ) {
        RaylDampingData & actDamp = regionRaylDamping_[actRegion];
        actIntDescrStiff->SetSecDestMat(DAMPING, actDamp.beta );
      }

      actIntDescrStiff->SetPtPdes(this, this);
      actIntDescrStiff->SetResults( results_[0], results_[0],
                                    actSDList, actSDList );

      assemble_->AddBiLinearForm( actIntDescrStiff );
      eqnMap_->AddResult( *results_[0], actSDList );
      
      // check for complex valued material parameter
      if( complexMatData_[actRegion] ) {
        Global::ComplexPart matType = Global::IMAG;
        FlatShellStiffInt * myIntC = new FlatShellStiffInt(&composite, false);
        myIntC->SetMatDataType(matType);
        BiLinFormContext * actIntDescrStiffC =
            new BiLinFormContext( myIntC, STIFFNESS);
        actIntDescrStiffC->SetPtPdes(this, this);
        actIntDescrStiffC->SetResults( results_[0], results_[0],
                                       actSDList, actSDList );
        actIntDescrStiffC->SetEntryType(matType);
        assemble_->AddBiLinearForm( actIntDescrStiffC );
      }

      // ==============  add mass ===========================================
      FlatShellMassInt * bilinearMass = new FlatShellMassInt(&composite, false);
      bilinearMass->SetPenaltyDof( penaltyDof );
      BiLinFormContext * actIntDescrMass = new BiLinFormContext(bilinearMass, MASS);

      //check for damping
      if ( dampingList_[actRegion] == RAYLEIGH 
          && complexMatData_[actRegion] == false ) {
        RaylDampingData & actDamp = regionRaylDamping_[actRegion];
        actIntDescrMass->SetSecDestMat( DAMPING, actDamp.alpha );
      }

      actIntDescrMass->SetPtPdes(this, this);
      actIntDescrMass->SetResults( results_[0], results_[0],
                                   actSDList, actSDList );

      assemble_->AddBiLinearForm( actIntDescrMass );
      eqnMap_->AddResult( *results_[0], actSDList );
    }
    
    //surface integrators
    //RHS-part
    DefinePressureIntegrators(pressSurf_, pressVals_, pressPhase_);
    
    // Add integrator for surface stresses
    std::map<RegionIdType,SurfStress>::iterator  stressIt;
    for( stressIt = surfStresses_.begin();
         stressIt != surfStresses_.end(); stressIt++ ) {
      shared_ptr<SurfElemList> surfElems( new SurfElemList(ptgrid_ ) );
      shared_ptr<ElemList> volElems( new ElemList(ptgrid_ ) );
      surfElems->SetRegion( stressIt->first );
      volElems->SetRegion( stressIt->second.region );

      SurfStress3DLinForm * stressInt =
        new SurfStress3DLinForm((*stressIt).second.stress, surfElems );
      LinearFormContext * stressContext =
        new LinearFormContext( stressInt );
      stressContext->SetPtPde(this);
      stressContext->SetResult( results_[0], volElems );
      assemble_->AddLinearForm( stressContext );
      eqnMap_->AddResult( *results_[0], volElems );

    }

    // =======================================================================
    // Integrators for NonConforming Interfaces
    // =======================================================================
    PtrParamNode ncIfaceListNode
        = param->Get("domain")->Get("ncInterfaceList", ParamNode::PASS);

    // Get index of LAGRANGE_MULT result, just in case... who knows...
    UInt lmResultIdx = 0;
    for(UInt i=0, n=results_.GetSize(); i<n; i++) {
      if(results_[i]->resultType == LAGRANGE_MULT) {
        lmResultIdx = i;
        break;
      }
    }
    LOG_DBG2(mechpde) << "NonMatching: Index of LAGRANGE_MULT result: "
                     << lmResultIdx;
    
    for( UInt i = 0, n = ncIFaces_.GetSize(); i < n; i++ ) {
      // get regionId of Lagrangian surface
      StdVector<std::string> keyVec, attrVec, valVec;
      std::string slaveSide;
      std::string ncIfaceName = ptgrid_->GetRegion().ToString(ncIFaces_[i]);

      if (!ncIfaceListNode) {
        EXCEPTION("No ncInterfaces defined in domain section.");
      }
      PtrParamNode curNciNode = ncIfaceListNode->GetByVal("ncInterface", "name",
                                                   ncIfaceName);
      slaveSide = curNciNode->Get("slaveSide")->As<std::string>();

      // Part 1: Define integrator M(u, Lambda) on
      //         non-conforming interface (master/slave side)
      LOG_DBG2(mechpde) << "NonMatching: Defining nonconforming integrator"
                        << " for M on interface '"
                        << ptgrid_->GetRegion().ToString(ncIFaces_[i]) << "'.";
      shared_ptr<ElemList> actNcList( new ElemList(ptgrid_ ) );
      actNcList->SetRegion( ncIFaces_[i] );

      NonConformingInt * ncInt =
        new NonConformingInt( dim_, isaxi_ );
//      MassInt * ncInt = new MassInt( -1.0, dim_, isaxi_ );

      NcBiLinFormContext * stiffIntDescr =
        new NcBiLinFormContext( ncInt , STIFFNESS );

      // Force assembling of M(u, Lambda)^T
      stiffIntDescr->SetCounterPart( true );

      stiffIntDescr->SetPtPdes(this, this);
      stiffIntDescr->SetResults( results_[0], results_[lmResultIdx],
                                 actNcList, actNcList );

      assemble_->AddBiLinearForm( stiffIntDescr );


      // Part 2: Define integrator D(u, Lambda) on
      //         Lagrangian surface (slave side)
      LOG_DBG2(mechpde) << "NonMatching: Defining mass integrator"
                        << " for D on interface '"
                        << ptgrid_->GetRegion().ToString(ncIFaces_[i]) << "'.";
      shared_ptr<SurfElemList> actSDList( new SurfElemList(ptgrid_ ) );
      actSDList->SetRegion( ptgrid_->GetRegion().Parse( slaveSide ) );

      // D(u, Lambda) has the form of a standard mass
      // integrator with factor 1.0
      MassInt * dMatInt = new MassInt( 1.0, dim_, isaxi_ );
      BiLinFormContext * dMatContext =
        new BiLinFormContext( dMatInt, STIFFNESS );

      // Force assembling of D(u, Lambda)^T
      dMatContext->SetCounterPart( true );
      dMatContext->SetPtPdes( this, this );
      dMatContext->SetResults( results_[0], results_[lmResultIdx],
                               actSDList, actSDList );

      assemble_->AddBiLinearForm( dMatContext );

      // Give result LAGRANGE_MULT to equation numbering class
      eqnMap_->AddResult( *results_[lmResultIdx], actSDList );
    }

    // Add integrators for region loads
    DefineRegionLoadIntegrators(regionLoads_);

    // Define Springs
    DefineSprings();

    // this is a very easy function, it checks if it is designed and adds itself in case
    ReadTestStrains();
  }
  
  void MechPDE::DefinePressureIntegrators(StdVector<shared_ptr<EntityList> >& pressSurf, StdVector<std::string>& pressVals, StdVector<std::string>& pressPhase, StdVector<LinearFormContext*>* linForms){
    for (UInt actSF = 0; actSF < pressSurf.GetSize(); actSF++) {

      LinearSurfForm * rhsSrcSurf = 
        new PressureLinForm(pressVals[actSF], pressPhase[actSF], 
                            subType_, isaxi_ );
      rhsSrcSurf->SetVoluInfo( materials_ );

      LinearFormContext * pressRhs = 
        new LinearFormContext( rhsSrcSurf );
      pressRhs->SetPtPde( this );
      pressRhs->SetResult( results_[0], pressSurf[actSF] );
      if(linForms != NULL){
        linForms->Push_back(pressRhs);
      }else{
        assemble_->AddLinearForm( pressRhs );
      }
      
      // Give entities and result to equation numbering class
      // and solution class
      eqnMap_->AddResult( *results_[0], pressSurf[actSF] );
    }
    
  }
  
  void MechPDE::DefineTestStrainIntegrator(const TestStrain test, StdVector<LinearFormContext*>* linForms)
  {
    Vector<double> vals = CalcTestStrainVector(test);

    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for(it = materials_.begin(); it != materials_.end(); it++)
    {
      // Set current region and material
      RegionIdType actRegion = it->first;
      BaseMaterial *actSDMat = it->second;

      shared_ptr<ElemList> actSDList(new ElemList(ptgrid_));
      actSDList->SetRegion( actRegion );

      //transform the type
      SubTensorType tensorType;
      String2Enum(subType_,tensorType);
      LinearForm * rhsStrain = new AddStrainRHSInt(actSDMat, vals, tensorType);

      LinearFormContext *linRhs = new LinearFormContext(rhsStrain);
      linRhs->SetPtPde(this);
      linRhs->SetResult(results_[0], actSDList);

      if(linForms != NULL)
      {
        linForms->Push_back(linRhs);
      }
      else
      {
        assemble_->AddLinearForm(linRhs);
      }
    }
  }
  
  Vector<Double> MechPDE::CalcTestStrainVector(TestStrain ts, bool reduced)
  {
    // the test strain vector is of size six with one value 1.0 and the others 0.0
    Vector<double> vals;

    if(dim_ == 2 && reduced)
    {
      assert(ts == MechPDE::X || ts == MechPDE::Y || ts == MechPDE::XY);
      vals.Resize(3, 0.0);
      vals[ts == MechPDE::XY ? MechPDE::Z : ts] = 1.0; // xy goes to the third element: z
    }
    else
    {
      vals.Resize(6, 0.0);
      vals[ts] = 1.0;
    }

    return vals;
  }


  void MechPDE::DefinePolarizationMatrixIntegrators(const Vector<Double> &vals,
      StdVector<LinearFormContext*> *linForms, const int num)
  {
    
    LinearForm * polRHSMech = new PiezoPolarizationMatrixMechRHSInt(vals, num);
    
    StdVector<std::string> region_names;
    domain->GetGrid()->GetRegionNames(region_names);
    LOG_DBG(mechpde) << "region names = " << region_names.ToString();
    // FIXME: hardcoded region name!!
    shared_ptr<EntityList> entlist = domain->GetGrid()->GetEntityList(EntityList::SURF_ELEM_LIST, 
        region_names[1], EntityList::NO_TYPE);

    LinearFormContext *linRhs = new LinearFormContext(polRHSMech);
    linRhs->SetPtPde(this);
    linRhs->SetResult(results_[0], entlist);
    
    if(linForms != NULL)
      linForms->Push_back(linRhs);
    else
      assemble_->AddLinearForm(linRhs);
  }
  
  void MechPDE::DefineRegionLoadIntegrators(std::map<RegionIdType, RegionLoad>& regionLoads, StdVector<LinearFormContext*>* linForms){
    VolForceInt * forceInt;
    std::map<RegionIdType, RegionLoad>::iterator loadIt = regionLoads.begin();
    for( loadIt = regionLoads.begin(); loadIt != regionLoads.end(); loadIt++ ) {
      forceInt = (*loadIt).second.GetIntegrator();

      // Create new element list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( loadIt->first );
      LinearFormContext * forceContext =
        new LinearFormContext( forceInt );
      forceContext->SetPtPde(this);
      forceContext->SetResult( results_[0], actSDList );
      if(linForms != NULL){
        linForms->Push_back(forceContext);
      }else{
        assemble_->AddLinearForm( forceContext );
      }

      //assemble_->AddRhsSrcIntegrator( forceInt, (*loadIt).first,
      //                                (*loadIt).second.dynamics, nonLin_ );
      (*loadIt).second.ToInfo(infoNode_->Get("regionLoad"));
    }
    
  }


  BaseForm *
  MechPDE::GetStiffIntegrator(BaseMaterial* actSDMat,
                              RegionIdType regionId,
                              bool reducedInt)
  {

    // Get region name
    std::string regionName = ptgrid_->GetRegion().ToString( regionId );

    BaseForm * bilinearStiff = NULL;

    // Check for FlatShellIntegrator, as this one is a special sub-type,
    // not handled by the SubTensorType
    if (subType_ == "flatShell" ) {
      FlatShellStiffInt * myInt = new FlatShellStiffInt(actSDMat, false);

      // Get region node
      PtrParamNode actRegionNode =
        myParam_->GetByVal( "region", "name", regionName );

      // Get thickness of region
      Double thickness = actRegionNode->Get("thickness")->As<Double>();
      myInt->SetThickness( thickness );

      // Get penalty value for drilling dof of region
      Double penaltyDof = actRegionNode->Get("penaltyDof")->As<Double>();
      myInt->SetPenaltyDof( penaltyDof );

      bilinearStiff = myInt;
      myInt = NULL;
    } else {

      //transform the type
      SubTensorType tensorType;
      String2Enum(subType_,tensorType);

      if( fracDamping_ == false ) {
        bilinearStiff = new linElastInt(actSDMat, tensorType);
      }
      else {
        Double dt = 0.0;
        mHandle_ =  domain->GetMathParser()->GetNewHandle();
        MathParser * mParser =  domain->GetMathParser();
        mParser->SetExpr( mHandle_, "dt" );
        dt = mParser->Eval( mHandle_ );
        bilinearStiff = new LinViscoElastInt(actSDMat,tensorType, "modifiedStiffness",dt );
      }
    } // if "flatShell"

    return bilinearStiff;
  }



  void MechPDE::DefineSolveStep()
  {
	  if(FSI_)
		  solveStep_ = new SolveStepMech(*this);
	  else
		  solveStep_ = new StdSolveStep(*this);

  }




  // ======================================================
  // COUPLING SECTION
  // ======================================================


  void MechPDE::InitCoupling(PDECoupling * Coupling)
  {

    isIterCoupled_ = true;
    ptCoupling_   = Coupling;

    for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
      {
        if (ptCoupling_->GetOutputQuantity(i) == MECH_DISPLACEMENT)
          {
            // Intialize the memory of the coupling values
            ptCoupling_->CreateCouplingVector(i,isComplex_);
          }

        if (ptCoupling_->GetOutputQuantity(i) == MECH_VELOCITY)
          {
            // Intialize the memory of the coupling values
            ptCoupling_->CreateCouplingVector(i,isComplex_);
          }

        if (ptCoupling_->GetOutputQuantity(i) == MECH_FORCE)
          {
            // Intialize the memory of the coupling values
            ptCoupling_->CreateCouplingVector(i,isComplex_);

            //now since we need a incremental formulation, initialize some necessary vectors
            isIncrFormulation_ = true;
          }
      }

  }


  void MechPDE::CalcOutputCoupling()
  {
    UInt dof = 0;
    SolutionType quantity;
    StdVector<UInt> * couplingnodes = NULL;
    StdVector<Elem*> * couplingElems = NULL;
    SingleVector * temp_values = NULL;
    Vector<Double> * values;
    SingleVector *tempDispValues = NULL;
    SingleVector *tempDispOldValues = NULL;
    StdVector<BaseMaterial*> * materials = NULL;
    StdVector<std::string> outputRegions;
    UInt interfaceDispCoupl = 0;
    UInt interfaceVelCoupl = 0;
    UInt interfaceForceCoupl = 0;
    bool foundDisp = false;
    bool foundVel = false;
    bool foundForce = false;

    if (useAitken_ == false)
      Info->PrintF( "RELAXATION", "Relaxation Factor = %e\n", displFac_);

    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == false)
      return;
    if (!FSI_)
    {
      // loop over all output coupling quantities
      for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
      {
        quantity = ptCoupling_->GetOutputQuantity(i);
        ptCoupling_->GetOutputValues(i, temp_values);
        values = dynamic_cast<Vector<Double>*>(temp_values);

        switch (ptCoupling_->GetOutputType(i))
        {
          case NODE:
            if (quantity == MECH_DISPLACEMENT)
            {
              ptCoupling_->GetOutputNodes(i, couplingnodes);
              sol_->NodeSolutionToCoupling(*values, *couplingnodes);
            }

            if (quantity == MECH_VELOCITY)
            {
              ptCoupling_->GetOutputNodes(i, couplingnodes);
              solDeriv1_.SetAlgSysVector(getDeriv(FIRST_DERIV));
              solDeriv1_.NodeSolutionToCoupling(*values, *couplingnodes);
            }

            if (quantity == MECH_FORCE)
            {
              ptCoupling_->GetOutputNodes(i, couplingnodes);
              ptCoupling_->GetOutputElements(i, couplingElems);
              ptCoupling_->GetOppositeMaterials(i, materials);
              dof = ptCoupling_->GetOutputDof(i);

              CalcAcousticCouplingRHS(couplingElems,
                  *materials,
                  *couplingnodes,
                  *values, dof);
            }
            break;

          case ELEM:
            EXCEPTION("No Element coupling output" );
        }
      }
    } else {
      // Outer loop over all OUTPUT coupling terms
      for (UInt i=0; i < ptCoupling_->GetNumOutputCouplings(); i++)
      {
        quantity = ptCoupling_->GetOutputQuantity(i);
        if (quantity == MECH_DISPLACEMENT)
        {
          interfaceDispCoupl = i;
          foundDisp = true;
        } else if (quantity == MECH_VELOCITY) {
          interfaceVelCoupl = i;
          foundVel = true;
        } else if (quantity == MECH_FORCE) {
          interfaceForceCoupl = i;
          foundForce = true;
        } else
          EXCEPTION("Could not handle coupling quantity: " << lexical_cast<std::string>(quantity));
      }

      if (foundDisp)
      {
        if (ptCoupling_->GetOutputType(interfaceDispCoupl) != NODE)
          EXCEPTION("MECH_DISPLACEMENT must be of type NODE");

        // OutPutCoupling for MechDisplacement
        quantity = ptCoupling_->GetOutputQuantity(interfaceDispCoupl);
        ptCoupling_->GetOutputValues(interfaceDispCoupl, tempDispValues);
        ptCoupling_->GetOutputOldValues(interfaceDispCoupl, tempDispOldValues);

        Vector<Double> DispValues    = dynamic_cast<Vector<Double>&>(*tempDispValues);
        Vector<Double> DispOldValues = dynamic_cast<Vector<Double>&>(*tempDispOldValues);

        ptCoupling_->GetOutputNodes(interfaceDispCoupl, couplingnodes);
        sol_->NodeSolutionToCoupling(DispValues, *couplingnodes);

        if ( analysistype_ == TRANSIENT ) 
        {

          fixedOmega_ = displFac_;
          Info->PrintF( pdename_," Fixed Relaxation Parameter= %e\n",fixedOmega_);

          ptCoupling_->GetOutputValues(interfaceDispCoupl, temp_values);

          Vector<Double> * values = dynamic_cast<Vector<Double>*>(temp_values);
          Vector<Double> auxValues, auxOldValues;

          //ptCoupling_->GetOutputNodes(interfaceDispCoupl, couplingnodes);
          sol_->NodeSolutionToCoupling(auxValues, *couplingnodes);

          Double auxValuesMax, auxValuesL1Norm;
          auxValuesMax = abs( auxValues[0]);
          auxValuesL1Norm = abs(auxValues[0]);
          for (UInt ii=1; ii<auxValues.GetSize(); ii++ )
          {
            auxValuesL1Norm += abs(auxValues[ii]);
            if (abs(auxValues[ii]) > auxValuesMax )
              auxValuesMax = abs(auxValues[ii]);
          }
          Info->PrintF( pdename_," Linf Norm of auxValues:= %e\n",auxValuesMax );
          Info->PrintF( pdename_," L1 Norm of auxValues:= %e\n",auxValuesL1Norm);
          Info->PrintF( pdename_," L2 Norm of auxValues:= %e\n",auxValues.NormL2());

          Vector<Double> gSol, naux1, naux2;
          sol_->GetAlgSysVector(gSol);

          if (firstTime_)
          {
            sol_tn_1_.SetAlgSysVector(getOld(TIMESTEP_1));
            sol_tn_1_.GetAlgSysVector(gSolOld_ );
            firstTime_ = false;
          }
          if (useAitken_) 
          {
            calcAitkenOmega(DispValues, DispOldValues);
            naux1 = gSol * aitkenOmega_;
            naux2 = gSolOld_ * (1.0 - aitkenOmega_);
            gSol = naux1 + naux2;
          } else {
            naux1 = gSol * fixedOmega_;
            naux2 = gSolOld_ * (1.0 - fixedOmega_);
            gSol = naux1 + naux2;
          }
          sol_->SetAlgSysVector(gSol);
          gSolOld_ = gSol;
          TS_alg_->Corrector(gSol);
          sol_->NodeSolutionToCoupling((*values), *couplingnodes);

          Double valuesMax, valuesL1Norm;
          valuesMax = abs((*values)[0]);
          valuesL1Norm = abs((*values)[0]);
          for (UInt ii=1; ii < (*values).GetSize(); ii++ )
          {
            valuesL1Norm += abs((*values)[ii]);
            if (abs((*values)[ii]) > valuesMax )
              valuesMax = abs((*values)[ii]);
          }
          Info->PrintF( pdename_,"\n Linf Norm of values:= %e\n",valuesMax );
          Info->PrintF( pdename_," L1 Norm of values:= %e\n",valuesL1Norm);
          Info->PrintF( pdename_," L2 Norm of values:= %e\n",(*values).NormL2());
        }
      } else {
        if ( analysistype_ == TRANSIENT )
        {
          Vector<Double> gSol;
          sol_->GetAlgSysVector(gSol);
          TS_alg_->Corrector(gSol);
        }
      }

      if (foundVel)
      {
        ptCoupling_->GetOutputValues(interfaceVelCoupl, temp_values);
        Vector<Double> * velValues = dynamic_cast<Vector<Double>*>(temp_values);

        if(ptCoupling_->GetOutputType(interfaceVelCoupl) != NODE)
          EXCEPTION("MECH_VELOCITY must be of type NODE");

        ptCoupling_->GetOutputNodes(interfaceVelCoupl, couplingnodes);
        solDeriv1_.SetAlgSysVector(getDeriv(FIRST_DERIV));
        solDeriv1_.NodeSolutionToCoupling((*velValues), *couplingnodes);
      }

      if (foundForce)
      {
        ptCoupling_->GetOutputValues(interfaceForceCoupl, temp_values);
        Vector<Double> * values = dynamic_cast<Vector<Double>*>(temp_values);

        if(ptCoupling_->GetOutputType(interfaceForceCoupl) != NODE)
          EXCEPTION("MECH_FORCE must be of type NODE");

        ptCoupling_->GetOutputNodes(interfaceForceCoupl, couplingnodes);
        ptCoupling_->GetOutputElements(interfaceForceCoupl, couplingElems);
        ptCoupling_->GetOppositeMaterials(interfaceForceCoupl, materials);
        dof = ptCoupling_->GetOutputDof(interfaceForceCoupl);
        CalcAcousticCouplingRHS(couplingElems,*materials,*couplingnodes,*values, dof);
      }
    }
  }


  void MechPDE::CalcAcousticCouplingRHS( StdVector<Elem*> * couplingElems,
                                         StdVector<BaseMaterial*> & materials,
                                         StdVector<UInt>& couplingNodes,
                                         Vector<Double> & elemCouplingSols,
                                         UInt couplingdof )
  {

    EXCEPTION( "Not working at the moment" );
    //     Matrix<Double> ptCoord, elemMat;
    //     Vector<Double> normal, sol, forceOnElem, nSol;
    //     Elem * ptVolElem;
    //     Double sign = 0.0;
    //     Double density = 0.0;
    //     Integer matIndex = -1;

    //     elemCouplingSols.Init();

    //     for (UInt actElem=0; actElem<couplingElems->GetSize(); actElem++)
    //       {
    //         // Perform cast from volume element to surface element, since
    //         // mech-acou coupling makes only sense on surface elements
    //         SurfElem * actCoupleElem =
    //           dynamic_cast<SurfElem*> ((*couplingElems)[actElem]);


    //         BaseFE * ptElem = actCoupleElem->ptElem;
    //         StdVector<UInt> & connecth = (*couplingElems)[actElem]->connect;
    //         GetElemCoords(connecth, ptCoord);

    //         // Try to find according region for first neighbouring volume
    //         // element of the surface element
    //         matIndex = subdoms_.Find(actCoupleElem->ptVolElem1->regionId);

    //         // If first volume element does not belong to acoustic PDE, try the
    //         // second one
    //         if ( matIndex == -1 ) {
    //           matIndex = subdoms_.Find(actCoupleElem->ptVolElem2->regionId);
    //           ptVolElem = actCoupleElem->ptVolElem2;
    //           sign = -1.0 * actCoupleElem->normalSign;
    //         } else {
    //           ptVolElem = actCoupleElem->ptVolElem1;
    //           sign = 1.0 * actCoupleElem->normalSign;
    //         }

    //         if ( matIndex == -1) {
    //           EXCEPTION( "MechPDE::CalcAcousticCouplingRHS: The two volume "
    //                    << "element neighbours of surface element Nr. "
    //                    << actCoupleElem->elemNum << " do not belong to my regions!" );
    //         }

    //         // Assign correct density
    // 	materials[actElem]->GetScalar(density,DENSITY,Global::REAL);

    //         // get correct density belonging to the the neighbouring element
    //         // in the fluid subdomain
    //         //density = (*couplingMaterials)[actElem]->GetDensity();

    //         BaseForm * bilinear_mass = new MassInt(ptElem, density, isaxi_);
    //         Matrix<Double> elemmat;
    //         bilinear_mass->CalcElementMatrix(ptCoord, elemmat);
    //         delete bilinear_mass;


    //         GetDerivSolVecOfElement(sol, connecth);
    //         nSol.Resize(connecth.GetSize());   // solution in normal direction

    //         // the normal vector points outwards of the mechanical domain
    //         // (see. Kaltenbacher, "Num. Sim. of Mech. Act. & Sens." chapter 8.2)
    //         ptgrid_->CalcSurfNormal( normal, *actCoupleElem );
    //         normal *= sign;


    //         for (UInt actNode=0; actNode < connecth.GetSize(); actNode++)
    //           for (UInt actDof=0; actDof<dofspernode_; actDof++)
    //             nSol[actNode] += sol[actDof + actNode*dofspernode_] * normal[actDof];


    //         Vector<Double> forceOnElem;
    //         forceOnElem= elemmat * nSol;

    //         for (UInt actNode=0; actNode<ptCoord.GetNumRows(); actNode++)
    //           {
    //             UInt nodePos = 0;

    //             while(connecth[actNode] != couplingNodes[nodePos] && nodePos < couplingNodes.GetSize())
    //               nodePos++;
    //             elemCouplingSols[nodePos] += forceOnElem[actNode];
    //           }
    //       }
  }



  bool MechPDE::HasOutput(SolutionType output)
  {

    if (output == MECH_DISPLACEMENT || output == MECH_VELOCITY || output == MECH_FORCE)
      return true;

    return false;
  }



  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================


  void MechPDE::InitTimeStepping()
  {

    // timestepping formulation
    PtrParamNode myLinSysNode = FindLinearSystem( pdename_ );

    // <system name="acoustic"/> exists
    if( myLinSysNode ) {

      std::string str = "";
      myLinSysNode->GetValue( "timeSteppingFormulation", str,  ParamNode::PASS );
      if ( str == "effMassMatrix" ) {
        effectiveMass_ = true;
        Info->PrintF( pdename_,
                      "      * effective mass matrix timestepping\n");
      }
      else if ( str == "diagMassMatrix" ) {
        diagMass_      = true;
        effectiveMass_ = true;
        Info->PrintF( pdename_,
                      "      * diagonal mass matrix in explicit timestepping\n");
      }
      else {
        effectiveMass_ = false;
        Info->PrintF( pdename_,
                      "      * effective stiffness matrix timestepping\n");
      }
    }

    PtrParamNode systemNode = FindLinearSystem(pdename_);
    
    if ( fracDamping_ == false ) {
      if ( effectiveMass_ == true ) {
        if ( diagMass_ == true ) {
          //explicit time stepping
          TS_alg_ = new NewmarkEffMass( algsys_, systemNode, true );
        }
        else {
          TS_alg_ = new NewmarkEffMass( algsys_, systemNode, false );
        }
      }
      else if ( effectiveMass_ == false ) {
        TS_alg_ = new Newmark( algsys_, systemNode );
      }
    }
    else {
      if (effectiveMass_ == false)
        TS_alg_ = new NewmarkFracDampMech( algsys_, pdeId_,
                                           eqnMap_, ptgrid_, this,
                                           subdoms_, dampingList_,
                                           systemNode);
      else
        EXCEPTION( "This needs to be implemented!" );
    }
  }




  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void MechPDE::ReadSpecialResults() {

     // ---------------------------
    //  Determine special results
    // ---------------------------
    PtrParamNode resultsNode = myParam_->Get("storeResults", ParamNode::PASS );
    if ( !resultsNode ) return;
    PtrParamNode volNode =
      resultsNode->GetByVal("surfRegionResult", "type", "volumeAboveDefSurf", ParamNode::PASS );
    if( !volNode ) return;
    ParamNodeList volListNodes;
    UInt saveBegin=0, saveEnd=2000000, saveInc=1;
    
    if(volNode->Has("surfRegionList"))
    {  
      volListNodes = volNode->Get("surfRegionList")->GetList( "surfRegion" );

      volNode->Get("surfRegionList")->GetValue("saveBegin", saveBegin );
      volNode->Get("surfRegionList")->GetValue("saveEnd", saveEnd );
      volNode->Get("surfRegionList")->GetValue("saveInc", saveInc );
    }
    
    if( volListNodes.GetSize() > 0 ) {
      ResultHandler * resHandler = domain->GetResultHandler();

      std::string dirName;
      volNode->GetValue( "dof", dirName );
      shared_ptr<ResultInfo> vol(new ResultInfo);
      vol->resultType = MECH_DEF_VOLUME;
      vol->dofNames = "";
      vol->unit = "m^3";
      vol->definedOn = ResultInfo::SURF_REGION;
      vol->entryType = ResultInfo::SCALAR;
      vol->fctType = shared_ptr<ConstFct>(new ConstFct() );

      shared_ptr<EntityList> actList;

      for( UInt i = 0; i < volListNodes.GetSize(); i++ ) {

        // get data from node
        std::string regionName, outputIdString;
        StdVector<std::string> outputIds;
        volListNodes[i]->GetValue( "name", regionName );
        volListNodes[i]->GetValue("outputIds", outputIdString );
        SplitStringList( outputIdString, outputIds, ',' );

        actList = ptgrid_->GetEntityList( EntityList::REGION_LIST,
                                          regionName, EntityList::REGION );
        shared_ptr<BaseResult> actSol;
        if( isComplex_ ) {
          actSol = shared_ptr<BaseResult>(new Result<Complex>);
        } else {
          actSol = shared_ptr<BaseResult>(new Result<Double>);
        }
        actSol->SetResultInfo( vol );
        actSol->SetEntityList( actList );
        resultLists_[vol].Push_back( actSol );
        volAboveDefSurfDir_[actList] = dirName;
        //! the result will be written to.
        resHandler->RegisterResult( actSol, sequenceStep_, saveBegin, saveInc, 
                                    saveEnd, outputIds, "", true, true );
      }
    }
  }

  void MechPDE::DefineAvailResults()
  {
    // Check for subType
    StdVector<std::string> dispDofNames, stressDofNames;
    usePlatePenaltyDof_ = false;

    if( subType_ == "3d" ) {
      dispDofNames = "x", "y", "z";
      stressDofNames = "xx", "yy", "zz", "yz", "xz", "xy";

    } else if( subType_ == "planeStrain" ) {
      dispDofNames = "x", "y";
      stressDofNames = "xx", "yy", "xy";

    } else if( subType_ == "planeStress" ) {
      dispDofNames = "x", "y";
      stressDofNames = "xx", "yy", "xy";

    } else if( subType_ == "axi" ) {
      dispDofNames = "r", "z";
      stressDofNames = "rr", "zz", "rz", "phiphi";

    } else if( subType_ == "flatShell" ) {

      // Job to do: check, if penalty dof is used
      ParamNodeList regionNodes = myParam_->Get("regionList")->GetChildren();

      for(UInt i = 0; i < regionNodes.GetSize(); i++ ) {
        Double val = regionNodes[i]->Get("penaltyDof")->As<Double>();
        if( val != 0.0 ) {
          usePlatePenaltyDof_ = true;
          break;
        }
      }

      dispDofNames = "x", "y", "z", "tx", "ty";
      if( usePlatePenaltyDof_ ) {
        dispDofNames.Push_back( "tz" );
      }
    }

    // === MECHANIC DISPLACEMENT ===
    shared_ptr<ResultInfo> disp(new ResultInfo);
    disp->resultType = MECH_DISPLACEMENT;
    disp->dofNames = dispDofNames;
    disp->unit = "m";
    if( subType_ != "flatShell" ) {
    disp->entryType = ResultInfo::VECTOR;
    } else {
      disp->entryType = ResultInfo::TENSOR;
    }

    // check if problem is lagrange or legendre
    std::string approxType = myParam_->Get("type")->As<std::string>();
    if ( approxType == "lagrange" ) {
      shared_ptr<AnsatzFct> fct(new LagrangeFct);
      disp->fctType = fct;
      disp->definedOn = ResultInfo::NODE;
    }else if(  approxType == "spectral" ) {
      UInt order = myParam_->Get("order")->As<UInt>();
      shared_ptr<SpectralFct> fct(new SpectralFct);
      disp->definedOn = ResultInfo::PFEM;
      fct->SetOrder(order);
      disp->fctType = fct;
    } else {
      // define Legendre type
      shared_ptr<LegendreFct> fct(new LegendreFct);
      if( myParam_->Get("isIsotropic")->As<bool>() ) {
        UInt order = myParam_->Get("order")->As<UInt>();
        fct->SetIsoOrder( order );
      } else {
        Matrix<UInt> orderMat;
        ParamTools::AsTensor<unsigned int>(myParam_->Get("anisotropic"), dim_, dim_, orderMat);
        fct->SetAnisoOrder( orderMat );
      }
      disp->fctType = fct;
      disp->definedOn = ResultInfo::PFEM;
    }

    results_.Push_back( disp );
    availResults_.insert( disp);

    if ( (analysistype_ == TRANSIENT) || (analysistype_ == HARMONIC))
    {
      // === MECHANIC VELOCITY ===
      shared_ptr<ResultInfo> vel(new ResultInfo);
      vel->resultType = MECH_VELOCITY;
      vel->dofNames = dispDofNames;
      vel->unit = "m/s";
      vel->entryType = disp->entryType;
      vel->definedOn = disp->definedOn;
      vel->fctType = disp->fctType;
      availResults_.insert( vel );

      // === MECHANIC ACCELERATION ===
      shared_ptr<ResultInfo> acc(new ResultInfo);
      acc->resultType = MECH_ACCELERATION;
      acc->dofNames = dispDofNames;
      acc->unit = "m/s^2";
      acc->entryType = disp->entryType;
      acc->definedOn = disp->definedOn;
      acc->fctType = disp->fctType;
      availResults_.insert( acc );
    }

    // === MECHANIC RHS ===
    shared_ptr<ResultInfo> rhs(new ResultInfo);
    rhs->resultType = MECH_RHS_LOAD;
    rhs->dofNames = dispDofNames;
    rhs->unit = "N";
    rhs->entryType = disp->entryType;
    rhs->definedOn = disp->definedOn;
    rhs->fctType = disp->fctType;
    availResults_.insert( rhs );

    // === MECHANIC STRESS ===
    shared_ptr<ResultInfo> stress(new ResultInfo);
    stress->resultType = MECH_STRESS;
    stress->dofNames = stressDofNames;
    stress->unit =  "N/m^2";
    stress->entryType = ResultInfo::TENSOR;
    stress->definedOn = ResultInfo::ELEMENT;
    stress->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( stress );

    // === MECHANIC VON MISES STRESS (yield criterion, a scalar value)===
    shared_ptr<ResultInfo> vonMises(new ResultInfo);
    vonMises->resultType = VON_MISES_STRESS;
    vonMises->dofNames = "";
    vonMises->unit =  "";
    vonMises->entryType = ResultInfo::SCALAR;
    vonMises->definedOn = ResultInfo::ELEMENT;
    vonMises->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( vonMises );


    // === MECHANIC VON MISES STRAIN (yield criterion, a scalar value)===
    shared_ptr<ResultInfo> vonMisesStrain(new ResultInfo);
    vonMisesStrain->resultType = VON_MISES_STRAIN;
    vonMisesStrain->dofNames = "";
    vonMisesStrain->unit =  "";
    vonMisesStrain->entryType = ResultInfo::SCALAR;
    vonMisesStrain->definedOn = ResultInfo::ELEMENT;
    vonMisesStrain->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( vonMisesStrain );



    // === MECHANIC STRAIN ===
    shared_ptr<ResultInfo> strain(new ResultInfo);
    strain->resultType = MECH_STRAIN;
    strain->dofNames = stressDofNames;
    strain->unit =  "";
    strain->entryType = ResultInfo::TENSOR;
    strain->definedOn = ResultInfo::ELEMENT;
    strain->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( strain );

    // === MECHANIC ENERGY ===
    shared_ptr<ResultInfo> energy(new ResultInfo);
    energy->resultType = MECH_ENERGY;
    energy->dofNames = "";
    energy->unit = "Ws";
    energy->entryType = ResultInfo::SCALAR;
    energy->definedOn = ResultInfo::REGION;
    energy->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( energy );

    // === LUMPED_MECHANIC DISPLACEMENT ===
    shared_ptr<ResultInfo> elem_disp(new ResultInfo);
    elem_disp->resultType = LUMPED_MECH_DISPLACEMENT;
    elem_disp->dofNames = dispDofNames;
    elem_disp->unit = "m";
    if(subType_ != "flatShell") {
      elem_disp->entryType = ResultInfo::VECTOR;
    } else {
      elem_disp->entryType = ResultInfo::TENSOR;
    }
    elem_disp->definedOn = ResultInfo::ELEMENT;
    elem_disp->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert(elem_disp);

    // === PSEUDO DENSITY for SIMP ===
    shared_ptr<ResultInfo> mechPD(new ResultInfo);
    mechPD->resultType = MECH_PSEUDO_DENSITY;
    mechPD->dofNames = "";
    mechPD->unit = "";
    mechPD->entryType = ResultInfo::SCALAR;
    mechPD->definedOn = ResultInfo::ELEMENT;
    mechPD->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( mechPD );

    shared_ptr<ResultInfo> pysicalPD(new ResultInfo);
    pysicalPD->resultType = PHYSICAL_PSEUDO_DENSITY;
    pysicalPD->dofNames = "";
    pysicalPD->unit = "";
    pysicalPD->entryType = ResultInfo::SCALAR;
    pysicalPD->definedOn = ResultInfo::ELEMENT;
    pysicalPD->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( pysicalPD );
    

    // === SHAPE offset ===
    shared_ptr<ResultInfo> shape(new ResultInfo);
    shape->resultType = MECH_SHAPE;
    shape->dofNames = dispDofNames;
    shape->unit = "m";
    shape->entryType = disp->entryType;
    shape->definedOn = disp->definedOn;
    shape->fctType = disp->fctType;
    availResults_.insert( shape );

    // === HOMOGENIZED_TENSOR ===
    // calculated as a special optimization result
    shared_ptr<ResultInfo> homogenTensor(new ResultInfo);
    homogenTensor->resultType = HOMOGENIZED_TENSOR;
    homogenTensor->dofNames = "";
    homogenTensor->unit = "";
    homogenTensor->entryType = ResultInfo::TENSOR;
    homogenTensor->definedOn = ResultInfo::REGION;
    homogenTensor->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert(homogenTensor);


    // gradDisplacement as a nodal property
    // TODO a method would reduce the copy and paste!
    shared_ptr<ResultInfo> grad_x_displ(new ResultInfo);
    grad_x_displ->resultType = GRAD_X_DISPLACEMENT;
    grad_x_displ->SetVectorDOFs(dim_, isaxi_);
    grad_x_displ->unit = "m/m";
    grad_x_displ->entryType = ResultInfo::VECTOR;
    grad_x_displ->definedOn = ResultInfo::NODE;
    grad_x_displ->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert(grad_x_displ);

    shared_ptr<ResultInfo> grad_y_displ(new ResultInfo);
    grad_y_displ->resultType = GRAD_Y_DISPLACEMENT;
    grad_y_displ->SetVectorDOFs(dim_, isaxi_);
    grad_y_displ->unit = "m/m";
    grad_y_displ->entryType = ResultInfo::VECTOR;
    grad_y_displ->definedOn = ResultInfo::NODE;
    grad_y_displ->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert(grad_y_displ);

    shared_ptr<ResultInfo> grad_z_displ(new ResultInfo);
    grad_z_displ->resultType = GRAD_Z_DISPLACEMENT;
    grad_z_displ->SetVectorDOFs(dim_, isaxi_);
    grad_z_displ->unit = "m/m";
    grad_z_displ->entryType = ResultInfo::VECTOR;
    grad_z_displ->definedOn = ResultInfo::NODE;
    grad_z_displ->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert(grad_z_displ);

    // === OPT_RESULT_1/2/3 ===
    // this is added via the optimization stuff in DesignSpace.

    // ===================================
    // Check for non-conforming interfaces
    // ===================================
    StdVector<std::string> ncIfaceNames, ncIfaceNamesForPDE;
    StdVector<RegionIdType> ncIfaceIds;

    LOG_DBG2(mechpde) << "NonMatching: Checking if nonconforming "
                      << "interfaces of PDE exist in domain.";

    PtrParamNode domainNCIfaceListNode;
    domainNCIfaceListNode = param->Get("domain")->Get("ncInterfaceList", ParamNode::PASS);

    if(domainNCIfaceListNode)
    {
      PtrParamNode ncInterfaceListNode = myParam_->Get("ncInterfaceList", ParamNode::PASS);
      ParamNodeList pdeNCIfaceNodes;

      if(ncInterfaceListNode)
      {
        pdeNCIfaceNodes = ncInterfaceListNode->GetList("ncInterface");

        for (UInt i = 0; i < pdeNCIfaceNodes.GetSize(); i++) {
          std::string pdeIfaceName = pdeNCIfaceNodes[i]->Get("name")->As<std::string>();
          std::string domainIfaceName;

          PtrParamNode domainIfaceNode = domainNCIfaceListNode->GetByVal("ncInterface",
              "name",
              pdeIfaceName,
              ParamNode::PASS);
          if(!domainIfaceNode)
          {
            LOG_DBG2(mechpde) << "NonMatching: Nonconforming "
                              << "interface '" << ncIfaceNames[i]
                              << "' does not exist in domain.";

            EXCEPTION( "ncInterface referenced from PDE not defined in domain!");
          }

          ncIfaceNamesForPDE.Push_back(pdeIfaceName);
        }
        ptgrid_->GetRegion().Parse( ncIfaceNamesForPDE, ncIfaceIds );

        for (UInt i = 0; i < ncIfaceIds.GetSize(); i++) {
          ncIFaces_.Push_back(ncIfaceIds[i]);
        }

        // In the case of the presence of non-conforming interfaces,
        // a second resultdof object has to be created, which describes the
        // Lagrange multiplier
        if( ncIFaces_.GetSize() > 0 ) {
          LOG_DBG2(mechpde) << "NonMatching: Defining new ResultDof "
                            << "Lagrange Multiplier (LM).";

          LOG_DBG3(mechpde) << "NonMatching: Lagrange Multiplier DOFs: ";
          StdVector<std::string> lmDofNames;
          for( UInt i=0, n=dispDofNames.GetSize(); i<n; i++ ) {
        	  lmDofNames.Push_back("LM_" + dispDofNames[i]);
              LOG_DBG3(mechpde) << "NonMatching: " << lmDofNames[i];
          }

          shared_ptr<ResultInfo> lagr ( new ResultInfo );
          lagr->resultType = LAGRANGE_MULT;
          lagr->dofNames = lmDofNames;
          lagr->fctType = results_[0]->fctType;
          lagr->definedOn = results_[0]->definedOn;
          results_.Push_back( lagr );
        }
      }

    }

  }

  template <class TYPE>
  void MechPDE::CalcLumpedDisplacement(shared_ptr<BaseResult> base_result)
  {
    Result<TYPE>& res = dynamic_cast<Result<TYPE>&>(*base_result);

    EntityIterator it = res.GetEntityList()->GetIterator();
    assert(it.GetType() == EntityList::ELEM_LIST);

    Vector<TYPE>& val = res.GetVector();
    val.Resize(res.GetEntityList()->GetSize() * dim_);

    shared_ptr<ResultInfo> resinfo = GetResultInfo(MECH_DISPLACEMENT);
    ElemList single_elem(domain->GetGrid());
    Vector<TYPE> elem_sol;

    // traverse the elements
    for(it.Begin(); !it.IsEnd(); it++)
    {
      const Elem* elem = it.GetElem();

      single_elem.SetElement(elem);

      GetSolVecOfElement(elem_sol, single_elem.GetIterator(), resinfo);
      assert(elem_sol.GetSize() == dim_ * elem->connect.GetSize());

      // traverse first over the dimension
      for(UInt d = 0; d < dim_; d++)
      {
        TYPE sum = 0.0;

        // traverse local node results
        for(UInt idx = d; idx < elem_sol.GetSize(); idx += dim_)
          sum += elem_sol[idx];

        val[it.GetPos() * dim_ + d] = sum / (TYPE) elem->connect.GetSize();
      }
    }
  }

  void MechPDE::CalcResults( shared_ptr<BaseResult> result )
  {
    SolutionType st = result->GetResultInfo()->resultType;

    switch(st)
    {
    case MECH_DISPLACEMENT:
      if( isComplex_ ) {
        ExtractResult<Complex>( result, sol_ );
      } else {
        ExtractResult<Double>( result, sol_ );
      }
      break;

    case LUMPED_MECH_DISPLACEMENT:
      if(isComplex_)
        CalcLumpedDisplacement<Complex>(result);
      else
        CalcLumpedDisplacement<Double>(result);
      break;

    case MECH_VELOCITY:
      ExtractDerivResult( result, FIRST_DERIV );
      break;

    case MECH_ACCELERATION:
      ExtractDerivResult( result, SECOND_DERIV );
      break;

    case MECH_RHS_LOAD:
      if( isComplex_ ) {
        ExtractRhsResult<Complex>( result, results_[0] );
      } else {
        ExtractRhsResult<Double>( result, results_[0] );
      }
      break;

    case MECH_STRESS:
    case VON_MISES_STRESS:
      if( isComplex_ )
      {
        CalcStressAndStrain<Complex>(result, MECH_STRESS);
        if(st == VON_MISES_STRESS)
          CalcVonMises<Complex>(result);
      }
      else
      {
        CalcStressAndStrain<Double>(result, MECH_STRESS);
        if(st == VON_MISES_STRESS)
          CalcVonMises<Double>(result);
      }
      break;


    case MECH_STRAIN:
    case VON_MISES_STRAIN:
      if( isComplex_ )
      {
        CalcStressAndStrain<Complex>(result, MECH_STRAIN);
        if(st == VON_MISES_STRAIN)
          CalcVonMises<Complex>(result);
      }
      else
      {
        CalcStressAndStrain<Double>(result, MECH_STRAIN);
        if(st == VON_MISES_STRAIN)
          CalcVonMises<Double>(result);
      }
      break;

    case MECH_ENERGY:
      if( isComplex_ ) {
        CalcEnergy<Complex>( result );
      } else {
        CalcEnergy<Double>( result );
      }
      break;

    case MECH_DEF_VOLUME:
      if( isComplex_ ) {
        ComputeVolDefSurf<Complex>( result );
      } else {
        ComputeVolDefSurf<Double>( result );
      }
      break;

    case MECH_PSEUDO_DENSITY:
    case PHYSICAL_PSEUDO_DENSITY:
      if(domain->GetErsatzMaterial(false) == NULL) // no excpetion
        result->Init();
      else
        domain->GetErsatzMaterial()->ExtractResults(result, isComplex_);
      break;

    case MECH_SHAPE:
      ExtractNodeOffset(result);
      break;

    case HOMOGENIZED_TENSOR:
      CalcHomogenizedTensor(result);
      break;

    // the actual case is given in the result info in result
    case OPT_RESULT_1:
    case OPT_RESULT_2:
    case OPT_RESULT_3:
    case OPT_RESULT_4:
    case OPT_RESULT_5:
    case OPT_RESULT_6:
    case OPT_RESULT_7:
    case OPT_RESULT_8:
    case OPT_RESULT_9:
      // design should work, this is checked in AvailabeResults()
      domain->GetErsatzMaterial()->ExtractResults(result, isComplex_);
      break;

    case GRAD_X_DISPLACEMENT:
      if(isComplex_) CalcGradSolution<Complex>(result, 1);
                else CalcGradSolution<Double>(result, 1);
      break;

    case GRAD_Y_DISPLACEMENT:
      if(isComplex_) CalcGradSolution<Complex>(result, 2);
                else CalcGradSolution<Double>(result, 2);
      break;

    case GRAD_Z_DISPLACEMENT:
      if(isComplex_) CalcGradSolution<Complex>(result, 3);
                else CalcGradSolution<Double>(result, 3);
      break;

    default:
      WARN( "Result type not computable by mechanic PDE" );
    }
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


  template <class TYPE>
  void MechPDE::CalcVonMises(shared_ptr<BaseResult> res)
  {
    // we assume the calculation of the stresses/strains to be already performed!
    // CalcStressAndStrain<TYPE>(res, ss);

    // we don't need the entity iterator, only the stresses.
    // copy the stress result, our result is scalar
    Result<TYPE> &  actRes = dynamic_cast<Result<TYPE>&>(*res);
    Vector<TYPE>  stresses_strains = actRes.GetVector(); // copy!
    Vector<TYPE>& v_m_s = actRes.GetVector();
    v_m_s.Resize(actRes.GetEntityList()->GetSize()); // make scalar

    assert(stresses_strains.GetSize() / v_m_s.GetSize() == stressDim_);

    // element stress matrix
    Vector<TYPE> stress_strain(stressDim_);

    assert(stressDim_ == 6  || stressDim_ == 3); // don't know what to do with axi!
    const Matrix<double>& m = GetVonMisesMatrix(stressDim_ == 6 ? 3 : 2);
    Vector<TYPE> tmp;

    LOG_DBG(mechpde) << "CalcVonMises: M=" << m.ToString();

    for(unsigned int e = 0, n = v_m_s.GetSize(); e < n; e++)
    {
       for(unsigned s = 0; s < stressDim_; s++)
         stress_strain[s] = stresses_strains[e * stressDim_ + s];

       tmp = m * stress_strain;
       TYPE squared = stress_strain.Inner(tmp);

       v_m_s[e] = std::sqrt(squared);
       LOG_DBG2(mechpde) << "CVMS: e=" << 0 << " stress_strain=" << stress_strain.ToString() << " M*stress_strain=" << tmp.ToString()
                         << " squared=" << squared << "-> " << v_m_s[e];
    }
  }



  template <class TYPE>
  void MechPDE::CalcStressAndStrain(shared_ptr<BaseResult> res, SolutionType ss)
  {
    assert(ss == MECH_STRESS || ss == MECH_STRAIN);

    Vector<Double> intPoint;
    Vector<TYPE> elem;
    elem.Resize(stressDim_); // same dimension for stress and strain
    elem.Init(0);

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
    MechStressStrain<TYPE> *stress_strain = new MechStressStrain<TYPE>(actSDMat,GetSubTensorType());
    stress_strain->SetAnsatzFct( results_[0]->fctType );
    // loop over elements
    for ( it.Begin(); !it.IsEnd(); it++ ) {
      it.GetElem()->ptElem->GetCoordMidPoint(intPoint);

      //set element solution
      Matrix<TYPE> elSol;
      sol_->GetElemSolutionAsMatrix(elSol, it);
      stress_strain->SetActElemSol(elSol);

      //calculates the element stress/strain
      stress_strain->SetIntPoint(intPoint);
      if(ss == MECH_STRAIN)
        stress_strain->CalcStrainVec(elem,1,it);
      else
        stress_strain->CalcStressVec(elem,1,it);

      stress_strain->UnsetIntPoint();
      for( UInt iDof = 0; iDof < stressDim_; iDof++ ) {
        actVal[it.GetPos()*stressDim_ + iDof] = elem[iDof];
      }
    }

    delete stress_strain;
  }

 // ======================================================
  // POSTPROCESSING
  // ======================================================
  template<class TYPE>
  void MechPDE::ComputeVolDefSurf( shared_ptr<BaseResult> vals ) {


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
  TYPE MechPDE::ComputeVolElem(BaseFE * surfEl, Matrix<Double>& surfCoord,
                               Vector<TYPE> disp) {


    TYPE elemVol, averageDis;
    UInt nrSurfNodes = surfEl->GetNumNodes();


    //compute average displacement
    averageDis = 0;
    for (UInt i=0; i<nrSurfNodes; i++) {
      averageDis += disp[i];
    }
    averageDis /= (Double)nrSurfNodes;

    //compute the deformed volume
    elemVol = averageDis * surfEl->CalcVolume(surfCoord,isaxi_);

    return elemVol;
  }



  // ********************************************************
  //   Query parameter object for information about prestressing
  // ********************************************************
  void MechPDE::ReadPreStressing() {
    
    // Check, if any prestressing boundary condition is present
    PtrParamNode bcsNode = myParam_->Get("bcsAndLoads", ParamNode::PASS );
    if( !bcsNode) return;

    // Get prestressing parameter nodes
    ParamNodeList stressNodes = bcsNode->GetList("preStress");

    std::string actRegionName, type;
    for (UInt k = 0; k < stressNodes.GetSize(); k++) {

      // obtain regions
      stressNodes[k]->GetValue( "region", actRegionName );
      RegionIdType actRegion = ptgrid_->GetRegion().Parse( actRegionName );

      // fetch data

      if ( stressNodes[k]->Has( "prescribedLHS" ) ) {
        preStressList_[actRegion] = "prescribedLHS";

        Matrix<Double> preStressMat;
        ParamTools::AsTensor<double>(stressNodes[k]->Get("prescribedLHS" )->Get("value"),
            1,6, preStressMat );
        // transform to vector
        Vector<Double> preStressVec;
        preStressVec.Resize( preStressMat.GetNumCols());
        for( UInt i=0; i<preStressMat.GetNumCols(); i++)
          preStressVec[i] = preStressMat[0][i];

        preStressVal_[actRegion] = preStressVec;
      }
      else if ( stressNodes[k]->Has( "computeLHS" ) ) {
        preStressList_[actRegion] = "computeLHS";
      }
    }
  }

  void MechPDE::ReadTestStrains()
  {
    // Check, if any prestressing boundary condition is present
    PtrParamNode bcsNode = myParam_->Get("bcsAndLoads", ParamNode::PASS );
    if(!bcsNode) return;

    ParamNodeList tsn = bcsNode->GetList("testStrain");

    for (UInt i = 0; i < tsn.GetSize(); i++)
    {
       TestStrain ts = testStrain.Parse(tsn[i]->Get("strain")->As<std::string>());
       if(dim_ == 2 && (ts == Z || ts == XZ || ts == YZ))
          throw Exception("test strain '" + testStrain.ToString(ts) + "' invalid for 2D");
       DefineTestStrainIntegrator(ts);
    }
  }

  void MechPDE::ReadSurfStress() {

    // try to get bcsAndLoads node
    PtrParamNode bcNode = myParam_->Get("bcsAndLoads", ParamNode::PASS);
    if( !bcNode )
      return;
    ParamNodeList stressNodes = bcNode->GetList("surfStress");


    // iterate over all surface stress definitions
    std::string surf, volume, phase;
    Matrix<Double> valMat;
    for( UInt i = 0; i < stressNodes.GetSize(); i++ ) {

      // fetch data
      stressNodes[i]->GetValue( "name", surf );
      stressNodes[i]->GetValue( "region", volume );
      stressNodes[i]->GetValue( "phase", phase );
      ParamTools::AsTensor<double>(stressNodes[i]->Get("value"),
                                         stressDim_, 1, valMat );

      // create new surface stress definition
      SurfStress actStress;

      actStress.surface = ptgrid_->GetRegion().Parse(surf);
      actStress.region = ptgrid_->GetRegion().Parse(volume);
      actStress.phase = phase[i];
      valMat.ConvertToVec_AppendRows( actStress.stress );

      // add surface stress definition
      RegionIdType regionId =  ptgrid_->GetRegion().Parse( surf );
      surfStresses_[regionId] = actStress;
    }

  }

  void MechPDE::ReadPressureLoads()
  {
    ReadPressureLoadsFromXML(myParam_->Get("bcsAndLoads", ParamNode::PASS), pressSurf_, pressVals_, pressPhase_);
  }
  
  void MechPDE::ReadPressureLoadsFromXML(PtrParamNode bcNode, StdVector<shared_ptr<EntityList> >& pressSurf, StdVector<std::string>& pressVals, StdVector<std::string>& pressPhase) {
    if( !bcNode )
      return;
    ParamNodeList presNodes = bcNode->GetList("pressure");

    // iterate over all pressure definitions
    std::string name, value, phase;
    for( UInt i = 0; i < presNodes.GetSize(); i++ ) {

      presNodes[i]->GetValue("name", name );
      presNodes[i]->GetValue("value", value );
      presNodes[i]->GetValue("phase", phase );

      pressSurf.Push_back( 
        ptgrid_->GetEntityList( EntityList::SURF_ELEM_LIST,
                                name, EntityList::REGION ) );
      pressVals.Push_back( value );
      pressPhase.Push_back( phase );

    }
  }


  void MechPDE::CalcHomogenizedTensor(shared_ptr<BaseResult> base_result)
  {
    // TODO: Either the calculation moves here ore the data interchange
    // with ErsatzMaterial is done nice. Now steal the data from ParamNode!!! :((

    ErsatzMaterial* em = domain->GetOptimization() != NULL ? dynamic_cast<ErsatzMaterial*>(domain->GetOptimization()) : NULL;
    if(em == NULL) EXCEPTION("the result 'homogenizedTensor' requires an appropriate optimization definition.");

    Matrix<double>& tensor = em->homogenizedTensor;

    Result<double>& res = dynamic_cast<Result<double>&>(*base_result);
    Vector<double>& vals = res.GetVector();
    const unsigned int trows(tensor.GetNumRows());
    const unsigned int tcols(tensor.GetNumCols());
    vals.Resize(trows * tcols);

    for(unsigned int r = 0; r < trows; ++r)
      for(unsigned int c = 0; c < tcols; ++c)
        vals[r * trows + c] = tensor[r][c];
  }

  template <class TYPE>
  void MechPDE::CalcEnergy( shared_ptr<BaseResult> vals )
  {

    Matrix<Double> elemmat;
    Vector<TYPE> help, eldisp;
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
        sol_->GetElemSolution(eldisp, it);
        help = elemmat * eldisp;
        energy += ( help * eldisp) * 0.5;
      }
      actVal[regionIt.GetPos()] = energy;
      delete bilinear_stiff;
    }
  }


  void MechPDE::RegisterFunctions() {

    typedef FctPointer<MechPDE> FCPT;
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
    pt.Push_back( new FCPT ( this, &MechPDE::ReadRegionLoads ) );
    name.Push_back( "setRegionLoad" );


    // Now register all functions with scripting
    for (UInt i = 0; i < pt.GetSize(); i++ ) {
      Script_RegisterFct(name[i], pt[i], a[i] );
    }

  }

  void MechPDE::ExtractNodeOffset(shared_ptr<BaseResult> result) {
    ResultInfo& actRes = *(result->GetResultInfo());
    UInt numDofs = actRes.dofNames.GetSize();
    EntityIterator it = result->GetEntityList()->GetIterator();
    Vector<double>& actSol = dynamic_cast<Result<double>&>(*result).GetVector();
    actSol.Resize(numDofs * result->GetEntityList()->GetSize());
    Grid* grd = domain->GetGrid();
    for(it.Begin(); !it.IsEnd(); it++){
      Point offset;
      grd->GetNodeOffset(it.GetNode()-1, offset);
      for(UInt iDof = 0; iDof < numDofs; iDof++){
        actSol[it.GetPos()*numDofs+iDof] = offset[iDof];
      }
    }
  }

  void MechPDE::calcAitkenOmega(const Vector<Double>& displTilde, const Vector<Double>& displPrevIter)
  {
    Vector<Double> aux1;

    if (iterCoupledCounter_ == 1)
    {
      if (aitkenOmegaPrevIter_ < displFac_)
      {
        aitkenOmegaPrevIter_ = displFac_;
        aitkenOmega_ = displFac_;
      }

      deltaTildeDisplPrevIter_  = displTilde;
      deltaTildeDisplPrevIter_ -= displPrevIter;
      return;
    }

    aux1  = displTilde - displPrevIter;
    aux1 -= deltaTildeDisplPrevIter_;

    Double aux2 = aux1 * deltaTildeDisplPrevIter_;
    Double aux3 = aux1 * aux1;
    deltaTildeDisplPrevIter_  = displTilde;
    deltaTildeDisplPrevIter_ -= displPrevIter;
    deltaTildeDisplPrevIter_ *= aitkenOmega_;

    if (aux3 == 0)
    {
      aitkenOmega_ = aitkenOmegaPrevIter_;
      EXCEPTION("division by zero: can not find relaxation parameter" << std::endl \
          << "Cause: previous iteration step and current iteration step do not show any difference" );
      return;
    }
    Double aux4 = aux2 / aux3;

    aitkenOmega_ = - aitkenOmegaPrevIter_ * aux4;

    Info->PrintF( pdename_," Relaxation Parameter according \
        to Aitken= %e\n",aitkenOmega_);

    aitkenOmegaPrevIter_ = aitkenOmega_;
  }


} // end namespace CoupledField
