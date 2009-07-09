// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "mechPDE.hh"

#include <sstream>
#include <iomanip>

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
#include "Driver/assemble.hh"
#include "newmark.hh"
#include "newmarkFracDampMech.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/ParamTools.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"
#include "DataInOut/resultHandler.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"
#include "Domain/ansatzFct.hh"
#include "Driver/stdSolveStep.hh"
#include "Driver/solveStepMech.hh"
#include "Utils/SmoothSpline.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/SIMP.hh"

#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh"
#endif

namespace CoupledField {

DECLARE_LOG(mechpde)
DEFINE_LOG(mechpde, "mechpde")


MechPDE::MechPDE(Grid * aptgrid, ParamNode* paramNode )
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

    firstTime_   = true;
    useAitken_   = false;
    //displFac_    = 0.3;
    displFac_ = myParam_->Get("fsiRelaxationParam")->AsDouble();
    FSI_=myParam_->Get("fsi")->AsBool();

    // ****************************
    // DETERMINE GEOMETRY
    // ****************************

    // Get problem geometry and PDE subtype
    myParam_->Get("subType", subType_ );

    std::string probGeo;
    param->Get("domain")->Get("geometryType", probGeo );

    // Set number of degrees of freedom and
    // ensure that subtype fits to problem geometry
    if ( subType_ == "3d" && probGeo == "3d" ) {
      stressDim_ = 6;
      Info->PrintF("", "=== 3D PROBLEM\n");
    }
    else if ( subType_ == "axi" && probGeo == "axi" ) {
      isaxi_ = true;
      stressDim_ = 4;
      Info->PrintF("", "=== AXISYSMMETRIC PROBLEM\n");
    }
    else if ( subType_ == "planeStrain" && probGeo == "plane" ) {
      stressDim_ = 3;
      Info->PrintF("", "=== PLANE STRAIN PROBLEM\n");
    }

    else if ( subType_ == "planeStress" && probGeo == "plane" ) {
      stressDim_ = 3;
      Info->PrintF("", "=== PLANE STRESS PROBLEM\n");
    }

    else if ( subType_ == "flatShell" ) {
      stressDim_ = 3;
      Info->PrintF("", "=== FLAT SHELL PROBLEM\n");
    }
    else {
      EXCEPTION( "Subtype '" <<  subType_ << "' of PDE '"
                 <<  pdename_ <<  "' does not fit to problem  geometry '"
                 << probGeo << "'"; );
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
    std::map<std::string, Double>      idDampFreq;
    std::map<std::string, Double>      idDampRatioDeltaF;

    // try to get dampingList
    ParamNode * dampListNode = myParam_->Get( "dampingList", false );
    if( dampListNode ) {

      // get specific damping nodes
      StdVector<ParamNode*> dampNodes = dampListNode->GetChildren();

      for( UInt i = 0; i < dampNodes.GetSize(); i++ ) {

        std::string dampString = dampNodes[i]->GetName();
        std::string actId = dampNodes[i]->Get("id")->AsString();

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
          dampNodes[i]->Get( "fracAlg", fracAlg, false );

          std::string interpol = "no";
          dampNodes[i]->Get( "interpolation", interpol, false );
          UInt fracMem = 1;
          dampNodes[i]->Get( "fracMemory", fracMem, false );

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
          Double rayleighDampingFreq, rayleighDampingRatioDeltaF;

          if( dampNodes[i]->Has("freq") ) {
            dampNodes[i]->Get( "freq", rayleighDampingFreq);
            idDampFreq[actId] = rayleighDampingFreq;
          }
          if( dampNodes[i]->Has("RatioDeltaF") ) {
            dampNodes[i]->Get( "RatioDeltaF", rayleighDampingRatioDeltaF);
            idDampRatioDeltaF[actId] = rayleighDampingRatioDeltaF;
          }
        }

        // store damping type string
        idDampType[actId] = actType;
      }
    }

    // Run over all region and set entry in "regionNonLinId"
    StdVector<ParamNode*> regionNodes =
      myParam_->Get("regionList")->GetChildren();

    RegionIdType actRegionId;
    std::string actRegionName, actDampingId;

    if( regionNodes.GetSize() > 0 ) {
      Info->PrintF( pdename_, "Damping in following region(s)\n" );
    }

    for (UInt k = 0; k < regionNodes.GetSize(); k++) {
      regionNodes[k]->Get( "name", actRegionName );
      regionNodes[k]->Get( "dampingId", actDampingId );
      if( actDampingId == "" )
        continue;

      actRegionId = ptgrid_->RegionNameToId( actRegionName );

      // Check actDampingId was already registerd
      if( idDampType.count( actDampingId ) == 0 ) {
        EXCEPTION( "Damping with id '" << actDampingId
                   << "' was not defined in 'dampingList'" );
      }

      dampingList_[actRegionId] = idDampType[actDampingId];
      if ( dampingList_[actRegionId] == RAYLEIGH ){
        Double dampFreq, ratioDeltaF;
        materials_[actRegionId]->GetScalar(dampFreq,RAYLEIGH_FREQUENCY,Global::REAL);
        ratioDeltaF = 0.01;

        // check, if deltaF and dampingFreq can be found in map
        if( idDampFreq.find(actDampingId) != idDampFreq.end() ){
          dampFreq = idDampFreq[actDampingId];
        }

        if( idDampRatioDeltaF.find(actDampingId) != idDampRatioDeltaF.end() ){
          ratioDeltaF = idDampRatioDeltaF[actDampingId];
        }

        materials_[actRegionId]->ComputeRayleighDamping(dampFreq,ratioDeltaF);
      }

      // Log to info file
      std::string dampString;
      Enum2String( dampingList_[actRegionId], dampString );
      Info->PrintF( pdename_, " %s: %s\n", actRegionName.c_str(),
                    dampString.c_str() );
    }

    // Check, if all entries are identical
    for ( UInt i = 1; i < dampingList_.size(); i++ ) {
      if ( dampingList_[subdoms_[i-1]] != dampingList_[subdoms_[i]] ) {
        identical = false;
        break;
      }
    }

    // Fractional damping can only be enabled, if all regions are damped
    // this way. Oterhwise an error is thrown.
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
    
    // check for prestraining
    ReadPreStraining();

    // read surface stress information
    ReadSurfStress();

    // read pressure information
    ReadPressureLoads();

  }


  void MechPDE::DefineSprings() {

    // try to get bcsAndLoads node
    ParamNode * bcNode = myParam_->Get("bcsAndLoads", false);
    if( !bcNode )
      return;

    // fetch parameter node specifying spring
    StdVector<ParamNode*> springNodes =
      bcNode->GetList("spring");

    // Iterate over all springs
    std::string name, dofName;
    Double massVal, dampVal, stiffVal;
    bool relStiff; // relative stiffness
    for( UInt i = 0; i < springNodes.GetSize(); i++ ) {

      // get data from node
      springNodes[i]->Get( "name", name );
      springNodes[i]->Get( "dof", dofName );
      springNodes[i]->Get( "massValue", massVal );
      springNodes[i]->Get( "dampingValue", dampVal );
      springNodes[i]->Get( "stiffnessValue", stiffVal );
      relStiff = springNodes[i]->Get("relStiffness")->AsBool(); // other Get is ambiguous

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
          new SingleEntryInt( GenStr(stiffVal),  dof, dim_ );
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
          new SingleEntryInt( GenStr(massVal), dof, dim_ );
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
          new SingleEntryInt( GenStr(dampVal), dof, dim_ );
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
    ParamNode * softListNode = myParam_->Get("softeningList", false );
    if( softListNode ) {

      // Get child elements and read data
      StdVector<ParamNode*> softNodes = softListNode->GetChildren();
      for( UInt i = 0; i < softNodes.GetSize(); i++ ) {
        type = softNodes[i]->GetName();
        softNodes[i]->Get( "id", id );
        idSoftTypeMap[id] = type;
      }

      if( softNodes.GetSize() ) {
        Info->PrintF( pdename_, "Applying softening for regions:\n" );
      }
    }

    // Now iterate over all regions and check for softening type
    StdVector<ParamNode*> regionNodes =
      myParam_->Get("regionList")->GetList("region");
    for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {

      // get region Name
      std::string regionName = regionNodes[i]->Get("name")->AsString();
      RegionIdType regionId = ptgrid_->RegionNameToId( regionName );

      // get softeningId of current region
      id = "";
      regionNodes[i]->Get( "softeningId", id, false );
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
    ParamNode * nonLinListNode = myParam_->Get("nonLinList", false );
    if( nonLinListNode) {

      // Get nonlinear types
      StdVector<ParamNode*> nonLinNodes = nonLinListNode->GetChildren();
      for( UInt i = 0; i < nonLinNodes.GetSize(); i++ ) {

        std::string actTypeString = nonLinNodes[i]->GetName();
        std::string actId = nonLinNodes[i]->Get("id")->AsString();

        NonLinType actType;
        String2Enum( actTypeString, actType );
        nonLinIdType_[actId] = actType;
      }
    }

    // Run over all regions and get entry in "regionNonLinId"
    StdVector<ParamNode*> regionNodes =
      myParam_->Get("regionList")->GetChildren();

    RegionIdType actRegionId;
    std::string actRegionName, actNonLinId;

    if( regionNodes.GetSize() > 0 ) {
      Info->PrintF( pdename_, "Non-linearity in following region(s)\n" );
    }
    for( UInt i = 0; i < regionNodes.GetSize(); i++ ) {

      // get data
      regionNodes[i]->Get( "name", actRegionName );
      regionNodes[i]->Get( "nonLinId", actNonLinId );

      if( actNonLinId == "" )
        continue;

      actRegionId = ptgrid_->RegionNameToId( actRegionName );

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
        if ( regionNodes[i]->Get("reducedInt")->AsString() == "yes" ) {
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
    BaseMaterial * actSDMat = NULL;


    // check for complex valued material parameter
    // (defined only in the direct piezo coupling case )
    bool complexMaterial = false;
    if( isDirectCoupled_ ) {
      ParamNode * piezoNode =
        param->Get( "sequenceStep", "index", GenStr(sequenceStep_) )
        ->Get("couplingList")->Get("direct")->Get("piezoDirect", false);
      if( piezoNode ) {
        ParamNode * matNode = NULL;
        if( piezoNode )
          matNode = piezoNode->Get("materialDataType", false);
        if( matNode )
          complexMaterial =
            matNode->Get("type")->AsString() == "imagMaterialParameter";
      }
    }


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
      actRegionName = ptgrid_->RegionIdToName( actRegion );


      //================= Check for Perfectly matched layers ====================//
      if ( dampingList_[actRegion] == PML ) {
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

        //damping factor
        Double dampPML;

        ParamNode * actRegionNode =
          myParam_->Get("regionList")->Get( "region", "name", actRegionName );

        std::string id = actRegionNode->Get("dampingId")->AsString();
        ParamNode * pmlNode = myParam_->Get("dampingList")->Get("pml", "id", id);
        ReadDataPML(dampingTypePML, inner, dampPML, pmlNode );

        GetPMLLayerData(inner, outer, actRegion);

        //====================================================================
        //	 stiffness integrator for PML
        //====================================================================

        std::string formsType = "laplaceInt";
        SubTensorType tensorType;
        String2Enum(subType_,tensorType);

        BaseForm * bilinearStiffPML =
          new MechPMLInt(formsType, actSDMat, dampingTypePML, dampPML, tensorType);

        bilinearStiffPML->SetPosPML(inner,outer);

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
        bilinearMassPML->SetPosPML(inner,outer);

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

        if ( regionNonLinType_[actRegion] != MATERIAL ) {
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
          if ( dampingList_[actRegion] == RAYLEIGH && complexMaterial == false ) {
            Double beta, measFreq;
            std::string fac;
            actSDMat->GetScalar(beta,RAYLEIGH_BETA,Global::REAL);
            actSDMat->GetScalar(measFreq,RAYLEIGH_FREQUENCY,Global::REAL);
            if( isComplex_ ) {
              // assemble string describing adjusted Rayleigh damping
              fac = GenStr( beta * measFreq) + "/ f";
            } else {
              fac = GenStr( beta );
            }
            actIntDescrStiff->SetSecDestMat(DAMPING, fac );
          }

          assemble_->AddBiLinearForm( actIntDescrStiff );

          // check for complex valued material parameter
          if( complexMaterial ) {
            BaseForm * bilinearStiffImag = GetStiffIntegrator(actSDMat,
                                                              actRegion);

            Global::ComplexPart matType = Global::IMAG;
            bilinearStiffImag->SetMatDataType(matType);

            //check  for softening!
            if ( regionSoftening_.count(actRegion) ) {
              std::cerr << "Applying softening for region " << actRegionName << std::endl;
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
        //check for prestressing
        //    if ( isPreStresss[

        //for prestressing
        //       for ( UInt preStr=0; preStr<preStressDomain_.GetSize(); preStr++ ) {
        //         if ( actRegion == preStressDomain_[preStr]) {
        //           Vector<Double> preStrVal(3);
        //           preStrVal[0] = preStressValX_[preStr];
        //           preStrVal[1] = preStressValY_[preStr];
        //           preStrVal[2] = preStressValZ_[preStr];

        //           BaseForm * bilinearPreStress;
        //           if (subType_ == "planeStrain")
        //             bilinearPreStress = new PreStressIntPlaneStrain(actSDMat, preStrVal);
        //           else if (subType_ == "axi")
        //             bilinearPreStress = new PreStressIntAxi(actSDMat, preStrVal);
        //           else if (subType_ == "3d")
        //             bilinearPreStress = new PreStressInt3D(actSDMat, preStrVal);
        //           else
        //             Info->Error("Unknown subtype in mech PDE! ",__FILE__,__LINE__);

        //           BiLinFormContext * actIntDescrPre =
        //             new BiLinFormContext(bilinearPreStress, STIFFNESS );
        //           actIntDescrPre->SetPtPdes(this, this);
        //           actIntDescrPre->SetResults( results_[0], results_[0],
        //                                       actSDList, actSDList );

        //           assemble_->AddBiLinearForm(actIntDescrPre);
        //         }
        //       }


      // ==============  add nonlinear stiffness ============================
      if (nonLin_ && regionNonLinType_[actRegion] == GEOMETRIC) {

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

      if (nonLin_ && regionNonLinType_[actRegion] == MATERIAL ) {

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
        if ( subType_ != "flatShell" ) {

          double density;
          actSDMat->GetScalar(density,DENSITY,Global::REAL);
          MassInt * bilinearMass  = new MassInt(density, dim_, isaxi_);
          if ( diagMass_ ) {
            // diagonal mass matrix
            bilinearMass->SetDiagMass();
          }

          actIntDescr =  new BiLinFormContext(bilinearMass, MASS );


        } else {

          FlatShellMassInt * bilinearMass = new FlatShellMassInt(actSDMat, false);

          // Obtain thickness and penalty dof
          ParamNode * actRegionNode =
            myParam_->Get( "region", "name", actRegionName );

          Double thickness = actRegionNode->Get("thickness")->AsDouble();
          bilinearMass->SetThickness( thickness );

          // Get penalty value for drilling dof of region
          Double penaltyDof = actRegionNode->Get("penaltyDof")->AsDouble();
          bilinearMass->SetPenaltyDof( penaltyDof );

          actIntDescr = new BiLinFormContext(bilinearMass, MASS);
        }

        actIntDescr->SetPtPdes(this, this);
        actIntDescr->SetResults( results_[0], results_[0],
                                 actSDList, actSDList );

        // Check for damping (mass part)
        if ( dampingList_[actRegion] == RAYLEIGH ) {
          Double alpha, measFreq;
          std::string fac;
          actSDMat->GetScalar(alpha,RAYLEIGH_ALPHA,Global::REAL);
          actSDMat->GetScalar(measFreq,RAYLEIGH_FREQUENCY,Global::REAL);
          if( isComplex_ ) {
            // assemble string describing adjusted Rayleigh damping
            fac = GenStr<Double>( alpha / measFreq ) + "* f";
          } else {
            fac = GenStr( alpha );
          }
          actIntDescr->SetSecDestMat( DAMPING, fac );
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
      actRegionName = ptgrid_->RegionIdToName( actRegion );
      ParamNode * actRegionNode =
        myParam_->Get("regionList")->Get( "region", "name", actRegionName );


      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );


      // Get penalty value for drilling dof of region
      Double penaltyDof = actRegionNode->Get("penaltyDof")->AsDouble();


      // ==============  add stiffness ===========================================

      // Now define integrator for this region
      FlatShellStiffInt * myInt = new FlatShellStiffInt(&composite, false);
      myInt->SetPenaltyDof( penaltyDof );
      BiLinFormContext * actIntDescrStiff =
        new BiLinFormContext( myInt, STIFFNESS);

      //check for damping
      if ( dampingList_[actRegion] == RAYLEIGH && complexMaterial == false ) {
        Double beta, measFreq;
        std::string fac;
        actSDMat->GetScalar(beta,RAYLEIGH_BETA,Global::REAL);
        actSDMat->GetScalar(measFreq,RAYLEIGH_FREQUENCY,Global::REAL);
        if( isComplex_ ) {
          // assemble string describing adjusted Rayleigh damping
          fac = GenStr( beta * measFreq) + "/ f";
        } else {
          fac = GenStr( beta );
        }
        actIntDescrStiff->SetSecDestMat(DAMPING, fac );
      }

      actIntDescrStiff->SetPtPdes(this, this);
      actIntDescrStiff->SetResults( results_[0], results_[0],
                                    actSDList, actSDList );

      assemble_->AddBiLinearForm( actIntDescrStiff );
      eqnMap_->AddResult( *results_[0], actSDList );

      // ==============  add mass ===========================================
      FlatShellMassInt * bilinearMass = new FlatShellMassInt(&composite, false);
      bilinearMass->SetPenaltyDof( penaltyDof );
      BiLinFormContext * actIntDescrMass = new BiLinFormContext(bilinearMass, MASS);

      //check for damping
      if ( dampingList_[actRegion] == RAYLEIGH && complexMaterial == false ) {
        Double alpha, measFreq;
        std::string fac;
        actSDMat->GetScalar(alpha,RAYLEIGH_ALPHA,Global::REAL);
        actSDMat->GetScalar(measFreq,RAYLEIGH_FREQUENCY,Global::REAL);
        if( isComplex_ ) {
          // assemble string describing adjusted Rayleigh damping
          fac = GenStr( alpha / measFreq ) + "* f";
        } else {
          fac = GenStr( alpha );
        }
        actIntDescrMass->SetSecDestMat( DAMPING, fac );
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
    
    // make sure not to call with empty vector!
    if(preStrainVal_.GetSize() > 0)
      DefineTestStrainIntegrators(preStrainVal_);
    
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
    ParamNode* ncIfaceListNode
        = param->Get("domain")->Get("ncInterfaceList", false);

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
      std::string ncIfaceName = ptgrid_->RegionIdToName(ncIFaces_[i]);

      if (!ncIfaceListNode) {
        EXCEPTION("No ncInterfaces defined in domain section.");
      }
      ParamNode* curNciNode = ncIfaceListNode->Get("ncInterface", "name",
                                                   ncIfaceName);
      slaveSide = curNciNode->Get("slaveSide")->AsString();

      // Part 1: Define integrator M(u, Lambda) on
      //         non-conforming interface (master/slave side)
      LOG_DBG2(mechpde) << "NonMatching: Defining nonconforming integrator"
                        << " for M on interface '"
                        << ptgrid_->RegionIdToName(ncIFaces_[i]) << "'.";
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
                        << ptgrid_->RegionIdToName(ncIFaces_[i]) << "'.";
      shared_ptr<SurfElemList> actSDList( new SurfElemList(ptgrid_ ) );
      actSDList->SetRegion( ptgrid_->RegionNameToId( slaveSide ) );

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
  }
  
  void MechPDE::DefinePressureIntegrators(StdVector<shared_ptr<EntityList> >& pressSurf, StdVector<std::string>& pressVals, StdVector<std::string>& pressPhase, std::set<LinearFormContext*>* linForms){
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
        linForms->insert(pressRhs);
      }else{
        assemble_->AddLinearForm( pressRhs );
      }
      
      // Give entities and result to equation numbering class
      // and solution class
      eqnMap_->AddResult( *results_[0], pressSurf[actSF] );
    }
    
  }
  
  void MechPDE::DefineTestStrainIntegrators(const Vector<Double> &vals, std::set<LinearFormContext*>* linForms)
  {
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
        linForms->insert(linRhs);
      }
      else
      {
        assemble_->AddLinearForm(linRhs);
      }
    }
  }
  
  void MechPDE::DefineRegionLoadIntegrators(std::map<RegionIdType, RegionLoad>& regionLoads, std::set<LinearFormContext*>* linForms){
    VolForceInt * forceInt;
    std::map<RegionIdType, RegionLoad>::iterator loadIt = regionLoads.begin();
    if (regionLoads.size() != 0 ) {
      (*loadIt).second.Print(true, pdename_ );
    }
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
        linForms->insert(forceContext);
      }else{
        assemble_->AddLinearForm( forceContext );
      }

      //assemble_->AddRhsSrcIntegrator( forceInt, (*loadIt).first,
      //                                (*loadIt).second.dynamics, nonLin_ );
      (*loadIt).second.Print(false, pdename_);

    }
    
  }


  BaseForm *
  MechPDE::GetStiffIntegrator(BaseMaterial* actSDMat,
                              RegionIdType regionId,
                              bool reducedInt)
  {

    // Get region name
    std::string regionName = ptgrid_->RegionIdToName( regionId );

    BaseForm * bilinearStiff = NULL;

    // Check for FlatShellIntegrator, as this one is a special sub-type,
    // not handled by the SubTensorType
    if (subType_ == "flatShell" ) {
      FlatShellStiffInt * myInt = new FlatShellStiffInt(actSDMat, false);

      // Get region node
      ParamNode * actRegionNode =
        myParam_->Get( "region", "name", regionName );

      // Get thickness of region
      Double thickness = actRegionNode->Get("thickness")->AsDouble();
      myInt->SetThickness( thickness );

      // Get penalty value for drilling dof of region
      Double penaltyDof = actRegionNode->Get("penaltyDof")->AsDouble();
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





  void MechPDE::CalcOutputCoupling() {
	  UInt dof = 0;
	  SolutionType quantity;
	  StdVector<UInt> * couplingnodes = NULL;
	  StdVector<Elem*> * couplingElems = NULL;
	  SingleVector * temp_values = NULL;
	  Vector<Double> * values;
	  SingleVector *tempDispValues=NULL;
	  SingleVector *tempDispOldValues=NULL;
	  StdVector<BaseMaterial*> * materials = NULL;
	  StdVector<std::string> outputRegions;
	  UInt interfaceDispCoupl=0, interfaceVelCoupl=0, interfaceForceCoupl=0;
	  bool foundDisp=false, foundVel=false, foundForce =false;

	  if(useAitken_==false)
		  Info->PrintF( "RELAXATION", "Relaxation Factor = %e\n",displFac_);

	  // at first, check if this PDE is iterative coupled
	  if (isIterCoupled_ == false)
		  return;
	  if(!FSI_){
		  // loop over all output coupling quantities
		  for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++)
		  {
			  quantity = ptCoupling_->GetOutputQuantity(i);
			  ptCoupling_->GetOutputValues(i, temp_values);

			  values = dynamic_cast<Vector<Double>*>(temp_values);

			  switch(ptCoupling_->GetOutputType(i))
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
					  solDeriv1_.SetAlgSysVector(getS1());
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

	  }
	  else{

		  // Outer loop over all OUTPUT coupling terms
		  for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++) {
			  quantity = ptCoupling_->GetOutputQuantity(i);
			  if(quantity == MECH_DISPLACEMENT){
				  interfaceDispCoupl=i;
				  foundDisp=true;
			  }
			  else if(quantity == MECH_VELOCITY){
				  interfaceVelCoupl=i;
				  foundVel=true;
			  }
			  else if(quantity == MECH_FORCE){
				  interfaceForceCoupl=i;
				  foundForce=true;
			  }
			  else
				  EXCEPTION("Could not handle coupling quantity: " << GenStr(quantity));
		  }

		  if(foundDisp){

			  if(ptCoupling_->GetOutputType(interfaceDispCoupl) != NODE)
				  EXCEPTION("MECH_DISPLACEMENT must be of type NODE");


			  // OutPutCoupling for MechDisplacement
			  quantity = ptCoupling_->GetOutputQuantity(interfaceDispCoupl);
			  ptCoupling_->GetOutputValues(interfaceDispCoupl, tempDispValues);
			  ptCoupling_->GetOutputOldValues(interfaceDispCoupl, tempDispOldValues);

			  Vector<Double> DispValues    = dynamic_cast<Vector<Double>&>(*tempDispValues);
			  Vector<Double> DispOldValues = dynamic_cast<Vector<Double>&>(*tempDispOldValues);

			  ptCoupling_->GetOutputNodes(interfaceDispCoupl, couplingnodes);
			  sol_->NodeSolutionToCoupling(DispValues, *couplingnodes);

			  if ( analysistype_ == TRANSIENT ){
				  if(useAitken_){
					  actDelta_ = DispOldValues-DispValues;

					  if(actDelta_.GetSize() != oldDelta_.GetSize()){
						  oldDelta_.Resize(actDelta_.GetSize());
						  oldDelta_.Init();
					  }

					  Vector<Double> aux1;
					  aux1 = (oldDelta_ - actDelta_);

					  Double aux2=aux1*actDelta_;
					  Double aux3=aux1*aux1;
					  Double aux4 = aux2/aux3;

					  aitkenMu_=(aitkenMu_)+((aitkenMu_-1.0)*aux4);
					  if (aitkenMu_ > 1.0-displFac_)
						  aitkenMu_=1.0-displFac_;
					  aitkenOmega_ = 1.0-aitkenMu_;
					  Info->PrintF( pdename_," unbounded relaxation Parameter according to Aitgken= %e\n",aitkenOmega_);

					  Double omegaMax=0.3;
					  if (iterCoupledCounter_>5)
						  omegaMax=1.0;
					  if (aitkenOmega_ > omegaMax)
						  aitkenOmega_=omegaMax;
					  oldDelta_ = actDelta_;
				  }

				  fixedOmega_ = displFac_;

				  Info->PrintF( pdename_," Relaxation Parameter according to Aitken= %e\n",aitkenOmega_);
				  Info->PrintF( pdename_," Fixed Relaxation Parameter= %e\n",fixedOmega_);


				  ptCoupling_->GetOutputValues(interfaceDispCoupl, temp_values);

				  Vector<Double> * values    = dynamic_cast<Vector<Double>*>(temp_values);
				  Vector<Double> auxValues, auxOldValues;

				  //ptCoupling_->GetOutputNodes(interfaceDispCoupl, couplingnodes);
				  sol_->NodeSolutionToCoupling(auxValues, *couplingnodes);

				  Double auxValuesMax, auxValuesL1Norm;
				  auxValuesMax=abs( auxValues[0]);
				  auxValuesL1Norm=abs(auxValues[0]);
				  for(UInt ii=1; ii<auxValues.GetSize(); ii++ ){
					  auxValuesL1Norm+=abs(auxValues[ii]);
					  if (abs(auxValues[ii])> auxValuesMax )
						  auxValuesMax=abs(auxValues[ii]);
				  }
				  Info->PrintF( pdename_," Linf Norm of auxValues:= %e\n",auxValuesMax );
				  Info->PrintF( pdename_," L1 Norm of auxValues:= %e\n",auxValuesL1Norm);
				  Info->PrintF( pdename_," L2 Norm of auxValues:= %e\n",auxValues.NormL2());

				  Vector<Double> gSol, naux1, naux2;
				  sol_->GetAlgSysVector(gSol );

				  if(firstTime_){
					  sol_tn_1_.SetAlgSysVector(getOld1());
					  sol_tn_1_.GetAlgSysVector(gSolOld_ );

					  firstTime_=false;
				  }
				  if(useAitken_ == true ){
					  naux1=gSol*aitkenOmega_;
					  naux2=gSolOld_*(1.0-aitkenOmega_);
					  gSol= naux1 + naux2;
				  }
				  else {
					  naux1=gSol*fixedOmega_;
					  naux2=gSolOld_*(1.0-fixedOmega_);
					  gSol=naux1 + naux2;
				  }
				  sol_->SetAlgSysVector(gSol );

				  gSolOld_= gSol;

				  TS_alg_->Corrector(gSol);

				  sol_->NodeSolutionToCoupling((*values), *couplingnodes);

				  Double valuesMax, valuesL1Norm;
				  valuesMax=abs((*values)[0]);
				  valuesL1Norm=abs((*values)[0]);
				  for(UInt ii=1; ii<(*values).GetSize(); ii++ ){
					  valuesL1Norm+=abs((*values)[ii]);
					  if (abs((*values)[ii])> valuesMax )
						  valuesMax=abs((*values)[ii]);
				  }
				  Info->PrintF( pdename_,"\n Linf Norm of values:= %e\n",valuesMax );
				  Info->PrintF( pdename_," L1 Norm of values:= %e\n",valuesL1Norm);
				  Info->PrintF( pdename_," L2 Norm of values:= %e\n",(*values).NormL2());

			  }
		  }
		  else {
			  if ( analysistype_ == TRANSIENT ){
				  Vector<Double> gSol;
				  sol_->GetAlgSysVector(gSol );

				  TS_alg_->Corrector(gSol);
			  }
		  }

		  if(foundVel){
			  ptCoupling_->GetOutputValues(interfaceVelCoupl, temp_values);

			  Vector<Double> * velValues    = dynamic_cast<Vector<Double>*>(temp_values);

			  if(ptCoupling_->GetOutputType(interfaceVelCoupl) != NODE)
				  EXCEPTION("MECH_VELOCITY must be of type NODE");

			  ptCoupling_->GetOutputNodes(interfaceVelCoupl, couplingnodes);
			  solDeriv1_.SetAlgSysVector(getS1());
			  solDeriv1_.NodeSolutionToCoupling((*velValues), *couplingnodes);
		  }



		  if (foundForce) {
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
    //           (*error) << "MechPDE::CalcAcousticCouplingRHS: The two volume "
    //                    << "element neighbours of surface element Nr. "
    //                    << actCoupleElem->elemNum << " do not belong to my regions!";
    //           Error( __FILE__, __LINE__ );
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
    ParamNode* myLinSysNode = FindLinearSystem( pdename_ );

    // <system name="acoustic"/> exists
    if( myLinSysNode ) {

      std::string str = "";
      myLinSysNode->Get( "timeSteppingFormulation", str,  false );
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

    if ( fracDamping_ == false ) {
      if ( effectiveMass_ == true ) {
        if ( diagMass_ == true ) {
          //explicit time stepping
          TS_alg_ = new NewmarkEffMass( algsys_, pdename_, true );
        }
        else {
          TS_alg_ = new NewmarkEffMass( algsys_, pdename_, false );
        }
      }
      else if ( effectiveMass_ == false ) {
        TS_alg_ = new Newmark( algsys_, pdename_ );
      }
    }
    else {
      if (effectiveMass_ == false)
        TS_alg_ = new NewmarkFracDampMech( algsys_, pdeId_,
                                           eqnMap_, ptgrid_, this,
                                           subdoms_, dampingList_ );
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
    ParamNode * resultsNode = myParam_->Get("storeResults", false );
    if ( !resultsNode ) return;
    ParamNode * volNode =
      resultsNode->Get("surfRegionResult", "type", "volumeAboveDefSurf", false );
    if( !volNode ) return;
    StdVector<ParamNode*> volListNodes =
      volNode->Get("surfRegionList")->GetList( "surfRegion" );

    UInt saveBegin, saveEnd, saveInc;
    volNode->Get("surfRegionList")->Get("saveBegin", saveBegin );
    volNode->Get("surfRegionList")->Get("saveEnd", saveEnd );
    volNode->Get("surfRegionList")->Get("saveInc", saveInc );

    if( volListNodes.GetSize() > 0 ) {
      ResultHandler * resHandler = domain->GetResultHandler();

      std::string dirName;
      volNode->Get( "dof", dirName );
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
        volListNodes[i]->Get( "name", regionName );
        volListNodes[i]->Get("outputIds", outputIdString );
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
        resHandler->RegisterResult( actSol, saveBegin, saveInc, saveEnd,
                                    outputIds, "", true, true );
      }
    }
  }

  void MechPDE::DefineAvailResults() {

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
      StdVector<ParamNode*> regionNodes = myParam_->Get("regionList")->GetChildren();

      for(UInt i = 0; i < regionNodes.GetSize(); i++ ) {
        Double val = regionNodes[i]->Get("penaltyDof")->AsDouble();
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
    std::string approxType = myParam_->Get("type")->AsString();
    if ( approxType == "lagrange" ) {
      shared_ptr<AnsatzFct> fct(new LagrangeFct);
      disp->fctType = fct;
      disp->definedOn = ResultInfo::NODE;
    }else if(  approxType == "spectral" ) {
      UInt order = myParam_->Get("order")->AsUInt();
      shared_ptr<SpectralFct> fct(new SpectralFct);
      disp->definedOn = ResultInfo::PFEM;
      fct->SetOrder(order);
      disp->fctType = fct;
    } else {

       // define Legendre type
       shared_ptr<LegendreFct> fct(new LegendreFct);
       if( myParam_->Get("isIsotropic")->AsBool() ) {
         UInt order = myParam_->Get("order")->AsUInt();
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

    // === PSEUDO DENSITY for simp ===
    shared_ptr<ResultInfo> pseudoDensity(new ResultInfo);
    pseudoDensity->resultType = MECH_PSEUDO_DENSITY;
    pseudoDensity->dofNames = "";
    pseudoDensity->unit = "";
    pseudoDensity->entryType = ResultInfo::SCALAR;
    pseudoDensity->definedOn = ResultInfo::ELEMENT;
    pseudoDensity->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( pseudoDensity );
    
    // === SHAPE offset ===
    shared_ptr<ResultInfo> shape(new ResultInfo);
    shape->resultType = MECH_SHAPE;
    shape->dofNames = dispDofNames;
    shape->unit = "m";
    shape->entryType = disp->entryType;
    shape->definedOn = disp->definedOn;
    shape->fctType = disp->fctType;
    availResults_.insert( shape );

    // === OPT_RESULT_1/2/3 ===
    // this is added via the optimization stuff in DesignSpace.

    // ===================================
    // Check for non-conforming interfaces
    // ===================================
    StdVector<std::string> ncIfaceNames, ncIfaceNamesForPDE;
    StdVector<RegionIdType> ncIfaceIds;

    LOG_DBG2(mechpde) << "NonMatching: Checking if nonconforming "
                      << "interfaces of PDE exist in domain.";

    ParamNode* domainNCIfaceListNode;
    domainNCIfaceListNode = param->Get("domain")->Get("ncInterfaceList", false);

    if(domainNCIfaceListNode)
    {
      ParamNode* ncInterfaceListNode =
        param->Get("sequenceStep", "index", GenStr(sequenceStep_) )
        ->Get("pdeList")->Get("mechanic")->Get("ncInterfaceList", false);
      StdVector<ParamNode*> pdeNCIfaceNodes;

      if(ncInterfaceListNode)
      {
        pdeNCIfaceNodes = ncInterfaceListNode->GetList("ncInterface");

        for (UInt i = 0; i < pdeNCIfaceNodes.GetSize(); i++) {
          std::string pdeIfaceName = pdeNCIfaceNodes[i]->Get("name")->AsString();
          std::string domainIfaceName;

          ParamNode* domainIfaceNode = domainNCIfaceListNode->Get("ncInterface",
              "name",
              pdeIfaceName,
              false);
          if(!domainIfaceNode)
          {
            LOG_DBG2(mechpde) << "NonMatching: Nonconforming "
                              << "interface '" << ncIfaceNames[i]
                              << "' does not exist in domain.";

            EXCEPTION( "ncInterface referenced from PDE not defined in domain!");
          }

          ncIfaceNamesForPDE.Push_back(pdeIfaceName);
        }
        ptgrid_->RegionNameToId( ncIfaceIds, ncIfaceNamesForPDE );

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

  void MechPDE::CalcResults( shared_ptr<BaseResult> result ) {

    switch (result->GetResultInfo()->resultType ) {

    case MECH_DISPLACEMENT:
      if( isComplex_ ) {
        ExtractResult<Complex>( result, sol_ );
      } else {
        ExtractResult<Double>( result, sol_ );
      }
      break;

    case MECH_VELOCITY:
      ExtractDerivResult( result, 1 );
      break;

    case MECH_ACCELERATION:
      ExtractDerivResult( result, 2 );
      break;

    case MECH_RHS_LOAD:
      if( isComplex_ ) {
        ExtractRhsResult<Complex>( result, results_[0] );
      } else {
        ExtractRhsResult<Double>( result, results_[0] );
      }
      break;

    case MECH_STRESS:
      if( isComplex_ ) {
        CalcStresses<Complex>( result );
      } else {
        CalcStresses<Double>( result );
      }
      break;

    case MECH_STRAIN:
      if( isComplex_ ) {
        CalcStrains<Complex>( result );
      } else {
        CalcStrains<Double>( result );
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
      if(domain->GetErsatzMaterial(false) == NULL) // no excpetion
        result->Init();
      else
        domain->GetErsatzMaterial()->ExtractResults(result, isComplex_);
      break;

    case MECH_SHAPE:
      ExtractNodeOffset(result);
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


    default:
      Warning( "Resulttype not computable by mechanic PDE",
               __FILE__, __LINE__ );
    }
  }


  template <class TYPE>
  void MechPDE::CalcStresses( shared_ptr<BaseResult> res ) {

    //get the correct bilinear form
    Vector<Double> intPoint;

    //transform the type
    SubTensorType type;
    String2Enum(subType_,type);

    Vector<TYPE> elemStress;
    elemStress.Resize(stressDim_);
    elemStress.Init(0);

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
  void MechPDE::CalcStrains( shared_ptr<BaseResult> res ) {
    //transform the subType of the pde
    SubTensorType type;
    String2Enum(subType_,type);

    Vector<Double> intPoint;
    Vector<TYPE> elemStrain;
    elemStrain.Resize(stressDim_);
    elemStrain.Init(0);

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
    MechStressStrain<TYPE> *strain = new MechStressStrain<TYPE>(actSDMat,type);

    // loop over elements
    for ( it.Begin(); !it.IsEnd(); it++ ) {
      it.GetElem()->ptElem->GetCoordMidPoint(intPoint);

      //set element solution
      Matrix<TYPE> elSol;
      sol_->GetElemSolutionAsMatrix(elSol, it);
      strain->SetActElemSol(elSol);

      //calculates the element strain
      strain->SetIntPoint(intPoint);
      strain->CalcStrainVec(elemStrain,1,it);
      strain->UnsetIntPoint();
      for( UInt iDof = 0; iDof < stressDim_; iDof++ ) {
        actVal[it.GetPos()*stressDim_ + iDof] = elemStrain[iDof];
      }
    }

    delete strain;
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



  // ********************************************************
  //   Query parameter object for information about prestressing
  // ********************************************************
  void MechPDE::ReadPreStressing() {


    // Check, if any prestressing boundary condition is present
    ParamNode * bcsNode = myParam_->Get("bcsAndLoads", false );
    if( !bcsNode) return;

    // Get prestressing parameter nodes
    StdVector<ParamNode*> stressNodes = bcsNode->GetList("preStress");

    for (UInt k = 0; k < stressNodes.GetSize(); k++) {

      // get current name,
      std::string actRegionName, actType;

      stressNodes[k]->Get( "type", actType );
      if ( actType != "RHS" ) {
        EXCEPTION("Prestressing of type '" << actType
                  << "' is not supported at the moment.");
      }

      stressNodes[k]->Get( "region", actRegionName );
      RegionIdType actRegion = ptgrid_->RegionNameToId( actRegionName );

      // fetch data
      std::string type;
      Vector<Double> stress(6);
      stressNodes[k]->Get( "type", type );
      stressNodes[k]->Get( "stress1", stress[0] );
      stressNodes[k]->Get( "stress2", stress[1] );
      stressNodes[k]->Get( "stress3", stress[2] );
      stressNodes[k]->Get( "stress4", stress[3] );
      stressNodes[k]->Get( "stress5", stress[4] );
      stressNodes[k]->Get( "stress6", stress[5] );

      preStressList_[actRegion] = "RHS";
      preStressVal_[actRegion] = stress;

    }
  }

  void MechPDE::ReadSurfStress() {

    // try to get bcsAndLoads node
    ParamNode * bcNode = myParam_->Get("bcsAndLoads", false);
    if( !bcNode )
      return;
    StdVector<ParamNode*> stressNodes = bcNode->GetList("surfStress");


    // iterate over all surface stress definitions
    std::string surf, volume, phase;
    Matrix<Double> valMat;
    for( UInt i = 0; i < stressNodes.GetSize(); i++ ) {

      // fetch data
      stressNodes[i]->Get( "name", surf );
      stressNodes[i]->Get( "region", volume );
      stressNodes[i]->Get( "phase", phase );
      ParamTools::AsTensor<double>(stressNodes[i]->Get("value"),
                                         stressDim_, 1, valMat );

      // create new surface stress definition
      SurfStress actStress;

      actStress.surface = ptgrid_->RegionNameToId(surf);
      actStress.region = ptgrid_->RegionNameToId(volume);
      actStress.phase = phase[i];
      valMat.ConvertToVec_AppendRows( actStress.stress );

      // add surface stress definition
      RegionIdType regionId =  ptgrid_->RegionNameToId( surf );
      surfStresses_[regionId] = actStress;
    }

  }
  
  // ********************************************************
  //   Query parameter object for information about prestraining
  // ********************************************************
  void MechPDE::ReadPreStrainingFromXML(ParamNode* bcsNode, Vector<Double> &vals)
  {
    // Check, if any prestraining boundary condition is present
    if(!bcsNode) return;

    // Get prestraining parameter nodes
    StdVector<ParamNode*> strainNodes = bcsNode->GetList("testStrain");
    if(strainNodes.IsEmpty()) return;
    if(strainNodes.GetSize() > 1)
    {
      Warning("Cannot handle more than one test strain per excitation, ignoring the rest",
                   __FILE__, __LINE__ );
    }
    
    // fetch data
    // for every call to this function, we only look at the first test strain
    vals.Resize(6, 0.0);
    strainNodes[0]->Get("strain1", vals[0]);
    strainNodes[0]->Get("strain2", vals[1]);
    strainNodes[0]->Get("strain3", vals[2]);
    strainNodes[0]->Get("strain4", vals[3]);
    strainNodes[0]->Get("strain5", vals[4]);
    strainNodes[0]->Get("strain6", vals[5]);
  }
  
  void MechPDE::ReadPreStraining()
  {
    ReadPreStrainingFromXML(myParam_->Get("bcsAndLoads", false), preStrainVal_);
  }

  void MechPDE::ReadPressureLoads()
  {
    ReadPressureLoadsFromXML(myParam_->Get("bcsAndLoads", false), pressSurf_, pressVals_, pressPhase_);
  }
  
  void MechPDE::ReadPressureLoadsFromXML(ParamNode* bcNode, StdVector<shared_ptr<EntityList> >& pressSurf, StdVector<std::string>& pressVals, StdVector<std::string>& pressPhase) {
    if( !bcNode )
      return;
    StdVector<ParamNode*> presNodes = bcNode->GetList("pressure");

    // iterate over all pressure definitions
    std::string name, value, phase;
    for( UInt i = 0; i < presNodes.GetSize(); i++ ) {

      presNodes[i]->Get("name", name );
      presNodes[i]->Get("value", value );
      presNodes[i]->Get("phase", phase );

      pressSurf.Push_back( 
        ptgrid_->GetEntityList( EntityList::SURF_ELEM_LIST,
                                name, EntityList::REGION ) );
      pressVals.Push_back( value );
      pressPhase.Push_back( phase );

    }
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
  

} // end namespace CoupledField
