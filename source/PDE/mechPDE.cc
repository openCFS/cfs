// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "mechPDE.hh"

#include <sstream>
#include <iomanip>

#include "Forms/forms_header.hh"
#include "Forms/linElastInt.hh"
#include "Forms/FlatShellStiffInt.hh"
#include "Forms/FlatShellMassInt.hh"
#include "Forms/massInt.hh"
#include "Forms/linPressureInt.hh"
#include "Forms/singleEntryInt.hh"
#include "Forms/linSurfStressInt.hh"
#include "Driver/assemble.hh"
#include "newmark.hh"
#include "newmarkFracDampMech.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"
#include "Domain/ansatzFct.hh"
#include "Driver/stdSolveStep.hh"
#include "Utils/SmoothSpline.hh"

#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh" 
#endif

namespace CoupledField {

  MechPDE::MechPDE(Grid * aptgrid, ParamNode* paramNode )
    :SinglePDE( aptgrid, paramNode ) {
    ENTER_FCN( "MechPDE::MechPDE" );

    pdename_          = "mechanic";
    pdematerialclass_ = MECHANIC;

    fracDamping_   = false;
    effectiveMass_ = false;
    nonLin_        = false;
    nonLinMaterial_= false;


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

    else if ( subType_ == "flatShell" && probGeo == "3d" ) {
      stressDim_ = 3;
      Info->PrintF("", "=== FLAT SHELL PROBLEM\n");
    }
    else
      {
        EXCEPTION( "Subtype '" <<  subType_ << "' of PDE '"
                   <<  pdename_ <<  "' does not fit to problem  geometry '"
                   << probGeo << "'"; );
      }

     // timestepping formulation
    std::string str = "";
    myParam_->Get( "timeSteppingFormulation", str,  false );
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

    //check for prestressing
    //    ReadPreStressing();

    // Register functions for scripting
    RegisterFunctions();
  }

  MechPDE::~MechPDE()
  {
    ENTER_FCN( "MechPDE::~MechPDE" );

  }

  void MechPDE::ReadDampingInformation( )
  {
    ENTER_FCN( "MechPDE::ReadDampingInformation" );
      
    fracMemory_ = 0;
    bool identical = true; // i.e. same type of damping for all regions
    Integer firstFrac=-1;

    
    // try to get dampingList
    ParamNode * dampListNode = myParam_->Get( "dampingList", false );
    if( !dampListNode) return;
    
    // get specific damping nodes
    StdVector<ParamNode*> dampNodes = dampListNode->GetChildren();
    std::map<std::string, DampingType> idDampType;
    std::map<std::string, Double>      idDampFreq;
    std::map<std::string, Double>      idDampRatioDeltaF;
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
        materials_[actRegionId]->GetScalar(dampFreq,RAYLEIGH_FREQUENCY,REAL);
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
    ENTER_FCN( "MechPDE::ReadSpecialBCs" );
    
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
    ENTER_FCN( "MechPDE::DefineSprings" );
    
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
    for( UInt i = 0; i < springNodes.GetSize(); i++ ) {

      // get data from node
      springNodes[i]->Get( "name", name );
      springNodes[i]->Get( "dof", dofName );
      springNodes[i]->Get( "massValue", massVal );
      springNodes[i]->Get( "dampingValue", dampVal );
      springNodes[i]->Get( "stiffnessValue", stiffVal );

      UInt dof = results_[0]->GetDofIndex( dofName );

      
      shared_ptr<NodeList> spNode (new NodeList(ptgrid_) );
      spNode->SetNamedNodes( name );
      
      // stiffness value
      if( stiffVal > EPS ) {
        SingleEntryInt * stiffInt = 
          new SingleEntryInt( stiffVal, dof, dim_ );
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
          new SingleEntryInt( massVal, dof, dim_ );
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
          new SingleEntryInt( dampVal, dof, dim_ );
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
    ENTER_FCN( "MechPDE::ReadSoftening" );

    // Check if softeningList node is present
    ParamNode * softListNode = myParam_->Get("softeningList", false );
    if( !softListNode ) return;

    // Get child elements and read data
    StdVector<ParamNode*> softNodes = softListNode->GetChildren();
    std::map<std::string, std::string> idSoftTypeMap;
    std::string type, id;
    for( UInt i = 0; i < softNodes.GetSize(); i++ ) {
      type = softNodes[i]->GetName();
      softNodes[i]->Get( "id", id );
      idSoftTypeMap[id] = type;
    }
    
    if( softNodes.GetSize() ) {
      Info->PrintF( pdename_, "Applying softening for regions:\n" );
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
    ENTER_FCN( "MechPDE::InitNonLin");
    
    nonLin_ = false;

    // ==============================================================
    // NOTE: Currently we can only treat geometric non-linearity and
    //       we assume that for a mechanic PDE all regions either
    //       are linear or non-linear!
    // ==============================================================

    // Check, if "nonLinList" is present
    ParamNode * nonLinListNode = myParam_->Get("nonLinList", false );
    if( !nonLinListNode) 
      return;

    // Get nonlinear types
    StdVector<ParamNode*> nonLinNodes = nonLinListNode->GetChildren();
    for( UInt i = 0; i < nonLinNodes.GetSize(); i++ ) {
      
      std::string actTypeString = nonLinNodes[i]->GetName();
      std::string actId = nonLinNodes[i]->Get("id")->AsString();

      NonLinType actType;
      String2Enum( actTypeString, actType );
      nonLinIdType_[actId] = actType;

      // check type
      if( actType == GEOMETRIC ) {
        nonLin_ = TRUE;
      }

      if( actType == MATERIAL ) {
        nonLin_ = true;
        nonLinMaterial_ = true;
      }
    }
    
    // Run over all region and set entry in "regionNonLinId"
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
      regionNonLinType_[actRegionId] = nonLinIdType_[actNonLinId];
      
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
    ENTER_FCN( "MechPDE::DefineIntegerators" );

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
      
                  
    //voulme integrators

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
      ParamNode * actRegionNode = 
        myParam_->Get( "region", "name", actRegionName );
      

      //================= Check for Perfectly matchec layers ====================//
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
	actSDMat->GetScalar(density,DENSITY,REAL);

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
            actSDMat->GetScalar(beta,RAYLEIGH_BETA,REAL);
            actSDMat->GetScalar(measFreq,RAYLEIGH_FREQUENCY,REAL);
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
            
            DataType matType = IMAG; 
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
          actSDMat->GetScalar(density,DENSITY,REAL);
          MassInt * bilinearMass  = new MassInt(density, dim_, isaxi_);
          if ( diagMass_ ) {
            // diagonal mass matrix
            bilinearMass->SetDiagMass();
          }
          
          actIntDescr =  new BiLinFormContext(bilinearMass, MASS );
          
          
        } else {
          
          FlatShellMassInt * bilinearMass = new FlatShellMassInt(actSDMat);
          
          // Obtain thickness and penalty dof
          Double thickness = actRegionNode->Get("thickenss")->AsDouble();
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
          actSDMat->GetScalar(alpha,RAYLEIGH_ALPHA,REAL);
          actSDMat->GetScalar(measFreq,RAYLEIGH_FREQUENCY,REAL);
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


    // Define Integrators for composite materials
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
        myParam_->Get( "region", "name", actRegionName );      


      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );
 
        
      // Get penalty value for drilling dof of region
      Double penaltyDof = actRegionNode->Get("penaltyDof")->AsDouble();


      // ==============  add stiffness ===========================================

      // Now define integrator for this region
      FlatShellStiffInt * myInt = new FlatShellStiffInt(&composite);
      myInt->SetPenaltyDof( penaltyDof );
      BiLinFormContext * actIntDescrStiff =
        new BiLinFormContext( myInt, STIFFNESS);
      
      //check for damping
      if ( dampingList_[actRegion] == RAYLEIGH && complexMaterial == FALSE ) {
        Double beta, measFreq;
        std::string fac;
        actSDMat->GetScalar(beta,RAYLEIGH_BETA,REAL);
        actSDMat->GetScalar(measFreq,RAYLEIGH_FREQUENCY,REAL);
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
      FlatShellMassInt * bilinearMass = new FlatShellMassInt(&composite);
      bilinearMass->SetPenaltyDof( penaltyDof );
      BiLinFormContext * actIntDescrMass = new BiLinFormContext(bilinearMass, MASS);
            
      //check for damping
      if ( dampingList_[actRegion] == RAYLEIGH && complexMaterial == FALSE ) {
        Double alpha, measFreq;
        std::string fac;
        actSDMat->GetScalar(alpha,RAYLEIGH_ALPHA,REAL);
        actSDMat->GetScalar(measFreq,RAYLEIGH_FREQUENCY,REAL);
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
    for (UInt actSF = 0; actSF < pressSurf_.GetSize(); actSF++) {

      // create new surface element list
      shared_ptr<SurfElemList> actPressSurf( new SurfElemList(ptgrid_ ) );
      shared_ptr<ElemList> actPressElem( new ElemList(ptgrid_ ) );
      actPressSurf->SetRegion( pressSurf_[actSF] );
      actPressElem->SetRegion( pressSurf_[actSF] );


      LinearSurfForm * rhsSrcSurf = 
        new PressureLinForm(pressVals_[actSF], pressPhase_[actSF], isaxi_ );
      rhsSrcSurf->SetVoluInfo( materials_ );

      LinearFormContext * pressRhs = 
        new LinearFormContext( rhsSrcSurf );
      pressRhs->SetPtPde( this );
      pressRhs->SetResult( results_[0], actPressSurf );
      assemble_->AddLinearForm( pressRhs );
      
      // Give entities and result to equation numbering class
      // and solution class
      eqnMap_->AddResult( *results_[0], actPressElem );
    }
    
    
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

    // Add integrators for region loads
    MechVolForceInt * forceInt;
    std::map<RegionIdType, RegionLoad>::iterator loadIt = regionLoads_.begin();
    if (regionLoads_.size() != 0 ) {
      (*loadIt).second.Print(true, pdename_ );
    }
    for( loadIt = regionLoads_.begin(); loadIt != regionLoads_.end(); loadIt++ ) {
      forceInt = (*loadIt).second.GetIntegrator();

      // Create new element list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( loadIt->first );
      LinearFormContext * forceContext = 
        new LinearFormContext( forceInt );
      forceContext->SetPtPde(this);
      forceContext->SetResult( results_[0], actSDList );
      assemble_->AddLinearForm( forceContext );

      //assemble_->AddRhsSrcIntegrator( forceInt, (*loadIt).first,
      //                                (*loadIt).second.dynamics, nonLin_ );
      (*loadIt).second.Print(false, pdename_);

    }

    // Define Springs
    DefineSprings();

  }


  BaseForm *
  MechPDE::GetStiffIntegrator(BaseMaterial* actSDMat, 
                              RegionIdType regionId,
                              bool reducedInt)
  {
    ENTER_FCN( "MechPDE::GetStiffIntegrator" );
    
    // Get region name
    std::string regionName = ptgrid_->RegionIdToName( regionId );

    // Get region node
    ParamNode * actRegionNode = 
      myParam_->Get( "region", "name", regionName );
    
    BaseForm * bilinearStiff = NULL;
     
    // Check for FlatShellIntegrator, as this one is a special sub-type,
    // not handled by the SubTensorType
    if (subType_ == "flatShell" ) {
      FlatShellStiffInt * myInt = new FlatShellStiffInt(actSDMat);
      
       
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
  
  
  void MechPDE::ReadRegionLoads( ) {
    ENTER_FCN ( "MechPDE::ReadRegionLoads" );
    
    StdVector<std::string> names, dofs, refCoord, type, phase;
    StdVector<std::string> tempNames, tempDofs,  tempPhase;
    StdVector<std::string>  tempRefCoord, tempType;
    StdVector<RegionIdType> regionIds;
    StdVector<UInt> vecComp;
    StdVector<std::string> loadVec, tempLoadVec, tempLoad(dim_);
    UInt locDof = 0;
    Integer index = -1;
    
    // Check, if function was called by a scripting command
#ifdef USE_SCRIPTING
    if ( messenger->IsEvaluating() == true ) {
      
      // obtain parameters from messenger object
      // Note: If scripting is used, only one region load
      // can be specified per call
      SCRIPT_GET( std::string,name);
      SCRIPT_GET( std::string, value );
      SCRIPT_GET( std::string, dof );
      SCRIPT_GET( std::string, refCoordSys );
      SCRIPT_GET( std::string, type );
      
      // Copy single entries into vectors
      tempNames.Push_back( name );
      tempLoadVec.Push_back( value );
      tempDofs.Push_back( dof );
      tempRefCoord.Push_back( refCoordSys );
      tempType.Push_back( type );
      
    } else {
#endif
      // obtain parameters from ParamHandler
      // Note: Here all region loads are read (in contrast
      // when called by an external script)
      
      // try to get bcsAndLoads node
      ParamNode * bcNode = myParam_->Get("bcsAndLoads", false);
      if( !bcNode )
        return;
      StdVector<ParamNode*> loadNodes = bcNode->GetList("regionLoad");

     
      for( UInt i = 0; i < loadNodes.GetSize(); i++ ) {
        
        ParamNode * actNode = loadNodes[i];

        tempNames.Push_back( actNode->Get("name")->AsString() );
        tempLoadVec.Push_back( actNode->Get("value")->AsString() );
        tempPhase.Push_back( actNode->Get("phase")->AsString() );
        tempDofs.Push_back( actNode->Get("dof")->AsString() );
        tempRefCoord.Push_back( actNode->Get("refCoordSys")->AsString() );
        tempType.Push_back( actNode->Get("type")->AsString() );
      }

#ifdef USE_SCRIPTING
    }
#endif

    // --- Common part for scripting and parameter file ---

    
    // Now sort the names and remove double entries
    for (UInt i = 0; i < tempNames.GetSize(); i++) {
      index = names.Find(tempNames[i]);
      if ( index == -1) {
        names.Push_back(tempNames[i]);
      }
    }

    // Convert region names to ID - vector
    ptgrid_->RegionNameToId(regionIds, names);
    
    
    // loop over all regions
    for (UInt i = 0; i < names.GetSize(); i++) {

      // get for each name all related entries for value,
      // dof, refCoordSys and type
      loadVec.Clear();
      phase.Clear();
      dofs.Clear();
      refCoord.Clear();
      type.Clear();
      for (UInt iEntry = 0; iEntry < tempNames.GetSize(); iEntry++ ) {
        if ( names[i] == tempNames[iEntry] ) {
          loadVec.Push_back(tempLoadVec[iEntry]);
          phase.Push_back(tempPhase[iEntry]);
          dofs.Push_back(tempDofs[iEntry]);
          refCoord.Push_back(tempRefCoord[iEntry]);
          type.Push_back(tempType[iEntry]);
        }
      }
      
      // check if all entries for  refCoord and type
      // are the same
      for (UInt k=0; k<refCoord.GetSize(); k++) {
        if ( refCoord[k] != refCoord[0] ||
             type[k] != type[0] ) {
          EXCEPTION( "MechPDE::DefineRegionLoads: The region load on region "
                     << names[i] << " has not for all dofs the same entry for "
                     << "refCoord or type (total/unit)!" );
        }
      }
      
      // Check if an entry already exists for this region
      RegionLoad * curLoad;
      
      std::map<RegionIdType, RegionLoad>::iterator it;
      it = regionLoads_.find( regionIds[i] );

      if ( it == regionLoads_.end() ) {
        regionLoads_.insert( std::map<RegionIdType, RegionLoad>::value_type( regionIds[i], 
                                                                             RegionLoad( dim_, isaxi_ ) ) );
      }
      it = regionLoads_.find( regionIds[i] );
      curLoad = & (*it).second;
      
      // -- Fill in the data we have so far --
      curLoad->name = ptgrid_->RegionIdToName( regionIds[i] );
      curLoad->phase = phase[0];

      if ( curLoad->refCoord != std::string() 
           && curLoad->refCoord != refCoord[0] ) {
        EXCEPTION( "Inconsistent definition of time data for regionLoads" );
      } else {
        curLoad->refCoord = refCoord[0];
      }

      if ( curLoad->volume < EPS ) {
        curLoad->volume = ptgrid_->CalcVolumeOfRegion(regionIds[i], isaxi_);
      }
      
      // now create local load vector
      for (UInt iDim=0; iDim < loadVec.GetSize(); iDim++) {
        locDof = domain->GetCoordSystem(refCoord[iDim])->
          GetVecComponent(dofs[iDim]);
        curLoad->value[locDof-1] = loadVec[iDim];
        curLoad->type = type[iDim];
      }
    }
  }
  
  
  void MechPDE::DefineSolveStep()
  {
    ENTER_FCN( "MechPDE::DefineSolveStep" );

    solveStep_ = new StdSolveStep(*this);
  }




  // ======================================================
  // COUPLING SECTION
  // ======================================================


  void MechPDE::InitCoupling(PDECoupling * Coupling)
  {
    ENTER_FCN( "MechPDE::InitCoupling" );
  
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
    ENTER_FCN( "MechPDE::CalcOutputCoupling" );

    UInt dof = 0;
    SolutionType quantity;
    StdVector<UInt> * couplingnodes = NULL;
    StdVector<Elem*> * couplingElems = NULL;
    CFSVector * temp_values = NULL;
    Vector<Double> * values;
    StdVector<BaseMaterial*> * materials = NULL;
    StdVector<std::string> outputRegions;
  
    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == false)
      return;

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


  void MechPDE::CalcAcousticCouplingRHS( StdVector<Elem*> * couplingElems, 
                                         StdVector<BaseMaterial*> & materials,
                                         StdVector<UInt>& couplingNodes,
                                         Vector<Double> & elemCouplingSols,
                                         UInt couplingdof )
  {
    ENTER_FCN( "MechPDE::CalcAcousticCouplingRHS" );

    EXCEPTION( "Not working at the moment" );
    //     Matrix<Double> ptCoord, elemMat;
    //     Vector<Double> normal, sol, forceOnElem, nSol;
    //     Elem * ptVolElem;
    //     Double sign = 0.0;
    //     Double density = 0.0;
    //     Integer matIndex = -1;

    //     elemCouplingSols.Init(0.0);
  
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
    // 	materials[actElem]->GetScalar(density,DENSITY,REAL);
        
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
      
    //         for (UInt actNode=0; actNode<ptCoord.GetSizeRow(); actNode++)
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
    ENTER_FCN( "MechPDE::HasOutput" );

    if (output == MECH_DISPLACEMENT || output == MECH_VELOCITY || output == MECH_FORCE)
      return true;

    return false;
  }



  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================


  void MechPDE :: InitTimeStepping()
  {
    ENTER_FCN( "MechPDE::InitTimeStepping" );


    if ( fracDamping_ == false ) {
      if ( effectiveMass_ == true ) {
        if ( diagMass_ == true ) {
	  //explicit time stepping
	  TS_alg_ = new NewmarkEffMass( algsys_, true );
	}
	else {
	  TS_alg_ = new NewmarkEffMass( algsys_ );
	}
      }
      else if ( effectiveMass_ == false ) {
        TS_alg_ = new Newmark( algsys_ );
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
    ENTER_FCN( "MechPDE::ReadSpecialResults" );

     // ---------------------------
    //  Determine special results
    // ---------------------------
    ParamNode * bcNode = myParam_->Get("bcsAndLoads", false );
    if ( !bcNode ) return;
    ParamNode * volNode = 
      bcNode->Get("surfRegionResult", "type",  "volumeAboveDefSurf", false );
    if( !volNode ) return;
    StdVector<ParamNode*> volListNodes = 
      volNode->Get("surfRegionList")->GetList( "surfRegion" );
     
    if( volListNodes.GetSize() > 0 ) {
      
      shared_ptr<ResultInfo> vol(new ResultInfo);
      vol->resultType = MECH_DEF_VOLUME;
      vol->dofNames = "";
      vol->unit = "m^3";
      vol->definedOn = ResultInfo::REGION;
      vol->entryType = ResultInfo::SCALAR;
      vol->fctType = shared_ptr<ConstFct>(new ConstFct() );
      
      std::string quantity;
      Enum2String( MECH_DEF_VOLUME, quantity );
      Info->PrintF( pdename_, " Computing '%s' for elements:\n", quantity.c_str() );   
      shared_ptr<EntityList> actList;

      for( UInt i = 0; i < volListNodes.GetSize(); i++ ) {

        // get data from node
        std::string regionName, dirName;
        volListNodes[i]->Get( "name", regionName );
        volListNodes[i]->Get( "dof", dirName );

        Info->PrintF( pdename_, " %s\n", regionName.c_str() );
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
      }
    }
  }

  void MechPDE::DefineAvailResults() {
    ENTER_FCN( "MechPDE::DefineAvailResults" );
    
    // Check for subType
    StdVector<std::string> dispDofNames, stressDofNames;
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
      dispDofNames = "x", "y", "z", "tx", "ty", "tz";
      Warning( "Not yet implemented" );
    }

    // === MECHANIC DISPLACEMENT ===
    shared_ptr<ResultInfo> disp(new ResultInfo);
    disp->resultType = MECH_DISPLACEMENT;
    disp->dofNames = dispDofNames;
    disp->unit = "m";
    disp->entryType = ResultInfo::VECTOR;

    // check if problem is lagrange or legendre
    std::string approxType = myParam_->Get("type")->AsString();
    if ( approxType == "lagrange" ) {
      shared_ptr<AnsatzFct> fct(new LagrangeFct);
      disp->fctType = fct;
      disp->definedOn = ResultInfo::NODE;
    } else {
      
       // define Legendre type
       shared_ptr<LegendreFct> fct(new LegendreFct);
       if( myParam_->Get("isIsotropic")->AsBool() ) {
         UInt order = myParam_->Get("order")->AsUInt();
         fct->SetIsoOrder( order );
       } else {
         Matrix<UInt> orderMat;
         myParam_->GetDim1xDim2Tensor("anisotropic", dim_, 1, orderMat);
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

    // === MECHANIC ENERGY ===
    shared_ptr<ResultInfo> energy(new ResultInfo);    
    energy->resultType = MECH_ENERGY;
    energy->dofNames = "";
    energy->unit = "Ws";
    energy->entryType = ResultInfo::SCALAR;
    energy->definedOn = ResultInfo::REGION;
    energy->fctType = shared_ptr<ConstFct>(new ConstFct() );
    availResults_.insert( energy );
  }

  void MechPDE::CalcResults( shared_ptr<BaseResult> result ) {
    ENTER_FCN( "MechPDE::CalcResults" );
    
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

    default:
      Warning( "Resulttype not computable by mechanic PDE",
               __FILE__, __LINE__ );
    }
  }




  template <class TYPE>
  void MechPDE::CalcStresses( shared_ptr<BaseResult> res ) {
    ENTER_FCN( "MechPDE::CalcStresses" );
    
    //get the correct bilinear form
    Vector<Double> intPoint;
  
    //transform the type
    SubTensorType type;
    String2Enum(subType_,type);

    Vector<TYPE> elemStress, sortedStress;
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




 // ======================================================
  // POSTPROCESSING  
  // ======================================================
  template<class TYPE>
  void MechPDE::ComputeVolDefSurf( shared_ptr<BaseResult> vals ) {

    ENTER_FCN( "MechPDE::ComputeVolDefSurf" );
    
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

    ENTER_FCN( "MechPDE::ComputeVolElem" );

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

    ENTER_FCN( "MechPDE::ReadPreStressing" );

    // Check, if any prestressing boundary condition is present
    ParamNode * bcsNode = myParam_->Get("bcsAndLoads", false );
    if( !bcsNode) return;

    // Get prestressing parameter nodes
    StdVector<ParamNode*> stressNodes = bcsNode->GetList("preStress");
    
    for (UInt k = 0; k < stressNodes.GetSize(); k++) {

      // get current name, 
      std::string actRegionName;
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

    //     StdVector<std::string> regionNames;

    //     keyVec = pdename_, "preStressing", "preStress", "name";
    //     attrVec = "", "", "tag";
    //     valVec  = "", "", "anyTag";

    //     params->GetList(keyVec, attrVec, valVec, regionNames );
    //     ptgrid_->RegionNameToId( preStressDomain_, regionNames );

    //     if ( preStressDomain_.GetSize() > 0 ) {

    //       Info->PrintF( pdename_,
    //                     " Found prestressing in the following regions:\n" );

    //       Double tmpDir;

    //       // Construct vectors for restricted search parameter
    //       StdVector<std::string> keyVec;
    //       StdVector<std::string> attrVec;
    //       StdVector<std::string> valVec;
    //       attrVec = "", "", "name";

    //       // for each prestress domain ...
    //       for ( UInt k = 0; k < preStressDomain_.GetSize(); k++ ) {

    //         // ... read direction of magnetisation
    //         valVec = "", "", regionNames[k];

    //         keyVec  = pdename_, "preStressing", "preStress", "orientX";
    //         params->Get( keyVec, attrVec, valVec, tmpDir);
    //         preStressValX_.Push_back( tmpDir);

    //         keyVec  = pdename_, "preStressing", "preStress", "orientY";
    //         params->Get( keyVec, attrVec, valVec, tmpDir );
    //         preStressValY_.Push_back( tmpDir );

    //         keyVec  = pdename_, "preStressing", "preStress", "orientZ";
    //         params->Get( keyVec, attrVec, valVec, tmpDir );
    //         preStressValZ_.Push_back( tmpDir );

    //         // ... report name to logfile
    //         Info->PrintF( pdename_, "%s\n", regionNames[k].c_str());
    //       }

    //    }


  }

  void MechPDE::ReadSurfStress() {
    ENTER_FCN( "MechPDE::ReadSurfStress" );
    
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
      stressNodes[i]->GetDim1xDim2Tensor( "value", stressDim_, 
                                          1, valMat );
        
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

  void MechPDE::ReadPressureLoads() {
    ENTER_FCN( "MechPDE::ReadPressureLoads" );

    // try to get bcsAndLoads node
    ParamNode * bcNode = myParam_->Get("bcsAndLoads", false);
    if( !bcNode )
      return;
    StdVector<ParamNode*> presNodes = bcNode->GetList("pressure");
    
    // iterate over all pressure definitions
    std::string name, value, phase;
    for( UInt i = 0; i < presNodes.GetSize(); i++ ) {
      
      presNodes[i]->Get("name", name );
      presNodes[i]->Get("value", value );
      presNodes[i]->Get("phase", phase );
      
      pressSurf_.Push_back( ptgrid_->RegionNameToId( name ) );
      pressVals_.Push_back( value );
      pressPhase_.Push_back( phase );

    }
  }


  template <class TYPE>
  void MechPDE::CalcEnergy( shared_ptr<BaseResult> vals )
  {
    ENTER_FCN( "MechPDE::CalcEnergy" );

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
    a.Last().RegisterParam( "refCoordSys", ArgList::STRING );
    a.Last().RegisterParam( "type", ArgList::STRING );
    pt.Push_back( new FCPT ( this, &MechPDE::ReadRegionLoads ) );
    name.Push_back( "setRegionLoad" );

    
    // Now register all functions with scripting 
    for (UInt i = 0; i < pt.GetSize(); i++ ) {
      Script_RegisterFct(name[i], pt[i], a[i] );
    }          

  }


  // ========================== RegionLoads ==========================
  MechPDE::RegionLoad::RegionLoad( UInt dim, bool isAxi ) {
    
    value.Resize( dim );
    value.Init( "0.0");
    
    isAxi_ = isAxi;
    volume = 0.0;
    
  }

  MechVolForceInt * MechPDE::RegionLoad::GetIntegrator() {
    MechVolForceInt * forceInt = new MechVolForceInt( value.GetSize(),
                                                      phase,
                                                      isAxi_);

    // Check, if type is "unit"
    bool isUnit;
    if ( type == "total" ) {
      isUnit = false;
    } else {
      isUnit = true;
    }
    forceInt->SetVolForceVector( value,
                                 domain->GetCoordSystem( refCoord),
                                 isUnit, volume );
    
    return forceInt;
    
    
  }

  void MechPDE::RegionLoad::Print( bool onlyHeader, std::string pdeName ) {
    std::ostringstream out;

    if (onlyHeader) {
      Info->PrintF(pdeName, "The following regions have a region load:\n\n");
      out.clear();
      out << std::setw(15) << "name" << " | " 
          << std::setw(15) << "refCoordSys" << " | "
          << std::setw(11) << "volume" << " | "
          << std::setw(6) << "type" << " | "
          << std::setw(11) << "load" <<std::endl;
      Info->PrintF(pdeName, out.str().c_str());
      out.str("");
      out << std::setw(90) << std::setfill('-') << "" 
          << std::setfill(' ') << std::endl;
      Info->PrintF(pdeName, out.str().c_str());
      out.str("");
    } else {
  
        
      // write logging information into info file
      for (UInt k = 0; k < value.GetSize(); k++ ) {
        out.str("");
        if ( k == 0) {
          out << std::setw(15) << name << " | " 
              << std::setw(15) << refCoord << " | "
              << std::setw(11) << volume << " | "
              << std::setw(6) << type << "|";
        } else {
          out << std::setw(15) << "" << " | " 
              << std::setw(15) << "" << " | "
              << std::setw(15) << "" << " | "
              << std::setw(11) << "" << " | "
              << std::setw(6) << "" << " | ";
            
        }
          
        out << std::setw(11) << value[k] << std::endl;
          
        Info->PrintF(pdeName,out.str().c_str());
      }
        
    }
  }
  
} // end namespace CoupledField
