#include "mechPDE.hh"

#include <sstream>
#include <iomanip>

#include "Forms/forms_header.hh"
#include "Forms/linElastInt.hh"
#include "Forms/massInt.hh"
#include "Forms/linPressureInt.hh"
#include "DataInOut/writeresults.hh"
#include "Driver/assemble.hh"
#include "newmark.hh"
#include "newmarkFracDampMech.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "Driver/solveStepMech.hh" 
#include "CoupledPDE/pdecoupling.hh"
#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"

#ifdef TCL_INTERFACE
#include "DataInOut/Scripting/cfsmessenger.hh" 
#endif

namespace CoupledField
{

  MechPDE::MechPDE(Grid * aptgrid, TimeFunc *aptTimeFunc, WriteResults *aptOut)
    :SinglePDE(aptgrid, aptOut, aptTimeFunc), 
     preStressVal_(0.0)

  {
    ENTER_FCN( "MechPDE::MechPDE" );

    pdename_          = "mechanic";
    pdematerialclass_ = MECHANIC;
    fracDamping_ = FALSE;

    // ****************************
    // DETERMINE GEOMETRY
    // ****************************

    // Get problem geometry and PDE subtype
    params->Get( "subType", subType_, pdename_ );
    std::string probGeo;
    params->Get( "type", probGeo, "geometry" );

    // Set number of degrees of freedom and
    // ensure that subtype fits to problem geometry
    if ( subType_ == "3d" && probGeo == "3d" ) {
      dofspernode_ = 3;
      stressDim_ = 6;
      Info->PrintF("", "=== 3D PROBLEM\n");
    }
    else if ( subType_ == "axi" && probGeo == "axi" ) {
      isaxi_ = TRUE;
      dofspernode_ = 2;
      stressDim_ = 4;
      Info->PrintF("", "=== AXISYSMMETRIC PROBLEM\n");
    }
    else if ( subType_ == "planeStrain" && probGeo == "plane" ) {
      dofspernode_ = 2;
      stressDim_ = 3;
      Info->PrintF("", "=== PLANE STRAIN PROBLEM\n");
    }
    else if ( subType_ == "planeStress" && probGeo == "plane" ) {
      dofspernode_ = 2;
      stressDim_ = 3;
      Info->PrintF("", "=== PLANE STRESS PROBLEM\n");
    }
    else
      {
        std::string errmsg = "Subtype '" + subType_;
        errmsg += "' of PDE '" + pdename_ + "' does not fit to problem ";
        errmsg += "geometry '" + probGeo + "'\n";
        Info->Error( errmsg, __FILE__, __LINE__ );
      }

    // =====================================================================
    // set solution information
    // =====================================================================
    solTypes_ = MECH_DISPLACEMENT;
    solDofs_ = dofspernode_;
    
    effectiveMass_ = params->IsSet( "effMass" );
    
    // *********************************
    //  Check for pressure loads
    // *********************************

    //check for pressure loads
    StdVector<std::string> regionNames;

    params->GetList( "name"    , regionNames , pdename_, "pressure" );
    params->GetList( "value"   , pressVals_ , pdename_, "pressure" );
    params->GetList( "dynamics", pressFnc_  , pdename_, "pressure" );
    
    ptgrid_->RegionNameToId( pressSurf_, regionNames );
    // Check consistency
    if ( pressSurf_.GetSize() != pressVals_.GetSize() )
      {
        std::string errmsg = "PressureLoads: ";
        errmsg += "#name = " + GenStr(pressSurf_.GetSize());
        errmsg += ", #value = " + GenStr(pressVals_.GetSize());
        errmsg += ", #dynamics = " + pressFnc_.GetSize() + '\n';
        Info->Error( errmsg, __FILE__, __LINE__ );
      }
    // append pressure Surface to surface region of this PDE
    for ( UInt i = 0; i < pressSurf_.GetSize(); i++ ) {
      surfdoms_.Push_back( pressSurf_[i] );
    }
    
    // We need not have as many function/filenames as pressureloads!
    for ( UInt k = pressFnc_.GetSize(); k < pressSurf_.GetSize(); k++ )
      {
        pressFnc_.Push_back( "none" );
      }

    //check for prestressing
    ReadPreStressing();

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
    Boolean identical = TRUE; // i.e. same type of damping for all regions
    Integer firstFrac=-1;

    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    for (UInt k = 0; k < subdoms_.GetSize(); k++) {
      
      RegionIdType actRegion = subdoms_[k];

      std::string actRegionName;
      actRegionName = ptgrid_->RegionIdToName( actRegion );
      keyVec = "mechanic" , "region" , "damping" , "type";
      attrVec= ""         , "name"   , "";
      valVec = ""         , actRegionName, "";
      StdVector<std::string> dampInfo;
      params->GetList( keyVec, attrVec, valVec, dampInfo);
      
      if ( dampInfo.IsEmpty() ) {
        dampingList_[actRegion] = NONE;
        Info->PrintF( pdename_, 
                      "No information specifying damping detected!\n" );
      }
      else if (dampInfo[0] == "no") {
        dampingList_[actRegion] = NONE;
        Info->PrintF( pdename_, 
                      "      * NO damping at all for region: %s\n",
                      actRegionName.c_str() );
      }
      else if (dampInfo[0] == "rayleigh") {
        dampingList_[actRegion] = RAYLEIGH;
        Info->PrintF( pdename_, 
                      "      * RAYLEIGH damping for region: %s\n",
                      actRegionName.c_str() );
      }
      else if (dampInfo[0] == "fractional") {
        dampingList_[actRegion] = FRACTIONAL;
        fracDamping_ = TRUE;
        Info->PrintF( pdename_, 
                      "      * FRACTIONAL damping for region: %s\n", 
                      actRegionName.c_str() );
        
        // Find first region containing fractional damping
        if ( firstFrac < 0 )
          firstFrac = k;
        
        keyVec = "mechanic" , "region" , "damping" , "fracAlg";
        StdVector<std::string> fracAlg;
        params->GetList( keyVec, attrVec, valVec, fracAlg );
        
        keyVec = "mechanic" , "region" , "damping" , "fracMemory";
        StdVector<UInt> fracMem;
        params->GetList( keyVec, attrVec, valVec, fracMem );
        
        keyVec = "mechanic" , "region" , "damping" , "interpolation";
        StdVector<std::string> interpol;
        params->GetList(  keyVec, attrVec, valVec, interpol );
        
        
        if( fracAlg.IsEmpty() || fracMem.IsEmpty() || interpol.IsEmpty() ) {
          (*error) << "Specify attributes fracAlg, fracMemory " 
                   << "and interpolation!";
          Error( __FILE__, __LINE__ ); 
        }
        else if  ( fracAlg[0] == "gl" ) {
          
          // Include fracAlg and interpolation info in dampingList
          Info->PrintF( "", "\t\t\t using Gruenwald-Letnikov algorithm,\n");
          if (interpol[0] == "no" )
            dampingList_[actRegion] = FRACTIONAL_GL;
          else {
            Error("Till now no interpolation is allowed in mechanics fractional damnping!",__FILE__,__LINE__);
          }
        }
        //         else if ( fracAlg[0] == "blank" ) {
          
        //           Info->PrintF( "", "\t\t\t using Blanks algorithm,\n");
        //           if (interpol[0] == "no" )
        //             dampingList_[actRegion] = FRACTIONAL_BLANK;
        //           else {
        //             dampingList_[actREgion] = FRACTIONAL_BLANK_INT;
        //             Info->PrintF("", 
        //                          "\t\t\t linear interpol. of single past values\n\n");
        //           }
        //         }
        
        // up to now take maximum of fracMemory
        if ( fracMem[0] > fracMemory_ )
          fracMemory_ = fracMem[0];
      }
      
    }
    
    // Check, if all entries are identical
    for ( UInt i = 1; i < dampingList_.size(); i++ ) {
      if ( dampingList_[subdoms_[i-1]] != dampingList_[subdoms_[i]] ) {
        identical = FALSE;
        break;
      }
    }


    // Fractional damping can only be enabled, if all regions are damped
    // this way. Oterhwise an error is thrown.
    if ( fracDamping_ == TRUE ) {
      
      if ( identical == TRUE ) {
        
        fracDamping_ = TRUE;
        Info->PrintF(pdename_, "Memory size for fractional damping  is: %d\n",
                     fracMemory_ );
      } else {
        
        (*error) << "Fractional damping can only be used if it is applied for "
                 << "ALL regions of the mechanical domain!\n"
                 << "Please check your parameter file!";
        Error( __FILE__, __LINE__ );
      }
        
        
    }
  }
  
  void MechPDE::ReadSpecialBCs() {
    ENTER_FCN( "MechPDE::ReadSpecialBCs" );
    
    ReadRegionLoads();

  }

  
  void MechPDE::InitNonLin()
  {
    ENTER_FCN( "MechPDE::InitNonLin");

    // ==============================================================
    // NOTE: Currently we can only treat geometric non-linearity and
    //       we assume that for a mechanic PDE all regions either
    //       are linear or non-linear!
    // ==============================================================
    StdVector<std::string> nonLinRegion;
    params->GetList( "nonLinear", nonLinRegion, pdename_, "region" );
    // Should not happen with validating parser, but beware!
    if ( nonLinRegion.GetSize() == 0 ) {
      nonLin_ = FALSE;
    }
    else {
      for ( UInt k = 1; k < nonLinRegion.GetSize(); k++ ) {
        if ( nonLinRegion[k] != nonLinRegion[0] ) {
          Info->Error( "Non-linearity should be the same for all regions!",
                       __FILE__, __LINE__ );
        }
      }
      nonLin_ = nonLinRegion[0] == "geo" ? TRUE : FALSE;
    }

    // In non-linear case determine type of line search strategy
    if ( nonLin_ == TRUE ) {

      Info->PrintF( pdename_, "Non-linearity in %d regions\n",
                    nonLinRegion.GetSize() );

      // type of line search
      params->Get( "type", lineSearch_, pdename_, "lineSearch" );

      if ( lineSearch_ == "no" ) {
        Info->PrintF( pdename_, "Performing no line search\n" );
      }
      else {
        Info->PrintF( pdename_, "Will perform line search\n" );
      }

    }

    // If no non-linearity we do not perform line search anyhow
    else {
      lineSearch_ = "no";
    }

    if( nonLin_ == TRUE )
      {
        // incremental stopping criterion
        params->Get( "incStopCrit", incStopCrit_, pdename_, "nonLinear" );

        // residual stopping criterion
        params->Get( "resStopCrit", residualStopCrit_, pdename_, "nonLinear" );
        
        // maximal number of NL-iterations
        params->Get("maxNumIters", nonLinMaxIter_, pdename_, "nonLinear");
      }

    // ------------------------------------------
    //   Get information on reduced integration
    // ------------------------------------------
    params->GetList( "reducedInt", reducedIntegration_, pdename_, "region" );

    if ( nonLin_ == TRUE ) {
      for ( UInt i = 0; i < reducedIntegration_.GetSize(); i++ ) {
        if ( reducedIntegration_[i] == "yes" ) {
          (*error) << "Currently we do not support non-linearity with "
                   << "reduced integration!";
          Error( __FILE__, __LINE__ );
        }
      }
    }
  }
  

  void MechPDE::DefineIntegrators()
  {
    ENTER_FCN( "MechPDE::DefineIntegerators" );

    // help variables for parameter checking
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;


    // check for complex valued material parameter
    Boolean complexMaterial = 
      params->HasValue( "type", "imagMaterialParameter", "materialDataType" );

     //voulme integrators
    for (UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++)
      {
	BaseMaterial * actSDMat = materialData_[actSD];

	// ==============  add "standard" stiffness ===========================
	BaseForm * bilinearStiff = GetStiffIntegrator(materialData_[actSD]);
	
	//check  for softening!
	RegionIdType actRegion = subdoms_[actSD];
	std::string actRegionName;
	actRegionName = ptgrid_->RegionIdToName( actRegion );
	keyVec = "mechanic" , "region" , "softening" , "type";
	attrVec= ""         , "name"   , "";
	valVec = ""         , actRegionName, "";
	StdVector<std::string> softeningInfo;
	params->GetList( keyVec, attrVec, valVec, softeningInfo);
	
	if (softeningInfo.GetSize() > 0)
	  bilinearStiff->SetSofteningModel(softeningInfo[0]);
	
	IntegratorDescriptor * actIntDescrStiff =
	  new IntegratorDescriptor(bilinearStiff, STIFFNESS);
	
	actIntDescrStiff->SetPDEIds(this, this);
	
	//check for damping
	if ( dampingList_[subdoms_[actSD]] == RAYLEIGH && complexMaterial == FALSE ) {
	  Double beta;
	  actSDMat->GetScalar(beta,RAYLEIGH_BETA,REAL);
	  actIntDescrStiff->SetSecondaryMat(DAMPING, beta,analysistype_);
	}
	
	assemble_->AddIntegrator(actIntDescrStiff, subdoms_[actSD]);

	// check for complex valued material parameter
	if( complexMaterial ) {
	  BaseForm * bilinearStiffImag = GetStiffIntegrator(materialData_[actSD]);

	  DataType matType = IMAG; 
	  bilinearStiffImag->SetMatDataType(matType);

	  if (softeningInfo.GetSize() > 0)
	    bilinearStiffImag->SetSofteningModel(softeningInfo[0]);
	
	  IntegratorDescriptor * actIntDescrStiffImag =
	    new IntegratorDescriptor(bilinearStiffImag, STIFFNESS);
	
	  actIntDescrStiffImag->SetPDEIds(this, this);
	  actIntDescrStiffImag->SetMatDataType(matType);

	  assemble_->AddIntegrator(actIntDescrStiffImag, subdoms_[actSD]);
	}


	//for prestressing
	for ( UInt preStr=0; preStr<preStressDomain_.GetSize(); preStr++ ) {
	  if ( subdoms_[actSD] == preStressDomain_[preStr]) {
	    Vector<Double> preStrVal(3);
	    preStrVal[0] = preStressValX_[preStr];
	    preStrVal[1] = preStressValY_[preStr];
	    preStrVal[2] = preStressValZ_[preStr];
	    
	    BaseForm * bilinearPreStress;
	    if (subType_ == "planeStrain")
	      bilinearPreStress = new PreStressIntPlaneStrain(actSDMat, preStrVal);
	    else if (subType_ == "axi")
	      bilinearPreStress = new PreStressIntAxi(actSDMat, preStrVal);
	    else if (subType_ == "3d")
	      bilinearPreStress = new PreStressInt3D(actSDMat, preStrVal);
	    else 
	      Info->Error("Unknown subtype in mech PDE! ",__FILE__,__LINE__);               
	    
	    IntegratorDescriptor * actIntDescrPre =
	      new IntegratorDescriptor(bilinearPreStress, STIFFNESS);
	    actIntDescrPre->SetPDEIds(this, this);
	    
	    assemble_->AddIntegrator(actIntDescrPre, subdoms_[actSD]);
	  }
	}
       

	// ==============  add nonlinear stiffness ============================
	if (nonLin_)
	  {
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
	    
	    IntegratorDescriptor * stiffNL1Descr = 
	      new IntegratorDescriptor(nLinPart1, STIFFNESS, nonLin_);
	    
	    stiffNL1Descr->SetPDEIds(this, this);
	    assemble_->AddIntegrator(stiffNL1Descr, subdoms_[actSD]);
	    
	    IntegratorDescriptor * stiffNL2Descr = 
	      new IntegratorDescriptor(nLinPart2, STIFFNESS, nonLin_);
	    
	    stiffNL2Descr->SetPDEIds(this, this);
	    assemble_->AddIntegrator(stiffNL2Descr, subdoms_[actSD]);
	    
	    // assemble prestress, if in config-file given
	    //    if (preStressVal_)
	    //      AssemblePreStressMat(ptEl, connect_PDE, ptCoord,
	    //      actMatData, elDisp);
	  }
	
	
	// ==============  add mass ===========================================
	double density;
	actSDMat->GetScalar(density,DENSITY,REAL);
	BaseForm * bilinearMass  = new MassInt(density, dofspernode_, isaxi_);
	
	IntegratorDescriptor * actIntDescr =
	  new IntegratorDescriptor(bilinearMass, MASS);
	actIntDescr->SetPDEIds(this, this);
	
	//check for damping (mass part)
	
	if ( dampingList_[subdoms_[actSD]] == RAYLEIGH ) {
	  Double alpha;
	  actSDMat->GetScalar(alpha,RAYLEIGH_ALPHA,REAL);
	  actIntDescr->SetSecondaryMat(DAMPING, alpha, analysistype_);
	}
	
	assemble_->AddIntegrator(actIntDescr, subdoms_[actSD]);
	
	
	// ==================== RHS ===========================================
	if (nonLin_)
	  {
	    BaseForm * rhsSource = new nLinMech_linFormInt(actSDMat, isaxi_);
	    assemble_->AddRhsIntegrator(rhsSource, subdoms_[actSD], nonLin_);
	  }
      }
    
  
   
    //surface integrators
    //RHS-part
    Boolean nonlin = FALSE;
    for (UInt actSF = 0; actSF < pressSurf_.GetSize(); actSF++) {
      LinearSurfForm * rhsSrcSurf = new PressureLinForm(pressVals_[actSF], isaxi_);
      rhsSrcSurf->SetVoluInfo( subdoms_, materialData_ );
      assemble_->AddRhsSrcSurfIntegrator(rhsSrcSurf, pressSurf_[actSF], pressFnc_[actSF],
					 nonlin);
    }
    
    
    // Add integrators for region loads
    MechVolForceInt * forceInt;
    std::map<RegionIdType, RegionLoad>::iterator it = regionLoads_.begin();
    (*it).second.Print(TRUE, pdename_ );
    for( it = regionLoads_.begin(); it != regionLoads_.end(); it++ ) {
      forceInt = (*it).second.GetIntegrator();
      assemble_->AddRhsSrcIntegrator( forceInt, (*it).first,
                                      (*it).second.dynamics, nonLin_ );
      (*it).second.Print(FALSE, pdename_);

    }

  }


  BaseForm *
  MechPDE::GetStiffIntegrator(BaseMaterial* actSDMat, Boolean reducedInt)
  {

    ENTER_FCN( "MechPDE::GetStiffIntegrator" );

   //transform the type
    SubTensorType tensorType;
    String2Enum(subType_,tensorType);
  
    BaseForm * bilinearStiff = NULL;

    if( fracDamping_ == FALSE ) {
      bilinearStiff = new linElastInt(actSDMat, tensorType);
    } 
    else {
      StdVector<std::string> keyVec, attrVec, valVec;
      keyVec = "transient", "firstDt";
      attrVec = "tag";
      valVec  =  bcSequenceTag_;
      Double dt = 0.0;
      params->Get(keyVec,attrVec,valVec,dt);
      bilinearStiff = new LinViscoElastInt(actSDMat,tensorType, "modifiedStiffness",dt );
    }
    
    return bilinearStiff;
  }
  
  
  void MechPDE::ReadRegionLoads( ) {
    ENTER_FCN ( "MechPDE::ReadRegionLoads" );
    
    StdVector<std::string> keyVec, attrVec, valVec;
    StdVector<std::string> names, dofs, dynamics, refCoord, type;
    StdVector<std::string> tempNames, tempDofs, tempDynamics;
    StdVector<std::string>  tempRefCoord, tempType;
    StdVector<RegionIdType> regionIds;
    StdVector<UInt> vecComp;
    StdVector<std::string> loadVec, tempLoadVec, tempLoad(dim_);
    UInt locDof = 0;
    Integer index = -1;
    
    // Check, if function was called by a scripting command
#ifdef TCL_INTERFACE
    if ( messenger->IsEvaluating() == TRUE ) {
      
      // obtain parameters from messenger object
      // Note: If scripting is used, only one region load
      // can be specified per call
      SCRIPT_GET( std::string,name);
      SCRIPT_GET( std::string, value );
      SCRIPT_GET( std::string, dof );
      SCRIPT_GET( std::string, refCoordSys );
      SCRIPT_GET( std::string, type );
      SCRIPT_GET( std::string, dynamics );
      
      // Copy single entries into vectors
      tempNames.Push_back( name );
      tempLoadVec.Push_back( value );
      if ( dynamics == "" ) {
        tempDynamics.Push_back( "none" );
      } else {
        tempDynamics.Push_back( dynamics );
      }
      tempDofs.Push_back( dof );
      tempRefCoord.Push_back( refCoordSys );
      tempType.Push_back( type );
      
    } else {
#endif
      // obtain parameters from ParamHandler
      // Note: Here all region loads are read (in contrast
      // when called by an external script)
      
      // get names of all regions with region loads
      keyVec = "mechanic", "bcsAndLoads", "regionLoad", "name";
      attrVec = "", "tag", "";
      valVec  = "", bcSequenceTag_, "";
      params->GetList(keyVec,attrVec,valVec,tempNames);
      
      // get value
      keyVec = "mechanic", "bcsAndLoads", "regionLoad", "value";
      params->GetList(keyVec, attrVec, valVec, tempLoadVec);
     
      // get dynamics
      keyVec = "mechanic", "bcsAndLoads", "regionLoad", "dynamics";
      params->GetList(keyVec, attrVec, valVec, tempDynamics);
     
      // get dofs
      keyVec = "mechanic", "bcsAndLoads", "regionLoad", "dof";
      params->GetList(keyVec, attrVec, valVec, tempDofs);
     
      // get coordinate system
      keyVec = "mechanic", "bcsAndLoads", "regionLoad", "refCoordSys";
      params->GetList(keyVec, attrVec, valVec, tempRefCoord);
     
      // get load type (total / unit)
      keyVec = "mechanic", "bcsAndLoads", "regionLoad", "type";
      params->GetList(keyVec, attrVec, valVec, tempType);
#ifdef TCL_INTERFACE
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

      // get for each name all related entries for value, dynamics,
      // dof, refCoordSys and type
      loadVec.Clear();
      dynamics.Clear();
      dofs.Clear();
      refCoord.Clear();
      type.Clear();
      for (UInt iEntry = 0; iEntry < tempNames.GetSize(); iEntry++ ) {
        if ( names[i] == tempNames[iEntry] ) {
          loadVec.Push_back(tempLoadVec[iEntry]);
          dynamics.Push_back(tempDynamics[iEntry]);
          dofs.Push_back(tempDofs[iEntry]);
          refCoord.Push_back(tempRefCoord[iEntry]);
          type.Push_back(tempType[iEntry]);
        }
      }
      
      // check if all entries for dynamics, refCoord and type
      // are the same
      for (UInt k=0; k<dynamics.GetSize(); k++) {
        if (dynamics[k] != dynamics[0] ||
            refCoord[k] != refCoord[0] ||
            type[k] != type[0] ) {
          (*error) << "MechPDE::DefineRegionLoads: The region load on region "
                   << names[i] << " has not for all dofs the same entry for "
                   << "dynamics, refCoord or type (total/unit)!";
          Error( __FILE__, __LINE__ );
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
      if ( curLoad->dynamics != std::string() 
           && curLoad->dynamics != dynamics[0] ) {
        Error( "Inconsistent definition of time data for regionLoads",
               __FILE__, __LINE__ );
      } else {
        curLoad->dynamics = dynamics[0];
      }

      if ( curLoad->refCoord != std::string() 
           && curLoad->refCoord != refCoord[0] ) {
        Error( "Inconsistent definition of time data for regionLoads",
               __FILE__, __LINE__ );
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

    solveStep_ = new SolveStepMech(*this);
  }



  // ======================================================
  // COUPLING SECTION
  // ======================================================


  void MechPDE::InitCoupling(PDECoupling * Coupling)
  {
    ENTER_FCN( "MechPDE::InitCoupling" );
  
    isIterCoupled_ = TRUE;
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
            isIncrFormulation_ = TRUE;
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
    if (isIterCoupled_ == FALSE)
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
            Error("No Element coupling output", __FILE__,__LINE__);
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

    Matrix<Double> ptCoord, elemMat;
    Vector<Double> normal, sol, forceOnElem, nSol;
    Elem * ptVolElem;
    Double sign = 0.0;
    Double density = 0.0;
    Integer matIndex = -1;

    elemCouplingSols.Init(0.0);
  
    for (UInt actElem=0; actElem<couplingElems->GetSize(); actElem++)
      {
        // Perform cast from volume element to surface element, since
        // mech-acou coupling makes only sense on surface elements
        SurfElem * actCoupleElem = 
          dynamic_cast<SurfElem*> ((*couplingElems)[actElem]);
        

        BaseFE * ptElem = actCoupleElem->ptElem;
        StdVector<UInt> & connecth = (*couplingElems)[actElem]->connect;
        GetElemCoords(connecth, ptCoord);
      
        // Try to find according region for first neighbouring volume
        // element of the surface element
        matIndex = subdoms_.Find(actCoupleElem->ptVolElem1->regionId);
        
        // If first volume element does not belong to acoustic PDE, try the
        // second one
        if ( matIndex == -1 ) {
          matIndex = subdoms_.Find(actCoupleElem->ptVolElem2->regionId);
          ptVolElem = actCoupleElem->ptVolElem2;
          sign = -1.0 * actCoupleElem->normalSign;
        } else {
          ptVolElem = actCoupleElem->ptVolElem1;
          sign = 1.0 * actCoupleElem->normalSign;
        }
        
        if ( matIndex == -1) {
          (*error) << "MechPDE::CalcAcousticCouplingRHS: The two volume "
                   << "element neighbours of surface element Nr. "
                   << actCoupleElem->elemNum << " do not belong to my regions!";
          Error( __FILE__, __LINE__ );
        }
      
        // Assign correct density
	materials[actElem]->GetScalar(density,DENSITY,REAL);
        
        // get correct density belonging to the the neighbouring element
        // in the fluid subdomain
        //density = (*couplingMaterials)[actElem]->GetDensity();
      
        BaseForm * bilinear_mass = new MassInt(ptElem, density, isaxi_);
        Matrix<Double> elemmat;
        bilinear_mass->CalcElementMatrix(ptCoord, elemmat);
        delete bilinear_mass;       

        
        GetDerivSolVecOfElement(sol, connecth);
        nSol.Resize(connecth.GetSize());   // solution in normal direction

        // the normal vector points outwards of the mechanical domain
        // (see. Kaltenbacher, "Num. Sim. of Mech. Act. & Sens." chapter 8.2)
        ptgrid_->CalcSurfNormal( normal, *actCoupleElem );
        normal *= sign;
      

        for (UInt actNode=0; actNode < connecth.GetSize(); actNode++)
          for (UInt actDof=0; actDof<dofspernode_; actDof++)
            nSol[actNode] += sol[actDof + actNode*dofspernode_] * normal[actDof];


        Vector<Double> forceOnElem;
        forceOnElem= elemmat * nSol;  
      
        for (UInt actNode=0; actNode<ptCoord.GetSizeRow(); actNode++)
          {
            UInt nodePos = 0;
          
            while(connecth[actNode] != couplingNodes[nodePos] && nodePos < couplingNodes.GetSize()) 
              nodePos++;
            elemCouplingSols[nodePos] += forceOnElem[actNode];
          }      
      }
  } 



  Boolean MechPDE::HasOutput(SolutionType output)
  {
    ENTER_FCN( "MechPDE::HasOutput" );

    if (output == MECH_DISPLACEMENT || output == MECH_VELOCITY || output == MECH_FORCE)
      return TRUE;

    return FALSE;
  }



  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================


  void MechPDE :: InitTimeStepping()
  {
    ENTER_FCN( "MechPDE::InitTimeStepping" );

    UInt rhsSize = eqnData_->GetNumEQNs() *
      eqnData_->GetNumDofsPerEQN();

    if ( fracDamping_ == FALSE ) { 
      if (effectiveMass_)  
        TS_alg_ = new NewmarkEffMass( algsys_, rhsSize );
      else
        TS_alg_ = new Newmark( algsys_, rhsSize );
    }
    else {
      if (effectiveMass_ == FALSE)
        TS_alg_ = new NewmarkFracDampMech( algsys_, rhsSize, pdeId_,
                                           eqnData_, ptgrid_, this, 
                                           subdoms_, dampingList_ );
      else
        Error("This needs to be implemented!",__FILE__,__LINE__);
    }


  }



  // ======================================================
  // POSTPROCESSING SECTION
  // ======================================================


  void MechPDE::WriteResultsInFile(const UInt kstep,
                                   const Double asteptime,
                                   UInt stepOffset,
                                   Double timeOffset)
  {
    ENTER_FCN( "MechPDE::WriteResultsInFile" );

    NodeStoreSol<Double> sol_der1Array, sol_der2Array;
    NodeStoreSol<Double> * solTransient;
    NodeStoreSol<Complex> * solHarmonic;

    Double actTime = asteptime + timeOffset;
    UInt actStep = kstep + stepOffset;

#ifdef WRITE_RHS
    NodeStoreSol<Double> rhs;
    rhs.SetNumSolutions(1);
    rhs.SetNumNodes(numPDENodes_);
    rhs.SetSolutionType(ACOU_RHSVAL);
    rhs.SetNumDofs(dim_);
    rhs.SetPtrEQNData(eqnData_, ptgrid_);
    rhs.Init(0.0);
    
    Double *ptRHS;
    algsys_->GetRHSVal( ptRHS );
    rhs.CopyFromAlgSysDataPointer(ptRHS);
    outFile_->WriteNodeSolutionTransient(rhs, actStep, actTime);
#endif

    if ( isComplex_ == FALSE ) {
      solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);
    
      if (saveSol_ == TRUE ) 
        outFile_->WriteNodeSolutionTransient(*solTransient, actStep, actTime);
    
      if (analysistype_== TRANSIENT) {
        if (saveDeriv1_ == TRUE)
          {
            solDeriv1_.SetAlgSysVector(getS1());
            outFile_->WriteNodeSolutionTransient(solDeriv1_, actStep, actTime);
          }
      
        if (saveDeriv2_ == TRUE)
          {
            solDeriv2_.SetAlgSysVector(getS2());
            outFile_->WriteNodeSolutionTransient(solDeriv2_, actStep, actTime);
          }
      }
    
      //element results
      if (calcStress_.GetSize() !=0 ) {
        ElemStoreSol<Double> & stressConverted = 
          dynamic_cast<ElemStoreSol<Double>&>(*stress_);
        outFile_->WriteElemSolutionTransient( stressConverted, actStep, actTime);
      }
    } else {
      solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);
      
      if (saveSol_ == TRUE )
        outFile_->WriteNodeSolutionHarmonic(*solHarmonic,  actStep,
                                            actTime, complexFormat_);
    

      //element results
      if (calcStress_.GetSize() !=0 ) {
        ElemStoreSol<Complex> & stressConverted = 
          dynamic_cast<ElemStoreSol<Complex>&>(*stress_);
        outFile_->WriteElemSolutionHarmonic( stressConverted, actStep, 
                                             actTime, complexFormat_);
      }
    }
  }


  void MechPDE::WriteHistoryInFile(const UInt kstep,
                                   const Double asteptime,
                                   UInt stepOffset,
                                   Double timeOffset)
  {
    ENTER_FCN( "MechPDE::WriteHistoryInFile" );

    NodeStoreSol<Double> sol_der1Array, sol_der2Array;
    NodeStoreSol<Double> * solTransient;
    NodeStoreSol<Complex> * solHarmonic;

    Double actTime = asteptime + timeOffset;
    UInt actStep = kstep + stepOffset;
    
    if ( isComplex_ == FALSE ) {
      solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);
      
      if (saveSolHist_ == TRUE)
        outFile_->WriteNodeHistoryTransient(*solTransient, actStep, actTime);
      
      if (analysistype_== TRANSIENT) {
        if (saveDeriv1Hist_ == TRUE) {
          solDeriv1_.SetAlgSysVector(getS1());
          outFile_->WriteNodeHistoryTransient(solDeriv1_, actStep, actTime);
        }
      }
      
      if (saveDeriv2Hist_ == TRUE) {
        solDeriv2_.SetAlgSysVector(getS2());
        outFile_->WriteNodeHistoryTransient(solDeriv2_, actStep, actTime);
      }
    } else  {
      solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);

      if (saveSolHist_ == TRUE)
        outFile_->WriteNodeHistoryHarmonic(*solHarmonic,  actStep, 
                                           actTime, complexFormat_);
    
    } 
  }

  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void MechPDE::ReadStoreResults() {

    ENTER_FCN( "MechPDE::ReadStoreResults" );

    // Construct vectors for restricted parameter search
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    StdVector<std::string> regionNames;
    std::string quantity;
  
    // -------------------------
    //  Determine nodal results
    // -------------------------
    StdVector<std::string> nodeValues;
    keyVec  = pdename_, "storeResults", "nodeResults", "region";
    attrVec = "", "", "type";  

    // --- Mechanic Displacement ---
    Enum2String(MECH_DISPLACEMENT, quantity);
    valVec = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, nodeValues);
    if (nodeValues.GetSize() > 0) {
      saveSol_ = TRUE;
      hasOutput_ = TRUE;
    }


    // --- Mechanic Velocity ---
    Enum2String(MECH_VELOCITY, quantity);
    valVec = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, nodeValues);
    if (nodeValues.GetSize() > 0) {
      saveDeriv1_ = TRUE;
      hasOutput_ = TRUE;
    
      // intialize corresponding storesolution object
      solDeriv1_.SetNumSolutions(1);
      solDeriv1_.SetNumNodes(numPDENodes_);
      solDeriv1_.SetSolutionType(MECH_VELOCITY);
      solDeriv1_.SetNumDofs(dim_);
      solDeriv1_.SetPtrEQNData(eqnData_, ptgrid_); 
      solDeriv1_.Init();
    }

    // --- Mechanic Acceleration ---
    Enum2String(MECH_ACCELERATION, quantity);
    valVec = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, nodeValues);
    if (nodeValues.GetSize() > 0) {
      saveDeriv2_ = TRUE;
      hasOutput_ = TRUE;
      
      // intialize corresponding storesolution object
      solDeriv2_.SetNumSolutions(1);
      solDeriv2_.SetNumNodes(numPDENodes_);
      solDeriv2_.SetSolutionType(MECH_ACCELERATION);
      solDeriv2_.SetNumDofs(dim_);
      solDeriv2_.SetPtrEQNData(eqnData_, ptgrid_);
      solDeriv2_.Init();
    }

    // ---------------------------
    //  Determine element results
    // ---------------------------
    StdVector<std::string> elemResults;
    keyVec  = pdename_, "storeResults", "elemResults", "region";
    attrVec = "", "", "type";  

    // --- Mechanic Stress ---
    Enum2String(MECH_STRESS, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, regionNames );
    ptgrid_->RegionNameToId( calcStress_, regionNames );
    
    // If the symbolic name is "all" compute electric field for all regions
    if ( calcStress_.GetSize() == 1 && calcStress_[0] == ALL_REGIONS ) {
      calcStress_ = subdoms_;
    }

    // Log to info file
    if ( calcStress_.GetSize() > 0 ) {
      hasOutput_ = TRUE;
      Info->PrintF( pdename_,
                    " Computing mechanical stress for regions:\n");
      for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", regionNames[k].c_str() );
      }

      if ( isComplex_ == FALSE ) {
        stress_ = new ElemStoreSol<Double>;
      } else {
        stress_ = new ElemStoreSol<Complex>;
      }
      stress_->SetNumSolutions(1);
      stress_->SetSolutionType(MECH_STRESS);
      stress_->SetNumElems(numElems_);
      stress_->SetNumDofs(6);
      stress_->SetPtrEQNData(eqnData_, ptgrid_);
      stress_->Init();
    }

    // --- Mechanic Energy ---
    Enum2String(MECH_ENERGY, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, regionNames);
    ptgrid_->RegionNameToId( calcEnergy_, regionNames );

    // If the symbolic name is "all" compute electric field for all regions
    if ( calcEnergy_.GetSize() == 1 && calcEnergy_[0] == ALL_REGIONS ) {
      calcEnergy_ = subdoms_;
    }

    // Log to info file
    if ( calcEnergy_.GetSize() > 0 ) {
      hasOutput_ = TRUE;
      Info->PrintF( pdename_,
                    " Computing mechanical Energy for regions:\n");
      for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", regionNames[k].c_str() );
      }
    }

    // -------------------------
    //  Determine nodal history
    // -------------------------
    StdVector<std::string> saveNodeHist; 
    keyVec  = pdename_, "storeResults", "nodeHistory", "saveNodes";
    attrVec = "", "", "type";

    // --- mechDisplacement ---
    Enum2String(MECH_DISPLACEMENT, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, saveNodeHist );
  
    if (saveNodeHist.GetSize() > 0) {
      saveSolHist_ = TRUE;
      hasOutput_ = TRUE;
      Info->PrintF( pdename_, " Saving mechDisplacement for Nodes:\n" );
      for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", saveNodeHist[k].c_str() );
      }
    }
  
    // --- mechVelocity ---
    Enum2String(MECH_VELOCITY, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, saveNodeHist );
  
    if (saveNodeHist.GetSize() > 0) {
      saveDeriv1Hist_ = TRUE;
      hasOutput_ = TRUE;
      Info->PrintF( pdename_, " Saving mechVelocity for Nodes:\n" );
      for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", saveNodeHist[k].c_str() );
      }
      // Check if solDeriv_1 was already set by the previous nodal
      // output check
      if ( saveDeriv1_ != TRUE ) {   
        solDeriv1_.SetNumSolutions(1);
        solDeriv1_.SetNumNodes(numPDENodes_);
        solDeriv1_.SetSolutionType(MECH_VELOCITY);
        solDeriv1_.SetNumDofs(dim_);
        solDeriv1_.SetPtrEQNData(eqnData_, ptgrid_); 
        solDeriv1_.Init();
      }
    }

    // --- mechAcceleration ---
    Enum2String(MECH_ACCELERATION, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, saveNodeHist );
  
    if (saveNodeHist.GetSize() > 0) {
      saveDeriv2Hist_ = TRUE;
      hasOutput_ = TRUE;
      Info->PrintF( pdename_, " Saving mechAcceleration for Nodes:\n" );
      for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", saveNodeHist[k].c_str() );
      }

      // Check if solDeriv_1 was already set by the previous nodal
      // output check
      if ( saveDeriv2_ != TRUE ) {
        solDeriv2_.SetNumSolutions(1);
        solDeriv2_.SetNumNodes(numPDENodes_);
        solDeriv2_.SetSolutionType(MECH_ACCELERATION);
        solDeriv2_.SetNumDofs(dim_);
        solDeriv2_.SetPtrEQNData(eqnData_, ptgrid_);
        solDeriv2_.Init();
      }
    }

    // ---------------------------
    //  Determine element history
    // ---------------------------
    StdVector<std::string> saveElemHist;
    keyVec  = pdename_, "storeResults", "elemHistory", "saveElems";
    attrVec = "", "", "";
    valVec = "", "", "";
    params->GetList(keyVec, attrVec, valVec, saveElemHist);

    if (saveElemHist.GetSize() >   0) {
      std::string errMsg = pdename_;
      errMsg += ": Saving history elements is not implemented yet!\n";
      errMsg += "Meanwhile you can use 'unvtool' to extract element data.";
      Error( errMsg.c_str(), __FILE__, __LINE__);
    }

    // ---------------------------
    //  Determine special results
    // ---------------------------
    params->GetList( "name", regionNames, pdename_, "volumeAboveDefSurf" );
    ptgrid_->RegionNameToId( volAboveDefSurfRegions_, regionNames );
    params->GetList( "dof", volAboveDefSurfDir_, pdename_, 
		     "volumeAboveDefSurf" );

  }


  // ************************************************************
  //   PostProcess
  // ************************************************************

  void MechPDE::PostProcess() {

    ENTER_FCN( "MechPDE::PostProcess" );

    //check for mechanical energy calculation
    if (calcEnergy_.GetSize() !=0 ) {
      if ( isComplex_ == FALSE) {
        CalcEnergy<Double>();
      } else {
        CalcEnergy<Complex>();
      }
    }
  
    if (calcStress_.GetSize() !=0 ) {
      if ( isComplex_ == FALSE) {
        CalcStresses<Double>();
      } else {
        CalcStresses<Complex>();
      }
    }

    // Compute volume above deformed surfaces
    if ( volAboveDefSurfRegions_.GetSize() > 0 ) {
      ComputeVolDefSurf( volAboveDefSurfRegions_, volAboveDefSurfDir_ );
    }

    // Last but no least trigger setting of BC from script file
#ifdef TCL_INTERFACE
    StdVector<std::string> context;
    context.Push_back( pdename_ );
    context.Push_back( GenStr(solveStep_->GetActStep() ) );
    
    if ( analysistype_ == TRANSIENT ||
         analysistype_ == STATIC ) {
      context.Push_back( GenStr(solveStep_->GetActTime() ) );
    } else {
      context.Push_back( GenStr(solveStep_->GetActFreq() ) );
    }
    messenger->TriggerEvent( CFSMessenger::CFS_PostProcess, 
                             context );
#endif

  }

  template <class TYPE>
  void MechPDE::CalcStresses() {
    ENTER_FCN( "MechPDE::CalcStresses" );
    
    //get the correct bilinear form
    Vector<Double> intPoint;
    // Resize solution arrays
  
    //transform the type
    SubTensorType type;
    String2Enum(subType_,type);

    Vector<TYPE> elemStress, sortedStress;
    elemStress.Resize(stressDim_);
    elemStress.Init(0);
    sortedStress.Resize(6);
    
    // loop over all subdomains
    for (UInt isd=0; isd<subdoms_.GetSize(); isd++) {
      
      BaseMaterial* actSDMat = materialData_[isd];
      MechStressStrain<TYPE> *stress = new MechStressStrain<TYPE>(actSDMat,type);
      
      // get vector of Elements of subdomains
      StdVector<Elem*> elemssd;     
      ptgrid_->GetVolElems( elemssd,subdoms_[isd] );
        
      // loop over elements of subdomain
      for (UInt iel=0; iel< elemssd.GetSize(); iel++) {
        UInt pdeElem = eqnData_->Mesh2PDEElem(elemssd[iel]->elemNum);
        elemssd[iel]->ptElem->GetCoordMidPoint(intPoint);
        //set element pointer
        BaseFE * ptEl = elemssd[iel]->ptElem;
        stress->SetElemPtr(ptEl);
          
        //set element solution        
        Matrix<TYPE> elSol;
        StdVector<UInt> connecth = elemssd[iel]->connect;
        sol_->GetElemSolutionAsMatrix(elSol, connecth);
        stress->SetActElemSol(elSol);
          
        //get coordinates of element
        Matrix<Double> ptCoord;
        GetElemCoords(connecth, ptCoord);
          
        Vector<TYPE> actStress;     
          
        //set the integration point
        stress->SetIntPoint(intPoint);

        //calculates the stress
        stress->CalcStressVec(elemStress,1,ptCoord);
        sortStresses(elemStress,sortedStress);
        stress_->SetElemResult(pdeElem-1, sortedStress);
      }

      delete stress;
    }
  }
  
  // ********************************************************
  //   Query parameter object for information about prestressing
  // ********************************************************
  void MechPDE::ReadPreStressing() {

    ENTER_FCN( "MechPDE::ReadPreStressing" );

    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    StdVector<std::string> regionNames;

    keyVec = pdename_, "preStressing", "preStress", "name";
    attrVec = "", "", "tag";
    valVec  = "", "", bcSequenceTag_;

    params->GetList(keyVec, attrVec, valVec, regionNames );
    ptgrid_->RegionNameToId( preStressDomain_, regionNames );

    if ( preStressDomain_.GetSize() > 0 ) {

      Info->PrintF( pdename_,
                    " Found prestressing in the following regions:\n" );

      Double tmpDir;

      // Construct vectors for restricted search parameter
      StdVector<std::string> keyVec;
      StdVector<std::string> attrVec;
      StdVector<std::string> valVec;
      attrVec = "", "", "name";

      // for each prestress domain ...
      for ( UInt k = 0; k < preStressDomain_.GetSize(); k++ ) {

        // ... read direction of magnetisation
        valVec = "", "", regionNames[k];

        keyVec  = pdename_, "preStressing", "preStress", "orientX";
        params->Get( keyVec, attrVec, valVec, tmpDir);
        preStressValX_.Push_back( tmpDir);

        keyVec  = pdename_, "preStressing", "preStress", "orientY";
        params->Get( keyVec, attrVec, valVec, tmpDir );
        preStressValY_.Push_back( tmpDir );

        keyVec  = pdename_, "preStressing", "preStress", "orientZ";
        params->Get( keyVec, attrVec, valVec, tmpDir );
        preStressValZ_.Push_back( tmpDir );

        // ... report name to logfile
        Info->PrintF( pdename_, "%s\n", regionNames[k].c_str());
      }
    }
  }

  template <class TYPE>
  void MechPDE::CalcEnergy()
  {
    ENTER_FCN( "MechPDE::CalcEnergy" );

    Matrix<Double> elemmat;  
    Matrix<Double> ptCoord;
    BaseFE         * ptElem;

    StdVector<UInt> connecth;
    Vector<TYPE> help;

    TYPE totalE = 0;

    UInt i, j;
    Vector<TYPE> energy(subdoms_.GetSize());

    for (i=0; i<subdoms_.GetSize(); i++) {
    
      StdVector<Elem*> elemssd;
      ptgrid_->GetVolElems( elemssd,subdoms_[i] );

      //get material
      BaseMaterial* actSDMat = materialData_[i];
    
      energy[i] = 0;
      for (j=0; j < elemssd.GetSize(); j++) {  
        ptElem=elemssd[j]->ptElem;
	BaseForm * bilinear_stiff = GetStiffIntegrator(actSDMat);

        connecth=elemssd[j]->connect;
        GetElemCoords(connecth, ptCoord);
        bilinear_stiff->SetElemPtr(ptElem);
        bilinear_stiff->CalcElementMatrix(ptCoord, elemmat);

        Vector<TYPE> eldisp;
        sol_->GetElemSolution(eldisp, connecth);
        help = elemmat * eldisp;
        energy[i] += help * eldisp;

        delete bilinear_stiff;      

      }  

      totalE += energy[i];

    }

    std::string resulttype = "Mechanic Deformation Energy";
    std::string unit = "(Ws)";
    std::string analysis;
    Double analysisVal;
    if ( analysistype_ == HARMONIC ) {
      analysis    = "Frequency:";
      analysisVal = solveStep_->GetActFreq();
    }
    else {
      analysis    = "Time:";
      analysisVal = solveStep_->GetActTime();
    }
    
    StdVector<std::string> regionNames;
    ptgrid_->RegionIdToName( regionNames, subdoms_ );
    Info->WriteResult(pdename_,  resulttype, regionNames, energy, unit,
                      analysis, analysisVal);

    StdVector<std::string> suball(1);
    Vector<TYPE> tmp(1);
    suball[0] = "Sum";
    tmp[0] = totalE;
    Info->WriteResult(pdename_,  resulttype, suball, tmp, unit,
                      analysis, analysisVal);
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
    a.Last().RegisterParam( "dynamics", ArgList::STRING );
    pt.Push_back( new FCPT ( this, &MechPDE::ReadRegionLoads ) );
    name.Push_back( "setRegionLoad" );

    
    // Now register all functions with scripting 
    for (UInt i = 0; i < pt.GetSize(); i++ ) {
      Script_RegisterFct(name[i], pt[i], a[i] );
    }          

  }


  // ========================== RegionLoads ==========================
  MechPDE::RegionLoad::RegionLoad( UInt dim, Boolean isAxi ) {
    
    value.Resize( dim );
    value.Init( "0.0");
    
    isAxi_ = isAxi;
    volume = 0.0;
    
  }

  MechVolForceInt * MechPDE::RegionLoad::GetIntegrator() {
    MechVolForceInt * forceInt = new MechVolForceInt( value.GetSize(),
                                                      isAxi_);

    // Check, if type is "unit"
    Boolean isUnit;
    if ( type == "total" ) {
      isUnit = FALSE;
    } else {
      isUnit = TRUE;
    }
    forceInt->SetVolForceVector( value,
                                 domain->GetCoordSystem( refCoord),
                                 isUnit, volume );
    
    return forceInt;
    
    
  }

  void MechPDE::RegionLoad::Print( Boolean onlyHeader, std::string pdeName ) {
      std::ostringstream out;

      if (onlyHeader) {
        Info->PrintF(pdeName, "The following regions have a region load:\n\n");
        out.clear();
        out << std::setw(15) << "name" << " | " 
            << std::setw(15) << "refCoordSys" << " | "
            << std::setw(15) << "dynamics" << " | "
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
                << std::setw(15) << dynamics << " | "
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
