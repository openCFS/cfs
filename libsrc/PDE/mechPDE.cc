#include "mechPDE.hh"

#include <sstream>
#include <iomanip>

#include "Forms/forms_header.hh"
#include "Forms/linElastInt.hh"
#include "Forms/FlatShellInt.hh"
#include "Forms/massInt.hh"
#include "Forms/linPressureInt.hh"
#include "Forms/singleEntryInt.hh"
#include "DataInOut/writeresults.hh"
#include "Driver/assemble.hh"
#include "newmark.hh"
#include "newmarkFracDampMech.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "CoupledPDE/pdecoupling.hh"
#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"
#include "Domain/ansatzFct.hh"
#include "Driver/stdSolveStep.hh"

#ifdef TCL_INTERFACE
#include "DataInOut/Scripting/cfsmessenger.hh" 
#endif

namespace CoupledField {

  MechPDE::MechPDE(Grid * aptgrid, TimeFunc *aptTimeFunc, WriteResults *aptOut)
    :SinglePDE(aptgrid, aptOut, aptTimeFunc) {
    ENTER_FCN( "MechPDE::MechPDE" );

    pdename_          = "mechanic";
    pdematerialclass_ = MECHANIC;
    fracDamping_ = false;

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
        std::string errmsg = "Subtype '" + subType_;
        errmsg += "' of PDE '" + pdename_ + "' does not fit to problem ";
        errmsg += "geometry '" + probGeo + "'\n";
        Info->Error( errmsg, __FILE__, __LINE__ );
      }

    // =====================================================================
    // set solution information
    // =====================================================================
    
    // Create new resultDof object
    shared_ptr<ResultDof> res1(new ResultDof);
    shared_ptr<AnsatzFct> fct(new LagrangeFct);
    res1->resultType = MECH_DISPLACEMENT;
    if( dim_ == 3 ) {
      res1->dofNames = "ux", "uy", "uz";
    } else {
      res1->dofNames = "ux", "uy";
    }
    res1->definedOn = ResultDof::NODE;
    res1->fctType = fct;
    results_.Push_back( res1 );
    
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
        fracDamping_ = true;
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

    //check for prestressing
    ReadPreStressing();

  }


  void MechPDE::DefineSprings() {
    ENTER_FCN( "MechPDE::DefineSprings" );

    // Read information from xml file
    StdVector<std::string> keyVec, attrVec, valVec;
    StdVector<std::string> springNames, springDofs;
    StdVector<Double> stiffVals, dampVals, massVals;

    attrVec = "", "tag", "";
    valVec = "", bcSequenceTag_, "";
    
    keyVec = pdename_, "bcsAndLoads", "spring", "name";
    params->GetList(keyVec, attrVec, valVec, springNames);

    keyVec = pdename_, "bcsAndLoads", "spring", "dof";
    params->GetList(keyVec, attrVec, valVec, springDofs);

    keyVec = pdename_, "bcsAndLoads", "spring", "massValue";
    params->GetList(keyVec, attrVec, valVec, massVals);

    keyVec = pdename_, "bcsAndLoads", "spring", "dampingValue";
    params->GetList(keyVec, attrVec, valVec, dampVals);

    keyVec = pdename_, "bcsAndLoads", "spring", "stiffnessValue";
    params->GetList(keyVec, attrVec, valVec, stiffVals);
    

    // Iterate over all springs
    for( UInt iSp = 0; iSp < springNames.GetSize(); iSp++ ) {

      UInt dof = results_[0]->GetDofIndex( springDofs[iSp] );

      shared_ptr<NodeList> spNode (new NodeList(ptgrid_) );
      spNode->SetNodes( springNames[iSp] );
      
      // stiffness value
      if( stiffVals[iSp] > EPS ) {
        SingleEntryInt * stiffInt = 
          new SingleEntryInt( stiffVals[iSp], dof, dim_ );
        BiLinFormContext * stiffIntContext = 
          new BiLinFormContext( stiffInt, STIFFNESS );
        stiffIntContext->SetPtPdes(this, this);
        stiffIntContext->SetResults( results_[0], results_[0],
                                     spNode, spNode );
        assemble_->AddBiLinearForm( stiffIntContext );
      }

      // mass Value
      if ( massVals[iSp] > EPS ) {
        SingleEntryInt * massInt = 
          new SingleEntryInt( massVals[iSp], dof, dim_ );
        BiLinFormContext * massIntContext = 
          new BiLinFormContext( massInt, MASS );
        massIntContext->SetPtPdes(this, this);
        massIntContext->SetResults( results_[0], results_[0],
                                    spNode, spNode );
        assemble_->AddBiLinearForm( massIntContext );
      }
      
      // damping value
      if( dampVals[iSp] > EPS ) {
        SingleEntryInt * dampInt = 
          new SingleEntryInt( dampVals[iSp], dof, dim_ );
        BiLinFormContext * dampIntContext = 
          new BiLinFormContext( dampInt, DAMPING );
        dampIntContext->SetPtPdes(this, this);
        dampIntContext->SetResults( results_[0], results_[0],
                                    spNode, spNode );
        assemble_->AddBiLinearForm( dampIntContext );
      }
    }
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
      nonLin_ = false;
    }
    else {
      for ( UInt k = 1; k < nonLinRegion.GetSize(); k++ ) {
        if ( nonLinRegion[k] != nonLinRegion[0] ) {
          Info->Error( "Non-linearity should be the same for all regions!",
                       __FILE__, __LINE__ );
        }
      }
      nonLin_ = nonLinRegion[0] == "geo" ? true : false;
    }

    // In non-linear case determine type of line search strategy
    if ( nonLin_ == true ) {

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

    if( nonLin_ == true )
      {
        // incremental stopping criterion
        params->Get( "incStopCrit", incStopCrit_, pdename_, "nonLinear" );

        // residual stopping criterion
        params->Get( "resStopCrit", residualStopCrit_, pdename_, "nonLinear" );
        
        // maximal number of NL-iterations
        params->Get("maxNumIters", nonLinMaxIter_, pdename_, "nonLinear");

	// perform logging?
	nonLinLogging_ = params->IsSet( "logging", pdename_, "nonLinear" );

      }

    // ------------------------------------------
    //   Get information on reduced integration
    // ------------------------------------------
    params->GetList( "reducedInt", reducedIntegration_, pdename_, "region" );

    if ( nonLin_ == true ) {
      for ( UInt i = 0; i < reducedIntegration_.GetSize(); i++ ) {
        if ( reducedIntegration_[i] == "yes" ) {
          (*error) << "Currently we do not support non-linearity with "
                   << "reduced integration!";
          Error( __FILE__, __LINE__ );
        }
      }
    }
  }
  

  void MechPDE::DefineIntegrators() {
    ENTER_FCN( "MechPDE::DefineIntegerators" );

    // help variables for parameter checking
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    RegionIdType actRegion;
    BaseMaterial * actSDMat = NULL;
    
    
    // check for complex valued material parameter
    bool complexMaterial = 
      params->HasValue( "type", "imagMaterialParameter", "materialDataType" );
    
    //voulme integrators

    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {

      // Get current region and material
      actRegion = it->first;
      actSDMat = it->second;

      // create new entity list
      shared_ptr<ElemList> actSDList( new ElemList(ptgrid_ ) );
      actSDList->SetRegion( actRegion );

      
      // ==============  add "standard" stiffness ===========================
      BaseForm * bilinearStiff = GetStiffIntegrator(actSDMat, actRegion);
      
      //check  for softening!
      std::string actRegionName;
      actRegionName = ptgrid_->RegionIdToName( actRegion );
      keyVec = "mechanic" , "region" , "softening" , "type";
      attrVec= ""         , "name"   , "";
      valVec = ""         , actRegionName, "";
      StdVector<std::string> softeningInfo;
      params->GetList( keyVec, attrVec, valVec, softeningInfo);
      
      if (softeningInfo.GetSize() > 0)
        bilinearStiff->SetSofteningModel(softeningInfo[0]);
      
      BiLinFormContext * actIntDescrStiff =
        new BiLinFormContext(bilinearStiff, STIFFNESS );
      
      actIntDescrStiff->SetPtPdes(this, this);
      actIntDescrStiff->SetResults( results_[0], results_[0],
                                    actSDList, actSDList );


      //check for damping
      if ( dampingList_[actRegion] == RAYLEIGH && complexMaterial == false ) {
        Double beta;
        actSDMat->GetScalar(beta,RAYLEIGH_BETA,REAL);
        actIntDescrStiff->SetSecDestMat(DAMPING, beta );
      }
      
      assemble_->AddBiLinearForm( actIntDescrStiff );
        
      // check for complex valued material parameter
      if( complexMaterial ) {
        BaseForm * bilinearStiffImag = GetStiffIntegrator(actSDMat, 
                                                          actRegion);
          
        DataType matType = IMAG; 
        bilinearStiffImag->SetMatDataType(matType);
          
        if (softeningInfo.GetSize() > 0)
          bilinearStiffImag->SetSofteningModel(softeningInfo[0]);
          
        BiLinFormContext * actIntDescrStiffImag =
          new BiLinFormContext(bilinearStiffImag, STIFFNESS );
          
        actIntDescrStiffImag->SetPtPdes(this, this);
        actIntDescrStiffImag->SetResults( results_[0], results_[0],
                                          actSDList, actSDList );
        actIntDescrStiffImag->SetMatDataType(matType);
          
        assemble_->AddBiLinearForm(actIntDescrStiffImag  );
      }
        
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
      if (nonLin_) {
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
	
	
      // ==============  add mass ===========================================
      double density;
      actSDMat->GetScalar(density,DENSITY,REAL);
      BaseForm * bilinearMass  = new MassInt(density, dim_, isaxi_);
	
      BiLinFormContext * actIntDescr =
        new BiLinFormContext(bilinearMass, MASS );
      actIntDescr->SetPtPdes(this, this);
      actIntDescr->SetResults( results_[0], results_[0],
                               actSDList, actSDList );
      
      //check for damping (mass part)
	
      if ( dampingList_[actRegion] == RAYLEIGH ) {
        Double alpha;
        actSDMat->GetScalar(alpha,RAYLEIGH_ALPHA,REAL);
        actIntDescr->SetSecDestMat(DAMPING, alpha);
      }
	
      assemble_->AddBiLinearForm(actIntDescr);
	
	
      // ==================== RHS ===========================================
      if (nonLin_) {
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
    
  
   


    //surface integrators
    //RHS-part
    for (UInt actSF = 0; actSF < pressSurf_.GetSize(); actSF++) {

      // create new surface element list
      shared_ptr<SurfElemList> actPressSurf( new SurfElemList(ptgrid_ ) );
      actPressSurf->SetRegion( pressSurf_[actSF] );
   


      LinearSurfForm * rhsSrcSurf = new PressureLinForm(pressVals_[actSF], isaxi_);
      rhsSrcSurf->SetVoluInfo( materials_ );

      LinearFormContext * pressRhs = 
        new LinearFormContext( rhsSrcSurf, 0.0, pressFnc_[actSF] );
      pressRhs->SetPtPde( this );
      pressRhs->SetResult( results_[0], actPressSurf );
      assemble_->AddLinearForm( pressRhs );
      
      // Give entities and result to equation numbering class
      // and solution class
      eqnMap_->AddResult( *results_[0], actPressSurf );
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
        new LinearFormContext( forceInt, 0.0, (*loadIt).second.dynamics );
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
    
    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    // Get region name
    std::string regionName = ptgrid_->RegionIdToName( regionId );
    attrVec = "", "", "name";
    valVec = "", "", regionName;
    
    BaseForm * bilinearStiff = NULL;
     
    // Check for FlatShellIntegrator, as this one is a special sub-type,
    // not handled by the SubTensorType
    if (subType_ == "flatShell" ) {
      FlatShellInt * myInt = new FlatShellInt(actSDMat);
       
      // Get thickness of region
      Double thickness;
      keyVec = "pdeList", "mechanic", "region", "thickness";
      params->Get( keyVec, attrVec, valVec, thickness );
      myInt->SetThickness( thickness );
      std::cerr << "Thickness of region " << regionName << " is " << thickness << std::endl;
       
      // Get penalty value for drilling dof of region
      Double penaltyDof;
      keyVec = "pdeList", "mechanic", "region", "penaltyDof";
      params->Get( keyVec, attrVec, valVec, penaltyDof );
      myInt->SetPenaltyDof( penaltyDof );
      std::cerr << "PenaltyDof of region " << regionName << " is " << penaltyDof << std::endl;
       
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
        StdVector<std::string> keyVec, attrVec, valVec;
        keyVec = "transient", "firstDt";
        attrVec = "tag";
        valVec  =  bcSequenceTag_;
        Double dt = 0.0;
        params->Get(keyVec,attrVec,valVec,dt);
        bilinearStiff = new LinViscoElastInt(actSDMat,tensorType, "modifiedStiffness",dt );
      }
    } // if "flatShell"

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
    if ( messenger->IsEvaluating() == true ) {
      
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

    Error( "Not working at the moment", __FILE__, __LINE__ );
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
      if (effectiveMass_)  
        TS_alg_ = new NewmarkEffMass( algsys_ );
      else
        TS_alg_ = new Newmark( algsys_);
    }
    else {
      if (effectiveMass_ == false)
        TS_alg_ = new NewmarkFracDampMech( algsys_, pdeId_,
                                           eqnMap_, ptgrid_, this, 
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
    rhs.SetPtrEQNData(eqnMap_.get(), ptgrid_);
    rhs.Init(0.0);
    
    Double *ptRHS;
    algsys_->GetRHSVal( ptRHS );
    rhs.CopyFromAlgSysDataPointer(ptRHS);
    outFile_->WriteNodeSolutionTransient(rhs, actStep, actTime);
#endif

    if ( isComplex_ == false ) {

      solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);
    
      if (saveSol_ == true ) 
        outFile_->WriteNodeSolutionTransient(*solTransient, actStep, actTime);
    
      if (analysistype_== TRANSIENT) {
        if (saveDeriv1_ == true)
          {
            solDeriv1_.SetAlgSysVector(getS1());
            outFile_->WriteNodeSolutionTransient(solDeriv1_, actStep, actTime);
          }
      
        if (saveDeriv2_ == true)
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
      
      if (saveSol_ == true )
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
    
    if ( isComplex_ == false ) {
      solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);
      
      if (saveSolHist_ == true)
        outFile_->WriteNodeHistoryTransient(*solTransient, actStep, actTime);
      
      if (analysistype_== TRANSIENT) {
        if (saveDeriv1Hist_ == true) {
          solDeriv1_.SetAlgSysVector(getS1());
          outFile_->WriteNodeHistoryTransient(solDeriv1_, actStep, actTime);
        }
      }
      
      if (saveDeriv2Hist_ == true) {
        solDeriv2_.SetAlgSysVector(getS2());
        outFile_->WriteNodeHistoryTransient(solDeriv2_, actStep, actTime);
      }
    } else  {
      solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);

      if (saveSolHist_ == true)
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
      saveSol_ = true;
      hasOutput_ = true;
    }


    // --- Mechanic Velocity ---
    Enum2String(MECH_VELOCITY, quantity);
    valVec = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, nodeValues);
    if (nodeValues.GetSize() > 0) {
      saveDeriv1_ = true;
      hasOutput_ = true;
    
      // intialize corresponding storesolution object
      solDeriv1_.SetNumSolutions(1);
      solDeriv1_.SetNumNodes(numPDENodes_);
      solDeriv1_.SetSolutionType(MECH_VELOCITY);
      solDeriv1_.SetNumDofs(dim_);
      solDeriv1_.SetPtrEQNData( eqnMap_.get(),ptgrid_); 
      solDeriv1_.SetRegions( subdoms_ );
      solDeriv1_.Init();
    }

    // --- Mechanic Acceleration ---
    Enum2String(MECH_ACCELERATION, quantity);
    valVec = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, nodeValues);
    if (nodeValues.GetSize() > 0) {
      saveDeriv2_ = true;
      hasOutput_ = true;
      
      // intialize corresponding storesolution object
      solDeriv2_.SetNumSolutions(1);
      solDeriv2_.SetNumNodes(numPDENodes_);
      solDeriv2_.SetSolutionType(MECH_ACCELERATION);
      solDeriv2_.SetNumDofs(dim_);
      solDeriv2_.SetPtrEQNData( eqnMap_.get(), ptgrid_);
      solDeriv2_.SetRegions( subdoms_ );
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
      hasOutput_ = true;
      Info->PrintF( pdename_,
                    " Computing mechanical stress for regions:\n");
      for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", regionNames[k].c_str() );
      }

      if ( isComplex_ == false ) {
        stress_ = new ElemStoreSol<Double>;
      } else {
        stress_ = new ElemStoreSol<Complex>;
      }
      stress_->SetNumSolutions(1);
      stress_->SetSolutionType(MECH_STRESS);
      stress_->SetNumElems(numElems_);
      stress_->SetNumDofs(6);
      stress_->SetPtrEQNData( eqnMap_.get(), ptgrid_);
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
      hasOutput_ = true;
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
      saveSolHist_ = true;
      hasOutput_ = true;
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
      saveDeriv1Hist_ = true;
      hasOutput_ = true;
      Info->PrintF( pdename_, " Saving mechVelocity for Nodes:\n" );
      for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", saveNodeHist[k].c_str() );
      }
      // Check if solDeriv_1 was already set by the previous nodal
      // output check
      if ( saveDeriv1_ != true ) {   
        solDeriv1_.SetNumSolutions(1);
        solDeriv1_.SetNumNodes(numPDENodes_);
        solDeriv1_.SetSolutionType(MECH_VELOCITY);
        solDeriv1_.SetNumDofs(dim_);
        solDeriv1_.SetPtrEQNData( eqnMap_.get(),ptgrid_); 
        solDeriv1_.Init();
      }
    }

    // --- mechAcceleration ---
    Enum2String(MECH_ACCELERATION, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, saveNodeHist );
  
    if (saveNodeHist.GetSize() > 0) {
      saveDeriv2Hist_ = true;
      hasOutput_ = true;
      Info->PrintF( pdename_, " Saving mechAcceleration for Nodes:\n" );
      for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", saveNodeHist[k].c_str() );
      }

      // Check if solDeriv_1 was already set by the previous nodal
      // output check
      if ( saveDeriv2_ != true ) {
        solDeriv2_.SetNumSolutions(1);
        solDeriv2_.SetNumNodes(numPDENodes_);
        solDeriv2_.SetSolutionType(MECH_ACCELERATION);
        solDeriv2_.SetNumDofs(dim_);
        solDeriv2_.SetPtrEQNData( eqnMap_.get(), ptgrid_);
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
      if ( isComplex_ == false) {
        CalcEnergy<Double>();
      } else {
        CalcEnergy<Complex>();
      }
    }
  
    if (calcStress_.GetSize() !=0 ) {
      if ( isComplex_ == false) {
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
    sortedStress.Init();
    
    // loop over all subdomains
    for (UInt isd=0; isd<subdoms_.GetSize(); isd++) {
      
      BaseMaterial* actSDMat = materials_[subdoms_[isd]];
      MechStressStrain<TYPE> *stress = new MechStressStrain<TYPE>(actSDMat,type);
      

      ElemList actSDList(ptgrid_ );
      actSDList.SetRegion( subdoms_[isd] );
      
      EntityIterator it = actSDList.GetIterator();
      for ( it.Begin(); !it.IsEnd(); it++ ) {
        
        UInt pdeElem = eqnMap_->Mesh2PdeElem(it.GetElem()->elemNum);
        it.GetElem()->ptElem->GetCoordMidPoint(intPoint);
          
        //set element solution        
        Matrix<TYPE> elSol;
        StdVector<UInt> connecth = it.GetElem()->connect;
        sol_->GetElemSolutionAsMatrix(elSol, connecth);
        stress->SetActElemSol(elSol);
          
        Vector<TYPE> actStress;     
          
        //set the integration point
        stress->SetIntPoint(intPoint);

        //calculates the stress
        stress->CalcStressVec(elemStress,1,it);
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

    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    for (UInt k = 0; k < subdoms_.GetSize(); k++) {

      RegionIdType actRegion = subdoms_[k];

      std::string actRegionName;
      actRegionName = ptgrid_->RegionIdToName( actRegion );
      keyVec = "mechanic" , "region" , "preStress" , "type";
      attrVec= ""         , "name"   , "";
      valVec = ""         , actRegionName, "";
      StdVector<std::string> stressInfo;

      params->GetList( keyVec, attrVec, valVec, stressInfo);
      
  
      if ( !stressInfo.IsEmpty() ) {
	preStressList_[actRegion] = "RHS";
	Vector<Double> stress(3);
	StdVector<Double> getStress;
	
	keyVec = "mechanic" , "region" , "preStress" , "stress1";
	params->GetList( keyVec, attrVec, valVec, getStress );
	stress[0] = getStress[0];
	  
	keyVec = "mechanic" , "region" , "preStress" , "stress2";
	params->GetList( keyVec, attrVec, valVec, getStress );      
	stress[1] = getStress[0];
	  
	keyVec = "mechanic" , "region" , "preStress" , "stress3";
	params->GetList( keyVec, attrVec, valVec, getStress );      
	stress[2] = getStress[0];

	preStressVal_[actRegion] = stress;
      }
    }

//     StdVector<std::string> keyVec;
//     StdVector<std::string> attrVec;
//     StdVector<std::string> valVec;
//     StdVector<std::string> regionNames;

//     keyVec = pdename_, "preStressing", "preStress", "name";
//     attrVec = "", "", "tag";
//     valVec  = "", "", bcSequenceTag_;

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

  template <class TYPE>
  void MechPDE::CalcEnergy()
  {
    ENTER_FCN( "MechPDE::CalcEnergy" );

    Matrix<Double> elemmat;  

    StdVector<UInt> connecth;
    Vector<TYPE> help;

    TYPE totalE = 0;

    UInt i;
    Vector<TYPE> energy(subdoms_.GetSize());

    for (i=0; i<subdoms_.GetSize(); i++) {
    
      //get material
      BaseMaterial* actSDMat = materials_[subdoms_[i]];
      energy[i] = 0.0;

      ElemList actSDList(ptgrid_ );
      actSDList.SetRegion( subdoms_[i] );
      
      EntityIterator it = actSDList.GetIterator();
      for ( it.Begin(); !it.IsEnd(); it++ ) {
        
	BaseForm * bilinear_stiff = GetStiffIntegrator(actSDMat, subdoms_[i]);

        connecth=it.GetElem()->connect;
        bilinear_stiff->CalcElementMatrix(elemmat,it,it);

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
  MechPDE::RegionLoad::RegionLoad( UInt dim, bool isAxi ) {
    
    value.Resize( dim );
    value.Init( "0.0");
    
    isAxi_ = isAxi;
    volume = 0.0;
    
  }

  MechVolForceInt * MechPDE::RegionLoad::GetIntegrator() {
    MechVolForceInt * forceInt = new MechVolForceInt( value.GetSize(),
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
