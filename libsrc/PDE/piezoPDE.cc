#include <stdio.h>

#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/WriteInfo.hh"
#include "Driver/assemble.hh"
#include "ParamIdent/piezoParamIdent.hh"
#include "newmark.hh"
#include "Elements/basefe.hh"
#include "blocknodeEQN.hh"
#include "scalarblockEQN.hh"
#include "Forms/forms_header.hh"
#include "Driver/solveStepPiezo.hh"

#include "piezoPDE.hh" 

namespace CoupledField {

  PiezoPDE::PiezoPDE( Grid *aptgrid, TimeFunc *aptTimeFunc,
                      WriteResults *aptOut ) :
    SinglePDE( aptgrid, aptOut, aptTimeFunc ) {

    ENTER_FCN( "PiezoPDE::PiezoPDE" );

    // =====================================================================
    // set solution information
    // =====================================================================

    // Set name of the PDE and its material class
    pdename_ = "piezo";
    pdematerialclass_ = "piezo";
    piezoMaterialType_ = REALMATERIALPARAMETER; // default

    if( params->HasValue( "type", "imagMaterialParameter", pdename_,
                          "materialDataType" ) ) {

      piezoMaterialType_ = IMAGMATERIALPARAMETER; 
      Info->PrintF( pdename_, "Using complex piezoMaterialData\n" );
      std::cout << "\n++ Be aware that you are about to consider "
                << "complex-valued material parameter!"
                << std::endl;
    }

    // Get problem geometry and PDE subtype
    params->Get( "subType", subType_, pdename_ );
    std::string probGeo;
    params->Get( "type", probGeo, "geometry" );

    // Check for equilibration of material parameters
    std::string equilibrate;
    params->Get( "equilibrate", equilibrate, "piezo" );
    if ( equilibrate == "yes" ) {
      Info->PrintF( pdename_, "Using EQUILIBRATION of material parameters!\n");
      std::cout << "\n++ Using equilibration of material parameters! Not "
                << "corrected for in results!"
                << std::endl;
    }

    // Set number of degrees of freedom and
    // ensure that subtype fits to problem geometry
    if ( subType_ == "3d" && probGeo == "3d" ) {
      dofspernode_ = 4;
      Info->PrintF("", "=== 3D PROBLEM\n");
    }
    else if ( subType_ == "axi" && probGeo == "axi" ) {
      isaxi_ = TRUE;
      dofspernode_ = 3;
      Info->PrintF("", "=== AXISYSMMETRIC PROBLEM\n");
    }
    else if ( subType_ == "planeStrain" && probGeo == "plane" ) {
      dofspernode_ = 3;
      Info->PrintF("", "=== PLANE STRAIN PROBLEM\n");
    }
    else {
      std::string errmsg = "Subtype '" + subType_;
      errmsg += "' of PDE '" + pdename_ + "' does not fit to problem ";
      errmsg += "geometry '" + probGeo + "'\n";
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    isHysteresis_ = FALSE;
    params->GetList( "nonLinear", nonLinType_, pdename_, "region" );
    
    for ( UInt k = 0; k < nonLinType_.GetSize(); k++ ) {
      if ( nonLinType_[k] == "hysteresis" ) {
        isHysteresis_ = TRUE;
        break;
      }
    }

    if( isHysteresis_ ) {

      // solution method
      params->Get( "method", nonLinMethod_, pdename_, "nonLinear" );
      
      // perform logging?
      nonLinLogging_ = params->IsSet( "logging", pdename_, "nonLinear" );

      // type of line search
      params->Get( "type", lineSearch_, pdename_, "lineSearch" );

      // incremental stopping criterion
      params->Get( "incStopCrit", incStopCrit_, pdename_, "nonLinear" );
      
      // residual stopping criterion
      params->Get( "resStopCrit", residualStopCrit_, pdename_, "nonLinear" );

      // maximal number of NL-iterations
      params->Get("maxNumIters", nonLinMaxIter_, pdename_, "nonLinear");
    }
 

    // =====================================================================
    // set solution information
    // =====================================================================
    solTypes_ = MECH_DISPLACEMENT, ELEC_POTENTIAL;
    solDofs_ = dofspernode_ - 1, 1;

    // Use effective mass approach?
    effectiveMass_ = params->IsSet( "effMass" );
    
    // check for pressure loads
    StdVector<std::string> regionNames;
    
    params->GetList( "name"    , regionNames, pdename_, "pressure" );
    params->GetList( "value"   , pressVals_ , pdename_, "pressure" );
    params->GetList( "phase"   , pressPhase_ , pdename_, "pressure" );
    params->GetList( "dynamics", pressFnc_  , pdename_, "pressure" );

    ptgrid_->RegionNameToId( pressSurf_, regionNames );
    
    // Check consistency
    if ( pressSurf_.GetSize() != pressVals_.GetSize() ) {
      std::string errmsg = "PressureLoads: ";
      errmsg += "#name = " + Info->GenStr(pressSurf_.GetSize());
      errmsg += ", #value = " + Info->GenStr(pressVals_.GetSize());
      errmsg += ", #dynamics = " + pressFnc_.GetSize() + '\n';
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    // We need not have as many function/filenames as pressureloads!
    for ( UInt k = pressFnc_.GetSize(); k < pressSurf_.GetSize(); k++ )
      {
        pressFnc_.Push_back( "none" );
      }
  }

  void PiezoPDE::ReadDampingInformation( )
  {
    ENTER_FCN( "PiezoPDE::ReadDampingInformation" );
    
    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    for (UInt k = 0; k < subdoms_.GetSize(); k++) {
      
      RegionIdType actRegion = subdoms_[k];

      std::string actRegionName;
      actRegionName = ptgrid_->RegionIdToName( actRegion );


      actRegionName = ptgrid_->RegionIdToName( actRegion );
      keyVec = "piezo" , "region" , "damping" , "type";
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


    }
    
  }

  // *********************
  //   DefineIntegrators
  // *********************
  void PiezoPDE::DefineIntegrators() {
    ENTER_FCN( "PiezoPDE::DefineIntegerators" );

    piezoMaterialType realMatParameter = REALMATERIALPARAMETER; 

    UInt posOfElectricPot;
    if ( subType_ == "3d" )
      posOfElectricPot = 4;
    else
      posOfElectricPot = 3;

    for ( UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++ ) {

      // ==============  add stiffness ======================================

      MaterialData actSDMat(materialData_[actSD]);

      // ==============  add "standard" stiffness ===========================

      BaseForm * bilinearStiff = GetStiffIntegrator(actSDMat);
      IntegratorDescriptor *actIntDescrStiff = 
        new IntegratorDescriptor(bilinearStiff, STIFFNESS);
      actIntDescrStiff->SetPDEIds(this, this);

      bilinearStiff->SetPiezoMaterialType(realMatParameter);
      actIntDescrStiff->SetPiezoMaterialType(realMatParameter);
      assemble_->AddIntegrator(actIntDescrStiff, subdoms_[actSD]);
      
      // check for complex-valued material parameter
      if (piezoMaterialType_ == IMAGMATERIALPARAMETER){
        BaseForm * bilinearStiffC = GetStiffIntegrator(actSDMat);
        IntegratorDescriptor *actComplexIntDescrStiff = 
          new IntegratorDescriptor(bilinearStiffC, STIFFNESS);
        actComplexIntDescrStiff->SetPDEIds(this, this);

        actComplexIntDescrStiff->SetPiezoMaterialType(piezoMaterialType_);
        bilinearStiffC->SetPiezoMaterialType(piezoMaterialType_);
        assemble_->AddIntegrator(actComplexIntDescrStiff, subdoms_[actSD]);
      }
      
      //check for damping
      if ( dampingList_[subdoms_[actSD]] == RAYLEIGH ) {
        Boolean isdamping = TRUE;
        Boolean reducedIntegration = FALSE; //is currently not supported
        BaseForm * dampStiff = 
          GetStiffIntegrator( actSDMat, reducedIntegration,isdamping );
        dampStiff->SetRaylDamping();

        IntegratorDescriptor *actIntDescrDamp =
          new IntegratorDescriptor(dampStiff, DAMPING);
        actIntDescrDamp->SetPDEIds(this, this);

        dampStiff->SetPiezoMaterialType(realMatParameter);
        actIntDescrDamp->SetPiezoMaterialType(realMatParameter);
        
        assemble_->AddIntegrator(actIntDescrDamp, subdoms_[actSD]);

        // check for complex-valued material parameter
        if (piezoMaterialType_ == IMAGMATERIALPARAMETER){
          Boolean isdamping = TRUE;
          Boolean reducedIntegration = FALSE; //is currently not supported
          BaseForm * dampStiffC = 
            GetStiffIntegrator( actSDMat,reducedIntegration, isdamping );
          IntegratorDescriptor *actComplexIntDescrDamp =
            new IntegratorDescriptor(dampStiffC, DAMPING);
          actComplexIntDescrDamp->SetPDEIds(this, this);

          dampStiffC->SetPiezoMaterialType(piezoMaterialType_);
          actComplexIntDescrDamp->SetPiezoMaterialType(piezoMaterialType_);
          assemble_->AddIntegrator(actComplexIntDescrDamp, subdoms_[actSD]);
        }
      }


      // ==============  add mass =============================================

      Double density = actSDMat.GetDensity();    
      BaseForm * bilinearMass  = new MassInt(density, dofspernode_,
                                             posOfElectricPot, isaxi_);
      bilinearMass->SetPiezoMaterialType(realMatParameter);

      IntegratorDescriptor * actIntDescrMass =
        new IntegratorDescriptor(bilinearMass, MASS);
      actIntDescrMass->SetPDEIds(this, this);

      //check for damping (mass part)
      if (dampingList_[subdoms_[actSD]] == RAYLEIGH ){    
        actIntDescrMass->SetSecondaryMat(DAMPING, actSDMat.GetDampingAlfa(),
                                         analysistype_);
      }
      actIntDescrMass->SetPiezoMaterialType(realMatParameter);
      assemble_->AddIntegrator(actIntDescrMass, subdoms_[actSD]);
      
      
    } // end for actSD ...
    

    //surface integrators
    //RHS-part
    if (analysistype_== HARMONIC || analysistype_==MULTIHARMONIC){
      UInt nonlin = 0;
      Vector<Complex> pressValsC_(pressVals_.GetSize());
            
      for (UInt actSF = 0; actSF < pressSurf_.GetSize(); actSF++) {
        BaseForm * rhsSrcSurf = new 
          PressureLinForm(pressVals_[actSF], isaxi_);
        rhsSrcSurf->SetDofZero(posOfElectricPot);
	
        assemble_->AddRhsSrcSurfIntegrator(rhsSrcSurf, pressSurf_[actSF], 
                                           pressFnc_[actSF],nonlin);
        assemble_->AddRhsSrcSurfIntegrator(rhsSrcSurf, pressSurf_[actSF], 
                                           pressPhase_[actSF],nonlin);
      }
	
    }
    else
      {

        Boolean nonlin = FALSE;
        for (UInt actSF = 0; actSF < pressSurf_.GetSize(); actSF++) {
          BaseForm * rhsSrcSurf = new 
            PressureLinForm(pressVals_[actSF], isaxi_);
          rhsSrcSurf->SetDofZero(posOfElectricPot);
          assemble_->AddRhsSrcSurfIntegrator(rhsSrcSurf, pressSurf_[actSF], 
                                             pressFnc_[actSF],nonlin);
        }
      }
  }


  // **********************
  //   GetStiffIntegrator
  // **********************
  BaseForm *PiezoPDE::GetStiffIntegrator( MaterialData& actSDMat,
                                          Boolean reducedInt,
                                          Boolean isdamping ) {

    ENTER_FCN( "PiezoPDE::GetStiffIntegrator" );
  
    BaseForm * bilinearStiff;
    if (subType_ == "planeStrain")
      bilinearStiff = new piezoPlainStrainInt(actSDMat,isdamping);

    else if (subType_ == "axi")      
      bilinearStiff = new piezoAxiInt(actSDMat,isdamping);
    
    else if (subType_ == "3d")
      bilinearStiff = new linPiezo3DInt(actSDMat,isdamping);

    else
      Error("Unknown subtype in mech PDE! ",__FILE__,__LINE__);

    return bilinearStiff;
  }

  void PiezoPDE::DefineSolveStep()
  {
    ENTER_FCN( "PiezoPDE::DefineSolveStep" );
    
    if ( isHysteresis_ ) {
      solveStep_ = new  SolveStepPiezo(*this);
    }
    else {
      solveStep_ = new  StdSolveStep(*this);
    }
  }


  // ======================================================
  // TRANSIENT SOLVING SECTION
  // ======================================================


  void PiezoPDE::InitTimeStepping() {
    ENTER_FCN( "PiezoPDE::InitTimeStepping" );
    UInt rhsSize = eqnData_->GetNumEQNs() *
      eqnData_->GetNumDofsPerEQN();

    if (effectiveMass_)
      TS_alg_ = new NewmarkEffMass( algsys_, rhsSize );
    else
      TS_alg_ = new Newmark( algsys_, rhsSize );
  }


  void PiezoPDE::WriteResultsInFile( const UInt kstep,
                                     const Double asteptime,
                                     UInt stepOffset,
                                     Double timeOffset ) {

    ENTER_FCN( "PiezoPDE::WriteResultsInFile" );

    NodeStoreSol<Double> * solTransient;
    NodeStoreSol<Complex> * solHarmonic;

    //     Double actTime = lasttimecalc_ + timeOffset;
    //     UInt actStep = laststepcalc_ + stepOffset;

    Double actTime = asteptime + timeOffset;
    UInt actStep = kstep + stepOffset;

    if (analysistype_ == STATIC ||
        analysistype_ == TRANSIENT) {
      solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);      

      if (saveSol_)
        outFile_->WriteNodeSolutionTransient(*solTransient, actStep, actTime);

      // Write derivatives
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
      if (calcEfield_.GetSize() !=0 ) 
        outFile_->WriteElemSolutionTransient(Efield_, actStep, actTime);
        
        
      if (calcStress_.GetSize() !=0 ) 
        outFile_->WriteElemSolutionTransient(stress_, actStep, actTime);
        
        
      //element results
      if (calcCharge_.GetSize() !=0 ) {
        outFile_->WriteElemSolutionTransient(charges_, actStep, actTime);
      } 
        
    } else if (analysistype_ == HARMONIC || analysistype_==MULTIHARMONIC) {
      
      if (saveSol_)
        {
          solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);
          outFile_->WriteNodeSolutionHarmonic(*solHarmonic,  actStep, 
                                              actTime, complexFormat_);
        }

      //element results
      if (calcCharge_.GetSize() !=0 ) {
        outFile_->WriteElemSolutionHarmonic(chargesComplex_, actStep,  
                                            actTime, complexFormat_);
      }
      if (calcStress_.GetSize() !=0 ) 
        outFile_->WriteElemSolutionHarmonic(stressComplex_, actStep,
                                            actTime, complexFormat_);
      if (calcEfield_.GetSize() !=0 ) 
        outFile_->WriteElemSolutionHarmonic(EfieldComplex_, actStep,
                                            actTime, complexFormat_);
    }
  }


  void PiezoPDE::WriteHistoryInFile( const UInt kstep,
                                     const Double asteptime,
                                     UInt stepOffset,
                                     Double timeOffset ) {

    ENTER_FCN( "PiezoPDE::WriteHIstoryInFile" );

    NodeStoreSol<Double> * solTransient;
    NodeStoreSol<Complex> * solHarmonic;

    Double actTime = asteptime + timeOffset;
    UInt actStep = kstep + stepOffset;

    if (analysistype_ == STATIC ||
        analysistype_ == TRANSIENT) {
      solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);      

      if (saveSolHist_) {
        outFile_->WriteNodeHistoryTransient(*solTransient, actStep, actTime);
      }
        
      // Write derivatives
      if (analysistype_== TRANSIENT) {
        if (saveDeriv1Hist_ == TRUE) {
          solDeriv1_.SetAlgSysVector(getS1());
          outFile_->WriteNodeHistoryTransient(solDeriv1_, actStep, actTime);
        }

        if (saveDeriv2Hist_ == TRUE){
          solDeriv2_.SetAlgSysVector(getS2());
          outFile_->WriteNodeHistoryTransient(solDeriv2_, actStep, actTime);
        }
      }
        
    } 
    else if (analysistype_ == HARMONIC || analysistype_==MULTIHARMONIC) {
      //history nodes
      if (saveSolHist_) {
        solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);
        outFile_->WriteNodeHistoryHarmonic(*solHarmonic, actStep,
                                           actTime, complexFormat_);
      }
    }
  }



  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void PiezoPDE::ReadStoreResults() {

    ENTER_FCN( "PiezoPDE::ReadStoreResults" );


    // Construct vectors for restricted parameter search
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    StdVector<std::string> regionNames;
    std::string quantity;
  
    // *****************************
    // Determine nodal results
    // ***************************** 
    StdVector<std::string> nodeValues;
    keyVec  = pdename_, "storeResults", "nodeResults", "region";
    attrVec = "", "", "type";  
  
    // --- Mechanic Displacement & Electric Potential ---
    Enum2String(MECH_DISPLACEMENT, quantity);
    valVec = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, nodeValues);
    if (nodeValues.GetSize() > 0) {
      saveSol_ = TRUE;
      hasOutput_ = TRUE;
    }
    Enum2String(ELEC_POTENTIAL, quantity);
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
      solDeriv1_.SetNumDofs(dofspernode_);
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
      solDeriv2_.SetNumDofs(dofspernode_);
      solDeriv2_.SetPtrEQNData(eqnData_, ptgrid_);
      solDeriv2_.Init();
    }


    // Construct vectors for restricted search parameter

    //  keyVec  = "piezo", "storeResults", "elemResults", "region";
    attrVec = "", "", "type";

    // *****************************
    // Determine element results
    // *****************************
    StdVector<std::string> elemResults;
    keyVec  = pdename_, "storeResults", "elemResults", "region";

    // --- Electric Field Intensity ---
    Enum2String(ELEC_FIELD_INTENSITY, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, regionNames );
    ptgrid_->RegionNameToId( calcEfield_, regionNames );

    // If the the symbolic name is "all" 
    // compute electric field for all regions
    if ( calcEfield_.GetSize() == 1 && calcEfield_[0] == ALL_REGIONS ) {
      calcEfield_ = subdoms_;
    }
    
    if ( calcEfield_.GetSize() > 0 ) {
      hasOutput_ = TRUE;
      Info->PrintF( pdename_, "Computing electric field for regions:\n" );

      for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
        Info->PrintF( pdename_, "%s\n", regionNames[k].c_str() );
      }
      Info->PrintF( "", "\n" );

      if ( analysistype_ == HARMONIC || analysistype_==MULTIHARMONIC) {
        EfieldComplex_.SetNumSolutions(1);
        EfieldComplex_.SetSolutionType(ELEC_FIELD_INTENSITY);
        EfieldComplex_.SetNumElems(numElems_);
        EfieldComplex_.SetNumDofs(dim_);
        EfieldComplex_.SetPtrEQNData(eqnData_, ptgrid_);
        EfieldComplex_.Init(); 
      }
      else {
        Efield_.SetNumSolutions(1);
        Efield_.SetSolutionType(ELEC_FIELD_INTENSITY);
        Efield_.SetNumElems(numElems_);
        Efield_.SetNumDofs(dim_);
        Efield_.SetPtrEQNData(eqnData_, ptgrid_);
        Efield_.Init(); 
      }
    }

    // --- Mechanic Stress ---
    Enum2String(MECH_STRESS, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, regionNames );
    ptgrid_->RegionNameToId ( calcStress_, regionNames );
  
    // If the symbolic name is "all" compute stresses for all regions
    if ( calcStress_.GetSize() == 1 && calcStress_[0] == ALL_REGIONS ) {
      calcStress_ = subdoms_;
    }

    // Log to info file
    if ( calcStress_.GetSize() > 0 ) {
      hasOutput_ = TRUE;
      Info->PrintF( pdename_,
                    "Computing mechanical stress for regions:\n" );
      for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
        Info->PrintF( pdename_, "%s\n", regionNames[k].c_str() );
      }
      Info->PrintF( "", "\n" );

      if( analysistype_ == HARMONIC || analysistype_==MULTIHARMONIC) {

        // Resize solution arrays
        stressComplex_.SetNumSolutions(1);
        stressComplex_.SetSolutionType(MECH_STRESS);
        stressComplex_.SetNumElems(numElems_);

        // We always store for six components (unverg-file-format as capa does
        stressComplex_.SetNumDofs(6);
        stressComplex_.SetPtrEQNData(eqnData_, ptgrid_);
        stressComplex_.Init(0);
      }

      else {

        // Resize solution arrays
        stress_.SetNumSolutions(1);
        stress_.SetSolutionType(MECH_STRESS);
        stress_.SetNumElems(numElems_);

        // We always store for six components (unverg-file-format as capa does
        stress_.SetNumDofs(6);
        stress_.SetPtrEQNData(eqnData_, ptgrid_);
        stress_.Init(0);

      }
    }
  
    // --- Electric Charges ---
    // check for charge computation
    params->GetList( "region", regionNames, pdename_, "charge" );
    ptgrid_->RegionNameToId( chargeNeighborRegion_, regionNames );

    params->GetList( "element", regionNames, pdename_, "charge" );
    ptgrid_->RegionNameToId( calcCharge_, regionNames );


    if ( calcCharge_.GetSize() > 0 ) {

      hasOutput_ = TRUE;
      Info->PrintF( pdename_, "Computing electric charges for regions:\n" );
      for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
        Info->PrintF( pdename_, "%s\n", regionNames[k].c_str() );
      }
      Info->PrintF( "", "\n" );


      std::cout<<"values for setting into charges_:"<<std::endl;
      
      std::cout<<"eqnData_:"<<std::endl;
      std::cout<<eqnData_<<std::endl;
      
      std::cout<<"numElems_:"<<std::endl;
      std::cout<<numElems_<<std::endl;

      // Resize solution arrays
      if( analysistype_ == HARMONIC || analysistype_==MULTIHARMONIC) {
        chargesComplex_.SetNumSolutions(1);
        chargesComplex_.SetSolutionType(ELEC_CHARGE);
        chargesComplex_.SetNumElems(numElems_);
        chargesComplex_.SetNumDofs(1);
        chargesComplex_.SetPtrEQNData(eqnData_, ptgrid_);
        chargesComplex_.Init();
      }
      else {
        charges_.SetNumSolutions(1);
        charges_.SetSolutionType(ELEC_CHARGE);
        charges_.SetNumElems(numElems_);
        charges_.SetNumDofs(1);
        charges_.SetPtrEQNData(eqnData_, ptgrid_);
        charges_.Init();
      }
    } 


    // --- Surface Deformation ---
    // Check for computation of
    // volume of deformed surface
    params->GetList( "name", regionNames, pdename_, "volumeAboveDefSurf" );
    ptgrid_->RegionNameToId( volAboveDefSurfRegions_, regionNames );

    params->GetList( "dof", volAboveDefSurfDir_, pdename_,
                     "volumeAboveDefSurf" );

    // *****************************
    // Determine nodal history
    // *****************************
    StdVector<std::string> saveNodeHist; 
    keyVec  = pdename_, "storeResults", "nodeHistory", "saveNodes";
    attrVec = "", "", "type";
  
    // --- mechDisplacement / elecPotential ---
    Enum2String(MECH_DISPLACEMENT, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, saveNodeHist );
  
    if (saveNodeHist.GetSize() > 0) {
      saveSolHist_ = TRUE;
      hasOutput_ = TRUE;
      Info->PrintF( pdename_, " Saving mechDisplacement for Nodes:\n" );

      for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
        Info->PrintF( pdename_, "%s\n", saveNodeHist[k].c_str() );
      }
    }

    Enum2String(ELEC_POTENTIAL, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, saveNodeHist );
  
    if (saveNodeHist.GetSize() > 0) {
      saveSolHist_ = TRUE;
      hasOutput_ = TRUE;
      Info->PrintF( pdename_, " Saving elctric potential for Nodes:\n" );

      for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
        Info->PrintF( pdename_, "%s\n", saveNodeHist[k].c_str() );
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
    }

    // --- mechAcceleration ---
    Enum2String(MECH_ACCELERATION, quantity);
    valVec  = "", "", quantity;
    params->GetList( keyVec, attrVec, valVec, saveNodeHist );
  
    if (saveNodeHist.GetSize() > 0) {
      saveDeriv1Hist_ = TRUE;
      hasOutput_ = TRUE;
      Info->PrintF( pdename_, " Saving mechAcceleration for Nodes:\n" );
      for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
        Info->PrintF( pdename_, " %s\n", saveNodeHist[k].c_str() );
      }
    } 

    // *****************************
    // Determine element history
    // *****************************
    StdVector<std::string> saveElemHist;
    keyVec  = pdename_, "storeResults", "elemHistory", "saveElems";
    attrVec = "", "", "";
    valVec = "", "", "";
    params->GetList(keyVec, attrVec, valVec, saveElemHist);

    if (saveElemHist.GetSize() > 0) {
      std::string errMsg = pdename_;
      errMsg += ": Saving history elements is not implemented yet!\n";
      errMsg += "Meanwhile you can use 'unvtool' to extract element data.";
      Error( errMsg.c_str(), __FILE__, __LINE__);
    }

  }

  // ************************************************************
  //   PostProcess
  // ************************************************************

  void PiezoPDE::PostProcess() {

    ENTER_FCN( "PiezoPDE::PostProcess" );

    // calc electric field
    if (calcEfield_.GetSize() !=0 ) {
      if (analysistype_==HARMONIC || analysistype_==MULTIHARMONIC)
        CalcComplexValuedEfield();
      else
        CalcEfield();
    }
  
    // calc stresses
    if (calcStress_.GetSize() !=0 ) {
      if (analysistype_ == HARMONIC || analysistype_==MULTIHARMONIC)
        CalcComplexValuedStress();
      else
        CalcStress();
    }
    

    //calc charges
    if (calcCharge_.GetSize() !=0 ) {
   
      if (analysistype_ == HARMONIC || analysistype_==MULTIHARMONIC)
        CalcComplexValuedCharges();
      else
        CalcCharges();
    }

    //
    if (volAboveDefSurfRegions_.GetSize() ) {
      ComputeVolDefSurf(volAboveDefSurfRegions_, volAboveDefSurfDir_);
    }

  }


  void PiezoPDE::CalcEfield(){
    ENTER_FCN( "PiezoPDE::CalcEfield" );
    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);

    GradientFieldOp<Double> * FieldOp 
      = new GradientFieldOp<Double>(ptgrid_, this, eqnData_,
                                    *solhelp, ELEC_POTENTIAL, 
                                    isaxi_);


    // ------ Calculation of the electric field ------

    Vector<Double> LCoord;
      
    StdVector<Elem*> elemssd;
    UInt counterElems=0;
    Vector<Double> TempE;
    UInt pdeElem;
  
    // loop over all subdomains
    for (UInt isd=0; isd<calcEfield_.GetSize(); isd++)
      {
        // get vector of Elem of subdomain with color: subdoms[isd]
        ptgrid_->GetVolElems(elemssd,calcEfield_[isd]);
      
        // loop over elements of subdomain
        for (UInt iel=0; iel< elemssd.GetSize(); iel++,counterElems++)
          {
            elemssd[iel]->ptElem->GetCoordMidPoint(LCoord);
            FieldOp->CalcElemGradField( TempE, elemssd[iel], LCoord, 1); 
            pdeElem = eqnData_->Mesh2PDEElem(elemssd[iel]->elemNum);
            Efield_.SetElemResult(pdeElem-1,TempE);
          }
      }
    delete FieldOp;
  }


  void PiezoPDE::CalcComplexValuedEfield(){
    ENTER_FCN( "PiezoPDE::CalcComplexValuedEfield" );
    NodeStoreSol<Complex> * solhelp = 
      dynamic_cast<NodeStoreSol<Complex>*>(sol_);
    
    GradientFieldOp<Complex> * FieldOp 
      = new GradientFieldOp<Complex>(ptgrid_, this, eqnData_,
                                     *solhelp, ELEC_POTENTIAL, 
                                     isaxi_);

    // ------ Calculation of the electric field ------

    Vector<Double> LCoord;
 
      
    StdVector<Elem*> elemssd;
    UInt counterElems=0;
    Vector<Complex> TempE;
    UInt pdeElem;
  
    // loop over all subdomains
    for (UInt isd=0; isd<calcEfield_.GetSize(); isd++)
      {
        // get vector of Elem of subdomain with color: subdoms[isd]
        ptgrid_->GetVolElems(elemssd,calcEfield_[isd]);
      
        // loop over elements of subdomain
        for (UInt iel=0; iel< elemssd.GetSize(); iel++,counterElems++)
          {
            elemssd[iel]->ptElem->GetCoordMidPoint(LCoord);
            FieldOp->CalcElemGradField(TempE, elemssd[iel], LCoord, 1); 
            pdeElem = eqnData_->Mesh2PDEElem(elemssd[iel]->elemNum);
            EfieldComplex_.SetElemResult(pdeElem-1,TempE);
          }
      }
    delete FieldOp;
  }


  void PiezoPDE::CalcStress(){
    ENTER_FCN( "PiezoPDE::CalcSress" );
  
    //get the correct bilinearform
    ShortInt  stressDim, elecDim;
    Vector<Double> intPoint;
  
    if (subType_ == "planeStrain") 
      {
        
        stressDim = 3;
        elecDim   = 2;
      }
  
    else if (subType_ == "axi") {
      stressDim = 4;
      elecDim   = 2;
    }
  
    else if (subType_ == "3d") {
      stressDim = 6;
      elecDim   = 3;
    }
  
    else 
      Info->Error("StressOp: Unknown subtype in mech PDE! "
                  ,__FILE__,__LINE__);  
  
  
    Vector<Double> elemElecStress, elemStress, sortedStress;
    elemElecStress.Resize(stressDim+elecDim);
    elemStress.Resize(stressDim);
    elemElecStress.Init(0);
    elemStress.Init(0);
    sortedStress.Resize(6);

    // loop over all subdomains
    for (UInt isd=0; isd<subdoms_.GetSize(); isd++) {
    
      MaterialData actSDMat(materialData_[isd]);
      linPiezoInt * stress;

      if (subType_ == "planeStrain")
        stress = new piezoPlainStrainInt(actSDMat);

      else if (subType_ == "axi")      
        stress = new piezoAxiInt(actSDMat);
    
      else if (subType_ == "3d")
        stress = new linPiezo3DInt(actSDMat);
    
      // get vector of Elements of subdomains
      StdVector<Elem*> elemssd;     
      ptgrid_->GetVolElems(elemssd,subdoms_[isd]);
    
      // loop over elements of subdomain
      for (UInt iel=0; iel< elemssd.GetSize(); iel++) {
        UInt pdeElem = eqnData_->Mesh2PDEElem(elemssd[iel]->elemNum);
        elemssd[iel]->ptElem->GetCoordMidPoint(intPoint);

        //set element pointer
        BaseFE * ptEl = elemssd[iel]->ptElem;
        stress->SetElemPtr(ptEl);
      
        //set element solution    
        Matrix<Double> elSol;
        StdVector<UInt> connecth = elemssd[iel]->connect;
        sol_->GetElemSolutionAsMatrix(elSol, connecth);
        stress->SetActElemSol(elSol);
      
        //get coordinates of element
        Matrix<Double> ptCoord;
        GetElemCoords(connecth, ptCoord);
      
        Vector<Double> actStress; 
      
        //set the integration point
        stress->SetIntPoint(intPoint);
      
        //calculates the stress
        stress->CalcStressVec(elemElecStress,1,ptCoord);
        elemStress = elemElecStress.Part(0,stressDim-1);
        sortStresses(elemStress,sortedStress);
        stress_.SetElemResult(pdeElem-1, sortedStress);
      }
      // Delete integrator again (Stressabbau ;-)
      delete stress;

    }
  }
  
  void PiezoPDE::CalcComplexValuedStress(){
    ENTER_FCN( "PiezoPDE::CalcComplexValuedStress" );
  
    //get the correct bilinearform
    ShortInt stressDim, elecDim;
    Vector<Double> intPoint;
  
    if (subType_ == "planeStrain") 
      {      
        stressDim = 3;
        elecDim   = 2;
      }
  
    else if (subType_ == "axi") {
      stressDim = 4;
      elecDim   = 2;
    }
  
    else if (subType_ == "3d") {
      stressDim = 6;
      elecDim   = 3;
    }
  
    else 
      Info->Error("StressOp: Unknown subtype in mech PDE! "
                  ,__FILE__,__LINE__);  
  
  
    Vector<Complex> elemElecStress, elemStress, sortedStress;
    elemElecStress.Resize(stressDim+elecDim);
    elemStress.Resize(stressDim);
    elemElecStress.Init(0);
    elemStress.Init(0);
    sortedStress.Resize(6);

    // loop over all subdomains
    linPiezoInt *stress;
    for (UInt isd=0; isd<subdoms_.GetSize(); isd++) {
    
      MaterialData actSDMat(materialData_[isd]);

      if (subType_ == "planeStrain")
        stress = new piezoPlainStrainInt(actSDMat);

      else if (subType_ == "axi")      
        stress = new piezoAxiInt(actSDMat);
    
      else if (subType_ == "3d")
        stress = new linPiezo3DInt(actSDMat);
    
      // get vector of Elements of subdomains
      StdVector<Elem*> elemssd;     
      ptgrid_->GetVolElems(elemssd,subdoms_[isd]);
    
      // loop over elements of subdomain
      for (UInt iel=0; iel< elemssd.GetSize(); iel++) {
        UInt pdeElem = eqnData_->Mesh2PDEElem(elemssd[iel]->elemNum);
        elemssd[iel]->ptElem->GetCoordMidPoint(intPoint);
            
        //set element pointer
        BaseFE * ptEl = elemssd[iel]->ptElem;
        stress->SetElemPtr(ptEl);
      
        //set element solution    
        Matrix<Complex> elSol;
        StdVector<UInt> connecth = elemssd[iel]->connect;
        sol_->GetElemSolutionAsMatrix(elSol, connecth);
        stress->SetActElemSol(elSol);
      
        //get coordinates of element
        Matrix<Double> ptCoord;
        GetElemCoords(connecth, ptCoord);
      
        Vector<Double> actStress; 
      
        //set the integration point
        stress->SetIntPoint(intPoint);
      
        //calculates the stress
  
        stress->CalcStressVec(elemElecStress,1,ptCoord);

        elemStress = elemElecStress.Part(0,stressDim-1);
        sortStresses(elemStress,sortedStress);
        stressComplex_.SetElemResult(pdeElem-1, sortedStress);

      }

      // Delete integrator again (Stressabbau ;-)
      delete stress;

    }
  }
  
  
  void PiezoPDE::CalcCharges(){
    ENTER_FCN( "PiezoPDE::CalcCharges" );
    
    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
    StdVector<SurfElem*> surfElems;
    Elem * ptVolElem;
    SurfElem * ptSurfElem;
    Vector<Double> lCoordSurf, lCoordVol, elemDField, Efield, normal;
    ElecChargeOp * chargeOp;
    BaseFE * ptSurfElemFE, * ptVolElemFE;
    Double elemNormalD = 0.0;
    Double charge = 0.0;
    Double Ecomp, Dcomp;
    UInt pdeElemNum = 0;
    Matrix<Double> elSol;
    Matrix<Double> ptCoord;
    Double normSign = 0.0;
    Integer regionIndex = 0;

    ShortInt stressDim, elecDim;
    Vector<Double> intPoint;
  

    GradientFieldOp<Double> * FieldOp;
    if (isHysteresis_) {
      FieldOp = new GradientFieldOp<Double>(ptgrid_, this, eqnData_,
                                            *solhelp, ELEC_POTENTIAL, 
                                            isaxi_);
    }

    if (subType_ == "planeStrain") {
      stressDim = 3;
      elecDim   = 2;
    } 
    else if (subType_ == "axi") {
      stressDim = 4;
      elecDim   = 2;
    } else if (subType_ == "3d") {
      stressDim = 6;
      elecDim   = 3;
    }
    else 
      Info->Error("StressOp: Unknown subtype in mech PDE! "
                  ,__FILE__,__LINE__);  
    
  
    Vector<Double> elemElecStress;
    elemElecStress.Resize(stressDim+elecDim);
    elemElecStress.Init(0);

    
    // Create vector with interpolation coordinate.
    // For simplicity we only evaluate the integral
    // in coordinate origin
    lCoordSurf.Resize(dim_-1);
    lCoordSurf.Init(0);

    //charge operator  
    chargeOp = new ElecChargeOp(ptgrid_, this, eqnData_, isaxi_);
                              
    Vector<Double> chargeSD(calcCharge_.GetSize());
    chargeSD.Init(0);
    Vector<Double> averageE(calcCharge_.GetSize());
    averageE.Init(0);
          
    // loop over all subdomains
    for (UInt iSD=0; iSD<calcCharge_.GetSize(); iSD++){
    
      // get surface and according volume elements
      ptgrid_->GetSurfElems( surfElems, calcCharge_[iSD] );

      // loop over all surface elements
      for (UInt iElem=0; iElem<surfElems.GetSize(); iElem++)
        {

          // Determine, which volume element is the right neighbour for the 
          // calculation

          std::cout<<"->chargeNeighborRegion_:"<<std::endl;
          std::cout<<chargeNeighborRegion_<<std::endl;

          if ( chargeNeighborRegion_.
               Find(surfElems[iElem]->ptVolElem1->regionId) != -1 ) {
            ptVolElem = surfElems[iElem]->ptVolElem1;
            normSign = -1.0;
          } else {
            ptVolElem = surfElems[iElem]->ptVolElem2;
            normSign = 1.0;
          }

          normSign *= (double) surfElems[iElem]->normalSign;

          ptSurfElemFE = surfElems[iElem]->ptElem;
          ptVolElemFE = ptVolElem->ptElem;
          const StdVector<UInt> & surfConnect = surfElems[iElem]->connect;
          const StdVector<UInt> & volConnect = ptVolElem->connect;

          // calculate volume integration coordinates from
          // surfe integration coordinates for evaluating the 
          // electric flux density on the surface of the volume
          // element
          ptSurfElemFE->GetCoordMidPoint(lCoordSurf);
          ptVolElemFE->GetLocalIntPoints4Surface(surfConnect, volConnect,
                                                 lCoordSurf, lCoordVol);          

          // Find correct material for volume element
          regionIndex = subdoms_.Find( ptVolElem->regionId );
          if ( regionIndex == -1 ) {
            (*error) << "PiezoPDE:CalcCharges:: The region with Name " 
                     << ptgrid_->RegionIdToName(ptVolElem->regionId)
                     << " of surface element Nr. " << ptSurfElem->elemNum
                     << "is not contained in my set of regions!.";
            Error( __FILE__, __LINE__ );
          }
          
          MaterialData actSDMat(materialData_[regionIndex]);
          linPiezoInt * stress;
           
          if (subType_ == "planeStrain")
            stress = new piezoPlainStrainInt(actSDMat);
           
          else if (subType_ == "axi")      
            stress = new piezoAxiInt(actSDMat);
           
          else if (subType_ == "3d")
            stress = new linPiezo3DInt(actSDMat);
      
          //set volume element
          stress->SetElemPtr(ptVolElemFE);
      
          //set element solution  
          sol_->GetElemSolutionAsMatrix(elSol, volConnect);
          stress->SetActElemSol(elSol);
    
          //set the integration point
          stress->SetIntPoint(lCoordVol);

          //get coordinates of element
          GetElemCoords(volConnect, ptCoord);

          // check for electric hysteresis
          if (isHysteresis_) {
            FieldOp->CalcElemGradField( Efield, ptVolElem, lCoordVol, 1);
            UInt comp =  materialData_[regionIndex].GetDirPol() - 1;
            Ecomp = Efield[comp];
            UInt elemNr = ptVolElem->elemNum;
	    
            Dcomp = solveStep_->GetHystOperator(regionIndex)->computeValue(Ecomp,elemNr); 
            //	    std::cerr << "nr=" << elemNr << " D=" << Dcomp << std::endl;
            //	    std::cerr << Ecomp << " " << Dcomp << std::endl;
            averageE[iSD] += Ecomp/surfElems.GetSize();
          }

          //calculates the stress
          Vector<Double> actElecD;

          if (isHysteresis_) {
            stress->CalcStressVec(elemElecStress,1,ptCoord, 0, Dcomp);
          }
          else {
            stress->CalcStressVec(elemElecStress,1,ptCoord);
          }

          actElecD = elemElecStress.Part(stressDim,elemElecStress.GetSize()-1);

          std::cout<<"actElecD:"<<std::endl;
          std::cout<<actElecD<<std::endl;

          //	  actElecD[1] = Dcomp;

          // Calc surface normal, which points from the surface element
          // into the volume
          ptgrid_->CalcSurfNormal(normal, *surfElems[iElem]);
          normal *= normSign;
          elemNormalD = normal * actElecD;

          std::cout<<"elemNormalD"<<std::endl;
          std::cout<<elemNormalD<<std::endl;

          std::cout<<"lCoordSurf"<<std::endl;
          std::cout<<lCoordSurf<<std::endl;


          chargeOp->CalcElemCharge(charge, surfElems[iElem], 
                                   lCoordSurf, elemNormalD);

          pdeElemNum = eqnData_->Mesh2PDEElem(ptVolElem->elemNum);
        
          // Create temporary vector, since SetElemResult only
          // can handle these
          Vector<Double> chargeVec(1);
          chargeVec[0] = charge;
          charges_.SetElemResult(pdeElemNum-1, chargeVec);
          chargeSD[iSD] += charge;   
          std::cout<<"charge:"<<std::endl;
          std::cout<<charge<<std::endl;     

          std::cout<<"Final charge = "<<std::endl;
          std::cout<<chargeSD[iSD]<<std::endl;
          getchar();

          // Delete integrator again (reeller Stressabbau !)
          delete stress;
        }
    }

    //print to info-file
    Info->PrintF(pdename_, "Computed Charges: ");
    Info->PrintVec(chargeSD);

    if (isHysteresis_) {
      Info->PrintF(pdename_, "Averaged Efield: ");
      Info->PrintVec(averageE);
    }

    //  Warning ("Charges are written to unv/gmv file, although capapost \
      // can not draw them", __FILE__, __LINE__);
  }



  void PiezoPDE::CalcComplexValuedCharges(){
    ENTER_FCN( "PiezoPDE::CalcComplexValuedCharges" );
    StdVector<SurfElem*> surfElems;
    Elem* ptVolElem;
    SurfElem * ptSurfElem;
    Matrix<Double> ptCoord;
    Matrix<Complex> elSol;
    Vector<Double> lCoordSurf, lCoordVol, normal;
    Vector<Complex> elemDField;
    ElecChargeOp * chargeOp;
    BaseFE * ptSurfElemFE, * ptVolElemFE;
    Complex   elemNormalD = Complex(0.0,0.0);
    Complex charge = Complex(0.0,0.0);
    UInt pdeElemNum = 0;
    Integer regionIndex = 0;
    Double normSign = 0.0;
    Vector<Complex> chargeSD;

    ShortInt stressDim, elecDim;
    Vector<Double> intPoint;
  
    if (subType_ == "planeStrain") 
      {
        stressDim = 3;
        elecDim   = 2;
      }
  
    else if (subType_ == "axi") {
      stressDim = 4;
      elecDim   = 2;
    }
  
    else if (subType_ == "3d") {
      stressDim = 6;
      elecDim   = 3;
    }
  
    else 
      Info->Error("StressOp: Unknown subtype in mech PDE! ",__FILE__,__LINE__);  
  
    Vector<Complex> elemElecStress;
    elemElecStress.Resize(stressDim+elecDim);
    elemElecStress.Init(0);
    
    // Create vector with interpolation coordinate.
    // For simplicity we only evaluate the integral
    // in coordinate origin

    lCoordSurf.Resize(dim_-1);
    lCoordSurf.Init(0);

    //charge operator  
    chargeOp = new ElecChargeOp(ptgrid_, this, eqnData_, isaxi_);

    chargeSD.Resize(calcCharge_.GetSize());
    chargeSD.Init(0);
                              
    // loop over all subdomains
    for (UInt iSD=0; iSD<calcCharge_.GetSize(); iSD++){    
    
      // get surface and acoording volume elements
      ptgrid_->GetSurfElems( surfElems, calcCharge_[iSD] );
      complexValuedCharge_.Resize(surfElems.GetSize());
    

      // loop over all surface elements
      for (UInt iElem=0; iElem<surfElems.GetSize(); iElem++)
        {
        
          // Determine, which volume element is the right neighbour for the 
          // calculation
          if ( chargeNeighborRegion_.
               Find(surfElems[iElem]->ptVolElem1->regionId) != -1 ) {
            ptVolElem = surfElems[iElem]->ptVolElem1;
            normSign = -1.0;
          } else {
            ptVolElem = surfElems[iElem]->ptVolElem2;
            normSign = 1.0;
          }
        
          normSign *= (double) surfElems[iElem]->normalSign;
      
          ptSurfElemFE = surfElems[iElem]->ptElem;
          ptVolElemFE = ptVolElem->ptElem;
          const StdVector<UInt> & surfConnect = surfElems[iElem]->connect;
          const StdVector<UInt> & volConnect = ptVolElem->connect;

          // calculate volume integration coordinates from
          // surfe integration coordinates for evaluating the 
          // electric flux density on the surface of the volume
          // element
          ptSurfElemFE->GetCoordMidPoint(lCoordSurf);
          ptVolElemFE->GetLocalIntPoints4Surface(surfConnect, volConnect,
                                                 lCoordSurf, lCoordVol);

          // Find correct material for volume element
          regionIndex = subdoms_.Find( ptVolElem->regionId );
          if ( regionIndex == -1 ) {
            (*error) << "PiezoPDE:CalcCharges:: The region with Name " 
                     << ptgrid_->RegionIdToName(ptVolElem->regionId)
                     << " of surface element Nr. " << ptSurfElem->elemNum
                     << "is not contained in my set of regions!.";
            Error( __FILE__, __LINE__ );
          }

          MaterialData actSDMat(materialData_[regionIndex]);
          linPiezoInt * stress;
        
          if (subType_ == "planeStrain")
            stress = new piezoPlainStrainInt(actSDMat);
        
          else if (subType_ == "axi")      
            stress = new piezoAxiInt(actSDMat);
        
          else if (subType_ == "3d")
            stress = new linPiezo3DInt(actSDMat);

          //set volume element
          stress->SetElemPtr(ptVolElemFE);

          //set element solution  
          sol_->GetElemSolutionAsMatrix(elSol, volConnect);

          //           std::cout<<"elSol in PiezoPDE:"<<std::endl;
          //           std::cout<<elSol<<std::endl;
          //           getchar();
          stress->SetActElemSol(elSol);

          //set the integration point
          stress->SetIntPoint(lCoordVol);

          //get coordinates of element
          GetElemCoords(volConnect, ptCoord);

          //calculates the stress
          Vector<Complex> actElecD;       
          stress->CalcStressVec(elemElecStress,1,ptCoord);

          actElecD = elemElecStress.Part(stressDim,elemElecStress.GetSize()-1);

          elemNormalD = Complex();

          // Calc surface normal, which points from the surface element
          // into the volume
          ptgrid_->CalcSurfNormal(normal, *surfElems[iElem]);
          normal *= normSign;
        
          // scalarProduct lCoordCal*actElecD
          for (UInt i=0; i<normal.GetSize(); i++)
            elemNormalD += normal[i] * actElecD[i];
        
          //lCoordVol.Mult(actElecD, elemNormalD);
          chargeOp->CalcElemCharge(charge, surfElems[iElem], 
                                   lCoordSurf, elemNormalD);
        
          pdeElemNum = eqnData_->Mesh2PDEElem(ptVolElem->elemNum);

        
          // Create temporar vector, since SetElemResult only
          // can handle these
          Vector<Complex> chargeVec(1);
          chargeVec[0] = charge;
          complexValuedCharge_[iElem]=charge;
          chargesComplex_.SetElemResult(pdeElemNum-1, chargeVec); 
          chargeSD[iSD] += charge;        

          // Delete integrator again (Stressabbau ;-)
          std::cout<<"Final charge " <<std::endl;
          std::cout<<chargeSD[iSD]<<std::endl;
          delete stress;
        }

    }
    //print to info-file
    Info->PrintF(pdename_, "Computed surface charge: ");
    Info->PrintVec(chargeSD);


  }

} // end of namespace
