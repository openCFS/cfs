#include <fstream>
#include <iostream>
#include <string>
#include <math.h>

#include "acousticPDE.hh" 
#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "Forms/forms_header.hh"
#include "Forms/abcInt.hh"
#include "Forms/linNeumannInt.hh"
#include "Estimator/spaceerror.hh"
#include "newmark.hh"
#include "newmarkFracDamp.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "PDE/scalarnodeEQN.hh"
#include "Utils/mathfunctions.hh"
#include "Utils/nodestoresol.hh"
#include "Driver/solveStepAcoustic.hh"
#include "Driver/solveStepAcousticBubble.hh"
#include "Driver/solveStepAcousticMechBubble.hh"
#include "CoupledPDE/pdecoupling.hh"

#ifdef TCL_INTERFACE
#include "DataInOut/Scripting/cfsmessenger.hh" 
#endif

#ifdef PARALLEL
#include <mpi.h>
#endif

namespace CoupledField {



  // =========================================================================
  // set solution information
  // =========================================================================
  AcousticPDE::AcousticPDE(Grid * aptgrid, TimeFunc *aptTimeFunc, 
                           WriteResults *aptOut)
    :SinglePDE(aptgrid,aptOut,aptTimeFunc) {

    ENTER_FCN( "AcousticPDE::AcousticPDE" );

    dofspernode_ = 1;
    solDofs_ = 1;
    pdename_          = "acoustic";
    pdematerialclass_ = FLUID;
    fracDamping_ = FALSE;

    isMechCoupled_ = FALSE;
    isNrbcCoupled_ = FALSE;
    saveRHSval_ = FALSE;
    saveRHSvalHist_ = FALSE;

    coarsealpha_ = 0.01; // solver parameter, see basePDE.cc

    // PDE formulation either in acoustic potential or pressure
    std::string str;
    params->Get( "formulation", str, "pdeList", pdename_ );
    String2Enum( str, formulation_ );
    str = "Using * " + str + " as state variable in formulation of PDE\n";
    Info->PrintF( pdename_, str.c_str() );

    // class NodeStoreSol will be initialized with acoustic potential
    if ( formulation_ ==  ACOU_PRESSURE) {
      solTypes_ = ACOU_PRESSURE;
    }
    else {
      solTypes_ = ACOU_POTENTIAL;
    }

    // timestepping formulation
    params->Get( "timeSteppingFormulation", str, "pdeList", pdename_ );
    if ( str == "effMassMatrix" ) {
      effectiveMass_ = TRUE;
      Info->PrintF( pdename_, 
                    "      * effective mass matrix timestepping\n");
    } 
    else {
      effectiveMass_ = FALSE;
      Info->PrintF( pdename_, 
                    "      * effective stiffness matrix timestepping\n");
    }

    // axisymmetric setup    
    isaxi_ = params->HasValue( "type", "axi", "geometry" );

    //To check if acoustic is coupled with nrbc and set isNrbcCoupled
    if( params->IsSet( "isCoupledNrbc","pdeList" ,"acoustic" ) ) {
      isNrbcCoupled_=TRUE;
    }

    // ===========================
    // DODO GridAdaption
    // check if <specialNodes> is given, to determine whether we to write them
    // to file,
    StdVector<std::string> vecKey, vecVals, vecAttr;
    StdVector<std::string> vecNodes, vecFiles;
    vecKey = pdename_, "storeResults", "specialNodes", "saveNodes";
    vecAttr = "", "", "";
    vecVals = "", "", "";
    params->GetList(vecKey, vecAttr, vecVals, vecNodes);
    
    // specialNodes is given => set flag
    if(vecNodes.GetSize() > 0) {
      m_bWriteSpecialBCs = TRUE;
      
      // get file to store specialNodes
      vecKey = pdename_, "storeResults", "specialNodes", "file";
      std::string strFile;
      params->Get(vecKey, strFile);
      
      // setup grid adaption object
      // it needs: the .data-file, the .xml-file and the grid
      m_pGridAdaption = new GridAdaption(strFile, aptgrid);
    }
    else {
      m_bWriteSpecialBCs = FALSE;
    }
    // DODO
    // ===========================

  }


  // *********************************************
  //   Check what type of damping should be used
  // *********************************************
  void AcousticPDE::ReadDampingInformation( ) {

    ENTER_FCN( "AcousticPDE::ReadDampingInformation" );

    fracMemory_ = 0;
    Integer firstFrac=-1;

    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    for (UInt k = 0; k < subdoms_.GetSize(); k++) {

      // Be aware, that the RegionId entries in subdoms_ need not be in order!
      // We follow the order in subdoms_, which has to be in acordance with
      //  AcousticPDE::DefineIntegrators()

      RegionIdType actRegion = subdoms_[k];
      std::string actRegionName;
      actRegionName = ptgrid_->RegionIdToName( actRegion );
      keyVec = "acoustic" , "region" , "damping" , "type";
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
                      actRegionName.c_str() );;
      }
      else if (dampInfo[0] == "rayleigh") {
        dampingList_[actRegion] = RAYLEIGH;
        Info->PrintF( pdename_, 
                      "      * RAYLEIGH damping for region: %s\n",
                      actRegionName.c_str() );
      }
      else if (dampInfo[0] == "thermoViscous") {
        dampingList_[actRegion] = THERMOVISCOUS;
        Info->PrintF( pdename_, 
                      "      * THERMOVISCOUS damping for region: %s\n",
                      actRegionName.c_str() );
      }
      else if (dampInfo[0] == "fractional") {
        fracDamping_ = TRUE;
        Info->PrintF( pdename_, 
                      "      * FRACTIONAL damping for region: %s\n",
                      actRegionName.c_str() );

        // Find first region containing fractional damping
        if ( firstFrac < 0 )
          firstFrac = k;

        // Gather additional information for fractional damping model
        keyVec = "acoustic" , "region" , "damping" , "fracAlg";
        StdVector<std::string> fracAlg;
        params->GetList( keyVec, attrVec, valVec, fracAlg );

        keyVec = "acoustic" , "region" , "damping" , "fracMemory";
        StdVector<UInt> fracMem;
        params->GetList( keyVec, attrVec, valVec, fracMem );

        keyVec = "acoustic" , "region" , "damping" , "interpolation";
        StdVector<std::string> interpol;
        params->GetList(  keyVec, attrVec, valVec, interpol );


        // Include fracAlg and interpolation info in dampingList
        if( fracAlg.IsEmpty() || fracMem.IsEmpty() || interpol.IsEmpty() ) {
          (*error) << "Specify attributes fracAlg, fracMemory " 
                   << "and interpolation for fractional damping model!";
          Error( __FILE__, __LINE__ ); 
        }
        else if  ( fracAlg[0] == "gl" ) {
          Info->PrintF( "", "\t\t\t using Gruenwald-Letnikov algorithm,\n");
          if (interpol[0] == "no" )
            dampingList_[actRegion] = FRACTIONAL_GL;
          else {
            dampingList_[actRegion] = FRACTIONAL_GL_INT;
            Info->PrintF("", 
                         "\t\t\t linear interpol. of single past values\n\n");
          }
        }
        else if ( fracAlg[0] == "blank" ) {
          Info->PrintF( "", "\t\t\t using Blanks algorithm,\n");
          if (interpol[0] == "no" )
            dampingList_[actRegion] = FRACTIONAL_BLANK;
          else {
            dampingList_[actRegion] = FRACTIONAL_BLANK_INT;
            Info->PrintF("", 
                         "\t\t\t linear interpol. of single past values\n\n");
          }
        }

        // up to now take maximum of specified fracMemory values
        if ( fracMem[0] > fracMemory_ )
          fracMemory_ = fracMem[0];
      }
    }

    if ( fracDamping_== TRUE ) {
      Info->PrintF(pdename_, "Memory size for fractional damping  is: %d\n",
                   fracMemory_ );
    }
  }

  void AcousticPDE::ReadSpecialBCs() {
    ENTER_FCN( "AcousticPDE::ReadSpecialBCs" );

    // ***************************************************************
    //   If no other damping type is specified and we have absorbing
    //   boundary conditions, then use ABCDAMP
    // ***************************************************************
    StdVector<std::string> auxVec;
    absorbingBCs_ = FALSE;
    params->GetList( "name", auxVec, pdename_, "absorbingBCs" );
    ptgrid_->RegionNameToId( absBCs_, auxVec );

    if ( absBCs_.GetSize() ) {
      absorbingBCs_ = TRUE;
      Info->PrintF( pdename_, "Apply Absorbing Boundary Conditions\n" );
    }
 
    // **************************************************
    //   Check what type of bubble model should be used
    // **************************************************

    // vectors for parameter handling
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    
    keyVec = pdename_, "bubbles", "bubbleType";
    attrVec = "", "tag";
    valVec =  "", bcSequenceTag_;
    params->GetList(keyVec, attrVec, valVec, auxVec);

    if ( auxVec.GetSize() == 1 ) {
      String2Enum( auxVec[0], bubbleDynType_ );
      //Set bubbledensity
      keyVec = pdename_, "bubbles", "bubbleNumDensity";
      params->Get(keyVec, attrVec, valVec, bubbleDensity_);
    }
    else if ( auxVec.GetSize() > 1 ) {
      Error("Specification of bubble type not unique in xml-file", __FILE__,
            __LINE__ );
    }
    else {
      bubbleDynType_ = NOBUBBLETYPE;
    }
  }


  // *************************************************************
  //   Check what type of nonlinear PDE formulation should be used
  // *************************************************************
  void  AcousticPDE::InitNonLin() {
    ENTER_FCN( "AcousticPDE::InitNonLin" );

    nonLin_ = FALSE; //declaration in StdPDE.hh
    StdVector<std::string> strVec;
    params->GetList( "nonLinear", strVec, pdename_, "region" );
    nonLinPDEName_.Resize(strVec.GetSize());

    for ( UInt k = 0; k < strVec.GetSize(); k++ ) {

      if ( strVec[k] == "no" )
        ;
      else if ( strVec[k] == "westervelt" ) {
        nonLin_ = TRUE;
        nonLinPDEName_[k] = WESTERVELT;
        Info->PrintF(pdename_, "      * Westervelt equation for region: %d\n",
                     k );
        if ( solTypes_ == ACOU_POTENTIAL )
          Error("Acoustic potential formulation not supported for \
Westervelt equation!" ,__FILE__,__LINE__);
      }
      else if ( strVec[k] == "kuznetsov" ) {
        nonLin_ = TRUE;
        nonLinPDEName_[k] = KUZNETSOV;
        Info->PrintF(pdename_, "      * Kuznetsov equation for region: %d\n",
                     k );
        if ( solTypes_ == ACOU_PRESSURE )
          Error("Acoustic pressure formulation not supported for \
Kuznetsov equation!" ,__FILE__,__LINE__);
      }
    }

    if( nonLin_ || bubbleDynType_ != NOBUBBLETYPE ) {
      // solution method
      params->Get("method", nonLinMethod_, pdename_, "nonLinear" );
      // perform logging?
      nonLinLogging_ = params->IsSet( "logging", pdename_, "nonLinear" );
      // type of line search
      params->Get("type", lineSearch_, pdename_, "lineSearch" );
      // incremental stopping criterion
      params->Get("incStopCrit", incStopCrit_, pdename_, "nonLinear" );
      // residual stopping criterion
      params->Get("resStopCrit", residualStopCrit_, pdename_, "nonLinear" );
      // maximal number of NL-iterations
      params->Get("maxNumIters", nonLinMaxIter_, pdename_, "nonLinear");
    }
  }


  void AcousticPDE::DefineIntegrators() {

    ENTER_FCN( "AcousticPDE::DefineIntegerators" );

    // help variables for parameter checking
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    Double density, compressibility, c0;
    Double coeffmass, coeffdamp;

    for (UInt actSD = 0; actSD < subdoms_.GetSize(); actSD++) {

      // ********************************************************************
      //   Attention:
      //   In AcousticPDE ALL integrators are multiplied with density!
      // ********************************************************************

      materials_[subdoms_[actSD]]->GetScalar(density,DENSITY,REAL);
      materials_[subdoms_[actSD]]->GetScalar(compressibility,ACOU_BULK_MODULUS,REAL);
      c0 = sqrt(compressibility/density);

      // if pde couples with mechanic, we have to multiply the density by -1
      if ( isMechCoupled_ == TRUE && formulation_ !=  ACOU_PRESSURE) {
        density *= -1.0;
      }

      // ********************************************************************
      //   Linear wave equation
      // ********************************************************************


      //do the stiffness part;
      IntegratorDescriptor * stiffIntDescr;

      //check  for PML!
      RegionIdType actRegion = subdoms_[actSD];
      std::string actRegionName;
      actRegionName = ptgrid_->RegionIdToName( actRegion );
      keyVec = "acoustic" , "region" , "pml" , "type";
      attrVec= ""         , "name"   , "";
      valVec = ""         , actRegionName, "";
      StdVector<std::string> pmlInfo;
      params->GetList( keyVec, attrVec, valVec, pmlInfo);

      if (pmlInfo.GetSize() > 0) {
        //read data for PML layer

        //type of PML damping
        std::string dampingTypePML;

        // inner / outer region
        Matrix<Double> inner;
        Matrix<Double> outer;

        //damping factor
        Double dampPML;

        ReadDataPML(dampingTypePML, inner, dampPML, actRegion);
        dampPML *= c0;

        GetPMLLayerData(inner, outer, actSD);

        //====================================================================
        //	 stiffness integrator for PML
        //====================================================================

        std::string formsType = "laplaceInt";

        //set real part
        BaseForm * bilinearStiffReal = 
          new PMLInt(formsType, density, dampingTypePML, dampPML, isaxi_);

        bilinearStiffReal->SetPosPML(inner,outer);
        DataType matType = REAL;
        bilinearStiffReal->SetMatDataType(matType);

        IntegratorDescriptor * stiffIntDescrReal = 
          new IntegratorDescriptor(bilinearStiffReal, STIFFNESS);

        stiffIntDescrReal->SetPDEIds(this, this);
        stiffIntDescrReal->SetMatDataType(matType);
        assemble_->AddIntegrator(stiffIntDescrReal, subdoms_[actSD]);

        //set imaginary part
        BaseForm * bilinearStiffImag = 
          new PMLInt(formsType, density, dampingTypePML, dampPML, isaxi_);

        bilinearStiffImag->SetPosPML(inner,outer);
        matType = IMAG;
        bilinearStiffImag->SetMatDataType(matType);

        IntegratorDescriptor * stiffIntDescrImag = 
          new IntegratorDescriptor(bilinearStiffImag, STIFFNESS);


        stiffIntDescrImag->SetPDEIds(this, this);
        stiffIntDescrImag->SetMatDataType(matType);
        assemble_->AddIntegrator(stiffIntDescrImag, subdoms_[actSD]);


        //====================================================================
        //	 mass integrator for PML
        //====================================================================
        formsType = "massInt";

        //	Double dampMass = (2.0/3.0)*density/c0;
        Double massFactor = density/(c0*c0);

        //set real part
        BaseForm * bilinearMassReal = 
          new PMLInt(formsType, massFactor, dampingTypePML, dampPML, isaxi_);

        bilinearMassReal->SetPosPML(inner,outer);
        matType = REAL;
        bilinearMassReal->SetMatDataType(matType);

        IntegratorDescriptor * massIntDescrReal = 
          new IntegratorDescriptor(bilinearMassReal, MASS);

        massIntDescrReal->SetPDEIds(this, this);
        massIntDescrReal->SetMatDataType(matType);
        assemble_->AddIntegrator(massIntDescrReal, subdoms_[actSD]);

        //set imaginary part
        BaseForm * bilinearMassImag = 
          new PMLInt(formsType, massFactor, dampingTypePML, dampPML, isaxi_);

        bilinearMassImag->SetPosPML(inner,outer);
        matType = IMAG;
        bilinearMassImag->SetMatDataType(matType);

        IntegratorDescriptor * massIntDescrImag = 
          new IntegratorDescriptor(bilinearMassImag, MASS);

        massIntDescrImag->SetPDEIds(this, this);
        massIntDescrImag->SetMatDataType(matType);
        assemble_->AddIntegrator(massIntDescrImag, subdoms_[actSD]);
      }

      else {
        // stiffness integrator 
        BaseForm * bilinearStiff = new LaplaceInt(density,isaxi_);        
        stiffIntDescr = new IntegratorDescriptor(bilinearStiff, STIFFNESS);

        stiffIntDescr->SetPDEIds(this, this);

        // mass integrator
        coeffmass = density / (c0*c0);

        BaseForm * bilinearMass  = new MassInt(coeffmass, dofspernode_, isaxi_);
        IntegratorDescriptor * massIntDescr = 
          new IntegratorDescriptor(bilinearMass, MASS);

        massIntDescr->SetPDEIds(this, this);
      

        // ********************************************************************
        //   Additional terms for damping
        // ********************************************************************

        if ( dampingList_.size() > 0 ) {
	  
          // We check, if damping has been specified for all regions.
          if ( dampingList_.size() != subdoms_.GetSize() ) {
            (*warning) << "Mismatch between dampingList_ and subdoms_!"
                       << "Size(dampingList_): " << dampingList_.size()
                       << "Size(subdoms_): " << subdoms_.GetSize();
            Warning(__FILE__, __LINE__);  
          }
	  
          if (dampingList_[subdoms_[actSD]] == RAYLEIGH) {
            // This works even after assemble_->AddIntegrator() is executed
            //   because of the pointers...
	    Double alpha, beta;
	    materials_[subdoms_[actSD]]->GetScalar(alpha,RAYLEIGH_ALPHA,REAL);
	    materials_[subdoms_[actSD]]->GetScalar(beta,RAYLEIGH_BETA,REAL);
    
            // stiffness part
            stiffIntDescr->SetSecondaryMat(DAMPING, beta, analysistype_);
	    
            // mass part
            massIntDescr->SetSecondaryMat(DAMPING, alpha, analysistype_);
          }
	  
	  
          else if ( dampingList_[subdoms_[actSD]] == THERMOVISCOUS ) {
	    Double viscousAlpha;
	    materials_[subdoms_[actSD]]->GetScalar(viscousAlpha, ACOU_ALPHA,REAL);

            coeffdamp  =  density * 2.0 * viscousAlpha * c0;
            BaseForm * bilinearStiff  = new LaplaceInt(coeffdamp, isaxi_);  
            IntegratorDescriptor * dampIntDescr = 
              new IntegratorDescriptor(bilinearStiff, DAMPING);
	    
            dampIntDescr->SetPDEIds(this, this);   
            assemble_->AddIntegrator(dampIntDescr, subdoms_[actSD]);
          }
	  
          else if ( dampingList_[subdoms_[actSD]] == FRACTIONAL_GL ||
                    dampingList_[subdoms_[actSD]] == FRACTIONAL_GL_INT ) {

	    Double fracAlpha, fracExp;
	    materials_[subdoms_[actSD]]->GetScalar(fracAlpha,ACOU_ALPHA,REAL);
	    materials_[subdoms_[actSD]]->GetScalar(fracExp,FRACTIONAL_EXPONENT,REAL);

            coeffdamp = - density * 2.0 * fracAlpha / c0 / sin((fracExp-1.0)*PI/2.0);

            BaseForm * bilinearDamp  = 
              new MassInt(coeffdamp, dofspernode_, isaxi_);
            bilinearDamp->SetFracDamping();

            // formulation using DAMPING matrix
            // adapt NewmarkFracDamp::Init and StdPDE::GetFracDampMatrixCoeff
            // IntegratorDescriptor * dampIntDescr = 
            //   new IntegratorDescriptor(bilinearDamp, DAMPING);
	    
            // two matrices formulation
            // added to STIFFNESS matrix because, because 
            //   matrix_factors[STIFFNESS] = 1.0
            IntegratorDescriptor * dampIntDescr = 
              new IntegratorDescriptor(bilinearDamp, STIFFNESS);
	    
            dampIntDescr->SetPDEIds(this, this);
            assemble_->AddIntegrator(dampIntDescr, subdoms_[actSD]);
          }
	  
          else if  ( dampingList_[subdoms_[actSD]] == FRACTIONAL_BLANK ||
                     dampingList_[subdoms_[actSD]] == FRACTIONAL_BLANK_INT ) {

	    Double fracAlpha, fracExp;
	    materials_[subdoms_[actSD]]->GetScalar(fracAlpha,ACOU_ALPHA,REAL);
	    materials_[subdoms_[actSD]]->GetScalar(fracExp,FRACTIONAL_EXPONENT,REAL);

            coeffdamp =  - density * 2.0 * fracAlpha / c0 / sin((fracExp-1.0)*PI/2.0);
            // prefactor of blank alg
            coeffdamp *= exp(-gammaln(1.0- (fracExp- 1.0)) ); 
            // weight factor of index 0
            coeffdamp *= 1.0/(1.0- (fracExp- 1.0));           
            BaseForm * bilinearDamp  = 
              new MassInt(coeffdamp, dofspernode_, isaxi_);
            bilinearDamp->SetFracDamping();
            // formulation using DAMPING matrix
            // adapt NewmarkFracDamp::Init and StdPDE::GetFracDampMatrixCoeff
            // IntegratorDescriptor * dampIntDescr = 
            //   new IntegratorDescriptor(bilinearDamp, DAMPING);
	    
            // two matrices formulation
            // added to STIFFNEss matrix because, because 
            //   matrix_factors[STIFFNESS] = 1.0
            IntegratorDescriptor * dampIntDescr = 
              new IntegratorDescriptor(bilinearDamp, STIFFNESS);
	    
            dampIntDescr->SetPDEIds(this, this);
            assemble_->AddIntegrator(dampIntDescr, subdoms_[actSD]);

          }
        }

        // Finally add the stiffness/mass integrators
        assemble_->AddIntegrator(stiffIntDescr, subdoms_[actSD]);
        assemble_->AddIntegrator(massIntDescr, subdoms_[actSD]); 
      }
    }
  
    // **********************************************************************
    //   von Neumann boundary condition
    // **********************************************************************
    StdVector<RegionIdType> IdVec;
    ptgrid_->RegionNameToId( IdVec, bcs_ni_ );
    for (UInt actSD = 0; actSD < bcs_ni_.GetSize(); actSD++) {

      // we assume the surface normal points out of domain,
      //  but we want to take deflections into the domain positive
      Double amplitude= -1.0 * val_ni_[actSD];

      //BaseForm *neumannBC = new VolumeSrcInt( amplitude, isaxi_ );
      LinearSurfForm *neumannBC = new LinNeumannInt( amplitude, DENSITY, isaxi_ );
      neumannBC->SetVoluInfo( materials_ );

      if (analysistype_ == TRANSIENT ) {
        assemble_->AddRhsSrcSurfIntegrator( neumannBC,
                                            IdVec[actSD],
                                            fncnames_ni_[actSD],
                                            nonLin_ );
      }
      else if (analysistype_ == HARMONIC ) {
        assemble_->AddRhsSrcSurfIntegrator( neumannBC,
                                            IdVec[actSD],
                                            bcs_ni_phase_[actSD],
                                            nonLin_ );
      }
    }

    // **********************************************************************
    //   surface-integration: Absorbing boundaries
    // **********************************************************************
    if ( absorbingBCs_ == TRUE) { // && analysistype_ != HARMONIC ) {
      for (UInt actSD = 0; actSD < absBCs_.GetSize(); actSD++) {
        SurfForm * bilinear_damp = new AbsorbingBCsInt(isaxi_);
        bilinear_damp->SetFirstVoluInfo(pdename_, materials_ );

        // In the case of acou-mech coupling we have to multiply the 
        // abc-Integrator matrix with -1
        if ( isMechCoupled_ == TRUE ) {
          bilinear_damp->SetFactor(-1.0);
        }
        // In the case of acou-nrbc coupling we have to multiply the 
        // abc-Integrator matrix with C0 and multiply it by nrbcBeta0
        if ( isNrbcCoupled_ == TRUE ) {
          bilinear_damp->SetFactor(sqrt( compressibility / density ));
        }
        IntegratorDescriptor * abcDescr = 
          new IntegratorDescriptor(bilinear_damp, DAMPING);
        abcDescr->SetPDEIds(this, this);     
        //	std::cout << "Add ABC" << std::endl; 
        assemble_->AddSurfIntegrator(abcDescr, absBCs_[actSD]);
      }
    }
  }

  void AcousticPDE::DefineSolveStep() {
    ENTER_FCN( "AcousticPDE::DefineSolveStep" );

    if( bubbleDynType_ != NOBUBBLETYPE ) {
      solveStep_ = new SolveStepAcousticBubble(*this, bubbleDynType_);
    }
    else {
      solveStep_ = new SolveStepAcoustic(*this);
    }
  }


  // ======================================================
  // TIME STEPPING SECTION
  // ======================================================

  void AcousticPDE::InitTimeStepping() {
    ENTER_FCN( "AcousticPDE::InitTimeStepping" );

    UInt rhsSize = eqnData_->GetNumEQNs() *
      eqnData_->GetNumDofsPerEQN();

    // this includes rayleigh and thermoviscous damping
    if ( fracDamping_ == FALSE ) { 
      if ( effectiveMass_ == FALSE ) {
        TS_alg_ = new Newmark( algsys_, rhsSize );
      }
      else {
        TS_alg_ = new NewmarkEffMass( algsys_, rhsSize );
      }
    }
    else {
      if ( effectiveMass_ == FALSE )
        TS_alg_ = new NewmarkFracDamp( algsys_, rhsSize, 
                                       pdeId_, eqnData_, ptgrid_, this, 
                                       subdoms_, dampingList_ );
      else
        Error("This needs to be implemented!",__FILE__,__LINE__);
    }
  
    // Needed for fractional damping model, see Assemble::AssembleMatrices
    assemble_->SetPDEPointer(this);

    // Needed for nonlinear wave equation, get memory for linear part of RHS
    if ( nonLin_ )
      RhsLinVal_.Resize(numPDENodes_);

  }



  // =========================================================================
  // COUPLING SECTION
  // =========================================================================

  void AcousticPDE::InitCoupling(PDECoupling * Coupling) {
    
    ENTER_FCN( "AcousticPDE::InitCoupling" );
    
    isIterCoupled_ = TRUE;
    ptCoupling_   = Coupling;

    const UInt numCouplings = ptCoupling_->GetNumOutputCouplings();

    StdVector<StdVector<UInt> > elemNodeToCouplingNode_tmp;
    elemNodeToCouplingNode_.Resize(numCouplings);
    
    for (UInt i = 0; i < numCouplings; i++) {

      if (ptCoupling_->GetOutputQuantity(i) == ACOU_FORCE) {
        // Intialize the memory of the coupling values
        ptCoupling_->CreateCouplingVector(i,isComplex_);

        // now since we need a incremental formulation, 
        //  initialize some necessary vectors
        isIncrFormulation_ = TRUE;
      }

      else if (ptCoupling_->GetOutputQuantity(i) == ACOU_POT_NRBC)    {
        // Intialize the memory of the coupling values
        ptCoupling_->CreateCouplingVector(i,isComplex_);
        
        // now since we need a incremental formulation, 
        //  initialize some necessary vectors
        isIncrFormulation_ = TRUE;
      }

      else if (ptCoupling_->GetOutputQuantity(i) == ACOU_POWERDENSITY) {
        // Intialize the memory of the coupling values
        // DIRTY HACK ALARM!
        // =================
        // In case of transient-transient coupling coupling vector is Double.
        // In case of harmonic-transient coupling we also want a Double vector!
        ptCoupling_->CreateCouplingVector(i,FALSE);

        //get the element-node to coupling node matching
        StdVector<std::string> couplRegions;
        ptCoupling_->GetOutputRegions(i, couplRegions);
        StdVector<RegionIdType> regionIds;
        ptgrid_->RegionNameToId( regionIds, couplRegions );

        //Get total number of coupling elements
        UInt totalCouplingElems = ptgrid_->GetNumElems( regionIds );
        
        elemNodeToCouplingNode_tmp.Clear();
        elemNodeToCouplingNode_tmp.Resize(totalCouplingElems);

        UInt offset = 0;
        for ( UInt reg = 0; reg < couplRegions.GetSize(); reg++ ) {

          // find subdomain index = SDidx
          Integer SDidx = subdoms_.Find( regionIds[reg] );
          if (SDidx==-1) {
            (*error) << "AcousticPDE: Coupling region is not within the "
                     << "subdomains of the PDE!";
            Error( __FILE__, __LINE__ );
          }

          // get elements belonging to subdomain
          StdVector<Elem*> elemssd;
          ptgrid_->GetVolElems(elemssd, subdoms_[SDidx]);

          StdVector<UInt> * couplingnodes = NULL;
          ptCoupling_->GetOutputNodes(i, couplingnodes);
          if ( couplingnodes == NULL ) {
            (*error) << "The pointer 'couplingnodes' is NULL!";
            Error( __FILE__, __LINE__ );
          }

          for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {
            StdVector<UInt> & connecth = elemssd[actEl]->connect;
            elemNodeToCouplingNode_tmp[offset+actEl].Resize(connecth.GetSize());

            for ( UInt elnode = 0; elnode < connecth.GetSize(); elnode++ ) {
              for (UInt cnode=0; cnode<(*couplingnodes).GetSize(); cnode++) {

                if (connecth[elnode] == (*couplingnodes)[cnode] ) {
                  elemNodeToCouplingNode_tmp[offset+actEl][elnode] = cnode;
                  break;
                }
              }
            }
          }
          //in the case, that we have more than one coupling region!
          offset = elemssd.GetSize();
        }
        elemNodeToCouplingNode_[i]  = elemNodeToCouplingNode_tmp;
      }

    }
  }
  

  void AcousticPDE::CalcOutputCoupling() {

    ENTER_FCN( "AcousticPDE::CalcOutputCoupling" );

    UInt dof;
    SolutionType quantity;
    StdVector<Elem*> * couplingElems = NULL;
    StdVector<UInt> * couplingNodes = NULL;
    CFSVector * temp_values = NULL;
    UInt regionCount = 0;
  
    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == FALSE)
      return;

    // loop over all output coupling quantities
    for (UInt i=0; i<ptCoupling_->GetNumOutputCouplings(); i++) {

      quantity = ptCoupling_->GetOutputQuantity(i);
      ptCoupling_->GetOutputValues(i, temp_values);

      // hard coded cast, since coupling is only possible with
      // real valued entries

      Vector<Double> * values = dynamic_cast<Vector<Double>*>(temp_values);

      switch(ptCoupling_->GetOutputType(i)) {
        
      case NODE:
        ptCoupling_->GetOutputNodes(i, couplingNodes);
        if (quantity == ACOU_FORCE) {
          
          ptCoupling_->GetOutputElements(i, couplingElems);
          dof = ptCoupling_->GetOutputDof(i);
          CalcMechCouplingRHS(couplingElems, *couplingNodes, *values, dof);
        }
        else if (quantity == ACOU_POWERDENSITY) {
          if ( isComplex_ == FALSE ) {
            CalcHeatCouplingRHS<Double>( *values, 
                                         elemNodeToCouplingNode_[regionCount],
                                         i, couplingNodes->GetSize() );
          } 
          else {
            CalcHeatCouplingRHS<Complex>(*values, 
                                         elemNodeToCouplingNode_[regionCount],
                                         i, couplingNodes->GetSize() );
          }        
          regionCount++;
        }
        else if (quantity == ACOU_POT_NRBC) {
          ptCoupling_->GetOutputElements(i, couplingElems);
          dof = ptCoupling_->GetOutputDof(i);
          
          // Here this call gives the values to phi0 in nrbcPDE
          CalcNRBCCouplingRHS(couplingElems, *couplingNodes,
                              *values, dof);                              
        }
        break;

      case ELEM:
        Error("No Element coupling output", __FILE__,__LINE__);
      }
    }
  }


  void AcousticPDE::
  CalcMechCouplingRHS( StdVector<Elem*> * couplingElems, 
                       StdVector<UInt> & couplingNodes,
                       Vector<Double>& elemCouplingSols,
                       UInt couplingdof ) {
    
    ENTER_FCN( "AcousticPDE::CalcMechCouplingRHS" );
    
    Double density = 0.0;
    Double sign = 0.0;
    Integer matIndex = -1;
    Elem * ptVolElem = NULL;
    Matrix<Double> ptCoord, elemMat;
    Vector<Double> sol, forceOnElem, normal;
    
    elemCouplingSols.Init(0.0);
    
    for (UInt actElem=0; actElem<couplingElems->GetSize(); actElem++) {
      
      // Perform cast from volume element to surface element, since
      // mech-acou coupling makes only sense on surface elements
      SurfElem * actCoupleElem = 
        dynamic_cast<SurfElem*> ((*couplingElems)[actElem]);

      if (actCoupleElem == NULL) {
        Error( "No elements found for coupling!", __FILE__, __LINE__ );
      }
      
      BaseFE * ptElem = actCoupleElem->ptElem;
      StdVector<UInt> & connecth = actCoupleElem->connect;
      GetElemCoords(connecth, ptCoord);
      
      // Try to find according region for first neighbouring volume
      // element of the surface element
      matIndex = subdoms_.Find(actCoupleElem->ptVolElem1->regionId);
      
      // If first volume element does not belong to acoustic PDE, try the
      // second one
      if ( matIndex == -1 ) {
        matIndex = subdoms_.Find(actCoupleElem->ptVolElem2->regionId);
        ptVolElem = actCoupleElem->ptVolElem2;
        sign = actCoupleElem->normalSign;
      } else {
        ptVolElem = actCoupleElem->ptVolElem1;
        sign = -1.0 * actCoupleElem->normalSign;
      }
      
      if ( matIndex == -1) {
        (*error) << "AcousticPDE::CalcMechCouplingRHS: The two volume "
                 << "element neighbours of surface element Nr. "
                 << actCoupleElem->elemNum << " do not belong to my regions!";
        Error( __FILE__, __LINE__ );
      }
      
      // Assign correct density
      materials_[subdoms_[matIndex]]->GetScalar(density,DENSITY,REAL);
      
      BaseForm * bilinear_mass = new MassInt(ptElem, density, isaxi_);
      bilinear_mass->CalcElementMatrix(ptCoord, elemMat);
      delete bilinear_mass;     
      
      GetDerivSolVecOfElement(sol, connecth);
      
      forceOnElem = elemMat * sol;
      
      // force has to be added on RHS with negative sign
      forceOnElem *= - 1.0;
      
      // the normal vector points outwards of the MECHANICAL domain
      // (see. Kaltenbacher, "Num. Sim. of Mechatr. Act. & Sens." chapter 8.2)
      ptgrid_->CalcSurfNormal(normal, *actCoupleElem);
      normal *= sign;

      for (UInt actNode=0; actNode<ptCoord.GetSizeCol(); actNode++) {
        UInt nodePos = 0;
          
        while(connecth[actNode] != couplingNodes[nodePos] &&
              nodePos < couplingNodes.GetSize()) {
          nodePos++;
        }
          
        for (UInt actDof=0; actDof < couplingdof ; actDof++) {
          elemCouplingSols[nodePos*couplingdof +actDof] += 
            forceOnElem[actNode] * normal[actDof];
        }
      }
    }
  }
  
  
  void AcousticPDE::
  CalcNRBCCouplingRHS( StdVector<Elem*> * couplingElems, 
                       StdVector<UInt> & couplingNodes,
                       Vector<Double>& elemCouplingSols,
                       UInt couplingdof ) {
    
    ENTER_FCN( "AcousticPDE::CalcNRBCCouplingRHS" );

    Double density = 0.0;
    Double sign = 0.0;
    Double coeff_Pmass, coeff_Qstiff;
    UInt nrbcMatType=0; // 1=Rmat, 2=pmat, 3=Qmat
    Integer matIndex = -1;
    Elem * ptVolElem = NULL;
    Matrix<Double> ptCoord, elemMat;
    Vector<Double> sol, deriv2sol, Qmat_x_sol, Pmat_x_2derSol, normal;
    
    elemCouplingSols.Init(0.0);
    
    for (UInt actElem=0; actElem<couplingElems->GetSize(); actElem++) {
      
      // Perform cast from volume element to surface element, since
      // nrbc-acou coupling makes only sense on surface elements
      SurfElem * actSurfCoupleElem = 
        dynamic_cast<SurfElem*> ((*couplingElems)[actElem]);

      if (actSurfCoupleElem == NULL) {
        Error( "No elements found for coupling!", __FILE__, __LINE__ );
      }
      
      BaseFE * ptElem = actSurfCoupleElem->ptElem;
      StdVector<UInt> & connecth = actSurfCoupleElem->connect;
      GetElemCoords(connecth, ptCoord);
      
      // Try to find according region for first neighbouring volume
      // element of the surface element
      matIndex = subdoms_.Find(actSurfCoupleElem->ptVolElem1->regionId);
      
      // If first volume element does not belong to acoustic PDE, try the
      // second one
      if ( matIndex == -1 ) {
        matIndex = subdoms_.Find(actSurfCoupleElem->ptVolElem2->regionId);
        ptVolElem = actSurfCoupleElem->ptVolElem2;
        //        sign = actSurfCoupleElem->normalSign;
      } else {
        ptVolElem = actSurfCoupleElem->ptVolElem1;
        //        sign = -1.0 * actSurfCoupleElem->normalSign;
      }
      
      if ( matIndex == -1) {
        (*error) << "AcousticPDE::CalcNRBCCouplingRHS: The two volume "
                 << "element neighbours of surface element Nr. "
                 << actSurfCoupleElem->elemNum << " do not belong to my regions!";
        Error( __FILE__, __LINE__ );
      }
      
      // Density set to 1 since for this mass integrator no fac is needed
      density = 1.0;
      
      // !!!!!!!!!!!!!!!USE PTELEM!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      // FIRST WE COMPUTE THE "Qmat_x_sol" RHS term
      // Standard "mass" matrix shape so using default constructor
      coeff_Qstiff=1.0;
      nrbcMatType=3;
      NrbcInt * bilinear_Qstiff  = 
        new NrbcInt(ptElem, coeff_Qstiff, nrbcMatType, isaxi_);
      
      // Now we get the Q matrix with tangential derivatives 
      bilinear_Qstiff->CalcElementMatrix(ptCoord, elemMat);
      
      //std::cout<<"%%%%%%----ACOUSTIC-PDE COUPLING OUTPUT----%%%%%"<<std::endl;
      // std::cout<<"Qmatrix= "<<std::endl;
      // std::cout<<elemMat<<std::endl;
      
      delete  bilinear_Qstiff;
      GetSolVecOfElement(sol, connecth);
      //          for(UInt k=0; k<sol.GetSize(); k++)
      //           if (abs(sol[k])<=5e-16)
      //             sol[k]=0;
        
      // std::cout<<"AcouSolution= "<<std::endl;
      // std::cout<<sol<<std::endl;
      Qmat_x_sol = elemMat * sol;
      
      // std::cout<<"Qmat_x_sol= "<<std::endl;
      // std::cout<<Qmat_x_sol<<std::endl;
      //  Q_x_sol has to be added on RHS of nrbcPDE1 with negative sign
      
      Qmat_x_sol *= - 1.0;
      //END OF "Qmat_x_sol" Computation   
      
      //NOW WE COMPUTE "Pmat_x_Deriv2Sol" RHS term 
      coeff_Pmass=1.0;
      nrbcMatType=2;
      NrbcInt * bilinear_Pmass  = 
        new NrbcInt(ptElem, coeff_Pmass, nrbcMatType, isaxi_);
      
      bilinear_Pmass->CalcElementMatrix(ptCoord, elemMat);
        
      //std::cout<<"Pmatrix= "<<std::endl;
      // std::cout<<elemMat<<std::endl;
      
      delete  bilinear_Pmass;
      
      GetDeriv2SolVecOfElement(deriv2sol, connecth);
      
      // std::cout<<"Acou2ndDerivSolution= "<<std::endl;
      // std::cout<<deriv2sol<<std::endl;
      
      Pmat_x_2derSol = elemMat * deriv2sol;
      // Computation of alpha factor (like in nrbcPDE)
      Double c0, alphaNRBC, Cj, density, compressibility;
      //actSD is 0 since there is only one acoustic domain
      materials_.begin()->second->GetScalar( density, DENSITY, REAL );
      density=1;
      materials_.begin()->second->GetScalar( compressibility, ACOU_BULK_MODULUS, REAL );
      c0 = sqrt(compressibility/density);
      Cj = 1.0;
      alphaNRBC = (1/(Cj*Cj) - 1/(c0*c0));
      Pmat_x_2derSol *= alphaNRBC;
      
      // std::cout<<"Pmat_x_2derSol= "<<std::endl;
      // std::cout<<Pmat_x_2derSol<<std::endl;
      
      // END OF "Pmat_x_2derSol" Computation   
      
      for (UInt actNode=0; actNode<ptCoord.GetSizeRow(); actNode++) {
        UInt nodePos = 0;
        while(connecth[actNode] != couplingNodes[nodePos] &&
              nodePos < couplingNodes.GetSize()) {
          nodePos++;
        }
        
        elemCouplingSols[nodePos*couplingdof] = 
          (Qmat_x_sol[actNode] + Pmat_x_2derSol[actNode]);
        
        // std::cout <<"elemCouplingSols["<<(nodePos*couplingdof)<<"]= "
        //           <<elemCouplingSols[nodePos*couplingdof]<<std::endl;
        
      }
      // std::cout<<"elemCouplingSols = "<<elemCouplingSols<<std::endl;
      
      //std::cout <<"%%%--END ACOUSTIC-PDE COUPLING OUTPUT INFO----%%%%"
      //          <<std::endl;                  
      
    }
    }

  
  template <class TYPE>
  void AcousticPDE::
  CalcHeatCouplingRHS(Vector<Double> & sourceValue, 
                      StdVector<StdVector<UInt> > & elemNodeToCouplingNode,
                      UInt actCoupling,
                      UInt numCouplingNodes) {


    ENTER_FCN( "AcousticPDE::CalcHeatCouplingRHS" );
    
    // get the coupling regions
    StdVector<std::string> couplRegions;
    ptCoupling_->GetOutputRegions(actCoupling, couplRegions);
    StdVector<RegionIdType> regionIds;
    ptgrid_->RegionNameToId( regionIds, couplRegions );

    // Operator for calculating energy density
    AcouPowerDensityOp<TYPE> *SourceOp = 
      new AcouPowerDensityOp<TYPE>( ptgrid_, this, eqnData_, isaxi_ );

    // initialize output vector
    sourceValue.Init();

    Vector<Double> elemPowerDensity;
    
    UInt offset = 0;
    for (UInt reg=0; reg<couplRegions.GetSize(); reg++) {
      
      // find subdomain index
      Integer SDidx = subdoms_.Find( regionIds[reg] );
      Double density;
      materials_[regionIds[reg]]->GetScalar(density,DENSITY,REAL);
      
      // get elements belonging to subdomain
      StdVector<Elem*> elemssd;
      ptgrid_->GetVolElems(elemssd, subdoms_[SDidx]);
      
      for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {

        SourceOp->CalcElemPD(elemPowerDensity, elemssd[actEl], density);

        // Add the element energy to the according coupling node
        StdVector<UInt> connecth = elemssd[actEl]->connect;
        for (UInt elnode=0; elnode<connecth.GetSize(); elnode++) {

          sourceValue[elemNodeToCouplingNode[actEl+offset][elnode]]
            += elemPowerDensity[elnode];
        }
      }
      
      //in the case, that we have more than one coupling region!
      offset = elemssd.GetSize();
    }
  }


  Boolean AcousticPDE::HasOutput(SolutionType output) {
    ENTER_FCN( "AcousticPDE::HasOutput" );
    if ((output == ACOU_FORCE) || (output == ACOU_POWERDENSITY)) {
      return TRUE;
    }
    if (output == ACOU_POT_NRBC) {
      return TRUE;
    }
    return FALSE;
  }



  // =========================================================================
  // POSTPROCESSING SECTION
  // =========================================================================

  void AcousticPDE::PostProcess() {
    ENTER_FCN( "AcousticPDE::PostProcess" );

    StdVector<Elem*> saveElems, temp;
    for ( UInt i = 0; i < saveElemForceHist_.GetSize(); i++ ) {

      ptgrid_->GetElemsByName(temp, saveElemForceHist_[i]);
      for( UInt j = 0; j < temp.GetSize(); j++ ) {
        saveElems.Push_back(temp[j]);
      }
    }

    // compute force from pressure/acoustic potential at surface element
    if (saveElemForceHist_.GetSize() != 0 && isComplex_ == FALSE ) {
      CalcForce(saveElems);
    }

    // compute pressure from acoustic potential
    if (calcElemPressure_.GetSize() != 0) {
      if ( isComplex_ == FALSE) {
        CalcElemPressure<Double>();
      } 
      else {
        CalcElemPressure<Complex>();
      }
    }
     

    // Last but no least trigger postprocessing fromt within script-file
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

  void AcousticPDE::CalcForce( StdVector<Elem*> & saveElems ) {  
    ENTER_FCN( "AcousticPDE::CalcForce" );

    Matrix<Double> ptCoord;
    Vector<Double> forceVec(saveElems.GetSize());
	
    for (UInt actEl=0; actEl < saveElems.GetSize(); actEl++) {

      // Perform cast from volume element to surface element,
      //  since calculation of force makes only sense on a surface
      SurfElem * actSaveElem = 
        dynamic_cast<SurfElem*> ((saveElems)[actEl]);

      if (actSaveElem == NULL) {
        Error( "No elements specified in storeResults section!"
               , __FILE__, __LINE__ );
      }

      BaseFE * ptElem = actSaveElem->ptElem;
      StdVector<UInt> connect = actSaveElem->connect;
      GetElemCoords(connect, ptCoord);
                
      Vector<Double> valueElem;
      if ( formulation_ == ACOU_POTENTIAL ) {

        Integer matIndex = -1;
        // Try to find according region for first neighbouring volume
        // element of the surface element
        matIndex = subdoms_.Find(actSaveElem->ptVolElem1->regionId);

        // If first volume element does not belong to acoustic PDE, try the
        // second one
        Elem * ptVolElem = NULL;
        if ( matIndex == -1 ) {
          matIndex = subdoms_.Find(actSaveElem->ptVolElem2->regionId);
          ptVolElem = actSaveElem->ptVolElem2;
        }
        else {
          ptVolElem = actSaveElem->ptVolElem1;
        }
		
        if ( matIndex == -1) {
          (*error) << "AcousticPDE::CalcForce: The two volume element"
                   << "neighbours of surface element no. "
                   << actSaveElem->elemNum << " don't belong to my region!";
          Error( __FILE__, __LINE__ );
        }
      
        // Assign correct density
        Double density;
	materials_[subdoms_[matIndex]]->GetScalar(density,DENSITY,REAL);
		
        // retrieve 1st derivative, since F = rho * dpsi/dt * A
        GetDerivSolVecOfElement(valueElem, connect);

        valueElem *= density;
      }
      else if ( formulation_ == ACOU_PRESSURE ) {

        // retrieve solution, since F = p * A
        GetSolVecOfElement(valueElem,connect);
      }

      const UInt nrIntPts= ptElem->GetNumIntPoints();
      const Vector<Double> & intWeights = ptElem->GetIntWeights();

      Vector<Double> forceElem(1); // is defined as vector, because 
      // the method SetElemResult can only handle Vector<Double>
      forceElem[0] = 0;
      Double jacDet;
      for (UInt actIntPt=1; actIntPt<=nrIntPts;  actIntPt++) {

        jacDet = ptElem->CalcJacobianDetAtIp(actIntPt, ptCoord);
        Vector<Double> shapeFnc;  
        ptElem -> GetShFncAtIp(shapeFnc, actIntPt);

        if (isaxi_) {
          Vector<Double> coordAtIP;
          coordAtIP = ptCoord * shapeFnc;
          forceElem[0] +=  (intWeights[actIntPt-1] * jacDet
                            * 2 * PI) * coordAtIP[0] * (valueElem * shapeFnc);
        }
        else {
          forceElem[0] +=  intWeights[actIntPt-1] * jacDet 
            * (shapeFnc * valueElem);
        }
      }

      // map element result back in global set of results
      UInt pdeElem;
      pdeElem = eqnData_->Mesh2PDEElem(actSaveElem->elemNum);
      acouForce_.SetElemResult(pdeElem-1,forceElem);

      forceVec[actEl] = forceElem[0];
    }

    sumAcouForce_ = 0;
    for(UInt k=0; k<saveElems.GetSize(); k++) {
      sumAcouForce_ += forceVec[k];
    }
  }

  template <class TYPE>
  void AcousticPDE::CalcElemPressure( ) {  
    ENTER_FCN( "AcousticPDE::CalcElemPressure" );

    // loop over all subdomains
    for (UInt actSD=0; actSD<calcElemPressure_.GetSize(); actSD++) {

      // get all elements belonging to subdomain
      StdVector<Elem*> elemsSD;
      ptgrid_->GetVolElems( elemsSD, calcElemPressure_[actSD] );

      // density of elements in subdomain
      Double density;
      materials_[subdoms_[actSD]]->GetScalar(density,DENSITY,REAL);

      // loop over all elements of subdomain
      for (UInt actEl=0; actEl< elemsSD.GetSize(); actEl++) {

        BaseFE * ptElem = elemsSD[actEl]->ptElem;
        const UInt nrIntPts= ptElem->GetNumIntPoints();
        const Vector<Double> & intWeights = ptElem->GetIntWeights();

        // get element coordinates
        StdVector<UInt> & connect = elemsSD[actEl]->connect;
        Matrix<Double> ptCoord;
        GetElemCoords(connect, ptCoord);

        // get shape function at center of the element
        Vector<Double> shapeFnc;
        Vector<Double> LCoord;
        ptElem -> GetCoordMidPoint(LCoord);
        ptElem -> GetShFnc(shapeFnc,LCoord);

        // retrieve 1st derivative and multiply with density, 
        //  since p = rho * dpsi/dt
        Vector<TYPE> valueElem;
        GetDerivSolVecOfElement(valueElem, connect);

        valueElem *= density;
        Vector<TYPE> pressureElem(1);
        pressureElem[0] = valueElem * shapeFnc;

        // map element result back in global set of results
        UInt pdeElem;
        pdeElem = eqnData_->Mesh2PDEElem(elemsSD[actEl]->elemNum);

        if ( isComplex_ == TRUE ) {
          ElemStoreSol<Complex> & pressure = 
            dynamic_cast<ElemStoreSol<Complex>&>(*acouPressure_);
          pressure.SetElemResult(pdeElem-1,pressureElem);
        }
        else {
          ElemStoreSol<Double> & pressure = 
            dynamic_cast<ElemStoreSol<Double>&>(*acouPressure_);
          pressure.SetElemResult(pdeElem-1,pressureElem);
        }

      }
    }
  }


  void AcousticPDE::WriteResultsInFile(const UInt kstep,
                                       const Double asteptime,
                                       UInt stepOffset,
                                       Double timeOffset) {
    ENTER_FCN( "AcousticPDE::WriteResultsInFile" );

#ifdef PARALLEL //only one thread should write the output
    int commrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
    if (!commrank) {
#endif
      NodeStoreSol<Double> solIm_mesh;
      NodeStoreSol<Double> * solTransient;
      NodeStoreSol<Complex> * solHarmonic;
      
      
      Double actTime = asteptime + timeOffset;
      UInt actStep = kstep + stepOffset;
      
      if (analysistype_==HARMONIC) {
        solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);

        // Please remember, that in harmonic case actTime means indeed
        //  the actual frequency
        if (saveSol_)
          outFile_->WriteNodeSolutionHarmonic(*solHarmonic,  actStep, 
                                              actTime, complexFormat_);

        if (calcElemPressure_.GetSize() != 0 ) {
          ElemStoreSol<Complex> & pressureConverted = 
            dynamic_cast<ElemStoreSol<Complex>&>(*acouPressure_);
          outFile_->WriteElemSolutionHarmonic(pressureConverted, actStep,
                                              actTime, complexFormat_);
        }
      }
      else {  
        solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);

        if (saveSol_){
          outFile_->WriteNodeSolutionTransient(*solTransient,actStep,actTime);
        }

        if (saveDeriv1_) {
          solDeriv1_.SetAlgSysVector(getS1()); 
          outFile_->WriteNodeSolutionTransient(solDeriv1_, actStep, actTime);
        }

        if (saveDeriv2_) {
          solDeriv2_.SetAlgSysVector(getS2());
          outFile_->WriteNodeSolutionTransient(solDeriv2_, actStep, actTime);
        }

        if(m_bWriteSpecialBCs) {
          m_pGridAdaption->Add2DataFile(*solTransient, actTime);
        }

        if (saveRHSval_){
          outFile_->WriteNodeSolutionTransient(rhs_, actStep, actTime);
          if (plotRHSVel_==TRUE)
            outFile_->WriteNodeSolutionTransient(rhs2_, actStep, actTime);
        }

        if (calcElemPressure_.GetSize() != 0 ) {
          ElemStoreSol<Double> & pressureConverted = 
            dynamic_cast<ElemStoreSol<Double>&>(*acouPressure_);
          outFile_->WriteElemSolutionTransient(pressureConverted,actStep,actTime);
        }

      }
 
      // ------- for bubble results ----------------------
      if(bubbleDynType_ != NOBUBBLETYPE &&
         (saveSol_ == TRUE || saveDeriv1_ == TRUE || saveDeriv2_ == TRUE)) {

        StdVector<Double> radius;
        StdVector<Double> velocity;

        if ( isMechCoupled_ == TRUE ) {

          radius = dynamic_cast<SolveStepAcousticMechBubble*>
            (solveStep_)->GetResultData("bubbleRadius");

          velocity = dynamic_cast<SolveStepAcousticMechBubble*>
            (solveStep_)->GetResultData("bubbleVelocity");
        }
        else {

          radius = dynamic_cast<SolveStepAcousticBubble*>
            (solveStep_)->GetResultData("bubbleRadius");

          velocity = dynamic_cast<SolveStepAcousticBubble*>
            (solveStep_)->GetResultData("bubbleVelocity");
        }

        ElemStoreSol<Double> bubbleResult;
    
        // Resize solution arrays
        bubbleResult.SetNumSolutions(1);
        bubbleResult.SetSolutionType(ELEC_FIELD_INTENSITY);
        bubbleResult.SetNumElems(numElems_);
          
        // dimension hard coded for .unv file!
        bubbleResult.SetNumDofs(3);  
        bubbleResult.SetPtrEQNData(eqnData_, ptgrid_);
        bubbleResult.Init();
          
        for (UInt el=0; el<numElems_; el++) {
          Vector<Double> result(3);
            
          result[0] = radius[el];
          result[1] = velocity[el];
          result[2] = (4.0/3.0)*PI*bubbleDensity_*radius[el]*radius[el]
            *radius[el];
          //        if (el == 90)
          //std::cerr<<actTime<<"   " <<el<< "   " << result[0] << "   " 
          // result[1] << "     " << result[2] << std::endl;
            
          bubbleResult.SetElemResult(el,result);
        }
          
        outFile_->WriteElemSolutionTransient(bubbleResult, actStep, actTime);
      }

#ifdef PARALLEL
    }//!commrank
#endif
  }

  void AcousticPDE::WriteHistoryInFile(const UInt kstep,
                                       const Double asteptime,
                                       UInt stepOffset,
                                       Double timeOffset) {
    ENTER_FCN( "AcousticPDE::WriteHistoryInFile" );

#ifdef PARALLEL //only one thread should write the output
    int commrank;
    MPI_Comm_rank(MPI_COMM_WORLD,&commrank);
    if (!commrank) {
#endif
      NodeStoreSol<Double> solIm_mesh;
      NodeStoreSol<Double> * solTransient;
      NodeStoreSol<Complex> * solHarmonic;
      
      
      Double actTime = asteptime + timeOffset;
      UInt actStep = kstep + stepOffset;
      
      if (analysistype_==HARMONIC) {
        solHarmonic = dynamic_cast<NodeStoreSol<Complex>*>(sol_);

        if (saveSolHist_)
          outFile_->WriteNodeHistoryHarmonic(*solHarmonic,  actStep, 
                                             actTime, complexFormat_);

        if (saveDeriv1Hist_) {
          // multiply solution with j * omega
        }
      }
      else {  
        solTransient = dynamic_cast<NodeStoreSol<Double>*>(sol_);

        if (saveSolHist_){
          outFile_->WriteNodeHistoryTransient(*solTransient, actStep,actTime);
        }

        if (saveDeriv1Hist_) {
          solDeriv1_.SetAlgSysVector(getS1()); 
          outFile_->WriteNodeHistoryTransient(solDeriv1_, actStep, actTime);
        }

        if (saveDeriv2Hist_){
          solDeriv2_.SetAlgSysVector(getS2());
          outFile_->WriteNodeHistoryTransient(solDeriv2_, actStep, actTime);
        }

        if (saveRHSvalHist_){
          outFile_->WriteNodeHistoryTransient(rhs_, actStep, actTime); 
          if (plotRHSVel_==TRUE)
            outFile_->WriteNodeHistoryTransient(rhs2_, actStep, actTime); 
        }

        if (saveElemForceHist_.GetSize() != 0) {
          outFile_->WriteElemHistoryTransient(acouForce_, actStep, actTime);

          // lazy mans way to produce output
          std::cerr << actTime << "   " << sumAcouForce_ << std::endl;
        }

        if (saveElemPressureHist_.GetSize() != 0) {
          ElemStoreSol<Double> & pressureConverted = 
            dynamic_cast<ElemStoreSol<Double>&>(*acouPressure_);
          outFile_->WriteElemHistoryTransient(pressureConverted, actStep, 
                                              actTime);
        }

      }
 
#ifdef PARALLEL
    }//!commrank
#endif
  }


  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void AcousticPDE::ReadStoreResults() {  
    ENTER_FCN( "AcousticPDE::ReadStoreResults" );
    
    // Construct vectors for restricted parameter search
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    std::string quantity;
    
    // *****************************
    // Determine nodal results
    // ***************************** 
    StdVector<std::string> nodeValues;
    keyVec  = pdename_, "storeResults", "nodeResults", "region";
    attrVec = "", "", "type";  

    if ( formulation_ == ACOU_POTENTIAL ) {
      // --- pressure ---
      Enum2String(ACOU_PRESSURE, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
        (*error) << "It makes no sense to have a PDE in acoustic potential "
                 << "and retrieve \n node results in pressure. "
                 << "Try element results!";
        Error( __FILE__, __LINE__ );
      }

      // --- acoustic potential ---
      Enum2String(ACOU_POTENTIAL, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
        saveSol_ = TRUE;
        hasOutput_ = TRUE;
      }

      // --- acoustic rhsval ---
      Enum2String(ACOU_RHSVAL, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
        saveRHSval_ = TRUE;
        hasOutput_ = TRUE;
      }


      // --- acoustic potential, 1. Deriv ---
      Enum2String(ACOU_POTENTIAL_DERIV_1, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0 && analysistype_ != HARMONIC ) {
        saveDeriv1_ = TRUE;
        hasOutput_ = TRUE;
        // intialize corresponding storesolution object
        solDeriv1_.SetNumSolutions(1);
        solDeriv1_.SetNumNodes(numPDENodes_);
        solDeriv1_.SetSolutionType(ACOU_POTENTIAL_DERIV_1);
        solDeriv1_.SetNumDofs(1);
        solDeriv1_.SetPtrEQNData(eqnData_, ptgrid_); 
        solDeriv1_.Init(0);
      }

      // --- acoustic potential, 2. Deriv ---
      Enum2String(ACOU_POTENTIAL_DERIV_2, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0 && analysistype_ != HARMONIC) {
        saveDeriv2_ = TRUE;
        hasOutput_ = TRUE;

        solDeriv2_.SetNumSolutions(1);
        solDeriv2_.SetNumNodes(numPDENodes_);
        solDeriv2_.SetSolutionType(ACOU_POTENTIAL_DERIV_2);
        solDeriv2_.SetNumDofs(1);
        solDeriv2_.SetPtrEQNData(eqnData_, ptgrid_); 
        solDeriv2_.Init(0);
      }
    }
    else if ( formulation_ == ACOU_PRESSURE ) {
      // --- acoustic potential ---
      Enum2String(ACOU_POTENTIAL, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
        (*error) << "It acoustic pressure formulation, \n" 
                 << "retrieving acoustic potential requires integration.\n"
                 << "So forget it!";
        Error(__FILE__,__LINE__);
      }

      // --- pressure ---
      Enum2String(ACOU_PRESSURE, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
        std::string warnMsg;
        warnMsg = "Due to the restrictions in the .unv file format, the ";
        warnMsg += "acoustic pressure is written as acoustic (fluid) ";
        warnMsg += "potential!";
        Warning(warnMsg.c_str(), __FILE__, __LINE__);

        saveSol_ = TRUE;
        hasOutput_ = TRUE;

        // Note: If errors occur using the acoustic pressure formulation,
        // just comment out the following three lines
        sol_->SetSolutionType(ACOU_PRESSURE);
        sol_->SetNumDofs(1);
        sol_->Init();

      }
    }
    
    
    // *****************************
    // Determine element results
    // *****************************
    StdVector<std::string> regionNames;
    keyVec  = pdename_, "storeResults", "elemResults", "region";
    attrVec = "", "", "type";

    if ( formulation_ == ACOU_POTENTIAL ) {
      // pressure
      Enum2String(ACOU_PRESSURE, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, regionNames);
      ptgrid_->RegionNameToId( calcElemPressure_, regionNames );

      // If the the symbolic name is "all" compute quatity for all regions
      if ( calcElemPressure_.GetSize() == 1 &&
           calcElemPressure_[0] == ALL_REGIONS ) {
        calcElemPressure_ = subdoms_;
      }

      if (calcElemPressure_.GetSize() > 0) {

        hasOutput_ = TRUE;
        Info->PrintF( pdename_, " Computing acoustic pressure for regions:\n" );
        for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
          Info->PrintF( pdename_, " %s\n", regionNames[k].c_str() );
        }

        if (analysistype_ == HARMONIC ) {
          acouPressure_ = new ElemStoreSol<Complex>;
        } 
        else {
          acouPressure_ = new ElemStoreSol<Double>;
        }

        acouPressure_->SetNumSolutions(1);
        acouPressure_->SetSolutionType(ACOU_PRESSURE);
        acouPressure_->SetNumElems(numElems_);
        acouPressure_->SetNumDofs(1);
        acouPressure_->SetPtrEQNData(eqnData_, ptgrid_);
        acouPressure_->Init();
      }
    } 


    // *****************************
    // Determine nodal history
    // *****************************
    StdVector<std::string> saveNodeHist; 
    keyVec  = pdename_, "storeResults", "nodeHistory", "saveNodes";
    attrVec = "", "", "type";

    if ( formulation_ == ACOU_POTENTIAL ) {
      // pressure
      Enum2String(ACOU_PRESSURE, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist);
      if (saveNodeHist.GetSize() > 0) {
        (*error) <<  "In acoustic potential formulation, \n"
                 << " the quantity pressure is a element result.";
        Error( __FILE__, __LINE__ );
      }

      // --- acoustic potential  ---
      Enum2String(ACOU_POTENTIAL, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
      if (saveNodeHist.GetSize() > 0) {
        saveSolHist_ = TRUE;
        hasOutput_ = TRUE;
        Info->PrintF( pdename_, "Saving acouPotential for Nodes:\n" );
        for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
          Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
        }
      }

      // --- acoustic potential, 1. Deriv ---
      Enum2String(ACOU_POTENTIAL_DERIV_1, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
      if (saveNodeHist.GetSize() > 0) {
        saveDeriv1Hist_ = TRUE;
        hasOutput_ = TRUE;
        Info->PrintF( pdename_, "Saving acouPotentialD1 for Nodes:\n" );
        for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
          Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
        }
        // Check if solDeriv_1 was already set by the previous nodal
        // output check
        if( saveDeriv1_ != TRUE ) {
          solDeriv1_.SetNumSolutions(1);
          solDeriv1_.SetNumNodes(numPDENodes_);
          solDeriv1_.SetSolutionType(ACOU_POTENTIAL_DERIV_1);
          solDeriv1_.SetNumDofs(1);
          solDeriv1_.SetPtrEQNData(eqnData_, ptgrid_); 
          solDeriv1_.Init(0);
        }
      }

      // --- acoustic potential, 2. Deriv ---
      Enum2String(ACOU_POTENTIAL_DERIV_2, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
      if (saveNodeHist.GetSize() > 0) {
        saveDeriv2Hist_ = TRUE;
        hasOutput_ = TRUE;
        Info->PrintF( pdename_, "Saving acouPotetentialD2 for Nodes:\n" );
        for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
          Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
        }
        // Check if solDeriv_2 was already set by the previous nodal
        // output check
        if( saveDeriv2_ != TRUE ) {
          solDeriv2_.SetNumSolutions(1);
          solDeriv2_.SetNumNodes(numPDENodes_);
          solDeriv2_.SetSolutionType(ACOU_POTENTIAL_DERIV_2);
          solDeriv2_.SetNumDofs(1);
          solDeriv2_.SetPtrEQNData(eqnData_, ptgrid_); 
          solDeriv2_.Init(0);
        }
      }

      // --- acoustic RHS ---
      Enum2String(ACOU_RHSVAL, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
      if (saveNodeHist.GetSize() > 0) {
        saveRHSvalHist_ = TRUE;
        hasOutput_ = TRUE;
        Info->PrintF( pdename_, "Saving acouRHSval for Nodes:\n" );
        for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
          Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );
        }
      }
    }

    else if ( formulation_ == ACOU_PRESSURE ) {
      // --- acoustic potential ---
      Enum2String(ACOU_POTENTIAL, quantity);
      valVec = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, nodeValues);
      if (nodeValues.GetSize() > 0) {
        (*error) << "It acoustic pressure formulation, \n" 
                 << "retrieving acoustic potential requires integration.\n"
                 << "So forget it!";
        Error(__FILE__,__LINE__);
      }

      // --- pressure ---
      Enum2String(ACOU_PRESSURE, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveNodeHist );
    
      if (saveNodeHist.GetSize() > 0) {
        saveSolHist_ = TRUE;
        hasOutput_ = TRUE;
        Info->PrintF( pdename_, "Saving acouPotential for Nodes:\n" );
        for ( UInt k = 0; k < saveNodeHist.GetSize(); k++ ) {
          Info->PrintF( pdename_, "  %s\n", saveNodeHist[k].c_str() );

          sol_->SetSolutionType(ACOU_PRESSURE);
          sol_->SetNumDofs(1);
          sol_->Init();
        }
      }
    }

    // *****************************
    // Determine element history
    // *****************************
    keyVec  = pdename_, "storeResults", "elemHistory", "saveElems";
    attrVec = "", "", "type";

    if ( (formulation_==ACOU_POTENTIAL) || (formulation_==ACOU_PRESSURE) ) {

      // --- force ---
      Enum2String(ACOU_FORCE, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveElemForceHist_ );

      if (saveElemForceHist_.GetSize() > 0) {
        hasOutput_ = TRUE;
        acouForce_.SetNumSolutions(1);
        acouForce_.SetSolutionType(ACOU_FORCE);
        acouForce_.SetNumElems(numElems_);
        acouForce_.SetNumDofs(1);
        acouForce_.SetPtrEQNData(eqnData_, ptgrid_);
        acouForce_.Init();

        Info->PrintF( pdename_, "Saving acouForce for Elements:\n" );
        for ( UInt k = 0; k < saveElemForceHist_.GetSize(); k++ ) {
          Info->PrintF( pdename_, "  %s\n", saveElemForceHist_[k].c_str() );
        }
      }
    }
    if ( formulation_==ACOU_POTENTIAL ) {

      // --- pressure ---
      Enum2String(ACOU_PRESSURE, quantity);
      valVec  = "", "", quantity;
      params->GetList( keyVec, attrVec, valVec, saveElemPressureHist_ );

      if (saveElemPressureHist_.GetSize() > 0) {
        hasOutput_ = TRUE;
        acouPressure_->SetNumSolutions(1);
        acouPressure_->SetSolutionType(ACOU_PRESSURE);
        acouPressure_->SetNumElems(numElems_);
        acouPressure_->SetNumDofs(1);
        acouPressure_->SetPtrEQNData(eqnData_, ptgrid_);
        acouPressure_->Init();

        Info->PrintF( pdename_, "Saving acouPressure for Elements:\n" );
        for ( UInt k = 0; k < saveElemPressureHist_.GetSize(); k++ ) {
          Info->PrintF( pdename_, "  %s\n", saveElemPressureHist_[k].c_str() );
        }
      }
    }


  }

  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void AcousticPDE::ReadDataPML(std::string& dampingTypePML, 
                                Matrix<Double>& inner, 
                                Double& dampPML, 
                                RegionIdType actRegion) {
  
    ENTER_FCN( "AcousticPDE::ReadDataPML" );

    inner.Resize(dim_,dim_);

    // help variables for parameter checking
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    StdVector<std::string> stringVal;

    // Construct vectors for restricted parameter search
    std::string actRegionName;
    actRegionName = ptgrid_->RegionIdToName( actRegion );
    keyVec = "acoustic" , "region" , "pml" , "type";
    attrVec= ""         , "name"   , "";
    valVec = ""         , actRegionName, "";
    params->GetList( keyVec, attrVec, valVec, stringVal);
    dampingTypePML = stringVal[0];

    StdVector<Double> val;
 
    //xMin
    keyVec = "acoustic" , "region" , "pml" , "xMin";
    params->GetList( keyVec, attrVec, valVec, val);
    inner[0][0] =  val[0];

    //xMax
    keyVec = "acoustic" , "region" , "pml" , "xMax";
    params->GetList( keyVec, attrVec, valVec, val);
    inner[1][0] =  val[0];

    //yMin
    keyVec = "acoustic" , "region" , "pml" , "yMin";
    params->GetList( keyVec, attrVec, valVec, val);
    inner[0][1] =  val[0];

    //yMax
    keyVec = "acoustic" , "region" , "pml" , "yMax";
    params->GetList( keyVec, attrVec, valVec, val);
    inner[1][1] =  val[0];

    if ( dim_ == 3 ) {
      //zMin
      keyVec = "acoustic" , "region" , "pml" , "zMin";
      params->GetList( keyVec, attrVec, valVec, val);
      inner[0][2] =  val[0];

      //zMax     
      keyVec = "acoustic" , "region" , "pml" , "zMax";
      params->GetList( keyVec, attrVec, valVec, val);
      inner[1][2] =  val[0];
    }

    keyVec = "acoustic" , "region" , "pml" , "dampFactor";
    params->GetList( keyVec, attrVec, valVec, val);
    dampPML = val[0];


  }

  // ***********************************************************************
  //   Obtain information on desired output quantities from parameter file
  // ***********************************************************************
  void AcousticPDE::GetPMLLayerData(Matrix<Double>& inner, 
                                    Matrix<Double>& outer,
                                    UInt actSD)  {  

    ENTER_FCN( "AcousticPDE::GetPMLLayerData" );

    outer.Resize(dim_,dim_);
    outer = inner;

    // inner/outer:   xmin  ymin  zmin
    //                xmax  ymax  zmax


    Double minXPML, minYPML, minZPML, maxXPML, maxYPML, maxZPML;

    StdVector<Elem*> elemssd;
    ptgrid_->GetElems(elemssd, subdoms_[actSD]);

    for (UInt actEl=0; actEl< elemssd.GetSize(); actEl++) {
      BaseFE * ptEl = elemssd[actEl]->ptElem;
      StdVector<UInt> connecth = elemssd[actEl]->connect;
                             
      Matrix<Double> ptCoord;
      GetElemCoords(connecth, ptCoord);
      for (UInt i=0; i< ptCoord.GetSizeCol(); i++) {
        //minXPML
        if ( ptCoord[0][i] < outer[0][0] )
          outer[0][0] = ptCoord[0][i];

        //minYPML
        if ( ptCoord[1][i] < outer[0][1] )
          outer[0][1] = ptCoord[1][i];

        if (inner.GetSizeCol() > 2 ) {
          //minZPML
          if ( ptCoord[2][i] < outer[0][2] )
            outer[0][2] = ptCoord[2][i];
        }

        //maxXPML
        if ( ptCoord[0][i] > outer[1][0] )
          outer[1][0] = ptCoord[0][i];

        //maxYPML
        if ( ptCoord[1][i] > outer[1][1] )
          outer[1][1] = ptCoord[1][i];

        if (inner.GetSizeCol() > 2 ) {
          //maxZPML
          if ( ptCoord[2][i] > outer[1][2] )
            outer[1][2] = ptCoord[2][i];
        }
      }

    }

    //   std::cout << "outer:\n" << outer << std::endl;

  }

} // end of namespace
