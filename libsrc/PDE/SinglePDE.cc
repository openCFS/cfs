#include "PDE/SinglePDE.hh"

// for coordinate handling
#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"

// header for Paramhandling
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/ParamHandling/CFSOLASParams.hh"

// header for Materialhandling
#include "DataInOut/LoadMaterialData.hh"
#include "DataInOut/LoadMaterialDataFile.hh"
#ifdef USE_DATABASE
#include "DataInOut/LoadMaterialDataDatabase.hh"
#endif

// header for equation numbering
#include "blocknodeEQN.hh"
#include "scalarblockEQN.hh"
#include "scalarnodeEQN.hh"

// header for Solvestep
#include "Driver/stdSolveStep.hh"

// header for iterative coupling
#include "CoupledPDE/pdecoupling.hh"

#ifdef PARALLEL
#include <mpi.h>
#endif

namespace CoupledField {


  SinglePDE::SinglePDE( Grid *aptgrid, WriteResults *aOutFile,
                        TimeFunc *aptTimeFunc )
    :  StdPDE(aptgrid, aOutFile, aptTimeFunc) {
  
    ENTER_FCN( "BasePDE::BasePDE" );
  
    nonLin_ = FALSE;
    incStopCrit_ = 1e-2;
    residualStopCrit_ = 1e-3;
  
    // =====================================================================
    // set file pointers
    // =====================================================================
    ptTimeFunc_ = aptTimeFunc;
    assemble_   = NULL;
    algsys_     = NULL;
    eqnData_    = NULL;
  
    // =====================================================================
    // set analysis parameters
    // =====================================================================
    complexFormat_ = AMPLITUDE_PHASE; // or REAL_IMAG
    couplingBCsCounter_ = 0;
    numDirichletBCs_ = 0;
    updateCouplingBCs_ = FALSE;
    dim_ = ptgrid_->GetDim();
    geoUpdate_ = FALSE;
    iterCoupledCounter_ = 0;
    effectiveMass_ = FALSE;

    // =====================================================================
    // set solver parameters
    // =====================================================================
    eps_         = 1.0e-8;
    dampiter_    = 1.0;
    maxnumiter_  = 500;
    coarsealpha_ = 0.1; // in acousticPDE.cc set to 0.01
  
    // savederiv1_ = FALSE;
    // savederiv2_ = FALSE;

    // =====================================================================
    // set postprocessing parameters
    // =====================================================================
    hasOutput_      = FALSE;
    saveSol_        = FALSE;
    saveDeriv1_     = FALSE;
    saveDeriv2_     = FALSE;
    saveSolHist_    = FALSE;
    saveDeriv1Hist_ = FALSE;
    saveDeriv2Hist_ = FALSE;


    // =====================================================================
    // set miscellaneous parameters
    // =====================================================================
    pdeId_ = NO_PDE_ID;
    isDirectCoupled_  = FALSE;
    m_bReadSpecialBCs = FALSE;
			
  }


  // **********************
  //   Default Destructor
  // **********************
  SinglePDE::~SinglePDE() {

    ENTER_FCN( "SinglePDE::~SinglePDE" );

    // Delete algebraic system only if
    // PDE is not direct coupled
    if ( isDirectCoupled_ == FALSE ) {
      DeleteAlgSys();
      delete solveStep_;
      delete TS_alg_;
    }


    delete assemble_;
    delete sol_;
    delete solVec_;
    delete eqnData_;

    delete[] materialData_;
  }


  void SinglePDE::Init( UInt bcSequenceIndex,
                        std::string  bcSequenceTag) {
    ENTER_FCN( "SinglePDE::Init()" );

    StdVector<RegionIdType> allIDs;

    bcSequenceIndex_ = bcSequenceIndex;
    bcSequenceTag_ = bcSequenceTag;
  
    // =====================================================================
    // get regions/subdomains for PDE
    // =====================================================================
    StdVector<std::string> regionNames, surfaceNames;
    StdVector<RegionIdType> surfIds;
    params->GetList( "name", regionNames, pdename_, "region" );
    ptgrid_->RegionNameToId( subdoms_, regionNames );

    params->GetList( "name", surfaceNames, pdename_, "surface" );
    ptgrid_->RegionNameToId( surfIds, surfaceNames );

    // Create vector of all IDs
    allIDs = subdoms_;
    for (UInt i = 0; i < surfIds.GetSize(); i++) {
      allIDs.Push_back(surfIds[i]);
    }


    Info->PrintF( pdename_, "The %s PDE lives on the following regions:\n",
                  pdename_.c_str());
    for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
      Info->PrintF( pdename_, "%s, index %i\n",
                    regionNames[k].c_str(), subdoms_[k] );
    }
    Info->PrintF( "", "\n" );

    // Generate a fitting algebraic system only if PDE is NOT
    // direct coupled
    if ( isDirectCoupled_ == FALSE ) {
	  algsys_ = new StandardSystem();
	}

    // Get parameter and report object of OLAS
    olasParams_ = algsys_->GetOLASParams();
    olasReport_ = algsys_->GetOLASReport();

    // Obtain unique pde identifier
    pdeId_ = algsys_->ObtainPDEId( pdename_ );

    // Determine, if this is a parallel run
    // and pass this information to OLAS
    bool parallel = false;

#ifdef PARALLEL

    // If more than one process is running,
    // then we consider this a parallel run
    int commsize;
    MPI_Comm_size( MPI_COMM_WORLD, &commsize );
    parallel = ( commsize > 1 );

#endif

    olasParams_->SetValue( "Parallel", parallel );

    // =====================================================================
    // Get type of analysis
    // =====================================================================
    
    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
    std::string stepString;
    std::string analysis;

    params->Get( "type", analysis, "analysis" );

    AnalysisType analysisHelp;
    StdVector<std::string> tags, analysisTypes, pdenames;
    String2Enum(analysis,analysisHelp);

    // NOTE: The concept of isAlwaysStatic bites with Direct Coupling
    //       and must be re-designed
    if ( isAlwaysStatic_ == TRUE &&
         analysisHelp == TRANSIENT ) {
      analysisHelp = STATIC;
    }

    if ( analysisHelp == STATIC ) {
      isComplex_ = FALSE;
      assemble_ = new StaticAssemble(algsys_, ptgrid_);
      analysistype_ = STATIC;
    }

    else if ( analysisHelp == TRANSIENT ) {
      isComplex_ = FALSE;
      assemble_ = new TransientAssemble(algsys_, ptgrid_);
      analysistype_ = TRANSIENT;
    }

    else if ( analysisHelp == TRANSIENT4SLICE ) {
      isComplex_ = FALSE;
      assemble_ = new TransientAssemble(algsys_, ptgrid_);
      analysistype_ = TRANSIENT;
    }

    else if ( analysis=="harmonic" || analysis == "paramIdent" ) {
      isComplex_ = TRUE;
      assemble_ = new HarmonicAssemble(algsys_, ptgrid_);
      analysistype_ = HARMONIC;
    }

    else if (analysis == "multiHarmonic"){
      isComplex_ = TRUE;
      assemble_ = new MHassemble(algsys_,ptgrid_);
      analysistype_ = MULTIHARMONIC;
    }

    else if ( analysisHelp == MULTI_SEQUENCE ) {

      stepString = Info->GenStr(bcSequenceIndex_);
      attrVec = "", "index", "type";
      valVec = "", stepString, pdename_;

      keyVec = "multiSequence", "step", "pde", "analysis";
      params->Get(keyVec, attrVec, valVec, analysis);


      String2Enum(analysis, analysistype_);

      if ( isAlwaysStatic_ == TRUE &&
           analysistype_ == TRANSIENT ) {
        analysistype_ = STATIC;
      }      

      if ( analysistype_ == STATIC ) {
        assemble_ = new StaticAssemble(algsys_, ptgrid_);
        isComplex_ = FALSE;
      }
      else if ( analysistype_ == TRANSIENT ) {
        isComplex_ = FALSE;
        assemble_ = new TransientAssemble(algsys_, ptgrid_);
      }
      else if ( analysis=="harmonic" ) {
        isComplex_ = TRUE;
        assemble_ = new HarmonicAssemble(algsys_, ptgrid_);
        analysistype_ = HARMONIC;
      }
      else if ( analysis=="multiHarmonic" ) {
        isComplex_ = TRUE;
        assemble_ = new MHassemble(algsys_, ptgrid_);
        analysistype_ = MULTIHARMONIC;
      }
      else {
        (*error) << "SinglePDE::Init: AnalysisType '" << analysis
                 << "' is not supported";
        Error( __FILE__, __LINE__ );
      }
    }
    else {
      (*error) << "SinglePDE::Init: AnalysisType '" << analysisHelp
               << "' is not supported";
      Error( __FILE__, __LINE__ );
    }
    
    // Give pde Id to assemble
    assemble_->SetPDEId ( pdeId_ );
    
    // Determine if solution is of complex type or not
    if ( analysistype_ == HARMONIC || analysistype_ == MULTIHARMONIC ) {
      sol_ = new NodeStoreSol<Complex>;
      solVec_ = new Vector<Complex>;
    }
    else {
      sol_ = new NodeStoreSol<Double>;
      solVec_ = new Vector<Double>;
    }

    // =====================================================================
    // initialize adaptivity
    // =====================================================================
    if( params->IsSet( "adaptspace" ) ) {
      params->Get( "tolerance_space_error", tolSpaceErr_ );
    }
    else {
      tolSpaceErr_ = 0;
    }
    
    // =====================================================================
    // read in damping information
    // =====================================================================
      ReadDampingInformation( ptgrid_ );
      
    // =====================================================================
    // read in boundary conditions
    // =====================================================================
    ReadBCs();
    ReadSpecialBCs();

    // =====================================================================
    // read in NonLinearities
    // =====================================================================
    InitNonLin();

    // =====================================================================
    // initialize EQN-object and Storeresults class
    // =====================================================================

    // What type of equation numbering does the user want?
    std::string typeOfNumbering;
    keyVec  = "linearSystems", "system", "setup", "eqnNumbering";
    attrVec = "", "name", "";
    valVec  = "", pdename_, "";
    params->Get( keyVec, attrVec, valVec, typeOfNumbering );

    // How do we want to treat inhomogeneous Dirichlet boundary conditions?
    bool usePenalty = true;
    {
      std::string aux;
      keyVec  = "linearSystems", "system", "setup", "idbcHandling";
      attrVec = "", "name", "";
      valVec  = "", pdename_, "";
      params->Get( keyVec, attrVec, valVec, aux );
      usePenalty = aux == "penalty" ? true : false;
      Info->PrintF( pdename_, "Treating IDBCs using '%s' approach\n",
                    aux.c_str() );
    }

    // Assemble a system matrix with scalar complex or double entries
    if ( typeOfNumbering == "scalar" ) {
      if ( dofspernode_ == 1 ) {
        eqnData_ = new ScalarNodeEQN( ptgrid_, allIDs, dofspernode_,
                                      !usePenalty );
      }
      else {
        eqnData_ = new ScalarBlockEQN( ptgrid_, allIDs, dofspernode_,
                                       !usePenalty );
      }
      Info->PrintF( pdename_, "Using scalar equation numbering\n" );
    }

    // Treat all dofs of a node together and assemble a system matrix with
    // small square matrices as entries
    else if ( typeOfNumbering == "block" ) {
      if ( dofspernode_ == 1 ) {
        (*warning) << "dofspernode = " << dofspernode_
                   << "so 'block' numbering identical to 'scalar'";
        Warning( __FILE__, __LINE__ );
        eqnData_ = new ScalarNodeEQN( ptgrid_, allIDs, dofspernode_,
                                      usePenalty );
        Info->PrintF( pdename_, "Using scalar equation numbering\n" );
      }
      else {
        eqnData_ = new BlockNodeEQN( ptgrid_, allIDs, dofspernode_,
                                     usePenalty );
        Info->PrintF( pdename_, "Using block equation numbering\n" );
      }
    }

    // Build in Dirichlet boundary conditions and compute the mapping
    // which relates nodes/dofs to equation numbers
    eqnData_->SetHomoDirichletBCs ( bcs_hd_, homDirichDof_  );
    eqnData_->SetInhomDirichletBCs( bcs_id_, inhomDirichDof_);
    SETPROFILE("Before Equation Mapping");
    eqnData_->CalcMapping();
    SETPROFILE("After Equation Mapping");

    // Report results to logfile
    Info->PrintF( pdename_, "Linear system will have %d equations\n\n",
                  eqnData_->GetNumEQNs() );

#ifdef DEBUG
    eqnData_->Print(*debug);
#endif

    numPDENodes_ = eqnData_->GetNumLocalNodes();
    numElems_ = eqnData_->GetNumLocalElems();

    // Initialize Storesolution class
    sol_->SetNumSolutions(solTypes_.GetSize());
    sol_->SetNumNodes(numPDENodes_);
    for (UInt iSol=0; iSol<solTypes_.GetSize(); iSol++) {
      sol_->SetSolutionType(solTypes_[iSol],iSol);
      sol_->SetNumDofs(solDofs_[iSol], solTypes_[iSol]);
    }
    sol_->SetPtrEQNData(eqnData_, ptgrid_);
    sol_->Init(); 

    SETPROFILE("Before Resizing StoreSol");
    solVec_->Resize( eqnData_->GetNumEQNs() * eqnData_->GetNumDofsPerEQN() );
    SETPROFILE("After Resizing StoreSol");


    // =====================================================================
    // initialize assemble object
    // =====================================================================
    assemble_->SetPtr2EQNData(eqnData_); 
    assemble_->SetPtr2TimeFnc(ptTimeFunc_);

    if (pdename_ == "piezo" || pdename_ == "mechanic" ) {
      assemble_->SetGeneralParams(pdename_, dofspernode_, 
                                  subdoms_, pressSurf_, bcSequenceTag_);
    }
    else {
      assemble_->SetGeneralParams(pdename_, dofspernode_, 
                                  subdoms_, surfdoms_, bcSequenceTag_);
    }

    if (isComplex_) {
      assemble_->SetMatrixEntryType( OLAS::COMPLEX );
    }
    else {
      assemble_->SetMatrixEntryType(OLAS::DOUBLE);
      assemble_->SetMatrixStorageType(OLAS::SPARSE_NONSYM);
    }

    assemble_->SetPtr2Sol(sol_);

    // =====================================================================
    // read in material data
    // =====================================================================
    ReadMaterialData();

    // =====================================================================
    // define the integrators for PDE
    // =====================================================================
    DefineIntegrators();

    // =====================================================================
    // define which solution types have to be saved
    // =====================================================================
    ReadStoreResults();

    // check, if any output is calculated at all
    if ( hasOutput_ == FALSE ) {
      (*warning) << "There was no output specified at all for PDE '"
                 << pdename_
                 << "'Please check your .xml-file, if this is really what "
                 << "you want!";
      Warning( __FILE__, __LINE__ );
    }


    // =====================================================================
    // Create time stepping algorithm
    // =====================================================================
     if ( analysistype_ == TRANSIENT && 
          isDirectCoupled_ == FALSE) {
       SETPROFILE("Before Definition of Timestepping");
      InitTimeStepping();
      SETPROFILE("After Definition of TimeStepping");
    }

    PreparePDE4Computation();

    //! Define step solution driver
    if ( isDirectCoupled_ == FALSE ) {
      DefineSolveStep();
    }

  }

  
  void SinglePDE::SaveSolution( const Double * ptSol, UInt size ) {
    ENTER_FCN( "SinglePDE::SaveSolutionPointer" );

    Vector<Double> & solHelp = dynamic_cast<Vector<Double>&>(*solVec_);

    solHelp.Resize(size);

    for ( UInt i = 0; i < size; i++ ) {
      solHelp[i] = ptSol[i];
    }

    sol_->SetAlgSysDataPointer( size, solHelp.GetPointer() );
            
  }

  void SinglePDE::SaveSolution( const Complex * ptSol, UInt size ) {
    ENTER_FCN( "SinglePDE::SaveSolutionPointer" );

    Vector<Complex> & solHelp = dynamic_cast<Vector<Complex>&>(*solVec_);

    solHelp.Resize(size);

    for ( UInt i = 0; i < size; i++ ) {
      solHelp[i] = ptSol[i];
    }

    sol_->SetAlgSysDataPointer( size, solHelp.GetPointer() );
            
  }

  
  void SinglePDE::WriteGeneralPDEdefines() {
    
    ENTER_FCN( "SinglePDE::WriteGeneralPDEdefines" );
    
    //BC-section
    for (UInt i=0; i< bcs_hd_.GetSize(); i++) {
      UInt dof;
      std::string doftype = bcs_hd_[i];
      if (dofspernode_ > 1) 
        dof = domain->GetCoordSystem()->GetVecComponent(homDirichDof_[i]);
      else
        dof = 1;
      
      Info->WriteHomBC(pdename_, bcs_hd_[i], dof);      
    }
    
    for (UInt i=0; i< bcs_id_.GetSize(); i++) {
      UInt dof;
      std::string doftype = bcs_id_[i];
      if (dofspernode_ > 1) {
        dof = domain->GetCoordSystem()->GetVecComponent(inhomDirichDof_[i]);
      }
      else {
        dof = 1;
      }

      Info->WriteInHomBC( pdename_, bcs_id_[i], val_id_[i], fncnames_id_[i],
                          dof );        
    }

    // Loads
    StdVector<std::string> loadDom = GetLoadDom();
    StdVector<std::string> loadDof = GetLoadDof();
    StdVector<Double> loadVals = GetLoadVals();
    StdVector<std::string> loadfncs = GetLoadFncs();

    for( UInt i = 0; i < loadDom.GetSize(); i++ ) {
      UInt dof;
      std::string doftype = loadDom[i];
      if (dofspernode_ > 1) {
        dof = domain->GetCoordSystem()->GetVecComponent(loadDof[i]);
      }
      else {
        dof = 1;
      }
      Info->WriteLoad(pdename_, loadDom[i], loadVals[i], loadfncs[i], dof);
    }
  }



  // **********
  //   SetBCs
  // **********
  void SinglePDE::SetBCs( const Double time ) {
    
    ENTER_FCN( "SinglePDE::SetBCs" );
    
    UInt dof;
    Double val, val_tfunc;
    
    StdVector<UInt> nodes;

    UInt bcNum = 0;
    Integer eqnNr; 
    UInt eqnDof;

    if ( isIterCoupled_ ) {
      bcNum = couplingBCsCounter_;
    }
    else {
      bcNum = 0;
    }


    // ---------------------------
    // HOMOGENEOUS DIRICHLET BC
    // ---------------------------

    // Notes:
    //
    // 1.) We only need the following part, if we perform dof-blocking
    //
    // 2.) Dof-blocking is currently only supported for real-valued problems

    // What type of equation numbering does the user want?
    std::string typeOfNumbering;
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    keyVec = "linearSystems", "system", "setup", "eqnNumbering";
    attrVec = "", "name", "";
    valVec  = "", pdename_, "";

    params->Get( keyVec, attrVec, valVec, typeOfNumbering );

    if ( typeOfNumbering == "block" ) {

      val = 0;
      for ( UInt i = 0; i < bcs_hd_.GetSize(); i++ )   {  
        dof = 1;
        if ( dofspernode_ == 1 ) {
          (*error) << "dofspernode_ = 1, but I assumed dof-blocking???";
          Error( __FILE__, __LINE__ );
        }
        std::string doftype = bcs_hd_[i];
        dof = domain->GetCoordSystem()->GetVecComponent(homDirichDof_[i]);
        ptgrid_->GetNodesByName( nodes, bcs_hd_[i] );
      
        for ( UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {
          eqnData_->Node2EQN( nodes[iNode], dof, eqnNr, eqnDof );
          if (eqnNr > 0) {
            if ( analysistype_ == HARMONIC ||
                 analysistype_ == MULTIHARMONIC) {
              (*error) << "Dof-blocking for complex systems is not "
                       << "supported, yet!";
              Error( __FILE__, __LINE__ );
            }
            else {
              algsys_->SetDirichlet( bcNum + 1, pdeId_, eqnNr, val, eqnDof );
            }
            bcNum++;
          }
        }
      }
    }

    // ---------------------------
    // INHOMOGENEOUS DIRICHLET BC
    // ---------------------------
    
    Double phase = 0.0;
    Double dirVal;

    for ( UInt i = 0; i < bcs_id_.GetSize(); i++ ) {
      dof = 1;
      if ( dofspernode_ > 1 ) {
        std::string doftype = bcs_id_[i]; 
        dof = domain->GetCoordSystem()->GetVecComponent(inhomDirichDof_[i]);
      }
      
      ptgrid_->GetNodesByName( nodes, bcs_id_[i] ); 
      
      // Get the correct time function value
      val_tfunc = 1.0;
      if ( ptTimeFunc_->GetmaxTimeFnc() > 0 &&
           (analysistype_ != HARMONIC || analysistype_ != MULTIHARMONIC) ) {
        val_tfunc=ptTimeFunc_->TimeFuncAtTime(time,fncnames_id_[i]);
      }

      val    =  val_id_[i] * val_tfunc;
      dirVal = val;

      for ( UInt iNode = 0; iNode < nodes.GetSize(); iNode++ ) {

        eqnData_->Node2EQN(nodes[iNode], dof, eqnNr, eqnDof);

	// DODO GridAdaption
	// interpolate acoustic pressure at current node?
	if(m_bReadSpecialBCs) {
	  Point<3> pt;
	  UInt ibcNode = nodes[iNode];
	  m_pGridAdaption->GetCoordinatesFromGridNode(ptgrid_, ibcNode, pt);
	  
	  Double x, y, z;
	  x = pt[0];
	  y = pt[1];
	  z = pt[2];
	  
	  // perform interpolation with previously set interpolation scheme
	  val = m_pGridAdaption->GetAt(x, y, z, time);
	  if(effectiveMass_) {
	    dirVal = val;
	  }
	}

        // Sanity check. This should not happen, but might appear
        // in the case that the same node/dof belongs to a region
        // with hom. and a region with inhom. Dirichlet BCs. This
        // problem was already encountered!
        if (eqnNr == 0) {

          std::cout << "Node | dof | eqnNr | eqnDof\n"
                    << nodes[iNode]  << " | "
                    << dof   << " | "
                    << eqnNr << " | "
                    << eqnDof << std::endl;

          (*error) << "Got eqn number 0 for inhom Dirichlet BC! "
                   << "Probably you have a node/dof that belongs to both "
                   << "a region with hom. and one with inhom. Dirichlet BCs."
                   << " Go check your .mesh file!";
          Error( __FILE__, __LINE__ );
        }

        // Transform Dirichlet boundary conditions for effmass-formulation
        if (effectiveMass_) {
          val = dirVal;
          val = TS_alg_->DirichletBC4EffMassMatrix(val,eqnNr);
        }

        // Case of complex-valued entries
        if (analysistype_ == HARMONIC || analysistype_ == MULTIHARMONIC) {

          phase = bcs_id_phase_[i];
          //Complex complexValue( val * std::cos( phase / 180 * PI ),
          //                      val * std::sin( phase / 180 * PI ) );
          Complex complexValue( val * cos( phase / 180 * PI ),
                                val * sin( phase / 180 * PI ) );
          algsys_->SetDirichlet( bcNum + 1, pdeId_, eqnNr, complexValue,
                                 eqnDof );
        }
        else {
          algsys_->SetDirichlet( bcNum + 1, pdeId_, eqnNr, val, eqnDof );
        }
        bcNum++;
      }
    }
  }

 
  void SinglePDE::SetTimeStepping(TimeStepping *timeStepping) {
    ENTER_FCN( "SinglePDE::SetTimeStepping" );
    TS_alg_ = timeStepping;
  }


  void SinglePDE::ReadBCs() {

    ENTER_FCN( "SinglePDE::ReadBCs" );

    // vectors for parameter handling
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;
 
    // Get names of node sets for homogeneous Dirichlet boundary conditions
    //   params->GetList( "name", bcs_hd_, pdename_, "dirichletHom"   );
    keyVec = pdename_, "bcsAndLoads", "dirichletHom", "name";
    attrVec = "", "tag", "";
    valVec = "", bcSequenceTag_, "";
    params->GetList(keyVec, attrVec, valVec, bcs_hd_);

    // Get names of node sets, values and filenames for inhomogenous
    // Dirichlet boundary conditions
    keyVec = pdename_, "bcsAndLoads", "dirichletInhom", "name";
    params->GetList(keyVec, attrVec, valVec, bcs_id_);

    keyVec = pdename_, "bcsAndLoads", "dirichletInhom", "value";
    params->GetList(keyVec, attrVec, valVec, val_id_);

    if (dofspernode_ > 1) {
      keyVec = pdename_, "bcsAndLoads", "dirichletHom", "dof";
      params->GetList(keyVec, attrVec, valVec, homDirichDof_);
      
      keyVec = pdename_, "bcsAndLoads", "dirichletInhom", "dof";
      params->GetList(keyVec, attrVec, valVec, inhomDirichDof_);
    }
    
    if (analysistype_ == TRANSIENT ||
        analysistype_ == STATIC) {
      keyVec = pdename_, "bcsAndLoads", "dirichletInhom", "dynamics";
      params->GetList(keyVec, attrVec, valVec, fncnames_id_);

      // DODO
      // determine filename for reading boundary conditions calculated before
      std::string strFileName = "none"; // default
      keyVec  = pdename_, "bcsAndLoads", "dirichletInhom", "fileName";
      attrVec = "", "", "";
      valVec  = "", "", "";
      StdVector<std::string> vecFiles;
      params->GetList(keyVec, attrVec, valVec, vecFiles);
      //params->Get(keyVec, strFileName);

      if(vecFiles.GetSize() > 0) {
	strFileName = vecFiles[0];
      }

      // if attribute "fileName" is set, we're in "post" mode, i.e. set flag
      // to use interpolation techqnique
      if(strFileName != "none") {
	m_bReadSpecialBCs = TRUE;

	// create new adaption object
	m_pGridAdaption = new GridAdaption();

	// determine no of nodes to use for interpolation
	keyVec = pdename_, "bcsAndLoads", "dirichletInhom", "noNodes";
	std::string strNoNodes;
	params->Get(keyVec, strNoNodes);

	// how many? atoi() strips percentage sign
	Integer nNodes = atoi(strNoNodes.c_str());

	// check if to use percentage, atoi() takes care of conversion
	bool bPercent = false;
	if(strNoNodes[strNoNodes.length()-1] == '%') {
	  bPercent = true;
	}
	m_pGridAdaption->SetInfluencingNodes(nNodes, bPercent);

	// default interpolation type is nearest neighbor
	GridAdaption::INTERPOLATION_TYPE it = GridAdaption::NNB;

	// now check for interpolation type, if given
	StdVector<std::string> vecInterpolType;
	std::string strInterpolType;
	keyVec = pdename_, "bcsAndLoads", "dirichletInhom", "interpolType";
	params->Get(keyVec, strInterpolType);

	// check interpolation type, if other than NNB
	if(strInterpolType == "IDW") {
	  it = GridAdaption::IDW;
	}
	else if(strInterpolType == "SHP") {
	  it = GridAdaption::SHP;
	}

	// now set type, default is NNB
	m_pGridAdaption->SetInterpolationType(it);

	// now parse given data file
	m_pGridAdaption->ReadFile(strFileName);
      }
      // DODO


    }
    else if (analysistype_ == HARMONIC||analysistype_ == MULTIHARMONIC) {
      keyVec = pdename_, "bcsAndLoads", "dirichletInhom", "phase";
      params->GetList(keyVec, attrVec, valVec, bcs_id_phase_);
    }

    // Check consistency
    if ( bcs_id_.GetSize() != val_id_.GetSize() ||
         fncnames_id_.GetSize() > bcs_id_.GetSize() ) {
      std::string errmsg = "dirichletInhom: ";
      errmsg += "#names of node sets = " + Info->GenStr(bcs_id_.GetSize());
      errmsg += ", #values = " + Info->GenStr(val_id_.GetSize());
      errmsg += ", #dynamics = " + fncnames_id_.GetSize() + '\n';
      Info->Error( errmsg, __FILE__, __LINE__ );
    }

    // We need not have as many function/filenames as boundary conditions!
    for ( UInt k = fncnames_id_.GetSize(); k < bcs_id_.GetSize(); k++ ) {
      fncnames_id_.Push_back( "none" );
    }
    
    // =====================================================================
    // if pde has more than one dof, initialize dof of boundary conditions
    // =====================================================================
    if ( dofspernode_ > 1 ) {
      if (bcs_hd_.GetSize() != homDirichDof_.GetSize()) {

        std::string errmsg = "Inconsistent definition of homogeneous ";
        errmsg += "Dirichlet Boundary Conditions\n";
        errmsg += " bcs_hd_.GetSize() = " + Info->GenStr( bcs_hd_.GetSize() );
        errmsg += "\n homDirichDof_.GetSize() = "
          + Info->GenStr( homDirichDof_.GetSize() ) + '\n';
        Info->Error( errmsg, __FILE__, __LINE__ );
      }
      if (bcs_id_.GetSize() != inhomDirichDof_.GetSize()) {
        std::string errmsg = "Inconsistent definition of inhomogeneous ";
        errmsg += "Dirichlet Boundary Conditions";
        Info->Error( errmsg, __FILE__, __LINE__ );
      }
    }

    // Calculate correct number of boundary conditions
    for ( UInt i = 0; i < bcs_hd_.GetSize(); i++ ) {
      numDirichletBCs_ += ptgrid_->GetNumNodes(bcs_hd_[i]);
    }
    
    for ( UInt i = 0; i < bcs_id_.GetSize(); i++) {
      numDirichletBCs_ += ptgrid_->GetNumNodes(bcs_id_[i]);
    }
  }


  void SinglePDE::ReadMaterialData() {
    ENTER_FCN( "SinglePDE::ReadMaterialData" );

    // -------------------
    // XMLPARAMS
    // -------------------

    std::string outformat="unverg";
    std::string matFileName;

    // Allocate space to hold material data for each subdomain of this PDE
    materialData_ = new MaterialData[subdoms_.GetSize()];

    // Get list of subdomains and materials
    StdVector< std::string > subdomName;
    StdVector< std::string > subdomMaterial;
    StdVector< RegionIdType > subdomId;
    params->GetList( "name", subdomName, "domain", "region" );
    params->GetList( "material", subdomMaterial, "domain", "region" );
    
    // Convert region names to Ids
    ptgrid_->RegionNameToId( subdomId, subdomName );

    params->Get("format", outformat, "output");

    if ( outformat != "database" ) {

      // Query name of file with material data
      params->Get( "file", matFileName, "materialData" );
    
      // Generate new material reader
      LoadMaterialDataFile loadMaterialFile( matFileName.c_str() );
      
      // Load material data for subdomains on which this PDE lives
      // from data file
      for( UInt i = 0; i < subdoms_.GetSize(); i++ ) {
        for( UInt k = 0; k <= subdomId.GetSize(); k++ ) {
          if( subdoms_[i] == subdomId[k] ){
            loadMaterialFile.GetMaterial( materialData_[i], subdomMaterial[k],
                                          pdematerialclass_ );
            break;
          }
        }
      }
    }
    else {

#ifdef USE_DATABASE
      LoadMaterialDataDatabase loadMaterialDB;
    
      // Load material data for subdomains on which this PDE lives
      // from data file
      for( UInt i = 0; i < subdoms_.GetSize(); i++ ) {
        for( UInt k = 0; k <= subdomName.GetSize(); k++ ) {
          if( subdoms_[i] == subdomName[k] ) {
            loadMaterialDB.GetMaterial( materialData_[i], subdomMaterial[k],
                                        pdematerialclass_ );
            break;
          }
        }
      }
#else  // No Database
      (*error) << "Re-compile with DATABASE = yes to get database support";
      Error( __FILE__, __LINE__ );
#endif
    }
  }


  UInt SinglePDE::GetNumRestraints() {
    
    ENTER_FCN( "SinglePDE::GetNumRestraints" );
    
    return numDirichletBCs_;
  }
  
  // ======================================================
  // GET /SET  METHODS
  // ======================================================

  //! Activate the direct coupling
  void SinglePDE::SetDirectCoupling () {
                           
    ENTER_FCN( "SinglePDE::SetDirectCoupling" );
    
    isDirectCoupled_ = TRUE;          
  }

  void SinglePDE::SetAlgebraicSystem( BaseSystem *algSys) {
    ENTER_FCN( "SinglePDE::SetAlgebraicSystem" );
    algsys_ = algSys;
  }

  void SinglePDE::SetSolveStep ( StdSolveStep *solveStep) {
    ENTER_FCN( "SinglePDE::SetSolveStep" );
    solveStep_ = solveStep;
  }


  // ======================================================
  // ALGSYS SECTION (SOLVER, ...) 
  // ======================================================
  void SinglePDE::DefineAlgSys() {

    ENTER_FCN( "StdPDE::DefineAlgSys" );


    // If PDE is not direct coupled then the PDE has to register
    // at the algebraic system and obtain an Id. 
    // Afterwards the matrix-graph has to be set up
    SETPROFILE("Before GraphSetupInit()");
    if ( isDirectCoupled_ == FALSE ) {

      // Set linear system parameters for OLAS
      ReadOlasParams( pdename_ );

      // Initialize the matrix graph object
      algsys_->GraphSetupInit(1);

      // obtain PDE identification tag from algebraic system
      algsys_->RegisterPDE( pdeId_, eqnData_->GetNumEQNs(),
                            eqnData_->GetNumLastFreeDof() );

      assemble_->SetPDEId( pdeId_ );

      solveStep_->SetPDEId( pdeId_ );
      
      
      // trigger the creation and assembly of the matrix graph
      assemble_->SetupMatrixGraph();
    
      // finish the assembly of the matrix graph
      algsys_->GraphSetupDone();
      
      // obtain reordering of the matrix graph
      // and give it to the EQN-object
      Integer * newOrder = algsys_->GetReordering( pdeId_ );
      eqnData_->ReorderMapping( newOrder );

    }

    // pass information about dofs, number of dirichlet equations
    // and constraints to the algebraic system
    algsys_->SetBlockSize( pdeId_, eqnData_->GetNumDofsPerEQN() );

    UInt numDir = GetNumRestraints() - eqnData_->GetNumDroppedDofs();
    algsys_->SetNumDirichletBCs(pdeId_, numDir );

    // create matrices and solver object, if PDE is not direct coupled
    if ( isDirectCoupled_ == FALSE ) {
      CreateMatrices_Solver();
    }

  }

  // ======================================================
  // METHODS FOR ASSEMBLING
  // ======================================================
  
  void  SinglePDE::AssembleMatrices() {
    assemble_->AssembleMatrices();
  }
  
  void  SinglePDE::AssembleSrcRHS(const Double time) {
    assemble_->AssembleSrcRHS(time);
  }
  
  void SinglePDE::AssembleNLRHS(const Double time) {
    assemble_->AssembleNLRHS(time);
  }
  
  void  SinglePDE::AssembleSprings(const Double time) {
    assemble_->AssembleSprings(time);
  }
  
  void  SinglePDE::InitNonLinMatrices() {
    assemble_->InitNonLinMatrices();
  }
  
  // constructs the matrix graph by providing to the algebraic system
  // the element connectivities
  void  SinglePDE::SetupMatrixGraph() {
    assemble_->SetupMatrixGraph();
  }

  void SinglePDE::SetFrequency(Double actFreq) {
    assemble_->SetFrequency(actFreq);
  }

  void SinglePDE::SetReassemble() {
    assemble_->SetReassemble();
  }

  
  // ======================================================
  // Adaptivity SECTION 
  // ======================================================

#ifdef ADAPTGRID
  void SinglePDE::RefineMesh( const Integer level) {

    ENTER_FCN( "SinglePDE::RefineMesh" );
  
    Integer          numChilds;
    Integer          numElems;
    StdVector<Elem*> elemssd;
    Double           theta_s;
    Double           coarseFactor;
    Double           tol4Elm;
    Integer          numRefLoops=0;
    Integer          numRefinements;
    Integer          iem;

    // get max num elements for the domain,on which we have the equation
    numElems=ptgrid_->GetMaxnumElem(level,subdoms_);

    // get pointer to array with elements
    ptgrid_->GetElemSD(elemssd,subdoms_[0],level);

    // if element is marked, then value of the array's element is equal 1,
    // else 0.
    markingElems_.Resize(numElems);

    if (!conf->ifget("safety_factor_for_space_adaptivity",theta_s)) {
      theta_s=1.0;
    }

    coarseFactor = 0.7;
    tol4Elm = theta_s*tolSpaceErr_;

    std::cout << " tolerance space error: " << tolSpaceErr_ << std::endl;
 
    // loop over elements
    for (iem=0; iem<numElems; iem++) {
      elemssd[iem]->refinementFlag = 0;
      if (errorMap_[iem]>tol4Elm) {
        elemssd[iem]->refinementFlag = 1; 

        numChilds=elemssd[iem]->ptElem->getNumChilds();

        numRefinements = defineRefinements(errorMap_[iem],tol4Elm,numChilds);

        if (numRefinements>numRefLoops)
          numRefLoops = numRefinements;

        markingElems_[iem]=numRefinements;
        elemssd[iem]->refinementNumber = numRefinements;
      }
      else {
        if (errorMap_[iem] < coarseFactor*tol4Elm) {
          elemssd[iem]->refinementFlag = 0;
        }
        else {
          elemssd[iem]->refinementFlag = 0;
        }
      }
    }
    
    // Fuehren die Verfeinerung durch
    std::cout << " number of refinement loops \n" << numRefLoops << std::endl;

    ptgrid_->Refine(numRefLoops);

  }


  Boolean SinglePDE::TestError(const Integer level) {

    ENTER_FCN( "SinglePDE::TestError" );

    if (!ptError_) ConstructorError();
  
    // Berechnung der Fehlerkarte
    Double            totalErr;
    Vector<Double>    solVec;
  
    solVec.Resize(sol_.size());
    int i;
    for (i=0; i<sol_.size(); i++)
      solVec[i]=sol_[0][i];
    //      sol_.toVector(solVec,1);
  
    ptError_->CalcErrorMap(solVec,subdoms_,ptgrid_,errorMap_,totalErr,level);
  
    std::cout << " space error: " << totalErr <<
      " tolerance: " << tolSpaceErr_ << std::endl;

    if (totalErr > tolSpaceErr_) return TRUE;
    else return FALSE;
  }


  //In this fnc we delete old pointer to Error-object & create new
  void SinglePDE::ConstructorError() {

    ENTER_FCN( "SinglePDE::ConstructorError" );

    if (ptError_) delete ptError_;
  
    ptError_=new SpaceErrorEstimator();
    ptError_->Init(this);
  }


  void SinglePDE::WriteErrorInfo(WriteResults * ptmeshes) {
    ptmeshes->WriteElemSolution(errorMap_,0,0,"ERR-errorMap");
    ptmeshes->WriteElemSolution(markingElems_,0,0,"ERR-markedElems");
  }

#endif // AdaptGrid

  // ======================================================
  // COUPLING SECTION
  // ======================================================
  
  void SinglePDE::CalcInputCoupling() {
    
    ENTER_FCN( "SinglePDE::CalcInputCoupling" );

    std::string errMsg;
    StdVector<UInt> * nodes;
    CFSVector * val;
    UInt pdeNode, eqnDof;
    Integer eqnNr;
    UInt couplingDof;
    Boolean clearCoords = TRUE;

    

    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == FALSE)
      return;

    // Reset counter for boundary conditions
    couplingBCsCounter_ = 0;
 

    // Outer loop over all INPUT coupling terms
    for (UInt i=0; i<ptCoupling_->GetNumInputCouplings(); i++) {

      //    ptCoupling_ = &ptCoupling_[i];
      ptCoupling_->GetInputValues(i, val);
      couplingDof = ptCoupling_->GetInputDof(i);
    
      // Up to now, Coupling is only possible with
      // Real valued solutions
      Vector<Double> const & help = dynamic_cast<Vector<Double>&>(*val);

      switch(ptCoupling_->GetInputType(i)) {
        
        // -------------------
        // COORDINATE COUPLING
        // -------------------
      case COORD: 
        
        // Set flag that the geometry has changed
        geoUpdate_ = TRUE;
        assemble_->SetNonlinGeo();

        ptCoupling_->GetInputNodes(i, nodes);

        // Resize + clear coordinate updates
        // only the first time
        if (clearCoords == TRUE) {
          deltCoords_.Resize(dim_, numPDENodes_);
          clearCoords = FALSE;
        }
          
        // set ptr of deltCoords to assembly-object
        assemble_->SetPtrDeltaCoordinates(&deltCoords_);
          
        for (UInt j=0; j<nodes->GetSize(); j++)
          for (UInt dof=0; dof<ptCoupling_->GetInputDof(i); dof++) {
            pdeNode = eqnData_->Mesh2PDENode((*nodes)[j]);
            deltCoords_(dof,pdeNode-1) = help[dof + j*dim_];

          }
        break;

        // -------------------
        // RHS COUPLING
        // -------------------
      case RHS:
        ptCoupling_->GetInputNodes(i, nodes);
          
        //for (UInt dof=0; dof<ptCoupling_->GetInputDof(i); dof++)
        for ( UInt dof = 0; dof < couplingDof; dof++ ) {
          for ( UInt j = 0; j < nodes->GetSize(); j++ ) {
            eqnData_->Node2EQN((*nodes)[j],dof+1,eqnNr,eqnDof);
            if ( eqnNr != 0 ) {
              algsys_->SetNodeRHS( help[ dof + couplingDof * j ], pdeId_,
                                   eqnNr, eqnDof );
            }
          }
        }
        break;

        // -----------------------
        // InhomDirichlet COUPLING
        // -----------------------
      case ID_BC:

        // How do we want to treat inhomogeneous Dirichlet
        // boundary conditions?
        {
          bool usePenalty = true;
          std::string aux;
          StdVector<std::string> keyVec;
          StdVector<std::string> attrVec;
          StdVector<std::string> valVec;
          keyVec  = "linearSystems", "system", "setup", "idbcHandling";
          attrVec = "", "name", "";
          valVec  = "", pdename_, "";
          params->Get( keyVec, attrVec, valVec, aux );
          usePenalty = aux == "penalty" ? true : false;
          Info->PrintF( pdename_, "Treating IDBCs using '%s' approach\n",
                        aux.c_str() );
          if ( usePenalty == false ) {
            (*error) << "Cannot use inhom. Dirichlet coupling together with "
                     << "IDBC elimination, since the equation numbering "
                     << "object does currently not have the information "
                     << "required to put those values at the end of the "
                     << "equation number interval! Please set idbcHandling="
                     << '"' << "penalty" << '"' << " in your xml-file";
            Error( __FILE__, __LINE__ );
          }
        }

        // Set flag that the boundary conditions have to be incorporated
        updateCouplingBCs_ = TRUE;

        ptCoupling_->GetInputNodes(i, nodes);

        for ( UInt dof = 0; dof < ptCoupling_->GetInputDof(i); dof++ ) {
          for ( UInt j = 0; j < nodes->GetSize();
                j++, couplingBCsCounter_++) {

            eqnData_->Node2EQN((*nodes)[j],dof+1,eqnNr,eqnDof);

            if (eqnNr==0) {
              Error( "The specified coupling node has no equation number",
                     __FILE__, __LINE__ );
            }

            algsys_->SetDirichlet( couplingBCsCounter_ + 1, pdeId_, eqnNr,
                                   help[dof+j*couplingDof], eqnDof );
          }
        }
        break;
          
      case MAT:
        Error( "Not implemented yet", __FILE__, __LINE__ );
        break;

      }
    }
  }
  
  void SinglePDE::ResetCoupling() {
    ENTER_FCN( "SinglePDE::ResetCoupling" );
    
    iterCoupledCounter_ = 0;
    
    
  }

} // end of namespace
