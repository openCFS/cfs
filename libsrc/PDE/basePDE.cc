#include <fstream>
#include <iostream>
#include <string>

#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/ParamHandling/CFSOLASParams.hh"
#include "DataInOut/WriteInfo.hh"
#include "Domain/elem.hh"
#include "Estimator/spaceerror.hh"
#include "basePDE.hh"
#include "scalarblockEQN.hh"
#include "scalarnodeEQN.hh"
#include "blocknodeEQN.hh"
#include "superblockEQN.hh"
#include "newmarkFracDamp.hh"


namespace CoupledField {

  BasePDE::BasePDE( Grid *aptgrid, BCs *aptBCs, FileType *aInFile,
                    WriteResults * aOutFile, TimeFunc *aptTimeFunc )
    : nonLin_(FALSE),
      incStopCrit_(1e-2), 
      residualStopCrit_(1e-3),
      firstTimeStepStatic_(TRUE),
      isaxi_(FALSE),
      isComplex_(FALSE),
      numPDENodes_(0),
      numElems_(0),
      sol_(NULL),
      isAlwaysStatic_(FALSE),
      dampingType_(NONE),
      needsDampingMatrix_(FALSE),
      isIncrFormulation_(FALSE),
      TS_alg_(NULL),
      materialData_(NULL),
      effectiveMass_(FALSE)
  {

    ENTER_FCN( "BasePDE::BasePDE" );

    nonLin_ = FALSE;
    incStopCrit_ = 1e-2;
    residualStopCrit_ = 1e-3;


    // =====================================================================
    // set file pointers
    // =====================================================================
    inFile_     = aInFile;
    outFile_    = aOutFile;
    ptTimeFunc_ = aptTimeFunc;
    ptgrid_     = aptgrid;
    ptBCs_      = aptBCs;
    assemble_   = NULL;
    algsys_     = NULL;
    eqnData_    = NULL;

    // =====================================================================
    // set analysis parameters
    // =====================================================================
    actlevel_ = 0;
    actFrequency_ = 0;
    complexFormat_ = AMPLITUDE_PHASE; // or REAL_IMAG
    couplingBCsCounter_ = 0;
    numDirichletBCs_ = 0;
    pdeIsCoupled_ = FALSE;
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
    numeqcoarse_ = 200;
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
  }

  
  void BasePDE::Init( Integer bcSequenceIndex,
                      std::string  bcSequenceTag ) {
    ENTER_FCN( "BasePDE::Init()" );

    bcSequenceIndex_ = bcSequenceIndex;
    bcSequenceTag_ = bcSequenceTag;

    // =====================================================================
    // get regions/subdomains for PDE
    // =====================================================================
    params->GetList( "name", subdoms_, pdename_, "region" );
    Info->PrintF( pdename_, "The %s PDE lives on the following regions:\n",
		  pdename_.c_str());
    for ( Integer k = 0; k < subdoms_.GetSize(); k++ ) {
      Info->PrintF( pdename_, " %s, index %d\n", subdoms_[k].c_str(), k );
    }
    Info->PrintF( "", "\n" );

    // Generate a fitting algebraic system
    algsys_ = new StandardSystem();

    // Get parameter and report object of OLAS
    olasParams_ = algsys_->GetOLASParams();
    olasReport_ = algsys_->GetOLASReport();

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
    
    if (analysisHelp == STATIC ||
        isAlwaysStatic_ == TRUE) {
      isComplex_ = FALSE;
      assemble_ = new StaticAssemble(algsys_, ptgrid_);
      analysistype_ = STATIC;
    }

    else if (analysisHelp == TRANSIENT) {
      isComplex_ = FALSE;
      assemble_ = new TransientAssemble(algsys_, ptgrid_);
      analysistype_ = TRANSIENT;
      laststepcalc_ = 1;
    }

    else if ( analysis=="harmonic" || analysis == "paramIdent" ) {
      isComplex_ = TRUE;
      assemble_ = new HarmonicAssemble(algsys_, ptgrid_);
      analysistype_ = HARMONIC;
    }

    else if ( analysisHelp == MULTI_SEQUENCE ) {

      stepString = Info->GenStr(bcSequenceIndex_);
      attrVec = "", "index", "type";
      valVec = "", stepString, pdename_;

      keyVec = "multiSequence", "step", "pde", "analysis";
      params->Get(keyVec, attrVec, valVec, analysis);

      String2Enum(analysis, analysistype_);

      if ( analysistype_ == STATIC ) {
        assemble_ = new StaticAssemble(algsys_, ptgrid_);
        isComplex_ = FALSE;
      }
      else if ( analysistype_ == TRANSIENT ) {
        isComplex_ = FALSE;
        assemble_ = new TransientAssemble(algsys_, ptgrid_);
        laststepcalc_ = 1;      
      }
      else {
        Error("BasePDE::Init: AnalysisType not supported", __FILE__, __LINE__);
      }
    }
    else {
      Error("Analysis Type not supported",__FILE__,__LINE__);
    }
    
    // Determine if solution is of complex type or not
    if ( analysistype_ == HARMONIC ) {
      sol_ = new NodeStoreSol<Complex>;
    }
    else {
      sol_ = new NodeStoreSol<Double>;
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
    // read in boundary conditions
    // =====================================================================
    ReadBCs();
    ReadSpecialBCs();
    numDirichletBCs_ += GetNumRestraints(actlevel_);

    // =====================================================================
    // read in NonLinearities
    // =====================================================================
    InitNonLin();

    // =====================================================================
    // initialize EQN-object and Storeresults class
    // =====================================================================

    // What type of equation numbering does the user want?
    keyVec = pdename_, "solver", "matrix", "eqnNumbering";
    std::string typeOfNumbering;
    params->Get( keyVec, typeOfNumbering );

    // Assemble a system matrix with scalar complex or double entries
    if ( typeOfNumbering == "scalar" ) {
      if ( dofspernode_ == 1 ) {
        eqnData_ = new ScalarNodeEQN( ptgrid_, ptBCs_, subdoms_, actlevel_,
                                      dofspernode_ );
      }
      else {
        eqnData_ = new ScalarBlockEQN( ptgrid_, ptBCs_, subdoms_, actlevel_,
                                       dofspernode_ );
      }
      Info->PrintF( pdename_, "Using scalar equation numbering\n" );
    }

    // Treat all dofs of a node together and assemble a system matrix with
    // small square matrices as entries
    else if ( typeOfNumbering == "block" ) {
      if ( dofspernode_ == 1 ) {
        Warning("dopspernode = 1, so 'block' numbering identical to 'scalar'");
        eqnData_ = new ScalarNodeEQN( ptgrid_, ptBCs_, subdoms_, actlevel_,
                                      dofspernode_ );
      }
      else {
        eqnData_ = new BlockNodeEQN( ptgrid_, ptBCs_, subdoms_, actlevel_,
                                     dofspernode_ );
        Info->PrintF( pdename_, "Using block equation numbering\n" );
      }
    }

    // Build in Dirichlet boundary conditions and compute the mapping
    // which relates nodes/dofs to equation numbers
    eqnData_->SetHomoDirichletBCs ( bcs_hd_, homDirichDof_  );
    eqnData_->SetInhomDirichletBCs( bcs_id_, inhomDirichDof_);
    eqnData_->CalcMapping();

    // Report results to logfile
    Info->PrintF( pdename_, "Linear system will have %d equations\n\n",
                  eqnData_->GetNumEQNs() );

    //eqnData_->Print(*debug);
    numPDENodes_ = eqnData_->GetNumLocalNodes();
    numElems_ = eqnData_->GetNumLocalElems();
    
    // Initialize Storesolution class
    sol_->SetNumSolutions(solTypes_.GetSize());
    sol_->SetNumNodes(numPDENodes_);
    for (Integer iSol=0; iSol<solTypes_.GetSize(); iSol++) {
      sol_->SetSolutionType(solTypes_[iSol],iSol);
      sol_->SetNumDofs(solDofs_[iSol], solTypes_[iSol]);
    }
    sol_->SetPtrEQNData(eqnData_, ptgrid_, actlevel_);
    sol_->Init(); 


    // =====================================================================
    // initialize assemble object
    // =====================================================================
    assemble_->SetPtr2EQNData(eqnData_); 
    assemble_->SetPtr2TimeFnc(ptTimeFunc_);

    if (pdename_ == "piezo" || pdename_ == "mechanic" ) {
      assemble_->SetGeneralParams(pdename_, dofspernode_, numPDENodes_, 
                                  subdoms_, pressSurf_, bcSequenceTag_);
    }
    else {
      assemble_->SetGeneralParams(pdename_, dofspernode_, numPDENodes_, 
                                  subdoms_, surfdoms_, bcSequenceTag_);
    }

    assemble_->SetGraphType(NODEGRAPH);

    if (isComplex_) {
      assemble_->SetMatrixEntryType( OLAS::COMPLEX );
    }
    else {
      assemble_->SetMatrixEntryType(OLAS::DOUBLE);
      assemble_->SetMatrixStorageType(OLAS::SPARSE_NONSYM);
    }

    assemble_->SetNumDirichlet(numDirichletBCs_);

    assemble_->SetPtrBCs(ptBCs_);
    assemble_->SetPtr2Sol(sol_);
    if (needsDampingMatrix_) {
      assemble_->NeedDampingMatrix();
    }

    // =====================================================================
    // read in material data
    // =====================================================================
    ReadMaterialData();

    // =====================================================================
    // define the integratos for PDE
    // =====================================================================
    DefineIntegrators(actlevel_);

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
    InitTimeStepping();

    PreparePDE4Computation();

    //!
    DefineSolveStep();
  }


  Integer BasePDE::GetBCDof( const std::string dofStartString ) {

    ENTER_FCN( "BasePDEPDE::GetBCDof" );

    Integer nrActDof = 0;

    if ( dofStartString == "ux" )
      nrActDof = 1;
    if ( dofStartString == "uy" )
      nrActDof = 2;
    if ( dofStartString == "uz" )
      nrActDof = 3;
    // HARD coded for piezo case
    if ( dofStartString == "ep" )
      nrActDof = dofspernode_;
    if ( nrActDof == 0 )
      Error( "Unknown dof-type in homog. BC; substring must start with ux, uy,"
             "uz or ep!!", __FILE__, __LINE__);

    return nrActDof;
  }


  void BasePDE::WriteGeneralPDEdefines() {

    ENTER_FCN( "BasePDE::WriteGeneralPDEdefines" );

    //BC-section
    for (Integer i=0; i< bcs_hd_.GetSize(); i++) {
      Integer dof;
      std::string doftype = bcs_hd_[i];
      if (dofspernode_ > 1) 
        dof = GetBCDof(homDirichDof_[i]);
      else
        dof = 1;

      Info->WriteHomBC(pdename_, bcs_hd_[i], dof);      
    }

    for (Integer i=0; i< bcs_id_.GetSize(); i++) {
      Integer dof;
      std::string doftype = bcs_id_[i];
      if (dofspernode_ > 1) {
        dof = GetBCDof(inhomDirichDof_[i]);
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

    for( int i = 0; i < loadDom.GetSize(); i++ ) {
      Integer dof;
      std::string doftype = loadDom[i];
      if (dofspernode_ > 1) {
        dof = GetBCDof(loadDof[i]);
      }
      else {
        dof = 1;
      }
      Info->WriteLoad(pdename_, loadDom[i], loadVals[i], loadfncs[i], dof);
    }
  }


  // initialize PDEs before iteration (done for each time step)
  void BasePDE::InitStepTransCoupled( Double asteptime ) {
    ENTER_FCN( "BasePDE::InitStepTransCoupled" );
    lasttimecalc_ = asteptime;
    iterCoupledCounter_ = 0;
  }


  void BasePDE::SetBCs( const Integer level,  const Double time ) {
    ENTER_FCN( "BasePDE::SetBCs" );

    Integer node, dof;
    Double val, val_tfunc;

    std::list<Integer> nodes;

    Integer i;
    Integer j;
    Integer eqnNr, eqnDof;

    if (pdeIsCoupled_) {
      j = couplingBCsCounter_;
    }
    else {
      j=0;
    }


    // ---------------------------
    // HOMOGENEOUS DIRICHLET BC
    // ---------------------------
    val = 0;
    for ( i = 0; i < bcs_hd_.GetSize(); i++ )   {  
      dof = 1;
      if ( dofspernode_ > 1 ) {
        std::string doftype = bcs_hd_[i];
        dof = GetBCDof(homDirichDof_[i]);
      }

      nodes = ptBCs_->GetNodesLevel(bcs_hd_[i]);
      std::list<Integer>::const_iterator p;
      for ( p = nodes.begin(); p != nodes.end(); p++ ) {
        node = *p;
        eqnData_->Node2EQN(node, dof, eqnNr, eqnDof);
        if (eqnNr > 0) {

          // Increment counter for BCs
          if ( analysistype_ == HARMONIC ) {

            // set real part 
            algsys_->SetDirichlet(j*2+1, eqnNr, val, eqnDof);

            // set imaginary part 
            algsys_->SetDirichlet(j*2+2, eqnNr, val, eqnDof+1);
          }
          else {
            algsys_->SetDirichlet(j+1, eqnNr ,val, eqnDof);
          }
          j++;
        }
      }
    }


    // ---------------------------
    // INHOMOGENEOUS DIRICHLET BC
    // ---------------------------

    Double phase = 0.0;
    Double dirVal;
    for ( i = 0; i < bcs_id_.GetSize(); i++ ) {
      dof = 1;
      if ( dofspernode_ > 1 ) {
        std::string doftype = bcs_id_[i]; 
        dof = GetBCDof(inhomDirichDof_[i]);
      }

      nodes = ptBCs_->GetNodesLevel(bcs_id_[i]); 

      //get the correct time function value
      val_tfunc = 1.0;
      if (ptTimeFunc_->GetmaxTimeFnc() > 0 ) {
        val_tfunc=ptTimeFunc_->TimeFuncAtTime(time,fncnames_id_[i]);
      }
      
      val    =  val_id_[i] * val_tfunc;
      dirVal = val;
      std::list<Integer>::const_iterator p;
      for ( p = nodes.begin(); p != nodes.end(); p++, j++ ) {
        node = *p;
        eqnData_->Node2EQN(node, dof, eqnNr, eqnDof);

	//transform Dirichlet boundary conditions for effmass-formulation
	if (effectiveMass_) {
	  val = dirVal;
	  val = TS_alg_->DirichletBC4EffMassMatrix(val,eqnNr);
	}

        if (analysistype_ == HARMONIC) {
          phase = bcs_id_phase_[i];

          // set real part
          algsys_->SetDirichlet(j*2+1, eqnNr, val * cos(phase/180*PI), eqnDof);

          // set imaginary part 
          algsys_->SetDirichlet(j*2+2, eqnNr, val * sin(phase/180*PI),
                                eqnDof+1);
        }
        else {
          algsys_->SetDirichlet(j+1, eqnNr, val, eqnDof);
        }
      }
    }
  }

  void BasePDE::GetMemento(PDEMemento & memento) {

    ENTER_FCN( "BasePDE::GetMemento" );

    // first get memento of coupling object
    if (pdeIsCoupled_) {
      ptCoupling_->GetMemento(memento.couplingMemento_);
    }

    // then write own data to PDEMemento
    memento.analysisType_ = analysistype_;

    if (memento.sol_)
      delete memento.sol_;

    if ( analysistype_ == STATIC || analysistype_ == TRANSIENT ) {

      // --- Real values --
      // Set solution
      memento.sol_ = new Vector<Double>;
      dynamic_cast<Vector<Double>&>(*(memento.sol_)) =
        dynamic_cast<NodeStoreSol<Double>&>(*(sol_)).GetAlgSysVector();

      if (analysistype_ == TRANSIENT) {
        Warning( "Currently only the first derivative is stored in the "
                 "memento object!", __FILE__, __LINE__ );

        // Set first derivative
        memento.solDeriv1_ = getS1();   
        
        // Set secondderivative
        memento.solDeriv2_ = getS2();
      }
    }
    else {

      // --- Complex values --      
      // Set solution
      memento.sol_ = new Vector<Complex>;
      dynamic_cast<Vector<Complex>&>(*(memento.sol_)) =
        dynamic_cast<NodeStoreSol<Complex>&>(*(sol_)).GetAlgSysVector();
    } 

    // now memento is initialized
    memento.isSet_ = TRUE;
  }


  void BasePDE::SetMemento( PDEMemento &memento ) {

    ENTER_FCN( "BasePDE::SetMemento" );
  
    // if there is no information in the menmento just leave
    if ( memento.isSet_ == FALSE ) {
      return;
    }
  
    if ( analysistype_ == STATIC || analysistype_ == TRANSIENT ) {

      // --- Real values --
      // Set solution
      dynamic_cast<NodeStoreSol<Double>&>(*(sol_)).SetAlgSysVector
        (dynamic_cast<Vector<Double>&>(*(memento.sol_)));
      
      // if previous step was transient and the current step is also
      // then give the time derivative to the timestepping algorithm
      if (analysistype_ == TRANSIENT
          && memento.analysisType_ == TRANSIENT) {
        TS_alg_->SetDeriv1(memento.solDeriv1_);
        TS_alg_->SetDeriv2(memento.solDeriv2_);
      }
    }

    else {

      // --- Complex values --      
      // Set solution
      dynamic_cast<NodeStoreSol<Complex>&>(*(sol_)).SetAlgSysVector
        (dynamic_cast<Vector<Complex>&>(*(memento.sol_)));
    }

    // If PDE is coupled, set the according coupling memento
    if (pdeIsCoupled_) {
      ptCoupling_->SetMemento(memento.couplingMemento_);
    }
  }


  const Vector<Double>& BasePDE::getS1() const {

    ENTER_FCN( "BasePDE::getS1" );
  
    if ( TS_alg_ != NULL ) {
      return TS_alg_->GetDeriv1();
    }
    else {
      (*error) << pdename_ << ":getS1: No timestepping defined for this PDE";
      Error( __FILE__, __LINE__ );
    }
  }


  const Vector<Double>& BasePDE::getS2() const {

    ENTER_FCN( "BasePDE::getS2" );
  
    if ( TS_alg_ != NULL ) {
      return TS_alg_->GetDeriv2();
    }
    else {
      (*error) << pdename_ << ":getS2: No timestepping defined for this PDE";
      Error( __FILE__, __LINE__ );
    }
  }


  void BasePDE::ReadBCs() {

    ENTER_FCN( "BasePDE::ReadBCs" );

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

    if (analysistype_ == TRANSIENT ||
        analysistype_ == STATIC) {
      keyVec = pdename_, "bcsAndLoads", "dirichletInhom", "dynamics";
      params->GetList(keyVec, attrVec, valVec, fncnames_id_);
    }
    else if (analysistype_ == HARMONIC) {
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
    for ( Integer k = fncnames_id_.GetSize(); k < bcs_id_.GetSize(); k++ ) {
      fncnames_id_.Push_back( "none" );
    }
    if (dofspernode_ > 1) {
      keyVec = pdename_, "bcsAndLoads", "dirichletHom", "dof";
      params->GetList(keyVec, attrVec, valVec, homDirichDof_);
   
      keyVec = pdename_, "bcsAndLoads", "dirichletInhom", "dof";
      params->GetList(keyVec, attrVec, valVec, inhomDirichDof_);
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
  }


  void BasePDE::ReadMaterialData() {
    ENTER_FCN( "BasePDE::ReadMaterialData" );

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
    params->GetList( "name", subdomName, "domain", "region" );
    params->GetList( "material", subdomMaterial, "domain", "region" );
  
    params->Get("format", outformat, "output");

    if (outformat!="database") {
      // Query name of file with material data
      params->Get( "file", matFileName, "materialData" );
    
      // Generate new material reader
      LoadMaterialDataFile loadMaterialFile( matFileName.c_str() );
      
      // Load material data for subdomains on which this PDE lives
      // from data file
      for( Integer i = 0; i < subdoms_.GetSize(); i++ ) {
        for( Integer k = 0; k <= subdomName.GetSize(); k++ ) {
          if( subdoms_[i] == subdomName[k] ){
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
      for( Integer i = 0; i < subdoms_.GetSize(); i++ ) {
        for( Integer k = 0; k <= subdomName.GetSize(); k++ ) {
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


  Integer BasePDE::GetNumRestraints( const Integer level ) {

    ENTER_FCN( "BasePDE::GetNumRestraints" );
    
    Integer res = 0;
    Integer i;

    for ( i = 0; i < bcs_hd_.GetSize(); i++ ) {
      res+=ptBCs_->GetNumNodesLevel(bcs_hd_[i]);
    }

    for (i=0; i<bcs_id_.GetSize(); i++) {
      res+=ptBCs_->GetNumNodesLevel(bcs_id_[i]);
    }

    return res;
  }


  // **********************
  //   Default Destructor
  // **********************
  BasePDE::~BasePDE() {

    ENTER_FCN( "BasePDE::~BasePDE" );
    // ATTENTION: Dummy value for as_id!!!!!!!!!!!!!!!!!!!!!!!!!!
    DeleteAlgSys(0);

    if (assemble_) {
      delete assemble_;
    }
  
    // ATTENTION: Dummy value for as_id!!!!!!!!!!!!!!!!!!!!!!!!!!
    if (sol_) {
      delete sol_;
    }

    if (eqnData_)
      delete eqnData_;

    if ( TS_alg_)
      delete TS_alg_;
  
    if ( materialData_ != NULL)
      delete[] materialData_;
  }


  // ======================================================
  // ALGSYS SECTION (SOLVER, ...) 
  // ======================================================
  void BasePDE::SetAlgSys() {

    ENTER_FCN( "BasePDE::SetAlgSys" );

    // Set parameter for solver and preconditioner

    (*cla) <<  "--- PDE: " << pdename_ << " ---" << std::endl;

    // Set parameters for OLAS
    std::string amExpert;
    params->Get( "override", amExpert, "expert" );
    CFSOLASParams::SetParams( pdename_, params, olasParams_,(amExpert=="yes"));

    // Set the graph type used for the system matrices
    assemble_->SetupMatrixGraph(eqnData_->GetNumEQNs());

    // Allocate the necessary matrices as well as solver and preconditioner
    CreateMatrices_Solver();
  }


  void BasePDE::CreateMatrices_Solver() {

    ENTER_FCN( "BasePDE::CreateMatrices_Solver" );

    // create and initialize matrices 
    assemble_->CreateMatrices();
    assemble_->InitMatrices();

    // create solver and preconditioner
    algsys_->CreateSolver();

    algsys_->CreatePrecond();

    // now reset AlgebraicSystem 
    algsys_->InitRHS();
    algsys_->InitSol();
  }


  // ======================================================
  // COUPLING SECTION
  // ======================================================

  void BasePDE::CalcInputCoupling() {

    ENTER_FCN( "BasePDE::CalcInputCoupling" );

    std::string errMsg;
    StdVector<Integer> * nodes;
    CFSVector * val;
    Integer pdeNode, eqnNr,eqnDof;
    Integer couplingDof;
    Boolean clearCoords = TRUE;

    // Reset counter for boundary conditions
    couplingBCsCounter_ = 0;
  
    // Outer loop over all INPUT coupling terms
    for (Integer i=0; i<ptCoupling_->GetNumInputCouplings(); i++) {

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
          
        for (Integer j=0; j<nodes->GetSize(); j++)
          for (Integer dof=0; dof<ptCoupling_->GetInputDof(i); dof++) {
            pdeNode = eqnData_->Mesh2PDENode((*nodes)[j]);
            if (pdeNode==-1) {
              errMsg =  pdename_;
              errMsg += "PDE: Coupling node Nr. ";
              errMsg += Info->GenStr((*nodes)[j]);
              errMsg += " is not in contained in list of my subdomains!";
              Error(errMsg.c_str(), __FILE__, __LINE__);
            }
            deltCoords_(dof,pdeNode-1) = help[dof + j*dim_];

          }
        break;

        // -------------------
        // RHS COUPLING
        // -------------------
      case RHS:
        ptCoupling_->GetInputNodes(i, nodes);
          
        //for (Integer dof=0; dof<ptCoupling_->GetInputDof(i); dof++)
        for ( Integer dof = 0; dof < couplingDof; dof++ ) {
          for ( Integer j = 0; j < nodes->GetSize(); j++ ) {
            eqnData_->Node2EQN((*nodes)[j],dof+1,eqnNr,eqnDof);
            if ( eqnNr != 0 ) {
              algsys_->SetNodeRHS(help[dof+couplingDof*j], eqnNr, eqnDof);
            }
          }
        }
        break;

        // -----------------------
        // InhomDirichlet COUPLING
        // -----------------------
      case ID_BC:
          
        // Set flag that the boundary conditions have to be incorporated
        updateCouplingBCs_ = TRUE;

        ptCoupling_->GetInputNodes(i, nodes);

        for ( Integer dof = 0; dof < ptCoupling_->GetInputDof(i); dof++ ) {
          for ( Integer j = 0; j < nodes->GetSize();
                j++, couplingBCsCounter_++) {

            eqnData_->Node2EQN((*nodes)[j],dof+1,eqnNr,eqnDof);

            if (eqnNr==0) {
              Error( "The specified coupling node has no equation number",
                     __FILE__, __LINE__ );
            }

            algsys_->SetDirichlet( couplingBCsCounter_ + 1, eqnNr,
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


  // ======================================================
  // GRID SECTION (Meshing, ...) 
  // ======================================================

  void BasePDE::GetElemCoords( const StdVector<Integer> connect, 
                               Matrix<Double> &coordMat, 
                               const Integer level ) {

    ENTER_FCN( "BasePDE::GetElemCoords" );

    Integer pdeNode;
    ptgrid_->GetCoordNodesElemMat(connect, coordMat, level);
  
    if (deltCoords_.GetSizeRow() != 0 && geoUpdate_ == TRUE) {
      for (Integer i=0; i<coordMat.GetSizeRow(); i++) {
        for (Integer j=0; j<coordMat.GetSizeCol(); j++) {
          pdeNode = eqnData_->Mesh2PDENode(connect[j]);
          //coordMat(i,j) += deltCoords_(i,mesh2PDENode_ [connect[j] - 1] - 1);
          coordMat(i,j) += deltCoords_(i, pdeNode - 1);
        }
      }
    }
  }


  // ======================================================
  // ADAPTIVITY SECTION 
  // ======================================================

#ifdef ADAPTGRID
  void BasePDE::RefineMesh( const Integer level ) {

    ENTER_FCN( "BasePDE::RefineMesh" );
  
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


  Boolean BasePDE::TestError(const Integer level) {

    ENTER_FCN( "BasePDE::TestError" );

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
  void BasePDE::ConstructorError() {

    ENTER_FCN( "BasePDE::ConstructorError" );

    if (ptError_) delete ptError_;
  
    ptError_=new SpaceErrorEstimator();
    ptError_->Init(this);
  }


  void BasePDE::WriteErrorInfo(WriteResults * ptmeshes) {
    ptmeshes->WriteElemSolution(errorMap_,0,0,"ERR-errorMap");
    ptmeshes->WriteElemSolution(markingElems_,0,0,"ERR-markedElems");
  }

#endif // AdaptGrid


  Double BasePDE::GetFracDampMatrixCoeff(Integer actSD) {

    ENTER_FCN( "BasePDE::GetFracDampMatrixCoeff" );

    Double coeff;

    // pre factor of fractional derivative (same for all algorithms)
    Double y  = materialData_[actSD].GetDampingBeta();
    Double dt = TS_alg_->GetTimeStep();

    coeff = exp(-(y-1.0)*log(dt));

    // needed for formulation with only MASS and STIFFNESS matrix
    // pre factor of Newmark time stepping scheme
    Double beta = TS_alg_->GetNewmarkBeta();
    coeff *= 1.0 / (beta*dt*dt);

    return coeff;
  }

  void BasePDE::GetSolVecOfElement( Vector<Double>& elemSol,
                                    StdVector<Integer>& connecth ) {

    ENTER_FCN( "BasePDE::GetSolVecOfElement" );

    // displacement of element nodes
    elemSol.Resize(dofspernode_ * connecth.GetSize());
    elemSol.Init(0);
    Integer eqnNr, eqnDof;
    Integer dofsPerEQN = eqnData_->GetNumDofsPerEQN();

    NodeStoreSol<Double> * solhelp = dynamic_cast<NodeStoreSol<Double>*>(sol_);
    Vector<Double> sol = solhelp->GetAlgSysVector();
  
    for(Integer actNode=0; actNode<connecth.GetSize(); actNode++) {
      for(Integer actDof=0; actDof < dofspernode_; actDof++) {
        eqnData_->Node2EQN(connecth[actNode],actDof+1,eqnNr,eqnDof);
        if (eqnNr!= 0) {
          elemSol[actDof + actNode*dofspernode_] =
            sol[eqnDof-1 + dofsPerEQN*(abs(eqnNr-1))];
        }
        else {
          elemSol[actDof + actNode*dofspernode_] = 0.0;
        }
      }
    }
  }


  void BasePDE::GetDerivSolVecOfElement(Vector<Double>& sol,
                                        StdVector<Integer>& connecth) {

    ENTER_FCN( "BasePDE::GetDerivSolVecOfElement" );

    // displacement of element nodes
    sol.Resize(dofspernode_ * connecth.GetSize());
    sol.Init(0);
    Integer eqnNr, eqnDof;
    Integer dofsPerEQN = eqnData_->GetNumDofsPerEQN();
  
    if (analysistype_ == TRANSIENT) {
      const Vector<Double> & sol_der1 = getS1();
    
      for( Integer actNode = 0; actNode < connecth.GetSize(); actNode++ ) {
        for( Integer actDof = 0; actDof < dofspernode_; actDof++ ) {
          eqnData_->Node2EQN(connecth[actNode],actDof+1,eqnNr,eqnDof);
          if ( eqnNr != 0 ) {
            sol[actDof + actNode*dofspernode_] =
              sol_der1[eqnDof-1 + dofsPerEQN*(abs(eqnNr-1))];
          }
          else {
            sol[actDof + actNode*dofspernode_] = 0.0;
          }
        }
      }
    }
  }

  void BasePDE::GetDeriv2SolVecOfElement( Vector<Double>& sol,
					  StdVector<Integer>& connecth ) {

    ENTER_FCN( "BasePDE::GetDeriv2SolVecOfElement" );

    // displacement of element nodes
    sol.Resize(dofspernode_ * connecth.GetSize());
    sol.Init(0);
    Integer eqnNr, eqnDof;
    Integer dofsPerEQN = eqnData_->GetNumDofsPerEQN();

    if (analysistype_ == TRANSIENT) {
      const Vector<Double> & sol_der2 = getS2();
    
      for( Integer actNode = 0; actNode < connecth.GetSize(); actNode++ ) {
	for( Integer actDof = 0; actDof < dofspernode_; actDof++ ) {
	  eqnData_->Node2EQN(connecth[actNode],actDof+1,eqnNr,eqnDof);
	  if (eqnNr!= 0) {
	    sol[actDof + actNode*dofspernode_] =
	      sol_der2[eqnDof-1 + dofsPerEQN*(abs(eqnNr-1))];
	  }
	  else {
	    sol[actDof + actNode*dofspernode_] = 0.0;
	  }
	}
      }
    }
  }


  void BasePDE::GetDerivSolOfElement( Matrix<Double>& sol,
                                        StdVector<Integer>& connect_PDE ) {

    ENTER_FCN( "BasePDE::GetDerivSolOfElement" );

    // displacement of element nodes
    sol.Resize(dofspernode_, connect_PDE.GetSize());

    const Vector<Double>& sol_der1 = getS1();  

    for( Integer actNode = 0; actNode < connect_PDE.GetSize(); actNode++ ) {
      for( Integer actDof = 0; actDof < dofspernode_; actDof++) {
	sol[actDof][actNode] =
	  sol_der1[actDof + dofspernode_*(connect_PDE[actNode]-1)];
      }
    }
  }


  void BasePDE::sortStresses(Vector<Double>& unsorted, Vector<Double>& sorted){

    ENTER_FCN( "BasePDE::SortStresses" );

    //soring according to capa (unv) notation
    //our notation is Voit: Txx Tyy Tzz Tyz Txz Txy

    if (subType_ == "3d") {
      //Txx Txy Tyy Txz Tyz Tzz
      sorted[0] = unsorted[0];
      sorted[1] = unsorted[5];
      sorted[2] = unsorted[1];
      sorted[3] = unsorted[4];
      sorted[4] = unsorted[3];
      sorted[5] = unsorted[2];
    }
    else if (subType_ == "axi") {

      //Tphiphi 0 Trr 0 Trz Tzz
      sorted[0] = unsorted[3];
      sorted[1] = 0.0;
      sorted[2] = unsorted[0];
      sorted[3] = 0.0;
      sorted[4] = unsorted[2];
      sorted[5] = unsorted[1];
    }
    else {
      //0 0 Tyy 0 Tyz Tzz
      sorted[0] = 0.0;
      sorted[1] = 0.0;
      sorted[2] = unsorted[0];
      sorted[3] = 0.0;
      sorted[4] = unsorted[2];
      sorted[5] = unsorted[1];
    }
  }


  void BasePDE::sortStresses( Vector<Complex> &unsorted,
                              Vector<Complex> &sorted ) {
    ENTER_FCN( "BasePDE::SortStresses" );

    //soring according to capa (unv) notation
    //our notation is Voit: Txx Tyy Tzz Tyz Txz Txy

    if (subType_ == "3d") {

      //Txx Txy Tyy Txz Tyz Tzz
      sorted[0] = unsorted[0];
      sorted[1] = unsorted[5];
      sorted[2] = unsorted[1];
      sorted[3] = unsorted[4];
      sorted[4] = unsorted[3];
      sorted[5] = unsorted[2];
    }

    else if (subType_ == "axi") {

      //Tphiphi 0 Trr 0 Trz Tzz
      sorted[0] = unsorted[3];
      sorted[1] = 0.0;
      sorted[2] = unsorted[0];
      sorted[3] = 0.0;
      sorted[4] = unsorted[2];
      sorted[5] = unsorted[1];
    }

    else {

      //0 0 Tyy 0 Tyz Tzz
      sorted[0] = 0.0;
      sorted[1] = 0.0;
      sorted[2] = unsorted[0];
      sorted[3] = 0.0;
      sorted[4] = unsorted[2];
      sorted[5] = unsorted[1];
    }
  }


  //stores an algsys_ vector into a StdVector and returns that L2-norm
  void BasePDE::StoreAlgsysToVec(Vector<Double>& vec, Double * pt) {

    ENTER_FCN( "BasePDE::StoreAlgsysToVec" );

    //const Integer numElems = numPDENodes_ * dofspernode_;
    Integer numElems = eqnData_->GetNumEQNs() * eqnData_->GetNumDofsPerEQN();

    for (Integer i=0; i<numElems; i++) {
      vec[i] = pt[i];
    }
  }

} // end of namespace

