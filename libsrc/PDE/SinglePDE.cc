#include "PDE/SinglePDE.hh"

// for coordinate handling
#include "Domain/domain.hh"
#include "Utils/coordSystem.hh"

// header for Paramhandling
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/ParamHandling/CFSOLASParams.hh"
#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"


// header for scripting
#ifdef TCL_INTERFACE
#include "DataInOut/Scripting/cfsmessenger.hh"
#endif

// header for Materialhandling
#include "DataInOut/MaterialHandler.hh"

// header for Solvestep and assemble
#include "Driver/stdSolveStep.hh"
#include "Driver/assemble.hh"
#include "Driver/basedriver.hh"

// header for iterative coupling
#include "CoupledPDE/pdecoupling.hh"

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "muParser.h"

namespace CoupledField {


  SinglePDE::SinglePDE( Grid *aptgrid, WriteResults *aOutFile,
                        TimeFunc *aptTimeFunc )
    :  StdPDE(aptgrid, aOutFile, aptTimeFunc) {
  
    ENTER_FCN( "BasePDE::BasePDE" );
  
    nonLin_ = false;
    incStopCrit_ = 1e-2;
    residualStopCrit_ = 1e-3;
  
    // =====================================================================
    // set file pointers
    // =====================================================================
    ptTimeFunc_ = aptTimeFunc;
    assemble_   = NULL;
    algsys_     = NULL;

    // =====================================================================
    // set analysis parameters
    // =====================================================================
    complexFormat_ = AMPLITUDE_PHASE; // or REAL_IMAG
    couplingBCsCounter_ = 0;
    updateCouplingBCs_ = false;
    dim_ = ptgrid_->GetDim();
    iterCoupledCounter_ = 0;
    effectiveMass_ = false;

    // =====================================================================
    // set postprocessing parameters
    // =====================================================================
    hasOutput_      = false;
    saveSol_        = false;
    saveDeriv1_     = false;
    saveDeriv2_     = false;
    saveSolHist_    = false;
    saveDeriv1Hist_ = false;
    saveDeriv2Hist_ = false;


    // =====================================================================
    // set miscellaneous parameters
    // =====================================================================
    pdeId_ = NO_PDE_ID;
    isDirectCoupled_  = false;
    numCouplingBcs_ = 0;

    // Obtain mathParser handler
    mHandle_ = domain->GetMathParser()->GetNewHandle();

    // Register functions for scripting
    RegisterFunctions();
  }


  // **************
  //   Destructor
  // **************
  SinglePDE::~SinglePDE() {

    ENTER_FCN( "SinglePDE::~SinglePDE" );

    // Delete algebraic system only if
    // PDE is not direct coupled
    if ( isDirectCoupled_ == false ) {
      if( needsAlgsys_ ) {
        delete algsys_;
      }
      delete solveStep_;
      delete TS_alg_;
      delete assemble_;
    }


    delete sol_;
    delete solVec_;


    std::map<RegionIdType, BaseMaterial*>::iterator it;
    for ( it = materials_.begin(); it != materials_.end(); it++ ) {
      delete it->second;
    }
    materials_.clear();

  }


  // ********
  //   Init
  // ********
  void SinglePDE::Init( UInt bcSequenceIndex,
                        std::string  bcSequenceTag ) {
    ENTER_FCN( "SinglePDE::Init()" );

    StdVector<RegionIdType> allIDs;

    bcSequenceIndex_ = bcSequenceIndex;
    bcSequenceTag_ = bcSequenceTag;
  
    // =====================================================================
    // get regions/subdomains for PDE
    // =====================================================================
    StdVector<std::string> regionNames, surfaceNames;
    params->GetList( "name", regionNames, pdename_, "region" );
    ptgrid_->RegionNameToId( subdoms_, regionNames );
    params->GetList( "name", surfaceNames, pdename_, "surface" );
    ptgrid_->RegionNameToId( surfdoms_, surfaceNames );

    // Create vector of all IDs
    allIDs = subdoms_;
    for (UInt i = 0; i < surfdoms_.GetSize(); i++) {
      allIDs.Push_back(surfdoms_[i]);
    }
    
    // output to info-file
    Info->PrintF( pdename_, "The %s PDE lives on the following regions:\n",
                  pdename_.c_str());
    for ( UInt k = 0; k < regionNames.GetSize(); k++ ) {
      Info->PrintF( pdename_, "%s, ID = %i\n",
                    regionNames[k].c_str(), subdoms_[k] );
    }
    Info->PrintF( "", "\n" );
    Info->PrintF( pdename_, "The %s PDE has the following surface regions:\n",
                  pdename_.c_str());
    for ( UInt k = 0; k < surfaceNames.GetSize(); k++ ) {
      Info->PrintF( pdename_, "%s, ID = %i\n",
                    surfaceNames[k].c_str(), surfdoms_[k] );
    }
    Info->PrintF( "", "\n" );

    // Generate a fitting algebraic system only if PDE is NOT
    // direct coupled
    if( needsAlgsys_ == true ) {
      if ( isDirectCoupled_ == false) {
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
    }
    
    // =====================================================================
    // Get type of analysis
    // =====================================================================
    
    // Construct vectors for restricted search parameter
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    analysistype_ = domain->GetDriver()->GetAnalysisType( pdename_ );

    // NOTE: The concept of isAlwaysStatic bites with Direct Coupling
    //       and must be re-designed
    if ( isAlwaysStatic_ == true &&
         analysistype_ == TRANSIENT ) {
      analysistype_ = STATIC;
    }

    switch( analysistype_ ) {
    case STATIC:
      isComplex_ = false;
      break;

    case TRANSIENT:
      isComplex_ = false;
      break;
      
    case HARMONIC:
      isComplex_ = true;
      break;
      
    case EIGENFREQUENCY:
      isComplex_ = false;
      break;
      
    default:
      
      (*error) << "SinglePDE::Init: AnalysisType '" << analysistype_
               << "' is not supported";
      Error( __FILE__, __LINE__ );
    }

    // Create new assemble class with according analysistype
    if( isDirectCoupled_ == false && needsAlgsys_ == true) {
      assemble_ = new Assemble( algsys_, analysistype_, ptTimeFunc_ );
    }
        
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
    ReadDampingInformation( );
      

    // =====================================================================
    // read in NonLinearities
    // =====================================================================
    InitNonLin();

    // =====================================================================
    // read in material data
    // =====================================================================
    ReadMaterialData();

    // =====================================================================
    // initialize EQN-object and Storeresults class
    // =====================================================================

    // Name of linear system depends of coupling type
    if ( isDirectCoupled_ == false ) {
      valVec  = "", pdename_, "";
    }
    else {
      valVec  = "", "direct", "";
    }

    // What type of equation numbering does the user want?
    std::string typeOfNumbering;
    keyVec  = "linearSystems", "system", "setup", "eqnNumbering";
    attrVec = "", "name", "";
    params->Get( keyVec, attrVec, valVec, typeOfNumbering );

    // How do we want to treat inhomogeneous Dirichlet boundary conditions?
    bool usePenalty = true;
    {
      std::string aux;
      keyVec  = "linearSystems", "system", "setup", "idbcHandling";
      attrVec = "", "name", "";
      params->Get( keyVec, attrVec, valVec, aux );
      usePenalty = aux == "penalty" ? true : false;
      Info->PrintF( pdename_, "Treating IDBCs using '%s' approach\n",
                    aux.c_str() );
    }
    
    // Create a new equation map
    eqnMap_ = shared_ptr<EqnMap>(new EqnMap( ptgrid_, pdeId_, !usePenalty ));


    // =====================================================================
    // read in boundary conditions
    // =====================================================================
    ReadBCs();
    ReadSpecialBCs();
    //
#ifdef TCL_INTERFACE
    // Trigger event for scripting 
    StdVector<std::string> args;
    args.Push_back( pdename_ );
    messenger->TriggerEvent( CFSMessenger::CFS_ReadBCs, args );
#endif   

    // =====================================================================
    // Create time stepping algorithm
    // =====================================================================
    if ( analysistype_ == TRANSIENT && 
         isDirectCoupled_ == false) {
      SETPROFILE("Before Definition of Timestepping");
      InitTimeStepping();
      SETPROFILE("After Definition of TimeStepping");
    }

    // =====================================================================
    // define the integrators for PDE and initialize eqn object
    // =====================================================================

    // Call initialization of (bi)linear integrators
    DefineIntegrators();

    // Print information about defined integrators
#ifdef DEBUG
    if( !isDirectCoupled_ && needsAlgsys_ == true ) {
      assemble_->PrintInfo( *debug );
    }
#endif

    // Finish equation mapping
    eqnMap_->SetHomoDirichletBCs ( hdBcs_ );
    eqnMap_->SetInhomDirichletBCs( idBcs_ );
    eqnMap_->SetConstraints( constraints_ );
    eqnMap_->Finalize();
    
     // Report results to logfile
    Info->PrintF( pdename_, "Linear system will have %d equations\n\n",
                  eqnMap_->GetNumEqns() );

    // Check if eqnmap should be printed onto screen
    if ( commandLine->GetShowEqnMap() == true ) {
      eqnMap_->Print( (*cla) );
    }

#ifdef DEBUG
    eqnMap_->Print(*debug);
#endif
    numPDENodes_ = eqnMap_->GetNumLocalNodes();
    numElems_ = eqnMap_->GetNumLocalElems();

    // Initialize Storesolution class
    sol_->SetNumSolutions(results_.GetSize());
    sol_->SetNumNodes(numPDENodes_);
    for (UInt iSol=0; iSol<results_.GetSize(); iSol++) {
      sol_->SetSolutionType(results_[iSol]->resultType,iSol);
      sol_->SetNumDofs( results_[iSol]->dofNames.GetSize(), 
                        results_[iSol]->resultType );
    }

    // Note: this is only a temporary solution
    sol_->SetPtrEQNData(eqnMap_.get(), ptgrid_);
    sol_->SetRegions( subdoms_ );
    sol_->Init(); 
    

    SETPROFILE("Before Resizing StoreSol");
    solVec_->Resize( eqnMap_->GetNumEqns() );
    if( isComplex_ ) {
      solVec_->Init( Complex(0.0,0.0) );
    } else {
      solVec_->Init( 0.0 );
    }
    SETPROFILE("After Resizing StoreSol");
    if ( analysistype_ == HARMONIC || analysistype_ == MULTIHARMONIC ) {
      sol_->SetAlgSysDataPointer(solVec_->GetSize(), 
                                 dynamic_cast<Vector<Complex>&>(*solVec_).GetPointer() );
    } else {
      sol_->SetAlgSysDataPointer(solVec_->GetSize(), 
                                 dynamic_cast<Vector<Double>&>(*solVec_).GetPointer() );

    }
    

    // =====================================================================
    // define which solution types have to be saved
    // =====================================================================
    ReadStoreResults();

#ifndef MpCCI
    // check, if any output is calculated at all
    if ( hasOutput_ == false ) {
      (*warning) << "There was no output specified at all for PDE '"
                 << pdename_
                 << "' Please check your .xml-file, if this is really what "
                 << "you want!";
      Warning( __FILE__, __LINE__ );
    }
#endif


    // Set information at sol_ object
    sol_->SetResult( results_[0] );


    PreparePDE4Computation();

    //! Define step solution driver
    if ( isDirectCoupled_ == false ) {
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
    
    // 1.) Homogeneous boundary condition
    Info-> WriteHomDirBC( pdename_, hdBcs_ );
    
    // 2.) Inhomogeneous boundary conditions
    Info->WriteInhomDirBC( pdename_, idBcs_ );

    // 3.) Inhom. Neumann boundary conditions
    Info->WriteInhomNeuBC( pdename_, inBcs_ );
    
    // 4.) Constraints
    Info->WriteConstraints( pdename_, constraints_ );

    // 5.) Loads
    Info->WriteLoad( pdename_, loads_ );
      }



  // **********
  //   SetBCs
  // **********
  void SinglePDE::SetBCs( const Double time ) {
    
     ENTER_FCN( "SinglePDE::SetBCs" );
  
     // Trigger setting of BC from script file
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
     messenger->TriggerEvent( CFSMessenger::CFS_SetBCs, context );
 #endif
  
     UInt dof;
     Double val, val_tfunc;
     StdVector<UInt> nodes;
     UInt bcNum = 0;
     Integer eqnNr; 
     Vector<Double> globCoord;
   
     // get global coordinate system and math parser
     CoordSystem * coosy = domain->GetCoordSystem();
     MathParser * parser = domain->GetMathParser();
  
     // set offset due to IDBC comming from input coupling
     if ( isIterCoupled_ ) {
       bcNum = couplingBCsCounter_;
     }
     else {
       bcNum = 0;
     }
     
     // ---------------------------
     // INHOMOGENEOUS DIRICHLET BC
     // ---------------------------
  
     Double phase = 0.0;
     StdVector<Double>  val_tfunc_vec, dirVal_vec;
     
     for ( UInt i = 0; i < idBcs_.GetSize(); i++ ) {
       
       // Get grip of actual idBC
       InhomDirichletBc const & actBc = *(idBcs_[i]);

       dof = actBc.dof;
       std::string const & dynamics = actBc.dynamics;
       
       // Get the correct time function value
       val_tfunc = 1.0;
       if ( ptTimeFunc_->GetmaxTimeFnc() > 0 &&
            (analysistype_ != HARMONIC || analysistype_ != MULTIHARMONIC) ) {
           val_tfunc=ptTimeFunc_->TimeFuncAtTime(time, dynamics);
       }
       
       // Get EntityIterator
       EntityIterator it = actBc.entities->GetIterator();

       for ( it.Begin(); !it.IsEnd(); it++ ) {

         eqnNr = eqnMap_->GetEqn( *actBc.result, it, dof );
      

         // If iterator points to a node, pass the current coordinate
         // to the parser
         if( it.GetType() == EntityList::NODE_LIST ) {
           
           // Get node coordinate
           ptgrid_->GetNodeCoordinate( globCoord, it.GetNode() );
           parser->SetCoordinates( mHandle_, *coosy, globCoord );
         }
         
         // Now evaluate value of IDBC
         parser->SetExpr( mHandle_, actBc.value );
         val = parser->Eval( mHandle_ ) * val_tfunc;

         // Sanity check. This should not happen, but might appear
         // in the case that the same node/dof belongs to a region
         // with hom. and a region with inhom. Dirichlet BCs. This
         // problem was already encountered!
         if (eqnNr == 0) {

           (*error) << "Got eqn number 0 for inhom Dirichlet BC! "
                    << "Probably you have a node/dof that belongs to both "
                    << "a region with hom. and one with inhom. Dirichlet BCs."
                    << " Go check your .mesh file!";
           Error( __FILE__, __LINE__ );
         }

         // Transform Dirichlet boundary conditions for effmass-formulation
         if (effectiveMass_) {
           val = TS_alg_->DirichletBC4EffMassMatrix(val,eqnNr);
         }
         
         // Case of complex-valued entries
         if (analysistype_ == HARMONIC || analysistype_ == MULTIHARMONIC) {

           parser->SetExpr( mHandle_, actBc.phase );
           phase = parser->Eval( mHandle_ );
           //Complex complexValue( val * std::cos( phase / 180 * PI ),
           //                      val * std::sin( phase / 180 * PI ) );
           Complex complexValue( val * cos( phase / 180 * PI ),
                                 val * sin( phase / 180 * PI ) );
           algsys_->SetDirichlet( bcNum + 1, pdeId_, eqnNr, complexValue,
                                  1 );
         }
         else {
 	  //	  std::cout << "IHDBC val=" << val << std::endl;
           algsys_->SetDirichlet( bcNum + 1, pdeId_, eqnNr, val, 1 );
         }
         bcNum++;
       }
     }
  }

 

  void SinglePDE::ReadBCs() {

    ENTER_FCN( "SinglePDE::ReadBCs" );

    // vectors for parameter handling
    StdVector<std::string> keyVec;
    StdVector<std::string> attrVec;
    StdVector<std::string> valVec;

    // =====================================================================
    // homogeneous Dirichlet BC
    // =====================================================================
    StdVector<std::string> names, dofs, entTypes;
    
    // Get names of node sets for homogeneous Dirichlet boundary conditions
    //   params->GetList( "name", bcs_hd_, pdename_, "dirichletHom"   );
    keyVec = pdename_, "bcsAndLoads", "dirichletHom", "name";
    attrVec = "", "tag", "";
    valVec = "", bcSequenceTag_, "";
    params->GetList(keyVec, attrVec, valVec, names);

    keyVec.Last() = "dof";
    params->GetList(keyVec, attrVec, valVec, dofs);

    keyVec.Last() = "entityType";
    params->GetList(keyVec, attrVec, valVec, entTypes);

    // Create homogeneous boundary conditions
    for( UInt i = 0; i < names.GetSize(); i++ ) {
      shared_ptr<HomDirichletBc> actBc ( new HomDirichletBc );
      shared_ptr<EntityList> actList =
        ptgrid_->GetEntityList( EntityList::StringToType(entTypes[i]), 
                                names[i] );
      actBc->entities = actList;
      actBc->result = results_[0];
      actBc->eqnMap = eqnMap_;
      if( dofs.GetSize() == 0 ) {
        actBc->dof = 1; 
      } else {
        actBc->dof = results_[0]->GetDofIndex( dofs[i] );
      }

      // add definition
      hdBcs_.Push_back( actBc );
    }

    // =====================================================================
    // inhomogeneous Dirichlet BC
    // =====================================================================
    StdVector<std::string> idName, idDof, idValue, idPhase, idDynamics;
    StdVector<std::string> idType;
    
    // Get names of node sets, values and filenames for inhomogenous
    // Dirichlet boundary conditions
    keyVec = pdename_, "bcsAndLoads", "dirichletInhom", "name";
    params->GetList(keyVec, attrVec, valVec, idName);

    keyVec.Last() = "entityType";
    params->GetList(keyVec, attrVec, valVec, idType);

    keyVec.Last() = "value";
    params->GetList(keyVec, attrVec, valVec, idValue);

    keyVec.Last() = "dof";
    params->GetList(keyVec, attrVec, valVec, idDof);
    
    keyVec.Last() = "dynamics";
    params->GetList(keyVec, attrVec, valVec, idDynamics);

    keyVec.Last() = "phase";
    params->GetList(keyVec, attrVec, valVec, idPhase);

    // Create inhomogeneous boundary conditions
    for( UInt i = 0; i < idName.GetSize(); i++ ) {
      shared_ptr<InhomDirichletBc> actBc ( new InhomDirichletBc );
      shared_ptr<EntityList> actList =
        ptgrid_->GetEntityList( EntityList::StringToType(idType[i]), 
                                idName[i] );
      actBc->entities = actList;
      actBc->result = results_[0];
      actBc->eqnMap = eqnMap_;
      if( idDof.GetSize() == 0 ) {
        actBc->dof = 1;
      } else {
        actBc->dof = results_[0]->GetDofIndex( idDof[i] );
      }
      actBc->value = idValue[i];
      actBc->phase = idPhase[i];
      actBc->dynamics = idDynamics[i];

      // add definition
      idBcs_.Push_back( actBc );
    }



    // =====================================================================
    // inhomogeneous Neumann BC
    // =====================================================================
    StdVector<std::string> inName, inDof, inValue, inPhase, inDynamics;
    StdVector<std::string> inType;
    
    // Get names of node sets, values and filenames for inhomogenous
    // Neumann boundary conditions
    keyVec = pdename_, "bcsAndLoads", "neumannInhom", "name";
    params->GetList(keyVec, attrVec, valVec, inName);

    keyVec.Last() = "entityType";
    params->GetList(keyVec, attrVec, valVec, inType);

    keyVec.Last() = "value";
    params->GetList(keyVec, attrVec, valVec, inValue);

    keyVec.Last() = "dof";
    params->GetList(keyVec, attrVec, valVec, inDof);
    
    keyVec.Last() = "dynamics";
    params->GetList(keyVec, attrVec, valVec, inDynamics);

    keyVec.Last() = "phase";
    params->GetList(keyVec, attrVec, valVec, inPhase);

    // Create inhomogeneous Neumann boundary conditions
    for( UInt i = 0; i < inName.GetSize(); i++ ) {
      shared_ptr<InhomNeumannBc> actBc ( new InhomNeumannBc );
      shared_ptr<EntityList> actList =
        ptgrid_->GetEntityList( EntityList::StringToType(inType[i]), 
                                inName[i] );
      actBc->entities = actList;
      actBc->result = results_[0];
      actBc->eqnMap = eqnMap_;
      if( inDof.GetSize() == 0 ) {
        actBc->dof = 1;
      } else {
        actBc->dof = results_[0]->GetDofIndex( inDof[i] );
      }
      actBc->value = inValue[i];
      actBc->phase = inPhase[i];
      actBc->dynamics = inDynamics[i];

      // add definition
      inBcs_.Push_back( actBc );
    }

    // =====================================================================
    // Constraint Conditions
    // =====================================================================
    StdVector<std::string> csName, csEntityType;


    // read name of constraints
    keyVec = pdename_, "bcsAndLoads", "constraint", "name";
    params->GetList( keyVec, attrVec, valVec, csName );

    keyVec.Last() = "entityType";
    params->GetList(keyVec, attrVec, valVec, csEntityType);

    // Create constraint condition
    for( UInt i = 0; i < csName.GetSize(); i++ ) {
      shared_ptr<Constraint> actBc ( new Constraint );
      shared_ptr<EntityList> actList =
        ptgrid_->GetEntityList( EntityList::StringToType(csEntityType[i]), 
                                csName[i] );
      actBc->masterEntities = actList;
      actBc->slaveEntities = actList;
      actBc->masterDof = 1;
      actBc->slaveDof = 1;
      actBc->result = results_[0];
      actBc->eqnMap = eqnMap_;

      // add definition
      constraints_.Push_back( actBc );
    }

    // =====================================================================
    // Load definitions
    // =====================================================================
    
    StdVector<std::string> name, dof, value, phase, dynamics, type;

    attrVec = "", "tag", "";
    valVec = "", bcSequenceTag_, "";
  
    keyVec = pdename_, "bcsAndLoads", "load", "name";
    params->GetList(keyVec, attrVec, valVec, name);

    keyVec.Last() = "entityType";
    params->GetList(keyVec, attrVec, valVec, type);

    keyVec.Last() = "dof";
    params->GetList(keyVec, attrVec, valVec, dof);

    keyVec.Last() =  "value";
    params->GetList(keyVec, attrVec, valVec, value);

    keyVec.Last() = "phase";
    params->GetList(keyVec, attrVec, valVec, phase);

    keyVec.Last() =  "dynamics";
    params->GetList(keyVec, attrVec, valVec, dynamics);

    for( UInt i=0; i<name.GetSize(); i++ ) {
    
      shared_ptr<LoadBc> actLoad( new LoadBc );
      shared_ptr<EntityList> actList =
        ptgrid_->GetEntityList( EntityList::StringToType(type[i]),
                                name[i] );

      actLoad->entities = actList;
      actLoad->result = results_[0];
      actLoad->eqnMap = eqnMap_;
      if( dof.GetSize() == 0 ) {
        actLoad->dof = 1;
      } else {
        actLoad->dof = 
          results_[0]->GetDofIndex(dof[i]);
      }
      actLoad->value = value[i];
      actLoad->phase = phase[i];
      actLoad->dynamics = dynamics[i];
      loads_.Push_back( actLoad);
    }
    assemble_->AddLoads( loads_ );    
  }
  
  
  void SinglePDE::ReadMaterialData() {
    ENTER_FCN( "SinglePDE::ReadMaterialData" );

    std::string outformat="unverg";
    std::string matFileName;


    // Get list of subdomains and materials
    StdVector< std::string > subdomName, keyVec, attrVec, valVec;
    StdVector< std::string > subdomMaterial, subdomComposite, subdomCoordSys;
    StdVector< RegionIdType > subdomId;
    params->GetList( "name", subdomName, "domain", "region" );
    params->GetList( "material", subdomMaterial, "domain", "region" );
    params->GetList( "composite", subdomComposite, "domain", "region" );
    params->GetList( "refCoordSys", subdomCoordSys, "domain", "region" );

    // Convert region names to Ids
    ptgrid_->RegionNameToId( subdomId, subdomName );
    params->Get("format", outformat, "output");

    // Obtain pointer to materialHandler
    MaterialHandler * matLoader = NULL;
    matLoader = domain->GetMaterialHandler();
    std::string actRegionName;      
    
    // -------------------
    // NORMAL MATERIALS
    // -------------------
    for( UInt i = 0; i < subdoms_.GetSize(); i++ ) {
      for( UInt k = 0; k < subdomId.GetSize(); k++ ) {
        if( subdoms_[i] == subdomId[k] ){
          
          // Check if region contains a material
          if ( subdomMaterial[k] != "" ) {
            
            // Be verbose
            actRegionName = ptgrid_->RegionIdToName( subdomId[k] );
            Info->PrintF( pdename_, "Material '%s' for region '%s' (ID = %d) "
                          "follows\n", subdomMaterial[k].c_str(),
                          actRegionName.c_str(), subdomId[k] );
            // Read data
            materials_[subdoms_[i]] = matLoader->
              LoadMaterial( subdomMaterial[k], pdematerialclass_ );

            // Check if coordinate system is present
            if( subdomCoordSys[k] != "" ) {
              CoordSystem * actCoosy = 
                domain->GetCoordSystem(subdomCoordSys[k]);
              materials_[subdoms_[i]]->SetCoordSys( actCoosy );
            }
            break;
          }
        }
      }
    }
    
    // -------------------
    // COMPOSITE MATERIALS
    // -------------------
    Double startHeight;
    StdVector<std::string> compMaterials;
    StdVector<Double> compOrient, compThick;
    
    std::ostringstream out;
    for( UInt i = 0; i < subdoms_.GetSize(); i++ ) {
      for( UInt k = 0; k < subdomId.GetSize(); k++ ) {
        if( subdoms_[i] == subdomId[k] ){
          
          // Check if region contains a material
          if ( subdomComposite[k] != "" ) {
            
            actRegionName = ptgrid_->RegionIdToName( subdomId[k] );
            out << "Composite material '" << subdomComposite[k] << "' for "
                << "region '" << actRegionName << "' (ID = " << subdomId[k]
                << ") follows:\n";
            Info->PrintF( pdename_, out.str().c_str());
            out.str("");
            
            // Read data from parameter file

            // StartHeight
            keyVec = "domain", "composite", "startHeight";
            attrVec = "", "name";
            valVec = "", subdomComposite[k];
            params->Get( keyVec, attrVec, valVec, startHeight );
            
            // materials of laminae
            attrVec = "", "name", "";
            valVec = "", subdomComposite[k], "";
            keyVec = "domain", "composite", "lamina", "material";
            params->GetList( keyVec, attrVec, valVec, compMaterials );
            
            // thickness of laminae
            keyVec = "domain", "composite", "lamina", "thickness";
            params->GetList( keyVec, attrVec, valVec, compThick );
            
            // orientation of laminae
            keyVec = "domain", "composite", "lamina", "orientation";
            params->GetList( keyVec, attrVec, valVec, compOrient );
            
            // Check, if all vectors have same length
            if ( compMaterials.GetSize() != compThick.GetSize() ||
                 compMaterials.GetSize() != compOrient.GetSize() ) {
              Error( "Inconsistent definition of composite material",
                     __FILE__, __LINE__ );
            }
            
            // Create new lamina and fill ine materials and thicknesses
            Composite & myMat = compositeMaterials_[subdomId[k]];
            UInt numLaminae = compMaterials.GetSize();
            myMat.name = subdomComposite[k];
            myMat.zStart = startHeight;
            myMat.thickness = compThick;
            myMat.orientation = compOrient;
            
            // Read single materials and print data
            myMat.materials.Resize( numLaminae );
            for ( UInt iLam = 0; iLam < numLaminae; iLam++ ) {
              
              // Print information
              out << " --- Lamina " << iLam+1 << ": thickness = " 
                  << compThick[iLam]
                  << " m, orientation = " << compOrient[iLam] << "° ---\n";
              Info->PrintF( pdename_, out.str().c_str());
              out.str("");
              
              myMat.materials[iLam] = matLoader->
                LoadMaterial( compMaterials[iLam], pdematerialclass_ );
            } // over  laminae
          }
        }
      }
    }
  }


  
  // ======================================================
  // GET /SET  METHODS
  // ======================================================

  //! Activate the direct coupling
  void SinglePDE::SetDirectCoupling () {
                           
    ENTER_FCN( "SinglePDE::SetDirectCoupling" );
    
    isDirectCoupled_ = true;          
  }



  // ======================================================
  // ALGSYS SECTION (SOLVER, ...) 
  // ======================================================
  void SinglePDE::DefineAlgSys() {

    ENTER_FCN( "StdPDE::DefineAlgSys" );

    // First check if the PDE needs an algebraic system at all
    if( needsAlgsys_ == false ) {
      return;
    }


    // If PDE is not direct coupled then the PDE has to register
    // at the algebraic system and obtain an Id. 
    // Afterwards the matrix-graph has to be set up
    SETPROFILE("Before GraphSetupInit()");
    if ( isDirectCoupled_ == false ) {

      // Set linear system parameters for OLAS
      ReadOlasParams( pdename_ );

      // Initialize the matrix graph object
      algsys_->GraphSetupInit(1);

      algsys_->RegisterPDE( pdeId_, eqnMap_->GetNumEqns(),
                            eqnMap_->GetNumLastFreeDof() );

      //assemble_->SetPDEId( pdeId_ );
      solveStep_->SetPDEId( pdeId_ );
      
      
      // trigger the creation and assembly of the matrix graph
      algsys_->AssembleInit( pdeId_, pdeId_, false );
      assemble_->SetupMatrixGraph(pdeId_, pdeId_);
      algsys_->AssembleDone( pdeId_, pdeId_, false );      
    
      // finish the assembly of the matrix graph
      algsys_->GraphSetupDone();
      
      // obtain reordering of the matrix graph and pass it to the EQN-object.
      Integer *newOrder = algsys_->GetReordering( pdeId_ );
      eqnMap_->ReorderMapping( &newOrder );
    }

    // pass information about dofs, number of dirichlet equations
    // and constraints to the algebraic system
    algsys_->SetBlockSize( pdeId_, 1 );

    UInt numDir = eqnMap_->GetNumInHomDirichletEqns() + numCouplingBcs_;
    algsys_->SetNumDirichletBCs(pdeId_, numDir );

    // create matrices and solver object, if PDE is not direct coupled
    if ( isDirectCoupled_ == false ) {
      CreateMatrices_Solver();
    }

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
    marlingElems_.Init();

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


  bool SinglePDE::TestError(const Integer level) {

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

    if (totalErr > tolSpaceErr_) return true;
    else return false;
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
    Integer eqnNr;
    UInt couplingDof;

    // Determine maximal allowed equation number for algebraic system
    Integer maxAllowedEqn = (Integer)algsys_->GetDimension();

    // at first, check if this PDE is iterative coupled
    if (isIterCoupled_ == false)
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

        ptCoupling_->GetInputNodes(i, nodes);
        ptgrid_->SetNodeOffset( *nodes , help );

        break;

        // -------------------
        // RHS COUPLING
        // -------------------
      case RHS:
        ptCoupling_->GetInputNodes(i, nodes);
          
        for ( UInt dof = 0; dof < couplingDof; dof++ ) {
          for ( UInt j = 0; j < nodes->GetSize(); j++ ) {
            eqnNr = eqnMap_->GetNodeEqn( (*nodes)[j], dof+1 );
            if ( eqnNr != 0 && eqnNr <= maxAllowedEqn ) {
              algsys_->SetNodeRHS( help[ dof + couplingDof * j ], pdeId_,
                                   eqnNr, 1 );
            }

#ifdef DEBUG
            else if ( eqnNr > maxAllowedEqn ) {
              (*debug) << "SinglePDE::CalcInputCoupling: "
                       << "(" << pdename_ << ") "
                       << "Refused to pass "
                       << "eqnNr = " << eqnNr << " to SetNodeRHS(), since "
                       << "it execeeds numLastFreeDof = " << maxAllowedEqn
                       << std::endl;
            }
            else if ( eqnNr == 0 ) {
              (*debug) << "SinglePDE::CalcInputCoupling: "
                       << "(" << pdename_ << ") "
                       << "Refused to pass "
                       << "eqnNr = " << eqnNr << " to SetNodeRHS(), since "
                       << "it is fixed by hom. Dirichlet BC" << std::endl;
            }
#endif

          }
        }
        break;

        // -----------------------
        // InhomDirichlet COUPLING
        // -----------------------
      case ID_BC:

        // How do we want to treat inhomogeneous Dirichlet boundary
        // conditions?
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
        updateCouplingBCs_ = true;

        ptCoupling_->GetInputNodes(i, nodes);

        for ( UInt dof = 0; dof < ptCoupling_->GetInputDof(i); dof++ ) {
          for ( UInt j = 0; j < nodes->GetSize();
                j++, couplingBCsCounter_++) {

            eqnNr = eqnMap_->GetNodeEqn( (*nodes)[j], dof+1 );

            if (eqnNr==0) {
              Error( "The specified coupling node has no equation number",
                     __FILE__, __LINE__ );
            }

            algsys_->SetDirichlet( couplingBCsCounter_ + 1, pdeId_, eqnNr,
                                   help[dof+j*couplingDof], 1 );
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

  void SinglePDE::SetIDBC( const std::string &name,
                           const std::string &dofString, 
                           const std::string &value, 
                           const std::string &phase,
                           const std::string &dynamics) {
    ENTER_FCN("SinglePDE::SetIDBC");

    Warning( "Not working yet!", __FILE__, __LINE__ );
  //   // Try ro find existing entry in IDBC vector
//     Integer index = -1;
//     UInt dofspernode = results_[0]->dofNames.GetSize();

    
//     for ( UInt i = 0; i < bcs_id_.GetSize(); i++ ) {

//       if ( bcs_id_[i] == name ) {
//         if ( dofspernode > 1 ) {
//           if ( inhomDirichDof_[i] == dofString ) {
//             index = i;
//             break;
//           }
//         } else {
//           index = i;
//           break;
//         }
//       }
//     }
    
//     // If entry was already defined, we can always change
//     // the value of the boundary condition
//     if ( index != -1 ) {
//       val_id_[index] = value;
//       if (analysistype_ == HARMONIC
//           ||analysistype_ == MULTIHARMONIC) {
//         bcs_id_phase_[index] = phase;
//       }
//     } else {
//       // Entry was not yet defined. We can only add a new
//       // idbc entry, if eqnData was not yet initialized / created
//       if ( eqnMap_.get() == NULL) {
//         bcs_id_.Push_back( name );
//         val_id_.Push_back( value );
//         if ( dofspernode > 1 ) {
//           inhomDirichDof_.Push_back( dofString );
//         }
//         if (analysistype_ == HARMONIC
//             ||analysistype_ == MULTIHARMONIC) {
//           bcs_id_phase_.Push_back( phase );
//         } else {
//           fncnames_id_.Push_back( dynamics );
//         }

//         // Re-calculate the correct number of boundary condition entries
//         for ( UInt i = 0; i < bcs_id_.GetSize(); i++) {
//           numDirichletBCs_ += ptgrid_->GetNumNodes(bcs_id_[i]);
//         }
        
//       } else {
//         (*error) << "No new ihom. Dirichlet BC can be defined " 
//                  << "at this point, as the equation numbering "
//                  << "was already performed!";
//         Error( __FILE__, __LINE__ );
//       }
//     } 

  }

  
  // ======================================================
  // SCRIPTING SECTION
  // ======================================================

  void SinglePDE::RegisterFunctions() {
    typedef FctPointer<SinglePDE> FCPT;
    StdVector<ArgList> a;
    StdVector<FCPT*> pt;
    StdVector<std::string> name;

    // --- IDBC ---
    a.Push_back();
    a.Last().RegisterParam("name", ArgList::STRING );
    a.Last().RegisterParam("dof", ArgList::STRING );
    a.Last().RegisterParam("value", ArgList::STRING );
    a.Last().RegisterParam("phase", ArgList::STRING );
    a.Last().RegisterParam("dynamics", ArgList::STRING );
    pt.Push_back( new FCPT( this, &SinglePDE::Wrap_IDBC) );
    name.Push_back( "setIdbc" );
                  

    // --- GetValue ---
    a.Push_back();
    a.Last().RegisterParam( "name", ArgList::STRING );
    pt.Push_back( new FCPT( this, &SinglePDE::Wrap_GetValue) );
    name.Push_back( "getValue" );
    
    // Now register all functions with scripting 
    for (UInt i = 0; i < pt.GetSize(); i++ ) {
      Script_RegisterFct(name[i], pt[i], a[i] );
    }
  }
  
  void SinglePDE::Wrap_IDBC() {
    SCRIPT_GET( std::string, name );
    SCRIPT_GET( std::string, dof );    
    SCRIPT_GET( std::string, value );
    SCRIPT_GET( std::string, phase );
    SCRIPT_GET( std::string, dynamics );
    if ( dynamics == "" ) {
      dynamics = "none";
    }
    SetIDBC(name, dof, value, phase, dynamics );
  }
  
  void SinglePDE::Wrap_GetValue() {
    SCRIPT_GET( std::string, name );
    StdVector<UInt> nodeNrs;
    ptgrid_->GetNodesByName( nodeNrs, name );
    for (UInt i=0; i<nodeNrs.GetSize(); i++) {
      Double val;
      sol_->Get(nodeNrs[i]-1,0,val);
      SCRIPT_RETVAL.Push_back(GenStr(val));
    }
  }
  
} // end of namespace
