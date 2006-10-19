#include <fstream>
#include <iostream>
#include <cstdio>

#include <string>

#include "definefiles.hh"
#include "DataInOut/AnsysFile/ansysfile.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "DataInOut/cfsRST/outCFS.hh"

#ifdef GSI
#include "DataInOut/GSI/outGSI.hh"
#endif

#ifdef GIDPOST
#include "DataInOut/GiD/outGiD.hh"
#endif

#include "DataInOut/PlainMaterialHandler.hh"
#include "DataInOut/XMLMaterialHandler.hh"

#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/ParamHandling/PlainXMLParamHandler.hh"
#include "DataInOut/ParamHandling/XMLParamHandler.hh"

#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh"
#endif

#ifdef USE_SCRIPTING_TCL
#include "DataInOut/Scripting/tcl/tcl-messenger.hh"
#endif

#ifdef USE_SCRIPTING_PY
#include "DataInOut/Scripting/python/py-messenger.hh"
#endif

#ifdef USE_DATABASE
#include "DataInOut/Database/outDB.hh"
#endif

// Maximal length of the trailing postfix of an auxilliary name,
// i.e. the length of the extension after the basename
#define MAXPOSTFIX 15

#ifdef PARALLEL
#define MAXRANK 4
#else
#define MAXRANK 0
#endif

#ifdef PARALLEL
#include <mpi.h>
#endif

namespace CoupledField
{

#ifdef MEMTRACE
  Double sumdmem;
  Double sumimem;
#endif




  // ===============
  //   Constructor
  // ===============
  DefineInOutFiles::DefineInOutFiles( const Char *name ) {

#ifdef PARALLEL
    // find out who I am and write to seperate files:
    int commrank;
    MPI_Comm comm = MPI_COMM_WORLD;
    MPI_Comm_rank(comm,&commrank);
    char *rank = new char[5];
    sprintf(rank,"_%d",commrank);
#endif

#ifdef PARALLEL
    // strcat( basename, rank );
#endif

    // Initialise internal pointers
    ptWriteResults_ = NULL;
    infileType_ = NULL;
    ptMaterialHandler_ = NULL;

  }
 

  // ==============
  //   Destructor
  // ==============
  DefineInOutFiles::~DefineInOutFiles() {

    ENTER_FCN( "DefineInOutFiles::~DefineInOutFiles" );

#ifdef DEBUG
    delete debug;
    debug = NULL;
#endif

#ifdef MEMTRACE
    delete memtrace;
    memtrace = NULL;
#endif

    // Delete pointer to OLAS report file
    delete cla;
    cla = NULL;

    // Delete internal pointers
    delete ptWriteResults_;
    ptWriteResults_ = NULL;

    delete infileType_;
    infileType_ = NULL;

    delete ptMaterialHandler_;
    ptMaterialHandler_ = NULL;

  }


  // ==============================
  //   Generate mesh file pointer
  // ==============================
  FileType* DefineInOutFiles::CreateMeshFileHandler() {

    ENTER_FCN( "DefineInOutFiles::CreateMeshFileHandler" );

    std::string meshFile = commandLine->GetMeshFile();
    if ( params->HasValue( "format", "mesh", "input" ) ) {
      infileType_ = new AnsysFile( meshFile.c_str() );
    }
    else {
      Error( "Wrong format for input file. Please, check your data!",
             __FILE__, __LINE__ );
    }

    return infileType_;
  }


  // ================================
  //   Generate output file pointer
  // ================================
  WriteResults*
  DefineInOutFiles::Create_ptWriteResults()
  {
    ENTER_FCN( "DefineInOutFiles::Create_ptWriteResults" );

    std::string simName = commandLine->GetSimName();

    std::string outformat;
    params->Get( "format", outformat, "output" );

    if ( outformat == "gmv" ) {
      ptWriteResults_= new WriteResultsGMV( simName.c_str() );
    }
    else if ( outformat == "unv" ) 
      ptWriteResults_= new WriteResultsUnverg( simName.c_str() ); 
    else if ( outformat == "cfsRST" ) 
      ptWriteResults_= new WriteResultsCFS( simName.c_str() ); 
#ifdef GSI
    else if ( outformat == "gsi" ) 
      ptWriteResults_= new WriteResultsGSI( simName.c_str() );
#endif

#ifdef GIDPOST
    else if ( outformat == "gid" ) 
      ptWriteResults_= new WriteResultsGiD( simName.c_str() );
#endif

#ifdef USE_DATABASE
    else if (outformat == "database") 
      ptWriteResults_= new WriteResultsDatabase( simName.c_str() );
#endif
    else {
      (*error) << "Output format '" << outformat << "' not recognised";
      Error( __FILE__, __LINE__ );
    }

    if (!ptWriteResults_) 
      Error("Can't open file for output results",__FILE__,__LINE__);
  
    return ptWriteResults_;
  }


  // ================================
  //   Generate output file pointer
  // ================================
  MaterialHandler *
  DefineInOutFiles::CreateMaterialHandler() {
    ENTER_FCN( "DefineInOutFiles::CreateMaterialHandler" );

    std::string fileName, format;    

    // Determine filename and format
    params->Get( "file", fileName, "materialData" );
    params->Get( "format", format, "materialData" );

    if ( format == "dat" ) {
      ptMaterialHandler_ = new PlainMaterialHandler( fileName );
    } else if ( format == "xml" ) {
      ptMaterialHandler_ = new XMLMaterialHandler( fileName );
    } else {
      (*error) << "CreateMaterialHandler: Format '" << format
               << "' is not recognized!";
      Error( __FILE__, __LINE__ );
    }
    return ptMaterialHandler_;

  }


  // ************
  //   OpenFile
  // ************
  void DefineInOutFiles::OpenFile( AuxFileType fileType ) {

    ENTER_FCN( "DefineInOutFiles::OpenFile" );

    std::string fileName;

    switch( fileType ) {


      // -------------
      //  TRACE_FILE
      // -------------
    case TRACE_FILE:

      fileName = commandLine->GetSimName() + ".trace";

#ifdef TRACE
      trace = new std::ofstream( fileName.c_str() );
      if ( trace == NULL ) {
        (*error) << "Failed to open file '" << fileName << "' for "
                 << "writing function trace information!";
        Error( __FILE__, __LINE__ );
      }
      break;
#else
      (*error) << "Executable was compiled without support for function "
               << "tracing! Refusing to open trace file '"
               << fileName << "'";
      Error( __FILE__, __LINE__ );
      break;
#endif


      // -------------
      //  DEBUG_FILE
      // -------------
    case DEBUG_FILE:

      fileName = commandLine->GetSimName() + ".deb";

#ifdef DEBUG
      debug = new std::ofstream( fileName.c_str() );
      if ( debug == NULL ) {
        (*error) << "Failed to open file '" << fileName << "' for "
                 << "writing debug information!";
        Error( __FILE__, __LINE__ );
      }
      break;
#else
      (*error) << "Executable was compiled without debug support!"
               << "Refusing to open debug file '" << fileName << "'";
      Error( __FILE__, __LINE__ );
      break;
#endif


      // ---------------
      //  MEMTRACE_FILE
      // ---------------
    case MEMTRACE_FILE:

      fileName = commandLine->GetSimName() + ".mem";

#ifdef MEMTRACE
      memtrace = new std::ofstream( fileName.c_str() );
      if ( memtrace == NULL ) {
        (*error) << "Failed to open file '" << fileName << "' for "
                 << "writing memory trace information!";
        Error( __FILE__, __LINE__ );
      }
      break;
#else
      (*error) << "Executable was compiled without support for memory "
               << "tracing! Refusing to open trace file '"
               << fileName << "'";
      Error( __FILE__, __LINE__ );
      break;
#endif


      // ------------
      //  OLAS_FILE
      // ------------
    case OLAS_FILE:

      fileName = commandLine->GetSimName() + ".las";
      cla = new std::ofstream( fileName.c_str() );
      if ( cla == NULL ) {
        (*error) << "Failed to open file '" << fileName << "' for "
                 << "writing status reports of OLAS, the linear algebra "
                 << "sub-system!";
        Error( __FILE__, __LINE__ );
      }
      break;

    }
  }

  CFSMessenger* DefineInOutFiles::CreateScriptMessenger( const std::string& fileName) {
    ENTER_FCN( "DefineInOutFiles::CreateScriptMessenger" );
    
    // check filename, if it is not empty
    if( fileName == "") {
      return new CFSMessenger();
    }

#ifdef USE_SCRIPTING_TCL
    else if( fileName.find( ".tcl") != std::string::npos ) {
      return new TCL_CFSMessenger();
    }
#endif
#ifdef USE_SCRIPTING_PY
    else if( fileName.find( ".py") != std::string::npos ) {
      return new PY_CFSMessenger();
    }
#endif
    else {
      std::stringstream msg;
      msg << "Could not determine script language of file '"
          << fileName << "'!";
      Error( msg.str().c_str(), __FILE__, __LINE__ );
    }
  }

} // end of namespace
