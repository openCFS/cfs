// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iostream>
#include <iomanip>

#include <def_use_xerces.hh>
#include <def_use_mpcci.hh>

#include "Utils/myclock.hh"
#include "DataInOut/DefineFiles/definefiles.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/MaterialHandler.hh"
#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "DataInOut/CommandLine/CommandLineHandlerSetting.hh"
#include "Domain/domain.hh"
#include "General/environment.hh"
#include "DataInOut/ParamHandling/SkeletonConf.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/Xerces.hh"
#include "DataInOut/resultHandler.hh"
#include "Utils/profiler.hh"

#include <def_use_mesh.hh>

#ifdef MpCCI

#if (MpCCI_RELEASE == 305)
#include <mpcci.h>
#else
#include <cci.h>
#endif

#endif

#ifdef GRIDLIB
#include "Domain/Gridlib/interface_gridlib.hh"
#endif

#ifdef USE_MESH
#include "DataInOut/SimInOut/AnsysFile/simInputMESH.hh"
#endif

#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh"
#endif


using namespace CoupledField;

#ifdef PARALLEL
#define STDOUT if (!commrank) std::cout
#else
#define STDOUT std::cout
#endif

int main( int argc, const char **argv ) {


  // Set segfault to false
  Exception::segfault_ = false;
  try 
  {
        
  // =========================================================================
  // TIMING
  // =========================================================================
  MyClock oClockTotal;
  oClockTotal.Start();

  // =========================================================================
  // CREATE INFORMATION OBJECT
  // =========================================================================

  // Create object for logging information
  Info = new WriteInfo();

  // Log program startup
  Info->PrintHeader();


  // =========================================================================
  // HANDLE COMMAND LINE PARAMETERS
  // =========================================================================
  commandLine = new CommandLineHandlerSetting( argc, argv );

# ifdef DEBUG
  // Get information about exception handling
  Exception::segfault_ = commandLine->GetForceSegFault();
#endif

  // =========================================================================
  // GENERATE OBJECT FOR HANDLING FILE-IO
  // =========================================================================
  Info->StartProgress( "Generating object for handling file-IO" );
  DefineInOutFiles FileHandler( commandLine->GetSimName().c_str() );
  Info->FinishProgress();
  

  // =========================================================================
  // GENERATE CENTRAL MESSENGER OBJECT
  // =========================================================================
#ifdef USE_SCRIPTING

  // Try to determine (optional) script file name
  std::string scriptFileName = commandLine->GetScriptFileName();
  
  messenger = FileHandler.CreateScriptMessenger( scriptFileName );
  
  // Check if script file was provided
  if ( scriptFileName != "" ) {
    std::stringstream msg;
    msg << "Activating Scripting using '" << scriptFileName << "'";
    Info->StartProgress( msg.str() );

    // Create new central messenger object (up to now only tcl available)
    messenger->ReadScriptFile(scriptFileName );
    Info->FinishProgress();    
  }

#endif


  // =========================================================================
  // ACTIVATE DEBUGGING STUFF
  // =========================================================================


  // Open file for debugging ouput
#ifdef DEBUG
  Info->StartProgress( "Opening file for writing debugging output" );
  FileHandler.OpenFile( DEBUG_FILE );
  Info->FinishProgress();
#endif

#ifdef MEMTRACE
  Info->StartProgress( "Opening file for tracing memory allocation" );
  FileHandler.OpenFile( MEMTRACE_FILE );
  Info->FinishProgress();
#endif

#ifdef PROFILING
  if ( commandLine->GetDoProfile() == true ) {
    Info->StartProgress( "Opening file for profiling output" );
  }
  else {
    Info->StartProgress( "Skipping generation of profiling output" );
  }
  profiler = new Profiler();
  SETPROFILE( "Begin of program" );
  Info->FinishProgress();
#endif

  // =========================================================================
  // SETUP OF MPCCI
  // =========================================================================

#ifdef MpCCI
  Info->StartProgress( "Setting up MpCCI interface" );

  CCI_Init_with_id_string( &argc, const_cast<char***>(&argv), "simulationcode2" );  
  //  CCI_Init( &argc, const_cast<char***>(&argv));

  Info->FinishProgress();
#endif

 // =========================================================================
  // WRITE SKELETON XML-FILE
  // =========================================================================

  // This block is used for writing a skeleton XML-File to make life easier
  // for the user. This also will import some information from the mesh-file.
  // Since normally all opening of files is handled in definefiles
  // the debug and trace file are opened separately here. Definefiles cannot
  // be used, since it would try to open other files also, that need not be
  // present????
  if ( commandLine->GetWriteSkeleton() == true ) {
 
#ifdef DEBUG
    FileHandler.OpenFile( DEBUG_FILE );
#endif

    // Hardcoded: Read mesh file
    std::string meshFile = commandLine->GetMeshFile();
    std::string simName = commandLine->GetSimName();
    if(meshFile == "")
      meshFile = simName + ".mesh";
    SimInput * ptInputfile = new SimInputMESH(meshFile, NULL);
    ptInputfile->InitModule();
    // class writing log-information
    SkeletonConf *ptskel = new SkeletonConf(ptInputfile);
    ptskel->WriteConf();
    delete ptskel;
    delete ptInputfile;
    delete Info;
    
    return 0;
  }


  // =========================================================================
  // HANDLE XML-FILE
  // =========================================================================

  // Generate parameter handler and pass address to global pointer
  std::string xmlFile = commandLine->GetParamFile();

  // Write information to command line
  std::stringstream msg;
  msg << "Reading parameters from file '" << xmlFile << "'";
  Info->StartProgress( msg.str() );
  
#ifdef USE_XERCES
  // this is the new param staff which replaces the old params - delete this comment finally
  param = NULL;           

  std::string schema = commandLine->GetSchemaPath();
  schema += "/CFS-Simulation/CFS.xsd";

  // Initialize our xerces dom parser to handle the cfs xml file
  Xerces* xerces = new Xerces(xmlFile, schema);

  // set the global ParamNode tree pointer
  param = xerces->CreateParamNodeInstance();

  // release the xerces ressources, param is not affected
  delete xerces;
  
#else
  EXCEPTION( "I am sorry to say, but CFS only can be compiled with XERCES-support");
#endif

  Info->FinishProgress();


  // =========================================================================
  // SETUP OF IO-STUFF
  // =========================================================================

  // read type of mesh-libraray
  std::string libmesh = "mesh";
  std::map<std::string, shared_ptr<SimInput> > inFiles;
  std::map<std::string, StdVector<shared_ptr<SimInput> > > gridInputs;
  ParamNode * meshNode = param->Get("input", false );
  if( meshNode )
    meshNode->Get("meshLibrary")->AsString();

  if ( libmesh != "structGrid" ) {
    // Generate mesh reader
    Info->StartProgress( "Generating mesh reader" );
    FileHandler.CreateSimInputFiles( inFiles, gridInputs );
    Info->FinishProgress();
  }

 

  // generate material handler
  MaterialHandler * ptMatHandler;
  Info->StartProgress( "Generating material reader");
  ptMatHandler = FileHandler.CreateMaterialHandler();
  Info->FinishProgress();
  
  Info->StartProgress( "Generating remaining output files" );

  // Open file for status reports by CFS++
  Info->CreateFile();

  // Create simulation output writer
  std::map<std::string, shared_ptr<SimOutput> > outFiles;
  FileHandler.CreateSimOutputFiles( outFiles );
  
  // Create resulthandler and pass the output files
  ResultHandler * ptHandler = 
    new ResultHandler( ResultHandler::EMBEDDED );
  std::map<std::string, shared_ptr<SimOutput> >::iterator outputIt;
  std::map<std::string, shared_ptr<SimInput> >::iterator inputIt;
  outputIt = outFiles.begin();
  inputIt = inFiles.begin();
  for( ; outputIt != outFiles.end(); outputIt++ ) {
    ptHandler->AddOutputDest( outputIt->second, outputIt->first );
  }
  for( ; inputIt != inFiles.end(); inputIt++ ) {
    ptHandler->AddInputReader( inputIt->second, inputIt->first );
  }
  
  
  // Log command line parameters
  std::ostream *myInfo = Info->GetInfoStreamPointer();
  commandLine->PrintParams( *myInfo, false );

  // Open file for status reports by OLAS
  FileHandler.OpenFile( OLAS_FILE );

  Info->FinishProgress();

  // =========================================================================
  // GENERATION OF DOMAIN OBJECT
  // =========================================================================
  SETPROFILE("Before Creation of Domain");
  domain = new  Domain( gridInputs, ptHandler, ptMatHandler );
  SETPROFILE("After Creation of Domain");

  // Create grid
  domain->CreateGrid();


  // =========================================================================
  // Only output of grid
  // =========================================================================
  if ( commandLine->GetPrintGrid() == true ) {
    STDOUT << "Printing grid to file " << myEndl << myEndl;
    domain->PrintGrid();
    delete domain;
    delete Info;
    delete commandLine; 
    delete ptHandler;
    return 0;
  }

   // Set up Problem
  domain->PostInit();
    
  // Solves the driver or optimization problem
  domain->SolveProblem();
    
  std::cout << std::endl;
  Info->StartProgress( "Finished solving the problem" );
  Info->FinishProgress();

  std::cout << std::endl;
  std::cout << std::endl;

#ifdef USE_SCRIPTING

    // Call intialization procedure
    StdVector<std::string> context;
    context.Push_back( commandLine->GetSimName() );
    messenger->TriggerEvent(CFSMessenger::CFS_Finish, context);
#endif

  Double wTime, cTime;
  oClockTotal.GetTime(wTime, cTime);
  std::cout << std::setw(70) << std::setfill('=') << "" << std::endl;
  std::cout << "TOTAL TIME:\t\t Wall clock: " << wTime
            << " s\t CPU: " << cTime << " s" << std::endl;
  std::cout << std::setw(70) << std::setfill('=') << "" 
            << std::endl << std::endl;


  // =========================================================================
  // CLEANUP PHASE
  // =========================================================================

#ifdef MpCCI
  CCI_Finalize();
#endif

  // Delete objects
  delete domain;
  delete Info;
  delete commandLine;
  delete ptHandler;

#ifdef USE_SCRIPTING
  delete messenger;
  messenger = NULL;
#endif

#ifdef PROFILING
  delete profiler;
  profiler = NULL;
#endif
  
  }  
  catch(std::exception& ex) 
  {
    Error(ex.what(), __FILE__, __LINE__);
  }


  // Delete global string streams
  delete error;
  delete warning;

  // Seems that everything went fine
  return 0;
}
