// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim:fenc=utf-8:ft=tcl:et:sw=2:ts=2:sts=2
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <iomanip>

#include <def_use_xerces.hh>
#include <def_use_mpcci.hh>
#include <boost/version.hpp>

#include "Utils/myclock.hh"
#include "DataInOut/DefineFiles/definefiles.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/MaterialHandler.hh"
#include "DataInOut/programOptions.hh"
#include "Domain/domain.hh"
#include "General/environment.hh"
#include "DataInOut/ParamHandling/SkeletonConf.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"
#include "DataInOut/ParamHandling/Xerces.hh"
#include "DataInOut/resultHandler.hh"
#include "DataInOut/coloredConsole.hh"
#include "Utils/profiler.hh"
#include <unistd.h>
#include <boost/date_time/posix_time/posix_time.hpp>

#include <def_use_mesh.hh>

#ifdef MpCCI

#if (MpCCI_RELEASE == 305)
#include <mpcci.h>
#else
#include <cci.h>
#endif

#endif

#ifdef USE_MESH
#include "DataInOut/SimInOut/AnsysFile/simInputMESH.hh"
#endif

#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh"
#endif


using namespace CoupledField;

int main( int argc, const char **argv ) {


  // Check if the boost version is >= 1.34, as the boost::filesystem::path
  // interface has changed heavily.
  // In version < 1.34 there is a name checking mechanism activated, which
  // doest not allow file/directory names containing characters like '+', ..
  // Therefore we need to replace this standard name checker.
#if ( BOOST_VERSION < 103400)
  boost::filesystem::path::default_name_check( boost::filesystem::native);
#endif

  // Set segfault to false
  Exception::segfault_ = false;
  try
  {

  BasePDE::SetEnums();

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


  // =========================================================================
  // HANDLE COMMAND LINE PARAMETERS
  // =========================================================================

  progOpts = new ProgramOptions(argc, argv);

  // Parse command line
  progOpts->ParseData();

  // Log program startup
  progOpts->GetHeaderString( std::cout );
  
  // Get information about exception handling
  Exception::segfault_ = progOpts->GetForceSegFault();

  // the new xml logging derived from the ParamNode
  info = new InfoNode(progOpts->GetSimName() + ".info.xml", "<?xml version=\"1.0\"?>");
  info->SetName("cfsInfo");
  info->Get("status")->SetValue("running"); // to be overwritten by "aborted" or "finished"

  // GENERATE OBJECT FOR HANDLING FILE-IO
  DefineInOutFiles FileHandler;

  // Print information about program start time and host
  using namespace boost::posix_time;
  using namespace boost::gregorian;
  std::string now = to_simple_string( second_clock::local_time() );
  char host[256];
  int ret = gethostname( host, 256 );

  // our calculation environment
  info->Get("calculation")->Get("environment")->Get("started")->SetValue(now);
  if(ret == 0)
    info->Get("calculation")->Get("environment")->Get("host")->SetValue(host);
  
  if(!progOpts->IsQuiet())
  {
    std::string out = "Simulation run started at " + now;
    if (ret==0)
      out += " on " + std::string( host ) + "\n";
    else
      out += "\n";
    std::cout << out << std::endl;
    Info->PrintF( "", out.c_str() );
  }


  // =========================================================================
  // GENERATE CENTRAL MESSENGER OBJECT
  // =========================================================================
#ifdef USE_SCRIPTING

  // Try to determine (optional) script file name
  std::string scriptFileName = progOpts->GetScriptFileStr();

  messenger = FileHandler.CreateScriptMessenger( scriptFileName );

  // Check if script file was provided
  if ( scriptFileName != "" ) {
    std::stringstream msg;
    msg << "Script file: '" << scriptFileName << "'";
    Info->StartProgress( msg.str() );

    // Create new central messenger object (up to now only tcl available)
    messenger->ReadScriptFile(scriptFileName );

    // Call initialization procedure
    StdVector<std::string> context;
    context.Push_back( progOpts->GetSimName() );
    messenger->TriggerEvent(CFSMessenger::CFS_Init, context);

    Info->FinishProgress();
  }

#endif


  // =========================================================================
  // ACTIVATE DEBUGGING STUFF
  // =========================================================================

#ifdef MEMTRACE
  Info->StartProgress( "Opening file for tracing memory allocation" );
  FileHandler.OpenFile( MEMTRACE_FILE );
  Info->FinishProgress();
#endif

#ifdef PROFILING
  if ( progOpts->GetDoProfile() == true ) {
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
  if ( progOpts->GetWriteSkeleton() == true )
  {
    // Hardcoded: Read mesh file
    std::string meshFile = progOpts->GetMeshFileStr();
    std::string simName = progOpts->GetSimName();
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
  std::string xmlFile = progOpts->GetParamFileStr();

  // Write information to command line
  std::stringstream msg;
  msg << "Reading parameter file '" << xmlFile << "'";
  Info->StartProgress( msg.str() );

#ifndef USE_XERCES
  EXCEPTION( "I am sorry to say, but CFS only can be compiled with XERCES-support");
#endif

  // this is the new param staff which replaces the old params - delete this comment finally
  param = NULL;

  std::string schema = progOpts->GetSchemaPathStr();
  schema += "/CFS-Simulation/CFS.xsd";

  // Initialize our xerces dom parser to handle the cfs xml file
  Xerces* xerces = new Xerces(xmlFile, schema);

  // set the global ParamNode tree pointer
  param = xerces->CreateParamNodeInstance();
  // save us in the info stuff, with defaults but no comments
  // todo: info->Get(InfoNode::HEADER)->Get("cfsSimulation")->SetValue(param);
  // release the xerces ressources, param is not affected
  delete xerces;

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
    FileHandler.CreateSimInputFiles( inFiles, gridInputs );
  }



  // generate material handler
  MaterialHandler * ptMatHandler;
  ptMatHandler = FileHandler.CreateMaterialHandler();

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
  progOpts->ToInfo(info->Get(InfoNode::HEADER)->Get("progOpts"));

  // Open file for status reports by OLAS
  FileHandler.OpenFile( OLAS_FILE );

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
  if ( progOpts->GetPrintGrid() == true ) {
    std::cout << "Printing grid to file " << std::endl << std::endl;
    domain->PrintGrid();
    delete domain;
    delete Info;
    delete progOpts;
    delete ptHandler;
    return 0;
  }

   // Set up Problem
  domain->PostInit();

  // Solves the driver or optimization problem
  domain->SolveProblem();

  Info->StartProgress( "Finished solving the problem" );
  Info->FinishProgress();

#ifdef USE_SCRIPTING

    // Call intialization procedure
    StdVector<std::string> context;
    context.Push_back( progOpts->GetSimName() );
    messenger->TriggerEvent(CFSMessenger::CFS_Finish, context);
#endif

  Double wTime, cTime;
  oClockTotal.GetTime(wTime, cTime);
  std::cout << std::endl
            << ">> Total time: wall clock: '" << wTime 
            << "s' CPU time: '" << cTime << "s'\n";
  InfoNode* in = info->Get(InfoNode::SUMMARY)->Get("runTime");
  in->Get("wall")->SetValue(wTime);
  in->Get("CPU")->SetValue(cTime);
  in->SetComment("in seconds");

  // =========================================================================
  // CLEANUP PHASE
  // =========================================================================

#ifdef MpCCI
  CCI_Finalize();
#endif

  // write the info object
  info->Get("status")->SetValue("finished"); // overwrite 'running'
  info->ToFile();
  delete info;
  info = NULL;

  // Delete objects
  delete param;
  delete domain;
  delete progOpts;
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
    if(!progOpts->IsQuiet())
    {
      std::cerr << std::endl << std::endl
                << "***********************************************************************"
                << std::endl << fg_red << " SIMULATION RUN FAILED!  -  CAUGHT EXCEPTION:" << fg_reset
                << std::endl << std::endl 
                << ex.what() << std::endl << std::endl
                << "***********************************************************************"
                << std::endl << std::endl;
    }
    else
    {
      std::cerr << std::endl << ">> Error: " << ex.what() << std::endl;
    }
    
    // Print error cause to info file
    if(info != NULL)
    {
      InfoNode* errorNode = info->Get(InfoNode::ERROR);
      errorNode->SetValue(ex.what());
      info->Get("status")->SetValue("aborted");
      info->ToFile();
    }
  }


  // Delete global string streams
  delete warning;
  delete Info; // might be used in catch!

  // Seems that everything went fine
  return 0;
}
