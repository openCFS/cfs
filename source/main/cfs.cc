// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim:fenc=utf-8:ft=tcl:et:sw=2:ts=2:sts=2
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#include <iomanip>

#include <def_use_xerces.hh>
#include <def_use_mpcci.hh>
#include <boost/version.hpp>

#include "main/cfs.hh"
#include "Utils/Timer.hh"
#include "DataInOut/DefineFiles/definefiles.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/MaterialHandler.hh"
#include "DataInOut/programOptions.hh"
#include "Domain/domain.hh"
#include "Domain/entityList.hh"
#include "General/environment.hh"
#include "DataInOut/ParamHandling/SkeletonConf.hh"
#include "DataInOut/ParamHandling/InfoNode.hh"
#include "DataInOut/ParamHandling/Xerces.hh"
#include "DataInOut/resultHandler.hh"
#include "DataInOut/coloredConsole.hh"
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
using namespace std;

int main(int argc, const char **argv)
{
  CFS* cfs = new CFS(argc, argv);
  int ret = cfs->Run();
  
  delete cfs; // so that we can delete
  return ret;
}

CFS::CFS(int argc, const char **argv)
{

  resultHandler = NULL;
  materialHandler = NULL;

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


  // Create object for logging information
  Info = new WriteInfo();

  // =========================================================================
  // HANDLE COMMAND LINE PARAMETERS
  // =========================================================================
  progOpts = new ProgramOptions(argc, argv);

  // Parse command line
  progOpts->ParseData();

  // Log program startup
  progOpts->GetHeaderString( cout );
  
  // Get information about exception handling
  Exception::segfault_ = progOpts->GetForceSegFault();

  // Set global Enums, the rest is set by the classes
  SetGlobalEnums();

  // the new xml logging derived from the ParamNode
  info = new InfoNode(progOpts->GetSimName() + ".info.xml", "<?xml version=\"1.0\"?>");
  info->SetName("cfsInfo");
  info->Get("status")->SetValue("running"); // to be overwritten by "aborted" or "finished"
  timer = info->Get(InfoNode::SUMMARY)->SetValue(new Timer());
  timer->Start(); // ignore that this is not the real beginning

  // Print information about program start time and host
  using namespace boost::posix_time;
  using namespace boost::gregorian;

  // our calculation environment
  InfoNode* env = info->Get(InfoNode::HEADER)->Get("environment");
  start_time_ = to_simple_string( second_clock::local_time() );
  env->Get("started")->SetValue(start_time_);

  char host[256];
  hostname_ = gethostname( host, 256 ) > 0 ? host : "";

  if(!hostname_.empty())
    env->Get("host")->SetValue(hostname_);
  


  // =========================================================================
  // GENERATE CENTRAL MESSENGER OBJECT
  // =========================================================================
#ifdef USE_SCRIPTING

  // Try to determine (optional) script file name
  string scriptFileName = progOpts->GetScriptFileStr();

  messenger = fileHandler.CreateScriptMessenger( scriptFileName );

  // Check if script file was provided
  if ( scriptFileName != "" ) {
    cout << "++ Script file: '" << scriptFileName << "'" << endl;

    // Create new central messenger object (up to now only tcl available)
    messenger->ReadScriptFile(scriptFileName );

    // Call initialization procedure
    StdVector<string> context;
    context.Push_back( progOpts->GetSimName() );
    messenger->TriggerEvent(CFSMessenger::CFS_Init, context);
  }

#endif

  // =========================================================================
  // ACTIVATE DEBUGGING STUFF
  // =========================================================================

#ifdef MEMTRACE
  fileHandler.OpenFile( MEMTRACE_FILE );
#endif

  // =========================================================================
  // SETUP OF MPCCI
  // =========================================================================

#ifdef MpCCI
  cout << "++ Setting up MpCCI interface" << endl;

  CCI_Init_with_id_string( &argc, const_cast<char***>(&argv), "simulationcode2" );
  //  CCI_Init( &argc, const_cast<char***>(&argv));
#endif

}


CFS::~CFS()
{
#ifdef MpCCI
  CCI_Finalize();
#endif

#ifdef USE_SCRIPTING
  delete messenger;
  messenger = NULL;
#endif

  // flush last information.
  info->ToFile();
  delete info;

  // Delete objects
  delete param;
  delete domain; // might write ersatz material file if <export save="finally"/> in optimization
  delete progOpts;
  delete resultHandler;

  // Delete global string streams
  delete Info;
}

int CFS::Run()
{
  try
  {
    if(progOpts->GetWriteSkeleton())
    {
      WriteXMLSkeleton();
      return 0;
    }

    if(!progOpts->IsQuiet())
    {
      cout << "Simulation run started at " << start_time_;
      if(!hostname_.empty()) cout << " on " << hostname_;
      cout << endl;
    }


    ReadXMLFile();
    SetupIO();

    domain = new Domain( gridInputs, resultHandler, materialHandler );
    // Create grid
    domain->CreateGrid();

    if(progOpts->GetPrintGrid())
      PrintGrid();
    else
      SolveProblem();

    timer->Stop();
    if(!progOpts->IsQuiet()) cout << endl; // conditional empty line
    cout << ">> Total time: wall clock: '" << timer->GetWallTime()
              << "s' CPU time: '" << timer->GetCPUTime() << "s'" << endl << endl;

    // write the info object
    info->Get("status")->SetValue("finished"); // overwrite 'running'
    info->Get(InfoNode::SUMMARY)->Get("memory/final")->SetValue(MemoryUsage(false));
    info->Get(InfoNode::SUMMARY)->Get("memory/peak")->SetValue(MemoryUsage(true));

    return 0;
  }
  catch(exception& ex)
  {
    if(!progOpts->IsQuiet())
    {
      cerr << endl << endl
                << "***********************************************************************"
                << endl << fg_red << " SIMULATION RUN FAILED!  -  CAUGHT EXCEPTION:" << fg_reset
                << endl << endl
                << ex.what() << endl << endl
                << "***********************************************************************"
                << endl << endl;
    }
    else
    {
      cerr << endl << ">> Error: " << ex.what() << endl;
    }

    // Print error cause to info file
    if(info != NULL)
    {
      InfoNode* errorNode = info->Get(InfoNode::ERROR);
      errorNode->SetValue(ex.what());
      info->Get("status")->SetValue("aborted");
      info->ToFile();
    }

    return 1;
  }
}

void CFS::SetGlobalEnums()
{
  SetEnvironmentEnums();
  BasePDE::SetEnums();
  EntityList::SetEnums();
}

void CFS::PrintGrid()
{
  cout << "Printing grid to file " << endl << endl;
  domain->PrintGrid();
}

void CFS::SolveProblem()
{
  // Set up Problem
 domain->PostInit();

 // Solves the driver or optimization problem
 domain->SolveProblem();

 cout << "++ Finished solving the problem" << endl;

#ifdef USE_SCRIPTING
   // Call intialization procedure
   StdVector<string> context;
   context.Push_back( progOpts->GetSimName() );
   messenger->TriggerEvent(CFSMessenger::CFS_Finish, context);
#endif
}

void CFS::Process()
{
  ReadXMLFile();
  SetupIO();

  domain = new Domain( gridInputs, resultHandler, materialHandler );
  // Create grid
  domain->CreateGrid();
}



void CFS::WriteXMLSkeleton()
{
  // This block is used for writing a skeleton XML-File to make life easier
  // for the user. This also will import some information from the mesh-file.
  // Since normally all opening of files is handled in definefiles
  // the debug and trace file are opened separately here. Definefiles cannot
  // be used, since it would try to open other files also, that need not be
  // present????
  // Hardcoded: Read mesh file
  string meshFile = progOpts->GetMeshFileStr();
  string simName = progOpts->GetSimName();
  if(meshFile == "")
    meshFile = simName + ".mesh";
  SimInput * ptInputfile = new SimInputMESH(meshFile, NULL);
  ptInputfile->InitModule();
  // class writing log-information
  SkeletonConf *ptskel = new SkeletonConf(ptInputfile);
  ptskel->WriteConf();
  delete ptskel;
  delete ptInputfile;
}


void CFS::ReadXMLFile()
{
  // Generate parameter handler and pass address to global pointer
  string xmlFile = progOpts->GetParamFileStr();

  // Write information to command line
  cout << "++ Reading parameter file '" + xmlFile + "'" << endl;

#ifndef USE_XERCES
  EXCEPTION( "I am sorry to say, but CFS only can be compiled with XERCES-support");
#endif

  // this is the new param staff which replaces the old params - delete this comment finally
  param = NULL;

  string schema = progOpts->GetSchemaPathStr();
  schema += "/CFS-Simulation/CFS.xsd";

  // Initialize our xerces dom parser to handle the cfs xml file
  Xerces* xerces = new Xerces(xmlFile, schema);

  // set the global ParamNode tree pointer
  param = xerces->CreateParamNodeInstance();
  // save us in the info stuff, with defaults but no comments
  // release the xerces ressources, param is not affected
  delete xerces;
}

void CFS::SetupIO()
{
  // read type of mesh-libraray
//  string libmesh = "mesh";
  map<string, shared_ptr<SimInput> > inFiles;

//  ParamNode * meshNode = param->Get("input", false );
//  if( meshNode )
//    meshNode->Get("meshLibrary", libmesh, false);

//  if ( libmesh != "structGrid" ) {
    // Generate mesh reader
  fileHandler.CreateSimInputFiles( inFiles, gridInputs );
//  }

  // generate material handler
  materialHandler = fileHandler.CreateMaterialHandler();

  // Open file for status reports by CFS++
  Info->CreateFile();

  // Create simulation output writer
  map<string, shared_ptr<SimOutput> > outFiles;
  fileHandler.CreateSimOutputFiles( outFiles );

  // Create resulthandler and pass the output files
  resultHandler = new ResultHandler( ResultHandler::EMBEDDED );
  map<string, shared_ptr<SimOutput> >::iterator outputIt;
  map<string, shared_ptr<SimInput> >::iterator inputIt;
  outputIt = outFiles.begin();
  inputIt = inFiles.begin();
  for( ; outputIt != outFiles.end(); outputIt++ ) {
    resultHandler->AddOutputDest( outputIt->second, outputIt->first );
  }
  for( ; inputIt != inFiles.end(); inputIt++ ) {
    resultHandler->AddInputReader( inputIt->second, inputIt->first );
  }


  // Log command line parameters
  progOpts->ToInfo(info->Get(InfoNode::HEADER)->Get("progOpts"));
  // log the optinal id/name/token/label from <cfsSimulation id="..">
  info->Get(InfoNode::HEADER)->Get("id")->SetValue(param->Get("id"));

  // Open file for status reports by OLAS
  fileHandler.OpenFile( OLAS_FILE );
}

