// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim:fenc=utf-8:ft=cpp:et:sw=2:ts=2:sts=2
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#include <iomanip>

#include <def_use_xerces.hh>
#include <boost/version.hpp>
#include <boost/asio/ip/host_name.hpp>

#include "main/CFS.hh"
#include "Utils/Timer.hh"
#include "DataInOut/DefineFiles/DefineInOutFiles.hh"
#include "DataInOut/ParamHandling/MaterialHandler.hh"
#include "DataInOut/ProgramOptions.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "General/Environment.hh"
#include "DataInOut/ParamHandling/SkeletonConf.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/Xerces.hh"
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/ColoredConsole.hh"
#include <unistd.h>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "DataInOut/Logging/LogConfigurator.hh"

#include <def_use_mesh.hh>

#ifdef USE_MESH
#include "DataInOut/SimInOut/AnsysFile/SimInputMESH.hh"
#endif

using namespace CoupledField;
using namespace std;


int main(int argc, const char **argv)
{
  CFS cfs(argc, argv);
  int ret = cfs.Run();
  return ret;
}

void PrintWarning(CoupledField::Exception& ex ) {
  
  // Print warning on command line
 std::string msg = ex.GetMsg();
 std::string fileName = ex.GetFileName();
 UInt lineNum = ex.GetLineNum();
 
  std::cerr << "\n "
      << fg_blue << "WARNING:" << fg_reset << "\n "
      << msg << endl;
  std::cerr << "\n(" << fileName << ", Line " 
            << lineNum  << ")\n\n";
  
  // Print warning also to info xml
  PtrParamNode warn = info->Get("warning",ParamNode::INSERT);
  warn->Get("lineNum")->SetValue(lineNum);
  warn->Get("fileName")->SetValue(fileName);
  warn->Get("message")->SetValue(msg);
}
    

CFS::CFS(int argc, const char **argv) :
    timer(new Timer())
{

  resultHandler = NULL;
  materialHandler = NULL;

  logConf_ = new LogConfigurator();
  
  // Set segfault to false
  Exception::segfault_ = false;


  // =========================================================================
  // HANDLE COMMAND LINE PARAMETERS
  // =========================================================================
  progOpts = new ProgramOptions(argc, argv);

  // Parse command line
  progOpts->ParseData();

  // Log program startup
  progOpts->GetHeaderString( cout );
  
  // Initialize logging class (read parameters from file if desired)
  logConf_->ParseLogConfFile();
  
  // Get information about exception handling
  Exception::segfault_ = progOpts->GetForceSegFault();
  
  // Register callback function with exception class for warning
  Exception::SetCallbackWarn(&PrintWarning);

  // Set global Enums, the rest is set by the classes
  SetGlobalEnums();

  // the new xml logging derived from the ParamNode
  info = PtrParamNode(new ParamNode(ParamNode::INSERT, ParamNode::ELEMENT ));
      //progOpts->GetSimName() + ".info.xml", "<?xml version=\"1.0\"?>");
  info->SetName("cfsInfo");
  info->Get("status")->SetValue("running"); // to be overwritten by "aborted" or "finished"
  info->Get(ParamNode::PN_SUMMARY)->Get("timer")->SetValue(timer);
  timer->Start(); // ignore that this is not the real beginning

  // Print information about program start time and host
  using namespace boost::posix_time;
  using namespace boost::gregorian;

  // our calculation environment
  PtrParamNode env = info->Get(ParamNode::PN_HEADER)->Get("environment");
  start_time_ = to_simple_string( second_clock::local_time() );
  env->Get("started")->SetValue(start_time_);
  
  hostname_ = boost::asio::ip::host_name();
  
//  char host[256];
//  hostname_ = gethostname( host, sizeof(host) ) > 0 ? host : "";

  if(!hostname_.empty())
    env->Get("host")->SetValue(hostname_);
  
  info->ToFile("", true);

}


CFS::~CFS()
{

  // flush last information.
  info->ToFile(std::string(), true);
  

  // Delete objects
  //delete param;
  delete domain; // might write ersatz material file if <export save="finally"/> in optimization
  delete progOpts;
  delete resultHandler;

  // delete some global objects because valgrind complains otherwise
  // does not really matter anyway...
  param.reset();
  info.reset();

  delete logConf_;
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
    
    cout << ">> Total time: wall clock: '";
    
    const int walltime((int) timer->GetWallTime());
    const int cputime((int) timer->GetCPUTime());

    if(walltime > 120) 
    {
      const int wallmin((int) (walltime / 60.0));
      const int cpumin((int) (cputime / 60.0));
      if(wallmin > 60)
      {
        cout << wallmin / 60 << "h " << (wallmin % 60) 
             << "m' CPU time: '" << cpumin / 60 << "h " << (cpumin % 60) << "m'"; 
      }
      else
      {
        cout << wallmin << "m " << (walltime % 60) 
             << "s' CPU time: '" << cpumin << "m " << (cputime % 60) << "s'"; 
      }
    }
    else
    {
      cout << walltime << "s' CPU time: '" 
           << cputime << "s'";
    }
    
    cout << endl << endl;
      
    // write the info object
    info->Get("status")->SetValue("finished"); // overwrite 'running'
    info->Get(ParamNode::PN_SUMMARY)->Get("memory/final")->SetValue(MemoryUsage(false));
    info->Get(ParamNode::PN_SUMMARY)->Get("memory/peak")->SetValue(MemoryUsage(true));

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
      PtrParamNode errorNode = info->Get(ParamNode::PN_ERROR);
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
  ElemShape::Initialize();
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
 
 using namespace boost::posix_time;
 using namespace boost::gregorian;
 cout << "\n++ Finished solving the problem at " 
     << to_simple_string( second_clock::local_time() ) 
     << endl;

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
  SimInput * ptInputfile = new SimInputMESH(meshFile, PtrParamNode());
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
  string schema = progOpts->GetSchemaPathStr();
  schema += "/CFS-Simulation/CFS.xsd";

  // Initialize our xerces dom parser to handle the cfs xml file
  Xerces* xerces = new Xerces(schema);
  xerces->SetFile(xmlFile);

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

//  PtrParamNode meshNode = param->Get("input", false );
//  if( meshNode )
//    meshNode->Get("meshLibrary", libmesh, false);

//  if ( libmesh != "structGrid" ) {
    // Generate mesh reader
  fileHandler.CreateSimInputFiles( inFiles, gridInputs );
//  }

  // generate material handler
  materialHandler = fileHandler.CreateMaterialHandler();

  // Create simulation output writer
  map<string, shared_ptr<SimOutput> > outFiles;
  map<string, string> outGridIds;
  fileHandler.CreateSimOutputFiles( outFiles, outGridIds );

  // Create resulthandler and pass the output files
  resultHandler = new ResultHandler( ResultHandler::EMBEDDED );
  map<string, shared_ptr<SimOutput> >::iterator outputIt;
  map<string, shared_ptr<SimInput> >::iterator inputIt;
  outputIt = outFiles.begin();
  inputIt = inFiles.begin();
  for( ; outputIt != outFiles.end(); outputIt++ ) {
    resultHandler->AddOutputDest( outputIt->second, 
                                  outputIt->first,
                                  outGridIds[outputIt->first] );
  }
  for( ; inputIt != inFiles.end(); inputIt++ ) {
    resultHandler->AddInputReader( inputIt->second, inputIt->first );
  }


  // Log command line parameters
  progOpts->ToInfo(info->Get(ParamNode::PN_HEADER)->Get("progOpts"));
  // log the optinal id/name/token/label from <cfsSimulation id="..">
  info->Get(ParamNode::PN_HEADER)->Get("id")->SetValue(param->Get("id"));
  // if requested five the problem file -> one can see the defaults then
  if(progOpts->DoDetailedInfo())
    info->Get(ParamNode::PN_HEADER)->Get("cfsSimulation")->SetValue(param);
  // Open file for status reports by OLAS
  fileHandler.OpenFile( OLAS_FILE );
}

