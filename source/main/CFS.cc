// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim:fenc=utf-8:ft=cpp:et:sw=2:ts=2:sts=2
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;


#include <iomanip>
#include <fstream>
#include <boost/version.hpp>
#include <boost/asio/ip/host_name.hpp>
#include <boost/exception/diagnostic_information.hpp>

#include "main/CFS.hh"
#include "Utils/Timer.hh"
#include "DataInOut/DefineInOutFiles.hh"
#include "DataInOut/SimState.hh"
#include "DataInOut/SimInOut/hdf5/SimOutputHDF5.hh"
#include "DataInOut/ParamHandling/MaterialHandler.hh"
#include "DataInOut/ProgramOptions.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "General/Environment.hh"
#include "DataInOut/ParamHandling/SkeletonConf.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/XmlReader.hh"
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/ColoredConsole.hh"
#include <unistd.h>
#include <boost/date_time/posix_time/posix_time.hpp>
#include "DataInOut/Logging/LogConfigurator.hh"

#include <def_use_mesh.hh>

#ifdef USE_MESH
#include "DataInOut/SimInOut/AnsysFile/SimInputMESH.hh"
#endif


#include "DataInOut/SimInOut/hdf5/SimInputHDF5.hh"
#include "PDE/SinglePDE.hh"

using namespace CoupledField;
using namespace std;
using namespace boost::posix_time;
using namespace boost::gregorian;


#ifdef __MINGW32__

#if 0
extern "C" {
int __security_cookie;
}

extern "C" void _fastcall __security_check_cookie(int i) {
//do nothing
}
#endif

// extern "C" void _chkstk() {
//do nothing
// }
extern "C" void _allmul() {
//do nothing
}
#endif //__MINGW32__

// Create global info node
PtrParamNode infoNode;

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
  PtrParamNode warn = infoNode->Get("warning",ParamNode::INSERT);
  warn->Get("lineNum")->SetValue(lineNum);
  warn->Get("fileName")->SetValue(fileName);
  warn->Get("message")->SetValue(msg);
}
    

CFS::CFS(int argc, const char **argv) :
    timer(new Timer())
{

  resultHandler = NULL;
  materialHandler = NULL;

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
  std::string confFile = progOpts->GetLogConfFileStr();
  logConf_ = new LogConfigurator(confFile);
  logConf_->ParseLogConfFile();
  
  // Get information about exception handling
  Exception::segfault_ = progOpts->GetForceSegFault();
  
  // Set global Enums, the rest is set by the classes
  SetGlobalEnums();

  // the new xml logging derived from the ParamNode
  infoNode = ParamNode::GenerateWriteNode("cfsInfo", progOpts->GetSimName() + ".info.xml", ParamNode::INSERT, true, true); // lazy write and add counters
  infoNode->Get("status")->SetValue("running"); // to be overwritten by "aborted" or "finished"
  infoNode->Get(ParamNode::SUMMARY)->Get("timer")->SetValue(timer);
  timer->Start(); // ignore that this is not the real beginning

  // Register callback function with exception class for warning
  Exception::SetCallbackWarn(&PrintWarning);
  
  // Print information about program start time and host
  using namespace boost::posix_time;
  using namespace boost::gregorian;

  // our calculation environment
  PtrParamNode env = infoNode->Get(ParamNode::HEADER)->Get("environment");
  start_time_ = to_simple_string( second_clock::local_time() );
  env->Get("started")->SetValue(start_time_);
  
  hostname_ = boost::asio::ip::host_name();
  
//  char host[256];
//  hostname_ = gethostname( host, sizeof(host) ) > 0 ? host : "";

  if(!hostname_.empty())
    env->Get("host")->SetValue(hostname_);
  
  infoNode->ToFile("", true);

}


CFS::~CFS()
{

  // flush last information.
  infoNode->ToFile(std::string(), true);


  delete resultHandler;
  resultHandler = NULL;
  
  // Delete objects
  //delete param;
  
  
  delete domain;
  domain = NULL;
  
  // might write ersatz material file if <export save="finally"/> in optimization
  delete progOpts;
  progOpts = NULL;
  
  if (simState) {
    simState->Finalize();
    simState.reset();
  }

  // delete some global objects because valgrind complains otherwise
  // does not really matter anyway...
  paramNode_.reset();
  infoNode.reset();

  delete logConf_;
  logConf_ = NULL;
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
    SetupIO(paramNode_);

    domain = new Domain( gridInputs, resultHandler, materialHandler, 
                         simState, paramNode_, infoNode );
    // Create grid
    domain->CreateGrid();

    if(progOpts->GetPrintGrid())
      PrintGrid();
    else{
      SolveProblem();
    }

    // wait for all drivers to be initialized before printing the math parser variables
    domain->GetMathParser()->ToInfo(infoNode->Get(ParamNode::HEADER)->Get("domain/globalMathParser"), MathParser::GLOB_HANDLER);

    timer->Stop();
    if(!progOpts->IsQuiet())
      cout << endl; // conditional empty line
    
    cout << ">> Total time: wall clock: '";
    
    int walltime = (int) timer->GetWallTime();
    double cputime = timer->GetCPUTime();

    if(walltime > 120) 
    {
      int wallmin((int) (walltime / 60.0));
      int cpumin((int) (cputime / 60.0));
      if(wallmin > 60)
        cout << wallmin / 60 << "h " << (wallmin % 60) << "m' CPU time: '" << cpumin / 60 << "h " << (cpumin % 60) << "m'";
      else
        cout << wallmin << "m " << (walltime % 60) << "s' CPU time: '" << cpumin << "m " << ((int) cputime % 60) << "s'";
    }
    else
    {
      cout << walltime << "s' CPU time: '" << cputime << "s'";
    }
    if(progOpts->IsQuiet())
      cout << " at " << to_simple_string(second_clock::local_time()) << endl;

    
    cout << endl << endl;
      
    // write the info object
    infoNode->Get("status")->SetValue("finished"); // overwrite 'running'
    infoNode->Get(ParamNode::SUMMARY)->SetComment("memory in MB");
    infoNode->Get(ParamNode::SUMMARY)->Get("memory/final")->SetValue(MemoryUsage(false)/1024.);
    infoNode->Get(ParamNode::SUMMARY)->Get("memory/peak")->SetValue(MemoryUsage(true)/1024.);

    return 0;
  }
  catch(const mu::ParserError& e)
  {
    cerr << endl << "mu::ParserError: " << e.GetMsg() << endl;
    return 1;
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
    if(infoNode != NULL)
    {
      PtrParamNode errorNode = infoNode->Get(ParamNode::FAIL);
      errorNode->SetValue(ex.what());
      infoNode->Get("status")->SetValue("aborted");
      infoNode->ToFile();
    }

    return 1;
  }
  catch (...)
  {
    cerr << "leftover exception caught:" << endl;
    cerr << boost::current_exception_diagnostic_information() << endl;
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
 
 if(!progOpts->IsQuiet())
   cout << "\n++ Finished solving the problem at " << to_simple_string(second_clock::local_time()) << endl;
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
  SimInput * ptInputfile = new SimInputMESH(meshFile, PtrParamNode(), infoNode);
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

  // Conditionally write information to command line
  if(!progOpts->IsQuiet())
    cout << "++ Reading parameter file '" + xmlFile + "'" << endl;

  // this is the new param stuff which replaces the old params - delete this comment finally
  string schema = progOpts->GetSchemaPathStr();
  schema += "/CFS-Simulation/CFS.xsd";

  // parse the problem xml file, validate and fill with defaults from schema
  // continue to work only with the ParamNode tree
  paramNode_ = XmlReader::ParseFile(xmlFile, schema);
}

void CFS::SetupIO(PtrParamNode rootNode )
{
  
  // Create structure for handling the simulation outputstate
  simState.reset(new SimState(false));
  
  // read meshes
  map<string, shared_ptr<SimInput> > inFiles;
  fileHandler.CreateSimInputFiles( rootNode, infoNode, inFiles, gridInputs );
  
  // generate material handler
  materialHandler = fileHandler.CreateMaterialHandler(rootNode );

  // Create simulation output writer
  map<string, shared_ptr<SimOutput> > outFiles;
  map<string, string> outGridIds;
  fileHandler.CreateSimOutputFiles( rootNode, infoNode, outFiles, outGridIds );

  // Create result handler and pass the output files
  resultHandler = new ResultHandler( paramNode_ );
  map<string, shared_ptr<SimOutput> >::iterator outputIt;
  map<string, shared_ptr<SimInput> >::iterator inputIt;
  outputIt = outFiles.begin();
  inputIt = inFiles.begin();
  
  // check for hdf5Reader
  shared_ptr<SimOutputHDF5> hdf5Writer;
  
  for( ; outputIt != outFiles.end(); outputIt++ ) {
    resultHandler->AddOutputDest( outputIt->second, 
                                  outputIt->first,
                                  outGridIds[outputIt->first] );

    // check if writer has hdf5 format
    // a dynamic cast leads for icc to error: cannot convert pointer to base class "CoupledField::SimOutput" to pointer to derived class "CoupledField::SimOutputHDF5" -- base class is virtual
    // in boost/smart_ptr/shared_ptr.hpp(805) itself :(
    if( typeid(*outputIt->second) == typeid(SimOutputHDF5) ) 
      hdf5Writer = boost::dynamic_pointer_cast<SimOutputHDF5>(outputIt->second);
  }
  for(; inputIt != inFiles.end(); inputIt++) 
    resultHandler->AddInputReader(inputIt->second, inputIt->first);
  
  // Pass hdf5 writer to simState class
  simState->SetOutputHdf5Writer(hdf5Writer);
  simState->SetMatParamFile(progOpts->GetParamFileStr(), materialHandler->GetFileName());
  
  // Log command line parameters
  progOpts->ToInfo(infoNode->Get(ParamNode::HEADER)->Get("progOpts"));
  
  // log the optinal id/name/token/label from <cfsSimulation id="..">
  infoNode->Get(ParamNode::HEADER)->Get("id")->SetValue(paramNode_->Get("id"));
  
  // if requested give the problem file -> one can see the defaults then
  if(progOpts->DoDetailedInfo())
    infoNode->Get(ParamNode::HEADER)->Get("cfsSimulation")->SetValue(paramNode_);
}

