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
#include "Utils/mathParser/mathParser.hh"
#include "DataInOut/DefineInOutFiles.hh"
#include "DataInOut/SimState.hh"
#include "DataInOut/SimInOut/hdf5/SimOutputHDF5.hh"
#include "DataInOut/ParamHandling/MaterialHandler.hh"
#include "DataInOut/ProgramOptions.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Driver/FormsContexts.hh"
#include "General/Environment.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/XmlReader.hh"
#include "DataInOut/ResultHandler.hh"
#include "DataInOut/ColoredConsole.hh"
#include <omp.h>
#include <sched.h>
#if not defined(WIN32) 
#  include <unistd.h>
#endif
#include <boost/date_time/posix_time/posix_time.hpp>
#include <def_use_petsc.hh>

#ifdef USE_PETSC
#include "petsc.h"
#include "OLAS/external/petsc/PETSCSolver.hh"
#endif




#include "DataInOut/SimInOut/hdf5/SimInputHDF5.hh"
#include "PDE/SinglePDE.hh"
#include "Utils/PythonKernel.hh"

#ifdef USE_MKL
#include <mkl_service.h>
#ifndef mkl_get_version
#define MKL_Get_Version MKLGetVersion
#define MKL_Free_Buffers MKL_FreeBuffers
#endif
#endif


using namespace CoupledField;
using namespace std;
using namespace boost::posix_time;
using namespace boost::gregorian;


// Create global info node
PtrParamNode infoNode;

#ifdef USE_PETSC
int main(int argc, const char **argv)
{
  PetscInitialize(NULL,NULL,PETSC_NULL,PETSC_NULL);
  int rank;
  int size;
  int ret =0;
  //find which is my rank
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  if(rank==0)
  {
    try
    {
      CFS cfs(argc, argv);  // possibly a std::exception
      ret = cfs.Run();      // catches all
    }
    catch(const std::exception& ex)
    {
      CFS::HandleException(ex);
      ret = 0;
    }
    //Send a Kill Tag to all workers before exiting the code
    if(size>1)
      for(rank = 1; rank < size; ++rank)
        MPI_Send(0, 0, MPI_INT, rank, DIETAG, MPI_COMM_WORLD);
  }
  else
  {
    PETSCWorker w(argc, argv);
    w.run();
  }
  PetscFinalize();
  return ret;
}
#else
int main(int argc, const char **argv)
{
  try
  {
    CFS cfs(argc, argv);
    int ret = cfs.Run();
    return ret;
  }
  catch(const exception& ex) // is catched in Run but can also come from constructor
  {
    CFS::HandleException(ex);
    return 1;
  }
}
#endif // USE_PETSC

void PrintWarning(const CoupledField::Exception& ex ) {
  
  // Print warning on command line
  std::string msg = ex.GetMsg();
  std::string fileName = ex.GetFileName();
  UInt lineNum = ex.GetLineNum();
 
  std::cerr << "\n " << fg_blue << "WARNING:" << fg_reset << "\n " << msg << endl;
  std::cerr << "\n(" << fileName << ", Line " << lineNum  << ")\n\n";
  
  // Print warning also to info xml
  PtrParamNode warn = infoNode->Get("warning",ParamNode::INSERT);
  warn->Get("lineNum")->SetValue(lineNum);
  warn->Get("fileName")->SetValue(fileName);
  warn->Get("message")->SetValue(msg);
}
    

CFS::CFS(int argc, const char **argv) :
  timer(new Timer())
{
  timer->Start();

  resultHandler = NULL;

  // Set segfault to false
  Exception::segfault_ = false;

  // =========================================================================
  // HANDLE COMMAND LINE PARAMETERS
  // =========================================================================
  progOpts = new ProgramOptions(argc, argv);

  // Parse command line, also initializes logging
  progOpts->ParseData();

  // Log program startup
  progOpts->PrintHeader(cout);

  // homogenize CFS_/OMP_/MKL_NUM_THREADS for mkl systems based on command line and environment setting
  // for openblas and netlib it is CFS_/OMP_/DUMMY_NUN_THREADS.
  // for Apple's accelerate it is CFS_/OMP_NUM_THREADS and VECLIB_MAXIMUM_THREADS
  SetNumberOfThreads(progOpts->GetNumThreads(), true, false); // yes, homogenize, no, don't be quiet
  
  // when we have the threads num set, we can print them
  progOpts->PrintNumThreads(cout, progOpts->IsQuiet());

  // Get information about exception handling
  Exception::segfault_ = progOpts->GetForceSegFault();
  
  // Set global Enums, the rest is set by the classes
  SetGlobalEnums();

  // the new xml logging derived from the ParamNode
  infoNode = ParamNode::GenerateWriteNode("cfsInfo", progOpts->GetSimName() + ".info.xml", ParamNode::INSERT, true, true); // lazy write and add counters
  infoNode->Get("status")->SetValue("running"); // to be overwritten by "aborted" or "finished"
  infoNode->Get(ParamNode::SUMMARY)->Get("timer")->SetValue(timer);

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
  
  delete domain;
  domain = NULL;
  
  // might write ersatz material file if <export save="finally"/> in optimization
  delete progOpts;
  progOpts = NULL;

  delete python;
  python = nullptr;
  
  if (simState) {
    simState->Finalize();
    simState.reset();
  }

  // delete some global objects because valgrind complains otherwise
  // does not really matter anyway...
  paramNode_.reset();
  infoNode.reset();

  #ifdef USE_MKL
  MKL_Free_Buffers();
  #endif
}

int CFS::Run()
{
  try
  {
    if(!progOpts->IsQuiet())
    {
      cout << "Simulation run started at " << start_time_;
      if(!hostname_.empty()) cout << " on " << hostname_;
      cout << endl;
    }

    ReadXMLFile(); // sets paramNode_

    // the python kernel only does stuff when compiled with USE_EMBEDDED_PYTHON and when we actually load a module (script)
    python = new PythonKernel(paramNode_->Get("python", ParamNode::PASS), infoNode);

    SetupIO(paramNode_);

    // =====================================================================================================================
    // ============================= precice start =========================================================================
#ifdef USE_PRECICE
    // according to the tutorial on how to couple to precice https://precice.org/couple-your-code-steering-methods.html,
    // we need rank and size of the "current thread". However, in openCFS, we are using
    // openMP. MPI is only used when compiling with PETSC
    int rank = 0;//omp_get_thread_num(); //0;
    int size = 1;
    // Initialize preCICE participant
    std::string participantName, configFileName;
    paramNode_->Get("fileFormats/preciceCoupling/participantName")->GetValue("name", participantName);
    paramNode_->Get("fileFormats/preciceCoupling/precice_configFile")->GetValue("file", configFileName);

    participant_ = new precice::Participant(participantName, configFileName, rank, size);

    participant_->initialize();
    // participant.advance(10);
    // participant.finalize();
#endif
    // =====================================================================================================================
    // ============================= precice end =========================================================================




    domain = new Domain( gridInputs, resultHandler, materialHandler,
        simState, paramNode_, infoNode );

    // Create grid
    domain->CreateGrid();

    // if a python function is registered, call it.
    python->CallHook(PythonKernel::POST_GRID);

    if(progOpts->GetPrintGrid())
      PrintGrid();
    else
      SolveProblem();

    python->CallHook(PythonKernel::POST_SOLVE_PROBLEM);


    // wait for all drivers to be initialized before printing the math parser variables
    domain->GetMathParser()->ToInfo(infoNode->Get(ParamNode::HEADER)->Get("domain/globalMathParser"), MathParser::GLOB_HANDLER);

    timer->Stop();

    cout << endl;
    cout << ">> Total wall-clock time: '" << Timer::GetTimeString(timer->GetWallTime())
                 << "' cpu time: '" << Timer::GetTimeString(timer->GetCPUTime())
                 << "' at " << to_simple_string(second_clock::local_time()) << endl << endl;

    // write the info object
    infoNode->Get("status")->SetValue("finished"); // overwrite 'running'

    infoNode->Get(ParamNode::SUMMARY)->SetComment("memory in MB");
    infoNode->Get(ParamNode::SUMMARY)->Get("memory/final")->SetValue(MemoryUsage(false)/1024);
    int peak = MemoryUsage(true); // macOS has no peak
    if(peak > 0)
      infoNode->Get(ParamNode::SUMMARY)->Get("memory/peak")->SetValue(peak/1024);

    return 0;
  }
  catch(const mu::Parser::exception_type& e)
  {
    cerr << endl << "mu::ParserError: '" << e.GetMsg() << "'"; // for expression '" << e.GetExpr() << "'" << endl;
    return 1;
  }
  catch(const exception& ex)
  {
    CFS::HandleException(ex);
    return 1;
  }
  catch (...)
  {
    cerr << "leftover exception caught:" << endl;
    cerr << boost::current_exception_diagnostic_information() << endl;
    return 1;
  }
}


void CFS::HandleException(const std::exception& ex)
{
  // long or brief output?
  if(progOpts != NULL && !progOpts->IsQuiet())
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

  cerr.flush();

  // Print error cause to info file
  if(infoNode != NULL)
  {
    PtrParamNode errorNode = infoNode->Get(ParamNode::FAIL);
    errorNode->SetValue(ex.what());
    infoNode->Get("status")->SetValue("aborted");
    infoNode->ToFile();
  }
}

void CFS::SetGlobalEnums()
{
  SetEnvironmentEnums();
  BasePDE::SetEnums();
  BiLinFormContext::SetEnums();
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

 python->CallHook(PythonKernel::POST_DOMAIN_INIT);

 // Solves the driver or optimization problem
 domain->SolveProblem();
 
 if(!progOpts->IsQuiet())
   cout << "\n++ Finished solving the problem at " << to_simple_string(second_clock::local_time()) << endl;
}

void CFS::ReadXMLFile()
{
  // Generate parameter handler and pass address to global pointer
  string xmlFile = progOpts->GetParamFileStr();

  // Conditionally write information to command line
  if(!progOpts->IsQuiet())
    cout << "++ Reading parameter file '" + xmlFile + "'" << endl;

  // this is the new param stuff which replaces the old params - delete this comment finally
  string schema = progOpts->GetSchemaPathStr() + "/CFS-Simulation/CFS.xsd";

  // parse the problem xml file, validate and fill with defaults from schema
  // continue to work only with the ParamNode tree
  paramNode_ = XmlReader::ParseFile(xmlFile, schema, "http://www.cfs++.org/simulation");

  // paramNode_->Dump();
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
  
  // log the optional id/name/token/label from <cfsSimulation id="..">
  std::string id = progOpts->GetId() != "" ? progOpts->GetId() : paramNode_->Get("id")->As<std::string>();
  infoNode->Get(ParamNode::HEADER)->Get("id")->SetValue(id);
  
  // additional log for all kind of information
  if (paramNode_->Has("info"))
    infoNode->Get(ParamNode::HEADER)->Get("info")->SetValue(paramNode_->Get("info"),false);
  // if requested give the problem file -> one can see the defaults then
  if(progOpts->DoDetailedInfo())
    infoNode->Get(ParamNode::HEADER)->Get("cfsSimulation")->SetValue(paramNode_,false);
}

