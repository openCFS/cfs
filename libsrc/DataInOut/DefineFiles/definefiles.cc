#include <fstream>
#include <iostream>
#include <string>

#include "definefiles.hh"
#include "DataInOut/AnsysFile/ansysfile.hh"
#include "DataInOut/conffile.hh"
#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/ParamHandling/PlainXMLParamHandler.hh"
#include "DataInOut/ParamHandling/XMLParamHandler.hh"

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


#ifdef PARALLEL
 // find out who I am and write to seperate files:
 int commrank;
 MPI_Comm comm = MPI_COMM_WORLD;
 MPI_Comm_rank(comm,&commrank);
 char *rank = new char[5];
 sprintf(rank,"_%d",commrank);
#endif


  // ===============
  //   Constructor
  // ===============
  DefineInOutFiles::DefineInOutFiles(const Char *name)
  {

    // Store the basename of the auxilliary files
    filename_ = new Char[strlen(name)+1];
    strcpy( filename_, name );

    // Generate array for names of auxilliary files
    Char *auxfile = new Char[strlen(name)+MAXPOSTFIX+MAXRANK+1];

    // Generate basename for files
    Char *basename = new Char[strlen(name)+MAXRANK+1];
    strcpy( basename, name );

#ifdef PARALLEL
    strcat( basename, rank );
#endif

#ifdef TRACE
    strcpy(auxfile, basename);
    trace = new std::ofstream(strcat(auxfile,".trace"));
    if (!trace) Error("Can't open trace-file");
#endif
 
#ifdef DEBUG
    strcpy(auxfile, basename);
    debug = new std::ofstream(strcat(auxfile,".deb"));
    if (!debug) Error("Can't open debug-file");
#endif

#ifdef MEMTRACE
    strcpy(auxfile, basename);
    CoupledField::memtrace = new std::ofstream(strcat(auxfile,".mem"));
    if (!memtrace) Error("Can't open memtrace-file");
#endif

    // File for reports of linear algebra system LAS
    strcpy(auxfile, basename);
    cla = new std::ofstream(strcat(auxfile,".las")); 
    if (!cla) Error("Can't open LAS++-file");


    strcpy(auxfile, basename);
    data=new std::ofstream(strcat(auxfile,".data"));
    if (!data) Error("Can't open data-file");

    // Generate configuration file object and pass address to global pointer
    // Note: This is the old-fashioned conffile format
    conf = new ConfFile(name);
    if (!conf) Error("Can't open conf-file");

    // Generate parameter handler and pass address to global pointer
    // Note: This is the new XML-based conffile format
    strcpy( auxfile, name );
    strcat( auxfile, ".xml" );

#ifdef USE_XERCES
    params = new XMLParamHandler( auxfile );
#else
    params = new PlainXMLParamHandler( auxfile );
#endif

    // Initialise internal pointers
    ptWriteResults_ = NULL;
    infileType_ = NULL;

    // ???
    flags = new Flags();

  }
 

  // ==============
  //   Destructor
  // ==============
  DefineInOutFiles ::~DefineInOutFiles()
  {
#ifdef TRACE
    (*trace) << "Entering DefineInOutFiles::~DefineInOutFiles" << std::endl;
#endif

#ifdef DEBUG
    delete debug;
#endif
 
    if (conf) delete conf;
    if (cla) delete cla;

    // Delete internal pointers
    if (ptWriteResults_) delete ptWriteResults_;
    if (infileType_) delete infileType_;
    delete[] filename_;

#ifdef TRACE
    delete trace;
#endif

 delete flags;

}

  FileType * DefineInOutFiles :: Create_ptFileType()
  {
    std::string informat="mesh";
    conf->ifget("format_input",informat);

    if  (informat=="mesh")
      {
	infileType_=new AnsysFile(filename_);
      }
    else
      {
	Error("Wrong format for input file. Please, check your data.",
	      __FILE__, __LINE__);
      }

    return infileType_;
  }

  WriteResults*
  DefineInOutFiles::Create_ptWriteResults(FileType *const aInFile)
  {
    std::string outformat="unverg";
    conf->ifget("format_output",outformat);
    
    Boolean withHistory=FALSE;
    Integer val;
    if (conf->ifget("history_node",val))
      if (val!=-1) 
	withHistory=TRUE;
  

    // save nodes may also be listed in config-command "save_nodes"
    std::vector<std::string> historyList; 
    conf->ifgetliststr("save_nodes", historyList);
    if (historyList.size())
      withHistory=TRUE;

    if (outformat=="gmv")
      ptWriteResults_=new WriteResultsGMV(filename_,withHistory, aInFile);
    else if (outformat=="unverg")
      ptWriteResults_=new WriteResultsUnverg(filename_, withHistory, aInFile);
    else
      Error("Wrong format for writing results. Please, check your data.",
	    __FILE__, __LINE__);
  
    if (!ptWriteResults_) 
      Error("Can't open file for output results",__FILE__,__LINE__);
  
    return ptWriteResults_;
  }

} // end of namespace
