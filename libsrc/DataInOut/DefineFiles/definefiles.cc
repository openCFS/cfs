#include <fstream>
#include <iostream>
#include <cstdio>
#include <string>

#include "definefiles.hh"
#include "DataInOut/AnsysFile/ansysfile.hh"
#include "DataInOut/ParamHandling/ConfFile.hh"
#include "DataInOut/WriteInfo.hh"
#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"
#include "DataInOut/GSI/outGSI.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"
#include "DataInOut/ParamHandling/PlainXMLParamHandler.hh"
#include "DataInOut/ParamHandling/XMLParamHandler.hh"

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
  DefineInOutFiles::DefineInOutFiles(const Char *name)
  {

#ifdef PARALLEL
 // find out who I am and write to seperate files:
 int commrank;
 MPI_Comm comm = MPI_COMM_WORLD;
 MPI_Comm_rank(comm,&commrank);
 char *rank = new char[5];
 sprintf(rank,"_%d",commrank);
 
#endif

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

    // *************************
    //   Handle Parameter File
    // *************************

#ifndef XMLPARAMS

    // Generate configuration file object and pass address to global pointer
    // Note: This is the old-fashioned conffile format
    Info->StartProgress("Reading in .conf-file");
    conf = new ConfFile(name);
    if (!conf) Error("Can't open conf-file");
    Info->FinishProgress();

#else

    conf = NULL;

    // Generate parameter handler and pass address to global pointer
    // Note: This is the new XML-based conffile format
    strcpy( auxfile, name );
    strcat( auxfile, ".xml" );

#ifdef USE_XERCES
    params = new XMLParamHandler( auxfile );
#else
    params = new PlainXMLParamHandler( auxfile );
#endif

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
    memtrace = new std::ofstream(strcat(auxfile,".mem"));
    if (!memtrace) Error("Can't open memtrace-file");
#endif

    // File for reports of linear algebra system LAS
    strcpy(auxfile, basename);
    cla = new std::ofstream(strcat(auxfile,".las")); 
    if (!cla) Error("Can't open LAS++-file");


    strcpy(auxfile, basename);
    data=new std::ofstream(strcat(auxfile,".data"));
    if (!data) Error("Can't open data-file");


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
    // Cannot trace into trace file once trace file has been deleted
    // so ENTER_FCN cannot write a leaving message for this destructor!!!
    // ENTER_FCN( "DefineInOutFiles::~DefineInOutFiles" );

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


  // ==============================
  //   Generate mesh file pointer
  // ==============================
  FileType * DefineInOutFiles :: Create_ptFileType()
  {
    ENTER_FCN( "DefineInOutFiles::Create_ptFileType" );

#ifndef XMLPARAMS
    std::string informat="mesh";
    conf->ifget("format_input",informat);

    if  (informat=="mesh")
      {
	infileType_=new AnsysFile(filename_);
      }
    else
      {
	Error( "Wrong format for input file. Please, check your data!",
	       __FILE__, __LINE__ );
      }
#else
    if ( params->HasValue( "format", "mesh", "input" ) )
      {
	infileType_ = new AnsysFile( filename_ );
      }
    else
      {
	Error( "Wrong format for input file. Please, check your data!",
	       __FILE__, __LINE__ );
      }
#endif

    return infileType_;
  }


  // ================================
  //   Generate output file pointer
  // ================================
  WriteResults*
  DefineInOutFiles::Create_ptWriteResults(FileType *const aInFile)
  {
    ENTER_FCN( "DefineInOutFiles::Create_ptWriteResults" );

#ifndef XMLPARAMS
    std::string outformat="unverg";
    conf->ifget("format_output",outformat);
    
#else
    std::string outformat;
    params->Get( "format", outformat, "output" );
#endif

#ifndef XMLPARAMS
    if (outformat=="gmv")
      ptWriteResults_=new WriteResultsGMV(filename_, aInFile);
    else if (outformat=="unverg")
      ptWriteResults_=new WriteResultsUnverg(filename_, aInFile); 
    else if (outformat=="gsi")
      ptWriteResults_=new WriteResultsGSI(filename_, aInFile);
#ifdef USE_DATABASE
    else if (outformat=="database")
      ptWriteResults_=new WriteResultsDatabase(filename_, aInFile);
#endif
    else
      Error("Wrong format for writing results. Please, check your data.",
	    __FILE__, __LINE__);
#else
    if ( outformat == "gmv" ) {
      ptWriteResults_= new WriteResultsGMV(filename_, aInFile);
    }
    else if ( outformat == "unv" ) 
      ptWriteResults_= new WriteResultsUnverg(filename_, aInFile); 
    else if ( outformat == "gsi" ) 
      ptWriteResults_= new WriteResultsGSI(filename_, aInFile);
#ifdef USE_DATABASE
    else if (outformat == "database") 
      ptWriteResults_= new WriteResultsDatabase(filename_, aInFile);
#endif
    else
      {
	std::string errmsg = "Output format '" + outformat;
	errmsg += "' not recognised";
	Info->Error( errmsg, __FILE__, __LINE__ );
      }
#endif

    if (!ptWriteResults_) 
      Error("Can't open file for output results",__FILE__,__LINE__);
  
    return ptWriteResults_;
  }

} // end of namespace
