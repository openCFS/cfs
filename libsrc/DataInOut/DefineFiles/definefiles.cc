#include <fstream>
#include <iostream>
#include <string>

#include "definefiles.hh"
#include "DataInOut/AnsysFile/ansysfile.hh"
#include "DataInOut/conffile.hh"
#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"

#ifdef PARALLEL
#include <mpi.h>
#endif

namespace CoupledField
{

#ifdef MEMTRACE
 Double sumdmem;
 Double sumimem;
#endif

DefineInOutFiles :: DefineInOutFiles(const Char * name)
{

 #ifdef PARALLEL
 // find out who I am and write to seperate files:
 int commrank;
 MPI_Comm comm = MPI_COMM_WORLD;
 MPI_Comm_rank(comm,&commrank);
 char *rank = new char[5];
 sprintf(rank,"_%d",commrank);
 #endif

 filename=new Char[100];



#ifdef TRACE
 strcpy(filename, name);
 #ifdef PARALLEL
 strcat(filename, rank);
 #endif
 trace=new std::ofstream(strcat(filename,".trace"));
 if (!trace) Error("Can't open trace-file");
#endif
 
#ifdef DEBUG
 strcpy(filename, name);
 #ifdef PARALLEL
 strcat(filename, rank);
 #endif
 debug=new std::ofstream(strcat(filename,".deb"));
 if (!debug) Error("Can't open debug-file");
#endif

#ifdef MEMTRACE
 #ifdef PARALLEL
 strcat(filename, rank);
 #endif
 strcpy(filename, name);
 CoupledField::memtrace=new std::ofstream(strcat(filename,".mem"));
 if (!memtrace) Error("Can't open memtrace-file");
#endif

 strcpy(filename, name);
  #ifdef PARALLEL
 strcat(filename, rank);
 #endif
 cla=new std::ofstream(strcat(filename,".las")); 
 if (!cla) Error("Can't open LAS++-file");


 strcpy(filename, name);
 data=new std::ofstream(strcat(filename,".data"));
 if (!data) Error("Can't open data-file");


 strcpy(filename, name);
 conf=new ConfFile(name);
 if (!conf) Error("Can't open conf-file");

 ptWriteResults=NULL;

 flags=new Flags();

}
 
DefineInOutFiles ::~DefineInOutFiles()
{
#ifdef TRACE
 (*trace) << "Entering DefineInOutFiles::~DefineInOutFiles" << std::endl;
#endif

#ifdef DEBUG
delete debug;
#endif
 
if (ptWriteResults) delete ptWriteResults;

if (conf) delete conf;
if (cla) delete cla;

if (infiletype) delete infiletype;

delete [] filename;

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
      infiletype=new AnsysFile(filename);
  else 
        Error("Wrong format for input file. Please, check your data.",__FILE__,__LINE__);

   return infiletype;
}

WriteResults * DefineInOutFiles :: Create_ptWriteResults(FileType * const aInFile)
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
    ptWriteResults=new WriteResultsGMV(filename,withHistory, aInFile);
  else if (outformat=="unverg")
    ptWriteResults=new WriteResultsUnverg(filename, withHistory, aInFile);
  else
    Error("Wrong format for writing results. Please, check your data.",__FILE__,__LINE__);
  
  if (!ptWriteResults) 
    Error("Can't open file for output results",__FILE__,__LINE__);
  
  return ptWriteResults;
}

} // end of namespace
