#include <fstream>
#include <iostream>
#include <string>

#include "definefiles.hh"
#include "DataInOut/AnsysFile/ansysfile.hh"
#include "DataInOut/conffile.hh"
#include "DataInOut/Unverg/outUnverg.hh"
#include "DataInOut/GMV/outGMV.hh"

namespace CoupledField
{

#ifdef MEMTRACE
 Double sumdmem;
 Double sumimem;
#endif

DefineInOutFiles :: DefineInOutFiles(const Char * name)
{

 filename=new Char[100];
 strcpy(filename, name);

#ifdef TRACE
 strcpy(filename, name);
 trace=new std::ofstream(strcat(filename,".trace"));
 if (!trace) Error("Can't open trace-file");
#endif
 
#ifdef DEBUG
 strcpy(filename, name);
 debug=new std::ofstream(strcat(filename,".deb"));
 if (!debug) Error("Can't open debug-file");
#endif

#ifdef MEMTRACE
 strcpy(filename, name);
 CoupledField::memtrace=new std::ofstream(strcat(filename,".mem"));
 if (!memtrace) Error("Can't open memtrace-file");
#endif

 strcpy(filename, name);
 cla=new std::ofstream(strcat(filename,".las")); 
 if (!cla) Error("Can't open LAS++-file");

 strcpy(filename, name);
 data=new std::ofstream(strcat(filename,".data"));
 if (!data) Error("Can't open data-file");

//  strcpy(filename, name);
//  infofile = new std::ofstream(strcat(filename,".info"));
//  if (!infofile) Error("Can't open info-file");
 

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

if (infofile) delete infofile;

if (conf) delete conf;
if (cla) delete cla;

if (infiletype) delete infiletype;

delete [] filename;

#ifdef TRACE
 delete trace;
#endif
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

WriteResults * DefineInOutFiles :: Create_ptWriteResults()
{
  std::string outformat="unverg";
  conf->ifget("format_output",outformat);

  Boolean withHistory=FALSE;
  Integer val;
  if (conf->ifget("history_node",val))
	if (val!=-1) withHistory=TRUE;

  if (outformat=="gmv")
    ptWriteResults=new WriteResultsGMV(filename,withHistory);
  else if (outformat=="unverg")
    ptWriteResults=new WriteResultsUnverg(filename,withHistory);
  else
    Error("Wrong format for writing results. Please, check your data.",__FILE__,__LINE__);
  
  if (!ptWriteResults) 
    Error("Can't open file for output results",__FILE__,__LINE__);
  
  return ptWriteResults;
}

} // end of namespace
