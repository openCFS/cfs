#include <fstream>
#include <iostream>
#include <string>

#include "definefiles.hh"
#include "datfile.hh"
#include "ansysfile.hh"
#include "conffile.hh"
#include "outUnverg.hh"
#include "outGMV.hh"

namespace CoupledField
{

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

 strcpy(filename, name);
 cla=new std::ofstream(strcat(filename,".cla"));
 if (!cla) Error("Can't open cla++-file");
 
 strcpy(filename, name);
 if (InfoPrint)
       {
         infofile = new std::ofstream(strcat(filename,".info"));
         if (!infofile) Error("Can't open info-file");
       }

 strcpy(filename, name); 

 conf=new ConfFile(name);
 if (!conf) Error("Can't open conf-file");

ptWriteResults=NULL;

}
 
DefineInOutFiles ::~DefineInOutFiles()
{
#ifdef TRACE
 (*trace) << "Entering DefineInOutFiles::~DefineInOutFiles" << std::endl;
#endif

#ifdef DEBUG
delete debug;
#endif
 
if (infofile) delete infofile;

if (conf) delete conf;
if (cla) delete cla;

if (ptWriteResults) delete ptWriteResults;

if (infiletype) delete infiletype;

 delete [] filename;

#ifdef TRACE
 delete trace;
#endif
}

FileType * DefineInOutFiles :: Create_ptFileType()
{
  std::string informat;
  conf->get("format_input",informat);

  if (informat=="dat")  
      infiletype=new DatFile(filename);
  else
  if  (informat=="mesh") 
      infiletype=new AnsysFile(filename);
  else 
        Error("Wrong format for input file. Please, check your data.",__FILE__,__LINE__);

   return infiletype;
}

WriteResults * DefineInOutFiles :: Create_ptWriteResults()
{
  std::string outformat;
  conf->get("format_output",outformat);

  if (outformat=="gmv") ptWriteResults=new WriteResultsGMV(filename);
  else 
    if (outformat=="unverg") ptWriteResults=new WriteResultsUnverg(filename);
      else
        Error("Wrong format for writing results. Please, check your data.",__FILE__,__LINE__);

  if (!ptWriteResults) 
    Error("Can't open file for output results",__FILE__,__LINE__);

  return ptWriteResults;
}

} // end of namespace
