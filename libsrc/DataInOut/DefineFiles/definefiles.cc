#include <fstream>
#include <iostream>
#include <string>

#include "definefiles.hh"
#include "datfile.hh"
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

// ptWriteResults2d=NULL;
// ptWriteResults3d=NULL;

}
 
DefineInOutFiles ::~ DefineInOutFiles()
{
#ifdef TRACE
 (*trace) << "Entering DefineInOutFiles::~DefineInOutFiles" << std::endl;
#endif
 
delete filename;
 
#ifdef TRACE
delete trace;
#endif
 
#ifdef DEBUG
delete debug;
#endif
 
if (InfoPrint) delete infofile;

if (conf) delete conf;
if (cla) delete cla;
}

FileType * DefineInOutFiles :: Create_ptFileType(Char * atype)
{
  if (!strcmp(atype, "-dat"))  
    {
      infiletype=new DatFile(filename);
    }
  else 
    {
      std::cerr << "ERROR: Sorry, we can't read files with type: "<< atype << std::endl;
      exit(-1);
    }
   return infiletype;
}

WriteResults<Point2D> * DefineInOutFiles :: Create_ptWriteResults2d()
{
  std::string outformat;
  conf->get("format_output",outformat);

  if (outformat=="gmv") ptWriteResults2d=new WriteResultsGMV<Point2D>(filename);
  else 
    if (outformat=="unverg") ptWriteResults2d=new WriteResultsUnverg<Point2D>(filename);
      else
        Error("Wrong format for writing results. Please, check your data.",__FILE__,__LINE__);

  if (!ptWriteResults2d) 
    Error("Can't open file for output results",__FILE__,__LINE__);

  return ptWriteResults2d;
}

WriteResults<Point3D> * DefineInOutFiles :: Create_ptWriteResults3d()
{
  std::string outformat;
  conf->get("format_output",outformat);

  if (outformat=="gmv") ptWriteResults3d=new WriteResultsGMV<Point3D>(filename);
  else
    if (outformat=="unverg") ptWriteResults3d=new WriteResultsUnverg<Point3D>(filename);
      else
        Error("Wrong format for writing results. Please, check your data.",__FILE__,__LINE__);

  if (!ptWriteResults3d)
      Error("Can't open file for output results",__FILE__,__LINE__);

  return ptWriteResults3d;
}

} // end of namespace
