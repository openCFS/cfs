//#include <string>
#include <fstream>
#include <iostream>

#include <general_head.hh>
#include <utils_head.hh>

#include "filetype.hh"
#include "datfile.hh"
#include "definefiles.hh"

namespace CoupledField
{

DefineInOutFiles :: DefineInOutFiles(const Char * name)
{
 
 filename=new Char[100];
 strcpy(filename, name);
 
#ifdef TRACE
 strcpy(filename, name);
 trace=new std::ofstream(strcat(filename,".trace"));
#endif
 
#ifdef DEBUG
 strcpy(filename, name);
 debug=new std::ofstream(strcat(filename,".deb"));
#endif
 
 strcpy(filename, name);
 if (InfoPrint) infofile = new std::ofstream(strcat(filename,".info"));

 strcpy(filename, name); 
}
 
DefineInOutFiles ::~ DefineInOutFiles()
{
 
delete filename;
 
#ifdef TRACE
delete trace;
#endif
 
#ifdef DEBUG
delete debug;
#endif
 
if (InfoPrint) delete infofile;
}

FileType * DefineInOutFiles :: Create_ptFileType(Char * atype)
{
  if (!strcmp(atype, "-dat"))  
    {
      infiletype=new DatFile(filename);
      return infiletype;
    }
  else 
    {
      std::cerr << "ERROR: Sorry, we can't read files with type: "<<atype << std::endl;
      exit(-1);
    }
}

} // end of namespace
