#include <fstream>
#include <iostream>
#include <string>

//#include <general_head.hh>
//#include <utils_head.hh>

//#include "filetype.hh"

#include "definefiles.hh"
#include "datfile.hh"

namespace CoupledField
{

DefineInOutFiles :: DefineInOutFiles(const Char * name)
{

 filename=new Char[100];
 strcpy(filename, name);

#ifdef TRACE
 mark
 std::strcpy(filename, name);
 trace=new std::ofstream(std::strcat(filename,".trace"));
 if (!trace) Error("Can't open trace-file");
#endif
 
#ifdef DEBUG
 std::strcpy(filename, name);
 debug=new std::ofstream(std::strcat(filename,".deb"));
 if (!debug) Error("Can't open debug-file");
#endif
 
 std::strcpy(filename, name);
 if (InfoPrint)
       {
         infofile = new std::ofstream(std::strcat(filename,".info"));
         if (!infofile) Error("Can't open info-file");
       }

 std::strcpy(filename, name); 
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
      std::cerr << "ERROR: Sorry, we can't read files with type: "<< atype << std::endl;
      exit(-1);
    }
}

} // end of namespace
