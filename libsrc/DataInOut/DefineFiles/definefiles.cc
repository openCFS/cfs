#include <fstream>
#include <iostream>
#include <string>

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
 if (InfoPrint)
       {
         infofile = new std::ofstream(strcat(filename,".info"));
         if (!infofile) Error("Can't open info-file");
       }

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
      std::cerr << "ERROR: Sorry, we can't read files with type: "<< atype << std::endl;
      exit(-1);
    }
}

} // end of namespace
