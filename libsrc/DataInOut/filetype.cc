#include <string>
#include <fstream.h>
#include <iostream.h>
#include <stdarg.h>

#include <general_head.hh>
#include <utils_head.hh>

#include "filetype.hh"

namespace CoupledField
{

FileType :: FileType(const Char * const afilename)
{
#ifdef TRACE
  (*trace) << "entering FileType::FileType" << endl;
#endif

  filename = new Char[100];
  strcpy(filename,afilename);
}
  
FileType :: ~FileType()
{
#ifdef TRACE
  (*trace) << "entering FileType::~FileType" << endl;
#endif

 delete [] filename ;
}

}
