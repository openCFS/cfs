//#include <string>
#include <fstream>
#include <iostream>
#include <stdarg.h>

#include <general_head.hh>
#include <utils_head.hh>

#include "filetype.hh"

namespace CoupledField
{

FileType :: FileType(const Char * const afilename)
{
#ifdef TRACE
  (*trace) << "entering FileType::FileType" << std::endl;
#endif

  filename = new Char[100];
  strcpy(filename,afilename);
}
  
FileType :: ~FileType()
{
#ifdef TRACE
  (*trace) << "entering FileType::~FileType" << std::endl;
#endif

 delete [] filename ;
}

}
