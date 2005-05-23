#include <fstream>
#include <iostream>

#include "filetype.hh"

namespace CoupledField
{


  FileType :: FileType(const Char * const afilename)
  {
    ENTER_FCN( "FileType::FileType" );

    filename = new Char[100];
    strcpy(filename,afilename);
  }
  
  FileType :: ~FileType()
  {
    ENTER_FCN( "FileType::~FileType" );

    delete [] filename ;
  }


}
