#include <fstream>
#include <iostream>
#include <string>
#include <stdarg.h>

#include "ansysfile.hh"

namespace CoupledField
{

AnsysFile :: AnsysFile(const Char * const afilename)
 :FileType(afilename)
{
#ifdef TRACE
  (*trace) << "entering AnsysFile::AnsysFile" << std::endl;
#endif

  infile.open(strcat(filename,".ansys"));
  if (!infile) {std::cerr << "ERROR(" << __FILE__ << " " << __LINE__ <<
                         ") Can't open " << filename << std::endl;
                exit(1);}

  infile.seekg(0, std::ios::end);
  pos_end=infile.tellg();

}
  
AnsysFile :: ~AnsysFile()
{
#ifdef TRACE
  (*trace) << "entering AnsysFile::~AnsysFile" << std::endl;
#endif

 infile.close() ;
}

}
