#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>

#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/filesystem/convenience.hpp"
#include "boost/filesystem/exception.hpp"
namespace fs = boost::filesystem;

#include "boost/algorithm/string/trim.hpp"

#include "sysparh.h"

#include "def_use_ansysrst.hh"

#include "General/Environment.hh"
#include "ResWrProtos.hh"

namespace CoupledField
{

void GetANSYSRevisionInfos(std::string& BinLibRev,
                           std::string& MachineId,                           
                           std::string& CompiledRev)
{
  char buffer[1024];
  std::string file_name = "prima.rst";
  int len;

  try 
  {
    // If file already exists, delete it.
    // ANSYS API can not handle this case.
    if(fs::exists(file_name))
    {
      fs::remove( file_name );
    }
  } catch (std::exception &ex)
  {
    std::cerr << ex.what() << std::endl;
  }

  len = file_name.length();
  sprintf(buffer, "%s", file_name.c_str());
  reswrtest_(buffer, &len, len);

  std::ifstream ifstr(file_name.c_str(), std::ios::binary);
  
  ifstr.seekg(44,std::ios::beg);
  ifstr.read(buffer, 4);
  buffer[4] = 0;
  BinLibRev = buffer;
  boost::algorithm::trim(BinLibRev);
  
  ifstr.seekg(52,std::ios::beg);
  ifstr.read(buffer, 12);
  buffer[12] = 0;
  MachineId = buffer;
  boost::algorithm::trim(MachineId);
  ifstr.close();

  CompiledRev = STR_VERSION;
  boost::algorithm::trim(CompiledRev);
  
  try 
  {
    // If file already exists, delete it.
    // ANSYS API can not handle this case.
    if(fs::exists(file_name))
    {
      fs::remove( file_name );
    }
  } catch (std::exception &ex)
  {
    std::cerr << ex.what() << std::endl;
  }

}
}

int main(int argc, char** argv) 
{
  std::string revision;
  std::string compiled_rev;
  std::string machine_id;

  CoupledField::GetANSYSRevisionInfos(revision,
                                      machine_id,
                                      compiled_rev);

  std::cout << revision << " "
            << machine_id << " compiled rev: "
            << compiled_rev
            << std::endl;
  

  return 0;  
}
