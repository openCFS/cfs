#include "BEM_PDE.hh"

#include <boost/lexical_cast.hpp>
#include <boost/filesystem.hpp>
#include <cmath>

#include "General/defs.hh"

// Bem_PDE::Bem_PDE()
// {
//   std::cout << "Bem_PDE::Bem_PDE()\n";
//   callNiHu_1();
// }

void Bem_PDE::callNiHu_1()
{
  // output pre
  const char* prefix_nihu = "++ NIHU >> ";
  const char* error_nihu = "++ NIHU_ERROR >> ";

  // Get the current working directory
  boost::filesystem::path current_directory = boost::filesystem::current_path();

  // Get the name of the current directory
  std::string directory = current_directory.stem().string();

  // Set the output directory based on the current directory name
  std::string output = directory + "_NiHu_output";

  
  if (chdir((directory + "_NiHu").c_str()) != 0)
  {
    std::cerr << error_nihu << "no directory " + directory + "_NiHu" << std::endl;
    return;
  }

  std::cout << prefix_nihu << "generate directory build" << std::endl;
  
  if (std::system("mkdir -p build") != 0)
  {
    std::cerr << error_nihu << "cannot generate directory build" << std::endl;
    return;
  }

  std::cout << prefix_nihu << "opening directory build" << std::endl;
  
  if (chdir("build") != 0)
  {
    std::cerr << error_nihu << "no directory build" << std::endl;
    return;
  }

  std::cout << prefix_nihu << "calling cmake" << std::endl;
  
  if (std::system("cmake ..") != 0)
  {
    std::cerr << error_nihu << "cannot call cmake properly" << std::endl;
    return;
  }

  std::cout << prefix_nihu << "building testcase" << std::endl;
  
  if (std::system("make") != 0)
  {
    std::cerr << error_nihu << "cannot execute make properly" << std::endl;
    return;
  }

  std::cout << prefix_nihu << "open result(executable)" << std::endl;
  
  if (std::system(("./" + output).c_str()) != 0)
  {
    std::cerr << error_nihu << "no executable found" << std::endl;
    return;
  }

  std::cout << prefix_nihu << "leaving directory NiHu/build" << std::endl;
  
  if (chdir("../..") != 0)
  {
    std::cerr << error_nihu << "possible segmentation fault" << std::endl;
    return;
  }
}