#include <iostream>
#include "xercesc/util/PlatformUtils.hpp"
#include "SPHandler.hh"

using namespace xercesc;

using namespace CoupledField;

int main(int argC, char* argV[])
{

  std::cout << std::endl;
  std::cout << " *********************************" << std::endl;
  std::cout << "   Parsing Parameter File:" << std::endl;
  std::cout << " *********************************" << std::endl;

  // Check and report input parameters
  if ( argC == 2 ) {
    std::cout << "\n --> Using input XML file: " << argV[1] << std::endl;
  }
  else {
    std::cout << "\n Please specify an XML file\n\n";
  }

  // Parse XML file
  SPHandler params( argV[1] );
  std::cout << " --> Parsing went fine\n\n";
  std::cout << " *********************************\n" << std::endl;
  return 0;

}
