#include <iostream>
#include <cstring>
#include "xercesc/util/PlatformUtils.hpp"
#include "SPHandler.hh"

using namespace xercesc;

using namespace CoupledField;

int main(int argC, char* argV[])
{

  char *xmlFile = NULL;
  bool showTree = false;

  // Check and report input parameters
  if ( argC == 2 ) {
    xmlFile = argV[1];
  }
  else if ( argC == 3 ) {
    if ( strcmp( argV[1], "-showTree" ) == 0 ) {
      showTree = true;
    }
    else {
      showTree = false;
      std::cout << "\n Ignoring command line option '" << argV[1] << "'\n"
		<< " Did you mean '-showTree'?" << std::endl;
    }
    xmlFile = argV[2];
  }
  else {
    std::cout << "\n Please call the TestParser using the following syntax:\n"
	      << "\n ./parse <name of xml file>\n\n"
	      << " or \n\n"
	      << " ./parse -showTree <name of xml file>\n\n";
    exit(1);
  }


  // ******************
  //   Parse XML file
  // ******************
  std::cout << "\n *********************************\n"
	    << "   Parsing Parameter File:\n"
            << " *********************************\n"
	    << "\n --> Using input XML file: " << xmlFile << std::endl;

  SPHandler params( xmlFile );

  std::cout << " --> Parsing went fine\n\n"
            << " *********************************\n" << std::endl;


  // *****************************
  //   Output DOM Tree to stdout
  // *****************************
  if ( showTree == true ) {
    std::cout << " *********************************\n"
	      << "   Exporting DOM Tree:\n"
	      << " *********************************\n\n";
    params.PrintTree();
  }
  return 0;

}
