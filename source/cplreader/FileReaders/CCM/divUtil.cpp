
#include <iostream>
#include <limits>
#include "divUtil.h"
#include <stdlib.h>	// For exit()


namespace CCMDualizer {

  void WaitReturn() {
    std::cout << "Press ENTER to continue...";
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  void Exit() {
    exit(0);
  }

}



