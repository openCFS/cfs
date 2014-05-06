#include <iostream>

#include "General/Exception.hh"

#include "ScatteredDataReaderCCM.hh"

using namespace std;

int main( int argc, char *argv[])
{
  CoupledField::ScatteredDataReaderCCM SCRCCM(argv[1], true);

  try 
  {
    SCRCCM.ReadInput();
    SCRCCM.WriteCellCenters();
  } catch ( CoupledField::Exception& ex ) 
  {
    std::cout << ex.what() << std::endl;
  }  
  
  return 0;
}
