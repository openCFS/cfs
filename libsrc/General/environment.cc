#include <iostream>
#include <fstream>
#include <string>

#include "environment.hh"


// Since OLAS uses a separate namespace for 
// writing out data, two different declarations
// have to be made
#ifdef USE_OLAS
namespace OutInfo{
  std::ostream * trace = NULL ;
  std::ostream * debug  = NULL;
  std::ostream * cla=NULL;
  std::ostream * memtrace=NULL;
  std::ostream * data = NULL;
}
#else
namespace CoupledField{
  std::ostream * trace = NULL ;
  std::ostream * debug  = NULL;
  std::ostream * cla=NULL;
  std::ostream * memtrace=NULL;
  std::ostream * data = NULL;
}
#endif



namespace CoupledField
{
  Boolean PrintGridOnly = FALSE;

  Flags * flags=NULL;

  BaseFE * ptQ    = NULL;
  BaseFE * ptQ2   = NULL;
  BaseFE * ptTet = NULL;
  BaseFE * ptL1 = NULL;
  BaseFE * ptL2 = NULL;
  BaseFE * ptTr1 = NULL;
  BaseFE * ptTr2 = NULL;
  BaseFE * ptHexa=NULL;
  BaseFE * ptPyra=NULL;
  BaseFE * ptWedge=NULL;

  WriteInfo * Info = NULL;
  ConfFile * conf           = NULL;
  BaseParamHandler * params = NULL;

  std::ostream& operator<< (std::ostream & outStr, std::vector<Double> xOut)
  {
    for (Integer i=0; i<xOut.size(); i++)
      outStr <<  " " << xOut[i];
    return outStr;
  }

  std::ostream& operator<< (std::ostream & outStr, std::vector<Integer> xOut)
  {
    for (Integer i=0; i<xOut.size(); i++)
      outStr <<  " " << xOut[i];
    return outStr;
  }
}
