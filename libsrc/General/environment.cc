#include <iostream>
#include <fstream>
#include <string>

#include "environment.hh"

namespace CoupledField
{
  std::ostream * trace = NULL ;
  std::ostream * debug  = NULL;
  std::ostream * cla=NULL;
  std::ostream * memtrace=NULL;
  std::ostream * data = NULL;

  Boolean PrintGridOnly = FALSE;

  Flags * flags=NULL;
  ConfFile * conf=NULL;

  BaseFE * ptQ   = NULL;
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
  
  std::ostream & operator << (std::ostream & out, const enum precond & type)
  {
    switch (type)
      {
      case Jacobi: out << " Jacobi "; break;
      case SSOR: out << " SSOR "; break;
      case LU: out << " LU "; break;
      case non: out << " non "; break;
      }
    return out;
  } // end of ostream


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
