#include <iostream>
#include <fstream>
#include <string>

#include "environment.hh"

namespace CoupledField
{
  std::ostream * trace = NULL ;
  std::ostream * debug  = NULL;
  std::ostream * infofile=NULL;
  std::ostream * cla=NULL;

  ConfFile * conf=NULL;

  Boolean InfoPrint=FALSE;

  BaseFE * ptQ   = NULL;
  BaseFE * ptTet = NULL;

  // BaseElem * ptQ = NULL;
  // BaseElem * ptTr = NULL;
  // BaseElem * ptTet = NULL;
  // BaseElem * ptL1=NULL;
  // BaseElem * ptHexa=NULL;

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
}
