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
