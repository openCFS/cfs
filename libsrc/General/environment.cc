#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>
#include <string>

#include "environment.hh"

namespace CoupledField
{
ostream * trace = NULL ;
ostream * debug  = NULL;
ostream * infofile=NULL;

Boolean InfoPrint=FALSE;

ostream & operator << (ostream & out, const enum precond & type)
{
switch (type)
{
case Jacobi: out << " Jacobi "; break;
case SSOR: out << " SSOR "; break;
case LU: out << " LU "; break;
case non: out << " non "; break;

return out;
}
} // end of ostream
}
