#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "basematrix.hh"
#include "sparsematrix.hh"
#include "complexmatrix.hh"

namespace CoupledField
{

ComplexMatrix :: ComplexMatrix()
  : BaseMatrix()
{
#ifdef TRACE
  (*trace) << "entering ComplexMatrix::ComplexMatrix" << endl;
#endif

}
  
ComplexMatrix :: ~ComplexMatrix()
{
#ifdef TRACE
  (*trace) << "entering ComplexMatrix::~ComplexMatrix" << endl;
#endif
}
}
