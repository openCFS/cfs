#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

BaseVector :: BaseVector()
{
#ifdef TRACE
  (*trace) << "entering BaseVector::BaseVector" << endl;
#endif

}
  
BaseVector :: ~BaseVector()
{
#ifdef TRACE
  (*trace) << "entering BaseVector::~BaseVector" << endl;
#endif
}

}
