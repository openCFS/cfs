#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

BaseTransfer :: BaseTransfer(Integer asize, Integer adof)
{
#ifdef TRACE
  (*trace) << "entering BaseTransfer::BaseTransfer" << endl;
#endif

  size = asize;
  dof  = adof;
}
  
BaseTransfer :: ~BaseTransfer()
{
#ifdef TRACE
  (*trace) << "entering BaseTransfer::~BaseTransfer" << endl;
#endif

}

}
