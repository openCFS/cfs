#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

BaseSmoother :: BaseSmoother(BaseParameter & param, Integer asize, Integer anumrhs, Integer adof)
{
#ifdef TRACE
  (*trace) << "entering BaseSmoother::BaseSmoother" << endl;
#endif

  size   = asize;
  numrhs = anumrhs;
  dof    = adof;

  numsmoothfor  = param.GetSmoothingStepFor();
  numsmoothback = param.GetSmoothingStepBack();
}
  
BaseSmoother :: ~BaseSmoother()
{
#ifdef TRACE
  (*trace) << "entering BaseSmoother::~BaseSmoother" << endl;
#endif
  
}
}
