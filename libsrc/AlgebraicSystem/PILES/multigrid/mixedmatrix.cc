#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

MixedMatrix :: MixedMatrix()
  : BaseMatrix()
{
#ifdef TRACE
  (*trace) << "entering MixedMatrix::MixedMatrix" << endl;
#endif

}
  
MixedMatrix :: ~MixedMatrix()
{
#ifdef TRACE
  (*trace) << "entering MixedMatrix::~MixedMatrix" << endl;
#endif
}

void MixedMatrix :: Mult(BaseVector & vec1, BaseVector & vec2, Double factor) const
{
}

void MixedMatrix :: Assemble(Double * v,Integer * p, Integer elemsize)
{
}

void MixedMatrix :: Print() const
{
}
}
