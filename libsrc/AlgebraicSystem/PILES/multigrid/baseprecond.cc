#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

BasePrecond :: BasePrecond()
{
#ifdef TRACE
  (*trace) << "entering BasePrecond::BasePrecond" << endl;
#endif

  offset = 15;

}
  
BasePrecond :: ~BasePrecond()
{
#ifdef TRACE
  (*trace) << "entering BasePrecond::~BasePrecond" << endl;
#endif

  delete mathier;
}

}
