#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

BaseParameter :: BaseParameter()
{
#ifdef TRACE
  (*trace) << "entering BaseParameter::BaseParameter" << endl;
#endif

}
  
BaseParameter :: ~BaseParameter()
{
#ifdef TRACE
  (*trace) << "entering BaseParameter::~BaseParameter" << endl;
#endif
}

void BaseParameter :: Set() 
{
#ifdef TRACE
  (*trace) << "entering BaseParameter::Set" << endl;
#endif

  eps           = 1e-8;
  epsint        = 0;
  epsmat        = 0;
  epsmach       = 1e-16;
  dampiter      = 0.8;
  dampsmooth    = 0.6;
  alpha         = 0.25;
  beta          = 0.0;
  penalty       = 1e8;
  precond       = 2;
  solver        = 2;
  coarsesystem  = 50000;
  coarsesolver  = 1;
  maxnumiter    = 300;
  coarsening    = 1;
  maxneighbour  = 200;
  transfer      = 1;
  norm          = 1;
  cycle         = 1;
  smoothingfor  = 1;
  smoothingback = 1;
  smoothingtype = 1;
  auxmat        = 0;
  spemat        = 0;
  outback       = 0;
}
}
