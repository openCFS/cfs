#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

BaseSolver :: BaseSolver()
{
#ifdef TRACE
  (*trace) << "entering BaseSolver::BaseSolver" << endl;
#endif

}
  
BaseSolver :: ~BaseSolver()
{
#ifdef TRACE
  (*trace) << "entering BaseSolver::~BaseSolver" << endl;
#endif
}

}
