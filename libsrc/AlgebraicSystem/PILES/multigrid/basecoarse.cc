#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

BaseCoarse :: BaseCoarse(Integer amaxnh)
{
#ifdef TRACE
  (*trace) << "entering BaseCoarse::BaseCoarse" << endl;
#endif

  maxnh = amaxnh;
}
  
BaseCoarse :: ~BaseCoarse()
{
#ifdef TRACE
  (*trace) << "entering BaseCoarse::~BaseCoarse" << endl;
#endif
}


void BaseCoarse :: Calc(Integer * rsw)
{
#ifdef TRACE
  (*trace) << "entering BaseCoarse::Calc" << endl;
#endif

  Integer actpoint;

  //startpoint = 1;

  while (startpoint != 0)
    {
      actpoint = startpoint;

      while (actpoint != 0)
	{
	  topology->SetCoarseGridPoint(actpoint);
	  topology->SetFineGridPoint(actpoint);
	  
	  actpoint = topology->GetNextCoarseGridPoint(actpoint, maxnh);
	}
      
      startpoint = topology->GetNextStartPoint();
    }

  topology->SetRSW(rsw);
  topology->CalcCoarseGraph();
}

void BaseCoarse :: Print() const
{
  topology->Print();
}
}
