#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

BaseMatrix :: BaseMatrix()
{
#ifdef TRACE
  (*trace) << "entering BaseMatrix::BaseMatrix" << endl;
#endif
  
  matfac = 1;
}
  
BaseMatrix :: ~BaseMatrix()
{
#ifdef TRACE
  (*trace) << "entering BaseMatrix::~BaseMatrix" << endl;
#endif

}

void BaseMatrix :: SetGraphRow(Integer end, Integer * apos, Integer row)
{
#ifdef TRACE
  (*trace) << "entering BaseGraph::SetGraphRow" << endl;
#endif  

  Integer i,ls;

  start[row] = start[row-1]+end;

  ls = start[row]-start[row-1];

  for (i=0; i<ls; i++)
    {
      pos[start[row-1]+i] = apos[i];
    }
}

Integer BaseMatrix :: GetGraphPos(Integer row, Integer col) const
{
#ifdef TRACE
  (*trace) << "entering BaseGraph::GetGraphPos" << endl;
#endif  

  Integer i,ls;

  ls = start[row]-start[row-1];

  for (i=0; i<ls; i++)
    {
      if (pos[start[row-1]+i] == col)
	{
	  return i+1;
	}
    }

  return 0;
}

void BaseMatrix :: SetCoarseGraph(BaseTopology * topology)
{
#ifdef TRACE
  (*trace) << "entering BaseGraph::SetCoarseGraph" << endl;
#endif  
  
  if (setgraph)
    {
      return;
    }

  Integer i,j,rs,s;
  Integer * p;

  for (i=0; i<size; i++)
    {
      rs = topology->GetCoarseRowSize(i+1);
      p  = topology->GetCoarsePos(i+1);
      s  = start[i];

      start[i+1] = s+rs+1;

      pos[s] = i+1;
      s++;

      for (j=0; j<rs; j++)
	{
	  pos[s] = p[j];
	  s++;
	}
    }

  setgraph = TRUE;
}

Integer BaseMatrix :: MaxRowSize() const
{
#ifdef TRACE
  (*trace) << "entering BaseGraph::MaxRowSize" << endl;
#endif  

  Integer i,maxrs;

  maxrs = 0;

  for (i=0; i<size; i++)
    {
      if ((start[i+1] - start[i]) > maxrs)
	{
	  maxrs = start[i+1] - start[i];
	}
    }

  return maxrs;
}

void BaseMatrix :: InitMatrix()
{
#ifdef TRACE
  (*trace) << "entering BaseGraph::InitMatrix" << endl;
#endif 

  Integer i;

  for (i=0; i<dof*dof*nne; i++)
    {
      val[i] = 0;
    }

  buildindir = FALSE;
}

void BaseMatrix :: CopyGraph(Integer * astart, Integer * apos)
{
#ifdef TRACE
  (*trace) << "entering BaseGraph::CopyGraph" << endl;
#endif 

  Integer i;

  for (i=0; i<=size; i++)
    {
      start[i] = astart[i];
    }

  for (i=0; i<nne; i++)
    {
      pos[i] = apos[i];
    }
}

}
