#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

BaseGraph :: BaseGraph(Integer asize)
{
#ifdef TRACE
  (*trace) << "entering BaseGraph::BaseGraph" << endl;
#endif

  size = asize;

  graph = new Graph[size];

  Integer i;

  for (i=0; i<size; i++)
    {
      graph[i].pos     = new Integer[50];
      graph[i].actsize = 0;
      graph[i].size    = 50;
    }

#ifdef MEMTRACE
  double dmb;
  double imb;

  dmb = 0;
  imb = 50*size*4./1e6;

  sumdmem += dmb;
  sumimem += imb;

  (*memtrace) << "+++ ALLOCATE MEMORY: double  BaseGraph        " << dmb << " MB" << endl;
  (*memtrace) << "+++ ALLOCATE MEMORY: integer BaseGraph        " << imb << " MB" << endl;
#endif

  //cm = new Integer[];
}
  
BaseGraph :: ~BaseGraph()
{
#ifdef TRACE
  (*trace) << "entering BaseGraph::~BaseGraph" << endl;
#endif

  Integer i;
  
  for (i=0; i<size; i++)
    {
      delete [] graph[i].pos;
    }
  
  delete graph;
}

void BaseGraph :: SetElementPos(Integer * connect, Integer elemsize)
{
#ifdef TRACE
  (*trace) << "entering BaseGraph::SetElementPos" << endl;
#endif

  Integer i,j,k,l,p1,p2;
  Boolean q;

  for (i=0; i<elemsize; i++)
    {
      p1 = -1;

      if (connect[i] > 0)
	{
	  p1 = connect[i]-1;
	}
      else if (connect[i] < 0)
	{
	  p1 = cm[abs(connect[i])]-1;
	}

      for (j=0; j<elemsize; j++)
	{
	  p2 = -1;

	  if (connect[j] > 0)
	    {
	      p2 = connect[j];
	    }
	  else if (connect[j] < 0)
	    {
	      p2 = cm[abs(connect[j])];
	    }

	  if (p1 != -1 && p2 != -1)
	    {
	      l = graph[p1].actsize;
	      k = 0;
	      q = TRUE;
	      
	      while (k < l && q == TRUE)
		{
		  if (graph[p1].pos[k] == p2)
		    {
		      q = FALSE;
		    }
		  
		  k++;
		}
	      
	      if (q == TRUE)
		{
		  graph[p1].pos[l] = p2;
		  graph[p1].actsize++;
		  
		  if (graph[p1].actsize >= graph[p1].size)
		    {
		      cout << endl;
		      cout << "Predefined value for neighbors in matrix graph too small! " << endl;
                      cout << " neeeded value is: " << p1 << endl; 
		      cout << endl;
                      exit (EXIT_FAILURE);
		    }
		}
	    }
	}
    }
}

void BaseGraph :: Create()
{
#ifdef TRACE
  (*trace) << "entering BaseGraph::Create" << endl;
#endif

  Integer i,j,k,l,rs;

  nne = 0;

  for (i=0; i<size; i++)
    {
      rs   = graph[i].actsize;
      nne += rs;

      for (j=0; j<rs; j++)
	{
	  if (graph[i].pos[j] == i+1)
	    {
	      k = j;
	    }
	}

      graph[i].pos[k] = graph[i].pos[0];
      graph[i].pos[0] = i+1;

      for (j=1; j<rs; j++)
	{
	  for (k=j+1; k<rs; k++)
	    {
	      if (graph[i].pos[j] > graph[i].pos[k])
		{
		  l = graph[i].pos[j];
		  graph[i].pos[j] = graph[i].pos[k];
		  graph[i].pos[k] = l;
		}
	    }
	}
    }
}

void BaseGraph :: Print() const
{
#ifdef TRACE
  (*trace) << "entering BaseGraph::Print" << endl;
#endif
  Integer i,j;

  (*test) << "printing auxiliary graph" << endl;

  for (i=0; i<size; i++)
    {
      (*test) << "line " << i+1 << endl;

      for (j=0; j<graph[i].actsize; j++)
	{
	  (*test) << graph[i].pos[j] << " ";
	}
      (*test) << endl;
    }
}

}
