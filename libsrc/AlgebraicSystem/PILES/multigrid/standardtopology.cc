#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

StandardTopology :: StandardTopology(Integer asize, Integer anne, Integer * pos, Integer * start)
  : BaseTopology(asize, anne, pos, start)
{
#ifdef TRACE
  (*trace) << "entering StandardTopology::StandardTopology" << endl;
#endif
  
}

StandardTopology :: ~StandardTopology()
{
#ifdef TRACE
  (*trace) << "entering StandardTopology::~StandardTopology" << endl;
#endif

}

void StandardTopology :: CalcCoarseGraph()
{
#ifdef TRACE
  (*trace) << "entering StandardTopology::CalcCoarseGraph" << endl;
#endif

  Integer i,j,k,l,ncn,mcn,nfn,maxnh,c1,c2;

  maxnh = MaxNumNeighbour();

  Integer * cnh   = new Integer[maxnh];
  Integer * hnh   = new Integer[maxnh];
  Integer * fnh   = new Integer[maxnh];
  crs             = new Integer[csize*offset];
  nrs             = new Integer[csize];

  for (i=0; i<csize; i++)
    {
      nrs[i] = 0;
    }
  
  for (i=0; i<csize*offset; i++)
    {
      crs[i] = 0;
    }

  /// build the graph - short connections & long connections!

  for (i=0; i<size; i++)
    {
      if (cnn[i] < 0)
	{
	  nfn = GetFineNeighbour(i+1,fnh);
	  ncn = GetStrongCoarseNeighbour(i+1, cnh);
	  
	  for (j=0; j<ncn; j++)
	    {
	      c1 = cnn[cnh[j]-1];

	      for (k=j+1; k<ncn; k++)
		{
		  c2 = cnn[cnh[k]-1];

		  Insert(c2, c1);
		}
	    }

	  for (j=0; j<ncn; j++)
	    {
	      c1 = cnn[cnh[j]-1];

	      for (k=0; k<nfn; k++)
		{
		  mcn = GetStrongCoarseNeighbour(fnh[k], hnh);
		  
		  for (l=0; l<mcn; l++)
		    {
		      if (hnh[l] != cnh[j])
			{
			  c2 = cnn[hnh[l]-1];
			 
			  Insert(c2,c1);
			}
		    }
		}
	    }
	}
    }

  /// count the cnne

  cnne = 0;

  for (i=0; i<csize; i++)
    {
      cnne += nrs[i];
    }

  cnne += csize;

  delete [] cnh;
  delete [] fnh;
  delete [] hnh;
}
}
