#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

RScalarTransfer :: RScalarTransfer(Integer asize, Integer * rsw, Boolean mem)
  : BaseTransfer(asize)
{
#ifdef TRACE
  (*trace) << "entering RScalarTransfer::RScalarTransfer" << endl;
#endif

  if (mem)
    {
      Integer i,nne;
      
      start = new Integer[size+1];
      
      start[0] = 0;
      nne      = 0;
      
      for (i=1; i<=size; i++)
	{
	  start[i] = start[i-1]+rsw[i-1];
	  nne += rsw[i-1];
	}
      
      val = new Double[nne];
      pos = new Integer[nne];

      calculated = FALSE;
    }
  else
    {
      val   = NULL;
      start = NULL;
      pos   = NULL;

      calculated = TRUE;
    }
}
  
RScalarTransfer :: ~RScalarTransfer()
{
#ifdef TRACE
  (*trace) << "entering RScalarTransfer::~RScalarTransfer" << endl;
#endif
}

void RScalarTransfer :: Calc(BaseTopology * topology)
{
#ifdef TRACE
  (*trace) << "entering RScalarTransfer::Calc" << endl;
#endif

  if (calculated)
    {
      return;
    }

  Integer i,j,ind,maxnh,ncn;
  Double weight;

  maxnh = topology->MaxNumNeighbour();

  Integer * cnh   = new Integer[maxnh];

  for (i=0; i<size; i++)
    {
      ind = topology->GetCNN(i+1);

      if (ind > 0)      // coarse grid
	{
	  val[start[i]] = 1;
	  pos[start[i]] = ind;
	}
      else if (ind < 0) // fine grid
	{
	  ncn = topology->GetStrongCoarseNeighbour(i+1, cnh);
	  
	  if (ncn != 0)
	    {
	      weight = 1./ncn;
	    }

	  for (j=0; j<ncn; j++)
	    {
	      val[start[i]+j] = weight;
	      pos[start[i]+j] = topology->GetCNN(cnh[j]);
	    }
	}
      else
	{
	  cout << i+1 << " " << ind << endl;
	  cout << "error in CoupledField::RScalarTransfer::Calc" << endl;
	  cout << "no specific grid point" << endl;
	  exit(1);
	}
    }

  calculated = TRUE;
}

void RScalarTransfer :: MultHh(BaseVector & uH, BaseVector & uh)
{
#ifdef TRACE
  (*trace) << "entering RScalarTransfer::MultHh" << endl;
#endif
  
  RealVector & uc = (RealVector &) uH;
  RealVector & uf = (RealVector &) uh;

  Integer i,j,q,rs;
  Double sum;

  for (i=0; i<size; i++)
    {
      q   = start[i];
      rs  = start[i+1] - q;
      sum = 0;

      for (j=0; j<rs; j++)
	{
	  sum += val[q]*uc.Get(pos[q]);
	  q++;
	}

      uf.Elem(i+1) += sum;
    }
}

void RScalarTransfer :: MulthH(BaseVector & uh, BaseVector & uH)
{
#ifdef TRACE
  (*trace) << "entering RScalarTransfer::MulthH" << endl;
#endif

  RealVector & uc = (RealVector &) uH;
  RealVector & uf = (RealVector &) uh;

  Integer i,j,q,rs;
  Double sum;

  for (i=0; i<size; i++)
    {
      q   = start[i];
      rs  = start[i+1] - q;
      sum = uf.Get(i+1);
      
      for (j=0; j<rs; j++)
	{
	  uc.Elem(pos[q]) += val[q]*sum;
	  q++;
	}
    }
}

void RScalarTransfer :: Print() const
{
#ifdef TRACE
  (*trace) << "entering RScalarTransfer::Print" << endl;
#endif
  
  Integer i,j,rs;

  for (i=0; i<size; i++)
    {
      rs = start[i+1] - start[i];

      cout << "row number " << i+1 << endl;

      for (j=0; j<rs; j++)
	{
	  cout << val[start[i]+j] << " ";
	}
      cout << endl;

      for (j=0; j<rs; j++)
	{
	  cout << pos[start[i]+j] << " ";
	}
      cout << endl;
    }

}

}
