#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

RBlockTransfer :: RBlockTransfer(Integer asize, Integer adof, Integer * rsw)
  : BaseTransfer(asize, adof)
{
#ifdef TRACE
  (*trace) << "entering RBlockTransfer::RBlockTransfer" << endl;
#endif

  Integer i,nne;
  
  start = new Integer[size+1];
  
  start[0] = 0;
  nne      = 0;
  
  for (i=1; i<=size; i++)
    {
      start[i] = start[i-1]+rsw[i-1];
      nne += rsw[i-1];
    }
  
  val = new Double[nne*dof*dof];
  pos = new Integer[nne];

  for (i=0; i<nne*dof*dof; i++)
    {
      val[i] = 0;
    }

  sum = new Double[dof];
  dm  = new DenseMatrix(dof,dof);

  calculated = FALSE;
}
  
RBlockTransfer :: ~RBlockTransfer()
{
#ifdef TRACE
  (*trace) << "entering RBlockTransfer::~RBlockTransfer" << endl;
#endif

  if (sum != NULL)
    {
      delete [] sum;
    }
}

void RBlockTransfer :: Calc(BaseTopology * topology)
{
#ifdef TRACE
  (*trace) << "entering RBlockTransfer::Calc" << endl;
#endif

  if (calculated)
    {
      return;
    }

  Integer i,j,k,p,ind,maxnh,ncn;
  Double weight;

  maxnh = topology->MaxNumNeighbour();

  Integer * cnh   = new Integer[maxnh];

  for (i=0; i<size; i++)
    {
      ind = topology->GetCNN(i+1);

      if (ind > 0)      // coarse grid
	{
	  p = start[i]*dof*dof;

	  for (k=0; k<dof; k++)
	    {
	      val[p+k*(dof+1)] = 1;
	    }

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
	      p = (start[i]+j)*dof*dof;
	      
	      for (k=0; k<dof; k++)
		{
		  val[p+k*(dof+1)] = weight;
		}

	      pos[start[i]+j] = topology->GetCNN(cnh[j]);
	    }
	}
      else
	{
	  cout << i+1 << " " << ind << endl;
	  cout << "error in CoupledField::RBlockTransfer::Calc" << endl;
	  cout << "no specific grid point" << endl;
	  exit(1);
	}
    }

  calculated = TRUE;

  delete [] cnh;
}

void RBlockTransfer :: MultHh(BaseVector & uH, BaseVector & uh)
{
#ifdef TRACE
  (*trace) << "entering RBlockTransfer::MultHh" << endl;
#endif
  
  RealVector & uc = (RealVector &) uH;
  RealVector & uf = (RealVector &) uh;

  Integer i,j,k,l,q,p,rs;

  DenseMatrix * x = new DenseMatrix(dof,dof);

  for (i=0; i<size; i++)
    {
      rs  = start[i+1] - start[i];

      for (k=0; k<dof; k++)
	{
	  sum[k] = 0;
	}

      for (j=0; j<rs; j++)
	{
	  p = (pos[start[i]+j]-1)*dof+1;
	  x = Get(i+1,j+1);

	  for (k=0; k<dof; k++)
	    {
	      for (l=0; l<dof; l++)
		{
		  sum[k] += x->Get(k+1,l+1)*uc.Get(p+l);
		}
	    }
	}

      p = i*dof+1;

      for (k=0; k<dof; k++)
	{
	  uf.Elem(p+k) += sum[k];
	}
    }
}

void RBlockTransfer :: MulthH(BaseVector & uh, BaseVector & uH)
{
#ifdef TRACE
  (*trace) << "entering RBlockTransfer::MulthH" << endl;
#endif

  RealVector & uc = (RealVector &) uH;
  RealVector & uf = (RealVector &) uh;

  Integer i,j,k,l,p,rs;

  DenseMatrix * x = new DenseMatrix(dof,dof);

  for (i=0; i<size; i++)
    {
      rs = start[i+1] - start[i];
      p  = i*dof+1;

      for (k=0; k<dof; k++)
	{
	  sum[k] = uf.Get(p+k);
	}

      for (j=0; j<rs; j++)
	{
	  p = (pos[start[i]+j]-1)*dof+1;
	  x = Get(i+1,j+1);

	  for (k=0; k<dof; k++)
	    {
	      for (l=0; l<dof; l++)
		{
		  uc.Elem(p+k) += x->Get(l+1,k+1)*sum[l];
		}
	    }
	}
    }
}

void RBlockTransfer :: Print() const
{
#ifdef TRACE
  (*trace) << "entering RBlockTransfer::Print" << endl;
#endif
  
  Integer i,j,p,rs;

  for (i=0; i<size; i++)
    {
      rs = start[i+1] - start[i];
      p  = start[i]*dof*dof;

      cout << "row number " << i+1 << endl;

      for (j=0; j<rs*dof*dof; j++)
	{
	  cout << val[p+j] << " ";
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
