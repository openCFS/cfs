#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>

#include "environment.hh"
#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

BaseTopology :: BaseTopology(Integer asize, Integer anne, Integer * pos, Integer * start)
{
#ifdef TRACE
  (*trace) << "entering BaseTopology::BaseTopology" << endl;
#endif
  
  size   = asize;
  nne    = anne;
  nh     = pos;
  snh    = start;
  csize  = 0;
  s1     = new Integer[nne];
  s2     = new Integer[nne];
  ss1    = new Integer[size];
  ss2    = new Integer[size];
  active = new Integer[size];
  cnn    = new Integer[size];
  offset = 30;

  Integer i,k;

  offset = 0;

  for (i=0; i<size; i++)
    {
      k = start[i+1] - start[i];

      if (k > offset)
	{
	  offset = k;
	}
    }

  offset *= 10;

  for (i=0; i<size; i++)
    {
      ss1[i]    = 0;
      ss2[i]    = 0;
      active[i] = 1;
      cnn[i]    = 0;
    }
#ifdef MEMTRACE
  double dmb;
  double imb;

  dmb = 0;
  imb = (4.*size+2*nne)*4/1e6;

  sumdmem += dmb;
  sumimem += imb;

  (*memtrace) << "+++ ALLOCATE MEMORY: double  BaseTopology     " << dmb << " MB" << endl;
  (*memtrace) << "+++ ALLOCATE MEMORY: integer BaseTopology     " << imb << " MB" << endl;
#endif

}

BaseTopology :: ~BaseTopology()
{
#ifdef TRACE
  (*trace) << "entering BaseTopology::~BaseTopology" << endl;
#endif
}

Boolean BaseTopology :: MemberS1(Integer p, Integer q)
{
#ifdef TRACE
  (*trace) << "entering BaseTopology::MemberS1" << endl;
#endif

  /// q member of S1(p)

  Integer i,rs;

  rs = ss1[p-1];

  for (i=0; i<rs; i++)
    {
      if (q == s1[snh[p-1]+i])
	{
	  return TRUE;
	}
    }

  return FALSE;
}

void BaseTopology :: UpdateS1(Integer p, Integer q)
{
#ifdef TRACE
  (*trace) << "entering BaseTopology::UpdateS1" << endl;
#endif
  /// insert q into S1(p)

  s1[snh[p-1]+ss1[p-1]] = q;
  ss1[p-1]++;
}

void BaseTopology :: UpdateS2(Integer p, Integer q)
{
  /// insert q into S2(p)

  s2[snh[p-1]+ss2[p-1]] = q;
  ss2[p-1]++;
}

Integer BaseTopology :: GetNextCoarseGridPoint(Integer p, Integer maxnh) const
{
#ifdef TRACE
  (*trace) << "entering BaseTopology::GetNextCoarseGridPoint" << endl;
#endif
  Integer i,j,m,n,n1,n2,c1,c2,max,maxsize;
  
  c2      = 0;
  max     = 0;
  m       = snh[p] - snh[p-1];
  maxsize = 0;
  
  for (i=1; i<m; i++)
    {
      n1    = nh[snh[p-1]+i];
      n     = snh[n1] - snh[n1-1];
      
      if (active[n1-1] == -1  ||  (active[n1-1] > 0 && active[n1-1] < 4 ) ) // active or defined nodes
	{
	  for (j=1; j<n; j++)
	    {
	      n2 = nh[snh[n1-1]+j];

	      if (cnn[n2-1] == 0  &&  (active[n2-1] >= 1 && active[n2-1] <= 2) )	// coarse nodes
		{
		  c1 = GetNumFine(n2);
		  c2 = ss2[n2-1];
		  
		  if ((c1+c2) > max)
		    {
		      max     = c1+c2;
		      maxsize = n2;
		      
		      if (max >= maxnh)
			{
			  return(maxsize);
			}
		    }
		}
	    }
	}
    }
  
  return(maxsize);
}

Integer BaseTopology :: GetNextStartPoint() const
{
#ifdef TRACE
  (*trace) << "entering BaseTopology::GetNextStartPoint" << endl;
#endif

  Integer i,p;

  for (i=1; i<=size; i++)
    {
      if (cnn[i-1] == 0)
	{
	  return i;
	}
    }

  return 0;
}

void BaseTopology :: SetFineGridPoint(Integer p)
{
#ifdef TRACE
  (*trace) << "entering BaseTopology::SetFineGridPoint" << endl;
#endif

  Integer i,rs;

  rs = ss2[p-1];

  for (i=0; i<rs; i++)
    {
      if (cnn[s2[snh[p-1]+i]-1] == 0)
	{
	  cnn[s2[snh[p-1]+i]-1] = -1;
	}
    }
}

void BaseTopology :: SetRSW(Integer * rsw)
{
#ifdef TRACE
  (*trace) << "entering BaseTopology::SetRSW" << endl;
#endif

  Integer i;

  for (i=0; i<size; i++)
    {
      if (cnn[i] > 0)
	{
	  rsw[i] = 1;
	}
      else
	{
	  rsw[i] = GetStrongNumCoarse(i+1);
	}
    }
}

Integer BaseTopology :: GetStrongNumCoarse(Integer p) const
{
#ifdef TRACE
  (*trace) << "entering BaseTopology::GetNumCoarse" << endl;
#endif

  Integer i,rs,c;

  rs = ss2[p-1];//snh[p]-snh[p-1];
  c  = 0;

  for (i=0; i<rs; i++)
    {
      if (cnn[s2[snh[p-1]+i]-1] > 0)
	{
	  c++;
	}
    }

  return c;
}

Integer BaseTopology :: GetNumCoarse(Integer p) const
{
#ifdef TRACE
  (*trace) << "entering BaseTopology::GetNumCoarse" << endl;
#endif

  Integer i,rs,c;

  rs = snh[p]-snh[p-1];
  c  = 0;

  for (i=0; i<rs; i++)
    {
      if (cnn[nh[snh[p-1]+i]-1] > 0)
	{
	  c++;
	}
    }

  return c;
}

Integer BaseTopology :: GetNumFine(Integer p) const
{
#ifdef TRACE
  (*trace) << "entering BaseTopology::GetNumFine" << endl;
#endif

  Integer i,rs,c;

  rs = snh[p]-snh[p-1];
  c  = 0;

  for (i=0; i<rs; i++)
    {
      if (cnn[nh[snh[p-1]+i]-1] < 0)
	{
	  c++;
	}
    }

  return c;
}

Integer BaseTopology :: MaxNumNeighbour() const
{
#ifdef TRACE
  (*trace) << "entering BaseTopology::MaxNumNeighbour" << endl;
#endif

  Integer i,maxnh;

  maxnh = 0;

  for (i=0; i<size; i++)
    {
      if ((snh[i+1]-snh[i]) > maxnh)
	{
	  maxnh = snh[i+1]-snh[i];
	}
    }

  return maxnh;
}

Integer BaseTopology ::GetStrongCoarseNeighbour(Integer p, Integer * cnh)
{
#ifdef TRACE
  (*trace) << "entering BaseTopology::GetStrongCoarseNeighbour" << endl;
#endif

  Integer i,rs,c;

  rs = ss2[p-1];//snh[p]-snh[p-1];
  c  = 0;

  for (i=0; i<rs; i++)
    {
      if (cnn[s2[snh[p-1]+i]-1] > 0)
	{
	  cnh[c] = s2[snh[p-1]+i];
	  c++;
	}
    }

  return c;
}

Integer BaseTopology ::GetCoarseNeighbour(Integer p, Integer * cnh)
{
#ifdef TRACE
  (*trace) << "entering BaseTopology::GetCoarseNeighbour" << endl;
#endif

  Integer i,rs,c;

  rs = snh[p]-snh[p-1];
  c  = 0;

  for (i=0; i<rs; i++)
    {
      if (cnn[nh[snh[p-1]+i]-1] > 0)
	{
	  cnh[c] = nh[snh[p-1]+i];
	  c++;
	}
    }

  return c;
}

Integer BaseTopology ::GetFineNeighbour(Integer p, Integer * fnh)
{
#ifdef TRACE
  (*trace) << "entering BaseTopology::GetFineNeighbour" << endl;
#endif

  Integer i,rs,c;

  rs = snh[p]-snh[p-1];
  c  = 0;

  for (i=0; i<rs; i++)
    {
      if (cnn[nh[snh[p-1]+i]-1] < 0)
	{
	  fnh[c] = nh[snh[p-1]+i];
	  c++;
	}
    }

  return c;
}

void BaseTopology :: Insert(Integer p, Integer q)
{
#ifdef TRACE
  (*trace) << "entering BaseTopology::Insert" << endl;
#endif

  Integer i,lower,upper;

  lower = (q-1)*offset;
  upper = (q-1)*offset+nrs[q-1];

  for (i=lower; i<upper; i++)
    {
      if (crs[i] == p)
	{
	  return;
	}
    }

  crs[upper] = p;
  nrs[q-1]++;

  upper      = (p-1)*offset+nrs[p-1];
  crs[upper] = q;
  nrs[p-1]++;
}

void BaseTopology :: Print() const
{
#ifdef TRACE
  (*trace) << "entering BaseTopology::Print" << endl;
#endif

  Integer i,j,rs;

  for (i=0; i<size; i++)
    {
      rs = snh[i+1]-snh[i];

      if (cnn[i] < 0)
	{
	  (*cla) << i+1 << " : " << "is fine grid node" << endl;
	}      
      else if (cnn[i] > 0)
	{
	  (*cla) << i+1 << " : " << "is coarse grid node " << cnn[i] << endl;
	}

      for (j=0; j<rs; j++)
	{
	  (*cla) << nh[snh[i]+j] << " ";
	}
      (*cla) << endl;

      rs = ss1[i];

      for (j=0; j<rs; j++)
	{
	  (*cla) << s1[snh[i]+j] << " ";
	}
      (*cla) << endl;

      rs = ss2[i];

      for (j=0; j<rs; j++)
	{
	  (*cla) << s2[snh[i]+j] << " ";
	}
      (*cla) << endl;
    }
}

////////////////////////////////////////////////////////////////////

LocalSet :: LocalSet(Integer asize, Integer adof)
{
#ifdef TRACE
  (*trace) << "entering Localset::LocalSet" << endl;
#endif
  
  size = asize;
  dof  = adof;

  v = new Double[size*dof*dof];
  p = new Integer[size];

  dm = new DenseMatrix(dof,dof);
}

LocalSet :: ~LocalSet()
{
#ifdef TRACE
  (*trace) << "entering LocalSet::~LocalSet" << endl;
#endif

  if (v != NULL)
    {
      delete [] v;
      delete [] p;
    }

  delete dm;
}

void LocalSet :: Init()
{
  Integer i;

  acts = 0;

  for (i=0; i<size*dof*dof; i++)
    {
      v[i] = 0;
    }

  for (i=0; i<size; i++)
    {
      p[i] = 0;
    }
}

void LocalSet :: Insert(Integer q)
{
#ifdef TRACE
  (*trace) << "entering LocalSet::Insert" << endl;
#endif

  Integer i;

  for (i=0; i<acts; i++)
    {
      if (p[i] == q)
	{
	  return;
	}
    }

  p[acts] = q;
  acts++;
}

void LocalSet :: Insert(Integer q, Double val)
{
#ifdef TRACE
  (*trace) << "entering LocalSet::Insert" << endl;
#endif

  Integer i;

  for (i=0; i<acts; i++)
    {
      if (p[i] == q)
	{
	  v[i] += val;
	  return;
	}
    }
}

void LocalSet :: Insert(Integer q, DenseMatrix * val)
{
#ifdef TRACE
  (*trace) << "entering LocalSet::Insert" << endl;
#endif

  Integer i,j,k,m,n;

  for (i=0; i<acts; i++)
    {
      if (p[i] == q)
	{
	  m = i*dof*dof;
	  n = 0;

	  for (j=1; j<=dof; j++)
	    {
	      for (k=1; k<=dof; k++)
		{
		  v[m+n] += val->Get(j,k);
		  n++;
		}
	    }

	  return;
	}
    }
}

void LocalSet :: Sort()
{
#ifdef TRACE
  (*trace) << "entering LocalSet::Sort" << endl;
#endif

  Integer i,j,q;

  for (i=0; i<acts; i++)
    {
      for (j=i+1; j<acts; j++)
	{
	  if (p[i] > p[j])
	    {
	      q    = p[i];
	      p[i] = p[j];
	      p[j] = q;
	    }
	}
    }
}

void LocalSet :: Print() const
{
#ifdef TRACE
  (*trace) << "entering LocalSet::Print" << endl;
#endif

  Integer i;

  (*cla) << "printing a local set" << endl;

  for (i=0; i<acts; i++)
    {
      (*cla) << p[i] << " ";
    }
  (*cla) << endl;

  for (i=0; i<acts*dof*dof; i++)
    {
      (*cla) << v[i] << " ";
    }
  (*cla) << endl;
}

}
