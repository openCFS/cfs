#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream.h>
#include <fstream.h>
#include <math.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

StandardCoarse :: StandardCoarse(BaseMatrix * amat, Integer amaxnh)
  : BaseCoarse(amaxnh)
{
#ifdef TRACE
  (*trace) << "entering StandardCoarse::StandardCoarse" << endl;
#endif
  
  Integer * pos, * start;

  mat = amat;
  fsize_sys = mat->GetSize();
  fnne_sys  = mat->GetNNE();
  pos       = mat->GetPosPointer();
  start     = mat->GetStartPointer();

  topology = new StandardTopology(fsize_sys, fnne_sys, pos, start);
}
  
StandardCoarse :: ~StandardCoarse()
{
#ifdef TRACE
  (*trace) << "entering StandardCoarse::~StandardCoarse" << endl;
#endif

  if (topology != NULL)
    {
      delete topology;
      topology = NULL;
    }
}

void StandardCoarse :: SetNeighbour(Double alpha, Double epsmat)
{
#ifdef TRACE
  (*trace) << "entering StandardCoarse::SetNeighbour" << endl;
#endif

  Integer i,j,k,p,rs,maxelem;
  Double v,minval,diaval;
  
  RScalarMatrix & a = (RScalarMatrix &) *mat;

  for (i=1; i<=fsize_sys; i++)
    {
      rs = a.GetRowSize(i);
      
      if (rs == 0)
	{
	  cout << "... matrix is singular exit in CoupledField::StandardCoarse::SetNeighbour" << endl;
	  cout << "... row number " << i << endl;
	  exit(1);
	}
      
      minval = a.MinRowVal(i);
      diaval = a.GetDiag(i);
      
      for (j=2; j<=rs; j++)
	{
	  v = a.Get(i,j);
	  p = a.GetMatrixPos(i,j);
	  
	  /// different coarsening functions!!!
	  
	  if (fabs(v) >= diaval*epsmat && v <= minval*alpha)
	    {
	      topology->UpdateS1(i,p);  /// p into S1(i)
	    }
	}
    }
      
  /// sort the local sets?????
  
  /// S(i,T) = s2
  
  startpoint = 0;
  maxelem    = 0;
  
  for (i=1; i<=fsize_sys; i++)
    {
      rs = mat->GetRowSize(i);
      
      for (j=2; j<=rs; j++)
	{
	  k = topology->GetNHNode(i,j);
	  
	  if (topology->MemberS1(i,k)) /// k member of S1(i)
	    {
	      topology->UpdateS2(i,k); /// k into S2(i)
	    }
	}  
      
      if ((topology->GetS2Size(i) > maxelem) && (topology->GetNHSize(i) != 0))
	{
	  startpoint = i;
	  maxelem    = topology->GetS2Size(i);
	}
    } 

  Integer numdir;

  numdir = a.GetNumDirichlet();

  for (i=1; i<=numdir; i++)
    {
      p = a.GetIndDirichlet(i);
      topology->SetDirichlet(p);
    }
}

void StandardCoarse :: CalcTopology()
{
#ifdef TRACE
  (*trace) << "entering StandardCoarse::CalcTopology" << endl;
#endif

}

void StandardCoarse :: SetSize()
{
#ifdef TRACE
  (*trace) << "entering StandardCoarse::SetSize" << endl;
#endif

  csize_sys = topology->GetCoarseSize();
  csize_aux = csize_sys;
  cnne_sys  = topology->GetCoarseNNE();
  cnne_aux  = cnne_sys;
}
}
