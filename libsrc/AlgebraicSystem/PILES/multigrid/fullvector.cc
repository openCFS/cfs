#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream.h>
#include <fstream.h>

#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{ 

RealVector :: RealVector(Integer asize, Integer anumrhs, Integer adof, Boolean mem)
  : BaseVector()
{
#ifdef TRACE
  (*trace) << "entering RealVector::RealVector" << endl;
#endif
  
  size   = asize;
  numrhs = anumrhs;
  dof    = adof;
  length = size*dof;

  if (mem)
    {
      val  = new Double[length*numrhs];

      Integer i;
      
      for (i=0; i<length; i++)
	{
	  val[i] = 0;
	}

      outback = FALSE;
    }
  else
    {
      val     = NULL;
      outback = TRUE;
    }

 
}
  
RealVector :: ~RealVector()
{
#ifdef TRACE
  (*trace) << "entering RealVector::~RealVector" << endl;
#endif

  if (val != NULL && !outback)
    {
      delete [] val;
    }
}

void RealVector :: Add(DenseVector & a, BaseVector &vec1, DenseVector & b, BaseVector &vec2)
{
  Integer i;

  RealVector & u   = (RealVector &) vec1;
  RealVector & v   = (RealVector &) vec2;
  DenseVector & g = (DenseVector &) a;
  DenseVector & h = (DenseVector &) b;

  for (i=0; i<length; i++)
    {
#ifdef SIMRHS
      Integer j;

      for (j=0; j<numrhs; j++)
	{
	  val[i*numrhs+j] = g.Get(j+1)*u.Get(i+1,j+1)+h.Get(j+1)*v.Get(i+1,j+1);
	}
#else
      val[i] = g.Get(1)*u.Get(i+1)+h.Get(1)*v.Get(i+1);
#endif
    }
}

void RealVector :: Add(DenseVector & a, BaseVector &vec)
{
  Integer i;

  RealVector & u  = (RealVector &) vec;
  DenseVector & g = (DenseVector &) a;

  for (i=0; i<length; i++)
    {
#ifdef SIMRHS
      Integer j;

      for (j=0; j<numrhs; j++)
	{
	  val[i*numrhs+j] += g.Get(j+1)*u.Get(i+1,j+1);
	}
#else
      val[i] += g.Get(1)*u.Get(i+1);
#endif
    }
}
  
void RealVector :: Add(BaseVector &vec)
{
  Integer i;

  RealVector & u = (RealVector &) vec;

  for (i=0; i<length; i++)
    {
#ifdef SIMRHS
      Integer j;

      for (j=0; j<numrhs; j++)
	{
	  val[i*numrhs+j] += u.Get(i+1,j+1);
	}
#else
     val[i] += u.Get(i+1); 
#endif
    }
}

DenseVector * RealVector :: Inner(BaseVector &vec) const
{
  Integer i;

  RealVector & u = (RealVector &) vec;
  
  DenseVector * sum = new DenseVector(numrhs);

  for (i=0; i<length; i++)
    {
#ifdef SIMRHS
      Integer j;

      for (j=0; j<numrhs; j++)
	{
	  sum.Elem(j+1) += val[i*numrhs+j]*u.Get(i+1,j+1);
	}
#else
      sum->Elem(1) += val[i]*u.Get(i+1);
#endif
    }
  
  return sum;
}

DenseVector & RealVector :: L2Norm() const
{
  Integer i,j;

  DenseVector sum(numrhs);

  for (i=0; i<length; i++)
    {
#ifdef SIMRHS
      for (j=0; j<numrhs; j++)
	{
	  sum.Elem(j+1) += val[i*numrhs+j]*val[i*numrhs+j];
	}
#else
      sum.Elem(1) += val[i]*val[i];
#endif
    }

  sum.Sqrt();

  return sum;
}

DenseVector & RealVector :: Max() const
{

}

void RealVector :: Rand()
{
  Integer i,j;

  srand48(1);
  
  for (i=0; i<length; i++)
    {
#ifdef SIMRHS
      for (j=0; j<numrhs; j++)
	{
	  val[i*numrhs+j] = drand48();
	}
#else
      val[i] = drand48();
#endif
    }
}

void RealVector :: Print() const
{
  Integer i;

  for (i=0; i<length; i++)
    {
      cout << val[i] << " ";
    }
  cout << endl;
}

void RealVector :: Assemble(Double * v, Integer * p, Integer elemsize)
{
  Integer i;

  for (i=0; i<elemsize; i++)
    {
      val[p[i]-1] += v[i];
    }
}

void RealVector :: Init()
{
  Integer i;

  for (i=0; i<length; i++)
    {
      val[i] = 0;
    }
}

void RealVector :: InitOne()
{
  Integer i;

  for (i=0; i<length; i++)
    {
      val[i] = 1;
    }
}



ComplexVector :: ComplexVector(Integer asize, Integer anumrhs, Integer adof, Boolean mem)
  : BaseVector()
{
#ifdef TRACE
  (*trace) << "entering ComplexVector::ComplexVector" << endl;
#endif

  size   = asize;
  numrhs = anumrhs;
  dof    = adof;
  length = size*dof;

  if (mem)
    {
      val  = new Double[2*length*numrhs];
      
      Integer i;
      
      for (i=0; i<2*length; i++)
	{
	  val[i] = 0;
	}
      
      outback = FALSE;
    }
  else
    {
      val = NULL;
      outback = TRUE;
    }
}
  
ComplexVector :: ~ComplexVector()
{
#ifdef TRACE
  (*trace) << "entering ComplexVector::~ComplexVector" << endl;
#endif
  
  if (val != NULL & !outback)
    {
      delete [] val;
    }
}

void ComplexVector :: Add(DenseVector & a, BaseVector &vec1, DenseVector & b, BaseVector &vec2)
{
  Integer i,j;

  ComplexVector & u = (ComplexVector &) vec1;
  ComplexVector & v = (ComplexVector &) vec2;
  DenseVector & g  = (DenseVector &) a;
  DenseVector & h  = (DenseVector &) b;

  Double re1,re2,im1,im2;

  j = 1;

  for (i=0; i<length; i++)
    {
      re1 = u.Get(j);
      im1 = u.Get(j+1);
      re2 = v.Get(j);
      im2 = v.Get(j+1);

      val[j-1] = g.Get(1)*re1-g.Get(2)*im1+h.Get(1)*re2-h.Get(2)*im2;
      val[j]   = g.Get(2)*re1+g.Get(1)*im1+h.Get(2)*re2+h.Get(1)*im2;

      j += 2;
    }
}

void ComplexVector :: Add(DenseVector & a, BaseVector &vec)
{
  Integer i,j;

  ComplexVector & u = (ComplexVector &) vec;
  DenseVector & g  = (DenseVector &) a;

  Double re,im;

  j = 1;

  for (i=0; i<length; i++)
    {
      re = u.Get(j);
      im = u.Get(j+1);

      val[j-1] += g.Get(1)*re-g.Get(2)*im;
      val[j]   += g.Get(2)*re+g.Get(1)*im;

      j += 2;
    }
}
  
void ComplexVector :: Add(BaseVector &vec)
{
  Integer i,j;

  ComplexVector & u = (ComplexVector &) vec;

  j = 1;

  for (i=0; i<length; i++)
    {
      val[j-1] += u.Get(j); 
      val[j]   += u.Get(j+1);
      
      j += 2;
    }
}

DenseVector * ComplexVector :: Inner(BaseVector &vec) const
{
  Integer i, j;

  ComplexVector & u = (ComplexVector &) vec;
  DenseVector * sum = new DenseVector(1);

  j = 1;

  for (i=0; i<size; i++)
    {
      sum->Elem(1) += val[j-1]*u.Get(j) - val[j]*u.Get(j+1);
      sum->Elem(2) += val[j-1]*u.Get(j+1) + val[j]*u.Get(j);

      j += 2;
    }

  return sum;
}

DenseVector & ComplexVector :: L2Norm() const
{
  Integer i,j;

  DenseVector sum(1);
  
  j = 0;

  for (i=0; i<length; i++)
    {
      sum.Elem(1) += val[j]*val[j] + val[j+1]*val[j+1];

      j += 2;
    }

  sum.Sqrt();

  return sum;
}

DenseVector & ComplexVector :: Max() const
{

}

void ComplexVector :: Rand()
{
  Integer i,j;

  srand48(1);
  
  j = 0;

  for (i=0; i<length; i++)
    {
      val[j]   = drand48();
      val[j+1] = drand48();

      j += 2;
    }
}

void ComplexVector :: Print() const
{
  Integer i,j;

  j = 0;

  for (i=0; i<length; i++)
    {
      (*test) << val[j] << " " << val[j+1] << ",";

      j += 2;
    }
  (*test) << endl;
}

void ComplexVector :: Assemble(Double * v, Integer * p, Integer elemsize)
{
  Integer i,j;

  for (i=0; i<elemsize; i++)
    {
      j = 2*p[i]-1;

      val[j-1] += v[i];
      val[j]   += v[i+elemsize];
    }
}

void ComplexVector :: Init()
{
  Integer i;

  for (i=0; i<2*length; i++)
    {
      val[i] = 0;
    }
}

}
