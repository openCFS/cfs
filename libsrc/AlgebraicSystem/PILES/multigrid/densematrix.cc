#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <iostream.h>
#include <fstream.h>

#include "environment.hh"
#include <general.hh>
#include "multigrid.hh"

namespace CoupledField
{

DenseMatrix :: DenseMatrix(Integer asize, Integer awidth, Boolean amem)
{
#ifdef TRACE
  (*trace) << "entering DenseMatrix::DenseMatrix" << endl;
#endif

  size = asize;
  width= awidth;
  mem  = amem;

  if (mem == TRUE)
    {
      val = new Double[size*width];
      fac = new Double[size*width];
      u   = new Double[size];
      v   = new Double[size];
      y   = new Double[size];
      
      Integer i;
      
      for (i=0; i<size*width; i++)
	{
	  val[i] = 0;
	  fac[i] = 0;
	}
    }
  else
    {
      val = NULL;
      fac = NULL;
      u   = NULL;
      v   = NULL;
      y   = NULL;
    }
}
  
DenseMatrix :: ~DenseMatrix()
{
#ifdef TRACE
  (*trace) << "entering DenseMatrix::~DenseMatrix" << endl;
#endif

  if (mem && val != NULL)
    {
      delete [] val;
      delete [] fac;
      delete [] u;
      delete [] v;
      delete [] y;
    }
}

void DenseMatrix :: InitMatrix()
{
  Integer i,j;

  for (i=0; i<size; i++)
    {
      for (j=0; j<size; j++)
	{
	  if (i == j)
	    {
	      val[i*size+j] = 5;
	    }
	  else
	    {
	      val[i*size+j] = -1;
	    }
	}
    }
}

void DenseMatrix :: Mult(DenseMatrix * a, DenseMatrix * b)
{
#ifdef TRACE
  (*trace) << "entering DenseMatrix::Mult" << endl;
#endif

  Integer i,j,k;
  Double sum;

  for (i=0; i<size; i++)
    {
      for (j=0; j<size; j++)
	{
	  sum = 0;

	  for (k=0; k<size; k++)
	    {
	      sum += (a->Get(i+1,k+1)*b->Get(k+1,j+1));
	    }

	  val[i*size+j] = sum;
	}
    }
}

void DenseMatrix :: MultT(DenseMatrix * a, DenseMatrix * b)
{
#ifdef TRACE
  (*trace) << "entering DenseMatrix::MultT" << endl;
#endif

  Integer i,j,k;
  Double sum;

  for (i=0; i<size; i++)
    {
      for (j=0; j<size; j++)
	{
	  sum = 0;

	  for (k=0; k<size; k++)
	    {
	      sum += (a->Get(k+1,i+1)*b->Get(k+1,j+1));
	    }

	  val[i*size+j] = sum;
	}
    }

}

Double DenseMatrix :: MaxNorm() const
{
#ifdef TRACE
  (*trace) << "entering DenseMatrix::~DenseMatrix" << endl;
#endif

  Integer i;
  Double max;

  max = val[0];

  for (i=1; i<size; i++)
    {
      if (max < val[i*size+i])
	{
	  max = val[i*size+i];
	}
    }

  return max;
}

DenseMatrix & DenseMatrix :: operator=(Double * x)
{
  val = x;
  
  return *this;
}

void DenseMatrix :: Inverse(DenseMatrix * a)
{
#ifdef TRACE
  (*trace) << "entering DenseMatrix::Inverse" << endl;
#endif

  Double det;

  switch (size)
    {
    case 1:
      val[0] = 1./a->Get(1,1);
      break;

    case 2:
      det = a->Get(1,1)*a->Get(2,2) - a->Get(1,2)*a->Get(2,1);
      val[0] = a->Get(2,2)/det;
      val[1] = -a->Get(1,2)/det;
      val[2] = -a->Get(2,1)/det;
      val[3] = a->Get(1,1)/det;
      break;

    case 3:
      det = a->Get(1,1)*a->Get(2,2)*a->Get(3,3)+
	a->Get(1,2)*a->Get(2,3)*a->Get(3,1)+
	a->Get(1,3)*a->Get(2,1)*a->Get(3,2)-
	a->Get(1,1)*a->Get(2,3)*a->Get(3,2)-
	a->Get(1,2)*a->Get(2,1)*a->Get(3,3)-
	a->Get(1,3)*a->Get(2,2)*a->Get(3,1);

      val[0] = (a->Get(2,2)*a->Get(3,3) - a->Get(2,3)*a->Get(3,2))/det;
      val[1] = (a->Get(1,3)*a->Get(3,2) - a->Get(1,2)*a->Get(3,3))/det;
      val[2] = (a->Get(1,2)*a->Get(2,3) - a->Get(1,3)*a->Get(2,2))/det;

      val[3] = (a->Get(2,3)*a->Get(3,1) - a->Get(2,1)*a->Get(3,3))/det;
      val[4] = (a->Get(1,1)*a->Get(3,3) - a->Get(1,3)*a->Get(3,1))/det;
      val[5] = (a->Get(1,3)*a->Get(2,1) - a->Get(1,1)*a->Get(2,3))/det;

      val[6] = (a->Get(1,2)*a->Get(3,2) - a->Get(2,2)*a->Get(3,1))/det;
      val[7] = (a->Get(1,2)*a->Get(3,1) - a->Get(1,1)*a->Get(3,2))/det;
      val[8] = (a->Get(1,1)*a->Get(2,2) - a->Get(1,2)*a->Get(2,1))/det;
      break;

    default:

      Integer i,j,k;
      Double h;

      /// factorization

      for (i=0; i<size; i++)
	{
	  for (j=i; j<size; j++)
	    {
	      h = 0;
	      
	      for (k=0; k<=i-1; k++)
		{
		  h += fac[i*size+k]*fac[k*size+j];
		}
	      
	      fac[i*size+j] = a->Get(i+1,j+1) - h;
	    }
	  
	  for (j=i+1; j<size; j++)
	    {
	      h = 0;
	      
	      for (k=0; k<=i-1; k++)
		{
		  h += fac[j*size+k]*fac[k*size+i];
		}
	      
	      fac[j*size+i] = (a->Get(i+1,j+1) - h)/fac[i*size+i];
	    }
	}

      /// calc the inverse matrix

      for (k=0; k<size; k++)
	{
	  v[k] = 1;
	  
	  y[0] = v[0];
	  
	  for (i=0; i<size; i++)
	    {
	      h = 0;
	      
	      for (j=0; j<=i-1; j++)
		{
		  h += fac[i*size+j]*y[j];
		}
	      
	      y[i] = v[i] - h;
	    }
	  
	  u[size-1] = y[size-1]/fac[size*size-1];
	  
	  for (i=size-2; i>=0; i--)
	    {
	      h = 0;
	      
	      for (j=i+1; j<size; j++)
		{
		  h += fac[i*size+j]*u[j];
		}
	      
	      u[i] = (y[i] - h)/fac[i*size+i];
	    }
	  
	  for (j=0; j<size; j++)
	    {
	      val[j*size+k] = u[j];
	      
	      v[j] = 0;
	      u[j] = 0;
	    }
	}
    }
}

void DenseMatrix :: Print() const
{
#ifdef TRACE
  (*trace) << "entering DenseMatrix::Print" << endl;
#endif

  Integer i,j;

  for (i=0; i<size; i++)
    {
      for (j=0; j<width; j++)
	{
	  (*cla) << val[i*size+j] << " ";
	}
      (*cla) << endl;
    }
}
}
