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
// DenseVector :: DenseVector()
// {
// #ifdef TRACE
//   (*trace) << "entering DenseVector::DenseVector" << endl;
// #endif

// }

// DenseVector :: ~DenseVector()
// {
// #ifdef TRACE
//   (*trace) << "entering DenseVector::~DenseVector" << endl;
// #endif
// }

//////////////////////////////// real ///////////////////////////////

DenseVector :: DenseVector(Integer asize)
{
#ifdef TRACE
  (*trace) << "entering RDenseVector::RDenseVector" << endl;
#endif

  size = asize;
  val  = new Double[size];

  Integer i;

  for (i=0; i<size; i++)
    {
      val[i] = 0;
    }
}

  
DenseVector :: ~DenseVector()
{
#ifdef TRACE
  (*trace) << "entering RDenseVector::~RDenseVector" << endl;
#endif

  delete val;
}

void DenseVector :: Sqrt()
{
  Integer i;

  for (i=0; i<size; i++)
    {
      val[i] = sqrt(val[i]);
    }
}

void DenseVector :: Abs()
{
  Integer i;

  for (i=0; i<size; i++)
    {
      val[i] = fabs(val[i]);
    }
}

DenseVector & DenseVector :: operator=(const Double src)
{
  Integer i;
  
  for (i=0; i<size; i++)   
    {
      val[i] = src;
    }
  
  return *this;
}

DenseVector & DenseVector :: operator=(const DenseVector &src)
{
  Integer i;
  
  for (i=1; i<=size; i++)
    {
      val[i-1] = src.Get(i);
    }

  return *this;
}

DenseVector & DenseVector :: operator*=(const Double x)
{
  Integer i;
  
  for (i=0; i<size; i++)   
    {
      val[i] *= x;
    }
  
  return * this;
}



DenseVector DenseVector :: operator-() const
{
  Integer i;

  DenseVector tmp(size);

  for (i=1; i<=size; i++)
    {
      tmp.Elem(i) = -val[i-1];
    }

 return tmp;
}  

void DenseVector :: Div(DenseVector &vec1, DenseVector &vec2)
{
  Integer i;

  for (i=1; i<=size; i++)
    {
      val[i-1] = vec1.Get(i)/vec2.Get(i);
    }
}

DenseVector DenseVector :: operator/(const DenseVector &src) const
{
  Integer i;

  DenseVector tmp(size);

  for (i=1; i<=size; i++)
    {
      tmp.Elem(i) = val[i-1]/src.Get(i);
    }

 return tmp;
} 

Boolean DenseVector :: operator<=(const DenseVector &x) const
{
  Boolean bb=TRUE;
  
  Integer i;

   for (i=1; i<=size; i++)
    {
      bb = ((val[i-1] <=  x.Get(i)) && bb);
    } 

  return bb;
  
}


Boolean DenseVector :: operator<=(const Double x) const
{
  Boolean bb=TRUE;

  Integer i;

   for (i=1; i<=size; i++)
    {
      bb = ((val[i-1] <=  x) && bb);
    } 

  return bb;
  
}

void DenseVector :: Print() const
{
  Integer i;

  for (i=0; i<size; i++)
    {
      cout << val[i] << " ";
    }
  cout << endl;
}
}
