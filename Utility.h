//  -*- C++ -*-
/***************************************************************************
    File        : Utility.h
    Description :

 ---------------------------------------------------------------------------
    Begin       : Sun Jan 27 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#ifndef UTILITY_H
#define UTILITY_H

#include <cmath>

namespace grd {

#define DATA_TYPE_CHAR       0
#define DATA_TYPE_SHORT      1
#define DATA_TYPE_FLOAT      2
#define DATA_TYPE_UNIFORM    3
#define DATA_TYPE_SCALED     4
#define DATA_TYPE_STRUCTURED 5

// *-
// for grid level 0
#define INIT_UNIT_CUBE_TETRA  1
#define INIT_UNIT_TETRA       2
#define INIT_UNIT_CUBE_HEXA   3


inline int
square(int a)
{
  return (a*a);
}


inline float
square(float a)
{
  return (a*a);
}


inline double
square(double a)
{
  return (a*a);
}


inline int
sym(int n, int m)
{
    int m2 = ((m << 1) - 2);
    int tp = ((n >= 0) ? (n % m2) : ((m2 - (-n % m2)) % m2) );
    return ((tp < m) ? tp : (m2 - tp));
}

inline double
length(double* v)
{
  double res = 0.0;
  for (int i = 0; i < 3; i++)
  {
    res += square(v[i]);
  }

  return sqrt(res);
}


inline double
distance(double* n1,double* n2)
{
  double res = 0.0;
  for (int i = 0; i < 3; i++)
  {
    res += square(n1[i] - n2[i]);
  }
  return sqrt(res);
}


inline void
crossProduct(double* res,double* v1,double* v2)
{
  res[0] = v1[1]*v2[2] - v1[2]*v2[1];
  res[1] = v1[2]*v2[0] - v1[0]*v2[2];
  res[2] = v1[0]*v2[1] - v1[1]*v2[0];
}

inline double
innerProduct(double* n1, double* n2)
{
  double res = 0.0;
  for (int i = 0; i < 3; i++)
  {
    res += n1[i]*n2[i];
  }

  return res;
}

} // namespace grd
#endif // UTILITY_H