// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

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

#include <cstdlib>
#include <iostream>
#include <string>
#include <cmath>

namespace grd {
/*
#define DATA_TYPE_CHAR       0
#define DATA_TYPE_SHORT      1
#define DATA_TYPE_FLOAT      2
#define DATA_TYPE_UNIFORM    3
#define DATA_TYPE_SCALED     4
#define DATA_TYPE_STRUCTURED 5
*/

  enum {
    DATA_TYPE_CHAR,
    DATA_TYPE_SHORT,
    DATA_TYPE_FLOAT,
    DATA_TYPE_UNIFORM,
    DATA_TYPE_SCALED,
    DATA_TYPE_STRUCTURED
  };


// *-
// for grid level 0
/*
#define INIT_UNIT_CUBE_TETRA  1
#define INIT_UNIT_TETRA       2
#define INIT_UNIT_CUBE_HEXA   3
*/

  enum {
    INIT_UNIT_CUBE_TETRA,
    INIT_UNIT_TETRA,
    INIT_UNIT_CUBE_HEXA
  };

// *-
// constants
const double ONE_THIRD = 1.0 / 3.0;
const double ONE_SIXTH = 1.0 / 6.0;

// Info Macros
inline void Info(const char* infoText)
{
  std::cerr << "Info: " << infoText << std::endl;
}

inline void Info(const char* infoText,const char* comment)
{
  std::cerr << "Info: "
            << infoText
            << std::endl
            << comment
            << std::endl;
}

inline void Info(const std::string& infoText)
{
  std::cerr << "Info:"
            << infoText
            << std::endl;
}

inline void Info(const std::string& infoText,const std::string& comment)
{
  std::cerr << "Info: "
            << infoText
            << std::endl
            << comment
            << std::endl;
}
// WARN Macros
inline void WARN(const char* warningText)
{
  std::cerr << "Grid WARNING: "
       << warningText
       << std::endl;
}

inline void WARN(const char* warningText,const char* comment)
{
  std::cerr << "Grid WARNING: "
       << warningText
       << std::endl
       << comment
       << std::endl;
}

inline void WARN(const std::string& warningText)
{
  std::cerr << "Grid WARNING: "
       << warningText
       << std::endl;
}

inline void WARN(const std::string& warningText,const std::string& comment)
{
  std::cerr << "Grid WARNING: "
       << warningText
       << std::endl
       << comment
       << std::endl;
}


// Error Macros
inline void Error(const char* errorText)
{
  std::cerr << "Gird ERROR: "
            << errorText
            << std::endl;
  exit(1);
}

inline void Error(const char* errorText, const char* comment)
{
  std::cerr << "Gird ERROR: "
            << errorText
            << std::endl
            << comment
            << std::endl;
  exit(1);
}
inline void Error(const std::string& errorText)
{
  std::cerr << "Gird ERROR: "
            << errorText
            << std::endl;
  exit(1);
}

inline void Error(const std::string& errorText, const std::string& comment)
{
  std::cerr << "Gird ERROR: "
            << errorText
            << std::endl
            << comment
            << std::endl;
  exit(1);
}

//=============================================================================
// Auxiliary math function
//=============================================================================
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


inline float
length(float* v)
{
  float res = 0.0;
  for (int i = 0; i < 3; i++)
  {
    res += square(v[i]);
  }

  return sqrt(res);
}

inline double
length(const double* v)
{
  double res = 0.0;
  for (int i = 0; i < 3; i++)
  {
    res += square(v[i]);
  }

  return sqrt(res);
}

inline double
lengthSquare(double v0,double v1,double v2)
{
  double res = 0.0;
  res = v0*v0+v1*v1+v2*v2;

  return res;
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
innerProduct(const double* n1,const double* n2)
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
