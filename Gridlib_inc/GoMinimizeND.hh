/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GOMINIMIZEND_HH
#define  GOMINIMIZEND_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbMath.hh"
#include "GoMinimize1D.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GoMinimizeND
{
public:
  typedef T (*Function)(const T*,void*);

  GoMinimizeND (int iDimensions, Function oF, int uiMaxLevel,
		int uiMaxBracket, int uiMaxIterations,
		void* pvUserData = NULL);

  ~GoMinimizeND ();

  INLINE int& maxLevel ();
  INLINE int& maxBracket ();
  INLINE void*& userData ();

  // find minimum on Cartesian-product domain
  void getMinimum (const T* afT0, const T* afT1,
		   const T* afTInitial, T* afTMin, T& rfFMin);

protected:
  int dimensions_;
  Function minFunction_;
  int maxIterations_;
  void* userData_;
  GoMinimize1D<T> oneDMinimizer_;
  T* directionStorage_;
  T** direction_;
  T* dConj_;
  T* dCurr_;
  T* tSave_;
  T* tCurr_;
  T fCurr_;
  T* lineArg_;

  void computeDomain (const T* afT0, const T* afT1, T& rfL0, T& rfL1);

  static T lineFunction (T fT, void* pvUserData);
};


#ifdef __GNUG__

#include "GoMinimizeND.in"
#include "GoMinimizeND.T"

#else

#pragma instantiate GoMinimizeND<float>
#pragma instantiate GoMinimizeND<double>

#ifndef OUTLINE
#include "GoMinimizeND.in"
#endif

#endif  // g++

#endif  // GOMINIMIZEND_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.2  2001/06/15 08:42:41  prkipfer
| resolved precision and sign flaw
|
| Revision 1.1  2001/01/03 15:30:12  prkipfer
| introduced intersection, linear system solver and minimizer
|
|
+---------------------------------------------------------------------*/
