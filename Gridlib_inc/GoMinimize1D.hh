/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GOMINIMIZE1D_HH
#define  GOMINIMIZE1D_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GbTypes.hh"
#include "GbMath.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GoMinimize1D
{
public:
  typedef T (*Function)(T,void*);

  GoMinimize1D (Function oF, int uiMaxLevel, int uiMaxBracket, void* pvUserData = NULL);

  INLINE ~GoMinimize1D ();

  INLINE int& maxLevel ();
  INLINE int& maxBracket ();
  INLINE void*& userData ();

  void getMinimum (T fT0, T fT1, T fTInitial, T& rfTMin, T& rfFMin);

  static const T EPSILON;

protected:
  Function minimizeFunction_;
  int maxLevel_, maxBracket_;
  T tMin_, fMin_;
  void* userData_;

  void getMinimum (T fT0, T fF0, T fT1, T fF1, int uiLevel);

  void getMinimum (T fT0, T fF0, T fTm, T fFm, T fT1, T fF1, int uiLevel);

  void getBracketedMinimum (T fT0, T fF0, T fTm, T fFm, T fT1, T fF1, int uiLevel);
};

#ifdef __GNUG__

#include "GoMinimize1D.in"
#include "GoMinimize1D.T"

#else

#pragma instantiate GoMinimize1D<float>
#pragma instantiate GoMinimize1D<double>

#ifndef OUTLINE
#include "GoMinimize1D.in"
#endif

#endif  // g++

#endif  // GOMINIMIZE1D_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.3  2001/06/15 08:42:11  prkipfer
| resolved precision and sign flaw
|
| Revision 1.2  2001/01/09 11:14:06  prkipfer
| precision epsilon is now static
|
| Revision 1.1  2001/01/03 15:30:11  prkipfer
| introduced intersection, linear system solver and minimizer
|
|
+---------------------------------------------------------------------*/
