// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.6  2001/08/16 16:53:23  prkipfer
| improved type safety for template parameter
|
| Revision 1.5  2001/06/19 16:32:27  prkipfer
| fixed typos and Linux compiler target
|
| Revision 1.4  2001/06/15 08:42:11  prkipfer
| resolved precision and sign flaw
|
| Revision 1.3  2001/02/13 10:59:14  prkipfer
| exchanged fixed numeric tolerances by machine dependent ones
|
| Revision 1.2  2001/01/09 11:14:06  prkipfer
| precision epsilon is now static
|
| Revision 1.1  2001/01/03 15:30:10  prkipfer
| introduced intersection, linear system solver and minimizer
|
|
|
+---------------------------------------------------------------------*/

#ifdef LINUX
#include <cassert>
#endif

template <class T>
const T
GoMinimize1D<T>::EPSILON = std::numeric_limits<T>::epsilon(); // max. 1e-08

template <class T>
GoMinimize1D<T>::GoMinimize1D (Function oF, int uiMaxLevel, int uiMaxBracket, void* pvUserData)
{
  assert( oF );
  minimizeFunction_ = oF;
  maxLevel_ = uiMaxLevel;
  maxBracket_ = uiMaxBracket;
  userData_ = pvUserData;
}

template <class T>
void 
GoMinimize1D<T>::getMinimum (T fT0, T fT1, T fTInitial, T& rfTMin, T& rfFMin)
{
  assert( fT0 <= fTInitial && fTInitial <= fT1 );

  tMin_ = std::numeric_limits<T>::max();
  fMin_ = std::numeric_limits<T>::max();

  T fF0 = minimizeFunction_(fT0,userData_);
  T fFInitial = minimizeFunction_(fTInitial,userData_);
  T fF1 = minimizeFunction_(fT1,userData_);

  getMinimum(fT0,fF0,fTInitial,fFInitial,fT1,fF1,maxLevel_);

  rfTMin = tMin_;
  rfFMin = fMin_;
}

template <class T>
void 
GoMinimize1D<T>::getMinimum (T fT0, T fF0, T fTm,
			     T fFm, T fT1, T fF1, int uiLevel)
{
  if ( fF0 < fMin_ ) {
    tMin_ = fT0;
    fMin_ = fF0;
  }

  if ( fFm < fMin_ ) {
    tMin_ = fTm;
    fMin_ = fFm;
  }

  if ( fF1 < fMin_ ) {
    tMin_ = fT1;
    fMin_ = fF1;
  }

  if ( uiLevel-- == 0 )
    return;

  if ( (fT1 - fTm)*(fF0 - fFm) > (fTm - fT0)*(fFm - fF1) ) {
    // quadratic fit has positive second derivative at midpoint

    if ( fF1 > fF0 ) {
      if ( fFm >= fF0 ) {
	// increasing, repeat on [t0,tm]
	getMinimum(fT0,fF0,fTm,fFm,uiLevel);
      }
      else {
	// not monotonic, have a bracket
	getBracketedMinimum(fT0,fF0,fTm,fFm,fT1,fF1,uiLevel);
      }
    }
    else if ( fF1 < fF0 ) {
      if ( fFm >= fF1 ) {
	// decreasing, repeat on [tm,t1]
	getMinimum(fTm,fFm,fT1,fF1,uiLevel);
      }
      else {
	// not monotonic, have a bracket
	getBracketedMinimum(fT0,fF0,fTm,fFm,fT1,fF1,uiLevel);
      }
    }
    else {
      // constant, repeat on [t0,tm] and [tm,t1]
      getMinimum(fT0,fF0,fTm,fFm,uiLevel);
      getMinimum(fTm,fFm,fT1,fF1,uiLevel);
    }
  }
  else {
    // quadratic fit has nonpositive second derivative at midpoint

    if ( fF1 > fF0 ) {
      // repeat on [t0,tm]
      getMinimum(fT0,fF0,fTm,fFm,uiLevel);
    }
    else if ( fF1 < fF0 ) {
      // repeat on [tm,t1]
      getMinimum(fTm,fFm,fT1,fF1,uiLevel);
    }
    else {
      // repeat on [t0,tm] and [tm,t1]
      getMinimum(fT0,fF0,fTm,fFm,uiLevel);
      getMinimum(fTm,fFm,fT1,fF1,uiLevel);
    }
  }
}

template <class T>
void 
GoMinimize1D<T>::getMinimum (T fT0, T fF0, T fT1, T fF1, int uiLevel)
{
  if ( fF0 < fMin_ ) {
    tMin_ = fT0;
    fMin_ = fF0;
  }

  if ( fF1 < fMin_ ) {
    tMin_ = fT1;
    fMin_ = fF1;
  }

  if ( uiLevel-- == 0 )
    return;

  T fTm = T(0.5)*(fT0+fT1);
  T fFm = minimizeFunction_(fTm,userData_);

  if ( fF0 - T(2.0)*fFm + fF1 > T(0.0) ) {
    // quadratic fit has positive second derivative at midpoint

    if ( fF1 > fF0 ) {
      if ( fFm >= fF0 ) {
	// increasing, repeat on [t0,tm]
	getMinimum(fT0,fF0,fTm,fFm,uiLevel);
      }
      else {
	// not monotonic, have a bracket
	getBracketedMinimum(fT0,fF0,fTm,fFm,fT1,fF1,uiLevel);
      }
    }
    else if ( fF1 < fF0 ) {
      if ( fFm >= fF1 ) {
	// decreasing, repeat on [tm,t1]
	getMinimum(fTm,fFm,fT1,fF1,uiLevel);
      }
      else {
	// not monotonic, have a bracket
	getBracketedMinimum(fT0,fF0,fTm,fFm,fT1,fF1,uiLevel);
      }
    }
    else {
      // constant, repeat on [t0,tm] and [tm,t1]
      getMinimum(fT0,fF0,fTm,fFm,uiLevel);
      getMinimum(fTm,fFm,fT1,fF1,uiLevel);
    }
  }
  else {
    // quadratic fit has nonpositive second derivative at midpoint

    if ( fF1 > fF0 ) {
      // repeat on [t0,tm]
      getMinimum(fT0,fF0,fTm,fFm,uiLevel);
    }
    else if ( fF1 < fF0 ) {
      // repeat on [tm,t1]
      getMinimum(fTm,fFm,fT1,fF1,uiLevel);
    }
    else {
      // repeat on [t0,tm] and [tm,t1]
      getMinimum(fT0,fF0,fTm,fFm,uiLevel);
      getMinimum(fTm,fFm,fT1,fF1,uiLevel);
    }
  }
}

template <class T>
void 
GoMinimize1D<T>::getBracketedMinimum (T fT0, T fF0,
				      T fTm, T fFm, 
				      T fT1, T fF1, int uiLevel)
{
  for (int uiI = 0; uiI < maxBracket_; uiI++) {
    // update minimum value
    if ( fFm < fMin_ ) {
      tMin_ = fTm;
      fMin_ = fFm;
    }

    // test for convergence
    const T fEps = GoMinimize1D<T>::EPSILON, fTol = GoMinimize1D<T>::EPSILON; // max 1e-8 und 1e-4
    if ( GbMath<T>::Abs(fT1-fT0) <= T(2.0)*fTol*GbMath<T>::Abs(fTm) + fEps )
      break;

    // compute vertex of interpolating parabola
    T fDT0 = fT0 - fTm, fDT1 = fT1 - fTm;
    T fDF0 = fF0 - fFm, fDF1 = fF1 - fFm;
    T fTmp0 = fDT0*fDF1, fTmp1 = fDT1*fDF0;
    T fDenom = fTmp1 - fTmp0;
    if ( GbMath<T>::Abs(fDenom) < fEps )
      return;

    T fTv = fTm + T(0.5)*(fDT1*fTmp1-fDT0*fTmp0)/fDenom;
    assert( fT0 <= fTv && fTv <= fT1 );
    T fFv = minimizeFunction_(fTv,userData_);

    if ( fTv < fTm ) {
      if ( fFv < fFm ) {
	fT1 = fTm;
	fF1 = fFm;
	fTm = fTv;
	fFm = fFv;
      }
      else {
	fT0 = fTv;
	fF0 = fFv;
      }
    }
    else if ( fTv > fTm ) {
      if ( fFv < fFm ) {
	fT0 = fTm;
	fF0 = fFm;
	fTm = fTv;
	fFm = fFv;
      }
      else {
	fT1 = fTv;
	fF1 = fFv;
      }
    }
    else {
      // vertex of parabola is already at middle sample point
      getMinimum(fT0,fF0,fTm,fFm,uiLevel);
      getMinimum(fTm,fFm,fT1,fF1,uiLevel);
    }
  }
}
