// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.5  2001/08/16 16:53:23  prkipfer
| improved type safety for template parameter
|
| Revision 1.4  2001/06/19 16:32:27  prkipfer
| fixed typos and Linux compiler target
|
| Revision 1.3  2001/06/15 08:42:42  prkipfer
| resolved precision and sign flaw
|
| Revision 1.2  2001/01/09 11:14:06  prkipfer
| precision epsilon is now static
|
| Revision 1.1  2001/01/03 15:30:12  prkipfer
| introduced intersection, linear system solver and minimizer
|
|
|
+---------------------------------------------------------------------*/

#ifdef LINUX
#include <cassert>
#endif

template <class T>
GoMinimizeND<T>::GoMinimizeND (int iDimensions, Function oF,
			       int uiMaxLevel, int uiMaxBracket,
			       int uiMaxIterations, void* pvUserData)
  : oneDMinimizer_(lineFunction,uiMaxLevel,uiMaxBracket)
{
  assert( iDimensions >= 1 && oF );

  dimensions_ = iDimensions;
  minFunction_ = oF;
  maxIterations_ = uiMaxIterations;
  userData_ = pvUserData;

  tCurr_ = new T[dimensions_];
  tSave_ = new T[dimensions_];
  directionStorage_ = new T[dimensions_*(dimensions_+1)];
  direction_ = new T*[dimensions_+1];
  for (int i = 0; i <= dimensions_; i++)
    direction_[i] = &directionStorage_[i*dimensions_];
  dConj_ = direction_[dimensions_];

  lineArg_ = new T[dimensions_];
}

template <class T>
GoMinimizeND<T>::~GoMinimizeND ()
{
  delete[] tCurr_;
  delete[] tSave_;
  delete[] directionStorage_;
  delete[] direction_;
  delete[] lineArg_;
}

template <class T>
void 
GoMinimizeND<T>::getMinimum (const T* afT0, const T* afT1,
			     const T* afTInitial, T* afTMin, T& rfFMin)
{
  // for 1D function callback
  oneDMinimizer_.userData() = this;

  // initial guess
  int uiQuantity = dimensions_*sizeof(T);
  fCurr_ = minFunction_(afTInitial,userData_);
  memcpy(tSave_,afTInitial,uiQuantity);
  memcpy(tCurr_,afTInitial,uiQuantity);

  // initialize direction set to the standard Euclidean basis
  memset(directionStorage_,0,uiQuantity*(dimensions_+1));
  int i, j;
  for (i = 0; i < dimensions_; i++)
    direction_[i][i] = T(1.0);

  T fL0, fL1, fLMin;
  for (int uiI = 0; uiI < maxIterations_; uiI++) {
    // find minimum in each direction and update current location
    for (i = 0; i < dimensions_; i++) {
      dCurr_ = direction_[i];
      computeDomain(afT0,afT1,fL0,fL1);
      oneDMinimizer_.getMinimum(fL0,fL1,T(0.0),fLMin,fCurr_);
      for (j = 0; j < dimensions_; j++)
	tCurr_[j] += fLMin*dCurr_[j];
    }

    // estimate a unit-length conjugate direction
    T fLength = T(0.0);
    for (i = 0; i < dimensions_; i++) {
      dConj_[i] = tCurr_[i] - tSave_[i];
      fLength += dConj_[i]*dConj_[i];
    }

    fLength = GbMath<T>::Sqrt(fLength);
    if ( fLength < GoMinimize1D<T>::EPSILON ) {
      // New position did not change significantly from old one.
      // Should there be a better convergence criterion here?
      break;
    }

    T fInvLength = T(1.0)/fLength;
    for (i = 0; i < dimensions_; i++)
      dConj_[i] *= fInvLength;

    // minimize in conjugate direction
    dCurr_ = dConj_;
    computeDomain(afT0,afT1,fL0,fL1);
    oneDMinimizer_.getMinimum(fL0,fL1,T(0.0),fLMin,fCurr_);
    for (j = 0; j < dimensions_; j++)
      tCurr_[j] += fLMin*dCurr_[j];

    // cycle the directions and add conjugate direction to set
    dConj_ = direction_[0];
    for (i = 0; i < dimensions_; i++)
      direction_[i] = direction_[i+1];

    // set parameters for next pass
    memcpy(tSave_,tCurr_,uiQuantity);
  }

  memcpy(afTMin,tCurr_,uiQuantity);
  rfFMin = fCurr_;
}

template <class T>
void
GoMinimizeND<T>::computeDomain (const T* afT0, const T* afT1,
				T& rfL0, T& rfL1)
{
  rfL0 = -std::numeric_limits<T>::max();
  rfL1 = +std::numeric_limits<T>::max();

  for (int i = 0; i < dimensions_; i++) {
    T fB0 = afT0[i] - tCurr_[i];
    T fB1 = afT1[i] - tCurr_[i];
    T fInv;
    if ( dCurr_[i] > T(0.0) ) {
      // valid t-interval is [B0,B1]
      fInv = T(1.0)/dCurr_[i];
      fB0 *= fInv;
      if ( fB0 > rfL0 )
	rfL0 = fB0;
      fB1 *= fInv;
      if ( fB1 < rfL1 )
	rfL1 = fB1;
    }
    else if ( dCurr_[i] < T(0.0) ) {
      // valid t-interval is [B1,B0]
      fInv = T(1.0)/dCurr_[i];
      fB0 *= fInv;
      if ( fB0 < rfL1 )
	rfL1 = fB0;
      fB1 *= fInv;
      if ( fB1 > rfL0 )
	rfL0 = fB1;
    }
  }

  // correction if numerical errors lead to values nearly zero
  if ( rfL0 > T(0.0) ) rfL0 = T(0.0);
  if ( rfL1 < T(0.0) ) rfL1 = T(0.0);
}

template <class T>
T 
GoMinimizeND<T>::lineFunction (T fT, void* pvUserData)
{
  GoMinimizeND* pkThis = (GoMinimizeND*) pvUserData;

  for (int i = 0; i < pkThis->dimensions_; i++) {
    pkThis->lineArg_[i] = pkThis->tCurr_[i] + fT*pkThis->dCurr_[i];
  }

  T fResult = pkThis->minFunction_(pkThis->lineArg_,pkThis->userData_);
  return fResult;
}
