// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/01/02 14:49:48  prkipfer
| introduced new classes
|
|
|
+---------------------------------------------------------------------*/

template <class T>
GbQuaternion<T>::GbQuaternion(T fW, T fX, T fY, T fZ)
{
  w_ = fW;
  x_ = fX;
  y_ = fY;
  z_ = fZ;
}

template <class T>
GbQuaternion<T>::GbQuaternion(const GbQuaternion<T>& rkQ)
{
  w_ = rkQ.w_;
  x_ = rkQ.x_;
  y_ = rkQ.y_;
  z_ = rkQ.z_;
}

template <class T>
void 
GbQuaternion<T>::fromRotationMatrix (const GbMatrix3<T>& kRot)
{
  // Algorithm in Ken Shoemake's article in 1987 SIGGRAPH course notes
  // article "Quaternion Calculus and Fast Animation".

  T fTrace = kRot[0][0]+kRot[1][1]+kRot[2][2];
  T fRoot;

  if ( fTrace > 0.0 ) {
    // |w| > 1/2, may as well choose w > 1/2
    fRoot = GbMath<T>::Sqrt(fTrace + 1.0);  // 2w
    w_ = 0.5*fRoot;
    fRoot = 0.5/fRoot;  // 1/(4w)
    x_ = (kRot[2][1]-kRot[1][2])*fRoot;
    y_ = (kRot[0][2]-kRot[2][0])*fRoot;
    z_ = (kRot[1][0]-kRot[0][1])*fRoot;
  }
  else {
    // |w| <= 1/2
    static int s_iNext[3] = { 1, 2, 0 };
    int i = 0;
    if ( kRot[1][1] > kRot[0][0] )
      i = 1;
    if ( kRot[2][2] > kRot[i][i] )
      i = 2;
    int j = s_iNext[i];
    int k = s_iNext[j];

    fRoot = GbMath<T>::Sqrt(kRot[i][i]-kRot[j][j]-kRot[k][k] + 1.0);
    T* apkQuat[3] = { &x_, &y_, &z_ };
    *apkQuat[i] = 0.5*fRoot;
    fRoot = 0.5/fRoot;
    w_ = (kRot[k][j]-kRot[j][k])*fRoot;
    *apkQuat[j] = (kRot[j][i]+kRot[i][j])*fRoot;
    *apkQuat[k] = (kRot[k][i]+kRot[i][k])*fRoot;
  }
}

template <class T>
void 
GbQuaternion<T>::toRotationMatrix (GbMatrix3<T>& kRot) const
{
  T fTx  = 2.0*x_;
  T fTy  = 2.0*y_;
  T fTz  = 2.0*z_;
  T fTwx = fTx*w_;
  T fTwy = fTy*w_;
  T fTwz = fTz*w_;
  T fTxx = fTx*x_;
  T fTxy = fTy*x_;
  T fTxz = fTz*x_;
  T fTyy = fTy*y_;
  T fTyz = fTz*y_;
  T fTzz = fTz*z_;
  
  kRot[0][0] = 1.0-(fTyy+fTzz);
  kRot[0][1] = fTxy-fTwz;
  kRot[0][2] = fTxz+fTwy;
  kRot[1][0] = fTxy+fTwz;
  kRot[1][1] = 1.0-(fTxx+fTzz);
  kRot[1][2] = fTyz-fTwx;
  kRot[2][0] = fTxz-fTwy;
  kRot[2][1] = fTyz+fTwx;
  kRot[2][2] = 1.0-(fTxx+fTyy);
}

template <class T>
void 
GbQuaternion<T>::fromAngleAxis (const T& rfAngle, const GbVec3<T>& rkAxis)
{
  // assert:  axis[] is unit length
  //
  // The quaternion representing the rotation is
  //   q = cos(A/2)+sin(A/2)*(x*i+y*j+z*k)

  T fHalfAngle = 0.5*rfAngle;
  T fSin = GbMath<T>::Sin(fHalfAngle);
  w_ = GbMath<T>::Cos(fHalfAngle);
  x_ = fSin*rkAxis[0];
  y_ = fSin*rkAxis[1];
  z_ = fSin*rkAxis[2];
}

template <class T>
void 
GbQuaternion<T>::toAngleAxis (T& rfAngle, GbVec3<T>& rkAxis) const
{
  // The quaternion representing the rotation is
  //   q = cos(A/2)+sin(A/2)*(x*i+y*j+z*k)

  T fSqrLength = x_*x_+y_*y_+z_*z_;
  if ( fSqrLength > 0.0 ) {
    rfAngle = 2.0*GbMath<T>::ACos(w_);
    T fInvLength = 1.0/GbMath<T>::Sqrt(fSqrLength);
    rkAxis[0] = x_*fInvLength;
    rkAxis[1] = y_*fInvLength;
    rkAxis[2] = z_*fInvLength;
  }
  else {
    // angle is 0 (mod 2*pi), so any axis will do
    rfAngle = 0.0;
    rkAxis[0] = 1.0;
    rkAxis[1] = 0.0;
    rkAxis[2] = 0.0;
  }
}

template <class T>
void 
GbQuaternion<T>::fromAxes (const GbVec3<T>* akAxis)
{
  GbMatrix3<T> kRot;

  for (int iCol = 0; iCol < 3; iCol++) {
    kRot[0][iCol] = akAxis[iCol][0];
    kRot[1][iCol] = akAxis[iCol][1];
    kRot[2][iCol] = akAxis[iCol][2];
  }

  fromRotationMatrix(kRot);
}

template <class T>
void 
GbQuaternion<T>::toAxes (GbVec3<T>* akAxis) const
{
  GbMatrix3<T> kRot;

  toRotationMatrix(kRot);

  for (int iCol = 0; iCol < 3; iCol++) {
    akAxis[iCol][0] = kRot[0][iCol];
    akAxis[iCol][1] = kRot[1][iCol];
    akAxis[iCol][2] = kRot[2][iCol];
  }
}

template <class T>
GbQuaternion<T> 
GbQuaternion<T>::exp () const
{
  // If q = A*(x*i+y*j+z*k) where (x,y,z) is unit length, then
  // exp(q) = cos(A)+sin(A)*(x*i+y*j+z*k).  If sin(A) is near zero,
  // use exp(q) = cos(A)+A*(x*i+y*j+z*k) since A/sin(A) has limit 1.

  T fAngle = GbMath<T>::Sqrt(x_*x_+y_*y_+z_*z_);
  T fSin = GbMath<T>::Sin(fAngle);

  GbQuaternion<T> kResult;
  kResult.w_ = GbMath<T>::Cos(fAngle);

  if ( GbMath<T>::Abs(fSin) >= GbQuaternion<T>::EPSILON ) {
    T fCoeff = fSin/fAngle;
    kResult.x_ = fCoeff*x_;
    kResult.y_ = fCoeff*y_;
    kResult.z_ = fCoeff*z_;
  }
  else {
    kResult.x_ = x_;
    kResult.y_ = y_;
    kResult.z_ = z_;
  }

  return kResult;
}

template <class T>
GbQuaternion<T> 
GbQuaternion<T>::log () const
{
  // If q = cos(A)+sin(A)*(x*i+y*j+z*k) where (x,y,z) is unit length, then
  // log(q) = A*(x*i+y*j+z*k).  If sin(A) is near zero, use log(q) =
  // sin(A)*(x*i+y*j+z*k) since sin(A)/A has limit 1.

  GbQuaternion<T> kResult;
  kResult.w_ = 0.0;

  if ( GbMath<T>::Abs(w_) < 1.0 ) {
    T fAngle = GbMath<T>::ACos(w_);
    T fSin = GbMath<T>::Sin(fAngle);
    if ( GbMath<T>::Abs(fSin) >= GbQuaternion<T>::EPSILON ) {
      T fCoeff = fAngle/fSin;
      kResult.x_ = fCoeff*x_;
      kResult.y_ = fCoeff*y_;
      kResult.z_ = fCoeff*z_;
      return kResult;
    }
  }

  kResult.x_ = x_;
  kResult.y_ = y_;
  kResult.z_ = z_;

  return kResult;
}

template <class T>
GbVec3<T> 
GbQuaternion<T>::operator* (const GbVec3<T>& rkVector) const
{
  // Given a vector u = (x0,y0,z0) and a unit length quaternion
  // q = <w,x,y,z>, the vector v = (x1,y1,z1) which represents the
  // rotation of u by q is v = q*u*q^{-1} where * indicates quaternion
  // multiplication and where u is treated as the quaternion <0,x0,y0,z0>.
  // Note that q^{-1} = <w,-x,-y,-z>, so no real work is required to
  // invert q.  Now
  //
  //   q*u*q^{-1} = q*<0,x0,y0,z0>*q^{-1}
  //     = q*(x0*i+y0*j+z0*k)*q^{-1}
  //     = x0*(q*i*q^{-1})+y0*(q*j*q^{-1})+z0*(q*k*q^{-1})
  //
  // As 3-vectors, q*i*q^{-1}, q*j*q^{-1}, and 2*k*q^{-1} are the columns
  // of the rotation matrix computed in GbQuaternion<T>::ToRotationMatrix.
  // The vector v is obtained as the product of that rotation matrix with
  // vector u.  As such, the quaternion representation of a rotation
  // matrix requires less space than the matrix and more time to compute
  // the rotated vector.  Typical space-time tradeoff...

  GbMatrix3<T> kRot;
  toRotationMatrix(kRot);
  return kRot*rkVector;
}

template <class T>
GbQuaternion<T> 
GbQuaternion<T>::slerp (T fT, const GbQuaternion<T>& rkP, const GbQuaternion<T>& rkQ)
{
  T fCos = rkP.dot(rkQ);
  T fAngle = GbMath<T>::ACos(fCos);

  if ( GbMath<T>::Abs(fAngle) < GbQuaternion<T>::EPSILON )
    return rkP;

  T fSin = GbMath<T>::Sin(fAngle);
  T fInvSin = 1.0/fSin;
  T fCoeff0 = GbMath<T>::Sin((1.0-fT)*fAngle)*fInvSin;
  T fCoeff1 = GbMath<T>::Sin(fT*fAngle)*fInvSin;
  return fCoeff0*rkP + fCoeff1*rkQ;
}

template <class T>
GbQuaternion<T> 
GbQuaternion<T>::slerpExtraSpins (T fT, const GbQuaternion<T>& rkP, const GbQuaternion<T>& rkQ, int iExtraSpins)
{
  T fCos = rkP.dot(rkQ);
  T fAngle = GbMath<T>::ACos(fCos);

  if ( GbMath<T>::Abs(fAngle) < GbQuaternion<T>::EPSILON )
    return rkP;

  T fSin = GbMath<T>::Sin(fAngle);
  T fPhase = GbMath<T>::PI*iExtraSpins*fT;
  T fInvSin = 1.0/fSin;
  T fCoeff0 = GbMath<T>::Sin((1.0-fT)*fAngle - fPhase)*fInvSin;
  T fCoeff1 = GbMath<T>::Sin(fT*fAngle + fPhase)*fInvSin;
  return fCoeff0*rkP + fCoeff1*rkQ;
}

template <class T>
void 
GbQuaternion<T>::intermediate (const GbQuaternion<T>& rkQ0, const GbQuaternion<T>& rkQ1, const GbQuaternion<T>& rkQ2,
			       GbQuaternion<T>& rkA, GbQuaternion<T>& rkB)
{
  // assert:  q0, q1, q2 are unit quaternions

  GbQuaternion<T> kQ0inv = rkQ0.unitInverse();
  GbQuaternion<T> kQ1inv = rkQ1.unitInverse();
  GbQuaternion<T> rkP0 = kQ0inv*rkQ1;
  GbQuaternion<T> rkP1 = kQ1inv*rkQ2;
  GbQuaternion<T> kArg = 0.25*(rkP0.log()-rkP1.log());
  GbQuaternion<T> kMinusArg = -kArg;

  rkA = rkQ1*kArg.exp();
  rkB = rkQ1*kMinusArg.exp();
}

template <class T>
GbQuaternion<T> 
GbQuaternion<T>::squad (T fT, const GbQuaternion<T>& rkP, const GbQuaternion<T>& rkA,
			const GbQuaternion<T>& rkB, const GbQuaternion<T>& rkQ)
{
  T fSlerpT = 2.0*fT*(1.0-fT);
  GbQuaternion<T> kSlerpP = slerp(fT,rkP,rkQ);
  GbQuaternion<T> kSlerpQ = slerp(fT,rkA,rkB);
  return slerp(fSlerpT,kSlerpP,kSlerpQ);
}
