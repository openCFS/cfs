// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.2  2002/03/21 14:58:58  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.2  2001/06/15 09:04:44  prkipfer
| added compiler hints and removed useless signedness
|
| Revision 1.1  2001/02/13 13:29:36  prkipfer
| introduced basic rendering tools
|
|
|
+---------------------------------------------------------------------*/

template <class T>
GrSpotLight<T>::GrSpotLight()
  : direction_(-GbVec3<T>::UNIT_Z)
{
  setAngle(GbMath<T>::PI);
  exponent_ = 0.0;
}

template <class T>
void 
GrSpotLight<T>::setAngle (T fAngle)
{
  angle_ = fAngle;
  cosSqr_ = GbMath<T>::Cos(angle_);
  cosSqr_ *= cosSqr_;
  sinSqr_ = GbMath<T>::Sin(angle_);
  sinSqr_ *= sinSqr_;
}

template <class T>
void 
GrSpotLight<T>::computeDiffuse (const GbMatrix3<T>& rkWorldRotate,
				const GbVec3<T>& rkWorldTranslate, T fWorldScale,
				const GbVec3<T>* akVertex, const GbVec3<T>* akNormal,
				int uiQuantity, const GbBool* abVisible, GrColor<T>* akDiffuse)
{
  // transform light position to model space of old mesh
  T fInvWScale = 1.0/fWorldScale;
  GbVec3<T> kDiff = location_ - rkWorldTranslate;
  GbVec3<T> kModelPos = (kDiff*rkWorldRotate)*fInvWScale;
  GbVec3<T> kModelDir = direction_*rkWorldRotate;

  // adjust diffuse color by light intensity
  GrColor<T> kAdjDiffuse = intensity_*diffuse_;

  for (int uiV = 0; uiV < uiQuantity; uiV++) {
    if ( abVisible[uiV] ) {
      kDiff = kModelPos - akVertex[uiV];
      T fDiffDotDir = kDiff.dot(kModelDir);
      if ( fDiffDotDir >= 0.0 )
	continue;

      T fLenSqr = kDiff.dot(kDiff);
      T fDiffDotDirSqr = fDiffDotDir*fDiffDotDir;
      T fAngleAttenuate = fDiffDotDirSqr - fLenSqr*cosSqr_;
      if ( fAngleAttenuate <= 0.0 )
	continue;

      T fDot = kDiff.dot(akNormal[uiV]);
      if ( fDot > 0.0 ) {
	T fNumer = fDot*fAngleAttenuate;
	T fDenom = fDiffDotDirSqr*sinSqr_*GbMath<T>::Sqrt(fLenSqr);
	akDiffuse[uiV] += (fNumer/fDenom)*kAdjDiffuse;
      }
    }
  }
}

template <class T>
void 
GrSpotLight<T>::computeSpecular (const GbMatrix3<T>& rkWorldRotate,
				 const GbVec3<T>& rkWorldTranslate, T fWorldScale,
				 const GbVec3<T>* akVertex, const GbVec3<T>* akNormal,
				 int uiQuantity, const GbBool* abVisible,
				 const GbVec3<T>& rkCameraModelLocation, GrColor<T>* akSpecular)
{
  // transform light position to model space of old mesh
  T fInvWScale = 1.0/fWorldScale;
  GbVec3<T> kDiff = location_ - rkWorldTranslate;
  GbVec3<T> kModelPos = (kDiff*rkWorldRotate)*fInvWScale;
  GbVec3<T> kModelDir = direction_*rkWorldRotate;

  // adjust diffuse color by light intensity
  GrColor<T> kAdjSpecular = intensity_*specular_;

  for (int uiV = 0; uiV < uiQuantity; uiV++) {
    if ( abVisible[uiV] ) {
      kDiff = kModelPos - akVertex[uiV];
      T fDiffDotDir = kDiff.dot(kModelDir);
      if ( fDiffDotDir >= 0.0 )
	continue;

      T fLenSqr = kDiff.dot(kDiff);
      T fDiffDotDirSqr = fDiffDotDir*fDiffDotDir;
      T fAngleAttenuate = fDiffDotDirSqr - fLenSqr*cosSqr_;
      if ( fAngleAttenuate <= 0.0 )
	continue;

      T fDot = kDiff.dot(akNormal[uiV]);
      GbVec3<T> kReflect = (T(2.0)*fDot)*akNormal[uiV] - kDiff;
      GbVec3<T> kViewDir = rkCameraModelLocation - akVertex[uiV];
      fDot = kViewDir.dot(kReflect);
      if ( fDot > 0.0 ) {
	T fNumer = fDot*fDot*fAngleAttenuate;
	T fDenom = fDiffDotDirSqr*sinSqr_ * GbMath<T>::Sqrt(fLenSqr)*kViewDir.dot(kViewDir);
	akSpecular[uiV] += (fNumer/fDenom)*kAdjSpecular;
      }
    }
  }
}

