// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.2  2001/06/15 09:04:44  prkipfer
| added compiler hints and removed useless signedness
|
| Revision 1.1  2001/02/13 13:29:35  prkipfer
| introduced basic rendering tools
|
|
|
+---------------------------------------------------------------------*/

template <class T>
GrPointLight<T>::GrPointLight ()
  : location_(GbVec3<T>::ZERO)
{
  debugmsg("created point light");
}

template <class T>
void 
GrPointLight<T>::computeDiffuse (const GbMatrix3<T>& rkWorldRotate,
				 const GbVec3<T>& rkWorldTranslate, T fWorldScale,
				 const GbVec3<T>* akVertex, const GbVec3<T>* akNormal,
				 int uiQuantity, const GbBool* abVisible, GrColor<T>* akDiffuse)
{
  // transform light position to model space of old mesh
  T fInvWScale = 1.0/fWorldScale;
  GbVec3<T> kDiff = location_ - rkWorldTranslate;
  GbVec3<T> kModelPos = (kDiff*rkWorldRotate)*fInvWScale;

  // adjust diffuse color by light intensity
  GrColor<T> kAdjDiffuse = intensity_*diffuse_;

  for (int uiV = 0; uiV < uiQuantity; uiV++) {
    if ( abVisible[uiV] ) {
      kDiff = kModelPos - akVertex[uiV];
      T fDot = kDiff.dot(akNormal[uiV]);
      if ( fDot > 0.0 ) {
	fDot /= kDiff.getNorm();
	akDiffuse[uiV] += fDot*kAdjDiffuse;
      }
    }
  }
}

template <class T>
void 
GrPointLight<T>::computeSpecular (const GbMatrix3<T>& rkWorldRotate,
				  const GbVec3<T>& rkWorldTranslate, T fWorldScale,
				  const GbVec3<T>* akVertex, const GbVec3<T>* akNormal,
				  int uiQuantity, const GbBool* abVisible,
				  const GbVec3<T>& rkCameraModelLocation, GrColor<T>* akSpecular)
{
  // transform light position to model space of old mesh
  T fInvWScale = 1.0/fWorldScale;
  GbVec3<T> kDiff = location_ - rkWorldTranslate;
  GbVec3<T> kModelPos = (kDiff*rkWorldRotate)*fInvWScale;

  // adjust specular color by light intensity
  GrColor<T> kAdjSpecular = intensity_*specular_;

  for (int uiV = 0; uiV < uiQuantity; uiV++) {
    if ( abVisible[uiV] ) {
      kDiff = kModelPos - akVertex[uiV];
      T fDot = kDiff.dot(akNormal[uiV]);
      GbVec3<T> kReflect = (T(2.0)*fDot)*akNormal[uiV] - kDiff;
      GbVec3<T> kViewDir = rkCameraModelLocation - akVertex[uiV];
      fDot = kViewDir.dot(kReflect);
      if ( fDot > 0.0 ) {
	T fVSL = kViewDir.dot(kViewDir);
	T fRSL = kReflect.dot(kReflect);
	T fDenom = GbMath<T>::Sqrt(fVSL*fRSL);
	akSpecular[uiV] += (fDot/fDenom)*kAdjSpecular;
      }
    }
  }
}

