// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.2  2001/06/15 09:02:30  prkipfer
| added compiler hints and removed useless signedness
|
| Revision 1.1  2001/02/13 13:29:30  prkipfer
| introduced basic rendering tools
|
|
|
+---------------------------------------------------------------------*/

template <class T>
GrDirectionalLight<T>::GrDirectionalLight()
  : direction_(GbVec3<T>::UNIT_Z)
{
  debugmsg("created directional light");
}

template <class T>
void 
GrDirectionalLight<T>::computeDiffuse (const GbMatrix3<T>& rkWorldRotate,
				       const GbVec3<T>& rkWorldTranslate, T fWorldScale,
				       const GbVec3<T>* akVertex, const GbVec3<T>* akNormal,
				       int uiQuantity, const GbBool* abVisible, GrColor<T>* akDiffuse)
{
  // transform light direction to model space of old mesh
  GbVec3<T> kModelDir = direction_*rkWorldRotate;

  // adjust diffuse color by light intensity
  GrColor<T> kAdjDiffuse = intensity_*diffuse_;

  for (int uiV = 0; uiV < uiQuantity; uiV++) {
    if ( abVisible[uiV] ) {
      T fDot = kModelDir.dot(akNormal[uiV]);
      if ( fDot < 0.0 )
	akDiffuse[uiV] -= fDot*kAdjDiffuse;
    }
  }
}

template <class T>
void 
GrDirectionalLight<T>::computeSpecular (const GbMatrix3<T>& rkWorldRotate,
					const GbVec3<T>& rkWorldTranslate, T fWorldScale,
					const GbVec3<T>* akVertex, const GbVec3<T>* akNormal,
					int uiQuantity, const GbBool* abVisible,
					const GbVec3<T>& rkCameraModelLocation, GrColor<T>* akSpecular)
{
  // transform light direction to model space of old mesh
  GbVec3<T> kModelDir = direction_*rkWorldRotate;

  // adjust specular color by light intensity
  GrColor<T> kAdjSpecular = intensity_*specular_;

  for (int uiV = 0; uiV < uiQuantity; uiV++) {
    if ( abVisible[uiV] ) {
      T fDot = kModelDir.dot(akNormal[uiV]);
      GbVec3<T> kReflect = kModelDir - (T(2.0)*fDot)*akNormal[uiV];
      GbVec3<T> kViewDir = rkCameraModelLocation - akVertex[uiV];
      fDot = kViewDir.dot(kReflect);
      if ( fDot > 0.0 ) {
	T fCoeff = fDot*fDot/kViewDir.dot(kViewDir);
	akSpecular[uiV] += fCoeff*kAdjSpecular;
      }
    }
  }
}

