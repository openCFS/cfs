// -*- C++ -*- automatisch in C++mode wechseln (emacs)
/*----------------------------------------------------------------------
|
|
| $Log$
| Revision 1.3  2002/04/03 11:17:07  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.2  2001/06/15 09:27:31  prkipfer
| removed useless signedness
|
| Revision 1.1  2001/03/20 09:50:18  prkipfer
| introduced image handling tool
|
|
|
+---------------------------------------------------------------------*/

#include "GbImageConvert.hh"

template <class T>
GbTImage<T>::GbTImage (int uiDimensions, int* auiBound, T* atData)
  : GbLattice(uiDimensions,auiBound)
{
  setData(atData);
}

template <class T>
GbTImage<T>::GbTImage (const GbTImage& rkImage)
  : GbLattice(rkImage)
{
  data_ = new T[iSize_];
  memcpy(data_,rkImage.data_,iSize_*sizeof(T));
}

template <class T>
GbTImage<T>::GbTImage (const char* acFilename)
{
  load(acFilename);
}

template <class T>
GbTImage<T>::GbTImage (int uiDimensions)
  : GbLattice(uiDimensions)
{
}

template <class T>
GbTImage<T>::~GbTImage ()
{
  delete[] data_;
}

template <class T>
GbBool 
GbTImage<T>::load (const char* acFilename)
{
  std::ifstream kIStr(acFilename,std::ios::in|std::ios::binary);
  if ( !kIStr ) {
    dimensions_ = 0;
    iSize_ = 0;
    iBound_ = NULL;
    offset_ = NULL;
    data_ = NULL;
    return false;
  }

  if ( !GbLattice::load(kIStr) ) {
    data_ = NULL;
    return false;
  }

  int uiRTTI, uiSizeOf;
  kIStr.read((char*)&uiRTTI,sizeof(int));
  kIStr.read((char*)&uiSizeOf,sizeof(int));

  data_ = new T[iSize_];
  if ( uiRTTI == T::getRTTI() ) {
    kIStr.read((char*)data_,iSize_*sizeof(T));
  }
  else {
    char* acFileData = new char[iSize_*uiSizeOf];
    kIStr.read(acFileData,iSize_*uiSizeOf);
    GbImageConvert(iSize_,uiRTTI,acFileData,T::getRTTI(),data_);
    delete[] acFileData;
  }

  return true;
}

template <class T>
GbBool 
GbTImage<T>::save (const char* acFilename) const
{
  std::ofstream kOStr(acFilename,std::ios::out|std::ios::binary);
  if ( !kOStr )
    return false;

  if ( !GbLattice::save(kOStr) )
    return false;

  int uiRTTI = T::getRTTI(), uiSizeOf = sizeof(T);
  kOStr.write((const char*)&uiRTTI,sizeof(int));
  kOStr.write((const char*)&uiSizeOf,sizeof(int));
  kOStr.write((const char*)data_,iSize_*sizeof(T));

  return true;
}
