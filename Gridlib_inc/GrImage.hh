/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GRIMAGE_HH
#define GRIMAGE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <map>

#include "GbTypes.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GrImage 
{
public:
  enum Type {
    IT_RGBA4444,
    IT_RGB888,
    IT_RGBA5551,
    IT_RGBA8888,
    IT_QUANTITY
  };

  // Construction and destruction.  GrImage accepts responsibility for
  // deleting the input array. The acImageName field is used as a
  // unique identifier for the image for purposes of sharing.  The caller of
  // the constructor may provided a name.  If not, the constructor generates
  // a unique name "imageN.gim" where N is the ID field. A
  // global map of images is maintained for sharing purposes.
  GrImage (Type eType, int uiWidth, int uiHeight, unsigned char* aucData, const char* acImageName = NULL);

  virtual ~GrImage ();

  // Streaming support.  The sharing system is automatically invoked by
  // these calls.  In Load, if an image corresponding to the filename is
  // already in memory, then that image is returned (i.e. shared).
  // Otherwise, a new image is created and returned.  The filename is used
  // as the image name.
  static GrImage* load (const char* acFilename);

  // This is intended to support saving procedurally generated images or
  // for utilities that convert to gridlib format from another format.  The filename
  // in this call does not replace the image name that might already exist.
  GbBool save (const char* acFilename);

  // Support for sharing images.  Note that setName is not
  // virtual, so subclasses hide it.  Do not circumvent this
  // function by explicitly calling GrImage::setName on a derived
  // object.
  void setName (const char* acName);

  INLINE Type getType () const;
  INLINE int getBytesPerPixel () const;

  INLINE int getWidth () const;
  INLINE int getHeight () const;
  INLINE int getNumPixels () const;

  INLINE unsigned char* getData () const;
  INLINE unsigned char* operator() (int uiIndex);

  class _GrImageInitTerm {
  public:
    _GrImageInitTerm () {
      GrImage::sharedImages_ = new std::map<GbKeyName,GrImage*>;
    }

    ~_GrImageInitTerm () {
      GrImage::sharedImages_->clear();
      delete GrImage::sharedImages_;
      GrImage::sharedImages_ = NULL;
    }
  };

protected:

  GrImage ();
  void setDefault ();

  Type type_;
  int width_, height_, wXh_;
  unsigned char* data_;

  static int bytesPerPixel_[IT_QUANTITY];

  // support for sharing images
//  friend class _GrImageInitTerm;
//  static const unsigned int mapSize_;
  char *name_;
  static std::map<GbKeyName,GrImage*>* sharedImages_;
  static unsigned int uiID_;

  static GbBool setAt (const char* acName, GrImage* pkImage);
  static GbBool getAt (const char* acName, GrImage*& rpkImage);
  static GbBool removeAt (const char* acName);
};

#ifndef OUTLINE
#include "GrImage.in"
#endif

#endif // GRIMAGE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.3  2001/06/19 16:30:21  prkipfer
| fixed typos and Linux compiler target
|
| Revision 1.2  2001/06/15 09:03:32  prkipfer
| switched to STL containers and introduced image sharing
|
| Revision 1.1  2001/02/13 13:29:32  prkipfer
| introduced basic rendering tools
|
|
+---------------------------------------------------------------------*/
