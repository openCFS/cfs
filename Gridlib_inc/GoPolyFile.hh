/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOPOLYFILE_HH
#define GOPOLYFILE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iomanip>

#include "GoFile.hh"
#include "GoFileRegistry.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GoPolyFileRegistration;

class GoPolyFile
  : public GoFile
{
public:
  //! Constructor
  GoPolyFile(const GbString&);
  virtual ~GoPolyFile();

//    virtual void write(const GoMesh& m, GbBool smooth, int level=0) const;
  virtual void read(std::vector<GoMesh*>& mv) const;

  static GoFile* create(const GbString& filename);
  static const GbString className_;
  
  //! This operator pretty prints info about the file
//  friend std::ostream& operator<<(std::ostream&, const GoOffFile&);

private:
  //! hide copy constructor
  GoPolyFile(const GoPolyFile&);

  static const GoPolyFileRegistration registration_;
};

extern GoFileRegistry FILE_IO_REGISTRY;

class GoPolyFileRegistration
{
public:
  GoPolyFileRegistration() {
//    debugmsg("registering "<<GoPolyFile::className_);
    FILE_IO_REGISTRY.registerCreator(GoPolyFile::className_,
				     GoPolyFile::create);
    FILE_IO_REGISTRY.registerExtension(GbString("poly"),
				       GoPolyFile::className_);
    FILE_IO_REGISTRY.registerDescription(GoPolyFile::className_,
					 GbString("Polygonal Surface File (*.poly)"));
  }
};


//#ifndef OUTLINE
//#include "GoPolyFile.in"
//#endif

#endif // GOPOLYFILE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.2  2001/02/13 11:27:26  prkipfer
| adapted to new GoFile::read method
|
| Revision 1.1  2000/10/30 13:14:06  prkipfer
| initial revision
|
|
+---------------------------------------------------------------------*/
