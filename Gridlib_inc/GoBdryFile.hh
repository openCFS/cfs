/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOBDRYFILE_HH
#define GOBDRYFILE_HH

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

class GoBdryFileRegistration;

class GoBdryFile
  : public GoFile
{
public:
  //! Constructor
  GoBdryFile(const GbString&);
  virtual ~GoBdryFile();

//    virtual void write(const GoMesh& m, GbBool smooth, int level=0) const;
  virtual void read(GoMesh& m) const;

  static GoFile* create(const GbString& filename);
  static const GbString className_;
  
  //! This operator pretty prints info about the file
//  friend std::ostream& operator<<(std::ostream&, const GoOffFile&);

private:
  //! hide copy constructor
  GoBdryFile(const GoBdryFile&);

  static const GoBdryFileRegistration registration_;
};

extern GoFileRegistry FILE_IO_REGISTRY;

class GoBdryFileRegistration
{
public:
  GoBdryFileRegistration() {
//    debugmsg("registering "<<GoBdryFile::className_);
    FILE_IO_REGISTRY.registerCreator(GoBdryFile::className_,
				     GoBdryFile::create);
    FILE_IO_REGISTRY.registerExtension(GbString("bdry"),
				       GoBdryFile::className_);
    FILE_IO_REGISTRY.registerDescription(GoBdryFile::className_,
					 GbString("Boundary Geometry File (*.bdry)"));
  }
};


//#ifndef OUTLINE
//#include "GoBdryFile.in"
//#endif

#endif // GOBDRYFILE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.1  2000/10/23 09:19:12  prkipfer
| initial version
|
|
+---------------------------------------------------------------------*/
