/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOUNVFILE_HH
#define GOUNVFILE_HH

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

class GoUnvFileRegistration;

class GoUnvFile
  : public GoFile
{
public:
  //! Constructor
  GoUnvFile(const GbString&);
  virtual ~GoUnvFile();

  virtual void read(std::vector<GoMesh*>& mv) const;

  static GoFile* create(const GbString& filename);
  static const GbString className_;
  
  //! This operator pretty prints info about the file
//  friend std::ostream& operator<<(std::ostream&, const GoUnvFile&);

private:
  //! hide copy constructor
  GoUnvFile(const GoUnvFile&);

  static const GoUnvFileRegistration registration_;
};

extern GoFileRegistry FILE_IO_REGISTRY;

class GoUnvFileRegistration
{
public:
  GoUnvFileRegistration() {
//    debugmsg("registering "<<GoUnvFile::className_);
    FILE_IO_REGISTRY.registerCreator(GoUnvFile::className_,
				     GoUnvFile::create);
    FILE_IO_REGISTRY.registerExtension(GbString("unv"),
				       GoUnvFile::className_);
    FILE_IO_REGISTRY.registerDescription(GoUnvFile::className_,
					 GbString("Universal File Format (*.unv)"));
  }
};


//#ifndef OUTLINE
//#include "GoUnvFile.in"
//#endif

#endif // GOUNVFILE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.2  2001/02/13 11:23:21  prkipfer
| adapted to new GoFile::read method
|
| Revision 1.1  2000/12/06 14:04:57  prkipfer
| added UNV File Type
|
|
+---------------------------------------------------------------------*/
