/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOPMFILE_HH
#define GOPMFILE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iomanip>

#include "GoFile.hh"
#include "GoFileRegistry.hh"
#include "GbEdgeCollapse.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GoPMFileRegistration;

class GoPMFile
  : public GoFile
{
public:
  //! Constructor
  GoPMFile(const GbString&);
  virtual ~GoPMFile();

  virtual void read(std::vector<GoMesh*>& mv) const;
  virtual void readHEC(std::vector<GbHalfEdgeCollapse *> &hecList, std::vector<int> &levelNum) const;

  virtual void write(const GoMesh& m, GbBool smooth, int level=0) const;

  static GoFile* create(const GbString& filename);
  static const GbString className_;
  
  //! This operator pretty prints info about the file
//  friend std::ostream& operator<<(std::ostream&, const GoPMFile&);

private:
  //! hide copy constructor
  GoPMFile(const GoPMFile&);

  static const GoPMFileRegistration registration_;

};

extern GoFileRegistry FILE_IO_REGISTRY;

class GoPMFileRegistration
{
public:
  GoPMFileRegistration() {
//    debugmsg("registering "<<GoPMFile::className_);
    FILE_IO_REGISTRY.registerCreator(GoPMFile::className_,
				     GoPMFile::create);
    FILE_IO_REGISTRY.registerExtension(GbString("pm"),
				       GoPMFile::className_);
    FILE_IO_REGISTRY.registerDescription(GoPMFile::className_,
					 GbString("Progressive Mesh (*.pm)"));
  }
};


//#ifndef OUTLINE
//#include "GoPMFile.in"
//#endif

#endif // GOPMFILE_HH
/*----------------------------------------------------------------------
|
|
|
+---------------------------------------------------------------------*/
