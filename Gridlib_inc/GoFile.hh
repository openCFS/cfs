/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOFILE_HH
#define GOFILE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iostream>
#include <fstream>
#include <vector>

#include "GbTypes.hh"
#include "GoVertex.hh"
#include "GoGeometryElement.hh"
#include "GoMesh.hh"
#include "GbMeshFunctions.hh"
#include "GbEdgeCollapse.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GoFile
{
public:
  //! Constructor
  virtual ~GoFile();

  static GoFile* guess(const GbString& s);
  static GbString identify(const GbString& s);
  static GbString listTypes();

  virtual GbBool has(GbContents f);
  virtual int content();

  virtual void write(const GoMesh& m, GbBool smooth, int level=0) const;
  virtual void write(std::vector<GoMesh*>& mv) const;
  virtual void read(std::vector<GoMesh*>& mv) const;

  virtual void readHEC(std::vector<GbHalfEdgeCollapse *> &hecList) const;

  //! there is no create routine - a instance of this class is nonsense anyway

  //! This operator pretty prints info about the file
//  friend std::ostream& operator<<(std::ostream&, const GoFile&);

protected:
  GbString filename_;
  mutable int contents_;

  //! protect constructor
  GoFile(const GbString&);

  //! hide copy constructor
  GoFile(const GoFile&);
};


//#ifndef OUTLINE
//#include "GoFile.in"
//#endif

#endif // GOFILE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.7  2001/08/20 09:15:51  prkipfer
| made content_ mutable for parsers that dont know the file content in advance
|
| Revision 1.6  2001/05/03 10:21:39  uflabsik
| hec reading added
|
| Revision 1.5  2001/02/13 11:35:17  prkipfer
| introduced methods to read and write multiple meshes from one file
|
| Revision 1.4  2000/11/07 17:27:47  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.3  2000/10/23 09:33:05  prkipfer
| added possibility to save subdivision levels
|
| Revision 1.2  2000/07/25 09:10:17  prkipfer
| new content enum
|
| Revision 1.1  2000/07/03 16:15:05  prkipfer
| extended IO subsystem, created class registry, added tetrahedron
|
|
+---------------------------------------------------------------------*/
