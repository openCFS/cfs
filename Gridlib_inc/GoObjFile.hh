/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOOBJFILE_HH
#define GOOBJFILE_HH

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

class GoObjFileRegistration;

class GoObjFile
  : public GoFile
{
public:
  //! Constructor
  GoObjFile(const GbString&);
  virtual ~GoObjFile();

  virtual void write(const GoMesh& m, GbBool smooth, int level=0) const;
  virtual void read(std::vector<GoMesh*>& mv) const;

  static GoFile* create(const GbString& filename);
  static const GbString className_;
  
  //! This operator pretty prints info about the file
//  friend std::ostream& operator<<(std::ostream&, const GoObjFile&);

private:
  //! hide copy constructor
  GoObjFile(const GoObjFile&);

  struct AddFace {
    AddFace(std::ofstream& f) : fh(f) {}
    void operator()(GoGeometryElement<float> *t) {
      fh<<"f ";
      for (int i=0; i<t->getNumVertices();i++) 
	fh<<t->getVertex(i)->getId()+1<<" ";
      fh<<std::endl;
    }
    std::ofstream& fh;
  };

  struct AddCoordinate {
    AddCoordinate(std::ofstream& f) : fh(f) {}
    void operator()(GoVertex<float> *t) {
      float x,y,z;
      t->getPosition(x,y,z);
      fh<<"v "<<std::setprecision(8)<<x<<" "<<std::setprecision(8)<<y<<" "<<std::setprecision(8)<<z<<std::endl;
    }
    std::ofstream& fh;
  };

  struct AddVertexNormal {
    AddVertexNormal(std::ofstream& f) : fh(f) {}
    void operator()(GoVertex<float> *t) {
      float x,y,z;
      t->getNormal(x,y,z);
      fh<<"vn "<<std::setprecision(8)<<x<<" "<<std::setprecision(8)<<y<<" "<<std::setprecision(8)<<z<<std::endl;
    }
    std::ofstream& fh;
  };

  static const GoObjFileRegistration registration_;
};

extern GoFileRegistry FILE_IO_REGISTRY;

class GoObjFileRegistration
{
public:
  GoObjFileRegistration() {
//    debugmsg("registering "<<GoObjFile::className_);
    FILE_IO_REGISTRY.registerCreator(GoObjFile::className_,
				     GoObjFile::create);
    FILE_IO_REGISTRY.registerExtension(GbString("obj"),
				       GoObjFile::className_);
    FILE_IO_REGISTRY.registerDescription(GoObjFile::className_,
					 GbString("Wavefront (*.obj)"));
  }
};


//#ifndef OUTLINE
//#include "GoObjFile.in"
//#endif

#endif // GOOBJFILE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.6  2001/02/13 11:28:41  prkipfer
| adapted to new GoFile::read method
|
| Revision 1.5  2000/11/28 16:39:37  uflabsik
| bug in obj format loader fixed
|
| Revision 1.4  2000/11/07 17:27:49  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.3  2000/10/23 09:33:06  prkipfer
| added possibility to save subdivision levels
|
| Revision 1.2  2000/07/17 11:21:37  prkipfer
| corrected floating point output precision
|
| Revision 1.1  2000/07/03 16:15:07  prkipfer
| extended IO subsystem, created class registry, added tetrahedron
|
|
+---------------------------------------------------------------------*/
