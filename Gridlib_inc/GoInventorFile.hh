/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOINVENTORFILE_HH
#define GOINVENTORFILE_HH

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

class GoInventorFileRegistration;

class GoInventorFile
  : public GoFile
{
public:
  //! Constructor
  GoInventorFile(const GbString&);
  virtual ~GoInventorFile();

  virtual void write(const GoMesh& m, GbBool smooth, int level=0) const;

  static GoFile* create(const GbString& filename);
  static const GbString className_;
  
  //! This operator pretty prints info about the file
//  friend std::ostream& operator<<(std::ostream&, const GoInventorFile&);

private:
  //! hide copy constructor
  GoInventorFile(const GoInventorFile&);

  struct AddFace {
    AddFace(std::ofstream& f) : fh(f) {}
    void operator()(GoGeometryElement<float> *t) {
      fh<<"\t";
      for (int i=0; i<t->getNumVertices();i++) 
	fh<<t->getVertex(i)->getId()<<", ";
      fh<<"-1,"<<std::endl;
    }
    std::ofstream& fh;
  };

  struct AddFaceNormal {
    AddFaceNormal(std::ofstream& f) : fh(f) {}
    void operator()(GoGeometryElement<float> *t) {
      GbVec3<float> n = t->getNormal();
      fh<<"\t"<<std::setprecision(8)<<n[0]<<" "<<std::setprecision(8)<<n[1]<<" "<<std::setprecision(8)<<n[2]<<","<<std::endl;
    }
    std::ofstream& fh;
  };

  struct AddCoordinate {
    AddCoordinate(std::ofstream& f) : fh(f) {}
    void operator()(GoVertex<float> *t) {
      float x,y,z;
      t->getPosition(x,y,z);
      fh<<"\t"<<std::setprecision(8)<<x<<" "<<std::setprecision(8)<<y<<" "<<std::setprecision(8)<<z<<","<<std::endl;
    }
    std::ofstream& fh;
  };

  struct AddVertexNormal {
    AddVertexNormal(std::ofstream& f) : fh(f) {}
    void operator()(GoVertex<float> *t) {
      float x,y,z;
      t->getNormal(x,y,z);
      fh<<"\t"<<std::setprecision(8)<<x<<" "<<std::setprecision(8)<<y<<" "<<std::setprecision(8)<<z<<","<<std::endl;
    }
    std::ofstream& fh;
  };

  static const GoInventorFileRegistration registration_;
};

extern GoFileRegistry FILE_IO_REGISTRY;

class GoInventorFileRegistration
{
public:
  GoInventorFileRegistration() {
//    debugmsg("registering "<<GoInventorFile::className_);
    FILE_IO_REGISTRY.registerCreator(GoInventorFile::className_,
				     GoInventorFile::create);
    FILE_IO_REGISTRY.registerExtension(GbString("iv"),
				       GoInventorFile::className_);
    FILE_IO_REGISTRY.registerDescription(GoInventorFile::className_,
					 GbString("OpenInventor (*.iv)"));
  }
};


//#ifndef OUTLINE
//#include "GoInventorFile.in"
//#endif

#endif // GOINVENTORFILE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.5  2000/11/07 17:27:49  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.4  2000/10/23 09:33:06  prkipfer
| added possibility to save subdivision levels
|
| Revision 1.3  2000/07/17 11:20:10  prkipfer
| corrected floating point output precision
|
| Revision 1.2  2000/07/03 16:15:06  prkipfer
| extended IO subsystem, created class registry, added tetrahedron
|
| Revision 1.1  2000/06/16 13:33:06  prkipfer
| just starting this subsystem...
|
|
+---------------------------------------------------------------------*/
