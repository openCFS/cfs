/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOUTMFILE_HH
#define GOUTMFILE_HH

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

class GoUTMFileRegistration;

class GoUTMFile
  : public GoFile
{
public:
  //! Constructor
  GoUTMFile(const GbString&);
  virtual ~GoUTMFile();

  virtual void read(std::vector<GoMesh*>& mv) const;
  virtual void write(const GoMesh& m, GbBool smooth, int level=0) const;

  static GoFile* create(const GbString& filename);
  static const GbString className_;
  
  //! This operator pretty prints info about the file
//  friend std::ostream& operator<<(std::ostream&, const GoUTMFile&);

private:
  //! hide copy constructor
  GoUTMFile(const GoUTMFile&);

  struct AddTet {
    AddTet(std::ofstream& f) : fh(f) {}
    void operator()(GoGeometryElement<float> *t) {
      for (int i=0; i<t->getNumVertices();i++) 
	fh<<t->getVertex(i)->getId()<<" ";
      fh<<std::endl;
    }
    std::ofstream& fh;
  };

  struct UnmarkAll {
    UnmarkAll() {}
    void operator()(GoGeometryElement<float> *t) {
      t->delFlag(MARK_G);
    }
  };

  struct TestAddTri {
    TestAddTri(GbVertexFlag v, std::ofstream& f) : flag(v), fh(f) {}
    void operator()(GoGeometryElement<float> *t) {
      if (! t->testFlag(MARK_G)) {
	int n=0;
	// count number of surface triangles in tetrahedron
	for (int i=0; i<t->getNumVertices();i++) {
	  GoVertex<float> *v = t->getVertex(i);
	  if (v->testFlag(BOUNDARY_V) && v->testFlag(flag)) {
	    t->setFlag(MARK_G);
	    ++n;
	  }
	}
	if (n==3) {
	  for (int i=0; i<t->getNumVertices();i++) {
	    GoVertex<float> *v = t->getVertex(i);
	    if (v->testFlag(BOUNDARY_V) && v->testFlag(flag)) {
	      fh<<v->getId()<<" ";
	    }
	  }
	  fh<<std::endl;
	}
	else if (n==4) {
	  // obviously they are all boundary vertices
	  for (int i=0; i<3;i++)
	    fh<<t->getVertex(i)->getId()<<" ";
	  fh<<std::endl;
	  for (int i=1; i<4;i++)
	    fh<<t->getVertex(i)->getId()<<" ";
	  fh<<std::endl;
	}
      }
    }
    GbVertexFlag flag;
    std::ofstream& fh;
  };

  struct AddCoordinate {
    AddCoordinate(std::ofstream& f) : fh(f) {}
    void operator()(GoVertex<float> *t) {
      float x,y,z;
      t->getPosition(x,y,z);
      fh<<std::fixed<<std::setprecision(8)<<x<<" "<<std::fixed<<std::setprecision(8)<<y<<" "<<std::fixed<<std::setprecision(8)<<z<<std::endl;
    }
    std::ofstream& fh;
  };

  struct CheckFlag {
    CheckFlag(GbVertexFlag f, unsigned int* c) : flag(f), cnt(c) {}
    void operator()(GoVertex<float> *t) {
      if (t->testFlag(BOUNDARY_V) && t->testFlag(flag)) ++(*cnt);
    }
    GbVertexFlag flag;
    unsigned int *cnt;
  };

  static const GoUTMFileRegistration registration_;
};

extern GoFileRegistry FILE_IO_REGISTRY;

class GoUTMFileRegistration
{
public:
  GoUTMFileRegistration() {
//    debugmsg("registering "<<GoUTMFile::className_);
    FILE_IO_REGISTRY.registerCreator(GoUTMFile::className_,
				     GoUTMFile::create);
    FILE_IO_REGISTRY.registerExtension(GbString("utmASCII"),
				       GoUTMFile::className_);
    FILE_IO_REGISTRY.registerDescription(GoUTMFile::className_,
					 GbString("Udo Tremel Fileformat (*.utmASCII)"));
  }
};


//#ifndef OUTLINE
//#include "GoUTMFile.in"
//#endif

#endif // GOUTMFILE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.3  2001/04/17 12:39:43  prkipfer
| added write method
|
| Revision 1.2  2001/02/13 11:26:46  prkipfer
| adapted to new GoFile::read method
|
| Revision 1.1  2000/12/12 10:04:43  prkipfer
| introduced UTM file format
|
|
+---------------------------------------------------------------------*/
