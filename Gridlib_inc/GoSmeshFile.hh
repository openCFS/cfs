/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOSMESHFILE_HH
#define GOSMESHFILE_HH

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

class GoSmeshFileRegistration;

class GoSmeshFile
  : public GoFile
{
public:
  //! Constructor
  GoSmeshFile(const GbString&);
  virtual ~GoSmeshFile();

  virtual void read(std::vector<GoMesh*>& mv) const;
  virtual void write(const GoMesh& m, GbBool smooth, int level=0) const;

  static GoFile* create(const GbString& filename);
  static const GbString className_;
  
  //! This operator pretty prints info about the file
//  friend std::ostream& operator<<(std::ostream&, const GoOffFile&);

private:
  //! hide copy constructor
  GoSmeshFile(const GoSmeshFile&);

  struct AddFace {
    AddFace(std::ofstream& f) : fh(f) {}
    void operator()(GoGeometryElement<float> *t) {
      fh<<t->getNumVertices()<<" ";
      for (int i=0; i<t->getNumVertices();i++) 
	fh<<t->getVertex(i)->getId()<<" ";
      fh<<std::endl;
    }
    std::ofstream& fh;
  };

  struct AddCoordinate {
    AddCoordinate(std::ofstream& f) : fh(f) {}
    void operator()(GoVertex<float> *t) {
      float x,y,z;
      t->getPosition(x,y,z);
      fh<<std::setprecision(8)<<x<<" "<<std::setprecision(8)<<y<<" "<<std::setprecision(8)<<z<<std::endl;
    }
    std::ofstream& fh;
  };


  static const GoSmeshFileRegistration registration_;
};

extern GoFileRegistry FILE_IO_REGISTRY;

class GoSmeshFileRegistration
{
public:
  GoSmeshFileRegistration() {
//    debugmsg("registering "<<GoSmeshFile::className_);
    FILE_IO_REGISTRY.registerCreator(GoSmeshFile::className_,
				     GoSmeshFile::create);
    FILE_IO_REGISTRY.registerExtension(GbString("smesh"),
				       GoSmeshFile::className_);
    FILE_IO_REGISTRY.registerDescription(GoSmeshFile::className_,
					 GbString("Surface Mesh File (*.smesh)"));
  }
};


//#ifndef OUTLINE
//#include "GoSmeshFile.in"
//#endif

#endif // GOSMESHFILE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.2  2002/03/21 14:58:57  elena
| new: changes in dat-file for reading tetrahedral (bugs in element connection)
|
| Revision 1.4  2001/02/13 11:27:08  prkipfer
| adapted to new GoFile::read method
|
| Revision 1.3  2000/11/07 17:27:50  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.2  2000/11/07 09:57:28  uflabsik
| write routine added
|
| Revision 1.1  2000/10/23 09:19:11  prkipfer
| initial version
|
|
+---------------------------------------------------------------------*/
