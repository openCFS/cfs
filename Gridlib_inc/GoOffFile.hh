/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOOFFFILE_HH
#define GOOFFFILE_HH

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

class GoOffFileRegistration;

class GoOffFile
  : public GoFile
{
public:
  //! Constructor
  GoOffFile(const GbString&);
  virtual ~GoOffFile();

  virtual void write(const GoMesh& m, GbBool smooth, int level=0) const;
  virtual void read(std::vector<GoMesh*>& mv) const;

  static GoFile* create(const GbString& filename);
  static const GbString className_;
  
  //! This operator pretty prints info about the file
//  friend std::ostream& operator<<(std::ostream&, const GoOffFile&);

private:
  //! hide copy constructor
  GoOffFile(const GoOffFile&);

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

  static const GoOffFileRegistration registration_;
};

extern GoFileRegistry FILE_IO_REGISTRY;

class GoOffFileRegistration
{
public:
  GoOffFileRegistration() {
//    debugmsg("registering "<<GoOffFile::className_);
    FILE_IO_REGISTRY.registerCreator(GoOffFile::className_,
				     GoOffFile::create);
    FILE_IO_REGISTRY.registerExtension(GbString("off"),
				       GoOffFile::className_);
    FILE_IO_REGISTRY.registerDescription(GoOffFile::className_,
					 GbString("Object File Format (*.off)"));
  }
};


//#ifndef OUTLINE
//#include "GoOffFile.in"
//#endif

#endif // GOOFFFILE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:57  elena
| new: dir Gridlib_inc
|
| Revision 1.5  2001/02/13 11:28:10  prkipfer
| adapted to new GoFile::read method
|
| Revision 1.4  2000/11/07 17:27:50  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.3  2000/10/23 09:33:07  prkipfer
| added possibility to save subdivision levels
|
| Revision 1.2  2000/07/17 11:22:04  prkipfer
| corrected floating point output precision
|
| Revision 1.1  2000/07/03 16:15:08  prkipfer
| extended IO subsystem, created class registry, added tetrahedron
|
|
+---------------------------------------------------------------------*/
