/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GOFKRFILE_HH
#define GOFKRFILE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <iomanip>

#include "GoFile.hh"
#include "GoFileRegistry.hh"

#include "GoFlowVertex.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GoFkrFileRegistration;

class GoFkrFile
  : public GoFile
{
public:
  //! Constructor
  GoFkrFile(const GbString&);
  virtual ~GoFkrFile();

  virtual void write(const GoMesh& m, GbBool smooth, int level=0) const;
  virtual void read(std::vector<GoMesh*>& mv) const;
  void readBin(std::vector<GoMesh*>& mv) const;

  static GoFile* create(const GbString& filename);
  static const GbString className_;
  
  //! This operator pretty prints info about the file
//  friend std::ostream& operator<<(std::ostream&, const GoFkrFile&);

private:
  //! hide copy constructor
  GoFkrFile(const GoFkrFile&);

  struct AddTet {
    AddTet(std::ofstream& f) : fh(f) {}
    void operator()(GoGeometryElement<float> *t) {
      for (int i=0; i<t->getNumVertices();i++) 
	fh<<t->getVertex(i)->getId()+1<<" ";
      fh<<std::endl;
    }
    std::ofstream& fh;
  };

  struct AddCoordinate {
    AddCoordinate(std::ofstream& f) : fh(f) {}
    void operator()(GoVertex<float> *t) {
      float x,y,z,d,e;
#ifdef HAS_DYN_CAST
      GoFlowVertex *vert = dynamic_cast<GoFlowVertex*>(t);
#else
      GoFlowVertex *vert = static_cast<GoFlowVertex*>(t);
#endif
      vert->getPosition(x,y,z);
      d=vert->getDensity();
      e=vert->getEnergy();
      GbVec3<float> m=vert->getMomentum();
      fh<<std::fixed<<std::setprecision(8)<<x<<" "<<std::fixed<<std::setprecision(8)<<y<<" "<<std::fixed<<std::setprecision(8)<<z<<" ";
      fh<<std::fixed<<std::setprecision(8)<<d<<" ";
      fh<<std::fixed<<std::setprecision(8)<<m[0]<<" "<<std::fixed<<std::setprecision(8)<<m[1]<<" "<<std::fixed<<std::setprecision(8)<<m[2]<<" ";
      fh<<std::fixed<<std::setprecision(8)<<e<<std::endl;
    }
    std::ofstream& fh;
  };

  static const GoFkrFileRegistration registration_;
};

extern GoFileRegistry FILE_IO_REGISTRY;

class GoFkrFileRegistration
{
public:
  GoFkrFileRegistration() {
//    debugmsg("registering "<<GoFkrFile::className_);
    FILE_IO_REGISTRY.registerCreator(GoFkrFile::className_,
				     GoFkrFile::create);
    FILE_IO_REGISTRY.registerExtension(GbString("fkr"),
				       GoFkrFile::className_);
    FILE_IO_REGISTRY.registerDescription(GoFkrFile::className_,
					 GbString("Frank Reck File (*.fkr)"));
  }
};


//#ifndef OUTLINE
//#include "GoFkrFile.in"
//#endif

#endif // GOFKRFILE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.6  2001/12/11 12:42:30  prkipfer
| fixes for KCC compiler on new PCs
|
| Revision 1.5  2001/06/18 11:06:18  prkipfer
| update binary format reading routine
|
| Revision 1.4  2001/04/17 12:38:06  prkipfer
| output is now fixed precision
|
| Revision 1.3  2001/02/13 11:30:32  prkipfer
| adapted to new GoFile::read method
|
| Revision 1.2  2001/01/02 14:26:23  prkipfer
| added binary file reader
|
| Revision 1.1  2000/11/07 17:27:48  prkipfer
| introduced external polymorphism for vertices
|
|
+---------------------------------------------------------------------*/
