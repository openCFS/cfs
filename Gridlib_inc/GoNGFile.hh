/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/
#ifndef GONGFILE_HH
#define GONGFILE_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include "GoFile.hh"
#include "GoFileRegistry.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

class GoNGFileRegistration;

class GoNGFile
  : public GoFile
{
public:
  //! Constructor
  GoNGFile(const GbString&);
  virtual ~GoNGFile();

  virtual void read(std::vector<GoMesh*>& mv) const;
  virtual void write(const GoMesh& m, GbBool smooth, int level=0) const;

  static GoFile* create(const GbString& filename);
  static const GbString className_;
  
  //! This operator pretty prints info about the file
//  friend std::ostream& operator<<(std::ostream&, const GoNGFile&);

private:
  //! hide copy constructor
  GoNGFile(const GoNGFile&);

  struct AddTetra {
    AddTetra(std::ofstream& f) : fh(f) {}
    void operator()(GoGeometryElement<float> *t) {
      fh<<t->getPartition()+1<<" 4 ";
      fh<<t->getVertex(0)->getId()+1<<" ";
      fh<<t->getVertex(2)->getId()+1<<" ";
      fh<<t->getVertex(1)->getId()+1<<" ";
      fh<<t->getVertex(3)->getId()+1<<" ";
      fh<<std::endl;
    }
    std::ofstream& fh;
  };

  struct AddCoordinate {
    AddCoordinate(std::ofstream& f) : fh(f) {}
    void operator()(GoVertex<float> *t) {
      float x,y,z;
      t->getPosition(x,y,z);
      fh.setf(std::ios::fixed, std::ios::floatfield);
      //fh<<std::setprecision(8)<<x<<" "<<std::setprecision(8)<<y<<" "<<std::setprecision(8)<<z<<std::endl;
      fh<<x<<" "<<y<<" "<<z<<std::endl;
      
    }
    std::ofstream& fh;
  };
  


  static const GoNGFileRegistration registration_;
};

extern GoFileRegistry FILE_IO_REGISTRY;

class GoNGFileRegistration
{
public:
  GoNGFileRegistration() {
//    debugmsg("registering "<<GoNGFile::className_);
    FILE_IO_REGISTRY.registerCreator(GoNGFile::className_,
				     GoNGFile::create);
    FILE_IO_REGISTRY.registerExtension(GbString("vol"),
				       GoNGFile::className_);
    FILE_IO_REGISTRY.registerExtension(GbString("N"),
				       GoNGFile::className_);
    FILE_IO_REGISTRY.registerDescription(GoNGFile::className_,
					 GbString("NetGen Output (*.vol)"));
  }
};


//#ifndef OUTLINE
//#include "GoNGFile.in"
//#endif

#endif // GONGFILE_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.3  2002/04/03 11:17:08  elena
| new: changes in CFS++ on SGI and new function for gridlib
|
| Revision 1.8  2001/09/12 09:34:41  prkipfer
| fixed FORTRAN count bug and adapted to new GoMesh
|
| Revision 1.7  2001/02/13 11:29:25  prkipfer
| adapted to new GoFile::read method
|
| Revision 1.6  2000/11/15 11:55:53  uflabsik
| changed order of vertices in write routine
|
| Revision 1.5  2000/11/07 17:27:49  prkipfer
| introduced external polymorphism for vertices
|
| Revision 1.4  2000/11/07 09:57:06  uflabsik
| registry entry added for .N-files
|
| Revision 1.3  2000/10/23 09:33:06  prkipfer
| added possibility to save subdivision levels
|
| Revision 1.2  2000/09/03 14:04:15  uflabsik
| output inserted
|
| Revision 1.1  2000/07/03 16:15:06  prkipfer
| extended IO subsystem, created class registry, added tetrahedron
|
|
+---------------------------------------------------------------------*/
