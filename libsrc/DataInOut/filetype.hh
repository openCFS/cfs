#ifndef FILE_FILETYPE_2001
#define FILE_FILETYPE_2001

#include "vector.hh"
#include <list>
#include <vector>

namespace CoupledField
{

  //! base class for reading initial data
  /*! 
    Base class for reading input data. Their functions are virtual due to handle the different types of input files and to hide their features from developer of code.
  */

struct NodeRestraint;
struct Elem;

class FileType
{

public:

  //!
  FileType(const Char * const afilename);

  //!
  virtual ~FileType();

  //!
  virtual void ReadMaxnumnodes(Integer & maxnumnodes)=0;  

  //!
  virtual void ReadMaxnumelem(Integer & maxnumelem)=0;

  //!
  virtual void ReadCoordinate(Point3D * const InitNodalCo,                                                     const Integer maxnumNodes)=0;

  //!
  virtual void ReadCoordinate(Point2D * const InitNodalCo,
			      const Integer maxnumNodes)=0;

  //!
  virtual void ReadBoundRestr(std::list<NodeRestraint> & restr, Integer & numberRestr)=0;

  //!
  virtual void ReadElems(std::vector<Elem> & allelems)=0;
 
  //!
  virtual Integer ReadDim()=0;  

protected:

  //! name of input file
  Char * filename;

};

}


#endif // FILE_FILETYPE
