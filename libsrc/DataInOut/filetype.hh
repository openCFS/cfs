#ifndef FILE_FILETYPE_2001
#define FILE_FILETYPE_2001

#include "vector.hh"
#include <list>
#include <vector>
#include <string>

#ifdef ADAPTGRID
#include "Element.h"  
#endif

namespace CoupledField
{

  //! base class for reading initial data
  /*! 
    Base class for reading input data. Their functions are virtual due to handle the different types of input files and to hide their features from developer of code.
  */

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
  virtual void ReadCoordinate(Point<3> * const InitNodalCo,                                                     const Integer maxnumNodes)=0;

  //!
  virtual void ReadCoordinate(Point<2> * const InitNodalCo,
			      const Integer maxnumNodes)=0;

  //!
  virtual void ReadBCs(std::list<Integer> * bcs, std::vector<std::string> levels)=0;  

  //!
  virtual void ReadEl(std::vector<Elem*> * elems, const std::vector<std::string>
sd)
 { Error(" not implemented",__FILE__,__LINE__);}

//! read 1D element. we cause it directly when we set BCs
  virtual void ReadEl1d(std::vector<Elem*> * allelems, const std::vector<std::string> sd)
 { Error(" not implemented",__FILE__,__LINE__);}  

#ifdef ADAPTGRID
  virtual void ReadGrid_RG(std::vector<grd::Element*> & elems, std::vector<grd::Vertex*> * vertex, const std::vector<std::string> sd)
   { Error(" not implemented",__FILE__,__LINE__);}

  virtual void ReadBCs_GridRG(std::vector<Integer> & idBCs,std::vector<Integer> &colorBCs)
   { Error(" not implemented",__FILE__,__LINE__);}
#endif
 
  //!
  virtual Integer ReadDim()=0;  

protected:

  //! name of input file
  Char * filename;

};

}


#endif // FILE_FILETYPE
