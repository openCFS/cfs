#ifndef FILE_FILETYPE_2001
#define FILE_FILETYPE_2001

#include "vector.hh"
#include <list>

namespace CoupledField
{

  //! base class for reading initial data
  /*! 
    Base class for reading input data. Its function is virtual due to handle the different types of input files and to hide their features from developer of code.
  */

struct NodeRestraint;

class FileType
{

public:

  //! enum for reading boundary condition
  enum nameBound{ numdofs, numconstr, numrestr, numloads, resistors,
		  numspring, bembdry, numflux, numrad, numpress, ncurrdens,
		  regl_itmx, reglr_use, endBound}; // 13   

  //! enum for reading data about general analisys
  enum nameGAn{soltype, analtype, numnode, numgroup, restart, inactdofs,
	       circuit, deactDf, masstype, damptype, nooptimiz, endGAnal}; //11

  //! enum for redaing general info about element
  enum nameGElem {numelem, ielemtyp, isubtype, ielemsave, maxnode,
		  nonlinear, form1, form2, endGElem}; // 8 

  //!
  FileType(const Char * const afilename);

  //!
  virtual ~FileType();

  //!
  virtual void ReadMaxnumnodes(Integer & maxnumnodes)=0;  

  //!
  virtual void ReadMaxnumelem(Integer & maxnumelem)=0;

  //!
  virtual void ReadNumberNodesPerElem(Integer &)=0;

  //!
  virtual void ReadElemConnectionGH(const Integer maxelem, Integer * connect, const Integer maxnode, const Integer numelemgr, const Integer startposinarrayconn)=0;

  //!
  virtual void ReadCoordinate(Point3D * const InitNodalCo,                                                     const Integer maxnumNodes)=0;

  //!
  virtual void ReadCoordinate(Point2D * const InitNodalCo,
			      const Integer maxnumNodes)=0;

  //!
  virtual void ReadBoundRestr(std::list<NodeRestraint> & restr, Integer & numberRestr)=0;

protected:

  //! name of input file
  Char * filename;
};

}


#endif // FILE_FILETYPE
