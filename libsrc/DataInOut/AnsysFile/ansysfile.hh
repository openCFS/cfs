#ifndef FILE_ANSYSFILE_2002
#define FILE_ANSYSFILE_2002

#include "filetype.hh"

namespace CoupledField
{

  //! base class for reading initial data
  /*! 
    Class, that is derived from class FileType for reading input data, which is produced by Ansys interface. 
  */

class AnsysFile: virtual public FileType 
{

public:

  //!
  AnsysFile(const Char * const afilename);

  //!
  virtual ~AnsysFile();

  //!
  virtual void ReadMaxnumnodes(Integer & );

  //!
  virtual void ReadMaxnumelem(Integer & );

 void ReadNumberNodesPerElem(Integer & nnodesperelem);

  //!
  virtual void ReadElemConnectionGH(const Integer maxelem, Integer * connect, const Integer maxnode, const Integer numelemgr, const Integer startposinarrayconn);

  //!
  virtual void ReadCoordinate(Point3D * const NodesCoord,                                                     const Integer maxnumNodes);

  //!
  virtual void ReadCoordinate(Point2D * const NodesCoord,
			      const Integer maxnumNodes);

  //!
  virtual void ReadBoundRestr(std::list<NodeRestraint> & restr, Integer & numberRestr);

protected:
  //! 
  std::ifstream infile;
  
private:

  std::string::size_type pos_end;

  // dimension of problem
  Integer dim_;

  //!
  void ReadDim();

  // get position after line with seekexp and comments lines
  void getPosLine(const std::string seekexp, std::string::size_type & pos);

  // get position in line
void getPosition(const std::string seekexp, std::string::size_type & pos);

  // transform string for level of boundary condition in number
  Integer TransformInDof(const std::string type_bc);

  // read number of nodes for boundary condition
  void ReadMaxnumnodesbc(Integer & nbc);

};

}

#endif // FILE_FILETYPE
