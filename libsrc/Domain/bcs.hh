#ifndef FILE_BOUNDCOND_2001
#define FILE_BOUNDCOND_2001

#include "filetype.hh"

namespace CoupledField
{

/// contains information about calculation domain and according to different meshes create different grids

//! struct for Restraints:

struct NodeRestraint
{
  Integer nodalnum;
  Integer dof;
  Integer numfunc;
  Double factor;

 int operator<(const NodeRestraint & t);
};

class BCs
{
public:
  //!
  BCs(FileType * const InFile);

  //!
  virtual ~BCs();

  //!
  void ReadBCs();

  //!
  Integer GetNumRestraints(const Integer level) 
  { return numrestr_[level]; }

  //!
  void GetRestraints(std::list<NodeRestraint> & arestr, const Integer level)
  { arestr=restr_[level]; }
 
protected:

private:
  //!
  Integer numNeumann[20];

  //!    
  Integer numConstraints[20];

  //!
  FileType* InFile_;

  //!
  std::list<NodeRestraint> restr_[20];

  //! number of nodes with restraints
  Integer numrestr_[20];
};

}

#endif // FILE_BCS
