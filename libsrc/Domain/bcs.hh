#ifndef FILE_BOUNDCOND_2001
#define FILE_BOUNDCOND_2001

#define NUMLEVELGRID 20

#include "filetype.hh"

namespace CoupledField
{
/// class BCs contains information about calculation domain and according to different meshes create different grids

//! struct for Restraints:
struct NodeRestraint
{
  Integer nodalnum;
  enum TypeBCs dof;

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
  Integer numNeumann[NUMLEVELGRID];

  //!    
  Integer numConstraints[NUMLEVELGRID];

  //!
  FileType* InFile_;

  //!
  std::list<NodeRestraint> restr_[NUMLEVELGRID];

  //! number of nodes with restraints
  Integer numrestr_[NUMLEVELGRID];
};

}

#endif // FILE_BCS
