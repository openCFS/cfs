#ifndef FILE_BCS_2001
#define FILE_BCS_2001

#include "filetype.hh"

namespace CoupledField
{

/// contains information about calculation domain and according to different meshes create different grids

//! struct for Restraints:
struct Restraints
{
  Integer num;
  Integer ** info;
  Double * val;
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
  Integer GetNumDirichlet(Integer level) 
  { return numDirichlet[level]; }

  //!
  Integer GetDirichletNode(Integer pos, Integer level)
  { return restr[level].info[pos][0];}

  //!
  Double GetDirichletVal(Integer pos, Integer level)
  { return restr[level].val[pos]; }

  Integer GetDirichletDof(Integer pos, Integer level)
  { return restr[level].info[pos][1];}

  Integer GetDirichletTfunc(Integer pos, Integer level)
  { return restr[level].info[pos][2];}

protected:

private:
  //!
  Integer level;

  //!
  Integer numDirichlet[20];

  //!
  Integer numNeumann[20];

  //!    
  Integer numConstraints[20];

  //!
  FileType *InFile;

  //!
  Restraints restr[20];

};

}

#endif // FILE_BCS
