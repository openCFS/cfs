#ifndef FILE_BOUNDCOND_2001
#define FILE_BOUNDCOND_2001

#define NUMLEVELGRID 20

#include "filetype.hh"
#include "grid.hh"

namespace CoupledField
{
/// class BCs contains information about calculation domain and according to different meshes create different grids

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
  void Update(Grid * ptgrid);  

  //!
  std::list<Integer> GetNodesLevel(const std::string level, const Integer lev);  
  //!
  Integer GetNumNodesLevel(const std::string level, const Integer lev);

protected:

private:

   std::vector<std::string> levels_;

   std::list<Integer> * bcs_[NUMLEVELGRID];

  //!
   FileType* InFile_;
};

}

#endif // FILE_BCS
