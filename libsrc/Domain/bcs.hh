#ifndef FILE_BOUNDCOND_2001
#define FILE_BOUNDCOND_2001

#define NUMLEVELGRID 20

#include <DataInOut/filetype.hh>
#include "grid.hh"

struct Elem;

namespace CoupledField
{
//! class BCs contains information about boundary and interface data

class BCs
{
public:
  //! constructor
  /*!
    \param InFile pointer to input file
  */
  BCs(FileType * const InFile);

  //! destructor
  virtual ~BCs();

  //! reads boundary conditions out of file
  void ReadBCs();

  //! 
  void Update(Grid * ptgrid);  

  //! returns all bc-nodes belonging to the requested level
  /*!
    \param level name of boundary type
    \param lev index for hierarchy in multilevel method
  */
  std::list<Integer> GetNodesLevel(const std::string level, const Integer lev=-1);  

  //! returns all 1D elements belonging to the requested interface
  /*!
    \param color name of interface
    \param lev index for hierarchy in multilevel method
  */
  std::vector<Elem*> BCs::getEdgesBC(const std::string color, const Integer lev=-1);

//! returns all 2D elements belonging to the requested interface
  /*!
    \param color name of interface
    \param lev index for hierarchy in multilevel method
  */
  std::vector<Elem*> BCs::getFacesBC(const std::string color, const Integer lev=-1);

  //! get number of nodes belonging to specified level and hierarchy index
  /*!
    \param level name of boundary type
    \param lev index for hierarchy in multilevel method
  */
  Integer GetNumNodesLevel(const std::string level, const Integer lev=-1);

  //! get 2D-neighbors elems for set of 1D-elements
  /*!
    \param color (input) name of boundary elements 
    \param lev (input)  index for hierarchy in multilevel method
  */
  std::vector<Elem*> getNeighElemsForSurfaces(const std::string color, const Integer lev=-1);

  //! 
  //! prints boundary data to screen
  void BCs :: printBCs(const Integer alevel=-1);

protected:

private:

  std::vector<std::string> levels_;            //!< stores all names of boundary nodes
  std::vector<std::string> color_edges_;       //!< stores all names of 1D interfaces
   std::vector<std::string> color_faces_;       //!< stores all names of 2D interfaces
  std::list<Integer> * bcs_[NUMLEVELGRID];     //!< stores all boundary nodes
  std::vector<Elem*> * bcsEdges_[NUMLEVELGRID]; //!< stores all 1D elements along interfaces
   std::vector<Elem*> * bcsFaces_[NUMLEVELGRID]; //!< stores all 2D elements along interfaces
  //! color of neighbors elements. it is got from conf-file
  std::vector<std::string> color_neighelems_;
  //! array with 2D neighbors elements for 1D;
  std::vector<Elem*> * bcsNeighElems_[NUMLEVELGRID];

  FileType* InFile_; //!< pointer to input file
  Integer toplevel_; //!< stores current top level
};

}

#endif // FILE_BCS
