#ifndef FILE_ELEC2DPDE_2001
#define FILE_ELEC2DPDE_2001

#include "elecPDE.hh"

 
namespace CoupledField
{

  //! Class for electrostatic 2D (no adaptivity)
  /*! 
    This class is derived from class ElecPDE. It is used for solving electrostatic equation.  
  */

class Elec2dPDE: public ElecPDE
{
public:
    
  //!  Constructor. here we read integration parameters
  /*!
    \param aGrid pointer to grid
    \param aBCs pointer to Boundary condition object
    \param aInFile pointer to class FileType. input data.
    \param aOutFile  pointer to class WriteResults. output data.
    \param aTimeFunc pointer to class TimeFunc
  */
  Elec2dPDE(Grid *aGrid, BCs *aBCs, TimeFunc *aTimeFunc, FileType *aInFile, WriteResults *aOutFile );

  //!  Deconstructor
  virtual ~Elec2dPDE();



  //! set up the element matrices and hands it over to AlgSys
  /*!
    \param level (input) level of Grid
  */
  virtual void SetupMatrices(const Integer level);

};

} // end of namespace
#endif
