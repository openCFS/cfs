#ifndef FILE_ELEC3DPDE_2002
#define FILE_ELEC3DPDE_2002

#include "elecPDE.hh" 

namespace CoupledField
{

  //! Class for electrostatic equation in 3D (no adaptivity)
  /*! 
    This class is derived from class ElecPDE. It is used for solving electrostatic equation in 3D. 
  */

class Elec3dPDE: public ElecPDE
{
public:

  //! Constructor. here we read integration parameters
  /*!
    \param 
    \param aGrid pointer to grid
    \param aBCs pointer to Boundary condition object
    \param aMatFile pointer to class Material. material data.
    \param aGrid pointer to class Grid
    \param aInFile pointer to class FileType. input data.
    \param aOutFile  pointer to class WriteResults. output data.
    \param aTimeFunc pointer to class TimeFunc
  */
  Elec3dPDE(Grid * aptgrid, BCs *aptbcs, Material *ptMaterial, TimeFunc *aptTimeFunc, 
	    FileType *aptFileType, WriteResults *aptOut);

  //! Deconstructor
  virtual ~Elec3dPDE();

  //! setup element matrices and hands it over to AlgSys
  /*!
    \param level level of grid
   */
  void SetupMatrices(const Integer level=0);


};

} // end of namespace
#endif
