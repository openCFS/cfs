#ifndef FILE_ACOUSTIC3DPDE_2003
#define FILE_ACOUSTIC3DPDE_2003

#include "acousticPDE.hh"

namespace CoupledField
{

  //! Class for acoustic equation (no adaptivity)
  /*! 
    This class is derived from class AcousticPDE. It is used for solving acoustic equation in 3D on one time step.
  */

class Acoustic3dPDE: public AcousticPDE
{
public:

  //!  Constructor. here we read integration parameters
  /*!
    \param aGrid pointer to grid
    \param aBCs pointer to Boundary condition object
    \param aMatFile pointer to class Material. material data.
    \param aInFile pointer to class FileType. input data.
    \param aOutFile  pointer to class WriteResults. output data.
    \param aTimeFunc pointer to class TimeFunc
  */
  Acoustic3dPDE(Grid *aGrid, BCs *aBCs, Material *aMatFile, TimeFunc *aTimeFunc, FileType *aInFile, 
	    WriteResults *aOutFile );

  //! Deconstructor
  virtual ~Acoustic3dPDE();

  //!  setup element matrices and hand it over to AlgebraicSystem 
  /*!
    \param level level of grid
   */
  void SetupMatrices(const Integer level);


protected:

  //! some preliminary actions for calculation of RHS. they are the same for all steps. 
  virtual void preComputeRHS()
  { 
    Error("preComputeRHS not implemenmted in Acoustic3D",__FILE__,__LINE__);
  }

};



} // end of namespace
#endif
