#ifndef FILE_MECH2DPDE
#define FILE_MECH2DPDE

#include "mechPDE.hh"

 
namespace CoupledField
{

  //! Class for mechanic equation (no adaptivity)
  /*! 
    This class is derived from class MechPDE. 
    It is used for solving mechanic equation on one time step in 2D.  
  */

class PlainStrainPDE : public MechPDE
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
  PlainStrainPDE(Grid *aGrid, BCs *aBCs, TimeFunc *aTimeFunc, FileType *aInFile, 
	    WriteResults *aOutFile );

  //!  Deconstructor
  virtual ~PlainStrainPDE() {;};

  //! specify type of system matrix for AlgebraicSystem
  /*!
    \param level (input) level of Grid
  */
  virtual void SetupMatrices(const Integer level);

    //! set boundary condition
  /*!
    \param level level of grid
    \param update indicator: do we update boundary condition in algebraic system ot set new
    \param atimestep time step of claculation
  */
  virtual void SetBCs(const Integer level, const Integer update, const Double atimestep);

protected:

};

} // end of namespace
#endif
