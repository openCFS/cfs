#ifndef FILE_ACOUSTIC2DPDE_2001
#define FILE_ACOUSTIC2DPDE_2001

#include "acousticPDE.hh"

 
namespace CoupledField
{

  //! Class for acoustic equation (no adaptivity)
  /*! 
    This class is derived from class AcousticPDE. It is used for solving acoustic equation on one time step in 2D.  
  */

class Acoustic2dPDE: public AcousticPDE
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
  Acoustic2dPDE(Grid *aGrid, BCs *aBCs, Material *aMatFile, TimeFunc *aTimeFunc, FileType *aInFile, 
	    WriteResults *aOutFile );

  //!  Deconstructor
  virtual ~Acoustic2dPDE() {;};

  //! specify type of system matrix for AlgebraicSystem
  /*!
    \param level (input) level of Grid
  */
  virtual void SetupMatrices(const Integer level);

  //! compute rhs
  /*!
    \param atime time of calculation
  */
  virtual void ComputeRHS(const Double atime);

  
protected:
  //! some preliminary actions for calculation of RHS. they are the same for all steps. 
  virtual void preComputeRHS();

  //! function for RHS
  Integer arg_rhs_;
  pfn1var ptRHSFnc_;
  std::vector<Double> directionFnc_;

  Boolean SetRHSFnc; //!< Indicator: is there RHS function
  std::vector<std::string> rhs_surfaces_; //!< list of surfaces, on which we have force
};

} // end of namespace
#endif
