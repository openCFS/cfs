#ifndef FILE_ELECPDE_2002
#define FILE_ELECPDE_2002

#include "basepde.hh" 

namespace CoupledField
{

  //! Class for electrostatic equation in 3D (no adaptivity)
  /*! 
    This class is derived from class BasePDE. It is used for solving electrostatic equation in 3D. 
  */

class ElecPDE: public BasePDE
{
public:

  //! Constructor. here we read integration parameters
  /*!
    \param 
    \param aGrid pointer to grid
    \param aBCs pointer to Boundary condition object
    \param aGrid pointer to class Grid
    \param aInFile pointer to class FileType. input data.
    \param aOutFile  pointer to class WriteResults. output data.
    \param aTimeFunc pointer to class TimeFunc
  */
  ElecPDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, WriteResults *aptOut);

  //! Deconstructor
  virtual ~ElecPDE();

  //! define discrete PDE
  virtual void DiscreteParamsPDE();

  //! set information for algebraic system about PDE. set matrix factors.
  void SetMatrixFactors();

  //! setup element matrices for AlgebraicSystem for assembling procedure
  /*!
    \param level level of grid
   */
  void SetupMatrices(const Integer level=0);
  


  //! solve one step for static problem 
  /*!
    \param ptBCs pointer to class with data about boundary condition
    \param level level of grid
  */
  virtual void SolveStepStatic(const Integer level);

  virtual void SolveStepTrans(const Integer kstep, const Double steptime, const Integer level, 
			      const Boolean updatesysmat)
  { 
    Error("Makes no sense for Electrostatics to perform transient step",__FILE__,__LINE__);
  }

  //! write results in file
  virtual void WriteResultsInFile();

  //! return pointer to vector with solution
  virtual const Vector<Double> & getS() const { return sol_;}

  //! return size of solution
  virtual Integer getSize() const { return size_;}

protected:

  Vector<Double> sol_;  //!< store solution,
  Integer size_;        //! size of solution (number of equations)

};

} // end of namespace
#endif
