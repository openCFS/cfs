#ifndef FILE_SMOOTHPDE
#define FILE_SMOOTHPDE

#include "basepde.hh"

 
namespace CoupledField
{

  //! Class for mechanic equation (no adaptivity)
  /*! 
    This class is derived from class BasePDE. 
    It is used for solving mechanic equation on one time step.  
  */

class SmoothPDE: public BasePDE
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
  SmoothPDE(Grid *aGrid, BCs *aBCs, TimeFunc *aTimeFunc, FileType *aInFile, WriteResults *aOutFile );

  //!  Deconstructor
  virtual ~SmoothPDE() {;};

  //! define all (bilinearform) integrators needed for this pde
  virtual void DefineIntegrators(const Integer level);

  //! initalize PDE coupling
  virtual void InitCoupling(PDECoupling * Coupling);

  //! perform ..
  virtual void PreStepStatic(const Integer level);

  //!
  virtual void StepStaticNonLin(const Integer level);

  virtual void SolveStepTrans(const Integer kstep, const Double steptime, const Integer level, 
			      const Boolean updatesysmat)
  {Error("Currently not available",__FILE__,__LINE__);}

  //! 
  virtual void PostStepStatic(const Integer level);

  //! calculate coupling terms
  virtual void CalcOutputCoupling();
  
  //! write results in file
   virtual void WriteResultsInFile();

   //! returns if PDE can compute the quantity
  virtual Boolean HasOutput(std::string output);
  
protected:

  Integer size_;        //!< total number of unknowns (equations)

private:

  // defines subtype of mechanic PDE: plainStrain, 3d, ...
  std::string subType_;
  
  Integer GetNrBCDof (const std::string & dofStartString);

  Integer GetBCDof(const std::string dofString);

  std::string method_;

  Boolean firstTurn_;

  Vector<Double> factor_;

};

} // end of namespace
#endif
