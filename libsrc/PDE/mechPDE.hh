#ifndef FILE_BASEMECHPDE
#define FILE_BASEMECHPDE

#include "basepde.hh"

 
namespace CoupledField
{

  //! Class for mechanic equation (no adaptivity)
  /*! 
    This class is derived from class BasePDE. 
    It is used for solving mechanic equation on one time step.  
  */

class MechPDE: public BasePDE
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
  MechPDE(Grid *aGrid, BCs *aBCs, TimeFunc *aTimeFunc, FileType *aInFile, WriteResults *aOutFile );

  //!  Deconstructor
  virtual ~MechPDE() {;};

  //! define discrete PDE
  virtual void DiscreteParamsPDE();

  //! set information for algebraic system about PDE. set matrix factors
  virtual void SetMatrixFactors();

 //! initalize PDE coupling
  virtual void InitCoupling(PDECoupling * Coupling);
  
  //! specify type of system matrix for AlgebraicSystem
  /*!
    \param level (input) level of Grid
  */
  virtual void SetupMatrices(const Integer level);


    //! set boundary condition
  /*!
    \param level level of grid
    \param update indicator: do we update boundary condition in algebraic system or set new
    \param atimestep time step of calculation
  */
  virtual void SetBCs(const Integer level, const Integer update, const Double atimestep);


  //! compute rhs
  /*!
    \param atime time of calculation
  */
  virtual void ComputeRHS(const Double atime) {;};
  
  virtual void SolveStepStatic(const Integer level);
  

  //! solve one step for transient problem 
  /*!
    \param kstep number of calculating step
    \param steptime time of calculation
    \param level level of grid
    \param updatesysmat indicator: need we to update algebraic system. it is used for adaptive procedure in space
  */
  virtual void SolveStepTrans(const Integer kstep, const Double steptime, const Integer level, 
			      const Boolean updatesysmat)
  { 
    Error("Currently not available",__FILE__,__LINE__);
  }

  //! calculate coupling terms
  virtual void CalcOutputCoupling();
  
  //! write results in file
   virtual void WriteResultsInFile();

  //! returns if PDE can compute the quantity
  virtual Boolean HasOutput(std::string output);

  //! Assemble mass part
  void AssembleMass(BaseFE * ptEl, Vector<Integer>& connect_PDE,  Matrix<Double>& ptCoord, Double density);

  //! Assemble stiffness part
  void AssembleStiffness(BaseFE * ptEl, Vector<Integer>& connect_PDE,  Matrix<Double>& ptCoord, MaterialData& actMatData);
  

protected:
  
  
  Integer size_;        //!< total number of unknowns (equations)


private:
  // defines subtype of mechanic PDE: plainStrain, 3d, ...
  std::string subType_;

  Integer GetNrBCDof (const std::string & dofStartString);

  /// flag for nonlinear calculations
  Boolean nonLin_;
  
  
};

} // end of namespace
#endif
