#ifndef FILE_MAGEDGEPDE_2002
#define FILE_MAGEDGEPDE_2002

#include "basepde.hh" 

namespace CoupledField
{

  //! Class for electrostatic equation in 3D (no adaptivity)
  /*! 
    This class is derived from class BasePDE. It is used for solving electrostatic equation in 3D. 
  */

class MagEdgePDE: public BasePDE
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
  MagEdgePDE(Grid * aptgrid, BCs *aptbcs, TimeFunc *aptTimeFunc, FileType *aptFileType, WriteResults *aptOut);

  //! Deconstructor
  virtual ~MagEdgePDE();

  //! define discrete PDE
  virtual void DiscreteParamsPDE();

  //! set information for algebraic system about PDE. set matrix factors.
  void SetMatrixFactors();


  //! set up the edge graph
  void SetupEdgeGraph();

  //! define algebraic system 
  /*!
    \param AS_sysid id of PDE in algebraic system
   */
  virtual void SetAlgSys(const Integer AS_sysid);

  //! Create the matrices and Solver as well as Preconditioner
  virtual void CreateMatrices_Solver();

  //! compute the number of edge) dirichlets
  void EvalNumDirichlet();

  //! get the edge numbers and theis signs (from AlgSys)
  /*!
    \param nodes node numbers of the element
    \param epos integer array conatining the edge numbers
    \param esign integer array containing the edge signs
  */
  void GetEdgeNumber(Integer *nodes, std::vector<Integer>& epos, 
		     std::vector<Integer>& esign);

  //! setup element matrices for AlgebraicSystem for assembling procedure
  /*!
    \param level level of grid
   */
  void SetupMatrices(const Integer level=0);
  
  //! set boundary condition
  /*!
    \param level level of grid
    \param update indicator: do we update boundary condition in algebraic system ot set new
    \param atimestep time step of claculation
  */
  virtual void SetBCs(const Integer level, const Integer update, const Double atimestep);

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

  /// calc element quantities
  void PostProcess(const Integer level);
  

protected:
  /// setup source term
  void SetupRHS(const Integer level);
  


  Vector<Double> sol_;  //!< store solution,
  Integer size_;        //!< size of solution (number of edges)
  Integer NumElems_;    //!< number of elements belonging to PDE
  Integer * EdgeDir_;   //!< stores the Dirichlet-edges
  Integer numEdgedir_;  //!< number of Dirichlet edges
  Vector<Double>* bField_; //!< vector of magnetic field induction
  
};

} // end of namespace
#endif
