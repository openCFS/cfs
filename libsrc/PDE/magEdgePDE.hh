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
		     std::vector<Integer>& esign, BaseFE * ptElem);

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

  //! return pointer to (real) vector with solution
  virtual const Vector<Double> & getS() const { return solRe_;}

  //! return pointer to imaginary vector with  solution
  virtual const Vector<Double> & getSIm() const { return solIm_;}

  //! return size of solution
  virtual Integer getSize() const { return size_;}

  /// calc element quantities
  void PostProcess(const Integer level);

  /// solve one harmonic step
  virtual void SolveStepHarmonic(const Integer level);
  
  

protected:
  /// setup source term
  void SetupRHS(const Integer level);

  /// reads all data in the config-file belonging to coils
  void MagEdgePDE::ReadCoils();
  
  /// correct the sign in the element matrix due to orientation of edges
  /*!
    \param elemmat Element matrix
    \param esign   Vector of edge orientations
  */
  void CorrectEdgeDir(Matrix<Double> & elemmat, std::vector<Integer>& esign);
  


  Vector<Double> solRe_;      //!< store real part of solution,
  Vector<Double> solIm_;      //!< store imaginary part of solution
  Integer size_;              //!< size of solution (number of edges)
  Integer NumElems_;          //!< number of elements belonging to PDE
  Integer * EdgeDir_;         //!< stores the Dirichlet-edges
  Integer numEdgedir_;        //!< number of Dirichlet edges
  Vector<Double>* bFieldRe_;  //!< real part of vector of magnetic field induction
  Vector<Double>* bFieldIm_;  //!< imaginary part of vector of magnetic field induction
  Vector<Double>* magVecPotRe_; //!< real part of vector of magnetic magnetic vector potential
  Vector<Double>* magVecPotIm_; //!< imaginary part of vector of magnetic magnetic vector potential
  std::vector <std::string> coilDomain_;  //!< name of all subdomains containing coils
  Double freq_;               //!< excitation frequency for harmonic analysis
  /// parameters necessary to describe coils
  struct coilDefStruct
  {
    Integer iDir;
    Double  current;
    Double  coilArea;
    Double currentPhase;
    std::vector<Double> coilMidPt;
  };
  

  std::vector<struct coilDefStruct> coilDef_; //!< vector of paramters describing coils
  
};

} // end of namespace
#endif
