#ifdef NEWBASEPDE


#ifndef FILE_ELECPDE_NEW
#define FILE_ELECPDE_NEW

#include "newBasePDE.hh" 

namespace CoupledField
{

  //! Class for electrostatic equation in 3D (no adaptivity)
  /*! 
    This class is derived from class BasePDE. It is used for solving electrostatic equation in 3D. 
  */

class ElecPDE : public BasePDE
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
  virtual ~ElecPDE(){};



  //! define all (bilinearform) integrators needed for this pde
  virtual void DefineIntegrators(const Integer level);


  //! reset 
  virtual void Reset();

   //! return size of solution
  virtual Integer getSize() const 
  { return numPDENodes_*dofspernode_;}

// ======================================================
// SOLUTION SECTION
// ======================================================

  //!
  virtual void StepStaticNonLin(const Integer level, const Double aTime=0);

  //!
  virtual void PreStepStatic(const Integer level);

  //!
  virtual void PostStepStatic(const Integer level);


  //! initialize time stepping: nothing to do in electrostatics!
  virtual void InitTimeStepping(const Double dt){;};

  //!
  virtual void SolveStepTrans(const Integer kstep, const Double asteptime,
			      const Integer level, const Boolean updatesysmat)
  {SolveStepStatic(level, asteptime);};

  //!
  virtual void StepTransLin(const Integer level)
  {StepStaticLin(level);};

  //!
  virtual void PreStepTrans(const Integer level)
  {PreStepStatic(level);};

  //!
  virtual void PostStepTrans(const Integer level)
  {PostStepStatic(level);};

  //!
  virtual void StepTransNonLin(const Integer level)
  {StepStaticNonLin(level);};

// ======================================================
// POSTPROCESSING SECTION
// ======================================================


  //! do PostProcessing step
  virtual void PostProcess(const Integer level);

  //! write results in file
  virtual void WriteResultsInFile();

  //! computes the electric energy for each subdomain
  void CalcEnergy();


  //! callculates nodal forces
  void CalcNodeForce(Array<Double> & force, 
		     std::vector<Integer> & nodes, 
		     std::vector<Elem*> & elems,
		     std::vector<std::vector<ShortInt> > & isBoundaryNode,
		     std::vector<std::vector<Integer> > & elemNodeToCouplingNode);


  //! GET SOLUTION AT ALL NODES OF AN ELEMENT
  void GetSolOfElement( Vector<Double>& elpot, Vector<Integer>& connect_PDE);


// ======================================================
// COUPLING SECTION
// ======================================================


  //! initalize PDE coupling
  virtual void InitCoupling(PDECoupling * Coupling);
  
  
  //! calculate coupling terms
  virtual void CalcOutputCoupling();


  //! returns if PDE can compute the quantity
  virtual Boolean HasOutput(std::string output);

  

protected:
  Array<Double> E_;  //!< conatins elecric field
  
  // ---- Electric Force variables ---
  Array<Double> Force_;        //!< stores Electric force of each element
  std::vector<std::vector<Elem*> > F_Interface_; //!<vector of vectors conaining Elements with acting force
  std::vector<std::vector<std::vector<ShortInt> > > isBoundaryNode_; //!< vector containing flag array for element boundary nodes
  std::vector<std::vector<std::vector<Integer> > > elemNodeToCouplingNode_; //!< assigns each coupling element node the according Coupling Node number
  std::vector<std::vector<Integer> > numBoundaryNodes_;               //!< contains number of surface nodes per element

  //postprocessing
  std::vector<std::string> calcEfield_;  //!< contains the subdomains, on which the electric field is computed
  std::vector<std::string> calcEnergy_;  //!< contains the subdomains, on which the electric energy is computed


 // for check: own solver
  Boolean SolverCFS_; //<! parameter indicator: TRUE, if you want to use Solver CFS. reading from config-file
  Matrix<Double> sysmat_;
  Vector<Double> vecrhs_;

 
};

} // end of namespace
#endif

#endif //#ifdef NEWBASEPDE
