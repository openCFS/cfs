#ifndef FILE_ACOUSTICPDE_2001
#define FILE_ACOUSTICPDE_2001

#include "newBasePDE.hh"

 
namespace CoupledField
{

  //! Class for acoustic equation (no adaptivity)
  /*! 
    This class is derived from class BasePDE. It is used for solving acoustic equation on one time step.  
  */

class AcousticPDE: public BasePDE
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
  AcousticPDE(Grid *aGrid, BCs *aBCs, TimeFunc *aTimeFunc, FileType *aInFile, 
	      WriteResults *aOutFile );

  //!  Deconstructor
  virtual ~AcousticPDE() {;};


  //! define all (bilinearform) integrators needed for this pde
  virtual void DefineIntegrators(const Integer level);

  //! prepare for correct time stepping
  /*!
    \param dt time step
  */
  virtual void InitTimeStepping(const Double dt);

  //! write results in file
   virtual void WriteResultsInFile();

  //!  return pointer to vector with first derivative of solution
  virtual const Array<Double>& getS1() const { return TS_alg_->GetDeriv1();}

  //! return pointer to vector with second derivative of solution
  virtual const Array<Double>& getS2() const { return TS_alg_->GetDeriv2();}

  //! return size of solution
  virtual Integer getSize() const 
  { return size_;}

#ifdef ADAPTGRID
  //! test error of computation
  virtual Boolean TestError(const Integer level);
#endif


  // ======================================================
  // COUPLING SECTION
  // ======================================================
  
  //! initalize PDE coupling
  void InitCoupling(PDECoupling * Coupling);;
  
  //! calculate coupling terms
  void CalcOutputCoupling();

  //! returns if PDE can compute the quantity
  virtual Boolean HasOutput(std::string output);
  
  //! calculate the vector of coupling forces to the mechanical PDE
  void CalcMechCouplingRHS(std::vector<Elem*> * couplingElems, 
			   std::vector<Integer> & couplingNodes,
			   std::vector<MaterialData*> * couplingMaterials,
			   Array<Double>& elemCouplingSols,
			   Integer couplingdof);
  

protected:

  Double freq_;   //!< excitation frequency for harmonic analysis
  Array<Double> solIm_; //!< stores the imaginary part in case of harmonic analysis

  //  Double lasttimecalc_;  //!< Last time on which we have calculated solution
  //  Integer laststepcalc_; //!< Number of last timestep on which we have calculated our solution

  Integer size_; //!< total number of unknowns (equations)

  Boolean with_absBCs_; //!< Indicator for absorbing boundary conditions 
  std::vector<std::string> bnd_absBCs_;   //!< list of bnds( for absorbing BCs)

  DampingType damping_type_; //!< specifies the type of damping model (see environment.hh)
  Boolean with_fracdamping_; //!< attenuation according to power law
  Integer frac_memory_;      //!< number of old time steps to be saved

  //General dimension of problem
  Integer dim_;
  
};

} // end of namespace
#endif
