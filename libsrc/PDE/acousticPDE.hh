#ifndef FILE_ACOUSTICPDE_2001
#define FILE_ACOUSTICPDE_2001

#include "basePDE.hh"

 
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
				
  //! write results in file
  //! \param stepOffset offset for starting (time)step
  //! \param timeOffset offset for starting time  
  virtual void WriteResultsInFile(Integer stepOffset = 0,
				  Double timeOffset = 0.0);

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
  virtual Boolean HasOutput(SolutionType output);
  
  //! calculate the vector of coupling forces to the mechanical PDE
  void CalcMechCouplingRHS(StdVector<Elem*> * couplingElems, 
			   StdVector<Integer> & couplingNodes,
			   StdVector<MaterialData*> * couplingMaterials,
			   Vector<Double>& elemCouplingSols,
			   Integer couplingdof,
			   StdVector<Elem*> * neighbours);
  


protected:

  //! Init the time stepping
  void InitTimeStepping();

  // Double freq_;   //!< excitation frequency for harmonic analysis
  NodeStoreSol<Double> sol_der1Array_, sol_der2Array_;

  Integer size_; //!< total number of unknowns (equations)

  StdVector<std::string> absBCs_;   //!< list of boundaries( for absorbing BCs)
  Boolean absorbingBCs_;               //!< switch for absorbing boundary conditions
  //DampingType dampingType_; //!< specifies the type of damping model (see environment.hh)
  Integer fracMemory_;      //!< number of old time steps to be saved

  // Postprocessing results

  //! contains 1. derivative of acoustic potential
  NodeStoreSol<Double> solDeriv1_;
  
  //! contains 2. derivative of acoustic potential
  NodeStoreSol<Double> solDeriv2_;

private:

#ifdef XMLPARAMS
    //! Obtain information on desired output quantities from parameter file

    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    //! The output quantities currently supported by the acoustics PDE are
    //! given in the following table. Here 'Keyword' and 'Result Type' refer
    //! to the XML parameter file, while 'Class Attribute' refers to the
    //! internal attribute of the AcousticPDE class that is set, if the keyword
    //! is specified.\n
    //! \todo Specification of ReadStoreResults for AcousticPDE!!!
    void ReadStoreResults();
#endif

};

} // end of namespace
#endif
