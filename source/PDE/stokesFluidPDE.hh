#ifndef FILE_NEWBASESTOKESFLUIDPDE
#define FILE_NEWBASESTOKESFLUIDPDE

#include "SinglePDE.hh"
 
namespace CoupledField
{

  //! Forward class declarations
  class BaseForm;


  //! Class for mechanic equation (no adaptivity)
  class StokesFluidPDE: public SinglePDE
  {

  public:

    //!  Constructor. here we read integration parameters
    /*!
      \param aGrid pointer to grid
      \param aOutFile  pointer to class WriteResults. output data.
    */
    StokesFluidPDE( Grid* aGrid, ParamNode* paramNode );

    //!  Deconstructor
    virtual ~StokesFluidPDE();

    //! Initialize NonLinearities
    virtual void InitNonLin();

    //! define all (bilinearform) integrators needed for this pde
    virtual void DefineIntegrators( );

    //! define the SoltionStep-Driver
    virtual void DefineSolveStep();

    /// returns a stiffness integrator appropriate to the actual problem (e.g.3D)
    BaseForm * GetStiffIntegrator(Double density, Double dynamicViscosity);
  

    // ======================================================
    // COUPLING SECTION
    // ======================================================
  
    //! initalize PDE coupling
    virtual void InitCoupling(PDECoupling * Coupling);

    //! Fill in input coupling terms
    virtual void CalcInputCoupling();
  
    //! calculate coupling terms
    virtual void CalcOutputCoupling();
  
    //! returns if PDE can compute the quantity
    virtual bool HasOutput(SolutionType output);

    //! do PostProcessing step
    void CalcResults( shared_ptr<BaseResult> result );


  protected:
 
    // Define avilable result types
    void DefineAvailResults();
  
    //! Init the time stepping
    void InitTimeStepping();

  private:

    /// returns the vector of the stokesFluid pressure solution
    ///belonging to all nodes of the actual element
    void GetPresSolVecOfElement(Vector<Double>& sol, StdVector<UInt>& connect_PDE);
    

    //! calculate the vector of coupling forces to the mechanical PDE
    void CalcMechCouplingRHS( StdVector<Elem*> * couplingElems, 
                              StdVector<UInt> & couplingNodes,
                              Vector<Double>& elemCouplingSols,
                              UInt couplingdof );
                              
  

    /// stores an algsys_ vector into a StdVector and returns that L2-norm
    void StoreAlgsysToVec(Vector<Double>& vec, Double * pt);

    /// returns that L2-norm of an algsys vector
    Double AlgsysL2Norm(Double * pt);
  
    /// returns the solution matrix belonging to all nodes of the actual element
    void GetSolOfElement( Matrix<Double>& elDisp, StdVector<UInt>& connect_PDE);

    /// flag for reduced Integration for each subdomain
    StdVector<std::string> reducedIntegration_;

    //! Contains the regions above which the deformed volume is computed
    StdVector<RegionIdType> volAboveDefSurfRegions_;

    //! Contains the directions for which the deformed volume is computed
    StdVector<std::string> volAboveDefSurfDir_;

    // ========================
    // coupling
    // ========================
    //! assigns each coupling element node the according Coupling Node number
    StdVector<StdVector<StdVector<UInt> > > elemNodeToCouplingNode_; 

    //! true, if solution should be written to result file
    bool saveSolVel_;
    bool saveSolPres_;
    bool saveSolVort_;


  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class StokesFluidPDE
  //! 
  //! \purpose 
  //! This class defines the stokes fluid field PDE.
  //! 
  //! \collab 
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 

#endif

} // end of namespace
#endif

