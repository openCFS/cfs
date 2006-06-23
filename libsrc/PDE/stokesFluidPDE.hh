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
      \param aTimeFunc pointer to class TimeFunc
    */
    StokesFluidPDE(Grid *aGrid, TimeFunc *aTimeFunc, WriteResults *aOutFile );

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


    // ======================================================
    // POSTPROC SECTION
    // ======================================================

    //! write results in file
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    virtual void WriteResultsInFile(const UInt kstep = 0,
                                    const Double asteptime = 0.0,
                                    UInt stepOffset = 0,
                                    Double timeOffset = 0.0);

    //! write history results in file
    //! \param stepOffset offset for starting (time)step
    //! \param timeOffset offset for starting time  
    void WriteHistoryInFile(const UInt kstep,
                            const Double asteptime,
                            UInt stepOffset = 0,
                            Double timeOffset = 0.0);

    //! do PostProcessing step
    virtual void PostProcess( );

  protected:
 
    
    //! Obtain information on desired output quantities from parameter file
  
    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    //! The output quantities currently supported by the mechanics PDE are
    //! given in the following table. Here 'Keyword' and 'Result Type' refer
    //! to the XML parameter file, while 'Class Attribute' refers to the
    //! internal attribute of the StokesFluidPDE class that is set, if the keyword
    //! is specified.\n\n
    //! <table border="1">
    //!   <tr>
    //!     <td><b>Keyword</b></td>
    //!     <td><b>Result Type</b></td>
    //!     <td><b>Class Attribute</b></td>
    //!   </tr>
    //!   <tr>
    //!     <td>displacement</td>
    //!     <td>nodeResults</td>
    //!     <td>savesol_</td>
    //!   </tr>
    //!   <tr>
    //!     <td>velocity</td>
    //!     <td>nodeResults</td>
    //!     <td>savederiv_</td>
    //!   </tr>
    //!   <tr>
    //!     <td>acceleration</td>
    //!     <td>nodeResults</td>
    //!     <td>savederiv2_</td>
    //!   </tr>
    //! </table>
    void ReadStoreResults();
  
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

