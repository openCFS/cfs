#ifndef FILE_NEWBASEMECHPDE
#define FILE_NEWBASEMECHPDE

#include "SinglePDE.hh"
 
namespace CoupledField
{



  //! Class for mechanic equation (no adaptivity)
  class MechPDE: public SinglePDE
  {

  public:

    //!  Constructor. here we read integration parameters
    /*!
      \param aGrid pointer to grid
      \param aOutFile  pointer to class WriteResults. output data.
      \param aTimeFunc pointer to class TimeFunc
    */
    MechPDE(Grid *aGrid, TimeFunc *aTimeFunc, WriteResults *aOutFile );

    //!  Deconstructor
    virtual ~MechPDE();

    //! read in damping information, see SinglePDE.cc  and SinglePDE.hh
    void ReadDampingInformation();

    //! Initialize NonLinearities
    virtual void InitNonLin();

    //! define all (bilinearform) integrators needed for this pde
    virtual void DefineIntegrators( );

    //! define the SoltionStep-Driver
    virtual void DefineSolveStep();

    /// returns a stiffness integrator appropriate to the actual problem (e.g.3D)
    BaseForm * GetStiffIntegrator(MaterialData& actSDMat,
                                  Boolean reducedInt=FALSE);
  

    // ======================================================
    // COUPLING SECTION
    // ======================================================
  
    //! initalize PDE coupling
    virtual void InitCoupling(PDECoupling * Coupling);

    //! calculate coupling terms
    virtual void CalcOutputCoupling();
  
    //! returns if PDE can compute the quantity
    virtual Boolean HasOutput(SolutionType output);

    /// setup source term
    void SetupRHS( );
  

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

  
    //computes mechanical deformation energy
    template <class TYPE>
    void CalcEnergy();

    //computes mechanical stresses
    template <class TYPE>
    void CalcStresses();
    
    //! Obtain information on desired output quantities from parameter file
  
    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    //! The output quantities currently supported by the mechanics PDE are
    //! given in the following table. Here 'Keyword' and 'Result Type' refer
    //! to the XML parameter file, while 'Class Attribute' refers to the
    //! internal attribute of the MechPDE class that is set, if the keyword
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

    // calc rhs coupling to acoustic pde
    // void CalcAcousticCouplingRHS(StdVector<Elem*> * couplingElems,
    // Vector<Double>& forceOnElem);
  
    /// calc rhs coupling to acoustic pde
    void CalcAcousticCouplingRHS( StdVector<Elem*> * couplingElems, 
                                  StdVector<MaterialData*> & materials, 
                                  StdVector<UInt>& couplingNodes,
                                  Vector<Double> & forceOnElem,
                                  UInt couplingdof );
                               
  

    /// does a line search and returns the optimal residual norm
    Double LineSearch(Vector<Double>& solIncrement, Vector<Double>& actSol, 
                      Double& etaLineSearch, Boolean trans=FALSE);


    /// Write nonlin iteration norms to the cla-file
    void WriteClaNlNorms( const UInt iterationCounter,
                          const Double residualL2Norm,
                          const Double extForcesL2Norm, const Double residualErr, 
                          const Double solIncrL2Norm, const Double actSolL2Norm, 
                          const Double incrementalErr );
  

    //! read in the domains with prestressing
    void ReadPreStressing();

    //! read in volume sources
    void DefineRegionLoads();

    /// calculates matrices D^_ and D^__ (see Hughes p. 217) for reduced integration
    void CalcReducedMat(MaterialData& lambdaMat, MaterialData& mueMat,
                        MaterialData& mat);

    /// stores an algsys_ vector into a StdVector and returns that L2-norm
    void StoreAlgsysToVec(Vector<Double>& vec, Double * pt);

    /// returns that L2-norm of an algsys vector
    Double AlgsysL2Norm(Double * pt);
  
    /// flag for reduced Integration for each subdomain
    StdVector<std::string> reducedIntegration_;

    /// returns the solution matrix belonging to all nodes of the actual element
    void GetSolOfElement( Matrix<Double>& elDisp, StdVector<UInt>& connect_PDE);

    //! Number of dimension for stresses
    UInt stressDim_;
    
    //! Flag indicating the use of fractional damping
    Boolean fracDamping_;

    //!
    StdVector<RegionIdType> preStressDomain_;
    StdVector<Double> preStressValX_; //! orientation in x
    StdVector<Double> preStressValY_; //! orientation in y
    StdVector<Double> preStressValZ_; //! orientation in z


    /// value of prestress
    Double preStressVal_;

    /// direction of prestress
    Directions preStressDir_;
  
    /// material data for reduced integration
    MaterialData * lambdaMat;

    /// material data for reduced integration
    MaterialData * mueMat;

    /// external forces (for nonlin simulations)
    Vector<Double> extForces_;

    //@{ \name Attributes related to post-processing

    //! Contains elemen stresses
    BaseElemStoreSol * stress_;

    //! Contains the subdomains on which the stress is computed
    StdVector<RegionIdType> calcStress_;

    //! Contains the subdomains on which the energy is computed
    StdVector<RegionIdType> calcEnergy_;

    //! Contains mechanic velocity
    NodeStoreSol<Double> solDeriv1_;
  
    //! Contains mechanic acceleration
    NodeStoreSol<Double> solDeriv2_;

    //! Contains the regions above which the deformed volume is computed
    StdVector<RegionIdType> volAboveDefSurfRegions_;

    //! Contains the directions for which the deformed volume is computed
    StdVector<std::string> volAboveDefSurfDir_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class MechPDE
  //! 
  //! \purpose 
  //! This class defines the mechanical field PDE and the according 
  //! postprocessing methods.
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

