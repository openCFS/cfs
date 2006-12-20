#ifndef FILE_BASEPAIRCOUPLING_HH
#define FILE_BASEPAIRCOUPLING_HH

#include <set>

#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include "Utils/vector.hh"
#include "PDE/eqnMap.hh"

namespace CoupledField
{

  // forward class declaration
  class SinglePDE;
  class Assemble;
  class MHassemble;
  class BaseSystem;
  class Grid;
  class BaseMaterial;
  class BaseNodeStoreSol;
  class Composite;
  
  //! Base class for pairwise direct coupling of two pdes
  class BasePairCoupling
  {
  public:

    // public typedefs
    typedef StdVector<shared_ptr<ResultDof> > ResultList;

    //! Destructor
    virtual ~BasePairCoupling();
    
    //! Initialization method
    virtual void Init(UInt sequenceStep = 0,
                      std::string  bcSequenceTag = "anyTag");

    //! Trigger calculation of postprocessing results
    virtual void PostProcess() = 0;

    //! Write solutions of postprocessing into gmv/unv ... files
    virtual void WriteResultsInFile(const UInt kstep,
                             const Double asteptime,
                             UInt stepOffset,
                             Double timeOffse) = 0;

    //! initialize nonlinearities
    virtual void InitNonLin() = 0;
    
    // ======================================================
    // GET / SET METHODS
    // ======================================================
    
    //! Return coupling name
    std::string GetName() 
    { return couplingName_; }

    //! Set pointer to algsys
    void SetAlgSys( BaseSystem *algSys)
    { algsys_ = algSys;}

    //! Set pointer to assemble class
    void SetAssemble( Assemble* assemble ) 
    { assemble_ = assemble; }

    //! Return pointer to first PDE
    SinglePDE* GetPde1() {
      return pde1_;
    }

    //! Return pointer to second PDE

    SinglePDE* GetPde2()
    { return pde2_;}

    //! Get material data object
    std::map<RegionIdType, BaseMaterial*>  getPDEMaterialData()
    {return materials_;};


    //! Return identifier of first PDE
    PdeIdType GetPdeId1();

    //! Return identifier of second PDE
    PdeIdType GetPdeId2();

    bool nonLin_;           //!< flag for nonlinear calculations
    bool nonLinMaterial_;           //!< flag for nonlinear material calculations
    
    void SetNonLinearity(bool nonLin){
      nonLin_=nonLin;};
    
    void SetMaterialNonLinearity(bool nonLin){
      nonLinMaterial_=nonLin;};


  protected:

    //! Constructor
    BasePairCoupling( SinglePDE *pde1, SinglePDE *pde2 );

    //! Definition of the (bi)linear forms
    virtual void DefineIntegrators() = 0;

    //! Obtain information on desired output quantities from parameter file
    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    virtual void ReadStoreResults() = 0;

    //! read material data
    virtual void ReadMaterialData();

    // =====================================================
    // Miscellaneous
    // =====================================================

    BaseNodeStoreSol * sol_;    //!< solution

    bool isaxi_;             //!< true: axisymmetric problem

    bool geoUpdate_;        //!< flag for geometric update


    // -----------------------------------------------------------------------
    // NON_LINEARITY
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes connected to nonlinearity
    //    bool nonLin_;           //!< flag for nonlinear calculations
    Double incStopCrit_;       //!< stopping criterion for incremental error
    Double residualStopCrit_;  //!< stopping criterion for residual error
    UInt nonLinMaxIter_;    //!< maximal number of NL-iterations
    std::string nonLinMethod_; //!< method for handling the non-linearity
    bool nonLinLogging_;    //!< log progress of non-linear iterations
    bool isHysteresis_;     //!< flag for hysteresis

    std::string lineSearch_;   //!< switch for lineSearch
    StdVector<NonLinPDE> nonLinPDEName_;//!< some PDEs carry a name (->acoustics!)
    //@}

    // ======================================================
    // DATA SECTION
    // ======================================================

    // --- Boundary Conditions ---
    
    //@{
    //! \name Attributes connected to the handling of boundary conditions

    //! tag of current set of boundary conditions

    //! tag of current set of boundary conditions. For a multiSequence-simulation
    //! this id determines, which set of boundary conditions is applied.
    std::string bcSequenceTag_;

    //! index of current set of boundary conditions. For a multiSequence-simulation
    //! this index determines, which set of boundary conditions is applied.
    UInt bcSequenceIndex_;
    //@}

    //! Type of current analysis
    AnalysisType analysisType_;

    //! true, if solution should be written to result file
    bool saveSol_;

    //! true, if first derivative of solution should be written to result file
    bool saveDeriv1_;

    //! true, if second derivative of solution should be written to result file
    bool saveDeriv2_;

    //! contains element results of complex valued charge 
    Vector<Complex> complexValuedCharge_;

    //! contains element results of complex valued Efield 
    Vector<Complex> complexValuedEfield_;

    CFSVector * solVec_;        //! needed in iterative coupled computation 

    //! vector containing solutiontypes of PDE
    StdVector<SolutionType> solTypes_;

    //! vector containgin dofs of solutiontypes
    StdVector<UInt> solDofs_;

    // -----------------------------------------------------------------------
    // Material data
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes handling info on material data

    //! Maps regions and (simple) materials
    std::map<RegionIdType, BaseMaterial*> materials_;    

    //! Maps regions and composite materials
    std::map<RegionIdType, Composite> compositeMaterials_;

    //! material class
    MaterialClass materialClass_;
  
    //@}
    
    //! (Volume) regions of coupling object
    StdVector<RegionIdType> subdoms_;

    //! Surface regions of coupling object
    StdVector<RegionIdType> surfRegions_;

    //! Name of coupling
    std::string couplingName_;

    //! Set with matrix types
    std::set<FEMatrixType> matrixTypes_;

    //! Pointer to first pde
    SinglePDE *pde1_;
    
    //! Pointer to second pde
    SinglePDE *pde2_;
    
    //! Pointer to assemble object
    Assemble *assemble_;

    //! Pointer to grid object
    Grid * ptGrid_;

    //! Pointer to algebraic system
    BaseSystem * algsys_;

    //! Pointer to equation map of first PDE
    shared_ptr<EqnMap> eqnMap1_;

    //! Pointer to equation map of second PDE
    shared_ptr<EqnMap> eqnMap2_;

    //! Results of first PDE
    ResultList results1_;

    //! Results of second PDE
    ResultList results2_;

    // -----------------------------------------------------------------------
    // Geometry & node numbering
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes related to geometry and node numbering
    UInt dofspernode_;  //!< number of unknowns per node
    UInt numPDENodes_;  //!< number of nodes in subdomains
    UInt numElems_;     //!< number of elements in subdomains
    UInt dim_;                  //!< space dimension of pde
  
    //! defines subtype of mechanic PDE: plainStrain, 3d, ...
    std::string subType_;
    //@}


  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class BasePairCoupling
  //! 
  //! \purpose This class serves as base class for all pair-wise coupling
  //! objects, like piezo- or mech-acoustic coupling. This class and its
  //! derived classes have to perform similar tasks like a SinglePDE, e.g.
  //! - getting the region for the object
  //! - defining couple integrators
  //! - creating an assemble object
  //! 
  //! \collab Objects of this class are instantiated by the class Domain.
  //! Afterwards they are passed to an instance of DirectCoupldPDE, which 
  //! can have arbitrarely many of them.
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused
  //! 
  //! \improve 
  //! - There should be only one assemble object for all direct-coupled 
  //!   SinglePDEs, so that this class won't have to create an own
  //!

#endif

} // end of namespace

#endif
