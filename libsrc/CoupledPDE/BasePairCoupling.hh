#ifndef FILE_BASEPAIRCOUPLING_HH
#define FILE_BASEPAIRCOUPLING_HH

#include <set>

#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include "PDE/nodeEQN.hh"

namespace CoupledField
{

  // forward class declaration
  class SinglePDE;
  class Assemble;
  class MHassemble;
  class BaseSystem;
  class Grid;
  class MaterialData;
  class BaseNodeStoreSol;

  //! Base class for pairwise direct coupling of two pdes
  class BasePairCoupling
  {
  public:
    
    //! Destructor
    virtual ~BasePairCoupling();
    
    //! Initialization method
    virtual void Init(UInt sequenceStep = 0,
                      std::string  bcSequenceTag = "anyTag");

    //! Trigger calculation of postprocessing results
    virtual void PostProcess() = 0;
    
    // ======================================================
    // GET / SET METHODS
    // ======================================================
    
    //! Return coupling name
    std::string GetName() 
    { return couplingName_; }

    //! Set pointer to algsys
    void SetAlgSys( BaseSystem *algSys)
    { algsys_ = algSys;}

    //! Return pointer to first PDE
    SinglePDE* GetPde1() {
      return pde1_;
    }

    //! Return pointer to second PDE

    SinglePDE* GetPde2()
    { return pde2_;}

    //! computes the coordinates of an element including the delta
    //! \param connect (input) global node numbers of element
    //! \param ptCoord (output) coordinates of the element nodes
    //!                (nrNodes \f$\times\f$ spaceDim);
    virtual void GetElemCoords(const StdVector< UInt > connect,
                               Matrix< Double > &coordMat );


    //! Return identifier of first PDE
    PdeIdType GetPdeId1();

    //! Return identifier of second PDE
    PdeIdType GetPdeId2();

    
    // ======================================================
    // METHODS FOR ASSEMBLING
    // ======================================================

    //! constructes the matrix graph by providing to the algebraic system the element connectivities
    void SetupMatrixGraph();

    //! specify type of system matrix for AlgebraicSystem
    void AssembleMatrices();
    
    //! setup source term
    void AssembleSrcRHS(const Double time = 0.0);
    
    //!  assemble a nonlinear RHS part
    void AssembleNLRHS(const Double time = 0.0);

    //!  assemble a spring into the system matrix
    void AssembleSprings(const Double time = 0.0);

    //! sets the actual frequency (just needed for harmonic analysis)
    void SetFrequency(Double actFreq);
    
    //! trigger the reassembling of the matrices
    void SetReassemble();

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

    Integer isaxi_;             //!< TRUE: axisymmetric problem

    Matrix<Double> deltCoords_;    //!< offset to grid coordinates

    Boolean geoUpdate_;        //!< flag for geometric update

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

    //! TRUE, if solution should be written to result file
    Boolean saveSol_;

    //! TRUE, if first derivative of solution should be written to result file
    Boolean saveDeriv1_;

    //! TRUE, if second derivative of solution should be written to result file
    Boolean saveDeriv2_;

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

    //! material data structure
    MaterialData *materialData_;     
    
    //! material class
    std::string materialClass_;
  
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

    NodeEQN * eqnData_;         //!< equation handling

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
