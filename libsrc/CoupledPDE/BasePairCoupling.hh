#ifndef FILE_BASEPAIRCOUPLING_HH
#define FILE_BASEPAIRCOUPLING_HH

#include <set>

#include "General/environment.hh"
#include "Utils/StdVector.hh"

namespace CoupledField
{

  // forward class declaration
  class SinglePDE;
  class Assemble;
  class MHassemble;
  class BaseSystem;
  class Grid;
  class MaterialData;

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
    SinglePDE* GetPde1()
    { return pde1_;}

    //! Return pointer to second PDE
    SinglePDE* GetPde2()
    { return pde2_;}
    
    //! Get types of needed matrices (sysmtem, stiffness,..)
    void GetMatrixTypes( std::set<FEMatrixType> &matTypes);

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
    Assemble * assemble_;

    //! Pointer to grid object
    Grid * ptGrid_;

    //! Pointer to algebraic system
    BaseSystem * algsys_;


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
