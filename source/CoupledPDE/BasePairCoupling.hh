// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_BASEPAIRCOUPLING_HH
#define FILE_BASEPAIRCOUPLING_HH

#include <set>

#include "General/Environment.hh"
#include "Utils/StdVector.hh"
#include "MatVec/Vector.hh"
#include "Domain/Results/BaseResults.hh"
#include "PDE/StdPDE.hh"
#include "FeBasis/FeFunctions.hh"

namespace CoupledField
{

  // forward class declaration
  class StdPDE;
  class SinglePDE;
  class Assemble;
  class AlgebraicSys;
  class Grid;
  class BaseMaterial;
  struct Composite;
  class FunctionDescription;
  class BaseBDBInt;
  class BaseFieldFunctor;
  class ResultFunctor;
  class SimState;

  //! Base class for pairwise direct coupling of two pdes
  
  //! This class serves as base class for all pair-wise coupling
  //! objects, like piezo- or mech-acoustic coupling. This class and its
  //! derived classes have to perform similar tasks like a SinglePDE, e.g.
  //! - getting the region for the object
  //! - defining couple integrators
  //! - creating an assemble object
  //!
  //! Objects of this class are instantiated by the class Domain.
  //! Afterwards they are passed to an instance of DirectCoupldPDE, which
  //! can have arbitrarely many of them.
  class BasePairCoupling
  {
  public:

    // public typedefs
    typedef std::set<shared_ptr<ResultInfo> > ResultSet;
    typedef StdVector<shared_ptr<BaseResult> > ResultList;
    typedef std::map<shared_ptr<ResultInfo> , ResultList > ResultMap;
    typedef StdVector<shared_ptr<ResultInfo> > ResultInfoList;


    //! Destructor
    virtual ~BasePairCoupling();

    //! Initialization method
    virtual void Init( UInt sequenceStep);
    
    //! 2nd part of initialization
    virtual void FinalizeInit();

    //! specify the time stepping
    virtual void InitTimeStepping() {;};

    //! Write solutions of postprocessing into unv ... files
    void WriteResultsInFile( const UInt kstep,
                             const Double asteptime );

    //! initialize nonlinearities
    virtual void InitNonLin() {};

    //! Do last caluclation step and cleanup
    virtual void Finalize() {};

    // ======================================================
    // GET / SET METHODS
    // ======================================================

    /** @copydoc BasePDE::GetName() */
    const std::string& GetName() const { return couplingName_; }

    //! Set pointer to algsys
    void SetAlgSys( AlgebraicSys *algSys)
    { algsys_ = algSys;}

    //! Get available results
    ResultSet GetAvailResults() const {
      return availResults_;
    }

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


    //! Ger ParamNode of coupling object
    PtrParamNode GetParamNode() { return myParam_; }

    bool nonLin_;             //!< flag for nonlinear calculations

    void SetNonLinearity(bool nonLin){
      nonLin_=nonLin;};

    //! Pointer to solution strategy object
    shared_ptr<SolStrategy> solStrat_;

    //! Map Storing FeSpaces for each unknown of PDE
    std::map<SolutionType, shared_ptr<BaseFeFunction> > feFunctions_;
    
    //! Map storing the feFunctions of the RHS
    std::map<SolutionType, shared_ptr<BaseFeFunction> > rhsFeFunctions_;

  protected:

    //! Constructor
    BasePairCoupling( SinglePDE *pde1, SinglePDE *pde2,
                      PtrParamNode paramNode, PtrParamNode infoNode,
                      shared_ptr<SimState> simState,
                      Domain* domain);

    //! Definition of the (bi)linear forms
    virtual void DefineIntegrators() = 0;

    //! define the FE-spaces and functions
    void DefineFeFunctions();

    //! Create FeSpaces according to formulation
    //! needs to be implemented in the child classes if a lagrange multiplyer is needed
    virtual void CreateFeSpaces( const std::string&  formulation, PtrParamNode infoNode,
                                 std::map<SolutionType, shared_ptr<FeSpace> >& crSpace ) {};

    //! define primary results (new unknowns, just living on the coupling!)
    virtual void DefinePrimaryResults() {};

    //! Define post-processing results
    //! This method defines the post-processing results, which are computed
    //! mostly using the bilienarform of the main problem.
    //! Thus, this method gets called after the integrators are defined
    virtual void DefinePostProcResults() {};

    //! define all computable results
    virtual void DefineAvailResults() {};

    //! Define a field result
    void DefineFieldResult( PtrCoefFct coef, shared_ptr<ResultInfo> res );
    
    //! Obtain information on desired output quantities from parameter file
    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    void ReadStoreResults();

    //! Read special store results information
    virtual void ReadSpecialResults(){};

    //! read material data
    virtual void ReadMaterialData();
    
    virtual void ReadSensorArrayResults();

    // =====================================================
    // Miscellaneous
    // =====================================================

    
    PtrParamNode infoNode_; // from constructor()

    bool isaxi_;             //!< true: axisymmetric problem

    bool geoUpdate_ = false;        //!< flag for geometric update

    //! Dimension of problem
    UInt dim_;
    
    // ======================================================
    // DATA SECTION
    // ======================================================


    // --- Boundary Conditions ---

    //@{
    //! \name Attributes connected to the handling of boundary conditions

    //! index of current set of boundary conditions. For a multiSequence-simulation
    //! this index determines, which set of boundary conditions is applied.
    UInt sequenceStep_  = 0;
    //@}

    //! set containing the types of possible solutions
    ResultSet availResults_;

    //! map containing the result types and their lists
    ResultMap resultLists_;

    //! map containing the results and their lists for history results
    ResultMap histResultLists_;

    //! Type of current analysis
    BasePDE::AnalysisType analysisType_ = BasePDE::NO_ANALYSIS;

    //! Flag indicating use of complex values
    bool isComplex_;
    
    //! Parameter node of direct coupling section
    PtrParamNode myParam_;

    // -----------------------------------------------------------------------
    // Material data
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes handling info on material data

    //! Maps regions and (simple) materials
    std::map<RegionIdType, BaseMaterial*> materials_;
    
    //! use of complex material data per region
    std::map<RegionIdType,bool> complexMatData_;

    //! Maps regions and composite materials
    std::map<RegionIdType, Composite> compositeMaterials_;

    //! material class
    MaterialClass materialClass_ = NO_CLASS;

    //@}

    
    
    //! Map for storing the primary ADB integrators of the problem

    //! This map stores the primary BDB integrators, which can be used for 
    //! calculating spatial derivatives, fluxes and energy.
    std::map<RegionIdType, BaseBDBInt*> bdbInts_;
    
    // for the case that the coupling is not symmetric (i.e. in case of hysteresis
    // with delta formulation) we need to store the counterpart in an additional map
    std::map<RegionIdType, BaseBDBInt*> bdbIntsCounterpart_;
    
    // flag reporting if the additional map bdbIntsCounterpart_ was used
    // true: map has to be considered
    // false: map does not have to be considered
    bool considerCounterpart_;

    //! Map storing functors for calculating general results
    std::map<SolutionType, shared_ptr<ResultFunctor> > resultFunctors_;

    //! Store field coefficient functions
    std::map<SolutionType, PtrCoefFct > fieldCoefs_;
    
    //! (Volume) regions of coupling object
    StdVector<RegionIdType> subdoms_;

    //! Entity lists for current coupling object
    StdVector<shared_ptr<EntityList> > entityLists_;

    //! ncInterface regions of coupling object
    StdVector<RegionIdType> ncInterfaceIds_;

    //! ncInterface information
    StdVector< StdPDE::NcInterfaceInfo > ncInterfaces_;

    //! Name of coupling
    std::string couplingName_;
    
    //! Pointer for saving the internal simulation state
    shared_ptr<SimState> simState_;

    //! Set with matrix types
    std::set<FEMatrixType> matrixTypes_;

    //! Pointer to first pde
    SinglePDE *pde1_;

    //! Pointer to second pde
    SinglePDE *pde2_;

    //! Pointer to assemble object
    Assemble *assemble_ = NULL;

    //! Pointer to grid object
    Grid * ptGrid_;

    //! Pointer to domain object
    Domain* domain_;
    
    //! Pointer to algebraic system
    AlgebraicSys * algsys_;

    ////! Pointer to equation map of first PDE
    //StdVector<FunctionDescription> feFcts1_;

    ////! Pointer to equation map of second PDE
    //StdVector<FunctionDescription> feFcts2_;

    //! ResultInfos of first PDE
    ResultInfoList results1_;

    //! ResultInfos of second PDE
    ResultInfoList results2_;
    //! Helper struct for interpolating field variables at arbitrary points
    struct FieldAtPoints {

      //! Physical Quantity
      shared_ptr<ResultInfo> resultInfo;
     
      //! Filename where points get written to
      std::string fileName;

      //! Format output file as CSV (comma separated values)
      bool csv;

      //! Delimiter for CSV fields
      char delim;

      //! Pointer to coordinate system
      CoordSystem * coordSys;
      
      //! Vector with elements 
      StdVector<const Elem*> elems;
      
      //! Vector with local points
      StdVector<LocPoint> locPoints;

      //! Flux values for each point
      SingleVector* field;
    };
    
    StdVector<FieldAtPoints> sensors_;
    //@}

  };


} // end of namespace

#endif
