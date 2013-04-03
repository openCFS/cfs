#ifndef FILE_STDPDE
#define FILE_STDPDE

#include <fstream>
#include "PDE/BasePDE.hh"

#include <set>

#include "Materials/Composite.hh"
#include "FeBasis/FeFunctions.hh"
#include "FeBasis/FeSpace.hh"

namespace CoupledField {


  // forward class declarations
  class BasePairCoupling;
  class WriteResults;
  class BaseNodeStoreSol;
  class StdSolveStep;
  class ParamNode;
  class BiotSavart;
  class BaseFeFunction;
  
  //! Base class for all single-field and direct-coupled problems

  class StdPDE : public BasePDE
  {
  
  public:

    // public typedefs
    typedef StdVector<shared_ptr<ResultInfo> > ResultInfoList;
    //! typedefs for result handling
    typedef std::set<shared_ptr<ResultInfo> > ResultSet;

    //! Typedef for nonconforming stuff Nietsche or Mortar
    typedef enum{
      NONE,
      MORTAR,
      NITSCHE
    }NcCouplingType;
    static Enum<NcCouplingType> ncCouplingType_;

    //! Virtual destructor
    virtual ~StdPDE();
    
    // ======================================================
    // ALGSYS SECTION (SOLVER, ...)
    // ======================================================
  
    //! Define the algebraic system 
    virtual void DefineAlgSys();

    //! Create the matrices and Solver as well as Preconditioner
    virtual void CreateMatrices_Solver();
    

    // ======================================================
    // GET/SET METHODS
    // ======================================================

    //! Returns the feFunction which holds a result related to the specified solutionType
    virtual shared_ptr<BaseFeFunction> GetFeFunction( SolutionType solType);
    
    //! Return all solution FeFunctions
    virtual  
    std::map<SolutionType, shared_ptr<BaseFeFunction> > GetFeFunctions( ) {
      return feFunctions_;
    }
    
    //! Return all Rhs FeFunctions
    virtual 
    std::map<SolutionType, shared_ptr<BaseFeFunction> > GetRhsFeFunctions() {
      return rhsFeFunctions_;
    }
    
    //! Return pointer to the SolveStep object
    BaseSolveStep * GetSolveStep();
    
    //! Return vector with resultInfo types
    ResultInfoList& GetResultInfos() { return results_;}

    //! Return result for given solutionType

    //! Returns the resultInfo related to the specified solutionType
    virtual shared_ptr<ResultInfo> GetResultInfo( SolutionType solType );


    virtual AnalysisType GetAnalysisType() {
      return analysistype_;
    }
  
    //! returns if PDE can compute the quantity
    virtual bool HasOutput(SolutionType output)
    {
      EXCEPTION("not implemented");
      return false;
    }
  
    virtual bool GetIsaxi() {
      return isaxi_;
    }
    
    //! Pass boundary conditions to the algebraic system 
    virtual void SetBCs() = 0;
    
    //! Return pointer to paramNode of current pde
    PtrParamNode GetParamNode() { return myParam_; }
    
    //! Return pointer to infoNode of curent pde
    PtrParamNode GetInfoNode() { return myInfo_; }

    //! Init the time stepping
    virtual void InitTimeStepping()
    {EXCEPTION("InitTimeStepping not implemented");};

    //! Get couplings object
    virtual StdVector<BasePairCoupling*>* GetCouplingsObject() 
    {
      return NULL;
    };

    // ======================================================
    // COMMUNICATION ROUTINES FOR PARAMETER IDENTIFICATION
    // ======================================================

    //@{

    //! Return list with material definition for each region
    std::map<RegionIdType, BaseMaterial*>  GetMaterialData()
    {return materials_;};
    
    //! Return assemble class, which olds all integrators
    Assemble * GetAssemble(){return assemble_;}

    //! Return all regions of the PDE
    StdVector<RegionIdType> GetRegions() {
      return regions_;}

    //! Return pointer to algebraic system
    AlgebraicSys * GetAlgSys(){return algsys_;}
    
    //! Return pointer to grid the PDE is defined on
    Grid * GetGrid() {return ptGrid_;}

    //! Set if PDE is nonlinear
    virtual void SetNonLinearity(bool nonLin){
      nonLin_=nonLin;};

    // St if PDE is nonlinear (material dependency)
    virtual void SetMaterialNonLinearity(bool nonLin){
      nonLinMaterial_=nonLin;};
    //@}

    //@{
    //!  Get functions concerning nonlinearity
    bool IsNonLin() 
    { return nonLin_;};

    bool IsNonLinMaterial() 
    { return nonLinMaterial_;};

    bool IsHysteresis() 
    { return isHysteresis_;};

    bool IsIterCoupled() 
    { return isIterCoupled_;};

    std::map<RegionIdType, StdVector<NonLinType> >& GetNonLinRegionTypes() 
    { return regionNonLinTypes_;};

    //! Return material class
    MaterialClass GetMaterialClass() const { return pdematerialclass_; }
    //@}
      
  protected:
  
    //! Constructor
    /*!
      \param aptgrid pointer to grid
    */
    StdPDE(Grid *aptgrid, PtrParamNode paramNode, PtrParamNode infoNode, 
           shared_ptr<SimState> simState,  Domain* domain);
  
    //! private copy constructor
    StdPDE & operator= (const StdPDE & myPDE) {
      EXCEPTION("Not implemented");

      // The following line is only to satisfy the compiler
      return *this;
    };

    //! retruns, if PDE needs to store previous solution or not
    virtual bool NeedsPrevSol() {
      return  needSolPrev_;
    }

    // ======================================================
    // DATA SECTION
    // ======================================================
  
    // -----------------------------------------------------------------------
    // Geometry & node numbering
    // -----------------------------------------------------------------------
    //@{
    //! \name Attributes related to geometry and node numbering

    //! Pointer to grid object
    Grid * ptGrid_;

    //! Vector containing all regions the PDE is defined on
    StdVector<RegionIdType> regions_;

    //@}

    //!  This is our pde info node. To be set/overwritten in each PDE! 
    PtrParamNode infoNode_; 
    
    //! Parameter node for OLAS
    PtrParamNode olasNode_;

    // -----------------------------------------------------------------------
    // Geometry & node numbering
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes related to geometry and node numbering 
    
    //! Holds the PDE-specific subType of the PDE (e.g. planeStrain for mech)
    std::string subType_;
    
    //@}

    // -----------------------------------------------------------------------
    // Boundary conditions
    // -----------------------------------------------------------------------    

    //@{
    //! \name boundary conditions

    //! Associate the xml-name of hom. Dirichlet Bc with the SolutionType
    std::map<SolutionType, std::string> hdbcSolNameMap_;
    
    //! Associate the xml-name of inhom. Dirichlet Bc with SolutionType
    std::map<SolutionType, std::string> idbcSolNameMap_;
        
    //! Number of additional in. Dirichlet boundary equations due to coupling
    UInt numCouplingBcs_;

    //@}

    // -----------------------------------------------------------------------
    // Non-linearity
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes connected to nonlinearity
    bool nonLin_;           //!< flag for nonlinear calculations
    bool nonLinMaterial_;           //!< flag for nonlinear material calculations
    bool isHysteresis_;     //!< flag for hysteresis

    //! map for each region the type of nonlinearity
    std::map<RegionIdType, StdVector<NonLinType> > regionNonLinTypes_;

    //! map for each nonlinearity the id
    std::map<std::string, NonLinType> nonLinTypes_;

    // type of nonlinear algorithm (e.g., Newton)
    NonLinMethodType nonLinMethod_;

    //@}


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
    MaterialClass pdematerialclass_;   

    //@}


    // -----------------------------------------------------------------------
    // PDE coupling
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes connected to handling PDE coupling
    bool isIterCoupled_;        //!< PDE couples with others
    Vector<Double> matParam_;      //!< change to material parameter
    bool updateCouplingBCs_ ;  //!< flag if coupling BC were already set
  
    //! nodes at which coupling terms are calculated
    std::list<UInt> couplingNodes;    
  
    //! elements at which coupling terms are calculated
    StdVector<Elem*> couplingElements;
  
    //! iteration counter for coupled PDE solution process
    UInt iterCoupledCounter_;
    //@}


    // -----------------------------------------------------------------------
    // Time stepping
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes connected to time stepping
    bool diagMass_;           //!< use of diagonal mass matrix in explicit time stepping
    //@}



    // -----------------------------------------------------------------------
    // Miscellaneous
    // -----------------------------------------------------------------------

    //@{
    //! \name Miscellaneous attributes
    //! Vector containing the results calculated by this PDE
    //! OBSOLETE
    ResultInfoList results_;
    
    //! Set containing the types of possible results
    ResultSet availResults_;
    //!

    //! flag indicating if this PDE needs the algebraic system
    bool needsAlgsys_;

    AnalysisType analysistype_; //!< analysis type
    bool isAlwaysStatic_;    //!< flag for static PDEs (like electrostatic)
    UInt dim_;                  //!< space dimension of pde
    bool isaxi_;             //!< true: axisymmetric problem
    bool isComplex_;         //!< true, if some part of PDE is complex (Material, solution)
    bool needSolPrev_;          //! true, if solution at time step n has also to bve stored

    //! list of damping types for all regions
    std::map<RegionIdType,DampingType> dampingList_;
    
    //! use of complex material data per region
    std::map<RegionIdType,bool> complexMatData_;

    //! checks, if we have for the coupling a incremental solution
    bool isIncrFormulation_;    
    
    //! if yes, PDE is computed on deformed geometry
    bool updatedLagrangeForm_;

    //! Pointer to object of analysis (Static, Trans, Harm or Eig)
    Assemble * assemble_;

    //! pointer to SolveStep classes
    StdSolveStep * solveStep_;
    
    AlgebraicSys * algsys_;      //!< pointer to algebraic system
    
    //! Pointer to solution strategy object
    shared_ptr<SolStrategy> solStrat_;
  
    /** This is the node for linear system responsible for this pde. */
    PtrParamNode olasInfo_;
    
    //@}

    //! Map Storing FeSpaces for each unknown of PDE
    std::map<SolutionType, shared_ptr<BaseFeFunction> > feFunctions_;
    
    //! Map storing the feFunctions of the previous step
    std::map<SolutionType, shared_ptr<BaseFeFunction> > prevFeFunctions_;

    //! Map storing time derivatives of FeFunctions
    std::map<SolutionType, shared_ptr<BaseFeFunction> > timeDerivFeFunctions_;
    
    //! Map Storing the time derivative order of the specific result
    std::map<SolutionType, UInt> timeDerivOrder_;
    
    //! Map associating the primary time derivative results with the primary one
    std::map<SolutionType, SolutionType> timeDerivPrimaryResults_;
    
    //! Map storing the feFunctions of the RHS
    std::map<SolutionType, shared_ptr<BaseFeFunction> > rhsFeFunctions_;

    //! vector containing regionIds of non-conforming interfaces
    StdVector<RegionIdType> ncIFaces_;

    //! map storing for each ncIface the nitsche NMGformulation factor
    std::map<RegionIdType,Double> nitscheFactors_;

    //! Type of non-matching formulation
    std::map<RegionIdType,NcCouplingType> ncTypes_;

  }; // class StdPDE

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class StdPDE
  //! 
  //! \purpose 
  //! This class serves as base class for all single-field and 
  //! direct coupled problems. The idea is, that an iterative coupled
  //! pde can have one or more objects of StdPDE, which are then
  //! solved iteratively.
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
