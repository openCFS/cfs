#ifndef FILE_STDPDE
#define FILE_STDPDE

#include <fstream>
#include "PDE/BasePDE.hh"

#include <set>

#include "Materials/Composite.hh"
#include "FeBasis/FeFunctions.hh"
#include "FeBasis/FeSpace.hh"
#include "Domain/CoefFunction/CoefFunctionMulti.hh"
#include <Domain/CoefFunction/CoefFunctionMaterialModel.hh>

namespace CoupledField {
  
  
  // forward class declarations
  class BasePairCoupling;
  class WriteResults;
  class StdSolveStep;

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
    
    //! Virtual destructor
    virtual ~StdPDE();
    
    // ======================================================
    // ALGSYS SECTION (SOLVER, ...)
    // ======================================================
    
    //! Define the algebraic system 
    virtual void DefineAlgSys();
    
    //! Create the matrices and Solver as well as Preconditioner
    virtual void CreateMatrices_Solver();
    
    /** call this, e.g. in the case of ProgramOptions::DoEquationMapping() which is cfs -M */
    void CreateEquationMapFile();
    
    // ======================================================
    // GET/SET METHODS
    // ======================================================
    
    //! Returns the feFunction which holds a result related to the specified solutionType
    shared_ptr<BaseFeFunction> GetFeFunction(SolutionType solType);
    
    //! Return all solution FeFunctions
    std::map<SolutionType, shared_ptr<BaseFeFunction> >& GetFeFunctions() { return feFunctions_; }
    
    //! Return all Rhs FeFunctions
    std::map<SolutionType, shared_ptr<BaseFeFunction> >& GetRhsFeFunctions() { return rhsFeFunctions_;  }
    
    //! Return pointer to the SolveStep object
    BaseSolveStep* GetSolveStep();
    
    //! Return vector with resultInfo types
    ResultInfoList& GetResultInfos() { return results_;}
    
    //! Return result for given solutionType
    
    //! Returns the resultInfo related to the specified solutionType
    virtual shared_ptr<ResultInfo> GetResultInfo( SolutionType solType );
    
    shared_ptr<CoefFunctionMaterialModel<Complex>> GetModelCoef(){ return matModelCoef_; };
    
    virtual AnalysisType GetAnalysisType() const {
      return analysistype_;
    }
    
    bool HasComplexMatData(RegionIdType actRegion) {
      return complexMatData_[actRegion]; }
    
    //! returns if PDE can compute the quantity
    virtual bool HasOutput(SolutionType output)
    {
      EXCEPTION("not implemented");
      return false;
    }
    
    virtual bool GetIsaxi() {
      return isaxi_;
    }
    
    //! Set geometry information
    virtual void SetGeomInfo();
    
    //! Set special RHS values
    virtual void SetRhsValues();
    
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
    std::map<RegionIdType, BaseMaterial*>&  GetMaterialData() { return materials_; }
    
    //! Return assemble class, which holds all integrators
    Assemble* GetAssemble() { return assemble_; }
    
    //! Return all regions of the PDE
    StdVector<RegionIdType> GetRegions() { return regions_; }
    
    //! Return pointer to algebraic system
    AlgebraicSys* GetAlgSys() { return algsys_; }

    //! Return pointer to grid the PDE is defined on
    Grid* GetGrid() { return ptGrid_; }
    
    /** Give the damping type by region.
     * @return NONE if no damping in map! */
    DampingType GetDamping(RegionIdType reg_id) const;
    
    //! Set if PDE is nonlinear
    virtual void SetNonLinearity(bool nonLin){ nonLin_=nonLin; };

    // St if PDE is nonlinear (material dependency)
    virtual void SetMaterialNonLinearity(bool nonLin){ nonLinMaterial_=nonLin; };
    //@}

    //@{
    //!  Get functions concerning nonlinearity
    bool IsNonLin()
    { return nonLin_;};

    bool IsNonLinMaterial()
    { return nonLinMaterial_;};

    bool IsTotalNonLinFormulation()
    { return nonLinTotalFormulation_;};

    bool IsHysteresis()
    { return isHysteresis_;};

    bool IsNonLinAndNonHyst()
    { return nonLinNonHyst_;};

    bool IsIterCoupled()
    { return isIterCoupled_;};

    shared_ptr<CoefFunctionMulti> GetHystCoefs(){
      return hysteresisCoefs_;
    }

    void InitHystCoefs(){
      EXCEPTION("Not implemented in base class");
    }

    std::map<RegionIdType, StdVector<NonLinType> >& GetNonLinRegionTypes()
    { return regionNonLinTypes_;};

    //! Return material class
    MaterialClass GetMaterialClass() const { return pdematerialclass_; }
    //@}

    /** Shortcut for the DOF names */
    StdVector<std::string>& GetDofNames(SolutionType st) {
      return feFunctions_[st]->GetResultInfo()->dofNames;
    }

    void TestInversionOfHystOperator(PtrParamNode testNode);

    void EstimateCurrentSlopeForHysteresis(Double steppingLength, Double scaling);

    void CheckSaturationOfHystOperators(Double& lastTSSatAvg, Double& lastItSatAvg, Double& curItSatAvg,
          Double& oppositeDirAsTSAvg, Double& oppositeDirAsItAvg);

    /*
     *
     */
    void SetFlagInCoefFncHyst(std::string flagName, Integer newState);
    void SetDoubleFlagInCoefFncHyst(std::string flagName, Double newState);

    /*
     * when dealing with Hysteresis using StdSolveStep, we need to set/adjust parts of
     * the underlying CoefFunctionHyst
     * However, stdSolveStep cannot directtly access these CoefFunctions but has to go
     * over the PDEs.
     * As basically a PDEs that use Hysteresis will need the same set of functions, it
     * makes sense to directly include them in the base class.
     */
    /*!
     * SetPreviousHystVals -> store input and output values from last iteration
     */
    void SetPreviousHystVals(bool setNextToLastTS = false, bool forceMemoryLock = false);

    bool MaterialTensorsHystDependent();

    virtual void FinalizeAfterTimeStep() {
      EXCEPTION("FinalizeAfterTimeStep has to be implemented for specific PDE");
    }

    virtual void FinilizeBeforTimeStep() {
      ;
    }

    //! Enum for type of nonconforming coupling (Nitsche or Mortar)
    typedef enum {
      NC_NONE,
        NC_MORTAR,
        NC_NITSCHE
    } NcCouplingType;
    static Enum<NcCouplingType> ncCouplingType_;
    static EnumTuple ncTypeTuples_[3];

    //! Enum for type of ansatz functions for Lagrange Multiplier (Mortar only)
    typedef enum {
      LM_STANDARD,
        LM_DUAL_DISCONTINUOUS,
        LM_DUAL_CUBIC
    } LagrangeMultType;
    static Enum<LagrangeMultType> lmType_;
    static EnumTuple lmTypeTuples_[3];

    //! Struct for collecting all formulation-specific options on an NcInterface
    struct NcInterfaceInfo {
      UInt              interfaceId;
      NcCouplingType    type;
      LagrangeMultType  lagrangeMultType;
      Double            nitscheFactor;
      Double            nitscheFactorDamp;
      Double            layerThickness; // parameter for thin layer formulation with non conforming interface condition
      string            layerMaterial; // material for thin layer formulation with non conforming interface condition
      bool              thinLayer;
      bool              crossPointHandling;
      bool              movingMortarForm;
    };
        

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
    
    //! Reads ncInterfaces defined in the XML file
    virtual void ReadNcInterfaces() = 0;
    
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
    
    //! Associate the xml-name of inhom. Dirichlet Bc with SolutionType - 1 timederivative
    std::map<SolutionType, std::string> idbcSolNameMapD1_;
    
    //! Associate the xml-name of inhom. Dirichlet Bc with SolutionType - 2 timederivative
    std::map<SolutionType, std::string> idbcSolNameMapD2_;
    
    //@}
    
    // -----------------------------------------------------------------------
    // Non-linearity
    // -----------------------------------------------------------------------
    
    //@{
    //! \name Attributes connected to nonlinearity
    bool nonLin_;           //!< flag for nonlinear calculations
    bool nonLinMaterial_;           //!< flag for nonlinear material calculations
    bool isAlwaysStatic_;
    bool nonLinTotalFormulation_;   //!< flag for total or incremental NL formulation
    // note: not all regions have to have hysteretic material behavior
    // if there are some regions which are hysteretic and others have other nonlinearities
    // isHysteresis_ and nonLinNonHyst_ are both true
    bool isHysteresis_;
    bool nonLinNonHyst_;
    
    std::string modelName_;

    bool matDepend_;        //!< flag for material dependencies
    
    //! map for each region the type of nonlinearity
    std::map<RegionIdType, StdVector<NonLinType> > regionNonLinTypes_;
    
    //! map for each region the type of material dependency
    std::map<RegionIdType, StdVector<NonLinType> > regionMatDepTypes_;
    
    //! map: for each nonlinearity the id
    std::map<std::string, NonLinType> nonLinTypes_;
    
    //! map: for each material dependency the id
    std::map<std::string, NonLinType> matDepTypes_;
    
    // type of nonlinear algorithm (e.g., Newton)
    NonLinMethodType nonLinMethod_;
    
    // Coefficient function for material model
    shared_ptr<CoefFunctionMaterialModel<Complex>> matModelCoef_;
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
    UInt dim_;                  //!< space dimension of pde
    bool isaxi_;                //!< true: axisymmetric problem
    bool isComplex_;            //!< true, if some part of PDE is complex (Material, solution)
    
    //! list of damping types for all regions
    std::map<RegionIdType,DampingType> dampingList_;
    
    //! use of complex material data per region
    std::map<RegionIdType,bool> complexMatData_;
    
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
    
    //! This map stores the hysteresis coefFunctions
    shared_ptr<CoefFunctionMulti> hysteresisCoefs_;
    
    //! This map stores the irreversible strains obtained by hysteresis coefFunctions
    shared_ptr<CoefFunctionMulti> hysteresisStrain_;
    
  }; // class StdPDE
  
#ifdef DOXYGEN_DETAILED_DOC
  
  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================
  StdPDE
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
