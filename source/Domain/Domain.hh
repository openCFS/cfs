#ifndef DOMAIN_HH
#define DOMAIN_HH

#include <map>
#include "Utils/StdVector.hh"
#include "General/Environment.hh"  
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/SimState.hh"

namespace CoupledField
{

  //! Forward class declarations
  class BasePDE;
  class BaseMaterial;
  class StdPDE;
  class SinglePDE;
  class IterCoupledPDE;
  class DirectCoupledPDE;
  class Grid;
  class FileType;
  class CoordSystem;
  class MaterialHandler;
  class DesignElement;
  class SingleDriver;
  class MultiSequenceDriver;
  class SimInput;
  class ResultHandler;
  
  class SimState;
  class BaseDriver;
  class MathParser;
  class Optimization;
  class DesignSpace;
  struct Elem;
  

  //! This class defines the computational domain.

  //! This class contains information about the geometrical domain,
  //! creates the (coupled) PDE objects, manages the handling of
  //! time functions (used for boundary and load conditions) and
  //! holds the pointers to the input mesh-file and the output data-file.
  class Domain
  {
  public:
    
    //! Constructor. call PostInit() afterwards!
    /*!
      \param aptFileType (input) input file (mesh-data)
      \param isParentDomain If true, this object is considered the primary
                            domain, which will log output to the console.
    */
    Domain( std::map<std::string, StdVector<shared_ptr<SimInput> > >& gridInputs,
            ResultHandler * handler, MaterialHandler * ptMat,
            shared_ptr<SimState> simState, PtrParamNode xmlNode,
            PtrParamNode infoNode,
            bool isParentDomain = true);
    
    //! Destructor
    virtual ~Domain();
    
    //! Trigger output of the grid
    void PrintGrid( );
    
    //! Dumps interesting information for the developer
    void Dump();
    
    //! Return parameter root node of current domain
    PtrParamNode GetParamRoot() {
      return param_;
    }
    
    //! Return info root node of current domain
    PtrParamNode GetInfoRoot() {
      return info_;
    }
    
    // ======================================================
    // INIT AND UPDATE ROUTINES
    // ======================================================
    
    //@{
    //! \name Methods for initialization and update

    //! Create Grid and grid related objects (coordinate systems etc.)
    void CreateGrid();

    //! Create PDE objects
    //! \param sequenceStep step index in MultiSequenceSimulation
    //! \param infoNode infoNode for adding information
    void CreatePDEs(UInt sequenceStep, PtrParamNode infoNode);
    
    /** Alternative to CreatePDEs for optimization with the MultiSequenceDriver
     * Up to now only implemented for single pdes, not for coupled ones
     * @see MultiSequenceDriver::keep_
     * @see MultiSequenceDriver::keptPDEs_ */
    void RestorePDEs(StdVector<SinglePDE*>& pdes);

    //! Initialize all PDEs which are previously created
    //! \param sequenceStep step index in MultiSequenceSimulation
    void InitPDEs( UInt sequenceStep );

    /** reset pdes
     * @param keep if false the pointers are deleted if true they are to be kept in the multi sequence driver */
    void ResetPDEs(bool keep);
    
    //! Set the grids and their IDs from external
    
    //! This methods allows on the pas a map of pre-created
    //! grids. This is useful in a multisequence analysis, where
    //! both sequence steps share the same grid map.  
    void SetGridMap( const std::map<std::string, Grid* >& gridMap );

    //@}

    // ======================================================
    // SET / GETTER METHODS
    // ======================================================

  
    //@{
    //! \name Methods for setting / getting data

    //! Get pointer to basePDE

    //! If only one PDE is defined, this method returns the pointer to it.
    //! In the iterative coupled case, the pointer to the coupled PDE is 
    //! returned
    BasePDE* GetBasePDE();

    //! Get pointer to StdPDE by name
    StdPDE* GetStdPDE(const std::string pdename);

    /** Get pointer to SinglePDE by name.
     * @param throw_exception shall an exception be thrown if the name does not exist
     * @return the pde or NULL if !throw_exception */
    SinglePDE* GetSinglePDE(const std::string pdename, bool throw_exception = true);

    StdVector<SinglePDE*> GetSinglePDEs() const { return ptSinglePde_; }

    //! Get driver object
    BaseDriver* GetDriver();

    /** This extracts the SingleDriver, which is NULL for not initialized 
     * MultiSequenceDrivers(). */ 
    void SetDriver( BaseDriver * driver );

    /** Get driver object. Note, that this might be NULL for not initialized MultiSequenceDriver()! */
    SingleDriver* GetSingleDriver() { return ptSingleDriver_; }

    /** Gets the multi sequence driver or NULL if we have none */
    MultiSequenceDriver* GetMultiSequenceDriver() { return multiSequenceDriver_; }

    //! Get pointer to CoupledPDE
    DirectCoupledPDE* GetDirectCoupledPDE() { return ptDirectCoupledPde_.GetSize() > 0 ? ptDirectCoupledPde_[0] : NULL; }
    
    //! Get map for all registered grids and their reader
    std::map<std::string, Grid* >  GetGridMap() {  return gridMap_;  }

    //! Get pointer to input-file
      //    FileType * GetInFile(){ return InFile_;}

    //! Get pointer to result handling object
    ResultHandler * GetResultHandler() { return resultHandler_; }

    //! Get pointer to material handler
    MaterialHandler * GetMaterialHandler() {return ptMatHandler_; }

    //! Get pointer to simulation state object
    shared_ptr<SimState> GetSimState() {return simState_; }
    
    /** Get pointer to grid object
    * @param id "default" would work */
    Grid* GetGrid(const std::string& id);

    /** works only with a single grid */
    Grid* GetGrid()
    {
      //assert(gridMap_.size() == 1);
      return gridMap_.begin()->second;
    }


    /** return the name of the registered coordinate systems */
    StdVector<std::string> GetCoordSystems() const;

    //! Return local coordinate system by name
    CoordSystem* GetCoordSystem(const std::string& name = std::string("default"));

    //! Return Math Parser object for evaluating math expressions
    MathParser* GetMathParser() { return mathParser_; }

    /** Returns the optimization
     *  @return null if there is none */
    Optimization* GetOptimization() { return optimization_; };

    /** Sets the optimization from outside, like the driver */
    void SetOptimization(Optimization* optimization) { this->optimization_ = optimization; };

    /** E.g. the MechPDE needs it in CalcResuls() to write pseudo densities. */
    DesignSpace* GetDesign() { return designSpace_; }

    bool HasDesign() const { return designSpace_ != NULL; }

    /** This is set by optimization which holds the data (in a derved form) or when we do loadErsatzMaterial or -x
     * Is is also reset here by the optimization destructor.
     * @param ersatzMaterial pointer to a data set. NULL to reset, such that ~Domain() doesn't delete it.
     * @param regionId the region for the ersatz material */
     void SetDesign(DesignSpace* data) { this->designSpace_ = data; }

    /** The post init does more advanced stuff like reading the ersatz material.
     * For this purpose the constructor needs to be finished. 
     * @excpetion checks for error, therefore this is a void method */
    void PostInit( UInt sequenceStep = 1 );

    /** solves the problem, either the "driver" or the optimization problem.
     * Suerly you have to call PostInit() first!*/
    void SolveProblem(); 

    //@}

    /** e.g. coordinate systems
     * Is also used by DensityFile which can enforce the defaults (min_x, max_x, ...) */
    void ToInfo(PtrParamNode info, bool force_default = false);

    /** has any of the single pdes perdiodic boundary conditions set? */
    bool HasPerdiodicBC() const;

<<<<<<< HEAD
    //! Get the generic result index (used in PostProc.cc) to keep track of
    //! generic result definitions
    int GetGenericResultIndex() {
      return genResId_;
    }

    //! Increment the generic result index
    void IncrementGenericResultIndex() {
      genResId_ += 1;
    }


=======
>>>>>>> 31b1c3fad (up-date)
  protected:

  private:
  
    // ======================================================
    // CREATION METHODS
    // ======================================================

    //@{
    //! \name Methods for creating various objects
  
    //! Create driver object
    
    //! Create the SinglePDE objects
    
    //! Create the SinglePDE objects
    //! \param sequenceStep step index in MultiSequenceSimulation
    //! \param infoNode infoNode for adding information
    void CreateSinglePDEs( UInt sequenceStep, PtrParamNode infoNode );

    //! Initialize direct coupled pde(s)

    //! Initialize direct coupled pde(s)
    //! \param sequenceStep step index in MultiSequenceSimulation
    //! \param infoNode infoNode for adding information
    void CreateDirectCoupledPDEs( UInt sequenceStep, PtrParamNode infoNode );

    //! Initialize iterative coupled pde

    //! Initialize iterative coupled pde
    //! \param sequenceStep step index in MultiSequenceSimulation
    //! \param infoNode infoNode for adding information
    void CreateIterCoupledPDE( UInt sequenceStep, PtrParamNode infoNode );

    // LUCA ON
    void CreateBemPDE();
    // bool use_nihu_;
    // LUCA OFF

    //! Initialize local coordinate systems as read in from the parameter file
    void CreateCoordinateSystems();
    
    //! Register variables of in element <variableList/>
    void RegisterVariables();
    
    //! Read in a single grid
    void ReadGrid(const std::string & gridId,
                  const StdVector< shared_ptr<SimInput> > & inputs,
                  bool isAxi, Double depth2d);
    //@}
  
    // ======================================================
    // DATA SECTION
    // ======================================================
  
    //@{
    //! \name Data about (coupled) PDEs

    //! ParamNode of xml file
    PtrParamNode param_;
    
    //! Info node
    PtrParamNode info_;
    
    //! Number of Single PDEs
    UInt numSinglePde_;

    //! Number of DirectCoupled PDEs
    UInt numDirectCoupledPde_;

    //! Number of StdPDEs which can couple iteratively

    //! Holds the number of StdPDEs, which can be coupled iteratively. 
    //! This means that SinglePDEs, which are already coupled directly, are
    //! not counted, since each SinglePDE is either directly OR iteratively
    //! coupled.
    UInt numIterCoupledStdPde_;
  
    //! Pointers to SinglePDEs
    StdVector<SinglePDE*> ptSinglePde_;

    //! Pointers to DirectCoupledPDEs
    StdVector<DirectCoupledPDE*> ptDirectCoupledPde_;

    //! Pointer to iterative coupled PDE
    IterCoupledPDE * ptIterCoupledPde_;

    //! Direct coupling status of each SinglePDE

    //! Flagfield, which indicates, if a SinglePDE is direct coupled
    //! or not
    std::map<SinglePDE*,bool> isDirectCoupled_;

    //@}

    /** Pointer to SingleDriver. Note, that this might be NULL in 
     * the case of MultiSequenceDriver - at least before Init() */
    SingleDriver * ptSingleDriver_;

    /** Only set if we really have a MultiSequcenceDriver */
    MultiSequenceDriver* multiSequenceDriver_;

    //! Map with grid ids and related grid objects
    std::map<std::string, Grid* > gridMap_;

    //! Map with grid ids an related input sources
    std::map<std::string, StdVector<shared_ptr<SimInput> > > gridInputs_;
    
    //! Pointer to object managing results
    ResultHandler * resultHandler_;

    //! Pointer to material handler
    MaterialHandler * ptMatHandler_;
    
    //! Pointer to simulation state object
    shared_ptr<SimState> simState_;

    //! Mapping between name and coordinate system pointer
    std::map<std::string, CoordSystem*> coordSys_;

    /** an optinal optimizer */
    Optimization* optimization_ = NULL;

    /** The ersatz material pointer is set be the domain or it points
     * to optimization data */
    DesignSpace* designSpace_ = NULL;


    //! Mathematic parser object
    MathParser * mathParser_;

    //! dimension of the problem
    UInt dim_;
    
    //! flag if domain uses pre-initialized grid map
    bool useExternalGridMap_;
    
    //! flag if object is main domain and output can be logged to console
    bool isParentDomain_;

    //! Generic result index
    //! For iteratively coupled PDEs we can define arbitrary generic results
    //! Since each PDE can define such results, we have to keep track of
    //! those using actual coefFunctions since they require a valid Enum
    //! Since these Enums can not be created on the fly, we have to keep
    //! track of them here and assign them individually to avoid redefinition */
    int genResId_ = 0;
  };

}

#endif // FILE_DOMAIN
