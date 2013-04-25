#ifndef DOMAIN_HH
#define DOMAIN_HH

#include <map>

#include "Utils/StdVector.hh"
#include "Utils/mathParser/mathParser.hh"

#include "Driver/BaseDriver.hh"


namespace CoupledField
{

  //! Forward class declarations
  class BasePDE;
  class BaseForm;
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
  class ParamNode;
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
      \param allowOutput If true, not output is logged to stdout
    */
    Domain( std::map<std::string, StdVector<shared_ptr<SimInput> > >& gridInputs,
            ResultHandler * handler, MaterialHandler * ptMat,
            shared_ptr<SimState> simState, PtrParamNode xmlNode,
            PtrParamNode infoNode,
            bool allowOutput = true);
    
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
    void CreatePDEs( UInt sequenceStep, PtrParamNode infoNode );
    
    //! Initialize all PDEs which are previously created
    //! \param sequenceStep step index in MultiSequenceSimulation
    void InitPDEs( UInt sequenceStep );

    //! Delete pointer to PDEs and create them new
    void ResetPDEs();
    
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
    BasePDE * GetBasePDE();

    //! Get pointer to StdPDE by name
    StdPDE * GetStdPDE(const std::string pdename);

    /** Get pointer to SinglePDE by name.
     * @param throw_exception shall an exception be thrown if the name does not exist
     * @return the pde or NULL if !throw_exception */
    SinglePDE * GetSinglePDE(const std::string pdename, bool throw_exception = true);

    //! Get driver object
    BaseDriver* GetDriver();

    /** This extracts the SingleDriver, which is NULL for not initialized 
     * MultiSequenceDrivers(). */ 
    void SetDriver( BaseDriver * driver );

    /** Get driver object. Note, that this might be NULL for not initialized 
     * MultiSequenceDriver()! */
    SingleDriver * GetSingleDriver() { return ptSingleDriver_; }

    //! Get pointer to CoupledPDE
    DirectCoupledPDE* GetDirectCoupledPDE()
    {  if (ptDirectCoupledPde_.GetSize() > 0)
      return ptDirectCoupledPde_[0]; 
    else 
      return NULL;
    }
    
    //! Get map for all registered grids and their reader
    std::map<std::string, Grid* >  GetGridMap() {
      return gridMap_;
    }

    //! Get pointer to input-file
      //    FileType * GetInFile(){ return InFile_;}

    //! Get pointer to result handling object
    ResultHandler * GetResultHandler() { return resultHandler_; }

    //! Get pointer to material handler
    MaterialHandler * GetMaterialHandler() {return ptMatHandler_; }

    //! Get pointer to simulation state object
    shared_ptr<SimState> GetSimState() {return simState_; }
    
    //! Get pointer to grid object
    Grid * GetGrid( const std::string& id = "default" );

    //! Return local coordinate system by name
    CoordSystem * GetCoordSystem( const std::string & name 
                                  = std::string("default") );

    //! Return Math Parser object for evaluating math expressions
    MathParser * GetMathParser() { return mathParser_; }

    /** The post init does more advancec stuff like reading the ersatz material.
     * For this purpose the constructor needs to be finished. 
     * @excpetion checks for error, thefore this is a void method */
    void PostInit( UInt sequenceStep = 0 );

    /** solves the problem, either the "driver" or the optimization problem.
     * Suerly you have to call PostInit() first!*/
    void SolveProblem(); 

    //@}


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


    //! Initialize local coordinate systems as read in from the parameter file
    void CreateCoordinateSystems();
    
    //! Register variables of in element <variableList/>
    void RegisterVariables();
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

    //! Mapping between name and coordinate sysem pointer
    std::map<std::string, CoordSystem*> coordSys_;

    //! Mathematic parser object
    MathParser * mathParser_;

    //! dimension of the problem
    UInt dim_;
    
    //! flag if domain uses pre-initialized grid map
    bool useExternalGridMap_;
    
    //! flag if logging output should be enabled
    bool output_;
  };

}

#endif // FILE_DOMAIN
