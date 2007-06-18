// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_DOMAIN_2001
#define FILE_DOMAIN_2001

#include <map>

#include "Utils/StdVector.hh"
#include "Driver/basedriver.hh"
#ifndef INTEGLIB
#include "Utils/mathParser/mathParser.hh"
#endif

namespace CoupledField
{

  //! Forward class declarations
  class BasePDE;
  class BaseForm;
  class StdPDE;
  class SinglePDE;
  class IterCoupledPDE;
  class PDECoupling;
  class DirectCoupledPDE;
  class Grid;
  class FileType;
  class CoordSystem;
  class MaterialHandler;
  class Optimization;
  class DesignSpace;
  class SingleDriver;
  class MultiSequenceDriver;
  class SimInput;
  class ResultHandler;
  class ParamNode;



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
    */
    Domain( std::map<std::string, StdVector<shared_ptr<SimInput> > >& gridInputs,
            ResultHandler * handler,
            MaterialHandler * ptMat );
    
    //! Destructor
    virtual ~Domain();
    
    //! Trigger output of the grid
    void PrintGrid( );
    
    /** dumps intersting information for the developer */
    void Dump();
    
    // ======================================================
    // INIT AND UPDATE ROUTINES
    // ======================================================
    
    //@{
    //! \name Methods for initialization and update

    //! Create Grid and grid related objects (coordinate systems etc.)
    void CreateGrid();

    //! Create PDE objects
    //! \param sequenceStep step index in MultiSequenceSimulation
    void CreatePDEs( UInt sequenceStep );
    
    //! Initialize all PDEs which are previously created
    //! \param sequenceStep step index in MultiSequenceSimulation
    void InitPDEs( UInt sequenceStep );

    //! Delete pointer to PDEs and create them new
    void ResetPDEs();

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

    //! Get pointer to SinglePDE by name
    SinglePDE * GetSinglePDE(const std::string pdename);

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
    { return ptDirectCoupledPde_[0]; };

    //! Get pointer to input-file
      //    FileType * GetInFile(){ return InFile_;}

    //! Get pointer to result handling object
    ResultHandler * GetResultHandler() { return resultHandler_; }

    //! Get pointer to material handler
    MaterialHandler * GetMaterialHandler() {return ptMatHandler_; }

    //! Get pointer to grid object
    Grid * GetGrid( const std::string& id = "default" );

    //! Return local coordinate system by name
    CoordSystem * GetCoordSystem( const std::string & name 
                                  = std::string("default") );

    /** For some optimization (SIMP) or normal simulation with external
     * ersatz material file the forms get the density value by this
     * function. Any penaltisation (SIMP) is applied!
     * @param from eleme elemNum and regionId are used
     * @param form the actual bilinearform, this is used to find the power for SIMP
     * @return -1.0 if nothing available for the elment, otherwise density^power */
    Double GetErsatzMaterial(const Elem* elem, const BaseForm* form);

    /** This is set by optimization which holds the data (in a derved form). It
     * is also reset here by the optimization destructor.
     * @param ersatzMaterial pointer to a data set. NULL to reset, such that ~Domain() doesn't delete it.
     * @param regionId the region for the ersatMaterial */
     void SetErsatzMaterial(DesignSpace* data)
     {
         ersatzMaterial = data;
     }
     
     /** E.g. the MechPDE needs it in CalcResuls() to write pseudo densities. 
      * @return  NULL but an exception if not set and not silent*/
     DesignSpace* GetErsatzMaterial(bool throw_exception = true); 
      

    /** Returns the optimization
     *  @return null if there is none */
    Optimization* GetOptimization() { return optimization_; };
    
    /** Sets the optimization from outside, like the driver */
    void SetOptimization(Optimization* optimization) { this->optimization_ = optimization; };  

#ifndef INTEGLIB
    //! Return Math Parser object for evaluating math expressions
    MathParser * GetMathParser() { return &mathParser_; }
#endif

    /** The post init does more advancec stuff like reading the ersatz material.
     * For this purpose the constructor needs to be finished. 
     * @excpetion checks for error, thefore this is a void method */
    void PostInit();

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
    void CreateSinglePDEs( UInt sequenceStep );
  
    //! Initialize direct coupled pde(s)

    //! Initialize direct coupled pde(s)
    //! \param sequenceStep step index in MultiSequenceSimulation
    void CreateDirectCoupledPDEs( UInt sequenceStep );

    //! Initialize iterative coupled pde

    //! Initialize iterative coupled pde
    //! \param sequenceStep step index in MultiSequenceSimulation
    void CreateIterCoupledPDE( UInt sequenceStep );


    //! Initialize local coordinate systems as read in from the parameter file
    void CreateCoordinateSystems();
    //@}
  
    /** Reads the 'loadErsatzMaterial' node. Check before that it exists and does 
     * not occur concurrently with an optimization */
    void ReadErsatzMaterial(ParamNode* pn); 
  
    // ======================================================
    // DATA SECTION
    // ======================================================
  
    //@{
    //! \name Data about (coupled) PDEs

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

    //! Pointer to coupling objects

    //! Pointer to the object which hold the information for
    //! pairwise iterative coupling of two StdPDEs
    StdVector<PDECoupling*> couplings_; //!<pointer to coupling objects

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

    //! Mapping between name and coordinate sysem pointer
    std::map<std::string, CoordSystem*> coordSys_;

    /** an optinal optimizer */
    Optimization* optimization_;

    /** The ersatz material pointer is set be the domain or it points 
     * to optimization data */    
    DesignSpace* ersatzMaterial;

#ifndef INTEGLIB
    //! Mathematic parser object
    MathParser mathParser_;
#endif

    //! dimension of the problem
    UInt dim_;
  };

}

#endif // FILE_DOMAIN
