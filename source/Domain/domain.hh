// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_DOMAIN_2001
#define FILE_DOMAIN_2001

#include <map>

#include "Utils/StdVector.hh"
#ifndef INTEGLIB
#include "Utils/mathParser/mathParser.hh"
#endif

namespace CoupledField
{

  //! Forward class declarations
  class BasePDE;
  class StdPDE;
  class SinglePDE;
  class IterCoupledPDE;
  class PDECoupling;
  class DirectCoupledPDE;
  class Grid;
  class FileType;
  class CoordSystem;
  class MaterialHandler;
  class BaseDriver;
  class SingleDriver;
  class SimInput;
  class ResultHandler;

  //! This class defines the computational domain.

  //! This class contains information about the geometrical domain,
  //! creates the (coupled) PDE objects, manages the handling of
  //! time functions (used for boundary and load conditions) and
  //! holds the pointers to the input mesh-file and the output data-file.
  class Domain
  {
  public:
    
    //! Constructor
    /*!
      \param aptFileType (input) input file (mesh-data)
    */
    Domain(SimInput * const aptMeshIOModule, ResultHandler * handler,
           MaterialHandler * ptMat );
    
    //! Destructor
    virtual ~Domain();
    
    //! Trigger output of the grid
    void PrintGrid( );
    
    
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

    //! Set driver object
    void SetDriver( BaseDriver * driver );

    //! Get driver object
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
    Grid * GetGrid(){ return ptgrid_;}

    //! Return local coordinate system by name
    CoordSystem * GetCoordSystem( const std::string & name 
                                  = std::string("default") );

#ifndef INTEGLIB
    //! Return Math Parser object for evaluating math expressions
    MathParser * GetMathParser() { return &mathParser_; }
#endif

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


    //! Initialize local coordinate systems

    //! Initialize local coordinate systems as read in from the parameter file
    void CreateCoordinateSystems();
    //@}
  
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

    //! Pointer to grid object
    Grid * ptgrid_;

    //! Pointer to SingleDriver
    SingleDriver * ptSingleDriver_;

    //! Pointer to object handling input file (mesh data)
    SimInput *simInput_;

    //! Pointer to object managing results
    ResultHandler * resultHandler_;

    //! Pointer to material handler
    MaterialHandler * ptMatHandler_;

    //! Mapping between name and coordinate sysem pointer
    std::map<std::string, CoordSystem*> coordSys_;

#ifndef INTEGLIB
    //! Mathematic parser object
    MathParser mathParser_;
#endif

    //! dimension of the problem
    UInt dim_;

  };

}

#endif // FILE_DOMAIN
