#ifndef FILE_DOMAIN_2001
#define FILE_DOMAIN_2001

#include <map>

#include "Utils/StdVector.hh"

namespace CoupledField
{

  //! Forward class declarations
  class TimeFunc;
  class BasePDE;
  class StdPDE;
  class SinglePDE;
  class IterCoupledPDE;
  class PDECoupling;
  class DirectCoupledPDE;
  class Grid;
  class FileType;
  class WriteResults;

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
      \param ptOut (input) output file
      \param aptTimeFunc (input) time function data base
    */
    Domain(FileType * const aptFileType, WriteResults * ptOut, 
           TimeFunc * aptTimeFunc);
    
    //! Destructor
    virtual ~Domain();
    
    //! Trigger output of the grid
    void PrintGrid( );
    
    
    // ======================================================
    // INIT AND UPDATE ROUTINES
    // ======================================================
    
    //@{
    //! \name Methods for initialization and update

    //! Initialize all PDEs
    //! \param pdes vector of pointers to pdes
    //! \param sequenceStep step index in MultiSequenceSimulation
    //! \param tags tags for each PDE 
    void InitPDEs(StdVector<std::string> &pdeNames,
                  UInt sequenceStep,
                  StdVector<std::string> tags);

    //! Delete pointer to PDEs and create them new
    void ResetPDEs();

    //! Update algebraic system and bcs after refinement of the mesh
    void Update( );

    //! Update algebraic system in case of new mesh
    void UpdateAlgSys( );

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

    //! Get pointer to input-file
    FileType * GetInFile(){ return InFile_;}

    //! Get pointer to output-file
    WriteResults * GetOutFile(){ return OutFile_;}

    //! Get pointer to grid object
    Grid * GetGrid(){ return ptgrid_;}

    //! Get time pointer to function
    TimeFunc* GetTimeFncPointer() {return ptTimeFunc_;}

    //@}


  protected:

  private:
  
    // ======================================================
    // CREATION METHODS
    // ======================================================

    //@{
    //! \name Methods for creating (coupled) PDE objects
  
    //! Create the SinglePDE objects
    
    //! Create the SinglePDE objects
    //! \param pdeNames Vector of names of PDEs
    void CreateSinglePDEs(StdVector<std::string> & pdeNames);
  
    //! Initialize direct coupled pde(s)

    //! Initialize direct coupled pde(s)
    //! \param sequenceTag Tag specifying the current coupling section
    void CreateDirectCoupledPDEs(StdVector<std::string> & sequenceTags);

    //! Initialize iterative coupled pde

    //! Initialize iterative coupled pde
    //! \param sequenceTag Tag specifying the current coupling section
    void CreateIterCoupledPDE(StdVector<std::string> & sequenceTags);
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
    std::map<SinglePDE*,Boolean> isDirectCoupled_;

    //@}

    //! Pointer to grid object
    Grid * ptgrid_;

    //! Pointer to object handling time functions
    TimeFunc * ptTimeFunc_;

    //! Pointer to object handling input file (mesh data)
    FileType *InFile_;

    //!  Pointer to object handling output file 
    WriteResults * OutFile_;

  };

}

#endif // FILE_DOMAIN
