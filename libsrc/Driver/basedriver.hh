#ifndef FILE_BASEDRIVER_2001
#define FILE_BASEDRIVER_2001

#include "General/environment.hh"
#include "Utils/tools.hh"

namespace CoupledField
{

  // forward class declarations
  class Domain;
  class WriteResults;

  //! Base class for driving classes where we implemented time-stepping
  class BaseDriver
  {
  public:

    //! Constructor
    //! \param adomain pointer to class Domain
    BaseDriver(Domain * adomain);

    //! Destructor
    virtual ~BaseDriver();
    
    //! Initialization method
    virtual void Init() = 0;

    //! Main method for solvin the problem
    virtual void SolveProblem()=0;

    //! Return current analysistype

    //! Returns the current analysistype. 
    //! \param pdeName Name of the pdename in case there is a coexistence
    //!                of two different analysistypes in one analysisstep
    //!                (e.g. transient-harmonic)
    virtual AnalysisType GetAnalysisType( const std::string& pdename) { 
      return analysis_; }
  
  protected:
    
    //! type of analysis
    AnalysisType analysis_;

    //! pointer to class Domain
    Domain * ptdomain_;

    //! --------------------- stuff for computation with adaptivity
    //! for printing a sequence of files in dir meshes in gmv-format
    WriteResults * ptMeshes_;

    //! counter of meshes for printing meshes
    UInt nummeshes_;  

    //! print mesh in special file. this method is used in adaptive procedure for space.
    void PrintSeqMeshes();

    //! auxiliary function for computation with adaptivity: open files for printing sequence of refined meshes with error map 
    bool printMeshesOrNot();
  
  private:
  
  };

}

#endif // FILE_DRIVER
