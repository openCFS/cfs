// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

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
    BaseDriver();

    //! Destructor
    virtual ~BaseDriver();
    
    //! Initialization method
    virtual void Init() = 0;

    //! Main method for solvin the problem
    virtual void SolveProblem()=0;

    //! Return current analysistype

    //! Returns the current analysistype. 
    virtual AnalysisType GetAnalysisType( ) { 
      return analysis_; }
    
    //! Return current (multi)sequenceStep
    virtual UInt GetActSequenceStep();
    
    //! Return current time / frequency step of simulation
    virtual UInt GetActStep ( const std::string& pdename ) = 0;
  
  protected:
    
    //! type of analysis
    AnalysisType analysis_;

    //! --------------------- stuff for computation with adaptivity
    //! for printing a sequence of files in dir meshes in gmv-format

    //! counter of meshes for printing meshes
    UInt nummeshes_;  

    //! current analysis step in a multiSequence analysis
    UInt actSequenceStep_;

    //! print mesh in special file. this method is used in adaptive procedure for space.
    void PrintSeqMeshes();

    //! auxiliary function for computation with adaptivity: open files for printing sequence of refined meshes with error map 
    bool printMeshesOrNot();
  
  private:
  
  };

}

#endif // FILE_DRIVER
