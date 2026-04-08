// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_SIMOUTPUT_HH
#define FILE_CFS_SIMOUTPUT_HH

#include <set>

#include <filesystem>
#include <fstream>
namespace fs = std::filesystem;

#include "Domain/Results/BaseResults.hh"
#include "Domain/Results/ResultInfo.hh"
#include "PDE/BasePDE.hh"

namespace CoupledField {

  //! Forward class declaration
  
  

  //! Base class for writing result objects to file
  class SimOutput {

  public:
    
    //! Type for associating a result name the related BaseResult objects
    typedef std::map<std::string, StdVector<shared_ptr<BaseResult> > > ResultMapType;

    //! Define capabilities of writing out certain information
    typedef enum {
      NONE,          /*!< No specific output data*/
      MESH,          /*!< Mesh information present*/
      MESH_RESULTS,  /*!< Results defined spatial entities (nodes, elements)*/
      HISTORY,       /*!< Non spatially resolved results (energy, power etc.)
                          or just specified for some few nodes, elements*/
      USERDATA,      /*!< Information about user (environment, build status for
                          program etc.)*/  
      DATABASE       /*!< Detailed internal data, for re-creating FeFunctions
                          and the internal state of PDEs */
    } Capability;

    //! Constructor
    SimOutput( const std::string& fileName, PtrParamNode outputNode,
               PtrParamNode infoNode, bool isRestart );

    //! Destructor
    virtual ~SimOutput();

    //! Get name of output file format(unv, gid, hdf5)
    const std::string& GetName() { return formatName_; };

    //! Initialize
    virtual void Init( Grid* ptGrid,
                       bool printGridOnly ) = 0;

    //! Get capabilites of interface
    const std::set<Capability>& GetCapabilities() const 
    { return capabilities_; }

    //! Begin multisequence step
    virtual void BeginMultiSequenceStep( UInt step,
                                         BasePDE::AnalysisType type,
                                         UInt numSteps) {};
    
    //! Register result (within one multisequence step)
    virtual void RegisterResult( shared_ptr<BaseResult> sol,
                                 UInt saveBegin, UInt saveInc,
                                 UInt saveEnd,
                                 bool isHistory ) {};

    //! Begin single analysis step
    virtual void BeginStep( UInt stepNum, Double stepVal ) {};

    //! Add result to current step
    virtual void AddResult( shared_ptr<BaseResult> sol ) {};

    //! End single analysis step
    virtual void FinishStep( ) {};

    //! End multisequence step
    virtual void FinishMultiSequenceStep( ) {};

    //! Finalize the output
    virtual void Finalize() {};

    //! If this is a streaming class, default false
    virtual bool IsStreaming() {return false;}

    //! Fill global vector
    template <class TYPE>
    static void FillGlobalVec(Vector<TYPE>& gSol, 
                              const StdVector<shared_ptr<BaseResult> > & solList,
                              ResultInfo::EntityUnknownType entityType );

    /** a public available timer for measuring output times */
    shared_ptr<Timer> timer; 
  protected:
    
    //! Get from complex number the angle in degree with lower limit
    inline Double CPhase( const Complex& c ) const {
      return (std::abs(c.imag()) > 1e-16) ?                   
        std::atan2(c.imag(),c.real() )*180/M_PI : 
         ( c.real() < 0.0 ) ? 180 : 0 ; 
    }
    
    /** checks the result info for nodes and elements. if not such a WARN is printed!
     * @return true if nodes and elements */
    bool ValidateNodesAndElements(ResultInfo& actInfo);

    //! Name of output format
    std::string formatName_;

    //! FileName
    std::string fileName_;

    //! Name of output directory
    fs::path dirName_;

    //! Capabilities of output class
    std::set<Capability> capabilities_;
    
    //! Flag if result file is for a restarted simulation
    bool isRestart_ = false;

    //! Grid class
    Grid* ptGrid_ = nullptr;
    
    //! Parameter node for current output class
    PtrParamNode myParam_;
    
    //! Parameter node for current output class
    PtrParamNode myInfo_;

    //! Current multisequence step in analysis
    UInt actMSStep_ = 0;
    
    //! Current step in analysis (time/frequency step)
    UInt actStep_ = 0;

    //! Current step value ( time / frequency )
    Double actStepVal_ = 0.0;
  };

}

#endif
