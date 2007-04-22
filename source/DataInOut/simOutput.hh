// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_CFS_SIMOUTPUT_HH
#define FILE_CFS_SIMOUTPUT_HH


#include "Utils/result.hh"
#include "Domain/resultInfo.hh"

namespace CoupledField {

  //! Forward class declaration
  class ParamNode;

  //! Base class for writing result objects to file
  class SimOutput {

  public:
    
    //! Type for associating a result name the related BaseResult objects
    typedef std::map<std::string, 
                     StdVector<shared_ptr<BaseResult> > > ResultMapType;

    //! Define capabilities of writing out certain information
    typedef enum {NONE, MESH, MESH_RESULTS, HISTORY} Capability;

    //! Constructor
    SimOutput( const std::string& fileName, ParamNode * outputNode );

    //! Destructor
    virtual ~SimOutput();

    //! Get name of output file format(unv, gid, hdf5)
    const std::string& GetName() { return formatName_; };

    //! Initialize
    virtual void Init( Grid* ptGrid) { ptGrid_ = ptGrid; }

    //! Get capabilites of interface
    const std::set<Capability>& GetCapabilities() const 
    { return capabilities_; }

    //! Write grid
    virtual void WriteGrid() = 0;

    //! Begin multisequence step
    virtual void BeginMultiSequenceStep( UInt step,
                                         AnalysisType type,
                                         UInt numSteps) {};
    
    //! Register result (within one multisequence step)
    virtual void RegisterResult( shared_ptr<BaseResult> sol,
                                 UInt saveBegin, UInt saveInc,
                                 UInt saveEnd ) {};

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

  protected:
    
    //! Fill global vector
    template <class TYPE>
    void FillGlobalVec( Vector<TYPE>& gSol, 
                        const StdVector<shared_ptr<BaseResult> > & solList,
                        ResultInfo::EntityUnknownType entityType );

    //! Get from complex number the angle in degree with lower limit
    inline Double CPhase( const Complex& c ) const {
      return (std::abs(c.imag()) > 1e-16) ?                   
        std::atan2(c.imag(),c.real() )*180/PI : 0.0;
    }

    //! Name of output format
    std::string formatName_;

    //! FileName
    std::string fileName_;

    //! Name of output directory
    std::string dirName_;

    //! Capabilities of output class
    std::set<Capability> capabilities_;

    //! Grid class
    Grid* ptGrid_;

    //! Parameter node for current output class
    ParamNode * myParam_;

    //! Current step in analysis (time/frequency step)
    UInt actStep_;

    //! Current step value ( time / frequency )
    Double actStepVal_;

  };

}

#endif
