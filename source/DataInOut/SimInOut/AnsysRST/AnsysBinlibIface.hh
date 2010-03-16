#ifndef ANSYS_BINLIB_IFACE_H
#define ANSYS_BINLIB_IFACE_H

#include <DataInOut/simOutput.hh>

#include "Utils/DynamicObject.hh"

// AnsysBinlibIface class -- this class should be in a shared library, so it
// can be linked with both the program, and dynamically loadable library

// AnsysBinlibIface.h
namespace CoupledField
{

  class AnsysBinlibIface : public DynamicObject
  {
  public:
    AnsysBinlibIface(void (*delObj)(void*));
    virtual ~AnsysBinlibIface(void);

    virtual void GetANSYSRevisionInfos(std::string& BinLibRev,
                                       std::string& MachineId,
                                       std::string& CompiledRev,
                                       const std::string& TestFN) = 0;

    //! Constructor
    virtual void SetFNAndOutputNode( const std::string& fileName,
                                     const std::string& formatName,
                                     const std::string& dirName,
                                     const std::string& pathSep,
                                     PtrParamNode outputNode ) = 0;
    //! Initialize
    virtual void Init( Grid* ptGrid,
                       bool printGridOnly ) = 0;

    //! Begin multisequence step
    virtual void BeginMultiSequenceStep( UInt step,
                                         BasePDE::AnalysisType type,
                                         UInt numSteps) {};

    virtual void SetResultFileName( std::string fn ) = 0;

    
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

  };
  
}

#endif
