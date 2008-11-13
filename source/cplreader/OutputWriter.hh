// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_OUTPUTWRITER_2008
#define FILE_OUTPUTWRITER_2008

#include <vector>
#include <map>
#include <string>

#include "General/environment.hh"
#include "FlowDataTypes.hh"

namespace CoupledField
{
  class CouplingHandler;
  
  class OutputWriter {
  public:
    OutputWriter();
    
    virtual ~OutputWriter();
    
    // Init output writer
    virtual void Init(int argc, char** argv,
                      const std::string& outputPath) = 0;

    // Write nodal coordinates
    virtual void WriteNodalCoords(const std::vector<Double>& coords) = 0;

    // Write element connectivity/topology
    virtual void WriteTopology(UInt maxNumElemNodes,
                               const std::vector<UInt>& elemTypes,
                               const std::vector<UInt>& elems) = 0;

    // Write regions
    virtual void WriteRegion(std::string regionName,
                             const std::vector<UInt>& regionElems,
                             const std::vector<UInt>& regionNodes,
                             UInt regionDim) = 0;

    // Write node groups
    virtual void WriteNodeGroups(const std::map< std::string, std::vector<UInt> >&
                                 nodeGroups) = 0;

    virtual void WriteElemGroups(const std::map< std::string, std::vector<UInt> >&
                                 elemGroups,
                                 const std::map< std::string, std::vector<UInt> >&
                                 groupNodes,
                                 const std::map< std::string, UInt >&
                                 groupDims) = 0;

    virtual void WriteFlowData(CouplingHandler* cplHandler,
                               UInt actRegion,
                               const std::vector<std::string>& outputFields) = 0;

    virtual void WriteUserData(const std::map<std::string, std::string>& userData) = 0;
    
    // Set dimension of output mesh
    void SetDim(UInt dim) { dim_ = dim; };

    virtual void BeginResults(const std::vector<FlowDataType>* flowData) {};
    
    virtual void BeginStep(UInt stepNum, Double stepValue);
    virtual void EndStep();

    virtual void EndResults() {};
    
  protected:
    std::map< std::string, std::vector<UInt> > stepNums_;
    std::map< std::string, std::vector<Double> > stepValues_;
    
    UInt actStepNum_;
    Double actStepValue_;
    std::vector<std::string> actStepResults_;
    UInt dim_;
    const std::vector<FlowDataType>* flowData_;
    bool doneWriting_;
  };

  typedef std::vector< shared_ptr<OutputWriter> >
  OutputWriterVectorType;
  
} // end of namespace
#endif
