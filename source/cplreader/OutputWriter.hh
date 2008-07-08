#ifndef FILE_OUTPUTWRITER_2008
#define FILE_OUTPUTWRITER_2008

#include <vector>
#include <map>
#include <string>

#include "General/environment.hh"
#include "FlowDataTypes.hh"

namespace CoupledField
{

class OutputWriter {
public:
  OutputWriter();

  virtual ~OutputWriter();

  // Init output writer
  virtual void Init(int argc, char** argv) = 0;

  // Write nodal coordinates
  virtual void WriteNodalCoords(const std::vector<Double>& coords) = 0;

  // Write element connectivity/topology
  virtual void WriteTopology(UInt maxNumElemNodes,
                                const std::vector<FEType> elemTypes,
                                const std::vector<UInt> elems) = 0;

  // Set dimension of output mesh
  void SetDim(UInt dim) { dim_ = dim; };

  // Write regions
  virtual void WriteRegion(std::string regionName,
                              const std::vector<UInt> regionElems,
                              const std::vector<UInt> regionNodes) = 0;

  // Write node groups
  virtual void WriteNodeGroups(const std::vector<std::string> groupNames,
                          const std::vector< std::vector<UInt> > groupNodes) = 0;

  virtual void WriteElemGroups(const std::vector<std::string> groupNames,
                          const std::vector< std::vector<UInt> > groupElems,
                          const std::vector< std::vector<UInt> > groupNodes) = 0;

  virtual void BeginStep(UInt stepNum, Double stepValue);
  virtual void EndStep();

  virtual void WriteFlowData(std::vector<FlowDataType> flowData,
                       std::vector<bool> actParts) = 0;

  virtual void WriteUserData(const std::map<std::string, std::string>& userData) = 0;

protected:
  std::map< std::string, std::vector<UInt> > stepNums_;
  std::map< std::string, std::vector<Double> > stepValues_;

  UInt actStepNum_;
  Double actStepValue_;
  std::vector<std::string> actStepResults_;
  UInt dim_;
};

} // end of namespace
#endif
