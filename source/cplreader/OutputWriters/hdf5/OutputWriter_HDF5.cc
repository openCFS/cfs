#include "OutputWriter_HDF5.hh"

namespace CoupledField
{

  OutputWriter_HDF5::OutputWriter_HDF5() {

  }

  OutputWriter_HDF5::~OutputWriter_HDF5() {

  }

  void OutputWriter_HDF5::Init(int argc, char** argv) {

  }

  void OutputWriter_HDF5::WriteNodalCoords(const std::vector<Double>& coords) {

  }

  void OutputWriter_HDF5::WriteTopology(UInt maxNumElemNodes,
                                        const std::vector<FEType> elemTypes,
                                        const std::vector<UInt> elems) {

  }

  void OutputWriter_HDF5::WriteRegion(std::string regionName,
                                      const std::vector<UInt> regionElems,
                                      const std::vector<UInt> regionNodes) {

  }

  void OutputWriter_HDF5::WriteNodeGroups(const std::vector<std::string> groupNames,
                                          const std::vector< std::vector<UInt> > groupNodes) {

  }

  void OutputWriter_HDF5::WriteElemGroups(const std::vector<std::string> groupNames,
                                          const std::vector< std::vector<UInt> > groupElems,
                                          const std::vector< std::vector<UInt> > groupNodes) {

  }

  void OutputWriter_HDF5::WriteFlowData(std::vector<FlowDataType> flowData,
                                        std::vector<bool> actParts) {

  }

  void OutputWriter_HDF5::WriteUserData(const std::map<std::string, std::string>& userData) {

  }

  void OutputWriter_HDF5::CheckOpenObjects() {

  }

  void OutputWriter_HDF5::WriteFileInfoHeader() {

  }

  void OutputWriter_HDF5::InitHDF5() {

  }

  void OutputWriter_HDF5::WriteRegion(const H5::Group& meshGroup,
                                      const std::vector<UInt>& nodes,
                                      const std::vector<UInt>& elems,
                                      const UInt dim,
                                      const std::string name) {

  }

  void OutputWriter_HDF5::InitResultsGroup() {

  }

  void OutputWriter_HDF5::WriteResultDescriptions(UInt numSteps,
                                                  const std::vector<FlowDataType>& outputFields,
                                                  const std::vector<UInt> stepNumbers,
                                                  const std::vector<Double> stepValues) {

  }

  void OutputWriter_HDF5::WriteResults( H5::Group& resultGroup,
                                        std::vector<Double>& resultVals,
                                        const UInt numDOFs,
                                        const bool isImag ) {
  }

  void OutputWriter_HDF5::CreateExternalFile(UInt timeStep) {

  }

  void OutputWriter_HDF5::WriteStringToUserData(const std::string& dSetName,
                                                const std::string& str) {

  }

  void OutputWriter_HDF5::WriteNodeGroups(const H5::Group& meshGroup) {

  }

  void OutputWriter_HDF5::WriteElemGroups(const H5::Group& meshGroup) {

  }
}
