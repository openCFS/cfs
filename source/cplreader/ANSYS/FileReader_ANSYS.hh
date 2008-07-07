#ifndef FILE_FILEREADER_ANSYS_2008
#define FILE_FILEREADER_ANSYS_2008

#include <cplreaderdefs.hh>
#include "FileReader.hh"

namespace CoupledField
{

  class FileReader_ANSYS : public FileReader
  {
  public:

    //! Constructor
    FileReader_ANSYS(const std::string& name,
                     const UInt dim,
                     const UInt numFiles,
                     const UInt startIndex);

    //! Deconstructor
    virtual ~FileReader_ANSYS();

    // Determine type of reader necessary for reading given files.
    static std::string GetReaderType();

    //! get node coordinates from the corresponding file
    virtual void ReadNodalCoords(std::vector<Double> & NODECOORD);


    //! get topology information from the corresponding topology file
    virtual void ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
                                 std::vector<UInt> & elemTypes);

    //! get nodal values from the corresponding fluid datafile the new way
    virtual void ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                                 const std::vector<bool>& activeParts,
                                 const UInt timeStepIdx);

    virtual std::string GetRegionName(const UInt regionIdx);

    virtual void GetRegionElements(std::vector<UInt> & regionElements,
                                   const UInt regionIdx);

  protected:

    void OpenFile(std::string fn);
    bool ReadChunk();
    bool GetNextLine(std::string& line);

    virtual void ResortNodes(std::vector<UInt>& elemNodes);
    virtual void DegenerateElement(const FEType elemTypeIn,
                                   FEType& elemTypeOut,
                                   std::vector<UInt>& elemNodes);

    // Checks if connectivities of two elements match.
    bool CompareElements(const std::vector<UInt>& elemNodes1,
                         const std::vector<UInt>& elemNodes2);

    UInt chunkSize_;
    UInt fSize_;
    bool strict_;
    bool degen_;
    bool everyThingRead_;

    std::ifstream inFile_;

    std::vector<std::string> regionNames_;
    std::vector<std::string> lines_;
    std::vector<std::string>::iterator lineIt_, lineEnd_;

    // Nodal coordinates for all nodes
    std::vector<Double> nodalCoords_;

    // Topology for per "original" element number
    std::map<UInt, std::vector<UInt> > topology_;

    // Element type per "original" element number
    std::map<UInt, UInt > elemTypes_;

    // Vector of "original" element numbers per region
    std::map<UInt, std::vector<UInt> > elemNumsOrig_;

    // Map from original node number to new node number
    std::map<UInt, UInt > nodeNumsMap_;

    // Map from original element number to new element number
    std::map<UInt, UInt > elemNumsMap_;

    // Map from original element number to region index
    std::map<UInt, UInt > elemRegionMap_;

    // Topology for per "original" element number
    std::vector<std::string> fileNames_;

  };


} // end of namespace
#endif
