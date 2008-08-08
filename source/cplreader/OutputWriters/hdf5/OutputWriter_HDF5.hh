#ifndef FILE_OUTPUTWRITER_HDF5_2008
#define FILE_OUTPUTWRITER_HDF5_2008

#include "H5Cpp.h"

#include "cplreader/OutputWriter.hh"

namespace CoupledField
{

  //! define CFS-HDF5 file format version
#define CFS_HDF5_FORMAT_MAJOR 0
#define CFS_HDF5_FORMAT_MINOR 9


  class OutputWriter_HDF5 : public OutputWriter {
  public:
    OutputWriter_HDF5();

    virtual ~OutputWriter_HDF5();

    // Init output writer
    virtual void Init(int argc, char** argv,
                      const std::string& outputPath);

    // Write nodal coordinates
    virtual void WriteNodalCoords(const std::vector<Double>& coords);

    // Write element connectivity/topology
    virtual void WriteTopology(UInt maxNumElemNodes,
                               const std::vector<UInt>& elemTypes,
                               const std::vector<UInt>& elems);

    // Write regions
    virtual void WriteRegion(std::string regionName,
                             const std::vector<UInt>& regionElems,
                             const std::vector<UInt>& regionNodes,
                             UInt regionDim);

    // Write node groups
    virtual void WriteNodeGroups(const std::map< std::string, std::vector<UInt> >&
                                 nodeGroups);

    virtual void WriteElemGroups(const std::map< std::string, std::vector<UInt> >&
                                 elemGroups,
                                 const std::map< std::string, std::vector<UInt> >&
                                 groupNodes,
                                 const std::map< std::string, UInt >&
                                 groupDims);
    

    virtual void WriteFlowData(UInt actRegion,
                               const std::vector<std::string>& outputFields);

    virtual void WriteUserData(const std::map<std::string, std::string>& userData);

    virtual void BeginResults(const std::vector<FlowDataType>* flowData);
    
    virtual void BeginStep(UInt stepNum, Double stepValue);
    virtual void EndStep();

    virtual void EndResults();

  private:
    // HDF5 specific member functions
    void CheckOpenObjects();
    void WriteFileInfoHeader();
    void InitHDF5();
    void InitResultsGroup();
    void WriteResultDescriptions();
    void WriteResults( H5::Group& resultGroup,
                       const std::vector<Double>& resultVals,
                       const UInt numDOFs,
                       const bool isImag );
    void CreateExternalFile(UInt timeStep);
    void WriteStringToUserData(const std::string& dSetName,
                               const std::string& str);

  private:
    // Directory name for storing HDF5 files
    std::string hdf5DirName_;

    // HDF5 specific members
    //! Dataset property list (used for chunking and compression )
    H5::DSetCreatPropList dPropList_;

    //! Main file containing grid and meta information
    H5::H5File mainFile_;

    //! In case we use
    H5::H5File currStepFile_;

    //! Main / Root Group
    H5::Group mainGroup_;

    //! Mesh Group
    H5::Group meshGroup_;

    //! Group for results
    H5::Group resultsGroup_;

    //! Group for user data
    H5::Group userDataGroup_;

    //! Group for mesh results
    H5::Group meshResultsGroup_;

    //! Group for current analysis step for mesh results
    H5::Group currMeshStepGroup_;

    std::vector< std::string > regionNames_;
  };

} // end of namespace
#endif
