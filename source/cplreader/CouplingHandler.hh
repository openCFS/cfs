#ifndef FILE_COUPLINGHANDLER_2006
#define FILE_COUPLINGHANDLER_2006

#include "H5Cpp.h"

#include <def_cplreader.hh>
#include "FileReader.hh"
#include "OutputWriter.hh"

namespace CoupledField
{

//! define CFS-HDF5 file format version
#define CFS_HDF5_FORMAT_MAJOR 0
#define CFS_HDF5_FORMAT_MINOR 9

  class ElemIntegr;

  //! Class for mesh coupling
  /*!
    This class handles the subroutines calls concerning MpCCI for coupling the fluid and acoustic computations.
  */

  class CouplingHandler
  {
  public:

    //!
    CouplingHandler(shared_ptr<FileReader> ptFileReader,
                      std::vector< shared_ptr<OutputWriter> >& outputWriters);

    //!
    virtual ~CouplingHandler();

    // Perform initialization
    void Init(int argc, char *argv[]);

    //! Hand over mesh from reader to writer
    void ConvertMesh();

    //! Performs the coupled computation phase
    void Couple();

    // Perform last operations
    void Finish();

  private:

    void CalculateAcouSrcs(const int partitionIdx,
                           FlowDataType& flowData);

    // Throw away unneccessary entries in vectors
    void ShrinkNodalVector(const UInt partitionIdx,
                           const UInt numDOFs,
                           const std::vector<Double>& input,
                           std::vector<Double>& output);

    void CollectElementNodes(const std::vector<UInt>& elems,
                             std::vector<UInt>& nodes,
                             UInt& dim);

    // MpCCI specific member functions
    int ElemTypes2MpCCI(FEType et);

    // HDF5 specific member functions
    void CheckOpenObjects();
    void WriteFileInfoHeader();
    void InitHDF5();
    void WriteRegion(const H5::Group& meshGroup,
                     const std::vector<UInt>& nodes,
                     const std::vector<UInt>& elems,
                     const UInt dim,
                     const std::string name);
    void InitResultsGroup();
    void WriteResultDescriptions(UInt numSteps,
                                 const std::vector<FlowDataType>& outputFields,
                                 const std::vector<UInt> stepNumbers,
                                 const std::vector<Double> stepValues);
    void WriteResults( H5::Group& resultGroup,
                       std::vector<Double>& resultVals,
                       const UInt numDOFs,
                       const bool isImag );
    void CreateExternalFile(UInt timeStep);
    void WriteStringToUserData(const std::string& dSetName,
                               const std::string& str);

    void WriteNodeGroups(const H5::Group& meshGroup);
    void WriteElemGroups(const H5::Group& meshGroup);

    shared_ptr<FileReader>  ptFileReader_;
    std::vector<Double> nodalCoords_;
    std::vector<UInt> topology_;
    std::vector<UInt> elemTypes_;

    std::map<UInt, std::vector<UInt> > regionElems_;
    std::map<UInt, std::map<UInt, UInt> > regionNodeIndices_;
    std::vector<UInt> numRegionNodes_;

    UInt dim_;

    std::map<std::string, std::vector<UInt> > nodeGroups_;
    std::map<std::string, std::vector<UInt> > elemGroups_;

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

    //!MpCCI
    int meshId_;
    int GlobalDim_;
    int MpCCIprocess_;
    int quantityId1_;
    int quantityDim1_;
    int quantityId2_;
    int quantityDim2_;

    std::map<UInt, ElemIntegr *> ptElemIntegr_;
    std::vector<std::string> outputFields_;
    std::vector<bool> activeParts_;

    std::vector< shared_ptr<OutputWriter> > outputWriters_;
  };

} // end of namespace
#endif
