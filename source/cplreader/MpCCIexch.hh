#ifndef FILE_MPCCIEXCH_2006
#define FILE_MPCCIEXCH_2006

#include "H5Cpp.h"

#include <cplreaderdefs.hh>
#include "filereader.hh"

namespace CoupledField
{

//! define CFS-HDF5 file format version
#define CFS_HDF5_FORMAT_MAJOR 0
#define CFS_HDF5_FORMAT_MINOR 9

  //! Class for mesh coupling
  /*! 
    This class handles the subroutines calls concerning MpCCI for coupling the fluid and acoustic computations.
  */

  class MpCCIExchangeCPLR
  {
  public:

    //!
    MpCCIExchangeCPLR(FileReader * ptFileReader);

    //!
    virtual ~MpCCIExchangeCPLR();

    // Perform initialization
    void Init(int argc, char *argv[]);
    
    //! Reorganizing grid info for MpCCi and hand over to MpCCI
    void PutExchangeGrid2MpCCI();

    //! Performs the coupled computation phase
    void Couple();

    // Perform last operations
    void Finish();

  private:

    void CalculateAcouSrcs(const int partitionIdx,
                           FlowDataType& flowData);
    

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
                                 const std::vector<std::string>& outputFields,
                                 const std::vector<UInt> stepNumbers,
                                 const std::vector<Double> stepValues,
                                 const std::vector<std::string> regions);
    void WriteResults( H5::Group& resultGroup,
                       std::vector<Double>& resultVals,
                       const UInt numDOFs,
                       const bool isImag );
    void CreateExternalFile(UInt timeStep);
    void WriteStringToUserData(const std::string& dSetName, 
                               const std::string& str);

    FileReader *  ptFileReader_;
    std::map<UInt, std::vector<Double> > NodalCoords_;
    std::map<UInt, std::vector<UInt> > Topology_;

    std::map<UInt, std::vector<UInt> > numNodesPerElem_;
    std::map<UInt, std::vector<UInt> > elemTypes_;

    std::map<UInt, std::map<UInt, UInt> > nodesInPartition_;

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

  };

} // end of namespace
#endif
