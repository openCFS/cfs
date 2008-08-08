#ifndef FILE_FILEREADER_2006
#define FILE_FILEREADER_2006

#include <fstream>
#include <vector>

#include <def_cplreader.hh>
#include "FlowDataTypes.hh"

//! Base class for reading topology and nodal data from results of fluid computations
/*!
  This class contains methods for the reading information
  from a fluid topology file and from the fluid results files.
  It reads and couples flow results from a set of files (or set of
  region files coming from a blocked structured fluid computation).
*/

namespace CoupledField
{

  class FileReader
  {
  public:

    //! Constructor
    FileReader(const std::string& name, const UInt dim, const UInt numFiles);

    //! Deconstructor
    virtual ~FileReader();

    virtual void Init() = 0;

    //! get node coordinates from the corresponding file
    virtual void ReadNodalCoords(std::vector<Double> & NODECOORD) = 0;


    //! get topology information from the corresponding topology file
    virtual void ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
                                 std::vector<UInt> & elemTypes) = 0;

    //! get elements in a certain region.
    virtual void GetRegionElements(std::vector<UInt> & regionElements,
                                   const UInt regionIdx)
    {
      std::cerr << "GetRegionElements not implemented" << std::endl;
    }

    //! get nodal values from the corresponding fluid datafile the new way
    virtual void ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                                 const std::vector<bool>& activeParts,
                                 const UInt timeStepIdx)
    {
      std::cerr << "ReadNodalValues (FlowDataType& nodalFlowData, const UInt timeStepIdx) " <<
        "not implemented!\n\n";

      exit(1);
    }

    //! Return number of regions
    virtual UInt GetNumRegions() { return numRegions_;}

    //! Return number number of timestep files
    virtual UInt GetNumFiles() { return numSteps_;}

    //! Return maximum number of nodes
    virtual UInt GetNumNodes(const UInt regionIdx) {
      return numNodesPerRegion_[regionIdx];
    }

    //! Return maximum number of nodes
    virtual UInt GetNumElems(const UInt regionIdx) {
      return numElemsPerRegion_[regionIdx];
    }

    //! return dimension of grid
    virtual UInt GetDim() { return dim_;}

    //! Get node groups
    virtual void GetNodeGroups(std::map<std::string,
                                        std::vector<UInt> >& nodeGroups);

    //! Get element groups
    virtual void GetElemGroups(std::map<std::string,
                                        std::vector<UInt> >& elemGroups);

    //! return time step value in seconds
    virtual double GetTimeStep(UInt stepNumber);

    //! return name of region
    virtual std::string GetRegionName(const UInt regionIdx);

    //! get user data from file reader
    virtual void GetUserData(std::map<std::string, std::string>& userData) {};

    UInt GetMaxNumElemNodes() {
      return maxNumElemNodes_;
    }

    std::string GetPreferredOutputPath() {
      return preferredOutputPath_;
    }

  protected:


    //! ifstream which are associated with the fluid files
//    std::ifstream flowdatafile;
    std::string name_;
    std::string baseName_;

    //! Number of regions
    UInt numRegions_;

    //! dimension of grid
    UInt dim_;

    //! Number of nodes per region
    std::vector<UInt> numNodesPerRegion_;

    //! Number of elements per Region
    std::vector<UInt> numElemsPerRegion_;

//    UInt numResults_;

    //! Number of time steps
    UInt numSteps_;

    // Maximum number of nodes per element
    UInt maxNumElemNodes_;

    // Preferred path for output of HDF5 files.
    std::string  preferredOutputPath_;

    std::map<SolutionType, bool> requiredResults_;
  };


} // end of namespace
#endif
