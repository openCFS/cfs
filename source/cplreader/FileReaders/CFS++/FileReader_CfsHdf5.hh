#ifndef FILE_FILEREADER_CFSDF5_2010
#define FILE_FILEREADER_CFSDF5_2010

#include "def_cplreader.hh"
#include "cplreader/FileReader.hh"
#include "hdf5Reader.h"
#include "hdf5Common.h"

namespace CoupledField
{

  class FileReader_CfsHdf5 : public FileReader
  {
  public:

    //! Constructor
    FileReader_CfsHdf5(const std::string& name,
                     const UInt dim,
                     const UInt numFiles,
                     const UInt startIndex);

    //! Deconstructor
    virtual ~FileReader_CfsHdf5();

    virtual void Init();

    virtual UInt GetNumRegions();

    void ReadNodalCoords(std::vector<Double>& nodalCoords);

    //! return name of region
    virtual std::string GetRegionName(const UInt regionIdx);

    //! get topology information from the corresponding topology file
    virtual void ReadTopology(std::vector<UInt> & connectivities,
                              std::vector<UInt> & elemTypes);

    //! get elements in a certain region.
    virtual void GetRegionElements(std::vector<UInt> & regionElements,
                                   const UInt regionIdx);

    //! get nodal values from the corresponding fluid datafile the new way
    virtual void ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                                 const std::vector<bool>& activeParts,
                                 const UInt timeStepIdx);

    //! return time step value in seconds
    virtual double GetTimeStep(UInt stepNumber);

    //! Correct the numbering of the nodes. This is needed if nodes are tossed
    //! away, e.g. with the reduce_elementOrder_ flag
    virtual void CorrectNumbering(std::vector<Double>& nodalCoords,
                                  std::vector<UInt>& connectivities, \
                                  std::vector<UInt>& elemTypes);

    //! Reduce the order of the result vectors. This may only be called if
    //! CorrectNumbering has already been called since a map from reduced grid
    //! to original grid is necessary.
    //! @param nodalFlowData The FlowData which has all the results which should
    //! be reduced
    //! @param regionNodes The nodes in the reduced region. The results are
    //! store per each region that is why it is necessary to know which nodes
    //! are in which region. Furthermore this is why CorrectNumbering() needs to
    //! be called which creates a map from the global original grid to the
    //! global reduced grid. Additional requirements are calculated and
    //! determined internaly.
    void ReduceOrderOfNodalValues(std::vector<FlowDataType>& nodalFlowData, \
        const std::map<std::string, std::vector<UInt> >& regionNodes);

  private:
    H5CFS::Hdf5Reader hdf5Reader_;
    std::vector<std::string> regionNames_;
    std::map<unsigned int, double> timeStepValues_;
    UInt reduce_elementOrder_;
    std::vector<UInt> elemTypes_;
    std::map<UInt,UInt> nodesMap_;
    std::vector<UInt> nodesMapAsVec_;
    std::map<RegionIdType, std::map<UInt, UInt> > origRegToNewReg_;

    /** The method creates a map from position in original vector to position in
     * reduced vector. To be more precise:
     * There are different node numbering vector, for the original grid there is 
     * the whole node numbering (global nodes) and vectors wich have the node
     * numbers for each region. The same goes for the new reduced grid. It is
     * now necessesary to find out which regional node is to be found were in
     * both grids. This is necessary to copy the correct solution data from one
     * to the other vector 
     * @param resultMap The map which will have the wanted map
     * @param nodesMapAsVec This is a map from reduced grid to the original grid
     * @param origRegNodes The node numbers of the original grid
     * @param regionNodes The node numbers of the reduced grid
     **/
    void mapRegionalOrigToReduce(std::map<UInt, UInt>& resultMap,
                                 const std::vector<UInt>& nodesMapAsVec,
                                 const std::vector<UInt>& origRegNodes,
                                 const std::vector<UInt>& regionNodes)
    {
      for (UInt i = 0; i < regionNodes.size(); ++i)
      {
        const UInt& currentRegNode = regionNodes[i];
        const UInt& currentGlobalNode = nodesMapAsVec[currentRegNode];
        std::vector<UInt>::const_iterator origRegNodesIter = \
               std::find(origRegNodes.begin(), origRegNodes.end(), currentGlobalNode);
        UInt origRegNodePosition = origRegNodesIter - origRegNodes.begin();
        resultMap[i] = origRegNodePosition;
      }
    }
  };


} // end of namespace
#endif
