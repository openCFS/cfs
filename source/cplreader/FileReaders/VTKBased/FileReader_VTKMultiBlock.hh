// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_FILEREADER_VTKMULTIBLOCK_2011
#define FILE_FILEREADER_VTKMULTIBLOCK_2011

#include "def_cplreader.hh"
#include "cplreader/FileReader.hh"

class vtkUnstructuredGrid;
class vtkMultiBlockDataSetAlgorithm;

namespace CoupledField
{

  class FileReader_VTKMultiBlock : public FileReader
  {
  public:

    //! Constructor
    FileReader_VTKMultiBlock(const std::string& name,
                        const UInt dim,
                        const UInt numFiles);

    //! Deconstructor
    virtual ~FileReader_VTKMultiBlock();

    virtual void Init();

    //! get node coordinates from the corresponding file
    virtual void ReadNodalCoords(std::vector<Double> & NODECOORD);


    //! get topology information from the corresponding topology file
    virtual void ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
                              std::vector<UInt> & elemTypes);

    /* get nodal values from the corresponding fluid datafile the new way */
    void ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                         const std::vector<bool>& activeParts,
                         const UInt timeStepIdx) = 0;

    /* get element values from the corresponding fluid datafile and interpolate to nodes */
    virtual void ReadElemValues(std::vector<FlowDataType>& nodalFlowData,
                                const std::vector<bool>& activeParts,
                                const UInt timeStepIdx){
      EXCEPTION("Not implemented in base class!");
    };

    virtual Double GetTimeStep(UInt stepNumber);
    virtual std::string GetRegionName(const UInt partitionIdx);

    //! get user data from file reader
    virtual void GetUserData(std::map<std::string, std::string>& userData) = 0;

    virtual void GetRegionElements(std::vector<UInt> & regionElements,
                                   const UInt regionIdx);

  protected:
    virtual Elem::FEType VTKCellTypeToFEType(UInt cellType);
    virtual void InitElemNodeMapping();

    virtual void CreateReader() = 0;
    virtual void EnableRegions() = 0;
    virtual void GetTimeValues() = 0;
    virtual void SetTimeValue(Double val) = 0;

    /* Actual time step values */
    std::vector<Double> timeValues_;
    std::vector< std::map<UInt, UInt> > nodeMap_;
    std::vector< std::vector<UInt> > actualRegionNodes_;
    std::vector< std::map< UInt, UInt > > actualRegionNodesToNewIdx_;

    std::map< UInt, std::vector<UInt> > unstrucElemNodeMapping_;
    std::map< UInt, std::vector<UInt> > uniformElemNodeMapping_;

    // Pointer for storing the reader
    vtkMultiBlockDataSetAlgorithm* reader_;

    // VTK dataset for merging nodes in order to get global node numbering
    vtkUnstructuredGrid* pts_;

    UInt numElems_;
  };


} // end of namespace
#endif
