// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_FILEREADER_OPENFOAM_2007
#define FILE_FILEREADER_OPENFOAM_2007

#include <def_cplreader.hh>
#include "cplreader/FileReader.hh"

class vtkOpenFOAMReader;
class vtkUnstructuredGrid;

namespace CoupledField
{

  class FileReader_OPENFOAM : public FileReader
  {
  public:

    //! Constructor
    FileReader_OPENFOAM(const std::string& name,
                        const UInt dim,
                        const UInt numFiles);

    //! Deconstructor
    virtual ~FileReader_OPENFOAM();

    virtual void Init();

    //! get node coordinates from the corresponding file
    virtual void ReadNodalCoords(std::vector<Double> & NODECOORD);


    //! get topology information from the corresponding topology file
    virtual void ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
                              std::vector<UInt> & elemTypes);

    /* get nodal values from the corresponding fluid datafile the new way */
    void ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                         const std::vector<bool>& activeParts,
                         const UInt timeStepIdx);

    virtual Double GetTimeStep(UInt stepNumber);
    virtual std::string GetRegionName(const UInt partitionIdx);

    //! get user data from file reader
    virtual void GetUserData(std::map<std::string, std::string>& userData);

    virtual void GetRegionElements(std::vector<UInt> & regionElements,
                                   const UInt regionIdx);

  protected:
      Elem::FEType VTKCellTypeToFEType(UInt cellType);

      /* Actual time step values */
      std::vector<Double> timeValues_;
      vtkOpenFOAMReader* reader_;
      std::vector< std::map<UInt, UInt> > nodeMap_;
      std::vector< std::vector<UInt> > actualRegionNodes_;
    vtkUnstructuredGrid* pts_;

      UInt numElems_;
  };


} // end of namespace
#endif
