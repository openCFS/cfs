#ifndef FILE_FILEREADER_CFSDF5_2010
#define FILE_FILEREADER_CFSDF5_2010

#include <def_cplreader.hh>
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

  private:
    H5CFS::Hdf5Reader hdf5Reader_;
    std::vector<std::string> regionNames_;
    std::map<unsigned int, double> timeStepValues_;
    UInt firststep_;
  };


} // end of namespace
#endif
