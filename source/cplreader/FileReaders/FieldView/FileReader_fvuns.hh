// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_FILEREADER_fvuns_2010
#define FILE_FILEREADER_fvuns_2010

#include "def_cplreader.hh"
#include "cplreader/FileReader.hh"

namespace CoupledField
{

  class FileReader_fvuns : public FileReader
  {
  public:

    //! Constructor
    FileReader_fvuns(const std::string& name,
                     const UInt dim,
                     const UInt numFiles,
                     const UInt startIndex);

    //! Deconstructor
    virtual ~FileReader_fvuns();

    virtual void Init();

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

    //! Get node groups
    virtual void GetNodeGroups(std::map<std::string,
                                        std::vector<UInt> >& nodeGroups);

  protected:
    std::vector<SolutionType> solutionTypes_;
          std::vector<UInt> dofIndices_;
  };
  

} // end of namespace
#endif
