// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_FILEREADER_GENGRIDS_2008
#define FILE_FILEREADER_GENGRIDS_2008

#include "def_cplreader.hh"
#include "cplreader/FileReader.hh"
#include "Utils/mathParser/mathParser.hh"

namespace CoupledField
{

  class FileReader_GENGRIDS : public FileReader
  {
  public:

    //! Constructor
    FileReader_GENGRIDS(const std::string& name,
                     const UInt dim,
                     const UInt numFiles,
                     const UInt startIndex);

    //! Deconstructor
    virtual ~FileReader_GENGRIDS();

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

    // Nodal coordinates for all nodes
    std::vector<Double> nodalCoords_;

    // Topology for per "original" element number
    std::vector<UInt> topology_;

    // Element type per "original" element number
    std::vector< UInt > elemTypes_;

    std::map<std::string, std::vector<UInt> > regionElems_;
    
    std::map<std::string, std::vector<UInt> > nodeGroups_;

    std::map<std::string, std::vector<std::string> > fieldExprs_;

    MathParser mathParser_;
  };


} // end of namespace
#endif
