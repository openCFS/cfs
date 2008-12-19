// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_FILEREADER_FASTEST_2006
#define FILE_FILEREADER_FASTEST_2006

#include <def_cplreader.hh>
#include "cplreader/FileReader.hh"

namespace CoupledField
{

  class FileReader_FASTEST : public FileReader
  {
  public:

    //! Constructor
    FileReader_FASTEST(const std::string& name,
                       const UInt dim,
                       const UInt numFiles,
                       const UInt startIndex);

    //! Deconstructor
    virtual ~FileReader_FASTEST();

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

    virtual void GetRegionElements(std::vector<UInt> & regionElements,
                                   const UInt regionIdx);

  protected:

    UInt startIndex_;
    std::string partFmtStr_;
    std::string timeFmtStr_;
    std::vector<Integer> dataColumns_;
    std::vector<SolutionType> solutionTypes_;
    std::vector<UInt> dofIndices_;
    std::vector<UInt> elsize_;
    UInt numResults_;

    std::ifstream inFile_;

  };


} // end of namespace
#endif
