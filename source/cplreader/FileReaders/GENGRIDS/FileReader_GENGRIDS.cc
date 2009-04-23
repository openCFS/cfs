// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "SphericalShellGenerator.hh"
#include "RectilinearUniformGridGenerator.hh"

#include "FileReader_GENGRIDS.hh"

namespace CoupledField
{

  FileReader_GENGRIDS::FileReader_GENGRIDS(const std::string& name,
                                     const UInt dim,
                                     const UInt numFiles,
                                     const UInt startIndex) :
    FileReader(name, dim, numFiles)
  {
    numSteps_ = 0;
  }

  FileReader_GENGRIDS::~FileReader_GENGRIDS()
  {
  }

  void FileReader_GENGRIDS::Init()
  {
    SphericalShellGenerator shellGen;
    RectilinearUniformGridGenerator ugGen;

    ugGen.GenUniformGrid(nodalCoords_,
                         topology_,
                         elemTypes_,
                         maxNumElemNodes_,
                         regionElems_,
                         nodeGroups_,
                         RectilinearUniformGridGenerator::SERENDIPITY);

 #if 0    
    shellGen.GenSphericalShell(nodalCoords_,
                               topology_,
                               elemTypes_,
                               maxNumElemNodes_,
                               regionElems_,
                               nodeGroups_);
 #endif

    numRegions_ = regionElems_.size();
    dim_ = 3;
    numNodesPerRegion_.resize(numRegions_);
    numElemsPerRegion_.resize(numRegions_);

    for(UInt i=0; i<numRegions_; i++) 
    {
      numNodesPerRegion_[i] = 1;
      numElemsPerRegion_[i] = elemTypes_.size() / 6;
    }
  }
  
  void FileReader_GENGRIDS::ReadNodalCoords(std::vector<Double> & NODECOORD)
  {
    NODECOORD = nodalCoords_;
  }

  void FileReader_GENGRIDS::ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
                                      std::vector<UInt> & elemTypes)
  {
    TOPOLOGYDATA = topology_;
    elemTypes = elemTypes_;
  }

  void FileReader_GENGRIDS::GetRegionElements(std::vector<UInt> & regionElements,
                                           const UInt regionIdx)
  {
    std::map<std::string, std::vector<UInt> >::iterator it;
    it = regionElems_.begin();
    
    for(UInt i=0; i<regionIdx; i++, it++);

    regionElements = it->second;
  }

  std::string FileReader_GENGRIDS::GetRegionName(const UInt regionIdx)
  {
    std::map<std::string, std::vector<UInt> >::iterator it;
    it = regionElems_.begin();
    
    for(UInt i=0; i<regionIdx; i++, it++);

    return it->first;    
  }

  void FileReader_GENGRIDS::GetNodeGroups(std::map<std::string,
                                          std::vector<UInt> >& nodeGroups) 
  {
    nodeGroups = nodeGroups_;
  }
}

