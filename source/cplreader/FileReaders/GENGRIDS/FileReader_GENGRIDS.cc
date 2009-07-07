// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "Settings.hh"

#include "SphericalShellGenerator.hh"
#include "RectilinearUniformGridGenerator.hh"
#include "ReferenceElementGenerator.hh"

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
    ReferenceElementGenerator refGen;
    
    Settings& settings = Settings::Instance();

    ParamNode *param = settings.GetParamNode();

    ParamNode* typeNode = NULL;

    if(param) 
    {
      typeNode = param->Get("type");
      if(typeNode->Get("name")->AsString() == "GENGRIDS") 
      {
        std::string subType = "RectilinearUniformGridGenerator";
        //        ParamNode* subTypeNode = typeNode->Get("subType", false);
        typeNode->Get("subType", subType, false);

        if(subType == "RectilinearUniformGridGenerator") 
        {
          ugGen.GenUniformGrid(nodalCoords_,
                               topology_,
                               elemTypes_,
                               maxNumElemNodes_,
                               regionElems_,
                               nodeGroups_,
                               RectilinearUniformGridGenerator::LAGRANGE);
        }

        if(subType == "ReferenceElementGenerator") 
        {
          std::string elemTypeName = "HEXA8";
          ParamNode* elemTypeNode = typeNode->Get("elementType", false);
          
          if(elemTypeNode) {
            elemTypeNode->Get("name", elemTypeName, false);
          }
          
          refGen.GenGrid(nodalCoords_,
                         topology_,
                         elemTypes_,
                         maxNumElemNodes_,
                         regionElems_,
                         nodeGroups_,
                         Elem::feType.Parse(elemTypeName));
        }

        if(subType == "SphericalShellGenerator") 
        {
          shellGen.GenSphericalShell(nodalCoords_,
                                     topology_,
                                     elemTypes_,
                                     maxNumElemNodes_,
                                     regionElems_,
                                     nodeGroups_);
        }
      }
    }
    

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

