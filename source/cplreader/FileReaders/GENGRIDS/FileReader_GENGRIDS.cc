// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "Settings.hh"

#include "SphericalShellGenerator.hh"
#include "RectilinearUniformGridGenerator.hh"
#include "ReferenceElementGenerator.hh"
#include "Utils/trivialCartesianCoordSys.hh"

#include "FileReader_GENGRIDS.hh"

namespace CoupledField
{

  FileReader_GENGRIDS::FileReader_GENGRIDS(const std::string& name,
                                     const UInt dim,
                                     const UInt numFiles,
                                     const UInt startIndex) :
    FileReader(name, dim, numFiles)
  {
    std::cout << "FileReader_GENGRIDS " << name << " " << dim << " " << numFiles << std::endl;
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
    ParamNode* fieldsNode = NULL;

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

        fieldsNode = typeNode->Get("fields", false);
        if(fieldsNode) 
        {
          StdVector<ParamNode*> fields = fieldsNode->GetList("field");

          for(UInt i=0; i<fields.GetSize(); i++) 
          {
            std::string fieldName = "default";
            //            ParamNode* elemTypeNode = typeNode->Get("elementType", false);
            fields[i]->Get("name", fieldName, false);
            std::cout << "field: " << fieldName << std::endl;

            StdVector<ParamNode*> components = fields[i]->GetList("component");
            fieldExprs_[fieldName].resize(components.GetSize());
            
            for(UInt j=0; j<components.GetSize(); j++) 
            {
              std::string compAxis = "x";
              std::string expr;
              
              components[j]->Get("axis", compAxis, false);
              expr = components[j]->AsString();

              for(UInt l=0; l<expr.length(); l++) 
              {
                if(expr[l] == '\n' ||
                   expr[l] == '\r' ||
                   expr[l] == '\t')
                  expr[l] = ' ';
              }
            
              if(compAxis == "x")
                fieldExprs_[fieldName][0] = expr;
              if(compAxis == "y")
                fieldExprs_[fieldName][1] = expr;
              if(compAxis == "z")
                fieldExprs_[fieldName][2] = expr;
              
              std::cout << "component: " << compAxis << " " << expr << std::endl;
              
            }
          }

          if(!fields.GetSize())
            numSteps_ = 0;            
        }
        else 
        {
          numSteps_ = 0;
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

  void FileReader_GENGRIDS::ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                                            const std::vector<bool>& activeParts,
                                            const UInt timeStepIdx)
  {
    UInt regionIdx, nRegions;
    Settings& settings = Settings::Instance();
    Double timeStep = timeStepIdx * settings.GetDouble("timestep");
    TrivialCartesianCoordSystem coordSys("cplreaderCoordSys", NULL, NULL);
    Vector<Double> globCoord(3);

    std::map<std::string, std::vector<UInt> >::iterator it;
    it = regionElems_.begin();

    MathParser::HandleType handle = mathParser_.GetNewHandle();
    
    for( regionIdx = 0, nRegions = nodalFlowData.size(); regionIdx < nRegions; regionIdx++, it++ ) 
    {
      if(!activeParts[regionIdx])
        continue;

      std::set<UInt> regionNodeSet;
      
      for( UInt el=0, n = it->second.size(); el < n; el++) 
      {
        UInt elemIdx = it->second[el]-1;
        regionNodeSet.insert(&topology_[maxNumElemNodes_*elemIdx],
                             &topology_[maxNumElemNodes_*(elemIdx+1)]);
      }

      UInt numRegionNodes = regionNodeSet.size();
      
      std::map<std::string, std::vector<std::string> >::iterator exprIt, exprEnd;
      exprIt = fieldExprs_.begin();
      exprEnd = fieldExprs_.end();

      for( ; exprIt != exprEnd; exprIt++) 
      {  
        SolutionType solType;
        if(exprIt->first == "fluidMechPressure") {
          solType = FLUIDMECH_VELOCITY;
        } else if(exprIt->first == "fluidMechVelocity") {
          solType = FLUIDMECH_VELOCITY;
        } else if(exprIt->first == "acouRhsLoad") {
          solType = ACOU_RHS_LOAD;
        } else {
          EXCEPTION("Solution type '" << exprIt->first << "' not supported!");
        }

        UInt numComps = exprIt->second.size();
        
        nodalFlowData[regionIdx][solType].resultName = exprIt->first;
        nodalFlowData[regionIdx][solType].isActive = true;
        nodalFlowData[regionIdx][solType].definedOn = ResultInfo::NODE; // nodes
        switch(exprIt->second.size()) 
        {
        case 1:
          nodalFlowData[regionIdx][solType].entryType = ResultInfo::SCALAR;
        case 3:
          nodalFlowData[regionIdx][solType].entryType = ResultInfo::VECTOR;
        }
        nodalFlowData[regionIdx][solType].dofNames.resize(exprIt->second.size());
        nodalFlowData[regionIdx][solType].data.resize(numRegionNodes * numComps);
        
        
        for(UInt comp=0; comp<numComps; comp++) 
        {  
          mathParser_.SetExpr(handle, exprIt->second[comp]);
          mathParser_.SetValue(handle, "t", timeStep);
          mathParser_.SetValue(handle, "t_idx", timeStepIdx);

          std::set<UInt>::iterator nsIt, nsEnd;
          nsIt= regionNodeSet.begin();
          nsEnd= regionNodeSet.end();

          for(UInt nodeIdx = 0; nsIt != nsEnd; nsIt++, nodeIdx++) 
          {
            if(!*nsIt) 
            {
              nodeIdx--;
              continue;
            }
            
            //            std::cout << "nodeIdx " << nodeIdx << " *nsIt " << *nsIt << std::endl;
            globCoord[0] = nodalCoords_[(*nsIt-1)*3+0];
            globCoord[1] = nodalCoords_[(*nsIt-1)*3+1];
            globCoord[2] = nodalCoords_[(*nsIt-1)*3+2];
            
            mathParser_.SetCoordinates(handle, coordSys, globCoord);
            nodalFlowData[regionIdx][solType].data[nodeIdx*numComps+comp] =  mathParser_.Eval(handle);
            //            std::cout << "data " << nodalFlowData[regionIdx][solType].data[nodeIdx*numComps+comp] << std::endl;
          }
          
        }
      }
      
    }
  }

}

