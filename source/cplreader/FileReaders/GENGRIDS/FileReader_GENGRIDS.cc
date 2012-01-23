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
    FileReader(name, dim, numFiles, startIndex)
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

    PtrParamNode param = settings.GetParamNode();

    PtrParamNode typeNode;
    PtrParamNode fieldsNode;

    if(param) 
    {
      typeNode = param->Get("type");
      if(typeNode->Get("name")->As<std::string>() == "GENGRIDS") 
      {
        std::string subType = "RectilinearUniformGridGenerator";
        //        PtrParamNode subTypeNode = typeNode->Get("subType", false);
        typeNode->GetValue("subType", subType, ParamNode::PASS);

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
          PtrParamNode elemTypeNode = typeNode->Get("elementType", ParamNode::PASS);
          
          if(elemTypeNode) {
            elemTypeNode->GetValue("name", elemTypeName, ParamNode::PASS);
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

        fieldsNode = typeNode->Get("fields", ParamNode::PASS);
        if(fieldsNode) 
        {
          ParamNodeList fields = fieldsNode->GetList("field");

          for(UInt i=0; i<fields.GetSize(); i++) 
          {
            std::string fieldName = "default";
            //            PtrParamNode elemTypeNode = typeNode->Get("elementType", false);
            fields[i]->GetValue("name", fieldName, ParamNode::PASS);
            std::cout << "field: " << fieldName << std::endl;

            ParamNodeList components = fields[i]->GetList("component");
            fieldExprs_[fieldName].resize(components.GetSize());
            
            for(UInt j=0; j<components.GetSize(); j++) 
            {
              std::string compAxis = "x";
              std::string expr;
              
              components[j]->GetValue("axis", compAxis, ParamNode::PASS);
              expr = components[j]->As<std::string>();

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
    TrivialCartesianCoordSystem coordSys("cplreaderCoordSys", NULL, PtrParamNode());
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
        } else if(exprIt->first == "acouLambRhs") {
          solType = ACOU_LAMB_RHS;
        } else if(exprIt->first == "acouLambVec") {
          solType = FLUIDMECH_VELOCITY;
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
        
        
        if(exprIt->second[0] == "vortexPair"){
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
            globCoord[0] = nodalCoords_[(*nsIt-1)*3+0];
            globCoord[1] = nodalCoords_[(*nsIt-1)*3+1];
            globCoord[2] = nodalCoords_[(*nsIt-1)*3+2];
            std::vector<Double> nodeVel;
            //commented out clean implementation
            //if( solType == ACOU_LAMB_VEC){
            //  ComputeVortexPairSource(globCoord, timeStep,nodeVel);
            //}else if(solType == FLUIDMECH_VELOCITY){
             ComputeVortexPairVelocity(globCoord, timeStep,nodeVel);
            //}else{
            //  EXCEPTION("Solution type '" << solType << "' not supported for analytical caluculation!");
            //}
            nodalFlowData[regionIdx][solType].data[nodeIdx*numComps+0] =  nodeVel[0];
            nodalFlowData[regionIdx][solType].data[nodeIdx*numComps+1] =  nodeVel[1];
            nodalFlowData[regionIdx][solType].data[nodeIdx*numComps+2] =  nodeVel[2];
          }
        }else{
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
  inline void FileReader_GENGRIDS::ComputeVortexPairSource(Vector<Double> coord, Double t, std::vector<Double> & velos){
    //global coordinates
    Double x,y,gamma,omega,eFac;

    x = coord[0];
    y = coord[1];
    gamma = 1.00531;
    omega = 0.08;
    eFac = 0.5;


    velos.resize(3);

    Double exp1 = exp(-eFac*(pow((x+cos(omega*t)),2) + pow((y+sin(omega*t)),2)));
    Double exp2 = exp(-eFac*(pow((x-cos(omega*t)),2) + pow((y-sin(omega*t)),2)));

    Double constant = pow(gamma,2)/(8*pow(PI,2))*(exp1-exp2);

    velos[0] = constant * cos(omega*t);
    velos[1] = constant * sin(omega*t);
    velos[2] = 0;
  }
  
  inline void FileReader_GENGRIDS::ComputeVortexPairVelocity(Vector<Double> coord, Double t, std::vector<Double> & velos){
    //global coordinates
    Double x,y;

    //user parameters
    Double bulkModulus,rho0,r0,gamma,rc,modelRad;
    
    x = coord[0];
    y = coord[1];
    //we know that this is a 2D probolem so we ignore the 3rd coord!
    
    bulkModulus = 1.0;
    rho0        = 1.0;
    r0          = 1.0;
    gamma       = 1.00531;
    rc          = 0.1;
    modelRad    = 0.15;

    //Compute remaining quantities
    //Double c0 = sqrt(bulkModulus/rho0);
    Double omega = gamma/(4.0*PI*r0*r0);
    std::complex<Double> Im(0.0,1.0);

    //compute radius from center
    //Double r = sqrt(x*x+y*y);
    Double arg = omega*t;

    std::complex<Double> z( x, y);
    std::complex<Double> b( r0 * cos(arg), r0*sin(arg));
    std::complex<Double> distance(pow(z,2) - pow(b,2));

    //compute complex value of velocity components
    std::complex<double> w_z(gamma*z/(PI*Im*distance));
    //compute complex derivative of velocity components
    std::complex<double> w_zz(-gamma*(pow(z,2) + pow(b,2)) / (PI*Im*pow(distance,2)));

    Double x_vort = 0;
    Double y_vort = 0;
    velos.resize(3);
    if( (abs(z-b) <= modelRad) || (abs(z+b) <= modelRad) ){
      if(abs(z-b) <= modelRad){
        y_vort = y - imag(b);
        x_vort = x - real(b);
      }
      if(abs(z+b) <= modelRad){
        y_vort = y + imag(b);
        x_vort = x + real(b);
      }
      velos[0] = -1.35 * gamma/(2*PI) * y_vort/(rc*rc + x_vort*x_vort + y_vort*y_vort);
      velos[1] =  1.35 * gamma/(2*PI) * x_vort/(rc*rc + x_vort*x_vort + y_vort*y_vort);
      velos[2] =  0;
    }else{
      velos[0] = real(w_z);
      velos[1] = -imag(w_z);
      velos[2] =  0;
    }

    //Double sigma = r0;
    //Double constant = gamma/(8*PI*r0*sigma*sigma);
    //constant *= (exp(-1.0/(2.0*pow(sigma,2)) * (pow(x+r0*cos(omega*t),2) + pow(y+r0*sin(omega*t),2)) ) -
    //             exp(-1.0/(2.0*pow(sigma,2)) * (pow(x-r0*cos(omega*t),2) + pow(y-r0*sin(omega*t),2)) ) );
    //velos[0] = constant * cos(omega*t); 
    //velos[1] = constant * sin(omega*t);
    //now compute the lambvector
    //std::vector<Double> velTmp = velos;
    //Matrix<Double> veloDeriv(3,3);
    //veloDeriv.Init();
    //veloDeriv[0][0] = real(w_zz);
    //veloDeriv[1][0] = -imag(w_zz);
    //veloDeriv[0][1] = veloDeriv[1][0];
    //veloDeriv[1][1] = -veloDeriv[0][0];

    //if( (abs(z-b) <= modelRad) || (abs(z+b) <= modelRad) ){
    //  velos[0] = 0.0;
    //  velos[1] = 0.0;
    //  velos[2] = 0.0;
    //}else{
    //  velos[0] = (veloDeriv[1][0] - veloDeriv[0][1]) * velTmp[1];
    //  velos[1] = (veloDeriv[0][1] - veloDeriv[1][0]) * velTmp[0];
    //  velos[2] = 0.0;
    //}
    //std::cout << veloDeriv << std::endl;
    return;
  }
  

}

