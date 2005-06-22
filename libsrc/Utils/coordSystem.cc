#include "coordSystem.hh"

#include "DataInOut/ParamHandling/BaseParamHandler.hh"

namespace CoupledField{


  CoordSystem::CoordSystem(const std::string & name,
                           Grid * ptGrid) {
    ENTER_FCN("CoordSystem::CoordSystem");

    name_ = name;
    ptGrid_ = ptGrid;

    dim_ = ptGrid_->GetDim();
    
  }

  CoordSystem::~CoordSystem(){
    ENTER_FCN("CoordSystem::~CoordSystem");
  }

  void CoordSystem::GetPoint(Vector<Double> & vec,
                             StdVector<std::string> & keyVec,
                             StdVector<std::string> & attrVec,
                             StdVector<std::string> & valVec 
                             ) {
    
    ENTER_FCN( "CoordSystem::GetVector");
    
   
    std::string nodeName;
    StdVector<std::string> mKeyVec, coordNames;
    StdVector<UInt> nodes;
    coordNames = "x", "y", "z";

    // check, if node is given by name and eventually get it
    // from the grid object
    attrVec.Push_back(std::string());
    valVec.Push_back(std::string());    
    mKeyVec = keyVec;
    mKeyVec.Push_back("node");

    params->Get(mKeyVec, attrVec, valVec, nodeName);
    std::cerr << "name of the node was " << nodeName << std::endl;

    vec.Resize(dim_);
    if ( nodeName != "none" ) {
      ptGrid_->GetNodesByName(nodes,nodeName);

      // check if more than one node is defined by this name
      if ( nodes.GetSize() != 1 ) {
        (*error) << "CoordinateSystem: There are more than 1 nodes defined "
                 << "defined by '" << nodeName << "'.\nTherefore it is "
                 << "impossible to determine ONE coordinate location. Please "
                 << "ensure, that only ONE node is defined by this name!";
        Error( __FILE__, __LINE__ );
      }

      ptGrid_->GetNodeCoordinate(vec, nodes[0]);

    } else {

      // if no node name was given, read in global x,y and z coordinate
      for (UInt i=0; i<vec.GetSize(); i++) {
        mKeyVec = keyVec;
        mKeyVec.Push_back(coordNames[i]);    
        params->Get(mKeyVec, attrVec, valVec, vec[i]);
      }
    }
      
    std::cerr << "vec = \n" << vec << std::endl;
  }


} // end of namespace
