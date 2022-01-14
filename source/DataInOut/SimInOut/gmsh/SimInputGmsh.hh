// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_GMSHREADER_2009
#define FILE_GMSHREADER_2009

#include <map>

// #include <boost/bimap/set_of.hpp>
#include <boost/bimap.hpp>

#include <DataInOut/SimInput.hh>

namespace CoupledField {

  class CoordSystem;
  
  /**
   **/
  class SimInputGmsh : public SimInput
  {
  private:

  public:
    SimInputGmsh(std::string fileName, PtrParamNode inputNode,
                 PtrParamNode infoNode);
    virtual ~SimInputGmsh();

    virtual void InitModule();

    virtual void ReadMesh(Grid *mi);
    
    //! Return dimension of the mesh
    virtual UInt GetDim();
    
    virtual void GetNumMultiSequenceSteps( std::map<UInt, BasePDE::AnalysisType>& analysis,
                                           std::map<UInt, UInt>& numSteps,
                                           bool isHistory );
  
  private:

    void ElemType2FEType(UInt eType, Elem::FEType& feType, UInt& numNodes);
    void InitElemNodeMap();
    void GetPhysEnt2RegionMapFromXML();
    void LinearizeElem(Elem::FEType* elemType);
    void TransformNodes(CoordSystem& coordSys, double scaleFac);

    struct NodeStruct 
    {
      UInt nodeId;
      Double x,y,z;
    };

    StdVector<Double> nodeCoords_;
    std::map<UInt, UInt> nodeNumMap_;
    std::vector< Elem::FEType > elementTypes_;
    std::vector< UInt > elementPhysicsTypes_;
    std::vector< UInt > connectivity_;
    UInt numNodesPerElem_;
    UInt dim_;
    std::map<UInt, std::string> physEntities2RegionNames_;
    std::map<UInt, std::string> physEntities2NamedNodes_;
    std::map<UInt, std::string> physEntities2NamedElems_;
    std::map<std::string, bool> linearizeRegions_;
    bool readOnlySomeRegions_;
    std::string coordSysId_;
    Double scaleFac_;
  }; 

} 

#endif
