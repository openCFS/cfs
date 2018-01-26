#ifndef CAPA_IF_H
#define CAPA_IF_H

#include "def_use_unv.hh"

#define MAX_SETS 1000


#ifdef USE_ADOLC
//typedef adouble Double;
typedef adtl::adouble Double;
#else
typedef double Double;
#endif


struct SetInfo {
  Double time;  // absolute time
  int idx;      // index to data array
  int ndata;    // number of values per node/element
};

struct GDataInfo {
  int numtype55; // number of different 55 datasets
  int numtype56; // number of different 56 datasets
  int numtype58; // number of different 58 datasets
  int types55[MAX_SETS]; // integer coding of 55 data types
  int types56[MAX_SETS]; // integer coding of 56 data types
  int types58[MAX_SETS]; // integer coding of 58 data types
  int n55[MAX_SETS]; // contains for each different node-result type the number of sets
  int e56[MAX_SETS]; // contains for each different element-result type the number of sets
  SetInfo *Nsetinfo[MAX_SETS]; // stores for each type of node-result data set all its set information
  SetInfo *Esetinfo[MAX_SETS]; // stores for each type of element-result data set all its set information
  SetInfo *Hsetinfo[MAX_SETS];
};

class CapaInterfaceC {
public:
  enum NodeDataTypes {
    NoNodeData=0,
    Displacements, Velocities, Accelerations,
    ElectricPot, ElectricPotDeriv1, ElectricPotDeriv2,
    VelocityPot, VelocityPotDeriv1, VelocityPotDeriv2,
    MagVectorPot, MagVectorPotDeriv1, MagVectorPotDeriv2,
    ElectricSPot, ElectricSPotDeriv1, ElectricSPotDeriv2,
    MagScalarPot, MagScalarPotDeriv1, MagScalarPotDeriv2,
    ElectricCharge,
    NormalVelocity, SingleLayerPot, MultiLayerPot,
    Temperature
  };
  
  enum ElemDataTypes {
    NoElemData=0,
    Strains,
    Stresses,
    ElectricField,
    MagneticField,
    VelocityField,
    EddyCurrent
  };

  enum DataTypes {
    NoData=0,
    StaticData,
    TimeData,
    CwData,
    ComplexCwData,
    ModalData
  };


protected:
  long          nNodes;  // number of nodes
  long          nElems;  // number of elements
  long          nDim; // dimension of the mesh
  DataTypes     dataType; //type of data: transient, harmonic, etc.
  long          nNodeDataSets; //number of node data sets (set 55)
  Double        nodeXData[MAX_SETS]; //time stamp of the node data sets
  long          nodeDataSetsNumData[MAX_SETS]; 
  long          nElemDataSets; //number of element data sets (set 56)
  long          elemDataSetsNumData[MAX_SETS];
  Double        elemXData[MAX_SETS];
  long          numHistTypes;
  long          histDataSetsNumData[MAX_SETS];
  NodeDataTypes nodeDataSets[MAX_SETS]; //type of node data (displ., velocity, etc.); coded as integers
  ElemDataTypes elemDataSets[MAX_SETS];
  NodeDataTypes nodeHistTypes[MAX_SETS]; // UNIQUE types of node data (NOT per datatset)

  void ClipNodeIndex(long& iNode) const ;
  void ClipElemIndex(long& iElem) const ;
  void ClipDataIndex(long& iData, long nData) const ;
  void ClipNodeSet(long& iData) const ;
  void ClipElemSet(long& iData) const ;

public:
  CapaInterfaceC();
  ~CapaInterfaceC();

  int GetAnalysisType() const { return dataType; }
  
  long GetNumOfNodes() const { return nNodes; }
  long GetNumOfElements() const { return nElems; }
  long GetDimension() const { return nDim; }

  long GetNumOfNodeDataSets() const {return nNodeDataSets; }
  long GetNodeDataSetNumData(long set) const;
  const NodeDataTypes* GetNodeDataSets() const {return nodeDataSets; }
  
  long GetNumOfElemDataSets() const {return nElemDataSets; }
  const ElemDataTypes* GetElemDataSets() const {return elemDataSets; }
  long GetElemDataSetNumData(long set) const;
  
  int ReadUniversalfile(const char *fileName);
  const Double* GetNodeXData(void) const {return nodeXData;}  
  const Double* GetElemXData(void) const {return elemXData;}
  void GetHistXData(int &numSteps, Double *stepVals) const;
  long GetHistNumData(long type) const;
  
  void GetPos(long iNode, Double* pos) const;
  int GetMaxElemNodes(void) const;
  void GetElemNodes(long iElem, long& numNodes, long* nodes) const;
  void GetElemColor(long iElem, long& elemcolor) const;
  void GetElemType(long iElem, long& elemtype) const;
  int GetNodeData(long iNode, long iData, Double* d1, Double* d2) const;
  int GetNodeHistData(long iNode, int type, long step,
                      Double *d1, Double *d2) const;
  int GetElemData(long iElem, long iData, Double* d1, Double* d2) const;

  Double GetNodeDataTime(int idataset) const {return nodeXData[idataset];}

  void GetDataInfo(GDataInfo & datainfo) const;
};

#endif

/// Local Variables:
/// mode: C++
/// c-basic-offset: 2
/// End:
