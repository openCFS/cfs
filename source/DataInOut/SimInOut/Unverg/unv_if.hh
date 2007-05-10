#ifndef CAPA_IF_H
#define CAPA_IF_H

#include <def_use_unv.hh>

#define MAX_SETS 1000

struct SetInfo {
  double time;  //absolute time
  int idx;      //index to data array
  int ndata;    //number of values per node/element
};

struct GDataInfo {
  int numtype55; //number of different 55 data sets
  int numtype56; //number of different 56 data sets
  int types55[MAX_SETS]; // integer coding of 55 data types
  int types56[MAX_SETS]; //integer coding of 56 data types
  int n55[MAX_SETS];  //contains for each different node-result type the number of sets
  int e56[MAX_SETS]; //contains for each different element-result type the number of sets
  SetInfo *Nsetinfo[MAX_SETS]; //stores for each type of node-result data set all its set information
  SetInfo *Esetinfo[MAX_SETS]; //stores for each type of element-result data set all its set information
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
    NoElemData,
    Strains,
    Stresses,
    ElectricField,
    MagneticField,
    VelocityField,
    EddyCurrent
  };

  enum DataTypes {
    NoData,
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
  double        nodeXData[MAX_SETS]; //time stemp of the node data sets
  long          nodeDataSetsNumData[MAX_SETS]; 
  long          nElemDataSets; //number of element data sets (set 56)
  long          elemDataSetsNumData[MAX_SETS];
  double        elemXData[MAX_SETS];
  NodeDataTypes nodeDataSets[MAX_SETS]; //type of node data (displ., velocity, etc.); coded as integers
  ElemDataTypes elemDataSets[MAX_SETS];

  void ClipNodeIndex(long& iNode);
  void ClipElemIndex(long& iElem);
  void ClipDataIndex(long& iData, long nData);
  void ClipNodeSet(long& iData);
  void ClipElemSet(long& iData);

public:
  CapaInterfaceC();
  ~CapaInterfaceC();

    
  long GetNumOfNodes() const { return nNodes; }
  long GetNumOfElements() const { return nElems; }
  long GetDimension() const { return nDim; }

  long GetNumOfNodeDataSets() const {return nNodeDataSets; }
  long GetNodeDataSetNumData(long set);
  const NodeDataTypes* GetNodeDataSets() const {return nodeDataSets; }
  
  long GetNumOfElemDataSets() const {return nElemDataSets; }
  const ElemDataTypes* GetElemDataSets() const {return elemDataSets; }
  long GetElemDataSetNumData(long set);
  
  int ReadUniversalfile(const char *fileName);
  const double* GetNodeXData(void) {return nodeXData;}  
  const double* GetElemXData(void) {return elemXData;}  
  void GetPos(long iNode, double* pos);
  int GetMaxElemNodes(void);
  void GetElemNodes(long iElem, long& numNodes, long* nodes);
  void GetElemColor(long iElem, long& elemcolor);
  void GetElemType(long iElem, long& elemtype);
  int GetNodeData(long iNode, long iData, double* d1, double* d2);
  int GetElemData(long iElem, long iData, double* d1, double* d2);

  double GetNodeDataTime(int idataset) {return nodeXData[idataset];}

  void GetDataInfo(GDataInfo & datainfo);
};

#endif

/// Local Variables:
/// mode: C++
/// c-basic-offset: 2
/// End:
