#ifndef CAPA_IF_H
#define CAPA_IF_H

#define MAX_SETS 1000

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
      VelocityField
      };

  enum DataTypes {
      NoData,
      TimeData,
      CwData,
      ComplexCwData,
      ModalData
      };
    
protected:
  long          nNodes;
  long          nElems;
  DataTypes     dataType;
  long          nNodeDataSets;
  double        nodeXData[MAX_SETS];
  long          nodeDataSetsNumData[MAX_SETS];
  long          nElemDataSets;
  long          elemDataSetsNumData[MAX_SETS];
  double        elemXData[MAX_SETS];
  NodeDataTypes nodeDataSets[MAX_SETS];
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
  int GetNodeData(long iNode, long iData, double* d1, double* d2);
  int GetElemData(long iElem, long iData, double* d1, double* d2);

};

#endif



