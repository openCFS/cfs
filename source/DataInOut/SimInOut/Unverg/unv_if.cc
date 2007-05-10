#include <cstdio>
#include <string>
#include <cstdlib>
#include <iostream>

#include "unv_if.hh"
#include "unv_dat.hh"
//#include "bas_tool.h"

const char *nodeDataTypesStr[30] = {
  "NoNodeData",
  "Displacements", "Velocities", "Accelerations",
  "ElectricPot", "ElectricPotDeriv1", "ElectricPotDeriv2",
  "VelocityPot", "VelocityPotDeriv1", "VelocityPotDeriv2",
  "MagVectorPot", "MagVectorPotDeriv1", "MagVectorPotDeriv2",
  "ElectricSPot", "ElectricSPotDeriv1", "ElectricSPotDeriv2",
  "MagScalarPot", "MagScalarPotDeriv1", "MagScalarPotDeriv2",
  "ElectricCharge",
  "NormalVelocity, SingleLayerPot, MultiLayerPot",
  "Temperature"
};

const char *elemDataTypesStr[10] = {
  "NoElemData",
  "Strains",
  "Stresses",
  "ElectricField",
  "MagneticField",
  "VelocityField",
  "EddyCurrent"
};

const char* dataTypesStr[]= {
  "NoData",
  "TimeData",
  "CwData",
  "ComplexCwData",
  "ModalData"
};

CapaInterfaceC::CapaInterfaceC()
{
  nNodes = 0;
  nElems = 0;
  nNodeDataSets = 0;
  nElemDataSets = 0;
  int i;
  for (i=0; i<MAX_SETS; i++) {
    nodeDataSets[i] = NoNodeData;
    nodeDataSetsNumData[i] = 0;
  }
  for (i=0; i<MAX_SETS; i++) {
    elemDataSets[i] = NoElemData;
    elemDataSetsNumData[i] = 0;
  }
}

CapaInterfaceC::~CapaInterfaceC()
{
}

void CapaInterfaceC::ClipNodeIndex(long& iNode)
{
  if (iNode > nNodes)
    iNode = nNodes;
  if (iNode<1)
    iNode =1;
}

void CapaInterfaceC::ClipElemIndex(long& iElem)
{
  if (iElem > nElems)
    iElem = nElems;
  if (iElem<1)
    iElem =1;
}

void CapaInterfaceC::ClipDataIndex(long& iData, long nData) {
  if (iData > nData-1)
    iData = nData-1;
  if (iData<0)
    iData =0;
}

void CapaInterfaceC::ClipElemSet(long& iData) {
  iData = (iData > nElemDataSets-1) ? nElemDataSets-1 : iData;
  iData = (iData < 0) ? 0 : iData;
}

void CapaInterfaceC::ClipNodeSet(long& iData) {
  iData = (iData > nNodeDataSets-1) ? nNodeDataSets-1 : iData;
  iData = (iData < 0) ? 0 : iData;
}

long CapaInterfaceC::GetNodeDataSetNumData(long set) {
  ClipNodeSet(set);
  return nodeDataSetsNumData[set];
}

long CapaInterfaceC::GetElemDataSetNumData(long set) {
  ClipElemSet(set);
  return elemDataSetsNumData[set];
}

//
// liest Universal-File und fuellt Datenstruktur
// aktualisiert nNodes, nElem
//
int CapaInterfaceC::ReadUniversalfile(const char *fileName)
{
  //    if (!existFile(fileName))
  //        Fatal("File '%s' does not exist",fileName);

  char fname[1024];
  strncpy(fname,fileName,sizeof(fname));
  UNV_FILE=fname;

  Widget w=0;
  if(ReadUniversalFile(w) < 0)
    return -1;

  nNodes = N_NODES;
  nElems = N_ELEMENTS;
  nDim = DIM;
  nNodeDataSets = N_SETS55;

  int i;

    
  for (i=0; i<nNodeDataSets ; i++) {
    nodeDataSetsNumData[i] = SETS55[i].header.n_data_val_per_node;
    if (strstr(SETS55[i].header.id[0],"time"))
      dataType = TimeData;
    else if (strstr(SETS55[i].header.id[0],"mode no"))
      dataType = ModalData;
    else if (strstr(SETS55[i].header.id[0],"cw"))
      dataType = CwData;
    else if (strstr(SETS55[i].header.id[0],"complex"))
      dataType = ComplexCwData;

    nodeXData[i]=0.0;
    switch (dataType) {
      char* cp;
    case NoData:
      break;
    case TimeData:
      cp=strstr(SETS55[i].header.id[0],"time");
      if (cp) {
        cp+=4;
        nodeXData[i]=atof(cp);
      }
      break;
    case ModalData:
    case ComplexCwData:
    case CwData:
      cp=strstr(SETS55[i].header.id[0],"frequency");
      if (cp) {
        cp+=9;
        nodeXData[i]=atof(cp);
      }
    }

    char *setdata=SETS55[i].header.id[0];

    if (strstr(setdata,"displacement"))
      nodeDataSets[i]=Displacements;
    else if (strstr(setdata,"velocities"))
      nodeDataSets[i]=Velocities;
    else if (strstr(setdata,"accelerations"))
      nodeDataSets[i]=Accelerations;
    else if (strstr(setdata,"electric potential")) {
      if (strstr(setdata,"1st deriv"))
        nodeDataSets[i]=ElectricPotDeriv1;
      else if (strstr(setdata,"2nd deriv"))
        nodeDataSets[i]=ElectricPotDeriv2;
      else 
        nodeDataSets[i]=ElectricPot;
    }
    else if (strstr(setdata,"fluid potential")) {
      if (strstr(setdata,"1st deriv"))
        nodeDataSets[i]=VelocityPotDeriv1;
      else if (strstr(setdata,"2nd deriv"))
        nodeDataSets[i]=VelocityPotDeriv2;
      else if (strstr(setdata,"norm."))
        nodeDataSets[i]=NormalVelocity;
      else 
        nodeDataSets[i]=VelocityPot;
    }
    else if (strstr(setdata,"mag. vector potential")) {
      if (strstr(setdata,"1st deriv"))
        nodeDataSets[i]=MagVectorPotDeriv1;
      else if (strstr(setdata,"2nd deriv"))
        nodeDataSets[i]=MagVectorPotDeriv2;
      else 
        nodeDataSets[i]=MagVectorPot;
    }
    else if (strstr(setdata,"electric scalar potential")) {
      if (strstr(setdata,"1st deriv"))
        nodeDataSets[i]=ElectricSPotDeriv1;
      else if (strstr(setdata,"2nd deriv"))
        nodeDataSets[i]=ElectricSPotDeriv2;
      else 
        nodeDataSets[i]=ElectricSPot;
    }
    else if (strstr(setdata,"mag. scalar potential")) {
      if (strstr(setdata,"1st deriv"))
        nodeDataSets[i]=MagScalarPotDeriv1;
      else if (strstr(setdata,"2nd deriv"))
        nodeDataSets[i]=MagScalarPotDeriv2;
      else 
        nodeDataSets[i]=MagScalarPot;
    }
    else if (strstr(setdata,"electric charge"))
      nodeDataSets[i]=ElectricCharge;
    else if (strstr(setdata,"single layer potential"))
      nodeDataSets[i]=SingleLayerPot;
    else if (strstr(setdata,"double layer potential"))
      nodeDataSets[i]=MultiLayerPot;
    else if (strstr(setdata,"temperature"))
      nodeDataSets[i]=Temperature;
    else
      nodeDataSets[i] = NoNodeData;
        
    //      printf("iNodeData=%d, dataType=%s, type=%s, x=%e,\n\t header='%s'\n",
    //             i,nodeDataTypesStr[nodeDataSets[i]], dataTypesStr[dataType], nodeXData[i],
    //             SETS55[i].header.id[0]);
  }
  for (i=nNodeDataSets; i<MAX_SETS; i++)
    nodeDataSets[i] = NoNodeData;
    
  nElemDataSets= N_SETS56;

  for (i=0; i<nElemDataSets ; i++) {

    if (strstr(SETS56[i].header.id[0],"time"))
      dataType = TimeData;
    else if (strstr(SETS56[i].header.id[0],"mode no"))
      dataType = ModalData;
    else if (strstr(SETS56[i].header.id[0],"cw"))
      dataType = CwData;
    else if (strstr(SETS56[i].header.id[0],"complex"))
      dataType = ComplexCwData;

    elemDataSetsNumData[i] = SETS56[i].header.n_data_val_per_element_pos;
    elemXData[i]=0.0;
    switch (dataType) {
      char* cp;
    case NoData:
      break;
    case TimeData:
      cp=strstr(SETS56[i].header.id[0],"at time");
      if (cp) {
        cp+=7;
        elemXData[i]=atof(cp);
      }
      break;
    case ModalData:
    case ComplexCwData:
    case CwData:
      cp=strstr(SETS56[i].header.id[0],"frequency");
      if (cp) {
        cp+=9;
        elemXData[i]=atof(cp);
      }
    }

    char *setdata=SETS56[i].header.id[0];

    if (strstr(setdata,"stress"))
      elemDataSets[i]=Stresses;
    else if (strstr(setdata,"strain"))
      elemDataSets[i]=Strains;
    else if (strstr(setdata,"electric field"))
      elemDataSets[i]=ElectricField;
    else if (strstr(setdata,"velocity-gradient"))
      elemDataSets[i]=VelocityField;
    else if (strstr(setdata,"mag. flux density"))
      elemDataSets[i]=MagneticField;
    else if (strstr(setdata,"eddy current"))
      elemDataSets[i]=EddyCurrent;

    //      printf("iElemData=%d, dataType=%s, type=%s, x=%e,\n\t header='%s'\n",
    //             i,elemDataTypesStr[elemDataSets[i]], 
    //             dataTypesStr[dataType], elemXData[i], SETS56[i].header.id[0]);
  }

  for (i=nElemDataSets; i<MAX_SETS; i++)
    elemDataSets[i] = NoElemData;

  return 0;
}


 
//
// Knotengeometrie holen
//
void CapaInterfaceC::GetPos(long iNode, double* pos)
{
  ClipNodeIndex(iNode);
  pos[0]=NODES[iNode].x1;
  pos[1]=NODES[iNode].x2;
  pos[2]=NODES[iNode].x3;
}

//
// groesste Anzahl von Knoten in irgendeinem Element ermitteln
//
int CapaInterfaceC::GetMaxElemNodes(void)
{
  int maxNumNodes=0;
  for(int iElem=0; iElem<N_ELEMENTS; iElem++)
    if (maxNumNodes<ELEMENTS[iElem].n_nodes)
      maxNumNodes = ELEMENTS[iElem].n_nodes;
  return maxNumNodes;
}

//
// Element-Info holen
//
void CapaInterfaceC::GetElemNodes(long iElem, long& numNodes, long* nodes)
{
  ClipElemIndex(iElem);
  element_data *elem = &ELEMENTS[iElem];
  numNodes=elem->n_nodes;
  for(int i=0; i<numNodes; i++)
    nodes[i]=elem->points[i];
}

void CapaInterfaceC::GetElemColor(long iElem, long& elemcolor)
{
  ClipElemIndex(iElem);
  element_data *elem = &ELEMENTS[iElem];
  elemcolor=elem->color;
}

void CapaInterfaceC::GetElemType(long iElem, long& elemtype)
{
  ClipElemIndex(iElem);
  element_data *elem = &ELEMENTS[iElem];
  elemtype=elem->fe_type;
}

//
// Knotendaten holen
//
int CapaInterfaceC::GetNodeData(long iNode, long dataIndex,
                                double* data1, double* data2) {

  int iData;
  ClipNodeIndex(iNode);
  ClipNodeSet(dataIndex);

  long iDataSet = dataIndex, nData=nodeDataSetsNumData[iDataSet];

  switch (dataType) {
  case TimeData:
    for (iData=0; iData<nData; iData++) {
      data1[iData]=SETS55[iDataSet].dat[iNode].data[iData];
    }
    break;
  case ModalData:
  case CwData:
    for (iData=0; iData<nData; iData++)
      data1[iData]=SETS55[iDataSet].dat[iNode].data[iData];
    break;
  case ComplexCwData:
    for (iData=0; iData<nData; iData++) {
      data1[iData]=SETS55[iDataSet].dat[iNode].data[iData];
    }
  }
  return 0;
}

//
// Elementdaten holen
//
int CapaInterfaceC::GetElemData(long iElem, long dataIndex,
                                double* data1, double* data2) {

  int iData;
  ClipElemIndex(iElem);
  ClipElemSet(dataIndex);

  long iDataSet=dataIndex, nData=elemDataSetsNumData[iDataSet];

  switch (dataType) {
  case TimeData:
    for (iData=0; iData<nData; iData++) {
      data1[iData]=SETS56[iDataSet].dat[iElem].data[iData];
    }
    break;
  case ModalData:
  case CwData:
    for (iData=0; iData<nData; iData++)
      data1[iData]=SETS56[iDataSet].dat[iElem].data[iData];
    break;
  case ComplexCwData:
    for (iData=0; iData<nData; iData++) {
      data1[iData]=SETS56[iDataSet].dat[iElem].data[iData];
    }
  }
  return 0;
}


void CapaInterfaceC::GetDataInfo(GDataInfo & datainfo)
{
  int set,k;
  int ntypes=0;
  int doadd = 0;

  for( set=0; set<MAX_SETS; set++ ) {
    datainfo.types55[set]=0;
    datainfo.types56[set]=0;
    datainfo.n55[set]=0;
    datainfo.e56[set]=0;
  }


  //--------------------------Node Data--------------------------
  
  //get all different node data types
  for (set=0; set<nNodeDataSets; set++) {
    doadd = 1;
    for (k=0; k<=ntypes; k++)
      if (datainfo.types55[k]==nodeDataSets[set]) doadd=0;
          
    if (doadd==1) 
    {	    
      datainfo.types55[ntypes] = nodeDataSets[set];
      ntypes++;
    }
  }
  
  datainfo.numtype55=ntypes;

  //get for each node data type the number of sets
  int nset;
  for (k=0; k<ntypes; k++) {
    for (set=0; set<nNodeDataSets; set++)  {
      if (datainfo.types55[k]==nodeDataSets[set]) datainfo.n55[k]++;
    }
  }
  //get for each node set its info: time and index to solution
  for (k=0; k<ntypes; k++)
  {
    nset = 0;
    datainfo.Nsetinfo[k] = new SetInfo[datainfo.n55[k]];
    for (set=0; set<nNodeDataSets; set++) 
    {
      if (datainfo.types55[k]==nodeDataSets[set]) 
      {
        //std::cout << "node: unv_if: " << set << "; nset " << nset << "; " <<GetNodeDataSetNumData(set) << std::endl;
        datainfo.Nsetinfo[k][nset].time = GetNodeDataTime(set);
        datainfo.Nsetinfo[k][nset].idx  = set;
        datainfo.Nsetinfo[k][nset].ndata = GetNodeDataSetNumData(set);
        //std::cout << "node: *unv_if: " << set << "; nset " << datainfo.Nsetinfo[k][nset].idx << "; " << datainfo.Nsetinfo[k][nset].ndata << std::endl;
        nset++;
      }
    }
  }
  
  //--------------------------Element Data--------------------------
  
  int etypes=0;
  //get all different element data types
  for (set=0; set<nElemDataSets; set++) {
    doadd = 1;
    for (k=0; k<=etypes; k++)
      if (datainfo.types56[k]==elemDataSets[set]) doadd=0;
          
    if (doadd==1) 
    {	    
      datainfo.types56[etypes] = elemDataSets[set];
      etypes++;
    }
  }
  
  datainfo.numtype56=etypes;
  
  //get for each node data type the number of sets
  int eset;
  for (k=0; k<etypes; k++) {
    for (set=0; set<nElemDataSets; set++) {
      if (datainfo.types56[k]==elemDataSets[set]) datainfo.e56[k]++;
    }
  }
  //get for each element set its info: time and index to solution
  for (k=0; k<etypes; k++)
  {
    eset = 0;
    datainfo.Esetinfo[k] = new SetInfo[datainfo.e56[k]];
    for (set=0; set<nElemDataSets; set++) 
    {
      if (datainfo.types56[k]==elemDataSets[set]) 
      {
        //std::cout << "elem: unv_if: " << set << "; eset " << eset << "; " <<GetElemDataSetNumData(set) << std::endl;
        datainfo.Esetinfo[k][eset].time = GetNodeDataTime(set);
        datainfo.Esetinfo[k][eset].idx  = set;
        datainfo.Esetinfo[k][eset].ndata = GetElemDataSetNumData(set);
        //std::cout << "elem: *unv_if: " << set << "; eset " << datainfo.Esetinfo[k][eset].idx << "; " << datainfo.Esetinfo[k][eset].ndata << std::endl;
        eset++;
      }
    }
  }
          
}

/// Local Variables:
/// mode: C++
/// c-basic-offset: 2
/// End:
