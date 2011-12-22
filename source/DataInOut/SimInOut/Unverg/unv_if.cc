#include <cstdio>
#include <cstdlib>
#include <cstring>

#include <string>
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
  for ( int i=0; i<MAX_SETS; ++i ) {
    nodeDataSets[i] = NoNodeData;
    nodeDataSetsNumData[i] = 0;
    elemDataSets[i] = NoElemData;
    elemDataSetsNumData[i] = 0;
  }
}

CapaInterfaceC::~CapaInterfaceC()
{
}

void CapaInterfaceC::ClipNodeIndex(long& iNode) const {
  if (iNode > nNodes)
    iNode = nNodes;
  if (iNode<1)
    iNode =1;
}

void CapaInterfaceC::ClipElemIndex(long& iElem) const {
  if (iElem > nElems)
    iElem = nElems;
  if (iElem<1)
    iElem =1;
}

void CapaInterfaceC::ClipDataIndex(long& iData, long nData) const {
  if (iData > nData-1)
    iData = nData-1;
  if (iData<0)
    iData =0;
}

void CapaInterfaceC::ClipElemSet(long& iData) const {
  iData = (iData > nElemDataSets-1) ? nElemDataSets-1 : iData;
  iData = (iData < 0) ? 0 : iData;
}

void CapaInterfaceC::ClipNodeSet(long& iData) const {
  iData = (iData > nNodeDataSets-1) ? nNodeDataSets-1 : iData;
  iData = (iData < 0) ? 0 : iData;
}

long CapaInterfaceC::GetNodeDataSetNumData(long set) const {
  ClipNodeSet(set);
  return nodeDataSetsNumData[set];
}

long CapaInterfaceC::GetElemDataSetNumData(long set) const {
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

    
  for (i=0; i<nNodeDataSets ; ++i) {
    nodeDataSetsNumData[i] = SETS55[i].header.n_data_val_per_node;
    // dataType can be determined from header
    dataType=NoData;
    switch (SETS55[i].header.analysis_type) {
      case 1: // static
        dataType = StaticData;
        break;
      case 2: // "normal" mode ?
      case 5: // frequency response
        if (SETS55[i].header.data_type == 5) // complex data
          dataType = ComplexCwData;
        else if (std::strstr(SETS55[i].header.id[0],"mode no"))
          dataType = ModalData;
        else
          dataType = CwData;
        break;
      case 4: // transient
        dataType = TimeData;
        break;
      default: // fall back to CAPA algorithm in any other case
        if (std::strstr(SETS55[i].header.id[0],"time"))
          dataType = TimeData;
        else if (std::strstr(SETS55[i].header.id[0],"mode no"))
          dataType = ModalData;
        else if (std::strstr(SETS55[i].header.id[0],"cw") || 
                 std::strstr(SETS55[i].header.id[0],"_AM") ||
                 std::strstr(SETS55[i].header.id[0],"_PH") ||
                 std::strstr(SETS55[i].header.id[0],"_RE") ||
                 std::strstr(SETS55[i].header.id[0],"_IM"))
          dataType = CwData;
        else if (std::strstr(SETS55[i].header.id[0],"complex"))
          dataType = ComplexCwData;
        break;
    }

    nodeXData[i]=0.0;
    // the time/frequency value is stored in 1st type specific real parameter
    if (dataType != NoData)
      nodeXData[i] = SETS55[i].header.type_spec_double_par[0];
    /*switch (dataType) {
      char* cp;
    case NoData:
      break;
    case TimeData:
      cp=std::strstr(SETS55[i].header.id[0],"time");
      if (cp) {
        cp+=4;
        nodeXData[i]=atof(cp);
      }
      break;
    case ModalData:
    case ComplexCwData:
    case CwData:
      cp=std::strstr(SETS55[i].header.id[0],"frequency");
      if (cp) {
        cp+=9;
        nodeXData[i]=atof(cp);
      }
    }*/

    char *setdata=SETS55[i].header.id[0];
    switch (SETS55[i].header.spec_data_type) {
      case  2: // stress
      case  3: // strain
      case  4: // element force
      case  6: // heat flux
      case  7: // strain energy
      case  9: // reaction force
      case 10: // kinetic energy
      case 13: // strain energy density
      case 14: // kinetic energy density
      case 15: // hydro-static pressure
      case 16: // heat gradient
      case 18: // coefficient of pressure
        nodeDataSets[i]=NoNodeData;
        break;
      case 5: // temperature
        nodeDataSets[i]=Temperature;
        break;
      case 8: // displacement
        nodeDataSets[i]=Displacements;
        break;
      case 11: // velocity
        nodeDataSets[i]=Velocities;
        break;
      case 12: // acceleration
        nodeDataSets[i]=Accelerations;
        break;
      default: // default to CAPA descriptions
        if (std::strstr(setdata,"displacement"))
          nodeDataSets[i]=Displacements;
        else if (std::strstr(setdata,"velocities"))
          nodeDataSets[i]=Velocities;
        else if (std::strstr(setdata,"accelerations"))
          nodeDataSets[i]=Accelerations;
        else if (std::strstr(setdata,"electric potential")) {
          if (std::strstr(setdata,"1st deriv"))
            nodeDataSets[i]=ElectricPotDeriv1;
          else if (std::strstr(setdata,"2nd deriv"))
            nodeDataSets[i]=ElectricPotDeriv2;
          else 
            nodeDataSets[i]=ElectricPot;
        }
        else if (std::strstr(setdata,"fluid potential")) {
          if (std::strstr(setdata,"1st deriv"))
            nodeDataSets[i]=VelocityPotDeriv1;
          else if (std::strstr(setdata,"2nd deriv"))
            nodeDataSets[i]=VelocityPotDeriv2;
          else if (std::strstr(setdata,"norm."))
            nodeDataSets[i]=NormalVelocity;
          else 
            nodeDataSets[i]=VelocityPot;
        }
        else if (std::strstr(setdata,"APRES")) {
          nodeDataSets[i]=VelocityPot;
        }
        else if (std::strstr(setdata,"mag. vector potential")) {
          if (std::strstr(setdata,"1st deriv"))
            nodeDataSets[i]=MagVectorPotDeriv1;
          else if (std::strstr(setdata,"2nd deriv"))
            nodeDataSets[i]=MagVectorPotDeriv2;
          else 
            nodeDataSets[i]=MagVectorPot;
        }
        else if (std::strstr(setdata,"electric scalar potential")) {
          if (std::strstr(setdata,"1st deriv"))
            nodeDataSets[i]=ElectricSPotDeriv1;
          else if (std::strstr(setdata,"2nd deriv"))
            nodeDataSets[i]=ElectricSPotDeriv2;
          else 
            nodeDataSets[i]=ElectricSPot;
        }
        else if (std::strstr(setdata,"mag. scalar potential")) {
          if (std::strstr(setdata,"1st deriv"))
            nodeDataSets[i]=MagScalarPotDeriv1;
          else if (std::strstr(setdata,"2nd deriv"))
            nodeDataSets[i]=MagScalarPotDeriv2;
          else 
            nodeDataSets[i]=MagScalarPot;
        }
        else if (std::strstr(setdata,"electric charge"))
          nodeDataSets[i]=ElectricCharge;
        else if (std::strstr(setdata,"single layer potential"))
          nodeDataSets[i]=SingleLayerPot;
        else if (std::strstr(setdata,"double layer potential"))
          nodeDataSets[i]=MultiLayerPot;
        else if (std::strstr(setdata,"temperature"))
          nodeDataSets[i]=Temperature;
        else
          nodeDataSets[i] = NoNodeData;
        break;
    }
        
    //      printf("iNodeData=%d, dataType=%s, type=%s, x=%e,\n\t header='%s'\n",
    //             i,nodeDataTypesStr[nodeDataSets[i]], dataTypesStr[dataType], nodeXData[i],
    //             SETS55[i].header.id[0]);
  }
  for (i=nNodeDataSets; i<MAX_SETS; ++i)
    nodeDataSets[i] = NoNodeData;
    
  nElemDataSets= N_SETS56;

  for (i=0; i<nElemDataSets ; i++) {

    if (std::strstr(SETS56[i].header.id[0],"time"))
      dataType = TimeData;
    else if (std::strstr(SETS56[i].header.id[0],"mode no"))
      dataType = ModalData;
    else if (std::strstr(SETS56[i].header.id[0],"cw"))
      dataType = CwData;
    else if (std::strstr(SETS56[i].header.id[0],"complex"))
      dataType = ComplexCwData;
    else
      dataType = StaticData;

    elemDataSetsNumData[i] = SETS56[i].header.n_data_val_per_element_pos;
    elemXData[i]=0.0;
    switch (dataType) {
      char* cp;
    case NoData:
    case StaticData:
      break;
    case TimeData:
      cp=std::strstr(SETS56[i].header.id[0],"at time");
      if (cp) {
        cp+=7;
        elemXData[i]=atof(cp);
      }
      break;
    case ModalData:
    case ComplexCwData:
    case CwData:
      cp=std::strstr(SETS56[i].header.id[0],"frequency");
      if (cp) {
        cp+=9;
        elemXData[i]=atof(cp);
      }
    }

    char *setdata=SETS56[i].header.id[0];

    if (std::strstr(setdata,"stress"))
      elemDataSets[i]=Stresses;
    else if (std::strstr(setdata,"strain"))
      elemDataSets[i]=Strains;
    else if (std::strstr(setdata,"electric field"))
      elemDataSets[i]=ElectricField;
    else if (std::strstr(setdata,"velocity-gradient"))
      elemDataSets[i]=VelocityField;
    else if (std::strstr(setdata,"mag. flux density"))
      elemDataSets[i]=MagneticField;
    else if (std::strstr(setdata,"eddy current"))
      elemDataSets[i]=EddyCurrent;

    //      printf("iElemData=%d, dataType=%s, type=%s, x=%e,\n\t header='%s'\n",
    //             i,elemDataTypesStr[elemDataSets[i]], 
    //             dataTypesStr[dataType], elemXData[i], SETS56[i].header.id[0]);
  }

  for (i=nElemDataSets; i<MAX_SETS; i++)
    elemDataSets[i] = NoElemData;

  if(NODES58 != NULL) 
  {
    numHistTypes = 0;
    for ( i=0; i<N_NODES; ++i ) {
      for ( int j=0; j<NODES58[i].num_datasets; ++j ) {
      
        switch (NODES58[i].sets[j].header.fnc_type) {
        case 1:  // time response
          dataType = TimeData;
          break;
        case 2:  // auto spectrum
        case 3:  // cross spectrum
        case 4:  // frequency response
        case 12: // spectrum
        case 24: // shock response spectrum
          if (NODES58[i].sets[j].header.data_type > 4)
            dataType = ComplexCwData;
          else
            dataType = CwData;
          break;
        case 22: // Eigenvalue
        case 23: // Eigenvector
          dataType = ModalData;
          break;
        default: // try to determine from x-axis
          switch (NODES58[i].sets[j].header.axis_data[0].data_type) {
          case 17: // time
            dataType = TimeData;
            break;
          case 18: // frequency
            if (NODES58[i].sets[j].header.data_type > 4)
              dataType = ComplexCwData;
            else
              dataType = CwData;
            break;
          }
          break;
        }
      
        switch (NODES58[i].sets[j].header.axis_data[1].data_type) {
        case 5: // temperature
          NODES58[i].sets[j].capaDataType = Temperature;
          break;
        case 8: // displacement
          NODES58[i].sets[j].capaDataType = Displacements;
          break;
        case 11: // velocity
          NODES58[i].sets[j].capaDataType = Velocities;
          break;
        case 12: // acceleration
          NODES58[i].sets[j].capaDataType = Accelerations;
          break;
        default:
          NODES58[i].sets[j].capaDataType = NoNodeData;
          break;
        }
      
        int k;
        for ( k=0; k<numHistTypes; ++k) {
          if (nodeHistTypes[k] == NODES58[i].sets[j].capaDataType)
            break;
        }
        if ( k == numHistTypes ) {
          nodeHistTypes[k] = (NodeDataTypes) NODES58[i].sets[j].capaDataType;
          histDataSetsNumData[k]
            = (NODES58[i].sets[j].header.res_dir == 0 ? 1 : 3);
          ++numHistTypes;
        }
      }
    }
  }
  

  return 0;
}

void CapaInterfaceC::GetHistXData(int &numSteps, double *stepVals) const {
  int i, j;

  if(NODES58 == NULL) 
  {
    numSteps = 0;
    return;
  }
  
  for ( i=0; i<N_NODES; ++i ) {
    if (NODES58[i].num_datasets > 0) {
      if ((numSteps <= 0)
          || (numSteps > NODES58[i].sets[0].header.num_values))
        numSteps = NODES58[i].sets[0].header.num_values;
      if (stepVals != NULL) {
        for ( j=0; j<numSteps; ++j ) {
          stepVals[j] = NODES58[i].sets[0].data[j].step_val;
        }
      }
      // We just take the steps of the first dataset, that we find.
      // If the other datasets differ from that, we are in trouble.
      break;
    }
  }
}

long CapaInterfaceC::GetHistNumData(long type) const {
  if ((type < 0) || (type >= numHistTypes))
    return 0;
  
  return histDataSetsNumData[type];
}

//
// Knotengeometrie holen
//
void CapaInterfaceC::GetPos(long iNode, double* pos) const {
  ClipNodeIndex(iNode);
  pos[0]=NODES[iNode].x1;
  pos[1]=NODES[iNode].x2;
  pos[2]=NODES[iNode].x3;
}

//
// groesste Anzahl von Knoten in irgendeinem Element ermitteln
//
int CapaInterfaceC::GetMaxElemNodes(void) const {
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
const {
  ClipElemIndex(iElem);
  element_data *elem = &ELEMENTS[iElem];
  numNodes=elem->n_nodes;
  for(int i=0; i<numNodes; i++)
    nodes[i]=elem->points[i];
}

void CapaInterfaceC::GetElemColor(long iElem, long& elemcolor) const {
  ClipElemIndex(iElem);
  element_data *elem = &ELEMENTS[iElem];
  elemcolor=elem->color;
}

void CapaInterfaceC::GetElemType(long iElem, long& elemtype) const {
  ClipElemIndex(iElem);
  element_data *elem = &ELEMENTS[iElem];
  elemtype=elem->fe_type;
}

//
// Knotendaten holen
//
int CapaInterfaceC::GetNodeData(long iNode, long dataIndex,
                                double* data1, double* data2) const {

  int iData;
  ClipNodeIndex(iNode);
  ClipNodeSet(dataIndex);

  long iDataSet = dataIndex, nData=nodeDataSetsNumData[iDataSet];

  switch (dataType) {
    case StaticData:
    case TimeData:
    case ModalData:
    case CwData:
      for (iData=0; iData<nData; ++iData)
        data1[iData]=SETS55[iDataSet].dat[iNode].data[iData];
      break;
    case ComplexCwData:
      for (iData=0; iData<nData; ++iData) {
        data1[iData]=SETS55[iDataSet].dat[iNode].data[2*iData];
        data2[iData]=SETS55[iDataSet].dat[iNode].data[2*iData+1];
      }
      break;
    case NoData:
      break;
  }
  return 0;
}

int CapaInterfaceC::GetNodeHistData(long iNode, int type, long step,
                                    double *data1, double *data2) const {
  int i;
  
  ClipNodeIndex(iNode);
  iNode -= 1;
  
  for ( i=0; i<6; ++i ) {
    data1[i] = 0.0;
    data2[i] = 0.0;
  }

  for ( i=0; i<NODES58[iNode].num_datasets; ++i ) {
    if ( NODES58[iNode].sets[i].capaDataType == type) {
      
      if ((step < 1) || (step > NODES58[iNode].sets[i].header.num_values))
        return 1;
      
      int iDof = labs(NODES58[iNode].sets[i].header.res_dir);
      if (iDof != 0)
        iDof -= 1;
      
      double sign = (NODES58[iNode].sets[i].header.res_dir < 0 ? -1.0 : 1.0);
      
      data1[iDof] = sign * NODES58[iNode].sets[i].data[step-1].real;
      
      if (NODES58[iNode].sets[i].header.data_type >= 5)
        data2[iDof] = sign * NODES58[iNode].sets[i].data[step-1].imag;
    }
  }
  
  return i == 0 ? 1 : 0;
}

//
// Elementdaten holen
//
int CapaInterfaceC::GetElemData(long iElem, long dataIndex,
                                double* data1, double* data2) const {

  int iData;
  ClipElemIndex(iElem);
  ClipElemSet(dataIndex);

  long iDataSet=dataIndex, nData=elemDataSetsNumData[iDataSet];

  switch (dataType) {
    case StaticData:
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
      break;
    case NoData:
      break;
  }
  return 0;
}


void CapaInterfaceC::GetDataInfo(GDataInfo & datainfo) const {
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

  //=============================================================
  // node histories (dataset 58)
  //=============================================================
    
  datainfo.numtype58 = numHistTypes;
  for ( k=0; k<numHistTypes; ++k ) {
    datainfo.types58[k] = nodeHistTypes[k];
  }
  
}

/// Local Variables:
/// mode: C++
/// c-basic-offset: 2
/// End:
