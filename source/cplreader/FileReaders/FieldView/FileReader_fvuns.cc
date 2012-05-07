// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "Settings.hh"

#include "FileReader_fvuns.hh"
#include "readfvuns.h"
#include "readfvuns.cpp"
#include "fv_reader_tags.h"

//Filling Start

#include "stdlib.h" // For exit()
#include <vector>
#include <map>
#include <iostream>
#include <iomanip>
#include <string>





//Manipulation Start
#include <fstream>
#include "stdio.h"
#include <iomanip>
#include <sstream>
#include "sys/stat.h"


#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/filesystem/convenience.hpp"
#include "boost/filesystem/exception.hpp"
namespace fs=boost::filesystem;

#include "Domain/resultInfo.hh"

#include "cplreader/Settings.hh"

#include "def_cplreader.hh"
#include "cplreader/FileReader.hh"

//Manipulation End





namespace CoupledField
{

FileReader_fvuns::FileReader_fvuns(const std::string& name,
    const UInt dim,
    const UInt numFiles,
    const UInt startIndex) :
    FileReader(name, dim, numFiles, startIndex)
    {
  std::cout << "FileReader_fvuns " << name << " " << dim << " " << numFiles << std::endl;
    }

FileReader_fvuns::~FileReader_fvuns()
{
}

void FileReader_fvuns::Init()
{
  Settings& settings = Settings::Instance();
  std::stringstream fullFileName;
  fullFileName <<  settings.GetString("basedir") << "/" << name_ << "/" << name_ << ".fvuns";
  file_name=fullFileName.str().c_str();

  Initreadfvuns();
  numRegions_=1;   //alle Teilgebiete werden zu einer Region zusammengefasst
  dim_=3;
  numNodesPerRegion_.resize(numRegions_);

  int allnodes=0;
  int allhex=0;
  int alltet=0;
  int allpri=0;
  int allpyr=0;
  std::cout<<"number of grids: "<<num_grids<<std::endl;
  for(int i=0; i<num_grids; i++)
  {
    allnodes+=parts[i].nnodes;
    allhex+=parts[i].nhex;
    alltet+=parts[i].ntet;
    allpri+=parts[i].npri;
    allpyr+=parts[i].npyr;
  }

  numNodesPerRegion_[0]=(UInt) allnodes;



  numElemsPerRegion_.resize(numRegions_);

  numElemsPerRegion_[0]=(UInt) (allhex+alltet+allpri+allpyr); //nur eine globale Region

  if (settings.GetInt("numsteps"))
    numSteps_ = (UInt) settings.GetInt("numsteps");
  { 
    if(allhex!=0)
      maxNumElemNodes_=8;
    else
    {
      if(allpri!=0)
        maxNumElemNodes_=6;
      else
      {
        if(allpyr!=0)
          maxNumElemNodes_=5;
        else
        {
          if(alltet!=0)
            maxNumElemNodes_=4;
        }
      }
    }

  }
}

static std::vector<unsigned int> IDcount;  
void FileReader_fvuns::ReadNodalCoords(std::vector<Double> & NODECOORD)
{
  std::cout<<"~~~~~~~~~~~~~~~~ReadNodalCoords~~~~~~~~~~~~~~~"<<std::endl;


  int allnodes=0;
  int nodecount=0;
  for(int i=0; i<num_grids; i++)
    allnodes+=parts[i].nnodes;
  NODECOORD.resize(allnodes*3);
  int r;
  for(int gw=0; gw<num_grids;gw++)
  {
    r=0; 
    for(int i=nodecount; i<parts[gw].nnodes+nodecount; i++)   //globale Knotenliste wird aus lokalen Knotenlisten erstellt
    {
      NODECOORD[i*3+0]=parts[gw].xyz[parts[gw].nnodes*0+r];
      NODECOORD[i*3+1]=parts[gw].xyz[parts[gw].nnodes*1+r];
      NODECOORD[i*3+2]=parts[gw].xyz[parts[gw].nnodes*2+r];
      r++;
    }
    std::cout<<NODECOORD[15*3+0]<<" "<<NODECOORD[15*3+1]<<" "<<NODECOORD[15*3+2]<<std::endl;
    std::cout<<NODECOORD[17*3+0]<<" "<<NODECOORD[17*3+1]<<" "<<NODECOORD[17*3+2]<<std::endl;
    std::cout<<NODECOORD[42*3+0]<<" "<<NODECOORD[42*3+1]<<" "<<NODECOORD[42*3+2]<<std::endl;
    std::cout<<NODECOORD[18*3+0]<<" "<<NODECOORD[18*3+1]<<" "<<NODECOORD[18*3+2]<<std::endl;



    nodecount+=parts[gw].nnodes;
  }
}

void FileReader_fvuns::ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
    std::vector<UInt> & elemTypes)
{ 
  std::cout<<"~~~~~~~~~~~~~~~~ReadTopology~~~~~~~~~~~~~~~"<<std::endl;
  std::vector <int> Topologydata;
  int nodecount=0;
  int nhex_=0;
  int npyr_=0;
  int npri_=0;
  int ntet_=0;
  // int VertsPerCell; // TODO: Unused variable VertsPerCell
  // int Cells=0; // TODO: Unused variable Cells


  int allCells=0;
  for(int i=0; i<num_grids;i++)
  {
    allCells+=parts[i].nhex+parts[i].ntet+parts[i].npyr+parts[i].npri;
  }


  for(int i=0; i<num_grids;i++)
  {
    nhex_+=parts[i].nhex;
    ntet_+=parts[i].ntet;
    npyr_+=parts[i].npyr;
    npri_+=parts[i].npri;
    // if(nhex_!=0) 
    // {
      // VertsPerCell=8;
      // Cells=parts[i].nhex;
    // }
    // if(ntet_!=0) 
    // {
      // VertsPerCell=4;
      // Cells=parts[i].ntet;
    // }
    // if(npyr_!=0) 
    // {
      // VertsPerCell=5;
      // Cells=parts[i].npyr;
    // }
    // if(npri_!=0) 
    // {
      // VertsPerCell=6;
      // Cells=parts[i].npri;
    // }

    if(parts[i].nhex!=0)
    {
      for(int j=0; j<8*parts[i].nhex; j++)
        Topologydata.push_back(parts[i].hex_ids[j]+nodecount+1); //+1, weil VertID bei 0 beginnt!
      for(int j=0; j< parts[i].nhex; j++)
      {
        elemTypes.push_back(Elem::HEXA8);

      }
    }
    if(parts[i].ntet!=0)
    {
      for(int j=0; j<4*parts[i].ntet; j++)
      { if(j<4)
        std::cout<<"Nodes: "<< parts[i].tet_ids[j]+nodecount<<std::endl;
      Topologydata.push_back(parts[i].tet_ids[j]+nodecount+1);} //+1, weil VertID bei 0 beginnt!
      for(int j=0; j< parts[i].ntet; j++)
        elemTypes.push_back(Elem::TET4);
      std::cout<<" tets"<<parts[i].ntet<<std::endl;

    }
    if(parts[i].npyr!=0)
    {
      for(int j=0; j<5*parts[i].npyr; j++)
        Topologydata.push_back(parts[i].pyr_ids[j]+nodecount+1); //+1, weil VertID bei 0 beginnt!
      for(int j=0; j< parts[i].npyr; j++)
        elemTypes.push_back(Elem::PYRA5);

    }
    if(parts[i].npri!=0)
    {
      for(int j=0; j<6*parts[i].npri; j++)
        Topologydata.push_back(parts[i].pri_ids[j]+nodecount+1); //+1, weil VertID bei 0 beginnt!
      for(int j=0; j< parts[i].npri; j++)
        elemTypes.push_back(Elem::WEDGE6);
      std::cout<<" pris"<<parts[i].npri<<std::endl;
    }

    nodecount+=parts[i].nnodes;
  }


  TOPOLOGYDATA.resize(elemTypes.size()*maxNumElemNodes_,0);  //Elementliste wird umsortiert von Fieldview- auf cplReaderelementformat
  int topowalker=0, bigtopowalker=0;
  std::cout<<"size Elementypes: "<<elemTypes.size()<<std::endl;

  for(unsigned int k=0; k<elemTypes.size(); k++)
  {

    if(elemTypes[k]==Elem::HEXA8)
    {
      TOPOLOGYDATA[0+bigtopowalker]=Topologydata[0+topowalker];
      TOPOLOGYDATA[1+bigtopowalker]=Topologydata[1+topowalker];
      TOPOLOGYDATA[2+bigtopowalker]=Topologydata[3+topowalker];
      TOPOLOGYDATA[3+bigtopowalker]=Topologydata[2+topowalker];
      TOPOLOGYDATA[4+bigtopowalker]=Topologydata[4+topowalker];
      TOPOLOGYDATA[5+bigtopowalker]=Topologydata[5+topowalker];
      TOPOLOGYDATA[6+bigtopowalker]=Topologydata[7+topowalker];
      TOPOLOGYDATA[7+bigtopowalker]=Topologydata[6+topowalker];
      topowalker+=8;
    }
    if(elemTypes[k]==Elem::TET4)
    {
      TOPOLOGYDATA[0+bigtopowalker]=Topologydata[1+topowalker];
      TOPOLOGYDATA[1+bigtopowalker]=Topologydata[2+topowalker];
      TOPOLOGYDATA[2+bigtopowalker]=Topologydata[3+topowalker];
      TOPOLOGYDATA[3+bigtopowalker]=Topologydata[0+topowalker];
      topowalker+=4;
    }

    if(elemTypes[k]==Elem::PYRA5)
    {
      TOPOLOGYDATA[0+bigtopowalker]=Topologydata[0+topowalker];
      TOPOLOGYDATA[1+bigtopowalker]=Topologydata[1+topowalker];
      TOPOLOGYDATA[2+bigtopowalker]=Topologydata[2+topowalker];
      TOPOLOGYDATA[3+bigtopowalker]=Topologydata[3+topowalker];
      TOPOLOGYDATA[4+bigtopowalker]=Topologydata[4+topowalker];
      topowalker+=5;
    }

    if(elemTypes[k]==Elem::WEDGE6)
    {
      TOPOLOGYDATA[0+bigtopowalker]=Topologydata[0+topowalker];
      TOPOLOGYDATA[1+bigtopowalker]=Topologydata[3+topowalker];
      TOPOLOGYDATA[2+bigtopowalker]=Topologydata[5+topowalker];
      TOPOLOGYDATA[3+bigtopowalker]=Topologydata[1+topowalker];
      TOPOLOGYDATA[4+bigtopowalker]=Topologydata[2+topowalker];
      TOPOLOGYDATA[5+bigtopowalker]=Topologydata[4+topowalker];
      topowalker+=6;

    }
    bigtopowalker+=maxNumElemNodes_;

  }



}

void FileReader_fvuns::GetRegionElements(std::vector<UInt> & regionElements,
    const UInt regionIdx)
{

  UInt elemOffset = 0;

  for(UInt i=0; i < regionIdx; i++)
    elemOffset += numElemsPerRegion_[i];

  regionElements.resize(numElemsPerRegion_[regionIdx]);
  for (UInt wuz=0; wuz<numElemsPerRegion_[regionIdx];wuz++)
    regionElements[wuz]=elemOffset+wuz+1;



}

std::string FileReader_fvuns::GetRegionName(const UInt regionIdx)
{
  std::vector <std::string>RegionName;
  RegionName.push_back("innerMesh");  //nur eine globale Region

  return RegionName[regionIdx];

}

void FileReader_fvuns::GetNodeGroups(std::map<std::string,
    std::vector<UInt> >& nodeGroups) 
{
}

void FileReader_fvuns::ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
    const std::vector<bool>& activeParts,
    const UInt timeStepIdx)
{

  std::cout<<"~~~~~~~~~~~~~~~~ReadNodalValues~~~~~~~~~~~~~~~"<<std::endl;


  Settings& settings = Settings::Instance();


  cleanup();      //alte Werte entfernen                                              
  std::stringstream fullFileName;
  fullFileName <<  settings.GetString("basedir") << "/" << name_ << "/" << name_ << timeStepIdx<<".fvuns";  //nächster Zeitschritt, auf Dateikonvention achten
  //name_timeStepIdx, z.B. FFS2D_0, FFS2D_1, FFS2D_2...
  file_name=fullFileName.str().c_str();

  Initreadfvuns();  //neue Werte beziehen




  std::vector <Double> NodalValues;
  int allnodes=0;
  int nodecount=0;
  for(int i=0; i<num_grids; i++)
    allnodes+=parts[i].nnodes;

  NodalValues.resize(allnodes*4);

  for(int i=0;i<num_grids;i++)
  {
    int r=0;
    for(int j=nodecount;j<parts[i].nnodes+nodecount;j++)  //aus lokalen Wertelisten globale Werteliste erstellen
    {
      NodalValues[j]=parts[i].vars[r];
      NodalValues[j+allnodes]=parts[i].vars[r+parts[i].nnodes];
      NodalValues[j+(allnodes*2)]=parts[i].vars[r+(parts[i].nnodes*2)];
      NodalValues[j+(allnodes*3)]=parts[i].vars[r+(parts[i].nnodes*3)];
      r++;

    }
    nodecount+=parts[i].nnodes;
  }
  int nandetektor=0;
  for(unsigned int i=0; i< NodalValues.size(); i++)   //NAN / INF ausschliessen
    if(isinf(NodalValues[i]) || isnan(NodalValues[i]))
    {
      NodalValues[i]=0;
      nandetektor+=1;
    }

  if(nandetektor!=0)
    std::cout<<nandetektor<<"-times NAN or INF detected in fluidMechVelocity!"<<std::endl;



  int dim=3;       
  UInt numDOFs=0;
  UInt actRegion=0;
  if(timeStepIdx==0)
  {
    FlowDataType& fd = nodalFlowData[0];//wie in OpenFoam wird nur die erste Region betrachtet

    FlowDataPartStruct* fdps = NULL;


    fdps = &fd[ACOU_RHS_LOAD];
    fdps->isActive = true;
    fdps->definedOn = ResultInfo::NODE; // nodes
    fdps->entryType = ResultInfo::SCALAR;
    fdps->dofNames.push_back("-");
    fdps->unit = MapSolTypeToUnit(ACOU_RHS_LOAD);
    fdps->resultName = SolutionTypeEnum.ToString(ACOU_RHS_LOAD);

    fdps = &fd[FLUIDMECH_PRESSURE];
    fdps->isActive = true;
    fdps->definedOn = ResultInfo::NODE; // nodes
    fdps->entryType = ResultInfo::SCALAR;
    fdps->dofNames.push_back("-");
    fdps->unit = MapSolTypeToUnit(FLUIDMECH_PRESSURE);
    fdps->resultName = SolutionTypeEnum.ToString(FLUIDMECH_PRESSURE);


    fdps = &fd[FLUIDMECH_VELOCITY];
    fdps->isActive = true;
    if(fdps->dofNames.empty())
    {    
      fdps->definedOn = ResultInfo::NODE; // nodes
      fdps->entryType = ResultInfo::VECTOR;

      fdps->dofNames.push_back("x");
      fdps->dofNames.push_back("y");
      if(dim==3)
        fdps->dofNames.push_back("z");
      fdps->unit = MapSolTypeToUnit(FLUIDMECH_VELOCITY);
      fdps->resultName = SolutionTypeEnum.ToString(FLUIDMECH_VELOCITY);
    }     
    numDOFs = fdps->dofNames.size();






    fdps->data.resize(numDOFs*numNodesPerRegion_[actRegion]); 
  }      
  solutionTypes_.push_back(ACOU_RHS_LOAD);
  solutionTypes_.push_back(FLUIDMECH_VELOCITY);
  solutionTypes_.push_back(FLUIDMECH_VELOCITY);
  solutionTypes_.push_back(FLUIDMECH_VELOCITY);
  solutionTypes_.push_back(FLUIDMECH_PRESSURE);

  dofIndices_.push_back(0);
  dofIndices_.push_back(0);
  dofIndices_.push_back(1);
  dofIndices_.push_back(2);
  dofIndices_.push_back(0);


  int NNODES=allnodes;




  for(UInt j=0; j<solutionTypes_.size(); j++)
  {
    if(j==0)
    {
      for(int i=0; i < NNODES; i++)

      {
        numDOFs = nodalFlowData[actRegion][solutionTypes_[j]].dofNames.size();
        std::vector<Double>& data = nodalFlowData[actRegion][solutionTypes_[j]].data;
        data.resize(numDOFs*numNodesPerRegion_[actRegion]);

        /*fdps->*/data[i*numDOFs+dofIndices_[j]] = 0;
      }
    }
    if(j==1)
    {//std::cout<<"Velo0"<<std::endl;
      for(int i=0; i < NNODES; i++)
      {
        numDOFs = nodalFlowData[actRegion][solutionTypes_[j]].dofNames.size();
        std::vector<Double>& data = nodalFlowData[actRegion][solutionTypes_[j]].data;
        data.resize(numDOFs*numNodesPerRegion_[actRegion]);

        /*fdps->*/data[i*numDOFs+dofIndices_[j]] = NodalValues[i+NNODES*0];
      }
    }
    if(j==2)
    {//std::cout<<"Velo1"<<std::endl;
      for(int i=0; i < NNODES; i++)
      {
        numDOFs = nodalFlowData[actRegion][solutionTypes_[j]].dofNames.size();
        std::vector<Double>& data = nodalFlowData[actRegion][solutionTypes_[j]].data;
        data.resize(numDOFs*numNodesPerRegion_[actRegion]);

        /*fdps->*/data[i*numDOFs+dofIndices_[j]] = NodalValues[i+NNODES*1];
      }
    }
    if(j==3)
    {//std::cout<<"Velo2"<<std::endl;
      for(int i=0; i < NNODES; i++)
      {
        numDOFs = nodalFlowData[actRegion][solutionTypes_[j]].dofNames.size();
        std::vector<Double>& data = nodalFlowData[actRegion][solutionTypes_[j]].data;
        data.resize(numDOFs*numNodesPerRegion_[actRegion]);

        /*fdps->*/data[i*numDOFs+dofIndices_[j]] = NodalValues[i+NNODES*2];
      }
    }
    if(j==4)
    {//std::cout<<"Pressure"<<std::endl;
      for(int i=0; i < NNODES; i++)
      {
        numDOFs = nodalFlowData[actRegion][solutionTypes_[j]].dofNames.size();
        std::vector<Double>& data = nodalFlowData[actRegion][solutionTypes_[j]].data;
        data.resize(numDOFs*numNodesPerRegion_[actRegion]);

        /*fdps->*/data[i*numDOFs+dofIndices_[j]] = NodalValues[i+NNODES*3];
      }
    }
  }


}
}



