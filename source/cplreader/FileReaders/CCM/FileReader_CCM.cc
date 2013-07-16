// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <string>
#include <iostream>
#include <fstream>
#include "stdio.h"
#include <iomanip>
#include <vector>
#include <algorithm>
#include <sstream>

#include "boost/tokenizer.hpp"

#include "Settings.hh"
#include "cplreader/FileReader.hh"
#include "FileReader_CCM.hh"
#include "Geometrics.hh"

#include "General/exception.hh"
#include "Settings.hh"
#include "General/environment.hh"
#include "integlib/elemIntegr.hh"

#include "FlowDataTypes.hh"
#include "OutputWriter.hh"
#include "CouplingHandler.hh"

#include "Meshes.hh"
#include "DataFile.hh"
#include "MeshFile.hh"
#include "MeshPreparator.hh"
#include "DualMesher.hh"
#include "MeshSplitter.hh"
#include "FlowDataLoader.hh"


using namespace CCM;

namespace CoupledField
{

	FileReader_CCM::FileReader_CCM(const std::string& name,
                         const UInt dim,
                         const UInt numFiles,
                         const UInt startIndex):
                         FileReader(name, dim, numFiles, startIndex)
  {
  }

  FileReader_CCM::~FileReader_CCM()
  {
  }


  void FileReader_CCM::Init() {
    std::cout << "Initialzing dual mesh Star-CCM+ File Reader" << std::endl;
    Settings& settings = Settings::Instance();
    firstStep_ = settings.GetInt("firststep");
    stepOffset_ = settings.GetInt("stepoffset");
    digits_ = settings.GetInt("digits");
    if (digits_ < 0) {
      digits_ = 0;
    }
    velocityName = settings.GetString("VelocityName");
    std::cout << "Name of Velocity field: " << velocityName << std::endl;
    pressureName = settings.GetString("PressureName");
    std::cout << "Name of Pressure field: " << pressureName << std::endl;
    divLHTName = settings.GetString("DivLHTName");
    std::cout << "Name of Divergence of Lighthill Tensor field: " << divLHTName << std::endl;
    laplacePName = settings.GetString("LaplacePName");
    std::cout << "Name of Laplace of Pressure field: " << laplacePName << std::endl;
    std::cout << "Reading Steps: [ " << (firstStep_ + stepOffset_) << " , " << (firstStep_ + stepOffset_ + numSteps_ - 1) <<  " ] " << std::endl;
    std::cout << "Saving as Steps: [ " << firstStep_ << " , " << (firstStep_ + numSteps_ - 1) <<  " ] " << std::endl;
    
    // opening mesh file  
    std::stringstream ss;
	  ss << baseName_ << "mesh.ccm";
	  std::string inFileName = ss.str();
    MeshFile meshFile(inFileName);
    
    // preparing dualization
    MeshPreparator preparator;
    DualizableMesh cfdMesh;
    std::vector<ConsecutiveMap> valueConMaps;
    preparator.PrepareMesh(&meshFile, cfdMesh, valueConMaps);
    meshFile.Close();
    //cfdMesh.Clear();
    //WaitReturn();
    
    // dualizing mesh
    DualizedMesh dualMesh;
    DualMesher* dualizer = new DualMesher();
    dualizer->DualizeMesh(cfdMesh, dualMesh, true);
    cfdMesh.Clear();
    
    // splitting mesh
    MeshSplitter* splitter = new MeshSplitter();
    splitter->SplitMesh(dualMesh, cplmesh_);
    dualMesh.Clear();
    
    // applying results to members
    maxNumElemNodes_ = cplmesh_.maxNumElemNodes;
    numRegions_ = cplmesh_.numRegions;
    numElemsPerRegion_.swap(cplmesh_.numElemsPerRegion);
    numNodesPerRegion_.swap(cplmesh_.numNodesPerRegion);
    dim_ = 3;
    
    //WaitReturn();
    
    loader = new FlowDataLoader(dualMesh.vertexCount, valueConMaps);
    std::cout << "Initialization of dual mesh Star-CCM+ File Reader finished" << std::endl;
  }

  //! get node coordinates from the corresponding file
  void FileReader_CCM::ReadNodalCoords(std::vector<Double> & NODECOORD) {
    NODECOORD.swap(cplmesh_.NODECOORD);
  }


  //! get topology information from the corresponding topology file
  void FileReader_CCM::ReadTopology(std::vector<UInt> & TOPOLOGYDATA,
                               std::vector<UInt> & elemTypes) {
    TOPOLOGYDATA.swap(cplmesh_.TOPOLOGYDATA);
    elemTypes.swap(cplmesh_.elemTypes);
  }

  //! get nodal values from the corresponding fluid datafile the new way
  void FileReader_CCM::ReadNodalValues(std::vector<FlowDataType>& nodalFlowData,
                               const std::vector<bool>& activeParts,
                               const UInt timeStepIdx) {
    std::stringstream thisFileName;
    uint actualStep = timeStepIdx + firstStep_ + stepOffset_;
    std::stringstream ss;
    ss << actualStep;
    std::string numberstr = ss.str();
    int zeros = (int) digits_ - (int) numberstr.length();
    if (zeros < 0) {
      zeros = 0;
    }
    thisFileName << baseName_;
    for (int i = 0; i < zeros; i++) {
      thisFileName << "0";
    }
    thisFileName << (actualStep) << ".ccm";
    
    DataFile dataFile(thisFileName.str(), false);
    FlowDataType& fdt = nodalFlowData[0];
    if (requiredResults_[FLUIDMECH_VELOCITY]) {
      FlowDataPartStruct* fdps = &fdt[FLUIDMECH_VELOCITY];
      loader->LoadVector(fdps, dataFile, velocityName, FLUIDMECH_VELOCITY);
    }
    if (requiredResults_[FLUIDMECH_PRESSURE]) {
      FlowDataPartStruct* fdps = &fdt[FLUIDMECH_PRESSURE];
      loader->LoadScalar(fdps, dataFile, pressureName, FLUIDMECH_PRESSURE);
    }
    if (requiredResults_[FLUIDMECH_PRESSURE_DERIV_2]) {
      FlowDataPartStruct* fdps = &fdt[FLUIDMECH_PRESSURE_DERIV_2];
      loader->LoadScalar(fdps, dataFile, laplacePName, FLUIDMECH_PRESSURE_DERIV_2);
    }
    if (requiredResults_[FLUIDMECH_DIV_LH_T]) {
      FlowDataPartStruct* fdps = &fdt[FLUIDMECH_DIV_LH_T];
      loader->LoadVector(fdps, dataFile, divLHTName, FLUIDMECH_DIV_LH_T);
    }
    dataFile.Close();
  }

  void FileReader_CCM::GetRegionElements(std::vector<UInt> & regionElements,
                                 const UInt regionIdx) {
    regionElements.swap(cplmesh_.regionElements);
  }

} // end of namespace
