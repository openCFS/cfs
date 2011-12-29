<<<<<<< .working
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/filesystem/convenience.hpp"
#include "boost/filesystem/exception.hpp"
namespace fs=boost::filesystem;

#include "def_use_mpcci.hh"

#if (MpCCI_RELEASE == 305)
#include "mpcci.h"
#endif

#include "General/exception.hh"
#include "Settings.hh"
#include "General/environment.hh"
#include "DataInOut/SimInOut/hdf5/hdf5io.hh"
#include "cplreader/CouplingHandler.hh"

#include "OutputWriter_MpCCI.hh"

namespace CoupledField
{

  OutputWriter_MpCCI::OutputWriter_MpCCI() {

  }

  OutputWriter_MpCCI::~OutputWriter_MpCCI() {
    Settings& settings = Settings::Instance();
    if(settings.GetString("coupling") == "MpCCI") {
      CCI_Finalize();
    }
  }

  void OutputWriter_MpCCI::Init(int argc, char** argv,
                               const std::string& outputPath)
  {
    if(settings.GetString("coupling") == "MpCCI") {
      quantityId1_ = 7;
      quantityDim1_ = 1;

      quantityId2_ = 8;
      quantityDim2_ = 3;

      //Define specific coupling description for cplreader side
      meshId_    = 1;
      GlobalDim_ = 3;

      CCI_Init_with_id_string(&argc, &argv, "simulationcode1");
    }
  }

  void OutputWriter_MpCCI::WriteNodalCoords(const std::vector<Double>& coords) {
    // write the dimension of the grid.
    H5IO::WriteAttribute( meshGroup_, "Dimension", dim_ );

    // ================
    //  Node Locations
    // ================
    UInt nNodes = coords.size()/3;
    H5::Group nodeGroup;
    try {
      nodeGroup = meshGroup_.createGroup( "Nodes" );
    } H5_CATCH( "Could not create node group" );

    H5IO::WriteAttribute( nodeGroup, "NumNodes", nNodes );

    H5IO::Write2DArray( nodeGroup, "Coordinates", nNodes,
                        3, &coords[0], dPropList_ );

    nodeGroup.close();
  }

  void OutputWriter_MpCCI::WriteTopology(UInt maxNumElemNodes,
                                        const std::vector<UInt>& elemTypes,
                                        const std::vector<UInt>& elems)
  {
  }

  void OutputWriter_MpCCI::WriteRegion(std::string regionName,
                                      const std::vector<UInt>& regionElems,
                                      const std::vector<UInt>& regionNodes,
                                      UInt regionDim)
  {
    if(settings.GetString("coupling") == "MpCCI") {
      // MpCCI part
      CCI_Def_partition(meshId_, actRegion+1);
      std::cout << "Define nodes" << std::endl;
      
      /*
        Uint numNodes = ptFileReader_->GetNumNodes(actRegion);
        std::vector<Double> mpcciCoords;
        mpcciCoords.resize(dim*numNodes);
        for(Uint i=0; i<numNodes; i++) {
        for(Uint d=0; d<dim; d++) {
        mpcciCoords[i*dim+d] = nodalCoords_[actRegion][i*3+d];
        }
        }
      */
      
      //define the nodes
      CCI_Def_nodes(meshId_, actRegion+1, 3,
                    (int)ptFileReader_->GetNumNodes(actRegion), 0, NULL,
                    CCI_DOUBLE, &nodalCoords_[actRegion][0]);
      
      std::cout << "Define elements" << std::endl;
      
      // Convert CFS++ element types to MpCCI element types
      std::vector<int> mpcciElemTypes;
      std::vector<int>::iterator etIt, etEnd;
      mpcciElemTypes.resize(elemTypes_[actRegion].size());
      std::copy(elemTypes_[actRegion].begin(),
                elemTypes_[actRegion].end(),
                mpcciElemTypes.begin());
      etIt = mpcciElemTypes.begin();
      etEnd = mpcciElemTypes.end();
      for( ; etIt != etEnd; etIt++)
      {
        *etIt = ElemTypes2MpCCI((FEType) *etIt);
      }
      
      //define the elements
      CCI_Def_elems(meshId_, actRegion+1, nElems, 0, NULL,
                    nElems, (int*) &mpcciElemTypes[0], (int*) &numNodesPerElem_[actRegion][0],
                    (int*) &topology_[actRegion][0]);
    }
  }

  void OutputWriter_MpCCI::WriteNodeGroups(const std::map< std::string, std::vector<UInt> >&
                                          nodeGroups)
  {
    if(settings.GetString("coupling") == "MpCCI") {
      //Close the definition phase; contact detection.
      //i take part in the coupling
      std::cout << "Closing MpCCI Setup" << std::endl;
      MpCCIprocess_ = 1;
      CCI_Close_setup(MpCCIprocess_);
    }
  }

  void OutputWriter_MpCCI::WriteElemGroups(const std::map< std::string, std::vector<UInt> >&
                                          elemGroups,
                                          const std::map< std::string, std::vector<UInt> >&
                                          elemGroupNodes,
                                          const std::map< std::string, UInt >&
                                          groupDims)
  {
  }

  void OutputWriter_MpCCI::WriteFlowData(CouplingHandler* cplHandler,
                                         UInt actRegion,
                                         const std::vector<std::string>& outputFields)
  {
    // Send fields for current partition to MpCCI
    if(settings.GetString("coupling") == "MpCCI") {
      UInt numNodesPart = ptFileReader_->GetNumNodes(actRegion);
      UInt numDOFs;
      static std::vector<Double> acouRhsField;
      static std::vector<Double> velField;
      
      // acouSrc scalar
      acouRhsField.resize(numNodesPart * quantityDim1_);
      
      fdIt = flowData[actRegion].find(ACOU_RHS_LOAD);
      
      // Copy data from partition struct to a dummy vector.
      if(fdIt != flowData[actRegion].end())
      {
        FlowDataPartStruct& fdps = flowData[actRegion][ACOU_RHS_LOAD];
        std::copy(fdps.data.begin(), fdps.data.end(),
                  acouRhsField.begin());
      }
      
      CCI_Put_nodes( meshId_, actRegion+1, quantityId1_,
                     quantityDim1_, numNodesPart, 0, NULL,
                     CCI_DOUBLE, &acouRhsField[0]);
      
      // velocity u, v, w
      velField.resize(numNodesPart * quantityDim2_);
      
      fdIt = flowData[actRegion].find(FLUIDMECH_VELOCITY);
      
      // Copy data from partition struct to a dummy vector.
      if(fdIt != flowData[actRegion].end())
      {
        FlowDataPartStruct& fdps = flowData[actRegion][FLUIDMECH_VELOCITY];
        numDOFs = fdps.dofNames.size();
        
        for(UInt i=0; i<numNodesPart; i++)
        {
          UInt idx1 = numNodesPart * numDOFs;
          UInt idx2 = numNodesPart * quantityDim2_;
          
          for(UInt j=0; j<numDOFs; j++)
            velField[idx2 + j] = fdps.data[idx1 + j];
        }
      }
      
      CCI_Put_nodes( meshId_, actRegion+1, quantityId2_,
                     quantityDim2_, numNodesPart, 0, NULL,
                     CCI_DOUBLE, &velField[0]);
    }
  }
  
  void OutputWriter_MpCCI::WriteUserData(const std::map<std::string, std::string>& userData)
  {
  }

  void OutputWriter_MpCCI::BeginResults(const std::vector<FlowDataType>* flowData)
  {
    // MpCCI status variables
    int globalConvergence = CCI_CONTINUE;
    int myConvergence     = CCI_CONTINUE;
  }
  
  void OutputWriter_MpCCI::BeginStep(UInt stepNum, Double stepValue) {
    OutputWriter::BeginStep(stepNum, stepValue);
    Settings& settings = Settings::Instance();
    UInt externalFiles = settings.GetInt("extfiles");

    // Check, if step group is already open and create it if not
    if( currMeshStepGroup_.getId() <= 0 ) {
      std::stringstream stepName;
      stepName << "/Results/Mesh/MultiStep_1/Step_" << stepNum;
      
      // Create new step group.
      try {
        currMeshStepGroup_= mainGroup_.createGroup(stepName.str());
        H5IO::WriteAttribute( currMeshStepGroup_, "StepValue", actStepValue_ );
        
        if(externalFiles)
          CreateExternalFile(stepNum);
      } H5_CATCH( "Can not create dataset for step " << actStepNum_ );
    }
  }

  void OutputWriter_MpCCI::EndStep() {
    OutputWriter::EndStep();

    if(settings.GetString("coupling") == "MpCCI") {
      //Send values via MpCCI
      int nQuantityIds=0;
      int quantityIds[] = {0, 0};
      int nLocalMeshIds=1;
      int localMeshIds[] = {1};
      int remoteCodeId = 2;
      int comm = CCI_COMM_RCODE[remoteCodeId];
      
      quantityIds[nQuantityIds] = quantityId1_;
      nQuantityIds++;
      quantityIds[nQuantityIds] = quantityId2_;
      nQuantityIds++;
      
      CCI_Send(nQuantityIds, quantityIds, nLocalMeshIds, localMeshIds, comm);
      CCI_Check_convergence(myConvergence,&globalConvergence,CCI_ANY_CODE);
      std::cout<<"CCI_CONTINUE_VALUE"<<CCI_CONTINUE<<std::endl;
    }
  }

  void OutputWriter_MpCCI::EndResults()
  {
  }

  int OutputWriter_MpCCI::ElemTypes2MpCCI(FEType et)
  {
    std::string elemTypeName;
    switch(et)
    {
    case ET_LINE2:
      return CCI_ELEM_LINE;
      break;

    case ET_TRIA3:
      return CCI_ELEM_TRIANGLE;
      break;

    case ET_QUAD4:
      return CCI_ELEM_QUAD;
      break;

    case ET_QUAD8:
      return CCI_ELEM_QUAD8;
      break;

    case ET_TET4:
      return CCI_ELEM_TETRAHEDRON;
      break;

    case ET_WEDGE6:
      return CCI_ELEM_PRISM;
      break;

    case ET_HEXA8:
      return CCI_ELEM_HEXAHEDRON;
      break;

    case ET_PYRA5:
      return CCI_ELEM_PYRAMID;
      break;

    default:
      Enum2String(et, elemTypeName);
      EXCEPTION("Element type " << elemTypeName
                << " cannot be used with MpCCI.");
    }
    return 0;
  }

}
=======
// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/filesystem/operations.hpp"
#include "boost/filesystem/path.hpp"
#include "boost/filesystem/convenience.hpp"
#include "boost/filesystem/exception.hpp"
namespace fs=boost::filesystem;

#include "def_use_mpcci.hh"

#if (MpCCI_RELEASE == 305)
#include "mpcci.h"
#endif

#include "General/exception.hh"
#include "Settings.hh"
#include "General/environment.hh"
#include "DataInOut/SimInOut/hdf5/hdf5io.hh"
#include "cplreader/CouplingHandler.hh"

#include "OutputWriter_MpCCI.hh"

namespace CoupledField
{

  OutputWriter_MpCCI::OutputWriter_MpCCI() {

  }

  OutputWriter_MpCCI::~OutputWriter_MpCCI() {
    Settings& settings = Settings::Instance();
    if(settings.GetString("coupling") == "MpCCI") {
      CCI_Finalize();
    }
  }

  void OutputWriter_MpCCI::Init(int argc, char** argv,
                               const std::string& outputPath)
  {
    if(settings.GetString("coupling") == "MpCCI") {
      quantityId1_ = 7;
      quantityDim1_ = 1;

      quantityId2_ = 8;
      quantityDim2_ = 3;

      //Define specific coupling description for cplreader side
      meshId_    = 1;
      GlobalDim_ = 3;

      CCI_Init_with_id_string(&argc, &argv, "simulationcode1");
    }
  }

  void OutputWriter_MpCCI::WriteNodalCoords(const std::vector<Double>& coords) {
    // write the dimension of the grid.
    H5IO::WriteAttribute( meshGroup_, "Dimension", dim_ );

    // ================
    //  Node Locations
    // ================
    UInt nNodes = coords.size()/3;
    H5::Group nodeGroup;
    try {
      nodeGroup = meshGroup_.createGroup( "Nodes" );
    } H5_CATCH( "Could not create node group" );

    H5IO::WriteAttribute( nodeGroup, "NumNodes", nNodes );

    H5IO::Write2DArray( nodeGroup, "Coordinates", nNodes,
                        3, &coords[0], dPropList_ );

    nodeGroup.close();
  }

  void OutputWriter_MpCCI::WriteTopology(UInt maxNumElemNodes,
                                        const std::vector<UInt>& elemTypes,
                                        const std::vector<UInt>& elems)
  {
  }

  void OutputWriter_MpCCI::WriteRegion(std::string regionName,
                                      const std::vector<UInt>& regionElems,
                                      const std::vector<UInt>& regionNodes,
                                      UInt regionDim)
  {
    if(settings.GetString("coupling") == "MpCCI") {
      // MpCCI part
      CCI_Def_partition(meshId_, actRegion+1);
      std::cout << "Define nodes" << std::endl;
      
      /*
        Uint numNodes = ptFileReader_->GetNumNodes(actRegion);
        std::vector<Double> mpcciCoords;
        mpcciCoords.resize(dim*numNodes);
        for(Uint i=0; i<numNodes; i++) {
        for(Uint d=0; d<dim; d++) {
        mpcciCoords[i*dim+d] = nodalCoords_[actRegion][i*3+d];
        }
        }
      */
      
      //define the nodes
      CCI_Def_nodes(meshId_, actRegion+1, 3,
                    (int)ptFileReader_->GetNumNodes(actRegion), 0, NULL,
                    CCI_DOUBLE, &nodalCoords_[actRegion][0]);
      
      std::cout << "Define elements" << std::endl;
      
      // Convert CFS++ element types to MpCCI element types
      std::vector<int> mpcciElemTypes;
      std::vector<int>::iterator etIt, etEnd;
      mpcciElemTypes.resize(elemTypes_[actRegion].size());
      std::copy(elemTypes_[actRegion].begin(),
                elemTypes_[actRegion].end(),
                mpcciElemTypes.begin());
      etIt = mpcciElemTypes.begin();
      etEnd = mpcciElemTypes.end();
      for( ; etIt != etEnd; etIt++)
      {
        *etIt = ElemTypes2MpCCI((FEType) *etIt);
      }
      
      //define the elements
      CCI_Def_elems(meshId_, actRegion+1, nElems, 0, NULL,
                    nElems, (int*) &mpcciElemTypes[0], (int*) &numNodesPerElem_[actRegion][0],
                    (int*) &topology_[actRegion][0]);
    }
  }

  void OutputWriter_MpCCI::WriteNodeGroups(const std::map< std::string, std::vector<UInt> >&
                                          nodeGroups)
  {
    if(settings.GetString("coupling") == "MpCCI") {
      //Close the definition phase; contact detection.
      //i take part in the coupling
      std::cout << "Closing MpCCI Setup" << std::endl;
      MpCCIprocess_ = 1;
      CCI_Close_setup(MpCCIprocess_);
    }
  }

  void OutputWriter_MpCCI::WriteElemGroups(const std::map< std::string, std::vector<UInt> >&
                                          elemGroups,
                                          const std::map< std::string, std::vector<UInt> >&
                                          elemGroupNodes,
                                          const std::map< std::string, UInt >&
                                          groupDims)
  {
  }

  void OutputWriter_MpCCI::WriteFlowData(CouplingHandler* cplHandler,
                                         UInt actRegion,
                                         const std::vector<std::string>& outputFields)
  {
    // Send fields for current partition to MpCCI
    if(settings.GetString("coupling") == "MpCCI") {
      UInt numNodesPart = ptFileReader_->GetNumNodes(actRegion);
      UInt numDOFs;
      static std::vector<Double> acouRhsField;
      static std::vector<Double> velField;
      
      // acouSrc scalar
      acouRhsField.resize(numNodesPart * quantityDim1_);
      
      fdIt = flowData[actRegion].find(ACOU_RHS_LOAD);
      
      // Copy data from partition struct to a dummy vector.
      if(fdIt != flowData[actRegion].end())
      {
        FlowDataPartStruct& fdps = flowData[actRegion][ACOU_RHS_LOAD];
        std::copy(fdps.data.begin(), fdps.data.end(),
                  acouRhsField.begin());
      }
      
      CCI_Put_nodes( meshId_, actRegion+1, quantityId1_,
                     quantityDim1_, numNodesPart, 0, NULL,
                     CCI_DOUBLE, &acouRhsField[0]);
      
      // velocity u, v, w
      velField.resize(numNodesPart * quantityDim2_);
      
      fdIt = flowData[actRegion].find(FLUIDMECH_VELOCITY);
      
      // Copy data from partition struct to a dummy vector.
      if(fdIt != flowData[actRegion].end())
      {
        FlowDataPartStruct& fdps = flowData[actRegion][FLUIDMECH_VELOCITY];
        numDOFs = fdps.dofNames.size();
        
        for(UInt i=0; i<numNodesPart; i++)
        {
          UInt idx1 = numNodesPart * numDOFs;
          UInt idx2 = numNodesPart * quantityDim2_;
          
          for(UInt j=0; j<numDOFs; j++)
            velField[idx2 + j] = fdps.data[idx1 + j];
        }
      }
      
      CCI_Put_nodes( meshId_, actRegion+1, quantityId2_,
                     quantityDim2_, numNodesPart, 0, NULL,
                     CCI_DOUBLE, &velField[0]);
    }
  }
  
  void OutputWriter_MpCCI::WriteUserData(const std::map<std::string, std::string>& userData)
  {
  }

  void OutputWriter_MpCCI::BeginResults(const std::vector<FlowDataType>* flowData)
  {
    // MpCCI status variables
    int globalConvergence = CCI_CONTINUE;
    int myConvergence     = CCI_CONTINUE;
  }
  
  void OutputWriter_MpCCI::BeginStep(UInt stepNum, Double stepValue) {
    OutputWriter::BeginStep(stepNum, stepValue);
    Settings& settings = Settings::Instance();
    UInt externalFiles = settings.GetInt("extfiles");

    // Check, if step group is already open and create it if not
    if( currMeshStepGroup_.getId() <= 0 ) {
      std::stringstream stepName;
      stepName << "/Results/Mesh/MultiStep_1/Step_" << stepNum;
      
      // Create new step group.
      try {
        currMeshStepGroup_= mainGroup_.createGroup(stepName.str());
        H5IO::WriteAttribute( currMeshStepGroup_, "StepValue", actStepValue_ );
        
        if(externalFiles)
          CreateExternalFile(stepNum);
      } H5_CATCH( "Can not create dataset for step " << actStepNum_ );
    }
  }

  void OutputWriter_MpCCI::EndStep() {
    OutputWriter::EndStep();

    if(settings.GetString("coupling") == "MpCCI") {
      //Send values via MpCCI
      int nQuantityIds=0;
      int quantityIds[] = {0, 0};
      int nLocalMeshIds=1;
      int localMeshIds[] = {1};
      int remoteCodeId = 2;
      int comm = CCI_COMM_RCODE[remoteCodeId];
      
      quantityIds[nQuantityIds] = quantityId1_;
      nQuantityIds++;
      quantityIds[nQuantityIds] = quantityId2_;
      nQuantityIds++;
      
      CCI_Send(nQuantityIds, quantityIds, nLocalMeshIds, localMeshIds, comm);
      CCI_Check_convergence(myConvergence,&globalConvergence,CCI_ANY_CODE);
      std::cout<<"CCI_CONTINUE_VALUE"<<CCI_CONTINUE<<std::endl;
    }
  }

  void OutputWriter_MpCCI::EndResults()
  {
  }

  int OutputWriter_MpCCI::ElemTypes2MpCCI(FEType et)
  {
    std::string elemTypeName;
    switch(et)
    {
    case ET_LINE2:
      return CCI_ELEM_LINE;
      break;

    case ET_TRIA3:
      return CCI_ELEM_TRIANGLE;
      break;

    case ET_QUAD4:
      return CCI_ELEM_QUAD;
      break;

    case ET_QUAD8:
      return CCI_ELEM_QUAD8;
      break;

    case ET_TET4:
      return CCI_ELEM_TETRAHEDRON;
      break;

    case ET_WEDGE6:
      return CCI_ELEM_PRISM;
      break;

    case ET_HEXA8:
      return CCI_ELEM_HEXAHEDRON;
      break;

    case ET_PYRA5:
      return CCI_ELEM_PYRAMID;
      break;

    default:
      Enum2String(et, elemTypeName);
      EXCEPTION("Element type " << elemTypeName
                << " cannot be used with MpCCI.");
    }
    return 0;
  }

}
>>>>>>> .merge-right.r8654
