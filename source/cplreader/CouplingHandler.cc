#include <fstream>
#include <iostream>
#include <iterator>
#include <set>
#include <string>
#include <algorithm>
#include <cmath>

#include <boost/tokenizer.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/date_time/posix_time/posix_time.hpp>
namespace fs=boost::filesystem;

#include <def_use_mpcci.hh>
#include <def_cfs_stats.hh>

#if (MpCCI_RELEASE == 305)
#include <mpcci.h>
#endif

#include "General/exception.hh"
#include "Settings.hh"
#include "General/environment.hh"
#include "DataInOut/SimInOut/hdf5/hdf5io.hh"
#include "integlib/elemIntegr.hh"

#include "FlowDataTypes.hh"
#include "OutputWriter.hh"
#include "CouplingHandler.hh"

#define H5_EXCEPTION(STR, EX)                   \
  EXCEPTION( STR, EX.getCDetailMsg() );

#define H5_CATCH( STR )                                 \
  catch (H5::Exception& h5Ex ) {                        \
    EXCEPTION( STR << ":\n" << h5Ex.getCDetailMsg() );  \
  }




namespace CoupledField
{

  CouplingHandler::CouplingHandler(shared_ptr<FileReader> ptFileReader,
                                       std::vector< shared_ptr<OutputWriter> >& outputWriters)
  {
    ptFileReader_ = ptFileReader;
    outputWriters_ = outputWriters;

    if(!ptFileReader_)
      EXCEPTION("Invalid pointer to file reader!");
  }

  CouplingHandler::~CouplingHandler()
  {

    if( mainGroup_.getLocId() <= 0 )
      return;

    // check for open groups, datasets etc.
    if (mainFile_.getObjCount( H5F_OBJ_DATASET |
                               H5F_OBJ_GROUP |
                               H5F_OBJ_DATATYPE | H5F_OBJ_ATTR) > 0 ) {
      std::cerr << "There are still objects open in the hdf5 file\n\n";
      CheckOpenObjects();
    }

    mainFile_.close();

#ifdef MpCCI
    Settings& settings = Settings::Instance();
    if(settings.GetString("coupling") == "MpCCI") {
      CCI_Finalize();
    }
#endif // MpCCI
  }

  void CouplingHandler::Init(int argc, char *argv[])
  {
    Settings& settings = Settings::Instance();

    ptFileReader_->Init();

    dim_ = settings.GetInt("dim");

    InitHDF5();
    WriteFileInfoHeader();

#ifdef MpCCI
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
#endif // MpCCI

    typedef boost::tokenizer< boost::char_separator<char> > Tok;
    boost::char_separator<char> sep(";| ");

    // Initialize vector with output fields
    Tok tokenizer(settings.GetString("outputfields"), sep);
    std::copy(tokenizer.begin(), tokenizer.end(),
              std::back_inserter(outputFields_));

    // Initialize vector with active regions
    Tok actp(settings.GetString("activeparts"), sep);
    Tok::iterator tit, tend;
    tit = actp.begin();
    tend = actp.end();

    activeParts_.resize(ptFileReader_->GetNumRegions());
    if(*tit == "all")
      std::fill(activeParts_.begin(), activeParts_.end(), true);
    else
    {
      std::stringstream sstr;
      UInt partIdx;

      for( ; tit != tend; tit++)
      {
        sstr.clear(); sstr.str("");
        sstr << *tit;
        sstr >> partIdx;

        if(partIdx > 0 && partIdx <= activeParts_.size())
          activeParts_[partIdx-1] = true;
      }
    }

    if(settings.GetInt("verbose"))
    {
      for(UInt i=0; i < activeParts_.size(); i++)
        std::cout << "Partition " << (i+1) << " active: "
                  << activeParts_[i] << std::endl;
    }

    // Initialize element integrators for source term calculation
    ptElemIntegr_[ET_LINE2]  = new ElemIntegr(ET_LINE2);
    ptElemIntegr_[ET_TRIA3]  = new ElemIntegr(ET_TRIA3);
    ptElemIntegr_[ET_QUAD4]  = new ElemIntegr(ET_QUAD4);
    ptElemIntegr_[ET_TET4]   = new ElemIntegr(ET_TET4);
    ptElemIntegr_[ET_WEDGE6] = new ElemIntegr(ET_WEDGE6);
    ptElemIntegr_[ET_PYRA5]  = new ElemIntegr(ET_PYRA5);
    ptElemIntegr_[ET_HEXA8]  = new ElemIntegr(ET_HEXA8);

  }

  void CouplingHandler::ConvertMesh()
  {
    std::cout << "========================================"
              << "========================================"
              << std::endl;
    std::cout << "                        "
              << "Entering Mesh Conversion Phase!" << std::endl;
    std::cout << "========================================"
              << "========================================"
              << std::endl;


    Settings& settings = Settings::Instance();
    std::cout << "Reading " << settings.GetString("type") << " mesh...\n";
    uint32_t dim = settings.GetInt("dim");
    int numRegions = ptFileReader_->GetNumRegions();
    std::vector<UInt>::iterator it, it2, end;
    std::set<UInt> regionNodeSet;
    std::vector<UInt> regionNodes;
    UInt maxNumElemNodes = 0;
    std::vector<Integer> feTypes; // type array for all elements (numElems)
    std::vector< UInt > numElemsOfDim ( 3 );
    UInt numElemTypes = sizeof(ELEM_DIM) / sizeof(UInt);
    std::map<FEType, UInt> numElemsOfType;
    std::vector<std::string> regionNames;

    numRegionNodes_.resize(numRegions);

    // First read everything into internal buffers
    ptFileReader_->ReadNodalCoords(nodalCoords_);
    // scale the nodal coordinates
    const UInt sizeNodCoords = nodalCoords_.size();
    std::stringstream geomstr;
    geomstr << settings.GetString("geomscale");
    Double geomScaleX, geomScaleY, geomScaleZ;
    geomstr >> geomScaleX >> geomScaleY >> geomScaleZ;
    if (geomScaleX != 1.0)
    {
      UInt i = 0;
      while (i < sizeNodCoords)
      {
        nodalCoords_[i] *= geomScaleX;
        i += 3;
      }
    }
    if (geomScaleY != 1.0)
    {
      UInt i = 1;
      while (i < sizeNodCoords)
      {
        nodalCoords_[i] *= geomScaleY;
        i += 3;
      }
    }
    if (geomScaleZ != 1.0)
    {
      UInt i = 2;
      while (i < sizeNodCoords)
      {
        nodalCoords_[i] *= geomScaleZ;
        i += 3;
      }
    }
    // <-- end scaling
    ptFileReader_->ReadTopology(topology_,
                                elemTypes_);

    // Determine the maximum number of element nodes
    maxNumElemNodes = ptFileReader_->GetMaxNumElemNodes();

    std::cout << "Writing nodal coords... ";

    // write the dimension of the grid.
    H5IO::WriteAttribute( meshGroup_, "Dimension", dim );

    // ================
    //  Node Locations
    // ================
    UInt nNodes = nodalCoords_.size()/3;
    H5::Group nodeGroup;
    try {
      nodeGroup = meshGroup_.createGroup( "Nodes" );
    } H5_CATCH( "Could not create node group" );

    H5IO::WriteAttribute( nodeGroup, "NumNodes", nNodes );

    H5IO::Write2DArray( nodeGroup, "Coordinates", nNodes,
                        3, &nodalCoords_[0], dPropList_ );

    nodeGroup.close();

    std::cout << "done.\nWriting element connectivity... ";

    // =====================
    //  Element definitions
    // =====================
    UInt nElems = elemTypes_.size();

    H5::Group elemGroup;
    try{
      elemGroup = meshGroup_.createGroup("Elements");
    } H5_CATCH( "Could not create element group" );

    // write connectivity
    H5IO::Write2DArray( elemGroup, "Connectivity", nElems,
                        maxNumElemNodes, &topology_[0], dPropList_ );

    // write element types
    H5IO::Write1DArray( elemGroup, "Types", nElems,
                        &elemTypes_[0], dPropList_ );

    std::cout << "done.\nWriting grid meta info... ";

    // ==========================
    //  Grid Meta Information
    // ==========================

    H5IO::WriteAttribute( elemGroup, "NumElems", nElems );

    it = elemTypes_.begin();
    end = elemTypes_.end();
    UInt quadrElems = 0;

    for( ; it != end; it++ )
    {
      numElemsOfDim[ ELEM_DIM[*it]-1 ]++;
      quadrElems &= ELEM_QUADRATIC[*it];
      numElemsOfType[(FEType)*it]++;
    }

    // This has still to to be checked
    H5IO::WriteAttribute( elemGroup, "QuadraticElems",
                          quadrElems );

    // number of elements per dimension
    for(UInt i=0; i<3; i++) {
      std::stringstream attrName;
      attrName << "Num" << (i+1) << "DElems";
      H5IO::WriteAttribute( elemGroup, attrName.str(), numElemsOfDim[i] );
    }

    // number of elements per type
    for(UInt i=0; i<numElemTypes; i++) {
      std::stringstream attrName;
      attrName << "Num_" << ELEM_TYPE_NAMES[i].substr(3);
      H5IO::WriteAttribute( elemGroup, attrName.str(),
                            numElemsOfType[(FEType)i] );
    }

    std::cout << "done.\n";

    // close element group
    elemGroup.close();

    std::cout << "Reading mesh done.\nConverting mesh...\n";

    for (int actRegion=0; actRegion<numRegions; actRegion++)
    {
      // Fill vector with region names
      regionNames.push_back(ptFileReader_->GetRegionName(actRegion));

      // Clear temporary containers
      regionNodeSet.clear();
      regionNodes.clear();

      std::cout << "Writing region " << (*regionNames.rbegin())
                << "... ";

      // Get all element numbers in current region and order them
      ptFileReader_->GetRegionElements(regionElems_[actRegion], actRegion);
      std::sort(regionElems_[actRegion].begin(), regionElems_[actRegion].end());

      // Put all nodes in a partition into a set to get an ordered list
      UInt idx = 0;
      UInt regionDim = 0;

      CollectElementNodes(regionElems_[actRegion], regionNodes, regionDim);

      WriteRegion(meshGroup_, regionNodes, regionElems_[actRegion], regionDim,
                  regionNames[actRegion]);

      it = regionNodes.begin();
      end = regionNodes.end();

      numRegionNodes_[actRegion] = std::distance(it, end);

      for ( idx=0; it != end; it++, idx++ )
      {
        regionNodeIndices_[actRegion][*it] = idx;
      }


#ifdef MpCCI
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
#endif // MpCCI

      std::cout << "done.\n";
      if(settings.GetInt("verbose"))
        std::cout << "Number of nodes in region " << regionNames[actRegion]
                  << " " << numRegionNodes_[actRegion] << std::endl;

    }

    // Read node and element groups
    ptFileReader_->GetNodeGroups(nodeGroups_);
    ptFileReader_->GetElemGroups(elemGroups_);

    // create new group for entity groups
    H5::Group groupsGroup;
    try{
      groupsGroup = meshGroup_.createGroup("Groups");
    } H5_CATCH( "Could not create 'Groups' group." );

    WriteNodeGroups( groupsGroup );
    WriteElemGroups( groupsGroup );
    groupsGroup.close();

    // close meshGroup
    meshGroup_.close();

#ifdef MpCCI
    if(settings.GetString("coupling") == "MpCCI") {
      //Close the definition phase; contact detection.
      //i take part in the coupling
      std::cout << "Closing MpCCI Setup" << std::endl;
      MpCCIprocess_ = 1;
      CCI_Close_setup(MpCCIprocess_);
    }
#endif // MpCCI

    std::cout << "Converting mesh done.\n";

    std::cout << "========================================"
              << "========================================"
              << std::endl;
    std::cout << "                        "
              << "Leaving Mesh Conversion Phase!" << std::endl;
    std::cout << "========================================"
              << "========================================"
              << std::endl;
  }

  void CouplingHandler::Couple()
  {
    Settings& settings = Settings::Instance();
    bool calcSrc = settings.GetInt("calcsrc");
    UInt counter = 0;
    Double stepVal = 0;
    UInt numFiles = ptFileReader_->GetNumFiles();
    UInt numRegions = ptFileReader_->GetNumRegions();
    std::vector<UInt> timeStepNumbers;
    std::vector<Double> timeStepValues;
    UInt actRegion;
    UInt stepNum = 0;
    std::vector<std::string> regionNames;
    bool externalFiles = settings.GetInt("extfiles");
    std::vector<FlowDataType> flowData(numRegions);
    bool readOK = true;

#ifdef MpCCI
    // MpCCI status variables
    int globalConvergence = CCI_CONTINUE;
    int myConvergence     = CCI_CONTINUE;
#endif // MpCCI


    std::cout << "========================================"
              << "========================================"
              << std::endl;
    std::cout << "                        "
              << "Entering Time Data Conversion Phase!" << std::endl;
    std::cout << "                        ";
    std::cout << "Number of transient files: " << numFiles << std::endl;
    std::cout << "========================================"
              << "========================================"
              << std::endl;

    // Initialize results tree in HDF5 file
    InitResultsGroup();


    // Fill vector with region names
    for (actRegion = 0; actRegion<numRegions; actRegion++)
      regionNames.push_back(ptFileReader_->GetRegionName(actRegion));

    while ( ( counter < numFiles ) && readOK)
    {
      stepVal = ptFileReader_->GetTimeStep(counter);
      stepNum = counter + 1;
      timeStepValues.push_back(stepVal);
      timeStepNumbers.push_back(stepNum);

      std::cout << "----------------------------------------"
                << "----------------------------------------"
                << std::endl;
      std::cout << "                        "
                << "Step Number: " << stepNum << "; ";
      std::cout << " Time Step: " << stepVal << std::endl;
      std::cout << "----------------------------------------"
                << "----------------------------------------"
                << std::endl;


      // Check, if step group is already open and create it if not
      if( currMeshStepGroup_.getId() <= 0 ) {
        std::stringstream stepName;
        stepName << "/Results/Mesh/MultiStep_1/Step_" << stepNum;

        // Create new step group.
        try {
          currMeshStepGroup_= mainGroup_.createGroup(stepName.str());
          H5IO::WriteAttribute( currMeshStepGroup_, "StepValue", stepVal );

          if(externalFiles)
            CreateExternalFile(stepNum);
        } H5_CATCH( "Can not create dataset for step " << stepNum );
      }

      // Read nodal values for all partitions
      try
      {
        ptFileReader_->ReadNodalValues(flowData, activeParts_, counter);
        // scale the nodal values
        const UInt sizeFlowData = flowData.size();
        std::stringstream velstr;
        velstr << settings.GetString("velscale");
        Double velScaleX, velScaleY, velScaleZ;
        velstr >> velScaleX >> velScaleY >> velScaleZ;
        if (velScaleX != 1.0)
        {
          velScale_(flowData, velScaleX, 0);
        }
        if (velScaleY != 1.0)
        {
          velScale_(flowData, velScaleY, 1);
        }
        if (velScaleZ != 1.0)
        {
          velScale_(flowData, velScaleZ, 2);
        }
        // <-- end scaling velocity

        // Override the setting of --outprec for CFX
        if(settings.GetString("type") == "CFX" && settings.GetInt("floatDataset"))
          settings.SetString("outprec", "single");

        readOK = true;
      } catch (std::exception& ex)
      {
        std::cerr << "CAUGHT EXCEPTION while trying to read nodal values:"
                  << std::endl << ex.what() << std::endl
                  << "Exiting read time values loop...";

        readOK = false;
        continue;
      }


      for (actRegion = 0; actRegion<numRegions && readOK; actRegion++)
      {
        std::string groupName;
        H5::Group currResultGroup;
        FlowDataType::iterator fdIt, fdEnd;
        UInt numDOFs;

        if(!activeParts_[actRegion])
        {
          FlowDataType::iterator fIt, fEnd;
          fIt = flowData[actRegion].begin();
          fEnd = flowData[actRegion].end();

          for( ; fIt != fEnd; fIt++ ) {
            fIt->second.isActive = false;
          }

          continue;
        }

        // If the user requests the calculation of the Lighthill
        // source term, follow his order!
        if(calcSrc)
        {
          // We need fluidMechVelocity for Lighthill source term.
          // This must be adapted for other source term formulations!!
          if ( flowData[actRegion].find(FLUIDMECH_VELOCITY)
               != flowData[actRegion].end() )
          {
            flowData[actRegion][FLUIDMECH_VELOCITY].isActive = true;
          }

          CalculateAcouSrcs(actRegion, flowData[actRegion]);
        }

        // Send fields for current partition to MpCCI
#ifdef MpCCI
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
#endif // MpCCI

        // Iterate over all flow datasets for this partition
        // and write them to HDF5
        fdIt = flowData[actRegion].begin();
        fdEnd = flowData[actRegion].end();
        for( ; fdIt != fdEnd; fdIt++)
        {
          FlowDataPartStruct& fdps = fdIt->second;
          // SolutionType st = fdIt->first;

          // Only write required datasets
          if(!fdps.isActive)
            continue;

          if(*outputFields_.begin() != "all")
          {
            if(std::find(outputFields_.begin(),
                         outputFields_.end(),
                         fdps.resultName) == outputFields_.end())
            {
              fdps.isActive = false;
              continue;
            }
          }

          // Create result groups in HDF5 file and write result.
          try {
            // Open or create subgroup for result
            try {
              currResultGroup = currMeshStepGroup_.openGroup( fdps.resultName );
            } catch (H5::GroupIException& ex)
            {
              currResultGroup = currMeshStepGroup_.createGroup( fdps.resultName );
            }

            // Create subgroup for region
            groupName = regionNames[actRegion];
            currResultGroup = currResultGroup.createGroup(groupName);

            // Create subgroup for Nodes
            groupName = H5IO::MapUnknownTypeAsString(fdps.definedOn);
            currResultGroup = currResultGroup.createGroup(groupName);

            // Write result dataset
            numDOFs = fdps.dofNames.size();
            std::cout << "Writing result: " << fdps.resultName
                      << " on region " << regionNames[actRegion] << "... ";

            std::vector<Double> output;
            ShrinkNodalVector(actRegion, numDOFs, fdps.data, output);

            WriteResults(currResultGroup, output, numDOFs, false);
            std::cout << "done." << std::endl;

            // Close nodes subgroup
            currResultGroup.close();
          } H5_CATCH("Failed to write result '" << fdps.resultName
                     << "' in step " << stepNum);
        }

      }//end of for


#ifdef MpCCI
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
#endif // MpCCI

      // Close current step group in HDF5 file
      try {
        currMeshStepGroup_.close();
      } H5_CATCH( "Could close current step group" );

      // increment time step counter
      counter++;
    }//end of while

    // Finally write out result descriptions. This must be done in
    // the end since we do not know how many steps there are in advance.
    WriteResultDescriptions( counter, flowData, timeStepNumbers,
                             timeStepValues );

    std::cout << "========================================"
              << "========================================"
              << std::endl;
    std::cout << "                        "
              << "Leaving Time Data Conversion Phase!" << std::endl;
    std::cout << "========================================"
              << "========================================"
              << std::endl;
  }

  void CouplingHandler::Finish()
  {
    Settings& settings = Settings::Instance();

    // Fetch custom data from file reader and write it to the UserData
    // section of the HDF5 file.
    std::map<std::string, std::string> userData;
    std::map<std::string, std::string>::const_iterator udIt, udEnd;

    ptFileReader_->GetUserData(userData);
    settings.DumpXML(userData["settings"]);
    udIt = userData.begin();
    udEnd = userData.end();
    for( ; udIt != udEnd; udIt++ )
      WriteStringToUserData(udIt->first, udIt->second);

    // Close main HDF5 group and finish
    try {
      mainGroup_.close();
    } H5_CATCH( "Could not close main group" );

    // Delete element integrators
    std::map<UInt, ElemIntegr *>::iterator it, end;
    it = ptElemIntegr_.begin();
    end = ptElemIntegr_.end();
    for( ; it != end; it++)
      delete it->second;
  }

  void CouplingHandler::CheckOpenObjects() {
    std::vector<UInt> types;
    std::vector<std::string> typeNames;
    hid_t* ids;

    types.push_back(H5F_OBJ_DATASET); typeNames.push_back("Dataset");
    types.push_back(H5F_OBJ_GROUP); typeNames.push_back("Group");
    types.push_back(H5F_OBJ_DATATYPE); typeNames.push_back("DataType");
    types.push_back(H5F_OBJ_ATTR); typeNames.push_back("Attribute");

    // check for open groups, datasets etc.
    std::cerr << "Number of open objects:\n"
              << "--------------------------\n";

    for(UInt t=0; t<types.size(); t++)
    {
      UInt numObjs = mainFile_.getObjCount(types[t]);
      std::cerr << typeNames[t] << "s: "<<  numObjs << std::endl;

      ids = new hid_t[numObjs];
      mainFile_.getObjIDs(types[t], numObjs, ids);

      for(UInt i=0; i<numObjs; i++)
      {
        H5::DataSet ds;
        H5::Group group;
        H5::DataType dt;
        H5::Attribute attr;

        switch(types[t])
        {
        case H5F_OBJ_DATASET:
          ds.setId((ids[i]));
          std::cerr << "  " << ds.fromClass() << std::endl;
          break;
        case H5F_OBJ_GROUP:
          group.setId((ids[i]));
          for(UInt idx=0; idx < group.getNumObjs(); idx++)
            std::cerr << "  subgroup " << group.getObjnameByIdx(idx) << std::endl;
          break;
        case H5F_OBJ_DATATYPE:
          dt.setId((ids[i]));
          std::cerr << "  " << dt.fromClass() << std::endl;
          break;
        case H5F_OBJ_ATTR:
          attr.setId((ids[i]));
          std::cerr << "  " << attr.fromClass() << std::endl;
          break;
        }
      }

      delete[] ids;
    }
  }

  void CouplingHandler::CollectElementNodes(const std::vector<UInt>& elems,
                                              std::vector<UInt>& nodes,
                                              UInt& dim)
  {
    std::vector<UInt>::const_iterator it, end;
    std::set<UInt> nodeSet;

    UInt topoIdx = 0;
    UInt idx = 0;
    UInt numElemNodes;
    FEType et;
    UInt maxNumElemNodes = ptFileReader_->GetMaxNumElemNodes();


    dim = 0;
    it = elems.begin();
    end = elems.end();

    for ( ; it != end; it++ )
    {
      idx = (*it-1);
      topoIdx = idx * maxNumElemNodes;
      et = (FEType) elemTypes_[idx];
      numElemNodes = NUM_ELEM_NODES[ et ];

      nodeSet.insert(&topology_[topoIdx],
                     &topology_[topoIdx+numElemNodes]);

      dim = dim < ELEM_DIM[et] ? ELEM_DIM[et] : dim;
    }

    std::copy(nodeSet.begin(), nodeSet.end(),
              std::back_inserter(nodes));
  }

  void CouplingHandler::InitHDF5() {
    Settings& settings = Settings::Instance();
    std::string pathsep;
    std::string fileName  = settings.GetString("name");
    std::ostringstream strBuffer;

    // Set compression and chunk size parameters for HDF5
    UInt compressionLevel = 9;
    UInt maxChunkSize = 4096;
    dPropList_ = H5::DSetCreatPropList::DEFAULT;
    dPropList_.setLayout( H5D_CHUNKED );
    dPropList_.setDeflate( compressionLevel );
    H5IO::SetMaxChunkSize( maxChunkSize );

    // Do not print HDF5 exceptions by default
    H5::Exception::dontPrint();

    hdf5DirName_ = ptFileReader_->GetPreferredOutputPath();

    // concatenate output file name
    try {
      fs::create_directory( hdf5DirName_ );
      pathsep = fs::path("/").native_directory_string();
    } catch (std::exception &ex) {
      EXCEPTION(ex.what());
    }

    strBuffer.clear();
    strBuffer.str("");
    strBuffer << hdf5DirName_ << pathsep << fileName << ".h5";
    fileName = strBuffer.str();

    // create main file and obtain main group
    try {
      mainFile_ = H5::H5File (fileName, H5F_ACC_TRUNC );
    } H5_CATCH( "Could not create hdf5 file '" << fileName << "' : " );

    mainGroup_ = mainFile_.openGroup( "/" );
    meshGroup_ = mainGroup_.createGroup( "Mesh" );
    mainGroup_.createGroup( "FileInfo" ).close();

    mainGroup_.createGroup( "UserData" );
    mainGroup_.createGroup( "Results" );
  }

  void CouplingHandler::WriteFileInfoHeader() {
    Settings& settings = Settings::Instance();
    H5::Group infoGroup;
    try {
      infoGroup = mainGroup_.openGroup( "FileInfo" );
    } H5_CATCH( "Could not open group for FileInfo" );

    // write file version
    std::stringstream version;
    version << CFS_HDF5_FORMAT_MAJOR << "." << CFS_HDF5_FORMAT_MINOR;
    std::string versionString = version.str();
    H5IO::Write1DArray( infoGroup, "Version", 1, &versionString, dPropList_ );

    // write date / time information
    using namespace boost::posix_time;
    using namespace boost::gregorian;
    std::string now = to_simple_string( second_clock::local_time() );
    H5IO::Write1DArray( infoGroup, "Date", 1, &now, dPropList_ );

    // write creator
    std::stringstream creator;
    creator << "cplreader " << CFS_VERSION << " ( "
            << CFS_SUBVERSION_REPOS
            << " rev. " << CFS_SUBVERSION_REV
            << ", date " << CFS_CONF_DATE
            << " )";
    std::string creatorString = creator.str();
    H5IO::Write1DArray( infoGroup, "Creator", 1, &creatorString, dPropList_ );

    // create group for content
    StdVector<Integer> content;
    content.Push_back( H5IO::MapCapabilityType( SimOutput::USERDATA ) );
    if(!settings.GetInt("justinit"))
      content.Push_back( H5IO::MapCapabilityType( SimOutput::MESH ) );
    if(!settings.GetInt("justmesh") && !settings.GetInt("justinit"))
      content.Push_back( H5IO::MapCapabilityType( SimOutput::MESH_RESULTS ) );
    H5IO::Write1DArray( infoGroup, "Content", content.GetSize(),
                        &content[0], dPropList_ );

  }

  void CouplingHandler::WriteRegion(const H5::Group& meshGroup,
                                      const std::vector<UInt>& nodes,
                                      const std::vector<UInt>& elems,
                                      const UInt dim,
                                      const std::string name)
  {
    H5::Group regionListGroup;

    // Try to open Regions group. If it does not exist try to create it.
    try{
      regionListGroup = meshGroup.openGroup("Regions");
    } catch(H5::GroupIException& ex) {;
      // create region group
      try{
        regionListGroup = meshGroup.createGroup("Regions");
      } H5_CATCH( "Could not create region group" );
    }

    // create new region group
    H5::Group actRegionGroup;
    try {
      actRegionGroup = regionListGroup.createGroup(name);
    } H5_CATCH( "Could not create region group for region '"
                << name << "'" );
    H5IO::WriteAttribute( actRegionGroup, "Dimension",
                          dim );

    // create new node group
    H5IO::Write1DArray<UInt>( actRegionGroup, "Nodes",
                              nodes.size(),
                              (const UInt*)&nodes[0], dPropList_ );

    // create new element group
    H5IO::Write1DArray( actRegionGroup,
                        "Elements",
                        elems.size(),
                        (const Integer*)&elems[0],
                        dPropList_);

    try{
      actRegionGroup.close();
      regionListGroup.close();
    } H5_CATCH( "Could not close region group" );
  }

  void CouplingHandler::WriteNodeGroups(const H5::Group& meshGroup) {
    H5::Group myGroup;
    std::string nodesName;
    std::map<std::string, std::vector<UInt> >::const_iterator it, end;
    it = nodeGroups_.begin();
    end = nodeGroups_.end();

    for(; it != end; it++ ) {
      nodesName = it->first;

      std::cout << "Writing node group " << nodesName << "... ";

      // try to open group with given name
      try {
        myGroup = meshGroup.openGroup( nodesName );
      } catch (H5::Exception& h5Ex ) {
        myGroup = meshGroup.createGroup( nodesName );
      }
      H5IO::WriteAttribute( myGroup, "Dimension", 0 );
      H5IO::Write1DArray( myGroup, "Nodes",
                          it->second.size(), &it->second[0], dPropList_ );

      // close nodes array of current group
      myGroup.close();

      std::cout << "done." << std::endl;
    }
  }

  void CouplingHandler::WriteElemGroups(const H5::Group& meshGroup) {
    H5::Group myGroup;
    std::vector< UInt > elemNodes;
    std::string elemsName;
    UInt dim;
    std::map<std::string, std::vector<UInt> >::const_iterator it, end;
    it = elemGroups_.begin();
    end = elemGroups_.end();

    for(; it != end; it++ ) {
      elemsName = it->first;
      std::cout << "Writing element group " << elemsName << "... ";

      CollectElementNodes(it->second, elemNodes, dim);

      try {
        myGroup = meshGroup.openGroup( elemsName );
      } catch (H5::Exception& h5Ex ) {
        myGroup = meshGroup.createGroup( elemsName );
      }
      H5IO::WriteAttribute( myGroup, "Dimension", dim );
      H5IO::Write1DArray( myGroup, "Elements",
                          it->second.size(), &it->second[0],
                          dPropList_);

      // Write nodes of element group
      H5IO::Write1DArray( myGroup, "Nodes",
                          elemNodes.size(), &elemNodes[0],
                          dPropList_);

      // close nodes array of current group
      myGroup.close();

      std::cout << "done." << std::endl;
    }
  }

  void CouplingHandler::InitResultsGroup()
  {
    Settings& settings = Settings::Instance();
    UInt externalFiles = settings.GetInt("extfiles");
    std::string analysisType = "transient";

    try {
      resultsGroup_ = mainGroup_.openGroup("Results");
    } H5_CATCH( "Could open group for results" );

    try {
      meshResultsGroup_ = resultsGroup_.createGroup("Mesh");
      H5IO::WriteAttribute( meshResultsGroup_, "ExternalFiles", externalFiles );

      H5::Group msGroup = meshResultsGroup_.createGroup("MultiStep_1");

      H5IO::WriteAttribute( msGroup, "AnalysisType", analysisType );

      msGroup.createGroup("ResultDescription").close();
      msGroup.close();
      meshResultsGroup_.close();
    } H5_CATCH( "Could init group for results" );

  }

  void CouplingHandler::WriteResultDescriptions(UInt numSteps,
      const std::vector<FlowDataType>& outputFields,
      const std::vector<UInt> stepNumbers,
      const std::vector<Double> stepValues)
  {
    std::string resultName, unit;
    UInt definedOn, numDofs, entryType;
    std::vector<std::string> dofNames;
    FlowDataType::const_iterator it, end;
    H5::Group resultDescGroup;
    H5::Group msGroup;
    std::vector<std::string> resultRegions;

    try {
      // open the group for the result description datasets.
      msGroup = mainGroup_.openGroup("/Results/Mesh/MultiStep_1");
      resultDescGroup = msGroup.openGroup("ResultDescription");

      H5IO::WriteAttribute( msGroup, "LastStepNum", *stepNumbers.rbegin() );
      H5IO::WriteAttribute( msGroup, "LastStepValue", *stepValues.rbegin() );
      msGroup.close();

    } H5_CATCH( "Could not open result description group" );

    // Extract active regionsc
    UInt numRegions = ptFileReader_->GetNumRegions();
    for(UInt i=0; i<numRegions; i++)
    {
      it = outputFields[i].begin();
      end = outputFields[i].end();

      for( ; it != end; it++ ) {
        if(it->second.isActive)
          resultRegions.push_back(ptFileReader_->GetRegionName(i));
        break;
      }
    }

    // Maybe the user made some wrong specifications...
    if(resultRegions.empty())
    {
      std::cerr << "No result description has been written!" << std::endl;
      return;
    }

    // Write result descriptions
    it = outputFields[0].begin();
    end = outputFields[0].end();

    for( ; it != end; it++ ) {
      if(!it->second.isActive)
        continue;

      resultName = it->second.resultName;
      definedOn = H5IO::MapUnknownType(it->second.definedOn);
      dofNames = it->second.dofNames;
      unit = it->second.unit;
      numDofs = dofNames.size();
      entryType = H5IO::MapEntryType(it->second.entryType);

      try {
        // === Second version: Separate datasets for each entry
        H5::Group actGroup = resultDescGroup.createGroup(resultName);

        H5IO::Write1DArray( actGroup, "DefinedOn", 1, &definedOn, dPropList_ );
        H5IO::Write1DArray( actGroup, "EntityNames", resultRegions.size(),
            &resultRegions[0], dPropList_ );
        H5IO::Write1DArray( actGroup, "NumDOFs", 1, &numDofs, dPropList_ );
        H5IO::Write1DArray( actGroup, "DOFNames", dofNames.size(),
            &dofNames[0], dPropList_ );
        H5IO::Write1DArray( actGroup, "EntryType", 1, &entryType, dPropList_ );
        H5IO::Write1DArray( actGroup, "Unit", 1, &unit, dPropList_ );

        H5IO::Write1DArray<Double>( actGroup, "StepValues",
            stepValues.size(), &stepValues[0], dPropList_ );
        H5IO::Write1DArray<UInt>( actGroup, "StepNumbers",
            stepNumbers.size(), &stepNumbers[0], dPropList_ );

        actGroup.close();

      } H5_CATCH( "Could not write result description for result '"
          << resultName << "'" );
    }

    try {
      resultDescGroup.close();
    } H5_CATCH( "Could not close result description group" );

  }

  void CouplingHandler::WriteResults( H5::Group& resultGroup,
                                        std::vector<Double>& resultVals,
                                        const UInt numDOFs,
                                        const bool isImag ) {
    Settings& settings = Settings::Instance();
    std::string outputPrecision = settings.GetString("outprec");

    // create dataset with related name
    std::string name;
    if( !isImag )
      name = "Real";
    else
      name = "Imag";

    UInt numEntities = (UInt) resultVals.size() / numDOFs;

    if(outputPrecision == "double")
    {
      H5IO::Write2DArray( resultGroup, name,
                          numEntities, numDOFs, &resultVals[0],
                          dPropList_ );
    }
    else
    {
      std::vector<Float> floatResultVals(resultVals.size());

      std::copy(resultVals.begin(), resultVals.end(),
                floatResultVals.begin());

      H5IO::Write2DArray( resultGroup, name,
                          numEntities, numDOFs, &floatResultVals[0],
                          dPropList_ );
    }
  }

  void CouplingHandler::CreateExternalFile(UInt timeStep) {
    Settings& settings = Settings::Instance();
    std::string fileName = settings.GetString("name");
    std::stringstream fName, masterGroup;
    std::string pathsep, fn;

    try {

      // open external file
      pathsep = fs::path("/").native_directory_string();

      fName << fileName << "_ms" << 1 << "_step"
            << timeStep << ".h5";
      fn = fName.str();
      fName.str("");
      fName << hdf5DirName_ << pathsep << fn;
      currStepFile_ = H5::H5File(fName.str(), H5F_ACC_TRUNC);

      // Write reference to external file to main file
      H5IO::WriteAttribute( currMeshStepGroup_, "ExtHDF5FileName", fn );

      // set current step group to external file
      currMeshStepGroup_.close();
      currMeshStepGroup_ = currStepFile_.openGroup("/");

      // Store reference to master file in external file
      fName.str("");
      fName << fileName << ".h5";
      fn = fName.str();
      H5IO::WriteAttribute( currMeshStepGroup_, "MasterHDF5FileName", fn );

      // Store reference to main group in external file
      masterGroup << "/Results/Mesh/MultiStep_" << 1
                  << "/Step_" << timeStep;
      H5IO::WriteAttribute( currMeshStepGroup_, "MasterGroup", masterGroup.str() );
    } H5_CATCH( "Could not create external file" );
  }


  void CouplingHandler::CalculateAcouSrcs(const int regionIdx,
                                            FlowDataType& flowData)
  {
    Settings& settings = Settings::Instance();

    std::string regionName = ptFileReader_->GetRegionName(regionIdx);

    if(flowData.find(FLUIDMECH_VELOCITY) == flowData.end())
    {
      std::cerr << "Cannot calculate Lighthill sources on " << regionName
                << " since no velocity field is available!" << std::endl;
      return;
    }

    // Acquire reference to velocity field
    FlowDataPartStruct& fdps = flowData[FLUIDMECH_VELOCITY];
    std::vector<Double>& velField = fdps.data;

    if(!fdps.isActive)
    {
      std::cerr << "Will not calculate Lighthill sources on " << regionName
                << " since velocity field is not active!" << std::endl;
      return;
    }

    std::cout << "Calculating Lighthill sources on " << regionName << " ";
    Double density = settings.GetDouble("density");


    // Init Lighthill structures
    FlowDataPartStruct& fdps2 = flowData[ACOU_RHS_LOAD];
    fdps2.isActive = true; // all partitions have results
    fdps2.definedOn = ResultInfo::NODE; // nodes
    if(fdps2.dofNames.empty())
      fdps2.dofNames.push_back("-");
    fdps2.unit = MapSolTypeToUnit(ACOU_RHS_LOAD);
    fdps2.resultName = "acouRhsLoad";
    fdps2.data.resize(numRegionNodes_[regionIdx]);
    fdps2.entryType = ResultInfo::SCALAR;
    std::vector<Double>& acouRhsField = fdps2.data;

    // Fill acouRhsLoad field with zeros
    std::fill(acouRhsField.begin(), acouRhsField.end(), 0);

    int nElems = ptFileReader_->GetNumElems(regionIdx);

    Matrix<Double> coordMat;
    Matrix<Double> nodaldTijdxj;
    Matrix<Double> nodalVel;
    Vector<Double> elemVec;

    FEType elemType;
    UInt numElemNodes;
    UInt elemDim;
    UInt elemIdx;
    UInt maxNENodes = ptFileReader_->GetMaxNumElemNodes();
    UInt nodeNum;

    for( int i=0; i<nElems; i++)
    {
      elemIdx = regionElems_[regionIdx][i] - 1;
      elemType = (FEType) elemTypes_[elemIdx];
      numElemNodes = NUM_ELEM_NODES[elemType];
      elemDim = ELEM_DIM[elemType];

      // Just calculate sources for volume elements!
      if(elemDim < dim_)
        continue;

      coordMat.Resize(elemDim, numElemNodes);
      nodaldTijdxj.Resize(elemDim, numElemNodes);
      nodalVel.Resize(elemDim, numElemNodes);

      for( UInt n=0; n<numElemNodes; n++)
      {
        nodeNum = topology_[elemIdx * maxNENodes + n];
        UInt topoIdx = (nodeNum - 1) * 3;
        UInt velIdx = regionNodeIndices_[regionIdx][nodeNum] * dim_;

        for( UInt d=0; d<elemDim; d++)
        {
          coordMat[d][n] = nodalCoords_[topoIdx+d];
          nodalVel[d][n] = velField[velIdx+d];
        }
      }

      try {
        ptElemIntegr_[elemType]->PerformIntegration( coordMat, nodaldTijdxj,
                                                     nodalVel, elemVec, density);
      } catch (CoupledField::Exception &ex)
      {
        std::cerr << "Warning: An Exception occurred during source term "
                  << "computation:\nElement " << i+1 << " of partition "
                  << ptFileReader_->GetRegionName(regionIdx) << std::endl;

        std::cerr << ex.what()<< std::endl;

        if(settings.GetInt("verbose")) {
          UInt oldPrec = std::cerr.precision(8);

          std::cerr << "Corner coords:\n";
          for (UInt iCol = 0, numCols = coordMat.GetSizeCol();
               iCol < numCols; ++iCol) {
            for (UInt iRow = 0, numRows = coordMat.GetSizeRow();
                 iRow < numRows; ++iRow) {
              std::cerr << "\t" << coordMat[iRow][iCol];
            }
            std::cerr << std::endl;
          }

          std::cerr.precision(oldPrec);
          std::cerr << ex.what();
        }

        std::cerr << "Setting contribution to acousrc to zero!\n\n";
        elemVec.Resize(numElemNodes);
        elemVec.Init(0.0);
      }

      // Add contributions of all element nodes
      for( UInt n=0; n<numElemNodes; n++)
      {
        UInt idx;
        nodeNum = topology_[elemIdx * maxNENodes + n];
        idx = regionNodeIndices_[regionIdx][nodeNum];

#ifdef DEBUG
        if (std::isnan(elemVec[n]) || std::isinf(elemVec[n])) {
          EXCEPTION("Source term calculated on element " << i+1
                    << " is Inf or Nan.");
        }
#endif

        acouRhsField[idx] -= elemVec[n];
      }
    }

    std::cout << "done." << std::endl;
  }

  void CouplingHandler::ShrinkNodalVector(const UInt partitionIdx,
                                            const UInt numDOFs,
                                            const std::vector<Double>& input,
                                            std::vector<Double>& output)
  {
    UInt numInputNodes = input.size() / numDOFs;
    UInt numNodes = numRegionNodes_[partitionIdx];
    std::map<UInt, UInt>::const_iterator it, end;
    UInt idxInput = 0;
    UInt idxOutput = 0;

    if(numInputNodes == numNodes)
    {
      output = input;
      return;
    }

    it = regionNodeIndices_[partitionIdx].begin();
    end = regionNodeIndices_[partitionIdx].end();
    output.resize(numNodes * numDOFs);

    for( ; it != end; it++ )
    {
      idxInput = (it->first - 1) * numDOFs;
      idxOutput = it->second * numDOFs;

      for(UInt dof=0; dof < numDOFs; dof++ )
        output[idxOutput+dof] = input[idxInput+dof];
    }
  }

  void CouplingHandler::WriteStringToUserData(const std::string& dSetName,
                                                const std::string& str) {
    H5::Group userDataGroup;

    // If it does not exist, create Group for Data.
    try {
      userDataGroup = mainGroup_.openGroup("UserData");
    } H5_CATCH( "Can not write meta information to hdf5 file" );

    H5IO::Write1DArray( userDataGroup, dSetName, 1, &str, dPropList_ );
    userDataGroup.close();
  }

  int CouplingHandler::ElemTypes2MpCCI(FEType et)
  {
#ifdef MpCCI
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
#endif // MpCCI
    return 0;
  }

} // end of namespace
