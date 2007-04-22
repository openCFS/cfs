// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <algorithm>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/convenience.hpp>
#include <boost/filesystem/exception.hpp>
#include <boost/algorithm/string/trim.hpp>

#include "DataInOut/CommandLine/BaseCommandLineHandler.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/exception.hh"
#include "simOutputXMDF.hh"

namespace fs = boost::filesystem;

namespace CoupledField {

#define CHECK_HDF5_ERROR(HID_T, STR) { if(HID_T < 0) {                  \
      std::ostringstream ostr;                                          \
      ostr << STR;                                                      \
      Exception ex(NULL, __FILE__, __LINE__, ostr.str().c_str()); } }

  SimOutputXMDF::SimOutputXMDF(std::string fileName, ParamNode * inputNode) :
    SimOutput(fileName, inputNode) 
  {
    ENTER_FCN( "SimOutputXMDF::XMDF" );
    fileName_ = fileName;
    fileId_ = -1;
    formatName_ = "hdf5";
    dirName_ = "simoutput_hdf5";

    capabilities_.insert( MESH );
    capabilities_.insert( MESH_RESULTS );

    gridWritten_ = false;
    externalFiles_ = false;

    // Change defaults according to XML file
    myParam_->Get("externalFiles", externalFiles_, false);

    // Do not print HDF5 exceptions by default
    H5::Exception::dontPrint();
  }


  SimOutputXMDF::~SimOutputXMDF() {
    ENTER_FCN( "SimOutputXMDF::~XMDF" );
    if(fileId_ > 0)
      xfCloseFile(fileId_);
  }

  void SimOutputXMDF::BeginMultiSequenceStep( UInt step,
                                              AnalysisType type,
                                              UInt numSteps  ) 
  {
    std::stringstream msName;
    H5::Attribute attrib;
    H5::DataSpace space;
    std::string analysisType;

    // If it does not exist, create Group for Data.
    try 
    {
      dataGroup_ = mainRoot_.openGroup("Data");
      volDataGroup_ = dataGroup_.openGroup("Volume");
    }
    catch (H5::Exception& h5ex)
    {
      try 
      {
        dataGroup_ = mainRoot_.createGroup("Data");
        volDataGroup_ = dataGroup_.createGroup("Volume");
      }
      catch (H5::Exception& h5ex)
      {
        EXCEPTION(h5ex.getCDetailMsg());
      }
    }

    // Assemble name of multistep
    msName << "Multistep " << step;

    try
    {
      // Create new multistep group.
      currMSGroup_ = volDataGroup_.createGroup(msName.str());

      // Add the attribute AnalysisType to multisequence step.
      attrib = currMSGroup_.createAttribute("AnalysisType",
                                            H5::StrType(H5::PredType::C_S1,
                                                        32),
                                            space);
      Enum2String(type, analysisType);
      currAnalysisType_ = type;
      analysisType.resize(32);
      attrib.write(H5::StrType(H5::PredType::C_S1, analysisType.size()),
                   analysisType);

      // Add a group for the attribute description datasets.
      currAttrDescGroup_ = currMSGroup_.createGroup("Attribute Description");

      WriteAttribDescriptions();

    }
    catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg());
    }

    currMS_ = step;
  }

  void SimOutputXMDF::RegisterResult( shared_ptr<BaseResult> sol,
                                      UInt saveBegin, UInt saveInc,
                                      UInt saveEnd ) 
  {
    registeredResults_[sol->GetResultInfo()->resultName].push_back(sol);
  }

  void SimOutputXMDF::BeginStep( UInt stepNum, Double stepVal ) 
  {
    std::stringstream stepName;
    H5::Attribute attrib;
    H5::DataSpace space;
    Double timeValue;
    Double freqValue;

    switch(currAnalysisType_)
    {
    case STATIC:

    case TRANSIENT:
      timeValue = stepVal;
      freqValue = 0;
      break;

    case HARMONIC:

    case EIGENFREQUENCY:
      timeValue = 0;
      freqValue = stepVal;
      break;

    case TRANSIENTHARMONIC:
      EXCEPTION("Don't know how to treat TRANSIENTHARMONIC in SimOutputXMDF.");
      break;

    case MULTI_SEQUENCE:
      EXCEPTION("Don't know how to treat MULTI_SEQUENCE in SimOutputXMDF.");
      break;
    }


    // Assemble name of multistep
    stepName << "Step " << stepNum;
    currStep_ = stepNum;

    try
    {
      // Create new step group.
      currStepGroup_ = currMSGroup_.createGroup(stepName.str());

      // Add the attribute StepFrequency to step.
      attrib = currStepGroup_.createAttribute("StepFrequency",
                                            H5::PredType::NATIVE_DOUBLE,
                                            space);
      attrib.write(H5::PredType::IEEE_F64LE,
                   &freqValue);
      attrib.close();

      // Add the attribute StepTime to step.
      attrib = currStepGroup_.createAttribute("StepTime",
                                            H5::PredType::NATIVE_DOUBLE,
                                            space);
      attrib.write(H5::PredType::IEEE_F64LE,
                   &timeValue);
      attrib.close();

      if(externalFiles_)
        CreateExternalFile();
    }
    catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg());
    }
  }

  void SimOutputXMDF::AddResult( shared_ptr<BaseResult> sol ) 
  {
    std::string regionName = sol->GetEntityList()->GetName();
    std::string resultName = (*sol->GetResultInfo()).resultName;
    std::stringstream regionStr;
    UInt regionNum;
    UInt numDOFs;
    std::vector< Vector<Double>* > results;
    Vector<Double> complexVec1;
    Vector<Double> complexVec2;
    std::vector<std::string> resultNames;


    if( sol->GetEntryType() == EntryType::DOUBLE ) {

      Vector<Double> & resultVec = dynamic_cast<Result<Double>&>
                                   (*sol).GetVector();

      numDOFs = (*sol->GetResultInfo()).dofNames.GetSize();

      results.push_back(&resultVec);
      resultNames.push_back(resultName);
    } else {
      Vector<Complex> & resultVec = dynamic_cast<Result<Complex>&>
                                    (*sol).GetVector();

      numDOFs = (*sol->GetResultInfo()).dofNames.GetSize();
      UInt numEntities = resultVec.GetSize() / numDOFs;

      complexVec1.Resize(resultVec.GetSize());
      complexVec2.Resize(resultVec.GetSize());

      for(UInt i=0; i<numEntities; i++)
      {
        complexVec1[i] = resultVec[i].real();
        complexVec2[i] = resultVec[i].imag();        
      }

      results.push_back(&complexVec1);
      results.push_back(&complexVec2);

      std::stringstream resName;
      if((*sol->GetResultInfo()).complexFormat == REAL_IMAG)
      {
        resName << resultName << "_real";
        resultNames.push_back(resName.str());
        resName.str("");
        resName << resultName << "_imag";
        resultNames.push_back(resName.str());
      }
      else 
      {
        resName << resultName << "_ampl";
        resultNames.push_back(resName.str());
        resName.str("");
        resName << resultName << "_phase";
        resultNames.push_back(resName.str());
      }
    }


    regionNum = ptGrid_->RegionNameToId(regionName)+1;
    regionStr << "Region " << regionNum;

    WriteAttributes(resultNames, results, regionStr.str(), numDOFs);
  }

  void SimOutputXMDF::FinishStep( ) 
  {
    if(externalFiles_)
      currStepFile_.close();

    currStepGroup_.close();
  }

  void SimOutputXMDF::FinishMultiSequenceStep( ) 
  {
    registeredResults_.clear();
    currAttrDescGroup_.close();
    currMSGroup_.close();
  }

  void SimOutputXMDF::Finalize() 
  {
    std::vector<std::string> fileNames;
    std::vector<std::string> dataSetNames;
    std::ifstream fin;
    std::string dumpStr;

    fileNames.push_back(commandLine->GetParamFile());
    fileNames.push_back(param->Get("fileFormats")
                        ->Get("materialData")
                        ->Get("file")->AsString());

    dataSetNames.push_back("Parameter File");
    dataSetNames.push_back("Material File");
    
    for(UInt i=0; i<fileNames.size(); i++)
    {
      fin.open( fileNames[i].c_str(), std::ios::binary );
    
      if(fin.fail())
        EXCEPTION("Cannot open file '" << fileNames[i]
                  <<"' to dump into HDF5!");

      // seek to the end of the file
      fin.seekg (0, std::ios::end);
      UInt numBytes = fin.tellg();
      fin.seekg (0, std::ios::beg);

      std::string str;
      str.resize(numBytes);
      fin.read(&str[0], numBytes);
      WriteStringToUserData(dataSetNames[i], str);
      fin.close();
    }

    commandLine->GetDumpString(dumpStr);
    WriteStringToUserData("Program Stats", dumpStr);
  }

  void SimOutputXMDF::InitModule()
  {
    std::string pathsep;
    std::string filename;
    std::ostringstream strBuffer;

    try 
    {
      fs::create_directory( dirName_ );
      pathsep = fs::path("/").native_directory_string();

    } catch (std::exception &ex)
    {
      EXCEPTION(ex.what());
    }

    strBuffer << dirName_ << pathsep << fileName_ << ".h5";
    filename = strBuffer.str();

    xid  xGroupId = -1;
    Integer    status;
    std::string pathToMesh = "Mesh";
    std::stringstream strBuf;

    // Create the file and the mesh group
    status = xfCreateFile(filename.c_str(), &fileId_, true);
    strBuf << "Could not create XMDF file " << fileName_;
    CHECK_HDF5_ERROR(status, strBuf.str());

    status = xfCreateGroupForMesh(fileId_, pathToMesh.c_str(), &xGroupId);
    strBuf << "Could not create mesh group '" << pathToMesh << "'";
    CHECK_HDF5_ERROR(status, strBuf.str());

    try
    {
      H5::Group group(xGroupId);
      mainRoot_ = group.openGroup("/");
    }
    catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg());
    }
  }

  void SimOutputXMDF::WriteGrid()
  {
    Integer fg_nElems;
    Integer fg_nNodes;
    std::vector<Integer> fg_ElemTypes;
    std::vector<Double> fg_XNodeLocs, fg_YNodeLocs, fg_ZNodeLocs;
    std::vector<Integer> fg_NodesInElem;
    Integer    status;
    Integer    fg_Compression = 1;
    std::stringstream strBuf;
    H5::Group mGroup;
    H5::Group elemGroup;
    xid meshGroupId;
    H5::DataSpace space;
    H5::Attribute attrib;
    UInt dim;
    UInt numElemTypes = sizeof(ELEM_DIM) / sizeof(UInt);

    if(!gridWritten_)
      InitModule();
    else
      return;

    try
    {
      mGroup = mainRoot_.openGroup("Mesh");
    }
    catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg());
    }

    // Write the dimension of the grid.
    try
    {
      attrib = mGroup.createAttribute("Dimension",
                                      H5::PredType::STD_U32LE,
                                      space);
      dim = ptGrid_->GetDim();
      attrib.write(H5::PredType::NATIVE_UINT32, &dim);
      attrib.close();
    }
    catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg());
    }

    // Get number of elements and nodes from grid.
    fg_nElems = ptGrid_->GetNumElems();
    fg_nNodes = ptGrid_->GetNumNodes();

    // Set the number of elements and nodes
    meshGroupId = mGroup.getLocId();
    status = xfSetNumberOfElements(meshGroupId, fg_nElems);
    if (status >= 0) {
      status = xfSetNumberOfNodes(meshGroupId, fg_nNodes);
    }

    if (status < 0) {
      EXCEPTION("Unable to write number of elements and nodes.");
    }

    // NodeLocations
    fg_XNodeLocs.resize(fg_nNodes);
    fg_YNodeLocs.resize(fg_nNodes);
    fg_ZNodeLocs.resize(fg_nNodes);

    for (UInt i = 0; i < fg_nNodes; i++) {
      Point p;

      ptGrid_->GetNodeCoordinate(p, i+1);

      fg_XNodeLocs[i] = p[0];
      fg_YNodeLocs[i] = p[1];
      fg_ZNodeLocs[i] = p[2];
    }

    // Write node locations
    status = xfWriteXNodeLocations(meshGroupId,
                                   fg_nNodes,
                                   &fg_XNodeLocs[0],
                                   fg_Compression);
    if (status >= 0) {
      status = xfWriteYNodeLocations(meshGroupId,
                                     fg_nNodes,
                                     &fg_YNodeLocs[0]);
      if (status >= 0) {
        status = xfWriteZNodeLocations(meshGroupId,
                                       fg_nNodes,
                                       &fg_ZNodeLocs[0]);
      }
    }

    if (status < 0) {
      EXCEPTION("Unable to write nodes.");
    }

    fg_ElemTypes.resize(fg_nElems);

    std::vector<UInt> connect;
    UInt elemNum, maxNumNodes,numNodes;
    UInt numRegions;
    FEType eType;
    RegionIdType region;
    std::vector< UInt > numElemsOfType;
    std::vector< UInt > numElemsOfDim;
    UInt quadraticElems = 0;

    maxNumNodes = ptGrid_->GetMaxNumNodesPerElem();
    connect.resize(maxNumNodes);
    numRegions = ptGrid_->GetNumRegions();
    numElemsOfType.resize(numElemTypes);
    numElemsOfDim.resize(3);

    // Build arrays for writing data
    for( UInt i=0; i < fg_nElems; i++ )
    {
      elemNum = i+1;
      std::fill(connect.begin(), connect.end(), 0);
      ptGrid_->GetElemData(elemNum, eType, region, &connect[0]);

      numElemsOfType[eType]++;
      numElemsOfDim[ELEM_DIM[eType]-1]++;

      numNodes = NUM_ELEM_NODES[eType];
      fg_ElemTypes[i] = ElemType2XMDFElemType(eType);

      switch(eType)
      {
      case ET_LINE3:
      case ET_TRIA6:
      case ET_QUAD8:
      case ET_QUAD9:
        quadraticElems = 1;
        ReorderConnectivity(eType, true, &connect[0], &connect[0]);
        break;

      case ET_TET10:
      case ET_HEXA20:
      case ET_HEXA27:
      case ET_PYRA13:
      case ET_WEDGE15:
        quadraticElems = 1;
        break;

      default:
        break;
      }

      fg_NodesInElem.insert(fg_NodesInElem.end(),
                            connect.begin(),
                            connect.end());
    }

    // Write element types
    status = xfWriteElemTypes(meshGroupId,
                              fg_nElems,
                              (int*) &fg_ElemTypes[0],
                              fg_Compression);
    if (status < 0) {
      EXCEPTION("Unable to write element types.");
    }

    // Write element connectivity
    status = xfWriteElemNodeIds(meshGroupId,
                                fg_nElems,
                                maxNumNodes,
                                (int*) &fg_NodesInElem[0],
                                fg_Compression);
    if (status < 0) {
      EXCEPTION("Unable to write element connectivity.");
    }

    // Write statistics about grid.
    try
    {
      elemGroup = mGroup.openGroup("Elements");

      attrib = elemGroup.createAttribute("QuadraticElems",
                                         H5::PredType::STD_U32LE,
                                         space);
      attrib.write(H5::PredType::NATIVE_UINT32, &quadraticElems);
      attrib.close();

      for(UInt i=0; i<3; i++)
      {
        std::stringstream attrName;
        attrName << "Num" << (i+1) << "DElems";
        attrib = elemGroup.createAttribute(attrName.str(),
                                           H5::PredType::STD_U32LE,
                                           space);
        attrib.write(H5::PredType::NATIVE_UINT32, &numElemsOfDim[i]);
        attrib.close();
      }

      for(UInt i=0; i<numElemTypes; i++)
      {
        std::stringstream attrName;
        attrName << "NUM_" << ELEM_TYPE_NAMES[i].substr(3);
        attrib = elemGroup.createAttribute(attrName.str(),
                                           H5::PredType::STD_U32LE,
                                           space);
        attrib.write(H5::PredType::NATIVE_UINT32, &numElemsOfType[i]);
        attrib.close();
      }

    }
    catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg());
    }

    try 
    {
      WriteRegions(mGroup);
      WriteNamedNodes(mGroup);
      WriteNamedElems(mGroup);
    }
    catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg());
    }

    gridWritten_ = true;
  }



  FEType SimOutputXMDF::XMDFElemType2ElemType( const Integer type ) {

    switch (type) {
    case 100:   return ET_LINE2;   // 1D linear
    case 101:   return ET_LINE3;   // 1D quadratic
    case 200:   return ET_TRIA3;   // 2D linear triangle
    case 201:   return ET_TRIA6;   // 2D quadratic triangle
    case 210:   return ET_QUAD4;   // 2D linear quadrilateral
    case 211:   return ET_QUAD8;   // 2D quadratic quadrilateral
    case 212:   return ET_QUAD9;   // 2D quadr. quad. with center node
    case 300:   return ET_TET4;    // 3D linear tetrahedron
    case 30010: return ET_TET10;   // 3D quadratic tetrahedron
    case 310:   return ET_WEDGE6;  // 3D linear prism
    case 31010: return ET_WEDGE15; // 3D quadratic prism
    case 320:   return ET_HEXA8;   // 3D linear hexahedron
    case 32010: return ET_HEXA20;  // 3D quadratic hexahedron
    case 32011: return ET_HEXA27;  // 3D quadr. hexa. with center nodes
    case 330:   return ET_PYRA5;   // 3D linear pyramid
    case 33010: return ET_PYRA13;  // 3D quadratic pyramid
    }

    return ET_UNDEF;
    // This place should never be reached!
  }

  Integer SimOutputXMDF::ElemType2XMDFElemType( const FEType type )
  {
    switch (type) {
    case ET_LINE2:   return 100;   // 1D linear
    case ET_LINE3:   return 101;   // 1D quadratic
    case ET_TRIA3:   return 200;   // 2D linear triangle
    case ET_TRIA6:   return 201;   // 2D quadratic triangle
    case ET_QUAD4:   return 210;   // 2D linear quadrilateral
    case ET_QUAD8:   return 211;   // 2D quadratic quadrilateral
    case ET_QUAD9:   return 212;   // 2D quadr. quad. with center node
    case ET_TET4:    return 300;   // 3D linear tetrahedron
    case ET_TET10:   return 30010; // 3D quadratic tetrahedron
    case ET_WEDGE6:  return 310;   // 3D linear prism
    case ET_WEDGE15: return 31010; // 3D quadratic prism
    case ET_HEXA8:   return 320;   // 3D linear hexahedron
    case ET_HEXA20:  return 32010; // 3D quadratic hexahedron
    case ET_HEXA27:  return 32011; // 3D quadr. hexa. with center nodes
    case ET_PYRA5:   return 330;   // 3D linear pyramid
    case ET_PYRA13:  return 33010; // 3D quadratic pyramid
    default: return ET_UNDEF;
    }

    // This place should never be reached!
    return -1;
  }

  void SimOutputXMDF::ReorderConnectivity( const Integer eType,
                                           const bool toXMDF,
                                           const UInt* in,
                                           UInt* out) {
    static Integer toXMDFIdxsLine[] = {0, 2, 1};
    static Integer fromXMDFIdxsLine[] = {0, 2, 1};
    static Integer toXMDFIdxsTria[] = {0, 3, 1, 4, 2, 5};
    static Integer fromXMDFIdxsTria[] = {0, 2, 4, 1, 3, 5};
    static Integer toXMDFIdxsQuad[] = {0, 4, 1, 5, 2, 6, 3, 7, 8};
    static Integer fromXMDFIdxsQuad[] = {0, 2, 4, 6, 1, 3, 5, 7, 8};

    std::vector<UInt> tmp;
    Integer* it;

    UInt numNodes = NUM_ELEM_NODES[eType];
    tmp.resize(numNodes);
    memcpy(&tmp[0], in, numNodes*sizeof(UInt));

    switch (eType) {
    case ET_LINE3: // 1D quadratic line
      if(toXMDF)
        it = toXMDFIdxsLine;
      else
        it = fromXMDFIdxsLine;
      break;
    case ET_TRIA6: // 2D quadratic triangle
      if(toXMDF)
        it = toXMDFIdxsTria;
      else
        it = fromXMDFIdxsTria;
      break;
    case ET_QUAD8: // 2D quadratic quadrilateral
    case ET_QUAD9: // 2D quadratic quadrilateral with center node
      if(toXMDF)
        it = toXMDFIdxsQuad;
      else
        it = fromXMDFIdxsQuad;
      break;
    }

    for(UInt i=0; i<numNodes; i++, it++)
    {
      out[i] = tmp[*it];
    }
  }

  void SimOutputXMDF::WriteRegions(const H5::Group& meshGroup)
  {
    hid_t status;
    H5::Group regionGroup, regionElemGroup, regionNodeGroup;
    std::vector< std::string > regionNames;
    std::vector< UInt > regionDims;
    std::vector< std::vector<UInt> > regionElems;
    std::vector< StdVector<UInt> > regionNodes;
    StdVector<Elem*> elems;
    std::vector< region_desc_type > regionDescs;
    std::vector<RegionIdType> surfRegionIds, volRegionIds;
    int fg_Compression = 1;
    UInt dim;
    UInt idx;
    hsize_t numRegions;

    try
    {
      regionGroup = meshGroup.createGroup("Regions");
      regionElemGroup = regionGroup.createGroup("Elements");
      regionNodeGroup = regionGroup.createGroup("Nodes");
    } catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg());
    }

    numRegions = ptGrid_->GetNumRegions();
    if(!numRegions)
      return;

    dim = ptGrid_->GetDim();
    ptGrid_->GetVolRegionIds(volRegionIds);
    ptGrid_->GetSurfRegionIds(surfRegionIds);
    ptGrid_->GetRegionNames(regionNames);
    regionDims.resize(numRegions);
    regionElems.resize(numRegions);
    regionNodes.resize(numRegions);

    for(UInt i=0; i<surfRegionIds.size(); i++)
    {
      idx = surfRegionIds[i];
      regionDims[idx] = dim-1;

      ptGrid_->GetElems(elems, idx);
      regionElems[idx].resize(elems.GetSize());
      for(UInt j=0; j<elems.GetSize(); j++)
      {
        regionElems[idx][j] = elems[j]->elemNum;
      }

      ptGrid_->GetNodesByRegion(regionNodes[idx], idx);

    }

    for(UInt i=0; i<volRegionIds.size(); i++)
    {
      idx = volRegionIds[i];
      regionDims[idx] = dim;

      ptGrid_->GetElems(elems, idx);
      regionElems[idx].resize(elems.GetSize());
      for(UInt j=0; j<elems.GetSize(); j++)
      {
        regionElems[idx][j] = elems[j]->elemNum;
      }

      ptGrid_->GetNodesByRegion(regionNodes[idx], idx);
    }

    regionDescs.resize(numRegions);

    for(UInt i=0; i<numRegions; i++)
    {
      strncpy(regionDescs[i].name,
              regionNames[i].c_str(),
              sizeof(regionDescs[i].name));
      regionDescs[i].dim = regionDims[i];

      status = xfWritePropertyInt(regionElemGroup.getLocId(),
                                  regionNames[i].c_str(),
                                  regionElems[i].size(),
                                  (Integer*)&regionElems[i][0],
                                  fg_Compression);
      CHECK_HDF5_ERROR(status, "Unable to write region elements.");
      status = xfWritePropertyInt(regionNodeGroup.getLocId(),
                                  regionNames[i].c_str(),
                                  regionNodes[i].GetSize(),
                                  (Integer*)&regionNodes[i][0],
                                  fg_Compression);
      CHECK_HDF5_ERROR(status, "Unable to write region nodes.");
    }


    try
    {
      H5::DataSpace space( 1, &numRegions );

      // Create a memory datatype for region_desc_type
      H5::CompType compTypeMem( sizeof(region_desc_type) );

      compTypeMem.insertMember( "Name", HOFFSET(region_desc_type, name),
                                 H5::StrType(H5::PredType::C_S1,
                                             sizeof(regionDescs[0].name)));

      compTypeMem.insertMember( "Dimension", HOFFSET(region_desc_type, dim),
                                 H5::PredType::NATIVE_UINT32);

      H5::CompType compTypeFile( sizeof(region_desc_type) );
      compTypeFile.insertMember( "Name", HOFFSET(region_desc_type, name),
                                  H5::StrType(H5::PredType::C_S1,
                                              sizeof(regionDescs[0].name)));
      compTypeFile.insertMember( "Dimension", HOFFSET(region_desc_type, dim),
                                  H5::PredType::STD_U32LE);

      H5::DataSet dataset = regionGroup.createDataSet("Description",
                                                      compTypeFile,
                                                      space);

      dataset.write(&regionDescs[0], compTypeMem);
    } catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg());
    }
  }

  void SimOutputXMDF::WriteNamedNodes(const H5::Group& meshGroup)
  {
    H5::Group nNodeGroup, namedNodeGroup;
    std::vector< named_entity_desc_type > namedNodeDescs;
    std::vector< UInt > nodes;
    std::vector<std::string> nodeNames;

    hid_t status;
    int fg_Compression = 1;
    hsize_t numNamedNodes;

    ptGrid_->GetListNodeNames(nodeNames);
    numNamedNodes = nodeNames.size();

    if(!numNamedNodes)
      return;

    try
    {
      namedNodeGroup = meshGroup.createGroup("Named Nodes");
      nNodeGroup = namedNodeGroup.createGroup("Nodes");
    } catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg());
    }

    namedNodeDescs.resize(numNamedNodes);

    for(UInt i=0; i<numNamedNodes; i++)
    {
      strncpy(namedNodeDescs[i].name,
              nodeNames[i].c_str(),
              sizeof(namedNodeDescs[i].name));

      ptGrid_->GetNodesByName(nodes, nodeNames[i]);

      status = xfWritePropertyInt(nNodeGroup.getLocId(),
                                  nodeNames[i].c_str(),
                                  nodes.size(),
                                  (int*) &nodes[0],
                                  fg_Compression);
      CHECK_HDF5_ERROR(status, "Unable to write named nodes.");
    }

    try
    {
      H5::DataSpace space( 1, &numNamedNodes );

      // Create a memory datatype for region_desc_type
      H5::CompType comptype( sizeof(named_entity_desc_type) );

      comptype.insertMember( "Name", HOFFSET(region_desc_type, name),
                             H5::StrType(H5::PredType::C_S1,
                                         sizeof(namedNodeDescs[0].name)));

      H5::DataSet dataset = namedNodeGroup.createDataSet("Description",
                                                         comptype,
                                                         space);

      dataset.write(&namedNodeDescs[0], comptype);
    } catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg());
    }
  }

  void SimOutputXMDF::WriteNamedElems(const H5::Group& meshGroup)
  {
    H5::Group nElemGroup, namedElemGroup;
    std::vector< named_entity_desc_type > namedElemDescs;
    std::vector< UInt > elems;
    StdVector< Elem* > elemPtrs;
    std::vector<std::string> elemNames;

    hid_t status;
    int fg_Compression = 1;
    hsize_t numNamedElems;

    ptGrid_->GetListElemNames(elemNames);
    numNamedElems = elemNames.size();

    if(!numNamedElems)
      return;

    try
    {
      namedElemGroup = meshGroup.createGroup("Named Elements");
      nElemGroup = namedElemGroup.createGroup("Elements");
    } catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg());
    }

    namedElemDescs.resize(numNamedElems);

    for(UInt i=0; i<numNamedElems; i++) {
      strncpy(namedElemDescs[i].name,
              elemNames[i].c_str(),
              sizeof(namedElemDescs[i].name));

      ptGrid_->GetElemsByName(elemPtrs, elemNames[i]);

      elems.resize(elemPtrs.GetSize());
      for(UInt j=0; j<elemPtrs.GetSize(); j++)
      {
        elems[j] = elemPtrs[j]->elemNum;
      }

      status = xfWritePropertyInt(nElemGroup.getLocId(),
                                  elemNames[i].c_str(),
                                  elems.size(),
                                  (int*) &elems[0],
                                  fg_Compression);
      CHECK_HDF5_ERROR(status, "Unable to write named elems.");
    }

    try {
      H5::DataSpace space( 1, &numNamedElems );

      // Create a memory datatype for region_desc_type
      H5::CompType comptype( sizeof(named_entity_desc_type) );

      comptype.insertMember( "Name", HOFFSET(region_desc_type, name),
                             H5::StrType(H5::PredType::C_S1,
                                         sizeof(namedElemDescs[0].name)));

      H5::DataSet dataset = namedElemGroup.createDataSet("Description",
                                                         comptype,
                                                         space);

      dataset.write(&namedElemDescs[0], comptype);
    } catch (H5::Exception& h5ex) {
      EXCEPTION(h5ex.getCDetailMsg());
    }
  }

  void SimOutputXMDF::WriteAttribDescriptions()
  {
    std::string resultName;
    UInt definedOn;
    UInt numDOFs;
    UInt entryType;
    hsize_t numRegions;
    std::string unit;
    std::string dofNames;
    std::vector<UInt> regions;
    RegionIdType regionId;
    ResDescType::const_iterator it, end;
    std::vector< boost::shared_ptr<BaseResult> >::const_iterator solIt, solEnd;
    boost::shared_ptr<ResultInfo> resInfo, actResInfo;

    it = registeredResults_.begin();
    end = registeredResults_.end();

    for( ; it != end; it++ )
    {
      resultName = it->first;

      solIt = it->second.begin();
      solEnd = it->second.end();

      resInfo = (*solIt)->GetResultInfo();
      numDOFs = resInfo->dofNames.GetSize();
      unit = resInfo->unit;

      switch(resInfo->definedOn)
      {
        case ResultInfo::NODE:
          definedOn = 1;
          break;
        case ResultInfo::EDGE:
          definedOn = 2;
          break;
        case ResultInfo::FACE:
          definedOn = 3;
          break;
        case ResultInfo::ELEMENT:
          definedOn = 4;
          break;
        case ResultInfo::SURF_ELEM:
          definedOn = 5;
          break;
        case ResultInfo::PFEM:
          definedOn = 6;
          break;
        case ResultInfo::REGION:
          definedOn = 7;
          break;
        case ResultInfo::SURF_REGION:
          definedOn = 8;
          break;
        case ResultInfo::NODELIST:
          definedOn = 9;
          break;
        case ResultInfo::COIL:
          definedOn = 10;
          break;
        case ResultInfo::FREE:
          definedOn = 11;
          break;
      }

      switch(resInfo->entryType)
      {
        case ResultInfo::UNKNOWN:
          entryType = 0;
          break;
        case ResultInfo::SCALAR:
          entryType = 1;
          break;
        case ResultInfo::VECTOR:
          entryType = 3;
          break;
        case ResultInfo::TENSOR:
          entryType = 6;
          break;
        case ResultInfo::STRING:
          entryType = 32;
          break;
      }

      UInt i;
      for(i = 0; i<numDOFs-1; i++)
      {
        dofNames.append(resInfo->dofNames[i]);
        dofNames.append(" ");
      }
      dofNames.append(resInfo->dofNames[i]);

      dofNames.resize(32);
      unit.resize(32);

      // Generate list of regions for the current result.
      regions.clear();
      for( ; solIt != solEnd; solIt++ )
      {
        actResInfo = (*solIt)->GetResultInfo();

        regionId = ptGrid_->RegionNameToId(
            (*solIt)->GetEntityList()->GetName());
        regions.push_back(++regionId);
      }

      // Reset solIt to beginning of result vector
      solIt = it->second.begin();

      try
      {
        hsize_t dims;
        dims=1;
        H5::DataSpace space( 1, &dims );
        numRegions = regions.size();

        UInt dtSize = sizeof(definedOn) + numRegions*sizeof(UInt) +
            sizeof(numDOFs) + dofNames.size() + sizeof(entryType) +
            unit.size();
        UInt offset = 0;

        // Create a file datatype for the attribute description
        H5::CompType compTypeFile( (size_t) dtSize );

        compTypeFile.insertMember( "DefinedOn",
                                   offset,
                                   H5::PredType::STD_U32LE);
        offset += sizeof(definedOn);

        compTypeFile.insertMember( "Regions",
                                   offset,
                                   H5::ArrayType(H5::PredType::STD_U32LE,
                                       1, &numRegions) );
        offset += numRegions*sizeof(UInt);

        compTypeFile.insertMember( "NumDOFs",
                                   offset,
                                   H5::PredType::STD_U32LE);
        offset += sizeof(numDOFs);

        compTypeFile.insertMember( "DOFNames",
                                   offset,
                                   H5::StrType(H5::PredType::C_S1,
                                               dofNames.size()));
        offset += dofNames.size();

        compTypeFile.insertMember( "EntryType",
                                   offset,
                                   H5::PredType::STD_U32LE);
        offset += sizeof(entryType);

        compTypeFile.insertMember( "Unit", offset,
                                   H5::StrType(H5::PredType::C_S1,
                                               unit.size()));

        // Create memory datatypes for attribute description
        H5::CompType compTypeMem1( (size_t) sizeof(definedOn) );
        compTypeMem1.insertMember( "DefinedOn",
                                   0,
                                   H5::PredType::NATIVE_UINT32);

        H5::CompType compTypeMem2( (size_t) numRegions*sizeof(UInt) );
        compTypeMem2.insertMember( "Regions",
                                   0, 
                                   H5::ArrayType(H5::PredType::NATIVE_UINT32,
                                       1, &numRegions) );

        H5::CompType compTypeMem3( (size_t) sizeof(numDOFs) );
        compTypeMem3.insertMember( "NumDOFs",
                                   0,
                                   H5::PredType::NATIVE_UINT32);

        H5::CompType compTypeMem4( (size_t) dofNames.size() );
        compTypeMem4.insertMember( "DOFNames",
                                   0,
                                   H5::StrType(H5::PredType::C_S1,
                                               dofNames.size()));

        H5::CompType compTypeMem5( (size_t) sizeof(entryType) );
        compTypeMem5.insertMember( "EntryType",
                                   0,
                                   H5::PredType::NATIVE_UINT32);

        H5::CompType compTypeMem6( (size_t) unit.size() );
        compTypeMem6.insertMember( "Unit",
                                   0,
                                   H5::StrType(H5::PredType::C_S1,
                                               unit.size()));

        std::vector< std::string > resNames;

        if( (*solIt)->GetEntryType() == EntryType::DOUBLE ) {
          resNames.push_back(resultName);
        }
        else
        {
          std::string resName;

          if(resInfo->complexFormat == REAL_IMAG)
          {
            resName = resultName; resName.append("_real");
            resNames.push_back(resName);
            resName = resultName; resName.append("_imag");
            resNames.push_back(resName);
          }
          else
          {
            resName = resultName; resName.append("_ampl");
            resNames.push_back(resName);
            resName = resultName; resName.append("_phase");
            resNames.push_back(resName);
          }
        }

        for(i=0; i<resNames.size(); i++)
        {
          H5::DataSet dataset = currAttrDescGroup_.createDataSet(resNames[i],
              compTypeFile,
              space);

          dataset.write(&definedOn, compTypeMem1);
          dataset.write(&regions[0], compTypeMem2);
          dataset.write(&numDOFs, compTypeMem3);
          dataset.write(dofNames.c_str(), compTypeMem4);
          dataset.write(&entryType, compTypeMem5);
          dataset.write(unit.c_str(), compTypeMem6);
        }

      } catch (H5::Exception& h5ex2)
      {
        EXCEPTION(h5ex2.getCDetailMsg());
      }
    }
  }

  void SimOutputXMDF::WriteAttributes(const std::vector<std::string>& resultNames,
                                      std::vector< Vector<Double>* >& results,
                                      const std::string& regionStr,
                                      const UInt numDOFs)
  {
    H5::Group resultGroup;
    UInt numResults = resultNames.size();

    for(UInt i=0; i<numResults; i++)
    {
      UInt numEntities = results[i]->GetSize() / numDOFs;

      try
      {
        // Try to open result group.
        resultGroup = currStepGroup_.openGroup(resultNames[i]);
      }
      catch (H5::Exception& h5ex)
      {
        try
        {
          resultGroup = currStepGroup_.createGroup(resultNames[i]);
        }
        catch (H5::Exception& h5ex)
        {
          EXCEPTION(h5ex.getCDetailMsg());
        }
      }

      try
      {
        // dim sizes of ds (on disk)
        hsize_t fdim[] = {numEntities, numDOFs}; 
        H5::DataSpace fspace( 2, fdim );

        // Create dataset and write it into the file.
        H5::DataSet dataset;
        dataset = resultGroup.createDataSet( regionStr,
                                             H5::PredType::NATIVE_DOUBLE,
                                             fspace );
        dataset.write( &(*results[i])[0], H5::PredType::IEEE_F64LE);

        dataset.close();
        fspace.close();
      }
      catch (H5::Exception& h5ex)
      {
        EXCEPTION(h5ex.getCDetailMsg());
      }
    }

    resultGroup.close();
  }

  void SimOutputXMDF::CreateExternalFile()
  {
    H5::Attribute attrib;
    H5::DataSpace space;
    std::stringstream fName;
    std::string pathsep;
    std::string fn;
    std::stringstream masterGroup;
    std::string mg;

    try 
    {

      pathsep = fs::path("/").native_directory_string();

      fName << fileName_ << "_ms" << currMS_ << "_step"
          << currStep_ << ".h5";
      fn = fName.str();
      fName.str("");
      fName << dirName_ << pathsep << fn;

      currStepFile_ = H5::H5File(fName.str(), H5F_ACC_TRUNC);

      // Add the attribute StepFrequency to step.
      attrib = currStepGroup_.createAttribute("ExtHDF5FileName",
                                              H5::StrType(H5::PredType::C_S1,
                                                  fn.size()),
                                              space);
      attrib.write(H5::StrType(H5::PredType::C_S1,
                   fn.size()),
      fn.c_str());
      attrib.close();
      currStepGroup_.close();

      currStepGroup_ = currStepFile_.openGroup("/");

      fName.str("");
      fName << fileName_ << ".h5";
      fn = fName.str();

      // Add the attribute StepFrequency to step.
      attrib = currStepGroup_.createAttribute("MasterHDF5FileName",
                                              H5::StrType(H5::PredType::C_S1,
                                                  fn.size()),
                                              space);
      attrib.write(H5::StrType(H5::PredType::C_S1,
                   fn.size()),
      fn.c_str());
      attrib.close();

      masterGroup << "/Data/Volume/Multistep " << currMS_
          << "/Step " << currStep_;
      mg = masterGroup.str();

      // Add the attribute StepFrequency to step.
      attrib = currStepGroup_.createAttribute("MasterGroup",
                                              H5::StrType(H5::PredType::C_S1,
                                                  mg.size()),
                                              space);
      attrib.write(H5::StrType(H5::PredType::C_S1,
                   mg.size()),
      mg.c_str());
      attrib.close();
    }
    catch (H5::Exception& h5ex)
    {
      EXCEPTION(h5ex.getCDetailMsg());
    }
  }

  void SimOutputXMDF::WriteStringToUserData(const std::string& dSetName, 
                                            const std::string& str)
  {
    H5::Group userDataGroup;
    
    // If it does not exist, create Group for Data.
    try 
    {
      userDataGroup = dataGroup_.openGroup("User Data");
    }
    catch (H5::Exception& h5ex)
    {
      try 
      {
        userDataGroup = dataGroup_.createGroup("User Data");
      }
      catch (H5::Exception& h5ex)
      {
        EXCEPTION(h5ex.getCDetailMsg());
      }
    }

    H5::DataSpace space;

    H5::DataSet dataset;
    dataset = userDataGroup.createDataSet(dSetName,
                                          H5::StrType(H5::PredType::C_S1,
                                                      str.length()),
                                          space);
    
    dataset.write(str.c_str(), H5::StrType(H5::PredType::C_S1, str.length()));
    
    userDataGroup.close();
  }
}

