// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <cstdio>

#include <string>

// Include headers which define what types of in/output files CFS++ supports
#include <def_use_mesh.hh>
#include <def_use_gidpost.hh>
#include <def_use_hdf5.hh>
#include <def_use_gmv.hh>
#include <def_use_gmsh.hh>
#include <def_use_unv.hh>
#include <def_use_ansysrst.hh>
#include <def_use_comsol.hh>
#include <def_use_cgns.hh>
#include <def_use_ensight.hh>

#include "DefineInOutFiles.hh"

#include "DataInOut/SimInOut/AnsysCDB/SimInputCDB.hh"

#ifdef USE_MESH
#include "DataInOut/SimInOut/AnsysFile/SimInputMESH.hh"
#include "DataInOut/SimInOut/internalMesh/InternalMesh.hh"
#endif

#ifdef USE_GMV
#include "DataInOut/SimInOut/gmv/SimInputGMV.hh"
#include "DataInOut/SimInOut/gmv/SimOutGMV.hh"
#endif

#ifdef USE_GMSH
#include "DataInOut/SimInOut/gmsh/SimInputGmsh.hh"
#include "DataInOut/SimInOut/gmsh/SimOutputGmsh.hh"
#include "DataInOut/SimInOut/gmsh/SimOutputParsed.hh"
#endif

#ifdef USE_HDF5
// HDF5 readers and writers
#include "DataInOut/SimInOut/hdf5/SimInputHDF5.hh"
#include "DataInOut/SimInOut/hdf5/SimOutputHDF5.hh"
#endif

#include "DataInOut/SimInOut/RefElems/SimInputRefElems.hh"

#ifdef USE_GIDPOST
#include "DataInOut/SimInOut/GiD/SimOutGiD.hh"
#endif

#ifdef USE_UNV
#include "DataInOut/SimInOut/Unverg/SimInputUnv.hh"
#include "DataInOut/SimInOut/Unverg/SimOutputUnv.hh"
#endif

#ifdef USE_ANSYSRST
#include "DataInOut/SimInOut/AnsysRST/SimOutputRST.hh"
#endif

#ifdef USE_COMSOL
#include "DataInOut/SimInOut/COMSOL/SimInputMPHTXT.hh"
#endif

#ifdef USE_ENSIGHT
#include "DataInOut/SimInOut/VTKBased/Ensight/SimInputEnsight.hh"
#endif

#ifdef USE_CGNS
#include "DataInOut/SimInOut/CGNS/SimInputCGNS.hh"
#include "DataInOut/SimInOut/CGNS/SimOutputCGNS.hh"
#endif

#include "DataInOut/SimInOut/TextOutput/TextSimOutput.hh"
#include "DataInOut/SimInOut/InfoResultOutput/SimOutputInfo.hh"
#ifndef __MINGW32__
#ifndef __INTEL_COMPILER
#include "DataInOut/SimInOut/Streaming/SimOutputStreaming.hh"
#endif
#endif

#include "DataInOut/ParamHandling/XMLMaterialHandler.hh"

#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{

  // ===============
  //   Constructor
  // ===============
  DefineInOutFiles::DefineInOutFiles() {
    // Initialise internal pointers
    ptMaterialHandler_ = NULL;
  }
  
  // ==============
  //   Destructor
  // ==============
  DefineInOutFiles::~DefineInOutFiles()
  {
    delete ptMaterialHandler_;
    ptMaterialHandler_ = NULL;
  }
  
// ==============================
//   Generate mesh file pointer
// ==============================
void DefineInOutFiles::CreateSimInputFiles(PtrParamNode rootNode,
                                           PtrParamNode infoNode,
                                           std::map<std::string, shared_ptr<
    SimInput> >& inFiles, std::map<std::string,
    StdVector<shared_ptr<SimInput> > >& gridInputs)
{

  std::string meshFile = progOpts->GetMeshFileStr();
  std::string simName = progOpts->GetSimName();
  std::string fileName = "default";
  std::string actId, actGridId;

  // resest map
  inFiles.clear();
  gridInputs.clear();

  std::string informat = "mesh";
  StdVector<PtrParamNode> inputOptionNodes;
  PtrParamNode inputNode = rootNode->Get("fileFormats") ->Get("input",ParamNode::INSERT);
  inputOptionNodes = inputNode->GetChildren();

  // if no reader is defined explictly, create implicit one
  if (!inputOptionNodes.GetSize())
  {
    actId = "default";
    actGridId = "default";
    if (meshFile.empty())
      meshFile = simName + ".mesh";
    
    PtrParamNode meshNode = inputNode->Get("mesh", ParamNode::INSERT);
    meshNode->GetValue("id", actId, ParamNode::INSERT);
    meshNode->GetValue("gridId", actGridId, ParamNode::INSERT);
    
    if ( meshFile.find(".h5", meshFile.length()-4) != std::string::npos ) {
      inFiles[actId] = shared_ptr<SimInput>( new SimInputHDF5(meshFile, PtrParamNode(new ParamNode()), infoNode) );
    } else if ( meshFile.find(".msh", meshFile.length()-5) != std::string::npos ){
      inFiles[actId] = shared_ptr<SimInput>( new SimInputGmsh(meshFile, PtrParamNode(new ParamNode()), infoNode) );
    } else {
      inFiles[actId] = shared_ptr<SimInput>( new SimInputMESH(meshFile, PtrParamNode(), infoNode) );
    }
    gridInputs[actGridId].Push_back(inFiles[actId]);
    return;
  }

  for (UInt i = 0; i < inputOptionNodes.GetSize(); i++)
  {

    // fetch format and id of output class
    PtrParamNode actNode = inputOptionNodes[i];
    informat = actNode->GetName();
    actId = actNode->Get("id")->As<std::string>();
    actNode->GetValue("fileName", fileName, ParamNode::PASS);
    actNode->GetValue("gridId", actGridId, ParamNode::EX);

    if (i == 0)
    {
      if ((meshFile.empty()) && (fileName != "default"))
        meshFile = fileName;
    }
    else
      meshFile = fileName;

    // ensure, that id is unique
    if (inFiles.find(actId) != inFiles.end())
    {
      EXCEPTION( "Id '" << actId
          << "' for input format '"<< informat
          << "' was already found!\n"
          << "Please ensure, that the ids for the "
          << "input entries are unique!" );
    }

    inFiles[actId] = CreateSingleInputFileObject(meshFile,simName,actNode,infoNode);
    // relate gridId with input reader id
    gridInputs[actGridId].Push_back(inFiles[actId]);
  }
}

// ================================
//   Generate output file pointer
// ================================
void DefineInOutFiles::
CreateSimOutputFiles(PtrParamNode rootNode,
                     PtrParamNode infoNode,
                     std::map<std::string, shared_ptr<SimOutput> >& out,
                     std::map<std::string, std::string> & gridIds )
{

  // resest map
  out.clear();

  std::string simName = progOpts->GetSimName();
  
  // check for restart
  bool restart = progOpts->GetRestart();

  // get list of output formats
  PtrParamNode outNode = rootNode->Get("fileFormats")->Get("output", ParamNode::PASS);

  if (!outNode)
  {
    WARN("There was no output writer specified at all");
  }
  ParamNodeList formatNodes = outNode->GetChildren();

  // Check if only  one reader per format has been defined  and that the given
  // ids are unique.
  std::set<std::string> formatSet;
  std::set<std::string> idSet;
  std::string actFormat, actId, hdf5Id;
  for (UInt i = 0; i < formatNodes.GetSize(); i++)
  {
    // fetch format and id of output class
    PtrParamNode actNode = formatNodes[i];
    actFormat = actNode->GetName();
    actId = actNode->Get("id")->As<std::string>();

    // Note: In general, we should ensure, that output writers exist only once
    // (especially hdf5, gmv ,etc.). But for text-writers, it makes intentionally
    // sense to have several ones (e.g. for collecting over frequency or over
    // space). We should therefore restrict this check for format like 
    // gmv and hdf5.
//    // ensure, that format is unique
//    if (formatSet.find(actFormat) != formatSet.end())
//    {
//      EXCEPTION( "Several output tags for format '" << actFormat
//                 << "' were found!\n"
//                 << "Please ensure, that the output tags for the "
//                 << "formats are unique!" );
//    } else 
//    {
//      formatSet.insert(actFormat);
//    }

    if(actFormat == "hdf5")
      hdf5Id = actId;

    // ensure, that id is unique
    if (idSet.find(actId) != idSet.end())
    {
      EXCEPTION( "Id '" << actId
          << "' for output format '"<< actFormat
          << "' was already found!\n"
          << "Please ensure, that the ids for the "
          << "output entries are unique!" );
    } else 
    {
      idSet.insert(actId);
    }
    
  }
  
  // iterate over all found files
  for (UInt i = 0; i < formatNodes.GetSize(); i++)
  {
    // fetch format and id of output class
    PtrParamNode actNode = formatNodes[i];
    actFormat = actNode->GetName();
    actId = actNode->Get("id")->As<std::string>();

    // Read in the gridId
    std::string gridId = "default";
    actNode->GetValue( "gridId", gridId, ParamNode::PASS );
    gridIds[actId] = gridId;
    
    if(!progOpts->IsQuiet())
      std::cout << "++ Creating " << actFormat << " writer with ID '" << actId << "'" << std::endl;
    out[actId]  = CreateSingleOutputFileObject(simName,actNode,infoNode,restart);

  } // loop over reader nodes
}

// ==================================
//   Generate material file handler
// ==================================
MaterialHandler *
DefineInOutFiles::CreateMaterialHandler(PtrParamNode rootNode )
{

  std::string fileName = "mat.dat";
  std::string format = "dat";

  // Determine filename and format
  PtrParamNode matNode = 
      rootNode->Get("fileFormats")->Get("materialData", ParamNode::PASS);
  if (matNode)
  {
    matNode->GetValue("file", fileName);
    matNode->GetValue("format", format);
  }

  if (format == "dat")
  {
    //ptMaterialHandler_ = new PlainMaterialHandler( fileName );
    EXCEPTION("I am really sorry to tell you, but we have just "
        << "abandoned the dat-format ... FOREVER!");
  }
  else if (format == "xml")
  {
    XMLMaterialHandler * xmlHandler = new XMLMaterialHandler();
    xmlHandler->LoadFromFile(fileName);
    ptMaterialHandler_ = xmlHandler;
  }
  else
  {
    EXCEPTION( "CreateMaterialHandler: Format '" << format
        << "' is not recognized!" );
  }
  return ptMaterialHandler_;

}

shared_ptr<SimInput>  DefineInOutFiles::CreateSingleInputFileObject(std::string fName,
                                                                    std::string simName,
                                                                    PtrParamNode configNode,
                                                                    PtrParamNode infoNode){
  shared_ptr<SimInput> aInput;

  std::string fFormat = configNode->GetName();

  if (fFormat == "mesh")
  {
#ifdef USE_MESH
    if(fName.empty()){
      fName += simName + ".mesh";
    }
    aInput = shared_ptr<SimInput> (
        new SimInputMESH(fName, configNode, infoNode));
#else
    EXCEPTION( "No support for MESH input file format." );
#endif // USE_MESH
  }
  else if (fFormat == "cdb")
  {
    if(fName.empty()){
      fName += simName + ".cdb";
    }
    aInput = shared_ptr<SimInput> (
        new SimInputCDB(fName, configNode, infoNode));
  }
  else if (fFormat == "hdf5")
  {
#ifdef USE_HDF5
    if(fName.empty()){
      fName += simName + ".h5";
    }
    aInput = shared_ptr<SimInput> (
        new SimInputHDF5(fName, configNode, infoNode));
#else
    EXCEPTION( "No support for HDF5 input file format." );
#endif // USE_HDF5
  }
  else if (fFormat == "gmv")
  {
#ifdef USE_GMV
    if(fName.empty()){
      fName += simName + ".gmv";
    }
    aInput = shared_ptr<SimInput>(new SimInputGMV(fName,
        configNode, infoNode));
#else
    EXCEPTION( "No support for GMV input file format." );
#endif // USE_GMV
  }
  else if (fFormat == "refelem")
  {
    if(fName.empty()){
      fName += simName + ".refelem";
    }
    aInput = shared_ptr<SimInput>(new SimInputRefElems(fName,
        configNode, infoNode));
  }
  else if (fFormat == "gmsh")
  {
#ifdef USE_GMSH
    if(fName.empty()){
      fName += simName + ".msh";
    }
    aInput = shared_ptr<SimInput> (
        new SimInputGmsh(fName, configNode, infoNode));
#else
    EXCEPTION( "No support for GMSH input file format." );
#endif // USE_GMSH
  }
  else if (fFormat == "mphtxt")
  {
#ifdef USE_COMSOL
    if(fName.empty()){
      fName += simName + ".mphtxt";
    }
    aInput = shared_ptr<SimInput> (
        new SimInputMPHTXT(fName, configNode, infoNode));
#else
    EXCEPTION( "No support for Comsol .mphtxt input file format." );
#endif // USE_COMSOL
  }
  else if (fFormat == "cgns")
  {
#ifdef USE_CGNS
    if(fName.empty()){
      fName += simName + ".cgns";
    }
    aInput = shared_ptr<SimInput> (
        new SimInputCGNS(fName, configNode, infoNode));
#else
    EXCEPTION( "No support for CGNS .cgns input file format." );
#endif // USE_CGNS
  }
  else if (fFormat == "unv")
  {
#ifdef USE_UNV
    if(fName.empty()){
      fName += simName + ".unv";
    }
    aInput = shared_ptr<SimInput> (new SimInputUnv(fName, configNode, infoNode));
#else
    EXCEPTION( "No support for UNV input file format." );
#endif // USE_UNV
  }
  else if (fFormat == "ensight")
  {
#ifdef USE_ENSIGHT
    if(fName.empty()){
      fName += simName + ".case";
    }
    aInput = shared_ptr<SimInput> (new SimInputEnsight(fName, configNode, infoNode));
#else
    EXCEPTION( "No support for ENSIGHT Gold input file format." );
#endif // USE_UNV
  }
  else if (fFormat == "internal")
  {
#ifdef USE_MESH
    if(fName.empty()){
      fName += simName + ".mesh";
    }
    aInput =
        shared_ptr<SimInput> (new InternalMesh(fName, configNode, infoNode));
#else
    EXCEPTION( "No support for internalMesh input file format." );
#endif // USE_MESH
  }
  else
  {
    EXCEPTION( "Wrong format for input file. Please, check your data!" );
  }

  return aInput;

}

shared_ptr<SimOutput> DefineInOutFiles::CreateSingleOutputFileObject(std::string fName,
                                                                     PtrParamNode configNode,
                                                                     PtrParamNode infoNode,
                                                                     bool isRestart){
  shared_ptr<SimOutput> aOutput;
  std::string fFormat = configNode->GetName();
  if (fFormat == "unv")
  {
#ifdef USE_UNV
    aOutput = shared_ptr<SimOutput> (new SimOutputUnv(fName, configNode,
                                                      infoNode, isRestart));
#else
    EXCEPTION( "No support for UNV output file format." );
#endif
  }

  if (fFormat == "cgns")
  {
#ifdef USE_CGNS
    aOutput =   shared_ptr<SimOutput> (new SimOutputCGNS(fName, configNode,
                                                         infoNode, isRestart));
#else
    EXCEPTION( "No support for CGNS output file format." );
#endif
  }

  if (fFormat == "gid")
  {
#ifdef USE_GIDPOST
    aOutput =   shared_ptr<SimOutput> (new SimOutputGiD(fName, configNode,
                                                        infoNode, isRestart));
#else
    EXCEPTION( "No support for GiD output file format." );
#endif
  }

  if (fFormat == "gmsh")
  {
#ifdef USE_GMSH
    aOutput =  shared_ptr<SimOutput> (new SimOutputGmsh(fName, configNode,
                                                        infoNode, isRestart));
#else
    EXCEPTION( "No support for Gmsh output file format." );
#endif
  }


  if (fFormat == "gmshParsed")
  {
#ifdef USE_GMSH
    aOutput =  shared_ptr<SimOutput> (new SimOutputParsed(fName, configNode,
                                                          infoNode, isRestart));
#else
    EXCEPTION( "No support for Gmsh parsed output file format." );
#endif
  }

  if (fFormat == "gmv")
  {
#ifdef USE_GMV
    aOutput = shared_ptr<SimOutput> (new SimOutputGMV(fName, configNode,
                                                      infoNode, isRestart));
#else
    EXCEPTION( "No support for GMV output file format." );
#endif
  }

  if (fFormat == "hdf5")
  {
#ifdef USE_HDF5
    aOutput = shared_ptr<SimOutput> (new SimOutputHDF5(fName, configNode,
                                                       infoNode, isRestart));

#else
    EXCEPTION( "No support for HDF5 output file format." );
#endif
  }

  if (fFormat == "rst")
  {
#ifdef USE_ANSYSRST
    aOutput =
    shared_ptr<SimOutput>( new SimOutputRST( fName, configNode, infoNode, isRestart ) );
#else
    EXCEPTION( "No support for ANSYS RST output file format." );
#endif
  }

  if (fFormat == "text" || fFormat == "csv")
  {
    aOutput = shared_ptr<SimOutput> (new SimOutputText(fName, configNode,
                                                       infoNode, isRestart));
  }

  if (fFormat == "info")
  {
    aOutput = shared_ptr<SimOutput> (new SimOutputInfo(configNode,infoNode, isRestart));
  }

#ifndef __MINGW32__
#ifndef __INTEL_COMPILER
  if (fFormat == "streaming")
  {
    aOutput = shared_ptr<SimOutput> (new SimOutputStreaming(configNode, infoNode, isRestart));
  }
#endif
#endif

  return aOutput;
}


} // end of namespace
