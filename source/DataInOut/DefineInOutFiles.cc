// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <cstdio>

#include <string>
#include <boost/make_shared.hpp>

// Include headers which define what types of in/output files openCFS supports
#include <def_use_gidpost.hh>
#include <def_use_unv.hh>
#include <def_use_cgns.hh>
#include <def_use_ensight.hh>
#include <def_use_embedded_python.hh>

#include "DefineInOutFiles.hh"

#include "DataInOut/SimInOut/AnsysCDB/SimInputCDB.hh"
#include "DataInOut/SimInOut/AnsysFile/SimInputMESH.hh"
#include "DataInOut/SimInOut/internalMesh/InternalMesh.hh"

#include "DataInOut/SimInOut/gmsh/SimInputGmsh.hh"
#include "DataInOut/SimInOut/gmsh/SimOutputGmsh.hh"
#include "DataInOut/SimInOut/gmsh/SimOutputParsed.hh"

// HDF5 readers and writers
#include "DataInOut/SimInOut/hdf5/SimInputHDF5.hh"
#include "DataInOut/SimInOut/hdf5/SimOutputHDF5.hh"

#include "DataInOut/SimInOut/RefElems/SimInputRefElems.hh"

#ifdef USE_GIDPOST
#include "DataInOut/SimInOut/GiD/SimOutGiD.hh"
#endif

#include "DataInOut/SimInOut/Unverg/SimInputUnv.hh"
#include "DataInOut/SimInOut/Unverg/SimOutputUnv.hh"

#ifdef USE_ENSIGHT
#include "DataInOut/SimInOut/VTKBased/Ensight/SimInputEnsight.hh"
#endif

#ifdef USE_CGNS
#include "DataInOut/SimInOut/CGNS/SimInputCGNS.hh"
#include "DataInOut/SimInOut/CGNS/SimOutputCGNS.hh"
#endif

#ifdef USE_EMBEDDED_PYTHON
  #include "DataInOut/SimInOut/python/SimInputPython.hh"
#endif

#include "DataInOut/SimInOut/TextOutput/TextSimOutput.hh"
#include "DataInOut/SimInOut/InfoResultOutput/SimOutputInfo.hh"
#ifdef USE_STREAMING
  #include "DataInOut/SimInOut/Streaming/SimOutputStreaming.hh"
#endif
#include "DataInOut/ParamHandling/XMLMaterialHandler.hh"

#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

using std::string;

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
  DefineInOutFiles::~DefineInOutFiles() {}
  
// ==============================
//   Generate mesh file pointer
// ==============================
void DefineInOutFiles::CreateSimInputFiles(PtrParamNode rootNode,
                                           PtrParamNode infoNode,
                                           std::map<string, shared_ptr<SimInput> >& inFiles,
                                           std::map<string, StdVector<shared_ptr<SimInput> > >& gridInputs)
{
  fs::path mfp = progOpts->GetMeshFile();
  string meshFile = mfp.string();
  string simName = progOpts->GetSimName();
  string fileName = "default";
  string actId, actGridId;

  // resest map
  inFiles.clear();
  gridInputs.clear();

  string informat = "mesh";
  StdVector<PtrParamNode> inputOptionNodes;
  PtrParamNode inputNode = rootNode->Get("fileFormats") ->Get("input",ParamNode::INSERT);
  inputOptionNodes = inputNode->GetChildren();

  // if no reader is defined explicitly, create implicit one
  if (!inputOptionNodes.GetSize())
  {
    actId = "default";
    actGridId = "default";
    if (meshFile.empty())
      meshFile = simName + ".mesh";
    
    PtrParamNode meshNode = inputNode->Get("mesh", ParamNode::INSERT);
    meshNode->GetValue("id", actId, ParamNode::INSERT);
    meshNode->GetValue("gridId", actGridId, ParamNode::INSERT);
    
    // we assume inFiles was not set before
    string extension = mfp.extension().string();
    if(extension == ".h5" || extension == ".cfs")
      inFiles[actId] = shared_ptr<SimInput>(new SimInputHDF5(meshFile, PtrParamNode(new ParamNode()), infoNode));
    else if(extension == ".cdb")
      inFiles[actId] = shared_ptr<SimInput>(new SimInputCDB(meshFile, PtrParamNode(new ParamNode()), infoNode));
    else // even if this is not .mesh, we need a SimInput, otherwise it fails later - is rather stupid :(
      inFiles[actId] = shared_ptr<SimInput>(new SimInputMESH(meshFile, PtrParamNode(), infoNode));

    gridInputs[actGridId].Push_back(inFiles[actId]);
    return;
  }

  for (UInt i = 0; i < inputOptionNodes.GetSize(); i++)
  {

    // fetch format and id of output class
    PtrParamNode actNode = inputOptionNodes[i];
    informat = actNode->GetName();
    actId = actNode->Get("id")->As<string>();
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
                     std::map<string, shared_ptr<SimOutput> >& out,
                     std::map<string, string> & gridIds )
{

  // resest map
  out.clear();

  string simName = progOpts->GetSimName();
  
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
  std::set<string> formatSet;
  std::set<string> idSet;
  string actFormat, actId, hdf5Id;
  for (UInt i = 0; i < formatNodes.GetSize(); i++)
  {
    // fetch format and id of output class
    PtrParamNode actNode = formatNodes[i];
    actFormat = actNode->GetName();
    actId = actNode->Get("id")->As<string>();

    // Note: In general, we should ensure, that output writers exist only once
    // (especially hdf5, etc.). But for text-writers, it makes intentionally
    // sense to have several ones (e.g. for collecting over frequency or over
    // space). We should therefore restrict this check for format like hdf5.
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
    actId = actNode->Get("id")->As<string>();

    // Read in the gridId
    string gridId = "default";
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
shared_ptr<MaterialHandler>
DefineInOutFiles::CreateMaterialHandler(PtrParamNode rootNode )
{
  fs::path root = progOpts->ObtainCFSRootFromSystem();
  root.normalize();                                                                          // shall be save to remove when the depreciaton hurts
  string fileName = root.string() + "/share/xml/CFS-Material/Examples/MaterialDataBase.xml"; // shall work also on Windows
  string format = "xml";

  // Determine filename and format
  PtrParamNode matNode = 
      rootNode->Get("fileFormats")->Get("materialData", ParamNode::PASS);
  if (matNode)
  {
    matNode->GetValue("file", fileName);
    matNode->GetValue("format", format);
  }
  else
  {
    WARN("No materialData tag given in XML input. CFS is using the default material database in '" << fileName << "'");
  }

  if (format == "dat")
  {
    //ptMaterialHandler_ = new PlainMaterialHandler( fileName );
    EXCEPTION("I am really sorry to tell you, but we have just "
        << "abandoned the dat-format ... FOREVER!");
  }
  else if (format == "xml")
  {
    shared_ptr<XMLMaterialHandler> xmlHandler = boost::make_shared<XMLMaterialHandler>();
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

shared_ptr<SimInput>  DefineInOutFiles::CreateSingleInputFileObject(string fName,
                                                                    string simName,
                                                                    PtrParamNode configNode,
                                                                    PtrParamNode infoNode){
  shared_ptr<SimInput> aInput;

  string fFormat = configNode->GetName();

  if (fFormat == "mesh")
  {
    if(fName.empty()){
      fName += simName + ".mesh";
    }
    aInput = shared_ptr<SimInput>(new SimInputMESH(fName, configNode, infoNode));
  }
  else if (fFormat == "cdb")
  {
    if(fName.empty()){
      fName += simName + ".cdb";
    }
    aInput = shared_ptr<SimInput> (new SimInputCDB(fName, configNode, infoNode));
  }
  else if (fFormat == "hdf5")
  {
    if(fName.empty()){
      fName += simName + ".h5";
    }
    aInput = shared_ptr<SimInput>(new SimInputHDF5(fName, configNode, infoNode));
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
    if(fName.empty()){
      fName += simName + ".msh";
    }
    aInput = shared_ptr<SimInput> (new SimInputGmsh(fName, configNode, infoNode));
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
    if(fName.empty()){
      fName += simName + ".unv";
    }
    aInput = shared_ptr<SimInput> (new SimInputUnv(fName, configNode, infoNode));
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
#endif
  }
  else if (fFormat == "internal")
  {
    if(fName.empty()){
      fName += simName + ".mesh";
    }
    aInput = shared_ptr<SimInput>(new InternalMesh(fName, configNode, infoNode));
  }
  else if (fFormat == "python")
  {
#ifdef USE_EMBEDDED_PYTHON
    aInput = shared_ptr<SimInput> (new SimInputPython(fName, configNode, infoNode));
#else
    EXCEPTION( "Compile with USE_EMBEDDED_PYTHON for python mesh reader." );
#endif
  }
  else
  {
    EXCEPTION( "Wrong format for input file. Please, check your data!" );
  }

  return aInput;

}

shared_ptr<SimOutput> DefineInOutFiles::CreateSingleOutputFileObject(string fName,
                                                                     PtrParamNode configNode,
                                                                     PtrParamNode infoNode,
                                                                     bool isRestart){
  shared_ptr<SimOutput> aOutput;
  string fFormat = configNode->GetName();
  if (fFormat == "unv")
  {
    aOutput = shared_ptr<SimOutput> (new SimOutputUnv(fName, configNode, infoNode, isRestart));
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
    aOutput =  shared_ptr<SimOutput> (new SimOutputGmsh(fName, configNode, infoNode, isRestart));
  }

  if (fFormat == "gmshParsed")
  {
    aOutput =  shared_ptr<SimOutput> (new SimOutputParsed(fName, configNode, infoNode, isRestart));
  }

  if (fFormat == "hdf5")
  {
    aOutput = shared_ptr<SimOutput>(new SimOutputHDF5(fName,configNode,infoNode, isRestart));
  }

  if (fFormat == "rst")
  {
    // was SimOutputRST(), but this is no more in the code?!
    EXCEPTION( "No support for ANSYS RST output file format." );
  }

  if (fFormat == "text" || fFormat == "csv")
  {
    aOutput = shared_ptr<SimOutput> (new SimOutputText(fName, configNode, infoNode, isRestart));
  }

  if (fFormat == "info")
  {
    aOutput = shared_ptr<SimOutput> (new SimOutputInfo(configNode,infoNode, isRestart));
  }

  if (fFormat == "streaming")
  {
#ifdef USE_STREAMING
    aOutput = shared_ptr<SimOutput> (new SimOutputStreaming(configNode, infoNode, isRestart));
#else
   throw Exception("not compiled with USE_STREAMING");
#endif
  }

  return aOutput;
}


} // end of namespace
