// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <fstream>
#include <iostream>
#include <cstdio>

#include <string>

// Include headers which define what types
// of in/output files CFS++ should support
#include <def_use_mesh.hh>
#include <def_use_gidpost.hh>
#include <def_use_hdf5.hh>
#include <def_use_gmv.hh>
#include <def_use_gmsh.hh>
#include <def_use_unv.hh>
#include <def_use_ansysrst.hh>
#include <def_use_scripting.hh>
#include <def_use_tcl.hh>
#include <def_use_python.hh>

#include "definefiles.hh"
#include "DataInOut/WriteInfo.hh"

#ifdef USE_MESH
#include "DataInOut/SimInOut/AnsysFile/simInputMESH.hh"
#include "DataInOut/SimInOut/internalMesh/internalMesh.hh"
#endif

#ifdef USE_GMV_INPUT
#include "DataInOut/SimInOut/gmv/simInputGMV.hh"
#endif

#ifdef USE_GMV_OUTPUT
#include "DataInOut/SimInOut/gmv/simOutGMV.hh"
#endif

#ifdef USE_GMSH
#include "DataInOut/SimInOut/gmsh/simInputGmsh.hh"
#include "DataInOut/SimInOut/gmsh/simOutputGmsh.hh"
#endif

#ifdef USE_HDF5
// HDF5 readers and writers
#include "DataInOut/SimInOut/hdf5/simInputHDF5.hh"
#include "DataInOut/SimInOut/hdf5/simOutputHDF5.hh"

// XDMF writer
#include "DataInOut/SimInOut/xdmf/simOutputXDMF.hh"
#endif

#ifdef USE_GIDPOST
#include "DataInOut/SimInOut/GiD/simOutGiD.hh"
#endif

#ifdef USE_UNV
#include "DataInOut/SimInOut/Unverg/simInputUnv.hh"
#include "DataInOut/SimInOut/Unverg/simOutUnv.hh"
#endif

#ifdef USE_ANSYSRST
#include "DataInOut/SimInOut/AnsysRST/simOutputRST.hh"
#endif

#include "DataInOut/SimInOut/TextOutput/textSimOutput.hh"
#include "DataInOut/SimInOut/InfoResultOutput/SimOutputInfo.hh"
#include "DataInOut/SimInOut/Streaming/SimOutputStreaming.hh"

#include "DataInOut/XMLMaterialHandler.hh"

#include "DataInOut/programOptions.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

#ifdef USE_SCRIPTING
#include "DataInOut/Scripting/cfsmessenger.hh"
#endif

#ifdef USE_SCRIPTING_TCL
#include "DataInOut/Scripting/tcl/tcl-messenger.hh"
#endif

#ifdef USE_SCRIPTING_PYTHON
#include "DataInOut/Scripting/python/py-messenger.hh"
#endif

// Maximal length of the trailing postfix of an auxilliary name,
// i.e. the length of the extension after the basename
#define MAXPOSTFIX 15

namespace CoupledField
{

  // ===============
  //   Constructor
  // ===============
  DefineInOutFiles::DefineInOutFiles() {
    // Initialise internal pointers
    simInput_ = NULL;
    ptMaterialHandler_ = NULL;
  }
  
  // ==============
  //   Destructor
  // ==============
  DefineInOutFiles::~DefineInOutFiles()
  {
    
 #ifdef MEMTRACE
    delete memtrace;
    memtrace = NULL;
 #endif
    
    // Delete pointer to OLAS report file
    delete cla;
    cla = NULL;
    
    delete simInput_;
    simInput_ = NULL;
    
    delete ptMaterialHandler_;
    ptMaterialHandler_ = NULL;
    
  }
  
// ==============================
//   Generate mesh file pointer
// ==============================
void DefineInOutFiles::CreateSimInputFiles(std::map<std::string, shared_ptr<
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
  PtrParamNode inputNode = param->Get("fileFormats") ->Get("input",ParamNode::INSERT);
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
    
    inFiles[actId] = shared_ptr<SimInput> (new SimInputMESH(meshFile, PtrParamNode() ));
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

    if (informat == "mesh")
    {
#ifdef USE_MESH
      if (meshFile.empty())
        meshFile = simName + ".mesh";
      inFiles[actId] = shared_ptr<SimInput> (
          new SimInputMESH(meshFile, actNode));
#else
      EXCEPTION( "No support for MESH input file format." );
#endif // USE_MESH
    }
    else if (informat == "hdf5")
    {
#ifdef USE_HDF5
      if (meshFile.empty())
        meshFile = simName + ".h5";
      inFiles[actId] = shared_ptr<SimInput> (
          new SimInputHDF5(meshFile, actNode));
#else
      EXCEPTION( "No support for HDF5 input file format." );
#endif // USE_HDF5
    }
    else if (informat == "gmv")
    {
#ifdef USE_GMV_INPUT
      if(meshFile.empty())
      meshFile = simName + ".gmv";
      inFiles[actId] = shared_ptr<SimInput>(new SimInputGMV(meshFile,
              actNode));
#else
      EXCEPTION( "No support for GMV input file format." );
#endif // USE_GMV_INPUT
    }
    else if (informat == "gmsh")
    {
#ifdef USE_GMSH
      if (meshFile == "")
        meshFile = simName + ".msh";
      inFiles[actId] = shared_ptr<SimInput> (
          new SimInputGmsh(meshFile, actNode));
#else
      EXCEPTION( "No support for GMSH input file format." );
#endif // USE_GMSH
    }
    else if (informat == "unv")
    {
#ifdef USE_UNV
      if (meshFile.empty())
        meshFile = simName + ".unv";
      inFiles[actId]
          = shared_ptr<SimInput> (new SimInputUnv(meshFile, actNode));
#else
      EXCEPTION( "No support for UNV input file format." );
#endif // USE_UNV
    }
    else if (informat == "internal")
    {
#ifdef USE_MESH
      if (meshFile.empty())
        meshFile = simName + ".mesh";
      inFiles[actId] = 
          shared_ptr<SimInput> (new InternalMesh(meshFile, actNode));
#else
      EXCEPTION( "No support for internalMesh input file format." );
#endif // USE_MESH
    }
    else
    {
      EXCEPTION( "Wrong format for input file. Please, check your data!" );
    }

    // relate gridId with input reader id
    gridInputs[actGridId].Push_back(inFiles[actId]);
  }
}

// ================================
//   Generate output file pointer
// ================================
void DefineInOutFiles::CreateSimOutputFiles(std::map<std::string, shared_ptr<
    SimOutput> >& out)
{

  // resest map
  out.clear();

  std::string simName = progOpts->GetSimName();

  // get list of output formats
  PtrParamNode outNode = param->Get("fileFormats")->Get("output", ParamNode::PASS);

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
  
  // The HDF5 writer needs to be known to the XDMF writer.
  shared_ptr<SimOutput> hdf5Writer;
  
  // iterate over all found files
  for (UInt i = 0; i < formatNodes.GetSize(); i++)
  {
    // fetch format and id of output class
    PtrParamNode actNode = formatNodes[i];
    actFormat = actNode->GetName();
    actId = actNode->Get("id")->As<std::string>();

    if (actFormat == "unv")
    {
#ifdef USE_UNV
      out[actId] = shared_ptr<SimOutput> (new SimOutputUnv(simName, actNode));
#else
      EXCEPTION( "No support for UNV output file format." );
#endif
    }

    if (actFormat == "gid")
    {
#ifdef USE_GIDPOST
      out[actId] = shared_ptr<SimOutput> (new SimOutputGiD(simName, actNode));
      continue;
#else
      EXCEPTION( "No support for GiD output file format." );
#endif
    }

    if (actFormat == "gmsh")
    {
#ifdef USE_GMSH
      out[actId] = shared_ptr<SimOutput> (new SimOutputGmsh(simName, actNode));
      continue;
#else
      EXCEPTION( "No support for Gmsh output file format." );
#endif
    }

    if (actFormat == "gmv")
    {
#ifdef USE_GMV_OUTPUT
      out[actId] = shared_ptr<SimOutput> (new SimOutputGMV(simName, actNode));
#else
      EXCEPTION( "No support for GMV output file format." );
#endif
    }

    if (actFormat == "hdf5")
    {
#ifdef USE_HDF5
      if(!hdf5Writer) 
      {        
        hdf5Writer.reset(new SimOutputHDF5(simName, actNode));
        out[actId] = hdf5Writer;

        std::cout << "Creating HDF5 writer '" << actId << "'..." << std::endl;
      }
#else
      EXCEPTION( "No support for HDF5 output file format." );
#endif
    }

    if (actFormat == "xdmf")
    {
#ifdef USE_HDF5
      if(!hdf5Writer) 
      {        
        hdf5Writer.reset(new SimOutputHDF5(simName, actNode));

        if(hdf5Id == "")
          hdf5Id = actId + "_hdf5";

        out[hdf5Id] = hdf5Writer;

        std::cout << "Creating HDF5/XDMF writer '" << hdf5Id << "'..." << std::endl;
      }
      
      SimOutputXDMF* simOutXDMF = new SimOutputXDMF(simName, actNode);
      if(simOutXDMF) 
      {
        out[actId] = shared_ptr<SimOutput> (simOutXDMF);
        simOutXDMF->SetHDF5Writer(dynamic_cast<SimOutputHDF5*>( hdf5Writer.get() ), false);
      }
        
#else
      EXCEPTION( "No support for HDF5 output file format.\n"
                 << "Therefore it does not make sense to write XDMF output." );
#endif
    }

    if (actFormat == "rst")
    {
#ifdef USE_ANSYSRST
      out[actId] =
      shared_ptr<SimOutput>( new SimOutputRST( simName, actNode ) );
#else
      EXCEPTION( "No support for ANSYS RST output file format." );
#endif
    }

    if (actFormat == "text")
    {
      out[actId] = shared_ptr<SimOutput> (new SimOutputText(simName, actNode));
    }
    
    if (actFormat == "info")
    {
      out[actId] = shared_ptr<SimOutput> (new SimOutputInfo(actNode));
    }

    if (actFormat == "streaming")
    {
      out[actId] = shared_ptr<SimOutput> (new SimOutputStreaming(actNode));
    }

  }
}

// ==================================
//   Generate material file handler
// ==================================
MaterialHandler *
DefineInOutFiles::CreateMaterialHandler()
{

  std::string fileName = "mat.dat";
  std::string format = "dat";

  // Determine filename and format
  PtrParamNode matNode = 
      param->Get("fileFormats")->Get("materialData", ParamNode::PASS);
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
    ptMaterialHandler_ = new XMLMaterialHandler(fileName);
  }
  else
  {
    EXCEPTION( "CreateMaterialHandler: Format '" << format
        << "' is not recognized!" );
  }
  return ptMaterialHandler_;

}

// ************
//   OpenFile
// ************
void DefineInOutFiles::OpenFile(AuxFileType fileType)
{

  std::string fileName;

  switch (fileType)
  {

  // ---------------
  //  MEMTRACE_FILE
  // ---------------
  case MEMTRACE_FILE:

    fileName = progOpts->GetSimName() + ".mem";

#ifdef MEMTRACE
    memtrace = new std::ofstream(fileName.c_str());
    if (memtrace == NULL)
    {
      EXCEPTION( "Failed to open file '" << fileName << "' for "
          << "writing memory trace information!" );
    }
    break;
#else
    EXCEPTION( "Executable was compiled without support for memory "
        << "tracing! Refusing to open trace file '"
        << fileName << "'" );
    break;
#endif

    // ------------
    //  OLAS_FILE
    // ------------
  case OLAS_FILE:

    fileName = progOpts->GetSimName() + ".las";
    cla = new std::ofstream(fileName.c_str());
    if (cla == NULL)
    {
      EXCEPTION( "Failed to open file '" << fileName << "' for "
          << "writing status reports of OLAS, the linear algebra "
          << "sub-system!" );
    }
    break;

  }
}

CFSMessenger* DefineInOutFiles::CreateScriptMessenger(
    const std::string& fileName)
{

  CFSMessenger * messenger = NULL;

#ifdef USE_SCRIPTING
  // check filename, if it is not empty
  if (fileName == "")
  {
    messenger = new CFSMessenger();
  }

#ifdef USE_SCRIPTING_TCL
  else if (fileName.find(".tcl") != std::string::npos)
  {
    messenger = new TCL_CFSMessenger();
  }
#endif
#ifdef USE_SCRIPTING_PYTHON
  else if( fileName.find( ".py") != std::string::npos )
  {
    messenger = new PY_CFSMessenger();
  }
#endif
  else
  {
    EXCEPTION( "Could not determine script language of file '"
        << fileName << "'!" );
  }
#endif

  return messenger;

}

} // end of namespace
