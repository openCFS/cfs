// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// Include headers which define what types
// of in/output files CFS++ should support
#include <def_use_mesh.hh>
#include <def_use_gidpost.hh>
#include <def_use_hdf5.hh>
#include <def_use_gmsh.hh>
#include <def_use_gmv.hh>
#include <def_use_unv.hh>
#include <def_use_ansysrst.hh>
#include <def_use_comsol.hh>

#include <boost/tokenizer.hpp>
#include <boost/filesystem/convenience.hpp>
namespace fs = boost::filesystem;

#include "General/Environment.hh"
#include "DataInOut/SimInput.hh"
#include "DataInOut/SimOutput.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/Xerces.hh"

#ifdef USE_MESH
#include "DataInOut/SimInOut/AnsysFile/SimInputMESH.hh"
#endif

#ifdef USE_GMSH
#include "DataInOut/SimInOut/gmsh/SimInputGmsh.hh"
#include "DataInOut/SimInOut/gmsh/SimOutputGmsh.hh"
#endif

#ifdef USE_GMV
#include "DataInOut/SimInOut/gmv/SimInputGMV.hh"
#include "DataInOut/SimInOut/gmv/SimOutGMV.hh"
#endif

#ifdef USE_HDF5
#include "DataInOut/SimInOut/hdf5/SimInputHDF5.hh"
#include "DataInOut/SimInOut/hdf5/SimOutputHDF5.hh"

#include "DataInOut/SimInOut/xdmf/SimOutputXDMF.hh"
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

#include "DataInOut/SimInOut/TextOutput/TextSimOutput.hh"


#include "ParamsInit.hh"
#include "HelperFuncs.hh"
#include "WVT.hh"

namespace CFSTool {

  inline void setParamNode(PtrParamNode paramNode, std::string name, std::string value,
      ParamNodeList* children = NULL);

  // Here we print out the hard-coded node and element limit in response
  // to the --printLimits command line option.
  void PrintLimits(const std::string& type) 
  {
    bool typeCMake = (type == "cmake");

    std::cout << std::endl;
    
    if(typeCMake) { std::cout << "SET("; }
    std::cout << "SOFT_NODE_LIMIT " << SOFT_NODE_LIMIT;
    if(typeCMake) { std::cout << ")"; }
    std::cout << std::endl;

    if(typeCMake) { std::cout << "SET("; }
    std::cout << "HARD_NODE_LIMIT " << HARD_NODE_LIMIT;
    if(typeCMake) { std::cout << ")"; }
    std::cout << std::endl;

    if(typeCMake) { std::cout << "SET("; }
    std::cout << "SOFT_ELEM_LIMIT " << SOFT_ELEM_LIMIT;
    if(typeCMake) { std::cout << ")"; }
    std::cout << std::endl;

    if(typeCMake) { std::cout << "SET("; }
    std::cout << "HARD_ELEM_LIMIT " << HARD_ELEM_LIMIT;
    if(typeCMake) { std::cout << ")"; }
    std::cout << std::endl;

    std::cout << std::endl;
  }

  shared_ptr<SimInput> GetInputReader(const std::string& fileName,
                                      const PtrParamNode& param,
                                      const PtrParamNode& info)
  {
    // determine suffix of fileName
    shared_ptr<SimInput> reader;

    if (fileName.find( ".refelem") == std::string::npos &&
        !fs::exists( fileName ))
      EXCEPTION( "\nFile '" << fileName << "' does not exist!");

    if( fileName.find( ".mesh") != std::string::npos ) {
#ifdef USE_MESH
      reader = shared_ptr<SimInput>(new SimInputMESH( fileName, PtrParamNode(), info ) );
#else
      EXCEPTION( "No support for MESH input file format." );
#endif
    } else if( fileName.find( ".h5") != std::string::npos ) {
#ifdef USE_HDF5
      reader = shared_ptr<SimInput>(new SimInputHDF5(fileName, param, info));
#else
      EXCEPTION( "No support for HDF5 input file format." );
#endif
    } else if( fileName.find( ".refelem") != std::string::npos ) {
      reader = shared_ptr<SimInput>(new SimInputRefElems(fileName, param, info));
    } else if( fileName.find( ".msh") != std::string::npos ) {
#ifdef USE_GMSH
      PtrParamNode gmshNode(new ParamNode (ParamNode::EX, ParamNode::ELEMENT));
      reader = shared_ptr<SimInput>(new SimInputGmsh(fileName, gmshNode, info) );
#else  
      EXCEPTION( "No support for Gmsh input file format." );
#endif
    } else if( fileName.find( ".mphtxt") != std::string::npos ) {
#ifdef USE_COMSOL
      PtrParamNode comsolNode(new ParamNode (ParamNode::EX, ParamNode::ELEMENT));
      reader = shared_ptr<SimInput>(new SimInputMPHTXT(fileName, comsolNode, info) );
#else  
      EXCEPTION( "No support for Comsol .mphtxt input file format." );
#endif
    } else if( fileName.find( ".gmv") != std::string::npos ) {
#ifdef USE_GMV
      reader = shared_ptr<SimInput>(new SimInputGMV(fileName, param, info));
#else
      EXCEPTION( "No support for GMV input file format." );
#endif
    } else if( fileName.find( ".unv") != std::string::npos ||
        fileName.find( ".unverg") != std::string::npos ||
        fileName.find( ".unvref") != std::string::npos ) {
#ifdef USE_UNV
      reader = shared_ptr<SimInput>(new SimInputUnv( fileName, param, info ));
#else
      EXCEPTION( "No support for UNV input file format." );
#endif
    } else {
      EXCEPTION( "Found not suitable reader for file '" << fileName
          << "'" );
    }

    return reader;
  }

  shared_ptr<SimOutput> GetOutputWriter( const std::string& fileName,
                                         const PtrParamNode& param,
                                         const PtrParamNode& info ) {
    // determine suffix for fileName
    shared_ptr<SimOutput> writer;
    std::string baseName;
    
    // we do not have restart functionality in the cfstool
    bool restart = false;

    if( fileName.find( ".post") != std::string::npos ) {
#ifdef USE_GIDPOST
      baseName = std::string(fileName, 0, fileName.find(".post"));
      PtrParamNode gidNode( new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
      PtrParamNode binary (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
      binary->SetName( "binaryFormat");
      if( fileName.find( ".bin") != std::string::npos ) {
        binary->SetValue( "yes");
      } else {
        binary->SetValue( "false");
      }
      gidNode->AddChildNode( binary);
      writer = shared_ptr<SimOutput>( new SimOutputGiD( baseName, gidNode, 
                                                        info, restart ) );
#else
      EXCEPTION( "No support for GiD output file format." );
#endif
    } else if( fileName.find( ".gmv") != std::string::npos ) {
#ifdef USE_GMV
      baseName = std::string(fileName, 0, fileName.find(".gmv"));
      PtrParamNode gmvNode(new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
      PtrParamNode binary (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));

      binary->SetName("binaryFormat");
      binary->SetValue( "yes" );
      PtrParamNode fixedGrid (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
      fixedGrid->SetName("fixedGrid");
      fixedGrid->SetValue( "yes" );
      gmvNode->AddChildNode( binary);
      gmvNode->AddChildNode( fixedGrid );
      writer = shared_ptr<SimOutput>( new SimOutputGMV( baseName, gmvNode, 
                                                        info, restart ) );
#else
      EXCEPTION( "No support for GMV output file format." );
#endif
    } else if( fileName.find( ".msh") != std::string::npos ) {
#ifdef USE_GMSH
      baseName = std::string(fileName, 0, fileName.find(".msh"));
      PtrParamNode gmshNode(new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
      PtrParamNode binary (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
      binary->SetName("binaryFormat");
      binary->SetValue( "yes" );
      PtrParamNode bigEndian (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
      bigEndian->SetName("endianness");
      bigEndian->SetValue( "big" );
      gmshNode->AddChildNode(binary);
      gmshNode->AddChildNode(bigEndian);
      writer = shared_ptr<SimOutput>( new SimOutputGmsh( baseName, gmshNode, 
                                                         info, restart ) );
#else 
      EXCEPTION( "No support for GMsh output file format." );
#endif
    } else if(fileName.find( ".h5") != std::string::npos) {
#ifdef USE_HDF5
      baseName = std::string(fileName, 0, fileName.find(".h5"));
      PtrParamNode h5Node (new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
      PtrParamNode eFiles (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
      eFiles->SetName("externalFiles");
      eFiles->SetValue( "false" );
      h5Node->AddChildNode(eFiles);
      writer =  shared_ptr<SimOutput>( new SimOutputHDF5( baseName, h5Node, 
                                                          info, restart ) );
#else
      EXCEPTION( "No support for HDF5 output file format." );
#endif
    } else if(fileName.find( ".xmf") != std::string::npos) {
#ifdef USE_HDF5
      baseName = std::string(fileName, 0, fileName.find(".xmf"));
      PtrParamNode h5Node (new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
      PtrParamNode eFiles (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
      eFiles->SetName("externalFiles");
      eFiles->SetValue( "false" );
      h5Node->AddChildNode(eFiles);
      writer =  shared_ptr<SimOutput>( new SimOutputXDMF( baseName, h5Node, 
                                                          info, restart ) );
#else
      EXCEPTION( "No support for HDF5 output file format. Cannot write XDMF files." );
#endif
    } else if(fileName.find( ".rst") != std::string::npos) {
#ifdef USE_ANSYSRST
      baseName = std::string(fileName, 0, fileName.find(".rst"));
      PtrParamNode rstNode (new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
      writer =  shared_ptr<SimOutput>( new SimOutputRST( baseName, rstNode ) );
#else
      EXCEPTION( "No support for ANSYS .rst output file format." );
#endif
    } else if(fileName.find( ".unv") != std::string::npos) {
#ifdef USE_UNV
      baseName = std::string(fileName, 0, fileName.find(".unv"));
      PtrParamNode unvNode (new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
      if(fileName.find( ".unverg") != std::string::npos) 
      {
        PtrParamNode flavor (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
        flavor->SetName("flavor");
        flavor->SetValue( "CAPA" );
        unvNode->AddChildNode(flavor);
      }
      writer =  shared_ptr<SimOutput>( new SimOutputUnv( baseName, unvNode, 
                                                         info, restart ) );
#else
      EXCEPTION( "No support for IDEAS universal output file format." );
#endif
    } else {
      EXCEPTION( "Output format not supported!" );
    }
    // create output writer
    return writer;
  }

  PtrParamNode GetParamNodeFromHDF5( shared_ptr<SimInput>& inputHDF5,
                                     const std::string& xmlFile ) {

    // read XML and material file from lateral mode file
    SimInputHDF5* hdf5Reader = dynamic_cast<SimInputHDF5*>(inputHDF5.get());
    std::string xml;
    
    try { 
      hdf5Reader->DB_Init();

      if(xmlFile == "ParameterFile") 
      {
        hdf5Reader->DB_GetParamFileContent( xml );
      }
      else 
      {
        if(xmlFile == "ParameterFile") 
        {
          hdf5Reader->DB_GetParamFileContent( xml );
        }
        else 
        {
          hdf5Reader->DB_Close();
          EXCEPTION("XML file " << xmlFile << " not found in HDF5 database group.");
        }
      }
      hdf5Reader->DB_Close();
    } catch (Exception& ex) 
    {
      hdf5Reader->ReadStringFromUserData(xmlFile, xml);
    }

    PtrParamNode pNode;
    
    CoupledField::Xerces xerces;
    xerces.SetString(xml);
    pNode = xerces.CreateParamNodeInstance();

    return pNode;
  }


  double RadPhase( const Complex& c ) {
    return std::atan2(c.imag() ,c.real() );
  }

  bool CheckReaderCapabilities(const std::string& readerDescription,
                               shared_ptr<SimInput>& reader,
                               std::string& result) 
  {
    std::stringstream sstr;
    bool ret = false;
    std::string fileName = reader->GetFileName();
    
    ret = std::find( reader->GetCapabilities().begin(),
                      reader->GetCapabilities().end(),
                      SimInput::MESH_RESULTS )
      != reader->GetCapabilities().end();

    sstr << "Reader for " << readerDescription << " with input file '"
         << fileName << "' " <<
      (ret ? "can provide results." :  "CANNOT provide results!");

    result = sstr.str();

    return ret;
  }

  /** Initialize static Enums.
   * todo: do better once - Fabian */
  void InitEnums()
  {
    SetEnvironmentEnums();
    BasePDE::SetEnums();
    EntityList::SetEnums();
  }

  void SetFreeCoord(const PtrParamNode& param,
                    std::string coordSysId,
                    std::string node_name)
  {

    std::vector<std::string> freeCoord;
    typedef boost::tokenizer< boost::char_separator<char> > Tok;
    boost::char_separator<char> sep(";| ");

    // Initialize vector with output fields
    Tok tokenizer(param->Get("freeCoord")->As<std::string>(), sep);
    std::copy(tokenizer.begin(), tokenizer.end(),
              std::back_inserter(freeCoord));
    if (freeCoord.size() != 4)
    {
      EXCEPTION("Not enought arguments " << freeCoord.size() \
          << ". Need 4 arguments (comp, start, stop, inc')");
    }

    PtrParamNode parent;

    ParamNodeList childVec;
    childVec.Push_back(PtrParamNode(new ParamNode())); // comp
    childVec.Push_back(PtrParamNode(new ParamNode())); // start
    childVec.Push_back(PtrParamNode(new ParamNode())); // stop
    childVec.Push_back(PtrParamNode(new ParamNode())); // inc

    // create a tree (Paramnode):
    // create leafes (comp, start, stop, inc)
    setParamNode(childVec[0], "comp",  freeCoord[0]);
    setParamNode(childVec[1], "start", freeCoord[1]);
    setParamNode(childVec[2], "stop",  freeCoord[2]);
    setParamNode(childVec[3], "inc",   freeCoord[3]);
    

    // create node (freeCoord) with comp,start,stop,inc as children
    parent = PtrParamNode(new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
    setParamNode(parent, "freeCoord", "", &childVec);
    
    // create node (list) with freeCoord and coordSysId as children
    childVec.Clear();
    childVec.Push_back(PtrParamNode( new ParamNode())); // coordSysId
    childVec.Push_back(parent);
    setParamNode(childVec[0], "coordSysId", coordSysId);
    childVec.Push_back(PtrParamNode( new ParamNode())); // gridId
    childVec.Push_back(parent);
    setParamNode(childVec[1], "gridId", "default");
    parent = PtrParamNode( new ParamNode()); //list
    setParamNode(parent, "list", "", &childVec);

    // create leaf (name)
    childVec.Clear();
    childVec.Push_back(parent);
    parent = PtrParamNode( new ParamNode()); //nodes name=""
    setParamNode(parent, "name", node_name);
    childVec.Push_back(parent);


    // create node (nodes) with name and list as childer
    parent = PtrParamNode( new ParamNode()); //nodes
    setParamNode(parent, "nodes", "", &childVec);

    // create node (nodeList) with nodes as child
    childVec.Clear();
    childVec.Push_back(parent);
    parent = PtrParamNode( new ParamNode()); //nodeList
    setParamNode(parent, "nodeList", "", &childVec);

    // create node (domain) with nodeList as child
    childVec.Clear();
    childVec.Push_back(parent); //domain
    //parent = new ParamNode(false); //nodeList
    //setParamNode(param, "domain", "", &childVec);
    //childVec.Clear();

    // put the childs into param
    param->Get("domain",ParamNode::INSERT )->GetChildren().Push_back(parent);
  }

  inline void setParamNode(PtrParamNode paramNode, std::string name, std::string value,
      ParamNodeList* children  )
  {
    paramNode->SetName(name);
    if (name != "")
    {
      paramNode->SetValue(value);
    }
    if (children != NULL)
    {
      ParamNodeList &childTmp = paramNode->GetChildren();
      childTmp = *children;
    }
  }



} //Namespace CFSTool

