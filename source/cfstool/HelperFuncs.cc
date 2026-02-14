// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

// Include headers which define what types
// of in/output files openCFS should support
#include <def_use_gidpost.hh>
#include <def_use_unv.hh>
#include <def_use_cgns.hh>

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/filesystem/fstream.hpp>
namespace fs = boost::filesystem;

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/xml_parser.hpp>
#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/property_tree/json_parser.hpp>

#include "General/Environment.hh"
#include "DataInOut/SimInput.hh"
#include "DataInOut/SimOutput.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ParamHandling/XmlReader.hh"

#include "DataInOut/SimInOut/AnsysCDB/SimInputCDB.hh"
#include "DataInOut/SimInOut/AnsysFile/SimInputMESH.hh"

#include "DataInOut/SimInOut/gmsh/SimInputGmsh.hh"
#include "DataInOut/SimInOut/gmsh/SimOutputGmsh.hh"

#include "DataInOut/SimInOut/hdf5/SimInputHDF5.hh"
#include "DataInOut/SimInOut/hdf5/SimOutputHDF5.hh"
#include "DataInOut/SimInOut/RefElems/SimInputRefElems.hh"

#ifdef USE_GIDPOST
#include "DataInOut/SimInOut/GiD/SimOutGiD.hh"
#endif

#include "DataInOut/SimInOut/Unverg/SimInputUnv.hh"
#include "DataInOut/SimInOut/Unverg/SimOutputUnv.hh"

#ifdef USE_CGNS
#include "DataInOut/SimInOut/CGNS/SimInputCGNS.hh"
#include "DataInOut/SimInOut/CGNS/SimOutputCGNS.hh"
#endif

#include "DataInOut/SimInOut/TextOutput/TextSimOutput.hh"


#include "ParamsInit.hh"
#include "HelperFuncs.hh"
#include "WVT.hh"

using std::string;

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

  shared_ptr<SimInput> GetInputReader(const std::string& fileName, const PtrParamNode& param, const PtrParamNode& info)
  {
    // determine suffix of fileName
    shared_ptr<SimInput> reader;
    PtrParamNode inputNode = param->Get("fileFormats")->Get("input");
    PtrParamNode readerNode;


    if (fileName.find( ".refelem") == string::npos && !fs::exists( fileName ))
      EXCEPTION( "\nFile '" << fileName << "' does not exist!");

    if( fileName.find( ".mesh") != string::npos )
    {
      if(inputNode->Has("mesh"))
        readerNode = inputNode->Get("mesh");
      else
      {
        readerNode = PtrParamNode(new ParamNode());
        readerNode->SetName("mesh");
      }
        
      reader = shared_ptr<SimInput>(new SimInputMESH( fileName, readerNode, info ) );
    }
    else if(fileName.find( ".cdb") != string::npos || fileName.find( ".inp") != string::npos )
    {
      if(inputNode->Has("cdb"))
        readerNode = inputNode->Get("cdb");
      else
      {
        readerNode = PtrParamNode(new ParamNode());
        readerNode->SetName("cdb");
      }

      reader = shared_ptr<SimInput>(new SimInputCDB(fileName, readerNode, info));
    }
    else if(fileName.find(".h5") != string::npos || fileName.find(".cfs") != string::npos)
    {
      if(inputNode->Has("hdf5"))
        readerNode = inputNode->Get("hdf5");
      else
      {
        readerNode = PtrParamNode(new ParamNode());
        readerNode->SetName("hdf5");
      }
      reader = shared_ptr<SimInput>(new SimInputHDF5(fileName, readerNode, info));
    }
    else if( fileName.find(".refelem") != string::npos )
    {
      readerNode = PtrParamNode(new ParamNode());
      readerNode->SetName("refelem");
      reader = shared_ptr<SimInput>(new SimInputRefElems(fileName, param, info));
    }
    else if( fileName.find(".msh") != string::npos )
    {
      if(inputNode->Has("gmsh"))
        readerNode = inputNode->Get("gmsh");
      else
      {
        readerNode = PtrParamNode(new ParamNode());
        readerNode->SetName("gmsh");
      }
      reader = shared_ptr<SimInput>(new SimInputGmsh(fileName, readerNode, info) );
    }
    else if( fileName.find( ".cgns") != string::npos )
    {
    #ifdef USE_CGNS
      if(inputNode->Has("cgns"))
        readerNode = inputNode->Get("cgns");
      else
      {
        readerNode = PtrParamNode(new ParamNode());
        readerNode->SetName("cgns");
      }

      reader = shared_ptr<SimInput>(new SimInputCGNS(fileName, readerNode, info) );
    #else
      EXCEPTION( "No support for CGNS .cgns input file format." );
    #endif
    }
    else if(fileName.find( ".unv") != string::npos ||
            fileName.find( ".unverg") != string::npos ||
            fileName.find( ".unvref") != string::npos)
    {
    #ifdef USE_UNV
      if(inputNode->Has("unv"))
        readerNode = inputNode->Get("unv");
      else
      {
        readerNode = PtrParamNode(new ParamNode());
        readerNode->SetName("unv");

        PtrParamNode flavorNode(new ParamNode());
        readerNode->AddChildNode( flavorNode );
        flavorNode->SetName("flavor");
        flavorNode->SetValue( "ideas" );

        if(fileName.find( ".unverg") != string::npos ||
           fileName.find( ".unvref") != string::npos ) {
          flavorNode->SetValue( "capa" );
        }
      }

      reader = shared_ptr<SimInput>(new SimInputUnv( fileName, readerNode, info ));
    #else
      EXCEPTION( "No support for UNV input file format." );
    #endif
    }
    else
      EXCEPTION( "Found no suitable reader for file '" << fileName<< "'" );

    if(!inputNode->Has(readerNode->GetName()))
      inputNode->AddChildNode( readerNode );

    PtrParamNode fn (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
    fn->SetName("fileName");
    fn->SetValue( fileName );
    readerNode->AddChildNode(fn);

    return reader;
  }

  shared_ptr<SimOutput> GetOutputWriter( const std::string& fileName,
                                         const PtrParamNode& param,
                                         const PtrParamNode& info ) {
    // determine suffix for fileName
    shared_ptr<SimOutput> writer;
    std::string baseName;
    PtrParamNode outputNode = param->Get("fileFormats")->Get("output");
    PtrParamNode writerNode;
    
    // we do not have restart functionality in the cfstool
    bool restart = false;

    if( fileName.find( ".post") != string::npos ) {
#ifdef USE_GIDPOST
      baseName = std::string(fileName, 0, fileName.find(".post"));

      PtrParamNode binary (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
      binary->SetName( "binaryFormat");
      if( fileName.find( ".bin") != string::npos ) {
        binary->SetValue( "yes");
      } else {
        binary->SetValue( "false");
      }

      if(outputNode->Has("gid")) {
        writerNode = outputNode->Get("gid");
        if(writerNode->Has("binaryFormat")) {
          MergeParamNodes(binary, writerNode->Get("binaryFormat"));
        } else {
          writerNode->AddChildNode( binary);
        }        
      } else {
        writerNode = PtrParamNode(new ParamNode());
        writerNode->SetName("gid");
        writerNode->AddChildNode( binary);
      }

      writer = shared_ptr<SimOutput>( new SimOutputGiD( baseName, writerNode, 
                                                        info, restart ) );
#else
      EXCEPTION( "No support for GiD output file format." );
#endif
    }
    else if( fileName.find( ".msh") != string::npos )
    {
      baseName = std::string(fileName, 0, fileName.find(".msh"));
      PtrParamNode binary (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
      binary->SetName("binaryFormat");
      binary->SetValue( "yes" );
      PtrParamNode bigEndian (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
      bigEndian->SetName("endianness");
      bigEndian->SetValue( "big" );

      if(outputNode->Has("gmsh")) {
        writerNode = outputNode->Get("gmsh");
        if(!writerNode->Has("binaryFormat")) {
          writerNode->AddChildNode( binary );
        } 
        if(!writerNode->Has("endianness")) {
          writerNode->AddChildNode( bigEndian );
        } 
      } else {
        writerNode = PtrParamNode(new ParamNode());
        writerNode->SetName("gmsh");
        writerNode->AddChildNode( binary );
        writerNode->AddChildNode( bigEndian );
      }
      writer = shared_ptr<SimOutput>(new SimOutputGmsh( baseName, writerNode, info, restart));
    }
    else if(fileName.find(".h5") != string::npos || fileName.find(".cfs") != string::npos) {
      string ext = fileName.find(".h5") != string::npos ? ".h5" : ".cfs";
      baseName = std::string(fileName, 0, fileName.find(ext));
      PtrParamNode eFiles (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
      eFiles->SetName("externalFiles");
      eFiles->SetValue( "false" );

      if(outputNode->Has("hdf5")) {
        writerNode = outputNode->Get("hdf5");
        if(!writerNode->Has("externalFiles")) {
          writerNode->AddChildNode( eFiles );
        } 
      } else {
        writerNode = PtrParamNode(new ParamNode());
        writerNode->SetName("hdf5");
        writerNode->AddChildNode( eFiles );
      }
      writer = shared_ptr<SimOutput>(new SimOutputHDF5(baseName, writerNode, info, restart));
    } else if(fileName.find( ".unv") != string::npos) {
      baseName = std::string(fileName, 0, fileName.find(".unv"));

      if(outputNode->Has("unv")) {
        writerNode = outputNode->Get("unv");
      } else {
        writerNode = PtrParamNode(new ParamNode());
        writerNode->SetName("unv");
      }

      if(fileName.find( ".unverg") != string::npos)
      {
        PtrParamNode flavor (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
        flavor->SetName("flavor");
        flavor->SetValue( "CAPA" );
        
        if(writerNode->Has("flavor")) {
          MergeParamNodes(flavor, writerNode->Get("flavor"));
        } else {
          writerNode->AddChildNode( flavor );          
        }
      }

      writer = shared_ptr<SimOutput>(new SimOutputUnv(baseName, writerNode, info, restart));
    }
    else if(fileName.find( ".cgns") != string::npos) {
#ifdef USE_CGNS
      baseName = std::string(fileName, 0, fileName.find(".cgns"));

      if(outputNode->Has("cgns")) {
        writerNode = outputNode->Get("cgns");
      } else {
        writerNode = PtrParamNode(new ParamNode());
        writerNode->SetName("cgns");
      }

      writer =  shared_ptr<SimOutput>( new SimOutputCGNS( baseName, writerNode, 
                                                          info, restart ) );
#else
      EXCEPTION( "No support for CGNS output file format." );
#endif
    } else {
      EXCEPTION( "Output format not supported!" );
    }

    if(!outputNode->Has(writerNode->GetName())) {
      outputNode->AddChildNode( writerNode );
    }

    PtrParamNode fn (new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
    fn->SetName("fileName");
    fn->SetValue( fileName );
    writerNode->AddChildNode(fn);
    
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
        hdf5Reader->DB_Close();
        EXCEPTION("XML file " << xmlFile << " not found in HDF5 database group.");
      }
      hdf5Reader->DB_Close();
    } catch (Exception& ex) 
    {
      hdf5Reader->ReadStringFromUserData(xmlFile, xml);
    }

    PtrParamNode pNode = XmlReader::ParseString(xml);

    return pNode;
  }

  void WriteStringToHDF5UserData( shared_ptr<SimOutput>& outputHDF5,
                                  const std::string& dsetName,
                                  const std::string& str ) {
    SimOutputHDF5* hdf5Writer = dynamic_cast<SimOutputHDF5*>(outputHDF5.get());

    hdf5Writer->WriteStringToUserData( dsetName, str );
  }

  PtrParamNode ParamNodeFromPropertyTree( const std::string& propertyFile,
                                          const std::string& context ) {
    PtrParamNode pNode;
    std::stringstream sstr;
    // Create an empty property tree object
    using boost::property_tree::ptree;
    ptree pt;

    if( fs::exists(propertyFile) ) 
    {
      if( propertyFile.find(".xml") != string::npos) {
        try 
        {
          read_xml(propertyFile, pt);
        }  catch (boost::property_tree::xml_parser_error& xmlParseEx)
        {
          EXCEPTION("Could not parse XML file '" << propertyFile << "' for "
                    << context << "!" << std::endl
                    << "XML problem: " << xmlParseEx.what() << std::endl);
        }        
      } else if(propertyFile.find(".json") != string::npos) {
        try 
        {
          read_json(propertyFile, pt);
        }  catch (boost::property_tree::json_parser_error& jsonParseEx)
        {
          EXCEPTION("Could not parse JSON file '" << propertyFile << "' for "
                    << context << "!" << std::endl
                    << "JSON problem: " << jsonParseEx.what() << std::endl);
        }        
      }
    }
    else
    {
      std::string inputStr = propertyFile;
      boost::trim(inputStr);

      if(inputStr == "") 
      {
        return pNode;
      }
      
      try 
      {
        sstr << inputStr;
        read_xml(sstr, pt);
      } catch (boost::property_tree::xml_parser_error& xmlParseEx)
      {
        sstr.str(""); sstr.clear();
        sstr << inputStr;

        try 
        {
          read_json(sstr, pt);
        } catch (boost::property_tree::json_parser_error& jsonParseEx)
        {
          EXCEPTION("Could not parse input string neither as XML nor as JSON for "
                    << context << "!" << std::endl
                    << "XML problem: " << xmlParseEx.what() << std::endl
                    << "JSON problem: " << jsonParseEx.what() << std::endl);
        }        
      }

      sstr.str(""); sstr.clear();
    }
      
    // Write the property tree to a stream as XML.
    write_xml(sstr, pt);
    std::string str = sstr.str();

    // Create ParamNode tree from XML.
    pNode = XmlReader::ParseString(str);

    return pNode;
  }
  
  void MergeParamNodes(PtrParamNode src, PtrParamNode dest) 
  {
    if(src && dest && src->HasChildren()) 
    {
      ParamNodeList list = src->GetChildren();

      for(UInt node=0, nnodes=list.GetSize(); node < nnodes; node++) 
      {
        PtrParamNode newFormatNode = list[node];    
        PtrParamNode oldFormatNode;
        std::string format = newFormatNode->GetName();
        
        if( dest->Has(format) )
        {
          oldFormatNode = dest->Get(format);
        
          // Merge together informations
          ParamNodeList fList;
          if(newFormatNode->HasChildren()) fList = newFormatNode->GetChildren();
          for(UInt snode=0, nsnodes=fList.GetSize(); snode < nsnodes; snode++) 
          {
            std::string str = fList[snode]->GetName();
            
            if( oldFormatNode->Has(str) ) 
            {
              ParamNodeList oldList = oldFormatNode->GetChildren();
              for(UInt on=0, non=oldList.GetSize(); on < non; on++) 
              {
                if(oldList[on]->GetName() == str) 
                {
                  oldFormatNode->ReplaceChild(fList[snode], on);
                  break;
                }
              }
            }
            else
            {
              oldFormatNode->AddChildNode(fList[snode]);
            }
          }
        } 
        else
        {            
          dest->AddChildNode( list[node] );
        }
      }
    }
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

