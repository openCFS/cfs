// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <def_use_xerces.hh>

#ifdef USE_XERCES

#include <string>
#include <fstream>

#include <xercesc/util/PlatformUtils.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/algorithm/string.hpp>

#include "General/exception.hh"
#include "DataInOut/ParamHandling/Xerces.hh"


// we want to use the Xerces-C++ namespace
using namespace xercesc;
namespace fs = boost::filesystem;

namespace CoupledField
{

  Xerces::Xerces(const std::string& file, const std::string& schema)
  {
    parser_  = NULL;
    root_    = NULL;

    // create canonical path from native-representation of the
    // file and the schema path
    fs::path filePath = fs::path( file, fs::native );
    fs::path schemaPath = fs::path( schema, fs::native );

    if(!fs::exists(filePath))
        EXCEPTION("xml file " << file << " doesn't exist");

    if(!schema.empty() && !fs::exists(schemaPath))
        EXCEPTION("schema file " << schema << " doesn't exist");

    this->file_   = file;
    this->schema_ = schema;
  }

  /** this cannot be done in the constructor as the error handler takes "this" */
  void Xerces::Parse()
  {
     // Initialise the XML4C2 system
    try
    {
       XMLPlatformUtils::Initialize();
    }
    catch(const XMLException &event )
    {
        EXCEPTION("Error initializing xerces-c"
                  << XMLString::transcode(event.getMessage()));
    }

    parser_ = new XercesDOMParser();
    // skip whitespaces
    parser_->setIncludeIgnorableWhitespace(false);
    // do nice XML
    parser_->setDoNamespaces(true);

    //  Check if validation is desired and a schema file was provided
    if(!schema_.empty())
    {
      parser_->setValidationScheme(XercesDOMParser::Val_Always);
      parser_->setDoSchema(true);
      parser_->setValidationSchemaFullChecking(true);
      std::string completeSchema;
      completeSchema = "http://www.cfs++.org ";
      completeSchema += schema_;
      parser_->setExternalSchemaLocation(completeSchema.c_str());
    }
    else
    {
      parser_->setDoSchema(false);
      parser_->setValidationSchemaFullChecking(false);
    }

    // Have not yet understood what an entity reference node is,
    // but it seems we do not need them
    parser_->setCreateEntityReferenceNodes(false);

    // Attach our own error handler to the parser
    parser_->setErrorHandler(new EventHandler(this));

    try
    {
       // Parse and validate the XML file. This will generate the DOM tree.
       // Catch all exceptions that the parser could not pass to our error
       // handler.
       parser_->parse(file_.c_str());
    }
    catch(const XMLException &event)
    {
        EXCEPTION("Error parsing '" << file_ << "' -> '"
                  << event.getMessage() << "'");
    }
    catch(const DOMException &event )
    {
        const unsigned int maxChars = 2047;
        XMLCh errText[maxChars + 1];

        std::stringstream ss;
        ss << "DOM error in '" << file_ << "': DOMException.code = " << event.code;

        if(DOMImplementation::loadDOMExceptionMsg(event.code, errText, maxChars))
          ss << " Message is: " << errText;

        EXCEPTION(ss.str());
    }

    // Obtain and validate root element of document tree
    DOMDocument* doc = parser_->getDocument();
    DOMNodeList* children = doc->getChildNodes();

    // some final checking, cannot imagine a problem here
    if(children->getLength() != 1)
        EXCEPTION("document root has " << children->getLength()
                  << " childs, expceted 1");

    if(children->item(0)->getNodeType() != DOMNode::ELEMENT_NODE)
        EXCEPTION("root node type is " <<  children->item(0)->getNodeType());

    root_ = children->item(0);

  }

  Xerces::~Xerces()
  {

     if(parser_ != NULL && parser_->getErrorHandler() != NULL)
     {
        delete parser_->getErrorHandler();
     }

     if(parser_ != NULL)
     {
        delete parser_;
        parser_ = NULL;
     }



    // Shutdown platform dependend utilities
    XMLPlatformUtils::Terminate();

  }


  PtrParamNode Xerces::CreateParamNodeInstance()
  {
     // read the file, this cannot be done in the constructor as the error handler takes this object
     Parse();

     PtrParamNode out (new ParamNode(ParamNode::EX, ParamNode::ELEMENT ) );
     Fill(root_, out);
     return out;
  }

  void Xerces::Fill(DOMNode* node, PtrParamNode parent)
  {
    // determine if this is a valid node
    // std::cout << "node " << XMLString::transcode(node->getNodeName()) << " has type " << node->getNodeType();
    // std::cout << " value = " << (node->getNodeValue() != NULL ? XMLString::transcode(node->getNodeValue()) : "null") << std::endl;

    switch(node->getNodeType())
    {

    case DOMNode::TEXT_NODE:
    {
      // if we are a text node, we "are" the value of our parent.
      std::string temp(XMLString::transcode(node->getNodeValue()));
      boost::trim(temp);
      parent->SetValue(temp);
      return; // nothing else to do, we don not create a new ParamNode
    }
    case DOMNode::ELEMENT_NODE:
    case DOMNode::ATTRIBUTE_NODE:
      // this is the typical situation, we create a new ParamNode
      break;

    default:
      // comment type or something alike, don't do anything
      return;
    }
    // normally we create here a new element and add it to parent.
    // This is not the case when node is root_, then we set the properties of parent directly
    PtrParamNode pn;
    if(node != root_)
    {
      // create a new param node and set it as a new child at the father
      ParamNode::NodeType type = ParamNode::ELEMENT;
      if ( node->getNodeType() == DOMNode::ATTRIBUTE_NODE  ) 
        type = ParamNode::ATTRIBUTE;
      PtrParamNode newNode = 
          PtrParamNode(new ParamNode(ParamNode::EX, type));
      parent->AddChildNode( newNode );
      pn = newNode;
      //parent->GetChildren().Push_back(new ParamNode();
      // we work with the this just added element - here we avoid any
      // potential copy constructor issues
      pn = parent->GetChildren().Last();
    }
    else
    {
      // no new child created but we modify the parent directly
      pn = parent;
    }

    pn->SetName( XMLString::transcode( node->getNodeName() ) );

    // The value of an attribute or simple element is set by the text node children

    // first process attributes - map is NULL this is not DOMElement
    DOMNamedNodeMap* map = node->getAttributes();
    const unsigned int mlen = map == NULL ? 0 : map->getLength();
    for(unsigned int i = 0; i < mlen; ++i)
    {
      // recursive call
      Fill(map->item(i), pn);
    }
    
    // process childs (if there are any)
    DOMNode *child = node->getFirstChild();
    while(child != NULL)
    {
      Fill(child, pn);
      child = child->getNextSibling();
    }
  }

  Xerces::EventHandler::EventHandler(const Xerces* xerces)
  {
     this->xerces_ = xerces;
  }


  void Xerces::EventHandler::warning(const SAXParseException &event )
  {
    std::stringstream os;
    os << "WARN parsing the xml file'" << xerces_->file_ << "' in line "
       << event.getLineNumber() << ", column " << event.getColumnNumber() << std::endl
       << "-> '" << XMLString::transcode(event.getMessage()) << "'" << std::endl
       << " schema: '" << (!xerces_->schema_.empty() ? xerces_->schema_ : "<no-schema>") << "'";
    // killme use the new log stuff from Andi
    std::cerr << os.str() << std::endl;
  }


  void Xerces::EventHandler::error(const SAXParseException &event)
  {
    EXCEPTION("Error parsing the xml file' " << xerces_->file_
              << "' in line " << event.getLineNumber() << ", column "
              << event.getColumnNumber() << std::endl << "-> '"
              << XMLString::transcode(event.getMessage()) << "'"
              << std::endl << " schema: '"
              << (!xerces_->schema_.empty() ? xerces_->schema_ : "<no-schema>") << "'");

  }


  void Xerces::EventHandler::fatalError(const SAXParseException &event)
  {
     error(event);
  }

  void Xerces::EventHandler::resetErrors()
  {
  }



} // end of namespace

#endif
