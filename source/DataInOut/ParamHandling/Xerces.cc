#include <def_use_xerces.hh>

#include <string>
#include <fstream>

#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/TransService.hpp>
#include <xercesc/validators/datatype/InvalidDatatypeValueException.hpp>

#include <filesystem>
#include <boost/algorithm/string.hpp>
#include <boost/make_shared.hpp>

#include "General/Exception.hh"
#include "DataInOut/ParamHandling/Xerces.hh"


// we want to use the Xerces-C++ namespace
using namespace xercesc;
namespace fs = std::filesystem;

namespace CoupledField
{

  Xerces::Xerces(const std::string& schema, const std::string &url)
  {
    parser_  = NULL;
    root_    = NULL;
    buf_     = NULL;

    // create canonical path from native-representation of the
    // the schema path
    fs::path schemaPath = fs::absolute( fs::path( schema ) );

    if (!schema.empty() && !fs::exists(schemaPath)) {
        EXCEPTION("schema file " << schema << " doesn't exist");
    }

    schema_ = schema;
    schemaUrl_ = url;

    // Initialise the XML4C2 system
    try {
       XMLPlatformUtils::Initialize();
    }
    catch(const XMLException &event ) {
      EXCEPTION("Error initializing xerces-c" << Transcode(event.getMessage()));
    }
  }

  void Xerces::SetFile( const std::string& file ) 
  {
    // create canonical path from native-representation of the
    // file path
    fs::path filePath = fs::absolute( fs::path( file ) );

    if (!fs::exists(filePath)) {
      EXCEPTION("xml file " << file << " doesn't exist");
    }

    this->file_  = file;
    delete buf_;
    buf_ = NULL;
    
  }
  
  void Xerces::SetString( const std::string& text ) {
    this->buf_ = 
        new xercesc::MemBufInputSource((const XMLByte*)text.c_str(), 
                                       text.size(),"- memory xml -", false);
    this->file_ = "";
  }
  
  /** this cannot be done in the constructor as the error handler takes "this" */
  void Xerces::Parse()
  {

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
      std::string completeSchema = schemaUrl_ + " ";
      std::string schemaString = schema_;
      boost::replace_all( schemaString, " ", "%20");
      completeSchema += schemaString;
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
      // Parse and validate the XML file / buffer. This will generate the 
      // DOM tree. Catch all exceptions that the parser could not pass to 
      // our errorhandler.

      // Distinguish if we have to parse a file or the memory buffer
      if( this->file_ != "" )
        parser_->parse(file_.c_str()); // Parse File
      else
        parser_->parse(*buf_); // Parse Memory String
    }
    catch(const XMLException &event)
    {
        EXCEPTION("Error parsing '" << file_ << "' -> '" << event.getMessage() << "'");
    }
    catch(const DOMException &event )
    {
        const unsigned int maxChars = 2047;
        XMLCh errText[maxChars + 1];

        std::stringstream ss;
        ss << "DOM error in '" << file_ << "': DOMException.code = " << event.code;

        if(DOMImplementation::loadDOMExceptionMsg(event.code, errText, maxChars))
          ss << " Message is: " << Transcode(errText);

        EXCEPTION(ss.str());
    }

    // Obtain and validate root element of document tree
    DOMDocument* doc = parser_->getDocument();
    DOMNodeList* children = doc->getChildNodes();

    // some final checking, cannot imagine a problem here
    if(children->getLength() != 1)
        EXCEPTION("document root has " << children->getLength() << " children, expected 1");

    if(children->item(0)->getNodeType() != DOMNode::ELEMENT_NODE)
        EXCEPTION("root node type is " <<  children->item(0)->getNodeType());

    root_ = children->item(0);

  }

  Xerces::~Xerces()
  {
    if(parser_ != NULL && parser_->getErrorHandler() != NULL)
      delete parser_->getErrorHandler();

    if(parser_ != NULL)
    {
      delete parser_;
      parser_ = NULL;
    }

    if(buf_ != NULL)
    {
      delete parser_;
      buf_ = NULL;
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
    // std::cout << "node " << Transcode(node->getNodeName()) << " has type " << node->getNodeType();
    // std::cout << " value = " << (node->getNodeValue() != NULL ? Transcode(node->getNodeValue()) : "null") << std::endl;

    switch(node->getNodeType())
    {

    case DOMNode::TEXT_NODE:
    {
      // if we are a text node, we "are" the value of our parent.
      std::string temp = Transcode(node->getNodeValue());
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
      PtrParamNode newNode = boost::make_shared<ParamNode>(ParamNode::EX, type);
      parent->AddChildNode( newNode );
      pn = newNode;
      // we work with the this just added element - here we avoid any
      // potential copy constructor issues
      // pn = parent->GetChildren().Last();
    }
    else
    {
      // no new child created but we modify the parent directly
      pn = parent;
    }

    std::string nodename = Transcode(node->getNodeName());
    pn->SetName(nodename);

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

  std::string Xerces::Transcode(const XMLCh* input) 
  {
    std::string ret;
    bool transcoded = false;
    char *value = NULL;

    // Since the XMLString::transcode() method allways throws an exception
    // if problems occur, we enclose it in a try-catch block, in order to be
    // able to properly react to an error condition. If we cannot transcode
    // to the local codepage, we fall back to ASCII afterwards
    try 
    {
      
      value = XMLString::transcode(input);
      ret = value;
      XMLString::release(&value);
      transcoded = true;
      
      // The TranscodeToStr also throws exceptions all the time.
      // This cannot be prevented. 
#if 0
      TranscodeToStr tr(node->getNodeValue(), "ASCII");
      std::string ret((char*)tr.str(), tr.length());
#endif
      
      // We can also manually generate a transcoder for the local code page.
      // It also throws exceptions, which cannot be prevented. 
#if 0
      XMLLCPTranscoder* lcpt = XMLPlatformUtils::fgTransService->makeNewLCPTranscoder(XMLPlatformUtils::fgMemoryManager);
      
      XMLSize_t charsEaten = lcpt->calcRequiredSize(node->getNodeValue() );
      char* resultXMLString_Encoded;
      resultXMLString_Encoded = lcpt->transcode( node->getNodeValue() );
      
      resultXMLString_Encoded[charsEaten] = 0;
      ret = resultXMLString_Encoded;
      
      // release the memory
      delete lcpt;
#endif        
      
    } catch (XMLException& ex) 
    {
      XMLString::release(&value);

      value = XMLString::transcode(ex.getType());
      std::string exType = value;
      XMLString::release(&value);

      value = XMLString::transcode(ex.getMessage());
      std::string exMsg = value;
      XMLString::release(&value);
      
      WARN("Xerces " << exType << ": " << exMsg
           << std::endl << std::endl
           << "These problems usually occur if characters are present in the XML"
           << std::endl
           << "file which can not be represented in the local code page. Please"
           << std::endl
           << "check your XML file and locale settings. Under Linux these are"
           << std::endl
           << "determined by the LC_CTYPE and LANG environment variables."
           << std::endl
           << std::endl
           << "Falling back to ASCII transcoder and replacing invalid characters." );
    }

    // If the input string could not yet be transcoded, we fall back to an ASCII
    // transcoder and just replace all invalid characters.
    if(!transcoded)
    {
      // construct a transcoder to ASCII
      XMLTransService::Codes resCode;
      XMLTranscoder* inpEncTranscoder = NULL;
      inpEncTranscoder =XMLPlatformUtils::fgTransService->
                        makeNewTranscoderFor("ASCII", resCode, 16*1024, 
                                             XMLPlatformUtils::fgMemoryManager);
      
      // transcode the string into ASCII
      XMLSize_t charsEaten;
      char resultXMLString_Encoded[16*1024+4];
      inpEncTranscoder->transcodeTo(input,
                                    XMLString::stringLen(input),
                                    (XMLByte*) resultXMLString_Encoded,
                                    16*1024,
                                    charsEaten,
                                    XMLTranscoder::UnRep_RepChar );
      
      // release the memory
      delete inpEncTranscoder;
      
      resultXMLString_Encoded[charsEaten] = 0;
      ret = resultXMLString_Encoded;

      transcoded = true;
    }

    return ret;
  }
  
  Xerces::EventHandler::EventHandler(const Xerces* xerces)
  {
     this->xerces_ = xerces;
  }


  void Xerces::EventHandler::warning(const SAXParseException &event )
  {
    std::stringstream os;
    os << "WARN parsing the xml file '" << xerces_->file_ << "' in line "
       << event.getLineNumber() << ", column " << event.getColumnNumber() << std::endl
       << "-> '" << Xerces::Transcode(event.getMessage()) << "'" << std::endl
       << " schema: '" << (!xerces_->schema_.empty() ? xerces_->schema_ : "<no-schema>") << "'";
    // killme use the new log stuff from Andi
    std::cerr << os.str() << std::endl;
  }


  void Xerces::EventHandler::error(const SAXParseException &event)
  {
    std::stringstream ss;
    ss << "XML parsing error: '" << Xerces::Transcode(event.getMessage()) << "'"
       << std::endl << "file: '" << (xerces_->file_ == "" ? "no-file" : xerces_->file_) << "'"
       << " line: " << event.getLineNumber() << " column: " << event.getColumnNumber()
       << std::endl << "schema: '" << (!xerces_->schema_.empty() ? xerces_->schema_ : "no-schema") << "'";
    throw Exception(ss.str()); // the macro EXCEPTION adds this file number which is too much information for the user in this case

  }


  void Xerces::EventHandler::fatalError(const SAXParseException &event)
  {
    std::cout << "Xerces::EventHandler::fatalError\n";
    std::cout.flush();
    error(event);
  }

  void Xerces::EventHandler::resetErrors()
  {
  }



} // end of namespace

