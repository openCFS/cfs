// includes for xerces
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationLS.hpp>
#include <xercesc/dom/DOMWriter.hpp>
#include <xercesc/framework/StdOutFormatTarget.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/util/XMLUni.hpp>

// we want to use the Xerces-C++ namespace
using namespace xercesc;

// we want to use new with nothrow
#include <new>

#include <string>

#include "ErrorHandler.hh"
#include "SPHandler.hh"



namespace CoupledField {


  // **************************************************************************
  //   Constructors and Destructors
  // **************************************************************************


  // ===================================
  //   Constructor for SPHandler class
  // ===================================
  SPHandler::SPHandler( char *fname ): parser_(NULL) {

    // Initialise the XML4C2 system
    try {
      XMLPlatformUtils::Initialize();
    }
    catch( const XMLException &event ) {
      std::cerr << "\n The following error occured during initialisation of "
		<< "xerces-c!\n" << X2C(event.getMessage()) << std::endl;
      MYABORT;
    }

    // Create the parser, force validation using full schema checking with
    // namespaces and tell it to drop irrelevant white space
    parser_ = new XercesDOMParser();
    parser_->setValidationScheme(XercesDOMParser::Val_Always);
    parser_->setDoNamespaces(true);
    parser_->setDoSchema(true);
    parser_->setValidationSchemaFullChecking(true);
    parser_->setIncludeIgnorableWhitespace(false);

    // We may separate the schema file from the instance file
    std::string mySchema = "http://www.cfs++.org ";
    mySchema += XMLSCHEMA;
    mySchema += "/CFS.xsd";
    parser_->setExternalSchemaLocation( mySchema.c_str() );

    // Have not yet understood what an entity reference node is, but it seems
    // we do not need them
    parser_->setCreateEntityReferenceNodes(false);

    // Attach our own error handler to the parser
    ErrorHandler *errHandler = new ErrorHandler();
    parser_->setErrorHandler(errHandler);

    // Parse and validate the XML file. This will generate the DOM tree. Catch
    // all exceptions that the parser could not pass to our error handler.
    try {
      parser_->parse( fname );
    }
    catch( const XMLException &event ) {
      std::cerr << "\n The following XML error was encountered during "
		<< "parsing:\n"
		<< X2C( event.getMessage() ) << std::endl;
      MYABORT;
    }
    catch( const DOMException &event ) {
      const unsigned int maxChars = 2047;
      XMLCh errText[maxChars + 1];

      std::cerr << "\nDOM Error during parsing! DOMException code is:\n"
		<< event.code << std::endl;

      if( DOMImplementation::loadDOMExceptionMsg(event.code, errText,
						 maxChars) ) {
	std::cerr << "Message is: " << X2C(errText) << std::endl;
      }

      MYABORT;
    }
    catch(...) {
      std::cerr << "\n An unknown error occurred during parsing!\n "
		<< " All I can say is that it was neither an XMLException "
		<< " nor a DOMException." << std::endl;
      MYABORT;
    }

    // Obtain and validate root element of document tree
    DOMDocument *doc = parser_->getDocument();
    DOMNodeList *children = doc->getChildNodes();
    if ( children->getLength() != 1 ) {
      std::cerr << " Document root has more than one child!\n"
		<< " We are in real trouble here!\n"
		<< " Exiting ...\n\n";
      MYABORT;
    }
    if ( children->item(0)->getNodeType() != DOMNode::ELEMENT_NODE ) {
      std::cerr << " Root node is not of type DOMNode::ELEMENT_NODE!\n"
		<< " We are in real trouble here!\n"
		<< " Exiting ...\n\n";
      MYABORT;
    }
    rootelem_ = (DOMElement *)(children->item(0)); 
  }


  // ==================================
  //   Destructor for SPHandler class
  // ==================================
  SPHandler::~SPHandler() {
    if ( parser_ != NULL ) {
      delete parser_;
    }
  }


  // ===========================================
  //   Default constructor for SPHandler class
  // ===========================================
  SPHandler::SPHandler() {}

}
