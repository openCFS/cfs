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

// #include <xercesc/util/NameIdPool.hpp>
// #include <xercesc/framework/XMLValidator.hpp>

#include <xercesc/validators/schema/SchemaValidator.hpp>
// #include <xercesc/validators/common/ContentSpecNode.hpp>
// #include <xercesc/validators/schema/SchemaSymbols.hpp>
// #include <xercesc/util/OutOfMemoryException.hpp>

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
    parser_->setCreateSchemaInfo(true);
    parser_->setCreateCommentNodes(false);

    // We may separate the schema file from the instance file
    std::string mySchema = "http://www.cfs++.org ";
    mySchema += XMLSCHEMA;
    mySchema += "CFS.xsd";
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

    // Experimental
    // tryStuff();
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


  // =====================================
  //   Print DOM tree to standard output
  // =====================================
  void SPHandler::PrintTree() {

    try {

      // Obtain a DOMImplementationLS implementation from the registry and
      // use this factory to generate an instance of DOMWriter
      DOMImplementation *impl =
	DOMImplementationRegistry::getDOMImplementation( C2X( "LS" ) );
      DOMWriter *theSerializer =
	((DOMImplementationLS*)impl)->createDOMWriter();

      // Plug in a format target to receive the resultant
      // XML stream from the serializer.
      //
      // StdOutFormatTarget prints the resultant XML stream
      // to stdout once it receives any thing from the serializer.
      XMLFormatTarget *myFormTarget = new StdOutFormatTarget();

      // Query parser for the DOM representation
      DOMNode *doc = parser_->getDocument();

      // Set some features
      DOMWriterSetFeature( theSerializer, XMLUni::fgDOMWRTFormatPrettyPrint,
			   true, true );
      DOMWriterSetFeature( theSerializer,
			   XMLUni::fgDOMWRTDiscardDefaultContent, false, true);

      // Report what comes next
      std::cout	<< "\n"
		<< " --------------------------------------------------------"
		<< "----------------------\n"
		<< "   The following is a serialisation of the DOM tree "
		<< "generated by the parser.\n   It includes all values "
		<< "for optional attributes that were not specified in\n"
		<< "   the XML file and were added by the parser based on the"
		<< " Schema definitions\n   used for validation.\n\n"
		<< "   Note:\n\n"
		<< "   The W3C Schema specifications demand that only "
		<< " default values for un-\n   supplied attributes"
		<< " are added by the parser. Optional elements are not\n"
		<< "   added to the tree and default values are only"
		<< " supplied for empty but\n   not for missing elements.\n"
		<< " --------------------------------------------------------"
		<< "----------------------\n"
		<< std::endl;

      // Perform the serialization through DOMWriter::writeNode();
      theSerializer->writeNode( myFormTarget, *doc );

      // Clean up
      delete theSerializer;

    }
    catch (XMLException& e) {
      std::cerr << "Some weird shit happened!\n" << std::endl;
    }
  }


  // ================================
  //   Set feature of the DOMWriter
  // ================================
  bool SPHandler::DOMWriterSetFeature( DOMWriter *serialiser,
				       const XMLCh *const feature,
				       bool featureVal,
				       bool shouldHave ) {

    if ( featureVal == true ) {
      std::cout << " Enabling feature '";
    }
    else {
      std::cout << " Disabling feature '";
    }
    std::cout << X2C( feature ) << "'." << std::endl;

    // Check if serialiser supports setting the feature
    // to the specified value
    bool supportsFeature = serialiser->canSetFeature( feature, featureVal );
    if ( supportsFeature == false ) {
      if ( shouldHave == true ) {
	std::cout << "\n ERROR: The serialiser does not support enabling "
		  << "feature '" << X2C(feature) << "'\n\n";
	exit(-1);
      }
      else {
	std::cout << "\n WARNING: The serialiser does not support enabling "
		  << "feature '" << X2C(feature) << "'\n\n";
      }
    }

    // The serialiser supports it, so enable/disable the feature
    else {
      serialiser->setFeature( feature, featureVal );
    }

    // Return, whether the serialiser supported setting
    // the feature
    return supportsFeature;
  }

#include "tryStuff.cc"

}
