#ifdef USE_XERCES

#include <iostream>

#include "XMLParserErrorHandler.hh"

#include "xercesc/util/XercesDefs.hpp"
#include "xercesc/sax/ErrorHandler.hpp"
#include "xercesc/sax/SAXParseException.hpp"

// what to do in case of an emergency
// #define XMLPARSERABORT ERROR( "", __FILE__, __LINE__ );
#define XMLPARSERABORT exit(1)

namespace CoupledField
{
 
  // ========================
  //   Handler for warnings
  // ========================
  void XMLParserErrorHandler::warning(const xercesc::SAXParseException &event){
    std::cerr << std::endl
	      << "\n Sorry to say, but the parser seems to be displeased "
	      << "with the input file.\n"
	      << " Please check line " << event.getLineNumber()
	      << ", column " << event.getColumnNumber() << '\n'
	      << " It issued the following warning:\n\n "
	      << xercesc::XMLString::transcode( event.getMessage() )
	      << "\n\n"
              << " XML parsers used the Schema: http://www.cfs++.org\n "
              << XMLSCHEMA <<  "/CFS.xsd\n" << std::endl;
    XMLPARSERABORT;
  }


  // =============================
  //   Handler for simple errors
  // =============================
  void XMLParserErrorHandler::error(const xercesc::SAXParseException &event){
    std::cerr << "\n Oh, we are in trouble here ...\n"
	      << " The parser seems to have encountered an error"
	      << " in the input file.\n"
	      << " Please check line " << event.getLineNumber()
	      << ", column " << event.getColumnNumber()
	      << "\n The only hint I can give is the following:\n\n "
	      << xercesc::XMLString::transcode( event.getMessage() )
	      << "\n\n"
              << " XML parsers used the Schema: http://www.cfs++.org\n "
              << XMLSCHEMA <<  "/CFS.xsd\n" << std::endl;

    XMLPARSERABORT;
  }


  // ============================
  //   Handler for fatal errors
  // ============================
  void XMLParserErrorHandler::fatalError(const xercesc::SAXParseException
					 &event) {
    std::cerr << "\n Oh, oh! We are in major trouble here ...\n"
	      << " The parser seems to have encountered a fatal error in the"
	      << " input file!\n"
              << " Please check line " << event.getLineNumber()
	      << ", column " << event.getColumnNumber()
	      << "\n This is what it was complaining about:\n\n "
	      << xercesc::XMLString::transcode( event.getMessage() )
	      << "\n\n"
              << " XML parsers used the Schema: http://www.cfs++.org\n "
              << XMLSCHEMA <<  "/CFS.xsd\n" << std::endl;
    XMLPARSERABORT;
  }


  // =================================
  //  Only provided for completeness
  // =================================
  void XMLParserErrorHandler::resetErrors(){};

}

#endif
