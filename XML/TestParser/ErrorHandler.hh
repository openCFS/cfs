#ifndef FILE_ErrorHandler
#define FILE_ErrorHandler


#include <xercesc/util/XercesDefs.hpp>
#include <xercesc/sax/ErrorHandler.hpp>
#include <xercesc/sax/SAXParseException.hpp>
#include <iostream>


// what to do in case of an emergency
#define MYABORT exit(1)


// indicate using Xerces-C++ namespace in general
// XERCES_CPP_NAMESPACE_USE
// using namespace xercesc;

namespace CoupledField {

  //! Class for handling errors encountered by the parser

  //! This class provides methods for handling errors encountered by the
  //! parser. It is an implementation of the ErrorHandler class of xerces.
  //! If an object of this class is registered with the DOM parser the latter
  //! will use the former as an interface instead of throwing an exception.\n\n
  //! Currently all the class provides is outputting a more or less nicely
  //! formatted error message and terminate program execution. This is done for
  //! all three severities of errors: warning, error, fatal error.
  class ErrorHandler: public xercesc::ErrorHandler {

  public:

    //! Constructor
    ErrorHandler(){};

    //! Destructor
    ~ErrorHandler(){};

    //! Handler for warnings
    void warning( const xercesc::SAXParseException &event ) {
      std::cerr << std::endl
		<< " Sorry to say, but the parser seems to be displeased "
		<< "with the input file." << std::endl
		<< " Please check line " << event.getLineNumber()
		<< ", column " << event.getColumnNumber() << std::endl
		<< " It issued the following warning:"
		<< std::endl << std::endl << " "
		<< xercesc::XMLString::transcode( event.getMessage() )
		<< std::endl << std::endl;
      MYABORT;
    }

    //! Handler for simple errors
    void error( const xercesc::SAXParseException &event ) {
      std::cerr << std::endl << " Oh, we are in trouble here ..." << std::endl
		<< " The parser seems to have encountered an error"
		<< " in the input file." << std::endl
		<< " Please check line " << event.getLineNumber()
		<< ", column " << event.getColumnNumber()	<< std::endl
		<< " The only hint I can give is the following:"
		<< std::endl << std::endl << " "
		<< xercesc::XMLString::transcode( event.getMessage() )
		<< std::endl << std::endl;
      MYABORT;
    }

    //! Handler for fatal errors
    void fatalError(const xercesc::SAXParseException &event ) {
      std::cerr << std::endl << " Oh, oh! We are in major trouble here ..."
		<< std::endl
		<< " The parser seems to have encountered a fatal error in the"
		<< " input file!" << std::endl
		<< " Please check line " << event.getLineNumber()
		<< ", column " << event.getColumnNumber() << std::endl
		<< " This is what it was complaining about:"
		<< std::endl << std::endl << " "
		<< xercesc::XMLString::transcode( event.getMessage() )
		<< std::endl << std::endl;
      MYABORT;
    }

    //! Only provided for completeness
    void resetErrors(){};

  };

}

#endif
