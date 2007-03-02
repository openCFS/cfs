#include <def_use_xerces.hh>

#ifdef USE_XERCES

#ifndef FILE_CFSParserErrorHandler
#define FILE_CFSParserErrorHandler

#include <xercesc/util/XercesDefs.hpp>
#include <xercesc/sax/ErrorHandler.hpp>
#include <xercesc/sax/SAXParseException.hpp>

namespace CoupledField
{

  // Forward declaration of classes
  class xercesc::SAXParseException;

  //! Class for handling errors encountered by the parser

  //! This class provides methods for handling errors encountered by the
  //! XML parser of Xerces. It is an implementation of the ErrorHandler class
  //! of Xerces. If an object of this class is registered with the DOM parser
  //! the latter will use the former as an interface instead of throwing an
  //! exception.
  //!
  //! Currently all the class provides is to output a more or less nicely
  //! formatted error message and terminate program execution. This is done for
  //! all three severities of errors: warning, error, fatal error.
  class XMLParserErrorHandler: public xercesc::ErrorHandler
  {

  public:

    //! Handler for warnings
    void warning( const xercesc::SAXParseException &event );

    //! Handler for simple errors
    void error( const xercesc::SAXParseException &event );

    //! Handler for fatal errors
    void fatalError( const xercesc::SAXParseException &event );

    //! Only provided for completeness
    void resetErrors();

  };

}

#endif

#endif
