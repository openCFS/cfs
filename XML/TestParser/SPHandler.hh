#ifndef FILE_SPHANDLER
#define FILE_SPHANDLER

#include <iostream>
#include <vector>

#include "xercesc/util/XercesDefs.hpp"
#include "xercesc/sax/ErrorHandler.hpp"
#include "xercesc/util/XMLString.hpp"
#include "xercesc/parsers/XercesDOMParser.hpp"

namespace CoupledField {

  //! Class for validated parsing of an XML file

  //! This class is intended as test parser for an XML parameter file
  //! containing parameters for a CFS++ simulaton. The parameters are read from
  //! an XML file whose grammar is specified by a given schema and stored in
  //! an abstract tree structure according to the DOM model.
  //!
  //! Apache's xerces library is employed for parsing the input file,
  //! generating the DOM tree and accessing the information contained in its
  //! nodes. The parser will validate the input file against the rules given
  //! by the schema description.
  class SPHandler {

  public:

    //! Constructor

    //! The constructor expects as input the name of the input file containing
    //! the steering parameters. The file has to be in XML format. The
    //! constructor will parse and validate the file and generate an internal
    //! DOM tree containing the parameters.
    SPHandler( char *fname );

    //! Default destructor

    //! The destructor is deep in the sense that it will delete the parser
    //! and the DOM tree it owns.
    ~SPHandler();

  private:

    //! Default constructor

    //! The default constructor is private in order to disallow its use. We
    //! require a file name in the constructor for passing to the XM Lparser.
    SPHandler();

    //! Parser object

    //! This is a pointer to the XercesDOMParser object. The parser is the
    //! owner of the DOMDocument it generates. Thus, deleting the parser will
    //! also delete the document, i.e. all information on the steering
    //! parameters.
    xercesc::XercesDOMParser *parser_;

    //! Root element of document tree

    //! This is a pointer to the root element of the DOM document tree. It is
    //! not the be mistaken for the root node of the document. The former is
    //! a child of the latter.
    xercesc::DOMElement *rootelem_;
    // ************************************************************************
    //   Private Auxilliary Methods: Conversion Routines for Strings
    // ************************************************************************

    //@{
    //! \name Conversion routines for string representations:

    //! Transcodes a string to native code-page

    //! The C++ DOM of xerces uses the plain, null-terminated (XMLCh *) utf-16
    //! strings as the String type. This method converts such a string into a
    //! standard char array. The returned buffer is dynamically allocated and
    //! it is the responsibility of the caller to delete it when not longer
    //! needed.
    char* X2C( const XMLCh *const toTranscode ) {
      return xercesc::XMLString::transcode( toTranscode );
    };

    //! Transcodes a string from native code-page

    //! The C++ DOM of xerces uses the plain, null-terminated (XMLCh *) utf-16
    //! strings as the String type. This method converts a standard char array
    //! into such a string. The returned buffer is dynamically allocated and
    //! it is the responsibility of the caller to delete it when not longer
    //! needed.
    XMLCh* C2X( const char *const toTranscode ) {
      return xercesc::XMLString::transcode( toTranscode );
    };

    //! Transcodes an STL string to native code-page

    //! The C++ DOM of xerces uses the plain, null-terminated (XMLCh *) utf-16
    //! strings as the String type. This method converts such a string into a
    //! standard STL string.
    std::string& X2S( const XMLCh *const toTranscode ) {
      std::string *stlstring = new std::string( X2C( toTranscode ) );
      return *stlstring;
    };

    //! Transcodes an STL string from native code-page

    //! The C++ DOM of xerces uses the plain, null-terminated (XMLCh *) utf-16
    //! strings as the String type. This method converts an STL string into
    //! such a string. The returned buffer is dynamically allocated and
    //! it is the responsibility of the caller to delete it when not longer
    //! needed.
    const XMLCh* S2X( const std::string toTranscode ) {
      return xercesc::XMLString::transcode( toTranscode.c_str() );
    };

  };

}

#endif
