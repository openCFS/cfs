#ifndef FILE_SPHANDLER
#define FILE_SPHANDLER

#include <iostream>
#include <vector>

#include "xercesc/util/XercesDefs.hpp"
#include "xercesc/sax/ErrorHandler.hpp"
#include "xercesc/util/XMLString.hpp"
#include "xercesc/parsers/XercesDOMParser.hpp"
#include "xercesc/dom/DOMImplementation.hpp"
#include "xercesc/dom/DOMImplementationLS.hpp"
#include "xercesc/dom/DOMWriter.hpp"

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

    //! Report parameters to the specified output stream

    //! This method can be used to serialise the DOM tree the parser has
    //! generated and write it to a specified output stream. Thus, we are
    //! able to see the whole tree including all optional parameters whose
    //! values were set by the Schema specification and not in the input
    //! XML file.
    //! \note
    //! - The current implementation of this method makes use of several
    //!   experimental featuers of xerces that might change in future releases!
    //! - Currently the stream parameter is ignored and output is sent to
    //!   standard output.
    void PrintTree();

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

    // ************************************************************************
    //   Private Auxilliary Methods: Diverse
    // ************************************************************************

    //! Turn on/off a feature of the DOMWriter

    //! This method can be used to turn a feature of the DOMWriter / serialiser
    //! on or off. The specifications demand that an application should query
    //! a serialiser first, before attempting to set a certain feature, or be
    //! prepared to catch the resulting exception, if it cannot. Thus, this
    //! method will first query the serisalier, whether it does support
    //! setting the feature to the specified value. If the serialiser does,
    //! then the method will set the feature and return true. If the serialsier
    //! does not supported it, then the behaviour depends on the value of the
    //! shouldHave argument. Setting shouldHave to true, indicates that the
    //! application relies on the serialiser having the desired feature
    //! capability. If the latter is not the case, an error will be issued. If
    //! shouldHave is false, then only a warning will be issued, if the
    //! serialiser does not support this feature capability and false will be
    //! returned.
    //! \note For a list of features that every serialiser should support see
    //!       the Xerces documentation of the DOMWriter class or the Document
    //!       object model specification.
    //! \param serialiser  pointer to the DOMWriter / serialiser
    //! \param feature     feature specification encoded as constant XMLCh
    //!                    array; feature names are static attributes of the
    //!                    xercesc::XMLUni class
    //! \param featureVal  value to which the feature should be set
    //! \param shouldHave  boolean signalling that the application relies on
    //!                    the serialiser having the specified feature
    //!                    capability
    //! \return a boolean signaling whether the serialiser supports setting the
    //!         specified feature to the specified value.
    bool DOMWriterSetFeature( xercesc::DOMWriter *serialiser,
			      const XMLCh *const feature,
			      bool featureVal,
			      bool shouldHave );
  };

}

#endif
