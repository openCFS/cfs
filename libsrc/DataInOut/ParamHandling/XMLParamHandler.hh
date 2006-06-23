#ifdef USE_XERCES

#ifndef FILE_XMLPARAMHANDLER
#define FILE_XMLPARAMHANDLER

// includes for Xerces
#include "xercesc/util/PlatformUtils.hpp"
#include "xercesc/dom/DOM.hpp"
#include "xercesc/dom/DOMImplementation.hpp"
#include "xercesc/dom/DOMImplementationLS.hpp"
#include "xercesc/dom/DOMWriter.hpp"
#include "xercesc/framework/StdOutFormatTarget.hpp"
#include "xercesc/framework/LocalFileFormatTarget.hpp"
#include "xercesc/parsers/XercesDOMParser.hpp"
#include "xercesc/util/XMLUni.hpp"


namespace CoupledField {

  // Forward declaration of classes
  class BaseParamHandler;

  //! Class for handling steering parameters for CFS++ simulations

  //! This class is intended for handling the input parameters steering a
  //! CFS++ simulation. The parameters are read from an XML file whose
  //! grammar is specified by a given schema and stored in an abstract tree
  //! structure according to the DOM model.
  //!
  //! Apache's xerces library is employed for parsing the input file,
  //! generating the DOM tree and accessing the information contained in its
  //! nodes. The parser will validate the input file against the rules given
  //! by the schema description.
  //!
  //! The class specifies several methods for accessing the values of steering
  //! parameters. These hide the underlying XML/DOM structure.
  //! \todo Check the class for memory leaks. We might introduce some leaks by
  //! the string conversion routines. Also it is not yet clear, whether we are
  //! expected to free the results from getElementsByTagName() and the like,
  //! or if this is even forbidden.
  //! \note With respect to function tracing the policy of this class is the
  //! following. All public methods use ENTER_FCN, while all other methods,
  //! besides those defined in the header itself, use ENTER_IFCN. The latter
  //! do not appear in the trace log at all.
  class XMLParamHandler : public BaseParamHandler {

  public:

    //! Constructor

    //! The constructor expects as input the name of the input file containing
    //! the steering parameters. The file has to be in XML format. The
    //! constructor will parse and validate the file and generate an internal
    //! DOM tree containing the parameters.
    XMLParamHandler( const char * fname, 
                     const char * schema = NULL,
                     const char * defaultFile = NULL );

    //! Default destructor

    //! The destructor is deep in the sense that it will delete the parser
    //! and the DOM tree it owns.
    ~XMLParamHandler();


    // =======================================================================
    // Simple query methods without side constraints
    // =======================================================================

    //@{
    //! \name Simple query methods without side constraints

    void Get( const std::string key, std::string &value,
              const std::string section="",
              const std::string subsection="" );

    void Get( const std::string key, Integer &value,
              const std::string section="", const std::string subsection="" );

    void Get( const std::string key, UInt &value,
              const std::string section="", const std::string subsection="" );

    void Get( const std::string key, Double &value,
              const std::string section="", const std::string subsection="" );

    void GetList( const std::string key,
                  StdVector<std::string> &list,
                  const std::string section = "",
                  const std::string subsection = "" );

    void GetList( const std::string key, StdVector<Double> &list,
                  const std::string section = "",
                  const std::string subsection = "" );

    void GetList( const std::string key, StdVector<Integer> &list,
                  const std::string section = "",
                  const std::string subsection = "" );

    void GetList( const std::string key, StdVector<UInt> &list,
                  const std::string section = "",
                  const std::string subsection = "" );

    void Get( const StdVector<std::string> &keyVec,
              std::string &value );

    void Get( const StdVector<std::string> &keyVec,
              Double &value );

    void Get( const StdVector<std::string> &keyVec,
              Integer &value );

    void Get( const StdVector<std::string> &keyVec,
              UInt &value );
    //@}


    // =======================================================================
    // Query methods allowing side-constraints and expecting vectors as input
    // =======================================================================

    //@{
    //! \name Query methods allowing side-constraints
    //! These query methods expect the search parameters to be specified in
    //! vector format. The search parameters have the following meaning:\n\n
    //! <center>
    //! <table border="1" width="80%" cellpadding="10">
    //! <tr>
    //! <td>keyVec</td>
    //! <td>a vector containing keywords, i.e. names of elements. These names
    //! specify the 'path' leading from the root of the parameter tree to the
    //! desired element. The last entry of this vector, however, can also be
    //! the name of an attribute.</td>
    //! </tr>
    //! <tr>
    //! <td>attrVec</td>
    //! <td>This vector contains names of attributes, which will be used to
    //! restrict the search space.</td>
    //! </tr>
    //! <tr>
    //! <td>valVec</td>
    //! <td>This vector contains the values of the attributes used for
    //! restricting the search.</td>
    //! </tr>
    //! </table>
    //! </center>
    //! \n\n
    //! When trying to locate the desired parameter value(s) in the tree, an
    //! element is only considered to be a match, if its name fits to the
    //! name given in keyVec and it has an attribute with the name given in
    //! the attrVec vector, whose value matches the one specified in the
    //! valVec vector. For this test the k-th entries in the vectors belong
    //! together.\n\n
    //! If a match is found for the specified set of search parameters, the
    //! query method will return the value of the element or attribute
    //! specified by the last entry in the keyVec vector. The Get methods
    //! expect this final match to be unique, while the GetList methods will
    //! also return a (possibly empty) StdVector of results, if there is no
    //! or multiple matches.\n\n
    //! The following is an example set of search parameters. We are trying
    //! to determine here the analysis type for the mechanic PDE in second
    //! step of a multi-sequence analysis:
    //! \f[\mbox{keyVec} = \left(\begin{array}{c}
    //! \mbox{"cfsSimulation"} \\
    //! \mbox{"multiSequence"} \\
    //! \mbox{"step"} \\
    //! \mbox{"pde"} \\
    //! \mbox{"analysis"}
    //! \end{array}\right)
    //! \mbox{attrVec} = \left(\begin{array}{c}
    //! \mbox{""} \\
    //! \mbox{""} \\
    //! \mbox{"index"} \\
    //! \mbox{"type"} \\
    //! \end{array}\right)
    //! \mbox{valVec} = \left(\begin{array}{c}
    //! \mbox{""} \\
    //! \mbox{""} \\
    //! \mbox{"2"} \\
    //! \mbox{"mechanic"} \\
    //! \end{array}\right)
    //! \f]
    //! Note that the first two entries of the attrVec vector are empty
    //! strings. This indicates that the elements matching the first two
    //! entries of the keyVec vector will not be tested for having attributes
    //! of a certain name and value.

    //! Perform search and return values of matches as strings.
    void GetList( const StdVector<std::string> &keyVec,
                  const StdVector<std::string> &attrVec,
                  const StdVector<std::string> &valVec,
                  StdVector<std::string> &list );

    //! Perform search and return values of matches as Doubles.
    void GetList( const StdVector<std::string> &keyVec,
                  const StdVector<std::string> &attrVec,
                  const StdVector<std::string> &valVec,
                  StdVector<Double> &list );

    //! Perform search and return values of matches as Integers.
    void GetList( const StdVector<std::string> &keyVec,
                  const StdVector<std::string> &attrVec,
                  const StdVector<std::string> &valVec,
                  StdVector<Integer> &list );

    //! Perform search and return values of matches as unsigned Integers.
    void GetList( const StdVector<std::string> &keyVec,
                  const StdVector<std::string> &attrVec,
                  const StdVector<std::string> &valVec,
                  StdVector<UInt> &list );

    //! Search for a unique match and return value as string.
    void Get( const StdVector<std::string> &keyVec,
              const StdVector<std::string> &attrVec,
              const StdVector<std::string> &valVec,
              std::string &value );

    //! Search for a unique match and return value as Double.
    void Get( const StdVector<std::string> &keyVec,
              const StdVector<std::string> &attrVec,
              const StdVector<std::string> &valVec,
              Double &value );

    //! Search for a unique match and return value as Integer.
    void Get( const StdVector<std::string> &keyVec,
              const StdVector<std::string> &attrVec,
              const StdVector<std::string> &valVec,
              Integer &value );

    //! Search for a unique match and return value as unsigned Integer.
    void Get( const StdVector<std::string> &keyVec,
              const StdVector<std::string> &attrVec,
              const StdVector<std::string> &valVec,
              UInt &value );

    //! Perform search and return values of matches as Doubles.
    void GetDim1xDim2Tensor( const StdVector<std::string> &keyVec,
                             const StdVector<std::string> &attrVec,
                             const StdVector<std::string> &valVec,
                             const unsigned int &dim1,
                             const unsigned int &dim2,
                             Matrix<Double> &matr );

    //@}


    // =======================================================================
    // Specialised query methods
    // =======================================================================

    //@{
    //! \name Specialised query methods
    void GetPDEList( StdVector<std::string> &list );

    void GetIterCoupledPDEList( StdVector<std::string> &list,
                                const std::string sequenceTag = "");
    
    void GetDirectCouplingList( StdVector<std::string> &list,
                                const std::string sequenceTag = "");   

    void GetCoilList( StdVector<std::string> &list,
                      const std::string pde = "" ,
                      const std::string sequenceTag = "");

    void GetCoilType( std::string &coilType, const std::string coilName,
                      const std::string pde ="" );

    bool IsSet( const std::string key,
                   const std::string section = "",
                   const std::string subsection = "" );


    bool ContainElem( StdVector<std::string> &keyVec,
                         StdVector<std::string> &attrVec,
                         StdVector<std::string> &valVec );

    bool HasValue( const std::string key,
                      const std::string value,
                      const std::string section = "",
                      const std::string subsection = "" );
    //@}

  private:

    //! Default Constructor

    //! The default constructor is private in order to avoid a meaningless
    //! instantiation of an object of this type, which is not associated
    //! with a parameter file and a DOM tree.
    XMLParamHandler();

    //! Parsing of XML file

    //! This method generates a parser and parses the XML input
    //! file with filename xmlFile. A pointer to the generated parser is
    //! assigned to the parser input parameter and a pointer to the root
    //! element of the DOM tree representing the contents of the XML file
    //! is returned. If no schema is specified, no validation is performed.
    xercesc::DOMElement* ParseFile( xercesc::XercesDOMParser **parser,
                                    const char *xmlFile,
                                    const char *schema = NULL );

    //! Print search parameters to an output stream

    //! This method will print the search parameters given in the three input
    //! vectors to the output stream specified by myStream. The method
    //! currently uses the classical C routine fprintf to generate formatted
    //! output.
    void PrintSearchParams( const StdVector<std::string> &key,
                            const StdVector<std::string> &attrib,
                            const StdVector<std::string> &value,
                            FILE *myStream );

    //! Auxilliary routine for generating vectors with search parameters

    //! This auxilliary routine can be used to generate for search parameters
    //! specified in the simple style by giving 'key', 'section' and
    //! 'subsection' arguments corresponding keyVec, attrVec and valVec
    //! vectors
    void GenerateSearchParams( const std::string key,
                               const std::string section,
                               const std::string subsection,
                               StdVector<std::string> &keyVec,
                               StdVector<std::string> &attrVec,
                               StdVector<std::string> &valVec );

    //! Report parameters to the specified output stream

    //! This method can be used to serialise the specified DOM tree and write
    //! it to standard output.
    //! \note
    //! - We will see the whole tree including all optional parameters whose
    //!   values were set by the Schema specification and not in the input
    //!   XML file.
    //! - The current implementation of this method makes use of several
    //!   experimental features of xerces that might change in future
    //!   releases!
    void PrintTree( xercesc::DOMDocument *doc );

    //! Turn on/off a feature of the DOMWriter

    //! This method can be used to turn a feature of the DOMWriter/serialiser
    //! on or off. The specifications demand that an application should query
    //! a serialiser first, before attempting to set a certain feature, or be
    //! prepared to catch the resulting exception, if it cannot. Thus, this
    //! method will first query the serialiser, whether it does support
    //! setting the feature to the specified value. If the serialiser does,
    //! then the method will set the feature and return true. If the
    //! serialiser does not supported it, then the behaviour depends on the
    //! value of the shouldHave argument. %Setting shouldHave to true,
    //! indicates that the application relies on the serialiser having the
    //! desired feature capability. If the latter is not the case, an error
    //! will be issued. If shouldHave is false, then only a warning will be
    //! issued, if the serialiser does not support this feature capability and
    //! false will be returned.
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
    //! \return a boolean signalling whether the serialiser supports setting
    //!         the specified feature to the specified value.
    bool DOMWriterSetFeature( xercesc::DOMWriter *serialiser,
			      const XMLCh *const feature,
			      bool featureVal,
			      bool shouldHave );

    //! Auxilliary routine for clearing input vector

    //! This routine is used to empty the given input vector if it does not
    //! have zero length. If this happends and beVerbose_ is true, a warning
    //! message will be logged in case this happened accidentially.
    template <class T>
    void ClearVector( StdVector<T> &myVec ) {
      if ( myVec.IsEmpty() != true ) {
        if ( beVerbose_ == true ) {
          std::string msg = "Warning input vector was not empty!\n";
          msg += "Contents have been erased!";
          Info->Warning( msg );
        }
        myVec.Clear();
      }
    }


    // ***********************************************************************
    //   Private Auxilliary Methods: Problem Handlers
    // ***********************************************************************

    //@{
    //! \name Routines for handling problems and default parameters

    //! Error handler in case we have no match at all

    //! This method is called by those query routines, e.g. the Get methods,
    //! which expect to obtain a unqiue match, if no match is found in the
    //! parameter tree. The method will try to obtain the value form the
    //! tree containing the default parameters by calling CheckForDefault.
    //! If this, too, fails, the method will call the NoMatchErrorReporter.
    void NoMatchHandler( const StdVector<std::string> &keyVec,
                         const StdVector<std::string> &attrVec,
                         const StdVector<std::string> &valVec,
                         std::string &defaultValue );

    //! Error handler in case we have no unique match

    //! This method is called by those query routines, e.g. the Get methods,
    //! which expect to obtain a unqiue match, if a match is found but it is
    //! not unique. The error handler will issue a corresponding error message
    //! and terminate programm execution.
    void MultipleMatchHandler( const StdVector<std::string> &keyVec,
                               const StdVector<std::string> &attrVec,
                               const StdVector<std::string> &valVec,
                               const unsigned int nmatches );

    //! Error reporter for the case that we have no match at all

    //! This method is called by the methods handling the no match case, if
    //! no default value was found either. The method assembles a
    //! corresponding error message and sends this to the error method of
    //! the WriteInfo class.
    void NoMatchErrorReporter( const StdVector<std::string> &keyVec,
                               const StdVector<std::string> &attrVec,
                               const StdVector<std::string> &valVec );

    //! Warning reporter for the case that we have no match at all

    //! This method is called by the methods handling the no match case, if
    //! no default value was found either. The method assembles a
    //! corresponding warning message and sends this to the
    //! the WriteInfo class.
    void NoMatchWarningReporter( const StdVector<std::string> &keyVec,
                               const StdVector<std::string> &attrVec,
                               const StdVector<std::string> &valVec );

    //! Determine default value for a parameter

    //! This method will try to determine, whether for a given set of search
    //! parameters a default value exists. If this is the case, then the
    //! method will return true and defaultValue will contain the default
    //! value of the parameter in string format. If there is no default
    //! parameter specified, then the method will return false and
    //! defaultValue will be set to the empty string.
    bool CheckForDefault(  const StdVector<std::string> &keyVec,
                              const StdVector<std::string> &attrVec,
                              const StdVector<std::string> &valVec,
                              std::string &defaultValue );

    //@}


    // ***********************************************************************
    //   Private Auxilliary Methods: Conversion Routines for Strings
    // ***********************************************************************

    //@{
    //! \name Conversion routines for string representations:

    //! Transcodes a string to native code-page

    //! The C++ DOM of xerces uses the plain, null-terminated (XMLCh *) utf-16
    //! strings as the String type. This method converts such a string into a
    //! standard Char array. The returned buffer is dynamically allocated and
    //! it is the responsibility of the caller to delete it using the FreeC()
    //! method when not longer needed.
    Char* X2C( const XMLCh *const toTranscode ) {
      return xercesc::XMLString::transcode( toTranscode );
    };

    //! Transcodes a string from native code-page

    //! The C++ DOM of xerces uses the plain, null-terminated (XMLCh *) utf-16
    //! strings as the String type. This method converts a standard Char array
    //! into such a string. The returned buffer is dynamically allocated and
    //! it is the responsibility of the caller to delete it using the FreeX()
    //! method when not longer needed.
    XMLCh* C2X( const Char *const toTranscode ) {
      return xercesc::XMLString::transcode( toTranscode );
    };

    //! Transcodes an STL string to native code-page

    //! The C++ DOM of xerces uses the plain, null-terminated (XMLCh *) utf-16
    //! strings as the String type. This method converts such a string into a
    //! standard STL string. The returned string is dynamically allocated and
    //! it is the responsibility of the caller to delete it using the FreeS()
    //! method when not longer needed.
    std::string* X2S( const XMLCh *const toTranscode ) {
      char *auxString = X2C( toTranscode );
      std::string *stlstring = new std::string( auxString );
      FreeC( &auxString );
      return stlstring;
    };

    //! Transcodes an STL string from native code-page

    //! The C++ DOM of xerces uses the plain, null-terminated (XMLCh *) utf-16
    //! strings as the String type. This method converts an STL string into
    //! such a string. The returned buffer is dynamically allocated and
    //! it is the responsibility of the caller to delete it using the FreeX()
    //! method when not longer needed.
    XMLCh* S2X( const std::string toTranscode ) {
      return xercesc::XMLString::transcode( toTranscode.c_str() );
    };

    //! Free a dynamically generated XML string

    //! This method can be used to delete the memory dynamically allocated
    //! for an XML string (XMLCh array) generated by the C2X() or S2X()
    //! method. Note that since the initial memory allocation was performed
    //! by a memory manager within Xerces we better not do a simple delete.
    //! On return the pointer is set to NULL.
    void FreeX( XMLCh **string ) {
      xercesc::XMLString::release ( string );
    }

    //! Free a dynamically generated string from native code-page

    //! This method can be used to delete the memory dynamically allocated
    //! for a string (char *array) generated by the X2C() method.
    //! Note that since the initial memory allocation was performed by a
    //! memory manager within Xerces we better not do a simple delete.
    //! On return the pointer is set to NULL.
    void FreeC( char **string ) {
      xercesc::XMLString::release ( string );
    }

    //! Free a dynamically generated STL string

    //! This method can be used to delete the memory dynamically allocated
    //! for an STL string generated by the X2S() method.
    void FreeS( std::string **string ) {
      delete (*string);
      (*string) = NULL;
    }

    //! Transcodes an STL string to Double

    //! The method converts an STL string to a Double value by the help of the
    //! atof function.
    Double String2Double( const std::string toTranscode ) {
      return (Double)atof( toTranscode.c_str() );
    };

    //! Transcodes an STL string to Integer

    //! The method converts an STL string to an Integer value by the help of
    //! the atoi function.
    Integer String2Integer( const std::string toTranscode ) {
      return (Integer)atoi( toTranscode.c_str() );
    };
    
    //! Transcodes an STL string to UInt
    
    //! The method converts an STL string to an unsigned Integer value by the 
    //! help of the atoi function.
    Integer String2UInt( const std::string toTranscode ) {
      UInt retVal = 0;
      int aux = 0;
      aux = atoi( toTranscode.c_str() );
      if ( aux >= 0 ) {
        retVal = (UInt)aux;
      }
      else {
        (*error) << "String2UInt: Cannot convert negative value " << aux
                 << " to unsigned integral type UInt";
        Error( __FILE__, __LINE__ );
      }
      return retVal;
    };

    //@}


    // ***********************************************************************
    //   Private Auxilliary Methods: Conversion Routines for Nodes
    // ***********************************************************************

    //@{
    //! \name Conversion routines for nodes:

    //! Convert a DOMNode to a DOMElement

    //! This method will convert a pointer to a DOMNode to a pointer to a
    //! DOMElement. Since DOMElement is derived from DOMNode this actually is
    //! a down-cast. An error will be issued, when the type of the node is
    //! not ELEMENT_NODE.
    xercesc::DOMElement* Node2Elem( xercesc::DOMNode *node );

    //! Convert a DOMNode to a DOMAttr

    //! This method will convert a pointer to a DOMNode to a pointer to a
    //! DOMAttr. Since DOMAttr is derived from DOMNode this actually is
    //! a down-cast. An error will be issued, when the type of the node is
    //! not ATTRIBUTE_NODE.
    xercesc::DOMAttr* Node2Attr( xercesc::DOMNode *node );

    //! Convert element to string

    //! This method takes a DOMNode and converts is value to a standard
    //! STL string. For this to work the type of the node must be TEXT_NODE,
    //! otherwise an error is issued.
    void Node2String( xercesc::DOMNode *textnode, std::string &value );

    //@}


    // ***********************************************************************
    //   Private Auxilliary Methods: Search Routines
    // ***********************************************************************

    //@{
    //! \name Search routines:

    // Search (restricted) tree for attributes matching keyword 

    //! Find all elements matching keys and side-constraints

    //! This is the central search routine of the class. It works in a
    //! recursive fashion. On every level of the recursion all elements in the
    //! current collection of sub-trees are determined, which match the
    //! keyword for this level and have a certain attribute with a prescribed
    //! value. These elements form the root nodes of the subtrees for the
    //! search on the next level. The routine returns a pointer to a vector
    //! containing the elements found on the lowest level of the recursion. If
    //! there are no matches on this or the previous levels, the return vector
    //! will have zero length.
    //! \param keys     Vector of keywords
    //! \param attribs  Vector of attributes used as side-constraints
    //! \param aValues  Vector of attributes' values used as side-constraints
    //! \param treetop  Root node of the subtree in which to search
    //! \param curdepth Current level in the recurive search
    StdVector<xercesc::DOMElement *>*
    FindMatchingElements( StdVector<std::string> &keys,
                          StdVector<std::string> &attribs,
                          StdVector<std::string> &aValues,
                          xercesc::DOMElement *treetop,
                          unsigned int curdepth );

    //! Find all attributes matching keys and side-constraints

    //! This method will determine the value of all attributes which have a
    //! prescribed name. The path to the elements that are tested for having
    //! the prescribed attribute is specified by given the certain keywords
    //! and the corresponding side-constraints. The method makes use of
    //! FindMatchingElements for determining these elements. The method
    //! returns a pointer to a vector of the found attributes.
    //! If no attributes were found, the vector has zero length.
    //! \param attr_key Name of desired attribute
    //! \param keys     Vector of keywords
    //! \param attribs  Vector of attributes used as side-constraints
    //! \param aValues  Vector of attributes' values used as side-constraints
    //! \param treetop  Root node of the subtree in which to search
    StdVector<xercesc::DOMAttr *>*
    FindMatchingAttributes( std::string attr_key,
                            StdVector<std::string> &keys,
                            StdVector<std::string> &attribs,
                            StdVector<std::string> &aValues,
                            xercesc::DOMElement *treetop );

    //! Find values of all elements and attributes matching search criteria

    //! This method will determine all elements and attributes matching the
    //! search criteria specified by the given keywords and side-constraints.
    //! It makes use of the FindMatchingElements and FindMatchingAttributes
    //! methods for this. If there are matches for both elements and
    //! attributes an error will be reported and program execution terminated.
    //! The method will extract the values of the found elements or attributes
    //! and return them in a vecotr of strings. If there are no matches, this
    //! vector will have zero length.
    //! \param key      Vector of keywords
    //! \param attrib   Vector of attributes used as side-constraints
    //! \param aValue   Vector of attributes' values used as side-constraints
    //! \param list     Vector containing the found values
    //! \param treeTop  Root node of the subtree in which to search
    void FindAllMatches( const StdVector<std::string> &key,
                         const StdVector<std::string> &attrib,
                         const StdVector<std::string> &aValue,
                         StdVector<std::string> &list,
                         xercesc::DOMElement *treeTop );


    // Central auxilliary search method for get functions

    // This is the central auxilliary search method used by the public get
    // methods. It will try to find the specified keyword in the parameter
    // tree and returning the values of the respective elements or attributes
    // in a vector of strings. The search can be restricted to certain
    // subtrees by specifying keywords for section and subsection. The method
    // will return an empty vector, if there is no match at all. It will issue
    // an error message, if there are matches for both, elements and
    // attributes. The final argument to the method is a pointer to a vector.
    // If this pointer is not the NULL pointer on input (which is the default)
    // then on return it will point to a vector containing the found elements.
    // \param key        Keyword for element tag or attribute name
    // \param list       Vector of strings (output)
    // \param section    Name of a section in which to look for keyword
    // \param subsection Name of a subsection in which to look for keyword
    // \param treeTop    Top node of search tree
    // \param elemlist   Vector containing the found elements
    //                   (output/optional)
    //@}


    // ***********************************************************************
    //   Private Auxilliary Methods: Obtain / Test values
    // ***********************************************************************

    //@{
    //! \name Routines to obtain/test the value of an element or attribute

    //! Obtain the value of an element that is a text node

    //! This method will return the value of an element as character array.
    //! Note that the method will only succeed if the element possesses
    //! one child, which is a text node. If this is not the case, an error
    //! will be issued. 
    //! Note that it is the caller's responsibility to free the dynamically
    //! allocate array via the FreeC() method.
    void GetElementValue( xercesc::DOMElement *elem,
                             StdVector<std::string> &values );

    //! Obtain the value of an elements attribute

    //! This method will determine the value of the specified attribute of a
    //! given element as standard string. If the element does not posses an
    //! attribute with the specified keyword, the return string is empty and
    //! the method returns false.
    //! \param element  The element
    //! \param keyword  Name of element's attribute
    //! \param attrVal  Value of element's attribute
    bool GetElementAttribute( xercesc::DOMElement* element,
                                 const std::string keyword,
                                 std::string &attrVal );

    //! Test value of an attribute

    //! This method can be used to test, if the attribute of a certain element
    //! has a specified value. If the element does not have an attribute with
    //! the specified name, an error will be reported, if the parameter
    //! failIfNoAttrib is 'true', otherwise the method will simply return
    //! 'false'.
    //! \param elem           Element whose attribute is to be tested
    //! \param attribute      Name of the attribute of interest
    //! \param value          Value to compare atrribute's value against
    //! \param failIfNoAttrib Determines behaviour, if element does not have
    //!                       an attribute with specified name
    //! \note In the case that the attribute we are looking for is the 'tag'
    //!       attribute used in a 'multiSequence' analysis, we do not go for
    //!       a one-to-one match of the attribute against the prescribed
    //!       value, but instead call MatchesTag to determine, whether this
    //!       is a matching element.
    bool AttribHasValue( xercesc::DOMElement* elem,
                         const std::string attribute,
                         const std::string value,
                         bool failIfNoAttrib = true );

    //! Test whether a refTag is contained in a tag

    //! This method is important in the context of multi-sequence analysis.
    //! In this type of analysis in each PDE of a step receives a unique
    //! refTag. In order to find the parameters for this PDE in this step
    //! the refTag is compared to the tag attribute of certain other parameter
    //! sections. There are three possible cases
    //! <center>
    //! <table border="1" width="80%" cellpadding="10">
    //! <tr><td align="center"><b>case</b></td>
    //! <td align="center"><b>meaning</b></td>
    //! </tr>
    //! <tr><td align="center">refTag = tag</td>
    //! <td>The tagged parameter section is used for this step only</td>
    //! </tr>
    //! <tr><td align="center">tag = anyTag</td>
    //! <td>The parameter section is used for all steps of the multi-sequence
    //!     analysis</td>
    //! </tr>
    //! <tr><td align="center">tag is a comma-separated list of tags and
    //!     refTag is one of them</td>
    //! <td>The tagged parameter section is used for this step among others
    //! </td>
    //! </tr>
    //! </table>
    //! </center>
    //! \n This method will compare tag and refTag and according to the above
    //! table decide, if the parameter section to which the tag belongs is
    //! to be used in the corresponding step. In the latter case it will
    //! return 'true'.
    bool MatchesTag( const std::string tag, const std::string refTag );

    //@}


    // ***********************************************************************
    //   Private class attributes
    // ***********************************************************************

    //! Parser object

    //! This is a pointer to the XercesDOMParser object for the XML parameter
    //! file. The parser is the owner of the DOMDocument it generates. Thus,
    //! deleting the parser will also delete the document, i.e. all
    //! information on the steering parameters.
    xercesc::XercesDOMParser *parser_;

    //! Parser object for defaults

    //! This is a pointer to the XercesDOMParser object for the XML file
    //! connect to specification of default values. The parser is the owner of
    //! the DOMDocument it generates. Thus, deleting the parser will also
    //! delete the document, i.e. all information on the steering parameters.
    xercesc::XercesDOMParser *parserDefaults_;

    //! Root element of document tree

    //! This is a pointer to the root element of the DOM document tree. It is
    //! not the be mistaken for the root node of the document. The former is
    //! a child of the latter.
    xercesc::DOMElement *rootElem_;

    //! Root element of document tree for defaults

    //! This is a pointer to the root element of the DOM document tree
    //! containing the default values. It is not the be mistaken for the root
    //! node of the document. The former is a child of the latter.
    xercesc::DOMElement *rootElemDefaults_;

    //! Turns on/off verbosity of class methods

    //! Setting this class attribute to true will turn on additional status
    //! messages by the methods of the class making them more verbose. This
    //! is intended for debugging purposes only. The attribute is assigned
    //! in the constructor. It is set to 'true', if the macro DEBUG_SPHANDLER
    //! is defined and to 'false' otherwise.
    bool beVerbose_;

    //! Flag indicating if a secondary default-xml file should be used
    bool useDefaults_;

  };

}

#endif

#endif
