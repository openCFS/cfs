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

namespace CoupledField
{

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
  //! the string conersion routines. Also it is not yet clear, whether we are
  //! expected to free the results from getElementsByTagName() and the like, or
  //! if this is even forbidden.
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
    XMLParamHandler( const char *fname );

    //! Default destructor

    //! The destructor is deep in the sense that it will delete the parser
    //! and the DOM tree it owns.
    ~XMLParamHandler();


    // ========================================================================
    // Simple query methods without side constraints
    // ========================================================================

    //@{
    //! \name Simple query methods without side constraints

    void Get( const std::string key, std::string &value,
	      const std::string section="",
	      const std::string subsection="" );

    void Get( const std::string key, Integer &value,
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

    void Get( const StdVector<std::string> &keyVec,
	      std::string &value );

    void Get( const StdVector<std::string> &keyVec,
	      Double &value );

    void Get( const StdVector<std::string> &keyVec,
	      Integer &value );
    //@}


    // ========================================================================
    // Query methods allowing side-constraints and expecting vectors as input
    // ========================================================================

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
    //! The following is an example set of search parameters. We are trying to
    //! determine here the analysis type for the mechanic PDE in second step of
    //! a multi-sequence analysis:
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
    //@}


    // ========================================================================
    // Specialised query methods
    // ========================================================================

    //@{
    //! \name Specialised query methods
    void GetPDEList( StdVector<std::string> &list );

    void GetIterCoupledPDEList( StdVector<std::string> &list);
    
    void GetCoilList( StdVector<std::string> &list,
		      const std::string pde = "" );

    void GetCoilType( std::string &coilType, const std::string coilName,
		      const std::string pde ="" );

    Boolean IsSet( const std::string key,
		   const std::string section = "",
		   const std::string subsection = "" );

    Boolean HasValue( const std::string key,
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

    //! This method generates a validating parser and parses the XML input
    //! file with filename xmlFile. A pointer to the generated parser is
    //! assigned to the parser input parameter and a pointer to the root
    //! element of the DOM tree representing the contents of the XML file
    //! is returned.
    xercesc::DOMElement* ParseFile( xercesc::XercesDOMParser **parser,
				    const char *xmlFile );

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
    //! 'subsection' arguments corresponding keyVec, attrVec and valVec vectors
    void GenerateSearchParams( const std::string key,
			       const std::string section,
			       const std::string subsection,
			       StdVector<std::string> &keyVec,
			       StdVector<std::string> &attrVec,
			       StdVector<std::string> &valVec );


    // ************************************************************************
    //   Private Auxilliary Methods: Problem Handlers
    // ************************************************************************

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
    //! no default value was found either. The method assembles a corresponding
    //! error message and sends this to the error method of the WriteInfo
    //! class.
    void NoMatchErrorReporter( const StdVector<std::string> &keyVec,
			       const StdVector<std::string> &attrVec,
			       const StdVector<std::string> &valVec );

    //! Determine default value for a parameter

    //! This method will try to determine, whether for a given set of search
    //! parameters a default value exists. If this is the case, then the method
    //! will return TRUE and defaultValue will contain the default value of the
    //! parameter in string format. If there is no default parameter specified,
    //! then the method will return FALSE and defaultValue will be set to the
    //! empty string.
    Boolean CheckForDefault(  const StdVector<std::string> &keyVec,
			      const StdVector<std::string> &attrVec,
			      const StdVector<std::string> &valVec,
			      std::string &defaultValue );

    //@}


    // ************************************************************************
    //   Private Auxilliary Methods: Conversion Routines for Strings
    // ************************************************************************

    //@{
    //! \name Conversion routines for string representations:

    //! Transcodes a string to native code-page

    //! The C++ DOM of xerces uses the plain, null-terminated (XMLCh *) utf-16
    //! strings as the String type. This method converts such a string into a
    //! standard Char array. The returned buffer is dynamically allocated and
    //! it is the responsibility of the caller to delete it when not longer
    //! needed.
    Char* X2C( const XMLCh *const toTranscode ) {
      return xercesc::XMLString::transcode( toTranscode );
    };

    //! Transcodes a string from native code-page

    //! The C++ DOM of xerces uses the plain, null-terminated (XMLCh *) utf-16
    //! strings as the String type. This method converts a standard Char array
    //! into such a string. The returned buffer is dynamically allocated and
    //! it is the responsibility of the caller to delete it when not longer
    //! needed.
    XMLCh* C2X( const Char *const toTranscode ) {
      return xercesc::XMLString::transcode( toTranscode );
    };

    //! Transcodes an STL string to native code-page

    //! The C++ DOM of xerces uses the plain, null-terminated (XMLCh *) utf-16
    //! strings as the String type. This method converts such a string into a
    //! standard STL string. The returned string is dynamically allocated and
    //! it is the responsibility of the caller to delete it when not longer
    //! needed.
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
    //@}


    // ************************************************************************
    //   Private Auxilliary Methods: Conversion Routines for Nodes
    // ************************************************************************

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


    // ************************************************************************
    //   Private Auxilliary Methods: Search Routines
    // ************************************************************************

    //@{
    //! \name Search routines:

    // Search (restricted) tree for attributes matching keyword 

    // The method will try to find all attributes whose name matches the
    // specified keyword in the parameter tree returning a vector of pointers
    // to the matching attributes. The search can be restricted to subtrees,
    // called sections, subsections and so on, which themselves match certain
    // keywords. The method works in a recursive fashion. If no match at all
    // is found, an empty vector is returned. The method relies on
    // FindMatchingElements.
    // \note The method requires to specify at least a section keyword. This
    //       is taken as identifier for the element, whose attributes we are
    //       looking for.
    // \param key      Keyword describing the desired attribute
    // \param keys     A vector of strings, containing the keywords
    //                 restricting the search tree. The order of the keywords
    //                 is section, subsection, ...
    // \param treeTop  Top node  of search tree.
    // \param elemlist Optionally we can return a vector of the elements for
    //                 which there have been matches.

    // Search (restricted) tree for elements matching keyword 

    // The method will try to find all elements whose tag matches the
    // specified keyword in the parameter tree returning a vector of pointers
    // to the matching elements. The search can be restricted to subtrees,
    // called sections, subsections and so on, which themselves match certain
    // keywords. The method works in a recursive fashion. If no match at all
    // is found, a NULL pointer will be returned.
    // \param keys       A vector of strings, containing the keywords. The
    //                   Order of the keywords is section, subsection, ...
    //                   up to the final keyword for the element tag.
    // \param treeTop    Top node of the subtree that is currently searched
    //                   through in the recursion.
    // \param curdepth   Current depth of recursion. This is used to identitfy
    //                   the current keyword.

    //! Under development
    StdVector<xercesc::DOMElement *>*
    FindMatchingElements( StdVector<std::string> &keys,
			  StdVector<std::string> &attribs,
			  StdVector<std::string> &aValues,
			  xercesc::DOMElement *treetop,
			  unsigned int curdepth );

    //! Under development
    StdVector<xercesc::DOMAttr *>*
    FindMatchingAttributes( std::string attr_key,
			    StdVector<std::string> &keys,
			    StdVector<std::string> &attribs,
			    StdVector<std::string> &aValues,
			    xercesc::DOMElement *treetop );

    //! Under development
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


    // ************************************************************************
    //   Private Auxilliary Methods: Obtain / Test values
    // ************************************************************************

    //@{
    //! \name Routines to obtain/test the value of an element or attribute

    //! Obtain the value of an element that is a text node

    //! This method will return the value of an element as character array.
    //! Note that the method will only succeed if the element possesses
    //! one child, which is a text node. If this is not the case, an error
    //! will be issued.
    Char* GetElementValue( xercesc::DOMElement *elem );

    //! Obtain the value of an elements attribute

    //! This method will determine the value of the specified attribute of a
    //! given element as standard string. If the element does not posses an
    //! attribute with the specified keyword, the return string is empty and
    //! the method returns FALSE.
    //! \param element  The element
    //! \param keyword  Name of element's attribute
    //! \param attrVal  Value of element's attribute
    Boolean GetElementAttribute( xercesc::DOMElement* element,
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
    bool AttribHasValue( xercesc::DOMElement* elem,
			 const std::string attribute,
			 const std::string value,
			 bool failIfNoAttrib = true );

    //@}


    // ************************************************************************
    //   Private class attributes
    // ************************************************************************

    //! Parser object

    //! This is a pointer to the XercesDOMParser object for the XML parameter
    //! file. The parser is the owner of the DOMDocument it generates. Thus,
    //! deleting the parser will also delete the document, i.e. all information
    //! on the steering parameters.
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

    //! Location of main Schema file for validated parsing

    //! We use a validating parser. Thus, we must specify a Schema against
    //! which the input xml-file can be tested. This attribute stores the
    //! location of the main Schema definition file. Its value is derived from
    //! the value of the XMLSCHEMA macro as
    //! <b>http://www.cfs++.org XMLSCHEMA/CFS.xsd</b>.
    std::string cfsSchema_;

    //! Location of XML file for handling of default parameters

    //! In the current implementation we specify default values for parameters
    //! in the main Schema description. However, due to the fact that default
    //! values for optional elements are only incorporated into the parsing
    //! result, if the element appears and is empty, but not, if the element
    //! does not appear at all, we make use of a sort of skeleton xml file
    //! for forcing appearance of the default values. This attribute stores
    //! the location of this defaults file which is parsed to generate the
    //! default parameter tree, see also rootElemDefaults_. The value of this
    //! attribute is derived from the macro XMLSCHEMA as
    //! <b>XMLSCHEMA/Defaults/CFS++Defaults.xml</b>.
    std::string cfsDefaults_;

  };

}

#endif

#endif
