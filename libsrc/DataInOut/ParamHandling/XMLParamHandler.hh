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
  //! \todo Fix the small memory leaks introduced by the string conersion
  //! routines. Probably we also must free all results from
  //! getElementsByTagName() and the like.
  //! \note With respect to function tracing the policy of this class is the
  //! following. All public methods use ENTER_FCN, while all other methods,
  //! besides those defined in the header itself, use ENTER_IFCN. The latter
  //! do not appear in the trace log at all.
  class XMLParamHandler : public BaseParamHandler
  {
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

    //! Get string-value for a keyword

    //! The method will try to find the specified keyword in the parameter tree
    //! returning the found value as a string. The search can be restricted to
    //! a certain section and subsection.
    //! \param key        Keyword
    //! \param value      String (output)
    //! \param section    Name of a section in which to look for keyword
    //!                   (optional)
    //! \param subsection Name of a subsection in which to look for keyword
    //!                   (optional)
    void Get( const std::string key, std::string &value,
	      const std::string section="",
	      const std::string subsection="" );

    //! Get int-value for a keyword

    //! The method will try to find the specified keyword in the parameter tree
    //! returning the found value as an Integer. The search can be restricted
    //! to a certain section and subsection.
    //! \param key        Keyword
    //! \param value      Integer (output)
    //! \param section    Name of a section in which to look for keyword
    //!                   (optional)
    //! \param subsection Name of a subsection in which to look for keyword
    //!                   (optional)
    void Get( const std::string key, Integer &value,
	      const std::string section="", const std::string subsection="" );

    //! Get double-value for a keyword

    //! The method will try to find the specified keyword in the parameter tree
    //! returning the found value as a Double. The search can be restricted
    //! to a certain section and subsection.
    //! \param key        Keyword
    //! \param value      Double (output)
    //! \param section    Name of a section in which to look for keyword
    //!                   (optional)
    //! \param subsection Name of a subsection in which to look for keyword
    //!                   (optional)
    void Get( const std::string key, Double &value,
	      const std::string section="", const std::string subsection="" );

    //! Get string-value for a element with certain attribute

    //! The method will try to find the specified keyword in the parameter tree
    //! returning the found value as a string, if the corresponding element
    //! has an attribute of the specified name which has the specified value.
    //! The search can be restricted to a certain section and subsection.
    //! If the applyToElem input argument is set to false, then the
    //! attribute/value test is applied not to the element itself, but to its
    //! parent.
    //! \param key          Keyword
    //! \param value        String (output)
    //! \param attribute    Name of attribute of element
    //! \param aValue       Value to test attribute's value against
    //! \param applyToElem  Specifies which element to test
    //! \param section      Name of a section in which to look for keyword
    //!                     (optional)
    //! \param subsection   Name of a subsection in which to look for keyword
    //!                     (optional)
    void CGet( const std::string key,
	       std::string &value,
	       const std::string attribute,
	       const std::string aValue,
	       Integer applyToElem,
	       const std::string section,
	       const std::string subsection );

    //! Get Double-value for a element with certain attribute

    //! The method will try to find the specified keyword in the parameter tree
    //! returning the found value as a Double, if the corresponding element
    //! has an attribute of the specified name which has the specified value.
    //! The search can be restricted to a certain section and subsection.
    //! If the applyToElem input argument is set to false, then the
    //! attribute/value test is applied not to the element itself, but to its
    //! parent.
    //! \param key          Keyword
    //! \param value        Double (output)
    //! \param attribute    Name of attribute of element
    //! \param aValue       Value to test attribute's value against
    //! \param applyToElem  Specifies which element to test
    //! \param section      Name of a section in which to look for keyword
    //!                     (optional)
    //! \param subsection   Name of a subsection in which to look for keyword
    //!                     (optional)
    void CGet( const std::string key,
	       Double &value,
	       const std::string attribute,
	       const std::string aValue,
	       Integer applyToElem,
	       const std::string section,
	       const std::string subsection );

    //! Get Integer-value for a element with certain attribute

    //! The method will try to find the specified keyword in the parameter tree
    //! returning the found value as a Integer, if the corresponding element
    //! has an attribute of the specified name which has the specified value.
    //! The search can be restricted to a certain section and subsection.
    //! If the applyToElem input argument is set to false, then the
    //! attribute/value test is applied not to the element itself, but to its
    //! parent.
    //! \param key          Keyword
    //! \param value        Integer (output)
    //! \param attribute    Name of attribute of element
    //! \param aValue       Value to test attribute's value against
    //! \param applyToElem  Specifies which element to test
    //! \param section      Name of a section in which to look for keyword
    //!                     (optional)
    //! \param subsection   Name of a subsection in which to look for keyword
    //!                     (optional)
    void CGet( const std::string key,
	       Integer &value,
	       const std::string attribute,
	       const std::string aValue,
	       Integer applyToElem,
	       const std::string section,
	       const std::string subsection );

    //! Get a list of strings matching a keyword

    //! The method will try to find the specified keyword in the parameter
    //! tree and will return the values of the respective elements or
    //! attributes in a vector of strings. The search can be restricted to
    //! certain subtrees by specifying keywords for section and subsection.
    //! The method will return an empty vector, if there is no match at all.
    //! It will issue an error message, if there are matches for both, elements
    //! and attributes.
    //! \param key        Keyword
    //! \param list       Vector of strings (output)
    //! \param section    Name of a section in which to look for keyword
    //!                   (optional)
    //! \param subsection Name of a subsection in which to look for keyword
    //!                   (optional)
    void GetList( const std::string key,
		  std::vector<std::string> &list,
		  const std::string section = "",
		  const std::string subsection = "" );

    //! Get a list of Doubles matching a keyword

    //! The method will try to find the specified keyword in the parameter
    //! tree and will return the values of the respective elements or
    //! attributes in a vector of Doubles. The search can be restricted to
    //! certain subtrees by specifying keywords for section and subsection.
    //! The method will return an empty vector, if there is no match at all.
    //! It will issue an error message, if there are matches for both, elements
    //! and attributes.
    //! \param key        Keyword
    //! \param list       Vector of Doubles (output)
    //! \param section    Name of a section in which to look for keyword
    //!                   (optional)
    //! \param subsection Name of a subsection in which to look for keyword
    //!                   (optional)
    void GetList( const std::string key, std::vector<Double> &list,
		  const std::string section = "",
		  const std::string subsection = "" );

    //! Get a list of strings for keyword and elements with certain attribute

    //! The method will try to find the specified keyword in the parameter
    //! tree. Once found, it tests, whether the corresponding elements have
    //! a specified value for a specified attribute. For these elements it
    //! will return the values of the respective elements or of the attributes
    //! matching the keyword. The search can be restricted to certain subtrees
    //! by specifying keywords for section and subsection.
    //! The method will return an empty vector, if there is no match at all.
    //! It will issue an error message, if there are matches for both, elements
    //! and attributes, or if one of the found elements does not have an
    //! attribute of the specified type. If the applyToElem input argument is
    //! set to false, then the attribute/value test is applied not to the
    //! element itself, but to its parent.
    //! \param key          Keyword
    //! \param list         Vector of strings (output)
    //! \param attribute    Name of attribute of element
    //! \param value        Value to test attribute's value against
    //! \param applyToElem  Specifies which element to test
    //! \param section      Name of a section in which to look for keyword
    //!                     (optional)
    //! \param subsection   Name of a subsection in which to look for keyword
    //!                     (optional)
    void CGetList( const std::string key,
		   std::vector<std::string> &list,
		   const std::string attribute,
		   const std::string value,
		   Integer applyToElem,
		   const std::string section,
		   const std::string subsection );

    //! Obtain list of PDEs defined in parameter file

    //! This method queries the parameter object for a list of all PDEs defined
    //! in the parameter file. The list is returned as a vector of standard
    //! strings.
    void GetPDEList( std::vector<std::string> &list );

    //! Obtain list of coils defined in parameter file

    //! This method queries the parameter object for a list of all coils
    //! defined in the parameter file. The list is returned as a vector of
    //! standard strings. The method will return an empty vector, if there are
    //! no matches. By specifying the optional pde input parameter the search
    //! can be restricted to a certain PDE entry in the pdeList section.
    void GetCoilList( std::vector<std::string> &list,
		      const std::string pde = "" );

    //! Obtain the type of a certain coils

    //! This method queries the parameter object for the type of a coils.
    //! The desired coil is specified via the coilName input argument.
    //! By specifying the optional pde input parameter the search can be
    //! restricted to a certain PDE entry in the pdeList section.
    void GetCoilType( std::string &coilType, const std::string coilName,
		      const std::string pde ="" );

    //! Query the on/off status of a flag/switch

    //! The method will search the parameter tree for the parameter matching
    //! the keyword. If the keyword is found, it will be compared to the token
    //! "yes". If the keyword is found and its value matches, then the method
    //! will return TRUE. If the keyword is not found or it does not match,
    //! it will return FALSE. The search can be restricted to a subtree, by
    //! specifying a section and subsection argument.
    //! \param key        Keyword
    //! \param section    Name of a section in which to look for keyword
    //!                   (optional)
    //! \param subsection Name of a subsection in which to look for keyword
    //!                   (optional)
    virtual Boolean IsSet( const std::string key,
			   const std::string section = "",
			   const std::string subsection = "" );

    //! Query whether a parameter has a certain value

    //! The method will search the parameter tree for the parameter matching
    //! the keyword. If the keyword is found, its value will be compared to the
    //! value given by the value input parameter. If the keyword is found and
    //! its value matches, then the method will return TRUE. If the keyword
    //! is not found (and there is also no default), or it does not match,
    //! it will return FALSE. The search can be restricted to a subtree, by
    //! specifying a section and subsection
    //! argument.
    //! \param key        Keyword
    //! \param section    Name of a section in which to look for keyword
    //!                   (optional)
    //! \param subsection Name of a subsection in which to look for keyword
    //!                   (optional)
    //! \param value      String against which to test value of parameter
    //!                   (optional)
    virtual Boolean HasValue( const std::string key,
			      const std::string value,
			      const std::string section = "",
			      const std::string subsection = "" );

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


    // ************************************************************************
    //   Private Auxilliary Methods: Problem Handlers
    // ************************************************************************

    //@{
    //! \name Routines for handling problems and default parameters

    //! Error handler in case we have no unique match

    //! This method is called by those get routines, which expect to obtain a
    //! unqiue match, if a match is found but it is not unique. The error
    //! handler will issue a corresponding error message and terminate
    //! programm execution.
    void MultipleMatchHandler( const std::string key,
			       const std::string section,
			       const std::string subsection,
			       const unsigned int nmatches );

    //! Error handler in case we have no match at all

    //! This method is called by those get routines, which expect to obtain a
    //! unqiue match, if no match is found and we are looking for a string
    //! value.
    void NoMatchHandler( std::string &value, const std::string key,
			 const std::string section,
			 const std::string subsection );

    //! Error handler in case we have no match at all

    //! This method is called by those get routines, which expect to obtain a
    //! unqiue match, if no match is found and we are looking for a Double
    //! value.
    void NoMatchHandler( Double &value, const std::string key,
			 const std::string section,
			 const std::string subsection );

    //! Error handler in case we have no match at all

    //! This method is called by those get routines, which expect to obtain a
    //! unqiue match, if no match is found and we are looking for an Integer
    //! value.
    void NoMatchHandler( Integer &value, const std::string key,
			 const std::string section,
			 const std::string subsection );

    //! Error handler in case we have no match at all

    //! This method is called by the CGet routine, which expects to obtain a
    //! unqiue match, if no match is found and we are looking for a string
    //! value.
    void NoMatchHandler( std::string &value,
			 const std::string key,
			 const std::string attribute,
			 const std::string aValue,
			 Integer applyToElem,
			 const std::string section,
			 const std::string subsection );

    //! Error handler in case we have no match at all

    //! This method is called by the CGet routine, which expects to obtain a
    //! unqiue match, if no match is found and we are looking for an Integer
    //! value.
    void NoMatchHandler( Integer &value,
			 const std::string key,
			 const std::string attribute,
			 const std::string aValue,
			 Integer applyToElem,
			 const std::string section,
			 const std::string subsection );

    //! Error handler in case we have no match at all

    //! This method is called by the CGet routine, which expects to obtain a
    //! unqiue match, if no match is found and we are looking for a Double
    //! value.
    void NoMatchHandler( Double &value,
			 const std::string key,
			 const std::string attribute,
			 const std::string aValue,
			 Integer applyToElem,
			 const std::string section,
			 const std::string subsection );

    //! Error reporter for the case that we have no match at all

    //! This method is called by the methods handling the no match case, if
    //! no default value was found either. The method assembles a corresponding
    //! error message and sends this to the error method of the WriteInfo
    //! class.
    void NoMatchErrorReporter( const std::string key,
			       const std::string section,
			       const std::string subsection );

    //! Determine default value for a parameter

    //! This method will try to determine, whether for a given parameter key
    //! a default value exists. If this is the case, then the method will
    //! return TRUE and defaultValue will contain the default value of the
    //! parameter in string format. If there is no default parameter specified,
    //! then the method will return FALSE and defaultValue will be set to the
    //! empty string.
    //! \param defaultValue (output) Default value of parameter, if it exists
    //! \param key          (input)  Keyword describing parameter
    //! \param section      (input)  Name of a section in which to look for
    //!                              keyword (optional)
    //! \param subsection   (input)  Name of a subsection in which to look for
    //!                              keyword (optional)
    //! \param constrained  (input)  Shall we perform a search with side
    //!                              constraints.
    //! \param attribute    (input)  Name of attribute for constraint
    //! \param aValue       (input)  Value of attribute for constraint
    //! \param applyToElem  (input)  Level of element for attribute/value test
    //!                              in tree above actual element
    //! \param rootElem     (input)  root element for search tree
    Boolean CheckForDefault( std::string &defaultValue, const std::string key,
			     const std::string section,
			     const std::string subsection,
			     bool constrained = false,
			     const std::string attribute = "",
			     const std::string aValue = "",
			     Integer applyToElem = 0 );

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

    //! Search (restricted) tree for attributes matching keyword 

    //! The method will try to find all attributes whose name matches the
    //! specified keyword in the parameter tree returning a vector of pointers
    //! to the matching attributes. The search can be restricted to subtrees,
    //! called sections, subsections and so on, which themselves match certain
    //! keywords. The method works in a recursive fashion. If no match at all
    //! is found, an empty vector is returned. The method relies on
    //! FindMatchingElements.
    //! \note The method requires to specify at least a section keyword. This
    //!       is taken as identifier for the element, whose attributes we are
    //!       looking for.
    //! \param key      Keyword describing the desired attribute
    //! \param keys     A vector of strings, containing the keywords
    //!                 restricting the search tree. The order of the keywords
    //!                 is section, subsection, ...
    //! \param treeTop  Top node  of search tree.
    //! \param elemlist Optionally we can return a vector of the elements for
    //!                 which there have been matches.
    std::vector<xercesc::DOMAttr*>*
    FindMatchingAttributes( std::string key, std::vector<std::string> &keys,
			    xercesc::DOMElement *treeTop,
			    std::vector<xercesc::DOMElement*> *elemlist=NULL );

    //! Search (restricted) tree for elements matching keyword 

    //! The method will try to find all elements whose tag matches the
    //! specified keyword in the parameter tree returning a vector of pointers
    //! to the matching elements. The search can be restricted to subtrees,
    //! called sections, subsections and so on, which themselves match certain
    //! keywords. The method works in a recursive fashion. If no match at all
    //! is found, a NULL pointer will be returned.
    //! \param keys       A vector of strings, containing the keywords. The
    //!                   Order of the keywords is section, subsection, ...
    //!                   up to the final keyword for the element tag.
    //! \param treeTop    Top node of the subtree that is currently searched
    //!                   through in the recursion.
    //! \param curdepth   Current depth of recursion. This is used to identitfy
    //!                   the current keyword.
    std::vector<xercesc::DOMElement*>*
    FindMatchingElements( std::vector<std::string> &keys,
			  xercesc::DOMElement *treeTop,
			  unsigned int curdepth );

    //! Central auxilliary search method for get functions

    //! This is the central auxilliary search method used by the public get
    //! methods. It will try to find the specified keyword in the parameter
    //! tree and returning the values of the respective elements or attributes
    //! in a vector of strings. The search can be restricted to certain
    //! subtrees by specifying keywords for section and subsection. The method
    //! will return an empty vector, if there is no match at all. It will issue
    //! an error message, if there are matches for both, elements and
    //! attributes. The final argument to the method is a pointer to a vector.
    //! If this pointer is not the NULL pointer on input (which is the default)
    //! then on return it will point to a vector containing the found elements.
    //! \param key        Keyword for element tag or attribute name
    //! \param list       Vector of strings (output)
    //! \param section    Name of a section in which to look for keyword
    //! \param subsection Name of a subsection in which to look for keyword
    //! \param treeTop    Top node of search tree
    //! \param elemlist   Vector containing the found elements
    //!                   (output/optional)
    void FindAllMatches( const std::string key,
			 std::vector<std::string> &list,
			 const std::string section,
			 const std::string subsection,
			 xercesc::DOMElement *treeTop,
			 std::vector<xercesc::DOMElement*> *elemlist=NULL );

    //! Get a list of strings for keyword and elements with certain attribute

    //! The method will try to find the specified keyword in the parameter
    //! tree. Once found, it tests, whether the corresponding elements have
    //! a specified value for a specified attribute. For these elements it
    //! will return the values of the respective elements or of the attributes
    //! matching the keyword. The search can be restricted to certain subtrees
    //! by specifying keywords for section and subsection.
    //! The method will return an empty vector, if there is no match at all.
    //! It will issue an error message, if there are matches for both, elements
    //! and attributes, or if one of the found elements does not have an
    //! attribute of the specified type. If the applyToElem input argument is
    //! set to false, then the attribute/value test is applied not to the
    //! element itself, but to its parent.
    //! \param key          Keyword
    //! \param list         Vector of strings (output)
    //! \param attribute    Name of attribute of element
    //! \param value        Value to test attribute's value against
    //! \param applyToElem  Specifies which element to test
    //! \param section      Name of a section in which to look for keyword
    //!                     (optional)
    //! \param subsection   Name of a subsection in which to look for keyword
    //!                     (optional)
    void CFindAllMatches( const std::string key,
			  std::vector<std::string> &list,
			  const std::string attribute,
			  const std::string value,
			  Integer applyToElem,
			  const std::string section,
			  const std::string subsection,
			  xercesc::DOMElement *rootElem );

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
    //! the specified name, an error will be reported.
    //! \param elem      Element whose attribute is to be tested
    //! \param attribute Name of the attribute of interest
    //! \param value     Value to compare atrribute's value against
    bool AttribHasValue( xercesc::DOMElement* elem,
			 const std::string attribute,
			 const std::string value );

    //@}


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
  };

}

#endif

#endif
