#ifndef FILE_PLAIN_XML_PARAMHANDLER
#define FILE_PLAIN_XML_PARAMHANDLER

#include "General/environment.hh"
#include "Utils/StdVector.hh"
#include <fstream>


namespace CoupledField
{

  // Forward declaration of classes
  class BaseParamHandler;

  //! Class for handling steering parameters for CFS++ simulations

  //! This class is intended for handling the input parameters steering a
  //! CFS++ simulation. The parameters are read from an XML file. 
  //!
  //! The XML file must follow the grammar defined in the standard CFS++
  //! Schema definition for parameter files. However, the underlying parser
  //! is non-validating.
  //!
  //! This class is intended as a replacement for the XMLParamHandler class
  //! for situations (read architectures and compilers) where the latter cannot
  //! be used due to its dependance on Apache's xerces library.
  class PlainXMLParamHandler : public BaseParamHandler
  {

  public:

    //! Constructor

    //! The constructor expects as input the name of the input file containing
    //! the steering parameters.
    PlainXMLParamHandler( Char *fname );

    //! Deconstructor
    ~PlainXMLParamHandler();

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
	      const std::string section="",
	      const std::string subsection="" );

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
	      const std::string section="",
	      const std::string subsection="" );

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
		  StdVector<std::string> &list,
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
    void GetList( const std::string key,
		  StdVector<Double> &list,
		  const std::string section = "",
		  const std::string subsection = "" );

    //! Obtain list of PDEs defined in parameter file

    //! This method queries the parameter object for a list of all PDEs defined
    //! in the parameter file. The list is returned as a vector of standard
    //! strings.
    void GetPDEList( StdVector<std::string> &list );

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
    Boolean IsSet( const std::string key,
		   const std::string section = "",
		   const std::string subsection = "" );

    //! Query whether a parameter has a certain value

    //! The method will search the parameter tree for the parameter matching
    //! the keyword. If the keyword is found, its value will be compared to the
    //! value given by the value input parameter. If the keyword is found and
    //! its value matches, then the method will return TRUE. If the keyword
    //! is not found, or it does not match, it will return FALSE. The search
    //! can be restricted to a subtree, by specifying a section and subsection
    //! argument.
    //! \param key        Keyword
    //! \param section    Name of a section in which to look for keyword
    //!                   (optional)
    //! \param subsection Name of a subsection in which to look for keyword
    //!                   (optional)
    //! \param value      String against which to test value of parameter
    //!                   (optional)
    Boolean HasValue( const std::string key, const std::string value,
		      const std::string section = "",
		      const std::string subsection = "" );

  private:

    //! Shortcut for standard type
    typedef std::string::size_type sType;

    //! XML-file with input parameters
    std::ifstream infile;

    //! Final position in file. Used instead of std::eof
    sType pos_end_;

    //! Search element tree for keyword

    //! This method will search the parameter tree for the element matching the
    //! keyword. If it is found, the position of the next line will be
    //! returned, otherwise, std::string::npos will be returned. By default,
    //! we start the search from the begining of the file.
    sType getposElem( const std::string keyword,
				       const sType
				       startpos = 0 );

    //! Auxillary method.

    //! This is an auxiliary method. It extracts the name of the PDE from the
    //! input string unpeeled and returns it in the output string peeled.
    void peel( const std::string unpeeled, std::string & peeled);

    //! This method checks, whether the buffer contains a tag
    Boolean is_tag( const std::string buffer);

    //! To find the final and the initial positions for required elements
    //! with keys: section, subsection, etc.
    void FindPosElems( StdVector<std::string> keys, int level,
		       StdVector<sType> & s_elems,
		       StdVector<sType> & e_elems,
		       StdVector<sType> s_section,
		       StdVector<sType> e_section);
    
    //! To find the final and the initial positions for required elements
    //! with keys: section, subsection, etc., where the last key is an
    //! attribute.
    void FindPosAttrs( StdVector<std::string> keys, int level,
		       StdVector<sType> & s_elems,
		       StdVector<sType> & e_elems,
		       StdVector<sType> s_section,
		       StdVector<sType> e_section);

    //! To find the final and initial positions of required elements with the
    //! key
    void getElems( const std::string key,
		   StdVector<sType> & s_elems,
		   StdVector<sType> & e_elems,
		   StdVector<sType> s_section,
		   StdVector<sType> e_section,
		   sType start=0,
		   const int type = 1 );

    //! To find the final and initial positions of required elements with the
    //! key
    void getAttr(const std::string attr,
		 const std::string element,
		 StdVector<sType> & s_elems,
		 StdVector<sType> & e_elems,
		 StdVector<sType> s_section,
		 StdVector<sType> e_section,
		 sType start=0);
    
    //! Find position of the input keyword key in the file
    sType findPos( const std::string key, sType start );

    //! Read values of an element
    void ReadValuesElem( StdVector<std::string> & list,
			 StdVector<sType> s_elems,
			 StdVector<sType> e_elems );
    
    //! Read attributes of an element
    void ReadAttrsElem( StdVector<std::string> & list,
			StdVector<sType> s_elems,
			StdVector<sType> e_elems );

    //! Read the position of the end of an element
    sType findEndPosElementType2( sType start );

  };

}

#endif
