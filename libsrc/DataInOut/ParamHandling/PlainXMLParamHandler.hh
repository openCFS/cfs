#ifndef FILE_PLAIN_XML_PARAMHANDLER
#define FILE_PLAIN_XML_PARAMHANDLER

#include "General/environment.hh"

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

    //! Obtain list of PDEs defined in parameter file

    //! This method queries the parameter object for a list of all PDEs defined
    //! in the parameter file. The list is returned as a vector of standard
    //! strings.
    void GetPDEList( std::vector<std::string> &list );

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
  };

}

#endif
