#ifndef FILE_BASEPPARAMHANDLER
#define FILE_BASEPPARAMHANDLER

// Required for the CFS own data types
#include "General/environment.hh"

namespace CoupledField
{


  //! Base class for parameter handlers.

  //! This is the base class from which all parameter handlers are derived.
  //! It is a pure virtual class which defines the public query methods a
  //! parameter handler for CFS++ must provide. These are also the only
  //! methods which other CFS classes should call.
  class BaseParamHandler
  {

  public:

    //! Virtual default destructor

    //! The default destructor is made virtual
    //! for inheritance purposes.
    virtual ~BaseParamHandler(){};

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
    virtual void Get( const std::string key, std::string &value,
		      const std::string section="",
		      const std::string subsection="" ) = 0;

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
    virtual void Get( const std::string key, Integer &value,
		      const std::string section="",
		      const std::string subsection="" ) = 0;

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
    virtual void Get( const std::string key, Double &value,
		      const std::string section="",
		      const std::string subsection="" ) = 0;

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
    virtual void GetList( const std::string key,
			  std::vector<std::string> &list,
			  const std::string section = "",
			  const std::string subsection = "" ) = 0;

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
    virtual void GetList( const std::string key,
			  std::vector<Double> &list,
			  const std::string section = "",
			  const std::string subsection = "" ) = 0;

    //! Obtain list of PDEs defined in parameter file

    //! This method queries the parameter object for a list of all PDEs defined
    //! in the parameter file. The list is returned as a vector of standard
    //! strings. The method will return an empty vector, if there are no
    //! matches.
    virtual void GetPDEList( std::vector<std::string> &list ) = 0;

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
			   const std::string subsection = "" ) = 0;

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
    virtual Boolean HasValue( const std::string key,
			      const std::string value,
			      const std::string section = "",
			      const std::string subsection = "" ) = 0;
  };

}

#endif
