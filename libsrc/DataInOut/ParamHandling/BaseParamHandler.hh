#ifndef FILE_BASEPPARAMHANDLER
#define FILE_BASEPPARAMHANDLER

// Required for the CFS own data types
#include "General/environment.hh"
#include "Utils/StdVector.hh"

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

    // ========================================================================
    // SIMPLE QUERY METHODS WITHOUT SIDE CONSTRAINTS
    // ========================================================================

    //@{
    //! \name Simple query methods without side constraints
    //! These are simple query methods that do not allow to restrict the search
    //! path by side constraints and only allow at most three keywords to
    //! specify the path to the desired parameter value.

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

    //! Get unsigned int-value for a keyword

    //! The method will try to find the specified keyword in the parameter tree
    //! returning the found value as an UInt. The search can be restricted
    //! to a certain section and subsection.
    //! \param key        Keyword
    //! \param value      UInt (output)
    //! \param section    Name of a section in which to look for keyword
    //!                   (optional)
    //! \param subsection Name of a subsection in which to look for keyword
    //!                   (optional)
    virtual void Get( const std::string key, UInt &value,
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
                          StdVector<std::string> &list,
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
    virtual void GetList( const std::string key, StdVector<Double> &list,
                          const std::string section = "",
                          const std::string subsection = "" ) = 0;

    //! Get a list of Integers matching a keyword

    //! The method will try to find the specified keyword in the parameter
    //! tree and will return the values of the respective elements or
    //! attributes in a vector of Integers. The search can be restricted to
    //! certain subtrees by specifying keywords for section and subsection.
    //! The method will return an empty vector, if there is no match at all.
    //! It will issue an error message, if there are matches for both, elements
    //! and attributes.
    //! \param key        Keyword
    //! \param list       Vector of Integer (output)
    //! \param section    Name of a section in which to look for keyword
    //!                   (optional)
    //! \param subsection Name of a subsection in which to look for keyword
    //!                   (optional)
    virtual void GetList( const std::string key, StdVector<Integer> &list,
                          const std::string section = "",
                          const std::string subsection = "" ) = 0;

    //! Get a list of unsigned Integers matching a keyword

    //! The method will try to find the specified keyword in the parameter
    //! tree and will return the values of the respective elements or
    //! attributes in a vector of unsigned Integers. The search can be 
    //! restricted to certain subtrees by specifying keywords for section 
    //! and subsection.
    //! The method will return an empty vector, if there is no match at all.
    //! It will issue an error message, if there are matches for both, elements
    //! and attributes.
    //! \param key        Keyword
    //! \param list       Vector of UInt (output)
    //! \param section    Name of a section in which to look for keyword
    //!                   (optional)
    //! \param subsection Name of a subsection in which to look for keyword
    //!                   (optional)
    virtual void GetList( const std::string key, StdVector<UInt> &list,
                          const std::string section = "",
                          const std::string subsection = "" ) = 0;

    //! Search for a unique match and return value as string.

    //! The method searches the parameter file for an element or attribute and
    //! returns the value as string, if a unique match is found. Otherwise an
    //! error is reported and program execution terminated.
    //! \param keyVec     (input) vector of keywords describing path to
    //!                           desired parameter
    //! \param value      (output) value of desired parameter
    virtual void Get( const StdVector<std::string> &keyVec,
                      std::string &value ) = 0;

    //! Search for a unique match and return value as Double.

    //! The method searches the parameter file for an element or attribute and
    //! returns the value as Double, if a unique match is found. Otherwise an
    //! error is reported and program execution terminated.
    //! \param keyVec     (input) vector of keywords describing path to
    //!                           desired parameter
    //! \param value      (output) value of desired parameter
    virtual void Get( const StdVector<std::string> &keyVec,
                      Double &value ) = 0;

    //! Search for a unique match and return value as Integer.

    //! The method searches the parameter file for an element or attribute and
    //! returns the value as Integer, if a unique match is found. Otherwise an
    //! error is reported and program execution terminated.
    //! \param keyVec     (input) vector of keywords describing path to
    //!                           desired parameter
    //! \param value      (output) value of desired parameter
    virtual void Get( const StdVector<std::string> &keyVec,
                      Integer &value ) = 0;

    //! Search for a unique match and return value as unsigned Integer.

    //! The method searches the parameter file for an element or attribute and
    //! returns the value as UInt, if a unique match is found. Otherwise an
    //! error is reported and program execution terminated.
    //! \param keyVec     (input) vector of keywords describing path to
    //!                           desired parameter
    //! \param value      (output) value of desired parameter
    virtual void Get( const StdVector<std::string> &keyVec,
                      UInt &value ) = 0;
    //@}


    // ========================================================================
    // QUERY METHODS INCLUDING SIDE CONSTRAINTS
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
    virtual void GetList( const StdVector<std::string> &keyVec,
                          const StdVector<std::string> &attrVec,
                          const StdVector<std::string> &valVec,
                          StdVector<std::string> &list ) = 0;

    //! Perform search and return values of matches as Doubles.
    virtual void GetList( const StdVector<std::string> &keyVec,
                          const StdVector<std::string> &attrVec,
                          const StdVector<std::string> &valVec,
                          StdVector<Double> &list ) = 0;

    //! Perform search and return values of matches as Integers.
    virtual void GetList( const StdVector<std::string> &keyVec,
                          const StdVector<std::string> &attrVec,
                          const StdVector<std::string> &valVec,
                          StdVector<Integer> &list ) = 0;

    //! Perform search and return values of matches as unsigned Integers.
    virtual void GetList( const StdVector<std::string> &keyVec,
                          const StdVector<std::string> &attrVec,
                          const StdVector<std::string> &valVec,
                          StdVector<UInt> &list ) = 0;

    //! Search for a unique match and return value as string.
    virtual void Get( const StdVector<std::string> &keyVec,
                      const StdVector<std::string> &attrVec,
                      const StdVector<std::string> &valVec,
                      std::string &value ) = 0;

    //! Search for a unique match and return value as Double.
    virtual void Get( const StdVector<std::string> &keyVec,
                      const StdVector<std::string> &attrVec,
                      const StdVector<std::string> &valVec,
                      Double &value ) = 0;

    //! Search for a unique match and return value as Integer.
    virtual void Get( const StdVector<std::string> &keyVec,
                      const StdVector<std::string> &attrVec,
                      const StdVector<std::string> &valVec,
                      Integer &value ) = 0;

    //! Search for a unique match and return value as unsigned Integer.
    virtual void Get( const StdVector<std::string> &keyVec,
                      const StdVector<std::string> &attrVec,
                      const StdVector<std::string> &valVec,
                      UInt &value ) = 0;
    //@}


    // ========================================================================
    // SPECIALISED QUERY METHODS
    // ========================================================================

    //@{
    //! \name Specialised query methods

    //! Obtain list of PDEs defined in parameter file

    //! This method queries the parameter object for a list of all PDEs defined
    //! in the parameter file. The list is returned as a vector of standard
    //! strings. The method will return an empty vector, if there are no
    //! matches.
    virtual void GetPDEList( StdVector<std::string> &list ) = 0;

    
    //! Obtain list of iterative coupled PDEs defined in parameter file

    //! This methode queries the parameter object for a list of all PDEs
    //! which are iterative coupled. 
    //! \param list        List with name of all iterative coupling pdes
    //! \param sequenceTag Symbolic tag for multi-sequence analysis (optional)
    virtual void GetIterCoupledPDEList( StdVector<std::string> &list,
                                        const std::string sequenceTag = "")= 0;

    //! Obtain list with names of pairwise direct couplings

    //! This methode queries the parameter object for a list of all names of 
    //! pairwise direct couplings (e.g. piezo, acouMech)
    //! \param list        List with names of all pairwise direct couplings
    //! \param sequenceTag Symbolic tag for multi-sequence analysis (optional)
    virtual void GetDirectCouplingList( StdVector<std::string> &list,
                                        const std::string sequenceTag = "")= 0;
    
    //! Obtain list of coils defined in parameter file

    //! This method queries the parameter object for a list of all coils
    //! defined in the parameter file. The list is returned as a vector of
    //! standard strings. The method will return an empty vector, if there are
    //! no matches. By specifying the optional pde input parameter the search
    //! can be restricted to a certain PDE entry in the pdeList section and
    //! with the additional sequenceTag it can be restricted to a certain step
    //! in a multiSequence.
    virtual void GetCoilList( StdVector<std::string> &list,
                              const std::string pde ="",
                              const std::string sequenceTag = "") = 0;

    //! Obtain the type of a certain coils

    //! This method queries the parameter object for the type of a coils.
    //! The desired coil is specified via the coilName input argument.
    //! By specifying the optional pde input parameter the search can be
    //! restricted to a certain PDE entry in the pdeList section.
    virtual void GetCoilType( std::string &coilType,
                              const std::string coilName,
                              const std::string pde ="" ) = 0;

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
