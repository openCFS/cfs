// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef CFS_SCRIPTABLE_HH
#define CFS_SCRIPTABLE_HH

#include <map>
#include <ostream>
#include <string>
#include <utility>

#include "General/exception.hh"
#include "Utils/StdVector.hh"

namespace CoupledField
{

  class ArgList;
  //! Forward class declaration
  class BaseFctPointer;

  // Base class for implementing an interface to a central scripting instance.
  class Scriptable {

  public: 

    //! TypeDef for named function pointer map
    typedef std::map<std::string, BaseFctPointer*> FctPtMap;
   
    //! TypeDef for named argument list map
    typedef std::map<std::string, ArgList> ArgListMap;
    

    //! Destructor
    virtual ~Scriptable();

    //! Register a new function for scripting

    //! Registers a new function for scripting under a given name. 
    //! The pointer to the function has to be passed in form of a 
    //! BaseFctPointer and the needed arguments as an argList object.
    //! \param name Name of the function in the target (script) language
    //! \param fctPointer Function pointer object to 'void f()'-function 
    //! \param argList Argument list which encapsulated the needed parameters 
    //!        of the function
    //! \note The pointer to BaseFctPointer gets deleted in the end by 
    //! Scriptable, so do not try delete it outside this class!
    void Script_RegisterFct( const std::string name,
                             BaseFctPointer  * fctPointer,
                             const ArgList & argList );


    //! Get the arguments of the function currently executing

    //! Returns the argument list of the function which is
    //! currently executed. If no function is executed, NULL is
    //! returned.
    ArgList * Script_GetArgList( );

    
    //! Returns true, if currently a scripting function is evaluating
    bool Script_IsExecuting() { return isExecuting_; }

    
    //! Central method vor evaluating a given scripting command

    //! This method evaluates the given arguments, beginning from an offset 
    //! prescribed by argOffset. If it is successful, it returns true and the
    //! as a vector of strings.
    //! \param args Vector of arguments in string format to be evaluated
    //! \param argOffset Offset for starting position in args vector
    //! \param retVal Vector of return values in string format
    //! \return true, if evaluation was successful
    bool Script_Eval( const StdVector<std::string> & args,
                         unsigned int & argOffset,
                         StdVector<std::string> & retVal);

    //! Get list of all available commands of this object

    //! This method returns a list of all available scripting commands
    //! offered by the particular class.
    //! \param commands Vector of available commands of this class
    void Script_GetCommands( StdVector<std::string> & commands);
    
    //! Get last error message

    //! In case of an error this method returns the error message and
    //! resets it
    //! \param errMsg Error message
    static void Script_GetError( std::string & errMsg );
    
  protected: 

    //! Wrapper for getting all scripting commands of concrete class
    void Wrap_GetCommands();
    
    //! Private constructor (class is abstract and must be inherited)
    Scriptable();
    
    //! Stream for error messages
    static std::stringstream errMsg_;

    //! Flag, if script is executing at the moment
    bool isExecuting_;

    //! Map for assigning function-names to pointers
    FctPtMap fctPointers_;
    
    //! Map for assigning each function-name an argumentList
    ArgListMap argLists_;

    //! Pointer to current argument list
    ArgList * currentArgs_;

  private:
    
  };


  // =========================================================================
  //    Useful Defines for Parameter Handling
  // =========================================================================
 
  //! Short form for getting a previously registered parameter of ArgList
#define SCRIPT_GET(TYPE, NAME)                  \
  TYPE NAME;                                    \
  Script_GetArgList()->GetParam( #NAME, NAME);

  //! Short form for getting a reference to the return values of ArgList
#define SCRIPT_RETVAL                           \
  Script_GetArgList()->GetRetVal()

  // =========================================================================
  //    Argument Handling of Functions
  // =========================================================================

  //! This class defines a general definition of a argument / parameter list
  class ArgList{
  public:

    //! Define type of arguments
    typedef enum { STRING, DOUBLE, UINT, INT,  
                   //VEC_UINT, VEC_INT, VEC_DOUBLE, -- killme TODO: Fabian 
                   STDVEC_STR, STDVEC_DOUBLE, STDVEC_UINT, STDVEC_INT } ParamType;

    //! Register a new paramter, which can be later set/get
    void RegisterParam( const std::string name, 
                        ParamType paramType );

    //! Set all parameters at once
    void SetParams( const StdVector<std::string> params ); 

    //! Get vector with return values
    StdVector<std::string> & GetRetVal();

    //! Return number of mandatory arguments
    unsigned int GetNumParams();

    //! Get Signature of argument list
    void GetSignature( std::string & signature );

    //! Conversion function ParamType -> string
    std::string ParamType2String( const ParamType type );
    
  private:
    
    //! Vector for order of parameters
    StdVector<std::string> orderedParams_;

    //! Map for defining the type of each parameter
    std::map<std::string, ParamType> paramMap_;

    //! Vector with return values
    StdVector<std::string> returnVals_;

#define DEFINE_PARAM(TYPE_ENUM, TYPE)                                   \
    private:                                                            \
    std::map<std::string, TYPE > TYPE_ENUM ## pool_;                    \
  public:                                                               \
        void SetParam( const std::string key, TYPE & param ) {          \
          if ( paramMap_.find(key) == paramMap_.end() ) {               \
            EXCEPTION( "Parameter '" << key << "' is not registered!" );\
          }                                                             \
          TYPE_ENUM ## pool_[key] = param; }                            \
        void GetParam( const std::string key, TYPE & out ) const {      \
          if ( paramMap_.find(key) == paramMap_.end() ) {               \
            EXCEPTION("Parameter '" << key << "' is not registered!");  \
          }                                                             \
          std::map<std::string, TYPE >::const_iterator it;              \
          it = TYPE_ENUM ## pool_.find(key);                            \
          if ( it ==  TYPE_ENUM ## pool_.end() ) {                      \
            EXCEPTION("Parameter '" << key << "' is not set!" );        \
          }                                                             \
          out = (*it).second;                                           \
        }                                                                   

    
    DEFINE_PARAM(STRING,std::string);
    DEFINE_PARAM(DOUBLE,double);
    DEFINE_PARAM(UINT,unsigned int);
    DEFINE_PARAM(INT,int);
    DEFINE_PARAM(STDVEC_STR,StdVector<std::string>);
    DEFINE_PARAM(STDVEC_DOUBLE,StdVector<double>);
    DEFINE_PARAM(STDVEC_UINT,StdVector<unsigned int>);
    DEFINE_PARAM(STDVEC_INT,StdVector<int>);
    
#undef DEFINE_PARAM

  };
  
  // =========================================================================
  //    General Pointer Class to "void f()"-type class member functions
  // =========================================================================
  
  //! Base class for pointer to 'void Class:f()'-type functions
  class BaseFctPointer {
  public:
    
    //! Constructor
    BaseFctPointer() {;}

    //! Destructor
    virtual ~BaseFctPointer() {;}
     
    //! Executes the function
    virtual void Call() = 0;
  };
  
  //! Templated class for specialized pointers to class member functions
  template<class T>
  class FctPointer : public BaseFctPointer {

  public:

    //! Define void function pointer of class T
    typedef void (T::*ptFct)();
    
    //! Constructor
    
    //! Constructor with reference to an object and a particular 
    //! member function
    //! \param obj Pointer to class instance / object
    //! \param function Pointer to 'void Class::f()' function
    FctPointer( T * obj, ptFct function ) 
      : obj_(obj), function_(function) {}

    //! Destructor
    virtual ~FctPointer() {}

    //! Execute Command
    void Call() { (obj_->*function_)(); }

  private:
    
    //! Pointer to object of class
    T * obj_;

    //! Pointer to class member function
    ptFct function_;
  };
  
  
  

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of class Scriptable
  // =========================================================================

  //! \class Scriptable
  //! 
  //! \purpose 
  //! This class serves as a generic interface for all scriptable functions.
  //! It acts like the receiver (client) of the messages sent by CFSMessenger
  //! and its derived classes.
  //! Every Class, which offers scriptable functions has to inherit this class
  //! and register its functions via the RegisterFunctions method with the 
  //! related parameters needed.
  //! 
  //! \collab 
  //! - As stated above, Scriptable receives messages from the generic 
  //! object CFSMessenger.
  //! - For storing pointers to the functions to be called, it uses the 
  //! BaseFctPointer class.
  //! - For storing the arguments needed for the scriptbale functions, it 
  //! utilizes the ArgList class.
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! - Optional function parameters with default values
  //! - Handling of real list and vector elements (up to now every list
  //! provided by the scripting-language ist treated as a string.
  //! 

  // =========================================================================
  //     Detailed description of class ArgList
  // =========================================================================

  //! \class ArgList
  //! 
  //! \purpose 
  //! This class defines an abstract interface to a general parameter / 
  //! argument list. Therefore it offers methods to register arguments
  //! with name and Type (defined by Arglist::ParamType), set the values
  //! and return the signature of the argument list. 
  //! The input parameters (always std::string) are automatically 
  //! converted to the specified type.
  //! 
  //! \collab 
  //! This class is mainly used by Scriptable for mapping associating
  //! scriptable functions with their needed arguments.
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve

  //! 
  
  // =========================================================================
  //     Detailed description of class BaseFctPointer
  // =========================================================================

  //! \class BaseFctPointer 
  //! 
  //! \purpose 
  //! This class is needed to wrap a function pointer to a general Class-
  //! member function of type 'void Class::f()'. Since the Class is part of the
  //! function pointer signature, the real pointer is stored in a derived class.
  //! By calling the method Call(), the function gets actually called.
  //! 
  //! \collab 
  //! This class is mainly used by Scriptable, which uses it to map function 
  //! names (used in the scripting-language) to C++ Class-member functions.
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve

  //! 
#endif

  

} // end of namespace

#endif
