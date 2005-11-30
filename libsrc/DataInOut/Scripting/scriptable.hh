#ifndef CFS_SCRIPTABLE_HH
#define CFS_SCRIPTABLE_HH

#include <string>
#include "Utils/StdVector.hh"

namespace CoupledField
{

  
  // Base class for implementing an interface to a central scripting instance.
  class Scriptable {

  public: 
    
    //! Destructor
    virtual ~Scriptable();

    //! Central method vor evaluating a given scripting command

    //! This method evaluates the given arguments, beginning from an offset 
    //! prescribed by argOffset. If it is successful, it returns TRUE and the
    //! as a vector of strings.
    //! \param args Vector of arguments in string format to be evaluated
    //! \param argOffset Offset for starting position in args vector
    //! \param retVal Vector of return values in string format
    //! \return TRUE, if evaluation was successful
    virtual Boolean Script_Eval( const StdVector<std::string> & args,
                                 UInt & argOffset,
                                 StdVector<std::string> & retVal) = 0;

    //! Get list of all available commands of this object

    //! This method returns a list of all available scripting commands
    //! offered by the particular class.
    //! \param commands Vector of available commands of this class
    //! \param argOffset Offset for start position in args vector
    virtual void Script_GetCommands( StdVector<std::string> & commands,
                                     UInt & argOffset);
    
    //! Get last error message

    //! In case of an error this method returns the error message and
    //! resets it
    //! \param errMsg Erro message
    static void Script_GetError( std::string & errMsg );
    
  protected: 
    
    //! Private constructor, since this class is merely an interface class
    Scriptable();
    
    //! Stream for error messages
    static std::stringstream errMsg_;

  private:
    
  };

} // end of namespace

#endif
