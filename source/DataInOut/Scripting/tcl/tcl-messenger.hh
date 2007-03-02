#ifndef CFS_TCL_MESSENGER_HH
#define CFS_TCL_MESSENGER_HH

#include <tcl.h>
#include <map>

#include "DataInOut/Scripting/cfsmessenger.hh"


namespace CoupledField {

  //! Central instance of message dispatcher for CFS++
  class TCL_CFSMessenger : public CFSMessenger {
    
  public: 

    //! Constructor
    TCL_CFSMessenger( );
    
    //! Destructor
    virtual ~TCL_CFSMessenger();

    //! Trigger reading from script-file
    //! \param fileName Script file to be evaluated
    void ReadScriptFile ( const std::string & fileName );

    //! This method triggers the call of a related event procedure.
    //! \param event Type of event function to be called
    //! \param context Additional read only parameters to be passed to the 
    //!  function
    //! \return true, if function exists and is successfully executed
    bool TriggerEvent( const EventType event, 
                          const StdVector<std::string> & context);
    
    //! Trigger the writing of a warning message
    virtual void Warning( const Char * msg, const Char * const filename,
                          const UInt numline);
    
    //! Trigger the abortion of the program with a given error message
    virtual void Error( const Char * msg, const Char * const filename,
                        const UInt numline);

    // ===================================================
    // T C L  -  I N T E R F A C E 
    // ===================================================

    //! TLC version of eval function of base class
    static int TCL_CFSEval( ClientData clientdata, Tcl_Interp *interp,
                            int argc, const char *argv[]);

  protected:
    
    //! Register function with external scripting language
    void RegisterFunctions();

    //! Triggers the registration of dummy procedures for each event
    void RegisterEvents();
    
    //! Pointer to instance of tcl interpreter
    Tcl_Interp * tcl_;

    //! Complete filename of script file to be evaluated
    std::string fileName_;

    //! Map from event enumeration type to string representation
    std::map<EventType, std::string> eventNames_;
    

    //! Store current parameters for tracing purpose
    static StdVector<std::string> curParams_;

    //! Store current Event for tracing purpose
    std::string curEvent_;
  };
  
  
} // end of namespace

#endif
 

