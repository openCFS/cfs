#ifndef CFS_TCL_MESSENGER_HH
#define CFs_TCL_MESSENGER_HH

#include <tcl.h>
#include <map>

#include "cfsmessenger.hh"


namespace CoupledField {

  //! Central instance of message dispatcher for CFS++
  class TCL_CFSMessenger : public CFSMessenger {
    
  public: 

    //! Constructor
    //! \param scriptFileName Script file to be evaluated
    TCL_CFSMessenger( const std::string & scriptFileName );
    
    //! Destructor
    virtual ~TCL_CFSMessenger();


    //! This method triggers the call of a related event procedure.
    //! \param event Type of event function to be called
    //! \param context Additional read only parameters to be passes to the 
    //!  function
    //! \return TRUE, if function exists and is successfully executed
    Boolean TriggerEvent( const EventType event, 
                          const StdVector<std::string> & context);
    
    // ===================================================
    // T C L  -  I N T E R F A C E 
    // ===================================================

    //! TLC version of set function of base class
    static int TCL_Set(ClientData clientdata, Tcl_Interp *interp,
                       int argc, const char *argv[]);

    //! TCL version of get function of base class
    static int TCL_Get(ClientData clientdata, Tcl_Interp *interp,
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
    
  };
  
  
} // end of namespace

#endif
 

