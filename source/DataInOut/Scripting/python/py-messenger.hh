// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef CFS_PYTHON_MESSENGER_HH
#define CFS_PYTHON_MESSENGER_HH

#ifdef _POSIX_C_SOURCE
#undef _POSIX_C_SOURCE
#endif
#ifdef _XOPEN_SOURCE
#undef _XOPEN_SOURCE
#endif

#include <Python.h>

#include "DataInOut/Scripting/CFSMessenger.hh"


namespace CoupledField {

  //! Python instance of messenger object
  class PY_CFSMessenger : public CFSMessenger {
    
  public: 

    //! Constructor
    PY_CFSMessenger();

    //! Destructor
    virtual ~PY_CFSMessenger();

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
    virtual void Warning( const char * msg, const char * const filename,
                          const UInt numline);
    
    //! Trigger the abortion of the program with a given error message
    virtual void Error( const char * msg, const char * const filename,
                        const UInt numline);
    
    // ===================================================
    // P Y T H O N  -  I N T E R F A C E 
    // ===================================================
    
    //! Python version of eval function of base class
    static PyObject *
    PY_CFSEval( PyObject *self, PyObject *args);

    //! Initialization structure for function binding
    static PyMethodDef CFSMethods[];

  protected:
    
    //! Triggers the registration of dummy procedures for each event
    void RegisterEvents();

    //! Complete filename of script file to be evaluated
    std::string fileName_;

    //! Map from event enumeration type to string representation
    std::map<EventType, std::string> eventNames_;
    

    //! Store current parameters for tracing purpose
    static StdVector<std::string> curParams_;

    //! Store current Event for tracing purpose
    std::string curEvent_;

    //! Global dictionary
    PyObject * mainDict_;

  };
  
  //! Python specific module initialization function
  //PyMODINIT_FUNC initcfs();

}

#endif
