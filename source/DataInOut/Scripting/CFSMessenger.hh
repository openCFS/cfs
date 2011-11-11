// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef CFS_MESSENGER_HH
#define CFS_MESSENGER_HH

#include <map>

#include "Domain/Domain.hh"
#include "General/Environment.hh"
#include "Utils/StdVector.hh"
#include "Scriptable.hh"
namespace CoupledField {
  
  //! Central instance of message dispatcher for CFS++
  class CFSMessenger {
    
  public: 

    //! Define all available eventtypes
    typedef enum {CFS_Init, CFS_PdeInit, CFS_ReadBCs, CFS_AssembleMat,
      CFS_AssembleRhs, CFS_SetBCs, CFS_CalcResults,
      CFS_Finish} EventType;
    
    //! Constructor
    CFSMessenger();

    //! Destructor
    virtual ~CFSMessenger();

    //! Trigger reading from script-file
    virtual void ReadScriptFile ( const std::string & fileName );
    
    //! Call event procedure in target interface language

    //! This method triggers the call of a related event procedure.
    //! \param event Type of event function to be called
    //! \param context Additional read only parameters to be passes to the 
    //!  function
    //! \return true, if function exists and is successfully executed
    virtual bool TriggerEvent( const EventType event, 
                                  const StdVector<std::string> & context);
    
    //! Returns true, if a script command is evaluated at the moment
    bool IsEvaluating() {return isEvaluating_;}
    
    //! Trigger the writing of a warning message
    virtual void Warning( const char * msg, const char * const filename,
                          const UInt numline);
    
    //! Trigger the abortion of the program with a given error message
    virtual void Error( const char * msg, const char * const filename,
                        const UInt numline);
    
    // ===================================================
    // INTERFACE FUNCTIONS PROVIDED BY CFS++
    // ===================================================
    
    //! Evaluate expression
    static bool CFSEval( const StdVector<std::string> & args,
                            StdVector<std::string> & retVal );
    
  protected:
    
    //! Error Message of central parsing function
    static std::string errMsg_;
    
    //! Flag indicating if a scripting command is evaluated at the moment
    bool isEvaluating_;
    
    //! Map containing for each event the number of parameters to be passed
    std::map<EventType, UInt> eventNumParams_;
  };
  
} // end of namespace

#endif


