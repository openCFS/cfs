#ifndef FILE_DEBUGGER_HH
#define FILE_DEBUGGER_HH

#ifdef TRACE

#include <iostream>

// Depending on whether we use LAS or OLAS the trace stream pointer is in
// different namespaces
#ifdef USE_OLAS
namespace OutInfo {
  extern std::ostream * trace;
}
#else
namespace CoupledField {
  extern std::ostream * trace;
}
#endif

//! Stream to which trace information is sent
#define TRACESTREAM trace

//! Macro for specifying depth of indentation in function trace log
#define TRACE_INDENT "    "

namespace OutInfo {

  // For LAS import trace from namespace
#ifndef USE_OLAS
  using CoupledField::trace;
#endif


  // ========================================================================
  // FCNTRACELISTELEM CLASS
  // ========================================================================

  //! Class for generating a doubly linked list of called functions

  //! This class performs the actual logging of function trace information
  //! to the trace stream. Each object stores a pointer to its predecessor
  //! (i.e. the function which called this one) and its successor (i.e. the
  //! function it called itself). The latter information is dynamic and
  //! changes, once a called subfunction has terminated. It is used by the
  //! FcnTraceHandler class.
  class FcnTraceListElem{

  public:

    //! Constructor
    FcnTraceListElem( char *name, int depth ){
      this->name_ = name;
      fcnDepth_ = depth;
    
    if (fcnDepth_<=TRACE && 
	TRACESTREAM != NULL)
      {
	for (int i = 0; i < fcnDepth_; i++) {
	  (*TRACESTREAM) << TRACE_INDENT;
	}
	(*TRACESTREAM) << "entering function " << name_ << std::endl;
      }
    }
    
    //! Default destructor

    //! The default destructor is responsible for issuing a "leaving
    //! function" message to the trace stream object.
    ~FcnTraceListElem(){
      if (fcnDepth_<=TRACE &&
	  TRACESTREAM != NULL)
	{
	  for (int i = 0; i < fcnDepth_; i++ ) {
	    (*TRACESTREAM) << TRACE_INDENT;
	  }
	  (*TRACESTREAM) << "leaving function " << name_ << std::endl;
	}
      fcnDepth_ = 0;
      name_ = NULL;
    }
    
    FcnTraceListElem *caller_; //!< Link to FcnTrace object for predecessor
    FcnTraceListElem *called_; //!< Link to FcnTrace object for successor

  private:
    int fcnDepth_;
    char *name_;

  };


  // ========================================================================
  // FCNTRACEHANDLER CLASS
  // ========================================================================

  //! Central class object for fuction tracing.

  //! This static class is taking care of generating function trace
  //! information. If will only exist, if the TRACE macro is defined
  //! during compilation.
  class FcnTraceHandler {

  public:

    //! Init method
    static void Init(){
      fcnTraceDepth_ = 0;
      foo_ = NULL;
    }

    //! Handle startup of a function.

    //! This method performs all things needed for tracing at the startup
    //! of a method.
    static void EnterFcn(char *name){
      if ( foo_ == NULL ){
	foo_ = new FcnTraceListElem( name, fcnTraceDepth_++ );
	foo_->caller_ = NULL;
	foo_->called_ = NULL;
      }
      else{
	foo_->called_ = new FcnTraceListElem( name, fcnTraceDepth_++ );
	foo_->called_->caller_ = foo_;
	foo_ = foo_->called_;
      }
    }

    //! Handle termination of a function.

    //! This method performs all things needed for tracing when a function
    //! is left.
    static void LeaveFcn(){
      FcnTraceListElem *tmp;
      if (foo_){
	tmp = foo_->caller_;
	delete foo_;
	foo_ = tmp;
	fcnTraceDepth_--;
      }
    }

    //! This method writes a complete call trace
    static void Dump(){
      while(foo_) LeaveFcn();
    }

  private:
    static FcnTraceListElem *foo_;
    static unsigned int fcnTraceDepth_;
  };


  // ========================================================================
  // FCNTRACEOBJLOCAL CLASS
  // ========================================================================

  //! When the TRACE macro is defined and the ENTER_FCN macro is called in
  //! a function, then a local object of this class will be instantiated.
  //! Its constructor and destructor trigger function tracing via the
  //! static FcnTraceHandler class.
  class FcnTraceObjLocal {
  public:

    //! Constructor
    FcnTraceObjLocal(char *name){
      FcnTraceHandler::EnterFcn(name);
    }

    //! Destructor
    ~FcnTraceObjLocal(){
      FcnTraceHandler::LeaveFcn();
    }

  };


} // namespace

#endif // TRACE

#endif // FILE_DEBUGGER_HH
