#ifndef TRACING_HH
#define TRACING_HH

#include <def_build_type_options.hh>

#ifdef TRACE

#include <iostream>

// The trace stream pointer is also used by OLAS and therefore belongs to the
// joint CFS++ / OLAS namespace for IO-objects
namespace OutInfo {
  extern std::ostream *trace;
}

//! Stream to which trace information is sent
#define TRACESTREAM trace

//! Macro for specifying depth of indentation in function trace log
#define TRACE_INDENT "    "


namespace OutInfo {


  //! Private unsigned integer data type for function tracing
  typedef unsigned long TraceInt;


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
  class FcnTraceListElem {

  public:

    //! Constructor
    FcnTraceListElem( char *name, TraceInt depth, TraceInt limit ) {

      name_     = name;
      fcnDepth_ = depth;
      limit_    = limit;
    
      if ( fcnDepth_ <= limit_ && TRACESTREAM != NULL ) {
        for ( TraceInt i = 0; i < fcnDepth_; i++ ) {
          (*TRACESTREAM) << TRACE_INDENT;
        }
        (*TRACESTREAM) << "entering function " << name_ << std::endl;
      }
    }

    //! Destructor

    //! The destructor is responsible for issuing a "leaving
    //! function" message to the trace stream object.
    ~FcnTraceListElem() {
      if ( fcnDepth_ <= limit_ && TRACESTREAM != NULL ) {
        for ( TraceInt i = 0; i < fcnDepth_; i++ ) {
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
    TraceInt fcnDepth_;
    TraceInt limit_;
    char *name_;

  };


  // ========================================================================
  // FCNTRACEHANDLER CLASS
  // ========================================================================

  //! Central class object for function tracing.

  //! This static class is taking care of generating function trace
  //! information. If will only exist, if the TRACE macro is defined
  //! during compilation.
  class FcnTraceHandler {

  public:

    //! Init method
    static void Init() {
      fcnTraceDepth_ = 0;
      foo_ = NULL;
    }

    //! Handle startup of a function.

    //! This method performs all things needed for tracing at the startup
    //! of a method.
    static void EnterFcn( char *name ) {
      if ( fcnTraceDepth_ < fcnTraceDepthLimit_ ) {
        if ( foo_ == NULL ){
          foo_ = new FcnTraceListElem( name, fcnTraceDepth_++,
                                       fcnTraceDepthLimit_ );
          foo_->caller_ = NULL;
          foo_->called_ = NULL;
        }
        else{
          foo_->called_ = new FcnTraceListElem( name, fcnTraceDepth_++,
                                                fcnTraceDepthLimit_ );
          foo_->called_->caller_ = foo_;
          foo_ = foo_->called_;
        }
      }
    }

    //! Handle termination of a function.

    //! This method performs all things needed for tracing when a function
    //! is left.
    static void LeaveFcn() {
      FcnTraceListElem *tmp;
      if ( foo_ ) {
        tmp = foo_->caller_;
        delete foo_;
        foo_ = tmp;
        fcnTraceDepth_--;
      }
    }

    //! This method writes a complete call trace
    static void Dump() {
      while( foo_ ) LeaveFcn();
    }

    //! Set maximum tracing depths
    static void SetMaxTraceDepth( unsigned int limit ) {
      fcnTraceDepthLimit_ = limit;
    }

  private:
    static FcnTraceListElem *foo_;
    static TraceInt fcnTraceDepth_;
    static TraceInt fcnTraceDepthLimit_;
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
    FcnTraceObjLocal( char *name ) {
      FcnTraceHandler::EnterFcn( name );
    }

    //! Destructor
    ~FcnTraceObjLocal() {
      FcnTraceHandler::LeaveFcn();
    }

  };


}

#endif

#endif
