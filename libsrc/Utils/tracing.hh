#ifndef FILE_DEBUGGER_HH
#define FILE_DEBUGGER_HH

#ifdef TRACE

#include <iostream>


//! Macro for specifying depth of indentation in function trace log
#define TRACE_INDENT "  "


namespace CoupledField {

  //! Function tracing is written to this file stream
  extern std::ostream *trace;


  // ========================================================================
  // FCN_TRACE CLASS
  // ========================================================================

  //! Class for generating a doubly linked list of called functions

  //! This class performs the actual logging of function trace information
  //! to the trace stream. Each object stores a pointer to its predecessor
  //! (i.e. the function which called this one) and its successor (i.e. the
  //! function it called itself). The latter information is dynamic and
  //! changes, once a called subfunction has terminated.
  class FcnTrace{

  public:

    //! Constructor
    FcnTrace( char *name, int depth ){
      this->name_ = name;
      fcn_depth_ = depth;
    
    if (fcn_depth_<=TRACE){
	for (int i=0;i<fcn_depth_;i++)
	  (*trace) << TRACE_INDENT;
	(*trace) << "entering function " << name_ << std::endl;
      }
    }

    //! Default destructor

    //! The default destructor is responsible for issuing a "leaving
    //! function" message to the trace stream object.
    ~FcnTrace(){
      if (fcn_depth_<=TRACE){
	for (int i=0;i<fcn_depth_;i++)
	  (*trace) << TRACE_INDENT;
	(*trace) << "leaving function " << name_ << std::endl;
      }
      fcn_depth_=0;
      name_=NULL;
    }

    FcnTrace *caller_; //!< Link to FcnTrace object for predecessor
    FcnTrace *called_; //!< Link to FcnTrace object for successor

  private:
    int fcn_depth_;
    char *name_;

  };


  // ========================================================================
  // DEBUGGER CLASS
  // ========================================================================

  //! Central class object for fuction tracing.

  //! This static class is taking care of generating function trace
  //! information. If will only exist, if the TRACE macro is defined
  //! during compilation.
  class Debugger {

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
	foo_ = new FcnTrace( name, fcnTraceDepth_++ );
	foo_->caller_ = NULL;
	foo_->called_ = NULL;
      }
      else{
	foo_->called_ = new FcnTrace( name, fcnTraceDepth_++ );
	foo_->called_->caller_ = foo_;
	foo_ = foo_->called_;
      }
    }

    //! Handle termination of a function.

    //! This method performs all things needed for tracing when a function
    //! is left.
    static void LeaveFcn(){
      FcnTrace *tmp;
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
    static FcnTrace *foo_;
    static unsigned int fcnTraceDepth_;
  };


  // ========================================================================
  // FCNOBJ CLASS
  // ========================================================================

  //! When the TRACE macro is defined and the ENTER_FCN macro is called in
  //! a function, then a local object of this class will be instantiated.
  //! Its constructor and destructor trigger function tracing via the
  //! static Debugger class.
  class FcnObj {
  public:

    //! Constructor
    FcnObj(char *name){
      Debugger::EnterFcn(name);
    }

    //! Destructor
    ~FcnObj(){
      Debugger::LeaveFcn();
    }

  };


} // namespace

#endif // TRACE

#endif // FILE_DEBUGGER_HH
