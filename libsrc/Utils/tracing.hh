#ifndef FILE_DEBUGGER_HH
#define FILE_DEBUGGER_HH

#ifdef TRACE

#include <iostream>

// namespace OutInfo{
//   extern std::ostream *trace;
//}

namespace CoupledField {

//! function tracing is written to this file stream
  extern std::ostream *trace;
//using OutInfo::trace;

//! this class is a list entry for function tracing
  class fcn_trace{
  public:
    fcn_trace(char *name,int depth){
      this->name_=name;
      fcn_depth_=depth;
    
    if (fcn_depth_<=TRACE){
	for (int i=0;i<fcn_depth_;i++)
	  (*trace) << "\t";
	(*trace) << "entering function " << name_ << std::endl;
      }//if
    }//fcn_trace
  
    ~fcn_trace(){
      if (fcn_depth_<=TRACE){
	for (int i=0;i<fcn_depth_;i++)
	  (*trace) << "\t";
	(*trace) << "leaving function " << name_ << std::endl;
      }//if
      fcn_depth_=0;
      name_=NULL;
    }//~
  public:
    int fcn_depth_;

    fcn_trace *caller_, *called_;

  protected:
    char *name_;
  };

//! this class is static. It controls function tracing
//! if TRACE is defined
  class Debugger {

  public:
    static void Init(){
      fcn_trace_depth_=0;
      foo_ = NULL;
    }

    //!enter a function
    static void EnterFcn(char *name){
      if (!foo_){
	foo_=new fcn_trace(name,fcn_trace_depth_++);
	foo_->caller_=NULL;
	foo_->called_=NULL;
      }
      else{
	foo_->called_=new fcn_trace(name,fcn_trace_depth_++);
	foo_->called_->caller_=foo_;
	foo_=foo_->called_;
      }
    }
    //!leave a function
    static void LeaveFcn(){
      fcn_trace *tmp;
      if (foo_){
	tmp=foo_->caller_;
	delete foo_;
	foo_=tmp;
	fcn_trace_depth_--;
      }//if
    }//leave

    //!write a complete call trace
    static void Dump(){
      while(foo_) LeaveFcn();
    }

  protected:
    static fcn_trace *foo_;
    static int fcn_trace_depth_;
  };


//! when an instance of this class is generated,
//! it tells the static Debugger that a function has been entered.
//! upon destruction, it tells the Debugger that a function has been left
  class FcnObj {
public:

  FcnObj(char *name){
  Debugger::EnterFcn(name);
  }
  
  ~FcnObj(){Debugger::LeaveFcn();}

  };


} // namespace

#endif // TRACE

#endif // FILE_DEBUGGER_HH
