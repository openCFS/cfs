// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_PROFILER_HH
#define OLAS_PROFILER_HH

// Include build type options header containing the #defines
#include <def_build_type_options.hh>

#ifdef PROFILING

#include <map>

namespace CoupledField {


  //! this class is a map entry for profiling info
  class fcn_prof{

  public:
    //!
    fcn_prof(){
      nops_   = 0;
      ncalls_ = 0;
      time_   = 0;
    }
     
    //! set number of operations in the function
    void setnops( long ops ) {
      nops_ = ops;
    }
    
    //! add a new statistical value for the execution time
    void addtime( double time ) {
      time_ += time;
      ncalls_++;
    }
    
    //! how many operations are performed in the function?
    long nops_;

    //! how often has the function been called?
    long ncalls_;

    //! what is the total (real) time spent in the function (in microsecs)
    double time_;
  };

  //! define a type for storing function profiling info
  typedef std::map<char*,fcn_prof> ProfMap;

  //! this class provides simple profiling functionality

  /*! the profiler provides a static function which can be called
    at the end of functions implementing "number crunching" algorithms
    to obtain the profiling report, another static function is available
	  
    \note the profiling utility is only active if the macro PROFILING
    is defined, i.e by passing a -DPROFILING flag to the compiler
  */
  class Profiler {

  public:

    friend class FcnProf;

    //! called after a function with profiling
    static void Update(char *name,int ops, double timediff);

    //! writes a profiling report to the stream *cla
    static void WriteReport();
		
    //! get current real time in microseconds
    static double GetRealTime();
		
  protected:
    //! store profiling data for different functions  
    static ProfMap profile_;

  };


  //! instances of this class are created at the beginning
  //! of functions which are to be profiled. At the end 
  //! of the function, the real time elapsed inside the 
  //! function is computed and passed to the Profiler
  class FcnProf {

  public:
    FcnProf(char* name, int nops);
		
    ~FcnProf();
		
  private: 

    char *name_;	
    long nops_;
    double stime_,etime_;
  };

} // namespace

#endif // PROFILE

#endif // OLAS_PROFILER_HH

