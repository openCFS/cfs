// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "utils/profiler.hh"

#ifdef PROFILING

#include <iostream>
#include <iomanip>
#include <map>
#include <sys/time.h>
#include "utils/environment.hh"
//#include "graph/baseparallel.hh"

namespace OLAS {

  //! called after a function with profiling
  void Profiler::Update(char *name,int ops, double timediff) {
    profile_[name].setnops(ops);
    profile_[name].addtime(timediff);
  }

  void Profiler::WriteReport() {
		
    char *name;
    int ncalls;
    double nops,NOPS;
    double totaltime;
    double perf;
	
    (*cla) << std::endl;
    (*cla) << "###################      PROFILING REPORT      ################"
	   << std::endl;
    (*cla) << 	std::setiosflags(std::ios::left) 
	   << std::setw(24) << "function" 
	   << std::setw(8) 	<< "calls"
	   << std::setw(16) << "real time"
	   << std::setw(16) << "performance" << std::endl;
    (*cla) << "###############################################################"
	   << std::endl;
		
    for( ProfMap::iterator i = profile_.begin(); i != profile_.end(); i++ ) {
      name = i->first;
      ncalls = i->second.ncalls_;
      nops = (double)(i->second.nops_)*1e-6;

      NOPS=nops;

      totaltime = i->second.time_;
      perf = (double)NOPS;
      perf *= ncalls;
      perf /= totaltime;
	
      (*cla)<< std::setw(24) << name
	    << std::setw(8)  << ncalls
	    << std::setw(16) << totaltime
	    << std::setw(16) << perf
	    << std::endl;
			   	   
      // (*cla) << " ticks : " << totaltime << " ops in fcn: " << 
      // nops << std::endl;
    }
			
    (*cla) << "###############################################################"
	   << std::endl;

  }//WriteReport
		
  double Profiler::GetRealTime() {
    struct timeval tv;
    double t;
    gettimeofday(&tv,NULL);
    t = tv.tv_sec +  tv.tv_usec*1.0e-6;
    return t;
  }


  FcnProf::FcnProf( char* name, int nops ) {
    stime_ = Profiler::GetRealTime();
    name_=name;
    nops_=nops;
  }
		
  FcnProf::~FcnProf() {
	//// commented out the barrier on 6.Nov 2004 (Uwe Fabricius)
	// The profiling class FcnProf is also instanciated in member
	// functions of purely serial classes, e.g. Vector, CRS_Matrix.
	// These classes might be instantiated in a parallel program
	// only on one process or on a subset of the processes in
	// OLAS_MPI_COMM_WORLD, what would result in a dead lock.
	// Skipping the barrier influences the profiling results, but
	// as a first step it removes at least these dead locks.
//  OLAS_MPI_Barrier(OLAS_MPI_COMM_WORLD);

    etime_ = Profiler::GetRealTime();
    Profiler::Update(name_,nops_,etime_-stime_);
  }

  ProfMap Profiler::profile_;

} // namespace

#endif


