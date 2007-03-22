// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef CONFIGURATION_HH_
#define CONFIGURATION_HH_

#include <olas_use_pardiso.hh>
#include <olas_use_lapack.hh>
#include <olas_use_ilupack.hh>
#include <olas_use_metis.hh>


namespace OLAS 
{
  /** This class is meant to be a static instance in global. It
   * is a frontend for the configuration, i.e. the compiler switches. */
  class Configuration
  {
     public:
        /** Pardiso is a state if the art direct solver. It is free for academic use
         * but needs to be registered per host, per user, per year. It is also useful
         * for the iterative precond/solver ilupack */
        bool HasPardiso(){
           #ifdef USE_PARDISO
              return true;
           #else 
              return false;
           #endif   
        };   
        
        /** Ilupack is a set of state of the art preconditioners and includes also
         * solvers. Is free for academic use. Can (should) make use of pardiso and mc64 */
        bool HasIlupack() {
           #ifdef USE_ILUPACK
              return true;
           #else 
              return false;
           #endif   
        };
 
        
        bool HasMetis() {
           #ifdef USE_METIS
              return true;
           #else 
              return false;
           #endif   
        };

        /** MC64 is a solver package which can optionally be used by ilupack */
        bool HasMC64() {
           #ifdef USE_MC64
              return true;
           #else 
              return false;
           #endif   
        };
        
  }; 

} // end of namespace 

#endif /*CONFIGURATION_HH_*/
