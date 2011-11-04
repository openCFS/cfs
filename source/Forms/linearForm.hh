// =====================================================================================
// 
//       Filename:  rhsIntegrator.hh
// 
//    Description:  This is the base Class for all Right hand side Integrators. A rhs
//                  integrator is defined via the structure of the integrand. These have 
//                  the most general form
//                  \int_K {\cal B} [d] \vec{U} \ \text{d} K
//                  Here, U can be a given quantity like the Lighthill sources for
//                  aeroacoustics or the solution of the PDE e.g. for postprocessing
//                  puposes. Right now we propose three derived classes
//                  1. rhsBUInt.hh - For integrants without material dependency
//                  2. rhsBDUInt.hh - For integrands with material dependency
//                  3. SingeEntryInt.hh - For the direct specificatin of the right hand
//                                        side without any kind of integration
//                  
// 
//        Version:  1.0
//        Created:  11/02/2011 10:02:29 AM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#ifndef FILE_LINEARFORM_
#define FILE_LINEARFORM_

namespace CoupledField{

  // forward class declaration
  class FeSpace;

  class LinearForm{
    public:
      LinearForm(){

      }

      ~LinearForm(){

      }

      virtual void CalcElemVector(Vector<Complex> elemVec,EntityIterator& ent){
        EXCEPTION("RhsIntegrator::CalcElementVector called this may not happen!");
      }

      virtual void CalcElemVector(Vector<Double> elemVec,EntityIterator& ent){
        EXCEPTION("RhsIntegrator::CalcElementVector called this may not happen!");
      }

      std::string GetName(){
        return name_;
      }

      virtual bool IsComplex()=0;

      virtual void SetFeSpace(shared_ptr<FeSpace> feSpace ){
        ptFeSpace_ = feSpace;
      }

    protected:
      //! name of linearform
      std::string name_;

      //! pointer to finite element space 1
      shared_ptr<FeSpace> ptFeSpace_;
  };
}

#endif
