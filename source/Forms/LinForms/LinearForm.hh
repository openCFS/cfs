// =====================================================================================
// 
//       Filename:  LinearForm.hh
// 
//    Description:  This is the base Class for all Right hand side Integrators. A rhs
//                  integrator is defined via the structure of the integrand. These have 
//                  the most general form
//                  \int_K {\cal B} [d] \vec{U} \ \text{d} K
//                  Here, U can be a given quantity like the Lighthill sources for
//                  aeroacoustics or the solution of the PDE e.g. for postprocessing
//                  puposes. Right now we propose three derived classes
//                  1. BUInt.hh - For integrants without material dependency
//                  2. BDUInt.hh - For integrands with material dependency
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

#include "Domain/CoefFunction/CoefFunction.hh"
#include "FeBasis/FeSpace.hh"
#include "Domain/Domain.hh"


namespace CoupledField{

  // forward class declaration
  class FeSpace;

  class LinearForm{
    public:
      LinearForm(bool coordUpdate = false ){
        coordUpdate_ = coordUpdate;
        isSolDependent_ = false;
        extractReal_ = false;
        harm_ = 0;
      }

      //! Copy constructor
      LinearForm(const LinearForm& right ){
        this->coordUpdate_ = right.coordUpdate_;
        this->isSolDependent_ = right.isSolDependent_;
        this->name_ = right.name_;
        this->harm_ = right.harm_;
        this->ptFeSpace_ = right.ptFeSpace_;
        this->intScheme_ = right.intScheme_;
        this->extractReal_ = right.extractReal_;
        this->normalizeToVol_ = false;
      }

      //! Make deep copy of object pointer
      virtual LinearForm* Clone()=0;

      virtual ~LinearForm(){

      }

      virtual void CalcElemVector(Vector<Complex>& elemVec,EntityIterator& ent){
        EXCEPTION("RhsIntegrator::CalcElementVector called this may not happen!");
      }

      virtual void CalcElemVector(Vector<Double>& elemVec,EntityIterator& ent){
        EXCEPTION("RhsIntegrator::CalcElementVector called this may not happen!");
      }

      void SetName(const std::string& name ){
         name_ = name;
      }
      
      void SetHarm(const UInt& harm ){
         harm_ = harm;
      }

      UInt GetHarm() const {
         return harm_;
      }

      const std::string& GetName() const {
        return name_;
      }

      /** some debug output for logging */
      virtual const std::string ToString() { return GetName(); }

      virtual bool IsComplex() const =0;

      virtual void SetFeSpace(shared_ptr<FeSpace> feSpace ){
        ptFeSpace_ = feSpace;
        intScheme_ = ptFeSpace_->GetIntScheme();
      }
      
      //! Return if element matrix is solution dependend
      virtual bool IsSolDependent() { return isSolDependent_;}
      
      //! Return if inearform uses updated Lagrangian formulation
      bool IsCoordUpdate() { return coordUpdate_; }
      
      /** set coordUpdate, required by shape opt */
      void SetCoordUpdate(bool val) { coordUpdate_ = val; }

      //! Set explicit that linearform is solution dependent
      virtual void SetSolDependent() {isSolDependent_ = true;}

      bool IsExtractReal() const {
    	  return extractReal_;
      }

      virtual void SetScalVec(Vector<Double>& vec ) {
        scalVec_ = vec; 
      }

      virtual void SetNormalizeToVol() {
        normalizeToVol_ = true;
      }

    protected:
      //! name of linearform
      std::string name_;

      //! for a multiharmonic analysis, we need to have the information, which
      //! harmonic the rhs entry belongs to
      UInt harm_;

      //! flag for use of updated Lagrangian formulation
      bool coordUpdate_; 
      
      //! flag if linearform is solution dependent
      bool isSolDependent_;

      //! normalize RHS to volume
      bool normalizeToVol_;

      //! pointer to finite element space 1
      shared_ptr<FeSpace> ptFeSpace_;
      
      //! point to integration scheme
      shared_ptr<IntScheme> intScheme_;

      //! if true, only the real part of a complex quantity is considered,
      //! needed e.g. when taking the result of a complex analysis (energy) as
      //! input of a static simulation (source for example HeatPDE)
      bool extractReal_;

      //! set scaling for coefVec
      Vector<Double> scalVec_;
  };
}

#endif
