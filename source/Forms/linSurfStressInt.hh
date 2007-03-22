// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_LIN_SURF_STRESS_INT
#define FILE_LIN_SURF_STRESS_INT

#include "Forms/linSurfForm.hh"

namespace CoupledField {
  
  //! Class for calculating surface stress (traction) load vectors in 3D
  class SurfStress3DLinForm : public LinearForm {
    
  public:
    
    //! Standard constructor
    SurfStress3DLinForm( Vector<Double> stressVec, 
                         shared_ptr<SurfElemList> surfElems );
    
    //! Destructor
    virtual ~SurfStress3DLinForm();
    
    //! Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & elemVec,
                         EntityIterator& ent );
    
  private:
    
    //! Global stress vector in Voigt notation
    Vector<Double> stress_;

    //! List which maps the surface elements to volume ones
    std::map<const Elem*,const SurfElem*> surfElems_;
    
  };

}// end of namespace

#endif
