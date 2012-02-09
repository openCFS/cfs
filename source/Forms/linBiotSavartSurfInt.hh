// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_LIN_BIOT_SAVART_SURF_INT
#define FILE_LIN_BIOT_SAVART_SURF_INT

#include <string>

#include "Forms/linSurfForm.hh"
#include "General/defs.hh"

namespace CoupledField {
  

  //! Forward class declarations
  class BiotSavart;
class EntityIterator;
template <class TYPE> class Vector;
  
  //! Calculates the Neumann-like RHS part of the Biot-Savart exciation
  
  //! This class calculates the 
  class LinBiotSavartSurfInt : public LinearSurfForm {
    
  public:
    
    //! Standard constructor
    LinBiotSavartSurfInt( shared_ptr<BiotSavart> bisa, 
                          bool isaxi );
    
    //! Destructor
    virtual ~LinBiotSavartSurfInt();
    
    //! Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & elemVec,
                         EntityIterator& ent );

    //! Set additional multiplicatove factor
    void SetFactor(const std::string& value);
    
  private:

    //! Pointer to Biot-Savart object
    shared_ptr<BiotSavart> biotSavart_;
    
  };

}// end of namespace

#endif

  
