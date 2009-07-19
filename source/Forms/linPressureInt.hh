// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_LIN_PRESSURE_INT
#define FILE_LIN_PRESSURE_INT

#include "Forms/linSurfForm.hh"

namespace CoupledField {
  
  //! Class for calculating Pressure load vector on surface elements
  class PressureLinForm : public LinearSurfForm {
    
  public:
    
    //! Standard constructor
    PressureLinForm( const std::string& value, 
                     const std::string& phase,
                     const std::string& subType,
                     bool isaxi );
    
    //! Destructor
    virtual ~PressureLinForm();
    
    //! Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & elemVec,
                         EntityIterator& ent );

    //! Calculation of vector of right hand side 
    void CalcElemVector( Vector<Complex> & elemVec,
                         EntityIterator& ent );

    void SetValue(const std::string& value) { value_ = value; }
    
    //! Calculapressure factor, this is called from shapeopt
    double GetPressureFactor(const SurfElem* elem);

  protected:
    
  private:

    //! internal method, which performs the integration part
    void PrepareElemVec( Vector<Double>& elemVec, 
                         EntityIterator& ent );
    
    //! value string
    std::string value_;

    //! phase string
    std::string phase_;
    
    //! subtype of integrator
    std::string subType_;


  };

}// end of namespace

#endif

  
