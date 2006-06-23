#ifndef FILE_LIN_PRESSURE_INT
#define FILE_LIN_PRESSURE_INT

#include "Forms/linSurfForm.hh"

namespace CoupledField {
  
  //! Class for calculating Pressure load vector on surface elements
  class PressureLinForm : public LinearSurfForm {
    
  public:
    
    //! Standard constructor
    PressureLinForm( Double aVal, bool isaxi );
    
    //! Destructor
    virtual ~PressureLinForm();
    
    //! Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & elemVec,
                         EntityIterator& ent );
    
  protected:
    
  private:
    
    /// factor of load
    Double multiplier_;
    
  };

}// end of namespace

#endif

  
