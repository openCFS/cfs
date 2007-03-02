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
                     bool isaxi );
    
    //! Destructor
    virtual ~PressureLinForm();
    
    //! Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & elemVec,
                         EntityIterator& ent );

    //! Calculation of vector of right hand side 
    void CalcElemVector( Vector<Complex> & elemVec,
                         EntityIterator& ent );
    
  protected:
    
  private:

    //! internal method, which performs the integration part
    void PrepareElemVec( Vector<Double>& elemVec, 
                         EntityIterator& ent );
    
    //! value string
    std::string value_;

    //! phase string
    std::string phase_;

  };

}// end of namespace

#endif

  
