#ifndef FILE_HEAT_NEUMANN_INT
#define FILE_HEAT_NEUMANN_INT

#include "Forms/linSurfForm.hh"

namespace CoupledField 
{
  
  //! Class implementing an inhomogeneous Neumann integrator for acoustic field
  class HeatNeumannInt : public LinearSurfForm {

  public:
    
    //! Standard constructor
    HeatNeumannInt( Double factor, Boolean isaxi );
    
    //! Destructor
    ~HeatNeumannInt();

    /// Calculation of vector of right hand side 
    void CalcElemVector( Matrix<Double>& ptCoord, Vector<Double> & elemVec );
    
  protected:

  };

} // end of namespace
#endif
