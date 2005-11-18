#ifndef FILE_ACOU_NEUMANN_INT
#define FILE_ACOU_NEUMANN_INT

#include "Forms/linSurfForm.hh"

namespace CoupledField 
{
  
  //! Class implementing an inhomogeneous Neumann integrator for acoustic field
  class AcouNeumannInt : public LinearSurfForm {

  public:
    
    //! Standard constructor
    AcouNeumannInt( Double factor, Boolean isaxi );
    
    //! Destructor
    ~AcouNeumannInt();

    /// Calculation of vector of right hand side 
    void CalcElemVector( Matrix<Double>& ptCoord, Vector<Double> & elemVec );
    
  protected:

  };

} // end of namespace
#endif
