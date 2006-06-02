#ifndef FILE_ACOU_NEUMANN_INT
#define FILE_ACOU_NEUMANN_INT

#include "Forms/linSurfForm.hh"

namespace CoupledField 
{
  
  //! Class implementing an inhomogeneous Neumann integrator for acoustic field
  class LinNeumannInt : public LinearSurfForm {

  public:
    
    //! Standard constructor
    LinNeumannInt( Double amplitude, MaterialType materialParam, Boolean isaxi );
    
    //! Destructor
    ~LinNeumannInt();

    /// Calculation of vector of right hand side 
    void CalcElemVector( Matrix<Double>& ptCoord, Vector<Double> & elemVec );
    
  protected:

    Double amplitude_;
    MaterialType materialParam_;

  };

} // end of namespace
#endif
