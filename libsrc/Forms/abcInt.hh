#ifndef File_ABC_MECH_INT
#define File_ABC_MECH_INT

#include "baseForm.hh"

namespace CoupledField {

  //! Class for calculation surface integrals of absorbing boundary
  //! elements
  class AbsorbingBCsInt : public SurfForm {

  public:
    
    //! Default constructor
    //! \param isAxi flag indicating axi-symmetrix geometry
    AbsorbingBCsInt( Boolean isAxi );

    //! Default destructor
    virtual ~AbsorbingBCsInt();

    //! Set an additional factor for multiplying the form
    void SetFactor(Double factor) {
      factor_ = factor;
    }

    //! Calculate elementwise matrix

    //! This method calculates the element matrix for absorbing boundary
    //! conditions of first order. 
    void CalcElementMatrix(Matrix<Double>& ptCoord, 
                           Matrix<Double> & StiffMat);

  private:
    
    //! Additional multiplicative factor
    Double factor_;
    
  };



} // end of namespace

#endif
