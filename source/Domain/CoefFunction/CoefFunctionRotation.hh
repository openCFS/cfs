/*
 * CoefFunctionRotation.hh
 *
 *  Created on: 02.12.2025
 *      Author: Michael Hofmarcher and Mark Socaciu Lendvai
 */

#ifndef COEFFUNCTIONROTATION
#define COEFFUNCTIONROTATION

#include <boost/shared_ptr.hpp>
#include "CoefFunction.hh"

namespace CoupledField {
// forward class declaration
// ============================================================================
//  Coef Function Up Coth
// ============================================================================
//! Coefficient function for alligning tensor fields with a vector field

  class CoefFunctionRotation : public CoefFunction , public boost::enable_shared_from_this<CoefFunctionRotation>
  {

  public:

    //! Constructor
    CoefFunctionRotation(PtrCoefFct matCoeff, PtrCoefFct vec1Coeff, PtrCoefFct vec2Coeff);

    //! Destructor
    virtual ~CoefFunctionRotation(){ }

    virtual string GetName() const override
    { 
      return "CoefFunctionRotation";
    }

    //! Return rotated tensor at integration point
    void GetTensor(Matrix<Double>& tensor, const LocPointMapped& lpm ) override;

    //! \copydoc CoefFunction::GetTensorSize
    virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const override{
      matCoeff_->GetTensorSize(numRows, numCols);
    }

  protected:

    // material tensor coefficient function
    PtrCoefFct matCoeff_;

    // reference vector field for rotation
    PtrCoefFct vec1Coeff_;

    // target vector field for rotation
    PtrCoefFct vec2Coeff_;

  };

}


#endif /* COEFFUNCTIONROTATION */
