/*
 * CoefFunctionRotation.hh
 *
 *  Created on: 02.12.2025
 *      Author: Michael Hofmarcher and Mark Socaciu Lendvai
 */

#ifndef COEFFUNCTIONROTATION
#define COEFFUNCTIONROTATION

#include <memory>
#include "CoefFunction.hh"

namespace CoupledField {

// ============================================================================
//  Coef Function Rotation
// ============================================================================
//! Rotates a material tensor coefficient function using two vector fields.
//!
//! The rotation procedure works as follows:
//!  1. The material tensor is assumed to be defined in the x-y-z coordinate
//!     system (i.e. the 1-2-3 axes in the mat-XML correspond to x-y-z).
//!  2. A rotation matrix is computed as the shortest rotation that
//!     maps \c refVec onto \c targetVec.
//!  3. The rotation matrix is applied to the material tensor via
//!     \f$ T' = R \, T \, R^\top \f$.
//!
//! \note This procedure is only sufficient for transversely isotropic material
//!       tensors where a single principal direction fully determines the tensor.
//!
//! Both 2D and 3D tensors are supported. The dimension of the vector fields
//! must match the tensor dimension (2D tensor → 2D vectors, 3D tensor → 3D
//! vectors).

  class CoefFunctionRotation : public CoefFunction , public enable_shared_from_this<CoefFunctionRotation>
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

    //! Material tensor coefficient function (defined in x-y-z / 1-2-3 system)
    PtrCoefFct matCoeff_;

    //! Reference vector field (principal direction as given in the mat-XML)
    PtrCoefFct vec1Coeff_;

    //! Target vector field (desired orientation in the simulation domain)
    PtrCoefFct vec2Coeff_;

  };

}


#endif /* COEFFUNCTIONROTATION */
