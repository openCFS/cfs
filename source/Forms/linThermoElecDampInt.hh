#ifndef FILE_LINTHERMOELECDAMPINT
#define FILE_LINTHERMOELECDAMPINT

#include <Materials/baseMaterial.hh>
#include <General/environment.hh>

#include "Elements/basefe.hh"
#include "PDE/SinglePDE.hh"

#include "baseForm.hh"


namespace CoupledField
{

  /// Class for calculation  element stiffness matrix of a laplace operator
  class LinThermoElecDampInt : public BaseForm
  {
  public:

    /// Constructor
    LinThermoElecDampInt(BaseMaterial* matData, SubTensorType type);

    ///
    ~LinThermoElecDampInt();

    //void setNegate(bool negate){flagNegate_ = negate;}

    //!  Compute matrix \f$ A \f$ at given integration point.
    void calcAMat( Matrix<Double> &aMat, UInt ip,
                   const Matrix<Double> &ptCoord );

    //!  Compute matrix \f$ D \f$ at given integration point.
    void calcDMat( Matrix<Double> &dMat );

    //!  Compute matrix \f$ B \f$ at given integration point.
    void CalcBMat( Matrix<Double> &bMat, UInt ip,
                   const Matrix<Double> &ptCoord );

    //! Compute element matrix associated to damping integral form
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1,
                            EntityIterator& ent2 );
    // 	QUERIES


    /*!   Query dimension of the matrix \f$ D \f$

    This method returns the dimensions of the data-matrix \f$ D \f$.
    Which corresponds to pyroelectric constants.
    In the case of 3D, a pyroelectric constants matrix has
    \f$ \mbox{dim}D=3\times 1 \f$
    In the case of 2D, a pyroelectric constants matrix has
    \f$ \mbox{dim}D=2\times 1 \f$ */
    void getDimD( UInt nRows, UInt nCols ) {
      //nRows = matDimRow_;
      //nCols = matDimCol_;
    };

    /*!   Query number of degrees of freedom for first physical quantity

    This method can be used to query the number of degrees of freedom at
    any node of a finite element for the physical quantity associated to
    the base functions whose derivatives form the matrix \f$ A \f$.
    In the case of 2D plane strain piezoelectric coupling the first
    physical quantity is the mechanical displacements.
    \return 2 */
    UInt getNumDofsA() {
      return numDofsA_;
    }

    /*!   Query number of degrees of freedom for second physical quantity

    This method can be used to query the number of degrees of freedom at
    any node of a finite element for the physical quantity associated to
    the base functions whose derivatives form the matrix \f$ B \f$.
    In the case of 2D plain strain piezoelectric coupling the second
    physical quantity is the electric potential.
    \return 1 */
    UInt getNumDofsB() {
      return numDofsB_;
    }

    //! Query material type for \f$ D \f$ tensor.
    MaterialType getDMaterialType();

  private:

    UInt numDofsA_;
    UInt numDofsB_;

    // dimension of the D matrix
    UInt matDimRow_;
    UInt matDimCol_;

		// To keep the last heat solution
		Vector< Double > teta_;

    PtrParamNode pn_;

    //flag to set -K
    //bool flagNegate_;


  };


}

#endif // FILE_LINTHERMOELECDAMPINT
