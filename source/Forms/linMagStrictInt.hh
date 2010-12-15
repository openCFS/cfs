// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_LIN_MAGMECH_COUPLING
#define FILE_LIN_MAGMECH_COUPLING


#include "Elements/basefe.hh"
#include "Forms/adbInt.hh"
#include "Materials/baseMaterial.hh"

namespace CoupledField {

//! Class describing the  \f$ ADB \f$ operator for magmech coupling 

//! The main objective of this class is to implement the pure virtual
//! methods of the ADBInt parent class for the case of a linear
//! magnetostrictive coupling operator 
class LinMagStrictInt : public ADBInt {


public:

  // =======================================================================
  // CONSTRUCTION and DESTRUCTION
  // =======================================================================

  //@{ \name Construction and destruction

  //! Default constructor
  LinMagStrictInt(BaseMaterial* matData, SubTensorType type); 

  //! Destructor
  ~LinMagStrictInt() {
  };

  //@}

  /*!
  Compute the matrix \f$ A \f$ of the \f$ ADB \f$ operator
  \param aMat    (output) computed matrix \f$ A \f$
  \param ip      (input)  number of integration point
  \param ptCoord (input)  matrix containing co-ordinates of all
  integration points 
  */
  void calcAMat( Matrix<Double> &aMat, UInt ip,
                 const Matrix<Double> &ptCoord );

  /*!    Compute the matrix \f$ B \f$ of the \f$ ADB \f$ operator
  \param bMat    (output) computed matrix \f$ B \f$
  \param ip      (input)  number of integration point
  \param ptCoord (input)  matrix containing co-ordinates of all
  integration points
  */
  void calcBMat( Matrix<Double> &bMat, UInt ip,
                 const Matrix<Double> &ptCoord );

  /*!   Compute the data-matrix \f$ D \f$
  The method computes the matrix \f$ D \f$  of the piezoelectric coupling
  properties of the element's material. In the 3D setting the latter
  is given by
  \f$[
  D = \left( e^T \right) \in \mathcal{R}^{6\times 3}.
  \f$]
  where \f$ e^T \f$ is the local tensor of piezoelectric coupling.*/
  virtual void calcDMat(Matrix<Double> &dMat, const Elem* elem);

  /*!   Query dimension of the matrix \f$ D \f$

  This method returns the dimensions of the data-matrix \f$ D \f$.
  In the case of 2D plane strain piezoelectric coupling we have
  \f$ \mbox{dim}D=3\times 2 \f$ */
  virtual void getDimD( UInt nRows, UInt nCols ) {
    nRows = matDimRow_;
    nCols = matDimCol_;
  };

  /*!   Query number of degrees of freedom for first physical quantity

  This method can be used to query the number of degrees of freedom at
  any node of a finite element for the physical quantity associated to
  the base functions whose derivatives form the matrix \f$ A \f$.
  In the case of 2D plane strain piezoelectric coupling the first
  physical quantity is the mechanical displacements.
  \return 2 */
  virtual UInt getNumDofsA() {
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

  /*!  Query material type for \f$ D \f$ tensor */
  MaterialType getDMaterialType() { return MAGNETOSTRICTION_TENSOR; }

private:

  //! Number of unknowns associated with differential operator A
  UInt numDofsA_;

  //! Number of unknowns associated with differential operator A
  UInt numDofsB_;

  //! Number of rows of material tensor
  UInt matDimRow_;

  //! Number of columns of material tensor
  UInt matDimCol_;

};

}

#endif




