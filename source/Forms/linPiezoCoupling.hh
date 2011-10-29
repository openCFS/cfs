// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_LIN_PIEZO_COUPLING
#define FILE_LIN_PIEZO_COUPLING


#include "Elements/basefe.hh"
//#include "Forms/adbInt.hh"
#include "Materials/baseMaterial.hh"

#include <Utils/ApproxData.hh>
#include <Forms/bdbInt.hh>
//#include <Forms/gradfieldop.hh>



namespace CoupledField {

  //! Class describing the  \f$ ADB \f$ operator for piezoelectric coupling 

  //! The main objective of this class is to implement the pure virtual
  //! methods of the ADBInt parent class for the case of a linear
  //! piezoelectric coupling operator 
  class linPiezoCoupling {


  public:

    // =======================================================================
    // CONSTRUCTION and DESTRUCTION
    // =======================================================================

    //@{ \name Construction and destruction

    //! Default constructor
    linPiezoCoupling(BaseMaterial* matData, SubTensorType type); 

    //! Destructor
    ~linPiezoCoupling() {
    };

    //@}
/*!
     Compute the matrix \f$ A \f$ of the \f$ ADB \f$ operator

     The method computes the matrix \f$ A \f$ of the \f$ ADB \f$ operator
     given in the case of 3D piezoelectric coupling by
     \f$[
     A = \left( \begin{array}{*{6}{c}}
     \displaystyle \frac{\partial N_u}{\partial x} &
     \displaystyle 0 &
     \displaystyle 0 &
     \displaystyle 0 &
     \displaystyle \frac{\partial N_u}{\partial z} &
     \displaystyle \frac{\partial N_u}{\partial y} \\[3ex]
     \displaystyle 0 &
     \displaystyle \frac{\partial N_u}{\partial y} &
     \displaystyle  0 &
     \displaystyle \frac{\partial N_u}{\partial z} &
     \displaystyle  0 &
     \displaystyle \frac{\partial N_u}{\partial x}\\[3ex]
     \displaystyle 0 &
     \displaystyle 0 &
     \displaystyle \frac{\partial N_u}{\partial z} &
     \displaystyle \frac{\partial N_u}{\partial y} &
     \displaystyle \frac{\partial N_u}{\partial x} &
     \displaystyle 0
     \end{array}\right)
     \f$]
     where \f$ N_u \f$ are the Finite Element ansatz functions for the
     mechanical displacement part. To be more precise the above matrix
     \f$ A \f$ is computed for the given integration point ip and every FE
     ansatz function belonging to a node of the element. These partial
     matrices are appended one after another to form the output matrix
     aMat in the following fashion
     \f$[
     \mbox{aMat} = \left(
     \begin{array}{c}
     \displaystyle A(p_1) \\[3ex]
     \displaystyle A(p_2) \\[3ex]
     \displaystyle \vdots \\[3ex]
     \displaystyle A(p_N)
     \end{array} \right)
     \in \mathcal{R}^{3N\times 6}
     \f$]
     where \f$ p_k, k=1,\ldots,N \f$ are the nodes of the element.
     \param aMat    (output) computed matrix \f$ A \f$
     \param ip      (input)  number of integration point
     \param ptCoord (input)  matrix containing co-ordinates of all
                             integration points */
    void calcAMat( Matrix<Double> &aMat, UInt ip,
                   const Matrix<Double> &ptCoord );

   /*!    Compute the matrix \f$ B \f$ of the \f$ ADB \f$ operator

      The method computes the matrix \f$ B \f$ of the \f$ ADB \f$ operator
      given in the case of 3D piezoelectric coupling by
      \f$[
      B = \left(\begin{array}{c}
      \displaystyle\frac{\partial N_\Phi}{\partial x} \\[3ex]
      \displaystyle\frac{\partial N_\Phi}{\partial y} \\[3ex]
      \displaystyle\frac{\partial N_\Phi}{\partial z}
      \end{array}\right)
      \f$]
      where \f$ N_\Phi \f$ are the Finite Element ansatz functions for the
      electric potential.
      To be more precise the above matrix \f$ B \f$ is computed for the given
      integration point ip and every FE ansatz function belonging to a node
      of the element.
      These partial matrices are appended one after another to form the
      output matrix bMat in the following fashion
      \f$[
      \mbox{bMat} = \left( B(p_1), B(p_2), \ldots, B(p_N) \right)
      \in \mathcal{R}^{3\times N}
      \f$]
      where \f$ p_k, k=1,\ldots,N \f$ are the nodes of the element.
      @see BaseForm::CalcBMat() */
    void CalcBMat(Matrix<Double> &bMat, UInt ip, const Matrix<Double> &ptCoord);

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
    MaterialType getDMaterialType() { return PIEZO_TENSOR; }



  private:

    UInt numDofsA_;
    UInt numDofsB_;

    UInt matDimRow_;
    UInt matDimCol_;

    SubTensorType sunTensorType_;

  };

}

#endif
