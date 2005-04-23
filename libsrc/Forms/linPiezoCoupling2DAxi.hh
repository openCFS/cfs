#ifndef FILE_LIN_PIEZO_COUPLING_2D_AXI
#define FILE_LIN_PIEZO_COUPLING_2D_AXI


#include "Elements/basefe.hh"
#include "Forms/adbInt.hh"
#include "DataInOut/MaterialData.hh"


namespace CoupledField {

  //! Class describing ADB operator for piezoelectric axi-symmetric coupling

  //! The main objective of this class is to implement the pure virtual
  //! methods of the ADBInt parent class for the case of a linear
  //! piezoelectric coupling operator for a 2D axi-symmetric simulation.
  //! \note The class currently assumes that the computational 2D area lies
  //! in the yz-plane with the z-axis being the axis of symmetry.
  class linPiezoCoupling2DAxi : public ADBInt {


  public:

    // =======================================================================
    // CONSTRUCTION and DESTRUCTION
    // =======================================================================

    //@{ \name Construction and destruction

    //! Default constructor
    linPiezoCoupling2DAxi() {
      ENTER_FCN( "linPiezoCoupling2DAxi::linPiezoCoupling2DAxi" );
      isaxi_ = TRUE;
    }

    //! Destructor
    ~linPiezoCoupling2DAxi() {
      ENTER_FCN( "linPiezoCoupling2DAxi::~linPiezoCoupling2DAxi" );
    };

    //@}


    //! Compute the matrix \f$A\f$ of the \f$ADB\f$ operator

    //! The method computes the matrix \f$A\f$ of the \f$ADB\f$ operator
    //! which, in the case of 2D axi-symmetric piezoelectric coupling,
    //! is given by
    //! \f[
    //! A = \left( \begin{array}{*{4}{c}}
    //! \displaystyle \frac{\partial N_u}{\partial r} &
    //! \displaystyle  0 &
    //! \displaystyle \frac{\partial N_u}{\partial z} &
    //! \displaystyle \frac{1}{r} \\[3ex]
    //! \displaystyle 0 &
    //! \displaystyle \frac{\partial N_u}{\partial z} &
    //! \displaystyle \frac{\partial N_u}{\partial r} &
    //! \displaystyle 0
    //! \end{array}\right)
    //! \f]
    //! where \f$N_u\f$ are the Finite Element ansatz functions for the
    //! mechanical displacement part. To be more precise the above matrix
    //! \f$A\f$ is computed for the given integration point ip and every FE
    //! ansatz function belonging to a node of the element. These partial
    //! matrices are appended one after another to form the output matrix
    //! aMat in the following fashion
    //! \f[
    //! \mbox{aMat} = \left(
    //! \begin{array}{c}
    //! \displaystyle A(p_1) \\[3ex]
    //! \displaystyle A(p_2) \\[3ex]
    //! \displaystyle \vdots \\[3ex]
    //! \displaystyle A(p_N)
    //! \end{array} \right)
    //! \in \mathcal{R}^{2N\times 4}
    //! \f]
    //! where \f$p_k, k=1,\ldots,N\f$ are the nodes of the element.
    //! \param aMat    (output) computed matrix \f$A\f$
    //! \param ip      (input)  number of integration point
    //! \param ptCoord (input)  matrix containing co-ordinates of all
    //!                         integration points
    void calcAMat( Matrix<Double> &aMat, Integer ip,
                   const Matrix<Double> &ptCoord );

    //! Compute the matrix \f$B\f$ of the \f$ADB\f$ operator

    //! The method computes the matrix \f$B\f$ of the \f$ADB\f$ operator
    //! which, in the case of 2D axi-symmetric piezoelectric coupling,
    //! is given  by
    //! \f[
    //! B = \left(\begin{array}{c}
    //! \displaystyle\frac{\partial N_\Phi}{\partial r} \\[3ex]
    //! \displaystyle\frac{\partial N_\Phi}{\partial z}
    //! \end{array}\right)
    //! \f]
    //! where \f$N_\Phi\f$ are the Finite Element ansatz functions for the
    //! electric potential.
    //! To be more precise the above matrix \f$B\f$ is computed for the given
    //! integration point ip and every FE ansatz function belonging to a node
    //! of the element.
    //! These partial matrices are appended one after another to form the
    //! output matrix bMat in the following fashion
    //! \f[
    //! \mbox{bMat} = \left( B(p_1), B(p_2), \ldots, B(p_N) \right)
    //! \in \mathcal{R}^{2\times N}
    //! \f]
    //! where \f$p_k, k=1,\ldots,N\f$ are the nodes of the element.
    //! \param bMat    (output) computed matrix \f$B\f$
    //! \param ip      (input)  number of integration point
    //! \param ptCoord (input)  matrix containing co-ordinates of all
    //!                         integration points
    void calcBMat( Matrix<Double> &bMat, Integer ip,
                   const Matrix<Double> &ptCoord );

    //! Compute the data-matrix \f$D\f$

    //! The method computes the matrix D of the piezoelectric coupling
    //! properties of the element's material. In the 2D setting we assume that
    //! the computational 2D area lies in the yz-plane with the z-axis being
    //! the axis of symmetry. Thus D is given by
    //! \f[
    //! D = \left(
    //! \begin{array}{cc}
    //! \displaystyle d_{28} & \displaystyle d_{29} \\
    //! \displaystyle d_{38} & \displaystyle d_{39} \\
    //! \displaystyle d_{48} & \displaystyle d_{49} \\
    //! \displaystyle d_{18} & \displaystyle d_{19} \\
    //! \end{array} \right)
    //! \f]
    //! where \f$\displaystyle d=(d_{ij})\f$ is the full local 9 x 9 tensor
    //! of piezoelectric coupling in 3D. If the full tensor is given by
    //! \f[
    //! D_{\mbox{full}} = \left(
    //! \begin{array}{ccc}
    //!            0         &           0          & \displaystyle e_{31} \\
    //!            0         &           0          & \displaystyle e_{32} \\
    //!            0         &           0          & \displaystyle e_{33} \\
    //!            0         & \displaystyle e_{34} &           0          \\
    //! \displaystyle e_{35} &           0          &           0          \\
    //!            0         &           0          &           0          \\
    //! \end{array} \right)
    //! \qquad\longrightarrow\qquad
    //! D = \left(
    //! \begin{array}{cc}
    //!             0        & \displaystyle e_{32} \\
    //!             0        & \displaystyle e_{33} \\
    //! \displaystyle e_{34} &           0          \\
    //!             0        & \displaystyle e_{31} \\
    //! \end{array} \right)
    //! \enspace.
    //! \f]
    void calcDMat( Matrix<Double> &dMat );

    //! Query dimension of the matrix \f$D\f$

    //! This method returns the dimensions of the data-matrix \f$D\f$.
    //! In the case of 2D axi-symmetric piezoelectric coupling we have
    //! \f$\mbox{dim}D=4\times 2\f$
    void getDimD( Integer nRows, Integer nCols ) {
      ENTER_IFCN( "linPiezoCoupling2DAxi::getDimD" );
      nRows = 4;
      nCols = 2;
    };

    //! Query number of degrees of freedom for first physical quantity

    //! This method can be used to query the number of degrees of freedom at
    //! any node of a finite element for the physical quantity associated to
    //! the base functions whose derivatives form the matrix \f$A\f$.
    //! In the case of 2D axi-symmetric piezoelectric coupling the first
    //! physical quantity is the mechanical displacements.
    //! \return 2
    Integer getNumDofsA() {
      ENTER_IFCN( "linPiezoCoupling2DAxi::getNumDofsA" );
      return 2;
    }

    //! Query number of degrees of freedom for second physical quantity

    //! This method can be used to query the number of degrees of freedom at
    //! any node of a finite element for the physical quantity associated to
    //! the base functions whose derivatives form the matrix \f$B\f$.
    //! In the case of 2D axi-symmetric piezoelectric coupling the second
    //! physical quantity is the electric potential.
    //! \return 1
    Integer getNumDofsB() {
      ENTER_IFCN( "linPiezoCoupling2DAxi::getNumDofsB" );
      return 1;
    }

  };

}

#endif
