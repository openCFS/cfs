#ifndef FILE_LIN_PIEZO_COUPLING_2D_PLANESTRAIN
#define FILE_LIN_PIEZO_COUPLING_2D_PLANESTRAIN


#include "Elements/basefe.hh"
#include "Forms/adbInt.hh"
#include "DataInOut/MaterialData.hh"


namespace CoupledField {

  //! Class describing ADB operator for piezoelectric plane strain coupling

  //! The main objective of this class is to implement the pure virtual
  //! methods of the ADBInt parent class for the case of a linear
  //! piezoelectric coupling operator for a 2D plane strain simulation.
  //! \note The class currently assumes that the computational 2D area lies
  //! in the yz-plane with the x-axis being the direction of large extension
  //! that can be neglected.
  class linPiezoCoupling2DPlaneStrain : public ADBInt {


  public:

    // =======================================================================
    // CONSTRUCTION and DESTRUCTION
    // =======================================================================

    //@{ \name Construction and destruction

    //! Default constructor
    linPiezoCoupling2DPlaneStrain() {
      ENTER_FCN( "linPiezoCoupling2DPlaneStrain::"
                 "linPiezoCoupling2DPlaneStrain" );
    }

    //! Destructor
    ~linPiezoCoupling2DPlaneStrain() {
      ENTER_FCN( "linPiezoCoupling2DPlaneStrain::"
                 "~linPiezoCoupling2DPlaneStrain" );
    };

    //@}


    //! Compute the matrix \f$A\f$ of the \f$ADB\f$ operator

    //! The method computes the matrix \f$A\f$ of the \f$ADB\f$ operator
    //! which, in the case of 2D PlaneStrain-symmetric piezoelectric coupling,
    //! is given by
    //! \f[
    //! A = \left( \begin{array}{*{3}{c}}
    //! \displaystyle \frac{\partial N_u}{\partial y} &
    //! \displaystyle 0                               &
    //! \displaystyle \frac{\partial N_u}{\partial z} \\[3ex]
    //! \displaystyle  0 &
    //! \displaystyle \frac{\partial N_u}{\partial z} &
    //! \displaystyle \frac{\partial N_u}{\partial y}
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
    void calcAMat( Matrix<Double> &aMat, UInt ip,
                   const Matrix<Double> &ptCoord ) {

      ENTER_FCN( "linPiezoCoupling2DPlaneStrain::calcAMat" );

      // Obtain info on problem sizes
      const UInt numNodes = ptelem->GetNumNodes();
      const UInt nDofMech = 2;
      const UInt nRowsD   = 3;

      // Set correct size of matrix A and initialise with zeros
      aMat.Resize( numNodes * nDofMech, nRowsD );

      // Get local shape functions and their derivatives with respect to
      // global coords (format for the latter: nrNodes x spaceDim)
      Matrix<Double> xiDx;
      Vector<Double> ShpFncAtIp;
      ptelem->GetGlobDerivShFncAtIp( xiDx, ip, ptCoord );
      ptelem->GetShFncAtIp( ShpFncAtIp, ip );

      // The matrix aMat can be seen as a numNodes x 1 block-vector.
      // The k-th entry of this block vector corresponds to the matrix
      // A of the ADB product evaluated at the k-th node of the finite
      // element. We assemble aMat in a top-down fashion and assume that
      // the first coordinate equals y and the second z.
      UInt actInd = 0;
      UInt actNode;

      for( actNode = 0; actNode < numNodes; actNode++ ) {

        // 1st row of sub-matrix A(actNode)
        aMat[actInd][0] = xiDx[actNode][0];   // dN/dy
        aMat[actInd][2] = xiDx[actNode][1];   // dN/dz
        actInd++;

        // 2nd row of sub-matrix A(actNode)
        aMat[actInd][1] = xiDx[actNode][1];    // dN/dz
        aMat[actInd][2] = xiDx[actNode][0];    // dN/dy
        actInd++;
      }
    }

    //! Compute the matrix \f$B\f$ of the \f$ADB\f$ operator

    //! The method computes the matrix \f$B\f$ of the \f$ADB\f$ operator
    //! which, in the case of 2D plane strain piezoelectric coupling,
    //! is given  by
    //! \f[
    //! B = \left(\begin{array}{c}
    //! \displaystyle\frac{\partial N_\Phi}{\partial y}\\[3ex]
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
    void calcBMat( Matrix<Double> &bMat, UInt ip,
                   const Matrix<Double> &ptCoord ) {

      ENTER_FCN( "linPiezoCoupling2DPlaneStrain::calcBMat" );

      // Obtain info on number of element's nodes
      const UInt numNodes = ptelem->GetNumNodes();

      // Set correct size of matrix B and initialise with zeros
      bMat.Resize( 2, numNodes );

      // Get derivatives of local shape functions with respect to global
      // coords (format: nrNodes x spaceDim)
      Matrix<Double> xiDx;
      ptelem->GetGlobDerivShFncAtIp( xiDx, ip, ptCoord );

      // The matrix bMat can be seen as a 1 x numNodes block-vector.
      // The k-th entry of this block vector corresponds to the matrix
      // B of the ADB product evaluated at the k-th node of the finite
      // element. We assume that the first coordinate equals y and the
      // second z.
      for( UInt actNode = 0; actNode < numNodes; actNode++ ) {
        bMat[0][actNode] = xiDx[actNode][0];
        bMat[1][actNode] = xiDx[actNode][1];
      }
    }

    //! Compute the data-matrix \f$D\f$

    //! The method computes the matrix D of the piezoelectric coupling
    //! properties of the element's material. In the 2D setting we assume that
    //! the computational 2D area lies in the yz-plane with the x-axis being
    //! direction of large extension that is neglected. Thus D is given by
    //! \f[
    //! D = \left(
    //! \begin{array}{cc}
    //! \displaystyle d_{18} & \displaystyle d_{19} \\
    //! \displaystyle d_{28} & \displaystyle d_{29} \\
    //! \displaystyle d_{38} & \displaystyle d_{39} \\
    //! \displaystyle d_{48} & \displaystyle d_{49} \\
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
    //! \displaystyle e_{34} &           0
    //! \end{array} \right)
    //! \enspace.
    //! \f]
    void calcDMat( Matrix<Double> &dMat ) {
      ENTER_FCN( "linPiezoCoupling2DPlaneStrain::calcDMat" );
      dMat.Resize( 3, 2 );
      Matrix<Double> *matMatrix = ptMaterial->GetMatrix();

      dMat[0][0] = (*matMatrix)[1][7];
      dMat[1][0] = (*matMatrix)[2][7];
      dMat[2][0] = (*matMatrix)[3][7];

      dMat[0][1] = (*matMatrix)[1][8];
      dMat[1][1] = (*matMatrix)[2][8];
      dMat[2][1] = (*matMatrix)[3][8];
    }

    //! Query dimension of the matrix \f$D\f$

    //! This method returns the dimensions of the data-matrix \f$D\f$.
    //! In the case of 2D plane strain piezoelectric coupling we have
    //! \f$\mbox{dim}D=3\times 2\f$
    void getDimD( UInt nRows, UInt nCols ) {
      ENTER_IFCN( "linPiezoCoupling2DPlaneStrain::getDimD" );
      nRows = 3;
      nCols = 2;
    };

    //! Query number of degrees of freedom for first physical quantity

    //! This method can be used to query the number of degrees of freedom at
    //! any node of a finite element for the physical quantity associated to
    //! the base functions whose derivatives form the matrix \f$A\f$.
    //! In the case of 2D plane strain piezoelectric coupling the first
    //! physical quantity is the mechanical displacements.
    //! \return 2
    UInt getNumDofsA() {
      ENTER_IFCN( "linPiezoCoupling2DPlaneStrain::getNumDofsA" );
      return 2;
    }

    //! Query number of degrees of freedom for second physical quantity

    //! This method can be used to query the number of degrees of freedom at
    //! any node of a finite element for the physical quantity associated to
    //! the base functions whose derivatives form the matrix \f$B\f$.
    //! In the case of 2D plain strain piezoelectric coupling the second
    //! physical quantity is the electric potential.
    //! \return 1
    UInt getNumDofsB() {
      ENTER_IFCN( "linPiezoCoupling2DPlaneStrain::getNumDofsB" );
      return 1;
    }

  };

}

#endif
