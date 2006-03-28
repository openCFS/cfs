#ifndef FILE_LIN_ELEC_INT_3D
#define FILE_LIN_ELEC_INT_3D

namespace CoupledField {


  //! Class describing the general BDB operator for 3D electrostatic

  //! The main objective of this class is to implement the pure vitual
  //! methods of the BDBInt parent class for the case of a linear 
  //! electrostatic 3D simulation.
  class linElecInt3D : public BDBInt {
    
  public:

    // =======================================================================
    // CONSTRUCTION and DESTRUCTION
    // =======================================================================

    //@{ \name Construction and destruction

    //! Default constructor
    linElecInt3D() {
      ENTER_FCN( "linElecInt3D::linElecInt3D" );
      
      factor_ = 0.0;
    }

    //! Constructor with material data
    linElecInt3D( BaseMaterial* matData)  :
      BDBInt(matData)
    {
      ENTER_FCN( "linElecInt3D::linElecInt3D" );

      factor_ = 0.0;
    }
    

    //! Destructor
    ~linElecInt3D() {
      ENTER_FCN( "linElecInt3D::~linElecInt3D" );

    }
    //@}
    
    // =======================================================================
    // CALCULATION 
    // =======================================================================
    
    //@{ \name Calculation Methods

    //! Compute the matrix \f$B\f$ of the \f$BDB\f$ operator

    //! The method computes the matrix \f$B\f$ of the \f$BDB\f$ operator
    //! which, in the case of a 3D electrostatic field,
    //! is given  by
    //! \f[
    //! B = \left(\begin{array}{c}
    //! \displaystyle\frac{\partial N_\Phi}{\partial x}\\[3ex]
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
    //! \in \mathcal{R}^{3 \times N}
    //! \f]
    //! where \f$p_k, k=1,\ldots,N\f$ are the nodes of the element.
    //! \param bMat    (output) computed matrix \f$B\f$
    //! \param ip      (input)  number of integration point
    //! \param ptCoord (input)  matrix containing co-ordinates of all
    //!                         integration points
    void calcBMat( Matrix<Double> &bMat, UInt ip,
                   Matrix<Double> &ptCoord ) {

      ENTER_FCN( "linElecInt3D::calcBMat" );

      // Obtain info on number of element's nodes
      const UInt numNodes = ptelem->GetNumNodes();

      // Set correct size of matrix B and initialise with zeros
      bMat.Resize( 3, numNodes );

      // Get derivatives of local shape functions with respect to global
      // coords (format: nrNodes x spaceDim)
      Matrix<Double> xiDx;
      ptelem->GetGlobDerivShFncAtIp( xiDx, ip, ptCoord );

      // The matrix bMat can be seen as a 3 x numNodes block-vector.
      // The k-th entry of this block vector corresponds to the matrix
      // B of the BDB product evaluated at the k-th node of the finite
      // element. 
      for( UInt actNode = 0; actNode < numNodes; actNode++ ) {
        bMat[0][actNode] = xiDx[actNode][0];
        bMat[1][actNode] = xiDx[actNode][1];
        bMat[2][actNode] = xiDx[actNode][2];
      }
    }
    
    //! Compute the data-matrix \f$D\f$

    //! The method computes the matrix D of the electrostatic
    //! properties of the element's material. In the 3D setting the latter
    //! is given by
    //! \f[
    //! D = \left( \epsilon \right) \in \mathcal{R}^{3\times 3}.
    //! \f]
    //! where \f$\epsilon\f$ is the local tensor of dielectric constants
    void calcDMat( Matrix<Double> &dMat ) {
      ENTER_FCN( "linElecInt3D::calcDMat" );
      
      ptMaterial->GetTensor(dMat,ELEC_PERMITTIVITY,REAL);
      dMat *= factor_;
    }
    
    //! Returns dimension of D matrix
    UInt getDimD() {
      ENTER_IFCN( "linElecInt3D::getDimD" );
      return 3;
    }

    //! Returns nr. of degrees of freedom
    UInt getNrDofs() {
      ENTER_IFCN( "linElecInt3D::getNrDofs" );
      return 1;
    }
    

    //! Sets a multiplicative factor for element matrix
    void SetFactor( const Double factor ) {
      factor_ = factor;
    }
    
    //@}
    
  protected:
    
    //! multiplicative factor for material
    Double factor_;
    
  };

} // end of namespace

#endif

