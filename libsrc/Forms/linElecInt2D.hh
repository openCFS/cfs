#ifndef FILE_LIN_ELEC_INT_2D
#define FILE_LIN_ELEC_INT_2D

namespace CoupledField {

  //! Class describing the general BDB operator for 2D electrostatic

  //! The main objective of this class is to implement the pure vitual
  //! methods of the BDBInt parent class for the case of a linear 
  //! electrostatic 2D simulation.
  //! \Note: The class asssumes that the computational 2D area lies in the
  //! yz-plane.
  class linElecInt2D : public BDBInt {
    
  public:

    // =======================================================================
    // CONSTRUCTION and DESTRUCTION
    // =======================================================================

    //@{ \name Construction and destruction

    //! Default constructor
    linElecInt2D( Boolean isAxi) {
      ENTER_FCN( "linElecInt2D::linElecInt2D" );
      
      isaxi_ = isAxi;
      factor_ = 0.0;
    }

    //! Constructor with material data
    linElecInt2D( MaterialData & matData, Boolean isAxi)  :
      BDBInt( matData )
    {
      ENTER_FCN( "linElecInt2D::linElecInt2D" );

      factor_ = 0.0;
      isaxi_ = isAxi;
    }
    

    //! Destructor
    ~linElecInt2D() {
      ENTER_FCN( "linElecInt2D::~linElecInt2D" );

    }
    //@}
    
    // =======================================================================
    // CALCULATION 
    // =======================================================================
    
    //@{ \name Calculation Methods

     //! Compute the matrix \f$B\f$ of the \f$BDB\f$ operator

    //! The method computes the matrix \f$B\f$ of the \f$BDB\f$ operator
    //! which, in the case of a 2D electrostatic field,
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
    void calcBMat( Matrix<Double> &bMat, Integer ip,
                   Matrix<Double> &ptCoord ) {

      ENTER_FCN( "linElecInt2D::calcBMat" );

      // Obtain info on number of element's nodes
      const Integer numNodes = ptelem->GetNumNodes();

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
      for( Integer actNode = 0; actNode < numNodes; actNode++ ) {
        bMat[0][actNode] = xiDx[actNode][0];
	bMat[1][actNode] = xiDx[actNode][1];
      }
    }
    
    //! Compute the data-matrix \f$D\f$
    
    //! The method computes the matrix D of the electrostatic
    //! properties of the element's material. 
    //! In the 2D setting we assume that
    //! the computational 2D area lies in the yz-plane 
    //! Therefore D looks like
    //! \f[
    //! D = \left(
    //! \begin{array}{cc}
    //! \displaystyle \epsilon_{11} & \displaystyle \epsilon_{12} \\
    //! \displaystyle \epsilon_{21} & \displaystyle \epsilon_{22} \\
    //! \end{array} \right)
    //! \f]
    //! where \f$\epsilon\f$ is the local tensor of dielectric constants.
    void calcDMat( Matrix<Double> &dMat ) {
      ENTER_FCN( "linElecInt2D::calcDMat" );
      dMat.Resize( 2, 2 );
      Matrix<Double> *matMatrix = ptMaterial->GetMatrix();
      
      // copy electric part of material matrix, which 
      // is the lower-right sub-diagonal block
      // d[8-9][8-9]
      Integer startRow = 7;
      Integer startCol = 7;
      for( Integer i = 0; i < 2; i++ ) {
	for ( Integer j = 0; j < 2; j++ ) {
	  dMat[i][j] = factor_ *
	    (*matMatrix)[startRow+i][startCol+j];
	}
      }
    }
    
    //! Returns dimension of D matrix
    Integer getDimD() {
      ENTER_IFCN( "linElecInt2D::getDimD" );
      return 2;
    }

    //! Returns nr. of degrees of freedom
    Integer getNrDofs() {
      ENTER_IFCN( "linElecInt2D::getNrDofs" );
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

