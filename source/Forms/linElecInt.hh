#ifndef FILE_LIN_ELEC_INT
#define FILE_LIN_ELEC_INT

#include <Elements/basefe.hh>
#include <Forms/bdbInt.hh>
#include <Materials/baseMaterial.hh>
#include <General/environment.hh>

namespace CoupledField {


  //! Class describing the general BDB operator for electrostatic

  //! The main objective of this class is to implement the pure vitual
  //! methods of the BDBInt parent class for the case of a linear 
  //! electrostatic simulation.
  class linElecInt : public BDBInt {
    
  public:

    // =======================================================================
    // CONSTRUCTION and DESTRUCTION
    // =======================================================================

    //@{ \name Construction and destruction

    //! Constructor with material data
    linElecInt( BaseMaterial* matData, SubTensorType type = FULL,
                bool geoUpdate = false );

    //! Destructor
    ~linElecInt() {
      ENTER_FCN( "linElecInt::~linElecInt" );

    }
    //@}
    
    // =======================================================================
    // CALCULATION 
    // =======================================================================
    //! Compute the matrix \f$B\f$ of the \f$BDB\f$ operator
    //! \param bMat    (output) computed matrix \f$B\f$
    //! \param ip      (input)  number of integration point
    //! \param ptCoord (input)  matrix containing co-ordinates of all
    //!                         integration points
    void calcBMat( Matrix<Double> &bMat, UInt ip,
                   Matrix<Double> &ptCoord );
    
    //! Compute the data-matrix \f$D\f$
    void calcDMat( Matrix<Double> &dMat );
    
    //! Returns dimension of D matrix
    UInt getDimD() {
      ENTER_IFCN( "linElecInt::getDimD" );
      return dim_;
    }

    //! Query material type for \f$D\f$ tensor
    MaterialType getDMaterialType() { return ELEC_PERMITTIVITY; }

    //! Returns nr. of degrees of freedom
    UInt getNrDofs() {
      ENTER_IFCN( "linElecInt::getNrDofs" );
      return 1;
    }

    //! Sets a multiplicative factor for element matrix
    void SetFactor( const std::string& factor );
    
    //@}
    
  protected:
    
    //! dimension
    UInt dim_;

  };

} // end of namespace

#endif

