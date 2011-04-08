#ifndef FILE_LIN_HEATCOND_INT
#define FILE_LIN_HEATCOND_INT

#include <Elements/basefe.hh>
#include <Forms/bdbInt.hh>
#include <Materials/baseMaterial.hh>
#include <General/environment.hh>

namespace CoupledField {


  //! Class describing the general BDB operator for Heat Conduction PDE

  //! The main objective of this class is to implement the pure vitual
  //! methods of the BDBInt parent class for the case of a linear 
  //! heat conduction simulation.
  class linHeatCondInt : public BDBInt {
    
  public:

    // =======================================================================
    // CONSTRUCTION and DESTRUCTION
    // =======================================================================

    //@{ \name Construction and destruction

    //! Constructor with material data
    linHeatCondInt( BaseMaterial* matData, SubTensorType type = FULL,
                bool geoUpdate = false );

    //! Destructor
    ~linHeatCondInt() {

    }
    //@}

    /** @see BaseForm::CalcBMat() */
    void CalcBMat( Matrix<Double> &bMat, UInt ip, const Matrix<Double> &ptCoord );
    
    //! Compute the data-matrix \f$D\f$
    void calcDMat( Matrix<Double> &dMat );
    
    //! Returns dimension of D matrix
    UInt getDimD() {
      return dim_;
    }

    //! Query material type for \f$D\f$ tensor
    MaterialType getDMaterialType() { return HEAT_CONDUCTIVITY_TENSOR; }

    //! Returns nr. of degrees of freedom
    UInt getNrDofs() {
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

