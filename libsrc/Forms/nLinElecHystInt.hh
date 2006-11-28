#ifndef FILE_NLIN_ELEC_HYST_INT
#define FILE_NLIN_ELEC_HYST_INT

#include <Elements/basefe.hh>
#include <Forms/bdbInt.hh>
#include <Forms/gradfieldop.hh>
#include <Materials/baseMaterial.hh>
#include <General/environment.hh>

namespace CoupledField {


  //! Class describing the general BDB operator for electrostatic

  //! The main objective of this class is to implement the pure vitual
  //! methods of the BDBInt parent class for the case of a hysteretic 
  //! electrostatic simulation.
  class nlinElecHystInt : public BDBInt {
    
  public:

    // =======================================================================
    // CONSTRUCTION and DESTRUCTION
    // =======================================================================

    //@{ \name Construction and destruction

    //! Constructor with material data
    nlinElecHystInt( BaseMaterial* matData, SubTensorType type = FULL,
                bool geoUpdate = false );

    //! Destructor
    ~nlinElecHystInt() {
      ENTER_FCN( "nlinElecHystInt::~nlinElecHystInt" );

    }
    //@}
    
    // =======================================================================
    // CALCULATION 
    // =======================================================================

    //! Compute element matrix associated to ADB form
    void CalcElementMatrix( Matrix<Double>& elemMat,
			    EntityIterator& ent1, 
			    EntityIterator& ent2 );

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
      ENTER_IFCN( "nlinElecHystInt::getDimD" );
      return dim_;
    }

    //! Returns nr. of degrees of freedom
    UInt getNrDofs() {
      ENTER_IFCN( "linElecInt::getNrDofs" );
      return 1;
    }
    

    //! Sets a multiplicative factor for element matrix
    void SetFactor( const Double factor ) {
      factor_ = factor;
    }
    
    //@}

    //! set objects for computation of E-field
    void Set4Hyst(Grid * ptGrid, 
		  StdPDE* ptPDE,
		  shared_ptr<EqnMap> eqnMap,
		  shared_ptr<ResultDof> result);

    
  protected:
    
    //! dimension
    UInt dim_;

    //! multiplicative factor for material
    Double factor_;

    /// scalar electric potential of all nodes of actual element
    Vector<Double> elemPot_;

    //! electric field operator
    GradientFieldOp<Double> * EfieldOp_;

    //! diferential electric permittivity dD/dE
    Double diffEpsVal_;    

    //! direction of ploarization
    UInt dirP_;
  };

} // end of namespace

#endif

