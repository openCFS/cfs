// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_LIN_GRAD_BDB_INT
#define FILE_LIN_GRAD_BDB_INT

#include <Elements/basefe.hh>
#include <Forms/bdbInt.hh>
#include <Materials/baseMaterial.hh>
#include <General/environment.hh>

namespace CoupledField {


  //! Specialized BDB-Integrator with gradient operator as B-differential operator

  //! This class implements the BDB-integrator, where the B-operator is the 
  //! gradient operator and the D-tensor can be set variable.
  //! It is used e.g. for the electrostatic problem, as well as the magnetic
  //! scalar problem for calculating the stiffness matrix.
  class linGradBDBInt : public BDBInt {
    
  public:

    // =======================================================================
    // CONSTRUCTION and DESTRUCTION
    // =======================================================================

    //@{ \name Construction and destruction

    //! Constructor with material data
    linGradBDBInt( BaseMaterial* matData,
                   MaterialType matType,
                   SubTensorType type = FULL,
                   bool geoUpdate = false );

    //! Destructor
    ~linGradBDBInt() {

    }
    //@}
    
    /** @see CalcBMat() */
    void CalcBMat( Matrix<Double> &bMat, UInt ip,
                   const Matrix<Double> &ptCoord );
    
    //! Compute the data-matrix \f$D\f$
    void calcDMat( Matrix<Double> &dMat, const Elem* elem);
    
    //! Returns dimension of D matrix
    UInt getDimD() {
      return dim_;
    }

    //! Query material type for \f$D\f$ tensor
    MaterialType getDMaterialType() { return matType_; }

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
    
    //! material type of D-matrix/tensor
    MaterialType matType_;

  };

} // end of namespace

#endif

