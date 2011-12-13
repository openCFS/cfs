// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NLIN_ELEC_HYST_INT
#define FILE_NLIN_ELEC_HYST_INT

#include <string>

#include "Forms/bdbInt.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "MatVec/vector.hh"

namespace CoupledField {


  //! Class describing the general BDB operator for electrostatic

  //! The main objective of this class is to implement the pure vitual
  //! methods of the BDBInt parent class for the case of a hysteretic 
  //! electrostatic simulation.
class BaseMaterial;
class EntityIterator;
class EqnMap;
class Grid;
class StdPDE;
struct ResultInfo;
template <class TYPE> class GradientFieldOp;
template <class TYPE> class Matrix;

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

    }
    //@}
    
    // =======================================================================
    // CALCULATION 
    // =======================================================================

    //! Compute element matrix associated to ADB form
    void CalcElementMatrix( Matrix<Double>& elemMat,
			    EntityIterator& ent1, 
			    EntityIterator& ent2 );

    /** @see BaseForm::CalcBMat() */
    void CalcBMat( Matrix<Double> &bMat, UInt ip,
                   const Matrix<Double> &ptCoord );
    
    //! Compute the data-matrix \f$D\f$
    void calcDMat( Matrix<Double> &dMat );
    
    //! Returns dimension of D matrix
    UInt getDimD() {
      return dim_;
    }

    //! Query material type for \f$D\f$ tensor
    MaterialType getDMaterialType() { return ELEC_PERMITTIVITY; }


    //! Returns nr. of degrees of freedom
    UInt getNrDofs() {
      return 1;
    }
        
    //@}

    //! set objects for computation of E-field
    void Set4Hyst(Grid * ptGrid, 
		  StdPDE* ptPDE,
		  shared_ptr<EqnMap> eqnMap,
		  shared_ptr<ResultInfo> result);

    //! Sets a multiplicative factor for element matrix
    void SetFactor( const std::string& factor );
      
  protected:
    
    //! dimension
    UInt dim_;

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

