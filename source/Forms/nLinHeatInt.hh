#ifndef FILE_NLIN_HEAT_INT
#define FILE_NLIN_HEAT_INT

#include "Forms/linHeatCondInt.hh"
#include "Forms/linearForm.hh"

#include "General/defs.hh"
#include "General/environment.hh"
#include "MatVec/vector.hh"

namespace CoupledField {


  //! heat conduction: nonlinear stiffness integrator 
  //! (heat conductivity is a function of temperature)
class ApproxData;
class BaseMaterial;
class EntityIterator;
template <class TYPE> class Matrix;

  class nlinHeatStiffInt : public linHeatCondInt {
    
  public:

    // =======================================================================
    // CONSTRUCTION and DESTRUCTION
    // =======================================================================

    //@{ \name Construction and destruction

    //! Constructor with material data
    nlinHeatStiffInt( BaseMaterial* matData, SubTensorType type = FULL,
                bool geoUpdate = false );

    //! Destructor
    ~nlinHeatStiffInt() {
    }
    //@}


    /// Calculation of stiffmess matrix
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );
    
  protected:
    
  private:
    
    ApproxData *nlinFnc_;      //!< pointer to approximation object
    Vector<Double> temperature_;    //!< magnetic vector potential at nodes

  };



  /// heat conduction: nonlinear mass integrator 
  /// (heat capacity is a function of temperature)
  class nlinHeatMassInt : public linHeatCondInt {
    
  public:

    // =======================================================================
    // CONSTRUCTION and DESTRUCTION
    // =======================================================================

    //@{ \name Construction and destruction

    //! Constructor with material data
    nlinHeatMassInt( BaseMaterial* matData, SubTensorType type = FULL,
                bool geoUpdate = false );

    //! Destructor
    ~nlinHeatMassInt() {
    }
    //@}


    /// Calculation of stiffmess matrix
    void CalcElementMatrix( Matrix<Double>& elemMat,
                            EntityIterator& ent1, 
                            EntityIterator& ent2 );
    
  protected:
    
  private:
    
    ApproxData *nlinFnc_;      //!< pointer to approximation object
    Vector<Double> temperature_;    //!< magnetic vector potential at nodes
    Double density_;

  };


  // =============================================================================
  // nonlinear RHS
  // =============================================================================

  /// class for calculation of right-hand-side of nonlinear heat PDE
  class nLinHeat_linFormInt : public LinearForm
  {
  public:

    /// constructor
    nLinHeat_linFormInt( BaseMaterial* matData, bool axi=false, 
                         bool coordUpdate = false );

    /// destructor
    virtual ~ nLinHeat_linFormInt();

    /// Calculation of vector of right hand side 
    void CalcElemVector( Vector<Double> & result,
                         EntityIterator& ent );
  
  protected:

    ApproxData *nlinFnc_;
  };

} // end of namespace

#endif

