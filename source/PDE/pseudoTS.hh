#ifndef FILE_PSEUDOTS_2007
#define FILE_PSEUDOTS_2007

#include "General/defs.hh"
#include "MatVec/vector.hh"
#include "timestepping.hh"

namespace CoupledField {
class BaseSystem;
}  // namespace CoupledField

namespace CoupledField
{

  //! class for time stepping of hyperbolic PDE: method is PseudoTS

  class PseudoTS: public TimeStepping {
  public:
    //! constructor
    //! \param algebraicsystem pointer to algebraic system 
    PseudoTS( BaseSystem * algebraicsystem );
    
    //! destructor
    virtual ~PseudoTS();
  
    //! initilization
    //! \param rhsSize total number of entries in the rhs vector
    void Init( Double dt, UInt rhsSize );

    //! perform predictor step
    void Predictor(Vector<Double>& solold);

    //! perform predictor step
    void RelaxedPredictor(Vector<Double>& solold, 
                          Vector<Double>& solpredRelaxed, 
                          Vector<Double>& solderiv1predRelaxed, 
                          Double omega){;};

    //! perform corrector step
    void Corrector(Vector<Double>& solnew);

    //! perform an update to RHS
    void UpdateRHS(){;};

    //! perform an update to RHS with actual solution (for nonlin calculation)
    void UpdateRHS(Vector<Double>& actSol){;};  

    //! perform calculations at end of timestep (set the new timestep)
    void AdvanceTimestep(Vector<Double>& solnew);

  private:
    //! predictor for nodal solution
    Vector<Double> solold_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class PseudoTS
  //! 
  //! \purpose 
  //! 
  //! \collab 
  //! 
  //! \implement 
  //! 
  //! \status In use
  //! 
  //! \unused 
  //! 
  //! \improve
  //! 

#endif

}

#endif // FILE_PSEUDOTS
