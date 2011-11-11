#ifndef FILE_PSEUDOTS_2007
#define FILE_PSEUDOTS_2007

#include "TimeStepping.hh"

namespace CoupledField
{

  //! class for time stepping of hyperbolic PDE: method is PseudoTS

  class PseudoTS: public TimeStepping {
  public:
    //! constructor
    //! \param algebraicsystem pointer to algebraic system 
    PseudoTS( AlgebraicSys * algebraicsystem );
    
    //! destructor
    virtual ~PseudoTS();
  
    //! initilization
    //! \param rhsSize total number of entries in the rhs vector
    void Init( Double dt, UInt rhsSize );

    //! perform predictor step
    void Predictor(SBM_Vector& solold);

    //! perform predictor step
    void RelaxedPredictor(SBM_Vector& solold, 
                          SBM_Vector& solpredRelaxed, 
                          SBM_Vector& solderiv1predRelaxed, 
                          Double omega){;};

    //! perform corrector step
    void Corrector(SBM_Vector& solnew);

    //! perform an update to RHS
    void UpdateRHS(){;};

    //! perform an update to RHS with actual solution (for nonlin calculation)
    void UpdateRHS(SBM_Vector& actSol){;};  

    //! perform calculations at end of timestep (set the new timestep)
    void AdvanceTimestep(SBM_Vector& solnew);

  private:
    //! predictor for nodal solution
    SBM_Vector solold_;

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
