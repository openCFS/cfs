#ifndef FILE_BDF2_2007
#define FILE_BDF2_2004

#include <General/environment.hh>
#include "timestepping.hh"


namespace CoupledField
{

  //! class for time stepping of parabolic PDE: method is Trapezoidal

  class Bdf2: public TimeStepping
  {
  public:
    //! constructor
    //! \param algebraicsystem pointer to algebraic system
    Bdf2(  AlgebraicSys * algebraicsystem );

    //! destructor
    virtual ~Bdf2();

    //! initilization
    //! \param rhsSIze total number of entries in the rhs vector
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
    void UpdateRHS();

    //! perform calculations at end of timestep (set the new timestep)
    void AdvanceTimestep(SBM_Vector& solnew);

    //! Substract Stiffness from RHS
    virtual void SubstractStiffnessFromRHS(SBM_Vector& actSol)
    {;};

  private:

    //! compute parameters for multiplication
    void CalcParameters(Double dt);
    bool firstTime_;
    SBM_Vector coeffMass_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class
  // =========================================================================

  //! \class Bdf2
  //!
  //! \purpose
  //!
  //! \collab
  //!
  //! \implement
  //!
  //! \status In use.
  //!
  //! \unused
  //!
  //! \improve
  //!

#endif

}

#endif // FILE_BDF2
