// ================================================================================================
/*!
 *       \file     ContactLinInt.hh
 *       \brief    Reference-gap force of penalty contact, for one side of a contact pair.
 *
 *  ------------------------------------------------------------------------------------------
 *  WHAT THIS COMPUTES AND WHY IT IS SEPARATE FROM THE STIFFNESS
 *  ------------------------------------------------------------------------------------------
 *
 *  Split the normal gap into its reference part and its displacement part,
 *
 *      g_N = g_0 + (u_s - u_m).n = g_0 + n^T N_D u        N_D = [ N_s , -N_m ]
 *
 *  Then the penalty contact residual is, for a fixed active set, exactly
 *
 *      r_c = -eps_N INT N_D^T n g_N dG
 *          = -eps_N g_0 INT N_D^T n dG   -   eps_N INT N_D^T n n^T N_D dG . u
 *            \________ THIS CLASS _______/   \________ ContactABInt _______/
 *
 *  StdSolveStep::StepStaticNonLin() builds the residual as  f - K(u).u  from the non-Newton
 *  matrices, so the K_c part is produced automatically by registering ContactABInt into the
 *  non-Newton STIFFNESS set. The g_0 part is a plain right-hand-side vector and has to be
 *  supplied here, there is no matrix that can generate a term independent of u.
 *
 *  g_0 vanishes only when the two surfaces are coincident at u = 0. Any assembly with a real
 *  clearance or interference has g_0 != 0, and then THIS TERM CARRIES THE ENTIRE INITIAL LOAD.
 *
 *  ------------------------------------------------------------------------------------------
 *  BLOCK STRUCTURE
 *  ------------------------------------------------------------------------------------------
 *
 *  N_D^T = [ N_s^T ; -N_m^T ], so the two sides get opposite signs:
 *
 *      secondary rows:  -eps_N w g_0 (N_s^T n)
 *      primary rows:    +eps_N w g_0 (N_m^T n)
 *
 *  One instance produces one side, selected by the constructor flag, matching the
 *  SurfaceLinFormContext that wraps it. The sign convention is the same rowSign that
 *  ContactABInt uses for its rows: +1 on the secondary side, -1 on the primary side.
 *
 *  ------------------------------------------------------------------------------------------
 *  MODE_FRICTION -- THE WHOLE TANGENTIAL RESIDUAL, NOT JUST A CONSTANT PART
 *  ------------------------------------------------------------------------------------------
 *
 *  The normal split above works because the penalty term is exactly quadratic in u, so all but
 *  a constant of it can be manufactured as -K_c.u from a non-Newton stiffness. The tangential
 *  term is not quadratic, on a slipping point the traction is mu*p*e, nonlinear in u through
 *  both the pressure and the direction, so there is nothing for a stiffness to manufacture.
 *
 *  MODE_FRICTION therefore carries the COMPLETE friction residual
 *
 *      r_T = -INT N_D^T t_T dG
 *
 *  with t_T the traction the return map produces (stick or slip), and ContactABInt's
 *  MODE_FRICTION supplies the matching consistent tangent as a NEWTON form.
 */
//================================================================================================

#ifndef FILE_CFS_CONTACTLININT_HH
#define FILE_CFS_CONTACTLININT_HH

#include "Forms/LinForms/LinearForm.hh"
#include "Domain/Mesh/NcInterfaces/ContactInterface.hh"

namespace CoupledField {

class ContactLinInt : public LinearForm {

  public:

    enum Mode {
      //! The constant reference-gap force -(lambda + eps_N g_0) INT N_D^T n dG.
      MODE_GAP_FORCE,
      //! The complete tangential (friction) residual -INT N_D^T t_T dG.
      MODE_FRICTION
    };

    ContactLinInt(ContactInterface* iface, bool rowIsPrimary,
                  Mode mode = MODE_GAP_FORCE);

    ContactLinInt(const ContactLinInt& right);

    virtual ~ContactLinInt() {}

    virtual ContactLinInt* Clone() { return new ContactLinInt(*this); }

    virtual bool IsComplex() const { return false; }

    //! Element vector of this side on one contact segment.
    void CalcElemVector(Vector<Double>& elemVec, EntityIterator& ent) override;

    //! The active set depends on the solution, so this form must be re-assembled every
    //! Newton iteration - which is what puts it into Assemble::AssembleNonLinRHS().
    bool IsSolDependent() override { return true; }

  private:

    //! the contact pair supplying the integration points
    ContactInterface* iface_;

    //! which body this instance assembles onto
    bool rowIsPrimary_;

    //! reference-gap force or friction residual
    Mode mode_;

    Matrix<Double> bMat_;
    Vector<Double> slip_, tracT_;
};

} // namespace CoupledField

#endif /* FILE_CFS_CONTACTLININT_HH */
