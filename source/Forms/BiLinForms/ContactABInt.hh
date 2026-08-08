// ================================================================================================
/*!
 *       \file     ContactABInt.hh
 *       \brief    Penalty contact stiffness for one coupling block of a contact pair.
 *
 *  ------------------------------------------------------------------------------------------
 *  WHY THIS IS A STIFFNESS AND NOT A NEWTON FORM
 *  ------------------------------------------------------------------------------------------
 *
 *  With the gap split into its reference part and its displacement part,
 *
 *      g_N = g_0 + (u_s - u_m).n
 *
 *  the penalty contact contribution to the residual is *exactly quadratic* in u for a fixed
 *  active set:
 *
 *      r_c = -eps_N g_0 INT N_D^T n dG   -   eps_N INT N_D^T n n^T N_D dG . u
 *            \_____ constant force _____/     \_________ K_c . u __________/
 *
 *  with N_D the operator mapping the combined DOF vector to (u_s - u_m).
 *
 *  StdSolveStep::StepStaticNonLin() forms the residual as f - K(u).u using the NON-Newton
 *  matrices, and adds the Newton-flagged forms to the tangent only. Putting K_c into the
 *  non-Newton STIFFNESS set therefore produces BOTH the -K_c.u residual term AND the correct
 *  tangent contribution, because -dr_c/du = K_c exactly. A separate Newton form would
 *  double-count.
 *
 *  The only nonlinearity left is the active set, so the form is marked solution dependent and
 *  gets re-assembled every Newton iteration.
 *
 *  The constant force is NOT this class's business: it is a right-hand-side vector and lives in
 *  ContactLinInt (Forms/LinForms). It vanishes only when the two surfaces are coincident at
 *  u = 0; for any assembly with a real clearance or interference it carries the whole initial
 *  load, so the two classes must always be registered together. They share
 *  CalcContactInterpolationMatrix() (Forms/ContactPointInterpolation.hh) and the same
 *  active-set test (ContactInterface::EvalCurrentGap) so that the residual and the tangent
 *  cannot describe different discrete gaps.
 *
 *  ------------------------------------------------------------------------------------------
 *  BLOCK STRUCTURE
 *  ------------------------------------------------------------------------------------------
 *
 *  N_D = [ N_s , -N_m ], so K_c has four blocks whose signs follow from that:
 *
 *      (s,s): +eps w (N_s^T n)(n^T N_s)      (s,m): -eps w (N_s^T n)(n^T N_m)
 *      (m,s): -eps w (N_m^T n)(n^T N_s)      (m,m): +eps w (N_m^T n)(n^T N_m)
 *
 *  One instance of this class produces one block, selected by the CouplingDirection it is
 *  constructed with, matching the SurfaceBiLinFormContext that wraps it.
 *
 *  ------------------------------------------------------------------------------------------
 *  MODE_STABILIZATION
 *  ------------------------------------------------------------------------------------------
 *
 *  A body supported ONLY through a contact interface, starting SEPARATED, has no stiffness at
 *  all in the direction it must travel: at u = 0 no point is active, so K_c is empty and the
 *  body's normal rigid-body mode is unconstrained. The Newton step is then whatever the linear
 *  solver makes of a numerically singular system, and the iteration does not so much diverge as
 *  wander.
 *  MODE_STABILIZATION adds a small artificial normal stiffness
 *
 *      K_stab = eps_stab INT N_D^T n n^T N_D dG        eps_stab = stabilization * E / h
 *
 *  on the points that are PROJECTED BUT NOT ACTIVE. It removes the singular direction, so every Newton
 *  step is bounded and the iterate
 *  walks towards contact instead of flailing. Once a point closes it becomes active, gets the
 *  real eps_N stiffness from MODE_CONTACT, and is skipped here, the stabilisation switches
 *  itself off as contact establishes, with no schedule to tune.
 *
 *  IT MUST BE A NEWTON FORM, and this is the whole point of the design:
 *
 *  StdSolveStep builds the residual as f - K.u from the NON-Newton matrices and adds the
 *  Newton-flagged ones to the tangent ONLY. A stabilisation stiffness in the non-Newton set
 *  would therefore contribute -K_stab.u to the residual, i.e. a fictitious force, and would
 *  change the converged solution. As a Newton form it touches the tangent alone, so
 *  the residual, and hence the solution it defines, is provably identical with and without it.
 *
 *  ------------------------------------------------------------------------------------------
 *  MODE_FRICTION
 *  ------------------------------------------------------------------------------------------
 *
 *  The normal penalty term is exactly quadratic in u for a fixed active set, which is what lets
 *  it live in the non-Newton STIFFNESS set and generate its own residual as -K_c.u. The
 *  tangential term is NOT on a slipping point
 *
 *      t_T = mu p e ,      e = t_T^trial/||t_T^trial|| ,      p = -t_N
 *
 *  and both the magnitude (through p, hence through the normal gap) and the direction e depend
 *  on u nonlinearly. There is no quadratic form to hide it in.
 *
 *  So friction uses the OTHER, general split:
 *
 *      residual   r_T = -INT N_D^T t_T dG        ContactLinInt, MODE_FRICTION, sol dependent
 *      tangent    K_T = +INT N_D^T D N_D dG      THIS CLASS, MODE_FRICTION, NEWTON form
 *
 *  with D = dt_T/d(u_s - u_m). That is uniformly correct for both branches, so both are handled
 *  here rather than splitting stick (which happens to be quadratic) away from slip. The stick
 *  branch could have gone into the non-Newton set; keeping the two branches in one place is
 *  worth more than the assembly it would save, because a point crossing the cone would otherwise
 *  have to migrate between two differently-flagged integrators within one Newton iteration.
 *
 *  The two branches of D
 *
 *      STICK   D = eps_T P_T                                       symmetric, PSD
 *      SLIP    D = -mu eps_N e n^T  +  (mu p/||s||) (P_T - e e^T)   NON-SYMMETRIC
 *
 *  with P_T = I - n n^T, s the elastic slip and e = s/||s||. The first slip term is the
 *  derivative of the Coulomb limit with respect to the normal gap, and it is the reason
 *  frictional contact has a non-symmetric tangent: e lies in the tangent plane and n does not,
 *  so e n^T has no symmetric counterpart anywhere in the system. openCFS handles this without
 *  configuration -- BiLinearForm defaults to isSymmetric_ = false, so the contact blocks already
 *  select a sparseNonSym matrix and a direct LU solver, verified on the frictionless testcases.
 *
 *  In 2D the second slip term vanishes identically: the tangent plane is one-dimensional, so e
 *  spans it and P_T = e e^T exactly. A 2D slip tangent is therefore purely the non-symmetric
 *  normal-coupling term, which is worth knowing before reading a 2D convergence history.
 * */
//================================================================================================

#ifndef FILE_CFS_CONTACTABINT_HH
#define FILE_CFS_CONTACTABINT_HH

#include "Forms/BiLinForms/BiLinearForm.hh"
#include "Domain/Mesh/NcInterfaces/ContactInterface.hh"

namespace CoupledField {

class BaseFeFunction;
class ContactABInt : public BiLinearForm {

  public:

    enum Mode {
      MODE_CONTACT,
      MODE_STABILIZATION,
      MODE_FRICTION
    };

    ContactABInt(ContactInterface* iface, BiLinearForm::CouplingDirection direction,
                 Mode mode = MODE_CONTACT);

    ContactABInt(const ContactABInt& right);

    virtual ~ContactABInt() {}

    virtual ContactABInt* Clone() { return new ContactABInt(*this); }

    virtual bool IsComplex() const { return false; }

    //! Element matrix for this block on one contact segment.
    void CalcElementMatrix(Matrix<Double>& elemMat,
                           EntityIterator& ent1,
                           EntityIterator& ent2) override;

    void SetFeSpace(shared_ptr<FeSpace> feSpace1, shared_ptr<FeSpace> feSpace2) override;

    void SetFeSpace(shared_ptr<FeSpace> feSpace) override { SetFeSpace(feSpace, feSpace); }

    bool IsSolDependent() override { return true; }

  private:

    //! the contact pair supplying the integration points
    ContactInterface* iface_;

    //! which block this instance computes
    BiLinearForm::CouplingDirection direction_;

    //! contact stiffness or separated-start stabilisation
    Mode mode_;

    //! scratch, kept as members to avoid reallocating per integration point
    Matrix<Double> bMatRow_, bMatCol_;
    Vector<Double> nRow_, nCol_;
    //! friction scratch: the material tangent D = dt_T/d(u_s - u_m), and the slip vectors
    Matrix<Double> dMat_;
    Vector<Double> slip_, tracT_;
};

} // namespace CoupledField

#endif /* FILE_CFS_CONTACTABINT_HH */
