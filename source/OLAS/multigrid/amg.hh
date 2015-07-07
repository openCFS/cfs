// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id$ */

#ifndef OLAS_AMG_HH
#define OLAS_AMG_HH
/**********************************************************/

#include "multigrid/hierarchylevel.hh"

namespace CoupledField {
/**********************************************************/

//! AMG solver or preconditioner

/*! This class encapsulates a AMG matrix hierarchy and other needed
 *  stuff for an AMG used as solver or preconditioner. For usage create
 *  an AMGSolver object, call Setup(), then Cycle().
 *  <br><br>
 *  <b>Settings:</b><br>
 *  There are some settings, the behaviour of the AMG components can be
 *  be controlled with. These parameters are passed to the constructor
 *  as a pointer to an OLAS_Params object. If a NULL pointer is passed,
 *  some default values are used (see below). The settings, that are
 *  specified in this object are copied into a temporary object of type
 *  HierarchyLevel::Settings, which is passed during the setup phase
 *  through all hierarchy levels. This copy prevents the frequent usage
 *  of the Get?Value routines of OLAS_Params, that might involve
 *  performance loss in heavy use.
 *  \note All settings for serial AMG are a subset of the settings
 *        for parallel AMG and therefore also control parallel AMG.
 *  <br><br>
 *  <b>List of Settings</b><br>
 *  <table>
 *     <tr>
 *        <td><b>ID string in OLAS_Params</b></td>
 *        <td><b>member in HierarchyLevel::Settings</b></td>
 *        <td><b>default value (params = NULL)</b></td>
 *        <td><b>description</b></td>
 *     </tr><tr>
 *        <td>AMG_MaxCoarseDependency</td>
 *        <td>maxCoarseDepend</td>
 *        <td>10</td>
 *        <td>Setup phase, C-F-splitting. The maximal number of coarse
 *            neighbours, an F-point can get due to their strong influence
 *            on that point. This does not mean, that this will be the
 *            number of C-neighbours, that will maximally surround an
 *            F-point, but the number of C-neighbours, that will cause
 *            an F-point to be judged as <i>well interpolated</i>. Thus
 *            a further neighbour will not be made a C-point due to its
 *            dependency to the named point, but it could still get C-point
 *            because of other reasons. (used in
 *            Topology::CalcCoarseFineSplitting). </td>
 *     </tr><tr>
 *        <td>AMG_MinSystemSize</td>
 *        <td>minSystemSize</td>
 *        <td>20</td>
 *        <td>Setup phase, C-F-splitting. The minimal size for the
 *            coarsest system. If the coarsening process yields a number
 *            of C-points for the next coarser system that falls below
 *            this limit, the coarsening result is rejected and the
 *            corresponding hierarchy level is not created. The present
 *            level gets the coarsest. (used in HierarchyLevel::Setup)</td>
 *     </tr><tr>
 *        <td>AMG_NumPreSmoothing</td>
 *        <td>numPreSmoothings</td>
 *        <td>1</td>
 *        <td>Solution Phase. Number of pre-smoothing steps in the AMG
 *            solution cycle (used in HierarchyLevel::Cycle). </td>
 *     </tr><tr>
 *        <td>AMG_NumPostSmoothing</td>
 *        <td>numPostSmoothings</td>
 *        <td>1</td>
 *        <td>Solution Phase. Number of post-smoothing steps in the AMG
 *            solution cycle (used in HierarchyLevel::Cycle). </td>
 *     </tr><tr>
 *        <td>AMG_CycleParameter</td>
 *        <td>CycleParameter</td>
 *        <td>1</td>
 *        <td>Solution phase. The recursion parameter for the solution
 *            cycle, which is the number of coarse corrections, applied
 *            <b>on each level except the top and the bottom level</b>.
 *            A value of 1 for example causes a common V-cycle, 2 a
 *            common W-cycle. Note that excluding the top level from
 *            this multiple coarse grid corrections causes the parameter
 *            value 2 to generate a real W-cycle, not a WW-cycle.
 *            (used in HierarchyLevel::Cycle) </td>
 *     </tr><tr>
 *        <td>AMG_NumBadCoarsenings</td>
 *        <td>numBadCoarsenings</td>
 *        <td>5</td>
 *        <td>Setup phase, coarsening The recursive creation of coarser
 *            levels is interrupted, if the coarsening gets slow. "Getting
 *            slow" is defined as "the ratio of coarse to fine grid size
 *            was greater than \c AMG_BadCoarseningRatio in
 *            \c AMG_NumBadCoarsenings+1 contiguous coarsenings".</td>
 *     </tr><tr>
 *        <td>AMG_CoarseningAlpha</td>
 *        <td>alpha</td>
 *        <td>0.15</td>
 *        <td>Setup phase, coarsening. The "alpha" in the classical
 *            Ruge-St�ben setup phase. It controls the judgement
 *            whether a dependency is strong or weak: \f$ (i,j) \in S
 *            \Leftrightarrow  d(i,j) \ge \alpha \f$, with the dependency
 *            strength of point i from point j \f$ \displaystyle d(i,j)
 *            := \frac{|a_{ij}|}{\max\{|a_{ik}| : k \ne i\}} \f$ and matrix
 *            entries \f$ a_{ij} \f$. (used in
 *            Topology::CreateDependencyGraphs, see also class Topology).
 *            </td>
 *     </tr><tr>
 *        <td>AMG_StrongDiagRatio</td>
 *        <td>strongDiagRatio</td>
 *        <td>0.0</td>
 *        <td>Setup phase, coarsening. Additional necessary condition
 *            for a dependency to be judged as strong. The alpha
 *            (AMG_CoarseningAlpha) is related to the size of an
 *            offdiagonal entry compared to the other offdiagonal
 *            entries. But it must be also be large in comparison with
 *            the corresponding diagonal entry. AMG_StrongDiagRatio
 *            represents this ratio: \f$ (i,j) \in S \Rightarrow
 *            (AMG\_StrongDiagRatio \cdot a_{ii}) \le |a_{ij}| \f$.
 *            (used in Topology::CreateDependencyGraphs, see also class
 *            Topology). </td>
 *     </tr><tr>
 *        <td>AMG_ForceFineRatio</td>
 *        <td>forceFineRatio</td>
 *        <td>1e-10</td>
 *        <td>Setup phase, coarsening. If a row of the system matrix
 *            is strongly diagonal dominant, it can be forced to get
 *            an F-point. AMG_ForceFineRatio is the limit value for
 *            the ratio \f$ \frac{\sum_j |a_{ij}|}{a_{ii}} \f$. If this
 *            ratio falls below AMG_ForceFineRatio, the point is forced
 *            immediately to get an F-point. This is in particular useful
 *            for example in systems arrising from a PDE discretization
 *            in which Dirichlet boundary conditions are treated using
 *            penalty terms. In this case all Dirichlet points are forced
 *            to get F-points. (used in Topology::CreateDependencyGraphs,
 *            parameter diag_dominance) </td>
 *     </tr><tr>
 *        <td>AMG_BadCoarseningRatio</td>
 *        <td>badCoarseningRatio</td>
 *        <td>0.7</td>
 *        <td>Setup phase, coarsening. See setting \c AMG_NumBadCoarsenings.
 *        </td>
 *     </tr><tr>
 *        <td>AMG_BoundaryCoarsening</td>
 *        <td>coarsenBoundariesFirst</td>
 *        <td>false</td>
 *        <td>Setup phase of \b parallel AMG. The normal, per process
 *            purely serial coarsening works fine, but gets inefficient
 *            on coarser levels or even completely sticks. If we
 *            coarsen the boundaries first, also coarsening \b across
 *            boundaries, this effect can be switched of. Unfortunately
 *            this coarsening is sightly slower and the resulting AMG
 *            hierarchy is worse in terms of convergence. Therefore
 *            this flags is set to \b false per default, but it is
 *            activated automatically on deeper levels, if the
 *            coarsening rate reaches the value of setting
 *            AMG_BadCoarseningRatio. Setting this flag to \b true
 *            would activate the boundary coarsening per default on
 *            all levels, which seems to make sense only for experimental
 *            purposes.
 *            </td>
 *     </tr><tr>
 *        <td>AMG_DirectSolver</td>
 *        <td>directSolver</td>
 *        <td>NOSOLVER</td>
 *        <td>Setup phase, solution phase. Selects the direct solver
 *            that will be used to solve the system on the coarsest
 *            level. Use the values of the enumeration SolverType for this
 *            setting. Supported values:<br>
 *            - \c NOSOLVER: the smoother is used (expect poor performance)
 *            - \c LAPACK_LU: Lapack LU decomposition (compile code with
 *              preprocessor flag USE_LAPACK).
 *            - \c PARDISO: Pardiso LU decomposition (compile code with
 *              preprocessor flag USE_PARDISO).
 *            - \c LDL_SOLVER: LDL decomposition implemented in OLAS </td>
 *     </tr>
 *     </tr><tr>
 *        <td>AMG_logging</td>
 *        <td>logging</td>
 *        <td>false</td>
 *        <td>Setup phase, solution phase. If set to \c true, AMG will
 *            print some logging output to (*cla) during setup and
 *            solution cycles.</td>
 *     </tr>
 *  </table>
 *  <br><br>
 *  For developing settings and preprocessor flags see files ppflags.hh
 *  and multiGrid.hh.
 */
template <typename T>
class AMGSolver
{
    public:

        //! constructor with optional setting parameter

        /*! Constructor for the AMGSolver. As optional argument a
         *  pointer to an object of type OLAS_Params may be passed.
         *  If done so, all AMG parameters will be extracted from this
         *  object. If a NULL pointer (== default) is passed, some
         *  internal default values will be used.
         *  \param params pointer to an object of type OLAS_Param
         */
        AMGSolver( OLAS_Params* params = NULL );
        ~AMGSolver();

        //! set up the AMG solver

        /*! Sets up the AMG solver.
         *  \param sys_matrix pointer to the system matrix of the problem
         *  \param aux_matrix pointer to an optional auxiliary matrix
         *         for the problem
         */
        bool Setup( const CRS_Matrix<T>* const sys_matrix,
                    const CRS_Matrix<T>* const aux_matrix = NULL );

        //! resets the object to the state after creation
        void Reset();

        //! one AMG cycle
        void Cycle( const Vector<T>& rhs,
                          Vector<T>& sol ) const;

        //! prints out some solver data
        std::ostream& Print( std::ostream& out ) const;

    protected:

        OLAS_Params *const  Parameters_;
        HierarchyLevel<T>  *Hierarchy_;

        bool prepared_,
             logging_;
};

/**********************************************************/
} // namespace CoupledField

#endif // OLAS_AMG_HH
