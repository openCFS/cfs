// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id$ */

#ifndef OLAS_AMG_HH
#define OLAS_AMG_HH
/**********************************************************/
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "MatVec/CRS_Matrix.hh"

#include "OLAS/multigrid/hierarchylevel.hh"

#ifndef  LOGLINE
#define  LOGLINE  " --------------------------------------\n"
#endif

namespace CoupledField {
/**********************************************************/

//! AMG solver or preconditioner

/*! This class encapsulates a AMG matrix hierarchy and other needed
 *  stuff for an AMG used as solver or preconditioner. For usage create
 *  an AMGSolver object, call SetupMG(), then Cycle().
 *  <br><br>
 *  <b>Settings:</b><br>
 *  There are some settings, the behaviour of the AMG components can be
 *  be controlled with. These parameters are passed to the constructor
 *  as a pointer to an PtrParamNode object. If a NULL pointer is passed,
 *  some default values are used (see below). The settings, that are
 *  specified in this object are copied into a temporary object of type
 *  HierarchyLevel::Settings, which is passed during the setup phase
 *  through all hierarchy levels. This copy prevents the frequent usage
 *  of the GetValue routines of PtrParamNode, that might involve
 *  performance loss in heavy use.
 *  <br><br>
 *  <b>List of Settings</b><br>
 *  <table>
 *     <tr>
 *        <td>maxCoarseDepend</td>
 *        <td>6</td>
 *        <td>Setup phase, Ruge-Stueben C-F-splitting. The maximal number of coarse
 *            neighbours, an F-point can get due to their strong influence
 *            on that point. This does not mean, that this will be the
 *            number of C-neighbours, that will maximally surround an
 *            F-point, but the number of C-neighbours, that will cause
 *            an F-point to be judged as <i>well interpolated</i>. Thus
 *            a further neighbour will not be made a C-point due to its
 *            dependency to the named point, but it could still get C-point
 *            because of other reasons. </td>
 *     </tr><tr>
 *        <td>minSystemSize</td>
 *        <td>200</td>
 *        <td>Setup phase, C-F-splitting. The minimal size for the
 *            coarsest system. If the coarsening process yields a number
 *            of C-points for the next coarser system that falls below
 *            this limit, the coarsening result is rejected and the
 *            corresponding hierarchy level is not created. The present
 *            level gets the coarsest.</td>
 *     </tr><tr>
 *        <td>numPreSmoothings</td>
 *        <td>3</td>
 *        <td>Solution Phase. Number of pre-smoothing steps in the AMG
 *            solution cycle. Smoother type is chosen via parameter SmootherType </td>
 *     </tr><tr>
 *        <td>numPostSmoothings</td>
 *        <td>3</td>
 *        <td>Solution Phase. Number of post-smoothing steps in the AMG
 *            solution cycle. Smoother type is chosen via parameter SmootherType </td>
 *     </tr><tr>
 *        <td>CycleParameter</td>
 *        <td>1</td>
 *        <td>Solution phase. The recursion parameter for the solution
 *            cycle, which is the number of coarse corrections, applied
 *            <b>on each level except the top and the bottom level</b>.
 *            A value of 1 for example causes a common V-cycle, 2 a
 *            common W-cycle. Note that excluding the top level from
 *            this multiple coarse grid corrections causes the parameter
 *            value 2 to generate a real W-cycle, not a WW-cycle. </td>
 *     </tr><tr>
 *        <td>numBadCoarsenings</td>
 *        <td>5</td>
 *        <td>Setup phase, coarsening The recursive creation of coarser
 *            levels is interrupted, if the coarsening gets slow. "Getting
 *            slow" is defined as "the ratio of coarse to fine grid size
 *            was greater than \c numBadCoarsenings in
 *            \c numBadCoarsenings + 1 contiguous coarsenings".</td>
 *     </tr><tr>
 *        <td>alpha</td>
 *        <td>0.1</td>
 *        <td>Setup phase, coarsening. The "alpha" in the classical
 *            Ruge-Stueben setup phase. It controls the judgement
 *            whether a dependency is strong or weak: \f$ (i,j) \in S
 *            \Leftrightarrow  d(i,j) \ge \alpha \f$, with the dependency
 *            strength of point i from point j \f$ \displaystyle d(i,j)
 *            := \frac{|a_{ij}|}{\max\{|a_{ik}| : k \ne i\}} \f$ and matrix
 *            entries \f$ a_{ij} \f$.
 *            </td>
 *     </tr><tr>
 *        <td>strongDiagRatio</td>
 *        <td>0.0</td>
 *        <td>Setup phase, coarsening. Additional necessary condition
 *            for a dependency to be judged as strong. The alpha
 *            (AMG_CoarseningAlpha) is related to the size of an
 *            offdiagonal entry compared to the other offdiagonal
 *            entries. But it must be also be large in comparison with
 *            the corresponding diagonal entry. AMG_StrongDiagRatio
 *            represents this ratio: \f$ (i,j) \in S \Rightarrow
 *            (AMG\_StrongDiagRatio \cdot a_{ii}) \le |a_{ij}| \f$. </td>
 *     </tr><tr>
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
 *            to get F-points. </td>
 *     </tr><tr>
 *        <td>directSolver</td>
 *        <td>BaseSolver::PARDISO_SOLVER</td>
 *        <td>Setup phase, solution phase. Selects the direct solver
 *            that will be used to solve the system on the coarsest
 *            level. Currently only PARDISO is available </td>
 *     </tr>
 *  </table>
 *  <br><br>
 */
template <typename T>
class AMGSolver
{
public:

  //! constructor with optional setting parameter

  /*! Constructor for the AMGSolver. As optional argument a
   *  pointer to an object of type PtrParamNode may be passed.
   *  If done so, all AMG parameters will be extracted from this
   *  object. If a NULL pointer (== default) is passed, some
   *  internal default values will be used.
   *  \param params pointer to an object of type PtrParamNode
   */
  AMGSolver( PtrParamNode& xml, PtrParamNode& olasInfo  )
    : Hierarchy_( NULL ),
      prepared_( false ),
      logging_( true )
    {
      Parameters_ = xml;
      Info_ = olasInfo;

    };

  ~AMGSolver(){
    Reset();
  };

  struct EdgeGeom{
    StdVector<Integer> eNodes; // edge nodes
    Double length;
  };

  //! set up the AMG solver

  /*! Sets up the AMG solver.
   *  \param sys_matrix pointer to the system matrix of the problem
   *  \param aux_matrix pointer to the auxiliary matrix
   */
  bool Setup( const CRS_Matrix<T>* const sys_matrix,
              const CRS_Matrix<Double>* const aux_matrix,
              AMGType amgType,
              const StdVector< StdVector< Integer> >& edgeIndNode,
              const StdVector<Integer>& nodeNumIndex){

    std::cout << "++ Performing Setup of AMG-preconditioner "<<std::endl;

    // reset the whole object
    Reset();


    Hierarchy_ = new HierarchyLevel<T>();

    // set logging state
    Hierarchy_->SetLogging( logging_ );
    // insert the matrices into the top level
    Hierarchy_->InsertSysMatrix( sys_matrix );
    Hierarchy_->InsertAuxMatrix( aux_matrix );
    Hierarchy_->SetAMGtype(amgType);
    if(amgType == EDGE){
      Hierarchy_->SetEdgeIndNode(edgeIndNode);
      Hierarchy_->SetNodeNumIndex(nodeNumIndex);
    }

    // create a settings object for the hierarchy
    typename HierarchyLevel<T>::Settings settings( Parameters_, Info_ );

    // setup
    if(amgType != EDGE){
      if( false == Hierarchy_->Setup(&settings) ) {
        Reset();
        return false;
      }
    }else{
      if( false == Hierarchy_->SetupEdge(&settings) ) {
              Reset();
              return false;
            }
    }
    return this->prepared_ = true;
  };


  //! resets the object to the state after creation
  void Reset(){

    // Remove system and auxiliary matrix from top level
    // of hierarchy, since they are not intended to be owned
    // by this object.
    if( Hierarchy_ ) {
      Hierarchy_->UnhookSysMatrix();
      Hierarchy_->UnhookAuxMatrix();
    }

    delete Hierarchy_;
    Hierarchy_ = NULL;

    this->prepared_ = false;
  };

  //! one AMG cycle
  void Cycle( const Vector<T>& rhs,
      Vector<T>& sol ){

              if( !prepared_ ) {
                  WARN( "called AMGSolver::Cycle at a non prepared AMG" );
                  return;
              }

              typename HierarchyLevel<T>::Settings settings( Parameters_, Info_ );

              //if( logging_ ) (*cla)<<" AMG: starting cycle\n";
              // We need this dirty cast to a non-constant HierarchyLevel,
              // since the cycle will change the HierarchyLevel's temporary
              // vectors.
              static_cast<HierarchyLevel<T>*>(Hierarchy_)->Cycle( rhs, sol,
                                   settings.numPreSmoothings,
                                   settings.numPostSmoothings,
                                   settings.CycleParameter );

              //if( logging_ )  (*cla) <<" AMG: finished cycle\n";
  };

  //! prints out some solver data
  std::ostream& Print( std::ostream& out ){
    out << LOGLINE" AMG solver:"<< std::endl;

    if( prepared_ ) {
        Hierarchy_->Print( out );
    } else {
        out << " AMG is not prepared" << std::endl;
    }

    out << LOGLINE;
    return out;
  };


protected:

  PtrParamNode  Parameters_;
  PtrParamNode Info_;
  HierarchyLevel<T>  *Hierarchy_;

  bool prepared_, logging_;
};

/**********************************************************/
} // namespace CoupledField


#endif // OLAS_AMG_HH
