// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id$ */

#ifndef OLAS_HIERARCHYLEVEL_HH
#define OLAS_HIERARCHYLEVEL_HH

/**********************************************************/

#include "matvec/matvec.hh"
#include "multigrid/topology.hh"
#include "multigrid/smoother.hh"
#include "multigrid/transfer.hh"

namespace CoupledField {
/**********************************************************/
class BaseSolver;

//! base class for one level in a multigrid hierarchy

/*! This class encapsulates all data needed on one level
 *  in a multigrid hierarchy, e.g.
 *   - system matrix and eventual auxiliary matrix
 *   - transfer operator
 *   - smoother
 *   - some vectors, that will be needed temporarily and are created
 *     for performance reasons in the setup phase
 *   - the topology, containing graphs of strong dependencies and the
 *     coarse-fine splitting of the level
 *  To create a AMG hierarchy, create an object of type HierarchyLevel,
 *  insert a matrix and call HierarchyLevel::Setup. Setup will prepare
 *  the whole AMG scheme for this level and produce the underlying levels
 *  recusively until the minimal system size (-> class AMGSolver) is
 *  reached.
 */
template <typename T>
class HierarchyLevel
{
    public:

        //! entry type of the matrices (e.g. tiny matrices)
        typedef typename AssocType<T>::T_Mtype T_Mtype;
        //! entry type of the vectors (e.g. tiny vectors)
        typedef typename AssocType<T>::T_Vtype T_Vtype;
        //! scalar type (e.g. double, even if T_Mtype is a tiny matrix)
        typedef typename AssocType<T>::T_Stype T_Stype;

        //! class for passing settings through all levels

        /*! The class OLAS_Params offers a comfortable way to store settings
         *  mapped to identification strings. But inside the routings the
         *  parameters are needed in it might be inefficient to get them each
         *  time from a demand to OLAS_Params. In particular parameters must
         *  be passed between the AMG levels, since each level works
         *  autonomiously. A pointer to an object of type class Settings
         *  can be passed to the top level and is then passed to all the
         *  underlying ones. So each parameter can be seen on each level
         *  without copying the parameters or passing them to the whole bunch
         *  of functions that are called recursively through the levels.
         */
        class Settings { public:
            Settings( OLAS_Params* const params = NULL );
            std::ostream& Print( std::ostream& out ) const;
            Integer    maxCoarseDepend,  //!< maximal number of interpolation points
                       minSystemSize,    //!< minimal size of coarsest system
                       numPreSmoothings, //!< number of presmoothing steps
                       numPostSmoothings,//!< number of postsmoothing steps
                       CycleParameter,   //!< recursion parameter for the cycle type
                       numBadCoarsenings, //!< number of tolerated bad coarsenings
                       gamma;
            Double     alpha,            //!< coarsening threshold
                       strongDiagRatio,  //!< second coarsening parameter
                       forceFineRatio,   //!< forces F-points, if diagonal entry dominates
                       badCoarseningRatio; //!< defines a "bad coarsening"
            bool       keepTopology,     //!< flag to keep topology after setup
                       logging;          //!< flag for logging into (*cla)
            SolverType directSolver;     //!< direct solver on the coarsest level
            AMGInterpolationType InterpolationType; //!< type of interpolation
            AMGSmootherType      SmootherType; //!< smoother type
            Double     smootherOmega; //!< damping parameter for smoother
        };

        //! constructor
        HierarchyLevel();
        //! destructor
        ~HierarchyLevel();

        //! returns a pointer to the next coarser level
        inline const HierarchyLevel<T>* GetNextLevel() const;
        //! returns the level ID
        inline Integer GetLevelID() const;

        //! Resets the hole object to the state after a constructor call
        void Reset();

        //! inserts the system matrix

        /*! Inserts the system matrix into the hierarchy level object.
         *  An already inserted matrix will be destroyed. After this
         *  call the passed matrix object belongs to this class, i.e.
         *  it will be destroyed in the destructor or in the a call of
         *  Reset(). To regain the ownership of the matrix call
         *  UnhookSysMatrix().
         */
        void InsertSysMatrix( const CRS_Matrix<T>* const sysmat );

        //! inserts the auxiliary matrix
        
        /*! Inserts the auxiliary matrix into the hierarchy level object.
         *  An already inserted matrix will be destroyed. After this
         *  call the passed matrix object belongs to this class, i.e.
         *  it will be destroyed in the destructor or in the a call of
         *  Reset(). To regain the ownership of the matrix call
         *  UnhookSysMatrix().
         */
        void InsertAuxMatrix( const CRS_Matrix<T>* const auxmat );

        //! returns a pointer to the system matrix

        /*! Returns a pointer to the system matrix of the AMG level.
         *  <b>Note</b> that the CRS_Matrix inside a AMG level might
         *  not conform the layout conventions for CRS matrices, which
         *  fixes the diagonal element beeing placed at leading position
         *  and the offdiagonal elements beeing sorted with ascending
         *  column indices in each row. The matrices on an AMG level
         *  might not conform the second property, but it is granted
         *  that the diagonal element is the first one in each row.
         */
        inline const CRS_Matrix<T> *GetSysMatrixPtr() const;

        //! returns a pointer to the auxiliary matrix

        /*! Returns a pointer to the auxiliary matrix of the AMG level.
         *  <b>Note</b> that the CRS_Matrix inside a AMG level might
         *  not conform the layout conventions for CRS matrices, which
         *  fixes the diagonal element beeing placed at leading position
         *  and the offdiagonal elements beeing sorted with ascending
         *  column indices in each row. The matrices on an AMG level
         *  might not conform the second property, but it is granted
         *  that the diagonal element is the first one in each row.
         */
        inline const CRS_Matrix<T> *GetAuxMatrixPtr() const;

        //! unhooks the system matrix

        /*! Unhooks the system matrix from the hierarchy level object.
         *  As a matrix, that was inserted with InsertSysMatrix(..),
         *  is owned then by the hierarchy level object, it is also
         *  destroyed in the destructor or in function Reset(). This
         *  function returns a pointer to the matrix and deprives
         *  the hierarchy object the control over the matrix.
         */
        CRS_Matrix<T> *UnhookSysMatrix();

        //! unhooks the auxiliary matrix

        /*! Unhooks the auxiliary matrix from the hierarchy level object.
         *  As a matrix, that was inserted with InsertAuxMatrix(..),
         *  is owned then by the hierarchy level object, it is also
         *  destroyed in the destructor or in function Reset(). This
         *  function returns a pointer to the matrix and deprives
         *  the hierarchy object the control over the matrix.
         */
        CRS_Matrix<T> *UnhookAuxMatrix();

        //! setup for this level

        /*! Setup for this AMG-Level. Also the underlying, coarser
         *  levels will be set up recursively.
         *  \param settings pointer to an object of type
         *         HierarchyLevel<T>::Settings
         *  \param numBadCoarsenings number of bad coarsenings on
         *         higher levels (see class AMGSolver, setting
         *         \c AMG_NumBadCoarsenings.
         */        
        bool Setup( Settings* const settings,
                    Integer         numBadCoarsenings = 0 );

        //! one AMG cycle
        void Cycle( const Vector<T>& rhs,
                          Vector<T>& sol,
                    const Integer    num_presmoothings  = 1,
                    const Integer    num_postsmoothings = 1,
                    const Integer    cycle_parameter    = 1  );

        //! returns the operator complexity \f$\sum_{i=0}^N\frac{nnz(A_i}{A}\f$
        Double GetOperatorComplexity() const;

        //! switches the logging on and off
        void SetLogging( const bool logging = true ) {
            logging_ = logging;
#ifdef PROFILE_HIERARCHYLEVEL
            if( !logging && !GetLevelID() ) {
                Warning( "AMG: PROFILING implies logging, which "
                         "is explicitly switched off.\n      "
                         "Setting ignored and reactivated logging",
                         __FILE__, __LINE__ );
                
            }
            logging_ = true;
#endif
        }

        //! prints this and all underlying levels
        std::ostream& Print( std::ostream& out,
                             const bool    color = false ) const;

    protected:

        //! protected: only other levels can set the ID
        HierarchyLevel( const Integer level_id );

        //! creates and prepares the topology object
        bool SetupTopology( const Double alpha,
                            const Double strongDiagRatio,
                            const Double forceFineRatio,
                            bool *const  penaltyFlags = NULL );

        //! creates and prepares the smoother
        bool SetupSmoother( const Settings* const settings,
                            bool *const penaltyFlags = NULL );

        //! creates and prepares the direct solver (on the coarsest level)
        bool SetupDirectSolver( const Settings *const settings );

        //! creates the temporary data (e.g. vectors, ..)
        bool SetupAuxData( const Integer sizeh, const Integer sizeH );

#ifdef AMG_EXPORT_COMPONENTS
        void ExportComponents();
#endif

#ifdef AMG_TRACE_CF_SPLITTING
        bool TraceCFSplitting( AMG_CF_TRACE_TYPE *trace    = NULL,
                               Integer           *indexMap = NULL ) const;
#endif

        //! level number in the hierarchy
        const Integer        LevelID_;
        //! system matrix
        CRS_Matrix<T>       *SysMatrix_,
        //! auxiliary matrix
                            *AuxMatrix_;
        //! system matrix in different format (e.g. for the direct solver)
        BaseMatrix          *DirSysMatrix_;
        //! pattern pool (used, if \c DirSysMatrix_ is a SCRS_Matrix)
        PatternPool         *dirSysMatPPool_;
        //! topology
        Topology<T>         *Topology_;
        //! transfer operator
        TransferOperator<T> *Transfer_;
        //! smoother
        Smoother<T>         *Smoother_;
        //! direct solver for the coarses level
        BaseSolver          *directSolver_;
        //! Pointer to parameter object for direct solver
        OLAS_Params *directSolverParams_;

        //! temporary vector of size nh
        Vector<T>           *vh1_,
        //! first temporary vector of size nH
                            *vH1_,
        //! second temporary vector of size nH
                            *vH2_;            

        //! pointer to the next coarser level
        HierarchyLevel<T>   *nextLevel_;

        //! if true, the class gets verbose to (*cla)
        bool                 logging_,
                             deleteDirSysMatrix_;
};

/**********************************************************/
} // namespace CoupledField

namespace std {
template <typename T> std::ostream& operator <<
( std::ostream& out, const OLAS::HierarchyLevel<T>& level );
}

/**********************************************************/

#endif // OLAS_HIERARCHYLEVEL_HH
