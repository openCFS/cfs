#ifndef OLAS_TOPOLOGY_HH
#define OLAS_TOPOLOGY_HH

/**********************************************************/

#include  "multigrid/ppflags.hh"
#include  "multigrid/depgraph.hh"

namespace CoupledField {
/**********************************************************/
#define  TOPOLOGY_COARSE_GRAPH_STRETCH_FACTOR  5
/**********************************************************/
// Basically in method GetNextCoarsePoint method GetCoarseImportance
// is invoked redundantly for some points. This can be avoided by
// keeping track of the points, for which the coarse importance has
// already been calculated within one call of GetNextCoarsePoint.
// This strategy reduces the time of the coarsening process by a
// factor of about 0.8, but requires the additional memory for the
// tracking array.
// The flag TOPOLOGY_AVOID_REDUNDANT_IMPORTANCE_CALCULATION activates
// this feature.
#define  TOPOLOGY_AVOID_REDUNDANT_IMPORTANCE_CALCULATION
/**********************************************************/

//! class Topology keeps the dependencies of an AMG hierarchy level.

/*! class Topology keeps the dependency graphs \f$ S \f$ and \f$ S^T\f$
 *  of an AMG hierarchy level.
 *  <br>
 *  \f$ (i,j) \in S \Leftrightarrow  d(i,j) \ge \alpha \f$ <br> with
 *  the dependency strength of point i from point j
 *  <br>
 *  \f$ \displaystyle d(i,j) := \frac{|a_{ij}|}{max\{|a_{ik}| : k \ne
 *  i\}} \f$
 *  <br>
 *  and matrix entries \f$ a_{ij} \f$.
 *  <br>
 *  \f$ (i,j) \in S^T \Leftrightarrow  d(j,i) \in S \f$ <br>
 */
template <typename T>
class Topology
{
    public:
        //! enumeration constants for the array <code>CoarseIndex_</code>.
        //! <b>Note</b> that <code>CoarseIndex_</code> does not contain the
        //! constant <code>COARSE</code> for coarse nodes, but their index
        //! > 0 in the coarse system. <code>COARSE</code> is more ore less
        //! a dummy constant, or usable in comparison like "<code>< COARSE
        //! </code>".
        enum { DIRICHLET_FINE = -2, FINE = -1, UNDEFINED = 0, COARSE };

        //! simple constructor
        Topology();

        //! simple constructor followed by a call of CreateDependencyGraphs

        /*! This constructor calls CreateDependencyGraphs. For more
         *  details about the parameters see the description of method
         *  Topology::CreateDependencyGraphs.
         */
        Topology( const CRS_Matrix<T>& matrix,
                  const Double         alpha,
                  const Double         tolerance,
                  const Double         diag_dominance );
        //! destructor
        ~Topology();

        //! returns the dependency graph \f$ S \f$
        inline const DependencyGraph<T>& GetS() const;
        //! returns the dependency graph \f$ S^T \f$
        inline const DependencyGraph<T>& GetST() const;
        //! returns the size of the fine system
        inline Integer GetSizeh() const;
        //! returns the size of the coarse system
        inline Integer GetSizeH() const;
        //! returns a pointer to <code>CoarseIndex_</code>
        inline const Integer* GetCoarseFineSplitting() const;
        //! returns true, if point [i] is an F-point (splitting must be present!)
        inline bool IsFPoint( const int i ) const;
        //! returns true, if point [i] is a C-point (splitting must be present!)
        inline bool IsCPoint( const int i ) const;
        //! returns the coarse index of point \c [i]
        inline Integer GetCoarseIndex( const Integer i ) const;
        
        

        //! creates the graphs of strong dependencies

        /*! Creates the graphs \f$ S \f$ and \f$ S^T \f$ of strong
         *  dependencies. In addition to this the first candidate (node)
         *  for the coarse system is evaluated and explicit Dirichlet
         *  nodes are processed. 
         *  Returns the index of the first coarse node.
         *  \param matrix system matrix
         *  \param alpha the coarsening parameter \f$ \alpha \f$ in the
         *         Ruge-St�ben algorithm for the evaluation of strong
         *         dependencies
         *  \param tolerance additional evaluation criterion for strong
         *         dependencies: \f$ (i,j) \in S \Rightarrow (tolerance
         *         \cdot a_{ii}) \le |a_{ij}| \f$.
         *  \param diag_dominance a value \f$\gamma\f$ < 1 that forces
         *         point \f$i\f$ to get an F-points, if \f$ \gamma \cdot
         *         a_{ii} > \sum_j |a_{ij}| \f$. This is usefull in,
         *         particular, to force F-points, if Dirichlet boundary
         *         conditions are treated via penalty factors.
         *  \return index of the first node, that should be set as coarse
         *          node in a coarsening routine. This will be the point
         *          with maximal number of other points, that strongly
         *          depend on this point.
         */
        Integer CreateDependencyGraphs( const CRS_Matrix<T>& matrix,
                                        const Double         alpha,
                                        const Double         tolerance,
                                        const Double         diag_dominance
#ifdef AMG_DIRICHLET_MIXED_SMOOTHING
                                      , bool *const dirichlet_flags = NULL
#endif
                                       );

        //! calculates the C-F-splitting and returns the number of coarse points
        Integer CalcCoarseFineSplitting( const Integer max_dependency,
                                         const Integer gamma );
        
        //! returns \f$ |N_i \cap F| \f$ == # of fine grid neighbours
        inline Integer GetNumFineNeighbours( const Integer i ) const;
        
        //! calculates the graphs, needed for the Galerkin product

        /*!
         */
        void CalcGalerkinGraphs( DependencyGraph<T>& graph_AHT,
                                 DependencyGraph<T>& graph_VT,
                                 const Integer stretch_factor =
                                       TOPOLOGY_COARSE_GRAPH_STRETCH_FACTOR ) const;
        
        //! <code>sizes[i]</code> <- # nodes, node[i] is interpolated from

        /*! Writes at <code>sizes[i]</code> the number of nodes, node [i]
         *  is interpolated from: <br>
         *  <code>sizes[i]</code> = \f$ i \in C \: ? \: 1 \: :
         *  |S_i \cap C| \f$ <br>
         *  So an already present C-F-splitting (in <code>CoarseIndex_
         *  </code>) is prerequisite for the call of this function. If you
         *  need only the return value, call the overloaded function version
         *  without a parameter.
         *  \param sizes <code>Integer</code> array of size \f$ n_h \f$,
         *               one-based indexed (use <code>NEWARRAY()</code>)
         *  \return total number of points needed for interpolation ( ==
         *          \f$ \sum_{i=1}^n size[i] \f$ )
         */
        Integer GetNumInterpolatingPoints( Integer* const sizes ) const;

        //! returns the number of points, point [i] is interpolated from
        Integer GetNumInterpolatingPoints( const Integer i ) const;
        
        //! returns the number of points, needed for interpolation

        /*! Returns the total number of points, needed for interpolation.
         *  This, for example, can be used as number of non-zero entries
         *  in the prolongation or restriction operator. This function
         *  returns the same value as GetNumInterpolatingPoints(Integer* const).
         */
        Integer GetNumInterpolatingPoints() const;

        //! Returns \f$ |S^T_i| + |S^T_i \cap F| \f$ (also in setup phase)
        Integer GetCoarseImportance( const Integer i ) const;

        //! <code>sizes[i]</code> <- # points, that are interpolated from point[i]

        /*! This function writes for each C-point, the number of F-points,
         *  that are interpolated from it, into the passed array \c sizes.
         *  Therefor \c sizes is addressed with \b coarse indices. So it
         *  must have a lenght equal to the number of coarse points.
         */
        Integer GetNumInterpolatedPoints( Integer* const sizes ) const;

        
        //! writes strongly influencing C-neigbours into the array <code>neighbours</code>

        /*! Writes all points \f$ i \in C \f$ with \f$ (p,i) \in S\f$
         *  into the passed array <code>neighbours</code>. The first
         *  is written at <code>neigbours[0]</code>, the number of written
         *  points is returned. <code>neighbours</code> must be large
         *  enough to keep all entries.
         *  \param p the point you want the neighbours of
         *  \param neighbours array, into that the neigbours are written
         *  \return number of written neighbours
         */
        Integer WriteStrongCoarseNeighbours( const Integer        p,
                                                   Integer* const neighbours ) const;

        //!
        Integer WriteCoarseNeighbours( const Integer        p,
                                             Integer* const c_neighbours ) const;

        //! writes all fine grid neigbours into the array <code>neighbours</code>

        /*! Writes all points \f$ i \in F \f$ with \f$ a_{pi} \ne 0 \f$
         *  into the passed array <code>neighbours</code>. The first
         *  is written at <code>neigbours[0]</code>, the number of written
         *  points is returned. <code>neighbours</code> must be large
         *  enough to keep all entries.
         *  \param p the point you want the neighbours of
         *  \param f_neighbours array, into that the neigbours are written
         *  \return number of written neighbours
         */
        Integer WriteFineNeighbours( const Integer        p,
                                           Integer* const f_neighbours ) const;

        //! returns the next point, that can be put into C

        /*! Returns the next point, that can be put into C, starting
         *  the search for this next C-point at the passed point p.
         *  This start point p is assumed to be a C-point, too. Therefore
         *  the algorithm skips one layer of dependencies, i.e. N(p) and
         *  checks all points j in N(p). The value of each j in N(p) for
         *  getting a C-point is calculated as \f$ v_j := |N_j \cup F| +
         *  |S^T_j| \f$. If \f$ v_j \ge \f$ max_dependency, j is
         *  immediately taken as new C-point. Otherwise the j with the
         *  maximal \f$ v_j \f$ is returned.
         */
        Integer GetNextCoarsePoint( const Integer p,
                                    const Integer max_dependency );

        //! sets node [i] as a coarse grid point, in short: [i] -> C

        /*! Makes node [i] a coarse grid node. This also implies that all
         *  points, that strongly depend on this point, can be put into F.
         *  So the function has the following effect: <br>
         *  \f$ \quad \quad i \in C \f$ and \f$ j \in F \quad \forall \:
         *  (i,j) \in S^T \f$
         *  \param i index of node, that should be put into C
         */
        inline void SetCoarsePoint( const Integer i );
        //! sets node [i] as a Dirichlet node (in particular into F)
        inline void SetDirichlet( const Integer i );

        //! resets the object to the state after creation with the standard constructor
        void Reset();

    protected:

#ifdef TOPOLOGY_IMPORT_CF_SPLITTING
        Integer ImportCFSplitting();
#endif
#ifdef TOPOLOGY_EXPORT_CF_SPLITTING
        void ExportCFSplitting();
#endif

        //! pointer to the row starts of the matrix
        const Integer *NhStartIndex_;
        //! pointer to column indices of the matrix
        const Integer *NhEdges_;

        //! graph \f$S_h\f$ of strong influences
        DependencyGraph<T>  S_;
        //! graph \f$S_h^T\f$ of strong dependencies
        DependencyGraph<T>  ST_;

        //! index map \f$i_h \mapsto i_H\f$
        Integer *CoarseIndex_;
        //! index of the next coarse node
        Integer  startCoarsePoint_;
        //! number of points in the fine system
        Integer  Size_h_;
        //! number of points in the coarse system
        Integer  Size_H_;

#ifdef TOPOLOGY_AVOID_REDUNDANT_IMPORTANCE_CALCULATION
        //! mrgl
        Integer *importanceKnown_;
#endif
};

/**********************************************************/
} // namespace CoupledField

#endif // OLAS_TOPOLOGY_HH
