#ifndef OLAS_SLOAN_HH
#define OLAS_SLOAN_HH

#include "OLAS/graph/BaseOrdering.hh"

namespace CoupledField {

  //! Sloan ordering algorithm for a matrix graph

  //! This class implements the %Sloan re-ordering algorithm for a matrix
  //! graph. The latter belongs to the class of profile minimisation
  //! approaches. For a detailed description see:
  //! <em>S.W. %Sloan, A Fortran Program for Profile and Wavefront Reduction,
  //! IJNME, Vol. 28, pp. 2651-2679, 1989</em>
  class Sloan : public BaseOrdering {

  public:

    //! Constructor

    //! This is the parameterised constructor of the class
    //! \param agraph pointer to original graph of matrix
    //! \param order contains for global node i the mapped index
    //! \param asize  number of elements (nodes) in the graph
    Sloan(StdVector<std::vector<UInt>>& agraph, StdVector<UInt>& order );

    //! Destructor
    virtual ~Sloan();

    //! renumber the nodes for minimum profile: main method
    void LabelGraph();

    //! compute a pair of nodes with maximum distance in the graph
    /*! goal: maximum on levels and on each level minimum on nodes belonging
      to this level
    */
    /*! \param ls list of nodes in rooted level structure
      \param xls indes to ls
      \param hlevel working array
      \param snode starting node, at which the labeling will start
      \param nc number of nodes in rooted level structure 
    */
    void PseudoDiameter(Integer *ls, Integer *xls, Integer *hlevel, 
                        Integer &snode, Integer &nc);

    //! computes rooted level structure for node = snode
    /*! \param snode starting node of rooted level structure (also called
      root)
      \param maxwidth maximum allowed number of nodes belonging to one level
      \param ls containes the rooted level structure
      \param xls index array to ls
      \param sdepth number of levels in rooted level structure
      \param width maximum number of nodes belonging to one level
    */
    void RootedLevel(Integer snode, Integer maxwidth, Integer *ls,
                     Integer *xls, Integer &sdepth, Integer &width);

    //! performs the renumbering
    /*! \param nc number of nodes to be reordered
      \param snode node at which numbering starts
      \param lstnum count of nodes, which have already be numbered
      \param q working array
      \param p working array
    */
    void NumberNodes(Integer nc, Integer snode, Integer &lstnum, Integer* q, 
                     Integer* p);

    //! orders a list in ascending sequence of their keys (using insertion
    //! sort)
    /*! \param nl length of the list
      \param list list to be sorted in ascending sequence of key
      \param key contains the keys
    */
    void SortList(Integer nl, Integer* list, Integer *key);

    //!
    void CalcProfile();

    //! Computes the total profile of the matrix graph

    //! This method can be used to compute both profiles, that of the old
    //! graph and that of the re-ordered graph. We return the profile sizes
    //! as floating point values to avoid problems with limited size of
    //! integral types.
    void GetProfile( Double &oldprof, Double &newprof );

    //! @{ \name Attributes for profile computation
    //! The attribute are related to computing and storing the profile
    //! of the matrix graph before and after reordering. Since the size
    //! of the profile (especially the un-reordered one) may easily exceed
    //! UINT_MAX, i.e. the maximal value of the standard unsigned
    //! integral data type in OLAS, we store the profiles as tuples
    //! (Mult, Rem). The actual profile size is then given by
    //! Mult * UINT_MAX + Rem.
    UInt profOldRem_;      //!< profile before renumbering
    UInt profNewRem_;      //!< profile after renumbering
    UInt profOldMult_;     //!< profile before renumbering
    UInt profNewMult_;     //!< profile after renumbering
    //@}
  };

} // namespace

#endif // OLAS_Sloan_HH
