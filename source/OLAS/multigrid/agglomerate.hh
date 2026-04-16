#ifndef OLAS_AGGLOMERATE_HH
#define OLAS_AGGLOMERATE_HH

/**********************************************************/

#include  "OLAS/multigrid/ppflags.hh"

#include <unordered_map>
#include <iostream>
#include <fstream>

#include "MatVec/opdefs.hh"

namespace CoupledField {

//! Class Agglomerate keeps the dependencies of an AMG hierarchy level.

/*! class Agglomerate keeps the dependency graphs of an edge-AMG hierarchy level.
 *  This means for every node we have to define the edges, connected to this node
 *  and also further neighbours, as specified in variable xxxxxxxxx
 */
template <typename T>
class Agglomerate
{
    public:

        //! enumeration constants for the array <code>CoarseIndex_</code>.
        //! <b>Note</b> that <code>CoarseIndex_</code> does not contain the
        //! constant <code>COARSE</code> for coarse nodes, but their index
        //! > 0 in the coarse system. <code>COARSE</code> is more ore less
        //! a dummy constant, or usable in comparison like "<code>< COARSE
        //! </code>".
        enum { UNDEFINED = -3, DIRICHLET_FINE = -2, FINE = -1, COARSE };


        //! simple constructor
        Agglomerate();

        //! simple constructor followed by a call of CreateAgglomerates

        /*! This constructor calls CreateAgglomertes.*/
        Agglomerate( const CRS_Matrix<Double>& auxMat);

        //! destructor
        ~Agglomerate();

        //! creates the graphs for the agglomerates

        Integer CreateAgglomerateGraphs( const CRS_Matrix<Double>& auxMat);

        inline void InsertEdgeList(const StdVector< StdVector< Integer> >& edgeIndeNode);

        inline void InsertNodeIndex(const StdVector< Integer>& nodeNumIndex);

        //! Define the coarse edges, based on the previously defined agglomerates
        UInt CreateCoarseEdges(const CRS_Matrix<Double>& coarseAuxMat);

        //! returns the whole agglomerate-construct
        inline const StdVector< StdVector< Integer> > GetAgglomerates()  const;
        //! returns the size of the fine system
        inline UInt GetSizeh()  const;
        //! returns the size of the coarse system
        inline UInt GetSizeH()  const;
        //! returns a pointer to <code>cInd_</code>
        inline const StdVector<Integer> GetCoarseFineSplitting()  const;
        //! returns true, if point [i] is an F-point (splitting must be present!)
        inline bool IsFPoint( const int i ) const;
        //! returns true, if point [i] is a C-point (splitting must be present!)
        inline bool IsCPoint( const int i )  const;
        //! returns the coarse index of point \c [i]
        inline Integer GetCoarseIndex( const Integer i ) const ;
        //! actually the inverse of GetCoarseIndex
        inline Integer GetIndexOfCoarseNode( const Integer i ) const;
        //! returns the agglomerate, which belongs to node i
        inline StdVector<Integer> GetAgglomerateOfNode( const Integer i ) const ;
        //! returns the agglomerate-number, which belongs to node i
        inline UInt GetAgglomerateNumOfNode( const Integer i ) const ;
        // for every C-point, get the number of F-points
        Integer GetNumInterpolatedPoints() const;
        //! if index i is valid, return the real nodeNum of coarse node i
        inline Integer GetRealNodeNumOfCIndex( const Integer i ) const;
        //! number of coarse edges...size of coarse system matrix
        inline UInt GetNumCoarseEdges( ) const;
        inline UInt GetNumFineEdges() const;
        //! get edge i, consisting of real nodeNum i and j
        inline StdVector<Integer> GetCoarseEdgeReal( const int i ) const;
        //! get edge i, consisting of node Index i and j
        inline StdVector<Integer> GetCoarseEdgeInd( const int i ) const;

        inline StdVector<Integer> GetCNodeNumIndex() const;
        inline Integer GetIndexOfRealNode( UInt i ) const;

        //! resets the object to the state after creation with the standard constructor
        void Reset();

        //! for every coarse node get the real node-number
        void CreateCoarseRealNodeNum();

    protected:

        //! check the orientation of every coarse-edge
        void OrientUniqueCEdges();

        void FillIndexNodeNum();

        static bool compareEdge0(const StdVector< Integer >& a, const StdVector< Integer >& b);
        static bool compareEdge1(const StdVector< Integer >& a, const StdVector< Integer >& b);

        //! pointer to the row starts of the matrix
        const UInt *NhStartIndex_;
        //! pointer to column indices of the matrix
        const UInt *NhEdges_;
        //! index of vector is edge in system matrix and contains
        //! nodes of the edge (in correct order)
        StdVector< StdVector< Integer> > edgeIndNode_;
        //! stores the real node number (from edgeIndNode_) for
        //! every index, this means index i in the auxiliary matrix
        //! has the real nodeNumber nodeNumIndex_[i]
        StdVector<Integer> nodeNumIndex_;
        //! inverse of nodeNumIndex_, this means value of
        //! indexNodeNum_[i]  is the index in the system matrix
        std::unordered_map<UInt, UInt> indexNodeNum_;
        StdVector< StdVector< Integer> > Agglomerates_;
        //! index map from i_h to i_H
        StdVector<Integer> CoarseIndex_;
        //! get the index in CoarseIndex_ of coarse node i
        StdVector<Integer> IndexOfCoarse_;
        //! store agglomerate-number for every node
        StdVector<UInt> nodeAgg_;
        //! store the real nodeNum for every coarse index
        StdVector<Integer> coarseRealNN_;
        //! number of points in the fine system
        UInt  Size_h_;
        UInt numCEdges_;
        //! number of points in the coarse system
        UInt  Size_H_;
        //! maximum agglomerate size
        UInt maxNumA_;
        //! minimal agglomerate size
        UInt minNumA_;
        //! coarse edges, from node i to node j...real nodeNums (correctly oriented)
        StdVector< StdVector< Integer> > cEdgesReal_;
        //! coarse edges, from node i to node j... indices in aux (correctly oriented)
        StdVector< StdVector< Integer> > cEdgesInd_;

        StdVector< StdVector< Integer> > coincidences_;


};

/**********************************************************/
} // namespace CoupledField

#endif // OLAS_AGGLOMERATE_HH
