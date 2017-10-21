#include "OLAS/multigrid/topology.hh"


namespace CoupledField {


template <typename T>
Topology<T>::Topology()
    : NhStartIndex_(NULL),
      NhEdges_(NULL),
      CoarseIndex_(0),
      startCoarsePoint_(-1),
      Size_h_(0),
      Size_H_(0)
#ifdef TOPOLOGY_AVOID_REDUNDANT_IMPORTANCE_CALCULATION
     ,importanceKnown_(0x0)
#endif
{
}

/**********************************************************/

template <typename T>
Topology<T>::Topology( const CRS_Matrix<T>& auxMat,
                       const Double         alpha,
                       const Double         tolerance,
                       const Double         diag_dominance )
    : NhStartIndex_(NULL),
      NhEdges_(NULL),
      CoarseIndex_(0),
      startCoarsePoint_(-1),
      Size_h_(0),
      Size_H_(0)
#ifdef TOPOLOGY_AVOID_REDUNDANT_IMPORTANCE_CALCULATION
     ,importanceKnown_(0x0)
#endif
{

    CreateDependencyGraphs( auxMat, alpha, tolerance, diag_dominance );
}

/**********************************************************/

template <typename T>
Topology<T>::~Topology()
{
    
    Reset();
}

/**********************************************************/

template <typename T>
inline const DependencyGraph<T>& Topology<T>::GetS() const
{
    return S_;
}

/**********************************************************/

template <typename T>
inline const DependencyGraph<T>& Topology<T>::GetST() const
{
    return ST_;
}

/**********************************************************/

template <typename T>
inline UInt Topology<T>::GetSizeh() const
{
    return Size_h_;
}

/**********************************************************/

template <typename T>
inline UInt Topology<T>::GetSizeH() const
{
    return Size_H_;
}

/**********************************************************/

template <typename T>
inline const Integer* Topology<T>::GetCoarseFineSplitting() const
{
    return CoarseIndex_;
}

/**********************************************************/

template <typename T>
inline bool Topology<T>::IsFPoint( const int i ) const
{
    return CoarseIndex_[i] == FINE;
}

/**********************************************************/

template <typename T>
inline bool Topology<T>::IsCPoint( const int i ) const
{
    return CoarseIndex_[i] >= COARSE;
}

/**********************************************************/

template <typename T>
inline Integer Topology<T>::GetCoarseIndex( const Integer i ) const
{
    return CoarseIndex_[i];
}


/**********************************************************/

template <typename T>
Integer Topology<T>::
CreateDependencyGraphs( const CRS_Matrix<T>& auxMat,
                        const Double         alpha,
                        const Double         tolerance,
                        const Double         diag_dominance)
{
    
    // loop indices, consistent for the whole function
    Integer i,  // i = row index in auxMat, in [1,n]
            ij; // ij = absolute index in pCol and pDat for an entry in row i

    // pointers for direct matrix access
    const UInt *const pRow = auxMat.GetRowPointer();
    const UInt *const pCol = auxMat.GetColPointer();
    const T *const pDat = auxMat.GetDataPointer();

    /////////////////////////////////
    // prepare some data structures
    /////////////////////////////////

    Size_h_ = auxMat.GetNumRows();
    
    // prepare dependency graph S
    S_.Create( auxMat,  // build it from matrix "auxMat"
               true,    // initialize it as empty graph
               false ); // use the matrix' start index array (do not copy)

    // prepare dependency graph ST
    // Eploiting the knowledge about the symmetric non-zero structure
    // of the matrix we generate also ST from the problem matrix. ST
    // will use the start-array of matrix, too.
    ST_.Create( auxMat,  // build it from matrix "auxMat"
                true,    // initialize it as empty graph
                false ); // use the matrix' start index array (do not copy)


    // allocate the array for the C-F-splitting and the coarse indices
    delete [] ( CoarseIndex_ );  CoarseIndex_  = NULL;
    NEWARRAY( CoarseIndex_, Integer, auxMat.GetNumRows() );
    for( i = 0; i < (Integer)auxMat.GetNumRows(); i++ ) {
        CoarseIndex_[i] = (Integer)UNDEFINED;
    }

    ///////////////////////////////
    // create the graphs S and ST
    ///////////////////////////////

    // inspect auxMat row by row
    for( i = 0; i < (Integer)auxMat.GetNumRows(); i++ ) {
        // get number of non-zero entries in row [i]
        UInt RowLength = auxMat.GetRowSize( i );
        if( RowLength == 0 ) {
            EXCEPTION( "topology.cc: auxMat is singular");
        }
        // get diagonal entry
        T		diag       = pDat[pRow[i]],
                maxOffdiag = 0, // maximal offdiagonal entry (abs. value)
                sumOffdiag = 0, // sum of offdiagonal entries (abs. values)
                Aij        = 0; // an arbitrary temporary matrix entry
        // get maximal offdiagonal entry and rowsum (without diagonal)
        for( ij = (Integer)pRow[i]+1; ij < (Integer)pRow[i+1]; ij++ ) {
            if( maxOffdiag < fabs(pDat[ij]) )  maxOffdiag = fabs(pDat[ij]);
            sumOffdiag += fabs(pDat[ij]);
        }
        // points with an extremely dominant diagonal entry are
        // immediately forced to get F-points, that are treated
        // like Dirichlet points
        if( sumOffdiag < diag_dominance * diag ) {
            CoarseIndex_[i] = (Integer)DIRICHLET_FINE;
        // otherwise continue the usual building of the graphs
        } else {
            // check all dependencies in row [i]
            for( ij = (Integer)pRow[i]/*+1*/; ij < (Integer)pRow[i+1]; ij++ ) {
                Aij = pDat[ij];
                // (i,j) in S, if A_ij sufficently large in relation to
                // diagonal entry A_ii AND i depends strongly on j
                if(    fabs(Aij) >= tolerance * diag
                    && fabs(Aij) >  alpha * maxOffdiag
                    && i != (Integer)pCol[ij]) {
                    // add (i,j) to S (no double insertion possible)
                    S_.AddEdge( i, pCol[ij] );
                    // add (j,i) to ST (no double insertion possible)
                    ST_.AddEdge( pCol[ij], i );
                    // LAS: In LAS ST_ has been just a copy of S_. A I
                    //      cannot understand the sense of this copy, it
                    //      is implemented differently here.
                }
            }
        }
    }

/*
    S_.Print(std::cout);
    ST_.Print(std::cout);
    std::ofstream fileoutS("logS_.txt");
    S_.Print(fileoutS);
    std::ofstream fileoutST("logST_.txt");
    ST_.Print(fileoutST);
*/
    ////////////////////////////////////////////////////////////
    // process explicit and implicit Dirichlet Values and
    // evaluate the first node that should become a coarse node
    ////////////////////////////////////////////////////////////

    // maximal number of strong dependencies for one node
    Integer maxNumStrongDep = 0;

    for( i = 0, startCoarsePoint_ = 0; i < (Integer)auxMat.GetNumRows(); i++ ) {
        // The first coarse node will be the one that influences the
        // most other nodes strongly, i.e. i with maximal |ST(i)|
        // LAS: LAS additionally (AND) checked, if the number of entries
        //      in row [i] of the matrix is != 0. I do not understand this
        //      condition yet, but I port it anyway.
        if( maxNumStrongDep < ST_.GetNumEdges(i) && auxMat.GetRowSize(i) ) {
            maxNumStrongDep = ST_.GetNumEdges(i);
            startCoarsePoint_ = i;
        }
        /// Nodes with an extremly dominant diagonal entry were marked
        /// in the loop above.
        if( CoarseIndex_[i] == DIRICHLET_FINE ||
        /// Nodes without any dependency from other nodes are treated as
        /// Dirichlet nodes and put in F
            auxMat.GetRowSize(i) == 1 ) {
            // LAS: LAS called BaseTopology::SetDirichlet(i) at this
            //      position. To rebuild this call in detail, we would
            //      have to call S_.RemoveAllEdges(i), too. This can be
            //      omitted here due to the following reasons:
            /////
            // Here we know:
            //  (i) In case (coarseIndex_[i] == DIRICHLET_FINE) not a
            //      single edge has been added to graph S, since we did
            //      not even enter the loop, in which edges S(i,..) are
            //      set.
            // (ii) In case (auxMat.GetRowSize(i) == 1) we know from exactly
            //      this condition that row [i] has not got any offdiagonal
            //      entries, and therefore no edges S(i,..) either.  
            // only remove them in ST_
            ST_.RemoveAllEdges( i );
            // In both cases we set the splitting number to FINE.
            CoarseIndex_[i] = (Integer)FINE;
        }
    }

    // LAS:  build in the explicit Dirichlet values
    // OLAS: we do not have direct access to Dirichlet information

    // set the pointers for the neigbourhood relations
    NhStartIndex_ = pRow;
    NhEdges_      = pCol;

    return startCoarsePoint_;
}

/**********************************************************/

template <typename T>
Integer Topology<T>::
CalcCoarseFineSplitting( const Integer max_dependency,
                         const Integer gamma )
{

#ifdef TOPOLOGY_AVOID_REDUNDANT_IMPORTANCE_CALCULATION
    // create array to prevent redundant evaluation of "the
    // coarse importance" of possible C-points
    NEWARRAY( importanceKnown_, Integer, Size_h_ );
    for( Integer i = 0; i < (Integer)Size_h_; i++ )  importanceKnown_[i] = -3;
#endif
    
    // get the first coarse point (should be stored in
    // startCoarsePoint_ at this time)
    Integer coarsepoint      = GetNextCoarsePoint( 0, max_dependency ),
            lowest_undefined = 0;

    // while there are undefined points left
    while( coarsepoint < (Integer)Size_h_ ) {
        // while there are some candidates for coarse points around
        while( coarsepoint == UNDEFINED || coarsepoint > -1 ) {
            // set the coarse point (note that the surrounding points
            // are put into F, as long as the are UNDEFINED, and
            // increments the coarse point counter Size_H_)
            SetCoarsePoint( coarsepoint );
            // get the next possible coarse point
            coarsepoint = GetNextCoarsePoint( coarsepoint, max_dependency );

        }
        // we are stuck in this chain of coarse points
        // -> search the next UNDEFINED point, starting at point 0
        for( coarsepoint = lowest_undefined;
             coarsepoint < (Integer)Size_h_; coarsepoint++ )
        {
            if( CoarseIndex_[coarsepoint] == UNDEFINED ) {
                lowest_undefined = coarsepoint;
                break;
            }
        }
        // if we could not find an UNDEFINED point, coarsepoint will
        // be (Size_h_ + 1) -> we leave the loop
    }

    // eventual second pass
    if( gamma >= 0 ) {
        for( Integer i = 0; i < (Integer)Size_h_; i++ ) {
            if( IsFPoint(i)                       &&
                GetNumInterpolatingPoints(i) == 1 &&
                S_.GetNumEdges(i)            >  1    ) {
                Integer nn = 0;
                for( Integer ij = (Integer)NhStartIndex_[i]+0; //ROK: was +1
                             ij < (Integer)NhStartIndex_[i+1]; ij++ ) {
                    const Integer j = NhEdges_[ij];
                    if( IsFPoint(j)                       &&
                        GetNumInterpolatingPoints(j) == 1 &&
                        S_.GetNumEdges(j)            >  1    ) {
                        if( ++nn >= gamma )  break;
                    }
                }
                // since ALL points are either FINE or COARSE
                // at this point, we can simply set the coarse
                // index instead of calling SetCoarsePoint
                if( nn >= gamma )  CoarseIndex_[i] = ++Size_H_;
            }
        }
    }

#ifdef TOPOLOGY_AVOID_REDUNDANT_IMPORTANCE_CALCULATION
    delete [] ( importanceKnown_ );  importanceKnown_  = NULL;
#endif

    return Size_H_;
}

/**********************************************************/

#define ADD_EDGES_SORTED

#ifdef ADD_EDGES_SORTED
#define  _ADD_EDGE_  AddEdgeSorted
#else
#define  _ADD_EDGE_  AddEdgeSavely
#endif
/**********************************************************/

template <typename T>
void Topology<T>::CalcGalerkinGraphs( DependencyGraph<T>& graph_AHT,
                                      DependencyGraph<T>& graph_VT,
                                      const Integer stretch_factor ) const
{

    graph_AHT.Reset();
    graph_VT.Reset();

    // Calculate the maximal row length in the matrix of the coarse
    // system as maximal row length of the system matrix multiplied
    // with a constant factor. This might turn out as a dangerous
    // heuristic, but LAS used it successfully, so let's see.
    Integer maxNodeSize = 0;
    for( Integer i = 1; i <= Size_h_; i++ ) {
        if( maxNodeSize < NhStartIndex_[i+1] - NhStartIndex_[i] ) {
            maxNodeSize = NhStartIndex_[i+1] - NhStartIndex_[i];
        }
    }
    // arrays for buffering some graph edges
    //   CoNS_i == strong Coarse Neighbours of point i (== in S(i) n C)
    Integer *const CoNS_i = new Integer [maxNodeSize];
    Integer *const CoNS_k = new Integer [maxNodeSize];
    // stretch the maximal number of connections
    maxNodeSize *= stretch_factor;

    // create the empty graph for VT and the one for AHT, already
    // initialized with "diagonal edges"
    if( false == graph_AHT.CreateWithDiagonals(Size_H_, maxNodeSize) ) {
        EXCEPTION( "topology.cc: could not create graph for A_H"); }
    if( false == graph_VT.Create(Size_H_, maxNodeSize) ) {
      EXCEPTION( "topology.cc: could not create graph for VT"); }

    // the variables nomenclature:
    //   i,j,k : fine grid indices
    //   xC    : index of point x in the coarse system
    //   nCSx  : number of points in C, that strongly influence point x
    //   Nx    : pointer to points, that influence point x; N(x)
    //   nNx   : number of points, that influence point x; |N(x)|
    //   ix    : running index belonging to point x

    // i runs over all points in the grid
    for( Integer i = 1; i <= Size_h_; i++ ) {
        // iC = coarse index of point i (i still might be F-point!)
        const Integer iC = CoarseIndex_[i];
        // nCSi = |j in (S(i) n C)|
        Integer nCSi = 0;
        // i in C
        // the edges, inserted in both cases of the following
        // if-statement, correspond exactly to the ones, that would
        // be set inside the loop over ik, if we would permit a point
        // being in it own neighbourhood (i in N(i)). Taking it out of
        // that loop is due to performance aspects. (NOTE: we have set
        // all edges (iC,iC) in AHT at creation time of AHT)
        if( iC >= COARSE ) {
            graph_VT._ADD_EDGE_( iC, i );
            // here we should also insert (iC,iC) -> graph_AHT, but this
            // graph was already created including all "diagonals"
        } else {
            // nCSi = |j in (S(i) n C)|
            // note that CoNS_i is only needed in the whole algorithm,
            // if i is F-point; therefore it is OK to fill it only in
            // this case, i.e. at this position
            nCSi = WriteStrongCoarseNeighbours( i, CoNS_i );
            // As we will never need the fine grid indices of S(i) n C,
            // we map them all at once to C-indices
            for( Integer il = 0; il < nCSi; il++ ) {
                CoNS_i[il] = CoarseIndex_[CoNS_i[il]];
                graph_VT._ADD_EDGE_( CoNS_i[il], i );
            }
        }
        // Ni = N(i), nNi = |N(i)|
        // (by the way: (j in N(i)) <==> (A_ij != 0) AND i != j)
        const Integer* const  Ni = NhEdges_ + NhStartIndex_[i] + 1;
        const Integer        nNi = NhStartIndex_[i+1] - NhStartIndex_[i] - 1;
        // k in N(i), ( <==> A_ik != 0 )
        for( Integer ik = 0; ik < nNi; ik++ ) {
            const Integer k = Ni[ik];
            // k is F-point
            if( IsFPoint(k) ) {
                // nCSk = |S(k) n C|
                const Integer nCSk = WriteStrongCoarseNeighbours( k, CoNS_k );
                // l in (S(k) n C)
                for( Integer il = 0; il < nCSk; il++ ) {
                    const Integer lC = CoarseIndex_[CoNS_k[il]];
                    graph_VT._ADD_EDGE_( lC, i );
                    // i is F-point
                    if( iC == FINE ) {
                        // j in (S(i) n C)
                        for( Integer ij = 0; ij < nCSi; ij++ ) {
                            // "long" C-F-F-C connection: jC-iF-kF-lC
                            // TRANSPOSE edge!
                            graph_AHT._ADD_EDGE_( lC, CoNS_i[ij] );
                        }
                    // i is C-point
                    } else {
                        // "short" C-F-C connection: iC-kF-lC
                        // TRANSPOSE edge!
                        graph_AHT._ADD_EDGE_( lC, iC );
                    }
                }
            // k is C-point
            } else {
                const Integer kC = CoarseIndex_[Ni[ik]];
                graph_VT._ADD_EDGE_( kC, i );
                // i is F-point
                if( iC == FINE ) {
                    // j in (S(i) n C)
                    for( Integer ij = 0; ij < nCSi; ij++ ) {
                        // "short" C-F-C connection: jC-iF-kC
                        // TRANSPOSE edge!
                        graph_AHT._ADD_EDGE_( kC, CoNS_i[ij] );
                    }
                // i is C-point
                } else {
                    // "direct" C-C connection: iC-kC
                    // TRANSPOSE edge!
                    graph_AHT._ADD_EDGE_( kC, iC );
                }
            }
        } // k in N(i)
    } // i

    // sort the graphs
#ifndef ADD_EDGES_SORTED
    graph_VT.Sort();
    graph_AHT.Sort();
#endif

    delete [] CoNS_i;
    delete [] CoNS_k;

}

/**********************************************************/
#undef ADD_EDGES_SORTED
#undef _ADD_EDGE_

/**********************************************************/

template <typename T>
inline Integer Topology<T>::GetNumFineNeighbours( const Integer i ) const
{
    const Integer* const Neighbourhood   = NhEdges_ + NhStartIndex_[i];
    const Integer        nNeighbours     = NhStartIndex_[i+1] -
                                           NhStartIndex_[i];
    Integer              nFineNeighbours = 0;

    // start with index 1, since 0 would point i itself
    for( Integer ij = 1; ij < nNeighbours; ij++ ) {
        if( CoarseIndex_[Neighbourhood[ij]] == FINE )  nFineNeighbours++;
    }

    return nFineNeighbours;
}

/**********************************************************/

template <typename T>
Integer Topology<T>::GetNumInterpolatingPoints( Integer* const sizes ) const
{

    Integer totalNumPoints = 0;
    
    // run over all points
    for( Integer i = 0; i < S_.GetNumNodes(); i++ ) {
        // coarse points are only interpolated from themselves
        if( CoarseIndex_[i] >= COARSE ) {
            sizes[i] = 1;
        // fine points (we assume that there do not exist any undefined
        // points) are interpolated from their coarse neigbourhood, they
        // strongly depend from
        } else {
                  Integer        size   = 0;
            const Integer* const edges  = S_.GetEdges( i );
            const Integer        nedges = S_.GetNumEdges( i );
            // inkrement size for each edge to a coarse node
            for( Integer ij = 0; ij < nedges; ij++ ) {
                if( CoarseIndex_[edges[ij]] >= COARSE )  size++;
            }
            sizes[i] = size;
        }
        totalNumPoints += sizes[i];
    }

    return totalNumPoints;
}

/**********************************************************/

template <typename T>
Integer Topology<T>::GetNumInterpolatingPoints( const Integer i ) const
{

    Integer NumPoints = 0;

    if( CoarseIndex_[i] >= COARSE ) {
        return 1;
    } else {
        const Integer* const edges  = S_.GetEdges( i );
        const Integer        nedges = S_.GetNumEdges( i );
        for( Integer ij = 0; ij < nedges; ij++ ) {
            if( CoarseIndex_[edges[ij]] >= COARSE )  NumPoints++;
        }
    }

    return NumPoints;
}

/**********************************************************/

template <typename T>
Integer Topology<T>::GetNumInterpolatingPoints() const
{

    Integer nPoints = 0;
    
    // run over all points
    for( Integer i = 0; i < S_.GetNumNodes(); i++ ) {
        // coarse points are only interpolated from themselves
        if( CoarseIndex_[i] >= COARSE ) {
            nPoints++;
        // fine points (we assume that there do not exist any undefined
        // points) are interpolated from their coarse neigbourhood, they
        // strongly depend from
        } else {
            if( CoarseIndex_[i] >= COARSE ) {
                nPoints++;
            } else {
                      Integer        size   = 0;
                const Integer* const edges  = S_.GetEdges( i );
                const Integer        nedges = S_.GetNumEdges( i );
                // inkrement size for each edge to a coarse node
                for( Integer ij = 0; ij < nedges; ij++ ) {
                    if( CoarseIndex_[edges[ij]] >= COARSE )  nPoints++;
                }
            }
        }
    }

    return nPoints;
}

/**********************************************************/

template <typename T>
Integer Topology<T>::
GetCoarseImportance( Integer const i ) const
{

    const Integer        nedges = ST_.GetNumEdges( i );
    const Integer* const  edges = ST_.GetEdges( i );
    // initialize the result with |ST(i)|
    Integer value = nedges;
    for( Integer ij = 0; ij < nedges; ij++ ) {
        // increment value for each j in (ST(i) n F)
        // NOTE that we cannot use method IsFPoint(i) here, since
        //      IsFPoint(i) prints a warning, if i is UNDEFINED.
        //      But since this method is called in setup phase
        //      i == UNDEFINED is completely normal.
        if( CoarseIndex_[edges[ij]] == FINE )  value++;
    }

    return value;
}

/**********************************************************/

template <typename T>
Integer Topology<T>::
GetNumInterpolatedPoints( StdVector<UInt>& sizes ) const
{
    sizes.Resize( this->GetSizeh() );
    sizes.Init(0);

    Integer totalNumPoints = 0;

    // run over all fine grid points
    for( Integer i = 0; i < (Integer)GetSizeh(); i++ ) {
        // only coarse points interpolate other points
        if( IsCPoint(i) ) {
            // initialize the counter with 1, because every coarse
            // point interpolates at least one point, namely itself
                  Integer        size   = 1;
            const Integer* const edges  = ST_.GetEdges( i );
            const Integer        nedges = ST_.GetNumEdges( i );
            // increment counter for each edge to an F-node
            for( Integer ij = 0; ij < nedges; ij++ ) {
                if( IsFPoint(edges[ij]) )  size++;
            }
            totalNumPoints += size;
            sizes[CoarseIndex_[i]] = size;
        }
    }

    return totalNumPoints;
}

/**********************************************************/

template <typename T>
Integer Topology<T>::
WriteStrongCoarseNeighbours( const Integer        p,
                                   Integer* const neighbours ) const
{
          Integer        nSCN = 0;
    const Integer        nSN  = S_.GetNumEdges( p );
    const Integer* const  SN  = S_.GetEdges( p );

    for( Integer i = 0; i < nSN; i++ ) {
        if( CoarseIndex_[SN[i]] > UNDEFINED ) {
            neighbours[nSCN++] = SN[i];
        }
    }

    return nSCN;
}

/**********************************************************/

template <typename T>
Integer Topology<T>::
WriteCoarseNeighbours( const Integer        p,
                             Integer* const c_neighbours ) const
{
          Integer        nCN = 0;
    const Integer        nN  = NhStartIndex_[p+1] - NhStartIndex_[p];
    const Integer* const  N  = NhEdges_ + NhStartIndex_[p];

    for( Integer i = 1; i < nN; i++ ) {
        if( CoarseIndex_[N[i]] >= COARSE ) {
            c_neighbours[nCN++] = N[i];
        }
    }

    return nCN;
}

/**********************************************************/

template <typename T>
Integer Topology<T>::
WriteFineNeighbours( const Integer        p,
                           Integer* const f_neighbours ) const
{
          Integer        nFN = 0;
    const Integer        nN  = NhStartIndex_[p+1] - NhStartIndex_[p];
    const Integer* const  N  = NhEdges_ + NhStartIndex_[p];

    for( Integer i = 1; i < nN; i++ ) {
        if( CoarseIndex_[N[i]] == FINE ) {
            f_neighbours[nFN++] = N[i];
        }
    }

    return nFN;
}

/**********************************************************/

template <typename T>
Integer Topology<T>::GetNextCoarsePoint( const Integer p,
                                         const Integer max_dependency )
{
    // if there is already a next coarse point ready selected (for
    // example after the call of CreateDependencyGraphs), return this
    // point and erase it from the class data
    if( startCoarsePoint_ > -1 &&
        CoarseIndex_[startCoarsePoint_] == UNDEFINED )
    {
        const Integer point = startCoarsePoint_;
        startCoarsePoint_ = -1;
        return point;
    // otherwise we must search the next coarse grid point.
    // some notes about the names of the indices:
    //   N_x   : neighbourhood of point x
    //   nN_x  : |N_x| == number of nodes in N_x
    //   iNx   : loop index for the i-th node in N_x
    //   ST_x  : { i | (x,i) in ST }
    //   nST_x : |ST_x|
    } else {

        const Integer        nST_p = ST_.GetNumEdges(p);
        const Integer* const  ST_p = ST_.GetEdges(p);
        Integer max_depend = 0,
                nextcoarse = -1;
        // walk through ST(p)
        for( Integer iSTp = 0; iSTp < nST_p; iSTp++ ) {
            // i = real index of node number nN_p in N(p)
            const Integer  i = ST_p[iSTp];
            const Integer  nN_i = NhStartIndex_[i+1] - NhStartIndex_[i];
            const UInt* const  N_i = NhEdges_ + NhStartIndex_[i];
            // walk through N(i)
            for( UInt iNi = 0; iNi < (UInt)nN_i; iNi++ ) {
                const UInt j = N_i[iNi];
                // if we found an UNDEFINED point, evaluate it

//#ifdef TOPOLOGY_AVOID_REDUNDANT_IMPORTANCE_CALCULATION
//std::cout<<    CoarseIndex_[j]<<std::endl;
//std::cout<<    importanceKnown_[j]<<std::endl;
//                if( CoarseIndex_[j] == UNDEFINED &&
                    // only check the "coarse importance" of point j, if
                    // we have not calculated it in this function call
//                    importanceKnown_[j] < Size_H_ ) {
//                    importanceKnown_[j] = Size_H_;
//#else
                if( CoarseIndex_[j] == UNDEFINED ) {
//#endif
                    // cImportance = |ST(j)| + |ST(j) n F|
                    const Integer cImportance = GetCoarseImportance( j );
                    // the measure for the value of point j as C-point is
                    // |N(j) n F| + |ST(j)|
                    if( cImportance > max_depend ) {
                        max_depend = cImportance;
                        // if the value exceeds the passed maximum,
                        // take j immediately as C-point
                        if( max_depend >= max_dependency )  return j;
                        // else continue the search for the point with
                        // the best value
                        else  nextcoarse = j;
                    }
                }
            }
        }
        
        return nextcoarse;
    }

}

/**********************************************************/

template <typename T>
inline void Topology<T>::SetCoarsePoint( const Integer i )
{
    
    // put node [i] into C
    CoarseIndex_[i] = Size_H_;
    Size_H_++;
    // fix all undefined points, that are strongly influenced by point
    // [i], as new points in F
    const Integer        nStronglyInfluenced = ST_.GetNumEdges( i );
    const Integer* const  StronglyInfluenced = ST_.GetEdges( i );
    for( Integer ij = 0; ij < nStronglyInfluenced; ij++ ) {
        if( CoarseIndex_[StronglyInfluenced[ij]] == UNDEFINED ) {
            CoarseIndex_[StronglyInfluenced[ij]] = (Integer)FINE;
        }
    }
}

/**********************************************************/

template <typename T>
inline void Topology<T>::SetDirichlet( const Integer i )
{
    
    S_.RemoveAllEdges( i );
    ST_.RemoveAllEdges( i );
    CoarseIndex_[i] = (Integer)FINE;
}

/**********************************************************/

template <typename T>
void Topology<T>::Reset()
{
    
    S_.Reset();
    ST_.Reset();

    // reset neighbourhood pointers to NULL
    NhStartIndex_ = NULL;
    NhEdges_      = NULL;
    // delete array of C-F-Splitting
    delete [] ( CoarseIndex_ );  CoarseIndex_  = NULL;
    CoarseIndex_ = NULL;
    
    startCoarsePoint_ = -1;
    Size_h_ = 0;
    Size_H_ = 0;
}


} // namespace CoupledField

/**********************************************************/
#ifdef DEBUG_TO_CERR
#undef debug
#endif // DEBUG_TO_CERR
/**********************************************************/
