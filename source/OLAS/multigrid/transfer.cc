/* $Id$ */

#include "multigrid/transfer.hh"
#include "multigrid/prematrix.hh"
    
/**********************************************************/
#ifdef DEBUG_TO_CERR
#ifndef DEBUG_TRANSFEROPERATOR
#define DEBUG_TRANSFEROPERATOR
#endif // DEBUG_TRANSFEROPERATOR
#define  debug  &std::cerr
#endif // DEBUG_TO_CERR

#ifdef PROFILING
#ifdef PROFILE_TRANSFEROPERATOR
#include "utils/utils.hh"
#endif
#endif
/**********************************************************/
#ifdef TRANSFER_OPERATOR_IMPORT_GALERKIN_PRODUCT
#include "exporttools.hh"
// If the import prefix for the Galerkin product is not yet
// defined, but there is a global import prefix, use the global one.
#ifndef TRANSFER_OPERATOR_GALERKIN_PRODUCT_IMPORT_PREFIX
#ifdef  AMG_IMPORT_PREFIX
#define TRANSFER_OPERATOR_GALERKIN_PRODUCT_IMPORT_PREFIX  AMG_IMPORT_PREFIX
#else  // AMG_IMPORT_PREFIX
#define TRANSFER_OPERATOR_GALERKIN_PRODUCT_IMPORT_PREFIX ""
#endif // AMG_IMPORT_PREFIX
#endif // TRANSFER_OPERATOR_GALERKIN_PRODUCT_IMPORT_PREFIX
#endif // TRANSFER_OPERATOR_IMPORT_GALERKIN_PRODUCT

#ifdef TRANSFER_OPERATOR_IMPORT_INTERPOLATION
// If the import prefix for the interpolation operator is not yet
// defined, but there is a global import prefix, use the global one.
#ifndef TRANSFER_OPERATOR_INTERPOLATION_IMPORT_PREFIX
#ifdef  AMG_IMPORT_PREFIX
#define TRANSFER_OPERATOR_INTERPOLATION_IMPORT_PREFIX  AMG_IMPORT_PREFIX
#else  // AMG_IMPORT_PREFIX
#define TRANSFER_OPERATOR_INTERPOLATION_IMPORT_PREFIX ""
#endif // AMG_IMPORT_PREFIX
#endif // TRANSFER_OPERATOR_INTERPOLATION_IMPORT_PREFIX
#endif // TRANSFER_OPERATOR_IMPORT_INTERPOLATION
/**********************************************************/

namespace OLAS {
/**********************************************************/

template <typename T>
TransferOperator<T>::TransferOperator()
                   : Prolongation_(NULL),
                     Restriction_(NULL)
{
    ENTER_FCN("TransferOperator::TransferOperator()");
}

/**********************************************************/

template <typename T>
TransferOperator<T>::~TransferOperator()
{
    ENTER_FCN("TransferOperator::~TransferOperator()");
    
    Reset();
}

/**********************************************************/

template <typename T>
void TransferOperator<T>::Reset()
{
    ENTER_FCN("TransferOperator::Reset()");

    delete Prolongation_;  Prolongation_ = NULL;
    delete Restriction_;   Restriction_  = NULL;
}

/**********************************************************
 * Note: The restriction operator built by this algorithm must
 * have sorted (with respect to their column indices) rows. If this
 * fact should change once, the Galerkin product must be adapted.
 **********************************************************/

template <typename T>
bool TransferOperator<T>::
CreateOperators( const CRS_Matrix<T>& matrix,
                 const Topology<T>&   topology,
                 const AMGInterpolationType itype,
                 const bool build_interpolation )
{
    ENTER_FCN("TransferOperator::CreateOperators");

    ////////////////////////////////////////////////////////
    // eventual import of interpolation
    ////////////////////////////////////////////////////////

#ifdef TRANSFER_OPERATOR_IMPORT_INTERPOLATION
    ImportInterpolation( topology.GetSizeh(), topology.GetSizeH() );
    return true;
#else // TRANSFER_OPERATOR_IMPORT_INTERPOLATION


    ////////////////////////////////////////////////////////
    // time measurement and debugging checks
    ////////////////////////////////////////////////////////

#ifdef PROFILE_TRANSFEROPERATOR
    const Double t1 = AMG_GET_REAL_TIME
#endif

#ifdef DEBUG_TRANSFEROPERATOR
    if( matrix.GetNcols() != topology.GetSizeh() ) {
        (*debug) << "ERROR: # columns of the matrix == " << matrix.GetNcols()
                 << std::endl << "       != size_h in topology   == "
                 << topology.GetSizeh() << std::endl << "must match!"
                 << std::endl;
        Error( "Dimensions of matrix and topology do not match",
                __FILE__, __LINE__ );
    }
#endif

    Reset();

    ////////////////////////////////////////////////////////
    // create specific interpolation
    ////////////////////////////////////////////////////////

    bool creationResult = false;

    switch( itype ) {
        // constant interpolation
        case AMG_INTERPOLATION_CONSTANT:
            creationResult =
            CreateOperatorsConstant( matrix, topology );
            break;
        // simple weighted average interpolation
        case AMG_INTERPOLATION_SIMPLE_WEIGHTED:
            creationResult =
            CreateOperatorsSimpleWeighted( matrix, topology );
            break;
        //
        case AMG_INTERPOLATION_SMOOTHED_SCALING:
            creationResult =
            CreateOperatorsSmoothedScaling( matrix, topology );
            break;
        // developing version of interpolation
        case AMG_INTERPOLATION_DEVELOP:
            creationResult =
            CreateOperatorsDevelop( matrix, topology );
            break;
        //
        default:
            Error( __FILE__, __LINE__,
                   "TransferOperator::CreateOperators: non "
                   "supported constant %d as interpolation type\n",
                   itype );
            break;
    }

    ////////////////////////////////////////////////////////
    // debug checks
    ////////////////////////////////////////////////////////

#ifdef DEBUG_TRANSFEROPERATOR
    creationResult = CheckOperators();
#endif
    
#ifdef PROFILE_TRANSFEROPERATOR
    const Double t2 = AMG_GET_REAL_TIME
    (*cla) << " AMG: created transfer operators: "
           << (t2-t1)<<" seconds\n";
#endif

    return creationResult;

#endif // TRANSFER_OPERATOR_IMPORT_INTERPOLATION
}

/**********************************************************
 * Note: The restriction operator built by this algorithm has
 * sorted (with respect to their column indices) rows. If this
 * fact should change once, the Galerkin product must be adapted.
 **********************************************************/

template <typename T>
bool TransferOperator<T>::
CreateOperatorsConstant( const CRS_Matrix<T>& matrix,
                         const Topology<T>&   topology )
{
    ENTER_FCN("TransferOperator::CreateOperatorsConstant");

    Integer  nNonZeros  = 0,
            *RowLengths = NULL;
    NewArray( RowLengths, Integer, topology.GetSizeH() );
    // get number of interpolated points for each C-point
    nNonZeros = topology.GetNumInterpolatedPoints( RowLengths );

    // create the restriction and prolongation matrix object
    Restriction_ =
    New CRS_Matrix<T>( topology.GetSizeH(),
                             topology.GetSizeh(),
                             nNonZeros );
    Prolongation_ =
    New CRS_Matrix<T>( topology.GetSizeh(),
                             topology.GetSizeH(),
                             nNonZeros );
    // from now on nNonZeros is the index of the next free index
    // in pDatP and pColP (see below)
    ((Integer*)Prolongation_->GetRowPointer())[1] = nNonZeros = 1;

    // get some pointers for performant and dirty (CAST!) access
    Integer* const pRowR = (Integer*)Restriction_->GetRowPointer();
    Integer* const pColR = (Integer*)Restriction_->GetColPointer();
    T* const pDatR = (T*)Restriction_->GetDataPointer();
    // pointers for interpolation
    Integer* const pRowP = (Integer*)Prolongation_->GetRowPointer();
    Integer* const pColP = (Integer*)Prolongation_->GetColPointer();
    T* const pDatP = (T*)Prolongation_->GetDataPointer();
    // pointer to the splitting and fine-coarse mapping
    const Integer* const CoarseIndex = topology.GetCoarseFineSplitting();

    // set the row starts and initialize the array RowLengths (after
    // this loop we will use it to keep at position [i] the position of
    // the next free entry in row [i] of the restriction matrix
    pRowR[1] = 1;
    for( Integer i = 1; i <= Restriction_->GetNrows(); i++ ) {
        pRowR[i+1]    = pRowR[i] + RowLengths[i];
        // change the functionality of RowLengths[i] (see above)
        RowLengths[i] = pRowR[i];
    }
    // for the interpolation matrix the row length needs not to be fixed
    // now, because it will be filled row by row.

    // fill restriction matrix //
    // Run over all points in the fine systems. Note that this means
    // running over the COLUMNS of the restriction matrix. This "direction"
    // is chosen, since the interpolation weights can be reused inside
    // one column.
    for( Integer i = 1; i <= topology.GetSizeh(); i++ ) {
        // if point [i] is C-point, in column [i] there is only one
        // entry, since C-points are simply injected with weight 1.0
        const Integer iC = CoarseIndex[i];
        if( iC >= Topology<T>::COARSE ) {
            // R(iC,i) = 1
            pColR[RowLengths[iC]  ] = i;
            pDatR[RowLengths[iC]++] = (T)1.0;
            // P(i,iC) = 1
            pColP[nNonZeros  ] = iC;
            pDatP[nNonZeros++] = (T)1.0;
        } else {
                  Integer   denominator = 0;
            const Integer        nEdges = topology.GetS().GetNumEdges( i );
            const Integer* const  Edges = topology.GetS().GetEdges( i );
            // build interpolation for row [i] in P
            for( Integer ij = 0; ij < nEdges; ij++ ) {
                const Integer  j = Edges[ij],
                              jC = CoarseIndex[j];
                if( jC >= Topology<T>::COARSE ) {
                    pColP[nNonZeros++] = jC;
                    denominator++;
                }
            }
            if( denominator ) {
                const T weight = opType<T>::invert(
                                          static_cast<T>(denominator));
                for( Integer ij = pRowP[i]; ij < nNonZeros; ij++ ) {
                    const Integer jC = pColP[ij];
#ifdef DEBUG_TRANSFEROPERATOR
                        if( RowLengths[jC] >= pRowR[jC+1] ) {
                            (*debug) << "ERROR: TransferOperator<T>::CreateOperators:"
                                     << std::endl << "number of weights exceeds row "
                                        "length " << (pRowR[jC+1] < pRowR[jC])
                                     << " in row " << jC << std::endl;
                            Error( "TransferOperator::CreateOperators: "
                                   "-> see debug file", __FILE__, __LINE__ );
                        }
#endif
                    // set the interpolation weight ...
                    pDatP[ij] = weight;
                    // ... and set it in the restriction, too.
                    pColR[RowLengths[jC]  ] = i;
                    pDatR[RowLengths[jC]++] = weight;
                }
            } // if( denominator )
        }
        // terminate row i of P by start index of the next one
        pRowP[i+1] = nNonZeros;
    }

    DeleteArray( RowLengths );

    return true;
}

/**********************************************************
 * Note: The restriction operator built by this algorithm has
 * sorted (with respect to their column indices) rows. If this
 * fact should change once, the Galerkin product must be adapted.
 **********************************************************/

template <typename T>
bool TransferOperator<T>::
CreateOperatorsSimpleWeighted( const CRS_Matrix<T>& matrix,
                               const Topology<T>&   topology )
{
    ENTER_FCN("TransferOperator::CreateOperators_new");

    Integer  nNonZeros  = 0,
            *RowLengths = NULL;
    NewArray( RowLengths, Integer, topology.GetSizeH() );
    // get number of interpolated points for each C-point
    nNonZeros = topology.GetNumInterpolatedPoints( RowLengths );

    // create restriction matrix object
    Restriction_ =
    New CRS_Matrix<T>( topology.GetSizeH(),
                             topology.GetSizeh(),
                             nNonZeros );
    // create prolongation matrix object
    Prolongation_ =
    New CRS_Matrix<T>( topology.GetSizeh(),
                             topology.GetSizeH(),
                             nNonZeros );
    // from now on nNonZeros is the index of the next free index
    // in pDatP and pColP (see below)
    ((Integer*)Prolongation_->GetRowPointer())[1] = nNonZeros = 1;

    // get some pointers for performant and dirty (CAST!) access
    // pointers for the problem matrix
    const Integer* const pRowA = matrix.GetRowPointer();
    const Integer* const pColA = matrix.GetColPointer();
    const T* const pDatA = matrix.GetDataPointer();
    // pointers for the restriction
    Integer* const pRowR = (Integer*)Restriction_->GetRowPointer();
    Integer* const pColR = (Integer*)Restriction_->GetColPointer();
    T* const pDatR = (T*)Restriction_->GetDataPointer();
    // pointers for interpolation
    Integer* const pRowP = (Integer*)Prolongation_->GetRowPointer();
    Integer* const pColP = (Integer*)Prolongation_->GetColPointer();
    T* const pDatP = (T*)Prolongation_->GetDataPointer();
    // pointer to the splitting and fine-coarse mapping
    const Integer* const CoarseIndex = topology.GetCoarseFineSplitting();

    // set the row starts and initialize the array RowLengths (after
    // this loop we will use it to keep at position [i] the position of
    // the next free entry in row [i] of the restriction matrix
    pRowR[1] = 1;
    for( Integer i = 1; i <= Restriction_->GetNrows(); i++ ) {
        pRowR[i+1]    = pRowR[i] + RowLengths[i];
        // change the functionality of RowLengths[i] (see above)
        RowLengths[i] = pRowR[i];
    }
    // for the interpolation matrix the row length needs not to be fixed
    // now, because it will be filled row by row.
    
    for( Integer i = 1; i <= topology.GetSizeh(); i++ ) {
        // if point [i] is C-point, in column [i] there is only one
        // entry, since C-points are simply injected with weight 1.0
        const Integer iC = CoarseIndex[i];
        if( iC >= Topology<T>::COARSE ) {
            // R(iC,i) = 1
            pColR[RowLengths[iC]  ] = i;
            pDatR[RowLengths[iC]++] = (T)1.0;
            // P(i,iC) = 1
            pColP[nNonZeros  ] = iC;
            pDatP[nNonZeros++] = (T)1.0;
        } else {
            T         denominator = (T)0.0;
            const Integer        nEdges = topology.GetS().GetNumEdges( i );
            const Integer* const  Edges = topology.GetS().GetEdges( i );
                  Integer           ijA = pRowA[i];
            // build interpolation for row [i] in P
            for( Integer ij = 0; ij < nEdges; ij++ ) {
                const Integer  j = Edges[ij],
                              jC = CoarseIndex[j];
                if( jC >= Topology<T>::COARSE ) {
                    // search the corresponding entry in A
                    while( ijA < pRowA[i+1] ) {
                        if( pColA[ijA] == j ) {
                            pColP[nNonZeros  ]  = jC;
                            denominator        +=
                            pDatP[nNonZeros++]  = Abs(pDatA[ijA]);
                            break;
                        }
                        ijA++;
                    }
                }
            }
            if( denominator ) {
                const T factor = opType<T>::invert(denominator);
                for( Integer ij = pRowP[i]; ij < nNonZeros; ij++ ) {
                    const Integer jC = pColP[ij];
                    // finish the interpolation weight ...
                    pDatP[ij] *= factor;
                    // ... and set it in the restriction, too.
                    pColR[RowLengths[jC]  ] = i;
                    pDatR[RowLengths[jC]++] = pDatP[ij];
                }
            }
        }
        // terminate row i of P by start index of the next one
        pRowP[i+1] = nNonZeros;
    }

    DeleteArray( RowLengths );
    
    return true;
}

/**********************************************************
 * Note: The restriction operator built by this algorithm must
 * have sorted (with respect to their column indices) rows.
 * Otherwise the Galerkin product must be adapted.
 **********************************************************/

template <typename T>
bool TransferOperator<T>::
CreateOperatorsSmoothedScaling( const CRS_Matrix<T>& matrix,
                                const Topology<T>&   topology )
{
    //////////////////////////////////////////////
    // Calculate scaling factors
    //////////////////////////////////////////////

    T *v = 0x0;
    NewArray( v, T, topology.GetSizeh() );

    const Integer* const pRowA = matrix.GetRowPointer();
    const Integer* const pColA = matrix.GetColPointer();
    const T* const pDatA = matrix.GetDataPointer();
    const Integer        sizeh = topology.GetSizeh();

    // emulate one Gauss-Seidel iteration with zero right hand side
    for( Integer i = 1; i <= sizeh; i++ )  v[i] = (T)1.0;
    for( Integer i = 1; i <= sizeh; i++ ) {
        T vi = 0;
        for( Integer ij = pRowA[i]+1; ij < pRowA[i+1]; ij++ ) {
            vi -= pDatA[ij] * v[pColA[ij]];
        }
        v[i] = vi / pDatA[pRowA[i]];
    }

    //////////////////////////////////////////////
    // Here follows the exact code copy of CreateOperatorsConstant
    // augmented by the scaling with v[i].
    // It would be also possible to create first the operators by
    // a call of CreateOperatorsConstant and apply the scaling
    // afterwards, but this version is slightly faster.
    //////////////////////////////////////////////

    Integer  nNonZeros  = 0,
            *RowLengths = NULL;
    NewArray( RowLengths, Integer, topology.GetSizeH() );
    // get number of interpolated points for each C-point
    nNonZeros = topology.GetNumInterpolatedPoints( RowLengths );

    // create the restriction and prolongation matrix object
    Restriction_ =
    New CRS_Matrix<T>( topology.GetSizeH(),
                             topology.GetSizeh(),
                             nNonZeros );
    Prolongation_ =
    New CRS_Matrix<T>( topology.GetSizeh(),
                             topology.GetSizeH(),
                             nNonZeros );
    // from now on nNonZeros is the index of the next free index
    // in pDatP and pColP (see below)
    ((Integer*)Prolongation_->GetRowPointer())[1] = nNonZeros = 1;

    // get some pointers for performant and dirty (CAST!) access
    Integer* const pRowR = (Integer*)Restriction_->GetRowPointer();
    Integer* const pColR = (Integer*)Restriction_->GetColPointer();
    T* const pDatR = (T*)Restriction_->GetDataPointer();
    // pointers for interpolation
    Integer* const pRowP = (Integer*)Prolongation_->GetRowPointer();
    Integer* const pColP = (Integer*)Prolongation_->GetColPointer();
    T* const pDatP = (T*)Prolongation_->GetDataPointer();
    // pointer to the splitting and fine-coarse mapping
    const Integer* const CoarseIndex = topology.GetCoarseFineSplitting();

    // set the row starts and initialize the array RowLengths (after
    // this loop we will use it to keep at position [i] the position of
    // the next free entry in row [i] of the restriction matrix
    pRowR[1] = 1;
    for( Integer i = 1; i <= Restriction_->GetNrows(); i++ ) {
        pRowR[i+1]    = pRowR[i] + RowLengths[i];
        // change the functionality of RowLengths[i] (see above)
        RowLengths[i] = pRowR[i];
    }
    // for the interpolation matrix the row length needs not to be fixed
    // now, because it will be filled row by row.

    // fill restriction matrix //
    // Run over all points in the fine systems. Note that this means
    // running over the COLUMNS of the restriction matrix. This "direction"
    // is chosen, since the interpolation weights can be reused inside
    // one column.
    for( Integer i = 1; i <= topology.GetSizeh(); i++ ) {
        // if point [i] is C-point, in column [i] there is only one
        // entry, since C-points are simply injected with weight 1.0
        const Integer iC = CoarseIndex[i];
        if( iC >= Topology<T>::COARSE ) {
            // R(iC,i) = 1
            pColR[RowLengths[iC]  ] = i;
            pDatR[RowLengths[iC]++] = v[i];
            // P(i,iC) = 1
            pColP[nNonZeros  ] = iC;
            pDatP[nNonZeros++] = v[i];
        } else {
                  Integer   denominator = 0;
            const Integer        nEdges = topology.GetS().GetNumEdges( i );
            const Integer* const  Edges = topology.GetS().GetEdges( i );
            // build interpolation for row [i] in P
            for( Integer ij = 0; ij < nEdges; ij++ ) {
                const Integer  j = Edges[ij],
                              jC = CoarseIndex[j];
                if( jC >= Topology<T>::COARSE ) {
                    pColP[nNonZeros++] = jC;
                    denominator++;
                }
            }
            if( denominator ) {
                const T weight = v[i] / denominator;
                for( Integer ij = pRowP[i]; ij < nNonZeros; ij++ ) {
                    const Integer jC = pColP[ij];
#ifdef DEBUG_TRANSFEROPERATOR
                        if( RowLengths[jC] >= pRowR[jC+1] ) {
                            (*debug) << "ERROR: TransferOperator<T>::CreateOperators:"
                                     << std::endl << "number of weights exceeds row "
                                        "length " << (pRowR[jC+1] < pRowR[jC])
                                     << " in row " << jC << std::endl;
                            Error( "TransferOperator::CreateOperators: "
                                   "-> see debug file", __FILE__, __LINE__ );
                        }
#endif
                    // set the interpolation weight ...
                    pDatP[ij] = weight;
                    // ... and set it in the restriction, too.
                    pColR[RowLengths[jC]  ] = i;
                    pDatR[RowLengths[jC]++] = weight;
                }
            } // if( denominator )
        }
        // terminate row i of P by start index of the next one
        pRowP[i+1] = nNonZeros;
    }

    DeleteArray( RowLengths );
    DeleteArray( v );

    return true;
}

/**********************************************************
 * Note: The restriction operator built by this algorithm must
 * have sorted (with respect to their column indices) rows.
 * Otherwise the Galerkin product must be adapted.
 **********************************************************/

template <typename T>
bool TransferOperator<T>::
CreateOperatorsDevelop( const CRS_Matrix<T>& matrix,
                        const Topology<T>&   topology )
{
    Error( "TransferOperator::CreateOperatorsDevelop: currently "
           "not implemented\n", __FILE__, __LINE__ );
    
    return false;
}

/**********************************************************/

template <typename T>
void TransferOperator<T>::Prolongate( const Vector<T>& v_H,
                                            Vector<T>& v_h,
                                      const bool       add  ) const
{
    ENTER_FCN("TransferOperator::Prolongate");

#ifdef DEBUG_TRANSFEROPERATOR
    // Since the restriction operator is always build, its absence
    // is a sure sign that the transfer operator has not been set
    // up, and therefore there is no prolongation matrix, either.
    if( Restriction_ == NULL ) {
        Error( "cannot prolongate vector, transfer operator not"
               " yet set up", __FILE__, __LINE__ );
    }
    // if there is a prolongation matrix, we want to use it
    if( Prolongation_ != NULL ) {
        if( Prolongation_->GetNcols() != v_H.GetSize() ||
            Prolongation_->GetNrows() != v_h.GetSize()    ) {
            (*debug) << "ERROR: dimensions of vectors and prolongation "
                        "operator do not match" << std::endl << "       "
                        "prolongation : " << Prolongation_->GetNrows()
                     << " x " << Prolongation_->GetNcols() << std::endl
                     << "       v_h         : " << v_h.GetSize() << std::endl
                     << "       v_H         : " << v_H.GetSize() << std::endl;
            Error( "dimensions for prolongation do not match",
                    __FILE__, __LINE__ );
        }
    // without prolongation matrix we use the restriction
    } else {
        if( Restriction_->GetNcols() != v_h.GetSize() ||
            Restriction_->GetNrows() != v_H.GetSize()    ) {
            (*debug) << "ERROR: dimensions of vectors and restriction "
                        "operator do not match" << std::endl << "       "
                        "restriction : " << Restriction_->GetNrows()
                     << " x " << Restriction_->GetNcols() << std::endl
                     << "       v_h         : " << v_h.GetSize() << std::endl
                     << "       v_H         : " << v_H.GetSize() << std::endl;
            Error( "dimensions for prolongation do not match",
                    __FILE__, __LINE__ );
        }
    }
#endif

    if( Prolongation_ != NULL ) {
        if( add )  Prolongation_->MultAdd( v_H, v_h );
        else       Prolongation_->Mult( v_H, v_h );
    } else {
        if( add )  Restriction_->MultTAdd( v_H, v_h );
        else       Restriction_->MultT( v_H, v_h );
    }
}

/**********************************************************/

template <typename T>
void TransferOperator<T>::Restrict( const Vector<T>& v_h,
                                          Vector<T>& v_H,
                                    const bool       add  ) const
{
    ENTER_FCN("TransferOperator::Restrict");

#ifdef DEBUG_TRANSFEROPERATOR
    if( Restriction_ == NULL ) {
        Error( "cannot restrict vector, restriction operator not"
               " yet built", __FILE__, __LINE__ );
    }
    if( Restriction_->GetNcols() != v_h.GetSize() ||
        Restriction_->GetNrows() != v_H.GetSize()    ) {
        (*debug) << "ERROR: dimensions of vectors and restriction "
                    "operator do not match" << std::endl << "       "
                    "restriction : " << Restriction_->GetNrows()
                 << " x " << Restriction_->GetNcols() << std::endl
                 << "       v_h         : " << v_h.GetSize() << std::endl
                 << "       v_H         : " << v_H.GetSize() << std::endl;
        Error( "dimensions for restriction do not match",
                __FILE__, __LINE__ );
    }
#endif

    if( add )  Restriction_->MultAdd( v_h, v_H );
    else       Restriction_->Mult( v_h, v_H );
}

/**********************************************************
 * This function computes the coarse system matrix AH via the Galerkin
 * product AH = R * Ah * P, where R is the restriction matrix, Ah is the
 * system matrix, and P is the prolongation matrix, with P = R^T (== R
 * R transposed.
 **********************************************************
 *
 * In an arbitrary matrix product A*B, with A and be stored in the CRS
 * format there is the problem, that we must run over columns of B. So
 * It would be much better to have B stored as B^T. Since we have stored
 * R in compressed row storage format, we multiply Ah * P^T. So AH is
 * computed as
 *                          n
 *                         ---
 *               AH[i,j] = >   R[i,k] * AH[k,l] * R[j,l]
 *                         ---
 *                       k,l = 1
 *
 * The resulting vector of the sum over l of AH[k,l] * R[j,l] is stored
 * in a temporary, also compressed, vector V for each j. So the following
 * expression describes the algorithm in a better way
 *
 *             n                                 n
 *            ---                               ---
 *  AH[i,j] = >   R[i,k] * Vj[k] , with Vj[k] = >   AH[k,l] * R[j,l]
 *            ---                               ---
 *           k = 1                             l = 1
 *
 * We exploit the fact that the rows of Ah are sorted. As a result the
 * entries of Vj and the rows of AH are sorted, too. We must only note
 * that the diagonal entries of Ah and AH are and must be stored at the
 * first position in the row. NOTE: The restriction matrix has not got
 * this property and also does not need it.
 * As the resulting matrix is also sparse, he passed graph of the coarse
 * matrix tells us, which entries must be computed. Therefore this graph
 * must be sorted to (no special positions of "diagonal edges). Since
 * we can reuse the temporary vector Vj for all valid rows i in AH, the
 * outer loop of this calculation must run over j. To avoid, that we
 * have to seek the j-th column for all i in the graph of AH, we do not
 * pass the graph of AH, but its transpose.
 * But also by knowing, which of the entries k in Vj are valid, can avoid
 * a lot of useless work. Therefore we pass the graph of (Ah * P)^T, too.
 *
 **********************************************************/

template <typename T>
void TransferOperator<T>::
GalerkinProduct(       CRS_Matrix<T>**           a_H,
                 const CRS_Matrix<T>&            a_h,
                 const DependencyGraph<T>& graph_AHT,
                 const DependencyGraph<T>& graph_VT  ) const
{
    ENTER_FCN("TransferOperator::GalerkinProduct");
    
#ifdef PROFILE_TRANSFEROPERATOR
    Double t1 = AMG_GET_REAL_TIME
#endif

#ifdef DEBUG_TRANSFEROPERATOR
    // check, if Ah is sorted correctly
    for( int rr = 1; rr <= a_h.GetNrows(); rr++ ) {
        // diagonal entry at leading position?
        if( a_h.GetColPointer()[a_h.GetRowPointer()[rr]] != rr ) {
            (*debug) << "ERROR: TransferOperator<T>::GalerkinProduct:"
                     << std::endl << "       A_h: in row ["<<rr<<"] "
                        "the first entry not diagonal (col ["
                     << a_h.GetColPointer()[a_h.GetRowPointer()[rr]+1]
                     << "]" << std::endl;
            Error( "TransferOperator::GalerkinProduct: non-diagonal "
                   "entry at leading position -> details in debug file",
                   __FILE__, __LINE__ );
        }
        // offdiagonal part sorted?
        for( int icc = a_h.GetRowPointer()[rr]+2;
                 icc < a_h.GetRowPointer()[rr+1]; icc++ ) {
            if( a_h.GetColPointer()[icc] < a_h.GetColPointer()[icc-1] ) {
                (*debug) << "ERROR: TransferOperator<T>::GalerkinProduct:"
                         << std::endl << "       A_h: row [" << rr << "] is "
                            "not sorted: col [" << a_h.GetColPointer()[icc]
                         << "] follows col[" << a_h.GetColPointer()[icc-1]
                         << "]" << std::endl;
                Error( "TransferOperator::GalerkinProduct: rows "
                       "of matrix A_h are not correctly sorted -> "
                       "details in debug file", __FILE__, __LINE__ );
            }
        }
    }
    // check, if the restriction matrix is sorted
    for( int rr = 1; rr <= Restriction_->GetNrows(); rr++ ) {
        for( int icc = Restriction_->GetRowPointer()[rr]+1;
                 icc < Restriction_->GetRowPointer()[rr+1]; icc++ ) {
            if( Restriction_->GetColPointer()[icc]   <
                Restriction_->GetColPointer()[icc-1]   ) {
                (*debug) << "ERROR: TransferOperator<T>::GalerkinProduct:"
                         << std::endl << "       R: row [" << rr << "] is "
                            "not sorted: col [" << Restriction_->GetColPointer()[icc]
                         << "] follows col[" << Restriction_->GetColPointer()[icc-1]
                         << "]" << std::endl;
                Error( "TransferOperator::GalerkinProduct: rows "
                       "of matrix A_h are not correctly sorted -> "
                       "details in debug file", __FILE__, __LINE__ );
            }
        }
    }
    // check, if the graph is sorted
    for( int rr = 1; rr <= graph_AHT.GetNumNodes(); rr++ ) {
        for( int ee = 0; ee < graph_AHT.GetNumEdges(rr)-1; ee++ ) {
            if( graph_AHT.GetEdges(rr)[ee] > graph_AHT.GetEdges(rr)[ee+1] ) {
                (*debug) << "ERROR: TransferOperator<T>::GalerkinProduct:"
                         << std::endl << "       edges ["<<ee<<"] and ["
                         << (ee+1)<<"] of node ["<<rr<<"] not sorted\n";
                Error( "TransferOperator::GalerkinProduct: graph "
                       "not sorted", __FILE__, __LINE__ );
                                     
            }
        }
    }
#endif

  ///////////////////////////////////////
  // prepare temporary line buffer of V^T
    T *V       = NULL; // array containing the data of V^TT (j-th row)
    Integer  lengthV = 0;    // maximal edges at a node in graph_VT
    // set lengthV as maximal number of edges over all nodes in graph_VT
    for( Integer i = 1; i <= graph_VT.GetNumNodes(); i++ ) {
        if( lengthV < graph_VT.GetNumEdges(i) ) {
            lengthV = graph_VT.GetNumEdges( i );
        }
    }
    // create the array, but zero-based indexed (this matches with
    // the also zero-based return value of graph_VT.GetEdges(..))
    NewArray( V, T, lengthV );  V++;

 ///////////////////////////////////
 // create object for coarse matrix
    // first calc the line lengths of AH (note that we have the graph of AH^T)
    Integer *RowSizeAH = NULL;
    NewArray( RowSizeAH, Integer, graph_AHT.GetNumNodes() );
    for( Integer i = 1; i <= graph_AHT.GetNumNodes(); i++ )  RowSizeAH[i] = 0;
    for( Integer i = 1; i <= graph_AHT.GetNumNodes(); i++ ) {
        for( Integer ij = 0; ij < graph_AHT.GetNumEdges(i); ij++ ) {
            RowSizeAH[graph_AHT.GetEdges(i)[ij]]++;
        }
    }
    // create matrix object
    if( *a_H )  delete *a_H;
    *a_H = New CRS_Matrix<T>( graph_AHT.GetNumNodes(),   // # rows
                              graph_AHT.GetNumNodes(),   // # columns
                              graph_AHT.GetNumEdges() ); // # non-zeros
    // set row lengths of the coarse matrix
    for( Integer i = 1; i <= graph_AHT.GetNumNodes(); i++ ) {
        (*a_H)->SetRowSize( i, RowSizeAH[i] );
        // Reset the row size, because from now on we will use this
        // array to count the already set entries in the rows of AH.
        // Set to 1, because it should point to the first offdiagonal
        // entry, since the diagonal entries are set separatly.
        RowSizeAH[i] = 1;
    }

 ////////////////////////////////////////////
 // get some pointers for performant access
    const Integer* const pRow_Ah = a_h.GetRowPointer();
    const Integer* const pCol_Ah = a_h.GetColPointer();
    const T* const pDat_Ah = a_h.GetDataPointer();
    // pointers to the restriction matrix
    const Integer* const pRow_R = Restriction_->GetRowPointer();
    const Integer* const pCol_R = Restriction_->GetColPointer();
    const T* const pDat_R = Restriction_->GetDataPointer();
    // pointers to the coarse matrix
    Integer* const pRow_AH = (Integer* const)(*a_H)->GetRowPointer();
    Integer* const pCol_AH = (Integer* const)(*a_H)->GetColPointer();
    T* const pDat_AH = (T* const)(*a_H)->GetDataPointer();

 //////////////////////////
 // build Galerkin product
    // Nomenclatur of variables and indices:
    //   i,j,k,l : according to the comment given above this function
    //   ix      : running index, belonging to the index x (i,j,k,l)
    //   ixY     : running index, belonging to the index x (i,j,k,l)
    //             in the matrix Y 
    //   V       : the temporary vector Vj
    //   kV      : pointer to the column indices of Vj, simply points
    //             to the edges of node j in the graph of V^T
    // 

    // the outer loop runs over the columns of a_H
    for( Integer j = 1; j <= (*a_H)->GetNcols(); j++ ) {
        // build temporary vector V. Note that we need not store the
        // column indices or the number of its entries, since these
        // data are already contained in graph_VT. Nevertheless we store
        // the vector's size, because is heavily used as loop boundary.
        const Integer        sizeV = graph_VT.GetNumEdges(j);
        const Integer* const    kV = graph_VT.GetEdges(j);
        // Vj[k] = sum over l for fix k, with (j,k) in graph_VT
        for( Integer ik = 0; ik < sizeV; ik++ ) {
            const Integer k = kV[ik];
            Integer ilA = pRow_Ah[k]+1, // offdiagonals of k-th row in Ah
                    ilR = pRow_R[j];    // j-th row in R
            V[ik] = 0.0;
            while( ilA < pRow_Ah[k+1] && ilR < pRow_R[j+1] ) {
                // l == k => R_jl must be multiplied with the diagonal
                //           entry A_kk, of which the position is special
                if( k == pCol_R[ilR] ) {
                    V[ik] += pDat_Ah[pRow_Ah[k]] * pDat_R[ilR++];
                // the column indices match => multiply A_kl * R_jl
                } else if( pCol_Ah[ilA] == pCol_R[ilR] ) {
                    V[ik] += pDat_Ah[ilA++] * pDat_R[ilR++];
                // if one of the indices "l" is smaller than the other,
                // increment its running index
                } else if( pCol_Ah[ilA] < pCol_R[ilR] ) {
                    ilA++;
                // <==> ( pCol_Ah[ilA] > pCol_R[ilR] )
                } else {
                    ilR++;
                }
            }
            // if we stopped due to the end of the row of A_h, we
            // might have forgotten the diagonal element -> further
            // search in the rest of the SORTED j-th row of R.
            while( ilR < pRow_R[j+1] && pCol_R[ilR] <= k ) {
                if( pCol_R[ilR] == k ) {
                    V[ik] += pDat_Ah[pRow_Ah[k]] * pDat_R[ilR];
                    break;
                } else {
                    ilR++;
                }
            }
        }
        // finished building auxiliary vector Vj, now multiply R_i * Vj
        // for all rows i in a_H. (here we need not care about diagonal
        // entries, because they neither are stored at a special position
        // in R, nor in V)
        for( Integer ii = 0; ii < graph_AHT.GetNumEdges(j); ii++ ) {
            const Integer i = graph_AHT.GetEdges(j)[ii];
            Integer ikV = 0,
                    ikR = pRow_R[i];
            T AH_ij = 0.0;
            while( ikV < sizeV && ikR < pRow_R[i+1] ) {
                // matching indices
                if( kV[ikV] == pCol_R[ikR] ) {
                    AH_ij += V[ikV++] * pDat_R[ikR++];
                } else if( kV[ikV] < pCol_R[ikR] ) {
                    ikV++;
                } else {
                    ikR++;
                }
            }
            // set column index, diagonal entries must be treated specially
            if( i == j ) {
                pDat_AH[pRow_AH[i]] = AH_ij;
                pCol_AH[pRow_AH[i]] = j;
            } else {
                pDat_AH[pRow_AH[i]+RowSizeAH[i]] = AH_ij;
                pCol_AH[pRow_AH[i]+RowSizeAH[i]] = j;
                RowSizeAH[i]++;
            }
        }
    } // j

    DeleteArray( RowSizeAH );
    V--;  DeleteArray( V );

#ifdef PROFILE_TRANSFEROPERATOR
    Double t2 = AMG_GET_REAL_TIME
    (*cla) << " AMG: Galerkin product: "<<(t2-t1)<<" seconds\n";
#endif
}

/**********************************************************/

template <typename T>
void TransferOperator<T>::
GalerkinProduct(       CRS_Matrix<T>** a_H,
                 const CRS_Matrix<T>&  a_h  ) const
{
#ifdef DEBUG_TRANSFEROPERATOR
    if( Prolongation_ == NULL ) {
        Error( "TransferOperator::GalerkinProduct: this variant of "
               "GalerkingProduct needs the interpolation matrix built "
               "explicitly\n", __FILE__, __LINE__ );
    }
#endif

#ifdef PROFILE_TRANSFEROPERATOR
    Double t1 = AMG_GET_REAL_TIME
#endif

    // get some pointers for performant access
    const Integer* const pRowA = a_h.GetRowPointer();
    const Integer* const pColA = a_h.GetColPointer();
    const T* const pDatA = a_h.GetDataPointer();
    // pointers for the restriction matrix
    const Integer* const pRowR = Restriction_->GetRowPointer();
    const Integer* const pColR = Restriction_->GetColPointer();
    const T* const pDatR = Restriction_->GetDataPointer();
    // pointers for the interpolation matrix
    const Integer* const pRowP = Prolongation_->GetRowPointer();
    const Integer* const pColP = Prolongation_->GetColPointer();
    const T* const pDatP = Prolongation_->GetDataPointer();
    // pointer to the splitting and fine-coarse mapping

    // create the object for the pre-matrix
    PreMatrix<T> preAH( Restriction_->GetNrows() );

    /////////////////////////////////////////////////
    //              n                              //
    //             ---                             //
    // AH[iC,jC] = >   R[iC,k] * AH[k,l] * P[l,jC] //
    //             ---                             //
    //           k,l = 1                           //
    /////////////////////////////////////////////////

    // for i in [1,nH]
    for( Integer iC = 1; iC <= Restriction_->GetNrows(); iC++ ) {
        // for k in [1,nh] (only non-zero entries)
        for( Integer ik = pRowR[iC]; ik < pRowR[iC+1]; ik++ ) {
            const Integer k     = pColR[ik];
            const T R_iCk = pDatR[ik];
            // for l in [i,nh] (only non-zero entries)
            for( Integer il = pRowA[k]; il < pRowA[k+1]; il++ ) {
                const Integer l          = pColA[il];
                // R_iCk_A_kl = R[iC][k] * Ah[k][l]
                const T R_iCk_A_kl = R_iCk * pDatA[il];
                for( Integer ijC = pRowP[l]; ijC < pRowP[l+1]; ijC++ ) {
                    // AH[iC][jC] += R[iC][k] * Ah[k][l] * P[l][jC]
                    preAH.AddToEntry( iC, pColP[ijC], R_iCk_A_kl * pDatP[ijC] );
                }
            }
        }
    }

    // For developing purposes it is possible to import an external
    // Galerkin product and overwrite the calculated.
#ifdef TRANSFER_OPERATOR_IMPORT_GALERKIN_PRODUCT
    char filename[500];
    sprintf( filename,
             TRANSFER_OPERATOR_GALERKIN_PRODUCT_IMPORT_PREFIX "A_%d",
             Restriction_->GetNrows() );
    std::cout
    << "\033[37m importing Galerkin product from \033[0;1m\""
    << filename<<"\"\033[0m\n";
    preAH.ImportASCII( filename );
#endif // TRANSFER_OPERATOR_IMPORT_GALERKIN_PRODUCT

    // sort the prematrix
    // NOTE: If we need once the resulting matrix to conform all layout
    //       specifications (additionally sorted offdiagonal entries),
    //       call Sort() instead of SortDiag() and set the matrixes
    //       layout flag by calling SetLayoutFlag().
    preAH.SortDiagonal();
    // convert the prematrix into the CRS format
    *a_H = New CRS_Matrix<T>( Restriction_->GetNrows(),
                              Restriction_->GetNrows(),
                              preAH.GetNumNonzeros()   );
    for( Integer i = 1, position = 1; i <= Restriction_->GetNrows(); i++ ) {
        const Integer rowlength = preAH.GetRowSize( i );
        // copy column indices
        memcpy( (*a_H)->GetColPointer() + position,
                preAH.GetRowCols(i)+1, rowlength * sizeof(Integer) );
        // copy data
        memcpy( (*a_H)->GetDataPointer() + position,
                preAH.GetRowData(i)+1, rowlength * sizeof(T) );
        // set row length
        (*a_H)->SetRowSize( i, rowlength );
        // update the position
        position += rowlength;
    }

#ifdef PROFILE_TRANSFEROPERATOR
    Double t2 = AMG_GET_REAL_TIME
    (*cla) << " AMG: Galerkin product: "<<(t2-t1)<<" seconds\n";
#endif

}

/**********************************************************/
#ifdef TRANSFER_OPERATOR_IMPORT_INTERPOLATION

template <typename T>
void TransferOperator<T>::ImportInterpolation( const Integer size_h,
                                               const Integer size_H )
{
    ENTER_FCN("TransferOperator::ImportInterpolation");

    char filename[500];
    sprintf( filename,
             TRANSFER_OPERATOR_INTERPOLATION_IMPORT_PREFIX "P_%d",
             size_h );
    std::cout
    << "\033[37m importing interpolation from \033[0;1m\""
    << filename<<"\"\033[0m\n";
    
    PreMatrix<T> preP, preR;
    // read interpolation matrix
    if( ! preP.ImportASCII(filename) ) {
        Error( "TransferOperator::ImportInterpolation: could not import "
               "interpolation matrix.", __FILE__, __LINE__ );
    }
    // If the interpolation is imported, we create the restriction
    // in any case.
    preR.SetNumRows( size_H );
    for( Integer i = 1; i <= preP.GetNrows(); i++ ) {
        for( Integer ij = 1; ij <= preP.GetRowSize(i); ij++ ) {
            preR.SetEntry( preP.GetRowCols(i)[ij], i,
                           preP.GetRowData(i)[ij]     );
        }
    }
    preP.Sort();
    preR.Sort();
    // convert the pre-matrices
    delete Prolongation_;
    delete Restriction_;
    Prolongation_ = New CRS_Matrix<T>( size_h, size_H, preP );
    Restriction_  = New CRS_Matrix<T>( size_H, size_h, preR );
    if( !Prolongation_ || !Restriction_ ) {
        Error( "TransferOperator::ImportInterpolation: could not convert "
               "imported pre-matrices into CRS format", __FILE__, __LINE__ );
    }
}

#endif // TRANSFER_OPERATOR_IMPORT_INTERPOLATION
/**********************************************************/
#ifdef DEBUG_TRANSFEROPERATOR

template <typename T>
bool TransferOperator<T>::CheckOperators() const
{
    bool transpositionOK = true,
         rowSumOK        = true;

    if( Prolongation_ && Restriction_ ) {
        // get some pointers for performant and dirty (CAST!) access
        const Integer* const pRowR = Restriction_->GetRowPointer();
        const Integer* const pColR = Restriction_->GetColPointer();
        const T* const pDatR = Restriction_->GetDataPointer();
        // pointers for interpolation
        const Integer* const pRowP = Prolongation_->GetRowPointer();
        const Integer* const pColP = Prolongation_->GetColPointer();
        const T* const pDatP = Prolongation_->GetDataPointer();

        ////////////////////////////////////////////////////
        // check, whether P is really the transpose of R
        ////////////////////////////////////////////////////
        // inclusion R in P
        for( Integer i = 1; i <= Restriction_->GetNrows(); i++ ) {
            for( Integer ijR = pRowR[i]; ijR < pRowR[i+1]; ijR++ ) {
                const Integer j   = pColR[ijR];
                const T Rij = pDatR[ijR];
                      Integer ijP = 0;
                for( ijP = pRowP[j]; ijP < pRowP[j+1]; ijP++ ) {
                    if( pColP[ijP] == i ) {
                        if( pDatP[ijP] != Rij ) {
                            (*debug)
                            << "ERROR: TransferOperator::CheckOperators:\n"
                               " -> R("<<i<<","<<j<<")["<<ijR<<"] = "<<Rij
                            <<" and P("<<j<<","<<i<<")["<<ijP<<"] = "<<pDatP[ijP]
                            << " differ."<<std::endl;
                            transpositionOK = false;
                        }
                        break;
                    }
                }
                if( ijP >= pRowP[j+1] ) {
                    (*debug)
                    << "ERROR: TransferOperator::CheckOperators:\n -> "
                       "R("<<i<<","<<j<<")["<<ijR<<"] = "<<Rij<<", but "
                       "P("<<j<<","<<i<<") == 0 (not existing)"
                    << std::endl;
                    transpositionOK = false;
                }
            }
        }
        // inclusion P in R
        for( Integer i = 1; i <= Prolongation_->GetNrows(); i++ ) {
            for( Integer ijP = pRowP[i]; ijP < pRowP[i+1]; ijP++ ) {
                const Integer j   = pColP[ijP];
                const T Pij = pDatP[ijP];
                      Integer ijR = 0;
                for( ijR = pRowR[j]; ijR < pRowR[j+1]; ijR++ ) {
                    // The corresponding entries are present in both matrices.
                    // This means that we have already compared them.
                    if( pColR[ijR] == i )  break;
                }
                if( ijR >= pRowR[j+1] ) {
                    (*debug)
                    << "ERROR: TransferOperator::CheckOperators:\n -> "
                       "P("<<i<<","<<j<<")["<<ijP<<"] = "<<Pij<<", but "
                       "R("<<j<<","<<i<<") == 0"<<std::endl;
                    transpositionOK = false;
                }
            }
        }
        if( ! transpositionOK ) {
            Warning( "TransferOperator::CheckOperators: P != R'\n "
                     "for more details inspect the debug file",
                     __FILE__, __LINE__ );
        }

        ////////////////////////////////////////////////////
        // check the prolongation's property to interpolate
        // constant functions exactly
        ////////////////////////////////////////////////////
        if( false &&
            Prolongation_ )  {
            for( Integer i = 1; i <= Prolongation_->GetNrows(); i++ ) {
                if( pRowP[i] < pRowP[i+1] ) {
                    T sum = static_cast<T>(0.0);
                    for( Integer ij = pRowP[i]; ij < pRowP[i+1]; ij++ ) {
                        sum += pDatP[ij];
                    }
                    if( Abs(sum - static_cast<T>(1.0)) > 1e-10 ) {
                        (*debug)
                        << "ERROR: TransferOperator::CheckOperators:\n"
                           " -> sum_i { P("<<i<<",i) } = "<<sum
                        << " != 1.0"<<std::endl;
                        rowSumOK = false;
                    }
                }
            }
        }
        if( ! rowSumOK ) {
            Warning( "TransferOperator::CheckOperators: row sum of "
                     "P is not 1.0\n for more details inspect the debug "
                     "file", __FILE__, __LINE__ );
        }
    }
    return transpositionOK && rowSumOK;
}

#endif // DEBUG_TRANSFEROPERATOR
/**********************************************************/
} // namespace OLAS

/**********************************************************/
#ifdef DEBUG_TO_CERR
#undef debug
#endif // DEBUG_TO_CERR
/**********************************************************/


/**********************************************************
 * Creation of the interpolation matrix only.
 * Just for the case we will need it once
 **********************************************************

template <typename T>
bool TransferOperator<T>::
CreateOperators( const CRS_Matrix<T>& matrix,
                 const Topology<T>&   topology )
{
    ENTER_FCN("TransferOperator::CreateOperators");

#ifdef DEBUG_TRANSFEROPERATOR
    if( matrix.GetNcols() != topology.GetSizeh() ) {
        (*debug) << "ERROR: # columns of the matrix == " << matrix.GetNcols()
                 << std::endl << "       != size_h in topology   == "
                 << topology.GetSizeh() << std::endl << "must match!"
                 << std::endl;
        Error( "Dimensions of matrix and topology do not match",
                __FILE__, __LINE__ );
    }
#endif

    Reset();

    Integer  nNonZeros  = 0,
            *RowLengths = NULL;
    NewArray( RowLengths, Integer, matrix.GetNcols() );
    // get number of interpolating points for each point
    nNonZeros = topology.GetNumInterpolatingPoints( RowLengths );

    // create the prolongation matrix object
    Prolongation_ = New CRS_Matrix<T>( topology.GetSizeh(),
                                       topology.GetSizeH(),
                                       nNonZeros );

    // get some pointers for performant and dirty (CAST!) access
    Integer* const pRow = (Integer* const)Prolongation_->GetRowPointer();
    Integer* const pCol = (Integer* const)Prolongation_->GetColPointer();
    T* const pDat = (T* const)Prolongation_->GetDataPointer();
    // pointer to the splitting and fine-coarse mapping
    const Integer* const CoarseIndex = topology.GetCoarseFineSplitting();
    Integer RowStart = 1;

    // fill prolongation row by row
    for( Integer i = 1; i <= topology.GetSizeh(); i++ ) {
        // set row start
        pRow[i] = RowStart;
        // set interpolation weights
        if( CoarseIndex[i] >= Topology<T>::COARSE ) {
            pCol[RowStart] = CoarseIndex[i];
            pDat[RowStart] = (T)1.0;
        } else {
#ifdef DEBUG_TRANSFEROPERATOR
            if( RowLengths[i] == 0 ) {
                (*debug) << "ERROR: F-point [" << i << "] does not have any "
                            "interpolating C-neighbours" << std::endl;
                Error( "TransferOperator::CreateOperators: found F-point "
                       "without an interpolating C-neighbour" );
            }
#endif
            const T              weight = opType<T>::invert( (T)RowLengths[i] );
            const Integer        nEdges = topology.GetS().GetNumEdges( i );
            const Integer* const  Edges = topology.GetS().GetEdges( i );
            for( Integer ij = 0, iedge = RowStart; ij < nEdges; ij++ ) {
                // i.e. point i is coarse point
                if( CoarseIndex[Edges[ij]] >= Topology<T>::COARSE ) {
#ifdef DEBUG_TRANSFEROPERATOR
                    if( iedge >= RowStart + RowLengths[i] ) {
                        (*debug) << "ERROR: TransferOperator<T>::CreateOperators:"
                                 << std::endl << "number of weights exceeds row "
                                    "length " << RowLengths[i] << " in row " << i
                                 << std::endl;
                        Warning( "TransferOperator::CreateOperators: "
                                 "-> see debug file" );
                    }
#endif
                    pCol[iedge] = CoarseIndex[Edges[ij]];
                    pDat[iedge] = weight;
                    iedge++;
                }
            }
        }
        // increment row start for the next row
        RowStart += RowLengths[i];
    }
    pRow[topology.GetSizeh()+1] = RowStart;

    DeleteArray( RowLengths );
    
    return true;
}

 **********************************************************/
