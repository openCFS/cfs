#include "OLAS/multigrid/transfer.hh"
#include "OLAS/multigrid/hierarchylevel.hh"
#include "OLAS/multigrid/amg.hh"
#include "OLAS/multigrid/prematrix.hh"

#include <assert.h>
#include <limits>

namespace CoupledField {

/**********************************************************/

template <typename T>
TransferOperator<T>::TransferOperator()
                   : Prolongation_(NULL),
                     Restriction_(NULL),
                     AuxProlongation_(NULL),
                     AuxRestriction_(NULL)
{
}

/**********************************************************/

template <typename T>
TransferOperator<T>::~TransferOperator()
{
    Reset();
}

/**********************************************************/

template <typename T>
void TransferOperator<T>::Reset()
{

    delete Prolongation_;  Prolongation_ = NULL;
    delete Restriction_;   Restriction_  = NULL;
    delete AuxProlongation_;  AuxProlongation_ = NULL;
    delete AuxRestriction_;   AuxRestriction_  = NULL;
}


template <typename T>
bool TransferOperator<T>::
CreateOperators( const CRS_Matrix<T>& sysMatrix,
                 const CRS_Matrix<Double>& auxMatrix,
                 const Topology<Double>&   topology,
                 const Agglomerate<Double>& agglomerates,
                 const AMGInterpolationType itype,
                 const AMGType amgType,
                 const bool build_interpolation )
{
    Reset();

    ////////////////////////////////////////////////////////
    // create specific interpolation
    ////////////////////////////////////////////////////////

    bool creationResult = false;

    switch( itype ) {
        // constant interpolation
        // choose amg-type
        case AMG_INTERPOLATION_CONSTANT:
          switch(amgType){
          case SCALAR:
            creationResult = CreateOperatorsConstantScalar( sysMatrix, auxMatrix, topology );
          break;

          case VECTORIAL:
            creationResult = CreateOperatorsConstantVectorial( sysMatrix, auxMatrix, topology );
          break;

          case EDGE:
            EXCEPTION("transfer.cc: edge-version shouldn't be handled here!");
          break;
          }
        break;
        // simple weighted average interpolation
        case AMG_INTERPOLATION_SIMPLE_WEIGHTED:
          EXCEPTION("AMG: simple-weighted interpolation not yet implemented");
            break;
        //
        case AMG_INTERPOLATION_SMOOTHED_SCALING:
          EXCEPTION("AMG: smoothed scaling interpolation not yet implemented");
            break;
        // developing version of interpolation
        case AMG_INTERPOLATION_DEVELOP:
            EXCEPTION("AMG: develop interpolation not yet implemented");
            break;

        default:
            EXCEPTION("topology.cc: TransferOperator::CreateOperators: non "
                   "supported constant %d as interpolation type\n");
            break;
    }


    return creationResult;

}


/************************************************************************/

template <typename T>
void TransferOperator<T>::GetCoarseEdgeIndNode(StdVector< StdVector< Integer> >& r)
{
  r = cEdgeIndNode_;
}
/************************************************************************/

template <typename T>
bool TransferOperator<T>::
CreateProlongationOperatorEdgeAux(const CRS_Matrix<Double>& auxMatrix,
                       const Agglomerate<Double>& agglomerates)
{
    Reset();


    Integer  nNonZeros  = 0;
    // get number of interpolated points for each C-point
    nNonZeros = agglomerates.GetNumInterpolatedPoints();

    //sanity check
//TODO temporarily deactivated, sometimes nNonZeros is 2 higher than GetSizeh...CHECK THIS
    //if((UInt)nNonZeros != agglomerates.GetSizeh()){
    //  EXCEPTION("transfer.cc: number of points in agglomerates doesn't match size of all nodes!! \n"
    //      "Implementation error");
    //}


    // create the restriction and prolongation matrix object for the AUXILIARY-matrix
    AuxRestriction_ = new CRS_Matrix<Double>( agglomerates.GetSizeH(), agglomerates.GetSizeh(), nNonZeros );
    // we can not set the sys-prolongation matrix because
    // we don't know the nnz entries yet ...

    // row pointers for R and P
    StdVector<UInt> pRowR;
    pRowR.Resize( agglomerates.GetSizeH() + 1, 0 );
    // column pointers for R and P
    StdVector<UInt> pColR;
    pColR.Resize( nNonZeros, 0 );
    // data pointer for R and P
    StdVector<Double> pDatR;
    pDatR.Resize( nNonZeros, 0.0 );


    // pointer to the splitting and fine-coarse mapping
    const StdVector<Integer> CoarseIndex = agglomerates.GetCoarseFineSplitting();

    /******************** Version P = R^T *********************/
    /* so we only calculate the restriction matrix
     * and define the prolongation as the transpose restriction
     */

    // fill row pointer
    nNonZeros = 0;
    for(UInt i  = 0; i < agglomerates.GetSizeh() ; ++i){
      if(agglomerates.IsCPoint(i)){
        // cast is safe
        UInt cI = (UInt)agglomerates.GetCoarseIndex(i);
        StdVector<Integer> t = agglomerates.GetAgglomerateOfNode(i);
        nNonZeros += t.GetSize();
        pRowR[cI + 1] = nNonZeros;
      }
    }

    // the following algorithm constructs R
    for( Integer i = 0; i < (Integer)agglomerates.GetSizeh(); ++i ) {
      if( agglomerates.IsCPoint(i) ){
        // get the coarse index of node i, cast is safe, it already is a C-point
        // this is now the row index
        UInt cI = (UInt)agglomerates.GetCoarseIndex(i);
        //get the agglomerate, point i belongs to
        StdVector<Integer> agg = agglomerates.GetAgglomerateOfNode(i);

        //sanity check
        if( (Integer)cI != agglomerates.GetCoarseIndex(agg[0]) ) EXCEPTION("transfer.cc: Coarse-index problem");
        if( agg.GetSize() != (pRowR[cI + 1] - pRowR[cI])) EXCEPTION("transfer.cc: Coarse-index problem");
        // every index in the agglomerate is now assigned a 1
        UInt it = 0;
        for(UInt j = pRowR[cI]; j < pRowR[cI + 1]; ++j){
          pColR[j] = agg[it];
          pDatR[j] = 1.0;
          it++;
        }
      }
    }

    AuxRestriction_->SetSparsityPatternData( pRowR, pColR, pDatR );

    return true;
}





/************************************************************************/

template <typename T>
bool TransferOperator<T>::
CreateProlongationOperatorEdgeSys(const CRS_Matrix<T>& SysMatrix,
                        const StdVector< StdVector< Integer> >& edgeIndNode,
                        const StdVector< Integer>& nodeNumIndex,
                        const Agglomerate<Double>& A)
{
  // sizeH is is size of the system matrix...numRows of prolongation matrix
  UInt sizeH = A.GetNumCoarseEdges();
  UInt sizeh = edgeIndNode.GetSize();
  cEdgeIndNode_.Resize(sizeH);
  // row pointers for R and P
  StdVector<UInt> pRowR;
  pRowR.Resize( sizeH + 1, 0 );
  // column pointers for R and P
  StdVector<UInt> pColR;
  // data pointer for R and P
  StdVector<Double> pDatR;

  StdVector< StdVector<Integer> > eIndRNode (sizeh);
  StdVector< StdVector<Integer> > eAggNum (sizeh);
  // just a performance improvement...outsorced from the O(~N^2) loop
  StdVector<Integer> eIN (edgeIndNode.GetSize() * 2);
  for(UInt j = 0; j < sizeh; ++j){
	  const StdVector<Integer>& e = edgeIndNode[j];

	  StdVector<Integer> & etmp = eIndRNode[j];
	  StdVector<Integer> & eA = eAggNum[j];
	  etmp.Resize(2);
	  eA.Resize(2);
	  etmp[0] = A.GetIndexOfRealNode(e[1]);
	  etmp[1] = A.GetIndexOfRealNode(e[0]);
	  eA[0] = A.GetAgglomerateNumOfNode(etmp[0]);
	  eA[1] = A.GetAgglomerateNumOfNode(etmp[1]);
    if( (eA[0] == -1) || (eA[1] == -1) ){
      EXCEPTION("coarse index problem");
    }
  }

  UInt count = 0;
  StdVector<Integer> cEindSwap (2);
  for(UInt i = 0; i < sizeH; ++i){
    // get edge i, nodes of these edges are coarse-nodes
    const StdVector<Integer>& cEind = A.GetCoarseEdgeInd(i);
    const StdVector<Integer>& cEreal = A.GetCoarseEdgeReal(i);
    StdVector<Integer>& cE = cEdgeIndNode_[i];
    cE.Resize(2,0);
    // now swap the nodes and search again
    cEindSwap = cEind;
    cEindSwap.Swap(0, 1);
    for(UInt j = 0; j < sizeh; ++j){
      if( eAggNum[j] == cEind ){
        pColR.Push_back(j);
        pDatR.Push_back(1.0);
        cE[0] = cEreal[0];
        cE[1] = cEreal[1];
        count++;
      }
      if( eAggNum[j] == cEindSwap ){
        pColR.Push_back(j);
        pDatR.Push_back(-1.0);
        cE[0] = cEreal[1];
        cE[1] = cEreal[0];
        count++;
      }
    }
    pRowR[i + 1] = count;
  }
  Restriction_ = new CRS_Matrix<Double>( sizeH, sizeh, count );
  Restriction_->SetSparsityPatternData( pRowR, pColR, pDatR );

  return true;
}


/***************************************************************************/
template <typename T>
bool TransferOperator<T>::
CreateOperatorsConstantScalar( const CRS_Matrix<T>& sysMatrix,
                               const CRS_Matrix<Double>& auxMatrix,
                               const Topology<Double>&   topology )
{

  Integer  nNonZeros  = 0;
  StdVector<UInt> RowLengths;
  // get number of interpolated points for each C-point
  nNonZeros = topology.GetNumInterpolatedPoints( RowLengths );
  // create the restriction and prolongation matrix object
  Prolongation_ = new CRS_Matrix<Double>( topology.GetSizeh(), topology.GetSizeH(), nNonZeros );
  AuxProlongation_ = new CRS_Matrix<Double>( topology.GetSizeh(), topology.GetSizeH(), nNonZeros );


  // row pointers for R and P
  StdVector<UInt> pRowP;
  pRowP.Resize( topology.GetSizeh() + 1, 0 );
  // column pointers for R and P
  StdVector<UInt> pColP;
  pColP.Resize( nNonZeros, 0 );
  // data pointer for R and P
  StdVector<Double> pDatP;
  pDatP.Resize( nNonZeros, 0.0 );


  // pointer to the splitting and fine-coarse mapping
  const Integer* const CoarseIndex = topology.GetCoarseFineSplitting();

  /********************* Version P = R^T **********************/
  /* so we only calculate the prolongation (interpolation) matrix
   * and define the restriction as the transpose prolongation
   */

  nNonZeros = 0;
  for( Integer i = 0; i < (Integer)topology.GetSizeh(); i++ ) {
    // if point [i] is C-point, in column [i] there is only one
    // entry, since C-points are simply injected with weight 1.0
    const Integer iC = CoarseIndex[i];
    if( iC >= Topology<Double>::COARSE ) {
      pColP[nNonZeros  ] = iC;
      pDatP[nNonZeros] = 1.0;
      nNonZeros++;
    } else {
      Integer   denominator = 0;
      const Integer        nEdges = topology.GetS().GetNumEdges( i );
      const Integer* const  Edges = topology.GetS().GetEdges( i );
      // build interpolation for row [i] in P
      for( Integer ij = 0; ij < nEdges; ij++ ) {
        const Integer  j = Edges[ij],
            jC = CoarseIndex[j];
        if( jC >= Topology<Double>::COARSE ) {
          pColP[nNonZeros] = jC;
          nNonZeros++;
          denominator++;
        }
      }
      if( denominator ) {
        const Double weight = OpType<Double>::invert(static_cast<Double>(denominator));

        for( Integer ij = pRowP[i]; ij < nNonZeros; ij++ ) {
          // set the interpolation weight ...
          pDatP[ij] = weight;
        }
      } // if( denominator )
    }
    // terminate row i of P by start index of the next one
    pRowP[i+1] = nNonZeros;
  }

  Prolongation_->SetSparsityPatternData( pRowP, pColP, pDatP );
  AuxProlongation_->SetSparsityPatternData( pRowP, pColP, pDatP );


  return true;
}

/***************************************************************************/

template <typename T>
bool TransferOperator<T>::
CreateOperatorsConstantVectorial( const CRS_Matrix<T>& sysMatrix,
                                  const CRS_Matrix<Double>& auxMatrix,
                                  const Topology<Double>&   topology )
{
  Integer  nNonZeros  = 0;
  StdVector<UInt> RowLengths;
  // get number of interpolated points for each C-point
  nNonZeros = topology.GetNumInterpolatedPoints( RowLengths );
  // create the restriction and prolongation matrix object for the AUXILIARY-matrix
  AuxProlongation_ = new CRS_Matrix<Double>( topology.GetSizeh(), topology.GetSizeH(), nNonZeros );
  // we can not set the sys-prolongation matrix because
  // we don't know the nnz entries yet ...

  // row pointers for R and P
  StdVector<UInt> pRowP;
  pRowP.Resize( topology.GetSizeh() + 1, 0 );
  // column pointers for R and P
  StdVector<UInt> pColP;
  pColP.Resize( nNonZeros, 0 );
  // data pointer for R and P
  StdVector<Double> pDatP;
  pDatP.Resize( nNonZeros, 0.0 );


  // pointer to the splitting and fine-coarse mapping
  const Integer* const CoarseIndex = topology.GetCoarseFineSplitting();

  /********************* Version P = R^T **********************/
  /* so we only calculate the prolongation (interpolation) matrix
   * and define the restriction as the transpose prolongation
   */
  nNonZeros = 0;
  for( Integer i = 0; i < (Integer)topology.GetSizeh(); i++ ) {
    // if point [i] is C-point, in column [i] there is only one
    // entry, since C-points are simply injected with weight 1.0
    const Integer iC = CoarseIndex[i];
    if( iC >= Topology<Double>::COARSE ) {
      pColP[nNonZeros  ] = iC;
      pDatP[nNonZeros] = 1.0;
      nNonZeros++;
    } else {
      Integer   denominator = 0;
      const Integer        nEdges = topology.GetS().GetNumEdges( i );
      const Integer* const  Edges = topology.GetS().GetEdges( i );
      // build interpolation for row [i] in P
      for( Integer ij = 0; ij < nEdges; ij++ ) {
        const Integer  j = Edges[ij],
            jC = CoarseIndex[j];
        if( jC >= Topology<Double>::COARSE ) {
          pColP[nNonZeros] = jC;
          nNonZeros++;
          denominator++;
        }
      }
      if( denominator ) {
        const Double weight = OpType<Double>::invert(static_cast<Double>(denominator));

        for( Integer ij = pRowP[i]; ij < nNonZeros; ij++ ) {
          // set the interpolation weight ...
          pDatP[ij] = weight;
        }
      } // if( denominator )
    }
    // terminate row i of P by start index of the next one
    pRowP[i+1] = nNonZeros;
  }


  AuxProlongation_->SetSparsityPatternData( pRowP, pColP, pDatP );


  /*********** Now build the system prolongation matrix **********/
  // Example
  /* P for AUXILIARY MATRIX
   *                 _    _
   *                | 1 0  |
   *                | 1 0  |
   *                |_0 1 _|
   *
   *                becomes
   *P for SYSTEM MATRIX (for 2D)
   *                _        _
   *               |  1 0 0 0 |
   *               |  0 1 0 0 |
   *               |  1 0 0 0 |
   *               |  0 1 0 0 |
   *               |  0 0 1 0 |
   *               |_ 0 0 0 1_|
   */

  UInt sizeSys = sysMatrix.GetNumRows();
  UInt sizeAux = auxMatrix.GetNumRows();

  //find out dimension of vector dof's
  //sanity check
  if(sizeSys % sizeAux != 0) EXCEPTION("AMG: Ah and Bh have wrong dimension!");
  UInt dim = sizeSys / sizeAux;

  // row pointers for P_sys
  StdVector<UInt> pRowPsys;
  pRowPsys.Resize( dim * topology.GetSizeh() + 1, 0 );
  // column pointers for P_sys
  StdVector<UInt> pColPsys;
  pColPsys.Resize( dim * nNonZeros, 0 );
  // data pointer for P_sys
  StdVector<Double> pDatPsys;
  pDatPsys.Resize( dim * nNonZeros, 0.0 );

  //over every row of P-aux
  UInt cP = 0;
  UInt rowS = 0;
  for(UInt i = 0; i < sizeAux; ++i){
    //every entry in aux corresponds to a dim x dim entry in sys
    UInt r = 0;
    //rowS = rowS + pRowP[i+1] - pRowP[i];
    for(UInt rIter = dim * i; rIter < dim * (i+1); ++rIter){

      for(UInt j = pRowP[i]; j < pRowP[i+1]; ++j){
        UInt pSysCol = dim * pColP[j] + r;
        pColPsys[ cP ] = pSysCol;
        pDatPsys[ cP ] = pDatP[j];
        cP++;
        rowS++;
      }
      pRowPsys[dim * i+r + 1] = rowS;
      r++;
    }
  }

  Prolongation_ = new CRS_Matrix<Double>( dim * topology.GetSizeh(),
                                        dim * topology.GetSizeH(),
                                        dim * nNonZeros );


  Prolongation_->SetSparsityPatternData(pRowPsys, pColPsys, pDatPsys);
  return true;
}

/***************************************************************************/



template <>
inline void TransferOperator<Double>::Prolongate( const Vector<Double>& v_H,
                                            Vector<Double>& v_h,
                                      const bool       add  ) const
{

    if( Prolongation_ != NULL ) {
        if( add )  Prolongation_->MultAdd( v_H, v_h );
        else       Prolongation_->Mult( v_H, v_h );
    } else {
        if( add )  Restriction_->MultTAdd( v_H, v_h );
        else       Restriction_->MultT( v_H, v_h );
    }

}

template <>
inline void TransferOperator<Complex>::Prolongate( const Vector<Complex>& v_H,
                                            Vector<Complex>& v_h,
                                      const bool       add  ) const
{

    if( Prolongation_ != NULL ) {
        if( add )  Prolongation_->MultAdd_type( v_H, v_h );
        else       Prolongation_->Mult_type( v_H, v_h );
    } else {
        if( add )  Restriction_->MultTAdd_type( v_H, v_h );
        else       Restriction_->MultT_type( v_H, v_h );
    }

}
/**********************************************************/
template <>
inline void TransferOperator<Double>::Restrict( const Vector<Double>& v_h,
                                          Vector<Double>& v_H,
                                    const bool       add  ) const
{
  if( Restriction_ != NULL ) {
      if( add )  Restriction_->MultAdd( v_h, v_H );
      else       Restriction_->Mult( v_h, v_H );
  } else {
      if( add )  Prolongation_->MultTAdd( v_h, v_H );
      else       Prolongation_->MultT( v_h, v_H );
  }
}

template <>
inline void TransferOperator<Complex>::Restrict( const Vector<Complex>& v_h,
                                          Vector<Complex>& v_H,
                                    const bool       add  ) const
{
  if( Restriction_ != NULL ) {
      if( add )  Restriction_->MultAdd_type( v_h, v_H );
      else       Restriction_->Mult_type( v_h, v_H );
  } else {
      if( add )  Prolongation_->MultTAdd_type( v_h, v_H );
      else       Prolongation_->MultT_type( v_h, v_H );
  }
}






/* To avoid constantly repeating the part of code that checks inbound SparseBLAS functions' status,
     use macro CALL_AND_CHECK_STATUS */
#define CALL_AND_CHECK_STATUS(function, error_message) do { \
		if(function != SPARSE_STATUS_SUCCESS)             \
		{                                                 \
			printf(error_message); fflush(0);                 \
		}                                                 \
} while(0)


template <typename T>
void TransferOperator<T>::
GalerkinProduct(StdVector<UInt>& A_H_rP,
                StdVector<UInt>& A_H_cP,
                StdVector<T>& A_H_dP,
                StdVector<UInt>& B_H_rP,
                StdVector<UInt>& B_H_cP,
                StdVector<Double>& B_H_dP,
                const CRS_Matrix<T>&  A_h,
                const CRS_Matrix<Double>& B_h,
                const AMGType amgType)
{

  CRS_Matrix<T>* TMP1 = new CRS_Matrix<T>();
  TMP1->MultTriple_MKL(A_h, *Prolongation_,  A_H_rP, A_H_cP, A_H_dP, 1, false);

  CRS_Matrix<Double>* TMP2 = new CRS_Matrix<Double>();
  TMP2->MultTriple_MKL(B_h, *AuxProlongation_, B_H_rP, B_H_cP, B_H_dP, 1, false);

  delete TMP1; TMP1 = NULL;
  delete TMP2; TMP2 = NULL;

}


/***************************************************************************/
template <typename T>
void TransferOperator<T>::
GalerkinProductEdgeAux(StdVector<UInt>& B_H_rP,
                      StdVector<UInt>& B_H_cP,
                      StdVector<Double>& B_H_dP,
                      const CRS_Matrix<Double>& B_h)
{

  CRS_Matrix<Double>* TMP = new CRS_Matrix<Double>();
  TMP->MultTriple_MKL(B_h, *AuxRestriction_, B_H_rP, B_H_cP, B_H_dP, 2, false);

  delete TMP; TMP = NULL;

}


/************** Specialization for real valued system matrix *******/
template <typename T>
void TransferOperator<T>::
GalerkinProductEdgeSys(StdVector<UInt>& A_H_rP,
                      StdVector<UInt>& A_H_cP,
                      StdVector<Double>& A_H_dP,
                      const CRS_Matrix<Double>& A_h)
{

  CRS_Matrix<Double>* TMP = new CRS_Matrix<Double>();
  TMP->MultTriple_MKL(A_h, *Restriction_, A_H_rP, A_H_cP, A_H_dP, 2, false);

  delete TMP; TMP = NULL;
}

/************** Specialization for complex valued system matrix *******/
template <typename T>
void TransferOperator<T>::
GalerkinProductEdgeSys(StdVector<UInt>& A_H_rP,
                        StdVector<UInt>& A_H_cP,
                        StdVector<Complex>& A_H_dP,
                        const CRS_Matrix<Complex>& A_h)
{

  UInt nnz = A_h.GetNnz();
  const Complex *pC = A_h.GetDataPointer();
  const UInt *cPC = A_h.GetColPointer();
  const UInt *rPC = A_h.GetRowPointer();

  double *re = new double[nnz];
  double *im = new double[nnz];
  for(UInt i = 0; i < nnz; ++i){
    re[i] = pC[i].real();
    im[i] = pC[i].imag();
  }

  CRS_Matrix<Double>* R = new CRS_Matrix<Double>(A_h.GetNumRows(), A_h.GetNumCols(), A_h.GetNnz(), rPC, cPC, re);
  CRS_Matrix<Double>* I = new CRS_Matrix<Double>(A_h.GetNumRows(), A_h.GetNumCols(), A_h.GetNnz(), rPC, cPC, im);

  CRS_Matrix<Double>* TMP_R = new CRS_Matrix<Double>();
  StdVector<Double>A_H_dP_Re;
  StdVector<Double>A_H_dP_Im;
  TMP_R->MultTriple_MKL(*R, *Restriction_, A_H_rP, A_H_cP, A_H_dP_Re, 2, false);
  TMP_R->MultTriple_MKL(*I, *Restriction_, A_H_rP, A_H_cP, A_H_dP_Im, 2, false);

  A_H_dP.Resize(A_H_dP_Re.GetSize());
  for (UInt i = 0; i < A_H_dP_Re.GetSize(); ++i) {
      A_H_dP[i].real(A_H_dP_Re[i]);
      A_H_dP[i].imag(A_H_dP_Im[i]);
  }

  delete [] re;
  delete [] im;
  delete R; R = NULL;
  delete I; I = NULL;
  delete TMP_R; TMP_R = NULL;

}

#undef CALL_AND_CHECK_STATUS

} // namespace CoupledField
