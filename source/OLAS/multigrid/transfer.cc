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
                        const CRS_Matrix<Double>& coarseAuxMat,
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

  // just a performance improvement...outsorced from the O(~N^2) loop
  StdVector<Integer> eIN (edgeIndNode.GetSize() * 2);
  for(UInt i = 0; i < edgeIndNode.GetSize(); ++i){
	  const StdVector<Integer>& e = edgeIndNode[i];
	  eIN[2 * i] = e[0];
	  eIN[2 * i + 1] = e[1];
  }

  UInt count = 0;
  StdVector<Integer> eInd (2);
  UInt aggN1, aggN2;
  //TODO the following O(~N^2) is way too expensive
  for(UInt i = 0; i < sizeH; ++i){
    // get edge i, nodes of these edges are coarse-nodes
    const StdVector<Integer>& cEind = A.GetCoarseEdgeInd(i);
    const StdVector<Integer>& cEreal = A.GetCoarseEdgeReal(i);
    StdVector<Integer>& cE = cEdgeIndNode_[i];
    cE.Resize(2,0);
    for(UInt j = 0; j < sizeh; ++j){
      //const StdVector<Integer>& e = edgeIndNode[j];
      //eInd[1] = A.GetIndexOfRealNode(e[1]);
      //eInd[0] = A.GetIndexOfRealNode(e[0]);
      eInd[1] = A.GetIndexOfRealNode(eIN[2 * j]);
      eInd[0] = A.GetIndexOfRealNode(eIN[2 * j + 1]);
      if( (eInd[0] == -1) || (eInd[1] == -1) ){
        EXCEPTION("coarse index problem");
      }
      aggN1 = A.GetAgglomerateNumOfNode(eInd[0]);
      aggN2 = A.GetAgglomerateNumOfNode(eInd[1]);
      if( ((Integer)aggN1 == cEind[0]) && ((Integer)aggN2 == cEind[1]) ){
        pColR.Push_back(j);
        pDatR.Push_back(1.0);
        cE[0] = cEreal[0];
        cE[1] = cEreal[1];
        count++;
      }
      if( ((Integer)aggN1 == cEind[1]) && ((Integer)aggN2 == cEind[0]) ){
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



template <typename T>
void TransferOperator<T>::Prolongate( const Vector<T>& v_H,
                                            Vector<T>& v_h,
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

/**********************************************************/

template <typename T>
void TransferOperator<T>::Restrict( const Vector<T>& v_h,
                                          Vector<T>& v_H,
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
ExportCSRMatrix(sparse_matrix_t AH,
				  int &rows,
				  int &cols,
				  int *&pointerB_1,
				  int *&pointerE_1,
				  int *&solCols,
				  double *&values){

	sparse_index_base_t solbase;
	CALL_AND_CHECK_STATUS(mkl_sparse_d_export_csr( AH, &solbase, &rows, &cols,
			&pointerB_1, &pointerE_1, &solCols, &values),
			"Error after MKL_SPARSE_D_EXPORT_CSR, Export of AH \n");
}

template <typename T>
void TransferOperator<T>::
ExportCSRMatrix(sparse_matrix_t AH,
				  int &rows,
				  int &cols,
				  int *&pointerB_1,
				  int *&pointerE_1,
				  int *&solCols,
				  std::complex<double> *&values){

	EXCEPTION("transfer.cc: Complex Galerkin product not yet handled");

	//sparse_index_base_t solbase;
	//CALL_AND_CHECK_STATUS(mkl_sparse_c_export_csr( AH, &solbase, &rows, &cols,
	//		&pointerB_1, &pointerE_1, &solCols, &values),
	//		"Error after MKL_SPARSE_D_EXPORT_CSR, Export of AH \n");
}


template <typename T>
void TransferOperator<T>::CreateCSRMatrix(sparse_matrix_t& Ah,
		  const sparse_index_base_t& t,
		  const int &rowsAh,
		  const int &colsAh,
		  int *rPAh,
		  int *cPAh,
		  double *dPAh){


	  CALL_AND_CHECK_STATUS(mkl_sparse_d_create_csr( &Ah, t, rowsAh, colsAh,
	      rPAh, rPAh+1, cPAh, dPAh ),
	      "Error after MKL_SPARSE_D_CREATE_CSR, Ah-matrix \n");
}



template <typename T>
void TransferOperator<T>::CreateCSRMatrix(sparse_matrix_t &Ah,
		  const sparse_index_base_t& t,
		  const int &rowsAh,
		  const int &colsAh,
		  int *rPAh,
		  int *cPAh,
		  std::complex<double> *dPAh){

	EXCEPTION("transfer.cc: Complex Galerkin product not yet handled");


	  //CALL_AND_CHECK_STATUS(mkl_sparse_c_create_csr( &Ah, t, rowsAh, colsAh,
	  //    rPAh, rPAh+1, cPAh, dPAh ),
	  //    "Error after MKL_SPARSE_D_CREATE_CSR, Ah-matrix \n");
}

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

  /***************************************************************
   ******************** 1) Coarse System matrix ******************
   ***************************************************************/
  sparse_matrix_t P = NULL;
  sparse_index_base_t t = SPARSE_INDEX_BASE_ZERO;
  MKL_INT rowsProlo = Prolongation_->GetNumRows();
  MKL_INT colsProlo = Prolongation_->GetNumCols();
  MKL_INT *rPProlo = (Integer*)Prolongation_->GetRowPointer();
  MKL_INT *cPProlo = (Integer*)Prolongation_->GetColPointer();
  double *dPProlo = Prolongation_->GetDataPointer();

  sparse_matrix_t Ah = NULL;
  MKL_INT rowsAh = A_h.GetNumRows();
  MKL_INT colsAh = A_h.GetNumCols();
  MKL_INT *rPAh = (Integer*) A_h.GetRowPointer();
  MKL_INT *cPAh = (Integer*) A_h.GetColPointer();
  T *dPAh = (T*) A_h.GetDataPointer();


  /******** Convert our CRS_Matrices into MKL-sparse csr-matrix *********/
  CALL_AND_CHECK_STATUS(mkl_sparse_d_create_csr( &P, t, rowsProlo, colsProlo,
      rPProlo, rPProlo+1, cPProlo, dPProlo ),
      "Error after MKL_SPARSE_D_CREATE_CSR, Prolongation-matrix \n");

  // specialization between complex and double handled in own function
  CreateCSRMatrix(Ah, t, rowsAh, colsAh, rPAh, cPAh, dPAh );


  /*************** Perform AH = P * Ah * P^T ***********************/
  sparse_matrix_t tempAH = NULL;
  sparse_operation_t transA = SPARSE_OPERATION_NON_TRANSPOSE;
  CALL_AND_CHECK_STATUS(mkl_sparse_spmm(transA, Ah, P, &tempAH),
      "Error after MKL_SPARSE_SPMM, tempAH = Ah * P \n");

  sparse_matrix_t AH = NULL;
  sparse_operation_t transAH = SPARSE_OPERATION_TRANSPOSE;
  CALL_AND_CHECK_STATUS(mkl_sparse_spmm(transAH, tempAH, P, &AH),
      "Error after MKL_SPARSE_SPMM, AH = P^T * tempAH \n");


  /*************** Export AH in MKL-csr format ***********************/

  MKL_INT rows = -1,
		  cols = -1;
  MKL_INT *solCols = NULL;
  T *values = NULL;


  MKL_INT *pointerB_1 = NULL, *pointerE_1 = NULL;

  // specialization between complex and double handled in own function
  ExportCSRMatrix(AH, rows, cols, pointerB_1, pointerE_1, solCols, values);


  A_H_rP.Resize(rows + 1, 0);
  // I'm sure there's a method in MKL to find nnz, but I couldn't find it
  UInt nnz;
  nnz = 0;
  for(Integer i = 0; i < rows; i++ ){
    nnz += pointerE_1[i] - pointerB_1[i];
    A_H_rP[i + 1] = nnz;
  }


  A_H_cP.Resize(nnz, 0);
  A_H_dP.Resize(nnz, 0.0);
  for(UInt i = 0; i < (UInt)rows; ++i){
    for(UInt j = (UInt)A_H_rP[i]; j < (UInt)A_H_rP[i+1]; ++j){
      A_H_cP[j] = (UInt)solCols[j];
      A_H_dP[j] = values[j];
    }
  }


  /*********** Release handles and deallocate memory ************/
  if( mkl_sparse_destroy( P ) != SPARSE_STATUS_SUCCESS)
  { printf(" Error after MKL_SPARSE_DESTROY, P \n");fflush(0); }

  if( mkl_sparse_destroy( Ah ) != SPARSE_STATUS_SUCCESS)
  { printf(" Error after MKL_SPARSE_DESTROY, Ah \n");fflush(0); }

  if( mkl_sparse_destroy( tempAH ) != SPARSE_STATUS_SUCCESS)
  { printf(" Error after MKL_SPARSE_DESTROY, tempAH \n");fflush(0); }

  if( mkl_sparse_destroy( AH ) != SPARSE_STATUS_SUCCESS)
  { EXCEPTION("Error in the MKL sparse matrix-matrix multiplication\n"
      "try to increase the size of the coarse system");}




  //if(amgType != SCALAR){
  /*****************************************************************/
  /********************* 2) Coarse Auxiliary matrix ****************/
  /*****************************************************************/
  sparse_matrix_t PB = NULL;
  sparse_index_base_t tB = SPARSE_INDEX_BASE_ZERO;
  MKL_INT rowsProloB = AuxProlongation_->GetNumRows();
  MKL_INT colsProloB = AuxProlongation_->GetNumCols();
  MKL_INT *rPProloB = (Integer*)AuxProlongation_->GetRowPointer();
  MKL_INT *cPProloB = (Integer*)AuxProlongation_->GetColPointer();
  double *dPProloB = AuxProlongation_->GetDataPointer();

  sparse_matrix_t Bh = NULL;
  MKL_INT rowsBh = B_h.GetNumRows();
  MKL_INT colsBh = B_h.GetNumCols();
  MKL_INT *rPBh = (Integer*) B_h.GetRowPointer();
  MKL_INT *cPBh = (Integer*) B_h.GetColPointer();
  double *dPBh = (double*) B_h.GetDataPointer();


  /********* Convert our CRS_Matrices into MKL-sparse csr-matrix **********/
  CALL_AND_CHECK_STATUS(mkl_sparse_d_create_csr( &PB, tB, rowsProloB, colsProloB,
      rPProloB, rPProloB+1, cPProloB, dPProloB ),
      "Error after MKL_SPARSE_D_CREATE_CSR, Aux-Prolongation-matrix \n");

  CALL_AND_CHECK_STATUS(mkl_sparse_d_create_csr( &Bh, tB, rowsBh, colsBh,
      rPBh, rPBh+1, cPBh, dPBh ),
      "Error after MKL_SPARSE_D_CREATE_CSR, Bh-matrix \n");



  /**************** Perform BH = P * Bh * P^T ************************/
  sparse_matrix_t tempBH = NULL;
  sparse_operation_t transB = SPARSE_OPERATION_NON_TRANSPOSE;
  CALL_AND_CHECK_STATUS(mkl_sparse_spmm(transB, Bh, PB, &tempBH),
      "Error after MKL_SPARSE_SPMM, tempBH = Bh * P \n");

  sparse_matrix_t BH = NULL;
  sparse_operation_t transBH = SPARSE_OPERATION_TRANSPOSE;
  CALL_AND_CHECK_STATUS(mkl_sparse_spmm(transBH, tempBH, PB, &BH),
      "Error after MKL_SPARSE_SPMM, BH = P^T * tempBH \n");


  /**************** Export BH in MKL-csr format ************************/
  sparse_index_base_t solbaseB;
  MKL_INT rowsB, colsB;
  MKL_INT *solColsB;
  double *valuesB = NULL;
  MKL_INT *pointerB_1B, *pointerE_1B;
  CALL_AND_CHECK_STATUS(mkl_sparse_d_export_csr( BH, &solbaseB, &rowsB, &colsB,
      &pointerB_1B, &pointerE_1B, &solColsB, &valuesB),
      "Error after MKL_SPARSE_D_EXPORT_CSR, Export of BH \n");



  B_H_rP.Resize(rowsB + 1, 0);
  // I'm sure there's a method in MKL to find nnz, but I couldn't find it
  nnz = 0;
  for(Integer i = 0; i < rowsB; i++ ){
    nnz += pointerE_1B[i] - pointerB_1B[i];
    B_H_rP[i + 1] = nnz;
  }


  B_H_cP.Resize(nnz, 0);
  B_H_dP.Resize(nnz, 0.0);
  for(UInt i = 0; i < (UInt)rowsB; ++i){
    for(UInt j = B_H_rP[i]; j < B_H_rP[i+1]; ++j){
      B_H_cP[j] = (UInt)solColsB[j];
      B_H_dP[j] = valuesB[j];
    }
  }

  /************ Release handles and deallocate memory *************/
  if( mkl_sparse_destroy( PB ) != SPARSE_STATUS_SUCCESS)
  { printf(" Error after MKL_SPARSE_DESTROY, P \n");fflush(0); }

  if( mkl_sparse_destroy( Bh ) != SPARSE_STATUS_SUCCESS)
  { printf(" Error after MKL_SPARSE_DESTROY, Bh \n");fflush(0); }

  if( mkl_sparse_destroy( tempBH ) != SPARSE_STATUS_SUCCESS)
  { printf(" Error after MKL_SPARSE_DESTROY, tempBH \n");fflush(0); }

  if( mkl_sparse_destroy( BH ) != SPARSE_STATUS_SUCCESS)
  { EXCEPTION("Error in the MKL sparse matrix-matrix multiplication\n"
      "try to increase the size of the coarse system");
    /*printf(" Error after MKL_SPARSE_DESTROY, BH \n");fflush(0); */}

  //}

}


/***************************************************************************/
template <typename T>
void TransferOperator<T>::
GalerkinProductEdgeAux(StdVector<UInt>& B_H_rP,
                      StdVector<UInt>& B_H_cP,
                      StdVector<Double>& B_H_dP,
                      const CRS_Matrix<Double>& B_h)
{
  //if(amgType != SCALAR){
  /*****************************************************************/
  /*********************  Coarse Auxiliary matrix ******************/
  /*****************************************************************/
  sparse_matrix_t PB = NULL;
  sparse_index_base_t tB = SPARSE_INDEX_BASE_ZERO;
  MKL_INT rowsProloB = AuxRestriction_->GetNumRows();
  MKL_INT colsProloB = AuxRestriction_->GetNumCols();
  MKL_INT *rPProloB = (Integer*)AuxRestriction_->GetRowPointer();
  MKL_INT *cPProloB = (Integer*)AuxRestriction_->GetColPointer();
  double *dPProloB = AuxRestriction_->GetDataPointer();

  sparse_matrix_t Bh = NULL;
  MKL_INT rowsBh = B_h.GetNumRows();
  MKL_INT colsBh = B_h.GetNumCols();
  MKL_INT *rPBh = (Integer*) B_h.GetRowPointer();
  MKL_INT *cPBh = (Integer*) B_h.GetColPointer();
  double *dPBh = (double*) B_h.GetDataPointer();

  /********* Convert our CRS_Matrices into MKL-sparse csr-matrix **********/
  CALL_AND_CHECK_STATUS(mkl_sparse_d_create_csr( &PB, tB, rowsProloB, colsProloB,
      rPProloB, rPProloB+1, cPProloB, dPProloB ),
      "Error after MKL_SPARSE_D_CREATE_CSR, Aux-Prolongation-matrix \n");

  CALL_AND_CHECK_STATUS(mkl_sparse_d_create_csr( &Bh, tB, rowsBh, colsBh,
      rPBh, rPBh+1, cPBh, dPBh ),
      "Error after MKL_SPARSE_D_CREATE_CSR, Bh-matrix \n");



  /**************** Perform BH = P * Bh * P^T ************************/
  sparse_matrix_t tempBH = NULL;
  sparse_operation_t transB = SPARSE_OPERATION_NON_TRANSPOSE;
  CALL_AND_CHECK_STATUS(mkl_sparse_spmm(transB, PB, Bh, &tempBH),
      "Error after MKL_SPARSE_SPMM, tempBH = R * Bh \n");

  sparse_matrix_t tempBHT = NULL;
  sparse_operation_t t = SPARSE_OPERATION_TRANSPOSE;
  CALL_AND_CHECK_STATUS(mkl_sparse_convert_csr(tempBH, t, &tempBHT),
        "Error after MKL_SPARSE_CONVERT_CSR, tempBHT = tempBH^T \n");


  sparse_matrix_t BH = NULL;
  sparse_operation_t transBH = SPARSE_OPERATION_NON_TRANSPOSE;
  CALL_AND_CHECK_STATUS(mkl_sparse_spmm(transBH, PB, tempBHT, &BH),
      "Error after MKL_SPARSE_SPMM, BH = tempBH * R^T \n");


  /**************** Export BH in MKL-csr format ************************/
  sparse_index_base_t solbaseB;
  MKL_INT rowsB, colsB;
  MKL_INT *solColsB;
  double *valuesB = NULL;
  MKL_INT *pointerB_1B, *pointerE_1B;
  CALL_AND_CHECK_STATUS(mkl_sparse_d_export_csr( BH, &solbaseB, &rowsB, &colsB,
      &pointerB_1B, &pointerE_1B, &solColsB, &valuesB),
      "Error after MKL_SPARSE_D_EXPORT_CSR, Export of BH \n");


  B_H_rP.Resize(rowsB + 1, 0);
  // I'm sure there's a method in MKL to find nnz, but I couldn't find it
  UInt nnz = 0;
  for(Integer i = 0; i < rowsB; i++ ){
    nnz += pointerE_1B[i] - pointerB_1B[i];
    B_H_rP[i + 1] = nnz;
  }

  B_H_cP.Resize(nnz, 0);
  B_H_dP.Resize(nnz, 0.0);
  for(UInt i = 0; i < (UInt)rowsB; ++i){
    for(UInt j = B_H_rP[i]; j < B_H_rP[i+1]; ++j){
      B_H_cP[j] = (UInt)solColsB[j];
      B_H_dP[j] = valuesB[j];
    }
  }


  /************ Release handles and deallocate memory *************/
  if( mkl_sparse_destroy( PB ) != SPARSE_STATUS_SUCCESS)
  { printf(" Error after MKL_SPARSE_DESTROY, P \n");fflush(0); }

  if( mkl_sparse_destroy( Bh ) != SPARSE_STATUS_SUCCESS)
  { printf(" Error after MKL_SPARSE_DESTROY, Bh \n");fflush(0); }

  if( mkl_sparse_destroy( tempBH ) != SPARSE_STATUS_SUCCESS)
  { printf(" Error after MKL_SPARSE_DESTROY, tempBH \n");fflush(0); }

  if( mkl_sparse_destroy( tempBHT ) != SPARSE_STATUS_SUCCESS)
  { printf(" Error after MKL_SPARSE_DESTROY, tempBHT \n");fflush(0); }

  if( mkl_sparse_destroy( BH ) != SPARSE_STATUS_SUCCESS)
  { EXCEPTION("Error in the MKL sparse matrix-matrix multiplication\n"
      "try to increase the size of the coarse system");
    /*printf(" Error after MKL_SPARSE_DESTROY, BH \n");fflush(0); */}

  //}

}


/***************************************************************/
template <typename T>
void TransferOperator<T>::
GalerkinProductEdgeSys(StdVector<UInt>& A_H_rP,
                        StdVector<UInt>& A_H_cP,
                        StdVector<T>& A_H_dP,
                        const CRS_Matrix<T>& A_h)
{
  //if(amgType != SCALAR){
  /*****************************************************************/
  /*********************  Coarse System matrix *********************/
  /*****************************************************************/
  sparse_matrix_t PB = NULL;
  sparse_index_base_t tB = SPARSE_INDEX_BASE_ZERO;
  MKL_INT rowsProloB = Restriction_->GetNumRows();
  MKL_INT colsProloB = Restriction_->GetNumCols();
  MKL_INT *rPProloB = (Integer*)Restriction_->GetRowPointer();
  MKL_INT *cPProloB = (Integer*)Restriction_->GetColPointer();
  double *dPProloB = Restriction_->GetDataPointer();

  sparse_matrix_t Ah = NULL;
  MKL_INT rowsBh = A_h.GetNumRows();
  MKL_INT colsBh = A_h.GetNumCols();
  MKL_INT *rPBh = (Integer*) A_h.GetRowPointer();
  MKL_INT *cPBh = (Integer*) A_h.GetColPointer();
  double *dPBh = (double*) A_h.GetDataPointer();


  /********* Convert our CRS_Matrices into MKL-sparse csr-matrix **********/
  CALL_AND_CHECK_STATUS(mkl_sparse_d_create_csr( &PB, tB, rowsProloB, colsProloB,
      rPProloB, rPProloB+1, cPProloB, dPProloB ),
      "Error after MKL_SPARSE_D_CREATE_CSR, Aux-Prolongation-matrix \n");

  CALL_AND_CHECK_STATUS(mkl_sparse_d_create_csr( &Ah, tB, rowsBh, colsBh,
      rPBh, rPBh+1, cPBh, dPBh ),
      "Error after MKL_SPARSE_D_CREATE_CSR, Ah-matrix \n");



  /**************** Perform AH = P * Ah * P^T ************************/
  sparse_matrix_t tempBH = NULL;
  sparse_operation_t transB = SPARSE_OPERATION_NON_TRANSPOSE;
  CALL_AND_CHECK_STATUS(mkl_sparse_spmm(transB, PB, Ah, &tempBH),
      "Error after MKL_SPARSE_SPMM, tempBH = R * Bh \n");

  sparse_matrix_t tempBHT = NULL;
  sparse_operation_t t = SPARSE_OPERATION_TRANSPOSE;
  CALL_AND_CHECK_STATUS(mkl_sparse_convert_csr(tempBH, t, &tempBHT),
        "Error after MKL_SPARSE_CONVERT_CSR, tempBHT = tempBH^T \n");


  sparse_matrix_t AH = NULL;
  sparse_operation_t transBH = SPARSE_OPERATION_NON_TRANSPOSE;
  CALL_AND_CHECK_STATUS(mkl_sparse_spmm(transBH, PB, tempBHT, &AH),
      "Error after MKL_SPARSE_SPMM, BH = tempBH * R^T \n");


  /**************** Export AH in MKL-csr format ************************/
  //sparse_index_base_t solbaseB;
  MKL_INT rowsB, colsB;
  MKL_INT *solColsB;
  double *valuesB = NULL;
  MKL_INT *pointerB_1B, *pointerE_1B;

  ExportCSRMatrix(AH, rowsB, colsB, pointerB_1B, pointerE_1B, solColsB, valuesB);

  //CALL_AND_CHECK_STATUS(mkl_sparse_d_export_csr( AH, &solbaseB, &rowsB, &colsB,
  //    &pointerB_1B, &pointerE_1B, &solColsB, &valuesB),
  //    "Error after MKL_SPARSE_D_EXPORT_CSR, Export of AH \n");


  A_H_rP.Resize(rowsB + 1, 0);
  // I'm sure there's a method in MKL to find nnz, but I couldn't find it
  UInt nnz = 0;
  for(Integer i = 0; i < rowsB; i++ ){
    nnz += pointerE_1B[i] - pointerB_1B[i];
    A_H_rP[i + 1] = nnz;
  }

  A_H_cP.Resize(nnz, 0);
  A_H_dP.Resize(nnz, 0.0);
  for(UInt i = 0; i < (UInt)rowsB; ++i){
    for(UInt j = A_H_rP[i]; j < A_H_rP[i+1]; ++j){
      A_H_cP[j] = (UInt)solColsB[j];
      A_H_dP[j] = valuesB[j];
    }
  }


  /************ Release handles and deallocate memory *************/
  if( mkl_sparse_destroy( PB ) != SPARSE_STATUS_SUCCESS)
  { printf(" Error after MKL_SPARSE_DESTROY, P \n");fflush(0); }

  if( mkl_sparse_destroy( Ah ) != SPARSE_STATUS_SUCCESS)
  { printf(" Error after MKL_SPARSE_DESTROY, Bh \n");fflush(0); }

  if( mkl_sparse_destroy( tempBH ) != SPARSE_STATUS_SUCCESS)
  { printf(" Error after MKL_SPARSE_DESTROY, tempBH \n");fflush(0); }

  if( mkl_sparse_destroy( tempBHT ) != SPARSE_STATUS_SUCCESS)
  { printf(" Error after MKL_SPARSE_DESTROY, tempBHT \n");fflush(0); }

  if( mkl_sparse_destroy( AH ) != SPARSE_STATUS_SUCCESS)
  { EXCEPTION("Error in the MKL sparse matrix-matrix multiplication\n"
      "try to increase the size of the coarse system");}

}
#undef CALL_AND_CHECK_STATUS

} // namespace CoupledField
