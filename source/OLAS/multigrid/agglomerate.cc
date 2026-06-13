#include "OLAS/multigrid/agglomerate.hh"


namespace CoupledField {


template <typename T>
Agglomerate<T>::Agglomerate()
    : NhStartIndex_(NULL),
      NhEdges_(NULL),
      Size_h_(0),
      Size_H_(0),
      maxNumA_(20),
      minNumA_(5)
      {}

/**********************************************************/

template <typename T>
Agglomerate<T>::Agglomerate( const CRS_Matrix<Double>& auxMat)
    : NhStartIndex_(NULL),
      NhEdges_(NULL),
      Size_h_(0),
      Size_H_(0),
      maxNumA_(20),
      minNumA_(5)
{

    CreateAgglomerateGraphs( auxMat);
}

/**********************************************************/

template <typename T>
Agglomerate<T>::~Agglomerate()
{
    
    Reset();
}

/**********************************************************/

template <typename T>
inline void Agglomerate<T>::InsertEdgeList(const StdVector< StdVector< Integer> >& edgeIndNode)
{
  edgeIndNode_ = edgeIndNode;
}

/**********************************************************/

template <typename T>
inline void Agglomerate<T>::InsertNodeIndex(const StdVector< Integer>& nodeNumIndex)
{
  nodeNumIndex_ = nodeNumIndex;

  this->FillIndexNodeNum();
}

/**********************************************************/

template <typename T>
inline const StdVector< StdVector< Integer> > Agglomerate<T>::GetAgglomerates() const
{
    return Agglomerates_;
}

/**********************************************************/

template <typename T>
inline UInt Agglomerate<T>::GetSizeh() const
{
    return Size_h_;
}

/**********************************************************/

template <typename T>
inline UInt Agglomerate<T>::GetSizeH() const
{
    return Size_H_;
}

/**********************************************************/

template <typename T>
inline const StdVector<Integer> Agglomerate<T>::GetCoarseFineSplitting() const
{
    return CoarseIndex_;
}

/**********************************************************/

template <typename T>
inline bool Agglomerate<T>::IsFPoint( const int i ) const
{
    return CoarseIndex_[i] == FINE;
}

/**********************************************************/

template <typename T>
inline bool Agglomerate<T>::IsCPoint( const int i ) const
{
    return CoarseIndex_[i] >= COARSE;
}

/**********************************************************/

template <typename T>
inline Integer Agglomerate<T>::GetCoarseIndex( const Integer i ) const
{
    return CoarseIndex_[i];
}


/**********************************************************/

template <typename T>
inline Integer Agglomerate<T>::GetIndexOfCoarseNode( const Integer i ) const
{
  return IndexOfCoarse_[i];
}

/**********************************************************/

template <typename T>
inline StdVector<Integer> Agglomerate<T>::GetAgglomerateOfNode( const Integer i ) const
{
  // which agglomerate belongs to node i?
  UInt aggNum = nodeAgg_[i];
  return Agglomerates_[(Integer)aggNum];
}

/**********************************************************/

template <typename T>
inline UInt Agglomerate<T>::GetAgglomerateNumOfNode( const Integer i ) const
{
  return nodeAgg_[i];
}
/**********************************************************/

template <typename T>
inline Integer Agglomerate<T>::GetRealNodeNumOfCIndex( const Integer i ) const
{
  if(i != -1){
    return coarseRealNN_[i];
  }else{
    EXCEPTION("GetRealNodeNumOfCIndex: agglomerate error");
    return -1;
  }
}
/**********************************************************/

template <typename T>
inline StdVector<Integer> Agglomerate<T>::GetCNodeNumIndex() const
{
  return coarseRealNN_;

}

/**********************************************************/

template <typename T>
inline UInt Agglomerate<T>::GetNumCoarseEdges() const
{
    return numCEdges_;
}

/**********************************************************/

template <typename T>
inline UInt Agglomerate<T>::GetNumFineEdges() const
{
    return edgeIndNode_.GetSize();
}

/**********************************************************/

template <typename T>
inline StdVector<Integer> Agglomerate<T>::GetCoarseEdgeReal( const int i ) const
{
  return cEdgesReal_[i];
}

/**********************************************************/

template <typename T>
inline StdVector<Integer> Agglomerate<T>::GetCoarseEdgeInd( const int i ) const
{
  return cEdgesInd_[i];
}

/**********************************************************/

template <typename T>
inline Integer Agglomerate<T>::GetIndexOfRealNode( UInt i ) const
{

  Integer t = -1;
  t = indexNodeNum_.at(i);
  return t;
}

/**********************************************************/

template <typename T>
Integer Agglomerate<T>::
GetNumInterpolatedPoints() const
{
    Integer totalNumPoints = 0;
    // run over all coarse grid points
    for( Integer i = 0; i < (Integer)GetSizeH(); i++ ) {
      const StdVector<Integer>& agg = Agglomerates_[i];
      if(agg.GetSize() <= 0)EXCEPTION("Empty agglomerate, decrease minNumNodes");
      totalNumPoints += agg.GetSize();
    }
    return totalNumPoints;
}


/**********************************************************/

template <typename T>
Integer Agglomerate<T>::CreateAgglomerateGraphs( const CRS_Matrix<Double>& auxMat)
{
    // loop indices, consistent for the whole function
    UInt i,  // i = row index in auxMat, in [1,n]
         ij; // ij = absolute index in pCol and pDat for an entry in row i

    // pointers for direct matrix access
    const UInt *const pRow = auxMat.GetRowPointer();
    const UInt *const pCol = auxMat.GetColPointer();

    if(auxMat.GetCurrentLayout() != auxMat.LEX_DIAG_FIRST){
      EXCEPTION("AMG: Agglomerate needs Auxiliary auxMat to be LEX_DIAG_FIRST");
    }

    Size_h_ = auxMat.GetNumRows();
    
    /*********************************************************************/
    /********* Step 1) Define neighbourhoods for every node **************/
    /*********************************************************************/
    StdVector< StdVector< UInt> > Neigh;
    Neigh.Resize(Size_h_);
    // loop over every node and store it's neighbourhood
    for(i = 0; i < Size_h_; ++i){
      StdVector<UInt>& a = Neigh[i];
      // loop over all direct neighbours
      for( ij = pRow[i]+1; ij < pRow[i+1]; ij++ ) {
        a.Push_back(pCol[ij]);
      }

      // also store the second layer (neighbours of direct neighbours)
      UInt numN = 0;
      UInt s = a.GetSize();
      for(UInt nI = 0; nI < s; ++nI){
        // neighbour-point
        UInt n = a[nI];
        // loop over all direct neighbours of this neighbour-point
        for( ij = pRow[n]+1; ij < pRow[n+1]; ij++ ) {
          // don't add the node twice
          if(!a.Contains(pCol[ij])){
            a.Push_back(pCol[ij]);
            numN++;
          }
        }
      }

      // add a final third layer
      s = a.GetSize();
      for(UInt nI = numN; nI < s; ++nI){
        // neighbour-point
        UInt n = a[nI];
        // loop over all direct neighbours of this neighbour-point
        for( ij = pRow[n]+1; ij < pRow[n+1]; ij++ ) {
          // don't add the node twice
          if(!a.Contains(pCol[ij]))a.Push_back(pCol[ij]);
        }
      }
    }

    /*********************************************************************/
    /********* Step 2) Define the agglomerates ***************************/
    /*********************************************************************/
    // we neeed a temporary construct, because some entries will
    // be deleted and we would get some indexing troubles
    StdVector< StdVector< Integer> > tempAgglomerates;

    // list of possible root points (COARSE) and their connected points in the
    // agglomerates (FINE)
    // all nodes are initialized as undefined in the first step
    CoarseIndex_.Resize(Size_h_, (Integer)UNDEFINED);

    // stores agglomerate-number for every node
    nodeAgg_.Resize(Size_h_, 0);

    // number of undefined points
    UInt undefN = CoarseIndex_.GetSize();

    StdVector<UInt> smallAgg;
    UInt iter = 0;
    StdVector<Integer>::iterator test;
    while( undefN > 0){
      // pick an undefined point, no matter which one
      test = std::find(CoarseIndex_.Begin(), CoarseIndex_.End(), (Integer)UNDEFINED);
      UInt actInd = test - CoarseIndex_.Begin();
      //UInt actInd = CoarseIndex_.Find((Integer)UNDEFINED);
      // this is a root-point, mark it as COARSE
      CoarseIndex_[actInd] = (Integer)COARSE;
      nodeAgg_[actInd] = iter;
      undefN--;
      StdVector<Integer> agg;
      // this root point will be the first entry of every agglomerate
      agg.Push_back(actInd);
      // get the neighbours of node actInd.
      // it can actually not happen that there is only a single point in n
      const StdVector<UInt>& n = Neigh[actInd];

      // only add undefined nodes to aggregate
      for( i = 0; i < n.GetSize(); ++i){
        if(CoarseIndex_[n[i]] == (Integer)UNDEFINED){
          agg.Push_back(n[i]);
          CoarseIndex_[n[i]] = (Integer)FINE;
          nodeAgg_[n[i]] = iter;
          undefN--;
        }
      }
      tempAgglomerates.Push_back(agg);

      // mark the agglomerate if its size is smaller than minNumA_
      if(agg.GetSize() < minNumA_) smallAgg.Push_back(iter);
      iter++;
    }

    // Now it probably happened that some agglomerates are smaller than minNumA_
    // therefore we merge the nodes of the small ones to the larger ones, until
    // the resulting agglomerate will not be larger than maxNumA_, if that happens,
    // merge them with another agglomerate

    // vector has 1 entry if agglomerate gets deleted, else it's 0
    StdVector<UInt> deleteAgg;
    deleteAgg.Resize(tempAgglomerates.GetSize(), 0);
//TODO Clean this up
    for(UInt i = 0; i < smallAgg.GetSize(); ++i){
      // get the whole small agglomerate
      StdVector<Integer> sAgg = tempAgglomerates[smallAgg[i]];
      // don't forget the root-point, it also belongs to the agglomerate
      // and it becomes a fine-node
      CoarseIndex_[sAgg[0]] = (Integer)FINE;

      // loop over each node of the small agglomerate and add it to
      // one of the agglomerates of its neighbours
      for(UInt j = 0; j < sAgg.GetSize(); ++j){
        // get neighbours of node smallAgg[j]
        const StdVector<UInt>& n = Neigh[sAgg[j]];
        // it can happen that the following if-clause isn't met for any node in
        // n, therefore we add the smallAgg entry to the smalles of the >maxNumA_
        UInt backupAggNumNode = nodeAgg_[n[0]];
        bool added = false;
        for(UInt neighIt = 0; neighIt < n.GetSize(); ++neighIt){
          // get neighIt's agglomerate-number
          UInt aggNumNode = nodeAgg_[n[neighIt]];
          StdVector<Integer>& tAgg = tempAgglomerates[aggNumNode];
          // add the node from the small-agglomerate to the neighbour-ones
          // only if the merged-size isn't larger than maxNumA_ and smaller
          // than minNumA_
          if( (tAgg.GetSize() < maxNumA_) &&
              (tAgg.GetSize() > minNumA_) &&
              !smallAgg.Contains(aggNumNode) ){
            tAgg.Push_back(sAgg[j]);
            nodeAgg_[sAgg[j]] = aggNumNode;
            added = true;
            break;
          }
          Integer bA = -1;
          Integer backupA = -1;
          StdVector<Integer>& backupAgg = tempAgglomerates[backupAggNumNode];
          if(tAgg.GetSize() > maxNumA_) bA = tAgg.GetSize() -  maxNumA_;
          if(backupAgg.GetSize() > maxNumA_){
            backupA = backupAgg.GetSize() - maxNumA_;
          }
          if( (bA < backupA) && !smallAgg.Contains(aggNumNode)) backupAggNumNode = aggNumNode;
        }
        if(added == false){
          tempAgglomerates[backupAggNumNode].Push_back(sAgg[j]);
          nodeAgg_[sAgg[j]] = backupAggNumNode;
        }

      }
      // store the agglomerate to be deleted
      deleteAgg[smallAgg[i]] = 1;
    }
    // define the coarse-indices for every coarse-point
    Integer cInd = 0;
    for(UInt i = 0; i < CoarseIndex_.GetSize(); ++i){
      if( CoarseIndex_[i] >= (Integer)COARSE ){
        CoarseIndex_[i] = cInd;
        IndexOfCoarse_.Push_back(i);
        cInd++;
      }
    }

    // fill the final Agglomerate_ graph
    // Index of every agglomerate will be the coarse node index of
    // the first node of the agglomerate
    Agglomerates_.Resize(cInd);
    for(UInt i = 0; i < tempAgglomerates.GetSize(); ++i){
      if(deleteAgg[i] != 1){
        const StdVector<Integer>& t = tempAgglomerates[i];
        Agglomerates_[CoarseIndex_[t[0]]] = t;
      }
    }

    nodeAgg_.Resize(Size_h_, 0);
    for(UInt i = 0; i < Size_h_; ++i){
      for(UInt aggIt = 0; aggIt < Agglomerates_.GetSize(); ++aggIt){
        StdVector<Integer>& a = Agglomerates_[aggIt];
              if( a.Contains(i) ) nodeAgg_[i] = CoarseIndex_[a[0]];
      }
    }

    Size_H_ = Agglomerates_.GetSize();

    this->CreateCoarseRealNodeNum();

    return Size_H_;
}

/**********************************************************/



template <typename T>
UInt Agglomerate<T>::CreateCoarseEdges(const CRS_Matrix<Double>& coarseAuxMat){
  /*
   * Idea: every edge in the auxiliary-matrix represents an edge between
   * two nodes. But these two nodes can lie in the same aggregate, therefore
   * they don't form a coarse-grid edge. If however the nodes belong to different
   * aggregates, there is a coarse grid edge...an entry in the coarse
   *  system matrix.
   *
   *  The plan is to loop over every edge in the coarse auxiliary-matrix
   *  and check if the nodes belong to different aggregates. If this is
   *  the case, insert a coarse edge
   */


  if(coarseAuxMat.GetCurrentLayout() != coarseAuxMat.LEX_DIAG_FIRST){
    EXCEPTION("coarse auxiliary matrix is not stored as LEX_DIAG_FIRST...Implementation error");
  }

  if( (cEdgesReal_.GetSize() != 0) || (cEdgesInd_.GetSize() != 0) ){
	  cEdgesReal_.Clear(false);
	  cEdgesInd_.Clear(false);
  }


  UInt sizeBH = coarseAuxMat.GetNumRows();

  const UInt *const pRowBH = coarseAuxMat.GetRowPointer();
  const UInt *const pColBH = coarseAuxMat.GetColPointer();

  for(UInt i = 0; i < sizeBH; ++i){
    // get the agglomerate
    Integer c = this->GetIndexOfCoarseNode(i);
    const StdVector<Integer>& agg = this->GetAgglomerateOfNode(c);
    if(agg[0] != c) EXCEPTION("agglomerate-error in prolongation operator");

    Integer realCnumRoot = this->GetRealNodeNumOfCIndex(i);

    // loop over off-diagonals and check if they lie in the same agglomerate
    // pRowBH[i]+1 since we don't want the diagonal
    for(UInt j = pRowBH[i] + 1; j < pRowBH[i + 1]; ++j){
      Integer cI = this->GetIndexOfCoarseNode(pColBH[j]);
      if( !agg.Contains(cI) ){
        StdVector<Integer> cEreal, cEindex;
        cEreal.Resize(2,0);
        cEreal[0] = realCnumRoot;
        cEindex.Resize(2,0);
        cEindex[0] = i;
        Integer realCnum = this->GetRealNodeNumOfCIndex(pColBH[j]);
        cEreal[1] = realCnum;
        cEindex[1] = pColBH[j];
        if( cEreal[1] < cEreal[0] ) cEreal.Swap(0,1);
        cEdgesReal_.Push_back(cEreal);
        cEdgesInd_.Push_back(cEindex);
      }
    }
  }

  // delete twin-edges and re-orient wrong coarse edges
  this->OrientUniqueCEdges();

  numCEdges_ = cEdgesReal_.GetSize();

  return numCEdges_;
}

/**********************************************************/

template <typename T>
bool Agglomerate<T>::compareEdge0(const StdVector< Integer>& a, const StdVector< Integer>& b){
  return (a[0] < b[0]);
}

template <typename T>
bool Agglomerate<T>::compareEdge1(const StdVector< Integer>& a, const StdVector< Integer>& b){
  return (a[1] < b[1]);
}

/**********************************************************/



template <typename T>
void Agglomerate<T>::OrientUniqueCEdges(){
  for(UInt i = 0; i < cEdgesReal_.GetSize(); ++i){
    StdVector<Integer>& cEreal = cEdgesReal_[i];
    StdVector<Integer>& cEindex = cEdgesInd_[i];
    if(cEreal[0] > cEreal[1]){
      cEreal.Swap(0,1);
      cEindex.Swap(0,1);
    }
  }



  // clean twins (equivalent of Matlab's unique()
//TODO with C++11 there is a more elegant way, using lambda-functions!!!
  std::sort(cEdgesReal_.Begin(), cEdgesReal_.End(), compareEdge1);


  StdVector<UInt> found;
  for(UInt i = 0; i < cEdgesReal_.GetSize(); ++i){
    StdVector<Integer>& c = cEdgesReal_[i];
    for(UInt j = i + 1 ; j < cEdgesReal_.GetSize(); ++j){
      StdVector<Integer>& ct = cEdgesReal_[j];
      if( c == ct) {
        found.Push_back(j);
        break;
      }
    }
  }

  StdVector< StdVector< Integer> > temp;
  for(UInt i = 0; i < cEdgesReal_.GetSize(); ++i){
    if( !found.Contains(i) ) temp.Push_back(cEdgesReal_[i]);
  }

  std::sort(temp.Begin(), temp.End(), compareEdge1);


  // now write back into a cleared cEdgesReal_ and cEdgesInd_ vector
  cEdgesReal_.Clear(false);
  cEdgesInd_.Clear(false);
  cEdgesReal_.Resize(temp.GetSize());
  cEdgesInd_.Resize(temp.GetSize());
  for(UInt i = 0; i  < temp.GetSize(); ++i){
    StdVector<Integer>& t = temp[i];
    StdVector<Integer>& cR = cEdgesReal_[i];
    StdVector<Integer>& cI = cEdgesInd_[i];
    cR.Resize(2);
    cI.Resize(2);
    for(UInt j = 0; j < 2; ++j){
      cR[j] = t[j];
      Integer cn = this->GetIndexOfRealNode(t[j]);
      Integer cnt = this->GetCoarseIndex(cn);
      if(cnt != -1){
        cI[j] = cnt;
      }else{
        EXCEPTION("Coarse Edge Error");
      }

    }
  }
}
/**********************************************************/



template <typename T>
void Agglomerate<T>::CreateCoarseRealNodeNum(){
  // we have Size_H_ coarse nodes
  coarseRealNN_.Resize(Size_H_, -1);
  for(UInt i = 0; i < Agglomerates_.GetSize(); ++i){
    StdVector<Integer>& a = Agglomerates_[i];
    UInt cI = this->GetCoarseIndex(a[0]);
    // get the real node number of cI
    coarseRealNN_[cI] = nodeNumIndex_[a[0]];
  }

}

/**********************************************************/

template <typename T>
void Agglomerate<T>::FillIndexNodeNum(){

  for(UInt i = 0; i < nodeNumIndex_.GetSize(); ++i){
    indexNodeNum_[ nodeNumIndex_[i] ] = i;
  }


}

/**********************************************************/

template <typename T>
void Agglomerate<T>::Reset()
{
    // reset neighbourhood pointers to NULL
    NhStartIndex_ = NULL;
    NhEdges_      = NULL;
    // delete array of C-F-Splitting
    //delete [] ( CoarseIndex_ );  CoarseIndex_  = NULL;
    //CoarseIndex_ = NULL;
    
    Size_h_ = 0;
    Size_H_ = 0;
}


} // namespace CoupledField
