// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id: AFWsmoother.cc 15673 2017-05-12 19:46:39Z kroppert $ */

#include "OLAS/multigrid/AFWsmoother.hh"
#include <unordered_map>

namespace CoupledField {
/**********************************************************/

template <typename T>
AFWSmoother<T>::AFWSmoother()
    : Size_( 0 ),
      Omega_( 1.0 )
{
}

/**********************************************************/

template <typename T>
AFWSmoother<T>::~AFWSmoother()
{
    Reset();
}

/**********************************************************/

template <typename T>
bool AFWSmoother<T>::Setup( const CRS_Matrix<T>& SrcMat )
{

  this->ExtractPatches(SrcMat);

  // set flag for the prepared state
  return this->prepared_ = true;
}


/**********************************************************/
template <typename T>
void AFWSmoother<T>::
CreatePatches(const CRS_Matrix<Double>& AuxMatrix,
               const StdVector< StdVector< Integer> >& edgeIndNode,
               const StdVector<Integer> nodeNumIndex){

  SizeNodes_ = AuxMatrix.GetNumRows();
  // size of the system matrix is equal to number of edges
  Size_ = edgeIndNode.GetSize();
  Patches_.Resize(SizeNodes_);

  // now we have to loop over every node and store the connected edges (indices in sysMat)
  // NOTE: index of nodeNumIndex is index of the node in the auxiliary-matrix and
  // value of nodeNumIndex is the real node-number, used in edgeIndNode

  // "invert" edgeIndNode, so that we get indices in edgeIndNode for every node
  std::unordered_map<Integer, StdVector<UInt> > edgesOfNode;
  for(UInt i = 0; i < edgeIndNode.GetSize(); ++i ){
    StdVector<Integer> nodes = edgeIndNode[i];
    StdVector<UInt>& t1 = edgesOfNode[nodes[0]];
    StdVector<UInt>& t2 = edgesOfNode[nodes[1]];
    t1.Push_back( i );
    t2.Push_back( i );
  }

  // loop over every node and look for this node in edgeIndNode
  //UInt eIndNodeSize = edgeIndNode.GetSize();
  for(UInt n = 0; n < nodeNumIndex.GetSize(); ++n ){
    Integer nNum = nodeNumIndex[n]; // real nodeNum
    // connected edges of node nNum
    StdVector<UInt>::iterator eIt = edgesOfNode[nNum].Begin();
    while( eIt != edgesOfNode[nNum].End() ){
      Patches_[n].Push_back(*eIt);
      ++eIt;
    }
  }
}

/**********************************************************/


template <typename T>
void AFWSmoother<T>::
ExtractPatches( const CRS_Matrix<T>& SysMat){
  ExtMat_.Resize(SizeNodes_);
  InvExtMat_.Resize(SizeNodes_);

  if(SysMat.GetCurrentLayout() != SysMat.LEX_DIAG_FIRST){
    EXCEPTION("AFW-Smoother: system matrix expected to be LEX_DIAG_FIRST");
  }

  pRow_ = SysMat.GetRowPointer();
  pCol_ = SysMat.GetColPointer();
  pDat_ = SysMat.GetDataPointer();


  T entry;
  // theoretically we construct an extraction matrix R and
  // perform R*K*R^T, but we do it the following way
  for(UInt n = 0; n < SizeNodes_; ++n){
    // vector containing the edge numbers (indices
    // in the system matrix), connected to node n (index in auxiliary-matrix)
    StdVector<UInt>& eNums = Patches_[n];
    UInt s = eNums.GetSize();

    Matrix<T> sM;
    sM.Resize(s,s);
    for(UInt i = 0; i < s; ++i) {
      for(UInt j = i; j < s  ; ++j){
        if( eNums[i] == eNums[j] ){
          for(UInt p = pRow_[eNums[i]]; p < pRow_[eNums[i]+1]; ++p){
		  //handle diagonal extra, since CRS_Matrix return some wrong results..
            if(eNums[i] == pCol_[p]){
              sM[i][i] = pDat_[p];
              break;
            }
          }
        }else{
          if(SysMat.HasMatrixEntry(eNums[i], eNums[j], entry)){
            sM[i][j] = entry;
            if( i != j) sM[j][i] = entry;
          }else{
            sM[i][j] = 0.0;
            sM[j][i] = 0.0;
          }
      }
      }
    }
    ExtMat_[n] = sM;
    // also store the inverse patch-matrix
    //Matrix<T> t;
    //sM.Invert(t);
    sM.Invert_Lapack();
    InvExtMat_[n] = sM;
  }

}

/**********************************************************/



template <typename T>
void AFWSmoother<T>::
Step( const CRS_Matrix<T>&                  matrix,
      const Vector<T>&                      rhs,
            Vector<T>&                      sol,
      const typename Smoother<T>::Direction direction,
      const bool                            force_setup )
{

    // call Setup, if object is not prepared or
    // preparation forced by force_setup == true
    if( !(this->IsPrepared()) || force_setup ) {
        if( false == Setup(matrix) )  return;
    }

    Vector<T> tmp1, tmp2, nSol, exRhs, exSol;
    T row = 0.0;

    UInt pS = 0;
    //NOTE: SizeNodes_ is the same as Patches_.GetSize()
    for(UInt n = 0; n < SizeNodes_; ++n){
      StdVector<UInt>& p = Patches_[n];
      pS = p.GetSize();

      //extract rhs and sol
      exRhs = rhs.GetEntries(p);
      exSol = sol.GetEntries(p);

      //multiply only columns of Kh with u_sol, which are
      //changed, according to Patches_[n]
      tmp1.Resize(pS);


      for(UInt i = 0; i < pS; ++i){
        row = matrix.MultColumnWithVec(p[i], sol);
        tmp1[i] = exRhs[i] - row;
      }

      tmp2 = InvExtMat_[n] * tmp1;
      nSol = exSol + tmp2;

      //insert it back into sol
      for(UInt i = 0; i < pS; ++i){
        sol[p[i]] = nSol[i];
      }
    }

}

/**********************************************************/

template <typename T>
void AFWSmoother<T>::Reset()
{
    Size_            =    0;         // reset the size of the LES
    Omega_           =  1.0;         // reset damping factor

    // call Reset() of base class
    Smoother<T>::Reset();
}

/**********************************************************/
} // namespace CoupledField
