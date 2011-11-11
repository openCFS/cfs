// =====================================================================================
// 
//       Filename:  identityOperator.hh
// 
//    Description:  This class implements the identity operator just
//                  consisting of the elements shape functions evaulated at the 
//                  given local coordinate
//                      / N_1x  0    0   N_2_x .. \
//                  b = |  0   N_1y  0    0    .. |
//                      \  0    0   N_1z  0    .. /
//
//                  for the basic shape functions N_1x = N_1y = N_1z
//                  for scalar unknowns the operator is just a vector
// 
//        Version:  1.0
//        Created:  10/04/2011 09:10:00 AM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#ifndef IDENTITYOP_HH
#define IDENTITYOP_HH
#include "BaseBOperator.hh"
#include "FeBasis/HCurl/HCurlElems.hh"

namespace CoupledField{
  template<class FE, class TYPE>
  class IdentityOperator : public BaseBOperator<FE,TYPE>{
    public:
      IdentityOperator(){
        return;
      }

      ~IdentityOperator(){
        return;
      }

      virtual void CalcOpMat(Matrix<Double> & bMat,LocPointMapped& lp, BaseFE* ptFe );

      virtual void CalcOpMatTransposed(Matrix<Double> & bMat,LocPointMapped& lp, BaseFE* ptFe );

      //avoid reimplementation of complex operator by making the bas class function
      //available
      using BaseBOperator<FE,TYPE>::CalcOpMat;

      using BaseBOperator<FE,TYPE>::CalcOpMatTransposed;


    protected:

  };
  template< class TYPE>
  class IdentityOperator<FeHCurl,TYPE> : public BaseBOperator<FeHCurl,TYPE>{
    public:
      IdentityOperator(){
        return;
      }

      ~IdentityOperator(){
        return;
      }

      virtual void CalcOpMat(Matrix<Double> & bMat,LocPointMapped& lp, BaseFE* ptFe ){
        FeHCurl* fe = static_cast<FeHCurl*>(ptFe);
        fe->GetShFnc( bMat, lp, lp.shapeMap->GetElem(), 0);
      }

      virtual void CalcOpMatTransposed(Matrix<Double> & bMat,LocPointMapped& lp, BaseFE* ptFe ){
        FeHCurl* fe = static_cast<FeHCurl*>(ptFe);
        Matrix<Double> xiDx;
        fe->GetShFnc( xiDx, lp, lp.shapeMap->GetElem(), 0);
        xiDx.Transpose(bMat);
      }

    protected:

  };

  template<class FE,class TYPE>
  class ScaledByEdgeIdentityOperator : public IdentityOperator<FE,TYPE>{
    public:
      ScaledByEdgeIdentityOperator(){
        return;
      }

      ~ScaledByEdgeIdentityOperator(){
        return;
      }

      virtual void CalcOpMat(Matrix<Double> & bMat,LocPointMapped& lp, BaseFE* ptFe ){
        IdentityOperator<FeHCurl,TYPE>::CalcOpMat(bMat,lp,ptFe);
        //scale By Edge
        Double minE,maxE;
        lp.shapeMap->GetMaxMinEdgeLength(maxE,minE);
        bMat /= maxE;
      }

      virtual void CalcOpMatTransposed(Matrix<Double> & bMat,LocPointMapped& lp, BaseFE* ptFe ){
        IdentityOperator<FeHCurl,TYPE>::CalcOpMatTransposed(bMat,lp,ptFe);
        //scale By Edge
        Double minE,maxE;
        lp.shapeMap->GetMaxMinEdgeLength(maxE,minE);
        bMat /= maxE;
      }

  };

  template<class FE, class TYPE>
  void IdentityOperator<FE,TYPE>::CalcOpMat(Matrix<Double> & bMat,LocPointMapped& lp, BaseFE* ptFe){
    UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( this->dim_, numFncs * this->dim_ );
    bMat.Init();

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Vector<Double> S;
    FE *fe = (static_cast<FE*>(ptFe));
    for(UInt d = 0; d < this->dim_ ; d ++){
      fe->GetShFnc( S, lp.lp, lp.shapeMap->GetElem() , d );
      for(UInt sh = 0; sh < S.GetSize(); sh ++){
        bMat[d][sh*(this->dim_) + d] = S[sh];
      }
    }

  }

  template<class FE, class TYPE>
  void IdentityOperator<FE,TYPE>::CalcOpMatTransposed(Matrix<Double> & bMat,LocPointMapped& lp, BaseFE* ptFe){
    UInt numFncs = ptFe->GetNumFncs();
    // Set correct size of matrix B and initialise with zeros
    bMat.Resize( numFncs * this->dim_ , this->dim_ );
    bMat.Init();

    // Get derivatives of local shape functions with respect to global
    // coords (format: nrNodes x spaceDim)
    Vector<Double> S;
    FE *fe = (static_cast<FE*>(ptFe));
    for(UInt d = 0; d < this->dim_ ; d ++){
      fe->GetShFnc( S, lp.lp, lp.shapeMap->GetElem() , d );
      for(UInt sh = 0; sh < S.GetSize(); sh ++){
        bMat[sh*(this->dim_) + d][d] = S[sh];
      }
    }

  }
}
#endif
