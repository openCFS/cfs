// =====================================================================================
// 
//       Filename:  baseOperator.hh
// 
//    Description:  Basic Class for differential operators. These Classes are passed as a 
//                  template parameter to the forms object. Thereby a PDE rather specifies
//                  an operator than a specific form.
//                  The operator Matrix itself is always a Double.
//                  Are there complex valued functions?
// 
//        Version:  1.0
//        Created:  09/30/2011 10:26:33 AM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#ifndef BaseBOperator_HH
#define BaseBOperator_HH


#include <string>
#include "MatVec/matrix.hh"

namespace CoupledField{
  template<class FE, class DATA_TYPE>
  class BaseBOperator{
    public:
      BaseBOperator(){

      }

      ~BaseBOperator(){

      }

      virtual void CalcOpMat(Matrix<Double> & bMat,LocPointMapped& lp, BaseFE* ptFe ) = 0;

      virtual void CalcOpMatTransposed(Matrix<Double> & bMat,LocPointMapped& lp, BaseFE* ptFe ) = 0;

      virtual void CalcOpMat(Matrix<Complex> & bMat,LocPointMapped& lp, BaseFE* ptFe ){
        Matrix<Double> realMat;
        this->CalcOpMat(realMat,lp,ptFe);
        UInt nrow = realMat.GetNumRows();
        UInt ncol = realMat.GetNumCols();
        bMat.Resize(nrow,ncol);
        bMat.Init();
        for(UInt i=0;i<nrow;++i){
          for(UInt j=0;j<ncol;++j){
            bMat[i][j] = Complex(1.0,0.0) * realMat[i][j];
          }
        }
      }

      virtual void CalcOpMatTransposed(Matrix<Complex> & bMat,LocPointMapped& lp, BaseFE* ptFe ){
          Matrix<Double> realMat;
          this->CalcOpMatTransposed(realMat,lp,ptFe);
          UInt nrow = realMat.GetNumRows();
          UInt ncol = realMat.GetNumCols();
          bMat.Resize(nrow,ncol);
          bMat.Init();
          for(UInt i=0;i<nrow;++i){
            for(UInt j=0;j<ncol;++j){
              bMat[i][j] = Complex(1.0,0.0) * realMat[i][j];
            }
          }
        }

      virtual void ApplyOp(Vector<DATA_TYPE>& retVec,
                            LocPointMapped& lp,
                            BaseFE* ptFe,
                            const Vector<DATA_TYPE>& solVec ){
        Matrix<DATA_TYPE> bOp;
        CalcOpMat(bOp,lp,ptFe);
        retVec = bOp * solVec;

      }

      virtual void ApplyOpTranspose(Vector<DATA_TYPE>& retVec,
                                        LocPointMapped& lp,
                                        BaseFE* ptFe,
                                        const Vector<DATA_TYPE>& solVec ){
        Matrix<DATA_TYPE> bOp;
        CalcOpMat(bOp,lp,ptFe);
        retVec = Transpose(bOp) * solVec;

      }

      UInt GetDiffOrder(){
        return diffOrder_;
      }

      std::string GetName(){
        return name_;
      }

      virtual void TransformJacDet(DATA_TYPE & jacLoc,LocPointMapped & lp, BaseFE* ptFe){
        return;
      }

      void SetSolDim(UInt dim){
        dim_ = dim;
      }

    protected:
      UInt diffOrder_;
      std::string name_;
      UInt dim_;
  };
}
#endif
