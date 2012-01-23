// =====================================================================================
// 
//       Filename:  curlOp.hh
// 
//    Description:  This operator implements the curl of the shape function
//                  Its operator-matrix has the structure
//                      / d N_1z/dy - d N_1y/dz ...\
//                  B = | d N_1x/dz - d N_1z/dx ...|
//                      \ d N_1y/dx - d N_1x/dy .../
//                      
//        Version:  1.0
//        Created:  10/28/2011 05:35:08 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================

#include "BaseBOperator.hh"

namespace CoupledField{
  template<class FE, UInt D, class TYPE>
  class CurlOperator : public BaseBOperator<FE,TYPE>{
    public:
    
    CurlOperator();

      ~CurlOperator();

      virtual void CalcOpMat(Matrix<Double> & bMat,
                             const LocPointMapped& lp, BaseFE* ptFe ){
        EXCEPTION("Curl Operator not implemented for general functions yet")
      }

      virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                       const LocPointMapped& lp, BaseFE* ptFe ){
        EXCEPTION("Curl Operator not implemented for general functions yet")
      }

      //avoid reimplementation of complex operator by making the bas class function
      //available, template caused
      using BaseBOperator<FE,TYPE>::CalcOpMat;

      using BaseBOperator<FE,TYPE>::CalcOpMatTransposed;

    protected:

  };

  //!Specialized class for HCurl elements
  template<UInt D, class TYPE >
  class CurlOperator<FeHCurl,D,TYPE> : public BaseBOperator<FeHCurl,TYPE>{
    public:
    
    // ------------------
    //  STATIC CONSTANTS 
    // ------------------
    //@{ 
    //! \name Static constants

    //! Order of differentiation
    static const UInt ORDER_DIFF = 1;

    //! Number of components of the problem (scalar, vector)
    static const UInt DIM_DOF = 1;

    //! Dimension of the underlying domain / space
    static const UInt DIM_SPACE = D;

    //! Dimension of the finite element
    static const UInt DIM_ELEM = D;

    //! Dimension of the related material 
    static const UInt DIM_D_MAT = D; 
    //@}
    
      CurlOperator(){

      }

      ~CurlOperator(){

      }

      virtual void CalcOpMat(Matrix<Double> & bMat,
                             const LocPointMapped& lp, BaseFE* ptFe ){
        FeHCurl *fe = (static_cast<FeHCurl*>(ptFe));
        fe->GetCurlShFnc(bMat, lp, lp.shapeMap->GetElem(), 1);
      }

      virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                       const LocPointMapped& lp, BaseFE* ptFe ){
        FeHCurl *fe = (static_cast<FeHCurl*>(ptFe));
        Matrix<Double> xiDx;
        fe->GetCurlShFnc(xiDx, lp, lp.shapeMap->GetElem(), 1);
        xiDx.Transpose(bMat);
      }

      //avoid reimplementation of complex operator by making the bas class function
      //available, template caused
      using BaseBOperator<FeHCurl,TYPE>::CalcOpMat;

      using BaseBOperator<FeHCurl,TYPE>::CalcOpMatTransposed;

    protected:

  };

  //!Specialized class for curl elements
  template<UInt D, class TYPE>
  class ScaledByEdgeCurlOperator : public CurlOperator<FeHCurl,D,TYPE>{
    public:
    
      ScaledByEdgeCurlOperator(){

      }

      ~ScaledByEdgeCurlOperator(){

      }

      virtual void CalcOpMat(Matrix<Double> & bMat,
                             const LocPointMapped& lp, BaseFE* ptFe ){
        CurlOperator<FeHCurl,D,TYPE>::CalcOpMat(bMat,lp,ptFe);
        Double minE,maxE;
        lp.shapeMap->GetMaxMinEdgeLength(maxE,minE);
        bMat /= maxE;
      }

      virtual void CalcOpMatTransposed(Matrix<Double> & bMat,
                                       const LocPointMapped& lp, BaseFE* ptFe ){
        CurlOperator<FeHCurl,D,TYPE>::CalcOpMatTransposed(bMat,lp,ptFe);
        Double minE,maxE;
        lp.shapeMap->GetMaxMinEdgeLength(maxE,minE);
        bMat /= maxE;
      }

      //avoid reimplementation of complex operator by making the bas class function
      //available, template caused
      using BaseBOperator<FeHCurl,TYPE>::CalcOpMat;

      using BaseBOperator<FeHCurl,TYPE>::CalcOpMatTransposed;

    protected:

  };
}
