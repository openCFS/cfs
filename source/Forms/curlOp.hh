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

#include "baseOperator.hh"

namespace CoupledField{
  template<class FE, class TYPE>
  class CurlOperator : public BaseBOperator<FE,TYPE>{
    public:
    CurlOperator();

      ~CurlOperator();

      virtual void CalcOpMat(Matrix<Double> & bMat,LocPointMapped& lp, BaseFE* ptFe ){
        EXCEPTION(__FILE__ << __LINE__ << "Curl Operator not implemented for general functions yet")
      }

      //avoid reimplementation of complex operator by making the bas class function
      //available, template caused
      using BaseBOperator<FE,TYPE>::CalcOpMat;

    protected:

  };

  //!Specialized class for curl elements
  template< class TYPE>
  class CurlOperator<FeHCurl,TYPE> : public BaseBOperator<FeHCurl,TYPE>{
    public:
      CurlOperator(){

      }

      ~CurlOperator(){

      }

      virtual void CalcOpMat(Matrix<Double> & bMat,LocPointMapped& lp, BaseFE* ptFe ){
        FeHCurl *fe = (static_cast<FeHCurl*>(ptFe));
        fe->GetCurlShFnc(bMat, lp, lp.shapeMap->GetElem(), 1);
      }

      //avoid reimplementation of complex operator by making the bas class function
      //available, template caused
      using BaseBOperator<FeHCurl,TYPE>::CalcOpMat;

    protected:

  };

  //!Specialized class for curl elements
  template<class FE, class TYPE>
  class ScaledByEdgeCurlOperator : public CurlOperator<FeHCurl,TYPE>{
    public:
      ScaledByEdgeCurlOperator(){

      }

      ~ScaledByEdgeCurlOperator(){

      }

      virtual void CalcOpMat(Matrix<Double> & bMat,LocPointMapped& lp, BaseFE* ptFe ){
        CurlOperator<FeHCurl,TYPE>::CalcOpMat(bMat,lp,ptFe);
        Double minE,maxE;
        lp.shapeMap->GetMaxMinEdgeLength(maxE,minE);
        bMat /= maxE;
      }

      //avoid reimplementation of complex operator by making the bas class function
      //available, template caused
      using BaseBOperator<FE,TYPE>::CalcOpMat;

    protected:

  };
}
