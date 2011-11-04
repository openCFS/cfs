// =====================================================================================
// 
//       Filename:  coefFunctionExpression.hh
// 
//    Description:  This coefficient funciton handles coefficients based on math parser
//                  expressions. Therefore, the set Matrix/scalar/vector function expect
//                  a string rather then a Double or Complex
// 
//        Version:  1.0
//        Created:  11/03/2011 12:17:37 PM
//       Revision:  none
//       Compiler:  g++
// 
//         Author:  Andreas Hueppe (AHU), andreas.hueppe@uni-klu.ac.at
//        Company:  Universitaet Klagenfurt
// 
// =====================================================================================


#ifndef COEFFUNCTIONEXPRESSION_HH
#define COEFFUNCTIONEXPRESSION_HH

#include "coefFunction.hh"

#include "General/environment.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "Domain/elem.hh"
#include "Elements/elemShapeMap.hh"
#include "Utils/mathParser/mathParser.hh"

namespace CoupledField{

template<class DATA_TYPE>
class CoefFunctionExpression : public CoefFunction<DATA_TYPE>{
  public:
    CoefFunctionExpression():
      mHandle_(domain->GetMathParser()->GetNewHandle(true))
    {

    }

    ~CoefFunctionExpression(){
      domain->GetMathParser()->ReleaseHandle(mHandle_);
    }

    void GetTensor(Matrix<DATA_TYPE>& CoefMat, LocPointMapped lp, const Elem* elem){
      Excpetion("CoefFunctionExpression::GetMatrix not implemented");
    }

    void GetVector(Matrix<DATA_TYPE>& CoefScalar, LocPointMapped lp, const Elem* elem){
      Excpetion("CoefFunctionExpression::GetVector not implemented");
    }

    void GetScalar(Matrix<DATA_TYPE>& CoefScalar, LocPointMapped lp, const Elem* elem){
      Excpetion("CoefFunctionExpression::GetScalar not implemented");
    }

    void SetTensor(Matrix<std::string> val){
      assert((this->mType_ == UNDEF) || (this->mType_ == TENSOR) );
      this->coefMat_ = val;
      this->mType = CoefFunction<DATA_TYPE>::TENSOR;
    }

    void SetVector(Vector<std::string> val){
      assert((this->mType_ == UNDEF) || (this->mType_ == VECTOR) );
      this->coefVec_ = val;
      this->mType = CoefFunction<DATA_TYPE>::VECTOR;
    }

    void SetScalar(std::string val){
      assert((this->mType_ == UNDEF) || (this->mType_ == SCALAR) );
      this->coefScalar_ = val;
      this->mType = CoefFunction<DATA_TYPE>::SCALAR;;
    }

  protected:
    Matrix<std::string> coefMat_;

    Vector<std::string> coefVec_;

    std::string coefScalar_;

    MathParser::HandleType mHandle_;


};

template<>
class CoefFunctionExpression<Double>{
  public:
    void GetVector(Matrix<Double>& CoefVec, LocPointMapped lp, const Elem* elem){
      //evaluate class variable, vector of strings using math paser and return it
      //in future implementations we want to distiguish between expressions that depend
      //only on time or frequency and expressions which include a spatial dependency.
      //The latter have to be evaluated every time while the first have to be evaluated only
      //at the beginning of the time step

      //FIRST IMPLEMENTATION
      //we assume everything has to be evaluated everytime....
      assert(this->mType_ == VECTOR);
      MathParser * parser = domain->GetMathParser();
      Vector<Double> pointCoord;

      CoefVec.Resize(this->coefVec_.GetSize());
      CoefVec.Init();

      lp.shapeMap->Local2Global(pointCoord,lp.lp);
      for(UInt i = 0; i<coefVec_.GetSize();++i){
        parser->SetExpr( mHandle_, coefVec_[i] );
        CoefVec[i] = parser->Eval( mHandle_ );
      }
    }
};


template<>
class CoefFunctionExpression<Complex>{
  public:
    void GetVector(Matrix<Complex>& CoefVec, LocPointMapped lp, const Elem* elem){
      //This is the complex version. Ok right now, mathparser only supports
      //Real valued expression. Therefore we implement this extra. In the future
      //we hope that the mathParserX will enable us to automatically evaluate also
      //real valued expressions as complex numbers....

      //FIRST IMPLEMENTATION
      //we assume everything has to be evaluated everytime....
      assert(this->mType_ == VECTOR);
      MathParser * parser = domain->GetMathParser();
      Vector<Double> pointCoord;
      UInt vSize = this->coefVec_.GetSize();

      CoefVec.Resize(vSize);
      CoefVec.Init();

      lp.shapeMap->Local2Global(pointCoord,lp.lp);
      parser->SetCoordinates( mHandle_, this->coordSys_, pointCoord );

      for(UInt i = 0; i<vSize;++i){
        parser->SetExpr( mHandle_, coefVec_[i] );
        CoefVec[i] = Complex(parser->Eval( mHandle_ ));
      }
    }
};

}
