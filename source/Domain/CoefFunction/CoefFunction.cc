
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefFunctionConst.hh"
#include "Domain/CoefFunction/CoefFunctionExpression.hh"
#include "Domain/CoefFunction/CoefFunctionTimeFreq.hh"
namespace CoupledField{

//! Generate scalar-valued coefficient function
shared_ptr<CoefFunction> 
CoefFunction::Generate( Global::ComplexPart format, 
                        const std::string& realVal, 
                        const std::string& imagVal  ) {

  MathParser* mp = domain->GetMathParser();
  MathParser::HandleType handle = mp->GetNewHandle();

  shared_ptr<CoefFunction> ret;

  bool depTime = ExprDependsOnTimeFreq(realVal);
  bool depSpace = ExprDependsOnSpace(realVal);
  if( format == Global::REAL) {
    // === REAL CASE ===
    if( depSpace ) {
      // --- a) general case: expression 
      shared_ptr<CoefFunctionExpression<Double> > c 
      (new CoefFunctionExpression<Double>());
      c->SetScalar(realVal);
      ret = c;
    } else if (depTime) {
      // --- b) time / freq dependency 
      shared_ptr<CoefFunctionTimeFreq<Double> > c 
      (new CoefFunctionTimeFreq<Double>());
      c->SetScalar(realVal);
      ret = c;
    } else { 
      // --- c) constant ---
      mp->SetExpr(handle,realVal);
      shared_ptr<CoefFunctionConst<Double> > c 
      (new CoefFunctionConst<Double>());
      c->SetScalar(mp->Eval(handle));
      ret = c;
    }      
  } else {
    
    // === COMPLEX CASE ===
    assert(imagVal != std::string(""));
    depTime |= ExprDependsOnTimeFreq(imagVal);
    depSpace |= ExprDependsOnSpace(imagVal);
    if( depSpace ) {
      // --- a) variable ---
      shared_ptr<CoefFunctionExpression<Complex> > c 
      (new CoefFunctionExpression<Complex>());
      c->SetScalar(realVal,imagVal);
      ret = c;
    } else if( depTime ) {
      // --- b) time / freq dependency 
      shared_ptr<CoefFunctionTimeFreq<Complex> > c 
      (new CoefFunctionTimeFreq<Complex>());
      c->SetScalar(realVal,imagVal);
      ret = c;
    } else {
      // --- c) constant ---
      mp->SetExpr(handle,realVal);
      Double real  = mp->Eval(handle);
      mp->SetExpr(handle,imagVal);
      Double imag  = mp->Eval(handle);
      shared_ptr<CoefFunctionConst<Complex> > c 
      (new CoefFunctionConst<Complex>());
      c->SetScalar(Complex(real, imag));
      ret = c;
    }
  }
  mp->ReleaseHandle(handle);
  return ret;

}

//! Generate vector-valued coefficient function
shared_ptr<CoefFunction> 
CoefFunction::Generate( Global::ComplexPart format, 
                        const StdVector<std::string>& realVal, 
                        const StdVector<std::string>& imagVal ) {
  
  MathParser* mp = domain->GetMathParser();
  MathParser::HandleType handle = mp->GetNewHandle();

  shared_ptr<CoefFunction> ret;

  bool depTime = false;
  bool depSpace = false;
  
  // Loop over all entries
  for( UInt i = 0; i < realVal.GetSize(); ++i ) {
    depTime  |= ExprDependsOnTimeFreq(realVal[i]);
    depSpace |= ExprDependsOnSpace(realVal[i]);
  }
  if( format == Global::REAL) {

    // === REAL CASE ===
    if( depSpace ) {
      // --- a) general case: expression 
      shared_ptr<CoefFunctionExpression<Double> > c 
      (new CoefFunctionExpression<Double>());
      c->SetVector(realVal);
      ret = c;
    } else if (depTime) {
      // --- b) time / freq dependency 
      shared_ptr<CoefFunctionTimeFreq<Double> > c 
      (new CoefFunctionTimeFreq<Double>());
      c->SetVector(realVal);
      ret = c;
    } else {
      // --- c) constant ---
      mp->SetExpr(handle,realVal.Serialize(','));
      shared_ptr<CoefFunctionConst<Double> > c 
      (new CoefFunctionConst<Double>());
      Vector<Double> realVec;
      mp->EvalVector(handle, realVec);
      c->SetVector(realVec);
      ret = c;
    }

  } else {
    
    // === COMPLEX CASE ===
    assert(imagVal.GetSize() > 0);
    assert(imagVal.GetSize() == realVal.GetSize());
    
    for( UInt i = 0; i < imagVal.GetSize(); ++i ) {
      depTime  |= ExprDependsOnTimeFreq(imagVal[i]);
      depSpace |= ExprDependsOnSpace(imagVal[i]);
    }
    if( depSpace ) {
      // --- a) general case: expression
      shared_ptr<CoefFunctionExpression<Complex> > c 
      (new CoefFunctionExpression<Complex>());
      c->SetVector(realVal,imagVal);
      ret = c;
    } else if (depTime) {
      // --- b) time / freq dependency
      shared_ptr<CoefFunctionTimeFreq<Complex> > c 
      (new CoefFunctionTimeFreq<Complex>());
      c->SetVector(realVal,imagVal);
      ret = c;
    } else {
      shared_ptr<CoefFunctionConst<Complex> > c 
            (new CoefFunctionConst<Complex>());
      Vector<Double> realVec, imagVec;
      mp->SetExpr(handle,realVal.Serialize(','));
      mp->EvalVector(handle, realVec);
      mp->SetExpr(handle,imagVal.Serialize(','));
      mp->EvalVector(handle, imagVec);
      
      // assemble complex vector
      Vector<Complex> cVec(realVec.GetSize());
      cVec.SetPart(Global::REAL, realVec );
      cVec.SetPart(Global::IMAG, imagVec );
      c->SetVector(cVec);
      ret = c;
    } 
  }
  mp->ReleaseHandle(handle);
  return ret;

}


//! Generate tensor-valued coefficient function
shared_ptr<CoefFunction> 
CoefFunction::Generate( Global::ComplexPart format,
                        UInt numRows, UInt numCols,
                        const StdVector<std::string>& realVal,
                        const StdVector<std::string>& imagVal ) {
  MathParser* mp = domain->GetMathParser();
   MathParser::HandleType handle = mp->GetNewHandle();

   shared_ptr<CoefFunction> ret;

   bool depTime = false;
   bool depSpace = false;
   // Loop over all entries
   for( UInt i = 0; i < realVal.GetSize(); ++i ) {
     depTime  |= ExprDependsOnTimeFreq(realVal[i]);
     depSpace |= ExprDependsOnSpace(realVal[i]);
   }
   if( format == Global::REAL) {
     
     // === REAL CASE ===
     if( depSpace ) {
       // --- a) general case: expression 
       shared_ptr<CoefFunctionExpression<Double> > c 
       (new CoefFunctionExpression<Double>());
       c->SetTensor(realVal, numRows, numCols );
       ret = c;
     } else if (depTime) {
       // --- b) time / freq dependency
       shared_ptr<CoefFunctionTimeFreq<Double> > c 
       (new CoefFunctionTimeFreq<Double>());
       c->SetTensor(realVal, numRows, numCols );
       ret = c;
     } else {
       // --- c) constant ---
       mp->SetExpr(handle,realVal.Serialize(','));
       shared_ptr<CoefFunctionConst<Double> > c 
       (new CoefFunctionConst<Double>());
       Matrix<Double> mat;
       mp->EvalMatrix(handle, mat, numRows, numCols );
       c->SetTensor(mat);
       ret = c; 
     }
   } else {
     
     // === COMPLEX CASE ===
     assert(imagVal.GetSize() > 0);
     assert(imagVal.GetSize() == realVal.GetSize());

     for( UInt i = 0; i < imagVal.GetSize(); ++i ) {
       depTime  |= ExprDependsOnTimeFreq(imagVal[i]);
       depSpace |= ExprDependsOnSpace(imagVal[i]);
     }
     if( depSpace ) {
       // --- a) general case: expression
       shared_ptr<CoefFunctionExpression<Complex> > c 
       (new CoefFunctionExpression<Complex>());
       c->SetTensor(realVal, imagVal, numRows, numCols );
       ret = c;
     } else if (depTime) {
       // --- b) time / freq dependency
       shared_ptr<CoefFunctionTimeFreq<Complex> > c 
       (new CoefFunctionTimeFreq<Complex>());
       c->SetTensor(realVal, imagVal, numRows, numCols );
       ret = c;
     } else {
       Matrix<Double> realPart, imagPart;
       Matrix<Complex> cMat(numRows,numCols);
       mp->SetExpr(handle,realVal.Serialize(','));
       mp->EvalMatrix(handle, realPart, numRows, numCols );

       mp->SetExpr(handle,imagVal.Serialize(','));
       mp->EvalMatrix(handle, imagPart, numRows, numCols );

       // assemble complex matrix
       cMat.SetPart(Global::REAL, realPart);
       cMat.SetPart(Global::IMAG, imagPart);
       shared_ptr<CoefFunctionConst<Complex> > c 
       (new CoefFunctionConst<Complex>());
       c->SetTensor(cMat);
       ret = c;
     }
   }
   mp->ReleaseHandle(handle);


   return ret;

}
bool CoefFunction::ExprDependsOnTimeFreq(const std::string& expr) {
  MathParser* mp = domain->GetMathParser();
  MathParser::HandleType handle = mp->GetNewHandle();
  mp->SetExpr(handle, expr);
  bool depends = false;
  if ( mp->IsExprVariable(handle, "t") ||
       mp->IsExprVariable(handle, "f") ){
    depends = true;
  }
  
  mp->ReleaseHandle(handle);
  return depends;
}

//! Returns true, if expression depends on space
bool CoefFunction::ExprDependsOnSpace(const std::string& expr) {
  MathParser* mp = domain->GetMathParser();
  MathParser::HandleType handle = mp->GetNewHandle();
  mp->SetExpr(handle, expr);
  bool depends = false;
  if ( mp->IsExprVariable(handle, "x") ||
       mp->IsExprVariable(handle, "y") ||
       mp->IsExprVariable(handle, "z") ){
    depends = true;
  }
  
  mp->ReleaseHandle(handle);
  return depends;
}

// ************************************************************************
// ENUM INITIALIZATION
// ************************************************************************

// Definition of coefficient function dimension type
static EnumTuple dimTypeTuples[] = 
{
 EnumTuple(CoefFunction::NO_DIM,  "NO_DIM"), 
 EnumTuple(CoefFunction::SCALAR,  "SCALAR"),
 EnumTuple(CoefFunction::VECTOR,  "VECTOR"),
 EnumTuple(CoefFunction::TENSOR,  "TENSOR")
};

Enum<CoefFunction::CoefDimType> CoefFunction::CoefDimType_ = \
    Enum<CoefFunction::CoefDimType>("Dimension of CoefFunction",
                                    sizeof(dimTypeTuples) / sizeof(EnumTuple),
                                    dimTypeTuples);

// Definition of coefficient function dependency type
static EnumTuple coefDependTuples[] = 
{
 EnumTuple(CoefFunction::CONST,    "CONST"), 
 EnumTuple(CoefFunction::TIMEFREQ, "TIMEFREQ"),
 EnumTuple(CoefFunction::GENERAL,  "GENERAL")
};

Enum<CoefFunction::CoefDependType>CoefFunction::CoefDependType_ = \
    Enum<CoefFunction::CoefDependType>("Dependency of CoefFunction",
                                       sizeof(coefDependTuples) / sizeof(EnumTuple),
                                       coefDependTuples);



}
