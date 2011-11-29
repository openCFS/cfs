
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefFunctionConst.hh"
#include "Domain/CoefFunction/CoefFunctionExpression.hh"
namespace CoupledField{

//! Generate scalar-valued coefficient function
shared_ptr<CoefFunction> 
CoefFunction::Generate( Global::ComplexPart format, 
                        const std::string& realVal, 
                        const std::string& imagVal  ) {

  MathParser* mp = domain->GetMathParser();
  MathParser::HandleType handle = mp->GetNewHandle();

  shared_ptr<CoefFunction> ret;

  bool isVar = false;
  isVar |= ExprDependsOnTimeFreq(realVal);
  isVar |= ExprDependsOnSpace(realVal);
  if( format == Global::REAL) {
    
    // === REAL CASE ===
    if( !isVar ) {
      // --- a) constant ---
      mp->SetExpr(handle,realVal);
      shared_ptr<CoefFunctionConst<Double> > c 
      (new CoefFunctionConst<Double>());
      c->SetScalar(mp->Eval(handle));
      ret = c;
    } else {
      // --- b) variable ---
      shared_ptr<CoefFunctionExpression<Double> > c 
      (new CoefFunctionExpression<Double>());
      c->SetScalar(realVal);
      ret = c;
    }

  } else {
    
    // === COMPLEX CASE ===
    assert(imagVal != std::string(""));
    isVar |= ExprDependsOnTimeFreq(imagVal);
    isVar |= ExprDependsOnSpace(imagVal);
    if( !isVar ) {
      // --- a) constant ---
      mp->SetExpr(handle,realVal);
      Double real  = mp->Eval(handle);
      mp->SetExpr(handle,imagVal);
      Double imag  = mp->Eval(handle);
      shared_ptr<CoefFunctionConst<Complex> > c 
      (new CoefFunctionConst<Complex>());
      c->SetScalar(Complex(real, imag));
      ret = c;
    } else {
      // --- b) variable ---
      shared_ptr<CoefFunctionExpression<Complex> > c 
      (new CoefFunctionExpression<Complex>());
      c->SetScalar(realVal,imagVal);
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

  bool isVar = false;
  
  // Loop over all entries
  for( UInt i = 0; i < realVal.GetSize(); ++i ) {
    isVar |= ExprDependsOnTimeFreq(realVal[i]);
    isVar |= ExprDependsOnSpace(realVal[i]);
  }
  if( format == Global::REAL) {
    
    // === REAL CASE ===
    if( !isVar ) {
      // --- a) constant ---
      mp->SetExpr(handle,realVal.Serialize(','));
      shared_ptr<CoefFunctionConst<Double> > c 
      (new CoefFunctionConst<Double>());
      Vector<Double> realVec;
      mp->EvalVector(handle, realVec);
      c->SetVector(realVec);
      ret = c;
    } else {
      // --- b) variable ---
      shared_ptr<CoefFunctionExpression<Double> > c 
      (new CoefFunctionExpression<Double>());
      c->SetVector(realVal);
      ret = c;
    }

  } else {
    
    // === COMPLEX CASE ===
    assert(imagVal.GetSize() > 0);
    for( UInt i = 0; i < imagVal.GetSize(); ++i ) {
      isVar |= ExprDependsOnTimeFreq(imagVal[i]);
      isVar |= ExprDependsOnSpace(imagVal[i]);
    }
    if( !isVar ) {
      // --- a) constant ---
      mp->SetExpr(handle,realVal.Serialize(','));
      Double real  = mp->Eval(handle);
      mp->SetExpr(handle,imagVal.Serialize(','));
      Double imag  = mp->Eval(handle);
      shared_ptr<CoefFunctionConst<Complex> > c 
      (new CoefFunctionConst<Complex>());
      c->SetScalar(Complex(real, imag));
      ret = c;
    } else {
      // --- b) variable ---
      shared_ptr<CoefFunctionExpression<Complex> > c 
      (new CoefFunctionExpression<Complex>());
      c->SetVector(realVal,imagVal);
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

   bool isVar = false;
   
   // Loop over all entries
   for( UInt i = 0; i < realVal.GetSize(); ++i ) {
     isVar |= ExprDependsOnTimeFreq(realVal[i]);
     isVar |= ExprDependsOnSpace(realVal[i]);
   }
   if( format == Global::REAL) {
     
     // === REAL CASE ===
     if( !isVar ) {
       // --- a) constant ---
       mp->SetExpr(handle,realVal.Serialize(','));
       shared_ptr<CoefFunctionConst<Double> > c 
       (new CoefFunctionConst<Double>());
       Matrix<Double> mat;
       mp->EvalMatrix(handle, mat, numRows, numCols );
       c->SetTensor(mat);
       ret = c; 
     } else {
       // --- b) variable ---
       shared_ptr<CoefFunctionExpression<Double> > c 
       (new CoefFunctionExpression<Double>());
       c->SetTensor(realVal, numRows, numCols );
       ret = c;
     }

   } else {
     
     // === COMPLEX CASE ===
     assert(imagVal.GetSize() > 0);
     for( UInt i = 0; i < imagVal.GetSize(); ++i ) {
       isVar |= ExprDependsOnTimeFreq(imagVal[i]);
       isVar |= ExprDependsOnSpace(imagVal[i]);
     }
     if( !isVar ) {
       // --- a) constant ---
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
     } else {
       // --- b) variable ---
       shared_ptr<CoefFunctionExpression<Complex> > c 
       (new CoefFunctionExpression<Complex>());
       c->SetTensor(realVal, imagVal, numRows, numCols );
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
}
