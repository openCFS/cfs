
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/CoefFunction/CoefFunctionConst.hh"
#include "Domain/CoefFunction/CoefFunctionExpression.hh"
#include "Domain/CoefFunction/CoefFunctionTimeFreq.hh"
#include "Domain/CoefFunction/CoefFunctionCompound.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "Utils/mathParser/mathParser.hh"

namespace CoupledField{

//! Get the maximum CoefFunction dependency type
CoefFunction::CoefDependType CoefFunction::GetMaxCoefDependType(CoefFunction::CoefDependType typeA, CoefFunction::CoefDependType typeB) {
  if (typeA == CoefFunction::SOLUTION || typeB == CoefFunction::SOLUTION) {
    return CoefFunction::SOLUTION;
  }
  if (typeA == CoefFunction::GENERAL || typeB == CoefFunction::GENERAL) {
    return CoefFunction::GENERAL;
  }
  if (typeA == CoefFunction::SPACE) {
    if (typeB == CoefFunction::TIMEFREQ) {
      return CoefFunction::GENERAL;
    }
    return CoefFunction::SPACE;
  } else if (typeA == CoefFunction::TIMEFREQ) {
    if (typeB == CoefFunction::SPACE) {
      return CoefFunction::GENERAL;
    }
    return CoefFunction::TIMEFREQ;
  }
  return typeB;
}  
  
//! Generate scalar-valued coefficient function
PtrCoefFct 
CoefFunction::Generate( MathParser * mp,
                        Global::ComplexPart format, 
                        const std::string& realVal, 
                        const std::string& imagVal  ) {

  unsigned int handle = mp->GetNewHandle(true);


  PtrCoefFct ret;

  bool depTime = ExprDependsOnTimeFreq(mp, realVal);
  bool depSpace = ExprDependsOnSpace(mp, realVal);
  if( format == Global::REAL) {
    // === REAL CASE ===
    if( depSpace ) {
      // --- a) general case: expression 
      shared_ptr<CoefFunctionExpression<Double> > c 
      (new CoefFunctionExpression<Double>(mp));
      c->SetScalar(realVal);
      ret = c;
    } else if (depTime) {
      // --- b) time / freq dependency 
      shared_ptr<CoefFunctionTimeFreq<Double> > c 
      (new CoefFunctionTimeFreq<Double>(mp));
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
    depTime |= ExprDependsOnTimeFreq(mp, imagVal);
    depSpace |= ExprDependsOnSpace(mp, imagVal);
    if( depSpace ) {
      // --- a) variable ---
      shared_ptr<CoefFunctionExpression<Complex> > c 
      (new CoefFunctionExpression<Complex>(mp));
      c->SetScalar(realVal,imagVal);
      ret = c;
    } else if( depTime ) {
      // --- b) time / freq dependency 
      shared_ptr<CoefFunctionTimeFreq<Complex> > c 
      (new CoefFunctionTimeFreq<Complex>(mp));
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

//! Generate vector-valued coefficient function from string
PtrCoefFct 
CoefFunction::Generate( MathParser * mp,
                        Global::ComplexPart format, 
                        const StdVector<std::string>& realVal, 
                        const StdVector<std::string>& imagVal ) {
  
  unsigned int handle = mp->GetNewHandle(true);

  PtrCoefFct ret;

  bool depTime = ExprDependsOnTimeFreq(mp, realVal);
  bool depSpace = ExprDependsOnSpace(mp, realVal);

  if( format == Global::REAL) {

    // === REAL CASE ===
    if( depSpace ) {
      // --- a) general case: expression 
      shared_ptr<CoefFunctionExpression<Double> > c 
      (new CoefFunctionExpression<Double>(mp));
      c->SetVector(realVal);
      ret = c;
    } else if (depTime) {
      // --- b) time / freq dependency 
      shared_ptr<CoefFunctionTimeFreq<Double> > c 
      (new CoefFunctionTimeFreq<Double>(mp));
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
    
    depTime |= ExprDependsOnTimeFreq(mp, imagVal);
    depSpace |= ExprDependsOnSpace(mp, imagVal);
    
    if( depSpace ) {
      // --- a) general case: expression Matrix<MAT_DATA_TYPE> dMat_;
      shared_ptr<CoefFunctionExpression<Complex> > c 
      (new CoefFunctionExpression<Complex>(mp));
      c->SetVector(realVal,imagVal);
      ret = c;
    } else if (depTime) {
      // --- b) time / freq dependency
      shared_ptr<CoefFunctionTimeFreq<Complex> > c 
      (new CoefFunctionTimeFreq<Complex>(mp));
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

//! Generate vector-valued coefficient function from CoefFunctions
PtrCoefFct 
CoefFunction::Generate( MathParser * mp, 
                        Global::ComplexPart format, 
                        const StdVector<PtrCoefFct>& realVal, 
                        const StdVector<PtrCoefFct>& imagVal ) {
  
  PtrCoefFct ret;
  // Check, if all entries are analytical
  bool isAnalytical = false;
  for( UInt i = 0; i < realVal.GetSize(); ++i ) {
    isAnalytical |= realVal[i]->IsAnalytic();
    // also ensure, that every CoefFunction is of scalar type
    if( realVal[i]->GetDimType() != CoefFunction::SCALAR )  {
      EXCEPTION( "Vector valued expression can be only composed of "
          << "scalar entries!" );
    }
  }
  for( UInt i = 0; i < imagVal.GetSize(); ++i ) {
    isAnalytical |= imagVal[i]->IsAnalytic();
    // also ensure, that every CoefFunction is of scalar type
    if( imagVal[i]->GetDimType() != CoefFunction::SCALAR )  {
      EXCEPTION( "Vector valued expression can be only composed of "
          << "scalar entries!" );
    }
  }
  
  if( isAnalytical ) {
    // use string generation method to generate the vector valued CoefFct
    StdVector<std::string> rStrVals, iStrVals;
    std::string empty;
    rStrVals.Resize( realVal.GetSize() );
    iStrVals.Resize( imagVal.GetSize() );
    for( UInt i = 0; i < realVal.GetSize(); ++i ) {
      shared_ptr<CoefFunctionAnalytic> tmp;
      tmp = dynamic_pointer_cast<CoefFunctionAnalytic>(realVal[i]);
      tmp->GetStrScalar(rStrVals[i], empty);
    }
    for( UInt i = 0; i < imagVal.GetSize(); ++i ) {
      shared_ptr<CoefFunctionAnalytic> tmp;
      tmp = dynamic_pointer_cast<CoefFunctionAnalytic>(realVal[i]);
      tmp->GetStrScalar(iStrVals[i], empty);
    }
    ret = Generate(mp, format, rStrVals, iStrVals);
  } else {
    EXCEPTION( "Vector-valued expression can only be composed of "
        << "analytical expressions currently!" );

  }
  return ret;
}


//! Generate tensor-valued coefficient function from strings
PtrCoefFct 
CoefFunction::Generate( MathParser * mp,
                        Global::ComplexPart format,
                        UInt numRows, UInt numCols,
                        const StdVector<std::string>& realVal,
                        const StdVector<std::string>& imagVal ) {
   unsigned int handle = mp->GetNewHandle(true);

   PtrCoefFct ret;

   bool depTime = ExprDependsOnTimeFreq(mp, realVal);
   bool depSpace = ExprDependsOnSpace(mp, realVal);

   if( format == Global::REAL) {
     
     // === REAL CASE ===
     if( depSpace ) {
       // --- a) general case: expression 
       shared_ptr<CoefFunctionExpression<Double> > c 
       (new CoefFunctionExpression<Double>(mp));
       c->SetTensor(realVal, numRows, numCols );
       ret = c;
     } else if (depTime) {
       // --- b) time / freq dependency
       shared_ptr<CoefFunctionTimeFreq<Double> > c 
       (new CoefFunctionTimeFreq<Double>(mp));
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

     depTime |= ExprDependsOnTimeFreq(mp, imagVal);
     depSpace |= ExprDependsOnSpace(mp, imagVal);

     if( depSpace ) {
       // --- a) general case: expression
       shared_ptr<CoefFunctionExpression<Complex> > c 
       (new CoefFunctionExpression<Complex>(mp));
       c->SetTensor(realVal, imagVal, numRows, numCols );
       ret = c;
     } else if (depTime) {
       // --- b) time / freq dependency
       shared_ptr<CoefFunctionTimeFreq<Complex> > c 
       (new CoefFunctionTimeFreq<Complex>(mp));
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

//! Generate tensor-valued coefficient function from CoefFunctions
PtrCoefFct 
CoefFunction::Generate( MathParser * mp,
                        Global::ComplexPart format,
                        UInt numRows, UInt numCols,
                        const StdVector<PtrCoefFct>& realVal,
                        const StdVector<PtrCoefFct>& imagVal )
{
  PtrCoefFct ret;
  UInt realSize = realVal.GetSize();
  UInt imagSize = imagVal.GetSize();

  // Check, if all entries are analytical
  bool isAnalytical = false;
  for ( UInt i = 0; i < realSize; ++i ) {
    if (!realVal[i]) continue;
    isAnalytical |= realVal[i]->IsAnalytic();
    // also ensure, that every CoefFunction is of scalar type
    if( realVal[i]->GetDimType() != CoefFunction::SCALAR )  {
      EXCEPTION( "Vector valued expression can be only composed of "
          << "scalar entries!" );
    }
  }
  for ( UInt i = 0; i < imagSize; ++i ) {
    if (!imagVal[i]) continue;
    isAnalytical |= imagVal[i]->IsAnalytic();
    // also ensure, that every CoefFunction is of scalar type
    if ( imagVal[i]->GetDimType() != CoefFunction::SCALAR )  {
      EXCEPTION( "Vector valued expression can be only composed of "
          << "scalar entries!" );
    }
  }

  if ( isAnalytical ) {
    // use string generation method to generate the vector valued CoefFct
    StdVector<std::string> rStrVals(realSize), iStrVals(imagSize);
    std::string empty;
    shared_ptr<CoefFunctionAnalytic> tmp;

    for ( UInt i = 0; i < realSize; ++i ) {
      if (realVal[i]) {
        tmp = dynamic_pointer_cast<CoefFunctionAnalytic>(realVal[i]);
        tmp->GetStrScalar(rStrVals[i], empty);
      }
      else {
        rStrVals[i] = "0.0";
      }
    }
    for ( UInt i = 0; i < imagSize; ++i ) {
      if (imagVal[i]) {
        tmp = dynamic_pointer_cast<CoefFunctionAnalytic>(imagVal[i]);
        tmp->GetStrScalar(iStrVals[i], empty);
      }
      else {
        iStrVals[i] = "0.0";
      }
    }

    if (imagSize == 0 && format == Global::COMPLEX) {
      iStrVals.Resize(realSize, "0.0");
    }

    ret = Generate(mp, format, numRows, numCols, rStrVals, iStrVals);
  }
  else {
    EXCEPTION( "Vector-valued expression can only be composed of "
        << "analytical expressions currently!" );
  }

  return ret;
}

//! Generate tensor-valued coefficient function from CoefFunctions
PtrCoefFct
CoefFunction::Generate( MathParser * mp,
                        Global::ComplexPart format,
                        UInt numRows, UInt numCols,
                        const StdVector<PtrCoefFct>& scalars)
{
  PtrCoefFct ret;
  UInt vecSize = scalars.GetSize();

  // Check, if all entries are analytical
  bool isAnalytical = false;
  for ( UInt i = 0; i < vecSize; ++i ) {
    if (!scalars[i]) continue;
    isAnalytical |= scalars[i]->IsAnalytic();
    // also ensure, that every CoefFunction is of scalar type
    if( scalars[i]->GetDimType() != CoefFunction::SCALAR )  {
      EXCEPTION( "Vector valued expression can be only composed of "
          << "scalar entries!" );
    }
  }

  if ( isAnalytical ) {
    // use string generation method to generate the vector valued CoefFct
    StdVector<std::string> rStrVals(vecSize), iStrVals(vecSize);
    shared_ptr<CoefFunctionAnalytic> tmp;

    for ( UInt i = 0; i < vecSize; ++i ) {
      if (scalars[i]) {
        tmp = dynamic_pointer_cast<CoefFunctionAnalytic>(scalars[i]);
        tmp->GetStrScalar(rStrVals[i], iStrVals[i]);
      }
      else {
        rStrVals[i] = "0.0";
        iStrVals[i] = "0.0";
      }
    }

    ret = Generate(mp, format, numRows, numCols, rStrVals, iStrVals);
  }
  else {
    EXCEPTION( "Vector-valued expression can only be composed of "
        << "analytical expressions currently!" );
  }

  return ret;
}

PtrCoefFct CoefFunction::Generate( MathParser * mp,
                                   Global::ComplexPart type,
                                   const CoefXpr& xpr ) {
  
  PtrCoefFct ret;
  // check if expression is analytical
  if( xpr.IsAnalytical( ) ) {
    
    // ==== ANALYTICAL COEFFICIENT FUNCTION ===
    // Use normal factory method to generate real or complex valued function
    Global::ComplexPart part = type;
    if( xpr.GetDimType() == CoefFunction::SCALAR ) {
      std::string real, imag;
      xpr.GetScalarXpr( real, imag );
      ret = Generate( mp, part, real, imag );
    } else if (xpr.GetDimType() == CoefFunction::VECTOR ) {
      StdVector<std::string> real, imag;
      xpr.GetVectorXpr( real, imag );
      ret = Generate( mp, part, real, imag );
    } else if (xpr.GetDimType() == CoefFunction::TENSOR ) {
      StdVector<std::string> real, imag;
      UInt numRows, numCols;
      xpr.GetTensorXpr( numRows, numCols, real, imag );
      ret = Generate( mp, part, numRows, numCols, real, imag );
    } else {
        EXCEPTION( "Unknown dimtype of coefficient function" );
    }
  } else {
    
    // ==== GENERAL COEFFICIENT FUNCTION ===
    // In this case we have to generate a compound coefficient function, as the
    // single arguments involved can not be represented in closed expression form.
    Global::ComplexPart part = type;
    std::map<std::string, PtrCoefFct > vars;
    xpr.GetArgs( vars );
    if( part == Global::REAL) {
      // ===== REAL-VALUED PART ====
      shared_ptr<CoefFunctionCompound<Double> > cf (new CoefFunctionCompound<Double>(mp));
      if( xpr.GetDimType() == CoefFunction::SCALAR ) {
        std::string real, imag;
        xpr.GetScalarXpr( real, imag );
        cf->SetScalar( real, vars );

      } else if (xpr.GetDimType() == CoefFunction::VECTOR ) {
        StdVector<std::string> real, imag;
        xpr.GetVectorXpr( real, imag );
        cf->SetVector( real, vars );

      } else if (xpr.GetDimType() == CoefFunction::TENSOR ) {
        UInt numRows, numCols;
        StdVector<std::string> real, imag;
        xpr.GetTensorXpr( numRows, numCols, real, imag );
        cf->SetTensor( real, numRows, numCols, vars );
      } else {
        EXCEPTION( "Unknown dimtype of coefficient function" );
      }
      ret = cf;
    } else {
      // ===== COMPLEX VALUED PART ====
      shared_ptr<CoefFunctionCompound<Complex> > cf (new CoefFunctionCompound<Complex>(mp));
      if( xpr.GetDimType() == CoefFunction::SCALAR ) {
        std::string real, imag;
        xpr.GetScalarXpr( real, imag );
        cf->SetScalar( real, imag, vars );

      } else if (xpr.GetDimType() == CoefFunction::VECTOR ) {
        StdVector<std::string> real, imag;
        xpr.GetVectorXpr( real, imag );
        cf->SetVector( real, imag, vars );

      } else if (xpr.GetDimType() == CoefFunction::TENSOR ) {
        UInt numRows, numCols;
        StdVector<std::string> real, imag;
        xpr.GetTensorXpr( numRows, numCols, real, imag );
        cf->SetTensor( real, imag, numRows, numCols, vars );
      } else {
        EXCEPTION( "Unknown dimtype of coefficient function" );
      }
      ret = cf;
    }
    // Determine the dependType of the returned coefFunction
    ret->dependType_ = xpr.GetDependency(); 
  }
  
  // Check for coordinate system
  if( xpr.GetCoordinateSystem() ) {
    ret->SetCoordinateSystem( xpr.GetCoordinateSystem() );
  }
  return ret;
}



bool CoefFunction::ExprDependsOnTimeFreq(MathParser* mp, const std::string& expr) {
  unsigned int handle = mp->GetNewHandle(true);
  mp->SetExpr(handle, expr);
  bool depends = false;
  if ( mp->IsExprVariable(handle, "t") ||
       mp->IsExprVariable(handle, "f") ){
    depends = true;
  }
  
  mp->ReleaseHandle(handle);
  return depends;
}

bool CoefFunction::ExprDependsOnTimeFreq(MathParser* mp, const StdVector<std::string>& expr) {
  bool dep = false;
  for( UInt i = 0; i < expr.GetSize(); ++i ) {
    dep  |= ExprDependsOnTimeFreq(mp, expr[i]);
  }
  return dep;
}

bool CoefFunction::ExprDependsOnSpace(MathParser* mp, const std::string& expr) {
  unsigned int handle = mp->GetNewHandle(true);
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

bool CoefFunction::ExprDependsOnSpace(MathParser* mp, const StdVector<std::string>& expr) {
  bool dep = false;
  for( UInt i = 0; i < expr.GetSize(); ++i ) {
    dep  |= ExprDependsOnSpace(mp, expr[i]);
  }
  return dep;
}

void  CoefFunction::GenScalCompNames( std::string& realVar, 
                                      std::string& imagVar, 
                                      const std::string& prefix,
                                      PtrCoefFct cf ) {
  assert( cf->GetDimType() == SCALAR );
  realVar = prefix + "_R";
  //if( cf->IsComplex() )
    imagVar = prefix + "_I";
}

void CoefFunction::GenVecCompNames( StdVector<std::string>& realVar, 
                                    StdVector<std::string>& imagVar, 
                                    const std::string& prefix,
                                    PtrCoefFct cf ) {
  
  assert( cf->GetDimType() == VECTOR );
  UInt numEntries = cf->GetVecSize();
  realVar.Resize( numEntries );
  for( UInt i = 0; i < numEntries; ++i ) {
    realVar[i] = prefix + "_" + lexical_cast<std::string>(i) + "_R";
  }
  if( cf->IsComplex() ) {
    imagVar.Resize( numEntries );
    for( UInt i = 0; i < numEntries; ++i ) {
      imagVar[i] = prefix + "_" + lexical_cast<std::string>(i) + "_I";
    }
  }
}

void CoefFunction::GenTensorCompNames( StdVector<std::string>& realVar,
                                       StdVector<std::string>& imagVar,
                                       const std::string& prefix,
                                       PtrCoefFct cf ) {
  assert( cf->GetDimType() == TENSOR );
  UInt numRows, numCols;
  cf->GetTensorSize( numRows, numCols );
  UInt numEntries = numRows * numCols;
  UInt pos = 0;
  realVar.Resize( numEntries );
  for( UInt i = 0; i < numRows; ++i ) {
    for( UInt j = 0; j < numCols; ++j ) {
      realVar[pos++] = prefix + "_" + lexical_cast<std::string>(i)  
                       + lexical_cast<std::string>(j) + "_R";
    }
  }
  if( cf->IsComplex() ) {
    imagVar.Resize( numEntries );
    pos = 0;
    for( UInt i = 0; i < numRows; ++i ) {
      for( UInt j = 0; j < numCols; ++j ) {
        imagVar[pos++] = prefix + "_" + lexical_cast<std::string>(i)  
                         + lexical_cast<std::string>(j) + "_I";
      }
    }
  }

}

void CoefFunction::GetVectorValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                            StdVector< Vector<Double> >& values, 
                                            Grid* ptGrid,
                                            const StdVector<shared_ptr<EntityList> >& srcEntities,
                                            bool updatedGeo) {

  // loop over all coordinates
  UInt numEntries = globCoord.GetSize();
  values.Resize( numEntries );
  StdVector< LocPoint > localCoords;
  StdVector< const Elem* >  elems;

  // Get mapping
  ptGrid->GetElemsAtGlobalCoords( globCoord, localCoords, elems, srcEntities, 1e-3, 1e-2, true, updatedGeo );


  LocPointMapped lpm;
  for( UInt i = 0; i < numEntries; ++i ) {
    if (elems[i] == NULL ) {
      WARN("Could not find suitable element for position " << globCoord[i]  
           << ". Setting values to zero for this point" );
    } else {
      LocPoint & lp = localCoords[i];
      const Elem * ptElem = elems[i];
      // shared_ptr<ElemShapeMap>  esm = ptGrid->GetElemShapeMap( ptElem );
      // shared_ptr<ElemShapeMap> esm(ptElem->ptrShapeMap);
      shared_ptr<ElemShapeMap> esm = (ptElem)->GetElemShapeMap(ptGrid);
      lpm.Set( lp, esm, 0);
      this->GetVector( values[i], lpm);
    }
  }
}
void CoefFunction::GetVectorValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                            StdVector< Vector<Complex> >& values, 
                                            Grid* ptGrid,
                                            const StdVector<shared_ptr<EntityList> >& srcEntities,
                                            bool updatedGeo) {
  // loop over all coordinates
  UInt numEntries = globCoord.GetSize();
  values.Resize( numEntries );
  StdVector< LocPoint > localCoords;
  StdVector< const Elem* >  elems;

  // Get mapping
  ptGrid->GetElemsAtGlobalCoords( globCoord, localCoords, elems, srcEntities, 1e-3, 1e-2, true, updatedGeo );

  LocPointMapped lpm;
  for( UInt i = 0; i < numEntries; ++i ) {
    LocPoint & lp = localCoords[i];
    // Note: we can not be sure that an element was always found
    if (elems[i] == NULL ) {
      WARN("Could not find suitable element for position " << globCoord[i]  
           << ". Setting values to zero for this point" );
    } else {
      const Elem * ptElem = elems[i];
      // shared_ptr<ElemShapeMap>  esm = ptGrid->GetElemShapeMap( ptElem );
      // shared_ptr<ElemShapeMap>  esm(ptElem->ptrShapeMap);
      shared_ptr<ElemShapeMap> esm = (ptElem)->GetElemShapeMap(ptGrid);
      lpm.Set( lp, esm, 0);
      this->GetVector( values[i], lpm);
    }
  }
}

void CoefFunction::GetScalarValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                            StdVector< Double >& values,
                                            Grid* ptGrid,
                                            const StdVector<shared_ptr<EntityList> >& srcEntities) {
  // loop over all coordinates
  UInt numEntries = globCoord.GetSize();
  values.Resize( numEntries );
  StdVector< LocPoint > localCoords;
  StdVector< const Elem* >  elems;

  // Get mapping
  ptGrid->GetElemsAtGlobalCoords( globCoord, localCoords, elems, srcEntities );

  LocPointMapped lpm;
  for( UInt i = 0; i < numEntries; ++i ) {
    if (elems[i] == NULL ) {
      WARN("Could not find suitable element for position " << globCoord[i]  
           << ". Setting values to zero for this point" );
    } else {
      LocPoint & lp = localCoords[i];
      const Elem * ptElem = elems[i];
      // shared_ptr<ElemShapeMap>  esm = ptGrid->GetElemShapeMap( ptElem );
      // shared_ptr<ElemShapeMap> esm(ptElem->ptrShapeMap);
      shared_ptr<ElemShapeMap> esm = (ptElem)->GetElemShapeMap(ptGrid);
      lpm.Set( lp, esm, 0);
      this->GetScalar( values[i], lpm);
    }
  }
}
void CoefFunction::GetScalarValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                            StdVector< Complex >& values,
                                            Grid* ptGrid,
                                            const StdVector<shared_ptr<EntityList> >& srcEntities) {
  // loop over all coordinates
  UInt numEntries = globCoord.GetSize();
  values.Resize( numEntries );
  StdVector< LocPoint > localCoords;
  StdVector< const Elem* >  elems;

  // Get mapping
  ptGrid->GetElemsAtGlobalCoords( globCoord, localCoords, elems, srcEntities );

  LocPointMapped lpm;
  for( UInt i = 0; i < numEntries; ++i ) {
    if (elems[i] == NULL ) {
      WARN("Could not find suitable element for position " << globCoord[i]  
           << ". Setting values to zero for this point" );
    } else {
      LocPoint & lp = localCoords[i];
      const Elem * ptElem = elems[i];
      // shared_ptr<ElemShapeMap>  esm = ptGrid->GetElemShapeMap( ptElem );
      // shared_ptr<ElemShapeMap> esm(ptElem->ptrShapeMap);
      shared_ptr<ElemShapeMap> esm = (ptElem)->GetElemShapeMap(ptGrid);
      lpm.Set( lp, esm, 0);
      this->GetScalar( values[i], lpm);
    }
  }
}

template<class TYPE>
CoefFunctionConst<TYPE>* CoefFunction::AsConst(bool throw_exception)
{
  CoefFunctionConst<TYPE>* ret = dynamic_cast<CoefFunctionConst<TYPE>* >(this);

  if(ret == NULL && throw_exception)
    EXCEPTION("cannot cast coef" << ToString() << "to CoefFunctionConst.")

  return ret;
}

//! Rotates a vector from the local to the global coordinate system
template<typename TYPE>
void CoefFunction::TransformVectorByCoordSys(Vector<TYPE> &outVec,
                                             const Vector<TYPE> &inVec,
                                             const Vector<Double> &point)
{
  if (coordSys_ && coordSys_->GetName() != "default") {
    // rotate the local vector back to global coordinates
    this->coordSys_->Local2GlobalVector(outVec, inVec, point);
  } else {
    outVec = inVec;
  }
}

//! Rotates a Vector from the local to the global coordinate system
template<typename TYPE>
void CoefFunction::TransformVectorByCoordSys(Vector<TYPE> &outVec,
                                             const Vector<TYPE> &inVec,
                                             const LocPointMapped &lpm)
{
  if (coordSys_ && coordSys_->GetName() != "default") {
    // rotate the local vector back to global coordinates
    Vector<Double> point;
    lpm.shapeMap->Local2Global(point, lpm.lp);
    this->coordSys_->Local2GlobalVector(outVec, inVec, point);
  } else {
    outVec = inVec;
  }
}

// Rotates a tensor from the local to the global coordinate system
template<typename TYPE>
void CoefFunction::TransformTensorByCoordSys(Matrix<TYPE> &outMat,
                                          const Matrix<TYPE> &inMat,
                                          const Vector<Double> &point)
{
  // Just return the original tensor if there is no coordinate system
  if (!coordSys_ || coordSys_->GetName() == "default") {
    outMat = inMat;
    return;
  }

  Matrix<Double> rotMat;
  coordSys_->GetFullGlobRotationMatrix(rotMat, point);
  inMat.PerformRotation(rotMat, outMat);
}

// Rotates a tensor from the local to the global coordinate system
template<typename TYPE>
void CoefFunction::TransformTensorByCoordSys(Matrix<TYPE> &outMat,
                                          const Matrix<TYPE> &inMat,
                                          const LocPointMapped &lpm)
{
  // Just return the original tensor if there is no coordinate system
  if (!coordSys_ || coordSys_->GetName() == "default") {
    outMat = inMat;
    return;
  }

  Matrix<Double> rotMat;
  coordSys_->GetFullGlobRotationMatrix(rotMat, lpm);
  inMat.PerformRotation(rotMat, outMat);
}

// ************************************************************************
// EXPLICIT TEMPLATE INSTANTIATION
// ************************************************************************
template void CoefFunction::TransformVectorByCoordSys(Vector<Double> &outVec,
                                                      const Vector<Double> &inVec,
                                                      const Vector<Double> &point);
template void CoefFunction::TransformVectorByCoordSys(Vector<Complex> &outVec,
                                                      const Vector<Complex> &inVec,
                                                      const Vector<Double> &point);
template void CoefFunction::TransformVectorByCoordSys(Vector<Double> &outVec,
                                                      const Vector<Double> &inVec,
                                                      const LocPointMapped &lpm);
template void CoefFunction::TransformVectorByCoordSys(Vector<Complex> &outVec,
                                                      const Vector<Complex> &inVec,
                                                      const LocPointMapped &lpm);
template void CoefFunction::TransformTensorByCoordSys(Matrix<Double> &outMat,
                                                      const Matrix<Double> &inMat,
                                                      const Vector<Double> &point);
template void CoefFunction::TransformTensorByCoordSys(Matrix<Complex> &outMat,
                                                      const Matrix<Complex> &inMat,
                                                      const Vector<Double> &point);
template void CoefFunction::TransformTensorByCoordSys(Matrix<Double> &outMat,
                                                      const Matrix<Double> &inMat,
                                                      const LocPointMapped &lpm);
template void CoefFunction::TransformTensorByCoordSys(Matrix<Complex> &outMat,
                                                      const Matrix<Complex> &inMat,
                                                      const LocPointMapped &lpm);

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

Enum<CoefFunction::CoefDimType> CoefFunction::coefDimType = \
    Enum<CoefFunction::CoefDimType>("Dimension of CoefFunction",
                                    sizeof(dimTypeTuples) / sizeof(EnumTuple),
                                    dimTypeTuples);

// Definition of coefficient function dependency type
static EnumTuple coefDependTuples[] = 
{
 EnumTuple(CoefFunction::CONSTANT, "CONSTANT"), 
 EnumTuple(CoefFunction::TIMEFREQ, "TIMEFREQ"),
 EnumTuple(CoefFunction::SPACE, "SPACE"),
 EnumTuple(CoefFunction::GENERAL,  "GENERAL"),
 EnumTuple(CoefFunction::SOLUTION, "SOLUTION")
};

Enum<CoefFunction::CoefDependType>CoefFunction::coefDependType = \
    Enum<CoefFunction::CoefDependType>("Dependency of CoefFunction",
                                       sizeof(coefDependTuples) / sizeof(EnumTuple),
                                       coefDependTuples);



#ifdef PERFORM_COMPOUND_TEST
void CoefCompoundTest() {
  {
//    PtrCoefFct a =
//        CoefFunction::Generate(Global::COMPLEX, "5.0", "2.0" );
//    //      CoefFunction::Generate(Global::REAL, "5.0");
//    PtrCoefFct b =
//        CoefFunction::Generate(Global::COMPLEX, "3.0", "4.0" );
//
//    CoefXprBinOp x(a, b, CoefXpr::OP_DIV );
//
//    PtrCoefFct ret = 
//        CoefFunction::Generate(Global::COMPLEX, CoefXprBinOp (a, b, CoefXpr::OP_DIV ));
//
//    std::string real, imag;
//    x.GetScalarXpr( real, imag );
//    std::cerr << a->ToString() << " * " << b->ToString()
//           << " = (" << real << ", " << imag << ")\n";
//    std::cerr << "It evaluates to " << ret->ToString();
//
//
//
//
//    PtrCoefFct ret2;
//
//    //PtrCoefFct temp;
//    //temp = CoefFunction::Generate(Global::COMPLEX, CoefXprBinOp(a,a,CoefXpr::OP_MULT) );
//
//
//
//    //std::cerr << "a-squared is " << temp->ToString();
//
//    //oefXprBinOp xpr ("1.0", temp , CoefXpr::OP_DIV );
//    // calculate sqrt( abs( 1/a^2) )
//    ret2 = CoefFunction::Generate(
//        Global::REAL,
//        CoefXprUnaryOp(
//            CoefXprUnaryOp(
//                CoefXprBinOp("1.0", 
//                             CoefXprBinOp(a,a,CoefXpr::OP_MULT),
//                             CoefXpr::OP_DIV),
//                             CoefXpr::OP_NORM ),
//                             CoefXpr::OP_SQRT ) );
//    std::cerr << "Wave velocity evaluates to " << ret2->ToString() << std::endl;
//    std::cerr << "\n\n";
  }
  
  // ========================================================================
  //   VECTOR TEST
  // ========================================================================
  {
    StdVector<std::string> aValsReal, aValsImag, bValsReal, bValsImag;
    aValsReal = "1", "1.0", "2.5" ;
    aValsImag = "0.0", "1.0", "4.0";
    
    bValsReal = "8", "6.0", "3.5";
    bValsImag = "4.0", "0.0", "1.0";
    Domain  * myDom = domain;
    
    
    PtrCoefFct a = 
        CoefFunction::Generate(Global::COMPLEX, aValsReal, aValsImag );
    PtrCoefFct b = 
            CoefFunction::Generate(Global::COMPLEX, bValsReal, bValsImag );
    
    CoefXprBinOp xpr(a,b,CoefXpr::OP_CROSS); 
    
    PtrCoefFct res = 
        CoefFunction::Generate(Global::COMPLEX, xpr );
    
    std::cerr << "result is " << res->ToString() << std::endl;

    // generate dummy local point mapped for first element to see the 
    // toolchain working
    const Elem * ptEl = myDom->GetGrid()->GetElem( 1 );
    // shared_ptr<ElemShapeMap> esm = 
    //     myDom->GetGrid()->GetElemShapeMap( ptEl, false );
    // shared_ptr<ElemShapeMap> esm(ptEl->ptrShapeMap);
    shared_ptr<ElemShapeMap> esm = (ptEl)->GetElemShapeMap(myDom->GetGrid());
    LocPoint lp;
    Vector<Double> midPoint(3);
    midPoint.Init();
    lp.coord = midPoint;
    LocPointMapped lpm;
    lpm.Set( lp, esm);

    //    Double scal;
    //    res->GetScalar( scal, lpm );
    //    std::cerr << "scalar double result is " << scal << std::endl;
    Vector<Complex> out;
    res->GetVector( out, lpm );
    std::cerr << "vector valued result is " << out << std::endl;
    
  }
  // ========================================================================
   //   TENSOR TEST
   // ========================================================================
   { 
     
     std::cerr << "=============\n";
     std::cerr << " TENSOR TEST \n";
     std::cerr << "=============\n";
    
     StdVector<std::string> aValsReal, aValsImag, bValsReal, bValsImag;
     aValsReal = "1", "1.0", "2.5",   "8.5";
     aValsImag = "0.0", "1.0", "4.0", "-2.0";
     
     //aValsReal = "1", "1.0";
     //aValsImag = "0.0", "1.0";

     bValsReal = "8",    "6.0", "-3.5", "3.0";
     bValsImag = "-4.0", "0.0",  "1.0", "7.0";
     PtrCoefFct a = CoefFunction::Generate(Global::COMPLEX, 2,2, aValsReal, aValsImag);
     PtrCoefFct b = CoefFunction::Generate(Global::COMPLEX, 2,2, bValsReal, bValsImag);
     
     CoefXprBinOp xpr(a,b,CoefXpr::OP_MULT);
     PtrCoefFct res = 
         CoefFunction::Generate(Global::COMPLEX, xpr );

     std::cerr << "result is " << res->ToString() << std::endl;
     // generate dummy local point mapped for first element to see the 
     // toolchain working
     const Elem * ptEl = myDom->GetGrid()->GetElem( 1 );
    //  shared_ptr<ElemShapeMap> esm = 
    //      myDom->GetGrid()->GetElemShapeMap( ptEl, false );
    // shared_ptr<ElemShapeMap> esm(ptEl->ptrShapeMap);
    shared_ptr<ElemShapeMap> esm = (ptEl)->GetElemShapeMap(myDom->GetGrid());
     LocPoint lp;
     Vector<Double> midPoint(3);
     midPoint.Init();
     lp.coord = midPoint;
     LocPointMapped lpm;
     lpm.Set( lp, esm);

     //    Double scal;
     //    res->GetScalar( scal, lpm );
     //    std::cerr << "scalar double result is " << scal << std::endl;
     Matrix<Complex> out;
     res->GetTensor( out, lpm );
     std::cerr << "tensor valued result is " << out << std::endl;

   }
//  // ========================================================================
//  //   MATHPARSER EXTERNAL VARIABLE TEST
//  // ========================================================================
//  {
//    MathParser * mp = myDom->GetMathParser();
//    unsigned int h = mp->GetNewHandle();
//    
//    Double var = 42.0;
//    mp->RegisterExternalVar( h, "var", &var );
//    mp->SetExpr(h, "var *2");
//    std::cerr << mp->Eval(h) << std::endl;
//    
//    
//    
//  }
//  
  
//  // ========================================================================
//  //   COMPOUND TEST
//  // ========================================================================
//  {
//    PtrCoefFct a = 
//              CoefFunction::Generate(Global::COMPLEX, "2.0", "1.0");
//    PtrCoefFct b = 
//           CoefFunction::Generate(Global::COMPLEX, "3.0", "4.0");
//    StdVector<std::string> cValsReal, cValsImag;
//    cValsReal = "5", "1.0", "2.5";
//    cValsImag = "2", "4.0", "1.5";
//    PtrCoefFct c = 
//     CoefFunction::Generate(Global::COMPLEX, cValsReal, cValsImag );
//    
//    // CoefXprBinOp xpr (a,b, CoefXpr::OP_MULT);
//    CoefXprVecScalOp xpr(c,a, CoefXpr::OP_MULT );
//    
//    PtrCoefFct res = CoefFunction::Generate( Global::COMPLEX, xpr );
//    std::cerr << "res is " << res->ToString();
//    
//    const Elem * ptEl = myDom->GetGrid()->GetElem( 1 );
//    shared_ptr<ElemShapeMap> esm = 
//            myDom->GetGrid()->GetElemShapeMap( ptEl, false );
//
//    // generate dummy local point mapped for first element to see the 
//    // toolchain working
//    LocPoint lp;
//    Vector<Double> midPoint(3);
//    midPoint.Init();
//    lp.coord = midPoint;
//    LocPointMapped lpm;
//    lpm.Set( lp, esm);
//    
////    Double scal;
////    res->GetScalar( scal, lpm );
////    std::cerr << "scalar double result is " << scal << std::endl;
//    Vector<Complex> out;
//    res->GetVector( out, lpm );
//    std::cerr << "vector valued result is " << out << std::endl;
//    
//    
    
//  }
  
}
#endif

// template instantiation stuff
template CoefFunctionConst<double>* CoefFunction::AsConst<double>(bool throw_exception);
template CoefFunctionConst<std::complex<double> >* CoefFunction::AsConst<std::complex<double> >(bool throw_exception);

} // end of namespace
