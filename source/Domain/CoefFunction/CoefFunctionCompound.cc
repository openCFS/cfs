#include "CoefFunctionCompound.hh"

#include "Domain/CoordinateSystems/CoordSystem.hh"
#include "Utils/mathParser/mathParser.hh"

namespace CoupledField {

// ===========================================================================
//  REAL VALUED COEFFICIENT FUNCTION
// ===========================================================================

CoefFunctionCompound<Double>::CoefFunctionCompound(MathParser * mp) {
  dependType_ = CoefFunction::GENERAL;
  isAnalytic_ = false;
  isComplex_ = false;
  numRows_ = 0;
  numCols_ = 0;
  parser_ = NULL;

  parser_ = mp;
  handle_ = parser_->GetNewHandle( true );

  // always store default coordinate system
  coordSysDefault_ = domain->GetCoordSystem();
}

CoefFunctionCompound<Double>::~CoefFunctionCompound() {
  parser_->ReleaseHandle( handle_ );
}

void CoefFunctionCompound<Double>::
GetTensor( Matrix<Double>& coefMat, const LocPointMapped& lpm ) {
  assert(this->dimType_ == CoefFunction::TENSOR);
  Vector<Double> pointCoord;
  Matrix<Double> locMatrix;

  // First, obtain global coordinates of current point, register it at the mathParser
  // and compute local coefficient vector
  lpm.shapeMap->Local2Global(pointCoord,lpm.lp);
  parser_->SetCoordinates( handle_, *(this->coordSysDefault_), pointCoord);

  // update internal variables
  UpdateXpr( lpm );

  parser_->EvalMatrix( handle_, locMatrix,
                         this->numRows_, this->numCols_ );

  // Rotate material, if coordinate system is not the global one
  TransformTensorByCoordSys(coefMat, locMatrix, lpm);
}

void CoefFunctionCompound<Double>::
GetVector( Vector<Double>& coefVec, const LocPointMapped& lpm ) {
  assert(this->dimType_ == VECTOR );

  Vector<Double> pointCoord;

  // First, obtain global coordinates of current point,
  lpm.GetGlobal(pointCoord);

  // register it at the mathParser and compute local coefficient vector
  parser_->SetCoordinates( handle_, *(this->coordSysDefault_), pointCoord);

  // update internal variables
  UpdateXpr(lpm);

  Vector<Double> locVector;
  parser_->EvalVector( handle_, locVector );

  if( this->coordSys_ ) {
    // Afterwards rotate the local vector back to global coordinates
    this->coordSys_->Local2GlobalVector( coefVec, locVector, pointCoord );
  } else {
    coefVec = locVector;
  }
}

void CoefFunctionCompound<Double>::
GetScalar( Double& coefScalar, const LocPointMapped& lpm ) {
  assert(this->dimType_ == SCALAR);

  // First, obtain global coordinates of current point and  register it at the mathParser
  Vector<Double> pointCoord;
  assert(lpm.shapeMap != nullptr);
  lpm.shapeMap->Local2Global(pointCoord,lpm.lp);
  parser_->SetCoordinates( handle_, *(this->coordSysDefault_), pointCoord);

  // update internal variables
  UpdateXpr( lpm );
  coefScalar = parser_->Eval( handle_ );
}

void CoefFunctionCompound<Double>::
SetScalar( const std::string& expr,
           std::map<std::string, PtrCoefFct >& vars ) {
  dimType_ = SCALAR;
  
  // loop over all coefficient functions and register them
  std::map<std::string, PtrCoefFct >::iterator it = vars.begin();
  for( ; it != vars.end(); ++it ) {
    RegisterCoefFct( it->first, it->second );
  }
  
  // set expression with mathparser
  expr_ = expr;
  parser_->SetExpr( handle_, expr );
}

void CoefFunctionCompound<Double>::
SetVector( StdVector<std::string>& expr, 
           std::map<std::string, PtrCoefFct >& vars ) {
  dimType_ = VECTOR;

  // loop over all coefficient functions and register them
  std::map<std::string, PtrCoefFct >::iterator it = vars.begin();
  for( ; it != vars.end(); ++it ) {
    RegisterCoefFct( it->first, it->second );
  }

  // set expression with mathparser
  expr_ = expr.Serialize(',');
  parser_->SetExpr( handle_, expr_ );
}

void CoefFunctionCompound<Double>::
SetTensor( StdVector<std::string>& expr,
           UInt numRows, UInt numCols,
           std::map<std::string, PtrCoefFct >& vars ) {
  dimType_ = TENSOR;

  // security check: make sure that nRows x nCols = val.size()
  if( numCols * numRows != expr.GetSize() ) {
    EXCEPTION("Tensor is supposed to have size " << numRows << " x " 
              << numCols << ", but only " << expr.GetSize() 
              << " entries have been provided!");
  }
  this->numRows_ = numRows;
  this->numCols_ = numCols;

  // loop over all coefficient functions and register them
  std::map<std::string, PtrCoefFct >::iterator it = vars.begin();
  for( ; it != vars.end(); ++it ) {
    RegisterCoefFct( it->first, it->second );
  }

  // set expression with mathparser
  expr_ = expr.Serialize(',');
  parser_->SetExpr( handle_, expr_ );
}

UInt CoefFunctionCompound<Double>::GetVecSize() const {
  assert(this->dimType_ == VECTOR );
  return parser_->GetNumExprs(handle_);
}

void CoefFunctionCompound<Double>::GetTensorSize( UInt& numRows, UInt& numCols ) const {
  assert(this->dimType_ == TENSOR );
  numRows = numRows_;
  numCols = numCols_;
}

std::string CoefFunctionCompound<Double>::ToString() const {
  std::stringstream out;
  out << "CoefFunctionCompound" << std::endl;
  out << "expression: '" << expr_ << std::endl;
  out << "registered variables:" << std::endl;
  std::map<std::string, PtrCoefFct>::const_iterator it;
  it = coefs_.begin();
  for( ; it != coefs_.end(); ++it ) {
    out << "\tvariable: " << it->first << std::endl;
    out << "\ttype: " << coefDimType.ToString(it->second->GetDimType() ) << std::endl;
    out << "\tvalue: " << it->second->ToString() << std::endl;
  }
  
  return out.str();
}

void CoefFunctionCompound<Double>::
RegisterCoefFct( const std::string& name,
                 PtrCoefFct& coef ) {
#ifdef USE_OPENMP
  if(omp_get_num_threads() > 1 ){
    EXCEPTION("Calling from parallel region. Not allowed.")
  }
#endif

  // make sure, that coefficient function with the same name is not already
  // registered
  if( coefs_.find( name ) != coefs_.end() ) {
    EXCEPTION( "A coefficient function with name '" << name 
               << "' was already registered!" );
  }
  
  // adjust dependency type of own coefficient function
  dependType_ = std::max(this->GetDependency(), 
                         coef->GetDependency());
  
  // query type of coefficient
  CoefDimType dim = coef->GetDimType();
  coefs_[name] = coef;
  coefDimTypes_[name] = dim;
  
  if( dim == SCALAR ) {
    
    // query variables names
    std::string real, imag;
    CoefFunction::GenScalCompNames( real, imag, name, coef );
#pragma omp parallel num_threads(CFS_NUM_THREADS)
    {
#pragma omp critical (CoefFunctionCompound_Double)
      {
        // insert value
        parser_->RegisterExternalVar(handle_, real,  &(scalVars_[name]) );
      }
    }
    
  } else if( dim == VECTOR ) {
  
    StdVector<std::string> real, imag;
    CoefFunction::GenVecCompNames( real, imag, name, coef );
#pragma omp parallel num_threads(CFS_NUM_THREADS)
    {
#pragma omp critical (CoefFunctionCompound_Double)
      {
        vecVars_[name].Resize(real.GetSize());
        for( UInt i = 0; i < real.GetSize(); ++ i ) {
          parser_->RegisterExternalVar(handle_, real[i],  &(vecVars_[name][i]) );
        }
      }
    }
    
  } else if( dim == TENSOR ) {
    StdVector<std::string> real, imag;
    CoefFunction::GenTensorCompNames( real, imag, name, coef );
#pragma omp parallel num_threads(CFS_NUM_THREADS)
    {
#pragma omp critical (CoefFunctionCompound_Double)
      {
        tensorVars_[name].Resize(numRows_, numCols_);
        UInt pos = 0;
        for( UInt i = 0; i < numRows_; ++i ) {
          for( UInt j = 0; j < numCols_; ++j ) {
            parser_->RegisterExternalVar(handle_, real[pos],  &(tensorVars_[name][i][j]) );
            pos++;
          }
        }
      }
    }
  } else {
    EXCEPTION( "Unknown dimension type of coefficient function '" << name << "'")
  }
}

void CoefFunctionCompound<Double>::
UpdateXpr( const LocPointMapped& lpm ) {
  //should be fine if parallel or not...

  // loop over all registered coefficients and update their values
  std::map<std::string, PtrCoefFct >::iterator it = coefs_.begin();
  for( ; it != coefs_.end(); ++it ) {
    const std::string & name = it->first;
    // query dimType
    CoefDimType dim = coefDimTypes_[name];
    if( dim == SCALAR ) {
      it->second->GetScalar( scalVars_[name], lpm );
    } else if( dim == VECTOR ) {
      // This is a dirty hack for buckling optimization.
      // In that case, we have two pdes, first linear elasticity (real),
      // then buckling (complex eigenvalue problem). In the optimization
      // (Excitation::Apply) we take the solution of the first pde,
      // convert it to complex and inject it in the second one to get the
      // stresses which we need to assemble the buckling pde. For assembly
      // the values have to be real again.
      if (it->second->IsComplex()) {
        Vector<Complex> vec;
        it->second->GetVector( vec, lpm );
        vecVars_[name] = vec.GetPart(Global::REAL);
      } else {
        it->second->GetVector( vecVars_[name], lpm );
      }
    } else if( dim == TENSOR ) {
      it->second->GetTensor( tensorVars_[name], lpm );
    } else {
      EXCEPTION( "Unknown dimension type of coefficient function '" << name << "'")
    }
  }
}


// ===========================================================================
//  Complex/Double VALUED COEFFICIENT FUNCTION
// ===========================================================================
void CoefFunctionCompound<Double>::CreateDivOperator(UInt spaceDim, UInt dofDim){

  if(spaceDim != dofDim)
    EXCEPTION("CoefFunctionCompound<Double>: Divergence need vectorial data!");

  if(spaceDim == 2)
    this->myOperator_.reset(new ScalarDivergenceOperator<FeH1,2,Double>());
  else if(spaceDim == 3)
    this->myOperator_.reset(new ScalarDivergenceOperator<FeH1,3,Double>());
}

void CoefFunctionCompound<Complex>::CreateDivOperator(UInt spaceDim, UInt dofDim){

  if(spaceDim != dofDim)
    EXCEPTION("CoefFunctionCompound<Complex>: Divergence need vectorial data!");

  if(spaceDim == 2)
    this->myOperator_.reset(new ScalarDivergenceOperator<FeH1,2,Complex>());
  else if(spaceDim == 3)
    this->myOperator_.reset(new ScalarDivergenceOperator<FeH1,3,Complex>());
}
// ===========================================================================
//  Complex VALUED COEFFICIENT FUNCTION
// ===========================================================================

CoefFunctionCompound<Complex>::CoefFunctionCompound(MathParser * mp) {
  dependType_ = CoefFunction::GENERAL;
  isAnalytic_ = false;
  isComplex_ = true;

  numRows_ = 0;
  numCols_ = 0;
  parser_ = NULL;

  parser_ = mp;
  handleReal_ = parser_->GetNewHandle( true );
  handleImag_ = parser_->GetNewHandle( true );

  // always store default coordinate system
  coordSysDefault_ = domain->GetCoordSystem();
}

CoefFunctionCompound<Complex>::~CoefFunctionCompound() {
  parser_->ReleaseHandle( handleReal_ );
  parser_->ReleaseHandle( handleImag_ );
}

void CoefFunctionCompound<Complex>::
GetTensor( Matrix<Complex>& coefMat, const LocPointMapped& lpm ) {
  assert(this->dimType_ == CoefFunction::TENSOR);
  Vector<Double> pointCoord;
  Matrix<Complex> locMatrix (this->numRows_, this->numCols_ );
  Matrix<Double> temp;

  // First, obtain global coordinates of current point, register it at the mathParser
  // and compute local coefficient vector
  lpm.shapeMap->Local2Global(pointCoord,lpm.lp);
  parser_->SetCoordinates( handleReal_, *(this->coordSysDefault_), pointCoord);
  parser_->SetCoordinates( handleImag_, *(this->coordSysDefault_), pointCoord);

  // update internal variables
  UpdateXpr( lpm );

  parser_->EvalMatrix( handleReal_, temp,
                       this->numRows_, this->numCols_ );
  locMatrix.SetPart(Global::REAL, temp);
  parser_->EvalMatrix( handleImag_, temp,
                       this->numRows_, this->numCols_ );
  locMatrix.SetPart(Global::IMAG, temp);

  // Rotate material, if coordinate system is not the global one
  TransformTensorByCoordSys(coefMat, locMatrix, lpm);
}

void CoefFunctionCompound<Complex>::
GetVector( Vector<Complex>& coefVec, const LocPointMapped& lpm ) {
  assert(this->dimType_ == VECTOR );
  Vector<Double> temp, pointCoord;
  Vector<Complex> locVector;

  // First, obtain global coordinates of current point
  lpm.shapeMap->Local2Global(pointCoord,lpm.lp);
  parser_->SetCoordinates( handleReal_, *(this->coordSysDefault_), pointCoord);
  parser_->SetCoordinates( handleImag_, *(this->coordSysDefault_), pointCoord);

  // update internal variables
  UpdateXpr( lpm );

  parser_->EvalVector( handleReal_, temp );
  locVector.Resize(temp.GetSize());
  locVector.SetPart(Global::REAL, temp);
  parser_->EvalVector( handleImag_, temp );
  locVector.SetPart(Global::IMAG, temp);
  if( this->coordSys_ ) {
    // Afterwards rotate the local vector back to global coordinates
    this->coordSys_->Local2GlobalVector( coefVec, locVector, pointCoord );
  } else {
    coefVec = locVector;
  }
}

void CoefFunctionCompound<Complex>::
GetScalar( Complex& coefScalar, const LocPointMapped& lpm ) {
  assert(this->dimType_ == SCALAR);
  Double real, imag;

  // First, obtain global coordinates of current point
  Vector<Double> pointCoord;
  lpm.shapeMap->Local2Global(pointCoord,lpm.lp);
  parser_->SetCoordinates( handleReal_, *(this->coordSysDefault_), pointCoord);
  parser_->SetCoordinates( handleImag_, *(this->coordSysDefault_), pointCoord);

  // update internal variables
  UpdateXpr( lpm );

  real = parser_->Eval( handleReal_ );
  imag = parser_->Eval( handleImag_ );
  coefScalar = Complex(real, imag);
}

void CoefFunctionCompound<Complex>::
SetScalar( std::string& exprReal,
           std::string& exprImag,
           std::map<std::string, PtrCoefFct >& vars ) {
  
  dimType_ = SCALAR;
  
  // loop over all coefficient functions and register them
  std::map<std::string, PtrCoefFct >::iterator it = vars.begin();
  for( ; it != vars.end(); ++it ) {
    RegisterCoefFct( it->first, it->second );
  }
  
  // set expression with mathparser
  exprReal_ = exprReal;
  exprImag_ = exprImag;
  parser_->SetExpr( handleReal_, exprReal_ );
  parser_->SetExpr( handleImag_, exprImag_ );
}

void CoefFunctionCompound<Complex>::
SetVector( StdVector<std::string>& exprReal, 
           StdVector<std::string>& exprImag,
           std::map<std::string, PtrCoefFct >& vars ) {
  dimType_ = VECTOR;

  // loop over all coefficient functions and register them
  std::map<std::string, PtrCoefFct >::iterator it = vars.begin();
  for( ; it != vars.end(); ++it ) {
    RegisterCoefFct( it->first, it->second );
  }

  // set expression with mathparser
  exprReal_ = exprReal.Serialize(',');
  exprImag_ = exprImag.Serialize(',');
  parser_->SetExpr( handleReal_, exprReal_ );
  parser_->SetExpr( handleImag_, exprImag_ );
}

void CoefFunctionCompound<Complex>::
SetTensor( StdVector<std::string>& exprReal,
           StdVector<std::string>& exprImag,
           UInt numRows, UInt numCols,
           std::map<std::string, PtrCoefFct >& vars ) {
  dimType_ = TENSOR;

  // security check: make sure that nRows x nCols = val.size()
  if( numCols * numRows != exprReal.GetSize() ) {
    EXCEPTION("Tensor is supposed to have size " << numRows << " x " 
              << numCols << ", but only " << exprReal.GetSize() 
              << " entries have been provided!");
  }
  this->numRows_ = numRows;
  this->numCols_ = numCols;

  // loop over all coefficient functions and register them
  std::map<std::string, PtrCoefFct >::iterator it = vars.begin();
  for( ; it != vars.end(); ++it ) {
    RegisterCoefFct( it->first, it->second );
  }

  // set expression with mathparser
  exprReal_ = exprReal.Serialize(',');
  exprImag_ = exprImag.Serialize(',');
  parser_->SetExpr( handleReal_, exprReal_ );
  parser_->SetExpr( handleImag_, exprImag_ );
}

UInt CoefFunctionCompound<Complex>::GetVecSize() const {
  assert(this->dimType_ == VECTOR );
  return parser_->GetNumExprs(handleReal_);
}

void CoefFunctionCompound<Complex>::GetTensorSize( UInt& numRows, UInt& numCols ) const {
  assert(this->dimType_ == TENSOR );
  numRows = numRows_;
  numCols = numCols_;
}

std::string CoefFunctionCompound<Complex>::ToString() const {
  std::stringstream out;
  out << "CoefFunctionCompound" << std::endl;
  out << "expression (real): '" << exprReal_ << std::endl;
  out << "expression (imag): '" << exprImag_ << std::endl;
  out << "registered variables:" << std::endl;
  std::map<std::string, PtrCoefFct>::const_iterator it;
  it = coefs_.begin();
  for( ; it != coefs_.end(); ++it ) {
    out << "\tvariable: " << it->first << std::endl;
    out << "\ttype: " 
        << coefDimType.ToString(it->second->GetDimType() ) << std::endl;
    out << "\tvalue: " << it->second->ToString() << std::endl;
  }

  return out.str();
}

void CoefFunctionCompound<Complex>::
RegisterCoefFct( const std::string& name,
                 PtrCoefFct& coef ) {
#ifdef USE_OPENMP
  if(omp_get_num_threads() > 1 ){
    EXCEPTION("Calling from parallel region. Not allowed.")
  }
#endif

  // make sure, that coefficient function with the same name is not already
  // registered
  if( coefs_.find( name ) != coefs_.end() ) {
    EXCEPTION( "A coefficient function with name '" << name 
               << "' was already registered!" );
  }
  
  // query type of coefficient
  CoefDimType dim = coef->GetDimType();
  coefs_[name] = coef;
  coefDimTypes_[name] = dim;
  
  if( dim == SCALAR ) {
    
    // query variables names
    std::string real, imag;
    CoefFunction::GenScalCompNames( real, imag, name, coef );
#pragma omp parallel num_threads(CFS_NUM_THREADS)
    {
#pragma omp critical (CoefFunctionCompound_Double)
      {
    parser_->RegisterExternalVar(handleReal_, real,  &(scalVarsReal_[name]) );
    parser_->RegisterExternalVar(handleReal_, imag,  &(scalVarsImag_[name]) );
    parser_->RegisterExternalVar(handleImag_, real,  &(scalVarsReal_[name]) );
    parser_->RegisterExternalVar(handleImag_, imag,  &(scalVarsImag_[name]) );
      }
    }
    
  } else if( dim == VECTOR ) {
  
    StdVector<std::string> real, imag;
    CoefFunction::GenVecCompNames( real, imag, name, coef );
#pragma omp parallel num_threads(CFS_NUM_THREADS)
    {
#pragma omp critical (CoefFunctionCompound_Double)
      {
    vecVarsReal_[name].Resize(real.GetSize());
    vecVarsImag_[name].Resize(imag.GetSize());
    
    for( UInt i = 0; i < real.GetSize(); ++ i ) {
      parser_->RegisterExternalVar(handleReal_, real[i],  &(vecVarsReal_[name][i]) );
      parser_->RegisterExternalVar(handleImag_, real[i],  &(vecVarsReal_[name][i]) );
    }
    for( UInt i = 0; i < imag.GetSize(); ++ i ) {

      parser_->RegisterExternalVar(handleReal_, imag[i],  &(vecVarsImag_[name][i]) );
      parser_->RegisterExternalVar(handleImag_, imag[i],  &(vecVarsImag_[name][i]) );

    }
      }
    }
  } else if( dim == TENSOR ) {
    StdVector<std::string> real, imag;
    CoefFunction::GenTensorCompNames( real, imag, name, coef );
#pragma omp parallel num_threads(CFS_NUM_THREADS)
    {
#pragma omp critical (CoefFunctionCompound_Double)
      {
    tensorVarsReal_[name].Resize(numRows_, numCols_);
    tensorVarsImag_[name].Resize(numRows_, numCols_);
    
    UInt pos = 0;
    for( UInt i = 0; i < numRows_; ++i ) {
      for( UInt j = 0; j < numCols_; ++j ) {

        parser_->RegisterExternalVar(handleReal_, real[pos],  
                                     &(tensorVarsReal_[name][i][j]) );
        parser_->RegisterExternalVar(handleImag_, real[pos],  
                                     &(tensorVarsReal_[name][i][j]) );

        pos++;
      }
    }
    pos = 0;
    for( UInt i = 0; i < numRows_; ++i ) {
      for( UInt j = 0; j < numCols_; ++j ) {
        parser_->RegisterExternalVar(handleReal_, imag[pos],  
                                     &(tensorVarsImag_[name][i][j]) );
        parser_->RegisterExternalVar(handleImag_, imag[pos],  
                                     &(tensorVarsImag_[name][i][j]) );
        pos++;
      }
    }
      }
    }
  } else {
    EXCEPTION( "Unknown dimension type of coefficient function '" << name << "'")
  }
}

void CoefFunctionCompound<Complex>::
UpdateXpr( const LocPointMapped& lpm ) {

//#ifdef USE_OPENMP
//  if(omp_get_num_threads() == 1 && CFS_NUM_THREADS>1){
//    WARN("Calling from serial region. May be dangerous")
//  }
//#endif

  // loop over all registered coefficients and update their values
  std::map<std::string, PtrCoefFct >::iterator it = coefs_.begin();
  for( ; it != coefs_.end(); ++it ) {
    const std::string & name = it->first;
    // query dimType
    CoefDimType dim = coefDimTypes_[name];

    if( dim == SCALAR ) {
      Complex temp;
      it->second->GetScalar( temp, lpm );
      scalVarsReal_[name] = temp.real();
      scalVarsImag_[name] = temp.imag();
//      std::cerr << "setting " << name << " to (" << temp.real()
//          << ", " << temp.imag()<< ")\n";
    } else if( dim == VECTOR ) {
      Vector<Complex> temp;
      it->second->GetVector( temp, lpm );
      vecVarsReal_[name] = temp.GetPart(Global::REAL);
      vecVarsImag_[name] = temp.GetPart(Global::IMAG);
//      std::cerr << "setting " << name << " to \n\tREAL:"
//          << (temp.GetPart(Global::REAL)).ToString()
//          << "\n\tIMAG: " << (temp.GetPart(Global::IMAG)).ToString() << "\n";
    } else if( dim == TENSOR ) {
      Matrix<Complex> temp;
      it->second->GetTensor( temp, lpm );
      tensorVarsReal_[name] = temp.GetPart(Global::REAL);
      tensorVarsImag_[name] = temp.GetPart(Global::IMAG);
    } else {
      EXCEPTION( "Unknown dimension type of coefficient function '" << name << "'")
    }

  }
}

} // end of namespace

