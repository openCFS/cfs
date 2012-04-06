#include "CoefXpr.hh"

#include "boost/algorithm/string/erase.hpp"

namespace CoupledField {


// initialize static variable
UInt CoefXpr::numVarNames_ = 0;

// ==========================================================================
//  GENERAL STATIC METHODS
// ==========================================================================

UInt CoefXpr::GetNumOperands(OpType op ) {
  switch(op) {
    case NOOP:
      return 0;
      break;
      
    // UNARY FUNCTIONS
    case OP_NORM:
    case OP_SQRT:
      return 1;
      break;
        
    //BINARY FUNCTIONS
    case OP_ADD:
    case OP_SUB:
    case OP_MULT:
    case OP_DIV:
    case OP_CROSS:
    case OP_POW:
      return 2;
      break;
      
    default:
      EXCEPTION("Case not implemented");
      break;
  }
}

std::string CoefXpr::OpToString( OpType op ) {
  switch(op) {
    case NOOP:
      return "";
      break;
    case OP_ADD:
      return "+";
      break;
    case OP_SUB:
      return "-";
      break;
    case OP_MULT:
      return "*";
      break;
    case OP_DIV:
      return "/";
      break;
    default:
      EXCEPTION("No string representation for operand type");
      break;
  }

}

CoefFunction::CoefDimType CoefXpr::GetDimType( PtrCoefFct a,
                                               OpType op ) {

  CoefFunction::CoefDimType dim = CoefFunction::NO_DIM;
  switch( a->GetDimType() ) {
    case CoefFunction::SCALAR:
      // any operation applied to a scalar will return again a scalar
      dim = CoefFunction::SCALAR;
      break;

    case CoefFunction::VECTOR:
      switch( op ) {
        case OP_NORM:
          dim = CoefFunction::SCALAR;
          break;
        default:
          dim = CoefFunction::VECTOR;
          break;
      }  
      break;

    case CoefFunction::TENSOR:
      switch( op ) {
        case OP_NORM:
          dim = CoefFunction::SCALAR;
          break;
        default:
          dim = CoefFunction::TENSOR;
          break;
      }  
      break;

    default:
      EXCEPTION( "Could not determine dimensionality" );
      break;
  }
  return dim;

}

CoefFunction::CoefDimType CoefXpr::GetDimType( PtrCoefFct a,
                                               PtrCoefFct b,
                                               OpType op ) {
  
  CoefFunction::CoefDimType dim = CoefFunction::NO_DIM;
  
  if( a->GetDimType() == b->GetDimType() ) {
    // 1) Both arguments have same dimension
    switch( a->GetDimType() ) {
      
       case CoefFunction::SCALAR:
         // SCALAR-SCALAR CASE
         dim = CoefFunction::SCALAR;
         break;

       case CoefFunction::VECTOR:
         switch( op ) {
           case OP_ADD:
           case OP_SUB:
           case OP_CROSS:
             dim = CoefFunction::VECTOR;
             break;
           case OP_MULT:
             // scalar product of two vectors
             dim = CoefFunction::SCALAR;
             break;
           case OP_DIV:
             EXCEPTION( "Division of two vectors not defined");
             break;
           default:
             EXCEPTION( "Operation for vector type not defined" );
             break;
         }  
         break;
         
       case CoefFunction::TENSOR:
         switch( op ) {
           case OP_ADD:
           case OP_SUB:
           case OP_MULT:
             dim = CoefFunction::TENSOR;
             break;
           case OP_DIV:
             EXCEPTION( "Division of two tensors not defined" );
             break;
           default:
             EXCEPTION( "Operation for vector type not defined" );
             break;
         }  
         break;
               
       default:
         EXCEPTION( "Could not determine dimensionality" );
         break;
    }
  } else {
    // 2) Both arguments have different dimension
    // case: vec - scalar
    if( a->GetDimType() == CoefFunction::VECTOR &&  
        b->GetDimType() == CoefFunction::SCALAR ) {
      dim = CoefFunction::VECTOR;
    }
    // case: tensor - scalar
    else if( a->GetDimType() == CoefFunction::TENSOR &&  
            b->GetDimType() == CoefFunction::SCALAR ) {
      dim = CoefFunction::TENSOR;
    }
    // case: tensor - vec
    else if( a->GetDimType() == CoefFunction::TENSOR &&  
             b->GetDimType() == CoefFunction::VECTOR ) {
      dim = CoefFunction::VECTOR;
    } else {
      EXCEPTION( "Can not apply a function on arguments of type "
          << CoefFunction::CoefDimType_.ToString(a->GetDimType())
          << " and " << CoefFunction::CoefDimType_.ToString(a->GetDimType()) );
    }
    
  }
  
  return dim;
}

bool CoefXpr::IsComplex( PtrCoefFct a,
                         OpType op ) {
  if( op == OP_NORM ) {
    return false;
  } else {
    return a->IsComplex();
  }
}

bool CoefXpr::IsComplex( PtrCoefFct a,
                         PtrCoefFct b,
                         OpType op ) {
  return (a->IsComplex() || b->IsComplex() );
}

void CoefXpr::ApplyUnaryFunc( std::string& retReal, const std::string& argReal,
                              OpType op ) {
  StdVector<std::string> args;
  switch( op ) {
    case OP_NORM:
      args = "( abs(", argReal, ") )";
      retReal = args.Serialize(' ');
      break;
    
    case OP_SQRT:
      args = "( sqrt(", argReal, ") )";
      retReal = args.Serialize(' ');
      break;
    
    default:
      EXCEPTION(" Unknown operand type '" << OpToString(op) 
                << "' to ApplyUnaryFunc" );
      break;
  }
}

void CoefXpr::ApplyUnaryFunc( std::string& retReal, std::string& retImag, 
                              const std::string& argReal,
                              const std::string& argImag,
                              OpType op ) {
  StdVector<std::string> args;
  switch( op ) {
    case OP_NORM:
      args = "( sqrt(", argReal, "*", argReal, "+", argImag, "*", argImag, ") )";
      retReal = args.Serialize(' ');
      retImag = "0.0";
      break;
    
    case OP_SQRT:
      EXCEPTION( "Complex square root not implemented" );
      break;
    
    default:
      EXCEPTION(" Unknown operand type '" << OpToString(op) 
                << "' to ApplyUnaryFunc" );
      break;
  }

}

void CoefXpr::ApplyBinaryFunc( std::string& retReal, 
                              const std::string& arg1Real,
                              const std::string& arg2Real,
                              OpType op ) {
  StdVector<std::string> args;
  switch( op ) {
    case OP_ADD:
    case OP_SUB:
    case OP_MULT:
    case OP_DIV:
      args = "(", arg1Real, ")", OpToString(op), "(", arg2Real, ")";
      retReal = args.Serialize(' ');
      break;
    
    case OP_POW:
      args = "(", arg1Real, "^", arg2Real, ")";
      break;
    
    default:
    EXCEPTION(" Unknown operand type '" << OpToString(op) 
              << "' to ApplyBinaryFunc");
    break;
  }

  // in the end, remove all whitespaces
  //boost::erase_all(retReal, " ");
}

void CoefXpr::ApplyBinaryFunc( std::string& retReal, std::string& retImag, 
                               const std::string& arg1Real,
                               const std::string& arg2Real,
                               const std::string& arg1Imag,
                               const std::string& arg2Imag,
                               OpType op ) {
  StdVector<std::string> args;
  std::string denom;
  switch( op ) {
    case OP_ADD:
    case OP_SUB:
      args = "(", arg1Real, ")", OpToString(op), "(", arg2Real, ")";
      retReal = args.Serialize(' ');
      args = "(", arg1Imag, ")", OpToString(op), "(", arg2Imag, ")";
      retImag = args.Serialize(' ');
      break;

    case OP_MULT:
      args = "(", arg1Real, "*", arg2Real, " - ", arg1Imag, "*", arg2Imag, ")";
      retReal = args.Serialize(' ');
      args = "(", arg1Real, "*", arg2Imag, " + ", arg1Imag, "*", arg2Real, ")";
      retImag = args.Serialize(' ');
      break;

    case OP_DIV:
      args = "( (", arg2Real, "*", arg2Real, ") + (", arg2Imag, "*", arg2Imag, ") )";
      denom = args.Serialize(' ');

      args = "(", arg1Real, "*", arg2Real, "+", arg1Imag, "*", arg2Imag, ")/", denom;
      retReal = args.Serialize(' ');

      args = "(", arg1Imag, "*", arg2Real, "-", arg1Real, "*", arg2Imag, ") / ", denom;
      retImag = args.Serialize(' ');
      break;
      
    case OP_POW:
      EXCEPTION( "Complex POW operation not implemented" )
      break;
      
    default:
      EXCEPTION("Missing operand type");
      break;
  }

  // in the end, remove all whitespaces
  //boost::erase_all(retReal, " ");
  //boost::erase_all(retImag, " ");
  
}

std::string CoefXpr::GetUniqueVarName() {
  
  std::string ret;
  // we use all 26 characters of the alphabet
  UInt base = 26;
  
  // offset for ASCII character (97 = 'a')
  UInt offset = 97;
  UInt rem = 0;
  UInt val = numVarNames_;
  bool cont = true;
  while (cont) {
    // get current remainder
    rem = val % base;
    ret = char(rem+offset) + ret;
    val /= base;
    if( val == 0)
      cont = false;
  }
  
  // increment number of generate variable names
  numVarNames_++;

  // add prefix
  ret= "_"+ret;
  
  return ret;  
}


// ==========================================================================
//  UNARY OPERATION EXPRESSIONS
// ==========================================================================
void CoefXprUnaryOp::Init( PtrCoefFct a, 
                           CoefXpr::OpType op ) {
  if( CoefXpr::GetNumOperands(op) != 1  ) {
     EXCEPTION( "Operand must have exactly 1 arguments" );
   }
   
   dimType_ = CoefXpr::GetDimType( a, op );
   isAnalytical_ = a->IsAnalytic();
   isComplex_ = CoefXpr::IsComplex(a, op );
   a_ = a;
   aName_ = CoefXpr::GetUniqueVarName();
   op_ = op;
}
  
CoefXprUnaryOp::CoefXprUnaryOp( PtrCoefFct a, 
                                CoefXpr::OpType op ) : CoefXpr() {
 Init(a, op);
}

CoefXprUnaryOp::CoefXprUnaryOp( const std::string& a,
                                CoefXpr::OpType op ) : CoefXpr() {
  PtrCoefFct temp  
  = CoefFunction::Generate( Global::REAL, a );
  Init( temp, op );
}

CoefXprUnaryOp::CoefXprUnaryOp( const CoefXpr& a,
                                CoefXpr::OpType op ) : CoefXpr() {
  Global::ComplexPart part = a.IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct temp  
    = CoefFunction::Generate( part, a );
  Init( temp, op );
}

void CoefXprUnaryOp::GetScalarXpr( std::string& real, std::string& imag ) const {
  
  // Switch depending on type of argument
  if( a_->GetDimType() == CoefFunction::SCALAR ) {
    // -------------
    //  SCALAR CASE
    // -------------
    std::string aR, aI;
    
    if( isAnalytical_) {
      CoefFunctionAnalytic & coefA = 
              dynamic_cast<CoefFunctionAnalytic&>(*a_);
      coefA.GetStrScalar( aR, aI );
    } else {
      CoefFunction::GenScalCompNames(aR, aI, aName_, a_);
    }

    if( ! a_->IsComplex() ) {
      CoefXpr::ApplyUnaryFunc( real, aR, op_ );
      imag = "0.0";
    } else {
      CoefXpr::ApplyUnaryFunc( real, imag, aR, aI, op_ );
    }
  }
  // -------------
  //  VECTOR CASE
  // -------------
  else if( a_->GetDimType() == CoefFunction::VECTOR ) {
    
    // only allowed vec->scalar operation is NORM
    if( op_ == OP_NORM ) {
      StdVector<std::string> aR, aI;
      if( isAnalytical_) {
        CoefFunctionAnalytic & coefA = 
            dynamic_cast<CoefFunctionAnalytic&>(*a_);
        coefA.GetStrVector( aR, aI );
      } else {
        CoefFunction::GenVecCompNames(aR, aI, aName_, a_);
      }

      StdVector<std::string> args;
      if( ! a_->IsComplex() ) {
        args = "( sqrt(";
        std::string tmp;
        for( UInt i = 0; i < aR.GetSize(); ++ i ) {
          CoefXpr::ApplyBinaryFunc( tmp, aR[i], aR[i], OP_MULT );
          args.Push_back(tmp);
          args.Push_back("+");
        }
        args.Push_back(") )");
        imag = "0.0";
      } else {
        EXCEPTION("Norm of complex valued vector not implemented yet");
      }
    }
  }
}

void CoefXprUnaryOp::GetVectorXpr( StdVector<std::string>& real, 
                                 StdVector<std::string>& imag ) const {
  
  EXCEPTION( "No vector valued unary function available") 
}

void CoefXprUnaryOp::GetTensorXpr( UInt& numRows, UInt& numCols,
                                 StdVector<std::string>& real, 
                                 StdVector<std::string>& imag ) const {
  EXCEPTION( "No tensor valued unary function available")
}

void CoefXprUnaryOp::
GetArgs( std::map<std::string, PtrCoefFct >& vars ) const {
  vars[aName_] = a_;
}

// ==========================================================================
//  BINARY OPERATION EXPRESSIONS
// ==========================================================================
void CoefXprBinOp::Init( PtrCoefFct a, 
                         PtrCoefFct b,
                         CoefXpr::OpType op ) {
  
  if( a->GetDimType() != b->GetDimType() ) {
     EXCEPTION( "Both arguments have same dimension" );
   }
   if( GetNumOperands(op) != 2  ) {
     EXCEPTION( "Operand must have exactly 2 operands" );
   }
   
   dimType_ = GetDimType( a, b, op);
   isAnalytical_ = a->IsAnalytic() && b->IsAnalytic();
   isComplex_ = a->IsComplex() || b->IsComplex();
   a_ = a;
   b_ = b;
   aName_ = CoefXpr::GetUniqueVarName();
   bName_ = CoefXpr::GetUniqueVarName();
   
   op_ = op;
}

CoefXprBinOp::CoefXprBinOp( PtrCoefFct a, 
                            PtrCoefFct b,
                            CoefXpr::OpType op ) : CoefXpr() {
  Init( a, b, op );
}

CoefXprBinOp::CoefXprBinOp( PtrCoefFct a, 
                            const CoefXpr& b,
                            CoefXpr::OpType op ) : CoefXpr() {
  PtrCoefFct temp 
  = CoefFunction::Generate( Global::REAL, b );
  Init( a, temp, op );
}

CoefXprBinOp::CoefXprBinOp( PtrCoefFct a, 
                            const std::string& b,
                            CoefXpr::OpType op ) : CoefXpr() {
  PtrCoefFct temp  
  = CoefFunction::Generate( Global::REAL, b );
  Init( a, temp, op );
}

CoefXprBinOp::CoefXprBinOp( const std::string& a, 
                            PtrCoefFct b, 
                            CoefXpr::OpType op ) : CoefXpr() {
  
  PtrCoefFct temp  
  = CoefFunction::Generate( Global::REAL, a );
  Init( temp, b, op );
  
}

CoefXprBinOp::CoefXprBinOp( const std::string& a, 
                            const CoefXpr& b, 
                            CoefXpr::OpType op ) : CoefXpr() {
  PtrCoefFct tempA  
   = CoefFunction::Generate( Global::REAL, a );
  Global::ComplexPart part = b.IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct tempB  
    = CoefFunction::Generate( part, b );
  Init( tempA, tempB, op );
}

void CoefXprBinOp::GetScalarXpr( std::string& real, std::string& imag ) const {
  
  assert((this->dimType_ == CoefFunction::NO_DIM) || 
         (this->dimType_ == CoefFunction::SCALAR) );

  if( a_->GetDimType() == CoefFunction::SCALAR  &&
      b_->GetDimType() == CoefFunction::SCALAR ) {
    
    // --------------------
    //  SCALAR-SCALAR CASE
    // --------------------
    
    std::string aR, aI, bR, bI;
    if( isAnalytical_) {
      CoefFunctionAnalytic & coefA = 
          dynamic_cast<CoefFunctionAnalytic&>(*a_);
      CoefFunctionAnalytic & coefB = 
          dynamic_cast<CoefFunctionAnalytic&>(*b_);
      coefA.GetStrScalar( aR, aI );
      coefB.GetStrScalar( bR, bI );
    } else {
      CoefFunction::GenScalCompNames(aR, aI, aName_, a_);
      CoefFunction::GenScalCompNames(bR, bI, bName_, b_);
    }

    if( !isComplex_ ) {
      CoefXpr::ApplyBinaryFunc( real, aR, bR, op_ );
      imag = "0.0";
    } else {
      CoefXpr::ApplyBinaryFunc( real, imag, aR, bR, aI, bI, op_ );
    }
  } else if( a_->GetDimType() == CoefFunction::VECTOR  &&
             b_->GetDimType() == CoefFunction::VECTOR ) {  
    // --------------------
    //  VECTOR-VECTOR CASE
    // --------------------
    
    // Only known (vec,vec)->scalar operation is scalar product
    if( op_ == OP_MULT ) {
      StdVector<std::string> aR, aI, bR, bI;
      if( isAnalytical_) {
        CoefFunctionAnalytic & coefA = 
            dynamic_cast<CoefFunctionAnalytic&>(*a_);
        CoefFunctionAnalytic & coefB = 
            dynamic_cast<CoefFunctionAnalytic&>(*b_);
        coefA.GetStrVector( aR, aI );
        coefB.GetStrVector( bR, bI );
      } else {
        CoefFunction::GenVecCompNames(aR, aI, aName_, a_);
        CoefFunction::GenVecCompNames(bR, bI, bName_, b_);
      }
      UInt numEntries = aR.GetSize();
      
      if( !isComplex_ ) {
        std::string temp;
        // treat first entry
        CoefXpr::ApplyBinaryFunc( real, aR[0], bR[0], OP_MULT );
        for( UInt i = 1; i < numEntries; ++ i ) {
          CoefXpr::ApplyBinaryFunc( temp, aR[i], bR[i], op_ );
          real += "+" + temp;
        }
        imag = "0.0";
      } else {
        std::string tempReal, tempImag;
        // treat first entry
        CoefXpr::ApplyBinaryFunc( real, imag, aR[0], bR[0], 
                                  aI[0], bI[0], OP_MULT );
        for( UInt i = 1; i < numEntries; ++ i ) {
          CoefXpr::ApplyBinaryFunc( tempReal, tempImag , 
                                    aR[i], bR[i],aI[i], bI[i], OP_MULT );
          real += "+" + tempReal;
          imag += "+" + tempImag;
        }
      }
        
    } else {
      EXCEPTION("The only allowed binary vector-vector operation "
                << "is OP_MULT" ); 
    } // if OP_MULT
  } else {
    EXCEPTION("Arguments must be either both of scalar or vector type.")
  }
}

void CoefXprBinOp::GetVectorXpr( StdVector<std::string>& real, 
                                 StdVector<std::string>& imag ) const {
  assert((this->dimType_ == CoefFunction::NO_DIM) || 
         (this->dimType_ == CoefFunction::VECTOR) );
  
  // Possible operations:
  // +
  // -
  // cross

  if( a_->GetDimType() == CoefFunction::VECTOR  &&
      b_->GetDimType() == CoefFunction::VECTOR ) {
    StdVector<std::string> aR, aI, bR, bI;
    if( isAnalytical_) {
      CoefFunctionAnalytic & coefA = 
          dynamic_cast<CoefFunctionAnalytic&>(*a_);
      CoefFunctionAnalytic & coefB = 
          dynamic_cast<CoefFunctionAnalytic&>(*b_);
      coefA.GetStrVector( aR, aI );
      coefB.GetStrVector( bR, bI );
    } else {
      CoefFunction::GenVecCompNames(aR, aI, aName_, a_);
      CoefFunction::GenVecCompNames(bR, bI, bName_, b_);
    }
    UInt numEntries = aR.GetSize();
    real.Resize( numEntries );
    imag.Resize( numEntries );

    // switch depending on operation type
    // === OP_ADD / OP_SUB ===
    if( op_ == OP_ADD || op_ == OP_SUB ) {
      if( !isComplex_ ) {
        for( UInt i = 0; i < numEntries; ++ i ) {
          CoefXpr::ApplyBinaryFunc( real[i], aR[i], bR[i], op_ );
          imag[i] = "0.0";
        }
      } else {
        for( UInt i = 0; i < numEntries; ++ i ) {
          CoefXpr::ApplyBinaryFunc( real[i], imag[i], aR[i], bR[i],
                                    aI[i], bI[i], op_ );
        }
      }
    }
    
    // === OP_CROSS ===
    else if( op_ == OP_CROSS ) {
      
      // switch, depending on dimension
      if( numEntries == 1 ) {
        
        // 1) simple product ob both entries
        if(!isComplex_) {
          CoefXpr::ApplyBinaryFunc( real[0], aR[0], bR[0], OP_MULT );
        } else {
          CoefXpr::ApplyBinaryFunc( real[0], imag[0], aR[0], bR[0],
                                    aI[0], bI[0], OP_MULT );
        }
       } else if( numEntries == 3) {
         
        // 2) full 3D cross product
         if(!isComplex_) {
           std::string temp;
           // 1st row
           CoefXpr::ApplyBinaryFunc( temp, aR[1], bR[2], OP_MULT );
           real[0] = temp + " - ";
           CoefXpr::ApplyBinaryFunc( temp, aR[2], bR[1], OP_MULT );
           real[0] += temp;
           // 2nd row
           CoefXpr::ApplyBinaryFunc( temp, aR[2], bR[0], OP_MULT );
           real[1] = temp + " - ";
           CoefXpr::ApplyBinaryFunc( temp, aR[0], bR[2], OP_MULT );
           real[1] += temp;
           // 3rd row
           CoefXpr::ApplyBinaryFunc( temp, aR[0], bR[1], OP_MULT );
           real[2] = temp + " - ";
           CoefXpr::ApplyBinaryFunc( temp, aR[1], bR[0], OP_MULT );
           real[2] += temp;
         } else {
           std::string tempR, tempI;
           // 1st row
           CoefXpr::ApplyBinaryFunc( tempR, tempI, 
                                     aR[1], bR[2], 
                                     aI[1], bI[2], OP_MULT );
           real[0] = tempR + " - ";
           imag[0] = tempI + " - ";
           CoefXpr::ApplyBinaryFunc( tempR, tempI, 
                                     aR[2], bR[1], 
                                     aI[2], bI[1], OP_MULT );
           real[0] += tempR;
           imag[0] += tempI;
           // 2nd row
           CoefXpr::ApplyBinaryFunc( tempR, tempI, 
                                     aR[2], bR[0], 
                                     aI[2], bI[0], OP_MULT );
           real[1] = tempR + " - ";
           imag[1] = tempI + " - ";
           CoefXpr::ApplyBinaryFunc( tempR, tempI, 
                                     aR[0], bR[2], 
                                     aI[0], bI[2], OP_MULT );
           real[1] += tempR;
           imag[1] += tempI;
           // 3rd row
           CoefXpr::ApplyBinaryFunc( tempR, tempI, 
                                     aR[0], bR[1], 
                                     aI[0], bI[1], OP_MULT );
           real[2] = tempR + " - ";
           imag[2] = tempI + " - ";
           CoefXpr::ApplyBinaryFunc( tempR, tempI, 
                                     aR[1], bR[0], 
                                     aI[1], bI[0], OP_MULT );
           real[2] += tempR;
           imag[2] += tempI;
         }
         
      } else 
        EXCEPTION( "Vector / cross product only defined in 1D and 3D");

    } else {
      EXCEPTION( "The only allowed (vector,vector)->vector operations "
                << "are OP_ADD, OP_SUB and OP_CROSS" );
    }
    
  } else {
    EXCEPTION("Arguments must be both of vector type.")
  }
  
  
}

void CoefXprBinOp::GetTensorXpr( UInt& numRows, UInt& numCols,
                                 StdVector<std::string>& real, 
                                 StdVector<std::string>& imag ) const {
  assert((this->dimType_ == CoefFunction::NO_DIM) || 
         (this->dimType_ == CoefFunction::TENSOR) );
   
  // Possible operations:
  // +
  // -
  // *
  if( a_->GetDimType() == CoefFunction::TENSOR  &&
      b_->GetDimType() == CoefFunction::TENSOR ) {

    StdVector<std::string> aR, aI, bR, bI;
     UInt numRowsA, numRowsB, numColsA, numColsB;
     
     if( isAnalytical_) {
       CoefFunctionAnalytic & coefA = 
           dynamic_cast<CoefFunctionAnalytic&>(*a_);
       CoefFunctionAnalytic & coefB = 
           dynamic_cast<CoefFunctionAnalytic&>(*b_);
       coefA.GetStrTensor( numRowsA, numColsA, aR, aI );
       coefB.GetStrTensor( numRowsB, numColsB, bR, bI );
     } else {
       CoefFunction::GenTensorCompNames(aR, aI, aName_, a_);
       CoefFunction::GenTensorCompNames(bR, bI, bName_, b_);
       a_->GetTensorSize(numRowsA, numColsA);
       b_->GetTensorSize(numRowsB, numColsB);
     }

     // switch depending on operand type
     if( op_ == OP_ADD || op_ == OP_SUB ) {

       // ensure that both tensors have the same size
       if( numRowsA != numRowsB || numColsA != numColsB ) {
         EXCEPTION( "Both tensors must have the same dimension");
       }
       numRows = numRowsA;
       numCols = numColsA;
       UInt numEntries = numRows * numRows;
       real.Resize( numEntries );
       imag.Resize( numEntries );
       if( !isComplex_ ) {
         for( UInt i = 0; i < numEntries; ++ i ) {
           CoefXpr::ApplyBinaryFunc( real[i], aR[i], bR[i], op_ );
           imag[i] = "0.0";
         }
       } else {
         for( UInt i = 0; i < numEntries; ++ i ) {
           CoefXpr::ApplyBinaryFunc( real[i], imag[i], aR[i], bR[i],
                                     aI[i], bI[i], op_ );
         }
       }
     } else  if( op_ == OP_MULT ) {
       
       // ensure that both tensors can be "multiplied"
       if( numColsA != numRowsB  ) {
         EXCEPTION( "Tensors can not be multiplied: First tensor is of size "
             << numRowsA << " x " << numColsA << " and second one is of size "
             << numRowsB << " x " << numColsB  );
       }
       numRows = numRowsA;
       numCols = numColsB;
       real.Resize( numRows * numCols );
       imag.Resize( numRows * numCols );

       UInt posA, posB, posRet;
       posRet = 0;
       if( !isComplex_ ) {
         for( UInt i = 0; i < numRows; ++i ) {
           for( UInt j = 0; j < numCols; ++j ) {
             std::string temp1, temp2;
             posA = i*numCols;
             posB = j;
             CoefXpr::ApplyBinaryFunc( temp1, aR[posA], bR[posB], OP_MULT );
             for( UInt k = 1; k < numColsA; ++k ) {
               posA ++;
               posB += numCols;
               CoefXpr::ApplyBinaryFunc( temp2, aR[posA], bR[posB], OP_MULT );
               temp1 += " + " + temp2;
             } // numColsA
             real[posRet++] = temp1;
           } // numCols
         } // numRows
         imag.Init("0.0");
       } else {
         for( UInt i = 0; i < numRows; ++i ) {
           for( UInt j = 0; j < numCols; ++j ) {
             std::string temp1R, temp1I, temp2R, temp2I;
             posA = i*numCols;
             posB = j;
             CoefXpr::ApplyBinaryFunc( temp1R, temp1I, aR[posA], bR[posB], 
                                       aI[posA], bI[posB], OP_MULT );
             for( UInt k = 1; k < numColsA; ++k ) {
               posA ++;
               posB += numCols;
               CoefXpr::ApplyBinaryFunc( temp2R, temp2I, aR[posA], bR[posB], 
                                         aI[posA], bI[posB], OP_MULT );
               temp1R += " + " + temp2R;
               temp1I += " + " + temp2I;
             } // numColsA
             real[posRet] = temp1R;
             imag[posRet] = temp1I;
             posRet++;
           } // numCols
         } // numRows
       }
       
       
     } else {
       EXCEPTION( "The only allowed (tensor,tensor)->tensor operations "
                << "are OP_ADD, OP_SUB and OP_MULT");
     }
  } else {
    EXCEPTION("Arguments must be both of tensor type.")
  }
}

void CoefXprBinOp::
GetArgs( std::map<std::string, PtrCoefFct >& vars ) const {
  vars[aName_] = a_;
  vars[bName_] = b_;
}

// ==========================================================================
//  VECTOR-SCALAR OPERATION EXPRESSIONS
// ==========================================================================
void CoefXprVecScalOp::Init( PtrCoefFct a, 
                             PtrCoefFct b,
                             CoefXpr::OpType op ) {

  // check dimensionality
  if( a->GetDimType() != CoefFunction::VECTOR ) {
    EXCEPTION( "First argument must be of type VECTOR" );
  }

  if( b->GetDimType() != CoefFunction::SCALAR ) {
    EXCEPTION( "Second argument must be of type SCALAR" );
  }

  if( GetNumOperands(op) != 2  ) {
    EXCEPTION( "Operand must have exactly 2 operands" );
  }

  dimType_ = a->GetDimType();
  isAnalytical_ = a->IsAnalytic() && b->IsAnalytic();
  isComplex_ = a->IsComplex() || b->IsComplex();
  a_ = a;
  b_ = b;
  aName_ = CoefXpr::GetUniqueVarName();
  bName_ = CoefXpr::GetUniqueVarName();
  op_ = op;
}

CoefXprVecScalOp::CoefXprVecScalOp( PtrCoefFct a, 
                                    PtrCoefFct b,
                                    CoefXpr::OpType op ) : CoefXpr() {
  Init( a, b, op );
}

CoefXprVecScalOp::CoefXprVecScalOp( PtrCoefFct a, 
                                    const std::string& b,
                                    CoefXpr::OpType op ) : CoefXpr() {
  PtrCoefFct temp  
  = CoefFunction::Generate( Global::REAL, b );
  Init( a, temp, op );
}

CoefXprVecScalOp::CoefXprVecScalOp( PtrCoefFct a, 
                                    const CoefXpr& b, 
                                    CoefXpr::OpType op ) : CoefXpr() {

  PtrCoefFct temp  
  = CoefFunction::Generate( Global::REAL, b );
  Init( a, temp, op );
  
}


void CoefXprVecScalOp::GetScalarXpr( std::string& real, std::string& imag ) const {

  EXCEPTION( "Vector-Scalar-Operation provides no scalar results")
}

void CoefXprVecScalOp::GetVectorXpr( StdVector<std::string>& real, 
                                  StdVector<std::string>& imag ) const {
  
  assert((this->dimType_ == CoefFunction::NO_DIM) || 
         (this->dimType_ == CoefFunction::VECTOR) );
  
  StdVector<std::string> aR, aI;
  std::string bR, bI;
  
  if( isAnalytical_) {
    CoefFunctionAnalytic & coefA = 
        dynamic_cast<CoefFunctionAnalytic&>(*a_);
    CoefFunctionAnalytic & coefB = 
        dynamic_cast<CoefFunctionAnalytic&>(*b_);
    coefA.GetStrVector( aR, aI );
    coefB.GetStrScalar( bR, bI );
  } else {
    CoefFunction::GenVecCompNames(aR, aI, aName_, a_);
    CoefFunction::GenScalCompNames(bR, bI, bName_, b_);
  }
  
  UInt numEntries = aR.GetSize();
  real.Resize( numEntries );
  imag.Resize( numEntries );
  if( !isComplex_ ) {
    for( UInt i = 0; i < numEntries; ++ i ) {
      CoefXpr::ApplyBinaryFunc( real[i], aR[i], bR, op_ );
      imag[i] = "0.0";
    }
  } else {
    for( UInt i = 0; i < numEntries; ++ i ) {
      CoefXpr::ApplyBinaryFunc( real[i], imag[i], aR[i], bR,
                                aI[i], bI, op_ );
    }
  }
}

void CoefXprVecScalOp::GetTensorXpr( UInt& numRows, UInt& numCols,
                                     StdVector<std::string>& real, 
                                     StdVector<std::string>& imag ) const {

  EXCEPTION( "Vector-Scalar-Operation provides no tensor results")
}

void CoefXprVecScalOp::
GetArgs( std::map<std::string, PtrCoefFct >& vars ) const {
  vars[aName_] = a_;
  vars[bName_] = b_;
}


// ==========================================================================
//  TENSOR-SCALAR OPERATION EXPRESSIONS
// ==========================================================================
void CoefXprTensScalOp::Init( PtrCoefFct a, 
                             PtrCoefFct b,
                             CoefXpr::OpType op ) {

  // check dimensionality
  if( a->GetDimType() != CoefFunction::TENSOR ) {
    EXCEPTION( "First argument must be of type TENSOR" );
  }

  if( b->GetDimType() != CoefFunction::SCALAR ) {
    EXCEPTION( "Second argument must be of type SCALAR" );
  }

  if( GetNumOperands(op) != 2  ) {
    EXCEPTION( "Operand must have exactly 2 operands" );
  }

  dimType_ = a->GetDimType();
  isAnalytical_ = a->IsAnalytic() && b->IsAnalytic();
  isComplex_ = a->IsComplex() || b->IsComplex();
  a_ = a;
  b_ = b;
  aName_ = CoefXpr::GetUniqueVarName();
  bName_ = CoefXpr::GetUniqueVarName();
  op_ = op;
}

CoefXprTensScalOp::CoefXprTensScalOp( PtrCoefFct a, 
                                      PtrCoefFct b,
                                      CoefXpr::OpType op ) : CoefXpr() {
  Init( a, b, op );
}

CoefXprTensScalOp::CoefXprTensScalOp( PtrCoefFct a, 
                                      const std::string& b,
                                      CoefXpr::OpType op ) : CoefXpr() {
  PtrCoefFct temp  
  = CoefFunction::Generate( Global::REAL, b );
  Init( a, temp, op );
}

CoefXprTensScalOp::CoefXprTensScalOp( PtrCoefFct a, 
                                      const CoefXpr& b, 
                                      CoefXpr::OpType op ) : CoefXpr() {
  PtrCoefFct temp  
  = CoefFunction::Generate( Global::REAL, b );
  Init( a, temp, op );
}


void CoefXprTensScalOp::GetScalarXpr( std::string& real, std::string& imag ) const {

  EXCEPTION( "Tensor-Scalar-Operation provides no scalar results")
}

void CoefXprTensScalOp::GetVectorXpr( StdVector<std::string>& real, 
                                  StdVector<std::string>& imag ) const {
  EXCEPTION( "Tensor-Scalar-Operation provides no vector results")
 
}

void CoefXprTensScalOp::GetTensorXpr( UInt& numRows, UInt& numCols,
                                     StdVector<std::string>& real, 
                                     StdVector<std::string>& imag ) const {

  assert((this->dimType_ == CoefFunction::NO_DIM) || 
         (this->dimType_ == CoefFunction::TENSOR) );

  StdVector<std::string> aR, aI;
  std::string bR, bI;
  
  if( isAnalytical_) {
    CoefFunctionAnalytic & coefA = 
        dynamic_cast<CoefFunctionAnalytic&>(*a_);
    CoefFunctionAnalytic & coefB = 
        dynamic_cast<CoefFunctionAnalytic&>(*b_);
    coefA.GetStrTensor( numRows, numCols, aR, aI );
    coefB.GetStrScalar( bR, bI );
  } else {
    CoefFunction::GenTensorCompNames(aR, aI, aName_, a_);
    CoefFunction::GenScalCompNames(bR, bI, bName_, b_);
    a_->GetTensorSize(numRows, numCols);
  }

  UInt numEntries = aR.GetSize();
  real.Resize( numEntries );
  imag.Resize( numEntries );
  
  if( !isComplex_ ) {
    for( UInt i = 0; i < numEntries; ++ i ) {
      CoefXpr::ApplyBinaryFunc( real[i], aR[i], bR, op_ );
      imag[i] = "0.0";
    }
  } else {
    for( UInt i = 0; i < numEntries; ++ i ) {
      CoefXpr::ApplyBinaryFunc( real[i], imag[i], aR[i], bR,
                                aI[i], bI, op_ );
    }
  }
}

void CoefXprTensScalOp::
GetArgs( std::map<std::string, PtrCoefFct >& vars ) const {
  vars[aName_] = a_;
  vars[bName_] = b_;
}

} // end of namespace
