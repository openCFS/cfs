#include "CoefXpr.hh"

#include "boost/algorithm/string/erase.hpp"

namespace CoupledField {


// define local helper function for wrapping string expressions
inline std::string B(const std::string& xpr ) {
  return "(" + xpr + ")";
}

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
        
    // BINARY FUNCTIONS
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
  return 0;
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
  return "";

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

void CoefXpr::Tranpose( UInt nRows, UInt nCols, 
                        const StdVector<std::string>& in,
                        StdVector<std::string>& trans ) {
  trans.Resize( in.GetSize() );
  
  for( UInt iRow = 0; iRow < nRows; ++iRow ) {
    for( UInt iCol = 0; iCol < nCols; ++iCol ) {  
      trans[iCol*nRows+iRow] = in[iRow*nCols + iCol];
    }
  }
}


void CoefXpr::ApplyUnaryFunc( std::string& retReal, const std::string& argReal,
                              OpType op ) {
  StdVector<std::string> args;
  switch( op ) {
    case OP_NORM:
      args =  "abs(", B(argReal), ")";
      retReal = B(args.Serialize(' '));
      break;
    
    case OP_SQRT:
      args = "sqrt(", B(argReal), ")";
      retReal = B(args.Serialize(' '));
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
      args = "sqrt(", B(argReal), "*", B(argReal), "+", B(argImag), "*", B(argImag), ")";
      retReal = B(args.Serialize(' '));
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
      args =  B(arg1Real), OpToString(op), B(arg2Real);
      retReal = B(args.Serialize(' ') );
      break;
    
    case OP_POW:
      args = B(arg1Real), "^", B(arg2Real);
      retReal = B(args.Serialize(' ' ));
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
      args = B(arg1Real), OpToString(op), B(arg2Real);
      retReal = B(args.Serialize(' ') );
      args = B(arg1Imag), OpToString(op), B(arg2Imag);
      retImag = B(args.Serialize(' '));
      break;

    case OP_MULT:
      args = "(", B(arg1Real), "*", B(arg2Real), ") - (", B(arg1Imag), "*", B(arg2Imag), ")";
      retReal = B(args.Serialize(' '));
      args = "(", B(arg1Real), "*", B(arg2Imag), ") + (", B(arg1Imag), "*", B(arg2Real), ")";
      retImag = B(args.Serialize(' '));
      break;

    case OP_DIV:
      args = "(", B(arg2Real), "*", B(arg2Real), ") + (", B(arg2Imag), "*", B(arg2Imag), ")";
      denom = B(args.Serialize(' '));

      args = "((", B(arg1Real), "*", B(arg2Real), ")+(", B(arg1Imag), "*", B(arg2Imag), "))/(", denom,")";
      retReal = B(args.Serialize(' '));

      args = "((", B(arg1Imag), "*", B(arg2Real), ")-(", B(arg1Real), "*", B(arg2Imag), ")) / (", denom,")";
      retImag = B(args.Serialize(' '));
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

void CoefXpr::ApplyTernaryFunc( std::string& retReal, 
                                const std::string& arg1Real,
                                const std::string& arg2Real,
                                const std::string& arg3Real,
                                OpType op1, OpType op2 ) {
  // first, compute 
  // tmp = op1(arg2, arg3)
  std::string tmp;
  CoefXpr::ApplyBinaryFunc( tmp, arg2Real, arg3Real, op2 );
  
  // then compute ret = op2(arg1, tmp)
  CoefXpr::ApplyBinaryFunc( retReal, arg1Real, tmp, op1 );
}

void CoefXpr::ApplyTernaryFunc( std::string& retReal, std::string& retImag, 
                                const std::string& arg1Real,
                                const std::string& arg2Real,
                                const std::string& arg3Real,
                                const std::string& arg1Imag,
                                const std::string& arg2Imag,
                                const std::string& arg3Imag,
                                OpType op1, OpType op2 ) {
  // first, compute 
  // tmp = op1(arg2, arg3)
  std::string tmpR, tmpI;
  CoefXpr::ApplyBinaryFunc( tmpR, tmpI, arg2Real, arg3Real, 
                            arg2Imag, arg3Imag, op2 );

  // then compute ret = op2(arg1, tmp)
  CoefXpr::ApplyBinaryFunc( retReal, retImag, arg1Real, tmpR, 
                            arg1Imag, tmpI, op1 );
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
  
  PtrCoefFct temp = CoefFunction::Generate( Global::REAL, a );
  Init( temp, op );
}

CoefXprUnaryOp::CoefXprUnaryOp( const CoefXpr& a,
                                CoefXpr::OpType op ) : CoefXpr() {
  
  Global::ComplexPart part = a.IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct temp = CoefFunction::Generate( part, a );
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
  
  Global::ComplexPart part = b.IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct temp  = CoefFunction::Generate( part, b );
  Init( a, temp, op );
}

CoefXprBinOp::CoefXprBinOp( PtrCoefFct a, 
                            const std::string& b,
                            CoefXpr::OpType op ) : CoefXpr() {
  
  Global::ComplexPart part = a->IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct temp = CoefFunction::Generate( part, b );
  Init( a, temp, op );
}

CoefXprBinOp::CoefXprBinOp( const std::string& a, 
                            PtrCoefFct b, 
                            CoefXpr::OpType op ) : CoefXpr() {
  
  Global::ComplexPart part = b->IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct temp = CoefFunction::Generate( part, a );
  Init( temp, b, op );
  
}

CoefXprBinOp::CoefXprBinOp( const std::string& a, 
                            const CoefXpr& b, 
                            CoefXpr::OpType op ) : CoefXpr() {
  Global::ComplexPart part = b.IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct tempA  = CoefFunction::Generate( part, a );
  PtrCoefFct tempB  = CoefFunction::Generate( part, b );
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
      // put brackets around expression
      real = Bracket(real);
      imag = Bracket(imag);
        
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
       UInt numEntries = numRows * numCols;
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
             posA = i*numColsA;
             posB = j;
             CoefXpr::ApplyBinaryFunc( temp1, aR[posA], bR[posB], OP_MULT );
             for( UInt k = 1; k < numColsA; ++k ) {
               posA ++;
               posB += numCols;
               CoefXpr::ApplyBinaryFunc( temp2, aR[posA], bR[posB], OP_MULT );
               temp1 += " + " + temp2;
             } // numColsA
             real[posRet++] = Bracket(temp1);
           } // numCols
         } // numRows
         imag.Init("0.0");
       } else {
         for( UInt i = 0; i < numRows; ++i ) {
           for( UInt j = 0; j < numCols; ++j ) {
             std::string temp1R, temp1I, temp2R, temp2I;
             posA = i*numColsA;
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
             real[posRet] = Bracket(temp1R);
             imag[posRet] = Bracket(temp1I);
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
  
  Global::ComplexPart part = a->IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct temp = CoefFunction::Generate( part, b );
  Init( a, temp, op );
}

CoefXprVecScalOp::CoefXprVecScalOp( PtrCoefFct a, 
                                    const CoefXpr& b, 
                                    CoefXpr::OpType op ) : CoefXpr() {

  Global::ComplexPart part;
  if( a->IsComplex() || b.IsComplex() ) {
    part = Global::COMPLEX;
  } else {
    part = Global::REAL;
  }
  PtrCoefFct temp = CoefFunction::Generate( part, b );
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
  
  Global::ComplexPart part = a->IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct temp = CoefFunction::Generate( part, b );
  Init( a, temp, op );
}

CoefXprTensScalOp::CoefXprTensScalOp( PtrCoefFct a, 
                                      const CoefXpr& b, 
                                      CoefXpr::OpType op ) : CoefXpr() {
  
  Global::ComplexPart part;
  if( a->IsComplex() || b.IsComplex() ) {
    part =  Global::COMPLEX;
  } else {
    part = Global::REAL;
  }
  PtrCoefFct temp = CoefFunction::Generate( part, b );
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

// ==========================================================================
//  TENSOR REPRESENTATIONS (MECHANIC)
// ==========================================================================
void CoefXprMechSubTensor::Init( PtrCoefFct a ) {
  // check dimensionality
  if( a->GetDimType() != CoefFunction::TENSOR ) {
    EXCEPTION( "Argument must be of type TENSOR" );
  }
  
  // ensure that dimension is a full 6x6 tensor
  UInt numRowsA, numColsA;
  a->GetTensorSize(numRowsA, numColsA);
  if( numRowsA != 6 || numColsA != 6 ) {
    EXCEPTION( "Tensor must have dimension 6 x 6 " );
  }
  
  dimType_ = CoefFunction::TENSOR;
  isAnalytical_ = a->IsAnalytic();
  isComplex_ = a->IsComplex();
  a_ = a;
  aName_ = CoefXpr::GetUniqueVarName();
  tensorType_ = NO_TENSOR;
  transposed_ = false;
}

CoefXprMechSubTensor::CoefXprMechSubTensor( PtrCoefFct a ) {
  Init(a);
}
CoefXprMechSubTensor::CoefXprMechSubTensor( const CoefXpr& a)  {
  
  Global::ComplexPart part = a.IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct temp = CoefFunction::Generate( part, a );  
}

void CoefXprMechSubTensor::SetSubTensorType(SubTensorType subType,
                                            bool transpose ) {
  tensorType_ = subType;
  transposed_ = transpose;
}


void CoefXprMechSubTensor::GetTensorXpr( UInt& numRows, UInt& numCols,
                                     StdVector<std::string>& real, 
                                     StdVector<std::string>& imag ) const {
  // ensure, that a subtype is set
  if( tensorType_ == NO_TENSOR ) {
    EXCEPTION( "No tensor sub type was set");
  }

  StdVector<std::string> aR, aI;
  UInt numRowsA, numColsA;

  if( isAnalytical_) {
    CoefFunctionAnalytic & coefA = 
        dynamic_cast<CoefFunctionAnalytic&>(*a_);
    coefA.GetStrTensor( numRowsA, numColsA, aR, aI );
  } else {
    CoefFunction::GenTensorCompNames(aR, aI, aName_, a_);
    a_->GetTensorSize(numRowsA, numColsA);
  }

  // switch depending on subTensorType
  if( tensorType_ == AXI ) {
    // resize to 4x4 tensor
    numCols = 4;
    numRows = 4;
    real.Resize( numRows * numCols );
    imag.Resize( numRows * numCols );
    imag.Init("0.0");
    UInt rowPtr[] = {1,2,6,3};
    for(UInt i = 0; i < numRows; ++i ) {
      for(UInt j = 0; j < numCols; ++j ) {
        real[i*numCols+j] = aR[(rowPtr[i]-1)*6 + (rowPtr[j]-1)];
        if(isComplex_ )
          imag[i*numCols+j] = aI[(rowPtr[i]-1)*6 + (rowPtr[j]-1)];
      }
    }
    
  } else if( tensorType_ == PLANE_STRAIN ) {
    // resize to 4x4 tensor
    numCols = 3;
    numRows = 3;
    real.Resize( numRows * numCols );
    imag.Resize( numRows * numCols );
    imag.Init("0.0");
    UInt rowPtr[] = {1,2,6,3};
    for(UInt i = 0; i < numRows; ++i ) {
      for(UInt j = 0; j < numCols; ++j ) {
        real[i*numCols+j] = aR[(rowPtr[i]-1)*6 + (rowPtr[j]-1)];
        if(isComplex_ )
          imag[i*numCols+j] = aI[(rowPtr[i]-1)*6 + (rowPtr[j]-1)];
      }
    }
  } else if( tensorType_ == PLANE_STRESS)  {
    numCols = 3;
    numRows = 3;
    real.Resize( numRows * numCols );
    imag.Resize( numRows * numCols );
    imag.Init("0.0");
    
    // 1st part: compute same as plane strain representation
    UInt rowPtr[] = {1,2,6};
    for(UInt i = 0; i < numRows; ++i ) {
      for(UInt j = 0; j < numCols; ++j ) {
        real[i*numCols+j] = aR[(rowPtr[i]-1)*6 + (rowPtr[j]-1)];
        if(isComplex_ )
          imag[i*numCols+j] = aI[(rowPtr[i]-1)*6 + (rowPtr[j]-1)];
      }
    }
    
    // 2nd part: subtract entries
    std::string& c22R = aR[6*2+2]; //c[2][2] - real
    std::string& c22I = aI[6*2+2]; //c[2][2] - imag
    
    std::string& c02R = aR[6*0+2]; //c[2][0] - real
    std::string& c02I = aI[6*0+2]; //c[2][0] - imag
    
    std::string& c12R = aR[6*1+2]; //c[1][2] - real
    std::string& c12I = aI[6*1+2]; //c[1][2] - imag
    
    std::string& c20R = aR[6*2+0]; //c[2][0] - real
    std::string& c20I = aI[6*2+0]; //c[2][0] - imag
    
    std::string& c21R = aR[6*2+1]; //c[2][1] - real
    std::string& c21I = aI[6*2+1]; //c[2][1] - imag
    
    std::string tmpR, tmpI;
    
    // modify c[0][0]
    CoefXpr::ApplyTernaryFunc( tmpR, tmpI, c20R, c02R, c22R, 
                               c20I, c02I, c22I, OP_MULT, OP_DIV );
    real[3*0+0] += "- (" + tmpR + ")";
    imag[3*0+0] += "- (" + tmpI + ")";
    // modify c[0][1]
    CoefXpr::ApplyTernaryFunc( tmpR, tmpI, c21R, c02R, c22R, 
                               c21I, c02I, c22I, OP_MULT, OP_DIV );
    real[3*0+1] += "- (" + tmpR + ")";
    imag[3*0+1] += "- (" + tmpI + ")";
    // modify c[1][0]
    CoefXpr::ApplyTernaryFunc( tmpR, tmpI, c20R, c12R, c22R, 
                               c20I, c12I, c22I, OP_MULT, OP_DIV );
    real[3*1+0] += "- (" + tmpR + ")";
    imag[3*1+0] += "- (" + tmpI + ")";
    // modify c[1][1]
    CoefXpr::ApplyTernaryFunc( tmpR, tmpI, c21R, c12R, c22R, 
                               c21I, c12I, c22I, OP_MULT, OP_DIV );
    real[3*1+1] += "- (" + tmpR + ")";
    imag[3*1+1] += "- (" + tmpI + ")";
    
    if( !isComplex_ ) {
      imag.Init("0.0");
    }
    
  } else if( tensorType_ == FULL ) {
    numRows = 6;
    numCols = 6;
    real = aR;
    imag = aI;
    if( !isComplex_ )
      imag.Init("0.0");
    
  } else {
    EXCEPTION( "Desired sub-tensor type not known" );
  }
  
  // In the end, put everything in brackets
  for( UInt i = 0; i < numRows*numCols; ++i ) {
    real[i] = Bracket(real[i]);
    imag[i] = Bracket(imag[i]);
  }
  
  // interchange variables, if transposed tensor is requested
  if( transposed_ ) {
    StdVector<std::string> tmp;
    CoefXpr::Tranpose( numRows, numCols, real, tmp );
    real = tmp;
    CoefXpr::Tranpose( numRows, numCols, imag, tmp );
    imag = tmp;
  }
 
}

void CoefXprMechSubTensor::GetArgs( std::map<std::string, PtrCoefFct > & vars ) const {
  vars[aName_] = a_;
}

// ==========================================================================
//  TENSOR REPRESENTATIONS (GENERAL MATERIAL)
// ==========================================================================
void CoefXprSubTensor::Init( PtrCoefFct a ) {
  // check dimensionality
  if( a->GetDimType() != CoefFunction::TENSOR ) {
    EXCEPTION( "Argument must be of type TENSOR" );
  }
  
  // ensure that tensor is either a 3x3 or a 6x3 tensor
  UInt numRowsA, numColsA;
  a->GetTensorSize(numRowsA, numColsA);
  if( !((numRowsA == 3 && numColsA == 6) || 
        (numRowsA == 3 && numColsA == 3 ) ) ) {
    EXCEPTION( "Tensor must have eiter dimension 3 x 3 or 3 x 6" );
  }
  
  dimType_ = CoefFunction::TENSOR;
  isAnalytical_ = a->IsAnalytic();
  isComplex_ = a->IsComplex();
  a_ = a;
  aName_ = CoefXpr::GetUniqueVarName();
  tensorType_ = NO_TENSOR;
  transposed_ = false;
}

CoefXprSubTensor::CoefXprSubTensor( PtrCoefFct a ) {
  Init(a);
}
CoefXprSubTensor::CoefXprSubTensor( const CoefXpr& a)  {
  
  Global::ComplexPart part = a.IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct temp = CoefFunction::Generate( part, a );  
}

void CoefXprSubTensor::SetSubTensorType(SubTensorType subType,
                                            bool transpose ) {
  tensorType_ = subType;
  transposed_ = transpose;
}


void CoefXprSubTensor::GetTensorXpr( UInt& numRows, UInt& numCols,
                                     StdVector<std::string>& real, 
                                     StdVector<std::string>& imag ) const {
  // ensure, that a subtype is set
  if( tensorType_ == NO_TENSOR ) {
    EXCEPTION( "No tensor sub type was set");
  }

  StdVector<std::string> aR, aI;
  UInt numRowsA, numColsA;

  if( isAnalytical_) {
    CoefFunctionAnalytic & coefA = 
        dynamic_cast<CoefFunctionAnalytic&>(*a_);
    coefA.GetStrTensor( numRowsA, numColsA, aR, aI );
  } else {
    CoefFunction::GenTensorCompNames(aR, aI, aName_, a_);
    a_->GetTensorSize(numRowsA, numColsA);
  }
  
  
  // switch depending on size
  if( numRowsA == 3 && numColsA == 3) {
    // --------------------
    //  3 x 3 TENSOR 
    // --------------------
    if( tensorType_ == FULL ) {
      numRows = 3;
      numCols = 3;
      real = aR;
      imag = aI;
      if( !isComplex_ )
        imag.Init("0.0");
    } else {
      // in all other cases, we just consider the 2x2 submatrix
      numRows = 2;
      numCols = 2;
      real.Resize(4);
      imag.Resize(4);
      imag.Init("0.0");
      real[2*0+0] = aR[3*0+0];
      real[2*0+1] = aR[3*0+1];
      real[2*1+0] = aR[3*1+0];
      real[2*1+1] = aR[3*1+1];
      if( isComplex_ ) {
        imag[2*0+0] = aI[3*0+0];
        imag[2*0+1] = aI[3*0+1];
        imag[2*1+0] = aI[3*1+0];
        imag[2*1+1] = aI[3*1+1];
      }
    }
    
    
  } else if( numRowsA == 3 && numColsA == 6 ) {
    // --------------------
    //  6 x 3 TENSOR 
    // --------------------
    if( tensorType_ == FULL ) {
      numRows = 6;
      numCols = 3;
      real = aR;
      imag = aI;
      if( !isComplex_ )
        imag.Init("0.0");
    } else if( tensorType_ == AXI ) {
      numRows = 2;
      numCols = 4;
      real.Resize(numRows * numCols );
      imag.Resize(numRows * numCols );
      imag.Init("0.0");
      real[4*0+0] = aR[6*0+0];
      real[4*0+1] = aR[6*0+1];
      real[4*0+2] = aR[6*0+5];
      real[4*0+3] = aR[6*0+2];
      real[4*1+0] = aR[6*1+0];
      real[4*1+1] = aR[6*1+1];
      real[4*1+2] = aR[6*1+5];
      real[4*1+3] = aR[6*1+2];
      if( isComplex_ ) {
        imag[4*0+0] = aI[6*0+0];
        imag[4*0+1] = aI[6*0+1];
        imag[4*0+2] = aI[6*0+5];
        imag[4*0+3] = aI[6*0+2];
        imag[4*1+0] = aI[6*1+0];
        imag[4*1+1] = aI[6*1+1];
        imag[4*1+2] = aI[6*1+5];
        imag[4*1+3] = aI[6*1+2];
      }
    } else if( tensorType_ == PLANE_STRAIN ||
        tensorType_ == PLANE_STRESS ) {
      numRows = 2;
      numCols = 3;
      real.Resize(numRows * numCols );
      imag.Resize(numRows * numCols );
      imag.Init("0.0");
      real[3*0+0] = aR[6*0+0];
      real[3*0+1] = aR[6*0+1];
      real[3*0+2] = aR[6*0+5];
      real[3*1+0] = aR[6*1+0];
      real[3*1+1] = aR[6*1+1];
      real[3*1+2] = aR[6*1+5];
      if( isComplex_ ) {
        imag[3*0+0] = aI[6*0+0];
        imag[3*0+1] = aI[6*0+1];
        imag[3*0+2] = aI[6*0+5];
        imag[3*1+0] = aI[6*1+0];
        imag[3*1+1] = aI[6*1+1];
        imag[3*1+2] = aI[6*1+5];
      }
    } // PLANE STRAIN
  }  else {
    EXCEPTION( "Unsupported tensor dimension");
  }// DIMENSIONALITY

  // In the end, put everything in brackets
  for( UInt i = 0; i < numRows*numCols; ++i ) {
    real[i] = Bracket(real[i]);
    imag[i] = Bracket(imag[i]);
  }
  
  // interchange variables, if transposed tensor is requested
  if( transposed_ ) {
    StdVector<std::string> tmp;
    CoefXpr::Tranpose( numRows, numCols, real, tmp );
    real = tmp;
    CoefXpr::Tranpose( numRows, numCols, imag, tmp );
    imag = tmp;
    std::swap(numRows, numCols);
  }
 
}

void CoefXprSubTensor::GetArgs( std::map<std::string, PtrCoefFct > & vars ) const {
  vars[aName_] = a_;
}

} // end of namespace
