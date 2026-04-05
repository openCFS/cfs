#include "CoefXpr.hh"
#include "Utils/ToolsFull.hh"
#include "boost/algorithm/string/erase.hpp"

namespace CoupledField {


// define local helper function for wrapping string expressions
inline std::string B(const std::string& xpr ) {
  return "(" + xpr + ")";
}

// helper function to determine if function is zero
const std::string zero1 = boost::lexical_cast<std::string>(0.0);
const std::string zero2 = "0.0";
const std::string zero3 = "0";
const std::string zero4 = B(zero1);
const std::string zero5 = B(zero2);
const std::string zero6 = B(zero3);

inline bool IsZero( const std::string& arg ) {
  return (arg == zero1 || arg == zero2 || arg == zero3 || 
          arg == zero4 || arg == zero5 || arg == zero6);
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
    case OP_SQRT_NEGATIVE:
    case OP_TRACE:
    case OP_INV:
    case OP_DET:
    case OP_RE:
    case OP_IM:
      return 1;
      break;
        
    // BINARY FUNCTIONS
    case OP_ADD:
    case OP_SUB:
    case OP_MULT:
    case OP_MULT_COMP:
    case OP_MULT_CONJ:
    case OP_MULT_VOIGT_TENSOR_VEC:
    case OP_MULT_VOIGT_TENSOR_VEC_CONJ:
    case OP_MULT_TENSOR:
    case OP_DIV:
    case OP_CROSS:
    case OP_CROSS_AXI:
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
    case OP_MULT_CONJ:
      return "*";
      break;
    case OP_DIV:
      return "/";
      break;
    case OP_MULT_TENSOR:
      return "\u2297";
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
        case OP_DET:
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
           case OP_MULT_COMP:
           case OP_CROSS_AXI:
             dim = CoefFunction::VECTOR;
             break;
           case OP_MULT_VOIGT_TENSOR_VEC:
           case OP_MULT_VOIGT_TENSOR_VEC_CONJ:
             // the first "tensor" is assumed to be symmetric, so we 
             // assume Voigt notation, so technically both arguments are "vectors"
             dim = CoefFunction::VECTOR;
             break;
           case OP_MULT:
           case OP_MULT_CONJ:
             // scalar product of two vectors
             dim = CoefFunction::SCALAR;
             break;
           case OP_DIV:
             EXCEPTION( "Division of two vectors not defined");
             break;
           case OP_MULT_TENSOR:
             dim = CoefFunction::TENSOR;
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
           case OP_MULT_CONJ:
           case OP_MULT_TENSOR:
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
    if( (a->GetDimType() == CoefFunction::VECTOR &&  
        b->GetDimType() == CoefFunction::SCALAR) ||
        (a->GetDimType() == CoefFunction::SCALAR &&  
         b->GetDimType() == CoefFunction::VECTOR)
        ) {
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
          << CoefFunction::coefDimType.ToString(a->GetDimType())
          << " and " << CoefFunction::coefDimType.ToString(b->GetDimType()) );
    }
    
  }
  
  return dim;
}

bool CoefXpr::IsComplex( PtrCoefFct a,
                         OpType op ) {
  if( (op == OP_NORM) || (op == OP_RE) || (op == OP_IM) ) {
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

void CoefXpr::Transpose( UInt nRows, UInt nCols, 
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
      
    case OP_SQRT_NEGATIVE:
      args = "-sqrt(", B(argReal), ")";
      retReal = B(args.Serialize(' '));
      break;

    case OP_INV:
      args = "1.0/(", B(argReal), ")";
      retReal = B(args.Serialize(' '));
      break;

    case OP_RE:
      args = B(argReal);
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
  StdVector<std::string> args_imag;
  switch( op ) {
    case OP_NORM:
      args = "sqrt(", B(argReal), "*", B(argReal), "+", B(argImag), "*", B(argImag), ")";
      retReal = B(args.Serialize(' '));
      retImag = "0.0";
      break;
    
    case OP_RE:
      args = B(argReal);
      retReal = B(args.Serialize(' '));
      retImag = "0.0";
      break;

    case OP_IM:
      args = B(argImag);
      retReal = B(args.Serialize(' '));
      retImag = "0.0";
      break;

    case OP_SQRT:
      // sqrt(z) = sqrt(a+i*b) = +-[ sqrt( (|z|+a)/2 )  +  i * b/|b| * sqrt( (|z|-a)/2 ) ]
      if (IsZero(argImag)) {
      args = "sqrt(", B(argReal), ")";
      args_imag = "0.0";
      }
      else {
        args = "sqrt((sqrt(", B(argReal), "*", B(argReal), "+", B(argImag), "*", B(argImag), ")+", B(argReal), ")/2)";
        args_imag = "", B(argImag), "/abs(", B(argImag), ")*sqrt((sqrt(", B(argReal), "*", B(argReal), "+", B(argImag), "*", B(argImag), ")-", B(argReal), ")/2)";
      }
      retReal = B(args.Serialize(' '));
      retImag = B(args_imag.Serialize(' '));
      break;
      
    case OP_SQRT_NEGATIVE:
      if (IsZero(argImag)) {
      args = "-sqrt(", B(argReal), ")";
      args_imag = "0.0";
      }
      else {
        args = "-sqrt((sqrt(", B(argReal), "*", B(argReal), "+", B(argImag), "*", B(argImag), ")+", B(argReal), ")/2)";
        args_imag = "-", B(argImag), "/abs(", B(argImag), ")*sqrt((sqrt(", B(argReal), "*", B(argReal), "+", B(argImag), "*", B(argImag), ")-", B(argReal), ")/2)";
      }
      retReal = B(args.Serialize(' '));
      retImag = B(args_imag.Serialize(' '));
      break;

    case OP_INV:
      // Simply apply binary function
      ApplyBinaryFunc( retReal, retImag, "1.0", argReal, "0.0" , argImag, OP_DIV );
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
    case OP_MULT_CONJ:
    case OP_DIV:
      if( IsZero(arg1Real) && IsZero(arg2Real) ) {
        retReal = zero2;
      } else {
        args =  B(arg1Real), OpToString(op), B(arg2Real);
        retReal = B(args.Serialize(' ') );
      }
      break;
    case OP_POW:
      args = B(arg1Real), "^", B(arg2Real);
      retReal = B(args.Serialize(' ' ));
      break;
    case OP_MULT_TENSOR:
      if( IsZero(arg1Real) && IsZero(arg2Real) ) {
        retReal = zero2;
      } else {
        args =  B(arg1Real), "*" , B(arg2Real);
        retReal = B(args.Serialize(' ') );
      }
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
  
  bool z1Real = IsZero(arg1Real);
  bool z2Real = IsZero(arg2Real);
  bool z1Imag = IsZero(arg1Imag);
  bool z2Imag = IsZero(arg2Imag);

  
  
  switch( op ) {
    case OP_ADD:
    case OP_SUB:
      if( z1Real && z2Real ) {
        retReal = zero2;
      } else {
        args = B(arg1Real), OpToString(op), B(arg2Real);
        retReal = B(args.Serialize(' ') );
      }
      if( z1Imag && z2Imag ) {
        retImag = zero2;
      } else {
        args = B(arg1Imag), OpToString(op), B(arg2Imag);
        retImag = B(args.Serialize(' '));
      }
      break;

    case OP_MULT:
    case OP_MULT_TENSOR:
      // if one of the two operands is zero, just take the shortcut
      if( ( z1Real && z1Imag ) || (z2Real && z2Imag) ) { 
        retReal = zero2;
        retImag = zero2;
      } else {
        if( z1Real || z2Real ) {
          args = "- (", B(arg1Imag), "*", B(arg2Imag), ")";
        } else if (z1Imag || z2Imag ) {
          args = "(", B(arg1Real), "*", B(arg2Real), ")";
        } else {
          args = "(", B(arg1Real), "*", B(arg2Real), ") - (", B(arg1Imag), "*", B(arg2Imag), ")";
        }
        retReal = B(args.Serialize(' '));
        args.Clear();
        if( z1Real || z2Imag ) {
          args = "(", B(arg1Imag), "*", B(arg2Real), ")";
        } else if ( z1Imag || z2Real ) {
          args = "(", B(arg1Real), "*", B(arg2Imag), ")";
        } else {
          args = "(", B(arg1Real), "*", B(arg2Imag), ") + (", B(arg1Imag), "*", B(arg2Real), ")";
        }
        retImag = B(args.Serialize(' '));
      }
      break;
    case OP_MULT_CONJ:
      // if one of the two operands is zero, just take the shortcut
      if( ( z1Real && z1Imag ) || (z2Real && z2Imag) ) { 
        retReal = zero2;
        retImag = zero2;
      } else {
        if( z1Real || z2Real ) {
          args = "(", B(arg1Imag), "*", B(arg2Imag), ")";
        } else if (z1Imag || z2Imag ) {
          args = "(", B(arg1Real), "*", B(arg2Real), ")";
        } else {
          args = "(", B(arg1Real), "*", B(arg2Real), ") + (", B(arg1Imag), "*", B(arg2Imag), ")";
        }
        retReal = B(args.Serialize(' '));
        args.Clear();
        if( z1Real || z2Imag ) {
          args = "(", B(arg1Imag), "*", B(arg2Real), ")";
        } else if ( z1Imag || z2Real ) {
          args = "-(", B(arg1Real), "*", B(arg2Imag), ")";
        } else {
          args = "-(", B(arg1Real), "*", B(arg2Imag), ") + (", B(arg1Imag), "*", B(arg2Real), ")";
        }
        retImag = B(args.Serialize(' '));
      }
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
   dependType_ = a->GetDependency();
   coordSys_ = a->GetCoordinateSystem();
   isComplex_ = CoefXpr::IsComplex(a, op );
   a_ = a;
   aName_ = CoefXpr::GetUniqueVarName();
   op_ = op;
}
  
CoefXprUnaryOp::CoefXprUnaryOp( MathParser * mp,
                                PtrCoefFct a, 
                                CoefXpr::OpType op ) : CoefXpr(mp) {
 Init(a, op);
}

CoefXprUnaryOp::CoefXprUnaryOp( MathParser * mp,
                                const std::string& a,
                                CoefXpr::OpType op ) : CoefXpr(mp) {
  
  PtrCoefFct temp = CoefFunction::Generate( mp_, Global::REAL, a );
  Init( temp, op );
}

CoefXprUnaryOp::CoefXprUnaryOp( MathParser * mp,
                                const CoefXpr& a,
                                CoefXpr::OpType op ) : CoefXpr(mp) {
  
  Global::ComplexPart part = a.IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct temp = CoefFunction::Generate( mp_, part, a );
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
        CoefXpr::ApplyBinaryFunc( tmp, aR[0], aR[0], OP_MULT );
        args.Push_back(tmp);
        for( UInt i = 1; i < aR.GetSize(); ++ i ) {
          CoefXpr::ApplyBinaryFunc( tmp, aR[i], aR[i], OP_MULT );
          args.Push_back("+");
          args.Push_back(tmp);
        }
        args.Push_back(") )");
        imag = "0.0";
        real = args.Serialize(' ');
      } else {
        args = "( sqrt(";
        // Although we create temporary return variables for real- and imaginary part, 
        // we only get a real-part in case of multiplication with conjugate transposed.
        std::string tmpR, tmpI;
        CoefXpr::ApplyBinaryFunc( tmpR, tmpI,  aR[0], aR[0], aI[0], aI[0], OP_MULT_CONJ );
        args.Push_back(tmpR);
        for( UInt i = 1; i < aR.GetSize(); ++ i ) {
          CoefXpr::ApplyBinaryFunc( tmpR, tmpI, aR[i], aR[i], aI[i], aI[i], OP_MULT_CONJ );
          args.Push_back("+");
          args.Push_back(tmpR);
        }
        args.Push_back(") )");
        imag = "0.0";
        real = args.Serialize(' ');
        }
      } else {
        EXCEPTION( "Unary operation of type " << op_ 
                   << " not supported for vector type functions!" );
      }
    }
  // -------------
  //  TENSOR CASE
  // -------------
  else if( a_->GetDimType() == CoefFunction::TENSOR ) {
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
    if( op_ == OP_DET ) { 
      if( numRowsA != numColsA)  {
        EXCEPTION( "Can not calculate determinant of tensor of size "
            << numRowsA << "x"<< numColsA );
      }
      if( numRowsA == 1 ) {
        real = B(aR[0]);
        imag = B(aI[0]);
      } else if( numRowsA == 2 ) {
        std::string tmpR, tmpI;
        if( !isComplex_ ) {
          ApplyBinaryFunc( tmpR, aR[0], aR[3], OP_MULT );
          real = tmpR + "-";
          ApplyBinaryFunc( tmpR, aR[1], aR[2], OP_MULT );
          real += tmpR;
          real = B(real);
        } else {
          ApplyBinaryFunc( tmpR, tmpI, aR[0], aI[0], aR[3], aI[3], OP_MULT );
          real = tmpR + "-";
          imag = tmpI + "-";
          ApplyBinaryFunc( tmpR, tmpI, aR[1], aI[1], aR[2], aI[2], OP_MULT );
          real += tmpR;
          imag += tmpI;
          real = B(real);
          imag = B(imag);
        }
      } else if( numRowsA == 3 ) {
        std::string tmpR, tmpI;
        if( !isComplex_ ) {
          ApplyTernaryFunc( tmpR, aR[0], aR[4], aR[8], OP_MULT, OP_MULT);
          real = tmpR;
          ApplyTernaryFunc( tmpR, aR[1], aR[5], aR[6], OP_MULT, OP_MULT);
          real += " + " + tmpR;
          ApplyTernaryFunc( tmpR, aR[2], aR[3], aR[7], OP_MULT, OP_MULT);
          real += " + " + tmpR;

          ApplyTernaryFunc( tmpR, aR[2], aR[4], aR[6], OP_MULT, OP_MULT);
          real += " - " + tmpR;
          ApplyTernaryFunc( tmpR, aR[1], aR[3], aR[8], OP_MULT, OP_MULT);
          real += " - " + tmpR;
          ApplyTernaryFunc( tmpR, aR[0], aR[5], aR[7], OP_MULT, OP_MULT);
          real += " - " + tmpR;
          real = B(real);
        } else {
          ApplyTernaryFunc( tmpR, tmpI, aR[0], aI[0], aR[4], aI[4],
                            aR[8], aI[8], OP_MULT, OP_MULT);
          real = tmpR;
          imag = tmpI;

          ApplyTernaryFunc( tmpR, tmpI, aR[1], aI[1], aR[5], aI[5],
                            aR[6], aI[6], OP_MULT, OP_MULT);
          real += " + " + tmpR;
          imag += " + " + tmpI;

          ApplyTernaryFunc( tmpR, tmpI, aR[2], aI[2], aR[3], aI[3],
                            aR[7], aI[7], OP_MULT, OP_MULT);
          real += " + " + tmpR;
          imag += " + " + tmpI;
          ApplyTernaryFunc( tmpR, tmpI, aR[2], aI[2], aR[4], aI[4],
                            aR[6], aI[6], OP_MULT, OP_MULT);
          real += " - " + tmpR;
          imag += " - " + tmpI;

          ApplyTernaryFunc( tmpR, tmpI, aR[1], aI[1], aR[3], aI[3],
                            aR[8], aI[8], OP_MULT, OP_MULT);

          real += " - " + tmpR;
          imag += " - " + tmpI;
          ApplyTernaryFunc( tmpR, tmpI, aR[0], aI[0], aR[5], aI[5],
                            aR[7], aI[7], OP_MULT, OP_MULT);
          real += " - " + tmpR;
          imag += " - " + tmpI;
          real = B(real);
          imag = B(imag);
        }
      } else {
        EXCEPTION( "General tensor inversion only implemented up to size 3x3");
      }

    } else {
      EXCEPTION( "Unary operation of type " << op_ 
                 << " not supported for matrix type functions!" );
    }
  }// if tensor
}

void CoefXprUnaryOp::GetVectorXpr( StdVector<std::string>& real, 
                                 StdVector<std::string>& imag ) const {
  

  // Switch depending on type of argument
  if( a_->GetDimType() == CoefFunction::SCALAR ) {
    EXCEPTION("CoefXprUnaryOp::GetVectorXpr SCALAR case not implemented!");
  }else if( a_->GetDimType() == CoefFunction::VECTOR ) {
    StdVector<std::string> aR, aI;
    UInt sizeA;

    if( isAnalytical_) {
      EXCEPTION("CoefXprUnaryOp::GetVectorXpr No vector valued analytic expression implemented!");
    } else {
      CoefFunction::GenVecCompNames(aR, aI, aName_, a_);
      sizeA = a_->GetVecSize();
      if( op_ == OP_RE) {
        real.Resize(sizeA);
        imag.Resize(sizeA);
        for(UInt i = 0; i < sizeA; ++i){
          ApplyUnaryFunc( real[i], imag[i], aR[i], aI[i], OP_RE );
        }
      } else if ( op_ == OP_IM ) {
        imag.Resize(sizeA);
        real.Resize(sizeA);
        for(UInt i = 0; i < sizeA; ++i){
          ApplyUnaryFunc( real[i], imag[i], aR[i], aI[i], OP_IM );
        }
      } else {
        EXCEPTION("CoefXprUnaryOp::GetVectorXpr This vector valued expression is not implemented yet!");
      }
    }
  }else if( a_->GetDimType() == CoefFunction::TENSOR ) {
    EXCEPTION("CoefXprUnaryOp::GetVectorXpr TENSOR case not implemented!");
  }
}

void CoefXprUnaryOp::GetTensorXpr( UInt& numRows, UInt& numCols,
                                   StdVector<std::string>& real, 
                                   StdVector<std::string>& imag ) const {
  assert((this->dimType_ == CoefFunction::NO_DIM) || 
         (this->dimType_ == CoefFunction::TENSOR) );

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
  if( op_ == OP_INV) {
    // Choose algorithm depending on size of "tensor"
    if( numRowsA != numColsA)  {
      EXCEPTION( "Tensor is of size " << numRowsA << "x"
                 << numColsA << " and can not be inverted! ");
    }
    // Obtain Jacobian detmerinant
    std::string detR, detI;
    CoefXprUnaryOp detOp = CoefXprUnaryOp( mp_, a_, OP_DET );
    detOp.GetScalarXpr( detR, detI );


    if( numRowsA == 1 ) {
      real.Resize(1);
      imag.Resize(1);
      ApplyUnaryFunc( real[0], imag[0], aR[0], aI[0], OP_INV );

    } else if( numRowsA == 2 ) {
      EXCEPTION( "Inversion of 2x2 matrix not fully implemented yet");
      //      real.Resize(4);
      //      std::string det;
      //      if( !isComplex_ ) {
      //        real[0] = aR[0];
      //        real[3] = aR[3];
      //        real[1] = B("-"+aR[1]);
      //        real[2] = B("-"+aR[2]);
      //      } else {
      //        imag.Resize(4);
      //      }
      //      

    } else if( numRowsA == 3 ) {
      EXCEPTION( "Inversion of 3x3 matrix not fully implemented yet");
    } else {
      EXCEPTION( "General tensor inversion only implemented up to size 3x3");
    }
  } else {
    EXCEPTION( "Unary tensor operation " << op_ << " not implemented" );
  }
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
  
   if( GetNumOperands(op) != 2  ) {
     EXCEPTION( "Operand must have exactly 2 operands" );
   }
   
   dimType_ = GetDimType( a, b, op);
   isAnalytical_ = a->IsAnalytic() && b->IsAnalytic();
   isComplex_ = a->IsComplex() || b->IsComplex();
   dependType_ = std::max(a->GetDependency(), 
                          b->GetDependency());
   
   // check for distinct coordinate systems
   CoordSystem* aCoord = a->GetCoordinateSystem();
   CoordSystem* bCoord = b->GetCoordinateSystem();
   if( aCoord != NULL || bCoord != NULL) {
     if (aCoord != NULL && bCoord == NULL) {
       coordSys_ = aCoord;
     } else if (aCoord == NULL && bCoord != NULL) {
       coordSys_ = bCoord;
     } else if (aCoord == bCoord) {
       coordSys_ = aCoord;
     } else {
       EXCEPTION("Case of two distinct coordinate systemw within "\
                 "one expression is not tested yet");
     }
   }
   
   a_ = a;
   b_ = b;
   aName_ = CoefXpr::GetUniqueVarName();
   bName_ = CoefXpr::GetUniqueVarName();
   
   op_ = op;
}

CoefXprBinOp::CoefXprBinOp( MathParser * mp,
                            PtrCoefFct a, 
                            PtrCoefFct b,
                            CoefXpr::OpType op ) : CoefXpr(mp) {
  Init( a, b, op );
}

CoefXprBinOp::CoefXprBinOp( MathParser * mp,
                            PtrCoefFct a, 
                            const CoefXpr& b,
                            CoefXpr::OpType op ) : CoefXpr(mp) {
  
  Global::ComplexPart part = b.IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct temp  = CoefFunction::Generate( mp_, part, b );
  Init( a, temp, op );
}

CoefXprBinOp::CoefXprBinOp( MathParser * mp,
                            const CoefXpr& a,
                            PtrCoefFct b, 
                            CoefXpr::OpType op ) : CoefXpr(mp) {
  
  Global::ComplexPart part = a.IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct temp  = CoefFunction::Generate( mp_, part, a );
  Init( temp, b, op );
}


CoefXprBinOp::CoefXprBinOp( MathParser * mp,
                            PtrCoefFct a, 
                            const std::string& b,
                            CoefXpr::OpType op ) : CoefXpr(mp) {
  
  Global::ComplexPart part = a->IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct temp = CoefFunction::Generate( mp_, part, b );
  Init( a, temp, op );
}

CoefXprBinOp::CoefXprBinOp( MathParser * mp,
                            const std::string& a, 
                            PtrCoefFct b, 
                            CoefXpr::OpType op ) : CoefXpr(mp) {
  
  Global::ComplexPart part = b->IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct temp = CoefFunction::Generate( mp_, part, a );
  Init( temp, b, op );
  
}

CoefXprBinOp::CoefXprBinOp( MathParser * mp,
                            const std::string& a, 
                            const CoefXpr& b, 
                            CoefXpr::OpType op ) : CoefXpr(mp) {
  Global::ComplexPart part = b.IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct tempA  = CoefFunction::Generate( mp_, part, a );
  PtrCoefFct tempB  = CoefFunction::Generate( mp_, part, b );
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
    if( op_ == OP_MULT || op_ == OP_MULT_CONJ ) {
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
        CoefXpr::ApplyBinaryFunc( real, aR[0], bR[0], op_ );
        for( UInt i = 1; i < numEntries; ++ i ) {
          CoefXpr::ApplyBinaryFunc( temp, aR[i], bR[i], op_ );
          real += "+" + temp;
        }
        imag = "0.0";
      } else {
        std::string tempReal, tempImag;
        // treat first entry
        CoefXpr::ApplyBinaryFunc( real, imag, aR[0], bR[0], 
                                  aI[0], bI[0], op_ );
        for( UInt i = 1; i < numEntries; ++ i ) {
          CoefXpr::ApplyBinaryFunc( tempReal, tempImag , 
                                    aR[i], bR[i],aI[i], bI[i], op_ );
          real += "+" + tempReal;
          imag += "+" + tempImag;
        }
      }
      // put brackets around expression
      real = Bracket(real);
      imag = Bracket(imag);
        
    } else {
      EXCEPTION("The only allowed binary vector-vector operation "
                << "is OP_MULT and OP_MULT_CONJ" ); 
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
  // *

  if( a_->GetDimType() == CoefFunction::VECTOR  &&
      b_->GetDimType() == CoefFunction::VECTOR ) {
    // --------------------
    //  VECTOR-VECTOR CASE
    // --------------------
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
    if( op_ == OP_ADD || op_ == OP_SUB || op_ == OP_MULT_COMP ) {
      CoefXpr::OpType op;
      if (op_ == OP_ADD || op_ == OP_SUB) {
    	  op = op_;
      } else {
    	  op = OP_MULT;
      }
      if( !isComplex_ ) {
        for( UInt i = 0; i < numEntries; ++ i ) {
          CoefXpr::ApplyBinaryFunc( real[i], aR[i], bR[i], op);
          imag[i] = "0.0";
        }
      } else {
        for( UInt i = 0; i < numEntries; ++ i ) {
          CoefXpr::ApplyBinaryFunc( real[i], imag[i], aR[i], bR[i],
                                    aI[i], bI[i], op );
        }
      }
    }
    
    // === OP_CROSS ===
    else if( op_ == OP_CROSS ) {
      
      // switch, depending on dimension
      if( numEntries == 1 ) {
        // 1) reduced cross product in 2D
        if(!isComplex_) {
          real.Resize(2);
          // 1st row
          std::string temp;
          CoefXpr::ApplyBinaryFunc( temp, aR[0], bR[1], OP_MULT );
          real[0] = "-" + B(temp);
          // 2nd row
          CoefXpr::ApplyBinaryFunc( temp, aR[0], bR[0], OP_MULT );
          real[1] = temp;
        } else {
          real.Resize(2);
          imag.Resize(2);
          std::string tempR, tempI;
          // 1st row
          CoefXpr::ApplyBinaryFunc( tempR, tempI,  
                                    aR[0], bR[1],
                                    aI[0], bI[1], OP_MULT );
          real[0] = "-" + B(tempR);
          imag[0] = "-" + B(tempI);
          // 2nd row
          CoefXpr::ApplyBinaryFunc( tempR, tempI,  
                                    aR[0], bR[0],
                                    aI[0], bI[0], OP_MULT );
          real[1] = tempR;
          imag[1] = tempI;
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
         
       } else { 
         EXCEPTION( "Vector / cross product only defined in 1D and 3D");
       }
      // === OP_CROSS_AXI ===
    } else if( op_ == OP_CROSS_AXI ) {

      // switch, depending on dimension
      if( numEntries == 1 ) {
        // this case occurs e.g. in the 2D / axisymmetric case, when
        // the first operand is pointing into the third direction e.g.
        //
        //     (0  )    (B_x)   (-J_z * B_y)
        // F = (0  ) x  (B_y) = ( J_z * B_x)
        //     (J_z)    ( 0 )   (     0    )
        // 
        // Here, H is just given as a scalar coefficient function
        // 1) reduced cross product
        if(!isComplex_) {
          real.Resize(2);
          // 1st row
          std::string temp;
          CoefXpr::ApplyBinaryFunc( temp, aR[0], bR[1], OP_MULT );
          real[0] = temp;
          // 2nd row
          CoefXpr::ApplyBinaryFunc( temp, aR[0], bR[0], OP_MULT );
          real[1] = "-"+B(temp);
        } else {
          real.Resize(2);
          imag.Resize(2);
          std::string tempR, tempI;
          // 1st row
          CoefXpr::ApplyBinaryFunc( tempR, tempI,  
                                    aR[0], bR[1],
                                    aI[0], bI[1], OP_MULT );
          real[0] = tempR;
          imag[0] = tempI;
          // 2nd row
          CoefXpr::ApplyBinaryFunc( tempR, tempI,  
                                    aR[0], bR[0],
                                    aI[0], bI[0], OP_MULT );
          real[1] = "-"+B(tempR);
          imag[1] = "-"+B(tempI);
        }
      } else {
        EXCEPTION("Wrong dimension");
      }
     
    }
    // === OP_MULT_TENSOR_VEC
    else if( op_ == OP_MULT_VOIGT_TENSOR_VEC ||
             op_ == OP_MULT_VOIGT_TENSOR_VEC_CONJ ) {
      
      // Here we assume that the first "vector" represents a symmetric
      // tensor in Voigt Notation and the second argument is a vector.
      // We allow three cases:
      // a) FULL: (3x3) tensor * vector with 3 components
      // b) PLANE: (2x2) tensor * vector with 2 components
      // c) AXI: (3 x 3 )tensor * vector with 2 components. 
      //          Attention: Here only the upper 2x2 block is considered
      //          for multiplication (circumferential components can
      //          not be considered)
      
      // Attention: We have to re-size the result to the number
      //            of vector components
      real.Resize( bR.GetSize() );
      imag.Resize( bR.GetSize() );
      
      CoefXpr::OpType scalOp;
      if( op_ == OP_MULT_VOIGT_TENSOR_VEC ) {
        scalOp = OP_MULT;
      } else {
        scalOp = OP_MULT_CONJ;
      }
      
      if( aR.GetSize() == 6 && bR.GetSize() == 3 ) {
        // ----------------
        //  3D / FULL case
        // ----------------
        // tensor [t] = (t_xx, t_yy, t_zz, t_yz, t_xz, t_xy)
        // 1st row
        if(!isComplex_) {
          std::string temp;
          // 1st row
          CoefXpr::ApplyBinaryFunc( temp, aR[0], bR[0], scalOp );
          real[0] = temp + " + ";
          CoefXpr::ApplyBinaryFunc( temp, aR[5], bR[1], scalOp );
          real[0] += temp + " + ";
          CoefXpr::ApplyBinaryFunc( temp, aR[4], bR[2], scalOp );
          real[0] += temp;
          // 2nd row
          CoefXpr::ApplyBinaryFunc( temp, aR[5], bR[0], scalOp );
          real[1] = temp + " + ";
          CoefXpr::ApplyBinaryFunc( temp, aR[1], bR[1], scalOp );
          real[1] += temp + " + ";
          CoefXpr::ApplyBinaryFunc( temp, aR[3], bR[2], scalOp );
          real[1] += temp;
          // 3rd row
          CoefXpr::ApplyBinaryFunc( temp, aR[4], bR[0], scalOp );
          real[2] = temp + " + ";
          CoefXpr::ApplyBinaryFunc( temp, aR[3], bR[1], scalOp );
          real[2] += temp + " + ";
          CoefXpr::ApplyBinaryFunc( temp, aR[2], bR[2], scalOp );
          real[2] += temp;
        } else {
          std::string tempR, tempI;
          // 1st row
          CoefXpr::ApplyBinaryFunc( tempR, tempI, 
                                    aR[0], bR[0], 
                                    aI[0], bI[0], scalOp );
          real[0] = tempR + " + ";
          imag[0] = tempI + " + ";
          CoefXpr::ApplyBinaryFunc( tempR, tempI, 
                                    aR[5], bR[1], 
                                    aI[5], bI[1], scalOp );
          real[0] += tempR + " + ";
          imag[0] += tempI + " + ";
          CoefXpr::ApplyBinaryFunc( tempR, tempI, 
                                    aR[4], bR[2], 
                                    aI[4], bI[2], scalOp );
          real[0] += tempR;
          imag[0] += tempI;
          // 2nd row
          CoefXpr::ApplyBinaryFunc( tempR, tempI, 
                                    aR[5], bR[0], 
                                    aI[5], bI[0], scalOp );
          real[1] = tempR + " + ";
          imag[1] = tempI + " + ";
          CoefXpr::ApplyBinaryFunc( tempR, tempI, 
                                    aR[1], bR[1], 
                                    aI[1], bI[1], scalOp );
          real[1] += tempR + " + ";
          imag[1] += tempI + " + ";
          CoefXpr::ApplyBinaryFunc( tempR, tempI, 
                                    aR[3], bR[2], 
                                    aI[3], bI[2], scalOp );
          real[1] += tempR;
          imag[1] += tempI;
          // 3rd row
          CoefXpr::ApplyBinaryFunc( tempR, tempI, 
                                    aR[4], bR[0], 
                                    aI[4], bI[0], scalOp );
          real[2] = tempR + " + ";
          imag[2] = tempI + " + ";
          CoefXpr::ApplyBinaryFunc( tempR, tempI, 
                                    aR[3], bR[1], 
                                    aI[3], bI[1], scalOp );
          real[2] += tempR + " + ";
          imag[2] += tempI + " + ";
          CoefXpr::ApplyBinaryFunc( tempR, tempI, 
                                    aR[2], bR[2], 
                                    aI[2], bI[2], scalOp );
          real[2] += tempR;
          imag[2] += tempI;
        }
      } else if ( (aR.GetSize() == 3 && bR.GetSize() == 2) ||  
                  (aR.GetSize() == 4 && bR.GetSize() == 2) ) {
        // ----------------
        //  2D / AXI case
        // ----------------
        // PLANE case: tensor [t] = (t_xx, t_yy, t_xy) 
        // AXI case: tensor [t] = (t_rr, t_zz, t_rz, r_phiphi)
        // 1st row
        if(!isComplex_) {
          std::string temp;
          // 1st row
          CoefXpr::ApplyBinaryFunc( temp, aR[0], bR[0], scalOp );
          real[0] = temp + " + ";
          CoefXpr::ApplyBinaryFunc( temp, aR[2], bR[1], scalOp );
          real[0] += temp;
          // 2nd row
          CoefXpr::ApplyBinaryFunc( temp, aR[2], bR[0], scalOp );
          real[1] += temp + " + ";
          CoefXpr::ApplyBinaryFunc( temp, aR[1], bR[1], scalOp );
          real[1] += temp;
        } else {
          std::string tempR, tempI;
          // 1st row
          CoefXpr::ApplyBinaryFunc( tempR, tempI, 
                                    aR[0], bR[0], 
                                    aI[0], bI[0], scalOp );
          real[0] = tempR + " + ";
          imag[0] = tempI + " + ";
          CoefXpr::ApplyBinaryFunc( tempR, tempI, 
                                    aR[2], bR[1], 
                                    aI[2], bI[1], scalOp );
          real[0] += tempR;
          imag[0] += tempI;
          // 2nd row
          CoefXpr::ApplyBinaryFunc( tempR, tempI, 
                                    aR[2], bR[0], 
                                    aI[2], bI[0], scalOp );
          real[1] = tempR + " + ";
          imag[1] = tempI + " + ";
          CoefXpr::ApplyBinaryFunc( tempR, tempI, 
                                    aR[1], bR[1], 
                                    aI[1], bI[1], scalOp );
          real[1] += tempR;
          imag[1] += tempI;
        }
      } else {
        EXCEPTION( "Dimension of arguments wrong." )
      }
      
      
    }
    else {
      EXCEPTION( "The only allowed (vector,vector)->vector operations "
                << "are OP_ADD, OP_SUB, OP_MULT_COMP and OP_CROSS" );
    }
  } else if( a_->GetDimType() == CoefFunction::SCALAR  &&
              b_->GetDimType() == CoefFunction::VECTOR ) {
    // --------------------
    //  SCALAR-VECTOR CASE
    // --------------------
      std::string aR, aI; 
      StdVector<std::string>  bR, bI;
       if( isAnalytical_) {
         CoefFunctionAnalytic & coefA = 
             dynamic_cast<CoefFunctionAnalytic&>(*a_);
         CoefFunctionAnalytic & coefB = 
             dynamic_cast<CoefFunctionAnalytic&>(*b_);
         coefA.GetStrScalar( aR, aI );
         coefB.GetStrVector( bR, bI );
       } else {
         CoefFunction::GenScalCompNames(aR, aI, aName_, a_);
         CoefFunction::GenVecCompNames(bR, bI, bName_, b_);
       }
       UInt numEntries = bR.GetSize();
       real.Resize( numEntries );
       imag.Resize( numEntries );

       // switch depending on operation type
       // === OP_ADD / OP_SUB ===
       if( op_ == OP_ADD || op_ == OP_SUB ||  
           op_ == OP_SUB || op_ == OP_MULT || 
           op_ == OP_MULT_CONJ ) {
         if( !isComplex_ ) {
           for( UInt i = 0; i < numEntries; ++ i ) {
             CoefXpr::ApplyBinaryFunc( real[i], aR, bR[i], op_ );
             imag[i] = "0.0";
           }
         } else {
           for( UInt i = 0; i < numEntries; ++ i ) {
             CoefXpr::ApplyBinaryFunc( real[i], imag[i], aR, bR[i],
                                       aI, bI[i], op_ );
           }
         }
    } else {
      EXCEPTION( "The only allowed (scalar,vector)->vector operations "
          << "are OP_ADD, OP_SUB, OP_MULT, OP_MULT_CONJ and OP_DIV" );
    }
  }
  else if( a_->GetDimType() == CoefFunction::TENSOR  &&
          b_->GetDimType() == CoefFunction::VECTOR ) {
      // --------------------
      //  TENSOR-VECTOR CASE
      // --------------------
      if ( op_ == OP_MULT ) { // compute multiplication
          // check size
          UInt ai,aj;
          a_->GetTensorSize(ai,aj);
          UInt bi = b_->GetVecSize();
          if (aj != bi) {
              EXCEPTION( "Size must be compatible" )
          }
          // compute
          StdVector<std::string> aR, aI;
          StdVector<std::string> bR, bI;

          UInt numRows, numCols; // size of a

          if( isAnalytical_) {
              CoefFunctionAnalytic & coefA =
                      dynamic_cast<CoefFunctionAnalytic&>(*a_);
              CoefFunctionAnalytic & coefB =
                      dynamic_cast<CoefFunctionAnalytic&>(*b_);
              coefA.GetStrTensor( numRows, numCols, aR, aI );
              coefB.GetStrVector( bR, bI );
          } else {
              CoefFunction::GenTensorCompNames(aR, aI, aName_, a_);
              CoefFunction::GenVecCompNames(bR, bI, bName_, b_);
              a_->GetTensorSize(numRows, numCols);
          }
          // output size
          real.Resize( numRows );
          imag.Resize( numRows );

          if( !isComplex_ ) {
              for( UInt i = 0; i < numRows; ++ i ) {
                  std::string sum;
                  CoefXpr::ApplyBinaryFunc(sum,aR[i*numCols],bR[0],CoefXpr::OP_MULT);
                  for (UInt j = 1; j < numCols; ++ j) {
                      CoefXpr::ApplyTernaryFunc( sum, sum, aR[i*numCols+j], bR[j], CoefXpr::OP_ADD, CoefXpr::OP_MULT);
                  }
                  real[i] = sum;
                  imag[i] = "0.0";
              }
          }
          else {
              for( UInt i = 0; i < numRows; ++ i ) {
                  std::string sumR;
                  std::string sumI;
                  CoefXpr::ApplyBinaryFunc(sumR,sumI,aR[i*numCols],bR[0],aI[i*numCols],bI[0],CoefXpr::OP_MULT);
                  for (UInt j = 1; j < numCols; ++ j) {
                      CoefXpr::ApplyTernaryFunc( sumR, sumI, sumR, aR[i*numCols+j], bR[j], sumI, aI[i*numCols+j], bI[j], CoefXpr::OP_ADD, CoefXpr::OP_MULT);
                  }
                  real[i] = sumR;
                  imag[i] = sumI;
              }
          }
      } else {
          EXCEPTION("only tesor x vector possible")
      }
  }
  else {
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
     } else  if( op_ == OP_MULT || op_ == OP_MULT_CONJ ) {
       
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
             CoefXpr::ApplyBinaryFunc( temp1, aR[posA], bR[posB], op_ );
             for( UInt k = 1; k < numColsA; ++k ) {
               posA ++;
               posB += numCols;
               CoefXpr::ApplyBinaryFunc( temp2, aR[posA], bR[posB], op_ );
               if( !IsZero(temp2)) {
                 temp1 += " + " + temp2;
               }
             } // numColsA
             if( IsZero(temp1)) {
               real[posRet++] = zero2;
             } else {
               real[posRet++] = Bracket(temp1);
             }
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
                                       aI[posA], bI[posB], op_ );
             for( UInt k = 1; k < numColsA; ++k ) {
               posA ++;
               posB += numCols;
               CoefXpr::ApplyBinaryFunc( temp2R, temp2I, aR[posA], bR[posB], 
                                         aI[posA], bI[posB], op_ );
               if( !IsZero(temp2R)) {
                 temp1R += " + " + temp2R;
               }
               if( !IsZero(temp2I)) {
                 temp1I += " + " + temp2I;
               }
             } // numColsA
             if( IsZero(temp1R) ) {
               real[posRet] = zero2;
             } else {
               real[posRet] = Bracket(temp1R);
             }
             if( IsZero(temp1I) ) {
               imag[posRet] = zero2;
             } else {
               imag[posRet] = Bracket(temp1I);  
             }
             
             posRet++;
           } // numCols
         } // numRows
       }
       
       
     } else if( op_ == OP_MULT_TENSOR ){
       numRows = numRowsA * numRowsB;
       numCols = numColsA * numColsB;
       real.Resize( numRows * numCols );
       imag.Resize( numRows * numCols );

       if(!isComplex_){
         std::string temp1;
         for(UInt iA = 0; iA < aR.GetSize(); iA++){
           std::string aReal = aR[iA];
           for(UInt iB = 0; iB < bR.GetSize(); iB++){
             temp1 = "";
             UInt rowA = iA / numColsA;
             UInt colA = iA % numColsA;
             UInt rowB = iB / numColsB;
             UInt colB = iB % numColsB;

             UInt targetRow = rowA*numRowsB+rowB;
             UInt targetCol = colA*numColsB+colB;

             CoefXpr::ApplyBinaryFunc( temp1, aReal, bR[iB], op_ );
             if( IsZero(temp1)) {
               real[targetRow*numRows+targetCol] = zero2;
             } else {
               real[targetRow*numRows+targetCol] = Bracket(temp1);
             }
           }
         }
       }else{
         std::string temp1,temp2,temp3,temp4;
         for(UInt iA = 0; iA < aR.GetSize(); iA++){
           std::string aReal = aR[iA];
           std::string aImag = aI[iA];
           for(UInt iB = 0; iB < bR.GetSize(); iB++){
             temp1 = "";
             UInt rowA = iA / numColsA;
             UInt colA = iA % numColsA;
             UInt rowB = iB / numColsB;
             UInt colB = iB % numColsB;

             UInt targetRow = rowA*numRowsB+rowB;
             UInt targetCol = colA*numColsB+colB;

             CoefXpr::ApplyBinaryFunc( temp1, aReal, bR[iB], op_ );
             CoefXpr::ApplyBinaryFunc( temp2, aImag, bI[iB], op_ );
             CoefXpr::ApplyBinaryFunc( temp3, aReal, bI[iB], op_ );
             CoefXpr::ApplyBinaryFunc( temp4, aImag, bR[iB], op_ );
             if( IsZero(temp1) && IsZero(temp2) ) {
               real[targetRow*numRows+targetCol] = zero2;
             } else {
               real[targetRow*numRows+targetCol] = Bracket(temp1+"-"+temp2);
             }
             if(IsZero(temp3) && IsZero(temp4) ){
               imag[targetRow*numRows+targetCol] = zero2;
             }else{
               imag[targetRow*numRows+targetCol] = Bracket(temp3+"+"+temp4);
             }
           }
         }
      }
     }else {
       EXCEPTION( "The only allowed (tensor,tensor)->tensor operations "
                << "are OP_ADD, OP_SUB, OP_MULT, OP_MULT_CONJ and OP_MULT_TENSOR");
     }
  } else if( a_->GetDimType() == CoefFunction::VECTOR  &&
             b_->GetDimType() == CoefFunction::VECTOR ) {
    StdVector<std::string> aR, aI, bR, bI;
    UInt sizeA, sizeB;
    if( isAnalytical_) {
      CoefFunctionAnalytic & coefA =
          dynamic_cast<CoefFunctionAnalytic&>(*a_);
      CoefFunctionAnalytic & coefB =
          dynamic_cast<CoefFunctionAnalytic&>(*b_);
      coefA.GetStrVector( aR, aI );
      coefB.GetStrVector( bR, bI );
      sizeA = aR.GetSize();
      sizeB = bR.GetSize();
    } else {
      sizeA = a_->GetVecSize();
      sizeB = b_->GetVecSize();
    }
    numRows = sizeA;
    numCols = sizeB;
    real.Resize( numRows * numCols );
    imag.Resize( numRows * numCols );

    if(!isComplex_){
      std::string temp1;
      for(UInt iA = 0; iA < aR.GetSize(); iA++){
        std::string aReal = aR[iA];
        for(UInt iB = 0; iB < bR.GetSize(); iB++){
          temp1 = "";
          CoefXpr::ApplyBinaryFunc( temp1, aReal, bR[iB], op_ );
          if( IsZero(temp1)) {
            real[iA*numRows+iB] = zero2;
          } else {
            real[iA*numRows+iB] = Bracket(temp1);
          }
        }
      }
    }else{
      EXCEPTION("Dyadic multiplication not available for complex valued vectors.")
    }
  }else {

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
  dependType_ = std::max(a->GetDependency(), 
                         b->GetDependency());
  
  // check for distinct coordinate systems
  CoordSystem* aCoord = a->GetCoordinateSystem();
  CoordSystem* bCoord = b->GetCoordinateSystem();
  if( aCoord != NULL || bCoord != NULL) {
    if (aCoord != NULL && bCoord == NULL) {
      coordSys_ = aCoord;
    } else if (aCoord == NULL && bCoord != NULL) {
      coordSys_ = bCoord;
    } else if (aCoord == bCoord) {
      coordSys_ = aCoord;
    } else {
      EXCEPTION("Case of two distinct coordinate systems within "\
                "one expression is not tested yet");
    }
  }
  
  isComplex_ = a->IsComplex() || b->IsComplex();
  a_ = a;
  b_ = b;
  aName_ = CoefXpr::GetUniqueVarName();
  bName_ = CoefXpr::GetUniqueVarName();
  op_ = op;
}

CoefXprVecScalOp::CoefXprVecScalOp( MathParser * mp,
                                    PtrCoefFct a, 
                                    PtrCoefFct b,
                                    CoefXpr::OpType op ) : CoefXpr(mp) {
  Init( a, b, op );
}

CoefXprVecScalOp::CoefXprVecScalOp( MathParser * mp,
                                    PtrCoefFct a, 
                                    const std::string& b,
                                    CoefXpr::OpType op ) : CoefXpr(mp) {
  
  Global::ComplexPart part = a->IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct temp = CoefFunction::Generate( mp_, part, b );
  Init( a, temp, op );
}

CoefXprVecScalOp::CoefXprVecScalOp( MathParser * mp,
                                    PtrCoefFct a, 
                                    const CoefXpr& b, 
                                    CoefXpr::OpType op ) : CoefXpr(mp) {

  Global::ComplexPart part;
  if( a->IsComplex() || b.IsComplex() ) {
    part = Global::COMPLEX;
  } else {
    part = Global::REAL;
  }
  PtrCoefFct temp = CoefFunction::Generate( mp_, part, b );
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
  dependType_ = std::max(a->GetDependency(), 
                         b->GetDependency());
  
  // check for distinct coordinate systems
  CoordSystem* aCoord = a->GetCoordinateSystem();
  CoordSystem* bCoord = b->GetCoordinateSystem();
  if( aCoord != NULL || bCoord != NULL) {
    if (aCoord != NULL && bCoord == NULL) {
      coordSys_ = aCoord;
    } else if (aCoord == NULL && bCoord != NULL) {
      coordSys_ = bCoord;
    } else if (aCoord == bCoord) {
      coordSys_ = aCoord;
    } else {
      EXCEPTION("Case of two distinct coordinate systemw within "\
                "one expression is not tested yet");
    }
  }
  isComplex_ = a->IsComplex() || b->IsComplex();
  a_ = a;
  b_ = b;
  aName_ = CoefXpr::GetUniqueVarName();
  bName_ = CoefXpr::GetUniqueVarName();
  op_ = op;
}

CoefXprTensScalOp::CoefXprTensScalOp( MathParser * mp,
                                      PtrCoefFct a, 
                                      PtrCoefFct b,
                                      CoefXpr::OpType op ) : CoefXpr(mp) {
  Init( a, b, op );
}

CoefXprTensScalOp::CoefXprTensScalOp( MathParser * mp,
                                      PtrCoefFct a, 
                                      const std::string& b,
                                      CoefXpr::OpType op ) : CoefXpr(mp) {
  
  Global::ComplexPart part = a->IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct temp = CoefFunction::Generate( mp_, part, b );
  Init( a, temp, op );
}

CoefXprTensScalOp::CoefXprTensScalOp( MathParser * mp,
                                      PtrCoefFct a, 
                                      const CoefXpr& b, 
                                      CoefXpr::OpType op ) : CoefXpr(mp) {
  
  Global::ComplexPart part;
  if( a->IsComplex() || b.IsComplex() ) {
    part =  Global::COMPLEX;
  } else {
    part = Global::REAL;
  }
  PtrCoefFct temp = CoefFunction::Generate( mp_, part, b );
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
  /*if( numRowsA != 6 || numColsA != 6 ) {
    EXCEPTION( "Tensor must have dimension 6 x 6 " );
  }*/
  
  dimType_ = CoefFunction::TENSOR;
  isAnalytical_ = a->IsAnalytic();
  dependType_ = a->GetDependency();
  coordSys_ = a->GetCoordinateSystem();
  isComplex_ = a->IsComplex();
  a_ = a;
  aName_ = CoefXpr::GetUniqueVarName();
  tensorType_ = NO_TENSOR;
  transposed_ = false;
}

CoefXprMechSubTensor::CoefXprMechSubTensor( MathParser * mp, PtrCoefFct a )
: CoefXpr(mp) {
  Init(a);
}
CoefXprMechSubTensor::CoefXprMechSubTensor( MathParser * mp, const CoefXpr& a) 
: CoefXpr(mp){
  transposed_=false;
  tensorType_ = NO_TENSOR;
  Global::ComplexPart part = a.IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct temp = CoefFunction::Generate( mp_, part, a );  
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
  // TODO: check if material is tricilinic with symmetry plane 1-2, otherwise reduction is not allowed!

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
  } else if( tensorType_ == PLANE_STRESS && numRowsA == 3 && numColsA == 3)  {
    //TODO: 2D Tensor can't be given in mat.xml, dirty fix, only implemented for PLANE_STRESS
      numCols = 3;
      numRows = 3;
      real.Resize( numRows * numCols );
      imag.Resize( numRows * numCols );
      imag.Init("0.0");

      // 1st part: compute same as plane strain representation
      UInt rowPtr[] = {1,2,3};
      for(UInt i = 0; i < numRows; ++i ) {
        for(UInt j = 0; j < numCols; ++j ) {
          real[i*numCols+j] = aR[(rowPtr[i]-1)*3 + (rowPtr[j]-1)];
          if(isComplex_ )
            imag[i*numCols+j] = aI[(rowPtr[i]-1)*3 + (rowPtr[j]-1)];
        }
      }
      if( !isComplex_ ) {
        imag.Init("0.0");
      }
  } else if( tensorType_ == PLANE_STRESS && numRowsA != 3 && numColsA != 3)  {
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
    CoefXpr::Transpose( numRows, numCols, real, tmp );
    real = tmp;
    CoefXpr::Transpose( numRows, numCols, imag, tmp );
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
        (numRowsA == 3 && numColsA == 3 ) ||
        (numRowsA == 6 && numColsA == 6 ) ) ) {
    EXCEPTION( "Tensor must have either dimension 3 x 3, 3 x 6 or 6 x 6" );
  }
  
  dimType_ = CoefFunction::TENSOR;
  isAnalytical_ = a->IsAnalytic();
  isComplex_ = a->IsComplex();
  dependType_ = a->GetDependency();
  a_ = a;
  aName_ = CoefXpr::GetUniqueVarName();
  tensorType_ = NO_TENSOR;
  transposed_ = false;
}

CoefXprSubTensor::CoefXprSubTensor(MathParser * mp, PtrCoefFct a ) 
: CoefXpr(mp) {
  Init(a);
}
CoefXprSubTensor::CoefXprSubTensor( MathParser * mp, const CoefXpr& a) 
: CoefXpr(mp) {
  SetSubTensorType(NO_TENSOR,false);
  Global::ComplexPart part = a.IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct temp = CoefFunction::Generate( mp_, part, a );  
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
    
  } else if( numRowsA == 6 && numColsA == 6 ) {
    // --------------------
    //  6 x 6 TENSOR 
    // --------------------
    if( tensorType_ == FULL ) {
      numRows = 6;
      numCols = 6;
      real = aR;
      imag = aI;
      if( !isComplex_ )
        imag.Init("0.0");
    } else {
      // in all other cases, we just consider the 4x4 submatrix
      numRows = 3;
      numCols = 3;
      real.Resize(9);
      imag.Resize(9);
      UInt indices[] = {0,1,3};
      imag.Init("0.0");
      for( UInt i = 0; i < numRows; ++i ) {
        for( UInt j = 0; j < numCols; ++j ) {
          real[i*numCols+j] = aR[indices[i]*6+indices[j]];
        }
      }
      if( isComplex_ ) {
        for( UInt i = 0; i < numRows; ++i ) {
          for( UInt j = 0; j < numCols; ++j ) {
            imag[i*numCols+j] = aI[indices[i]*6+indices[j]];
          }
        }
      }
    }
  } else if( numRowsA == 3 && numColsA == 6 ) {
    // --------------------
    //  6 x 3 TENSOR 
    // --------------------
    if( tensorType_ == FULL ) {
      numRows = 3;
      numCols = 6;
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
    CoefXpr::Transpose( numRows, numCols, real, tmp );
    real = tmp;
    CoefXpr::Transpose( numRows, numCols, imag, tmp );
    imag = tmp;
    std::swap(numRows, numCols);
  }
 
}

void CoefXprSubTensor::GetArgs( std::map<std::string, PtrCoefFct > & vars ) const {
  vars[aName_] = a_;
}

// ==========================================================================
//  VECTOR REPRESENTATIONS (MECHANIC)
// ==========================================================================
void CoefXprMechSubVector::Init( PtrCoefFct a ) {
  // check dimensionality
  if( a->GetDimType() != CoefFunction::VECTOR ) {
    EXCEPTION( "Argument must be of type VECTOR" );
  }

  // ensure that dimension is a full 6x vectorr
  UInt size = a->GetVecSize();
  if( size != 6  ) {
    EXCEPTION( "Vector must have length 6" );
  }

  dimType_ = CoefFunction::VECTOR;
  isAnalytical_ = a->IsAnalytic();
  dependType_ = a->GetDependency();
  coordSys_ = a->GetCoordinateSystem();
  isComplex_ = a->IsComplex();
  a_ = a;
  aName_ = CoefXpr::GetUniqueVarName();
  tensorType_ = NO_TENSOR;
}

CoefXprMechSubVector::CoefXprMechSubVector( MathParser * mp, PtrCoefFct a )
: CoefXpr(mp) {
  Init(a);
}
CoefXprMechSubVector::CoefXprMechSubVector( MathParser * mp, const CoefXpr& a)
: CoefXpr(mp){
  SetSubTensorType(NO_TENSOR);
  Global::ComplexPart part = a.IsComplex() ? Global::COMPLEX : Global::REAL;
  PtrCoefFct temp = CoefFunction::Generate( mp_, part, a );
}

void CoefXprMechSubVector::SetSubTensorType(SubTensorType subType ) {
  tensorType_ = subType;
}


void CoefXprMechSubVector::GetVectorXpr( StdVector<std::string>& real,
                                     StdVector<std::string>& imag ) const {
  // ensure, that a subtype is set
  if( tensorType_ == NO_TENSOR ) {
    EXCEPTION( "No tensor sub type was set");
  }

  StdVector<std::string> aR, aI;
  UInt size;

  if( isAnalytical_) {
    CoefFunctionAnalytic & coefA =
        dynamic_cast<CoefFunctionAnalytic&>(*a_);
    coefA.GetStrVector( aR, aI );
    size = coefA.GetVecSize();
  } else {
    CoefFunction::GenVecCompNames(aR, aI, aName_, a_);
    size = a_->GetVecSize();
  }

  // switch depending on subTensorType
  if( tensorType_ == AXI ) {
    // resize to 4x Vector
    size = 4;
    real.Resize( size );
    imag.Resize( size );
    imag.Init("0.0");
    UInt rowPtr[] = {1,2,6,3};
    for(UInt i = 0; i < size; ++i ) {
        real[i] = aR[(rowPtr[i]-1)];
        if(isComplex_ )
          imag[i] = aI[(rowPtr[i]-1)];
    }
    // TODO: check if a_5 & a_4 are zero: this is required for the reduction to be allowed!
  } else if( (tensorType_ == PLANE_STRAIN) | (tensorType_ == PLANE_STRESS) ) {
    // resize to 3x Vector
    size = 3;
    real.Resize( size );
    imag.Resize( size );
    imag.Init("0.0");
    UInt rowPtr[] = {1,2,6};
    for(UInt i = 0; i < size; ++i ) {
        real[i] = aR[(rowPtr[i]-1)];
        if(isComplex_ )
          imag[i] = aI[(rowPtr[i]-1)];
    }
    // TODO: check if a_5 & a_4 are zero: this is required for the reduction to be allowed!
  } else if( tensorType_ == FULL ) {
    size = 6;
    real = aR;
    imag = aI;
    if( !isComplex_ )
      imag.Init("0.0");

  } else {
    EXCEPTION( "Desired sub-tensor type not known" );
  }

  // In the end, put everything in brackets
  for( UInt i = 0; i < size; ++i ) {
    real[i] = Bracket(real[i]);
    imag[i] = Bracket(imag[i]);
  }

}

void CoefXprMechSubVector::GetArgs( std::map<std::string, PtrCoefFct > & vars ) const {
  vars[aName_] = a_;
}

} // end of namespace
