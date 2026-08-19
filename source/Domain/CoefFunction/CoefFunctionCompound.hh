// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef COEFFUNCTION_COMPOUND_HH
#define COEFFUNCTION_COMPOUND_HH

#include "CoefFunction.hh"
#include "Forms/Operators/IdentityOperator.hh"
#include "Forms/Operators/DivOperator.hh"
#include "Utils/ThreadLocalStorage.hh"
#include "Utils/mathParser/mathParser.hh"
//#include "FeBasis/FeFunctions.hh"
#include "FeBasis/H1/FeSpaceH1Nodal.hh"

namespace CoupledField {

//! Coefficient class, which is formed by a compound of several others

//! This coefficient class can calculate arbitrary expressions with other 
//! CoefFunction(s) as arguments. The expression itself is written as 
//! MathParser expression. The arguments of the expression get registered
//! by a map with the variable name as key and the pointer to the CoefFunction
//! as value. 
//! Internally, a MathParser instance is generated to store the scalar / vector
//! / tensor expression and all argument CoefFunctions are registered as 
//! external variables to the MathParser. For every argument CoefFunction, 
//! internally a //! scalar / vector / tensor is created and the entrie(s) are 
//! registered as external variables to the MathParser.
//!  
//! Upon access via CoefFunctionCompound::GetScalar,
//! CoefFunctionCompound::GetVector or CoefFunctionCompound::GetTensor,
//! every argument (= CoefFunction) gets first evaluated and stored into its
//! scalar / vector / tensor variable (thus the variables within the MathParser
//! get updated "transparently") and then the own expression is evaluated.
//!
//! To assign variable names to the single scalar / vector / tensor entries,
//! the methods CoefFunction::GenScalCompNames, CoefFunction::GenVecCompNames
//! and CoefFunction::GenTensorCompNames are used.
//! 
//! In the following example, a vector Coefficient function a and scalar 
//! CoefFunction b get multiplied:
//!
//! \code
//! StdVector<std::string> aVals;
//! aVals = "2.0", "3.0";
//! PtrCoefFct a = 
//! CoefFunction::Generate(Global::REAL, aVals );
//! PtrCoefFct b = 
//! CoefFunction::Generate(Global::REAL, "4.0" );
//! 
//! CoefFunctionCompound c;
//! std::map<std::string, PtrCoefFct> vars;
//! vars["a"] = a;
//! vars["b"] = b;
//! 
//! // Generate expression for vector*scalar multiplication. 
//! // The single vector components are split by ','. Indices of vectors / scalars are 0-based. 
//! // The real part is denoted by the suffix "_R".
//! // This is the same convention as in the methods CoefFunction::GenScalCompNames()
//! std::string vecXpr = "(a_0_R * b_R), (a_1_R * b_R);"
//! c.SetVector( vecXpr, vars);
//! \endcode
//! 
//! \note Objects of this class normally get created using CoefXpr classes to
//!       form the mathematical expression and then get generated using
//!       the method CoefFunction::Generate(Global::ComplexPart, const CoefXpr&).
template<class DATA_TYPE>
class CoefFunctionCompound : public CoefFunction {
};

// ===========================================================================
//  REAL VALUED COEFFICIENT FUNCTION
// ===========================================================================
//! Coefficient class, which is formed as a compound of several others - real-valued
template<>
class CoefFunctionCompound<Double> : public CoefFunction {
public:

  //! Constructor
  CoefFunctionCompound(MathParser * mp);

  //! Destructor
  virtual ~CoefFunctionCompound();

  virtual string GetName() const { return "CoefFunctionCompound_double"; }

  //! \copydoc CoefFunction::GetTensor
  void GetTensor( Matrix<Double>& coefMat, const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetVector
  void GetVector( Vector<Double>& coefVec, const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetScalar
  void GetScalar( Double& coefScalar, const LocPointMapped& lpm );

  //! Set Scalar valued expression
  void SetScalar( const std::string& expr,
                  std::map<std::string, PtrCoefFct >& vars );

  //! Set Vector valued expression
  void SetVector( StdVector<std::string>& expr, 
                  std::map<std::string, PtrCoefFct >& vars );

  //! Set Tensor valued expression
  void SetTensor( StdVector<std::string>& expr,
                  UInt numRows, UInt numCols,
                  std::map<std::string, PtrCoefFct >& vars );

  //! Return reference to coefs_ map. coefs_ can be changed!
  std::map<std::string, PtrCoefFct>& GetCoefFcts() { return coefs_; }

  //! \copydoc CoefFunction::GetVecSize
  virtual UInt GetVecSize() const;

  //! \copydoc CoefFunction::GetTensorSize
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const;
     
  //! \copydoc CoefFunction::ToString
  std::string ToString() const;

  //! \copydoc CoefFunction::SetDerivativeOperation
  virtual void SetDerivativeOperation(CoefDerivativeType type, UInt gDim, UInt dDim){
    this->derivType_ = type;

    //make some checks here!
    switch(dimType_){
    case SCALAR:
      //only NONE is valid right now
      //if extended to gradient, this would be fine too
      if(type==VECTOR_DIVERGENCE){
        EXCEPTION("CoefFunctionExpression: VECTOR_DIVERGENCE is not a valid operator for scalar coefFunction");
      }
      break;
    case VECTOR:
      //this is fine in all cases right now
      if(type==VECTOR_DIVERGENCE){
        //change dim type to scalar
        this->dimType_ = SCALAR;
//          UInt gDim = this->domain_->GetGrid()->GetDim();
//          UInt dDim = this->resultInfo_->dofNames.GetSize();
        this->CreateDivOperator(gDim,dDim);
      }
      break;
    case TENSOR:
      if(type==VECTOR_DIVERGENCE){
        EXCEPTION("CoefFunctionExpression: VECTOR_DIVERGENCE is not a valid operator for tensor coefFunction");
      }
      break;
    default:
      break;
    }
    return;
  }

protected:
  
  //! Register coefficient function
  void RegisterCoefFct( const std::string& name,
                        PtrCoefFct& coef );
  
  //! Update expressions
  void UpdateXpr( const LocPointMapped& lpm );

  //!Creates a divergence operator in case we want the divergence of an
  //! external vector field
  void CreateDivOperator(UInt spaceDim, UInt dofDim);

  //! BOperator to map solutions to arbitrary points
  //! Right now, hardcoded identity operator
  shared_ptr<BaseBOperator > myOperator_;

  //! Number of rows of tensor
  UInt numRows_;

  //! Number of columns of tensor
  UInt numCols_;
  
  //! Expression containing the result
  std::string expr_;
  
  //! Mathparser handle
  unsigned int handle_;
  
  //! Mathparser object
  MathParser* parser_;
  
  //! Default coordinate system 
  CoordSystem* coordSysDefault_;
  
  //! Map variable names to coefficient functions
  std::map<std::string, PtrCoefFct > coefs_;

  //! Map variable names to dimensionality types
  std::map<std::string, CoefDimType > coefDimTypes_;
  
  //! Map variable names to scalar valued variables
  TLMap<std::string, Double> scalVars_;

  //! Map variable names to vector valued variables
  TLMap<std::string, Vector<Double> > vecVars_;

  //! Map variable names to tensor valued variables
  TLMap<std::string, Matrix<Double> > tensorVars_;
};

// ===========================================================================
//  Complex VALUED COEFFICIENT FUNCTION
// ===========================================================================
//! Coefficient class, which is formed as a compound of several others - complex-valued
template<>
class CoefFunctionCompound<Complex> : public CoefFunction {
public:

  //! Constructor
  CoefFunctionCompound(MathParser * mp);

  //! Destructor
  virtual ~CoefFunctionCompound();

  virtual string GetName() const { return "CoefFunctionCompound_complex"; }

  //! \copydoc CoefFunction::GetTensor
  void GetTensor( Matrix<Complex>& coefMat, const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetVector
  void GetVector( Vector<Complex>& coefVec, const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetScalar
  void GetScalar( Complex& coefScalar, const LocPointMapped& lpm );

  //! Set Scalar valued expression
  void SetScalar( std::string& realExpr,
                  std::string& imagExpr,
                  std::map<std::string, PtrCoefFct >& vars );

  //! Set Vector valued expression
  void SetVector( StdVector<std::string>& realExpr,
                  StdVector<std::string>& imagExpr,
                  std::map<std::string, PtrCoefFct >& vars );

  //! Set Tensor valued expression
  void SetTensor( StdVector<std::string>& realExpr,
                  StdVector<std::string>& imagExpr,
                  UInt numRows, UInt numCols,
                  std::map<std::string, PtrCoefFct >& vars );

  //! Return reference to coefs_ map. coefs_ can be changed!
  std::map<std::string, PtrCoefFct>& GetCoefFcts() { return coefs_; }

  //! \copydoc CoefFunction::GetVecSize
  virtual UInt GetVecSize() const;

  //! \copydoc CoefFunction::GetTensorSize
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const;
     
  //! \copydoc CoefFunction::ToString
  std::string ToString() const;
  
        //! \copydoc CoefFunction::SetDerivativeOperation
    virtual void SetDerivativeOperation(CoefDerivativeType type, UInt gDim, UInt dDim){
      this->derivType_ = type;

      //make some checks here!
      switch(dimType_){
      case SCALAR:
        //only NONE is valid right now
        //if extended to gradient, this would be fine too
        if(type==VECTOR_DIVERGENCE){
          EXCEPTION("CoefFunctionExpression: VECTOR_DIVERGENCE is not a valid operator for scalar coefFunction");
        }
        break;
      case VECTOR:
        //this is fine in all cases right now
        if(type==VECTOR_DIVERGENCE){
          //change dim type to scalar
          this->dimType_ = SCALAR;
//          UInt gDim = this->domain_->GetGrid()->GetDim();
//          UInt dDim = this->resultInfo_->dofNames.GetSize();
          this->CreateDivOperator(gDim,dDim);
        }
        break;
      case TENSOR:
        if(type==VECTOR_DIVERGENCE){
          EXCEPTION("CoefFunctionExpression: VECTOR_DIVERGENCE is not a valid operator for tensor coefFunction");
        }
        break;
      default:
        break;
      }
      return;
    }

protected:
  
  //! Register coefficient function
  void RegisterCoefFct( const std::string& name,
                        PtrCoefFct& coef );
  
  //! Update expressions
  void UpdateXpr( const LocPointMapped& lpm );


  //!Creates a divergence operator in case we want the divergence of an
  //! external vector field
  void CreateDivOperator(UInt spaceDim, UInt dofDim);

  //! BOperator to map solutions to arbitrary points
  //! Right now, hardcoded identity operator
  shared_ptr<BaseBOperator > myOperator_;

  //! Number of rows of tensor
  UInt numRows_;

  //! Number of columns of tensor
  UInt numCols_;
  
  //! Expression containing the real-valued results
  std::string exprReal_;
  
  //! Expression containing the imag-valued results
  std::string exprImag_;
  
  //! Mathparser handle (real part)
  unsigned int handleReal_;
  
  //! Mathparser handle (imag part)
  unsigned int handleImag_;
  
  //! Mathparser object
  MathParser* parser_;
  
  //! Default coordinate system 
  CoordSystem* coordSysDefault_;
  
  //! Map variable names to coefficient functions
  std::map<std::string, PtrCoefFct > coefs_;

  //! Map variable names to dimensionality types
  std::map<std::string, CoefDimType > coefDimTypes_;
  
  //=================================================
  // Make this thread safe as we relate to the address
  //=================================================

  //! Map variable names to scalar valued variables (real part)
  TLMap<std::string, Double> scalVarsReal_;
  
  //! Map variable names to scalar valued variables (imag part)
  TLMap<std::string, Double> scalVarsImag_;

  //! Map variable names to vector valued variables (real part)
  TLMap<std::string, Vector<Double> > vecVarsReal_;
  
  //! Map variable names to vector valued variables (imag part)
  TLMap<std::string, Vector<Double> > vecVarsImag_;

  //! Map variable names to tensor valued variables (real part)
  TLMap<std::string, Matrix<Double> > tensorVarsReal_;

  //! Map variable names to tensor valued variables (imag part )
  TLMap<std::string, Matrix<Double> > tensorVarsImag_;
};


} // end of namespace

#endif
