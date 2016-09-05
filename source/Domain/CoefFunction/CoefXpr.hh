#ifndef COEF_VECTOR_XPR_HH
#define COEF_VECTOR_XPR_HH

#include "CoefFunction.hh"

namespace CoupledField {

//! Base class for describing expressions based on coefficient functions

//! This is the base class for describing general mathematical expressions,
//! which are composed of CoefFunction objects.
//! They can be converted to normal CoefFunction objects using the method
//! CoefFunction::Generate( Global::ComplexPart, const CoefXpr&).
//! Similar to a normal mathematical formula, which consists of operations,
//! applied to arguments, we establish a similar structure for coefficient 
//! functions.
//! Thus, this class offers e.g. operations ( CoefXpr::OpType), which 
//! expresses the type of operation to be performed. It mostly gets passed
//! as arguments to the derived classes, which describe, e.g. operations with
//! one argument ( CoefXprUnaryOp), two arguments (CoefXprBinaryOp)
//! or binary operations on vectors-scalar ( CoefXprVecScalOp).
//!
//! The elementary expressions can be composed by passing them as arguments
//! to other coefficient expression.
//! 
//! When querying the final expression (CoefXpr::GetScalarXpr, 
//! CoefXpr::GetVectorXpr, CoefXpr::GetTensorXpr), we have to distinguish, 
//! if the resulting expression
//! can be written in a closed string form (e.g. all operands involved are 
//! themselves analytical expressions, CoefXpr::IsAnalytical) or if 
//! it involves non-analytical expressions (e.g. sampled data or a 
//! nonlinear curve):
//! 
//!   - If the expression is analytic, the resulting string / string vector
//!     can be directly used by the CoefFunction::Generate() method,
//!     to build again an analytic expression. This is the most efficient way.
//! 
//!   - If the expression is not analytic (e.g. one argument is a nonlinear
//!     function), the result can only be used in a  CoefFunctionCompound.
//!     In this case, the expression returned by the Get-methods are just
//!     formulas, which reference the single entries of the arguments.
//!     Consider e.g. the case of an addition of two vector-valued
//!     coefficient functions a and b (each having two entries). In this case,
//!     the returned expression from  CoefXprBinaryOp::GetVectorXpr would
//!     be the component-wise formula:
//!       (a_0 + b_0), (a_1 + b_1)
//!     This formula can be evaluated in an MathParser expression, if the variables
//!     (a_0, a_1, b_0, b_1) get registered. This can be accomplished by
//!     getting the arguments of an expression using  CoefXpr::GetArgs
//!     and using the method  CoefFunction
//!     The "arguments" a and b can be retrieved using e.g. the method
//!     CoefFunction::GenVecCompNames, which returns a list of component names
//!     for a given vector-valued coefficient function.
class CoefXpr {
public:
  
    //! Constructor
    CoefXpr(MathParser * mp) {
      isAnalytical_ = false;
      isComplex_ = false;
      dependType_ = CoefFunction::CONSTANT;
      coordSys_ = NULL;
      dimType_ = CoefFunction::NO_DIM;
      mp_ = mp;
    }
    
    //! Destructor
    virtual ~CoefXpr() {
    }
    
    //! Get dimension of expression
    CoefFunction::CoefDimType GetDimType() const {
      return dimType_;
    }
    
    //! Query, if expression is analytical
    
    //! This method returns, if the composed expression is still analytical, 
    //! i.e. it can be written in a closed formula (e.g. all variables are 
    //! linear) or if not (e.g. a sampled function or nonlinear approximation).
    //! In the first case, this expression can be efficiently converted in
    //! one of the classes  CoefFunctionConst, CoefFunctionTimeFreq
    //! or CoefFunctionExpression. 
    //! Otherwise, only a CoefFunctionCompound can be generated.
    bool IsAnalytical() const {
      return isAnalytical_;
    }
    
    //! Return dependency of CoefFunction
    CoefFunction::CoefDependType GetDependency() const {
      return dependType_;
    }
    
    //! Query, if expression is complex-valued
    bool IsComplex() const {
      return isComplex_;
    }
    
    //! Get associated coordinate system
    CoordSystem* GetCoordinateSystem() const {
      return coordSys_;
    }
    
    //! Get MathParser instance
    MathParser * GetMathParser() {
      return mp_;
    }
  
  // --------------------------------------------------------------------------
  // GENERAL STATIC METHODS
  // --------------------------------------------------------------------------
  //! Typedef for defined mathematical operations
  typedef enum {
    NOOP,                     /*!< No operation */
    OP_ADD,                   /*!< Binary + operation */
    OP_SUB,                   /*!< Binary - operation */
    OP_MULT,                  /*!< Binary * operation (scal-scal, scalar-vector)
                                   or inner product (vector-vector) */
    OP_MULT_CONJ,             /*!< Binary * operation (scal-scal, scalar-vector), conjugated */
    OP_MULT_VOIGT_TENSOR_VEC, /*!< Binary * operation (tensor-vector, Voigt case) */
    OP_MULT_VOIGT_TENSOR_VEC_CONJ, /*!< Binary * operation (tensor-vector, Voigt case), conjugated */
    OP_MULT_TENSOR,           /*!< Binary * operation (tensor-tensor) e.g. dyad in case of vector-vector */
    OP_DIV,                   /*!< Binary / operation */
    OP_CROSS,                 /*!< Binary x operation (cross product 3D and 2D) */
    OP_CROSS_AXI,             /*!< Binary x operation (axisymmetric cross product) */
    OP_POW,                   /*!< Binary x^y operation */
    OP_NORM,                  /*!< Unary L2-Norm operation */
    OP_SQRT,                   /*!< Unary square root operation */
    OP_TRACE,                  /*!< Unary trace of tensor operation  */
    OP_INV,                   /*!< Unary inversion operation */
    OP_DET                    /*!< Unary operation, returning the determinant of a square tensor*/
  } OpType;
  
  //! Get number of operands for OpType
  static UInt GetNumOperands( OpType op );
  
  //! Get string representation of operator
  static std::string OpToString( OpType op );
  
  //! Return dimensionality of result of unary expression
  static CoefFunction::CoefDimType GetDimType( PtrCoefFct a,
                                               OpType op ); 
  
  //! Return dimensionality of result of unary expression
  static CoefFunction::CoefDimType GetDimType( PtrCoefFct a,
                                               PtrCoefFct b,
                                               OpType op ); 

  //! Return if result of unary expression is complex
  static bool IsComplex( PtrCoefFct a,
                         OpType op ); 

  //! Return if result of binary expression is complex
  static bool IsComplex( PtrCoefFct a,
                         PtrCoefFct b,
                         OpType op ); 
  
  //! Helper method to create the transposed of a tensor / matrix
  
  //! This method returns the transposed of the original vector. The parameters
  //! nCols and Rows refer to the constant original matrix / tensor.
  static void Transpose( UInt nRows, UInt nCols, 
                        const StdVector<std::string>& in,
                        StdVector<std::string>& trans );
                        

  //! Apply real-valued unary function on scalar
  static void ApplyUnaryFunc( std::string& retReal, const std::string& argReal,
                              OpType op );

  //! Apply complex-valued unary function on scalar
  static void ApplyUnaryFunc( std::string& retReal, std::string& retImag, 
                              const std::string& argReal,
                              const std::string& argImag,
                              OpType op );

  //! Apply real-valued binary function on scalar
  static void ApplyBinaryFunc( std::string& retReal, 
                               const std::string& arg1Real,
                               const std::string& arg2Real,
                               OpType op );

  //! Apply complex-valued binary function on scalar
  static void ApplyBinaryFunc( std::string& retReal, std::string& retImag, 
                               const std::string& arg1Real,
                               const std::string& arg2Real,
                               const std::string& arg1Imag,
                               const std::string& arg2Imag,
                               OpType op );
  
  //! Apply real-valued ternary function on scalar
  
  //! Calculates the expression  op1(arg1, op2(arg2, arg3) )
  static void ApplyTernaryFunc( std::string& retReal, 
                                const std::string& arg1Real,
                                const std::string& arg2Real,
                                const std::string& arg3Real,
                                OpType op1, OpType op2 );

  //! Apply complex-valued ternary function on scalar
  
  //! Calculates the expression  op1(arg1, op2(arg2, arg3) )
  static void ApplyTernaryFunc( std::string& retReal, std::string& retImag, 
                                const std::string& arg1Real,
                                const std::string& arg2Real,
                                const std::string& arg3Real,
                                const std::string& arg1Imag,
                                const std::string& arg2Imag,
                                const std::string& arg3Imag,
                                OpType op1, OpType op2 );

  

  // --------------------------------------------------------------------------
  // Case 1: Expression evaluates to an analytical expressions
  // --------------------------------------------------------------------------

  //! Get scalar expression
  virtual void GetScalarXpr( std::string& real, std::string& imag ) const {
    EXCEPTION( "Not implemented here");
  }

  //! Get vector expression
  virtual void GetVectorXpr( StdVector<std::string>& real, 
                             StdVector<std::string>& imag ) const {
    EXCEPTION( "Not implemented here");
  }

  //! Get tensor expression
  virtual void GetTensorXpr( UInt& numRows, UInt& numCols,
                             StdVector<std::string>& real, 
                             StdVector<std::string>& imag ) const {
    EXCEPTION( "Not implemented here");
  }

  // --------------------------------------------------------------------------
  // Case 2: General expression, which is not analytical
  // --------------------------------------------------------------------------
  
  //! Return coefficient functions and their internal variable name 
  
  //! This method returns a map, containing as key the internal variable
  //! name and as value the CoefFunction. 
  //! If e.g. the non-analytic result of a scalar expression is "a+b",
  //! this method returns a map with two entries: One entry with 
  //! key "a" with the coefficient function representing the first argument
  //! and one entry with key "b", representing the coefficient function for
  //! the second argument.
  virtual void GetArgs( std::map<std::string, PtrCoefFct >& vars ) const = 0;
  
protected:
  
  //! Obtain a unique variable name
  
  //! This method returns on each call a unique variable name, which can be 
  //! used in expressions. This is necessary, as we can not re-use variable
  //! names in nested expressions.
  static std::string GetUniqueVarName();
  
  //! Dimension of expression
  CoefFunction::CoefDimType dimType_;
  
  //! Flag, if expression evaluates to an analytical expression
  bool isAnalytical_;
  
  //! Flag, if expression is complex-valued
  bool isComplex_;
  
  //! Dependency type of the coefficient function
  CoefFunction::CoefDependType dependType_;
  
  //! Pointer to coordinate system
  CoordSystem* coordSys_;
  
  //! Pointer to MathParser
  MathParser * mp_;
  
private:
  
  //! Internal counter for used generated variable names
  static UInt numVarNames_;
  
};

// --------------------------------------------------------------------------
//  COEFFICIENT EXPRESSIONS WITH ONE ARGUMENT (UNARY)
// --------------------------------------------------------------------------

//! Unary operation applied to a CoefFunction

//! This class wraps all expressions, which need just one CoefFunction as 
//! as argument. This are for example the operations 
//!   - CoefXpr::OP_NORM
//!   - CoefXpr::OP_SQRT
class CoefXprUnaryOp : public CoefXpr {
public:
  
  //! Constructor
  CoefXprUnaryOp( MathParser * mp,
                  PtrCoefFct a,
                  CoefXpr::OpType op );

  //! Constructor for string representation of coefficient function
  CoefXprUnaryOp( MathParser * mp,
                  const std::string& a,
                  CoefXpr::OpType op );

  //! Constructor for string and coefficient function
  CoefXprUnaryOp( MathParser * mp,
                  const CoefXpr& a,
                  CoefXpr::OpType op );

  //! Get scalar expression
  void GetScalarXpr( std::string& real, std::string& imag ) const;

  //! Get vector expression
  void GetVectorXpr( StdVector<std::string>& real, 
                     StdVector<std::string>& imag ) const;

  //! Get tensor expression
  void GetTensorXpr( UInt& numRows, UInt& numCols,
                     StdVector<std::string>& real, 
                     StdVector<std::string>& imag ) const;

  //! \copydoc CoefXpr::GetArgs
  void GetArgs( std::map<std::string, PtrCoefFct > & vars ) const;
  
  //! Get operand type
  CoefXpr::OpType GetOpType() const {
    return op_;
  }
  
protected:
  
  //! Private initialization
   void Init( PtrCoefFct a, 
              CoefXpr::OpType op );
  
  //! Coefficient function representing the first (and only) argument
  PtrCoefFct a_;
  
  //! Variable name of the first argument
  std::string aName_;
  
  //! Operator
  CoefXpr::OpType op_;

};

// --------------------------------------------------------------------------
//  COEFFICIENT EXPRESSIONS WITH TWO ARGUMENTS (BINARY)
// --------------------------------------------------------------------------

//! Binary operation between two coefficient functions

//! This class wraps up all expressions, which are composed of
//! one binary operation and its two arguments. Note, that both
//! arguments must be of same dimension (e.g. scalar-scalar, vector-vector).
//! Otherwise the classes CoefXprVecScalOp or CoefXprTensScalOp
//! must be used.
//!
//! To generate e.g. a binary expression "coef1 / coef2", write:
//! \code
//! CoefXprBinOp xpr(coef1, coef2, CoefXpr::OP_DIV); 
//! \endcode

class CoefXprBinOp : public CoefXpr {
public:
  //! Constructor
  CoefXprBinOp( MathParser * mp,
                PtrCoefFct a, 
                PtrCoefFct b,
                CoefXpr::OpType op );
  
  //! Constructor for coefficient and xpr
  CoefXprBinOp( MathParser * mp,
                PtrCoefFct a, 
                const CoefXpr& b,
                CoefXpr::OpType op );
  
  //! Constructor for xpr and coefficient
  CoefXprBinOp( MathParser * mp,
                const CoefXpr& a,
                PtrCoefFct b, 
                CoefXpr::OpType op );

  //! Constructor for coefficient function, string
  CoefXprBinOp( MathParser * mp,
                PtrCoefFct a, 
                const std::string& b,
                CoefXpr::OpType op );

  //! Constructor for string and coefficient function
  CoefXprBinOp( MathParser * mp,
                const std::string& a, 
                PtrCoefFct b, 
                CoefXpr::OpType op );
  
  //! Constructor for string and coefficient function
  CoefXprBinOp( MathParser * mp,
                const std::string& a, 
                const CoefXpr& b, 
                CoefXpr::OpType op );
  
  //! Get scalar expression
  void GetScalarXpr( std::string& real, std::string& imag ) const;

  //! Get vector expression
  void GetVectorXpr( StdVector<std::string>& real, 
                     StdVector<std::string>& imag ) const;

  //! Get tensor expression
  void GetTensorXpr( UInt& numRows, UInt& numCols,
                     StdVector<std::string>& real, 
                     StdVector<std::string>& imag ) const;

  //! \copydoc CoefXpr::GetArgs
  void GetArgs( std::map<std::string, PtrCoefFct >& vars ) const ;
  
  //! Get operand type
  CoefXpr::OpType GetOpType() const {
    return op_;
  }
  
protected:
  
  //! Private initialization
  void Init( PtrCoefFct a, 
             PtrCoefFct b,
             CoefXpr::OpType op );
  
  //! Coefficient function A
  PtrCoefFct a_;
  
  //! Coefficient function B
  PtrCoefFct b_;
  
  //! Variable name of the first argument
  std::string aName_;

  //! Variable name of the second argument
  std::string bName_;

  //! Operator
  CoefXpr::OpType op_;

};

// --------------------------------------------------------------------------
//  VECTOR-SCALAR COEFFICIENT EXPRESSIONS WITH TWO ARGUMENTS (BINARY)
// --------------------------------------------------------------------------

//! Binary operation between vector and scalar coefficients

//! This class wraps up binary expressions, which take as a first argument
//! a vector-valued CoefFunction and as second argument a scalar CoefFunction.
//!
//! The only operations allowed are CoefXpr::OP_MULT and CoefXpr::OP_DIV.
class CoefXprVecScalOp : public CoefXpr {
public:
  //! Constructor
  CoefXprVecScalOp( MathParser * mp,
                    PtrCoefFct a, 
                    PtrCoefFct b,
                    CoefXpr::OpType op );

  //! Constructor for coefficient function, string
  CoefXprVecScalOp( MathParser * mp,
                    PtrCoefFct a, 
                    const std::string& b,
                    CoefXpr::OpType op );

  //! Constructor for string and coefficient function
  CoefXprVecScalOp( MathParser * mp,
                    PtrCoefFct a, 
                    const CoefXpr& b, 
                    CoefXpr::OpType op );
  
  //! Get scalar expression
  void GetScalarXpr( std::string& real, std::string& imag ) const;

  //! Get vector expression
  void GetVectorXpr( StdVector<std::string>& real, 
                     StdVector<std::string>& imag ) const;

  //! Get tensor expression
  void GetTensorXpr( UInt& numRows, UInt& numCols,
                     StdVector<std::string>& real, 
                     StdVector<std::string>& imag ) const;

  //! \copydoc CoefXpr::GetArgs
  void GetArgs( std::map<std::string, PtrCoefFct >& vars ) const;

  //! Get operand type
  CoefXpr::OpType GetOpType() const {
    return op_;
  }

protected:

  //! Private initialization
  void Init( PtrCoefFct a, 
             PtrCoefFct b,
             CoefXpr::OpType op );
  
  //! Coefficient function A
  PtrCoefFct a_;
  
  //! Coefficient function B
  PtrCoefFct b_;
  
  //! Variable name of the first argument
  std::string aName_;

  //! Variable name of the second argument
  std::string bName_;
  
  //! Operator
  CoefXpr::OpType op_;

};

// --------------------------------------------------------------------------
//  TENSOR-SCALAR COEFFICIENT EXPRESSIONS WITH TWO ARGUMENTS (BINARY)
// --------------------------------------------------------------------------

//! Binary operation between tensor and scalar coefficients

//! This class wraps up binary expressions, which take as a first argument
//! a tensor-valued CoefFunction and as second argument a scalar CoefFunction.
//!
//! The only operations allowed are CoefXpr::OP_MULT and CoefXpr::OP_DIV.
class CoefXprTensScalOp : public CoefXpr {
public:
  //! Constructor
  CoefXprTensScalOp( MathParser * mp,
                     PtrCoefFct a, 
                     PtrCoefFct b,
                     CoefXpr::OpType op );

  //! Constructor for coefficient function, string
  CoefXprTensScalOp( MathParser * mp,
                     PtrCoefFct a, 
                     const std::string& b,
                     CoefXpr::OpType op );

  //! Constructor for string and coefficient function
  CoefXprTensScalOp( MathParser * mp,
                     PtrCoefFct a, 
                     const CoefXpr& b, 
                     CoefXpr::OpType op );

  //! Get scalar expression
  void GetScalarXpr( std::string& real, std::string& imag ) const;

  //! Get vector expression
  void GetVectorXpr( StdVector<std::string>& real, 
                     StdVector<std::string>& imag ) const;

  //! Get tensor expression
  void GetTensorXpr( UInt& numRows, UInt& numCols,
                     StdVector<std::string>& real, 
                     StdVector<std::string>& imag ) const;

  //! \copydoc CoefXpr::GetArgs
  void GetArgs( std::map<std::string, PtrCoefFct >& vars ) const;

  //! Get operand type
  CoefXpr::OpType GetOpType() const {
    return op_;
  }

protected:

  //! Private initialization
  void Init( PtrCoefFct a, 
             PtrCoefFct b,
             CoefXpr::OpType op );
  
  //! Coefficient function A
  PtrCoefFct a_;
  
  //! Coefficient function B
  PtrCoefFct b_;
  
  //! Variable name of the first argument
  std::string aName_;

  //! Variable name of the second argument
  std::string bName_;
  
  //! Operator
  CoefXpr::OpType op_;
};

// --------------------------------------------------------------------------
//  TENSOR REPRESENTATION (MECHANIC MATERIAL)
// --------------------------------------------------------------------------
//! Models the sub-tensor representation of mechanical stiffness

//! This class represents the different tensor representations of the 
//! mechanical stiffness tensor. It takes the original, full tensor and
//! returns the sub-tensor in general form 
class CoefXprMechSubTensor : public CoefXpr {
  
public:
  
  //! Constructor
  CoefXprMechSubTensor( MathParser * mp,
                        PtrCoefFct a );

  //! Constructor for given coefficient function
  CoefXprMechSubTensor( MathParser * mp,
                        const CoefXpr& a) ;
  
  //! Set given subType and if tensor should be transposed
  void SetSubTensorType(SubTensorType subType, bool transposed );
  

  //! Get tensor expression
  void GetTensorXpr( UInt& numRows, UInt& numCols,
                     StdVector<std::string>& real, 
                     StdVector<std::string>& imag ) const;

  //! \copydoc CoefXpr::GetArgs
  void GetArgs( std::map<std::string, PtrCoefFct > & vars ) const;
   
protected:
   
   //! Private initialization
   void Init( PtrCoefFct );
   
   //! Coefficient function representing the original tensor
   PtrCoefFct a_;
   
   //! Variable name of the first argument
   std::string aName_;
   
   //! Subtensor type
   SubTensorType tensorType_;
   
   //! Flag if transposed of tensor should be used
   bool transposed_;
};

// --------------------------------------------------------------------------
//  TENSOR REPRESENTATION (GENERAL MATERIAL)
// --------------------------------------------------------------------------
//! Models the sub-tensor representation of general tensors

//! This class represents the different tensor representations of the 
//! general material tensors in 3 dimensions. It takes the original, full tensor 
//! and returns the sub-tensor in general form  
class CoefXprSubTensor : public CoefXpr {
  
public:
  
  //! Constructor
  CoefXprSubTensor( MathParser * mp,
                    PtrCoefFct a );

  //! Constructor for given coefficient function
  CoefXprSubTensor( MathParser * mp,
                    const CoefXpr& a) ;
  
  //! Set given subType and if tensor should be transposed
  void SetSubTensorType(SubTensorType subType, bool transposed );
  

  //! Get tensor expression
  void GetTensorXpr( UInt& numRows, UInt& numCols,
                     StdVector<std::string>& real, 
                     StdVector<std::string>& imag ) const;

  //! \copydoc CoefXpr::GetArgs
  void GetArgs( std::map<std::string, PtrCoefFct > & vars ) const;
   
protected:
   
   //! Private initialization
   void Init( PtrCoefFct );
   
   //! Coefficient function representing the original tensor
   PtrCoefFct a_;
   
   //! Variable name of the first argument
   std::string aName_;
   
   //! Subtensor type
   SubTensorType tensorType_;
   
   //! Flag if transposed of tensor should be used
   bool transposed_;
};


// --------------------------------------------------------------------------
//  TENSOR REPRESENTATION (Voigt Notation)
// --------------------------------------------------------------------------
//! Models the sub-tensor representation of mechanic tensors in Voigt notation

//!
class CoefXprMechSubVector : public CoefXpr {

public:

  //! Constructor
  CoefXprMechSubVector( MathParser * mp,
                        PtrCoefFct a );

  //! Constructor for given coefficient function
  CoefXprMechSubVector( MathParser * mp,
                        const CoefXpr& a) ;

  //! Set given subType and if tensor should be transposed
  void SetSubTensorType(SubTensorType subType );


  //! Get vector expression
  void GetVectorXpr( StdVector<std::string>& real,
                     StdVector<std::string>& imag ) const;

  //! \copydoc CoefXpr::GetArgs
  void GetArgs( std::map<std::string, PtrCoefFct > & vars ) const;

protected:

   //! Private initialization
   void Init( PtrCoefFct );

   //! Coefficient function representing the original tensor
   PtrCoefFct a_;

   //! Variable name of the first argument
   std::string aName_;

   //! Subtensor type
   SubTensorType tensorType_;

};

}
#endif // header guard
