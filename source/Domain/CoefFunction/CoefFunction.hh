// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunction.hh
 *       \brief    Base class for describing coefficients
 *
 *       \date     30/10/2011
 *       \author   Andreas Hueppe
 */
//================================================================================================

#ifndef COEFFUNCTION_HH
#define COEFFUNCTION_HH

#include "General/Environment.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/ElemMapping/ElemShapeMap.hh"
#include "Domain/ElemMapping/EntityLists.hh"
namespace CoupledField{


// forward class declarations
class CoordSystem;

//! This is the base class for describing coefficients

//!   It is used in the Integrator classes to obtain the material tensor
//!   (D-Matrix) as well as in the B-Operator classes to realize e.g. a
//!   change in the derivatives. 
//!
//!   Example1: consider the PML
//!   Coefficient function provides a material matrix multiplied by the
//!   complex jacobian. Moreover the same coefficient function is passed to
//!   the operator to provide the damping parameters for the "strechted"
//!   derivatives.
//!
//!   Example2: Consider geometric non-linear mechanics
//!   The CoefFunction can provide the solution at the current integration
//!   point to enable the b-operator to compute its matrix of derivatives
//!
//!   This class is just the interface class right now, the following
//!   structure should be realized:
//!   \li \c coefFunctionAnalytic
//!     This class descibes the quantity by an analytic expression this
//!     includes (non-linear) materials, analytic flow fields, etc.
//!   \li \c coefFunctionGrid
//!     This class provides the Information based on information read in
//!     from a grid. E.g. in the case of (interpolated) aeroacoustic source terms 
//!   \li \c coefFucnctionSol
//!     This function provides the information based on the current
//!     solution of the problem
class CoefFunction{
public:

  // ========================
  //  ENUMS / TYPEDEFS
  // ========================
  //@{ \name Public type definitions

  //! Dimension of coefficient functions
  typedef enum{ 
    NO_DIM,   /*!< Uninitialized */
    SCALAR,   /*!< Scalar entry (cardinality 0)*/
    VECTOR,   /*!< Vector entry (cardinality 1)*/
    TENSOR    /*!< Tensor entry (cardinality 2)*/ 
  } CoefDimType;
  static Enum<CoefDimType> CoefDimType_;
  
  //! Dependency of coefficient function
  typedef enum{ 
    CONST,         /*!< No dependency on space or time */
    TIMEFREQ,      /*!< Only depending on time / frequency, not space */
    GENERAL        /*!< General dependency (spatial and / or time / freq) */
  } CoefDependType;
  static Enum<CoefDependType> CoefDependType_;
  
  // ========================
  //  FACTORY METHODS
  // ========================
  //@{ \name Factory Methods

  //! Generate scalar-valued coefficient function
  
  //! This method generates a scalar-valued coefficient function, either
  //! real- or complex valued.
  //! The method internally investigates the expression and generates
  //! the appropriate coefficient function: If the expression evaluates
  //! to a constant, a CoefFunctionConst-object is created, otherwise
  //! the more general (but more expensive) CoefFunctionExpression.
  //! \param type If COMPLEX, a complex valued CoefFunction is generated;
  //!               If REAL, a real-valued CoefFunction is generated
  //! \param realVal Real-part of the CoefFunction
  //! \param imagVal Imag-part of the CoefFunction (optional)
  static shared_ptr<CoefFunction> 
  Generate( Global::ComplexPart type, 
            const std::string& realVal, 
            const std::string& imagVal = "" );

  //! Generate vector-valued coefficient function
    
  //! This method generates a vector-valued coefficient function, either
  //! real- or complex valued.
  //! The method internally investigates the expression and generates
  //! the appropriate coefficient function: If the expression evaluates
  //! to a constant, a CoefFunctionConst-object is created, otherwise
  //! the more general (but more expensive) CoefFunctionExpression.
  //! \param type If COMPLEX, a complex valued CoefFunction is generated;
  //!               If REAL, a real-valued CoefFunction is generated
  //! \param realVal Real-part vector of the CoefFunction
  //! \param imagVal Imag-part vector of the CoefFunction (optional)
  static shared_ptr<CoefFunction> 
  Generate( Global::ComplexPart type, 
            const StdVector<std::string>& realVal, 
            const StdVector<std::string>& imagVal = StdVector<std::string>() );


  //! Generate tensor-valued coefficient function

  //! This method generates a tensor-valued coefficient function, either
  //! real- or complex valued.
  //! The method internally investigates the expression and generates
  //! the appropriate coefficient function: If the expression evaluates
  //! to a constant, a CoefFunctionConst-object is created, otherwise
  //! the more general (but more expensive) CoefFunctionExpression.
  //! \param type If COMPLEX, a complex valued CoefFunction is generated;
  //!             If REAL, a real-valued CoefFunction is generated
  //! \param numRows Number of rows of the tensor function
  //! \param numCols Number of columns of the tensor function
  //! \param realVal Real-part tensor of the CoefFunction
  //! \param imagVal Imag-part tensor of the CoefFunction (optional)
  static shared_ptr<CoefFunction> 
  Generate( Global::ComplexPart type,
            UInt numRows, UInt numCols,
            const StdVector<std::string>& realVal,
            const StdVector<std::string>& imagVal = StdVector<std::string>() );
  //@}

  // ========================
  //  CONSTRUCTOR / DESTRUCTOR
  // ========================
  //@{ \name Constructor / Destructor

  //! Constructor
  CoefFunction(){
    dimType_ = NO_DIM;
    dependType_ = CONST;
    
    // by default, the coefficients do not
    // depend on any coordiante system
    this->coordSys_ = NULL;
  }

  //! Destructor
  virtual ~CoefFunction(){
    ;
  }
  //@}
  
  // ========================
  //  ACCESS METHODS
  // ========================
  //@{ \name Access Methods

  //! Return real-valued tensor at integration point
  virtual void GetTensor(Matrix<Double>& tensor, 
                         const LocPointMapped& lpm ) {
    EXCEPTION("CoefFunction::GetTensor<Double> called: This may not happen");
  }

  //! Return real-valued vector at integration point
  virtual void GetVector(Vector<Double>& vec, 
                         const LocPointMapped& lpm ) {
    EXCEPTION("CoefFunction::GetVector<Double> called: This may not happen");
  }

  //! Return real-valued scalar at integration point
  virtual void GetScalar(Double& scal, 
                         const LocPointMapped& lpm ) {
    EXCEPTION("CoefFunction::GetScalar<Double> called: This may not happen");
  }

  //! Return complex-valued tensor at integration point
  virtual void GetTensor(Matrix<Complex>& tensor, 
                         const LocPointMapped& lpm ) {
    EXCEPTION("CoefFunction::GetTensor<Complex> called: This may not happen");
  }

  //! Return complex-valued vector at integration point
  virtual void GetVector(Vector<Complex>& vec, 
                         const LocPointMapped& lpm ) {
    EXCEPTION("CoefFunction::GetVector<Complex> called: This may not happen");
  }

  //! Return complex-valued scalar at integration point
  virtual void GetScalar(Complex& scalar, 
                         const LocPointMapped& lpm ) {
    EXCEPTION("CoefFunction::GetScalar<Complex> called: This may not happen");
  }
  //@}


  //! Return associated coordinate system
  void SetCoordinateSystem(CoordSystem* cSys){
    coordSys_ = cSys;
  }
  
  //! Return dependency of CoefFunction
  CoefDependType GetDependency() {
    return dependType_;
  }

  //! Return type of entry (scalar, vector, tensor)
  virtual CoefDimType GetDimType(){
    return dimType_;
  }
  
  //! Returm, if coeffunction is complex
  virtual bool IsComplex() {
    EXCEPTION("Method not properly overwritten");
    return false;
  }

  virtual void AddEntities(shared_ptr<EntityList> ent){
    if(!entities_.Contains(ent)){
      entities_.Push_back(ent);
    }else{
      WARN("entity list " << ent->GetName() << " already contained in CoefFunction")
    }
  }

  virtual void UpdateCoefs(Double tfParam){
    WARN(__FILE__ << __LINE__ << " Update Coefs has no functionality yet");
  }

  //! Dump coefficient function to string 
  virtual std::string ToString() const {
    EXCEPTION("CoefFuncion: ToString() not properly overwritten");
  }


protected:

  // ========================
  //  HELPER METHODS
  // ========================
  //@{ \name Helper methods
  
  //! Returns true, if expression depends on time / freq
  static bool ExprDependsOnTimeFreq(const std::string& expr);
  
  //! Returns true, if expression depends on space
  static bool ExprDependsOnSpace(const std::string& expr);
  //@}
  
  //TODO: CHANGE THIS TO SHARED POINTER
  // i.e. change the daomin to hold a shared_ptr to the coordinate systems!
  CoordSystem* coordSys_;

  //! Dimension of coefficient function (scalar, vector, tensor)
  CoefDimType dimType_;
  
  //! Dependency type of the coefficient function
  CoefDependType dependType_;

  //! Support of the CoefFunction. Only needed for grid/solution results
  StdVector<shared_ptr<EntityList> > entities_;
};

}
#endif
