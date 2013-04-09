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

#include <boost/utility.hpp>

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
class CoefXpr;
class CoefFunction;
class FeSpace;

//! This is the base class for describing coefficients

//!   It is used in the Integrator classes to obtain the material tensor
//!   (D-Matrix) as well as in the B-Operator classes to realize e.g. a
//!   change in the derivatives. 
//!
//!   Example1: consider the PML
//!   Coefficient function provides a material matrix multiplied by the
//!   complex Jacobian. Moreover the same coefficient function is passed to
//!   the operator to provide the damping parameters for the "stretched"
//!   derivatives.
//!
//!   Example2: Consider geometric non-linear mechanics
//!   The CoefFunction can provide the solution at the current integration
//!   point to enable the b-operator to compute its matrix of derivatives
//!
//!   This class is just the interface class right now, the following
//!   structure should be realized:
//!   \li \c coefFunctionAnalytic
//!     This class describes the quantity by an analytic expression this
//!     includes (non-linear) materials, analytic flow fields, etc.
//!   \li \c coefFunctionGrid
//!     This class provides the Information based on information read in
//!     from a grid. E.g. in the case of (interpolated) aeroacoustic source terms 
//!   \li \c coefFunctionSol
//!     This function provides the information based on the current
//!     solution of the problem
//!
//! \note As derived classes of the CoefFunction can be very complex and hide
//!       even whole FeFunctions on other grids, we generally disallow copying
//!       and assigning CoefFunction objects. 
//! 
class CoefFunction : public boost::noncopyable  {
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
    CONSTANT,         /*!< No dependency on space or time */
    TIMEFREQ,      /*!< Only depending on time / frequency, not space */
    GENERAL,       /*!< General dependency (spatial and / or time / freq) */
    SOLUTION       /*!< Dependency on another FeFunction */
  } CoefDependType;
  static Enum<CoefDependType> CoefDependType_;
  
  // ========================
  //  FACTORY METHODS
  // ========================
  //@{ \name Factory Methods

  //! Generate scalar-valued coefficient function from string(s)
  
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
  static PtrCoefFct 
  Generate( MathParser * mp,
            Global::ComplexPart type, 
            const std::string& realVal, 
            const std::string& imagVal = "0" );

  //! Generate vector-valued coefficient function from strings
    
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
  static PtrCoefFct 
  Generate( MathParser * mp,
            Global::ComplexPart type, 
            const StdVector<std::string>& realVal, 
            const StdVector<std::string>& imagVal = StdVector<std::string>() );


  //! Generate vector-valued coefficient function from scalar CoefFunctions
  static PtrCoefFct 
  Generate( MathParser * mp,
            Global::ComplexPart type, 
            const StdVector<PtrCoefFct>& realVal, 
            const StdVector<PtrCoefFct>& imagVal = StdVector<PtrCoefFct>() );


  //! Generate tensor-valued coefficient function from strings

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
  static PtrCoefFct 
  Generate( MathParser * mp,
            Global::ComplexPart type,
            UInt numRows, UInt numCols,
            const StdVector<std::string>& realVal,
            const StdVector<std::string>& imagVal = StdVector<std::string>() );
  
  //! Generate tensor-valued coefficient function from scalar CoefFunctions
  static PtrCoefFct 
  Generate( MathParser * mp,
            Global::ComplexPart type,
            UInt numRows, UInt numCols,
            const StdVector<PtrCoefFct>& realVal,
            const StdVector<PtrCoefFct>& imagVal = StdVector<PtrCoefFct>() );


  //! Generate coefficient function from coefficient expression

  //! This method generates a new coefficient function based on an expression.
  //! The method investigates the expression and tries to generate the most
  //! efficient coefficient representation, i.e. if the expression evaluates
  //! to a constant, a CoefFunctionConst/Analytical/TimeFreq is generated.
  //! In all other cases, a compound coefficient function is generated. 
  static PtrCoefFct Generate( MathParser * mp, 
                              Global::ComplexPart type,
                              const CoefXpr&  xpr );
  //@}

  // ==========================
  //  CONSTRUCTOR / DESTRUCTOR
  // ==========================
  //@{ \name Constructor / Destructor

  //! Constructor
  CoefFunction(){
    dimType_ = NO_DIM;
    dependType_ = CONSTANT;
    isAnalytic_ = false;
    isComplex_ = false;
    
    // by default, the coefficients do not
    // depend on any coordinate system
    this->coordSys_ = NULL;
  }

  //! Destructor
  virtual ~CoefFunction(){
    ;
  }
  
  //! Get an instance of the object with only a given ComplexPart
  
  //! This method allows to retrieve just the real / imaginary part of a 
  //! complex valued coefficient function if applicable. In case the 
  //! real-part of an already real-valued coefficient function is requested,
  //! a pointer to the same instance is returned
  virtual PtrCoefFct GetComplexPart( Global::ComplexPart part ) {
    EXCEPTION( "Method CoefFunction::GetComplexPart not properly overwritten" );
    return PtrCoefFct();
  }
  
  //@}
  
  // ========================
  //  ACCESS METHODS
  // ========================
  //@{ \name Access Methods

  //! Return real-valued tensor at integration point
  virtual void GetTensor(Matrix<Double>& tensor, 
                         const LocPointMapped& lpm ) {
    EXCEPTION( "CoefFunction::GetTensor<Double> called: This may not happen. "
        << "Most likely this method is called with a complex-valued "
        << "CoefFunction object." );
  }

  //! Return real-valued vector at integration point
  virtual void GetVector(Vector<Double>& vec, 
                         const LocPointMapped& lpm ) {
    EXCEPTION( "CoefFunction::GetVector<Double> called: This may not happen " 
        << "Most likely this method is called with a complex-valued "
        << "CoefFunction object." );
  }

  //! Return real-valued scalar at integration point
  virtual void GetScalar(Double& scal, 
                         const LocPointMapped& lpm ) {
    EXCEPTION( "CoefFunction::GetScalar<Double> called: This may not happen. " 
        << "Most likely this method is called with a complex-valued "
        << "CoefFunction object." );
  }

  //! Return complex-valued tensor at integration point
  virtual void GetTensor(Matrix<Complex>& tensor, 
                         const LocPointMapped& lpm ) {
    
    // Provide default implementation in the base class, which returns
    // just the double values as real-valued complex matrix
    Matrix<Double> temp;
    GetTensor( temp, lpm);
    tensor.Resize(temp.GetNumRows(), temp.GetNumCols() );
    tensor.SetPart(Global::REAL, temp);
  }

  //! Return complex-valued vector at integration point
  virtual void GetVector(Vector<Complex>& vec, 
                         const LocPointMapped& lpm ) {
    // Provide default implementation in the base class, which returns
     // just the double values as real-valued complex vector
     Vector<Double> temp;
     GetVector( temp, lpm);
     vec.Resize(temp.GetSize() );
     vec.SetPart(Global::REAL, temp);
  }

  //! Return complex-valued scalar at integration point
  virtual void GetScalar(Complex& scalar, 
                         const LocPointMapped& lpm ) {
    // Provide default implementation in the base class, which returns
    // just the double value as real part of a complex
    Double temp;
    GetScalar( temp, lpm);
     scalar = Complex(temp, 0.0);
  }
  
  //! Return size of vector in case coefficient function is a vector
  virtual UInt GetVecSize() const {
    EXCEPTION( "CoefFunction::GetVecSize: Not overwritten");
    return 0;
  }
  
  //! Return row and columns size of tensor if coefficient function is a tensor
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
      EXCEPTION( "CoefFunction::GetVecSize: Not overwritten");
  }
  
  //@}


  //! Set associated coordinate system
  virtual void SetCoordinateSystem(CoordSystem* cSys){
    coordSys_ = cSys;
  }
  
  //! Get associated coordinate system
  CoordSystem* GetCoordinateSystem(){
    return coordSys_;
  }
  
  //! Return dependency of CoefFunction
  CoefDependType GetDependency() {
    return dependType_;
  }

  //! Return type of entry (scalar, vector, tensor)
  virtual CoefDimType GetDimType() const{
    return dimType_;
  }
  
  //! Return if coeffunction is zero
  
  //! Returns, if the coefficient function is zero. In general, this
  //! is true for all but the CoefFunctionConst class.
  virtual bool IsZero() const {
    return false;
  }
  
  //! Return if coefficient function is analytic
  
  //! This method can be used to query, if a coefficient function can be
  //! described by a closed formula (= string). 
  virtual bool IsAnalytic() const {
    return isAnalytic_;
  }
  
  //! Return if coeffunction is complex
  virtual bool IsComplex() const {
    return isComplex_;
  }

  virtual void AddEntityList(shared_ptr<EntityList> ent){
    if(!entities_.Contains(ent)){
      entities_.Push_back(ent);
    }else{
      WARN("entity list " << ent->GetName() << " already contained in CoefFunction")
    }
  }

  virtual StdVector<shared_ptr<EntityList> > GetEntityList(){
    return entities_;
  }

  //! Dump coefficient function to string 
  virtual std::string ToString() const {
    EXCEPTION("CoefFuncion: ToString() not properly overwritten");
    return "";
  }
  
  // ======================================================================
  //  Helper methods for generating variable names of coefficient function
  // ======================================================================
  //@{ \name Generate variable names for a given coefficient function
  
  //! Generate vector of variable names for all components of a vector
  static void GenScalCompNames( std::string& realVar, 
                                std::string& imagVar, 
                                const std::string& prefix,
                                PtrCoefFct cf );
  
  //! Generate vector of variable names for all components of a vector
  static void GenVecCompNames( StdVector<std::string>& realVar, 
                               StdVector<std::string>& imagVar, 
                               const std::string& prefix,
                               PtrCoefFct cf );
  
  //! Generate vector of variable names for all components of a tensor
  static void GenTensorCompNames( StdVector<std::string>& realVar,
                                  StdVector<std::string>& imagVar,
                                  const std::string& prefix,
                                  PtrCoefFct cf );

  //@}

  // ======================================================================
  //  Some specialized methods for data from external grids
  //  This bearks to some amount the class interface but is done to improve speed
  // ======================================================================
  //@{ \name External Data interfaces

  //! Map Conservative to FeFunction Vector
  virtual void MapConservative( shared_ptr<FeSpace> targetSpace,
                                    Vector<Double>& feFncVec){
    EXCEPTION("This coefficient function does not support Conservative Mapping");
  }

  //! Map Conservative to FeFunction Vector
  virtual void MapConservative( shared_ptr<FeSpace> targetSpace,
                                    Vector<Complex>& feFncVec){
    EXCEPTION("This coefficient function does not support Conservative Mapping");
  }

  //! Determine if coefFunction has conservative mapping
  virtual bool IsConservative(){
    return false;
  }

  //! Set if coefFunction has conservative mapping usually done by PDE
  //! for Coeffunction grid, this overrides the user settings
  virtual void SetConservative(bool value){
    return;
  }


  //! Give Values at global coordinate locations
  virtual void GetVectorValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                             StdVector< Vector<Double> >& values){
    EXCEPTION("This coefficient function does not support evaluation at global coordinates");
  }

  //! Give Values at global coordinate locations
  virtual void GetVectorValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                             StdVector< Vector<Complex> >& values){
    EXCEPTION("This coefficient function does not support evaluation at global coordinates");
  }

  //! Give Values at global coordinate locations
  virtual void GetScalarValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                             StdVector< Double >& values){
    EXCEPTION("This coefficient function does not support evaluation at global coordinates");
  }

  //! Give Values at global coordinate locations
  virtual void GetScalarValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                             StdVector< Complex >& values){
    EXCEPTION("This coefficient function does not support evaluation at global coordinates");
  }


  //@}
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
  // i.e. change the domain to hold a shared_ptr to the coordinate systems!
  CoordSystem* coordSys_;

  //! Dimension of coefficient function (scalar, vector, tensor)
  CoefDimType dimType_;
  
  //! Dependency type of the coefficient function
  CoefDependType dependType_;
  
  //! Flag, if coefficient function is analytic (= can be represented as string)
  bool isAnalytic_;
  
  //! Flag, if coefficient function is complex-valued
  bool isComplex_;

  //! Support of the CoefFunction. Only needed for grid/solution results
  StdVector<shared_ptr<EntityList> > entities_;
};


//! Subclass for analytical coefficient functions

//! This class represents all analytical coefficient functions,
//! which can be represented by a formula (constant, analytical expression).
//! Thus they can deliver a string representation of their scalar, vector
//! or tensor value.
class CoefFunctionAnalytic : public CoefFunction {
public:

  //! Constructor
  CoefFunctionAnalytic()
  : CoefFunction () {
    isAnalytic_ = true;
  }

  //! Destructor
  virtual ~CoefFunctionAnalytic() {}


  //! Get scalar expression
  virtual void GetStrScalar( std::string& real, std::string& imag ) {
    EXCEPTION( "Not implemented here");
  }

  //! Get vector expression
  virtual void GetStrVector( StdVector<std::string>& real, 
                             StdVector<std::string>& imag ) {
    EXCEPTION( "Not implemented here");
  }

  //! Get tensor expression
  virtual void GetStrTensor( UInt& numRows, UInt& numCols,
                             StdVector<std::string>& real, 
                             StdVector<std::string>& imag ) {
    EXCEPTION( "Not implemented here");
  }

};

}
#endif
