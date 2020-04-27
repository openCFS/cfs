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

#include "def_use_openmp.hh"
#ifdef USE_OPENMP
#include <omp.h>
#endif

namespace CoupledField{


// forward class declarations
class CoordSystem;
class CoefXpr;
class CoefFunction;
template < typename TYPE > class CoefFunctionConst;
class FeSpace;
class BaseFeFunction;

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
    CONSTANT,      /*!< No dependency on space or time */
    TIMEFREQ,      /*!< Only depending on time / frequency, not space */
    SPACE,         /*!< Only depending on space */
    GENERAL,       /*!< General dependency on space and time /freq */
    SOLUTION       /*!< Dependency on another FeFunction */
  } CoefDependType;
  static Enum<CoefDependType> CoefDependType_;
  
  //! Dependency of coefficient function
  typedef enum{
	NOINVERS,
    INVSOURCE,         /*!< Invserse scheme: source data */
    INVMEASURE         /*!< Inverse scheme: measured data */
  } CoefInverseType;
  static Enum<CoefInverseType> CoefInverseType_;

  //! Dependency of coefficient function
  typedef enum{
	NOINFORMATION,
    FEBASIS,         /*!< Invserse scheme: source data */
    DELTA         /*!< Inverse scheme: measured data */
  } CoefInverseSourceApprox;
  static Enum<CoefInverseSourceApprox> CoefInverseSourceApprox_;

  //! Modifications of coefficient function
  typedef enum{
    NONE,              /*!< Default interpolation of data*/
    VECTOR_DIVERGENCE  /*!< Return divergence of vector valued CoefFuncton when called with getScalar*/
  } CoefDerivativeType;
  static Enum<CoefDerivativeType> CoefDerivativeType_;
  
  //! Get the maximum CoefFunction dependency type
  static CoefDependType GetMaxCoefDependType(CoefDependType typeA, CoefDependType typeB);

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
  
  //! Generate tensor-valued coefficient function from scalar CoefFunctions,
  //! separate for real and imaginary parts
  static PtrCoefFct
  Generate( MathParser * mp,
            Global::ComplexPart type,
            UInt numRows, UInt numCols,
            const StdVector<PtrCoefFct>& realVal,
            const StdVector<PtrCoefFct>& imagVal );

  //! Generate tensor-valued coefficient function from scalar CoefFunctions
  static PtrCoefFct 
  Generate( MathParser * mp,
            Global::ComplexPart type,
            UInt numRows, UInt numCols,
            const StdVector<PtrCoefFct>& scalars );


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
    if(!IsSerialRegion()){
      EXCEPTION("Constructor of CoefFunction is called from parallel region which is not safe.")
    }
    dimType_ = NO_DIM;
    dependType_ = CONSTANT;
    isAnalytic_ = false;
    isComplex_ = false;
    supportDerivative_ = false;
    derivType_ = NONE;
    inverseType_ = NOINVERS;
    
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
  virtual void GetTensor(Matrix<Double>& tensor, const LocPointMapped& lpm ) {
    EXCEPTION( "CoefFunction::GetTensor<Double> called: This may not happen. "
        << "Most likely this method is called with a complex-valued CoefFunction object." );
  }

  //! Return real-valued vector at integration point
  virtual void GetVector(Vector<Double>& vec, const LocPointMapped& lpm ) {
    EXCEPTION( "CoefFunction::GetVector<Double> called: This may not happen "
        << "Most likely this method is called with a complex-valued CoefFunction object." );
  }

  //! Return real-valued element averaged value
  virtual void GetAvgElemValue(Double & vec, const Elem* elem) {
    EXCEPTION( "CoefFunction::GetAvgElemValue<Double> not implemented in base class" );
  }


  //! Return real-valued scalar at integration point
  virtual void GetScalar(Double& scal, const LocPointMapped& lpm ) {
    EXCEPTION( "CoefFunction::GetScalar<Double> called: This may not happen. " 
        << "Most likely this method is called with a complex-valued "
        << "CoefFunction object." );
  }

  //! Return complex-valued tensor at integration point
  virtual void GetTensor(Matrix<Complex>& tensor, const LocPointMapped& lpm ) {
    // Provide default implementation in the base class, which returns
    // just the double values as real-valued complex matrix
    Matrix<Double> temp;
    GetTensor( temp, lpm);
    tensor.Resize(temp.GetNumRows(), temp.GetNumCols() );
    tensor.SetPart(Global::REAL, temp);
  }

  //! Return complex-valued vector at integration point
  virtual void GetVector(Vector<Complex>& vec, const LocPointMapped& lpm ) {
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

  virtual void Init(shared_ptr<BaseFeFunction> feFct,
            shared_ptr<FeSpace> feSpc,
            const StdVector<RegionIdType>& regions,
            std::map<RegionIdType, BaseMaterial*>& materials,
            Grid* ptGrid,
            PtrCoefFct magFluxCoef,
            const UInt& N,
            const UInt& M,
            const Double& baseFreq,
            const UInt& nFFT){
    EXCEPTION("Not implemented here in base class");
    return;
  }

  //! Set associated coordinate system
  virtual void SetCoordinateSystem(CoordSystem* cSys){
    coordSys_ = cSys;
  }
  
  //!
  virtual void SetFeFunction( shared_ptr<BaseFeFunction> fct1, SolutionType solType) {
	  feFunctions_[solType] = fct1;
  }

  //! Get associated coordinate system
  CoordSystem* GetCoordinateSystem(){
    return coordSys_;
  }
  
  //! Return dependency of CoefFunction
  CoefDependType GetDependency() {
    return dependType_;
  }

  //! Return dependency of CoefFunction
  CoefInverseType GetInverseType() {
    return inverseType_;
  }

  //! Set type of approximation for source type
  void SetInverseSourceApproxType( CoefInverseSourceApprox type )  {
	  approxSourceType_ = type;
  }

  //! Return dependency of CoefFunction
  CoefInverseSourceApprox GetInverseSourceApproxType() {
    return approxSourceType_;
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

  /** helper for using simplified CoefFunctionConst interface GetScalar(), ... without lpm.
   * For performance reason don't use shared pointers, hence don't store the pointer!
   * Does NOT return a const version but a CoefFunctionCons cast!
   * To be used like double val = coef->AsConst<double>->GetScalar();
   * @return NULL if the type is wrong (is just a dynamic cast) */
  template<class TYPE>
  CoefFunctionConst<TYPE>* AsConst(bool throw_exception = false);

  //! stes the coefFnc as active (just used for inverse source identififcation)
  virtual void SetActive(bool val) {
    isActive_ = val;
  }


  //! Dump coefficient function to string 
  virtual std::string ToString() const {
    EXCEPTION("CoefFuncion: ToString() not properly overwritten");
    return "";
  }
  
  //! sets the derivative modification to the coefFunction
  virtual void SetDerivativeOperation(CoefDerivativeType type){
    EXCEPTION("CoefFunction: This CoefFunction does not support derivatives");
    return;
  }
  //! sets the derivative modification to the coefFunction
  virtual void SetDerivativeOperation(CoefDerivativeType type, UInt gDim, UInt dDim){
    EXCEPTION("CoefFunction: This CoefFunction does not support derivatives");
    return;
  }
  //! computes the optimality condition
  virtual void ComputeOptCondition(Double& optAmp, Double& optPhase) {
	  EXCEPTION("CoefFuncion::ComputeOptCondition not implemented");
   }

  //! computes the L2 norm of error
  virtual void ComputeDiff2Meas( Double& error ) {
	  EXCEPTION("CoefFuncion::ComputeDiff2Meas not implemented");
  }

  //! set all parameters for inverse scheme
  virtual void SetInverseParam( Double& alpha, Double& beta, Double& rho, Double& qExp,
		                        Double& freq, std::string fileNameMeasdata,
								std::string logLevel, Double& scalVal) {
 	  EXCEPTION("CoefFuncion::SetInverseParam not implemented");
   }

  //! set all parameters for inverse scheme
  virtual void ChangeInverseParam( Double& alpha, Double& beta, Double& rho) {
   	  EXCEPTION("CoefFuncion::ChangeInverseParam not implemented");
  }

  //! update the source values (amplitude and phase)
  virtual void UpdateSource(Double& stepLength, bool lineSearch, bool scaleBack=false) {
	  EXCEPTION("CoefFuncion::UpdateSource not implemented");
  }

  //! computes the L2 norm of error
  virtual void ComputeTikh(Double& funcVal, Double& resSquared) {
	  EXCEPTION("CoefFuncion::ComputeTikh not implemented");
  }

  //! compute square of L2-norm of measured pressure at mic-positions
  virtual void ComputeMeasL2squared( Double& vaL2 ) {
	  EXCEPTION("CoefFuncion::ComputeMeasL2squared not implemented");
  }

  virtual void SetApproxSourceDelta() {
	  EXCEPTION("CoefFuncion::SetApproxSourceDelta not implemented");
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


  //@{  
  //! Return vectorial values at global coordinate locations
  
  //! This method allows to get the values at several global coordinate 
  //! locations at once. This can be very efficient for simple expressions
  //! (constant, analytical) or rather costly (e.g. for element-discerete values
  //! involving a global-local transformation.
  //!
  //! In the base class, the most general approach is implemented, i.e. we
  //! map every coordinate to element local coordinates and call the related
  //! GetScalar() method, which can be rather costly.
  //! In this case it is advisable to override this method in the derived 
  //! CoefFunction class.
  virtual void GetVectorValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                        StdVector< Vector<Double> >& values,
                                        Grid* ptGrid,
                                        const StdVector<shared_ptr<EntityList> >& srcEntities =
                                        StdVector<shared_ptr<EntityList> >() );
  virtual void GetVectorValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                        StdVector< Vector<Complex> >& values, 
                                        Grid* ptGrid,
                                        const StdVector<shared_ptr<EntityList> >& srcEntities =
                                        StdVector<shared_ptr<EntityList> >() );
  //@}

  //@{
  //! Return scalar values at global coordinate locations
  
  //! This method allows to get the values at several global coordinate 
  //! locations at once. This can be very efficient for simple expressions
  //! (constant, analytical) or rather costly (e.g. for element-discerete values
  //! involving a global-local transformation.
  //!
  //! In the base class, the most general approach is implemented, i.e. we
  //! map every coordinate to element local coordinates and call the related
  //! GetScalar() method, which can be rather costly.
  //! In this case it is advisable to override this method in the derived 
  //! CoefFunction class.
  virtual void GetScalarValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                        StdVector< Double >& values, 
                                        Grid* ptGrid,
                                        const StdVector<shared_ptr<EntityList> >& srcEntities =
                                        StdVector<shared_ptr<EntityList> >() );
  virtual void GetScalarValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                        StdVector< Complex >& values, 
                                        Grid* ptGrid,
                                        const StdVector<shared_ptr<EntityList> >& srcEntities =
                                        StdVector<shared_ptr<EntityList> >() );
  //@}

  //@{
  //! Return scalar values at global coordinate locations

  //! This method allows to get the values at several global coordinate
  //! locations at once. This can be very efficient for simple expressions
  //! (constant, analytical) or rather costly (e.g. for element-discerete values
  //! involving a global-local transformation.
  //!
  //! In the base class, the most general approach is implemented, i.e. we
  //! map every coordinate to element local coordinates and call the related
  //! GetScalar() method, which can be rather costly.
  //! In this case it is advisable to override this method in the derived
  //! CoefFunction class.
  virtual void GetTensorValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                        StdVector< Matrix<Double> >& values,
                                        Grid* ptGrid,
                                        const StdVector<shared_ptr<EntityList> >& srcEntities =
                                        StdVector<shared_ptr<EntityList> >() ) {
    Exception("GetTensorValuesAtCoords<Double> not implemented in base class");
  }

  virtual void GetTensorValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                        StdVector< Matrix<Complex> >& values,
                                        Grid* ptGrid,
                                        const StdVector<shared_ptr<EntityList> >& srcEntities =
                                        StdVector<shared_ptr<EntityList> >() ) {
    Exception("GetTensorValuesAtCoords<Complex> not implemented in base class");
  }

  //! Needed for harmonic balancing CoefFunctionHarmBalance
  virtual shared_ptr<CoefFunction> GenerateMatCoefFnc(const UInt& iRegion,
                                                      const std::string& name,
                                                      const bool nonLin,
                                                      shared_ptr<ElemList> actSDList){
    EXCEPTION("Not implemented in base class");
  }

  //! Needed for harmonic balancing CoefFunctionHarmBalance
  virtual void RegisterElemsInRegion(shared_ptr<ElemList> actSDList,
                                     const UInt& iRegion){
    EXCEPTION("Not implemented in base class");
  }


  //! Functions needed for Hystersis
  virtual void GetCouplTensorSize(UInt& numRows, UInt& numCols){
    EXCEPTION("GetCouplTensorSize not available");
  }

  virtual void ComputeVector(Vector<Double>& outputVector,const LocPointMapped& lpm, int timeLevel, int baseSign, std::string vectorName, bool onBoundary, bool usedAsRHSload ){
    EXCEPTION("ComputeVector not available");
  }

  virtual void ComputeTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm,
          std::string tensorName, std::string implementationVersion, bool transposed, bool rotate, bool useAbs, bool lockPrecomputationAndDeltaMat ){
    EXCEPTION("ComputeTensor not available");
  }

  virtual void PrecomputeMaterialTensorForInverison(){
    EXCEPTION("PrecomputeMaterialTensorForInverison not available");
  }

  virtual void SetPreviousHystVals(bool lastTS = false, bool forceMemoryLock = false){
	  EXCEPTION("SetPreviousHystVals not available");
  }

  virtual void TestInversion(PtrParamNode testNode){
    EXCEPTION( "Not implemented in base class");
  }
  virtual void SetFlag(std::string flagName, Integer intState){
    EXCEPTION( "Not implemented in base class");
  }
  virtual void SetDoubleFlag(std::string flagName, Double intState){
    EXCEPTION( "Not implemented in base class");
  }
  
  virtual bool useStrainForm(){
    EXCEPTION("Not implemented in base class");
  }
  
  virtual int GetStrainForm(){
    EXCEPTION("Not implemented in base class");
  }

  virtual void SetStrainForm(int intState){
    EXCEPTION("Not implemented in base class");
  }

  virtual std::string getPDEName(){
    EXCEPTION("Not implemented in base class");
  }
  
  virtual bool deltaMatActive(){
    EXCEPTION("Not implemented in base class");
  }
  
  virtual bool couplingTensorSet(){
    EXCEPTION("Not implemented in base class");
  }
  virtual int GetDeltaForm(){
    EXCEPTION("Not implemented in base class");
  }
  
  virtual int GetTimeLevel(std::string EntitiyType){
    EXCEPTION("Not implemented in base class");
  }
  
  virtual shared_ptr<CoefFunction> GenerateMatCoefFnc(std::string tensorName ){
    EXCEPTION("Not implemented in base class");
  }
  virtual shared_ptr<CoefFunction> GenerateRHSCoefFnc(std::string vectorName, bool onBoundary = false){
    EXCEPTION("Not implemented in base class");
  }
  virtual shared_ptr<CoefFunction> GenerateRHSCoefFnc(std::string vectorName, shared_ptr<CoefFunction> coefFunctionToBeIncluded){
    EXCEPTION("Not implemented in base class");
  }

  virtual shared_ptr<CoefFunction> GenerateOutputCoefFnc(std::string ResultName){
    EXCEPTION("Not implemented in base class");
  }
  
  virtual void ActiveOneShotSlopeEstimation(Double steppingLength, Double scaling){
    EXCEPTION( "Not implemented in base class");
  }

  virtual void checkSaturationStateAllElements(Double& lastTSSatAvg, Double& lastItSatAvg, Double& curItSatAvg,
      Double& oppositeDirAsTSAvg, Double& oppositeDirAsItAvg){
    EXCEPTION( "Not implemented in base class");
  }

  virtual Vector<Double> GetIrreversibleStrains(const LocPointMapped& Originallpm, int timeLevel){
    EXCEPTION( "Not implemented in base class");
  }

  virtual Matrix<Double> ConvertFromVoigtToTensor(Vector<Double> Si_voigt){
    EXCEPTION( "Not implemented in base class");
  }

  virtual Matrix<Double> GetIrreversibleStrainTensor(const LocPointMapped& Originallpm, int timeLevel){
    EXCEPTION( "Not implemented in base class");
  }

  virtual Vector<Double> GetPrecomputedOutputOfHysteresisOperator(const LocPointMapped& lpm, int timeLevel, bool forStrain){
    EXCEPTION( "Not implemented in base class");
  }

  virtual Vector<Double> GetPrecomputedInputToHysteresisOperator(const LocPointMapped& lpm, int timeLevel){
    EXCEPTION( "Not implemented in base class");
  }

  virtual void SetFPMaterialTensors(Integer intState){
    EXCEPTION( "Not implemented in base class");
  }
  virtual UInt GetFPMaterialState(){
    EXCEPTION( "Not implemented in base class");
  }
  virtual UInt GetFPMaterialStateRHS(){
    EXCEPTION( "Not implemented in base class");
  }
  virtual Matrix<Double> GetFPMaterialTensor(const LocPointMapped& OriginalLPM){
    EXCEPTION( "Not implemented in base class");
  }
  virtual Vector<Double> GetFPCorrectionVector(const LocPointMapped& OriginalLPM, Integer timeLevel){
    EXCEPTION( "Not implemented in base class");
  }

  virtual bool anyMatrixForLocalInversionRequiresComputation(){
    EXCEPTION( "Not implemented in base class");
  }

  virtual void getMatrixForLocalInversion(const LocPointMapped& Originallpm, Matrix<Double>& matrixForInversion, Matrix<Double>& matrixForInversionInverse){
    EXCEPTION( "Not implemented in base class");
  }

  virtual void setMatrixForLocalInversion(Matrix<Double> matrixForInversion, Matrix<Double> matrixForInversionInverse, UInt storageIdx, bool reuse){
    EXCEPTION( "Not implemented in base class");
  }

  virtual void getLPMMaps(std::map<UInt, LocPointMapped >& allLPM, std::map<UInt, LocPointMapped >& midpointLPM){
    EXCEPTION( "Not implemented in base class");
  }

  virtual void GetScaledAndRotatedCouplingTensor(const LocPointMapped& lpm, Matrix<Double>& couplTensor, Matrix<Double>& rotatedCouplTensor, int timeLevel,
  bool rotate=true){
    EXCEPTION( "Not implemented in base class");
  }
  
  virtual Matrix<Double> GetDeltaMat(const LocPointMapped& Originallpm, int timelevel_new, int timelevel_old, bool useStrains, bool useAbs,
      std::string implementationVersion){
    EXCEPTION( "Not implemented in base class");
  }
  
  virtual void SetElastAndCouplTensor(PtrCoefFct elastTensor, PtrCoefFct couplTensor){
    EXCEPTION( "Not implemented in base class");
  }
  
  virtual void AddAdditionalSDList(shared_ptr<EntityList> actSDList, RegionIdType volReg, bool onSurface){
    EXCEPTION( "Not implemented in base class");
  }
  
  virtual Double GetOutputSaturation(){
    EXCEPTION( "Not implemented in base class");
  }


  //! set volume regionId being the correct neighbor of a surface region id
  virtual void SetVolNeighborRegionId(RegionIdType surfId, RegionIdType volId) {
	  neighborRegionId_[surfId] = volId;
  }

  //! return volume regionId being the correct neighbor of a surface region id
  virtual RegionIdType GetVolNeighborRegionId(RegionIdType surfId) {
	  return neighborRegionId_[surfId];
  }
  //@}

  //@}

  // ========================
  //  HELPER METHODS
  // ========================
  //@{ \name Helper methods
  
  //! in case of parallel execution we can ask if we are in parallel
  inline bool IsSerialRegion(){
#ifdef USE_OPENMP
    return (omp_get_num_threads()==1);
#else
    return true;
#endif
  }

  //! Returns true, if expression depends on time / freq
  static bool ExprDependsOnTimeFreq(MathParser* mp, const std::string& expr);
  
  //! Returns true, if one of the given expressions depends on time / freq
  static bool ExprDependsOnTimeFreq(MathParser* mp, const StdVector<std::string>& expr);
  
  //! Returns true, if one of the given expressions depends on space
  static bool ExprDependsOnSpace(MathParser* mp, const std::string& expr);
  
  //! Returns true, if expression depends on space
  static bool ExprDependsOnSpace(MathParser* mp, const StdVector<std::string>& expr);

protected:

  //! Rotates a vector from the local to the global coordinate system
  template<typename TYPE>
  void TransformVectorByCoordSys(Vector<TYPE> &outVec,
                              const Vector<TYPE> &inVec,
                              const Vector<Double> &point);

  //! Rotates a Vector from the local to the global coordinate system
  template<typename TYPE>
  void TransformVectorByCoordSys(Vector<TYPE> &outVec,
                              const Vector<TYPE> &inVec,
                              const LocPointMapped &lpm);

  //! Rotates a tensor from the local to the global coordinate system
  template<typename TYPE>
  void TransformTensorByCoordSys(Matrix<TYPE> &outMat,
                              const Matrix<TYPE> &inMat,
                              const Vector<Double> &point);

  //! Rotates a tensor from the local to the global coordinate system
  template<typename TYPE>
  void TransformTensorByCoordSys(Matrix<TYPE> &outMat,
                              const Matrix<TYPE> &inMat,
                              const LocPointMapped &lpm);

  //@}

  //TODO: CHANGE THIS TO SHARED POINTER
  // i.e. change the domain to hold a shared_ptr to the coordinate systems!
  CoordSystem* coordSys_;

  //! Dimension of coefficient function (scalar, vector, tensor)
  CoefDimType dimType_;
  
  //! Dependency type of the coefficient function
  CoefDependType dependType_;
  
  //! storing the derivative type of the CoefFunction
  CoefDerivativeType derivType_;

  //! storing the type for inverse scheme
  CoefInverseType inverseType_;

  //! how the source term is approximated
  CoefInverseSourceApprox inverseApproxType_;

  //! Flag, if coefficient function is analytic (= can be represented as string)
  bool isAnalytic_;
  
  //! Flag, if coefficient function is complex-valued
  bool isComplex_;

  //! Flag indicating if the CoefFunction supports derivatives
  bool supportDerivative_;

  //! only needed for hystersis
  StdPDE* linkedPDE_;

  //! for each surface Region id we have to correct neighbor volume region id
  std::map<RegionIdType, RegionIdType> neighborRegionId_;

  //! Map Storing FeSpaces for each unknown of PDE
  std::map<SolutionType, shared_ptr<BaseFeFunction> > feFunctions_;

  //! sets the rhsFnc active
  bool isActive_;

  //! approximate source terms with delta functions
  CoefInverseSourceApprox approxSourceType_;
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

  //! \copydoc CoefFunction::SetDerivativeOperation
  virtual void SetDerivativeOperation(CoefDerivativeType type){
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
        //PAY ATTENTION: In case of a derivative, the coefFunction is
        // no longer analytic due to the current implementation!
        this->isAnalytic_ = false;
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
};

}
#endif
