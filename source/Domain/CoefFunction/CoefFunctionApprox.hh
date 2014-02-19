#ifndef FILE_COEFFUNCTION_APPROX_HH
#define FILE_COEFFUNCTION_APPROX_HH

#include <boost/tr1/type_traits.hpp>
#include "CoefFunction.hh"
#include "FeBasis/FeFunctions.hh"

namespace CoupledField {

// forward clas declaration
class ApproxData;
class BaseBOperator;
class Grid;
class FeFunctions;
//class LocalPointMapped;


// ============================================================================
//  ISOTROPIC VERSIONS
// ============================================================================
//! Provide a coefficient for approximated sample data (scalar)

//! This class encapsulates the approximation of sampled data as used e.g. 
//! for nonlinear methods. It internally
//! utilizes the ApproxData class for approximating nonlinear curves.
//! \note This class only works for real-valued scalar data.
class CoefFunctionApprox : public CoefFunction{
public:

  //! Constructor
  CoefFunctionApprox();

  //! Destructor
  virtual ~CoefFunctionApprox();
  
  //! Initialize with data
  void Init( Double coefScalar, ApproxData * nLinFnc,
             PtrCoefFct dependCoef  );

  
  //! \see CoefFunction::GetScalar
  void GetScalar(Double& coefScalar, 
                 const LocPointMapped& lpm );
    

  //! \see CoefFunction::ToString
  std::string ToString() const;

protected:
  
  //! Constant initial value of the curve
  Double coefScalar_;
  
  //! Class for function approximation
  ApproxData * nLinFnc_;
  
  //! Coefficient function which this one depends one
  PtrCoefFct dependCoef_;
};

// ============================================================================
//  Coef Function Composite
// ============================================================================
//! Provide a coefficient for approximated sample data (scalar)

//! This class works similar to the CoefFunctionApprox 
//! It has the following additional features:
//!  - it can hold `n` dependant coef functions. e.g. when a non-linearity depends on two solutions (e.g. heat and elec potential)
//!  - it allows the definition of regions (called terminals) where the value of a given dependency is evaluated
//! 
//! \note This class only works for real-valued scalar data.

class CoefFunctionComposite : public CoefFunction{
public:

  //! Constructor
  CoefFunctionComposite();

  //! Destructor
  virtual ~CoefFunctionComposite() {}
  
  //! \see CoefFunction::GetScalar
  virtual void GetScalar(Double& coefScalar, 
                 const LocPointMapped& lpm ) {EXCEPTION("not implemented");}
    
  void SetRegion(TerminalConnector tc, RegionIdType reg);
  void SetDependCoef(NonLinType nl, PtrCoefFct dep );

  //! Initialize with data
  void Init( Double coefScalar, ApproxData * nLinFnc  );

  //! Helper function to return the average value at the defined terminal/region \param tc
  //! for the NonLinType \param nl for a point \param lpm
  Double GetAvgTerminalValue(TerminalConnector tc, NonLinType nl, const LocPointMapped & lpm);

  //! \see CoefFunction::ToString
  std::string ToString() const {return "Composite Coefficient function";}

protected:

  //! region Ids of the terminals
  std::map<TerminalConnector, RegionIdType> terminals_;

  //! number of regions
  UInt nRegions_;

  //! Class for function approximation
  ApproxData * nLinFnc_;
  
  //! Constant initial value of the curve
  Double coefScalar_;
  
  //! Coefficient function which this one depends one
  std::map<NonLinType, PtrCoefFct> dependCoefs_;

  //! number of depending coef fcts
  UInt nDepCoefs_;
};

// ============================================================================

// ============================================================================
//  Coef Function Bipole
// ============================================================================
//! Provide an approximated value for a bipolar (two terminals) which depends on
//! the difference of the dependant coef function on the two terminals

class CoefFunctionBipole : public CoefFunctionComposite {
public:
  //! Constructor
  CoefFunctionBipole() {}

  //! Destructor
  virtual ~CoefFunctionBipole() {}

  //! \see CoefFunction::GetScalar
  void GetScalar(Double& coefScalar, 
                 const LocPointMapped& lpm );
    
  //! \see CoefFunction::ToString
  std::string ToString() const {return "Composite Coefficient function for temperature dependent bipole";}
};

// ============================================================================

// ============================================================================
//  Coef Function Heat Bipole
// ============================================================================
//! Provide an approximated value for a bipolar (two terminals) which depends on
//! the difference of the dependant coef function on the two terminals
//! _and_ on the temperature (second dependant coef function)
//! This class requires that the nonlinear function depends on two variables. 
//! (Currently only supported by BiLinApproximate)
//! 

class CoefFunctionHeatBipole : public CoefFunctionComposite {
public:
  //! Constructor
  CoefFunctionHeatBipole() {}

  //! Destructor
  virtual ~CoefFunctionHeatBipole() {}

  //! \see CoefFunction::GetScalar
  void GetScalar(Double& coefScalar, 
                 const LocPointMapped& lpm );
    
  //! \see CoefFunction::ToString
  std::string ToString() const {return "Composite Coefficient function for temperature dependent bipole";}
};


// ============================================================================
//! Provide a coefficient for approximated derivative of sample data (scalar)

//! This class encapsulates the approximation of the derivative of sampled 
//! data, e.g. as needed for the nonlinear Newton methods. It internally
//! utilizes the ApproxData class for approximating nonlinear curves.
//! \note This class only works for real-valued scalar data.
class CoefFunctionApproxDeriv : public CoefFunction{
public:

  //! Constructor
  CoefFunctionApproxDeriv();

  //! Destructor
  virtual ~CoefFunctionApproxDeriv();
  
  //! Initialize with data
  void Init( ApproxData * nLinFnc,
             UInt dimDMat,
             PtrCoefFct dependCoef );

  //! \see CoefFunction::GetTensor
  void GetTensor(Matrix<Double>& coefMat, 
                 const LocPointMapped& lpm );

  //! \see CoefFunction::ToString
  std::string ToString() const;

protected:
  
  //! Dimension of the D-matrix
  UInt dimDMat_;
  
  //! Class for function approximation
  ApproxData * nLinFnc_;
  
  //! Coefficient function which this one depends one
  PtrCoefFct dependCoef_;
};

// ============================================================================
//  ANISOTROPIC VERSIONS
// ============================================================================
//! Provide a coefficient for approximated anisotropic sampled data (scalar)

//! This class encapsulates the approximation of sampled data as used e.g. 
//! for nonlinear methods. It internally
//! utilizes the ApproxData class for approximating nonlinear curves.
//! In addition it maps the 
//! \note This class only works for real-valued scalar data.
class CoefFunctionApproxAniso : public CoefFunction{
public:

  //! Constructor
  CoefFunctionApproxAniso();

  //! Destructor
  virtual ~CoefFunctionApproxAniso();

  //! Initialize with data
  void Init( Double coefScalar, 
             StdVector<ApproxData*>  nLinFnc,
             StdVector<Double> angles,
             PtrCoefFct dependCoef );

  //! \see CoefFunction::GetScalar
  void GetScalar(Double& coefScalar, 
                 const LocPointMapped& lpm );

  //! \see CoefFunction::ToString
  std::string ToString() const;

protected:
  
  //! Constant initial value of the curve
  Double coefScalar_;
  
  //! Vector containing the approximations of the curves
  StdVector<ApproxData* > nLinFnc_;
  
  //! Vector containing the approximations of the curve
  StdVector<Double> angles_;
  
  //! Coefficient function which this one depends one
  PtrCoefFct dependCoef_;
};

//! Provide a coefficient for approximated anisotropic derivative of sample data (scalar)

//! This class encapsulates the approximation of the derivative of sampled 
//! data, e.g. as needed for the nonlinear Newton methods. It internally
//! utilizes the ApproxData class for approximating nonlinear curves.
//! \note This class only works for real-valued scalar data.
class CoefFunctionApproxDerivAniso : public CoefFunction{
public:

  //! Constructor
  CoefFunctionApproxDerivAniso();

  //! Destructor
  virtual ~CoefFunctionApproxDerivAniso();
  
  //! Initialize with data
  void Init( StdVector<ApproxData*>  nLinFnc,
             StdVector<Double> angles,
             UInt dimDMat,
             PtrCoefFct dependCoef );

  //! \see CoefFunction::GetTensor
  void GetTensor(Matrix<Double>& coefMat, 
                 const LocPointMapped& lpm );

  //! \see CoefFunction::ToString
  std::string ToString() const;

protected:
  
  //! Dimension of the D-matrix
  UInt dimDMat_;
  
  //! Vector containing the approximations of the curves
  StdVector<ApproxData* > nLinFnc_;
  
  //! Vector containing the approximations of the curve
  StdVector<Double> angles_;
  
  //! Coefficient function which this one depends one
  PtrCoefFct dependCoef_;
};

}
#endif
