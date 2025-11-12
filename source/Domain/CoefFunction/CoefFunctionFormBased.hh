// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef COEF_FUNCTION_FORM_BASED_HH
#define COEF_FUNCTION_FORM_BASED_HH

#include "CoefFunction.hh"
#include "FeBasis/BaseFE.hh"
#include "FeBasis/FeFunctions.hh"
#include "Forms/Operators/BaseBOperator.hh"
#include "Optimization/Design/DesignMaterial.hh"
#include "Optimization/Context.hh"

namespace CoupledField  {

//! Forward class declarations
class BaseBDBInt;
class FeSpace;

// ==========================================================================
//  FORM BASED COEFFICIENT FUNCTION
// ==========================================================================

//! Base class for coefficient functions based on (Bi)Linearforms

//! This class serves as base for all coefficient functions, where the value
//! is calculated using a bilinearform. Examples are the flux quantity,
//! the energy density or the derivative of a potential function.
//! Thus this class mainly serves as an interface to store bilinearforms
//! associated with regions.
class CoefFunctionFormBased : public CoefFunction {
public:
  
  //! Constructor
  CoefFunctionFormBased(  );
  
  //! Destructor
  virtual ~CoefFunctionFormBased();
  
  virtual string GetName() const { return "CoefFunctionFormBased"; }

  //! Set integrator for specific region
  virtual void AddIntegrator(BaseBDBInt* form, RegionIdType region);

  //! Set name request for specific integrator
  void SetIntegratorName(const std::string& integratorName) {
    if( !integratorName_.empty() ) {
      EXCEPTION("Implementation error: Overwriting already defined integrator " << integratorName << " with integrator " << integratorName_ << " for post-processing result makes no sense. Please fix the implementation by using SetIntegratorName() correctly (only once!).");
    }
    integratorName_ = integratorName;
  }
  
  //! Return type of entry (scalar, vector, tensor)
  virtual CoefDimType GetDimType() const { return dimType_;  }

  std::map<RegionIdType, BaseBDBInt* > GetForms() const { return forms_.Mine(); }

protected:

  //! Store bilinearform for each region
  CfsTLS< std::map<RegionIdType, BaseBDBInt* > > forms_;

  //! Store the name of the requestet integrator type
  std::string integratorName_;

};


// ==========================================================================
//  COEFFICIENT FUNCTION BASED ON DIFFERENTIAL OPERATOR
// ==========================================================================
//! Coefficient function, wrapping the application of a differential operator  

//! This class computes the spatial derivative of the primary unknown by
//! applying the B-Operator of the related BDB-class (e..g the grad or curl
//! of a scalar / vector potential). The BDB-bilinearforms /
//! B-Operators have to get passed to this class for every region the result 
//! might get calculated at.
template<class TYPE>
class CoefFunctionBOp : public CoefFunctionFormBased {
public:

  //! Constructor 
  CoefFunctionBOp( shared_ptr<BaseFeFunction> feFct,
                   shared_ptr<ResultInfo> info,
                   TYPE factor = 1.0 );
  
  virtual string GetName() const { return "CoefFunctionBOp"; }


  //! Set integrator for specific region
  virtual void AddIntegrator( BaseBDBInt* form,  
                              RegionIdType region );
  
  //! Pass directly a B-operator
  void AddBOperator( BaseBOperator* bOp,
                     RegionIdType region,
                     std::string integratorName="" );

  //! Destructor
  virtual ~CoefFunctionBOp();
  
  // ========================
  //  ACCESS METHODS
  // ========================
  //@{ \name Access Methods

  //! \copydoc CoefFunction::GetTensor
  virtual void GetTensor(Matrix<TYPE>& coefMat,
                           const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetVector
  virtual void GetVector(Vector<TYPE>& coefVec,
                         const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetScalar
  virtual void GetScalar(TYPE& coefScal,
                           const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetVecSize
  virtual UInt GetVecSize() const;

  //! \copydoc CoefFunction::GetTensorSize
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
    EXCEPTION("This class defined coefficients of vector type only." );
  }
  //@}

  //! \copydoc CoefFunction::ToString
  virtual std::string ToString() const;
  
protected:

  //! Differential operator for each region (not thread relevant)
  std::map<RegionIdType, BaseBOperator* > bOps_;
 
  //! FeFunction containing the coefficients
  shared_ptr<FeFunction<TYPE> > feFct_;

  //! Result info object of result to be calculated (not thread relevant)
  shared_ptr<ResultInfo> res_;
  
  //! Point to FeSpace (not thread relevant)
  shared_ptr<FeSpace> feSpace_;
  
  //! Solution of element (thread relevant)
  CfsTLS< Vector<TYPE> > elemSol_;
  
  //! Operator matrix
  CfsTLS< Matrix<TYPE> > bMat_;
  
  //! Additional factor
  TYPE factor_;
};

// ==========================================================================
//  COEFFICIENT FUNCTION BASED ON FLUX COMPUTATION OF A BDB-FORM
// ==========================================================================
//! Coefficient function, defined by the flux of a FeFunction

//! This class computes the flux of the primary unknown by
//! applying the dB-Operator of the related BDB-class. The BDB-bilinearforms
//! have to get passed to this class for every region the result might get
//! calculated at.
//! \tparam TYPE Real-valued or Complex valued instantiation
//! \tparam TRANS Apply transposed version of dB / Ad matrix
template<class TYPE, bool TRANS = false >
class CoefFunctionFlux : public CoefFunctionFormBased {
public:

  //! Constructor
  CoefFunctionFlux( shared_ptr<BaseFeFunction> feFct,
                    shared_ptr<ResultInfo> info,
                    TYPE factor = 1.0 );
  //! Destructor
  virtual ~CoefFunctionFlux() {}; 

  virtual string GetName() const { return "CoefFunctionFlux"; }

  //! Set integrator for specific region
  virtual void AddIntegrator( BaseBDBInt* form,  
                              RegionIdType region );

  // ========================
  //  ACCESS METHODS
  // ========================
  //@{ \name Access Methods

  //! \copydoc CoefFunction::GetTensor
  virtual void GetTensor(Matrix<TYPE>& coefMat,
                           const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetVector
  virtual void GetVector(Vector<TYPE>& coefVec,
                         const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetScalar
  virtual void GetScalar(TYPE& coefScal,
                           const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetVecSize
  virtual UInt GetVecSize() const;

  //! \copydoc CoefFunction::GetTensorSize
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
    EXCEPTION("This class defined coefficients of vector type only." );
  }
  //@}

  //! \copydoc CoefFunction::ToString
  virtual std::string ToString() const;

protected:

  //! FeFunction containing the coefficients
  shared_ptr<FeFunction<TYPE> > feFct_;

  //! Pointer to result info of desired result
  shared_ptr<ResultInfo> res_;
  
  //! Solution of element
  CfsTLS< Vector<TYPE> > elemSol_;

  //! Additional factor
  TYPE factor_;
};


// ==========================================================================
//  COEFFICIENT FUNCTION EIGEN
// ==========================================================================
//! Coefficient function that calculates the eigenvalues of another CoefFunction.

//! This CoefFunction was designed for the purpose of calculating principal stress and principal strain
//! Only works if the handled CoefFct (e.g. CoefFctFlux for stresses) calls a GetVector(), other methods are not implemented.

class CoefFunctionEigen : public CoefFunctionFormBased {
public:

  //! Constructor
  CoefFunctionEigen( shared_ptr<BaseFeFunction> feFct,
                    shared_ptr<ResultInfo> info,
                    PtrCoefFct stressCoef,
                    Double factor = 1.0 );
  //! Destructor
  virtual ~CoefFunctionEigen();

  virtual string GetName() const { return "CoefFunctionEigen"; }


  // ========================
  //  ACCESS METHODS
  // ========================
  //@{ \name Access Methods

  //! \copydoc CoefFunction::GetVector
  virtual void GetVector(Vector<Double>& coefVec, const LocPointMapped& lpm);

  //! \copydoc CoefFunction::GetVecSize
  virtual UInt GetVecSize() const;

  //! Calculates Eigenvector and Eigenvalue from a stress- or strain coefVec. For principal stresses and strain.
  void GetEigenFromCoefVec(Vector<Double> &solVec);

  //! \copydoc CoefFunction::GetTensorSize
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
    EXCEPTION("This class defined coefficients of vector type only." );
  }

  //! \copydoc CoefFunction::ToString
  virtual std::string ToString() const;

protected:

  //! FeFunction containing the coefficients
  shared_ptr<FeFunction<Double> > feFct_;

  //! Pointer to result info of desired result
  shared_ptr<ResultInfo> res_;

  //! Additional factor
  Double factor_;

  //! Coefficient Function of derived solution (e.g. strain or stress coefficients)
  PtrCoefFct stressCoef_;

};

// ==========================================================================
//  COEFFICIENT FUNCTION BASED ON KERNEL OF BDB-INTEGRATOR
// ==========================================================================
//! Coefficient function, calculated by the kernel of a BdB-integrator

//! This coefficient functions is computed by applying the kernel of a
//! BdB-integrator to the element solution. This can be used e.g. to calculate
//! the energy density.
template<class TYPE>
class CoefFunctionBdBKernel : public CoefFunctionFormBased
{
public:

  //! Constructor
  CoefFunctionBdBKernel(shared_ptr<BaseFeFunction> feFct,  TYPE factor = 1.0);
  //! Destructor
  virtual ~CoefFunctionBdBKernel();

  virtual string GetName() const { return "CoefFunctionBdBKernel"; }

  // ========================
  //  ACCESS METHODS
  // ========================
  //@{ \name Access Methods

  //! \copydoc CoefFunction::GetScalar
  virtual void GetScalar( TYPE& coefScal, const LocPointMapped& lpm);

  //! \copydoc CoefFunction::GetVecSize
  virtual UInt GetVecSize() const {
    EXCEPTION("This class defines coefficients of scalar type only.");
    return 0;
  }

  //! \copydoc CoefFunction::GetTensorSize
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const {
    EXCEPTION("This class defines coefficients of scalar type only.");
  }

  //@}

  //! \copydoc CoefFunction::ToString
  virtual std::string ToString() const;

protected:

  //! FeFunction containing the coefficients
  shared_ptr<FeFunction<TYPE> > feFct_;

  //! Pointer to result info of desired result
  shared_ptr<ResultInfo> res_;
  
  //! Kernel of element matrix
  CfsTLS< Matrix<TYPE> > kernel_;
  
  //! Kernel of element matrix (always real-valued)
  CfsTLS< Matrix<Double> > kernelR_;
  
  //! Solution of element
  CfsTLS< Vector<TYPE> > elemSol_;

  //! Additional factor
  TYPE factor_;

};


/** Calculates the dyadic product of strain vs. strain
 * This is required for external topology gradient evaluation for Bloch mode analysis (Nazarov).
 * The class is a modification of CoefFunctionBdBKernel */
template<class TYPE>
class CoefFunctionDyadicStrain : public CoefFunctionFormBased
{
public:
  CoefFunctionDyadicStrain(shared_ptr<BaseFeFunction> feFct);

  virtual ~CoefFunctionDyadicStrain();

  virtual string GetName() const { return "CoefFunctionDyadicStrain"; }

  //! \copydoc CoefFunction::GetTensorSize
  virtual void GetTensorSize(unsigned int& numRows, unsigned int& numCols ) const;

  virtual void GetTensor(Matrix<TYPE>& tensor, const LocPointMapped& lpm);

  //! \copydoc CoefFunction::ToString
  virtual std::string ToString() const;

protected:

  //! FeFunction containing the coefficients
  shared_ptr<FeFunction<TYPE> > feFct_;

  //! Solution of element
  CfsTLS< Vector<TYPE> > elemSol_;
};


/** Calculates the dyadic product of strain vs. strain
 * This is required for external topology gradient evaluation for Bloch mode analysis (Nazarov).
 * The class is a modification of CoefFunctionBdBKernel */
template<class TYPE>
class CoefFunctionQuadSol : public CoefFunctionFormBased
{
public:
  CoefFunctionQuadSol(shared_ptr<BaseFeFunction> feFct);

  virtual ~CoefFunctionQuadSol();

  virtual string GetName() const { return "CoefFunctionQuadSol"; }

  virtual void GetScalar(TYPE& coefScal, const LocPointMapped& lpm);

  //! \copydoc CoefFunction::ToString
  virtual std::string ToString() const;

protected:

  //! FeFunction containing the coefficients
  shared_ptr<FeFunction<TYPE> > feFct_;

  //! Solution of element
  CfsTLS< Vector<TYPE> > elemSol_;

};

/**
 * Calculates scalar values (e.g. pressure) from lattice Boltzmann(LBM) particle distribution function values
 */
template<class TYPE> 
class CoefFunctionLBM : public CoefFunctionFormBased
{
public:
  CoefFunctionLBM(LatticeBoltzmannPDE* lbm, shared_ptr<BaseFeFunction> feFct,shared_ptr<ResultInfo> resInfo);

  virtual ~CoefFunctionLBM();

  virtual string GetName() const { return "CoefFunctionLBM"; }


  virtual void GetScalar(TYPE& coefScal, const LocPointMapped& lpm);

  virtual void GetVector(Vector<TYPE>& vec, const LocPointMapped& lpm);

  //! \copydoc CoefFunction::ToString
  virtual std::string ToString() const;

protected:

  //! FeFunction containing the coefficients
//  shared_ptr<FeFunction<TYPE> > feFct_;

  //! Result name
  SolutionType resType_;

  //! Pointer to LBM object. We need this to call functions like CalcDensity() etc.
  LatticeBoltzmannPDE* lbm_;
};

/** Simply returns the homogenized tensor. Does not work by using the existing CoefFunction :( */
template<class TYPE, App::Type APP>
class CoefFunctionHomogenization : public CoefFunctionFormBased
{
public:
  CoefFunctionHomogenization(shared_ptr<BaseFeFunction> feFct, MaterialTensorNotation notation = NO_NOTATION);

  virtual ~CoefFunctionHomogenization();

  virtual string GetName() const { return "CoefFunctionHomogenization"; }

  //! \copydoc CoefFunction::GetTensorSize
  unsigned int GetVecSize() const;

  void GetVector(Vector<TYPE>& vec, const LocPointMapped& lpm);

  void GetTensorSize(unsigned int& numRows, unsigned int& numCols ) const;

  void GetTensor(Matrix<TYPE>& tensor, const LocPointMapped& lpm);

  void GetScalar(TYPE& coefScal, const LocPointMapped& lpm);

  //! \copydoc CoefFunction::ToString
  virtual std::string ToString() const;

private:

  //! FeFunction containing the coefficients
  shared_ptr<FeFunction<TYPE> > feFct_;

  MaterialTensorNotation notation_;
};



} // end of namespace
#endif
