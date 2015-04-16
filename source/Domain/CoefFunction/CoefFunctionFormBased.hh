// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef COEF_FUNCTION_FORM_BASED_HH
#define COEF_FUNCTION_FORM_BASED_HH

#include <boost/tr1/type_traits.hpp>


#include "CoefFunction.hh"
#include "FeBasis/BaseFE.hh"
#include "FeBasis/FeFunctions.hh"
#include "Forms/Operators/BaseBOperator.hh"


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
  
  //! Set integrator for specific region
  virtual void AddIntegrator(BaseBDBInt* form, RegionIdType region);
  
  //! Return type of entry (scalar, vector, tensor)
  virtual CoefDimType GetDimType() const { return dimType_;  }

protected:
  
  //! Store bilinearform for each region
  std::map<RegionIdType, BaseBDBInt* > forms_;
  
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
  
  //! Set integrator for specific region
  virtual void AddIntegrator( BaseBDBInt* form,  
                              RegionIdType region );
  
  //! Pass directly a B-operator
  void AddBOperator( BaseBOperator* bOp,
                     RegionIdType region );

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

  //! Differential operator for each region
  std::map<RegionIdType, BaseBOperator* > bOps_;
 
  //! FeFunction containing the coefficients
  shared_ptr<FeFunction<TYPE> > feFct_;

  //! Result info object of result to be calculated
  shared_ptr<ResultInfo> res_;
  
  //! Point to FeSpace
  shared_ptr<FeSpace> feSpace_;
  
  //! Solution of element
  Vector<TYPE> elemSol_;
  
  //! Operator matrix
  Matrix<TYPE> bMat_;
  
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
  virtual ~CoefFunctionFlux();

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
  Vector<TYPE> elemSol;

  //! Additional factor
  TYPE factor_;
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
  Matrix<TYPE> kernel_;
  
  //! Kernel of element matrix (always real-valued)
  Matrix<Double> kernelR_;
  
  //! Solution of element
  Vector<TYPE> elemSol_;

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

  //! \copydoc CoefFunction::GetTensorSize
  virtual void GetTensorSize(unsigned int& numRows, unsigned int& numCols ) const;

  virtual void GetTensor(Matrix<TYPE>& tensor, const LocPointMapped& lpm);

  //! \copydoc CoefFunction::ToString
  virtual std::string ToString() const;

protected:

  //! FeFunction containing the coefficients
  shared_ptr<FeFunction<TYPE> > feFct_;

  //! Solution of element
  Vector<TYPE> elemSol_;
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

  virtual void GetScalar(TYPE& coefScal, const LocPointMapped& lpm);

  //! \copydoc CoefFunction::ToString
  virtual std::string ToString() const;

protected:

  //! FeFunction containing the coefficients
  shared_ptr<FeFunction<TYPE> > feFct_;

  //! Solution of element
  Vector<TYPE> elemSol_;
};




} // end of namespace
#endif
