#ifndef CFS_RESULT_FUNCTOR_HH
#define CFS_RESULT_FUNCTOR_HH


#include "Domain/Results/BaseResults.hh"
#include "General/Environment.hh"
#include "Domain/ElemMapping/ElemShapeMap.hh"
#include "Domain/ElemMapping/EntityLists.hh"
#include "Forms/BiLinForms/BiLinearForm.hh"
#include "Forms/BiLinForms/BDBInt.hh"

namespace CoupledField {

//! Class for evaluating a result

//! This class encapsulates the calculation of results.
//! In the most general case, it provides only the abstract method to
//! fill a BaseResult object, as needed by the ResultHandler. In the derived
//! classes, also field values can be evaluated at arbitrary local points 
//! within elements.
//! The idea is that the ResultFunctor objects get instantiated within the
//! PDEs, which has the information, how objects have to be generated.
class ResultFunctor {
public:

  //! Constructor
  ResultFunctor ( shared_ptr<ResultInfo> info ) {
    derivType_ = NOT_DEFINED;
    resultInfo_ = info;
    dim_ = info->dofNames.GetSize();
  }
  
  virtual ~ResultFunctor() {
  }
  
  //! Typedef for spatial type of results
  //! PRIMARY: denotes primary field variables or their time-derivative
  //! VOL_FIELD: denotes fields in the volume space
  //! SURF_FIELD: denotes fields defined by normals
  //! INTEGRATED: denotes values, which get integrated
  typedef enum {NOT_DEFINED, PRIMARY, VOL_FIELD, SURF_FIELD, INTEGRATED} ResultDerivType;
  
  //! Typedef for integration accuracy
  //! MIDPOINT: only the element midpoints is used for integration of results
  //! FULL: use full accuracy, i.e. all integration points for result integration
  typedef enum {NO_ACCURACY, MIDPOINT, FULL} IntegAccuracy;
  
  //! Evaluate result for complete entity list
  virtual void EvalResult(shared_ptr<BaseResult> res ) = 0;
  
  //! Obtain resultinfo
  shared_ptr<ResultInfo> GetResultInfo() {return resultInfo_;}
  
  //! Add BDB-Integrator
  void AddIntegrator( BaseBDBInt* integrator, RegionIdType region) {
    forms_[region] = integrator;
  }
  
  //! perform some finial steps
  virtual void Finalize() {;}

  //! Obtain derivType
  ResultDerivType GetDerivType() {return derivType_;}

  //! Return Coefficient function
  virtual PtrCoefFct GetCoefFct() {
    EXCEPTION("Base Class ResultFunctor has no Coef-Function");
  }

protected:
    
  //! Type of result (primary, field, integrated etc.)
  ResultDerivType derivType_;
  
  //! Information about quantity to be calculated
  shared_ptr<ResultInfo> resultInfo_;
  
  //! Dimension of the result to be calculated
  UInt dim_;
  
  //! List of integrators
  std::map<RegionIdType, BaseBDBInt*> forms_;

  //! Store coefficient function
  PtrCoefFct coef_;
};

// --------------------------------------------------------------------------
//  RESULT BASED ON COEFFICIENT FUNCTION
// --------------------------------------------------------------------------
//! Functor for purely evaluating a coefficient function
template<class TYPE>
class FieldCoefFunctor : public ResultFunctor {
public:

  //! Constructor
  FieldCoefFunctor( PtrCoefFct coef,
                    shared_ptr<ResultInfo> inf );

  //! Destructor
  virtual ~FieldCoefFunctor();

  //! Return Coefficient function
  virtual PtrCoefFct GetCoefFct() {
    return coef_;
  }

  //! Evaluate field for complete result
  virtual void EvalResult( shared_ptr<BaseResult> res );

  //! Evaluate field at local point
  void GetVector(Vector<TYPE>& vec, 
                 const LocPointMapped& lpm);
};

// --------------------------------------------------------------------------
//  QUANTITIES DERIVED BY SURFACE / VOLUME INTEGRATION
// --------------------------------------------------------------------------

//! Calculate the result by integration of a coefficient function

//! This class takes a coefficient function and integrates it. This can be 
//! used to calculate e.g. the total energy from the energy density or the
//! total weight using the density.
//! 
//! In theory we could have an arbitrary integration order, but this
//! feature can only be implemented, if the definition of the integration order
//! is handled more flexible in the FeSpace classes.
template<class TYPE>
class ResultFunctorIntegrate : public ResultFunctor {
public:
  
  //! Constructor
  ResultFunctorIntegrate( PtrCoefFct coef,
                          shared_ptr<BaseFeFunction> feFct,
                          shared_ptr<ResultInfo> inf );

  //! Destructor
  virtual ~ResultFunctorIntegrate();

  //! Evaluate result for complete entity list
  virtual void EvalResult(shared_ptr<BaseResult> res );

  //! Return Coefficient function
  virtual PtrCoefFct GetCoefFct() {
    return coef_;
  }
  void SetAveraged(bool average) {
    average_ = average;
  }

private:

  //! Pointer to FeFunction
  shared_ptr<BaseFeFunction> feFct_;
  bool average_;
};


// --------------------------------------------------------------------------
//  Calculate the result by integration and sums up to total force
// --------------------------------------------------------------------------

template<class TYPE>
class ResultFunctorVWP : public ResultFunctor {
public:

  //! Constructor
  ResultFunctorVWP( PtrCoefFct coef,
                    shared_ptr<BaseFeFunction> feFct,
                    shared_ptr<ResultInfo> inf,
					Grid* ptGrid);

  //! Destructor
  virtual ~ResultFunctorVWP();

  //! Evaluate result for complete entity list
  virtual void EvalResult(shared_ptr<BaseResult> res );

  //! Return Coefficient function
  virtual PtrCoefFct GetCoefFct() {
    return coef_;
  }

private:

  //! Calculate element force
  //! \param F              (output) Array containing nodal forces
  //!                                (dim x nodes) of each element
  //! \param ptElem         (input)  Pointer to element
  //! \param dim            (input)  number of dofs = dim
  //! \param IsBoundaryNode (input)  contains 1, if corresponding node is a
  //!                                boundary node, otherwise 0
  void CalcElemElecForce(Matrix<Double>& Force, const EntityIterator nameIt,
		                 const Elem * ptElement, const StdVector<ShortInt> & IsBoundaryNode);

  //! Calculates the expression \f[ \frac{\delta \vert J \vert}{\delta r} /f]
  //! \param J (input) Jacobian matrix
  //! \param J_dr (input) derivative of Jacobian matrix in r-direction
  Double CalcDetJDr(const Matrix<Double> &J, const Matrix<Double> &dJ_dr);

  //! Pointer to FeFunction
  shared_ptr<BaseFeFunction> feFct_;

  //! Pointer to grid
  Grid* ptGrid_;
};

// --------------------------------------------------------------------------
//  FIELDS BASED ON ENERGY
// --------------------------------------------------------------------------
//! Functor for calculating the total energy within a list of elements

//! This class computes the energy density of the primary unknown using the
//! the kernel of the BDB-integrator. For many physical fields, the energy
//! can be calculated as 1/2 * sol^T * elemMat * sol. This is how this functor
//! calculates the energy.
//!  The BDB-bilinearforms have to get passed to this class for every region 
//! the result might get  calculated at.
template<class TYPE>
class EnergyResultFunctor : public ResultFunctor {
public:
  //! Constructor
  //! \param feFct Pointer to FeFunction
  //! \param inf ResultInfo of the final result
  //! \param factor Multiplicative factor for the element matrix (mostly 0.5)
  EnergyResultFunctor(shared_ptr<BaseFeFunction> feFct,
                      shared_ptr<ResultInfo> inf,
                      TYPE factor );

  //! Destructor
  virtual ~EnergyResultFunctor();
  
  //! Set integration accuracy
  void SetIntegAccuracy( IntegAccuracy acc );

  //! Evaluate result for complete entity list
  virtual void EvalResult(shared_ptr<BaseResult> res );


protected:

  //! FeFunction containing the coefficients, on which the functor computes
  shared_ptr<FeFunction<TYPE> > feFct_;

  //! Accuracy for integration
  IntegAccuracy accuracy_;

  //! Additional factor
  TYPE factor_;
  
  //! Solution of element
  Vector<TYPE> elemSol, temp;

  //! Element matrix (real-valued case)
  Matrix<Double> elemMatR;
  
  //! Element matrix (complex-valued case)
  Matrix<Complex> elemMatC;

  //! Mapped local points
  LocPointMapped lpm;

};
    
} // end of namespace

#endif
