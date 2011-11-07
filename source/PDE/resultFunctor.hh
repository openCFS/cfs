// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef CFS_RESULT_FUNCTOR_HH
#define CFS_RESULT_FUNCTOR_HH

#include "General/environment.hh"
#include "Elements/elemShapeMap.hh"
#include "Domain/entityList.hh"
#include "Forms/biLinearForm.hh"

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
  ResultFunctor ( shared_ptr<BaseFeFunction> feFct,
                  shared_ptr<ResultInfo> info ) {
    derivType_ = NOT_DEFINED;
    resultInfo_ = info;
    ptGrid_ = feFct->GetGrid();
    dim_ = info->dofNames.GetSize();
  }
  
  virtual ~ResultFunctor() {}
  //! Typedef for spatial type of results
  //! PRIMARY: denotes primary field variables or their time-derivative
  //! VOL_FIELD: denotes fields in the volume space
  //! SURF_FIELD: dentes fields defined by normals
  //! INTEGRATED: denotes valus, which get integrated
  typedef enum {NOT_DEFINED, PRIMARY, VOL_FIELD, SURF_FIELD, INTEGRATED} ResultDerivType;
  
  //! Evaluate result for complete entity list
  virtual void EvalResult(shared_ptr<BaseResult> res ) = 0;
  
  //! Obtain resultinfo
  shared_ptr<ResultInfo> GetResultInfo() {return resultInfo_;}
  
  //! Add BDB-Integrator
  void AddIntegrator( BaseBDBInt* integrator, RegionIdType region) {
    forms_[region] = integrator;
  }
  
  //! Obtain derivType
  ResultDerivType GetDerivType() {return derivType_;}

protected:
  
  //! Type of result (primary, field, integrated etc.)
  ResultDerivType derivType_;
  
  //! Information about quantity to be calculated
  shared_ptr<ResultInfo> resultInfo_;
  
  //! Pointer to grid
  Grid * ptGrid_;

  //! Dimension of the result to be calcualted
  UInt dim_;
  
  //! List of integrators
  std::map<RegionIdType, BaseBDBInt*> forms_;
};

// --------------------------------------------------------------------------
//  FIELD RESULTS 
// --------------------------------------------------------------------------
//! Base class for evaluating spatially dependent fields
class BaseFieldFunctor : public ResultFunctor {
public:
  
  //! Constructor    
  BaseFieldFunctor( shared_ptr<BaseFeFunction> feFct,
                    shared_ptr<ResultInfo> info  ) 
  : ResultFunctor(feFct, info) {
    
    derivType_ = VOL_FIELD;
  }
  //! Destructor
  virtual ~BaseFieldFunctor() {}

protected:

};
// --------------------------------------------------------------------------
//  FIELD RESULTS (TEMPLATIZED)
// --------------------------------------------------------------------------
//! Templatized result functor for real/complex valued quantities
template<class TYPE>
class FieldFunctor : public BaseFieldFunctor {
public:

  //! Constructor
  FieldFunctor( shared_ptr<BaseFeFunction> feFct,
                shared_ptr<ResultInfo> info  ) 
  : BaseFieldFunctor( feFct, info) {  
    feFct_ = dynamic_pointer_cast<FeFunction<TYPE> >(feFct);
  }
  
  //! Destructor
  ~FieldFunctor() {}
  
  //! Evaluate field for complete result
  virtual void EvalResult( shared_ptr<BaseResult> res ) {

    Result<TYPE>& actSol = static_cast<Result<TYPE>& >(*res);
    EntityIterator it = actSol.GetEntityList()->GetIterator();
    Vector<TYPE>& vec = actSol.GetVector();
    Vector<TYPE> tempField;
    vec.Resize( it.GetSize() * this->dim_ );

    // loop over elements
    for ( it.Begin(); !it.IsEnd(); it++ ) {
      const Elem * el = it.GetElem();
      LocPoint lp = Elem::shapes[el->type].midPointCoord;
      this->Eval(tempField, el, lp );
      // loop over dofs
      for(UInt iDim = 0; iDim < dim_; iDim++ ) {
        vec[it.GetPos()*dim_ + iDim] = tempField[iDim];
      }
    }
  }

  //! Evaluate field at one local point within a given element
  virtual void Eval(Vector<TYPE>& vec, const Elem* el, LocPoint lp) = 0;

protected:
  
  //! FeFunction containing the coefficients, on which the functor compuates
  shared_ptr<FeFunction<TYPE> > feFct_;
  
};

// --------------------------------------------------------------------------
//  FIELDS BASED ON PRIMARY RESULT
// --------------------------------------------------------------------------
// ... to be implemented (interpolation of "nodal" values to arbitrary points")

// --------------------------------------------------------------------------
//  FIELDS BASED ON DERIVATIVE (GRADIENT, CURL)
// --------------------------------------------------------------------------

//! Functor for calculating fields based on the spatial derivative 

//! This class computes the spatial derivative of the primary unknown by
//! applying the B-Ooperator of the related BDB-class. The BDB-bilinearforms
//! have to get passed to this class for every region the result might get
//! calculated at.
template<class TYPE>
class DiffFieldFunctor : public FieldFunctor<TYPE> {
public:
  
  //! Constructor
  DiffFieldFunctor( shared_ptr<BaseFeFunction> feFct,
                    shared_ptr<ResultInfo> inf ) :
                     FieldFunctor<TYPE>( feFct, inf) {
  }

  //! Destructor
  virtual ~DiffFieldFunctor() {}

  //! Evaluate field at local point
  virtual void Eval(Vector<TYPE>& vec, const Elem* el, LocPoint lp) {
    vec.Resize(this->ptGrid_->GetDim());
    shared_ptr<ElemShapeMap> esm = this->ptGrid_->GetElemShapeMap( el, true );
    lpm.Set( lp, esm );
    this->feFct_->GetElemSolution( elemSol, el);
    this->forms_[el->regionId]->ApplyBMat(vec, elemSol, el, lpm );
  }
  
protected:
  
  //! Solution of element
  Vector<TYPE> elemSol;
  
  //! Mapped local points
  LocPointMapped lpm;
};

// --------------------------------------------------------------------------
//  FIELDS BASED ON FLUX (dB-Intergator)
// -------------------------------------------------------------------------
//! Functor for calculating fields based on the flux

//! This class computes the flux  of the primary unknown by
//! applying the dB-Ooperator of the related BDB-class. The BDB-bilinearforms
//! have to get passed to this class for every region the result might get
//! calculated at.
template<class TYPE>
class FluxFieldFunctor : public FieldFunctor<TYPE> {
public:
  
  //! Constructor
  FluxFieldFunctor(shared_ptr<BaseFeFunction> feFct,
                   shared_ptr<ResultInfo> inf ) :
                     FieldFunctor<TYPE>( feFct, inf) {
  }
  
  //! Destructor
  virtual ~FluxFieldFunctor() {}
  
  //! Evaluate field at local point
  virtual void Eval(Vector<TYPE>& vec, const Elem* el, LocPoint lp) {
    vec.Resize(this->ptGrid_->GetDim());
    shared_ptr<ElemShapeMap> esm = this->ptGrid_->GetElemShapeMap( el, true );
    lpm.Set( lp, esm );
    this->feFct_->GetElemSolution( elemSol, el);
    BaseBDBInt* bdb = this->forms_[el->regionId]; 
    bdb->ApplydBMat(vec, elemSol, el, lpm );
  }
  
protected:
  
  //! Solution of element
  Vector<TYPE> elemSol;

  //! Mapped local points
  LocPointMapped lpm;

};

// --------------------------------------------------------------------------
//  FIELDS BASED ON ENERGY DENSITY (KERNEL OF BDB-INTEGRATOR)
// --------------------------------------------------------------------------
//! Functor for calculating energy density using kernel of BDB-integrator

//! This class computes the energy density of the primary unknown using the
//! the kernel of the BDB-integrator. For many physical fields, the energy
//! can be calculated as 1/2 * sol^T * elemMat * sol. This is how this functor
//! calculates the energy.
//!  The BDB-bilinearforms have to get passed to this class for every region 
//! the result might get  calculated at.
template<class TYPE>
class EnergyDensFieldFunctor : public FieldFunctor<TYPE> {
public:

  //! Constructor
  EnergyDensFieldFunctor(shared_ptr<BaseFeFunction> feFct,
                         shared_ptr<ResultInfo> inf ) :
                           FieldFunctor<TYPE>( feFct, inf) {
  }
  
  //! Destructor
  ~EnergyDensFieldFunctor() {}

  //! Evaluate field at local point
  virtual void Eval(Vector<TYPE>& vec, const Elem* el, LocPoint lp) {
    vec.Resize(1);
    shared_ptr<ElemShapeMap> esm = this->ptGrid_->GetElemShapeMap( el, true );
    lpm.Set( lp, esm );
    this->forms_[el->regionId]->CalcKernel(elemMat, el, lpm);
    
    // energy density is 1/2 * elemSol^T * kernel * elemSol
    this->feFct_->GetElemSolution( elemSol, el);
    temp = elemMat * elemSol;
    vec[0] = (temp * elemSol) * 0.5;
  }
  
protected:
  //! Solution of element
  Vector<TYPE> elemSol, temp;

  //! Element matrix
  Matrix<TYPE> elemMat;

  //! Mapped local points
  LocPointMapped lpm;
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
  EnergyResultFunctor(shared_ptr<BaseFeFunction> feFct,
                      shared_ptr<ResultInfo> inf ) :
                        ResultFunctor( feFct, inf) {
    feFct_ = dynamic_pointer_cast<FeFunction<TYPE> >(feFct);
    derivType_ = INTEGRATED;
  }
  
  //! Destructor
  ~EnergyResultFunctor() {}

  //! Evaluate result for complete entity list
  virtual void EvalResult(shared_ptr<BaseResult> res ) {
    Result<TYPE>& actSol = static_cast<Result<TYPE>& >(*res);
    EntityIterator regionIt = actSol.GetEntityList()->GetIterator();
    Vector<TYPE>& vec = actSol.GetVector();
    vec.Resize( regionIt.GetSize() );

    // Loop over regions
    for( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ ) {
      ElemList actSDList(ptGrid_);
      actSDList.SetRegion( regionIt.GetRegion() );
      EntityIterator elemIt = actSDList.GetIterator();

      TYPE tempEnergy = 0.0;
      // loop over elements
      for ( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {
        const Elem * el = elemIt.GetElem();
        forms_[el->regionId]->CalcElementMatrix(elemMat, 
                                                elemIt, elemIt);

        // energy density is 1/2 * elemSol^T * kernel * elemSol
        this->feFct_->GetElemSolution( elemSol, el);
        temp = elemMat * elemSol;
        tempEnergy += (temp * elemSol) * 0.5;
      }
      vec[regionIt.GetPos()] = tempEnergy;
    }
  }
protected:
  
  //! FeFunction containing the coefficients, on which the functor compuates
  shared_ptr<FeFunction<TYPE> > feFct_;

  //! Solution of element
  Vector<TYPE> elemSol, temp;

  //! Element matrix
  Matrix<TYPE> elemMat;

  //! Mapped local points
  LocPointMapped lpm;

};
// --------------------------------------------------------------------------
//  FIELDS BASED ON NORMAL FLUX (SURFACE)
// --------------------------------------------------------------------------
// The operators on surface elements need additional information about the
// correct "side" and "normal orientation". The utilize BOperators like
// IdNormalOp, GradNormalOp or CurlNormalOp
// ... to be implemented

// --------------------------------------------------------------------------
//  QUANTITIES DERIVED BY VOLUME INTEGRATION
// --------------------------------------------------------------------------

// Here we would have a META-functor, which takes another FieldFuctor and
// calcultes the integral of it.
// This can be used e.g. for a more flexible approach to calculate the energy:
// If we use the energy density functor, we can have a different accuracy /
// integration scheme then used for the normal calculation of the element matrix


} // end of namespace


#endif
