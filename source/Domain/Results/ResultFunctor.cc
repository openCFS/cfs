#include "ResultFunctor.hh"

namespace CoupledField {


// --------------------------------------------------------------------------
//  RESULT BASED ON COEFFICIENT FUNCTION
// --------------------------------------------------------------------------
template<class TYPE> FieldCoefFunctor<TYPE>::
  FieldCoefFunctor( PtrCoefFct coef,
                    shared_ptr<ResultInfo> inf ) :
                     ResultFunctor(  inf) {
    coef_ = coef;
  }

template<class TYPE> FieldCoefFunctor<TYPE>::
~FieldCoefFunctor() {
}

  
template<class TYPE> void FieldCoefFunctor<TYPE>::
EvalResult( shared_ptr<BaseResult> res ) {

  EntityList::ListType entityListType = res->GetEntityList()->GetType();

  // check, if the node list is a 
  if( entityListType == EntityList::NODE_LIST 
      && typeid(*coef_) == typeid(FeFunction<TYPE>) ) {
    FeFunction<TYPE> & feFct= 
        dynamic_cast<FeFunction<TYPE>&> (*coef_);
    feFct.ExtractResult(res);
    return;
  }
  
  Result<TYPE>& actSol = static_cast<Result<TYPE>& >(*res);
  EntityIterator it = actSol.GetEntityList()->GetIterator();
  Vector<TYPE>& vec = actSol.GetVector();
  Vector<TYPE> tempField;
  vec.Resize( it.GetSize() * this->dim_ );

  switch( entityListType ) 
  {
  case EntityList::ELEM_LIST:
    // loop over elements
    for ( it.Begin(); !it.IsEnd(); it++ ) {
      const Elem * el = it.GetElem();
      LocPoint lp = Elem::shapes[el->type].midPointCoord;
      LocPointMapped lpm;
      shared_ptr<ElemShapeMap> esm = 
        this->ptGrid_->GetElemShapeMap( el, true );
      lpm.Set( lp, esm );
      this->GetVector(tempField, lpm );
      // loop over dofs
      for(UInt iDim = 0; iDim < dim_; iDim++ ) {
        vec[it.GetPos()*dim_ + iDim] = tempField[iDim];
      }
    }
    break;

  case EntityList::NODE_LIST:
    {
      WARN("Evaluation of nodal values is very inefficient at the moment.")

      StdVector<Vector<double> > globCoords;
      StdVector< LocPoint > localCoords;
      StdVector< const Elem* > elems;
      Vector<Double> coord;

      for ( it.Begin(); !it.IsEnd(); it++ ) {
        UInt node = it.GetNode();

        ptGrid_->GetNodeCoordinate(coord, node, true );
        
        globCoords.Push_back(coord);
      }
      
      ptGrid_->GetElemsAtGlobalCoords(globCoords, localCoords, elems);

      UInt numElems = elems.GetSize();

      if(numElems != globCoords.GetSize()) 
      {
        EXCEPTION("Element and node vectors have different size. "
                  << "Cannot perform evaluation of node results.");
      }
      
      // loop over elems
      for ( UInt i=0; i < numElems; i++ ) {
        const Elem * el = elems[i];
        LocPoint& lp = localCoords[i];
        LocPointMapped lpm;
        shared_ptr<ElemShapeMap> esm = 
          this->ptGrid_->GetElemShapeMap( el, true );
        lpm.Set( lp, esm );
        this->GetVector(tempField, lpm );
        // loop over dofs
        for(UInt iDim = 0; iDim < dim_; iDim++ ) {
          vec[i*dim_ + iDim] = tempField[iDim];
        }
      }

    }    
    break;
  default:
    EXCEPTION("Cannot extract result for entity list type '" <<
              EntityList::listType.ToString(entityListType) << "'.");    
    break;
    
  }
  
}

template<class TYPE> void FieldCoefFunctor<TYPE>::
GetVector(Vector<TYPE>& vec, 
          const LocPointMapped& lpm) {
  switch( coef_->GetDimType()) {
    case CoefFunction::VECTOR:
      coef_->GetVector( vec, lpm );
      break;
    case CoefFunction::SCALAR:
      vec.Resize(1);
      coef_->GetScalar( vec[0], lpm );
      break;
    case CoefFunction::TENSOR:
      EXCEPTION("Not implemented");
      break;
    default:
      EXCEPTION("Missing case statement");
      break;
  }
}

template class FieldCoefFunctor<Double>;
template class FieldCoefFunctor<Complex>;

// --------------------------------------------------------------------------
//  FIELDS BASED ON ENERGY
// --------------------------------------------------------------------------

template<class TYPE> EnergyResultFunctor<TYPE>::
EnergyResultFunctor(shared_ptr<BaseFeFunction> feFct,
                    shared_ptr<ResultInfo> inf ) :
                    ResultFunctor( inf) {
  feFct_ = dynamic_pointer_cast<FeFunction<TYPE> >(feFct);
  derivType_ = INTEGRATED;
}
  
template<class TYPE> EnergyResultFunctor<TYPE>::
  ~EnergyResultFunctor() {
}

template<class TYPE> void EnergyResultFunctor<TYPE>::
EvalResult(shared_ptr<BaseResult> res ) {
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

template class EnergyResultFunctor<Complex>;
template class EnergyResultFunctor<Double>;
// --------------------------------------------------------------------------
//  QUANTITIES DERIVED BY SURFACE / VOLUME INTEGRATION
// --------------------------------------------------------------------------

template<class TYPE> ResultFunctorIntegrate<TYPE>::
ResultFunctorIntegrate( PtrCoefFct coef,
                        shared_ptr<BaseFeFunction> feFct,
                        shared_ptr<ResultInfo> inf ) :
                        ResultFunctor( inf) {
  derivType_ = INTEGRATED;
  coef_ = coef;
  feFct_ = feFct;
}
    
template<class TYPE> ResultFunctorIntegrate<TYPE>::
~ResultFunctorIntegrate() {
  
}
  

template<class TYPE> void ResultFunctorIntegrate<TYPE>:: 
 EvalResult(shared_ptr<BaseResult> res ) {
  Result<TYPE>& actSol = static_cast<Result<TYPE>& >(*res);
  EntityIterator regionIt = actSol.GetEntityList()->GetIterator();
  shared_ptr<FeSpace> feSpace = feFct_->GetFeSpace();
  shared_ptr<IntScheme> intScheme = feSpace->GetIntScheme();
  Vector<TYPE>& vec = actSol.GetVector();
  UInt numDofs =actSol.GetResultInfo()->dofNames.GetSize();
  vec.Resize( regionIt.GetSize() * numDofs );

  
  // switch depending on type of coefficient function:
  
  if( coef_->GetDimType() == CoefFunction::SCALAR ) {


    // ---------------
    //  SCALAR RESULT
    // ---------------
    // Loop over regions
    for( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ ) {
      ElemList actSDList(ptGrid_);
      actSDList.SetRegion( regionIt.GetRegion() );
      EntityIterator elemIt = actSDList.GetIterator();

      TYPE tempVal = 0.0;
      // loop over elements
      for ( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {
        const Elem * el = elemIt.GetElem();


        // Obtain FE element from feSpace and integration scheme
        IntegOrder order;
        IntScheme::IntegMethod method;

        feSpace->GetFe( elemIt, method, order );
        // Get shape map from grid
        shared_ptr<ElemShapeMap> esm = 
            domain->GetGrid()->GetElemShapeMap( el, true );

        // Get integration points
        StdVector<LocPoint> intPoints;
        StdVector<Double> weights;
        intScheme->GetIntPoints( Elem::GetShapeType(el->type), method, order, 
                                 intPoints, weights );

        // Loop over all integration points
        tempVal = 0.0;
        LocPointMapped lpm;
        TYPE elemVal = 0.0;
        for( UInt i = 0; i < intPoints.GetSize(); i++  ) {

          // Calculate for each integration point the LocPointMapped
          lpm.Set( intPoints[i], esm );
          coef_->GetScalar(tempVal, lpm );
          elemVal += tempVal * lpm.jacDet * weights[i];
        } // loop integration points
        vec[regionIt.GetPos()] += elemVal;
      } // loop elements
    } // loop regions
    
  } else if( coef_->GetDimType() == CoefFunction::VECTOR ) {
    // ---------------
    //  VECTOR RESULT
    // ---------------
    // Loop over regions
    UInt pos = 0;
    for( regionIt.Begin(); !regionIt.IsEnd(); regionIt++ ) {
      ElemList actSDList(ptGrid_);
      actSDList.SetRegion( regionIt.GetRegion() );
      EntityIterator elemIt = actSDList.GetIterator();

      Vector<TYPE> tempVal;
      // loop over elements
      for ( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {
        const Elem * el = elemIt.GetElem();


        // Obtain FE element from feSpace and integration scheme
        IntegOrder order;
        IntScheme::IntegMethod method;

        feSpace->GetFe( elemIt, method, order );
        // Get shape map from grid
        shared_ptr<ElemShapeMap> esm = 
            domain->GetGrid()->GetElemShapeMap( el, true );

        // Get integration points
        StdVector<LocPoint> intPoints;
        StdVector<Double> weights;
        intScheme->GetIntPoints( Elem::GetShapeType(el->type), method, order, 
                                 intPoints, weights );

        // Loop over all integration points
        LocPointMapped lpm;
        Vector<TYPE> elemVal(numDofs);
        elemVal.Init();
        for( UInt i = 0; i < intPoints.GetSize(); i++  ) {

          // Calculate for each integration point the LocPointMapped
          lpm.Set( intPoints[i], esm );
          coef_->GetVector(tempVal, lpm );
          elemVal += tempVal * (lpm.jacDet * weights[i]);
        } // loop integration points
        for( UInt jDof = 0; jDof < numDofs; ++jDof ) {
          vec[pos+jDof] += elemVal[jDof];
        }
      } // loop elements
      pos+= numDofs;
    } // loop regions
  } else {
    EXCEPTION("Can only integrate scalar or vector results");
  }

    
}

template class ResultFunctorIntegrate<Complex>;
template class ResultFunctorIntegrate<Double>;
} // end of namespace
