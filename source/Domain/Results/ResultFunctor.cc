#include "ResultFunctor.hh"

namespace CoupledField {


// --------------------------------------------------------------------------
//  FIELDS BASED ON COEFFICIENT FUNCTION
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

  // check for (combination of node list and fefunction) or coil list
  // since the coil list is used with the FeSpaceConst which does not have elements
  if( ( entityListType == EntityList::NODE_LIST
      && typeid(*coef_) == typeid(FeFunction<TYPE>) ) ||
      (entityListType == EntityList::COIL_LIST) ) {
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
  case EntityList::SURF_ELEM_LIST:
    // loop over elements
    for ( it.Begin(); !it.IsEnd(); it++ ) {
      const Elem * el = it.GetElem();
      LocPoint lp = Elem::shapes[el->type].midPointCoord;
      LocPointMapped lpm;
      shared_ptr<ElemShapeMap> esm = 
        it.GetGrid()->GetElemShapeMap( el, true );
      lpm.Set( lp, esm, 0.0 );
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

        it.GetGrid()->GetNodeCoordinate(coord, node, true );
        
        globCoords.Push_back(coord);
      }
      StdVector<shared_ptr<EntityList> > lists(1);
      lists[0] = res->GetEntityList();
      it.GetGrid()->GetElemsAtGlobalCoords(globCoords, localCoords, 
                                           elems, lists);

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
          it.GetGrid()->GetElemShapeMap( el, true );
        lpm.Set( lp, esm, 0.0 );
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
                    shared_ptr<ResultInfo> inf,
                    TYPE factor) :
                    ResultFunctor( inf) {
  feFct_ = dynamic_pointer_cast<FeFunction<TYPE> >(feFct);
  factor_ = factor;
  derivType_ = INTEGRATED;
  
  // default: use full integration order 
  accuracy_ = FULL;
}
  
template<class TYPE> EnergyResultFunctor<TYPE>::
  ~EnergyResultFunctor() {
}

template<class TYPE> void EnergyResultFunctor<TYPE>::
SetIntegAccuracy( IntegAccuracy acc ){
  accuracy_ = acc;
}

template<class TYPE> void EnergyResultFunctor<TYPE>::
EvalResult(shared_ptr<BaseResult> res ) {
  EXCEPTION("General implementation not available");
}

template<> void EnergyResultFunctor<Double>::
EvalResult(shared_ptr<BaseResult> res ) {
  Result<Double>& actSol = static_cast<Result<Double>& >(*res);
  EntityIterator nameIt = actSol.GetEntityList()->GetIterator();
  Vector<Double>& vec = actSol.GetVector();
  vec.Resize( nameIt.GetSize() );

  // Loop over regions
  for( nameIt.Begin(); !nameIt.IsEnd(); nameIt++ ) {
    shared_ptr<EntityList> actSDList = 
        nameIt.GetGrid()->GetEntityList( EntityList::ELEM_LIST, nameIt.GetName() );
    EntityIterator elemIt = actSDList->GetIterator();

    Double tempEnergy = 0.0;
    // loop over elements
    for ( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {
      const Elem * el = elemIt.GetElem();
      // energy density is factor_ * elemSol^T * kernel * elemSol
      this->feFct_->GetElemSolution( elemSol, el);
      if( accuracy_ == FULL ) {
        // ==================
        //  FULL INTEGRATION
        // ==================
        forms_[el->regionId]->CalcElementMatrix(elemMatR, 
                                                elemIt, elemIt);
        temp = elemMatR * elemSol;
        tempEnergy += (temp * elemSol) * factor_;
        
      } else if ( accuracy_ == MIDPOINT ) {
        // =====================
        //  REDUCED INTEGRATION
        // =====================
        // Obtain FE element from feSpace and integration scheme
        IntegOrder order;
        // Use at least 2nd order integration, as the kernel
        // is composed by a product of 2 functions.
        order.SetIsoOrder(2);
        shared_ptr<FeSpace> feSpace = feFct_->GetFeSpace();
        shared_ptr<IntScheme> intScheme = feSpace->GetIntScheme();
        // Get shape map from grid
        shared_ptr<ElemShapeMap> esm = 
            elemIt.GetGrid()->GetElemShapeMap( el, true );

        // Get integration points
        StdVector<LocPoint> intPoints;
        StdVector<Double> weights;
        intScheme->GetIntPoints( Elem::GetShapeType(el->type), IntScheme::GAUSS, order, 
                                 intPoints, weights );
        // Loop over all integration points
        LocPointMapped lpm;
        for( UInt i = 0; i < intPoints.GetSize(); i++  ) {
          // Calculate for each integration point the LocPointMapped
          lpm.Set( intPoints[i], esm, weights[i] );
          forms_[el->regionId]->CalcKernel(elemMatR, lpm );
          temp = elemMatR * elemSol;
          tempEnergy += (temp * elemSol) * factor_ * lpm.jacDet * weights[i]; 
        } // loop integration points
        
      } else {
        EXCEPTION("No valid integration method defined");
      }
    } // loop elements
    vec[nameIt.GetPos()] = tempEnergy;
  }
}

template<> void EnergyResultFunctor<Complex>::
EvalResult(shared_ptr<BaseResult> res ) {
  Result<Complex>& actSol = static_cast<Result<Complex>& >(*res);
  EntityIterator nameIt = actSol.GetEntityList()->GetIterator();
  Vector<Complex>& vec = actSol.GetVector();
  vec.Resize( nameIt.GetSize() );

  // Loop over regions
  for( nameIt.Begin(); !nameIt.IsEnd(); nameIt++ ) {
    shared_ptr<EntityList> actSDList = 
        nameIt.GetGrid()->GetEntityList( EntityList::ELEM_LIST, nameIt.GetName() );
    EntityIterator elemIt = actSDList->GetIterator();

    Complex tempEnergy = 0.0;
    // loop over elements
    for ( elemIt.Begin(); !elemIt.IsEnd(); elemIt++ ) {
      const Elem * el = elemIt.GetElem();

      // energy density is 1/2 * elemSol^T * kernel * elemSol
      this->feFct_->GetElemSolution( elemSol, el);

      
      if( accuracy_ == FULL ) {
        // ==================
        //  FULL INTEGRATION
        // ==================

        // in the complex valued case, we can have either
        // real-valued matrices or complex ones.
        if( forms_[el->regionId]->IsComplex() )  {
          forms_[el->regionId]->CalcElementMatrix(elemMatC, 
                                                  elemIt, elemIt);
          temp = elemMatC * elemSol;
        } else {
          forms_[el->regionId]->CalcElementMatrix(elemMatR, 
                                                  elemIt, elemIt);
          temp = elemMatR * elemSol;
        }
        tempEnergy += (temp * Conj(elemSol) ) * factor_;
      }  else if( accuracy_ == MIDPOINT ) {

        // =====================
        //  REDUCED INTEGRATION
        // =====================

        // Obtain FE element from feSpace and integration scheme
        IntegOrder order;
        // Use at least 2nd order integration, as the kernel
        // is composed by a product of 2 functions.
        order.SetIsoOrder(2);
        shared_ptr<FeSpace> feSpace = feFct_->GetFeSpace();
        shared_ptr<IntScheme> intScheme = feSpace->GetIntScheme();
        // Get shape map from grid
        shared_ptr<ElemShapeMap> esm = 
            elemIt.GetGrid()->GetElemShapeMap( el, true );

        // Get integration points
        StdVector<LocPoint> intPoints;
        StdVector<Double> weights;
        intScheme->GetIntPoints( Elem::GetShapeType(el->type), IntScheme::GAUSS, order, 
                                 intPoints, weights );
        // Loop over all integration points
        LocPointMapped lpm;
        for( UInt i = 0; i < intPoints.GetSize(); i++  ) {
          //std::cerr << "i = " << i << ", point = " << intPoints[i] << ", weight = " << weights[i] << std::endl;
          // Calculate for each integration point the LocPointMapped
          lpm.Set( intPoints[i], esm, weights[i] );
          if( forms_[el->regionId]->IsComplex() )  {
            forms_[el->regionId]->CalcKernel(elemMatC, lpm );
            temp = elemMatC * elemSol;
          } else {
            forms_[el->regionId]->CalcKernel(elemMatR, lpm );
            temp = elemMatR * elemSol;
          }
          tempEnergy += (temp * Conj(elemSol) ) * factor_ * lpm.jacDet * weights[i]; 
        } // loop integration points
      } else {
        EXCEPTION("No valid integration method defined");
      }
    } // loop elements

    vec[nameIt.GetPos()] = tempEnergy;
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
  EntityIterator nameIt = actSol.GetEntityList()->GetIterator();
  shared_ptr<FeSpace> feSpace = feFct_->GetFeSpace();
  shared_ptr<IntScheme> intScheme = feSpace->GetIntScheme();
  Vector<TYPE>& vec = actSol.GetVector();
  UInt numDofs =actSol.GetResultInfo()->dofNames.GetSize();
  vec.Resize( nameIt.GetSize() * numDofs );
  vec.Init();
  
  // when integrating a coil result, the FeSpaceConst is not used but
  // the space of the magnetic unknowns, which is why we need to know
  // if this is a coil list in order to get the elements
  bool isCoil = nameIt.GetType() == EntityList::COIL_LIST;

  // switch depending on type of coefficient function:
  
  if( coef_->GetDimType() == CoefFunction::SCALAR ) {


    // ---------------
    //  SCALAR RESULT
    // ---------------
    // Loop over names (= regions / surface regions / named elements)
    for( nameIt.Begin(); !nameIt.IsEnd(); nameIt++ ) {
      shared_ptr<EntityList> actSDList;
      if( !isCoil ){
        actSDList =
          nameIt.GetGrid()->GetEntityList( EntityList::ELEM_LIST, nameIt.GetName() );
      } else {
        actSDList = nameIt.GetCoil()->GetElems();
      }
      EntityIterator elemIt = actSDList->GetIterator();

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
            elemIt.GetGrid()->GetElemShapeMap( el, true );

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
          lpm.Set( intPoints[i], esm, weights[i] );
          coef_->GetScalar(tempVal, lpm );
          elemVal += tempVal * lpm.jacDet * weights[i];
        } // loop integration points
        vec[nameIt.GetPos()] += elemVal;
      } // loop elements
    } // loop regions
    
  } else if( coef_->GetDimType() == CoefFunction::VECTOR ) {
    // ---------------
    //  VECTOR RESULT
    // ---------------
    // Loop over regions
    UInt pos = 0;
    for( nameIt.Begin(); !nameIt.IsEnd(); nameIt++ ) {
      shared_ptr<EntityList> actSDList;
      if( !isCoil ){
        actSDList =
          nameIt.GetGrid()->GetEntityList( EntityList::ELEM_LIST, nameIt.GetName() );
      } else {
        actSDList = nameIt.GetCoil()->GetElems();
      }
      EntityIterator elemIt = actSDList->GetIterator();

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
            elemIt.GetGrid()->GetElemShapeMap( el, true );

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
          lpm.Set( intPoints[i], esm, weights[i] );
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
