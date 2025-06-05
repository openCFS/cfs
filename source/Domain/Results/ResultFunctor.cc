#include "ResultFunctor.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Domain/Domain.hh"
#include "Domain/Results/ResultInfo.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "FeBasis/HCurl/HCurlElems.hh"

namespace CoupledField {

DEFINE_LOG(resfunc, "resultFunctor")


// --------------------------------------------------------------------------
//  FIELDS BASED ON COEFFICIENT FUNCTION
// --------------------------------------------------------------------------
template<class TYPE> FieldCoefFunctor<TYPE>::
  FieldCoefFunctor( PtrCoefFct coef,
                    shared_ptr<ResultInfo> inf ) :
                     ResultFunctor(  inf) {
    coef_ = coef;
  }

template<class TYPE> FieldCoefFunctor<TYPE>::~FieldCoefFunctor() {
}

  
template<class TYPE> void FieldCoefFunctor<TYPE>::EvalResult( shared_ptr<BaseResult> res )
{
  EntityList::ListType entityListType = res->GetEntityList()->GetType();
  bool updatedGeo = true; //TODO: let the PDE set this flag
  // optimization results are generated in DesignSpace(). This includes complicated ones like opt_result_*
  if(res->GetResultInfo()->fromOptimization)
  {
    if(domain->HasDesign())
      domain->GetDesign()->ExtractResults<TYPE>(res);
    else {
      Vector<TYPE>* data = dynamic_cast<Vector<TYPE>* >(res->GetSingleVector());
      data->Resize(res->GetEntityList()->GetSize(), 0.0);
    }
    return;
  }

  LOG_DBG(resfunc) << "ER " << res->ToString();

  // TODO element based interpolation from elemResults to nodeResults very slow
  // check for (combination of node list and fefunction) or coil list
  // since the coil list is used with the FeSpaceConst which does not have elements
  if( ( entityListType == EntityList::NODE_LIST && typeid(*coef_) == typeid(FeFunction<TYPE>) )
      || (entityListType == EntityList::COIL_LIST) )
  {
    FeFunction<TYPE> & feFct=  dynamic_cast<FeFunction<TYPE>&> (*coef_);
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
    for(it.Begin(); !it.IsEnd(); it++)
    {
      const Elem * el = it.GetElem();
      LocPoint lp = Elem::shapes[el->type].midPointCoord;
      LocPointMapped lpm;
      shared_ptr<ElemShapeMap> esm = it.GetGrid()->GetElemShapeMap(el, updatedGeo);
      lpm.Set(lp, esm, 0.0);
      this->GetVector(tempField, lpm );
      // loop over dofs
      for(UInt iDim = 0; iDim < dim_; iDim++ )
        vec[it.GetPos()*dim_ + iDim] = tempField[iDim];
    }
    break;

  case EntityList::NODE_LIST:
    {
      // WARN("Evaluation of nodal values is very inefficient at the moment.")

      StdVector<Vector<double> > globCoords;
      StdVector< LocPoint > localCoords;
      StdVector< const Elem* > elems;
      Vector<Double> coord;

      for ( it.Begin(); !it.IsEnd(); it++ ) {
        UInt node = it.GetNode();

        it.GetGrid()->GetNodeCoordinate(coord, node, updatedGeo );
        
        globCoords.Push_back(coord);
      }
      StdVector<shared_ptr<EntityList> > lists(1);
      lists[0] = res->GetEntityList();
      it.GetGrid()->GetElemsAtGlobalCoords(globCoords, localCoords, 
                                           elems, lists, 1e-3, 1e-2, true, updatedGeo);

      UInt numElems = elems.GetSize();


      LOG_DBG2(resfunc) << "ER NL numElemes=" << numElems << " globCoords " << globCoords.ToString();


      if(numElems != globCoords.GetSize()) 
        EXCEPTION("Element and node vectors have different size. Cannot perform evaluation of node results.");
      
      // loop over elems
      for ( UInt i=0; i < numElems; i++ ) {
        const Elem* el = elems[i];
        LocPoint& lp = localCoords[i];
        LocPointMapped lpm;
        shared_ptr<ElemShapeMap> esm = it.GetGrid()->GetElemShapeMap( el, updatedGeo );
        lpm.Set( lp, esm, 0.0 );
        this->GetVector(tempField, lpm );

        LOG_DBG3(resfunc) << "ER NT e=" << el->elemNum << " temField=" << tempField.ToString();

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

template<class TYPE> void FieldCoefFunctor<TYPE>::GetVector(Vector<TYPE>& vec, const LocPointMapped& lpm)
{
  LOG_DBG(resfunc) << "FCF::GV "<< coef_->ToString() << " dim=" << coef_->GetDimType();

  switch(coef_->GetDimType())
  {
    case CoefFunction::VECTOR:
      coef_->GetVector( vec, lpm );
      break;
    case CoefFunction::SCALAR:
      vec.Resize(1);
      coef_->GetScalar( vec[0], lpm );
      break;
    case CoefFunction::TENSOR:
      {
        Matrix<TYPE> tmp;
        coef_->GetTensor(tmp, lpm);
        if(resultInfo_->dofNames.GetSize() == tmp.GetNumRows() * tmp.GetNumCols())
          tmp.ConvertToVec_AppendRows(vec);
        else
          tmp.ConvertToVec_UpperTriangular(vec);
        break;
      }
    default:
      EXCEPTION("Missing case statement");
      break;
  }
  LOG_DBG2(resfunc) << "FCF::GV vec=" << vec.ToString();
}

template class FieldCoefFunctor<Double>;
template class FieldCoefFunctor<Complex>;

// --------------------------------------------------------------------------
//  FIELDS BASED ON ENERGY
// --------------------------------------------------------------------------

template<class TYPE> EnergyResultFunctor<TYPE>::EnergyResultFunctor(shared_ptr<BaseFeFunction> feFct, shared_ptr<ResultInfo> inf, TYPE factor)
                      : ResultFunctor( inf)
{
  feFct_ = dynamic_pointer_cast<FeFunction<TYPE> >(feFct);
  factor_ = factor;
  derivType_ = INTEGRATED;
  
  // default: use full integration order 
  accuracy_ = FULL;
}
  
template<class TYPE> EnergyResultFunctor<TYPE>:: ~EnergyResultFunctor() {
}

template<class TYPE> void EnergyResultFunctor<TYPE>::SetIntegAccuracy( IntegAccuracy acc ){
  accuracy_ = acc;
}

template<class TYPE> void EnergyResultFunctor<TYPE>::EvalResult(shared_ptr<BaseResult> res ) {
  EXCEPTION("General implementation not available");
}

template<> void EnergyResultFunctor<Double>::EvalResult(shared_ptr<BaseResult> res ) {
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
        forms_[el->regionId]->CalcElementMatrix(elemMatR, elemIt, elemIt);
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
        intScheme->GetIntPoints(Elem::GetShapeType(el->type), IntScheme::GAUSS, order, intPoints, weights);
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

template<> void EnergyResultFunctor<Complex>::EvalResult(shared_ptr<BaseResult> res ) {
  Result<Complex>& actSol = static_cast<Result<Complex>& >(*res);
  EntityIterator nameIt = actSol.GetEntityList()->GetIterator();
  Vector<Complex>& vec = actSol.GetVector();
  vec.Resize( nameIt.GetSize() );

  LOG_DBG2(resfunc) << "ERF<C>:ER " << res->GetResultInfo()->resultName << " acc=" << accuracy_ << " factor=" << factor_;

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
          forms_[el->regionId]->CalcElementMatrix(elemMatC, elemIt, elemIt);
          temp = elemMatC * elemSol;
        } else {
          forms_[el->regionId]->CalcElementMatrix(elemMatR, elemIt, elemIt);
          temp = elemMatR * elemSol;
        }
        tempEnergy += (temp.Inner(elemSol) ) * factor_;
        LOG_DBG3(resfunc) << "ERF<C>:ER e=" << el->elemNum << " sol=" << elemSol.ToString() << " spf=" << (temp.Inner(elemSol) ) * factor_ << " -> " << tempEnergy;

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
         tempEnergy += temp.Inner(elemSol) * factor_ * lpm.jacDet * weights[i];
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

template<class TYPE> ResultFunctorIntegrate<TYPE>::ResultFunctorIntegrate( PtrCoefFct coef,
                        shared_ptr<BaseFeFunction> feFct,
                        shared_ptr<ResultInfo> inf) :
                        ResultFunctor( inf) {
  derivType_ = INTEGRATED;
  coef_ = coef;
  feFct_ = feFct;
  average_ = false;
}

    
template<class TYPE> ResultFunctorIntegrate<TYPE>::~ResultFunctorIntegrate() {
  
}
  
template<class TYPE> void ResultFunctorIntegrate<TYPE>::EvalResult(shared_ptr<BaseResult> res ) {
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

  // hande the different cases commonly instead of copy & paste hell!
  TYPE tempValScal = 0.0;
  TYPE elemValScal = 0.0;
  TYPE elemVolume = 0.0;
  Vector<TYPE> tempValVec;
  Vector<TYPE> elemValVec(numDofs);
  Matrix<TYPE> tempValMat;
  Matrix<TYPE> elemValMat;

  // Loop over names (= regions / surface regions / named elements)
  UInt pos = 0;
  for(nameIt.Begin(); !nameIt.IsEnd(); nameIt++)
  {
    shared_ptr<EntityList> actSDList = isCoil ? nameIt.GetCoil()->GetElems()
                                              : nameIt.GetGrid()->GetEntityList( EntityList::ELEM_LIST, nameIt.GetName());

    EntityIterator elemIt = actSDList->GetIterator();
    TYPE regionVolume = 0.0;
    // loop over elements
    for (elemIt.Begin(); !elemIt.IsEnd(); elemIt++)
    {
      switch(coef_->GetDimType())
      {
      case CoefFunction::SCALAR:
        tempValScal = 0.0;
        elemValScal = 0.0;
        break;
      case CoefFunction::VECTOR:
        tempValVec.Init();
        elemValVec.Init();
        break;
      case CoefFunction::TENSOR:
        tempValMat.Init();
        elemValMat.Init();
        break;
      default:
        assert(false);
      }

      const Elem* el = elemIt.GetElem();

      // Obtain FE element from feSpace and integration scheme
      IntegOrder order;
      IntScheme::IntegMethod method;

      feSpace->GetFe(elemIt, method, order);
      // Get shape map from grid
      shared_ptr<ElemShapeMap> esm =  elemIt.GetGrid()->GetElemShapeMap( el, true );

      // Get integration points
      StdVector<LocPoint> intPoints;
      StdVector<Double> weights;
      intScheme->GetIntPoints( Elem::GetShapeType(el->type), method, order, intPoints, weights );

      // Loop over all integration points
      LocPointMapped lpm;
      elemVolume =0;
      for(UInt i = 0; i < intPoints.GetSize(); i++)
      {
        // Calculate for each integration point the LocPointMapped
        lpm.Set(intPoints[i], esm, weights[i]);

        switch(coef_->GetDimType())
        {
        case CoefFunction::SCALAR:
          coef_->GetScalar(tempValScal, lpm );
          elemValScal += tempValScal * lpm.jacDet * weights[i];
          break;
        case CoefFunction::VECTOR:
          coef_->GetVector(tempValVec, lpm );
          elemValVec += tempValVec * (lpm.jacDet * weights[i]);
          break;
        case CoefFunction::TENSOR:
          coef_->GetTensor(tempValMat, lpm);
          if(i == 0)
            elemValMat.Resize(tempValMat); //
          elemValMat += tempValMat * (lpm.jacDet * weights[i]);
          break;
        default:
          assert(false);
        }
        if (average_) { // we need to compute the volume
          elemVolume += lpm.jacDet * weights[i];
        }
      } // loop integration points

      switch(coef_->GetDimType())
      {
      case CoefFunction::SCALAR:
        vec[nameIt.GetPos()] += elemValScal;
        break;
      case CoefFunction::VECTOR:
        for(unsigned jDof = 0; jDof < numDofs; ++jDof )
          vec[pos+jDof] += elemValVec[jDof];
        break;
      case CoefFunction::TENSOR:
      {
        unsigned int rows = tempValMat.GetNumRows();
        unsigned int cols = tempValMat.GetNumCols();

        LOG_DBG(resfunc) << "RFI::ER TENSOR res=" << res->ToString() << " nD=" << numDofs << " r=" << rows << " c=" << cols;

        if(numDofs == rows * cols)
          elemValMat.ConvertToVec_AppendRows(tempValVec);
        else
          elemValMat.ConvertToVec_UpperTriangular(tempValVec);

        assert(tempValVec.GetSize() == numDofs);

        // feSpace->GetFe(elemIt, method, order);
        // Get shape map from grid
        // shared_ptr<ElemShapeMap> esm = elemIt.GetGrid()->GetElemShapeMap(el, true);

        for(unsigned jDof = 0; jDof < numDofs; ++jDof)
          vec[pos+jDof] += tempValVec[jDof];

        break;
      }

      default:
        assert(false);
      }
      regionVolume += elemVolume;
    } // loop elements
    if (average_) { // divide by total volume
      switch(coef_->GetDimType())
      {
      case CoefFunction::SCALAR:
        vec[nameIt.GetPos()] /= regionVolume;
        break;
      default: // = CoefFunction::VECTOR | CoefFunction::TENSOR
        for(unsigned jDof = 0; jDof < numDofs; ++jDof )
          vec[pos+jDof] /= regionVolume;
        break;
      }
    }
    pos+= numDofs;
  } // loop regions
}

template class ResultFunctorIntegrate<Complex>;
template class ResultFunctorIntegrate<Double>;

// --------------------------------------------------------------------------
//   Calculate the result by integration and sums up to total force: VWP
// --------------------------------------------------------------------------

template<class FE, class DATA_TYPE>
ResultFunctorVWP<FE, DATA_TYPE>::ResultFunctorVWP(shared_ptr< CoefFunctionSurfVWP<FE, DATA_TYPE> > coef,
                                             shared_ptr<ResultInfo> inf)
: ResultFunctor(inf)
{
  derivType_ = INTEGRATED;
  coef_      = coef;
  surfCoef_  = coef;
}

template<class FE, class DATA_TYPE> 
ResultFunctorVWP<FE, DATA_TYPE>::~ResultFunctorVWP() {

}

template<class FE, class DATA_TYPE>
void ResultFunctorVWP<FE,DATA_TYPE>::EvalResult(shared_ptr<BaseResult> res ) {
  Result<DATA_TYPE>& actSol = static_cast< Result<DATA_TYPE>& >(*res);
  Vector<DATA_TYPE>& vec = actSol.GetVector();
  Vector<DATA_TYPE> force;
  EntityIterator nameIt = actSol.GetEntityList()->GetIterator();

  vec.Resize( nameIt.GetSize() * dim_ );
  vec.Init();

  // Loop over names (= regions / surface regions / named elements)
  for (nameIt.Begin(); !nameIt.IsEnd(); nameIt++)  {
    surfCoef_->GetTotalForce(nameIt.GetName(), force);
    for (UInt dof = 0; dof < dim_; ++dof) {
      vec[nameIt.GetPos()*dim_+dof] = force[dof];
    }
  }
}

template class ResultFunctorVWP<FeH1,Double>;
template class ResultFunctorVWP<FeHCurl,Double>;
template class ResultFunctorVWP<FeH1,Complex>;
template class ResultFunctorVWP<FeHCurl,Complex>;

} // end of namespace
