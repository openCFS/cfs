#include "ResultFunctor.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Domain/Domain.hh"
#include "Domain/Results/ResultInfo.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

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
      // shared_ptr<ElemShapeMap> esm = it.GetGrid()->GetElemShapeMap(el, updatedGeo);
      // shared_ptr<ElemShapeMap> esm(el->ptrShapeMap);
      shared_ptr<ElemShapeMap> esm = (el)->GetElemShapeMap(it.GetGrid(), updatedGeo);
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
        // shared_ptr<ElemShapeMap> esm = it.GetGrid()->GetElemShapeMap( el, updatedGeo );
        // shared_ptr<ElemShapeMap> esm(el->ptrShapeMap);
        shared_ptr<ElemShapeMap> esm = (el)->GetElemShapeMap(it.GetGrid(), updatedGeo);
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
        // shared_ptr<ElemShapeMap> esm =
        //     elemIt.GetGrid()->GetElemShapeMap( el, true );
        // shared_ptr<ElemShapeMap> esm(el->ptrShapeMap);
        shared_ptr<ElemShapeMap> esm = (el)->GetElemShapeMap(elemIt.GetGrid(), true);

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
        // shared_ptr<ElemShapeMap> esm =
        //     elemIt.GetGrid()->GetElemShapeMap( el, true );
        // shared_ptr<ElemShapeMap> esm(el->ptrShapeMap);
        shared_ptr<ElemShapeMap> esm = (el)->GetElemShapeMap(elemIt.GetGrid(), true);

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
      // shared_ptr<ElemShapeMap> esm =  elemIt.GetGrid()->GetElemShapeMap( el, true );
      // shared_ptr<ElemShapeMap> esm(el->ptrShapeMap);
      shared_ptr<ElemShapeMap> esm = (el)->GetElemShapeMap(elemIt.GetGrid(), true);

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
//   Calculate the result by integration and sums up to total force
// --------------------------------------------------------------------------

template<class TYPE> ResultFunctorVWP<TYPE>::ResultFunctorVWP( PtrCoefFct coef,
                        shared_ptr<BaseFeFunction> feFct,
                        shared_ptr<ResultInfo> inf,
						Grid* ptGrid) :
                        ResultFunctor( inf) {
  derivType_ = INTEGRATED;
  coef_      = coef;
  feFct_     = feFct;
  ptGrid_    = ptGrid;
}

template<class TYPE> ResultFunctorVWP<TYPE>::~ResultFunctorVWP() {

}

template<class TYPE> void ResultFunctorVWP<TYPE>::EvalResult(shared_ptr<BaseResult> res ) {
  Result<TYPE>& actSol = static_cast<Result<TYPE>& >(*res);
  EntityIterator nameIt = actSol.GetEntityList()->GetIterator();

  // get coefFunction
  PtrCoefFct coef = GetCoefFct();

  StdVector<RegionIdType> neighborIds(1);
  Vector<TYPE>& vec = actSol.GetVector();
  vec.Resize( nameIt.GetSize() * dim_ );
  vec.Init();

  UInt pos = 0;
  // Loop over names (= regions / surface regions / named elements)
  for(nameIt.Begin(); !nameIt.IsEnd(); nameIt++)  {
    shared_ptr<EntityList> actSDList = ptGrid_->GetEntityList( EntityList::SURF_ELEM_LIST,
    		                                                            nameIt.GetName());

    //get correct volume neighbor region id
    neighborIds[0] = coef->GetVolNeighborRegionId(actSDList->GetRegion());

    // get nodes belonging to the surface elements
    StdVector<UInt> surfNodeList;
    nameIt.GetGrid()->GetNodesByRegion(surfNodeList, ptGrid_->GetRegionId(nameIt.GetName()) );

    // get volume elements next to nodes
    StdVector<const Elem*> elemList;
    ptGrid_->GetElemsNextToNodes( elemList, surfNodeList, neighborIds);

    //get memory
    StdVector<StdVector<ShortInt> > isBoundaryNode, elemNodeToCouplingNode;
    isBoundaryNode.Resize(elemList.GetSize());
    isBoundaryNode.Init();
    elemNodeToCouplingNode.Resize(elemList.GetSize());
    elemNodeToCouplingNode.Init();


    for (UInt ielem=0; ielem<elemList.GetSize(); ielem++) {
    	isBoundaryNode[ielem].Resize(elemList[ielem]->connect.GetSize());
        isBoundaryNode[ielem].Init();
        elemNodeToCouplingNode[ielem].Resize(elemList[ielem]->connect.GetSize());
        elemNodeToCouplingNode[ielem].Init();

        // Determine Boundary Nodes
        for (UInt ielemnode=0; ielemnode<isBoundaryNode[ielem].GetSize(); ielemnode++) {
          for (UInt inodes=0; inodes<surfNodeList.GetSize(); inodes++) {
        	  if (elemList[ielem]->connect[ielemnode] == surfNodeList[inodes] ) {
        		  isBoundaryNode[ielem][ielemnode] = 1;
        		  elemNodeToCouplingNode[ielem][ielemnode] = inodes;
        		  break;
            } // end if
          }
        }
    }

    EntityIterator elemIt = actSDList->GetIterator();

    Vector<Double> force(surfNodeList.GetSize()*dim_);
    force.Init();
    for (UInt ielem=0; ielem<elemList.GetSize(); ielem++) {
       	Matrix<Double> Force;
       	CalcElemElecForce( Force, elemIt, elemList[ielem], isBoundaryNode[ielem] );

        // Add the element force to the according coupling node
        for (UInt ielemnode=0; ielemnode<elemList[ielem]->connect.GetSize(); ielemnode++) {
        	for( UInt idim=0; idim<dim_; idim++) {
        		force[elemNodeToCouplingNode[ielem][ielemnode]*dim_+idim] += Force(ielemnode,idim);
        	}
        }
    } // end elements!

    Vector<Double> totalForce(dim_);
    totalForce.Init();

    for (UInt i=0; i<surfNodeList.GetSize(); i++)
    	for (UInt dim=0; dim<dim_; dim++)
    		totalForce[dim] += force[i*dim_+dim];

    //compute total force
    for(unsigned jDof = 0; jDof < dim_; jDof++ )
    	vec[pos+jDof] += totalForce[jDof];

    pos+= dim_;
  }
}


template<class TYPE> void ResultFunctorVWP<TYPE>::CalcElemElecForce(Matrix<Double>& Force,
		                                          const EntityIterator nameIt,
		                                          const Elem* ptElement,
                                                  const StdVector<ShortInt> & IsBoundaryNode) {


    Vector<Double> field;
    Matrix<Double> dJ_dr, CornerCoords,J;
    Double DetdJ_dr;

    StdVector<UInt>  connectivity = ptElement->connect;
    UInt numNodes  = connectivity.GetSize();

    // Obtain FE element from feSpace and integration scheme
    IntegOrder order;
    order.SetIsoOrder(2);
    //IntScheme::IntegMethod method;

    shared_ptr<FeSpace> feSpace = feFct_->GetFeSpace();
    // Get shape map from grid
    // shared_ptr<ElemShapeMap> esm =  ptGrid_->GetElemShapeMap( ptElement, true );
    // shared_ptr<ElemShapeMap> esm(ptElement->ptrShapeMap);
    shared_ptr<ElemShapeMap> esm = (ptElement)->GetElemShapeMap(ptGrid_, true);

    // Get integration points
    StdVector<LocPoint> intPoints;
    StdVector<Double> weights;
    shared_ptr<IntScheme> intScheme = feSpace->GetIntScheme();
    intScheme->GetIntPoints( Elem::GetShapeType(ptElement->type), IntScheme::GAUSS, order,
                                     intPoints, weights );
    //intScheme->GetIntPoints( Elem::GetShapeType(ptElement->type), method, order, intPoints, weights );

    Force.Resize(IsBoundaryNode.GetSize(), dim_);
    Force.Init();

    // Loop over all integration points
    LocPointMapped lpm;
    lpm.SetCheckJacobi(false);
    for(UInt i = 0; i < intPoints.GetSize(); i++) {
    	// Calculate for each integration point the LocPointMapped
    	lpm.Set(intPoints[i], esm, weights[i]);
    	//lpm.SetCheckJacobi(false);
    	// get field vector scaled by square of material parameter
    	coef_->GetVector(field, lpm );

    	Matrix<Double> J = lpm.jac;
    	Matrix<Double> Jinv = lpm.jacInv;

    	Matrix<Double> JinvT; Jinv.Transpose(JinvT);
    	Double Jdet = lpm.jacDet;
    	Matrix<Double> SpecCornerCoords( dim_, numNodes );
    	for (UInt nNode=0; nNode<IsBoundaryNode.GetSize(); nNode++) {
    		// loop over all dimension
    		for (UInt idim=0; idim<dim_; idim++) {
    			// form SpecCornerCoords-Array
    			SpecCornerCoords.Init();
    			if (IsBoundaryNode[nNode] == 1)
    					SpecCornerCoords[idim][nNode] = 1;
    			else
    				break;

    			// calculate dJ_dr
    			lpm.Set(intPoints[i], esm, weights[i], SpecCornerCoords);
    			Matrix<Double> dJ_dr = lpm.jac;

    			// calculate dJ_dr
    			DetdJ_dr = CalcDetJDr(J, dJ_dr);


    			Matrix<Double> J_r_Trans;
    			dJ_dr.Transpose(J_r_Trans);
    			dJ_dr = J_r_Trans;

    			Double field2 = field.Inner(field);

    			Vector<Double> dJdrTimesE(dim_); dJdrTimesE.Init();
    			dJdrTimesE = dJ_dr * field;

    			Vector<Double> JInvTimesdJdr(dim_); JInvTimesdJdr.Init();
    			JInvTimesdJdr = JinvT * dJdrTimesE;

    			Force(nNode,idim) -=  ( field.Inner(JInvTimesdJdr)*Jdet
    					              - 0.5 *  field2 * DetdJ_dr ) * weights[i];

    		} // loop over dimension
    	} // loop over boundary nodes
    } // loop over integration points

}


template<class TYPE> Double ResultFunctorVWP<TYPE>::CalcDetJDr(Matrix<Double> &J,
		                                                       Matrix<Double> &dJ_dr ) {
  Double det;

  if (J.GetNumRows() == 2) {
	  det = dJ_dr[0][0]*J[1][1]+dJ_dr[1][1]*J[0][0]-dJ_dr[0][1]*J[1][0]-dJ_dr[1][0]*J[0][1];
  }
  else {
	  det =  dJ_dr[0][0] * (J[1][1]*J[2][2] - J[1][2]*J[2][1])
          -  dJ_dr[0][1] * (J[1][0]*J[2][2] - J[1][2]*J[2][0])
          +  dJ_dr[0][2] * (J[1][0]*J[2][1] - J[1][1]*J[2][0])
          -  dJ_dr[1][0] * (J[0][1]*J[2][2] - J[0][2]*J[2][1])
          +  dJ_dr[1][1] * (J[0][0]*J[2][2] - J[0][2]*J[2][0])
          -  dJ_dr[1][2] * (J[0][0]*J[2][1] - J[0][1]*J[2][0])
          +  dJ_dr[2][0] * (J[0][1]*J[1][2] - J[0][2]*J[1][1])
          -  dJ_dr[2][1] * (J[0][0]*J[1][2] - J[0][2]*J[1][0])
          +  dJ_dr[2][2] * (J[0][0]*J[1][1] - J[0][1]*J[1][0]);
  }
  return det;
}

template class ResultFunctorVWP<Complex>;
template class ResultFunctorVWP<Double>;
} // end of namespace
