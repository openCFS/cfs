// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionGridNodalDefault.cc 
 *       \brief    Implementation File
 *
 *       \date     May. 12, 2015
 *       \author   Manfred Kaltenbacher
 */
//================================================================================================

#include <def_expl_templ_inst.hh>
#include <complex>

#include "CoefFunctionGridNodalSource.hh"
#include "FeBasis/FeSpace.hh"

namespace CoupledField{


template<typename DATA_TYPE>
CoefFunctionGridNodalSource<DATA_TYPE>::CoefFunctionGridNodalSource(Domain* ptDomain,
        PtrParamNode configNode,PtrParamNode curInfo)
		:CoefFunctionGridNodal<DATA_TYPE>(ptDomain, configNode) {

  //====================================================
  // Determine information about source grid and result
  //====================================================

  this->inputId_ = "default";
  this->gridId_ = "default";
  this->curInterpType_ = CoefFunctionGrid::NO_INTERPOLATION;
  this->conservativeReady_ = false;
  this->srcIsSurface_ = false;

  this->extDataInfo_ = curInfo->Get("defaultGrid",ParamNode::APPEND);
  this->extDataInfo_->Get("interpolation")->Get("type")->SetValue("noInterpolation");

  std::string factorString = this->factorFnc_->ToString();
  this->extDataInfo_->Get("interpolation")->Get("factor",ParamNode::INSERT)->SetValue(factorString);

  this->srcGrid_ = this->domain_->GetGrid();

  //lets determine the destination region and set it to our source regions
  std::string destreg = configNode->GetParent()->GetParent()->Get("name")->As<std::string>();

  //define regularization parameter
  //beta_ = 0.5;
  isDataReadFromFile_ = false;

  //obtain entitylist from grid and add it
  nodeListSource_ = this->srcGrid_->GetEntityList(EntityList::NODE_LIST,destreg);
  //this->AddEntityList( nodeListSource_ );

  //resize the source, source-incremental and source-save vectors
  UInt numNodes = nodeListSource_->GetSize();
  sourceAmp_.Resize( numNodes );
  sourceAmp_.Init();
  sourcePhi_.Resize( numNodes );
  sourcePhi_.Init();
  sourceAmpDelta_.Resize( numNodes );
  sourceAmpDelta_.Init();
  sourcePhiDelta_.Resize( numNodes );
  sourcePhiDelta_.Init();
  sourceAmpSave_.Resize( numNodes );
  sourceAmpSave_.Init();
  sourcePhiSave_.Resize( numNodes );
  sourcePhiSave_.Init();

  this->DetermineResult(this->inputId_,this->aSeqStep_);
  this->dimDof_ = this->resultInfo_->dofNames.GetSize();
  // Determine which steps are available
  this->domain_->GetResultHandler()->GetStepValues(this->inputId_,this->aSeqStep_,this->resultInfo_,this->stepValueMap_,false);
  this->curStep_ = this->stepValueMap_.begin()->first;
  this->curTStep_ = this->stepValueMap_.begin()->second;


  //check for data of inverse scheme
  std::string inverseString;
  configNode->GetValue("dependtype",inverseString);

  if ( inverseString == "INVSOURCE") {
	  this->inverseType_ = CoefFunction::INVSOURCE;
	  std::cout << "Generate SOURCE " << std::endl;
	  std::cout << "Do Source" << "  numNodes: " << numNodes << std::endl;
  }
  else {
	  std::cout << "Generate MEASURE " << std::endl;
	  this->inverseType_ = CoefFunction::INVMEASURE;
	  configNode->GetValue("measureNodes",inverseString);
	  std::cout << "Name: " << inverseString << std::endl;

	  this->srcGrid_->GetNodesByName( measNodes_, inverseString );
	  for ( UInt i=0; i<measNodes_.GetSize(); i++)
		  std::cout << "Node: " << measNodes_[i] << std::endl;


	  //set number of equations
	 // this->numEqns_ = measNodes_.GetSize();
  }

  //GridcoefFunctions are always general!
  this->dependType_ = CoefFunction::GENERAL;
}

// ========================
//  ACCESS METHODS
// ========================
template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::GetTensor(Matrix<DATA_TYPE>& CoefMat,
                       const LocPointMapped& lpm ){
  EXCEPTION("TENSORIAL RESULTS NOT SUPPORTED FOR INTERPOLATION FROM EXTERNAL GRIDS");
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::GetVector(Vector<DATA_TYPE>& CoefMat,
                       const LocPointMapped& lpm ){
  //no tensors right now
  assert(this->dimType_ != CoefFunction::TENSOR);

  //cover the case of nc_surfElems
  const Elem* sourceElem = NULL;

  this->UpdateSolution();
   //ok we always take the volume element
  //because even if we are coping here with surfaces, e.g. for boundary conditions,
  //the volume element connectivity should give also the boundary nodes and
  //we do not cover special, nonconforming cases here...
  if(lpm.isSurface && !this->srcIsSurface_)
    sourceElem = lpm.lpmVol->ptEl;
  else
    sourceElem = lpm.ptEl;

  Vector<DATA_TYPE> elemSol;
  if(this->dimType_ == CoefFunction::SCALAR)
    CoefMat.Resize(1);
  else if (this->dimType_ == CoefFunction::VECTOR)
    CoefMat.Resize(this->dimDof_);
  
  if(!sourceElem){
    EXCEPTION("Could not determine source element.")
  }

  this->GetElemSolution( elemSol, sourceElem->elemNum);

  shared_ptr<ElemShapeMap> esm = this->srcGrid_->GetElemShapeMap( sourceElem, true );
  BaseFE * ptFe = esm->GetBaseFE();
  if(lpm.isSurface)
    this->myOperator_->ApplyOp(CoefMat,(*lpm.lpmVol),ptFe,elemSol);
  else
    this->myOperator_->ApplyOp(CoefMat,lpm,ptFe,elemSol);



}

template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::GetScalar(DATA_TYPE& CoefMat,
                       const LocPointMapped& lpm ){
  //no tensors right now
  assert(this->dimType_ != CoefFunction::TENSOR);
  //this is really simple. we just take the nodal result and
  //interpolate it to lpm

  this->UpdateSolution();

  //cover the case of nc_surfElems
  const Elem* sourceElem = NULL;

  if(lpm.isSurface && !this->srcIsSurface_)
    sourceElem = lpm.lpmVol->ptEl;
  else
    sourceElem = lpm.ptEl;
  Vector<DATA_TYPE> elemSol;
  Vector<DATA_TYPE> ptSol(1);
  ptSol.Init();

  sourceElem = lpm.ptEl;
  this->GetElemSolution( elemSol, sourceElem->elemNum);
  shared_ptr<ElemShapeMap> esm = this->srcGrid_->GetElemShapeMap( sourceElem, true );
  BaseFE * ptFe = esm->GetBaseFE();
  if(lpm.isSurface)
    this->myOperator_->ApplyOp(ptSol,(*lpm.lpmVol),ptFe,elemSol);
  else
    this->myOperator_->ApplyOp(ptSol,lpm,ptFe,elemSol);

  CoefMat = ptSol[0];
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::AddEntityList(shared_ptr<EntityList> ent){
    //this function may not be called if eqnMapping has already
    //been performed because there will be exactly one srcRegion to one target 
    //region and both are identical
    if(this->eqnMapComplete_ && !this->entities_.Contains(ent))
      EXCEPTION("Call to AddEntities after eqnMapping has been performed. This should not happen!");
    
    //first we check if its a region list.CoefFunctionGrid
    //otherwise this will be shit
    if(ent->GetDefineType() != EntityList::REGION)
      EXCEPTION("CoefFunctionGridNodalSource can only be used with regions");

    if(!this->entities_.Contains(ent)){
      this->entities_.Push_back(ent);
    }else{
      WARN("entity list " << ent->GetName() << " already contained in CoefFunction")
    }

    //in this special coeffunction,
    //the sourceRegion is equal to the destination
    this->srcRegions_.insert(ent->GetName());

    if(ent->GetType() == EntityList::SURF_ELEM_LIST){
      this->srcIsSurface_ = true;
    }


    //====================================================
    // Create Data structures for easy solution access
    //====================================================                                          
    //in this special class we start with the equation mapping right
    //after the first call to get entities as this method should not be called twice!
    this->MapEqns();
    //read in the first solution
    this->ReadSolution(this->stepValueMap_.begin()->first,this->solVec_);
    this->extDataInfo_->Get("RegionList")->Get("SourceRegion")->Get("name")->SetValue(ent->GetName());
  }

template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::GetScalarValuesAtPoints( const StdVector<Vector<Double> >  & points,
                                                                       StdVector<DATA_TYPE >  & vals){
  assert(this->dimType_ == CoefFunction::SCALAR);

  vals.Resize(points.GetSize());
  vals.Init();
  StdVector< const Elem* > elements;
  StdVector<LocPoint> locals;
  this->GetElemsForPoints(points,elements,locals);
  shared_ptr<ElemShapeMap> esm;
  LocPointMapped lpm;
  for(UInt curPoint=0;curPoint < points.GetSize();++curPoint){
    const Elem* sourceElem = elements[curPoint];
    esm = this->srcGrid_->GetElemShapeMap( sourceElem, true );
    LocPoint lp = locals[curPoint];
    lpm.Set(lp,esm,1.0);
    this->GetScalar(vals[curPoint],lpm);
  }
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::GetVectorValuesAtPoints( const StdVector<Vector<Double> >  & points,
                                                                      StdVector<Vector<DATA_TYPE> >  & vals){
  vals.Resize(points.GetSize(),Vector<DATA_TYPE>(this->dimDof_));

  vals.Init();
  StdVector< const Elem* >  elements;
  StdVector<LocPoint>  locals;
  this->GetElemsForPoints(points,elements,locals);
  shared_ptr<ElemShapeMap> esm;
  LocPointMapped lpm;
  for(UInt curPoint=0;curPoint < points.GetSize();++curPoint){
    const Elem* sourceElem = elements[curPoint];
    esm = this->srcGrid_->GetElemShapeMap( sourceElem, true );
    LocPoint lp = locals[curPoint];
    lpm.Set(lp,esm,1.0);
    this->GetVector(vals[curPoint],lpm);
  }
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::GetTensorValuesAtPoints( const StdVector<Vector<Double> >  & points,
                                                                       StdVector<Matrix<DATA_TYPE> >  & vals){

  assert(this->dimType_ == CoefFunction::TENSOR);
  EXCEPTION("GetTensorValuesAtPoints is not implemented");
  //vals.Resize(points.GetSize(),Matrix<DATA_TYPE>(this->numRows_,this->numCols_));
  //vals.Init();
  //StdVector< const Elem* > elements;
  //StdVector<LocPoint> locals;
  //this->GetElemsForPoints(points,elements,locals);
  //shared_ptr<ElemShapeMap> esm;
  //LocPointMapped lpm;
  //for(UInt curPoint=0;curPoint < points.GetSize();++curPoint){
  //  const Elem* sourceElem = elements[curPoint];
  //  esm = this->srcGrid_->GetElemShapeMap( sourceElem, true );
  //  LocPoint lp = locals[curPoint];
  //  lpm.Set(lp,esm,1.0);
  //  this->GetTensor(vals[curPoint],lpm);
  //}
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::MapConservative( shared_ptr<FeSpace> targetSpace,
                                                              Vector<DATA_TYPE>& feFncVec) {

  //if the targetSpace is also using grid ordering, we just take the nodal values and map them to the
  //target equations
  //if the target space has arbitrary order, this does not really make sense, so we issue an exception
  //the sources then need more some sort of  distribution than a clean mapping
  //this is the old problem with conservative interpolation. what happens if the target grid is finer than the source grid?
  std::set<std::string>::iterator regIter = this->srcRegions_.begin();
  for( ; regIter != this->srcRegions_.end(); ++regIter) {
    RegionIdType curId = this->srcGrid_->GetRegion().Parse(*regIter);
    FeSpace::MappingType curtype = targetSpace->GetMapType(curId);
    if(curtype!=FeSpace::GRID){
      EXCEPTION("You are trying to map a nodal grid to a polynomial space. This is not supported!")
    }
  }

  //IDEA> create vector of pairs and associate each eqnNumber of feFncVec a index of own solvec
  //a little more memory but very efficient updates
  if(!this->conservativeReady_){
    BuildNodeIdxAssoc(targetSpace);
    this->conservativeReady_ = true;
  }

  if ( this->isActive_ ) {
	  if (  this->inverseType_ == CoefFunction::INVMEASURE ) {
		 // std::cout << "\n DO RHS INVERSE_MEASURE" << std::endl;
		  //First, fetch the solutions at the measurement position
		  //and store them in measVec_
		  if ( !isDataReadFromFile_ ) {
			  this->UpdateSolution();
			  measVec_ = this->solVec_;
			  isDataReadFromFile_ = true;
		  }

		  //reset RHS to zero:
		  this->solVec_.Init();

		  NodeList* nodeList = new NodeList(this->srcGrid_);
		  nodeList->SetNodes(measNodes_);

		  //StdVector<Integer> measEqs(measNodes_.GetSize());
		  EntityIterator entIt = nodeList->GetIterator();

		  Vector<DATA_TYPE> actPDEsol(measNodes_.GetSize());
		  UInt k=0;
		  while( !entIt.IsEnd() ) {
			  //equation number not needed
//			  StdVector<Integer> eq;
//			  targetSpace->GetEqns(eq, entIt);
//			  measEqs[k]   = eq[0];

			  //value
			  Vector<DATA_TYPE> result;
			  this->feFunctions_[ACOU_PRESSURE]->GetEntitySolution( result, entIt );
			  actPDEsol[k] = result[0];

			  entIt++;
			  k++;
		  }

		  UInt idx=0;
		  for(UInt i=0;i<this->fctSolAssoc_.GetSize();++i) {
		      const std::pair<UInt,UInt> & curP = this->fctSolAssoc_[i];
		      if ( curP.second > 0 ) {
		    	  if ( isMeasuredNode_[i] ) {
		    		  Complex val = actPDEsol[idx] - measVec_[i];
		    	  	  this->solVec_[i] = -1.0*std::conj( val );
		    	  	  std::cout << "RHS-MEASURE: " <<  this->solVec_[i] << std::endl;
		    	  	  idx++;
		    	  }
		      }
		  }

//		  for ( UInt i=0; i<measNodes_.GetSize(); i++ ) {
//			  std::cout << "SolPrev: " << actPDEsol[i] << "  Meas: " << measVec_[i] << std::endl;
//			  Complex val = actPDEsol[i] - measVec_[i];
//			  this->solVec_[i] =  	std::conj( val );
//			  std::cout << "RHS-MEASURE: " <<  this->solVec_[i] << std::endl;
//		  }
	  }
	  else if (this->inverseType_ == CoefFunction::INVSOURCE ) {
		  //std::cout << "\n DO RHS INVERS_SOURCE" << std::endl;
		  //reset RHS
		  this->solVec_.Init();

		  //get actual solution of PDE at the nodes of source region
		  Vector<DATA_TYPE> actPDEsol(nodeListSource_->GetSize());
		  EntityIterator entIt = nodeListSource_->GetIterator();
		  UInt k=0;
		  while( !entIt.IsEnd() ) {
			  Vector<DATA_TYPE> result;
			  this->feFunctions_[ACOU_PRESSURE]->GetEntitySolution( result, entIt );
			  actPDEsol[k] = result[0];
			  //std::cout << "INVERSE_SOURCE solPrev: " << actPDEsol[k] << std::endl;
			  entIt++;
			  k++;
		  }
		  //ComputeSourceData(actPDEsol);
		  for ( UInt i=0; i<actPDEsol.GetSize(); i++ ) {
			  Complex val( sourceAmp_[i]*std::cos(sourcePhi_[i]), sourceAmp_[i]*std::sin(sourcePhi_[i]) );
			  this->solVec_[i] = val;
			  //std::cout << "RHS-Source: " << val << std::endl;
		  }
	  }
  }
  else {
	  //set RHS nodal source to zero
	  this->solVec_.Init();
//	  if ( this->inverseType_ == CoefFunction::INVSOURCE  )
//		  std::cout << "INVERSESOURCE:  set zero" << std::endl;
//	  else
//		  std::cout << "INVERSEMEASURE:  set zero" << std::endl;
  }


  for(UInt i=0;i<this->fctSolAssoc_.GetSize();++i){
    const std::pair<UInt,UInt> & curP = this->fctSolAssoc_[i];
    feFncVec[curP.first] += this->solVec_[curP.second];
    //std::cout << "Global Eq: " << curP.first << "Local Eq: " << curP.second  << std::endl;
  }
  //std::cout << "RHS-VEC:\n " << feFncVec << std::endl;

  //reset that the rhsFnc is active
  this->isActive_ = false;
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::BuildNodeIdxAssoc(shared_ptr<FeSpace> targetSpace){
	//loop over the entitylist and obtain for each element the equation numbers
	std::set<std::string>::iterator regIter = this->srcRegions_.begin();
	for( ; regIter != this->srcRegions_.end(); ++regIter) {
		shared_ptr<NodeList> actSDList( new NodeList(this->srcGrid_ ) );
		RegionIdType curId = this->srcGrid_->GetRegion().Parse(*regIter);
		actSDList->SetNodesOfRegion( curId );
		EntityIterator ents = actSDList->GetIterator();
		this->fctSolAssoc_.Reserve(this->fctSolAssoc_.GetSize()+actSDList->GetSize());
		while(!ents.IsEnd()){
			//obtain eqn and node number
			std::pair<UInt,UInt> curP;
			StdVector<Integer> eqns;
			targetSpace->GetEqns(eqns,ents);
			UInt curNode = ents.GetNode();
			//std::cout << "Node: " << curNode << std::endl;
			UInt solIdx = this->nodeIdxMap_[curNode];
			//safety check
			if(this->eqnNumbers_[solIdx].GetSize() != eqns.GetSize()){
				WARN("Detected incompatibility during external data mapping...")
        		ents++;
				continue;
			}

			for(UInt i=0; i<eqns.GetSize(); i++){
				bool nodeFound = false;
				if (  this->inverseType_ == CoefFunction::INVMEASURE ) {
					//check if node belongs to measured nodes
					for (UInt k=0; k<measNodes_.GetSize(); k++ ) {
						if ( curNode == measNodes_[k] )
							nodeFound = true;
					}
				}
				if(eqns[i] > 0){
					curP.first = eqns[i]-1;

					if ( this->inverseType_ == CoefFunction::INVMEASURE ) {
						if ( nodeFound )
							isMeasuredNode_.Push_back(true);
						else {
							isMeasuredNode_.Push_back(false);
							//std::cout << "Node does not belong to measured node" << std::endl;
						}
					}
					curP.second = this->eqnNumbers_[solIdx][i];
					//std::cout << "EqGlobal: " << curP.first << "  EqLocal: " << curP.second << std::endl;
					this->fctSolAssoc_.Push_back(curP);
				}
			}
			ents++;
		}
	}
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::GetElemsForPoints(const StdVector<Vector<Double> >  & points,
                                                                  StdVector< const Elem* > & elements,
                                                                  StdVector<LocPoint> & locals){
  //build up set of source regions
  std::set<std::string>::iterator regIter = this->srcRegions_.begin();
  StdVector<shared_ptr<EntityList> > lists; 
  for( ; regIter != this->srcRegions_.end(); ++regIter) {
    RegionIdType curId = this->srcGrid_->GetRegion().Parse(*regIter);
    shared_ptr<ElemList> newList(new ElemList(this->srcGrid_));
    newList->SetRegion(curId);
    lists.Push_back(newList);
  }

  //this is unfortunate but we need to search even for the defaultGridCase...
  this->srcGrid_->GetElemsAtGlobalCoords( points,
                                          locals,
                                          elements,
                                          lists);

}

template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::ComputeSourceData(Vector<DATA_TYPE>& actPDEsol) {

	//Vector<Double> amp(actPDEsol.GetSize());
	for (UInt i=0; i<actPDEsol.GetSize(); i++ ) {
		//amplitude
		Double val;
		Complex jphi(0,sourcePhi_[i]); //j*phi
		Complex valC = actPDEsol[i] * std::exp(jphi);
		val  = std::pow( (std::abs ( valC.real() ) / (alpha_ * qExp_)), 1.0/(qExp_-1) );
		if ( valC.real() < 0.0 )
			val *= -1.0;
		sourceAmp_[i] = val;

		//phase
		Double error = 1e3;
		Double fnc, fncDeriv;
		Double psi = sourcePhi_[i];
		Complex sol = actPDEsol[i];
		Double psiPrev = psi;
		//std::cout << "\n StartPhi: " << psi << "  Sol: " << sol << " Src: " << sourceAmp_[i] << std::endl;
		UInt iter = 1;
		while ( error > 1e-3 && iter < 100) {
			fnc = 2.0*beta_*psiPrev + sourceAmp_[i] * ( sol.real()*std::sin(psiPrev) + sol.imag()*std::cos(psiPrev) );
			fncDeriv = 2.0*beta_ + sourceAmp_[i] * ( -sol.real()*std::cos(psiPrev) + sol.imag()*std::sin(psiPrev) );
			//std::cout << "fnc: " <<  fnc << "  fncderiv: " << fncDeriv << std::endl;
			psi -= 0.5*fnc/fncDeriv;
			//std::cout << "New psi: " << psi << std::endl;
//			if ( psi < -PI/2.0)
//				psi = -PI/2.0;
//			else if (psi > PI/2.0 )
//				psi = PI/2.0;

			error = std::abs(psi - psiPrev);
			psiPrev = psi;
			//std::cout << "PhiAct: " << psi << "  ErrorPhi: " << error << std::endl;
			iter++;
		}
		if ( iter > 100 )
			EXCEPTION("ComputeSourceData: No convergence for computation of phase!");

		sourcePhi_[i] = psi;
		//std::cout << "Idx: " << i << "  Amp: " << sourceAmp_[i] << "  Phi: " << sourcePhi_[i] << std::endl;
	}
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::ComputeOptCondition(Double& valAmp, Double& valPhi) {

	Vector<DATA_TYPE> actPDEsol(nodeListSource_->GetSize());
	EntityIterator entIt = nodeListSource_->GetIterator();

	UInt k=0;
	while( !entIt.IsEnd() ) {
		Vector<DATA_TYPE> result;
		this->feFunctions_[ACOU_PRESSURE]->GetEntitySolution( result, entIt );
		actPDEsol[k] = result[0];
		//std::cout << "INVERSE_SOURCE sol: " << actPDEsol[k] << std::endl;
		entIt++;
		k++;
	}

	//compute optCondition
	valAmp = 0.0;
	valPhi = 0.0;
	for (UInt i=0; i<actPDEsol.GetSize(); i++ ) {
		//Complex sol = actPDEsol[i];
		Double psi  = sourcePhi_[i];
		//std::cout << "\n sol= " << sol << "  Phi=" << psi << std::endl;

		//amplitude
		Double val;
		val = alpha_* qExp_* std::pow( std::abs( sourceAmp_[i]), qExp_-1.0 );
		if ( sourceAmp_[i] < 0.0 )
			val *= -1.0;

		Complex jphi(0,sourcePhi_[i]); //j*phi
		Complex valC = actPDEsol[i] * std::exp(jphi);
		val -= valC.real();

//		val -=  ( sol.real()*std::cos(psi) - sol.imag()*std::sin(psi) );

		//store derivative
		sourceAmpDelta_[i] = val;
		valAmp += val*val;

		//phase
		val = 2*beta_*psi + sourceAmp_[i]*valC.imag(); //( sol.real()*std::sin(psi) + sol.imag()*std::cos(psi) );
		sourcePhiDelta_[i] = val;
		valPhi +=  val*val;
	}
//	optAmp = std::sqrt(valAmp);
//	optPhase = std::sqrt(valPhi);
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::ComputeDiff2Meas(Double& error) {

	NodeList* nodeList = new NodeList(this->srcGrid_);
	nodeList->SetNodes(measNodes_);

	//StdVector<Integer> measEqs(measNodes_.GetSize());
	EntityIterator entIt = nodeList->GetIterator();

	Vector<DATA_TYPE> actPDEsol(measNodes_.GetSize());
	UInt k=0;
	while( !entIt.IsEnd() ) {
		//value
		Vector<DATA_TYPE> result;
		this->feFunctions_[ACOU_PRESSURE]->GetEntitySolution( result, entIt );
		actPDEsol[k] = result[0];

		entIt++;
		k++;
	}

	error = 0.0;
	UInt idx=0;
	for(UInt i=0;i<this->fctSolAssoc_.GetSize();++i) {
		const std::pair<UInt,UInt> & curP = this->fctSolAssoc_[i];
		if ( curP.second > 0 ) {
			if ( isMeasuredNode_[i] ) {
				Complex val = actPDEsol[idx] - measVec_[i];
				std::cout << "Nr: " << idx << "  PDE: " << actPDEsol[idx] << " Meas:" << measVec_[i]
						<< "  Diff: " << val << std::endl;
				error += std::abs( val ) * std::abs( val ) ;
				idx++;
			}
		}
	}
//	std::cout << std::endl;
//	error = std::sqrt( error ); // (Double) idx );
}



template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::UpdateSource( Double& stepLength, bool lineSearch ) {
	//! update the source values (amplitude and phase)

	// save source
	if ( lineSearch ) {
		for ( UInt i=0; i<sourceAmp_.GetSize(); i++ ) {
			sourceAmp_[i] = sourceAmpSave_[i] - stepLength * sourceAmpDelta_[i];
			sourcePhi_[i] = sourcePhiSave_[i] - stepLength * sourcePhiDelta_[i];
		}
	}
	else {
		sourceAmpSave_ = sourceAmp_;
		sourcePhiSave_ = sourcePhi_;
//		for ( UInt i=0; i<sourceAmp_.GetSize(); i++ ) {
//			std::cout << "SourceAmp: " << sourceAmpSave_[i] << "  SourcePhi: " << sourcePhiSave_[i] << std::endl;
//		}
//		std::cout << std::endl;
	}
}


//! compute Tikhonov function
template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::ComputeTikh( Double& funcVal, Double& resSquared ) {

	funcVal = resSquared;
	for ( UInt i=0; i<sourceAmp_.GetSize(); i++ ) {
		funcVal += alpha_*std::pow( std::abs(sourceAmp_[i]), qExp_ ) +	beta_*std::pow(std::abs(sourcePhi_[i]), 2);
	}
}


}

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
//  template class CoefFunctionGridNodalSource<Double>;
  template class CoefFunctionGridNodalSource<Complex>;
#endif
