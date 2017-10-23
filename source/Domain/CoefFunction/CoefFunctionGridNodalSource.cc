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
#include <iostream>
#include <fstream>

#include "CoefFunctionGridNodalSource.hh"
#include "Domain/CoefFunction/CoefXpr.hh"
#include "FeBasis/FeSpace.hh"
#include "Materials/AcousticMaterial.hh"

namespace CoupledField{


template<typename DATA_TYPE>
CoefFunctionGridNodalSource<DATA_TYPE>::CoefFunctionGridNodalSource(Domain* ptDomain,
        PtrParamNode configNode,PtrParamNode curInfo, shared_ptr<RegionList> regions)
		:CoefFunctionGridNodal<DATA_TYPE>(ptDomain, configNode, regions) {

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
  this->srcGrid_ = this->domain_->GetGrid();

  //lets determine the destination region and set it to our source regions
  std::string destreg = configNode->GetParent()->GetParent()->Get("name")->As<std::string>();

  //define regularization parameter
  //beta_ = 0.5;
  isDataReadFromFile_ = false;

  //obtain all nodes
  nodeListSource_ = this->srcGrid_->GetEntityList(EntityList::NODE_LIST,destreg);

  //resize the source, source-incremental and source-save vectors
  UInt numNodes = nodeListSource_->GetSize();
  std::cout<< "NumNodes: " << numNodes << std::endl;
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

  typeCoeff_ = inverseString;

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
  }

  //Gid CoefFunctions are always general!
  this->dependType_ = CoefFunction::GENERAL;


  // add Region list
  AddEntityList(regions);

  //define the operator
  this->CreateOperator(this->srcGrid_->GetDim(),this->dimDof_);
}


template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::SetInverseParam( Double& alpha, Double& beta,
		Double& qExp, Double& freq, std::string fileName) {

	alpha_ = alpha;
	beta_  = beta;
	qExp_  = qExp;
	freq_ = freq;

	if ( this->inverseType_ == CoefFunction::INVMEASURE ) {
		std::cout << "Read measured data from FILE" << std::endl;
		std::string line;
		std::ifstream myfile ( fileName.c_str() );

		Integer nodeNr;
		Double realPart, imagPart;

		//get number of measurement positions
		UInt numMeas=0;
		if ( myfile.is_open() ) {
			while ( getline (myfile,line) ) {
				numMeas++;
			}
			myfile.close();
		}
		else
			EXCEPTION("Could not open file including measurement data");


		//resize
		measNodes_.Resize(numMeas);
		readMeasVec_.Resize(numMeas);

		//reopen and read measurement data
		std::ifstream myfile1 ( fileName.c_str() );
		UInt idx=0;
		if ( myfile1.is_open() ) {
			while ( getline (myfile1,line) ) {
				//std::cout << line << '\n';
				std::sscanf(line.c_str(), "%d%lf%lf", &nodeNr, &realPart, &imagPart);
				measNodes_[idx] = nodeNr;
				readMeasVec_[idx] = Complex(realPart,imagPart);
				idx++;
			}
			myfile1.close();
		}
		else
			EXCEPTION("Could not open file including measurement data");

		for ( UInt i=0; i<measNodes_.GetSize(); i++)
			std::cout << "Node: " << measNodes_[i] << "  Value: " <<  readMeasVec_[i] << std::endl;

	}

	scalingHesse_ = 1.0;
	//compute the scaling for the source terms
	if ( this->inverseType_ == CoefFunction::INVSOURCE &&
			this->approxSourceType_ != CoefFunction::DELTA) {

		Double meanSOS = 340.0; //fix speed of sound
		Double waveNumber = 2*M_PI*freq_ / std::pow( meanSOS, 2.0 );;
		Double factorPropStiff = std::pow( meanElemVol_, 1.0/(Double)this->srcGrid_->GetDim() );
		Double factorPropMass  = waveNumber * waveNumber * meanElemVol_;

		if ( factorPropMass > factorPropStiff ) {
			scalingHesse_ = std::pow( waveNumber, 4.0 );
		}
		else {
			scalingHesse_ = std::pow( factorPropStiff, -4.0 );
		}
	}

	std::cout << "SCALING: " << scalingHesse_  << std::endl;
}


template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::AddEntityList(shared_ptr<EntityList> ent){
    //this function may not be called if eqnMapping has already
    //been performed because there will be exactly one srcRegion to one target 
    //region and both are identical

	if(!this->entities_.Contains(ent)){
      this->entities_.Push_back(ent);
    }else{
      WARN("entity list " << ent->GetName() << " already contained in CoefFunction")
    }

    this->srcRegions_.insert(ent->GetName());
    std::cout << "SRC Region: " << ent->GetName() << std::endl;

    //compute mean element volume
    if ( this->inverseType_ == CoefFunction::INVSOURCE ) {
    	//get all elements belonging to region
    	meanElemVol_ = 0.0;
    	UInt numEl=0;
    	std::set<std::string>::iterator regIter = this->srcRegions_.begin();
    	for( ; regIter != this->srcRegions_.end(); ++regIter) {
    	    RegionIdType curId = this->srcGrid_->GetRegion().Parse(*regIter);
    	    StdVector<Elem*> elems;
    	    this->srcGrid_->GetElems(elems, curId);
    	    for ( UInt i=0; i<elems.GetSize(); i++) {
    	    	numEl++;
    	    	const Elem* ptElem = elems[i];
    	    	shared_ptr<ElemShapeMap> esm = this->srcGrid_->GetElemShapeMap( ptElem, true );
    	    	meanElemVol_ += esm->CalcVolume();
    	    }
    	}
    	meanElemVol_ /= (Double)numEl;
    	std::cout << "Mean element volume: " <<  meanElemVol_  << std::endl;
    }

    //====================================================
    // Create Data structures for easy solution access
    //====================================================                                          
    //in this special class we start with the equation mapping right
    //after the first call to get entities as this method should not be called twice!
    this->MapEqns();
    this->InitSolVec();

    //read in the first solution: no longer necessary, since we obtain the measure
    //data from the input file "measDataFilename" specified in xml-file
    //if ( this->inverseType_ == CoefFunction::INVMEASURE ) {
    // 	this->ReadSolution(this->stepValueMap_.begin()->first,this->solVec_);
    //}
    //this->extDataInfo_->Get("RegionList")->Get("SourceRegion")->Get("name")->SetValue(ent->GetName());
  }


template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::MapConservative( shared_ptr<FeSpace> targetSpace,
                                                              Vector<DATA_TYPE>& feFncVec) {

  std::set<std::string>::iterator regIter = this->srcRegions_.begin();
  for( ; regIter != this->srcRegions_.end(); ++regIter) {
    RegionIdType curId = this->srcGrid_->GetRegion().Parse(*regIter);
    FeSpace::MappingType curtype = targetSpace->GetMapType(curId);
    if(curtype!=FeSpace::GRID){
      EXCEPTION("You are trying to map a nodal grid to a polynomial space. This is not supported!")
    }
  }

  //IDEA: create vector of pairs and associate each eqnNumber of feFncVec a index of own solvec
  //a little more memory but very efficient updates
  if(!this->conservativeReady_){
    BuildNodeIdxAssoc(targetSpace);
    this->conservativeReady_ = true;
  }

  if ( this->isActive_ ) {
	  if (  this->inverseType_ == CoefFunction::INVMEASURE ) {
		  //First, fetch the solutions at the measurement position
		  //and store them in measVec_
		  if ( !isDataReadFromFile_ ) {

			  this->UpdateSolution();
			  measVec_ = this->solVec_;
			  //std::cout << "Meas vec set!!" << std::endl;
			  isDataReadFromFile_ = true;

			  //==================NEW=============================================
			  //set measured data read from file to measVec_!
			  UInt idx = 0;
			  for(UInt i=0;i<this->fctSolAssoc_.GetSize();++i) {
				  const std::pair<UInt,UInt> & curP = this->fctSolAssoc_[i];
				  if ( curP.second > 0 ) {
					  if ( isMeasuredNode_[i] ) {
						  //std::cout << "Pos: " << i << "  Value: " << measVec_[i] << std::endl;
						  measVec_[i] = readMeasVec_[idx];
						  idx++;
					  }
				  }
			  }
		  }
		  //reset RHS to zero:
		  this->solVec_.Init();

		  NodeList* nodeList = new NodeList(this->srcGrid_);
		  nodeList->SetNodes(measNodes_);

		  //StdVector<Integer> measEqs(measNodes_.GetSize());
		  EntityIterator entIt = nodeList->GetIterator();

		  Vector<DATA_TYPE> actPDEsol(measNodes_.GetSize());
		  actPDEsol.Init();

		  UInt k=0;
		  while( !entIt.IsEnd() ) {
			  //value
			  Vector<DATA_TYPE> result;
			  result.Init();
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
		    	  	  this->solVec_[i] = 1.0*std::conj( val );
		    	  	  //std::cout << "RHS-ADJ: " <<  this->solVec_[i] << std::endl;
		    	  	  idx++;
		    	  }
		      }
		  }
	  }
	  else if (this->inverseType_ == CoefFunction::INVSOURCE ) {
		  //reset RHS
		  this->solVec_.Init();

		  //get actual solution of PDE at the nodes of source region
		  Vector<DATA_TYPE> actPDEsol(nodeListSource_->GetSize());
		  actPDEsol.Init();
		  EntityIterator entIt = nodeListSource_->GetIterator();
		  UInt k=0;
		  while( !entIt.IsEnd() ) {
			  Vector<DATA_TYPE> result;
			  result.Init();
			  this->feFunctions_[ACOU_PRESSURE]->GetEntitySolution( result, entIt );
			  actPDEsol[k] = result[0];
			  entIt++;
			  k++;
		  }
		  //ComputeSourceData(actPDEsol);
		  for ( UInt i=0; i<actPDEsol.GetSize(); i++ ) {
			  Complex val( sourceAmp_[i]*std::cos(sourcePhi_[i]), sourceAmp_[i]*std::sin(sourcePhi_[i]) );
			  this->solVec_[i] = val;
		  }
	  }
  }
  else {
	  //set solution vector to zero
	  this->solVec_.Init();
  }


  for(UInt i=0;i<this->fctSolAssoc_.GetSize();++i){
    const std::pair<UInt,UInt> & curP = this->fctSolAssoc_[i];
    feFncVec[curP.first] += this->solVec_[curP.second];
  }

  //reset that the rhsFnc is active
  if ( this->approxSourceType_ == CoefFunction::DELTA )
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
						if ( nodeFound ) {
							isMeasuredNode_.Push_back(true);
							//std::cout << "Node is MeasureNode" << std::endl;
						}
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
void CoefFunctionGridNodalSource<DATA_TYPE>::ComputeSourceData(Vector<DATA_TYPE>& actPDEsol) {

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
		result.Init();
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
		Double deltaAmp;
		deltaAmp = alpha_* qExp_* std::pow( std::abs( sourceAmp_[i]), qExp_-1.0 );
		if ( sourceAmp_[i] < 0.0 )
			deltaAmp *= -1.0;

		Complex jphi(0,sourcePhi_[i]); //j*phi
		Complex valC = actPDEsol[i] * std::exp(jphi);
		deltaAmp += valC.real() * scalingHesse_;

		//store derivative
		sourceAmpDelta_[i] = deltaAmp;

		//std::cout << "ApmVal: " << val << std::endl;
		valAmp += deltaAmp*deltaAmp;

		//phase
		Double deltaPsi;
		deltaPsi = 2*beta_*psi - sourceAmp_[i]*valC.imag()
						+ 2*psi / ( (M_PI/2 + psi) * (M_PI/2-psi) );

		sourcePhiDelta_[i] = deltaPsi;
		valPhi +=  deltaPsi * deltaPsi;
	}
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
		result.Init();
		this->feFunctions_[ACOU_PRESSURE]->GetEntitySolution( result, entIt );
		actPDEsol[k] = result[0];
		//std::cout << "ActSol: " << result[0] << std::endl;
		entIt++;
		k++;
	}

	error = 0.0;
	Double meanMeasP = 0.0;
	UInt idx=0;
	for(UInt i=0;i<this->fctSolAssoc_.GetSize();++i) {
		const std::pair<UInt,UInt> & curP = this->fctSolAssoc_[i];
		if ( curP.second > 0 ) {
			if ( isMeasuredNode_[i] ) {
				meanMeasP += std::abs( measVec_[i]) * std::abs( measVec_[i]) ;
				Complex val = actPDEsol[idx] - measVec_[i];
				std::cout << "Nr: " << idx << "  PDE: " << actPDEsol[idx] << " Meas:" << measVec_[i]
						<< "  Diff: " << val << std::endl;
				error += std::abs( val ) * std::abs( val ) ;
				idx++;
			}
		}
	}

    std::cout << "\n Relative L2 in %: " << 100*std::sqrt( error ) / std::sqrt( meanMeasP ) << std::endl;
//	error = std::sqrt( error ); // (Double) idx );
}



template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::UpdateSource( Double& stepLength, bool lineSearch ) {

	if ( lineSearch ) {
		for ( UInt i=0; i<sourceAmp_.GetSize(); i++ ) {
			sourceAmp_[i] = sourceAmpSave_[i] - stepLength * sourceAmpDelta_[i];
			sourcePhi_[i] = sourcePhiSave_[i] - stepLength * sourcePhiDelta_[i];
		}
	}
	else {
		sourceAmpSave_ = sourceAmp_;
		sourcePhiSave_ = sourcePhi_;
	}
}


//! compute Tikhonov function
template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::ComputeTikh( Double& funcVal, Double& resSquared,
		                                                  bool adjustAlpha, bool adjustBeta ) {

	Double valAmp = 0;
	Double valPhi = 0;

	if ( adjustAlpha )
		alpha_ = 1.0;
	if ( adjustBeta )
		beta_ = 1.0;

	for ( UInt i=0; i<sourceAmp_.GetSize(); i++ ) {
		valAmp += alpha_*std::pow( std::abs(sourceAmp_[i]), qExp_ );
		valPhi += beta_*std::pow(std::abs(sourcePhi_[i]), 2);
		//funcVal += alpha_*std::pow( std::abs(sourceAmp_[i]), qExp_ ) +	beta_*std::pow(std::abs(sourcePhi_[i]), 2);
	}

	if ( adjustAlpha ) {
		alpha_ = resSquared / valAmp *0.1;
		valAmp *= alpha_;
		std::cout << "Adjusted alpha: " << alpha_ << std::endl;
	}

	if ( adjustBeta) {
		if ( valPhi > 1e-4 )
			beta_ = resSquared / valPhi *0.1;
		valPhi *= beta_;
		std::cout << "Adjusted beta: " << beta_ << std::endl;
	}
	std::cout << "valL2: " << resSquared << "  valAmP: " <<  valAmp << "  valPhi: " << valPhi << std::endl;
	funcVal = resSquared + valAmp + valPhi;
}

template<typename DATA_TYPE>
void CoefFunctionGridNodalSource<DATA_TYPE>::ComputeMeasL2squared(Double& valL2 ) {

	Complex val(0,0);
	for(UInt i=0;i<this->fctSolAssoc_.GetSize();++i) {
		const std::pair<UInt,UInt> & curP = this->fctSolAssoc_[i];
		if ( curP.second > 0 ) {
			if ( isMeasuredNode_[i] )
				val += measVec_[i] * measVec_[i] ;
		}
	}
	valL2 = std::abs( val );
}


}

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
//  template class CoefFunctionGridNodalSource<Double>;
  template class CoefFunctionGridNodalSource<Complex>;
#endif
