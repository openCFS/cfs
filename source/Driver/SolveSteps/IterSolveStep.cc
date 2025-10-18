#include "IterSolveStep.hh"

#include <iomanip>

#include "MatVec/BaseMatrix.hh"
#include "PDE/StdPDE.hh"
#include "PDE/SinglePDE.hh"
#include "CoupledPDE/IterCoupledPDE.hh"
#include "CoupledPDE/DirectCoupledPDE.hh"
#include "Domain/CoefFunction/CoefFunctionAccumulator.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Driver/TimeSchemes/BaseTimeScheme.hh"
#include "Driver/AnalysisID.hh"
#include "Driver/BaseDriver.hh"
namespace CoupledField
{


// declare logging stream
DEFINE_LOG(itersolvestep, "itersolvestep")

// ======================================================================
//  Classes for Convergence Criterions
// ======================================================================

  ConvCriterion::ConvCriterion( NormType type, Double value ) {
    normType_ = type;
    finalNorm_ = value;
  }
  
  Double ConvCriterion::CalcNorm( Double newVal, Double oldVal ) {
    Double ret = 0.0;
    /*
     * // Tested, but did not work
    if((std::abs(newVal) > 1e30) || (std::abs(oldVal) > 1e30)){
	    std::cout << "Abs(newVal)= "<<std::abs(newVal) << "; Abs(oldVal)= "<<std::abs(oldVal) << std::endl;
	    std::cout << "Input out of bounds, reset norm to 1 to avoid NaN!" << std::endl;
	    return 1.0; 
    }
    */
    Double delta = std::abs( std::abs(newVal) - std::abs(oldVal) );

    switch( normType_) {
      case NO_NORM:
        break;
      case L2ABS:
        ret = std::abs(newVal);
        break;
      case L2REL:

        //if (std::abs(newVal) > 1.0 ) { 
        if (std::abs(newVal) > 0.0 ) {
          ret = delta / std::abs(newVal);
        } else {
          ret = delta;
        }
        break;
        
    }

    return ret;
    
    
      /* FOR comparison from OLD TRUNK  
    // Calculate difference
    delta->Add(1.0, val, -1.0, oldval );

//    Vector<Double> & val_vec = dynamic_cast<Vector<Double>& >(val);
//    Vector<Double> & oldval_vec = dynamic_cast<Vector<Double>& >(oldval);
//    delta = val_vec - oldval_vec;

    switch (normtype)
      {
      case NO_NORM:
        return 0;
        break;
      
      case L2ABS:
        norm = delta->NormL2();
        break;

      case L2REL:
        valNorm2 =  val.NormL2();
        if (valNorm2 > 0) 				<---- this is the interesting part here!
          norm = delta->NormL2() / valNorm2;
        else
          norm = delta->NormL2();

        break;
      }
*/
    
    
  }

  // Definition of norm types
  static EnumTuple normTypeTuples[] = 
  {
   EnumTuple(ConvCriterion::NO_NORM, "no"),                                    
   EnumTuple(ConvCriterion::L2REL, "rel"),
   EnumTuple(ConvCriterion::L2ABS, "abs")
  };
  Enum<ConvCriterion::NormType> ConvCriterion::NormTypeEnum = 
      Enum<ConvCriterion::NormType> ("Types of boundary conditions",
         sizeof(normTypeTuples) / sizeof(EnumTuple),
         normTypeTuples);


  ConvCriterionAccu::ConvCriterionAccu(NormType type, Double value, bool overrideNumInt )
  : ConvCriterion(type, value) {
    actNorm_ = 0.0;
    oldNorm_ = 0.0;
    overrideNumInt_ = false;
  }
  
  ConvCriterionAccu::~ConvCriterionAccu() {
    
  }
  
  void ConvCriterionAccu::AddCoefFct( shared_ptr<EntityList> list,
                                      shared_ptr<CoefFunctionAccumulator> coefFct ) {
    
   // change the integration style if necessary
   coefFct->SetIntegrateFlag(overrideNumInt_);
   
   coefs_[list] = coefFct;
   oldNorm_ = 0.0;
  }
  
  void ConvCriterionAccu::ResetValues() {
    oldNorm_ = 0.0;
    actNorm_ = 0.0;
    
  }
  
  void ConvCriterionAccu::StartSampling() {
    
    std::map<shared_ptr<EntityList>, shared_ptr<CoefFunctionAccumulator> >::iterator it;
    it = coefs_.begin();
    
    // Remember old norm
    oldNorm_ = actNorm_;
    
    // Re-set all oceffunction accumulator
    for( ; it != coefs_.end(); ++it ) {
      it->second->ResetSampling();
    }
  }
  
  void ConvCriterionAccu::StopSampling(){
    std::map<shared_ptr<EntityList>, shared_ptr<CoefFunctionAccumulator> >::iterator it;
    it = coefs_.begin();
    actNorm_ = 0.0;
    // Re-set all oceffunction accumulator
    for( ; it != coefs_.end(); ++it ) {
      Double tmp = it->second->GetNorm();
      actNorm_ +=  tmp * tmp; 
    }
    actNorm_ = sqrt(actNorm_);
  }
  
  Double ConvCriterionAccu::GetNorm() {
    return CalcNorm(actNorm_, oldNorm_);
  }
  
  bool ConvCriterionAccu::Converged() {
    if( normType_ == ConvCriterion::NO_NORM ) {
      return true;
    } else {
      return CalcNorm( actNorm_, oldNorm_) <= finalNorm_;
    }
  }
  
  StdVector<shared_ptr<EntityList> > ConvCriterionAccu::GetSupport() {
    StdVector<shared_ptr<EntityList> > ret;
    ret.Reserve( coefs_.size());
    std::map<shared_ptr<EntityList>, shared_ptr<CoefFunctionAccumulator> >::iterator it;
    for( ; it != coefs_.end(); ++it ) {
      ret.Push_back(it->first);
    }
    return ret;
  }
  // ------------------------------------------------------------------------
  
  ConvCriterionDisplacement::ConvCriterionDisplacement( NormType type, 
                                                        Double value )
  : ConvCriterion(type, value) {
    actNorm_ = 0.0;
    oldNorm_ = 0.0;
    SetNormFlag(true);
  }
  
  ConvCriterionDisplacement::~ConvCriterionDisplacement() {
    
  }
  
  void ConvCriterionDisplacement::SetDispFct( shared_ptr<FeFunction<Double> > disp) {
    assert(disp);
    disp_ = disp;
  }

  void ConvCriterionDisplacement::SetDispFctComplex( shared_ptr<FeFunction<Complex> > disp) {
    assert(disp);
    dispComplex_ = disp;
  }
  
  void ConvCriterionDisplacement::SetVelFct( shared_ptr<FeFunction<Double> > vel) {
    assert(vel);
    vel_ = vel;
  }

  void ConvCriterionDisplacement::SetAccFct( shared_ptr<FeFunction<Double> > acc) {
    assert(acc);
    acc_ = acc;
  }

  void ConvCriterionDisplacement::SetNormFlag( bool justNorm) {
    justNorm_ = justNorm;
  }  
  
  
  void ConvCriterionDisplacement::AddRegion(RegionIdType region ) {
    updatedRegions_.insert(region);
  }
  
  
  void ConvCriterionDisplacement::ResetValues() {
      oldNorm_ = 0.0;
      actNorm_ = 0.0;
      
    }
  
  void ConvCriterionDisplacement::StartSampling() {
    // update grid to current values
      // Check if displacement fefunction is set
    // if we just want to calculate the norm, we need no geometry update
    if( !(disp_ || dispComplex_) || justNorm_)
        return;
/*    
    double currentNorm = CalcNorm(actNorm_, oldNorm_);
    if(currentNorm >= 1.0){
	   std::cout << "Difference to last timestep too large (currentNorm = " << currentNorm << " )" << std::endl;
	   std::cout << "Skip mesh smoothing this time! " << std::endl;
	   return;
   } 
*/	    
      Grid * ptGrid;
      if( isComplex_ ) {
        ptGrid = dispComplex_->GetGrid();
      } else {
        ptGrid = disp_->GetGrid();
      }
      const UInt dim = ptGrid->GetDim();
      // Grid for vel and acc
      //Grid * ptGridVel = vel_->GetGrid();
      //Grid * ptGridAcc = acc_->GetGrid();
      // Loop over all regions of FeFunction
      shared_ptr<EntityList> nodes;
      std::set<RegionIdType> dispRegions;
      if( isComplex_ ) {
        dispRegions = dispComplex_->GetRegions();

        // get the potential reference node and calcualte the phase offset
        if( refNodeEnabled_ ) {
          shared_ptr<EntityList> refNodeList;
          refNodeList = ptGrid->GetEntityList(EntityList::NODE_LIST, refNodeName_);
          if( refNodeList->GetSize() != 1 ) {
            EXCEPTION("The reference node name specified either does not exist or contains more than one node!");
          }
          EntityIterator refNodeIt = refNodeList->GetIterator();
          Vector<Complex> refNodeOffset(dim);
          dispComplex_->GetEntitySolution(refNodeOffset, refNodeIt);

          // calculate the phase and phase correction
          Vector<double> refNodePhase(dim);
          for( UInt ii=0; ii<refNodePhase.GetSize(); ii++) {
            refNodePhase[ii] = std::arg(refNodeOffset[ii]);
          }
          phaseCorrection_ = -refNodePhase[refNodeDOF_]+phaseOffset_/180.0*3.14159265358979323846;
          // this complex value can be used to directly multiply it with the result, leading to a zero phase for the specified DOF
          // the rest will be projected by using the real value
          phaseCorrMult_ = std::polar(1.0,phaseCorrection_);

          // calculate the actual offset of the reference node and give feedback to the user
          Vector<Complex> refNodeOffsetCorr(dim);
          refNodeOffsetCorr = refNodeOffset*phaseCorrMult_;
          std::cout << "Found reference node " << refNodeName_ << "!" << std::endl;
          std::cout << "Initial position: " << refNodeOffset.ToString() << std::endl;
          std::cout << "Reference calculated based on DOF " << refNodeDOF_ << std::endl;
          std::cout << "Corrected position: " << refNodeOffsetCorr.ToString() << std::endl;
        }
      } else {
        dispRegions = disp_->GetRegions();
      }

      // dummy to compute the offset of the current node and DOF based on the phase correction multiplier
      Complex curOffset;
      
      std::set<RegionIdType>::const_iterator regionIt = updatedRegions_.begin();

      for( ; regionIt != updatedRegions_.end(); regionIt++ ) {

        std::string regionName = ptGrid->GetRegion().ToString(*regionIt);

        // check if this region is contained in the displacement function as well
        if( dispRegions.find(*regionIt) == dispRegions.end() ) {
          WARN( "Can not perform geometry update on region"
              << regionName << ", as there are no displacement defined on it!");
        }

        nodes = ptGrid->GetEntityList(EntityList::NODE_LIST, regionName);
        EntityIterator nodeIt = nodes->GetIterator();

        // Loop over all nodes
        Vector<Double> offset(dim), totalOffset(nodes->GetSize() * dim );
        StdVector<UInt> nodeNums(nodes->GetSize());
        UInt pos = 0;
        for( ; !nodeIt.IsEnd(); nodeIt++, pos++ ) {

          nodeNums[pos] = nodeIt.GetNode();
          // aquire nodal solution
          if( isComplex_ ) {
            Vector<Complex> offsetComplex(dim);
            dispComplex_->GetEntitySolution(offsetComplex, nodeIt);
            for( UInt i=0; i<offsetComplex.GetSize(); i++ ) {
              curOffset = offsetComplex[i]*phaseCorrMult_;
              offset[i] = curOffset.real();
            }
          } else {
            disp_->GetEntitySolution(offset, nodeIt);
          }

          UInt offsetPos = pos*dim;
          for( UInt iDim = 0; iDim < dim; ++iDim ) {
            totalOffset[offsetPos+iDim] = offset[iDim];
          }

        }
        // Pass total array
        ptGrid->SetNodeOffset(nodeNums, totalOffset);
        // We do the same for the velocity and acceleration since they have their own FeFunction
        //ptGridVel->SetNodeOffset(nodeNums, totalOffset);
        //ptGridAcc->SetNodeOffset(nodeNums, totalOffset);

      }
      // update nc interfaces if existing
      ptGrid->MoveNcInterfaces();
      //ptGridVel->MoveNcInterfaces();
      //ptGridAcc->MoveNcInterfaces();
  }
  
  void ConvCriterionDisplacement::StopSampling() {

    // if no displacement is set, just leave
    if( !(disp_ || dispComplex_))
      return;

    oldNorm_ = actNorm_;
    actNorm_ = 0.0;

    // Calculate norm of total displacement
    Grid * ptGrid;
    if( isComplex_ ) {
      ptGrid = dispComplex_->GetGrid();
    } else {
      ptGrid = disp_->GetGrid();
    }
    
    // update nc interfaces if existing
    //ptGrid->MoveNcInterfaces();
    
    const UInt dim = ptGrid->GetDim();
    // Loop over all regions of FeFunction
    shared_ptr<EntityList> nodes;
    std::set<RegionIdType> dispRegions;
    if( isComplex_ ) {
      dispRegions = dispComplex_->GetRegions();
    } else {
      dispRegions = disp_->GetRegions();
    }
    
    // dummy to compute the offset of the current node and DOF based on the phase correction multiplier
    Complex curOffset;
    
    std::set<RegionIdType>::const_iterator regionIt = updatedRegions_.begin();

    for( ; regionIt != updatedRegions_.end(); regionIt++ ) {
      std::string regionName = ptGrid->GetRegion().ToString(*regionIt);

      // check if this region is contained in the displacement function as well
      if( dispRegions.find(*regionIt) == dispRegions.end() ) {
        WARN( "Can not perform geometry update on region"
            << regionName << ", as there are no displacement defined on it!");
      }

      nodes = ptGrid->GetEntityList(EntityList::NODE_LIST, regionName);
      EntityIterator nodeIt = nodes->GetIterator();

      // Loop over all nodes
      Vector<Double> offset(dim);
      StdVector<UInt> nodeNums(nodes->GetSize());
      UInt pos = 0;
      for( ; !nodeIt.IsEnd(); nodeIt++, pos++ ) {

        nodeNums[pos] = nodeIt.GetNode();
        // aquire nodal solution
        if( isComplex_ ) {
          Vector<Complex> offsetComplex(dim);
          dispComplex_->GetEntitySolution(offsetComplex, nodeIt);
          for( UInt i=0; i<offsetComplex.GetSize(); i++ ) {
            curOffset = offsetComplex[i]*phaseCorrMult_;
            offset[i] = curOffset.real();
          }
        } else {
          disp_->GetEntitySolution(offset, nodeIt);
        }
        
        for( UInt iDim = 0; iDim < dim; ++iDim ) {
          actNorm_ +=  offset[iDim] * offset[iDim];
        }

      } // loop: nodes
    } // loop: regions
    
    // take square root 
    actNorm_ = sqrt(actNorm_);
      
  }

  Double ConvCriterionDisplacement::GetNorm() {
    return CalcNorm( actNorm_, oldNorm_);
  
  }
      
  bool ConvCriterionDisplacement::Converged() {
    if( normType_ == ConvCriterion::NO_NORM ) {
      return true;
    } else {
      return CalcNorm( actNorm_, oldNorm_) <= finalNorm_;
    }
  }
  
  StdVector<shared_ptr<EntityList> > ConvCriterionDisplacement::GetSupport() {
    StdVector<shared_ptr<EntityList> > ret;
    return ret;
  }


  
// ======================================================================
  //! Derived class for step-wise solving of iterative coupled StdPDEs

  IterSolveStep::IterSolveStep(IterCoupledPDE &apde, PtrParamNode node,
                               PtrParamNode infoNode) 
    : BaseSolveStep(),
      rPDE_(apde)
  {
    param_ = node;
    info_ = infoNode;
    startStep_ = 1;
    isFinalized_ = false;
    nonLinLogging_ = false;
    stopOnDivergence_ =  false;
    maxiter_ = 0;
    actAnalysisType_ = BasePDE::NO_ANALYSIS;
    // for mechPDE if only the norm shall converge, but no geometry change shall be executed
    justNorm_ = false;
    
    // use user defined PDE order
    customReorderPDE_  = false;

    // Initialize solution map
    solutionMap_[MAG_FORCE_LORENTZ_DENSITY] = MAG_FORCE_LORENTZ;
    solutionMap_[MAG_FORCE_MAXWELL_DENSITY] = MAG_FORCE_MAXWELL;
    solutionMap_[ELEC_POWER_DENSITY] = ELEC_POWER;
  }
  
    
   

  IterSolveStep::~IterSolveStep()
  {
  }


  void IterSolveStep::Init() {

    LOG_TRACE(itersolvestep) << "Initializing iterative coupling";
    
    // fetch "convergence" node
    PtrParamNode convParamNode; 
    if( param_ )
      convParamNode = param_->Get("convergence", ParamNode::PASS );
    convNode_ = info_->Get("convergence");

    // get maximum number of iterations (optional)
    maxiter_ = 1;
    if( convParamNode )
      convParamNode->GetValue( "maxNumIters", maxiter_, ParamNode::PASS );
    convNode_->Get("maxNumIters")->SetValue(maxiter_);

    // query logging flag
    nonLinLogging_ = true;
    if( convParamNode ) 
      convParamNode->GetValue( "logging", nonLinLogging_, ParamNode::PASS );
    convNode_->Get("logging")->SetValue(nonLinLogging_);
    
    // query divergence behavior
    stopOnDivergence_ = true;
    if( convParamNode ) {
      convParamNode->GetValue( "stopOnDivergence", stopOnDivergence_, ParamNode::PASS );
    }
    convNode_->Get("stopOnDivergence")->SetValue(stopOnDivergence_);

    justNorm_ = true;
    if( convParamNode ) {
      convParamNode->GetValue( "justNorm", justNorm_, ParamNode::PASS );
    } 
    convNode_->Get("justNorm")->SetValue(justNorm_);

    // check for custom PDE order
    if( param_ ) {
      // use the prescribed order of the PDEs
      PDEorder_ = param_->Get("PDEorder")->As<std::string>();
      if ( !PDEorder_.empty() ) {
        customReorderPDE_ = true;
      }
      endWithCoupledPDEs_ = param_->Get("endWithCoupledPDEs")->As<bool>();
    }

    // 1) Check for general convergence criterions
    if( convParamNode && convParamNode->Has("quantity") ) {
      LOG_TRACE(itersolvestep) << "Checking convergence criterions";
      ParamNodeList q = convParamNode->GetList("quantity");
      for( UInt i = 0; i < q.GetSize(); ++i ) {
        std::string quantity = q[i]->Get("name")->As<std::string>();

        // For generic result, a simple SolutionTypeEnum.Parse(quantity) might fail because
        // they are not defined yet. Hence, we check if a postProcList is given and if we
        // can find and define the quantity
        // check if we can parse the quantity (by passing the false flag, we just get an
        // invalid solution type but not an error)
        // we do this before checking for generic postProcessing results (although they might not be defined already)
        // since they COULD be defined already be a PDE read before this one
        SolutionType solType = SolutionTypeEnum.TryParse(quantity, INVALID_SOLUTION_TYPE);
        if( solType == INVALID_SOLUTION_TYPE ) {
          // direct conversion was not successful, check if we can add the result based on postProcessing results

          // get the postProcList and access all names of type function
          PtrParamNode ppListNode = param_->GetParent()->GetParent()->Get("postProcList",ParamNode::PASS);
          std::string funcName;

          if( ppListNode ) {
            ParamNodeList ppListNodeChildren = ppListNode->GetChildren();

            // loop over all postProc definitions and check, if a function is defined
            for( UInt i = 0; i < ppListNodeChildren.GetSize(); i++ ) {

              // get all children of one postProc definition
              ParamNodeList ppNodeChildren = ppListNodeChildren[i]->GetChildren();
              
              // we only consider the function type postProcResults
              ParamNodeList postProcFuncs = ParamNode::GetListByChild(ppNodeChildren, "function");

              // loop over all function-type children and check their name
              for( UInt u = 0; u < postProcFuncs.GetSize(); u++ ) {
                postProcFuncs[u]->GetValue("resultName",funcName);

                if( funcName==quantity ) {
                  // the quantity we need is defined as a postProcResult but not yet parseable
                  // add it to the environment and get a usable solution type
                  AddGenericSolution(quantity, this->rPDE_.GetDomain());
                }
              }
            }
          }
          // finally, try to parse the quantity again
          // if it does not work, the quantity is not defined (neither pre-defined nor via the postProcList)
          solType = SolutionTypeEnum.Parse(quantity);
        }
        
        Double norm = q[i]->Get("value")->As<Double>();
        ConvCriterion::NormType type = 
            ConvCriterion::NormTypeEnum.Parse(q[i]->Get("normType")
                                              ->As<std::string>());
        bool overrideNumInt = q[i]->Get("avoidIntegratedNorm")->As<bool>();
        LOG_DBG3(itersolvestep) << "\tQuantity: " << quantity;
        LOG_DBG3(itersolvestep) << "\tNorm:     " << norm;
        LOG_DBG3(itersolvestep) << "\tNormType: " << ConvCriterion::NormTypeEnum.ToString(type);
        LOG_DBG3(itersolvestep) << "\tavoidIntegratedNorm: " << overrideNumInt;
        // Create new convergence criterion, depending on solution type
        shared_ptr<ConvCriterion> crit;
        // test: if no geometry update is present, calculate mech displacement as normal criterion

        if( solType == MECH_DISPLACEMENT ) {
          LOG_DBG3(itersolvestep) << "\t=> Creating special displacement convergence criterion";
          shared_ptr<ConvCriterionDisplacement> accu (new ConvCriterionDisplacement(type, norm));
          accu->SetNormFlag(justNorm_);
          crit  = accu;
        } else if( solType == SMOOTH_DISPLACEMENT ) {
          LOG_DBG3(itersolvestep) << "\t=> Creating special displacement convergence criterion (smoothPDE)";
          shared_ptr<ConvCriterionDisplacement> accu (new ConvCriterionDisplacement(type, norm));
          accu->SetNormFlag(justNorm_);
          crit  = accu;
        } else {
          LOG_DBG3(itersolvestep) << "\t=> Creating general accumulated convergence criterion";
          shared_ptr<ConvCriterionAccu> accu (new ConvCriterionAccu(type, norm, overrideNumInt));
          crit = accu;
        }
      criterions_[solType] = crit;
      } // loop: quantities
    }
  }
  
  
  void IterSolveStep::TriggerFinalize() {
    if( !isFinalized_) {
      LOG_DBG(itersolvestep) << "Calling ::Finalize() from TriggerFinalize()";
      Finalize();
    } 
  }
  
  void IterSolveStep::Finalize() {
    LOG_TRACE(itersolvestep) << "Finalizing iterative coupled solve step";
    
    // 1) Check for updated geometry
    if( param_->Has("geometryUpdate") ) {
      ParamNodeList regionNodes = param_->Get("geometryUpdate")->GetList("region");

      bool enableRefNode = false;
      PtrParamNode geomUpdateParamNode = param_->Get("geometryUpdate", ParamNode::PASS );
      geomUpdateParamNode->GetValue( "EnableRefNode", enableRefNode, ParamNode::PASS );

      std::string refNodeName;
      double phaseOffset = 0.0;
      UInt refNodeDOF = 0;
      if( enableRefNode ) {
        // search for the named node, get its ID and use it as a reference
        refNodeName = geomUpdateParamNode->Get("RefNodeName")->As<std::string>();
        phaseOffset = geomUpdateParamNode->Get("RefNodePhaseOffset")->As<double>();
        refNodeDOF = geomUpdateParamNode->Get("RefNodeDOFNumber")->As<UInt>();
        // we get the actual node number later on when we have access to the grid pointer
      }

      if( regionNodes.GetSize() > 0 ) {

        // check for presence of mechanical PDE
        shared_ptr<FeFunction<Double> > disp;
        shared_ptr<FeFunction<Complex> > dispComplex;
        //shared_ptr<FeFunction<Double> > vel;
        //shared_ptr<FeFunction<Double> > acc;
        // check for presence of smooth PDE
        shared_ptr<FeFunction<Double> > dispSmooth;
        shared_ptr<FeFunction<Complex> > dispSmoothComplex;
        //shared_ptr<FeFunction<Double> > velSmooth;
        UInt numSinglePDEs = rPDE_.singlePDEs_.GetSize();
        for( UInt i = 0; i < numSinglePDEs; ++i ) {
          // since we have to differentiate between mech- and smooth displacement, we include it in the lopp
          // furthermore, we have to distinguish which region uses which criterion
          SinglePDE * ptPde = rPDE_.singlePDEs_[i]; 
          if( ptPde->GetName() == "mechanic" ) {
            // enable harmonic case
            if ( ptPde->IsComplex() ) {
              dispComplex = dynamic_pointer_cast<FeFunction<Complex> >
                            (ptPde->GetFeFunction(MECH_DISPLACEMENT));
            } else {
              disp = dynamic_pointer_cast<FeFunction<Double> >
                            (ptPde->GetFeFunction(MECH_DISPLACEMENT));
            }
            LOG_DBG(itersolvestep) << "=> Found MECH_DISPLACEMENT as coupling quantity";
//            vel = dynamic_pointer_cast<FeFunction<Double> >
//                        (ptPde->GetFeFunction(MECH_VELOCITY));
//            acc = dynamic_pointer_cast<FeFunction<Double> >
//                        (ptPde->GetFeFunction(MECH_ACCELERATION));

            // Check if convergence criterion for mechanic is present
            shared_ptr<ConvCriterionDisplacement> convDisp;
            if( criterions_.find(MECH_DISPLACEMENT) != criterions_.end()) {
              convDisp =
                  dynamic_pointer_cast<ConvCriterionDisplacement>(criterions_[MECH_DISPLACEMENT]);
            } else {
              convDisp.reset(new ConvCriterionDisplacement(ConvCriterion::NO_NORM, 0.0));
            }

            if ( ptPde->IsComplex() ) {
              WARN("You specified a geometry update for a harmonic simulation. Only the amplitude of the deformation will be used!");
              convDisp->SetDispFctComplex( dispComplex );
              convDisp->SetIsComplex(true);
              convDisp->SetEnableRefNode(enableRefNode);
              convDisp->SetRefNodeName(refNodeName);
              convDisp->SetRefNodeDOF(refNodeDOF);
              convDisp->SetPhaseOffset(phaseOffset);
            } else {
              convDisp->SetDispFct( disp );
            }
            // We set the vel and acc as well since they have their own FeFunction and need the geometry update too
            //convDispSmooth->SetVelFct( velSmooth );
            //convDisp->SetAccFct( acc );
            Grid * ptGrid;
            if ( ptPde->IsComplex() ) {
              ptGrid = dispComplex->GetGrid();
            } else {
              ptGrid = disp->GetGrid();
            }

            LOG_DBG(itersolvestep) << "Performing geometry update on the following regions:";
            // Read in all regions, which have geometric update and check if they are present in the mechPDE
            for( UInt i = 0; i < regionNodes.GetSize(); ++i ) {
              std::string regionName = regionNodes[i]->Get("name")->As<std::string>();
              RegionIdType regionId = ptGrid->GetRegion().Parse(regionName);
              StdVector<RegionIdType> regionsPde;
              regionsPde = ptPde->GetRegions();
              if( std::find(regionsPde.Begin(), regionsPde.End(), regionId )!= regionsPde.End() ) {
                // the region is defined for this PDE, set the convergence criterion
                convDisp->AddRegion( regionId );
                LOG_DBG(itersolvestep) << "\t Region for mechPDE: "<< regionName;
              }
            }

          } else if( ptPde->GetName() == "smooth" ) {
            // enable harmonic case
            if ( ptPde->IsComplex() ) {
              dispSmoothComplex = dynamic_pointer_cast<FeFunction<Complex> >
                            (ptPde->GetFeFunction(SMOOTH_DISPLACEMENT));
            } else {
              dispSmooth = dynamic_pointer_cast<FeFunction<Double> >
                            (ptPde->GetFeFunction(SMOOTH_DISPLACEMENT));
            }
            LOG_DBG(itersolvestep) << "=> Found SMOOTH_DISPLACEMENT as coupling quantity";
//            velSmooth = dynamic_pointer_cast<FeFunction<Double> >
//                        (ptPde->GetFeFunction(SMOOTH_VELOCITY));

            // Check if convergence criterion for mechanic is present
            shared_ptr<ConvCriterionDisplacement> convDispSmooth;
            if( criterions_.find(SMOOTH_DISPLACEMENT) != criterions_.end()) {
              convDispSmooth =
                  dynamic_pointer_cast<ConvCriterionDisplacement>(criterions_[SMOOTH_DISPLACEMENT]);
            } else {
              convDispSmooth.reset(new ConvCriterionDisplacement(ConvCriterion::NO_NORM, 0.0));
            }

            if ( ptPde->IsComplex() ) {
              WARN("You specified a geometry update for a harmonic simulation. Only the amplitude of the deformation will be used!");
              convDispSmooth->SetDispFctComplex( dispSmoothComplex );
              convDispSmooth->SetIsComplex(true);
              convDispSmooth->SetEnableRefNode(enableRefNode);
              convDispSmooth->SetRefNodeName(refNodeName);
              convDispSmooth->SetRefNodeDOF(refNodeDOF);
              convDispSmooth->SetPhaseOffset(phaseOffset);
            } else {
              convDispSmooth->SetDispFct( dispSmooth );
            }
            // We set the vel and acc as well since they have their own FeFunction and need the geometry update too
            //convDispSmooth->SetVelFct( velSmooth );
            //convDisp->SetAccFct( acc );
            Grid * ptGrid;
            if ( ptPde->IsComplex() ) {
              ptGrid = dispSmoothComplex->GetGrid();
            } else {
              ptGrid = dispSmooth->GetGrid();
            }

            LOG_DBG(itersolvestep) << "Performing geometry update on the following regions:";
            // Read in all regions, which have geometric update and check if they are present in the mechPDE
            for( UInt i = 0; i < regionNodes.GetSize(); ++i ) {
              std::string regionName = regionNodes[i]->Get("name")->As<std::string>();
              RegionIdType regionId = ptGrid->GetRegion().Parse(regionName);
              StdVector<RegionIdType> regionsPde;
              regionsPde = ptPde->GetRegions();
              if( std::find(regionsPde.Begin(), regionsPde.End(), regionId )!= regionsPde.End() ) {
                // the region is defined for this PDE, set the convergence criterion
                convDispSmooth->AddRegion( regionId );
                LOG_DBG(itersolvestep) << "\t Region for smoothPDE: "<< regionName;
              }
            }

          }

        }
        if(!disp && !dispComplex && !dispSmooth && !dispSmoothComplex) {
          WARN( "No geometry updated will performed, as no mechanical "
              << "physic is defined");
        } else {
//          // Check if convergence criterion for mechanic is present
//          shared_ptr<ConvCriterionDisplacement> convDisp;
//          if( criterions_.find(MECH_DISPLACEMENT) != criterions_.end()) {
//            convDisp =
//                dynamic_pointer_cast<ConvCriterionDisplacement>(criterions_[MECH_DISPLACEMENT]);
//          } else if( criterions_.find(SMOOTH_DISPLACEMENT) != criterions_.end()) {
//            convDisp =
//                dynamic_pointer_cast<ConvCriterionDisplacement>(criterions_[SMOOTH_DISPLACEMENT]);
//          } else {
//            convDisp.reset(new ConvCriterionDisplacement(ConvCriterion::NO_NORM, 0.0));
//          }
//
//          convDisp->SetDispFct( disp );
//          // We set the vel and acc as well since they have their own FeFunction and need the geometry update too
//          convDisp->SetVelFct( vel );
//          //convDisp->SetAccFct( acc );
//          Grid * ptGrid = disp->GetGrid();
//
//          LOG_DBG(itersolvestep) << "Performing geometry update on the following regions:";
//          // Read in all regions, which have geometric update
//          for( UInt i = 0; i < regionNodes.GetSize(); ++i ) {
//            std::string regionName = regionNodes[i]->Get("name")->As<std::string>();
//            RegionIdType regionId = ptGrid->GetRegion().Parse(regionName);
//            convDisp->AddRegion( regionId );
//            LOG_DBG(itersolvestep) << "\t "<< regionName;
//          }
        }

      }
    }
    
    // 2) Loop over all convergence criterions, which are not explicitly "set"
    // and try to get CoefFunction from a related PDE.
    // ... to be implemented in a further step
    
    // 3) Resort the PDE order 
    ResortPDEOrder();
    isFinalized_ = true;
  }
  
  void IterSolveStep::ResortPDEOrder() {
    LOG_TRACE(itersolvestep) << "Resorting PDE order";

    // Collect all uncoupled SinglePDes
    std::set<SinglePDE*> uncoupledPdes;
    uncoupledPdes.insert(rPDE_.singlePDEs_.Begin(), 
                         rPDE_.singlePDEs_.End() );

    // loop over all coupled SinglePDEs and remove involved
    // SinglePDEs
    for( UInt i = 0; i < rPDE_.coupledPDEs_.GetSize(); ++i ) {
      StdVector<SinglePDE*> sPdes = rPDE_.coupledPDEs_[i]->GetSinglePDEs();
      for( UInt j = 0; j < sPdes.GetSize(); ++j ) {
        SinglePDE * singlePde = sPdes[j];
        if( uncoupledPdes.find(singlePde) != uncoupledPdes.end() ) {
          uncoupledPdes.erase(singlePde);
        }
      }
    }
    
    // Now assemble the list of StdPDEs. Currently we pursue the 
    // following strategy:
    // 1) We start by all uncoupled Pdes
    // 2) Add coupled Pdes in the end
    rPDE_.numPDEs_ = uncoupledPdes.size() + rPDE_.coupledPDEs_.GetSize();

    rPDE_.PDEs_.Reserve( rPDE_.numPDEs_ );



    if ( customReorderPDE_ ) {
      // use the prescribed order of the PDEs
      typedef boost::tokenizer<boost::char_separator<char> > Tok;
      boost::char_separator<char> sep(";| ");
      Tok t(PDEorder_, sep);
      Tok::iterator iter, end;
      iter = t.begin();
      end = t.end();

      if( !endWithCoupledPDEs_ ) {
        // append the coupled PDEs
        for( UInt i = 0; i < rPDE_.coupledPDEs_.GetSize(); ++i ) {
          rPDE_.PDEs_.Push_back( rPDE_.coupledPDEs_[i] );
        }
      }

      std::set<SinglePDE*>::iterator     it = uncoupledPdes.begin();
      for( ; iter != end; iter++) {
        // loop over all specified PDEs and append them
        it = uncoupledPdes.begin();
        for( ; it != uncoupledPdes.end(); ++it ) {
          if( (*it)->GetName() == *iter ) {
            // we found a match, append the PDE
            rPDE_.PDEs_.Push_back( *it );
            break;
          } else {
            if( iter==end ) {
              // PDE does not exist or the name is wrong
              EXCEPTION("PDE " << (*iter) << " does not exist." );
            }
          }
        }
      }
      // now we loop over the PDEs and append those not contained in the specified PDE order
      // re-init iterators
      it = uncoupledPdes.begin();
      bool skipPDE;
      for( ; it != uncoupledPdes.end(); ++it ) {
        skipPDE = 0;
        iter = t.begin();
        for( ; iter != end; iter++) {
          if( (*it)->GetName() == *iter ) {
            // we found the match again and have to skip it
            skipPDE = 1;
          }
        }
        if( !skipPDE ) {
          rPDE_.PDEs_.Push_back( *it );
        }
      }

      if( endWithCoupledPDEs_ ) {
        // append the coupled PDEs
        for( UInt i = 0; i < rPDE_.coupledPDEs_.GetSize(); ++i ) {
          rPDE_.PDEs_.Push_back( rPDE_.coupledPDEs_[i] );
        }
      }
    } else {
      // use the classical reordering scheme
      std::set<SinglePDE*>::iterator it = uncoupledPdes.begin();
      // remember mechanic PDE if present
      SinglePDE * mechPDE = NULL;
      SinglePDE * heatPDE = NULL;
      SinglePDE * smoothPDE = NULL;
      for( ; it != uncoupledPdes.end(); ++it ) {
        if( (*it)->GetName() == "mechanic" ) {
          mechPDE = *it;
        }
        else if ( (*it)->GetName() == "heatConduction" ) {
          heatPDE = *it;
        }
        else if ( (*it)->GetName() == "smooth" ) {
                 smoothPDE = *it;
               }
        else {
          rPDE_.PDEs_.Push_back( *it );
        }
      }

      if ( heatPDE )
        rPDE_.PDEs_.Push_back( heatPDE );

      if( mechPDE )
        rPDE_.PDEs_.Push_back(mechPDE);

      if( smoothPDE )
        rPDE_.PDEs_.Push_back(smoothPDE);

      for( UInt i = 0; i < rPDE_.coupledPDEs_.GetSize(); ++i ) {
        rPDE_.PDEs_.Push_back( rPDE_.coupledPDEs_[i] );
      }
    }
    
    // In the end print final ordering:
    LOG_DBG(itersolvestep) << "Final ordering of PDEs:";
    for( UInt i = 0; i < rPDE_.numPDEs_; ++i ) {
      LOG_DBG(itersolvestep) << "\t" << i+1 << ": " 
          <<  rPDE_.PDEs_[i]->GetName();
    }
  }
  
  PtrCoefFct IterSolveStep::GetCouplingCoefFct( SolutionType type,
                                                shared_ptr<EntityList>  list,
                                                const std::string& pdeName,
                                                bool& updateGeo ) {
    LOG_TRACE(itersolvestep) << "Returning Coupling CoefFct for quantity '"
        << SolutionTypeEnum.ToString(type) << "' of PDE '" << pdeName 
        << "' on entityList '" << list->GetName() << "'";

    // finalize, if not yet done
    
    if( !isFinalized_) {
      LOG_DBG(itersolvestep) << "Calling ::Finalize()";
      Finalize();
    } 
     PtrCoefFct ret, coef;
     
     
    // Initial implementation: Directly access CoefFct of PDE
    // Later we use the interface, which keeps track of the norm of the coupling quantity
    for( UInt i = 0; i < rPDE_.singlePDEs_.GetSize(); ++i ) {
  
      if( rPDE_.singlePDEs_[i]->GetName() == pdeName ) {

        coef = rPDE_.singlePDEs_[i]->GetCoefFct(type);

        updateGeo = rPDE_.singlePDEs_[i]->IsUpdatedGeo();
        if( coef ) {
          LOG_DBG(itersolvestep) << "\t=> Found quantity, updateGeo: " << 
                      (updateGeo ? "yes" : "no");
          break;
        }
      }
    }
    if( !coef ) {
      EXCEPTION( "Could not return coupling quantity '" 
        << SolutionTypeEnum.ToString(type) << "' for Physic '"
        << pdeName << "' on entityList '" << list->GetName() << "'");
    }

    // wrap the return coefficient function in a CoefFunctionAccumulator
    shared_ptr<CoefFunctionAccumulator> acc(new CoefFunctionAccumulator(coef, true));
    
    // If this is a density quantity (e.g. force density), a possible convergence
    // criterion might be defined in terms of the absolute force (e.g. force), so we initially try
    // to find the derived value
    SolutionType mappedType = type;
    if( solutionMap_.find(type) != solutionMap_.end() ) {
      mappedType = solutionMap_[type]; 
      LOG_DBG(itersolvestep) << "\tRe-map solution type  to " <<
          SolutionTypeEnum.ToString(mappedType);
    }

    // Check, if there was a convergence criterion defined for this quantity.
    // In this case, we add it to the ConvergenceCriterion instance. If we use the displacement
    // as Dirichlet BC, we have to skip it here since its already defined for the whole region
    if ( SolutionTypeEnum.ToString(mappedType) == "mechDisplacement" ) {
      LOG_DBG(itersolvestep) << "\tQuantity is associated to convergence criterion but was skipped (mechDisplacement)";
    } else if ( SolutionTypeEnum.ToString(mappedType) == "smoothDisplacement" ) {
      LOG_DBG(itersolvestep) << "\tQuantity is associated to convergence criterion but was skipped (smoothDisplacement)";
    } else {
      if( criterions_.find(mappedType) != criterions_.end() ) {

        LOG_DBG(itersolvestep) << "\tQuantity is associated to convergence criterion";
        // add accumulated coefficient function to list
        shared_ptr<ConvCriterionAccu> c
        = dynamic_pointer_cast<ConvCriterionAccu>(criterions_[mappedType]);
        c->AddCoefFct( list, acc );
      }
    }
    
    return acc;
  }

  void IterSolveStep::GetUpdateGeoForPDE( SolutionType type,
                                                shared_ptr<EntityList>  list,
                                                const std::string& pdeName,
                                                bool& updateGeo ) {
    LOG_TRACE(itersolvestep) << "Returning updateGeo "
         << " of PDE '" << pdeName 
        << "' on entityList '" << list->GetName() << "'";
     
    // Loop over PDEs and return the flag
    for( UInt i = 0; i < rPDE_.singlePDEs_.GetSize(); ++i ) {
      if( rPDE_.singlePDEs_[i]->GetName() == pdeName ) {
        updateGeo = rPDE_.singlePDEs_[i]->IsUpdatedGeo();
      }
    }
  }
  
  // ========================================================================
  //  STATIC COUPLED ITERATION
  // ========================================================================

  void IterSolveStep::SolveStepStatic() {
    LOG_TRACE(itersolvestep) << "----------------------";
    LOG_TRACE(itersolvestep) << " Solving static step  ";
    LOG_TRACE(itersolvestep) << "----------------------\n";
    
    UInt iter = 0;
    bool normsReached = false;
    std::map<SolutionType, shared_ptr<ConvCriterion> >::iterator convIt;

    PtrParamNode actNode = convNode_->Get("step",ParamNode::APPEND); 
    actNode->Get("number")->SetValue(actStep_);
    if (nonLinLogging_) {
      for( convIt = criterions_.begin(); 
          convIt != criterions_.end(); ++convIt ) {
        std::string quantityName = SolutionTypeEnum.ToString(convIt->first);
        actNode->Get("quantity")->Get("name")->SetValue(quantityName);
      }
    }
    
    // Create a stream for logging the convergence. It will be used in the end,
    // if no convergence was achieved to report the status of the lastest iteration.
    std::stringstream msg;
    UInt width[4] = {20, 10, 15, 15}; // define widths for convergence output
    while (iter < maxiter_ &&  (! normsReached)) {
      LOG_DBG(itersolvestep) << "\n";
      LOG_DBG(itersolvestep) << "=== Iteration #" << iter+1 << "===";

      // --------------------------------------
      //  1) Re-Set all convergence criterions
      // --------------------------------------
      LOG_DBG(itersolvestep) << "Calling StartSampling for criterions";
      for( convIt = criterions_.begin(); 
          convIt != criterions_.end(); ++convIt ) {
        LOG_DBG3(itersolvestep) << "\t" << SolutionTypeEnum.ToString(convIt->first);
        convIt->second->StartSampling();
      }

      // -----------------------------------
      //  2) Calculate Sinlge PDEs
      // -----------------------------------
      for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {

        LOG_DBG(itersolvestep) << "Processing PDE '" << 
            rPDE_.PDEs_[i]->GetName() << "'";

        rPDE_.PDEs_[i]->GetSolveStep()->SetActTime(actTime_);
        rPDE_.PDEs_[i]->GetSolveStep()->SetActStep(actStep_);
        rPDE_.PDEs_[i]->GetSolveStep()->SetCouplingIter(iter);
        rPDE_.PDEs_[i]->GetSolveStep()->PreStepStatic();
        rPDE_.PDEs_[i]->GetSolveStep()->SolveStepStatic();
        rPDE_.PDEs_[i]->GetSolveStep()->PostStepStatic();
      } // end of for-loop


      // -----------------------------------
      //  3) Compute Coupling Criterions
      // -----------------------------------
      normsReached = true;
      msg.str(""); // clear logging information stream
      msg << std::setw(width[0]) << "Quantity"
          << std::setw(width[1]) << "Converged"
          << std::setw(width[2]) << "Norm"
          << std::setw(width[3]) << "Goal" << std::endl
          << std::setw(width[0]+width[1]+width[2]+width[3]);
      msg << std::setfill('-') << "" << std::setfill(' ') << std::endl;
      
      LOG_DBG(itersolvestep) << "Calling StopSampling for criterions";
      for( convIt = criterions_.begin(); 
          convIt != criterions_.end(); ++convIt ) {
        LOG_DBG3(itersolvestep) << "\t" << SolutionTypeEnum.ToString(convIt->first);
        convIt->second->StopSampling();
        

        // Obtain norm
        Double norm = convIt->second->GetNorm();
        normsReached &= convIt->second->Converged();
        if (nonLinLogging_) {
          std::string quantityName = SolutionTypeEnum.ToString(convIt->first);
          PtrParamNode itNode =actNode->GetByVal("quantity","name", quantityName)
                      ->Get("iteration",ParamNode::APPEND);
          itNode->Get("count")->SetValue(iter+1);
          itNode->Get("norm")->SetValue(norm);
          itNode->Get("converged")->SetValue(convIt->second->Converged());
          
          // put information also to message stream
          msg << std::setw(width[0]) << quantityName 
              << std::setw(width[1]) << (convIt->second->Converged() ? "yes" : "no")
              << std::setw(width[2]) << std::setiosflags(std::ios::scientific) << norm
              << std::setw(width[3]) << std::setiosflags(std::ios::scientific) << convIt->second->GetFinalNorm()
              << std::endl;
        }
      }
      iter++;
    } // end of while-loop
    
    LOG_DBG(itersolvestep) << "Finished loop in " << iter << " iterations";
    
    actNode->Get("numIters")->SetValue(iter);
    if (iter >= maxiter_ && !normsReached) {
      actNode->Get("converged")->SetValue(normsReached);
      if( stopOnDivergence_ ) {
        EXCEPTION("Iterative PDE coupling did not converge\n\n" << msg.str());
      } else {
        WARN("Iterative PDE coupling did not converge\n\n" << msg.str());
      }
    } else {
      actNode->Get("converged")->SetValue(normsReached); 
    }
    
    // now we are converged and can compute any postprocessing-quantities
    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++)
      rPDE_.PDEs_[i]->GetSolveStep()->PostStepStatic();
  }


  //----------------------- TRANSIENT-----------------------------------------
  
  void IterSolveStep::InitTimeStepping() {
    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {
      rPDE_.PDEs_[i]->GetSolveStep()->InitTimeStepping();
    }
  }
  
  void IterSolveStep::PreStepTrans()
  {
  }
  

  void IterSolveStep::SolveStepTrans() {

    LOG_TRACE(itersolvestep) << "--------------------------------------"; 
    LOG_TRACE(itersolvestep) <<" Solving transient step " << actStep_ 
                             << ", t = " << actTime_;
    LOG_TRACE(itersolvestep) << "--------------------------------------\n";
    
    UInt iter = 0;
    bool normsReached = false;
    std::map<SolutionType, shared_ptr<ConvCriterion> >::iterator convIt;


    PtrParamNode actNode = convNode_->Get("step",ParamNode::APPEND); 
    actNode->Get("number")->SetValue(actStep_);
    if (nonLinLogging_) {

      for( convIt = criterions_.begin(); 
          convIt != criterions_.end(); ++convIt ) {
        convIt->second->ResetValues();
        std::string quantityName = SolutionTypeEnum.ToString(convIt->first);
        actNode->Get("quantity")->Get("name")->SetValue(quantityName);
      }
    }
    // Create a stream for logging the convergence. It will be used in the end,
    // if no convergence was achieved to report the status of the lastest iteration.
    std::stringstream msg;
    UInt width[4] = {20, 10, 15, 15}; // define widths for convergence output

    while ( iter < maxiter_ &&  (! normsReached) ) {
	    	    
      LOG_DBG(itersolvestep) << "\n";
      LOG_DBG(itersolvestep) << "=== Iteration #" << iter+1 << "===";

      //std::cout << "=== Iteration #" << iter+1 << "===" << std::endl;

      // --------------------------------------
      //  1) Re-Set all convergence criterions
      // --------------------------------------
      LOG_DBG(itersolvestep) << "Calling StartSampling for criterions";
      for( convIt = criterions_.begin(); 
          convIt != criterions_.end(); ++convIt ) {
        convIt->second->StartSampling();
        LOG_DBG3(itersolvestep) << "\t" << SolutionTypeEnum.ToString(convIt->first);
      }

      // -----------------------------------
      //  2) Calculate Single PDEs
      // -----------------------------------
      for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {

        LOG_DBG(itersolvestep) << "Processing PDE '" << 
            rPDE_.PDEs_[i]->GetName() << "'";

        // Set the iteration counter and the PDE name for the exportLinSys output-stream
        Domain* domain = rPDE_.PDEs_[i]->GetDomain();
        AnalysisID& id = domain->GetDriver()->GetAnalysisId();
        id.coupleIter = iter;
        id.pdeName = rPDE_.PDEs_[i]->GetName();

        rPDE_.PDEs_[i]->GetSolveStep()->SetActTime(actTime_);
        rPDE_.PDEs_[i]->GetSolveStep()->SetActStep(actStep_);
        rPDE_.PDEs_[i]->GetSolveStep()->SetCouplingIter(iter);
        rPDE_.PDEs_[i]->GetSolveStep()->PreStepTrans();
        rPDE_.PDEs_[i]->GetSolveStep()->SolveStepTrans();
        rPDE_.PDEs_[i]->GetSolveStep()->PostStepTrans();
      } // end of for-loop

      // -----------------------------------
      //  3) Compute Coupling Criterions 
      // -----------------------------------
      normsReached = true;
      msg.str(""); // clear logging information stream
      msg << std::setw(width[0]) << "Quantity"
               << std::setw(width[1]) << "Converged"
               << std::setw(width[2]) << "Norm"
               << std::setw(width[3]) << "Goal" << std::endl
               << std::setw(width[0]+width[1]+width[2]+width[3]);
      msg << std::setfill('-') << "" << std::setfill(' ') << std::endl;
      
      LOG_DBG(itersolvestep) << "Calling StopSampling for criterions";
      for( convIt = criterions_.begin(); 
          convIt != criterions_.end(); ++convIt ) {
        LOG_DBG3(itersolvestep) << "\t" << SolutionTypeEnum.ToString(convIt->first);
            
        convIt->second->StopSampling();
        

        // Obtain norm
        Double norm = convIt->second->GetNorm();
        normsReached &= convIt->second->Converged();
        if (nonLinLogging_) {
          std::string quantityName = SolutionTypeEnum.ToString(convIt->first);
          PtrParamNode itNode =actNode->GetByVal("quantity","name", quantityName)
                  ->Get("iteration",ParamNode::APPEND);
          itNode->Get("count")->SetValue(iter+1);
          itNode->Get("norm")->SetValue(norm);
          itNode->Get("converged")->SetValue(convIt->second->Converged());

          // put information also to message stream
          msg << std::setw(width[0]) << quantityName 
              << std::setw(width[1]) << (convIt->second->Converged() ? "yes" : "no")
              << std::setw(width[2]) << std::setiosflags(std::ios::scientific) << norm
              << std::setw(width[3]) << std::setiosflags(std::ios::scientific) << convIt->second->GetFinalNorm()
              << std::endl;

          //std::cout << "Quantity " << quantityName << " :" << norm << std::endl;
        }
      }
      // reset the glmVec to the initial copy if we did not converge
      // here we only trigger the event, the reset will happen later during the function call
      // we do this so that other PDEs still see the newly calcualted values, only when the calculation of the new sub-step for the PDE is done, we reset the glmVector
      for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {
        std::map<SolutionType, shared_ptr<BaseFeFunction> > feFunctions;
        feFunctions = rPDE_.PDEs_[i]->GetFeFunctions();
        std::map<SolutionType, shared_ptr<BaseFeFunction> >::iterator fncIt;
        for(fncIt = feFunctions.begin();fncIt != feFunctions.end(); ++fncIt){
          fncIt->second->GetTimeScheme()->ResetGlmVector();
        }
      }
      iter++;
    } // end of while-loop
   
   /* 
           for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {

        LOG_DBG(itersolvestep) << "Processing PDE '" << 
            rPDE_.PDEs_[i]->GetName() << "'";

            std::cout << "Processing PDE " << rPDE_.PDEs_[i]->GetName() << std::endl;

        rPDE_.PDEs_[i]->GetSolveStep()->SetActTime(actTime_);
        rPDE_.PDEs_[i]->GetSolveStep()->SetActStep(actStep_);
        rPDE_.PDEs_[i]->GetSolveStep()->SetCouplingIter(iter);
        rPDE_.PDEs_[i]->GetSolveStep()->PreStepStatic();
       if(rPDE_.PDEs_[i]->GetName() == "mechanic"){
	 }else{
        rPDE_.PDEs_[i]->GetSolveStep()->SolveStepStatic(analysis_id);
  }
        rPDE_.PDEs_[i]->GetSolveStep()->PostStepStatic();
      } // end of for-loop
*/
    
    
    LOG_DBG(itersolvestep) << "Finished loop in " << iter << " iterations";

    actNode->Get("numIters")->SetValue(iter);
    if (iter >= maxiter_ && !normsReached) {

      actNode->Get("converged")->SetValue(normsReached);
      if( stopOnDivergence_ ) {
        EXCEPTION("Iterative PDE coupling did not converge\n\n" << msg.str());
      } else {
        WARN("Iterative PDE coupling did not converge\n\n" << msg.str());
      }
    } else {
      actNode->Get("converged")->SetValue(normsReached); 
    }
    
    // now we are converged and can compute any postprocessing-quantities
    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++)
      rPDE_.PDEs_[i]->GetSolveStep()->PostStepTrans();
  } 

  //----------------------- HARMONIC---------------------------------------
  void IterSolveStep::PreStepHarmonic(){
    }

  void IterSolveStep::SolveStepHarmonic()
  {
    LOG_TRACE(itersolvestep) << "--------------------------------------";
    LOG_TRACE(itersolvestep) <<" Solving harmonic step " << actStep_
                             << ", f = " << actFreq_;
    LOG_TRACE(itersolvestep) << "--------------------------------------\n";

    UInt iter = 0;
    bool normsReached = false;
    std::map<SolutionType, shared_ptr<ConvCriterion> >::iterator convIt;


    PtrParamNode actNode = convNode_->Get("step",ParamNode::APPEND);
    actNode->Get("number")->SetValue(actStep_);
    if (nonLinLogging_) {

      for( convIt = criterions_.begin();
          convIt != criterions_.end(); ++convIt ) {
        convIt->second->ResetValues();
        std::string quantityName = SolutionTypeEnum.ToString(convIt->first);
        actNode->Get("quantity")->Get("name")->SetValue(quantityName);
      }
    }
    // Create a stream for logging the convergence. It will be used in the end,
    // if no convergence was achieved to report the status of the lastest iteration.
    std::stringstream msg;
    UInt width[4] = {20, 10, 15, 15}; // define widths for convergence output

    while ( iter < maxiter_ &&  (! normsReached) ) {

      LOG_DBG(itersolvestep) << "\n";
      LOG_DBG(itersolvestep) << "=== Iteration #" << iter+1 << "===";

      //std::cout << "=== Iteration #" << iter+1 << "===" << std::endl;

      // --------------------------------------
      //  1) Re-Set all convergence criterions
      // --------------------------------------
      LOG_DBG(itersolvestep) << "Calling StartSampling for criterions";
      for( convIt = criterions_.begin();
          convIt != criterions_.end(); ++convIt ) {
        convIt->second->StartSampling();
        LOG_DBG3(itersolvestep) << "\t" << SolutionTypeEnum.ToString(convIt->first);
      }

      // -----------------------------------
      //  2) Calculate Single PDEs
      // -----------------------------------
      for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {

        LOG_DBG(itersolvestep) << "Processing PDE '" <<
            rPDE_.PDEs_[i]->GetName() << "'";

        rPDE_.PDEs_[i]->GetSolveStep()->SetActFreq(actFreq_);
        rPDE_.PDEs_[i]->GetSolveStep()->SetActStep(actStep_);
        rPDE_.PDEs_[i]->GetSolveStep()->SetCouplingIter(iter);
        rPDE_.PDEs_[i]->GetSolveStep()->PreStepHarmonic();
        rPDE_.PDEs_[i]->GetSolveStep()->SolveStepHarmonic();
        rPDE_.PDEs_[i]->GetSolveStep()->PostStepHarmonic();
      } // end of for-loop

      // -----------------------------------
      //  3) Compute Coupling Criterions
      // -----------------------------------
      normsReached = true;
      msg.str(""); // clear logging information stream
      msg << std::setw(width[0]) << "Quantity"
               << std::setw(width[1]) << "Converged"
               << std::setw(width[2]) << "Norm"
               << std::setw(width[3]) << "Goal" << std::endl
               << std::setw(width[0]+width[1]+width[2]+width[3]);
      msg << std::setfill('-') << "" << std::setfill(' ') << std::endl;

      LOG_DBG(itersolvestep) << "Calling StopSampling for criterions";
      for( convIt = criterions_.begin();
          convIt != criterions_.end(); ++convIt ) {
        LOG_DBG3(itersolvestep) << "\t" << SolutionTypeEnum.ToString(convIt->first);

        convIt->second->StopSampling();


        // Obtain norm
        Double norm = convIt->second->GetNorm();
        normsReached &= convIt->second->Converged();
        std::string quantityName = SolutionTypeEnum.ToString(convIt->first);
        LOG_DBG3(itersolvestep) << "Quantity " << quantityName << " :" << norm << std::endl;
        if (nonLinLogging_) {
          PtrParamNode itNode =actNode->GetByVal("quantity","name", quantityName)
                  ->Get("iteration",ParamNode::APPEND);
          itNode->Get("count")->SetValue(iter+1);
          itNode->Get("norm")->SetValue(norm);
          itNode->Get("converged")->SetValue(convIt->second->Converged());

          // put information also to message stream
          msg << std::setw(width[0]) << quantityName
              << std::setw(width[1]) << (convIt->second->Converged() ? "yes" : "no")
              << std::setw(width[2]) << std::setiosflags(std::ios::scientific) << norm
              << std::setw(width[3]) << std::setiosflags(std::ios::scientific) << convIt->second->GetFinalNorm()
              << std::endl;
        }
      }
      iter++;
    } // end of while-loop

    LOG_DBG(itersolvestep) << "Finished loop in " << iter << " iterations";

    actNode->Get("numIters")->SetValue(iter);
    if (iter >= maxiter_ && !normsReached) {

      actNode->Get("converged")->SetValue(normsReached);
      if( stopOnDivergence_ ) {
        EXCEPTION("Iterative PDE coupling did not converge\n\n" << msg.str());
      } else {
        WARN("Iterative PDE coupling did not converge\n\n" << msg.str());
      }
    } else {
      actNode->Get("converged")->SetValue(normsReached);
    }

    // now we are converged and can compute any postprocessing-quantities
    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++)
      rPDE_.PDEs_[i]->GetSolveStep()->PostStepHarmonic();
  }


  void IterSolveStep::SetActTime( const Double actTime )
  {

    actTime_ = actTime;

    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {
      actAnalysisType_ = rPDE_.PDEs_[i]->GetAnalysisType();

      if ( actAnalysisType_ == BasePDE::TRANSIENT )
        rPDE_.PDEs_[i]->GetSolveStep()->SetActTime(actTime);
    }
  }

  void IterSolveStep::SetActFreq( const Double actFreq )
  {

    actFreq_ = actFreq;
    
    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {
      actAnalysisType_ = rPDE_.PDEs_[i]->GetAnalysisType();

      if ( actAnalysisType_ == BasePDE::HARMONIC )
        rPDE_.PDEs_[i]->GetSolveStep()->SetActFreq(actFreq);
    }
  }

  void IterSolveStep::SetActStep( const UInt actStep )
  {
    
    actStep_ = actStep;

    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {
      actAnalysisType_ = rPDE_.PDEs_[i]->GetAnalysisType();

      if ( actAnalysisType_ == BasePDE::TRANSIENT )
        rPDE_.PDEs_[i]->GetSolveStep()->SetActStep(actStep);
      // Dirty Hack!!
      else if ( actAnalysisType_ == BasePDE::HARMONIC )
        rPDE_.PDEs_[i]->GetSolveStep()->SetActStep(1);
    }
  }


  void IterSolveStep::SetNumTimeSteps( UInt numTimeStep)
  {

    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {
      actAnalysisType_ = rPDE_.PDEs_[i]->GetAnalysisType();

      if ( actAnalysisType_ == BasePDE::TRANSIENT ) {
        rPDE_.PDEs_[i]->GetSolveStep()->SetNumTimeSteps(numTimeStep);
        numTimeStep_=numTimeStep;
      }
    }
  }

  void IterSolveStep::SetStartStep( const UInt startStep )
  {
    
    for (UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {
      actAnalysisType_ = rPDE_.PDEs_[i]->GetAnalysisType();

      if ( actAnalysisType_ == BasePDE::TRANSIENT )
        rPDE_.PDEs_[i]->GetSolveStep()->SetStartStep(startStep);
    }
  }

  void IterSolveStep::SetupGetRidOfZerosActive() {
    for(UInt i=0; i<rPDE_.PDEs_.GetSize(); i++) {
      rPDE_.PDEs_[i]->GetSolveStep()->SetupGetRidOfZerosActive();
    }
  }
    
} // end of namespace
