// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
//================================================================================================
/*!
 *       \file     CoefFunctionGridNodal.cc 
 *       \brief    Implementation of nodal grid inteprolation 
 *
 *       \date     11/21/2012
 *       \author   Andreas Hueppe
 */
//================================================================================================

#include <def_expl_templ_inst.hh>

#include "CoefFunctionGridNodal.hh"
#include "Driver/SolveSteps/BaseSolveStep.hh"

#include "FeBasis/FeSpace.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include <boost/lexical_cast.hpp>
#include <boost/tr1/type_traits.hpp>
#include <set>
#include <string>
#include <limits>

namespace CoupledField{

// declare class specific logging stream
DECLARE_LOG(coeffunctiongridnodal);
DEFINE_LOG(coeffunctiongridnodal, "coefFunctionGridNodal");

template<class DATA_TYPE>
  CoefFunctionGridNodal<DATA_TYPE>::CoefFunctionGridNodal(Domain* ptDomain,
                                                          PtrParamNode configNode,
                                                          shared_ptr<RegionList> regions)
                                   :CoefFunctionGrid(ptDomain, configNode, regions){
    lastStepUpdate_ = std::numeric_limits<uint>::max();
    stepNumberInterpolationA_ = std::numeric_limits<uint>::max();
    stepNumberInterpolationB_ = std::numeric_limits<uint>::max();
    initializedSpaceFactors_ = false;
    initializedRegionNodes_ = false;
    initializedUsedRegionNodesForSum_ = false;
    initializedRegionNodeCoordinates_ = false;
#ifdef USE_OPENMP
    omp_init_lock(&updateSolutionLock_);
#endif
    
    //Set sequence Step according to XML definition
    this->aSeqStep_ = configNode->Get("sequenceStep")->As<UInt>();
    
    std::string dependType = configNode->Get("dependtype")->As<std::string>();
    if (dependType == "CONST" || dependType == "SPACE") {
      this->dependType_ = CoefFunction::SPACE;
    } else {
      this->dependType_ = CoefFunction::GENERAL;
    }

    // get some boolean setting values
    this->snapToCFSStep_ = configNode->Get("snapToCFSTimeStep")->As<bool>();
    this->verboseSum_ = configNode->Get("verboseSum")->As<bool>();
    this->verboseTimeFreqFactor_ = configNode->Get("verboseTimeFreqFactor")->As<bool>();
    
    

    // get the solution quantity
    solName_ = configNode->Get("quantity")->As<std::string>();
    solType_ = SolutionTypeEnum.Parse(solName_);
    isComplex_ =  std::tr1::is_same<DATA_TYPE,Complex>::value;
    this->numEqns_ = 0;
    
    mp_ = domain_->GetMathParser();
    mHandleStep_ = mp_->GetNewHandle(true);
    std::string var;
    if(isComplex_)
      var = "f";
    else
      var = "t";
    mp_->SetExpr(mHandleStep_, var);
    this->InitGlobalFactorFunctions(configNode);
    //if(this->srcGrid_->GetDim() == 3 && this->domain_->GetGrid()->GetDim()==2){
      //TODO: Set 3D coodinate system to global factor
   //   WARN("3D->2D not supported.");
      //factorFnc_->SetCoordinateSystem();
   // }
  }

  //! Destructor
  template<class DATA_TYPE>
  CoefFunctionGridNodal<DATA_TYPE>::~CoefFunctionGridNodal(){
#ifdef USE_OPENMP
    omp_destroy_lock(&updateSolutionLock_);
#endif
  }
  
  //! function to create on string of global factors from multiple strings
  std::string OneFactorStringFromMultipleFactors(StdVector<std::string>& factors) {
    UInt size = factors.GetSize();
    if (size > 1) {
      std::stringstream ss;
      for (UInt i = 0; i < size; ++i) {
        ss << "(" << factors[i] << ")";
        if (i < size - 1) {
          ss << "*";
        }
      }
      return ss.str();
    } else if (size == 1) {
      return factors[0];
    }
    return std::string("");
  }
  
  template<class DATA_TYPE>
  shared_ptr<CoefFunction> CoefFunctionGridNodal<DATA_TYPE>::CreateFactorFunction(std::string factorString){
    shared_ptr<CoefFunction> factorFnc;
    if(isComplex_) {
      factorFnc = CoefFunction::Generate(mp_,Global::COMPLEX,factorString);
    } else {
      factorFnc = CoefFunction::Generate(mp_,Global::REAL,factorString);
    }
    return factorFnc;
  }
  
  template<class DATA_TYPE>
  CoefFunction::CoefDependType CoefFunctionGridNodal<DATA_TYPE>::ReadFactorFunctionStringDependency(std::string factorString){
    shared_ptr<CoefFunction> function = this->CreateFactorFunction(factorString);
    return function->GetDependency();
  }
  
  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::WriteGlobalFactorsToXML(PtrParamNode configNode){
    ParamNodeList subNodes = configNode->GetList("globalFactor");
    UInt size = subNodes.GetSize();
    if (size > 1) {
      for (UInt i = 0; i < subNodes.GetSize(); i++) {
        std::string factorString = subNodes[i]->As<std::string>();
        this->extDataInfo_->Get("interpolation")->Get("factors")->Get("factor",ParamNode::APPEND)->Get("value")->SetValue(factorString);
      }
    } else if (size == 1) {
      std::string factorString = subNodes[0]->As<std::string>();
      this->extDataInfo_->Get("interpolation")->Get("factor",ParamNode::APPEND)->SetValue(factorString);
    }
  }

  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::InitGlobalFactorFunctions(PtrParamNode configNode){
    // collect strings of given factors according to their type
    StdVector<std::string> constantFactors;
    StdVector<std::string> timeFreqFactors;
    StdVector<std::string> spaceFactors;
    StdVector<std::string> generalFactors;

    // test each string and detect its type
    ParamNodeList subNodes = configNode->GetList("globalFactor");
    for (UInt i = 0; i < subNodes.GetSize(); i++) {
      std::string factorString = subNodes[i]->As<std::string>();
      CoefFunction::CoefDependType type = ReadFactorFunctionStringDependency(factorString);
      if (type == CoefFunction::CONSTANT) {
        constantFactors.Push_back(factorString);
      } else if (type == CoefFunction::TIMEFREQ) {
        timeFreqFactors.Push_back(factorString);
      } else if (type == CoefFunction::SPACE) {
        spaceFactors.Push_back(factorString);
      } else if (type == CoefFunction::GENERAL) {
        generalFactors.Push_back(factorString);
      }
    }
    hasFactor_ = subNodes.GetSize() > 0;
    
    // preoptimize factors not dependend on time / frequency
    // afterwards only maximum one of these factor remains
    hasSpaceFactor_ = false;
    hasConstantFactor_ = false;
    if (spaceFactors.GetSize() > 0 && constantFactors.GetSize() > 0) {
      for (UInt i = 0; i < constantFactors.GetSize(); i++) {
        spaceFactors.Push_back(constantFactors[i]);
      }
      hasSpaceFactor_ = true;
    } else if (spaceFactors.GetSize() > 0)  {
      hasSpaceFactor_ = true;
    } else if (constantFactors.GetSize() > 0) {
      hasConstantFactor_ = true;
    }
    if (hasSpaceFactor_) {
      std::string factorString = OneFactorStringFromMultipleFactors(spaceFactors);
      spaceFactorFunction_ = CreateFactorFunction(factorString);
      this->dependType_ = GetMaxCoefDependType(this->dependType_, spaceFactorFunction_->GetDependency());
    } else if (hasConstantFactor_) {
      std::string factorString = OneFactorStringFromMultipleFactors(constantFactors);
      shared_ptr<CoefFunction> factorFnc = CreateFactorFunction(factorString);
      LocPointMapped lpm;
      factorFnc->GetScalar(constantFactor_,lpm);
    }
    
    // create factor depending solely on time / frequency
    hasTimeFreqFactor_ = timeFreqFactors.GetSize() > 0;
    if (hasTimeFreqFactor_) {
      std::string factorString = OneFactorStringFromMultipleFactors(timeFreqFactors);
      timeFreqFactorFunction_ = CreateFactorFunction(factorString);
      this->dependType_ = GetMaxCoefDependType(this->dependType_, timeFreqFactorFunction_->GetDependency());
    }
    hasGeneralFactor_ = generalFactors.GetSize() > 0;
    if (hasGeneralFactor_) {
      std::string factorString = OneFactorStringFromMultipleFactors(generalFactors);
      generalFactorFunction_ = CreateFactorFunction(factorString);
      this->dependType_ = GetMaxCoefDependType(this->dependType_, generalFactorFunction_->GetDependency());
    }
  }
  
  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::GetTensor(Matrix<DATA_TYPE>& CoefMat,
                                                   const LocPointMapped& lpm ){
    EXCEPTION("GetTensor is not implemented here"); 
  }

  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::GetVector(Vector<DATA_TYPE>& CoefMat,
                                                   const LocPointMapped& lpm ){
    EXCEPTION("GetVector is not implemented here"); 
  }

  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::GetScalar(DATA_TYPE& CoefMat,
                                                   const LocPointMapped& lpm ){
    EXCEPTION("GetScalar is not implemented here"); 
  }


  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::CreateOperator(UInt spaceDim, UInt dofDim){
    if(dofDim == 1){
      if(spaceDim == 2)
        this->myOperator_.reset(new IdentityOperator<FeH1,2,1,DATA_TYPE>());
      else if(spaceDim == 3)
        this->myOperator_.reset(new IdentityOperator<FeH1,3,1,DATA_TYPE>());
    }else if(dofDim == 2){
      if(spaceDim == 2)
        this->myOperator_.reset(new IdentityOperator<FeH1,2,2,DATA_TYPE>());
      else if(spaceDim == 3)
        this->myOperator_.reset(new IdentityOperator<FeH1,3,2,DATA_TYPE>());
    }else if(dofDim == 3){
      if(spaceDim == 2)
        this->myOperator_.reset(new IdentityOperator<FeH1,2,3,DATA_TYPE>());
      else if(spaceDim == 3)
        this->myOperator_.reset(new IdentityOperator<FeH1,3,3,DATA_TYPE>());
    }
  }

  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::CreateDivOperator(UInt spaceDim, UInt dofDim){
    if(spaceDim != dofDim)
      EXCEPTION("CoefFunctionGridNodal<DATA_TYPE>: Divergence need vectorial data!");

    if(spaceDim == 2)
      this->myOperator_.reset(new ScalarDivergenceOperator<FeH1,2,DATA_TYPE>());
    else if(spaceDim == 3)
      this->myOperator_.reset(new ScalarDivergenceOperator<FeH1,3,DATA_TYPE>());
  }

  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::InitSolVec(){
    numNodes_ = srcGrid_->GetNumNodes();
    numEqns_ = (numNodes_ + 1) * dimDof_;
    solVec_.Resize(numEqns_);
    solVec_.Init(0.0);
  }

  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::GetElemSolution(Vector<DATA_TYPE> & sol, UInt eNum){
    StdVector<UInt> connect;
    srcGrid_->GetElemNodes(connect,eNum);
    sol.Resize(connect.GetSize()*dimDof_);
    if (dimDof_ > 1) {
      for(UInt aNode=0;aNode < connect.GetSize();aNode++){
        UInt node = connect[aNode];
        for(UInt d = 0; d<dimDof_;++d){
          sol[aNode*dimDof_+d] = solVec_[node*dimDof_+d];
        }
      }
    } else {
      for(UInt aNode=0;aNode < connect.GetSize();aNode++){
        UInt node = connect[aNode];
        sol[aNode] = solVec_[node];
      }
    }
  }

  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::InitRegionNodes(){
    if (initializedRegionNodes_) {
      return;
    }
    regionNodes_.Resize(srcRegions_.size());
    std::set<std::string>::iterator regIter = srcRegions_.begin();
    for( UInt i = 0; regIter != srcRegions_.end(); ++i,++regIter) {
      srcGrid_->GetNodesByRegion(regionNodes_[i], srcGrid_->GetRegion().Parse(*regIter));
    }
    initializedRegionNodes_ = true;
  }
  
  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::InitRegionNodeCoordinates(){
    if (initializedRegionNodeCoordinates_) {
      return;
    }
    this->InitRegionNodes();
    regionNodeCoordinates_.Resize(srcRegions_.size());
    std::set<std::string>::iterator regIter = srcRegions_.begin();
    for( UInt i = 0; regIter != srcRegions_.end(); ++i,++regIter) {
      StdVector<UInt>& nodeNums = regionNodes_[i];
      UInt size = nodeNums.GetSize();
      StdVector<Vector<Double> >& CoordVec = regionNodeCoordinates_[i];
      CoordVec.Resize(size);
#pragma omp parallel for      
      for( UInt aN=0;aN<size;aN++){
        srcGrid_->GetNodeCoordinate(CoordVec[aN],nodeNums[aN],true);
      }
    }
    initializedRegionNodeCoordinates_ = true;
  }
  
  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::ClearRegionNodeCoordinates(){
    initializedRegionNodeCoordinates_ = false;
    regionNodeCoordinates_ = StdVector<StdVector<Vector<Double> > >();
  }
  
  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::InitSpaceFactor(){
    if (initializedSpaceFactors_) {
      return;
    }
    initializedSpaceFactors_ = true;
    if (!hasSpaceFactor_) {
      return;
    }
    this->InitRegionNodes();
    this->InitRegionNodeCoordinates();
    spaceFactor_.Resize(srcRegions_.size());
    std::set<std::string>::iterator regIter = srcRegions_.begin();
    for( UInt i = 0; regIter != srcRegions_.end(); ++i,++regIter) {
      StdVector<UInt>& nodeNums = regionNodes_[i];
      StdVector<Vector<Double> >& CoordVec = regionNodeCoordinates_[i];
      
      StdVector< DATA_TYPE >& factors = spaceFactor_[i];
      factors.Resize(nodeNums.GetSize());
      spaceFactorFunction_->GetScalarValuesAtCoords(CoordVec,factors,this->domain_->GetGrid());
    }
    this->ClearRegionNodeCoordinates();
  }
  
  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::InitUsedRegionNodesForSum(){
    if (this->initializedUsedRegionNodesForSum_) {
      return;
    }
    this->InitRegionNodes();
    this->usedRegionNodesForSum_.Clear();
    std::vector<bool> allUsedNode(numNodes_ + 1, false);
    std::set<std::string>::iterator regIter = srcRegions_.begin();
    for( UInt i = 0; regIter != srcRegions_.end(); ++i,++regIter) {
      StdVector<UInt>& nodeNums = regionNodes_[i];
      const UInt size = nodeNums.GetSize();
      std::vector<bool> useNode(size,false);
      for (UInt iN = 0; iN < size; ++iN) {
        UInt node = nodeNums[iN];
        if (!allUsedNode[node]) {
          allUsedNode[node] = true;
          useNode[iN] = true;
        }
      }
      this->usedRegionNodesForSum_.Push_back(useNode);
    }
    this->useAllRegionNodesForSum_ = true;
    regIter = srcRegions_.begin();
    for( UInt i = 0; regIter != srcRegions_.end() && this->useAllRegionNodesForSum_; ++i,++regIter) {
      std::vector<bool>& useNode = this->usedRegionNodesForSum_[i];
      const UInt size = useNode.size();
      for (UInt n = 0; n < size && this->useAllRegionNodesForSum_; ++n) {
        this->useAllRegionNodesForSum_ = useNode[n];
      }
    }
    this->initializedUsedRegionNodesForSum_ = true;
  }
  
  template<class DATA_TYPE>
  template<bool useSpaceFactorA, bool useSpaceFactorB, bool useConstFactor, 
              bool countSum, bool countAllValues, UInt dimDof>
  void CoefFunctionGridNodal<DATA_TYPE>::CopyResult(Vector<DATA_TYPE>& sol,
                          Vector<DATA_TYPE>& res,
                          StdVector<UInt>& nodeNums, 
                          StdVector<DATA_TYPE>& spaceFactorA, 
                          StdVector<DATA_TYPE>& spaceFactorB,
                          const DATA_TYPE constFactor,
                          StdVector<DATA_TYPE>& sum, 
                          StdVector<DATA_TYPE>& factorSum,
                          std::vector<bool>& countNodes) {
    if (dimDof < 1) {
      return;
    }
    const UInt size = nodeNums.GetSize();
    const bool useSpaceFactor = useSpaceFactorA || useSpaceFactorB;
    const bool useOnlyConstFactor = useConstFactor && !useSpaceFactor;
#pragma omp parallel
    {
      StdVector<DATA_TYPE> tSum(dimDof);
      StdVector<DATA_TYPE> tFactorSum(dimDof);
      if (countSum) {
        tSum.Init(0.0);
        tFactorSum.Init(0.0);
      }
      DATA_TYPE tSum1 = 0.0;
      DATA_TYPE tFactorSum1 = 0.0;
      DATA_TYPE value;
      DATA_TYPE factor = 1.0;
      UInt node;
      UInt sourceIndex;
      UInt targetIndex;
#pragma omp for
      for (UInt i = 0; i < size; ++i) {
        if (useSpaceFactorA) {
          if (useSpaceFactorB) {
            if (useConstFactor) {
              factor = spaceFactorA[i] * spaceFactorB[i] * constFactor;
            } else {
              factor = spaceFactorA[i] * spaceFactorB[i];
            }
          } else {
            if (useConstFactor) {
              factor = spaceFactorA[i] * constFactor;
            } else {
              factor = spaceFactorA[i];
            }
          }
        } else if (useSpaceFactorB) {
          if (useConstFactor) {
            factor = spaceFactorB[i] * constFactor;
          } else {
            factor = spaceFactorB[i];
          }
        }
        node = nodeNums[i];
        if (dimDof > 1) {
          sourceIndex = i * dimDof;
          targetIndex = node * dimDof;
          for (UInt d = 0; d < dimDof; ++d) {
            value = res[sourceIndex + d];
            if (countSum) {
              if (countAllValues || countNodes[i]) {
                tSum[d] += value;
              }
            }
            if (useOnlyConstFactor) {
              value *= constFactor;
            } else if (useSpaceFactor) {
              value *= factor;
            }
            sol[targetIndex + d] = value;
            if (countSum && useSpaceFactor) {
              if (countAllValues || countNodes[i]) {
                tFactorSum[d] += value;
              }
            }
          }
        } else { 
          // here we define the special case of dimDof == 1, because we do not exactly know 
          // the optimization behaviour of the compiler
          value = res[i];
          if (countSum) {
            if (countAllValues || countNodes[i]) {
              tSum1 += value;
            }
          }
          if (useOnlyConstFactor) {
            value *= constFactor;
          } else if (useSpaceFactor) {
            value *= factor;
          }
          sol[node] = value;
          if (countSum && useSpaceFactor) {
            if (countAllValues || countNodes[i]) {
              tFactorSum1 += value;
            }
          }
        }
      }
#pragma omp critical
      if (countSum && dimDof == 1) {
        tSum[0] += tSum1;
        tFactorSum[0] += tFactorSum1;
      }
      {
        if (countSum) {
          for (UInt d = 0; d < dimDof; ++d) {
            sum[d] += tSum[d];
            if (useSpaceFactor) {
              factorSum[d] += tFactorSum[d];
            } else {
              if (useOnlyConstFactor) {
                factorSum[d] += tSum[d] * constFactor;
              } else {
                factorSum[d] += tSum[d];
              }
            }
          }
        }
      }
    }
  }
  
  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::ReadSolution(UInt step, Double stepValue, Vector<DATA_TYPE> & sol){
    Double actValue = this->mp_->Eval(mHandleStep_);
    std::string varName = this->isComplex_ ? "f" : "t";
    this->mp_->SetValue(MathParser::GLOB_HANDLER, varName, stepValue);
    this->ReadSolution(step, sol);
    this->mp_->SetValue(MathParser::GLOB_HANDLER, varName, actValue);
  }
  
  template<class DATA_TYPE>
  void CoefFunctionGridNodal<DATA_TYPE>::ReadSolution(UInt step,Vector<DATA_TYPE> & sol){
    // evaulate space factors and time/frequency factors
    this->InitSpaceFactor();
    bool useConstantFactor = this->hasConstantFactor_;
    DATA_TYPE constantFactor = this->constantFactor_;
    if (this->hasTimeFreqFactor_) {
      DATA_TYPE timeFreqFactor;
      LocPointMapped lpm;
      this->timeFreqFactorFunction_->GetScalar(timeFreqFactor,lpm);
      if (this->verboseTimeFreqFactor_) {
        std::cout << "Global time/frequency factor at " << this->allSrcRegionNames_ << ": " << timeFreqFactor << std::endl;
      }
      if (useConstantFactor) {
        constantFactor *= timeFreqFactor;
      } else {
        useConstantFactor = true;
        constantFactor = timeFreqFactor;
      }
    }
    // data for summing loaded values
    StdVector<DATA_TYPE> loadedSum(dimDof_);
    loadedSum.Init(0.0);
    StdVector<DATA_TYPE> factorSum(dimDof_);
    factorSum.Init(0.0);
    
    std::set<std::string>::iterator regIter = srcRegions_.begin();
    sol.Resize(numEqns_);
    sol.Init(0.0);
    this->InitRegionNodes();
    for( UInt i = 0; regIter != srcRegions_.end(); ++i,++regIter) {
      // get result and cast it to required type
      shared_ptr<BaseResult> Bres = domain_->GetResultHandler()->GetResult( this->inputId_, this->aSeqStep_ , step , this->solType_, *regIter );
      Vector<DATA_TYPE> resVec;
      try{
        Result<DATA_TYPE>* myResult = dynamic_cast<Result<DATA_TYPE>* >(Bres.get());
        resVec =   myResult->GetVector();
      }catch(...){
        EXCEPTION("Cannot cast to desired vector type. Are you trying to load real data into a harmonic computation?");
      }
      // get node numbers
      StdVector<UInt>& nodeNums = regionNodes_[i];
      
      const UInt size = nodeNums.GetSize();
      // evaulate general factors
      StdVector< DATA_TYPE > generalFactor;
      if (hasGeneralFactor_) {
        generalFactor.Resize(size);
        this->InitRegionNodeCoordinates();
        StdVector<Vector<Double> >& CoordVec = regionNodeCoordinates_[i];
        generalFactorFunction_->GetScalarValuesAtCoords(CoordVec,generalFactor,this->domain_->GetGrid());
      }
      StdVector< DATA_TYPE >& spaceFactor = this->hasSpaceFactor_ ? this->spaceFactor_[i] : generalFactor;
      
      if (verboseSum_) {
        this->InitUsedRegionNodesForSum();
      }
      std::vector<bool> dummy;
      std::vector<bool>& countNodes = this->verboseSum_ ? this->usedRegionNodesForSum_[i] : dummy;
      
      // get function for copying the result vector
      typename CoefFunctionGridNodal<DATA_TYPE>::CopyResultFunction::Ptr crf = GetCopyResultFunction(
        this->hasGeneralFactor_, this->hasSpaceFactor_, useConstantFactor, 
        this->verboseSum_, this->useAllRegionNodesForSum_, this->dimDof_);
      
      // copying the result vector
      crf(sol,resVec,nodeNums,generalFactor,spaceFactor,constantFactor,loadedSum,factorSum,countNodes);
    }
    if (this->verboseSum_) {
      std::cout << "Sum of " << this-> solName_ <<" read at " << allSrcRegionNames_ << ": ";
      for (UInt d = 0; d < this->dimDof_; ++d) {
        std::cout << loadedSum[d];
        if (d < this->dimDof_ - 1) {
          std::cout << ",";
        }
      }
      std::cout << std::endl;
      if (hasFactor_) {
        std::cout << "Sum of " << this-> solName_ <<" read at " << allSrcRegionNames_ << " after applying global factors: ";
        for (UInt d = 0; d < this->dimDof_; ++d) {
          std::cout << factorSum[d];
          if (d < this->dimDof_ - 1) {
            std::cout << ",";
          }
        }
        std::cout << std::endl;
      }
    }
  }

  template<class DATA_TYPE>
  bool CoefFunctionGridNodal<DATA_TYPE>::UpdateSolution(){
    const UInt stepnumber = this->domain_->GetBasePDE()->GetSolveStep()->GetActStep();
#ifdef USE_OPENMP
    omp_set_lock(&updateSolutionLock_);
#endif
    if (stepnumber == this->lastStepUpdate_) {
#ifdef USE_OPENMP
      omp_unset_lock(&updateSolutionLock_);
#endif
      return false;
    }
    if (this->dependType_ == CoefFunction::SPACE) {
      if (this->lastStepUpdate_ != std::numeric_limits<uint>::max()) {
        this->lastStepUpdate_ = stepnumber;
#ifdef USE_OPENMP
        omp_unset_lock(&updateSolutionLock_);
#endif
        return false;
      }
    }
    if (this->snapToCFSStep_) {
      this->ReadSolution(stepnumber,this->solVec_);
    } else {
      Double stepValue;
      if (this->isComplex_) {
        stepValue = this->domain_->GetBasePDE()->GetSolveStep()->GetActFreq();
      } else {
        stepValue = this->domain_->GetBasePDE()->GetSolveStep()->GetActTime();
      }
      std::map<UInt,Double>::iterator stepIter = stepValueMap_.begin();
      UInt preStepNumber = stepIter->first;
      Double preStepValue = stepIter->second;
      UInt postStepNumber = stepIter->first;
      Double postStepValue = stepIter->second;
      bool notEnd = true;
      if (stepValue > postStepValue) {
        while (stepIter->second < stepValue && notEnd) {
          preStepNumber = postStepNumber;
          preStepValue = postStepValue;
          ++stepIter;
          if (stepIter != stepValueMap_.end()) {
            postStepNumber = stepIter->first;
            postStepValue = stepIter->second;
          } else {
            notEnd = false;
          }
        }
        if (stepValue > postStepValue) {
          preStepNumber = postStepNumber;
          preStepValue = postStepValue;
        }
      }
      if (preStepNumber == postStepNumber) {
        this->ReadSolution(preStepNumber,this->solVec_);
      } else {
        const Double factorPre = (postStepValue - stepValue) / (postStepValue - preStepValue);
        const Double factorPost = (stepValue - preStepValue) / (postStepValue - preStepValue);
        
        bool Apre = preStepNumber != this->stepNumberInterpolationB_;
        const Double facA = Apre ? factorPre : factorPost;
        const Double facB = Apre ? factorPost : factorPre;
        const UInt newA = Apre ? preStepNumber : postStepNumber;
        const Double newAValue = Apre ? preStepValue : postStepValue;
        const UInt newB = Apre ? postStepNumber : preStepNumber;
        const Double newBValue = Apre ? postStepValue : preStepValue;
        
        if (this->stepNumberInterpolationA_ != newA) {
          if (this->stepNumberInterpolationA_ == std::numeric_limits<uint>::max()) {
            this->solVecInterpolationA_.Resize(this->numEqns_);
            this->solVecInterpolationA_.Init(0.0);
          }
          this->ReadSolution(newA,newAValue,this->solVecInterpolationA_);
          this->stepNumberInterpolationA_ = newA;
        }
        if (this->stepNumberInterpolationB_ != newB) {
          if (this->stepNumberInterpolationB_ == std::numeric_limits<uint>::max()) {
            this->solVecInterpolationB_.Resize(this->numEqns_);
            this->solVecInterpolationB_.Init(0.0);
          }
          this->ReadSolution(newB,newBValue,this->solVecInterpolationB_);
          this->stepNumberInterpolationB_ = newB;
        }
        
#pragma omp parallel for
        for(UInt i=0;i<numEqns_;i++){
          this->solVec_[i] = facA * this->solVecInterpolationA_[i] + facB *  this->solVecInterpolationB_[i];
        }
      }
    }
    this->lastStepUpdate_ = stepnumber;
#ifdef USE_OPENMP
    omp_unset_lock(&updateSolutionLock_);
#endif
    return true;
  }

#ifdef EXPLICIT_TEMPLATE_INSTANTIATION
  template class CoefFunctionGridNodal<Double>;
  template class CoefFunctionGridNodal<Complex>;
#endif

} // namespace CoupledField
