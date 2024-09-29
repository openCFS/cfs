// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
//================================================================================================
/*!
 *       \file     CoefFunctionGridElem.cc 
 *       \brief    Implementation of element grid inteprolation 
 *
 *       \date     29/09/2024
 *       \author   Dominik Mayrhofer
 */
//================================================================================================


#include "Utils/mathParser/mathParser.hh"
#include "CoefFunctionGridElem.hh"
#include "Driver/SolveSteps/BaseSolveStep.hh"
#include "FeBasis/FeSpace.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include <boost/lexical_cast.hpp>
#include <set>
#include <string>
#include <limits>

namespace CoupledField{

// declare class specific logging stream
DEFINE_LOG(coeffunctiongridelem, "CoefFunctionGridElem")

template<class DATA_TYPE>
  CoefFunctionGridElem<DATA_TYPE>::CoefFunctionGridElem(Domain* ptDomain,
                                                          PtrParamNode configNode,
                                                          shared_ptr<RegionList> regions)
                                   :CoefFunctionGrid(ptDomain, configNode, regions){
    lastStepUpdate_ = std::numeric_limits<unsigned int>::max();
    stepNumberInterpolationA_ = std::numeric_limits<unsigned int>::max();
    stepNumberInterpolationB_ = std::numeric_limits<unsigned int>::max();
    initializedSpaceFactors_ = false;
    initializedRegionElements_ = false;
    initializedRegionElementCoordinates_ = false;
    initializedConstantInput_ = false;
    numElements_ = 0;
    hasSpaceFactor_ = false;
    hasConstantFactor_ = false;
    hasTimeFreqFactor_ = false;
    useAllRegionElementsForSum_ = false;
    hasFactor_ = false;
    hasGeneralFactor_ = false;
    eqnMapComplete_ = false;
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
    CoefDependType inputDependType = this->dependType_;

    // get some boolean setting values
    this->snapToCFSStep_ = configNode->Get("snapToCFSTimeStep")->As<bool>();
    this->verboseSum_ = configNode->Get("verboseSum")->As<bool>();
    this->verboseTimeFreqFactor_ = configNode->Get("verboseTimeFreqFactor")->As<bool>();
    
    // get the solution quantity
    solName_ = configNode->Get("quantity")->As<std::string>();
    solType_ = SolutionTypeEnum.Parse(solName_);
    isComplex_ =  std::is_same<DATA_TYPE,Complex>::value;
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
    
    this->hasSpaceInputButTimeFactor_ = inputDependType == CoefFunction::SPACE && this->dependType_ == CoefFunction::GENERAL;
    
    //if(this->srcGrid_->GetDim() == 3 && this->domain_->GetGrid()->GetDim()==2){
      //TODO: Set 3D coodinate system to global factor
   //   WARN("3D->2D not supported.");
      //factorFnc_->SetCoordinateSystem();
   // }
  }

  //! Destructor
  template<class DATA_TYPE>
  CoefFunctionGridElem<DATA_TYPE>::~CoefFunctionGridElem(){
#ifdef USE_OPENMP
    omp_destroy_lock(&updateSolutionLock_);
#endif
  }
  
  template<class DATA_TYPE>
  std::string CoefFunctionGridElem<DATA_TYPE>::OneFactorStringFromMultipleFactors(StdVector<std::string>& factors) {
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
  shared_ptr<CoefFunction> CoefFunctionGridElem<DATA_TYPE>::CreateFactorFunction(std::string factorString){
    shared_ptr<CoefFunction> factorFnc;
    if(isComplex_) {
      factorFnc = CoefFunction::Generate(mp_,Global::COMPLEX,factorString);
    } else {
      factorFnc = CoefFunction::Generate(mp_,Global::REAL,factorString);
    }
    return factorFnc;
  }
  
  template<class DATA_TYPE>
  CoefFunction::CoefDependType CoefFunctionGridElem<DATA_TYPE>::ReadFactorFunctionStringDependency(std::string factorString){
    shared_ptr<CoefFunction> function = this->CreateFactorFunction(factorString);
    return function->GetDependency();
  }
  
  template<class DATA_TYPE>
  void CoefFunctionGridElem<DATA_TYPE>::WriteGlobalFactorsToXML(PtrParamNode configNode){
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
  void CoefFunctionGridElem<DATA_TYPE>::InitGlobalFactorFunctions(PtrParamNode configNode){
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
  void CoefFunctionGridElem<DATA_TYPE>::GetTensor(Matrix<DATA_TYPE>& CoefMat,
                                                   const LocPointMapped& lpm ){
    EXCEPTION("CoefFunctionGridElem: GetTensor is not implemented here"); 
  }

  template<class DATA_TYPE>
  void CoefFunctionGridElem<DATA_TYPE>::GetVector(Vector<DATA_TYPE>& CoefMat,
                                                   const LocPointMapped& lpm ){
    EXCEPTION("CoefFunctionGridElem: GetVector is not implemented here"); 
  }

  template<class DATA_TYPE>
  void CoefFunctionGridElem<DATA_TYPE>::GetScalar(DATA_TYPE& CoefMat,
                                                   const LocPointMapped& lpm ){
    EXCEPTION("CoefFunctionGridElem: GetScalar is not implemented here"); 
  }

  template<class DATA_TYPE>
  void CoefFunctionGridElem<DATA_TYPE>::InitSolVec(){
    numElements_ = srcGrid_->GetNumElems();
    solVec_.Resize(numElements_*dimDof_);
    solVec_.Init(0.0);
  }

  template<class DATA_TYPE>
  void CoefFunctionGridElem<DATA_TYPE>::GetElemSolution(Vector<DATA_TYPE> & sol, UInt eNum){
    sol.Resize(dimDof_);
    if (dimDof_ > 1) {
      for(UInt d = 0; d<dimDof_;++d){
        sol[d] = solVec_[eNum*dimDof_+d];
      }
    } else {
      sol[0] = solVec_[eNum];
    }
  }

  template<class DATA_TYPE>
  void CoefFunctionGridElem<DATA_TYPE>::InitRegionElements(){
    if (initializedRegionElements_) {
      return;
    }
    regionElements_.Resize(srcRegions_.size());
    std::set<std::string>::iterator regIter = srcRegions_.begin();
    for( UInt i = 0; regIter != srcRegions_.end(); ++i,++regIter) {
      srcGrid_->GetElementsByRegion(regionElements_[i], srcGrid_->GetRegion().Parse(*regIter));
    }
    initializedRegionElements_ = true;
  }
  
  template<class DATA_TYPE>
  void CoefFunctionGridElem<DATA_TYPE>::InitRegionElementCoordinates(){
    if (initializedRegionElementCoordinates_) {
      return;
    }
    this->InitRegionElements();
    regionElementCoordinates_.Resize(srcRegions_.size());
    std::set<std::string>::iterator regIter = srcRegions_.begin();
    for( UInt i = 0; regIter != srcRegions_.end(); ++i,++regIter) {
      StdVector<UInt>& elementNums = regionElements_[i];
      UInt size = elementNums.GetSize();
      StdVector<Vector<Double> >& CoordVec = regionElementCoordinates_[i];
      CoordVec.Resize(size);
#pragma omp parallel for      
      for( Integer aN=0;aN< (Integer) size;aN++){
        srcGrid_->GetElemCentroid(CoordVec[aN],elementNums[aN],true);
      }
    }
    initializedRegionElementCoordinates_ = true;
  }
  
  template<class DATA_TYPE>
  void CoefFunctionGridElem<DATA_TYPE>::ClearRegionElementCoordinates(){
    initializedRegionElementCoordinates_ = false;
    regionElementCoordinates_ = StdVector<StdVector<Vector<Double> > >();
  }
  
  template<class DATA_TYPE>
  void CoefFunctionGridElem<DATA_TYPE>::InitSpaceFactor(){
    if (initializedSpaceFactors_) {
      return;
    }
    initializedSpaceFactors_ = true;
    if (!hasSpaceFactor_) {
      return;
    }
    this->InitRegionElements();
    this->InitRegionElementCoordinates();
    spaceFactor_.Resize(srcRegions_.size());
    std::set<std::string>::iterator regIter = srcRegions_.begin();
    for( UInt i = 0; regIter != srcRegions_.end(); ++i,++regIter) {
      StdVector<UInt>& elementNums = regionElements_[i];
      StdVector<Vector<Double> >& CoordVec = regionElementCoordinates_[i];
      
      StdVector< DATA_TYPE >& factors = spaceFactor_[i];
      factors.Resize(elementNums.GetSize());
      spaceFactorFunction_->GetScalarValuesAtCoords(CoordVec,factors,this->domain_->GetGrid());
    }
    this->ClearRegionElementCoordinates();
  }

  template<class DATA_TYPE>
  void CoefFunctionGridElem<DATA_TYPE>::InitUsedRegionElementsForSum(){
    if (this->initializedUsedRegionElementsForSum_) {
      return;
    }
    this->InitRegionElements();
    this->usedRegionElementsForSum_.Clear();
    std::vector<bool> allUsedElement(numElements_ + 1, false);
    std::set<std::string>::iterator regIter = srcRegions_.begin();
    for( UInt i = 0; regIter != srcRegions_.end(); ++i,++regIter) {
      StdVector<UInt>& elementNums = regionElements_[i];
      const UInt size = elementNums.GetSize();
      std::vector<bool> useElement(size,false);
      for (UInt iN = 0; iN < size; ++iN) {
        UInt elem = elementNums[iN];
        if (!allUsedElement[elem]) {
          allUsedElement[elem] = true;
          useElement[iN] = true;
        }
      }
      this->usedRegionElementsForSum_.Push_back(useElement);
    }
    this->useAllRegionElementsForSum_ = true;
    regIter = srcRegions_.begin();
    for( UInt i = 0; regIter != srcRegions_.end() && this->useAllRegionElementsForSum_; ++i,++regIter) {
      std::vector<bool>& useElement = this->usedRegionElementsForSum_[i];
      const UInt size = useElement.size();
      for (UInt n = 0; n < size && this->useAllRegionElementsForSum_; ++n) {
        this->useAllRegionElementsForSum_ = useElement[n];
      }
    }
    this->initializedUsedRegionElementsForSum_ = true;
  }
  
  template<class DATA_TYPE>
  template<bool useSpaceFactorA, bool useSpaceFactorB, bool useConstFactor, 
              bool countSum, bool countAllValues, UInt dimDof>
  void CoefFunctionGridElem<DATA_TYPE>::CopyResult(Vector<DATA_TYPE>& sol,
                          Vector<DATA_TYPE>& res,
                          StdVector<UInt>& elementNums, 
                          StdVector<DATA_TYPE>& spaceFactorA, 
                          StdVector<DATA_TYPE>& spaceFactorB,
                          const DATA_TYPE constFactor,
                          StdVector<DATA_TYPE>& sum, 
                          StdVector<DATA_TYPE>& factorSum,
                          std::vector<bool>& countElements) {
    if (dimDof < 1) {
      return;
    }
    const UInt size = elementNums.GetSize();
    const bool useSpaceFactor = useSpaceFactorA || useSpaceFactorB;
    const bool useOnlyConstFactor = useConstFactor && !useSpaceFactor;
#pragma omp parallel
    {
      StdVector<DATA_TYPE> tSum;
      StdVector<DATA_TYPE> tFactorSum;
      if (countSum) {
        tSum.Resize(dimDof);
        tSum.Init(0.0);
        tFactorSum.Resize(dimDof);
        tFactorSum.Init(0.0);
      }
      DATA_TYPE tSum1 = 0.0;
      DATA_TYPE tFactorSum1 = 0.0;
      DATA_TYPE value;
      DATA_TYPE factor = 1.0;
      UInt element;
      UInt sourceIndex;
      UInt targetIndex;
#pragma omp for
      for (Integer i = 0; i < (Integer) size; ++i) {
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
        element = elementNums[i];
        if (dimDof > 1) {
          sourceIndex = i * dimDof;
          targetIndex = element * dimDof;
          for (UInt d = 0; d < dimDof; ++d) {
            value = res[sourceIndex + d];
            if (countSum) {
              if (countAllValues || countElements[i]) {
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
              if (countAllValues || countElements[i]) {
                tFactorSum[d] += value;
              }
            }
          }
        } else { 
          // here we define the special case of dimDof == 1, because we do not exactly know 
          // the optimization behaviour of the compiler
          value = res[i];
          if (countSum) {
            if (countAllValues || countElements[i]) {
              tSum1 += value;
            }
          }
          if (useOnlyConstFactor) {
            value *= constFactor;
          } else if (useSpaceFactor) {
            value *= factor;
          }
          sol[element] = value;
          if (countSum && useSpaceFactor) {
            if (countAllValues || countElements[i]) {
              tFactorSum1 += value;
            }
          }
        }
      }
      if (countSum && dimDof == 1) {
        tSum[0] += tSum1;
        tFactorSum[0] += tFactorSum1;
      }
#pragma omp critical
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
  void CoefFunctionGridElem<DATA_TYPE>::ReadSolution(UInt step, Double stepValue, Vector<DATA_TYPE> & sol){
    Double actValue = this->mp_->Eval(mHandleStep_);
    std::string varName = this->isComplex_ ? "f" : "t";
    this->mp_->SetValue(MathParser::GLOB_HANDLER, varName, stepValue);
    this->ReadSolution(step, sol);
    this->mp_->SetValue(MathParser::GLOB_HANDLER, varName, actValue);
  }
  
  template<class DATA_TYPE>
  void CoefFunctionGridElem<DATA_TYPE>::ReadSolution(UInt step,Vector<DATA_TYPE> & sol){
    if (srcRegions_.size() < 1) {
      return;
    }
    // evaulate space factors 
    this->InitSpaceFactor();
    
    if (this->hasSpaceInputButTimeFactor_ && !this->initializedConstantInput_ && (constantInput_.GetSize() == 0)) {
      bool cachedHasTimeFreqFactor = this->hasTimeFreqFactor_;
      this->hasTimeFreqFactor_ = false;
      bool cachedHasGeneralFactor = this->hasGeneralFactor_;
      this->hasGeneralFactor_ = false;
      
      constantInput_.Resize(srcRegions_.size());
      this->ReadSolution(stepValueMap_.begin()->first, sol);
      
      this->hasConstantFactor_ = false;
      this->hasTimeFreqFactor_ = cachedHasTimeFreqFactor;
      this->hasSpaceFactor_ = false;
      this->hasGeneralFactor_ = cachedHasGeneralFactor;
      
      this->initializedConstantInput_ = true;
    }
    
    // evaluate time/frequency factors
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
    StdVector<DATA_TYPE> loadedSum;
    StdVector<DATA_TYPE> factorSum;
    if (this->verboseSum_) {
      loadedSum.Resize(dimDof_);
      loadedSum.Init(0.0);
      factorSum.Resize(dimDof_);
      factorSum.Init(0.0);
    }
    
    // get function for copying the result vector
    typename CoefFunctionGridElem<DATA_TYPE>::CopyResultFunction::Ptr crf = GetCopyResultFunction(
      this->hasGeneralFactor_, this->hasSpaceFactor_, useConstantFactor, 
      this->verboseSum_, this->useAllRegionElementsForSum_, this->dimDof_);
      
    std::set<std::string>::iterator regIter = srcRegions_.begin();
    if (sol.GetSize() != numElements_*dimDof_) {
      sol.Resize(numElements_*dimDof_);
      sol.Init(0.0);
    }
    this->InitRegionElements();
    bool loadValues = true;
    if (this->hasSpaceInputButTimeFactor_ && this->initializedConstantInput_) {
      loadValues = false;
    }
    for( UInt i = 0; regIter != srcRegions_.end(); ++i,++regIter) {
      // get result and cast it to required type
      shared_ptr<BaseResult> Bres;
      Result<DATA_TYPE> dummyResult; //to eliminate compiler warning 
      Result<DATA_TYPE>* myResult = &dummyResult;
      if (loadValues) {
        Bres = domain_->GetResultHandler()->GetResult( this->inputId_, this->aSeqStep_ , step , this->solType_, *regIter );
        try{
          myResult = dynamic_cast<Result<DATA_TYPE>* >(Bres.get());
        }catch(...){
          EXCEPTION("Cannot cast to desired vector type. Are you trying to load real data into a harmonic computation?");
        }
      }
      Vector<DATA_TYPE>& resVec = loadValues ?  myResult->GetVector() : this->constantInput_[i];
      // get element numbers
      StdVector<UInt>& elementNums = regionElements_[i];
      
      const UInt size = elementNums.GetSize();
      // evaulate general factors
      StdVector< DATA_TYPE > generalFactor;
      if (hasGeneralFactor_) {
        generalFactor.Resize(size);
        this->InitRegionElementCoordinates();
        StdVector<Vector<Double> >& CoordVec = regionElementCoordinates_[i];
        generalFactorFunction_->GetScalarValuesAtCoords(CoordVec,generalFactor,this->domain_->GetGrid());
      }
      StdVector< DATA_TYPE >& spaceFactor = this->hasSpaceFactor_ ? this->spaceFactor_[i] : generalFactor;
      
      if (verboseSum_) {
        this->InitUsedRegionElementsForSum();
      }
      std::vector<bool> dummy;
      std::vector<bool>& countElements = this->verboseSum_ ? this->usedRegionElementsForSum_[i] : dummy;
      
      // copying the result vector
      if (this->hasSpaceInputButTimeFactor_ && !this->initializedConstantInput_) {
        this->constantInput_[i].Resize(resVec.GetSize());
        StdVector<UInt> dummyElementNums(size);
#pragma omp parallel for
        for (Integer ii = 0; ii < (Integer) size; ii++) {
          dummyElementNums[ii] = ii;
        }
        crf(this->constantInput_[i],resVec,dummyElementNums,generalFactor,spaceFactor,constantFactor,loadedSum,factorSum,countElements);
      } else {
        crf(sol,resVec,elementNums,generalFactor,spaceFactor,constantFactor,loadedSum,factorSum,countElements);
      }
    }
    
    // verbose of the loaded sum
    bool verboseLoadedSum = this->verboseSum_;
    if (verboseLoadedSum && this->hasSpaceInputButTimeFactor_) {
      verboseLoadedSum = !this->initializedConstantInput_;
    }
    if (verboseLoadedSum) {
      std::cout << "Sum of " << this-> solName_ <<" read at " << allSrcRegionNames_ << ": ";
      for (UInt d = 0; d < this->dimDof_; ++d) {
        std::cout << loadedSum[d];
        if (d < this->dimDof_ - 1) {
          std::cout << ",";
        }
      }
      std::cout << std::endl;
    }
    
    // verbose of the sum after applying factors
    bool verboseFactorSum = this->verboseSum_ && this->hasFactor_;
    if (verboseFactorSum && this->hasSpaceInputButTimeFactor_) {
      verboseFactorSum = this->initializedConstantInput_;
    }
    if (verboseFactorSum) {
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
  bool CoefFunctionGridElem<DATA_TYPE>::UpdateSolution(){
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
      bool updated = this->lastStepUpdate_ == std::numeric_limits<unsigned int>::max();
      if (updated) {
        this->ReadSolution(stepValueMap_.begin()->first,this->solVec_);
      }
      this->lastStepUpdate_ = stepnumber;
#ifdef USE_OPENMP
      omp_unset_lock(&updateSolutionLock_);
#endif
      return updated;
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
        while (stepIter != stepValueMap_.end() && stepIter->second < stepValue && notEnd){
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
          if (this->stepNumberInterpolationA_ == std::numeric_limits<unsigned int>::max()) {
            this->solVecInterpolationA_.Resize(this->numEqns_);
            this->solVecInterpolationA_.Init(0.0);
          }
          this->ReadSolution(newA,newAValue,this->solVecInterpolationA_);
          this->stepNumberInterpolationA_ = newA;
        }
        if (this->stepNumberInterpolationB_ != newB) {
          if (this->stepNumberInterpolationB_ == std::numeric_limits<unsigned int>::max()) {
            this->solVecInterpolationB_.Resize(this->numEqns_);
            this->solVecInterpolationB_.Init(0.0);
          }
          this->ReadSolution(newB,newBValue,this->solVecInterpolationB_);
          this->stepNumberInterpolationB_ = newB;
        }
        
#pragma omp parallel for
        for(Integer i=0;i< (Integer) numEqns_;i++){
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


  template class CoefFunctionGridElem<Double>;
  template class CoefFunctionGridElem<Complex>;

} // namespace CoupledField
