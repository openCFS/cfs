#include "CoefFunctionHyst.hh"

// classes for function / spline approximation
#include "Materials/Models/Preisach.hh"
#include "Materials/Models/VectorPreisachSutor.hh"
#include "Materials/Models/VectorPreisachMayergoyz.hh"
#include "FeBasis/FeFunctions.hh"
#include "FeBasis/FeSpace.hh"
#include "Forms/Operators/BaseBOperator.hh"
#include "Domain/CoefFunction/CoefFunctionFormBased.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include <cmath>
#include "Utils/Interpolate1D.hh"
#include "Utils/Timer.hh"

namespace CoupledField {
	DECLARE_LOG(coeffcthyst)
	DEFINE_LOG(coeffcthyst, "coeffcthyst")
  DECLARE_LOG(coeffcthysthelper)
	DEFINE_LOG(coeffcthysthelper, "coeffcthysthelper")
  DECLARE_LOG(coeffcthystdeltamat)
	DEFINE_LOG(coeffcthystdeltamat, "coeffcthystdeltamat")
  
  void CoefFunctionHelper::ComputeVector(Vector<Double>& outputVector,const LocPointMapped& lpm, int timeLevel, int baseSign, 
          std::string vectorName, bool onBoundary ){
    LOG_DBG(coeffcthysthelper) << "+++++ Coef Function Hyst Vector - Get " << vectorName <<" ++++++";
    //std::cout << "Coef Function Hyst RHS Load - Get Vector" << std::endl;
    // return vector that can be put onto rhs of system; 
    // the sign of term shall be such, that it can be added with + in the
    // corresponding pdes
    // Example: in piezo systems we have (among others) -int_Omega (BN)^T*P dOmega 
    //          on rhs, so we would return -P here
    
    Double specificSign;
    
    if(tensorsInitialized_ == false){
      EXCEPTION("Linear tensors have to be initialized by now!");
    }
    
    //    if(hystCoefFunction_->anyMatrixForLocalInversionRequiresComputation() && ((vectorName == "MagPolarization")||(vectorName == "MagMagnetization")) ){
    //      // for calculation of magnetic magnetization, we need to know the material tensor mu/nu
    //      // as M = nu*P = mu^-1*P
    //      // the problem is, that mu itself depends on P in case of magstrict coupling
    //      // and that P requires mu to be computed from B
    //      // > specify mu/nu with previous value of P, i.e. 
    //      //      nu = nu(P_k-1) at iteration k
    //      // to ensure that all integration points/all evaluations during a timestep-iteration
    //      // uses the same state of nu/mu, we precompute this value for all evluation points the first time
    //      // we call ComputeVector or ComputeTensor
    //      PrecomputeMaterialTensorForInverison();
    //    }
    
    //std::cout << "Timelevel: " << timeLevel << " (0 = current, 1 = last it, -1 = last ts)" << std::endl;
    //      if(onBoundary){
    //        if(timeLevelRHS == -1){
    //          //        std::cout << "Get BC wrt to the last ts" << std::endl;
    //        } else if (timeLevelRHS != 0){
    //          //        std::cout << "Get BC wrt to the last "<< timeLevelRHS << "th it" << std::endl;
    //        } else {
    //          //        std::cout << "Put current value of hyst to BC" << std::endl;
    //        }
    //      } else {
    //        if(timeLevelRHS == -1){
    //          //        std::cout << "Get RHS Load wrt to the last ts" << std::endl;
    //        } else if (timeLevelRHS != 0){
    //          //        std::cout << "Get RHS Load wrt to the last "<< timeLevelRHS << "th it" << std::endl;
    //        } else {
    //          //        std::cout << "Put current value of hyst to RHS" << std::endl;
    //        }
    //      }
    
    if(vectorName == "ElecPolarization"){
      //std::cout << "ElecPolarization was requested" << std::endl;
      // rhs (electrostatics, w = testfunction): 
      // + int_Volume nabla(w)^T P dOmega
      // - int_Surface w P^T n dGamma
      // > specificSign = +
      specificSign = 1.0;
      
      outputVector = Vector<Double>(GetVecSize());
      outputVector.Init();
      outputVector = hystCoefFunction_->GetPrecomputedOutputOfHysteresisOperator(lpm,timeLevel, false);
      
      //std::cout << "BaseSign = " << baseSign << std::endl;
      
      outputVector.ScalarMult(specificSign*baseSign);
      
    } else if ((vectorName == "IrrStressForMechPDE") || (vectorName == "IrrStressesPiezo_VectorForm")
            || (vectorName == "IrrStressesMagstrict_VectorForm")){
      /*
       * IrrStressForMechPDE should be the same for piezo and magstrict case
       * piezo: sigma = cS - eE
       * magstrict: sigma = cS - hB
       */
      //      std::cout << vectorName << " was requested" << std::endl;
      // v = mech testfunction
      // on rhs:  - int_Volume (Bv)T[c]S_irr dOmega
      // > basically what is left when S = Bu is decomposed into S_r + S_irr
      //    -div(sigma) > -div([c]S_r) > -div([c](Bu - S_irr)) > -div([c]S_irr) on rhs
      //      but due to integration by parts > div(v)[c]S_irr
      // NOTE: we actually compute the irreversible part of the stress tensor here
      // > we can use this for output, too
      // > in the output case, we need +1 sign however
      if(vectorName == "IrrStressForMechPDE"){
        specificSign = 1.0;
      } else {
        specificSign = 1.0;
      }
      
      Vector<Double> S_irr = hystCoefFunction_->GetIrreversibleStrains(lpm, timeLevel);
      // c is not solution dependet here > take it directly
      outputVector = Vector<Double>(S_irr.GetSize());
      outputVector.Init();
      LOG_DBG(coeffcthysthelper) << "S_irr " << S_irr.ToString();
      LOG_DBG(coeffcthysthelper) << "elastTensor_: " << elastTensor_.ToString();
      elastTensor_.Mult(S_irr,outputVector);
      LOG_DBG(coeffcthysthelper) << "elastTensor_: " << elastTensor_.ToString();
      LOG_DBG(coeffcthysthelper) << "elastTensor_.Mult(S_irr,outputVector): " << outputVector.ToString();
      outputVector.ScalarMult(specificSign*baseSign);
      //      std::cout << "outputVector: " << outputVector.ToString() << std::endl;
    } else if ((vectorName == "CoupledIrrStrainForElecPDE")||(vectorName == "CoupledIrrStrainsPiezo")){
      //      std::cout << vectorName << " was requested" << std::endl;
      // w = elec testfunction
      // on lhs: div( [e](Bu - S_irr ) > -(Bw)T [e] Bu + (Bw)T [e] S_irr
      //  x -1 to get symmetry of coupling matrix
      // (Bw)T [e] Bu - (Bw)T [e] S_irr
      // on rhs:  +(Bw)T[e]S_irr
      // > return [e]S_irr * specific sign
      // note that e might be a function of P itself > use correct timelevel > current
      specificSign = 1.0;
      Matrix<Double> couplTensor = Matrix<Double>(1,1); // will be resized in compute tensor function
      std::string implementationVersion = "none"; // no deltaMat!
      std::string tensorName = "CouplingMechToElec";
      bool transposed = false; // we need e not e^T
      bool rotate = true; // allow rotation if directtion of P changed
      bool useAbs = false; // only for deltaMat
      bool lockPrecomputationAndDeltaMat = false; //!! if true we use old value, but we want current value (timelevel = 0)!
      
      ComputeTensor(couplTensor, lpm, tensorName, implementationVersion, transposed, rotate, useAbs, lockPrecomputationAndDeltaMat );
      
      Vector<Double> S_irr = hystCoefFunction_->GetIrreversibleStrains(lpm, timeLevel);
      // c is not solution dependet here > take it directly
      outputVector = Vector<Double>(GetVecSize());
      outputVector.Init();
      couplTensor.Mult(S_irr,outputVector);
      
      outputVector.ScalarMult(specificSign*baseSign);
     
    } else if ((vectorName == "CoupledIrrStrainForMagPDE")||(vectorName == "CoupledIrrStrainsMagstrict")){
      //      std::cout << vectorName << " was requested" << std::endl;
      // w = mag testfunction
      // on lhs: rot( -[h](Bu - S_irr) > -(rot w)T [h]Bu + (rot w)T [h]S_irr
      // > note: Greens identity for rot does not lead to change in sign in front of
      //          volume term
      //  x -1 to get symmetry of coupling matrix (done in mag PDE)
      // (rot w)T [h] Bu - (rot w)T [h]S_irr
      // on rhs:  +(rot w)T [h]S_irr
      // > return [h]S_irr * specific sign
      // note that h might be a function of P itself > use correct timelevel > current
      specificSign = 1.0;
      Matrix<Double> couplTensor = Matrix<Double>(1,1); // will be resized in compute tensor function
      std::string implementationVersion = "none"; // no deltaMat!
      std::string tensorName = "CouplingMechToMag";
      bool transposed = false; // we need h not h^T
      bool rotate = true; // allow rotation if directtion of P changed
      bool useAbs = false; // only for deltaMat
      bool lockPrecomputationAndDeltaMat = false; //!! if true we use old value, but we want current value (timelevel = 0)!
      
      ComputeTensor(couplTensor, lpm, tensorName, implementationVersion, transposed, rotate, useAbs, lockPrecomputationAndDeltaMat );
      
      Vector<Double> S_irr = hystCoefFunction_->GetIrreversibleStrains(lpm, timeLevel);
      // c is not solution dependet here > take it directly
      outputVector = Vector<Double>(GetVecSize());
      outputVector.Init();
      couplTensor.Mult(S_irr,outputVector);
      
      outputVector.ScalarMult(specificSign*baseSign);
      
      //      std::cout << "outputVector: " << outputVector.ToString() << std::endl;
    } else if(vectorName == "MagPolarization"){
      //      if(hystCoefFunction_->matrixForLocalInversionRequiresComputation(lpm)){
      //        Matrix<Double> tmpTensor = Matrix<Double>(GetVecSize(),GetVecSize());
      //        ComputeTensor(tmpTensor, lpm, 0, "MAG_RELUCTIVITY_FOR_INVERSION", false, false, false);
      //      }
      
      //std::cout << "MagPolarization was requested" << std::endl;
      // rhs (magnetics, w = testfunction): 
      // + int_Volume rot(w)^T M dOmega
      // > specificSign = +
      specificSign = 1.0;
      
      outputVector = Vector<Double>(GetVecSize());
      outputVector.Init();
      // note: getOutputOfHystOperator requires the material tensor for inversion to be set
      // therefore, PrecomputeMaterialTensorForInverison is triggered when MagPolarization is requested
      outputVector = hystCoefFunction_->GetPrecomputedOutputOfHysteresisOperator(lpm,timeLevel,false);
      
      outputVector.ScalarMult(specificSign*baseSign);
//      std::cout << "MagPolarization: " << outputVector.ToString() << std::endl;
    } else if(vectorName == "MagFieldIntensityHyst"){
      //std::cout << "MagMagnetization was requested" << std::endl;
      // rhs (magnetics, w = testfunction): 
      // + int_Volume rot(w)^T M dOmega
      // > specificSign = +
      specificSign = 1.0;
      
      outputVector = Vector<Double>(GetVecSize());
      outputVector.Init();
      // get polarization J_pol
      outputVector = hystCoefFunction_->GetPrecomputedInputToHysteresisOperator(lpm,timeLevel);
      
    } else if(vectorName == "MagMagnetization"){
      //std::cout << "MagMagnetization was requested" << std::endl;
      // rhs (magnetics, w = testfunction): 
      // + int_Volume rot(w)^T M dOmega
      // > specificSign = +
      specificSign = 1.0;
      
      outputVector = Vector<Double>(GetVecSize());
      outputVector.Init();
      // get polarization J_pol
      outputVector = hystCoefFunction_->GetPrecomputedOutputOfHysteresisOperator(lpm,timeLevel,false);
      //std::cout << "Output of hyst operator: " << outputVector.ToString() << std::endl;
      
      // get the matrix that was used to compute P from B to scale P to M
      Matrix<Double> mu = Matrix<Double>(GetVecSize(),GetVecSize());
      Matrix<Double> nu = Matrix<Double>(GetVecSize(),GetVecSize());
      hystCoefFunction_->getMatrixForLocalInversion(lpm, mu, nu);
      LOG_DBG(coeffcthysthelper) << "Used material tensor during inversion (mu):" << mu.ToString();
      LOG_DBG(coeffcthysthelper) << "Used material tensor for computation of magnetization (nu):" << nu.ToString();
      
      // scale with the actual nu
      Vector<Double> tmpVec = Vector<Double>(GetVecSize());
      nu.Mult(outputVector,tmpVec);
      //std::cout << "tmpVec: " << tmpVec.ToString() << std::endl;
      outputVector = tmpVec;
      
      //std::cout << "Return: " << outputVector.ToString() << std::endl;
      outputVector.ScalarMult(specificSign*baseSign);
      //outputVector.ScalarMult(specificSign*baseSign*795774.7155);
    }
    else {
      //std::cout << "Vector " << vectorName << " requested, but not yet implemented " << std::endl;
      EXCEPTION("Vector not implemented yet");
    }
    
    //std::cout << "return: " << std::endl;
    //std::cout << outputVector.ToString() << std::endl;
  }
  
  void CoefFunctionHyst::EvaluateRotDirectionForVectorExtension(){
    if(!POL_useExtension_){
      // this function only is relevant for scalar preisach model with vector extension
      return;
    }
    /*
     * Iterate over all scalar hyst operators and and evaluate the rotation direction
     * with the latest value of the flux or field intensity (E_B_) (i.e. with element solution)
     * note that the storageIdx = operatorIdx corresponds to the midpoint of the
     * element if extended evaluation is used
     */
    if(POL_setWithFlux_){
      
      Vector<Double> latestFlux;
      for(UInt i = 0; i < numHystOperators_; i++){
        latestFlux = E_B_[i];
        hyst_->UpdateRotationStateWithFluxDensity(latestFlux,matrixForInversion_[i],i);
        hyst_->EvaluateRotationState(i);
      }
    } else {
      Vector<Double> latestField;
      for(UInt i = 0; i < numHystOperators_; i++){
        latestField = E_B_[i];
        hyst_->UpdateRotationStateWithFieldIntensity(latestField,i);
        hyst_->EvaluateRotationState(i);
      }
    }
  }
  
  void CoefFunctionHelper::PrecomputeMaterialTensorForInverison(){
    LOG_DBG(coeffcthysthelper) << "CoefFunctionHelper::PrecomputeMaterialTensorForInverison()";
    
    /*
     * Iterate over all integration points used/managed by the hyst operator and
     * for each of these LPM set the matrix for iversion, i.e. mu/nu in magnetic case
     */
    std::map<UInt, LocPointMapped > allLPM;
    std::map<UInt, LocPointMapped > midpointLPM;
    hystCoefFunction_->getLPMMaps(allLPM, midpointLPM);
    
    UInt dim = hystCoefFunction_->GetVecSize();    
    UInt size = allLPM.size();
    
#ifdef USE_OPENMP
#pragma omp parallel num_threads(CFS_NUM_THREADS)
    {
      Matrix<Double> nu = Matrix<Double>(dim,dim);  
      Matrix<Double> mu = Matrix<Double>(dim,dim);
      nu.Init();
      mu.Init();
      // we will make use of the computeTensor function to avoid double code
      // however, to make sure that we do neither get a deltaMat nor that we
      // call this function here recursively, we have to adapt some flags
      bool lockPrecomputationAndDeltaMat = true;
      bool reuse = false;
      
      //std::cout << "Running parallel" << std::endl;
      UInt numT = CFS_NUM_THREADS;
      UInt aThread = omp_get_thread_num();
      //std::cout << "My number: " << aThread << std::endl;
      UInt chunksize = std::floor(size/numT);
      
      std::map<UInt, LocPointMapped >::iterator LPMit;
      std::map<UInt, LocPointMapped >::iterator LPMitStart = allLPM.begin();
      std::map<UInt, LocPointMapped >::iterator LPMitEnd;
      
      std::advance(LPMitStart, aThread*chunksize);
      
      if(aThread == numT-1){
        LPMitEnd = allLPM.end();
      } else {
        LPMitEnd = LPMitStart;
        std::advance(LPMitEnd, chunksize);
      }  
      
#pragma omp barrier
      for(LPMit = LPMitStart; LPMit != LPMitEnd; LPMit++){
        if((isCoupled_ == false)&&(materialTensorsSetOnce_)){
          LOG_DBG(coeffcthysthelper) << "Single field case > small signal tensor cannot depend on coupling";
          reuse = true;
          // calling setMatrixForLocalInversion with reuse = true will simply set the flag
          // matrixForInversionWasSet_ to true
          hystCoefFunction_->setMatrixForLocalInversion(mu, nu, LPMit->first, reuse);
        } else if((strainForm_ == -1 )&&(materialTensorsSetOnce_)){
          LOG_DBG(coeffcthysthelper) << "Small signal tensor shall not be used; mu/nu independent on M/P; reuse same value";
          reuse = true;
          // calling setMatrixForLocalInversion with reuse = true will simply set the flag
          // matrixForInversionWasSet_ to true
          hystCoefFunction_->setMatrixForLocalInversion(mu, nu, LPMit->first, reuse);
        } else {
          LOG_DBG(coeffcthysthelper) << "Pre-Compute material tensor for local inversion";
          // now we basically compute the reluctivity at each LPM but not with the
          // current value of P (as its computation would require the tensor) but with
          // the previous value > by setting lockPrecomputationAndDeltaMat to true, three options are set
          //    a) no delta Matrix
          //    b) no call to this function (would otherwise end in endless recursion)
          //    c) evaluate all tensors at last iteration
          ComputeTensor(nu, LPMit->second, "Reluctivity", "none", false, true, false, lockPrecomputationAndDeltaMat);
          //ComputeTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm, , 
          //std::string tensorName, std::string implementationVersion, bool transposed, bool rotate, bool useAbs, bool lockPrecomputationAndDeltaMat )
          
          // now we have nu but for inversion we use mu > compute inverse
          nu.Invert(mu);
          
          hystCoefFunction_->setMatrixForLocalInversion(mu, nu, LPMit->first, reuse);
          LOG_DBG(coeffcthysthelper) << "Computed material tensor for storageIdx="<<LPMit->first<<": "<<mu.ToString();
        }
      }
#pragma omp barrier
    }
#else
    {
      Matrix<Double> nu = Matrix<Double>(dim,dim);  
      Matrix<Double> mu = Matrix<Double>(dim,dim);
      nu.Init();
      mu.Init();
      // we will make use of the computeTensor function to avoid double code
      // however, to make sure that we do neither get a deltaMat nor that we
      // call this function here recursively, we have to adapt some flags
      bool lockPrecomputationAndDeltaMat = true;
      bool reuse = false;
      
      //std::cout << "Running serial" << std::endl;
      std::map<UInt, LocPointMapped >::iterator LPMit;
      for(LPMit = allLPM.begin(); LPMit != allLPM.end(); LPMit++){
        if((COUPLED_inXMLFile_ == false)&&(materialTensorsSetOnce_)){
          LOG_DBG(coeffcthysthelper) << "Single field case > small signal tensor cannot depend on coupling";
          reuse = true;
          // calling setMatrixForLocalInversion with reuse = true will simply set the flag
          // matrixForInversionWasSet_ to true
          hystCoefFunction_->setMatrixForLocalInversion(mu, nu, LPMit->first, reuse);
        } else if((strainForm_ == -1 )&&(materialTensorsSetOnce_)){
          LOG_DBG(coeffcthysthelper) << "Small signal tensor shall not be used; mu/nu independent on M/P; reuse same value";
          reuse = true;
          // calling setMatrixForLocalInversion with reuse = true will simply set the flag
          // matrixForInversionWasSet_ to true
          hystCoefFunction_->setMatrixForLocalInversion(mu, nu, LPMit->first, reuse);
        } else {
          LOG_DBG(coeffcthysthelper) << "Pre-Compute material tensor for local inversion";
          // now we basically compute the reluctivity at each LPM but not with the
          // current value of P (as its computation would require the tensor) but with
          // the previous value > by setting lockPrecomputationAndDeltaMat to true, three options are set
          //    a) no delta Matrix
          //    b) no call to this function (would otherwise end in endless recursion)
          //    c) evaluate all tensors at last iteration
          ComputeTensor(nu, LPMit->second, "Reluctivity", "none", false, true, false, lockPrecomputationAndDeltaMat);
          //ComputeTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm, int timeLevel_to_diff, 
          //std::string tensorName, std::string implementationVersion, bool transposed, bool rotate, bool useAbs, bool lockPrecomputationAndDeltaMat )
          
          // now we have nu but for inversion we use mu > compute inverse
          nu.Invert(mu);
          
          hystCoefFunction_->setMatrixForLocalInversion(mu, nu, LPMit->first, reuse);
          LOG_DBG(coeffcthysthelper) << "Computed material tensor for storageIdx="<<LPMit->first<<": "<<mu.ToString();
        }
      }
    }
#endif
    
    materialTensorsSetOnce_ = true;
  }
  
  void CoefFunctionHelper::ComputeTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm, 
          std::string tensorName, std::string implementationVersion, bool transposed, bool rotate, bool useAbs, bool lockPrecomputationAndDeltaMat ){
    LOG_DBG(coeffcthysthelper) << "+++++ Coef Function Hyst Mat - Get " << tensorName <<" ++++++";
    
    /*
     * IMPORTANT NOTE TO SIGN CONVENTION:
     *  all tensors will be returned with a "+" sign in front of the actual base tensor
     *  even if the later equation system requires a "-" sign
     *  Example mag-strict
     *    T = cS - hB = c(Bu) - h(rotA)
     *    H =  -hS + nuB = -h(Bu) + nu(rotA)
     *  > this function returns +c,+h,+nu even though +c,-h,+nu will be used later
     *  Example pizeo
     *    T = cS - eE = c(Bu) + e(gradV)
     *    D = eS + epsE = e(Bu) - eps(gradV) 
     *  > this function returns +c,+e,+eps event though +c,+e,+e,-eps will be used later
     * 
     *  If delta matrices are involved, their sign will be chose such, that the former point
     *  holds.
     *  Example mag-strict including deltaMat 
     *    nuMod = nu - dM/dB + h*dS/dB
     *    hMod = h - c*dS/dB
     *  > according to formula, we would need c*dS/dB on lhs but due to the "-" that needs
     *    to be multiplied to h, we have to subtract this contribution here
     *  Exmaple piezo including deltaMat
     *    epsMod = eps + dP/dE - e*dS/dE
     *    eMod = e + c*dS/dE
     *  > according to formula, we need c*dS/dE on lhs; as e is positve (when usign E = -gradV)
     *    we simply add c*dS/dE
     * 
     * CONSEQUENCE: the actual sign for the tensors has to be implemented in the corresponding classes
     *  e.g. in elecPDE we have to assemble -eps (is already the case for piezos) and in magstrictPDE
     *  we have to assemble the coupling terms with -1
     * 
     */
    //    std::cout << "+++++ Coef Function Hyst Mat - Compute Tensor 1 ++++++" << std::endl;
    UInt numCols, numRows;
    // GetTensorSize will return the size of the actual tensor that shall
    // be returned, i.e. if it is the transposed of the couplTensor_
    // numCols = numRows(coupTensor_)
    // > NEW: this function can now be called via CoefFunctionHelper::PrecomputeMaterialTensorForInverison()
    //        and thus from CoefFunctionRHSload, CoefFunctionOutput and so on
    //        therewith, getTensor size of these functions will be called which lead to errors as this was not 
    //        meant to happen
    GetTensorSize(numRows,numCols,tensorName);
    Matrix<Double> tmp;
    
    outputTensor.Resize(numRows,numCols);
    outputTensor.Init();
    bool alreadyTransposed = false;
    /*
     * THE NEXT TWO FUNCTIONS LEAD TO ISSUES WITH RACE CONDITIONS
     */
    if(tensorsInitialized_ == false){
      EXCEPTION("Linear tensors have to be initialized by now!");
    }
    //      if(tensorsInitialized_ == false){
    //        LOG_DBG(coeffcthysthelper) << "Initialize linear tensors";
    //        InitLinearTensors(lpm);
    //      }
    //      std::cout << "+++++ Coef Function Hyst Mat - Compute Tensor 2 ++++++" << std::endl;
    //      if(hystCoefFunction_->anyMatrixForLocalInversionRequiresComputation() && 
    //              ((tensorName == "Reluctivity")||(tensorName == "CouplingMechToMag")||(tensorName == "CouplingMagToMech")) ){
    //        LOG_DBG(coeffcthysthelper) << "Set matrix for inversion";
    //        // for calculation of magnetic magnetization, we need to know the material tensor mu/nu
    //        // as M = nu*P = mu^-1*P
    //        // the problem is, that mu itself depends on P in case of magstrict coupling
    //        // and that P requires mu to be computed from B
    //        // > specify mu/nu with previous value of P, i.e. 
    //        //      nu = nu(P_k-1) at iteration k
    //        // to ensure that all integration points/all evaluations during a timestep-iteration
    //        // uses the same state of nu/mu, we precompute this value for all evluation points the first time
    //        // we call ComputeVector or ComputeTensor
    //        if(lockPrecomputationAndDeltaMat == false){
    //          LOG_DBG(coeffcthysthelper) << "PrecomputeMaterialTensorForInverison unlocked";
    //          PrecomputeMaterialTensorForInverison();
    //        }
    //      }
    //      std::cout << "+++++ Coef Function Hyst Mat - Compute Tensor 3 ++++++" << std::endl;
    if(transposed){
      tmp = Matrix<Double>(numCols,numRows);
    } else {
      tmp = Matrix<Double>(numRows,numCols);
    }
    tmp.Init();
    
    int deltaForm_ = hystCoefFunction_->GetDeltaForm();
    bool deltaFormActive_ = hystCoefFunction_->deltaMatActive();
    if( (implementationVersion == "none") || (lockPrecomputationAndDeltaMat == true) ){
      deltaFormActive_ = false;
    }
    //      std::cout << "+++++ Coef Function Hyst Mat - Compute Tensor 4 ++++++" << std::endl;
    // always compute deltaMat from current value to a previous one (timeLevel_to_diff)
    // except for the precomputation case (which calls this function with the flag 
    // lockPrecomputationAndDeltaMat = true)
    int timelevel_cur;
    if(lockPrecomputationAndDeltaMat == true){
      timelevel_cur = 1; // value from previous iteration
    } else {
      timelevel_cur = 0; // current value
    }
    //    int timeLevel_to_diff = GetTimeLevel("DeltaMat");
    
    // new flags
    bool deltaMat_Pol_active = true;
    bool deltaMat_Strain_active = true;
    bool deltaMat_Coupling_active = true;
    
    int timeLevelDeltaMat_Pol = GetTimeLevel("DeltaMatPol");
    int timeLevelDeltaMat_Strain = GetTimeLevel("DeltaMatStrain");
    int timeLevelDeltaMat_Coupling = GetTimeLevel("DeltaMatCoupling");
    
    if( (timeLevelDeltaMat_Pol == 0)||(deltaFormActive_==false) ){
      deltaMat_Pol_active = false;
    }
    if( (timeLevelDeltaMat_Strain == 0)||(deltaFormActive_==false) ){
      deltaMat_Strain_active = false;
    }
    if( (timeLevelDeltaMat_Coupling == 0)||(deltaFormActive_==false) ){
      deltaMat_Coupling_active = false;
    }
    
    if( (tensorName == "IrrStressesPiezo_TensorForm")||(tensorName == "IrrStressesMagstrict_TensorForm") ){
      Vector<Double> Si_vector = hystCoefFunction_->GetIrreversibleStrains(lpm, 0);
      Vector<Double> tmpVector = Vector<Double>(Si_vector.GetSize());
      elastTensor_.Mult(Si_vector,tmpVector);      
      tmp = hystCoefFunction_->ConvertFromVoigtToTensor(tmpVector);
      
    } else if((tensorName == "IrrStrainsPiezo_TensorForm")||(tensorName == "IrrStrainsMagstrict_TensorForm")){
      // always compute the current timelevel > 0
      // > tensor form is only used for output; otherwise use vector form
      tmp = hystCoefFunction_->GetIrreversibleStrainTensor(lpm, 0);
    } else if(tensorName == "Permittivity"){
      //        std::cout << "Get Permittivity" << std::endl;
      /*
       * The following cases are possible:
       *  I. pure electrostatics:
       *     a) no-deltaMat > return epsS
       *     b) deltamat    > return epsS + dP/dE
       *  II. coupled piezo-electric case
       *     a) e-form as basis
       *        return epsS + dP/dE - e*dS/dE
       *     b) d-form as basis
       *        return epsT - d [c] d^T + dP/dE - d [c] dS/dE
       */
      
      // in case of piezoelectricity, we have to check for e-form or d-form
      // in case of e-form, fieldTensor is seens (and set to be) epsS which needs
      // no further treatment
      Matrix<Double> e_scaled;
      Matrix<Double> tmp2 = Matrix<Double>(numRows,numCols);
      if(strainForm_ == 1){
        //std::cout << "Get Permittivity - Compute epsS from epsT" << std::endl;
        // d-form shall be used as basis; fieldTensor such represents epsT
        // to get the required epsS we have to compute
        //  epsS = epsT - d(P)[c]d(P)^T
        // where d(P) is the scaled and rotated coupling tensor in d-form
        Matrix<Double> d_scaled, d_scaled_transposed;
        
        hystCoefFunction_->GetScaledAndRotatedCouplingTensor(lpm,couplTensor_,d_scaled,timelevel_cur,rotate);
        
        // e = d*c
        UInt numRows,numCols;
        numRows = d_scaled.GetNumRows();
        numCols = d_scaled.GetNumCols();
        e_scaled = Matrix<Double>(numRows,numCols);
        d_scaled.Mult(elastTensor_,e_scaled);
        d_scaled_transposed = Matrix<Double>(numCols,numRows);
        d_scaled.Transpose(d_scaled_transposed);
        // tmp2 = e*d^T = dcd^T
        e_scaled.Mult(d_scaled_transposed,tmp2);
        
        // tmp = epsT
        tmp = fieldTensor_;
        
        // tmp = epsT - d*c*d^T
        tmp.Add(-1.0,tmp2);
        
      } else {
        //std::cout << "Get Permittivity - Take epsS directly" << std::endl;
        // uncoupled case or e-form
        // just take epsS directly
        tmp = fieldTensor_;
        // also get e in case we use deltaMat for strains, too
        hystCoefFunction_->GetScaledAndRotatedCouplingTensor(lpm,couplTensor_,e_scaled,timelevel_cur,rotate);  
      }      
      // check if deltaMat shall be added
      // here we have two flags:
      //  deltaForm_ is used to indicate if we are using a delta formulation in general
      //  deltaFormActive_ indicates if we actually need this deltaMatrix for the current
      //    evaluation (the issue is, that we need a deltaMatrix on the lhs and a non-deltaMatrix
      //    on the rhs
      //      if(deltaFormActive_ && (deltaForm_ != 0) && (lockPrecomputationAndDeltaMat == false) ) {
      if(deltaMat_Pol_active){
        // std::cout << "Get Permittivity - Compute DeltaMatrix" << std::endl;
        //        std::cout << "Compute DeltaMatrix" << std::endl;
        bool useStrains = false;
        
        // deltaMat will be computed using the current value (timelevel_cur = 0)
        // and the value at timelevel (-1 > last ts; +1 > last iteration)
        Matrix<Double> deltaMat = hystCoefFunction_->GetDeltaMat(lpm, timelevel_cur, timeLevelDeltaMat_Pol, useStrains, useAbs,implementationVersion);
        //        std::cout << "deltaMat elec = " << deltaMat.ToString() << std::endl;
        tmp.Add(1.0,deltaMat);
        //std::cout << "DeltaMat: " << deltaMat.ToString() << std::endl;
        if( (strainForm_ != -1) && (deltaMat_Strain_active) ){
          //std::cout << "Get Permittivity - Add dS/dE" << std::endl;
          // coupled case
          // here we have to add -e*dS/dE in addition
          useStrains = true;
          Matrix<Double> deltaMat_strains = hystCoefFunction_->GetDeltaMat(lpm, timelevel_cur, timeLevelDeltaMat_Strain, useStrains, useAbs,implementationVersion);
          //          std::cout << "deltaMat strains = " << deltaMat_strains.ToString() << std::endl;
          e_scaled.Mult(deltaMat_strains,tmp2);
          
          tmp.Add(-1.0,tmp2);
          //std::cout << "DeltaMatStrain: " << deltaMat_strains.ToString() << std::endl;
        }
        
      } 
      //          else if (deltaForm_ != 0){
      //                  std::cout << "deltaFormSet_ == true, but deltaMatNotActive!" << std::endl;
      //        } else {
      //                  std::cout << "deltaFormSet_ == false" << std::endl;
      //        }
      
    } else if ((tensorName == "CouplingMechToElec")||(tensorName == "CouplingMechToMag")){
      //      std::cout << "ComputeMechToElec" << std::endl;
      Matrix<Double> rotatedCouplTensor;
      
      if(strainForm_ == -1){
        // small signal tensor defined, but shall not be used > return empty matrix/vector
        rotatedCouplTensor.Init();
        tmp = rotatedCouplTensor;
      } else if(strainForm_ == 1){
        // use d-form/g-form as basis, i.e. the followings steps have to be applied
        // 1. scale and rotate d/g
        hystCoefFunction_->GetScaledAndRotatedCouplingTensor(lpm,couplTensor_,rotatedCouplTensor,timelevel_cur,rotate);
        // 2. compute e/h from d/g
        // e = d*c ; h = g*c
        rotatedCouplTensor.Mult(elastTensor_,tmp);
      } else {
        // use e-form/h-form as basis, i.e. the following step has to be done
        // 1. acale and rotate e
        hystCoefFunction_->GetScaledAndRotatedCouplingTensor(lpm,couplTensor_,rotatedCouplTensor,timelevel_cur,rotate);
        tmp = rotatedCouplTensor;
      }
      
    } else if ((tensorName == "CouplingElecToMech")||(tensorName == "CouplingMagToMech")){ 
      // this is basically the same tensor as mechToElec except for the case of
      // deltaFormulation; in the later case, we have to add c*dS/dE    
      //      std::cout << "ComputeElecToMech" << std::endl;
      Matrix<Double> rotatedCouplTensor;
      
      if(strainForm_ == -1){
        // small signal tensor defined, but shall not be used > return empty matrix/vector
        rotatedCouplTensor.Init();
        tmp = rotatedCouplTensor;
      } else if(strainForm_ == 1){
        // use d-form as basis, i.e. the followings steps have to be applied
        // 1. scale and rotate d/g
        hystCoefFunction_->GetScaledAndRotatedCouplingTensor(lpm,couplTensor_,rotatedCouplTensor,timelevel_cur,rotate);
        // 2. compute e/h from d/g
        // e = d*c ; h = g*c
        rotatedCouplTensor.Mult(elastTensor_,tmp);
      } else {
        // use e-form/h-form as basis, i.e. the following step has to be done
        // 1. acale and rotate e/h
        hystCoefFunction_->GetScaledAndRotatedCouplingTensor(lpm,couplTensor_,rotatedCouplTensor,timelevel_cur,rotate);
        tmp = rotatedCouplTensor;
      }
      
      //      if(deltaFormActive_ && (deltaForm_ != 0) && (lockPrecomputationAndDeltaMat == false) ) {
      if( (deltaMat_Coupling_active)&&(strainForm_ != -1) ){
        LOG_DBG(coeffcthysthelper) << "Compute DeltaMatrix";
        //        std::cout << "Compute DeltaMatrix" << std::endl;
        bool useStrains = true;
        
        // deltaMat will be computed using the current value (timelevel_cur = 0)
        // and the value at timelevel (-1 > last ts; +1 > last iteration)
        Matrix<Double> deltaMat = hystCoefFunction_->GetDeltaMat(lpm, timelevel_cur, timeLevelDeltaMat_Coupling, useStrains, useAbs,implementationVersion);
        //        std::cout << "deltaMat " << deltaMat.ToString() << std::endl;
        UInt numRows,numCols;
        numRows = tmp.GetNumRows();
        numCols = tmp.GetNumCols();
        Matrix<Double> tmp2 = Matrix<Double>(numCols,numRows);
        // tmp2 = dS/dE*c / dS/dB*c
        // > sollte eigentlich c*dS/dx sein
        //deltaMat.Mult(elastTensor_,tmp2);
        elastTensor_.Mult(deltaMat,tmp2);
        
        // tmp = e_scaled_rotated + dS/dE*c / -dS/dB*c
        // > does not work! e is 2x3 / 3x6 but c*dS/dE is 3x2 / 6x3
        // > reason: in pde we need e^T + c*dS/dE instead of (e + c*dS/dE)^T
        // > transpose first, then add deltaMat
        // > transpose directly into outputTensor
        if(transposed){
          //std::cout << "Perform transpose" << std::endl;
          tmp.Transpose(outputTensor);
          alreadyTransposed = true;
          //std::cout << "Transposed tensor: " << tmp.ToString() << std::endl;
        } else {
          EXCEPTION("We should require transpose here!");
        }
        
        if(tensorName == "CouplingMagToMech"){
          // NOTE: we need +dS/dB on lhs exactly like in piezo-case; however,
          // the coupling term here gets multiplied by -1 in the magnetostrictive case
          // as we need -h + dS/dB*c; therewith dS/dB has to be subtracted here to get the correct +1 sign later
          outputTensor.Add(-1.0,tmp2);
        } else {
          outputTensor.Add(1.0,tmp2);
        }
      }  
    } else if(tensorName == "Reluctivity"){
      //std::cout << "Get Reluctivity" << std::endl;
      /*
       * The following cases are possible:
       *  I. pure magnetics:
       *     a) no-deltaMat > return nuS
       *     b) deltamat    > return nuS - dM/dB
       *  II. coupled magnetostrictive case
       *     a) h-form as basis
       *        return nuS -dM/dB + h*dS/dB
       *     b) g-form as basis
       *        return nuT + g [c] g^T - dM/dB + g [c] dS/dB
       */
      
      // in case of magnetostriction, we have to check for h-form or g-form
      // in case of h-form, fieldTensor is seen (and set to be) nuS which needs
      // no further treatment
      Matrix<Double> h_scaled;
      Matrix<Double> tmp2;
      if(strainForm_ == 1){
        //std::cout << "Get Reluctivity - Compute nuS from nuT" << std::endl;
        // g-form shall be used as basis; fieldTensor such represents nuT
        // to get the required nuS we have to compute
        // (note that compared to piezos, we have a + here)
        //  nuS = nuT + g(P)[c]g(P)^T
        // where g(P) is the scaled and rotated coupling tensor in g-form
        Matrix<Double> g_scaled, g_scaled_transposed;
        
        hystCoefFunction_->GetScaledAndRotatedCouplingTensor(lpm,couplTensor_,g_scaled,timelevel_cur,rotate);
        
        // h = g*c
        g_scaled.Mult(elastTensor_,h_scaled);
        h_scaled.Transpose(g_scaled_transposed);
        // tmp2 = h*g^T = gcg^T
        h_scaled.Mult(g_scaled_transposed,tmp2);
        
        // tmp = nuT
        tmp = fieldTensor_;
        
        // tmp = nuT + g*c*g^T
        tmp.Add(1.0,tmp2);
        
      } else {
        //std::cout << "Get Reluctivity - Take nuS directly" << std::endl;
        // uncoupled case or h-form
        // just take nuS directly
        tmp = fieldTensor_;
      }
      
      // check if deltaMat shall be added
      // here we have two flags:
      //  deltaForm_ is used to indicate if we are using a delta formulation in general
      //  deltaFormActive_ indicates if we actually need this deltaMatrix for the current
      //    evaluation (the issue is, that we need a deltaMatrix on the lhs and a non-deltaMatrix
      //    on the rhs
      if(deltaMat_Pol_active){
        //      if(deltaFormActive_ && (deltaForm_ != 0) && (lockPrecomputationAndDeltaMat == false) ) {
        //std::cout << "Get Reluctivity - Compute DeltaMatrix" << std::endl;
        //        std::cout << "Compute DeltaMatrix" << std::endl;
        bool useStrains = false;
        
        // deltaMat will be computed using the current value (timelevel_cur = 0)
        // and the value at timelevel (-1 > last ts; +1 > last iteration)
        Matrix<Double> deltaMat = hystCoefFunction_->GetDeltaMat(lpm, timelevel_cur, timeLevelDeltaMat_Pol, useStrains, useAbs,implementationVersion);
        
        // note: deltaMat will compute dP/dB, but we want to have -dM/dB
        // i.e. add -nu*deltaMat to tmp
        // in other words we want to have:
        // dH = nu*(Identiy - deltaMat)*dB
        Matrix<Double> tmp3;
        tmp3.Resize(numRows,numCols);
        tmp.Mult(deltaMat,tmp3);
        
        tmp.Add(-1.0,tmp3);
        
        //tmp.Add(-795774.7155,deltaMat);
        //std::cout << "DeltaMat (dP/dB): " << deltaMat.ToString() << std::endl;
        if( (strainForm_ != -1) && (deltaMat_Strain_active) ){
          //        if(strainForm_ != -1){
          //std::cout << "Get Reluctivity - Add dS/dB" << std::endl;
          // coupled case
          // here we have to add +g[c] dS/dB in addition
          useStrains = true;
          Matrix<Double> deltaMat_strains = hystCoefFunction_->GetDeltaMat(lpm, timelevel_cur, timeLevelDeltaMat_Strain, useStrains, useAbs,implementationVersion);
          h_scaled.Mult(deltaMat_strains,tmp2);
          
          tmp.Add(1.0,tmp2);
          //std::cout << "DeltaMatStrain: " << deltaMat_strains.ToString() << std::endl;
        }
        
      } else if (deltaForm_ != 0){
        //        std::cout << "deltaFormSet_ == true, but deltaMatNotActive!" << std::endl;
      } else {
        //        std::cout << "deltaFormSet_ == false" << std::endl;
      }
      
    } else {
      //std::cout << "Tensor " << tensorName << " requested, but not yet implemented " << std::endl;
      EXCEPTION("Tensor not implemented yet");
    }
    
    if(alreadyTransposed == false){
      if(transposed){
        //std::cout << "Perform transpose" << std::endl;
        tmp.Transpose(outputTensor);
        //std::cout << "Transposed tensor: " << tmp.ToString() << std::endl;
      } else {        
        outputTensor = tmp;
      }
    } 
    
    LOG_DBG(coeffcthysthelper) << "Computed material tensor:" << outputTensor.ToString();
    //      std::cout << "+++++ Coef Function Hyst Mat - Compute Tensor END ++++++" << std::endl;
    //std::cout << "Return the following tensor: " << outputTensor.ToString() << std::endl;
    //    std::cout << "Computed tensor - " << tensorName << " - " << std::endl;
    //    std::cout << outputTensor.ToString() << std::endl;        
  }
  
	CoefFunctionHyst::CoefFunctionHyst(BaseMaterial * const material,
          shared_ptr<ElemList> actSDList,
          PtrCoefFct dependency1,
          SubTensorType tensorType, MaterialType matType, shared_ptr<FeSpace> ptFeSpace) 
  : CoefFunction() {
    
		dependCoef1Surf_ = NULL;
		Init(material, actSDList, dependency1, tensorType, matType, ptFeSpace);
	}
  
	CoefFunctionHyst::CoefFunctionHyst(BaseMaterial * const material,
          shared_ptr<ElemList> actSDList,
          PtrCoefFct dependency1, PtrCoefFct dependCoef1Surf,
          SubTensorType tensorType, MaterialType matType, shared_ptr<FeSpace> ptFeSpace) 
  : CoefFunction() {
    
    //		WARN("Currently we support only single-input hysteresis operators! The second dependency is only allowed if it is the surface version of the original dependency!")
		dependCoef1Surf_ = dependCoef1Surf;
		Init(material, actSDList, dependency1, tensorType, matType, ptFeSpace);
	}
  
	void CoefFunctionHyst::Init(BaseMaterial * const material,
          shared_ptr<ElemList> actSDList,
          PtrCoefFct dependency1,
          SubTensorType tensorType,
          MaterialType matType,
          shared_ptr<FeSpace> ptFeSpace) {
    
		// this type of coefficient is nonlinear (i.e. solution dependent)
		dimType_ = VECTOR;
		dependType_ = SOLUTION;
		isAnalytic_ = false;
		isComplex_ = false;
		dependCoef1_ = dependency1;
    storageInitialized_ = false;    
		tensorType_ = tensorType;
    hystHelper_ = NULL;
    
    // dim_ is the dim of the output retrieved by GetVector
		// not tha same as Preisach_Dim that determines whether we use scalar or vector model
		dim_ = dependency1->GetVecSize();
    
    if(tensorType_ == AXI){
      EXCEPTION("Coef functions hyst does only support plane 2d and full 3d models (at the moment)")
    }
    
		matType_ = matType;
		material_ = material;
		ptFeSpace_ = ptFeSpace;
    
    ReadInMaterial(material);
    
		/*
     * for performance measurement
     */
		timer_ = new Timer();
		totalCallingCounter_ = 0;
		totalEvaluationCounter_ = 0;
		avgEvaluationTime_ = 0.0;
		totalEvaluationTime_ = 0.0;
    
    // for slope estimation
    // can be turned on in solvestep to allow for better initial step
    ES_includeStrains_ = false;
    ES_useAbs_ = false;
    ES_implementationVersion_ = "Division";
    ES_steppingLength_ = 1e-10;
    ES_scaling_ = 1.0;
    
		// for interaction with coefFunctionHystMat
		deltaMatActive_ = true;
		hystItself_ = PtrCoefFct(this);
    
    // per default the coef functtions for coupling are NULL
    // for coupled pdes, these functions have to be set BEFORE both coupled
    // and single field integrators get defined!
    elastTensorFct_ = NULL;
    couplTensorFct_ = NULL;
    
    /*
     * Default case: deltaMat towards last TS
     * > material tensors at current value (0)
     * > deltaMat towards last ts (-1)
     * > rhs hyst on last ts (-1)
     * > bc on current value (0)
     * > output on current value (0)
     */
    RUN_deltaForm_ = -1;
		timeLevel_Mat_ = 0; // i.e. time step for which to evaluate material tensors (eps, nu, ...)
		timeLevel_deltaMatPol_ = -1; // timestep to which delta is computed
    timeLevel_deltaMatStrain_ = 0; // no deltaMat for strain (deltaS/deltaE)
		timeLevel_deltaMatCoupling_ = 0; // no deltaMat for coupling ([e]*deltaS/deltaE)
    
		timeLevel_rhsPol_ = -1; // rhs term has to be last ts value as delta is towards that value
    timeLevel_rhsStrain_ = 0; // current value of S_irr
    timeLevel_rhsCoupling_ = 0;  // current value of [e]S_irr
		timeLevel_boundaryTerms_ = 0; // BC term has to use the last known, i.e. current value
		timeLevel_results_ = 0; // output should always be the current value
    forceCurrentTS_ = false;
    includeOldDelta_ = 0.0;
    
		/*
     * set initial values for runtime dependent parameter
     */
		// before anything can be computed, a system has to be assembled first -> 1
		RUN_evaluationPurpose_ = 1;
    
		// do not output switching state; that is very costly and should only be done
		// for debugging or special figure computation
		RUN_allowBMP_ = false;
    
    hystOperatorLocked_ = false;
    
		// is set to false, the hyst operator will keep its rotation state
		// if the initial state is unset, this initial rotation state will be 0 0
		// unless the flag is set to true, it will stay 0 0 and so the output of the
		// hyst operator will be 0 0, too
		// > at least during the first inputs, this flag should be true so that the
		//   rotation list gets initialized properly
		// > tests showed that this flag should always be true
    // 1.6.2018 removed as it actually always had to be true
    //		RUN_overwriteDirection_ = true;
		
		/*
     * set temporary values for xml dependent parameter
     */
		// currently, we have no direct handle of the xml input during the creation
		// of this object; the initialization of the XML dependent parameter is therefore
		// postponed to the first call of StdSolveStep (prior to that we do not need the
		// paramters either way)
		XMLParameterSet_ = false;
    
		XML_EvaluationDepth_ = 0;
		XML_performanceMeasurement_ = 0;
    
		/*
     * get elements and integration points
     */
		// set map: global to local element number
		EntityIterator it = actSDList->GetIterator();
		UInt iel = 0;
		UInt globalElNr;
		for (it.Begin(); !it.IsEnd(); it++, iel++) {
			globalElNr = it.GetElem()->elemNum;
			globalElem2Local_[globalElNr] = iel;
			globalElemOnSurf_[globalElNr] = false;
		}
    
		//store subdomain list of elements
		SDList_List_.push_back(actSDList);
		numElemSD_ = actSDList->GetSize();
    
		// pick out the first element (even though any of them would do as they share the
		// same material) and extract midpoint
		it.Begin();
		const Elem * el = it.GetElem();
		LocPoint lp = Elem::shapes[el->type].midPointCoord;
		LocPointMapped lpm;
    bool updateShapeMap = false;
		shared_ptr<ElemShapeMap> esm = it.GetGrid()->GetElemShapeMap(el, updateShapeMap);
		lpm.Set(lp, esm, 0.0);
    
		// get number of integration points
		IntegOrder order;
		IntScheme::IntegMethod method;
		StdVector<LocPoint> intPoints;
		StdVector<Double> weights;
    
		ptFeSpace_->GetFe(it, method, order);
		// store method and order for later on so that we can retrieve the
		// integration points without getting the fe
		// > this is useful for surface elements where we do not necessarily know the
		//   neighboring volume element (GetFe will fail in this case)
		IntegMethod_ = method;
		IntegOrder_ = order;
		
		ptFeSpace_->GetIntScheme()->GetIntPoints(Elem::GetShapeType(el->type), method, order,
            intPoints, weights);
    
		
		numIntegrationPoints_ = intPoints.GetSize();
    
		// create a map that gives to each integration point + to midpoint an index
		// we cannot have a map of loc points directly (not sortable) but each point
		// has also a number which we can use for sorting; midpoint gets index -1
    // > storageIdx of midpoint = operatorIdx 
		locPointIndices_.insert(std::pair<int, UInt>(-1, 0));
    
		for (UInt i = 0; i < numIntegrationPoints_; i++) {
			//std::cout << intPoints[i].number << std::endl;
			locPointIndices_.insert(std::pair<int, UInt>(intPoints[i].number, i + 1));
		}
    
		//std::cout << "NumIntegrationPoints: " << numIntegrationPoints_ << std::endl;
		//std::cout << "Integration points:" << std::endl;
		//std::cout << intPoints.ToString() << std::endl;
    
		//std::cout << "Mapping: " << std::endl;
		//std::map<int,UInt>::iterator mapit;
		//for(mapit=locPointIndices_.begin(); mapit!=locPointIndices_.end(); mapit++){
		//  std::cout << mapit->first << "; " << mapit->second << std::endl;
		//}
    
    
		PtrCoefFct matCoef = material_->GetTensorCoefFnc(matType_, tensorType_,
            Global::REAL, false);
    
		/*
     * POL_initialTensor_ is the small signal tensor from the mat file
     */
		//matCoef->GetTensor(POL_initialTensor_, lpm);
    
    needsInversion_ = false;
    POL_setWithFlux_ = false;
		// to calculate differential material properties, we need to know e0 / nu0
		if (material_->GetMaterialDatabaseName() == "Electrostatic") {
			rev_mat_fac_ = 8.854187817e-12; //eps0
			PDEName_ = "Electrostatic";
      
      // no longer used > CoefFunctionHystMat, CoefFunctionHystOuptut etc have their own
      // values for this; 
      // small signal tensor = permittivity
      // however, value is used for TestInversion function
      PtrCoefFct permittivity = material_->GetTensorCoefFnc(ELEC_PERMITTIVITY,tensorType_,
              Global::REAL, false);
      
      POL_eps_mu_SmallSignal_ = Matrix<Double>(dim_, dim_);
      permittivity->GetTensor(POL_eps_mu_SmallSignal_, lpm);
      
		} else if (material_->GetMaterialDatabaseName() == "Electromagnetics") {
			rev_mat_fac_ = 795774.7155; //nu0
			PDEName_ = "Electromagnetics";
      needsInversion_ = true;
      POL_setWithFlux_ = true; // only for scalar model with vecttor extension
      // no longer used > CoefFunctionHystMat, CoefFunctionHystOuptut etc have their own
      // values for this; 
      // however, value is used for TestInversion function
      // small signal tensor = permeability
      PtrCoefFct permeability = material_->GetTensorCoefFnc(MAG_PERMEABILITY,tensorType_,
              Global::REAL, false);
      
      POL_eps_mu_SmallSignal_ = Matrix<Double>(dim_, dim_);
      permeability->GetTensor(POL_eps_mu_SmallSignal_, lpm);
      
		} else {
			EXCEPTION("Currently only Electrostatics and Electromagnetics are supported");
		}
    
    //    std::cout << "StrainForm: " << strainForm << std::endl;
    //    std::cout << "COUPLED_useStrainForm_: " << COUPLED_useStrainForm_ << std::endl;
    //    
    //    
	}
  
  void CoefFunctionHyst::ReadAndSetWeights(BaseMaterial* const material, bool setForStrains){
    
    ParameterPreisachWeights paramSet = ParameterPreisachWeights(); 
    
    // use same offset as in XMLMaterialHandler.cc
    int enumOffset = 0;
    if(setForStrains){
      enumOffset = 100;
    }
    
    std::string usedHystModel;
    material->GetScalar(usedHystModel, MaterialType(HYST_MODEL+enumOffset));
    
    int numRows;
    material->GetScalar(numRows, MaterialType(PREISACH_WEIGHTS_DIM+enumOffset));
    paramSet.numRows_ = UInt(numRows);
    int weightTypeInt;
    material->GetScalar(weightTypeInt, MaterialType(PREISACH_WEIGHTS_TYPE+enumOffset));
    
    if(weightTypeInt == 0){
      paramSet.weightType_ = "Constant";
      material->GetScalar(paramSet.constWeight_, MaterialType(PREISACH_WEIGHTS_CONSTVALUE+enumOffset), Global::REAL);
    } else if(weightTypeInt == 1){
      paramSet.weightType_ = "muDat";
      int forHalfRange; // script for determining muDat parameter uses a Preisach model that is normalized
      // to [-0.5,0.5]; here we use [-1,1] as range instead
      // > parameter have to be adapted accordingly to fit different range
      // > if flag forHalfRange == true; this adapteion has to be done here; otherwise adaption is assumed to
      //    be done by user 
      material->GetScalar(forHalfRange, MaterialType(PREISACH_WEIGHTS_MUDAT_PARAMSFORHALFRANGE+enumOffset));
      
      material->GetScalar(paramSet.muDat_A_, MaterialType(PREISACH_WEIGHTS_MUDAT_A+enumOffset), Global::REAL);
      material->GetScalar(paramSet.muDat_sigma1_, MaterialType(PREISACH_WEIGHTS_MUDAT_SIGMA+enumOffset), Global::REAL);
      material->GetScalar(paramSet.muDat_h1_, MaterialType(PREISACH_WEIGHTS_MUDAT_H+enumOffset), Global::REAL);
      material->GetScalar(paramSet.muDat_eta_, MaterialType(PREISACH_WEIGHTS_MUDAT_ETA+enumOffset), Global::REAL);
      
      if(forHalfRange){
        /*
         * A_cfs = A_script/4
         * eta_cfs = eta_script
         * h_cfs = 2*h_script
         * sigma_cfs = sigma_script/2
         */
        paramSet.muDat_A_ /= 4;
        paramSet.muDat_sigma1_ /= 2;
        paramSet.muDat_h1_ *= 2;
      }

    } else if(weightTypeInt == 2){
      paramSet.weightType_ = "muDatExtended";
      int forHalfRange; // script for determining muDat parameter uses a Preisach model that is normalized
      // to [-0.5,0.5]; here we use [-1,1] as range instead
      // > parameter have to be adapted accordingly to fit different range
      // > if flag forHalfRange == true; this adapteion has to be done here; otherwise adaption is assumed to
      //    be done by user 
      material->GetScalar(forHalfRange, MaterialType(PREISACH_WEIGHTS_MUDAT_PARAMSFORHALFRANGE+enumOffset));
      
      //std::cout << "MuDatExtended" << std::endl;
      material->GetScalar(paramSet.muDat_A_, MaterialType(PREISACH_WEIGHTS_MUDAT_A+enumOffset), Global::REAL);
      material->GetScalar(paramSet.muDat_sigma1_, MaterialType(PREISACH_WEIGHTS_MUDAT_SIGMA+enumOffset), Global::REAL);
      material->GetScalar(paramSet.muDat_sigma2_, MaterialType(PREISACH_WEIGHTS_MUDAT_SIGMA2+enumOffset), Global::REAL);
      material->GetScalar(paramSet.muDat_h1_, MaterialType(PREISACH_WEIGHTS_MUDAT_H+enumOffset), Global::REAL);
      material->GetScalar(paramSet.muDat_h2_, MaterialType(PREISACH_WEIGHTS_MUDAT_H2+enumOffset), Global::REAL);
      material->GetScalar(paramSet.muDat_eta_, MaterialType(PREISACH_WEIGHTS_MUDAT_ETA+enumOffset), Global::REAL);
            
      if(forHalfRange){
        /*
         * A_cfs = A_script/4
         * eta_cfs = eta_script
         * h_cfs = 2*h_script
         * h2_cfs = 2*h2_script
         * sigma_cfs = sigma_script/2
         * sigma2_cfs = sigma2_script/2
         */
        paramSet.muDat_A_ /= 4;
        paramSet.muDat_sigma1_ /= 2;
        paramSet.muDat_h1_ *= 2;
        paramSet.muDat_sigma2_ /= 2;
        paramSet.muDat_h2_ *= 2;
      }
      
    } else if(weightTypeInt == 3){
      paramSet.weightType_ = "givenTensor";
      material->GetTensor(paramSet.weightTensor_, MaterialType(PREISACH_WEIGHTS_TENSOR+enumOffset), Global::REAL);
      //std::cout << "Found weights: " << paramSet.weightTensor_.ToString() << std::endl;
    } else {
      EXCEPTION("Weight type unknown");
    }
    
    // NEW: add additional anhysteretic curve to Preisach models
    material->GetScalar(paramSet.anhysteretic_a_ , MaterialType(PREISACH_WEIGHTS_ANHYST_A+enumOffset), Global::REAL);
    material->GetScalar(paramSet.anhysteretic_b_ , MaterialType(PREISACH_WEIGHTS_ANHYST_B+enumOffset), Global::REAL);
    material->GetScalar(paramSet.anhysteretic_c_ , MaterialType(PREISACH_WEIGHTS_ANHYST_C+enumOffset), Global::REAL);
    material->GetScalar(paramSet.anhysteretic_cInAtan_ , MaterialType(PREISACH_WEIGHTS_ANHYST_CINATAN+enumOffset));
    
//          std::cout << "anhystA: " << paramSet.anhysteretic_a_ << std::endl;
//      std::cout << "anhystB: " << paramSet.anhysteretic_b_ << std::endl;
//      std::cout << "anhystC: " << paramSet.anhysteretic_c_ << std::endl;
//    
    int anhystForHalfRange;
    // script for determining muDat parameter uses a Preisach model that is normalized
    // to [-0.5,0.5]; here we use [-1,1] as range instead
    // > parameter have to be adapted accordingly to fit different range
    // > if flag forHalfRange == true; this adapteion has to be done here; otherwise adaption is assumed to
    //    be done by user 
    material->GetScalar(anhystForHalfRange, MaterialType(PREISACH_WEIGHTS_ANHYST_PARAMSFORHALFRANGE+enumOffset));
    
    if(anhystForHalfRange){
      /*
       * a_cfs = 2*a_script (script assume anhystpart to have max ampl of 0.5 instead of 1)
       * b_cfs = b_script/2 (script multiplies b with e_norm in range [-0.5,0.5])
       * c_cfs = 2*c_script 
       * (if c is part of atan, i.e. a*atan(b*(e+c)) we have to scale by 2 as we have doubled range for e, too
       *  if c is not part of atan, i.e. a*atan(b*e) + c*e, c has to be scaled like a, i.e. by factor 2 
       */
      paramSet.anhysteretic_a_ *= 2;
      paramSet.anhysteretic_b_ /= 2;
      paramSet.anhysteretic_c_ *= 2;
    }

    int anhystOnlyInt = 0;
    material->GetScalar(anhystOnlyInt , MaterialType(PREISACH_WEIGHTS_ANHYST_ONLY+enumOffset));
    if(anhystOnlyInt == 1){
      paramSet.anhystOnly_ = true;
    } else {
      paramSet.anhystOnly_ = false;
    }
    
    int anhystPartCountsTowardsOutputSat = 0;
    material->GetScalar(anhystPartCountsTowardsOutputSat, MaterialType(PREISACH_WEIGHTS_ANHYSTCOUNTINGTOOUTPUTSAT+enumOffset));
    if(anhystPartCountsTowardsOutputSat == 1){
      paramSet.anhystCountingToOutputSat_ = true;
    } else {
      paramSet.anhystCountingToOutputSat_ = false;
    }
    
    
    // compute preisach weights (for scalar and vector sutor case first; vector mayergoyz gets special treatment later)
    paramSet.weightTensor_ = evaluatePreisachWeights(&paramSet);
    
    /*
     * Important note: ScalarPreisach model currently does not go over the triangle
     * -1 <= beta <= alpha <= 1 but over the whole square
     * -1 <= beta <= 1, -1 <= alpha <= 1 and divides the result by 2
     * > if we just fill up the triangle part of the Preisach plane with weights, we
     *    get wrong results
     * > easy fix: mirror weights along diagonal alpha = beta
     * > nicer fix: make sure that ScalarPreisach only uses the trianlge
     *            > requires more work as one needs to make sure that subtraction for
     *              partially overlapped elements still works
     *            > use easy fix initially
     */
    for(UInt i = 0; i < paramSet.numRows_; i++){
      // note: diagonal element does not need mirroring
      // >  k < i instead of k <= i
      for(UInt k = 0; k < i; k++){
        paramSet.weightTensor_[k][i] = paramSet.weightTensor_[i][k];
      }
    }
    
    // VERY IMPORTANT: make sure that \int_alpha,beta paramSet.PreisachWeights_ = 1 
    // NOTE: this is not true for mayergoyz adapted weights! for those, we have to calculate the
    // intergal over the scalar weights and divide by this term (if possible!)
    // NOTE2: for constant scalar weights, the transformed vector weights will all be 0.75* the scalar value (0f 0.5)
    // NOTE3: go over full range for alpha and beta, then divide by 2 later on
    Double intOverWeights = 0.0;
    for(UInt i = 0; i < paramSet.numRows_; i++){
      for(UInt k = 0; k < paramSet.numRows_; k++){
        intOverWeights += paramSet.weightTensor_[i][k];
      }
    }
    
    if(intOverWeights == 0){
      // special case: if all Preisach weights are 0, we solve only for the anhyst part
      paramSet.anhystOnly_ = true;
    } else {
      // problem: scaling to 1 does work for scalar model and vector model by sutor but
      // not for mayergoyz model; here the integral over all weights must be smaller
      // than one such that the sum over multiple directions adds up to correct values
      // (for const. weight for example, the weights for the mayergoyz model will be 3/4 
      // of that const. value > integral sum will only be 3/4, too)
      //
      // solution: determine scaling factor for standard weights, then apply this factor
      // to the transformed weights
      paramSet.anhystOnly_ = false;
      Double elemArea = 2.0/(paramSet.numRows_)*2.0/(paramSet.numRows_);
      intOverWeights *= elemArea/2; // factor of 1/2 as we have to integrate over the triangle area
      //std::cout << "intOverWeights = " << intOverWeights << std::endl;
      bool scalingRequired = true;
      
      if(usedHystModel == "vectorPreisach_Mayergoyz"){
        int isIsotropic = 1;
        material->GetScalar(isIsotropic, MaterialType(PREISACH_MAYERGOYZ_ISOTROPIC+enumOffset));
        
        if( (dim_ != 2) || (isIsotropic == 0)){
          EXCEPTION("Mayergoyz vector model currently only implemented for 2d isotropic materials");
        } else {
          POL_weightsAlreadyAdapted_ = 0;
          material->GetScalar(POL_weightsAlreadyAdapted_, MaterialType(PREISACH_WEIGHTS_FOR_MAYERGOYZ_VECTOR+enumOffset));
          if(POL_weightsAlreadyAdapted_ == 1){
            // here we have to assume that the weights given via the material file
            // already have been adapted (and scaled correctly)
            // > take weights as they are and do not further scale them!
            scalingRequired = false;
          } else {
            // compute transformation of weights
            paramSet.weightTensor_ = transformPreisachWeightsForIsotropicVectorCase(&paramSet);
            
            // transformPreisachWeightsForIsotropicVectorCase does only compute
            // the upper triangle > mirror
            for(UInt i = 0; i < paramSet.numRows_; i++){
              // note: diagonal element does not need mirroring
              // >  k < i instead of k <= i
              for(UInt k = 0; k < i; k++){
                paramSet.weightTensor_[k][i] = paramSet.weightTensor_[i][k];
              }
            }
            
            scalingRequired = true;
          } 
          // check if weights are symmetric w.r.t. alpha = -beta
          // this property derives from vector model (see "Mathematical Models of Hysteresis and their Applications" - Mayergoyz  p.164 eq(3.58) )
          // make weights symmetric
          //std::cout << "Weights before forcing symmetry: " << paramSet.PreisachWeights_.ToString() << std::endl;
          for(UInt i = 0; i < paramSet.numRows_; i++){
            for(UInt k = 0; k < paramSet.numRows_-i; k++){
              // iterate over triangle -1 < alpha < 1; -1 < beta < -alpha
              // in indices: 0 < i < numRows; 0 < k < numRows-i
              paramSet.weightTensor_[i][k] = paramSet.weightTensor_[paramSet.numRows_-k-1][paramSet.numRows_-i-1];
            }
          }
          //std::cout << "Weights after forcing symmetry: " << paramSet.PreisachWeights_.ToString() << std::endl;
        }
      }
      
      // store integral over weights; 
      // for our implementation we want to have the integral over all weights to be 1 which is why we scale
      // the weights; the implementation that derives the muDat parameter, however, allows the integral over
      // the weights to be different from 1, especially in case of anhysteretic parts being present; in the later
      // case, the saturation value outputSat refers to the sum of preisach in sat and anhyst part in sat;
      // to obtain the same curves we have two options:
      // a) allow int over weights to be different from 1, too
      // b) define saturation value of preisach alone to be outputSat - anhystPart in sat
      // > currently using option b
      // > this option can be turned of in mat.xml by setting anhystIncludedInSatValue to false
      //    by doing so, outputSat will be the saturation value of the preisach operator AND the scaling factor
      //    for the anhyst part
      paramSet.intOverWeights_ = intOverWeights;
//      std::cout << "intOverWeights_: " << intOverWeights << std::endl;
      if(scalingRequired){
        for(UInt i = 0; i < paramSet.numRows_; i++){
          for(UInt k = 0; k < paramSet.numRows_; k++){
            paramSet.weightTensor_[i][k] /= intOverWeights;
          }
        }
      }
      
      // get anhyst part at positive saturation to determine the saturation value of preisach operator alone
      // > the setting is done in operatorParams afterwards
      paramSet.anhystAtSat_normalized_ = Hysteresis::evalAnhystPart_normalized_atSaturation(paramSet.anhysteretic_a_,
              paramSet.anhysteretic_b_,paramSet.anhysteretic_c_,paramSet.anhysteretic_cInAtan_);
      
      

    }
    
    //    // check for negative weights
    //    for(UInt i = 0; i < paramSet.numRows_; i++){
    //      for(UInt k = 0; k < paramSet.numRows_; k++){
    //        if(paramSet.weightTensor_[i][k] < 0){
    //          std::cout << "i,k: " << i << ", " << k << std::endl;
    //          std::cout << "weight[i][k]: " << paramSet.weightTensor_[i][k] << std::endl;
    //        }
    //      }
    //    }
    
    //		/*
    //     * paramSet.freeFieldTensor_ is the tensor to be returned, when
    //     * RUN_tensorToReturn_ = 2
    //     * and the tensor to be added to deltaMat if
    //     * RUN_tensorToAdd_ = 2
    //     * (only ONE needed per PDE)
    //     */
    //		paramSet.freeFieldTensor_ = Matrix<Double>(dim_, dim_);
    //		paramSet.freeFieldTensor_.Init();
    //		for (UInt i = 0; i < dim_; i++) {
    //			paramSet.freeFieldTensor_[i][i] = rev_mat_fac_;
    //		}
    if(setForStrains){
      STRAIN_weightParams_ = paramSet;
    } else {
      POL_weightParams_ = paramSet;
    }
  }
  
  void CoefFunctionHyst::ReadAndSetParamsForHystOperator(BaseMaterial* const material, bool setForStrains){
//    std::cout << "ReadAndSetParamsForHystOperator" << std::endl;
    ParameterPreisachOperators paramSet = ParameterPreisachOperators(); 
    
    // use same offset as in XMLMaterialHandler.cc
    int enumOffset = 0;
    if(setForStrains){
      enumOffset = 100;
    }
    
    std::string dimTypeStr;
		material->GetScalar(dimTypeStr, MaterialType(PREISACH_DIM+enumOffset));
    
    if (dimTypeStr == "SCALAR") {
			paramSet.methodType_ = 0;
		} else if (dimTypeStr == "VECTOR") {
			paramSet.methodType_ = 1;
		}
    
    material->GetScalar(paramSet.methodName_, MaterialType(HYST_MODEL+enumOffset));
//    std::cout << "paramSet.methodName_: " << paramSet.methodName_ << std::endl;
    /*
     * Get model specific paramter
     */
    // initialize parameter for vector and scalar model just to make sure that
    // no random values from heap are used (even if some parameter are only needed
    // in scalar case and some only in vector case)
    // SCALAR ONLY
    paramSet.fixDirection_ = Vector<Double>(dim_);
    paramSet.fixDirection_.Init();
    paramSet.fixDirection_[0] = 1.0;
    
    // VECTOR ONLY (SUTOR/Mayergoyz isotropic)
    paramSet.rotResistance_ = 1.0;
    paramSet.angularDistance_ = 0.0;
    paramSet.angularResolution_ = 1e-4;
    paramSet.amplitudeResolution_ = 1e-17;
    paramSet.angularClipping_ = 0;
    paramSet.evalVersion_ = 2;   
    paramSet.isClassical_ = false;
    paramSet.scaleUpToSaturation_ = true;
    paramSet.printOut_ = false;
    paramSet.bmpResolution_ = 1000;
    paramSet.numDirections_ = 11;
    paramSet.outputClipping_ = 0;
    
    // USED FOR BOTH MODELS
    paramSet.fieldsAlignedAboveSat_ = true;
		paramSet.hystOutputRestrictedToSat_ = true;
    paramSet.hasInverseModel_ = true;
    
    material->GetScalar(paramSet.inputSat_, MaterialType(X_SATURATION+enumOffset), Global::REAL);
    material->GetScalar(paramSet.outputSat_, MaterialType(Y_SATURATION+enumOffset), Global::REAL);
    
    if(paramSet.methodName_ == "scalarPreisach"){
			//get direction
      // NEW: allow arbitrary direction
      paramSet.fixDirection_.Init();
      
			material->GetScalar(paramSet.fixDirection_[0], MaterialType(P_DIRECTION_X+enumOffset), Global::REAL);
			material->GetScalar(paramSet.fixDirection_[1], MaterialType(P_DIRECTION_Y+enumOffset), Global::REAL);
			if (dim_ == 3) {
				material->GetScalar(paramSet.fixDirection_[2], MaterialType(P_DIRECTION_Z+enumOffset), Global::REAL);
			}
      
      if(paramSet.fixDirection_.NormL2() == 0){
        WARN("Zero direction specified; taking default = x-direction");
        paramSet.fixDirection_[0] = 1.0;
      } else {
        paramSet.fixDirection_[0]/=paramSet.fixDirection_.NormL2();
        paramSet.fixDirection_[1]/=paramSet.fixDirection_.NormL2();
        if (dim_ == 3) {
          paramSet.fixDirection_[2]/=paramSet.fixDirection_.NormL2();
        }
      }
    } else if(paramSet.methodName_ == "vectorPreisach_Sutor"){
      paramSet.hasInverseModel_ = false;
      material->GetScalar(paramSet.rotResistance_, MaterialType(ROT_RESISTANCE+enumOffset), Global::REAL);
			material->GetScalar(paramSet.angularDistance_, MaterialType(ANG_DISTANCE+enumOffset), Global::REAL);
			material->GetScalar(paramSet.angularResolution_, MaterialType(ANG_RESOLUTION+enumOffset), Global::REAL);
			material->GetScalar(paramSet.angularClipping_, MaterialType(ANG_CLIPPING+enumOffset), Global::REAL);
			material->GetScalar(paramSet.amplitudeResolution_, MaterialType(AMP_RESOLUTION+enumOffset), Global::REAL);
      
			int printOut;
			int bmpResolution;
      if(!setForStrains){
        material->GetScalar(printOut, MaterialType(PRINT_PREISACH+enumOffset));
        material->GetScalar(bmpResolution, MaterialType(PRINT_PREISACH_RESOLUTION+enumOffset));
        /*
         * if printOut > 0 -> activate output of overlaid switching and rotation state; output every printOut timestep
         * bmpResolution -> number of pixels (std = 1000)
         */
        paramSet.printOut_ = (UInt) printOut;
        paramSet.bmpResolution_ = (UInt) bmpResolution;
      }
      
			int eval;
			material->GetScalar(eval, MaterialType(EVAL_VERSION+enumOffset));
      
			paramSet.evalVersion_ = (UInt) eval;
      
//      std::cout << "blub" << std::endl;
      
      int scaleToSat = 1;
      material->GetScalar(scaleToSat, MaterialType(SCALETOSAT+enumOffset));
      
//      std::cout << "Scale to Sat: " << scaleToSat << std::endl;
      
      if(scaleToSat == 1){
        paramSet.scaleUpToSaturation_ = true;
      } else {
        paramSet.scaleUpToSaturation_ = false;
      }
      
    } else if(paramSet.methodName_ == "vectorPreisach_Mayergoyz"){
      int isIsotropic = 1;
      material->GetScalar(isIsotropic, MaterialType(PREISACH_MAYERGOYZ_ISOTROPIC+enumOffset));
      if( (dim_ != 2) || (isIsotropic == 0)){
        EXCEPTION("Mayergoyz vector model currently only implemented for 2d isotropic materials");
      }
      paramSet.fieldsAlignedAboveSat_ = false;
      paramSet.hasInverseModel_ = false;
      paramSet.isIsotropic_ = true;
      
      paramSet.outputClipping_ = 0;
      material->GetScalar(paramSet.outputClipping_, MaterialType(PREISACH_MAYERGOYZ_CLIPOUTPUT+enumOffset));
      
			if(paramSet.outputClipping_ == 0){
				paramSet.hystOutputRestrictedToSat_ = false;
			}
			
      paramSet.numDirections_ = 0;
      material->GetScalar(paramSet.numDirections_, MaterialType(PREISACH_MAYERGOYZ_NUM_DIR+enumOffset));
      
    } else {
      std::stringstream exceptionMSG;
      exceptionMSG << paramSet.methodName_ << " is not available as hysteresis model";
      EXCEPTION(exceptionMSG.str());
    }
    
    if(setForStrains){
      STRAIN_operatorParams_ = paramSet;
    } else {
      POL_operatorParams_ = paramSet;
    }
//    std::cout << "POL_operatorParams_.methodName_: " << POL_operatorParams_.methodName_ << std::endl;  
  }
  
  
  void CoefFunctionHyst::ReadInMaterial(BaseMaterial* const material){
    
    bool forStrain = false;
    
    /*
     * Get weights first, then operator parameters
     */
    ReadAndSetWeights(material, forStrain);
    ReadAndSetParamsForHystOperator(material, forStrain);
    
    if(POL_weightParams_.anhystCountingToOutputSat_){
      POL_operatorParams_.preisachSat_ = POL_operatorParams_.outputSat_*(1.0 - POL_weightParams_.anhystAtSat_normalized_);
    } else {
      POL_operatorParams_.preisachSat_ = POL_operatorParams_.outputSat_;
    }

//    std::cout << "anhystCountingToOutputSat_: " << POL_weightParams_.anhystCountingToOutputSat_ << std::endl;
//    std::cout << "anhystAtSat_normalized_: " << POL_weightParams_.anhystAtSat_normalized_ << std::endl;
//    std::cout << "preisachSat_: " << POL_operatorParams_.preisachSat_ << std::endl;
//    std::cout << "outputSat_: " << POL_operatorParams_.outputSat_ << std::endl;
    
    /*
     * Set inversion for vector models
     * > only needed for polarization 
     */  
    inversionSet_ = false;
    if ( (POL_operatorParams_.methodName_ == "vectorPreisach_Sutor") || (POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") ){
      // inversion via LM > set parameter
      if (material_->GetMaterialDatabaseName() == "Electromagnetics") {
        inversionSet_ = true;
        Integer invMat, invMethod, maxNumReg, maxNumLS;
        
        material->GetScalar(invMat, MAX_NUM_IT_HYST_INV);
        LM_inversion_.maxNumIts = (UInt) invMat;
        
        material->GetScalar(invMethod, VEC_HYST_INV_METHOD);
        LM_inversion_.inversionMethod = (UInt) invMethod;
        
        material->GetScalar(LM_inversion_.tolH, RES_TOL_H_HYST_INV, Global::REAL);
        material->GetScalar(LM_inversion_.tolB, RES_TOL_B_HYST_INV, Global::REAL);
        
        material->GetScalar(maxNumReg, MAX_NUM_REG_IT_HYST_INV);
        LM_inversion_.maxNumRegIts = (UInt) maxNumReg;
        material->GetScalar(LM_inversion_.alphaRegStart, ALPHA_REG_HYST_INV, Global::REAL);
        material->GetScalar(LM_inversion_.alphaRegMin, ALPHA_REG_MIN_HYST_INV, Global::REAL);
        material->GetScalar(LM_inversion_.alphaRegMax, ALPHA_REG_MAX_HYST_INV, Global::REAL);  

        material->GetScalar(LM_inversion_.trustLow, TRUST_LOW_HYST_INV, Global::REAL);
        material->GetScalar(LM_inversion_.trustMid, TRUST_MID_HYST_INV, Global::REAL);
        material->GetScalar(LM_inversion_.trustHigh, TRUST_HIGH_HYST_INV, Global::REAL);  
        
        material->GetScalar(maxNumLS, MAX_NUM_LS_IT_HYST_INV);
        LM_inversion_.maxNumLSIts = (UInt) maxNumLS;
        material->GetScalar(LM_inversion_.alphaLSMin, ALPHA_LS_MIN_HYST_INV, Global::REAL);
        material->GetScalar(LM_inversion_.alphaLSMax, ALPHA_LS_MAX_HYST_INV, Global::REAL);  
        
        material->GetScalar(LM_inversion_.jacImplementation, JAC_IMPLEMENTATION_HYST_INV, Global::REAL);
        material->GetScalar(LM_inversion_.jacRes, JAC_RESOLUTION_HYST_INV, Global::REAL);
        int stopLSAtLocalMin = 0;
        material->GetScalar(stopLSAtLocalMin, STOP_INV_LS_AT_LOCAL_MIN); 
        if(stopLSAtLocalMin == 1){
          LM_inversion_.stopLineSearchAtLocalMin = true;
        } else {
          LM_inversion_.stopLineSearchAtLocalMin = false;
        }
        
        material->GetScalar(LM_inversion_.projLM_mu, HYST_INV_PROJLM_MU, Global::REAL);
        material->GetScalar(LM_inversion_.projLM_rho, HYST_INV_PROJLM_RHO, Global::REAL);
        material->GetScalar(LM_inversion_.projLM_beta, HYST_INV_PROJLM_BETA, Global::REAL);
        material->GetScalar(LM_inversion_.projLM_sigma, HYST_INV_PROJLM_SIGMA, Global::REAL);
        material->GetScalar(LM_inversion_.projLM_gamma, HYST_INV_PROJLM_GAMMA, Global::REAL);
        material->GetScalar(LM_inversion_.projLM_tau, HYST_INV_PROJLM_TAU, Global::REAL);
        material->GetScalar(LM_inversion_.projLM_c, HYST_INV_PROJLM_C, Global::REAL);
        material->GetScalar(LM_inversion_.projLM_p, HYST_INV_PROJLM_P, Global::REAL);
        
      } 
    } else if(POL_operatorParams_.methodName_ == "scalarPreisach"){
      inversionSet_ = true; // scalar model needs no additional parameter for inversion, it is always ready
    }
    
    /*
     * Read in initial input > has to be done AFTER operator for polarization
     * has been read in; initial state only given for polarization!
     * > input however will be applied to strain operator, too!
     */
    POL_initialInput_ = Vector<Double>(dim_);
    
    material->GetScalar(POL_initialInput_[0], INITIAL_STATE_X, Global::REAL);
    material->GetScalar(POL_initialInput_[1], INITIAL_STATE_Y, Global::REAL);
    if (dim_ == 3) {
      material->GetScalar(POL_initialInput_[2], INITIAL_STATE_Z, Global::REAL);
    }
    
//    std::cout << "Initial Input: " << POL_initialInput_.ToString() << std::endl;
    
    bool scaleBySaturation = false;
    int scaleBySaturationInt = 0;
    material->GetScalar(scaleBySaturationInt, PREISACH_SCALEINITIALSTATE);
    if(scaleBySaturationInt != 0){
      scaleBySaturation = true;
    }
    
    POL_prescribeInitialOutput_ = false;
    int prescribeOutputInt = 0;
    material->GetScalar(prescribeOutputInt, PREISACH_PRESCRIBEOUTPUT);
    if(prescribeOutputInt != 0){
      POL_prescribeInitialOutput_ = true;
    }
    
    if(POL_prescribeInitialOutput_ && scaleBySaturation){
      POL_initialInput_.ScalarMult(POL_operatorParams_.outputSat_);
    } else if(!POL_prescribeInitialOutput_ && scaleBySaturation){
      POL_initialInput_.ScalarMult(POL_operatorParams_.inputSat_);
    }
    
//    std::cout << "Initial Input after scaling: " << POL_initialInput_.ToString() << std::endl;
    
    POL_initial_.inputVector = POL_initialInput_;
    POL_initial_.prescribeOutput = POL_prescribeInitialOutput_;
    if(POL_initialInput_.NormL2() > 1e-16){
      POL_initial_.useInitialInput = true;
    } else {
      POL_initial_.useInitialInput = false;
    }
    POL_initial_.scaleBySaturation = scaleBySaturation;
    
    
    int isCoupled = 0;
    int reusePolarizationForStrains = 1;
    material->GetScalar(isCoupled, HYST_COUPLING_DEFINED);
    
    CouplingParams_  = ParameterIrrStrainsAndCoupling(); 
    
    if(isCoupled == 0){
      CouplingParams_.ownHystOperator_ = false;
      CouplingParams_.couplingDefined_inMatFile_ = false;
      CouplingParams_.useStrainForm_ = false;
      
    } else  if(isCoupled == 1){
      CouplingParams_.couplingDefined_inMatFile_ = true;
      
      /*
       * Get additional operators for coupling case
       */
      int strainForm;
      material->GetScalar(strainForm, HYST_STRAIN_FORM);
      
      /*
       * -1 : not coupled at all
       *  0 : coupled e-form/h-form (piezo/magstrict)
       *  1 : coupled d-form (piezo)
       *  2 : coupled g-form (magstrict)
       */
      CouplingParams_.strainForm_ = strainForm;
      if(strainForm > 0){
        CouplingParams_.useStrainForm_ = true;
      } else {
        CouplingParams_.useStrainForm_ = false;
      }
      
      /*
       * Input for piezoelectric / magnetostrictive setups
       */  
      material->GetScalar(CouplingParams_.ci_size_, HYST_IRRSTRAIN_CI_SIZE);
      material->GetTensor(CouplingParams_.ci_, HYST_IRRSTRAIN_CI, Global::REAL);
      material->GetScalar(CouplingParams_.c1_, HYST_IRRSTRAIN_C1, Global::REAL);
      material->GetScalar(CouplingParams_.c2_, HYST_IRRSTRAIN_C2, Global::REAL);
      material->GetScalar(CouplingParams_.c3_, HYST_IRRSTRAIN_C3, Global::REAL);
      material->GetScalar(CouplingParams_.d0_, HYST_IRRSTRAIN_D0, Global::REAL);
      material->GetScalar(CouplingParams_.d1_, HYST_IRRSTRAIN_D1, Global::REAL);
      material->GetScalar(CouplingParams_.irrStrainForm_, HYST_IRRSTRAINS);
      material->GetScalar(CouplingParams_.paramsForHalfRange_, HYST_IRRSTRAIN_PARAMSFORHALFRANGE);
      material->GetScalar(CouplingParams_.scaleTosSat_, HYST_IRRSTRAIN_SCALETOSAT);
      material->GetScalar(CouplingParams_.sSat_, S_SATURATION, Global::REAL);
      
      material->GetScalar(reusePolarizationForStrains, IRRSTRAIN_REUSE_P);
      
      if(reusePolarizationForStrains == 0){
        CouplingParams_.ownHystOperator_ = true;
        /*
         * Read in separate hyst operator for strains
         */
        forStrain = true;
        
        /*
         * Get weights first, then operator parameters
         */
        ReadAndSetWeights(material, forStrain);
        ReadAndSetParamsForHystOperator(material, forStrain);
        
        if(STRAIN_weightParams_.anhystCountingToOutputSat_){
          STRAIN_operatorParams_.preisachSat_ = STRAIN_operatorParams_.outputSat_*(1.0 - STRAIN_weightParams_.anhystAtSat_normalized_);
        } else {
          STRAIN_operatorParams_.preisachSat_ = STRAIN_operatorParams_.outputSat_;
        }
        
      } else {
        CouplingParams_.ownHystOperator_ = false;
      }
      
    } else {
      EXCEPTION("isCoupled should be 0 or 1")
    }
    
    //    EXCEPTION("Stop here");
  }
  
  
	CoefFunctionHyst::~CoefFunctionHyst() {
    
		delete timer_;
    
    if(storageInitialized_){
      delete[] E_B_;
      delete[] P_J_;
      delete[] E_H_;
      
      delete[] E_B_lastIt_;
      delete[] P_J_lastIt_;
      delete[] E_H_lastIt_;
      
      delete[] E_B_lastTS_;
      delete[] P_J_lastTS_;
      delete[] E_H_lastTS_;
      
      delete[] Si_;
      delete[] Si_lastTS_;
      delete[] Si_lastIt_;
      
      delete[] deltaMat_;
      delete[] deltaMatPrev_;
      delete[] deltaMatStrain_;
      delete[] deltaMatStrainPrev_;
      delete[] matrixForInversion_;
      delete[] matrixForInversionInverted_;
      delete[] rotatedCouplingTensor_;
      
      delete[] requiresReeval_;
      delete[] Si_requiresReeval_;
      delete[] deltaMat_requiresReeval_;
      delete[] deltaMatStrain_requiresReeval_;
      delete[] rotatedCouplingTensor_requiresReeval_;
      
      delete[] takeEstimatedSlope_;
    }
	}
  
	void CoefFunctionHyst::InitStorage() {
		/*
     * this function sets up the actual storage vectors and the the hysteretic
     * objects
     * it has to be called AFTER the XML dependent parameter are known
     * as the size of the storage depends on the evaluation depth (standard,
     * extended, full) which is retrieved from the non-linear parameters (see
     * stdSolveStep::ReadNonLinData)
     */
    if(storageInitialized_){
      return;
    }
    
		// 1. determine the number of required storages
		if (XML_EvaluationDepth_ == 1) {
			// standard evaluation
			// > one hyst operator per element
			// > only one storage for each element
			numHystOperators_ = numElemSD_;
			numStorageEntries_ = numElemSD_;
		} else if (XML_EvaluationDepth_ == 2) {
			// extended evaluation
			// > one hyst operator per element
			// > each integration point (+ midpoint) needs one storage
			numHystOperators_ = numElemSD_;
			numStorageEntries_ = numElemSD_ * (numIntegrationPoints_ + 1);
		} else if (XML_EvaluationDepth_ == 3) {
			// full evaluation
			// > each integration point (+ midpoint) get a hyst operator
			// > each integration point (+ midpoint) needs one storage
			numHystOperators_ = numElemSD_ * (numIntegrationPoints_ + 1);
			numStorageEntries_ = numElemSD_ * (numIntegrationPoints_ + 1);
		} else {
			EXCEPTION("Evaluation depth < 1 or > 3 not allowed")
		}
    
		//std::cout << "XML_EvaluationDepth_: " << XML_EvaluationDepth_ << std::endl;
    //		std::cout << "numIntegrationPoints_: " << numIntegrationPoints_ << std::endl;
    //		std::cout << "numStorageEntries_: " << numStorageEntries_ << std::endl;
    //		std::cout << "numHystOperators_: " << numHystOperators_ << std::endl;
//    std::cout << "CoefFunction::CoefDimType(POL_operatorParams_.methodType_): " << CoefFunction::CoefDimType(POL_operatorParams_.methodType_) << std::endl;
		if (POL_operatorParams_.methodType_ == 0) {
      
      bool isVirgin = true;
      POL_useExtension_ = false;
      //      int useExtensionInt;
      //      material_->GetScalar(useExtensionInt, SCALPREISACH_USE_EXT);
      //      if(useExtensionInt == 1){
      //        POL_useExtension_ = true;
      //      }
      if(POL_useExtension_){
        EXCEPTION("No longer available");
        //        material_->GetScalar(POL_operatorParams_.rotResistance_, ROT_RESISTANCE, Global::REAL);
        //        material_->GetScalar(POL_operatorParams_.angularDistance_, ANG_DISTANCE, Global::REAL);
        //        
        //        hyst_ = new ExtendedPreisach(numHystOperators_, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_, POL_weightParams_.weightTensor_, 
        //                POL_operatorParams_.rotResistance_, POL_operatorParams_.angularDistance_, dim_, isVirgin, POL_weightParams_.anhysteretic_a_, POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_, POL_weightParams_.anhystOnly_);
        //        
        //        // set initial direction
        //        POL_initialInput_ = Vector<Double>(dim_);
        //        
        //        material_->GetScalar(POL_initialInput_[0], INITIAL_STATE_X, Global::REAL);
        //        material_->GetScalar(POL_initialInput_[1], INITIAL_STATE_Y, Global::REAL);
        //        if (dim_ == 3) {
        //          material_->GetScalar(POL_initialInput_[2], INITIAL_STATE_Z, Global::REAL);
        //        }
        //        
        //        Vector<Double> initialInputForRotstate = POL_initialInput_;
        //        // if no value is given for initial state, set it to dirP
        //        if (initialInputForRotstate.NormL2() < 1e-16) {
        //          initialInputForRotstate.Init();
        //          if(POL_setWithFlux_){
        //            // update rotation states with flux quantity 
        //            // POL_initialInput_ = (POL_operatorParams_.outputSat_*Identity + POL_operatorParams_.inputSat_*POL_eps_mu_SmallSignal_)*POL_operatorParams_.fixDirection_
        //            // > flux with value just at saturation
        //            // > but: saturation alone would not suffice! we need to consdier anhyst parts, too
        //            POL_eps_mu_SmallSignal_.Mult(POL_operatorParams_.fixDirection_,initialInputForRotstate);
        //            initialInputForRotstate.ScalarMult(POL_operatorParams_.inputSat_);
        //            Double anhystPart = hyst_->evalAnhystPart_normalized(1.0);
        //            initialInputForRotstate.Add((anhystPart+1.0)*POL_operatorParams_.outputSat_,POL_operatorParams_.fixDirection_);
        //          } else {
        //            // update rotation states with field intensity
        //            // > scaling with x saturation sufficient
        //            initialInputForRotstate.Add(POL_operatorParams_.inputSat_,POL_operatorParams_.fixDirection_);
        //          }
        //        }
        //        
        //        if(POL_setWithFlux_){
        //          for(UInt i = 0; i < numHystOperators_; i++){
        //            hyst_->UpdateRotationStateWithFluxDensity(initialInputForRotstate,POL_eps_mu_SmallSignal_,i);
        //            hyst_->EvaluateRotationState(i);
        //          } 
        //        } else {
        //          for(UInt i = 0; i < numHystOperators_; i++){
        //            hyst_->UpdateRotationStateWithFieldIntensity(POL_initialInput_,i);
        //            hyst_->EvaluateRotationState(i);
        //          }
        //        }
        
      } else {
        bool ignoreAnhyst = false;
        hyst_ = new Preisach(numHystOperators_, POL_operatorParams_, POL_weightParams_, isVirgin, ignoreAnhyst);
//                (numHystOperators_, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_, POL_weightParams_.weightTensor_, isVirgin, 
//                POL_weightParams_.anhysteretic_a_, POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
      }
      
      POL_operatorParams_.hasInverseModel_ = false;
      // used during testing of hyst operator so it is set here
      LM_inversion_.tolB = 1e-10;
      LM_inversion_.tolH = 1e-8; // criterion as used in Preisach.cc for inversion; note: here for actual pol; in Preisach.cc scaled by XSatuated as it calcs with normalized values
                    
		} else if (POL_operatorParams_.methodName_ == "vectorPreisach_Sutor") {
      POL_operatorParams_.hasInverseModel_ = false;
      
			// this flag is not used currently;
			// if we have an initial state, isVirgin should be false, but as already set,
			// it is not used at the momement
			bool isVirgin = true;
      
			if (POL_operatorParams_.evalVersion_ == 1) {
				POL_operatorParams_.isClassical_ = true; // original vector preisach model -> sutor2012
        
				hyst_ = new VectorPreisachSutor_ListApproach(numHystOperators_, POL_operatorParams_, POL_weightParams_, dim_, isVirgin);
//                (numHystOperators_, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_,
//                POL_weightParams_.weightTensor_, POL_operatorParams_.rotResistance_, dim_, isVirgin,
//                POL_operatorParams_.isClassical_, POL_operatorParams_.scaleUpToSaturation_,
//                POL_operatorParams_.angularDistance_,POL_operatorParams_.angularResolution_,
//                POL_weightParams_.anhysteretic_a_, POL_weightParams_.anhysteretic_b_, 
//                POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
			} else if (POL_operatorParams_.evalVersion_ == 2) {
				POL_operatorParams_.isClassical_ = false; // revised vector preisach model -> sutor2015
        
				hyst_ = new VectorPreisachSutor_ListApproach(numHystOperators_, POL_operatorParams_, POL_weightParams_, dim_, isVirgin);
//                (numHystOperators_, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_,
//                POL_weightParams_.weightTensor_, POL_operatorParams_.rotResistance_, dim_, isVirgin,
//                POL_operatorParams_.isClassical_, POL_operatorParams_.scaleUpToSaturation_,
//                POL_operatorParams_.angularDistance_,POL_operatorParams_.angularResolution_,
//                POL_weightParams_.anhysteretic_a_, POL_weightParams_.anhysteretic_b_, 
//                POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
			} else if (POL_operatorParams_.evalVersion_ == 10) {
				POL_operatorParams_.isClassical_ = true; // original vector preisach model -> sutor2015; matrix based implementation
        
				hyst_ = new VectorPreisachSutor_MatrixApproach(numHystOperators_, POL_operatorParams_, POL_weightParams_, dim_, isVirgin);
//                (numHystOperators_, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_,
//                POL_weightParams_.weightTensor_, POL_operatorParams_.rotResistance_, dim_, isVirgin,
//                POL_operatorParams_.isClassical_, POL_operatorParams_.scaleUpToSaturation_,
//                POL_operatorParams_.angularDistance_,POL_operatorParams_.angularResolution_,
//                POL_weightParams_.anhysteretic_a_, POL_weightParams_.anhysteretic_b_, 
//                POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
			} else if (POL_operatorParams_.evalVersion_ == 20) {
				POL_operatorParams_.isClassical_ = false; // revised vector preisach model -> sutor2015; matrix based implementation
        
				hyst_ = new VectorPreisachSutor_MatrixApproach(numHystOperators_, POL_operatorParams_, POL_weightParams_, dim_, isVirgin);
//                (numHystOperators_, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_,
//                POL_weightParams_.weightTensor_, POL_operatorParams_.rotResistance_, dim_, isVirgin,
//                POL_operatorParams_.isClassical_, POL_operatorParams_.scaleUpToSaturation_,
//                POL_operatorParams_.angularDistance_,POL_operatorParams_.angularResolution_,
//                POL_weightParams_.anhysteretic_a_, POL_weightParams_.anhysteretic_b_, 
//                POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
			} else {
				EXCEPTION("POL_operatorParams_.evalVersion_ has to be one of the following: \n "
                "1: classical vector model (sutor2012) \n"
                "2: revised vector model (sutor2015) [DEFAULT] \n"
                "10: classical vector model (sutor2012) - Matrix implementation, only for reference \n"
                "20: revised vector model (sutor2015) - Matrix implementation, only for reference \n")
			}
      
      if( (POL_operatorParams_.scaleUpToSaturation_ == false)&&(POL_operatorParams_.isClassical_ == false) ){
        POL_operatorParams_.fieldsAlignedAboveSat_ = false;
      } 
      
		}
    else if (POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") {
      // basically a scalar model in multiple directions
      // isotropic case: all scalar models are equal (same weights etc)
      // anisotropic case: each model different; choice of directions matters; weights are harder to obtain
      POL_operatorParams_.fieldsAlignedAboveSat_ = false;
      bool isVirgin = true;
      
      /*
       * IMPORTANT REMARK:
       *  > although the Mayergoyz model is based on the scalar models 
       *     we are not alloewed to directly apply the preisach parameter for
       *     the scalar case (i.e. the weights, the anhyst parameter and so on)
       *  > make sure that the passed parameter are already transformed correctly
       *      > see constructor above
       */
      hyst_ = new VectorPreisachMayergoyz(numHystOperators_, POL_operatorParams_, POL_weightParams_, dim_, isVirgin);

//      hyst_ = new VectorPreisachMayergoyz(numHystOperators_, POL_operatorParams_.numDirections_, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_, 
//              POL_weightParams_.weightTensor_,dim_,isVirgin,POL_weightParams_.anhysteretic_a_, POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,
//              POL_weightParams_.anhystOnly_,POL_operatorParams_.outputClipping_);
//      
      POL_operatorParams_.hasInverseModel_ = false;      
    } else {
      EXCEPTION("Unknown hyst model");
    }
    
    if ( (POL_operatorParams_.methodName_ == "vectorPreisach_Sutor") || (POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") ){
      // inversion via LM > set parameter
      if (material_->GetMaterialDatabaseName() == "Electromagnetics") {
        hyst_->SetParamsForInversion(LM_inversion_);
      }
    }
    
    if(POL_operatorParams_.angularClipping_ != 0){
      WARN("Angular clipping currently not used; parameter will be ignored");
    }
    
    if(CouplingParams_.ownHystOperator_ && CouplingParams_.couplingDefined_inMatFile_){
//      std::cout << "Use own hyst operator for strain" << std::endl;
      bool isVirgin = true;
      
      if (STRAIN_operatorParams_.methodType_ == 0) {
        bool ignoreAnhystPart = false;
        hystStrain_ = new Preisach(numHystOperators_, STRAIN_operatorParams_, STRAIN_weightParams_, isVirgin, ignoreAnhystPart);
//                (numHystOperators_, STRAIN_operatorParams_.inputSat_, STRAIN_operatorParams_.outputSat_, 
//                STRAIN_weightParams_.weightTensor_, isVirgin, 
//                STRAIN_weightParams_.anhysteretic_a_, STRAIN_weightParams_.anhysteretic_b_, STRAIN_weightParams_.anhysteretic_c_,
//                STRAIN_weightParams_.anhystOnly_);
        
      } else if (POL_operatorParams_.methodName_ == "vectorPreisach_Sutor") {
        
        if (STRAIN_operatorParams_.evalVersion_ == 1) {
          STRAIN_operatorParams_.isClassical_ = true; // original vector preisach model -> sutor2012
          
          hystStrain_ = new VectorPreisachSutor_ListApproach(numHystOperators_, STRAIN_operatorParams_, STRAIN_weightParams_, dim_, isVirgin);
//                  (numHystOperators_, 
//                  STRAIN_operatorParams_.inputSat_, STRAIN_operatorParams_.outputSat_,
//                  STRAIN_weightParams_.weightTensor_, STRAIN_operatorParams_.rotResistance_, dim_, isVirgin,
//                  STRAIN_operatorParams_.isClassical_, STRAIN_operatorParams_.scaleUpToSaturation_,
//                  STRAIN_operatorParams_.angularDistance_,
//                  STRAIN_operatorParams_.angularResolution_,
//                  STRAIN_weightParams_.anhysteretic_a_, STRAIN_weightParams_.anhysteretic_b_, 
//                  STRAIN_weightParams_.anhysteretic_c_,STRAIN_weightParams_.anhystOnly_);
        } else if (STRAIN_operatorParams_.evalVersion_ == 2) {
          STRAIN_operatorParams_.isClassical_ = false; // revised vector preisach model -> sutor2015
          
          hystStrain_ = new VectorPreisachSutor_ListApproach(numHystOperators_, STRAIN_operatorParams_, STRAIN_weightParams_, dim_, isVirgin);
//                  (numHystOperators_, 
//                  STRAIN_operatorParams_.inputSat_, STRAIN_operatorParams_.outputSat_,
//                  STRAIN_weightParams_.weightTensor_, STRAIN_operatorParams_.rotResistance_, dim_, isVirgin,
//                  STRAIN_operatorParams_.isClassical_, STRAIN_operatorParams_.scaleUpToSaturation_,
//                  STRAIN_operatorParams_.angularDistance_,
//                  STRAIN_operatorParams_.angularResolution_,
//                  STRAIN_weightParams_.anhysteretic_a_, STRAIN_weightParams_.anhysteretic_b_, 
//                  STRAIN_weightParams_.anhysteretic_c_,STRAIN_weightParams_.anhystOnly_);
        } else if (STRAIN_operatorParams_.evalVersion_ == 10) {
          STRAIN_operatorParams_.isClassical_ = true; // original vector preisach model -> sutor2015; matrix based implementation
          
          hystStrain_ = new VectorPreisachSutor_MatrixApproach(numHystOperators_, STRAIN_operatorParams_, STRAIN_weightParams_, dim_, isVirgin);
//                  (numHystOperators_, 
//                  STRAIN_operatorParams_.inputSat_, STRAIN_operatorParams_.outputSat_,
//                  STRAIN_weightParams_.weightTensor_, STRAIN_operatorParams_.rotResistance_, dim_, isVirgin,
//                  STRAIN_operatorParams_.isClassical_, STRAIN_operatorParams_.scaleUpToSaturation_,
//                  STRAIN_operatorParams_.angularDistance_,
//                  STRAIN_operatorParams_.angularResolution_,
//                  STRAIN_weightParams_.anhysteretic_a_, STRAIN_weightParams_.anhysteretic_b_, 
//                  STRAIN_weightParams_.anhysteretic_c_,STRAIN_weightParams_.anhystOnly_);
        } else if (STRAIN_operatorParams_.evalVersion_ == 20) {
          STRAIN_operatorParams_.isClassical_ = false; // revised vector preisach model -> sutor2015; matrix based implementation
          
          hystStrain_ = new VectorPreisachSutor_MatrixApproach(numHystOperators_, STRAIN_operatorParams_, STRAIN_weightParams_, dim_, isVirgin);
//                  (numHystOperators_, 
//                  STRAIN_operatorParams_.inputSat_, STRAIN_operatorParams_.outputSat_,
//                  STRAIN_weightParams_.weightTensor_, STRAIN_operatorParams_.rotResistance_, dim_, isVirgin,
//                  STRAIN_operatorParams_.isClassical_, STRAIN_operatorParams_.scaleUpToSaturation_,
//                  STRAIN_operatorParams_.angularDistance_,
//                  STRAIN_operatorParams_.angularResolution_,
//                  STRAIN_weightParams_.anhysteretic_a_, STRAIN_weightParams_.anhysteretic_b_, 
//                  STRAIN_weightParams_.anhysteretic_c_,STRAIN_weightParams_.anhystOnly_);
        } else {
          EXCEPTION("STRAIN_operatorParams_.evalVersion_ has to be one of the following: \n "
                  "1: classical vector model (sutor2012) \n"
                  "2: revised vector model (sutor2015) [DEFAULT] \n"
                  "10: classical vector model (sutor2012) - Matrix implementation, only for reference \n"
                  "20: revised vector model (sutor2015) - Matrix implementation, only for reference \n")
        }
        if( (STRAIN_operatorParams_.scaleUpToSaturation_ == false)&&(STRAIN_operatorParams_.isClassical_ == false) ){
          STRAIN_operatorParams_.fieldsAlignedAboveSat_ = false;
        }  
      }
      else if (STRAIN_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") {
        // basically a scalar model in multiple directions
        // isotropic case: all scalar models are equal (same weights etc)
        // anisotropic case: each model different; choice of directions matters; weights are harder to obtain
        STRAIN_operatorParams_.fieldsAlignedAboveSat_ = false;
        
        /*
         * IMPORTANT REMARK:
         *  > although the Mayergoyz model is based on the scalar models 
         *     we are not alloewed to directly apply the preisach parameter for
         *     the scalar case (i.e. the weights, the anhyst parameter and so on)
         *  > make sure that the passed parameter are already transformed correctly
         *      > see constructor above
         */
        hystStrain_ = new VectorPreisachMayergoyz(numHystOperators_, STRAIN_operatorParams_, STRAIN_weightParams_, dim_, isVirgin);
//                (numHystOperators_, STRAIN_operatorParams_.numDirections_, 
//                STRAIN_operatorParams_.inputSat_, STRAIN_operatorParams_.outputSat_, 
//                STRAIN_weightParams_.weightTensor_,dim_,isVirgin,
//                STRAIN_weightParams_.anhysteretic_a_, STRAIN_weightParams_.anhysteretic_b_, STRAIN_weightParams_.anhysteretic_c_,
//                STRAIN_weightParams_.anhystOnly_,STRAIN_operatorParams_.outputClipping_);
        
      } else {
        EXCEPTION("Unknown hyst model");
      }
      
    } 
//    else {
//      std::cout << "Reuse pol to compute strain" << std::endl;
//    }
    
    /*
     * Initial state
     */
    // initial input that shall be feeded to the vectorPreisach operator
    
    Vector<Double> initial_E_H = Vector<Double>(dim_);
    initial_E_H.Init();
    Vector<Double> initial_E_B = Vector<Double>(dim_);
    initial_E_B.Init();
    Vector<Double> initial_P_J = Vector<Double>(dim_);
    initial_P_J.Init();
    Vector<Double> P_J_forStrains = Vector<Double>(dim_);
    P_J_forStrains.Init();
    
    if(POL_initial_.useInitialInput == true){
      //      WARN("Currently the treatment of initial states is not working properly! \n"
      //                "Depending on the selected evaluation approach, the initial state of \n"
      //                "the hysteresis operator will act as an excitation from the first iteration on or not. \n"
      //                "Furthermore, the initial state has to fit to the boundary condition \n"
      //                "(i.e. fluxParallel with initial state standing perpendicular on the boundary will not work).\n"
      //                "The only save usage is in context of the debugging fixPoint iteration \n"
      //                "(i.e. output of hysteresis operator does not couple back).")
      
      if(POL_initial_.prescribeOutput == false){
        // directly feed given input to hyst operator
        initial_E_H = POL_initialInput_;
      } else {       
        // output is prescribed (better said the total flux quantity is given, i.e.
        // irreversible + reversible + anhyst parts shall in sum be equal to the given vector
        // > feed this vector to inversion and retrieve actual input
        // > note: at the very start of the program, all hyst operators will be in virgin state
        // thus it suffices to compute the input just for one hyst operator
        initial_E_B = POL_initial_.inputVector;
        int successCode = 0;
        bool overwrite = false;
        UInt operatorIdx = 0;

        if (POL_operatorParams_.methodType_ == 0) {       
          Double muForInversion;
          Vector<Double> tmp = Vector<Double>(dim_);
          POL_eps_mu_SmallSignal_.Mult(POL_operatorParams_.fixDirection_,tmp);
          POL_operatorParams_.fixDirection_.Inner(tmp,muForInversion);
          
          Double scalInput;
          POL_operatorParams_.fixDirection_.Inner(initial_E_B,scalInput);
          Double scalOutput = hyst_->computeInputAndUpdate(scalInput,muForInversion, operatorIdx, overwrite, successCode);
          
          initial_E_H.Init();
          initial_E_H.Add(scalOutput,POL_operatorParams_.fixDirection_);
        } else {
          initial_E_H = hyst_->computeInput_vec(initial_E_B,operatorIdx,POL_eps_mu_SmallSignal_,
                  POL_operatorParams_.fieldsAlignedAboveSat_,POL_operatorParams_.hystOutputRestrictedToSat_,successCode);
        }
      }
      
      // now iterate over ALL hyst operator and evaluate them using initial_E_H_ to get initial_P_J_ AND
      // to set the initial memory accordingly
      bool overwrite = true;
      int successCode = 0;
      bool debugOut = false;
      for (UInt k = 0; k < numHystOperators_; k++) {
        /*
         * although all hyst operators should produce the same output, we have to
         * call the evaluation for each one of them to get its internal memory set
         */
        if (POL_operatorParams_.methodType_ == SCALAR) {    
          Double scalInput;

          POL_operatorParams_.fixDirection_.Inner(initial_E_H,scalInput);
          
          Double scalOutput = hyst_->computeValueAndUpdate(scalInput,k, overwrite,successCode);
          
          initial_P_J.Init();
          initial_P_J.Add(scalOutput,POL_operatorParams_.fixDirection_);
          
          if(CouplingParams_.ownHystOperator_){
            Double scalOutputForStrain = hystStrain_->computeValueAndUpdate(scalInput,k, overwrite,successCode);
            P_J_forStrains.Init();
            P_J_forStrains.Add(scalOutputForStrain,STRAIN_operatorParams_.fixDirection_);
          } else {
            P_J_forStrains = initial_P_J;
          }
          
        } else {
          initial_P_J = hyst_->computeValue_vec(initial_E_H, k, overwrite, debugOut, successCode);
          
          if(CouplingParams_.ownHystOperator_){
            P_J_forStrains = hystStrain_->computeValue_vec(initial_E_H, k, overwrite, debugOut, successCode);
          } else {
            P_J_forStrains = initial_P_J;
          }
        }
      }
    }
    
//    std::cout << "Initial output polarization: " << initial_P_J.ToString() << std::endl;
    
    Vector<Double> zeroVec = Vector<Double>(dim_);
    zeroVec.Init();
    Vector<Double> initial_Si = ComputeIrreversibleStrains(P_J_forStrains, initial_E_H, zeroVec);
    
//    std::cout << "Initial output strain: " << initial_Si.ToString() << std::endl;
		/*
     * finally initialize storage
     * NEW: use same storage for vector and scalar model
     * (even though the scalar model could use smaller storage, too)
     * NEW: create storage for each hystOperator or for each integration point
     * depending on evaluation depth
     */
    // axi case not supported!
    UInt strainDim = 3;
    if(dim_ == 3){
      strainDim = 6;
    }
    Vector<Double> zeroStrainVec = Vector<Double>(strainDim);
		zeroStrainVec.Init();
    
		E_B_ = new Vector<Double>[numStorageEntries_];
		P_J_ = new Vector<Double>[numStorageEntries_];
    P_J_forStrains_ = new Vector<Double>[numStorageEntries_];
		E_H_ = new Vector<Double>[numStorageEntries_];
    
		E_B_lastIt_ = new Vector<Double>[numStorageEntries_];
		P_J_lastIt_ = new Vector<Double>[numStorageEntries_];
		E_H_lastIt_ = new Vector<Double>[numStorageEntries_];
    
		E_B_lastTS_ = new Vector<Double>[numStorageEntries_];
		P_J_lastTS_ = new Vector<Double>[numStorageEntries_];
    P_J_forStrains_lastTS_ = new Vector<Double>[numStorageEntries_];
		E_H_lastTS_ = new Vector<Double>[numStorageEntries_];
    
    Si_ = new Vector<Double>[numStorageEntries_];
    Si_lastTS_ = new Vector<Double>[numStorageEntries_];
    Si_lastIt_ = new Vector<Double>[numStorageEntries_];
    
		deltaMat_ = new Matrix<Double>[numStorageEntries_];
    deltaMatPrev_ = new Matrix<Double>[numStorageEntries_];
    
    deltaMatStrain_ = new Matrix<Double>[numStorageEntries_];
    deltaMatStrainPrev_ = new Matrix<Double>[numStorageEntries_];
    
    matrixForInversion_ = new Matrix<Double>[numStorageEntries_];
    matrixForInversionInverted_ = new Matrix<Double>[numStorageEntries_];
    matrixForInversionWasSet_ = new bool[numStorageEntries_];
    
		requiresReeval_ = new bool[numStorageEntries_];
    Si_requiresReeval_ = new bool[numStorageEntries_];
    deltaMat_requiresReeval_ = new bool[numStorageEntries_];
    deltaMatStrain_requiresReeval_ = new bool[numStorageEntries_];
    
		rotatedCouplingTensor_ = new Matrix<Double>[numStorageEntries_];
		rotatedCouplingTensor_requiresReeval_ = new bool[numStorageEntries_];
    lastUsedTimeLevelForRotation_ = Vector<int>(numStorageEntries_);
    
    takeEstimatedSlope_ = new bool[numStorageEntries_];
		//TODO: coupling tensor einlesen und rotatedCouplingTensor entsprechend initialisieren
    
		Vector<Double> slope = Vector<Double>(dim_);
		slope.Init();
    
		for (UInt k = 0; k < numStorageEntries_; k++) {
			requiresReeval_[k] = true;
      Si_requiresReeval_[k] = true;
      deltaMat_requiresReeval_[k] = true;
      deltaMatStrain_requiresReeval_[k] = true;
      rotatedCouplingTensor_requiresReeval_[k] = true;
      lastUsedTimeLevelForRotation_[k] = -1;
      takeEstimatedSlope_[k] = false;
      // create with wrong size first; during first evaluation, set correct size
      rotatedCouplingTensor_[k] = Matrix<Double>(1,1);
      
      P_J_[k] = initial_P_J;
      P_J_forStrains_[k] = P_J_forStrains;
      E_B_[k] = initial_E_B;
      E_H_[k] = initial_E_H;
      
			E_B_lastIt_[k] = zeroVec;
			E_B_lastTS_[k] = zeroVec;
      
			E_H_lastIt_[k] = zeroVec;
			E_H_lastTS_[k] = zeroVec;
      
			P_J_lastIt_[k] = zeroVec;
			P_J_lastTS_[k] = zeroVec;
      P_J_forStrains_lastTS_[k] = zeroVec;
      
      Si_lastTS_[k] = zeroStrainVec;
      Si_lastIt_[k] = zeroStrainVec;
      Si_[k] = initial_Si;
      
      deltaMat_[k] = Matrix<Double>(dim_,dim_);
      deltaMat_[k].Init();// = POL_initialTensor_;
      
      deltaMatPrev_[k] = Matrix<Double>(dim_,dim_);
      deltaMatPrev_[k].Init();// = POL_initialTensor_;
      
      UInt numRows;
      if(dim_ == 3){
        // 6x3 tensor
        numRows = 6;
      } else {
        // 3x2 tensor
        numRows = 3;
      }
      deltaMatStrain_[k] = Matrix<Double>(numRows,dim_);
      deltaMatStrain_[k].Init();
      
      deltaMatStrainPrev_[k] = Matrix<Double>(numRows,dim_);
      deltaMatStrainPrev_[k].Init();
      
  		if (POL_initialInput_.NormL2() > 1e-16) {
        // TODO: compute first deltaMat or start without it
				// create a first deltaMatrix    
				//CreateDeltaMatrix(POL_initialInput_, POL_initialOutput_[k], deltaMat_[k], evalMethodName, k, intoSat, outofSat, satToSat, POL_initialInput_);
			}
      // for inversion we need mu in magnetic case (not nu!)
      // TODO: init matrixForInversion with 0 (and let CoefFundtionHystMatrial etc initialize with actual values
      //        > this has to be done in those additinal functions as this function here does not have the required tensors for the coupled case
      matrixForInversion_[k] = POL_eps_mu_SmallSignal_;//.Init();
      matrixForInversionInverted_[k] = Matrix<Double>(dim_,dim_);
      matrixForInversionInverted_[k].Init();
      matrixForInversionWasSet_[k] = false;
		}
    
    InitLPMMaps();
    storageInitialized_ = true;
	}
  
	
	void CoefFunctionHyst::AddAdditionalSDList(shared_ptr<EntityList> actSDList, RegionIdType volReg, bool isSurface){
		if(storageInitialized_ == true){
			EXCEPTION("Storage was already initialized! Cannot add further SDLists");
		}
		//std::cout << "Add additional SD list to HystOperator" << std::endl;
		SDList_List_.push_back(actSDList);
		volRegions_.insert(volReg);
    
		EntityIterator it = actSDList->GetIterator();
		
		UInt globalElNr;
		for (it.Begin(); !it.IsEnd(); it++) {
			globalElNr = it.GetElem()->elemNum;
			//std::cout << "Add element with #"<<globalElNr << " to list" << std::endl;
			globalElem2Local_[globalElNr] = numElemSD_;
			globalElemOnSurf_[globalElNr] = isSurface;
			numElemSD_++;
		}
		
		// add integration points to locPointIndices_
		
		// pick out the first element (even though any of them would do as they share the
		// same material) and extract midpoint
		it.Begin();
		const Elem * el = it.GetElem();
		
		//std::cout << "Elem Type: " << el->type << std::endl;
		LocPoint lp = Elem::shapes[el->type].midPointCoord;
		
		//std::cout << "Midpoint coords: " << Elem::shapes[el->type].midPointCoord.ToString() << std::endl;
		LocPointMapped lpm;
    bool updateShapeMap = false;
		shared_ptr<ElemShapeMap> esm = it.GetGrid()->GetElemShapeMap(el, updateShapeMap);
		//std::cout << "0" << std::endl;
		lpm.Set(lp, esm, 0.0);
		//std::cout << "1" << std::endl;
		// get number of integration points
		
		StdVector<LocPoint> intPoints;
		StdVector<Double> weights;
		// here we use the information about the integration method and order that we retrieved in the
		// constructor
		ptFeSpace_->GetIntScheme()->GetIntPoints(Elem::GetShapeType(el->type), IntegMethod_, IntegOrder_,
            intPoints, weights);
    
		UInt numIntegrationPointsLocal = intPoints.GetSize();
		if(numIntegrationPointsLocal > numIntegrationPoints_){
			// should not be the case as the number of integration points on surf elements usually is smaller than the
			// number of integration points on the volume
			//std::cout << "On the newly added element, more integration points are defined as on the already stored ones." << std::endl;
			numIntegrationPoints_ = numIntegrationPointsLocal;
		}
		
		// create a map that gives to each integration point + to midpoint an index
		// we cannot have a map of loc points directly (not sortable) but each point
		// has also a number which we can use for sorting; midpoint gets index -1
		// midpoint already inserted
		
		//locPointIndices_.insert(std::pair<int, UInt>(-1, 0));
		
		for (UInt i = 0; i < numIntegrationPointsLocal; i++) {
			//std::cout << intPoints[i].number << std::endl;
			//std::cout << "Add integration point pair (" << intPoints[i].number << ", " << i+1 << ")" << std::endl;
			locPointIndices_.insert(std::pair<int, UInt>(intPoints[i].number, i + 1));
		}
	}
	
  void CoefFunctionHyst::InitLPMMaps(){
    // This function has to be called after the storage was initialized (i.e. when all SDlists have been added)
    // > iterate over all elements and store each lpm (integration points + midpoint for full and extended, midpoint for standard)
    //    into a map with the storageIdx as key
    // > additionally store all midpoints in a separate map, again with storageIdx as key
    bool storeOnlyMidpoints = true;
    if(XML_EvaluationDepth_ > 1){
      storeOnlyMidpoints = false;
    }
    UInt storageIdx,operatorIdx;     
 		
		/*
     * 1. iterate over all elements
     */
    for (std::list<shared_ptr<EntityList> >::iterator listIt=SDList_List_.begin(); listIt != SDList_List_.end(); ++listIt) {
      shared_ptr<EntityList> SDList_ = *listIt;
			
      EntityIterator it = SDList_->GetIterator();
      LocPoint lp;
      LocPointMapped lpm;
      const Elem * el;
      shared_ptr<ElemShapeMap> esm;
      Vector<Double> zeroVec = Vector<Double>(dim_);
      zeroVec.Init();
      LocPointMapped actualLPM;
			
      for (it.Begin(); !it.IsEnd(); it++) {
        /*
         * 1.0 get element and shape map
         */
        el = it.GetElem();
        bool updateShapeMap = false;
        esm = it.GetGrid()->GetElemShapeMap(el, updateShapeMap);
				//std::cout << "ElemType: " << el->type << std::endl;
        
        if (storeOnlyMidpoints == false) {
          // Get lpm for integration points
          IntegOrder order;
          IntScheme::IntegMethod method;
          StdVector<LocPoint> intPoints;
          StdVector<Double> weights;
					
          ptFeSpace_->GetFe(it, method, order);
					
          ptFeSpace_->GetIntScheme()->GetIntPoints(Elem::GetShapeType(el->type), method, order,
                  intPoints, weights);
					
          // Loop over all integration points
          const UInt numIntPts = intPoints.GetSize();
          for (UInt i = 0; i < numIntPts; i++) {
						
            // Calculate for each integration point the LocPointMapped
            lpm.Set(intPoints[i], esm, weights[i]);
            
            
            //   preprocess LPM to get correct operator and storage indices
            PreprocessLPM(lpm, actualLPM, operatorIdx, storageIdx);
						//std::cout << "Integration point NR: " << i << std::endl;
						//std::cout << "On boundary? " << onBoundary << std::endl;
            
            // store ONLY the ORIGINAL lpm
            // even thought this one needs to be postprocessed again 
            // > reason for not storing the actual lpm: PreprocessLPM has to be called either way to get operatorIdx	   
            LOG_TRACE(coeffcthyst) << "Insert pair "<<storageIdx << " / " << lpm.lp.coord.ToString() << " into allLPMmap_";
            allLPMmap_[storageIdx]=lpm;
            
            // new: create mapping between operator and storage
            bool isMidpoint = false;
            MapStorageIndexToOperator(operatorIdx, storageIdx, isMidpoint);
            
          }
        }
        
        // get lpm for midpoint
        lp = Elem::shapes[el->type].midPointCoord;
        lpm.Set(lp, esm, 0.0);
        
				//std::cout << "Midpoint " << std::endl;
				//std::cout << "On boundary? " << onBoundary << std::endl;
        //      // preprocess LPM to get correct operator and storage indices
        bool forceMidpoint = true;
        PreprocessLPM(lpm, actualLPM, operatorIdx, storageIdx,forceMidpoint);
        LOG_TRACE(coeffcthyst) << "Insert pair "<<storageIdx << " / " << lpm.lp.coord.ToString() << " into allLPMmap_";
        allLPMmap_[storageIdx]=lpm;
        LOG_TRACE(coeffcthyst) << "Insert pair "<<storageIdx << " / " << lpm.lp.coord.ToString() << " into midpointLPMmap_";
        midpointLPMmap_[storageIdx]=lpm;
        
        // new: create mapping between operator and storage
        bool isMidpoint = true;
        // note: storageIdx for midpoint is always the one with the smallest number
        // compared to the other integration points
        // > as inner map is sorted by storageIdx, the first entry will always
        //   be the midpoint
        MapStorageIndexToOperator(operatorIdx, storageIdx, isMidpoint);
        
      }
    }
    
    LOG_TRACE(coeffcthyst) << "Mapping between operatorIndices and storageIndices";
    std::map< UInt, std::map< UInt, bool > >::iterator outerIt;
    std::map< UInt, bool >::iterator innerIt;
    for(outerIt = operatorToStorage_.begin(); outerIt != operatorToStorage_.end(); outerIt++){
      LOG_TRACE(coeffcthyst) << "Operator Idx - Storage Idx - isMidpoint";
      
      for(innerIt = (outerIt->second).begin(); innerIt != (outerIt->second).end(); innerIt++){
        LOG_TRACE(coeffcthyst) << outerIt->first << " - " << innerIt->first << " - " << innerIt->second;
      }
    }
    
  }
  
	bool CoefFunctionHyst::PreprocessLPM(const LocPointMapped& lpmInput,
          LocPointMapped& lpmOutput, UInt& operatorIdx, UInt& storageIdx, bool forceMidpoint) {
		/*
     *	Input:
     *		LocPointMapped lpmInput coming from numerical integration
     *
     *	Output:
     *		LocPointMapped lpmOutput specifying the position at which the hysteresis operator
     *		shall be evaluated
     *
     *		UInt operatorIdx reffering to the hysteresis operator that has to be evaluated
     *
     *		UInt storageIdx specifying the position in the storage arrays
     */
		const Elem * el = lpmInput.ptEl;
    //std::cout << "Preprocess LPM" << std::endl;
    //std::cout << "Input LPM on surface? " << lpmInput.isSurface << std::endl;
    //std::cout << "Corresponding element: " << el->ToString() << std::endl;
    //std::cout << "Coordinates of LP: " << lpmInput.lp.coord.ToString() << std::endl;
    
		/*
     * standard evaluation: midpointOnly ALWAYS used
     * extended evaluation: midpointOnly = true for output computation and
     *                                      during SetPreviousHystValues
     * full evaluation: midpointOnly = true only for output computation
     *
     * > use function EvaluateAtMidpointOnly() to determine the flag
     * 
     * for setprevioushystvalues, we might want to evaluate BOTH at integration points
     * and the midpoint; however the function EvaluateAtMidpointOnly would only give us
     * the actualy integrations points by then; by the additional input parameter
     * forceMidpoint, we can still retrieve the midpoint, though
     */
    if(forceMidpoint == false){
      forceMidpoint = EvaluateAtMidpointOnly();
    }
    
		if (forceMidpoint) {
      //std::cout << "Evaluate only at midpoint" << std::endl;
			/*
       * get solution at midpoint of element (elemSolution)
       */
			LocPoint lp = Elem::shapes[el->type].midPointCoord;
			shared_ptr<ElemShapeMap> esm = lpmInput.shapeMap;
      
			lpmOutput.Set(lp, esm, 0.0);
      
		} else {
      //std::cout << "Evaluate at integration point" << std::endl;
			/*
       * get solution at actual integration point
       */
			lpmOutput = lpmInput;
		}
    
		UInt idxElem = globalElem2Local_[el->elemNum];
    
		bool onBoundary = globalElemOnSurf_[el->elemNum];
    
    if(onBoundary){
      // set surface information for lpm
      lpmOutput.SetSurfInfo( volRegions_ );
    }
    
    //    if(onBoundary){
    //      std::cout << "LPM on boundary" << std::endl;
    //      
    //      std::cout << "lpmInput.isSurface: " << lpmInput.isSurface << std::endl;
    //      std::cout << "lpmInput.normal.ToString(): " << lpmInput.normal.ToString() << std::endl;
    //      std::cout << "lpmOutput.isSurface: " << lpmOutput.isSurface << std::endl;
    //      std::cout << "lpmOutput.normal.ToString(): " << lpmOutput.normal.ToString() << std::endl;
    //    }
    
		UInt idxPoint = locPointIndices_[lpmOutput.lp.number];
    //std::cout << "Global Indices: " << std::endl;
    //std::cout << "ElementI Number: " << el->elemNum << std::endl;
    //std::cout << "LPM Number: " << lpmOutput.lp.number << std::endl;		
		
		//std::cout << "Mapped indices: " << std::endl;
    //std::cout << "ElementIndex: " << idxElem << std::endl;
    //std::cout << "PointIndex: " << idxPoint << std::endl;
    
		// storage index = index where to store results, input, ...
		// operator index = index of hyst operator to evaluate
		if (XML_EvaluationDepth_ == 1) {
			// standard evaluation
			// storageIndex = operatorIndex = elementIndex
			storageIdx = idxElem;
			operatorIdx = idxElem;
		} else if (XML_EvaluationDepth_ == 2) {
			// extended evaluation
			// operatorIndex = elementIndex (only one operator per element)
			// storageIndex = combinedIndex = elementIndex*(numIntegrationPoints+1)+pointIndex
			operatorIdx = idxElem;
			storageIdx = idxElem * (numIntegrationPoints_ + 1) + idxPoint;
		} else {
			// full evaluation
			// one operator and one storage per integration point
			// operatorIndex = storageIndex = combinedIndex
			operatorIdx = idxElem * (numIntegrationPoints_ + 1) + idxPoint;
			storageIdx = operatorIdx;
		}
		
		return onBoundary;
	}
  
  Matrix<Double> CoefFunctionHyst::EstimateSlope(const LocPointMapped& Originallpm, bool includeStrains, bool useAbs, std::string implementationVersion,
          Double steppingLength, Double scaling){
    
    Matrix<Double> estimatedSlope;
    
    UInt numCols = dim_;
    UInt numRows = dim_;
    estimatedSlope = Matrix<Double>(numRows,numCols);
    
    UInt operatorIdx, storageIdx;
		LocPointMapped actualLPM;
    
		bool onBoundary = PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx);
	  // deltaMat can only be computed on volume elements, not on boundaries
		//bool onBoundary = false;
		if(onBoundary == true){
			EXCEPTION("Slope estimation only for volume elements.");
		}
    
    UInt currentTimelevel = 0;
    Vector<Double> E_B_current = RetrieveLPMSolution(actualLPM, storageIdx, currentTimelevel, onBoundary);
    
    // approximate slope via deltaMatrix approach
    // for the estimaten, step into 1,1,1 direction by some small value
    Double steppingDistance = steppingLength*POL_operatorParams_.inputSat_;
    Vector<Double> step = Vector<Double>(dim_);
    step.Init(steppingDistance);
    
    Vector<Double> E_B_stepping = Vector<Double>(dim_);
    E_B_stepping.Init();
    E_B_stepping.Add(1.0,E_B_current,1.0,step);
    
    Vector<Double> E_B_diff = Vector<Double>(dim_);
    E_B_diff.Init();
    E_B_diff.Add(1.0,E_B_stepping,-1.0,E_B_current);
    
    Vector<Double> P_J_current = GetPrecomputedOutputOfHysteresisOperator(Originallpm,currentTimelevel,false);
    bool forceMemoryLock = true;
    bool forceMemoryWrite = false;
    Vector<Double> P_J_stepping = CalcOutputOfHysteresisOperator(E_B_stepping,operatorIdx,storageIdx,forceMemoryLock,forceMemoryWrite,false);
    
    Vector<Double> P_J_diff = Vector<Double>(dim_);
    P_J_diff.Init();
    P_J_diff.Add(1.0,P_J_stepping,-1.0,P_J_current);
    
    P_J_diff.ScalarMult(scaling);
    
    estimatedSlope = CalcDeltaMat(E_B_diff, P_J_diff, useAbs,implementationVersion);
    
    LOG_TRACE(coeffcthyst) << "Estimated Slope / Delta Mat (storageIdx="<<storageIdx << "): " << estimatedSlope.ToString();
    
    if(includeStrains){
      Matrix<Double> estimatedSlopeStrains;
      if(dim_ == 3){
        // 6x3 tensor
        numRows = 6;
      } else {
        // 3x2 tensor
        numRows = 3;
      }
      estimatedSlopeStrains = Matrix<Double>(numRows,numCols);
      estimatedSlopeStrains.Init();
      // TODO:
      // 1. compute Irreversible strains from P_J_
      // 2. get deltaMat for strains
      // 3. multiply deltaMat for strains with coupling tensors
      // 4. add multiplied deltaMat to deltaMat from polarization
      
      WARN("Strain computation not yet added to slope estimation.");
    } 
    
    return estimatedSlope;
  }
  
  Matrix<Double> CoefFunctionHyst::CalcDeltaMat(Vector<Double> E_B_diff, Vector<Double> P_J_diff, bool useAbs,std::string implementationVersion, Double cuttingTol){
    //    std::cout << "CalcDeltaMat" << std::endl;
    Matrix<Double> deltaMat;
    UInt numCols = dim_;
    UInt numRows = dim_;
    
    Double absTol = 1e-14;//1e-16;
    Double relTol = 1e-14;//1e-16;
    
    deltaMat = Matrix<Double>(numRows,numCols);
    deltaMat.Init();
    
    Double E_B_diff_norm = E_B_diff.NormL2();
    
    if(E_B_diff_norm <= absTol){
      LOG_DBG(coeffcthystdeltamat) << "E_B_diff_norm <= " << absTol;
      LOG_DBG(coeffcthystdeltamat) << "> set deltaMat to 0";
      // variation of solution is very small
      // return zero matrix
      return deltaMat;
    }
    
    for(UInt i = 0; i < P_J_diff.GetSize(); i++){
      // if cuttingTol is negative (default), nothing happens
      // else: cut down nominato
      // note: s_diff is passed by value, so it is safe to edit vector directly
      LOG_DBG(coeffcthystdeltamat) << "P_J_diff["<<i<<"]: " << P_J_diff[i];
      LOG_DBG(coeffcthystdeltamat) << "cuttingTol: " << cuttingTol;
      LOG_DBG(coeffcthystdeltamat) << "E_B_diff["<<i<<"]: " << E_B_diff[i];
      if(abs(P_J_diff[i]) < cuttingTol){
        LOG_DBG(coeffcthystdeltamat) << "P_J_diff["<<i<<"] <= " << cuttingTol;
        LOG_DBG(coeffcthystdeltamat) << "> set P_J_diff["<<i<<"] to zero";
        //P_J_diff[i] = 0.0;
      }
    }
    
    //std::cout << "Difference Vector: " << P_J_diff.ToString() << std::endl;
    if(implementationVersion == "Division"){
      if(dim_ == 2){
        
        if(abs(E_B_diff[0])/E_B_diff_norm > relTol){
          // division can be used
          deltaMat[0][0] = P_J_diff[0]/E_B_diff[0];
        } else if(abs(E_B_diff[1])/E_B_diff_norm > relTol) {
          // distribute entry to non-diagonal entries
          deltaMat[0][1] = P_J_diff[0]/E_B_diff[1];
          deltaMat[1][0] = P_J_diff[0]/E_B_diff[1];
        } 
        //          else {
        //            // all entries are too small (relatively > return zero matrx
        //            // leave entry zero
        //          }
        
        if(abs(E_B_diff[1])/E_B_diff_norm > relTol){
          deltaMat[1][1] = P_J_diff[1]/E_B_diff[1];
        } else if(abs(E_B_diff[0])/E_B_diff_norm > relTol) {
          // distribute entry to non-diagonal entries
          deltaMat[0][1] = P_J_diff[1]/E_B_diff[0];
          deltaMat[1][0] = P_J_diff[1]/E_B_diff[0];
        } 
        //          else {
        //            // all entries are too small (relatively > return zero matrx
        //            // leave entry zero
        //          }
        
      } else {
        Double sum,diff;
        sum = E_B_diff[1]+E_B_diff[2];
        diff = E_B_diff[1]-E_B_diff[2];
        
        if(abs(E_B_diff[0])/E_B_diff_norm > relTol){
          // division can be used
          deltaMat[0][0] = P_J_diff[0]/E_B_diff[0];
        } else if(abs(sum)/E_B_diff_norm > relTol) {
          // distribute entry to non-diagonal entries
          deltaMat[0][1] = P_J_diff[0]/sum;
          deltaMat[1][0] = P_J_diff[0]/sum;
          deltaMat[0][2] = P_J_diff[0]/sum;
          deltaMat[2][0] = P_J_diff[0]/sum;
        } else if(abs(diff)/E_B_diff_norm > relTol) {
          // distribute entry to non-diagonal entries
          deltaMat[0][1] = P_J_diff[0]/diff;
          deltaMat[1][0] = P_J_diff[0]/diff;
          deltaMat[0][2] = -P_J_diff[0]/diff;
          deltaMat[2][0] = -P_J_diff[0]/diff;
        } 
        //          else {
        //            // all entries are too small (relatively > return zero matrx
        //            // leave entry zero
        //          }
        
        sum = E_B_diff[0]+E_B_diff[2];
        diff = E_B_diff[0]-E_B_diff[2];
        
        if(abs(E_B_diff[1])/E_B_diff_norm > relTol){
          // division can be used
          deltaMat[1][1] = P_J_diff[1]/E_B_diff[1];
        } else if(abs(sum)/E_B_diff_norm > relTol) {
          // distribute entry to non-diagonal entries
          deltaMat[0][1] = P_J_diff[1]/sum;
          deltaMat[1][0] = P_J_diff[1]/sum;
          deltaMat[1][2] = P_J_diff[1]/sum;
          deltaMat[2][1] = P_J_diff[1]/sum;
        } else if(abs(diff)/E_B_diff_norm > relTol) {
          // distribute entry to non-diagonal entries
          deltaMat[0][1] = P_J_diff[1]/diff;
          deltaMat[1][0] = P_J_diff[1]/diff;
          deltaMat[1][2] = -P_J_diff[1]/diff;
          deltaMat[2][1] = -P_J_diff[1]/diff;
        } 
        //          else {
        //            // all entries are too small (relatively > return zero matrx
        //            // leave entry zero
        //          }
        
        sum = E_B_diff[0]+E_B_diff[1];
        diff = E_B_diff[0]-E_B_diff[1];
        
        if(abs(E_B_diff[2])/E_B_diff_norm > relTol){
          // division can be used
          deltaMat[2][2] = P_J_diff[2]/E_B_diff[2];
        } else if(abs(sum)/E_B_diff_norm > relTol) {
          // distribute entry to non-diagonal entries
          deltaMat[0][2] = P_J_diff[2]/sum;
          deltaMat[2][0] = P_J_diff[2]/sum;
          deltaMat[1][2] = P_J_diff[2]/sum;
          deltaMat[2][1] = P_J_diff[2]/sum;
        } else if(abs(diff)/E_B_diff_norm > relTol) {
          // distribute entry to non-diagonal entries
          deltaMat[0][2] = P_J_diff[2]/diff;
          deltaMat[2][0] = P_J_diff[2]/diff;
          deltaMat[1][2] = -P_J_diff[2]/diff;
          deltaMat[2][1] = -P_J_diff[2]/diff;
        } 
        //          else {
        //            // all entries are too small (relatively > return zero matrx
        //            // leave entry zero
        //          }                   
      }
    } else {
      EXCEPTION("Delta mat version not implemented yet");
    }
    
    return deltaMat;
  };
  
  Matrix<Double> CoefFunctionHyst::CalcDeltaMatStrains(Vector<Double> E_B_diff, Vector<Double> S_diff, bool useAbs, std::string implementationVersion, Double cuttingTol){
    
    //    std::cout << "CalcDeltaMatStrains" << std::endl;
    Matrix<Double> deltaMat;
    UInt numCols = dim_;
    UInt numRows;
    
    Double absTol = 1e-16;
    Double relTol = 1e-16;
    
    if(dim_ == 3){
      // 6x3 tensor
      numRows = 6;
    } else {
      // 3x2 tensor
      numRows = 3;
    }
    deltaMat = Matrix<Double>(numRows,numCols);
    deltaMat.Init();
    
    Double E_B_diff_norm = E_B_diff.NormL2();
    
    if(E_B_diff_norm <= absTol){
      //      std::cout << "E_B_diff_norm <= absTol" << std::endl;
      // variation of solution is very small
      // return zero matrix
      return deltaMat;
    }
    
    for(UInt i = 0; i < S_diff.GetSize(); i++){
      // if cuttingTol is negative (default), nothing happens
      // else: cut down nominato
      // note: s_diff is passed by value, so it is safe to edit vector directly
      if(abs(S_diff[i]) < cuttingTol){
        S_diff[i] = 0.0;
      }
    }
    
    if(implementationVersion == "Division"){
      //      std::cout << "Division" << std::endl;
      if(dim_ == 2){
        
        if(abs(E_B_diff[0])/E_B_diff_norm > relTol){
          // division can be used
          deltaMat[0][0] = S_diff[0]/E_B_diff[0];
        } else if(abs(E_B_diff[1])/E_B_diff_norm > relTol) {
          // distribute entry to non-diagonal entries
          deltaMat[0][1] = S_diff[0]/E_B_diff[1];
          deltaMat[1][0] = S_diff[0]/E_B_diff[1];
        } 
        //          else {
        //            // all entries are too small (relatively > return zero matrx
        //            // leave entry zero
        //          }
        
        if(abs(E_B_diff[1])/E_B_diff_norm > relTol){
          deltaMat[1][1] = S_diff[1]/E_B_diff[1];
        } else if(abs(E_B_diff[0])/E_B_diff_norm > relTol) {
          // distribute entry to non-diagonal entries
          deltaMat[0][1] = S_diff[1]/E_B_diff[0];
          deltaMat[1][0] = S_diff[1]/E_B_diff[0];
        } 
        //          else {
        //            // all entries are too small (relatively > return zero matrx
        //            // leave entry zero
        //          }
        Double sum,diff;
        sum = E_B_diff[0]+E_B_diff[1];
        diff = E_B_diff[0]-E_B_diff[1];
        
        if(abs(sum)/E_B_diff_norm > relTol) {
          // distribute entry to non-diagonal entries
          deltaMat[2][0] = S_diff[1]/sum;
          deltaMat[2][1] = S_diff[1]/sum;
        } else if(abs(diff)/E_B_diff_norm > relTol) {
          // distribute entry to non-diagonal entries
          deltaMat[2][0] = S_diff[1]/diff;
          deltaMat[2][1] = -S_diff[1]/diff;
        } 
        
      } else {
        Double sum,diff;
        sum = E_B_diff[1]+E_B_diff[2];
        diff = E_B_diff[1]-E_B_diff[2];
        
        if(abs(E_B_diff[0])/E_B_diff_norm > relTol){
          // division can be used
          deltaMat[0][0] = S_diff[0]/E_B_diff[0];
        } else if(abs(sum)/E_B_diff_norm > relTol) {
          // distribute entry to non-diagonal entries
          deltaMat[0][1] = S_diff[0]/sum;
          deltaMat[1][0] = S_diff[0]/sum;
          deltaMat[0][2] = S_diff[0]/sum;
          deltaMat[2][0] = S_diff[0]/sum;
        } else if(abs(diff)/E_B_diff_norm > relTol) {
          // distribute entry to non-diagonal entries
          deltaMat[0][1] = S_diff[0]/diff;
          deltaMat[1][0] = S_diff[0]/diff;
          deltaMat[0][2] = -S_diff[0]/diff;
          deltaMat[2][0] = -S_diff[0]/diff;
        } 
        //          else {
        //            // all entries are too small (relatively > return zero matrx
        //            // leave entry zero
        //          }
        
        sum = E_B_diff[0]+E_B_diff[2];
        diff = E_B_diff[0]-E_B_diff[2];
        
        if(abs(E_B_diff[1])/E_B_diff_norm > relTol){
          // division can be used
          deltaMat[1][1] = S_diff[1]/E_B_diff[1];
        } else if(abs(sum)/E_B_diff_norm > relTol) {
          // distribute entry to non-diagonal entries
          deltaMat[0][1] = S_diff[1]/sum;
          deltaMat[1][0] = S_diff[1]/sum;
          deltaMat[1][2] = S_diff[1]/sum;
          deltaMat[2][1] = S_diff[1]/sum;
        } else if(abs(diff)/E_B_diff_norm > relTol) {
          // distribute entry to non-diagonal entries
          deltaMat[0][1] = S_diff[1]/diff;
          deltaMat[1][0] = S_diff[1]/diff;
          deltaMat[1][2] = -S_diff[1]/diff;
          deltaMat[2][1] = -S_diff[1]/diff;
        } 
        //          else {
        //            // all entries are too small (relatively > return zero matrx
        //            // leave entry zero
        //          }
        
        sum = E_B_diff[0]+E_B_diff[1];
        diff = E_B_diff[0]-E_B_diff[1];
        
        if(abs(E_B_diff[2])/E_B_diff_norm > relTol){
          // division can be used
          deltaMat[2][2] = S_diff[2]/E_B_diff[2];
        } else if(abs(sum)/E_B_diff_norm > relTol) {
          // distribute entry to non-diagonal entries
          deltaMat[0][2] = S_diff[2]/sum;
          deltaMat[2][0] = S_diff[2]/sum;
          deltaMat[1][2] = S_diff[2]/sum;
          deltaMat[2][1] = S_diff[2]/sum;
        } else if(abs(diff)/E_B_diff_norm > relTol) {
          // distribute entry to non-diagonal entries
          deltaMat[0][2] = S_diff[2]/diff;
          deltaMat[2][0] = S_diff[2]/diff;
          deltaMat[1][2] = -S_diff[2]/diff;
          deltaMat[2][1] = -S_diff[2]/diff;
        } 
        //          else {
        //            // all entries are too small (relatively > return zero matrx
        //            // leave entry zero
        //          }     
        
        sum = E_B_diff[1]+E_B_diff[2];
        diff = E_B_diff[1]-E_B_diff[2];
        
        if(abs(sum)/E_B_diff_norm > relTol) {
          // ds_zy/(de_y + de_z)
          deltaMat[3][1] = S_diff[3]/sum;
          deltaMat[3][2] = S_diff[3]/sum;
        } else if(abs(diff)/E_B_diff_norm > relTol) {
          // ds_zy/(de_y + de_z)
          deltaMat[3][1] = S_diff[3]/diff;
          deltaMat[3][2] = -S_diff[3]/diff;
        } else if(abs(E_B_diff[0])/E_B_diff_norm > relTol) {
          deltaMat[3][0] = S_diff[3]/E_B_diff[0];
        }
        
        sum = E_B_diff[0]+E_B_diff[2];
        diff = E_B_diff[0]-E_B_diff[2];
        
        if(abs(sum)/E_B_diff_norm > relTol) {
          // ds_zy/(de_y + de_z)
          deltaMat[4][0] = S_diff[4]/sum;
          deltaMat[4][2] = S_diff[4]/sum;
        } else if(abs(diff)/E_B_diff_norm > relTol) {
          // ds_zy/(de_y + de_z)
          deltaMat[4][0] = S_diff[4]/diff;
          deltaMat[4][2] = -S_diff[4]/diff;
        } else if(abs(E_B_diff[1])/E_B_diff_norm > relTol) {
          deltaMat[4][1] = S_diff[4]/E_B_diff[1];
        }
        
        sum = E_B_diff[0]+E_B_diff[1];
        diff = E_B_diff[0]-E_B_diff[1];
        
        if(abs(sum)/E_B_diff_norm > relTol) {
          // ds_zy/(de_y + de_z)
          deltaMat[5][0] = S_diff[5]/sum;
          deltaMat[5][1] = S_diff[5]/sum;
        } else if(abs(diff)/E_B_diff_norm > relTol) {
          // ds_zy/(de_y + de_z)
          deltaMat[5][0] = S_diff[5]/diff;
          deltaMat[5][1] = -S_diff[5]/diff;
        } else if(abs(E_B_diff[2])/E_B_diff_norm > relTol) {
          deltaMat[5][2] = S_diff[5]/E_B_diff[2];
        } 
        
      }
    } else {
      EXCEPTION("Delta mat version not implemented yet");
    }
    return deltaMat;
  };
  
  
  Matrix<Double> CoefFunctionHyst::GetDeltaMat(const LocPointMapped& Originallpm, int timelevel_new, int timelevel_old, bool useStrains, bool useAbs,
          std::string implementationVersion){
    
    
    LOG_DBG(coeffcthystdeltamat) << "GetDeltaMat";
    //LOG_DBG(coeffcthystdeltamat) << "Compute delta between timelevel " << timelevel_new << " and time level " << timelevel_old << "(-2 = deactivated; -1 = lastTS; 0 = current IT; 1 = lastIT)";
    
    if(tensorType_ == AXI){
      EXCEPTION("GetDeltaMat only implemented for 2d plane and full 3d setups");
    }
    
    Matrix<Double> deltaMat;
    UInt numCols = dim_;
    UInt numRows;
    
    if(timelevel_new == timelevel_old){
      WARN("DeltaMat was requested with timelevel_new = timelevel_old, i.e. delta would be 0; this should not happen"); 
    }
    
    UInt operatorIdx, storageIdx;
		LocPointMapped actualLPM;
    
		bool onBoundary = PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx);
	  // deltaMat can only be computed on volume elements, not on boundaries
		//bool onBoundary = false;
		if(onBoundary == true){
			EXCEPTION("DeltaMat not defined on boundary");
		}
		
    if(takeEstimatedSlope_[storageIdx]){
      takeEstimatedSlope_[storageIdx] = false;
      return EstimateSlope(Originallpm, ES_includeStrains_, ES_useAbs_, ES_implementationVersion_, ES_steppingLength_, ES_scaling_);
    }
    
    
    if(useStrains){
      //      std::cout << "GetDeltaMat - for strains" << std::endl;
      if(deltaMatStrain_requiresReeval_[storageIdx] == false){
        //        std::cout << "NO reeval" << std::endl;
        return deltaMatStrain_[storageIdx];
      }
      
      if(dim_ == 3){
        numRows = 6;
      } else {
        numRows = 3;
      }
      deltaMat = Matrix<Double>(numRows,numCols);
      
      /*
       * Compute matrix dM, such that
       * dS = dM*dE (or dB)
       * 
       * with dS = Strains(timeLevel1)-Strains(timeLevel2)
       *      dE = E(timeLevel1)-E(timeLevel2)
       */
      
      Vector<Double> E_B_new = RetrieveLPMSolution(actualLPM, storageIdx, timelevel_new, onBoundary);
      Vector<Double> E_B_old = RetrieveLPMSolution(actualLPM, storageIdx, timelevel_old, onBoundary);
      
      Vector<Double> E_B_diff = E_B_new;
      E_B_diff -= E_B_old;
      
      Vector<Double> S_new = GetIrreversibleStrains(Originallpm,timelevel_new);
      Vector<Double> S_old = GetIrreversibleStrains(Originallpm,timelevel_old);
      Vector<Double> S_diff = S_new;
      S_diff -= S_old;
      
      // cutting does not work so well either
      Double cuttingTol = 1e-5*S_old.NormL2();
      
      Matrix<Double> deltaTMP = CalcDeltaMatStrains(E_B_diff, S_diff, useAbs, implementationVersion,cuttingTol);
      
      // Blending with old deltaMat does not work as expected
      if( (includeOldDelta_ > 0)&&(includeOldDelta_ <= 1)){
        deltaMat.Init();
        deltaMat.Add(1.0-includeOldDelta_,deltaTMP);
        deltaMat.Add(includeOldDelta_,deltaMatStrainPrev_[storageIdx]);
      } else {
        deltaMat = deltaTMP;
      }
      
      deltaMatStrain_[storageIdx] = deltaMat;
      deltaMatStrain_requiresReeval_[storageIdx] = false;
      
    } else {
      
      if(deltaMat_requiresReeval_[storageIdx] == false){
        LOG_DBG(coeffcthyst) << "Reuse old deltaMat";
        //std::cout << "Reuse old deltaMat" << std::endl;
        return deltaMat_[storageIdx];
      }
      
      /*
       * Compute matrix dM, such that
       * dP = dM*dE (or dB)
       * 
       * with dP = Polarization(timeLevel1)-Polarizatoin(timeLevel2)
       *      dE = E(timeLevel1)-E(timeLevel2)
       */
      numRows = dim_;
      deltaMat = Matrix<Double>(numRows,numCols);
      
      Vector<Double> E_B_new = RetrieveLPMSolution(actualLPM, storageIdx, timelevel_new, onBoundary);
      Vector<Double> E_B_old = RetrieveLPMSolution(actualLPM, storageIdx, timelevel_old, onBoundary);
      
      LOG_DBG(coeffcthystdeltamat) << "Old solution of system (solution at timelevel = " << timelevel_old << "): " << E_B_old.ToString();
      LOG_DBG(coeffcthystdeltamat) << "Current solution of system (solution at timelevel = " << timelevel_new << "): " << E_B_new.ToString();
      
      // does not work well if true
      bool cut_E_B_toSat = !true;
      if(cut_E_B_toSat){
        
        Double satValue;
        if(PDEName_ == "Electromagnetics"){
          // will not  be exactly correct when taking the norm but gives beter estimate than skipping it
          satValue = POL_operatorParams_.outputSat_ + POL_eps_mu_SmallSignal_.NormL2()*POL_operatorParams_.inputSat_;
        } else {
          satValue = POL_operatorParams_.inputSat_;
        }
        
        if(E_B_new.NormL2() > satValue){
          E_B_new.ScalarMult(satValue/E_B_new.NormL2());
        }
        if(E_B_old.NormL2() > satValue){
          E_B_old.ScalarMult(satValue/E_B_old.NormL2());
        }
        
        LOG_DBG(coeffcthystdeltamat) << "After scaling to saturation: ";
        LOG_DBG(coeffcthystdeltamat) << "Old solution of system (solution at timelevel = " << timelevel_old << "): " << E_B_old.ToString();
        LOG_DBG(coeffcthystdeltamat) << "Current solution of system (solution at timelevel = " << timelevel_new << "): " << E_B_new.ToString();
        
      }
      
      
      //std::cout << "Old solution (solution at timelevel = " << timelevel_old << "): " << E_B_old.ToString() << std::endl;
      //std::cout << "Current solution (solution at timelevel = " << timelevel_new << "): " << E_B_new.ToString() << std::endl;
      
      Vector<Double> E_B_diff = E_B_new;
      E_B_diff -= E_B_old;
      
      Vector<Double> P_J_new = GetPrecomputedOutputOfHysteresisOperator(Originallpm,timelevel_new,false);
      Vector<Double> P_J_old = GetPrecomputedOutputOfHysteresisOperator(Originallpm,timelevel_old,false);
      
      LOG_DBG(coeffcthystdeltamat) << "Old output of hyst operator: " << P_J_old.ToString();
      LOG_DBG(coeffcthystdeltamat) << "Current output of hyst operator: " << P_J_new.ToString();
      
      //std::cout << "Old output of hystoperator (solution at timelevel = " << timelevel_old << "): " << P_J_old.ToString() << std::endl;
      //std::cout << "Current output of hystoperator (solution at timelevel = " << timelevel_new << "): " << P_J_new.ToString() << std::endl;
      
      Vector<Double> P_J_diff = P_J_new;
      P_J_diff -= P_J_old;
      
      Double cuttingTol = 1e-5*P_J_old.NormL2();
      
      Matrix<Double> deltaTMP = CalcDeltaMat(E_B_diff, P_J_diff, useAbs, implementationVersion,cuttingTol);
      if( (includeOldDelta_ > 0)&&(includeOldDelta_ <= 1)){
        deltaMat.Init();
        deltaMat.Add(1.0-includeOldDelta_,deltaTMP);
        deltaMat.Add(includeOldDelta_,deltaMatPrev_[storageIdx]);
      } else {
        deltaMat = deltaTMP;
      }
      
      deltaMat_[storageIdx] = deltaMat;
      deltaMat_requiresReeval_[storageIdx] = false;
      
    }
    LOG_DBG(coeffcthystdeltamat) << "Computed deltaMat for storageIDX " << storageIdx << ": " << deltaMat.ToString();
    return deltaMat;
    
  }
  
  Matrix<Double> CoefFunctionHyst::GetIrreversibleStrainTensor(const LocPointMapped& Originallpm, int timeLevel) {
    
    Vector<Double> Si_voigt = GetIrreversibleStrains(Originallpm, timeLevel);
    return ConvertFromVoigtToTensor(Si_voigt);
  }
  Matrix<Double> CoefFunctionHyst::ConvertFromVoigtToTensor(Vector<Double> Si_voigt){
    Matrix<Double> Si_tensor;
    if(Si_voigt.GetSize() == 3){
      Si_tensor = Matrix<Double>(2,2);
    } else {
      Si_tensor = Matrix<Double>(3,3);
    }
    // transform voigt notation to tensor
    // TODO: check implementation of [c]
    //      > do we have Si = sxx,syy,szz,2szy,2sxz,2sxy or do we have Si = sxx,syy,szz,szy,sxz,sxy ?
    if(dim_ == 2){
      Si_tensor[0][0] = Si_voigt[0];
      Si_tensor[1][1] = Si_voigt[1];
      Si_tensor[0][1] = Si_voigt[2];
      Si_tensor[1][0] = Si_voigt[2];
    } else {
      Si_tensor[0][0] = Si_voigt[0];
      Si_tensor[1][1] = Si_voigt[1];
      Si_tensor[2][2] = Si_voigt[2];
      Si_tensor[1][2] = Si_voigt[3];
      Si_tensor[2][1] = Si_voigt[3];
      Si_tensor[0][2] = Si_voigt[4];
      Si_tensor[2][0] = Si_voigt[4];
      Si_tensor[0][1] = Si_voigt[5];
      Si_tensor[1][0] = Si_voigt[5];
    }
    return Si_tensor;
  }
  
  // compute irreversible strains Si
  // Note: Return Vector in VoigtNotation!
  Vector<Double> CoefFunctionHyst::GetIrreversibleStrains(const LocPointMapped& Originallpm, int timeLevel) {
    
    if(tensorType_ == AXI){
      EXCEPTION("ComputeIrreversibleStrains only implemented for 2d plane and full 3d setups");
    }
    /*
     * Model irreversible strains by polynomial of polarization/output of hyst operator
     * ( see Numerical Simulation of Mechatronic Sensors and Actuators 3rd Edition p.386
     * 
     * Si = 3/2*(beta0 + beta1*|P| + beta2*|P|^2 + ... + betan*|P|^n)*(dirP*dirP^T - 1/3[I])
     * 
     */
    
    // Initial step: 
    // check timeLevel
    // only if timeLevel == 0 (current value) and reevaluation is required
    // we evaluate Si; otherwise reuse value
    UInt operatorIdx, storageIdx;
		LocPointMapped actualLPM;
    
		PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx);
    
		if (timeLevel == -1) {
			// return value from last time step
			return Si_lastTS_[storageIdx];
		} else if (timeLevel == 1) {
			// return value from last iteration
			return Si_lastIt_[storageIdx];
		} else {
			// get current state; here we have to check if a reevaluation is needed
			// note: requiresReeval_ has one entry for each storage, not for each operator 
      //			if (false == Si_requiresReeval_[storageIdx]) {
      //				std::cout << "no reeval needed as flag is false!" << std::endl;
      // return current value
      return Si_[storageIdx];
      //			}
		}
    //    
    //    // computation required
    //    // get polarizaton / output of hyst operator first
    //    Vector<Double> P = GetPrecomputedOutputOfHysteresisOperator(Originallpm, timeLevel);
    //    // for case that P is zero, we also obtain the last value in order to compute the last known direction of P
    //    // if this one is 0, too then dir = zeroVec
    //    Vector<Double> Pold = GetPrecomputedOutputOfHysteresisOperator(Originallpm, -1);
    //    
    //    Vector<Double> S_irr = ComputeIrreversibleStrains(P,Pold);
    //    
    //    // flag gets reset at end of each iteration (so the stored value can beu
    //    // used during one iteration by several terms/function calls)
    //    Si_requiresReeval_[storageIdx] = false;
    //    Si_[storageIdx] = S_irr;
    //    
    //    return S_irr;
  }
  
  Vector<Double> CoefFunctionHyst::ComputeIrreversibleStrains(Vector<Double> P, Vector<Double> E_H, Vector<Double> Pold) {
    
    // setup storage
    Matrix<Double> Si_tensor = Matrix<Double>(dim_,dim_);
    Si_tensor.Init();
    
    UInt sizeVoigt;
    if(dim_ == 3){
      sizeVoigt = 6;
    } else {
      // rememeber: axi not supported
      sizeVoigt = 3;
    }
    
    Vector<Double> Si_voigt = Vector<Double>(sizeVoigt);
    Si_voigt.Init(0.0);
    Double scalarStrain = 0.0;
    
    //    std::cout << "P for strain: " << P.ToString() << std::endl;
    //    std::cout << "Norm of P: " << P.NormL2() << std::endl;
    //    
    Vector<Double>dirP = Vector<Double>(dim_);
    dirP.Init();
    if(P.NormL2() != 0){
      dirP.Add(1.0/P.NormL2(),P);
    } else if (Pold.NormL2() != 0){
      dirP.Add(1.0/Pold.NormL2(),Pold);
    }
    // else: current and previous state are both 0 > 0 direction
    
    if (CouplingParams_.irrStrainForm_ == 0){
//        Model irreversible strains as described by Felix Wolf in his PHD-Thesis
//        S = c1 + |p(e) + c2| + (e - 0.5)*c3
//        p,e = normalized polarization,electric field; normalized to range [-0.5,0.5]
      Double pNormalized = P.NormL2();
      Double eNormalized = E_H.Inner(dirP);
      Double pMax,eMax;
      Double toDiff;
      if(CouplingParams_.paramsForHalfRange_ == 1){
        // normalize p,e to twice the sat value
        pNormalized /= (2*POL_operatorParams_.outputSat_);
        eNormalized /= (2*POL_operatorParams_.inputSat_);
        pMax = 0.5;
        eMax = 0.5;
        toDiff = 0.5;
      } else {
        // normalize p,e, to sat value
        pNormalized /= (POL_operatorParams_.outputSat_);
        eNormalized /= (POL_operatorParams_.inputSat_);
        pMax = 1;
        eMax = 1;
        toDiff = 1.0; // as e is now in range [-1 to 1] subtract 1 instead of 0.5
      }
      scalarStrain = CouplingParams_.c1_ + abs(pNormalized + CouplingParams_.c2_)
              + CouplingParams_.c3_*(eNormalized - toDiff);
      
      // false by default for this model
      if(CouplingParams_.scaleTosSat_ == 1){
        Double maxStrain = CouplingParams_.c1_ + abs(pMax + CouplingParams_.c2_)
              + CouplingParams_.c3_*(eMax - toDiff);
        scalarStrain *= (CouplingParams_.sSat_/maxStrain);
      }
      
    } else if (CouplingParams_.irrStrainForm_ == 1){
//        Model irreversible strains as described used M. Loeffler in his muDat identification script
//        (arising from his Diplomarbeit)
//        S = strainSat*( sum_i=0^N c_i*p(e)^i + d0*e + d1*e^2 )
//        note: p,e = normalized polarization,electric field; normalized to range [-0.5,0.5]; sum starts at 0
      Double pNormalized = P.NormL2();
      Double eNormalized = E_H.Inner(dirP);
      if(CouplingParams_.paramsForHalfRange_ == 1){
        // normalize p,e to twice the sat value
        pNormalized /= (2*POL_operatorParams_.outputSat_);
        eNormalized /= (2*POL_operatorParams_.inputSat_);
      } else {
        // normalize p,e, to sat value
        pNormalized /= (POL_operatorParams_.outputSat_);
        eNormalized /= (POL_operatorParams_.inputSat_);
      }
      
      // use Horner scheme
      // start with beta n
      scalarStrain = CouplingParams_.ci_[0][CouplingParams_.ci_size_-1];
      for(int i = CouplingParams_.ci_size_-2; i >= 0; i--){
        scalarStrain = scalarStrain*pNormalized + CouplingParams_.ci_[0][i];
      }
      scalarStrain += CouplingParams_.d0_*eNormalized;
      scalarStrain += CouplingParams_.d1_*eNormalized*eNormalized;

      // true by default for this model
      if(CouplingParams_.scaleTosSat_ == 1){
        scalarStrain *= CouplingParams_.sSat_;
      }
    } else if (CouplingParams_.irrStrainForm_ == 2){
//        Model irreversible strains as described by Manfred in his book
//        S = sum_i=1^N c_i*P(E)^i 
//        note: P,E are polarization and electric field but not normalized! sum starts at 1 (no c0 for offset!)
            
      Double pNorm = P.NormL2();
      // use Horner scheme
      // start with beta n
      scalarStrain = CouplingParams_.ci_[0][CouplingParams_.ci_size_-1];
      Double maxStrain = CouplingParams_.ci_[0][CouplingParams_.ci_size_-1];
      for(int i = CouplingParams_.ci_size_-2; i >= -1; i--){
        scalarStrain = scalarStrain*pNorm;
        maxStrain = maxStrain*POL_operatorParams_.outputSat_;
        // here ci array ends with c1 (and has to be multiplied by P, too
        // to do so, go over rest of loop but do not load [-1] index 
        if(i != -1){
          scalarStrain += CouplingParams_.ci_[0][i];
          maxStrain += CouplingParams_.ci_[0][i];
        }
      }

      // false by default for this model
      if(CouplingParams_.scaleTosSat_ == 1){
        scalarStrain *= (CouplingParams_.sSat_/maxStrain);
      }
    } else {
      return Si_voigt;
      //EXCEPTION("Implementation of irrStrains not known");
    }
    
    // old version; does not fit to parameter determination script and was therefore changed to version above
//    
//    if (CouplingParams_.irrStrainForm_ == 0){
//      // Model strains as discribed in "Generalisiertes Preisach-Modell für die Simulation und Kompensation
//      // der Hysterese piezokeramischer Wandler" - PHD Thesis, Felix Wolf
//      // Original form:
//      //   S = c1 + abs( Hyst(E) + c2) + (E - ESat)/ESat * c3
//      // Modified and used form:
//      //   S = strainSat*(c1 + abs( P(E).NormL2()/Psat + c2) + (inner(E,dirP) - ESat)/ESat * c3) (1D case)
//      // for 2D and 3D extend similar as in case of polynomial approximation
//      // i.e.,
//      // [S] = 3/2*S*(dirP dirP^T - 1/3[I])
//      scalarStrain = CouplingParams_.c1_ + abs(P.NormL2()/POL_operatorParams_.outputSat_ + CouplingParams_.c2_)
//              + CouplingParams_.c3_*(E_H.Inner(dirP) - POL_operatorParams_.inputSat_)/POL_operatorParams_.inputSat_;
//      scalarStrain *= COUPLED_sSat_;
//    } else if (CouplingParams_.irrStrainForm_ == 1){
//      // Model strains (similar) as discribed in "Numerical Simulation of Mechatronic Sensors and Actuators"
//      // Original form:
//      //  [S] = 3/2*(beta0 + beta1*|P| + beta2*|P|^2 + ... + betan*|P|^n)*(dirP*dirP^T - 1/3[I])
//      // Modified and used form:
//      //   S = strainSat*(beta0 + beta1 * P(E)/Psat + beta2 * (P(E)/Psat)^2 + beta3 * (P(E)/Psat)^3 ...)/sum(beta_i)  (1D case)
//      //  [S] = 3/2*S*(dirP dirP^T - 1/3[I])
//      Double normP = P.NormL2();
//      
//      // NEW approach: compute Si via
//      // Si = Ssat*3/2*(beta0_ + beta1_*|P|/Psat + beta2_*(|P|/Psat)^2 + ... + betan_*(|P|/Psat)^n)*(dirP*dirP^T - 1/3[I])/sumOfBetas
//      // Ssat is a new parameter 
//      // betai_ = beta/Psat^i will be read from input instead of beta
//      // sumOfBetas = (beta0_ + beta1_ + beta2_ + ... + betan_)
//      // why the change? to allow for easier adaption to amplitude Ssat
//      normP = normP/POL_operatorParams_.outputSat_;  
//      // use Horner scheme
//      // start with beta n
//      scalarStrain = COUPLED_betaCoefs_[0][COUPLED_dim_beta_-1];
//      Double betaSum = COUPLED_betaCoefs_[0][COUPLED_dim_beta_-1];
//      for(int i = COUPLED_dim_beta_-2; i >= 0; i--){
//        scalarStrain = scalarStrain*normP + COUPLED_betaCoefs_[0][i];
//        betaSum += COUPLED_betaCoefs_[0][i];
//      }
//      if(betaSum != 0){
//        scalarStrain = COUPLED_sSat_*scalarStrain/betaSum;
//      }
//    } else {
//      return Si_voigt;
//    }
    
    // bring to 2D/3D form
    Matrix<Double> negeye = Matrix<Double>(dim_,dim_);
    for(UInt i = 0; i < dim_; i++){
      negeye[i][i] = -1.0;
    }
    
    Matrix<Double>dyadic = Matrix<Double>(dim_,dim_);
    dyadic.DyadicMult(dirP);
    dyadic.Add(1.0/3.0,negeye);
    
    Si_tensor.Add(1.5*scalarStrain,dyadic);
    
    // transform matrix to voigt notation
    // TODO: check implementation of [c]
    //      > do we have Si = sxx,syy,szz,2szy,2sxz,2sxy or do we have Si = sxx,syy,szz,szy,sxz,sxy ?
    if(dim_ == 2){
      Si_voigt[0] = Si_tensor[0][0];
      Si_voigt[1] = Si_tensor[1][1];
      Si_voigt[2] = Si_tensor[0][1];
    } else {
      Si_voigt[0] = Si_tensor[0][0];
      Si_voigt[1] = Si_tensor[1][1];
      Si_voigt[2] = Si_tensor[2][2];
      Si_voigt[3] = Si_tensor[1][2];
      Si_voigt[4] = Si_tensor[0][2];
      Si_voigt[5] = Si_tensor[0][1];
    }
    //    std::cout << "Si_voigt = " << Si_voigt.ToString() << std::endl;
    
    return Si_voigt;   
  }
  
	Vector<Double> CoefFunctionHyst::GetPrecomputedOutputOfHysteresisOperator(const LocPointMapped& Originallpm, 
          int timeLevel, bool forStrain) {
    LOG_TRACE(coeffcthyst) << "GetPrecomputedOutputOfHysteresisOperator for timelevel " << timeLevel;
    //std::cout << "Get Output of Hyst Operator" << std::endl;
    //std::cout << "Timelevel: " << timeLevel << " (0 = current, -1 = lastTS, +1 = lastIt)" << std::endl;
		/*
     * Input:
     *  LocPointMapped specifying the integration point where the state of the hyst operator
     *  shall be obtained
     *
     *  UInt timeLevel:
     *		-1: get value from last ts
     *    0: get current value
     *    1: get value from last iteration
     *
     * Flags:
     *  bool invert_
     *    true: use inverse hyst operator if possible; otherwise perform iterative inversion
     *    false: use standard hyst operator
     *
     * Return value:
     *  Vector defining the state of the hyst operator
     */
    
		assert(XMLParameterSet_);
		totalCallingCounter_++;
    
		/*
     * 1. preprocess lpm to get information about the actual lpm where we have to
     *  evaluate the hyst operator, the index of the operator and the index of the storage array
     */
		UInt operatorIdx, storageIdx;
		LocPointMapped actualLPM;
    
    //		bool onBoundary =     PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx);
    
    PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx);
    
		/*std::cout << "timeLevel " << timeLevel << std::endl;
     std::cout << "storageIdx: " << storageIdx << std::endl;
     std::cout << "size of arrays: " << P_J_lastTS_->GetSize() << std::endl;
     */
    /*
     * 2. check time level and need for evaluation
     */    
		if(timeLevel == -2){
      Vector<Double> zeroVec = Vector<Double>(dim_);
			zeroVec.Init();
			return zeroVec;
		} else if (timeLevel == -1) {
			LOG_TRACE(coeffcthyst) << "> prevTS value";
      // return value from last time step
      //std::cout << "Get value from last TS: " << P_J_lastTS_[storageIdx].ToString() << std::endl;
      if(forStrain){
        return P_J_forStrains_lastTS_[storageIdx];
      } else {
        return P_J_lastTS_[storageIdx];
      }
		} else if (timeLevel == 1) {
      LOG_TRACE(coeffcthyst) << "> prevIt value";
			// return value from last iteration
      if(forStrain){
        EXCEPTION("Polarization for strains for last iteration should not be called");
      } else {
        return P_J_lastIt_[storageIdx];
      }
		} else {
      if(forStrain){
        return P_J_forStrains_[storageIdx];
      } else {
        return P_J_[storageIdx];
      }
      
      //      OLD
      //      
      //      LOG_TRACE(coeffcthyst) << "> current value";
      //      if(hystOperatorLocked_){
      //        // no evaluation > return current state
      //				LOG_TRACE(coeffcthyst) << "Output: hystOperatorLocked_ locked";
      //        return P_J_[storageIdx];
      //      }
      //      
      //			// get current state; here we have to check if a reevaluation is needed
      //			// note: requiresReeval_ has one entry for each storage, not for each operator
      //			if (false == requiresReeval_[storageIdx]) {
      //        LOG_TRACE(coeffcthyst) << "> NO reeval; take stored value";
      //				//std::cout << "no reeval needed as flag is false!" << std::endl;
      //				// return current value
      //				return P_J_[storageIdx];
      //			}
		}
    //    LOG_TRACE(coeffcthyst) << "> reeval needed";
    //		/*
    //     * 3. Evaluate hysteresis operator
    //     */
    //		// first get input
    //		Vector<Double> input = RetrieveInputToHysteresisOperator(actualLPM, operatorIdx, storageIdx, onBoundary);
    //    
    //		// then compute output
    //		return CalcOutputOfHysteresisOperator(input, operatorIdx, storageIdx);
    
	}
  
	Vector<Double> CoefFunctionHyst::RetrieveLPMSolution(LocPointMapped& actualLPM, UInt storageIdx, int timeLevel, bool onBoundary) {
    
		if (timeLevel == -1) {
			// get solution from last ts
			return E_B_lastTS_[storageIdx];
		} else if (timeLevel == 1) {
			// get solution from last iteration
			return E_B_lastIt_[storageIdx];
		} else {
			// get current input from system
			Vector<Double> LPMSolution;
			if(onBoundary){
				if(dependCoef1Surf_ != NULL){
					LPMSolution = Vector<Double>(dependCoef1Surf_->GetVecSize());
					dependCoef1Surf_->GetVector(LPMSolution, actualLPM);
				} else {
					EXCEPTION("LPM solution on boundary requested, but coefFunction on boundary not set.")
				}
				//std::cout << "Retrieved input ON boundary: " << LPMSolution.ToString() << std::endl;
				
				// for fieldParallel boundary conditions, En should be zero, Pn not
				// (En should be zero due to NeumannBC but this does not seem to work, so 
				// we might have to subtract normal part)
				bool forceTangential = true;
        
        if(forceTangential){
          
          Vector<Double> normalDirection = actualLPM.normal;
          //          std::cout << "Stored normal vector: " << normalDirection.ToString() << std::endl;
          Double normalProjection = LPMSolution.Inner(normalDirection);
          //          std::cout << "LPMSolution from system: " << LPMSolution.ToString() << std::endl;
          LPMSolution.Add(-normalProjection,normalDirection);
          //          std::cout << "LPMSolution after cutting normal direction: " << LPMSolution.ToString() << std::endl;
        }
        
			} else {
				LPMSolution = Vector<Double>(dependCoef1_->GetVecSize());
				dependCoef1_->GetVector(LPMSolution, actualLPM);
			}
			
      //			if(POL_operatorParams_.angularClipping_ > 0){
      //				/*
      //				 * clipping
      //				 */
      //				ClipDirection(LPMSolution);
      //			}
			
			return LPMSolution;
		}
	}
  
	Vector<Double> CoefFunctionHyst::RetrieveInputToHysteresisOperator(LocPointMapped& actualLPM, UInt operatorIdx, UInt storageIdx, bool onBoundary) {
    LOG_TRACE(coeffcthyst) << "RetrieveInputToHysteresisOperator";
    //std::cout << "RetrieveInputToHysteresisOperator" << std::endl;
		// get current solution first
		// as we want the current input only, set timeLevel to 0
    
    /*
     * NOTE/TODO:
     *  in case of magnetostrictive coupling, we actually would have to solve the nonlinear
     *  problem
     *    H = -[h](S - Si(H)) + nu(B - P(H))
     *  
     * in case of vector model (where we use Levenberg Marquart) we could actually build in
     * this whole function (which would drastically increase the computational effort)
     * in case of scalar model, we can only invert for 
     *      nu(B - P(H)) = H
     * 
     * furthermore, nu depends on P itself, such that we have even more trouble 
     * 
     * --> currently, we only get H_partial such that
     *         H_partial = nu(P_old) (B - P(H_partial))
     *    
     */
    if(hystOperatorLocked_){
      // no evaluation > return current state
			LOG_TRACE(coeffcthyst) << "Input: hystOperatorLocked_ locked";
      return E_H_[storageIdx];
    }
    
    
		UInt timeLevel = 0;
		Vector<Double> curLPMSolution = RetrieveLPMSolution(actualLPM, storageIdx, timeLevel, onBoundary);
//    std::cout << "curLPMSolution: " << curLPMSolution.ToString() << std::endl;
		if (needsInversion_ && (POL_operatorParams_.hasInverseModel_ == false) ) {
      LOG_TRACE(coeffcthyst) << "Inversion needed";
      //std::cout << "Inversion needed" << std::endl;
      // we have an inverse field (e.g. magnetics) but no inverse model (i.e. VectorPreisach)
      //
			// check if solution did change since last time > if not, hope that input / output pair of hyst
			// operator can be reused
			Vector<Double> diff = curLPMSolution;
			Vector<Double> toDiff = E_B_[storageIdx];
			diff -= toDiff;
      
      //std::cout << "StorageIDX: " << storageIdx << std::endl;
      //std::cout << "current Solution (input): " << curLPMSolution.ToString() << std::endl;
      //std::cout << "old Solution (loaded): " << E_B_[storageIdx] << std::endl;
      //std::cout << "difference: " << diff.ToString() << std::endl;
      //std::cout << "diff.NormL2(): " << diff.NormL2() << std::endl;
      Vector<Double> retrievedInput;
			if (diff.NormL2() < POL_operatorParams_.amplitudeResolution_) {
        LOG_TRACE(coeffcthyst) << "> difference " << diff.NormL2() << " below tolerance > reuse";
        //std::cout << "(Nearly) no difference in amplitude to previous element solution" << std::endl;
        retrievedInput = E_H_[storageIdx];
				//}
			} else {
        LOG_TRACE(coeffcthyst) << "> calculate new input";
        //bool overwriteMemory = OverwriteHystMemory();
        // for input computation never overwrite storage
        // we basically search in the current state for a fitting solution
        // then, during the step calcoutput we evaluate the hyst operator with
        // the retrieved input and here we actually overwrite the memory, if necessary
        bool overwriteMemory = false;
        
        if(matrixForInversionWasSet_[storageIdx] == false){
          EXCEPTION("Material matrix for inversion was not set yet! This has to be done in the calling helper classes!");
        } else {
          LOG_TRACE(coeffcthyst) << "Matrix for inversion (storage="<<storageIdx<<"): "<<matrixForInversion_[storageIdx].ToString();
        }
        if(POL_operatorParams_.methodType_ == 0){
          retrievedInput = Vector<Double>(dim_);
					
          //					retrievedInput[POL_operatorParams_.fixDirection_] = hyst_->computeInputAndUpdate(curLPMSolution[POL_operatorParams_.fixDirection_], 
          //						E_B_lastIt_[storageIdx][POL_operatorParams_.fixDirection_], E_H_lastIt_[storageIdx][POL_operatorParams_.fixDirection_],P_J_lastIt_[storageIdx][POL_operatorParams_.fixDirection_],
          //						operatorIdx, POL_eps_mu_SmallSignal_[POL_operatorParams_.fixDirection_][POL_operatorParams_.fixDirection_], overwriteMemory);
          // NEW: POL_operatorParams_.fixDirection_ is a vector now
          // > to get the correct scalar input data, we have to project onto this direction vector
          Double muForInversion;
          Vector<Double> tmp = Vector<Double>(dim_);
          matrixForInversion_[storageIdx].Mult(POL_operatorParams_.fixDirection_,tmp);
          POL_operatorParams_.fixDirection_.Inner(tmp,muForInversion);
          
          Double scalInput;
          POL_operatorParams_.fixDirection_.Inner(curLPMSolution,scalInput);
          int successFlag = 0;
          Double scalOutput = hyst_->computeInputAndUpdate(scalInput,muForInversion, operatorIdx, overwriteMemory, successFlag);
          
          retrievedInput.Init();
          retrievedInput.Add(scalOutput,POL_operatorParams_.fixDirection_);
          
          /*
           * NEW:
           *  as we use retrievedInput later to output H directly, we have to add the missing components
           *  in the scalar case
           * so far we have H in direction dirP; to get H perpendicular to that direction, we compute the
           * following:
           * 1) B_perpendicular = B - dirP*(dirP.inner(B)
           * 2) H_perpendicular = nu*B_perpendicualr
           * 3) add H_perpendicular to H
           */
          Vector<Double> curLPMSolution_perpendicular = Vector<Double>(dim_);
          curLPMSolution_perpendicular.Init();
          curLPMSolution_perpendicular.Add(1.0,curLPMSolution,-scalInput,POL_operatorParams_.fixDirection_);
          Vector<Double> retrievedInput_perpendicular = Vector<Double>(dim_);
          matrixForInversionInverted_[storageIdx].Mult(curLPMSolution_perpendicular,retrievedInput_perpendicular);
          retrievedInput.Add(1.0,retrievedInput_perpendicular);
          
          //          retrievedInput[POL_operatorParams_.fixDirection_] = hyst_->computeInputAndUpdate(curLPMSolution[POL_operatorParams_.fixDirection_], 
          //						muForInversion, operatorIdx, overwriteMemory);
        } else {
          // Idea: use Levenberg Marquart to estimate  inversion of hyst operator
          //          retrievedInput = hyst_->computeInput_vec(curLPMSolution, operatorIdx, POL_eps_mu_SmallSignal_,
          //                  alphaLinesearch_, overwriteMemory, RUN_overwriteDirection_ );
					
          //					retrievedInput = hyst_->computeInput_vec(curLPMSolution, E_B_lastIt_[storageIdx], E_H_lastIt_[storageIdx],
          //						P_J_lastIt_[storageIdx],operatorIdx,POL_eps_mu_SmallSignal_,RUN_overwriteDirection_);
          //          
          int successFlag = 0;

          retrievedInput = hyst_->computeInput_vec(curLPMSolution,operatorIdx,matrixForInversion_[storageIdx],
                  POL_operatorParams_.fieldsAlignedAboveSat_,POL_operatorParams_.hystOutputRestrictedToSat_,successFlag);
          
        }
        //std::cout << "Computed input: " << retrievedInput.ToString() << std::endl;
        // we should not store the value here but first after evaluating the hyst operator in calcOutput function
        // otherwise we will never get to the point where P_J_ gets overwritten
        // and thus we get no real stepping at all
        //E_H_[storageIdx] = retrievedInput;
			}
      // store current solution for later checks
      E_B_[storageIdx] = curLPMSolution;
      
      //std::cout << "Returned input: " << retrievedInput << std::endl;
      //      if(POL_operatorParams_.angularClipping_ > 0){
      //				/*
      //				 * clipping; note: LPMSolution is already clipped in RetrieveLPMSolution
      //				 */
      //				ClipDirection(retrievedInput);
      //			}
			
			return retrievedInput;
      
		} else {
      // store current solution for later checks
      E_B_[storageIdx] = curLPMSolution;
      //std::cout << "Returned input: " << curLPMSolution << std::endl;
			// for electrostatic / piezoelectric material, we can use the element solution (E) as input
			return curLPMSolution;
		} 
	}
  
	void CoefFunctionHyst::ClipDirection(Vector<Double>& targetVector){
    /*
     * EDIT: 12.4.2018
     *  Clipping does not work as expected; 
     *  Reason:
     *    If we clip the hyst output (P), we shrink the space of possible solutions.
     *    Consider the magnetic case with B pointing in 45 degree angle
     *    Due to clipping, P just might be allowed to have 44 or 46 degree so the difference
     *    in B has to be compensated by H. Assume that P is at 44 degree. Now H will point 
     *    towards 130+ degree to make up for the remaining part of B. Due to small mu, H
     *    will be very large and thus, will have a huge influence on P, eventually turing it
     *    to 46 degree. Now H has to change to the opposite direction to compensate for B.
     *    This will turn P to 44 degree again etc.
     *  Possible new purpose:
     *    Magnetics: clip B > smaller target solution space; maybe easier to find a solution via inversion
     *    Electrostatics: clip E > hyst operator will ignore tiny rotations in E, might help in terms of robustness
     *                              but could also backfire as fine tuning might be harder 
     * 
     * 
     * 
     * (Former name: clipNewRotationDirection, taken from VecPreisach)
     * New purpose March/April 2018
     * The general clipping idea might be ok, but it does not help to clip the rotation states
     * directly
     *  > reason: due to the weighting via the switching state, we still get a continuous
     *            resolution, even if the single rotations states were clipped
     *  > new purpose: instead of clipping the single rotation states, we clip the final output
     *            of the hyst operator!
     *  > further improvings/changes:
     *      1. 3d case implemented
     *      2. angularClipping is treated as resolution in DEGREE (not rad!)
     * 
     * Addon: shifted to CoefFunctionHyst; both input and output of hyst operator
     *				can be clipped
     * 
     * New function April 2017
     * Idea: Restrict the range of possible rotation directions to fixed angular steps
     * (i.e. x rad)
     * Reason: Due to numerical pollutions, we might encounter rotation directions like
     * (1.00000000,-1e-9) or (0.99999999,1e-10). These deviations from the input direction
     * (1.0,0) seem negligble but will lead to actual problems when evaluting deltaMatrices
     * of the form deltaP/deltaE (as deltaE normally is similarly small). Due to this, the
     * resulting deltaMatrices may have huge permittivities/permeabilites in direction where
     * normally no value should be. This leads to serious convergence issues.
     *
     */
    
    /*
     * The following two clipping steps will be applied: 
     * 1. transform to circular/spherical coordinates; restrict angles (in deg) to angularClipping
     * 2. check for special angles (90,180) degree and restrict vector further
     *      (as e.g. sin(pi) != 0 due to numerics)
     */
    LOG_TRACE(coeffcthyst) << "Clip direction of target vector to next full " << POL_operatorParams_.angularClipping_ << " degree";
    LOG_DBG(coeffcthyst) << "Original vector: " << targetVector.ToString();
    if(dim_ == 2){
      /*
       * use polar/circular coordinates
       * x = r cos(alpha)
       * y = r sin(alpha)
       */
      Double radius, tmp, alphaDeg;
      radius = targetVector.NormL2();
      
      if(targetVector[0] == 0){
        if(targetVector[1] > 0){
          alphaDeg = 90.0; 
        } else {
          alphaDeg = -90.0;
        }
      } else if(targetVector[0] > 0){
        tmp = atan2(targetVector[1],targetVector[0]);
        alphaDeg = tmp*180/M_PI;
      } else {
        tmp = -asin(targetVector[1]/radius);
        alphaDeg = tmp*180/M_PI + 180.0;
      }
      LOG_DBG(coeffcthyst) << "Circular/Polar coordinates (r,alpha): " << radius << "," << alphaDeg;
      
      // apply clipping
      tmp = alphaDeg/POL_operatorParams_.angularClipping_;
      tmp = round(tmp);
      alphaDeg = tmp*POL_operatorParams_.angularClipping_;
      
      LOG_DBG(coeffcthyst) << "Circular/Polar coordinates after clipping (r,alpha): " << radius << "," << alphaDeg;
      
      // now rebuild output vector with clipped coordinates
      targetVector[0] = radius*cos(alphaDeg/180*M_PI);
      targetVector[1] = radius*sin(alphaDeg/180*M_PI);
      
      LOG_DBG(coeffcthyst) << "Rebuild vector after clipping: " << targetVector.ToString();
      
      // finally check for special cases i.e. 90/180 deg (which can not perfectly be reproduced by computing cos()/sin()
      if( abs(alphaDeg - 90) < POL_operatorParams_.angularClipping_/1000.0 ){
        // alpha = 90 > positive y axis
        targetVector[0] = 0.0;
        targetVector[1] = radius;
      } else if( abs(alphaDeg + 90) < POL_operatorParams_.angularClipping_/1000.0 ){
        // alpha = -90 > negative y axis
        targetVector[0] = 0.0;
        targetVector[1] = -radius;
      } else if( (abs(alphaDeg - 180) < POL_operatorParams_.angularClipping_/1000.0) || (abs(alphaDeg + 180) < POL_operatorParams_.angularClipping_/1000.0) ){
        // alpha = +/- 180 deg > negative x axis
        targetVector[0] = -radius;
        targetVector[1] = 0.0;
      }
      
      LOG_DBG(coeffcthyst) << "Rebuild vector after further treatment: " << targetVector.ToString();
    } else {
      /*
       * use spherical coordinates
       * x = r sin(theta) cos(phi)
       * y = r sin(theta) sin(phi)
       * z = r cos(theta)
       */
      Double radius, tmp, phiDeg, thetaDeg;
      radius = targetVector.NormL2();
      
      tmp = acos(targetVector[2]/radius);
      thetaDeg = tmp*180/M_PI;
      tmp = atan2(targetVector[1],targetVector[0]);
      phiDeg = tmp*180/M_PI;
      
      LOG_DBG(coeffcthyst) << "Spherical coordinates (r,theta,phi): " << radius << "," << thetaDeg << "," << phiDeg;
      
      // apply clipping
      tmp = thetaDeg/POL_operatorParams_.angularClipping_;
      tmp = round(tmp);
      thetaDeg = tmp*POL_operatorParams_.angularClipping_;
      
      tmp = phiDeg/POL_operatorParams_.angularClipping_;
      tmp = round(tmp);
      phiDeg = tmp*POL_operatorParams_.angularClipping_;
      
      LOG_DBG(coeffcthyst) << "Spherical coordinates after clipping (r,theta,phi): " << radius << "," << thetaDeg << "," << phiDeg;
      
      targetVector[0] = radius*sin(thetaDeg/180*M_PI)*cos(phiDeg/180*M_PI);
      targetVector[1] = radius*sin(thetaDeg/180*M_PI)*sin(phiDeg/180*M_PI);
      targetVector[2] = radius*cos(thetaDeg/180*M_PI);
      
      LOG_DBG(coeffcthyst) << "Rebuild vector after clipping: " << targetVector.ToString();
      
      // finally check for special cases i.e. 90/180 deg (which can not perfectly be reproduced by computing cos()/sin()
      if( (abs(thetaDeg - 90) < POL_operatorParams_.angularClipping_/1000.0 ) || (abs(thetaDeg + 90) < POL_operatorParams_.angularClipping_/1000.0 ) ){
        // theta = +/- 90 > z = 0
        targetVector[2] = 0.0;
      } else if( (abs(thetaDeg - 180) < POL_operatorParams_.angularClipping_/1000.0) || (abs(thetaDeg + 180) < POL_operatorParams_.angularClipping_/1000.0) ){
        // theta = +/- 180 deg > negative z-axis
        targetVector[0] = 0.0;
        targetVector[1] = 0.0;
        targetVector[2] = -radius;
      }
      if( (abs(phiDeg - 90) < POL_operatorParams_.angularClipping_/1000.0 ) || (abs(phiDeg + 90) < POL_operatorParams_.angularClipping_/1000.0 ) ){
        // phi = +/- 90 > x = 0
        targetVector[0] = 0.0;
      } else if( (abs(phiDeg - 180) < POL_operatorParams_.angularClipping_/1000.0) || (abs(phiDeg + 180) < POL_operatorParams_.angularClipping_/1000.0) ){
        // phi = +/- 180 deg > y = 0
        targetVector[1] = 0.0;
      }
      
      LOG_DBG(coeffcthyst) << "Rebuild vector after further treatment: " << targetVector.ToString();
      
    }
  }
	
  Vector<Double> CoefFunctionHyst::GetPrecomputedInputToHysteresisOperator(const LocPointMapped& Originallpm, int timeLevel){
    
		UInt operatorIdx, storageIdx;
		LocPointMapped actualLPM;
    
    PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx);
    
    if(timeLevel == -1){
      return E_H_lastTS_[storageIdx];
    } else if(timeLevel == -1){
      return E_H_lastIt_[storageIdx];
    } else {
      return E_H_[storageIdx];
    }
  }
  
  void CoefFunctionHyst::EvaluateHystOperatorsInt(Integer intFlag){
    
    bool setMatForInversion = false;
    bool resetSolToZeroFirst = false;
    bool overwriteMemory = false;
    bool setRotStateForVecExtension = false;
    
    if(abs(intFlag) >= 1){
      setMatForInversion = true;
    }
    if(abs(intFlag) == 2){
      resetSolToZeroFirst = true;
    }
    if(abs(intFlag) == 3){
      setRotStateForVecExtension = true;
    }
    
    if(intFlag < 0){
      overwriteMemory = true;
    }
    
    CoefFunctionHyst::EvaluateHystOperators(setMatForInversion, resetSolToZeroFirst, overwriteMemory, setRotStateForVecExtension);
    
  }
  
  void CoefFunctionHyst::EvaluateHystOperators(bool setMatForInversion, bool resetSolToZeroFirst, bool overwriteMemory, bool setRotStateForVecExtension){
    
    /*
     * NEW function added April 16th 2018
     * Idea/background:
     *  There are two central functions that need to be executed before the actual evaluation
     *  of the hyst operatros can take place
     *  a) hystHelper->InitLinearTensors 
     *      > has to be executed once during the first time step
     *  b) hystHelper->PrecomputeMaterialTensorForInverison()
     *      > has to be exectued for magnetic case every time that solution changes
     *        (like evaluation of hyst operators itself)
     *  Till now, these two functions above were executed during the first call to
     *  hystHelper->GetVector or hystHelper->GetTensor. This works well for single
     *  threads but for parallel execution segfaults appear sometimes which leads to
     *  the assumption that there are some race conditions present.
     * 
     * > Solution approach:
     *  Evaluate these two critial functions as well as all hyst operators BEFORE
     *  GetVector or GetTensor get called and set the flags ...requiresReeval_ to false
     *  such that GetVector and GetTensor (which call GetOutputOfHystOperator) will use this
     *  precalculated values;
     * 
     *  Advantages:
     *  - Execution done only once for each solution vector
     *  - Execution done at one central place
     *  - Execution can be parallel for non-critial parts 
     * 
     */
    
    /*
     * Create helper class if not done yet;
     * > this helper class is normally created alongside the material or rhs
     *    classes during setup of the pdes
     * > for coupled case, we have to make sure, that couplTensorFct and elastTensorFct were
     *    set during setup of the coupling
     */
    if(hystHelper_ == NULL){
      
      PtrCoefFct fieldTensor;
      if (material_->GetMaterialDatabaseName() == "Electrostatic") {
        fieldTensor = material_->GetTensorCoefFnc( ELEC_PERMITTIVITY,tensorType_,Global::REAL);
      } else if (material_->GetMaterialDatabaseName() == "Electromagnetics") {
        fieldTensor = material_->GetTensorCoefFnc( MAG_RELUCTIVITY,tensorType_,Global::REAL);
      }
      hystHelper_ = new CoefFunctionHelper(fieldTensor,elastTensorFct_,couplTensorFct_, hystItself_);
    }
    
    /*
     * Init lienar Tensors; as we assume these tensors to be constant for each
     * region, we can take any lpm for init
     */
    if(allLPMmap_.empty()){
      EXCEPTION("Map of LPM was not set yet");
    }
    LocPointMapped anyLPM = allLPMmap_.begin()->second;
    hystHelper_->InitLinearTensors(anyLPM);
    
    /*
     * PrecomputeMaterialTensorForInverison() will access the polarization state
     * of the PREVIOUS ITERATION, i.e. the values stored in P_J_lastIt_ in order
     * to evaluate the matrix for inversion;
     * if we reset sol to 0 during solvestephyst, we might think about setting these values
     * to 0, too
     */
    if(resetSolToZeroFirst){
      for(UInt i = 0; i < numStorageEntries_; i++){
        P_J_lastIt_[i].Init();
      }
    }
    
    /*
     * Compute matrices which are needed for local inversion
     * (only for magnetics due to inversion of Preisach operator)
     * this function uses the last known values of the hystOperator in case
     * that the matrix depends on it (only for coupled case where nu can
     * be dependent on h which depends on hysteresis output)
     */
    if((setMatForInversion)&&(needsInversion_)){
      hystHelper_->PrecomputeMaterialTensorForInverison();
    }
    
    //    if(setRotStateForVecExtension){
    //      EvaluateRotDirectionForVectorExtension();
    //    }
    
    if(resetSolToZeroFirst){
      for(UInt i = 0; i < numStorageEntries_; i++){
        LOG_DBG(coeffcthyst) << "matrixForInversion_[ " << i << " ] = " << matrixForInversion_[i];
      }
    }
    /*
     * Evaluate all hyst operators at corresponding storage entries
     *  > standard evaluation: one hyst operator for each midpoint; evaluation only at midpoint
     *  > extended evaluation: one hyst operator for each midpoint; evaluation at each integration point + midpoint
     *  > full evaluation: one hyst operator for each integration point + midpoint; evaluation at each integration point + midpoint
     */
    PrecomputePolarization(overwriteMemory);
    
    if(resetSolToZeroFirst){
      for(UInt i = 0; i < numStorageEntries_; i++){
        LOG_DBG(coeffcthyst) << "E_H_[ " << i << " ] = " << E_H_[i];
        LOG_DBG(coeffcthyst) << "E_B_[ " << i << " ] = " << E_B_[i];
        LOG_DBG(coeffcthyst) << "P_J_[ " << i << " ] = " << P_J_[i];
      }
    }
    /*
     * -1 : not coupled at all
     *  0 : coupled e-form/h-form (piezo/magstrict)
     *  1 : coupled d-form (piezo)
     *  2 : coupled g-form (magstrict)
     */
    if(CouplingParams_.strainForm_ != -1){
//      std::cout << "Coupled" << std::endl;
      PrecomputeIrrStrains();
      bool rotate = true;
      Matrix<Double> baseTensor = hystHelper_->GetBaseCouplingTensor();
      PrecomputeScaledAndRotatedCouplingTensor(baseTensor,rotate);
    }
    
    //    
    //#ifdef USE_OPENMP
    //#pragma omp parallel num_threads(CFS_NUM_THREADS)
    //    {
    //      std::map<UInt, std::map<UInt, bool> >::iterator itStart = operatorToStorage_.begin();
    //      std::map<UInt, std::map<UInt, bool> >::iterator itEnd;
    //      std::map<UInt, std::map<UInt, bool> >::iterator it;
    //      UInt size = operatorToStorage_.size();
    //      
    //      //std::cout << "Running parallel" << std::endl;
    //      UInt numT = CFS_NUM_THREADS;
    //      UInt aThread = omp_get_thread_num();
    //      //std::cout << "My number: " << aThread << std::endl;
    //      UInt chunksize = std::floor(size/numT);
    //      
    //      std::advance(itStart, aThread*chunksize);
    //      
    //      if(aThread == numT-1){
    //        itEnd = operatorToStorage_.end();
    //      } else {
    //        itEnd = itStart;
    //        std::advance(itEnd, chunksize);
    //      }  
    //      //UInt cnt = 0;
    //      //UInt cntIn = 0;
    //#pragma omp barrier
    //      for(it = itStart; it != itEnd; it++){
    //        //      UInt operatorIdx = it->first;
    //        //cnt++;
    //        std::map<UInt, bool>::iterator innerIt;
    //        for(innerIt = it->second.begin(); innerIt != it->second.end(); innerIt++){
    //          //cntIn++;
    //          UInt storageIdx = innerIt->first;
    //          //        bool isMidpoint = innerIt->second;
    //          LocPointMapped curLoc = allLPMmap_[storageIdx];
    //          GetPrecomputedOutputOfHysteresisOperator(curLoc, 0);
    //        }
    //      } 
    //#pragma omp barrier
    //      //std::cout << "My number " << aThread << "; did " << cnt << " operators with total " << cntIn << " evals" << std::endl;
    //    }
    //#else
    //    {
    //    std::map<UInt, std::map<UInt, bool> >::iterator itStart = operatorToStorage_.begin();
    //    std::map<UInt, std::map<UInt, bool> >::iterator itEnd = operatorToStorage_.end();
    //    std::map<UInt, std::map<UInt, bool> >::iterator it;
    //    
    //    for(it = itStart; it != itEnd; it++){
    ////      UInt operatorIdx = it->first;
    //      std::map<UInt, bool>::iterator innerIt;
    //      for(innerIt = it->second.begin(); innerIt != it->second.end(); innerIt++){
    //        UInt storageIdx = innerIt->first;
    ////        bool isMidpoint = innerIt->second;
    //        LocPointMapped curLoc = allLPMmap_[storageIdx];
    //        GetPrecomputedOutputOfHysteresisOperator(curLoc, 0);
    //      }
    //    } 
    //#endif
  }
  
  void CoefFunctionHyst::PrecomputeIrrStrains(){
    LOG_DBG(coeffcthyst) << "CoefFunctionHyst::PrecomputeIrrStrains()";
    /*
     * Iterate over all lpm
     * for each lpm extract the precomputed values of P (current), P (lastTimestep)
     * use those values to calculate the irreversible strains and store them
     */
    
    bool forStrain = true;
    UInt size = allLPMmap_.size();
    UInt strainSize = 3;
    if(dim_ == 3){
      strainSize = 6;
    }
    
#ifdef USE_OPENMP
#pragma omp parallel num_threads(CFS_NUM_THREADS)
    {     
      Vector<Double> P = Vector<Double>(dim_);
      Vector<Double> E_H = Vector<Double>(dim_);
      Vector<Double> Pold = Vector<Double>(dim_);
      Vector<Double> S_irr = Vector<Double>(strainSize);
      //std::cout << "Running parallel" << std::endl;
      UInt numT = CFS_NUM_THREADS;
      UInt aThread = omp_get_thread_num();
      //std::cout << "My number: " << aThread << std::endl;
      UInt chunksize = std::floor(size/numT);
      
      std::map<UInt, LocPointMapped >::iterator LPMit;
      std::map<UInt, LocPointMapped >::iterator LPMitStart = allLPMmap_.begin();
      std::map<UInt, LocPointMapped >::iterator LPMitEnd;
      
      std::advance(LPMitStart, aThread*chunksize);
      
      if(aThread == numT-1){
        LPMitEnd = allLPMmap_.end();
      } else {
        LPMitEnd = LPMitStart;
        std::advance(LPMitEnd, chunksize);
      }  
      
#pragma omp barrier
      for(LPMit = LPMitStart; LPMit != LPMitEnd; LPMit++){
        // computation required
        // get current polarizaton / output of hyst operator (timelevel = 0)
        P = GetPrecomputedOutputOfHysteresisOperator(LPMit->second, 0, forStrain);
        // for case that P is zero, we also obtain the last value in order to compute the last known direction of P
        // if this one is 0, too then dir = zeroVec
        Pold = GetPrecomputedOutputOfHysteresisOperator(LPMit->second, -1, forStrain);
        
        E_H = GetPrecomputedInputToHysteresisOperator(LPMit->second, 0);
        
        S_irr = ComputeIrreversibleStrains(P,E_H,Pold);
        
        // flag gets reset at end of each iteration (so the stored value can beu
        // used during one iteration by several terms/function calls)
        Si_requiresReeval_[LPMit->first] = false;
        Si_[LPMit->first] = S_irr;
      }
#pragma omp barrier
         }
#else
    {
      Vector<Double> P = Vector<Double>(dim_);
      Vector<Double> E_H = Vector<Double>(dim_);
      Vector<Double> Pold = Vector<Double>(dim_);
      Vector<Double> S_irr = Vector<Double>(strainSize);
      //std::cout << "Running serial" << std::endl;
      
      std::map<UInt, LocPointMapped >::iterator LPMit;
      
      for(LPMit = allLPMmap_.begin(); LPMit != allLPMmap_.end(); LPMit++){
        // computation required
        // get current polarizaton / output of hyst operator (timelevel = 0)
        P = GetPrecomputedOutputOfHysteresisOperator(LPMit->second, 0, forStrain);
        // for case that P is zero, we also obtain the last value in order to compute the last known direction of P
        // if this one is 0, too then dir = zeroVec
        Pold = GetPrecomputedOutputOfHysteresisOperator(LPMit->second, -1, forStrain);
        
        E_H = GetPrecomputedInputToHysteresisOperator(LPMit->second, 0);
        
        S_irr = ComputeIrreversibleStrains(P,E_H,Pold);
        
        // flag gets reset at end of each iteration (so the stored value can beu
        // used during one iteration by several terms/function calls)
        Si_requiresReeval_[LPMit->first] = false;
        Si_[LPMit->first] = S_irr;
      }
    }
#endif
  };
  
  void CoefFunctionHyst::PrecomputeScaledAndRotatedCouplingTensor(Matrix<Double> baseTensor, bool rotate){
    LOG_DBG(coeffcthysthelper) << "CoefFunctionHelper::PrecomputeScaledAndRotatedCouplingTensor()";
    
    /*
     * Iterate over all integration points used/managed by the hyst operator and
     * for each of these LPM set the matrix for iversion, i.e. mu/nu in magnetic case
     */    
    
    UInt size = allLPMmap_.size();
    
    // obtain correct sizes first
    UInt numRows, numCols;
    numCols = baseTensor.GetNumCols();
    numRows = baseTensor.GetNumRows();
    
    if(numCols == 4){
      // axi case: 2x4 matrix
      WARN("Rotation for axi case not implemented yet. No rotation will be peformed");
      rotate = false;
    }
    
#ifdef USE_OPENMP
#pragma omp parallel num_threads(CFS_NUM_THREADS)
    {
      Matrix<Double> scaledCouplTensor = Matrix<Double>(numRows,numCols);
      Matrix<Double> rotatedCouplTensor = Matrix<Double>(numRows,numCols);
      Matrix<Double> R;
      Matrix<Double> RT;
      
      Vector<Double> P = Vector<Double>(dim_);
      
      // get current polarization (elec or mag)
      Vector<Double> dirP = Vector<Double>(dim_);
      Double alpha, beta, gamma, scaling;
      //      std::cout << "Current polarization vector " << P.ToString() << std::endl;
      //      
      // calculate scaling
      
      //std::cout << "Running parallel" << std::endl;
      UInt numT = CFS_NUM_THREADS;
      UInt aThread = omp_get_thread_num();
      //std::cout << "My number: " << aThread << std::endl;
      UInt chunksize = std::floor(size/numT);
      
      std::map<UInt, LocPointMapped >::iterator LPMit;
      std::map<UInt, LocPointMapped >::iterator LPMitStart = allLPMmap_.begin();
      std::map<UInt, LocPointMapped >::iterator LPMitEnd;
      
      std::advance(LPMitStart, aThread*chunksize);
      
      if(aThread == numT-1){
        LPMitEnd = allLPMmap_.end();
      } else {
        LPMitEnd = LPMitStart;
        std::advance(LPMitEnd, chunksize);
      }  
      
#pragma omp barrier
      for(LPMit = LPMitStart; LPMit != LPMitEnd; LPMit++){
        
        // get current value of P
        Vector<Double> P = GetPrecomputedOutputOfHysteresisOperator(LPMit->second, 0, false);
        scaling = P.NormL2()/POL_operatorParams_.outputSat_;
        
        scaledCouplTensor = baseTensor*scaling;
        
        //      std::cout << "Scalaed coupling tensor " << scaledCouplTensor.ToString() << std::endl;
        //      
        if(rotate == true){
          if(P.NormL2() == 0){
            //        std::cout << "Polarization is zero > perform no rotation" << std::endl;
            rotatedCouplTensor = scaledCouplTensor;
            
          } else {
            dirP = P / P.NormL2();
            
            if (numCols == 6){
              // 3d plane case
              // Idea: Create rotation matrix, that maps the current direction of P onto z-axis
              //       as the z-axis is the default polarization axis in the mat file
              //       To obtained the desired behavior, we have to rotate the z-axis onto P however.
              //       This can be done by taking the transposed rotation matrix.
              
              // 1. rotate around z-axis by angle gamma such that P lies in z-y plane
              // gamma = angle between z-y plane and dirP
              gamma = std::atan2(dirP[0],dirP[1]);
              
              // 2. rotate around x-axis by angle alpha such that P lies on top of z-axis
              // alpha = angle between x-y plane and z
              alpha = std::atan2(std::sqrt(dirP[1]*dirP[1]+dirP[0]*dirP[0]),dirP[2]);
              
              // no rotation arouind y-axis needed > beta = 0
              // WARNING: this whole procedure is only valid if coupling tensor is at least
              // transverse isotropic! Otherwise we have to figure out on which axis to rotate
              // the transverse directions
              beta = 0.0;
              R = Compute3DRotationMatrix(alpha, beta, gamma);
              
              // take transpose matrix to revert rotation (i.e. rotate z-axis onto dirP)
              RT.Resize(3,3);
              R.Transpose(RT);
              
              //              assert(rotatedCouplTensor.GetNumRows() == 3);
              //              assert(rotatedCouplTensor.GetNumCols() == 6);
              
              scaledCouplTensor.PerformRotation(RT,rotatedCouplTensor);
              
            } else {
              // 2d plane strain or plane stress case
              // Important remark:
              //  for 2d plane strain and stress (and somehow also for axi) the coupling tensor
              //  will be rotated by default by alpha = -90 and gamma = -90.
              //  This results in the following mapping of the coordinate axis:
              //    z > y, y > x, x > z
              //  I.e. the material is rotated such that the default polarization is in +y direction.
              //  We thus have two possibilities to align the material tensor to the 2d polarization
              //  direction:
              //  a) get original 3x6 tensor and rotate that tensor according to the 3d version above
              //      then cut out subtensor
              //  b) take cut out subtensor and perform rotation directly in 2d
              //      > Version b done in the following
              
              // std::atan2(dirP[0],dirP[1]) would rotate towards the y-axis
              // we want however the y-axis to rotate, so that we take -gamma instead
              gamma = -std::atan2(dirP[0],dirP[1]);
              
              R = Compute2DRotationMatrix(gamma);
              
              //              assert(rotatedCouplTensor.GetNumRows() == 2);
              //              assert(rotatedCouplTensor.GetNumCols() == 3);
              scaledCouplTensor.PerformRotation(R,rotatedCouplTensor);
            }
          }
        } else {
          // no rotation > take scale tensor
          rotatedCouplTensor = scaledCouplTensor;
        }
        
        rotatedCouplingTensor_requiresReeval_[LPMit->first] = false;
        rotatedCouplingTensor_[LPMit->first] = rotatedCouplTensor;
        //        lastUsedTimeLevelForRotation_[storageIdx] = timeLevel;
      }
#pragma omp barrier
    }
#else
    {
      Matrix<Double> scaledCouplTensor = Matrix<Double>(numRows,numCols);
      Matrix<Double> rotatedCouplTensor = Matrix<Double>(numRows,numCols);
      Matrix<Double> R;
      Matrix<Double> RT;
      
      Vector<Double> P = Vector<Double>(dim);
      
      // get current polarization (elec or mag)
      Vector<Double> dirP = Vector<Double>(dim);
      Double alpha, beta, gamma, scaling;
      //      std::cout << "Current polarization vector " << P.ToString() << std::endl;
      //      
      // calculate scaling
      
      //std::cout << "Running serial" << std::endl;      
      std::map<UInt, LocPointMapped >::iterator LPMit;
      std::map<UInt, LocPointMapped >::iterator LPMitStart = allLPMmap_.begin();
      std::map<UInt, LocPointMapped >::iterator LPMitEnd = allLPMmap_.end();
      
#pragma omp barrier
      for(LPMit = LPMitStart; LPMit != LPMitEnd; LPMit++){
        
        // get current value of P
        Vector<Double> P = GetPrecomputedOutputOfHysteresisOperator(LPMit->second, 0);
        scaling = P.NormL2()/POL_operatorParams_.outputSat_;
        
        scaledCouplTensor = baseTensor*scaling;
        
        //      std::cout << "Scalaed coupling tensor " << scaledCouplTensor.ToString() << std::endl;
        //      
        if(rotate == true){
          if(P.NormL2() == 0){
            //        std::cout << "Polarization is zero > perform no rotation" << std::endl;
            rotatedCouplTensor = scaledCouplTensor;
            
          } else {
            dirP = P / P.NormL2();
            
            if (numCols == 6){
              // 3d plane case
              // Idea: Create rotation matrix, that maps the current direction of P onto z-axis
              //       as the z-axis is the default polarization axis in the mat file
              //       To obtained the desired behavior, we have to rotate the z-axis onto P however.
              //       This can be done by taking the transposed rotation matrix.
              
              // 1. rotate around z-axis by angle gamma such that P lies in z-y plane
              // gamma = angle between z-y plane and dirP
              gamma = std::atan2(dirP[0],dirP[1]);
              
              // 2. rotate around x-axis by angle alpha such that P lies on top of z-axis
              // alpha = angle between x-y plane and z
              alpha = std::atan2(std::sqrt(dirP[1]*dirP[1]+dirP[0]*dirP[0]),dirP[2]);
              
              // no rotation arouind y-axis needed > beta = 0
              // WARNING: this whole procedure is only valid if coupling tensor is at least
              // transverse isotropic! Otherwise we have to figure out on which axis to rotate
              // the transverse directions
              beta = 0.0;
              R = Compute3DRotationMatrix(alpha, beta, gamma);
              
              // take transpose matrix to revert rotation (i.e. rotate z-axis onto dirP)
              RT.Resize(3,3);
              R.Transpose(RT);
              
              //              assert(rotatedCouplTensor.GetNumRows() == 3);
              //              assert(rotatedCouplTensor.GetNumCols() == 6);
              
              scaledCouplTensor.PerformRotation(RT,rotatedCouplTensor);
              
            } else {
              // 2d plane strain or plane stress case
              // Important remark:
              //  for 2d plane strain and stress (and somehow also for axi) the coupling tensor
              //  will be rotated by default by alpha = -90 and gamma = -90.
              //  This results in the following mapping of the coordinate axis:
              //    z > y, y > x, x > z
              //  I.e. the material is rotated such that the default polarization is in +y direction.
              //  We thus have two possibilities to align the material tensor to the 2d polarization
              //  direction:
              //  a) get original 3x6 tensor and rotate that tensor according to the 3d version above
              //      then cut out subtensor
              //  b) take cut out subtensor and perform rotation directly in 2d
              //      > Version b done in the following
              
              // std::atan2(dirP[0],dirP[1]) would rotate towards the y-axis
              // we want however the y-axis to rotate, so that we take -gamma instead
              gamma = -std::atan2(dirP[0],dirP[1]);
              
              R = Compute2DRotationMatrix(gamma);
              
              //              assert(rotatedCouplTensor.GetNumRows() == 2);
              //              assert(rotatedCouplTensor.GetNumCols() == 3);
              scaledCouplTensor.PerformRotation(R,rotatedCouplTensor);
            }
          }
        } else {
          // no rotation > take scale tensor
          rotatedCouplTensor = scaledCouplTensor;
        }
        
        rotatedCouplingTensor_requiresReeval_[LPMit->first] = false;
        rotatedCouplingTensor_[LPMit->first] = rotatedCouplTensor;
        //        lastUsedTimeLevelForRotation_[storageIdx] = timeLevel;
      }
#pragma omp barrier
    }
#endif
    
  };
	
	Vector<Double> CoefFunctionHyst::CalcOutputOfHysteresisOperator(Vector<Double> inputToHystOperator, UInt operatorIdx, UInt storageIdx, 
          bool forceMemoryLock, bool forceMemoryWrite, bool useOperatorForStrain) {
    
    //std::cout << "CalcOutputOfHysteresisOperator" << std::endl;
    //std::cout << "InputToHystOperator: " << inputToHystOperator.ToString() << std::endl;
    
		Vector<Double> outputOfHystOperator;
    
		/*
     *  overwrite the hysteresis memory
     *          or
     *  work on temporal storage
     */
    // during standard evaluation, the helperfunction OverwriteHystMemory()
    // will return the correct flag; however, for setting the previous hyst value
    // it might be that we are forced to write or lock the memory
    bool overwriteMemory;
    if(forceMemoryLock == true){
      overwriteMemory = false;
    } else if(forceMemoryWrite == true){
      overwriteMemory = true;
    } else {
      overwriteMemory = OverwriteHystMemory();
    }
    
		/*
     * Check if a reevaluation really is required by comparing the actual input to
     * the last used input
     */
		Vector<Double> diff = inputToHystOperator;
		Vector<Double> toDiff = E_H_[storageIdx];
		diff -= toDiff;
    
    if(hystOperatorLocked_){
      // no evaluation > return current state
			LOG_TRACE(coeffcthyst) << "Calc output: hystOperatorLocked_ locked";
      //			std::cout << "hystOperatorLocked_ locked" << std::endl;
      if(useOperatorForStrain){
        return P_J_forStrains_[storageIdx];
      } else {
        return P_J_[storageIdx];
      }
    }
    
		if ((diff.NormL2() < POL_operatorParams_.amplitudeResolution_)&&(overwriteMemory == false)) {
      //	if (XML_textOutputLevel_ == 2) {
      //std::cout << "(Nearly) no difference in amplitude to previous input" << std::endl;
      //std::cout << "Input: " << inputToHystOperator.ToString() << std::endl;
      //std::cout << "Old: " << toDiff.ToString() << std::endl;
      //std::cout << "Apmplitude resolution: " << POL_operatorParams_.amplitudeResolution_ << std::endl;
      //	}
      //std::cout << "return P_J_[storageIdx]: " << P_J_[storageIdx].ToString() << std::endl;
      if(useOperatorForStrain){
        return P_J_forStrains_[storageIdx];
      } else {
        return P_J_[storageIdx];
      }
		}
    //std::cout << "Evaluate Hysteresis operator" << std::endl;
		/*
     * Call hystersis operator
     */
		totalEvaluationCounter_++;
    
    ParameterPreisachOperators operatorParams;
    if(useOperatorForStrain){
      operatorParams = STRAIN_operatorParams_;
    } else {
      operatorParams = POL_operatorParams_;
    }
    
		if (operatorParams.methodType_ == 0) {
      
			outputOfHystOperator = Vector<Double>(dim_);
			outputOfHystOperator.Init();
      Vector<Double> curDir;
      
			if (XML_performanceMeasurement_ == 1) {
				timer_->Start();
			}
      
      //      if(POL_useExtension_){
      //        // get latest rotation state from extension
      //        curDir = hyst_->getRotationDirection(operatorIdx);
      //      } else {
      curDir = operatorParams.fixDirection_;
      //      }
      //std::cout << "curDir: " << curDir.ToString() << std::endl;
      //std::cout << "Used input for SCALAR model: " << inputToHystOperator[POL_operatorParams_.fixDirection_] << std::endl;
      //      outputOfHystOperator[POL_operatorParams_.fixDirection_] = hyst_->computeValueAndUpdate(inputToHystOperator[POL_operatorParams_.fixDirection_], 
      //				E_H_[storageIdx][POL_operatorParams_.fixDirection_],operatorIdx, overwriteMemory);
      
      // NEW: POL_operatorParams_.fixDirection_ is a vector now
      // > to get the correct scalar input data, we have to project onto this direction vector
      Double scalInput;
      Double scalOutput;
      curDir.Inner(inputToHystOperator,scalInput);
      int successFlag = 0;
      
      if(useOperatorForStrain){
        scalOutput = hystStrain_->computeValueAndUpdate(scalInput,operatorIdx, overwriteMemory,successFlag);
      } else {
        scalOutput = hyst_->computeValueAndUpdate(scalInput,operatorIdx, overwriteMemory,successFlag);
      }
      
      outputOfHystOperator.Init();
      outputOfHystOperator.Add(scalOutput,curDir);
      
			//outputOfHystOperator[POL_operatorParams_.fixDirection_] = hyst_->computeValueAndUpdate(inputToHystOperator[POL_operatorParams_.fixDirection_],operatorIdx, overwriteMemory);
			
			if (XML_performanceMeasurement_ == 1) {
				timer_->Stop();
				//        std::cout << "Scalar Preisach operator" << std::endl;
				//        std::cout << "Total time (" << totalEvaluationCounter_ << "calls): " << timer_->GetCPUTime() << std::endl;
				//        std::cout << "Avg Time: " << timer_->GetCPUTime()/totalEvaluationCounter_ << std::endl;
				totalEvaluationTime_ = timer_->GetCPUTime();
				avgEvaluationTime_ = timer_->GetCPUTime() / totalEvaluationCounter_;
			}
      
		} else {
			if (XML_performanceMeasurement_ == 1) {
				timer_->Start();
			}
      //std::cout << "Comput hyst output using Vector model " << std::endl;
			//std::cout << "OperatorIdx: " << operatorIdx << std::endl;
      int successFlag = 0;
      
      if(useOperatorForStrain){
        outputOfHystOperator = hyst_->computeValue_vec(inputToHystOperator, operatorIdx, overwriteMemory, 
                STRAIN_operatorParams_.fieldsAlignedAboveSat_, successFlag);
      } else {
        outputOfHystOperator = hyst_->computeValue_vec(inputToHystOperator, operatorIdx, overwriteMemory, 
                POL_operatorParams_.fieldsAlignedAboveSat_, successFlag);
      }

			if (XML_performanceMeasurement_ == 1) {
				timer_->Stop();
				//        std::cout << "Vector Preisach operator" << std::endl;
				//        std::cout << "Total time (" << totalEvaluationCounter_ << "calls): " << timer_->GetCPUTime() << std::endl;
				//        std::cout << "Avg Time: " << timer_->GetCPUTime()/totalEvaluationCounter_ << std::endl;
				totalEvaluationTime_ = timer_->GetCPUTime();
				avgEvaluationTime_ = timer_->GetCPUTime() / totalEvaluationCounter_;
			}
      
			/*!
       * Print out of Hysteresis state
       * ONLY FOR DEBUGGING REASONS!
       * USE CAREFULLY! MASSIVE IMAGE OUTPUT!
       * only for vector case
       */
			static UInt cnt = 0;
			static UInt firstIdx = 1;
      
      if(!useOperatorForStrain){
        if ((POL_operatorParams_.printOut_ > 0) && (RUN_allowBMP_ == true)) {
          if ((cnt % POL_operatorParams_.printOut_ == 0)&&(operatorIdx == firstIdx)) {
            //std::cout << "Outputting bmp" << std::endl;
            std::stringstream filenamebuf;
            filenamebuf << "Switch_Elem" << firstIdx << "_Step" << std::setfill('0') << std::setw(5) << cnt << "_v" << POL_operatorParams_.evalVersion_ << "_numRows" << POL_numRows_ << ".bmp";
            hyst_->switchingStateToBmp(POL_operatorParams_.bmpResolution_, filenamebuf.str(), operatorIdx, true);
          }
          
          if (operatorIdx == firstIdx) {
            cnt++;
            /*
             * disable output until reset
             * -> otherwise we would get two images for each timestep if P and D are computed
             */
            RUN_allowBMP_ = false;
          }
        }
      }
    }
    //    
    //		if(POL_operatorParams_.angularClipping_ > 0){
    //			/*
    //			 * clipping; note: inputToHystOperator is already clipped
    //			 */
    //			ClipDirection(outputOfHystOperator);
    //		}
		
    if(useOperatorForStrain){
      P_J_forStrains_[storageIdx] = outputOfHystOperator;
    } else {
      //std::cout << "Store input/output combination for later usage" << std::endl;
      E_H_[storageIdx] = inputToHystOperator;
      //std::cout << "P_J_[storageIdx] before setting: " << P_J_[storageIdx].ToString() << std::endl;
      P_J_[storageIdx] = outputOfHystOperator;
    }
    //std::cout << "P_J_[storageIdx] after setting: " << P_J_[storageIdx].ToString() << std::endl;
    
    //std::cout << "Computed output: " << outputOfHystOperator.ToString() << std::endl;
		// this flag has to be reset by hand!
		requiresReeval_[storageIdx] = false;
    
		return outputOfHystOperator;
    
	}
  
  void CoefFunctionHyst::ForceRemanence(){
		LOG_TRACE(coeffcthyst) << "ForceRemanence";
    /*
     * Iterate over all operators and evaluate them at 0
     * > at this time also set the last input to the hyst operator and set arrays
     */
    Vector<Double> output;
    Vector<Double> zeroVec = Vector<Double>(dim_);
    zeroVec.Init();
    
    bool forceMemoryLock = true;
    UInt operatorIdx, storageIdx;
    
    LocPointMapped actualLPM;
    std::map<UInt, LocPointMapped>::iterator listIt;
    std::map<UInt, LocPointMapped> currentMap;
		currentMap = allLPMmap_;
		
    for(listIt = currentMap.begin(); listIt!= currentMap.end(); ++listIt){
      LOG_TRACE(coeffcthyst) << "Storage idx: " << listIt->first;
      bool forceMidpoint;
      
      if(midpointLPMmap_.find(listIt->first) != midpointLPMmap_.end()){
        forceMidpoint = true;
      } else {
        forceMidpoint = false;
      }
      
      PreprocessLPM(listIt->second, actualLPM, operatorIdx, storageIdx, forceMidpoint);
      if(listIt->first != storageIdx){
        EXCEPTION("StorageIdx in map != storageIdx as obtained by PreprocessLPM!");
      }
      
      E_H_[storageIdx] = zeroVec;
      
      // then compute output
      output = CalcOutputOfHysteresisOperator(zeroVec, operatorIdx, storageIdx, forceMemoryLock, false);
			LOG_TRACE(coeffcthyst) << "Computed output: " << output.ToString();
			
      P_J_[storageIdx] = output;
      if (needsInversion_ && (POL_operatorParams_.hasInverseModel_ == false) ){
        // magnetics > H = 0, B = J
        E_B_[storageIdx] = output;
      } else {
        // electrostatics > E = 0, D = P
        E_B_[storageIdx] = zeroVec;
      }
    }
  }
  
  //  void EvaluateAllHystOperators(bool forceMemoryLock){
  //    
  //    // Unfortunately, there are some unresolved race conditions when during the assembly
  //    // process; most probably sources of this problem are the functions
  //    // InitLinearTensors and PrecomputeMaterialTensorForInversion
  //    // To be save of such race conditions, we will manually take care of all evaluations
  //    // of the hyst operators in the following way:
  //    // 1. At the beginning of the very first timestep, we
  //    //    a) Set all linear tensors
  //    //    b) Setup the matrix for local inversion of the Preisach operator
  //    //    c) evaluate each hyst-operator with the initial solution (most probably 0) and
  //    //        the previously set matrix for local inversion
  //    // 2. Every time after the solution vector was updated, i.e. during linesearch,
  //    //    after linesearch, after a possible reset and after convergence; only in
  //    //    the last case, we allow setting of memory
  //    
  //    
  //    
  //  }
  
  void CoefFunctionHyst::PrecomputePolarization(bool overwriteMemory){
    
    UInt RUN_evaluationPurpose_bak = RUN_evaluationPurpose_;
    RUN_evaluationPurpose_ = 1;
    
#ifdef USE_OPENMP
#pragma omp parallel num_threads(CFS_NUM_THREADS)
    {
      std::map<UInt, std::map<UInt, bool> >::iterator itStart = operatorToStorage_.begin();
      std::map<UInt, std::map<UInt, bool> >::iterator itEnd;
      std::map<UInt, std::map<UInt, bool> >::iterator it;
      UInt size = operatorToStorage_.size();
      
      //std::cout << "Running parallel" << std::endl;
      UInt numT = CFS_NUM_THREADS;
      UInt aThread = omp_get_thread_num();
      //std::cout << "My number: " << aThread << std::endl;
      UInt chunksize = std::floor(size/numT);
      
      std::advance(itStart, aThread*chunksize);
      
      if(aThread == numT-1){
        itEnd = operatorToStorage_.end();
      } else {
        itEnd = itStart;
        std::advance(itEnd, chunksize);
      }  
      
      Vector<Double> solution = Vector<Double>(dim_);
      Vector<Double> input = Vector<Double>(dim_);
      Vector<Double> output = Vector<Double>(dim_);
      
      // RUN_evaluationPurpose_ only needed if we want to check
      // a) if we want to evaluate at the midpoint or at every point
      // b) if we want to set hyst memory or not
      // > b is irrelevant here as we set the overwrite flag per hand
      // > a would be irrelevant, too as we just go over whole map of LPM
      //    but the function PreProcessLPM calls EvaluateAtMidPointOnly which
      //    needs the correct evaluationpurpose
      // > set purpose to 1 (assemble/eval)
      
      bool onBoundary;
      bool forceMemoryLock;
      bool forceMemoryWrite;
      
#pragma omp barrier
      for(it = itStart; it != itEnd; it++){
        UInt operatorIdx = it->first;
        UInt cnt = 0;
        
        std::map<UInt, bool>::iterator innerIt;
        for(innerIt = it->second.begin(); innerIt != it->second.end(); innerIt++){
          /*
           * NOTE:
           * for standard and full evaluation, the length of the inner list should
           *    always be 1
           *  > standard eval: 1 operator per element, evaluation only at midpoint 
           *  > full eval: 1 operato per int. point + 1 for midpoint, evaluation at corrsp. point
           * 
           * for extended evaluation, the length should be larger than 1 
           *  > extended eval: 1 operator per element, eval at midpoint and integration points
           *      (hyst memory may only be written at midpoint)
           */
          UInt storageIdx = innerIt->first;
          bool isMidpoint = innerIt->second;
          if(isMidpoint && (cnt != 0)){
            EXCEPTION("Midpoint should always be the first storage entry assigned to each operator (if at all");
          }
          cnt++;
          
          if(overwriteMemory){
            forceMemoryLock = false;
            forceMemoryWrite = true;
          } else {
            forceMemoryLock = true;
            forceMemoryWrite = false;
          }
          
          if(!isMidpoint && (XML_EvaluationDepth_ == 2)){
            // extended evaluation: only overwrite at midpoint!
            forceMemoryLock = true;
            forceMemoryWrite = false;
          }
          
          LocPointMapped curLPM = allLPMmap_[storageIdx];
          LocPointMapped actualLPM;
          
          UInt operatorIdxTMP, storageIdxTMP;
          onBoundary = PreprocessLPM(curLPM, actualLPM, operatorIdxTMP, storageIdxTMP, isMidpoint);
          
          if( (operatorIdxTMP != operatorIdx) || (storageIdxTMP != storageIdx)){
            EXCEPTION("Storage or operator index does not fit");
          }         
          // get input
          input = RetrieveInputToHysteresisOperator(actualLPM, operatorIdx, storageIdx,onBoundary);
          // RetrieveInputToHysteresisOperator will store the actual solution directly to E_B_
          
          // compute additional hyst operator for strains first
          // reason: CalcOutputOfHysteresisOperator sets E_H_ if last flag is false; 
          //  during the next call input will be compared to E_H_
          // > if output is computed for standard hyst operator first, the call for the strain operator will
          //  not evaluate anything as E_H_ == input
          if(CouplingParams_.ownHystOperator_ && COUPLED_inXMLFile_ ){
            output = CalcOutputOfHysteresisOperator(input, operatorIdx, storageIdx, forceMemoryLock, forceMemoryWrite, true); 
          }
          
          // then compute output
          output = CalcOutputOfHysteresisOperator(input, operatorIdx, storageIdx, forceMemoryLock, forceMemoryWrite, false);  
        }
      }
#pragma omp barrier 
    }
    
#else
    Vector<Double> solution = Vector<Double>(dim_);
		Vector<Double> input = Vector<Double>(dim_);
		Vector<Double> output = Vector<Double>(dim_);
    
    // RUN_evaluationPurpose_ only needed if we want to check
    // a) if we want to evaluate at the midpoint or at every point
    // b) if we want to set hyst memory or not
    // > b is irrelevant here as we set the overwrite flag per hand
    // > a would be irrelevant, too as we just go over whole map of LPM
    //    but the function PreProcessLPM calls EvaluateAtMidPointOnly which
    //    needs the correct evaluationpurpose
    // > set purpose to 1 (assemble/eval)
    
		bool onBoundary;
    bool forceMemoryLock;
    bool forceMemoryWrite;
    // NEW APPROACH compared to the SetPreviousHystValues approach
    // 1. iterate over all operators
    // 2. iterate over all assigned evaluation points
    //  > standard: only midpoint; overwrite = overwriteMemory
    //  > extended: midpoint first with overwrite = overwriteMemory, 
    //              then integration points with overwrite = false!
    //  > full: one point per operator; overwrite = overwriteMemory
    std::map<UInt, std::map<UInt, bool> >::iterator itStart = operatorToStorage_.begin();
    std::map<UInt, std::map<UInt, bool> >::iterator itEnd = operatorToStorage_.end();
    std::map<UInt, std::map<UInt, bool> >::iterator it;
    
    for(it = itStart; it != itEnd; it++){
      UInt operatorIdx = it->first;
      UInt cnt = 0;
      
      std::map<UInt, bool>::iterator innerIt;
      for(innerIt = it->second.begin(); innerIt != it->second.end(); innerIt++){
        /*
         * NOTE:
         * for standard and full evaluation, the length of the inner list should
         *    always be 1
         *  > standard eval: 1 operator per element, evaluation only at midpoint 
         *  > full eval: 1 operato per int. point + 1 for midpoint, evaluation at corrsp. point
         * 
         * for extended evaluation, the length should be larger than 1 
         *  > extended eval: 1 operator per element, eval at midpoint and integration points
         *      (hyst memory may only be written at midpoint)
         */
        UInt storageIdx = innerIt->first;
        bool isMidpoint = innerIt->second;
        if(isMidpoint && (cnt != 0)){
          EXCEPTION("Midpoint should always be the first storage entry assigned to each operator (if at all");
        }
        cnt++;
        
        if(overwriteMemory){
          forceMemoryLock = false;
          forceMemoryWrite = true;
        } else {
          forceMemoryLock = true;
          forceMemoryWrite = false;
        }
        
        if(!isMidpoint && (XML_EvaluationDepth_ == 2)){
          // extended evaluation: only overwrite at midpoint!
          forceMemoryLock = true;
          forceMemoryWrite = false;
        }
        
        LocPointMapped curLPM = allLPMmap_[storageIdx];
        LocPointMapped actualLPM;
        
        UInt operatorIdxTMP, storageIdxTMP;
        onBoundary = PreprocessLPM(curLPM, actualLPM, operatorIdxTMP, storageIdxTMP, isMidpoint);
        
        if( (operatorIdxTMP != operatorIdx) || (storageIdxTMP != storageIdx)){
          EXCEPTION("Storage or operator index does not fit");
        }         
        // get input
        input = RetrieveInputToHysteresisOperator(actualLPM, operatorIdx, storageIdx,onBoundary);
        // RetrieveInputToHysteresisOperator will store the actual solution directly to E_B_
        
        // compute additional hyst operator for strains first
        // reason: CalcOutputOfHysteresisOperator sets E_H_ if last flag is false; 
        //  during the next call input will be compared to E_H_
        // > if output is computed for standard hyst operator first, the call for the strain operator will
        //  not evaluate anything as E_H_ == input
        if(CouplingParams_.ownHystOperator_ && COUPLED_inXMLFile_ ){
          output = CalcOutputOfHysteresisOperator(input, operatorIdx, storageIdx, forceMemoryLock, forceMemoryWrite, true); 
        }
        
        // then compute output
        output = CalcOutputOfHysteresisOperator(input, operatorIdx, storageIdx, forceMemoryLock, forceMemoryWrite, false); 
      }
    }
#endif
    
    RUN_evaluationPurpose_ = RUN_evaluationPurpose_bak;
  }
  
  void CoefFunctionHyst::SetPreviousHystVals(bool setLastTS) {
    /*
     * NEW: just copy values from current array to lastIT or lastTS array
     */
    std::map<UInt, std::map<UInt, bool> >::iterator itStart = operatorToStorage_.begin();
    std::map<UInt, std::map<UInt, bool> >::iterator itEnd = operatorToStorage_.end();
    std::map<UInt, std::map<UInt, bool> >::iterator it;
    
    for(it = itStart; it != itEnd; it++){
      std::map<UInt, bool>::iterator innerIt;
      for(innerIt = it->second.begin(); innerIt != it->second.end(); innerIt++){
        UInt storageIdx = innerIt->first;
        
        if (setLastTS) {
          E_B_lastTS_[storageIdx] = E_B_[storageIdx];
          P_J_lastTS_[storageIdx] = P_J_[storageIdx];
          P_J_forStrains_lastTS_[storageIdx] = P_J_forStrains_[storageIdx];
          E_H_lastTS_[storageIdx] = E_H_[storageIdx];
          Si_lastTS_[storageIdx] = Si_[storageIdx];
        } else {
          E_B_lastIt_[storageIdx] = E_B_[storageIdx];
          P_J_lastIt_[storageIdx] = P_J_[storageIdx];
          E_H_lastIt_[storageIdx] = E_H_[storageIdx];
          Si_lastIt_[storageIdx] = Si_[storageIdx];
        }
        deltaMatPrev_[storageIdx] = deltaMat_[storageIdx];
        deltaMatStrainPrev_[storageIdx] = deltaMatStrain_[storageIdx];
      }
    }
  }
  
  /*
   * OLD VERSION
   */
	void CoefFunctionHyst::SetPreviousHystValsOLD(bool setLastTS, bool forceMemoryLock) {
    
    //		if(XML_textOutputLevel_ == 2){
    //      std::cout << "++ SetPreviousHystVals ++" << std::endl;
    //      std::cout << "lastTS? " << setLastTS << std::endl;
    //		}
    
		/*
     * Function to backup input/output pair as well as current deltaMatrix;
     * This function is called at the beginning of each iteration BEFORE the
     * system gets solved. Therewith, the backup always store the values of
     * the last iteration. These values are needed for the computation of the
     * new deltaMatrix. Evaluation of the hysteresis operator is done onto a temp.
     * copy, i.e. the actual memory is not changed).
     *
     * If setLastTS_ is true, the current state is instead stored to
     * the arrays ..._lastTS_ and the evaluation is performed on the actual storage
     * (i.e. the hysteresis operator is set). This call should only be performed
     * after the end of a timestep.
     *
     * Depending on the XML_EvaluationDepth_ this function does the following steps:
     * depth = 1 or 2 (standard or extended):
     *  iterate over all elements
     *    extract element index
     *    extract solution and input to hysteresis operator at midpoint of element
     *    evaluate hysteresis operator with the extracted input
     *    store extracted solution, extracted input and evaluated output of hystoperator
     *     and store current state of deltaMatrix (no reevaluation)
     *
     * depth = 3 (full evalution)
     *  iterate over all elements
     *    extract element index
     *    get integration points for element
     *
     *    iterate over all integration points
     *      compute corresponding index (unique combination of element index and integration point index)
     *      extract solution and input to hystoperator at integration point
     *      evaluate hyst operator with the determined index
     *      store extracted solution, input and output and store current state of
     *        deltaMattrix
     *
     *    additionally: get coordinates of element midpoint (is normally no integration point)
     *    extract solution and input at midpoint
     *    evaluate hyst operator with index of midpoint (again a unique combination of element index + point index)
     *    store extracted solution, ...
     *
     * Note: the additional evaluation at the midpoint is needed for the output computation as we
     *  evaluate the hyst operator at the midpoint; therefore we need to know the state at the midpoint
     *
     */
    
		Vector<Double> solution = Vector<Double>(dim_);
		Vector<Double> input = Vector<Double>(dim_);
		Vector<Double> output = Vector<Double>(dim_);
		UInt storageIdx;
		UInt operatorIdx;
    
		/*
     * 0.1 backup current evaluationPurpose
     */
		UInt evalPurposeBackup = RUN_evaluationPurpose_;
    
		/*
     * 0.2 set evaluation purpose to
     *  2 if setLastTS is false
     *  3 if setLastTS is true
     *
     * by setting the evaluation purpose, we automatically will get the
     * right values for midpointOnly and OverwriteHystMemory from the
     * corresponding functions
     */
		if (setLastTS) {
			// overwriteMemory = true
			// for extended evaluation (where we also evaluate at each integration point)
			// we only overwrite memory for the midpoint (where the hystoperator is located)
			RUN_evaluationPurpose_ = 3;
		} else {
			// overwriteMemory = false
			RUN_evaluationPurpose_ = 2;
		}
    
		bool overwriteMemoryIntPoints = false;
		if (XML_EvaluationDepth_ == 3) {
      // only for full integration (where we have a hyst operator per integration
      // point) we overwrite the memory at each of these points
			overwriteMemoryIntPoints = true;
		}
		bool overwriteMemoryMidPoint = OverwriteHystMemory();
		bool midpointOnly = EvaluateAtMidpointOnly();
    
		if (forceMemoryLock) {
			overwriteMemoryIntPoints = false;
			overwriteMemoryMidPoint = false;
		}
    
    LOG_TRACE(coeffcthyst) << "Set previous hyst values";
    if(setLastTS){
      LOG_TRACE(coeffcthyst) << "For last timestep";
    } else {
      LOG_TRACE(coeffcthyst) << "For last iteration";
    }
    LOG_TRACE(coeffcthyst) << "Midpoint only? " << midpointOnly;
    LOG_TRACE(coeffcthyst) << "Overwrite hyst memory at midpoint? " << overwriteMemoryMidPoint;
    LOG_TRACE(coeffcthyst) << "Overwrite hyst memory at int. points? " << overwriteMemoryIntPoints;
		bool onBoundary;
    
    LocPointMapped actualLPM;
    std::map<UInt, LocPointMapped>::iterator listIt;
    std::map<UInt, LocPointMapped> currentMap;
		// NEW: iterate over maps with already stored LPM
    if (midpointOnly == false) {
      currentMap = allLPMmap_;
    } else {
      currentMap = midpointLPMmap_;
    }
    
    for(listIt = currentMap.begin(); listIt!= currentMap.end(); ++listIt){
      LOG_TRACE(coeffcthyst) << "Storage idx: " << listIt->first;
      bool curOverwrite;
      bool forceMidpoint;
      
      if(midpointLPMmap_.find(listIt->first) != midpointLPMmap_.end()){
        // lpm is a midpoint > overwrite flag for midpoint needed
        LOG_TRACE(coeffcthyst) << "Midpoint found";
        curOverwrite = overwriteMemoryMidPoint;
        forceMidpoint = true;
      } else {
        LOG_TRACE(coeffcthyst) << "Integration point found";
        curOverwrite = overwriteMemoryIntPoints;
        forceMidpoint = false;
      }
      
      onBoundary = PreprocessLPM(listIt->second, actualLPM, operatorIdx, storageIdx, forceMidpoint);
      if(listIt->first != storageIdx){
        
        EXCEPTION("StorageIdx in map != storageIdx as obtained by PreprocessLPM!");
      }
      
      // get input
      input = RetrieveInputToHysteresisOperator(actualLPM, operatorIdx, storageIdx,onBoundary);
      // new: retrieveinput will store the actual solution directly to E_B_
      solution = E_B_[storageIdx];
      
      // then compute output
      output = CalcOutputOfHysteresisOperator(input, operatorIdx, storageIdx, forceMemoryLock, curOverwrite);
      
      // store at storage index
      if (setLastTS) {
        //std::cout << "Set values for integration point " << storageIdx << std::endl;
        E_B_lastTS_[storageIdx] = solution;
        //std::cout << "P_J_lastTS_[storageIdx] before setting: " << P_J_lastTS_[storageIdx] << std::endl;
        P_J_lastTS_[storageIdx] = output;
        //std::cout << "P_J_lastTS_[storageIdx] after setting: " << P_J_lastTS_[storageIdx] << std::endl;
        E_H_lastTS_[storageIdx] = input;
      } else {
        E_B_lastIt_[storageIdx] = solution;
        P_J_lastIt_[storageIdx] = output;
        E_H_lastIt_[storageIdx] = input;
      }
      deltaMatPrev_[storageIdx] = deltaMat_[storageIdx];
      deltaMatStrainPrev_[storageIdx] = deltaMatStrain_[storageIdx];
    }
		/*
     * 2.0 restore old state of evaluationPurpose
     */
		RUN_evaluationPurpose_ = evalPurposeBackup;    
	}
  
	void CoefFunctionHyst::GetScalar(Double& outputScalar, const LocPointMapped& lpm) {
		EXCEPTION("GetScalar not implemented for coefFncHyst");
	}
  
	void CoefFunctionHyst::GetVector(Vector<Double>& outputVector, const LocPointMapped& lpm) {
    
		// return current state of hyst operator
		UInt timeLevel = 0;
		outputVector = GetPrecomputedOutputOfHysteresisOperator(lpm, timeLevel, false);
	}
  
	void CoefFunctionHyst::GetTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm) {
    EXCEPTION("GetTensor not implemented for coefFncHyst");
	}
  
  void CoefFunctionHyst::SetFlag(std::string flagName, Integer intState) {
    
    /*
     * Flags that (may) change during runtime
     */
    if(flagName == "evaluationPurpose"){
      RUN_evaluationPurpose_ = intState;
    }
    else if(flagName == "allowBMP"){
      RUN_allowBMP_ = bool(intState);
    } 
    else if(flagName == "includeOldDeltaInPercent"){
      includeOldDelta_ = (Double) intState/100.0;
    }
    else if(flagName == "forceRemanence"){
      ForceRemanence();
    }
    else if(flagName == "lockHystOperator"){
      hystOperatorLocked_ = bool(intState);
    }
    else if(flagName == "forceCurrentTS"){
      forceCurrentTS_ = bool(intState);
    } 
    else if( (flagName == "SetTimeLevelRHSHyst") || (flagName == "SetTimeLevelRHSPol") ){
      timeLevel_rhsPol_ = intState;
    }
    else if( flagName == "SetTimeLevelRHSStrain"){
      timeLevel_rhsStrain_ = intState;
    }
    else if( flagName == "SetTimeLevelRHSCoupling"){
      timeLevel_rhsCoupling_ = intState;
    }
    else if(flagName == "SetTimeLevelMaterial"){
      timeLevel_Mat_ = intState;
    }
    else if( (flagName == "SetTimeLevelDeltaMat") || (flagName == "SetTimeLevelDeltaMatPol") ){
      timeLevel_deltaMatPol_ = intState;
      if(timeLevel_deltaMatPol_ == -2){
        deltaMatActive_ = false;
      } else {
        deltaMatActive_ = true;
      }
    }
    else if(flagName == "SetTimeLevelDeltaMatStrain"){
      timeLevel_deltaMatStrain_ = intState;
      if(timeLevel_deltaMatStrain_ == -2){
        deltaMatActive_ = false;
      } else {
        deltaMatActive_ = true;
      }
    }
    else if(flagName == "SetTimeLevelDeltaMatCoupling"){
      timeLevel_deltaMatCoupling_ = intState;
      if(timeLevel_deltaMatCoupling_ == -2){
        deltaMatActive_ = false;
      } else {
        deltaMatActive_ = true;
      }
    }
    else if(flagName == "resetReeval"){
      //std::cout << "Reset reeval flag" << std::endl;
      for (UInt i = 0; i < numStorageEntries_; i++) {
        requiresReeval_[i] = true;
        Si_requiresReeval_[i] = true;
        deltaMat_requiresReeval_[i] = true;
        deltaMatStrain_requiresReeval_[i] = true;
        rotatedCouplingTensor_requiresReeval_[i] = true;
        lastUsedTimeLevelForRotation_[i] = -10; 
      }
    }
    else if(flagName == "allowSettingOfMatForLocalInversion"){
      // reset flag for matrix inversion, i.e. the next time we call SetMatrixForInversion, we will overwrite the current value
      for (UInt i = 0; i < numStorageEntries_; i++) {
        matrixForInversionWasSet_[i] = false;
      }
      //      // to set matrices for inversion, flag has to be 1
      //      EvaluateHystOperatorsInt(1);
      
    } else if(flagName == "EvaluateHystOperators"){
      EvaluateHystOperatorsInt(intState);
    } else if(flagName == "SetPreviousHystValues"){
      SetPreviousHystVals(bool (intState));
    }
    
    /*
     * Flags that are to be set from xml file
     */
    else if(flagName == "deltaForm"){
      // ==0 > no deltaMat 
      // !=0 > deltaMat gets computed
      RUN_deltaForm_ = intState;
    }
    else if(flagName == "evaluationDepth"){
      /*
       * 1 -> standard evaluation
       * 2 -> extended evaluation
       * 3 -> full evaluation
       *
       * -> details: see .hh file
       */
      if (intState < 1) {
        XML_EvaluationDepth_ = 1;
      } else if (intState > 3) {
        XML_EvaluationDepth_ = 3;
      } else {
        XML_EvaluationDepth_ = intState;
      }
      
      // the evaluation depth is crucial for the number of storage entries
      // > init only after evalDepth is set!
      InitStorage();
      XMLParameterSet_ = true;
    }
    else if(flagName == "measurePerformance"){
      XML_performanceMeasurement_ = intState; 
      hyst_->setFlags(XML_performanceMeasurement_);
    }    
    else {
      std::stringstream except;
      except << "Flag " << flagName << " unknown to CoefFncHyst.";
      EXCEPTION(except.str())
    } 
  }
  
	bool CoefFunctionHyst::EvaluateAtMidpointOnly() {
		/*
     *  XML_EvaluationDepth_   >    1 (standard)    2 (extended)                          3 (full)
     *  RUN_evaluationPurpose_ v
     *    1 (assemble)          midpoint            integration point                     integration point
     *    2 (store lastIt)      midpoint            integration point + midpoint          integration point + midpoint
     *    3 (store lastTS)      midpoint            integration point + midpoint          integration point + midpoint
     *    4 (output)            midpoint            midpoint                              midpoint
     *
     */
		if (XML_EvaluationDepth_ == 1) {
			return true;
		} else if (XML_EvaluationDepth_ == 2) {
			// extended
			// NEW: even for extended evaluation, we have to store values at each integration point
			// (to get correct deltaMatrices)
			// but in difference to full evaluation, we lock the hysteresis operator for all points except the midpoint!
			if (RUN_evaluationPurpose_ == 4) {
				return true;
			} else {
				return false;
			}
		} else {
			// full
			if (RUN_evaluationPurpose_ <= 3) {
				return false;
			} else {
				return true;
			}
		}
	}
  
	bool CoefFunctionHyst::OverwriteHystMemory() {
		/*
     *  RUN_evaluationPurpose_
     *    1 (assemble)                false
     *    2 (store lastIT)            false
     *    3 (store lastTS)            true
     *    4 (output)                  false
     *
     *  Motivation:
     *    During the nonlinear solution process, we oftentimes overestimate
     *    the input to the hysteresis operator and have to step back.
     *    These stepping (especially during linesearch) lead to permanent changes
     *    of the hysteresis memory which are (normally) not wanted.
     *    Therefore, the hysteresis memory should be locked, i.e. all changes
     *    are done to a temporal copy of the current state.
     *    This leads to a nonlinear but nonhysteretic stepping.
     *
     *    At the end of each timestep (where we assume to have found a valid solution)
     *    we finally unlock the hysteresis memory and evaluate the hysteresis operator
     *    once more at the last computed solution to bring the operator to the correct
     *    state.
     *
     *    Conclusion:
     *      during assemble and storage (of last iteration values) and during
     *      output computation, we lock the hysteresis memory
     *      (output computation would not be the problem as we would simply
     *      apply the same input as during the storage of the lastTS value)
     *
     *      for storing the lastTS values), we unlock it
     *
     *
     */
		if (RUN_evaluationPurpose_ == 3) {
			return true;
		} else {
			return false;
		}
    
	}
  
  void CoefFunctionHyst::CreatePeriodicTestSignal(std::string name, Double amplitudeScaling, Double numPeriods, UInt stepsPerPeriod, Vector<Double>& xVals, Vector<Double>& yVals){
    
    UInt totalSteps = UInt(stepsPerPeriod*numPeriods);
    xVals.Resize(totalSteps);
    yVals.Resize(totalSteps);
    xVals.Init();
    yVals.Init();
    
    if(name == "DecreasingSawtooth"){
      /*
       * Testcase: Put material into full saturation in y-direction; from that state on
       * apply a decreasing triangular signal in x-direction; for initial period hold
       * signla in y-direction; during first real period decrease to 1/2 of value;
       * then to 1/4 then to 1/8 and finally to 0
       */
      UInt numStepsFirstIncrease = stepsPerPeriod;
      Double incr = POL_operatorParams_.inputSat_/((Double) numStepsFirstIncrease-1);
      
      yVals[0] = -POL_operatorParams_.inputSat_;
      for(UInt i = 1; i < numStepsFirstIncrease; i++){
        xVals[i] = xVals[i-1] + incr; 
        yVals[i] = -POL_operatorParams_.inputSat_;
      }
      
      UInt cnt = numStepsFirstIncrease;
      UInt remainingSteps = totalSteps - numStepsFirstIncrease;
      UInt numFullPeriods = UInt(remainingSteps/stepsPerPeriod);
      remainingSteps = remainingSteps%stepsPerPeriod;
      
      incr = 2*POL_operatorParams_.inputSat_/((Double)(stepsPerPeriod-1));
      Double sign = -1.0;
      Double decrFactor = 1.0 - 1.0/((Double) numFullPeriods);
      
      for(UInt j = 0; j < numFullPeriods; j++){
        for(UInt i = 0; i < stepsPerPeriod; i++){
          xVals[cnt] = xVals[cnt-1] + sign*incr;
          
          if(j==0){
            yVals[cnt] = -POL_operatorParams_.inputSat_/2.0;
          } else if(j==1){
            yVals[cnt] = -POL_operatorParams_.inputSat_/4.0;
          } else if(j==2){
            yVals[cnt] = -POL_operatorParams_.inputSat_/8.0;
          }
          cnt++;
        }
        sign = -1.0*sign;
        incr = incr*decrFactor;
      }
      for(UInt i = 0; i < remainingSteps; i++){
        xVals[cnt] = xVals[cnt-1] + sign*incr;
        cnt++;
      }
    } else if(name == "DecreasingSine"){
      Double decrease = 1.0/totalSteps;
      
      for(UInt i = 0; i < totalSteps; i++){
        //xIn[i][0] = sin( (2*M_PI*i)/stepsPerPeriod );
        xVals[i] = POL_operatorParams_.inputSat_*(1.0 - decrease*i) * sin( (2*M_PI*i)/stepsPerPeriod );
        
        if(i >= 3*stepsPerPeriod){
          yVals[i] = -POL_operatorParams_.inputSat_/8.0;
        } else if(i >= 2*stepsPerPeriod){
          yVals[i] = -POL_operatorParams_.inputSat_/4.0;
        } else if(i >= stepsPerPeriod){
          yVals[i] = -POL_operatorParams_.inputSat_/2.0;
        } else if(i >= 0){
          yVals[i] = -POL_operatorParams_.inputSat_;
        }
      }
    } else if(name == "DecreasingRotation"){
      Double decrease = 1.0/totalSteps;
      
      for(UInt i = 0; i < totalSteps; i++){
        xVals[i] = POL_operatorParams_.inputSat_*(1.0 - decrease*i)*sin( (2*M_PI*i)/stepsPerPeriod );
        yVals[i] = POL_operatorParams_.inputSat_*(1.0 - decrease*i)*cos( (2*M_PI*i)/stepsPerPeriod );
      }
    } else if(name == "IncreasingRotation"){
      Double increase = 1.0/totalSteps;
      
      for(UInt i = 0; i < totalSteps; i++){
        xVals[i] = POL_operatorParams_.inputSat_*increase*i*sin( (2*M_PI*i)/stepsPerPeriod );
        yVals[i] = POL_operatorParams_.inputSat_*increase*i*cos( (2*M_PI*i)/stepsPerPeriod );
      }
    } else if(name == "Sine"){
      for(UInt i = 0; i < totalSteps; i++){
        xVals[i] = POL_operatorParams_.inputSat_*sin( (2*M_PI*i)/stepsPerPeriod );
        yVals[i] = POL_operatorParams_.inputSat_*0.25;
      }
    } else if(name == "Rotation"){
      for(UInt i = 0; i < totalSteps; i++){
        xVals[i] = POL_operatorParams_.inputSat_*sin( (2*M_PI*i)/stepsPerPeriod );
        yVals[i] = POL_operatorParams_.inputSat_*cos( (2*M_PI*i)/stepsPerPeriod );
      }
    } else if(name == "Forc"){
      Double decrease = 1.0/totalSteps;
      
      for(UInt i = 0; i < totalSteps; i++){
        //xIn[i][0] = sin( (2*M_PI*i)/stepsPerPeriod );
        xVals[i] = 1*POL_operatorParams_.inputSat_*(1.0 - decrease*i) * sin( (2*M_PI*i)/stepsPerPeriod ) - 1*POL_operatorParams_.inputSat_*decrease*i;
        yVals[i] = 0.3*POL_operatorParams_.inputSat_*(1.0 - decrease*i) * sin( 3*(2*M_PI*i)/stepsPerPeriod ) - 0.3*POL_operatorParams_.inputSat_*decrease*i;
      }
    } else {
      EXCEPTION("Signal not implemented yet");
    }
    
    xVals.ScalarMult(amplitudeScaling);
    yVals.ScalarMult(amplitudeScaling);
  }
  
  void CoefFunctionHyst::CreateNonPeriodicTestSignal(std::string name, Double amplitudeScaling, UInt numberOfSteps, Vector<Double>& xVals, Vector<Double>& yVals){
    
    if(name == "SelfDesigned"){
      UInt steps1,steps2,steps3,steps4,steps5,steps6,steps7;
      steps1=50;
      steps2=50;
      steps3=50;
      steps4=50;
      steps5=50;
      steps6=50;
      steps7=10;
      UInt totalSteps=steps1+steps2+steps3+steps4+steps5+steps6+steps7;
      
      xVals.Resize(totalSteps);
      xVals.Init();
      yVals.Resize(totalSteps);
      yVals.Init();
      Vector<Double> xScal = Vector<Double>(totalSteps);
      xScal.Init();
      Double delta;
      
      // 1. expontential increase from [10^-15 to 10^0[
      delta = 15.0/steps1;
      
      for(UInt i = 0; i < steps1; i++){
        xScal[i] = std::pow(10,-15.0+i*delta);
      }
      // 2. cos decrease towards -0.5[
      delta = -std::acos(-0.5)/steps2;
      
      for(UInt i = 0; i < steps2; i++){
        xScal[steps1+i] = std::cos(i*delta);
      }
      // 3. linear increase to 0.4
      delta = 0.9/steps3;
      
      for(UInt i = 0; i < steps3; i++){
        xScal[steps1+steps2+i] = -0.5+i*delta;
      }
      // 4. 1/x decrease towards -0.3
      Double Xstart = 1.0/(0.7 + 1.0/steps4);
      delta = (steps4 - Xstart)/steps4;
      
      for(UInt i = 0; i < steps4; i++){
        xScal[steps1+steps2+steps3+i] = 1.0/(Xstart+i*delta)-1.0/steps4-0.3;
      }
      // 5. logarithmic increase towards 0.2
      Xstart = std::exp(-3);
      Double Xend = std::exp(2);
      delta = (Xend-Xstart)/steps5;
      
      for(UInt i = 0; i < steps5; i++){
        xScal[steps1+steps2+steps3+steps4+i] = 0.1*std::log(Xstart+i*delta);
      }
      // 6. x^3 increase to 0.8
      Double denom = std::pow(steps6-1,3);
      for(UInt i = 0; i < steps6; i++){
        xScal[steps1+steps2+steps3+steps4+steps5+i] = 0.6*std::pow(i,3)/denom + 0.2;
      }
      
      // 7. hold signal at value 0.8
      for(UInt i = 0; i < steps7; i++){
        xScal[steps1+steps2+steps3+steps4+steps5+steps6+i] = 0.8;
      }
      
      //        std::ofstream testOutput;
      //        testOutput.open("GeneratedInput");
      //        for(UInt i = 0; i < totalSteps; i++){
      //          testOutput << i << " " << xScal[i] << std::endl;
      //        }
      //        testOutput.close();
      Vector<Double> dir = Vector<Double>(dim_);
      dir[0] = 1.0;
      dir[1] = 0.0;
      for(UInt i = 0; i < totalSteps; i++){
        dir[0] = (Double) (totalSteps-i)/totalSteps;
        dir[1] = (Double) i/totalSteps;
        
        xVals[i] = POL_operatorParams_.inputSat_*xScal[i]*dir[0];
        yVals[i] = POL_operatorParams_.inputSat_*xScal[i]*dir[1];
      }
    } else if(name == "RemDrop"){
      // saturate material with one step, then increase input signal perpendicular up to saturation
      // similar to "SatX-RemX-SatY" 
      xVals.Resize(numberOfSteps+1);
      yVals.Resize(numberOfSteps+1);
      xVals.Init();
      yVals.Init();
      
      Double maxAmplitude = POL_operatorParams_.inputSat_;
      xVals[0] = maxAmplitude;
      
      Double increase = 1.0/Double(numberOfSteps-1);
      for(UInt i = 0; i < numberOfSteps; i++){
        yVals[i+1] = maxAmplitude*i*increase;
      }   
    } else if(name == "SatX-RemX-SatY"){
      if(numberOfSteps < 8){
        numberOfSteps = 8;
      }
      if(numberOfSteps%4 != 0){
        numberOfSteps += numberOfSteps%4; 
      }
      
      Double increase = 4.0/numberOfSteps;
      
      xVals.Resize(numberOfSteps);
      yVals.Resize(numberOfSteps);
      xVals.Init();
      yVals.Init();
      
      Double maxAmplitude = POL_operatorParams_.inputSat_;
      
      for(UInt i = 0; i < numberOfSteps/4; i++){
        xVals[i] = maxAmplitude*i*increase;
      }
      for(UInt i = numberOfSteps/4; i < 2*numberOfSteps/4; i++){          
        xVals[i] = maxAmplitude - maxAmplitude*(i-numberOfSteps/4)*increase;
      }
      
      for(UInt i = 2*numberOfSteps/4; i < 3*numberOfSteps/4; i++){
        // linearly increase in y up to n*saturation
        yVals[i] = maxAmplitude*(i-2*numberOfSteps/4)*increase;
      }
      for(UInt i = 3*numberOfSteps/4; i < numberOfSteps; i++){
        // linearly decrease to 0
        yVals[i] = maxAmplitude - maxAmplitude*(i-3*numberOfSteps/4)*increase;
      }
    } else {
      EXCEPTION("Signal not implemented yet");
    }
    
    xVals.ScalarMult(amplitudeScaling);
    yVals.ScalarMult(amplitudeScaling);
    
  }
  
  void CoefFunctionHyst::WriteSignalToFile(std::string name, Vector<Double> xVals, Vector<Double> yVals){
    std::ofstream output;
    output.open(name);
    UInt totalSteps = xVals.GetSize();
    output << "# Step    x    y" << std::endl;
    for(UInt i = 0; i < totalSteps; i++){
      output << i << " " << xVals[i] << " " << yVals[i] << std::endl;
    }
    output.close();
  }
  
  void CoefFunctionHyst::TestHystOperatorWithSignal(std::string name, Vector<Double> xVals, Vector<Double> yVals, 
          bool testInversion, bool printStatistics, bool writeResultsToFile, 
          bool measurePerformance, std::string commonPerformanceFile, bool test1D, bool outputIrrStrains){
    
		/*
     * 0. Declare variables (there are alot)
     */
		std::ofstream statistics;
		std::ofstream results_x;
    std::ofstream results_p;
    std::ofstream results_s;
    std::ofstream results_xp;
    std::ofstream results_xps;
    std::ofstream results_y;
    std::ofstream angularResults_x;
    std::ofstream angularResults_p;
    std::ofstream angularResults_y;
		std::ofstream performance;
    std::ofstream orientationTowardsExcitation_p;
    std::ofstream commonPerfStream;
		
    std::string basedir = "./history/";
    
		std::stringstream statistics_name;
		statistics_name << basedir << name << "_statistics";
		
		std::stringstream results_name_x;
		results_name_x << basedir << name << "_results_x";		   
    
    std::stringstream results_name_p;
		results_name_p << basedir << name << "_results_p";		
    
    std::stringstream results_name_s;
		results_name_s << basedir << name << "_results_s";		
    
    std::stringstream results_name_y;
		results_name_y << basedir << name << "_results_y";		
		
    std::stringstream results_name_xp;
		results_name_xp << basedir << name << "_results_xp";	
    
    std::stringstream results_name_xps;
		results_name_xps << basedir << name << "_results_xps";	
    
    std::stringstream angularResults_name_x;
		angularResults_name_x << basedir << name << "_angularResults_x";		   
    
    std::stringstream angularResults_name_p;
		angularResults_name_p << basedir << name << "_angularResults_p";		
    
    std::stringstream angularResults_name_y;
		angularResults_name_y << basedir << name << "_angularResults_y";		
    
		std::stringstream performance_name;
		performance_name << basedir << name << "_performance";		  
		
    std::stringstream orientation_name;
		orientation_name << basedir << name << "_projectedCoords_p";	

    std::stringstream commonPerfStream_name;
		commonPerfStream_name << basedir << commonPerformanceFile;	
    
		if(writeResultsToFile){
			results_x.open(results_name_x.str());
      results_xp.open(results_name_xp.str());
      results_p.open(results_name_p.str());
      results_y.open(results_name_y.str());
      angularResults_x.open(angularResults_name_x.str());
      angularResults_p.open(angularResults_name_p.str());
      angularResults_y.open(angularResults_name_y.str());
      orientationTowardsExcitation_p.open(orientation_name.str());
		}
    
		if(printStatistics){
			statistics.open(statistics_name.str());
		}
		if(measurePerformance){
			performance.open(performance_name.str());
      if(commonPerformanceFile != "---"){
        commonPerfStream.open(commonPerfStream_name.str(),std::ios_base::app);
      }
		}
		
		UInt totalSteps = xVals.GetSize();
		UInt numFails = 0;
    UInt numHalfFails = 0; // fail error crit but pass residual crit
		int successFlagForward = -1;
		// forward:
		// -1 = fail
		//  0 = reuse value
		//  1 = anhyst only
		//  2 = eval on permanent storage
		//  3 = eval on temporal storage
		
		int successFlagBackward = -1;
		// backward:
		// -1 = fail
		//  0 = reuse value
		//  1 = anhyst only
		//  2 = bisection
		//  3-6 only for vector implementation using Levenberg Marquardt
		//  3 = reamnence
		//  4 = passed dut to error tolerance
		//  5 = passed due to tolerance wrt x
		//  6 = passed due to tolerance wrt y
		
		Vector<Double> xIn = Vector<Double>(dim_);
		Vector<Double> hIn = Vector<Double>(dim_);
		Vector<Double> yIn = Vector<Double>(dim_);
		
    Vector<Double> xInPrev = Vector<Double>(dim_);
		Vector<Double> xPrev = Vector<Double>(dim_);
		Vector<Double> hPrev = Vector<Double>(dim_);
		Vector<Double> yPrev = Vector<Double>(dim_);
		
		Vector<Double> xOut = Vector<Double>(dim_);
		Vector<Double> hOut = Vector<Double>(dim_);
    Vector<Double> hOutForStrains = Vector<Double>(dim_);
		Vector<Double> yOut = Vector<Double>(dim_);
		Vector<Double> hOutOld = Vector<Double>(dim_);
    Vector<Double> hOutOldForStrains = Vector<Double>(dim_);
    
    Vector<Double> xInBak = Vector<Double>(dim_);
		Vector<Double> xErr = Vector<Double>(dim_);
		Vector<Double> yErr = Vector<Double>(dim_);
		
		std::map< UInt, Vector<Double> > failedTests_xIn;
		std::map< UInt, Vector<Double> > failedTests_xOut;
		std::map< UInt, Vector<Double> > failedTests_yIn;
		std::map< UInt, Vector<Double> > failedTests_yOut;
		std::map< UInt, std::pair< bool, Double > > failedTests_wrtX;
		std::map< UInt, std::pair< bool, Double > > failedTests_wrtY;
		
		Vector<int> successCodeVectorForward = Vector<int>(totalSteps);
		Vector<int> successCodeVectorBackward = Vector<int>(totalSteps);
		
		Double minAlpha = 0;
		Double maxAlpha = 0;
		Double avgAlpha = 0;
		Double totalminAlpha = 1e10;
		Double totalmaxAlpha = -1e10;
		Double totalavgAlpha = 0;
		
		UInt numberOfLMIterations = 0;
		UInt numberOfLinesearchIterations = 0;
		UInt maxNumberOfLinesearchIterations = 0;
		
		// counter for forward evaluation
		UInt forwardFails = 0;
		UInt forwardReused = 0;
		UInt forwardAnhystOnly = 0;
		UInt forwardOverwrite = 0;
		UInt forwardTMP = 0;
		
		// counter for backward evaluation / inversion
		UInt LMFails = 0;
		UInt totalReused = 0;
		UInt totalAnhystOnly = 0;
		UInt totalBisection = 0;
		UInt totalRemanence = 0;
		UInt totalPassedErrorTol = 0;
		UInt totalPassedResTolX = 0;
		UInt totalPassedResTolY = 0;
		UInt totalNumberOfLMIterations = 0;
		UInt totalNumberOfLinesearchIterations = 0;
		UInt absmaxNumberOfLinesearchIterations = 0;
		
		Double startTime = 0.0;
    Double endTime = 0.0;
    Double evalTime = 0.0;
		Double backwardMaxEvalTime = -1.0;
		Double backwardTotalEvalTime = 0.0;
		Double forwardMaxEvalTime = -1.0;
		Double forwardTotalEvalTime = 0.0;
		
		UInt forwardEvalCounter = 0;
		UInt backwardEvalCounter = 0;
		Timer* forwardTimer;
		Timer* backwardTimer;
		
		if(printStatistics){
			statistics << "+++ STATISTICS +++" << std::endl;
			statistics << "TEST: " << name << std::endl;
			statistics << "MODEL: " << POL_operatorParams_.methodName_ << std::endl;
      statistics << "Error Criteria: " << std::endl;
			statistics << "- Residual wrt input (=x): " << LM_inversion_.tolH << std::endl;
      statistics << "- Residual wrt output (=y): " << LM_inversion_.tolB << std::endl;
		}
		
		/*
     * 1. Create temporal hyst operator; this should be done for each test to ensure 
     *		that operator is in initial state
     */
		bool isVirgin = true;
		bool vector = true;
		bool debugOut = false;
		
		if (POL_operatorParams_.methodName_ == "scalarPreisach") {
			
			POL_useExtension_ = false;
      //			int useExtensionInt;
      //			material_->GetScalar(useExtensionInt, SCALPREISACH_USE_EXT);
			if(POL_useExtension_){
				EXCEPTION("Extension not implemented for tests; remove completely as not working");
//				POL_useExtension_ = true;
//				
//				material_->GetScalar(POL_operatorParams_.rotResistance_, ROT_RESISTANCE, Global::REAL);
//				material_->GetScalar(POL_operatorParams_.angularDistance_, ANG_DISTANCE, Global::REAL);
//				
//				hystTMP = new ExtendedPreisach(1, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_, POL_weightParams_.weightTensor_, 
//                POL_operatorParams_.rotResistance_, POL_operatorParams_.angularDistance_, dim_, isVirgin, POL_weightParams_.anhysteretic_a_, 
//                POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
//				
//				// set initial direction
//				Vector<Double> initialInput = Vector<Double>(dim_);
//				POL_eps_mu_SmallSignal_.Mult(POL_operatorParams_.fixDirection_,initialInput);
//				initialInput.ScalarMult(POL_operatorParams_.inputSat_);
//				initialInput.Add(POL_operatorParams_.outputSat_,POL_operatorParams_.fixDirection_);
//				// initialInput = (POL_operatorParams_.outputSat_*Identity + POL_operatorParams_.inputSat_*POL_eps_mu_SmallSignal_)*POL_operatorParams_.fixDirection_
//				// > flux with value just at saturation
//				hystTMP->UpdateRotationState(initialInput,POL_eps_mu_SmallSignal_,0);
//				
			} else {
				vector = false;
        bool ignoreAnhystPart = false;
				hystTMP = new Preisach(1, POL_operatorParams_, POL_weightParams_, isVirgin, ignoreAnhystPart);
//                (1, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_, 
//                POL_weightParams_.weightTensor_, isVirgin, 
//                POL_weightParams_.anhysteretic_a_, POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,
//                POL_weightParams_.anhystOnly_);
      }
		} else if (POL_operatorParams_.methodName_ == "vectorPreisach_Sutor") {
			if(printStatistics){
				if (POL_operatorParams_.evalVersion_ == 1) {
					statistics << "- classical model, list based implementation" << std::endl;
				} else if (POL_operatorParams_.evalVersion_ == 2) {
					statistics << "- revised model, list based implementation" << std::endl;
				} else if (POL_operatorParams_.evalVersion_ == 10) {
					statistics << "- classical model, matrix based implementation" << std::endl;
				} else if (POL_operatorParams_.evalVersion_ == 20) {
					statistics << "- revised model, matrix based implementation" << std::endl;
				}
				
				statistics << "with rotational resistance = " << POL_operatorParams_.rotResistance_ << std::endl;
				if(POL_operatorParams_.isClassical_ == false){
					statistics << "and angular distance = " << POL_operatorParams_.angularDistance_ << std::endl;
				}
			}
			
			if (POL_operatorParams_.evalVersion_ == 1) {
				POL_operatorParams_.isClassical_ = true; // original vector preisach model -> sutor2012
				
				hystTMP = new VectorPreisachSutor_ListApproach(1, POL_operatorParams_, POL_weightParams_, dim_, isVirgin);               
//                (1, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_,
//                POL_weightParams_.weightTensor_, POL_operatorParams_.rotResistance_, dim_, isVirgin,
//                POL_operatorParams_.isClassical_, POL_operatorParams_.scaleUpToSaturation_,
//                POL_operatorParams_.angularDistance_,POL_operatorParams_.angularResolution_,POL_weightParams_.anhysteretic_a_, 
//                POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
			} else if (POL_operatorParams_.evalVersion_ == 2) {
				POL_operatorParams_.isClassical_ = false; // revised vector preisach model -> sutor2015
				
				hystTMP = new VectorPreisachSutor_ListApproach(1, POL_operatorParams_, POL_weightParams_, dim_, isVirgin);
//                (1, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_,
//                POL_weightParams_.weightTensor_, POL_operatorParams_.rotResistance_, dim_, isVirgin,
//                POL_operatorParams_.isClassical_, POL_operatorParams_.scaleUpToSaturation_,
//                POL_operatorParams_.angularDistance_,POL_operatorParams_.angularResolution_,POL_weightParams_.anhysteretic_a_, 
//                POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
			} else if (POL_operatorParams_.evalVersion_ == 10) {
				POL_operatorParams_.isClassical_ = true; // original vector preisach model -> sutor2015; matrix based implementation
				
				hystTMP = new VectorPreisachSutor_MatrixApproach(1, POL_operatorParams_, POL_weightParams_, dim_, isVirgin);
//                (1, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_,
//                POL_weightParams_.weightTensor_, POL_operatorParams_.rotResistance_, dim_, isVirgin,
//                POL_operatorParams_.isClassical_, POL_operatorParams_.scaleUpToSaturation_,
//                POL_operatorParams_.angularDistance_,POL_operatorParams_.angularResolution_,POL_weightParams_.anhysteretic_a_, 
//                POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
			} else if (POL_operatorParams_.evalVersion_ == 20) {
				POL_operatorParams_.isClassical_ = false; // revised vector preisach model -> sutor2015; matrix based implementation
				
				hystTMP = new VectorPreisachSutor_MatrixApproach(1, POL_operatorParams_, POL_weightParams_, dim_, isVirgin);
//                (1, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_,
//                POL_weightParams_.weightTensor_, POL_operatorParams_.rotResistance_, dim_, isVirgin,
//                POL_operatorParams_.isClassical_, POL_operatorParams_.scaleUpToSaturation_,
//                POL_operatorParams_.angularDistance_,POL_operatorParams_.angularResolution_,POL_weightParams_.anhysteretic_a_, 
//                POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
			} else {
				EXCEPTION("POL_operatorParams_.evalVersion_ has to be one of the following: \n "
                "1: classical vector model (sutor2012) \n"
                "2: revised vector model (sutor2015) [DEFAULT] \n"
                "10: classical vector model (sutor2012) - Matrix implementation, only for reference \n"
                "20: revised vector model (sutor2015) - Matrix implementation, only for reference \n")
			}
      hystTMP->SetParamsForInversion(LM_inversion_); 
      
//			hystTMP->SetParamsForInversion(LM_inversion_.inversionMethod, LM_inversion_.maxNumIts, LM_inversion_.maxNumLSIts,
//                LM_inversion_.tolH, LM_inversion_.tolB, 
//                LM_inversion_.jacRes, LM_inversion_.alphaLSStart,LM_inversion_.alphaLSMin,LM_inversion_.alphaLSMax,
//                LM_inversion_.stopLineSearchAtLocalMin,
//                POL_operatorParams_.angularClipping_);   
      
		} else if (POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") {			
      // basically a scalar model in multiple directions
      // isotropic case: all scalar models are equal (same weights etc)
      // anisotropic case: each model different; choice of directions matters; weights are harder to obtain
      int isIsotropic = 1;
      material_->GetScalar(isIsotropic, PREISACH_MAYERGOYZ_ISOTROPIC);
      if( (dim_ != 2) || (isIsotropic == 0)){
        EXCEPTION("Mayergoyz vector model currently only implemented for 2d isotropic materials");
      }
      
      statistics << "with num directions = " << POL_operatorParams_.numDirections_ << std::endl;
      statistics << "and output clipping mode = " << POL_operatorParams_.outputClipping_ << std::endl;
      
      /*
       * IMPORTANT REMARK:
       *  > although the Mayergoyz model is based on the scalar models 
       *     we are not allowed to directly apply the Preisach parameter for
       *     the scalar case (i.e. the weights, the anhyst parameter and so on)
       *  > make sure that the passed parameter are already transformed correctly
       *      > see constructor above
       */
      hystTMP = new VectorPreisachMayergoyz(1, POL_operatorParams_, POL_weightParams_, dim_, isVirgin);
//              (1, POL_operatorParams_.numDirections_, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_, 
//              POL_weightParams_.weightTensor_,dim_,isVirgin,
//              POL_weightParams_.anhysteretic_a_, POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,
//              POL_weightParams_.anhystOnly_,POL_operatorParams_.outputClipping_);
			
      hystTMP->SetParamsForInversion(LM_inversion_);
      
//			hystTMP->SetParamsForInversion(LM_inversion_.inversionMethod, LM_inversion_.maxNumIts, LM_inversion_.maxNumLSIts,
//                LM_inversion_.tolH, LM_inversion_.tolB, 
//                LM_inversion_.jacRes, LM_inversion_.alphaLSStart,LM_inversion_.alphaLSMin,LM_inversion_.alphaLSMax,
//                LM_inversion_.stopLineSearchAtLocalMin,
//                POL_operatorParams_.angularClipping_);   
		} else {
			EXCEPTION("Invalid model selected for inversion test");
		}
		
    if(CouplingParams_.ownHystOperator_ && CouplingParams_.couplingDefined_inMatFile_){
      if (STRAIN_operatorParams_.methodName_ == "scalarPreisach") {
        bool ignoreAnhystPart = false;
        hystStrainTMP = new Preisach(1, STRAIN_operatorParams_, STRAIN_weightParams_, isVirgin, ignoreAnhystPart);
//                (1, STRAIN_operatorParams_.inputSat_, STRAIN_operatorParams_.outputSat_, 
//                STRAIN_weightParams_.weightTensor_, isVirgin, 
//                STRAIN_weightParams_.anhysteretic_a_, STRAIN_weightParams_.anhysteretic_b_, STRAIN_weightParams_.anhysteretic_c_,
//                STRAIN_weightParams_.anhystOnly_);
        
      } else if (STRAIN_operatorParams_.methodName_ == "vectorPreisach_Sutor") {
        if (STRAIN_operatorParams_.evalVersion_ == 1) {
          STRAIN_operatorParams_.isClassical_ = true; // original vector preisach model -> sutor2012
          
          hystStrainTMP = new VectorPreisachSutor_ListApproach(1, STRAIN_operatorParams_, STRAIN_weightParams_, dim_, isVirgin);
//                  (1, STRAIN_operatorParams_.inputSat_, STRAIN_operatorParams_.outputSat_,
//                  STRAIN_weightParams_.weightTensor_, STRAIN_operatorParams_.rotResistance_, dim_, isVirgin,
//                  STRAIN_operatorParams_.isClassical_, STRAIN_operatorParams_.scaleUpToSaturation_,
//                  STRAIN_operatorParams_.angularDistance_,STRAIN_operatorParams_.angularResolution_,
//                  STRAIN_weightParams_.anhysteretic_a_, 
//                  STRAIN_weightParams_.anhysteretic_b_, STRAIN_weightParams_.anhysteretic_c_,STRAIN_weightParams_.anhystOnly_);
        } else if (STRAIN_operatorParams_.evalVersion_ == 2) {
          STRAIN_operatorParams_.isClassical_ = false; // revised vector preisach model -> sutor2015
          
          hystStrainTMP = new VectorPreisachSutor_ListApproach(1, STRAIN_operatorParams_, STRAIN_weightParams_, dim_, isVirgin);
//                  (1, STRAIN_operatorParams_.inputSat_, STRAIN_operatorParams_.outputSat_,
//                  STRAIN_weightParams_.weightTensor_, STRAIN_operatorParams_.rotResistance_, dim_, isVirgin,
//                  STRAIN_operatorParams_.isClassical_, STRAIN_operatorParams_.scaleUpToSaturation_,
//                  STRAIN_operatorParams_.angularDistance_,STRAIN_operatorParams_.angularResolution_,
//                  STRAIN_weightParams_.anhysteretic_a_, 
//                  STRAIN_weightParams_.anhysteretic_b_, STRAIN_weightParams_.anhysteretic_c_,STRAIN_weightParams_.anhystOnly_);
        } else if (STRAIN_operatorParams_.evalVersion_ == 10) {
          STRAIN_operatorParams_.isClassical_ = true; // original vector preisach model -> sutor2015; matrix based implementation
          
          hystStrainTMP = new VectorPreisachSutor_MatrixApproach(1, STRAIN_operatorParams_, STRAIN_weightParams_, dim_, isVirgin);
//                  (1, STRAIN_operatorParams_.inputSat_, STRAIN_operatorParams_.outputSat_,
//                  STRAIN_weightParams_.weightTensor_, STRAIN_operatorParams_.rotResistance_, dim_, isVirgin,
//                  STRAIN_operatorParams_.isClassical_, STRAIN_operatorParams_.scaleUpToSaturation_,
//                  STRAIN_operatorParams_.angularDistance_,STRAIN_operatorParams_.angularResolution_,
//                  STRAIN_weightParams_.anhysteretic_a_, 
//                  STRAIN_weightParams_.anhysteretic_b_, STRAIN_weightParams_.anhysteretic_c_,STRAIN_weightParams_.anhystOnly_);
        } else if (STRAIN_operatorParams_.evalVersion_ == 20) {
          STRAIN_operatorParams_.isClassical_ = false; // revised vector preisach model -> sutor2015; matrix based implementation
          
          hystStrainTMP = new VectorPreisachSutor_MatrixApproach(1, STRAIN_operatorParams_, STRAIN_weightParams_, dim_, isVirgin);
//                  (1, STRAIN_operatorParams_.inputSat_, STRAIN_operatorParams_.outputSat_,
//                  STRAIN_weightParams_.weightTensor_, STRAIN_operatorParams_.rotResistance_, dim_, isVirgin,
//                  STRAIN_operatorParams_.isClassical_, STRAIN_operatorParams_.scaleUpToSaturation_,
//                  STRAIN_operatorParams_.angularDistance_,STRAIN_operatorParams_.angularResolution_,
//                  STRAIN_weightParams_.anhysteretic_a_, 
//                  STRAIN_weightParams_.anhysteretic_b_, STRAIN_weightParams_.anhysteretic_c_,STRAIN_weightParams_.anhystOnly_);
        }
      } else if (STRAIN_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") {
        
          hystStrainTMP = new VectorPreisachMayergoyz(1, STRAIN_operatorParams_, STRAIN_weightParams_, dim_, isVirgin);
//                  (1, STRAIN_operatorParams_.numDirections_, 
//                STRAIN_operatorParams_.inputSat_, STRAIN_operatorParams_.outputSat_, 
//                STRAIN_weightParams_.weightTensor_,dim_,isVirgin,
//                STRAIN_weightParams_.anhysteretic_a_, STRAIN_weightParams_.anhysteretic_b_, STRAIN_weightParams_.anhysteretic_c_,
//                STRAIN_weightParams_.anhystOnly_,STRAIN_operatorParams_.outputClipping_);
                
      }
    }
    
    
    if(outputIrrStrains){
      results_s.open(results_name_s.str());
      results_s << "#Number S_irr_xx S_irr_yy S_irr_xy" << std::endl;
      
      results_xps.open(results_name_xps.str());
      results_xps << "#Number xIn yIn Px Py S_irr_xx S_irr_yy S_irr_xy" << std::endl;
    }
    
		if(printStatistics){
			statistics << "PARAMETER: " << std::endl;
			statistics << "- xSAT: " << POL_operatorParams_.inputSat_ << std::endl;
			statistics << "- pSAT: " << POL_operatorParams_.outputSat_ << std::endl;
			statistics << "- anhyst a: " << POL_weightParams_.anhysteretic_a_ << std::endl;
			statistics << "- anhyst b: " << POL_weightParams_.anhysteretic_b_ << std::endl;
			statistics << "- anhyst c: " << POL_weightParams_.anhysteretic_c_ << std::endl;
			statistics << "- only anhyst? " << POL_weightParams_.anhystOnly_ << std::endl;
			statistics << "PREISACH WEIGHTS: " << std::endl;
      
      if(POL_weightParams_.weightType_ == "Constant"){
        statistics << "> constant = " << POL_weightParams_.constWeight_ << std::endl;
      } else if(POL_weightParams_.weightType_ == "muDat"){
        statistics << "> muDat with " << std::endl;
        statistics << "- A = " << POL_weightParams_.muDat_A_ << std::endl;
        statistics << "- sigma = " << POL_weightParams_.muDat_sigma1_ << std::endl;
        statistics << "- h = " << POL_weightParams_.muDat_h1_ << std::endl;
        statistics << "- eta = " << POL_weightParams_.muDat_eta_ << std::endl;
      } else if(POL_weightParams_.weightType_ == "muDat"){
        statistics << "> extended muDat with " << std::endl;
        statistics << "- A = " << POL_weightParams_.muDat_A_ << std::endl;
        statistics << "- sigma1 = " << POL_weightParams_.muDat_sigma1_ << std::endl;
        statistics << "- sigma2 = " << POL_weightParams_.muDat_sigma2_ << std::endl;
        statistics << "- h1 = " << POL_weightParams_.muDat_h1_ << std::endl;
        statistics << "- h2 = " << POL_weightParams_.muDat_h2_ << std::endl;
        statistics << "- eta = " << POL_weightParams_.muDat_eta_ << std::endl;
      } else if(POL_weightParams_.weightType_ == "givenTensor"){
        statistics << "> given tensor = " << std::endl;
        statistics << POL_weightParams_.weightTensor_.ToString() << std::endl;
      }
      if (POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") {
        if(POL_weightsAlreadyAdapted_ == 1){
          statistics << "> weights directly used for Mayergoyz model" << std::endl;
        } else {
          statistics << "> weights were transferred for usage in Mayergoyz model" << std::endl;
        }
      }
      
//			
//			if ( (POL_operatorParams_.methodName_ == "vectorPreisach_Sutor")||(POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") ) {
//				statistics << "LEVENBERG-MARQUARDT: " << std::endl;
//				statistics << "- max number of iterations: " << LM_inversion_.maxNumIts << std::endl;
//				statistics << "- tolerance wrt x: " << LM_inversion_.tolH << std::endl;
//				statistics << "- tolerance wrt y: " << LM_inversion_.tolB << std::endl;
//				statistics << "- FD-resolution for Jacobian: " << LM_inversion_.jacRes << std::endl;
//				statistics << "LINESEARCH: " << std::endl;
//				statistics << "- alpha start: " << LM_inversion_.alphaLSStart << std::endl;
//				statistics << "- alpha min: " << LM_inversion_.alphaLSMin << std::endl;
//				statistics << "- alpha max: " << LM_inversion_.alphaLSMax << std::endl;		
//			}
		}
		
		/*
     * 2. Perform tests
     */
		forwardTimer = NULL;
		backwardTimer = NULL;
		if (measurePerformance) {
			forwardTimer = new Timer();
			backwardTimer = new Timer();
			performance << "+++ RUNTIMES +++" << std::endl;
			performance << "TEST: " << name << std::endl;
			performance << "MODEL: " << POL_operatorParams_.methodName_ << std::endl;
			if (POL_operatorParams_.methodName_ == "vectorPreisach_Sutor") {
				if (POL_operatorParams_.evalVersion_ == 1) {
					performance << "- classical model, list based implementation" << std::endl;
				} else if (POL_operatorParams_.evalVersion_ == 2) {
					performance << "- revised model, list based implementation" << std::endl;
				} else if (POL_operatorParams_.evalVersion_ == 10) {
					performance << "- classical model, matrix based implementation" << std::endl;
				} else if (POL_operatorParams_.evalVersion_ == 20) {
					performance << "- revised model, matrix based implementation" << std::endl;
				}
			}
			performance << "SINGLE TESTS: " << std::endl;

		}
    
		if(writeResultsToFile){
      //      std::stringstream results;
      //			results << "# +++ RESULTS +++" << std::endl;
      //			results << "# TEST: " << name << std::endl;
      //			results << "# MODEL: " << POL_operatorParams_.methodName_ << std::endl;
      //			if (POL_operatorParams_.methodName_ == "vectorPreisach_Sutor") {
      //				if (POL_operatorParams_.evalVersion_ == 1) {
      //					results << "# - classical model, list based implementation" << std::endl;
      //				} else if (POL_operatorParams_.evalVersion_ == 2) {
      //					results << "# - revised model, list based implementation" << std::endl;
      //				} else if (POL_operatorParams_.evalVersion_ == 10) {
      //					results << "# - classical model, matrix based implementation" << std::endl;
      //				} else if (POL_operatorParams_.evalVersion_ == 20) {
      //					results << "# - revised model, matrix based implementation" << std::endl;
      //				}
      //			}
      //			
      //      results_x << results.str();
      //      results_p << results.str();
      //      results_y << results.str();
      //      angularResults_x << results.str();
      //      angularResults_p << results.str();
      //      angularResults_y << results.str();
      
			if(testInversion){
        results_x << "#Number xIn[0] xIn[1] xOut[0] xOut[1]" << std::endl;
        results_x << "# > xIn[0],xIn[1] = x and y components of input; given by test signal" << std::endl;
        results_x << "# > xOut[0],xOut[1] = x and y components of inversion output" << std::endl;
        
        results_p << "#Number pIn[0] pIn[1] pOut[0] pOut[1]" << std::endl;
        results_p << "# > pIn[0],pIn[1] = x and y components of hyst operator computed from xIn" << std::endl;
        results_p << "# > pOut[0],pOut[1] = x and y components of hyst operator computed from xOut" << std::endl;
        
        results_y << "#Number yIn[0] yIn[1] yOut[0] yOut[1]" << std::endl;
				results_y << "# > yIn[0],yIn[1] = x and y components of output of hyst operator to xIn; used as input for inversion" << std::endl;
        results_y << "# > yOut[0],yOut[1] = x and y components computed from xOut" << std::endl;
        
        angularResults_x << "#Number abs(xIn) atan2(xIn[1]/xIn[0])*180/pi abs(xOut) atan2(xOut[1]/xOut[0])*180/pi" << std::endl;
        angularResults_x << "# > xIn[0],xIn[1] = x and y components of input; given by test signal" << std::endl;
        angularResults_x << "# > xOut[0],xOut[1] = x and y components of inversion output" << std::endl;
        
        angularResults_p << "#Number abs(pIn) atan2(pIn[1]/pIn[0])*180/pi abs(pOut) atan2(pOut[1]/pOut[0])*180/pi" << std::endl;
        angularResults_p << "# > pIn[0],pIn[1] = x and y components of hyst operator computed from xIn" << std::endl;
        angularResults_p << "# > pOut[0],pOut[1] = x and y components of hyst operator computed from xOut" << std::endl;
        
        angularResults_y << "#Number abs(yIn) atan2(yIn[1]/yIn[0])*180/pi abs(yOut) atan2(yOut[1]/yOut[0])*180/pi" << std::endl;
        angularResults_y << "# > yIn[0],yIn[1] = x and y components of output of hyst operator to xIn; used as input for inversion" << std::endl;
        angularResults_y << "# > yOut[0],yOut[1] = x and y components computed from xOut" << std::endl;
 
			} else {
        results_x << "#Number xIn[0] xIn[1]" << std::endl;
        results_x << "# > xIn[0],xIn[1] = x and y components of input; given by test signal" << std::endl;
        
        results_p << "#Number pOut[0] pOut[1]" << std::endl;
        results_p << "# > pOut[0],pOut[1] = x and y components of hyst operator computed from xOut" << std::endl;
        
        results_y << "#Number yOut[0] yOut[1]" << std::endl;
        results_y << "# > yOut[0],yOut[1] = x and y components computed from xOut" << std::endl;
        
        angularResults_x << "#Number abs(xIn) atan2(xIn[1]/xIn[0])*180/pi" << std::endl;
        angularResults_x << "# > xIn[0],xIn[1] = x and y components of input; given by test signal" << std::endl;
        
        angularResults_p << "#Number abs(pOut) atan2(pOut[1]/pOut[0])*180/pi" << std::endl;
        angularResults_p << "# > pOut[0],pOut[1] = x and y components of hyst operator computed from xOut" << std::endl;
        
        angularResults_y << "#Number abs(yOut) atan2(yOut[1]/yOut[0])*180/pi" << std::endl;
        angularResults_y << "# > yOut[0],yOut[1] = x and y components computed from xOut" << std::endl;  

			}
      
      orientationTowardsExcitation_p << "# Projetion of p_x and p_y to p_perpendicular and p_parallel" << std::endl;
      orientationTowardsExcitation_p << "# with p_parallel = p along excitation axis " << std::endl;
      orientationTowardsExcitation_p << "# and p_perpendicular = p perpendicualr to excitation axis " << std::endl;
      orientationTowardsExcitation_p << "# Columns: " << std::endl;
      orientationTowardsExcitation_p << "# 1. Testnumber " << std::endl;
      orientationTowardsExcitation_p << "# 2. Excitation amplitude " << std::endl;
      orientationTowardsExcitation_p << "# 3. Excitation angle towards x-axis " << std::endl;
      orientationTowardsExcitation_p << "# 4. p parallel to excitation " << std::endl;
      orientationTowardsExcitation_p << "# 5. p perpendicular to excitation " << std::endl;
		}
		
		
		bool overwriteMemory = false;
		// always overwrite Direction; deprecated flag; remove in future
    //		bool overwriteDirection = true;
		
		if( (printStatistics) ){
			std::cout << "##### STARTING TEST " << name << " #####" << std::endl;
		}
		hOutOld.Init();
    hOutOldForStrains.Init();
    
		for(UInt i = 0; i < totalSteps; i++){
			if( (i%10 == 0)&&(printStatistics) ){
				std::cout << "STEP NR " << i+1 << "/" << totalSteps << " #####" << std::endl;
			}
      performance << "Test " << i+1 << std::endl;
			
			// set previous state
			if(i > 0){
        xInPrev = xIn; 
				xPrev = xOut; // xOut = retrieved output from inversion; only set and used in case of testInversion
				hPrev = hOut; // hOut = output of hyst operator from forward computation with overwrite = true
				yPrev = yOut; // yOut = calculated y = h + eps*x from computation with overwrite = true
			} else {
        xInPrev.Init();
				xPrev.Init();
				hPrev.Init();
				yPrev.Init();
			}
			
			xIn.Init();
			xIn[0] = xVals[i];
			if(test1D == false){
				xIn[1] = yVals[i];
			}
			
      xInBak = xIn;
      
      if(testInversion){
				/*
         * Inversion test:
         *	1. evaluate forward hyst operator with given input; overwrite memory = false
         *  2. retrieve input by passing result of 2 to inverse hyst operator
         *			> measure backward
         *  3. pass retrieved input to forward hyst operator; overwrite memroy = true
         *			> measure forward time
         */
				
				// 1. forward; do not overwrite
				overwriteMemory = false;
				
				successFlagForward = -1;
				hIn.Init();
				if(vector){
          hIn = hystTMP->computeValue_vec(xIn, 0, overwriteMemory, debugOut, successFlagForward);
        } else {  
					hIn[0] = hystTMP->computeValueAndUpdate(xIn[0], 0, overwriteMemory, successFlagForward);
        }
				
				yIn.Init();
        yIn.Add(1.0,hIn);
        for(UInt j = 0; j < dim_; j++){
          yIn[j] += POL_eps_mu_SmallSignal_[j][j]*xIn[j];
        }
				
				// 2. backward; do not overwrite
				overwriteMemory = false;
        successFlagBackward = -1;
				numberOfLMIterations = 0;
				numberOfLinesearchIterations = 0;
				maxNumberOfLinesearchIterations = 0;
				
        if (measurePerformance) {
          backwardTimer->Start();
					startTime = backwardTimer->GetCPUTime();
        }
				
        if(vector){
          xOut = hystTMP->computeInput_vec_withStatistics(yIn, yPrev, xPrev, hPrev, 
                  0, POL_eps_mu_SmallSignal_, POL_operatorParams_.fieldsAlignedAboveSat_, POL_operatorParams_.hystOutputRestrictedToSat_,
                  numberOfLMIterations, numberOfLinesearchIterations, maxNumberOfLinesearchIterations,
                  successFlagBackward, minAlpha, maxAlpha, avgAlpha, xIn);	
        } else {	
					xOut[0] = hystTMP->computeInputAndUpdate(yIn[0], POL_eps_mu_SmallSignal_[0][0], 
                  0, overwriteMemory, successFlagBackward);
        }
				
        if (measurePerformance) {
          backwardTimer->Stop();
					
          // only count actual computations; if value gets reused, this should not count
          if(successFlagBackward != 0) {
            backwardEvalCounter++;
            
            endTime = backwardTimer->GetCPUTime();
            evalTime = endTime-startTime;
            if(evalTime > backwardMaxEvalTime){
              backwardMaxEvalTime = evalTime;
            }
            backwardTotalEvalTime += evalTime;
						
						performance << "- backward: " << evalTime << std::endl;
          } else {
            performance << "- backward: " << evalTime << "(reused)" << std::endl;
					}
        }
				
				if(successFlagBackward == -1){
          LMFails++;
        } else if(successFlagBackward == 0){
          totalReused++;
        } else if(successFlagBackward == 1){
          totalAnhystOnly++;
        } else if(successFlagBackward == 2){
          totalBisection++;
        } else if(successFlagBackward == 3){
          totalRemanence++;
        } else if(successFlagBackward == 4){
          totalPassedErrorTol++;
        } else if(successFlagBackward == 5){
          totalPassedResTolX++;
        } else if(successFlagBackward == 6){
          totalPassedResTolY++;
        }
        successCodeVectorBackward[i] = successFlagBackward;
				
				if(minAlpha < totalminAlpha){
          totalminAlpha = minAlpha;
        }
        if(maxAlpha > totalmaxAlpha){
          totalmaxAlpha = maxAlpha;
        }
        totalavgAlpha += avgAlpha;
        
				totalNumberOfLMIterations += numberOfLMIterations;
				totalNumberOfLinesearchIterations += numberOfLinesearchIterations;
				if(maxNumberOfLinesearchIterations > absmaxNumberOfLinesearchIterations){
					absmaxNumberOfLinesearchIterations = maxNumberOfLinesearchIterations;
				}
				
				/*
         * Check result
         */
				xErr.Init();
				xErr.Add(-1.0,xIn,1.0,xOut);
				// check below after eval of step 3
				
				// set xIn to xOut, so that during forward evaluation, xOut is used
				xIn = xOut;
        
			} // testInversion
			
			// 3. forward; overwrite this time and measure performance
      // NOTE: in case of inversion, we set the memory of the system with the RETRIEVED SOLUTION
      // i.e. if inversion fails, it will also affect the reference solution!
			overwriteMemory = true;
			
			successFlagForward = -1;
			hOut.Init();
      hOutForStrains.Init();
			
			if (measurePerformance) {
				forwardTimer->Start();
				startTime = forwardTimer->GetCPUTime();
			}
			
			if(vector){
				hOut = hystTMP->computeValue_vec(xIn, 0, overwriteMemory, debugOut, successFlagForward);
			} else {  
				hOut[0] = hystTMP->computeValueAndUpdate(xIn[0], 0, overwriteMemory, successFlagForward);
			}
			
      if(CouplingParams_.ownHystOperator_){
        if(vector){
          hOutForStrains = hystStrainTMP->computeValue_vec(xIn, 0, overwriteMemory, debugOut, successFlagForward);
        } else {  
          hOutForStrains[0] = hystStrainTMP->computeValueAndUpdate(xIn[0], 0, overwriteMemory, successFlagForward);
        }
      } else {
        hOutForStrains = hOut;
      }
      
			if (measurePerformance) {
				forwardTimer->Stop();
				
				// only count actual computations; if value gets reused, this should not count
				if(successFlagForward != 0) {
					forwardEvalCounter++;
					
					endTime = forwardTimer->GetCPUTime();
					evalTime = endTime-startTime;
					if(evalTime > forwardMaxEvalTime){
						forwardMaxEvalTime = evalTime;
					}
					forwardTotalEvalTime += evalTime;
					
					performance << "- forward: " << evalTime << std::endl;
				} else {
					performance << "- forward: " << evalTime << "(reused)" << std::endl;
				}
			}
			
			if(successFlagForward == -1){
				forwardFails++;
			} else if(successFlagForward == 0){
				forwardReused++;
			} else if(successFlagForward == 1){
				forwardAnhystOnly++;
			} else if(successFlagForward == 2){
				forwardOverwrite++;
			} else if(successFlagForward == 3){
				forwardTMP++; 
			}
			successCodeVectorForward[i] = successFlagForward;
			
			yOut.Init();
			yOut.Add(1.0,hOut);
			for(UInt j = 0; j < dim_; j++){
        // note that in case of inversion, xIn is the RETRIEVED input!
				yOut[j] += POL_eps_mu_SmallSignal_[j][j]*xIn[j];
			}
			
      if(outputIrrStrains){
        Vector<Double> S_irr = ComputeIrreversibleStrains(hOutForStrains,xIn,hOutOldForStrains);
        results_s << i+1 << " " << std::setprecision(9) << S_irr.ToString() << std::endl;
        results_xps << i+1 << " " << std::setprecision(9) << xIn.ToString() << " " << hOutForStrains.ToString() << " " << S_irr.ToString() << std::endl;
      }
      if(writeResultsToFile){
        results_xp << i+1 << " " << std::setprecision(9) << xIn.ToString() << " " << hOut.ToString() << std::endl;
        
        Double xAmpl = xIn.NormL2();
        Double xAngle = std::atan2(xIn[1],xIn[0])/M_PI*180;
        
        Vector<Double> xDir = Vector<Double>(dim_);
        xDir.Init();
        if(xAmpl != 0){
          xDir.Add(1.0/xAmpl,xIn);
        } else {
          Double xAmplOld = xInPrev.NormL2();
          if(xAmplOld != 0){
            xDir.Add(1.0/xAmplOld,xInPrev);
          } 
        }
        Double pParallel = hOut.Inner(xDir);
        Vector<Double> pPerpendicularVec = Vector<Double>(dim_);
        pPerpendicularVec.Add(1.0,hOut,-pParallel,xDir);
        Double pPerpendicular = pPerpendicularVec.NormL2();
        
        orientationTowardsExcitation_p << i+1 << " " << std::setprecision(9) << xAmpl << " " << xAngle << " " << pParallel << " " << pPerpendicular << std::endl;
      }
			/*
       * 4. Check result
       * > only for inversion test; for forward test, we do not have a
       *		value for comparison
       */
			bool failWRTX = false;
			bool failWRTY = false;
			if(testInversion){
				yErr.Init();
				yErr.Add(-1.0,yIn,1.0,yOut);
				
				if(xErr.NormL2() > LM_inversion_.tolH){
					failWRTX = true;
					failedTests_xIn[i] = xInBak;
					failedTests_xOut[i] = xOut;
					failedTests_yIn[i] = yIn;
					failedTests_yOut[i] = yOut;
				}
				
				if(yErr.NormL2() > LM_inversion_.tolB){
					failWRTY = true;
					failedTests_xIn[i] = xInBak;
					failedTests_xOut[i] = xOut;
					failedTests_yIn[i] = yIn;
					failedTests_yOut[i] = yOut;
				}	
				
				if( failWRTX || failWRTY ){
					failedTests_wrtX[i] = std::pair< bool, Double >(failWRTX, xErr.NormL2());
					failedTests_wrtY[i] = std::pair< bool, Double >(failWRTY, yErr.NormL2());
				}
        
        if( failWRTX && !failWRTY ){
          numHalfFails++;
        }
        
        if( failWRTX && failWRTY ){
          numFails++;
        }
			}
			
			/*
       * 5. Output and statistics
       */
      std::string delimiter = " ";
			if(writeResultsToFile){
        Double angle1 = 0;
        Double ampl1 = 0;
        Double angle2 = 0;
        Double ampl2 = 0;
				if(testInversion){
					//results << "# Number \t xIn[0] \t xIn[1] \t yIn[0] \t yIn[1] \t xOut[0] \t xOut[1] \t yOut[0] \t yOut[1]" << std::endl;
					results_x << i+1 << std::setprecision(9) << delimiter << xInBak[0] << delimiter << xInBak[1] << delimiter << xOut[0] << delimiter << xOut[1] << std::endl;
          ampl1 = xIn.NormL2();
          angle1 = atan2(xInBak[1],xInBak[0])*180.0/M_PI;
          ampl2 = xOut.NormL2();
          angle2 = atan2(xOut[1],xOut[0])*180.0/M_PI;
          angularResults_x << i+1 << std::setprecision(9) << delimiter << ampl1 << delimiter << angle1 << delimiter << ampl2 << delimiter << angle2 << std::endl;
          
          results_p << i+1 << std::setprecision(9) << delimiter << hIn[0] << delimiter << hIn[1] << delimiter << hOut[0] << delimiter << hOut[1] << std::endl;
          ampl1 = hIn.NormL2();
          angle1 = atan2(hIn[1],hIn[0])*180.0/M_PI;
          ampl2 = hOut.NormL2();
          angle2 = atan2(hOut[1],hOut[0])*180.0/M_PI;
          angularResults_p << i+1 << std::setprecision(9) << delimiter << ampl1 << delimiter << angle1 << delimiter << ampl2 << delimiter << angle2 << std::endl;
          
          results_y << i+1 << std::setprecision(9) << delimiter << yIn[0] << delimiter << yIn[1] << delimiter << yOut[0] << delimiter << yOut[1] << std::endl;
          ampl1 = yIn.NormL2();
          angle1 = atan2(yIn[1],yIn[0])*180.0/M_PI;
          ampl2 = yOut.NormL2();
          angle2 = atan2(yOut[1],yOut[0])*180.0/M_PI;
          angularResults_y << i+1 << std::setprecision(9) << delimiter << ampl1 << delimiter << angle1 << delimiter << ampl2 << delimiter << angle2 << std::endl;
				} else {
					//					results << "# Number \t xIn[0] \t xIn[1] \t yOut[0] \t yOut[1]" << std::endl;
          results_x << i+1 << std::setprecision(9) << delimiter << xInBak[0] << delimiter << xInBak[1] << std::endl;
          ampl1 = xIn.NormL2();
          angle1 = atan2(xInBak[1],xInBak[0])*180.0/M_PI;
          angularResults_x << i+1 << std::setprecision(9) << delimiter << ampl1 << delimiter << angle1 << std::endl;
          
          results_p << i+1 << std::setprecision(9) << delimiter << hOut[0] << delimiter << hOut[1] << std::endl;
          ampl1 = hOut.NormL2();
          angle1 = atan2(hOut[1],hOut[0])*180.0/M_PI;
          angularResults_p << i+1 << std::setprecision(9) << delimiter << ampl1 << delimiter << angle1 << std::endl;
          
          results_y << i+1 << std::setprecision(9) << delimiter << yOut[0] << delimiter << yOut[1] << std::endl; 
          ampl1 = yOut.NormL2();
          angle1 = atan2(yOut[1],yOut[0])*180.0/M_PI;
          angularResults_y << i+1 << std::setprecision(9) << delimiter << ampl1 << delimiter << angle1 << std::endl;
				}
			}
			
			LOG_TRACE(coeffcthyst) << "##### TARGET X-VECTOR #####";
			LOG_TRACE(coeffcthyst) << xInBak.ToString();
			if(testInversion){
				LOG_TRACE(coeffcthyst) << "##### RETRIEVED X-VECTOR #####";
				LOG_TRACE(coeffcthyst) << xOut.ToString();
				
				LOG_TRACE(coeffcthyst) << "##### ERROR VECTOR wrt X #####";
				LOG_TRACE(coeffcthyst) << xErr.ToString();
				LOG_TRACE(coeffcthyst) << "##### TARGET Y-VECTOR #####";
				LOG_TRACE(coeffcthyst) << yIn.ToString();
			}
			LOG_TRACE(coeffcthyst) << "##### RETRIEVED Y-VECTOR #####";
			LOG_TRACE(coeffcthyst) << yOut.ToString();
			if(testInversion){
				LOG_TRACE(coeffcthyst) << "##### ERROR VECTOR wrt Y #####";
				LOG_TRACE(coeffcthyst) << yErr.ToString();
			}
			hOutOld = hOut;
      hOutOldForStrains = hOutForStrains;
		}
		
		UInt LMcases = totalSteps-totalReused-totalBisection-totalAnhystOnly-totalRemanence;
		if(LMcases != 0){
			totalavgAlpha /= (Double) LMcases;
		}
		
		if(printStatistics){
			std::stringstream majorResults;
			majorResults << "############################# " <<	std::endl;	
			majorResults << "##### RESULTS FOR TEST " << name <<	std::endl;	
			majorResults << " " << totalSteps-numFails << " of " << totalSteps << " satisfied at least the failback criterion, i.e. |residual_Y| = |Y - mu_eps*X - P(X)| < " << LM_inversion_.tolB << std::endl;
			majorResults << " " << numHalfFails << " of " << totalSteps << " failed error criterion but passed due to failback criterion" << std::endl;
      majorResults << " " << numFails << " of " << totalSteps << " failed to satisfy even the failback criterion" << std::endl;
			
			// print major results to stdout
			std::cout << majorResults.str() << std::endl;
			
			statistics << majorResults.str() << std::endl;
			if( (numFails > 0)||(numHalfFails > 0) ){
				statistics << "## Detailed Statistics for failed tests: " << std::endl;
				
				std::map< UInt, std::pair< bool, Double > >::iterator mapItX;
				std::map< UInt, std::pair< bool, Double > >::iterator mapItY = failedTests_wrtY.begin();
				bool failX, failY;
				for(mapItX = failedTests_wrtX.begin(); mapItX != failedTests_wrtX.end(); mapItX++,mapItY++){
//          std::cout << "1" << std::endl;
					statistics << " Test Nr " << mapItX->first << std::endl;
					failX = mapItX->second.first;
					failY = mapItY->second.first;
//          std::cout << "2" << std::endl;
					bool xsolabovesat = failedTests_xIn[mapItX->first].NormL2()>POL_operatorParams_.inputSat_;
					bool xretabovesat = failedTests_xOut[mapItX->first].NormL2()>POL_operatorParams_.inputSat_;
					statistics << " > X_in: " << failedTests_xIn[mapItX->first].ToString()  << "; above sat? "<< xsolabovesat << std::endl;
					statistics << " > X_out: " << failedTests_xOut[mapItX->first].ToString() << "; above sat? "<< xretabovesat << std::endl;
					if(failX){
						statistics << " > FAIL - remaining error wrt X  " << mapItX->second.second << std::endl;
					} else {
						statistics << " > SUCCESS - remaining error wrt X  " << mapItX->second.second << std::endl;
					}
//					std::cout << "3" << std::endl;
					statistics << " > Y_in: " << failedTests_yIn[mapItX->first].ToString() << std::endl;
					statistics << " > Y_out: " << failedTests_yOut[mapItX->first].ToString() << std::endl;
					if(failY){
						statistics << " > FAIL - remaining error wrt Y  " << mapItY->second.second << std::endl;
					} else {
						statistics << " > SUCCESS - remaining error wrt Y  " << mapItY->second.second << std::endl;
					}
//					std::cout << "4" << std::endl;
					statistics << " - SuccessCode (forward): " << successCodeVectorForward[mapItX->first] << std::endl;
					statistics << " - SuccessCode (backward): " << successCodeVectorBackward[mapItX->first] << std::endl;
				}
			}
			
			statistics << "## Detailed Statistics on forward evaluations (excluding those inside of LM/during inversion): " << std::endl;
			statistics << " " << forwardReused << " of " << totalSteps << " reused old output as input change was below tolerance" << std::endl;
			statistics << " " << forwardAnhystOnly << " of " << totalSteps << " featured only anhysteretic parts" << std::endl;
			statistics << " " << forwardOverwrite << " of " << totalSteps << " evaluations were performed on permanent storage" << std::endl;
			statistics << " " << forwardTMP << " of " << totalSteps << " evaluations were performed on temporal storage" << std::endl;
			
			statistics << "## Detailed Statistics on inversion: " << std::endl;
			statistics << " " << totalReused << " of " << totalSteps << " reused old solution as DeltaY < " << LM_inversion_.tolB << std::endl;
			statistics << " " << totalAnhystOnly << " of " << totalSteps << " were solved using bisection for pure anhysteretic case" << std::endl;
			statistics << " " << totalBisection << " of " << totalSteps << " were solved using bisection / simple division" << std::endl;
			if(!vector){
				statistics << " " << totalPassedErrorTol << " of " << totalSteps << " tests passed due to error tolerance" << std::endl;
			}
			
			if(vector){
				statistics << " " << totalRemanence << " of " << totalSteps << " cases were in remanence" << std::endl;
				statistics << " " << LMcases << " of " << totalSteps << " were solved using Levenberg-Marquardt" << std::endl;
				if(LMcases != 0){
					statistics << "## Detailed Statistics for Levenberg-Marquardt: " << std::endl;
					statistics << " " << totalPassedErrorTol << " of " << LMcases << " tests passed due to error tolerance |JacT*Res| < " << LM_inversion_.tolH << std::endl;
					statistics << " " << totalPassedResTolX << " of " << LMcases << " tests passed due to |residual_X| = |X - mu_eps^-1*(Y - P(X)| < " << LM_inversion_.tolH << std::endl;
					statistics << " " << totalPassedResTolY << " of " << LMcases << " tests passed due to |residual_Y| = |Y - mu_eps*X - P(X)| < " << LM_inversion_.tolH << std::endl;
					statistics << " " << LMFails << " of " << LMcases << " failed the failback criterion (|residual_Y| < " << LM_inversion_.tolB << ")" << std::endl;
					statistics << " Total number of LM iterations: " << totalNumberOfLMIterations << std::endl;
					statistics << " Average number of LM iterations: " << (Double) totalNumberOfLMIterations/LMcases << std::endl;
					statistics << " Total number of Linesearch iterations: " << totalNumberOfLinesearchIterations << std::endl;
					statistics << " Average number of Linesearch iterations (per LM Iteration): " << (Double) totalNumberOfLinesearchIterations/totalNumberOfLMIterations << std::endl;
					statistics << " Maximal number of Linesearch iterations: " << absmaxNumberOfLinesearchIterations << std::endl;
					statistics << " Minimal alphaLS: " << totalminAlpha << std::endl;
					statistics << " Maximal alphaLS: " << totalmaxAlpha << std::endl;
					statistics << " Average alphaLS: " << totalavgAlpha << std::endl;
				}
			}
			statistics << "############################# " <<	std::endl;	
      statistics << std::endl;	
		}
		
		
		if (measurePerformance){
			Double forwardAvgEvalTime = forwardTotalEvalTime/forwardEvalCounter;
			performance << "## Runtime information: " << std::endl;
      performance << "# Runtime information for forward evaluations: " << std::endl;
			performance << "Total number of forward evaluations: " << forwardEvalCounter << std::endl;
			performance << "Total time for forward evaluations: " << forwardTotalEvalTime << std::endl;
			performance << "Average time for forward evaluations: " << forwardAvgEvalTime << std::endl;
			performance << "Maximal time for forward evaluations: " << forwardMaxEvalTime << std::endl;
      performance << "> Note: only evaluations on permanent storage counted, " << std::endl;
      performance << "   i.e. no reused values and no evaluations on temporal storage (as e.g. inside LM)" << std::endl;
			if(testInversion){	
				Double backwardAvgEvalTime = backwardTotalEvalTime/backwardEvalCounter;
        performance << "# Runtime information for backward evaluations: " << std::endl;
				performance << "Total number of backward evaluations: " << backwardEvalCounter << std::endl;
				performance << "Total time for backward evaluations: " << backwardTotalEvalTime << std::endl;
				performance << "Average time for backward evaluations: " << backwardAvgEvalTime << std::endl;
				performance << "Maximal time for backward evaluations: " << backwardMaxEvalTime << std::endl;
        performance << "# Results for grepping: " << std::endl;
        //        performance << "Flag for grep \t forwardAvgEvalTime (#evals) \t backwardAvgEvalTime (#evals) \t  numFails "<< std::endl;
        //				performance << "TOGREP: \t\t" << forwardAvgEvalTime << " (" << forwardEvalCounter << ") " << "\t " << backwardAvgEvalTime << " (" << backwardEvalCounter << ") " << "\t " << numFails << std::endl;
        /*
         * NOTE: we do measure forward time when computing inversion, too, but the results will be much different
         *        from the results without inversion; inversion seems to affect forward evaluation time, even though no storage will be overwritten
         *        by LM (or should not) and when inspecting the return flag of forward evaluation it is not 0, i.e. result is not reused;
         *        nevertheless the forward runtime will be as short as if results had been reused
         *      > not sure why this is the case but as a consequence, we measure forward time and backward time during separate runs
         */
        performance << "Flag for grep \t backwardAvgEvalTime (#evals) \t  numFails "<< std::endl;
				performance << "TOGREP: \t\t" << "\t " << backwardAvgEvalTime << " (" << backwardEvalCounter << ") " << "\t " << numFails << std::endl;
			} else {
        performance << "# Results for grepping: " << std::endl;
        performance << "Flag for grep \t forwardAvgEvalTime (forwardEvalCounter) " << std::endl;
        performance << "TOGREP: \t\t" << forwardAvgEvalTime << " (" << forwardEvalCounter << ") " << std::endl;
      }
			//statistics << "TOGREP: " << forwardAvgEvalTime << " (" << forwardEvalCounter << ") " << forwardTotalEvalTime << " " << backwardAvgEvalTime << " (" << backwardEvalCounter << ") " << " " << backwardTotalEvalTime << " " << numFails << std::endl;			performance << "############################# " <<	std::endl;	
      performance << std::endl;	
      
      
      if(commonPerformanceFile != "---"){
        if(testInversion){
          Double backwardAvgEvalTime = backwardTotalEvalTime/backwardEvalCounter;
          commonPerfStream << name << "\t\t Forward: " << forwardAvgEvalTime << " (Evals=" << forwardEvalCounter << ") " << "\t Backward: " << backwardAvgEvalTime << " (Evals=" << backwardEvalCounter << ") " << "\t InversionFails=" << numFails << std::endl;
        } else {
          commonPerfStream << name << "\t\t Forward: " << forwardAvgEvalTime << " (Evals=" << forwardEvalCounter << ") " << "\t Backward: 0.0 (Evals=0) \t InversionFails=0" << std::endl;
        }
      }
		}
    
		if(writeResultsToFile){
			results_x.close();
      results_xp.close();
      results_p.close();
      results_y.close();
      angularResults_x.close();
      angularResults_p.close();
      angularResults_y.close();
      orientationTowardsExcitation_p.close();
		}
    if(outputIrrStrains){
      results_s.close();
    }
    
		if(printStatistics){
			statistics.close();
		} 
		if(measurePerformance){
			performance.close();
      if(commonPerformanceFile != "---"){
        commonPerfStream.close();
      }
		}
		
  }
  
  void CoefFunctionHyst::TestInversion(PtrParamNode testNode){
    
    if(testNode == NULL){
      EXCEPTION("No input to inversion test provided");
    }
		if(storageInitialized_ == false){
      // memory has not been set yet > no tests possible
			EXCEPTION("Memory has not been initialized yet");
    } 
    
    bool testInversion = false;
    if(testNode->Has("TestInversion")){
      testNode->GetValue("TestInversion",testInversion,ParamNode::PASS);
    }
    
    if( (testInversion == true)&&(inversionSet_ == false) ){
      WARN("Test of inversion demanded but no inversion parameter for current material. Skip inversion test.");
      testInversion = false;      
    }
    
    bool stopAfterTests = false;
    if(testNode->Has("StopAfterTests")){
      testNode->GetValue("StopAfterTests",stopAfterTests,ParamNode::PASS);
    }
    bool printStatistics = false;
    if(testNode->Has("PrintStatistics")){
      testNode->GetValue("PrintStatistics",printStatistics,ParamNode::PASS);
    }
    bool writeResultsToFile = false;
    if(testNode->Has("WriteResultsToFile")){
      testNode->GetValue("WriteResultsToFile",writeResultsToFile,ParamNode::PASS);
    }
    bool writeInputToFile = false;
    if(testNode->Has("WriteInputToFile")){
      testNode->GetValue("WriteInputToFile",writeInputToFile,ParamNode::PASS);
    }
    bool measurePerformance = false;
    std::string commonPerfFile = "---";
    if(testNode->Has("MeasurePerformance")){
      PtrParamNode PerformanceNode = testNode->Get("MeasurePerformance");
      if(PerformanceNode->Has("activate")){
        PerformanceNode->GetValue("activate",measurePerformance,ParamNode::PASS);
      }
      if(PerformanceNode->Has("commonResultFile")){
        PerformanceNode->GetValue("commonResultFile",commonPerfFile,ParamNode::PASS);
      }
    }
    bool outputIrrStrains = false;
    if(testNode->Has("OutputIrrStrains")){
      testNode->GetValue("OutputIrrStrains",outputIrrStrains,ParamNode::PASS);
    }
    
    //    std::cout << "testInversion? " << testInversion << std::endl;
    //    std::cout << "test1D? " << test1D << std::endl;
    //    std::cout << "StopAfterTests? " << stopAfterTests << std::endl;
    //    std::cout << "printStatistics? " << printStatistics << std::endl;
    //    std::cout << "writeInputToFile? " << writeInputToFile << std::endl;
    //    std::cout << "writeResultsToFile? " << writeResultsToFile << std::endl;
    //    std::cout << "measurePerformance? " << measurePerformance << std::endl;
    
    PtrParamNode InputSignals = testNode->Get("InputSignals");
    
    Vector<Double> xVals = Vector<Double>(1);
    Vector<Double> yVals = Vector<Double>(1);
    
    Double amplitudeScaling = 1.0;
    Double numPeriods = 5.0;
    UInt stepsPerPeriod = 20;
    UInt numberOfSteps = 100;
    if(InputSignals->Has("Sine")){
      if(InputSignals->Get("Sine")->Has("AmplitudeScaling")){
        InputSignals->Get("Sine")->GetValue("AmplitudeScaling",amplitudeScaling,ParamNode::PASS);
      }
      if(InputSignals->Get("Sine")->Has("NumPeriods")){
        InputSignals->Get("Sine")->GetValue("NumPeriods",numPeriods,ParamNode::PASS);
      }
      if(InputSignals->Get("Sine")->Has("StepsPerPeriod")){
        InputSignals->Get("Sine")->GetValue("StepsPerPeriod",stepsPerPeriod,ParamNode::PASS);
      }
      bool test1D = false;
      if(InputSignals->Get("Sine")->Has("Test1D")){
        InputSignals->Get("Sine")->GetValue("Test1D",test1D,ParamNode::PASS);
      }
      std::string outputName = "---";
      if(InputSignals->Get("Sine")->Has("OutputName")){
        InputSignals->Get("Sine")->GetValue("OutputName",outputName,ParamNode::PASS);
      }
      if(outputName == "---"){
        outputName = "Sine";
      }
      
      CreatePeriodicTestSignal("Sine",amplitudeScaling,numPeriods,stepsPerPeriod,xVals,yVals);
      TestHystOperatorWithSignal(outputName,xVals,yVals,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains);
      if(writeInputToFile){
        WriteSignalToFile("Sine_input",xVals,yVals);
      }
    }
    if(InputSignals->Has("Rotation")){
      if(InputSignals->Get("Rotation")->Has("AmplitudeScaling")){
        InputSignals->Get("Rotation")->GetValue("AmplitudeScaling",amplitudeScaling,ParamNode::PASS);
      }
      if(InputSignals->Get("Rotation")->Has("NumPeriods")){
        InputSignals->Get("Rotation")->GetValue("NumPeriods",numPeriods,ParamNode::PASS);
      }
      if(InputSignals->Get("Rotation")->Has("StepsPerPeriod")){
        InputSignals->Get("Rotation")->GetValue("StepsPerPeriod",stepsPerPeriod,ParamNode::PASS);
      }
      bool test1D = false;
      if(InputSignals->Get("Rotation")->Has("Test1D")){
        InputSignals->Get("Rotation")->GetValue("Test1D",test1D,ParamNode::PASS);
      }
      std::string outputName = "---";
      if(InputSignals->Get("Rotation")->Has("OutputName")){
        InputSignals->Get("Rotation")->GetValue("OutputName",outputName,ParamNode::PASS);
      }
      if(outputName == "---"){
        outputName = "Rotation";
      }
      
      CreatePeriodicTestSignal("Rotation",amplitudeScaling,numPeriods,stepsPerPeriod,xVals,yVals);
      TestHystOperatorWithSignal(outputName,xVals,yVals,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains);
      if(writeInputToFile){
        WriteSignalToFile("Rotation_input",xVals,yVals);
      }
    }
    if(InputSignals->Has("DecreasingRotation")){
      if(InputSignals->Get("DecreasingRotation")->Has("AmplitudeScaling")){
        InputSignals->Get("DecreasingRotation")->GetValue("AmplitudeScaling",amplitudeScaling,ParamNode::PASS);
      }
      if(InputSignals->Get("DecreasingRotation")->Has("NumPeriods")){
        InputSignals->Get("DecreasingRotation")->GetValue("NumPeriods",numPeriods,ParamNode::PASS);
      }
      if(InputSignals->Get("DecreasingRotation")->Has("StepsPerPeriod")){
        InputSignals->Get("DecreasingRotation")->GetValue("StepsPerPeriod",stepsPerPeriod,ParamNode::PASS);
      }
      bool test1D = false;
      if(InputSignals->Get("DecreasingRotation")->Has("Test1D")){
        InputSignals->Get("DecreasingRotation")->GetValue("Test1D",test1D,ParamNode::PASS);
      }
      std::string outputName = "---";
      if(InputSignals->Get("DecreasingRotation")->Has("OutputName")){
        InputSignals->Get("DecreasingRotation")->GetValue("OutputName",outputName,ParamNode::PASS);
      }
      if(outputName == "---"){
        outputName = "DecreasingRotation";
      }
            
      CreatePeriodicTestSignal("DecreasingRotation",amplitudeScaling,numPeriods,stepsPerPeriod,xVals,yVals);
      TestHystOperatorWithSignal(outputName,xVals,yVals,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains);
      if(writeInputToFile){
        WriteSignalToFile("DecreasingRotation_input",xVals,yVals);
      }
    }
    if(InputSignals->Has("IncreasingRotation")){
      if(InputSignals->Get("IncreasingRotation")->Has("AmplitudeScaling")){
        InputSignals->Get("IncreasingRotation")->GetValue("AmplitudeScaling",amplitudeScaling,ParamNode::PASS);
      }
      if(InputSignals->Get("IncreasingRotation")->Has("NumPeriods")){
        InputSignals->Get("IncreasingRotation")->GetValue("NumPeriods",numPeriods,ParamNode::PASS);
      }
      if(InputSignals->Get("IncreasingRotation")->Has("StepsPerPeriod")){
        InputSignals->Get("IncreasingRotation")->GetValue("StepsPerPeriod",stepsPerPeriod,ParamNode::PASS);
      }
      bool test1D = false;
      if(InputSignals->Get("IncreasingRotation")->Has("Test1D")){
        InputSignals->Get("IncreasingRotation")->GetValue("Test1D",test1D,ParamNode::PASS);
      }
      std::string outputName = "---";
      if(InputSignals->Get("IncreasingRotation")->Has("OutputName")){
        InputSignals->Get("IncreasingRotation")->GetValue("OutputName",outputName,ParamNode::PASS);
      }
      if(outputName == "---"){
        outputName = "IncreasingRotation";
      }
          
      CreatePeriodicTestSignal("IncreasingRotation",amplitudeScaling,numPeriods,stepsPerPeriod,xVals,yVals);
      TestHystOperatorWithSignal(outputName,xVals,yVals,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains);
      if(writeInputToFile){
        WriteSignalToFile("IncreasingRotation_input",xVals,yVals);
      }
    }
    
    if(InputSignals->Has("DecreasingSine")){
      if(InputSignals->Get("DecreasingSine")->Has("AmplitudeScaling")){
        InputSignals->Get("DecreasingSine")->GetValue("AmplitudeScaling",amplitudeScaling,ParamNode::PASS);
      }
      if(InputSignals->Get("DecreasingSine")->Has("NumPeriods")){
        InputSignals->Get("DecreasingSine")->GetValue("NumPeriods",numPeriods,ParamNode::PASS);
      }
      if(InputSignals->Get("DecreasingSine")->Has("StepsPerPeriod")){
        InputSignals->Get("DecreasingSine")->GetValue("StepsPerPeriod",stepsPerPeriod,ParamNode::PASS);
      }
      bool test1D = false;
      if(InputSignals->Get("DecreasingSine")->Has("Test1D")){
        InputSignals->Get("DecreasingSine")->GetValue("Test1D",test1D,ParamNode::PASS);
      }
      std::string outputName = "---";
      if(InputSignals->Get("DecreasingSine")->Has("OutputName")){
        InputSignals->Get("DecreasingSine")->GetValue("OutputName",outputName,ParamNode::PASS);
      }
      if(outputName == "---"){
        outputName = "DecreasingSine";
      }
      
      CreatePeriodicTestSignal("DecreasingSine",amplitudeScaling,numPeriods,stepsPerPeriod,xVals,yVals);
      TestHystOperatorWithSignal(outputName,xVals,yVals,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains);
      if(writeInputToFile){
        WriteSignalToFile("DecreasingSine_input",xVals,yVals);
      }
    }
    if(InputSignals->Has("DecreasingSawtooth")){
      if(InputSignals->Get("DecreasingSawtooth")->Has("AmplitudeScaling")){
        InputSignals->Get("DecreasingSawtooth")->GetValue("AmplitudeScaling",amplitudeScaling,ParamNode::PASS);
      }
      if(InputSignals->Get("DecreasingSawtooth")->Has("NumPeriods")){
        InputSignals->Get("DecreasingSawtooth")->GetValue("NumPeriods",numPeriods,ParamNode::PASS);
      }
      if(InputSignals->Get("DecreasingSawtooth")->Has("StepsPerPeriod")){
        InputSignals->Get("DecreasingSawtooth")->GetValue("StepsPerPeriod",stepsPerPeriod,ParamNode::PASS);
      }
      bool test1D = false;
      if(InputSignals->Get("DecreasingSawtooth")->Has("Test1D")){
        InputSignals->Get("DecreasingSawtooth")->GetValue("Test1D",test1D,ParamNode::PASS);
      }
      std::string outputName = "---";
      if(InputSignals->Get("DecreasingSawtooth")->Has("OutputName")){
        InputSignals->Get("DecreasingSawtooth")->GetValue("OutputName",outputName,ParamNode::PASS);
      }
      if(outputName == "---"){
        outputName = "DecreasingSawtooth";
      }
      
      CreatePeriodicTestSignal("DecreasingSawtooth",amplitudeScaling,numPeriods,stepsPerPeriod,xVals,yVals);
      TestHystOperatorWithSignal(outputName,xVals,yVals,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains);
      if(writeInputToFile){
        WriteSignalToFile("DecreasingSawtooth_input",xVals,yVals);
      }
    }
    if(InputSignals->Has("Forc")){
      if(InputSignals->Get("Forc")->Has("AmplitudeScaling")){
        InputSignals->Get("Forc")->GetValue("AmplitudeScaling",amplitudeScaling,ParamNode::PASS);
      }
      if(InputSignals->Get("Forc")->Has("NumPeriods")){
        InputSignals->Get("Forc")->GetValue("NumPeriods",numPeriods,ParamNode::PASS);
      }
      if(InputSignals->Get("Forc")->Has("StepsPerPeriod")){
        InputSignals->Get("Forc")->GetValue("StepsPerPeriod",stepsPerPeriod,ParamNode::PASS);
      }
      bool test1D = false;
      if(InputSignals->Get("Forc")->Has("Test1D")){
        InputSignals->Get("Forc")->GetValue("Test1D",test1D,ParamNode::PASS);
      }
      std::string outputName = "---";
      if(InputSignals->Get("Forc")->Has("OutputName")){
        InputSignals->Get("Forc")->GetValue("OutputName",outputName,ParamNode::PASS);
      }
      if(outputName == "---"){
        outputName = "Forc";
      }
      
      CreatePeriodicTestSignal("Forc",amplitudeScaling,numPeriods,stepsPerPeriod,xVals,yVals);
      TestHystOperatorWithSignal(outputName,xVals,yVals,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains);
      if(writeInputToFile){
        WriteSignalToFile("Forc_input",xVals,yVals);
      }
    }
    if(InputSignals->Has("SelfDesigned")){
      if(InputSignals->Get("SelfDesigned")->Has("AmplitudeScaling")){
        InputSignals->Get("SelfDesigned")->GetValue("AmplitudeScaling",amplitudeScaling,ParamNode::PASS);
      }
      if(InputSignals->Get("SelfDesigned")->Has("NumSteps")){
        InputSignals->Get("SelfDesigned")->GetValue("NumSteps",numberOfSteps,ParamNode::PASS);
      }
      bool test1D = false;
      if(InputSignals->Get("SelfDesigned")->Has("Test1D")){
        InputSignals->Get("SelfDesigned")->GetValue("Test1D",test1D,ParamNode::PASS);
      }
      std::string outputName = "---";
      if(InputSignals->Get("SelfDesigned")->Has("OutputName")){
        InputSignals->Get("SelfDesigned")->GetValue("OutputName",outputName,ParamNode::PASS);
      }
      if(outputName == "---"){
        outputName = "SelfDesigned";
      }
      
      CreateNonPeriodicTestSignal("SelfDesigned",amplitudeScaling,numberOfSteps,xVals,yVals);
      TestHystOperatorWithSignal(outputName,xVals,yVals,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains);
      if(writeInputToFile){
        WriteSignalToFile("SelfDesigned_input",xVals,yVals);
      }
    }
    if(InputSignals->Has("SatX-RemX-SatY")){
      if(InputSignals->Get("SatX-RemX-SatY")->Has("AmplitudeScaling")){
        InputSignals->Get("SatX-RemX-SatY")->GetValue("AmplitudeScaling",amplitudeScaling,ParamNode::PASS);
      }
      if(InputSignals->Get("SatX-RemX-SatY")->Has("NumSteps")){
        InputSignals->Get("SatX-RemX-SatY")->GetValue("NumSteps",numberOfSteps,ParamNode::PASS);
      }
      bool test1D = false;
      if(InputSignals->Get("SatX-RemX-SatY")->Has("Test1D")){
        InputSignals->Get("SatX-RemX-SatY")->GetValue("Test1D",test1D,ParamNode::PASS);
      }
      std::string outputName = "---";
      if(InputSignals->Get("SatX-RemX-SatY")->Has("OutputName")){
        InputSignals->Get("SatX-RemX-SatY")->GetValue("OutputName",outputName,ParamNode::PASS);
      }
      if(outputName == "---"){
        outputName = "SatX-RemX-SatY";
      }
      
      CreateNonPeriodicTestSignal("SatX-RemX-SatY",amplitudeScaling,numberOfSteps,xVals,yVals);
      TestHystOperatorWithSignal(outputName,xVals,yVals,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains);
      if(writeInputToFile){
        WriteSignalToFile("SatX-RemX-SatY_input",xVals,yVals);
      }
    }
    
    if(InputSignals->Has("RemDrop")){
      if(InputSignals->Get("RemDrop")->Has("AmplitudeScaling")){
        InputSignals->Get("RemDrop")->GetValue("AmplitudeScaling",amplitudeScaling,ParamNode::PASS);
      }
      if(InputSignals->Get("RemDrop")->Has("NumSteps")){
        InputSignals->Get("RemDrop")->GetValue("NumSteps",numberOfSteps,ParamNode::PASS);
      }
      bool test1D = false;
      if(InputSignals->Get("RemDrop")->Has("Test1D")){
        InputSignals->Get("RemDrop")->GetValue("Test1D",test1D,ParamNode::PASS);
      }
      std::string outputName = "---";
      if(InputSignals->Get("RemDrop")->Has("OutputName")){
        InputSignals->Get("RemDrop")->GetValue("OutputName",outputName,ParamNode::PASS);
      }
      if(outputName == "---"){
        outputName = "RemDrop";
      }
      
      CreateNonPeriodicTestSignal("RemDrop",amplitudeScaling,numberOfSteps,xVals,yVals);
      TestHystOperatorWithSignal(outputName,xVals,yVals,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains);
      if(writeInputToFile){
        WriteSignalToFile("RemDrop_input",xVals,yVals);
      }
    }
    
    
    if(InputSignals->Has("ReadFromFile")){
      if(InputSignals->Get("ReadFromFile")->Has("AmplitudeScaling")){
        InputSignals->Get("ReadFromFile")->GetValue("AmplitudeScaling",amplitudeScaling,ParamNode::PASS);
      }
      std::string fileNameX = "None";
      if(InputSignals->Get("ReadFromFile")->Has("FileContaining_x_over_steps")){
        InputSignals->Get("ReadFromFile")->GetValue("FileContaining_x_over_steps",fileNameX,ParamNode::PASS);
      }
      std::string fileNameY = "None";
      if(InputSignals->Get("ReadFromFile")->Has("FileContaining_y_over_steps")){
        InputSignals->Get("ReadFromFile")->GetValue("FileContaining_y_over_steps",fileNameY,ParamNode::PASS);
      }
      bool test1D = false;
      if(InputSignals->Get("ReadFromFile")->Has("Test1D")){
        InputSignals->Get("ReadFromFile")->GetValue("Test1D",test1D,ParamNode::PASS);
      }
      
      Vector<Double> steps = Vector<Double>(1);
      Vector<Double> xVals = Vector<Double>(1);
      Vector<Double> yVals = Vector<Double>(1);
      std::stringstream combinedname;
      combinedname << "Test-";
      Interpolate1D reader = Interpolate1D();
      UInt maxSize = 1;
      if(fileNameX != "None"){
        reader.ReadFile(fileNameX.c_str(),steps,xVals);
        combinedname << fileNameX;
        if(xVals.GetSize() > maxSize){
          maxSize = xVals.GetSize();
        }
      }
      
      combinedname << "-";
      if(fileNameY != "None"){
        reader.ReadFile(fileNameY.c_str(),steps,yVals);
        combinedname << fileNameY;
        if(yVals.GetSize() > maxSize){
          maxSize = yVals.GetSize();
        }
      } 
      
      std::string outputName = "---";
      if(InputSignals->Get("ReadFromFile")->Has("OutputName")){
        InputSignals->Get("ReadFromFile")->GetValue("OutputName",outputName,ParamNode::PASS);
      }
      if(outputName == "---"){
        outputName = combinedname.str();
      }
      
      if(xVals.GetSize() < maxSize){
        xVals.Resize(maxSize);
      }
      if(yVals.GetSize() < maxSize){
        yVals.Resize(maxSize);
      }
      
      xVals.ScalarMult(amplitudeScaling);
      yVals.ScalarMult(amplitudeScaling);
      
      TestHystOperatorWithSignal(outputName,xVals,yVals,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains);
      
      combinedname << "_input";
      if(writeInputToFile){
        WriteSignalToFile(combinedname.str(),xVals,yVals);
      }
    }
    
    if(stopAfterTests){
      EXCEPTION("Stop after testing");
    }
    
    
    //     bool stopAfterInversionTests = false;
    //      bool printStatistics = false;
    //      bool writeResultsToFiles = false;
    //      std::string testFile;
    //      if( nonLinNode->Has("HYST_testInversion") ) {
    //        nonLinNode->Get( "HYST_testInversion")->GetValue( "Testnumber", testInvesion_, ParamNode::PASS );
    //        nonLinNode->Get( "HYST_testInversion")->GetValue( "StopAfterTests", stopAfterInversionTests, ParamNode::PASS );
    //        nonLinNode->Get( "HYST_testInversion")->GetValue( "PrintStatistics", printStatistics, ParamNode::PASS );
    //        nonLinNode->Get( "HYST_testInversion")->GetValue( "WriteResultsToFiles", writeResultsToFiles, ParamNode::PASS );
    //        nonLinNode->Get( "HYST_testInversion")->GetValue( "FileToRead", testFile, ParamNode::PASS );
    //        std::cout << "testFile: " << testFile << std::endl;
    //      } else {
    //        testInvesion_ = 0;
    //      }
    //      
    //      if(testFile != "none"){
    //        Vector<Double> xVals = Vector<Double>(1);
    //        Vector<Double> yVals = Vector<Double>(1);
    //        Interpolate1D reader = Interpolate1D();
    //        reader.ReadFile(testFile.c_str(),xVals,yVals);
    //        std::cout << "Read in data: " << std::endl;
    //        std::cout << xVals.ToString() << std::endl;
    //        std::cout << yVals.ToString() << std::endl;
    //      }
    
    //    
    //    OLD
    //    
    //  UInt testNumber = 0;
    //  bool abortAfterTests = false;
    ////  bool printStatistics = false;
    //  bool writeResultsToFiles = false;
    ////  bool testForwardOnly = false;
    //  //void CoefFunctionHyst::TestInversion(UInt testNumber, bool abortAfterTests, bool printStatistics, bool writeResultsToFiles, bool testForwardOnly){
    //    
    //    //std::cout << "POL_weightParams_.weightTensor_: " << POL_weightParams_.weightTensor_.ToString() << std::endl;
    //    
    //
    //    
    //    Matrix<Double> eps_mu = POL_eps_mu_SmallSignal_;
    //    
    //    std::stringstream generalInfo;
    //    generalInfo << std::endl;
    //    generalInfo << "############################" <<std::endl;
    //    generalInfo << "#### - INVERSION TEST - ####" <<std::endl;
    //    generalInfo << "############################" <<std::endl;
    //    generalInfo << std::endl;
    //    generalInfo << "## Used material parameter for all following tests: " << std::endl;
    //    generalInfo << "xSat: " << POL_operatorParams_.inputSat_ << std::endl;
    //    generalInfo << "ySat: " << POL_operatorParams_.outputSat_ << std::endl;
    //    generalInfo << "eps_mu: " << eps_mu.ToString() << std::endl;
    //    generalInfo << "## Error criteria: " << std::endl;
    //    if (POL_operatorParams_.methodType_ == VECTOR) {
    //      generalInfo << "residual tolerance wrt x: " << LM_inversion_.tolH << std::endl;
    //    }
    //    generalInfo << "residual tolerance wrt y: " << LM_inversion_.tolB << std::endl;
    //    generalInfo << std::endl;
    //    if (POL_operatorParams_.methodType_ == VECTOR) {
    //      generalInfo << "## Levenberg-Marquardt parameter set for inversion process: " << std::endl;
    //      generalInfo << "maximal number of iterations: " << LM_inversion_.maxNumIts << std::endl;
    //      generalInfo << "resolution for Jacobian: " << LM_inversion_.jacRes << std::endl;
    //      if(LM_inversion_.useTikhonov){
    //        generalInfo << "Tikhonov regularization applied" << std::endl;
    //      } else {
    //        generalInfo << "No Tikhonov regularization applied" << std::endl;
    //      }
    //      generalInfo << std::endl;
    //    }
    //    
    //    if(printStatistics){
    //      std::cout << generalInfo.str();
    //    } 
    //    
    //    LOG_TRACE(coeffcthyst) << generalInfo.str();
    //    
    //    // Tests inversion of VECTORPreisach Operator and returns 
    //    //  estimation for alpha (linesearch stepping paramater) that can be used
    //    //  for later computation
    //    // > can be used in actual simulations to setup an appropriate linesearch factor alpha
    //    //   to accelerate the actual computations later on!
    //    
    //    Vector<Double>* xIn;
    //    bool vector;
    //    bool isVirgin = true;
    //    
    //    UInt numCases = 11;
    //    UInt totalSteps;
    //    bool testAll = false;
    //    
    //    if(testNumber > numCases){
    //      testAll = true;
    //      LOG_TRACE(coeffcthyst) << "PERFORM ALL TESTS";
    //      if(printStatistics){
    //        std::cout << "PERFORM ALL TESTS" << std::endl;
    //      }
    //    } else {
    //      LOG_TRACE(coeffcthyst) << "PERFORM TEST #" << abs(testNumber);
    //      if(printStatistics){
    //        std::cout << "PERFORM TESTS #" << abs(testNumber) << std::endl;
    //      } 
    //    }
    //    
    //    for(UInt ca = 1; ca <= numCases; ca++){
    //      std::string testName;
    //			bool hystOutputRestrictedToSat = true;
    //      /*
    //       * I. Create testsignel
    //       */
    //      
    //      // recreate the hyst operator in each iteration!
    //      if (POL_operatorParams_.methodType_ == SCALAR) {
    //        vector = false;
    //        
    //        // create one temoporary ScaparPreisach operator with the same
    //        // paramter as the later used ones, except, we just need one entry   
    //        
    //        POL_useExtension_ = false;
    //        int useExtensionInt;
    //        material_->GetScalar(useExtensionInt, SCALPREISACH_USE_EXT);
    //        if(useExtensionInt == 1){
    //          POL_useExtension_ = true;
    //        }
    //        if(POL_useExtension_){
    //          material_->GetScalar(POL_operatorParams_.rotResistance_, ROT_RESISTANCE, Global::REAL);
    //          material_->GetScalar(POL_operatorParams_.angularDistance_, ANG_DISTANCE, Global::REAL);
    //          
    //          hystTMP = new ExtendedPreisach(1, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_, POL_weightParams_.weightTensor_, 
    //                  POL_operatorParams_.rotResistance_, POL_operatorParams_.angularDistance_, dim_, isVirgin, POL_weightParams_.anhysteretic_a_, POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
    //          
    //          // set initial direction
    //          Vector<Double> initialInput = Vector<Double>(dim_);
    //          POL_eps_mu_SmallSignal_.Mult(POL_operatorParams_.fixDirection_,initialInput);
    //          initialInput.ScalarMult(POL_operatorParams_.inputSat_);
    //          initialInput.Add(POL_operatorParams_.outputSat_,POL_operatorParams_.fixDirection_);
    //          // initialInput = (POL_operatorParams_.outputSat_*Identity + POL_operatorParams_.inputSat_*POL_eps_mu_SmallSignal_)*POL_operatorParams_.fixDirection_
    //          // > flux with value just at saturation
    //          hystTMP->UpdateRotationState(initialInput,POL_eps_mu_SmallSignal_,0);
    //          
    //        } else {
    //          hystTMP = new Preisach(1, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_, POL_weightParams_.weightTensor_, isVirgin, POL_weightParams_.anhysteretic_a_, POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
    //          //hyst_ = new Preisach(numHystOperators_, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_, POL_weightParams_.weightTensor_, isVirgin);
    //        }
    //        
    //        //hystTMP = new Preisach(1, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_, POL_weightParams_.weightTensor_, isVirgin);
    //        
    //				// normal preisach has worse resolution than vectorpreisach
    //				// > check for smaller tolerance
    //      } else if (POL_operatorParams_.methodName_ == "vectorPreisach_Sutor") {
    //        vector = true;
    //        // create one temporary VectorPreisach operator with the same
    //        // parameter as the later used ones, except, we just need it to store
    //        // entries for one single element / point
    //        if (POL_operatorParams_.evalVersion_ == 1) {
    //          POL_operatorParams_.isClassical_ = true; // original vector preisach model -> sutor2012
    //          
    //          hystTMP = new VectorPreisachSutor_ListApproach(1, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_,
    //                  POL_weightParams_.weightTensor_, POL_operatorParams_.rotResistance_, dim_, isVirgin,
    //                  POL_operatorParams_.isClassical_, POL_operatorParams_.angularDistance_,POL_operatorParams_.angularResolution_,POL_weightParams_.anhysteretic_a_, POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
    //        } else if (POL_operatorParams_.evalVersion_ == 2) {
    //          POL_operatorParams_.isClassical_ = false; // revised vector preisach model -> sutor2015
    //          
    //          hystTMP = new VectorPreisachSutor_ListApproach(1, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_,
    //                  POL_weightParams_.weightTensor_, POL_operatorParams_.rotResistance_, dim_, isVirgin,
    //                  POL_operatorParams_.isClassical_, POL_operatorParams_.angularDistance_,POL_operatorParams_.angularResolution_,POL_weightParams_.anhysteretic_a_, POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
    //        } else if (POL_operatorParams_.evalVersion_ == 10) {
    //          POL_operatorParams_.isClassical_ = true; // original vector preisach model -> sutor2015; matrix based implementation
    //          
    //          hystTMP = new VectorPreisachSutor_MatrixApproach(1, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_,
    //                  POL_weightParams_.weightTensor_, POL_operatorParams_.rotResistance_, dim_, isVirgin,
    //                  POL_operatorParams_.isClassical_, POL_operatorParams_.angularDistance_,POL_operatorParams_.angularResolution_,POL_weightParams_.anhysteretic_a_, POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
    //        } else if (POL_operatorParams_.evalVersion_ == 20) {
    //          POL_operatorParams_.isClassical_ = false; // revised vector preisach model -> sutor2015; matrix based implementation
    //          
    //          hystTMP = new VectorPreisachSutor_MatrixApproach(1, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_,
    //                  POL_weightParams_.weightTensor_, POL_operatorParams_.rotResistance_, dim_, isVirgin,
    //                  POL_operatorParams_.isClassical_, POL_operatorParams_.angularDistance_,POL_operatorParams_.angularResolution_,POL_weightParams_.anhysteretic_a_, POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
    //        } else {
    //          EXCEPTION("POL_operatorParams_.evalVersion_ has to be one of the following: \n "
    //                  "1: classical vector model (sutor2012) \n"
    //                  "2: revised vector model (sutor2015) [DEFAULT] \n"
    //                  "10: classical vector model (sutor2012) - Matrix implementation, only for reference \n"
    //                  "20: revised vector model (sutor2015) - Matrix implementation, only for reference \n")
    //        }
    //        hystTMP->SetParamsForInversion(LM_inversion_.maxNumIts, LM_inversion_.tolH, LM_inversion_.tolB, LM_inversion_.jacRes,
    //                LM_inversion_.useTikhonov, LM_inversion_.alphaLSStart,LM_inversion_.alphaLSMin,LM_inversion_.alphaLSMax,POL_operatorParams_.angularClipping_); 
    //      } else if (POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") {
    //        vector = true;
    //        // basically a scalar model in multiple directions
    //        // isotropic case: all scalar models are equal (same weights etc)
    //        // anisotropic case: each model different; choice of directions matters; weights are harder to obtain
    //        int isIsotropic = 1;
    //        material_->GetScalar(isIsotropic, PREISACH_MAYERGOYZ_ISOTROPIC);
    //        if( (dim_ != 2) || (isIsotropic == 0)){
    //          EXCEPTION("Mayergoyz vector model currently only implemented for 2d isotropic materials");
    //        }
    //        int numDirections = 0;
    //        material_->GetScalar(numDirections, PREISACH_MAYERGOYZ_NUM_DIR);
    //        
    //        int clipOutput = 0;
    //        material_->GetScalar(clipOutput, PREISACH_MAYERGOYZ_CLIPOUTPUT);
    //        
    //				if(clipOutput == 0){
    //					hystOutputRestrictedToSat = false;
    //				}
    //				
    //        /*
    //         * IMPORTANT REMARK:
    //         *  > although the Mayergoyz model is based on the scalar models 
    //         *     we are not alloewed to directly apply the preisach parameter for
    //         *     the scalar case (i.e. the weights, the anhyst parameter and so on)
    //         *  > make sure that the passed parameter are already transformed correctly
    //         *      > see constructor above
    //         */
    //        hystTMP = new VectorPreisachMayergoyz(1, numDirections, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_, 
    //                POL_weightParams_.weightTensor_,dim_,isVirgin,POL_weightParams_.anhysteretic_a_, POL_weightParams_.anhysteretic_b_, 
    //								POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_,clipOutput);
    //        
    //        hystTMP->SetParamsForInversion(LM_inversion_.maxNumIts, LM_inversion_.tolH, LM_inversion_.tolB, LM_inversion_.jacRes,
    //                LM_inversion_.useTikhonov, LM_inversion_.alphaLSStart,LM_inversion_.alphaLSMin,LM_inversion_.alphaLSMax,POL_operatorParams_.angularClipping_); 
    //      } else {
    //        EXCEPTION("Invalid model selected for inversion test");
    //      }
    //      
    //      bool fieldsAlignedAboveSat = true;
    //      if (POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz"){
    //        fieldsAlignedAboveSat = false;
    //      }
    //      
    //      
    //      if( (ca == 1)&&( (testNumber == 1) ||(testAll)) ) {
    //        if(printStatistics){
    //          std::cout << "##### TEST CASE " << ca << " - Self designed test signal" << std::endl;
    //        }
    //        testName = "Self designed test signal";
    //        LOG_TRACE(coeffcthyst) << "##### TEST CASE " << ca << " - Self designed test signal";
    //        /*
    //         * Testcase 1 --- (self-designed) benchmark signal
    //         * 
    //         * 1d-Signal shape:
    //         *  exp-increase, cos decrease, linear increase, 1/x decrease, log increase, x^3 increase, hold state
    //         * 
    //         * 2d/3d-Signal shape:
    //         *  Multiply 1d-Signal with chaning direction vector
    //         * 
    //         */
    //        
    //        UInt steps1,steps2,steps3,steps4,steps5,steps6,steps7;
    //        steps1=50;
    //        steps2=50;
    //        steps3=50;
    //        steps4=50;
    //        steps5=50;
    //        steps6=50;
    //        steps7=10;
    //        totalSteps=steps1+steps2+steps3+steps4+steps5+steps6+steps7;
    //        
    //        xIn = new Vector<Double>[totalSteps]; 
    //        Vector<Double> xScal = Vector<Double>(totalSteps);
    //        xScal.Init();
    //        Double delta;
    //        
    //        // 1. expontential increase from [10^-15 to 10^0[
    //        delta = 15.0/steps1;
    //        
    //        for(UInt i = 0; i < steps1; i++){
    //          xScal[i] = std::pow(10,-15.0+i*delta);
    //        }
    //        // 2. cos decrease towards -0.5[
    //        delta = -std::acos(-0.5)/steps2;
    //        
    //        for(UInt i = 0; i < steps2; i++){
    //          xScal[steps1+i] = std::cos(i*delta);
    //        }
    //        // 3. linear increase to 0.4
    //        delta = 0.9/steps3;
    //        
    //        for(UInt i = 0; i < steps3; i++){
    //          xScal[steps1+steps2+i] = -0.5+i*delta;
    //        }
    //        // 4. 1/x decrease towards -0.3
    //        Double Xstart = 1.0/(0.7 + 1.0/steps4);
    //        delta = (steps4 - Xstart)/steps4;
    //        
    //        for(UInt i = 0; i < steps4; i++){
    //          xScal[steps1+steps2+steps3+i] = 1.0/(Xstart+i*delta)-1.0/steps4-0.3;
    //        }
    //        // 5. logarithmic increase towards 0.2
    //        Xstart = std::exp(-3);
    //        Double Xend = std::exp(2);
    //        delta = (Xend-Xstart)/steps5;
    //        
    //        for(UInt i = 0; i < steps5; i++){
    //          xScal[steps1+steps2+steps3+steps4+i] = 0.1*std::log(Xstart+i*delta);
    //        }
    //        // 6. x^3 increase to 0.8
    //        Double denom = std::pow(steps6-1,3);
    //        for(UInt i = 0; i < steps6; i++){
    //          xScal[steps1+steps2+steps3+steps4+steps5+i] = 0.6*std::pow(i,3)/denom + 0.2;
    //        }
    //        
    //        // 7. hold signal at value 0.8
    //        for(UInt i = 0; i < steps7; i++){
    //          xScal[steps1+steps2+steps3+steps4+steps5+steps6+i] = 0.8;
    //        }
    //        
    //        //        std::ofstream testOutput;
    //        //        testOutput.open("GeneratedInput");
    //        //        for(UInt i = 0; i < totalSteps; i++){
    //        //          testOutput << i << " " << xScal[i] << std::endl;
    //        //        }
    //        //        testOutput.close();
    //        Vector<Double> dir = Vector<Double>(dim_);
    //        dir[0] = 1.0;
    //        dir[1] = 0.0;
    //        for(UInt i = 0; i < totalSteps; i++){
    //          xIn[i] = Vector<Double>(dim_);
    //          xIn[i].Init();
    //          dir[0] = (Double) (totalSteps-i)/totalSteps;
    //          if(vector){
    //            dir[1] = (Double) i/totalSteps;
    //          } else {
    //            dir[1] = 0.0;
    //          }
    //          xIn[i].Add(POL_operatorParams_.inputSat_*xScal[i],dir);
    //        }
    //        
    //      } else if( (ca == 2)&&( (testNumber == 2) ||(testAll)) ) {
    //        if(printStatistics){
    //          std::cout << "##### TEST CASE " << ca << " - Decreasing Sawtooth" << std::endl;
    //        }
    //        LOG_TRACE(coeffcthyst) << "##### TEST CASE " << ca << " - Decreasing Sawtooth";
    //        testName = "Decreasing Sawtooth";
    //        totalSteps = 200;
    //        xIn = new Vector<Double>[totalSteps]; 
    //        /*
    //         * Testcase: Put material into full saturation in y-direction; from that state on
    //         * apply a decreasing triangular signal in x-direction; for initial period hold
    //         * signla in y-direction; during first real period decrease to 1/2 of value;
    //         * then to 1/4 and finally to 0
    //         */
    //        xIn[0] = Vector<Double>(dim_);
    //        xIn[0].Init(0.0);
    //        xIn[0][1] = -POL_operatorParams_.inputSat_;
    //        
    //        UInt numStepsFirstIncrease = (UInt) ( (Double)totalSteps/8.0 ) +1;
    //        Double incr = POL_operatorParams_.inputSat_/((Double) numStepsFirstIncrease-1);
    //        
    //        for(UInt i = 1; i < numStepsFirstIncrease; i++){
    //          xIn[i] = Vector<Double>(dim_);
    //          xIn[i].Init(0.0);
    //          xIn[i][0] = xIn[i-1][0] + incr; 
    //          if(vector){
    //            xIn[i][1] = -POL_operatorParams_.inputSat_;
    //          }
    //        }
    //        
    //        UInt cnt = numStepsFirstIncrease;
    //        
    //        UInt remainingSteps = totalSteps - numStepsFirstIncrease;
    //        UInt stepsPerPeriod = 20;
    //        UInt numPeriods = remainingSteps/stepsPerPeriod;
    //        remainingSteps = remainingSteps%stepsPerPeriod;
    //        
    //        incr = 2*POL_operatorParams_.inputSat_/((Double)(stepsPerPeriod-1));
    //        Double sign = -1.0;
    //        Double decrFactor = 1.0 - 1.0/((Double) numPeriods);
    //        
    //        for(UInt j = 0; j < numPeriods; j++){
    //          for(UInt i = 0; i < stepsPerPeriod; i++){
    //            xIn[cnt] = Vector<Double>(dim_);
    //            xIn[cnt].Init(0.0);
    //            xIn[cnt][0] = xIn[cnt-1][0] + sign*incr;
    //            
    //            if(vector){
    //              if(j==0){
    //                xIn[cnt][1] = -POL_operatorParams_.inputSat_/2.0;
    //              } else if(j==1){
    //                xIn[cnt][1] = -POL_operatorParams_.inputSat_/4.0;
    //              }
    //            }
    //            
    //            cnt++;
    //          }
    //          sign = -1.0*sign;
    //          incr = incr*decrFactor;
    //        }
    //        for(UInt i = 0; i < remainingSteps; i++){
    //          xIn[cnt] = Vector<Double>(dim_);
    //          xIn[cnt].Init(0.0);
    //          xIn[cnt][0] = xIn[cnt-1][0] + sign*incr;
    //          cnt++;
    //        }
    //      } else if( (ca == 3)&&( (testNumber == 3) ||(testAll)) ) {
    //        if(printStatistics){
    //          std::cout << "##### TEST CASE " << ca << " - Decreasing Sin" << std::endl;
    //        }
    //        LOG_TRACE(coeffcthyst) << "##### TEST CASE " << ca << " - Decreasing Sin";
    //        testName = "Decreasing Sin";
    //        totalSteps = 200;
    //        Double numPeriods = 4;
    //        Double stepsPerPeriod = totalSteps/numPeriods;
    //        Double decrease = 1.0/totalSteps;
    //        xIn = new Vector<Double>[totalSteps]; 
    //        
    //        
    //        for(UInt i = 0; i < totalSteps; i++){
    //          xIn[i] = Vector<Double>(dim_);
    //          xIn[i].Init(0.0);
    //          //xIn[i][0] = sin( (2*M_PI*i)/stepsPerPeriod );
    //          xIn[i][0] = 1.2*POL_operatorParams_.inputSat_*(1.0 - decrease*i) * sin( (2*M_PI*i)/stepsPerPeriod );
    //          
    //          if(vector){
    //            if(i >= 3*stepsPerPeriod){
    //              xIn[i][1] = -POL_operatorParams_.inputSat_/8.0;
    //            } else if(i >= 2*stepsPerPeriod){
    //              xIn[i][1] = -POL_operatorParams_.inputSat_/4.0;
    //            } else if(i >= stepsPerPeriod){
    //              xIn[i][1] = -POL_operatorParams_.inputSat_/2.0;
    //            } else if(i >= 0){
    //              xIn[i][1] = -POL_operatorParams_.inputSat_;
    //            }
    //          }
    //        }
    //      } else if( (ca == 4)&&( (testNumber == 4) ||(testAll)) ) {
    //        if(printStatistics){
    //          std::cout << "##### TEST CASE " << ca << " - Decreasing Sin - extended" << std::endl;
    //        }
    //        LOG_TRACE(coeffcthyst) << "##### TEST CASE " << ca << " - Decreasing Sin - extended";
    //        testName = "Decreasing Sin - Extended";
    //        totalSteps = 500;
    //        Double numPeriods = 10;
    //        Double stepsPerPeriod = totalSteps/numPeriods;
    //        Double decrease = 1.0/totalSteps;
    //        xIn = new Vector<Double>[totalSteps]; 
    //        
    //        
    //        for(UInt i = 0; i < totalSteps; i++){
    //          xIn[i] = Vector<Double>(dim_);
    //          xIn[i].Init(0.0);
    //          //xIn[i][0] = sin( (2*M_PI*i)/stepsPerPeriod );
    //          xIn[i][0] = 0.6*POL_operatorParams_.inputSat_*(1.0 - decrease*i) * sin( (2*M_PI*i)/stepsPerPeriod );
    //          if(vector){
    //            // smaller amplitude, higher frequency
    //            xIn[i][1] = -0.4*POL_operatorParams_.inputSat_*(1.0 - decrease*i) * sin( 3*(2*M_PI*i)/stepsPerPeriod );
    //          }
    //        }
    //      } else if( (ca == 5)&&( (testNumber == 5) ||(testAll)) ) {
    //        if(printStatistics){
    //          std::cout << "##### TEST CASE " << ca << " - Decreasing Sin - extended, scaled down" << std::endl;
    //        }
    //        LOG_TRACE(coeffcthyst) << "##### TEST CASE " << ca << " - Decreasing Sin - extended, scaled down";
    //        testName = "Decreasing Sin - Extended";
    //        totalSteps = 500;
    //        Double numPeriods = 10;
    //        Double stepsPerPeriod = totalSteps/numPeriods;
    //        Double decrease = 1.0/totalSteps;
    //        xIn = new Vector<Double>[totalSteps]; 
    //        
    //        
    //        for(UInt i = 0; i < totalSteps; i++){
    //          xIn[i] = Vector<Double>(dim_);
    //          xIn[i].Init(0.0);
    //          //xIn[i][0] = sin( (2*M_PI*i)/stepsPerPeriod );
    //          xIn[i][0] = 0.006*POL_operatorParams_.inputSat_*(1.0 - decrease*i) * sin( (2*M_PI*i)/stepsPerPeriod );
    //          if(vector){
    //            // smaller amplitude, higher frequency
    //            xIn[i][1] = -0.004*POL_operatorParams_.inputSat_*(1.0 - decrease*i) * sin( 3*(2*M_PI*i)/stepsPerPeriod );
    //          }
    //        }
    //      } else if( (ca == 6)&&( (testNumber == 6) ||(testAll)) ) {
    //        if(printStatistics){
    //          std::cout << "##### TEST CASE " << ca << " - Decreasing Sin - extended, reverse" << std::endl;
    //        }
    //        LOG_TRACE(coeffcthyst) << "##### TEST CASE " << ca << " - Decreasing Sin - extended, reverse";
    //        testName = "Decreasing Sin - Extended";
    //        totalSteps = 500;
    //        Double numPeriods = 10;
    //        Double stepsPerPeriod = totalSteps/numPeriods;
    //        Double increase = 1.0/totalSteps;
    //        xIn = new Vector<Double>[totalSteps]; 
    //        
    //        
    //        for(UInt i = 0; i < totalSteps; i++){
    //          xIn[i] = Vector<Double>(dim_);
    //          xIn[i].Init(0.0);
    //          //xIn[i][0] = sin( (2*M_PI*i)/stepsPerPeriod );
    //          xIn[i][0] = 0.6*POL_operatorParams_.inputSat_*increase*i * sin( (2*M_PI*i)/stepsPerPeriod );
    //          if(vector){
    //            // smaller amplitude, higher frequency
    //            xIn[i][1] = -0.4*POL_operatorParams_.inputSat_*increase*i * sin( 3*(2*M_PI*i)/stepsPerPeriod );
    //          }
    //        }
    //      } else if( (ca == 7)&&( (testNumber == 7) ||(testAll)) ) {
    //        if(printStatistics){
    //          std::cout << "##### TEST CASE " << ca << " - Forc" << std::endl;
    //        }
    //        LOG_TRACE(coeffcthyst) << "##### TEST CASE " << ca << " - Forc";
    //        testName = "Forc";
    //        totalSteps = 500;
    //        Double numPeriods = 10;
    //        Double stepsPerPeriod = totalSteps/numPeriods;
    //        Double decrease = 1.0/totalSteps;
    //        xIn = new Vector<Double>[totalSteps]; 
    //        
    //        for(UInt i = 0; i < totalSteps; i++){
    //          xIn[i] = Vector<Double>(dim_);
    //          xIn[i].Init(0.0);
    //          //xIn[i][0] = sin( (2*M_PI*i)/stepsPerPeriod );
    //          xIn[i][0] = 1.2*POL_operatorParams_.inputSat_*(1.0 - decrease*i) * sin( (2*M_PI*i)/stepsPerPeriod ) - 1.2*POL_operatorParams_.inputSat_*decrease*i;
    //          if(vector){
    //            // smaller amplitude, higher frequency
    //            xIn[i][1] = 0.4*POL_operatorParams_.inputSat_*(1.0 - decrease*i) * sin( 3*(2*M_PI*i)/stepsPerPeriod ) - 0.4*POL_operatorParams_.inputSat_*decrease*i;
    //          }
    //        }
    //      } else if( (ca == 8)&&( (testNumber == 8) ||(testAll)) ) {
    //        if(printStatistics){
    //          std::cout << "##### TEST CASE " << ca << " - Saturation in x, Drop to rem, Saturation in y, Drop to rem" << std::endl;
    //        }
    //        LOG_TRACE(coeffcthyst) << "##### TEST CASE " << ca << " - SatX,SatY";
    //        testName = "SatX,SatY";
    //        totalSteps = 500;
    //				
    //				Double increase = 4.0/totalSteps;
    //				
    //        xIn = new Vector<Double>[totalSteps]; 
    //				Double maxAmplitude = 4.0*POL_operatorParams_.inputSat_;
    //        
    //        for(UInt i = 0; i < totalSteps/4; i++){
    //					// linearly increase in x up to n*saturation
    //          xIn[i] = Vector<Double>(dim_);
    //          xIn[i].Init(0.0);
    //					
    //					xIn[i][0] = maxAmplitude*i*increase;
    //        }
    //				for(UInt i = totalSteps/4; i < 2*totalSteps/4; i++){
    //					// linearly decrease to 0
    //					xIn[i] = Vector<Double>(dim_);
    //					xIn[i].Init(0.0);
    //					
    //					xIn[i][0] = maxAmplitude - maxAmplitude*(i-totalSteps/4)*increase;
    //				}
    //				
    //				for(UInt i = 2*totalSteps/4; i < 3*totalSteps/4; i++){
    //					// linearly increase in y up to n*saturation
    //          xIn[i] = Vector<Double>(dim_);
    //          xIn[i].Init(0.0);
    //					if(vector){
    //            xIn[i][1] = maxAmplitude*(i-2*totalSteps/4)*increase;
    //          }
    //        }
    //				for(UInt i = 3*totalSteps/4; i < totalSteps; i++){
    //					// linearly decrease to 0
    //					xIn[i] = Vector<Double>(dim_);
    //					xIn[i].Init(0.0);
    //					if(vector){
    //            xIn[i][1] = maxAmplitude - maxAmplitude*(i-3*totalSteps/4)*increase;
    //          }
    //				}
    //				
    //				
    //      } else if( (ca == 9)&&( (testNumber == 9) ||(testAll)) ) {
    //        if(printStatistics){
    //          std::cout << "##### TEST CASE " << ca << " - exactly Saturation in x, Drop to rem, exactly Saturation in y, Drop to rem" << std::endl;
    //        }
    //        LOG_TRACE(coeffcthyst) << "##### TEST CASE " << ca << " - SatX,SatY clipped";
    //        testName = "SatX,SatY clipped";
    //        totalSteps = 500;
    //				
    //				Double increase = 4.0/totalSteps;
    //				
    //        xIn = new Vector<Double>[totalSteps]; 
    //				Double maxAmplitude = 4.0*POL_operatorParams_.inputSat_;
    //        
    //        for(UInt i = 0; i < totalSteps/4; i++){
    //					// linearly increase in x up to n*saturation
    //          xIn[i] = Vector<Double>(dim_);
    //          xIn[i].Init(0.0);
    //          
    //					xIn[i][0] = maxAmplitude*i*increase;
    //          if(xIn[i][0] > POL_operatorParams_.inputSat_){
    //            xIn[i][0] = POL_operatorParams_.inputSat_;
    //          }
    //        }
    //				for(UInt i = totalSteps/4; i < 2*totalSteps/4; i++){
    //					// linearly decrease to 0
    //					xIn[i] = Vector<Double>(dim_);
    //					xIn[i].Init(0.0);
    //          
    //					xIn[i][0] = maxAmplitude - maxAmplitude*(i-totalSteps/4)*increase;
    //          if(xIn[i][0] > POL_operatorParams_.inputSat_){
    //            xIn[i][0] = POL_operatorParams_.inputSat_;
    //          }
    //				}
    //				
    //				for(UInt i = 2*totalSteps/4; i < 3*totalSteps/4; i++){
    //					// linearly increase in y up to n*saturation
    //          xIn[i] = Vector<Double>(dim_);
    //          xIn[i].Init(0.0);
    //          
    //          if(vector){
    //            xIn[i][1] = maxAmplitude*(i-2*totalSteps/4)*increase;
    //          }
    //          if(xIn[i][1] > POL_operatorParams_.inputSat_){
    //            xIn[i][1] = POL_operatorParams_.inputSat_;
    //          }
    //        }
    //				for(UInt i = 3*totalSteps/4; i < totalSteps; i++){
    //					// linearly decrease to 0
    //					xIn[i] = Vector<Double>(dim_);
    //					xIn[i].Init(0.0);
    //          
    //          if(vector){
    //            xIn[i][1] = maxAmplitude - maxAmplitude*(i-3*totalSteps/4)*increase;
    //          }
    //          if(xIn[i][1] > POL_operatorParams_.inputSat_){
    //            xIn[i][1] = POL_operatorParams_.inputSat_;
    //          }
    //				}
    //        
    //      } else if( (ca == 10)&&( (testNumber == 10) ||(testAll)) ) {
    //        if(printStatistics){
    //          std::cout << "##### TEST CASE " << ca << " - Forc, Short" << std::endl;
    //        }
    //        LOG_TRACE(coeffcthyst) << "##### TEST CASE " << ca << " - Forc, Short";
    //        testName = "Forc";
    //        totalSteps = 80;
    //        Double numPeriods = 2;
    //        Double stepsPerPeriod = totalSteps/numPeriods;
    //        Double decrease = 1.0/totalSteps;
    //        xIn = new Vector<Double>[totalSteps]; 
    //        
    //        for(UInt i = 0; i < totalSteps; i++){
    //          xIn[i] = Vector<Double>(dim_);
    //          xIn[i].Init(0.0);
    //          //xIn[i][0] = sin( (2*M_PI*i)/stepsPerPeriod );
    //          xIn[i][0] = 1.2*POL_operatorParams_.inputSat_*(1.0 - decrease*i) * sin( (2*M_PI*i)/stepsPerPeriod ) - 1.2*POL_operatorParams_.inputSat_*decrease*i;
    //          if(vector){
    //            // smaller amplitude, higher frequency
    //            xIn[i][1] = 0.4*POL_operatorParams_.inputSat_*(1.0 - decrease*i) * sin( 3*(2*M_PI*i)/stepsPerPeriod ) - 0.4*POL_operatorParams_.inputSat_*decrease*i;
    //          }
    //        }
    //      } else if( (ca == 11)&&( (testNumber == 11) ||(testAll)) ) {
    //        if(printStatistics){
    //          std::cout << "##### TEST CASE " << ca << " - 1.5*Saturation in x, Drop to rem, 1.5*Saturation in y, Drop to rem" << std::endl;
    //        }
    //        LOG_TRACE(coeffcthyst) << "##### TEST CASE " << ca << " - 1.5*SatX,1.5*SatY";
    //        testName = "1.5*SatX,1.5*SatY";
    //        totalSteps = 200;
    //				
    //				Double increase = 4.0/totalSteps;
    //				
    //        xIn = new Vector<Double>[totalSteps]; 
    //				Double maxAmplitude = 1.5*POL_operatorParams_.inputSat_;
    //        
    //        for(UInt i = 0; i < totalSteps/4; i++){
    //					// linearly increase in x up to n*saturation
    //          xIn[i] = Vector<Double>(dim_);
    //          xIn[i].Init(0.0);
    //					
    //					xIn[i][0] = maxAmplitude*i*increase;
    //        }
    //				for(UInt i = totalSteps/4; i < 2*totalSteps/4; i++){
    //					// linearly decrease to 0
    //					xIn[i] = Vector<Double>(dim_);
    //					xIn[i].Init(0.0);
    //					
    //					xIn[i][0] = maxAmplitude - maxAmplitude*(i-totalSteps/4)*increase;
    //				}
    //				
    //				for(UInt i = 2*totalSteps/4; i < 3*totalSteps/4; i++){
    //					// linearly increase in y up to n*saturation
    //          xIn[i] = Vector<Double>(dim_);
    //          xIn[i].Init(0.0);
    //					if(vector){
    //            xIn[i][1] = maxAmplitude*(i-2*totalSteps/4)*increase;
    //          }
    //        }
    //				for(UInt i = 3*totalSteps/4; i < totalSteps; i++){
    //					// linearly decrease to 0
    //					xIn[i] = Vector<Double>(dim_);
    //					xIn[i].Init(0.0);
    //					if(vector){
    //            xIn[i][1] = maxAmplitude - maxAmplitude*(i-3*totalSteps/4)*increase;
    //          }
    //				}
    //				
    //				
    //      } else {
    //        continue;
    //      }
    //      
    //			
    //			
    //			
    //			
    //			
    //			
    //      /*
    //       * II. Setup output streams
    //       */
    //      
    //      std::stringstream stringstream1;
    //      stringstream1 << "INVERSION_TEST_xIn_" << ca;
    //      std::string filenameInput = stringstream1.str();
    //      std::stringstream stringstream2;
    //      stringstream2 << "INVERSION_TEST_xRet_" << ca;
    //      std::string filenameOutput = stringstream2.str();
    //      
    //      std::stringstream stringstream4;
    //      stringstream4 << "INVERSION_TEST_xIn_hIn_" << ca;
    //      std::string filenameOutputHystIn = stringstream4.str();
    //      std::stringstream stringstream5;
    //      stringstream5 << "INVERSION_TEST_xRet_hRet_" << ca;
    //      std::string filenameOutputHystRet = stringstream5.str();
    //      
    //      std::stringstream stringstream6;
    //      stringstream6 << "LOG_" << ca;
    //      std::string filenameLog = stringstream6.str();
    //      
    //      
    //      //      std::ofstream testOutputSig;
    //      //      if(writeResultsToFiles){
    //      //        testOutputSig.open(filenameInput);
    //      //        testOutputSig << "## Results of TestInversion for test signal " << ca << " - " << testName << std::endl; 
    //      //        testOutputSig << "# Number of steps: " << totalSteps << std::endl;
    //      //        testOutputSig << "# eps_mu: " << eps_mu.ToString() << std::endl;
    //      //        testOutputSig << "# xSat: " << POL_operatorParams_.inputSat_ << std::endl;
    //      //        testOutputSig << "# ySat: " << POL_operatorParams_.outputSat_ << std::endl;
    //      //        testOutputSig << "# Starting value of inversion " << xIn[0][0] << " " << xIn[0][1] << std::endl;
    //      //        testOutputSig << "# " << std::endl;
    //      //        testOutputSig << "# Step    x_in[0]    x_in[1]" << std::endl;
    //      //        for(UInt i = 0; i < totalSteps; i++){
    //      //          testOutputSig << i << " " << xIn[i][0] << " " << xIn[i][1] << std::endl;
    //      //        }
    //      //        testOutputSig.close();
    //      //      }
    //      
    //      Vector<int> successCodeVector = Vector<int>(totalSteps);
    //      Vector<Double>* xRetrieved = new Vector<Double>[totalSteps];
    //      Vector<Double>* yRetrieved = new Vector<Double>[totalSteps];
    //      Vector<Double>* hIn= new Vector<Double>[totalSteps];
    //      Vector<Double>* hRetrieved = new Vector<Double>[totalSteps];
    //      Vector<Double>* yError = new Vector<Double>[totalSteps];
    //      Vector<Double>* xError = new Vector<Double>[totalSteps];
    //      Vector<Double>* yIn = new Vector<Double>[totalSteps]; 
    //      
    //      //std::ofstream testOutput2;
    //      std::ofstream testOutput3;
    //      std::ofstream testOutput4;
    //      std::ofstream testOutput5;
    //      
    //      if(writeResultsToFiles){
    //        //testOutput2.open(filenameOutput);
    //        testOutput3.open(filenameOutputHystIn);
    //        testOutput4.open(filenameOutputHystRet);
    //        
    //        //        testOutput2 << "## Results of TestInversion for test signal " << ca << " - " << testName << std::endl; 
    //        //        testOutput2 << "# Number of steps: " << totalSteps << std::endl;
    //        //        testOutput2 << "# eps_mu: " << eps_mu.ToString() << std::endl;
    //        //        testOutput2 << "# xSat: " << POL_operatorParams_.inputSat_ << std::endl;
    //        //        testOutput2 << "# ySat: " << POL_operatorParams_.outputSat_ << std::endl;
    //        //        testOutput2 << "# Starting value of inversion " << xIn[0][0] << " " << xIn[0][1] << std::endl;
    //        //        testOutput2 << "# " << std::endl;
    //        //        testOutput2 << "# Step    x_ret[0]    x_ret[1]" << std::endl;
    //        //        
    //        testOutput3 << "## Results of TestInversion for test signal " << ca << " - " << testName << std::endl; 
    //        testOutput3 << "# Number of steps: " << totalSteps << std::endl;
    //        testOutput3 << "# eps_mu: " << eps_mu.ToString() << std::endl;
    //        testOutput3 << "# xSat: " << POL_operatorParams_.inputSat_ << std::endl;
    //        testOutput3 << "# ySat: " << POL_operatorParams_.outputSat_ << std::endl;
    //        testOutput3 << "# Starting value of inversion " << xIn[0][0] << " " << xIn[0][1] << std::endl;
    //        testOutput3 << "# " << std::endl;
    //        testOutput3 << "# Step    x_in[0]    x_in[1]    h_in[0]    h_in[1]    y_in[0]    y_in[1]" << std::endl;
    //        
    //        testOutput4 << "## Results of TestInversion for test signal " << ca << " - " << testName << std::endl; 
    //        testOutput4 << "# Number of steps: " << totalSteps << std::endl;
    //        testOutput4 << "# eps_mu: " << eps_mu.ToString() << std::endl;
    //        testOutput4 << "# xSat: " << POL_operatorParams_.inputSat_ << std::endl;
    //        testOutput4 << "# ySat: " << POL_operatorParams_.outputSat_ << std::endl;
    //        testOutput4 << "# Starting value of inversion " << xIn[0][0] << " " << xIn[0][1] << std::endl;
    //        testOutput4 << "# " << std::endl;
    //        testOutput4 << "# Step    x_ret[0]    x_ret[1]    h_ret[0]    h_ret[1]    y_ret[0]    y_ret[1]" << std::endl;
    //      }
    //      
    //      
    //      /*
    //       * III. Start testing
    //       */
    //      
    //      std::map<UInt, Double> failedTests;
    //      
    //      Double minAlpha = 0;
    //      Double maxAlpha = 0;
    //      Double avgAlpha = 0;
    //      Double totalminAlpha = 1e10;
    //      Double totalmaxAlpha = -1e10;
    //      Double totalavgAlpha = 0;
    //      
    //			UInt numberOfLMIterations = 0;
    //			UInt numberOfLinesearchIterations = 0;
    //			UInt maxNumberOfLinesearchIterations = 0;
    //      
    //      // counter for forward evaluation
    //      UInt forwardFails = 0;
    //      UInt forwardReused = 0;
    //      UInt forwardAnhystOnly = 0;
    //      UInt forwardOverwrite = 0;
    //      UInt forwardTMP = 0;
    //
    //      // counter for backward evaluation / inversion
    //      UInt LMFails = 0;
    //      UInt totalReused = 0;
    //      UInt totalAnhystOnly = 0;
    //      UInt totalBisection = 0;
    //      UInt totalRemanence = 0;
    //      UInt totalPassedErrorTol = 0;
    //      UInt totalPassedResTolX = 0;
    //      UInt totalPassedResTolY = 0;
    //			UInt totalNumberOfLMIterations = 0;
    //			UInt totalNumberOfLinesearchIterations = 0;
    //			UInt absmaxNumberOfLinesearchIterations = 0;
    //
    //      // overwriteMemory: 
    //      // this was true all the time but I don't get why
    //      // when we evaluate the hyst operator we just want to know which value
    //      // of y we have to solve for but we do not want to set the material already
    //      // if we would do so, the material would already have the solution inside
    //      bool overwriteMemory = false;
    //      bool overwriteDirection = true;
    //      
    //      yIn[0] = Vector<Double>(dim_);
    //      yIn[0].Init(0.0);
    //      hIn[0] = Vector<Double>(dim_);
    //      hIn[0].Init(0.0);
    //      
    //      //      if(POL_operatorParams_.angularClipping_ > 0){
    //      //        /*
    //      //         * clip input to hyst operator; clip output of it
    //      //         */
    //      //        ClipDirection(xIn[0]);
    //      //      }
    //      
    //      Timer* forwardTimer = new Timer();
    //      Double forwardMaxEvalTime = 0;
    //      Double forwardAvgEvalTime = 0;
    //      Double forwardTotalEvalTime = 0;
    //      UInt forwardEvalCounter = 0;
    //      
    //      Timer* backwardTimer = new Timer();
    //      Double backwardMaxEvalTime = 0;
    //      Double backwardAvgEvalTime = 0;
    //      Double backwardTotalEvalTime = 0;
    //      UInt backwardEvalCounter = 0;
    //      Double startTime = 0;
    //      Double endTime = 0;
    //      Double evalTime = 0;
    //        
    //      int successFlagForward = -1;
    //      // forward:
    //      // -1 = fail
    //      //  0 = reuse value
    //      //  1 = anhyst only
    //      //  2 = eval on permanent storage
    //      //  3 = eval on temporal storage
    //      
    //      int successFlagBackward = -1;
    //      // backward:
    //      // -1 = fail
    //      //  0 = reuse value
    //      //  1 = anhyst only
    //      //  2 = bisection
    //      //  3-6 only for vector implementation using Levenberg Marquardt
    //      //  3 = reamnence
    //      //  4 = passed dut to error tolerance
    //      //  5 = passed due to tolerance wrt x
    //      //  6 = passed due to tolerance wrt y
    //      
    //      bool debugOut = false;
    //      
    ////      if (XML_performanceMeasurement_ == 3) {
    ////        startTime = forwardTimer->GetCPUTime();
    ////				forwardTimer->Start();
    ////			}
    ////
    ////      // Get polarization P
    ////      
    ////      if(vector){
    ////        hIn[0] = hystTMP->computeValue_vec(xIn[0], 0, overwriteMemory, overwriteDirection,debugOut,successFlagForward);
    ////      } else {
    ////				//hIn[0][0] = hystTMP->computeValueAndUpdate(xIn[0][0], 0, 0, overwriteMemory);
    ////        hIn[0][0] = hystTMP->computeValueAndUpdate(xIn[0][0], 0, overwriteMemory,successFlagForward);
    ////      }
    ////      
    ////      if (XML_performanceMeasurement_ == 3){
    ////        forwardTimer->Stop();
    ////      
    ////        // only count actual computations; if value gets reused, this should not count
    ////        if(successFlagForward != 0) {
    ////          forwardEvalCounter++;
    ////				
    ////          endTime = forwardTimer->GetCPUTime();
    ////          evalTime = endTime-startTime;
    ////          if(evalTime > forwardMaxEvalTime){
    ////            forwardMaxEvalTime = evalTime;
    ////          }
    ////          forwardTotalEvalTime += evalTime;
    ////        }
    ////			}
    ////
    ////      if(successFlagForward == -1){
    ////        forwardFails++;
    ////      } else if(successFlagForward == 0){
    ////        forwardReused++;
    ////      } else if(successFlagForward == 1){
    ////        forwardAnhystOnly++;
    ////      } else if(successFlagForward == 2){
    ////        forwardOverwrite++;
    ////      } else if(successFlagForward == 3){
    ////        forwardTMP++; 
    ////      }
    ////      
    ////      //      if(POL_operatorParams_.angularClipping_ > 0){
    ////      //        /*
    ////      //         * clip input to hyst operator; clip output of it
    ////      //         */
    ////      //        ClipDirection(hIn[0]);
    ////      //      }
    ////      
    ////      // D = eps*E + P
    ////      yIn[0].Add(1.0,hIn[0]);
    ////      for(UInt j = 0; j < dim_; j++){
    ////        yIn[0][j] += eps_mu[j][j]*xIn[0][j];
    ////      }
    ////      
    ////      xRetrieved[0] = Vector<Double>(dim_);
    ////      xRetrieved[0].Init(0.0);
    ////      
    ////			Vector<Double> zeroVec = Vector<Double>(dim_);
    ////			zeroVec.Init();
    ////      
    ////      if (XML_performanceMeasurement_ == 3) {
    ////        
    ////        startTime = backwardTimer->GetCPUTime();
    ////        backwardTimer->Start();
    ////      }
    ////      
    ////      if(vector){
    ////        //        if (POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz"){
    ////        //          xRetrieved[0] = hystTMP->computeInput_vec(yIn[0], 0, eps_mu,overwriteDirection);
    ////        //        } else {
    ////        // new: we have to pass the previous solution (here: 0)
    ////        // new2: pass arguments by value!
    ////        xRetrieved[0] = hystTMP->computeInput_vec_withStatistics(yIn[0], zeroVec, zeroVec, zeroVec, 0, eps_mu,
    ////                overwriteDirection, fieldsAlignedAboveSat, hystOutputRestrictedToSat, numberOfLMIterations, numberOfLinesearchIterations, 
    ////                maxNumberOfLinesearchIterations,successFlagBackward,minAlpha, maxAlpha, avgAlpha,xIn[0]);
    ////        //        }
    ////      } else {
    ////				overwriteMemory = false;
    ////        //        xRetrieved[0][0] = hystTMP->computeInputAndUpdate(yIn[0][0], zeroVec[0], zeroVec[0],
    ////        //						zeroVec[0], 0, eps_mu[0][0], overwriteMemory);
    ////				xRetrieved[0][0] = hystTMP->computeInputAndUpdate(yIn[0][0], eps_mu[0][0], 0, overwriteMemory,successFlagBackward);
    ////			}
    ////      
    ////      if (XML_performanceMeasurement_ == 3){
    ////        backwardTimer->Stop();
    ////      
    ////        // only count actual computations; if value gets reused, this should not count
    ////        if(successFlagBackward != 0) {
    ////          backwardEvalCounter++;
    ////				
    ////          endTime = backwardTimer->GetCPUTime();
    ////          evalTime = endTime-startTime;
    ////          if(evalTime > backwardMaxEvalTime){
    ////            backwardMaxEvalTime = evalTime;
    ////          }
    ////          backwardTotalEvalTime += evalTime;
    ////        }
    ////			}
    ////			
    ////      if(minAlpha < totalminAlpha){
    ////        totalminAlpha = minAlpha;
    ////      }
    ////      if(maxAlpha > totalmaxAlpha){
    ////        totalmaxAlpha = maxAlpha;
    ////      }
    ////      totalavgAlpha += avgAlpha;
    ////      
    ////      if(successFlagBackward == -1){
    ////        LMFails++;
    ////      } else if(successFlagBackward == 0){
    ////        totalReused++;
    ////      } else if(successFlagBackward == 1){
    ////        totalAnhystOnly++;
    ////      } else if(successFlagBackward == 2){
    ////        totalBisection++;
    ////      } else if(successFlagBackward == 3){
    ////        totalRemanence++;
    ////      } else if(successFlagBackward == 4){
    ////        totalPassedErrorTol++;
    ////      } else if(successFlagBackward == 5){
    ////        totalPassedResTolX++;
    ////      } else if(successFlagBackward == 6){
    ////        totalPassedResTolY++;
    ////      }
    ////      
    ////			totalNumberOfLMIterations += numberOfLMIterations;
    ////			totalNumberOfLinesearchIterations += numberOfLinesearchIterations;
    ////			if(maxNumberOfLinesearchIterations > absmaxNumberOfLinesearchIterations){
    ////				absmaxNumberOfLinesearchIterations = maxNumberOfLinesearchIterations;
    ////			}
    ////			
    ////      // evaluate Preisach OP to get the actual output; DO NOT OVERWRITE MEMORY
    ////      yRetrieved[0] = Vector<Double>(dim_);
    ////      yRetrieved[0].Init(0.0);
    ////      hRetrieved[0] = Vector<Double>(dim_);
    ////      hRetrieved[0].Init(0.0);
    ////      
    ////      xError[0] = Vector<Double>(dim_);
    ////      xError[0].Init(0.0);
    ////      
    ////      xError[0] = xRetrieved[0];
    ////      xError[0] -= xIn[0];
    ////      
    ////      yError[0] = Vector<Double>(dim_);
    ////      yError[0].Init(0.0);
    ////      
    ////      successFlagForward = -1;
    ////      if (XML_performanceMeasurement_ == 3) {
    ////				
    ////        startTime = forwardTimer->GetCPUTime();
    ////        forwardTimer->Start();
    ////			}
    ////      
    ////      overwriteMemory = true;
    ////      if(vector){
    ////        //				if (POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz"){
    ////        //					hRetrieved[0] = hystTMP->computeValue_vec(xIn[0], 0, overwriteMemory, overwriteDirection);
    ////        //				} else {
    ////        hRetrieved[0] = hystTMP->computeValue_vec(xRetrieved[0], 0, overwriteMemory, overwriteDirection, debugOut,successFlagForward);
    ////        //				}
    ////      } else {
    ////				//hRetrieved[0][0] = hystTMP->computeValueAndUpdate(xRetrieved[0][0], 0, 0, overwriteMemory);
    ////        hRetrieved[0][0] = hystTMP->computeValueAndUpdate(xRetrieved[0][0], 0, overwriteMemory,successFlagForward);
    ////      }
    ////      
    ////      if (XML_performanceMeasurement_ == 3){
    ////        forwardTimer->Stop();
    ////      
    ////        // only count actual computations; if value gets reused, this should not count
    ////        if(successFlagForward != 0) {
    ////          forwardEvalCounter++;
    ////				
    ////          endTime = forwardTimer->GetCPUTime();
    ////          evalTime = endTime-startTime;
    ////          if(evalTime > forwardMaxEvalTime){
    ////            forwardMaxEvalTime = evalTime;
    ////          }
    ////          forwardTotalEvalTime += evalTime;
    ////        }
    ////			}
    ////      
    ////      if(successFlagForward == -1){
    ////        forwardFails++;
    ////      } else if(successFlagForward == 0){
    ////        forwardReused++;
    ////      } else if(successFlagForward == 1){
    ////        forwardAnhystOnly++;
    ////      } else if(successFlagForward == 2){
    ////        forwardOverwrite++;
    ////      } else if(successFlagForward == 3){
    ////        forwardTMP++; 
    ////      }
    ////      
    ////      // D = eps*E + P
    ////      yRetrieved[0].Add(1.0,hRetrieved[0]);
    ////      for(UInt j = 0; j < dim_; j++){
    ////        yRetrieved[0][j] += eps_mu[j][j]*xRetrieved[0][j];
    ////      }
    ////      
    ////      yError[0] = yRetrieved[0];
    ////      yError[0] -= yIn[0];
    ////      
    ////      if(writeResultsToFiles){
    ////        //testOutput2 << std::setprecision(9) <<  "0   "<< xRetrieved[0][0] <<"    "<< xRetrieved[0][1]<<std::endl;
    ////        testOutput3 << std::setprecision(9) <<  "0   "<< xIn[0][0] <<"    "<< xIn[0][1]<<"    " << hIn[0][0] <<"    "<< hIn[0][1]<<"    " << yIn[0][0] <<"    "<< yIn[0][1]<<std::endl;
    ////        testOutput4 << std::setprecision(9) <<  "0   "<< xRetrieved[0][0] <<"    "<< xRetrieved[0][1]<<"    " << yRetrieved[0][0] <<"    "<< yRetrieved[0][1]<<std::endl;
    ////      }
    ////      //    Double inc = 2*POL_operatorParams_.inputSat_/( (Double) numTests);
    ////      
    //      
    //      
    //      UInt numFails = 0;
    //      for(UInt i = 0; i < totalSteps; i++){
    //        if( (i%10 == 0)&&(printStatistics) ){
    //          std::cout << "##### CASE " << ca << "; TEST NR " << i << "#####" << std::endl;
    //        }
    //        
    //        LOG_TRACE(coeffcthyst) << "##### CASE " << ca << "; TEST NR " << i << "#####";
    //        //      xIn[i] = Vector<Double>(dim_);
    //        //      xIn[i][0] = xIn[i-1][0] + inc;
    //        //      xIn[i][1] = xIn[i-1][0] + inc/2;
    //        overwriteMemory = false;
    //        
    //        yIn[i] = Vector<Double>(dim_);
    //        yIn[i].Init(0.0);
    //        hIn[i] = Vector<Double>(dim_);
    //        hIn[i].Init(0.0);
    //        
    //        successFlagForward = -1;
    //        if (XML_performanceMeasurement_ == 3) {
    //          startTime = forwardTimer->GetCPUTime();
    //          forwardTimer->Start();
    //          
    //        }
    //        evalTime = 0.0;
    //        Double innerTime = 0.0;
    //
    //        //std::cout << "Forward" << std::endl;
    //        if(vector){
    //          //          std::cout << "Start forward" << std::endl;
    //          hIn[i] = hystTMP->computeValue_vecMeasure(xIn[i], 0, overwriteMemory, overwriteDirection, debugOut, successFlagForward, innerTime);
    //          //          std::cout << "End forward" << std::endl;
    //        } else {  
    //					hIn[i][0] = hystTMP->computeValueAndUpdateMeasure(xIn[i][0], 0, overwriteMemory, successFlagForward, innerTime);
    //          //hIn[i][0] = hystTMP->computeValueAndUpdate(xIn[i][0], xIn[i-1][0], 0, overwriteMemory);
    //          //std::cout << ", x=" << xIn[i][0] << ", flag="<< successFlagForward << ", time=" << innerTime << std::endl;
    //        }
    //
    //        if (XML_performanceMeasurement_ == 3){
    //          forwardTimer->Stop();
    //          
    //          // only count actual computations; if value gets reused, this should not count
    //          if(successFlagForward != 0) {
    //            forwardEvalCounter++;
    //            
    //            endTime = forwardTimer->GetCPUTime();
    //            evalTime = endTime-startTime;
    //            if(evalTime > forwardMaxEvalTime){
    //              forwardMaxEvalTime = evalTime;
    //            }
    //            forwardTotalEvalTime += evalTime;
    //          }
    //        }
    //     
    ////        std::cout << "Success Flag: " << successFlagForward << std::endl;
    ////        std::cout << "evalTime: " << evalTime << std::endl;
    ////        std::cout << "innerTime: " << innerTime << std::endl;
    ////        std::cout << "forwardTotalEvalTime: " << forwardTotalEvalTime << std::endl;
    ////        
    ////        if(i == 10)
    ////          EXCEPTION("stop here");
    //        
    //        if(successFlagForward == -1){
    //          forwardFails++;
    //        } else if(successFlagForward == 0){
    //          forwardReused++;
    //        } else if(successFlagForward == 1){
    //          forwardAnhystOnly++;
    //        } else if(successFlagForward == 2){
    //          forwardOverwrite++;
    //        } else if(successFlagForward == 3){
    //          forwardTMP++; 
    //        }
    //        
    //        yIn[i].Add(1.0,hIn[i]);
    //        for(UInt j = 0; j < dim_; j++){
    //          yIn[i][j] += eps_mu[j][j]*xIn[i][j];
    //        }
    //        
    //				numberOfLMIterations = 0;
    //				numberOfLinesearchIterations = 0;
    //				maxNumberOfLinesearchIterations = 0;
    //        
    //        overwriteMemory = false;
    //        xRetrieved[i] = Vector<Double>(dim_);
    //        xRetrieved[i].Init(0.0);
    //        
    //        successFlagBackward = -1;
    //        if (XML_performanceMeasurement_ == 3) {
    //          startTime = backwardTimer->GetCPUTime();
    //          backwardTimer->Start();
    //          
    //        }
    //        
    //        if(vector){
    //          //          std::cout << "Start backward" << std::endl;          
    //          //          if (POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz"){
    //          //            xRetrieved[i] = hystTMP->computeInput_vec(yIn[i],0, eps_mu, overwriteDirection);
    //          //          } else {
    //          xRetrieved[i] = hystTMP->computeInput_vec_withStatistics(yIn[i], yRetrieved[i-1], xRetrieved[i-1], hRetrieved[i-1], 
    //                  0, eps_mu, overwriteDirection, fieldsAlignedAboveSat, hystOutputRestrictedToSat,
    //									numberOfLMIterations, numberOfLinesearchIterations, 
    //                  maxNumberOfLinesearchIterations,successFlagBackward,minAlpha, maxAlpha, avgAlpha,xIn[i]);	
    //          //          }
    //          //            std::cout << "End backward" << std::endl;
    //        } else {
    //          //          xRetrieved[i][0] = hystTMP->computeInputAndUpdate(yIn[i][0], yRetrieved[i-1][0], xRetrieved[i-1][0],
    //          //						hRetrieved[i-1][0], 0, eps_mu[0][0], overwriteMemory);		
    //					xRetrieved[i][0] = hystTMP->computeInputAndUpdate(yIn[i][0], eps_mu[0][0], 0, overwriteMemory,successFlagBackward);
    ////          successFlagBackward = 7;
    //        }
    //        // debugging
    //        //xRetrieved[i] = xIn[i];
    //        
    //        if (XML_performanceMeasurement_ == 3){
    //          backwardTimer->Stop();
    //          
    //          // only count actual computations; if value gets reused, this should not count
    //          if(successFlagBackward != 0) {
    //            backwardEvalCounter++;
    //            
    //            endTime = backwardTimer->GetCPUTime();
    //            evalTime = endTime-startTime;
    //            if(evalTime > backwardMaxEvalTime){
    //              backwardMaxEvalTime = evalTime;
    //            }
    //            backwardTotalEvalTime += evalTime;
    //          }
    //        }
    //        
    //        if(successFlagBackward == -1){
    //          LMFails++;
    //        } else if(successFlagBackward == 0){
    //          totalReused++;
    //        } else if(successFlagBackward == 1){
    //          totalAnhystOnly++;
    //        } else if(successFlagBackward == 2){
    //          totalBisection++;
    //        } else if(successFlagBackward == 3){
    //          totalRemanence++;
    //        } else if(successFlagBackward == 4){
    //          totalPassedErrorTol++;
    //        } else if(successFlagBackward == 5){
    //          totalPassedResTolX++;
    //        } else if(successFlagBackward == 6){
    //          totalPassedResTolY++;
    //        }
    //								
    //        successCodeVector[i] = successFlagBackward;
    //        
    //        if(minAlpha < totalminAlpha){
    //          totalminAlpha = minAlpha;
    //        }
    //        if(maxAlpha > totalmaxAlpha){
    //          totalmaxAlpha = maxAlpha;
    //        }
    //        totalavgAlpha += avgAlpha;
    //        
    //				totalNumberOfLMIterations += numberOfLMIterations;
    //				totalNumberOfLinesearchIterations += numberOfLinesearchIterations;
    //				if(maxNumberOfLinesearchIterations > absmaxNumberOfLinesearchIterations){
    //          //	std::cout << "maxNumberOfLinesearchIterations: " << maxNumberOfLinesearchIterations << std::endl;
    //          //	std::cout << "absmaxNumberOfLinesearchIterations: " << absmaxNumberOfLinesearchIterations << std::endl;
    //					absmaxNumberOfLinesearchIterations = maxNumberOfLinesearchIterations;
    //				}
    //        
    //        
    //        yError[i] = Vector<Double>(dim_);
    //        yError[i].Init(0.0);
    //        
    //        yRetrieved[i] = Vector<Double>(dim_);
    //        yRetrieved[i].Init(0.0);
    //        
    //        hRetrieved[i] = Vector<Double>(dim_);
    //        hRetrieved[i].Init(0.0);
    //        
    //        overwriteMemory = true;
    //        
    //        successFlagForward = -1;
    //        if (XML_performanceMeasurement_ == 3) {
    //          startTime = forwardTimer->GetCPUTime();
    //          forwardTimer->Start();
    //          
    //        }
    //
    //        if(vector){
    //          //					if (POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz"){
    //          //						hRetrieved[i] = hystTMP->computeValue_vec(xIn[i], 0, overwriteMemory, overwriteDirection);
    //          //					} else {
    //          hRetrieved[i] = hystTMP->computeValue_vec(xRetrieved[i], 0, overwriteMemory, overwriteDirection, debugOut, successFlagForward);
    //          
    //          //					}
    //          //hRetrieved[i] = hystTMP->computeValue_vec(xRetrieved[i], 0, overwriteMemory, overwriteDirection);
    //        } else {
    //          hRetrieved[i][0] = hystTMP->computeValueAndUpdate(xRetrieved[i][0], 0, overwriteMemory, successFlagForward);
    //          
    //          //hRetrieved[i][0] = hystTMP->computeValueAndUpdate(xRetrieved[i][0], xRetrieved[i-1][0], 0, overwriteMemory);
    //        }
    //
    //        if (XML_performanceMeasurement_ == 3){
    //          forwardTimer->Stop();
    //          
    //          // only count actual computations; if value gets reused, this should not count
    //          if(successFlagForward != 0) {
    //            forwardEvalCounter++;
    //            
    //            endTime = forwardTimer->GetCPUTime();
    //            evalTime = endTime-startTime;
    //            if(evalTime > forwardMaxEvalTime){
    //              forwardMaxEvalTime = evalTime;
    //            }
    //            forwardTotalEvalTime += evalTime;
    //          }
    //        }
    //        
    //        if(successFlagForward == -1){
    //          forwardFails++;
    //        } else if(successFlagForward == 0){
    //          forwardReused++;
    //        } else if(successFlagForward == 1){
    //          forwardAnhystOnly++;
    //        } else if(successFlagForward == 2){
    //          forwardOverwrite++;
    //        } else if(successFlagForward == 3){
    //          forwardTMP++; 
    //        }
    //        
    //        eps_mu.Mult(xRetrieved[i],yRetrieved[i]);
    //        yRetrieved[i].Add(hRetrieved[i]);
    //        
    //        yError[i] = yRetrieved[i];
    //        yError[i].Add(-1.0,yIn[i]);
    //        
    //        xError[i] = Vector<Double>(dim_);
    //        xError[i].Init(0.0);
    //        
    //        xError[i] = xRetrieved[i];
    //        xError[i] -= xIn[i];
    //        
    //        if(writeResultsToFiles){
    //          //testOutput2 << i << std::setprecision(9) << "   "<< xRetrieved[i][0] <<"    "<< xRetrieved[i][1]<<std::endl;
    //          testOutput3 << i << std::setprecision(9) << "   "<< xIn[i][0] <<"    "<< xIn[i][1]<< "   "<< hIn[i][0] <<"    "<< hIn[i][1]<< "   "<< yIn[i][0] <<"    "<< yIn[i][1]<<std::endl;
    //          testOutput4 << i << std::setprecision(9) << "   "<< xRetrieved[i][0] <<"    "<< xRetrieved[i][1]<< "   "<< hRetrieved[i][0] <<"    "<< hRetrieved[i][1]<< "   "<< yRetrieved[i][0] <<"    "<< yRetrieved[i][1]<<std::endl;
    //        }
    //        
    //        LOG_TRACE(coeffcthyst) << "##### TARGET X-VECTOR #####";
    //        LOG_TRACE(coeffcthyst) << xIn[i].ToString();
    //        LOG_TRACE(coeffcthyst) << "##### RETRIEVED X-VECTOR #####";
    //        LOG_TRACE(coeffcthyst) << xRetrieved[i].ToString();
    //        
    //        LOG_TRACE(coeffcthyst) << "##### ERROR VECTOR wrt X #####";
    //        LOG_TRACE(coeffcthyst) << xError[i].ToString();
    //        
    //        LOG_TRACE(coeffcthyst) << "##### HYST_OP EVALUATED WITH TARGET X-VECTOR #####";
    //        LOG_TRACE(coeffcthyst) << hIn[i].ToString();
    //        
    //        LOG_TRACE(coeffcthyst) << "##### HYST_OP EVALUATED WITH RETRIEVED X-VECTOR #####";
    //        LOG_TRACE(coeffcthyst) << hRetrieved[i].ToString();
    //        
    //        LOG_TRACE(coeffcthyst) << "##### TARGET Y-VECTOR #####";
    //        LOG_TRACE(coeffcthyst) << yIn[i].ToString();
    //        LOG_TRACE(coeffcthyst) << "##### RETRIEVED Y-VECTOR #####";
    //        LOG_TRACE(coeffcthyst) << yRetrieved[i].ToString();
    //        
    //        LOG_TRACE(coeffcthyst) << "##### ERROR VECTOR wrt Y #####";
    //        LOG_TRACE(coeffcthyst) << yError[i].ToString();
    //        
    //        if(yError[i].NormL2() > LM_inversion_.tolB){
    //          numFails++;
    //          failedTests[i] = yError[i].NormL2();
    //        } 
    //      }
    //      
    //      UInt LMcases = totalSteps-totalReused-totalBisection-totalAnhystOnly-totalRemanence;
    //      if(LMcases != 0){
    //        totalavgAlpha /= (Double) LMcases;
    //      }
    //      
    //      if(writeResultsToFiles){
    //        //testOutput2.close();
    //        testOutput3.close();
    //        testOutput4.close();
    //      }
    //      
    //      std::stringstream statistics;
    //      statistics << std::endl;	
    //      statistics << "############################# " <<	std::endl;	
    //      statistics << "##### RESULTS FOR TEST CASE " << ca <<	std::endl;	
    //      statistics << " " << totalSteps-numFails << " of " << totalSteps << " satisfied at least the failback criterion, i.e. |residual_Y| = |Y - mu_eps*X - P(X)| < " << LM_inversion_.tolB << std::endl;
    //      statistics << " " << numFails << " of " << totalSteps << " failed to satisfy the failback criterion" << std::endl;
    //      if(numFails != 0){
    //        statistics << "## Detailed Statistics for failed tests: " << std::endl;
    //        std::map<UInt, Double>::iterator mapIt;
    //        for(mapIt = failedTests.begin(); mapIt != failedTests.end(); mapIt++){
    //          statistics << " Test Nr " << mapIt->first << std::endl;
    //          statistics << " - Remaining |residual_Y|: " << mapIt->second << std::endl;
    //          statistics << " - Y_in: " << yIn[mapIt->first].ToString() << std::endl;
    //          statistics << " - Y_ret: " << yRetrieved[mapIt->first].ToString() << std::endl;
    //          statistics << " - h_sol: " << hIn[mapIt->first].ToString() << std::endl;
    //          statistics << " - h_ret: " << hRetrieved[mapIt->first].ToString() << std::endl;
    //          bool xsolabovesat = xIn[mapIt->first].NormL2()>POL_operatorParams_.inputSat_;
    //          bool xretabovesat = xRetrieved[mapIt->first].NormL2()>POL_operatorParams_.inputSat_;
    //          statistics << " - X_sol: " << xIn[mapIt->first].ToString()  << "; above sat? "<< xsolabovesat << std::endl;
    //          statistics << " - X_ret: " << xRetrieved[mapIt->first].ToString() << "; above sat? "<< xretabovesat << std::endl;
    //          statistics << " - SuccessCode: " << successCodeVector[mapIt->first] << std::endl;
    //        }
    //      }
    //      
    //      statistics << "## Detailed Statistics on forward evaluations (excluding those inside of LM/during inversion): " << std::endl;
    //      statistics << " " << forwardReused << " of " << 2*totalSteps << " reused old output as input change was below tolerance" << std::endl;
    //      statistics << " " << forwardAnhystOnly << " of " << 2*totalSteps << " featured only anhysteretic parts" << std::endl;
    //      statistics << " " << forwardOverwrite << " of " << 2*totalSteps << " evaluations were performed on permanent storage" << std::endl;
    //      statistics << " " << forwardTMP << " of " << 2*totalSteps << " evaluations were performed on temporal storage" << std::endl;
    //      
    //      statistics << "## Detailed Statistics on inversion: " << std::endl;
    //      statistics << " " << totalReused << " of " << totalSteps << " reused old solution as DeltaY < " << LM_inversion_.tolB << std::endl;
    //      statistics << " " << totalAnhystOnly << " of " << totalSteps << " were solved using bisection for pure anhysteretic case" << std::endl;
    //      statistics << " " << totalBisection << " of " << totalSteps << " were solved using bisection / simple division" << std::endl;
    //      if(!vector){
    //        statistics << " " << totalPassedErrorTol << " of " << totalSteps << " tests passed due to error tolerance" << std::endl;
    //      }
    //      
    //      if(vector){
    //        statistics << " " << totalRemanence << " of " << totalSteps << " cases were in remanence" << std::endl;
    //        statistics << " " << LMcases << " of " << totalSteps << " were solved using Levenberg-Marquardt" << std::endl;
    //        if(LMcases != 0){
    //          statistics << "## Detailed Statistics Levenberg-Marquardt Statistics: " << std::endl;
    //          statistics << " " << totalPassedErrorTol << " of " << LMcases << " tests passed due to error tolerance |JacT*Res| < " << LM_inversion_.tolH << std::endl;
    //          statistics << " " << totalPassedResTolX << " of " << LMcases << " tests passed due to |residual_X| = |X - mu_eps^-1*(Y - P(X)| < " << LM_inversion_.tolH << std::endl;
    //          statistics << " " << totalPassedResTolY << " of " << LMcases << " tests passed due to |residual_Y| = |Y - mu_eps*X - P(X)| < " << LM_inversion_.tolH << std::endl;
    //          statistics << " " << LMFails << " of " << LMcases << " failed the failback criterion (|residual_Y| < " << LM_inversion_.tolB << ")" << std::endl;
    //          statistics << " Total number of LM iterations: " << totalNumberOfLMIterations << std::endl;
    //          statistics << " Average number of LM iterations: " << (Double) totalNumberOfLMIterations/LMcases << std::endl;
    //          statistics << " Total number of Linesearch iterations: " << totalNumberOfLinesearchIterations << std::endl;
    //          statistics << " Average number of Linesearch iterations (per LM Iteration): " << (Double) totalNumberOfLinesearchIterations/totalNumberOfLMIterations << std::endl;
    //          statistics << " Maximal number of Linesearch iterations: " << absmaxNumberOfLinesearchIterations << std::endl;
    //          statistics << " Minimal alphaLS: " << totalminAlpha << std::endl;
    //          statistics << " Maximal alphaLS: " << totalmaxAlpha << std::endl;
    //          statistics << " Average alphaLS: " << totalavgAlpha << std::endl;
    //        }
    //      }
    //      
    //      if (XML_performanceMeasurement_ == 3) {
    //        
    //        backwardAvgEvalTime = backwardTotalEvalTime/backwardEvalCounter;
    //        forwardAvgEvalTime = forwardTotalEvalTime/forwardEvalCounter;
    //        statistics << "## Runtime information: " << std::endl;
    //        statistics << "Total number of forward evaluations (excluding those inside LM): " << forwardEvalCounter << std::endl;
    //        statistics << "Total time for forward evaluations (excluding those inside LM): " << forwardTotalEvalTime << std::endl;
    //        statistics << "Average time for forward evaluations (excluding those inside LM): " << forwardAvgEvalTime << std::endl;
    //        statistics << "Maximal time for forward evaluations (excluding those inside LM): " << forwardMaxEvalTime << std::endl;
    //        statistics << "Total number of backward evaluations: " << backwardEvalCounter << std::endl;
    //        statistics << "Total time for backward evaluations: " << backwardTotalEvalTime << std::endl;
    //        statistics << "Average time for backward evaluations: " << backwardAvgEvalTime << std::endl;
    //        statistics << "Maximal time for backward evaluations: " << backwardMaxEvalTime << std::endl;
    //        //statistics << "TOGREP: " << forwardAvgEvalTime << " (" << forwardEvalCounter << ") " << forwardTotalEvalTime << " " << backwardAvgEvalTime << " (" << backwardEvalCounter << ") " << " " << backwardTotalEvalTime << " " << numFails << std::endl;
    //        statistics << "TOGREP: \t\t" << forwardAvgEvalTime << " (" << forwardEvalCounter << ") " << "\t " << backwardAvgEvalTime << " (" << backwardEvalCounter << ") " << "\t " << numFails << std::endl;
    //			}
    //      statistics << "############################# " <<	std::endl;	
    //      statistics << std::endl;	
    //      if(printStatistics){
    //        std::cout << statistics.str();
    //			}
    //      if(writeResultsToFiles){
    //        testOutput5.open(filenameLog);
    //        testOutput5 << statistics.str();
    //        testOutput5.close();
    //      }
    //			
    //      LOG_TRACE(coeffcthyst) << statistics.str();
    //      
    //    }
    //    if(abortAfterTests){
    //      EXCEPTION("STOP AFTER TESTS");
    //    }
  }
  
  
  std::string CoefFunctionHyst::ToString() const {
    
    static UInt initialOutput = 0;
    std::ostringstream oss;
    
    if (initialOutput == 0) {
      oss << "+++ CoefFunctionHyst +++\n";
      oss << "++ General information from FE context: \n"
              << "  coefFunctionHyst created for pde " << PDEName_ << "\n"
              << "  coefFunctionHyst depends on " << dependCoef1_->ToString() << "\n";
      oss << "\n";
      
      oss << "++ Parameter extracted from input: \n"
              << "+ From material file: \n"
              << "  initial/smallfield tensor " << POL_eps_mu_SmallSignal_.ToString() << "\n"
              << "  xSaturation " << POL_operatorParams_.inputSat_ << "\n"
              << "  ySaturation " << POL_operatorParams_.outputSat_ << "\n";
      
      if (POL_operatorParams_.methodType_ == 0) {
        oss << "  dirP " << POL_operatorParams_.fixDirection_.ToString() << "\n";
      } else {
        oss << "  rotRes " << POL_operatorParams_.rotResistance_ << "\n"
                << "  angResistance " << POL_operatorParams_.angularDistance_ << "\n"
                << "  angClipping " << POL_operatorParams_.angularClipping_ << "\n"
                << "  angResolution " << POL_operatorParams_.angularResolution_ << "\n"
                << "  amplitudeResolution " << POL_operatorParams_.amplitudeResolution_ << "\n"
                << "  evalVersion " << POL_operatorParams_.evalVersion_ << "\n";
        oss << "  initial input to hyst operator " << POL_initialInput_.ToString() << "\n";
        //        oss << "  initial output of hyst operator (exemplary for element 0) " << POL_initialOutput_[0].ToString() << "\n";
        oss << "  printOut " << POL_operatorParams_.printOut_ << "\n";
        if (POL_operatorParams_.printOut_) {
          oss << "  bmpResolution " << POL_operatorParams_.bmpResolution_ << "\n";
        }
      }
      oss << "\n";
      
      oss << "+ From xml file: \n"
              << "  evaluationDepth " << XML_EvaluationDepth_ << "\n"
              << "  performanceMeasurement " << XML_performanceMeasurement_ << "\n"
              << "\n";
    }
    
    if (initialOutput > 0) {
      oss << "++ Current state of runtime parameter: \n"
              << "  evaluationPurpose " << RUN_evaluationPurpose_ << "\n"
              << "  allowBMP " << RUN_allowBMP_ << "\n";
      
      if (XML_performanceMeasurement_ == 1) {
        oss << "++ Performance measurements:\n "
                << "  Total number of calls to EvaluateHysteresisOperator: " << totalCallingCounter_ << "\n"
                << "  Total number of actual evaluations: " << totalEvaluationCounter_ << "\n"
                << "  Total evaluation time for Hysteresis operator: " << totalEvaluationTime_ << "\n"
                << "  Average evaluation time for Hysteresis operator: " << avgEvaluationTime_ << "\n";
        if (XML_performanceMeasurement_ == 2) {
          oss << hyst_->runtimeToString();
        }
      }
    }
    initialOutput++;
    
    return oss.str();
  }
  
}// namespace
