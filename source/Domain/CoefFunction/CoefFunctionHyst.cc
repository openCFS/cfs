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

#include "DataInOut/SimOutput.hh"

namespace CoupledField {
    /*
   * Logging system
   *
   * LOG_TRACE = Calls to / steps inside main functions / most important results
   * LOG_TRACE2 = Steps inside sub functions
   * LOG_DBG = Global debugging info like convergence results
   * LOG_DBG2 = Specialized debugging stuff like state of vectors
   * LOG_DBG3 = Very specialized debugging stuff e.g. for checking input parameter
   */
	DEFINE_LOG(coeffunctionhyst_main, "coeffunctionhyst_main")

  // important subfunctions get own logger
	DEFINE_LOG(coeffunctionhyst_helper, "coeffunctionhyst_helper")
	DEFINE_LOG(coeffunctionhyst_deltamat, "coeffunctionhyst_deltamat")

  // forward definition for static members
  std::map< std::string, TracedData > CoefFunctionHyst::tracedOperatorData_;
  std::mutex CoefFunctionHyst::tracedDataMutex_;

  void CoefFunctionHyst::ComputeVector(Vector<Double>& outputVector,const LocPointMapped& lpm, int timeLevel, int baseSign,
          std::string vectorName, bool onBoundary, bool usedAsRHSload ){
    LOG_DBG(coeffunctionhyst_helper) << "+++++ Coef Function Hyst Vector - Get " << vectorName <<" ++++++";
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
      outputVector = GetPrecomputedOutputOfHysteresisOperator(lpm,timeLevel, false);

      if((GetFPMaterialStateRHS() != 0)&&(usedAsRHSload)){
        Vector<Double> correctionVector = GetFPCorrectionVector(lpm, timeLevel);
        outputVector.Add(-1.0,correctionVector);
      }
      //std::cout << "BaseSign = " << baseSign << std::endl;

      outputVector.ScalarMult(specificSign*baseSign);

    } else if ((vectorName == "IrrStressForMechPDE") || (vectorName == "IrrStressesPiezo_VectorForm")
            || (vectorName == "IrrStressesMagstrict_VectorForm")){
//      std::cout << vectorName << " was requested" << std::endl;
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

      Vector<Double> S_irr = GetIrreversibleStrains(lpm, timeLevel);
      // c is not solution dependet here > take it directly
      outputVector = Vector<Double>(S_irr.GetSize());
      outputVector.Init();
      LOG_DBG(coeffunctionhyst_helper) << "S_irr " << S_irr.ToString();
      LOG_DBG(coeffunctionhyst_helper) << "elastTensor_: " << elastTensor_.ToString();
      elastTensor_.Mult(S_irr,outputVector);
      LOG_DBG(coeffunctionhyst_helper) << "elastTensor_: " << elastTensor_.ToString();
      LOG_DBG(coeffunctionhyst_helper) << "elastTensor_.Mult(S_irr,outputVector): " << outputVector.ToString();
      outputVector.ScalarMult(specificSign*baseSign);
      //      std::cout << "outputVector: " << outputVector.ToString() << std::endl;
    } else if ((vectorName == "CoupledIrrStrainForElecPDE")||(vectorName == "CoupledIrrStrainsPiezo")){
//            std::cout << vectorName << " was requested" << std::endl;
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

      Vector<Double> S_irr = GetIrreversibleStrains(lpm, timeLevel);
      // c is not solution dependet here > take it directly
      outputVector = Vector<Double>(GetVecSize());
      outputVector.Init();


//      std::cout << "e(P) = " << std::endl;
//      std::cout << couplTensor.ToString() << std::endl;
//      std::cout << "S_irr = " << std::endl;
//      std::cout << S_irr.ToString() << std::endl;
      couplTensor.Mult(S_irr,outputVector);

//      std::cout << "specificSign = " << specificSign << std::endl;
//      std::cout << "baseSign = " << baseSign << std::endl;
      outputVector.ScalarMult(specificSign*baseSign);

    } else if ((vectorName == "CoupledIrrStrainForMagPDE")||(vectorName == "CoupledIrrStrainsMagstrict")){
//            std::cout << vectorName << " was requested" << std::endl;
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

      Vector<Double> S_irr = GetIrreversibleStrains(lpm, timeLevel);
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
      outputVector = GetPrecomputedOutputOfHysteresisOperator(lpm,timeLevel,false);

//      if((GetFPMaterialStateRHS() != 0)&&(usedAsRHSload)){
//        WARN("Mag polarization should not be on rhs; use magnetization instead");
//      }
      //std::cout << "BaseSign = " << baseSign << std::endl;

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
      outputVector = GetPrecomputedInputToHysteresisOperator(lpm,timeLevel);

    } else if(vectorName == "MagMagnetization"){
//      std::cout << "MagMagnetization was requested" << std::endl;
      // rhs (magnetics, w = testfunction):
      // + int_Volume rot(w)^T M dOmega
      // > specificSign = +
      specificSign = 1.0;

      outputVector = Vector<Double>(GetVecSize());
      outputVector.Init();
      // get polarization J_pol
      outputVector = GetPrecomputedOutputOfHysteresisOperator(lpm,timeLevel,false);
//      std::cout << "Output of hyst operator: " << outputVector.ToString() << std::endl;

      // get the matrix that was used to compute P from B to scale P to M
      Matrix<Double> mu = Matrix<Double>(GetVecSize(),GetVecSize());
      Matrix<Double> nu = Matrix<Double>(GetVecSize(),GetVecSize());
      getMatrixForLocalInversion(lpm, mu, nu);
      LOG_DBG(coeffunctionhyst_helper) << "Used material tensor during inversion (mu):" << mu.ToString();
      LOG_DBG(coeffunctionhyst_helper) << "Used material tensor for computation of magnetization (nu):" << nu.ToString();

      //std::cout << "BaseSign = " << baseSign << std::endl;

      // scale with the actual nu
      Vector<Double> tmpVec = Vector<Double>(GetVecSize());
      nu.Mult(outputVector,tmpVec);
//      std::cout << "tmpVec: " << tmpVec.ToString() << std::endl;
      outputVector = tmpVec;

      // correction is already in terms of H, so no further scaling with nu!
      if((GetFPMaterialStateRHS() != 0)&&(usedAsRHSload)){
        Vector<Double> correctionVector = GetFPCorrectionVector(lpm, timeLevel);
        outputVector.Add(-1.0,correctionVector);
      }
//
//      std::cout << "Return: " << outputVector.ToString() << std::endl;
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

  void CoefFunctionHyst::ComputeTensor(Matrix<Double>& outputTensor, const LocPointMapped& lpm,
          std::string tensorName, std::string implementationVersion, bool transposed, bool rotate, bool useAbs, bool lockPrecomputationAndDeltaMat ){
//    LOG_DBG(coeffunctionhyst_helper) << "+++++ Coef Function Hyst Mat - Get " << tensorName <<" ++++++";

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
    GetTensorSize(numRows,numCols);//,tensorName);
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
    //        LOG_DBG(coeffunctionhyst_helper) << "Initialize linear tensors";
    //        InitLinearTensors(lpm);
    //      }
    //      std::cout << "+++++ Coef Function Hyst Mat - Compute Tensor 2 ++++++" << std::endl;
    //      if(hystCoefFunction_->anyMatrixForLocalInversionRequiresComputation() &&
    //              ((tensorName == "Reluctivity")||(tensorName == "CouplingMechToMag")||(tensorName == "CouplingMagToMech")) ){
    //        LOG_DBG(coeffunctionhyst_helper) << "Set matrix for inversion";
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
    //          LOG_DBG(coeffunctionhyst_helper) << "PrecomputeMaterialTensorForInverison unlocked";
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

    int deltaForm_ = GetDeltaForm();
    bool deltaFormActive_ = deltaMatActive();
    if( (implementationVersion == "none") || (lockPrecomputationAndDeltaMat == true) ){
      deltaFormActive_ = false;
    }
    //      std::cout << "+++++ Coef Function Hyst Mat - Compute Tensor 4 ++++++" << std::endl;
    // always compute deltaMat from current value to a previous one (timeLevel_to_diff)
    // except for the precomputation case (which calls this function with the flag
    // lockPrecomputationAndDeltaMat = true)
    int timelevel_cur;
//    if(lockPrecomputationAndDeltaMat == true){
//      // no longer required: in old version, we computed the rotated/scaled coupling tensors here which require
//      // the polarization/output of hyst operator which itself might require the material tensors and so on;
//      // to avoid troubles, we simply used the old polarization at this point to avoid trouble
//      //
//      // currently however we precompute polarization first and then the coupling tensors before coming to this function
//      // in that sense, we do not even use the timelevel_cur flag anymore for getting the coupling tensor (as we just
//      // take the last precomputed state)
//      // > timelevel_cur not relevant here
//      timelevel_cur = 1; // value from previous iteration
//    } else {
      timelevel_cur = 0; // current value
//    }
    //    int timeLevel_to_diff = GetTimeLevel("DeltaMat");

      /*
       * -- NOTE 1.7.2019 --
       * in case of non-isotropic materials, like transversalisotropic piezos, we do not only have
       * to rotate the coupling tensor towards the polarization direction but also the other material
       * tensors!
       * per default, the tensors from the material file are assumed for a polarization along the z-axis;
       * in 2d those tensors are automatically rotated such that z -> y; however if the polarization changes
       * during runtime, we have to adapt to the new directions!
       *
       * TODO:
       * 1) add check for anisotropic material tensors (e.g., eps and c in mechanics)
       * 2) if tensors are isotropic > just take them; otherwise rotate them accordingly
       *
       */


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
      Vector<Double> Si_vector = GetIrreversibleStrains(lpm, 0);
      Vector<Double> tmpVector = Vector<Double>(Si_vector.GetSize());
      elastTensor_.Mult(Si_vector,tmpVector);
      tmp = ConvertFromVoigtToTensor(tmpVector);

    } else if((tensorName == "IrrStrainsPiezo_TensorForm")||(tensorName == "IrrStrainsMagstrict_TensorForm")){
      // always compute the current timelevel > 0
      // > tensor form is only used for output; otherwise use vector form
      tmp = GetIrreversibleStrainTensor(lpm, 0);
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
      // in case of e-form, fieldTensor is seen (and set to be) epsS which needs
      // no further treatment
      Matrix<Double> e_scaled;
      Matrix<Double> tmp2 = Matrix<Double>(numRows,numCols);
      if(strainForm_ == 1){
//        std::cout << "Get Permittivity - Compute epsS from epsT" << std::endl;
        // d-form shall be used as basis; fieldTensor such represents epsT
        // to get the required epsS we have to compute
        //  epsS = epsT - d(P)[c]d(P)^T
        // where d(P) is the scaled and rotated coupling tensor in d-form
        Matrix<Double> d_scaled, d_scaled_transposed;

        GetScaledAndRotatedCouplingTensor(lpm,couplTensor_,d_scaled,timelevel_cur,rotate);

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
        if(GetFPMaterialState() == 0){
          tmp = eps_nu_base_;
        } else {
          tmp = GetFPMaterialTensor(lpm);
        }

        // tmp = epsT - d*c*d^T
        tmp.Add(-1.0,tmp2);

      } else {
//        std::cout << "Get Permittivity - Take epsS directly" << std::endl;
        // uncoupled case or e-form
        // just take epsS directly
        if(GetFPMaterialState() == 0){
          tmp = eps_nu_base_;
        } else {
//          std::cout << "Get Permittivity - GetFPMaterialTensor" << std::endl;
          tmp = GetFPMaterialTensor(lpm);
        }
//        std::cout << "Obtained tensor: " << tmp.ToString() << std::endl;
        // also get e in case we use deltaMat for strains, too
        GetScaledAndRotatedCouplingTensor(lpm,couplTensor_,e_scaled,timelevel_cur,rotate);
      }
      // check if deltaMat shall be added
      // here we have two flags:
      //  deltaForm_ is used to indicate if we are using a delta formulation in general
      //  deltaFormActive_ indicates if we actually need this deltaMatrix for the current
      //    evaluation (the issue is, that we need a deltaMatrix on the lhs and a non-deltaMatrix
      //    on the rhs
      //      if(deltaFormActive_ && (deltaForm_ != 0) && (lockPrecomputationAndDeltaMat == false) ) {
      if(deltaMat_Pol_active){
//         std::cout << "Get Permittivity - Compute DeltaMatrix" << std::endl;
        //        std::cout << "Compute DeltaMatrix" << std::endl;
        bool useStrains = false;

        // deltaMat will be computed using the current value (timelevel_cur = 0)
        // and the value at timelevel (-1 > last ts; +1 > last iteration)
        Matrix<Double> deltaMat = GetDeltaMat(lpm, timelevel_cur, timeLevelDeltaMat_Pol, useStrains, useAbs,implementationVersion);
        //        std::cout << "deltaMat elec = " << deltaMat.ToString() << std::endl;
        tmp.Add(1.0,deltaMat);
//        std::cout << "DeltaMat: " << deltaMat.ToString() << std::endl;
        if( (strainForm_ != -2) && (deltaMat_Strain_active) ){
//          std::cout << "Get Permittivity - Add dS/dE" << std::endl;
          // coupled case
          // here we have to add -e*dS/dE in addition
          useStrains = true;
          Matrix<Double> deltaMat_strains = GetDeltaMat(lpm, timelevel_cur, timeLevelDeltaMat_Strain, useStrains, useAbs,implementationVersion);
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
        rotatedCouplTensor = couplTensor_;
        rotatedCouplTensor.Init();
        tmp = rotatedCouplTensor;
      } else if(strainForm_ == 1){
        // use d-form/g-form as basis, i.e. the followings steps have to be applied
        // 1. scale and rotate d/g
        GetScaledAndRotatedCouplingTensor(lpm,couplTensor_,rotatedCouplTensor,timelevel_cur,rotate);
        // 2. compute e/h from d/g
        // e = d*c ; h = g*c
        rotatedCouplTensor.Mult(elastTensor_,tmp);
      } else {
        // use e-form/h-form as basis, i.e. the following step has to be done
        // 1. acale and rotate e
        GetScaledAndRotatedCouplingTensor(lpm,couplTensor_,rotatedCouplTensor,timelevel_cur,rotate);
        tmp = rotatedCouplTensor;
      }

    } else if ((tensorName == "CouplingElecToMech")||(tensorName == "CouplingMagToMech")){
      // this is basically the same tensor as mechToElec except for the case of
      // deltaFormulation; in the later case, we have to add c*dS/dE
      //      std::cout << "ComputeElecToMech" << std::endl;
      Matrix<Double> rotatedCouplTensor;

      if(strainForm_ == -1){
        // small signal tensor defined, but shall not be used > return empty matrix/vector
        rotatedCouplTensor = couplTensor_;
        rotatedCouplTensor.Init();
        tmp = rotatedCouplTensor;
      } else if(strainForm_ == 1){
        // use d-form as basis, i.e. the followings steps have to be applied
        // 1. scale and rotate d/g
        GetScaledAndRotatedCouplingTensor(lpm,couplTensor_,rotatedCouplTensor,timelevel_cur,rotate);
        // 2. compute e/h from d/g
        // e = d*c ; h = g*c
        rotatedCouplTensor.Mult(elastTensor_,tmp);
      } else {
        // use e-form/h-form as basis, i.e. the following step has to be done
        // 1. acale and rotate e/h
        GetScaledAndRotatedCouplingTensor(lpm,couplTensor_,rotatedCouplTensor,timelevel_cur,rotate);
        tmp = rotatedCouplTensor;
      }

      //      if(deltaFormActive_ && (deltaForm_ != 0) && (lockPrecomputationAndDeltaMat == false) ) {
      if( (deltaMat_Coupling_active)&&(strainForm_ != -2) ){
        LOG_DBG(coeffunctionhyst_helper) << "Compute DeltaMatrix";
        //        std::cout << "Compute DeltaMatrix" << std::endl;
        bool useStrains = true;

        // deltaMat will be computed using the current value (timelevel_cur = 0)
        // and the value at timelevel (-1 > last ts; +1 > last iteration)
        Matrix<Double> deltaMat = GetDeltaMat(lpm, timelevel_cur, timeLevelDeltaMat_Coupling, useStrains, useAbs,implementationVersion);
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
    } else if( (tensorName == "Reluctivity") || (tensorName == "LinearReluctivity") ){
//      std::cout << "Get Reluctivity" << std::endl;
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

        GetScaledAndRotatedCouplingTensor(lpm,couplTensor_,g_scaled,timelevel_cur,rotate);

        // h = g*c
        g_scaled.Mult(elastTensor_,h_scaled);
        h_scaled.Transpose(g_scaled_transposed);
        // tmp2 = h*g^T = gcg^T
        h_scaled.Mult(g_scaled_transposed,tmp2);

        // tmp = nuT
        if(GetFPMaterialState() == 0){
          tmp = eps_nu_base_;
        } else {
          tmp = GetFPMaterialTensor(lpm);
        }

//        tmp = eps_nu_base_;

        // tmp = nuT + g*c*g^T
        tmp.Add(1.0,tmp2);

      } else {
        //std::cout << "Get Reluctivity - Take nuS directly" << std::endl;
        // uncoupled case or h-form
        // just take nuS directly
        if(GetFPMaterialState() == 0){
          tmp = eps_nu_base_;
        } else {
          tmp = GetFPMaterialTensor(lpm);
        }
      }

      // check if deltaMat shall be added
      // here we have two flags:
      //  deltaForm_ is used to indicate if we are using a delta formulation in general
      //  deltaFormActive_ indicates if we actually need this deltaMatrix for the current
      //    evaluation (the issue is, that we need a deltaMatrix on the lhs and a non-deltaMatrix
      //    on the rhs
      if( (deltaMat_Pol_active)&&(tensorName != "LinearReluctivity") ){

        bool useStrains = false;
        bool oldVersion = false;
        bool newVersion = !oldVersion;

        // for testing > compute both but return only the new one
        bool testing = !true;
        bool miniOutput = false;
        UInt operatorIdx, storageIdx,idxElem;
        Matrix<Double> tmp_BAK;
        if(testing){
          LocPointMapped actualLPM;
          PreprocessLPM(lpm, actualLPM, operatorIdx, storageIdx,idxElem);
          if(operatorIdx == 1){
            oldVersion = true;
            newVersion = true;
            miniOutput = true;

            if(miniOutput){
              std::cout << "+++ DeltaMat/Jacobian Reluctivity for operator "<<operatorIdx<<" +++" << std::endl;
              std::cout << "timelevel_cur = " << timelevel_cur << std::endl;
              std::cout << "timeLevelDeltaMat_Pol = " << timeLevelDeltaMat_Pol << std::endl;
            }
          }
        }

        if(oldVersion && newVersion){
          tmp_BAK = Matrix<Double>(numRows,numCols);
          tmp_BAK = tmp;
        }

        if(oldVersion){
          // deltaMat will be computed using the current value (timelevel_cur = 0)
          // and the value at timelevel (-1 > last ts; +1 > last iteration)
          // note: here we do not have to pass a flag forceMidpoint as GetDeltaMat checks internally
          Matrix<Double> deltaMat = GetDeltaMat(lpm, timelevel_cur, timeLevelDeltaMat_Pol, useStrains, useAbs,implementationVersion);

          // note: deltaMat will compute dP/dB, but we want to have -dM/dB
          // i.e. add -nu*deltaMat to tmp
          // in other words we want to have:
          // dH = nu*(Identiy - deltaMat)*dB
          Matrix<Double> tmp3 = Matrix<Double>(numRows,numCols);
          tmp.Mult(deltaMat,tmp3);
          tmp.Add(-1.0,tmp3);

          if(miniOutput){
            std::cout << "OLD IMPLEMENTATION" << std::endl;
            std::cout << tmp.ToString() << std::endl;
          }
        }

        if(oldVersion && newVersion){
          tmp = tmp_BAK;
          deltaMat_requiresReeval_[storageIdx] = true;
        }

        if(newVersion){
          Matrix<Double> tmpInverted = Matrix<Double>(dim_,dim_);
          tmp.Invert(tmpInverted);
          // note: here we do not have to pass a flag forceMidpoint as GetMaterialRelation checks internally
          tmp = GetMaterialRelation(lpm, timelevel_cur, timeLevelDeltaMat_Pol, tmpInverted, "Reluctivity", false);

          if(miniOutput){
            std::cout << "NEW IMPLEMENTATION" << std::endl;
            std::cout << tmp.ToString() << std::endl;
          }
        }

        //tmp.Add(-795774.7155,deltaMat);
        //std::cout << "DeltaMat (dP/dB): " << deltaMat.ToString() << std::endl;
        if( (strainForm_ != -2) && (deltaMat_Strain_active) ){
          //        if(strainForm_ != -1){
          //std::cout << "Get Reluctivity - Add dS/dB" << std::endl;
          // coupled case
          // here we have to add +g[c] dS/dB in addition
          useStrains = true;
          // note: here we do not have to pass a flag forceMidpoint as GetDeltaMat checks internally
          Matrix<Double> deltaMat_strains = GetDeltaMat(lpm, timelevel_cur, timeLevelDeltaMat_Strain, useStrains, useAbs,implementationVersion);
          h_scaled.Mult(deltaMat_strains,tmp2);

          tmp.Add(1.0,tmp2);
          //std::cout << "DeltaMatStrain: " << deltaMat_strains.ToString() << std::endl;
        }

      } else if (deltaForm_ != 0){
        //        std::cout << "deltaFormSet_ == true, but deltaMatNotActive!" << std::endl;
      } else {
        //        std::cout << "deltaFormSet_ == false" << std::endl;
      }
//      std::cout << tmp.ToString() << std::endl;
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

    LOG_DBG(coeffunctionhyst_helper) << "Computed material tensor:" << outputTensor.ToString();
    //      std::cout << "+++++ Coef Function Hyst Mat - Compute Tensor END ++++++" << std::endl;
    //std::cout << "Return the following tensor: " << outputTensor.ToString() << std::endl;
//        std::cout << "Computed tensor - " << tensorName << " - " << std::endl;
//        std::cout << outputTensor.ToString() << std::endl;
  }

  void CoefFunctionHyst::PrecomputeMaterialTensorForInverison(){
    LOG_DBG(coeffunctionhyst_helper) << "CoefFunctionHelper::PrecomputeMaterialTensorForInverison()";

    /*
     * Iterate over all integration points used/managed by the hyst operator and
     * for each of these LPM set the matrix for iversion, i.e. mu/nu in magnetic case
     */
    std::map<UInt, LocPointMapped > allLPM;
    std::map<UInt, LocPointMapped > midpointLPM;
    getLPMMaps(allLPM, midpointLPM);

    UInt dim = GetVecSize();

    UInt size = allLPMmap_.size();
    UInt chunksize;
    UInt numT;
    DetermineChunkSizeAndThreads(size, chunksize, numT);


#ifdef USE_OPENMP
#pragma omp parallel num_threads(numT)
    {
//      UInt size = allLPM.size();
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
//      UInt numT = CFS_NUM_THREADS;
      UInt aThread = omp_get_thread_num();
      //std::cout << "My number: " << aThread << std::endl;
//      UInt chunksize = std::floor(size/numT);

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
          LOG_DBG(coeffunctionhyst_helper) << "Single field case > small signal tensor cannot depend on coupling";
          reuse = true;
          // calling setMatrixForLocalInversion with reuse = true will simply set the flag
          // eps_mu_local_Set_ to true
          setMatrixForLocalInversion(mu, nu, LPMit->first, reuse);
        } else if((strainForm_ <= -1 )&&(materialTensorsSetOnce_)){
          LOG_DBG(coeffunctionhyst_helper) << "Small signal tensor shall not be used; mu/nu independent on M/P; reuse same value";
          reuse = true;
          // calling setMatrixForLocalInversion with reuse = true will simply set the flag
          // eps_mu_local_Set_ to true
          setMatrixForLocalInversion(mu, nu, LPMit->first, reuse);
        } else {
          LOG_DBG(coeffunctionhyst_helper) << "Pre-Compute material tensor for local inversion";
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

          setMatrixForLocalInversion(mu, nu, LPMit->first, reuse);
          LOG_DBG(coeffunctionhyst_helper) << "Computed material tensor for storageIdx="<<LPMit->first<<": "<<mu.ToString();
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
          LOG_DBG(coeffunctionhyst_helper) << "Single field case > small signal tensor cannot depend on coupling";
          reuse = true;
          // calling setMatrixForLocalInversion with reuse = true will simply set the flag
          // eps_mu_local_Set_ to true
          setMatrixForLocalInversion(mu, nu, LPMit->first, reuse);
        } else if((strainForm_ == -1 )&&(materialTensorsSetOnce_)){
          LOG_DBG(coeffunctionhyst_helper) << "Small signal tensor shall not be used; mu/nu independent on M/P; reuse same value";
          reuse = true;
          // calling setMatrixForLocalInversion with reuse = true will simply set the flag
          // eps_mu_local_Set_ to true
          setMatrixForLocalInversion(mu, nu, LPMit->first, reuse);
        } else {
          LOG_DBG(coeffunctionhyst_helper) << "Pre-Compute material tensor for local inversion";
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

          setMatrixForLocalInversion(mu, nu, LPMit->first, reuse);
          LOG_DBG(coeffunctionhyst_helper) << "Computed material tensor for storageIdx="<<LPMit->first<<": "<<mu.ToString();
        }
      }
    }
#endif

    materialTensorsSetOnce_ = true;
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
    formerHystHelperParamsSet_ = false;
    hystModelTraced_ = false;
    globalFPset_ = false;

    lastTSStorageInitialized_ = false;
    lastItStorageInitialized_ = false;
    backupStorageInitialized_ = false;
    estimatorStorageInitialized_ = false;
    strainRequired_ = false;

    DEBUG_infoOutput_ = false;
    DEBUG_debugOutput_ = false;
    DEBUG_printCallingOrder_ = false;

    // dim_ is the dim of the output retrieved by GetVector
		// not tha same as HYSTERESIS_DIM_ that determines whether we use scalar or vector model
		dim_ = dependency1->GetVecSize();

    if(tensorType_ == AXI){
      EXCEPTION("Coef functions hyst does only support plane 2d and full 3d models (at the moment)")
    }

    hyst_ = NULL;
    hystStrain_ = NULL;
    
		matType_ = matType;
		material_ = material;
		ptFeSpace_ = ptFeSpace;

    // default values used for tracing until nov 2020; 
    // worked acceptably well for most cases but failed for tracing extended Mayergoyz model
    // can now be set via mat.xml
    // > init must be done before call to readInMaterial!
    trace_forceRetracing_ = false; // false: load trace from if possible; true: retrace (but only once for each material!)
    trace_forceCentral_ = false;
    trace_JacResolution_ = 1e-5;
    
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

    // new flag for localized fixed-point methods > see header for comments
    evalJacAtMidpointOnly_ = false;
    
		// do not output switching state; that is very costly and should only be done
		// for debugging or special figure computation
		RUN_allowBMP_ = false;

    hystOperatorLocked_ = false;

    // default values
    useFPH_ = false;
    skipStorage_ = false;
    
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
		if (material_->GetClass() == ELECTROSTATIC) {
			rev_mat_fac_ = 8.854187817e-12; //eps0
			PDEName_ = "Electrostatic";

      // no longer used > CoefFunctionHystMat, CoefFunctionHystOuptut etc have their own
      // values for this;
      // small signal tensor = permittivity
      // however, value is used for TestInversion function
      PtrCoefFct permittivity = material_->GetTensorCoefFnc(ELEC_PERMITTIVITY_TENSOR,tensorType_,
              Global::REAL, false);

      eps_nu_base_ = Matrix<Double>(dim_, dim_);
      permittivity->GetTensor(eps_nu_base_, lpm);

      // in electrostatics we only need eps
      eps_mu_base_ = eps_nu_base_;

      PtrCoefFct permittivityFULL = material_->GetTensorCoefFnc(ELEC_PERMITTIVITY_TENSOR,FULL,
              Global::REAL, false);
      
      eps_nu_baseFULL_ = Matrix<Double>(3, 3);
      permittivityFULL->GetTensor(eps_nu_baseFULL_, lpm);

      // in electrostatics we only need eps
      eps_mu_baseFULL_ = eps_nu_baseFULL_;

	  }
		else if (material_->GetClass() == ELECTROMAGNETIC) {
			rev_mat_fac_ = 795774.7155; //nu0
			PDEName_ = "Electromagnetics";
      needsInversion_ = true;
      POL_setWithFlux_ = true; // only for scalar model with vecttor extension
      // no longer used > CoefFunctionHystMat, CoefFunctionHystOuptut etc have their own
      // values for this;
      // however, value is used for TestInversion function
      // small signal tensor = permeability
      PtrCoefFct permeability = material_->GetTensorCoefFnc(MAG_PERMEABILITY_TENSOR,tensorType_,
              Global::REAL, false);

      eps_mu_base_ = Matrix<Double>(dim_, dim_);
      permeability->GetTensor(eps_mu_base_, lpm);

      PtrCoefFct reluctivity = material_->GetTensorCoefFnc(MAG_RELUCTIVITY_TENSOR,tensorType_,
			Global::REAL, false);
      
      eps_mu_base_ = Matrix<Double>(dim_, dim_);
      permeability->GetTensor(eps_mu_base_, lpm);

      eps_nu_base_ = Matrix<Double>(dim_, dim_);
      reluctivity->GetTensor(eps_nu_base_, lpm);
      
      // for inversion test with dimOperator > dim_
      PtrCoefFct permeabilityFULL = material_->GetTensorCoefFnc(MAG_PERMEABILITY_TENSOR,FULL,
              Global::REAL, false);
      PtrCoefFct reluctivityFULL = material_->GetTensorCoefFnc(MAG_RELUCTIVITY_TENSOR,FULL,
              Global::REAL, false);
      eps_mu_baseFULL_ = Matrix<Double>(3, 3);
      permeabilityFULL->GetTensor(eps_mu_baseFULL_, lpm);
      
      eps_nu_baseFULL_ = Matrix<Double>(3, 3);
      reluctivityFULL->GetTensor(eps_nu_baseFULL_, lpm);
		} else {
			EXCEPTION("Currently only Electrostatics and Electromagnetics are supported");
		}
    //    std::cout << "StrainForm: " << strainForm << std::endl;
    //    std::cout << "COUPLED_useStrainForm_: " << COUPLED_useStrainForm_ << std::endl;
    //
    //

    bool testVectorVersionOfScalarModel = false;
    if(testVectorVersionOfScalarModel){
      UInt numTests = 100;
      Hysteresis* hystTEST = new Preisach(1, POL_operatorParams_, POL_weightParams_, true, false);
      Vector<Double> testInput = Vector<Double>(dim_);
      Vector<Double> testOutput = Vector<Double>(dim_);
      Vector<Double> testOutput2 = Vector<Double>(dim_);
      Vector<Double> diff = Vector<Double>(dim_);
      Double testInputScal, testOutputScal;
      int successFlag = 0;
      bool overwrite = false;

      for(UInt i = 0; i < numTests; i++){
        testInput.Init();
        testInput[0] = POL_operatorParams_.inputSat_*( 2.0*static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - 1.0 );
        testInput[1] = POL_operatorParams_.inputSat_*( 2.0*static_cast <float> (rand()) / static_cast <float> (RAND_MAX) - 1.0 );

        // evaluate scalar model as usual by computing a scalar input first and
        // calling computeValueAndUpdate afterwards
        POL_operatorParams_.fixDirection_.Inner(testInput,testInputScal);
        testOutputScal = hystTEST->computeValueAndUpdate(testInputScal, 0, overwrite, successFlag);

        testOutput.Init();
        testOutput.Add(testOutputScal,POL_operatorParams_.fixDirection_);

        // now use vector version
        testOutput2 = hystTEST->computeValue_vec(testInput,0,overwrite,false,successFlag);

        diff.Add(1.0,testOutput,-1.0,testOutput2);

        std::cout << "Input to scalar hyst operator: " << testInput.ToString() << std::endl;
        std::cout << "Output of scalar using old approach: " << testOutput.ToString() << std::endl;
        std::cout << "Output of scalar using new approach: " << testOutput2.ToString() << std::endl;
        std::cout << "Output difference: " << diff.ToString() << std::endl;
      }
      EXCEPTION("Stop after test");
    }

    bool traceHystModelForTesting = true;
    bool stopAfterTracing = !true;
    if(!hystModelTraced_ && traceHystModelForTesting){
      UInt baseSteps = 100;
      // get some additinal info if we are already traverse loop
      Double negCoercivity = 0.0;
      Double maxPolarization = 0.0;

//      TraceHystOperator(baseSteps, maxSlopeGlobal_, minSlopeGlobal_, negCoercivity, maxPolarization);
      TraceHystOperatorVector(baseSteps, maxSlopeGlobal_, minSlopeGlobal_, negCoercivity, maxPolarization);
      hystModelTraced_ = true;
      if(stopAfterTracing){
        exit(0);
      }
    }

    bool testBasicMatrixRoutines = !true;
    if(testBasicMatrixRoutines){
      Matrix<Double> testMatrix = Matrix<Double>(3,3);
      testMatrix[0][0] = 1;
      testMatrix[1][0] = 2;
      testMatrix[0][1] = 3;
      testMatrix[2][0] = 4;
      testMatrix[0][2] = 5;
      testMatrix[1][1] = 6;
      testMatrix[2][1] = 7;
      testMatrix[1][2] = 8;
      testMatrix[2][2] = 9;

      Double max = testMatrix.GetMax();
      Double min = testMatrix.GetMin();

      std::cout << "Testmatrix 1: " << testMatrix.ToString() << std::endl;
      std::cout << "max: " << max << std::endl;
      std::cout << "min: " << min << std::endl;

      testMatrix[0][0] = -1;
      testMatrix[1][0] = -2;
      testMatrix[0][1] = -3;
      testMatrix[2][0] = 4;
      testMatrix[0][2] = 5;
      testMatrix[1][1] = -6;
      testMatrix[2][1] = 7;
      testMatrix[1][2] = 8;
      testMatrix[2][2] = -9;

      max = testMatrix.GetMax();
      min = testMatrix.GetMin();

      std::cout << "Testmatrix 2: " << testMatrix.ToString() << std::endl;
      std::cout << "max: " << max << std::endl;
      std::cout << "min: " << min << std::endl;

      testMatrix[0][0] = rand()%100-50.0;
      testMatrix[1][0] = rand()%100-50.0;
      testMatrix[0][1] = rand()%100-50.0;
      testMatrix[2][0] = rand()%100-50.0;
      testMatrix[0][2] = rand()%100-50.0;
      testMatrix[1][1] = rand()%100-50.0;
      testMatrix[2][1] = rand()%100-50.0;
      testMatrix[1][2] = rand()%100-50.0;
      testMatrix[2][2] = rand()%100-50.0;

      max = testMatrix.GetMax();
      min = testMatrix.GetMin();

      std::cout << "Testmatrix 2: " << testMatrix.ToString() << std::endl;
      std::cout << "max: " << max << std::endl;
      std::cout << "min: " << min << std::endl;

      Matrix<Double> diag;
      testMatrix.GetDiagInMatrix(diag);
      max = diag.GetMax();
      min = diag.GetMin();

      std::cout << "diag of Testmatrix 2: " << diag.ToString() << std::endl;
      std::cout << "max(diag): " << max << std::endl;
      std::cout << "min(diag): " << min << std::endl;


      exit(0);
    }

	}

  void CoefFunctionHyst::ReadAndSetWeights(BaseMaterial* const material, bool setForStrains, int forcedPreisachResolutionForTests){
    ParameterPreisachWeights paramSet = ParameterPreisachWeights();

    // use same offset as in XMLMaterialHandler.cc
    int enumOffset = 0;
    if(setForStrains){
      enumOffset = 100;
    }

    std::string usedHystModel;
    material->GetString(usedHystModel, MaterialType(HYST_MODEL+enumOffset));

    int isPreisachType_Int = 1;
    material->GetScalar(isPreisachType_Int, MaterialType(HYST_TYPE_IS_PREISACH+enumOffset));

    if(isPreisachType_Int == 0){
      std::cout << "Used hysteresis model is not a Preisach-based model; no Preisach weights will be read in" << std::endl;
      return;
    }

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
        /* MUDAT and EXTENDED MUDAT
         * mu = A/(1 + pow( pow((alpha+beta+h1)*sigma1,2) + pow((alpha-beta-h2)*sigma2,2),eta))
         * 
         * In the parameter determination script, the Preisach plane is set to the range -0.5 to +0.5.
         * To use these paramater for the implemented models going from -1 to +1, we have to scale appropriately.
         * A_cfs = A_script/4
         * eta_cfs = eta_script
         * h1_cfs = 2*h1_script
         * h2_cfs = 2*h2_script
         * sigma1_cfs = sigma1_script/2
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
    } else if(weightTypeInt == 4){
      paramSet.weightType_ = "muLorentz";
      // NOTE: we use the same storage as for muDat to save enums
      int forHalfRange;
      material->GetScalar(forHalfRange, MaterialType(PREISACH_WEIGHTS_MUDAT_PARAMSFORHALFRANGE+enumOffset));
      material->GetScalar(paramSet.muLorentz_A_, MaterialType(PREISACH_WEIGHTS_MUDAT_A+enumOffset), Global::REAL);
      material->GetScalar(paramSet.muLorentz_sigma1_, MaterialType(PREISACH_WEIGHTS_MUDAT_SIGMA+enumOffset), Global::REAL);
      material->GetScalar(paramSet.muLorentz_h1_, MaterialType(PREISACH_WEIGHTS_MUDAT_H+enumOffset), Global::REAL);
      if(forHalfRange){
        /*
         * MULORENTZ
         * mu = A /(1 + ((beta+h)/h*sigma)^2) * A /(1 + ((alpha-h)/h*sigma)^2)
         * 
         * In the parameter determination script, the Preisach plane is set to the range -0.5 to +0.5.
         * To use these paramater for the implemented models going from -1 to +1, we have to scale appropriately.
         * 
         * A_cfs = A_script/2 > only factor two compared to muDat but this is due to using B^2 in formula
         * h_cfs = 2*h_script
         * sigma_cfs = sigma_script
         */
        paramSet.muLorentz_A_ /= 2;
        //paramSet.muLorentz_sigma1_ /= 1;
        paramSet.muLorentz_h1_ *= 2;
      }
      // in order to use the same full function (with h1,h2,sigma1 and sigma2) for evaluation, set h2 and sigma2 to
      // h1 and sigma1 respectively
      paramSet.muLorentz_sigma2_ = paramSet.muLorentz_sigma1_;
      paramSet.muLorentz_h2_ = paramSet.muLorentz_h1_;
    } else if(weightTypeInt == 5){
      paramSet.weightType_ = "muLorentzExtended";
      // NOTE: we use the same storage as for muDat to save enums
      int forHalfRange;
      material->GetScalar(forHalfRange, MaterialType(PREISACH_WEIGHTS_MUDAT_PARAMSFORHALFRANGE+enumOffset));
      material->GetScalar(paramSet.muLorentz_A_, MaterialType(PREISACH_WEIGHTS_MUDAT_A+enumOffset), Global::REAL);
      material->GetScalar(paramSet.muLorentz_sigma1_, MaterialType(PREISACH_WEIGHTS_MUDAT_SIGMA+enumOffset), Global::REAL);
      material->GetScalar(paramSet.muLorentz_sigma2_, MaterialType(PREISACH_WEIGHTS_MUDAT_SIGMA2+enumOffset), Global::REAL);
      material->GetScalar(paramSet.muLorentz_h1_, MaterialType(PREISACH_WEIGHTS_MUDAT_H+enumOffset), Global::REAL);
      material->GetScalar(paramSet.muLorentz_h2_, MaterialType(PREISACH_WEIGHTS_MUDAT_H2+enumOffset), Global::REAL);
      
      if(forHalfRange){
        /*
         * MULORENTZ and EXTENDED MULORENTZ
         * mu = A /(1 + ((beta+h1)/h1*sigma1)^2) * A /(1 + ((alpha-h2)/h2*sigma2)^2)
         * 
         * In the parameter determination script, the Preisach plane is set to the range -0.5 to +0.5.
         * To use these paramater for the implemented models going from -1 to +1, we have to scale appropriately.
         * 
         * A_cfs = A_script/2 > only factor two compared to muDat but this is due to using B^2 in formula
         * h_cfs = 2*h_script
         * sigma_cfs = sigma_script
         */
        paramSet.muLorentz_A_ /= 2;
        //paramSet.muLorentz_sigma1_ /= 1;
        paramSet.muLorentz_h1_ *= 2;
        //paramSet.muLorentz_sigma2_ /= 1;
        paramSet.muLorentz_h2_ *= 2;
      }
      
    }  else {
      EXCEPTION("Weight type unknown");
    }

    // NEW: add additional anhysteretic curve to Preisach models
    material->GetScalar(paramSet.anhysteretic_a_ , MaterialType(PREISACH_WEIGHTS_ANHYST_A+enumOffset), Global::REAL);
    material->GetScalar(paramSet.anhysteretic_b_ , MaterialType(PREISACH_WEIGHTS_ANHYST_B+enumOffset), Global::REAL);
    material->GetScalar(paramSet.anhysteretic_c_ , MaterialType(PREISACH_WEIGHTS_ANHYST_C+enumOffset), Global::REAL);
    material->GetScalar(paramSet.anhysteretic_d_ , MaterialType(PREISACH_WEIGHTS_ANHYST_D+enumOffset), Global::REAL);

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
       *  a*atan(b*e+d) + c*e
       * 
       * a_cfs = 2*a_script (script assume anhystpart to have max ampl of 0.5 instead of 1)
       * b_cfs = b_script/2 (script multiplies b with e_norm in range [-0.5,0.5])
       * c_cfs = c_script*2 (for scaling y to +1,-1) /2 (as e is double compared to script)
       * d_cfs = d_script
       */
      paramSet.anhysteretic_a_ *= 2;
      paramSet.anhysteretic_b_ /= 2;
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

    bool useTMPstorage = false;
    if(forcedPreisachResolutionForTests > 0 ){
      if( (paramSet.weightType_ == "givenTensor")&&(UInt(forcedPreisachResolutionForTests) != paramSet.numRows_) ){
        WARN("Test-Option forcedPreisachResolutionForTests is not compatible with weight type givenTensor.");
      } else {
        paramSet.numRows_ = UInt(forcedPreisachResolutionForTests);
        useTMPstorage = true;
      }
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
//    std::cout << "intOverWeights: " << intOverWeights << std::endl;
    
    
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
//      paramSet.anhystOnly_ = false;
      Double elemArea = 2.0/(paramSet.numRows_)*2.0/(paramSet.numRows_);
      intOverWeights *= elemArea/2; // factor of 1/2 as we have to integrate over the triangle area
      //std::cout << "intOverWeights = " << intOverWeights << std::endl;
      bool scalingRequired = true;

      if(usedHystModel == "vectorPreisach_Mayergoyz"){
        int isIsotropic = 1;
        material->GetScalar(isIsotropic, MaterialType(PREISACH_MAYERGOYZ_ISOTROPIC+enumOffset));

        if(isIsotropic == 0){
//        if( (dim_ != 2) || (isIsotropic == 0)){
          EXCEPTION("Mayergoyz vector model currently only implemented for isotropic materials");
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
//            std::cout << "transformPreisachWeightsForIsotropicVectorCase" << std::endl;
            paramSet.weightTensor_ = transformPreisachWeightsForIsotropicVectorCase(&paramSet);
            
            // transformPreisachWeightsForIsotropicVectorCase does only compute
            // the upper triangle > mirror along alpha = beta to obtain the full numRows_ x numRows_ tensor
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
//          std::cout << "Weights after forcing symmetry: " << paramSet.weightTensor_.ToString() << std::endl;
        }
      }

      // store integral over weights;
      // for our implementation we want to have the integral over all weights to be 1 which is why we scale
      // the weights; the implementation that derives the muDat parameter, however, allows the integral over
      // the weights to be different from 1, especially in case of anhysteretic parts being present; in the later
      // case, the saturation value outputSat refers to the sum of preisach in sat and anhyst part in sat;
      // to obtain the same curves we have two options:
      // a) allow int over weights to be different from 1, too
      // b) define saturation value of preisach alone to be outputSat - anhystPart[inputsat]
      // > currently using option b
      // > this option can be turned of in mat.xml by setting anhystIncludedInSatValue to false
      //    by doing so, outputSat will be the saturation value of the preisach operator AND the scaling factor
      //    for the anhyst part > where is the switch for that?
      // Note: in case of the Mayergoyz vector model, intOverWeights refers to the original scalar weights not to
      //    the transformed ones! the weights for the single scalar models inside the Mayergoyz model do not sum up to 1!
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
              paramSet.anhysteretic_b_,paramSet.anhysteretic_c_,paramSet.anhysteretic_d_);
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
    if(!useTMPstorage){
      if(setForStrains){
        STRAIN_weightParams_ = paramSet;
      } else {
        POL_weightParams_ = paramSet;
      }
    } else {
      if(setForStrains){
        STRAIN_weightParamsForTesting_ = paramSet;
      } else {
        POL_weightParamsForTesting_ = paramSet;
      }
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
		material->GetString(dimTypeStr, MaterialType(HYSTERESIS_DIM+enumOffset));

    if (dimTypeStr == "SCALAR") {
			paramSet.methodType_ = 0;
		} else if (dimTypeStr == "VECTOR") {
			paramSet.methodType_ = 1;
		}

    // print warnings (e.g., due to unsymmetric weights; usually set to false)
    paramSet.printWarnings_ = false;

    material->GetString(paramSet.methodName_, MaterialType(HYST_MODEL+enumOffset));
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
    // new: always set x-direction as default; will be used for test hystoperator
    //  if flag Test1D = true
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

    paramSet.startingAxisMG_ = Vector<Double>(dim_);
    paramSet.startingAxisMG_.Init();

    // USED FOR BOTH MODELS
    paramSet.fieldsAlignedAboveSat_ = true;
		paramSet.hystOutputRestrictedToSat_ = true;
    paramSet.hasInverseModel_ = true;

    // for debugging > should be put into mat file!
    paramSet.checkInversionResult_ = !true;

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
      Double dirNorm = paramSet.fixDirection_.NormL2();

      if(paramSet.fixDirection_.NormL2() == 0){
        WARN("Zero direction specified; taking default = x-direction");
        paramSet.fixDirection_[0] = 1.0;
      } else {
        paramSet.fixDirection_[0]/=dirNorm;
        paramSet.fixDirection_[1]/=dirNorm;
        if (dim_ == 3) {
          paramSet.fixDirection_[2]/=dirNorm;
        }
      }
      
      paramSet.fieldsAlignedAboveSat_ = true;
      paramSet.inputForAlignment_ = paramSet.inputSat_;
//      std::cout << "Scalar model - using direction (normalized): " << paramSet.fixDirection_.ToString() << std::endl;

    } else if(paramSet.methodName_ == "JilesAtherton"){
      Double dummyStorage;
      material->GetScalar(dummyStorage, MaterialType(JILES_TEST), Global::REAL);
      std::cout << "Test for JilesAtherton Parameters; found " << dummyStorage << std::endl;
      exit(0); // no model implemented yet
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

      int scaleToSat = 1;
      material->GetScalar(scaleToSat, MaterialType(SCALETOSAT+enumOffset));

//      std::cout << "Scale to Sat: " << scaleToSat << std::endl;

      if(scaleToSat == 1){
        paramSet.scaleUpToSaturation_ = true;
      } else {
        paramSet.scaleUpToSaturation_ = false;
      }
      
      /*
       * changes regarding treatment of revised model;
       * actually we have an alignment of input and output but not necessarily at input saturation directly
       * for kappa_revised >= 1, alignment is obtained at input saturation (angular lag is vanishing, whole Preisach plane filled in terms of 
       * rotational operator)
       * for kappa_revised (aka rotResistance_) < 1, alignment is obtained at input saturation times 1/kappa_revised (angular lag is
       * vanishing at input saturation; Preisach plane is filled completely at 1/kappa_revised
       * if scaleUpToSaturation_ is set to true, alignment is given at input saturation for all kappa_revised
       * -->  fieldsAlignedAboveSat_ should always be set to true, but to distinguish the different cases,
       * another parameter is used that defines the required input field (inputForAlignment_)
       */
      paramSet.fieldsAlignedAboveSat_ = true;
      if((paramSet.evalVersion_ == 1)||(paramSet.evalVersion_ == 10)){
        paramSet.isClassical_=true;
      } else {
        paramSet.isClassical_=false;
      }
//      std::cout << "1. - ReadAndSetParamsForHystOperator" << std::endl;
//      std::cout << "set value of paramSet.inputForAlignment_" << std::endl;
      if(paramSet.isClassical_==true){
//        std::cout << "Classic" << std::endl;
        paramSet.inputForAlignment_ = paramSet.inputSat_;
      } else {
//        std::cout << "Revised" << std::endl;
        if(paramSet.scaleUpToSaturation_ == true){
//          std::cout << "ForcedAlignment" << std::endl;
          paramSet.inputForAlignment_ = paramSet.inputSat_;
        } else {
          if(paramSet.rotResistance_ >= 1.0){
//            std::cout << "KappaRes=" << paramSet.rotResistance_ <<" >= 1" << std::endl;
            paramSet.inputForAlignment_ = paramSet.inputSat_;
          } else {
//            std::cout << "KappaRes=" << paramSet.rotResistance_ <<" < 1" << std::endl;
            paramSet.inputForAlignment_ = paramSet.inputSat_/paramSet.rotResistance_;
          }
        }
      }
      std::cout << "paramSet.inputSat_ = " << paramSet.inputSat_ << std::endl;
      std::cout << "paramSet.inputForAlignment_ = " << paramSet.inputForAlignment_ << std::endl;

      std::cout << "paramSet.fieldsAlignedAboveSat_? " << paramSet.fieldsAlignedAboveSat_ << std::endl;
      
    } else if(paramSet.methodName_ == "vectorPreisach_Mayergoyz"){
      int isIsotropic = 1;
      material->GetScalar(isIsotropic, MaterialType(PREISACH_MAYERGOYZ_ISOTROPIC+enumOffset));
      if(isIsotropic == 0){
//        if( (dim_ != 2) || (isIsotropic == 0)){
        EXCEPTION("Mayergoyz vector model currently only implemented for 2d isotropic materials");
      }
      paramSet.fieldsAlignedAboveSat_ = false;
      paramSet.inputForAlignment_ = std::numeric_limits<Double>::max();
      paramSet.hasInverseModel_ = false;
      paramSet.isIsotropic_ = true;

      paramSet.outputClipping_ = 0;
      material->GetScalar(paramSet.outputClipping_, MaterialType(PREISACH_MAYERGOYZ_CLIPOUTPUT+enumOffset));

			if(paramSet.outputClipping_ == 0){
				paramSet.hystOutputRestrictedToSat_ = false;
			}

      paramSet.numDirections_ = 0;
      material->GetScalar(paramSet.numDirections_, MaterialType(PREISACH_MAYERGOYZ_NUM_DIR+enumOffset));

			material->GetScalar(paramSet.startingAxisMG_[0], MaterialType(MAYERGOYZ_STARTAXIS_X+enumOffset), Global::REAL);
			material->GetScalar(paramSet.startingAxisMG_[1], MaterialType(MAYERGOYZ_STARTAXIS_Y+enumOffset), Global::REAL);
			if (dim_ == 3) {
				material->GetScalar(paramSet.startingAxisMG_[2], MaterialType(MAYERGOYZ_STARTAXIS_Z+enumOffset), Global::REAL);
			}
      
      material->GetScalar(paramSet.lossParam_a, MaterialType(MAYERGOYZ_LOSSPARAM_A+enumOffset), Global::REAL);
      material->GetScalar(paramSet.lossParam_b, MaterialType(MAYERGOYZ_LOSSPARAM_B+enumOffset), Global::REAL);

      int useAbsoluteValueOfdPhi_int;
      int normalizeXInExponentOfG_int;
      material->GetScalar(useAbsoluteValueOfdPhi_int, MaterialType(MAYERGOYZ_USEABSDPHI+enumOffset));
      material->GetScalar(normalizeXInExponentOfG_int, MaterialType(MAYERGOYZ_NORMALIZEXINEXP+enumOffset));
      material->GetScalar(paramSet.restrictionOfPsi_, MaterialType(MAYERGOYZ_RESTRICTIONOFPSI+enumOffset));
      if(useAbsoluteValueOfdPhi_int == 1){
        paramSet.useAbsoluteValueOfdPhi_ = true;
      } else {
        paramSet.useAbsoluteValueOfdPhi_ = false;
      }
      if(normalizeXInExponentOfG_int == 1){
        paramSet.normalizeXInExponentOfG_ = true;
      } else {
        paramSet.normalizeXInExponentOfG_ = false;
      }
      
      material->GetScalar(paramSet.scalingFactorXInExponent_, MaterialType(MAYERGOYZ_SCALINGOFXINEXP+enumOffset), Global::REAL);

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
     * 6.11.2020
     * Get general parameter required for tracing
     * motivation: tracing is one of the most critical parts for solving the hysteretic 
     * equation systems; if not done properly, the derived slopes - which in terms are
     * used to derive the Fixed-Point contraction factors - can be too large or too small
     * causing either a too slow convergence or no convergence at all
     * Most of the older testcases were successful for 
     * - forward differences > set by switching forceCentral in setFlag() function below
     * - JacResolution (min resolution for computation of jacobian) = 1e-5
     * However, this combination lead to serious problems with the Mayergoyz model so that
     * adaptions had to be made
     */
    int forceRetraceTmp; 
    int forceCentralTmp; 
    material->GetScalar(trace_JacResolution_, TRACE_JAC_RESOLUTION, Global::REAL);
    material->GetScalar(forceRetraceTmp, TRACE_FORCE_RETRACING);
    material->GetScalar(forceCentralTmp, TRACE_FORCE_CENTRALDIFF);
    
    if(forceRetraceTmp == 1){
      trace_forceRetracing_ = true;
    } else {
      trace_forceRetracing_ = false;
    }
    if(forceCentralTmp == 1){
      trace_forceCentral_ = true;
    } else {
      trace_forceCentral_ = false;
    }
    
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
//    if ( (POL_operatorParams_.methodName_ == "vectorPreisach_Sutor") || (POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") ){
//      // inversion via LM > set parameter
    // NEW: read in hysteresis inversion also for scalar model to allow different tests
    //      default inversion for scalar model should still be EverettBasedInversion, though!
      if (material_->GetClass() == ELECTROMAGNETIC) {
        inversionSet_ = true;
        Integer invMat, invMethod, maxNumReg, maxNumLS;
        material->GetScalar(invMat, MAX_NUM_IT_HYST_INV);
        InversionParams_.maxNumIts = (UInt) invMat;

        material->GetScalar(invMethod, VEC_HYST_INV_METHOD);
        InversionParams_.inversionMethod = static_cast<localInversionFlag>(invMethod);

        std::string usedMethod;
        Enum2String(InversionParams_.inversionMethod,usedMethod);
//        switch (InversionParams_.inversionMethod) {
//          case 0:
//            usedMethod = "Levenberg-Marquardt"; break;
//          case 1:
//            usedMethod = "Unregularized Newton"; break;
//          case 2:
//            usedMethod = "Jacobian-Free-Newton-Krylov"; break;
//          case 3:
//            usedMethod = "projected LM"; break;
//          case 4:
//            usedMethod = "Everett based (scalar model only)"; break;
//          case 5:
//            usedMethod = "Fixpoint"; break;
//          default:
//            usedMethod = "Undefined!";
//        }
        LOG_DBG(coeffunctionhyst_main) << "Defined inversion method: " << usedMethod;

        if(InversionParams_.inversionMethod == LOCAL_FIXPOINT){
          // local fp method
          // Trace hyst operator first
//          std::cout << "Setting parameter for local inversion - trace hyst operator for inversion via fp" << std::endl;
          if(!hystModelTraced_){
            UInt baseSteps = 100;
            // get some additinal info if we are already traverse loop
            Double negCoercivity = 0.0;
            Double maxPolarization = 0.0;

//          if (POL_operatorParams_.methodName_ == "scalarPreisach") {
//            TraceHystOperator(baseSteps, maxSlopeGlobal_, minSlopeGlobal_, negCoercivity, maxPolarization);
//          } else {
            TraceHystOperatorVector(baseSteps, maxSlopeGlobal_, minSlopeGlobal_, negCoercivity, maxPolarization);
//          }

//            //          TraceHystOperator(baseSteps, maxSlopeGlobal_, minSlopeGlobal_, negCoercivity, maxPolarization);
//            TraceHystOperatorVector(baseSteps, maxSlopeGlobal_, minSlopeGlobal_, negCoercivity, maxPolarization);
            hystModelTraced_ = true;
          }
          Double contractionFactorLocalInversion;
          material->GetScalar(contractionFactorLocalInversion, HYST_INV_FP_SAFETYFACTOR, Global::REAL);
          
          // the global system has its own safety factors (localFPFactor_C,globalFPFactor_C) see coeffnchyst.hh; 
          // note that the "local" of localFPFactor_C refers to the localized FP method, not
          // to the local inversion on element level!
          InversionParams_.safetyFactor_C = contractionFactorLocalInversion;
          InversionParams_.minSlopeHystOperator = minSlopeGlobal_;
          InversionParams_.maxSlopeHystOperator = maxSlopeGlobal_;
        }

        material->GetScalar(InversionParams_.tolH, RES_TOL_H_HYST_INV, Global::REAL);
        material->GetScalar(InversionParams_.tolB, RES_TOL_B_HYST_INV, Global::REAL);

        int tolH_useAsRelativeNormInt = 0;
        material->GetScalar(tolH_useAsRelativeNormInt, RES_TOL_H_HYST_INV_ISREL);
        if(tolH_useAsRelativeNormInt == 1){
          InversionParams_.tolH_useAsRelativeNorm = true;
        } else {
          InversionParams_.tolH_useAsRelativeNorm = false;
        }
        int tolB_useAsRelativeNormInt = 0;
        material->GetScalar(tolB_useAsRelativeNormInt, RES_TOL_B_HYST_INV_ISREL);
        if(tolB_useAsRelativeNormInt == 1){
          InversionParams_.tolB_useAsRelativeNorm = true;
        } else {
          InversionParams_.tolB_useAsRelativeNorm = false;
        }
        
        material->GetScalar(maxNumReg, MAX_NUM_REG_IT_HYST_INV);
        InversionParams_.maxNumRegIts = (UInt) maxNumReg;
        material->GetScalar(InversionParams_.alphaRegStart, ALPHA_REG_HYST_INV, Global::REAL);
        material->GetScalar(InversionParams_.alphaRegMin, ALPHA_REG_MIN_HYST_INV, Global::REAL);
        material->GetScalar(InversionParams_.alphaRegMax, ALPHA_REG_MAX_HYST_INV, Global::REAL);

        material->GetScalar(InversionParams_.trustLow, TRUST_LOW_HYST_INV, Global::REAL);
        material->GetScalar(InversionParams_.trustMid, TRUST_MID_HYST_INV, Global::REAL);
        material->GetScalar(InversionParams_.trustHigh, TRUST_HIGH_HYST_INV, Global::REAL);

        material->GetScalar(maxNumLS, MAX_NUM_LS_IT_HYST_INV);
        InversionParams_.maxNumLSIts = (UInt) maxNumLS;
        material->GetScalar(InversionParams_.alphaLSMin, ALPHA_LS_MIN_HYST_INV, Global::REAL);
        material->GetScalar(InversionParams_.alphaLSMax, ALPHA_LS_MAX_HYST_INV, Global::REAL);

        material->GetScalar(InversionParams_.jacImplementation, JAC_IMPLEMENTATION_HYST_INV);

        std::string usedJacImplementation;
        switch (InversionParams_.jacImplementation) {
          case -1:
            usedJacImplementation = "Jacobi-Free"; break;
          case 1:
            usedJacImplementation = "Central differences"; break;
          case 0:
            usedJacImplementation = "Forward/Backward"; break;
          case 2:
            usedJacImplementation = "Forward/Backward; possibly scaled diagonal"; break;
          default:
            usedJacImplementation = "Undefined!";
        }
        LOG_DBG(coeffunctionhyst_main) << "Defined method for computing jacobian matrix: " << usedJacImplementation;


        material->GetScalar(InversionParams_.jacRes, JAC_RESOLUTION_HYST_INV, Global::REAL);
        int stopLSAtLocalMin = 0;
        material->GetScalar(stopLSAtLocalMin, STOP_INV_LS_AT_LOCAL_MIN);
        if(stopLSAtLocalMin == 1){
          InversionParams_.stopLineSearchAtLocalMin = true;
        } else {
          InversionParams_.stopLineSearchAtLocalMin = false;
        }

        int printWarningsInt = 0;
        material->GetScalar(printWarningsInt, HYST_LOCAL_INVERSION_PRINT_WARNINGS);
        if(printWarningsInt == 1){
          InversionParams_.printWarnings = true;
        } else {
          InversionParams_.printWarnings = false;
        }

        material->GetScalar(InversionParams_.projLM_mu, HYST_INV_PROJLM_MU, Global::REAL);
        material->GetScalar(InversionParams_.projLM_rho, HYST_INV_PROJLM_RHO, Global::REAL);
        material->GetScalar(InversionParams_.projLM_beta, HYST_INV_PROJLM_BETA, Global::REAL);
        material->GetScalar(InversionParams_.projLM_sigma, HYST_INV_PROJLM_SIGMA, Global::REAL);
        material->GetScalar(InversionParams_.projLM_gamma, HYST_INV_PROJLM_GAMMA, Global::REAL);
        material->GetScalar(InversionParams_.projLM_tau, HYST_INV_PROJLM_TAU, Global::REAL);
        material->GetScalar(InversionParams_.projLM_c, HYST_INV_PROJLM_C, Global::REAL);
        material->GetScalar(InversionParams_.projLM_p, HYST_INV_PROJLM_P, Global::REAL);

      } else {
        InversionParams_.inversionMethod = LOCAL_NOTIMPLEMENTED; // undefined
      }
//    } else if(POL_operatorParams_.methodName_ == "scalarPreisach"){
//      inversionSet_ = true; // scalar model needs no additional parameter for inversion, it is always ready
//    }

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
      strainRequired_ = true;

      CouplingParams_.couplingDefined_inMatFile_ = true;

      /*
       * Get additional operators for coupling case
       */
      int strainForm;
      material->GetScalar(strainForm, HYST_STRAIN_FORM);

      /*
       * -2 : not coupled at all;
       * -1 : not coupled via small signal parameter; irreversible strains still present
       *  0 : coupled e-form/h-form (piezo/magstrict)
       *  1 : coupled d-form (piezo)
       *  2 : coupled g-form (magstrict)
       */

      CouplingParams_.strainForm_ = strainForm;
//      if(strainForm > 0){
      if(strainForm != -2){
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
      if(strainRequired_){
        delete[] Si_;
      }

      if(lastItStorageInitialized_){
        delete[] E_B_lastIt_;
        delete[] P_J_lastIt_;
        delete[] E_H_lastIt_;
        if(strainRequired_){
          delete[] Si_lastIt_;
        }
      }
      if(lastTSStorageInitialized_){
        delete[] E_B_lastTS_;
        delete[] P_J_lastTS_;
        delete[] E_H_lastTS_;
        if(strainRequired_){
          delete[] Si_lastTS_;
        }
      }
      if(estimatorStorageInitialized_){
        delete[] E_B_estimatorStep_;
        delete[] P_J_estimatorStep_;
        delete[] E_H_estimatorStep_;
        if(strainRequired_){
          delete[] Si_estimatorStep_;
        }
      }
      if(backupStorageInitialized_){
        delete[] E_B_BAK_;
        delete[] P_J_BAK_;
        delete[] E_H_BAK_;
        if(strainRequired_){
          delete[] Si_BAK_;
        }
      }

      delete[] deltaMat_;
      delete[] deltaMatPrev_;
      if(strainRequired_){
        delete[] deltaMatStrain_;
        delete[] deltaMatStrainPrev_;
      }
      delete[] eps_mu_local_;
      delete[] epsInv_nu_local_;
      delete[] eps_mu_local_Set_;
      delete[] rotatedCouplingTensor_;

      delete[] requiresReeval_;
      delete[] Si_requiresReeval_;
      delete[] deltaMat_requiresReeval_;
      delete[] deltaMatStrain_requiresReeval_;
      delete[] rotatedCouplingTensor_requiresReeval_;

      delete[] takeEstimatedSlope_;
      delete[] FPMaterialTensor_;
      
      if(hyst_ != NULL){
        delete hyst_;
      }
      if(hystStrain_ != NULL){
        delete hystStrain_;
      }
    }
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
        hyst_->UpdateRotationStateWithFluxDensity(latestFlux,eps_mu_local_[i],i);
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

  void CoefFunctionHyst::InitSpecificStorage(int storageLocation){

    if(storageInitialized_ == false){
      WARN("Normal storage has to be initialized first!");
      return;
    }

    LOG_DBG(coeffunctionhyst_main) << "Init specific storage - " << storageLocation;
//    std::cout << "Init specific storage - " << storageLocation << std::endl;

    Vector<Double> zeroVec = Vector<Double>(dim_);
    zeroVec.Init();
    Vector<Double> zeroVecStrain;
    if(dim_ == 2){
      zeroVecStrain = Vector<Double>(3);
    } else {
      zeroVecStrain = Vector<Double>(6);
    }
    zeroVecStrain.Init();

    if(storageLocation == 1){
      if(lastTSStorageInitialized_){
        return;
      }
      E_B_lastTS_ = new Vector<Double>[numStorageEntries_];
      P_J_lastTS_ = new Vector<Double>[numStorageEntries_];
      E_H_lastTS_ = new Vector<Double>[numStorageEntries_];
      if(strainRequired_){
        P_J_forStrains_lastTS_ = new Vector<Double>[numStorageEntries_];
        Si_lastTS_ = new Vector<Double>[numStorageEntries_];
      }

      for (UInt k = 0; k < numStorageEntries_; k++) {
        E_B_lastTS_[k] = zeroVec;
        E_H_lastTS_[k] = zeroVec;
        P_J_lastTS_[k] = zeroVec;
        if(strainRequired_){
          P_J_forStrains_lastTS_[k] = zeroVec;
          Si_lastTS_[k] = zeroVecStrain;
        }
      }

      // init storage for lastTS
      lastTSStorageInitialized_ = true;
    }

    if(storageLocation == 0){
      if(lastItStorageInitialized_){
        return;
      }
      E_B_lastIt_ = new Vector<Double>[numStorageEntries_];
      P_J_lastIt_ = new Vector<Double>[numStorageEntries_];
      E_H_lastIt_ = new Vector<Double>[numStorageEntries_];
      if(strainRequired_){
        Si_lastIt_ = new Vector<Double>[numStorageEntries_];
      }

      for (UInt k = 0; k < numStorageEntries_; k++) {
        E_B_lastIt_[k] = zeroVec;
        E_H_lastIt_[k] = zeroVec;
        P_J_lastIt_[k] = zeroVec;
        if(strainRequired_){
          Si_lastIt_[k] = zeroVecStrain;
        }
      }

      // init storage for lastIt
      lastItStorageInitialized_ = true;
    }

    if(storageLocation == 2){
      if(estimatorStorageInitialized_){
        return;
      }
      E_B_estimatorStep_ = new Vector<Double>[numStorageEntries_];
      P_J_estimatorStep_ = new Vector<Double>[numStorageEntries_];
      E_H_estimatorStep_ = new Vector<Double>[numStorageEntries_];
      if(strainRequired_){
        Si_estimatorStep_ = new Vector<Double>[numStorageEntries_];
      }

      for (UInt k = 0; k < numStorageEntries_; k++) {
        E_B_estimatorStep_[k] = zeroVec;
        P_J_estimatorStep_[k] = zeroVec;
        E_H_estimatorStep_[k] = zeroVec;
        if(strainRequired_){
          Si_estimatorStep_[k] = zeroVecStrain;
        }
      }

      // init storage for estimator
      estimatorStorageInitialized_ = true;
    }

    if(storageLocation == 3){
      if(backupStorageInitialized_){
        return;
      }
      E_B_BAK_ = new Vector<Double>[numStorageEntries_];
      P_J_BAK_ = new Vector<Double>[numStorageEntries_];
      E_H_BAK_ = new Vector<Double>[numStorageEntries_];
      if(strainRequired_){
        Si_BAK_ = new Vector<Double>[numStorageEntries_];
      }

      for (UInt k = 0; k < numStorageEntries_; k++) {
        E_B_BAK_[k] = zeroVec;
        P_J_BAK_[k] = zeroVec;
        E_H_BAK_[k] = zeroVec;
        if(strainRequired_){
          Si_BAK_[k] = zeroVecStrain;
        }
      }

      // init storage for backup
      backupStorageInitialized_ = true;
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

//    std::cout << "4. Init storage" << std::endl;
//    std::cout << "Check if new parameter inputForAlignment_ has been set" << std::endl;
//    std::cout << "POL_operatorParams.inputForAlignment_ = " << POL_operatorParams_.inputForAlignment_ << std::endl;
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
      // Scalar model
      bool isVirgin = true;
      POL_useExtension_ = false;
      //      int useExtensionInt;
      //      material_->GetScalar(useExtensionInt, SCALPREISACH_USE_EXT);
      //      if(useExtensionInt == 1){
      //        POL_useExtension_ = true;
      //      }
      if(POL_useExtension_){
        EXCEPTION("No longer available");
      } else {
        bool ignoreAnhyst = false;
        hyst_ = new Preisach(numHystOperators_, POL_operatorParams_, POL_weightParams_, isVirgin, ignoreAnhyst);
//                (numHystOperators_, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_, POL_weightParams_.weightTensor_, isVirgin,
//                POL_weightParams_.anhysteretic_a_, POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
      }

      POL_operatorParams_.hasInverseModel_ = false;
      // used during testing of hyst operator so it is set here
      // 4.11.2020
      // not a good idea! tolB and tolH are already set at this point
      // this just overwrites the given input
//      InversionParams_.tolB = 1e-10;
//      InversionParams_.tolH = 1e-8; // criterion as used in Preisach.cc for inversion; note: here for actual pol; in Preisach.cc scaled by XSatuated as it calcs with normalized values

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

      
      // setting of fieldsAlignedAboveSat_ is already treated in ReadAndSetParamsForHystOperator
      // the following lines are only for debugging purpose
//      bool useOldAlignmentCriterion = false;
//      if(useOldAlignmentCriterion){
//        std::cout << "POL_operatorParams_.fieldsAlignedAboveSat_ before checking old criterion: "<<POL_operatorParams_.fieldsAlignedAboveSat_<<std::endl;
//        if( (POL_operatorParams_.scaleUpToSaturation_ == false)&&(POL_operatorParams_.isClassical_ == false) ){
//          POL_operatorParams_.fieldsAlignedAboveSat_ = false;
//        }
//        std::cout << "POL_operatorParams_.fieldsAlignedAboveSat_ after checking old criterion: "<<POL_operatorParams_.fieldsAlignedAboveSat_<<std::endl;
//      }
		}
    else if (POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") {
      // basically a scalar model in multiple directions
      // isotropic case: all scalar models are equal (same weights etc)
      // anisotropic case: each model different; choice of directions matters; weights are harder to obtain
      // already treated in ReadAndSetParamsForHystOperator
      //      POL_operatorParams_.fieldsAlignedAboveSat_ = false;
      //      // set dummy value for new parameter inputForAlignment_
      //      POL_operatorParams_.inputForAlignment_ = numeric_limits<Double>::max();
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

//    if ( (POL_operatorParams_.methodName_ == "vectorPreisach_Sutor") || (POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") ){
//      // inversion via LM > set parameter
    //8.6.2020: Also set for scalar preisach model to pass user specified tolerances!
      if (material_->GetClass() == ELECTROMAGNETIC) {
        hyst_->SetParamsForInversion(InversionParams_);
      }
//    }

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
        // already treated in ReadAndSetParamsForHystOperator
        // compare with POL_operatorParams_ above
        //        if( (STRAIN_operatorParams_.scaleUpToSaturation_ == false)&&(STRAIN_operatorParams_.isClassical_ == false) ){
        //          STRAIN_operatorParams_.fieldsAlignedAboveSat_ = false;
        //          
        //        }
        
      }
      else if (STRAIN_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") {
        // basically a scalar model in multiple directions
        // isotropic case: all scalar models are equal (same weights etc)
        // anisotropic case: each model different; choice of directions matters; weights are harder to obtain
        // already treated in ReadAndSetParamsForHystOperator
        //        STRAIN_operatorParams_.fieldsAlignedAboveSat_ = false;
        //        // set dummy value for new parameter inputForAlignment_
        //        STRAIN_operatorParams_.inputForAlignment_ = numeric_limits<Double>::max();
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
          eps_mu_base_.Mult(POL_operatorParams_.fixDirection_,tmp);
          POL_operatorParams_.fixDirection_.Inner(tmp,muForInversion);

          Double scalInput;
          POL_operatorParams_.fixDirection_.Inner(initial_E_B,scalInput);
          Double scalOutput = hyst_->computeInputAndUpdate(scalInput,muForInversion, operatorIdx, overwrite, successCode);

          initial_E_H.Init();
          initial_E_H.Add(scalOutput,POL_operatorParams_.fixDirection_);
        } else {
          bool useEverett = false;
          bool overwrite = false;
          // new type: everett based inversion of mayergoyz vector model
          // as we want to overwrite the memory later on (see below flag overwrite = true)
          // during this initial phase, we have to pass it here as the mayergoyz model based on
          // inverted scalar models overwrites its memory directly in computeInput_vec
          if(InversionParams_.inversionMethod == 10){
            useEverett = true;
            overwrite = true;
          }
          
          initial_E_H = hyst_->computeInput_vec(initial_E_B,operatorIdx,eps_mu_base_,
                  POL_operatorParams_.fieldsAlignedAboveSat_,POL_operatorParams_.hystOutputRestrictedToSat_,successCode,useEverett,overwrite);
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
        if (POL_operatorParams_.methodType_ == 0) {
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
          if(InversionParams_.inversionMethod != LOCAL_NOTIMPLEMENTED){
            initial_P_J = hyst_->computeValue_vec(initial_E_H, k, overwrite, debugOut, successCode);
          } else {
            // in case of the everett based inversion, we are no longer allowed to call computeValue_vec as
            // the retrieved input from computeValue_vec will not return the original input (the inverse of a sum is
            // not the sum of the inverted terms!); instead we compute the output by subtracting mu*initial_E_H from 
            // initial initial_E_B
            initial_P_J.Init(); 
            eps_mu_base_.Mult(initial_E_H,initial_P_J);
            initial_P_J.Add(-1.0,initial_E_B);
            initial_P_J.ScalarMult(-1.0);
          }
          if(CouplingParams_.ownHystOperator_){
            if(InversionParams_.inversionMethod != LOCAL_NOTIMPLEMENTED){
              P_J_forStrains = hystStrain_->computeValue_vec(initial_E_H, k, overwrite, debugOut, successCode);
            } else {
              // here we need an initial input for hystStrain but this would require an initial strain to be given!
              // this is not case; BUT: as we never actually have to invert hystStrain_ (at least not at the moment)
              // we can use hystStrain_ in forward model as usual and just pass the retrieved input to it!
              P_J_forStrains = hystStrain_->computeValue_vec(initial_E_H, k, overwrite, debugOut, successCode);
            }
          } else {
            P_J_forStrains = initial_P_J;
          }
        }
      }
    }

//    std::cout << "Initial output polarization: " << initial_P_J.ToString() << std::endl;

    Vector<Double> zeroVec = Vector<Double>(dim_);
    zeroVec.Init();
    Vector<Double> initial_Si;
    if(strainRequired_){
      initial_Si = ComputeIrreversibleStrains(P_J_forStrains, initial_E_H, zeroVec);
    }
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
    E_H_ = new Vector<Double>[numStorageEntries_];
		P_J_ = new Vector<Double>[numStorageEntries_];
    if(strainRequired_){
      P_J_forStrains_ = new Vector<Double>[numStorageEntries_];
      Si_ = new Vector<Double>[numStorageEntries_];
      deltaMatStrain_ = new Matrix<Double>[numStorageEntries_];
      deltaMatStrainPrev_ = new Matrix<Double>[numStorageEntries_];
    }

    estimatorSet_ = false;

		deltaMat_ = new Matrix<Double>[numStorageEntries_];
    deltaMatPrev_ = new Matrix<Double>[numStorageEntries_];

    eps_mu_local_ = new Matrix<Double>[numStorageEntries_];
    epsInv_nu_local_ = new Matrix<Double>[numStorageEntries_];
    eps_mu_local_Set_ = new bool[numStorageEntries_];

		requiresReeval_ = new bool[numStorageEntries_];
    Si_requiresReeval_ = new bool[numStorageEntries_];
    deltaMat_requiresReeval_ = new bool[numStorageEntries_];
    deltaMatStrain_requiresReeval_ = new bool[numStorageEntries_];

		rotatedCouplingTensor_ = new Matrix<Double>[numStorageEntries_];
		rotatedCouplingTensor_requiresReeval_ = new bool[numStorageEntries_];
    lastUsedTimeLevelForRotation_ = Vector<int>(numStorageEntries_);

    takeEstimatedSlope_ = new bool[numStorageEntries_];
    FPMaterialTensor_ = new Matrix<Double>[numStorageEntries_];
    maxSlopeLocal_ = Vector<Double>(numStorageEntries_);
    minSlopeLocal_ = Vector<Double>(numStorageEntries_);
    muFPH_ = Vector<Double>(numStorageEntries_);

		//TODO: coupling tensor einlesen und rotatedCouplingTensor entsprechend initialisieren

		Vector<Double> slope = Vector<Double>(dim_);
		slope.Init();

    rev_mat_fac_MATRIX_ = Matrix<Double>(dim_,dim_);
    // use same tensor as for inversion, rhs evaluation and so on
//    rev_mat_fac_MATRIX_.Init();
//
//    for(UInt i = 0; i < dim_; i++){
//      rev_mat_fac_MATRIX_[i][i] = rev_mat_fac_;
//    }

    rev_mat_fac_MATRIX_ = eps_nu_base_;
//    std::cout << "nu_0 Tensor: " << rev_mat_fac_MATRIX_.ToString() << std::endl;

    FPMaterialTensorSet_ = 0;

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
      E_B_[k] = initial_E_B;
      E_H_[k] = initial_E_H;

      if(strainRequired_){
        P_J_forStrains_[k] = P_J_forStrains;
        Si_[k] = initial_Si;
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
      }

      deltaMat_[k] = Matrix<Double>(dim_,dim_);
      deltaMat_[k].Init();// = POL_initialTensor_;

      deltaMatPrev_[k] = Matrix<Double>(dim_,dim_);
      deltaMatPrev_[k].Init();// = POL_initialTensor_;

		if (POL_initialInput_.NormL2() > 1e-16) {
        // TODO: compute first deltaMat or start without it
				// create a first deltaMatrix
				//CreateDeltaMatrix(POL_initialInput_, POL_initialOutput_[k], deltaMat_[k], evalMethodName, k, intoSat, outofSat, satToSat, POL_initialInput_);
			}
      // for inversion we need mu in magnetic case (not nu!)
      // TODO: init matrixForInversion with 0 (and let CoefFundtionHystMatrial etc initialize with actual values
      //        > this has to be done in those additinal functions as this function here does not have the required tensors for the coupled case
      eps_mu_local_[k] = eps_mu_base_;//.Init();
      epsInv_nu_local_[k] = Matrix<Double>(dim_,dim_);
      epsInv_nu_local_[k].Init();
      eps_mu_local_Set_[k] = false;

      FPMaterialTensor_[k] = rev_mat_fac_MATRIX_;//.Init();
		}

    InitLPMMaps();
    storageInitialized_ = true;

#ifdef USE_OPENMP
    UInt numT = CFS_NUM_THREADS;
    if(numT > 1){
      LOG_DBG(coeffunctionhyst_main) << "Using OPENMP with " << numT << " threads; initialize storage directly";
      InitSpecificStorage(0);
      InitSpecificStorage(1);
      InitSpecificStorage(2);
      InitSpecificStorage(3);
    } else {
      LOG_DBG(coeffunctionhyst_main) << "Using OPENMP with only " << numT << " thread; initialize storage on demand";
    }
#else
    LOG_DBG(coeffunctionhyst_main) << "Running single-threaded; initialize storage on demand";
#endif
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
    UInt storageIdx,operatorIdx,idxElem;

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
            PreprocessLPM(lpm, actualLPM, operatorIdx, storageIdx,idxElem);
						//std::cout << "Integration point NR: " << i << std::endl;
						//std::cout << "On boundary? " << onBoundary << std::endl;

            // store ONLY the ORIGINAL lpm
            // even thought this one needs to be postprocessed again
            // > reason for not storing the actual lpm: PreprocessLPM has to be called either way to get operatorIdx
            LOG_TRACE(coeffunctionhyst_main) << "Insert pair "<<storageIdx << " / " << lpm.lp.coord.ToString() << " into allLPMmap_";
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
        PreprocessLPM(lpm, actualLPM, operatorIdx, storageIdx,idxElem,forceMidpoint);
        LOG_TRACE(coeffunctionhyst_main) << "Insert pair "<<storageIdx << " / " << lpm.lp.coord.ToString() << " into allLPMmap_";
        allLPMmap_[storageIdx]=lpm;
        LOG_TRACE(coeffunctionhyst_main) << "Insert pair "<<storageIdx << " / " << lpm.lp.coord.ToString() << " into midpointLPMmap_";
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

    LOG_TRACE(coeffunctionhyst_main) << "Mapping between operatorIndices and storageIndices";
    std::map< UInt, std::map< UInt, bool > >::iterator outerIt;
    std::map< UInt, bool >::iterator innerIt;
    for(outerIt = operatorToStorage_.begin(); outerIt != operatorToStorage_.end(); outerIt++){
      LOG_TRACE(coeffunctionhyst_main) << "Operator Idx - Storage Idx - isMidpoint";

      for(innerIt = (outerIt->second).begin(); innerIt != (outerIt->second).end(); innerIt++){
        LOG_TRACE(coeffunctionhyst_main) << outerIt->first << " - " << innerIt->first << " - " << innerIt->second;
      }
    }

  }

  /*
   * Extension 24.6.2020 - return index of element, too
   * - can be used e.g., in the estimation of min/max slopes for local fixpoint method; here a new flag is
   * introduced that triggers slope estimate only at the midpoint (and midpointidx = idxelem!); 
   * 
   */
	bool CoefFunctionHyst::PreprocessLPM(const LocPointMapped& lpmInput,
          LocPointMapped& lpmOutput, UInt& operatorIdx, UInt& storageIdx, UInt& idxElem, bool forceMidpoint) {
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
//      std::cout << "Evaluate only at midpoint" << std::endl;
			/*
       * get solution at midpoint of element (elemSolution)
       */
			LocPoint lp = Elem::shapes[el->type].midPointCoord;
			shared_ptr<ElemShapeMap> esm = lpmInput.shapeMap;

			lpmOutput.Set(lp, esm, 0.0);

		} else {
//      std::cout << "Evaluate at integration point" << std::endl;
			/*
       * get solution at actual integration point
       */
			lpmOutput = lpmInput;
		}

		idxElem = globalElem2Local_[el->elemNum];

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

  Matrix<Double> CoefFunctionHyst::EstimateSlope_Jacobian(const LocPointMapped& Originallpm, bool includeStrains, bool useAbs, std::string implementationVersion, bool forceMidpoint){

    Matrix<Double> estimatedSlope;

    UInt numCols = dim_;
    UInt numRows = dim_;
    estimatedSlope = Matrix<Double>(numRows,numCols);

    UInt operatorIdx, storageIdx,idxElem;
		LocPointMapped actualLPM;

		bool onBoundary = PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx,idxElem,forceMidpoint);
	  // deltaMat can only be computed on volume elements, not on boundaries
		//bool onBoundary = false;
		if(onBoundary == true){
			EXCEPTION("Slope estimation only for volume elements.");
		}

    UInt currentTimelevel = 0;
    Vector<Double> E_H_current = GetPrecomputedInputToHysteresisOperator(Originallpm, currentTimelevel,forceMidpoint);
    Vector<Double> E_B_fromSystem = RetrieveLPMSolution(actualLPM, storageIdx, currentTimelevel, onBoundary);
    Vector<Double> P_J_current = GetPrecomputedOutputOfHysteresisOperator(Originallpm,currentTimelevel,false,forceMidpoint);

    // NOTE: if inversion fails, P_J_ + E_H_ will not return E_B_ and thus Jacobian will have
    // negative entries on diagonal (which should not be the case)
    // Two solutions:
    // a) FIX INVERSION! (the only actual solution)
    // b) use E_B_cur = P_J_ + E_H_*mu instead of actaul E_B_current (workaround)

    Vector<Double> E_B_current = Vector<Double>(dim_);
    if(PDEName_ == "Electromagnetics"){
      // Magnetics: reconstruct B = mu*H + J
      E_B_current.Init();
      eps_mu_local_[storageIdx].Mult(E_H_current,E_B_current);
      E_B_current.Add(1.0,P_J_current);
    } else {
      E_B_current = E_H_current;
    }

//    std::cout << "E_H_current: " << E_H_current.ToString() << std::endl;
//    std::cout << "E_B_current: " << E_B_current.ToString() << std::endl;
//    std::cout << "P_J_current: " << P_J_current.ToString() << std::endl;
//    std::cout << "eps_mu_local_[storageIdx]: " << eps_mu_local_[storageIdx].ToString() << std::endl;
//
    bool miniOutput = false;
//    if(storageIdx == 1){
//      miniOutput = true;
//    }
    if(miniOutput){
      std::cout << "Estimate Slope using FD-Jacobian" << std::endl;
      std::cout << "currentTimelevel = " << currentTimelevel << std::endl;
      std::cout << "P_J_current = " << P_J_current.ToString() << std::endl;
      std::cout << "E_B_current (computed): " << E_B_current.ToString() << std::endl;
      std::cout << "E_H_current: " << E_H_current.ToString() << std::endl;
      if(estimatorSet_){
        std::cout << "Use estimator!" << std::endl;
        std::cout << "E_B_estimatorStep_[storageIdx] = " << E_B_estimatorStep_[storageIdx].ToString() << std::endl;
        std::cout << "E_H_estimatorStep_[storageIdx] = " << E_H_estimatorStep_[storageIdx].ToString() << std::endl;
      }
    }

    //
    /*
     * Approximate Jacobian J = [J_ik] of F using Finite Differences, i.e.
     *
     *  J_ik = [F_i(x + dx*e_K) - F_i(x)]/[dx]
     *
     * Remark:
     * According to "Iterative Methods for Linear and Nonlinear Equations" - C.T.Kelly - Siam 1995
     *  dx = sqrt(machinePrecision)*|x|
     * (see page 80)
     *
     * if |x| = 0: dx = sqrt(machinePrecision)
     *
     */
    Double scaling = 1e-7; // approx sqrt(1e-15)
    //Double steppingDistance = std::max(1e-6*POL_operatorParams_.inputSat_,1e-3*E_H_current.NormL2());
    Double steppingDistance = std::max(scaling,scaling*E_H_current.NormL2());

    Vector<Double> E_B_diff = Vector<Double>(dim_);
    Vector<Double> P_J_diff = Vector<Double>(dim_);
    Vector<Double> E_B_stepping = Vector<Double>(dim_);
    Vector<Double> E_H_stepping = Vector<Double>(dim_);
    Vector<Double> P_J_stepping;

    bool forceMemoryLock = true;
    bool forceMemoryWrite = false;

    Double sign = 1.0;

    if(miniOutput){
      if(estimatorSet_){
        std::cout << "Compute FD Jacobian with signs based on difference between estimated E_H_new and E_H" << std::endl;
      }
    }
    for(UInt i = 0; i < dim_; i++){

      /*
       * determine in which direction to step
       * a) step further into the current directtion
       * b) step backwards
       *
       * default: step forwards, i.e. use same sign = sgn(E_H_current[i])
       * better (if estimator step was done; see SolveStepHyst):
       *  step towards estimated value, i.e. sign = sgn(E_B_estimatorStep_[i] - E_B_current_[i])
       * Note: take E_B_ instead of E_H_ as these values come directly from the system
       */
      if(estimatorSet_){
        if(E_B_estimatorStep_[storageIdx][i] < E_B_current[i]){
          sign = -1.0;
        } else {
          sign = 1.0;
        }
        if(miniOutput){
          std::cout << "sign["<<i<<"]: "<<sign<<std::endl;
        }
      } else {
        if(E_H_current[i] < 0){
          sign = -1.0;
        } else {
          sign = 1.0;
        }
      }

//      std::cout << "Direction: " << i << std::endl;
      E_H_stepping = E_H_current;
      E_H_stepping[i] += sign*steppingDistance;

      P_J_stepping = CalcOutputOfHysteresisOperator(E_H_stepping,operatorIdx,storageIdx,forceMemoryLock,forceMemoryWrite,false);

      if(PDEName_ == "Electromagnetics"){
        // Magnetics: reconstruct B = mu*H + J
        E_B_stepping.Init();
        eps_mu_local_[storageIdx].Mult(E_H_stepping,E_B_stepping);
        E_B_stepping.Add(1.0,P_J_stepping);
      } else {
        // Electrostatics: E = E
        E_B_stepping = E_H_stepping;
      }

//      std::cout << "E_H_stepping: " << E_H_stepping.ToString() << std::endl;
//      std::cout << "E_B_stepping: " << E_B_stepping.ToString() << std::endl;
//      std::cout << "P_J_stepping: " << P_J_stepping.ToString() << std::endl;
//
//
      E_B_diff.Init();
      E_B_diff.Add(1.0,E_B_stepping,-1.0,E_B_current);

      P_J_diff.Init();
      P_J_diff.Add(1.0,P_J_stepping,-1.0,P_J_current);

      for(UInt j = 0; j < dim_; j++){
        estimatedSlope[j][i] = P_J_diff[j]/E_B_diff[i];
      }
    }

    if(miniOutput){
      if(estimatorSet_){
        LOG_DBG(coeffunctionhyst_main) << "Compute FD Jacobian with signs based on current E_H";
        // perform once more without estimator for comparison
        Matrix<Double> estimatedSlope2 = Matrix<Double>(numRows,numCols);
        for(UInt i = 0; i < dim_; i++){

          if(E_H_current[i] < 0){
            sign = -1.0;
          } else {
            sign = 1.0;
          }
          LOG_DBG2(coeffunctionhyst_main) << "sign["<<i<<"]: "<<sign;

          //      std::cout << "Direction: " << i << std::endl;
          E_H_stepping = E_H_current;
          E_H_stepping[i] += sign*steppingDistance;

          P_J_stepping = CalcOutputOfHysteresisOperator(E_H_stepping,operatorIdx,storageIdx,forceMemoryLock,forceMemoryWrite,false);

          if(PDEName_ == "Electromagnetics"){
            // Magnetics: reconstruct B = mu*H + J
            E_B_stepping.Init();
            eps_mu_local_[storageIdx].Mult(E_H_stepping,E_B_stepping);
            E_B_stepping.Add(1.0,P_J_stepping);
          } else {
            // Electrostatics: E = E
            E_B_stepping = E_H_stepping;
          }

          LOG_DBG(coeffunctionhyst_main) << "steppingDistance: " << steppingDistance;
          LOG_DBG(coeffunctionhyst_main) << "E_H_stepping: " << E_H_stepping.ToString();
          LOG_DBG(coeffunctionhyst_main) << "E_B_stepping: " << E_B_stepping.ToString();
          LOG_DBG(coeffunctionhyst_main) << "P_J_stepping: " << P_J_stepping.ToString();


          E_B_diff.Init();
          E_B_diff.Add(1.0,E_B_stepping,-1.0,E_B_current);

          P_J_diff.Init();
          P_J_diff.Add(1.0,P_J_stepping,-1.0,P_J_current);

          LOG_DBG(coeffunctionhyst_main) << "E_B_diff: " << E_B_diff.ToString();
          LOG_DBG(coeffunctionhyst_main) << "P_J_diff: " << P_J_diff.ToString();

          for(UInt j = 0; j < dim_; j++){
            estimatedSlope2[j][i] = P_J_diff[j]/E_B_diff[i];
          }
        }
        LOG_DBG2(coeffunctionhyst_main) << "Estimated Jacobian without estimator: " << estimatedSlope2.ToString();
        LOG_DBG2(coeffunctionhyst_main) << "> estimator is only used for estimating in which direction to perform FD!";
      }
      LOG_DBG(coeffunctionhyst_main) << "Estimated Jacobian: " << estimatedSlope.ToString();

    }


//    if(operatorIdx == 0){
//      std::cout << "Computed slope for operator " << operatorIdx << " / storage " << storageIdx << " : " << estimatedSlope.ToString() << std::endl;
//    }
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

    bool compareJacobians = !true;
    if(compareJacobians && miniOutput){
      std::cout << "---> Compare different Jacobian approximations <---" << std::endl;
      std::cout << "Common input parameter: " << std::endl;
      std::cout << "E_H: " << E_H_current.ToString() << std::endl;
      std::cout << "eps_mu: " << eps_mu_local_[storageIdx].ToString() << std::endl;

      Matrix<Double> jac = JacobianApproximationOfMaterialRelation(E_H_current, eps_mu_local_[storageIdx], operatorIdx, 1);
      jac = JacobianApproximationOfMaterialRelation(E_H_current, eps_mu_local_[storageIdx], operatorIdx, 2);
      jac = JacobianApproximationOfMaterialRelation(E_H_current, eps_mu_local_[storageIdx], operatorIdx, 3);
    }

    return estimatedSlope;
  }

  Matrix<Double> CoefFunctionHyst::EstimateSlope(const LocPointMapped& Originallpm, bool includeStrains, bool useAbs, std::string estimationType,
          std::string implementationVersion, Double steppingLength, Double scaling, bool forceMidpoint){

    if(estimationType == "Jacobian"){
      return EstimateSlope_Jacobian(Originallpm, includeStrains, useAbs, implementationVersion, forceMidpoint);//, steppingLength, scaling);
    } else if(estimationType == "DeltaMat"){
      return EstimateSlope_DeltaMat(Originallpm, includeStrains, useAbs, implementationVersion, forceMidpoint);//, steppingLength, scaling);
    } else {
      EXCEPTION("Slope estimation not known");
    }
  }


  Matrix<Double> CoefFunctionHyst::EstimateSlope_DeltaMat(const LocPointMapped& Originallpm, bool includeStrains, bool useAbs, std::string implementationVersion, bool forceMidpoint){

    Matrix<Double> estimatedSlope;

    UInt numCols = dim_;
    UInt numRows = dim_;
    estimatedSlope = Matrix<Double>(numRows,numCols);

    UInt operatorIdx, storageIdx,idxElem;
		LocPointMapped actualLPM;

		bool onBoundary = PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx,idxElem,forceMidpoint);
	  // deltaMat can only be computed on volume elements, not on boundaries
		//bool onBoundary = false;
		if(onBoundary == true){
			EXCEPTION("Slope estimation only for volume elements.");
		}

    UInt currentTimelevel = 0;
    Vector<Double> E_H_current = GetPrecomputedInputToHysteresisOperator(Originallpm, currentTimelevel,forceMidpoint);
    Vector<Double> E_B_current = RetrieveLPMSolution(actualLPM, storageIdx, currentTimelevel, onBoundary);
    Vector<Double> P_J_current = GetPrecomputedOutputOfHysteresisOperator(Originallpm,currentTimelevel,false,forceMidpoint);

    /*
     * Approximate Jacobian J = [J_ik] of F using DeltaMat like approach, i.e.
     *
     *      J_ii = [F*_i - F_i(x)]/[x* - x]
     *      J_ij = 0      if i neq j
     *
     * Remark:
     * According to "Iterative Methods for Linear and Nonlinear Equations" - C.T.Kelly - Siam 1995
     *  dx = sqrt(machinePrecision)*|x|
     * (see page 80)
     *
     * if |x| = 0: dx = sqrt(machinePrecision)
     *
     */
    Double scaling = 1e-7; // approx sqrt(1e-15)
    //Double steppingDistance = std::max(1e-6*POL_operatorParams_.inputSat_,1e-3*E_H_current.NormL2());
    Double steppingDistance = std::max(scaling,scaling*E_H_current.NormL2());

    /*
     * if estimator was computed, step towards it;
     * otherwise step into 1,1,1 direction
     */
    Vector<Double> step = Vector<Double>(dim_);

    if(estimatorSet_){
        step.Init();
        step.Add(steppingDistance,E_H_estimatorStep_[storageIdx],-steppingDistance,E_H_current);
      } else {
        step.Init(steppingDistance);
      }

    Vector<Double> E_H_stepping = Vector<Double>(dim_);
    Vector<Double> E_B_stepping = Vector<Double>(dim_);
    Vector<Double> P_J_stepping;

    E_H_stepping.Init();
    E_H_stepping.Add(1.0,E_H_current,1.0,step);

    // P_J has to be computed using E_H
    bool forceMemoryLock = true;
    bool forceMemoryWrite = false;
    P_J_stepping = CalcOutputOfHysteresisOperator(E_H_stepping,operatorIdx,storageIdx,forceMemoryLock,forceMemoryWrite,false);

    if(PDEName_ == "Electromagnetics"){
      // Magnetics: reconstruct B = mu*H + J
      E_B_stepping.Init();
      eps_mu_local_[storageIdx].Mult(E_H_stepping,E_B_stepping);
      E_B_stepping.Add(1.0,P_J_stepping);
    } else {
      // Electrostatics: E = E
      E_B_stepping = E_H_stepping;
    }

    Vector<Double> E_B_diff = Vector<Double>(dim_);
    E_B_diff.Init();
    E_B_diff.Add(1.0,E_B_stepping,-1.0,E_B_current);

    Vector<Double> P_J_diff = Vector<Double>(dim_);
    P_J_diff.Init();
    P_J_diff.Add(1.0,P_J_stepping,-1.0,P_J_current);

    P_J_diff.ScalarMult(scaling);

    estimatedSlope = CalcDeltaMat(E_B_diff, P_J_diff, useAbs,implementationVersion);

    LOG_TRACE(coeffunctionhyst_main) << "Estimated Slope / Delta Mat (storageIdx="<<storageIdx << "): " << estimatedSlope.ToString();

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

  Matrix<Double> CoefFunctionHyst::DeltaMatApproximationOfMaterialRelation(Vector<Double>& E_H_diff, Vector<Double>& P_J_diff,
          Matrix<Double>& eps_mu, bool needsInversion){
  /*
   * Approximate slope of current material law (see JacobianApproximationOfMaterialRelation) using deltaMatrix approach,
   * i.e. instead of evaluating the Jacobian around a given input, we compute the differential quotient by dividing
   * deltaD/deltaE (in electrostatics) and dH/dB (in magnetics)
   *
   * This is basically the same what was done so far in ComputeTensor for the premittivity and the reluctivity tensors
   * except that we do not compute the deltaMatrix of deltaP/deltaE and deltaM/deltaB but instead directly compute
   * deltaD/deltaE and deltaH/deltaB
   */

    // 1. start by getting the DeltaMatrix > use current implementation; seems to do well
    // only difference: we use E_H_diff instead of E_B_diff, i.e. we divide by the input difference even in magnetics
    // this will give us dJ/dH instead of the previously computed dJ/dB
    // however, this function returns
    // dH/dB = inv(dB/dH) = inv(mu + dJ/dH)
    // > dJ/dH works perfect!
    Matrix<Double> deltaMat = CalcDeltaMat(E_H_diff, P_J_diff, false, "Division", -1.0);

    // 2. add material tensor
    // this gives us dD/dE and dB/dH respectively
    deltaMat.Add(1.0,eps_mu);

    if(needsInversion){
      // magnetics
      Matrix<Double> deltaMat_inv = Matrix<Double>(dim_,dim_);
      deltaMat.Invert(deltaMat_inv);
      return deltaMat_inv;
    } else {
      return deltaMat;
    }

  }

    Matrix<Double> CoefFunctionHyst::JacobianApproximationOfMaterialRelation(Vector<Double>& E_H, Matrix<Double>& eps_mu,
          UInt operatorIdx, UInt version, Double steppingSign){

    /*
     * Compute approximation of material relation
     *
     * Electrostatics:
     *  D = eps*E + P(E)
     * > Jacobian_DE_ik = dD_i/dE_k = eps_ik + dP_i(E)/dE_k
     * > approximate dP_i(E)/dE_k with a Finite Difference approximation
     *
     * Magnetics:
     *  B = mu*H + P(H)
     * > Jacobian_BH_ik = dB_i/dH_k = mu_ik + dP_i(H)/dH_k
     * > dP_i(H)/dH_k could be approximated with Finite Differences
     * BUT: we need H instead of B
     *  H = nu*(B - P(H))
     * > Jacobian_HB_ik = dH_i/dB_k = nu_ik - (nu*dP(H))_i/dB_k
     *                  = nu_ik - sum_i (nu_ij*dP(H)_j)/dB_k
     *  -- assuming nu const -- = nu_ik - sum_i nu_ij*dP(H)_j/dB_k
     *  -- assuming nu const and diagonal -- = = nu_ik - nu_ik*dP(H)_i/dB_k
     *
     * > Option 1 (current version, assuming nu const):
     *    1. take H_current, compute H_shifted = H_current + H_stepping
     *       take P(H_current), compute P(H_shifted)
     *    2. compute DeltaP = P(H_shifted) - P(H_current)
     *    3. compute DeltaB = mu*H_stepping + DeltaP
     *    4. approximate dP(H)_j/dB_k via
     *           dP(H)_j/dB_k = DeltaP_j/DeltaB_k
     *    5. compute Jacobian_HB_ik = nu_ik - sum_i nu_ij*dP(H)_j/dB_k
     *
     *    > disadvantage: DeltaP fits to H_stepping but is divided by B_stepping
     *
     * > Option 2 (assuming nu const):
     *    1. take B_current, compute B_shifted = B_current + B_stepping
     *       compute P(B_current) and P(B_shifted) using two step approach (could be combined, too)
     *        a) compute input to hyst operator, i.e. find H such that B = mu*H + P(H)
     *        b) evaluate P(H)
     *    2. compute DeltaP = P(B_shifted) - P(B_current)
     *    3. approximate dP(H)_j/dB_k via
     *           dP(H)_j/dB_k = DeltaP_j/B_stepping_k
     *    4. compute Jacobian_HB_ik = nu_ik - sum_i nu_ij*dP(H)_j/dB_k
     *
     *    > disadvantage: Input inversion 1a) very costly AND the most critical part at all; if inversion fails,
     *                      the Jacobian will be completely off!
     *
     * > Option 3 (assuming nu const): DEFAULT
     *    1. take H_current, compute H_shifted = H_current + H_stepping
     *       take P(H_current), compute P(H_shifted)
     *    2. compute DeltaP = P(H_shifted) - P(H_current)
     *    3. approximate dP(H)_i/dH_k via
     *          dP(H)_i/dH_k = DeltaP_i/DeltaH_k
     *    4. compute Jacobian_BH_ik = mu_ik + dP_i(H)/dH_k
     *    5. compute Jacobian_HB = inverse(Jacobian_BH)
     *
     *    > disadvantage: no idea if this estimate of slope could be correct
     *                    > tests showed, that Option 3 and Option 2 lead to nearly the same results (if option 2 does
     *                        not fail); Option 1 has slightly different values on off-diagonal
     *                    > Option 3 seems to be the best working version!
     */
      
    if( (InversionParams_.inversionMethod == LOCAL_NOTIMPLEMENTED) && (inversionSet_) ){
      // NOTE: in electrostatics, the inversion parameter are not se as we need no local inversion!
      // in case of the everett based inversion, we are no longer allowed to call computeValue_vec as
      // the retrieved input from computeValue_vec will not return the original input (the inverse of a sum is
      // not the sum of the inverted terms!); instead we compute the output by subtracting mu*E_H from 
      // B; problem: we do neither know B nor P here!
      EXCEPTION("No prober inversion method defined (InversionParams_.inversionMethod == LOCAL_NOTIMPLEMENTED), although inversionSet_ == True!");
    }
      
    Double scaling = 1e-7; // sqrt(double precision)
    Double steppingDistance;
    Vector<Double> stepping = Vector<Double>(dim_);
    Vector<Double> P = Vector<Double>(dim_);
    Vector<Double> Pshifted = Vector<Double>(dim_);
    Matrix<Double> Jacobian = Matrix<Double>(dim_,dim_);

    int successFlag = 0;
    bool debugOutput = false;
    bool overwriteMemory = false;

    /*
     * Electrostatics
     */
    if(version == 0){
      steppingDistance = steppingSign*std::max(scaling,scaling*E_H.NormL2());
      P = hyst_->computeValue_vec(E_H, operatorIdx, overwriteMemory, debugOutput, successFlag);

      for(UInt k = 0; k < dim_; k++){
        stepping = E_H;
        stepping[k] += steppingDistance;

        Pshifted = hyst_->computeValue_vec(stepping, operatorIdx, overwriteMemory, debugOutput, successFlag);

        for(UInt i = 0; i < dim_; i++){
          Jacobian[i][k] = (Pshifted[i] - P[i])/steppingDistance;
        }
      }
      Jacobian.Add(1.0,eps_mu);
    }
    /*
     * Magnetics - Option 1
     */
    else if(version == 1){
//      std::cout << "Version 1 --- Compute dP/dB using H-stepping" << std::endl;

      Vector<Double> B = Vector<Double>(dim_);
      Vector<Double> Bshifted = Vector<Double>(dim_);
      Matrix<Double> nu = Matrix<Double>(dim_,dim_);
      Matrix<Double> tmp = Matrix<Double>(dim_,dim_);

      eps_mu.Invert(nu);

      // H-stepping
      steppingDistance = steppingSign*std::max(scaling,scaling*E_H.NormL2());
      P = hyst_->computeValue_vec(E_H, operatorIdx, overwriteMemory, debugOutput, successFlag);

      eps_mu.Mult(E_H,B);
      B.Add(1.0,P);

      for(UInt k = 0; k < dim_; k++){
        stepping = E_H;
        stepping[k] += steppingDistance;

        Pshifted = hyst_->computeValue_vec(stepping, operatorIdx, overwriteMemory, debugOutput, successFlag);

        eps_mu.Mult(stepping,Bshifted);
        Bshifted.Add(1.0,Pshifted);

        for(UInt i = 0; i < dim_; i++){
          tmp[i][k] = -(Pshifted[i] - P[i])/(Bshifted[k] - B[k]);
        }
        // tmp = 1.0 - dP/dB
        tmp[k][k] += 1.0;
      }
      // scale with nu
      nu.Mult(tmp,Jacobian);
    }
    /*
     * Magnetics - Option 2
     */
    else if(version == 2){
//      std::cout << "Version 2 --- Compute dP/dB using B-stepping" << std::endl;

      Vector<Double> B = Vector<Double>(dim_);
      Vector<Double> Hshifted = Vector<Double>(dim_);
      Vector<Double> Bshifted = Vector<Double>(dim_);
      Matrix<Double> nu = Matrix<Double>(dim_,dim_);
      Matrix<Double> tmp = Matrix<Double>(dim_,dim_);

      eps_mu.Invert(nu);

      P = hyst_->computeValue_vec(E_H, operatorIdx, overwriteMemory, debugOutput, successFlag);

      eps_mu.Mult(E_H,B);
      B.Add(1.0,P);

      LOG_DBG(coeffunctionhyst_main) << "P = " << P.ToString();
      LOG_DBG(coeffunctionhyst_main) << "B = " << B.ToString();

      // B-stepping
      steppingDistance = steppingSign*std::max(scaling,scaling*B.NormL2());

      for(UInt k = 0; k < dim_; k++){
        LOG_DBG2(coeffunctionhyst_main) << "k = " << k;

        stepping = B;
        stepping[k] += steppingDistance;

        Hshifted = hyst_->computeInput_vec(stepping,operatorIdx, eps_mu, POL_operatorParams_.fieldsAlignedAboveSat_,POL_operatorParams_.hystOutputRestrictedToSat_,successFlag);
        LOG_DBG2(coeffunctionhyst_main) << "Inversion success? " << successFlag;
        
        Pshifted = hyst_->computeValue_vec(Hshifted, operatorIdx, overwriteMemory, debugOutput, successFlag);

        eps_mu.Mult(Hshifted,Bshifted);
        Bshifted.Add(1.0,Pshifted);

        LOG_DBG2(coeffunctionhyst_main) << "Hshifted: " << Hshifted.ToString();
        LOG_DBG2(coeffunctionhyst_main) << "Pshifted: " << Pshifted.ToString();
        LOG_DBG2(coeffunctionhyst_main) << "Bshifted: " << Bshifted.ToString();
        LOG_DBG2(coeffunctionhyst_main) << "B + stepping*e_k = " << stepping.ToString();

        for(UInt i = 0; i < dim_; i++){
          tmp[i][k] = -(Pshifted[i] - P[i])/steppingDistance;
        }
        // tmp = 1.0 - dP/dB
        tmp[k][k] += 1.0;
      }

      LOG_DBG(coeffunctionhyst_main) << "Pre scaling: " << tmp.ToString();
      // scale with nu
      nu.Mult(tmp,Jacobian);
    }
    /*
     * Magnetics - Option 3
     */
    else if((version == 3)||(version == 4)){
//      std::cout << "Version 3 --- Compute dP/dH using H-stepping, then invert" << std::endl;
      bool output = false;
//      if(operatorIdx == 1){
//        output = true;
//      }
//

      Matrix<Double> tmp = Matrix<Double>(dim_,dim_);

      steppingDistance = steppingSign*std::max(scaling,scaling*E_H.NormL2());

      if(output){
        std::cout << "steppingDistance: " << steppingDistance << std::endl;
      }

      P = hyst_->computeValue_vec(E_H, operatorIdx, overwriteMemory, debugOutput, successFlag);

      if(output){
        std::cout << "P_base: " << P.ToString() << std::endl;
        std::cout << "E_H_base: " << E_H.ToString() << std::endl;
      }

      for(UInt k = 0; k < dim_; k++){
        stepping = E_H;
        // test: step towards larger values
//        if(stepping[k] > 0){
          stepping[k] += steppingDistance;
//        } else {
//          stepping[k] -= steppingDistance;
//        }
        Pshifted = hyst_->computeValue_vec(stepping, operatorIdx, overwriteMemory, debugOutput, successFlag);

        if(output){
          std::cout << "Pshifted ["<<k<<"]: " << Pshifted.ToString() << std::endl;
          std::cout << "E_H_shifted ["<<k<<"]: " << stepping.ToString() << std::endl;
        }

        // compute dP/dH
        for(UInt i = 0; i < dim_; i++){
          tmp[i][k] = (Pshifted[i] - P[i])/steppingDistance;
        }
      }
      if(output){
        std::cout << "dB/dH (ohne eps_mu): " << tmp.ToString() << std::endl;
      }

      // dB/dH
      tmp.Add(1.0,eps_mu);

      if(output){
        std::cout << "dB/dH (mit eps_mu): " << tmp.ToString() << std::endl;
      }

      if(version == 3){
        // reluctivity
        // dH/dB
        tmp.Invert(Jacobian);
        if(output){
          std::cout << "dH/dB: " << Jacobian.ToString() << std::endl;
        }
      } else {
        // permeability
        Jacobian = tmp;
        if(output){
          std::cout << "dB/dH: " << Jacobian.ToString() << std::endl;
        }
      }
    }
    else {
      EXCEPTION("Approximate Jacobian with flag version = "<<version<<" not allowed");
    }

//    std::cout << "> Jacobian approximation of material law: " << std::endl;
//    std::cout << Jacobian.ToString() << std::endl;

    return Jacobian;
  }

  Matrix<Double> CoefFunctionHyst::GetMaterialRelation(const LocPointMapped& Originallpm, int timelevel_new, int timelevel_old,
          Matrix<Double> eps_mu, std::string tensorName, bool prohibitUseOfStorage){

    // TODO: Implement axi case!
    if(tensorType_ == AXI){
      EXCEPTION("GetMaterialRelation only implemented for 2d plane and full 3d setups");
    }

    // define return matrix and get info about current evaluation point (originallpm)
    Matrix<Double> deltaMat;
    UInt operatorIdx, storageIdx,idxElem;
    LocPointMapped actualLPM;
    bool forceMidpoint = evalJacAtMidpointOnly_;
//    if(forceMidpoint) std::cout << "GetMaterialRelation - forceMidpoint" << std::endl;
    bool onBoundary = PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx,idxElem,forceMidpoint);
    // deltaMat can only be computed on volume elements, not on boundaries
    if(onBoundary == true){
      EXCEPTION("GetMaterialRelation not defined on boundary");
    }

    // check for special cases
    // a) timelevel_new = timelevel_old > no delta
    //      > reuse current state
    // b) timelevel_old = 2
    //      > reuse current state
    // c) timelevel_old = 3
    //      > instead of deltaMat, compute FiniteDifference-Jacobian arround the current solution
    bool reuseOldState = false;
    bool useJacobian = false;

    if(timelevel_old == 3){
      useJacobian = true;
    } else if((timelevel_new == timelevel_old)||(timelevel_old == 2)){
      /*
       * NEW: 11.10.2018 -> use this as trigger to relaod matrix
       */
      reuseOldState = true;
    }

    // prohibitUseOfStorage only for comparison with old version; without this flag,
    // the deltaMat approximation would already set deltaMat_requiresReeval_ to false and thus we would
    // just load that state without performing the actual new computation
    if( !prohibitUseOfStorage && ((deltaMat_requiresReeval_[storageIdx] == false)||(reuseOldState))){
      deltaMat = deltaMat_[storageIdx];
    } else {
      /*
       * Reluctivity/Permittivity
       */
      if(useJacobian){
        UInt version = 0; // for Permittivity
        Vector<Double> E_H_new = GetPrecomputedInputToHysteresisOperator(Originallpm, timelevel_new,forceMidpoint);
        if(tensorName == "Reluctivity"){
          version = 3; // works best so far
        } else if(tensorName == "Permeability"){
          version = 4;
        } 
        deltaMat = JacobianApproximationOfMaterialRelation(E_H_new, eps_mu, operatorIdx, version);
      } else {
        bool needsInversion = false;
        if(tensorName == "Reluctivity"){
          needsInversion = true;
        }

        Vector<Double> E_H_new = GetPrecomputedInputToHysteresisOperator(Originallpm, timelevel_new,forceMidpoint);
        Vector<Double> E_H_old = GetPrecomputedInputToHysteresisOperator(Originallpm, timelevel_old,forceMidpoint);

        Vector<Double> P_J_new = GetPrecomputedOutputOfHysteresisOperator(Originallpm,timelevel_new,false,forceMidpoint);
        Vector<Double> P_J_old = GetPrecomputedOutputOfHysteresisOperator(Originallpm,timelevel_old,false,forceMidpoint);

        Vector<Double> E_H_diff = Vector<Double>(dim_);
        E_H_diff.Add(1.0,E_H_new,-1.0,E_H_old);
        Vector<Double> P_J_diff = Vector<Double>(dim_);
        P_J_diff.Add(1.0,P_J_new,-1.0,P_J_old);

        deltaMat = DeltaMatApproximationOfMaterialRelation(E_H_diff, P_J_diff, eps_mu, needsInversion);

        bool compareToJacobian = false;
//        if(operatorIdx == 1){
//          compareToJacobian = true;
//        }

        if(compareToJacobian){
          Matrix<Double> jacobian = Matrix<Double>(dim_,dim_);
          jacobian = JacobianApproximationOfMaterialRelation(E_H_new, eps_mu, operatorIdx, 3);

          std::cout << "Timelevels old/new: " << timelevel_old << "/" << timelevel_new << std::endl;
          std::cout << "E_H_diff: " << E_H_diff.ToString() << std::endl;
          std::cout << "P_J_diff: " << P_J_diff.ToString() << std::endl;

          std::cout << "Computed deltaMat: " << std::endl;
          std::cout << deltaMat.ToString() << std::endl;
          std::cout << "Computed Jacobian: " << std::endl;
          std::cout << jacobian.ToString() << std::endl;
        }

      }
      // TODO: comment in, once everything goes over this function
      // kept out at the moment to check for compatibility
      if(!prohibitUseOfStorage){
        deltaMat_[storageIdx] = deltaMat;
        deltaMat_requiresReeval_[storageIdx] = false;
      }
    }

    return deltaMat;
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
      LOG_DBG(coeffunctionhyst_deltamat) << "E_B_diff_norm <= " << absTol;
      LOG_DBG(coeffunctionhyst_deltamat) << "> set deltaMat to 0";
      // variation of solution is very small
      // return zero matrix
//      std::cout << "set deltaMat to 0 as E_B_diff_norm is too small!" << std::endl;
      return deltaMat;
    }

    for(UInt i = 0; i < P_J_diff.GetSize(); i++){
      // if cuttingTol is negative (default), nothing happens
      // else: cut down nominato
      // note: s_diff is passed by value, so it is safe to edit vector directly
      LOG_DBG(coeffunctionhyst_deltamat) << "P_J_diff["<<i<<"]: " << P_J_diff[i];
      LOG_DBG(coeffunctionhyst_deltamat) << "cuttingTol: " << cuttingTol;
      LOG_DBG(coeffunctionhyst_deltamat) << "E_B_diff["<<i<<"]: " << E_B_diff[i];
      if(abs(P_J_diff[i]) < cuttingTol){
        LOG_DBG(coeffunctionhyst_deltamat) << "P_J_diff["<<i<<"] <= " << cuttingTol;
        LOG_DBG(coeffunctionhyst_deltamat) << "> set P_J_diff["<<i<<"] to zero";
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
  }

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
  }

  /*
   * Precompute/Preestimate Jacobian for Newtonlike schemes
   *
   * Idea:
   *  bundle different slope approximations together
   *
   * Notation:
   *  J = Jacobian of hyst operator
   *  F = hyst operator
   *  x = input to hyst operator
   *
   * Implemented approximations:
   *  1) FiniteDifferences-Jacobian
   *      J_ik = [F_i(x + dx*e_K) - F_i(x)]/[dx]
   *
   *  2) DeltaMat-Jacobian
   *      J_ii = [F*_i - F_i(x)]/[x* - x]
   *      J_ij = 0      if i neq j
   *  2a) towards estimator
   *        F* = F_estimated, x* = x_estimated
   *  2b) towards last timestep
   *        F* = F_lastTS, x* = x_lastTS
   *  2c) towards last iteration
   *        F* = F_lastIT, x* = x_lastIT
   *
   *  3) reuse previously computed value
   *      J = J_old
   */
//  bool CoefFunctionHyst::PrecomputeSlopes(const LocPointMapped& Originallpm, std::string implementationVersion,
//          bool strainVersion, bool forceReevaluation){
//
//    /*
//     * Initial checks
//     */
//    if(tensorType_ == AXI){
//      EXCEPTION("Slopes estimations only implemented for 2d plane and full 3d setups");
//    }
//    if((strainVersion == true)&&(implementationVersion == "FD-Jacobian")){
//      WARN("Currently, FD-Jacobian is not implemented for strains; taking DeltaMat-Jacobian instead");
//      implementationVersion = "DeltaMat-Jacobian";
//    }
//
//    /*
//     * Common variables
//     */
//    UInt numCols = dim_;
//    UInt numRows;
//
//    if(strainVersion){
//      if(dim_ == 3){
//        numRows = 6;
//      } else {
//        numRows = 3;
//      }
//    } else {
//      numRows = numCols;
//    }
//
//
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
//
//      Matrix<Double> deltaMat;
//
//      UInt numRows;
//
//      /*
//       * Initial step:
//       *  determine correct
//       */
//      UInt operatorIdx, storageIdx,idxElem;
//      LocPointMapped actualLPM;
//      bool onBoundary = PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx,idxElem);
//
//
//      Vector<Double> solution = Vector<Double>(dim_);
//      Vector<Double> input = Vector<Double>(dim_);
//      Vector<Double> output = Vector<Double>(dim_);
//
//      // RUN_evaluationPurpose_ only needed if we want to check
//      // a) if we want to evaluate at the midpoint or at every point
//      // b) if we want to set hyst memory or not
//      // > b is irrelevant here as we set the overwrite flag per hand
//      // > a would be irrelevant, too as we just go over whole map of LPM
//      //    but the function PreProcessLPM calls EvaluateAtMidPointOnly which
//      //    needs the correct evaluationpurpose
//      // > set purpose to 1 (assemble/eval)
//
//      bool onBoundary;
//      bool forceMemoryLock;
//      bool forceMemoryWrite;
//
//#pragma omp barrier
//      for(it = itStart; it != itEnd; it++){
//        UInt operatorIdx = it->first;
//        UInt cnt = 0;
//
//        std::map<UInt, bool>::iterator innerIt;
//        for(innerIt = it->second.begin(); innerIt != it->second.end(); innerIt++){
//          /*
//           * NOTE:
//           * for standard and full evaluation, the length of the inner list should
//           *    always be 1
//           *  > standard eval: 1 operator per element, evaluation only at midpoint
//           *  > full eval: 1 operato per int. point + 1 for midpoint, evaluation at corrsp. point
//           *
//           * for extended evaluation, the length should be larger than 1
//           *  > extended eval: 1 operator per element, eval at midpoint and integration points
//           *      (hyst memory may only be written at midpoint)
//           */
//          UInt storageIdx = innerIt->first;
//          bool isMidpoint = innerIt->second;
//          if(isMidpoint && (cnt != 0)){
//            EXCEPTION("Midpoint should always be the first storage entry assigned to each operator (if at all");
//          }
//          cnt++;
//
//          if(overwriteMemory){
//            forceMemoryLock = false;
//            forceMemoryWrite = true;
//          } else {
//            forceMemoryLock = true;
//            forceMemoryWrite = false;
//          }
//
//          if(!isMidpoint && (XML_EvaluationDepth_ == 2)){
//            // extended evaluation: only overwrite at midpoint!
//            forceMemoryLock = true;
//            forceMemoryWrite = false;
//          }
//
//          LocPointMapped curLPM = allLPMmap_[storageIdx];
//          LocPointMapped actualLPM;
//
//          UInt operatorIdxTMP, storageIdxTMP;
//          onBoundary = PreprocessLPM(curLPM, actualLPM, operatorIdxTMP, storageIdxTMP, isMidpoint);
//
//          if( (operatorIdxTMP != operatorIdx) || (storageIdxTMP != storageIdx)){
//            EXCEPTION("Storage or operator index does not fit");
//          }
//          // get input
//          input = RetrieveInputToHysteresisOperator(actualLPM, operatorIdx, storageIdx,onBoundary);
//          // RetrieveInputToHysteresisOperator will store the actual solution directly to E_B_
//
//          // compute additional hyst operator for strains first
//          // reason: CalcOutputOfHysteresisOperator sets E_H_ if last flag is false;
//          //  during the next call input will be compared to E_H_
//          // > if output is computed for standard hyst operator first, the call for the strain operator will
//          //  not evaluate anything as E_H_ == input
//          if(CouplingParams_.ownHystOperator_ && COUPLED_inXMLFile_ ){
//            output = CalcOutputOfHysteresisOperator(input, operatorIdx, storageIdx, forceMemoryLock, forceMemoryWrite, true);
//          }
//
//          // then compute output
//          output = CalcOutputOfHysteresisOperator(input, operatorIdx, storageIdx, forceMemoryLock, forceMemoryWrite, false);
//        }
//      }
//#pragma omp barrier
//    }
//
//
//
//
//
//    Matrix<Double> deltaMat;
//    UInt numCols = dim_;
//    UInt numRows;
//
//    /*
//     * Initial step:
//     *  determine correct
//     */
//    UInt operatorIdx, storageIdx;
//		LocPointMapped actualLPM;
//		bool onBoundary = PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx);
//
//    if(implementationVersion == "FD-Jacobian"){
//
//    } else if(implementationVersion == "ReuseOldState"){
//
//    } else {
//
//    }
//    if(implementationVersion == "DeltaMat-Estimator"){
//
//    }
//
//
//  }

  Matrix<Double> CoefFunctionHyst::GetDeltaMat(const LocPointMapped& Originallpm, int timelevel_new, int timelevel_old, bool useStrains, bool useAbs,
          std::string implementationVersion){

    LOG_DBG(coeffunctionhyst_deltamat) << "GetDeltaMat";
    //LOG_DBG(coeffunctionhyst_deltamat) << "Compute delta between timelevel " << timelevel_new << " and time level " << timelevel_old << "(-2 = deactivated; -1 = lastTS; 0 = current IT; 1 = lastIT)";

    // TODO: Implement axi case!
    if(tensorType_ == AXI){
      EXCEPTION("GetDeltaMat only implemented for 2d plane and full 3d setups");
    }
    //evalJacAtMidpointOnly_ = !true;
    bool forceMidpoint = evalJacAtMidpointOnly_;
//    if(forceMidpoint) std::cout << "GetDeltaMat - Force midpoint evaluation!" << std::endl;
    // define return matrix and get info about current evaluation point (originallpm)
    Matrix<Double> deltaMat;
    UInt operatorIdx, storageIdx,idxElem;
		LocPointMapped actualLPM;

		bool onBoundary = PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx,idxElem,forceMidpoint);
	  // deltaMat can only be computed on volume elements, not on boundaries
		if(onBoundary == true){
			EXCEPTION("DeltaMat not defined on boundary");
		}

    // for quick debugging
    bool miniOutput = false;
    bool compareDeltaMatrices = false;
//    if(storageIdx == 1){
//      miniOutput = true;
//    }

    // check for special cases
    // a) timelevel_new = timelevel_old > no delta
    //      > reuse current state
    // b) timelevel_old = 2
    //      > reuse current state
    // c) timelevel_old = 3
    //      > instead of deltaMat, compute FiniteDifference-Jacobian arround the current solution
    bool reuseOldState = false;
    bool useJacobian = false;

    if(timelevel_old == 3){
      useJacobian = true;
    } else if((timelevel_new == timelevel_old)||(timelevel_old == 2)){
      /*
       * NEW: 11.10.2018 -> use this as trigger to relaod matrix
       */
      reuseOldState = true;
    }

    if( ((deltaMatStrain_requiresReeval_[storageIdx] == false)||(reuseOldState)) && useStrains ){
      deltaMat = deltaMatStrain_[storageIdx];
    } else if( ((deltaMat_requiresReeval_[storageIdx] == false)||(reuseOldState)) && !useStrains ){
      if(miniOutput){
        std::cout << "Reuse previously computed deltaMat / slope approximation" << std::endl;
      }
      deltaMat = deltaMat_[storageIdx];
    } else {
      /*
       * NEW computation/evaluation required
       */
      UInt numCols = dim_;
      UInt numRows;

      if(useStrains){
        /*
         * Strain
         */

        // TODO: implement me!
        if(useJacobian){
          // use deltaMat instead!
          useJacobian = false;
//          EXCEPTION("FD-Jacobian currently not available for strains");
        }

        if(dim_ == 3){
          numRows = 6;
        } else {
          if(tensorType_ == AXI){
            numRows = 4;
          } else {
            numRows = 3;
          }
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

        Vector<Double> S_new = GetIrreversibleStrains(Originallpm,timelevel_new,forceMidpoint);
        Vector<Double> S_old = GetIrreversibleStrains(Originallpm,timelevel_old,forceMidpoint);
        Vector<Double> S_diff = S_new;
        S_diff -= S_old;

        // cutting does not work so well; remove in future; to disable > pass some negative value
        Double cuttingTol = -1.0; // 1e-5*S_old.NormL2();

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
        /*
         * Reluctivity/Permittivity
         */
        if(useJacobian){
          if(miniOutput){
            std::cout << "Estimate slope using a FD-Jacobian approximation around the current state" << std::endl;
          }
          deltaMat = EstimateSlope(Originallpm, ES_includeStrains_, ES_useAbs_, "Jacobian", ES_implementationVersion_, ES_steppingLength_, ES_scaling_,forceMidpoint);
        } else if(takeEstimatedSlope_[storageIdx]){
          takeEstimatedSlope_[storageIdx] = false;
          if(miniOutput){
            std::cout << "Estimate slope using a deltaMat approximation around the current state" << std::endl;
          }
          deltaMat = EstimateSlope(Originallpm, ES_includeStrains_, ES_useAbs_, "DeltaMat", ES_implementationVersion_, ES_steppingLength_, ES_scaling_,forceMidpoint);
        } else {

          if(compareDeltaMatrices && (storageIdx == 1)){
            LOG_DBG2(coeffunctionhyst_main) << "Compute deltaMat towards last TS and towards last IT for storage idx " << storageIdx;
            numRows = dim_;
            Matrix<Double> deltaMatTS = Matrix<Double>(numRows,numCols);
            Matrix<Double> deltaMatIT = Matrix<Double>(numRows,numCols);

            Vector<Double> E_B_cur = RetrieveLPMSolution(actualLPM, storageIdx, timelevel_new, onBoundary);
            Vector<Double> E_B_ts = RetrieveLPMSolution(actualLPM, storageIdx, -1, onBoundary);
            Vector<Double> E_B_it = RetrieveLPMSolution(actualLPM, storageIdx, 1, onBoundary);
            //
            //        std::cout << "E_B_ at " << std::endl;
            //        std::cout << "Current IT: " << E_B_cur.ToString() << std::endl;
            //        std::cout << "Previous IT: " << E_B_it.ToString() << std::endl;
            //        std::cout << "Previous TS: " << E_B_ts.ToString() << std::endl;

            Vector<Double> E_B_diff_it = E_B_cur;
            E_B_diff_it -= E_B_it;
            Vector<Double> E_B_diff_ts = E_B_cur;
            E_B_diff_ts -= E_B_ts;

            Vector<Double> P_J_cur = GetPrecomputedOutputOfHysteresisOperator(Originallpm,timelevel_new,false,forceMidpoint);
            Vector<Double> P_J_ts = GetPrecomputedOutputOfHysteresisOperator(Originallpm,-1,false,forceMidpoint);
            Vector<Double> P_J_it = GetPrecomputedOutputOfHysteresisOperator(Originallpm,1,false,forceMidpoint);

            //        std::cout << "P_J__ at " << std::endl;
            //        std::cout << "Current IT: " << P_J_cur.ToString() << std::endl;
            //        std::cout << "Previous IT: " << P_J_it.ToString() << std::endl;
            //        std::cout << "Previous TS: " << P_J_ts.ToString() << std::endl;
            //
            Vector<Double> P_J_diff_it = P_J_cur;
            P_J_diff_it -= P_J_it;
            Vector<Double> P_J_diff_ts = P_J_cur;
            P_J_diff_ts -= P_J_ts;

            Double cuttingTol = 1e-5*P_J_it.NormL2();
            deltaMatIT = CalcDeltaMat(E_B_diff_it, P_J_diff_it, useAbs, implementationVersion,cuttingTol);
            if(miniOutput){
              std::cout << "Computed deltaMat towards last IT: " << deltaMatIT.ToString() << std::endl;
            }
            cuttingTol = 1e-5*P_J_ts.NormL2();
            deltaMatTS = CalcDeltaMat(E_B_diff_ts, P_J_diff_ts, useAbs, implementationVersion,cuttingTol);
            if(miniOutput){
              std::cout << "Computed deltaMat towards last TS: " << deltaMatTS.ToString() << std::endl;
            }
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

          if(miniOutput){
            std::cout << "Calculated deltaMat approximation with" << std::endl;
            std::cout << "E_B_new (timelevel_new=" << timelevel_new <<"): " << E_B_new.ToString() << std::endl;
            std::cout << "E_B_old (timelevel_old=" << timelevel_old <<"): " << E_B_old.ToString() << std::endl;
          }

          LOG_DBG(coeffunctionhyst_deltamat) << "Old solution of system (solution at timelevel = " << timelevel_old << "): " << E_B_old.ToString();
          LOG_DBG(coeffunctionhyst_deltamat) << "Current solution of system (solution at timelevel = " << timelevel_new << "): " << E_B_new.ToString();

          //std::cout << "Old solution (solution at timelevel = " << timelevel_old << "): " << E_B_old.ToString() << std::endl;
          //std::cout << "Current solution (solution at timelevel = " << timelevel_new << "): " << E_B_new.ToString() << std::endl;

          Vector<Double> E_B_diff = E_B_new;
          E_B_diff -= E_B_old;

          Vector<Double> P_J_new = GetPrecomputedOutputOfHysteresisOperator(Originallpm,timelevel_new,false,forceMidpoint);
          Vector<Double> P_J_old = GetPrecomputedOutputOfHysteresisOperator(Originallpm,timelevel_old,false,forceMidpoint);

          LOG_DBG(coeffunctionhyst_deltamat) << "Old output of hyst operator: " << P_J_old.ToString();
          LOG_DBG(coeffunctionhyst_deltamat) << "Current output of hyst operator: " << P_J_new.ToString();

          //std::cout << "Old output of hystoperator (solution at timelevel = " << timelevel_old << "): " << P_J_old.ToString() << std::endl;
          //std::cout << "Current output of hystoperator (solution at timelevel = " << timelevel_new << "): " << P_J_new.ToString() << std::endl;

          Vector<Double> P_J_diff = P_J_new;
          P_J_diff -= P_J_old;

          // does not work well; remove in future
          Double cuttingTol = 1e-5*P_J_old.NormL2();

          Matrix<Double> deltaTMP = CalcDeltaMat(E_B_diff, P_J_diff, useAbs, implementationVersion,cuttingTol);
          if( (includeOldDelta_ > 0)&&(includeOldDelta_ <= 1)){
            deltaMat.Init();
            deltaMat.Add(1.0-includeOldDelta_,deltaTMP);
            deltaMat.Add(includeOldDelta_,deltaMatPrev_[storageIdx]);
          } else {
            deltaMat = deltaTMP;
          }
        }
        // TODO: move out, such that Jacobian will also be stored as deltaMat_
        // kept out at the moment to check for compatibility
        deltaMat_[storageIdx] = deltaMat;
        deltaMat_requiresReeval_[storageIdx] = false;
      }
    }

    if(miniOutput){
      std::cout << "deltaMat: " << deltaMat.ToString() << std::endl;
    }

    return deltaMat;
  }


//  Matrix<Double> CoefFunctionHyst::GetDeltaMat_preCleanup(const LocPointMapped& Originallpm, int timelevel_new, int timelevel_old, bool useStrains, bool useAbs,
//          std::string implementationVersion){
//
//    LOG_DBG(coeffunctionhyst_deltamat) << "GetDeltaMat";
//    //LOG_DBG(coeffunctionhyst_deltamat) << "Compute delta between timelevel " << timelevel_new << " and time level " << timelevel_old << "(-2 = deactivated; -1 = lastTS; 0 = current IT; 1 = lastIT)";
//
//    if(tensorType_ == AXI){
//      EXCEPTION("GetDeltaMat only implemented for 2d plane and full 3d setups");
//    }
//
//    Matrix<Double> deltaMat;
//    UInt numCols = dim_;
//    UInt numRows;
//
////    std::cout << "DeltaMat was requested with timelevel_new = " << timelevel_new << std::endl;
////    std::cout << "DeltaMat was requested with timelevel_old = " << timelevel_old << std::endl;
////
//    bool reuseOldState = false;
//    if((timelevel_new == timelevel_old)||(timelevel_old == 2)){
//      /*
//       * NEW: 11.10.2018 -> use this as trigger to relaod matrix
//       */
////      std::cout << "DeltaMat was requested with timelevel_new = " << timelevel_old << " > reuse old matrix" << std::endl;
//      reuseOldState = true;
//    }
//
//    UInt operatorIdx, storageIdx;
//		LocPointMapped actualLPM;
//
//		bool onBoundary = PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx);
//	  // deltaMat can only be computed on volume elements, not on boundaries
//		//bool onBoundary = false;
//		if(onBoundary == true){
//			EXCEPTION("DeltaMat not defined on boundary");
//		}
//
//    Matrix<Double> estimatedSlope;
//    // for testing
//    bool useJacobian = false;
//    bool miniOutput = false;
//
//    if(storageIdx == 1){
//      miniOutput = true;
//    }
//
//    if(timelevel_old == 3){
//      useJacobian = true;
//    }
//
//    if(useJacobian){
//       estimatedSlope = EstimateSlope(Originallpm, ES_includeStrains_, ES_useAbs_, "Jacobian", ES_implementationVersion_, ES_steppingLength_, ES_scaling_);
//    }
//    if(takeEstimatedSlope_[storageIdx]){
//      takeEstimatedSlope_[storageIdx] = false;
//      return EstimateSlope(Originallpm, ES_includeStrains_, ES_useAbs_, "DeltaMat", ES_implementationVersion_, ES_steppingLength_, ES_scaling_);
//    }
//
//
//    if(useStrains){
//      //      std::cout << "GetDeltaMat - for strains" << std::endl;
//      if((deltaMatStrain_requiresReeval_[storageIdx] == false)||(reuseOldState)){
//        //        std::cout << "NO reeval" << std::endl;
//        return deltaMatStrain_[storageIdx];
//      }
//
//      if(dim_ == 3){
//        numRows = 6;
//      } else {
//        numRows = 3;
//      }
//      deltaMat = Matrix<Double>(numRows,numCols);
//
//      /*
//       * Compute matrix dM, such that
//       * dS = dM*dE (or dB)
//       *
//       * with dS = Strains(timeLevel1)-Strains(timeLevel2)
//       *      dE = E(timeLevel1)-E(timeLevel2)
//       */
//
//      Vector<Double> E_B_new = RetrieveLPMSolution(actualLPM, storageIdx, timelevel_new, onBoundary);
//      Vector<Double> E_B_old = RetrieveLPMSolution(actualLPM, storageIdx, timelevel_old, onBoundary);
//
//      Vector<Double> E_B_diff = E_B_new;
//      E_B_diff -= E_B_old;
//
//      Vector<Double> S_new = GetIrreversibleStrains(Originallpm,timelevel_new);
//      Vector<Double> S_old = GetIrreversibleStrains(Originallpm,timelevel_old);
//      Vector<Double> S_diff = S_new;
//      S_diff -= S_old;
//
//      // cutting does not work so well either
//      Double cuttingTol = 1e-5*S_old.NormL2();
//
//      Matrix<Double> deltaTMP = CalcDeltaMatStrains(E_B_diff, S_diff, useAbs, implementationVersion,cuttingTol);
//
//      // Blending with old deltaMat does not work as expected
//      if( (includeOldDelta_ > 0)&&(includeOldDelta_ <= 1)){
//        deltaMat.Init();
//        deltaMat.Add(1.0-includeOldDelta_,deltaTMP);
//        deltaMat.Add(includeOldDelta_,deltaMatStrainPrev_[storageIdx]);
//      } else {
//        deltaMat = deltaTMP;
//      }
//
//      deltaMatStrain_[storageIdx] = deltaMat;
//      deltaMatStrain_requiresReeval_[storageIdx] = false;
//
//    } else {
//
//      if((deltaMat_requiresReeval_[storageIdx] == false)||(reuseOldState)){
//        LOG_DBG(coeffunctionhyst_main) << "Reuse old deltaMat";
////        std::cout << "Reuse old deltaMat" << std::endl;
//        return deltaMat_[storageIdx];
//      }
//
//      // compare deltaMat IT vs deltaMat TS
//      bool compareDeltaMatrices = !true;
//      if(compareDeltaMatrices && (storageIdx == 1)){
//
//        LOG_DBG2(coeffunctionhyst_main) << "Compute deltaMat towards last TS for storage idx " << storageIdx;
//        numRows = dim_;
//        Matrix<Double> deltaMatTS = Matrix<Double>(numRows,numCols);
//        Matrix<Double> deltaMatIT = Matrix<Double>(numRows,numCols);
//
//        Vector<Double> E_B_cur = RetrieveLPMSolution(actualLPM, storageIdx, timelevel_new, onBoundary);
//        Vector<Double> E_B_ts = RetrieveLPMSolution(actualLPM, storageIdx, -1, onBoundary);
//        Vector<Double> E_B_it = RetrieveLPMSolution(actualLPM, storageIdx, 1, onBoundary);
////
////        std::cout << "E_B_ at " << std::endl;
////        std::cout << "Current IT: " << E_B_cur.ToString() << std::endl;
////        std::cout << "Previous IT: " << E_B_it.ToString() << std::endl;
////        std::cout << "Previous TS: " << E_B_ts.ToString() << std::endl;
//
//        Vector<Double> E_B_diff_it = E_B_cur;
//        E_B_diff_it -= E_B_it;
//        Vector<Double> E_B_diff_ts = E_B_cur;
//        E_B_diff_ts -= E_B_ts;
//
//      Vector<Double> P_J_cur = GetPrecomputedOutputOfHysteresisOperator(Originallpm,timelevel_new,false);
//      Vector<Double> P_J_ts = GetPrecomputedOutputOfHysteresisOperator(Originallpm,-1,false);
//      Vector<Double> P_J_it = GetPrecomputedOutputOfHysteresisOperator(Originallpm,1,false);
//
////        std::cout << "P_J__ at " << std::endl;
////        std::cout << "Current IT: " << P_J_cur.ToString() << std::endl;
////        std::cout << "Previous IT: " << P_J_it.ToString() << std::endl;
////        std::cout << "Previous TS: " << P_J_ts.ToString() << std::endl;
////
//      Vector<Double> P_J_diff_it = P_J_cur;
//      P_J_diff_it -= P_J_it;
//          Vector<Double> P_J_diff_ts = P_J_cur;
//      P_J_diff_ts -= P_J_ts;
//
//      Double cuttingTol = 1e-5*P_J_it.NormL2();
//      deltaMatIT = CalcDeltaMat(E_B_diff_it, P_J_diff_it, useAbs, implementationVersion,cuttingTol);
////      std::cout << "Computed deltaMat towards last IT: " << deltaMatIT.ToString() << std::endl;
//
//      cuttingTol = 1e-5*P_J_ts.NormL2();
//      deltaMatTS = CalcDeltaMat(E_B_diff_ts, P_J_diff_ts, useAbs, implementationVersion,cuttingTol);
////      std::cout << "Computed deltaMat towards last TS: " << deltaMatTS.ToString() << std::endl;
//
//      }
//
//
//      /*
//       * Compute matrix dM, such that
//       * dP = dM*dE (or dB)
//       *
//       * with dP = Polarization(timeLevel1)-Polarizatoin(timeLevel2)
//       *      dE = E(timeLevel1)-E(timeLevel2)
//       */
//      numRows = dim_;
//      deltaMat = Matrix<Double>(numRows,numCols);
//
//      Vector<Double> E_B_new = RetrieveLPMSolution(actualLPM, storageIdx, timelevel_new, onBoundary);
//      Vector<Double> E_B_old = RetrieveLPMSolution(actualLPM, storageIdx, timelevel_old, onBoundary);
//
////      if(storageIdx == 1){
////        std::cout << "Timelevel new: " << timelevel_new << std::endl;
////        std::cout << "Timelevel_old: " << timelevel_old << std::endl;
////      }
//
//      LOG_DBG(coeffunctionhyst_deltamat) << "Old solution of system (solution at timelevel = " << timelevel_old << "): " << E_B_old.ToString();
//      LOG_DBG(coeffunctionhyst_deltamat) << "Current solution of system (solution at timelevel = " << timelevel_new << "): " << E_B_new.ToString();
//
//      // does not work well if true
//      bool cut_E_B_toSat = !true;
//      if(cut_E_B_toSat){
//
//        Double satValue;
//        if(PDEName_ == "Electromagnetics"){
//          // will not  be exactly correct when taking the norm but gives beter estimate than skipping it
//          satValue = POL_operatorParams_.outputSat_ + eps_nu_base_.NormL2()*POL_operatorParams_.inputSat_;
//        } else {
//          satValue = POL_operatorParams_.inputSat_;
//        }
//
//        if(E_B_new.NormL2() > satValue){
//          E_B_new.ScalarMult(satValue/E_B_new.NormL2());
//        }
//        if(E_B_old.NormL2() > satValue){
//          E_B_old.ScalarMult(satValue/E_B_old.NormL2());
//        }
//
//        LOG_DBG(coeffunctionhyst_deltamat) << "After scaling to saturation: ";
//        LOG_DBG(coeffunctionhyst_deltamat) << "Old solution of system (solution at timelevel = " << timelevel_old << "): " << E_B_old.ToString();
//        LOG_DBG(coeffunctionhyst_deltamat) << "Current solution of system (solution at timelevel = " << timelevel_new << "): " << E_B_new.ToString();
//
//      }
//
//
//      //std::cout << "Old solution (solution at timelevel = " << timelevel_old << "): " << E_B_old.ToString() << std::endl;
//      //std::cout << "Current solution (solution at timelevel = " << timelevel_new << "): " << E_B_new.ToString() << std::endl;
//
//      Vector<Double> E_B_diff = E_B_new;
//      E_B_diff -= E_B_old;
//
//      Vector<Double> P_J_new = GetPrecomputedOutputOfHysteresisOperator(Originallpm,timelevel_new,false);
//      Vector<Double> P_J_old = GetPrecomputedOutputOfHysteresisOperator(Originallpm,timelevel_old,false);
//
//      LOG_DBG(coeffunctionhyst_deltamat) << "Old output of hyst operator: " << P_J_old.ToString();
//      LOG_DBG(coeffunctionhyst_deltamat) << "Current output of hyst operator: " << P_J_new.ToString();
//
//      //std::cout << "Old output of hystoperator (solution at timelevel = " << timelevel_old << "): " << P_J_old.ToString() << std::endl;
//      //std::cout << "Current output of hystoperator (solution at timelevel = " << timelevel_new << "): " << P_J_new.ToString() << std::endl;
//
//      Vector<Double> P_J_diff = P_J_new;
//      P_J_diff -= P_J_old;
//
//      Double cuttingTol = 1e-5*P_J_old.NormL2();
//
//      Matrix<Double> deltaTMP = CalcDeltaMat(E_B_diff, P_J_diff, useAbs, implementationVersion,cuttingTol);
//      if( (includeOldDelta_ > 0)&&(includeOldDelta_ <= 1)){
//        deltaMat.Init();
//        deltaMat.Add(1.0-includeOldDelta_,deltaTMP);
//        deltaMat.Add(includeOldDelta_,deltaMatPrev_[storageIdx]);
//      } else {
//        deltaMat = deltaTMP;
//      }
//
//      deltaMat_[storageIdx] = deltaMat;
//      deltaMat_requiresReeval_[storageIdx] = false;
//
//
//    if(miniOutput){
//      std::cout << "Calculated deltaMat approximation with" << std::endl;
//      std::cout << "E_B_new (timelevel_new=" << timelevel_new <<"): " << E_B_new.ToString() << std::endl;
//      std::cout << "E_B_old (timelevel_old=" << timelevel_old <<"): " << E_B_old.ToString() << std::endl;
//    }
//
//    }
//    LOG_DBG(coeffunctionhyst_deltamat) << "Computed deltaMat for storageIDX " << storageIdx << ": " << deltaMat.ToString();
//
////    std::cout << "Storage Idx = " << storageIdx << std::endl;
////    std::cout << "Jacobian: " << estimatedSlope.ToString() << std::endl;
////    std::cout << "DeltaMat: " << deltaMat.ToString() << std::endl;
////
//    if(useJacobian){
//
//      if(miniOutput){
//        std::cout << "Return Jacobian approximation" << std::endl;
//        std::cout << "FD-Jacobian: " << estimatedSlope.ToString() << std::endl;
//      }
//      return estimatedSlope;
//    }
//    if(miniOutput){
//      std::cout << "Return deltaMat approximation" << std::endl;
//      std::cout << "deltaMat: " << deltaMat.ToString() << std::endl;
////      std::cout << "for comparison - FDJacobian: " << estimatedSlope.ToString() << std::endl;
//    }
//
//    return deltaMat;
//
//  }


  Matrix<Double> CoefFunctionHyst::GetIrreversibleStrainTensor(const LocPointMapped& Originallpm, int timeLevel, bool forceMidpoint) {

    Vector<Double> Si_voigt = GetIrreversibleStrains(Originallpm, timeLevel, forceMidpoint);
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
  Vector<Double> CoefFunctionHyst::GetIrreversibleStrains(const LocPointMapped& Originallpm, int timeLevel, bool forceMidpoint) {

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
    UInt operatorIdx, storageIdx,idxElem;
		LocPointMapped actualLPM;

		PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx,idxElem,forceMidpoint);

		if (timeLevel == -1) {
      if(lastTSStorageInitialized_ == false){
        InitSpecificStorage(1);
      }
			// return value from last time step
			return Si_lastTS_[storageIdx];
		} else if (timeLevel == 1) {
      if(lastItStorageInitialized_ == false){
        InitSpecificStorage(0);
      }
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

//    std::cout << "ComputeIrreversibleStrains; P = " << P.ToString() << "; E_H = " << E_H.ToString() << std::endl;

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
//      // Model strains as discribed in "Generalisiertes Preisach-Modell f��r die Simulation und Kompensation
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
          int timeLevel, bool forStrain, bool forceMidpoint) {
    LOG_TRACE(coeffunctionhyst_main) << "GetPrecomputedOutputOfHysteresisOperator for timelevel " << timeLevel;
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
		UInt operatorIdx, storageIdx,idxElem;
		LocPointMapped actualLPM;

    //		bool onBoundary =     PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx);

    PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx,idxElem,forceMidpoint);

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
      if(lastTSStorageInitialized_ == false){
        InitSpecificStorage(1);
      }

			LOG_TRACE(coeffunctionhyst_main) << "> prevTS value";
      // return value from last time step
      //std::cout << "Get value from last TS: " << P_J_lastTS_[storageIdx].ToString() << std::endl;
      if(forStrain){
        return P_J_forStrains_lastTS_[storageIdx];
      } else {
        return P_J_lastTS_[storageIdx];
      }
		} else if (timeLevel == 1) {
      if(lastItStorageInitialized_ == false){
        InitSpecificStorage(0);
      }

      LOG_TRACE(coeffunctionhyst_main) << "> prevIt value";
			// return value from last iteration
      if(forStrain){
        EXCEPTION("Polarization for strains for last iteration should not be called");
      } else {
//        if(storageIdx == 1){
//          std::cout << "Get P_J_lastIt_" << std::endl;
//        }
        return P_J_lastIt_[storageIdx];
      }
		} else {
      if(forStrain){
        return P_J_forStrains_[storageIdx];
      } else {
//        std::cout << "Get P_J_[storageIdx] " << P_J_[storageIdx].ToString() << std::endl;
        return P_J_[storageIdx];
      }

      //      OLD
      //
      //      LOG_TRACE(coeffunctionhyst_main) << "> current value";
      //      if(hystOperatorLocked_){
      //        // no evaluation > return current state
      //				LOG_TRACE(coeffunctionhyst_main) << "Output: hystOperatorLocked_ locked";
      //        return P_J_[storageIdx];
      //      }
      //
      //			// get current state; here we have to check if a reevaluation is needed
      //			// note: requiresReeval_ has one entry for each storage, not for each operator
      //			if (false == requiresReeval_[storageIdx]) {
      //        LOG_TRACE(coeffunctionhyst_main) << "> NO reeval; take stored value";
      //				//std::cout << "no reeval needed as flag is false!" << std::endl;
      //				// return current value
      //				return P_J_[storageIdx];
      //			}
		}
    //    LOG_TRACE(coeffunctionhyst_main) << "> reeval needed";
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
      if(lastTSStorageInitialized_ == false){
        InitSpecificStorage(1);
      }
			// get solution from last ts
			return E_B_lastTS_[storageIdx];
		} else if (timeLevel == 1) {
      if(lastItStorageInitialized_ == false){
        InitSpecificStorage(0);
      }
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

	Vector<Double> CoefFunctionHyst::RetrieveInputToHysteresisOperator(LocPointMapped& actualLPM, UInt operatorIdx, UInt storageIdx, bool onBoundary, bool overwriteMemory) {
    /*
     * New comments/description - January 2019
     *
     * Purpose of this function:
     *  This function returns the input-vector x to the forward hysteresis operator
     *  or y in case of an inverse hysterese operator (currently none implemented).
     *
     * Possible cases:
     *  Case A: Global FE system computes x or an adequate potential directly (e.g. electrostatics)
     *  > return element solution x
     *
     *  Case B: Global FE system computes y or an adequate potential directly where y = c*x + p (e.g. magnetics)
     *  and hysteresis operator p is a function of y, i.e. p = p(y) (currently no such model implemented in CFS)
     *  > return element solution y
     *
     *  Case C: Global FE system computes y or an adequate potential directly but hysteresis operator is a function
     *  of x (magnetic case, basically all implemented Preisach models)
     *  > subcase C1: iteratively compute x_k^j such that res_k^j = y_k - c*x_k^j + p(x_k^j) < tolerance
     *    > done in hysteresis classes via function compute_input_vec
     *    > note: j is the inner iteration counter; k the outer iteration counter (e.g. from solvestep)
     *    > relevant for standard Fixpoint methods (B-Version) and all Newton-like solution schemes for the global system
     *
     *  > subcase C2:
     *    1) estimate x_est from previous values
     *        x_est = (y_k - p(x_k-1))/c
     *    2) perform a local linesearch or take a fix appropriate relaxation parameter \omega to compute x_k
     *        x_k = (1-\omega)*x_k-1 + \omega x_est
     *
     *    > relevant for the so-called H-Version of the Fixpoint method
     *    > compared to the B-Version can we avoid the local inversion (can be teadious and expensiv) but due to the
     *        worse values of x_k we need more outer iterations
     *    > for a comparison between B and H Version see e.g.
     *        "Newton-Raphson method and fixed-point technique in finite element\ncomputation of magnetic field problems in media with hysteresis" - Saitz
     */


    LOG_TRACE(coeffunctionhyst_main) << "RetrieveInputToHysteresisOperator";
    /*
     * NOTE/TODO:
     *  in case of magnetostrictive coupling, we actually would have to solve the nonlinear
     *  problem
     *    H = -[h](S - Si(H)) + nu(B - P(H))
     *
     *  problem 1:
     *  > the inversion procedure based on finding H = \nu (B - P(H)) will not work directly
     *
     *  > solution approach:
     *      instead of feeding the local inversion method with B, we use B^mod instead with
     *      B^mod = B + \nu_inv [h](S - Si(H_old))
     *  > needs to be implemented yet
     *
     *  problem 2:
     *  > in magnetostrictive systems, nu can be a function of B if nu^T is given (we need nu^S which
     *    computes from nu^T, [h] and [c]) and [h] is scaled and rotated by P
     *  > inversion must not use some constant \nu but has to use the last knowon value of \nu
     *
     *  > solution approach:
     *      retrive \nu(P(H_k-1)) via the computeTensor function and pass this \nu (or its inverse) to the
     *      hyst operator
     *  > should already be implemented but needs check!
     */

    UInt timeLevel = 0;
    if( (needsInversion_ == false)||(POL_operatorParams_.hasInverseModel_ == true) ){
      /*
       * Forward model can be used directly > Case A
       * or
       * Backward model directly implemented > Case B
       *
       * Case A (electrostatics) > return E
       * Case B (magnetics) > return B
       *
       */
      if(hystOperatorLocked_){
        // no evaluation > return last stored input to hyst operator
        LOG_TRACE(coeffunctionhyst_main) << "Input: hystOperator locked";
        return E_B_[storageIdx];
      } else {
        // retrieve current system solution (E for electrostatics, B for magnetics)
        E_B_[storageIdx] = RetrieveLPMSolution(actualLPM, storageIdx, timeLevel, onBoundary);
        return E_B_[storageIdx];
      }
    } else {
      /*
       * Forward model used but system solution corresponds to output > Case C
       */
      if(hystOperatorLocked_){
        // no evaluation > return last stored input to hyst operator
        LOG_TRACE(coeffunctionhyst_main) << "Input: hystOperator locked";
//        std::cout << "Hyst operator locked" << std::endl;
        return E_H_[storageIdx];
      }

      Vector<Double> retrievedInput = Vector<Double>(dim_);
      Vector<Double> tmp = Vector<Double>(dim_);
      Vector<Double> curLPMSolution = RetrieveLPMSolution(actualLPM, storageIdx, timeLevel, onBoundary);

      /*
       * at this point we could subtract magnetostrictive parts as described in the comment above to
       * get an effective B^mod that is feeded into inversion or which is used for the estimation for the H-Version
       */

      if(useFPH_ == false){
        /*
         * Case C1 > perform local inversion
         */
        // check if solution did change since last time > if not, hope that input / output pair of hyst
        // operator can be reused
        Vector<Double> diff = curLPMSolution;
        Vector<Double> toDiff = E_B_[storageIdx];
        diff -= toDiff;

        if (diff.NormL2() < POL_operatorParams_.amplitudeResolution_) {
          LOG_TRACE(coeffunctionhyst_main) << "> difference " << diff.NormL2() << " below tolerance > reuse";
          retrievedInput = E_H_[storageIdx];
        } else {
          LOG_TRACE(coeffunctionhyst_main) << "> calculate new input";
          //bool overwriteMemory = OverwriteHystMemory();
          // for input computation never overwrite storage
          // we basically search in the current state for a fitting solution
          // then, during the step calcoutput we evaluate the hyst operator with
          // the retrieved input and here we actually overwrite the memory, if necessary
          /*
           * Update 7.5.2020
           * due to additional inversion routine for Mayergoyz model, the above note no longer holds
           * in all cses; when we use the Everett-based inversion for the Mayergoyz model (i.e., we define
           * the Mayergoyz model not in terms of forward models but in terms of inverted forward models),
           * we cannot use the forward evaluation anymore as the computed vector input will not return the
           * originally applied output (see VectorPreisachMayergoyz for more info); instead each inverted forward
           * model must overwrite its state directly which can be triggered by overwriteMemory = true
           */
//          bool overwriteMemory = false;
          if((InversionParams_.inversionMethod != LOCAL_NOTIMPLEMENTED)&&(overwriteMemory)){
            //              inversionMethod == 10 > everett for mayergoyz; only in this case overwriteMemroy might be true; 
            //              otherwise enforce overwriteMemory = false
            WARN("RetrieveInputToHysteresisOperator with flag overwriteMemory=true is only allowed for Mayergoyz vector model"
                    "using inverse scalar models!");
            overwriteMemory = false;
          }
          
          
          if(eps_mu_local_Set_[storageIdx] == false){
            EXCEPTION("Material matrix for inversion was not set yet! This has to be done in the calling helper classes!");
          } else {
            LOG_TRACE(coeffunctionhyst_main) << "Matrix for inversion (storage="<<storageIdx<<"): "<<eps_mu_local_[storageIdx].ToString();
          }
          if(POL_operatorParams_.methodType_ == 0){
            /*
             * Scalar hysteresis operator
             */
            if(InversionParams_.inversionMethod != 4){
              WARN("Inversion of scalar model should always use Everett based implementation (except for testing purposes); will ignore selected inversion method and use EverettBasedInversion instead.")
            }
            /*
             * Get material tensor (see comment above; matrixForInversion should contain the current value of mu = inv(nu^S))
             */
            Double muForInversion;
            eps_mu_local_[storageIdx].Mult(POL_operatorParams_.fixDirection_,tmp);
            POL_operatorParams_.fixDirection_.Inner(tmp,muForInversion);

            /*
             * Compute scalar input to hyst operator
             */
            Double scalInput, scalOutput;
            POL_operatorParams_.fixDirection_.Inner(curLPMSolution,scalOutput);
            int successFlag = 0;
            scalInput = hyst_->computeInputAndUpdate(scalOutput,muForInversion, operatorIdx, overwriteMemory, successFlag);

            retrievedInput = Vector<Double>(dim_);
            retrievedInput.Init();
            retrievedInput.Add(scalInput,POL_operatorParams_.fixDirection_);

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
            curLPMSolution_perpendicular.Add(1.0,curLPMSolution,-scalOutput,POL_operatorParams_.fixDirection_);
            Vector<Double> retrievedInput_perpendicular = Vector<Double>(dim_);
            // epsInv_nu_local_[storageIdx] = nu
            epsInv_nu_local_[storageIdx].Mult(curLPMSolution_perpendicular,retrievedInput_perpendicular);
            retrievedInput.Add(1.0,retrievedInput_perpendicular);

          } else {
            /*
             * Vector hysteresis model (or scalar model that is forced to be inverted using nonlinear algorithm)
             */
            int successFlag = 0;
            bool useEverett = false;
            
            if(InversionParams_.inversionMethod == 10){
              useEverett = true;
            }
            retrievedInput = hyst_->computeInput_vec(curLPMSolution,operatorIdx,eps_mu_local_[storageIdx],
                    POL_operatorParams_.fieldsAlignedAboveSat_,POL_operatorParams_.hystOutputRestrictedToSat_,successFlag, useEverett, overwriteMemory);

          }
        }
      } else {
        /*
         * Case C2 > estimate input only; input will not solve the material law exactly (or nearly exactly)
         *
         *    1) estimate x_est from previous values
         *        x_est = (y_k - p(x_k-1))/c
         *    2) perform a local linesearch or take a fix appropriate relaxation parameter \omega to compute x_k
         *        x_k = (1-\omega)*x_k-1 + \omega x_est
         */
        Vector<Double> H_estimate = Vector<Double>(dim_);
        // P_J_ should not have been update yet (the input for this update is currently evaluated), i.e.
        // P_J_ = hystOperator(H_k-1)
        tmp.Add(1.0,curLPMSolution,-1.0,P_J_[storageIdx]);
        // epsInv_nu_local_[storageIdx] = nu
        epsInv_nu_local_[storageIdx].Mult(tmp,H_estimate);

        bool miniOutput = false;
//        if(storageIdx == 1){
//          miniOutput = true;
//        }
        if(miniOutput){
          std::cout << "USE FPH" << std::endl;
          std::cout << "curLPMSolution = " << curLPMSolution.ToString() << std::endl;
          std::cout << "P_J_[storageIdx] = " << P_J_[storageIdx].ToString() << std::endl;
          std::cout << "epsInv_nu_local_[storageIdx] = " << epsInv_nu_local_[storageIdx].ToString() << std::endl;
          std::cout << "H_estimate = " << H_estimate.ToString() << std::endl;
        }

        // omega can either come from linesearch (not implemented here) or be fix; according to
        // "Anwendung der Methode der finiten Elemente und des Vektor-Preisach-Modells zur Berechnung ebener magnetostatischer Felder in hysteresebehafteten Medien" - Kurz et al.
        // omega should between 0 and 2/\mu_r_max
        // with mu_r_max*mu_0 = maxSlopeGlobal_ of forward hyst model
        //
        // TODO: maxSlopeGlobal_ has to be determined by tracing the hyst operator
        // > must not be done here as this function might be called in parallel and then it would get messed up
        // > beter: once useFPH is set to true, trace hyst operator
        Double mu_0 = eps_mu_local_[storageIdx].GetMax(); //4*M_PI*1e-7;
//        Double mu_FP = 0.5*(maxSlopeGlobal_ + minSlopeGlobal_);
        // used so far and works quite good
//        Double mu_FP = maxSlopeGlobal_;
        Double mu_FP = muFPH_[storageIdx];

//        Double mu_FP = (maxSlopeGlobal_+ minSlopeGlobal_)/2;
        Double omega;
        if(mu_FP != 0){
          omega = mu_0/mu_FP;
        } else {
          // mu FP not set yet
          omega = 0;
        }

//        if(operatorIdx == 1){
//          std::cout << "#### FP H ####" << std::endl;
////          std::cout << "maxSlopeGlobal_ = " << maxSlopeGlobal_ << std::endl;
////          std::cout << "omega = " << omega << std::endl;
//
//          Matrix<Double> FP_Tensor = Matrix<Double>(dim_,dim_);
//          Matrix<Double> baseTensor = Matrix<Double>(dim_,dim_);
//          Matrix<Double> tmpTensor = Matrix<Double>(dim_,dim_);
//          // we would need original lpm here, but for testing it could work
//          ComputeTensor(baseTensor, actualLPM, "Reluctivity", "DEFAULT", false, false, false, true );
//          baseTensor.Invert(tmpTensor);
//          FP_Tensor = GetMaterialRelation(actualLPM, 0, 3, tmpTensor, "Reluctivity",true);
//          std::cout << "estimated FP-Tesnor (nu) at current position: " << std::endl;
//          std::cout << FP_Tensor.ToString() << std::endl;
//          FP_Tensor.Invert(tmpTensor);
//          std::cout << "estimated FP-Tesnor (mu) at current position: " << std::endl;
//          std::cout << tmpTensor.ToString() << std::endl;
//          std::cout << "maxSlope in mu = " << tmpTensor.GetMax() << std::endl;
//          std::cout << "maxSlope (global estimate) = " << maxSlopeGlobal_ << std::endl;
//          if(tmpTensor.GetMax()-maxSlopeGlobal_ > 10*maxSlopeGlobal_){
//            std::cout << "## WARN!!!!: " << tmpTensor.GetMax()-maxSlopeGlobal_ << std::endl;
//          }
//        }

        retrievedInput.Add(1.0-omega,E_H_[storageIdx],omega,H_estimate);
        
        if(miniOutput){
          std::cout << "omega = " << omega << std::endl;
          std::cout << "E_H_[storageIdx] = " << E_H_[storageIdx].ToString() << std::endl;
          std::cout << "retrievedInput = " << retrievedInput.ToString() << std::endl;
        }
        
//        if(storageIdx == 1){
//          std::cout << "E_H_[storageIdx] = " << E_H_[storageIdx].ToString() << std::endl;
//          std::cout << "retrievedInput = " << retrievedInput.ToString() << std::endl;
//        }

      }
      //std::cout << "Computed input: " << retrievedInput.ToString() << std::endl;
      // we should not store the value here but first after evaluating the hyst operator in calcOutput function
      // otherwise we will never get to the point where P_J_ gets overwritten
      // and thus we get no real stepping at all
      //E_H_[storageIdx] = retrievedInput;

      // store current solution for later checks
      E_B_[storageIdx] = curLPMSolution;
      return retrievedInput;
    }
  }

  Vector<Double> CoefFunctionHyst::RetrieveInputToHysteresisOperatorOLD(LocPointMapped& actualLPM, UInt operatorIdx, UInt storageIdx, bool onBoundary) {
    EXCEPTION("RetrieveInputToHysteresisOperatorOLD no longer supported");
    
    if(hystOperatorLocked_){
      // no evaluation > return current state
			LOG_TRACE(coeffunctionhyst_main) << "Input: hystOperatorLocked_ locked";
      return E_H_[storageIdx];
    }

		UInt timeLevel = 0;

		Vector<Double> curLPMSolution = RetrieveLPMSolution(actualLPM, storageIdx, timeLevel, onBoundary);
//    std::cout << "curLPMSolution: " << curLPMSolution.ToString() << std::endl;
		if(needsInversion_ && (POL_operatorParams_.hasInverseModel_ == false) ) {

      LOG_TRACE(coeffunctionhyst_main) << "Inversion needed";
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
        LOG_TRACE(coeffunctionhyst_main) << "> difference " << diff.NormL2() << " below tolerance > reuse";
        //std::cout << "(Nearly) no difference in amplitude to previous element solution" << std::endl;
        retrievedInput = E_H_[storageIdx];
				//}
			} else {
        LOG_TRACE(coeffunctionhyst_main) << "> calculate new input";
        //bool overwriteMemory = OverwriteHystMemory();
        // for input computation never overwrite storage
        // we basically search in the current state for a fitting solution
        // then, during the step calcoutput we evaluate the hyst operator with
        // the retrieved input and here we actually overwrite the memory, if necessary
        bool overwriteMemory = false;

        if(eps_mu_local_Set_[storageIdx] == false){
          EXCEPTION("Material matrix for inversion was not set yet! This has to be done in the calling helper classes!");
        } else {
          LOG_TRACE(coeffunctionhyst_main) << "Matrix for inversion (storage="<<storageIdx<<"): "<<eps_mu_local_[storageIdx].ToString();
        }
        if(POL_operatorParams_.methodType_ == 0){
          if(InversionParams_.inversionMethod != 4){
            WARN("Inversion of scalar model should always use Everett based implementation (except for testing purposes); will ignore selected inversion method and use EverettBasedInversion instead.")
          }

          retrievedInput = Vector<Double>(dim_);

          //					retrievedInput[POL_operatorParams_.fixDirection_] = hyst_->computeInputAndUpdate(curLPMSolution[POL_operatorParams_.fixDirection_],
          //						E_B_lastIt_[storageIdx][POL_operatorParams_.fixDirection_], E_H_lastIt_[storageIdx][POL_operatorParams_.fixDirection_],P_J_lastIt_[storageIdx][POL_operatorParams_.fixDirection_],
          //						operatorIdx, eps_nu_base_[POL_operatorParams_.fixDirection_][POL_operatorParams_.fixDirection_], overwriteMemory);
          // NEW: POL_operatorParams_.fixDirection_ is a vector now
          // > to get the correct scalar input data, we have to project onto this direction vector
          Double muForInversion;
          Vector<Double> tmp = Vector<Double>(dim_);
          eps_mu_local_[storageIdx].Mult(POL_operatorParams_.fixDirection_,tmp);
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
          epsInv_nu_local_[storageIdx].Mult(curLPMSolution_perpendicular,retrievedInput_perpendicular);
          retrievedInput.Add(1.0,retrievedInput_perpendicular);

          //          retrievedInput[POL_operatorParams_.fixDirection_] = hyst_->computeInputAndUpdate(curLPMSolution[POL_operatorParams_.fixDirection_],
          //						muForInversion, operatorIdx, overwriteMemory);
        } else {
          // Idea: use Levenberg Marquart to estimate  inversion of hyst operator
          //          retrievedInput = hyst_->computeInput_vec(curLPMSolution, operatorIdx, eps_nu_base_,
          //                  alphaLinesearch_, overwriteMemory, RUN_overwriteDirection_ );

          //					retrievedInput = hyst_->computeInput_vec(curLPMSolution, E_B_lastIt_[storageIdx], E_H_lastIt_[storageIdx],
          //						P_J_lastIt_[storageIdx],operatorIdx,eps_nu_base_,RUN_overwriteDirection_);
          //
          int successFlag = 0;

          retrievedInput = hyst_->computeInput_vec(curLPMSolution,operatorIdx,eps_mu_local_[storageIdx],
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
    LOG_TRACE(coeffunctionhyst_main) << "Clip direction of target vector to next full " << POL_operatorParams_.angularClipping_ << " degree";
    LOG_DBG(coeffunctionhyst_main) << "Original vector: " << targetVector.ToString();
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
      LOG_DBG(coeffunctionhyst_main) << "Circular/Polar coordinates (r,alpha): " << radius << "," << alphaDeg;

      // apply clipping
      tmp = alphaDeg/POL_operatorParams_.angularClipping_;
      tmp = round(tmp);
      alphaDeg = tmp*POL_operatorParams_.angularClipping_;

      LOG_DBG(coeffunctionhyst_main) << "Circular/Polar coordinates after clipping (r,alpha): " << radius << "," << alphaDeg;

      // now rebuild output vector with clipped coordinates
      targetVector[0] = radius*cos(alphaDeg/180*M_PI);
      targetVector[1] = radius*sin(alphaDeg/180*M_PI);

      LOG_DBG(coeffunctionhyst_main) << "Rebuild vector after clipping: " << targetVector.ToString();

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

      LOG_DBG(coeffunctionhyst_main) << "Rebuild vector after further treatment: " << targetVector.ToString();
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

      LOG_DBG(coeffunctionhyst_main) << "Spherical coordinates (r,theta,phi): " << radius << "," << thetaDeg << "," << phiDeg;

      // apply clipping
      tmp = thetaDeg/POL_operatorParams_.angularClipping_;
      tmp = round(tmp);
      thetaDeg = tmp*POL_operatorParams_.angularClipping_;

      tmp = phiDeg/POL_operatorParams_.angularClipping_;
      tmp = round(tmp);
      phiDeg = tmp*POL_operatorParams_.angularClipping_;

      LOG_DBG(coeffunctionhyst_main) << "Spherical coordinates after clipping (r,theta,phi): " << radius << "," << thetaDeg << "," << phiDeg;

      targetVector[0] = radius*sin(thetaDeg/180*M_PI)*cos(phiDeg/180*M_PI);
      targetVector[1] = radius*sin(thetaDeg/180*M_PI)*sin(phiDeg/180*M_PI);
      targetVector[2] = radius*cos(thetaDeg/180*M_PI);

      LOG_DBG(coeffunctionhyst_main) << "Rebuild vector after clipping: " << targetVector.ToString();

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

      LOG_DBG(coeffunctionhyst_main) << "Rebuild vector after further treatment: " << targetVector.ToString();

    }
  }

  Vector<Double> CoefFunctionHyst::GetPrecomputedInputToHysteresisOperator(const LocPointMapped& Originallpm, int timeLevel, bool forceMidpoint){

		UInt operatorIdx, storageIdx,idxElem;
		LocPointMapped actualLPM;

    PreprocessLPM(Originallpm, actualLPM, operatorIdx, storageIdx,idxElem,forceMidpoint);

    if(timeLevel == -1){
      if(lastTSStorageInitialized_ == false){
        InitSpecificStorage(1);
      }
      return E_H_lastTS_[storageIdx];
    } else if(timeLevel == 1){
      if(lastItStorageInitialized_ == false){
        InitSpecificStorage(0);
      }
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
    if(formerHystHelperParamsSet_ == false){
//      std::cout << "Create hystHelper" << std::endl;
//      std::cout << "- where? EvaluateHystOperators" << std::endl;
      PtrCoefFct fieldTensor;
      if (material_->GetClass() == ELECTROSTATIC) {
        fieldTensor = material_->GetTensorCoefFnc( ELEC_PERMITTIVITY_TENSOR,tensorType_,Global::REAL);
      }
      else if (material_->GetClass() == ELECTROMAGNETIC) {
        fieldTensor = material_->GetTensorCoefFnc( MAG_RELUCTIVITY_TENSOR,tensorType_,Global::REAL);
      }

      // stepwise change to merged class of coefFncHyst and HystHelper
      SetParamsOfFormerHystHelper(fieldTensor,elastTensorFct_,couplTensorFct_);
      formerHystHelperParamsSet_ = true;
    }

    /*
     * Init lienar Tensors; as we assume these tensors to be constant for each
     * region, we can take any lpm for init
     */
    if(allLPMmap_.empty()){
      EXCEPTION("Map of LPM was not set yet; have you forgotten to set evaluation depth?");
    }
    LocPointMapped anyLPM = allLPMmap_.begin()->second;

    // call own init function
    InitLinearTensors(anyLPM);

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
      PrecomputeMaterialTensorForInverison();
    }

    //    if(setRotStateForVecExtension){
    //      EvaluateRotDirectionForVectorExtension();
    //    }

    if(resetSolToZeroFirst){
      for(UInt i = 0; i < numStorageEntries_; i++){
        LOG_DBG(coeffunctionhyst_main) << "eps_mu_local_[ " << i << " ] = " << eps_mu_local_[i];
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
        LOG_DBG(coeffunctionhyst_main) << "E_H_[ " << i << " ] = " << E_H_[i];
        LOG_DBG(coeffunctionhyst_main) << "E_B_[ " << i << " ] = " << E_B_[i];
        LOG_DBG(coeffunctionhyst_main) << "P_J_[ " << i << " ] = " << P_J_[i];
      }
    }
    /*
     * -2 : not coupled at all;
     * -1 : not coupled via small signal parameter; irreversible strains still present
     *  0 : coupled e-form/h-form (piezo/magstrict)
     *  1 : coupled d-form (piezo)
     *  2 : coupled g-form (magstrict)
     */
    if(CouplingParams_.strainForm_ != -2){
//      std::cout << "Coupled with strainForm_ = " << strainForm_ << std::endl;
      PrecomputeIrrStrains();
      bool rotate = true;
      Matrix<Double> baseTensor = GetBaseCouplingTensor();
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
    LOG_DBG(coeffunctionhyst_main) << "CoefFunctionHyst::PrecomputeIrrStrains()";
    /*
     * Iterate over all lpm
     * for each lpm extract the precomputed values of P (current), P (lastTimestep)
     * use those values to calculate the irreversible strains and store them
     */

    bool forStrain = false;

    if(CouplingParams_.ownHystOperator_ && COUPLED_inXMLFile_ ){
      // own hyst operator defined for strains > evaluate that one (or load its value)
      // otherwise use polarization operator
      forStrain = true;
    }

    UInt strainSize = 3;
    if(dim_ == 3){
      strainSize = 6;
    }

    UInt size = allLPMmap_.size();
    UInt chunksize;
    UInt numT;
    DetermineChunkSizeAndThreads(size, chunksize, numT);

#ifdef USE_OPENMP
#pragma omp parallel num_threads(numT)
    {
//      UInt size = allLPMmap_.size();
      Vector<Double> P = Vector<Double>(dim_);
      Vector<Double> E_H = Vector<Double>(dim_);
      Vector<Double> Pold = Vector<Double>(dim_);
      Vector<Double> S_irr = Vector<Double>(strainSize);
      //std::cout << "Running parallel" << std::endl;
//      UInt numT = CFS_NUM_THREADS;
      UInt aThread = omp_get_thread_num();
      //std::cout << "My number: " << aThread << std::endl;
//      UInt chunksize = std::floor(size/numT);

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
  }

  void CoefFunctionHyst::PrecomputeScaledAndRotatedCouplingTensor(Matrix<Double> baseTensor, bool rotate){
    LOG_DBG(coeffunctionhyst_helper) << "CoefFunctionHelper::PrecomputeScaledAndRotatedCouplingTensor()";

    /*
     * Iterate over all integration points used/managed by the hyst operator and
     * for each of these LPM set the matrix for iversion, i.e. mu/nu in magnetic case
     */
    // obtain correct sizes first
    UInt numRows, numCols;
    numCols = baseTensor.GetNumCols();
    numRows = baseTensor.GetNumRows();

    if(numCols == 4){
      // axi case: 2x4 matrix
      WARN("Rotation for axi case not implemented yet. No rotation will be peformed");
      rotate = false;
    }

    UInt size = allLPMmap_.size();
    UInt chunksize;
    UInt numT;
    DetermineChunkSizeAndThreads(size, chunksize, numT);

#ifdef USE_OPENMP
#pragma omp parallel num_threads(numT)
    {
//      UInt size = allLPMmap_.size();
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
//      UInt numT = CFS_NUM_THREADS;
      UInt aThread = omp_get_thread_num();
      //std::cout << "My number: " << aThread << std::endl;
//      UInt chunksize = std::floor(size/numT);

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

        if(rotatedCouplingTensor_requiresReeval_[LPMit->first] == false){
//          std::cout << "No reeval of coupling tensor!" << std::endl;
          continue;
        }
//        std::cout << "Scale and rotate coupling tensor!" << std::endl;

        // get current value of P
        Vector<Double> P = GetPrecomputedOutputOfHysteresisOperator(LPMit->second, 0, false);
        scaling = P.NormL2()/POL_operatorParams_.outputSat_;

        scaledCouplTensor = baseTensor*scaling;
        rotatedCouplTensor = scaledCouplTensor;

        if(scaledCouplTensor.NormL2() != 0){
          //      std::cout << "Scalaed coupling tensor " << scaledCouplTensor.ToString() << std::endl;
          //
          if(rotate == true){
            //            if(P.NormL2() == 0){
            //              //        std::cout << "Polarization is zero > perform no rotation" << std::endl;
            //              rotatedCouplTensor = scaledCouplTensor;
            //
            //            } else {
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
          //          } else {
          //            // no rotation > take scaled tensor
          //            rotatedCouplTensor = scaledCouplTensor;
          //          }
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

      Vector<Double> P = Vector<Double>(dim_);

      // get current polarization (elec or mag)
      Vector<Double> dirP = Vector<Double>(dim_);
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

        if(rotatedCouplingTensor_requiresReeval_[LPMit->first] == false){
          continue;
        }

        // get current value of P
        Vector<Double> P = GetPrecomputedOutputOfHysteresisOperator(LPMit->second, 0, false);
        scaling = P.NormL2()/POL_operatorParams_.outputSat_;

        scaledCouplTensor = baseTensor*scaling;
        rotatedCouplTensor = scaledCouplTensor;

        if(scaledCouplTensor.NormL2() != 0){
          //      std::cout << "Scalaed coupling tensor " << scaledCouplTensor.ToString() << std::endl;
          //
          if(rotate == true){
            //            if(P.NormL2() == 0){
            //              //        std::cout << "Polarization is zero > perform no rotation" << std::endl;
            //              rotatedCouplTensor = scaledCouplTensor;
            //
            //            } else {
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
          //          } else {
          //            // no rotation > take scaled tensor
          //            rotatedCouplTensor = scaledCouplTensor;
          //          }
        }
        rotatedCouplingTensor_requiresReeval_[LPMit->first] = false;
        rotatedCouplingTensor_[LPMit->first] = rotatedCouplTensor;
        //        lastUsedTimeLevelForRotation_[storageIdx] = timeLevel;
      }

    //
    //      for(LPMit = LPMitStart; LPMit != LPMitEnd; LPMit++){
    //
    //        // get current value of P
    //        Vector<Double> P = GetPrecomputedOutputOfHysteresisOperator(LPMit->second, 0, false);
    //        scaling = P.NormL2()/POL_operatorParams_.outputSat_;
    //
    //        scaledCouplTensor = baseTensor*scaling;
    //
    //        if(scaledCouplTensor.NormL2() == 0){
    //          rotatedCouplTensor.Init();
    //        } else {
    //
    //          //      std::cout << "Scalaed coupling tensor " << scaledCouplTensor.ToString() << std::endl;
    //          //
    //          if(rotate == true){
    //            if(P.NormL2() == 0){
    //              //        std::cout << "Polarization is zero > perform no rotation" << std::endl;
    //              rotatedCouplTensor = scaledCouplTensor;
    //
    //            } else {
    //              dirP = P / P.NormL2();
    //
    //              if (numCols == 6){
    //                // 3d plane case
    //                // Idea: Create rotation matrix, that maps the current direction of P onto z-axis
    //                //       as the z-axis is the default polarization axis in the mat file
    //                //       To obtained the desired behavior, we have to rotate the z-axis onto P however.
    //                //       This can be done by taking the transposed rotation matrix.
    //
    //                // 1. rotate around z-axis by angle gamma such that P lies in z-y plane
    //                // gamma = angle between z-y plane and dirP
    //                gamma = std::atan2(dirP[0],dirP[1]);
    //
    //                // 2. rotate around x-axis by angle alpha such that P lies on top of z-axis
    //                // alpha = angle between x-y plane and z
    //                alpha = std::atan2(std::sqrt(dirP[1]*dirP[1]+dirP[0]*dirP[0]),dirP[2]);
    //
    //                // no rotation arouind y-axis needed > beta = 0
    //                // WARNING: this whole procedure is only valid if coupling tensor is at least
    //                // transverse isotropic! Otherwise we have to figure out on which axis to rotate
    //                // the transverse directions
    //                beta = 0.0;
    //                R = Compute3DRotationMatrix(alpha, beta, gamma);
    //
    //                // take transpose matrix to revert rotation (i.e. rotate z-axis onto dirP)
    //                RT.Resize(3,3);
    //                R.Transpose(RT);
    //
    //                //              assert(rotatedCouplTensor.GetNumRows() == 3);
    //                //              assert(rotatedCouplTensor.GetNumCols() == 6);
    //
    //                scaledCouplTensor.PerformRotation(RT,rotatedCouplTensor);
    //
    //              } else {
    //                // 2d plane strain or plane stress case
    //                // Important remark:
    //                //  for 2d plane strain and stress (and somehow also for axi) the coupling tensor
    //                //  will be rotated by default by alpha = -90 and gamma = -90.
    //                //  This results in the following mapping of the coordinate axis:
    //                //    z > y, y > x, x > z
    //                //  I.e. the material is rotated such that the default polarization is in +y direction.
    //                //  We thus have two possibilities to align the material tensor to the 2d polarization
    //                //  direction:
    //                //  a) get original 3x6 tensor and rotate that tensor according to the 3d version above
    //                //      then cut out subtensor
    //                //  b) take cut out subtensor and perform rotation directly in 2d
    //                //      > Version b done in the following
    //
    //                // std::atan2(dirP[0],dirP[1]) would rotate towards the y-axis
    //                // we want however the y-axis to rotate, so that we take -gamma instead
    //                gamma = -std::atan2(dirP[0],dirP[1]);
    //
    //                R = Compute2DRotationMatrix(gamma);
    //
    //                //              assert(rotatedCouplTensor.GetNumRows() == 2);
    //                //              assert(rotatedCouplTensor.GetNumCols() == 3);
    //                scaledCouplTensor.PerformRotation(R,rotatedCouplTensor);
    //              }
    //            }
    //          } else {
    //            // no rotation > take scale tensor
    //            rotatedCouplTensor = scaledCouplTensor;
    //          }
    //        }
    //        rotatedCouplingTensor_requiresReeval_[LPMit->first] = false;
    //        rotatedCouplingTensor_[LPMit->first] = rotatedCouplTensor;
    //        //        lastUsedTimeLevelForRotation_[storageIdx] = timeLevel;
    //      }
#pragma omp barrier
    }
#endif

  }

	Vector<Double> CoefFunctionHyst::CalcOutputOfHysteresisOperator(Vector<Double> inputToHystOperator, UInt operatorIdx, UInt storageIdx,
          bool forceMemoryLock, bool forceMemoryWrite, bool useOperatorForStrain) {

//    std::cout << "CalcOutputOfHysteresisOperator" << std::endl;
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
			LOG_TRACE(coeffunctionhyst_main) << "Calc output: hystOperatorLocked_ locked";
      			std::cout << "hystOperatorLocked_ locked" << std::endl;
      if(useOperatorForStrain){
        return P_J_forStrains_[storageIdx];
      } else {
        return P_J_[storageIdx];
      }
    }

    /*
     * NEW 18.4.2020:
     * do not check for diff in case of FP-H;
     * reason: in FP-H we are not allowed to set P_J_ and E_H_ during linesearches (see further below);
     * in this case diff might become zero if the same stepping is tested twice in a row; this would
     * return P_J_ which then would be the stored stated from after the previous outer iteration;
     * in short: the values cannot fit in this case!
     */
    
		if ((!useFPH_) && ((diff.NormL2() < POL_operatorParams_.amplitudeResolution_)&&(overwriteMemory == false))) {
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
//    std::cout << "Evaluate Hysteresis operator" << std::endl;
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
      bool debugOutput = false;

      if(useOperatorForStrain){
        outputOfHystOperator = hystStrain_->computeValue_vec(inputToHystOperator, operatorIdx, overwriteMemory,
                debugOutput, successFlag);
      } else {
        outputOfHystOperator = hyst_->computeValue_vec(inputToHystOperator, operatorIdx, overwriteMemory,
                debugOutput, successFlag);
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
            filenamebuf << "Switch_Elem" << firstIdx << "_Step" << std::setfill('0') << std::setw(5) << cnt << "_v" << POL_operatorParams_.evalVersion_ << "_numRows" << POL_weightParams_.numRows_ << ".bmp";
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
    bool miniOutput = false;
//    if(storageIdx == 1){
//      miniOutput = true;
//    }
    /*
     * NEW 18.4.2020:
     * unlike all B-based version, the H-versions (FixpointH_Global and FixpointH_Localized)
     * utilize E_H_ and P_J_ for actual computations during ALL evaluations;
     * this means, that we get an unwanted hysteretic behavior into the evaluation process;
     * this is especially problematic during linesearch, where multiple different steplengthes are tested;  
     * tests showed for example, that a stepping sequence 1, 0.5, 0.1, 1, 0.5, 0.1 lead to different residuals
     * for the first and second time although the solution vectors and the system were reset correctly after each trial
     * 
     * as mentioned, this came from the fact that the intermediate steps overwrote the storage and thus influenced the
     * eastimated input computed by the H-version; 
     * 
     * solution idea: add flage skipStorage_ that prohibits the storage from being overwritten; it is important however
     * that this flag is set to false at the end of each outer iteration (so after a valid steplength has been obtained)
     * as otherwise the H-version will not proceed
     * 
     * suggested treatment:
     * initialize flag with zero and keep flag at false until linesearch starts; set it to true in case of H-version;
     * after linesearch set flag to false again
     */
    if(skipStorage_ && useFPH_){
      if(miniOutput){
        std::cout << "Don NOT store E_H_ and P_J_!" << std::endl;
      }
    } else {
      if(useOperatorForStrain){
        P_J_forStrains_[storageIdx] = outputOfHystOperator;
      } else {
        if(miniOutput){
          std::cout << "Overwrite memory? " << overwriteMemory << std::endl;
          std::cout << "Store input/output combination for later usage" << std::endl;
          std::cout << "E_H_["<<storageIdx<<"] before setting: " << E_H_[storageIdx].ToString() << std::endl;
          std::cout << "P_J_[storageIdx] before setting: " << P_J_[storageIdx].ToString() << std::endl;
        }    
        E_H_[storageIdx] = inputToHystOperator;

        P_J_[storageIdx] = outputOfHystOperator;
        if(miniOutput){
          std::cout << "E_H_["<<storageIdx<<"] after setting: " << E_H_[storageIdx].ToString() << std::endl;
          std::cout << "P_J_[storageIdx] after setting: " << P_J_[storageIdx].ToString() << std::endl;
        }
      }
    }
    //std::cout << "P_J_[storageIdx] after setting: " << P_J_[storageIdx].ToString() << std::endl;

    //std::cout << "Computed output: " << outputOfHystOperator.ToString() << std::endl;
		// this flag has to be reset by hand!
		requiresReeval_[storageIdx] = false;

		return outputOfHystOperator;

	}

  void CoefFunctionHyst::ForceRemanence(){
		LOG_TRACE(coeffunctionhyst_main) << "ForceRemanence";
    /*
     * Iterate over all operators and evaluate them at 0
     * > at this time also set the last input to the hyst operator and set arrays
     */
    Vector<Double> output;
    Vector<Double> zeroVec = Vector<Double>(dim_);
    zeroVec.Init();

    bool forceMemoryLock = true;
    UInt operatorIdx, storageIdx,idxElem;

    LocPointMapped actualLPM;
    std::map<UInt, LocPointMapped>::iterator listIt;
    std::map<UInt, LocPointMapped> currentMap;
		currentMap = allLPMmap_;

    for(listIt = currentMap.begin(); listIt!= currentMap.end(); ++listIt){
      LOG_TRACE(coeffunctionhyst_main) << "Storage idx: " << listIt->first;
      bool forceMidpoint;

      if(midpointLPMmap_.find(listIt->first) != midpointLPMmap_.end()){
        forceMidpoint = true;
      } else {
        forceMidpoint = false;
      }

      PreprocessLPM(listIt->second, actualLPM, operatorIdx, storageIdx, idxElem,forceMidpoint);
      if(listIt->first != storageIdx){
        EXCEPTION("StorageIdx in map != storageIdx as obtained by PreprocessLPM!");
      }

      E_H_[storageIdx] = zeroVec;

      // then compute output
      output = CalcOutputOfHysteresisOperator(zeroVec, operatorIdx, storageIdx, forceMemoryLock, false);
			LOG_TRACE(coeffunctionhyst_main) << "Computed output: " << output.ToString();

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

  void CoefFunctionHyst::DetermineChunkSizeAndThreads(UInt size, UInt& chunksize, UInt& numThreads){
    /*
     * little helper to get a proper chunksize and number of threads for parallel computations with
     * openMP
     */
    UInt minChunksize = 2; // let each thread do at least some points
    numThreads = CFS_NUM_THREADS;

    if(minChunksize >= size){
      numThreads = 1;
      chunksize = size;
    } else {
      chunksize = std::floor(size/numThreads);
      //    std::cout << "Total number of operators to evaluate " << size << std::endl;
      if(chunksize < minChunksize){
        //      std::cout << "chunksize = " << chunksize << " < minimal chunksize = " << minChunksize << std::endl;
        numThreads = std::floor(size/minChunksize);
        chunksize = std::floor(size/numThreads);
      }
    }
  }

  void CoefFunctionHyst::PrecomputePolarization(bool overwriteMemory){
//    std::cout << "PrecomputePolarization" << std::endl;
    UInt RUN_evaluationPurpose_bak = RUN_evaluationPurpose_;
    RUN_evaluationPurpose_ = 1;
    UInt size = operatorToStorage_.size();
    UInt chunksize;
    UInt numT;
    DetermineChunkSizeAndThreads(size, chunksize, numT);

#ifdef USE_OPENMP
#pragma omp parallel num_threads(numT)
    {
      std::map<UInt, std::map<UInt, bool> >::iterator itStart = operatorToStorage_.begin();
      std::map<UInt, std::map<UInt, bool> >::iterator itEnd;
      std::map<UInt, std::map<UInt, bool> >::iterator it;
      //std::cout << "Running parallel" << std::endl;

      UInt aThread = omp_get_thread_num();
      std::advance(itStart, aThread*chunksize);

      if(aThread == numT-1){
        itEnd = operatorToStorage_.end();
      } else {
        itEnd = itStart;
        std::advance(itEnd, chunksize);
      }
//      #pragma omp critical(output)
//      {
//        std::cout << "My number: " << aThread << "; my workload range: "<< itStart->first << " to " << itEnd->first  <<std::endl;
//      }
//      #pragma omp barrier
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

      // barrier not required
//#pragma omp barrier
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

          UInt operatorIdxTMP, storageIdxTMP,idxElemTMP;
          onBoundary = PreprocessLPM(curLPM, actualLPM, operatorIdxTMP, storageIdxTMP, idxElemTMP,isMidpoint);

          if( (operatorIdxTMP != operatorIdx) || (storageIdxTMP != storageIdx)){
            EXCEPTION("Storage or operator index does not fit");
          }

//          std::cout << "My number: " << aThread << "; Working on operator: " << operatorIdx << "; on storage: " << storageIdx << std::endl;

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
//#pragma omp barrier
    } // end parallel section

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

        UInt operatorIdxTMP, storageIdxTMP,idxElemTMP;
        onBoundary = PreprocessLPM(curLPM, actualLPM, operatorIdxTMP, storageIdxTMP, idxElemTMP, isMidpoint);

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

  void CoefFunctionHyst::RestorePrevHystVals(UInt storageLocation) {
    //    std::cout << "RestorePrevHystVals" << std::endl;
    /*
     * NEW: just copy values from current array to lastIT or lastTS array
     * Addition 27.11.2018:
     *    three different storage locations
     *    0 = previousTimestep
     *    1 = previousIteration
     *    2 = estimation
     */
    std::map<UInt, std::map<UInt, bool> >::iterator itStart = operatorToStorage_.begin();
    std::map<UInt, std::map<UInt, bool> >::iterator itEnd = operatorToStorage_.end();
    std::map<UInt, std::map<UInt, bool> >::iterator it;

    for(it = itStart; it != itEnd; it++){
      std::map<UInt, bool>::iterator innerIt;
      for(innerIt = it->second.begin(); innerIt != it->second.end(); innerIt++){
        UInt storageIdx = innerIt->first;

        if (storageLocation == 1) {
          if(lastTSStorageInitialized_ == false){
            InitSpecificStorage(1);
          }
          // restore E_B_, P_J_, E_H_ and Si_ from prevous timestep
          E_B_[storageIdx] = E_B_lastTS_[storageIdx];
          P_J_[storageIdx] = P_J_lastTS_[storageIdx];
          P_J_forStrains_[storageIdx]= P_J_forStrains_lastTS_[storageIdx];
          E_H_[storageIdx] = E_H_lastTS_[storageIdx];
          if(strainRequired_){
            Si_[storageIdx] = Si_lastTS_[storageIdx];
          }
        } else if(storageLocation == 0){
          if(lastItStorageInitialized_ == false){
            InitSpecificStorage(0);
          }
          // restore E_B_, P_J_, E_H_ and Si_ from prevous timestep
          E_B_[storageIdx] = E_B_lastIt_[storageIdx];
          P_J_[storageIdx] = P_J_lastIt_[storageIdx];
          E_H_[storageIdx] = E_H_lastIt_[storageIdx];
          if(strainRequired_){
            Si_[storageIdx] = Si_lastIt_[storageIdx];
          }
        } else if(storageLocation == 3){
//          std::cout << "Restore from backup" << std::endl;
          if(backupStorageInitialized_ == false){
            InitSpecificStorage(3);
          }
          // restore E_B_, P_J_, E_H_ and Si_ from prevous timestep
          E_B_[storageIdx] = E_B_BAK_[storageIdx];
          P_J_[storageIdx] = P_J_BAK_[storageIdx];
          E_H_[storageIdx] = E_H_BAK_[storageIdx];
          if(strainRequired_){
            Si_[storageIdx] = Si_BAK_[storageIdx];
          }
        } else if(storageLocation == 2){
          WARN("Estimator step not meant to be used for restoring");
        } else {
          EXCEPTION("storageLocation must be 0 (lastIt) or 1 (lastTS)");
        }

        if(storageIdx == 1){
          LOG_DBG2(coeffunctionhyst_main) << "RESTORED PREV HYST VALUE";
          LOG_DBG2(coeffunctionhyst_main) << "storage location: " << storageLocation;
          LOG_DBG2(coeffunctionhyst_main) << "P_J_[storageIdx]: " << P_J_[storageIdx].ToString();
          LOG_DBG2(coeffunctionhyst_main) << "E_B_[storageIdx]: " << E_B_[storageIdx].ToString();
          LOG_DBG2(coeffunctionhyst_main) << "E_H_[storageIdx]: " << E_H_[storageIdx].ToString();
        }
      }
    }
  }

  void CoefFunctionHyst::SetPreviousHystVals(UInt storageLocation) {
    /*
     * NEW: just copy values from current array to lastIT or lastTS array
     * Addition 27.11.2018:
     *    three different storage locations
     *    0 = previousTimestep
     *    1 = previousIteration
     *    2 = estimation
     */
    std::map<UInt, std::map<UInt, bool> >::iterator itStart = operatorToStorage_.begin();
    std::map<UInt, std::map<UInt, bool> >::iterator itEnd = operatorToStorage_.end();
    std::map<UInt, std::map<UInt, bool> >::iterator it;

    for(it = itStart; it != itEnd; it++){
      std::map<UInt, bool>::iterator innerIt;
      for(innerIt = it->second.begin(); innerIt != it->second.end(); innerIt++){
        UInt storageIdx = innerIt->first;

        if(storageIdx == 1){
          LOG_DBG(coeffunctionhyst_main) << "SET PREV HYST VALUE";
          LOG_DBG(coeffunctionhyst_main) << "storage location: " << storageLocation;
          LOG_DBG(coeffunctionhyst_main) << "P_J_[storageIdx]: " << P_J_[storageIdx].ToString();
          LOG_DBG(coeffunctionhyst_main) << "E_B_[storageIdx]: " << E_B_[storageIdx].ToString();
          LOG_DBG(coeffunctionhyst_main) << "E_H_[storageIdx]: " << E_H_[storageIdx].ToString();
        }

        if (storageLocation == 1) {
          if(lastTSStorageInitialized_ == false){
            InitSpecificStorage(1);
          }
          E_B_lastTS_[storageIdx] = E_B_[storageIdx];
          P_J_lastTS_[storageIdx] = P_J_[storageIdx];
          E_H_lastTS_[storageIdx] = E_H_[storageIdx];
          if(strainRequired_){
            Si_lastTS_[storageIdx] = Si_[storageIdx];
            P_J_forStrains_lastTS_[storageIdx] = P_J_forStrains_[storageIdx];
          }
        } else if(storageLocation == 0){
          if(lastItStorageInitialized_ == false){
            InitSpecificStorage(0);
          }
          E_B_lastIt_[storageIdx] = E_B_[storageIdx];
          P_J_lastIt_[storageIdx] = P_J_[storageIdx];
          E_H_lastIt_[storageIdx] = E_H_[storageIdx];
          if(strainRequired_){
            Si_lastIt_[storageIdx] = Si_[storageIdx];
          }
        } else if(storageLocation == 2){
          if(estimatorStorageInitialized_ == false){
            InitSpecificStorage(2);
          }
          estimatorSet_ = true;
          E_B_estimatorStep_[storageIdx] = E_B_[storageIdx];
          P_J_estimatorStep_[storageIdx] = P_J_[storageIdx];
          E_H_estimatorStep_[storageIdx] = E_H_[storageIdx];
          if(strainRequired_){
            Si_estimatorStep_[storageIdx] = Si_[storageIdx];
          }
        } else if(storageLocation == 3){
//          std::cout << "Store to backup" << std::endl;
          if(backupStorageInitialized_ == false){
            InitSpecificStorage(3);
          }
          E_B_BAK_[storageIdx] = E_B_[storageIdx];
          P_J_BAK_[storageIdx] = P_J_[storageIdx];
          E_H_BAK_[storageIdx] = E_H_[storageIdx];
          if(strainRequired_){
            Si_BAK_[storageIdx] = Si_[storageIdx];
          }
        } else {
          EXCEPTION("storageLocation must be 0 (lastIt) or 1 (lastTS) or 2 (estimator");
        }
        deltaMatPrev_[storageIdx] = deltaMat_[storageIdx];
        if(strainRequired_){
          deltaMatStrainPrev_[storageIdx] = deltaMatStrain_[storageIdx];
        }
      }
    }
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

  void CoefFunctionHyst::SetDoubleFlag(std::string flagName, Double doubleState){
    if(flagName == "FPGlobalC"){
//      std::cout << "Set globalFPFactor_C to " << doubleState << std::endl;
      globalFPFactor_C = doubleState;
    }
    if(flagName == "FPLocalC"){
//      std::cout << "Set localFPFactor_C to " << doubleState << std::endl;
      localFPFactor_C = doubleState;
    }
  }

  void CoefFunctionHyst::SetFlag(std::string flagName, Integer intState) {

    /*
     * Flags that (may) change during runtime
     */
    if(flagName == "evaluationPurpose"){
      RUN_evaluationPurpose_ = intState;
    }
    else if(flagName == "evalJacAtMidpointOnly"){
      evalJacAtMidpointOnly_ = bool(intState);
    }
    else if(flagName == "allowBMP"){
      RUN_allowBMP_ = bool(intState);
    }
    else if(flagName == "debugOutput"){
      DEBUG_debugOutput_ = bool(intState);
    }
    else if(flagName == "infoOutput"){
      DEBUG_infoOutput_ = bool(intState);
    }
    else if(flagName == "printCallingOrder"){
      DEBUG_printCallingOrder_ = bool(intState);
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
//      std::cout << "forceCurrentTS was set to " << bool(intState) << std::endl;
      forceCurrentTS_ = bool(intState);
    }
    else if(flagName == "PrintTimeLevels"){
      std::cout << "timeLevel_Mat_ = " << timeLevel_Mat_ << std::endl;
      std::cout << "timeLevel_rhsPol_ = " << timeLevel_rhsPol_ << std::endl;
      std::cout << "timeLevel_rhsStrain_ = " << timeLevel_rhsStrain_ << std::endl;
      std::cout << "timeLevel_rhsCoupling_ = " << timeLevel_rhsCoupling_ << std::endl;
      std::cout << "timeLevel_deltaMatPol_ = " << timeLevel_deltaMatPol_ << std::endl;
      std::cout << "timeLevel_deltaMatStrain_ = " << timeLevel_deltaMatStrain_ << std::endl;
      std::cout << "timeLevel_deltaMatCoupling_ = " << timeLevel_deltaMatCoupling_ << std::endl;
      std::cout << "timeLevel_results_ = " << timeLevel_results_ << std::endl;
      std::cout << "timeLevel_boundaryTerms_ = " << timeLevel_boundaryTerms_ << std::endl;
    } else if (flagName == "TestJacobianApproximations"){
      TestJacobianApproximations();
    } else if (flagName == "SetFPH"){
      useFPH_ = bool(intState);
    } else if (flagName == "SkipStorage"){
      skipStorage_ = bool(intState);
    } else if (flagName == "FPApproach"){
      selectedFPApproach_ = fixpointFlag(intState);
//      std::cout << "Selected FPApproach: " << selectedFPApproach_ << std::endl;
    } else if( flagName == "EvalFPMaterialTensors"){
//      SetFPMaterialTensorsNEW(intState);
      SetFPMaterialTensorsTMP(intState);
    }
    else if( flagName == "SetFPCorrectionFlag"){
//      std::cout << "Set FP RHS Correction to " << intState << std::endl;
      FPMaterialRHSSet_ = UInt(intState);
    }
    else if( flagName == "SetFPMatrixFlag"){
//      std::cout << "Set FP System Matrix Flag to " << intState << std::endl;

      if((intState == 1)||(intState == 3)){
        // both global FP method (intState == 1) and H-Version (intState = 3) require information about
        // the maximal and minimal slope of the hysteresis curves; to estimate these quantities, we have to trace
        // the hysteresis model once
        if(!hystModelTraced_){
          UInt baseSteps = 100;
          // get some additinal info if we are already traverse loop
          Double negCoercivity = 0.0;
          Double maxPolarization = 0.0;

//          if (POL_operatorParams_.methodName_ == "scalarPreisach") {
//            TraceHystOperator(baseSteps, maxSlopeGlobal_, minSlopeGlobal_, negCoercivity, maxPolarization);
//          } else {
            TraceHystOperatorVector(baseSteps, maxSlopeGlobal_, minSlopeGlobal_, negCoercivity, maxPolarization);
//          }
          hystModelTraced_ = true;
        }
      }
      FPMaterialTensorSet_ = UInt(intState);
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

//      std::cout << "timeLevel_deltaMatPol_ was set to " << timeLevel_deltaMatPol_ << std::endl;

      // -2 means that deltaMat is completely deactivated
      // 0 would mean that this specific deltaMat is deactivated but another one might be active
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
    else if(flagName == "DeactivateDeltaMat"){
      deltaMatActive_ = false;
    }
    else if(flagName == "ActivateDeltaMat"){
      deltaMatActive_ = true;
    }
    else if(flagName == "resetEstimator"){
      estimatorSet_ = false;
    }
    else if(flagName == "resetGlobalFPTensors"){
      globalFPset_ = false;
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
    else if(flagName == "resetReevalCouplingTensor"){
      for (UInt i = 0; i < numStorageEntries_; i++) {
        rotatedCouplingTensor_requiresReeval_[i] = true;
      }
    }
    else if(flagName == "allowSettingOfMatForLocalInversion"){
      // reset flag for matrix inversion, i.e. the next time we call SetMatrixForInversion, we will overwrite the current value
      for (UInt i = 0; i < numStorageEntries_; i++) {
        eps_mu_local_Set_[i] = false;
      }
      //      // to set matrices for inversion, flag has to be 1
      //      EvaluateHystOperatorsInt(1);

    } else if(flagName == "EvaluateHystOperators"){
      EvaluateHystOperatorsInt(intState);
    } else if(flagName == "SetPreviousHystValues"){
      SetPreviousHystVals(intState);
    } else if(flagName == "RestorePreviousHystValues"){
      RestorePrevHystVals(intState);
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

  void CoefFunctionHyst::CreatePeriodicTestSignal(std::string name, Double amplitudeScaling, Double numPeriods, UInt stepsPerPeriod, Vector<Double>& xVals, Vector<Double>& yVals, Double initialAmplitude){

    UInt initialSteps = 0;
    if(initialAmplitude != std::numeric_limits<double>::infinity()){
      initialSteps = 1;
    }
    UInt totalSteps = UInt(stepsPerPeriod*numPeriods+initialSteps);
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
      if(initialSteps != 0){
        xVals[0] = initialAmplitude/amplitudeScaling*POL_operatorParams_.inputSat_*sin( (2*M_PI*0)/stepsPerPeriod );
        yVals[0] = initialAmplitude/amplitudeScaling*POL_operatorParams_.inputSat_*cos( (2*M_PI*0)/stepsPerPeriod );
      }
      for(UInt i = initialSteps; i < totalSteps; i++){
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

  void CoefFunctionHyst::WriteSignalToFile(std::string name, Vector<Double> xVals, Vector<Double> yVals, Vector<Double> zVals){
    std::ofstream output;
    output.open(name);
    UInt totalSteps = xVals.GetSize();
    output << "# Step    x    y    z" << std::endl;
    for(UInt i = 0; i < totalSteps; i++){
      output << i << " " << xVals[i] << " " << yVals[i] << " " << zVals[i] << std::endl;
    }
    output.close();
  }

  void CoefFunctionHyst::TestHystOperatorWithSignal(std::string name, Vector<Double> xVals, Vector<Double> yVals,
          Vector<Double> zVals, UInt dimHystOperator, int forcedPreisachResolution,
          bool testInversion, bool printStatistics, bool writeResultsToFile,
          bool measurePerformance, std::string commonPerformanceFile, bool test1D, bool outputIrrStrains,
          std::string nameTagForPerfFile,std::string additionalTag1,std::string additionalTag2,std::string additionalTag3){

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
    std::fstream commonPerfStream;

    std::stringstream basedir_name;
    std::stringstream forcedResolution_name;
    
    basedir_name << targetDirForResultFiles_;
    if(targetDirForResultFiles_.back() != '/'){
      basedir_name << "/";
    }
    //basedir_name << "./history_hystOperator/";
    
    try {
      fs::create_directory( basedir_name.str() );
    } catch (std::exception &ex) {
      EXCEPTION(ex.what());
    }
    
    forcedResolution_name << "";
//    std::string basedir = "./history_hystOperator/";

    std::string stepDirBack = "";
    if( forcedPreisachResolution > 0 ){
//      basedir_name << "forcedResolution-"<<forcedPreisachResolution<<"/";
//      stepDirBack = "../";
      forcedResolution_name << "-forcedRes-"<<forcedPreisachResolution;
    } else {
      forcedResolution_name << "-definedRes-"<<POL_weightParams_.numRows_;
    }
    std::string basedir = basedir_name.str();
    std::string forcedResolutionString = forcedResolution_name.str();

    try {
      fs::create_directory( basedir );
    } catch (std::exception &ex) {
      EXCEPTION(ex.what());
    }

		std::stringstream statistics_name;
		statistics_name << basedir << material_->GetName() << forcedResolutionString << name << "_statistics.txt";

		std::stringstream results_name_x;
		results_name_x << basedir << material_->GetName() << forcedResolutionString << name << "_results_x.txt";

    std::stringstream results_name_p;
		results_name_p << basedir << material_->GetName() << forcedResolutionString << name << "_results_p.txt";

    std::stringstream results_name_s;
		results_name_s << basedir << material_->GetName() << forcedResolutionString << name << "_results_s.txt";

    std::stringstream results_name_y;
		results_name_y << basedir << material_->GetName() << forcedResolutionString << name << "_results_y.txt";

    std::stringstream results_name_xp;
		results_name_xp << basedir << material_->GetName() << forcedResolutionString << name << "_results_xp.txt";

    std::stringstream results_name_xps;
		results_name_xps << basedir << material_->GetName() << forcedResolutionString << name << "_results_xps.txt";

    std::stringstream angularResults_name_x;
		angularResults_name_x << basedir << material_->GetName() << forcedResolutionString << name << "_angularResults_x.txt";

    std::stringstream angularResults_name_p;
		angularResults_name_p << basedir << material_->GetName() << forcedResolutionString << name << "_angularResults_p.txt";

    std::stringstream angularResults_name_y;
		angularResults_name_y << basedir << material_->GetName() << forcedResolutionString << name << "_angularResults_y.txt";

		std::stringstream performance_name;
		performance_name << basedir << material_->GetName() << forcedResolutionString << name << "_performance.txt";

    std::stringstream orientation_name;
		orientation_name << basedir << material_->GetName() << forcedResolutionString << name << "_projectedCoords_p";

    bool useOldNameTag=false;
    std::stringstream nameTagForPerfFile_name;
    nameTagForPerfFile_name << nameTagForPerfFile;
    if(useOldNameTag){
      std::stringstream nameTagForPerfFile_name;
      nameTagForPerfFile_name << nameTagForPerfFile << "-" << material_->GetName() << forcedResolutionString;
    }else{
      if(additionalTag1 != "---"){
        nameTagForPerfFile_name << "\t" << additionalTag1;
      }
      if(additionalTag2 != "---"){
        nameTagForPerfFile_name << "\t" << additionalTag2;
      }
      if(additionalTag3 != "---"){
        nameTagForPerfFile_name << "\t" << additionalTag3;
      }
    }
    nameTagForPerfFile = nameTagForPerfFile_name.str();
		        
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

//    writeResultsToInfoXML_
    
		if(printStatistics){
			statistics.open(statistics_name.str());
		}
    
    bool addHeaderToCommonPerFile = false;
    std::stringstream commonPerfStream_name;
		if(measurePerformance){
			performance.open(performance_name.str());
      if(commonPerformanceFile != "---"){
        commonPerfStream_name << basedir << stepDirBack << commonPerformanceFile;
        
        // read file first to check if it is empty
        commonPerfStream.open(commonPerfStream_name.str(),std::ios_base::in);
        if(commonPerfStream.is_open()){
          // file could be opened for reading, thus it must already exist
          if(commonPerfStream.peek() == std::ifstream::traits_type::eof()){
            // file is empty though > add header
            addHeaderToCommonPerFile = true;
//            std::cout << "Common performance file could be opened but is empty; add header" << std::endl;
          }
          commonPerfStream.close();
        } else {
          // file does not exist yet, opening it with std::ios_base::app should create the file
          // add header to new file
          addHeaderToCommonPerFile = true;
//          std::cout << "Common performance file could not be opened for reading" << std::endl;
        }
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
    //  3-6 only for vector jacobianImplementation using Levenberg Marquardt, Newton or fixpoint iteration
    //  3 = remanence
    //  4 = passed dut to error tolerance
    //  5 = passed due to tolerance wrt x
    //  6 = passed due to tolerance wrt y
    //  7 = passed fixpoint due to tolerance wrt to x
    //  8 = passed fixpoint due to tolerance wrt to y
		Vector<Double> xIn = Vector<Double>(dimHystOperator);
		Vector<Double> hIn = Vector<Double>(dimHystOperator);
		Vector<Double> yIn = Vector<Double>(dimHystOperator);

    Vector<Double> xInPrev = Vector<Double>(dimHystOperator);
		Vector<Double> xPrev = Vector<Double>(dimHystOperator);
		Vector<Double> hPrev = Vector<Double>(dimHystOperator);
		Vector<Double> yPrev = Vector<Double>(dimHystOperator);

		Vector<Double> xOut = Vector<Double>(dimHystOperator);
		Vector<Double> hOut = Vector<Double>(dimHystOperator);
    Vector<Double> hOutForStrains = Vector<Double>(dimHystOperator);
		Vector<Double> yOut = Vector<Double>(dimHystOperator);
		Vector<Double> hOutOld = Vector<Double>(dimHystOperator);
    Vector<Double> hOutOldForStrains = Vector<Double>(dimHystOperator);

    Vector<Double> xInBak = Vector<Double>(dimHystOperator);
		Vector<Double> xErr = Vector<Double>(dimHystOperator);
		Vector<Double> yErr = Vector<Double>(dimHystOperator);

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
		UInt forwardReused = 0;
		UInt forwardAnhystOnly = 0;
		UInt forwardOverwrite = 0;
		UInt forwardTMP = 0;

		// counter for backward evaluation / inversion
		UInt totalFails = 0;
		UInt totalReused = 0;
		UInt totalAnhystOnly = 0;
		UInt totalBisection = 0;
    UInt totalFixedPoint = 0;
    UInt totalNewtonBased = 0;
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
    UInt backwardReuseCounter = 0;
		Timer* forwardTimer = NULL;
		Timer* backwardTimer = NULL;
    hystTMP = NULL;
    hystStrainTMP = NULL;
    
    UInt precisionDigits = 12;//9
    string separatorString = " ";
    UInt minWidthOutputColumns = 10;

    std::string usedInversionMethod;
    Enum2String(InversionParams_.inversionMethod,usedInversionMethod);
    
		if(printStatistics){
			statistics << "+++ STATISTICS +++" << std::endl;
			statistics << "TEST: " << name << std::endl;
      statistics << "HYSTERESIS MODEL: " << POL_operatorParams_.methodName_ << std::endl;
      statistics << "LOCAL INVERSION METHOD: " << usedInversionMethod << std::endl;
      if(InversionParams_.inversionMethod == LOCAL_FIXPOINT){
        statistics << "FP-Contraction/Safety factor: " << InversionParams_.safetyFactor_C << std::endl;
        statistics << "Min/max slope from tracing: " << InversionParams_.minSlopeHystOperator << " / " << InversionParams_.maxSlopeHystOperator << std::endl;
      }      
      statistics << "Error Criteria: " << std::endl;
			statistics << "- Residual wrt input (=x): " << InversionParams_.tolH << std::endl;
      statistics << "-> Relative? " << InversionParams_.tolH_useAsRelativeNorm << std::endl;
      statistics << "- Residual wrt output (=y): " << InversionParams_.tolB << std::endl;
      statistics << "-> Relative? " << InversionParams_.tolB_useAsRelativeNorm << std::endl;
      statistics << "Maximal number of outer iterations: " << InversionParams_.maxNumIts << std::endl;
    }

		/*
     * 1. Create temporal hyst operator; this should be done for each test to ensure
     *		that operator is in initial state
     */
		bool vector = true;
		bool debugOut = false;
    bool forStrains = false;
    
    // used by trace hyst operator to guarantess that scalar model and tracing signal align (along x-axis)
    // not required here
    bool forceScalarDirection = false; 
    
    // resize
    if(dimHystOperator > dim_){
      POL_operatorParams_.fixDirection_.Resize(dimHystOperator);
      POL_operatorParams_.startingAxisMG_.Resize(dimHystOperator);
    }
    
    hystTMP = CoefFunctionHyst::getTemporalHystOperator(forStrains,forceScalarDirection,dimHystOperator,forcedPreisachResolution);
        
		if (POL_operatorParams_.methodName_ == "scalarPreisach") {
      vector = false;
      test1D = true;
      if(printStatistics){
        statistics << "- orientation of scalar model: " << POL_operatorParams_.fixDirection_.ToString() << std::endl;
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
		} else if (POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") {
      if(printStatistics){
        statistics << "with num directions = " << POL_operatorParams_.numDirections_ << std::endl;
        statistics << "and output clipping mode = " << POL_operatorParams_.outputClipping_ << std::endl;
        statistics << "Starting axis: " << POL_operatorParams_.startingAxisMG_.ToString() << std::endl;
      }
		} else {
			EXCEPTION("Invalid model selected for inversion test");
		}
    if((printStatistics)&&(vector)){
      statistics << "- dimension of vector model = " << dimHystOperator << std::endl;
    }
    
    if(CouplingParams_.ownHystOperator_ && CouplingParams_.couplingDefined_inMatFile_){
      if(printStatistics){
        statistics << "Separate HYSTERESIS MODEL for STRAINS: " << STRAIN_operatorParams_.methodName_ << std::endl;
      }
      forStrains = true;
      hystStrainTMP = CoefFunctionHyst::getTemporalHystOperator(forStrains,forceScalarDirection,dimHystOperator,forcedPreisachResolution);
      
      if (STRAIN_operatorParams_.methodName_ == "scalarPreisach") {
        if(printStatistics){
          statistics << "- orientation of scalar model: " << STRAIN_operatorParams_.fixDirection_.ToString() << std::endl;
        }
      } else if (STRAIN_operatorParams_.methodName_ == "vectorPreisach_Sutor") {
        if(printStatistics){
          if (STRAIN_operatorParams_.evalVersion_ == 1) {
            statistics << "- classical model, list based implementation" << std::endl;
          } else if (STRAIN_operatorParams_.evalVersion_ == 2) {
            statistics << "- revised model, list based implementation" << std::endl;
          } else if (STRAIN_operatorParams_.evalVersion_ == 10) {
            statistics << "- classical model, matrix based implementation" << std::endl;
          } else if (STRAIN_operatorParams_.evalVersion_ == 20) {
            statistics << "- revised model, matrix based implementation" << std::endl;
          }
          statistics << "with rotational resistance = " << STRAIN_operatorParams_.rotResistance_ << std::endl;
          if(STRAIN_operatorParams_.isClassical_ == false){
            statistics << "and angular distance = " << STRAIN_operatorParams_.angularDistance_ << std::endl;
          }
        }
      } else if (STRAIN_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") {
        if(printStatistics){
          statistics << "with num directions = " << STRAIN_operatorParams_.numDirections_ << std::endl;
          statistics << "and output clipping mode = " << STRAIN_operatorParams_.outputClipping_ << std::endl;
          statistics << "Starting axis: " << STRAIN_operatorParams_.startingAxisMG_.ToString() << std::endl;
        }
      } else {
        EXCEPTION("Invalid model selected for inversion test");
      }
    }

    if(outputIrrStrains){
      if(dimHystOperator == 2){
        results_s.open(results_name_s.str());
        results_s << "#Number S_irr_xx S_irr_yy S_irr_xy" << std::endl;

        results_xps.open(results_name_xps.str());
        results_xps << "#Number xIn yIn Px Py S_irr_xx S_irr_yy S_irr_xy" << std::endl;
      } else {
        results_s.open(results_name_s.str());
        results_s << "#Number S_irr_xx S_irr_yy S_irr_zz S_irr_yz S_irr_xz S_irr_xy" << std::endl;

        results_xps.open(results_name_xps.str());
        results_xps << "#Number xIn yIn zIn Px Py Pz S_irr_xx S_irr_yy S_irr_zz S_irr_yz S_irr_xz S_irr_xy" << std::endl;
      }
    }

		if(printStatistics){
			statistics << "PARAMETER: " << std::endl;
			statistics << "- xSAT: " << POL_operatorParams_.inputSat_ << std::endl;
			statistics << "- pSAT: " << POL_operatorParams_.outputSat_ << std::endl;
			statistics << "- anhyst a: " << POL_weightParams_.anhysteretic_a_ << std::endl;
			statistics << "- anhyst b: " << POL_weightParams_.anhysteretic_b_ << std::endl;
			statistics << "- anhyst c: " << POL_weightParams_.anhysteretic_c_ << std::endl;
      statistics << "- anhyst d: " << POL_weightParams_.anhysteretic_d_ << std::endl;
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
      } else if(POL_weightParams_.weightType_ == "muDatExtended"){
        statistics << "> extended muDat with " << std::endl;
        statistics << "- A = " << POL_weightParams_.muDat_A_ << std::endl;
        statistics << "- sigma1 = " << POL_weightParams_.muDat_sigma1_ << std::endl;
        statistics << "- sigma2 = " << POL_weightParams_.muDat_sigma2_ << std::endl;
        statistics << "- h1 = " << POL_weightParams_.muDat_h1_ << std::endl;
        statistics << "- h2 = " << POL_weightParams_.muDat_h2_ << std::endl;
        statistics << "- eta = " << POL_weightParams_.muDat_eta_ << std::endl;
      } else if(POL_weightParams_.weightType_ == "muLorentz"){
        statistics << "> Lorentzian weight function with " << std::endl;
        statistics << "- A = " << POL_weightParams_.muLorentz_A_ << std::endl;
        statistics << "- sigma = " << POL_weightParams_.muLorentz_sigma1_ << std::endl;
        statistics << "- h = " << POL_weightParams_.muLorentz_h1_ << std::endl;
      } else if(POL_weightParams_.weightType_ == "muLorentzExtended"){
        statistics << "> Lorentzian weight function with " << std::endl;
        statistics << "- A = " << POL_weightParams_.muLorentz_A_ << std::endl;
        statistics << "- sigma1 = " << POL_weightParams_.muLorentz_sigma1_ << std::endl;
        statistics << "- sigma2 = " << POL_weightParams_.muLorentz_sigma2_ << std::endl;
        statistics << "- h1 = " << POL_weightParams_.muLorentz_h1_ << std::endl;
        statistics << "- h2 = " << POL_weightParams_.muLorentz_h2_ << std::endl;
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
//				statistics << "- max number of iterations: " << InversionParams_.maxNumIts << std::endl;
//				statistics << "- tolerance wrt x: " << InversionParams_.tolH << std::endl;
//				statistics << "- tolerance wrt y: " << InversionParams_.tolB << std::endl;
//				statistics << "- FD-resolution for Jacobian: " << InversionParams_.jacRes << std::endl;
//				statistics << "LINESEARCH: " << std::endl;
//				statistics << "- alpha start: " << InversionParams_.alphaLSStart << std::endl;
//				statistics << "- alpha min: " << InversionParams_.alphaLSMin << std::endl;
//				statistics << "- alpha max: " << InversionParams_.alphaLSMax << std::endl;
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
        if(dimHystOperator == 2){
          results_x << "#Number xIn[0] xIn[1] xOut[0] xOut[1]" << std::endl;
          results_x << "# > xIn[0],xIn[1] = x and y components of input; given by test signal" << std::endl;
          results_x << "# > xOut[0],xOut[1] = x and y components of inversion output" << std::endl;

          results_p << "#Number pIn[0] pIn[1] pOut[0] pOut[1]" << std::endl;
          results_p << "# > pIn[0],pIn[1] = x and y components of hyst operator computed from xIn" << std::endl;
          results_p << "# > pOut[0],pOut[1] = x and y components of hyst operator computed from xOut" << std::endl;

          results_y << "#Number yIn[0] yIn[1] yOut[0] yOut[1]" << std::endl;
          results_y << "# > yIn[0],yIn[1] = x and y components of output of hyst operator to xIn; used as input for inversion" << std::endl;
          results_y << "# > yOut[0],yOut[1] = x and y components computed from xOut" << std::endl;

          angularResults_x << "#Number abs(xIn) atan2(xIn[1],xIn[0])*180/pi abs(xOut) atan2(xOut[1],xOut[0])*180/pi" << std::endl;
          angularResults_x << "# > xIn[0],xIn[1] = x and y components of input; given by test signal" << std::endl;
          angularResults_x << "# > xOut[0],xOut[1] = x and y components of inversion output" << std::endl;

          angularResults_p << "#Number abs(pIn) atan2(pIn[1],pIn[0])*180/pi abs(pOut) atan2(pOut[1],pOut[0])*180/pi" << std::endl;
          angularResults_p << "# > pIn[0],pIn[1] = x and y components of hyst operator computed from xIn" << std::endl;
          angularResults_p << "# > pOut[0],pOut[1] = x and y components of hyst operator computed from xOut" << std::endl;

          angularResults_y << "#Number abs(yIn) atan2(yIn[1],yIn[0])*180/pi abs(yOut) atan2(yOut[1],yOut[0])*180/pi" << std::endl;
          angularResults_y << "# > yIn[0],yIn[1] = x and y components of output of hyst operator to xIn; used as input for inversion" << std::endl;
          angularResults_y << "# > yOut[0],yOut[1] = x and y components computed from xOut" << std::endl;
        } else {
          results_x << "#Number xIn[0] xIn[1] xIn[2] xOut[0] xOut[1] xOut[2]" << std::endl;
          results_x << "# > xIn[0],xIn[1],xIn[2] = x,y and z components of input; given by test signal" << std::endl;
          results_x << "# > xOut[0],xOut[1],xOut[2] = x,y and z components of inversion output" << std::endl;

          results_p << "#Number pIn[0] pIn[1] pIn[2] pOut[0] pOut[1] pOut[2]" << std::endl;
          results_p << "# > pIn[0],pIn[1],pIn[2] = x,y and z components of hyst operator computed from xIn" << std::endl;
          results_p << "# > pOut[0],pOut[1],pOut[2] = x,y and z components of hyst operator computed from xOut" << std::endl;

          results_y << "#Number yIn[0] yIn[1] yIn[2] yOut[0] yOut[1] yOut[2]" << std::endl;
          results_y << "# > yIn[0],yIn[1],yIn[2] = x,y and z components of output of hyst operator to xIn; used as input for inversion" << std::endl;
          results_y << "# > yOut[0],yOut[1],yOut[2] = x,y and z components computed from xOut" << std::endl;

          angularResults_x << "#Number r_in = abs(xIn); phi_in = atan2(xIn[1],xIn[0])*180/pi; theta_in = atan2(r_in,xIn[2])*180/pi; r_out = abs(xOut); phi_out = atan2(xOut[1],xOut[0])*180/pi; theta_out = atan2(r_out,xOut[2])*180/pi" << std::endl;
          angularResults_x << "# > xIn[0],xIn[1],xIn[2] = x,y and z components of input; given by test signal" << std::endl;
          angularResults_x << "# > xOut[0],xOut[1],xOut[2] = x,y and z components of inversion output" << std::endl;

          angularResults_p << "#Number r_in = abs(pIn); phi_in = atan2(pIn[1],pIn[0])*180/pi; theta_in = atan2(r_in,pIn[2])*180/pi; r_out = abs(pOut); phi_out = atan2(pOut[1],pOut[0])*180/pi; theta_out = atan2(r_out,pOut[2])*180/pi" << std::endl;
          angularResults_p << "# > pIn[0],pIn[1],pIn[2] = x,y and z components of hyst operator computed from xIn" << std::endl;
          angularResults_p << "# > pOut[0],pOut[1],pOut[2] = x,y and z components of hyst operator computed from xOut" << std::endl;

          angularResults_y << "#Number r_in = abs(yIn); phi_in = atan2(yIn[1],yIn[0])*180/pi; theta_in = atan2(r_in,yIn[2])*180/pi; r_out = abs(yOut); phi_out = atan2(yOut[1],yOut[0])*180/pi; theta_out = atan2(r_out,yOut[2])*180/pi" << std::endl;
          angularResults_y << "# > yIn[0],yIn[1],yIn[2] = x,y and z components of output of hyst operator to xIn; used as input for inversion" << std::endl;
          angularResults_y << "# > yOut[0],yOut[1],yOut[2] = x,y and z components computed from xOut" << std::endl;
        }
			} else {
        if(dimHystOperator == 2){
          results_x << "#Number xIn[0] xIn[1]" << std::endl;
          results_x << "# > xIn[0],xIn[1] = x and y components of input; given by test signal" << std::endl;

          results_p << "#Number pOut[0] pOut[1]" << std::endl;
          results_p << "# > pOut[0],pOut[1] = x and y components of hyst operator computed from xOut" << std::endl;

          results_xp << "#Number\t xIn[0]\t xIn[1]\t pOut[0]\t pOut[1]" << std::endl;

          results_y << "#Number yOut[0] yOut[1]" << std::endl;
          results_y << "# > yOut[0],yOut[1] = x and y components computed from xOut" << std::endl;

          angularResults_x << "#Number abs(xIn) atan2(xIn[1]/xIn[0])*180/pi" << std::endl;
          angularResults_x << "# > xIn[0],xIn[1] = x and y components of input; given by test signal" << std::endl;

          angularResults_p << "#Number abs(pOut) atan2(pOut[1]/pOut[0])*180/pi" << std::endl;
          angularResults_p << "# > pOut[0],pOut[1] = x and y components of hyst operator computed from xOut" << std::endl;

          angularResults_y << "#Number abs(yOut) atan2(yOut[1]/yOut[0])*180/pi" << std::endl;
          angularResults_y << "# > yOut[0],yOut[1] = x and y components computed from xOut" << std::endl;
        } else {
          results_x << "#Number xIn[0] xIn[1] xIn[2]" << std::endl;
          results_x << "# > xIn[0],xIn[1],xIn[2] = x,y and z components of input; given by test signal" << std::endl;

          results_p << "#Number pOut[0] pOut[1] pOut[2]" << std::endl;
          results_p << "# > pOut[0],pOut[1],pOut[2] = x,y and z components of hyst operator computed from xOut" << std::endl;

          results_xp << "#Number\t xIn[0]\t xIn[1]\t xIn[2]\t pOut[0]\t pOut[1]\t pOut[2]" << std::endl;

          results_y << "#Number yOut[0] yOut[1] yOut[2]" << std::endl;
          results_y << "# > yOut[0],yOut[1],yOut[2] = x,y and z components computed from xOut" << std::endl;

          angularResults_x << "#Number r_in = abs(xIn); phi_in = atan2(xIn[1],xIn[0])*180/pi; theta_in = atan2(r_in,xIn[2])*180/pi" << std::endl;
          angularResults_x << "# > xIn[0],xIn[1],xIn[2] = x,y and z components of input; given by test signal" << std::endl;

          angularResults_p << "#Number r_out = abs(pOut); phi_out = atan2(pOut[1],pOut[0])*180/pi; theta_out = atan2(r_out,pOut[2])*180/pi" << std::endl;
          angularResults_p << "# > pOut[0],pOut[1],pOut[2] = x,y and z components of hyst operator computed from xOut" << std::endl;

          angularResults_y << "#Number r_out = abs(yOut); phi_out = atan2(yOut[1],yOut[0])*180/pi; theta_out = atan2(r_out,yOut[2])*180/pi" << std::endl;
          angularResults_y << "# > yOut[0],yOut[1],yOut[2] = x,y and z components computed from xOut" << std::endl;
        }
			}

      if(dimHystOperator == 2){
        orientationTowardsExcitation_p << "# Projetion of p_x and p_y to p_perpendicular and p_parallel" << std::endl;
        orientationTowardsExcitation_p << "# with p_parallel = p along excitation axis " << std::endl;
        orientationTowardsExcitation_p << "# and p_perpendicular = p perpendicualr to excitation axis " << std::endl;
        orientationTowardsExcitation_p << "# Columns: " << std::endl;
        orientationTowardsExcitation_p << "# 1. Testnumber " << std::endl;
        orientationTowardsExcitation_p << "# 2. Excitation amplitude: r_excite = abs(x_excite) " << std::endl;
        orientationTowardsExcitation_p << "# 3. Excitation angle towards x-axis: phi = atan2(x_ecite[1],x_ecite[0]) " << std::endl;
        orientationTowardsExcitation_p << "# 4. p parallel to excitation " << std::endl;
        orientationTowardsExcitation_p << "# 5. p perpendicular to excitation " << std::endl;
      } else {
        orientationTowardsExcitation_p << "# Projetion of p_x, p_y and p_z to p_perpendicular and p_parallel" << std::endl;
        orientationTowardsExcitation_p << "# with p_parallel = p along excitation axis " << std::endl;
        orientationTowardsExcitation_p << "# and p_perpendicular = p perpendicualr to excitation axis " << std::endl;
        orientationTowardsExcitation_p << "# Columns: " << std::endl;
        orientationTowardsExcitation_p << "# 1. Testnumber " << std::endl;
        orientationTowardsExcitation_p << "# 2. Excitation amplitude: r_excite = abs(x_excite) " << std::endl;
        orientationTowardsExcitation_p << "# 3. Excitation azimuth angle towards x-axis: phi = atan2(x_ecite[1],x_ecite[0]) " << std::endl;
        orientationTowardsExcitation_p << "# 4. Excitation inclination angle towards z-axis: theta = atan2(r_excite,x_ecite[2]) " << std::endl;
        orientationTowardsExcitation_p << "# 5. Component of p parallel to excitation " << std::endl;
        orientationTowardsExcitation_p << "# 6. Vector components of p perpendicular to excitation " << std::endl;
      }
		}


		bool overwriteMemory = false;
		// always overwrite Direction; deprecated flag; remove in future
    //		bool overwriteDirection = true;

		if( (printStatistics) ){
			std::cout << "##### STARTING TEST " << name << " #####" << std::endl;
		}
		hOutOld.Init();
    hOutOldForStrains.Init();

    Vector<Double> projectionDir = Vector<Double>(dimHystOperator);
    if(POL_operatorParams_.fixDirection_.NormL2() == 0){
//      std::cout << "No direction specified; using default x-direction for 1d test" << std::endl;
      // take x-direction as default case
      projectionDir.Init();
      projectionDir[0] = 1.0;
    } else {
      projectionDir = POL_operatorParams_.fixDirection_;
    }

    Matrix<Double> eps_mu_used = Matrix<Double>(dimHystOperator,dimHystOperator);
    if(dimHystOperator > dim_){
      eps_mu_used = eps_mu_baseFULL_;
    } else {
      eps_mu_used = eps_mu_base_;
    }
    Double eps_mu_scal;
    Vector<Double> tmp = Vector<Double>(dimHystOperator);

    eps_mu_used.Mult(projectionDir,tmp);
    tmp.Inner(projectionDir,eps_mu_scal);

//    std::cout << "eps_mu_scal: " << eps_mu_scal << std::endl;
//    std::cout << "eps_nu_base_[0][0]: " << eps_nu_base_[0][0] << std::endl;
//
    UInt outputEveryNthStep = 10; // 10
		for(UInt i = 0; i < totalSteps; i++){
			if( (i%outputEveryNthStep == 0)&&(printStatistics) ){
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
      xIn[1] = yVals[i];
      if(dimHystOperator == 3){
        xIn[2] = zVals[i];
      }

           // new usage of test1D
      // if test1D is set, project input onto POL_operatorParams_.fixDirection_
      // if this vector is empty, use x-axis as previously
			if(test1D == true){
//        std::cout << "Testing only along fixDirection = " << POL_operatorParams_.fixDirection_.ToString() << std::endl;
//        std::cout << "Demanded input = " << xIn.ToString() << std::endl;
        Double projection;
        projectionDir.Inner(xIn,projection);
        xIn.Init();
        xIn.Add(projection,projectionDir);
//        std::cout << "Used input = " << xIn.ToString() << std::endl;
			}

      xInBak = xIn;

      /*
       * Note 25.10.2020
       * Change:
       * - in case of inversion test, we have to measure the evaluation time during the
       *   initial evaluation, i.e., with the prescribed input
       * Reason:
       * - if inversion fails, the retrieved input might be completely different to the actual
       *   desired one; this can lead to completely different runtimes which are required for the
       *   models evaluation
       * - Example: during performance tests for the Mayergoyz model, the evaluation time for fixed-point
       *   inversion was much smaller than the evaluation time for Newton based inversion; this was caused
       *   by the fixed-point inversion which could not invert the hystoperator correctly; thus, during the
       *   evaluation with time-measurement, the Mayergoyz model became a much weaker input signal than it
       *   should have, which does not use the single scalar models to a full extend; thus it was cheaper
       */
      
      
      if(testInversion){
				/*
         * Inversion test:
         *	1. evaluate forward hyst operator with given input; overwrite memory = false
         *      > measure forward time
         *  2. retrieve input by passing result of 2 to inverse hyst operator
         *			> measure backward
         *  3. pass retrieved input to forward hyst operator; overwrite memroy = true
         *			CHANGED > do not measure forward time
         */

				// 1. forward; do not overwrite
				overwriteMemory = false;

				successFlagForward = -1;
        hIn.Init();
        
        if (measurePerformance) {
          forwardTimer->Start();
          startTime = forwardTimer->GetCPUTime();
        }

				hIn = hystTMP->computeValue_vec(xIn, 0, overwriteMemory, debugOut, successFlagForward);

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

				yIn.Init();
        yIn.Add(1.0,hIn);
        
        for(UInt j = 0; j < dimHystOperator; j++){
          yIn[j] += eps_mu_used[j][j]*xIn[j];
        }
        
//		std::cerr << "Input xIn = " << xIn.ToString(TS_PLAIN,separatorString,precisionDigits) << std::endl;
//        std::cerr << "Computed hIn = " << hIn.ToString(TS_PLAIN,separatorString,precisionDigits) << std::endl;
//        std::cerr << "Computed yIn = " << yIn.ToString(TS_PLAIN,separatorString,precisionDigits) << std::endl;

				// 2. backward; do not overwrite
				overwriteMemory = false;
        successFlagBackward = -1;
				numberOfLMIterations = 0;
				numberOfLinesearchIterations = 0;
				maxNumberOfLinesearchIterations = 0;

        Double scalIn = 0.0;
        projectionDir.Inner(yIn,scalIn);
//        Double scalIn = yIn[0];
        Double scalOut = 0;

        if (measurePerformance) {
          backwardTimer->Start();
					startTime = backwardTimer->GetCPUTime();
        }

        if((InversionParams_.inversionMethod != LOCAL_EVERETTBASED)&&(InversionParams_.inversionMethod != LOCAL_NOTIMPLEMENTED)){
          xOut = hystTMP->computeInput_vec_withStatistics(yIn, yPrev, xPrev, hPrev,
                  0, eps_mu_used, POL_operatorParams_.fieldsAlignedAboveSat_, POL_operatorParams_.hystOutputRestrictedToSat_,
                  numberOfLMIterations, numberOfLinesearchIterations, maxNumberOfLinesearchIterations,
                  successFlagBackward, minAlpha, maxAlpha, avgAlpha, xIn);
        } else if(InversionParams_.inversionMethod == LOCAL_NOTIMPLEMENTED) {
          EXCEPTION("Everett based inversion of the Mayergoyz model not supported yet as its output cannot be used in the forward model");
          // Problem: Throughout coefFunctionHyst, we need both, the forward and the backward/inverse Hysteresis model
          // the backward model is required, if B is known, but sometimes H is known or set, e.g., in testing of the hyst
          // operator or for the computation of the Jacobian where H is set and P needs to be retrieved; in these cases
          // we need the inverse of the backward model which unfortunately is NOT the forward Mayergoyz model!
          // the issue lies in the fact, that the sum over the inverted scalar models is not the same as the inverse of the
          // summed up forward scalar models; 
//          bool useEverett = true;
//          bool overwriteMemory = true;
//
//          xOut = hyst_->computeInput_vec(yIn,0,eps_mu_used,
//                    POL_operatorParams_.fieldsAlignedAboveSat_,POL_operatorParams_.hystOutputRestrictedToSat_,
//                    successFlagBackward, useEverett, overwriteMemory);
        } else {
//					xOut[0] = hystTMP->computeInputAndUpdate(yIn[0], eps_mu_base_[0][0],
//                  0, overwriteMemory, successFlagBackward);
          scalOut = hystTMP->computeInputAndUpdate(scalIn, eps_mu_scal,
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
            backwardReuseCounter++;
					}
        }

        if(InversionParams_.inversionMethod == 4){
          xOut.Init();
          xOut.Add(scalOut,projectionDir);
//          xOut[0] = scalOut;
        }

//        std::cout << "Computed xOut = " << xOut.ToString() << std::endl;
//
				if(successFlagBackward == -1){
          totalFails++;
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
        } else if((successFlagBackward == 5) || (successFlagBackward == 7)){
          totalPassedResTolX++;
          if(successFlagBackward == 7){
            totalFixedPoint++;
          } else {
            totalNewtonBased++;
          }
        } else if((successFlagBackward == 6) || (successFlagBackward == 8)){
          totalPassedResTolY++;
          if(successFlagBackward == 8){
            totalFixedPoint++;
          } else {
            totalNewtonBased++;
          }
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
      // Note 25.10.2020: If inversion is tested, measure forward time during the initial
      //  evaluation (step 1.) due to the same reason as in the NOTE below (failed inverison leads
      //  to wrong input which might be much faster or slower to be evaluated)
      // NOTE: in case of inversion, we set the memory of the system with the RETRIEVED SOLUTION
      // i.e. if inversion fails, it will also affect the reference solution!
			overwriteMemory = true;

			successFlagForward = -1;
			hOut.Init();
      hOutForStrains.Init();

      Double hinScal = 0.0;
      projectionDir.Inner(xIn,hinScal);
//      Double houtScal = 0.0;

      if ((POL_operatorParams_.printOut_ > 0)) {
        static UInt cntTS = 1;
        hystTMP->collectParallelProjections(true, cntTS);
        cntTS++;
      }
      
			if (measurePerformance && (testInversion == false)) {
				forwardTimer->Start();
				startTime = forwardTimer->GetCPUTime();
			}

//			if(vector){
				hOut = hystTMP->computeValue_vec(xIn, 0, overwriteMemory, debugOut, successFlagForward);
//			} else {
//				houtScal = hystTMP->computeValueAndUpdate(hinScal, 0, overwriteMemory, successFlagForward);
//			}

			if (measurePerformance && (testInversion == false)) {
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

      if ((POL_operatorParams_.printOut_ > 0)) {
        static UInt cnt = 1;
        bool appendToTxt = false;
        if(cnt > 1){
          appendToTxt = true;
        } 
        
        if (((cnt-1) % POL_operatorParams_.printOut_ == 0)) {
          std::stringstream filenamebuf;
          filenamebuf << "PreisachPlane_Step_" << std::setfill('0') << std::setw(5) << cnt << ".bmp";
//            filenamebuf << "PreisachPlane_Step_" << std::setfill('0') << std::setw(5) << cnt << "_v" << POL_operatorParams_.evalVersion_ << "_numRows" << POL_weightParams_.numRows_ << ".bmp";
          std::stringstream filenamebuf2;
          filenamebuf2 << "NestedList.txt";
          std::stringstream optionalHeaderStream;
          if(cnt == 1){
            optionalHeaderStream << "Parameter of vector Preisach model based on rotational operators: \n";
            if(POL_operatorParams_.isClassical_){
              optionalHeaderStream << "+++ version: classical \n";
              optionalHeaderStream << "+++ rotational resistance: " << POL_operatorParams_.rotResistance_ << "\n";
            } else {
              optionalHeaderStream << "+++ version: revised \n";
              optionalHeaderStream << "+++ rotational resistance: " << POL_operatorParams_.rotResistance_ << "\n";
              optionalHeaderStream << "+++ angular distance: " << POL_operatorParams_.angularDistance_ << "\n"; 
            }
            optionalHeaderStream << "\n";
          }
          optionalHeaderStream << "Step: " << std::setfill('0') << std::setw(5) << cnt << "\n";
          optionalHeaderStream << "+++ Normalized input vector: (" << xIn[0]/POL_operatorParams_.inputSat_ << ", " << xIn[1]/POL_operatorParams_.inputSat_ << ") \n";

          if(cnt == 1){
            std::cout << "--- Printing state of Preisach plane to " << filenamebuf.str() << std::endl;
            std::cout << "--- Printing nested list to " << filenamebuf2.str() << std::endl;            
          }

          hystTMP->switchingStateToBmp(POL_operatorParams_.bmpResolution_, filenamebuf.str(), 0, true);  
          hystTMP->rotationListToTxt(filenamebuf2.str(), 0, appendToTxt, optionalHeaderStream.str());
          // deactivate collection to save runtime
          hystTMP->collectParallelProjections(false, 0);
        }
        cnt++;
      }
        
//      if(!vector){
//        hOut.Init();
//        hOut.Add(houtScal,projectionDir);
//      }

      if(CouplingParams_.ownHystOperator_){
        if(vector){
          hOutForStrains = hystStrainTMP->computeValue_vec(xIn, 0, overwriteMemory, debugOut, successFlagForward);
        } else {
          hOutForStrains.Init();
          Double scalOut = hystStrainTMP->computeValueAndUpdate(hinScal, 0, overwriteMemory, successFlagForward);
          hOutForStrains.Add(scalOut,projectionDir);
        }
      } else {
        hOutForStrains = hOut;
      }

			if(successFlagForward == -1){
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
			for(UInt j = 0; j < dimHystOperator; j++){
        // note that in case of inversion, xIn is the RETRIEVED input!
				yOut[j] += eps_mu_used[j][j]*xIn[j];
			}
      
      if(outputIrrStrains){
        Vector<Double> S_irr = ComputeIrreversibleStrains(hOutForStrains,xIn,hOutOldForStrains);
        results_s << i+1 << separatorString << S_irr.ToString(TS_PLAIN,separatorString,precisionDigits) << std::endl;
        results_xps << i+1 << separatorString << xIn.ToString(TS_PLAIN,separatorString,precisionDigits) << separatorString;
        results_xps << hOutForStrains.ToString(TS_PLAIN,separatorString,precisionDigits) << separatorString;
        results_xps << S_irr.ToString(TS_PLAIN,separatorString,precisionDigits) << std::endl;
      }
      if(writeResultsToFile){
        
        if(writeResultsToInfoXML_){
          Double Px,Py,Pz,Pabs,Pangle1,Pangle2;
          Px = hOut[0];
          Py = hOut[1];
          Pabs = hOut.NormL2();
          Pangle1 = std::atan2(hOut[1],hOut[0])/M_PI*180;
          if(dimHystOperator == 3){
            Pz = hOut[2];
            Pangle2 = std::atan2(Pabs,hOut[2])/M_PI*180;
          } else {
            Pz = 0.0;
            Pangle2 = 0.0;
          }
          WriteTestStepToInfoXML(i+1,Px,Py,Pz,Pabs,Pangle1,Pangle2);
        }
        
        //results_xp << i+1 << " " << std::setprecision(precisionDigits) << xIn.ToString() << " " << hOut.ToString() << std::endl;
        //        results_xp << i+1 << separatorString << std::setprecision(precisionDigits) << xIn[0] << " " << xIn[1] << " " << hOut[0] << " " << hOut[1] << std::endl;
        results_xp << i+1 << separatorString << xIn.ToString(TS_PLAIN,separatorString,precisionDigits) << separatorString;
        results_xp << hOut.ToString(TS_PLAIN,separatorString,precisionDigits) << std::endl;
        
        Double xAmpl = xIn.NormL2();
        Double phi = std::atan2(xIn[1],xIn[0])/M_PI*180;
        Double theta = 0;
        if(dimHystOperator == 3){
          theta = std::atan2(xAmpl,xIn[2])/M_PI*180;
        }
        Vector<Double> xDir = Vector<Double>(dimHystOperator);
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
        Vector<Double> pPerpendicularVec = Vector<Double>(dimHystOperator);
        pPerpendicularVec.Add(1.0,hOut,-pParallel,xDir);
        Double pPerpendicular = pPerpendicularVec.NormL2();
        if(dimHystOperator == 2){
          orientationTowardsExcitation_p << i+1 << separatorString << std::setprecision(precisionDigits) << xAmpl << separatorString << phi;
          orientationTowardsExcitation_p << std::setprecision(precisionDigits) << separatorString << pParallel << separatorString << pPerpendicular << std::endl;
        } else {
          orientationTowardsExcitation_p << i+1 << separatorString << std::setprecision(precisionDigits) << xAmpl << separatorString << phi << separatorString << theta;
          orientationTowardsExcitation_p << std::setprecision(precisionDigits) << separatorString << pParallel << separatorString << pPerpendicularVec.ToString(TS_PLAIN,separatorString,precisionDigits) << std::endl;
        }
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

				if(xErr.NormL2() > InversionParams_.tolH){
					failWRTX = true;
					failedTests_xIn[i] = xInBak;
					failedTests_xOut[i] = xOut;
					failedTests_yIn[i] = yIn;
					failedTests_yOut[i] = yOut;
				}

				if(yErr.NormL2() > InversionParams_.tolB){
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
			if(writeResultsToFile){
        Double r_vec1 = 0;
        Double phi_vec1 = 0;
        Double theta_vec1 = 0;
        
        Double r_vec2 = 0;
        Double phi_vec2 = 0;
        Double theta_vec2 = 0;

				if(testInversion){				
					// input to hyst operator and retrieved input from inversion (xIn, xOut)
          results_x << i+1 << separatorString << xInBak.ToString(TS_PLAIN,separatorString,precisionDigits);
          results_x << separatorString << xOut.ToString(TS_PLAIN,separatorString,precisionDigits) << std::endl;
                  
          r_vec1 = xInBak.NormL2();
          phi_vec1 = atan2(xInBak[1],xInBak[0])*180.0/M_PI;
          
          angularResults_x << i+1 << std::setprecision(precisionDigits);
          angularResults_x << separatorString << r_vec1;
          angularResults_x << separatorString << phi_vec1;
          if(dimHystOperator == 3){
            theta_vec1 = atan2(r_vec1,xInBak[2])*180.0/M_PI;
            angularResults_x << separatorString << theta_vec1;
          }
          
          r_vec2 = xOut.NormL2();
          phi_vec2 = atan2(xOut[1],xOut[0])*180.0/M_PI;
          
          angularResults_x << separatorString << r_vec2;
          angularResults_x << separatorString << phi_vec2;
          if(dimHystOperator == 3){
            theta_vec2 = atan2(r_vec2,xOut[2])*180.0/M_PI;
            angularResults_x << separatorString << theta_vec2;
          }
          angularResults_x << std::endl;
          
          // polarization given and computed (hIn, hOut)
          results_p << i+1 << separatorString << hIn.ToString(TS_PLAIN,separatorString,precisionDigits);
          results_p << separatorString << hOut.ToString(TS_PLAIN,separatorString,precisionDigits) << std::endl;
                  
          r_vec1 = hIn.NormL2();
          phi_vec1 = atan2(hIn[1],hIn[0])*180.0/M_PI;
          
          angularResults_p << i+1 << std::setprecision(precisionDigits);
          angularResults_p << separatorString << r_vec1;
          angularResults_p << separatorString << phi_vec1;
          if(dimHystOperator == 3){
            theta_vec1 = atan2(r_vec1,hIn[2])*180.0/M_PI;
            angularResults_p << separatorString << theta_vec1;
          }
          
          r_vec2 = hOut.NormL2();
          phi_vec2 = atan2(hOut[1],hOut[0])*180.0/M_PI;
          
          angularResults_p << separatorString << r_vec2;
          angularResults_p << separatorString << phi_vec2;
          if(dimHystOperator == 3){
            theta_vec2 = atan2(r_vec2,hOut[2])*180.0/M_PI;
            angularResults_p << separatorString << theta_vec2;
          }
          angularResults_p << std::endl;
          
          // output given and computed (yIn, yOut)
          results_y << i+1 << separatorString << yIn.ToString(TS_PLAIN,separatorString,precisionDigits);
          results_y << separatorString << yOut.ToString(TS_PLAIN,separatorString,precisionDigits) << std::endl;
                  
          r_vec1 = yIn.NormL2();
          phi_vec1 = atan2(yIn[1],yIn[0])*180.0/M_PI;
          
          angularResults_y << i+1 << std::setprecision(precisionDigits);
          angularResults_y << separatorString << r_vec1;
          angularResults_y << separatorString << phi_vec1;
          if(dimHystOperator == 3){
            theta_vec1 = atan2(r_vec1,yIn[2])*180.0/M_PI;
            angularResults_y << separatorString << theta_vec1;
          }
          
          r_vec2 = yOut.NormL2();
          phi_vec2 = atan2(yOut[1],yOut[0])*180.0/M_PI;
          
          angularResults_y << separatorString << r_vec2;
          angularResults_y << separatorString << phi_vec2;
          if(dimHystOperator == 3){
            theta_vec2 = atan2(r_vec2,yOut[2])*180.0/M_PI;
            angularResults_y << separatorString << theta_vec2;
          }
          angularResults_y << std::endl;
          
				} else {
					// input to hyst operator (xIn)
          results_x << i+1 << separatorString << xInBak.ToString(TS_PLAIN,separatorString,precisionDigits);
          results_x << std::endl;
                  
          r_vec1 = xInBak.NormL2();
          phi_vec1 = atan2(xInBak[1],xInBak[0])*180.0/M_PI;
          
          angularResults_x << i+1 << std::setprecision(precisionDigits);
          angularResults_x << separatorString << r_vec1;
          angularResults_x << separatorString << phi_vec1;
          if(dimHystOperator == 3){
            theta_vec1 = atan2(r_vec1,xInBak[2])*180.0/M_PI;
            angularResults_x << separatorString << theta_vec1;
          }
          angularResults_x << std::endl;
          
          // polarization computed (hOut)
          results_p << i+1;
          results_p << separatorString << hOut.ToString(TS_PLAIN,separatorString,precisionDigits) << std::endl;
          
          angularResults_p << i+1 << std::setprecision(precisionDigits);
          r_vec2 = hOut.NormL2();
          phi_vec2 = atan2(hOut[1],hOut[0])*180.0/M_PI;
          
          angularResults_p << separatorString << r_vec2;
          angularResults_p << separatorString << phi_vec2;
          if(dimHystOperator == 3){
            theta_vec2 = atan2(r_vec2,hOut[2])*180.0/M_PI;
            angularResults_p << separatorString << theta_vec2;
          }
          angularResults_p << std::endl;
          
          // output computed (yOut)
          results_y << i+1;
          results_y << separatorString << yOut.ToString(TS_PLAIN,separatorString,precisionDigits) << std::endl;
                  
          angularResults_y << i+1 << std::setprecision(precisionDigits);        
          r_vec2 = yOut.NormL2();
          phi_vec2 = atan2(yOut[1],yOut[0])*180.0/M_PI;
          
          angularResults_y << separatorString << r_vec2;
          angularResults_y << separatorString << phi_vec2;
          if(dimHystOperator == 3){
            theta_vec2 = atan2(r_vec2,yOut[2])*180.0/M_PI;
            angularResults_y << separatorString << theta_vec2;
          }
          angularResults_y << std::endl;
				}
			}

			LOG_TRACE(coeffunctionhyst_main) << "##### TARGET X-VECTOR #####";
			LOG_TRACE(coeffunctionhyst_main) << xInBak.ToString(TS_PLAIN,separatorString,precisionDigits);
			if(testInversion){
				LOG_TRACE(coeffunctionhyst_main) << "##### RETRIEVED X-VECTOR #####";
				LOG_TRACE(coeffunctionhyst_main) << xOut.ToString(TS_PLAIN,separatorString,precisionDigits);

				LOG_TRACE(coeffunctionhyst_main) << "##### ERROR VECTOR wrt X #####";
				LOG_TRACE(coeffunctionhyst_main) << xErr.ToString(TS_PLAIN,separatorString,precisionDigits);
				LOG_TRACE(coeffunctionhyst_main) << "##### TARGET Y-VECTOR #####";
				LOG_TRACE(coeffunctionhyst_main) << yIn.ToString(TS_PLAIN,separatorString,precisionDigits);
			}
			LOG_TRACE(coeffunctionhyst_main) << "##### RETRIEVED Y-VECTOR #####";
			LOG_TRACE(coeffunctionhyst_main) << yOut.ToString(TS_PLAIN,separatorString,precisionDigits);
			if(testInversion){
				LOG_TRACE(coeffunctionhyst_main) << "##### ERROR VECTOR wrt Y #####";
				LOG_TRACE(coeffunctionhyst_main) << yErr.ToString(TS_PLAIN,separatorString,precisionDigits);
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
			majorResults << " " << totalSteps-numFails << " of " << totalSteps << " satisfied at least the failback criterion, i.e. |residual_Y| = |Y - mu_eps*X - P(X)| < " << InversionParams_.tolB << std::endl;
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
					statistics << " Test Nr " << mapItX->first + 1 << std::endl;
					failX = mapItX->second.first;
					failY = mapItY->second.first;
//          std::cout << "2" << std::endl;
					bool xsolabovesat = failedTests_xIn[mapItX->first].NormL2()>POL_operatorParams_.inputSat_;
					bool xretabovesat = failedTests_xOut[mapItX->first].NormL2()>POL_operatorParams_.inputSat_;
					statistics << " > X_in: " << failedTests_xIn[mapItX->first].ToString(TS_PLAIN,separatorString,precisionDigits)  << "; above sat? "<< xsolabovesat << std::endl;
					statistics << " > X_out: " << failedTests_xOut[mapItX->first].ToString(TS_PLAIN,separatorString,precisionDigits) << "; above sat? "<< xretabovesat << std::endl;
					if(failX){
						statistics << " > FAIL - remaining error wrt X  " << mapItX->second.second << std::endl;
					} else {
						statistics << " > SUCCESS - remaining error wrt X  " << mapItX->second.second << std::endl;
					}
//					std::cout << "3" << std::endl;
					statistics << " > Y_in: " << failedTests_yIn[mapItX->first].ToString(TS_PLAIN,separatorString,precisionDigits) << std::endl;
					statistics << " > Y_out: " << failedTests_yOut[mapItX->first].ToString(TS_PLAIN,separatorString,precisionDigits) << std::endl;
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
			statistics << " " << totalReused << " of " << totalSteps << " reused old solution as DeltaY < " << InversionParams_.tolB << std::endl;
			statistics << " " << totalAnhystOnly << " of " << totalSteps << " were solved using bisection for pure anhysteretic case" << std::endl;
			statistics << " " << totalBisection << " of " << totalSteps << " were solved using bisection / simple division" << std::endl;
			if(!vector){
				statistics << " " << totalPassedErrorTol << " of " << totalSteps << " tests passed due to error tolerance" << std::endl;
			}

			if(vector){
				statistics << " " << totalRemanence << " of " << totalSteps << " cases were in remanence" << std::endl;
				statistics << " " << totalFixedPoint << " of " << totalSteps << " cases were solved via local FP iteration" << std::endl;
        statistics << " " << totalNewtonBased << " of " << totalSteps << " cases were solved by a Newton based approach " << std::endl;
        //statistics << " " << LMcases << " of " << totalSteps << " were solved using Levenberg-Marquardt, Newton or Fixpoint" << std::endl;
				if(LMcases != 0){
					statistics << "## Detailed Statistics on local inversion: " << std::endl;
					statistics << " " << totalPassedErrorTol << " of " << LMcases << " tests passed due to error tolerance |JacT*Res| < " << InversionParams_.tolH << "(not checked for Fixpoint)" << std::endl;
					statistics << " " << totalPassedResTolX << " of " << LMcases << " tests passed due to |residual_X| = |X - mu_eps^-1*(Y - P(X)| < " << InversionParams_.tolH << std::endl;
					statistics << " " << totalPassedResTolY << " of " << LMcases << " tests passed due to |residual_Y| = |Y - mu_eps*X - P(X)| < " << InversionParams_.tolH << std::endl;
					statistics << " " << totalFails << " of " << LMcases << " failed the failback criterion (|residual_Y| < " << InversionParams_.tolB << ")" << std::endl;
					statistics << " Total number of outer iterations (LM,Newton,FP): " << totalNumberOfLMIterations << std::endl;
					statistics << " Average number of outer iterations (LM,Newton,FP): " << (Double) totalNumberOfLMIterations/LMcases << std::endl;
					statistics << " Total number of inner iterations (Linesearch inside LM,Newton): " << totalNumberOfLinesearchIterations << std::endl;
					statistics << " Average number of inner iterations (Linesearch inside LM,Newton): " << (Double) totalNumberOfLinesearchIterations/totalNumberOfLMIterations << std::endl;
					statistics << " Maximal number of inner iterations (Linesearch inside LM,Newton): " << absmaxNumberOfLinesearchIterations << std::endl;
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
        performance << "Used method for inversion: " << usedInversionMethod << std::endl;
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
        bool detailled = false;
        if(addHeaderToCommonPerFile){
          commonPerfStream << std::fixed << std::setw(minWidthOutputColumns) << std::setprecision(6) << "NAMETAG";
          commonPerfStream << "\t\t Average time for forward eval \t\t Number of forward evals \t---\t ";
          commonPerfStream << "INVERSIONMETHOD";
          commonPerfStream << "\t\t Average time for inversion \t\t Number of evals for inversion (excluding reused values) \t\t Number of failed inversions";
          commonPerfStream << std::endl;
        }
        std::string INVERSIONMETHOD = usedInversionMethod;
        Double backwardAvgEvalTime;
        if(testInversion){
          backwardAvgEvalTime = backwardTotalEvalTime/backwardEvalCounter;
        } else {
          backwardAvgEvalTime = 0;
          backwardTotalEvalTime = 0;
          backwardEvalCounter = 0;
        }
        if(detailled){
          commonPerfStream << std::fixed << std::setw(minWidthOutputColumns) << std::setprecision(6) << "Tag: " << nameTagForPerfFile;
          commonPerfStream << "\t\t Avg forward time: " << forwardAvgEvalTime*1000 << " ms (Evals=" << forwardEvalCounter << ") ";
          commonPerfStream << "\t---\t Inverison method: " << usedInversionMethod;
          commonPerfStream << "\t\t Avg backward time: " << backwardAvgEvalTime*1000 << " ms (Evals=" << backwardEvalCounter << ", Reused=" << backwardReuseCounter << ") ";
          commonPerfStream << "\t InversionFails=" << numFails << std::endl;
        } else {
          commonPerfStream << std::fixed << std::setw(minWidthOutputColumns) << std::setprecision(6) << nameTagForPerfFile;
          commonPerfStream << "\t\t" << forwardAvgEvalTime*1000 << " ms \t" << forwardEvalCounter;
          commonPerfStream << "\t---\t" << usedInversionMethod;
          commonPerfStream << "\t\t" << backwardAvgEvalTime*1000 << " ms \t" << backwardEvalCounter << "( + "<< backwardReuseCounter << " reused)";
          commonPerfStream << "\t" << numFails << std::endl;
//          commonPerfStream << std::fixed << std::setprecision(6) << nameTagForPerfFile << "\t" << forwardAvgEvalTime << "\t" << forwardEvalCounter << "\t" << backwardAvgEvalTime << "\t" << backwardEvalCounter << "\t" << numFails << std::endl;
        }
//        if(testInversion){
//          if(detailled){
//            commonPerfStream << std::fixed << std::setprecision(6) << nameTagForPerfFile << "\t\t Forward: " << forwardAvgEvalTime << " (Evals=" << forwardEvalCounter << ") " << "\t Backward: " << backwardAvgEvalTime << " (Evals=" << backwardEvalCounter << ") " << "\t InversionFails=" << numFails << std::endl;
//          } else {
//            commonPerfStream << std::fixed << std::setprecision(6) << nameTagForPerfFile << "\t" << forwardAvgEvalTime << "\t" << forwardEvalCounter << "\t" << backwardAvgEvalTime << "\t" << backwardEvalCounter << "\t" << numFails << std::endl;
//          }
//        } else {
//          if(detailled){
//            commonPerfStream << std::fixed << std::setprecision(6) << nameTagForPerfFile << "\t\t Forward: " << forwardAvgEvalTime << " (Evals=" << forwardEvalCounter << ") " << "\t Backward: NaN (Evals=NaN) \t InversionFails=NaN" << std::endl;
//          } else {
//            commonPerfStream << std::fixed << std::setprecision(6) << nameTagForPerfFile << "\t" << forwardAvgEvalTime << "\t" << forwardEvalCounter << "\t" << "NaN" << "\t" << "NaN" << "\t" << "NaN" << std::endl;
//          }
//        }
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

    /*
     * clean up
     */
    if(dimHystOperator > dim_){
      POL_operatorParams_.fixDirection_.Resize(dim_);
      POL_operatorParams_.startingAxisMG_.Resize(dim_);
    }
    if(hystTMP != NULL){
      delete hystTMP;
    }
    if(hystStrainTMP != NULL){
      delete hystStrainTMP;
    }
    if(forwardTimer != NULL){
      delete forwardTimer;
    }
    if(backwardTimer != NULL){
      delete backwardTimer;
    }
  }

  // Previous version
//    void CoefFunctionHyst::TestHystOperatorWithSignal(std::string name, Vector<Double> xVals, Vector<Double> yVals,
//          Vector<Double> zVals, UInt dimHystOperator, 
//          bool testInversion, bool printStatistics, bool writeResultsToFile,
//          bool measurePerformance, std::string commonPerformanceFile, bool test1D, bool outputIrrStrains,
//          std::string nameTagForPerfFile){
//
//		/*
//     * 0. Declare variables (there are alot)
//     */
//		std::ofstream statistics;
//		std::ofstream results_x;
//    std::ofstream results_p;
//    std::ofstream results_s;
//    std::ofstream results_xp;
//    std::ofstream results_xps;
//    std::ofstream results_y;
//    std::ofstream angularResults_x;
//    std::ofstream angularResults_p;
//    std::ofstream angularResults_y;
//		std::ofstream performance;
//    std::ofstream orientationTowardsExcitation_p;
//    std::ofstream commonPerfStream;
//
//    std::string basedir = "./history_hystOperator/";
//
//    try {
//      fs::create_directory( basedir );
//    } catch (std::exception &ex) {
//      EXCEPTION(ex.what());
//    }
//
//		std::stringstream statistics_name;
//		statistics_name << basedir << material_->GetName() << name << "_statistics";
//
//		std::stringstream results_name_x;
//		results_name_x << basedir << material_->GetName() << name << "_results_x";
//
//    std::stringstream results_name_p;
//		results_name_p << basedir << material_->GetName() << name << "_results_p";
//
//    std::stringstream results_name_s;
//		results_name_s << basedir << material_->GetName() << name << "_results_s";
//
//    std::stringstream results_name_y;
//		results_name_y << basedir << material_->GetName() << name << "_results_y";
//
//    std::stringstream results_name_xp;
//		results_name_xp << basedir << material_->GetName() << name << "_results_xp";
//
//    std::stringstream results_name_xps;
//		results_name_xps << basedir << material_->GetName() << name << "_results_xps";
//
//    std::stringstream angularResults_name_x;
//		angularResults_name_x << basedir << material_->GetName() << name << "_angularResults_x";
//
//    std::stringstream angularResults_name_p;
//		angularResults_name_p << basedir << material_->GetName() << name << "_angularResults_p";
//
//    std::stringstream angularResults_name_y;
//		angularResults_name_y << basedir << material_->GetName() << name << "_angularResults_y";
//
//		std::stringstream performance_name;
//		performance_name << basedir << material_->GetName() << name << "_performance";
//
//    std::stringstream orientation_name;
//		orientation_name << basedir << material_->GetName() << name << "_projectedCoords_p";
//
//    std::stringstream commonPerfStream_name;
//		commonPerfStream_name << basedir << material_->GetName() << commonPerformanceFile;
//
//		if(writeResultsToFile){
//			results_x.open(results_name_x.str());
//      results_xp.open(results_name_xp.str());
//      results_p.open(results_name_p.str());
//      results_y.open(results_name_y.str());
//      angularResults_x.open(angularResults_name_x.str());
//      angularResults_p.open(angularResults_name_p.str());
//      angularResults_y.open(angularResults_name_y.str());
//      orientationTowardsExcitation_p.open(orientation_name.str());
//		}
//
//		if(printStatistics){
//			statistics.open(statistics_name.str());
//		}
//		if(measurePerformance){
//			performance.open(performance_name.str());
//      if(commonPerformanceFile != "---"){
//        commonPerfStream.open(commonPerfStream_name.str(),std::ios_base::app);
//      }
//		}
//
//		UInt totalSteps = xVals.GetSize();
//		UInt numFails = 0;
//    UInt numHalfFails = 0; // fail error crit but pass residual crit
//		int successFlagForward = -1;
//		// forward:
//		// -1 = fail
//		//  0 = reuse value
//		//  1 = anhyst only
//		//  2 = eval on permanent storage
//		//  3 = eval on temporal storage
//
//		int successFlagBackward = -1;
//		// backward:
//		// -1 = fail
//		//  0 = reuse value
//		//  1 = anhyst only
//		//  2 = bisection
//		//  3-6 only for vector implementation using Levenberg Marquardt
//		//  3 = reamnence
//		//  4 = passed dut to error tolerance
//		//  5 = passed due to tolerance wrt x
//		//  6 = passed due to tolerance wrt y
//
//		Vector<Double> xIn = Vector<Double>(dim_);
//		Vector<Double> hIn = Vector<Double>(dim_);
//		Vector<Double> yIn = Vector<Double>(dim_);
//
//    Vector<Double> xInPrev = Vector<Double>(dim_);
//		Vector<Double> xPrev = Vector<Double>(dim_);
//		Vector<Double> hPrev = Vector<Double>(dim_);
//		Vector<Double> yPrev = Vector<Double>(dim_);
//
//		Vector<Double> xOut = Vector<Double>(dim_);
//		Vector<Double> hOut = Vector<Double>(dim_);
//    Vector<Double> hOutForStrains = Vector<Double>(dim_);
//		Vector<Double> yOut = Vector<Double>(dim_);
//		Vector<Double> hOutOld = Vector<Double>(dim_);
//    Vector<Double> hOutOldForStrains = Vector<Double>(dim_);
//
//    Vector<Double> xInBak = Vector<Double>(dim_);
//		Vector<Double> xErr = Vector<Double>(dim_);
//		Vector<Double> yErr = Vector<Double>(dim_);
//
//		std::map< UInt, Vector<Double> > failedTests_xIn;
//		std::map< UInt, Vector<Double> > failedTests_xOut;
//		std::map< UInt, Vector<Double> > failedTests_yIn;
//		std::map< UInt, Vector<Double> > failedTests_yOut;
//		std::map< UInt, std::pair< bool, Double > > failedTests_wrtX;
//		std::map< UInt, std::pair< bool, Double > > failedTests_wrtY;
//
//		Vector<int> successCodeVectorForward = Vector<int>(totalSteps);
//		Vector<int> successCodeVectorBackward = Vector<int>(totalSteps);
//
//		Double minAlpha = 0;
//		Double maxAlpha = 0;
//		Double avgAlpha = 0;
//		Double totalminAlpha = 1e10;
//		Double totalmaxAlpha = -1e10;
//		Double totalavgAlpha = 0;
//
//		UInt numberOfLMIterations = 0;
//		UInt numberOfLinesearchIterations = 0;
//		UInt maxNumberOfLinesearchIterations = 0;
//
//		// counter for forward evaluation
//		UInt forwardFails = 0;
//		UInt forwardReused = 0;
//		UInt forwardAnhystOnly = 0;
//		UInt forwardOverwrite = 0;
//		UInt forwardTMP = 0;
//
//		// counter for backward evaluation / inversion
//		UInt LMFails = 0;
//		UInt totalReused = 0;
//		UInt totalAnhystOnly = 0;
//		UInt totalBisection = 0;
//		UInt totalRemanence = 0;
//		UInt totalPassedErrorTol = 0;
//		UInt totalPassedResTolX = 0;
//		UInt totalPassedResTolY = 0;
//		UInt totalNumberOfLMIterations = 0;
//		UInt totalNumberOfLinesearchIterations = 0;
//		UInt absmaxNumberOfLinesearchIterations = 0;
//
//		Double startTime = 0.0;
//    Double endTime = 0.0;
//    Double evalTime = 0.0;
//		Double backwardMaxEvalTime = -1.0;
//		Double backwardTotalEvalTime = 0.0;
//		Double forwardMaxEvalTime = -1.0;
//		Double forwardTotalEvalTime = 0.0;
//
//		UInt forwardEvalCounter = 0;
//		UInt backwardEvalCounter = 0;
//		Timer* forwardTimer;
//		Timer* backwardTimer;
//
//    UInt precisionDigits = 12;//9
//    
//		if(printStatistics){
//			statistics << "+++ STATISTICS +++" << std::endl;
//			statistics << "TEST: " << name << std::endl;
//			statistics << "MODEL: " << POL_operatorParams_.methodName_ << std::endl;
//      statistics << "Error Criteria: " << std::endl;
//			statistics << "- Residual wrt input (=x): " << InversionParams_.tolH << std::endl;
//      statistics << "- Residual wrt output (=y): " << InversionParams_.tolB << std::endl;
//		}
//
//		/*
//     * 1. Create temporal hyst operator; this should be done for each test to ensure
//     *		that operator is in initial state
//     */
//		bool isVirgin = true;
//		bool vector = true;
//		bool debugOut = false;
//
//		if (POL_operatorParams_.methodName_ == "scalarPreisach") {
//
//			POL_useExtension_ = false;
//      //			int useExtensionInt;
//      //			material_->GetScalar(useExtensionInt, SCALPREISACH_USE_EXT);
//			if(POL_useExtension_){
//				EXCEPTION("Extension not implemented for tests; remove completely as not working");
////				POL_useExtension_ = true;
////
////				material_->GetScalar(POL_operatorParams_.rotResistance_, ROT_RESISTANCE, Global::REAL);
////				material_->GetScalar(POL_operatorParams_.angularDistance_, ANG_DISTANCE, Global::REAL);
////
////				hystTMP = new ExtendedPreisach(1, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_, POL_weightParams_.weightTensor_,
////                POL_operatorParams_.rotResistance_, POL_operatorParams_.angularDistance_, dim_, isVirgin, POL_weightParams_.anhysteretic_a_,
////                POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
////
////				// set initial direction
////				Vector<Double> initialInput = Vector<Double>(dim_);
////				eps_nu_base_.Mult(POL_operatorParams_.fixDirection_,initialInput);
////				initialInput.ScalarMult(POL_operatorParams_.inputSat_);
////				initialInput.Add(POL_operatorParams_.outputSat_,POL_operatorParams_.fixDirection_);
////				// initialInput = (POL_operatorParams_.outputSat_*Identity + POL_operatorParams_.inputSat_*eps_nu_base_)*POL_operatorParams_.fixDirection_
////				// > flux with value just at saturation
////				hystTMP->UpdateRotationState(initialInput,eps_nu_base_,0);
////
//			} else {
//				vector = false;
//        bool ignoreAnhystPart = false;
//				hystTMP = new Preisach(1, POL_operatorParams_, POL_weightParams_, isVirgin, ignoreAnhystPart);
////                (1, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_,
////                POL_weightParams_.weightTensor_, isVirgin,
////                POL_weightParams_.anhysteretic_a_, POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,
////                POL_weightParams_.anhystOnly_);
//        test1D = true;
//        if(InversionParams_.inversionMethod != 4){
////          std::cout << "Inversion of scalar model with Newton/LM selected; will restrict input to operator direction \n"
////                  "as otherwise inversion will not work" << std::endl;
//          hystTMP->SetParamsForInversion(InversionParams_);
//        }
//      }
//		} else if (POL_operatorParams_.methodName_ == "vectorPreisach_Sutor") {
//			if(printStatistics){
//				if (POL_operatorParams_.evalVersion_ == 1) {
//					statistics << "- classical model, list based implementation" << std::endl;
//				} else if (POL_operatorParams_.evalVersion_ == 2) {
//					statistics << "- revised model, list based implementation" << std::endl;
//				} else if (POL_operatorParams_.evalVersion_ == 10) {
//					statistics << "- classical model, matrix based implementation" << std::endl;
//				} else if (POL_operatorParams_.evalVersion_ == 20) {
//					statistics << "- revised model, matrix based implementation" << std::endl;
//				}
//
//				statistics << "with rotational resistance = " << POL_operatorParams_.rotResistance_ << std::endl;
//				if(POL_operatorParams_.isClassical_ == false){
//					statistics << "and angular distance = " << POL_operatorParams_.angularDistance_ << std::endl;
//				}
//			}
//
//			if (POL_operatorParams_.evalVersion_ == 1) {
//				POL_operatorParams_.isClassical_ = true; // original vector preisach model -> sutor2012
//
//				hystTMP = new VectorPreisachSutor_ListApproach(1, POL_operatorParams_, POL_weightParams_, dim_, isVirgin);
////                (1, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_,
////                POL_weightParams_.weightTensor_, POL_operatorParams_.rotResistance_, dim_, isVirgin,
////                POL_operatorParams_.isClassical_, POL_operatorParams_.scaleUpToSaturation_,
////                POL_operatorParams_.angularDistance_,POL_operatorParams_.angularResolution_,POL_weightParams_.anhysteretic_a_,
////                POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
//			} else if (POL_operatorParams_.evalVersion_ == 2) {
//				POL_operatorParams_.isClassical_ = false; // revised vector preisach model -> sutor2015
//
//				hystTMP = new VectorPreisachSutor_ListApproach(1, POL_operatorParams_, POL_weightParams_, dim_, isVirgin);
////                (1, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_,
////                POL_weightParams_.weightTensor_, POL_operatorParams_.rotResistance_, dim_, isVirgin,
////                POL_operatorParams_.isClassical_, POL_operatorParams_.scaleUpToSaturation_,
////                POL_operatorParams_.angularDistance_,POL_operatorParams_.angularResolution_,POL_weightParams_.anhysteretic_a_,
////                POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
//			} else if (POL_operatorParams_.evalVersion_ == 10) {
//				POL_operatorParams_.isClassical_ = true; // original vector preisach model -> sutor2015; matrix based implementation
//
//				hystTMP = new VectorPreisachSutor_MatrixApproach(1, POL_operatorParams_, POL_weightParams_, dim_, isVirgin);
////                (1, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_,
////                POL_weightParams_.weightTensor_, POL_operatorParams_.rotResistance_, dim_, isVirgin,
////                POL_operatorParams_.isClassical_, POL_operatorParams_.scaleUpToSaturation_,
////                POL_operatorParams_.angularDistance_,POL_operatorParams_.angularResolution_,POL_weightParams_.anhysteretic_a_,
////                POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
//			} else if (POL_operatorParams_.evalVersion_ == 20) {
//				POL_operatorParams_.isClassical_ = false; // revised vector preisach model -> sutor2015; matrix based implementation
//
//				hystTMP = new VectorPreisachSutor_MatrixApproach(1, POL_operatorParams_, POL_weightParams_, dim_, isVirgin);
////                (1, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_,
////                POL_weightParams_.weightTensor_, POL_operatorParams_.rotResistance_, dim_, isVirgin,
////                POL_operatorParams_.isClassical_, POL_operatorParams_.scaleUpToSaturation_,
////                POL_operatorParams_.angularDistance_,POL_operatorParams_.angularResolution_,POL_weightParams_.anhysteretic_a_,
////                POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,POL_weightParams_.anhystOnly_);
//			} else {
//				EXCEPTION("POL_operatorParams_.evalVersion_ has to be one of the following: \n "
//                "1: classical vector model (sutor2012) \n"
//                "2: revised vector model (sutor2015) [DEFAULT] \n"
//                "10: classical vector model (sutor2012) - Matrix implementation, only for reference \n"
//                "20: revised vector model (sutor2015) - Matrix implementation, only for reference \n")
//			}
//      hystTMP->SetParamsForInversion(InversionParams_);
//
////			hystTMP->SetParamsForInversion(InversionParams_.inversionMethod, InversionParams_.maxNumIts, InversionParams_.maxNumLSIts,
////                InversionParams_.tolH, InversionParams_.tolB,
////                InversionParams_.jacRes, InversionParams_.alphaLSStart,InversionParams_.alphaLSMin,InversionParams_.alphaLSMax,
////                InversionParams_.stopLineSearchAtLocalMin,
////                POL_operatorParams_.angularClipping_);
//
//		} else if (POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") {
//      // basically a scalar model in multiple directions
//      // isotropic case: all scalar models are equal (same weights etc)
//      // anisotropic case: each model different; choice of directions matters; weights are harder to obtain
//      int isIsotropic = 1;
//      material_->GetScalar(isIsotropic, PREISACH_MAYERGOYZ_ISOTROPIC);
//      if(isIsotropic == 0){
////        if( (dim_ != 2) || (isIsotropic == 0)){
//        EXCEPTION("Mayergoyz vector model currently only implemented for 2d isotropic materials");
//      }
//
//      statistics << "with num directions = " << POL_operatorParams_.numDirections_ << std::endl;
//      statistics << "and output clipping mode = " << POL_operatorParams_.outputClipping_ << std::endl;
//      statistics << "Starting axis: " << POL_operatorParams_.startingAxisMG_.ToString() << std::endl;
//
//      /*
//       * IMPORTANT REMARK:
//       *  > although the Mayergoyz model is based on the scalar models
//       *     we are not allowed to directly apply the Preisach parameter for
//       *     the scalar case (i.e. the weights, the anhyst parameter and so on)
//       *  > make sure that the passed parameter are already transformed correctly
//       *      > see constructor above
//       */
//      hystTMP = new VectorPreisachMayergoyz(1, POL_operatorParams_, POL_weightParams_, dim_, isVirgin);
////              (1, POL_operatorParams_.numDirections_, POL_operatorParams_.inputSat_, POL_operatorParams_.outputSat_,
////              POL_weightParams_.weightTensor_,dim_,isVirgin,
////              POL_weightParams_.anhysteretic_a_, POL_weightParams_.anhysteretic_b_, POL_weightParams_.anhysteretic_c_,
////              POL_weightParams_.anhystOnly_,POL_operatorParams_.outputClipping_);
//
//      hystTMP->SetParamsForInversion(InversionParams_);
//
////			hystTMP->SetParamsForInversion(InversionParams_.inversionMethod, InversionParams_.maxNumIts, InversionParams_.maxNumLSIts,
////                InversionParams_.tolH, InversionParams_.tolB,
////                InversionParams_.jacRes, InversionParams_.alphaLSStart,InversionParams_.alphaLSMin,InversionParams_.alphaLSMax,
////                InversionParams_.stopLineSearchAtLocalMin,
////                POL_operatorParams_.angularClipping_);
//		} else {
//			EXCEPTION("Invalid model selected for inversion test");
//		}
//
//    if(CouplingParams_.ownHystOperator_ && CouplingParams_.couplingDefined_inMatFile_){
//      if (STRAIN_operatorParams_.methodName_ == "scalarPreisach") {
//        bool ignoreAnhystPart = false;
//        hystStrainTMP = new Preisach(1, STRAIN_operatorParams_, STRAIN_weightParams_, isVirgin, ignoreAnhystPart);
////                (1, STRAIN_operatorParams_.inputSat_, STRAIN_operatorParams_.outputSat_,
////                STRAIN_weightParams_.weightTensor_, isVirgin,
////                STRAIN_weightParams_.anhysteretic_a_, STRAIN_weightParams_.anhysteretic_b_, STRAIN_weightParams_.anhysteretic_c_,
////                STRAIN_weightParams_.anhystOnly_);
//
//      } else if (STRAIN_operatorParams_.methodName_ == "vectorPreisach_Sutor") {
//        if (STRAIN_operatorParams_.evalVersion_ == 1) {
//          STRAIN_operatorParams_.isClassical_ = true; // original vector preisach model -> sutor2012
//
//          hystStrainTMP = new VectorPreisachSutor_ListApproach(1, STRAIN_operatorParams_, STRAIN_weightParams_, dim_, isVirgin);
////                  (1, STRAIN_operatorParams_.inputSat_, STRAIN_operatorParams_.outputSat_,
////                  STRAIN_weightParams_.weightTensor_, STRAIN_operatorParams_.rotResistance_, dim_, isVirgin,
////                  STRAIN_operatorParams_.isClassical_, STRAIN_operatorParams_.scaleUpToSaturation_,
////                  STRAIN_operatorParams_.angularDistance_,STRAIN_operatorParams_.angularResolution_,
////                  STRAIN_weightParams_.anhysteretic_a_,
////                  STRAIN_weightParams_.anhysteretic_b_, STRAIN_weightParams_.anhysteretic_c_,STRAIN_weightParams_.anhystOnly_);
//        } else if (STRAIN_operatorParams_.evalVersion_ == 2) {
//          STRAIN_operatorParams_.isClassical_ = false; // revised vector preisach model -> sutor2015
//
//          hystStrainTMP = new VectorPreisachSutor_ListApproach(1, STRAIN_operatorParams_, STRAIN_weightParams_, dim_, isVirgin);
////                  (1, STRAIN_operatorParams_.inputSat_, STRAIN_operatorParams_.outputSat_,
////                  STRAIN_weightParams_.weightTensor_, STRAIN_operatorParams_.rotResistance_, dim_, isVirgin,
////                  STRAIN_operatorParams_.isClassical_, STRAIN_operatorParams_.scaleUpToSaturation_,
////                  STRAIN_operatorParams_.angularDistance_,STRAIN_operatorParams_.angularResolution_,
////                  STRAIN_weightParams_.anhysteretic_a_,
////                  STRAIN_weightParams_.anhysteretic_b_, STRAIN_weightParams_.anhysteretic_c_,STRAIN_weightParams_.anhystOnly_);
//        } else if (STRAIN_operatorParams_.evalVersion_ == 10) {
//          STRAIN_operatorParams_.isClassical_ = true; // original vector preisach model -> sutor2015; matrix based implementation
//
//          hystStrainTMP = new VectorPreisachSutor_MatrixApproach(1, STRAIN_operatorParams_, STRAIN_weightParams_, dim_, isVirgin);
////                  (1, STRAIN_operatorParams_.inputSat_, STRAIN_operatorParams_.outputSat_,
////                  STRAIN_weightParams_.weightTensor_, STRAIN_operatorParams_.rotResistance_, dim_, isVirgin,
////                  STRAIN_operatorParams_.isClassical_, STRAIN_operatorParams_.scaleUpToSaturation_,
////                  STRAIN_operatorParams_.angularDistance_,STRAIN_operatorParams_.angularResolution_,
////                  STRAIN_weightParams_.anhysteretic_a_,
////                  STRAIN_weightParams_.anhysteretic_b_, STRAIN_weightParams_.anhysteretic_c_,STRAIN_weightParams_.anhystOnly_);
//        } else if (STRAIN_operatorParams_.evalVersion_ == 20) {
//          STRAIN_operatorParams_.isClassical_ = false; // revised vector preisach model -> sutor2015; matrix based implementation
//
//          hystStrainTMP = new VectorPreisachSutor_MatrixApproach(1, STRAIN_operatorParams_, STRAIN_weightParams_, dim_, isVirgin);
////                  (1, STRAIN_operatorParams_.inputSat_, STRAIN_operatorParams_.outputSat_,
////                  STRAIN_weightParams_.weightTensor_, STRAIN_operatorParams_.rotResistance_, dim_, isVirgin,
////                  STRAIN_operatorParams_.isClassical_, STRAIN_operatorParams_.scaleUpToSaturation_,
////                  STRAIN_operatorParams_.angularDistance_,STRAIN_operatorParams_.angularResolution_,
////                  STRAIN_weightParams_.anhysteretic_a_,
////                  STRAIN_weightParams_.anhysteretic_b_, STRAIN_weightParams_.anhysteretic_c_,STRAIN_weightParams_.anhystOnly_);
//        }
//      } else if (STRAIN_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") {
//
//          hystStrainTMP = new VectorPreisachMayergoyz(1, STRAIN_operatorParams_, STRAIN_weightParams_, dim_, isVirgin);
////                  (1, STRAIN_operatorParams_.numDirections_,
////                STRAIN_operatorParams_.inputSat_, STRAIN_operatorParams_.outputSat_,
////                STRAIN_weightParams_.weightTensor_,dim_,isVirgin,
////                STRAIN_weightParams_.anhysteretic_a_, STRAIN_weightParams_.anhysteretic_b_, STRAIN_weightParams_.anhysteretic_c_,
////                STRAIN_weightParams_.anhystOnly_,STRAIN_operatorParams_.outputClipping_);
//
//      }
//    }
//
//
//    if(outputIrrStrains){
//      results_s.open(results_name_s.str());
//      results_s << "#Number S_irr_xx S_irr_yy S_irr_xy" << std::endl;
//
//      results_xps.open(results_name_xps.str());
//      results_xps << "#Number xIn yIn Px Py S_irr_xx S_irr_yy S_irr_xy" << std::endl;
//    }
//
//		if(printStatistics){
//			statistics << "PARAMETER: " << std::endl;
//			statistics << "- xSAT: " << POL_operatorParams_.inputSat_ << std::endl;
//			statistics << "- pSAT: " << POL_operatorParams_.outputSat_ << std::endl;
//			statistics << "- anhyst a: " << POL_weightParams_.anhysteretic_a_ << std::endl;
//			statistics << "- anhyst b: " << POL_weightParams_.anhysteretic_b_ << std::endl;
//			statistics << "- anhyst c: " << POL_weightParams_.anhysteretic_c_ << std::endl;
//			statistics << "- only anhyst? " << POL_weightParams_.anhystOnly_ << std::endl;
//			statistics << "PREISACH WEIGHTS: " << std::endl;
//
//      if(POL_weightParams_.weightType_ == "Constant"){
//        statistics << "> constant = " << POL_weightParams_.constWeight_ << std::endl;
//      } else if(POL_weightParams_.weightType_ == "muDat"){
//        statistics << "> muDat with " << std::endl;
//        statistics << "- A = " << POL_weightParams_.muDat_A_ << std::endl;
//        statistics << "- sigma = " << POL_weightParams_.muDat_sigma1_ << std::endl;
//        statistics << "- h = " << POL_weightParams_.muDat_h1_ << std::endl;
//        statistics << "- eta = " << POL_weightParams_.muDat_eta_ << std::endl;
//      } else if(POL_weightParams_.weightType_ == "muDatExtended"){
//        statistics << "> extended muDat with " << std::endl;
//        statistics << "- A = " << POL_weightParams_.muDat_A_ << std::endl;
//        statistics << "- sigma1 = " << POL_weightParams_.muDat_sigma1_ << std::endl;
//        statistics << "- sigma2 = " << POL_weightParams_.muDat_sigma2_ << std::endl;
//        statistics << "- h1 = " << POL_weightParams_.muDat_h1_ << std::endl;
//        statistics << "- h2 = " << POL_weightParams_.muDat_h2_ << std::endl;
//        statistics << "- eta = " << POL_weightParams_.muDat_eta_ << std::endl;
//      } else if(POL_weightParams_.weightType_ == "muLorentz"){
//        statistics << "> Lorentzian weight function with " << std::endl;
//        statistics << "- A = " << POL_weightParams_.muLorentz_A_ << std::endl;
//        statistics << "- sigma1 = " << POL_weightParams_.muLorentz_sigma1_ << std::endl;
//        statistics << "- sigma2 = " << POL_weightParams_.muLorentz_sigma2_ << std::endl;
//        statistics << "- h1 = " << POL_weightParams_.muLorentz_h1_ << std::endl;
//        statistics << "- h2 = " << POL_weightParams_.muLorentz_h2_ << std::endl;
//      } else if(POL_weightParams_.weightType_ == "givenTensor"){
//        statistics << "> given tensor = " << std::endl;
//        statistics << POL_weightParams_.weightTensor_.ToString() << std::endl;
//      }
//      if (POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") {
//        if(POL_weightsAlreadyAdapted_ == 1){
//          statistics << "> weights directly used for Mayergoyz model" << std::endl;
//        } else {
//          statistics << "> weights were transferred for usage in Mayergoyz model" << std::endl;
//        }
//      }
//
////
////			if ( (POL_operatorParams_.methodName_ == "vectorPreisach_Sutor")||(POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") ) {
////				statistics << "LEVENBERG-MARQUARDT: " << std::endl;
////				statistics << "- max number of iterations: " << InversionParams_.maxNumIts << std::endl;
////				statistics << "- tolerance wrt x: " << InversionParams_.tolH << std::endl;
////				statistics << "- tolerance wrt y: " << InversionParams_.tolB << std::endl;
////				statistics << "- FD-resolution for Jacobian: " << InversionParams_.jacRes << std::endl;
////				statistics << "LINESEARCH: " << std::endl;
////				statistics << "- alpha start: " << InversionParams_.alphaLSStart << std::endl;
////				statistics << "- alpha min: " << InversionParams_.alphaLSMin << std::endl;
////				statistics << "- alpha max: " << InversionParams_.alphaLSMax << std::endl;
////			}
//		}
//
//		/*
//     * 2. Perform tests
//     */
//		forwardTimer = NULL;
//		backwardTimer = NULL;
//		if (measurePerformance) {
//			forwardTimer = new Timer();
//			backwardTimer = new Timer();
//			performance << "+++ RUNTIMES +++" << std::endl;
//			performance << "TEST: " << name << std::endl;
//			performance << "MODEL: " << POL_operatorParams_.methodName_ << std::endl;
//			if (POL_operatorParams_.methodName_ == "vectorPreisach_Sutor") {
//				if (POL_operatorParams_.evalVersion_ == 1) {
//					performance << "- classical model, list based implementation" << std::endl;
//				} else if (POL_operatorParams_.evalVersion_ == 2) {
//					performance << "- revised model, list based implementation" << std::endl;
//				} else if (POL_operatorParams_.evalVersion_ == 10) {
//					performance << "- classical model, matrix based implementation" << std::endl;
//				} else if (POL_operatorParams_.evalVersion_ == 20) {
//					performance << "- revised model, matrix based implementation" << std::endl;
//				}
//			}
//			performance << "SINGLE TESTS: " << std::endl;
//
//		}
//
//		if(writeResultsToFile){
//      //      std::stringstream results;
//      //			results << "# +++ RESULTS +++" << std::endl;
//      //			results << "# TEST: " << name << std::endl;
//      //			results << "# MODEL: " << POL_operatorParams_.methodName_ << std::endl;
//      //			if (POL_operatorParams_.methodName_ == "vectorPreisach_Sutor") {
//      //				if (POL_operatorParams_.evalVersion_ == 1) {
//      //					results << "# - classical model, list based implementation" << std::endl;
//      //				} else if (POL_operatorParams_.evalVersion_ == 2) {
//      //					results << "# - revised model, list based implementation" << std::endl;
//      //				} else if (POL_operatorParams_.evalVersion_ == 10) {
//      //					results << "# - classical model, matrix based implementation" << std::endl;
//      //				} else if (POL_operatorParams_.evalVersion_ == 20) {
//      //					results << "# - revised model, matrix based implementation" << std::endl;
//      //				}
//      //			}
//      //
//      //      results_x << results.str();
//      //      results_p << results.str();
//      //      results_y << results.str();
//      //      angularResults_x << results.str();
//      //      angularResults_p << results.str();
//      //      angularResults_y << results.str();
//
//			if(testInversion){
//        results_x << "#Number xIn[0] xIn[1] xOut[0] xOut[1]" << std::endl;
//        results_x << "# > xIn[0],xIn[1] = x and y components of input; given by test signal" << std::endl;
//        results_x << "# > xOut[0],xOut[1] = x and y components of inversion output" << std::endl;
//
//        results_p << "#Number pIn[0] pIn[1] pOut[0] pOut[1]" << std::endl;
//        results_p << "# > pIn[0],pIn[1] = x and y components of hyst operator computed from xIn" << std::endl;
//        results_p << "# > pOut[0],pOut[1] = x and y components of hyst operator computed from xOut" << std::endl;
//
//        results_y << "#Number yIn[0] yIn[1] yOut[0] yOut[1]" << std::endl;
//				results_y << "# > yIn[0],yIn[1] = x and y components of output of hyst operator to xIn; used as input for inversion" << std::endl;
//        results_y << "# > yOut[0],yOut[1] = x and y components computed from xOut" << std::endl;
//
//        angularResults_x << "#Number abs(xIn) atan2(xIn[1]/xIn[0])*180/pi abs(xOut) atan2(xOut[1]/xOut[0])*180/pi" << std::endl;
//        angularResults_x << "# > xIn[0],xIn[1] = x and y components of input; given by test signal" << std::endl;
//        angularResults_x << "# > xOut[0],xOut[1] = x and y components of inversion output" << std::endl;
//
//        angularResults_p << "#Number abs(pIn) atan2(pIn[1]/pIn[0])*180/pi abs(pOut) atan2(pOut[1]/pOut[0])*180/pi" << std::endl;
//        angularResults_p << "# > pIn[0],pIn[1] = x and y components of hyst operator computed from xIn" << std::endl;
//        angularResults_p << "# > pOut[0],pOut[1] = x and y components of hyst operator computed from xOut" << std::endl;
//
//        angularResults_y << "#Number abs(yIn) atan2(yIn[1]/yIn[0])*180/pi abs(yOut) atan2(yOut[1]/yOut[0])*180/pi" << std::endl;
//        angularResults_y << "# > yIn[0],yIn[1] = x and y components of output of hyst operator to xIn; used as input for inversion" << std::endl;
//        angularResults_y << "# > yOut[0],yOut[1] = x and y components computed from xOut" << std::endl;
//
//			} else {
//        results_x << "#Number xIn[0] xIn[1]" << std::endl;
//        results_x << "# > xIn[0],xIn[1] = x and y components of input; given by test signal" << std::endl;
//
//        results_p << "#Number pOut[0] pOut[1]" << std::endl;
//        results_p << "# > pOut[0],pOut[1] = x and y components of hyst operator computed from xOut" << std::endl;
//
//        results_xp << "#Number\t xIn[0]\t xIn[1]\t pOut[0]\t pOut[1]" << std::endl;
//        
//        results_y << "#Number yOut[0] yOut[1]" << std::endl;
//        results_y << "# > yOut[0],yOut[1] = x and y components computed from xOut" << std::endl;
//
//        angularResults_x << "#Number abs(xIn) atan2(xIn[1]/xIn[0])*180/pi" << std::endl;
//        angularResults_x << "# > xIn[0],xIn[1] = x and y components of input; given by test signal" << std::endl;
//
//        angularResults_p << "#Number abs(pOut) atan2(pOut[1]/pOut[0])*180/pi" << std::endl;
//        angularResults_p << "# > pOut[0],pOut[1] = x and y components of hyst operator computed from xOut" << std::endl;
//
//        angularResults_y << "#Number abs(yOut) atan2(yOut[1]/yOut[0])*180/pi" << std::endl;
//        angularResults_y << "# > yOut[0],yOut[1] = x and y components computed from xOut" << std::endl;
//
//			}
//
//      orientationTowardsExcitation_p << "# Projetion of p_x and p_y to p_perpendicular and p_parallel" << std::endl;
//      orientationTowardsExcitation_p << "# with p_parallel = p along excitation axis " << std::endl;
//      orientationTowardsExcitation_p << "# and p_perpendicular = p perpendicualr to excitation axis " << std::endl;
//      orientationTowardsExcitation_p << "# Columns: " << std::endl;
//      orientationTowardsExcitation_p << "# 1. Testnumber " << std::endl;
//      orientationTowardsExcitation_p << "# 2. Excitation amplitude " << std::endl;
//      orientationTowardsExcitation_p << "# 3. Excitation angle towards x-axis " << std::endl;
//      orientationTowardsExcitation_p << "# 4. p parallel to excitation " << std::endl;
//      orientationTowardsExcitation_p << "# 5. p perpendicular to excitation " << std::endl;
//		}
//
//
//		bool overwriteMemory = false;
//		// always overwrite Direction; deprecated flag; remove in future
//    //		bool overwriteDirection = true;
//
//		if( (printStatistics) ){
//			std::cout << "##### STARTING TEST " << name << " #####" << std::endl;
//		}
//		hOutOld.Init();
//    hOutOldForStrains.Init();
//
//    Vector<Double> projectionDir = Vector<Double>(dim_);
//    if(POL_operatorParams_.fixDirection_.NormL2() == 0){
////      std::cout << "No direction specified; using default x-direction for 1d test" << std::endl;
//      // take x-direction as default case
//      projectionDir.Init();
//      projectionDir[0] = 1.0;
//    } else {
//      projectionDir = POL_operatorParams_.fixDirection_;
//    }
//
//
//    Double eps_mu_scal;
//    Vector<Double> tmp = Vector<Double>(dim_);
//    eps_mu_base_.Mult(projectionDir,tmp);
//    tmp.Inner(projectionDir,eps_mu_scal);
//
////    std::cout << "eps_mu_scal: " << eps_mu_scal << std::endl;
////    std::cout << "eps_nu_base_[0][0]: " << eps_nu_base_[0][0] << std::endl;
////
//		for(UInt i = 0; i < totalSteps; i++){
//			if( (i%10 == 0)&&(printStatistics) ){
//				std::cout << "STEP NR " << i+1 << "/" << totalSteps << " #####" << std::endl;
//			}
//      performance << "Test " << i+1 << std::endl;
//
//			// set previous state
//			if(i > 0){
//        xInPrev = xIn;
//				xPrev = xOut; // xOut = retrieved output from inversion; only set and used in case of testInversion
//				hPrev = hOut; // hOut = output of hyst operator from forward computation with overwrite = true
//				yPrev = yOut; // yOut = calculated y = h + eps*x from computation with overwrite = true
//			} else {
//        xInPrev.Init();
//				xPrev.Init();
//				hPrev.Init();
//				yPrev.Init();
//			}
//
//			xIn.Init();
//			xIn[0] = xVals[i];
//      xIn[1] = yVals[i];
//
//           // new usage of test1D
//      // if test1D is set, project input onto POL_operatorParams_.fixDirection_
//      // if this vector is empty, use x-axis as previously
//			if(test1D == true){
////        std::cout << "Testing only along fixDirection = " << POL_operatorParams_.fixDirection_.ToString() << std::endl;
////        std::cout << "Demanded input = " << xIn.ToString() << std::endl;
//        Double projection;
//        projectionDir.Inner(xIn,projection);
//        xIn.Init();
//        xIn.Add(projection,projectionDir);
////        std::cout << "Used input = " << xIn.ToString() << std::endl;
//			}
//
//      xInBak = xIn;
//
//      if(testInversion){
//				/*
//         * Inversion test:
//         *	1. evaluate forward hyst operator with given input; overwrite memory = false
//         *  2. retrieve input by passing result of 2 to inverse hyst operator
//         *			> measure backward
//         *  3. pass retrieved input to forward hyst operator; overwrite memroy = true
//         *			> measure forward time
//         */
//
//				// 1. forward; do not overwrite
//				overwriteMemory = false;
//
//				successFlagForward = -1;
//				hIn.Init();
////        if(vector){
//          hIn = hystTMP->computeValue_vec(xIn, 0, overwriteMemory, debugOut, successFlagForward);
////        } else {
////          Double hscal = hystTMP->computeValueAndUpdate(xIn[0], 0, overwriteMemory, successFlagForward);
////          hIn.Add(hscal,projectionDir);
////        }
//
//				yIn.Init();
//        yIn.Add(1.0,hIn);
//        for(UInt j = 0; j < dim_; j++){
//          yIn[j] += eps_mu_base_[j][j]*xIn[j];
//        }
////
////        std::cout << "Computed hIn = " << hIn.ToString() << std::endl;
////        std::cout << "Computed yIn = " << yIn.ToString() << std::endl;
//
//				// 2. backward; do not overwrite
//				overwriteMemory = false;
//        successFlagBackward = -1;
//				numberOfLMIterations = 0;
//				numberOfLinesearchIterations = 0;
//				maxNumberOfLinesearchIterations = 0;
//
//        Double scalIn = 0.0;
//        projectionDir.Inner(yIn,scalIn);
////        Double scalIn = yIn[0];
//        Double scalOut = 0;
//
//        if (measurePerformance) {
//          backwardTimer->Start();
//					startTime = backwardTimer->GetCPUTime();
//        }
//
//        if(InversionParams_.inversionMethod != 4){
//          xOut = hystTMP->computeInput_vec_withStatistics(yIn, yPrev, xPrev, hPrev,
//                  0, eps_mu_base_, POL_operatorParams_.fieldsAlignedAboveSat_, POL_operatorParams_.hystOutputRestrictedToSat_,
//                  numberOfLMIterations, numberOfLinesearchIterations, maxNumberOfLinesearchIterations,
//                  successFlagBackward, minAlpha, maxAlpha, avgAlpha, xIn);
//        } else {
////					xOut[0] = hystTMP->computeInputAndUpdate(yIn[0], eps_mu_base_[0][0],
////                  0, overwriteMemory, successFlagBackward);
//          scalOut = hystTMP->computeInputAndUpdate(scalIn, eps_mu_scal,
//                  0, overwriteMemory, successFlagBackward);
//        }
//
//        if (measurePerformance) {
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
//
//						performance << "- backward: " << evalTime << std::endl;
//          } else {
//            performance << "- backward: " << evalTime << "(reused)" << std::endl;
//					}
//        }
//
//        if(InversionParams_.inversionMethod == 4){
//          xOut.Init();
//          xOut.Add(scalOut,projectionDir);
////          xOut[0] = scalOut;
//        }
//
////        std::cout << "Computed xOut = " << xOut.ToString() << std::endl;
////
//				if(successFlagBackward == -1){
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
//        successCodeVectorBackward[i] = successFlagBackward;
//
//				if(minAlpha < totalminAlpha){
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
//					absmaxNumberOfLinesearchIterations = maxNumberOfLinesearchIterations;
//				}
//
//				/*
//         * Check result
//         */
//				xErr.Init();
//				xErr.Add(-1.0,xIn,1.0,xOut);
//				// check below after eval of step 3
//
//				// set xIn to xOut, so that during forward evaluation, xOut is used
//				xIn = xOut;
//
//			} // testInversion
//
//			// 3. forward; overwrite this time and measure performance
//      // NOTE: in case of inversion, we set the memory of the system with the RETRIEVED SOLUTION
//      // i.e. if inversion fails, it will also affect the reference solution!
//			overwriteMemory = true;
//
//			successFlagForward = -1;
//			hOut.Init();
//      hOutForStrains.Init();
//
//      Double hinScal = 0.0;
//      projectionDir.Inner(xIn,hinScal);
////      Double houtScal = 0.0;
//
//      if ((POL_operatorParams_.printOut_ > 0)) {
//        static UInt cntTS = 1;
//        hystTMP->collectParallelProjections(true, cntTS);
//        cntTS++;
//      }
//      
//			if (measurePerformance) {
//				forwardTimer->Start();
//				startTime = forwardTimer->GetCPUTime();
//			}
//
////			if(vector){
//				hOut = hystTMP->computeValue_vec(xIn, 0, overwriteMemory, debugOut, successFlagForward);
////			} else {
////				houtScal = hystTMP->computeValueAndUpdate(hinScal, 0, overwriteMemory, successFlagForward);
////			}
//
//			if (measurePerformance) {
//				forwardTimer->Stop();
//
//				// only count actual computations; if value gets reused, this should not count
//				if(successFlagForward != 0) {
//					forwardEvalCounter++;
//
//					endTime = forwardTimer->GetCPUTime();
//					evalTime = endTime-startTime;
//					if(evalTime > forwardMaxEvalTime){
//						forwardMaxEvalTime = evalTime;
//					}
//					forwardTotalEvalTime += evalTime;
//
//					performance << "- forward: " << evalTime << std::endl;
//				} else {
//					performance << "- forward: " << evalTime << "(reused)" << std::endl;
//				}
//			}
//
//      
//
//      if ((POL_operatorParams_.printOut_ > 0)) {
//        static UInt cnt = 1;
//        bool appendToTxt = false;
//        if(cnt > 1){
//          appendToTxt = true;
//        } 
//        
//        if (((cnt-1) % POL_operatorParams_.printOut_ == 0)) {
//          std::stringstream filenamebuf;
//          filenamebuf << "PreisachPlane_Step_" << std::setfill('0') << std::setw(5) << cnt << ".bmp";
////            filenamebuf << "PreisachPlane_Step_" << std::setfill('0') << std::setw(5) << cnt << "_v" << POL_operatorParams_.evalVersion_ << "_numRows" << POL_weightParams_.numRows_ << ".bmp";
//          std::stringstream filenamebuf2;
//          filenamebuf2 << "NestedList.txt";
//          std::stringstream optionalHeaderStream;
//          if(cnt == 1){
//            optionalHeaderStream << "Parameter of vector Preisach model based on rotational operators: \n";
//            if(POL_operatorParams_.isClassical_){
//              optionalHeaderStream << "+++ version: classical \n";
//              optionalHeaderStream << "+++ rotational resistance: " << POL_operatorParams_.rotResistance_ << "\n";
//            } else {
//              optionalHeaderStream << "+++ version: revised \n";
//              optionalHeaderStream << "+++ rotational resistance: " << POL_operatorParams_.rotResistance_ << "\n";
//              optionalHeaderStream << "+++ angular distance: " << POL_operatorParams_.angularDistance_ << "\n"; 
//            }
//            optionalHeaderStream << "\n";
//          }
//          optionalHeaderStream << "Step: " << std::setfill('0') << std::setw(5) << cnt << "\n";
//          optionalHeaderStream << "+++ Normalized input vector: (" << xIn[0]/POL_operatorParams_.inputSat_ << ", " << xIn[1]/POL_operatorParams_.inputSat_ << ") \n";
//
//          if(cnt == 1){
//            std::cout << "--- Printing state of Preisach plane to " << filenamebuf.str() << std::endl;
//            std::cout << "--- Printing nested list to " << filenamebuf2.str() << std::endl;            
//          }
//
//          hystTMP->switchingStateToBmp(POL_operatorParams_.bmpResolution_, filenamebuf.str(), 0, true);  
//          hystTMP->rotationListToTxt(filenamebuf2.str(), 0, appendToTxt, optionalHeaderStream.str());
//          // deactivate collection to save runtime
//          hystTMP->collectParallelProjections(false, 0);
//        }
//        cnt++;
//      }
//        
////      if(!vector){
////        hOut.Init();
////        hOut.Add(houtScal,projectionDir);
////      }
//
//      if(CouplingParams_.ownHystOperator_){
//        if(vector){
//          hOutForStrains = hystStrainTMP->computeValue_vec(xIn, 0, overwriteMemory, debugOut, successFlagForward);
//        } else {
//          hOutForStrains.Init();
//          Double scalOut = hystStrainTMP->computeValueAndUpdate(hinScal, 0, overwriteMemory, successFlagForward);
//          hOutForStrains.Add(scalOut,projectionDir);
//        }
//      } else {
//        hOutForStrains = hOut;
//      }
//
//			if(successFlagForward == -1){
//				forwardFails++;
//			} else if(successFlagForward == 0){
//				forwardReused++;
//			} else if(successFlagForward == 1){
//				forwardAnhystOnly++;
//			} else if(successFlagForward == 2){
//				forwardOverwrite++;
//			} else if(successFlagForward == 3){
//				forwardTMP++;
//			}
//			successCodeVectorForward[i] = successFlagForward;
//
//			yOut.Init();
//			yOut.Add(1.0,hOut);
//			for(UInt j = 0; j < dim_; j++){
//        // note that in case of inversion, xIn is the RETRIEVED input!
//				yOut[j] += eps_mu_base_[j][j]*xIn[j];
//			}
//
//      if(outputIrrStrains){
//        Vector<Double> S_irr = ComputeIrreversibleStrains(hOutForStrains,xIn,hOutOldForStrains);
//        results_s << i+1 << " " << std::setprecision(precisionDigits) << S_irr.ToString() << std::endl;
//        results_xps << i+1 << " " << std::setprecision(precisionDigits) << xIn.ToString() << " " << hOutForStrains.ToString() << " " << S_irr.ToString() << std::endl;
//      }
//      if(writeResultsToFile){
//        //results_xp << i+1 << " " << std::setprecision(precisionDigits) << xIn.ToString() << " " << hOut.ToString() << std::endl;
//        results_xp << i+1 << " " << std::setprecision(precisionDigits) << xIn[0] << " " << xIn[1] << " " << hOut[0] << " " << hOut[1] << std::endl;
//        
//        Double xAmpl = xIn.NormL2();
//        Double xAngle = std::atan2(xIn[1],xIn[0])/M_PI*180;
//
//        Vector<Double> xDir = Vector<Double>(dim_);
//        xDir.Init();
//        if(xAmpl != 0){
//          xDir.Add(1.0/xAmpl,xIn);
//        } else {
//          Double xAmplOld = xInPrev.NormL2();
//          if(xAmplOld != 0){
//            xDir.Add(1.0/xAmplOld,xInPrev);
//          }
//        }
//        Double pParallel = hOut.Inner(xDir);
//        Vector<Double> pPerpendicularVec = Vector<Double>(dim_);
//        pPerpendicularVec.Add(1.0,hOut,-pParallel,xDir);
//        Double pPerpendicular = pPerpendicularVec.NormL2();
//
//        orientationTowardsExcitation_p << i+1 << " " << std::setprecision(precisionDigits) << xAmpl << " " << xAngle << " " << pParallel << " " << pPerpendicular << std::endl;
//      }
//			/*
//       * 4. Check result
//       * > only for inversion test; for forward test, we do not have a
//       *		value for comparison
//       */
//			bool failWRTX = false;
//			bool failWRTY = false;
//			if(testInversion){
//				yErr.Init();
//				yErr.Add(-1.0,yIn,1.0,yOut);
//
//				if(xErr.NormL2() > InversionParams_.tolH){
//					failWRTX = true;
//					failedTests_xIn[i] = xInBak;
//					failedTests_xOut[i] = xOut;
//					failedTests_yIn[i] = yIn;
//					failedTests_yOut[i] = yOut;
//				}
//
//				if(yErr.NormL2() > InversionParams_.tolB){
//					failWRTY = true;
//					failedTests_xIn[i] = xInBak;
//					failedTests_xOut[i] = xOut;
//					failedTests_yIn[i] = yIn;
//					failedTests_yOut[i] = yOut;
//				}
//
//				if( failWRTX || failWRTY ){
//					failedTests_wrtX[i] = std::pair< bool, Double >(failWRTX, xErr.NormL2());
//					failedTests_wrtY[i] = std::pair< bool, Double >(failWRTY, yErr.NormL2());
//				}
//
//        if( failWRTX && !failWRTY ){
//          numHalfFails++;
//        }
//
//        if( failWRTX && failWRTY ){
//          numFails++;
//        }
//			}
//
//			/*
//       * 5. Output and statistics
//       */
//      std::string delimiter = " ";
//			if(writeResultsToFile){
//        Double angle1 = 0;
//        Double ampl1 = 0;
//        Double angle2 = 0;
//        Double ampl2 = 0;
//				if(testInversion){
//					//results << "# Number \t xIn[0] \t xIn[1] \t yIn[0] \t yIn[1] \t xOut[0] \t xOut[1] \t yOut[0] \t yOut[1]" << std::endl;
//					results_x << i+1 << std::setprecision(precisionDigits) << delimiter << xInBak[0] << delimiter << xInBak[1] << delimiter << xOut[0] << delimiter << xOut[1] << std::endl;
//          ampl1 = xIn.NormL2();
//          angle1 = atan2(xInBak[1],xInBak[0])*180.0/M_PI;
//          ampl2 = xOut.NormL2();
//          angle2 = atan2(xOut[1],xOut[0])*180.0/M_PI;
//          angularResults_x << i+1 << std::setprecision(precisionDigits) << delimiter << ampl1 << delimiter << angle1 << delimiter << ampl2 << delimiter << angle2 << std::endl;
//
//          results_p << i+1 << std::setprecision(precisionDigits) << delimiter << hIn[0] << delimiter << hIn[1] << delimiter << hOut[0] << delimiter << hOut[1] << std::endl;
//          ampl1 = hIn.NormL2();
//          angle1 = atan2(hIn[1],hIn[0])*180.0/M_PI;
//          ampl2 = hOut.NormL2();
//          angle2 = atan2(hOut[1],hOut[0])*180.0/M_PI;
//          angularResults_p << i+1 << std::setprecision(precisionDigits) << delimiter << ampl1 << delimiter << angle1 << delimiter << ampl2 << delimiter << angle2 << std::endl;
//
//          results_y << i+1 << std::setprecision(precisionDigits) << delimiter << yIn[0] << delimiter << yIn[1] << delimiter << yOut[0] << delimiter << yOut[1] << std::endl;
//          ampl1 = yIn.NormL2();
//          angle1 = atan2(yIn[1],yIn[0])*180.0/M_PI;
//          ampl2 = yOut.NormL2();
//          angle2 = atan2(yOut[1],yOut[0])*180.0/M_PI;
//          angularResults_y << i+1 << std::setprecision(precisionDigits) << delimiter << ampl1 << delimiter << angle1 << delimiter << ampl2 << delimiter << angle2 << std::endl;
//				} else {
//					//					results << "# Number \t xIn[0] \t xIn[1] \t yOut[0] \t yOut[1]" << std::endl;
//          results_x << i+1 << std::setprecision(precisionDigits) << delimiter << xInBak[0] << delimiter << xInBak[1] << std::endl;
//          ampl1 = xIn.NormL2();
//          angle1 = atan2(xInBak[1],xInBak[0])*180.0/M_PI;
//          angularResults_x << i+1 << std::setprecision(precisionDigits) << delimiter << ampl1 << delimiter << angle1 << std::endl;
//
//          results_p << i+1 << std::setprecision(precisionDigits) << delimiter << hOut[0] << delimiter << hOut[1] << std::endl;
//          ampl1 = hOut.NormL2();
//          angle1 = atan2(hOut[1],hOut[0])*180.0/M_PI;
//          angularResults_p << i+1 << std::setprecision(precisionDigits) << delimiter << ampl1 << delimiter << angle1 << std::endl;
//
//          results_y << i+1 << std::setprecision(precisionDigits) << delimiter << yOut[0] << delimiter << yOut[1] << std::endl;
//          ampl1 = yOut.NormL2();
//          angle1 = atan2(yOut[1],yOut[0])*180.0/M_PI;
//          angularResults_y << i+1 << std::setprecision(precisionDigits) << delimiter << ampl1 << delimiter << angle1 << std::endl;
//				}
//			}
//
//			LOG_TRACE(coeffunctionhyst_main) << "##### TARGET X-VECTOR #####";
//			LOG_TRACE(coeffunctionhyst_main) << xInBak.ToString();
//			if(testInversion){
//				LOG_TRACE(coeffunctionhyst_main) << "##### RETRIEVED X-VECTOR #####";
//				LOG_TRACE(coeffunctionhyst_main) << xOut.ToString();
//
//				LOG_TRACE(coeffunctionhyst_main) << "##### ERROR VECTOR wrt X #####";
//				LOG_TRACE(coeffunctionhyst_main) << xErr.ToString();
//				LOG_TRACE(coeffunctionhyst_main) << "##### TARGET Y-VECTOR #####";
//				LOG_TRACE(coeffunctionhyst_main) << yIn.ToString();
//			}
//			LOG_TRACE(coeffunctionhyst_main) << "##### RETRIEVED Y-VECTOR #####";
//			LOG_TRACE(coeffunctionhyst_main) << yOut.ToString();
//			if(testInversion){
//				LOG_TRACE(coeffunctionhyst_main) << "##### ERROR VECTOR wrt Y #####";
//				LOG_TRACE(coeffunctionhyst_main) << yErr.ToString();
//			}
//			hOutOld = hOut;
//      hOutOldForStrains = hOutForStrains;
//		}
//
//		UInt LMcases = totalSteps-totalReused-totalBisection-totalAnhystOnly-totalRemanence;
//		if(LMcases != 0){
//			totalavgAlpha /= (Double) LMcases;
//		}
//
//		if(printStatistics){
//			std::stringstream majorResults;
//			majorResults << "############################# " <<	std::endl;
//			majorResults << "##### RESULTS FOR TEST " << name <<	std::endl;
//			majorResults << " " << totalSteps-numFails << " of " << totalSteps << " satisfied at least the failback criterion, i.e. |residual_Y| = |Y - mu_eps*X - P(X)| < " << InversionParams_.tolB << std::endl;
//			majorResults << " " << numHalfFails << " of " << totalSteps << " failed error criterion but passed due to failback criterion" << std::endl;
//      majorResults << " " << numFails << " of " << totalSteps << " failed to satisfy even the failback criterion" << std::endl;
//
//			// print major results to stdout
//			std::cout << majorResults.str() << std::endl;
//
//			statistics << majorResults.str() << std::endl;
//			if( (numFails > 0)||(numHalfFails > 0) ){
//				statistics << "## Detailed Statistics for failed tests: " << std::endl;
//
//				std::map< UInt, std::pair< bool, Double > >::iterator mapItX;
//				std::map< UInt, std::pair< bool, Double > >::iterator mapItY = failedTests_wrtY.begin();
//				bool failX, failY;
//				for(mapItX = failedTests_wrtX.begin(); mapItX != failedTests_wrtX.end(); mapItX++,mapItY++){
////          std::cout << "1" << std::endl;
//					statistics << " Test Nr " << mapItX->first << std::endl;
//					failX = mapItX->second.first;
//					failY = mapItY->second.first;
////          std::cout << "2" << std::endl;
//					bool xsolabovesat = failedTests_xIn[mapItX->first].NormL2()>POL_operatorParams_.inputSat_;
//					bool xretabovesat = failedTests_xOut[mapItX->first].NormL2()>POL_operatorParams_.inputSat_;
//					statistics << " > X_in: " << failedTests_xIn[mapItX->first].ToString()  << "; above sat? "<< xsolabovesat << std::endl;
//					statistics << " > X_out: " << failedTests_xOut[mapItX->first].ToString() << "; above sat? "<< xretabovesat << std::endl;
//					if(failX){
//						statistics << " > FAIL - remaining error wrt X  " << mapItX->second.second << std::endl;
//					} else {
//						statistics << " > SUCCESS - remaining error wrt X  " << mapItX->second.second << std::endl;
//					}
////					std::cout << "3" << std::endl;
//					statistics << " > Y_in: " << failedTests_yIn[mapItX->first].ToString() << std::endl;
//					statistics << " > Y_out: " << failedTests_yOut[mapItX->first].ToString() << std::endl;
//					if(failY){
//						statistics << " > FAIL - remaining error wrt Y  " << mapItY->second.second << std::endl;
//					} else {
//						statistics << " > SUCCESS - remaining error wrt Y  " << mapItY->second.second << std::endl;
//					}
////					std::cout << "4" << std::endl;
//					statistics << " - SuccessCode (forward): " << successCodeVectorForward[mapItX->first] << std::endl;
//					statistics << " - SuccessCode (backward): " << successCodeVectorBackward[mapItX->first] << std::endl;
//				}
//			}
//
//			statistics << "## Detailed Statistics on forward evaluations (excluding those inside of LM/during inversion): " << std::endl;
//			statistics << " " << forwardReused << " of " << totalSteps << " reused old output as input change was below tolerance" << std::endl;
//			statistics << " " << forwardAnhystOnly << " of " << totalSteps << " featured only anhysteretic parts" << std::endl;
//			statistics << " " << forwardOverwrite << " of " << totalSteps << " evaluations were performed on permanent storage" << std::endl;
//			statistics << " " << forwardTMP << " of " << totalSteps << " evaluations were performed on temporal storage" << std::endl;
//
//			statistics << "## Detailed Statistics on inversion: " << std::endl;
//			statistics << " " << totalReused << " of " << totalSteps << " reused old solution as DeltaY < " << InversionParams_.tolB << std::endl;
//			statistics << " " << totalAnhystOnly << " of " << totalSteps << " were solved using bisection for pure anhysteretic case" << std::endl;
//			statistics << " " << totalBisection << " of " << totalSteps << " were solved using bisection / simple division" << std::endl;
//			if(!vector){
//				statistics << " " << totalPassedErrorTol << " of " << totalSteps << " tests passed due to error tolerance" << std::endl;
//			}
//
//			if(vector){
//				statistics << " " << totalRemanence << " of " << totalSteps << " cases were in remanence" << std::endl;
//				statistics << " " << LMcases << " of " << totalSteps << " were solved using Levenberg-Marquardt" << std::endl;
//				if(LMcases != 0){
//					statistics << "## Detailed Statistics for Levenberg-Marquardt: " << std::endl;
//					statistics << " " << totalPassedErrorTol << " of " << LMcases << " tests passed due to error tolerance |JacT*Res| < " << InversionParams_.tolH << std::endl;
//					statistics << " " << totalPassedResTolX << " of " << LMcases << " tests passed due to |residual_X| = |X - mu_eps^-1*(Y - P(X)| < " << InversionParams_.tolH << std::endl;
//					statistics << " " << totalPassedResTolY << " of " << LMcases << " tests passed due to |residual_Y| = |Y - mu_eps*X - P(X)| < " << InversionParams_.tolH << std::endl;
//					statistics << " " << LMFails << " of " << LMcases << " failed the failback criterion (|residual_Y| < " << InversionParams_.tolB << ")" << std::endl;
//					statistics << " Total number of LM iterations: " << totalNumberOfLMIterations << std::endl;
//					statistics << " Average number of LM iterations: " << (Double) totalNumberOfLMIterations/LMcases << std::endl;
//					statistics << " Total number of Linesearch iterations: " << totalNumberOfLinesearchIterations << std::endl;
//					statistics << " Average number of Linesearch iterations (per LM Iteration): " << (Double) totalNumberOfLinesearchIterations/totalNumberOfLMIterations << std::endl;
//					statistics << " Maximal number of Linesearch iterations: " << absmaxNumberOfLinesearchIterations << std::endl;
//					statistics << " Minimal alphaLS: " << totalminAlpha << std::endl;
//					statistics << " Maximal alphaLS: " << totalmaxAlpha << std::endl;
//					statistics << " Average alphaLS: " << totalavgAlpha << std::endl;
//				}
//			}
//			statistics << "############################# " <<	std::endl;
//      statistics << std::endl;
//		}
//
//
//		if (measurePerformance){
//			Double forwardAvgEvalTime = forwardTotalEvalTime/forwardEvalCounter;
//			performance << "## Runtime information: " << std::endl;
//      performance << "# Runtime information for forward evaluations: " << std::endl;
//			performance << "Total number of forward evaluations: " << forwardEvalCounter << std::endl;
//			performance << "Total time for forward evaluations: " << forwardTotalEvalTime << std::endl;
//			performance << "Average time for forward evaluations: " << forwardAvgEvalTime << std::endl;
//			performance << "Maximal time for forward evaluations: " << forwardMaxEvalTime << std::endl;
//      performance << "> Note: only evaluations on permanent storage counted, " << std::endl;
//      performance << "   i.e. no reused values and no evaluations on temporal storage (as e.g. inside LM)" << std::endl;
//			if(testInversion){
//				Double backwardAvgEvalTime = backwardTotalEvalTime/backwardEvalCounter;
//        performance << "# Runtime information for backward evaluations: " << std::endl;
//				performance << "Total number of backward evaluations: " << backwardEvalCounter << std::endl;
//				performance << "Total time for backward evaluations: " << backwardTotalEvalTime << std::endl;
//				performance << "Average time for backward evaluations: " << backwardAvgEvalTime << std::endl;
//				performance << "Maximal time for backward evaluations: " << backwardMaxEvalTime << std::endl;
//        performance << "# Results for grepping: " << std::endl;
//        //        performance << "Flag for grep \t forwardAvgEvalTime (#evals) \t backwardAvgEvalTime (#evals) \t  numFails "<< std::endl;
//        //				performance << "TOGREP: \t\t" << forwardAvgEvalTime << " (" << forwardEvalCounter << ") " << "\t " << backwardAvgEvalTime << " (" << backwardEvalCounter << ") " << "\t " << numFails << std::endl;
//        /*
//         * NOTE: we do measure forward time when computing inversion, too, but the results will be much different
//         *        from the results without inversion; inversion seems to affect forward evaluation time, even though no storage will be overwritten
//         *        by LM (or should not) and when inspecting the return flag of forward evaluation it is not 0, i.e. result is not reused;
//         *        nevertheless the forward runtime will be as short as if results had been reused
//         *      > not sure why this is the case but as a consequence, we measure forward time and backward time during separate runs
//         */
//        performance << "Flag for grep \t backwardAvgEvalTime (#evals) \t  numFails "<< std::endl;
//				performance << "TOGREP: \t\t" << "\t " << backwardAvgEvalTime << " (" << backwardEvalCounter << ") " << "\t " << numFails << std::endl;
//			} else {
//        performance << "# Results for grepping: " << std::endl;
//        performance << "Flag for grep \t forwardAvgEvalTime (forwardEvalCounter) " << std::endl;
//        performance << "TOGREP: \t\t" << forwardAvgEvalTime << " (" << forwardEvalCounter << ") " << std::endl;
//      }
//			//statistics << "TOGREP: " << forwardAvgEvalTime << " (" << forwardEvalCounter << ") " << forwardTotalEvalTime << " " << backwardAvgEvalTime << " (" << backwardEvalCounter << ") " << " " << backwardTotalEvalTime << " " << numFails << std::endl;			performance << "############################# " <<	std::endl;
//      performance << std::endl;
//
//
//      if(commonPerformanceFile != "---"){
//        bool detailled = false;
//        if(testInversion){
//          Double backwardAvgEvalTime = backwardTotalEvalTime/backwardEvalCounter;
//
//          if(detailled){
//            commonPerfStream << std::fixed << std::setprecision(6) << nameTagForPerfFile << "\t\t Forward: " << forwardAvgEvalTime << " (Evals=" << forwardEvalCounter << ") " << "\t Backward: " << backwardAvgEvalTime << " (Evals=" << backwardEvalCounter << ") " << "\t InversionFails=" << numFails << std::endl;
//          } else {
//            commonPerfStream << std::fixed << std::setprecision(6) << nameTagForPerfFile << "\t" << forwardAvgEvalTime << "\t" << forwardEvalCounter << "\t" << backwardAvgEvalTime << "\t" << backwardEvalCounter << "\t" << numFails << std::endl;
//          }
//        } else {
//          if(detailled){
//            commonPerfStream << std::fixed << std::setprecision(6) << nameTagForPerfFile << "\t\t Forward: " << forwardAvgEvalTime << " (Evals=" << forwardEvalCounter << ") " << "\t Backward: 0.0 (Evals=0) \t InversionFails=0" << std::endl;
//          } else {
//            commonPerfStream << std::fixed << std::setprecision(6) << nameTagForPerfFile << "\t" << forwardAvgEvalTime << "\t" << forwardEvalCounter << "\t" << 0 << "\t" << 0 << "\t" << 0 << std::endl;
//          }
//        }
//      }
//		}
//
//		if(writeResultsToFile){
//			results_x.close();
//      results_xp.close();
//      results_p.close();
//      results_y.close();
//      angularResults_x.close();
//      angularResults_p.close();
//      angularResults_y.close();
//      orientationTowardsExcitation_p.close();
//		}
//    if(outputIrrStrains){
//      results_s.close();
//    }
//
//		if(printStatistics){
//			statistics.close();
//		}
//		if(measurePerformance){
//			performance.close();
//      if(commonPerformanceFile != "---"){
//        commonPerfStream.close();
//      }
//		}
//
//  }

  
  void CoefFunctionHyst::TestInversion(PtrParamNode testNode, PtrParamNode infoNode){

    if(testNode == NULL){
      EXCEPTION("No input to inversion test provided");
    }
		if(storageInitialized_ == false){
      // memory has not been set yet > no tests possible
			EXCEPTION("Memory has not been initialized yet");
    }

    bool use3dOperator = false;
    if(testNode->Has("Use3dOperator")){
      testNode->GetValue("Use3dOperator",use3dOperator,ParamNode::PASS);
    }
    UInt dimHystOperatorForTesting = dim_;
    if(use3dOperator){
      dimHystOperatorForTesting = 3;
    }
    
    int forcedPreisachResolution = -1;
    if(testNode->Has("PreisachPlaneResolution")){
      testNode->GetValue("PreisachPlaneResolution",forcedPreisachResolution,ParamNode::PASS);
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
    writeResultsToInfoXML_ = false;
    if(testNode->Has("writeResultsToInfoXML")){
      testNode->GetValue("writeResultsToInfoXML",writeResultsToInfoXML_,ParamNode::PASS);
      if(infoNode == NULL){
        EXCEPTION("Cannot print to info.xml - infoNode = NULL");
      }
      infoNodeLoc_ = infoNode;
    }
    
    bool writeResultsToFile = false;
    if(testNode->Has("WriteResultsToFile")){
      testNode->GetValue("WriteResultsToFile",writeResultsToFile,ParamNode::PASS);
    }
    bool writeInputToFile = false;
    if(testNode->Has("WriteInputToFile")){
      testNode->GetValue("WriteInputToFile",writeInputToFile,ParamNode::PASS);
    }
    
    targetDirForResultFiles_ = "history_hystOperator";
    if(testNode->Has("DirectoryForTestoutput")){
      testNode->GetValue("DirectoryForTestoutput",targetDirForResultFiles_,ParamNode::PASS);
    }
    
    bool measurePerformance = false;
    std::string commonPerfFile;
    std::string additionalTag1 = "---";
    std::string additionalTag2 = "---";
    std::string additionalTag3 = "---";
    if(testNode->Has("MeasurePerformance")){
      PtrParamNode PerformanceNode = testNode->Get("MeasurePerformance");
      if(PerformanceNode->Has("activate")){
        PerformanceNode->GetValue("activate",measurePerformance,ParamNode::PASS);
      }
      if(PerformanceNode->Has("commonResultFile")){
        PerformanceNode->GetValue("commonResultFile",commonPerfFile,ParamNode::PASS);
        commonPerfFile.erase(std::remove( commonPerfFile.begin(), commonPerfFile.end(), '\"' ),commonPerfFile.end());
      }
      if(PerformanceNode->Has("additionalTag1")){
        PerformanceNode->GetValue("additionalTag1",additionalTag1,ParamNode::PASS);
        additionalTag1.erase(std::remove( additionalTag1.begin(), additionalTag1.end(), '\"' ),additionalTag1.end());
      }
      if(PerformanceNode->Has("additionalTag2")){
        PerformanceNode->GetValue("additionalTag2",additionalTag2,ParamNode::PASS);
        additionalTag2.erase(std::remove( additionalTag2.begin(), additionalTag2.end(), '\"' ),additionalTag2.end());
      }
      if(PerformanceNode->Has("additionalTag3")){
        PerformanceNode->GetValue("additionalTag3",additionalTag3,ParamNode::PASS);
        additionalTag3.erase(std::remove( additionalTag3.begin(), additionalTag3.end(), '\"' ),additionalTag3.end());
      }
    } else {
      commonPerfFile = "---";
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
    Vector<Double> zVals = Vector<Double>(1);

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

      std::string nameTagForPerfFile = "---";
        if(InputSignals->Get("Sine")->Has("NameTagPerFile")){
        InputSignals->Get("Sine")->GetValue("NameTagPerFile",nameTagForPerfFile,ParamNode::PASS);
      }
      if(nameTagForPerfFile == "---"){
        nameTagForPerfFile = outputName;
      }

      CreatePeriodicTestSignal("Sine",amplitudeScaling,numPeriods,stepsPerPeriod,xVals,yVals);
      zVals = xVals;
      zVals.Init();
      
      /*
       * START ACTUAL TESTING ROUTINE
       */
      TestHystOperatorWithSignal(outputName,xVals,yVals,zVals,dimHystOperatorForTesting,forcedPreisachResolution,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains,nameTagForPerfFile,additionalTag1,additionalTag2,additionalTag3);
      
      if(writeInputToFile){
        WriteSignalToFile("Sine_input",xVals,yVals,zVals);
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

      Double initialAmplitude = 0.0;
      if(InputSignals->Get("Rotation")->Has("InitialAmplitudeScaling")){
        InputSignals->Get("Rotation")->GetValue("InitialAmplitudeScaling",initialAmplitude,ParamNode::PASS);
      }
      
      std::string outputName = "---";
      if(InputSignals->Get("Rotation")->Has("OutputName")){
        InputSignals->Get("Rotation")->GetValue("OutputName",outputName,ParamNode::PASS);
      }

      if(outputName == "---"){
        std::stringstream outputNameStream;
        outputNameStream << "Rot_At-"<< std::showpoint << amplitudeScaling <<"xSat_startAt-" << std::showpoint << initialAmplitude << "xSat";
        outputName = outputNameStream.str();
      }

      std::string nameTagForPerfFile = "---";
      if(InputSignals->Get("Rotation")->Has("NameTagPerFile")){
        InputSignals->Get("Rotation")->GetValue("NameTagPerFile",nameTagForPerfFile,ParamNode::PASS);
      }
      if(nameTagForPerfFile == "---"){
        nameTagForPerfFile = outputName;
      }

      CreatePeriodicTestSignal("Rotation",amplitudeScaling,numPeriods,stepsPerPeriod,xVals,yVals,initialAmplitude);
      zVals = xVals;
      zVals.Init();
      
      /*
       * START ACTUAL TESTING ROUTINE
       */
      TestHystOperatorWithSignal(outputName,xVals,yVals,zVals,dimHystOperatorForTesting,forcedPreisachResolution,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains,nameTagForPerfFile,additionalTag1,additionalTag2,additionalTag3);

      if(writeInputToFile){
        WriteSignalToFile("Rotation_input",xVals,yVals,zVals);
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

      std::string nameTagForPerfFile = "---";
        if(InputSignals->Get("DecreasingRotation")->Has("NameTagPerFile")){
        InputSignals->Get("DecreasingRotation")->GetValue("NameTagPerFile",nameTagForPerfFile,ParamNode::PASS);
      }
      if(nameTagForPerfFile == "---"){
        nameTagForPerfFile = outputName;
      }

      CreatePeriodicTestSignal("DecreasingRotation",amplitudeScaling,numPeriods,stepsPerPeriod,xVals,yVals);
      zVals = xVals;
      zVals.Init();
      
      /*
       * START ACTUAL TESTING ROUTINE
       */
      TestHystOperatorWithSignal(outputName,xVals,yVals,zVals,dimHystOperatorForTesting,forcedPreisachResolution,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains,nameTagForPerfFile,additionalTag1,additionalTag2,additionalTag3);
      
      if(writeInputToFile){
        WriteSignalToFile("DecreasingRotation_input",xVals,yVals,zVals);
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

      std::string nameTagForPerfFile = "---";
        if(InputSignals->Get("IncreasingRotation")->Has("NameTagPerFile")){
        InputSignals->Get("IncreasingRotation")->GetValue("NameTagPerFile",nameTagForPerfFile,ParamNode::PASS);
      }
      if(nameTagForPerfFile == "---"){
        nameTagForPerfFile = outputName;
      }

      CreatePeriodicTestSignal("IncreasingRotation",amplitudeScaling,numPeriods,stepsPerPeriod,xVals,yVals);
      zVals = xVals;
      zVals.Init();
      
      /*
       * START ACTUAL TESTING ROUTINE
       */
      TestHystOperatorWithSignal(outputName,xVals,yVals,zVals,dimHystOperatorForTesting,forcedPreisachResolution,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains,nameTagForPerfFile,additionalTag1,additionalTag2,additionalTag3);
      
      if(writeInputToFile){
        WriteSignalToFile("IncreasingRotation_input",xVals,yVals,zVals);
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

      std::string nameTagForPerfFile = "---";
      if(InputSignals->Get("DecreasingSine")->Has("NameTagPerFile")){
        InputSignals->Get("DecreasingSine")->GetValue("NameTagPerFile",nameTagForPerfFile,ParamNode::PASS);
      }
      if(nameTagForPerfFile == "---"){
        nameTagForPerfFile = outputName;
      }

      CreatePeriodicTestSignal("DecreasingSine",amplitudeScaling,numPeriods,stepsPerPeriod,xVals,yVals);
      zVals = xVals;
      zVals.Init();
      
      /*
       * START ACTUAL TESTING ROUTINE
       */
      TestHystOperatorWithSignal(outputName,xVals,yVals,zVals,dimHystOperatorForTesting,forcedPreisachResolution,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains,nameTagForPerfFile,additionalTag1,additionalTag2,additionalTag3);
      
      if(writeInputToFile){
        WriteSignalToFile("DecreasingSine_input",xVals,yVals,zVals);
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

      std::string nameTagForPerfFile = "---";
        if(InputSignals->Get("DecreasingSawtooth")->Has("NameTagPerFile")){
        InputSignals->Get("DecreasingSawtooth")->GetValue("NameTagPerFile",nameTagForPerfFile,ParamNode::PASS);
      }
      if(nameTagForPerfFile == "---"){
        nameTagForPerfFile = outputName;
      }

      CreatePeriodicTestSignal("DecreasingSawtooth",amplitudeScaling,numPeriods,stepsPerPeriod,xVals,yVals);
      zVals = xVals;
      zVals.Init();
      
      /*
       * START ACTUAL TESTING ROUTINE
       */
      TestHystOperatorWithSignal(outputName,xVals,yVals,zVals,dimHystOperatorForTesting,forcedPreisachResolution,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains,nameTagForPerfFile,additionalTag1,additionalTag2,additionalTag3);
      
      if(writeInputToFile){
        WriteSignalToFile("DecreasingSawtooth_input",xVals,yVals,zVals);
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

      std::string nameTagForPerfFile = "---";
        if(InputSignals->Get("Forc")->Has("NameTagPerFile")){
        InputSignals->Get("Forc")->GetValue("NameTagPerFile",nameTagForPerfFile,ParamNode::PASS);
      }
      if(nameTagForPerfFile == "---"){
        nameTagForPerfFile = outputName;
      }

      CreatePeriodicTestSignal("Forc",amplitudeScaling,numPeriods,stepsPerPeriod,xVals,yVals);
      zVals = xVals;
      zVals.Init();
      
      /*
       * START ACTUAL TESTING ROUTINE
       */
      TestHystOperatorWithSignal(outputName,xVals,yVals,zVals,dimHystOperatorForTesting,forcedPreisachResolution,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains,nameTagForPerfFile,additionalTag1,additionalTag2,additionalTag3);
      
      if(writeInputToFile){
        WriteSignalToFile("Forc_input",xVals,yVals,zVals);
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

      std::string nameTagForPerfFile = "---";
        if(InputSignals->Get("SelfDesigned")->Has("NameTagPerFile")){
        InputSignals->Get("SelfDesigned")->GetValue("NameTagPerFile",nameTagForPerfFile,ParamNode::PASS);
      }
      if(nameTagForPerfFile == "---"){
        nameTagForPerfFile = outputName;
      }

      CreateNonPeriodicTestSignal("SelfDesigned",amplitudeScaling,numberOfSteps,xVals,yVals);
      zVals = xVals;
      zVals.Init();
      
      /*
       * START ACTUAL TESTING ROUTINE
       */
      TestHystOperatorWithSignal(outputName,xVals,yVals,zVals,dimHystOperatorForTesting,forcedPreisachResolution,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains,nameTagForPerfFile,additionalTag1,additionalTag2,additionalTag3);
      
      if(writeInputToFile){
        WriteSignalToFile("SelfDesigned_input",xVals,yVals,zVals);
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

      std::string nameTagForPerfFile = "---";
        if(InputSignals->Get("SatX-RemX-SatY")->Has("NameTagPerFile")){
        InputSignals->Get("SatX-RemX-SatY")->GetValue("NameTagPerFile",nameTagForPerfFile,ParamNode::PASS);
      }
      if(nameTagForPerfFile == "---"){
        nameTagForPerfFile = outputName;
      }

      CreateNonPeriodicTestSignal("SatX-RemX-SatY",amplitudeScaling,numberOfSteps,xVals,yVals);
      zVals = xVals;
      zVals.Init();
      
      /*
       * START ACTUAL TESTING ROUTINE
       */
      TestHystOperatorWithSignal(outputName,xVals,yVals,zVals,dimHystOperatorForTesting,forcedPreisachResolution,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains,nameTagForPerfFile,additionalTag1,additionalTag2,additionalTag3);
      
      if(writeInputToFile){
        WriteSignalToFile("SatX-RemX-SatY_input",xVals,yVals,zVals);
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

      std::string nameTagForPerfFile = "---";
        if(InputSignals->Get("RemDrop")->Has("NameTagPerFile")){
        InputSignals->Get("RemDrop")->GetValue("NameTagPerFile",nameTagForPerfFile,ParamNode::PASS);
      }
      if(nameTagForPerfFile == "---"){
        nameTagForPerfFile = outputName;
      }

      CreateNonPeriodicTestSignal("RemDrop",amplitudeScaling,numberOfSteps,xVals,yVals);
      zVals = xVals;
      zVals.Init();
      
      /*
       * START ACTUAL TESTING ROUTINE
       */
      TestHystOperatorWithSignal(outputName,xVals,yVals,zVals,dimHystOperatorForTesting,forcedPreisachResolution,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains,nameTagForPerfFile,additionalTag1,additionalTag2,additionalTag3);
      
      if(writeInputToFile){
        WriteSignalToFile("RemDrop_input",xVals,yVals,zVals);
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
      std::string fileNameZ = "None";
      if(InputSignals->Get("ReadFromFile")->Has("FileContaining_z_over_steps")){
        InputSignals->Get("ReadFromFile")->GetValue("FileContaining_z_over_steps",fileNameZ,ParamNode::PASS);
      }
      bool test1D = false;
      if(InputSignals->Get("ReadFromFile")->Has("Test1D")){
        InputSignals->Get("ReadFromFile")->GetValue("Test1D",test1D,ParamNode::PASS);
      }

      Vector<Double> steps = Vector<Double>(1);
      
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

      if(fileNameY != "None"){
        reader.ReadFile(fileNameY.c_str(),steps,yVals);
        combinedname << "-";
        combinedname << fileNameY;
        if(yVals.GetSize() > maxSize){
          maxSize = yVals.GetSize();
        }
      }
      
      if(fileNameZ != "None"){
        reader.ReadFile(fileNameZ.c_str(),steps,zVals);
        combinedname << "-";
        combinedname << fileNameZ;
        if(dimHystOperatorForTesting < 3){
          std::cout << "3d input signal provided; scalar or 3d vector hysteresis operator will be used" << std::endl;
          dimHystOperatorForTesting = 3;
        }
        if(zVals.GetSize() > maxSize){
          maxSize = zVals.GetSize();
        }
      }

      std::string outputName = "---";
      if(InputSignals->Get("ReadFromFile")->Has("OutputName")){
        InputSignals->Get("ReadFromFile")->GetValue("OutputName",outputName,ParamNode::PASS);
      }
      if(outputName == "---"){
        outputName = combinedname.str();
      }

      std::string nameTagForPerfFile = "---";
        if(InputSignals->Get("ReadFromFile")->Has("NameTagPerFile")){
        InputSignals->Get("ReadFromFile")->GetValue("NameTagPerFile",nameTagForPerfFile,ParamNode::PASS);
      }
      if(nameTagForPerfFile == "---"){
        nameTagForPerfFile = outputName;
      }

      if(xVals.GetSize() < maxSize){
        xVals.Resize(maxSize);
      }
      if(yVals.GetSize() < maxSize){
        yVals.Resize(maxSize);
      }
      if(zVals.GetSize() < maxSize){
        zVals.Resize(maxSize);
      }

      xVals.ScalarMult(amplitudeScaling);
      yVals.ScalarMult(amplitudeScaling);
      zVals.ScalarMult(amplitudeScaling);
      
      /*
       * START ACTUAL TESTING ROUTINE
       */
      TestHystOperatorWithSignal(outputName,xVals,yVals,zVals,dimHystOperatorForTesting,forcedPreisachResolution,testInversion,printStatistics,writeResultsToFile,
              measurePerformance,commonPerfFile,test1D,outputIrrStrains,nameTagForPerfFile,additionalTag1,additionalTag2,additionalTag3);
      
      combinedname << "_input";
      if(writeInputToFile){
        WriteSignalToFile(combinedname.str(),xVals,yVals,zVals);
      }
    }

    if(stopAfterTests){
      std::cout << ("--- Stop after testing ---") << std::endl;
//      exit(0);
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
  }

  void CoefFunctionHyst::TestJacobianApproximations(){
		/*
     * 1. Create temporal hyst operator; this should be done for each test to ensure
     *		that operator is in initial state
     */
		bool isVirgin = true;
//		bool vector = true;

		if (POL_operatorParams_.methodName_ == "scalarPreisach") {

			POL_useExtension_ = false;
			if(POL_useExtension_){
				EXCEPTION("Extension not implemented for tests; remove completely as not working");
			} else {
//				vector = false;
        bool ignoreAnhystPart = false;
				hystTMP = new Preisach(1, POL_operatorParams_, POL_weightParams_, isVirgin, ignoreAnhystPart);

//        if(InversionParams_.inversionMethod != 4){
////          std::cout << "Inversion of scalar model with Newton/LM selected; will restrict input to operator direction \n"
////                  "as otherwise inversion will not work" << std::endl;
        // 8.6.2020: set parameter for inversion in order to pass user specified tolerances
          hystTMP->SetParamsForInversion(InversionParams_);
        }
//      }
		} else if (POL_operatorParams_.methodName_ == "vectorPreisach_Sutor") {
			if (POL_operatorParams_.evalVersion_ == 1) {
				POL_operatorParams_.isClassical_ = true; // original vector preisach model -> sutor2012

				hystTMP = new VectorPreisachSutor_ListApproach(1, POL_operatorParams_, POL_weightParams_, dim_, isVirgin);
			} else if (POL_operatorParams_.evalVersion_ == 2) {
				POL_operatorParams_.isClassical_ = false; // revised vector preisach model -> sutor2015

				hystTMP = new VectorPreisachSutor_ListApproach(1, POL_operatorParams_, POL_weightParams_, dim_, isVirgin);
			} else if (POL_operatorParams_.evalVersion_ == 10) {
				POL_operatorParams_.isClassical_ = true; // original vector preisach model -> sutor2015; matrix based implementation

				hystTMP = new VectorPreisachSutor_MatrixApproach(1, POL_operatorParams_, POL_weightParams_, dim_, isVirgin);
			} else if (POL_operatorParams_.evalVersion_ == 20) {
				POL_operatorParams_.isClassical_ = false; // revised vector preisach model -> sutor2015; matrix based implementation

				hystTMP = new VectorPreisachSutor_MatrixApproach(1, POL_operatorParams_, POL_weightParams_, dim_, isVirgin);
			} else {
				EXCEPTION("POL_operatorParams_.evalVersion_ has to be one of the following: \n "
                "1: classical vector model (sutor2012) \n"
                "2: revised vector model (sutor2015) [DEFAULT] \n"
                "10: classical vector model (sutor2012) - Matrix implementation, only for reference \n"
                "20: revised vector model (sutor2015) - Matrix implementation, only for reference \n")
			}
      hystTMP->SetParamsForInversion(InversionParams_);

		} else if (POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") {
      // basically a scalar model in multiple directions
      // isotropic case: all scalar models are equal (same weights etc)
      // anisotropic case: each model different; choice of directions matters; weights are harder to obtain
      int isIsotropic = 1;
      material_->GetScalar(isIsotropic, PREISACH_MAYERGOYZ_ISOTROPIC);
      if(isIsotropic == 0){
//        if( (dim_ != 2) || (isIsotropic == 0)){
        EXCEPTION("Mayergoyz vector model currently only implemented for 2d isotropic materials");
      }
      /*
       * IMPORTANT REMARK:
       *  > although the Mayergoyz model is based on the scalar models
       *     we are not allowed to directly apply the Preisach parameter for
       *     the scalar case (i.e. the weights, the anhyst parameter and so on)
       *  > make sure that the passed parameter are already transformed correctly
       *      > see constructor above
       */
      hystTMP = new VectorPreisachMayergoyz(1, POL_operatorParams_, POL_weightParams_, dim_, isVirgin);
      hystTMP->SetParamsForInversion(InversionParams_);

		} else {
			EXCEPTION("Invalid model selected for inversion test");
		}

    /*
     * 2. define some test vector
     */
    Vector<Double> xIn = Vector<Double>(dim_);
    Vector<Double> hIn = Vector<Double>(dim_);
    Vector<Double> yIn = Vector<Double>(dim_);
    xIn[0] = -1.0;
    xIn[1] = 9234;

    int dummyFlag = 0;
    hIn = hystTMP->computeValue_vec(xIn, 0, false, false, dummyFlag);

    yIn.Init();
    yIn.Add(1.0,hIn);
    for(UInt j = 0; j < dim_; j++){
      yIn[j] += eps_mu_base_[j][j]*xIn[j];
    }

    Matrix<Double> eps_mu_inv = Matrix<Double>(dim_,dim_);
    eps_mu_base_.Invert(eps_mu_inv);

    /*
     * 3. call actual test
     */
    hystTMP->CompareJacobianApproximations(xIn, yIn, eps_mu_inv, 0);
  }

  void CoefFunctionHyst::WriteTestStepToInfoXML(const UInt testNum,
          const Double Px, const Double Py, const Double Pz,
          const Double Pabs, const Double Pangle1, const Double Pangle2){

    // cp. EigenfrequencyDriver.cc
    PtrParamNode res;
    if(testNum <= 1){
      res = infoNodeLoc_->Get("result", ParamNode::APPEND);
    } else  {
      // create or take the last
      ParamNodeList l = infoNodeLoc_->GetList("result");
      res = l.IsEmpty() ? infoNodeLoc_->Get("result") : l.Last();
    }

    PtrParamNode testStep = res->Get("Teststep", ParamNode::APPEND);

    testStep->Get("nr")->SetValue(testNum); 
    testStep->Get("Px")->SetValue(Px,15);
    testStep->Get("Py")->SetValue(Py,15);
    testStep->Get("Pz")->SetValue(Pz,15);
    testStep->Get("Pabs")->SetValue(Pabs,15);
    testStep->Get("phi")->SetValue(Pangle1,15);
    testStep->Get("theta")->SetValue(Pangle2,15);
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
              << "  initial/smallfield tensor " << eps_mu_base_.ToString() << "\n"
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

  Hysteresis* CoefFunctionHyst::getTemporalHystOperator(bool forStrains,bool forceScalarDirection, UInt dimOfHystOperator, int forcedPreisachResolutionForTests, bool enforceContinuousEvaluation){
      /*
       * Create temporal hyst operator; e.g. for testing (-> TestHystOperator) or
       * for estimating the maximal and minimal slopes
       *
       * > Create this hyst operator with the same parameter set as the actually used operators in the system
       */

      Hysteresis* tmpHystOperator = NULL;
      bool isVirgin = true;

      // enable warnings
      POL_operatorParams_.printWarnings_ = true;
      STRAIN_operatorParams_.printWarnings_ = true;

      std::cout << "2. - getTemporalHystOperator" << std::endl;
    std::cout << "Check if new parameter inputForAlignment_ has been set" << std::endl;
    std::cout << "POL_operatorParams.inputForAlignment_ = " << POL_operatorParams_.inputForAlignment_ << std::endl;
      
      ParameterPreisachWeights weightParamsForTMPOperator;
      if(forcedPreisachResolutionForTests > 0){
        // get temporal weights with different resolution than specified in mat.xml
        ReadAndSetWeights(material_, forStrains, forcedPreisachResolutionForTests);
        if(forStrains){
          weightParamsForTMPOperator = STRAIN_weightParamsForTesting_;
        } else {
          weightParamsForTMPOperator = POL_weightParamsForTesting_;
        }
      } else {
        if(forStrains){
          weightParamsForTMPOperator = STRAIN_weightParams_;
        } else {
          weightParamsForTMPOperator = POL_weightParams_;
        }
      }
      
      if(!forStrains){

        if (POL_operatorParams_.methodName_ == "scalarPreisach") {
          POL_useExtension_ = false;

          if(POL_useExtension_){
            EXCEPTION("Extension not implemented for tests; remove completely as not working");
          } else {
            tmpHystOperator = new Preisach(1, POL_operatorParams_, weightParamsForTMPOperator, isVirgin, false);

            // important: for testing/tracing the scalarPreisach operator, we have to ensure, that the internal
            // direction is along the x-axis as this is the tracing direction!
            if(forceScalarDirection){
              Vector<Double> xDir = Vector<Double>(dimOfHystOperator);
              xDir.Init();
              xDir[0] = 1.0;
              tmpHystOperator->setFixDirection(xDir);
            }

//            if(InversionParams_.inversionMethod != 4){
//              //          std::cout << "Inversion of scalar model with Newton/LM selected; will restrict input to operator direction \n"
//              //                  "as otherwise inversion will not work" << std::endl;
            // 8.6.2020: set inversion params to pass user specified tolerances 
              tmpHystOperator->SetParamsForInversion(InversionParams_);
//            }
          }
        } else if (POL_operatorParams_.methodName_ == "vectorPreisach_Sutor") {
          if (POL_operatorParams_.evalVersion_ == 1) {
            POL_operatorParams_.isClassical_ = true; // original vector preisach model -> sutor2012

            tmpHystOperator = new VectorPreisachSutor_ListApproach(1, POL_operatorParams_, weightParamsForTMPOperator, dimOfHystOperator, isVirgin);

          } else if (POL_operatorParams_.evalVersion_ == 2) {
            POL_operatorParams_.isClassical_ = false; // revised vector preisach model -> sutor2015

            tmpHystOperator = new VectorPreisachSutor_ListApproach(1, POL_operatorParams_, weightParamsForTMPOperator, dimOfHystOperator, isVirgin);

          } else if (POL_operatorParams_.evalVersion_ == 10) {
            POL_operatorParams_.isClassical_ = true; // original vector preisach model -> sutor2015; matrix based implementation

            if(enforceContinuousEvaluation){
              std::cout << "Info: Using list-based implementation for tracing instead of matrix-based version to avoid uncrealistically high slopes resulting from the discrete steppings." << std::endl;
              tmpHystOperator = new VectorPreisachSutor_ListApproach(1, POL_operatorParams_, weightParamsForTMPOperator, dimOfHystOperator, isVirgin);     
            } else {
              WARN("Tracing the discrete matrix-based version of the vector model may cause unusable parameter for local fixpoint inversion. Please trace list-based implementation instead.");
              tmpHystOperator = new VectorPreisachSutor_MatrixApproach(1, POL_operatorParams_, weightParamsForTMPOperator, dimOfHystOperator, isVirgin);
            }
          } else if (POL_operatorParams_.evalVersion_ == 20) {
            POL_operatorParams_.isClassical_ = false; // revised vector preisach model -> sutor2015; matrix based implementation

            if(enforceContinuousEvaluation){
              std::cout << "Info: Using list-based implementation for tracing instead of matrix-based version to avoid uncrealistically high slopes resulting from the discrete steppings." << std::endl;
              tmpHystOperator = new VectorPreisachSutor_ListApproach(1, POL_operatorParams_, weightParamsForTMPOperator, dimOfHystOperator, isVirgin);
            } else {
              WARN("Tracing the discrete matrix-based version of the vector model may cause unusable parameter for local fixpoint inversion. Please trace list-based implementation instead.");
              tmpHystOperator = new VectorPreisachSutor_MatrixApproach(1, POL_operatorParams_, weightParamsForTMPOperator, dimOfHystOperator, isVirgin);
            }
          } else {
            EXCEPTION("POL_operatorParams_.evalVersion_ has to be one of the following: \n "
              "1: classical vector model (sutor2012) \n"
              "2: revised vector model (sutor2015) [DEFAULT] \n"
              "10: classical vector model (sutor2012) - Matrix implementation, only for reference \n"
              "20: revised vector model (sutor2015) - Matrix implementation, only for reference \n")
          }
          tmpHystOperator->SetParamsForInversion(InversionParams_);

        } else if (POL_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") {
          // basically a scalar model in multiple directions
          // isotropic case: all scalar models are equal (same weights etc)
          // anisotropic case: each model different; choice of directions matters; weights are harder to obtain
          int isIsotropic = 1;
          material_->GetScalar(isIsotropic, PREISACH_MAYERGOYZ_ISOTROPIC);
          if(isIsotropic == 0){
//            if( (dim_ != 2) || (isIsotropic == 0)){
            EXCEPTION("Mayergoyz vector model currently only implemented for isotropic materials");
          }

          /*
           * IMPORTANT REMARK:
           *  > although the Mayergoyz model is based on the scalar models
           *     we are not allowed to directly apply the Preisach parameter for
           *     the scalar case (i.e. the weights, the anhyst parameter and so on)
           *  > make sure that the passed parameter are already transformed correctly
           *      > see constructor above
           */
          tmpHystOperator = new VectorPreisachMayergoyz(1, POL_operatorParams_, weightParamsForTMPOperator, dimOfHystOperator, isVirgin);

          tmpHystOperator->SetParamsForInversion(InversionParams_);

        } else {
          EXCEPTION("Invalid model selected for inversion test");
        }
      } else {
        /*
         * these operators still compute the polarization (not the strains directly) but different parameters could be
         * set here; furthermore we do not set parameter for inversion
         */
        if (STRAIN_operatorParams_.methodName_ == "scalarPreisach") {
          bool ignoreAnhystPart = false;
          tmpHystOperator = new Preisach(1, STRAIN_operatorParams_, weightParamsForTMPOperator, isVirgin, ignoreAnhystPart);

          // important: for testing/tracing the scalarPreisach operator, we have to ensure, that the internal
          // direction is along the x-axis as this is the tracing direction!
          if(forceScalarDirection){
            Vector<Double>xDir = Vector<Double>(dimOfHystOperator);
            xDir.Init();
            xDir[0] = 1.0;
            tmpHystOperator->setFixDirection(xDir);
          }

        } else if (STRAIN_operatorParams_.methodName_ == "vectorPreisach_Sutor") {
          if (STRAIN_operatorParams_.evalVersion_ == 1) {
            STRAIN_operatorParams_.isClassical_ = true; // original vector preisach model -> sutor2012

            tmpHystOperator = new VectorPreisachSutor_ListApproach(1, STRAIN_operatorParams_, weightParamsForTMPOperator, dimOfHystOperator, isVirgin);

          } else if (STRAIN_operatorParams_.evalVersion_ == 2) {
            STRAIN_operatorParams_.isClassical_ = false; // revised vector preisach model -> sutor2015

            tmpHystOperator = new VectorPreisachSutor_ListApproach(1, STRAIN_operatorParams_, weightParamsForTMPOperator, dimOfHystOperator, isVirgin);

          } else if (STRAIN_operatorParams_.evalVersion_ == 10) {
            STRAIN_operatorParams_.isClassical_ = true; // original vector preisach model -> sutor2015; matrix based implementation
            if(enforceContinuousEvaluation){
              std::cout << "Info: Using list-based implementation for tracing instead of matrix-based version to avoid uncrealistically high slopes resulting from the discrete steppings." << std::endl;
              tmpHystOperator = new VectorPreisachSutor_ListApproach(1, STRAIN_operatorParams_, weightParamsForTMPOperator, dimOfHystOperator, isVirgin);     
            } else {
              WARN("Tracing the discrete matrix-based version of the vector model may cause unusable parameter for local fixpoint inversion. Please trace list-based implementation instead.");
              tmpHystOperator = new VectorPreisachSutor_MatrixApproach(1, STRAIN_operatorParams_, weightParamsForTMPOperator, dimOfHystOperator, isVirgin);
            }
//            tmpHystOperator = new VectorPreisachSutor_MatrixApproach(1, STRAIN_operatorParams_, weightParamsForTMPOperator, dimOfHystOperator, isVirgin);

          } else if (STRAIN_operatorParams_.evalVersion_ == 20) {
            STRAIN_operatorParams_.isClassical_ = false; // revised vector preisach model -> sutor2015; matrix based implementation

            if(enforceContinuousEvaluation){
              std::cout << "Info: Using list-based implementation for tracing instead of matrix-based version to avoid uncrealistically high slopes resulting from the discrete steppings." << std::endl;
              tmpHystOperator = new VectorPreisachSutor_ListApproach(1, STRAIN_operatorParams_, weightParamsForTMPOperator, dimOfHystOperator, isVirgin);
            } else {
              WARN("Tracing the discrete matrix-based version of the vector model may cause unusable parameter for local fixpoint inversion. Please trace list-based implementation instead.");
              tmpHystOperator = new VectorPreisachSutor_MatrixApproach(1, STRAIN_operatorParams_, weightParamsForTMPOperator, dimOfHystOperator, isVirgin);
            }
            //tmpHystOperator = new VectorPreisachSutor_MatrixApproach(1, STRAIN_operatorParams_, weightParamsForTMPOperator, dimOfHystOperator, isVirgin);

          }
        } else if (STRAIN_operatorParams_.methodName_ == "vectorPreisach_Mayergoyz") {

          tmpHystOperator = new VectorPreisachMayergoyz(1, STRAIN_operatorParams_, weightParamsForTMPOperator, dimOfHystOperator, isVirgin);

        }
      }

      // disable warnings again
      POL_operatorParams_.printWarnings_ = false;
      STRAIN_operatorParams_.printWarnings_ = false;

      if(tmpHystOperator == NULL){
        EXCEPTION("temporary hysteresis operator is NULL!")
      }
      
      return tmpHystOperator;
    }

  void CoefFunctionHyst::TraceHystOperator(UInt baseSteps, Double& maxSlope, Double& minSlope, Double& negCoercivity, Double& maxPolarization, bool dedicatedOperatorForStrains){
    EXCEPTION("TraceHystOperator (non-vector-version) Should not be needed/used anymore");
//
//    /*
//       * Get temporal hyst operator first
//       * > has to have same parameter as the actually used operators but must be temporal to not mess up with storage
//       */
//      bool forceScalarDirection = false;
//      if(dedicatedOperatorForStrains){
//        if(STRAIN_operatorParams_.methodName_ == "scalarPreisach"){
//          std::cout << "++ Tracing Dedicated SCALAR Strain Hysteresis Operator for Material " << material_->GetName() << std::endl;
//          forceScalarDirection = true;
//        } else {
//          std::cout << "++ Tracing Dedicated VECTOR Strain Hysteresis Operator for Material " << material_->GetName() << std::endl;
//          forceScalarDirection = false;
//        }
//      } else {
//        if(POL_operatorParams_.methodName_ == "scalarPreisach"){
//          std::cout << "++ Tracing Dedicated SCALAR Hysteresis Operator for Material " << material_->GetName() << std::endl;
//          forceScalarDirection = true;
//        } else {
//          std::cout << "++ Tracing Dedicated VECTOR Hysteresis Operator for Material " << material_->GetName() << std::endl;
//          forceScalarDirection = false;
//        }
//      }
//
//      Hysteresis* hystOperatorForTrace = getTemporalHystOperator(dedicatedOperatorForStrains,forceScalarDirection,dim_);
//      Timer* traceTimer = new Timer();
//      Double startTime = traceTimer->GetCPUTime();
//      traceTimer->Start();
//
//      /*
//       * Define signal to trace the relevant parts of the hyst loop
//       * 1) Virgin Curve to Saturation
//       * 2) Saturation to multiples of Saturation (for the case that Hyst model grows beyound hysteresis)
//       * 3) Saturation to remanence
//       * 4) Remanence to negative saturation
//       */
//      UInt virginSteps = 2*baseSteps;
//      UInt beyondSaturationSteps = baseSteps;
//      UInt saturationToRemanenceSteps = 2*baseSteps;
//      UInt remanenceToNegSaturationSteps = 5*baseSteps;
//      UInt totalSteps = virginSteps + beyondSaturationSteps + saturationToRemanenceSteps + remanenceToNegSaturationSteps;
//
//      Vector<Double> scalarInputs = Vector<Double>(totalSteps);
//      Vector<Double> scalarOutputs = Vector<Double>(totalSteps);
//      Vector<Double> estimatedSlope = Vector<Double>(totalSteps);
//
//      Double inputSat, outputSat;
//      if(dedicatedOperatorForStrains){
//        inputSat = STRAIN_operatorParams_.inputSat_;
//        outputSat = STRAIN_operatorParams_.outputSat_;
//      } else {
//        inputSat = POL_operatorParams_.inputSat_;
//        outputSat = POL_operatorParams_.outputSat_;
//      }
//
//
//      Vector<UInt> steps =  Vector<UInt>(4);
//      steps[0] = virginSteps;
//      steps[1] = beyondSaturationSteps;
//      steps[2] = saturationToRemanenceSteps;
//      steps[3] = remanenceToNegSaturationSteps;
//
//      Vector<UInt> positionOffsets = Vector<UInt>(4);
//      positionOffsets[0] = 0;
//      positionOffsets[1] = steps[0];
//      positionOffsets[2] = steps[0] + steps[1];
//      positionOffsets[3] = steps[0] + steps[1] + steps[2];
//
//      Double maxAmplitude = 10*inputSat;
//      Vector<Double> endAmplitudes = Vector<Double>(6);
//      endAmplitudes[0] = inputSat; // go to saturation
//      endAmplitudes[1] = maxAmplitude; // go beyond saturation
//      endAmplitudes[2] = 0; // go to remanence
//      endAmplitudes[3] = -inputSat; // go to negative saturation
//
//      for(UInt k = 0; k < 4; k++){
//        Double increment = 0.0;
//        Double startAmplitude = 0.0;
//        if(k != 0){
//          startAmplitude = endAmplitudes[k-1];
//        }
//        increment = (endAmplitudes[k]-startAmplitude)/(Double (steps[k]));
//        for(UInt i = 0; i < steps[k]; i++){
//          Double currentAmplitude = startAmplitude + (i+1)*increment; // we want to reach amplitudeEnd at maximal i
//          scalarInputs[positionOffsets[k] + i] = currentAmplitude;
//        }
//      }
//
//      /*
//       * Trace hyst operator (only in 1d direction; should be ok for estimate of slopes)
//       * > true for most cases, but for revised sutor model, we get problems!
//       * > here we have the angular distance which can completely wreck up slope!
//       * > additionally check for remdrop, i.e. at end when material is in remanence, drive it to saturation by
//       *   increasing field perpendicular to polarization
//       */
//      int successFlagForward = 0;
//      maxSlope = 0;
//      minSlope = 1e16;
//      // get some additinal info if we are already traverse loop
//      negCoercivity = 0.0;
//      maxPolarization = 0.0;
//
//      UInt idxMaxSlope, idxMinSlope;
//      Double fieldMaxPolarization;
//      Double fieldMaxSlope, fieldMinSlope;
//
//      Vector<Double> vecIn = Vector<Double>(dim_);
//      Vector<Double> vecOut = Vector<Double>(dim_);
//
//      for(UInt i = 0; i < totalSteps; i++){
//        vecIn.Init();
//        vecIn[0] = scalarInputs[i];
//
//        vecOut = hystOperatorForTrace->computeValue_vec(vecIn, 0, true, false, successFlagForward);
//        scalarOutputs[i] = vecOut[0];
//
//        if(i > 0){
//          estimatedSlope[i] = (scalarOutputs[i] - scalarOutputs[i-1])/(scalarInputs[i] - scalarInputs[i-1]);
//
//          if(estimatedSlope[i] < 0){
//            if(abs(estimatedSlope[i]) > 1e-16){
//              std::stringstream warnmsg;
//              warnmsg << "Negative slope (" << estimatedSlope[i]<<") detected between " << scalarInputs[i] << " and " << scalarInputs[i-1] << std::endl;
//              warnmsg << "Hysteresis model (Preisach) is assumed to be monoton (at least in the 1d case!)" << std::endl;
//              WARN(warnmsg.str());
//            }
//          }
//          if(estimatedSlope[i] > maxSlope){
//            maxSlope = estimatedSlope[i];
//            fieldMaxSlope = scalarInputs[i];
//            idxMaxSlope = i;
//          }
//          if(estimatedSlope[i] < minSlope){
//            minSlope = estimatedSlope[i];
//            fieldMinSlope = scalarInputs[i];
//            idxMinSlope = i;
//          }
//
//          if(scalarOutputs[i] > maxPolarization){
//            // only positive values checked here as we drive system beyond the positive polarization only
//            maxPolarization = scalarOutputs[i];
//            fieldMaxPolarization = scalarInputs[i];
//          }
//
//          if(scalarOutputs[i] == 0){
//            // excact hit; very unlikely
//            negCoercivity = scalarInputs[i];
//          }
//          if((scalarOutputs[i] < 0) && (scalarOutputs[i-1] > 0)){
//            // negative coercivity somewhere between the two inputs
//            // > note: due to the form of the tested signal, we just hit the negative coercivity here and only have to check for it
//            negCoercivity = (scalarInputs[i] + scalarInputs[i-1])/2.0;
//          }
//        }
//      }
//
//      traceTimer->Stop();
//      Double endTime = traceTimer->GetCPUTime();
//
//      bool testOutput = true;
//      if(testOutput){
////        std::cout << "Traced hyst operator for material " << material_->GetName() << std::endl;
//
//        std::ofstream tracedData;
//        std::stringstream tracedData_FILENAME;
//        tracedData_FILENAME << "Traced_Hysteresis_data_" << material_->GetName() << ".trace";
//
//        tracedData.open(tracedData_FILENAME.str());
//        tracedData << "### INFO " << std::endl;
//        tracedData << "# > hysteresis operator was tested with 1d input signal " << std::endl;
//        tracedData << "# > maximal slope  " << maxSlope << "\t (found at (idx,Field) = (" << idxMaxSlope << "," << fieldMaxSlope << ") )" << std::endl;
//        tracedData << "# > minimal slope  " << minSlope << "\t (found at (idx,Field) = (" << idxMinSlope << "," << fieldMinSlope << ") )" << std::endl;
//        tracedData << "# > maximal polarization  " << maxPolarization << "\t (found for Field " << fieldMaxPolarization << ")" << std::endl;
//        tracedData << "# > specified saturatin values for field and polarization  " << inputSat << ", " << outputSat << std::endl;
//        tracedData << "# > negative coercivity found for Field " << negCoercivity << std::endl;
//        tracedData << "# > number of traced positions " << totalSteps << std::endl;
//        tracedData << "# > requied runtime " << endTime-startTime << std::endl;
//        tracedData << "#" << std::endl;
//        tracedData << "# Traced data: " << std::endl;
//        tracedData << "#IDX   Field   Polarization" << std::endl;
//        for(UInt i = 0; i < totalSteps; i++){
//          tracedData << i << "\t" << scalarInputs[i] << "\t" << scalarOutputs[i] << std::endl;
//        }
//        tracedData.close();
//      }
//
//      /*
//       * free memory
//       */
//      delete hystOperatorForTrace;
//      delete traceTimer;
    }

  /*
   * Unfortunately, for revised sutor model, we cannot just rely on the slope along the excitation axis but must also
   * consider the slopes along the other axis, i.e. we have to calculate/estimate the Jacobian at each point along the
   * tracing curve
   * TODO: implement me!
   */
  void CoefFunctionHyst::TraceHystOperatorVector(UInt baseSteps, Double& maxSlope, Double& minSlope, Double& negCoercivity,
                                                          Double& maxPolarization, bool dedicatedOperatorForStrains){
//    std::cout << "++ Tracing Vector Hysteresis Operator for material " << material_->GetName() << std::endl;
    /*
       * Get temporal hyst operator first
       * > has to have same parameter as the actually used operators but must be temporal to not mess up with storage
       */
      bool startAtRemanence = true;
      UInt precisionDigits = 12;//+1 for output
      
      // Note 6.11.2020:
      // 1. the difference between centralFD and forwardFD is negligble for
      // most cases (only if jumps occur central differences show other results);
      // > no additional scaling factor required! 
      // > due to this additional factor, the testsuite examples were solved differently
      // 2. the scaling factor has significant influence on the derived Jacobian
      // > due to chaning it from 1e-5 to 1e-7 the max slope for the electrostatic testcases
      // changed by a factor of 4 which in turn leads to different convergence behavior
      /*
       * make scaling, additonalSafetyFactor and forceCentral to intput
       * arguments that are read from simulation xml file
       * > get read in via solve step and then passed to this function
       * > cmp. contraction factors for FP iteration 
       */

      // minimal stepping for FiniteDifference stepping
      Double scaling = trace_JacResolution_; //1e-5
      // edit 5.11.2020: 
      // central differences result in less peaks in the Jacobian which in turn leads to better convergence factors
      // for the fixedpoint methods (local and global); unfortunately, the computed slopes are sometimes too small
      // to allow for global convergence (as not all possible worst cases can be traced)
      // thus if we enforce central differences, we should increase the sefetyfactor a bit
      bool forceCentral = trace_forceCentral_;

//      if(forceCentral){
//        additionalSafetyFactor_ = 1.125;
//      } else {
//        additionalSafetyFactor_ = 1.00;
//      }
      /*
       * Note 26.6.2020 - there are some testcases that fail if traced data is loaded from file instead
       * of being traced on the fly. to ensure that this is not due different precisions, the directly traced and 
       * the loaded data have to be rounded/restricted to the same precision! 
       */
      
      bool forceScalarDirection = false;
      if(dedicatedOperatorForStrains){
        if(STRAIN_operatorParams_.methodName_ == "scalarPreisach"){
//          std::cout << "++ Tracing additional SCALAR Hysteresis Operator for Strains " << std::endl;
          forceScalarDirection = true;
        } else {
//          std::cout << "++ Tracing additional VECTOR Hysteresis Operator for Strains " << std::endl;
          forceScalarDirection = false;
        }
      } else {
        if(POL_operatorParams_.methodName_ == "scalarPreisach"){
//          std::cout << "++ Tracing SCALAR Hysteresis Operator for Material " << material_->GetName() << std::endl;
          forceScalarDirection = true;
        } else {
//          std::cout << "++ Tracing VECTOR Hysteresis Operator for Material " << material_->GetName() << std::endl;
          forceScalarDirection = false;
        }
      }

      std::stringstream tracedData_FILENAME;
      tracedData_FILENAME << material_->GetName();
      if(forceScalarDirection){
        tracedData_FILENAME << "_Scalar_HystOperator";
      } else {
        tracedData_FILENAME << "_Vector_HystOperator";
      }

      if(dedicatedOperatorForStrains){
        tracedData_FILENAME << "_forStrain";
      } else {
        tracedData_FILENAME << "_forPolarization";
      }
      if(startAtRemanence){
        tracedData_FILENAME << "_startWithRemdrop.trace";
      } else {
        tracedData_FILENAME << "_startAtZero.trace";
      }

      /*
       * check if this particular material has already been traced
       */
      TracedData currentMaterial;

      // search in interal array first > even if forceRetracing is true, we just trace every material just once
      // Use mutex for thread-safe access to static map during parallel assembly
      {
        std::lock_guard<std::mutex> lock(tracedDataMutex_);
        auto it = CoefFunctionHyst::tracedOperatorData_.find(tracedData_FILENAME.str());
        if(it != CoefFunctionHyst::tracedOperatorData_.end()){
          std::cout << "++ Material " << material_->GetName() << " has already been traced. Reuse values from internal array." << std::endl;
          currentMaterial = it->second;

          maxSlope = currentMaterial.maxSlope_;
          minSlope = currentMaterial.minSlope_;
          negCoercivity = currentMaterial.negCoercivity_;
          maxPolarization = currentMaterial.maxPolarization_;

          return;
        }
      }
      
      std::fstream tracedData;
      if(trace_forceRetracing_){
        std::cout << "++ Retracing of Material " << material_->GetName() << " has already been forced. Previous tracing-data will not be read from files." << std::endl;
      } else {
        // check if file with data already exists
        tracedData.open(tracedData_FILENAME.str(),std::ios_base::in);
        if(tracedData.is_open()){
          // material has already been traced; read data from file
          UInt numParamsToBeRead = 4;
          UInt numParamsRead = 0;
          size_t found;
          double tmp;
          std::string readLine;
          std::string dataField;
          std::string searchString;
          while(numParamsRead < numParamsToBeRead){
            if(std::getline(tracedData,readLine)){

              searchString = "###MaxSlope###";
              found = readLine.find(searchString); 
              if (found != string::npos){
                dataField = readLine.substr(found+searchString.length()+1);
                tmp = std::stod(dataField);
                maxSlope = tmp;
                numParamsRead++;
              } 

              searchString = "###MinSlope###";
              found = readLine.find(searchString); 
              if (found != string::npos){
                dataField = readLine.substr(found+searchString.length()+1);
                tmp = std::stod(dataField);
                minSlope = tmp;
                numParamsRead++;
              } 

              searchString = "###NegCoercivity###";
              found = readLine.find(searchString);
              if (found != string::npos){
                dataField = readLine.substr(found+searchString.length()+1);
                tmp = std::stod(dataField);
                negCoercivity = tmp;
                numParamsRead++;
              } 

              searchString = "###MaxPolarization###";
              found = readLine.find(searchString);
              if (found != string::npos){
                dataField = readLine.substr(found+searchString.length()+1);
                tmp = std::stod(dataField);
                maxPolarization = tmp;
                numParamsRead++;
              }           
            } else {
              break;
            }
          }
          tracedData.close();

          if(numParamsRead == numParamsToBeRead){
            /*
             * round data to ensure consistent results whether file was read or tracing was done on the fly
             */
  //          std::cout << "++ Material " << material_->GetName() << " has already been traced. All required parameter read from file." << std::endl;
  //          std::cout << "PRE-cutting" << std::endl;
  //          std::cout << std::scientific<< std::setprecision(precisionDigits+1) << "###MaxSlope### " << maxSlope << std::endl;
  //          std::cout << std::scientific<< std::setprecision(precisionDigits+1) << "###MinSlope### " << minSlope << std::endl;
  //          std::cout << std::scientific<< std::setprecision(precisionDigits+1) << "###NegCoercivity### " << negCoercivity << std::endl;
  //          std::cout << std::scientific<< std::setprecision(precisionDigits+1) << "###MaxPolarization### " << maxPolarization << std::endl;
  //          
            maxSlope = restrictPrecision(maxSlope,precisionDigits);
            minSlope = restrictPrecision(minSlope,precisionDigits);
            negCoercivity = restrictPrecision(negCoercivity,precisionDigits);
            maxPolarization = restrictPrecision(maxPolarization,precisionDigits);

            std::cout << "++ Material " << material_->GetName() << " has already been traced. All required parameter read from file." << std::endl;
            std::cout << std::scientific<< std::setprecision(precisionDigits+1) << "###MaxSlope### " << maxSlope << std::endl;
            std::cout << std::scientific<< std::setprecision(precisionDigits+1) << "###MinSlope### " << minSlope << std::endl;
            std::cout << std::scientific<< std::setprecision(precisionDigits+1) << "###NegCoercivity### " << negCoercivity << std::endl;
            std::cout << std::scientific<< std::setprecision(precisionDigits+1) << "###MaxPolarization### " << maxPolarization << std::endl;         

            // store for later usage with mutex for thread-safe access
            currentMaterial.maxSlope_ = maxSlope;
            currentMaterial.minSlope_ = minSlope;
            currentMaterial.negCoercivity_ = negCoercivity;
            currentMaterial.maxPolarization_ = maxPolarization;
            {
              std::lock_guard<std::mutex> lock(tracedDataMutex_);
              CoefFunctionHyst::tracedOperatorData_[tracedData_FILENAME.str()] = currentMaterial;
            }

            return;
          }
        }
      } 

      std::cout << "++ Tracing hysteresis operator for material " << material_->GetName() << std::endl;
      
      // if we trace the hyst operator with a discrete implementation, i.e., a matrix-based version of the hysteresis
      // models (should never be the case in practice due to inefficiency but for testing it is), the computed slopes
      // will be way larger than for the continuous evaluation using Everett-based functioning; this hugely impacts
      // the convergence behaviour for global and localized fixpoint methods and may easily result in non-convergence
      // whereas convergence would easily be doable in the continous case; therefore we force the tracing with continous
      // models if and hope that the retrieved parameter will work for the discrete models, too
      bool enforceContinuousEvalulation = true;
      int enforcedPreisachResolution = -1;
      Hysteresis* hystOperatorForTrace = getTemporalHystOperator(dedicatedOperatorForStrains,forceScalarDirection,dim_,enforcedPreisachResolution,enforceContinuousEvalulation);

//    /*
//       * Get temporal hyst operator first
//       * > has to have same parameter as the actually used operators but must be temporal to not mess up with storage
//       */
//      bool forceScalarDirection = true;
//      Hysteresis* hystOperatorForTrace = getTemporalHystOperator(dedicatedOperatorForStrains,forceScalarDirection);

      Timer* traceTimer = new Timer();
      Double startTime = traceTimer->GetCPUTime();
      traceTimer->Start();
      //std::cout << "++ Start " << std::endl;
      /*
       * Define signal to trace the relevant parts of the hyst loop
       * 1a) Virgin Curve to Saturation
       * 1b) Remanence in y; go to Saturation in x
       * 2) Saturation to multiples of Saturation (for the case that Hyst model grows beyound hysteresis)
       * 3) Saturation to remanence
       * 4) Remanence to negative saturation and beyond
       * 5) Back to positive saturation // NEW
       */
      UInt virginSteps = 2*baseSteps;
      UInt beyondSaturationSteps = baseSteps;
      UInt saturationToRemanenceSteps = 2*baseSteps;
      UInt remanenceToNegSaturationSteps = 5*baseSteps;
      UInt backToPosSaturationSteps = 1*baseSteps;
      UInt backToHalfNegSaturation = 1*baseSteps;
      UInt backToHalfPosSaturation = 1*baseSteps;
      UInt totalSteps = virginSteps + beyondSaturationSteps + saturationToRemanenceSteps + remanenceToNegSaturationSteps + backToPosSaturationSteps + backToHalfNegSaturation + backToHalfPosSaturation;

      Vector<Double> scalarInputs = Vector<Double>(totalSteps);
      Vector<Double> estimatedSlope = Vector<Double>(totalSteps);
      Vector<Double>* vectorOutputs = new Vector<Double>[totalSteps];
      Matrix<Double>* Jacobians = new Matrix<Double>[totalSteps];

      Double inputSat, outputSat;
      if(dedicatedOperatorForStrains){
        inputSat = STRAIN_operatorParams_.inputSat_;
        outputSat = STRAIN_operatorParams_.outputSat_;
      } else {
        inputSat = POL_operatorParams_.inputSat_;
        outputSat = POL_operatorParams_.outputSat_;
      }


      Vector<UInt> steps =  Vector<UInt>(7);
      steps[0] = virginSteps;
      steps[1] = beyondSaturationSteps;
      steps[2] = saturationToRemanenceSteps;
      steps[3] = remanenceToNegSaturationSteps;
      steps[4] = backToPosSaturationSteps;
      steps[5] = backToHalfNegSaturation;
      steps[6] = backToHalfPosSaturation;

      Vector<UInt> positionOffsets = Vector<UInt>(7);
      positionOffsets[0] = 0;
      positionOffsets[1] = steps[0];
      positionOffsets[2] = steps[0] + steps[1];
      positionOffsets[3] = steps[0] + steps[1] + steps[2];
      positionOffsets[4] = steps[0] + steps[1] + steps[2] + steps[3];
      positionOffsets[5] = steps[0] + steps[1] + steps[2] + steps[3] + steps[4];
      positionOffsets[6] = steps[0] + steps[1] + steps[2] + steps[3] + steps[4] + steps[5];

      Double maxAmplitude = 10*inputSat;
      Vector<Double> endAmplitudes = Vector<Double>(7);
      endAmplitudes[0] = inputSat; // go to saturation
      endAmplitudes[1] = maxAmplitude; // go beyond saturation
      endAmplitudes[2] = 0; // go to remanence
      endAmplitudes[3] = -maxAmplitude; // go beyond negative saturation
      endAmplitudes[4] = inputSat; // go back to pos saturation
      endAmplitudes[5] = -inputSat/2; // go back to half of neg saturation
      endAmplitudes[6] = inputSat/2; // go back to half of neg saturation

      for(UInt k = 0; k < 7; k++){
        Double increment = 0.0;
        Double startAmplitude = 0.0;
        if(k != 0){
          startAmplitude = endAmplitudes[k-1];
        }
        increment = (endAmplitudes[k]-startAmplitude)/(Double (steps[k]));
        for(UInt i = 0; i < steps[k]; i++){
          Double currentAmplitude = startAmplitude + (i+1)*increment; // we want to reach amplitudeEnd at maximal i
          scalarInputs[positionOffsets[k] + i] = currentAmplitude;
        }
      }

      /*
       * Trace hyst operator (only in 1d direction; should be ok for estimate of slopes)
       * > true for most cases, but for revised sutor model, we get problems!
       * > here we have the angular distance which can completely wreck up slope!
       * > additionally check for remdrop, i.e. at end when material is in remanence, drive it to saturation by
       *   increasing field perpendicular to polarization
       */
      int successFlagForward = 0;
      maxSlope = 0;
      minSlope = 1e16;
      // get some additinal info if we are already traverse loop
      negCoercivity = 0.0;
      maxPolarization = 0.0;

      UInt idxMaxSlope = 0;
      UInt idxMinSlope = 0;
      Double fieldMaxPolarization = 0;
      Double fieldMaxSlope = 0;
      Double fieldMinSlope = 0;

      Vector<Double> vecIn = Vector<Double>(dim_);
      Vector<Double> vecOut = Vector<Double>(dim_);

      Matrix<Double> tmp = Matrix<Double>(dim_,dim_);
      Vector<Double> Pshifted = Vector<Double>(dim_);
      Vector<Double> stepping = Vector<Double>(dim_);
      // for central differences; should be used if Jacobian is to be evaluated around the 0-point; here Jacobian
      // matrix might have extremely large entries compared to standard inputs!
      // > see note from VectorPreisachSutor:
      /*
       * Note 11-07-2019
       * the problem with 0 input only occurs if the the Preisach weights are unsymmetric and only as a consequence
       * of changing directions around 0-inputs
       * this occurs e.g., in the computation of the Jacobian around the point 0/0/0
       * here each evaluation direciton leads to a different rotation vector 1e-5/0/0 > 1/0/0; 0/1e-5/0 > 0/1/0 etc.
       * normally this would be no problem at any other point which is not (close to) 0/0/0 as we would have a valid
       * direction vector;
       * the suggested workaround thus is the following: whenever the Jacobian is to be evaluated (which is usually done
       * using forward or backward differences) and one of the points is close to 0/0/0 (e.g., 0/0/1e-10), we switch to
       * central differences instead
       */
      bool useCentral = false;
      Vector<Double> Pshifted_neg = Vector<Double>(dim_);
      Vector<Double> stepping_neg = Vector<Double>(dim_);

      Double steppingDistance;
      //std::cout << "++ Start with eval" << std::endl;
      if(startAtRemanence){
        // drive material above saturation in y
        vecIn.Init();
        vecIn[1] = maxAmplitude;
        vecOut = hystOperatorForTrace->computeValue_vec(vecIn, 0, true, false, successFlagForward);
        // then drop to remanence
        vecIn.Init();
        vecOut = hystOperatorForTrace->computeValue_vec(vecIn, 0, true, false, successFlagForward);
      }

      bool output = false;
      //std::cout << "++ Start loop" << std::endl;
      std::stringstream Jacwarnmsg;

      for(UInt i = 0; i < totalSteps; i++){
        vecIn.Init();
        // would be a workaround for peaking Jacobian issue around zero
//        if(scalarInputs[i] < scaling){
//          scalarInputs[i] = 1.0;
//        }
        vecIn[0] = scalarInputs[i];

//        if(((i >= 498)&&(i <= 500))||((i >= 1248)&&(i <= 1251))){
//          output = true;
//        } else {
//          output = false;
//        }

        if(output){
          std::cout << "++  " << i << std::endl;
          std::cout << "++  " << i << std::endl;
          std::cout << "++ Step " << i << std::endl;

          std::cout << "++++ PERMANENT STORAGE / BASE STEP  " << i << std::endl;
        }

        // bring material to next point first, then estimate Jacobian around this point
//        vecOut = hystOperatorForTrace->computeValue_vec(vecIn, 0, true, output, successFlagForward);
        bool debug_output = false;
        
//        output = true;
        Vector<Double> vecOutTMP;
        Vector<Double> vecOutTMPpost;
        if(output){
          // has to be evaluated before working on permanent storage
          // just for testing; do we get same output if working on tmp and perm storage?
          // 26.10.2020: works fine
          vecOutTMP = hystOperatorForTrace->computeValue_vec(vecIn, 0, false, debug_output, successFlagForward);
        }
        vecOut = hystOperatorForTrace->computeValue_vec(vecIn, 0, true, debug_output, successFlagForward);
        if(output){
          vecOutTMPpost = hystOperatorForTrace->computeValue_vec(vecIn, 0, false, debug_output, successFlagForward);
          std::cout << "In: = " << std::scientific<< std::setprecision(precisionDigits+1) << vecIn[0] << ", " << vecIn[1] << std::endl;
          std::cout << "Out (permStorage): = " << std::scientific<< std::setprecision(precisionDigits+1) << vecOut[0] << ", " << vecOut[1] << std::endl;
          std::cout << "Out (tmpStorage pre): = " << std::scientific<< std::setprecision(precisionDigits+1) << vecOutTMP[0] << ", " << vecOutTMP[1] << std::endl;
          std::cout << "Out (tmpStorage post): = " << std::scientific<< std::setprecision(precisionDigits+1) << vecOutTMPpost[0] << ", " << vecOutTMPpost[1] << std::endl;
          std::cout << "++  " << i << std::endl;
          std::cout << "++++ TEMPORAL STORAGE / JACOBI STEPS  " << i << std::endl;
        }
        output = false;
        
        vectorOutputs[i] = Vector<Double>(dim_);
        vectorOutputs[i] = vecOut;

        /*
         * Edit 25.10.2020
         * Sometimes the hysteresis operator shows strange behavior; evaluating 
         * the operator with the shifted input might jump way too much in both stepping directions;
         * however, the values in both shifted directions are fine in value, the problem seems to be
         * the value of vecOut; using central differences at such points leads to reasonable results for the
         * Jacobian; 
         * Why this is important: if fixed-point inversion shall be used, we take the maximal and minimal slope
         * of this tracing process as indicator for the local/global convergence factor to be used; if this value
         * is vastly overestimated, the convergence factor is way too small; in this case a convergence would work
         * but it will be so slow that the number of iterations will always exceed (example: the typical slopes are 
         * around 1e-5 for FeCo but at some peaks it goes up to 1e0)
         * Solution idea: Use forward backward differences as it is done currently; after setting up the Jacobian
         * check whether there is such a large jump in slopes; if so, repeat the evaluation with central differences
         * 26.10.2020: tryAgain works well but only if we reallly encounter spikes;
         * in some cases, however, the slopes creeped up continously; reason still
         * unknown; until then > forceCentral differences!
         */
        bool tryAgain = true;
        Double maxAllowedIncreaseFactor = 100.0;
        
        /*
         * estimate Jacobian by forward/backward OR central differences
         */
        useCentral = false; // per default to save runtime
        Jacobians[i] = Matrix<Double>(dim_,dim_);
        Matrix<Double> tmpJacAbs = Matrix<Double>(dim_,dim_);
        // for debugging; currently needed as output of hyst operator above (vecOut)
        // sometimes changes strongly if input is shifted slightly;
        // however, if input is slightly shifted into the other direction, too,
        // the difference between both shifed outputs is just fine; somehow vecOut
        // seems to be wrong > until error is found, take central differences
        // Example output: 
//        ++ Inner step 1 of outer step 913
//        UseCentral
//        Unshifted In: = -9.9358095600000e+05, 0.0000000000000e+00
//        Unshifted Out: = -2.9384986440569e+00, 6.4784209914851e-02
//        Shifted pos In: = -9.9358095600000e+05, 8.2800000000000e-05
//        Shifted pos Out: = -2.9384986440569e+00, 6.4981088529329e-02
//        Shifted neg In2: = -9.9358095600000e+05, -8.2800000000000e-05
//        Shifted neg Out2: = -2.9384986440569e+00, 6.4981088094017e-02
         // both shifted outputs lead to different values but the difference is only at 
         // the eigth-th digit; the diference to the unshifted input is however appearing
        // alread in the second digit > factor 1e6 higher!

        steppingDistance = std::max(scaling,scaling*vecIn.NormL2()/inputSat);
        if((vecIn.NormL2() < 1e-4*steppingDistance)){
          Jacwarnmsg << "# vecIn.NormL2() < 1e-4*steppingDistance (stepping distance = "<<steppingDistance<<") at step "<<i<<"; use central differences at this point." << std::endl;  
          useCentral = true;
        }
        if(forceCentral){
          useCentral = true;
        }

        while(tryAgain == true){
          for(UInt k = 0; k < dim_; k++){
//            output = true;
            if(output){
              std::cout << "++ Inner step " << k << " of outer step " << i << std::endl;
            }

            stepping = vecIn;
            stepping[k] += steppingDistance;

            debug_output = false;
  //          if((k == 0)&&(output)){
  //            debug_output = true;
  //          }

            // no overwrite here!
            Pshifted = hystOperatorForTrace->computeValue_vec(stepping, 0, false, debug_output, successFlagForward);

            Double usedSteppingDistance = steppingDistance;
            if(useCentral){
              if(output) std::cout << "UseCentral" << std::endl;
              usedSteppingDistance += steppingDistance;
              stepping_neg = vecIn;
              stepping_neg[k] -= steppingDistance;
              Pshifted_neg = hystOperatorForTrace->computeValue_vec(stepping_neg, 0, false, debug_output, successFlagForward);
            } else {
              if(output) std::cout << "UseForward" << std::endl;
              // to utilize the same structure further below, just assign vecOut and vecIn to Pshifted_neg and stepping_neg
              Pshifted_neg = vecOut;
            }

            if(output){
              std::cout << "Unshifted In: = " << std::scientific<< std::setprecision(precisionDigits+1) << vecIn[0] << ", " << vecIn[1] << std::endl;
              std::cout << "Unshifted Out: = " << std::scientific<< std::setprecision(precisionDigits+1) << vecOut[0] << ", " << vecOut[1] << std::endl;
              std::cout << "Shifted pos In: = " << std::scientific<< std::setprecision(precisionDigits+1) << stepping[0] << ", " << stepping[1] << std::endl;
              std::cout << "Shifted pos Out: = " << std::scientific<< std::setprecision(precisionDigits+1) << Pshifted[0] << ", " << Pshifted[1] << std::endl;
              if(useCentral){
                std::cout << "Shifted neg In2: = " << std::scientific<< std::setprecision(precisionDigits+1) << stepping_neg[0] << ", " << stepping_neg[1] << std::endl;
                std::cout << "Shifted neg Out2: = " << std::scientific<< std::setprecision(precisionDigits+1) << Pshifted_neg[0] << ", " << Pshifted_neg[1] << std::endl;
              }
            }
            output = false;

            // compute dP/dH
            for(UInt j = 0; j < dim_; j++){
              tmp[j][k] = (Pshifted[j] - Pshifted_neg[j])/usedSteppingDistance;
            }
          }
          
          if(i == 0){
            // the first Jacobian cannot be compared as there is no counterpart
            tryAgain = false;
          } else {
            // check tmp-Jacobian before inserting
            tmp.GetAbsValues(tmpJacAbs);
            Double tmpJacMaxVal = tmpJacAbs.GetMax();
          
            Jacobians[i-1].GetAbsValues(tmpJacAbs);
            Double tmpJacMaxValOld = tmpJacAbs.GetMax();
//            std::cout << "tmpJacMaxVal/tmpJacMaxValOld:" << tmpJacMaxVal/tmpJacMaxValOld << std::endl;
            if( (tmpJacMaxVal/tmpJacMaxValOld > maxAllowedIncreaseFactor) ){
              //|| (tmpJacMaxValOld/tmpJacMaxVal > maxAllowedIncreaseFactor) ){
              if(useCentral == false){
                Jacwarnmsg << "# Jump in Forward-Differences-Jacobian detected at step "<<i<<"! Try again with central differences." << std::endl;                
                std::cout<<"Jump in Forward-Differences-Jacobian detected at step "<<i<<"! Try again with central differences."<<std::endl;
                useCentral = true;
                tryAgain = true;
              } else {
                Jacwarnmsg << "# Jump in Central-Differences-Jacobian detected at step "<<i<<"!" <<std::endl;
                Jacwarnmsg << "# Current Jacobian: "<<tmp.ToString()<<std::endl;
                Jacwarnmsg << "# Previous Jacobian: "<<Jacobians[i-1].ToString()<<std::endl;
                Jacwarnmsg << "# Skip current Jacobian and replace it with the Jacobian from position "<<i-1<<std::endl;
                std::cout<<"Jump in Central-Differences-Jacobian detected at step "<<i<<"! Skip current Jacobian and replace it with the Jacobian from position "<<i-1<<std::endl;
                useCentral = false;
                tryAgain = false;
                tmp = Jacobians[i-1];
              }
            } else {
              useCentral = false;
              tryAgain = false;
            }
          }
        }
        Jacobians[i] = tmp;

        if(output){
          std::cout << "Jacobian:" << std::endl;
          std::cout << tmp.ToString();
          std::cout << std::endl;
        }

        if(i > 0){
          /*
           * estimated slope NOT used for determination of maximal and minimal slope;
           * just a debugging quantity!
           */
          estimatedSlope[i] = (vectorOutputs[i][0] - vectorOutputs[i-1][0])/(scalarInputs[i] - scalarInputs[i-1]);

//          if(output){
//            std::cout << "Previously estimated slope: " << estimatedSlope[i] << std::endl;
//            std::cout << "Jacobian at current point: " << tmp.ToString() << std::endl;
//            std::cout << "min(Jacobian): " << tmp.GetMin() << std::endl;
//            std::cout << "max(Jacobian): " << tmp.GetMax() << std::endl;
//          }

          if(estimatedSlope[i] < 0){
            std::stringstream warnmsg;
            if(abs(estimatedSlope[i]) > 1e-16){
              std::stringstream warnmsg;
              warnmsg << "Negative slope (" << estimatedSlope[i]<<") detected between " << scalarInputs[i] << " and " << scalarInputs[i-1] << std::endl;
              warnmsg << "Hysteresis model (Preisach) is assumed to be monoton (at least in the 1d case!)" << std::endl;
              WARN(warnmsg.str());
            }
//            warnmsg << "Negative slope (" << estimatedSlope[i]<<") detected between " << scalarInputs[i] << " and " << scalarInputs[i-1] << std::endl;
//            warnmsg << "Hysteresis model (Preisach) is assumed to be monoton (at least in the 1d case!)" << std::endl;
//            WARN(warnmsg.str());
          }
          if(abs(tmp.GetMax()) > maxSlope){
          //if(estimatedSlope[i] > maxSlope){
            maxSlope = abs(tmp.GetMax());
            //maxSlope = estimatedSlope[i];
            fieldMaxSlope = scalarInputs[i];
            idxMaxSlope = i;
          }
          if(abs(tmp.GetMin()) < minSlope){
//          if(estimatedSlope[i] < minSlope){
            minSlope = abs(tmp.GetMin());
//            minSlope = estimatedSlope[i];
            fieldMinSlope = scalarInputs[i];
            idxMinSlope = i;
          }

          if(vectorOutputs[i].Max() > maxPolarization){
            // only positive values checked here as we drive system beyond the positive polarization only
            maxPolarization = vectorOutputs[i].Max();
            fieldMaxPolarization = scalarInputs[i];
          }

          if(vectorOutputs[i][0] == 0){
            // excact hit; very unlikely
            negCoercivity = scalarInputs[i];
          }
          if((vectorOutputs[i][0] < 0) && (vectorOutputs[i-1][0] > 0)){
            // negative coercivity somewhere between the two inputs
            // > note: due to the form of the tested signal, we just hit the negative coercivity here and only have to check for it
            negCoercivity = (scalarInputs[i] + scalarInputs[i-1])/2.0;
          }
        }
      }

      traceTimer->Stop();
      Double endTime = traceTimer->GetCPUTime();

//      Double in = M_PI;
//      UInt precision = 9;
//      std::cout << std::scientific<< std::setprecision(precisionDigits+1) << "Pi up to " << precisionDigits+1 << " digits (due to output): " << in << std::endl;
//      
//      Double restrictedDouble = restrictPrecision(in, precision);
//      
//      std::cout << std::scientific<< std::setprecision(precisionDigits+1) << "Pi up to " << precision << " digits (due to function): " << restrictedDouble << std::endl;
//      
//      std::cout << "PRE-cutting" << std::endl;
//      std::cout << std::scientific<< std::setprecision(precisionDigits+1) << "###MaxSlope### " << maxSlope << std::endl;
//      std::cout << std::scientific<< std::setprecision(precisionDigits+1) << "###MinSlope### " << minSlope << std::endl;
//      std::cout << std::scientific<< std::setprecision(precisionDigits+1) << "###NegCoercivity### " << negCoercivity << std::endl;
//      std::cout << std::scientific<< std::setprecision(precisionDigits+1) << "###MaxPolarization### " << maxPolarization << std::endl;
//      

      /*
       * Restrict prescision via stringstreams
       */
      maxSlope = restrictPrecision(maxSlope,precisionDigits);
      minSlope = restrictPrecision(minSlope,precisionDigits);
      negCoercivity = restrictPrecision(negCoercivity,precisionDigits);
      maxPolarization = restrictPrecision(maxPolarization,precisionDigits);      

      // store for later usage with mutex for thread-safe access
      currentMaterial.maxSlope_ = maxSlope;
      currentMaterial.minSlope_ = minSlope;
      currentMaterial.negCoercivity_ = negCoercivity;
      currentMaterial.maxPolarization_ = maxPolarization;
      {
        std::lock_guard<std::mutex> lock(tracedDataMutex_);
        CoefFunctionHyst::tracedOperatorData_[tracedData_FILENAME.str()] = currentMaterial;
      }

      bool testOutput = true;
      if(testOutput){
        //std::cout << "Traced hyst operator for material " << material_->GetName() << std::endl;
        tracedData.open(tracedData_FILENAME.str(),std::ios_base::out);
        tracedData << "### INFO " << std::endl;
        tracedData << "# > hysteresis operator was tested with 1d input signal " << std::endl;
		if(enforcedPreisachResolution > 0){
			UInt preisachRes;
			if(dedicatedOperatorForStrains){
			  preisachRes = STRAIN_weightParamsForTesting_.numRows_;
			} else {
			  preisachRes = POL_weightParamsForTesting_.numRows_;
			}
			tracedData << "# > used resolution for Preisach plane as enforced in input.xml: " << preisachRes << " x " << preisachRes << std::endl;
        } else {
			UInt preisachRes;
			if(dedicatedOperatorForStrains){
			  preisachRes = STRAIN_weightParams_.numRows_;
			} else {
			  preisachRes = POL_weightParams_.numRows_;
			}
			tracedData << "# > used resolution for Preisach plane as given in mat.xml: " << preisachRes << " x " << preisachRes << std::endl;
		}
		if(dedicatedOperatorForStrains){
          tracedData << "# > hysteresis operator will be used ONLY for evaluation of irreversible strains " << std::endl;
        } else {
          tracedData << "# > hysteresis operator will be used for BOTH the evaluation of polarization and irreversible strains " << std::endl;
        }
        if(forceCentral){
          tracedData << "# > central differences were enforced" << std::endl;
        } else {
          tracedData << "# > forward differences were used whenever possible" << std::endl;
        }
        tracedData << "# > minimal stepping distance for Jacobian computation (parameter: trace_JacResolution_) = " << trace_JacResolution_ << std::endl;       
        tracedData << "# > maximal absolute entry of Jacobian (used as maxSlope)  " << std::scientific<< std::setprecision(precisionDigits+1)<< maxSlope << "\t (found at (idx,Field) = (" << idxMaxSlope << "," << fieldMaxSlope << ") )" << std::endl;
        tracedData << "# > minimal absolute entry of Jacobian (used as minSlope)  " << std::scientific<< std::setprecision(precisionDigits+1)<< minSlope << "\t (found at (idx,Field) = (" << idxMinSlope << "," << fieldMinSlope << ") )" << std::endl;
        tracedData << "# > maximal polarization  " << std::scientific<< std::setprecision(precisionDigits+1)<< maxPolarization << "\t (found for Field " << fieldMaxPolarization << ")" << std::endl;
        tracedData << "# > specified saturation values for field and polarization  " << std::scientific<< std::setprecision(precisionDigits+1)<< inputSat << ", " << outputSat << std::endl;
        tracedData << "# > negative coercivity found for Field " << std::scientific<< std::setprecision(precisionDigits+1)<< negCoercivity << std::endl;
        tracedData << "# > number of traced positions " << totalSteps << std::endl;
        tracedData << "# > required runtime " << std::scientific<< std::setprecision(precisionDigits+1)<< endTime-startTime << std::endl;
        tracedData << "# Warnings and errors during evaluation:" << std::endl;
        tracedData << "#" << std::endl;
        tracedData << Jacwarnmsg.str();
        tracedData << "#" << std::endl;
        tracedData << "# Parameter to be read-in/reused by CFS" << std::endl;
        tracedData << "###MaxSlope### " << std::scientific<< std::setprecision(precisionDigits+1)<< maxSlope << std::endl;
        tracedData << "###MinSlope### " << std::scientific<< std::setprecision(precisionDigits+1)<< minSlope << std::endl;
        tracedData << "###NegCoercivity### " << std::scientific<< std::setprecision(precisionDigits+1)<< negCoercivity << std::endl;
        tracedData << "###MaxPolarization### " << std::scientific<< std::setprecision(precisionDigits+1)<< maxPolarization << std::endl;
        tracedData << "#" << std::endl;
        tracedData << "# Traced data: " << std::endl;
        tracedData << "#IDX \t Field_x \t Pol_x \t Pol_y \t Jac_xx \t Jac_xy \t Jac_yx \t Jac_yy" << std::endl;
        for(UInt i = 0; i < totalSteps; i++){
          tracedData << i << "\t" << std::scientific<< std::setprecision(precisionDigits+1) << scalarInputs[i] << "\t" << vectorOutputs[i][0] << "\t" << vectorOutputs[i][1] << "\t" << Jacobians[i][0][0] << "\t" << Jacobians[i][0][1] << "\t" << Jacobians[i][1][0] << "\t" << Jacobians[i][1][1] << std::endl;
        }
        tracedData.close();
        std::cout << std::scientific<< std::setprecision(precisionDigits+1) << "###MaxSlope### " << maxSlope << std::endl;
        std::cout << std::scientific<< std::setprecision(precisionDigits+1) << "###MinSlope### " << minSlope << std::endl;
        std::cout << std::scientific<< std::setprecision(precisionDigits+1) << "###NegCoercivity### " << negCoercivity << std::endl;
        std::cout << std::scientific<< std::setprecision(precisionDigits+1) << "###MaxPolarization### " << maxPolarization << std::endl;
      }

      /*
       * free memory
       */
//      std::cout << "clean up" << std::endl;
      delete hystOperatorForTrace;
      delete traceTimer;
//      std::cout << "clean up 2" << std::endl;
      delete[] vectorOutputs;
      delete[] Jacobians;

    //  std::cout << "Stop test here; exit!" << std::endl;
//      exit(0);

    }


  
}// namespace
