#include "VectorPreisachMayergoyz.hh"

#include <iostream>
#include <fstream>
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Utils/Timer.hh"

namespace CoupledField
{
  class Preisach;

  DEFINE_LOG(vecpreisachmayergoyz, "vecpreisachmayergoyz")

  VectorPreisachMayergoyz::VectorPreisachMayergoyz(Integer numElem, Vector<Double> xSat, Vector<Double> ySat,
          Matrix<Double>* preisachWeight, UInt dim, bool isVirgin,
          Vector<Double> anhyst_A, Vector<Double> anhyst_B, Vector<Double> anhyst_C,
          bool anhystOnly, int clipOutput) : Hysteresis(numElem){

    EXCEPTION("Anisotropic case not yet implemented; Find way to compute weightings first!");

  }

  VectorPreisachMayergoyz::VectorPreisachMayergoyz(Integer numElem, ParameterPreisachOperators& operatorParams,
          ParameterPreisachWeights weightParams, UInt dim, bool isVirgin)
  : Hysteresis(numElem,operatorParams,weightParams){

//  (Integer numElem, UInt numDirections, Double xSat, Double ySat,
//          Matrix<Double>& preisachWeight, UInt dim, bool isVirgin,
//          Double anhyst_A, Double anhyst_B, Double anhyst_C, bool anhystOnly, int clipOutput) : Hysteresis(numElem){
    UInt numDirections = operatorParams.numDirections_;
    int clipOutput = operatorParams.outputClipping_;

    startingAxis_ = operatorParams.startingAxisMG_;

    /*
     * Idea of Mayergoyz model:
     *  Vectorial Hysteresis is obtained by summing up Scalar Hysteresis models into
     *  a (finite) number of directions; for the isotropic case, all these scalar functions
     *  are the same;
     * NOTE: the Preisach weights for the single scalar models do not correspond to the
     *  weights of a single preisach model but have to be adapted (see helper function)
     *
     */
    dim_ = dim;
//    if(dim != 2){
//      EXCEPTION("Mayergoyz vector model (currently) only implemented for 2d");
//      // for the 3d case the computation of the weightings has to be adapted; furtherone
//      // one has to specify many more directions as we have to integrate over a whole
//      // hemisphere instead of a hemicircle
//    }

    numDirections_ = numDirections;
    if(numDirections < 2){
      EXCEPTION("To obtain vector functionality, at least 2 directions are required");
    }

    if(dim_ == 3){
      // TODO: allow direct setting of theta and phi directions via mat file;
      Double sqrtRound = std::round(std::sqrt(numDirections));
      numDirectionsTheta_ = UInt(sqrtRound);
      numDirectionsPhi_ = UInt(numDirections/sqrtRound);
//      std::cout << "3D Isotropic Mayergoyz Vector Model with N_theta = " << numDirectionsTheta_ << " and N_phi = " << numDirectionsPhi_ << std::endl;
      numDirections_ = numDirectionsTheta_*numDirectionsPhi_;
      operatorParams.numDirections_ = numDirections_;
    }

    singleDirections_ = new Vector<Double>[numDirections_];

    if(startingAxis_.GetSize() != dim_){
      EXCEPTION("Length of starting axis does not match dimension of model");
    }

    if(startingAxis_.NormL2() < 1e-16){
//      std::cout << "Generate random starting axis" << std::endl;
      // generate random starting axis
      for(UInt i = 0; i < startingAxis_.GetSize(); i++){
        startingAxis_[i] = 2*((float) rand()) / (float) RAND_MAX-1.0;
      }
      Double axisNorm = startingAxis_.NormL2();
      if(axisNorm != 0){
        for(UInt i = 0; i < startingAxis_.GetSize(); i++){
          startingAxis_[i] /= axisNorm;
        }
      } else {
        startingAxis_[0] = 1.0; // default if nothing works
      }
    }

    if(dim_ == 2){
      // for 2D > only semicircle required due to symmetry
      deltaAngle_ = M_PI/numDirections;

      // accroding to "Magnetic Field Analysis of Electric Machines Taking Ferromagnetic Hysteresis into Account" by J. Saitz, p. 38
      // we should use randomly distributed starting angle to make up for the missing symmetry property due to discretization
      Double angleOffset = atan2(startingAxis_[1],startingAxis_[0]);

  //    Double startingAngle = -deltaAngle/2;
      Double currentAngle = 0.0;
  //    std::cout << "MG-model-Directions of scalar models: " << std::endl;
  //    std::cout << "Staring axis: " << startingAxis_.ToString() << std::endl;
  //    std::cout << "Staring angle: " << angleOffset << std::endl;

      // used coordinate transform:
      //  x = r cos(phi)
      //  y = r sin(phi)
      for(UInt i = 0; i < numDirections; i++){
        currentAngle = angleOffset + i*deltaAngle_;
        singleDirections_[i] = Vector<Double>(dim_);
        singleDirections_[i][0] = std::cos(currentAngle);
        singleDirections_[i][1] = std::sin(currentAngle);

  //      std::cout << "angle (" << i << "): " << currentAngle << std::endl;
  //      std::cout << "direction (" << i << "): " << singleDirections_[i].ToString() << std::endl;
      }
    } else {
      // for 3D > only one hemisphere required due to symmetry
      // NOTE: we just have to keep track of sin(theta) as integration requires sin(theta)dTheta; otherwise angles are
      //        not important for us
      deltaAnglePhi_ = 2*M_PI/numDirectionsPhi_;
      deltaAngleTheta_ = M_PI/(2.0*numDirectionsTheta_);
      sinTheta_ = Vector<Double>(numDirections_);

      Double angleOffsetPhi = atan2(startingAxis_[1],startingAxis_[0]);
      Double angleOffsetTheta = acos(startingAxis_[2]); // r is normalized to 1 > see above

      Double phi,theta;
      UInt cnt = 0;
      Double s,c;

      // used coordinate transform:
      // x = r cos(phi) sin(theta)
      // y = r sin(phi) sin(theta)
      // z = r cos(theta)
      for(UInt i = 0; i < numDirectionsTheta_; i++){
        theta = angleOffsetTheta + i*deltaAngleTheta_;
        c = std::cos(theta);
        s = std::sin(theta);
        sinTheta_[cnt] = s;

        for(UInt j = 0; j < numDirectionsPhi_; j++){
//          std::cout << "Index " << cnt << " of " << numDirections_ << std::endl;
          phi = angleOffsetPhi + j*deltaAnglePhi_;
          singleDirections_[cnt] = Vector<Double>(dim_);
          singleDirections_[cnt][0] = s*std::cos(phi);
          singleDirections_[cnt][1] = s*std::sin(phi);
          singleDirections_[cnt][2] = c;
          cnt++;
        }
      }
    }

    bool printDirections = !true;
    if(printDirections){
      std::cout << "Using Mayergoyz Vector Preisach Model in " << dim_ << "D with the following " << numDirections_ << " directions: " << std::endl;
      for(UInt i = 0; i < numDirections_; i++){
        std::cout << "direction (" << i << "): " << singleDirections_[i].ToString() << std::endl;
      }
    }
    clipOutput_ = clipOutput;
    isIsotropic_ = true;

    if(isIsotropic_){
      // add anhysteretic part directly to output and set saturation parameter accordingly
      // the single scalar models will not have anhysteretic parts
      // already done in Hysteresis base class
//      XSaturated_ = xSat;
//      PSaturated_ = ySat;
//      anhyst_A_ = anhyst_A;
//      anhyst_B_ = anhyst_B;
//      anhyst_C_ = anhyst_C;
//      anhystOnly_ = anhystOnly;

      /*
       * IMPORTANT REMARK:
       *  > although we use scalar models here, we are not alloewed to directly apply the preisach parameter for
       *     the scalar case (i.e. the weights, the anhyst parameter and so on)
       *  > make sure that the passed parameter are already transformed correctly
       *      > see CoefFunctionHyst
       */
      singlePreisachOperators_ = new Preisach*[numDirections];
      for(UInt i = 0; i < numDirections_; i++){
        // > set anhysteretic parameter to 0 here!
        // preisach models shall only return the pure hysteretic part; the anhysteretic part is added later in this class
        bool ignoreAnhystPart = true;
        singlePreisachOperators_[i] = new Preisach(numElem,operatorParams,weightParams, isVirgin, ignoreAnhystPart);
        singlePreisachOperators_[i]->setFixDirection(singleDirections_[i]);
      }

    } else {
      // in each single direction, we need specifiy values for
      // xSat, ySat, preisachWeights, anhyst_A, anhyst_B, anhyst_C
      // > unless we have a way to measure or compute these values, the anisotropic case is not available
      EXCEPTION("Only isotropic case implemented.");
    }

    bool testSatOutput = false;
    if(testSatOutput){
      std::cout << "hystSaturated: " << hystSaturated_ << std::endl;
      std::cout << "PSaturated: " << PSaturated_ << std::endl;
      std::cout << "anhyst_A_ " << anhyst_A_ << std::endl;
      std::cout << "anhyst_B_ " << anhyst_B_ << std::endl;
      std::cout << "anhyst_C_ " << anhyst_C_ << std::endl;

      Vector<Double> satInput = Vector<Double>(dim_);
      Vector<Double> vectorResult = Vector<Double>(dim_);
      std::cout << "direction [0]: " << singleDirections_[0].ToString() << std::endl;
      satInput.Init();
      satInput.Add(XSaturated_,singleDirections_[0]);

      int blub;
      Double scalarResult = singlePreisachOperators_[0]->computeValueAndUpdate(XSaturated_,0,false,blub);
      std::cout << "result of scalar model in direction [0]: " << scalarResult << std::endl;
      vectorResult.Init();
      vectorResult.Add(scalarResult,singleDirections_[0]);
      std::cout << "resulting vector in direction [0]: "<< vectorResult.ToString() << std::endl;

      vectorResult = this->computeValue_vec(satInput,0,false,false,blub,true);
      std::cout << "result of vector model in directon [0] (without anhyst part!): " << vectorResult.ToString() << std::endl;

      vectorResult = this->computeValue_vec(satInput,0,false,false,blub);
      std::cout << "result of vector model in directon [0] (with anhyst part!): " << vectorResult.ToString() << std::endl;

      /*
       * NOTE: we cannot get the same result here! reason: we compare the vector result with one of the internally stored
       *      scalar models; these scalar models use the adapted weights, however! for a correct comparison, we would have
       *      to compare the vector model (with adapted weights) to a scalar model with original weights!
       * > can only be done in coefFunctionHyst where the actual transformation of weights is done
       */
    }

    prevXVal_ = new Vector<Double>[numElem];
    prevHVal_ = new Vector<Double>[numElem];
    for(int k = 0; k < numElem; k++){
      prevXVal_[k] = Vector<Double>(dim_);
      prevXVal_[k].Init();
      prevHVal_[k] = Vector<Double>(dim_);
      prevHVal_[k].Init();
    }

  }

  VectorPreisachMayergoyz::~VectorPreisachMayergoyz(){
    delete [] singlePreisachOperators_;
    delete [] singleDirections_;

    delete[] prevXVal_;
    delete[] prevHVal_;
  }

  Vector<Double> VectorPreisachMayergoyz::computeValue_vec(Vector<Double>& xVal, Integer idx,
          bool overwrite,bool debugOutput,int& successFlag, bool skipAnhystPart){

        /*
     * Vectorial output = integral/sum over scalar models
     */
    Vector<Double> output = Vector<Double>(dim_);
    output.Init();

//    Vector<Double> tmp = Vector<Double>(dim_);
//    Vector<Double> currentDir;
//    Double scalarInput, scalarOutput;

    /*
     * Remarks to capping/clipping:
     * A) if input is only along one direction and the initial state is the virgin state,
     *    or has only components in this particular direction
     *    then, the output of the vector model will be equal to the scalar one for all input values
     *    if amplitude of input is clipped to saturation as noted above
     *    WITHOUT capping, equality only holds up to saturation as above saturation, more and more
     *      scalar models reach saturation, too
     * B) if input is only along one direction but during the past, different directions occured,
     *    then, the output of the vector model will not be equal to the scalar one due to remanent
     *    contributions of the old state
     * C) if capping is applied, remanent contributions cannot be erased completely if mutually orthogonal
     *    input is applied as due to the capping the input to the non-aligned directions will be scaled down
     *    by cos(angle) and thus can never reach saturation;
     *    WITHOUT capping, the remanent contribution can be erased if input goes above saturation
     *
     * From "Mathematical models of hysteresis and their applications" - Mayergoyz:
     * I) vector model should reduce to scalar model if input varies only in one direction; however, it is
     *      not stated if this equality has to hold above saturation
     * II) for mutually orthogonal input it MUST be possible to erase remamentn parts orthogonal to the current
     *      input if input amplitude goes to infinity; this implies that NO capping is done by Mayergoyz in his
     *      implementation;
     *
     * Conclusions:
     *  1. no capping of input
     *      > remanent parts can be erased
     *      > equality to scalar model only up to saturation
     *  2. idea: clip output
     *      > equaliy to scalar model beyond saturation but difference when leaving saturation again
     *        (as Preisach planes have been set inside the single scalar models)
     */

    successFlag = 0;
    int successFlagSingle = 0;
    Vector<Double> scalContribution = Vector<Double>(dim_);
    Vector<Double> testOutput = Vector<Double>(dim_);
    for(UInt i = 0; i < numDirections_; i++){
      scalContribution.Init();
      scalContribution = singlePreisachOperators_[i]->computeValue_vec(xVal,idx,overwrite,debugOutput,successFlagSingle);

      if(dim_ == 3){
        scalContribution *= sinTheta_[i];
      }

      output.Add(1.0,scalContribution);

//      // old version
//        currentDir = singleDirections_[i];
//        currentDir.Inner(xVal,scalarInput);
//        scalarOutput = singlePreisachOperators_[i]->computeValueAndUpdate(scalarInput,idx,overwrite,successFlagSingle);
//        output.Add(scalarOutput,currentDir);

      // just check if at least one of the scalar models required reevaluation
      if(successFlagSingle != 0){
        successFlag = successFlagSingle;
      }
    }

    //Double deltaAngle = M_PI/numDirections_;
    // for numerical integration we have to multiply by deltaAngle
    // and to average out correctly over the halfspace we have to multiply by 2.0/Pi
    // > in total, multiply by 2/numDir
		// Question: why do we average out in the first place?
    //output.ScalarMult(2.0/numDirections_);
    // actually M_PI/numDirections works better ...
//		output.ScalarMult(M_PI/numDirections_);
    if(dim_ == 2){
      output.ScalarMult(deltaAngle_);
    } else {
      // multiply with dTheta dPhi
      // Note: sin(theta) already considered above
      output.ScalarMult(deltaAngleTheta_*deltaAnglePhi_);
    }

    Vector<Double> dirInput = Vector<Double>(dim_);
    dirInput.Init();
    if((xVal.NormL2() != 0)){ //&&(clipOutput_ == 2)){
      dirInput.Add(1.0/xVal.NormL2(),xVal);
    }

    if(clipOutput_ == 1){
      // clip amplitude to saturation; works well if input only in 1d but
      // not so well if remanent parts perpendicular to input exist as those
      // will be scaled down, too
      // NOTE: here we have to scale by hystSaturated as this is the saturation
      //      value of the hysteretic part alone
      if(output.NormL2() > hystSaturated_){
        output.ScalarMult(hystSaturated_/output.NormL2());
      }
    } else if(clipOutput_ == 2){
      // > default
      // clip amplitude to saturation, but such that remanent part is not
      // affected; results seem to be more reasonable than unclipped and clipping 1
      // BIG DISADVANTAGE: hard to invert!!!!!
      Double projection = output.Inner(dirInput);
      if((abs(output.NormL2()) > hystSaturated_)&&(dirInput.NormL2() != 0)){
        output.Add(-projection,dirInput);
        Double normRemaining = output.NormL2();
        output.Add(std::sqrt(hystSaturated_*hystSaturated_-normRemaining*normRemaining),dirInput);
      }
    }
    if(skipAnhystPart == false){
      if(isIsotropic_){
        // for isotropic case, add anhystPart directly to output
        // make sure that the scalar models return no anhystPart in this case
        // > scale anhysteretic part by PSaturated_
        if(xVal.NormL2() != 0){
          Double amplitude = PSaturated_*evalAnhystPart_normalized(xVal.NormL2()/XSaturated_);
          output.Add(amplitude,dirInput);
        }
      }
    }
    if(overwrite == true){
      /*
       * store to arrays > as in VectorPreisachSutor
       * include anhyst part!
       */
      prevXVal_[idx] = xVal;
      prevHVal_[idx] = output;
    }

    return output;
  }
}
