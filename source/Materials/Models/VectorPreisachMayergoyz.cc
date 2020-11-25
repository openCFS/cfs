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
          Vector<Double> anhyst_A, Vector<Double> anhyst_B, Vector<Double> anhyst_C, Vector<Double> anhyst_D,
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

    numDirections_ = numDirections;
    if(numDirections_ < 2){
      EXCEPTION("To obtain vector functionality, at least 2 directions are required!");
    }

    if(dim_ == 3){
      // possible TODO: allow direct setting of theta and phi directions via mat file;
      // Edit: 3.6.2020 - Phi should be resolved with 4x as many directions as Theta
      // > until now the sqrt was taken 
      bool useOldDistribution=false;
      if(useOldDistribution){
        Double sqrtRound = std::round(std::sqrt(numDirections));
        numDirectionsTheta_ = UInt(sqrtRound);
        if(numDirectionsTheta_ < 2){
          EXCEPTION("To obtain vector functionality in 3d, at least 2 directions are required into the third dimension! Increase number of directions in mat.xml. "
                  "The number of directions into the third dimensions is computed via floor(sqrt(totalNumberDirections)).");
        }
        numDirectionsPhi_ = UInt(numDirections/sqrtRound);
      } else {
        UInt sqrtCeil = std::ceil(std::sqrt(numDirections/4.0));
        numDirectionsTheta_ = sqrtCeil;
        if(numDirectionsTheta_ < 2){
          EXCEPTION("To obtain vector functionality in 3d, at least 2 directions are required into the third dimension! Increase number of directions in mat.xml. "
                  "The number of directions into the third dimensions is computed via ceil(sqrt(totalNumberDirections/4)) so that resolution in both angular directions is equal.");
        }
        numDirectionsPhi_ = std::ceil(Double(numDirections)/Double(numDirectionsTheta_));
        //std::cout << "###NumTheta x NumPhi### " << numDirectionsTheta_ << " x " << numDirectionsPhi_ << std::endl;
      }
      if(numDirectionsPhi_ < 2){
        EXCEPTION("To obtain vector functionality in 3d, at least 2 directions are required into the each dimension! Increase number of directions in mat.xml.");
      }
      numDirections_ = numDirectionsTheta_*numDirectionsPhi_;
      numDirections = numDirections_; // just for the case that I missed it somewhere below to replace numDirections with numDirections_
    } 
//    else {
//      std::cout << "###NumPhi### " << numDirections_ << std::endl;
//    }
    operatorParams.numDirections_ = numDirections_;
    singleDirections_ = new Vector<Double>[numDirections_];

    if(startingAxis_.GetSize() != dim_){
      EXCEPTION("Length of starting axis does not match dimension of model");
    }

    if(startingAxis_.NormL2() < 1e-16){
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
      deltaAngle_ = M_PI/numDirections_;

      // accroding to "Magnetic Field Analysis of Electric Machines Taking Ferromagnetic Hysteresis into Account" by J. Saitz, p. 38
      // we should use randomly distributed starting angle to make up for the missing symmetry property due to discretization
      Double angleOffset = atan2(startingAxis_[1],startingAxis_[0]);

      Double currentAngle = 0.0;
  //    std::cout << "MG-model-Directions of scalar models: " << std::endl;
  //    std::cout << "Staring axis: " << startingAxis_.ToString() << std::endl;
  //    std::cout << "Staring angle: " << angleOffset << std::endl;

      // used coordinate transform:
      //  x = r cos(phi)
      //  y = r sin(phi)
      for(UInt i = 0; i < numDirections_; i++){
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
      // add anhysteretic part directly to output and set saturation parameter accordingly;
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
      singlePreisachOperators_ = new Preisach*[numDirections_];
      for(UInt i = 0; i < numDirections_; i++){
        // > set anhysteretic parameter to 0 here!
        // preisach models shall only return the pure hysteretic part; the anhysteretic part is added later in this class
        bool ignoreAnhystPart = true;
        singlePreisachOperators_[i] = new Preisach(numElem,operatorParams,weightParams, isVirgin, ignoreAnhystPart);
        singlePreisachOperators_[i]->setFixDirection(singleDirections_[i]);
      }

    } else {
      // anisotropic case:
      // in each single direction, we need specifiy values for
      // xSat, ySat, preisachWeights, anhyst_A, anhyst_B, anhyst_C
      // > unless we have a way to measure or compute these values, the anisotropic case is not available
      EXCEPTION("Only isotropic case implemented.");
    }

    // new for rotational loss computation
    // > extended Mayergoyz model
    // cf. Dlala - "Improving Loss Properties of the Mayergoyz Vector Hysteresis Model"
    lossParam_a_ = operatorParams.lossParam_a;
    lossParam_b_ = operatorParams.lossParam_b;
    improveRotLoss_ = false;
    if( (lossParam_a_ != 0)||(lossParam_b_ != 0) ){
      improveRotLoss_ = true;
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
    for(UInt i = 0; i < numDirections_; i++ ){
      delete singlePreisachOperators_[i];
    }
    delete [] singlePreisachOperators_;  
    delete [] singleDirections_;

    delete[] prevXVal_;
    delete[] prevHVal_;
  }

  Vector<Double> VectorPreisachMayergoyz::computeInput_vec_Everett(Vector<Double> yVal, Integer operatorIndex,
          Matrix<Double> mu, bool overwrite, int& successFlag){

//    computeValue_vec(Vector<Double>& xVal, Integer idx,
//          bool overwrite,bool debugOutput,int& successFlag, bool skipAnhystPart)
    
     /*
      * Added 7.5.2020
      * for reference see e.g., Dlala 2010 "Improving Loss Properties of the Mayergoyz Vector Hysteresis Model"
      * and Saitz 2001 "Magnetic Field Analysis of Electric Machines Taking Ferromagnetic Hysteresis into Account"
      * 
      * Idea: Directly build up an inverse form of the Mayergoyz vector model by replacing the sum over forward models
      * B_phi(H_phi) with backward models H_phi(B_phi). As we have no inverse Preisach model (currently) implemented, we
      * instead utilize the forward scalar model and invert it via the efficient Everett-based algorithm (see Preisach.cc).
      * Important note: This kind of evaluation will actually be able to retrieve a suitable H_vec from a given B_vec but
      * this H_vec fed into the forward algorithm (computeValue_vec) will not return the original B_vec!
      * Reason: The inverse of a sum is not the sum of its inverse terms (at least not in general). This means, if we want
      * to have the polarization P_vec, too, we may not just pass the computed H_vec to the forward evaluation. Luckily, it is 
      * even easier than that. We just have to subtracht mu time H_vec from B_vec to retrieve P_vec
      * 
      */     

    EXCEPTION("Everett based inversion of the Mayergoyz model not supported yet as its output cannot be used in the forward model");
    // Problem: Throughout coefFunctionHyst, we need both, the forward and the backward/inverse Hysteresis model
    // the backward model is required, if B is known, but sometimes H is known or set, e.g., in testing of the hyst
    // operator, for the computation of the Jacobian and for tracing of the hyst operator; in these cases H is set 
    // and P needs to be retrieved; to retrieve P from H we need the inverse of the backward model 
    // which unfortunately is NOT the forward Mayergoyz model!
    // the issue lies in the fact, that the sum over the inverted scalar models is not the same as the inverse of the
    // summed up forward scalar models; 
    
    /*
     * Vectorial computedInput = integral/sum over inverse scalar models
     */
    Vector<Double> computedInput = Vector<Double>(dim_);
    computedInput.Init();
    
    successFlag = 0;
    int successFlagSingle = 0;
    Vector<Double> scalContribution = Vector<Double>(dim_);
//    Vector<Double> testOutput = Vector<Double>(dim_);
    Vector<Double> currentDirection;
    Vector<Double> tmp = Vector<Double>(dim_);
    Double muForInversion;
    for(UInt i = 0; i < numDirections_; i++){     
      // get direction of current scalar model
      currentDirection = singleDirections_[i];
      
      // extract correct element of material tensor
      mu.Mult(currentDirection,tmp);
      currentDirection.Inner(tmp,muForInversion);
      
      /*
       * Compute scalar input to hyst operator
       */
      Double scalInput, scalOutput;
      currentDirection.Inner(yVal,scalOutput);
      scalInput = singlePreisachOperators_[i]->computeInputAndUpdate(scalOutput,muForInversion, operatorIndex, overwrite, successFlagSingle);

      scalContribution.Init();
      scalContribution.Add(scalInput,currentDirection);

      if(dim_ == 3){
        scalContribution *= sinTheta_[i];
      }

      computedInput.Add(1.0,scalContribution);

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
      computedInput.ScalarMult(deltaAngle_);
    } else {
      // multiply with dTheta dPhi
      // Note: sin(theta) already considered above
      computedInput.ScalarMult(deltaAngleTheta_*deltaAnglePhi_);
    }

    return computedInput;
  }
  
  Vector<Double> VectorPreisachMayergoyz::computeValue_vec(Vector<Double>& xVal, Integer idx,
          bool overwrite,bool debugOutput,int& successFlag, bool skipAnhystPart){

    /*
     * Vector output = integral/sum over scalar models
     */
    Vector<Double> output = Vector<Double>(dim_);
    output.Init();

    successFlag = 0;
    int successFlagSingle = 0;
    Vector<Double> scalContribution = Vector<Double>(dim_);
    
    // compute correction angle; 
    // if it is 0, scalInput is just the dotProduct betwen xVal and the direction
    // > Note that this angle does not depend on the direction of the scalar model, so we just have to evaluate it once
    Double angleCorrection = 0.0;
    if(improveRotLoss_ == true){
      // edit 19.10.2020: 
      // angleCorrection restricted to angle between prevHVal_ and xVal (previously restriction to 2pi)!
      // > see also comments in evaluateLagCorrectionAngle

      Double phaseLagMax = 2*M_PI;
      
      // restrict angle to angle between X and Y; as Y is not known yet, restrict to previous
      // output of hyst operator H instead (Y = mu*X + H)
      phaseLagMax = std::abs(xVal.InnerAngle(prevHVal_[idx]));
    
      angleCorrection = evaluateLagCorrectionAngle(xVal, prevXVal_[idx], phaseLagMax);
    }

    for(UInt i = 0; i < numDirections_; i++){
      scalContribution.Init();
      if(improveRotLoss_ == false){
        scalContribution = singlePreisachOperators_[i]->computeValue_vec(xVal,idx,overwrite,debugOutput,successFlagSingle);
      } else {
        Double scalInput;
        Double scalOutput;

        Double absProd = xVal.NormL2()*singleDirections_[i].NormL2();
        /*
         * Important edit 15.10.2020
         * by the old implementation we get the absolute value of the angle between
         * the two vectors; as long as we later use just the cos of this angle, it
         * makes no difference as cos is a symmetric function w.r.t. the angle
         * however, as we modify the angle by subtracting the lag angle, the sign
         * will play an important role!
         * > use new implementation in vector class
         */      
        Double baseAngle = 0;
        baseAngle = singleDirections_[i].InnerAngle(xVal);
        
        scalInput = absProd*std::cos(baseAngle - angleCorrection);
        scalOutput = singlePreisachOperators_[i]->computeValueAndUpdate(scalInput,idx, overwrite,successFlagSingle);
        scalContribution.Add(scalOutput,singleDirections_[i]);        
      }

      if(dim_ == 3){
        scalContribution *= sinTheta_[i];
      }

      output.Add(1.0,scalContribution);

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

  Double VectorPreisachMayergoyz::evaluateLagCorrectionAngle(Vector<Double>& xVal, Vector<Double>& prevXVal, Double maxAngle){   
    /*
     * Compute phase shift/angle that can be considered during the input projection of xVal onto the directions
     * of the scalar models to imporve the rotational loss properties of the vector model.
     * Source: Dlala - "Improving Loss Properties of the Mayergoyz Vector Hysteresis Model"
     * 
     * psi(t) = g(x)|x||dPhi|
     * g(x) = e^(-a|x| + b))
     * 
     * Note: 
     *  own modification: work with normalized xVal; if x = H, H might be in the order 1e5 or higher resulting
     *  in extremely large values; in the original paper, the authors only show the results for x = B where a
     *  typical B does not exceed 1-2T
     * 
     * Note 2 (19.06.2020, 20.10.2020):
     *  after matching 'a' and 'b' to rotloss curves, 'a' became negative; thus, the exponential function does not tend
     *  to 0 for large inputs but instead goes to infinity; as long as we stay in the range around saturation
     *  the normalized values as used above work but once we hit a critical input level, the exponential function will
     *  result in nan and from then on, computations are only non-sense;
     *  in general it seems to be non-sense if the computed correction is larger than the actual angle between the
     *  old and the new direction; the maximum lag should be this angle between old state and new state;
     *  > this aspect seems to be intended in the cited paper above (cf. Fig 10); the intention between the lag angle is
     *    to rotate the coordinate systems of the input or output in such a way, that input and output align
     *  > Consequence 1: all scalar directions have to be rotated by same value
     *  > Consequence 2: we have to make sure that all projections of the input vector on the single direction vectors
     *    follows the same turn; in the previous implementation only the absolute value of the angle was retrieved which
     *    then got modified by the correction/lag angle; in that way ALL angles reduced towards or increased away from
     *    the input instead of rotating all directions by the same value 
     *  > Consequence 3: the angle correction should be restricted to the actual angle between input and output; otherwise
     *    the angle correction could cause the formerly behind-lagging vector to preceed the other
     */

    Double absX = xVal.NormL2();
    Double absPrevX = prevXVal.NormL2();
    if(absPrevX <= 1e-16){
      return 0;
    }
   
    // 20.10.2020: according to the original source, the lag angle always takes the absolute value of dPhi, i.e.,
    // the absolute change in input-angle; the question is, if the computed correction angle will work correctly for
    // clockwise and counterclockwise rotations
    bool useAbsDPhi = true;
    Double absDPhi = 0;
    
//    // old way of computing the angle dPhi
//    Double dotProd = 0;
//    xVal.Inner(prevXVal,dotProd);
//    Double absProd = absX*absPrevX;
//    
//    Double tmp = 0;
//    if(absProd != 0){
//      tmp = dotProd/absProd;
//      if(tmp > 1.0){
//        tmp = 1.0;
//      } else if(tmp < -1.0){
//        tmp = -1.0;
//      }
//      absDPhi = (std::acos(tmp));
//    } else {
//      // if one of the vectors is zero (absProd = 0), assume alignment
//      absDPhi = 0;
//    }
//
//    std::cout << "-- absDPhi (old comp, in deg): " << absDPhi*180.0/M_PI << std::endl; 
   
    absDPhi = xVal.InnerAngle(prevXVal);
        
    if(useAbsDPhi){
//      std::cout << "Take absolute value of dPhi" << std::endl;
      absDPhi = std::abs(absDPhi);
    } 
//    else {
//      std::cout << "Take dPhi directly, not its absolute value" << std::endl;
//    }
    
    Double absXNormalized = absX/XSaturated_;
    Double gX = std::exp(-lossParam_a_*absXNormalized + lossParam_b_);
    // regarding note from 19.06.2020 above: limit angular correction (=maximal angular lag) to absDPhi
    Double psi = 0;

    // note: 20.6.2020 - be restricting the angluar lag to absDPhi we are again able to invert the mayergoyz
    // model appropriately as the jacobian no longer jumps; however, this restriction also leads to a non-sufficient
    // reduction of rotational losses
    int mode = 4;
    // 0: unrestricted; matching works but inversion fails
    // 1: restricted to absDPhi; inversion works but matching fails; rotloss is not reduced to 0
    // 2: restricted to absDPhi; negative value returned; same issue as 1
    // 3: restricted to 0,2pi; maybe this works > yes it does! it leads to the same rotational losses as mode 0 but
    //    keeps the model invertible with newton, as the jacobian no longer jumps;
    //    > matching to actual loss curves fails however
    // 4: restricted to angle between current input and previous output; best results so far
    // CONCLUSION: 
    //  - angular lag must be allowed to become larger than absXNormalized
    //  - angular lag must be allowed to become larger than 180 degree, i.e., it some of the scalar models
    //      must be allowed to flip their direction!
    //  - ONLY MODE 4 WORKING PROPERLY AT THE MOMEMENT!
    // Further note: if the starting value for Newton inversion is good enough, the problem with the exploding
    //  Jacobian does not appear for mode 0 either! Nevertheless, it is a risk to use mode 0.
    if(mode != 4){
      WARN("Lag correction for extended Mayergoyz model with mode "<<mode<<" requested. However, only mode 4 is working properly (at the moment).")
    }
    
    Double scaling = gX*absXNormalized;
    if(mode == 0){
      psi = scaling*absDPhi;
    } else if(mode == 1){
      if(scaling > 1.0){
        scaling = 1.0;
      }
      psi = scaling*absDPhi;
    } else if(mode == 2){
      if(scaling > 1.0){
        scaling = 1.0;
      }
      psi = -scaling*absDPhi;
    } else if(mode == 3){
      psi = scaling*absDPhi;
      
      if(psi > 2*M_PI){
        psi = 2*M_PI;
      }
    } else if(mode == 4){
      psi = scaling*absDPhi;
      
      if(psi > maxAngle){
        psi = maxAngle;
      }
    }
    return psi;
  }
}
