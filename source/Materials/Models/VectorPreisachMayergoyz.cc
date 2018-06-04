#include "VectorPreisachMayergoyz.hh"

#include <iostream>
#include <fstream>
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Utils/Timer.hh"

namespace CoupledField
{ 
  class Preisach;
  
  DECLARE_LOG(vecpreisachmayergoyz)
  DEFINE_LOG(vecpreisachmayergoyz, "vecpreisachmayergoyz")
  
  VectorPreisachMayergoyz::VectorPreisachMayergoyz(Integer numElem, Vector<Double> xSat, Vector<Double> ySat, 
          Matrix<Double>* preisachWeight, UInt dim, bool isVirgin,
          Vector<Double> anhyst_A, Vector<Double> anhyst_B, Vector<Double> anhyst_C, bool anhystOnly, int clipOutput) : Hysteresis(numElem){
    
    EXCEPTION("Anisotropic case not yet implemented; Find way to compute weightings first!");
    
  }
  
  VectorPreisachMayergoyz::VectorPreisachMayergoyz(Integer numElem, UInt numDirections, Double xSat, Double ySat, 
          Matrix<Double>& preisachWeight, UInt dim, bool isVirgin, 
          Double anhyst_A, Double anhyst_B, Double anhyst_C, bool anhystOnly, int clipOutput) : Hysteresis(numElem){
    
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
    if(dim != 2){
      EXCEPTION("Mayergoyz vector model only implemented for 2d");
      // for the 3d case the computation of the weightings has to be adapted; furtherone
      // one has to specify many more directions as we have to integrate over a whole
      // hemisphere instead of a hemicircle
    } 
    
    numDirections_ = numDirections;
    if(numDirections < 2){
      EXCEPTION("To obtain vector functionality, at least 2 directions are required");
    }
    
    singleDirections_ = new Vector<Double>[numDirections_];
    
    // accroding to "Magnetic Field Analysis of Electric Machines Taking Ferromagnetic Hysteresis into Account" by J. Saitz, p. 38
    // we should use randomly distributed starting angle to make up for the missing symmetry property due to discretization
    Double deltaAngle = M_PI/numDirections;
    Double startingAngle = -deltaAngle/2;
    Double currentAngle;
    for(UInt i = 0; i < numDirections; i++){
      currentAngle = startingAngle + i*deltaAngle;
      
      singleDirections_[i] = Vector<Double>(dim_);
      singleDirections_[i][0] = std::cos(currentAngle);
      singleDirections_[i][1] = std::sin(currentAngle);
    }
    
		// Inversion requires the following steps:
		// 1. decompose input vector Yin into the defined directions
		//			> Yin_k = Yin.Inner(singleDirections_[k])
		// 2. use inversion of scalar models to obtain Xout_k such that
		//			  Yin_k = HystForward(Xout_k) 
		// 3. combine Xout_k correctly to get the output vector Xout
		//			> Problem: Xout = sum_k Xout_k * singleDirections_[k] will not work
		//			> Xout.Inner(singleDirections_[k]) != Xout_k 
		//			> Solution:
		//			> build up Xout such, that Xout.Inner(singleDirections_[k]) = Xout_k 
		//			> to do so, have to solve the system:
		//				sum_j Xout_comp_j singleDirections_[j].Inner(singleDirections_[k]) = Xout_k 
		// > does not work! matrix nearly singular
//		std::cout << "Single directions: " << std::endl;
//		for(UInt i = 0; i < numDirections; i++){
//			std::cout << singleDirections_[i].ToString() << std::endl;
//		}
//		// tmpMatrix[j][k] = singleDirections_[j].Inner(singleDirections_[k])
//		Matrix<Double>tmpMatrix = Matrix<Double>(numDirections,numDirections);
//		Matrix<Double>tmpMatrix2 = Matrix<Double>(numDirections,numDirections);
//		for(UInt j = 0; j < numDirections; j++){
//			for(UInt k = 0; k < numDirections; k++){
//				tmpMatrix[j][k] = singleDirections_[j].Inner(singleDirections_[k]);
//			}
//		}
//		std::cout << "tmpmatrix = " << std::endl;
//		std::cout << tmpMatrix.ToString() << std::endl;
//		matrixForCoefComputation_ = Matrix<Double>(numDirections,numDirections);
//		
//		tmpMatrix.Invert(matrixForCoefComputation_);
////		matrixForCoefComputation_ = tmpMatrix;
////		matrixForCoefComputation_.Invert_Lapack();
//		std::cout << "matrixForCoefComputation_ = " << std::endl;
//		std::cout << matrixForCoefComputation_.ToString() << std::endl;
//		tmpMatrix.Mult(matrixForCoefComputation_,tmpMatrix2);
//		std::cout << "tmpmatrix*matrixForCoefComputation_ = " << std::endl;
//		std::cout << tmpMatrix2.ToString() << std::endl;
//		matrixForCoefComputation_ > inverse of tmpMatrix
		
    clipOutput_ = clipOutput;
    isIsotropic_ = true;
    
    if(isIsotropic_){
      // add anhysteretic part directly to output and set saturation parameter accordingly
      // the single scalar models will not have anhysteretic parts
      XSaturated_ = xSat;
      PSaturated_ = ySat;
      anhyst_A_ = anhyst_A;
      anhyst_B_ = anhyst_B;
      anhyst_C_ = anhyst_C;
      anhystOnly_ = anhystOnly;
      
      /*
       * IMPORTANT REMARK:
       *  > although we use scalar models here, we are not alloewed to directly apply the preisach parameter for
       *     the scalar case (i.e. the weights, the anhyst parameter and so on)
       *  > make sure that the passed parameter are already transformed correctly
       *      > see CoefFunctionHyst
       */
      singlePreisachOperators_ = new Preisach*[numDirections];
      for(UInt i = 0; i < numDirections_; i++){
        singlePreisachOperators_[i] = new Preisach(numElem, xSat, ySat, preisachWeight, isVirgin,0.0, 0.0, 0.0, false);
      }
      
    } else {
      // in each single direction, we need specifiy values for
      // xSat, ySat, preisachWeights, anhyst_A, anhyst_B, anhyst_C
      // > unless we have a way to measure or compute these values, the anisotropic case is not available
      EXCEPTION("Only isotropic case implemented.");
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
//    
//  Vector<Double> VectorPreisachMayergoyz::computeInput_vec(Vector<Double> Yin, Integer idx, Matrix<Double> eps_mu, bool overwriteDirectioncomputeValue){
//    /*
//     * TO FIND OUT: can we actually use the same scalar preisach models (i.e. with the same weights) for forward and backward computation?
//     * > not found so far in literatur
//     */
//    
//    /*
//     * Vectorial output = integral/sum over scalar models
//     */
//    Vector<Double> yCheck = Vector<Double>(dim_);
//    yCheck.Init();
//		
//		return yCheck;
////    Double yCheckScal;
////    
////    Vector<Double> output = Vector<Double>(dim_);
////    output.Init();
////    
////    
////    Vector<Double> currentDir;
////    Vector<Double> tmp = Vector<Double>(dim_);
////    Double scalarInput, scalarOutput, eps_mu_scal;
////    bool overwrite = false; // never overwrite here // overwriteDirection stemming from VectorPreisach to get same function header
////    for(UInt i = 0; i < numDirections_; i++){
////      currentDir = singleDirections_[i];
////      eps_mu.Mult(currentDir,tmp);
////      currentDir.Inner(tmp,eps_mu_scal);
////      
////      currentDir.Inner(Yin,scalarInput);
////      scalarOutput = singlePreisachOperators_[i]->computeInputAndUpdate(scalarInput,eps_mu_scal,idx,overwrite);
////      
////      yCheckScal = singlePreisachOperators_[i]->computeValueAndUpdate(scalarOutput,idx,false);
////      yCheckScal += eps_mu_scal*scalarOutput;
//////      std::cout << "Direction #" << i << std::endl;
//////      std::cout << "yIn = " << scalarInput << std::endl;
//////      std::cout << "yRet = " << yCheckScal << std::endl;
//////      std::cout << "xRet = " << scalarOutput << std::endl;
////      
////      output.Add(scalarOutput,currentDir);
////      yCheck.Add(yCheckScal,currentDir);
////    }
//////    std::cout << "output = " << output.ToString() << std::endl;
////    Vector<Double> tmpTest = Vector<Double>(dim_);
////    tmpTest.Init();
////    for(UInt i = 0; i < numDirections_; i++){
//////      std::cout << "xRet, extracted, dir " << i << " = " << singleDirections_[i].Inner(output) << std::endl;
////      tmpTest.Add(singleDirections_[i].Inner(output),singleDirections_[i]);
////    }
//////    std::cout << "tmpTest = " << tmpTest.ToString() << std::endl;
////    tmpTest.ScalarMult(2.0/numDirections_);
//////    std::cout << "tmpTest = " << tmpTest.ToString() << std::endl;
////    //Double deltaAngle = M_PI/numDirections_;
////    // for numerical integration we have to multiply by deltaAngle
////    // and to average out correctly over the halfspace we have to multiply by 2.0/Pi
////    // > in total, multiply by 2/numDir
////    output.ScalarMult(2.0/numDirections_);
////    yCheck.ScalarMult(2.0/numDirections_);
////
//////    std::cout << "yInVec = " << Yin.ToString() << std::endl;
//////    std::cout << "yRetVec = " << yCheck.ToString() << std::endl;
//////    std::cout << "output = " << output.ToString() << std::endl;
//////    
////    yCheck = computeValue_vec(output, idx, false, false, false);
//////    std::cout << "yRetVec from computeValue_vec = " << yCheck.ToString() << std::endl;
////    eps_mu.Mult(output,tmp);
////    yCheck.Add(tmp);
////    return output;
//  }
  Vector<Double> VectorPreisachMayergoyz::computeValue_vecMeasure(Vector<Double>& xVal, Integer idx, 
          bool overwrite,bool debugOutput,int& successFlag, Double& time){
    
    Timer* timer = new Timer();
    Double startTime = timer->GetCPUTime();
    timer->Start();
    
    Vector<Double> Yvec = computeValue_vec(xVal, idx, overwrite, debugOutput, successFlag);
    
    timer->Stop();
    Double endTime = timer->GetCPUTime();  
    time = endTime-startTime;
    
    return Yvec;
    
  }  
    
    
  Vector<Double> VectorPreisachMayergoyz::computeValue_vec(Vector<Double>& xVal, Integer idx, 
          bool overwrite,bool debugOutput,int& successFlag){
    
        /*
     * Vectorial output = integral/sum over scalar models
     */
    Vector<Double> output = Vector<Double>(dim_);
    output.Init();
    
    Vector<Double> currentDir;
//    Vector<Double> tmp = Vector<Double>(dim_);
    Double scalarInput, scalarOutput;
    
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
    for(UInt i = 0; i < numDirections_; i++){
      currentDir = singleDirections_[i];
      currentDir.Inner(xVal,scalarInput);
      scalarOutput = singlePreisachOperators_[i]->computeValueAndUpdate(scalarInput,idx,overwrite,successFlagSingle);
      output.Add(scalarOutput,currentDir);
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
		output.ScalarMult(M_PI/numDirections_);
    
    Vector<Double> dirInput = Vector<Double>(dim_);
    dirInput.Init();
    if((xVal.NormL2() != 0)&&(clipOutput_ == 2)){
      dirInput.Add(1.0/xVal.NormL2(),xVal);
    }
        
    if(clipOutput_ == 1){
      // clip amplitude to saturation; works well if input only in 1d but
      // not so well if remanent parts perpendicular to input exist as those
      // will be scaled down, too
      if(output.NormL2() > PSaturated_){
        output.ScalarMult(PSaturated_/output.NormL2());
      }
    } else if(clipOutput_ == 2){
      // > default
      // clip amplitude to saturation, but such that remanent part is not
      // affected; results seem to be more reasonable than unclipped and clipping 1  
      Double projection = output.Inner(dirInput);
      if((abs(output.NormL2()) > PSaturated_)&&(dirInput.NormL2() != 0)){
        output.Add(-projection,dirInput);
        Double normRemaining = output.NormL2();
        output.Add(std::sqrt(PSaturated_*PSaturated_-normRemaining*normRemaining),dirInput);
      }
    } 
    if(isIsotropic_){
      // for isotropic case, add anhystPart directly to output
      // make sure that the scalar models return no anhystPart in this case
      if(xVal.NormL2() != 0){
        Double amplitude = PSaturated_*evalAnhystPart_normalized(xVal.NormL2()/XSaturated_);
        output.Add(amplitude,dirInput); 
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
