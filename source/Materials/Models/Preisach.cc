#include "Preisach.hh"

#include <iostream>
#include <fstream>
#include "DataInOut/Logging/LogConfigurator.hh"


namespace CoupledField
{ 
  class Preisach;
  
  DECLARE_LOG(scalpreisach)
  DEFINE_LOG(scalpreisach, "scalpreisach")
  DECLARE_LOG(scalpreisachInversion)
  DEFINE_LOG(scalpreisachInversion, "scalpreisachInversion")
  
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
    Double startingAngle = 0.0;
    Double deltaAngle = M_PI/numDirections;
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
      YSaturated_ = ySat;
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
  
  Vector<Double> VectorPreisachMayergoyz::computeValue_vec(Vector<Double>& xVal, Integer idx, bool overwrite,bool overwriteDirection,bool debugOutput){
        /*
     * Vectorial output = integral/sum over scalar models
     */
    Vector<Double> output = Vector<Double>(dim_);
    output.Init();
    
    Vector<Double> currentDir;
//    Vector<Double> tmp = Vector<Double>(dim_);
    Double scalarInput, scalarOutput;
    
		// idea: cap xVal to Xsat BEFORE feeding it into the scalar models
		// background: inside each scalar model, we cap the scalar input to Xsat
		// here, each model gets (depending on the direction) only the projection of Xval
		// i.e. if Xval is above saturation, the projections will be capped which
		// leads to some inpute being capped whereas others are not; by this, we
		// grow above the point of saturation as due to the projection there will be
		// some scalar models which are first capped if xval is way above saturation
		Vector<Double> xValCapped = Vector<Double>(dim_);
		xValCapped = xVal;
		
		if(xVal.NormL2() > XSaturated_){
			xValCapped.ScalarMult(XSaturated_/xVal.NormL2());
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
    
    
    
    for(UInt i = 0; i < numDirections_; i++){
      currentDir = singleDirections_[i];
      currentDir.Inner(xVal,scalarInput);
			
      scalarOutput = singlePreisachOperators_[i]->computeValueAndUpdate(scalarInput,idx,overwrite);
      output.Add(scalarOutput,currentDir);
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
    if(xVal.NormL2() != 0){
      dirInput.Add(1.0/xVal.NormL2(),xVal);
    }
    
    if(clipOutput_ == 1){
      // clip amplitude to saturation; works well if input only in 1d but
      // not so well if remanent parts perpendicular to input exist as those
      // will be scaled down, too
      if(output.NormL2() > YSaturated_){
        output.ScalarMult(YSaturated_/output.NormL2());
      }
    } else if(clipOutput_ == 2){
      // > default
      // clip amplitude to saturation, but such that remanent part is not
      // affected; results seem to be more reasonable than unclipped and clipping 1
      Double projection = output.Inner(dirInput);
      if(abs(output.NormL2()) > YSaturated_){
        output.Add(-projection,dirInput);
        Double normRemaining = output.NormL2();
        output.Add(std::sqrt(YSaturated_*YSaturated_-normRemaining*normRemaining),dirInput);
      }
    }

    if(isIsotropic_){
      // for isotropic case, add anhystPart directly to output
      // make sure that the scalar models return no anhystPart in this case
      if(xVal.NormL2() != 0){
        Double amplitude = YSaturated_*evalAnhystPart_normalized(xVal.NormL2()/XSaturated_);
        output.Add(amplitude,dirInput); 
      } 
    }
    
    if(overwrite == true){
      /*
       * store to arrays > as in VectorPreisachv10
       * include anhyst part!
       */
      prevXVal_[idx] = xVal;
      prevHVal_[idx] = output;
    }
    
    return output;
  }
    
  ExtendedPreisach::ExtendedPreisach(Integer numElem, Double xSat, Double ySat, 
          Matrix<Double>& preisachWeight, Double rotationalResistance , Double angularDistance, UInt dim, bool isVirgin, 
          Double anhyst_A, Double anhyst_B, Double anhyst_C, bool anhystOnly):
  Preisach(numElem, xSat, ySat, preisachWeight, isVirgin, anhyst_A, anhyst_B, anhyst_C, anhystOnly){
    
    if(angularDistance != 0){
      WARN("angularDistance not yet used in ExtendedPreisach");
    }
    
    rotResistance_ = rotationalResistance;
    dim_ = dim;
    
    rotationStates_ = new std::map<Double,Vector<Double> >[numElem];
    currentDirection_ = new Vector<Double>[numElem];
  }
  
  ExtendedPreisach::~ExtendedPreisach(){
    delete [] rotationStates_;
    delete [] currentDirection_;
  }
  
  void ExtendedPreisach::UpdateRotationState(Vector<Double> flux_in, Matrix<Double> eps_mu, UInt idx){
    
    /*
     * 1. Decompose flux_in into amplitude and direction
     */
    Double ampl = flux_in.NormL2();
    if(ampl <= 0){
      // nothing to update
      return;
    }
    
    Vector<Double> dir = flux_in;
    dir.ScalarDiv(ampl);
    
    /*
     * 2. Compute actual setting value as for Vector Modell 2015er edition
     */
    // major difference: we use the flux quantity; full setting shall be achived
    // if material goes in saturation, i.e. we have to clip towards YSaturation_
    // however. YSaturation is the saturation value for P not for the flux quantity
    Vector<Double> satInCurDir = Vector<Double>(dim_);
    eps_mu.Mult(dir,satInCurDir);
    satInCurDir.ScalarMult(XSaturated_);
    satInCurDir.Add(YSaturated_,dir);
    // satInCurDir = (YSaturated_*Identity + XSaturated*eps_mu)*dir
    Double satValue = satInCurDir.NormL2();
    
    if(ampl >= satValue){
      ampl = 1.0;
    } else {
      ampl = ampl/satValue;
    }
    
    Double X_thres = ampl*rotResistance_;
    if(X_thres > 1.0){
      X_thres = 1.0;
    } else if (X_thres <= 0){
      return;
    }
    
    /*
     * 3. Update list using X_thres
     */
    // special case: list empty
    if(rotationStates_[idx].empty()){
      // directly insert X_thres, dir
      rotationStates_[idx][X_thres] = dir;
    } else {
      std::map<Double, Vector<Double> >::iterator insertPos;
      std::pair<std::map<Double, Vector<Double> >::iterator,bool> retVal; 
      /*
       * Advantage from taking a map:
       *  > map is sorted automatically
       *  > insert function will find the correct position (according to key)
       *      and will return its position
       *  > using this position we can easily wipe out rest according to wiping
       *      out rules (which are much easier for the rotation states as the 
       *      only rule is: remove all entries with key < X_thres, keep rest
       */
      retVal = rotationStates_[idx].insert(std::pair<Double, Vector<Double> >(X_thres, dir));
      insertPos = retVal.first;
      
      // if X_thres was inserted, insertPos point to the new entry
      // if X_thres was not inserted (due to its value already in the map), insertPos points to the 
      //  this prevously inserted element
      
      std::map<Double, Vector<Double> >::iterator startPos = rotationStates_[idx].begin();
      rotationStates_[idx].erase(startPos,insertPos); // insertPos is kept
      
    }
  }
  
  void ExtendedPreisach::EvaluateRotationState(UInt idx){
    // evaluate rotationStates_ and store in currentDirection_
    // iterate through map from back to front
    Vector<Double> dirVector = Vector<Double>(dim_);
    dirVector.Init();
    Double curWeight;
    
    if(!rotationStates_[idx].empty()){

      Double curVal;
      Double nextVal;
    
      std::map<Double, Vector<Double> >::reverse_iterator curPos;
      for(curPos = rotationStates_[idx].rbegin(); curPos != rotationStates_[idx].rend(); curPos++){
        if(curPos++ == rotationStates_[idx].rend()){
          // currently we are in the last entry 
          nextVal = 0.0;
        } else {
          nextVal = curPos->first;
        }
        curPos--;
        curVal = curPos->first;
        
        // each rotation state is weighted with the actual area that is assigned
        // to it in the Preisach-like plane
        // as in the SutorModel from 2015 this plane consists of triangles like
        /*
         * Rotation states form flipped L-shapes in the whole Preisach plane
         *
         *            S_U           alpha           T_U
         *      __ __ __ __ __ __ __ |_  __ __ __ __ __ __
         *     |                     |                  /
         *     |A0                   |                /
         *     |      _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ /_ ex1
         *     |     |   A1          |            /
         *     |     |     _ _ _ _ _ |_ _ _ _ _ /_ ex2
         *     |     |   |    A2     |        /
         *     |     |   |           |      /
         *     |     |   |        _ _| _ _/_ _ exN
         *     |     |   |      |    |  /
         *     |_____|___|______|_AN_|/____________________ beta
         *     |     |   |      |   /
         *     |     |   |      | /
         *     |     |   |      /
         *     |     |   |    /
         *     |     |   |  /
         *     |     |   |/
         *     |     |  /
         *     |     |/
         *     |    /
         *     |_ /
         *           T_L
         */
        // Unlike the actual vector model, these states are not furtehr sub
        // divided so that the area computation is very simple
        // > take outersquare - innersquare, divide by 2 (for triangles), divide
        //  by 2 to value between 0 and 1 (preisach plane goes from -1 to +1)
        curWeight = (curVal*curVal - nextVal*nextVal)/4.0;
        dirVector.Add(curWeight,curPos->second);
      }
    }
    currentDirection_[idx] = dirVector;  
  }
  
  Preisach::Preisach(Integer numElem, Double xSat, Double ySat, 
          Matrix<Double>& preisachWeight, bool isVirgin) 
  : Hysteresis(numElem,xSat,ySat,0.0, 0.0, 0.0,false){
    Preisach(numElem, xSat, ySat, preisachWeight, isVirgin, 0.0, 0.0, 0.0,false);
  }
  
  Preisach::Preisach(Integer numElem, Double xSat, Double ySat, 
          Matrix<Double>& preisachWeight, bool isVirgin, Double anhyst_A, Double anhyst_B, Double anhyst_C,bool anhystOnly) 
  : Hysteresis(numElem,xSat,ySat,anhyst_A,anhyst_B,anhyst_C,anhystOnly)
  {
    LOG_TRACE(scalpreisach) << "ScalarPreisach: Using Everett function";
    //tol_ = 1e-5;
    tol_ = 1e-15;
    
    if (xSat > 0 ) {
      XSaturated_  = xSat;
    }
    else {
      XSaturated_  = 1.0;
    }
       
    YSaturated_  = ySat;
    isVirgin_    = isVirgin;
    
    preisachWeights_ = preisachWeight;
    preisachSum_.Resize(numElem);
    preisachSum_.Init();
    StringLength_.Resize(numElem);
   
//    anhyst_A_ = anhyst_A;
//    anhyst_B_ = anhyst_B;
//    anhyst_C_ = anhyst_C;
    
    strings_     = new Vector<Double>[numElem];
    helpStrings_ = new Vector<Double>[numElem];
    minmaxtype_ = new Vector<Integer>[numElem];
    evaluatedEverettPixel_ = new Vector<Double>[numElem];
    maxStringLength_ = 100;

    for (Integer el=0; el<numElem; el++) {
      strings_[el].Resize(maxStringLength_);
			minmaxtype_[el].Resize(maxStringLength_);
			evaluatedEverettPixel_[el].Resize(maxStringLength_);
			
      helpStrings_[el].Resize(maxStringLength_+1);
      
      StringLength_[el] = 1;
      for ( UInt i=0; i<maxStringLength_; i++) {
        strings_[el][i] = 0.0;
				minmaxtype_[el][i] = 0; // neither min nor max
				evaluatedEverettPixel_[el][i] = 0.0;
      }
    }
    
    // for inversion
    previousXval_ = Vector<Double>(numElem);
    previousXval_.Init();
    
    previousPval_ = Vector<Double>(numElem);
    previousPval_.Init();
  }
  
  Preisach::~Preisach()
  {
    delete [] strings_;
    delete [] helpStrings_;
		delete [] evaluatedEverettPixel_;
		delete [] minmaxtype_;
  }
  
  Double Preisach::getValue( Integer idx ) 
  {
    return ( preisachSum_[idx]*YSaturated_ );
  }
  
    
  Double Preisach::bisect(Double dY,Double xMin,Double xMax, Double xFixed, Double eps_mu, Double tol){
		LOG_DBG(scalpreisachInversion) << "perform bisection for dY: " << dY;
//    std::cout << "perform bisection" << std::endl;
//    std::cout << "dY: " << dY << std::endl;
//    std::cout << "xMin: " << xMin << std::endl;
//    std::cout << "xMax: " << xMax << std::endl;
//    std::cout << "xFixed: " << xFixed << std::endl;
    
    /*
     *  Try to find xOpt, such that
     * 
     *  abs(dY - dP - dX_to_dY*dX) = min
     *  with dP = everettPixel(xFixed,xOpt)
     *  and  dX = xOpt-xFixed
     * 
     * if xFixed < -1, find xOpt, such that
     * 
     *  abs(dY) = abs(everettPixel(-xOpt,xOpt)) = min
     * 
     * > use bisection to find xOpt between xMin and xMax
     */ 
    UInt maxIter = 50;
    Double xMiddle;
    Double diffMin, diffMiddle;
    Double everettMin, everettMiddle;
    Double dX_to_dY = eps_mu*XSaturated_/YSaturated_;
    Double dX = 0.0;
    Double dAnhyst = 0.0;
    
    if(xFixed >= -1){
      LOG_DBG(scalpreisachInversion) << "perform bisection with fixed edge";
      //std::cout << "bisection with fixed edge" << std::endl;
      for(UInt i = 0; i<maxIter; i++){
        xMiddle = 0.5*(xMax+xMin);
        dX = xMiddle - xFixed;
        dAnhyst = evalAnhystPart_normalized(xMiddle) - evalAnhystPart_normalized(xFixed);
        // important: here we need a factor of 2.0 as we have an update
        // to the previous state (i.e. we have to revert from +1 to -1 or
        // vice versa)
        everettMiddle = 2.0*everettPixel(xFixed,xMiddle);
        // dY = dP + dX_to_dY * dX
        diffMiddle = dY - everettMiddle - dX*dX_to_dY - dAnhyst;
        
        if(abs(diffMiddle) < tol){
					LOG_DBG(scalpreisachInversion) << "Remaining diff (normalized): " << abs(diffMiddle);
          LOG_DBG(scalpreisachInversion) << "Remaining diff: " << YSaturated_*abs(diffMiddle);
          return xMiddle;
        }
        
        dX = xMin - xFixed;
        dAnhyst = evalAnhystPart_normalized(xMin) - evalAnhystPart_normalized(xFixed);
        
        // syntax for everettPixel: 1. parameter: older/previous entry, 2. parameter new/next entry
        everettMin = 2.0*everettPixel(xFixed,xMin);
        
        diffMin = dY - everettMin - dX*dX_to_dY - dAnhyst;
        
//				std::cout << "xMiddle: " << xMiddle << std::endl;
//        std::cout << "everettMiddle: " << everettMiddle << std::endl;
//        std::cout << "everettMin: " << everettMin << std::endl;
//        std::cout << "diffMiddle: " << diffMiddle << std::endl;
//        std::cout << "diffMin: " << diffMin << std::endl;
        if( ((diffMiddle < 0)&&(diffMin < 0)) || ((diffMiddle > 0)&&(diffMin > 0)) ){
          // search between middle and max
          xMin = xMiddle;
        } else {
          // search between min and middle
          xMax = xMiddle;
        }
      }
      LOG_DBG(scalpreisachInversion) << "Maxnumber of iteations reached ( case: xFixed >= -1 )";
			LOG_DBG(scalpreisachInversion) << "Remaining diff (normalized): " << abs(diffMiddle);
      LOG_DBG(scalpreisachInversion) << "Remaining diff: " << YSaturated_*abs(diffMiddle);
			return xMiddle;
    } else {
      LOG_DBG(scalpreisachInversion) << "perform bisection without fixed edge";
      // for the open edge case, we have to be careful as the everettPixel
      // inside an empty Preisach plane will always be positive for x in 0,1;
      // to get a negative dY, we have to search in -1,0 instead
      if(dY < 0){
        xMax = -xMax;
        xMin = -xMin;
      }
      for(UInt i = 0; i<maxIter; i++){
        //std::cout << "bisection with open edge" << std::endl;
        xMiddle = 0.5*(xMax+xMin);
        dX = xMiddle; // here we have dX towards 0
        dAnhyst = evalAnhystPart_normalized(xMiddle);
        // here we do not need a factor of 2.0 as we do not have to
        // revert a previous state
        everettMiddle = everettPixel(-xMiddle,xMiddle);

        diffMiddle = dY - everettMiddle - dX*dX_to_dY - dAnhyst;
        
        if(abs(diffMiddle) < tol){
          LOG_DBG(scalpreisachInversion) << "Remaining diff (normalized): " << abs(diffMiddle);
          LOG_DBG(scalpreisachInversion) << "Remaining diff: " << YSaturated_*abs(diffMiddle);
          return xMiddle;
        }
        
        // note that everettPixel always returns a positive value (the area)
        // and a sign that can be used to find out if this area is to be subtracted
        // or to be added
        dX = xMin; // here we have dX towards 0
        dAnhyst = evalAnhystPart_normalized(xMin);
        everettMin = everettPixel(-xMin,xMin);
        //everettMax = everettPixel(-xMax,xMax,sign);
				
        diffMin = dY - everettMin - dX*dX_to_dY - dAnhyst;
        //diffMin = dYabs - abs(everettMin + dX*dX_to_dY);
        //diffMax = dYabs - everettMax;
				
//        std::cout << "xMiddle: " << xMiddle << std::endl;
//        std::cout << "everettMiddle: " << everettMiddle << std::endl;
//        std::cout << "everettMin: " << everettMin << std::endl;
//        std::cout << "diffMiddle: " << diffMiddle << std::endl;
//        std::cout << "diffMin: " << diffMin << std::endl;
        if( ((diffMiddle < 0)&&(diffMin < 0)) || ((diffMiddle > 0)&&(diffMin > 0)) ){
          // search between middle and max
          xMin = xMiddle;
        } else {
          // search between min and middle
          xMax = xMiddle;
        }
      }
      LOG_DBG(scalpreisachInversion) << "Maxnumber of iteations reached ( case: xFixed < -1 )";
		  LOG_DBG(scalpreisachInversion) << "Remaining diff (normalized): " << abs(diffMiddle);
      LOG_DBG(scalpreisachInversion) << "Remaining diff: " << YSaturated_*abs(diffMiddle);
			return xMiddle;
    }
  }

  Double Preisach::computeInputAndUpdate(Double Yin, Double eps_mu, Integer idx, bool overwrite){
    /*
     * NEW implementation; this shall replace class SimplePreisachInv that
     * apparently is not working
     * 
     * 
     * Source: "Generalisiertes Preisach-Modell für die Simulation uund Kompensation
     *         piezokeramischer Aktoren" - Dissertation, Felix Wolf, p. 127ff
     */
    
		LOG_TRACE(scalpreisachInversion) << "Compute inverse of Preisach operator for Yin = " << Yin;
    LOG_TRACE(scalpreisachInversion) << "Index = " << idx;
		LOG_DBG(scalpreisachInversion) << "Compute inverse of Preisach operator for Yin_normalized = " << Yin/YSaturated_;
    /*
     * 0. Check if Input drives system into saturation
     *    > make sure to not only compare with YSaturated but also 
     *      include the contribution of XSaturated!
     *    > if eps_mu*XSaturated_ is not considered
     *      we could return XSaturated as Yin > YSaturated but
     *      when recomputing Yin via YSaturated + eps_mu*XSaturated
     *      we can get a significant difference to the actual Yin
     */
    Double Xout;
		Double x1,x2; // search interval for later bisection
		x1 = 0.0;
		x2 = 0.0;
    Double tol = 1e-10;
    // dX_to_dY is the factor which is needed to transfer a normalized dX to 
    // a normalized dY
    // dY_norm = dY/YSaturated = eps_mu*dX/YSaturated = eps_mu*dX*XSaturated/XSaturated/YSaturated
    //         = eps_mu*XSaturated/YSaturated * dX_norm
    Double dX_to_dY = eps_mu*XSaturated_/YSaturated_;
    UInt invcase = 0;
		UInt subcase = 0;
    //std::cout << "Inversion" << std::endl;
    //std::cout << "Yin: " << Yin << std::endl;
    //std::cout << "YSaturated: " << YSaturated_ << std::endl;
    //std::cout << "eps_mu: " << eps_mu << std::endl;
    // anhysteretic parameter are meant to be combined with saturated X > mult by 1.0 here
    Double anhystPartPosSat = YSaturated_*evalAnhystPart_normalized(1.0);
    Double anhystPartNegSat = YSaturated_*evalAnhystPart_normalized(-1.0);
    
    if(anhystOnly_ == true){
      // actual hyst operator is 0; we just solve for the anhyst part via bisection
      Double Xup, Xdown, Poffset;
      Xup = Yin/eps_mu;
      Xdown = 0.0;
      Poffset = 0.0;
      Xout = bisectForAnhyst(Yin, Xdown, Xup, Poffset, eps_mu, tol);
      //Xout = bisectForSaturation(Yin, eps_mu, tol, false);
      invcase = 9;
    } else {
      // we just invert the actual hysteresis operator, i.e we have to subtract all anhyst and reversible parts
      //  before starting the inversion
      if(Yin >= (YSaturated_ + eps_mu*XSaturated_ + anhystPartPosSat) ){
        LOG_DBG(scalpreisachInversion) << "Pos saturation";
        //std::cout << "pos saturation" << std::endl;
        /*
         * System is actually in positive saturation
         * 
         * yIn = Ysaturated_ + eps_mu*xOut + YSaturated_*(anhyst_A_*std::atan(anhyst_B_*xOut/Xsaturated_) + anhyst_C_*xOut/Xsaturated_)
         *     >= YSaturated_ + eps_mu*XSaturated_ + YSaturated_*(anhyst_A_*std::atan(anhyst_B_) + anhyst_C_);
         * 
         * > inversion without anhysteretic part
         * xOut = (yIn - Ysaturated_)/eps_mu
         * 
         * > with anhysteretic part, inversion more complicated!
         */
        if(anhystPartPosSat == 0){
          Xout = (Yin - YSaturated_)/eps_mu;
          invcase = 1;
        } else {
          Double Xup, Xdown, Poffset;
          Xup = (Yin - YSaturated_)/eps_mu;
          Xdown = XSaturated_;
          Poffset = YSaturated_;
          Xout = bisectForAnhyst(Yin, Xdown, Xup, Poffset, eps_mu, tol);
          //Xout = bisectForSaturation(Yin, eps_mu, tol, false);
          invcase = 11;
        }
      } else if (Yin <= (-YSaturated_ - eps_mu*XSaturated_ + anhystPartNegSat) ){
        LOG_DBG(scalpreisachInversion) << "Neg saturation";
        //std::cout << "neg saturation" << std::endl;
        /*
         * System is actually in negative saturation
         * 
         * yIn = -Ysaturated_ + eps_mu*xOut
         * 
         * xOut = (yIn + Ysaturated_)/eps_mu
         */
        if(anhystPartPosSat == 0){
          Xout = (Yin + YSaturated_)/eps_mu;
          invcase = 2;
        } else {
          Double Xup, Xdown, Poffset;
          Xup = (Yin + YSaturated_)/eps_mu;
          Xdown = -XSaturated_;
          Poffset = -YSaturated_;
          Xout = bisectForAnhyst(Yin, Xdown, Xup, Poffset, eps_mu, tol);
          //Xout = bisectForSaturation(Yin, eps_mu, tol, true);
          invcase = 21;
        }
      } else {
        invcase = 3;
        /*
         * 1. Compute difference between requested yVal and previously computed value
         */
        // factor for relating a change in x to a change in y for
        // the non-hysteretic part
        //  > dY = dP + dX_to_dY*dX
        // for dY, dP and dX all normalized

        // TODO:; should be normalized AFTER subtracting x
        // do not cap to +/-1 here as Y is not the polarization (that is actually capped!)
        // EDIT: 26.3.2018
        //  we have to cap, but not to +/-1 (that would be for polarization) but
        //  to +/- (1 + dX_to_dY*Xsaturated)
        // > reason: in the following we compute dY as the difference between oldthat can
        //           be represented by the Preisach plane is
        //           2*YSaturated_ + 2*eps_mu*XSaturated_ 
        //           (considering the reversible part);
        //           If we start in saturation, however, we have a previous |Y| > YSaturated.
        //           This leads to a dY that can be larger than 2*YSaturated_ + 2*eps_mu*XSaturated_.
        //           Even if this is not the case, the later finesearch will try to compensate this
        //           dY by adapting the Preisach plane. This works quite well actually, but is wrong
        //           apparently as the finesearch computes dY from YSaturated+eps_mu*XSaturated and not from the actual
        //           Yold.
        //            > restrict Yold to (YSaturated+eps_mu*XSaturated)/YSaturated
        //          

        // OLD without capping > works well unless we come from a saturated state
  //      Double Y_normalized = Yin/YSaturated_;
  //      Double P_old_normalized = previousPval_[idx];    
  //      Double Y_old_normalized = P_old_normalized+dX_to_dY*previousXval_[idx];
  //      LOG_DBG(scalpreisach) << "yOld: " << Y_old_normalized*YSaturated_;
  //      Double dY = Y_normalized - Y_old_normalized;

        // NEW with capping
        // NEW 25.4.2018 with anhyst part
        Double Y_normalized = Yin/YSaturated_;
        Double P_old_normalized = previousPval_[idx];    
        // IMPORTANT NOTE: previousXval_[idx] is normalized but NOT clipped > exactly what we need!
        // NOTE: previousPval_ contains already the anhysteretic part!
        Double Y_old_normalized = P_old_normalized + dX_to_dY*previousXval_[idx];// + evalAnhystPart_normalized(previousXval_[idx]);
        LOG_DBG(scalpreisachInversion) << "Y_old: " << Y_old_normalized*YSaturated_;
        LOG_DBG(scalpreisachInversion) << "Y_old_normalized: " << Y_old_normalized;
        if(Y_old_normalized > (1 + dX_to_dY + evalAnhystPart_normalized(1.0))){
          Y_old_normalized = 1 + dX_to_dY + evalAnhystPart_normalized(1.0);
        } else if (Y_old_normalized < -(1 + dX_to_dY)+evalAnhystPart_normalized(-1.0)){
          Y_old_normalized = -1 - dX_to_dY + evalAnhystPart_normalized(-1.0);
        }
        LOG_DBG(scalpreisachInversion) << "Y_old_normalized (clipped): " << Y_old_normalized;

        Double dY = Y_normalized - Y_old_normalized;


  //     dY = dY - dX_to_dY*previousXval_[idx];
        //std::cout << "Y_normalized - P_old_normalized - dX_to_dY*previousXval_[idx]: " << dY << std::endl;
        LOG_TRACE(scalpreisachInversion) << "Starting value for dY: " << dY*YSaturated_;
        LOG_DBG(scalpreisachInversion) << "Starting value for dY (normalized): " << dY;
        Integer minmaxcur = 0.0;
        if(dY > 0){
          // larger value than before > input has to be larger than
          // previous one (this is only true if all Preisach weights >= 0!)
          // > retrieved input will be a maximum
          minmaxcur = 1.0;
        } else if(dY < 0){
          minmaxcur = -1.0;
        }
        //std::cout << "minmaxcur: " << minmaxcur << std::endl;
        /*
         * 2. Check if difference to previous value is relevant or not
         */
        if(abs(dY) < 1e-16){
          // difference is negligible
          LOG_TRACE(scalpreisachInversion) << "take previous xvalue" << std::endl;
          invcase = 31;
          // attention: previousXval_ is normalized by XSaturated_
          Xout = previousXval_[idx];
        } else{
          invcase = 32;
          Vector<Double> &stringEl = strings_[idx];
          Vector<int> &minmaxEl = minmaxtype_[idx];
          Vector<Double> &everettEl = evaluatedEverettPixel_[idx];

          x1 = -1.0; // lower bound for max, left bound for min
          x2 = 1.0; // upper bound for max, right bound for min

          Double xfix; 
          UInt actLength = StringLength_[idx];

          /*
           * first step: list has to be cut down such that it ends with
           *						 the opposite minmaxtype (as it is easier to add
           *						 a max after a min than adding a max after a max)
           */
          if( (actLength >= 1)&&(minmaxcur == minmaxEl[actLength-1]) ){
            LOG_DBG(scalpreisachInversion) << "remove last entry of list: " << stringEl[actLength-1];
            //std::cout << "remove last entry of list: " << stringEl[actLength-1] << std::endl;
            // delete last entry (not actually here but just virtually)
            // the value of the corresponding everett pixel now increases
            // the difference dY of course
            // ATTENTION: if eps_mu*XSaturated_ comes into the range of 1e-3 or larger
            //            we might find an everett pixel that fits dY but the result
            //            will still be wrong as dY actually computes to be
            //              Y = everettPixel(X) + eps_mu*X
            // Conclusion: we have to subtract the influence of eps_mu*X from the
            //             searched area update
            Double dP = everettEl[actLength-1];
            //std::cout << "dP (=everett[actLength-1]): " << dP << std::endl;

            dY += dP;
            //std::cout << "dY+dP: " << dY << std::endl;

            // we search space for dY not dP, so we have to add the contribution of
            // the reversible part
            Double dX = stringEl[actLength-1];
            Double dAnhyst = evalAnhystPart_normalized(stringEl[actLength-1]);

            LOG_DBG(scalpreisachInversion) << "Invert Preisach - Add anhystPart " << evalAnhystPart_normalized(stringEl[actLength-1]);
            LOG_DBG(scalpreisachInversion) << "for X/norm = " << (stringEl[actLength-1]);

            if(actLength >= 2){
              LOG_DBG(scalpreisachInversion) << "remove a second entry: " << stringEl[actLength-2];
              dX -= stringEl[actLength-2];
              dAnhyst -= evalAnhystPart_normalized(stringEl[actLength-2]);
            }
            //std::cout << "dX = " << dX << std::endl;
            dY += dX*dX_to_dY + dAnhyst;
            //std::cout << "dY+dP+dX_to_dY*dX: " << dY << std::endl;

            if(minmaxcur > 0){
              // max> max was removed; lower bound can be increased to removed entry
              x1 = stringEl[actLength-1];
            } else {
              // min> min was removed; right bound can be decreased to removed entry
              x2 = stringEl[actLength-1];
            }
            LOG_DBG(scalpreisachInversion) << "Updated dY: " << dY;
            actLength--;
          }

          /*
           * second step: coarse search > find everettPixel that is large enought
           *							to hold dY
           */
          while(true){
            if(actLength <= 0){
              //std::cout << "list empty" << std::endl;
              // list empty

              // in bisection, search for x such that triangle
              // everett(-x,x) = dY
              // to indicate that we want to search for a triangle with
              // alpha and beta to be equal and variable, set xfix (the third parameter
              // to the bisect function to -2.0
              xfix = -2.0;
              // searching range between 0 and maximal value
              x1 = 0;
              x2 = 1;
              subcase = 0;
              break; // go to fine search
            } else {
              // here we need to check if our dY would lead to a new list entry
              // after the current one or if we have to step further backwards to find
              // a fitting place
              // Idea:
              //  everettEl[actLength-1] corresponds to the increment of the polarization, i.e. dP_old
              //  if abs(dP_new) < abs(dP_old), the new change in polarization (dP_new) could be
              //  achieved by inserting a new extremum to the end of the list
              //  if abs(dP_new) > abs(dP_old), the new input would wipe out older entries as the space is
              //  not large enougth; in this case, we go back by 2 positions (to get the same extremum again)
              //  and check again
              // Problem:
              //  aboves procedure only works for updates dP; however, we have dY = dP + eps*dX given
              //  if eps is negligble (i.e. for electrostatics with eps = eps0) this normally is no issue
              //  in case of magnetics and large values of dX, the contribution eps*dX might be significant, though
              //  > we would have to compute dP_new first but this would require the knowledge of dX_new (which of course is unknonw)
              //  > instead we compute dY_old = dP_old * eps*dX_old and compare dY_new with dY_old
              //    by aboves steps; at this point we can no longer compare the areas directly (i.e. abs(dP_new), abs(dP_old))
              //    this leads to problems like the following:
              //      material in saturation, dP (from 0 to saturation) = 1, eps*dX = 1 > dY_old = 2
              //      material shall now have a negative y-value, e.g. -0.1 > dY_new = -0.1-2 = -2.1
              //      > the value -0.1 is larger than neg. saturation (here -1) and thus should lead to a 
              //        a corresponding x_new that is a minimum after the previous maximum of +1 (=xSat)
              //      > if we check for dY_old and dY_new, we see however, that dY_new does not fit into the space
              //        of dY_old (which obviously is wrong)
              //    > solution idea: take the first everett entry twice in terms of space computation
              //              this seems legit as the full preisach plane can swtich from +1 to -1 and thus
              //              has an effetive space of 2
              //

              // the availableSpace consists of the space for dP, i.e. the size of the
              // everett pixel
              Double availSpace = everettEl[actLength-1];

              Double dX = stringEl[actLength-1];
              Double dAnhyst = evalAnhystPart_normalized(stringEl[actLength-1]);
  //            std::cout << "stringEl[actLength-1]: "<< stringEl[actLength-1] << std::endl;
  //            std::cout << "actLength: "<< actLength << std::endl;
              if(actLength >= 2){
                //std::cout << "stringEl[actLength-2]: "<< stringEl[actLength-2] << std::endl;
                dX -= stringEl[actLength-2];
                dAnhyst -= evalAnhystPart_normalized(stringEl[actLength-2]);
                subcase = 2;
              } else {
                // last entry; everettEl counts twice
                availSpace += everettEl[actLength-1];
                // we also have to take the difference towards the other edge of the triangle, so add last entry another time
                // dX += stringEl[actLength-1];
                // > instead: add dX*dX_to_dY to availspace instead of increasing dX!
                // dX is simply stringEL here
                availSpace += dX_to_dY*dX + dAnhyst;
                subcase = 1;
                //
              }
              availSpace += dX_to_dY*dX + dAnhyst;

  //            std::cout << "everettEl[actLength-1]: " << everettEl[actLength-1] << std::endl;
  //            std::cout << "dX_to_dY*dX: " << dX_to_dY*dX << std::endl;
  //            std::cout << "abs(everettEl[actLength-1] + dX_to_dY*dX): " << abs(everettEl[actLength-1] + dX_to_dY*dX) << std::endl;
  //            std::cout << "abs(dY): " << abs(dY) << std::endl;
              // if( abs(everettEl[actLength-1] + dX_to_dY*dX) >= abs(dY) ){
              if( abs(availSpace) >= abs(dY) ){
                // std::cout << "enough space" << std::endl;
                // > enough space

                // fix edge of the everett triangle
                xfix = stringEl[actLength-1];
                if(actLength >= 2){
                  // reduce searching range (if possible due to previous entries)
                  if(minmaxcur > 0){
                    // entry at position actLength-2 is a maximum again and
                    // marks the upper bound for the everett triangle
                    x2 = stringEl[actLength-2];
                  } else {
                    // entry at position actLength-2 is a minimum again and
                    // marks the left bound for the everett triangle
                    x1 = stringEl[actLength-2];
                  }
                }

                break; // go to fine search
              } else {
                // std::cout << "not enough space" << std::endl;

                // not enough space; decrease list by 2 (it has to end with the
                // same minmaxtype again!) and check again
                // the removed entries have to be added to dY (as the difference
                // to the last computed value increases by the removed entries)
                // and adapt x1,x2 if possible

                // everettEl = dP
                // dY = dP + dX_to_dY*dX
                // we actually try to find a fitting dP here, so we must subtract
                // the influence of dX
                dY += everettEl[actLength-1]+dX_to_dY*dX + dAnhyst;
                actLength--;
                if(actLength >= 1){
                  dX = stringEl[actLength-1];
                  dAnhyst = evalAnhystPart_normalized(stringEl[actLength-1]);
                  if(actLength >= 2){
                    dX -= stringEl[actLength-2];
                    dAnhyst -= evalAnhystPart_normalized(stringEl[actLength-2]);
                  }
                  dY += everettEl[actLength-1] + dX_to_dY*dX + dAnhyst;
                  if(minmaxcur > 0){
                    // max> max was removed; lower bound can be increased to removed entry
                    x1 = stringEl[actLength-1];
                  } else {
                    // min> min was removed; right bound can be decreased to removed entry
                    x2 = stringEl[actLength-1];
                  }
                  LOG_DBG(scalpreisachInversion) << "Updated dY: " << dY;
                  actLength--;
                }
              }
            }	
          } // while loop

          /*
           * thrid step: fine search via bisection (xOut in range -1,1)
           */
          Xout = bisect(dY,x1,x2,xfix,eps_mu,tol);
        } // reuse old value
        
        std::cout << "scale Xout with Saturation" << std::endl;
        Xout *= XSaturated_;
      } // pos/neg/no saturation
      // rescale to -xSat to +xSat
      //LOG_TRACE(scalpreisachInversion) << "Found Xout (normalized): " << Xout;
      LOG_TRACE(scalpreisachInversion) << "Found Xout: " << Xout;
      
    } // anhyst only

		/*
     * final step: if overwrite == true > compute forward step to
     *						 update the list (this was not done yet!)
     * 
     * > has to be done for ALL cases
     */	
    bool debug = true;
    if(overwrite || debug){
      LOG_DBG(scalpreisachInversion) << "overwrite? " << overwrite;   
      Double yRetrieved = computeValueAndUpdate( Xout, idx, overwrite );
      //std::cout << "Xout: " << Xout << std::endl;
      //std::cout << "pRetrieved: " << yRetrieved << std::endl;
			LOG_DBG(scalpreisachInversion) << "Found Xout: " << Xout;
      yRetrieved+=Xout*eps_mu;
      LOG_TRACE(scalpreisachInversion) << "yRequested: " << Yin;
      LOG_TRACE(scalpreisachInversion) << "Found Yout: " << yRetrieved;
      LOG_TRACE(scalpreisachInversion) << "InversionCase: " << invcase;
      LOG_TRACE(scalpreisachInversion) << "SubCase: " << subcase;
//      if(abs(yRetrieved-Yin) > tol){
//        LOG_TRACE(scalpreisachInversion) << "Difference: " << yRetrieved-Yin;
//        
//        Double P_old_normalized = previousPval_[idx];
//        Double Y_old_normalized = P_old_normalized+dX_to_dY*previousXval_[idx];
//        
//        LOG_TRACE(scalpreisachInversion) << "yOld: " << Y_old_normalized*YSaturated_;
//        LOG_TRACE(scalpreisachInversion) << "yRequested-yOld: " << Yin-Y_old_normalized*YSaturated_;
//        LOG_TRACE(scalpreisachInversion) << "yRequested-yOld (normalized): " << Yin-Y_old_normalized;
//        LOG_TRACE(scalpreisachInversion) << "x1, x2, xOut: " << x1 << ", " << x2 << ", " << Xout;
//      }
    } else {
      //store value for the case that we want to reuse it later
      // only needed if computeValueAndUpdate is not exectued with overwrite = true
      //NOTE: we must not execute this line before computeValueAndUpdate as this 
      // function would return without overwriting the hyst memory
      //previousXval_[idx] = Xout/XSaturated_;
    }
      		
		return Xout;
	}
  
  //  
  //  Double Preisach::computeValue(Double& Xin, Integer idx, bool overwrite) 
  //  {
  //
  //    /*
  //     * What is this function used for?
  //     * It does not update the list of minima and maxima but only evaluates
  //     * the everett function using the current list checking three cases:
  //     * a) new input wipes out complete list -> only one Everett pixel needed
  //     * b) new input replaces the currently last entry of list
  //     * c) new input attached to end of list
  //     *
  //     * What about cases in between?
  //     * No memory set.
  //     */
  //
  //    Vector<Double> &stringEl = strings_[idx];
  //    UInt& actLength = StringLength_[idx];
  //
  //    //normalize input
  //    Double newX, Yval;
  //    newX = normalizeInput(Xin);
  //
  //    Yval = 0.0;
  //    if ( abs(abs(newX)+eps_) > abs(stringEl[0]) || actLength == 0 ) {
  //      Yval =  everettPixel( -newX, newX );
  //    }
  //    else {
  //      Yval =  everettPixel(-stringEl[0],stringEl[0]);
  //      if ( abs(abs(newX)+eps_) > abs(stringEl[actLength-1]) ) {
  //        for ( UInt i=0; i<actLength-2; i++ ) {
  //          Yval +=  2.0*everettPixel(stringEl[i],stringEl[i+1]);
  //        }
  //        Yval +=  2.0*everettPixel(stringEl[actLength-2], newX);
  //      }
  //      else {
  //        for ( UInt i=0; i<actLength-1; i++ ) {
  //          Yval +=  2.0*everettPixel(stringEl[i],stringEl[i+1]);
  //        }
  //        Yval +=  2.0*everettPixel(stringEl[actLength-1], newX);
  //      }
  //    }
  //
  //    return ( Yval*YSaturated_ );
  //  }
  
  Double Preisach::computeValueAndUpdate( Double Xin, Integer idx,
          bool overwrite )  
  {
    
    //do the deletion
    Double newY = 0.0;
    if(anhystOnly_ == false){
      // update minmaxlist will add anhyst part 
      newY += updateMinMaxList(Xin, idx, overwrite);
    } else  {
      Double X_norm_unclipped = Xin / XSaturated_;
      LOG_DBG(scalpreisachInversion) << "Eval Preisach - Add anhystPart " << evalAnhystPart_normalized(X_norm_unclipped);
      LOG_DBG(scalpreisachInversion) << "for X/norm = " << (X_norm_unclipped);
      newY += evalAnhystPart_normalized(X_norm_unclipped);
    }
    return ( newY*YSaturated_ );
  }
  
  
  Double Preisach::updateMinMaxList(Double Xin, Integer idx, 
          bool overwrite )
  {
    LOG_TRACE(scalpreisachInversion) << "UpdateMinMaxList - Input: " << Xin;
    //std::cout << "UpdateMinMaxList - Input: " << Xin << std::endl;
    Double newY;
    
    //normalize input
    Double newX = normalizeAndClipInput(Xin);
    
    Vector<Double> &stringEl     = strings_[idx];
    Vector<Double> &helpStringEl = helpStrings_[idx];
		
    UInt& actLength = StringLength_[idx];
    UInt stringLength = actLength;
		
		// determine type of current input
		// only relevant if overwrite is true
    
    // NOTE: previousXval_[idx] = Xin/XSaturated_;
    // > previousXval_ is unclipped
    // > compare with clipped version of previousXval (better for saturation
    //    as > sat and >> sat will reuse max value
//    Double oldX = normalizeAndClipInput(previousXval_[idx]*XSaturated_);
//		Double diff = newX - oldX;
    
    
    Double diff = Xin-previousXval_[idx]*XSaturated_;
		int minmaxcur = 0;
		if(diff > 0){
			// new input is larger than last one > leads to a maximum
			minmaxcur = 1;
		} else if(diff < 0){
			// new input is smaller than last one > leads to a minimum
			minmaxcur = -1;
		} else if(diff == 0){
      // reuse old value but set previousXval anew
      // reason: above we compare the clipped values, Xin might be different from prevX
      LOG_TRACE(scalpreisachInversion) << "Reuse: " << preisachSum_[idx];
      previousXval_[idx] = Xin/XSaturated_;
			return preisachSum_[idx];
		}
    
    //    std::cout << "Element: " << idx << std::endl;
    //
    //    std::cout << "Input Xin: " << Xin << std::endl;
    
    //    std::cout << "Print out entries of min/max list before update" << std::endl;
    //    for ( UInt i=0; i<stringLength; i++ ) {
    //      std::cout << "index " << i << ": " << stringEl[i] << std::endl;
    //    }
    //    std::cout << "#############" << std::endl;
    
    //    std::cout << "Print out entries of helper min/max list before update" << std::endl;
    //    for ( UInt i=0; i<stringLength; i++ ) {
    //      std::cout << "index " << i << ": " << helpStringEl[i] << std::endl;
    //    }
    //    std::cout << "#############" << std::endl;
    
    if ( abs(newX) > abs(abs(stringEl[0]) - tol_) || stringLength == 0 ) {
      stringLength = 1;
      if ( overwrite ){
				stringEl[0]  = newX;
        // at the beginning we check if input became larger or smaller and
        // according to that, we set minmaxcur to -1 (for min) or +1 (for max)
        // this type of extremum is correct if a new entry is added to the list
        // but might be wrong if it wipes out the list (e.g. if we stay in neg
        // saturation)
        // > reset minmaxcur here
        if(newX < 0){
          minmaxcur = -1;
        } else {
          minmaxcur = 1;
        }
				minmaxtype_[idx][0] = minmaxcur;
			} else {
        helpStringEl[0] = newX;
			}
    }
    else { 
      if ( abs( newX - stringEl[stringLength-1] ) > tol_ ) {
        helpStringEl[0] = -stringEl[0];
        for ( UInt i=1; i<=stringLength; i++ ) {
          helpStringEl[i] = stringEl[i-1];
        }
        
        UInt k = 0;
        
        Double a = helpStringEl[stringLength-1];
        Double b = helpStringEl[stringLength];
        
        //      std::cout << "a=" << a << "  b=" << b << std::endl;
        
        while ( ( k<stringLength-1) && 
                ( ( newX<=std::min(a,b) ) || ( newX>=std::max(a,b) ) ) ) {
          k = k + 1;
          a = helpStringEl[stringLength-k-1];
          b = helpStringEl[stringLength-k];
        }
        //        std::cout << "Old string length: " << stringLength << std::endl;
        stringLength = stringLength - k + 1;
        //        std::cout << "New string length: " << stringLength << std::endl;
        
        if (overwrite ) {
          //check, if capacity of string-arrays is too less
          if ( stringLength > maxStringLength_ ) {
            //resize the string-arrays
            maxStringLength_ += (UInt) round( (Double)maxStringLength_ / 2.0 );
            stringEl.Resize(maxStringLength_);
            
            //store the resulting strings
            for ( UInt i=0; i<stringLength-1; i++ ) {
              stringEl[i] = helpStringEl[i+1];
            }
            
            // resize help-String-array 
            helpStringEl.Resize(maxStringLength_+1);
          }
          else {
            //store the resulting strings
            for ( UInt i=0; i<stringLength-1; i++ ) {
              stringEl[i] = helpStringEl[i+1];
            }
          }
        }
        
        if ( overwrite ) {
          //store the new input
          stringEl[stringLength-1] = newX;
					// we only have to store the minmax type of the newly added entry
					// the others are already set (or will be overwritten)
					minmaxtype_[idx][stringLength-1] = minmaxcur;
        }
        else {
          
          //          std::cout << "Print out entries of helper min/max list before resorting" << std::endl;
          //
          //           for ( UInt i=0; i<stringLength+1; i++ ) {
          //             std::cout << "index " << i << ": " << helpStringEl[i] << std::endl;
          //           }
          //           std::cout << "#############" << std::endl;
          
          //correct storage in helpString
          for ( UInt i=0; i<stringLength-1; i++ ) {
            helpStringEl[i] = helpStringEl[i+1];
          }
          
          //          std::cout << "Print out entries of helper min/max list after resorting" << std::endl;
          //
          //           for ( UInt i=0; i<stringLength+1; i++ ) {
          //             std::cout << "index " << i << ": " << helpStringEl[i] << std::endl;
          //           }
          //           std::cout << "#############" << std::endl;
          
          helpStringEl[stringLength-1] = newX; 
        }
      }
      else {
        /*
         * Note: We check if a new input is large/small enough to modify the
         * actual storage stringEl;
         * if this is the case, we update
         * either stringEl (if update = true) or helpStringEl (if update = false)
         * and everything works just fine
         *
         * if this is not the case, we do neither change stringEl nor helpStringEl
         * if update = true, this is no problem, as we evaluate stringEl which is
         * up to date (otherwise an update would have been needed)
         * if update = false, we evaluate helpStringEl. This can be a serious problem
         * if helpStringEl is not stringEL
         * -> if no update of list was performed, copy stringEl (which is up to date)
         * to helpStringEl
         */
        //  std::cout << "No update needed; copy current stringEl to helpStringEl" << std::endl;
        for ( UInt i=0; i<=stringLength-1; i++ ) {
          helpStringEl[i] = stringEl[i];
        }
      }
    }
    
    if ( overwrite ) {
      actLength = stringLength;
      
      //      std::cout << "Overwrite = true" << std::endl;
      //      std::cout << "Print out entries of min/max list" << std::endl;
      //      for ( UInt i=0; i<stringLength; i++ ) {
      //        std::cout << "index " << i << ": " << stringEl[i] << std::endl;
      //      }
      //      std::cout << "#############" << std::endl;
      
      //      std::cout << "Updated non-temporary min/max list" << std::endl;
      //      std::cout << "Print out entries of min/max list" << std::endl;
      //      for ( UInt i=0; i<actLength; i++ ) {
      //        std::cout << "index " << i << ": " << stringEl[i] << std::endl;
      //      }
      //      std::cout << "#############" << std::endl;
      
      //compute preisach-sum
      Double pixelToAdd = everettPixel(-stringEl[0],stringEl[0]);
			evaluatedEverettPixel_[idx][0] = pixelToAdd;
      preisachSum_[idx] = pixelToAdd;
      
      for ( UInt i=0; i<actLength-1; i++ ) {
        pixelToAdd = 2.0*everettPixel(stringEl[i],stringEl[i+1]);
				evaluatedEverettPixel_[idx][i+1] = pixelToAdd;
        preisachSum_[idx] += pixelToAdd;   
      }
      newY = preisachSum_[idx]; 
			
      // add optional anhysteretic part
      // > see e.g. "A preisach-based hysteresis model for magnetic and ferroelectric hysteresis" - A. Sutor 2010
      Double X_norm_unclipped = Xin / XSaturated_;
      LOG_DBG(scalpreisachInversion) << "Eval Preisach - Add anhystPart " << evalAnhystPart_normalized(X_norm_unclipped);
      LOG_DBG(scalpreisachInversion) << "for X/norm = " << (X_norm_unclipped);
      newY += evalAnhystPart_normalized(X_norm_unclipped);
      preisachSum_[idx] += evalAnhystPart_normalized(X_norm_unclipped);
      //newY += anhyst_A_*std::atan(anhyst_B_*X_norm_unclipped) + anhyst_C_*X_norm_unclipped;
      
			// store values for next evaluation
			// ONLY in case of overwrite
      // store the normalized values but Xin has to be unclipped!
      // (newX will be clipped to 1 if X > XSaturated)
			previousXval_[idx] = Xin/XSaturated_;
			previousPval_[idx] = newY;
			
    }
    else {
      /*
       * shouldn't we start at idx=1 here as helpStringEl[0] = -stringEl[0], helpStringEl[1] = stringEl[0]
       * by this we start with everettPixel(+stringEl[0],-stringEl[0])
       * then add 2*everettPixel(-stringEl[0],+stringEl[0])
       * ???
       * -> No, see line 181 -> elements get shifted to the right by 1
       *
       */
      
      //       std::cout << "Overwrite = false" << std::endl;
      //       std::cout << "Print out entries of helper min/max list" << std::endl;
      //
      //       for ( UInt i=0; i<stringLength; i++ ) {
      //         std::cout << "index " << i << ": " << helpStringEl[i] << std::endl;
      //       }
      //       std::cout << "#############" << std::endl;
      Double pixelToAdd = everettPixel(-helpStringEl[0], helpStringEl[0]);
      newY = pixelToAdd;
      for ( UInt i=0; i<stringLength-1; i++ ) {
				pixelToAdd = 2.0*everettPixel(helpStringEl[i],helpStringEl[i+1]);
        newY +=  pixelToAdd;
      }
      
      Double X_norm_unclipped = Xin / XSaturated_;
      LOG_DBG(scalpreisachInversion) << "Eval Preisach - Add anhystPart " << evalAnhystPart_normalized(X_norm_unclipped);
      LOG_DBG(scalpreisachInversion) << "for X/norm = " << (X_norm_unclipped);
      newY += evalAnhystPart_normalized(X_norm_unclipped);
      
      //newY += anhyst_A_*std::atan(anhyst_B_*X_norm_unclipped) + anhyst_C_*X_norm_unclipped;
    }
    LOG_TRACE(scalpreisachInversion) << "Computed new value: " << newY;
    //std::cout << "UpdateMinMaxList - Output: " << newY << std::endl;
    return newY;
  }
  
  
  Double Preisach::everettPixel(Double val1, Double val2)
  {
    
    Double X1 = std::max(val1,val2);
    Double X2 = std::min(val1,val2);
    
    UInt M = preisachWeights_.GetNumRows();
    Double delta = 2.0 / ( (Double) M );
    
    //  std::cout << "delta: " << delta << std::endl;
    //compute index for X1 (alpha)
    Integer idx1 = -1;
    Double alpha = -1.0;
    while ( alpha <= X1 ) {
      idx1++;
      alpha += delta;
    }
    
    if (alpha > 1.0) {
      alpha = 1.0;
      idx1 = M-1;
    }
    
    //compute index for X2 (beta)
    Integer idx2 = -1;
    Double  beta = -1.0;
    while ( beta <= X2 ) {
      idx2++;
      beta += delta;
    }
    beta -= delta;
    
    // X1 >= X2
    // -> mit alpha_0 = beta_0 = -1 und gleichem delta -> idx1 >= idx2
    
    //    std::cout << "idx1: " << idx1 << std::endl;
    //    std::cout << "idx2: " << idx2 << std::endl;
    
    Double area = 0.0;
    if ( idx1 >= 0 ) {
      UInt start = std::min(idx1,idx2);
      UInt stop  = std::max(idx1,idx2);
      
      // why do we iterate over the whole square and not over the triangle?
      for ( UInt i=start; i<=stop; i++ ) {
        for ( UInt j=start; j<=stop; j++ ) {
          area += preisachWeights_[i][j];
        }
      }
      
      area *= 0.5*delta*delta;
      
      //reduce the computed area
      Double diffX1 = alpha - X1;
      Double diffX2 = X2    - beta;
      Double minusArea;
      
      /*
       * The idea behind the minusArea is to handle input variation which are smaller than delta
       * Assume for example a min/max list with differences between mins and maxs < delta
       * In that case, the same preisach element would flip from +1 to -1 and the single steps would cancel
       * out in the sum. By reducing the weighting areas accordingly, we can also treat input lists of the described form
       * as now the summed up terms are weighted with different areas!
       * Commit out the section below and you will see, that it does not work anymore!
       */
      
      //check, if we are already on the diagonal!!
      if ( idx1 == idx2 ) {
        minusArea = (   diffX2 * (delta - 0.5*diffX2)
                + diffX1 * (delta - 0.5*diffX1)
                - diffX1*diffX2
                ) * preisachWeights_[idx1][idx2];
      }
      else {
        minusArea = ( (diffX1+diffX2 )*delta - diffX1*diffX2 )
                * preisachWeights_[idx1][idx2];
        
        Integer idx = idx1-1;
        while ( idx > idx2 ) {
          minusArea += diffX2 * delta * preisachWeights_[idx][idx2];
          idx--;
          //	  std::cout << "minusArea2=" << minusArea << std::endl;
        }
        minusArea += ( delta*diffX2 - 0.5*diffX2*diffX2 )
                * preisachWeights_[idx][idx2];
        
        idx = idx2 + 1;
        while ( idx < idx1 ) {
          minusArea += diffX1 * delta * preisachWeights_[idx1][idx];
          idx++;
        }
        minusArea += ( delta*diffX1 - 0.5*diffX1*diffX1 )
                * preisachWeights_[idx1][idx];
      }
      
      area -= minusArea;
    }
    
    //sgn-function
    if ( val2 < val1 ) {
      area *= -1.0;
    }
    
    return area;
  }
  
  // new name that better describes function
  // Xin gets not only normalized by dividing it with Xsaturated
  // but also gets clipped to +/-1 if Xsaturated is exceeded
  Double Preisach::normalizeAndClipInput(Double Xin)
  {
    
    Double Xout;
    
    if ( Xin > XSaturated_ ) {
      //saturation achieved!!
      Xout = 1.0;
    }
    else if ( Xin < -XSaturated_ ) {
      //saturation achieved!!
      Xout = -1.0;
    }
    else {
      //normalize input
      Xout = Xin / XSaturated_;
    }
    
    return Xout;
  }
  
}
