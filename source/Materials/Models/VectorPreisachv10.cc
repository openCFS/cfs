#include "VectorPreisachv10.hh"

#include <fstream>
#include <iostream>
#include <string>
#include <boost/algorithm/string.hpp>
#include "DataInOut/Logging/LogConfigurator.hh"
/* See VectorPreisachv10_ListApproach.hh for detailed description */

namespace CoupledField
{ 
  DECLARE_LOG(vecpreisach)
  DEFINE_LOG(vecpreisach, "vecpreisach")
  DECLARE_LOG(vecpreisachInversion)
  DEFINE_LOG(vecpreisachInversion, "vecpreisachInversion")

  /*
   * BASE CLASS FUNCTIONS
   */
  VectorPreisachv10::VectorPreisachv10(Integer numElem, Double xSat, Double ySat,
          Matrix<Double>& preisachWeight, Double rotationalResistance , UInt dim, bool isVirgin,
          bool classical, Double angularDistance, Double angResolution, Double anhystA, Double anhystB, Double anhystC, bool anhystOnly)
  : Hysteresis(numElem,xSat,ySat,anhystA,anhystB,anhystC,anhystOnly)
  {
    
    /*
     * Global quantities, i.e. the same for all FE elements of the same material
     */
    if (xSat > 0 ) {
      XSaturated_  = xSat;
    }
    else {
      XSaturated_  = 1.0;
    }
    
    YSaturated_  = ySat;
    
    isVirgin_    = isVirgin;
    
    // default: no measurements; no output; new mapping
    mappingVersion_ = 1;
    performanceMeasurement_ = 0;
    // TODO: change all occurances of debug output to logging system, then remove flag
    textOutputLevel_ = 0;
    
    numElem_ = numElem;
    dim_ = dim;
    
    preisachWeights_ = preisachWeight;
    
    tol_ = 1e-16;
    // restriction to halfspace not working well; better without
		restrictToHalfspace_ = !true;
     
    // get size of preisachWeights_
    UInt M = preisachWeights_.GetNumRows();
    UInt N = preisachWeights_.GetNumCols();
    
    if(M != N){
      EXCEPTION("Matrix preisachWeight has dim " << M << " x " << N << " and thus is not symmetric!");
    }
   
    numRows_ = M;
    
    isSymmetric_ = true;
    // check symmetry w.r.t. alpha = -beta
    for(UInt i = 0; i < numRows_; i++){
      for(UInt k = 0; k < numRows_-i; k++){
        // iterate over triangle -1 < alpha < 1; -1 < beta < -alpha
        // in indices: 0 < i < numRows; 0 < k < numRows-i
        if(abs(preisachWeights_[i][k]-preisachWeights_[numRows_-k-1][numRows_-i-1]) > tol_){
          isSymmetric_ = false;
          break;
        }
      }
    }
        
    // resolution of Preisach plane
    delta_ = 2.0/Double(numRows_);
    
    rotationalResistance_ = rotationalResistance;
    
    angularDistance_ = angularDistance;
		angResolution_ = angResolution;
		
    classical_ = classical;
    
    anhyst_A_ = anhystA;
    anhyst_B_ = anhystB;
    anhyst_C_ = anhystC;
    
    if(classical_){
      LOG_TRACE(vecpreisach) << "VectorPreisach: Using classical vector model (Sutor2012)";
    } else {
      LOG_TRACE(vecpreisach) << "VectorPreisach: Using revised vector model (Sutor2015)";
    }
    
		//std::cout << "NumElements: " << numElem_ << std::endl;
		
    /*
     * Local quantities, i.e. arrays storing different values for each FE element
     * > may only be written by compueValue_vec and only if overwrite = true
     */
    preisachSum_ = new Vector<Double>[numElem_];

    prevXVal_ = new Vector<Double>[numElem_];
    prevHVal_ = new Vector<Double>[numElem_];
    for(UInt k = 0; k < numElem_; k++){
      preisachSum_[k] = Vector<Double>(dim_);
      preisachSum_[k].Init();
      
      prevXVal_[k] = Vector<Double>(dim_);
      prevXVal_[k].Init();
      prevHVal_[k] = Vector<Double>(dim_);
      prevHVal_[k].Init();  
    }

    /*
     * Needed in context of the linesearch algorithm:
     * there we have to evaluate the hysteresis operator
     * for multiple value of eta; due to the memory property of the hysteresis operator
     * these test-steps will have a permanent effect on the output; to avoid this
     * we can use a flag which denotes if we want to overwrite or not;
     * in case of overwrite = false, we modify only a copy of the data structures and
     * store the final result in preisachSumTmp_l
     */
    preisachSumTmp_ = new Vector<Double>[numElem_];
    for(UInt k = 0; k < numElem_; k++){
      preisachSumTmp_[k] = Vector<Double>(dim_);
      preisachSumTmp_[k].Init();
    }
    
    // for inversion
    INV_maxIter_ = 35;
    INV_resTolH_ = 1e-12;
    INV_resTolB_ = 1e-12;
    INV_jacobiResolution_ = 1e-12;
    INV_useTikhonov_ = false;
    INV_alphaLSStart_ = 0.25;
    INV_angClipping_ = -1.0;
        
  }
  
  VectorPreisachv10::~VectorPreisachv10(){
    delete[] preisachSum_;
    delete[] preisachSumTmp_;
    
    delete[] prevXVal_;
    delete[] prevHVal_;
  }
//  
//  void VectorPreisachv10::ClipDirection(Vector<Double>& targetVector){
//    
//    if(INV_angClipping_ <= 0.0){
//      return; // no clipping
//    }
//    
//    /*
//     * (Former name: clipNewRotationDirection, taken from VecPreisach)
//     * New purpose March/April 2018
//     * The general clipping idea might be ok, but it does not help to clip the rotation states
//     * directly
//     *  > reason: due to the weighting via the switching state, we still get a continuous
//     *            resolution, even if the single rotations states were clipped
//     *  > new purpose: instead of clipping the single rotation states, we clip the final output
//     *            of the hyst operator!
//     *  > further improvings/changes:
//     *      1. 3d case implemented
//     *      2. angularClipping is treated as resolution in DEGREE (not rad!)
//     * 
//		 * Addon: shifted to CoefFunctionHyst; both input and output of hyst operator
//		 *				can be clipped
//     * 
//     * New function April 2017
//     * Idea: Restrict the range of possible rotation directions to fixed angular steps
//     * (i.e. x rad)
//     * Reason: Due to numerical pollutions, we might encounter rotation directions like
//     * (1.00000000,-1e-9) or (0.99999999,1e-10). These deviations from the input direction
//     * (1.0,0) seem negligble but will lead to actual problems when evaluting deltaMatrices
//     * of the form deltaP/deltaE (as deltaE normally is similarly small). Due to this, the
//     * resulting deltaMatrices may have huge permittivities/permeabilites in direction where
//     * normally no value should be. This leads to serious convergence issues.
//     *
//     */
//        
//    /*
//     * The following two clipping steps will be applied: 
//     * 1. transform to circular/spherical coordinates; restrict angles (in deg) to angularClipping
//     * 2. check for special angles (90,180) degree and restrict vector further
//     *      (as e.g. sin(pi) != 0 due to numerics)
//     */
////    LOG_TRACE(coeffcthyst) << "Clip direction of target vector to next full " << MAT_angClipping_ << " degree";
////    LOG_DBG(coeffcthyst) << "Original vector: " << targetVector.ToString();
//    if(dim_ == 2){
//      /*
//       * use polar/circular coordinates
//       * x = r cos(alpha)
//       * y = r sin(alpha)
//       */
//      Double radius, tmp, alphaDeg;
//      radius = targetVector.NormL2();
//      
//      if(targetVector[0] == 0){
//        if(targetVector[1] > 0){
//          alphaDeg = 90.0; 
//        } else {
//          alphaDeg = -90.0;
//        }
//      } else if(targetVector[0] > 0){
//        tmp = atan2(targetVector[1],targetVector[0]);
//        alphaDeg = tmp*180/M_PI;
//      } else {
//        tmp = -asin(targetVector[1]/radius);
//        alphaDeg = tmp*180/M_PI + 180.0;
//      }
////      LOG_DBG(coeffcthyst) << "Circular/Polar coordinates (r,alpha): " << radius << "," << alphaDeg;
//      
//      // apply clipping
//      tmp = alphaDeg/INV_angClipping_;
//      tmp = round(tmp);
//      alphaDeg = tmp*INV_angClipping_;
//      
////      LOG_DBG(coeffcthyst) << "Circular/Polar coordinates after clipping (r,alpha): " << radius << "," << alphaDeg;
//      
//      // now rebuild output vector with clipped coordinates
//      targetVector[0] = radius*cos(alphaDeg/180*M_PI);
//      targetVector[1] = radius*sin(alphaDeg/180*M_PI);
//      
////      LOG_DBG(coeffcthyst) << "Rebuild vector after clipping: " << targetVector.ToString();
//      
//      // finally check for special cases i.e. 90/180 deg (which can not perfectly be reproduced by computing cos()/sin()
//      if( abs(alphaDeg - 90) < INV_angClipping_/1000.0 ){
//        // alpha = 90 > positive y axis
//        targetVector[0] = 0.0;
//        targetVector[1] = radius;
//      } else if( abs(alphaDeg + 90) < INV_angClipping_/1000.0 ){
//        // alpha = -90 > negative y axis
//        targetVector[0] = 0.0;
//        targetVector[1] = -radius;
//      } else if( (abs(alphaDeg - 180) < INV_angClipping_/1000.0) || (abs(alphaDeg + 180) < INV_angClipping_/1000.0) ){
//        // alpha = +/- 180 deg > negative x axis
//        targetVector[0] = -radius;
//        targetVector[1] = 0.0;
//      }
//      
////      LOG_DBG(coeffcthyst) << "Rebuild vector after further treatment: " << targetVector.ToString();
//    } else {
//      /*
//       * use spherical coordinates
//       * x = r sin(theta) cos(phi)
//       * y = r sin(theta) sin(phi)
//       * z = r cos(theta)
//       */
//      Double radius, tmp, phiDeg, thetaDeg;
//      radius = targetVector.NormL2();
//      
//      tmp = acos(targetVector[2]/radius);
//      thetaDeg = tmp*180/M_PI;
//      tmp = atan2(targetVector[1],targetVector[0]);
//      phiDeg = tmp*180/M_PI;
//      
////      LOG_DBG(coeffcthyst) << "Spherical coordinates (r,theta,phi): " << radius << "," << thetaDeg << "," << phiDeg;
//      
//      // apply clipping
//      tmp = thetaDeg/INV_angClipping_;
//      tmp = round(tmp);
//      thetaDeg = tmp*INV_angClipping_;
//      
//      tmp = phiDeg/INV_angClipping_;
//      tmp = round(tmp);
//      phiDeg = tmp*INV_angClipping_;
//      
////      LOG_DBG(coeffcthyst) << "Spherical coordinates after clipping (r,theta,phi): " << radius << "," << thetaDeg << "," << phiDeg;
//      
//      targetVector[0] = radius*sin(thetaDeg/180*M_PI)*cos(phiDeg/180*M_PI);
//      targetVector[1] = radius*sin(thetaDeg/180*M_PI)*sin(phiDeg/180*M_PI);
//      targetVector[2] = radius*cos(thetaDeg/180*M_PI);
//      
////      LOG_DBG(coeffcthyst) << "Rebuild vector after clipping: " << targetVector.ToString();
//      
//      // finally check for special cases i.e. 90/180 deg (which can not perfectly be reproduced by computing cos()/sin()
//      if( (abs(thetaDeg - 90) < INV_angClipping_/1000.0 ) || (abs(thetaDeg + 90) < INV_angClipping_/1000.0 ) ){
//        // theta = +/- 90 > z = 0
//        targetVector[2] = 0.0;
//      } else if( (abs(thetaDeg - 180) < INV_angClipping_/1000.0) || (abs(thetaDeg + 180) < INV_angClipping_/1000.0) ){
//        // theta = +/- 180 deg > negative z-axis
//        targetVector[0] = 0.0;
//        targetVector[1] = 0.0;
//        targetVector[2] = -radius;
//      }
//      if( (abs(phiDeg - 90) < INV_angClipping_/1000.0 ) || (abs(phiDeg + 90) < INV_angClipping_/1000.0 ) ){
//        // phi = +/- 90 > x = 0
//        targetVector[0] = 0.0;
//      } else if( (abs(phiDeg - 180) < INV_angClipping_/1000.0) || (abs(phiDeg + 180) < INV_angClipping_/1000.0) ){
//        // phi = +/- 180 deg > y = 0
//        targetVector[1] = 0.0;
//      }
//
////      LOG_DBG(coeffcthyst) << "Rebuild vector after further treatment: " << targetVector.ToString();
//      
//    }
//  }
  
	Vector<Double> VectorPreisachv10::restrictToHalfspace(Vector<Double>& e_u_new){
		/*
		 * Idea: restrict e_u to halfspace y>0; inputs that point into the lower
		 *				halfspace y<0 are then represented by negative xPar values
		 *	> this might reduce the amount of rotation states in case of the field
		 *		switching direcion along an axis (e.g. from -x to x)
		 *			
		 */
		Vector<Double> e_restricted = Vector<Double>(dim_);
		e_restricted = e_u_new;
		if(e_u_new[1] < 0){
			// lower halfspace > rotate around z-axis by 180 degree > switch sign of x and y component
			// > works for 2d and 3d alike
			e_restricted[0] *= -1.0;
			e_restricted[1] *= -1.0;
		} else if(e_u_new[1] == 0){
			// here we are in the x-z plane
			// > restrict to positve x
			if(e_u_new[0] < 0){
				// negative x > rotate around y-axis by 180 degree
				e_restricted[0] *= -1.0;
				if(dim_ == 3){
					e_restricted[2] *= -1.0;
				}
			} else if(e_u_new[0] == 0){
				// x and y are negative > we are on the z-axis
				// > restrict to positve z
				if(dim_ == 3){
					if(e_restricted[2] < 0){
						e_restricted[2] *= -1.0;
					}
				}
			}
		}
		return e_restricted;
		
	}
	  
  Vector<Double> VectorPreisachv10::evaluateNewRotationDirection(Vector<Double>& e_u_new, Vector<Double>& e_u_old, Double xVal){
    /*
     * calculates the new rotation direction for overwritten rotation states according to the
     * revised vector model from 2015
     *
     * Idea: the old rotation state e_u_old is not simply overwritten, but instead it rotates in direction of e_u_new
     *       it will however, only rotate up to an angle deltaPhi which depends on the magnitude of the current input
     *
     * Further addition:
     *       an easy axis might be given via input file; in that case we allow the material to rotate easier into that diretion
     *       an harder into the perpendicular direction -> theory not yet created
     */
    if(classical_){
      return e_u_new;
    } else {
      /*
       * Evaluate the new rotation direction in the following steps:
       * 1. check for the angular distance delta_phi which shall be kept between the new input direction e_u_new and the
       *    rotation direction to be taken e_phi
       * 2. check angle alpha between e_u_new and the currently stored direction e_u_old
       * 3. if delta_phi = 0
       *      return e_u_new
       *    else if delta_phi < alpha
       *      return e_u_old
       *    else
       *      rotate e_u_old towards e_u_new but keep an angular distance of deltaphi
       *        (or rotate e_u_new by deltaphi towards e_u_old)
       */
      Double delta_phi, alpha;
      
      /*
       * 0. Check old rotation state: if it is zero (i.e. unset yet), rotate completely to new direction
       */
      if(e_u_old.NormL2() < tol_){
        return e_u_new;
      }
      /*
       * 1. calculate delta_phi based on (10) Sutor2015
       * Note: original formula uses 2*abs(xVal) due to normalization to range [-0.5,0.5] instead of [-1,1];
       * xVal not normalized yet
       *
       * angularDistance_ in deg
       */
      xVal /= XSaturated_;
      delta_phi = angularDistance_ * (1 - abs(xVal));
      //      std::cout << "delta_phi: " << delta_phi << std::endl;
      
      /*
       * 2. get angle between e_u_new and e_u_old
       * Note: e_u_new and e_u_old already have length 1
       * Note2: delta_phi is in degree, so alpha has to be in degree, too
       */
      Double tmp = e_u_old.Inner(e_u_new);
      //      std::cout << "e_u_old: " << e_u_old.ToString() << std::endl;
      //      std::cout << "e_u_new: " << e_u_new.ToString() << std::endl;
      //      std::cout << "e_u_old.Inner(e_u_new): " << tmp << std::endl;
      /*
       * due to rounding errors, the scalar product of two vectors of lenght 1
       * can be larger than 1 -> cut down to avoid NaN
       */
      if(tmp > 1.0){
        tmp = 1.0;
      } else if (tmp < -1.0){
        tmp = -1.0;
      }
      alpha = std::acos(tmp)*180/M_PI;
      
      //      std::cout << "alpha: " << alpha << std::endl;
      /*
       * NEW: 15.2.2017
       * If alpha (the angle between new and old state is close to zero,
       * we do not have to rotate at all as the rotation states are nearly the same)
       */
      // take much larger tolerance here (should we really go up to 1e-15?!)
      if(abs(alpha) < angResolution_){
        /*
         * take new state
         */
        return e_u_new;
      }
      
      /*
       * 3. calculate new rotation direction depending on delta_phi and alpha
       */
      if(delta_phi < angResolution_){
				//        std::cout << "delta_phi < angleTol > return e_u_new" << std::endl;
        /*
         * no resistance to rotation
         */
        return e_u_new;
      } else if(delta_phi >= alpha) {
				//        std::cout << "delta_phi >= angleTol > return e_u_old" << std::endl;
        /*
         * e_u_old is already closer than the resistance angle delta_phi that should remain
         */
        return e_u_old;
      } else {
        /*
         * construct new rotation vector e_phi
         * Note: e_u_new is rotated by -delta_phi towards e_u_old
         */
        Vector<Double> e_phi = Vector<Double>(dim_);
        Matrix<Double> rotMat = Matrix<Double>(dim_,dim_);
        Double c,s;
        
        if(dim_ == 2){
          /*
           * in 2d we need
           *  1. direction of rotation axis (+z or -z)
           *  2. a 2x2 rotation matrix describing a rotation of -delta_phi degree around rotation axis
           */
          
          /*
           * in 2d rotation axis is either +z or -z; find sign via cross product of e_u_old and e_u_new but take only 3rd entry
           * Note: we build the right-hand-system e_u_old, e_u_new, e_u_old x e_u_new
           */
          Double n3;
          n3 = e_u_old[0]*e_u_new[1] - e_u_old[1]*e_u_new[0];
          
          if(n3 != 0){
            n3 = n3/abs(n3);
          }
          
          //std::cout << "n3: " << n3 << std::endl;
          
          c = std::cos(-delta_phi/180*M_PI);
          s = std::sin(-delta_phi/180*M_PI);
          
          //          std::cout << "c: " << c << std::endl;
          //          std::cout << "s: " << s << std::endl;
          
          rotMat[0][0] = c;
          rotMat[0][1] = -n3*s;
          rotMat[1][0] = n3*s;
          rotMat[1][1] = c;
          
          e_phi = rotMat*e_u_new;
          
        } else {
          /*
           * in 3d we need
           *  1. rotation axis given by e_u_old x e_u_new
           *  2. a 3x3 rotation matrix describing a rotation of -delta_phi degree around rotation axis
           */
          Vector<Double> normal = Vector<Double>(dim_);
          e_u_old.CrossProduct(e_u_new,normal);
          
          if(normal.NormL2() != 0){
            normal = normal/normal.NormL2();
          }
          
          c = std::cos(-delta_phi/180*M_PI);
          s = std::sin(-delta_phi/180*M_PI);
          
          /*
           * rotation matrix for a rotation around 'normal';
           * Note: e_u_old, e_u_new and normal have to have the same origin
           * (which is fulfilled here as all are vectors starting at (0,0,0))
           */
          rotMat[0][0] = normal[1]*normal[1]*(1-c) + c;
          rotMat[0][1] = normal[1]*normal[2]*(1-c) - normal[3]*s;
          rotMat[0][2] = normal[1]*normal[3]*(1-c) + normal[2]*s;
          
          rotMat[1][0] = normal[2]*normal[1]*(1-c) + normal[3]*s;
          rotMat[1][1] = normal[2]*normal[2]*(1-c) + c;
          rotMat[1][2] = normal[2]*normal[3]*(1-c) - normal[1]*s;
          
          rotMat[2][0] = normal[3]*normal[1]*(1-c) - normal[2]*s;
          rotMat[2][1] = normal[3]*normal[2]*(1-c) + normal[1]*s;
          rotMat[2][2] = normal[3]*normal[3]*(1-c) + c;
          
          e_phi = rotMat*e_u_new;
        }
        return e_phi;
      }
    }
  }
  
  Vector<Double> VectorPreisachv10::computeAbsResidualX(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& hystVal, Matrix<Double> mu_inv){
    /*
     *    yVal = mu*xVal + hystVal(xVal)
     * 
     *  > res = xVal + hystVal(xVal)/mu - yVal/mu
     */
    Vector<Double> res = Vector<Double>(dim_);
    Vector<Double> tmp = Vector<Double>(dim_);
    tmp = hystVal;
    tmp.Add(-1.0,yVal);
    mu_inv.Mult(tmp,res);
    res.Add(1.0,xVal);
    
    //    for(UInt i = 0; i < dim_; i++){
    //      // assume that mu is diagonal and has entries != 0
    //      if(mu[i][i] != 0){
    //        res[i] += 1.0/mu[i][i]*hystVal[i];
    //        res[i] += -1.0/mu[i][i]*yVal[i];
    //      }
    //    }
    
    return res;
  }
  
  Vector<Double> VectorPreisachv10::computeResidual(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& hystVal, 
          Matrix<Double> mu, Matrix<Double> mu_inv, bool wrtX, bool relative){
    
    Vector<Double> ret = computeAbsResidualX(xVal, yVal, hystVal, mu_inv);
    if(wrtX){
      /*
       * wrtX = true
       * 
       * > resAbsX = xVal - yVal/mu + hystVal/mu
       * > resRelX = resAbsX/xVal.NormL2()
       */
      if(relative){
        Double xValNorm = xVal.NormL2();
        if(xValNorm != 0){
          ret.ScalarDiv(xValNorm);
        } else {
          WARN("relative residual wrt x requested, but xValNorm == 0; return abs value instead");
        }
      }
    } else {
      /*
       *  wrtX = false
       * 
       * > resAbsY = yVal - hystVal - mu*xVal = -mu*resAbsX
       * > resRelY = reaAbsY/yVal.NormL2()
       * 
       */
      Vector<Double> tmp = Vector<Double>(dim_);
      //ret = ret*(-mu);
      mu.Mult(ret,tmp);
      tmp.ScalarMult(-1.0);
      ret = tmp;
      //      for(UInt i = 0; i < dim_; i++){
      //        ret[i] = -1.0*ret[i]*mu[i][i];
      //      }
      
      if(relative){
        Double yValNorm = yVal.NormL2();
        if(yValNorm != 0){
          ret.ScalarDiv(yValNorm);
        } else {
          WARN("relative residual wrt y requested, but yValNorm == 0; return abs value instead");
        }
      }
    }
    
    return ret;
  }
  
  Matrix<Double> VectorPreisachv10::computeJacobianOfAbsResidualX(Vector<Double>& xVal, Vector<Double>& hystVal, 
          Matrix<Double> mu_inv, Integer operatorIdx, Double sign, UInt implementation, 
					bool overwriteMemory, bool overwriteDirection) {
    
    LOG_TRACE(vecpreisachInversion) << " --------- computeJacobianOfAbsResidualX --------- ";
    LOG_DBG(vecpreisach) << "VecPreisach::computeJacobianOfAbsResidualX";
    if(xVal.NormL2() >= XSaturated_){
      EXCEPTION("xVal.NormL2() >= XSaturated_");
    }
    
    Double deltaX = 0.0;
		Double scal = INV_jacobiResolution_;
    Double deltaXmin = scal*XSaturated_; 
   
//		Double deltaXmin = 1e-8*XSaturated_; 
//    Double scal = 1e-8;//1e-5;
//    
//		
//    bool overwriteMemory = false;
    Vector<Double> xShifted;
    Vector<Double> hystShifted;
    Vector<Double> deltaHyst;
    Matrix<Double> jac = Matrix<Double>(dim_,dim_);
    jac.Init();
    
    if(implementation == 0){
      LOG_DBG(vecpreisach) << "Use forward/backward differences for approximation of Jacobian";
      /*
       * Full Jacobian using forward/backward differences
       */
      for(UInt i = 0; i < dim_; i++){
        xShifted = xVal;
        
        //        if( xVal[i] < 0 ){
        //          deltaX = sign*std::min( -scal*XSaturated_, -deltaXmin );
        //        } else {
        //          deltaX = sign*std::max( scal*XSaturated_, deltaXmin );
        //        }
        if( xVal[i] < 0 ){
          deltaX = sign*std::min( scal*xVal[i], -deltaXmin );
        } else {
          deltaX = sign*std::max( scal*xVal[i], deltaXmin );
        }
        
        xShifted[i] += deltaX;
        
        Vector<Double> dXvec = xShifted;
        dXvec.Add(-1.0,xVal);
        LOG_TRACE(vecpreisachInversion) << " xVal " << xVal.ToString();
        LOG_TRACE(vecpreisachInversion) << " xShifted " << xShifted.ToString(); 
        LOG_TRACE(vecpreisachInversion) << " dXvec " << dXvec.ToString(); 
        hystShifted = computeValue_vec(xShifted, operatorIdx, overwriteMemory, overwriteDirection);
        /*
         * Compute Jacobian for residual wrt x
         * 
         * jac_ji = + delta_ji + mu_inv[j][:]*dhystVal/dxVal_i
         */   
        deltaHyst = hystShifted;
        deltaHyst.Add(-1.0,hystVal);
        
        LOG_TRACE(vecpreisachInversion) << " hystVal " << hystVal.ToString();
        LOG_TRACE(vecpreisachInversion) << " hystShifted " << hystShifted.ToString();
        LOG_TRACE(vecpreisachInversion) << " deltaHyst " << deltaHyst.ToString();
//        Vector<Double> hystValRecomputed = computeValue_vec(xVal, operatorIdx, overwriteMemory, overwriteDirection);
//        LOG_TRACE(vecpreisachInversion) << " hystValRecomputed " << hystValRecomputed.ToString();
        
//        Vector<Double> xOppShift = xShifted;
//        xOppShift[i] -= 2*deltaX;
//        dXvec = xOppShift;
//        dXvec.Add(-1.0,xVal);
//        LOG_TRACE(vecpreisachInversion) << " xOppShift " << xOppShift.ToString(); 
//        LOG_TRACE(vecpreisachInversion) << " dXvec " << dXvec.ToString(); 
//        Vector<Double> hystOppShifted = computeValue_vec(xOppShift, operatorIdx, overwriteMemory, overwriteDirection);
//        LOG_TRACE(vecpreisachInversion) << " hystOppShifted " << hystOppShifted.ToString();
        
        Vector<Double> curCol = Vector<Double>(dim_);
        mu_inv.Mult(deltaHyst,curCol);
        curCol.ScalarDiv(deltaX);
        
        jac[i][i] = 1.0;
        for(UInt j = 0; j < dim_; j++){
          jac[j][i] += curCol[j];
          
          //          //jac[j][i] += (hystShifted[j]-hyst[j])/deltaX/mu; 
          //          // assume matrix mu to be diagonal
          //          // > here we had an error by taking mu[i][i] instead of mu[j][j]
          //          if(mu[j][j]!=0){
          //            jac[j][i] += (hystShifted[j]-hystVal[j])/(xShifted[i]-xShifted_opp[i])/mu[j][j]; 
          //          }
        }
      }
      
      LOG_TRACE(vecpreisachInversion) << " Jacobian " << jac.ToString();
      Matrix<Double> jacT = Matrix<Double>(dim_,dim_);
      jac.Transpose(jacT);
      Matrix<Double> jacTjac = Matrix<Double>(dim_,dim_);
      jacT.Mult(jac,jacTjac);
      LOG_TRACE(vecpreisachInversion) << " jacTjac " << jacTjac.ToString();
      
      Double det;
      jacTjac.Determinant(det);
      LOG_TRACE(vecpreisachInversion) << " det(jacTjac) " << det;
//      if(abs(det) < 1e-14){
//        // repeat with debug output
//        for(UInt i = 0; i < dim_; i++){
//          xShifted = xVal;
//          
//          //        if( xVal[i] < 0 ){
//          //          deltaX = sign*std::min( -scal*XSaturated_, -deltaXmin );
//          //        } else {
//          //          deltaX = sign*std::max( scal*XSaturated_, deltaXmin );
//          //        }
//          if( xVal[i] < 0 ){
//            deltaX = sign*std::min( scal*xVal[i], -deltaXmin );
//          } else {
//            deltaX = sign*std::max( scal*xVal[i], deltaXmin );
//          }
//          
//          xShifted[i] += deltaX;
//          
//          Vector<Double> dXvec = xShifted;
//          dXvec.Add(-1.0,xVal);
//          LOG_TRACE(vecpreisachInversion) << " xVal " << xVal.ToString();
//          LOG_TRACE(vecpreisachInversion) << " xShifted " << xShifted.ToString(); 
//          LOG_TRACE(vecpreisachInversion) << " dXvec " << dXvec.ToString(); 
//          hystShifted = computeValue_vec(xShifted, operatorIdx, overwriteMemory, overwriteDirection,true);
//          /*
//           * Compute Jacobian for residual wrt x
//           * 
//           * jac_ji = + delta_ji + mu_inv[j][:]*dhystVal/dxVal_i
//           */   
//          deltaHyst = hystShifted;
//          deltaHyst.Add(-1.0,hystVal);
//          
//          LOG_TRACE(vecpreisachInversion) << " hystVal " << hystVal.ToString();
//          LOG_TRACE(vecpreisachInversion) << " hystShifted " << hystShifted.ToString();
//          LOG_TRACE(vecpreisachInversion) << " deltaHyst " << deltaHyst.ToString();
//          Vector<Double> hystValRecomputed = computeValue_vec(xVal, operatorIdx, overwriteMemory, overwriteDirection,true);
//          LOG_TRACE(vecpreisachInversion) << " hystValRecomputed " << hystValRecomputed.ToString();
//          
//          Vector<Double> xOppShift = xShifted;
//          xOppShift[i] -= 2*deltaX;
//          dXvec = xOppShift;
//          dXvec.Add(-1.0,xVal);
//          LOG_TRACE(vecpreisachInversion) << " xOppShift " << xOppShift.ToString(); 
//          LOG_TRACE(vecpreisachInversion) << " dXvec " << dXvec.ToString(); 
//          Vector<Double> hystOppShifted = computeValue_vec(xOppShift, operatorIdx, overwriteMemory, overwriteDirection,true);
//          LOG_TRACE(vecpreisachInversion) << " hystOppShifted " << hystOppShifted.ToString();
//          
//          Vector<Double> curCol = Vector<Double>(dim_);
//          mu_inv.Mult(deltaHyst,curCol);
//          curCol.ScalarDiv(deltaX);
//          
//          jac[i][i] = 1.0;
//          for(UInt j = 0; j < dim_; j++){
//            jac[j][i] += curCol[j];
//            
//            //          //jac[j][i] += (hystShifted[j]-hyst[j])/deltaX/mu; 
//            //          // assume matrix mu to be diagonal
//            //          // > here we had an error by taking mu[i][i] instead of mu[j][j]
//            //          if(mu[j][j]!=0){
//            //            jac[j][i] += (hystShifted[j]-hystVal[j])/(xShifted[i]-xShifted_opp[i])/mu[j][j]; 
//            //          }
//          }
//        }
//        
//      }
      
      
    } else if(implementation == 1){
      LOG_DBG(vecpreisach) << "Use central differences for approximation of Jacobian";
      /*
       * Full Jacobian using central differences
       */
      Vector<Double> xShifted_opp;
      Vector<Double> hystShifted_opp;
      
      for(UInt i = 0; i < dim_; i++){
        xShifted = xVal;
        
        //        if( xVal[i] < 0 ){
        //          deltaX = sign*std::min( -scal*XSaturated_, -deltaXmin );
        //        } else {
        //          deltaX = sign*std::max( scal*XSaturated_, deltaXmin );
        //        }
        if( xVal[i] < 0 ){
          deltaX = sign*std::min( scal*xVal[i], -deltaXmin );
        } else {
          deltaX = sign*std::max( scal*xVal[i], deltaXmin );
        }
        
        xShifted[i] += deltaX;
        hystShifted = computeValue_vec(xShifted, operatorIdx, overwriteMemory, overwriteDirection); 
        xShifted_opp = xVal;
        xShifted_opp[i] -= deltaX;
        hystShifted_opp = computeValue_vec(xShifted_opp, operatorIdx, overwriteMemory, overwriteDirection);
        
        /*
         * Compute Jacobian for residual wrt x
         * 
         * jac_ji = + delta_ji + mu_inv[j][:]*dhystVal/dxVal_i
         */   
        deltaHyst = hystShifted;
        deltaHyst.Add(-1.0,hystShifted_opp);
        
        Vector<Double> curCol = Vector<Double>(dim_);
        mu_inv.Mult(deltaHyst,curCol);
        curCol.ScalarDiv(2*deltaX);
        
        jac[i][i] = 1.0;
        for(UInt j = 0; j < dim_; j++){
          jac[j][i] += curCol[j];
          
          //          //jac[j][i] += (hystShifted[j]-hyst[j])/deltaX/mu; 
          //          // assume matrix mu to be diagonal
          //          // > here we had an error by taking mu[i][i] instead of mu[j][j]
          //          if(mu[j][j]!=0){
          //            jac[j][i] += (hystShifted[j]-hystVal[j])/(xShifted[i]-xShifted_opp[i])/mu[j][j]; 
          //          }
        }
      }
    } 
    
    LOG_DBG(vecpreisach) << "Retrieved Jacobian: " << jac.ToString();
    return jac;
  }
  
  Matrix<Double> VectorPreisachv10::computeJacobian(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& hyst, Vector<Double>& resX,
          Matrix<Double> mu, Matrix<Double> mu_inv, Integer operatorIdx, Double sign, bool wrtX, bool relative, 
					UInt implementation, bool overwriteMemory, bool overwriteDirection){

    Matrix<Double> jac = computeJacobianOfAbsResidualX(xVal, hyst, mu_inv, 
												operatorIdx, sign, implementation, overwriteMemory, overwriteDirection);
    
    if(wrtX){
      /*
       * wrtX = true
       * 
       * jacAbsX_ij = d resAbsX_i / d x_j
       *            = d ( x_i - y_i/mu + hyst_i/mu ) / d x_j
       *            ~ delta_ij + 1/mu (hyst(x+Dx_j*e_j) - hyst(x))_i / Dx_j 
       *            > computeJacobianOfAbsResidualX(xVal, hyst, mu, idElem, sign);
       * 
       * jacRelX_ij = d resRelX_i / d x_j
       *            
       *            = d (reaAbsX_i / xVal.NormL2()) / d x_j
       * 
       *            = ( xVal.NormL2()* d resAbsX_i / d x_j - resAbsX_i* d xVal.NormL2() / d x_j ) * (1/xVal.NormL2()^2)
       * 
       *            = jacAbsX_ij / xVal.NormL2() - resAbsX_i * xVal_j / (xVal.NormL2()^3) 
       * 
       *            = jacAbsX_ij / xVal.NormL2() - resRelX_i * xVal_j / (xVal.NormL2()^2) 
       * 
       */
      if(relative){
        Double xValNorm = xVal.NormL2();
        if(xValNorm != 0){
          jac = jac*(1.0/xValNorm);
          
          Double xValNormQuad = xValNorm*xValNorm;
          for(UInt i = 0; i < dim_; i++){
            for(UInt j = 0; j < dim_; j++){
              // resX has to be relative wrt x, too!
              jac[i][j] -= (resX[i]*xVal[j])/xValNormQuad;
            }
          }
        } else {
          WARN("Jacobian of relative residual wrt x requested, but xValNorm == 0; return Jacobian of abs value instead");
        }
      }
    } else {
      /*
       * wrtX = false
       * 
       * jacAbsY_ij = d resAbsY_i / d x_j
       *            = d (-mu*resAbsX_i) / d x_j
       *            = -mu * jacAbsX_ij
       * 
       * jacRelY_ij = d resRelY_i / d x_j
       *            = 1/yVal.normL2() * jacAbsY_ij
       *            
       */
      //jac = jac*(-mu);
      Matrix<Double> tmp = Matrix<Double>(dim_,dim_);
      tmp = jac;
      mu.Mult(tmp,jac);
      jac = jac*(-1.0);
      
      if(relative){
        Double yValNorm = yVal.NormL2();
        if(yValNorm != 0){
          jac = jac*(1.0/yValNorm);
        } else {
          WARN("relative residual wrt y requested, but yValNorm == 0; return abs value instead");
        }
      }
    }
    
    return jac;
  }
  
  bool VectorPreisachv10::checkConvergence(Vector<Double>& res, Matrix<Double>& jacT, Double& errorNorm, Double tol){
    // According to Dahmen&Reusken - Numerik partieller DFG
    // the residual is a non-sufficient condition
    // instead we have to check the norm of jacT*res
    Vector<Double> errorVec = Vector<Double>(dim_);
    jacT.Mult(res,errorVec); 
    errorNorm = errorVec.NormL2();
    
    //std::cout << "CheckConvergence" << std::endl;
    //std::cout << "resIn: " << res.ToString() << std::endl;
    //std::cout << "jacT: " << jacT.ToString() << std::endl;
    //std::cout << "ErrorVec: " << errorVec.ToString() << std::endl;
    
    if(errorNorm <= tol){
      return true;
    } else {
      return false;
    }
    
  }
  	
	Double VectorPreisachv10::computeRho(Vector<Double>& xNew, Vector<Double>& xUpdate, 
				Vector<Double>& res, Vector<Double>& resShifted, Matrix<Double>& jac){
		
    std::stringstream computeRhoTrace;
		// According to Dahmen & Reusken, we have to check
		// our increment by computing
		// 
		// rho = (||F(x)||^2 - ||F(x + increment)||^2) /
		//       (||F(x||^2 - ||F(x) + F'(x)*increment||^2)
		//
		// F = residual
		// F' = Jacobian
		Double resNorm = res.NormL2();
		Double resShiftedNorm = resShifted.NormL2();
		Vector<Double> tmp = Vector<Double>(dim_);
		jac.Mult(xUpdate,tmp);
		
		tmp = tmp + res;
		Double tmpNorm = tmp.NormL2();
		Double a = resNorm*resNorm - resShiftedNorm*resShiftedNorm;
		Double b = resNorm*resNorm - tmpNorm*tmpNorm;
		
    computeRhoTrace << "Nominator: ||F(x)||^2 - ||F(x + increment)||^2 = " << a << std::endl;
    computeRhoTrace << "Denominator: ||F(x||^2 - ||F(x) + F'(x)*increment||^2 = " << b << std::endl;
    
		/*
		 * Note: during various tests, rho reached values between
		 *       -200 to 2
		 * > if stepping is too small: rho << 0
		 * > if stepping is too large: rho ~ 1-2
		 * > consider this for creating an appropriate adation of alpha
     * 
     * > EDIT: according to Dahmen Reusken, rho < 0 is never acceptable
     *          and b should always be >= 0
		 */
		Double rho;
		if(b != 0){
			rho = a/b;
		} else {
      LOG_DBG(vecpreisach) << "Compute rho: b = 0";
			rho = a;
      if(a == 0){
        LOG_DBG(vecpreisach) << "Compute rho: a = 0, too!";
        // rho is actually undefined here (as increment is basically 0)
        // > alpha seems to be way too large (otherwise increment would not be 0)
        // > return rho >= 1 in order to trigger a much smaller alpha
        // > rho >= 10 > reset alpha
        rho = 10.0;
      }
		}
    bool output = false;
//    /*
//     * NEW: check if increment does lead to a usable change in residual at all
//     *      (a/b might be ok, but the actual change in residual (b) might be so small
//     *      that we do not step forward)
//     * 
//     */
//    if(b < 1e-25){
//      b = 0;
//    }
//    
//    if((b/(resNorm*resNorm) < 1e-12)&&(b != 0)){
//      computeRhoTrace << "b/(resNorm*resNorm) = " << b  <<"/"<< (resNorm*resNorm) << " = " << b/(resNorm*resNorm) << " < 1e-12" << std::endl; 
//      computeRhoTrace << " > residual did not really decrease > reset alpha to 1.0" << std::endl; 
//      rho = 10.0;
//      output = true;
//    }
//    
//    
    
    if(b < 0){
      // due to the model itsel, b should always be > 0 (which is the case except when it is 0 an numerical noise comes into play (< 1e-34)
      // however, during linesearch we force the solution to be below saturation; this forcing of the update might actually lead to b < 0
      computeRhoTrace << "b < 0 > SHOULD not happen" << std::endl; 
      output = !true;
    }
    if( (b == 0)&&(a != 0) ){
      computeRhoTrace << "b == 0 > a is taken directly, correct?" << std::endl; 
      output = true;
    }
    if(output == true){
      LOG_TRACE(vecpreisachInversion) << computeRhoTrace.str();
    }
  
		LOG_DBG(vecpreisach) << "Computed rho: " << rho;
		return rho;
	}
		
  Integer VectorPreisachv10::checkIncrementOLD(Vector<Double>& xNew, Vector<Double>& xUpdate, 
		Vector<Double>& res, Vector<Double>& resShifted, Matrix<Double>& jac, Double& alpha){
    
		/*
		 * NEW: change from bool to integer to allow to distinguish more cases
		 * success < 0 --> increment not accepted; alpha too small
		 * success = 0 --> alpha fits well; accept increment
		 * success > 0 --> alpha too large but increment acceptable
		 */
    LOG_DBG(vecpreisach) << "Check increment";
    Integer success = 0;
    if(xNew.NormL2() >= XSaturated_){
      // input yVal is not in saturation (otherwise we would not be here but
      // use the simple case)
      // > new value cannot be in saturation either
      // > current update will definitely not match
      // > problem here: if xNew is in saturation, further comvergence using
      //    the jacobian will not work as the slope will be such that we 
      //    cannot leave this (wrongly obtained) saturation
      // > discard update and try with different alpha again
      LOG_DBG(vecpreisach) << "Intermediate solution (" << xNew.ToString() << ") has gone into saturation > discard!";
      alpha = alpha*2.0;
      success = -1;
    } else {
      //std::cout << "Safe case: " << std::endl;
      //std::cout << "Check increment according to residual" << std::endl;
      Double rho = computeRho(xNew, xUpdate, res, resShifted, jac);
      
      //std::cout << "rho: " << rho << std::endl;
      //      std::cout << "current alpha: " << alpha << std::endl;
      // test different stepping; the furhter away we are, the stronger
      // we will adapt stepping alpha for next time
      UInt steppingVersion = 6;
      Double factor = 1;
      Double beta0;
		  Double beta1;
			
			if(steppingVersion == 0){
				// 0 < beta0 < beta1 < 1
				beta0 = 0.25;
				beta1 = 0.75;
				// WORKS best so far
        LOG_DBG(vecpreisach) << "Adapt alpha using discrete stepping";
        /*
         * Discrete stepping; works quite well so far
         */
        // increment not accepted; increase alpha in linesearch
        if(rho < -100*beta0){
          factor = 256.0;
        } else if(rho < -10*beta0){
          factor = 128.0;
        } else if(rho < -beta0){
          factor = 64.0;
        } else if(rho < 0){
          factor = 32.0;
        } else if(rho < beta0/8.0){
          factor = 16.0;
        } else if(rho < beta0/4.0){
          factor = 8.0;
        } else if(rho < beta0/2.0){
          factor = 4.0;
        } else if(rho < beta0/1.5){
          factor = 2.0;
        } else if(rho < beta0){
          factor = 1.5; 
        } 
        
        else if(rho > 1.8*beta1){
          factor = 1.0/32.0;
        } else if(rho > 1.6*beta1){
          factor = 1.0/16.0;
        } else if(rho > 1.4*beta1){
          factor = 1.0/8.0;
        } else if(rho > 1.2*beta1){
          factor = 1.0/4.0;
        } else if(rho > 1.1*beta1){
          factor = 1.0/2.0;
        } else if(rho > beta1){
          factor = 1.0/1.5;
        } 
      } else if(steppingVersion == 3){
        LOG_DBG(vecpreisach) << "Adapt alpha using continous stepping v2";
        /*
         * continuous stepping, longer range > works best so far
         */
        beta0 = 0.15;
        beta1 = 0.65;
        if(rho < -100000*beta0){
          factor = 2048;
        } else if(rho < -0.01*beta0){
          factor = pow(2,log10(-rho/beta0)+6);
        } else if(rho < 0.0001*beta0){
          factor = 16;
        } else if(rho <= beta0){
          factor = pow(2,-log10(rho/beta0));
        }
        
        else if(rho >= 3.0*beta1){
          factor = 1.0/2048;
        } else if(rho >= beta1){
          factor = pow(2,-(rho-beta1)/0.2);
        }
      } else if(steppingVersion == 5){
        // good
        beta0 = 0.15;
        beta1 = 0.65;

        if(rho < -100000*beta0){
          factor = pow(4,log10(10*(100000*beta0+beta0)/beta0 + 1));
        } else if(rho <= beta0){
          factor = pow(4,log10(10*(-rho+beta0)/beta0 + 1));
        } 
        
        else if(rho >= 3.0*beta1){
          factor = 1.0/2048;
        } else if(rho >= beta1){
          factor = pow(2,-(rho-beta1)/0.1);
        }
  
      } else if(steppingVersion == 6){
        // good, too
				beta0 = 0.15; // 0.25
				beta1 = 0.75;

        if(rho < -100000*beta0){
          factor = pow(4,log10(10*(100000*beta0+beta0)/beta0 + 1));
        } else if(rho <= beta0){
          factor = pow(4,log10(10*(-rho+beta0)/beta0 + 1));
        } 
        
        else if(rho >= 3.0*beta1){
          factor = pow(2,-(3.0*beta1-beta1)/0.2);
        } else if(rho >= beta1){
          factor = pow(2,-(rho-beta1)/0.2);
        }
  
      }
      
      if(rho >= 10){
        // reset alpha to 1
        alpha = 1.0;
        LOG_TRACE(vecpreisach) << "Reset alpha to 1.0 " << factor;
      } else {
        alpha = alpha*factor;
        LOG_TRACE(vecpreisach) << "factor " << factor;
      }
 
      if( (rho >= beta0)&&(rho <= beta1) ){
        // accept increment > according to literature we should accept all
        // increment that lead to rho > beta0
        success = 0;
      } else if (rho < beta0 ){
        // discard increment, alpha too small
        success = -1;
      } else {
				// accept increment, alpha too large
				success = 1;
			}
    }
   
    return success;
  }
  
  Integer VectorPreisachv10::checkIncrement(Vector<Double>& xNew, Vector<Double>& xUpdate, 
		Vector<Double>& res, Vector<Double>& resShifted, Matrix<Double>& jac, Double& alpha){
    /*
     * Inspired by "Levenberg-Marquart iteative regularization for th epulse-type impact-force reconstruction"
     * 
     * Trust region method with four regions for rho
     * 
     * rho = (||F(x)||^2 - ||F(x + increment)||^2) /
     *       (||F(x||^2 - ||F(x) + F'(x)*increment||^2)
     * 
     * rho < beta0 : discard update; increase alpha
     * beta0 <= rho < betaLow : accept update; increase alpha
     * betaLow <= rho < betaUp : accept update; keep alpha
     * betaUp <= rho : accept update; decrease alpha
     * 
     * Return value:
     * -1 : dicard update
     *  0 : keep update
     *  1 : keep update
     */
    Double beta0, betaLow, betaUp;
    beta0 = 0.1;
    betaLow = 0.25;
    betaUp = 0.75;
    Double factorUp = 2;
    Double factorDown = 0.5;
    
    LOG_DBG(vecpreisach) << "Check increment";
    Integer success = 0;
    if(xNew.NormL2() >= XSaturated_){
      // input yVal is not in saturation (otherwise we would not be here but
      // use the simple case)
      // > new value cannot be in saturation either
      // > current update will definitely not match
      // > problem here: if xNew is in saturation, further comvergence using
      //    the jacobian will not work as the slope will be such that we 
      //    cannot leave this (wrongly obtained) saturation
      // > discard update and try with different alpha again
      LOG_DBG(vecpreisach) << "Intermediate solution (" << xNew.ToString() << ") has gone into saturation > discard!";
      alpha = alpha*factorUp;
      success = -1;
    } else {
      //std::cout << "Safe case: " << std::endl;
      //std::cout << "Check increment according to residual" << std::endl;
      Double rho = computeRho(xNew, xUpdate, res, resShifted, jac);
      
      if(rho < beta0){
        success = -1;
        alpha = alpha*factorUp;
      } else if(rho < betaLow){
        success = 0;
        alpha = alpha*factorUp;
      } else if(rho < betaUp){
        success = 0;
      } else {
        success = 1;
        alpha = alpha*factorDown;
      }
    }
   
    return success;
  }
  
  
  
  
  bool VectorPreisachv10::performLinesearch(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& res, 
		Vector<Double>& xUpdate, Matrix<Double>& jac, Matrix<Double>& jacT, Matrix<Double> mu, Matrix<Double> mu_inv, 
		Integer operatorIdx, bool overwriteMemory, bool overwriteDirection,
		Double& alpha, Double alphaMin, Double alphaMax, bool wrtX, bool relative, UInt& numberOfIterations,Vector<Double>& xStart){
    
    LOG_TRACE(vecpreisachInversion) << " --------- START LINESEARCH --------- ";
    if(xVal.NormL2() >= XSaturated_){
      EXCEPTION("xInput to Linesearch already above saturation > must not be the case!");
    }
		LOG_DBG(vecpreisach) << "Old alpha " << alpha;
		
		UInt maxIter = 25;
    UInt itCnt = 0;
    
    Matrix<Double> matToInvert = Matrix<Double>(dim_,dim_);
    Matrix<Double> matInverted = Matrix<Double>(dim_,dim_);
    Matrix<Double> jacTjac = Matrix<Double>(dim_,dim_);
    jacT.Mult(jac,jacTjac);
    Vector<Double> jacTres_neg = Vector<Double>(dim_);
    Vector<Double> resNew = Vector<Double>(dim_);
    Vector<Double> xNew = Vector<Double>(dim_);
    Vector<Double> hystNew = Vector<Double>(dim_);
    
    jacT.Mult(res,jacTres_neg);
    jacTres_neg = jacTres_neg*(-1.0);
    
//    Double alphaMax = 256.0; //1e0;//1e1;//e10;
//    // too small not working, too large not working either e-18,e-12,e-4 > no; e-14,e-15 ok
//    Double alphaMin = 1.0/1024; //1e-10;//15; 
    
		/*
		 * new: success is now an integer
		 * <0 no success
		 * >0 success but alpha was too large
		 * =0 success alpha good
		 */
    Integer success = -1;
    bool discard = false;
    
    while(true){
      itCnt++;
			// for statistics
			numberOfIterations=itCnt;

      Double detMatToInvert;
      UInt cnt = 0;
      while(true){
        matToInvert = jacTjac;
        LOG_TRACE(vecpreisachInversion) << " jacTjac =  " << jacTjac.ToString();
        jacTjac.Determinant(detMatToInvert);
        LOG_TRACE(vecpreisachInversion) << " det(jacTjac) =  " << detMatToInvert;
        
        for(UInt i = 0; i < dim_; i++){
          matToInvert[i][i] += alpha*alpha;
        }
        matToInvert.Determinant(detMatToInvert);
        if(detMatToInvert != 0){
          break;
        } else {
          alpha = alpha*2;
        }
        cnt++;
        if(cnt > 200){
          EXCEPTION("Cannot get matrix for inversion that is invertible!");
        }
      }
      
      LOG_TRACE(vecpreisachInversion) << " alpha =  " << alpha;
      LOG_TRACE(vecpreisachInversion) << " matToInvert =  " << matToInvert.ToString();
      matToInvert.Determinant(detMatToInvert);
      LOG_TRACE(vecpreisachInversion) << " det(matToInvert) =  " << detMatToInvert;
      matToInvert.Invert(matInverted);
      LOG_TRACE(vecpreisachInversion) << " matInverted =  " << matInverted.ToString();
      
      if(INV_useTikhonov_){
        Vector<Double> tikhonovReg = Vector<Double>(dim_);
        tikhonovReg.Init();
        tikhonovReg.Add(alpha*alpha,xStart,-alpha*alpha,xVal);
        tikhonovReg.Add(1.0,jacTres_neg);
        matInverted.Mult(tikhonovReg,xUpdate);
      } else {
        matInverted.Mult(jacTres_neg,xUpdate);
      }
      assert(!xUpdate.ContainsNaN() && !xUpdate.ContainsInf());
      xNew.Init();
      
      if(xUpdate.NormL2() >= XSaturated_){
        // too large update
        // > set xNew back to saturation
        // use only 1/2 of it
        // > else oscilation between two near saturation states possible
        LOG_DBG(vecpreisachInversion) << "xUpdate > Saturation > scale down";
        LOG_DBG(vecpreisachInversion) <<  "xUpdate = " << xUpdate.ToString();
        xUpdate.ScalarDiv(2.0);
        LOG_DBG(vecpreisachInversion) <<  "xUpdate (after scaling down) = " << xUpdate.ToString();
        alpha = alpha*16;
      } 
      xNew.Add(1.0,xVal,1.0,xUpdate);
      
      if(xNew.NormL2() >= XSaturated_){
        // too large update
        // > set xNew back to saturation
        //LOG_TRACE(vecpreisachInversion) << "Reset xUpdate as x above saturation";
        xNew.ScalarMult(0.9*XSaturated_/xNew.NormL2());
        xUpdate = xNew;
        xUpdate -= xVal;
        alpha = alpha*2.0;
      }
  
      hystNew = computeValue_vec(xNew, operatorIdx, overwriteMemory, overwriteDirection);
      resNew = computeResidual(xNew,yVal,hystNew,mu,mu_inv,wrtX,relative);
        
      //std::cout << "Current iteration: " << itCnt << std::endl;
      //std::cout << "Alpha pre: " << alpha << std::endl;
      success = checkIncrement(xNew, xUpdate, res, resNew, jac, alpha);
      //std::cout << "Alpha post: " << alpha << std::endl;
      LOG_DBG(vecpreisach) << "New alpha " << alpha;
      LOG_DBG(vecpreisach) << "Success? " << success;
      if(alpha > alphaMax){
        LOG_DBG(vecpreisach) << "Maximal alpha reached > stop";
        // maximal alpha used; stop here (regardless of success)
        alpha = alphaMax;
        break;
      }
      if(alpha < alphaMin){
        LOG_DBG(vecpreisach) << "Minimal alpha reached > stop";
        // minimal alpha used; stop here
        alpha = alphaMin;
        break;
      }
	
      if(success >= 0){
        LOG_DBG(vecpreisach) << "Linesearch was successful after " << itCnt << " iterations";
        break;
      } else {
        if(itCnt >= maxIter){
          LOG_DBG(vecpreisach) << "Linesearch was not successful; Discard update.";
          discard = true;
          break;
        }
      }
    }

    return discard;
  }
  
	Vector<Double> VectorPreisachv10::computeInput_vec_withPrevStates(Vector<Double> yVal, Vector<Double> prevYval,
      Vector<Double> prevXval, Vector<Double> prevHystval, Integer operatorIndex, 
      Matrix<Double> mu, bool overwriteDirection){
		
		UInt totalNumberOfLMIterations=0;
		UInt totalNumberOfLinesearchIterations=0;
		UInt maximalNumberOfLinesearchIterations=0;
		UInt successCode=0;
    Double minAlpha,maxAlpha,avgAlpha;
		
		return computeInput_vec_withStatistics(yVal, prevYval, prevXval, prevHystval, 
			operatorIndex, mu, overwriteDirection, 
      totalNumberOfLMIterations, totalNumberOfLinesearchIterations, 
      maximalNumberOfLinesearchIterations, successCode, minAlpha,maxAlpha,avgAlpha);
	}
    
  Vector<Double> VectorPreisachv10::computeInput_vec_withStatistics(Vector<Double> yVal, Vector<Double> prevYval,
          Vector<Double> prevXval, Vector<Double> prevHystval, Integer operatorIndex, 
          Matrix<Double> mu, bool overwriteDirection, 
          UInt& totalNumberOfLMIterations, UInt& totalNumberOfLinesearchIterations, 
          UInt& maximalNumberOfLinesearchIterations, UInt& successCode, Double& minAlpha, Double& maxAlpha, Double& avgAlpha){
    
    assert(!yVal.ContainsNaN() && !yVal.ContainsInf());
		/*
     * IMPORTANT: do not pass yVal, prevYval, prevXval nor prevHystval as reference
     *						> we do not want to overwrite these values by accident!
     */
    LOG_TRACE(vecpreisachInversion) << " --------- START IVERSION --------- ";
    bool debug = true;
    std::stringstream traceMsg;
    if(debug){
      traceMsg << "VecPreisach::computeInput_vec for operator index: " << operatorIndex << std::endl;
      traceMsg << "Input value yVal: " << yVal.ToString() << std::endl;
      traceMsg << "Previous input value yVal_old: " << prevYval.ToString() << std::endl;
      traceMsg << "Previous hyst value hVal_old: " << prevHystval.ToString() << std::endl;
      traceMsg << "Previous retrieved input value xVal_old: " << prevXval.ToString() << std::endl;
      traceMsg << "Previous solution: Percentage of saturation (|xVal|/XSaturated): " << prevXval.NormL2()/XSaturated_ << std::endl;
      traceMsg << "Material Tensor: " << mu.ToString() << std::endl;
      traceMsg << "ySat: " << YSaturated_ << std::endl;
      traceMsg << "xSat: " << XSaturated_ << std::endl;
      traceMsg << "INV_maxIter_: " << INV_maxIter_ << std::endl;
      traceMsg << "INV_resTolH_: " << INV_resTolH_ << std::endl;
      traceMsg << "INV_resTolB_: " << INV_resTolB_ << std::endl;
      traceMsg << "INV_jacobiResolution_: " << INV_jacobiResolution_ << std::endl;
      traceMsg << "INV_useTikhonov_: " << INV_useTikhonov_ << std::endl;
      traceMsg << "INV_alphaLSStart_: " << INV_alphaLSStart_ << std::endl;     
      traceMsg << "INV_alphaLSMin_: " << INV_alphaLSMin_ << std::endl;  
      traceMsg << "INV_alphaLSMax_: " << INV_alphaLSMax_ << std::endl;  
    }
    
    // for statistics
    successCode = 0; 
    // 0: no success
    // 1: reused value
    // 2: simple approach
    // 3: remanence
    // 4: passed error tol
    // 5: passed res tol x
    // 6: passed res tol y
    // 9: anhyst part only
    
    totalNumberOfLMIterations = 0;
		totalNumberOfLinesearchIterations = 0;
		maximalNumberOfLinesearchIterations = 0;
    
    minAlpha = 1e16;
    maxAlpha = -1e16;
    avgAlpha = 0;
    
		// during inversion, never overwrite hyst-memory
		bool overwriteMemory = false;
    
    // return value
    Vector<Double> xVal = Vector<Double>(dim_);
    
    /*
     * Check if yVal is significantly different from previous value
     */
    Vector<Double> diff;
		diff = yVal;
    diff -= prevYval;
    
    if(diff.NormL2() < INV_resTolB_){
      traceMsg << "--A-- Inversion: Reuse old value" << std::endl;
      xVal = prevXval;
      LOG_DBG(vecpreisachInversion) << "Reused value xVal: " << xVal.ToString();
      successCode = 1;
      return xVal;
    }
		
    /*
     * Invert eps/mu for later usage
     */
    Matrix<Double> mu_inv = Matrix<Double>(dim_,dim_);
    Vector<Double> xTMP = Vector<Double>(dim_);
    mu.Invert(mu_inv);
    Double yNorm = yVal.NormL2();
    
    Vector<Double> yDir = Vector<Double>(dim_);
    yDir.Init();
    if(yNorm > 0){
      yDir = yVal;
      yDir.ScalarDiv(yNorm);
    } else {
      yDir[0] = 1.0; // default > take x-direction
      // only required to compute eps_mu
    }
    // get mu in current direction
    Vector<Double> tmp = Vector<Double>(dim_);
    Double eps_mu = 0.0;
    mu.Mult(yDir,tmp);
    yDir.Inner(tmp,eps_mu);

    if(anhystOnly_ == true){
      traceMsg << "--S-- Inversion: Special case, only anhysteretic part > solve by bisection" << std::endl;
      successCode = 9;
      xVal.Init();
      if(yNorm == 0){
        // anhyst part has no remanence
        return xVal;
      }
      // 
      Double tol = 1e-12;
      Double Xup, Xdown, Poffset, xScal;
      Xup = yNorm/eps_mu;
      Xdown = 0;
      Poffset = 0;
      xScal = bisectForAnhyst(yNorm, Xdown, Xup, Poffset, eps_mu, tol); 
      xVal.Add(xScal,yDir);
      
      return xVal;
    }    
    
    /*
     * Check if yVal is beyond saturation > easy case
     */
    /*
     * NOTE: when we are in saturation (or if we assume to be for checking)
     * we have y and x aligned, i.e. we can check for this case the same way
     * as we do for scalar problems
     * However, we only have to check for positive saturation as we check the
     * amplitude which always should be >= 0
     */
    if(yNorm > 0){
      Vector<Double> yDir = yVal;
      yDir.ScalarDiv(yNorm);
      
      // get mu in current direction
      Vector<Double> tmp = Vector<Double>(dim_);
      Double eps_mu = 0.0;
      mu.Mult(yDir,tmp);
      yDir.Inner(tmp,eps_mu);
            
      Double anhystPartPosSat = YSaturated_*evalAnhystPart_normalized(1.0);
      if(yNorm >= (YSaturated_ + eps_mu*XSaturated_ + anhystPartPosSat) ){
        // saturation detected > solve like in 1d case via bisection
        traceMsg << "--B-- Inversion: Saturation found" << std::endl;
        Double xScal = 0.0;
        
        if(anhystPartPosSat == 0){
          traceMsg << "--B1-- Inversion: Anhysteretic part zero > solve by simple division" << std::endl;
          xScal = (yNorm - YSaturated_)/eps_mu;
          xVal.Init();
          xVal.Add(xScal,yDir);
          successCode = 2;
          return xVal;
        } else {
          traceMsg << "--B2-- Inversion: Anhysteretic non-zero > solve by bisection" << std::endl;
          Double tol = 1e-12;
          Double Xup, Xdown, Poffset;
          Xup = (yNorm - YSaturated_)/eps_mu;
          Xdown = XSaturated_;
          Poffset = YSaturated_;
          xScal = bisectForAnhyst(yNorm, Xdown, Xup, Poffset, eps_mu, tol);
          xVal.Init();
          xVal.Add(xScal,yDir);
          successCode = 2;
          return xVal;
        }
      }
    }
    
    /*
     * Use Levenberg-Marquart algorithm as presented in Dahmen&Reusken - Numerik partialler DFG
     */   
     // LM
    Vector<Double> hystVal;
    Vector<Double> hystVal_rem;
    /*
     * Check for remanence
     * > apparently, LM seems to work fine except for the remanence case
     * i.e. yval = hystval, xval = 0
     * > check if xVal = 0 would be a proper solution
     */
    xTMP.Init();
    hystVal_rem = computeValue_vec(xTMP, operatorIndex, overwriteMemory, overwriteDirection);
    
    diff.Init();
    diff.Add(1.0,yVal,-1.0,hystVal_rem);
    
    LOG_DBG(vecpreisachInversion) << "Check difference between yVal and hystVal(0); remanence?";
    LOG_DBG(vecpreisachInversion) << "Diff Vector: " << diff.ToString();
    LOG_DBG(vecpreisachInversion) << "Norm of diff: " << diff.NormL2();
    
    if(diff.NormL2() < INV_resTolB_){
      traceMsg << "--C-- Inversion: Remanence detected" << std::endl;
      xVal = xTMP;
      LOG_DBG(vecpreisachInversion) << "Set xVal to 0: " << xVal.ToString();
      successCode = 3;
      //	compCase = 2;
    } else {
      traceMsg << "--D-- Inversion: Use Levenberg Marquart" << std::endl;
      
      // tolerance wrt y > 1e-10 or 1e-12 seems good > takes 2-3 its
      // only problem: y-x-loops look ugly as x can be quite off!
      //Double tolError = 1e-11;  
      // tolerance for reevalution
      
      UInt itCnt = 0;
      Double alpha = INV_alphaLSStart_;
      Double alphaMin = INV_alphaLSMin_;//1.0/256.0;
      Double alphaMax = INV_alphaLSMax_;//8192;//512.0;
      
      /*
       * only residual wrt X in abs form tested
       */
      bool wrtX = true;
      bool relError = !true;
      
      if (diff.NormL2() < INV_resTolB_){
        // we are quite close to remanence, but have not yet reached it
        // start from 0
        traceMsg << "--D-- Inversion: Start at 0" << std::endl;
        xVal = xTMP;
        //  startAtZero=true;
        //	compCase = 3;
      } else {
        //	compCase = 4;
        xVal = prevXval;
        traceMsg << "--D-- Inversion: Start at previous X value = " << xVal.ToString() << std::endl;
        // different approach:
        // to this point we know:
        // old and new Y value
        // old X value
        // old H value
        // H value at 0 (H_rem)
        // > by this we can check:
        // a) newY > oldY > X should increase, too
        // b) newY < oldY < X should decrease, too
        //      b1) H_rem > newY > X = 0 is still too large but maybe good starting point
        //      b2) H_rem < newY > X between old X value and 0
        // > note: the steps above (sort of bisect) would work well for scalar model; for
        //         vector model we might have a change in direction but at similiar amplitude
        //         such that we do not match aboves cases
        //         still, this ssems more reasonable than fp iteration that nearly always drives
        //         x way above saturation due to large mu_inv
        traceMsg << "Try to obtain better starting value" << std::endl;
        if(yVal.NormL2() > prevYval.NormL2()){
          traceMsg << "Y increased in norm > increase X ,too" << std::endl;
          // increase X a bit (does not work if X is zero, though)
          xVal = prevXval;
          
          if(xVal.NormL2() != 0){
            Double factor;
            if(prevYval.NormL2() != 0){
              factor = yVal.NormL2()/prevYval.NormL2();
            } else {
              factor = 1.15;
            }
            traceMsg << "Scale X by " << factor << std::endl;
            xVal.ScalarMult(factor);
          } else {
            traceMsg << "X is zero; scaling will not help; try a scaled version of yVal instead" << std::endl;
            xVal = yVal;
            xVal.ScalarMult(xVal.NormL2()*XSaturated_/YSaturated_);
          }
        } else {
          traceMsg << "Y decreased in norm > decrease X ,too" << std::endl;
          //            if(hystVal_rem.NormL2() > yVal.NormL2()){
          //              LOG_DBG(vecpreisach) << "Y smaller than remanence (norm-wise) > start at remanence";
          //              // start at 0
          //              xVal.Init();
          //            } else {
          //              LOG_DBG(vecpreisach) << "Y larger than remanence but smaller than current hyst-value (norm-wise) > start in between";
          // take midpoint between remanence and currenc xval
          xVal = prevXval;
          xVal.ScalarMult(0.85);
          //            }
        }
        
      }
      traceMsg << "starting Xval: " << xVal.ToString() << std::endl;
      assert(!xVal.ContainsNaN() && !xVal.ContainsInf());
      /*
       * In the followng, prev*val_ = value of previous iteration
       *                      *val_ = current value
       */
      
      // perform some initial steps with fp iteation
      /*
       * Note: one iteration seems fine; two might already move solution so
       *       far off that it does not recover
       * > 1 works ok; 
       * (mostly due to bringing x into sat and then resetting it
       */
      //        ClipDirection(xVal);
      hystVal = computeValue_vec(xVal, operatorIndex, overwriteMemory, overwriteDirection);
      LOG_DBG(vecpreisachInversion) << "starting Hval: " << hystVal.ToString();
      //hystVal = prevHystval_[idElem];
      
      /*
       * As y is not in saturation (that case was already ruled out, 
       * the actual xVal cannot be in saturation neither
       * > It seems that initial fp iteration only works due to the following capping;
       * FP with the original mu nearly always put the system over saturation; the following
       * reset will cap it below saturation and from this point on the LM method seems
       * to converge much better than using the previous value of x
       * TODOL
       * > check if fp works better if xsat/ysat is used instead of mu
       * > check if lm works better if reset is done to a smaller percentage of xsat (currently 0.9)
       */
      if(xVal.NormL2() >= XSaturated_){
        traceMsg << "Reset xVal as its value is above Xsaturation" << std::endl;
        // reduce amplitude to a value slightly below saturation
        // remember: we may not take yVal directly as this might be larger than
        // ysaturated (due to y = hyst + mu*x
        // instead compare to yTMP (as above)
        //Double percent = 0.9; //yTMP.NormL2()/YSaturated_;
        
        /*
         * IDEA: take direction of new Y and scale it
         *  > only if y != 0!
         *  if y == 0, take xDirecttion
         * 
         */
        Vector<Double> dir = Vector<Double>(dim_);
        if(yNorm > 0){
          dir = yVal;
          dir.ScalarDiv(yNorm);
        } else {
          dir = xVal;
          dir.ScalarDiv(xVal.NormL2());
        }

        xVal = dir;
        xVal.ScalarMult(0.90*XSaturated_);
        
        //xVal.ScalarMult(percent*XSaturated_/xVal.NormL2());
        // or reset to 0? > 0.9 works good; 
        //xVal.Init();
        //resetAfterSat=true;
        //          ClipDirection(xVal);
        hystVal = computeValue_vec(xVal, operatorIndex, overwriteMemory, overwriteDirection);
      }
      
      Vector<Double> xStart = xVal;
      
      Double sign = 1.0;
      bool successError = false;
      bool successX = false;
      bool successY = false;
      bool discardUpdate = false;
      
      Vector<Double> xUpdate = Vector<Double>(dim_);
      Vector<Double> res = Vector<Double>(dim_);
      Matrix<Double> jac = Matrix<Double>(dim_,dim_);
      Matrix<Double> jacT = Matrix<Double>(dim_,dim_);
      
      Double errorNorm;
      Double errorNormResX;
      Double errorNormResY;
      UInt implementation = 0;
      Vector<Double> bestSol = Vector<Double>(dim_);
      Double bestErrorNorm = 1e16;
      
      Vector<Double> bestSolRes = Vector<Double>(dim_);
      Double bestErrorNormRes = 1e16;
      
      /*
       * Start actual LM iteration
       */
      while(true){ 
        itCnt++;
        totalNumberOfLMIterations++;
        
        if(debug){
          traceMsg << "--OUTER ITERATION-- (" << itCnt << ")" << std::endl;
          traceMsg << "Current solution: " << xVal.ToString() << std::endl;
        }
        // do not override hyst memory here
        //          ClipDirection(xVal);
        hystVal = computeValue_vec(xVal, operatorIndex, overwriteMemory, overwriteDirection);
        
        res = computeResidual(xVal,yVal,hystVal,mu,mu_inv,wrtX,relError);
        jac = computeJacobian(xVal,yVal,hystVal,res,mu,mu_inv,operatorIndex,
                sign,wrtX,relError,implementation,overwriteMemory,overwriteDirection);
        jac.Transpose(jacT);
        
        if(debug){
          traceMsg << "Current hystVal: " << hystVal.ToString() << std::endl;
          traceMsg << "Current Residual (wrt x): " << res.ToString() << std::endl;
          traceMsg << "Current startvalue for alpha_LS: " << alpha << std::endl;
          traceMsg << "Current JacobiMatrix: " << jac.ToString() << std::endl;
          
          Vector<Double> jacTres = Vector<Double>(dim_);
          jacT.Mult(res,jacTres);
          traceMsg << "jacTres: " << jacTres.ToString() << std::endl;
          traceMsg << "tolForErrorVector: " << INV_resTolH_ << std::endl;
          traceMsg << "tolForResidual X: " << INV_resTolH_ << std::endl;
          traceMsg << "tolForResidual Y: " << INV_resTolB_ << std::endl;
        }
        
        successError = checkConvergence(res,jacT,errorNorm,INV_resTolH_);
        errorNormResX = res.NormL2();
        successX = (errorNormResX <= INV_resTolH_);
        successY = checkInversionOutput(xVal, yVal, mu, INV_resTolB_, errorNormResY,
                operatorIndex, overwriteMemory, overwriteDirection);
        
        if(successError){
          // main criterion > would be best if this one could be satisfied
          traceMsg << "Success! Error estimate |jacT*ResX| = " << errorNorm << " < " << INV_resTolH_ << std::endl;
          successCode = 4;
          break;
        } else if(successX){
          // failback; still top if this works
          traceMsg << "Success! Residual norm wrt X = |ResX| = " << errorNormResX << " < " << INV_resTolH_ << std::endl;
          successCode = 5;
          break;
        } else {
          if( (totalNumberOfLMIterations%10 == 0)||(itCnt >= INV_maxIter_) ){
            if(successY){
              // failback; might still have large res error in x buz might be the best we can find
              // check only each 10th iteration
              traceMsg << "Success! Residual norm wrt Y = |ResY| = " << errorNormResY << " < " << INV_resTolB_ << std::endl;
              successCode = 6;
              break;
            } else if (itCnt >= INV_maxIter_) {
              if(debug){
                traceMsg << "Max number of iterations reached." << std::endl;
                traceMsg << "Last found solution xFound = " << xVal.ToString() << std::endl;
                traceMsg << "p(xFound) = " << hystVal.ToString() << std::endl;
                traceMsg << "Remaining residual wrt x-value (xFound - nu*(yIn - p(xFound)): " << res.ToString() << std::endl;
                traceMsg << "Remaining error norm |jacT*ResX|: " << errorNorm << std::endl;
                traceMsg << "Remaining residual-norm wrt x: " << errorNormResX << std::endl;
                traceMsg << "Remaining residual-norm wrt y: " << errorNormResY << std::endl;
              }
              //                  std::cout << "Remaining error norm |jacT*ResX|: " << errorNorm << std::endl;
              //                  std::cout << "Remaining residual-norm wrt x: " << errorNormResX << std::endl;
              //                  std::cout << "Remaining residual-norm wrt y: " << errorNormResY << std::endl;
              successCode = 0;
              
              break;
            }
          }
        }
        
        if(errorNormResX < bestErrorNormRes){
          bestErrorNormRes = errorNormResX;
          bestSolRes = xVal;                
        }
        
        if(errorNorm < bestErrorNorm){
          bestErrorNorm = errorNorm;
          bestSol = xVal;                
        }
        
        LOG_DBG(vecpreisachInversion) << "Solution not appropriate yet > Perform linesearch";
        UInt numberOfIterations = 0;
        
        discardUpdate = performLinesearch(xVal, yVal, res, xUpdate, jac, jacT, mu, mu_inv, 
                operatorIndex, overwriteMemory, overwriteDirection, alpha, alphaMin, alphaMax, wrtX, relError, numberOfIterations,xStart);
        
        if(alpha < minAlpha){
          minAlpha = alpha;
        }
        if(alpha > maxAlpha){
          maxAlpha = alpha;
        }
        avgAlpha += alpha;
        
        traceMsg << "Computed update: " << xUpdate.ToString() << std::endl;
        traceMsg << "Discard update? " << discardUpdate << std::endl;
        
        totalNumberOfLinesearchIterations += numberOfIterations;
        
        if(numberOfIterations > maximalNumberOfLinesearchIterations){
          maximalNumberOfLinesearchIterations = numberOfIterations;
        }
        
        /*
         * New idea: start with tighter bounds for alpha but shift these limits
         * after each 10 iterations
         */
        if(totalNumberOfLMIterations%10 == 0){
          alphaMin = alphaMin/2.0;
          alphaMax = alphaMax*2.0;
          // for Tikhonov > set xStart to xVal
          xStart = bestSol; // own test
        } 
        
        LOG_DBG(vecpreisachInversion) << "Computed update: " << xUpdate.ToString();
        if(!discardUpdate){
          xVal = xVal+xUpdate;
        } else {
          LOG_DBG(vecpreisachInversion) << "Discard update";
        }
        sign = sign*(-1.0);
      }
      
      avgAlpha /= totalNumberOfLMIterations;
      
      if(xVal.NormL2() >= XSaturated_){
        EXCEPTION("LM lead xVal into saturation > must not be the case!");
      }
    }
    // end LM
    
    /*
     * set values for next time
     */
    if(overwriteMemory == true){
      EXCEPTION("Memory should not be overridden here");
    }
    
    if(debug && successCode==0){
      Double resYNorm = 0.0;
      checkInversionOutput(xVal, yVal, mu, INV_resTolB_, resYNorm, operatorIndex, overwriteMemory, overwriteDirection,true);
      LOG_TRACE(vecpreisachInversion) << traceMsg.str();
    }
    
    LOG_TRACE(vecpreisachInversion) << " --------- END IVERSION --------- ";
    
    return xVal;
  } 
  
	// returns true, if yTarget can be reached by mu*xComputed + hyst(xComputed)
	bool VectorPreisachv10::checkInversionOutput(Vector<Double>& xComputed, Vector<Double>& yTarget, 
          Matrix<Double>& mu, Double tol, Double& resYNorm, Integer operatorIdx, bool overwriteMemory, bool overwriteDirection, bool output){
		
		Vector<Double> yCheck = Vector<Double>(dim_);
    Vector<Double> hCheck = computeValue_vec(xComputed, operatorIdx, overwriteMemory, overwriteDirection);
    
		mu.Mult(xComputed,yCheck);
		yCheck.Add(hCheck);
		
		Vector<Double> diff;
		diff = yCheck;
		diff -= yTarget;
    resYNorm = diff.NormL2();
    if(output){
      LOG_TRACE(vecpreisachInversion) << "CheckInversionOutput - norm of residual wrt y-vector: " << diff.NormL2();
    }
		if(diff.NormL2() < tol){
      if(output){
        LOG_TRACE(vecpreisachInversion) << "Solution OK; |diff| = " << diff.NormL2() << " < " << tol;
      }
			return true;
		} else {
      if(output){
        LOG_TRACE(vecpreisachInversion) << "Solution NOT OK; |diff| = " << diff.NormL2() << " > " << tol;
        LOG_TRACE(vecpreisachInversion) << "yRetrieved = " << yCheck.ToString();
        LOG_TRACE(vecpreisachInversion) << "yIn = " << yTarget.ToString();
      }
			return false;
		}
	}
	
  /*
   * MATRIX BASED IMPLEMENTATION
   */
  VectorPreisachv10_MatrixApproach::VectorPreisachv10_MatrixApproach(Integer numElem, Double xSat, Double ySat,
          Matrix<Double>& preisachWeight, Double rotationalResistance , UInt dim, bool isVirgin,
          bool classical, Double angularDistance, Double angResolution, Double anhystA, Double anhystB, Double anhystC, bool anhystOnly)
  : VectorPreisachv10(numElem, xSat, ySat,
          preisachWeight, rotationalResistance , dim, isVirgin,
          classical, angularDistance, angResolution, anhystA, anhystB, anhystC, anhystOnly)
  {
    LOG_TRACE(vecpreisach) << "Using Matrix-based implementation";
    
    /*
     * Get storage for switchingStates and rotatationStates
     */
    switchingStates_ = new Matrix<Double>[numElem];
    rotationStateX_ = new Matrix<Double>[numElem];
    rotationStateY_ = new Matrix<Double>[numElem];
    rotationStateZ_ = new Matrix<Double>[numElem];
    
    /*
     * Initialize arrays/vectors/matrices for each element
     */
    for(UInt k = 0; k < (UInt) numElem; k++){
      switchingStates_[k] = Matrix<Double>(numRows_,numRows_);
      
      rotationStateX_[k] = Matrix<Double>(numRows_,numRows_);
      rotationStateX_[k].Init();
      
      rotationStateY_[k] = Matrix<Double>(numRows_,numRows_);
      rotationStateY_[k].Init();
      
      if(dim_ == 3){
        rotationStateZ_[k] = Matrix<Double>(numRows_,numRows_);
        rotationStateZ_[k].Init();
      }
    }
    
    /*
     * Initialize switchingStates
     */
    for(UInt k = 0; k < numElem_; k++){
      InitializeSwitchingState(k);
    }
    
  }
  
  VectorPreisachv10_MatrixApproach::~VectorPreisachv10_MatrixApproach(){
    delete[] switchingStates_;
    delete[] rotationStateX_;
    delete[] rotationStateY_;
    delete[] rotationStateZ_;
    
    delete updateMatricesTimer_;
    delete evaluateMatricesTimer_;
    delete copyToTemporalStorageTimer_;
    delete copyFromTemporalStorageTimer_;
  }
  
  void VectorPreisachv10_MatrixApproach::InitializeSwitchingState(UInt idElem){
    
    /*
     * split Preisach plane along diagonal alpha=-beta in a +1 and -1 part
     */
    for(UInt i = 0; i < numRows_; i++){
      
      for(UInt j = 0; j <= i; j++){
        
        if(i+j+1 == numRows_){
          // on diagonal alpha = -beta -> leave at 0
        } else if(i+j+1 < numRows_){
          // below diagonal alpha = -beta -> +1
          switchingStates_[idElem][i][j] = +1;
        } else {
          // above diagonal alpha = -beta -> -1
          switchingStates_[idElem][i][j] = -1;
        }
        if (i == j){
          // on alpha = beta divide value by 2
          switchingStates_[idElem][i][j] = switchingStates_[idElem][i][j]/2.0;;
        }
      }
    }
  }
  
  std::string VectorPreisachv10_MatrixApproach::runtimeToString(){
    std::ostringstream oss;
    oss << "--- VectorPreisach MatrixApproach ---\n";
    if(performanceMeasurement_){
      Double totalUpdateTime = updateMatricesTimer_->GetCPUTime();
      Double totalEvaluationTime = evaluateMatricesTimer_->GetCPUTime();
      Double totalCopyingTimeForward = copyToTemporalStorageTimer_->GetCPUTime();
      Double totalCopyingTimeBackward = copyFromTemporalStorageTimer_->GetCPUTime();
      
      oss << "-- Updating -- \n";
      oss << "  Total number of update steps: " << updateMatricesCounter_ << "\n";
      oss << "  Average time to update matrices: " << totalUpdateTime/updateMatricesCounter_ << "\n";
      
      oss << "-- Evaluation -- \n";
      oss << "  Total number of evaluation steps: " << evaluateMatricesCounter_ << "\n";
      oss << "  Percentage of evaluation on temporal storage: " << copyToTemporalStorageCounter_/evaluateMatricesCounter_*100 << "% \n";
      oss << "  Average time to evaluate nested list: " << totalEvaluationTime/evaluateMatricesCounter_ << "\n"
              << "  Average time to create temporal copy: " << totalCopyingTimeForward/copyToTemporalStorageCounter_ << "\n"
              << "  Average time to restore backup: " << totalCopyingTimeBackward/copyToTemporalStorageCounter_ << "\n";
    }
    return oss.str();
  }
  
  void VectorPreisachv10_MatrixApproach::UpdateSwitchingStates(Vector<Double>& u_in, UInt idElem){
    /*
     * Update switching states from current input u_in; update only for FE element with index idElem
     * -> has to be called AFTER UpdateRotationStates
     */
    Vector<Double> curState = Vector<Double>(dim_);
    Double alpha,betaNext,xPar;
    
    /*
     * iterate over switching states
     */
    for(UInt i = 0; i < numRows_; i++){
      alpha = -1 + i*delta_;
      
      for(UInt j = 0; j <= i; j++){
        betaNext = -1 + (j+1)*delta_;
        
        /*
         * get corresponding rotation state
         */
        curState[0] = rotationStateX_[idElem][i][j];
        curState[1] = rotationStateY_[idElem][i][j];
        if(dim_ == 3){
          curState[2] = rotationStateZ_[idElem][i][j];
        }
        
        /*
         * calcualte xPar, which is responsible for the setting process
         */
        xPar = u_in.Inner(curState);
        
        xPar /= XSaturated_;
        
        /*
         * check if update is needed
         */
        Double fac = 1.0;
        if(i == j){
          fac = 0.5;
        }
        
        if(xPar > alpha){
          switchingStates_[idElem][i][j] = fac;
        } else if (xPar < betaNext){
          switchingStates_[idElem][i][j] = -fac;
        }
      }
    }
  }
  
  void VectorPreisachv10_MatrixApproach::UpdateRotationStates(Double XThres, Double xVal, Vector<Double>& e_u_new, UInt idElem){
    /*
     * Update rotation states from current input u_in = xVal * e_u_new; change of rotation state indicated by xThres
     * Update only for FE element with index idElem
     */
    
    Vector<Double> curState = Vector<Double>(dim_);
    Vector<Double> newState;
    Double alpha,betaNext;
    
    /*
     * iterate over rotation states
     */
    for(UInt i = 0; i < numRows_; i++){
      alpha = -1 + i*delta_;
      
      for(UInt j = 0; j <= i; j++){
        betaNext = -1 + (j+1)*delta_;
        
        /*
         * check if update is necessary
         */
        if(classical_){
          /*
           * rotation state changes, if xThres > alpha or -xThres < beta+delta_, with alpha and beta being the
           * bottom left coordinates of each element
           * -> skip setting of rotation state, if XThres < alpha AND -xThres > beta+delta
           */
          if((XThres < alpha)&&(-XThres > betaNext)){
            continue;
          }
        } else {
          /*
           * rotation state changes, if xThres > alpha AND -xThres < beta+delta_
           * -> skip setting of rotation state, if xThres < alpha OR -xThres > beta+delta_
           */
          if((XThres < alpha)||(-XThres > betaNext)){
            continue;
          }
        }
        
        curState.Init();
        
        /*
         * extract current rotation state
         */
        curState[0] = rotationStateX_[idElem][i][j];
        curState[1] = rotationStateY_[idElem][i][j];
        if(dim_ == 3){
          curState[2] = rotationStateZ_[idElem][i][j];
        }
        
        /*
         * compute new rotation direction
         */
        newState = evaluateNewRotationDirection(e_u_new, curState, xVal);
        //newState = clipNewRotationDirection(newState);
				if(restrictToHalfspace_){
					newState = restrictToHalfspace(newState);
        }
        /*
         * store new state
         */
        rotationStateX_[idElem][i][j] = newState[0];
        rotationStateY_[idElem][i][j] = newState[1];
        if(dim_ == 3){
          rotationStateZ_[idElem][i][j] = newState[2];
        }
      }
    }
  }
  
  Vector<Double> VectorPreisachv10_MatrixApproach::computeValue_vec(Vector<Double>& u_in, Integer idElem, bool overwrite, bool overwriteDirection, bool debugOut){

    /*
     * Determine the current rotational threshold
     *
     * Note: u_in is no longer normalized, i.e. u_in.NormL2 > Xsaturated is possible!
     * in that case, we have to cap u_in.NormL2 to Xsaturated BEFORE computing X_thres
     * Otherwise, the Preisach plane may be filled further than it is intended by the model
     * (e.g. for k = 0.5 and the revised model, the maximal filling state of the Preisach plane
     * is 0.5, not 1)
     */
    Double X_thres;
    Double uNormTmp = u_in.NormL2();
    if(uNormTmp >= XSaturated_){
      uNormTmp = 1.0;
    } else {
      uNormTmp = uNormTmp/XSaturated_;
    }
    if(classical_){
      X_thres = std::pow(uNormTmp,rotationalResistance_);
    } else {
      X_thres = uNormTmp*rotationalResistance_;
    }
    
    /*
     * Get current direction
     */
    Vector<Double> e_u = Vector<Double>(u_in.GetSize());
    Double xVal = u_in.NormL2();
    
    if(xVal > tol_){
      e_u = u_in/xVal;
    } else {
      e_u.Init(0.0);
    }
    
    // compute anhysteretic part first   
    Double anhystPartScal = evalAnhystPart_normalized(xVal/XSaturated_);
    Vector<Double> anhystPart = Vector<Double>(dim_);
    anhystPart.Init();
    anhystPart.Add(anhystPartScal,e_u);
    
    /*
     * Storage for return values
     */
    Vector<Double> retVec = Vector<Double>(dim_);
    retVec.Init();
    
    if(anhystOnly_ == false){
      Matrix<Double> switchingStatesSingleElemBAK;
      Matrix<Double> rotationStateXSingleElemBAK;
      Matrix<Double> rotationStateYSingleElemBAK;
      Matrix<Double> rotationStateZSingleElemBAK;
      
      if(overwrite == false){
        
        if(performanceMeasurement_){
          copyToTemporalStorageCounter_++;
          copyToTemporalStorageTimer_->Start();
        }
        
        /*
         * get copies of the data structures
         */
        switchingStatesSingleElemBAK = switchingStates_[idElem];
        rotationStateXSingleElemBAK = rotationStateX_[idElem];
        rotationStateYSingleElemBAK = rotationStateY_[idElem];
        rotationStateZSingleElemBAK = rotationStateZ_[idElem];
        
        if(performanceMeasurement_){
          copyToTemporalStorageTimer_->Stop();
        }
        
      }
      
      if(performanceMeasurement_){
        updateMatricesCounter_++;
        updateMatricesTimer_->Start();
      }
      
      /*
       * Update rotation states
       */
      if(overwriteDirection){
        UpdateRotationStates(X_thres, xVal, e_u, idElem);
      }
      
      /*
       * Update switching states
       */
      UpdateSwitchingStates(u_in, idElem);
      
      if(performanceMeasurement_){
        updateMatricesTimer_->Stop();
      }
  
      /*
       * Storage for element value
       */
      Double x = 0;
      Double y = 0;
      Double z = 0;
      Double s = 0;
      
      if(performanceMeasurement_){
        evaluateMatricesCounter_++;
        evaluateMatricesTimer_->Start();
      }
      
      /*
       * iterate over matrices and multiply entries together
       */
      for(UInt i = 0; i < numRows_; i++){
        
        for(UInt j = 0; j <= i; j++){
          s = switchingStates_[idElem][i][j];
          s *= preisachWeights_[i][j];
          x += s*rotationStateX_[idElem][i][j];
          y += s*rotationStateY_[idElem][i][j];
        }
      }
      
      retVec[0] = x;
      retVec[1] = y;
      
      if(dim_ == 3){
        for(UInt i = 0; i < numRows_; i++){
          
          for(UInt j = 0; j <= i; j++){
            s = switchingStates_[idElem][i][j];
            s *= preisachWeights_[i][j];
            z += s*rotationStateZ_[idElem][i][j];
          }
        }
        retVec[2] = z;
      }
      
      if(performanceMeasurement_){
        evaluateMatricesTimer_->Stop();
      }
      
      /*
       * scale retVec with delta_^2; then add anhyst part; finally scale with Ysaturated
       */
      retVec.ScalarMult(delta_*delta_);
      
      if(overwrite == false){
        /*
         * reload backup
         */
        if(performanceMeasurement_){
          copyFromTemporalStorageTimer_->Start();
        }
        
        switchingStates_[idElem] = switchingStatesSingleElemBAK;
        rotationStateX_[idElem] = rotationStateXSingleElemBAK;
        rotationStateY_[idElem] = rotationStateYSingleElemBAK;
        rotationStateZ_[idElem] = rotationStateZSingleElemBAK;
        
        if(performanceMeasurement_){
          copyFromTemporalStorageTimer_->Stop();
        }
      }
    }
    
    retVec.Add(anhystPart);
    retVec.ScalarMult(YSaturated_);
    
    if(overwrite == false){
      /*
       * store to tmp array
       */
      preisachSumTmp_[idElem] = retVec;
     
      return preisachSumTmp_[idElem];
    } else {
      preisachSum_[idElem] = retVec;
      
      prevXVal_[idElem] = u_in;
      prevHVal_[idElem] = preisachSum_[idElem];
      
      return preisachSum_[idElem];
    }
  }
  
  void VectorPreisachv10_MatrixApproach::switchingStateToBmp(UInt numPixel, std::string filename, UInt idElem, bool overLayWithRotState)
  {
    /*
     * NEW: rotation state is evaluated along with the switching states if overLayWithRotState is true
     * in the old versions, the rotation state was evaluated separately although we have to iterate over the
     * rotation list to evaluate the switching lists
     */
    
    if(numPixel < 2){
      WARN("Image should have more than 2 x 2 pixel");
      return;
    }
    
    if(numPixel%2 != 0){
      WARN("Rounded number of pixel ("<<numPixel<<") to a multiple of 2 ("<<numPixel+1<<")");
      numPixel = numPixel + 1;
    }
    
    /*
     * Calculate upscaling factor
     */
    UInt upscaling = (UInt) std::floor(numPixel/numRows_);
    
    /*
     * now call output function of matrix
     */
    
    if(overLayWithRotState == true){
      UInt version = 2;
      
      if(version == 1){
        
        for(UInt comp = 0; comp < dim_; comp++){
          
          std::stringstream stream;
          std::string filename_new;
          if(comp == 0){
            stream << "x-" << filename;
            filename_new = stream.str();
            switchingStates_[idElem].matrix2Bmp(upscaling,filename_new,&rotationStateX_[idElem]);
          } else if(comp == 1){
            stream << "y-" << filename;
            filename_new = stream.str();
            switchingStates_[idElem].matrix2Bmp(upscaling,filename_new,&rotationStateY_[idElem]);
          } else {
            stream << "z-" << filename;
            filename_new = stream.str();
            switchingStates_[idElem].matrix2Bmp(upscaling,filename_new,&rotationStateZ_[idElem]);
          }
        }
        
      } else if(version == 2) {
        /*
         * New way of outputting matrix:
         *   encode rotation state in color: 0 = red; 120 = blue; 240 = green; angles in between colored as a mix
         *   encode switching state as sign and amplitude: negative switching -> angle + 180; absvalue < 1 -> scale final colorcombination
         *
         *   -> only for 2D rotstates, as the z component is not considered
         */
        
        std::stringstream stream;
        stream << "xy-" << filename;
        std::string filename_new = stream.str();
        
        switchingStates_[idElem].matrix2Bmp_v2(upscaling,filename_new,&rotationStateX_[idElem],&rotationStateY_[idElem]);
      }
    } else {
      switchingStates_[idElem].matrix2Bmp(upscaling,filename,NULL);
    }
  }
  
  /*
   * LIST BASED IMPLEMENTATIONl
   */
  VectorPreisachv10_ListApproach::VectorPreisachv10_ListApproach(Integer numElem, Double xSat, Double ySat,
          Matrix<Double>& preisachWeight, Double rotationalResistance , UInt dim, bool isVirgin,
          bool classical, Double angularDistance, Double angResolution, Double anhystA, Double anhystB, Double anhystC, bool anhystOnly)
  : VectorPreisachv10(numElem, xSat, ySat,
          preisachWeight, rotationalResistance , dim, isVirgin,
          classical, angularDistance, angResolution, anhystA, anhystB, anhystC, anhystOnly)
  {
    
    LOG_TRACE(vecpreisach) << "Using List-based implementation";
    
    globRotList_ = new std::list<RotListEntryv10>[numElem_];
    for(UInt k = 0; k < numElem_; k++){
      globRotList_[k] = std::list<RotListEntryv10>();
      Initialize_GlobalRotationList(globRotList_[k]);
    }
    
    //      std::cout << "Global rotlists after initialization: " << std::endl;
    //      std::list<RotListEntryv10>::iterator listIt;
    //      for(UInt k = 0; k < numElem_; k++){
    //        std::cout << "Element number: " << k << std::endl;
    //        for(listIt = globRotList_[k].begin(); listIt != globRotList_[k].end(); listIt++){
    //          std::cout << listIt->ToString() << std::endl;
    //        }
    //      }
    /*
     * lowerTriangleValue_ and lastEu_ only needed for classical_ model
     */
    lowerTriangleValue_ = 0.0;
    
    if(classical_){
      Evaluate_LowerTriangle();
    }
    
    lastEu_ = new Vector<Double>[numElem_];
    for(UInt k = 0; k < numElem_; k++){
      //lastEu_[k] = initDir;
      lastEu_[k] = Vector<Double>(dim_);
      lastEu_[k].Init();
    }
    
    /*
     * Timer and counter
     */
    updateNestedListCounter_ = 0;
    evaluateNestedListCounter_ = 0;
    copyToTemporalStorageCounter_ = 0;
    
    updateRotListTimer_ = new Timer();
    updateSwitchingListTimer_ = new Timer();
    simplifyRotListTimer_ = new Timer();
    simplifySwitchingListTimer_ = new Timer();
    evaluateNestedListTimer_ = new Timer();
    copyToTemporalStorageTimer_ = new Timer();
    
  }
  
  VectorPreisachv10_ListApproach::~VectorPreisachv10_ListApproach(){
    delete[] globRotList_;
    delete[] lastEu_;
    
    delete updateRotListTimer_;
    delete updateSwitchingListTimer_;
    delete simplifyRotListTimer_;
    delete simplifySwitchingListTimer_;
    delete copyToTemporalStorageTimer_;
    delete evaluateNestedListTimer_;
  }
  
  std::string VectorPreisachv10_ListApproach::runtimeToString(){
    std::ostringstream oss;
    oss << "--- VectorPreisach ListApproach ---\n";
    if(performanceMeasurement_){
      Double totalRotListUpdateTime = updateRotListTimer_->GetCPUTime();
      Double totalSwitchListUpdateTime = updateSwitchingListTimer_->GetCPUTime();
      
      Double totalRotListSimplifyTime = simplifyRotListTimer_->GetCPUTime();
      Double totalSwitchListSimplifyTime = simplifySwitchingListTimer_->GetCPUTime();
      
      Double totalEvalNestedTime = evaluateNestedListTimer_->GetCPUTime();
      Double totalCopyingTime = copyToTemporalStorageTimer_->GetCPUTime();
      
      oss << "-- Updating -- \n";
      oss << "  Total number of nested list updates: " << updateNestedListCounter_ << "\n";
      oss << "  Average time to update outer rotation list: " << totalRotListUpdateTime/updateNestedListCounter_ << "\n"
              << "  Average time to update inner switching lists: " << totalSwitchListUpdateTime/updateNestedListCounter_ << "\n";
      
      oss << "-- Simplification/merging -- \n";
      oss << "  Total number of simplification steps: " << updateNestedListCounter_ << "\n";
      oss << "  Average time to simplify outer list: " << totalRotListSimplifyTime/updateNestedListCounter_ << "\n"
              << "  Average time to simplify inner switching lists: " << totalSwitchListSimplifyTime/updateNestedListCounter_ << "\n";
      
      oss << "-- Evaluation -- \n";
      oss << "  Total number of evaluation steps: " << evaluateNestedListCounter_ << "\n";
      oss << "  Percentage of evaluation on temporal storage: " << ((Double) copyToTemporalStorageCounter_)/((Double) evaluateNestedListCounter_)*100 << "% \n";
      oss << "  Average time to evaluate nested list: " << totalEvalNestedTime/evaluateNestedListCounter_ << "\n"
              << "  Average time to create temporal copy: " << totalCopyingTime/copyToTemporalStorageCounter_ << "\n";
    }
    return oss.str();
  }
  
  
  void VectorPreisachv10_ListApproach::Initialize_GlobalRotationList(std::list<RotListEntryv10>& usedList){
    /*
     * Make sure that list is empty
     */
    usedList.clear();
    
    /*
     * Initialize list with a maximum of value 0 and 0-vector
     */
    Vector<Double> zeroVec = Vector<Double>(dim_);
    zeroVec.Init(0);
    
    /*
     * Create new entry for rotlist; the contained switching list is empty
     */
    std::list<ListEntryv10> switchingList = std::list<ListEntryv10>();
    RotListEntryv10 initEntry = RotListEntryv10(0.0, 0.0, zeroVec, switchingList, 0.0, false, false);
    
    /*
     * For classical model, we insert a minimum of value 0
     * -> evaluation will lead to an upper triangle which is completely -1 and a split upper square
     *    together with the +1 lower triangle we get a total state of 0 (assuming symmetric weights!
     */
    if(classical_ == true){
      initEntry.getListReference().push_back(ListEntryv10(0,true,false));
    }
    
    usedList.push_back(initEntry);
    //lastXpar_[idElem] = 0.0;
  }
  
  void VectorPreisachv10_ListApproach::Initialize_GlobalRotationListWithValues(std::list<RotListEntryv10>& usedList,Vector<Double>& initDir, Double initRotValue, Double initSwitchValue){
    /*
     * Make sure that list is empty
     */
    usedList.clear();
    
    /*
     * clamp initValue to range 0,1
     */
    if(initRotValue < 0){
      initRotValue = 0.0;
    } else if (initRotValue > 1.0){
      initRotValue = 1.0;
    }
    
    /*
     * Create new entry for rotlist; the contained switching list is empty
     */
    std::list<ListEntryv10> switchingList = std::list<ListEntryv10>();
    RotListEntryv10 initEntry = RotListEntryv10(initRotValue, 0.0, initDir, switchingList, 0.0, false, false);
    
    /*
     * Note: switching list is set with entries of xPar which
     * computes from inner(u_in/Xsaturated_,rotDirection)
     * thus it should only have values between 0 and 1 (assuming that input is restricted to
     * amplitudes Xsaturated_ (even if this is not the case, it would make no difference as
     * the evaluation is limited to that range
     */
    if(initSwitchValue <= 0){
      initSwitchValue = 0.0;
    } else if(initSwitchValue >= XSaturated_){
      initSwitchValue = 1.0;
    } else {
      initSwitchValue = initSwitchValue/XSaturated_;
    }
    
    /*
     * For classical model, we (always!)insert a minimum of value 0
     * -> evaluation will lead to an upper triangle which is completely -1 and a split upper square
     *    together with the +1 lower triangle we get a total state of 0 (assuming symmetric weights!
     */
    if(classical_ == true){
      initEntry.getListReference().push_back(ListEntryv10(0.0,true,false));
    }
    
    if(initSwitchValue != 0){
      /*
       * we do not need to insert a value of 0; this would be the empty case;
       * otherwise, we insert the value as a maximum and ADDITIONALLY
       * 0.0 as minimum (i.e. we assume that material is in a remanence case)
       */
      initEntry.getListReference().push_back(ListEntryv10(initSwitchValue,false,false));
      initEntry.getListReference().push_back(ListEntryv10(0.0,true,false));
    }
    
    usedList.push_back(initEntry);
    //lastXpar_[idElem] = 0.0;
  }
  
  
  void VectorPreisachv10_ListApproach::Update_GlobalRotationList(Double xThres, Double xVal, Vector<Double> e_u, std::list<RotListEntryv10>& usedList,bool overwriteDirection,bool debugOut){
    /*
     * function for updating the global rotation list with an entry pair xThres,e_u
     * furthermore, the magnitude xVal of the original input vector u_in is passed to update the switching lists
     *
     * Remarks and note compared to older version (which used to have a list for each element of the Preisach plane)
     * 1. the global rotation list describes only the behavior in the square region 0 <= alpha <= 1; -1 <= beta <= 0
     * 2. the global rotation list contains only min or max instead of pairs
     * 3. each entry of the rotation list stores a switching list
     *  -> rotation list entries cannot be wiped out as easily as normally, as we may not loose the switching lists
     *     (remember: in the original model, switching states and rotation states are stored independently and such
     *     switching states must be able to remain even if the rotation state is overwritten)
     *  -> solution approach:
     *      1. check which (if any) older entries would get overwritten
     *      2. insert new entry at the right spot
     *      3. adapt the lower value (defining the lower bound of the area) of the partially overwritten entry to xThres
     *      4. set the rotation state of all completely overwritten entries to e_u
     *          -> DO NOT DELETE ANYTHING
     *      5. update the switching lists of all entries based on the new rotation state
     *      6. use Simplify_GlobalRotationList to merge adjacent rotListEntries if possible
     */
    
    if(usedList.empty()){
      EXCEPTION("List may not be empty at this point");
    }
    
    /*
     * Step 1. Update rotation list by inserting new entries (if needed!)
     */
    std::list<RotListEntryv10>::iterator listIt;
    /*
     *  Note that the function list::insert adds the new entry before position listIt
     * (i.e. if listIt = list.end() it will be inserted before the end)
     */
    std::list<RotListEntryv10>::iterator insertPos = usedList.end();
    std::list<RotListEntryv10>::iterator listStart = usedList.begin();
    std::list<RotListEntryv10>::iterator listEnd = --(usedList.end());
    
    Double curVal;
    /*
     * lowerBound is important for RotListEntries
     * the smallest possible entry (the one which is take if an entry is appended to the end of the list)
     * is
     */
    Double lowerBound = 0;
    /*
     * Prevstate: rotation state of previous list entry
     * e_u_old: rotation state of current list entry but from last timestep
     */
    Vector<Double> prevState;
    Vector<Double> e_phi = Vector<Double>(dim_);
    Vector<Double> e_u_old;
    bool posFound = false;
    bool needsInsert = true;
    
    /*
     * cut down xThres to range of interest (0 to 1)
     * (as xThres is used to compare against alpha and alpha in the upper square goes from 0 to 1)
     *
     */
    if(xThres > 1.0){
      xThres = 1.0;
    }
    
    // for classical_ we have to add entries for xThres = 0.0 as this entry would influence
    // the result of the upper triangle!
    // -> already done during Initialize_GlobalRotationList -> no further insert needed
    if(xThres <= 0.0){
      xThres = 0.0;
      needsInsert = false;
    }
    
    //    std::cout << "XThres: " << xThres << std::endl;
    //    std::cout << "e_u: " << e_u.ToString() << std::endl;
    //
    //    std::cout << "##########################" << std::endl;
    //    std::cout << "GlobalRotationList pre insert" << std::endl;
    //
    //    for(listIt = usedList.begin(); listIt != usedList.end(); listIt++){
    //    std::cout << listIt->ToString() << std::endl;
    //    }
    
    if(performanceMeasurement_){
      updateRotListTimer_->Start();
      updateNestedListCounter_++;
    }
    bool listUpdated = false;
    int cntInner = 0;
    
    //	std::cout << "Overwrite direction? " << overwriteDirection << std::endl;
    if((overwriteDirection)&&(needsInsert == true)){
      /*
       * Update rotation states
       */
      for(listIt = usedList.begin(); listIt != usedList.end(); listIt++){
        cntInner++;
        curVal = listIt->getVal();
        
        /*
         * check if new input value (xThres) is larger than curVal
         */
        if((xThres >= curVal)&&(posFound == false)){
          posFound = true;
          
          if(xThres == curVal){
            /*
             * here we do not need to insert a new entry as the area
             * defined by xThres and nextEntry->getVal does not change
             * -> simply overwrite the rotation state
             */
            needsInsert = false;
          } else {
            /*
             * in that case, we would have to insert a new entry
             * however, we should check the rotation state of the previous entry first
             * -> if that rotation state is equal to e_u, the new entry would only define a
             * subarea with same rotation state
             * -> such subareas can be simplified (i.e. merged) later, but it is no bad idea to prohibit such
             * unnecessary inclusions
             *
             * Note: we still have to overwrite the rotation state of all following entries
             */
            if(listIt != listStart){
              listIt--;
              prevState = listIt->getVecReference();
              listIt++;

              needsInsert = false;
              /*
               * check if previous state has the same (or nearly) the same rotation state
               * in that case, we do not have to insert this entry, as it would be merged with
               * the previous one later on
               */
              for(UInt i = 0; i< e_u.GetSize();i++){
                if(abs(e_u[i]-prevState[i])>tol_){
                  needsInsert = true;
                  break;
                }
              }
            } else {
              needsInsert = true;
            }
            
            if(needsInsert == true){
              /*
               * mark position for later
               * Note that the function list::insert adds the new entry before position listIt
               * (i.e. if listIt = list.end() it will be inserted before the end)
               */
              insertPos = listIt;
              
              // ???
              //            if(insertPos == globRotList_[idElem].end()){
              //
              //            }
              
              /*
               * the new rotation area will be between xThresh and curVal
               */
              lowerBound = curVal;
            }
          }
        } // xThres >= curVal
        if(posFound == true){
          listUpdated = true;
          /*
           * we already have found a suitable position, i.e. all following entries would be
           * overwritten by the new entry; as we want to keep the switching list, we simply overwrite
           * the rotation state
           * (if the rotation state is already e_u, the flag rotHasChanged_ is not set to true!)
           */
          
          if(classical_){
            /*
             * classical model knows no angular distance -> full rotation is performed
             */
						if(restrictToHalfspace_){
							e_phi = restrictToHalfspace(e_phi);
						}
            listIt->setVec(e_phi);
          } else {
            /*
             * here we do not set the state to e_u directly but rotate it towards e_u;
             * Note that each previously stored rotation state gets rotated to a different e_phi!
             */
            e_u_old = listIt->getVecReference();
            e_phi = evaluateNewRotationDirection(e_u,e_u_old,xVal);
            //e_phi = clipNewRotationDirection(e_phi);
						// > clipping moved to ceofFunction hyst > input and output gets clipped and
						//    correct value are stored
						if(restrictToHalfspace_){
							e_phi = restrictToHalfspace(e_phi);
						}
            
            listIt->setVec(e_phi);
          }
        }
      } // loop
      
      if(posFound == false){
        /*
         * no suitable position was found (i.e. the current input is the smallest so far)
         * -> append to end of the list, but only if the previous entry does not have the same rotation state
         * (in that case we would create two adjacent areas with same rotation state; would be no big deal as
         * we later call the simplify list function, but this saves some unnecessary computations)
         */
        prevState = listEnd->getVecReference();
        
        needsInsert = false;
        
        Vector<Double> e_u_check;
        if(restrictToHalfspace_){
          e_u_check = restrictToHalfspace(e_u);
        } else {
          e_u_check = e_u;
        }
        
        for(UInt i = 0; i< e_u.GetSize();i++){
          if(abs(e_u_check[i]-prevState[i])>tol_){
            needsInsert = true;
            /*
             * insert pos will already be pointing to the end of the list in this case (see definition above)
             */
            break;
          }
        }
      }
      
      /*
       * finally, if we still need to insert an element -> do it
       *
       *  Note that the function list::insert adds the new entry before position listIt
       * (i.e. if listIt = list.end() it will be inserted before the end)
       */
      if(needsInsert == true){
        
        /*
         * Important notes:
         * 1. new entries inherit the switching list from the entry, which they (partially) overlay
         *    this is the entry prior to insertPos (if any!); otherwise the new list will be empty
         * 2. the previous entry (if any) gets a new lower bound (the current xThres)
         * 3. for Sutor2015: rotation state of new entry has to be created from rotating the direction of the (partially)
         *    overlaid entry towards the current direction
         */
        
        std::list<ListEntryv10> newList = std::list<ListEntryv10>();
        Double lastXpar = 0.0;
        UInt startCnt;
        bool wasWipedOut;
        if(insertPos != listStart){
          /*
           * get previous entry
           */
          insertPos--;
          newList = insertPos->getListCopy();
          lastXpar = insertPos->getLastLocalXpar();
          /*
           * e_u_old is the rotation state which will rotate towards e_u_new
           */
          e_u_old = insertPos->getVecReference();
          
          /*
           * as we inherit the switching list, we also inherit the wiped out property and the value of startCnt
           */
          wasWipedOut = insertPos->wasListWipedOut();
          startCnt = insertPos->getStartCnt();
          /*
           * due to the new entry, the previous one will be reduced in size
           * -> set lowerVal of that element (this will also set the flag isChanged_ of the element to true)
           */
          insertPos->setLowerVal(xThres);
          insertPos++;
        } else {
          if(classical_){
            /*
             * add minimum of value 0 -> the same as in initialize_globalRotationList
             * No dummy anymore!
             */
            newList.push_back(ListEntryv10(0,true,false));
            wasWipedOut = false;
            startCnt = 0;
          }
          /*
           * new previous rotation state available; e_u_old = 0-vector
           */
          e_u_old = Vector<Double>(dim_);
        }
        
        /*
         * Important: do not forget to insert rotated from evaluateNewRotationDirection instead of e_u!
         * As we insert at the end of the list, we overwrite (at least partially) the last rotentry
         * -> needed rotation direction of partially overlapped rotation state -> see above
         */
        e_phi = evaluateNewRotationDirection(e_u,e_u_old,xVal);
        //e_phi = clipNewRotationDirection(e_phi);
				if(restrictToHalfspace_){
					e_phi = restrictToHalfspace(e_phi);
				}
				
        usedList.insert(insertPos,RotListEntryv10(xThres,lowerBound,e_phi,newList,lastXpar,false,false,wasWipedOut,startCnt));
        listUpdated = true;
      }
    } // if overwriteDirection = true
    
    if(performanceMeasurement_){
      updateRotListTimer_->Stop();
      updateSwitchingListTimer_->Start();
    }
    
    /*
     * Step 2. For each entry in rotation list, update the switching list with the value xPar
     * -> this has to be done even if overwriteDirection = false
     */
    Double xPar;
    UInt updated;
    Vector<Double> rotState;
    /*
     * reset list iterator to the end; list may be extended!
     */
    listEnd = --(usedList.end());
    
    bool anySwitchingListUpdated = false;
    for(listIt = usedList.begin(); listIt != usedList.end(); listIt++){
      
      rotState = listIt->getVecReference();
      
      xPar = rotState.Inner(e_u)*xVal;// = rotState.Inner(u_in);
      /*
       * Normalize to Xsaturated
       */
      xPar /= XSaturated_;
      // do not set to zero!
      // this will lead to wrong results
//      if( abs(xPar) < 1e-15 ){
//        xPar = 0.0;
//      }
      
      /*
       * we need to pass lastXpar to list
       */
      bool isLastRotEntry = false;
      if(listIt == listEnd){
        /*
         * last entry of rotation list -> switching state for classical evaluation can extend into upper triangle area!
         *
         */
        isLastRotEntry = true;
      }
      
      Rectangle bbox = Rectangle(0,0,0,0);
      getBoundingBoxFromRotEntry(listIt, bbox, isLastRotEntry);
      
      updated = Update_SwitchingList(listIt->getListReference(),xPar,listIt->getLastLocalXpar(), bbox, listIt->wasListWipedOut(),isLastRotEntry);
      
      if(updated != 0){
        anySwitchingListUpdated = true;
      }
      
      //      if(updated != 0){
      //        std::cout << "UPDATE of switching list!" << std::endl;
      //      } else {
      //        std::cout << "NO UPDATE of switching list!" << std::endl;
      //      }
      
      if(updated == 1){
        /*
         * list was updated by the input was not strong enough to wipe out the diagonal splitting along alpha = -beta
         */
        listIt->setLastLocalXpar(xPar);
      } else if(updated == 2){
        /*
         * list was updated by the input and the diagonally split part was wiped out
         */
        listIt->setLastLocalXpar(xPar);
        listIt->setWipedOut(true);
        /*
         * reset startCnt (list was wiped out it was completely overwritten too
         * In that case set startCnt back to 0
         * (for exmplanation of startCnt see class RotListEntryv10)
         */
        listIt->setStartCnt(0);
      } else if(updated == 3){
        /*
         * list was updated by the input; no wipe out (i.e. still partially split by alpha = -beta) but first
         * element was replaced!
         */
        listIt->setLastLocalXpar(xPar);
        /*
         * reset startCnt
         */
        listIt->setStartCnt(0);
      }
      
      /*
       * else: no update was performed!
       */
    }
    
    //    std::cout << "##########################" << std::endl;
    //    std::cout << "GlobalRotationList pre merge and simplify" << std::endl;
    //
    //    for(listIt = usedList.begin(); listIt != usedList.end(); listIt++){
    //    std::cout << listIt->ToString() << std::endl;
    //    }
    
    if(performanceMeasurement_){
      updateSwitchingListTimer_->Stop();
      simplifyRotListTimer_->Start();
    }
    
    if(debugOut){
      LOG_TRACE(vecpreisachInversion) << "tmpList pre simplify of rotlist";
      std::list<RotListEntryv10>::iterator listIt;
      for(listIt = usedList.begin(); listIt != usedList.end(); listIt++){
        LOG_TRACE(vecpreisachInversion) << listIt->ToString();
      }
      LOG_TRACE(vecpreisachInversion) << "##########################";
    }	
    
    /*
     * Step 3. Merge adjacent rotList entries
     */
    bool rotListSimplified  = false;
    rotListSimplified = Simplify_GlobalRotationList(usedList);
    
    if(debugOut){
      LOG_TRACE(vecpreisachInversion) << "tmpList post simplify of rotlist";
      std::list<RotListEntryv10>::iterator listIt;
      for(listIt = usedList.begin(); listIt != usedList.end(); listIt++){
        LOG_TRACE(vecpreisachInversion) << listIt->ToString();
      }
      LOG_TRACE(vecpreisachInversion) << "##########################";
    }	
    
    
    if(performanceMeasurement_){
      simplifyRotListTimer_->Stop();
      simplifySwitchingListTimer_->Start();
    }
    
    /*
     * Step 4. Simplify local switching lists
     */
    bool anySwitchingListSimplified = Simplify_LocalSwitchingLists(usedList);
    
    if(debugOut){
      LOG_TRACE(vecpreisachInversion) << "tmpList post simplify of switchlist";
      std::list<RotListEntryv10>::iterator listIt;
      for(listIt = usedList.begin(); listIt != usedList.end(); listIt++){
        LOG_TRACE(vecpreisachInversion) << listIt->ToString();
      }
      LOG_TRACE(vecpreisachInversion) << "##########################";
    }	
    
    
    if(performanceMeasurement_){
      simplifySwitchingListTimer_->Stop();
    }
    
    if(debugOut){
      LOG_TRACE(vecpreisachInversion) << "GlobalRotlist updated? " << listUpdated;
      LOG_TRACE(vecpreisachInversion) << "Any switching list udpdated? " << anySwitchingListUpdated;
      LOG_TRACE(vecpreisachInversion) << "GlobalRotList simplified? " << rotListSimplified;
      LOG_TRACE(vecpreisachInversion) << "Any switching list simplified? " << anySwitchingListSimplified;
    }
    
    
    //    static int cnt = 0;
    //
    //    std::cout << "GlobalRotationList after merge and simplify -- step " << ++cnt << std::endl;
    //
    //    for(listIt = usedList.begin(); listIt != usedList.end(); listIt++){
    //    std::cout << listIt->ToString() << std::endl;
    //    }
    //    std::cout << "##########################" << std::endl;
    
  }
  
  UInt VectorPreisachv10_ListApproach::Update_SwitchingList(std::list<ListEntryv10>& list, Double newEntry, Double lastXpar, Rectangle boundingBox, bool wasWipedOut, bool lastRotEntry){
    /*
     * This function is used to update both the globalSwitching list as well as the
     * local switching list stored in each entry of the rotation list
     *
     * We can basically reuse most parts of the older UpdateList function (which used to update lists for
     * each element alpha,beta). Just consider the upper square (where all the local switching lists are working)
     * as one large element with alpha = 0, beta = -1 and the upper triangle (where the global switching list is
     * set) to be another large element with alpha = 0 and beta = 0; The size of these elements is of course not
     * delta_^2 anymore but 1^2.
     *
     *  NEW VERSION (7.8.2016):
     *    Current (v7) wiping out rule:
     *      1. min override all larger minima and all maxima which followed the deleted minima
     *      2. max override all smaller maxima and all minima which followed the deleted maxima
     *      3. minima cannot delete maxima and maxima cannot delete minima DIRECTLY
     *
     *    Resulting state of min/max list:
     *      1. list will be alternating list of minima and maxima
     *      2. maxima are descending; minima are ascending
     *      3. overall amplitude does not have to be decreasing, i.e. list could
     *          be 1,-0.4,0.8,-0.2 but not 0.8,-0.4,1,-0.2
     *
     *    Issues with this rule:
     *      Asuume that a list begins with minimum/maximum and a new maximum/minimum with larger absolute value
     *      is added --- like (-0.6,0.4) add (0.8) -> (-0.6,0.8) ---. If we now replace the first entry of the list
     *      according to aboves wiping out rule --- e.g. (-0.6,0.8) add (-0.7) -> (-0.7) ---, we might
     *      loose important information.
     *
     *                  -0.7                -0.7
     *               ____v____           ____v____
     *              |\   !     |         |\   !    |
     *              |__\_!_____|0.8      |  \_!    |-> note the step! it cannot recreated without 0.8
     *              |    \     |      -> |   ^!    |
     *              |    ! \   |         |    !    |
     *              |____!____\|         |____!____|
     *                     -0.6
     *
     *      Currently (v<8) we checked for cases as the one above and kept important maxima/minima if needed.
     *      E.g. we would have (-0.7,0.8) then. This is, however, not the only problem. During evaluation,
     *      a list like (-0.7,0.8) will not work with the used evaluation scheme, as that scheme assumes
     *      the first entry of the list to have the largest absolute value.
     *      Take the following example:
     *
     *           MinMax[1] = min
     *       ________v________
     *      | \      |    |   |
     *      |___\____|____|___| MinMax[0] = max
     *      |     \  |    |   |
     *      |_______\|____|___| MinMax[2] = max
     *      |        |\   |   |
     *      |________|__\_|___| MinMax[4] = max
     *      |        |    \   |
     *      |________|____| \_|
     *              MinMax[3] = min
     *
     *      Here, the first entry abs(MinMax[0]) is larger than abs(MinMax[1]) and thus allows
     *      to decompose the square into subareas like
     *
     *       area 0 (square!)
     *       _v_______________
     *      | \  | 1 |    |   |
     *      |___\|___| 3  |   |
     *      |   2    |    |5  |
     *      |________|____|   |
     *      |     4       |   | < std decomposition if starting with maximum
     *      |_____________|___|
     *      |        6        |
     *      |_________________|
     *
     *      If abs(MinMax[1]) is larger than abs(MinMax[0]), however, we get the following issues
     *
     *        here is MinMax[1] = min
     *       __v______________
     *      | \|1|        |   |
     *area0>|__|\|     3  |   |
     *no    |2 |          |5  |
     *square|__|__________|   | < std decomposition does not work
     *      |     4       |   |
     *      |_____________|___|
     *      |        6        |
     *      |_________________|
     *
     *      The decomposition above cannot be treated by the used evaluation anymore as
     *        - area0 is not square anymore (which is important as we than can leave it out assuming
     *          symmetric weights)
     *        - area1 is not completely -1
     *        - area2 is too small
     *        - area3 is no rectangle
     *      This issue can be solved by checking for the special case and than tweaking the decomposition.
     *
     *      However, there seems to be a much easier way.
     *
     *    New, simpler (original?) approach:
     *      Extend wiping out rule by a forth point:
     *        4. if abs(input) is larger than abs(first list entry) -> wipe out list completely and set input
     *           as first value
     *           -> a minimum may override a maximum and vice versa, but only in that special case
     *
     *       With respect to aboves list-example, we get:
     *        (-0.6,0.4) add (0.8) -> (0.8)
     *        (0.8) add (-0.7) -> (0.8,-0.7)
     *            *
     *                  -0.7                -0.7
     *               ____v____           ____v____
     *              |\   !     |         |\   !    |
     *              |__\_!_____|0.8      |__\_!____|0.8
     *              |    \     |      -> |   ^!    |       < state can now be recreated
     *              |    ! \   |         |    !    |
     *              |____!____\|         |____!____|
     *                     -0.6
     *
     *      And regarding the evaluation scheme:
     *
     *         original MinMax[0]      the former MinMax[1] is now MinMax[0]
     *        ____v____________          __v______________
     *       | \|1|        |   |        |_\|          |   |
     *       |__|\|     3  |   |        |  |     2    |   |
     *       |2 |          |5  |        |1 |          |4  |
     *       |__|__________|   |  ->    |__|__________|   | < std scheme for a list starting
     *       |     4       |   |        |     3       |   |   with a minimum!
     *       |_____________|___|        |_____________|___|
     *       |        6        |        |        5        |
     *       |_________________|        |_________________|
     *
     *      Note that the former starting value MinMax[0] is not needed for the transformed state.
     *
     *    IMPORTANT REMARKS:
     *      1. The new update/wiping-out rule is only relevant for the FIRST ENTRY of the list!
     *      2. The new treatment is only for the upper square!
     *          Reason: the upper triangle is not symmetric to alpha = -beta! Here it does not work
     *                and makes problems instead!
     *          Note: in the global sense it would be ok, but as we split along beta = 0 we have to
     *                distinguish.
     *
     *
     */
    
    // not needed anymore as we clip against bounding box
    //    Double alpha,beta,delta;
    //
    //    /*
    //     * New treatment: switching state is always over the whole Preisach plane; it will later get clipped to
    //     * fitting rotation states
    //     * in that case, the min-max list can have values from -1 to +1 for alpha and beta
    //     */
    //    alpha = -1.0;
    //    beta = -1.0;
    //    delta = 2.0;
    
    std::list<ListEntryv10>::iterator listIt;
    
    bool appendDirectly = false;
    if(list.empty()==true){
      appendDirectly = true;
    }
    
    int state;
    /*
     * check if the current input is a minimum or a maximum
     */
    if(abs(newEntry-lastXpar)<tol_){
      /*
       * no update -> return 0
       */
      return 0;
    } else if(newEntry > lastXpar){
      /*
       * maximum found -> 2
       */
      state = 2;
    } else {
      state = -2;
    }
    
    /*
     * get outer boundings of the corresponding rotation state
     */
    Double l,r,t,b;
    boundingBox.getBounds(l,r,t,b);
    
    //    std::cout << "xPar: " << newEntry << std::endl;
    //    std::cout << "lastXpar: " << lastXpar << std::endl;
    //    std::cout << "left/right: " << l << " / " << r << std::endl;
    //    std::cout << "isLastRotEntry: " << lastRotEntry << std::endl;
    
    /*
     * Check if lists gets completely overwritten (i.e. the saturated input value exceeds the bounds of the region
     */
    if(state > 0){
      /*
       * Element is to be set by maximum
       *
       * NEW: restrict to bounding box of rotation entry!
       * Background: Assume a rotation entry inside the bounding box l,r,t,b
       *  __ __ __ _t __ __ __ __
       * |               |       |
       * |               |       |
       * |               |       |
       * |               |       |
       * l               |       r
       * |               |       |
       * |__ __ __ __ __ |       |
       * |       rotstate        |
       * |                       |
       * |__ __ __ _b __ __ __ __|
       *
       * Stairs of the switching state contribute only, if the lines pass through the bounding box
       * -> if a line lies above t, it makes no difference if its value is >t or exactly t; in both cases the whole switching states
       *      is set to +1
       * -> in a similar way it behaves for r;
       * -> for l and b, it makes no difference, if input is smaller or exactly that value
       *
       *
       *  -> restrict to [b,t]
       */
      if(newEntry >= t){
        /*
         * current region will be set completely to +1
         * (note that in the case of local lists in the square, this applies only to the current rotation state!)
         */
        /*
         * empty list
         */
        list.clear();
        /*
         * cut value
         */
        newEntry = t;
        /*
         * add maximum to list
         */
        list.push_back(ListEntryv10(newEntry,false));
        /*
         * elements on alpha = -beta also loose their special state
         */
        return 2;
        
      } else if(newEntry <= b){
        /*
         * new entry will have no effect
         */
        return 0;
      } else if(appendDirectly){
        list.push_back(ListEntryv10(newEntry,false));
        /*
         * return code 3: list was not wiped, but first entry was set/replaced
         */
        return 3;
      } else if(classical_){
        /*
         * Works only on the L-shapes -> only for classical_
         */
        /*
         * for upper square, check value of lowerVal
         * if a maximum to be added has a value smaller than lowerVal (bottom of the rotation area), the resulting area will have
         * no effect at all onto the switching state
         * -> do not add to list
         * (similar to the check above with newEntry <= alpha)
         */
        if(newEntry <= b){
          return 0;
        }
      }
      /*
       * else -> see below
       */
    } else {
      /*
       * Element is to be set by a minimum -> restrict to [l,r]
       */
      if(newEntry < l){
        /*
         * Current element will be completely set to -1 (in case of MinMaxList), will be rotated into direction newVecEntry (in case of RotList)
         */
        /*
         * empty list
         */
        list.clear();
        /*
         * cut value
         */
        newEntry = l;
        /*
         * add minimum to list
         */
        list.push_back(ListEntryv10(newEntry,true));
        /*
         * elements on alpha = -beta also loose their special state
         */
        /*
         * return 2 -> list was updated and wiped
         */
        return 2;
        
      } else if(newEntry >= r){
        /*
         * new entry will have no effect
         */
        return 0;
      } else if(appendDirectly){
        list.push_back(ListEntryv10(newEntry,true));
        /*
         * return code 3: list was not wiped, but first entry was set/replaced
         */
        return 3;
      } else if(classical_){
        /*
         * Works only on the L-shapes -> only for classical_
         * -> make sure when calling this function, that it does set lower value to
         */
        /*
         * for upper square, check value of lowerVal
         * if a minimum to be added has a value larger than -lowerVal (right boundary), the resulting area will have
         * no effect at all onto the switching state
         * -> do not add to list
         * (similar to the check above with newEntry >= beta+delta)
         */
        if((newEntry >= r)&&(lastRotEntry == false)){
          /*
           * Note: Unlike the maxima (which are terminated by the lower bound or the beta axis),
           * the minima can NOW reach into the upper triangle region for the classical approach.
           * This is only the case for the last entry of the rot list. In that case we can allow minima
           * which do not stop at the lower value (which would be 0) but which can extend up to beta = +1
           */
          //NOTE: can only be equal r here, as the corresponding bounding box for the last entry goes up to +1
          return 0;
        }
      }
      /*
       * else -> see below
       */
    }
    
    /*
     * it is ok if list is empty at beginning of the function;
     * however, in that case, an entry is added to the list or the function
     * already returned
     */
    if(list.empty()==true){
      EXCEPTION("List is not allowed to be empty at this point!");
    }
    std::list<ListEntryv10>::iterator lastEntry = --(list.end());
    
    /*
     * if we came to this point, we have to iterate through the list and check if value shall be included
     */
    bool canBeInserted = false;
    bool canBeDeleted = false;
    
    /*
     * compare current input type with the extremum type of the last entry
     * if they are of opposite type (e.g. a maximum shall be inserted and last entry of list is a minimum)
     * the value can be added to the list regardless if a previous value is found which is to be replaced;
     * if they have the same extremum type, it cannot be appended to the list (in contrast to the case of rotList);
     * here we only insert if we find a fitting spot inside the list
     */
    if(state < 0){
      /*
       * new value is minimum
       */
      if(lastEntry->isMin() == false){
        canBeInserted = true;
      }
    } else if (state > 0){
      /*
       * new value is maximum
       */
      if(lastEntry->isMin() == true){
        canBeInserted = true;
      }
    }
    
    bool listMin;
    bool firstEntrySet = false;
    UInt cnt = 0;
    Double listVal;
    ListEntryv10 helperEntry = ListEntryv10(0,false);
    
    for(listIt = list.begin(); listIt != list.end(); listIt++){
      
      listMin = listIt->isMin();
      listVal = listIt->getVal();
      
      if(state == 2){
        /*
         * we iterate over MinMaxList and new entry is a maximum
         */
        if(!listMin){
          /*
           * compare value with value of list maximum
           */
          if(newEntry >= listVal){
            canBeInserted = true;
            canBeDeleted = true;
          }
        } else {
          /*
           * New in version 8
           * -> first entry of list can be overwritten by an extremum of different type if abs(newValue) >= abs(oldValue)!
           * Only valid for first entry!
           */
          if((cnt == 0)&&( (abs(newEntry)-abs(listVal)) >= 0)){
            canBeInserted = true;
            canBeDeleted = true;
          }
        }
        
      } else if(state == -2){
        /*
         * we iterate over MinMaxList and new entry is a minimum
         */
        if(listMin){
          /*
           * compare value with value of list minimum
           */
          if(newEntry <= listVal){
            canBeInserted = true;
            canBeDeleted = true;
          }
        } else {
          /*
           * New in version 8
           * -> first entry of list can be overwritten by an extremum of different type if abs(newValue) >= abs(oldValue)!
           * Only valid for first entry!
           */
          if((cnt == 0)&&( (abs(newEntry)-abs(listVal)) >= 0)){
            canBeInserted = true;
            canBeDeleted = true;
          }
        }
      }
      
      if(canBeDeleted){
        /*
         * delete all entries of the list starting with the current one
         */
        list.erase(listIt,list.end());
        
        /*
         * set flag denoting that the first entry of the list got deleted, too.
         * -> after this function, a new first entry is to be set
         * -> return 3 instead of return 1
         */
        if(cnt == 0){
          firstEntrySet = true;
        }
        
        break;
      }
      cnt++;
    }
    
    if(canBeInserted == true){
      
      if(state > 0){
        /*
         * insert maximum
         */
        list.push_back(ListEntryv10(newEntry,false));
      } else if(state < 0){
        /*
         * insert minimum
         */
        list.push_back(ListEntryv10(newEntry,true));
      }
      /*
       * mark element as updated
       */
      if(firstEntrySet == true){
        /*
         * return code 3: list was not wiped, but first entry was set/replaced
         * (only relevant for local swtiching lists!
         */
        return 3;
      } else {
        return 1;
      }
    }
    /*
     * if we end up here, we did not change the list -> return 0
     * (can be deleted would also change list, but that flag is only true, if canBeInserted is true, too)
     */
    return 0;
  }
  
  Vector<Double> VectorPreisachv10_ListApproach::computeValue_vec(Vector<Double>& u_in, Integer idElem, bool overwrite,bool overwriteDirection, bool debugOut){
    /*
     * Determine the current rotational threshold
     *
     * Note: u_in is no longer normalized, i.e. u_in.NormL2 > Xsaturated is possible!
     * in that case, we have to cap u_in.NormL2 to Xsaturated BEFORE computing X_thres
     * Otherwise, the Preisach plane may be filled further than it is intended by the model
     * (e.g. for k = 0.5 and the revised model, the maximal filling state of the Preisach plane
     * is 0.5, not 1)
     */
    Double X_thres;
    Double uNormTmp = u_in.NormL2();
        
    if(uNormTmp >= XSaturated_){
      uNormTmp = 1.0;
    } else {
      uNormTmp = uNormTmp/XSaturated_;
    }
    if(classical_){
      X_thres = std::pow(uNormTmp,rotationalResistance_);
    } else {
      X_thres = uNormTmp*rotationalResistance_;
    }
    
    /*
     * Get current direction
     */
    Vector<Double> e_u = Vector<Double>(u_in.GetSize());
    Double xVal = u_in.NormL2();
    
    //  std::cout << "xVal: " << xVal << std::endl;
    //	std::cout << "OverwriteDirection? " << overwriteDirection << std::endl;
    
    //if(xVal > tol_) //another tolerance?!
    if(xVal != 0){
      e_u = u_in/xVal;
    } else {
      // reuse old direction
      e_u = lastEu_[idElem];
      //e_u.Init(0.0);
    }
    
    // compute anhysteretic part first
    Double anhystPartScal = evalAnhystPart_normalized(xVal/XSaturated_);
    Vector<Double> anhystPart = Vector<Double>(dim_);
    anhystPart.Init();
    anhystPart.Add(anhystPartScal,e_u);
    
    // std::cout << "e_u: " << e_u.ToString() << std::endl;
    
    /*
     * set value of lastEu_ (only needed for classical_ model to get the rotation information for the lowerTriangle_)
     */
    if(overwriteDirection){
      lastEu_[idElem] = e_u;
    }
      
    /*
     * Storage for element value
     */
    Vector<Double> Yout = Vector<Double>(dim_);
    Yout.Init(0.0);
    
    if(anhystOnly_ == false){
      /*
       * Storage for return values
       */
      Vector<Double> retVec = Vector<Double>(dim_);
      retVec.Init();
      
      /*
       * Update and evaluate global rotation list
       */
      
      //    /*
      //     * check if copy works
      //     */
      //    std::cout << "GlobalRotationList pre updating " << std::endl;
      //
      //    std::list<RotListEntryv10>::iterator listIt;
      //    for(listIt = globRotList_[idElem].begin(); listIt != globRotList_[idElem].end(); listIt++){
      //    std::cout << listIt->ToString() << std::endl;
      //    }
      //    std::cout << "##########################" << std::endl;
      
      
      if(overwrite == true){
        /*
         * work on std data structure
         */
        if(debugOut){
          LOG_TRACE(vecpreisachInversion) << "Work on permanent storage" ;
        }
        Update_GlobalRotationList(X_thres, xVal, e_u, globRotList_[idElem],overwriteDirection,debugOut);
        //Update_GlobalRotationList(X_thres, xVal, e_u, globRotList_[idElem],true);
        
        /*
         * Evaluate_GlobalRotationList checks for each element if it was changed or not and
         * reevaluates only the ones that did
         */
        if(performanceMeasurement_){
          evaluateNestedListCounter_++;
          evaluateNestedListTimer_->Start();
        }
        
        Evaluate_GlobalRotationList(globRotList_[idElem], retVec);
        
        if(performanceMeasurement_){
          evaluateNestedListTimer_->Stop();
        }
        
      } else {
        /*
         * get copy of globRotList_[idElem]
         */
        if(debugOut){
          LOG_TRACE(vecpreisachInversion) << "Work on temporal storage" ;
        }
        
        if(performanceMeasurement_){
          copyToTemporalStorageCounter_++;
          copyToTemporalStorageTimer_->Start();
        }
        
        std::list<RotListEntryv10> tmpList = globRotList_[idElem];
        
        if(performanceMeasurement_){
          copyToTemporalStorageTimer_->Stop();
        }
        if(debugOut){
          LOG_TRACE(vecpreisachInversion) << "tmpList pre updating ";
          std::list<RotListEntryv10>::iterator listIt;
          for(listIt = tmpList.begin(); listIt != tmpList.end(); listIt++){
            LOG_TRACE(vecpreisachInversion) << listIt->ToString();
          }
          LOG_TRACE(vecpreisachInversion) << "##########################";
        }
        /*
         * then perform all updates on that temporal list only
         */
        //   std::cout << "Working only on temporal storage! " << std::endl;
        //Update_GlobalRotationList(X_thres, xVal, e_u, tmpList,true);
        Update_GlobalRotationList(X_thres, xVal, e_u, tmpList,overwriteDirection,debugOut);
        
        if(performanceMeasurement_){
          evaluateNestedListCounter_++;
          evaluateNestedListTimer_->Start();
        }
        
        Evaluate_GlobalRotationList(tmpList, retVec);
        
        if(performanceMeasurement_){
          evaluateNestedListTimer_->Stop();
        }
        
        //      std::cout << "retVec: " << retVec.ToString() << std::endl;
        if(debugOut){
          LOG_TRACE(vecpreisachInversion) << "tmpList post updating ";
          std::list<RotListEntryv10>::iterator listIt;
          for(listIt = tmpList.begin(); listIt != tmpList.end(); listIt++){
            LOG_TRACE(vecpreisachInversion) << listIt->ToString();
          }
          LOG_TRACE(vecpreisachInversion) << "##########################";
        }		
      }
      
      //    /*
      //     * check if copy works
      //     */
      //    if(idElem == 0){
      //    std::cout << "GlobalRotationList after updating and evaluation " << std::endl;
      //    //std::list<RotListEntryv10>::iterator listIt;
      //    for(listIt = globRotList_[idElem].begin(); listIt != globRotList_[idElem].end(); listIt++){
      //    std::cout << listIt->ToString() << std::endl;
      //    }
      //    std::cout << "##########################" << std::endl;
      //    }
      
      Yout += retVec;
      
      if(classical_){
        /*
         * Add value of lower triangle
         */
        Yout += e_u * lowerTriangleValue_;
      }
    }
    // add anhyst part
    Yout += anhystPart;
    
    //  std::cout << "###YOUT###################" << std::endl;
    //  std::cout << std::setprecision(16) << std::scientific << Yout.ToString() << std::endl;
    //  std::cout << "##########################" << std::endl;
    
    /*
     * in previous versions we scaled with delta_^2, too
     * this is no longer needed as we perform this step already in mapRectangleToPreisachWeights which is called
     * for the evaluation of all three parts
     */
    if(overwrite == true){
      preisachSum_[idElem] = Yout*YSaturated_;

      prevXVal_[idElem] = u_in;
      prevHVal_[idElem] = preisachSum_[idElem];
      
      return preisachSum_[idElem];
    } else {
      preisachSumTmp_[idElem] = Yout*YSaturated_;

      return preisachSumTmp_[idElem];
    }
  }
  
  Double VectorPreisachv10_ListApproach::clipRectangleToElement(Rectangle& source, UInt idAlpha, UInt idBeta, Double delta, bool isRotState){
    /*
     * Calculates the overlapping area of a rectangle (rectT,rectB,rectL,rectR) with the element defined by alphaId, betaId
     * This calculation does the following steps:
     *  1. calculate the overlapping rectangle via clipRectangles
     *  2. if success, check if lower left corner or upper right corner lie below diagonal
     *    -> if yes, cut down so that only the right corner is still below the diagonal
     *    -> subtract triangular area
     *
     * new flag: isRotState
     * -> if the rotation state is mapped to the helper matrix (in mapRectangleToHelperMatrix), we do not want to scale its value by 0.5 if its on the diagonal alpha = beta
     *    if we would do so, we would get the factor 0.5 twice in the final images which are created from the matrix (1. from switching state,
     *    2. from rotation state)
     * -> if isRotState == True, scale full entries on alpha = beta times 2 and for partially filled elements, skip
     *    the cutting of triangle parts below the diagonal
     *
     */
    
    if(delta <= 0){
      /*
       * default for clipping to Preisach elements
       */
      delta = delta_;
    }
    /*
     * else: use different delta for clipping to Bmp
     */
    
    Double elemLeft,elemBot;
    Double ovL,ovR,ovT,ovB;
    bool success;
    
    elemLeft = idBeta*delta - 1.0;
    elemBot = idAlpha*delta - 1.0;
    
    Rectangle elemRect = Rectangle(elemLeft,elemLeft+delta, elemBot+delta, elemBot);
    Rectangle target = Rectangle(0,0,0,0);
    
    success = elemRect.clipRectangles(source,target);
    
    if(!success){
      return 0.0;
    }
    
    target.getBounds(ovL,ovR,ovT,ovB);
    
    /*
     * check if element is a diagonal one, i.e. if it lies on alpha = beta -> idAlpha = idBeta
     */
    Double triang = 0.0;
    if(idAlpha == idBeta){
      /*
       * check if parts of the rectangle lie below the diagonal
       */
      if(ovL-elemLeft >= ovT-elemBot){
        /*
         * upper left corner below diagonal -> no overlap at all
         */
        return 0.0;
      }
      if(ovL-elemLeft >= ovB-elemBot){
        /*
         * bottom left corner below diagonal -> move bot up to diagonal -> ovB = ovL (remember alpha = beta here!)
         */
        ovB = ovL;
      }
      if(ovR-elemLeft >= ovT-elemBot){
        /*
         * upper right corner below diagonal -> move right to diagonal -> ovR = ovT
         */
        ovR = ovT;
      }
      
      if(ovR-elemLeft >= ovB-elemBot){
        /*
         * lower right corner below diagonal (only this triangular area is still below diagonal!)
         * -> subtract triangle
         */
        triang = 0.5*((ovR-ovB)*(ovR-ovB));
      }
    }
    
    if(isRotState == true){
      /*
       * do not subtract the part below the diagonal -> treat element as full element
       */
      triang = 0.0;
    }
    
    /*
     * calculate rectangular area
     * Note: execute this step after the cutting down steps above, so that the area below the diagonal is
     * at max a triangle (which value already is known)
     */
    Double area = (ovT-ovB)*(ovR-ovL);
    
    /*
     * subtract triangular part (if any)
     */
    area -= triang;
    
    return area;
  }
  
  void VectorPreisachv10_ListApproach::getBoundingBoxFromRotEntry(std::list<RotListEntryv10>::iterator rotListIt, Rectangle& rect, bool lastRotListEntryv10){
    /*
     * helper function returning a bounding box including all rotation areas which belong to a given rotListEntryv10
     *
     * Idea:
     *  when new values are added to the local switching list, these values can be clipped to this bounding box
     *
     */
    
    Double l,r,t,b;
    Double upperBound = rotListIt->getVal(); //exi
    Double lowerBound = rotListIt->getLowerVal(); //exi+1
    
    if(classical_ == true){
      /*
       * Rotation states form flipped L-shapes in S_U; last entry of list extends into T_U
       *
       *            S_U           alpha           T_U
       *      __ __ __ __ __ __ __ |_  __ __ __ __ __ __
       *     |     |   |      |    |                  /
       *     |A0   |A12|A22   |AN2 | AN2 (cont)     /
       * ex1_|_ _ _| _ |      |    |              /
       *     |A11      |      |    |            /
       * ex2_|_ _ _ _ _| _  _ |    |          /
       *     |A21             |    |        /
       *     |                |    |      /
       * exN_|_ _ _ _ _ _ _ _ | _ _| _ _/
       *     |AN1                  |  /
       *     |_____________________|/____________________ beta
       *
       *
       * Bounding boxes:
       *  if Entry is not last entry:
       *    l = -1; r = -lowerVal; t = +1; b = lowerVal
       *    (Bounding box includes all areas belong to this and to former rotListEntries
       *  if Entry is last entry:
       *    l = -1; r = +1; t = +1; b = 0
       *    (Bounding box includes all areas)
       *
       *
       */
      if(lastRotListEntryv10 == true){
        /*
         * last entry in list was found -> extend right up to +1
         */
        l = -1.0;
        r = 1.0;
        b = 0.0;
        t = 1.0;
        
      } else {
        /*
         * inner entry (also includes i=0)
         */
        l = -1.0;
        r = -lowerBound;
        b = lowerBound;
        t = 1.0;
      }
    } else {
      /*
       * Rotation states form flipped L-shapes in the whole Preisach plane
       *
       *            S_U           alpha           T_U
       *      __ __ __ __ __ __ __ |_  __ __ __ __ __ __
       *     |                     |                  /
       *     |A02                  |                /
       *     |_ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ /_ ex1
       *     |     |   A12         |            /
       *     |A01  |_ _  _ _ _ _ _ |_ _ _ _ _ /_ ex2
       *     |     |   |    A22    |        /
       *     |     |A11|           |      /
       *     |     |   |_ _ _ _ _ _| _ _/_ _ exN
       *     |     |   |      |    |  /
       *     |_____|___|_A21__|_AN_|/____________________ beta
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
       *
       *
       * Bounding boxes:
       *  l = -upperVal, r = upperVal, t = upperVal, b = -upperVal
       *
       */
      l = -upperBound;
      r = upperBound;
      b = -upperBound;
      t = upperBound;
    }
    
    rect.setBounds(l,r,t,b);
  }
  
  
  bool VectorPreisachv10_ListApproach::getRectanglesFromRotEntry(std::list<RotListEntryv10>::iterator rotListIt, Rectangle& rect1, Rectangle& rect2, bool lastRotListEntryv10){
    /*
     * -encapsulate the determination of rectangular rotation areas from a rot list entry
     * -return true if two non-zero rectangles are created; false otherwise
     * -lastRotListEntryv10 needed to check for last entry
     */
    
    bool twoAreas = true;
    Double l,r,t,b;
    Double l2,r2,t2,b2;
    Double upperBound = rotListIt->getVal(); //exi
    Double lowerBound = rotListIt->getLowerVal(); //exi+1
    
    if(classical_ == true){
      /*
       * Rotation states form flipped L-shapes in S_U; last entry of list extends into T_U
       *
       *            S_U           alpha           T_U
       *      __ __ __ __ __ __ __ |_  __ __ __ __ __ __
       *     |     |   |      |    |                  /
       *     |A0   |A12|A22   |AN2 | AN2 (cont)     /
       * ex1_|_ _ _| _ |      |    |              /
       *     |A11      |      |    |            /
       * ex2_|_ _ _ _ _| _  _ |    |          /
       *     |A21             |    |        /
       *     |                |    |      /
       * exN_|_ _ _ _ _ _ _ _ | _ _| _ _/
       *     |AN1                  |  /
       *     |_____________________|/____________________ beta
       *
       *     A0:
       *      L = -1, R = -ex1, B = ex1, T = 1
       *
       *     Ai1:
       *      L = -1, R = -exi+1, B = exi+1, T = exi
       *     Ai2:
       *      L = -exi, R = -exi+1, B = exi, T = 1
       *
       *     AN1:
       *      L = -1, R = exN, B = 0, T = exN
       *     AN2:
       *      L = -exN, R = 1, B = exN, T = 1
       *
       */
      if(upperBound >= 1){
        /*
         * A0 = 0; A11 = square; A12 = 0
         */
        twoAreas = false;
      }
      if(lastRotListEntryv10 == true){
        /*
         * AN -> last entry in list was found
         */
        /*
         * AN1
         */
        l = -1.0;
        r = upperBound;
        b = 0.0;
        t = upperBound;
        
        /*
         * AN2
         */
        l2 = -upperBound;
        r2 = 1;
        b2 = upperBound;
        t2 = 1;
      } else {
        /*
         * Ai (also includes i=0)
         */
        l = -1.0;
        r = -lowerBound;
        b = lowerBound;
        t = upperBound;
        
        l2 = -upperBound;
        r2 = -lowerBound;
        b2 = upperBound;
        t2 = 1.0;
      }
    } else {
      /*
       * Rotation states form flipped L-shapes in the whole Preisach plane
       *
       *            S_U           alpha           T_U
       *      __ __ __ __ __ __ __ |_  __ __ __ __ __ __
       *     |                     |                  /
       *     |A02                  |                /
       *     |_ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ /_ ex1
       *     |     |   A12         |            /
       *     |A01  |_ _  _ _ _ _ _ |_ _ _ _ _ /_ ex2
       *     |     |   |    A22    |        /
       *     |     |A11|           |      /
       *     |     |   |_ _ _ _ _ _| _ _/_ _ exN
       *     |     |   |      |    |  /
       *     |_____|___|_A21__|_AN_|/____________________ beta
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
       *
       *
       *     A01:
       *      L = -1, R = -ex1, B = -1, T = ex1
       *     A02:
       *      L = -1, R = 1, B = ex1, T = 1
       *
       *     Ai1:
       *      L = -exi, R = -exi+1, B = -exi, T = exi+1
       *     Ai2:
       *      L = -exi, R = exi, B = exi+1, T = exi
       *
       *     AN:
       *      L = -exN, R = exN, B = -exN, T = exN
       *
       */
      if(lastRotListEntryv10 == true){
        /*
         * AN -> last entry in list was found
         */
        twoAreas = false;
        l = -upperBound;
        r = upperBound;
        b = -upperBound;
        t = upperBound;
        
        l2 = 0;
        r2 = 0;
        b2 = 0;
        t2 = 0;
      } else {
        /*
         * Ai (also includes i=0)
         */
        l = -upperBound;
        r = -lowerBound;
        b = -upperBound;
        t = lowerBound;
        
        l2 = -upperBound;
        r2 = upperBound;
        b2 = lowerBound;
        t2 = upperBound;
      }
    }
    
    rect1.setBounds(l,r,t,b);
    rect2.setBounds(l2,r2,t2,b2);
    
    return twoAreas;
  }
  
  void VectorPreisachv10_ListApproach::Evaluate_GlobalRotationList(std::list<RotListEntryv10>& usedList, Vector<Double>& retVec){
    /*
     * Evaluates the weighted rotation state of
     *  a) upper square S_U and upper triangle T_U (classic)
     *  b) the whole Preisach plane (revised)
     *
     * Evaluation is done in the following way:
     *
     * for rotEntry in globalRotList:
     *   calculate two rectangular regions corresponding to the rotation area (forming a L)
     *
     *   for entry in localSwitchList:
     *     get rectangular bounds from entry in switch list
     *     clip rectangle with first rotation area
     *     clip result to Preisach plane
     *
     *     clip rectangle with second rotation area
     *     clip result to Preisach plane
     *
     *     set entry of switching list to not changed
     *     sum up both parts with appropriate signs
     *
     *   set entry of rotation list to not changed
     *   multiply sum with rotation state and add product to eval state
     *
     * save evaluated state in retVec
     *
     */
    
    /*
     * Init ret vec
     */
    if(retVec.GetSize() == 0){
      retVec = Vector<Double>(dim_);
    }
    retVec.Init(0);
    
    /*
     * iterators for local switching list
     */
    std::list<ListEntryv10>::iterator swListIt;
    std::list<ListEntryv10>::iterator swListStart;
    std::list<ListEntryv10>::iterator swListEnd;
    
    /*
     * iterators for global rotation list
     */
    std::list<RotListEntryv10>::iterator rotListIt;
    std::list<RotListEntryv10>::iterator rotListEnd = --(usedList.end());
    
    bool twoAreas = true;
    bool lastRotListEntryv10;
    
    Rectangle rotRect1 = Rectangle(0,0,0,0);
    Rectangle rotRect2 = Rectangle(0,0,0,0);
    Rectangle swRect = Rectangle(0,0,0,0);
    Rectangle overlapRect = Rectangle(0,0,0,0);
    
    Double sum = 0.0;
    Double tmp = 0.0;
    Double factor;
    UInt cnt = 0;
    UInt outercnt = 1;
    
    Vector<Double> rotState;
    
    for(rotListIt = usedList.begin(); rotListIt != usedList.end(); rotListIt++){
      
      if(rotListIt == rotListEnd){
        lastRotListEntryv10 = true;
      } else {
        lastRotListEntryv10 = false;
      }
      rotState = rotListIt->getVecReference();
      
      /*
       * check if reevaluation is needed at all
       */
      if(rotListIt->hasChanged() == false){
        //        std::cout << "Rotstate did not change > reuse old value" << std::endl;
        /*
         * neither switching list did change since last time (and weights did not change either)
         * nor did the lower bound of the rotation area
         * -> rotation area and all switching areas are the same as last time
         * -> overlap of switching areas and rotation area does not have to be computed again
         * -> take stored value
         */
        retVec += rotState*rotListIt->getLastEvalState();
      } else {
        /*
         * get rotation area(s)
         */
        
        swListStart = rotListIt->getListReference().begin();
        swListEnd = --(rotListIt->getListReference().end());
        
        twoAreas = getRectanglesFromRotEntry(rotListIt, rotRect1, rotRect2,lastRotListEntryv10);
        
        //        std::cout << "---------------------------" << std::endl;
        //        std::cout << "Rotlistentry " << outercnt << std::endl;
        //        Double ltmp,rtmp,ttmp,btmp;
        //        rotRect1.getBounds(ltmp,rtmp,ttmp,btmp);
        //        std::cout << "Rect1: " << ltmp << ", " << rtmp << ", " << ttmp << ", " << btmp << std::endl;
        //        rotRect2.getBounds(ltmp,rtmp,ttmp,btmp);
        //        std::cout << "Rect2: " << ltmp << ", " << rtmp << ", " << ttmp << ", " << btmp << std::endl;
        outercnt++;
        
        /*
         * reset counter (needed to check if we have the first entry of the switch list)
         * Regarding cnt:
         *  the absolute value of cnt is not relevant. We only have to check if it is 0 or not.
         *  In case that it is 0 we have to use the first switching entry twice (once for area with
         *  id = 0 and once for area with id = 1). Normally (and in older versions), we start at 0.
         *  However, the lists started to get very long with lots of entries which had no influence on
         *  the result (i.e. the switching areas did not overlap with the rotation area). To reduce the
         *  computational cost and to shorten the lists, the function simplify_switchinglist was introduced.
         *  This function checks directly after the merging step in update_globalrotationlist, all entries
         *  for a possible overlap. Entries at the beginning of the list, which do not lead to an overlap
         *  can be removed from the list (elements at the end which have no influence will not get inserted
         *  in the first place). The first entry of this shortened list does not correspond to area id 0
         *  anymore (as this is the special helper area). So we startCnt in rotListEntryv10 to mark that we
         *  deleted the first entry which originally corresponded to area 0. As all areas with id > 0 are
         *  calculated by the same scheme, the absolute value of cnt is not of interest and thus it is
         *  either 0 or 1.
         */
        cnt = rotListIt->getStartCnt();
        sum = 0.0;
        tmp = 0.0;
        bool success = false;
        
        //        std::cout << "RotRect1" << std::endl;
        //        std::cout << rotRect1.ToString() << std::endl;
        //
        //        std::cout << "RotRect2" << std::endl;
        //        std::cout << rotRect2.ToString() << std::endl;
        
        for(swListIt = rotListIt->getListReference().begin(); swListIt != rotListIt->getListReference().end(); ){
          
          /*
           * get rectangular bounds corresponding to entry of switching list (compare to Evaluate_GlobalSwitchingList)
           */
          factor = getRectangleFromSwitchingList(rotListIt->getListReference(),swListStart, swListIt, swListEnd,cnt, swRect);
          
          /*
           * clip resulting area against rotation area1
           */
          success = rotRect1.clipRectangles(swRect,overlapRect);
          
          //          std::cout << "Switching" << std::endl;
          //          std::cout << swRect.ToString() << std::endl;
          
          if(success){
            /*
             *  clip overlap to PreisachPlane and sum up
             */
            tmp = mapRectangleToPreisachWeights(overlapRect);
            
            //            std::cout << "Overlap1" << std::endl;
            //            std::cout << overlapRect.ToString() << std::endl;
            //            std::cout << "Factor: " << factor << std::endl;
            //            std::cout << "Value: " << tmp << std::endl;
            
            sum += factor*tmp;
          }
          
          if(twoAreas){
            /*
             * repeat the clipping steps for the second area
             */
            success = rotRect2.clipRectangles(swRect,overlapRect);
            
            if(success){
              /*
               *  clip overlap to PreisachPlane and sum up
               */
              tmp = mapRectangleToPreisachWeights(overlapRect);
              
              //             std::cout << "Overlap2" << std::endl;
              //             std::cout << overlapRect.ToString() << std::endl;
              //             std::cout << "Factor: " << factor << std::endl;
              //             std::cout << "Value: " << tmp << std::endl;
              
              sum += factor*tmp;
            }
          }
          
          /*
           * extra treatment for unsymmetric weights
           * here we have to overlap with the diagonally split area, too.
           * (area 0, area id = -1)
           *
           *      split area = area 0
           *       _v_______________
           *      | \  | 1 |    |   |
           *      |___\|___| 3  |   |
           *      |   2    |    |5  |
           *      |________|____|   |
           *      |     4       |   |
           *      |_____________|___|
           *      |        6        |
           *      |_________________|
           *
           * for symmetric weights, this area will always cancel out as positive and negative
           * part are equal. For unsymmetric weights, we have to evaluate, but only if this special
           * area has not been wiped out
           * Has only be checked once for each rotation state -> cnt == 0
           * (Note that cnt can start at values > 0, if the switching list was simplified.
           * In such a case, we can be sure that the split area does not overlap with a rotation state
           * as the areas 1 and 2 do not overlap either.)
           *
           */
          if(cnt == 0){
            
            if((isSymmetric_ == false)&&(rotListIt->wasListWipedOut() == false)){
              /*
               * set flag upperSplitSquare to true -> get area 0 instead of area 1
               */
              factor = getRectangleFromSwitchingList(rotListIt->getListReference(),swListStart, swListIt, swListEnd,cnt, swRect,true);
              
              if(factor != 0){
                // getRectangleFromSwitchingList returns the correct sign for
                // the addition (i.e. +1 or -1 depending on minima/maxima)
                // exception: area0, as this area is half +1 and half -1
                // for this special case, getRectangleFromSwitchingList returns 0
                // and the signs are considered in mapRectanglesToPreisachPlane
                EXCEPTION("Something got wrong here! Check function getRectangleBounds!");
              }
              
              /*
               * clip resulting area against rotation area1
               */
              success = rotRect1.clipRectangles(swRect,overlapRect);
              
              if(success){
                /*
                 *  clip overlap to PreisachPlane and sum up
                 *  -> NOTE: here we set the flag skipUpperDiagonal to true!
                 */
                tmp = mapRectangleToPreisachWeights(overlapRect,true);
                // factor is not needed here; mapRectangleToPreisachWeights with skipUpperDiagonal == true
                // will not only skip the diagonal entries, but also consider the right signs!
                sum += tmp;
              }
              
              if(twoAreas){
                /*
                 * repeat the clipping steps for the second area
                 */
                success = rotRect2.clipRectangles(swRect,overlapRect);
                
                if(success){
                  /*
                   *  clip overlap to PreisachPlane and sum up
                   *  -> NOTE: here we set the flag skipUpperDiagonal to true!
                   */
                  tmp = mapRectangleToPreisachWeights(overlapRect,true);
                  // factor is not needed here; mapRectangleToPreisachWeights with skipUpperDiagonal == true
                  // will not only skip the diagonal entries, but also consider the right signs!
                  sum += tmp;
                }
              }
            }
          }
          
          if(cnt > 0){
            /*
             * NOTE: area1 and area2 are both calculated using the first list entry, therefore we do not increase the iterator after
             * the first iteration
             */
            swListIt++;
          }
          cnt++;
        } // sw list
        
        /*
         * set rotation element to be unchanged
         * -> if neither the boundaries of the rotation area, nor the switching list change till next time
         * -> reuse value
         */
        rotListIt->setToUnchanged();
        rotListIt->setLastEvalState(sum);
        
        /*
         * add sum * rotState to retVec
         */
        //        std::cout << "Sum: " << sum << std::endl;
        //        std::cout << "rotState: " << rotState.ToString() << std::endl;
        retVec += rotState*sum;
        //        std::cout << "retVec: " << retVec.ToString() << std::endl;
      }
    } // rot list
  }
  
  void VectorPreisachv10_ListApproach::Evaluate_LowerTriangle(){
    /*
     * This function calculates the overall switching state of the lower triangle;
     * This value has to be multiplied with current rotation state e_u in each iteration to get
     * the overall contribution of the lowe triangle
     * -> this function has to be called only once during the initializatioh as the state is always +1
     */
    
    /*
     * use function mapRectangleToPreisachWeights for the evaluation
     * this function overlaps rectangles with the Preisach plane and intergrates over all weights in the
     * overlapping (parts below the diagonal alpha = beta will be cut away accordingly)
     * -> as overlapping area use the lower square -1 <= alpha <= 0, -1 <= beta <= 0
     */
    Rectangle rect = Rectangle(-1.0,0.0,0.0,-1.0);
    lowerTriangleValue_ = mapRectangleToPreisachWeights(rect);
  }
  
  Double VectorPreisachv10_ListApproach::mapRectangleToPreisachWeights(Rectangle& rect, bool skipUpperDiagonal){
    
    if(mappingVersion_ == 0){
      //if(textOutputLevel_ == 2){
      // std::cout << "Use OLD mapping method" << std::endl;
      //}
      return mapRectangleToPreisachWeightsOLD(rect, skipUpperDiagonal);
    } else {
      //if(textOutputLevel_ == 2){
      // std::cout << "Use NEW mapping method" << std::endl;
      //}
      return mapRectangleToPreisachWeightsNEW(rect, skipUpperDiagonal);
    }
  }
  Double VectorPreisachv10_ListApproach::mapRectangleToPreisachWeightsOLD(Rectangle& rect, bool skipUpperDiagonal){
    /*
     * Input: rectangle area described by its top (t), bottom (b), left (l) and right (r) boundary
     * Output: Sum over all PreisachWeights overlapped by the rectangle (partially overlapped elements
     * are added only partially!)
     *
     * skipUpperDiagonal:
     *  used for the calculation of the upper square area which is split along the diagonal alpha=-beta
     *  Normally, we can leave out this area, as positive and negative part cancel out. This is not true
     *  if the Preisach weights are no longer symmetric to alpha = -beta. In that case we have to sum up
     *  parts of both halves of this upper square (only the parts overlapping with a rotation state).
     *  To make calculation easier, we skip the diagonal entries on alpha = -beta. This is valid, as one
     *  element of the Preisach plane is always assumed to be symmetric, so that a single element which
     *  is split along its diagonal sums up to 0.
     *  Furthermore, skipUpperDiagonal will consider the correct signs, i.e. all entries above alpha = -beta will be subtracted
     *  all entries below will be added to the return value!
     */
    
    Double sum = 0.0;
    Double t,b,l,r;
    rect.getBounds(l,r,t,b);
    
    //   std::cout << "L,R,T,B: " << l << "," << r << "," << t << "," << b << "," << std::endl;
    
    /*
     * restrict input value to valid range
     * (Preisach plane just goes from -1 to +1 in alpha and beta)
     */
    t = std::min(t,1.0);
    b = std::max(b,-1.0);
    l = std::max(l,-1.0);
    r = std::min(r,1.0);
    
    //   std::cout << "L,R,T,B: " << l << "," << r << "," << t << "," << b << "," << std::endl;
    
    /*
     * check if rectangular has size != 0
     */
    if((t <= b)||(r <= l)){
      return sum;
    }
    
    /*
     * Lower triangluar part of Preisach plane (here 16 elements)
     * and overlapping rectangle
     *   _____ _____ _____ _____ _____
     *  |    .|.....|.....|.....|.   /
     *  | 1 xx | 2   | 3   | 4   |xx5/
     *  |___xx_|_____|_____|_____|/
     *  |   xx |     |     |    / xx
     *  | 6 xx | 7   | 8   | 9/   xx
     *  |___xx_|_____|_____|/     xx
     *  |   xx |     |    /       xx
     *  | a xx.|.b...|.c/.........xx
     *  |_____|_____|/
     *  |     |    /
     *  | d   | e/
     *  |_____|/
     *  |    /
     *  | f/
     *  |/
     *
     *  Element 7,8,9 are completely overlapped
     *  Element 6 is cut vertically
     *  Element 2,3,b,c are cut horizontally
     *  Element a is cut horizontally and vertically
     *  Element c,5 are cut horizontally, vertically and diagonally
     *
     *  Approach:
     *   Iterate over all (fully and partially overlapped) elements
     *    >Check if indices denote a fully overlapped element
     *      > if yes sum up value
     *      > else use clipToElement to determine the appropriate scaling value for weight, then sum up
     */
    
    /*
     * 1. get indices of all overlapped elements (partially and fully)
     * NOTE: element alpha = -1 and beta = -1 will have index 0,0!
     * -> add floor(1.0/delta_) != numRows/2
     */
    int rowMin =  std::floor(b/delta_) + std::floor(1.0/delta_);
    int rowMax = std::ceil(t/delta_) + std::floor(1.0/delta_)-1;
    
    int colMin =  std::floor(l/delta_) + std::floor(1.0/delta_);
    int colMax =  std::ceil(r/delta_) + std::floor(1.0/delta_)-1;
    
    /*
     * NEW: use integers instead of unsigned integer
     * REASON: colMaxFull became negative (which simply would indicate no fully overlapped elemets), but
     * unsigning took it to a very large positive value and thus -> wrong result
     */
    
    /*
     * 2. get indices of completely overlapped elements
     */
    int rowMinFull =  std::ceil(b/delta_) + std::floor(1.0/delta_);
    int rowMaxFull = std::floor(t/delta_) + std::floor(1.0/delta_)-1;
    
    int colMinFull =  std::ceil(l/delta_) + std::floor(1.0/delta_);
    int colMaxFull = std::floor(r/delta_) + std::floor(1.0/delta_)-1;
    
    //    std::cout << std::setprecision(12) << std::scientific << b << " " << t << " " << l << " " << r << std::endl;
    //    std::cout << std::setprecision(12) << std::scientific << floor(b) << " " << ceil(t) << " " << floor(l) << " " << ceil(r) << std::endl;
    //    std::cout << std::setprecision(12) << std::scientific << floor(b/delta_) << " " << ceil(t/delta_) << " " << floor(l) << " " << ceil(r) << std::endl;
    //    std::cout << rowMin << " " << rowMax << " " << colMin << " " << colMax << std::endl;
    //    std::cout << rowMinFull << " " << rowMaxFull << " " << colMinFull << " " << colMaxFull << std::endl;
    
    
    //    Matrix<Double> accessCounter = preisachWeights_;
    ////
    //    accessCounter.Init();
    
    //   std::cout << "AccessCounter (start): " << accessCounter.ToString() << std::endl;
    
    /*
     * Iterate over all elements
     */
    Double tmp = 0.0;
    
    if(skipUpperDiagonal == true){
      /*
       * special treatment:
       * 1. skip entries on alpha = -beta
       * 2. consider the correct signs, i.e. +1 for entries below alpha = -beta
       *    and -1 for entries above!
       */
      Double sign;
      for(int ii = rowMin; ii <= rowMax; ii++){
        /*
         * ensure, that j <= i -> other elements do not contribute
         */
        UInt i = (UInt) ii;
        
        for(int jj = colMin; jj <= std::min(ii,colMax); jj++){
          
          UInt j = (UInt) jj;
          
          if(numRows_-i-j-1 == 0){
            /*
             * alpha = -beta
             */
            continue;
          } else if (numRows_ < i+j+1){
            /*
             * above alpha = -beta
             * -> sign = -1
             */
            sign = -1.0;
          } else {
            /*
             * below diagonal
             * -> sign = +1
             */
            sign = 1.0;
          }
          
          tmp = preisachWeights_[i][j];
          
          /*
           * Check for an inner element
           */
          if((ii >= rowMinFull)&&(ii <= rowMaxFull)){
            if((jj >= colMinFull)&&(jj <= colMaxFull)){
              /*
               * check for diagonal element on alpha = beta
               * (should not happen here)
               */
              if(i == j){
                tmp = tmp/2.0;
              }
              
              /*
               * scale Preisach weights with area
               */
              tmp *= delta_*delta_;
            } else {
              /*
               * Partially overlapped element
               * -> clip rectangle to element
               * Note: it does not matter here, if we pass the full input area or the
               * subarea defined by the current pair i,j as the clip function will overlap
               * the element area with the provided input area and automatically find the
               * correct overlap; the clipToElement function will furthermore check for
               * areas below the diagonal alpha = beta and cut them from the overlapping area
               *
               * Note 2: here we do not have to scale with delta_^2 as the overlap is already a returning an
               * area <= delta_^2
               */
              tmp = tmp * clipRectangleToElement(rect, i, j);
            }
          } else {
            /*
             * Partially overlapped element
             * -> clip rectangle to element
             * Note: it does not matter here, if we pass the full input area or the
             * subarea defined by the current pair i,j as the clip function will overlap
             * the element area with the provided input area and automatically find the
             * correct overlap; the clipToElement function will furthermore check for
             * areas below the diagonal alpha = beta and cut them from the overlapping area
             *
             * Note 2: here we do not have to scale with delta_^2 as the overlap is already a returning an
             * area <= delta_^2
             */
            tmp = tmp * clipRectangleToElement(rect, i, j);
          }
          
          sum += sign*tmp;
        }
      }
      
    } else {
      /*
       * std treatment
       */
      /*
       * sum over all fully overlapped Elements
       */
      Double sumFullElems = 0.0;
      for(int ii = rowMin; ii <= rowMax; ii++){
        /*
         * ensure, that j <= i -> other elements do not contribute
         */
        UInt i = (UInt) ii;
        
        for(int jj = colMin; jj <= std::min(ii,colMax); jj++){
          UInt j = (UInt) jj;
          
          tmp = preisachWeights_[i][j];
          
          /*
           * Check for an inner element
           */
          if((ii >= rowMinFull)&&(ii <= rowMaxFull)){
            if((jj >= colMinFull)&&(jj <= colMaxFull)){
              /*
               * check for diagonal element on alpha = beta
               */
              if(i == j){
                tmp = tmp/2.0;
                // accessCounter[i][j] += 0.5;
              }
              //              else {
              //                accessCounter[i][j] += 1;
              //              }
              
              sumFullElems += tmp;
              
              
              /*
               * new: scale only at end of summation to save flops
               */
              /*
               * scale Preisach weights with area
               */
              //tmp *= delta_*delta_;
            } else {
              /*
               * Partially overlapped element
               * -> clip rectangle to element
               * Note: it does not matter here, if we pass the full input area or the
               * subarea defined by the current pair i,j as the clip function will overlap
               * the element area with the provided input area and automatically find the
               * correct overlap; the clipToElement function will furthermore check for
               * areas below the diagonal alpha = beta and cut them from the overlapping area
               *
               * Note 2: here we do not have to scale with delta_^2 as the overlap is already a returning an
               * area <= delta_^2
               */
              sum += tmp * clipRectangleToElement(rect, i, j);
              // accessCounter[i][j] += clipRectangleToElement(rect, i, j);
              //tmp = tmp * clipRectangleToElement(rect, i, j);
            }
          } else {
            /*
             * Partially overlapped element
             * -> clip rectangle to element
             * Note: it does not matter here, if we pass the full input area or the
             * subarea defined by the current pair i,j as the clip function will overlap
             * the element area with the provided input area and automatically find the
             * correct overlap; the clipToElement function will furthermore check for
             * areas below the diagonal alpha = beta and cut them from the overlapping area
             *
             * Note 2: here we do not have to scale with delta_^2 as the overlap is already a returning an
             * area <= delta_^2
             */
            sum += tmp * clipRectangleToElement(rect, i, j);
            // accessCounter[i][j] += clipRectangleToElement(rect, i, j);
            //tmp = tmp * clipRectangleToElement(rect, i, j);
          }
          
          /*
           * new: all partially filled elements get added directly to sum;
           * all others are added at end of loop
           */
          //sum += tmp;
        }
      }
      sum += sumFullElems*delta_*delta_;
    }
    
    //    std::cout << "rowMin / rowMinFull: " << rowMin << " / " << rowMinFull << std::endl;
    //    std::cout << "rowMax / rowMaxFull: " << rowMax << " / " << rowMaxFull << std::endl;
    //    std::cout << "colMin / colMinFull: " << colMin << " / " << colMinFull << std::endl;
    //    std::cout << "colMax / colMaxFull: " << colMax << " / " << colMaxFull << std::endl;
    //    std::cout << "AccessCounter (end): " << accessCounter.ToString() << std::endl;
    
    return sum;
  }
  
  //NEW
  Double VectorPreisachv10_ListApproach::mapRectangleToPreisachWeightsNEW(Rectangle& rect, bool skipUpperDiagonal){
    /*
     * Input: rectangle area described by its top (t), bottom (b), left (l) and right (r) boundary
     * Output: Sum over all PreisachWeights overlapped by the rectangle (partially overlapped elements
     * are added only partially!)
     *
     * skipUpperDiagonal:
     *  used for the calculation of the upper square area which is split along the diagonal alpha=-beta
     *  Normally, we can leave out this area, as positive and negative part cancel out. This is not true
     *  if the Preisach weights are no longer symmetric to alpha = -beta. In that case we have to sum up
     *  parts of both halves of this upper square (only the parts overlapping with a rotation state).
     *  To make calculation easier, we skip the diagonal entries on alpha = -beta. This is valid, as one
     *  element of the Preisach plane is always assumed to be symmetric, so that a single element which
     *  is split along its diagonal sums up to 0.
     *  Furthermore, skipUpperDiagonal will consider the correct signs, i.e. all entries above alpha = -beta will be subtracted
     *  all entries below will be added to the return value!
     */
    
    Double sum = 0.0;
    Double t,b,l,r;
    rect.getBounds(l,r,t,b);
    
    //   std::cout << "L,R,T,B: " << l << "," << r << "," << t << "," << b << "," << std::endl;
    
    /*
     * restrict input value to valid range
     * (Preisach plane just goes from -1 to +1 in alpha and beta)
     */
    t = std::min(t,1.0);
    b = std::max(b,-1.0);
    l = std::max(l,-1.0);
    r = std::min(r,1.0);
    
    //   std::cout << "L,R,T,B: " << l << "," << r << "," << t << "," << b << "," << std::endl;
    
    /*
     * check if rectangular has size != 0
     */
    if((t <= b)||(r <= l)){
      return sum;
    }
    
    /*
     * Lower triangluar part of Preisach plane (here 16 elements)
     * and overlapping rectangle
     *   _____ _____ _____ _____ _____
     *  |    .|.....|.....|.....|.   /
     *  | 1 ¦ | 2   | 3   | 4   |¦5/
     *  |___¦_|_____|_____|_____|/
     *  |   ¦ |     |     |    / ¦
     *  | 6 ¦ | 7   | 8   | 9/   ¦
     *  |___¦_|_____|_____|/     ¦
     *  |   ¦ |     |    /       ¦
     *  | a ¦.|.b...|.c/.........¦
     *  |_____|_____|/
     *  |     |    /
     *  | d   | e/
     *  |_____|/
     *  |    /
     *  | f/
     *  |/
     *
     *  Element 7,8,9 are completely overlapped
     *  Element 6 is cut vertically
     *  Element 2,3,b,c are cut horizontally
     *  Element a is cut horizontally and vertically
     *  Element c,5 are cut horizontally, vertically and diagonally
     *
     *  Approach:
     *   Iterate over all (fully and partially overlapped) elements
     *    >Check if indices denote a fully overlapped element
     *      > if yes sum up value
     *      > else use clipToElement to determine the appropriate scaling value for weight, then sum up
     */
    
    /*
     * 1. get indices of all overlapped elements (partially and fully)
     * NOTE: element alpha = -1 and beta = -1 will have index 0,0!
     * -> add floor(1.0/delta_) != numRows/2
     */
    int rowMin =  std::floor(b/delta_) + std::floor(1.0/delta_);
    int rowMax = std::ceil(t/delta_) + std::floor(1.0/delta_)-1;
    
    int colMin =  std::floor(l/delta_) + std::floor(1.0/delta_);
    int colMax =  std::ceil(r/delta_) + std::floor(1.0/delta_)-1;
    
    /*
     * NEW: use integers instead of unsigned integer
     * REASON: colMaxFull became negative (which simply would indicate no fully overlapped elemets), but
     * unsigning took it to a very large positive value and thus -> wrong result
     */
    
    /*
     * 2. get indices of completely overlapped elements
     */
    int rowMinFull =  std::ceil(b/delta_) + std::floor(1.0/delta_);
    int rowMaxFull = std::floor(t/delta_) + std::floor(1.0/delta_)-1;
    
    int colMinFull =  std::ceil(l/delta_) + std::floor(1.0/delta_);
    int colMaxFull = std::floor(r/delta_) + std::floor(1.0/delta_)-1;
    
    //    std::cout << std::setprecision(12) << std::scientific << b << " " << t << " " << l << " " << r << std::endl;
    //    std::cout << std::setprecision(12) << std::scientific << floor(b) << " " << ceil(t) << " " << floor(l) << " " << ceil(r) << std::endl;
    //    std::cout << std::setprecision(12) << std::scientific << floor(b/delta_) << " " << ceil(t/delta_) << " " << floor(l) << " " << ceil(r) << std::endl;
    //    std::cout << rowMin << " " << rowMax << " " << colMin << " " << colMax << std::endl;
    //    std::cout << rowMinFull << " " << rowMaxFull << " " << colMinFull << " " << colMaxFull << std::endl;
    
    /*
     * new (hopefully cheaper) approach
     *    instead of iterating over all elements and checking for each element
     *    if it is
     *      a) on the diagonal alpha = beta
     *      b) partially overlapped
     *    we do
     *      1) loop over all fully overlapped elements excluding diagonal
     *      2) loop over all fully overlapped elements on the diagonal
     *      3) iterate over all partially overlapped elements
     *        3a) bottom
     *        3b) left
     *        3c) top
     *        3d) right
     *
     *    for upper split upper part, we do
     *      1) loop over all fully overlapped elements excluding diagonal
     *      2) only fully overlapped diagonal elements
     *      3) iterate over all partially overlapped elements
     
     *    -> only relevant for case of unsymmetric weights!
     *    TODO: unsymmetric case has to be tested
     *
     *
     */
    
    Double tmp = 0.0;
    Double fullElementSum = 0.0;
    
    //    Matrix<Double> accessCounter = preisachWeights_;
    //
    //    accessCounter.Init();
    
    //std::cout << "AccessCounter (start): " << accessCounter.ToString() << std::endl;
    
    
    if(skipUpperDiagonal == true){
      /*
       * special treatment:
       * 1. skip entries on alpha = -beta
       * 2. consider the correct signs, i.e. +1 for entries below alpha = -beta
       *    and -1 for entries above!
       */
      
      /*
       * new approach
       * 1) only full elements without diagonal
       */
      for(int ii = rowMinFull; ii <= rowMaxFull; ii++){
        UInt i = (UInt) ii;
        
        // ii-1 -> diagonal not reached
        for(int jj = colMinFull; jj <= std::min(ii-1,colMaxFull); jj++){
          UInt j = (UInt) jj;
          
          if(numRows_-i-j-1 == 0){
            /*
             * alpha = -beta
             */
            continue;
          }
          
          if (numRows_ < i+j+1){
            /*
             * above alpha = -beta
             * -> sign = -1
             */
            fullElementSum -= preisachWeights_[i][j];
          } else {
            /*
             * below diagonal
             * -> sign = +1
             */
            fullElementSum += preisachWeights_[i][j];
          }
        }
      }
      
      /*
       * new approach
       * 2) only full diagonal elements
       */
      for(int ii = rowMinFull; ii <= rowMaxFull; ii++){
        UInt i = (UInt) ii;
        
        // ii-1 was already treated above; now check only if element on diagonal
        // colMaxFull == ii is needed, too
        if(colMaxFull >= ii){
          UInt j = i;
          
          if(numRows_-i-j-1 == 0){
            /*
             * alpha = -beta
             */
            continue;
          }
          
          if (numRows_ < i+j+1){
            /*
             * above alpha = -beta
             * -> sign = -1
             */
            fullElementSum -= preisachWeights_[i][j]/2.0;
          } else {
            /*
             * below diagonal
             * -> sign = +1
             */
            fullElementSum += preisachWeights_[i][j]/2.0;
          }
        }
      }
      
      // full elements have an area of delta_^2
      // including the ones on the diagonal as these have been scaled
      // by 2.0 already
      sum += delta_*delta_*fullElementSum;
      
      /*
       * partial bottom elements (including corners)
       */
      if(rowMinFull > rowMin){
        int ii = rowMin;
        UInt i = (UInt) rowMin;
        
        for(int jj = colMin; jj <= std::min(ii,colMax); jj++){
          UInt j = (UInt) jj;
          
          if(numRows_-i-j-1 == 0){
            /*
             * alpha = -beta
             */
            continue;
          }
          
          // clip to rectangles will return the
          // the partially overlapped area;
          // this value is already scaled by delta_^2
          // additonally, on alpha = beta, the part below the diagonal is
          // subtracted
          tmp = clipRectangleToElement(rect, i, j);
          
          if (numRows_ < i+j+1){
            /*
             * above alpha = -beta
             * -> sign = -1
             */
            sum -= tmp*preisachWeights_[i][j];
          } else {
            /*
             * below diagonal
             * -> sign = +1
             */
            sum += tmp*preisachWeights_[i][j];
          }
        }
      }
      
      /*
       * partial top elements (including corners)
       */
      if((rowMaxFull < rowMax)&&(rowMin != rowMax)){
        int ii = rowMax;
        UInt i = (UInt) rowMax;
        
        for(int jj = colMin; jj <= std::min(ii,colMax); jj++){
          UInt j = (UInt) jj;
          
          if(numRows_-i-j-1 == 0){
            /*
             * alpha = -beta
             */
            continue;
          }
          
          // clip to rectangles will return the
          // the partially overlapped area;
          // this value is already scaled by delta_^2
          tmp = clipRectangleToElement(rect, i, j);
          
          if (numRows_ < i+j+1){
            /*
             * above alpha = -beta
             * -> sign = -1
             */
            sum -= tmp*preisachWeights_[i][j];
          } else {
            /*
             * below diagonal
             * -> sign = +1
             */
            sum += tmp*preisachWeights_[i][j];
          }
        }
      }
      
      /*
       * partial right elements (excluding corners)
       * colIndex jj has to be
       * <= rowIndex
       *
       * start with right edge
       * > reason: if the first column is partially filled, we
       * have colMaxFull = -1; colMax = 1 > working
       * but colMinFull = colMin = 0 > not working
       */
      if(colMaxFull < colMax){
        int jj = colMax;
        UInt j = (UInt) colMax;
        
        for(int ii = std::max(rowMinFull,jj); ii <= rowMaxFull; ii++){
          UInt i = (UInt) ii;
          
          if(numRows_-i-j-1 == 0){
            /*
             * alpha = -beta
             */
            continue;
          }
          
          // clip to rectangles will return the
          // the partially overlapped area;
          // this value is already scaled by delta_^2
          tmp = clipRectangleToElement(rect, i, j);
          
          if (numRows_ < i+j+1){
            /*
             * above alpha = -beta
             * -> sign = -1
             */
            sum -= tmp*preisachWeights_[i][j];
          } else {
            /*
             * below diagonal
             * -> sign = +1
             */
            sum += tmp*preisachWeights_[i][j];
          }
        }
      }
      
      /*
       * partial left elements (excluding corners)
       * colIndex jj has to be
       * <= rowIndex
       * -> ii starts at jj, but only if it is larger than rowMinFull
       * (otherwise the corner could be included twice)
       * -> ii starts at jj, but only if it is larger than rowMinFull
       * (otherwise the corner could be included twice)
       */
      if((colMinFull > colMin)&&(colMin != colMax)){
        int jj = colMin;
        UInt j = (UInt) colMin;
        
        for(int ii = std::max(rowMinFull,jj); ii <= rowMaxFull; ii++){
          UInt i = (UInt) ii;
          
          if(numRows_-i-j-1 == 0){
            /*
             * alpha = -beta
             */
            continue;
          }
          
          // clip to rectangles will return the
          // the partially overlapped area;
          // this value is already scaled by delta_^2
          tmp = clipRectangleToElement(rect, i, j);
          
          if (numRows_ < i+j+1){
            /*
             * above alpha = -beta
             * -> sign = -1
             */
            sum -= tmp*preisachWeights_[i][j];
          } else {
            /*
             * below diagonal
             * -> sign = +1
             */
            sum += tmp*preisachWeights_[i][j];
          }
        }
      }
      
      
    } else {
      /*
       * standard treatment (i.e. no skipping of elements, no signing, just adding up)
       * Remark: the correct sign for the whole sum if obtained from function
       *         getRectangleFromSwitchingList during the evaluation in Evaluate_GlobalRotation
       */
      
      /*
       * new approach
       * 1) only full elements without diagonal
       */
      for(int ii = rowMinFull; ii <= rowMaxFull; ii++){
        UInt i = (UInt) ii;
        
        // ii-1 -> diagonal not reached
        for(int jj = colMinFull; jj <= std::min(ii-1,colMaxFull); jj++){
          UInt j = (UInt) jj;
          
          fullElementSum += preisachWeights_[i][j];
          // accessCounter[i][j] += 1;
        }
      }
      
      /*
       * new approach
       * 2) only full diagonal elements
       */
      for(int ii = rowMinFull; ii <= rowMaxFull; ii++){
        UInt i = (UInt) ii;
        
        // ii-1 was already treated above; now check only if element on diagonal
        // colMaxFull == ii is needed, too
        if(colMaxFull >= ii){
          UInt j = i;
          
          fullElementSum += preisachWeights_[i][j]/2.0;
          //  accessCounter[i][j] += 0.5;
        }
      }
      
      // full elements have an area of delta_^2
      // including the ones on the diagonal as these have been scaled
      // by 2.0 already
      sum += delta_*delta_*fullElementSum;
      
      /*
       * partial bottom elements (including corners)
       */
      if(rowMinFull > rowMin){
        int ii = rowMin;
        UInt i = (UInt) rowMin;
        
        for(int jj = colMin; jj <= std::min(ii,colMax); jj++){
          UInt j = (UInt) jj;
          
          // clip to rectangles will return the
          // the partially overlapped area;
          // this value is already scaled by delta_^2
          // additonally, on alpha = beta, the part below the diagonal is
          // subtracted
          sum += clipRectangleToElement(rect, i, j)*preisachWeights_[i][j];
          //  accessCounter[i][j] += clipRectangleToElement(rect, i, j);
        }
      }
      
      /*
       * partial top elements (including corners)
       * check if rowMin and rowMax are the same > in this case we only have
       * one partial row which was already treated above
       */
      if((rowMaxFull < rowMax)&&(rowMin != rowMax)){
        int ii = rowMax;
        UInt i = (UInt) rowMax;
        
        for(int jj = colMin; jj <= std::min(ii,colMax); jj++){
          UInt j = (UInt) jj;
          
          sum += clipRectangleToElement(rect, i, j)*preisachWeights_[i][j];
          //  accessCounter[i][j] += clipRectangleToElement(rect, i, j);
        }
      }
      
      /*
       * partial right elements (excluding corners)
       * colIndex jj has to be
       * <= rowIndex
       * -> ii starts at jj, but only if it is larger than rowMinFull
       * (otherwise the corner could be included twice)
       *
       * start with right edge
       * > reason: if the first column is partially filled, we
       * have colMaxFull = -1; colMax = 1 > working
       * but colMinFull = colMin = 0 > not working
       */
      if(colMaxFull < colMax){
        int jj = colMax;
        UInt j = (UInt) colMax;
        
        for(int ii = std::max(rowMinFull,jj); ii <= rowMaxFull; ii++){
          UInt i = (UInt) ii;
          
          sum += clipRectangleToElement(rect, i, j)*preisachWeights_[i][j];
          //  accessCounter[i][j] += clipRectangleToElement(rect, i, j);
        }
      }
      
      /*
       * partial left elements (excluding corners)
       * colIndex jj has to be
       * <= rowIndex
       * -> ii starts at jj, but only if it is larger than rowMinFull
       * (otherwise the corner could be included twice)
       *
       * check if colMin and colMax are the same > in this case we only have
       * one partial col which was already treated above
       */
      if((colMinFull > colMin)&&(colMin != colMax)){
        int jj = colMin;
        UInt j = (UInt) colMin;
        
        for(int ii = std::max(rowMinFull,jj); ii <= rowMaxFull; ii++){
          UInt i = (UInt) ii;
          
          sum += clipRectangleToElement(rect, i, j)*preisachWeights_[i][j];
          //  accessCounter[i][j] += clipRectangleToElement(rect, i, j);
        }
      }
    }
    
    //    std::cout << "rowMin / rowMinFull: " << rowMin << " / " << rowMinFull << std::endl;
    //    std::cout << "rowMax / rowMaxFull: " << rowMax << " / " << rowMaxFull << std::endl;
    //    std::cout << "colMin / colMinFull: " << colMin << " / " << colMinFull << std::endl;
    //    std::cout << "colMax / colMaxFull: " << colMax << " / " << colMaxFull << std::endl;
    //    std::cout << "AccessCounter (end): " << accessCounter.ToString() << std::endl;
    
    return sum;
  }
  
  Double VectorPreisachv10_ListApproach::getRectangleFromSwitchingList(std::list<ListEntryv10>& list,
          std::list<ListEntryv10>::iterator startIt, std::list<ListEntryv10>::iterator curIt, std::list<ListEntryv10>::iterator endIt,
          UInt idArea, Rectangle& rect, bool upperSplitSquare){
    /*
     * This functions is used to decompose a stair-case switching state into rectangular areas which are
     * easier to compute.
     * This functions will only be used for the switching lists
     *
     * upperSplitSquare:
     *  If true, the coordinates of the diagonally split upper square (area 0) will be returned)
     *  (cnt will be set to 0, to trigger the right calculations)
     *
     */
    /*
     * In older version we had to distinguish between upper square and upper triangle region when
     * calling this function. Depending on the region, different outer edges were set.
     * Actually, this was only needed for the upper triangle, as the resulting rectangular
     * bounds were used directly in the later evaluation. For the upper square, the resulting
     * bounds were clipped against the rotation area and thus possible extensions into other
     * regions would have been cut away.
     * In the new version, we clip always against the rotation state and thus we can use the
     * full Preisach plane as outer range. However, we now have to use the former TRIANGLE treatment in
     * each case, i.e. we have to cut the switching areas which go over the diagonal alpha = beta
     * down to triangle, so that the total shape of areas overlapping alpha = beta become squares.
     */
    
    Double alpha, beta, delta;
    alpha = -1.0;
    beta = -1.0;
    delta = 2.0;
    
    /*
     * Calculates the coordinates of a rectangle lying inside the area
     * (alpha, alpha+delta) x (beta, beta+delta)
     *
     * Input:
     *  list -> switching list
     *  startIt -> iterator pointing to the FIRST entry
     *  curIt -> iterator pointing to the CURRENT entry
     *  endIt -> iterator pointing to the LAST entry
     *  alpha,beta -> bottom left corner coordinates of outer area
     *  idArea -> index of area to be computed - 1 (i.e. area1 -> idArea = 0)
     *
     *  tRet,bRet -> y-coordinates of the top and bottom edge of the rectangle
     *  lRet,rRet -> x-coordinates of the left and right edge of the rectangle
     */
    
    /*
     * Splitting into splitting areas for e1 = max
     *
     *                         alpha
     *      __ __ __ __ __ __ __ |_  __ __ __ __ __ __
     *     |\    |          |       |               /
     *     |  \  |    A0    | A2    | A4          /
     *     |_ _ \|_ _ _ _ _ |       |           /___ e1
     *     |A1              |       |         /
     *     |_ _ _ _ _  _ _ _|_ _ __ |       /___ e3
     *     |A3                      |     /
     *     |_ _ _ _ _ _ _ _ _ _ _ _ |_ _/___ e5
     *     |A5                        /
     *     |                        /
     *     |                      /_|__________________ beta
     *     |                    /   |
     *     |                  /     |
     *     |                /       e4
     *     |              / |
     *     |            /   |
     *     |          /     e2
     *     |        /
     *     |     |/
     *     |    /|
     *     |_ /  |
     *          -e1
     *
     *     set e0 = -e1
     *
     *     negative areas (i = even)
     *       L = ei;             R = min(ei+2,1) > if ei+2 does not exist: ei+2 = 1
     *       B = max(ei+1,ei);   T = 1
     *             > limit extend of areas such, that the area below the diagonal is a triangle
     *             > if ei+1 does not exist: ei+1 = -1
     *
     *     positive areas (i = odd)
     *       L = -1;             R = min(ei+1,ei) > ensure that area below diagonal is triangular; if ei+1 does not exist: ei+1 = 1
     *       B = max(ei+2,-1);   T = ei
     *             > if ei+2 does not exist: ei+2 = -1
     *
     *     upper split area:
     *       L = -1; R = -e1; B = e1; T = 1
     *
     *
     * Splitting into splitting areas for e1 = min
     *
     *                         alpha
     *      __ __ __ __ __ __ __ |_  __ __ __ __ __ __
     *     |\    |     |           |                /
     *     |  \  | A1  |   A3      | A5           /
     *     |_ _ \|     |           |           _/___ -e1
     *     |A0   |     |           |          /
     *     |_ _ _|_ _  |           |        /___ e2
     *     |A2         |           |      /
     *     |_ _ ___ _ _|_ _  _ _ _ |   _/___ e4
     *     |A4                     |  /
     *     |                       |/
     *     |                      /|__________________ beta
     *     |                    /  |
     *     |                  /    |
     *     |                /      e5
     *     |              /
     *     |           |/
     *     |          /|
     *     |        /  |
     *     |     |/    e3
     *     |    /|
     *     |_ /  |
     *           e1
     *
     *     set e0 = -e1
     *
     *     negative areas (i = odd)
     *       L = ei;             R = min(ei+2,1) > if ei+2 does not exist: ei+2 = 1
     *       B = max(ei+1,ei);   T = 1
     *             > limit extend of areas such, that the area below the diagonal is a triangle
     *             (take A5 as example: as no e6 exists, A5 would reach down to the bottom of the Preisach plane (-1)
     *              by this, the area below the diagonal would be a trapezoid standing on its side
     *              by restricting to e5, we have only a triangle below the diagonal)
     *             > if ei+1 does not exist: ei+1 = -1
     *
     *     positive areas (i = even)
     *       L = -1;             R = min(ei+1,ei) > ensure that area below diagonal is triangular; if ei+1 does not exist: ei+1 = 1
     *       B = max(ei+2,-1);   T = ei
     *             > if ei+2 does not exist: ei+2 = -1
     *
     *     upper split area:
     *       L = -1; R = e1; B = -e1; T = 1
     *
     *
     * Note regarding different starting case (e1 = min or max):
     *   The rules for calculating positive and negative areas are independent of the type of e1.
     *   The only change is w.r.t. the indices. For e1 = min, positive areas have even indices and for
     *   e1 = max, they have odd indices. However, this is true for the type of ei, too. If list starts
     *   with minimum, all odd entries will be minima and if list starts with maximum, all odd entries will
     *   be maxima. Together with the calculation rules, we see that area Ai is positive, if ei is a maximum
     *   no matter what type e1 is. The calculation can thus be done without checking for indices or the
     *   type of e1. We only have to check, whether ei is a max or a min.
     *   Exception: Area0 is obtained from e1, too. To get this area, idArea has to be set to 0.
     *
     *
     * Note regarding initial values:
     *   When the classical model is used (version 2012), the lower triangle T_L gets not evaluated
     *   via switching lists as its switching state is always +1; to compensate for the +1 triangle, we
     *   should have an initial value of -1 for the upper triangle; this can easily be achieved by
     *   initializing the list with a max(or min) of value 0.
     *
     *
     */
    
    Double l,r,t,b;
    Double nextVal, nextnextVal;
    Double firstVal = startIt->getVal();
    Double curVal = curIt->getVal();
    bool firstMin = startIt->isMin();
    bool curMin = curIt->isMin();
    
    if(upperSplitSquare == true){
      
      l = -1;
      if(firstMin){
        r = firstVal;
        b = -firstVal;
      } else {
        r = -firstVal;
        b = firstVal;
      }
      t = 1;
      
      rect.setBounds(l,r,t,b);
      return 0.0; // value 0.0 not needed, but by this we can check if value was hit
    }
    
    if(idArea == 0){
      /*
       * use startIt regardless of the state of curIt!
       */
      if(firstMin){
        /*
         * List starts with minimum; get area 0 according to
         *
         *       L = -1;             R = e1
         *       B = max(e2,-1);     T = -e1
         *
         *       if e2 does not exist: e2 = -1
         *
         */
        if(startIt != endIt){
          startIt++;
          if(startIt->isMin() == firstMin){
            EXCEPTION("MinMaxList has to be alternating!")
          } else {
            nextVal = startIt->getVal();
          }
          startIt--;
        } else {
          nextVal = alpha;
        }
        
        l = beta;
        r = firstVal;
        b = nextVal; // will automatically be std::max(nextVal,alpha);
        t = -firstVal;
        
        /*
         * min-helper area is positive, so +1.0
         */
        rect.setBounds(l,r,t,b);
        return 1.0;
        
      } else {
        /*
         * List starts with maximum; get area 0 according to
         *
         *
         *       L = -e1;       R = min(e2,1)
         *       B = e1;        T = 1
         *
         *       if e2 does not exist: e2 = 1
         *
         */
        if(startIt != endIt){
          startIt++;
          if(startIt->isMin() == firstMin){
            EXCEPTION("MinMaxList has to be alternating!")
          } else {
            nextVal = startIt->getVal();
          }
          startIt--;
        } else {
          nextVal = beta+delta;
        }
        
        l = -firstVal;
        r = nextVal; // will automatically be std::min(nextVal,beta+delta);
        b = firstVal;
        t = alpha+delta;
        
        /*
         * max-helper area is negative, so -1.0
         */
        rect.setBounds(l,r,t,b);
        return -1.0;
      }
    } else {
      /*
       * area i with i > 0 -> use curIt
       */
      if(curMin){
        /*
         * entry is minimum; get area i according to
         *
         *       L = ei;             R = min(ei+2,1) > if ei+2 does not exist: ei+2 = 1
         *       B = max(ei+1,ei);   T = 1
         *             > limit extend of areas such, that the area below the diagonal is a triangle
         *             > if ei+1 does not exist: ei+1 = -1
         *
         */
        if(curIt != endIt){
          curIt++;
          
          if(curIt->isMin() == curMin){
            EXCEPTION("MinMaxList has to be alternating!")
          } else {
            nextVal = curIt->getVal();
          }
          
          if(curIt != endIt){
            curIt++;
            
            if(curIt->isMin() != curMin){
              EXCEPTION("MinMaxList has to be alternating!")
                      /*
                       * here we would expect the same entry type the the one of the current iteration
                       * e.g. a minimum after a maximum which followed a minimum
                       */
            } else {
              nextnextVal = curIt->getVal();
            }
            curIt--;
          } else {
            nextnextVal = beta+delta;
          }
          curIt--;
        } else {
          nextVal = alpha;
          nextnextVal = beta+delta;
        }
        
        l = curVal;
        r = nextnextVal; // will automatically be std::min(nextnextVal,beta+delta);
        b = std::max(curVal,nextVal);
        t = alpha+delta;
        
        /*
         * min area is negative, so -1.0
         */
        rect.setBounds(l,r,t,b);
        return -1.0;
        
      } else {
        /*
         * entry is maximum; get area i according to
         *
         *       L = -1;             R = min(ei+1,ei) > ensure that area below diagonal is triangular; if ei+1 does not exist: ei+1 = 1
         *       B = max(ei+2,-1);   T = ei
         *             > if ei+2 does not exist: ei+2 = -1
         *
         */
        
        if(curIt != endIt){
          curIt++;
          
          if(curIt->isMin() == curMin){
            EXCEPTION("MinMaxList has to be alternating!")
          } else {
            nextVal = curIt->getVal();
          }
          
          if(curIt != endIt){
            curIt++;
            
            if(curIt->isMin() != curMin){
              EXCEPTION("MinMaxList has to be alternating!")
                      /*
                       * here we would expect the same entry type the the one of the current iteration
                       * e.g. a minimum after a maximum which followed a minimum
                       */
            } else {
              nextnextVal = curIt->getVal();
            }
            curIt--;
          } else {
            nextnextVal = alpha;
          }
          curIt--;
        } else {
          nextVal = beta+delta;
          nextnextVal = alpha;
        }
        
        l = beta;
        r = std::min(curVal,nextVal);
        b = nextnextVal; // will automatically be std::max(nextnextVal,alpha);
        t = curVal;
        
        /*
         * max area is positive, so 1.0
         */
        rect.setBounds(l,r,t,b);
        return 1.0;
      }
    }
  }
  
  bool VectorPreisachv10_ListApproach::Simplify_LocalSwitchingLists(std::list<RotListEntryv10>& usedList){
    /*
     * This function iterates over the globalRotation list and checks each entry in the corresponding
     * local switching list for overlap with the two possible rotation areas.
     * This search is performed until at least one overlap is found. All switching entries BEFORE the
     * one which lead to the first overlap are removed from the list. If at least 1 entry got removed,
     * startCnt will be set to 1 (as the entry of the switching list which corresponds to the helper area
     * (id = 0) no longer is in the list).
     * Leave out all rotListEntries for which the flag isUpdated is false.
     */
    
    std::list<ListEntryv10>::iterator swListIt;
    std::list<ListEntryv10>::iterator firstToKeep;
    std::list<ListEntryv10>::iterator swListEnd; // = --(globSwitchList_[idElem].end());
    std::list<ListEntryv10>::iterator swListStart;
    
    /*
     * iterators for global rotation list
     */
    std::list<RotListEntryv10>::iterator rotListIt;
    std::list<RotListEntryv10>::iterator rotListEnd = --(usedList.end());
    
    bool switchingListSimplified = false;
    bool twoAreas,lastRotListEntry;
    
    Rectangle rotRect1 = Rectangle(0,0,0,0);
    Rectangle rotRect2 = Rectangle(0,0,0,0);
    Rectangle swRect = Rectangle(0,0,0,0);
    Rectangle overlapRect = Rectangle(0,0,0,0);
    
    UInt cnt;
    
    for(rotListIt = usedList.begin(); rotListIt != usedList.end(); rotListIt++){
      
      if(rotListIt == rotListEnd){
        lastRotListEntry = true;
      } else {
        lastRotListEntry = false;
      }
      /*
       * maybe the list was already simplified once or was inherited from a previous rotListState
       * in that case we have to make sure to not test for the wrong value
       */
      cnt = rotListIt->getStartCnt();
      
      /*
       * check if reevaluation is needed at all
       */
      if(rotListIt->hasChanged() == false){
        /*
         * neither switching list did change since last time (and weights did not change either)
         * nor did the lower bound of the rotation area
         * -> rotation area and all switching areas are the same as last time
         * -> overlap of switching areas and rotation area does not have to be computed again
         * -> take stored value
         */
        continue;
      } else {
        /*
         * get rotation area(s)
         */
        
        swListStart = rotListIt->getListReference().begin();
        swListEnd = --(rotListIt->getListReference().end());
        
        twoAreas = getRectanglesFromRotEntry(rotListIt, rotRect1, rotRect2,lastRotListEntry);
        
        firstToKeep = rotListIt->getListReference().begin();
        swListEnd = --(rotListIt->getListReference().end());
        
        bool success1 = false;
        bool success2 = false;
        bool gotSuccess = false;
        
        /*
         * here we do not increase the iterator directly (similar as we do it during evaluation)
         *
         * Reason:
         *  during evaluation, we have to use the first entry twice as it corresponds to area 0 and 1
         *  (as long as this first entry of the list really is the first one and not the first one
         *  after simplification).
         *  Here we do not want to calculate the resulting value but only want to check if there
         *  is an overlap of switching area and rotation area. As the first entry corresponds to two
         *  areas we have to check if one of them has a valid overlap.
         */
        for(swListIt = rotListIt->getListReference().begin(); swListIt != rotListIt->getListReference().end(); ){
          
          /*
           * get rectangular bounds corresponding to entry of switching list (compare to Evaluate_GlobalSwitchingList)
           */
          getRectangleFromSwitchingList(rotListIt->getListReference(),swListStart, swListIt, swListEnd,cnt, swRect);
          
          /*
           * clip resulting area against rotation area1
           */
          success1 = rotRect1.clipRectangles(swRect,overlapRect);
          
          if(twoAreas){
            /*
             * repeat the clipping steps for the second area
             */
            success2 = rotRect2.clipRectangles(swRect,overlapRect);
          }
          
          if((success1 == true)||(success2==true)){
            /*
             * the current switching area has an overlap with at least one of the two rotation areas
             * -> this entry has to be kept
             * -> end loop (as we suppose that all later entries have an overlap, too
             * This assumption is ok as only those entries are appended to the list which actually have
             * and effect, i.e. we only have to cut at the start of the list!
             */
            gotSuccess = true;
            break;
          }
          
          if(cnt > 0){
            
            /*
             * entry has no overlap, check next one
             * (start increasing only if cnt > 0 as list entry corrsponding to cnt == 0
             * has to be checked for area 0 and 1)
             */
            swListIt++;
            firstToKeep++;
          }
          cnt++;
        } // sw list
        
        if(firstToKeep != rotListIt->getListReference().begin()){
          /*
           * At least one entry has to be removed from the list
           */
          if(!gotSuccess){
            /*
             * NEW (30.6.2016)
             * check if we had an success at all!
             * for rotResistance < 1, it is possible, that there are rotation states
             * which have no overlap at all with switching areas (as xPar is too small).
             * In that case, this function will erase the whole switching list as there is
             * no overlap (ok so far). By setting the start cnt to value > 0 we do not
             * check the rotation state against the initially half split upper triangle
             * part (will only be checked for if start cnt = 0 -> see Evaluate_GlobalRotationList).
             * However, in the case of unsymmetric weights, we have to do just that!
             * -> if we had no success at all, keep starting counter at 0, so that the
             * later functions will overlap these rotation areas with the initial half
             * split upper square area.
             *
             */
            rotListIt->setStartCnt(0);
            //std::cout << "No overlap between switching list and rotation list. Will delete whole switching list! Overlap only with initially split area!" << std::endl;
            /*
             * furthermore, we are not allowed to delete the full list!
             * Reason: we need the first entry of the list, to determine the measure of the
             * initally split upper square part
             * -> leave first entry of switching list!
             */
            if(firstToKeep != ++(rotListIt->getListReference().begin())){
              rotListIt->getListReference().erase(++(rotListIt->getListReference().begin()),firstToKeep);
              switchingListSimplified = true;
            }
          } else {
            rotListIt->setStartCnt(1);
            rotListIt->getListReference().erase(rotListIt->getListReference().begin(),firstToKeep);
            switchingListSimplified = true;
          }
        } //else
        /*
         * whole list has to be kept
         * -> continue with next rotListElement
         */
      } // rot list has changed
    } // rot list
    
    return switchingListSimplified;
  }
  
  bool VectorPreisachv10_ListApproach::Simplify_GlobalRotationList(std::list<RotListEntryv10>& usedList){
    /*
     * New merging rule (applicable for all versions)
     *
     * Idea: due to new insertion rule for switching lists (entries get clipped to bounding box
     *       created from rotlist entries), adjacent rotlist entries very seldom have the same switching
     *       list, even if the new entries xPar would override both switching lists
     *
     *       Example (a > b):
     *        rotentry1: bounding box1 left = -a, right = 1, bot = -1, top = a; inner list: max = (a+b)/2; min = b
     *        rotentry2: bounding box2 left = -b, right = 1, bot = -1, top = b; inner list: max = b/2
     *
     *        new entry xPar = b:
     *          rotentry1: bounding box1 left = -a, right = 1, bot = -1, top = a; inner list: max = (a+b)/2;
     *          rotentry2: bounding box2 left = -b, right = 1, bot = -1, top = b; inner list: max = b
     *          -> cannot be merged
     *
     *        new entry xPar = (a+b)/2:
     *          rotentry1: bounding box1 left = -a, right = 1, bot = -1, top = a; inner list: max = (a+b)/2
     *          rotentry2: bounding box2 left = -b, right = 1, bot = -1, top = b; inner list: max = b -> because of clipping
     *                      (without clipping: inner list: max = (a+b)/2 -> could be merged)
     *          -> cannot be merged due to old merging rule in combination with new clipping rule
     *
     * New merging rule:
     *        entry_i_last shall be the last entry of the switching list switchinglist_i belonging to rotlist i
     *        entry_i+1_first shall be the first entry of the switching list switchinglist_i+1 belonging to rotlist i+1
     *
     *        if length(switchinglist_i+1) == 1    <- always the case if switchinglist_i+1 was completely overwritten
     *          if type(entry_i_last) == type(entry_i+1_first):
     *            if(entry_i_last == entry_i+1_first)   <- new list starts where old list ended
     *              OR
     *              if(type(entry_i_last == max)
     *                if(entry_i_last > entry_i+1_first) AND entry_i+1_first == top value of bounding box
     *                  -> merge
     *              else
     *                if(entry_i_last < entry_i+1_first) AND entry_i+1_first == left value of bounding box
     *                  -> merge
     *
     * Idea behind new rule:
     *       - merging in old version was more or less only possible if list of entry i+1 was completely overwritten;
     *        otherwise it is very seldom the case that both switching list can match
     *       - for the old rule, merging could only be done if the lists were the same for all entries;
     *        however, the first list can contain entries which do not affect the later rotation state
     *        (at least due to the new shape of the rotation entries coming from version 10); these entries
     *        will never be inserted into the switching list of the next following rotation state due to the
     *        new clipping and entries which was inherited will be sorted out in the simplify function for the
     *        switching list -> merging becomes nearly impossible for version 10
     *       - the new merging rule solves this problem as we already know that we can only merge if
     *        the switching list of entry i+1 is overwritten and that this will be the case no matter
     *        what entries the list of entry i has, as long as the last entry is a smaller or equal min
     *        /a larger of equal max
     *
     *       Example (a > b):
     *        rotentry1: bounding box1 left = -a, right = 1, bot = -1, top = a; inner list: max = (a+b)/2; min = b
     *        rotentry2: bounding box2 left = -b, right = 1, bot = -1, top = b; inner list: max = b/2
     *
     *        new entry xPar = b:
     *          rotentry1: bounding box1 left = -a, right = 1, bot = -1, top = a; inner list: max = (a+b)/2
     *          rotentry2: bounding box2 left = -b, right = 1, bot = -1, top = b; inner list: max = b
     *          -> can be merged now, as (a+b)/2 > b and b = top
     *
     *        new entry xPar = (a+b)/2:
     *          rotentry1: bounding box1 left = -a, right = 1, bot = -1, top = a; inner list: max = (a+b)/2
     *          rotentry2: bounding box2 left = -b, right = 1, bot = -1, top = b; inner list: max = b -> because of clipping
     *                      (without clipping: inner list: max = (a+b)/2 -> could be merged)
     *          -> can be merged now, as (a+b)/2 > b and b = top
     *
     */
    
    
    /*
     * iterate over rotation list and check if two adjacent entries have the same rotation state and the
     * same list of switching entries
     * This case can happen e.g. if an input overrides multiple older rotation areas at once and the
     * resulting value of xPar is large enough to overwrite the contained switching lists
     */
    bool rotListSimplified = false;
    if(usedList.size() < 2){
      /*
       * list has not enough entries to be merged together
       */
      return false;
    }
    
    std::list<RotListEntryv10>::iterator listIt;
    std::list<RotListEntryv10>::iterator nextListIt;
    std::list<RotListEntryv10>::iterator listEnd = --(usedList.end());
    
    for(listIt = usedList.begin(); listIt != usedList.end(); ){
      
      /*
       *  get next entry in list (if any)
       *  if not -> nothing more to do here
       */
      if(listIt != listEnd){
        nextListIt = listIt;
        nextListIt++;
        
        /*
         * check if entries can be merged together
         * (v2 = new merging rule; without v2 = classical merging)
         */
        if(listIt->canBeMergedWith_v2(*nextListIt,classical_)){
          /*
           * get lower bound of second element and set lower bound of first element to that value
           * (-> i.e. we extend the area of the first element)
           */
          Double low = nextListIt->getLowerVal();
          listIt->setLowerVal(low);
          /*
           * now we can delete the second entry
           */
          usedList.erase(nextListIt);
          
          /*
           * reset iterator to last entry
           * (normally the iterators should remain valid, but I do not trust them ...
           */
          listEnd = --(usedList.end());
          
          /*
           * do not increase iterator
           * -> it might be, that the by now extended first entry matches also the next one
           */
          rotListSimplified = true;
        } else {
          /*
           * increase iterator and check next pair
           */
          listIt++;
        }
      } else {
        break;
      }
    }
    
    return rotListSimplified;
  }
  
  void VectorPreisachv10_ListApproach::mapRectangleToHelperMatrix(Matrix<Double>& helper, Rectangle rect, Double factor, bool skipUpperDiagonal, bool isRotState){
    /*
     * similar function to mapRectangleToPreisachWeights
     * instead of summing up the Preisach weights, we set the entries in the helper matrix
     * -> only needed for output
     *
     * skipUpperDiagonal (see also mapRectangleToPreisachWeights):
     *  needed for the upper square part which is still split into a positive and a negative half
     *  (splitted along alpha = -beta)
     *  Go over this area and set all entries above alpha = -beta to negative values, all values below
     *  alpha = -beta to positive values and entries directly on the diagonal to 0.
     *
     * new flag: isRotState
     * -> if the rotation state is mapped to the helper matrix, we do not want to scale its value by 0.5 if its on the diagonal alpha = beta
     *    if we would do so, we would get the factor 0.5 twice in the final images which are created from the matrix (1. from switching state,
     *    2. from rotation state)
     * -> if isRotState == True, scale full entries on alpha = beta times 2 and for partially filled elements, skip
     *    the cutting of triangle parts below the diagonal
     *
     */
    
    Double l,r,t,b;
    rect.getBounds(l,r,t,b);
    
    /*
     * restrict input value to valid range
     * (Preisach plane just goes from -1 to +1 in alpha and beta)
     */
    t = std::min(t,1.0);
    b = std::max(b,-1.0);
    l = std::max(l,-1.0);
    r = std::min(r,1.0);
    
    /*
     * check if rectangular has size != 0
     */
    if((t <= b)||(r <= l)){
      return;
    }
    
    UInt numRows = helper.GetNumRows();
    
    Double delta = 2.0/numRows;
    
    /*
     * Lower triangluar part of Preisach plane (here 16 elements)
     * and overlapping rectangle
     *   _____ _____ _____ _____ _____
     *  |    .|.....|.....|.....|.   /
     *  | 1 xx | 2   | 3   | 4   |xx5/
     *  |___xx_|_____|_____|_____|/
     *  |   xx |     |     |    / xx
     *  | 6 xx | 7   | 8   | 9/   xx
     *  |___xx_|_____|_____|/     xx
     *  |   xx |     |    /       xx
     *  | a xx.|.b...|.c/.........xx
     *  |_____|_____|/
     *  |     |    /
     *  | d   | e/
     *  |_____|/
     *  |    /
     *  | f/
     *  |/
     *
     *  Element 7,8,9 are completely overlapped
     *  Element 6 is cut vertically
     *  Element 2,3,b,c are cut horizontally
     *  Element a is cut horizontally and vertically
     *  Element c,5 are cut horizontally, vertically and diagonally
     *
     *  Approach:
     *   Iterate over all (fully and partially overlapped) elements
     *    >Check if indices denote a fully overlapped element
     *      > if yes sum up value
     *      > else use clipToElement to determine the appropriate scaling value for weight, then sum up
     */
    
    /*
     * 1. get indices of all overlapped elements (partially and fully)
     * NOTE: element alpha = -1 and beta = -1 will have index 0,0!
     * -> add floor(1.0/delta_) != numRows/2
     */
    int rowMin =  std::floor(b/delta) + std::floor(1.0/delta);
    int rowMax = std::ceil(t/delta) + std::floor(1.0/delta)-1;
    
    int colMin =  std::floor(l/delta) + std::floor(1.0/delta);
    int colMax =  std::ceil(r/delta) + std::floor(1.0/delta)-1;
    
    /*
     * NEW: use integers instead of unsigned integer
     * REASON: colMaxFull became negative (which simply would indicate no fully overlapped elemets), but
     * unsigning took it to a very large positive value and thus -> wrong result
     */
    
    /*
     * 2. get indices of completely overlapped elements
     */
    int rowMinFull =  std::ceil(b/delta) + std::floor(1.0/delta);
    int rowMaxFull = std::floor(t/delta) + std::floor(1.0/delta)-1;
    
    int colMinFull =  std::ceil(l/delta) + std::floor(1.0/delta);
    int colMaxFull = std::floor(r/delta) + std::floor(1.0/delta)-1;
    
    
    if(skipUpperDiagonal == true){
      /*
       * special treatment:
       * 1. skip entries on alpha = -beta
       * 2. consider the right signs, i.e. +1 for entries below alpha = -beta
       *    and -1 for entries above!
       */
      Double sign;
      for(int ii = rowMin; ii <= rowMax; ii++){
        /*
         * ensure, that j <= i -> other elements do not contribute
         */
        UInt i = (UInt) ii;
        
        for(int jj = colMin; jj <= std::min(ii,colMax); jj++){
          
          UInt j = (UInt) jj;
          
          if(numRows-i-j-1 == 0){
            /*
             * alpha = -beta
             */
            helper[i][j] = 0.0;
            continue;
          } else if (numRows < i+j+1){
            /*
             * above alpha = -beta
             * -> sign = -1
             */
            sign = -1.0;
          } else {
            /*
             * below diagonal
             * -> sign = +1
             */
            sign = 1.0;
          }
          
          /*
           * Check for an inner element
           */
          if((ii >= rowMinFull)&&(ii <= rowMaxFull)){
            if((jj >= colMinFull)&&(jj <= colMaxFull)){
              /*
               * check for diagonal element on alpha = beta
               * (should not happen here)
               */
              if((i == j)&&(isRotState==false)){
                helper[i][j] = sign*0.5;
              }
              
              /*
               * scale Preisach weights with area
               */
              helper[i][j] = sign*1.0;
            } else {
              /*
               * Partially overlapped element
               * -> clip rectangle to element
               * Note: it does not matter here, if we pass the full input area or the
               * subarea defined by the current pair i,j as the clip function will overlap
               * the element area with the provided input area and automatically find the
               * correct overlap; the clipToElement function will furthermore check for
               * areas below the diagonal alpha = beta and cut them from the overlapping area
               *
               * Note 2: divide by delta^2 to normalize the values to range -1 to +1
               */
              /*
               * several element can contribute to the same helper entry -> sum up instead of setting
               */
              helper[i][j] += sign*clipRectangleToElement(rect, i, j, delta,isRotState)/(delta*delta);
            }
          } else {
            /*
             * Partially overlapped element
             * -> clip rectangle to element
             * Note: it does not matter here, if we pass the full input area or the
             * subarea defined by the current pair i,j as the clip function will overlap
             * the element area with the provided input area and automatically find the
             * correct overlap; the clipToElement function will furthermore check for
             * areas below the diagonal alpha = beta and cut them from the overlapping area
             *
             * Note 2: divide by delta^2 to normalize the values to range -1 to +1
             */
            helper[i][j] += sign*clipRectangleToElement(rect, i, j, delta,isRotState)/(delta*delta);
          }
        }
      }
      
    } else {
      /*
       * std treatment
       */
      
      /*
       * Iterate over all elements
       */
      for(int ii = rowMin; ii <= rowMax; ii++){
        /*
         * ensure, that j <= i -> other elements do not contribute
         */
        UInt i = (UInt) ii;
        
        for(int jj = colMin; jj <= std::min(ii,colMax); jj++){
          /*
           * Check for an inner element
           */
          UInt j = (UInt) jj;
          
          if((ii >= rowMinFull)&&(ii <= rowMaxFull)){
            if((jj >= colMinFull)&&(jj <= colMaxFull)){
              /*
               * check for diagonal element
               */
              if((i == j)&&(isRotState==false)){
                helper[i][j] = factor*0.5;
              } else {
                helper[i][j] = factor*1.0;
              }
            } else {
              /*
               * Partially overlapped element
               * -> clip rectangle to element
               * Note: it does not matter here, if we pass the full input area or the
               * subarea defined by the current pair i,j as the clip function will overlap
               * the element area with the provided input area and automatically find the
               * correct overlap; the clipToElement function will furthermore check for
               * areas below the diagonal alpha = beta and cut them from the overlapping area
               *
               * Note 2: divide by delta^2 to normalize the values to range -1 to +1
               */
              /*
               * several element can contribute to the same helper entry -> sum up instead of setting
               */
              helper[i][j] += factor*clipRectangleToElement(rect, i, j, delta,isRotState)/(delta*delta);
            }
          } else {
            /*
             * Partially overlapped element
             * -> clip rectangle to element
             * Note: it does not matter here, if we pass the full input area or the
             * subarea defined by the current pair i,j as the clip function will overlap
             * the element area with the provided input area and automatically find the
             * correct overlap; the clipToElement function will furthermore check for
             * areas below the diagonal alpha = beta and cut them from the overlapping area
             *
             * Note 2: divide by delta^2 to normalize the values to range -1 to +1
             */
            helper[i][j] += factor*clipRectangleToElement(rect, i, j, delta,isRotState)/(delta*delta);
          }
        }
      }
    }
  }
  
  
  void VectorPreisachv10_ListApproach::switchingStateToBmp(UInt numPixel, std::string filename, UInt idElem, bool overLayWithRotState)
  {
    /*
     * NEW: rotation state is evaluated along with the switching states if overLayWithRotState is true
     * in the old versions, the rotation state was evaluated separately although we have to iterate over the
     * rotation list to evaluate the switching lists
     */
    
    if(numPixel < 2){
      WARN("Image should have more than 2 x 2 pixel");
      return;
    }
    
    if(numPixel%2 != 0){
      WARN("Rounded number of pixel ("<<numPixel<<") to a multiple of 2 ("<<numPixel+1<<")");
      numPixel = numPixel + 1;
    }
    
    /*
     * create matrix needed to save switching state
     */
    Matrix<Double> helperMatrix = Matrix<Double>(numPixel,numPixel);
    helperMatrix.Init();
    
    /*
     * create matrices to store the rotation state
     */
    Matrix<Double> rotX, rotY, rotZ;
    if(overLayWithRotState){
      rotX = Matrix<Double>(numPixel,numPixel);
      rotY = Matrix<Double>(numPixel,numPixel);
      rotX.Init();
      rotY.Init();
      if(dim_ == 3){
        rotZ = Matrix<Double>(numPixel,numPixel);
        rotZ.Init();
      }
    }
    
    /*
     * Fill matrix / matrices
     * 1. if(classical_)
     *  i. clip lowerTriangle to helperMatrix
     *
     * 2. iterate over rotation list
     *  i. store rotation direction if needed
     *  ii. iterate over switching list
     *    a. determine overlapping areas and clip them against helperMatrix
     *    b. for unsymmetric weights include triagonal areas
     */
    
    /*
     * Part 1
     */
    if(classical_){
      for(UInt i = 0; i < numPixel/2; i++){
        for(UInt j = 0; j <= i; j++){
          if(i == j) helperMatrix[i][j] = 0.5;
          else helperMatrix[i][j] = 1.0;
        }
      }
      
      if(overLayWithRotState){
        for(UInt i = 0; i < numPixel/2; i++){
          for(UInt j = 0; j <= i; j++){
            rotX[i][j] = lastEu_[idElem][0];
          }
        }
        for(UInt i = 0; i < numPixel/2; i++){
          for(UInt j = 0; j <= i; j++){
            rotY[i][j] = lastEu_[idElem][1];
          }
        }
        if(dim_ == 3){
          for(UInt i = 0; i < numPixel/2; i++){
            for(UInt j = 0; j <= i; j++){
              rotZ[i][j] = lastEu_[idElem][2];
            }
          }
        }
      }
    }
    
    /*
     * Part 2
     */
    
    /*
     * iterators for local switching list
     */
    std::list<ListEntryv10>::iterator swListIt;
    std::list<ListEntryv10>::iterator swListStart;
    std::list<ListEntryv10>::iterator swListEnd;
    
    /*
     * iterators for global rotation list
     */
    std::list<RotListEntryv10>::iterator rotListIt;
    std::list<RotListEntryv10>::iterator rotListEnd = --(globRotList_[idElem].end());
    
    bool area0 = true;
    bool twoAreas = true;
    bool lastRotListEntryv10;
    Double upperBound;
    
    Rectangle rotRect1 = Rectangle(0,0,0,0);
    Rectangle rotRect2 = Rectangle(0,0,0,0);
    Rectangle swRect = Rectangle(0,0,0,0);
    Rectangle overlapRect = Rectangle(0,0,0,0);
    
    Vector<Double> curRotState;
    Double factor;
    UInt cnt = 0;
    
    for(rotListIt = globRotList_[idElem].begin(); rotListIt != globRotList_[idElem].end(); rotListIt++){
      
      if(rotListIt == rotListEnd){
        lastRotListEntryv10 = true;
      } else {
        lastRotListEntryv10 = false;
      }
      
      curRotState = rotListIt->getVecReference();
      upperBound = rotListIt->getVal(); //exi
      
      swListStart = rotListIt->getListReference().begin();
      swListEnd = --(rotListIt->getListReference().end());
      
      twoAreas = getRectanglesFromRotEntry(rotListIt, rotRect1, rotRect2,lastRotListEntryv10);
      
      if(twoAreas == false){
        /*
         * area0 is already an actually set area (i.e. it has a rotation state
         * -> no special treatment needed!
         */
        area0 = false;
      }
      
      if(area0 == true){
        /*
         * area0 has no rotation state but has a switching state that we want to output
         */
        
        if(classical_){
          /*
           * area 0 consists only of one square region
           */
          Rectangle area0 = Rectangle(-1.0,-upperBound,1.0,upperBound);
          
          /*
           * map rectangle to HelperMatrix and set flag skipUpperDiagonal to true
           * -> by setting the flag to true, the function will automatically assume that the area shall
           * be filled with +1 up to the splitting diagonal and -1 above it
           */
          mapRectangleToHelperMatrix(helperMatrix,area0,0,true);
        } else {
          /*
           * in the revised version, area 0 has the same L-kind shape as the other rotation states;
           * however, to get it split properly, we divide this area into three parts
           * a) square part which is split into +1 and -1 -> same as in classical model
           * b) left, vertical rectangle -> completely set to +1
           * c) upper, horizontal rectangle -> completely set to -1
           *
           */
          Rectangle area0_square = Rectangle(-1.0,-upperBound,1.0,upperBound);
          mapRectangleToHelperMatrix(helperMatrix,area0_square,0,true);
          
          Rectangle area0_left = Rectangle(-1.0,-upperBound,upperBound,-1.0);
          mapRectangleToHelperMatrix(helperMatrix,area0_left,1.0,false);
          
          Rectangle area0_top = Rectangle(-upperBound,1.0,1.0,upperBound);
          mapRectangleToHelperMatrix(helperMatrix,area0_top,-1.0,false);
        }
        
        area0 = false;
        
        /*
         * Note: we do not have to write to the matrix for the rotation states as we have no rotation state here
         */
      }
      
      /*
       * reset counter (needed to check if we have the first entry of the switch list)
       * Regarding cnt:
       *  the absolute value of cnt is not relevant. We only have to check if it is 0 or not.
       *  In case that it is 0 we have to use the first switching entry twice (once for area with
       *  id = 0 and once for area with id = 1). Normally (and in older versions), we start at 0.
       *  However, the lists started to get very long with lots of entries which had no influence on
       *  the result (i.e. the switching areas did not overlap with the rotation area). To reduce the
       *  computational cost and to shorten the lists, the function simplify_switchinglist was introduced.
       *  This function checks directly after the merging step in update_globalrotationlist, all entries
       *  for a possible overlap. Entries at the beginning of the list, which do not lead to an overlap
       *  can be removed from the list (elements at the end which have no influence will not get inserted
       *  in the first place). The first entry of this shortened list does not correspond to area id 0
       *  anymore (as this is the special helper area). So we startCnt in rotListEntryv10 to mark that we
       *  deleted the first entry which originally corresponded to area 0. As all areas with id > 0 are
       *  calculated by the same scheme, the absolute value of cnt is not of interest and thus it is
       *  either 0 or 1.
       */
      cnt = rotListIt->getStartCnt();
      
      bool success = false;
      
      for(swListIt = rotListIt->getListReference().begin(); swListIt != rotListIt->getListReference().end(); ){
        
        /*
         * get rectangular bounds corresponding to entry of switching list (compare to Evaluate_GlobalSwitchingList)
         */
        factor = getRectangleFromSwitchingList(rotListIt->getListReference(),swListStart, swListIt, swListEnd,cnt, swRect);
        
        /*
         * clip resulting area against rotation area1
         */
        success = rotRect1.clipRectangles(swRect,overlapRect);
        
        if(success){
          /*
           *  clip overlap to HelperMatrix (factor holds the value +1 or -1 and indicates how the matrix shall be filled)
           */
          mapRectangleToHelperMatrix(helperMatrix,overlapRect,factor);
        }
        
        if(twoAreas){
          /*
           * repeat the clipping steps for the second area
           */
          success = rotRect2.clipRectangles(swRect,overlapRect);
          
          if(success){
            /*
             *  clip overlap to HelperMatrix (factor holds the value +1 or -1 and indicates how the matrix shall be filled)
             */
            mapRectangleToHelperMatrix(helperMatrix,overlapRect,factor);
          }
        }
        
        /*
         * extra treatment for unsymmetric weights
         * here we have to overlap with the diagonally split area, too.
         * (area 0, area id = -1)
         *
         *      split area = area 0
         *       _v_______________
         *      | \  | 1 |    |   |
         *      |___\|___| 3  |   |
         *      |   2    |    |5  |
         *      |________|____|   |
         *      |     4       |   |
         *      |_____________|___|
         *      |        6        |
         *      |_________________|
         *
         * for symmetric weights, this area will always cancel out as positive and negative
         * part are equal. For unsymmetric weights, we have to evaluate, but only if this special
         * area has not been wiped out
         * Has only be checked once for each rotation state -> cnt == 0
         * (Note that cnt can start at values > 0, if the switching list was simplified.
         * In such a case, we can be sure that the split area does not overlap with a rotation state
         * as the areas 1 and 2 do not overlap either.)
         *
         */
        if(cnt == 0){
          
          /*
           * only difference to evaluation algorithm: do this also, if list was already wiped
           */
          if(isSymmetric_ == false){
            /*
             * set flag upperSplitSquare to true -> get area 0 instead of area 1
             */
            factor = getRectangleFromSwitchingList(rotListIt->getListReference(),swListStart, swListIt, swListEnd,cnt, swRect,true);
            
            if(factor != 0){
              EXCEPTION("Something got wrong here! Check function getRectangleBounds!");
            }
            
            /*
             * clip resulting area against rotation area1
             */
            success = rotRect1.clipRectangles(swRect,overlapRect);
            
            if(success){
              /*
               *  clip overlap to HelperMatrix (factor holds the value +1 or -1 and indicates how the matrix shall be filled)
               * -> NOTE: here we set the flag skipUpperDiagonal to true!
               */
              mapRectangleToHelperMatrix(helperMatrix,overlapRect,factor,true);
            }
            
            if(twoAreas){
              /*
               * repeat the clipping steps for the second area
               */
              success = rotRect2.clipRectangles(swRect,overlapRect);
              
              if(success){
                /*
                 *  clip overlap to HelperMatrix (factor holds the value +1 or -1 and indicates how the matrix shall be filled)
                 * -> NOTE: here we set the flag skipUpperDiagonal to true!
                 */
                mapRectangleToHelperMatrix(helperMatrix,overlapRect,factor,true);
              }
            }
          }
        }
        
        if(cnt > 0){
          /*
           * NOTE: area1 and area2 are both calculated using the first list entry, therefore we do not increase the iterator after
           * the first iteration
           */
          swListIt++;
        }
        cnt++;
      } // sw list
      
      if(overLayWithRotState){
        /*
         * write information about rotation state into matrices
         */
        Double currentEntry_x = curRotState[0];
        
        mapRectangleToHelperMatrix(rotX,rotRect1,currentEntry_x,false,true);
        
        if(twoAreas){
          mapRectangleToHelperMatrix(rotX,rotRect2,currentEntry_x,false,true);
        }
        
        Double currentEntry_y = curRotState[1];
        
        mapRectangleToHelperMatrix(rotY,rotRect1,currentEntry_y,false,true);
        
        if(twoAreas){
          mapRectangleToHelperMatrix(rotY,rotRect2,currentEntry_y,false,true);
        }
        
        if(dim_ == 3){
          Double currentEntry_z = curRotState[2];
          
          mapRectangleToHelperMatrix(rotZ,rotRect1,currentEntry_z,false,true);
          
          if(twoAreas){
            mapRectangleToHelperMatrix(rotZ,rotRect2,currentEntry_z,false,true);
          }
        }
      }
    } // rot list
    
    /*
     * now call output function of matrix
     */
    UInt upscaling = 1;
    
    if(overLayWithRotState == true){
      UInt version = 2;
      
      if(version == 1){
        
        for(UInt comp = 0; comp < dim_; comp++){
          
          std::stringstream stream;
          std::string filename_new;
          if(comp == 0){
            stream << "x-" << filename;
            filename_new = stream.str();
            helperMatrix.matrix2Bmp(upscaling,filename_new,&rotX);
          } else if(comp == 1){
            stream << "y-" << filename;
            filename_new = stream.str();
            helperMatrix.matrix2Bmp(upscaling,filename_new,&rotY);
          } else {
            stream << "z-" << filename;
            filename_new = stream.str();
            helperMatrix.matrix2Bmp(upscaling,filename_new,&rotZ);
          }
        }
        
      } else if(version == 2) {
        /*
         * New way of outputting matrix:
         *   encode rotation state in color: 0 = red; 120 = blue; 240 = green; angles in between colored as a mix
         *   encode switching state as sign and amplitude: negative switching -> angle + 180; absvalue < 1 -> scale final colorcombination
         *
         *   -> only for 2D rotstates, as the z component is not considered
         */
        
        std::stringstream stream;
        stream << "xy-" << filename;
        std::string filename_new = stream.str();
        
        helperMatrix.matrix2Bmp_v2(upscaling,filename_new,&rotX,&rotY);
        //
        //        std::cout << "switching: \n" << helperMatrix.ToString() << std::endl;
        //        std::cout << "rotx: \n " << rotX.ToString() << std::endl;
        //        std::cout << "roty: \n" << rotY.ToString() << std::endl;
      }
    } else {
      helperMatrix.matrix2Bmp(upscaling,filename,NULL);
    }
  }
  
}//end namespace
