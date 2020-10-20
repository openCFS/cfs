#ifndef FILE_VECPREISACH_M_2018
#define FILE_VECPREISACH_M_2018

#include "Hysteresis.hh"
#include "Preisach.hh"

#include <list>

#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"


namespace CoupledField {

  class VectorPreisachMayergoyz : public Hysteresis
  {
  public:
    // general constructor for anisotropic case, i.e. in each spatial direction
    // we have to specify xSaat,ySat as well as a set of fitting preisachWeights
    // the number of directions is specified by the length of xSat, ySat, etc
    VectorPreisachMayergoyz(Integer numElem, Vector<Double> xSat, Vector<Double> ySat,
      Matrix<Double>* preisachWeight, UInt dim, bool isVirgin,
      Vector<Double> anhyst_A, Vector<Double> anhyst_B, Vector<Double> anhyst_C, Vector<Double> anhyst_D,
      bool anhystOnly, int clipOutput);

    // constructor for isotropic case, i.e. same xSat, ySat and weights in all directions
    VectorPreisachMayergoyz(Integer numElem, ParameterPreisachOperators& operatorParams,
          ParameterPreisachWeights weightParams, UInt dim, bool isVirgin);

//      Integer numElem, UInt numDirections, Double xSat, Double ySat,
//      Matrix<Double>& preisachWeight, UInt dim, bool isVirgin,Double anhyst_A, Double anhyst_B, Double anhyst_C, bool anhystOnly, int clipOutput);
//
    virtual ~VectorPreisachMayergoyz();

    //! Try to compute input xVal to hyst operator, such that mu*xVal + H(xVal) = yVal
    // return usable input xVal
    /*
     * computeInput_vec is the one to be called in coefFunctionHyst
     * Exception: testInversion > here we use computeInput_vec_withStatistics
     */
    Vector<Double> computeInput_vec(Vector<Double> yVal, Integer operatorIndex,
      Matrix<Double> mu, bool fieldsAlignedAboveSat,
      bool hystOutputRestrictedToSat, int& successFlag, bool useEverett=false, bool overwrite=false){

      if(useEverett == true){
        return computeInput_vec_Everett(yVal, operatorIndex, mu, overwrite, successFlag);
      } else {
        Vector<Double> prevYval = Vector<Double>(dim_);
        mu.Mult(prevXVal_[operatorIndex],prevYval);
        prevYval.Add(1.0,prevHVal_[operatorIndex]);

        return computeInput_vec_withPrevStates(yVal, prevYval,
          prevXVal_[operatorIndex], prevHVal_[operatorIndex],
          operatorIndex, mu, fieldsAlignedAboveSat,
          hystOutputRestrictedToSat, successFlag);
      }
    }

    Vector<Double> computeInput_vec_Everett(Vector<Double> yVal, Integer operatorIndex,
          Matrix<Double> mu, bool overwrite, int& successFlag);
    
    Vector<Double> computeValue_vec(Vector<Double>& xVal, Integer idx, bool overwrite,
      bool debugOutput, int& successFlag, bool skipAnhystPart = false);

//    Vector<Double> computeValue_vecMeasure(Vector<Double>& xVal, Integer idx, bool overwrite,
//      bool debugOutput, int& successFlag, Double& time);

    void setFlags(UInt performanceFlag){
      ;
    };
    
    Double evaluateLagCorrectionAngle(Vector<Double>& xVal, Vector<Double>& prevXVal, Double maxAngle);

  private:
    // for 2D > only semicircle required due to symmetry
    // used coordinate transform:
    //  x = r cos(phi)
    //  y = r sin(phi)
    UInt numDirections_; // = numDirectionsPhi_ in 2D
    Double deltaAngle_;

    // for 3D > only one hemisphere required due to symmetry
    // used coordinate transform:
    // x = r cos(phi) sin(theta)
    // y = r sin(phi) sin(theta)
    // z = r cos(theta)
    UInt numDirectionsTheta_; // inner angle; from 0 to pi/2
    UInt numDirectionsPhi_; // outer angle from 0 to 2*pi

    Double deltaAnglePhi_;
    Double deltaAngleTheta_;
    Vector<Double> sinTheta_;

    int clipOutput_;
    bool isIsotropic_;

    Vector<Double> startingAxis_;
    Vector<Double>* singleDirections_;
    Preisach** singlePreisachOperators_;
    Matrix<Double> matrixForCoefComputation_;
    
    // additional parameter for fitting rotational loss
    // source: Dlala - "Improving Loss Properties of the Mayergoyz Vector Hysteresis Model"
    Double lossParam_a_;
    Double lossParam_b_;
    bool improveRotLoss_;
  };

} //end of namespace


#endif
