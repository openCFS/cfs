#ifndef FILE_PREISACH_2004
#define FILE_PREISACH_2004

#include "Hysteresis.hh"

#include <list>

#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"


namespace CoupledField {
  
  class Preisach : public Hysteresis
  {
  public:
    //! constructor
    Preisach(Integer numElem, Double xSat, Double ySat, 
      Matrix<Double>& preisachWeight, bool isVirgin);
    
    Preisach(Integer numElem, Double xSat, Double ySat, 
      Matrix<Double>& preisachWeight, bool isVirgin, Double anhyst_A, Double anhyst_B, Double anhyst_C, bool anhystOnly);
    
    //! default constructor
    Preisach();
    
    //! destructor
    virtual ~Preisach();
    
    //! actually never used
    //Double computeValue(Double& xVal, Integer idxElem, bool overwrite = true);
    
    Double computeInputAndUpdate(Double Yin, Double eps_mu, Integer idx, bool overwrite = true);
    
    //!computes for xVal a new output value and deletion rules are applied
    Double computeValueAndUpdate(Double xVal, Integer idxElem, 
    bool overwrite = true);
    
    //! returns the current output of the hyst-operator for element idxElem
    Double getValue(  Integer idxElem);
    
    //! returns the number of stored minima and maxima
    UInt getStringLength( Integer idxElem ) {
      return  StringLength_[ idxElem];
    };
    
    //! updates the list of minima and maxima due to new input
    Double updateMinMaxList(Double newX, Integer idxElem, 
    bool overwrite);
    
    //!
    void SetTimeStepVal(Double dt) 
    {;};
    
    //! normalizes the input to Xsaturated_
    Double normalizeAndClipInput(Double xInput);
    
    std::string runtimeToString(){
      return "No runtime information available for Scalar Preisach model";
    };
    
    void setFlags(UInt performanceFlag){
      ;
    };
    
  protected:
    
    //! computes  the everett function (area-integration for x1, x2)
    Double everettPixel(Double x1, Double x2);
    
//    inline Double evalAnhystPart_normalized(Double xNormalizedUnclipped){
//      // returns normalized anhysteretic part
//      return anhyst_A_*std::atan(anhyst_B_*xNormalizedUnclipped) + anhyst_C_*xNormalizedUnclipped;
//    }
    
    //Double bisectForSaturation(Double Yin, Double eps_mu, Double tol, bool negSaturation);
    
    Double bisect(Double dY,Double xMin,Double xMax, Double xFixed, Double eps_mu, Double tol);
    
    Double XSaturated_; //! saturation value for  input
    Double YSaturated_; //! saturation value for output
    
    /*
     * for optional anhysteretic parts
     */
//    Double anhyst_A_;
//    Double anhyst_B_;
//    Double anhyst_C_;
    
    bool isVirgin_; //! yes, if starting at zero
    
    Vector<Double> preisachSum_; //! output value of Preisach operator
    
    Vector<Double>* strings_; //! irreduceable minima and maxima
    Vector<Double>* helpStrings_; //! help array for string_
    Vector<Integer>* minmaxtype_; // stores for each entry of the min/max list if it is a minimum or not
    // -1 = min; 1 = max; 0 = initial state
    Vector<Double>* evaluatedEverettPixel_; // stores to each entry of the min/max list
    // the corresponding everett pixel (with sign)
    // i.e. everett(-strings[0],strings[0]), everett(strings[0],strings[1]), everett(strings[1],strings[2]), ...
    // note that we need one additional entry for the first min/max
    
    Vector<UInt> StringLength_; //! number of irreduceable minima and maxima
    UInt maxStringLength_; //! maximum allowd length for 
    
    Matrix<Double> preisachWeights_; //! preisach weight function
    
    Double tol_; //! accuracy parameter
    
    // previous input value X and polarization P
    // NOTE:
    //  X and P are normalized; P is clipped to saturation, X must not be clipped
    Vector<Double> previousXval_;
    Vector<Double> previousPval_;
  };
  
  /*
   * 13.4.2018 - Extension to Preisach model
   * 
   * Background and Idea:
   * 
   *  Scalar Preisach model computes the polarization P of a material into a fixed
   *  direction \vec{e}_P which is delivered via mat file.
   *  I.e. it computes
   *    P(\vec{E}\cdot\vec{e}_P) * \vec{e}_P
   * 
   *  The Vector Preisach model as presented by Sutor (2012, 2015; implemented in 
   *  class VectorPreisachv10) extends the scalar model by the following steps:
   *  1. Decompose \vec{E} into a direction \vec{e}_E and its amplitude \|E\|.
   *  2. Update an additional Preisach-like plane using \|E\| as setting value;
   *      each point in the Preisach-like plane stores a corresponding direction
   *      \vec{e}_E_k, depending on the setting history.
   *  3. FOR EACH stored direction vector \vec{e}_E_k / for each subarea k inside 
   *      this additional Preisach-like plane, store a separate Scalar Preisach 
   *      model that gets set by \vec{E}_current * \vec{e}_E_k; these Scalar
   *      Preisach models get evaluated and result in the weighting w_k for \vec{e}_E_k
   *  4. The vector output is computed via \sum_k w_k \vec{e}_E_k.
   * 
   *  Both the Scalar and the Vector model have certain pos and cons:
   *    Scalar          vs.           Vector
   *  + relatively fast           + allows direction changes
   *  + more robust**             + closer to physics / real world applications
   *  - direction fix             - much harder to invert
   *                              - problematic if input with large amplitude 
   *                                  oscillates in direction**
   * 
   * **: if D and P do not align in direction, the missing part D_rest stands
   *      perpendicular on P; D_rest is fully determined by eps0*E_rest which
   *      causes E_rest to become immense even if D_rest is rather small; the
   *      Vector model would immediately rotate towards E_rest due to its large
   *      value which most probably will lead to an overshoot, i.e. during the
   *      next step, E_rest will be very large again but in opposite direction
   *      > solution will never settle
   *      > this is the reason why you never should clip the direction of the
   *          Vector Modell as this would reduce the set of direction which 
   *          could be mapped by it > the Vector model has to be as continuous
   *          as possible!
   *      The Scalar model on the other hand does not care about E_rest at all
   *      as this part is not even passed to the model due to the projection of
   *      \vec{E} onto \vec{e}_P.
   * 
   * 
   * New idea - Extended Preisach
   *  To get the best of both worlds (stability and fastness of Scalar Preisach,
   *  flexibility of Vector Preisach) the following pseudo-approach might be helpful.
   * 
   *  1. Define initial direction \vec{e}_P^0 via mat file
   *  2. Compute \vec{P} as for standard Scalar Preisach model, i.e.
   *        vec{P} = ScalarPreisach( \vec{E} \cdot \vec{e}_P^0 ) * \vec{e}_P^0
   *  3. At the beginning of the next time step n, compute new direction 
   *      \vec{e}_P^n using the flux quantity \vec{D} as setting quantity in 
   *      combination with a Preisach-like plane like it is used in the Vector
   *      Preisach model.
   *  
   *  > Very similar to Vector Preisach model but with three major differences
   *    I. Rotation direction only updated once during each time step (or once
   *        in a while)
   *    II. Use flux quantity to set the rotation direction
   *    III. Rotation state and switching state are fully decoupled
   *        > Rotation states / Preisach-like plane gets set as in Vector model
   *          but instead of computing a scalar model for each subarea / for
   *          each stored direction vector to get the weights w_k for each direction
   *          we take the actual area inside the plane as weight;
   *        > the resulting final direction is used for the computation of a single
   *          Scalar Preisach model
   * 
   *  Motivation behind some important choices:
   *    a) why take flux quantity to determine the direction?
   *      > as described above for the Scalar model, the fixed direction during
   *        computation will most probably prohibit that D and P are aligned;
   *        this leads to very large values of E into the remaining perpendicular
   *        direction; the overall value of E is thus dominated by the part pointing
   *        perpendicular to the current direction of P and thus would lead to
   *        a drastic change in direction even if D is only slightly misaligned with
   *        P
   *      > if we use D to set the directio of P this can hopefully be avoided; it
   *        also makes more sense to take the flux quantity as it usually dominated
   *        by P
   *    b) if we just want to update the rotation state once in a while, why do we
   *        not use the Vector Model and lock its direction?
   *      > this was actually tested but the issue lies in the strong connection between
   *        switching states and rotation states; even if we fix the rotation states
   *        by keeping the Preisach-like plane fixed, the overall direction of the
   *        model will still change when the switching states do as those determine
   *        the weightings for the stored rotation directions
   *      > another reason was mentioned above: if P is not allowed to follow D, 
   *        E will become very large in the perpendicular direction to P; this
   *        will sooner or later cause P to fully rotate into wrong directions; by
   *        decoupling rotation and switching states, we can profit from the 
   *        projection of the input into the actual direction as it is the case in 
   *        the simple Scalar model
   *        
   *  LONG STORY SHORT - what's to be added
   *    1. add simple list of pairs<Double, Vector<Double> >
   *        > list stores rotation directions according to amplitude of corresponding
   *          flux vector, i.e. \vec{e}_D^k and \|D^k\|
   *        > list follows wiping out rules of classical Everett function
   *    2. add functions to update and evaluate the list
   * 
   *    All other functions can be taken from the Scalar model   
   * 
   */
  class ExtendedPreisach : public Preisach
  {
  public:
    ExtendedPreisach(Integer numElem, Double xSat, Double ySat, 
      Matrix<Double>& preisachWeight, Double rotationalResistance , Double angularDistance, UInt dim, bool isVirgin, 
      Double anhyst_A, Double anhyst_B, Double anhyst_C, bool anhystOnly);
    
    virtual ~ExtendedPreisach();
    
    void UpdateRotationState(Vector<Double> flux_in, Matrix<Double> eps_mu, UInt idx);
    
    void EvaluateRotationState(UInt idx);
    
    Vector<Double> getRotationDirectionAndUpdate(Vector<Double> flux_in, Matrix<Double> eps_mu, UInt idx){
      
      UpdateRotationState(flux_in, eps_mu, idx);
      
      EvaluateRotationState(idx);
      
      return currentDirection_[idx];
    }
    
    Vector<Double> getRotationDirection(UInt idx){
      
      return currentDirection_[idx];
    }
    
  private:
    Vector<Double>* dirVector_;
    Double rotResistance_;
    UInt dim_;
    
    // big advantage compared to scalar model: setting values are always > 0
    // and we do not have to distinguish between minima and maxima; the rotation
    // directions get inserted only according to the amplitude / threshold as in
    // Vector model of the corresponding input flux quantity
    // > map can be used
    std::map<Double,Vector<Double> >* rotationStates_; 
    Vector<Double>* currentDirection_;
    
  };
  
  class VectorPreisachMayergoyz : public Hysteresis
  {
  public:
    // general constructor for anisotropic case, i.e. in each spatial direction
    // we have to specify xSaat,ySat as well as a set of fitting preisachWeights
    // the number of directions is specified by the length of xSat, ySat, etc
    VectorPreisachMayergoyz(Integer numElem, Vector<Double> xSat, Vector<Double> ySat, 
          Matrix<Double>* preisachWeight, UInt dim, bool isVirgin,Vector<Double> anhyst_A, Vector<Double> anhyst_B, Vector<Double> anhyst_C, bool anhystOnly);
    
    // constructor for isotropic case, i.e. same xSat, ySat and weights in all directions
    VectorPreisachMayergoyz(Integer numElem, UInt numDirections, Double xSat, Double ySat, 
          Matrix<Double>& preisachWeight, UInt dim, bool isVirgin,Double anhyst_A, Double anhyst_B, Double anhyst_C, bool anhystOnly);

    virtual ~VectorPreisachMayergoyz();

    Vector<Double> computeInput_vec(Vector<Double> Yin, Integer idx, Matrix<Double> eps_mu, bool overwriteDirection = true);
    
    Vector<Double> computeValue_vec(Vector<Double>& xVal, Integer idx, bool overwrite = true,bool overwriteDirection = true,bool debugOutput = false);
    
    void setFlags(UInt performanceFlag){
      ;
    };
    
  private:
    UInt dim_;
    UInt numDirections_;
    
    bool clipOutput_;
    bool isIsotropic_;
    
    Vector<Double>* singleDirections_;
    Preisach** singlePreisachOperators_;
    Matrix<Double> matrixForCoefComputation_;
    
  };
  
} //end of namespace


#endif

