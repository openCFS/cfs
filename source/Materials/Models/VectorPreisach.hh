#ifndef FILE_VECPREISACH_2016
#define FILE_VECPREISACH_2016

#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"
#include "Hysteresis.hh"
//#include "Preisach.hh"

/*! \brief Class for Vector Preisach Model
 *
 * Basic concept taken from:
 *    Sutor et al. "An efficient vector Preisach hysteresis model based on a novel rotational operator"
 *    doi: http://dx.doi.org/10.1063/1.3672069
 *
 * Similarities to std. scalar Preisach:
 *  1. Array/Matrix which stores weights for normalized input values alpha/beta
 *  2. Array/Matrix which stores switching operators for normalized input values alpha/beta
 *
 * Extensions by vector model:
 *  1. Array/Matrix storing for each alpha/beta a vector denoting the elements/domains orientation
 *  2. Switching operator at position alpha, beta is changed, if xPar is larger alpha or smaller beta;
 *      xPar results from the scalar product of the input vector xIn with the direction of the element;
 *      xPar is therewith a function of alpha, beta, too! (not explicitly named in aboves publication;
 *       note that Fig 1. in the publication does not visualize this correctly!)
 *
 * Remarks:
 *  1. vector model can only resolve input variation up to steps of size delta_ = 2/M, with M being the number
 *     of rows of the weight matrix; smaller input variaton cannot be captured, which leads to problems if input
 *     amplitudes get small; in the std. scalar Preisach model, this problem is solved during the evaluation via
 *     the Everett function, as this function just integrates over triangular areas between two succeeding input extrema;
 *     the triangles are thereby not bounded to the elements of the Preisach plane, but can be subareas. In that way,
 *     smaller input variations can be captured;
 *     Problem: Everett function not possible for vector model as the rotation
 *     operators may change each step and thus require a recalculation over the complete Preisach plane
 *  2. one concept to relief the issue above is to allow the switching operators to take values between +1 and -1;
 *     the value depends thereby on the "filling state" of the element;
 *     Problem: the history of the filling state is not tracked, i.e. the single elements are free of hysteresis
 *  3. to resolve the issue of 2., the idea was to store a list of input extrema for one specific element (called the
 *     staircase element in the following); this staircase element was determined by the last global input maxima and
 *     minima; if input maxima smaller than the global maximum or minima larger than the global minimum occured, we got
 *     a staircasing shape inside this element; all other elements were of course also affected by the input min and max
 *     but they were either split along a horizontal line (input max) or along a vertical line (input min) or
 *     refilled completely
 *
 *     Example:
 *
 *         p = positive value
 *         n = negative value
 *         _ / | = element boundaries
 *
 *         1. global extremum = global min
 *       __________ _v________
 *      |     n    |n¦   n    |
 *      |----------|-¦--------| < 2. global extremum = global max
 *      |     p    |p¦   p    |
 *      |          | ¦        |
 *      |__________|_¦________|
 *      |          | ¦        |
 *      |     p    |p¦   p    |
 *      |          | ¦        |
 *      |          | ¦        |
 *      |_________ |_¦________|
 *
 *         1. local extremum = local min
 *       __________ _______v__
 *      |     n    |n¦   n !n |
 *      |----------|-¦-----!--|
 *      |     p    |p¦   p !n | < staircased element
 *      |          | ¦     !  |
 *      |__________|_¦_____!__|
 *      |          | ¦     !  |
 *      |     p    |p¦   p !n | < vertical splitted element
 *      |          | ¦     !  |
 *      |          | ¦     !  |
 *      |_________ |_¦_____!__|
 *
 *         2. local extremum = local max
 *       __________ __________
 *      |     n    |n¦   n !n |
 *      |----------|-¦-----!--|
 *      |     p    |p¦   p !n | < staircased element
 *      |~~~~~~~~~~|~¦~~~~~!~ |
 *      |_____p____|p¦___p_!p_|
 *      |          | ¦     !  |
 *      |     p    |p¦   p !p | < element got completely refilled with p
 *      |          | ¦     !  |
 *      |          | ¦     !  |
 *      |_________ |_¦_____!__|
 *
 *
 *  3 (cont.). if global extrema changed, the old staircased element got filled up completely and a new element became the
 *     staircased element
 *
 *     Problem: as the maxima and minima are in terms of xPar, they are different for each element;
 *              -> there are no global minima or maxima!
 *              -> each element can have a staircased subdivision!
 *
 *  4. NEW Concept:
 *      Store a list of min/max for EACH element (jep this is going to be costly and probably very inefficient);
 *      With these lists we can reconstruct the staircasing in each element if necessary.
 *        element = -1: entry1 = beta; min
 *                      entry2 = xPar, but only if xPar is a maximum (or a minimum with equal value -> obsolet)
 *        element = +1: entry1 = alpha+delta_; max
 *                      entry2 = xPar, but only if xPar is a minimum (or a maximum with equal value -> obsolet)
 *        element = 0 (for elements on the diagonal alpha = -beta):
 *          a) first entrance of element by minimum defined by value xPar = betaMin:
 *              entry1 = alpha+delta/2; max
 *              entry2 = betaMin; min
 *          b) first entrance of element by maximum defined by value xPar = alphaMax:
 *              entry1 = beta+delta/2; min
 *              entry2 = alphaMax; max
 *
 *        element = 0 (for elements on the diagonal alpha = beta)
 *        -> can only be 1 element, the alpha = beta and alpha = -beta one -> alpha = beta = 0)
 *        -> exists only for the preisach plane having an odd number of rows!
 *        -> here we have to split a triangle into a positive triangle, a negative triangle and a rect of pos or neg value
 *        -> the sum of these areas shall be 0
 *          a) first entrance of element by minimum defined by value xPar = betaMin:
 *              entry1 = alpha + delta * sqrt(2)/2; max -> tested in formula below, results in area = 0
 *              entry2 = betaMin; min
 *          b) first entrance by maximum defined by xPar = alphaMax:
 *              entry1 = beta + delta * (1 - sqrt(2)/2); min -> tested in formula below, results in area = 0
 *              entry2 = alphaMax; max
 *
 *      Calculation of value for each element:
 *
 *      1. full (i.e. quadratic) element:
 *      value = 0.0
 *      cnt = 0
 *      old_entry = 0
 *      for each entry in list:
 *        if cnt == 0 and entry is min:
 *          --> element is split vertically from beginning
 *          value = value + ((entry - beta) - (beta+delta_ - entry))/delta_
 *        else if cnt == 0 and entry is max:
 *          --> element is split horizontally from beginning
 *          value = value + ((entry - alpha) - (alpha+delta_ - entry))/delta_
 *        else if cnt != 0 and entry is min:
 *          --> subtract overlap of negative area with former positive area; subtract * 2 to switch from + to -
 *          value = value - 2 * ((old_entry - alpha) * (beta + delta_ - entry))/delta^2
 *        else if cnt != 0 and entry is max:
 *          --> add overlap of positive area with former negative area; add * 2 to switch form - to +
 *          value = value + 2 * ((entry - alpha) * (beta + delta_ - old_entry))/delta^ 2
 *        else
 *          --> should not occur
 *
 *        old_entry = entry
 *        cnt++
 *
 *      2. half (i.e. triangular) element -> for alpha = beta
 *      value = 0.0
 *      cnt = 0
 *      old_entry = 0
 *      for each entry in list:
 *        if cnt == 0 and entry is min:
 *          --> element is split vertically from beginnig, but only one half of the element (alpha >= beta) needed!
 *          --> still right half is neg, left half is pos
 *          value = value - 0.5 * ((beta+delta_ - entry)/delta_) * ((alpha+delta_ - entry)/delta_) // upper triangle
 *          value = value + ((entry - beta)/delta_) * ((alpha+delta_ - entry)/delta_)              // upper left rectangle
 *          value = value + 0.5 * ((entry - beta)/delta_) * ((entry - alpha)/delta_)               // lower triangle
 *        else if cnt == 0 and entry is max:
 *          --> element is split horizontally from beginning; upper part = -1, lower part = +1
 *          value = value - 0.5 * ((beta+delta_ - entry)/delta_) * ((alpha+delta_ - entry)/delta_) // upper triangle
 *          value = value - ((entry - beta)/delta_) * ((alpha+delta_ - entry)/delta_)              // upper left rectangle
 *          value = value + 0.5 * ((entry - beta)/delta_) * ((entry - alpha)/delta_)               // lower triangle
 *        else if cnt != 0 and entry is min:
 *          --> subtract overlap of negative area with former positive area; this time the overlap is a triangle
 *          --> note the difference to the case of square elements; there we had rectangular (not necessary quadratic) updates
 *          --> here we have trianular updates, where the triangular is a (quadratic) square separated at its diagonal
 *          value = value - 2 * 0.5*((old_entry - entry)/delta_)^2
 *        else if cnt != 0 and entry is max:
 *          --> add overlap of positive area with former negative area; again this area is triangular
 *          value = value + 2 * 0.5*((entry - old_entry)/delta_)^2
 *        else
 *          --> should not occur
 *
 * 5. EXTENSION:
 *      Due to the new concept, we are able to change the switching operators independently of the discretization delta
 *      of the Preisach plane. However, the rotation operators remained discrete. In case of small input signals this may
 *      lead to problems again as the elements on the line alpha = -beta are initialized with a rotation direction of (0,0,0).
 *      Only if the input becomes large enough to overcome the threshold of the elements on this line, these elements will
 *      contribute to the result. This switch will always occur as a jump in the ouput signals as suddenly a whole new
 *      element is added to the output sum.
 *
 *      Suggestion:
 *        Instead of rapidly switching the rotation operators from its old value to the new value, we make some sort of
 *        linear interpolation, thus allowing "microrotations":
 *
 *      Original version:
 *        rot =
 *        a  e_u       if uthres > alpha or uthres < beta
 *        b  rot_old   else
 *
 *                  alpha
 *           ________|__________
 *          |bbbb¦aaa|aaaaaaaa/
 *          |bbbb¦aaa|aaaaaa/
 *          |----+---|----/ uthres
 *          |aaaa¦aaa|aa/
 *          |aaaa¦aaa|/_________ beta
 *          |aaaa¦aa/
 *          |aaaa¦/
 *          |aaa/ -uthres
 *          |a/
 *          /
 *
 *      Suggested version:
 *      a if uthres > alpha+delta -> e_u
 *        else
 *      b   if -ut < beta -> e_u
 *          else
 *      c     if ut > alpha and ut < alpha+delta -> lin interpolation
 *            else
 *      d       if -ut > beta and -ut < beta+delta -> lin interpolation
 *      e       else rot_old
 *
 *
 *                  alpha
 *           ________|__________
 *          |eeee¦d!b|bbbbbbbb/     # -> e_u
 *          |eeee¦d!b|bbbbbb/
 *          |----+---|----/ uthres
 *          |cccc¦c!b|bb/
 *          |====¦=!b|/_________ beta
 *          |aaaa¦aa/
 *          |aaaa¦/
 *          |aaa/ -uthres
 *          |a/
 *          /
 *
 */
/*! \brief Conventions
 *
 *  1. All elements inside the Preisach plane are addressable by a pair of doubles (alpha, beta) with
 *        alpha, beta \in [-1,1-delta_].
 *  2. The pair (alpha, beta) is thereby the coordinate of the lower left corner of the element.
 *  3. The left and lower boundary are part of the element, right and upper boundary belong to
 *        the neighboring element.
 *
 */

namespace CoupledField {

  class ListEntryOld
  {
  public:

    ListEntryOld(Double value, Vector<Double> vector, bool isMin, bool isDummy = false){
      val_ = value;
      vec_ = vector;
      isMin_ = isMin;
      isDummy_ = isDummy;
    }
    ~ListEntryOld(){
    }

    bool isMin(){
      return isMin_;
    }

    bool isDummy(){
      return isDummy_;
    }

    void toggleMin(){
      isMin_ = !isMin_;
    }

    void setVal(Double newValue){
      val_ = newValue;
    }

    void setVec(Vector<Double> newVector){
      vec_ = newVector;
    }

    Double getVal(){
      return val_;
    }

    Vector<Double> getVec(){
      return vec_;
    }

    std::string ToString() const {

      std::ostringstream os;

      os << "Value: " << val_ << '\n';
      os << "Vector: " << vec_.ToString() << '\n';
      os << "IsMin?: " << bool(isMin_) << '\n';

      return os.str();
    }

  private:

    Double val_;
    Vector<Double> vec_;
    bool isMin_;
    bool isDummy_;

  };


  class VectorPreisach : public Hysteresis
  {
  public:
    //! constructor
    VectorPreisach(Integer numElem, Double xSat, Double ysat,
	     Matrix<Double>& preisachWeight, Double rotationalResistance , UInt dim, bool isVirgin,
	     bool isTesting = false, UInt evalVersion = 4);

    //! destructor
    virtual ~VectorPreisach();

    //! this function gets called from outside and calculates the output of the Preisach operator
    Vector<Double> computeValue_vec(Vector<Double>& xVal, Integer idElem);

    //! returns the current output of the hyst-operator for element idxElem
    Vector<Double> getValue_vec(Integer idElem);

    void switchingStateToBmp(UInt numPixel, std::string filename, UInt idElem, bool overLayWithRotState = false);

  private:

    //! initialize lists
    void Initialize_MinMaxList(std::list<ListEntryOld>& list, Double alpha, Double beta, UInt idElem, UInt idAlpha, UInt idBeta);

    void Initialize_RotList(std::list<ListEntryOld>& list, Double alpha, Double beta, UInt idElem, UInt idAlpha, UInt idBeta);

    // calculate overlapping area of two rectangles; returns true if overlap exists
    bool clipRectangles(Double t1, Double b1, Double l1, Double r1, Double t2, Double b2, Double l2, Double r2, Double& tRet, Double& bRet, Double& lRet, Double& rRet);

    Double getRectangleBounds(std::list<ListEntryOld>& list,std::list<ListEntryOld>::iterator startIt, std::list<ListEntryOld>::iterator endIt,
                            Double alpha, Double beta, UInt idElem, UInt idAlpha, UInt idBeta, UInt idArea,
                            bool isTriangle, bool isRot,
                            Double& tRet, Double& bRet, Double& lRet, Double& rRet);

    //! function used to update
    //!     1. list of rotation operators for the evaluation of the current orientation of each element
    //!     2. list of minima and maxima of xPar which gets used to calculate the value of the switching operator
    //!
    //! input
    //!     case 1: x = xthres, xVec = e_u, (alpha,beta) = coordinates of the bottom left corner of the current Preisach element
    //!             list = reference to preRotVecs_[idxElem][idAlpha][idBeta]
    //!     case 2: x = xPar (uncut), xVec = /, (alpha,beta) as above, list = reference to preMinMaxVals_[idxElem][idAlpha][idBeta]
    //!
    //! (idElem,idAlpha and idBeta are determined in the calling function)
    void UpdateList(std::list<ListEntryOld>& list, Double newEntry, Vector<Double> newVecEntry, Double alpha, Double beta, UInt idElem, UInt idAlpha, UInt idBeta, bool isRot);

    //! Write out list
    void outputList(std::list<ListEntryOld>& list,bool isRot);

    //! function used to calculate
    //!     1. averaged orientation of element with coordinates (alpha,beta) from updated list preRotVecs_[idxElem][idAlpha][idBeta]
    //!     2. averaged switching value of element with coords. (alpha,beta) from updated list preMinMaxVals_[idxElem][idAlpha][idBeta]
    void EvaluatePreisachElement(std::list<ListEntryOld>& list, Double& retVal, Vector<Double>& retVec,
                                 Double alpha, Double beta, UInt idElem, UInt idAlpha, UInt idBeta, bool isRot);
    void EvaluatePreisachElement_v2(std::list<ListEntryOld>& list, Double& retVal, Vector<Double>& retVec,
                                    Double alpha, Double beta, UInt idElem, UInt idAlpha, UInt idBeta, bool isRot);
    void EvaluatePreisachElement_v3(std::list<ListEntryOld>& listl, Double& retVal, Vector<Double>& retVec,
                                    Double alpha, Double beta, UInt idElem, UInt idAlpha, UInt idBeta, bool isRot);
    //! No retVal needed here, as retVec is already the fully computed element vector; isRot is not needed either as we treat both lists at the same time
    void EvaluatePreisachElement_v4(std::list<ListEntryOld>& rotList, std::list<ListEntryOld>& minMaxList,Double& retVal,Vector<Double>& retVec,
                                    Double alpha, Double beta,UInt idElem, UInt idAlpha, UInt idBeta);

    /*!
     * Global quantities, i.e. the same for all FE elements of the same material
     */
    Double Xsaturated_; //! saturation value for  input
    Double YSaturated_; //! saturation value for output

    Matrix<Double> preisachWeights_; //! preisach weight function
    Double rotationalResistance_; //! parameter describing the resistance of the domains to rotation

    bool isVirgin_; //! yes, if starting at zero; currently flag is unused

    UInt dim_; //! 2D or 3D (axi not tested at all!)
    UInt numElem_; //! total number of FE elements
    UInt numRows_; //! number of rows of the Preisach plane

    Double delta_; //! resolution of Preisach plane
    Double tol_; //! tolerance for all kind of comparisons

    /*!
     * Debugging quantities
     */
    Double evalVersion_;
    bool isTesting_; //! if true, rotationalWeigts will be initialized in +y direction and their update is deactivated
    Vector<Double> rotStateTesting_;

    /*!
     * Local quantities, i.e. arrays storing different values for each FE element
     */
    Vector<Double>* preisachSum_; //! output value of Preisach operator

    Matrix<Double>* oldXpar_; //! stores for each element of the Preisach plane the value of xPar from the last calculation step
    Matrix<Double>* switchingStates_; //! stores the current switching state
    Vector<Double>*** rotationStates_; //! stores the current rotation state
    Vector<Double>*** evaluatedState_; //! stores the completely evaluated element state here (for version 4)

    bool** wipedOut_; //! this array is for the elements on line alpha = -beta (center coordinates here)
    //! these elements need a special treatment as they are split along the diagonal alpha = -beta in +1/-1 which
    //! makes not only the calculation of the switching value but also the storage of the mins and maxs more complicated
    //! however, once the element on the diagonal has been set to +1 or -1 by a strong enough input, we can (and must) treat it as
    //! a std. element -> wipedOut_ stores exactly this information
    /*
     * NOTE: currently diagonal splitted elements not working perfectly for odd number of rows of the Preisach plane.
     *    reason: element with center 0,0 lies on alpha = beta and alpha = -beta then and such is not only split
     *    in half but in 4 triangular sections which need an even more sophisticated treatment
     *    -> give out a warning; error is hopefully not too large
     */

    //! for each FE element we need for each Preisach element a list of min and max values
    //! Important note: it seems not possible to use our own vector class in combination with std::vector
    //! -> Matrices will have no entries, segfaults, etc
    std::list<ListEntryOld>*** preMinMaxValues_;
    std::list<ListEntryOld>*** preRotValues_;

    //! flag checking if the min/max list has changed; if not, the switching value will also be the same so we do not
    //! need to recompute it
    bool*** preMinMaxChanged_;
    bool*** preRotChanged_;

    //! for optimized version
    //! according to the original vector model, the lower triangle part (alpha <= 0) always has a
    //! switching state of +1
    //! reason: the rotation state in the lower triangle will be set to e_u in each iteration, as xThres >= 0
    //! the scalar product of u_in with its direction e_u will give the length of u_in which has to be >= 0
    //! -> the lower triangluar area has always a switching state of +1 in each element, such that we just once
    //! have to sum up all weights in this region; this value is then multiplied with the current e_u vector in each
    //! iteration
    Double lowerTriangleValue_;


  };
} //end of namespace


#endif

