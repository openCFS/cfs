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

  class VectorPreisach : public Hysteresis
  {
  public:
    //! constructor
    VectorPreisach(Integer numElem, Double xSat, Double ysat,
	     Matrix<Double>& preisachWeight, Double rotationalResistance , UInt dim, bool isVirgin);

    //! destructor
    virtual ~VectorPreisach();

    //! update rotational Weights
    void updateRotationalWeights(Vector<Double>& Xin, Integer idx);

    //! update switching Operators; has to be called after updateRotationalWeights!
    void updateSwitchingOperators(Vector<Double>& Xin, Integer idx);

    //! evaluates the min/max list to compute the value of Preisach element alpha,beta of FE element idx
    Double evaluatePreisachElementValue(Double alpha,Double beta,UInt idx);

    //! this function gets called from outside and calculates the output of the Preisach operator
    Vector<Double> computeValue_vec(Vector<Double>& xVal, Integer idxElem);

    //! returns the current output of the hyst-operator for element idxElem
    Vector<Double> getValue_vec(Integer idxElem);

    //! Updates the min and max maps for microswitching
    void updateMinMax(Double xPar, Double alphaMax, Double betaMin, UInt idxElem);

  private:

    /*!
     * Global quantities, i.e. the same for all FE elements of the same material
     */
    Double Xsaturated_; //! saturation value for  input
    Double YSaturated_; //! saturation value for output

    Matrix<Double> preisachWeights_; //! preisach weight function
    Double rotationalResistance_; //! parameter describing the resistance of the domains to rotation

    bool isVirgin_; //! yes, if starting at zero; currently flag is unused

    UInt dim_; //! 2D or 3D (axi not tested at all!)
    UInt numRows_; //! number of rows of the Preisach plane

    Double delta_; //! resolution of Preisach plane

    bool isTesting_; //! if true, rotationalWeigts will be initialized in +y direction and their update is deactivated

    /*!
     * Local quantities, i.e. arrays storing different values for each FE element
     */
    Vector<Double>* preisachSum_; //! output value of Preisach operator
    Matrix<Double>* rotationalWeightsX_; //! rotational operator, x-component
    Matrix<Double>* rotationalWeightsY_; //! rotational operator, y-component
    Matrix<Double>* rotationalWeightsZ_; //! rotational operator, z-component ( for 3D )
    Matrix<Double>* switchingStates_; //! current state of the switching operators
    Matrix<Double>* oldXpar_; //! stores for each element of the preisach plane the value of xPar from the last calculation step

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

    bool* rotationalWeightsUpdated_; //! flag denoting if rotationalWeights have already been updated
    bool* switchingStatesUpdated_; //! flag denoting if switchingstates have already been updated

    //! for each FE element we need to know for each Preisach element if it was initialized with a min or a max
    bool*** firstMin_;

    //! for each FE element we need for each Preisach element a list of min and max values
    //! Importante note: it seems not possible to use our own vector class in combination with std::vector
    //! -> Matrices will have no entries, segfaults, etc
    std::list<Double>*** preMinMaxValues_;

    //! flag checking if the min/max list has changed; if not, the switching value will also be the same so we do not
    //! need to recompute it
    bool*** preMinMaxChanged_;

  };


} //end of namespace


#endif

