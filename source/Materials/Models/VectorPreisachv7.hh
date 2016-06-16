#ifndef FILE_VECPREISACHV7_2016
#define FILE_VECPREISACHV7_2016

#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"
#include "Hysteresis.hh"
//#include "Preisach.hh"

namespace CoupledField {

class ListEntry
  {
  /*
   * New version of ListEntry
   * as ListEntry is only used for the switching list, it does not have to store vectors anymore
   */
  public:

    ListEntry(Double value, bool isMin, bool isDummy = false){
      val_ = value;
      isMin_ = isMin;
      isDummy_ = isDummy;
      hasChanged_ = true;
    }
    ~ListEntry(){
    }

    bool isMin(){
      return isMin_;
    }

    bool isMinConst() const{
      return isMin_;
    }

    bool isDummy(){
      return isDummy_;
    }

    void toggleMin(){
      isMin_ = !isMin_;
    }

    void setVal(Double newValue){
      if(newValue != val_){
        hasChanged_ = true;
      }
      val_ = newValue;
    }

    Double getVal(){
      return val_;
    }

    Double getValConst() const{
      return val_;
    }

    bool operator==(ListEntry x) const{
      Double tol = 1e-16;

      if(abs(val_ - x.getVal()) > tol){
        //std::cout << "Diffval: " << abs(val_ - x.getVal()) << std::endl;
        return false;
      }
//      if(val_ != x.getVal()){
//        return false;
//      }
      if(isMin_ != x.isMin()){
        return false;
      }
      return true;
    }

    std::string ToString() const {

      std::ostringstream os;

      os << "Value: " << val_ << '\n';
      os << "IsMin?: " << bool(isMin_) << '\n';

      return os.str();
    }

    void setHasChanged(bool state){
      hasChanged_ = state;
    }

    bool hasChanged(){
      return hasChanged_;
    }

  protected:

    Double val_;
    bool isMin_;
    bool isDummy_;

    /*
     * marks if this entry has changed;
     * hasChanged_ becomes true, if the switching list gets updated;
     */
    bool hasChanged_;

  };

  class RotListEntry : public ListEntry
  {
  /*
   * Extension to normal list entries
   *
   * Each entry of the rotation list has to store a list of std ListEntries which are used to reconstruct the
   * switching state
   */
  public:

      /*
       * upperValue corresponds to the value of xThres which is used to set the entry
       * lowerValue is the value of the next entry (and thus needs to be updated)
       * -> using these two values, we can directly decide, if the area of the rotation state has changed or not
       * -> if the area did not change and the switching list did not change either, we do not have to recalculate
       *    the rotation state but only have to multiply with the rotation vector
       */
      RotListEntry(Double upperValue, Double lowerValue, Vector<Double> vector, std::list<ListEntry> switchingList, Double lastLocalXpar= 0.0,
                   bool isMin = false, bool isDummy = false, bool wasWipedOut=false, UInt startCnt = 0)
      : ListEntry(upperValue, isMin, isDummy)
      {
        // create empty list
        switchingList_ = switchingList;
        vec_ = vector;
        lowerVal_ = lowerValue;
        evaluatedVal_ = 0.0;
        rotHasChanged_ = true;
        lastLocalXpar_ = lastLocalXpar;
        wipedOut_ = wasWipedOut;
        startCnt_ = startCnt;
      }

      ~RotListEntry(){
      }

      /*
       * To keep the later rotation list as short as possible, we try to merge adjacent rotListEntries.
       * This is possible, if
       *  - they share the same rotation vector
       *  - they have the same switching list (e.g. after two adjacent regions were
       *     completely overwritten by a new input, which not only set the rotation vector
       *     but also set resulted in an xPar which overwrote the switching list)
       */
      bool canBeMergedWith(RotListEntry& rotEntry){

        Double tol = 1e-15;
        Vector<Double> tmp = rotEntry.getVecCopy();

        for(UInt i = 0; i< vec_.GetSize();i++){
          if(abs(vec_[i]-tmp[i])>tol){
            //std::cout << "Diffvec " << abs(vec_[i]-tmp[i]) << std::endl;
            return false;
          }
        }

        if(startCnt_ != rotEntry.getStartCnt()){
          /*
           * even though the switching lists may be the same,
           * one rotEntry starts evaluation at 0 the other one at >= 1.
           * But >= 1 means that the real first entry was deleted from the list
           * and thus the list wouldn't be equal
           */
          return false;
        }

//        if(!(vec_ == rotEntry.getVecCopy())){
//          return false;
//        }

        if(switchingList_ != rotEntry.getListReference()){
          return false;
        }

        return true;
      }

      std::list<ListEntry>& getListReference(){
        return switchingList_;
      }

      std::list<ListEntry> getListCopy(){
        return switchingList_;
      }

      void setList(std::list<ListEntry> switchingList){
        switchingList_ = switchingList;
      }

      std::string ToString() const {

        std::ostringstream os;

        os << "RotationEntry: " << '\n';
        os << "Upper/Lower Value: " << val_ << " / " << lowerVal_ << '\n';
        if(vec_.GetSize() == 3){
          os << "Vector: " << vec_[0] << " / " << vec_[1] << " / " << vec_[2] << '\n';
        } else {
          os << "Vector: " << vec_[0] << " / " << vec_[1] << '\n';
        }
        //os << "IsMin?: " << bool(isMin_) << '\n'; // is always 0
        os << "List of SwitchingEntries: " << '\n';
        os << "Value (isMin?): " << '\n';

        std::list<ListEntry>::const_iterator listIt;
        UInt cnt = 1;

        for(listIt = switchingList_.begin(); listIt != switchingList_.end(); listIt++){
          os << listIt->getValConst() << " (" << listIt->isMinConst() << ") -/- ";
          cnt++;
        }

        return os.str();
      }

      void setLowerVal(Double newValue){
        if(lowerVal_ != newValue){
          hasChanged_ = true;
        }
        lowerVal_ = newValue;
      }

      void setVec(Vector<Double> newVector){

        Double tol = 1e-15;
        for(UInt i = 0; i< vec_.GetSize();i++){
          if(abs(vec_[i]-newVector[i])>tol){
            rotHasChanged_ = true;
            break;
          }
        }

        vec_ = newVector;
      }

      Vector<Double>& getVecReference(){
        return vec_;
      }

      Vector<Double> getVecCopy(){
        return vec_;
      }

      Double getLowerVal(){
        return lowerVal_;
      }

      Double getLastLocalXpar(){
        return lastLocalXpar_;
      }

      void setLastLocalXpar(Double newLastLocalXpar){
        if(newLastLocalXpar != lastLocalXpar_){
          hasChanged_ = true;
        }
        lastLocalXpar_ = newLastLocalXpar;
      }

      void setToUnchanged(){
        /*
         * sets the value hasChanged_ of the element and of
         * each entry of the switching list to false
         * -> to be done after evaluating the list so that this state
         *    can be reused if neither list nor rotlistelem changes
         *    till next time
         */
        hasChanged_ = false;

        std::list<ListEntry>::iterator listIt;

         /*
          * Note that we can detect overall changes to the switchingList by checking its content!
          * Reason: the list changes if elements get
          *  - values changes -> setVal sets flag hasChanged_ to true
          *  - entries are deleted/overwritten by new ones -> new entries have automatically hasChanged_ = true
          *  - entries get appended to the list -> new entries have automatically hasChanged_ = true
          *
          *  -> the loop can only be passed, if the number of elements and their content remained the same
          */
         for(listIt = switchingList_.begin(); listIt != switchingList_.end(); listIt++){
           listIt->setHasChanged(false);
         }
      }

      bool hasChanged(){
        /*
         * check if either the rotation state entry has
         * changed its lower bound or if the switching list
         * has changed
         * NOTE: do not check for rotation state changes
         * -> this function should only help to determine if
         *    the overlap of switching areas and rotation area
         *    has changed
         * -> if not, we do not need to recalculate the overlap!
         */
        if(hasChanged_ == true){
          return true;
        }

        std::list<ListEntry>::iterator listIt;

        /*
         * Note that we can detect overall changes to the switchingList by checking its content!
         * Reason: the list changes if elements get
         *  - values changes -> setVal sets flag hasChanged_ to true
         *  - entries are deleted/overwritten by new ones -> new entries have automatically hasChanged_ = true
         *  - entries get appended to the list -> new entries have automatically hasChanged_ = true
         *
         *  -> the loop can only be passed, if the number of elements and their content remained the same
         */
        for(listIt = switchingList_.begin(); listIt != switchingList_.end(); listIt++){
          if(listIt->hasChanged() == true){
            return true;
          }
        }

        return false;
      }

      bool wasListWipedOut(){
        return wipedOut_;
      }

      void setWipedOut(bool state){
        if(state != wipedOut_){
          hasChanged_ = true;
        }
        wipedOut_ = state;
      }

      Double getLastEvalState(){
        return evaluatedVal_;
      }

      void setLastEvalState(Double newState){
        evaluatedVal_ = newState;
      }

      UInt getStartCnt(){
        return startCnt_;
      }

      void setStartCnt(UInt newVal){
        /*
         * Note: we do nt set flag hasChanged to true here
         *
         * Reason: hasChanged is used to determine if the overlapping state
         * might have changed, i.e. due to a new lower value or new entries in the switching list.
         * This startCnt_ is set only if unneeded elements were removed from the list, i.e. the
         * resulting switching state does (hopefully) not change and thus we do not have
         * to reevaluate the switching state.
         *
         * However, this value has to be checked during merging process!
         */
        startCnt_ = newVal;
      }

  private:
      std::list<ListEntry> switchingList_;

      /*
       * When evaluating the switching list we must iterate over the list and call the function
       * getRectangleBounds with a cnt specifying the area id (as there is a special treatment for
       * area with id = 0). However, to simplify and shorten the list, it can and will happen, that the
       * entries corresponding to area with id=0 get deleted. Therefore the cnt in the loop should not
       * start at 0 as this would lead to wrong results.
       * Therefore, add a start counter which is either 0 (meaning the value corresponding to area 0 is
       * still in the list) or 1 (meaning the value is no longer in the list).
       * This value should be reset to 0 each time the switching list got completely overwritten, i.e.
       * the first element of the list got replaced.
       */
      UInt startCnt_;

      /*
       * marks if this the rotation state has changed;
       */
      bool rotHasChanged_;

      /*
       * lower value defining the area
       * upper value is stored in val_
       */
      Double lowerVal_;

      /*
       * stores the evaluated state resulting from
       * overlapping switching list with rotation area and Preisach weights
       * the overall evaluated state results from this value * rotationState (= vec_)
       * -> this value needs only be recalculated if hasChanged_ == true
       */
      Double evaluatedVal_;

      /*
       * rotation state corresponding to this entry
       */
      Vector<Double> vec_;

      /*
       * stores the last value of xPar which was used to set the locally stored switching list
       */
      Double lastLocalXpar_;

      /*
       * this flag states, if the stored switching list reached the absolute min (beta = -1) or absolute max (alpha = +1)
       * if that is the case, the resulting switching areas do not have a subarea around alpha = -beta which is still split
       * along this diagonal;
       */
      bool wipedOut_;
  };


  class VectorPreisachv7 : public Hysteresis
  {
  /*!
   *
   * Completely reworked approach (v7)
   *
   * Concept:
   *  1. Divide Preisach plane into three regions
   *     - - - - - -
   *    | 1   |2  /
   *    |     | /
   *     - - -
   *    | 3 /
   *    | /
   *
   *  2. Treatment of region 1 (-1 <= beta <= 0; 0 <= alpha <= 1)
   *    Use a global list of RotListEntry to define the rotation state of this region. Each entry of this list has
   *    a pair xThres, e_u for the calculation of the rotation state. Additionally, each entry has a list of ListEntry.
   *    This second list is used to compute the switching state fitting to the rotation state. Using this stacked approach,
   *    the calculation of xPar can be done without the need to average the rotation state.
   *    The evulation of region 1 is done in the following steps:
   *
   *    1) update global rotation list with current input
   *    2) initialize evaluated state with 0
   *    3) for each entry in the global list
   *       i) get rotation vector and calculate xPar
   *       ii) update switching list with xPar
   *       iii) get x and y values defining the area which corresponds to the rotation state (see sketch below)
   *       iv) for each entry in switching list
   *          a) get x and y values defining the area corresponding to the current switching entry
   *          b) calculate overlap of switching area and rotation area
   *          c) map the overlap-area onto the Preisach plane; sum up all weights in the mapped part (cut if
   *              matrix entries are only partially overlapped)
   *          d) apply right sign depending on switching state being a min or a max
   *          e) multiply resulting summed up value with rotation state and add to evaluated state
   *
   *
   *       Possible rotation states:
   *         __ __ __ __ __ __ __ __ __ __
   *        |        |     |        |     |
   *        |  0     | 1   |  2     | 3   |
   *        |________|     |        |     |
   *        |              |        |     |
   *        |______________|        |     |
   *        |                       |     |
   *        |                       |     |
   *        |_______________________|     |
   *        |                             |
   *        |__ __ __ __ __ __ __ __ __ __|
   *
   *          -> areas are squares (0) or L-shaped (1,2,3)
   *
   *       Each rotation state has an own switching state which can have the form of a stair-case.
   *       Possible rotation state for area 1:
   *         __ __ __ __ __ __ __ __ __ __
   *        |        |                    |
   *        |  +1    |__            -1    |
   *        |           |__               |
   *        |              |              |
   *        |              |________      |
   *        |                       |     |
   *        |                       |     |
   *        |                       |     |
   *        |                       |     |
   *        |__ __ __ __ __ __ __ __|__ __|
   *
   *          -> determine the overlap switching area and rotation area to find correct value
   *
   *       Overlap of rotation state 1 and its switching state:
   *         __ __ __ __ __ __ __ __ __ __
   *        |        |xxxxx|              |
   *        |        |xxxxx|              |
   *        |________|..xxx|              |
   *        |..............|              |
   *        |..............|              |
   *        |                             |
   *        |                             |
   *        |                             |
   *        |                             |
   *        |__ __ __ __ __ __ __ __ __ __|
   *
   *          -> only dotted and "x"ed area count to final value
   *          -> before summing them together, map these two areas separately to the Preisach weights; then sum up
   *              the Preisach weights in the map area
   *
   *
   *  3. Treatment of region 2 (0 <= beta <= alpha; 0 <= alpha <= 1)
   *    According to the original model, this upper triangular region is set to the current value of e_u in each step
   *    as xThres >= 0 >= -beta for all inputs. Therewith, all parts of this region have the same rotation state at each given
   *    evaluation time and, thus, the value of xPar is the same for each point in that region.
   *    By this, the switching state of this region can be determined by a global list of non-wiped out minima and maxima
   *    like in the standard, scalar Preisach model.
   *    The evaluated state of this region is determined by evaluating the switching state, overlapping it with the Preisach
   *    weights and multiplying the resulting value with the current rotation state e_u.
   *
   *  4. Treatment of region 3 (-1 <= beta <= alpha; -1 <= alpha <= 0)
   *    From the original model we can directly derive, that the rotation state of this area will be set to e_u in each
   *    step as xThres >= 0 >= alpha for all inputs (compare to region 2). Additionally, the switching state is completely
   *    +1 for all times, as xPar = u_in dot e_u = |u_in| >= 0 >= alpha.
   *    To evaluate this region, we have to sum up all Preisach weights in the lower triangle (and maybe cut them at alpha = 0
   *    in half if the number of rows is odd) and multiply the value with e_u. As the switching state stays +1 for all times,
   *    the summation can be done during the initialization.l
   *
   */

  public:
    VectorPreisachv7(Integer numElem, Double xSat, Double ySat,
           Matrix<Double>& preisachWeight, Double rotationalResistance , UInt dim, bool isVirgin,
           bool isTesting = false, UInt evalVersion=7 );

    //virtual ~VectorPreisachv7();
    virtual ~VectorPreisachv7();

    //! this function gets called from outside and calculates the output of the Preisach operator
    Vector<Double> computeValue_vec(Vector<Double>& xVal, Integer idElem);

    //! returns the current output of the hyst-operator for element idxElem
    Vector<Double> getValue_vec( Integer idElem )
    {
      return ( preisachSum_[idElem] );
    }

    void switchingStateToBmp(UInt numPixel, std::string filename, UInt idElem, bool overLayWithRotState = false);
    void rotationStateToBmp(UInt numPixel, std::string filename, UInt idElem, UInt comp = 0, bool writeOutput = true, Matrix<Double>* helperMatrix = NULL);

  private:

    // calculate overlapping area of two rectangles; returns true if overlap exists
    bool clipRectangles(Double t1, Double b1, Double l1, Double r1, Double t2, Double b2, Double l2, Double r2, Double& tRet, Double& bRet, Double& lRet, Double& rRet);

    Double clipToElement(Double rectT, Double rectB, Double rectL, Double rectR, UInt idAlpha, UInt idBeta, Double delta = -1);

    Double mapRectangleToPreisachWeights(Double t, Double b, Double l, Double r, bool skipUpperDiagonal = false);

    //! only for output visualization; not needed during std calculation
    void mapRectangleToHelperMatrix(Matrix<Double>& helper, Double t, Double b, Double l, Double r, Double factor, bool skipUpperDiagonal = false);

    void Simplify_GlobalRotationList(UInt idElem);

    void Simplify_LocalSwitchingLists(UInt idElem);

    void Initialize_GlobalSwitchingList(UInt idElem);

    void Initialize_GlobalRotationList(UInt idElem);

    UInt Update_SwitchingList(std::list<ListEntry>& list, Double newEntry, Double lastXpar, bool isGlobal, bool wasWipedOut, Double lowerVal = 0);

    void Update_GlobalRotationList(Double xThres, Vector<Double> e_u, Vector<Double> u_in, UInt idElem);

    void Evaluate_GlobalSwitchingList(UInt idElem);

    void Evaluate_GlobalRotationList(UInt idElem, Vector<Double>& retVec);

    void Evaluate_LowerTriangle();

    Double getRectangleBounds(std::list<ListEntry>& list,std::list<ListEntry>::iterator startIt, std::list<ListEntry>::iterator endIt,
                              UInt idArea, bool isGlobal, Double& tRet, Double& bRet, Double& lRet, Double& rRet, bool upperSplitSquare = false);

    // NOT NEEDED YET
    // void outputList(std::list<ListEntry>& list,bool isRot);

    // NOT NEEDED
    //    bool GlobalSwitchingChanged(UInt idElem){
    //      /*
    //       * Iterate over global list of switching states and check each entry for a change;
    //       * the list needs to be reevaluated only, if at least one entry changed
    //       */
    //      std::list<ListEntry>::iterator listIt;
    //
    //      for(listIt = globSwitchList_[idElem].begin(); listIt != globSwitchList_[idElem].end(); listIt++){
    //        if(listIt->hasChanged() == true){
    //          return true;
    //        }
    //      }
    //
    //      return false;
    //    }

    /*!
     * Global quantities, i.e. the same for all FE elements of the same material
     */
    Double Xsaturated_; //! saturation value for  input
    Double YSaturated_; //! saturation value for output

    Matrix<Double> preisachWeights_; //! preisach weight function
    Double rotationalResistance_; //! parameter describing the resistance of the domains to rotation

    bool isVirgin_; //! yes, if starting at zero; currently flag is unused
    bool isSymmetric_; //! states if Preisach weights are symmetric w.r.t alpha = -beta

    UInt dim_; //! 2D or 3D (axi not tested at all!)
    UInt numElem_; //! total number of FE elements
    UInt numRows_; //! number of rows of the Preisach plane

    Double delta_; //! resolution of Preisach plane
    Double tol_; //! tolerance for all kind of comparisons

    Double lowerTriangleValue_;

    /*!
     * Local quantities, i.e. one for each FE elements of the same material
     */
    Vector<Double>* preisachSum_;

    std::list<RotListEntry>* globRotList_;
    std::list<ListEntry>* globSwitchList_;
    /*
     * lowerTriangleValue_ is the same for all FE elements using the same set of Preisach Weights
     * upperTriangleValue_ differs for each FE element
     */
    Vector<Double> upperTriangleValue_;
    /*
     * stores the last value of xPar which was used to update the globalSwitching list
     */
    Vector<Double> lastXpar_;

    /*
     * last rotation state; only used in output function rotationStateToBmp
     */
    Vector<Double>* lastEu_;

    /*!
     * Debugging quantities
     */
    Double evalVersion_;

    // CURRENTLY NOT IMPLEMENTED
    //bool isTesting_; //! if true, rotationalWeigts will be initialized in +y direction and their update is deactivated

  };
} // end of namespace

#endif

