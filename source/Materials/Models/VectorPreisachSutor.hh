#ifndef FILE_VECPREISACHV10_2016
#define FILE_VECPREISACHV10_2016

#include <list>
#include <iomanip>
#include "MatVec/Vector.hh"
#include "MatVec/Matrix.hh"
#include "Hysteresis.hh"
#include "Utils/Timer.hh"
//#include "Preisach.hh"

namespace CoupledField {

  class Rectangle
  {
  public:

    Rectangle(){
      l_ = 0.0;
      r_ = 0.0;
      t_ = 0.0;
      b_ = 0.0;
    }

    Rectangle(Double left, Double right, Double top, Double bot){
      l_ = left;
      r_ = right;
      t_ = top;
      b_ = bot;
    }

    void setBounds(Double left, Double right, Double top, Double bot){
      l_ = left;
      r_ = right;
      t_ = top;
      b_ = bot;
    }

    void getBounds(Double& left, Double& right, Double& top, Double& bot){
      left = l_;
      right = r_;
      top = t_;
      bot = b_;
    }

    std::string ToString() const {

      std::ostringstream os;

      os << std::setprecision(15) << "Left / Right: " << l_ << " / " << r_ << '\n';
      os << std::setprecision(15) << "Top / Bottom: " << t_ << " / " << b_ << '\n';

      return os.str();
    }


    bool clipRectangles(Rectangle& source, Rectangle& target){
      Double l,r,t,b;

      if( (this->t_ <= source.b_)||(source.t_ <= this->b_) ) {
        /*
         * top edge of one rectangle is smaller lies below or on top of the bottom edge of the other rectangle
         * -> no overlap area possible
         */
        return false;
      }
      if( (this->r_ <= source.l_)||(source.r_ <= this->l_) ){
        /*
         * right edge lies left or on top of the left edge of the other element
         * -> no overlap area possible
         */
        return false;
      }

      t = std::min(this->t_,source.t_);
      b = std::max(this->b_,source.b_);
      l = std::max(this->l_,source.l_);
      r = std::min(this->r_,source.r_);

      target = Rectangle(l,r,t,b);
      return true;
    }

    Double l_,r_,t_,b_;

  };


  class ListEntryv10
  {
    /*
     * New version of ListEntryv10
     * as ListEntryv10 is only used for the switching list, it does not have to store vectors anymore
     */
  public:

    ListEntryv10(Double value, bool isMin, bool isDummy = false,bool hasChanged = true){
      val_ = value;
      isMin_ = isMin;
      isDummy_ = isDummy;
      hasChanged_ = hasChanged;
    }
    ~ListEntryv10(){
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

    bool operator==(ListEntryv10 x) const{
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


  class RotListEntryv10 : public ListEntryv10
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
    RotListEntryv10(Double upperValue, Double lowerValue, Vector<Double> vector, std::list<ListEntryv10> switchingList, Double lastLocalXpar= 0.0,
      bool isMin = false, bool isDummy = false, bool wasWipedOut=false, UInt startCnt = 0)
      : ListEntryv10(upperValue, isMin, isDummy)
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

        critForVectorEquality_ = 1;
        tolForVectorEquality_ = 1e-14;
      }

      ~RotListEntryv10(){
      }

      /*
       * To keep the later rotation list as short as possible, we try to merge adjacent rotListEntries.
       * This is possible, if
       *  - they share the same rotation vector
       *  - they have the same switching list (e.g. after two adjacent regions were
       *     completely overwritten by a new input, which not only set the rotation vector
       *     but also set resulted in an xPar which overwrote the switching list)
       */
      bool canBeMergedWith(RotListEntryv10& rotEntry){

//        Double tol = 1e-15;
        Vector<Double> tmp = rotEntry.getVecCopy();

        bool vectorsEqual = Hysteresis::checkVectorEquality(vec_, tmp, critForVectorEquality_, tolForVectorEquality_);

        if(!vectorsEqual){
          return false;
        }

//        for(UInt i = 0; i< vec_.GetSize();i++){
//          if(abs(vec_[i]-tmp[i])>tol){
//            //std::cout << "Diffvec " << abs(vec_[i]-tmp[i]) << std::endl;
//            return false;
//          }
//        }

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

      /*
       * To keep the later rotation list as short as possible, we try to merge adjacent rotListEntries.
       * New merging rule:
       *        entry_i_last shall be the last entry of the switching list switchinglist_i belonging to rotlist i
       *        entry_i+1_first shall be the first entry of the switching list switchinglist_i+1 belonging to rotlist i+1
       *
       *        if length(switchinglist_i+1) == 1    <- always the case if switchinglist_i+1 was completely overwritten
       *          if type(entry_i_last) == type(entry_i+1_first)
       *            if(entry_i_last == entry_i+1_first)   <- new list starts where old list ended
       *              OR
       *              if(type(entry_i_last == max)
       *                if(entry_i_last > entry_i+1_first) AND entry_i+1_first == top value of bounding box
       *                  -> merge
       *              else
       *                if(entry_i_last < entry_i+1_first) AND entry_i+1_first == left value of bounding box
       *                  -> merge
       *
       */
      bool canBeMergedWith_v2(RotListEntryv10& rotEntry, bool classical){

        /*
         * this = current entry
         * rotEntry = next entry
         */

//        Double tol = 1e-15;
        Vector<Double> tmp = rotEntry.getVecCopy();

        bool vectorsEqual = Hysteresis::checkVectorEquality(vec_, tmp, critForVectorEquality_, tolForVectorEquality_);

        if(!vectorsEqual){
          return false;
        }

//        for(UInt i = 0; i< vec_.GetSize();i++){
//          if(abs(vec_[i]-tmp[i])>tol){
//            /*
//             * rotation direction does not coincide
//             */
//            return false;
//          }
//        }

        /*
         * new check for the new merging rule, which requires the second list to only contain one entry
         */
        if( rotEntry.getListReference().size() == 1 ){
          /*
           * check for new merging rule, which includes the old rule for the case that the second list
           * only has one entry
           * -> this merging rule is seldomly hit for the classical approach, as here the lists are longer
           * than one entry most of the time. The reason for this is, that the bounding boxes of the rotation
           * states are extending to the outer bounds at the left and top side so that overstanding entries are
           * not clipped down as it is the case for the revised model.
           */

          /*
           * check if current lists ends where the next list starts
           */
          Double lastEntry,firstEntry;
          bool lastEntryMin,firstEntryMin;

          lastEntry = switchingList_.back().getVal();
          lastEntryMin = switchingList_.back().isMin();
          firstEntry = rotEntry.getListReference().front().getVal();
          firstEntryMin = rotEntry.getListReference().front().isMin();

          if( lastEntryMin != firstEntryMin ){
            return false;
          }

          /*
           * switching list of current entry ends with same type of extremum as switching list of next entry begins
           */
          /*
           * check if values are the same
           */
          if(lastEntry == firstEntry){
            //std::cout << "Merged by new rule" << std::endl;
            return true;
          }

          /*
           * get left and top boundary of next rot entry
           */
          Double l,t;
          if(classical == true){
            l = -1.0;
            t = 1.0;
          } else {
            l = -rotEntry.getVal();
            t = rotEntry.getVal();
          }

          if(lastEntryMin){
            /*
             * if lastEntryMin == true -> check if lastEntry < firstEntry and firstEntry = left boundary of bounding box
             */
            // NEW 2018: <=
            if( (lastEntry <= firstEntry) && (firstEntry == l) ){
              //std::cout << "Merged by new rule" << std::endl;
              return true;
            }
          } else {
            /*
             * if lastEntryMin == false -> check if lastEntry > firstEntry and firstEntry = top boundary of bounding box
             */
            if( (lastEntry >= firstEntry) && (firstEntry == t) ){
              //std::cout << "Merged by new rule" << std::endl;
              return true;
            }
          }

        } else {
          /*
           * check the old merging rule (maybe we have identical lists, so why not merge?
           */
          if(switchingList_ == rotEntry.getListReference()){
            if(startCnt_ == rotEntry.getStartCnt()){
              /*
               * even though the switching lists may be the same,
               * one rotEntry starts evaluation at 0 the other one at >= 1.
               * But >= 1 means that the real first entry was deleted from the list
               * and thus the list wouldn't be equal
               */
              //std::cout << "Merged by old rule" << std::endl;
              return true;
            }
          }
        }

        return false;
      }

      std::list<ListEntryv10>& getListReference(){
        return switchingList_;
      }

      std::list<ListEntryv10> getListCopy(){
        return switchingList_;
      }

      void setList(std::list<ListEntryv10> switchingList){
        switchingList_ = switchingList;
      }

      std::string ToStringOLD() const {

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

        std::list<ListEntryv10>::const_iterator listIt;

        for(listIt = switchingList_.begin(); listIt != switchingList_.end(); listIt++){
          os << listIt->getValConst() << " (" << listIt->isMinConst() << ") -/- ";
        }

        return os.str();
      }

      std::string ToString(UInt cnt = 0) const {

        std::ostringstream os;

        /*
         * New format - 12.3.2020
         * 
         * cnt: x_thres (double), rotationState (vector), { [x_parallel (double), isMin (bool)]_1 [x_parallel (double), isMin (bool)]_2 ... [x_parallel (double), isMin (bool)]_n }
         * 
         */
        
        os << cnt << ": " << val_ << "; ("; 
        os << vec_[0] << ", " << vec_[1];
        if(vec_.GetSize() == 3){
          os << ", " << vec_[2];
        } 
        os  << "); {";
        
        std::list<ListEntryv10>::const_iterator listIt;
        
        for(listIt = switchingList_.begin(); listIt != switchingList_.end(); listIt++){
          os << "[" << listIt->getValConst() << "; ";
          if(listIt->isMinConst()){
            os << "min";
          } else {
            os << "max";
          }
          os << "]";
          if(listIt != switchingList_.end()){
            os << ";";
          }
        }
        os << "} \n";
        
        return os.str();
      }
      
      void setLowerVal(Double newValue){
        if(lowerVal_ != newValue){
          hasChanged_ = true;
        }
        /*
         * new addition: due to new treatment of upper triangle in the classical model, lowerVal_ is no longer the
         * right boundary of the bounding box in case of the last rotListEntry as the area corresponding to the
         * last rotListEntry is extended by the upper triangle AN2 (cont) plus a smaller triangle extending AN1.
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
         * When a new entry is added to the rotation list with xThresh = 0, we did not include it to the list
         * in older versions. In the new version, we have to add this entry to treat the changed rotation state
         * in AN2(cont) [which would become A(N+1)2]. This works perfectly fine except for one issue.
         * Due to saving computational resources, we try to reuse as many values as possible [which is why hasChanged_
         * was introduced]. If we consider the case of adding a new entry with value = zero, the old treatment fails.
         * The reason is, that the previous last entry [corresponding to AN1,AN2 and AN2(cont)] has a lower value of zero
         * already. If we insert a new entry with value 0, we set the lower value of entry N to 0, but this value is
         * already set, so the hasChanged_ flag stays false. Therewith, entry N will not be reevaluted so that the
         * loaded value still includes AN1, AN2 and AN2(cont). This is wrong as AN2(cont) as well as the extension of
         * AN1 into the upper triangle must be subtracted to get correct results.
         *
         * Short conclusion: check if newValue = 0, if yes -> set hasChanged_ to true
         */
        if(newValue == 0){
          hasChanged_ = true;
        }

        lowerVal_ = newValue;
      }

      void setVec(Vector<Double> newVector){

//        Double tol = 1e-15;
        rotHasChanged_ = !Hysteresis::checkVectorEquality(vec_, newVector, critForVectorEquality_, tolForVectorEquality_);

//        for(UInt i = 0; i< vec_.GetSize();i++){
//          if(abs(vec_[i]-newVector[i])>tol){
//            rotHasChanged_ = true;
//            break;
//          }
//        }

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

        std::list<ListEntryv10>::iterator listIt;

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

        std::list<ListEntryv10>::iterator listIt;

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
         * Note: we do not set flag hasChanged to true here
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
    std::list<ListEntryv10> switchingList_;

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

    UInt critForVectorEquality_;
    Double tolForVectorEquality_;
  };


  class VectorPreisachSutor : public Hysteresis
  {

  public:
    VectorPreisachSutor(Integer numElem, ParameterPreisachOperators operatorParams,
          ParameterPreisachWeights weightParams, UInt dim, bool isVirgin);

//      Integer numElem, Double xSat, Double ySat,
//      Matrix<Double>& preisachWeight, Double rotationalResistance , UInt dim, bool isVirgin,
//      bool classical, bool scaleUpToSaturation, Double angularDistance, Double angResolution,
//      Double anhystA, Double anhystB, Double anhystC, bool anhystOnly);

    virtual ~VectorPreisachSutor();

    //! this function gets called from outside and calculates the output of the Preisach operator
    virtual Vector<Double> computeValue_vec(Vector<Double>& xVal, Integer idElem, bool overwrite,
      bool debugOutput, int& successFlag, bool skipAnhystPart = false){
      EXCEPTION("Not implemented in base class");
    }

    //void ClipDirection(Vector<Double>& targetVector);

    //
    /*
     * get-function should no longer be used
     * background: coefFunctionHyst allows for "extended" evaluation, i.e.
     *  each element has one hysteresis operator (addresable by index operatorIdx)
     *  that gets evaluated at each integration point separately; the results
     *  of these evaluations is stored separately, too (addresable by index storageIdx)
     * > the hyst operator does not know about this additional storageIdx and such
     *   it cannot return the correct value when called by getValue_vec
     * > when comparisons to previous solutions are required (i.e. check for changes)
     *   this should be done in coefFunctionHyst
     * > exception: for inversion via LevenbergMarquart, we need to know the previous
     *   values at the storageIdx > those values should be passed to the compute input
     *   function
     */
//    //! returns the current output of the hyst-operator for element idxElem
//    Vector<Double> getValue_vec( Integer idElem, bool overwrite = true )
//    {
//      if(overwrite == false){
//        return (preisachSumTmp_[idElem]);
//      } else {
//        return ( preisachSum_[idElem] );
//      }
//    }



    //! Try to compute input xVal to hyst operator, such that mu*xVal + H(xVal) = yVal
    // return usable input xVal
    /*
     * computeInput_vec is the one to be called in coefFunctionHyst
     * Exception: testInversion > here we use computeInput_vec_withStatistics
     */
    Vector<Double> computeInput_vec(Vector<Double> yVal, Integer operatorIndex,
      Matrix<Double> mu, bool fieldsAlignedAboveSat, bool hystOutputRestrictedToSat,
      int& successFlag,bool useEverett, bool overwrite){

      Vector<Double> prevYval = Vector<Double>(dim_);
      mu.Mult(prevXVal_[operatorIndex],prevYval);
      prevYval.Add(1.0,prevHVal_[operatorIndex]);

      return computeInput_vec_withPrevStates(yVal, prevYval,
        prevXVal_[operatorIndex], prevHVal_[operatorIndex], operatorIndex,
        mu, fieldsAlignedAboveSat, hystOutputRestrictedToSat, successFlag);
    }

    void switchingStateToBmp(UInt numPixel, std::string filename, UInt idElem, bool overLayWithRotState){
      EXCEPTION("Not implemented in base class");
    }
    void rotationListToTxt(std::string filename, UInt idElem, bool append, std::string optionalHeader){
      EXCEPTION("Not implemented in base class");
    }
    
    void collectParallelProjections(bool switchOnOff, UInt cnt){    
      collectProjections_ = switchOnOff;
      tsCNTForProjections_ = cnt;
    }

    void setFlags(UInt performanceFlag){

      if(performanceFlag == 2){
        /*
         * start with 2
         * -> Reason: coefFunctionHyst uses performanceFlag with
         *      0 > no measurement
         *      1 > measurement but only in coefFunctionHyst
         *      2 > measurement also in Preisach operator
         */
        performanceMeasurement_ = true;
      } else {
        performanceMeasurement_ = false;
      }
    }

    std::string runtimeToString(){
      EXCEPTION("Not implemented in base class");
    }

  protected:

    Vector<Double> evaluateNewRotationDirection(Vector<Double>& e_u_new, Vector<Double>& e_u_old, Double xVal, UInt idElem);

    Vector<Double> restrictToHalfspace(Vector<Double>& e_u_new);

    /*!
     * Global quantities, i.e. the same for all FE elements of the same material
     */
//    Double XSaturated_; //! saturation value for  input
//    Double PSaturated_; //! saturation value for output
    Double maxOutputVal_; // actual maximial output if input is in or above saturation
    // usually = PSaturated, but not for rotResistance < 1 and revised model

    Matrix<Double> preisachWeights_; //! preisach weight function
    Double rotationalResistance_; //! parameter describing the resistance of the domains to rotation

    bool collectProjections_;
    std::list< std::string > listOfCollectedProjections_;
    UInt tsCNTForProjections_;
    
    bool testAngDistForClassical;
    bool isVirgin_; //! yes, if starting at zero; currently flag is unused
    bool isSymmetric_; //! states if Preisach weights are symmetric w.r.t alpha = -beta

//    UInt dim_; //! 2D or 3D (axi not tested at all!) // in bace class
    UInt numElem_; //! total number of FE elements
    UInt numRows_; //! number of rows of the Preisach plane

//    Double angResolution_; //! tolerance for resoultion of angular clipping; provided via mat.xml > no longer used
    Double delta_; //! resolution of Preisach plane
    Double tol_; //! tolerance for all kind of comparisons
    UInt critForVectorEquality_;
    Double tolForVectorEquality_;

    Double tolSwitchingEntry_;
    Double tolForWeightComparison_;

//    /*
//     * for optional anhysteretic parts
//     */
//    Double anhyst_A_;
//    Double anhyst_B_;
//    Double anhyst_C_;

    bool classical_; //! switch between classical evaluation (2012 version of vector model) or revised evaluation (2015 version)
    bool restrictToHalfspace_;
    bool scaleUpToSaturation_; // see computeValue_vec for details

    /*!
     * Local quantities, i.e. one for each FE elements of the same material
     * NOTE: these values may only be written by function
     *        computeValue_vec and ONLY if flag overwrite=true
     */
    Vector<Double>* preisachSum_;
    Vector<Double>* preisachSumTmp_;

    /*!
     * Quantities needed by the revised approach (Sutor2015)
     */
    Double angularDistance_; //! angular distance between old rotation state and new input direction which remains for small fields

    /*
     * runtime measurement
     */
    bool performanceMeasurement_;

    /*
     * can be used for output
     * O = no output
     * 1 = info
     * 2 = debug
     */
    UInt textOutputLevel_;


    /*
     * allows to switch between new and old mapping function
     * (mapRectangelTo...)
     *
     * 0 = old version
     * 1 = new version (default)
     */
    UInt mappingVersion_;

    /*
     * for revised version -> allow precomputation of deltaPhi as well as cos(deltaPhi), sin(deltaPhi)
     * for computation of new rotation state
     * > should be set for each fe element separately
     * > race conditions!
     */
    bool usePreComputedValue_;
    Vector<Double> deltaPhi_preComputed_;
    Vector<Double> cos_deltaPhi_preComputed_;
    Vector<Double> sin_deltaPhi_preComputed_;
  };

  class VectorPreisachSutor_MatrixApproach : public VectorPreisachSutor
  {
  public:
    VectorPreisachSutor_MatrixApproach(Integer numElem, ParameterPreisachOperators operatorParams,
          ParameterPreisachWeights weightParams, UInt dim, bool isVirgin);

//      Integer numElem, Double xSat, Double ySat,
//      Matrix<Double>& preisachWeight, Double rotationalResistance , UInt dim, bool isVirgin,
//      bool classical, bool scaleUpToSaturation, Double angularDistance, Double angResolution,
//      Double anhystA, Double anhystB, Double anhystC, bool anhystOnly);

    ~VectorPreisachSutor_MatrixApproach();

    //! this function gets called from outside and calculates the output of the Preisach operator
    Vector<Double> computeValue_vec(Vector<Double>& xVal, Integer idElem, bool overwrite,
      bool debugOut, int& successFlag, bool skipAnhystPart = false);

    void switchingStateToBmp(UInt numPixel, std::string filename, UInt idElem, bool overLayWithRotState);
    void rotationListToTxt(std::string filename, UInt idElem, bool append, std::string optionalHeader){
      std::cout << "+++ INFO: Matrix based implementation has no rotation list to output." << std::endl;
    }
        
    std::string runtimeToString();

  private:

    void InitializeSwitchingState(UInt idElem); //! initial splitting of Preisach plane into +1 and -1 part

    //! update function for rotation state
    void UpdateRotationStates(Double XThres, Double xVal, Vector<Double>& e_u_new, UInt idElem);

    //! update function for switching state
    void UpdateSwitchingStates(Vector<Double>& u_in, UInt idElem);

    Matrix<Double>* switchingStates_; //! stores the current switching state
    Matrix<Double>* rotationStateX_; //! stores the current rotation state in x
    Matrix<Double>* rotationStateY_; //! stores the current rotation state in y
    Matrix<Double>* rotationStateZ_; //! stores the current rotation state in z

    /*
     * Timer and counter
     */
    UInt updateMatricesCounter_;
    Timer* updateMatricesTimer_;

    UInt evaluateMatricesCounter_;
    Timer* evaluateMatricesTimer_;

    UInt copyToTemporalStorageCounter_;
    Timer* copyToTemporalStorageTimer_;
    Timer* copyFromTemporalStorageTimer_;

  };

  class VectorPreisachSutor_ListApproach : public VectorPreisachSutor
  {

  public:
    VectorPreisachSutor_ListApproach(Integer numElem, ParameterPreisachOperators operatorParams,
          ParameterPreisachWeights weightParams, UInt dim, bool isVirgin);

//      Integer numElem, Double xSat, Double ySat,
//      Matrix<Double>& preisachWeight, Double rotationalResistance , UInt dim, bool isVirgin,
//      bool classical, bool scaleUpToSaturation, Double angularDistance, Double angResolution,
//      Double anhystA=0.0, Double anhystB=0.0, Double anhystC=0.0, bool anhystOnly = false);

    virtual ~VectorPreisachSutor_ListApproach();

    //! this function gets called from outside and calculates the output of the Preisach operator
    Vector<Double> computeValue_vec(Vector<Double>& xVal, Integer idElem, bool overwrite,
      bool debugOut, int& successFlag, bool skipAnhystPart = false);

//    Vector<Double> computeValue_vecMeasure(Vector<Double>& u_in, Integer idElem, bool overwrite,
//          bool debugOut, int& successCode, Double& time);

    void switchingStateToBmp(UInt numPixel, std::string filename, UInt idElem, bool overLayWithRotState);
    void rotationListToTxt(std::string filename, UInt idElem, bool append, std::string optionalHeader);
    
    std::string runtimeToString();

  private:

    /*
     * for version 10 -> revised model
     */
    void Update_GlobalRotationList(Double xThres, Double xVal, Vector<Double> e_u,
      std::list<RotListEntryv10>& usedRotationList, UInt idElem, bool performSimplification, bool debutOutput = false);
    UInt Update_SwitchingList(std::list<ListEntryv10>& list, Double newEntry, Double lastXpar, Rectangle boundingBox, bool wasWipedOut, bool lastRotEntry);
    Double clipRectangleToElement(Rectangle& source, UInt idAlpha, UInt idBeta, Double delta = -1,bool isRotState = false);
    void getBoundingBoxFromRotEntry(std::list<RotListEntryv10>::iterator rotListIt, Rectangle& rect1, bool lastRotListEntryv10);
    bool getRectanglesFromRotEntry(std::list<RotListEntryv10>::iterator rotListIt, Rectangle& rect1, Rectangle& rect2, bool lastRotListEntryv10);
    void Evaluate_GlobalRotationList(std::list<RotListEntryv10>& usedRotationList, Vector<Double>& retVec, bool debugOutput = false);
    void Evaluate_LowerTriangle();
    Double mapRectangleToPreisachWeightsOLD(Rectangle& rect, bool skipUpperDiagonal = false);
    Double mapRectangleToPreisachWeightsNEW(Rectangle& rect, bool skipUpperDiagonal = false);
    Double mapRectangleToPreisachWeightsDEBUG(Rectangle& rect, bool skipUpperDiagonal = false);
    Double mapRectangleToPreisachWeights(Rectangle& rect, bool skipUpperDiagonal = false);
    Double getRectangleFromSwitchingList(std::list<ListEntryv10>& list,
    std::list<ListEntryv10>::iterator startIt, std::list<ListEntryv10>::iterator curIt, std::list<ListEntryv10>::iterator endIt,
    UInt idArea, Rectangle& rect, bool upperSplitSquare = false);
    bool Simplify_LocalSwitchingLists(std::list<RotListEntryv10>& usedRotationList);
    bool Simplify_GlobalRotationList(std::list<RotListEntryv10>& usedRotationList);
    void mapRectangleToHelperMatrix(Matrix<Double>& helper, Rectangle rect, Double factor, bool skipUpperDiagonal,bool isRotState);
    void Initialize_GlobalRotationList(std::list<RotListEntryv10>& usedRotationList);
    void Initialize_GlobalRotationListWithValues(std::list<RotListEntryv10>& usedList,Vector<Double>& initDir, Double initRotValue, Double initSwitchValue);


    /*!
     * Local quantities, i.e. one for each FE elements of the same material
     */
    std::list<RotListEntryv10>* globRotList_; //! outer list from which rotation and switching state get computed

    /*!
     * Quantities needed by the classical_ approach (Sutor2012)
     */
    Vector<Double>* lastEu_; //! last rotation state; only used in output function rotationStateToBmp

    Double lowerTriangleValue_; //! Integral over the Preisach weights in the lower triangular region

    /*
     * Timer and counter
     */
    UInt updateNestedListCounter_;
    Timer* updateRotListTimer_;
    Timer* updateSwitchingListTimer_;
    Timer* simplifyRotListTimer_;
    Timer* simplifySwitchingListTimer_;

    UInt evaluateNestedListCounter_;
    Timer* evaluateNestedListTimer_;

    UInt copyToTemporalStorageCounter_;
    Timer* copyToTemporalStorageTimer_;


  };
} // end of namespace

#endif
