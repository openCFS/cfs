#ifndef FILE_VECPREISACHV10_2016
#define FILE_VECPREISACHV10_2016

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
      
      os << "Left / Right: " << l_ << " / " << r_ << '\n';
      os << "Top / Bottom: " << t_ << " / " << b_ << '\n';
      
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
      
      t = fmin(this->t_,source.t_);
      b = fmax(this->b_,source.b_);
      l = fmax(this->l_,source.l_);
      r = fmin(this->r_,source.r_);
      
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
      
      if(fabs(val_ - x.getVal()) > tol){
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
        
        Double tol = 1e-15;
        Vector<Double> tmp = rotEntry.getVecCopy();
        
        for(UInt i = 0; i< vec_.GetSize();i++){
          if(fabs(vec_[i]-tmp[i])>tol){
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
        
        Double tol = 1e-15;
        Vector<Double> tmp = rotEntry.getVecCopy();
        
        for(UInt i = 0; i< vec_.GetSize();i++){
          if(fabs(vec_[i]-tmp[i])>tol){
            /*
             * rotation direction does not coincide
             */
            return false;
          }
        }
        
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
            if( (lastEntry < firstEntry) && (firstEntry == l) ){
              //std::cout << "Merged by new rule" << std::endl;
              return true;
            }
          } else {
            /*
             * if lastEntryMin == false -> check if lastEntry > firstEntry and firstEntry = top boundary of bounding box
             */
            if( (lastEntry > firstEntry) && (firstEntry == t) ){
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
        
        std::list<ListEntryv10>::const_iterator listIt;
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
        
        Double tol = 1e-15;
        for(UInt i = 0; i< vec_.GetSize();i++){
          if(fabs(vec_[i]-newVector[i])>tol){
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
  };
  
  
  class VectorPreisachv10 : public Hysteresis
  {
    
  public:
    VectorPreisachv10(Integer numElem, Double xSat, Double ySat,
      Matrix<Double>& preisachWeight, Double rotationalResistance , UInt dim, bool isVirgin,
      bool classical, Double angularDistance, Double angularClipping);
    
    virtual ~VectorPreisachv10();
    
    //! this function gets called from outside and calculates the output of the Preisach operator
    virtual Vector<Double> computeValue_vec(Vector<Double>& xVal, Integer idElem, bool overwrite = true,bool overwriteDirection = true){
      EXCEPTION("Not implemented in base class");
    }
    
    //! returns the current output of the hyst-operator for element idxElem
    Vector<Double> getValue_vec( Integer idElem, bool overwrite = true )
    {
      if(overwrite == false){
        return (preisachSumTmp_[idElem]);
      } else {
        return ( preisachSum_[idElem] );
      }
    }
    
    //! returns the current input of the hyst-operator for element idxElem
    Vector<Double> getInput_vec( Integer idElem )
    {
      return actXval_[idElem];
    }
    
    bool checkConvergence(Vector<Double>& res, Matrix<Double>& jacT, Double& errorNorm, Double tol){
      // According to Dahmen&Reusken - Numerik partieller DFG
      // the residual is a non-sufficient condition
      // instead we have to check the norm of jacT*res
      Vector<Double> errorVec = Vector<Double>(dim_);
      jacT.Mult(res,errorVec); 
      errorNorm = errorVec.NormL2();
      
      std::cout << "CheckConvergence" << std::endl;
      std::cout << "resIn: " << res.ToString() << std::endl;
      std::cout << "jacT: " << jacT.ToString() << std::endl;
      std::cout << "ErrorVec: " << errorVec.ToString() << std::endl;
      
      if(errorNorm <= tol){
        return true;
      } else {
        return false;
      }
      
    }
    
    bool checkIncrement(Vector<Double>& xUpdate, Vector<Double>& res, Vector<Double>& resShifted, Matrix<Double>& jac, Double& alpha){
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
      
      std::cout << "res: " << res.ToString() << std::endl;
      std::cout << "resShifted_: " << resShifted.ToString() << std::endl;
      std::cout << "xUpdate: " << xUpdate.ToString() << std::endl;
      std::cout << "jac: " << jac.ToString() << std::endl;
      
      tmp = tmp + res;
      Double tmpNorm = tmp.NormL2();
      Double a = resNorm*resNorm - resShiftedNorm*resShiftedNorm;
      Double b = resNorm*resNorm - tmpNorm*tmpNorm;
      std::cout << "resNorm*resNorm - resShiftedNorm*resShiftedNorm: " << a << std::endl;
      std::cout << "resNorm*resNorm - tmpNorm*tmpNorm: " << b  << std::endl;
      
      Double rho;
      if(b != 0){
        rho = a/b;
      } else {
        rho = a;
      }
      
      // 0 < beta0 < beta1 < 1
      Double beta0 = 0.25;
      Double beta1 = 0.75;
      
      std::cout << "rho: " << rho << std::endl;
      std::cout << "current alpha: " << alpha << std::endl;
      // test different stepping; the furhter away we are, the stronger
      // we will adapt stepping alpha for next time
      bool success = false;

      // increment not accepted; increase alpha in linesearch
      if(rho < -100*beta0){
        alpha = alpha*256.0;
      } else if(rho < -10*beta0){
	alpha = alpha*128.0;
      } else if(rho < -beta0){
	alpha = alpha*64.0;
      } else if(rho < 0){
        alpha = alpha*32.0;
      } else if(rho < beta0/8.0){
        alpha = alpha*16.0;
      } else if(rho < beta0/4.0){
        alpha = alpha*8.0;
      } else if(rho < beta0/2.0){
        alpha = alpha*4.0;
      } else if(rho < beta0/1.5){
	alpha = alpha*2.0;
      } else if(rho < beta0){
        alpha = alpha*1.5; 
      } 
      // increment accepted; decrease alpha for next linesearch call
      // usually rho does not become much larger than 1 (from theory 1 = max)#
      else if(rho > 1.8*beta1){
	alpha = alpha/16.0;
	success = true;
      } else if(rho > 1.4*beta1){
        alpha = alpha/8.0;
        success = true;
      } else if(rho > 1.2*beta1){
        alpha = alpha/4.0;
        success = true;
      } else if(rho > 1.1*beta1){
	alpha = alpha/2.0;
	success = true;
      } else if(rho > beta1){
        alpha = alpha/1.5;
        success = true;
      } else {
        // best case
        // increment ok; keep alpha as startvalue for next linesearch
        success = true;
      }
      std::cout << "new alpha: " << alpha << std::endl;
      return success;
    }
    
    Vector<Double> computeAbsResidualX(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& hystVal, Double mu){
      /*
       *    yVal = mu*xVal + hystVal(xVal)
       * 
       *  > res = xVal + hystVal(xVal)/mu - yVal/mu
       */
      Vector<Double> res = Vector<Double>(dim_);
      res = xVal;
      res.Add(1.0/mu,hystVal);
      res.Add(-1.0/mu,yVal);
      
      return res;
    }
        
    Vector<Double> computeResidual(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& hystVal, Double mu, Integer idElem, 
            bool wrtX, bool relative){

      Vector<Double> ret = computeAbsResidualX(xVal, yVal, hystVal, mu);
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
            ret = ret/xValNorm; 
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
        ret = ret*(-mu);
        if(relative){
          Double yValNorm = yVal.NormL2();
          if(yValNorm != 0){
            ret = ret/yValNorm; 
          } else {
            WARN("relative residual wrt y requested, but yValNorm == 0; return abs value instead");
          }
        }
      }
     
      return ret;
    }
    
    Matrix<Double> computeJacobianOfAbsResidualX(Vector<Double>& xVal, Vector<Double>& hyst, Double mu, Integer idElem, Double sign) {
      
      //    std::cout << "CompueJacobian" << std::endl;
//      Double deltaXmin = 1e-5;
//      Double scal = 1e-5;
      Double deltaXmin = 1e-10;
      Double scal = 1e-10;
      Double deltaX;
      
      bool overwrite = false;
      Vector<Double> xShifted;
      Vector<Double> hystShifted;
      
      Matrix<Double> jac = Matrix<Double>(dim_,dim_);
      
      std::cout << "xVal: " << xVal.ToString() << std::endl;
      std::cout << "hystVal: " << hyst.ToString() << std::endl;
      
      for(UInt i = 0; i < dim_; i++){
        xShifted = xVal;
        
        if(xVal.NormL2() >= XSaturated_){
          // last value is in saturation but yVal seems not (as we otherwise would
          // end up here);
          // we have to reduce xVal for the computation of Jacobian, as we otherwise
          // might get deltaY = hyst(xVal+delta) - hyst(xVal) = ySat-ySat = 0
          // 1. set xShifted to xVal > already donw
          // 2. scale xShifted by XSaturated_/xShifted.NormL2()
          xShifted *= (XSaturated_/xShifted.NormL2());
          
          sign = -1.0;
        } 
        
//        if( xVal[i] < 0 ){
//          deltaX = sign*fmin( scal*xVal[i], -deltaXmin );
//        } else {
//          deltaX = sign*fmax( scal*xVal[i], deltaXmin );
//        }
        
        if( xVal[i] < 0 ){
          deltaX = sign*fmin( -scal*XSaturated_, -deltaXmin );
        } else {
          deltaX = sign*fmax( scal*XSaturated_, deltaXmin );
        }

        xShifted[i] += deltaX;
        
        hystShifted = computeValue_vec(xShifted, idElem, overwrite);
        
        /*
         * Compute Jacobian for residual wrt x
         * 
         * jac_ji = -delta_ji - dhystVal_j/dxVal_i/mu
         */   
        jac[i][i] = 1.0;
        for(UInt j = 0; j < dim_; j++){
          jac[j][i] += (hystShifted[j]-hyst[j])/deltaX/mu; 
        }
      }
      return jac;
    }
   
    Matrix<Double> computeJacobian(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& hyst, Vector<Double>& resX,
            Double mu, Integer idElem, Double sign, bool wrtX, bool relative){
      
      Matrix<Double> jac = computeJacobianOfAbsResidualX(xVal, hyst, mu, idElem, sign);
      
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
        jac = jac*(-mu);
        
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
    
    Vector<Double> computeResidualOLD(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& hystVal, Double mu, Integer idElem, bool wrtX){  
      
      Vector<Double> res = Vector<Double>(dim_);
      res = yVal;
      res -= hystVal;
      
      std::cout << "Compute Residual: " << std::endl;
      std::cout << "xVal: " << xVal.ToString() << std::endl;
      std::cout << "yVal: " << yVal.ToString() << std::endl;
      std::cout << "hystVal: " << hystVal.ToString() << std::endl;

      if(wrtX){
        /*
         * Residual wrt x
         */
        // NOTE: hystVal = electric or magnetic POLARIZATION
        // yVal = mu*xVal + hystVal
        // > res = (yVal - hystVal)/mu - xVal
        // residual is computed wrt xVal
        // to get residual wrt yVal, we have to multiply residual by mu
        res /= mu;
        res -= xVal;
        
      } else {
        /*
         * Residual wrt y
         * 
         * NOTE: Error in x will be larger than if residual is reduced wrt to x
         * BUT: much faster convergence!
         * 
         * res = yVal - hystVal - mu*xVal
         */
         res.Add(-mu,xVal);
      }
      return res;
    } 
    
    Matrix<Double> computeJacobianOLD(Vector<Double>& xVal, Vector<Double>& hyst, Double mu, Integer idElem, Double sign,bool wrtX){
      
      //    std::cout << "CompueJacobian" << std::endl;
      Double deltaXmin = 1e-5;
      Double scal = 1e-5;
      Double deltaX;
      
      bool overwrite = false;
      Vector<Double> xShifted;
      Vector<Double> hystShifted;
      
      Matrix<Double> jac = Matrix<Double>(dim_,dim_);
      
      std::cout << "xVal: " << xVal.ToString() << std::endl;
      std::cout << "hystVal: " << hyst.ToString() << std::endl;
      
      for(UInt i = 0; i < dim_; i++){
        xShifted = xVal;
        
        if(xVal.NormL2() >= XSaturated_){
          // last value is in saturation but yVal seems not (as we otherwise would
          // end up here);
          // we have to reduce xVal for the computation of Jacobian, as we otherwise
          // might get deltaY = hyst(xVal+delta) - hyst(xVal) = ySat-ySat = 0
          // 1. set xShifted to xVal > already donw
          // 2. scale xShifted by XSaturated_/xShifted.NormL2()
          xShifted *= (XSaturated_/xShifted.NormL2());
          
          sign = -1.0;
        } 
        
//        if( xVal[i] < 0 ){
//          deltaX = sign*fmin( scal*xVal[i], -deltaXmin );
//        } else {
//          deltaX = sign*fmax( scal*xVal[i], deltaXmin );
//        }
        
        if( xVal[i] < 0 ){
          deltaX = sign*fmin( -scal*XSaturated_, -deltaXmin );
        } else {
          deltaX = sign*fmax( scal*XSaturated_, deltaXmin );
        }
        
        std::cout << "deltaX " << deltaX << std::endl;
        
        xShifted[i] += deltaX;
        
        hystShifted = computeValue_vec(xShifted, idElem, overwrite);
        
        if(wrtX){
          /*
           * Compute Jacobian for residual wrt x
           * 
           * jac_ji = -delta_ji - dhystVal_j/dxVal_i/mu
           */   
          jac[i][i] = -1.0;
          for(UInt j = 0; j < dim_; j++){
            jac[j][i] -= (hystShifted[j]-hyst[j])/deltaX/mu; 
          }
        } else {
          /*
           * Compute Jacobian for residual wrt y
           * 
           * jac_ji = -mu*delta_ji - dhystVal_j/dxVal_i
           */ 
          for(UInt j = 0; j < dim_; j++){
            jac[j][i] = -(hystShifted[j]-hyst[j])/deltaX; 
          }
          
          jac[i][i] -= mu;
        }
      }
      return jac;
    }  
    
    bool performLinesearch(Vector<Double>& xVal, Vector<Double>& yVal, Vector<Double>& res, Vector<Double>& xUpdate,
          Matrix<Double>& jac, Matrix<Double>& jacT, Double mu, Integer idElem,Double& alpha, 
          bool wrtX, bool relative){
      
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
      jacTres_neg = jacTres_neg*Double(-1.0);
      
      Double alphaMax = 1e8;
      Double alphaMin = 1e-16;
      
      bool success = false;
      bool discard = false;
      //    Double bestResidual = 1e15;
      //    Double bestAlpha = -1.0;
      //    Vector<Double> bestUpdate = Vector<Double>(dim_);
      
      while(true){
        itCnt++;
        
        matToInvert = jacTjac;
        for(UInt i = 0; i < dim_; i++){
          matToInvert[i][i] +=  alpha*alpha;
        }
        
        matToInvert.Invert(matInverted);
        matInverted.Mult(jacTres_neg,xUpdate);
        
        xNew = xVal;
        xNew += xUpdate;
        
        hystNew = computeValue_vec(xNew, idElem, false);
        resNew = computeResidual(xNew,yVal,hystNew,mu,idElem,wrtX,relative);

        //      if(resNew.NormL2() < bestResidual){
        //        bestResidual = resNew.NormL2();
        //        bestAlpha = alpha;
        //        bestUpdate = xUpdate;
        //      }
        
        std::cout << "--Check increment-- " << std::endl;
        std::cout << "Current alpha: " << alpha << std::endl;
        std::cout << "y,xOld,hystOld: " << std::endl;
        std::cout << yVal.ToString() << std::endl;
        std::cout << xVal.ToString() << std::endl;
        Vector<Double> hystOLD = computeValue_vec(xVal, idElem, false);
        std::cout << hystOLD.ToString() << std::endl;
        std::cout << "y,xOld,hystOld: " << std::endl;
        std::cout << yVal.ToString() << std::endl;
        std::cout << xNew.ToString() << std::endl;
        std::cout << hystNew.ToString() << std::endl;
        std::cout << "Current alpha: " << alpha << std::endl;
        
        success = checkIncrement(xUpdate, res, resNew, jac, alpha);
        
        if(alpha > alphaMax){
          // maximal alpha used; stop here (regardless of success
          alpha = alphaMax;
          break;
        }
        if(alpha < alphaMin){
          // minimal alpha used; stop here
          alpha = alphaMin;
          break;
        }
        
        if(success){
          break;
        } else {
          //        if((alpha == alphaMin) || (alpha == alphaMax)){
          //          // nothing more to achieve > stop
          //        }
          //        if(check == 2){
          //          std::cout << "Rho = 0 > current update has no effect; instead take best update so far" << std::endl;
          //          alpha = bestAlpha;
          //          xUpdate = bestUpdate;
          //          break;
          //        }
          
          if(itCnt >= maxIter){
            //          std::cout << "Max number of iterations used" << std::endl;
            //          std::cout << "Reuse best intermediate solution (the one with the smallest residual)" << std::endl;
            //          alpha = bestAlpha;
            //          xUpdate = bestUpdate;
            discard = true;
            break;
          }
        }
      }
      
      if(success){
        std::cout << "Linesearch was successful" << std::endl;
        std::cout << "newAlpha: " << alpha << std::endl;
      } else {
        std::cout << "Linesearch was NOT successful" << std::endl;
        std::cout << "newAlpha: " << alpha << std::endl;
        
      }
      if(discard){
        std::cout << "Rho < 0: Discard xUpdate: " << xUpdate.ToString() << std::endl;
      }
      
      return discard;
      
    }
    
    //! Try to compute input xVal to hyst operator, such that mu*xVal + H(xVal) = yVal
    // return usable input xVal
    Vector<Double> computeInput_vec(Vector<Double>& yVal, Integer idElem, Double mu, Double& alpha, 
          bool overwrite = true,bool overwriteDirection = true){
      
      
      Vector<Double> xVal = Vector<Double>(dim_);
      
      /*
       * Check if yVal is beyond saturation > easy case
       */
      Double yNorm = yVal.NormL2();
      
      std::cout << "yNorm: " << yNorm << std::endl;
      std::cout << "yNorm-YSaturated_: " << yNorm-YSaturated_ << std::endl;
      if(yNorm >= YSaturated_){
        std::cout << "Use simple approach" << std::endl;
        // Important consequences:
        // a) material is completely aligned with outer field (at least in Sutors model and therefore also x
        // b) mu (eps) adds an important contribution
        // c) xVal = (yVal - ySat+xSat*mu)/mu + xSat
        //         = (yVal - ySat)/mu
        Vector<Double> ySat = yVal;
        ySat *= (YSaturated_/yNorm);
        xVal = yVal;
        xVal -= ySat;
        xVal /= mu;
      } 
      /*
       * Use Levenberg-Marquart algorithm as presented in Dahmen&Reusken - Numerik partialler DFG
       */   
      else {
        std::cout << "Try LM" << std::endl;
        
        // tolerance wrt y > 1e-10 or 1e-12 seems good > takes 2-3 its
        // only problem: y-x-loops look ugly as x can be quite off!
        //Double tolError = 1e-11;  
        // tolerance for reevalution
        Double tolY = 1e-16;
        bool wrtX = true;
        bool relError = !true;
        
        // toleranace for relative error criterion
        Double tolError = 0.01;
        if(!relError){
          tolError = tolError*XSaturated_;
        }
//        if(wrtX){
//          tolError = 0.01;
//          // if we reduce wrt X use larger tol
//          // > error for y will be similar small
//          if(relError){
//            
//          } else {
//            tolError = 10;
//          }
//          //tolError /= mu;
//        } else {
//          if(relError){
//            tolError = 1e-11;
//          } else {
//            tolError = 1e-9;
//          }        
//        }
        
        UInt maxIter = 25;
        UInt itCnt = 0;
        
        // use last computed Xval as starting ppoint
        xVal = actXval_[idElem];
        
        Vector<Double> diff = yVal;
        diff -= actYval_[idElem];
        
        if(diff.NormL2() < tolY){
          std::cout << "Difference between requested yVal and previously computed yVal < " << tolY << std::endl;
          return xVal;
        }
        
        Vector<Double> xUpdate = Vector<Double>(dim_);
        
        // check if a starting value for the linesearch paramater alpha was given
        // if not: set some starting value (found out by testing to be quite ok)
        if(alpha < 0.0){
          if(wrtX){
            // larger alpha needed
            alpha = 10/sqrt(mu);;
          } else {
            // small alpha needed
            alpha = 0.1*sqrt(mu);
          }
        }
        
        Double sign = 1.0;
        bool success = false;
        bool discardUpdate = false;
        
        Vector<Double> res = Vector<Double>(dim_);
        Matrix<Double> jac = Matrix<Double>(dim_,dim_);
        Matrix<Double> jacT = Matrix<Double>(dim_,dim_);
        Vector<Double> hyst;
        Double errorNorm;
        
        while(true){
          itCnt++;
          // do not override here
          hyst = computeValue_vec(xVal, idElem, false);
          
          bool relative = relError;
          res = computeResidual(xVal,yVal,hyst,mu,idElem,wrtX,relative);
          jac = computeJacobian(xVal,yVal,hyst,res,mu,idElem,sign,wrtX,relative);
          jac.Transpose(jacT);
          
//          Double tolForCheck = tolError;
//          if(relError && (xVal.NormL2() != 0)){
//            // absError = |jacT*res|
//            // relError = |jacT*res|/|xVal|
//            // > to check for relError < tolError, we can instead check 
//            //   if absError < tolError*|xVal|
//            tolForCheck = tolError*xVal.NormL2();
//            std::cout << "Check absError against tol: " << tolForCheck << std::endl;
//          }
          
          if(relError) {
            std::cout << "Check relError against tol: " << tolError << std::endl;
          } else {
            std::cout << "Check absError against tol: " << tolError << std::endl;
          }

          success = checkConvergence(res,jacT,errorNorm,tolError);
          
          std::cout << "Actual errorCrit: " << errorNorm << std::endl;
          std::cout << "Actual residual ";
          if(wrtX){
            std::cout << "(wrtX) ";
          }
          std::cout << res.NormL2() << std::endl;
          
          if(success){
            std::cout << "+++++++++++++++++++++++++++" << std::endl;
            std::cout << "Inversion success after it " << itCnt << std::endl;
            //std::cout << "Last tol for convergence: " << tolForCheck << std::endl;
            if(wrtX){
              std::cout << "Remaining error-Norm wrt xVal: " << errorNorm << std::endl;
              Vector<Double> resY = res;
              resY *= mu;
              std::cout << "Remaining error-Norm wrt yVal: " << resY.NormL2() << std::endl;
            } else {
              std::cout << "Remaining error-Norm wrt yVal: " << errorNorm << std::endl;
              Vector<Double> resX = res;
              resX /= mu;
              std::cout << "Remaining error-Norm wrt xVal: " << resX.NormL2() << std::endl;
            }
            
            break;
          } else {
            if(itCnt >= maxIter){
              std::cout << "-------------------------------------------------------" << std::endl;
              std::cout << "Inversion could not find a solution after " << itCnt << " iterations" << std::endl;
              //std::cout << "Last tol for convergence: " << tolForCheck << std::endl;
              if(wrtX){
                std::cout << "Remaining error-Norm wrt xVal: " << errorNorm << std::endl;
                Vector<Double> resY = res;
                resY *= mu;
                std::cout << "Remaining error-Norm wrt yVal: " << resY.NormL2() << std::endl;
              } else {
                std::cout << "Remaining error-Norm wrt yVal: " << errorNorm << std::endl;
                Vector<Double> resX = res;
                resX /= mu;
                std::cout << "Remaining error-Norm wrt xVal: " << resX.NormL2() << std::endl;
              }
              
              break;
            }
          }
          
          //        std::cout << "currentValue of x: " << xVal.ToString() << std::endl;
          //        std::cout << "residual: " << res.ToString() << std::endl;
          
          std::cout << "Perform Linesearch; starting error: " << errorNorm << std::endl;
          discardUpdate = performLinesearch(xVal, yVal, res, xUpdate, jac, jacT, mu, idElem, alpha, wrtX, relative);
          //discardUpdate = performLinesearch(xVal,yVal,mu,idElem,res,jac,jacT,alpha,xUpdate,wrtX,relative);
          
          if(!discardUpdate){
            xVal = xVal+xUpdate;
          } else {
            break;
          }
          
          sign = sign*(-1.0);
          
        }
      }
      actYval_[idElem] = yVal;
      actXval_[idElem] = xVal;
      return xVal;
      
    }
    
    void switchingStateToBmp(UInt numPixel, std::string filename, UInt idElem, bool overLayWithRotState = false){
      EXCEPTION("Not implemented in base class");
    }
    
    void setFlags(UInt performanceFlag, UInt textOutputFlag, UInt mappingFlag){
      
      if(performanceFlag >= 2){
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
      
      if(textOutputFlag > 2){
        textOutputFlag = 2;
      }
      textOutputLevel_ = textOutputFlag;
      
      mappingVersion_ = mappingFlag;
      
    }
    
    std::string runtimeToString(){
      EXCEPTION("Not implemented in base class");
    }
    
  protected:
    Vector<Double> evaluateNewRotationDirection(Vector<Double>& e_u_new, Vector<Double>& e_u_old, Double xVal);
    
    Vector<Double> clipNewRotationDirection(Vector<Double>& e_u_new);
    
    /*!
     * Global quantities, i.e. the same for all FE elements of the same material
     */
    Double XSaturated_; //! saturation value for  input
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
    
    bool classical_; //! switch between classical evaluation (2012 version of vector model) or revised evaluation (2015 version)
    
    
    /*!
     * Local quantities, i.e. one for each FE elements of the same material
     */
    Vector<Double>* preisachSum_;
    Vector<Double>* preisachSumTmp_;
    
    /*
     * for inversion 
     */
    Vector<Double>* actXval_;
    Vector<Double>* actYval_;  
    
    /*!
     * Quantities needed by the revised approach (Sutor2015)
     */
    Double angularDistance_; //! angular distance between old rotation state and new input direction which remains for small fields
    
    /*
     * The input direction gets clipped to discrete angles inside the used coordinate system;
     * angularClipping = 0 -> no clipping
     * angularClipping > 0 -> minimal angular step that gets resolved; value in rad
     */
    Double angularClipping_;
    
    
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
    
  };
  
  class VectorPreisachv10_MatrixApproach : public VectorPreisachv10
  {
  public:
    VectorPreisachv10_MatrixApproach(Integer numElem, Double xSat, Double ySat,
      Matrix<Double>& preisachWeight, Double rotationalResistance , UInt dim, bool isVirgin,
      bool classical, Double angularDistance, Double angularClipping);
    
    ~VectorPreisachv10_MatrixApproach();
    
    //! this function gets called from outside and calculates the output of the Preisach operator
    Vector<Double> computeValue_vec(Vector<Double>& xVal, Integer idElem, bool overwrite = true,bool overwriteDirection=true);
    
    void switchingStateToBmp(UInt numPixel, std::string filename, UInt idElem, bool overLayWithRotState = false);
    
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
  
  class VectorPreisachv10_ListApproach : public VectorPreisachv10
  {
    
  public:
    VectorPreisachv10_ListApproach(Integer numElem, Double xSat, Double ySat,
      Matrix<Double>& preisachWeight, Double rotationalResistance , UInt dim, bool isVirgin,
      bool classical, Double angularDistance, Double angularClipping);
    
    virtual ~VectorPreisachv10_ListApproach();
    
    //! this function gets called from outside and calculates the output of the Preisach operator
    Vector<Double> computeValue_vec(Vector<Double>& xVal, Integer idElem, bool overwrite = true,bool overwriteDirection=true);
    
    void switchingStateToBmp(UInt numPixel, std::string filename, UInt idElem, bool overLayWithRotState = false);
    
    std::string runtimeToString();
    
  private:
    
    /*
     * for version 10 -> revised model
     */
    void Update_GlobalRotationList(Double xThres, Double xVal, Vector<Double> e_u, std::list<RotListEntryv10>& usedRotationList,bool overwriteDirection=true);
    UInt Update_SwitchingList(std::list<ListEntryv10>& list, Double newEntry, Double lastXpar, Rectangle boundingBox, bool wasWipedOut, bool lastRotEntry);
    Double clipRectangleToElement(Rectangle& source, UInt idAlpha, UInt idBeta, Double delta = -1,bool isRotState = false);
    void getBoundingBoxFromRotEntry(std::list<RotListEntryv10>::iterator rotListIt, Rectangle& rect1, bool lastRotListEntryv10);
    bool getRectanglesFromRotEntry(std::list<RotListEntryv10>::iterator rotListIt, Rectangle& rect1, Rectangle& rect2, bool lastRotListEntryv10);
    void Evaluate_GlobalRotationList(std::list<RotListEntryv10>& usedRotationList, Vector<Double>& retVec);
    void Evaluate_LowerTriangle();
    Double mapRectangleToPreisachWeightsOLD(Rectangle& rect, bool skipUpperDiagonal = false);
    Double mapRectangleToPreisachWeightsNEW(Rectangle& rect, bool skipUpperDiagonal = false);
    Double mapRectangleToPreisachWeights(Rectangle& rect, bool skipUpperDiagonal = false);
    Double getRectangleFromSwitchingList(std::list<ListEntryv10>& list,
    std::list<ListEntryv10>::iterator startIt, std::list<ListEntryv10>::iterator curIt, std::list<ListEntryv10>::iterator endIt,
    UInt idArea, Rectangle& rect, bool upperSplitSquare = false);
    void Simplify_LocalSwitchingLists(std::list<RotListEntryv10>& usedRotationList);
    void Simplify_GlobalRotationList(std::list<RotListEntryv10>& usedRotationList);
    void mapRectangleToHelperMatrix(Matrix<Double>& helper, Rectangle rect, Double factor, bool skipUpperDiagonal = false,bool isRotState = false);
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

