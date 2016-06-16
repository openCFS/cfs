#include "VectorPreisach.hh"

#include <iostream>
#include <fstream>

/* See VectorPreisach.hh for detailed description */

namespace CoupledField
{ 

  VectorPreisach::VectorPreisach(Integer numElem, Double xSat, Double ySat,
		       Matrix<Double>& preisachWeight, Double rotationalResistance, UInt dim,  bool isVirgin, bool isTesting, UInt evalVersion)
            : Hysteresis(numElem)
          //: Preisach(numElem, xSat, ySat, preisachWeight, isVirgin)
  {

    /*
     * Global quantities, i.e. the same for all FE elements of the same material
     */
    if (xSat > 0 ) {
      Xsaturated_  = xSat;
    }
    else {
      Xsaturated_  = 1.0;
    }

    YSaturated_  = ySat;

    isVirgin_    = isVirgin;

    numElem_ = numElem;
    dim_ = dim;

    preisachWeights_ = preisachWeight;

    // get size of preisachWeights_
    UInt M = preisachWeights_.GetNumRows();
    UInt N = preisachWeights_.GetNumCols();

    if(M != N){
      EXCEPTION("Matrix preisachWeight has dim " << M << " x " << N << " and thus is not symmetric!");
    }

    numRows_ = M;

    // fixed for 3,4,5,6
//    if(numRows_%2 == 1){
//      /*
//       * odd number of rows -> currently the element with center 0,0 is treated wrong as it is not treated as
//       * initially split along alpha = -beta
//       */
//      std::cout << "Preisach plane has an odd number of rows (M = " << numRows_ << ")." << std::endl
//      << "In that case the Preisach element with center (0,0) is not treated correctly right now." << std::endl
//      << "Consider to use an even number of rows instead." << std::endl;
//    }


    // in the optimized version, we do not have to reevaluate the switching state of the lower triangle (alpha <= 0)
    // as it always is +1
    // -> to caputre the value of the lower triangle, we just have to sum up all weights once
    // in each calculation step, the rotation state of the lower triangle is calculated by multiplying the stored value
    // with the current vector e_u
    lowerTriangleValue_ = 0.0;

    for(UInt i = 0; i < (UInt) (numRows_/2.0); i++){
      for(UInt j = 0; j <= i; j++){
        if(i == j){
          // elements on the diagonal alpha = beta count only as 0.5
          lowerTriangleValue_ += 0.5*preisachWeights_[i][j];
        } else {
          lowerTriangleValue_ += preisachWeights_[i][j];
        }
      }
    }

    // resolution of Preisach plane
    delta_ = 2.0/Double(numRows_);

    tol_ = 1e-16;

    rotationalResistance_ = rotationalResistance;

    /*
     * Debugging quantities
     */
    isTesting_ = isTesting; //!true;

    if(isTesting_ && (evalVersion_ >= 5)){
      std::cout << "Testcase is for eval version 5 and 6 may not match the one of the other versions." << std::endl;
      std::cout << "Reason: The testcase allows the lower triangle part of the rotation plane to oppose the input vector."
                << "Therewith, the switching states in this lower triangle are not longer fix at +1, which is assumed in version 5 and 6."
                << std::endl;
    }

    rotStateTesting_ = Vector<Double>(dim_);
    rotStateTesting_.Init();
    rotStateTesting_[1] = 1;

    evalVersion_ = evalVersion; //0 = classical; 1 = weighted rotations; 2 = weighed rotation compact; 3 = weighted rotation compact, using helper functions
                      //4 = combined evalutation of switches and rotation states
                      //5 = optimized verison of version 4; 6 = optimized version of version 3

//    std::cout << "IsTesting: " << isTesting_ << std::endl;
//    std::cout << "evalVersoin: " << evalVersion_ << std::endl;

    /*
     * Local quantities, i.e. arrays storing different values for each FE element
     */
    preisachSum_ = new Vector<Double>[numElem];

    oldXpar_ = new Matrix<Double>[numElem];
    switchingStates_ = new Matrix<Double>[numElem];
    rotationStates_ = new Vector<Double>**[numElem];
    evaluatedState_ = new Vector<Double>**[numElem];

    wipedOut_ = new bool*[numElem];

    preMinMaxChanged_ = new bool**[numElem];
    preRotChanged_ = new bool**[numElem];

    preMinMaxValues_ = new std::list<ListEntryOld>**[numElem];
    preRotValues_ = new std::list<ListEntryOld>**[numElem];

    /*
     * Initialize arrays/vectors/matrices for each element
     */
    for(UInt k = 0; k < (UInt) numElem; k++){

      preisachSum_[k] = Vector<Double>(dim);
      preisachSum_[k].Init(0.0);

      oldXpar_[k] = Matrix<Double>(numRows_,numRows_);
      oldXpar_[k].Init();

      switchingStates_[k] = Matrix<Double>(numRows_,numRows_);
      switchingStates_[k].Init();

      /*
       * For output via bmp: set lower triangle ov switching states to +1
       * (only for version 5 and 6 which precalculate the lower triangle)
       */
      if((evalVersion_ == 5)||(evalVersion_ == 6)){
        for(UInt i = 0; i < numRows_; i++){
          for(UInt j = 0; j <= i; j++){
            if(i+j+1 < numRows_){
              // below diagonal alpha = -beta -> +1
              switchingStates_[k][i][j] = +1;
            }
            if (i == j){
              // on alpha = beta divide value by 2
              switchingStates_[k][i][j] = switchingStates_[k][i][j]/2.0;;
            }
          }
        }
      }

      if(evalVersion_ == 0){
        /*
         * for the classical approach, we must initialize the switching states in the classical way
         * i.e. split Preisach plane along diagonal alpha=-beta in a +1 and -1 part
         * Note that for the list based versions (evalVersion > 0) the initialization is done
         * during the first evaluation step
         */
        for(UInt i = 0; i < numRows_; i++){
          for(UInt j = 0; j <= i; j++){
            if(i+j+1 == numRows_){
              // on diagonal alpha = -beta -> leave at 0
            } else if(i+j+1 < numRows_){
              // below diagonal alpha = -beta -> +1
              switchingStates_[k][i][j] = +1;
            } else {
              // above diagonal alpha = -beta -> -1
              switchingStates_[k][i][j] = -1;
            }
            if (i == j){
              // on alpha = beta divide value by 2
              switchingStates_[k][i][j] = switchingStates_[k][i][j]/2.0;;
            }
          }
        }
      }

      /*
       * we need one entry for each element with center points lying on alpha = -beta
       * -> numDiag = half (numRows_ = even)
       * -> numDiag = half + 1 (numRows_ = odd)
       *
       * Remark: use the BETA coordinate to address this boolean as we go from beta = -1(-> idx=0) to beta = 0(->idx=half-1)
       */
      UInt numDiag = (UInt) numRows_/2.0 + numRows_%2;

      wipedOut_[k] = new bool[numDiag];
      for (UInt i = 0; i < numDiag; i++){
        wipedOut_[k][i] = false;
      }

      if(numRows_%2 == 1){
         /*
          * TODO: currently we do not have a special treatment for the element lying on both alpha = beta and alpha = -beta
          *     (element is split into 4 triangles, two used, two not, ...)
          *     If this treatment has been implemented below, set flag to false.
          */
        wipedOut_[k][numDiag-1] = !true;

      }

      rotationStates_[k] = new Vector<Double>*[numRows_];
      for (UInt i = 0; i < numRows_; i++){
        rotationStates_[k][i] = new Vector<Double>[numRows_];

        for (UInt j = 0; j <= i; j++){
          rotationStates_[k][i][j] = Vector<Double>(dim_);
          rotationStates_[k][i][j].Init();
          if(isTesting_ == true){
            rotationStates_[k][i][j][1] = 1.0;
          }
        }
      }

      evaluatedState_[k] = new Vector<Double>*[numRows_];
      for (UInt i = 0; i < numRows_; i++){
        evaluatedState_[k][i] = new Vector<Double>[numRows_];

        for (UInt j = 0; j <= i; j++){
          evaluatedState_[k][i][j] = Vector<Double>(dim_);
          evaluatedState_[k][i][j].Init();
        }
      }

      preMinMaxChanged_[k] = new bool*[numRows_];
      for (UInt i = 0; i < numRows_; i++){
        preMinMaxChanged_[k][i] = new bool[numRows_];

        for (UInt j = 0; j <= i; j++){
          preMinMaxChanged_[k][i][j] = false;
        }
      }

      preRotChanged_[k] = new bool*[numRows_];
      for (UInt i = 0; i < numRows_; i++){
        preRotChanged_[k][i] = new bool[numRows_];

        for (UInt j = 0; j <= i; j++){
          preRotChanged_[k][i][j] = false;
        }
      }

      Double alpha,beta;
      alpha = -1;
      preMinMaxValues_[k] = new std::list<ListEntryOld>*[numRows_];
      for (UInt i = 0; i < numRows_; i++){

        beta = -1;
        preMinMaxValues_[k][i] = new std::list<ListEntryOld>[numRows_];
        for (UInt j = 0; j <= i; j++){

          preMinMaxValues_[k][i][j] = std::list<ListEntryOld>();
          Initialize_MinMaxList(preMinMaxValues_[k][i][j],alpha,beta,k,i,j);

          beta = beta+delta_;
        }
        alpha = alpha+delta_;
      }

      alpha = -1;
      preRotValues_[k] = new std::list<ListEntryOld>*[numRows_];
      for (UInt i = 0; i < numRows_; i++){

        beta = -1;
        preRotValues_[k][i] = new std::list<ListEntryOld>[numRows_];
        for (UInt j = 0; j <= i; j++){

          preRotValues_[k][i][j] = std::list<ListEntryOld>();
          Initialize_RotList(preRotValues_[k][i][j],alpha,beta,k,i,j);

          beta = beta+delta_;
        }
        alpha = alpha+delta_;
      }
    }

  } // end constructor

  VectorPreisach::~VectorPreisach()
  {
    /*
     * too bad that vector<vector<bool> > is not recommended
     * and vector<vector<vector<list> > > did not work ...
     */
    for(UInt k = 0; k < numElem_; k++){
      for (UInt i = 0; i < numRows_; i++){
        delete[] preMinMaxValues_[k][i];
        delete[] preRotValues_[k][i];
        delete[] preMinMaxChanged_[k][i];
        delete[] preRotChanged_[k][i];
        delete[] rotationStates_[k][i];
        delete[] evaluatedState_[k][i];
      }
      delete[] preMinMaxValues_[k];
      delete[] preRotValues_[k];
      delete[] preMinMaxChanged_[k];
      delete[] preRotChanged_[k];
      delete[] wipedOut_[k];
      delete[] rotationStates_[k];
      delete[] evaluatedState_[k];
    }
    delete[] wipedOut_;
    delete[] preMinMaxValues_;
    delete[] preRotValues_;
    delete[] preMinMaxChanged_;
    delete[] preRotChanged_;
    delete[] preisachSum_;
    delete[] oldXpar_;
    delete[] rotationStates_;
    delete[] switchingStates_;
    delete[] evaluatedState_;
  }

  void VectorPreisach::Initialize_MinMaxList(std::list<ListEntryOld>& list, Double alpha, Double beta, UInt idElem, UInt idAlpha, UInt idBeta){
    if(isVirgin_){
      /*
       * Initialize the MinMaxList such that the evaluated switching states result to be:
       *
       * numCols_ = numRows_ = 4            numCols_ = numRows_ = 3
       *  _ _ _ _                           _ _ _
       * |0 - - -|                         |0 - -|
       * |+ 0 -  |                         |+ 0  |
       * |+ +    |                         |+_ _ |
       * |+_ _ _ |
       *
       * note that only for an odd number of rows, we have an element with value 0 on
       * the diagonal alpha = beta!
       *
       * Notation/Convention:
       *              alpha+delta,beta+delta
       *   _____________.
       *  |             |         Switching state = 1  <--> Max with value alpha+delta
       *  |             |         Switching state = -1 <--> Min with value beta
       *  |             |         Switching state = 0  <--> Max with value alpha or Min with value beta+delta (insert both)
       *  |             |
       *  |             |
       *  ._____________|
       *  alpha,beta
       *
       *              alpha+delta,beta+delta
       *   _____________.
       *  |           / |         Switching state = 0.5  <--> Max with value alpha+delta
       *  |         /   |         Switching state = -0.5 <--> Min with value beta
       *  |       /     |         Switching state = 0  <--> Max with value alpha or Min with value beta+delta (insert both)
       *  |     /       |
       *  |   /         |
       *  ._/___________|
       *  alpha,beta
       *
       * –> same treatment for triangles and squares
       *
       */

      /*
       * MinMaxList needs only the value itself and no vector information
       */
      Vector<Double> dummy = Vector<Double>(dim_);
      dummy.Init();

      /*
       * Make sure list is empty
       */
      list.clear();

      if(idAlpha+idBeta+1 == numRows_){
        /*
         * Element on diagonal alpha = -beta -> initial value = 0
         */
        /*
         * NOTE: the evaluation of this element will return 0 as long as
         * the maximum is alpha and the minimum is beta+delta
         * When one of these entries gets overwritten, the other one needs to be
         * deleted, too. Otherwise, the evaluation of will still return 0 although this is not wanted.
         * For the case, that the first input is a maximum, the following minimum will be deleted automatically
         * and thus the evaluation will work.
         * For the case, that the first input is a minimum, the maximum will not be deleted and thus we get 0 as
         * result.
         *
         * SOLUTION: set flag isDummy to true and check for this flag during updateprocess. If dummy element is detected,
         * clear all elements starting with the dummy element.
         */
        /*
         * helper max
         */
        list.push_back(ListEntryOld(alpha,dummy,false,true));
        /*
         * helper min
         */
        list.push_back(ListEntryOld(beta+delta_,dummy,true,true));
      } else if(idAlpha+idBeta+1 < numRows_){
        /*
         * Element below diagonal alpha = -beta -> initial value = 1
         */
        list.push_back(ListEntryOld(alpha+delta_,dummy,false));
      } else {
        /*
         * Element above diagonal alpha = -beta -> initial value = -1
         */
        list.push_back(ListEntryOld(beta,dummy,true));
      }

      preMinMaxChanged_[idElem][idAlpha][idBeta] = true;

    } else {
      EXCEPTION("Currently hysteresis is only supported starting from virgin state")
    }
  }

  void VectorPreisach::Initialize_RotList(std::list<ListEntryOld>& list, Double alpha, Double beta, UInt idElem, UInt idAlpha, UInt idBeta){
    if(isVirgin_){
      /*
       * Initialize the RotList with initial 0 vector and depending on the positon with
       *  1. value = alpha -> for all elements below the diagonal alpha = -beta
       *  2. value = beta+delta -> for all elements above the diagonal alpha = -beta
       *  3. two entries, one with value = alpha the other with value = beta+delta for all elements on alpha = -beta
       *
       * Idea behind this initialization:
       *     value = alpha is the minimal maximum; value = beta+delta the maximal minimum in each element, i.e. it will
       *     be replaced by any minimum/maximum entering the element during later update processes; furthermore, these
       *     values will lead to weighting areas of size 0, so that they will not contribute to the final Preisach sum
       *     unless the element gets overwritten by an input
       *
       */

      /*
       * MinMaxList needs only the value itself and no vector information
       */
      Vector<Double> zeroVec = Vector<Double>(dim_);
      zeroVec.Init();

      if(isTesting_){
        /*
         * Make sure list is empty
         */
        list.clear();
        /*
         * Initialize in with testing vector instead
         */
        list.push_back(ListEntryOld(alpha+delta_,rotStateTesting_,false));
        preRotChanged_[idElem][idAlpha][idBeta] = true;
        return;
      }

      /*
       * Make sure list is empty
       */
      list.clear();

      if(idAlpha+idBeta+1 == numRows_){
        /*
         * Element on diagonal alpha = -beta -> add min and max
         * -> does not result in symmetric behavior!
         * -> minima with maximal value and maxima with minimal value may
         *    not remain in the list if at least one other element is appended;
         *    otherwise, the calculation of the area returns 0
         *
         * Solution approach:
         *
         */
        /*
         * helper max
         */
        list.push_back(ListEntryOld(alpha,zeroVec,false));
        /*
         * helper min
         */
        list.push_back(ListEntryOld(beta+delta_,zeroVec,true));
      } else if(idAlpha+idBeta+1 < numRows_){
        /*
         * Element below diagonal alpha = -beta -> add minimal maximum
         */
        list.push_back(ListEntryOld(alpha,zeroVec,false));
      } else {
        /*
         * Element above diagonal alpha = -beta -> add maximal minimum
         */
        list.push_back(ListEntryOld(beta+delta_,zeroVec,true));
      }

      preRotChanged_[idElem][idAlpha][idBeta] = true;

    } else {
      EXCEPTION("Currently hysteresis is only supported starting from virgin state")
    }
  }

  void VectorPreisach::outputList(std::list<ListEntryOld>& list,bool isRot){

    std::list<ListEntryOld>::iterator listIt;

    if(isRot){
      std::cout << "Rotation-List contains: " << std::endl;
    } else {
      std::cout << "Switching-List contains: " << std::endl;
    }
    UInt cnt = 1;
    for(listIt = list.begin(); listIt != list.end(); listIt++){
        std::cout << "Entry " << cnt << ":"  << std::endl;
        std::cout << listIt->ToString();
        cnt++;
    }
  }

  void VectorPreisach::UpdateList(std::list<ListEntryOld>& list, Double newEntry, Vector<Double> newVecEntry, Double alpha, Double beta, UInt idElem, UInt idAlpha, UInt idBeta, bool isRot){

    /*
     * Update either MinMaxList (isRot == false) or RotList (isRot == true) with set of newEntry,newVecEntry
     *
     * Major differences:
     *  1. MinMaxList is only interested in newEntry; newVecEntry is not relevant (at the moment)
     *  2. MinMaxList has to be alternating list of minima and maxima;
     *    ->new entries are only added if
     *    i) larger maximum is added -> delete all smaller or equal maxima and all minima in between
     *    ii) smaller minimum is added -> delete all larger or equal minima and all maxima in between
     *    iii) last entry is maximum and minimum (or vice versa) and the value does not fulfill i) or ii)
     *    -> no new entries are added, if
     *    iv) last entry is a minimum or maximum and the input value does not fulfill i) or ii), i.e. a larger min
     *       may not follow a smaller min directly without a maximum in between
     *    -> special case: for elements on alpha = -beta not all minima and maxima following a deleted entry are also deleted
     *                      -> see below
     *
     *  3. RotList has not to be an alternating list of minima and maxima and the values are set differently
     *    i) all elements with alpha <= 0 are completely set to the new value, i.e. the lists have only one entry which is updated
     *    ii) all elements with beta >= 0 are completely set to the new value, i.e. the lists have only one entry which is updated
     *    iii) all elements with alpha = -beta are updated with two values each time, a minimum and a maximum
     *      -> Convention: add maximum first
     *    iv) all elements below the diagonal alpha = -beta are set only by a series of consecutive maxima and elements above this diagonal
     *      are only set by a series of consecutive minima
     *      -> larger maxima delete smaller maxima, smaller minima delete larger minima, smaller maxima or larger minima are appended to list
     *         without an alternating extremum in between
     *
     * General notes:
     *  For the later computation, we require the entries in the list to be in the range alpha, alpha+delta (for max) and beta, beta+delta (for min)
     *  If the input value is out of these bounds -> cut it and isert the cutted value
     *
     *  Input:
     *    MinMaxList: newEntry = xPar resulting from current rotation state of element; newVecEntry = not used
     *    RotList: newEntry = xThresh, newVecEntry = e_u
     *
     */

    std::list<ListEntryOld>::iterator listIt;
    /*
     * notes
     *  1. list.end() points to the element behind the last entry in the list
     *  2. after initialization, list is not entry, so decreasing lastEntry is valid
     *
     */
    if(list.empty()==true){
      EXCEPTION("List is not allowed to be empty at this point!")
    }

    std::list<ListEntryOld>::iterator lastEntry = --(list.end());
    Double newEntryUncut = newEntry;

    int state;
    /*
     * check if the current input is a minimum or a maximum
     */
    if(isRot){
      /*
       * the list of rotation operators is updated by different rules as stated above
       */
      if(idAlpha+idBeta+1 == numRows_){
        /*
         * element on diagonal alpha = -beta -> insert pair of max and min -> state = 0
         */
        state = 0;
      } else if(idAlpha+idBeta+1 < numRows_){
        /*
         * element below diagonal alpha = beta -> value is determined by maxima only
         */
        state = 1;
      } else {
        state = -1;
        /*
         * newEntry is xThreshold in this case, which is always >= 0
         * -> according to Sutor et al. we have to set the region above the diagonal with -xThreshold
         * -> change sign
         */
        newEntry = -newEntry;
      }
    } else {
      /*
       * get last value to determine a change in input
       * Note that we cannot use the value of the last entry as this value has been restricted to the
       * element boundaries
       */
      Double oldEntry = oldXpar_[idElem][idAlpha][idBeta];

      if(abs(newEntry-oldEntry)<tol_){
        /*
         * Value stayed the (nearly) the same
         */
        return;
      } else if(newEntry > oldEntry){
        /*
         * Value got larger -> max
         */
        state = 2;
      } else {
        /*
         * Value got smaller -> min
         */
        state = -2;
      }
    }

//    if(!isRot){
//      std::cout << "Update" << std::endl;
//      std::cout << "idAlpha,idBeta: " << idAlpha << "," << idBeta << std::endl;
//      std::cout << "preRotChanged_[idElem][idAlpha][idBeta]: " << preRotChanged_[idElem][idAlpha][idBeta] << std::endl;
//      std::cout << "state: " << state << std::endl;
//      std::cout << "newEntry: " << newEntry << std::endl;
//      std::cout << "alpha, alpha+delta: " << alpha << "," << alpha+delta_ << std::endl;
//      std::cout << "beta, beta+delta: " << beta << "," << beta+delta_ << std::endl;
//    }

    if(state >= 0){
      /*
       * Element is to be set by maximum -> restrict to [alpha,alpha+delta_] (state = 1,2)
       * or
       * Element on diagonal alpha = -beta (state = 0)
       * -> restrict newEntry to range [alpha,alpha+delta_] / or -newValue to [beta,beta+delta_]
       *
       */
      if(newEntry >= alpha+delta_){
        /*
         * current element will be set completely to +1 (in case of MinMaxList), will be rotated into direction newVecEntry (in case of RotList)
         */
        /*
         * empty list
         */
        list.clear();
        /*
         * cut value
         */
        newEntry = alpha+delta_;

        if(state == 0){
          /*
           * rotation vector to be set for an element on alpha = -beta
           * -> add pair of min and max starting with max
           */
          list.push_back(ListEntryOld(newEntry,newVecEntry,false));
          list.push_back(ListEntryOld(-newEntry,newVecEntry,true));
        } else {
          /*
           * add maximum to list
           */
          list.push_back(ListEntryOld(newEntry,newVecEntry,false));
          /*
           * elements on alpha = -beta also loose their special state
           */
        }
        if( idAlpha+idBeta+1 == numRows_ ){
          wipedOut_[idElem][idBeta] = true;
        }

        if(isRot){
          preRotChanged_[idElem][idAlpha][idBeta] = true;
        } else {
          preMinMaxChanged_[idElem][idAlpha][idBeta] = true;
          oldXpar_[idElem][idAlpha][idBeta] = newEntryUncut;
        }

        return;
      } else if(newEntry <= alpha){
        /*
         * This value will have no influence on the element result -> leave function
         */
        //std::cout << "return!" << std::endl;
        return;
      }
    } else {
      /*
       * Element is to be set by a minimum -> restrict to [beta,beta+delta]
       */
      if(newEntry < beta){
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
        newEntry = beta;

        /*
         * add minimum to list
         */
        list.push_back(ListEntryOld(newEntry,newVecEntry,true));
        /*
         * elements on alpha = -beta also loose their special state
         */
        if( idAlpha+idBeta+1 == numRows_ ){
          wipedOut_[idElem][idBeta] = true;
        }

        if(isRot){
          preRotChanged_[idElem][idAlpha][idBeta] = true;
        } else {
          preMinMaxChanged_[idElem][idAlpha][idBeta] = true;
          oldXpar_[idElem][idAlpha][idBeta] = newEntryUncut;
        }
        return;
      } else if(newEntry >= beta+delta_){
        /*
         * This value will have no influence on the element result -> leave function
         */
        //std::cout << "return!" << std::endl;
        return;
      }
    }

    /*
     * if we came to this point, we have to iterate through the list and check if value shall be included
     */
    bool canBeInserted = false;
    bool canBeDeleted = false;
    bool insertHelper = false;

    if(isRot){
      /*
       * the rotation list is always updated if newEntry lies inside the element boundaries, i.e. larger minima can directly follow smaller ones
       * and smaller maxima can directly follow larger maxima
       */
      canBeInserted = true;
    } else {
      if( (idAlpha+idBeta+1 == numRows_)&&(wipedOut_[idElem][idBeta] == false) ){
        /*
         * for non-wiped out elements on the diagonal alpha = -beta we have to be little more careful
         *
         *  Case 1: xPar does not replace the first entry (x0) of the list but any other value later on
         *          -> std procedure
         *  Case 2: xPar replaces the first entry (x0) of the list; list is otherwise empty
         *          -> std procedure
         *  Case 3: xPar replaces the first entry, list is not empty (second element in list shall be x1)
         *          3a) xPar is min, xPar-beta <= alpha+delta-x1
         *              -> xPar replaces x0 and is so far to the left that it wipes out rest of the list -> std procedure
         *          3b) xPar is min, xPar-beta > alpha+delta-x1
         *              -> xPar is able to delete x0, but it cuts the diagonal below x1
         *              -> x1 will still have an influence on the result
         *
         *                   xPar                xPar
         *               ____v____           ____v____
         *              |\   !  ¦ |         |\   !    |
         *              |__\_!__¦_|<x1      |  \_!    |<x1 -> not wiped out!
         *              |    \  ¦ |      -> |   ^!    |    -> note the step! it cannot recreated without x1
         *              |    ! \¦ |         |    !    |
         *              |____!__¦\|         |____!____|
         *                     x0
         *
         *              -> save x1, delete list, add x1 as max, add xPar as min
         *
         *          3c) xPar is max, xPar-alpha >= beta+delta-x1 -> 3a)
         *          3d) xPar is max, xPar-alpha < beta+delta-x1 -> similar to 3b)
         *
         *          Note: for 3a it is really <= and for 3c >= -> checked
         *
         */
        if( (state == 2)&&(list.begin()->isMin() == false) ){
          /*
           * we want to insert a maximum and list starts with maximum -> take special care
           */
          state = 3;
        } else if( (state == -2)&&(list.begin()->isMin() == true) ){
          /*
           * we want to insert a minimum and list starts with minimum -> take special care
           */
          state = -3;
        } // else -> Case1 -> no special treatment necessary
      }

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
    }

    bool listMin;
    Double listVal;
    ListEntryOld helperEntry = ListEntryOld(0,Vector<Double>(dim_),false);

    for(listIt = list.begin(); listIt != list.end(); listIt++){

      /*
       * NEW: For elements on alpha = -beta, we start with a pair of min and max such that the resulting area = 0
       * However, if only the minimum is overwritten, the evaluation algorithm will return 0 due to the still remaining
       * initial maximum of value = alpha. The easiest solution is to remove the maximum, too.
       * To detect this case, check for the isDummy flag, which is only set for elements of the minmaxlist for alpha = -beta.
       *
       * NOTE: The rotation list does not have this issue, as all new values on alpha = -beta will be pairs of minima and maxima, too.
       */
      if(listIt->isDummy()){
          canBeInserted = true;
          canBeDeleted = true;
      }

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
        }
      } else if(state == 3){
        /*
         * special treatment -> see above
         */
        if(!listMin){
          /*
           * compare value with value of list maximum
           */
          if(newEntry >= listVal){
            /*
             * first entry of list is to be replaced
             * (Note: state == 3 can only occur on the first iteration as state is then set back to 2)
             */

            canBeInserted = true;
            canBeDeleted = true;
            /*
             * check if list is empty otherwise
             */
            if(listIt == lastEntry){
              /*
               * list only has one entry -> std procedure (Case 2 from aboves description)
               */
              insertHelper = false;
            } else {
              /*
               * check if value is large enough to wipe out the following value
               */
              listIt++;
              Double nextEntry = listIt->getVal();
              listIt--;

              if(newEntry-alpha >= beta+delta_-nextEntry){
                /*
                 * value is large enough -> std procedure (Case 3c)
                 */
                insertHelper = false;
              } else {
                /*
                 * value is not large enough; create a backup of next entry and insert it later
                 * as the new beginning of the list followed by the current entry (Case 3d)
                 */
                listIt++;
                helperEntry = *listIt;
                listIt--;
                insertHelper = true;
              }
            }
          } else {
            /*
             * Value is not to be inserted at first position -> change back to state = 2
             */
            state = 2;
          }
        }
      } else if(state == -3){
        /*
         * special treatment -> see above
         */
        if(listMin){
          /*
           * compare value with value of list minimum
           */
          if(newEntry <= listVal){
            /*
             * first entry of list is to be replaced
             * (Note: state == -3 can only occur on the first iteration as state is then set back to -2)
             */

            canBeInserted = true;
            canBeDeleted = true;
            /*
             * check if list is empty otherwise
             */
            if(listIt == lastEntry){
              /*
               * list only has one entry -> std procedure (Case 2 from aboves description)
               */
              insertHelper = false;
            } else {
              /*
               * check if value is large enough to wipe out the following value
               */
              listIt++;
              Double nextEntry = listIt->getVal();
              listIt--;

              if(newEntry-beta <= alpha+delta_-nextEntry){
                /*
                 * value is small enough -> std procedure (Case 3a)
                 */
                insertHelper = false;
              } else {
                /*
                 * value is not small enough; create a backup of next entry and insert it later
                 * as the new beginning of the list followed by the current entry (Case 3b)
                 */
                listIt++;
                helperEntry = *listIt;
                listIt--;
                insertHelper = true;
              }
            }
          } else {
            /*
             * Value is not to be inserted at first position -> change back to state = -2
             */
            state = -2;
          }
        }
      } else if(state == 1){
        /*
         * we iterate over RotList and are below the diagonal (alpha = -beta), i.e. we have to insert a max
         */
        if(!listMin){
          /*
           * compare with value of list maximum
           * Note: Even though these elements should only have max in their list, make sure that we really compare with one
           */
          if(newEntry >= listVal){
            canBeDeleted = true;
          }
        } else {
          EXCEPTION("A value below the diagonal alpha = -beta should not have minima in the rotList")
        }
      } else if(state == -1){
        /*
         * we iterate over RotList and are above the diagonal (alpha = -beta), i.e. we have to insert a min
         */
        if(listMin){
          if(newEntry <= listVal){
            canBeDeleted = true;
          }
        } else {
          EXCEPTION("A value above the diagonal alpha = -beta should not have maxima in the rotList")
        }
      } else if(state == 0){
        /*
         * we iterate over RotList and are on the diagonal (alpha = -beta), i.e. insert a max and a min (max first)
         */
        if(!listMin){
          if(newEntry >= listVal){
            canBeDeleted = true;
          } else {
            /*
             * we know that these elements always have pairs of min and max with same absolute value of the scalar entry
             * -> jump over minimum
             * Important: increase counter only if list is not to be deleted -> deletion stats at listIt!
             */
            listIt++;
          }
        }
      }

      if(canBeDeleted){
        /*
         * delete all entries of the list starting with the current one
         */
        list.erase(listIt,list.end());
        break;
      }
    }

    if(canBeInserted == true){
      if(insertHelper == true){
        list.push_back(helperEntry);
      }
      if(state == 0){
        /*
         * add pair of min and max
         */
        list.push_back(ListEntryOld(newEntry,newVecEntry,false));
        list.push_back(ListEntryOld(-newEntry,newVecEntry,true));
      } else if(state > 0){
        /*
         * insert maximum
         */
        list.push_back(ListEntryOld(newEntry,newVecEntry,false));
      } else if(state < 0){
        /*
         * insert minimum
         */
        list.push_back(ListEntryOld(newEntry,newVecEntry,true));
      }
      /*
       * mark element as updated
       *
       * Note: the flag will not be set to false in this function but at the end of the evaluation process
       *      by this, we assure that the evaluation process is triggered at least once
       */
      if(isRot){
        preRotChanged_[idElem][idAlpha][idBeta] = true;
      } else {
        preMinMaxChanged_[idElem][idAlpha][idBeta] = true;
        oldXpar_[idElem][idAlpha][idBeta] = newEntryUncut;
      }
    }
  }


  void VectorPreisach::EvaluatePreisachElement(std::list<ListEntryOld>& list, Double& retVal, Vector<Double>& retVec, Double alpha, Double beta, UInt idElem, UInt idAlpha, UInt idBeta, bool isRot){

    Double val = 0;
    Vector<Double> vec = Vector<Double>(dim_);
    vec.Init();

    if((isRot)&&(preRotChanged_[idElem][idAlpha][idBeta] == false)){
      /*
       * no change this last time -> return previous value
       */
      retVal = val;
      retVec = rotationStates_[idElem][idAlpha][idBeta];
      return;
    } else if((!isRot)&&(preMinMaxChanged_[idElem][idAlpha][idBeta] == false)){
      /*
       * no change this last time -> return previous value
       */
      retVal = switchingStates_[idElem][idAlpha][idBeta];
      retVec = vec;
      return;
    }

    if(list.empty()){
      EXCEPTION("List should not be empty at this point!");
    }

    Double update = 0.0;
    UInt cnt = 0;
    Double top = 0.0;
    Double bot = 0.0;
    Double left = 0.0;
    Double right = 0.0;

    bool isMinEntry;
    Double entry;
    Vector<Double> entryVec;

    std::list<ListEntryOld>::iterator listIt;
    /*
     * list is not empty, so we can get iterator to last element
     * note that .end() is pointing to the element past the last element
     */
    std::list<ListEntryOld>::iterator listEnd = --(list.end());

    /*
     * check if we have a full (i.e. square) element or if we are on the diagonal alpha = beta
     */
    if(idAlpha != idBeta){
      /*
       * square
       *
       */
      /*
       *  Switching operators (-> list is alternating sequence of min and max)
       *
       *
       *    element from (alpha,beta) to (alpha+delta,beta+delta)
       *    example input min max list (starting with max):
       *      max: alpha+3*delta/4
       *      min: beta+delta/2
       *      max: alpha+delta/2
       *      min: beta+delta/4
       *      max: alpha+delta/4
       *
       *            MinMax[1] = min
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
       *    resulting areas:
       *
       *       quadratic part which is still split by diagonal -> only for elements on alpha = -beta
       *       (note that upper and lower part cancel out in sum, so that we can omit this part)
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
       *      TODO:
       *      OBSERVATION 10.5.16:
       *      MinMax[1] may be smaller than beta* ! -> area 1 is negative and overlaps area3!!!
       *
       *      RESULTING PROBLEMS:
       *      - As area1 will be negative (as left > right), it was not counted as area1 was supposed to be > tol
       *      - Area0 is no square anymore, i.e. its value does not cancel out anymore!
       *
       *     here is MinMax[1] = min
       *       __v______________
       *      | \|1|        |   |
       *area0>|__|\|     3  |   |
       *no    |2 |          |5  |
       *square|__|__________|   |
       *      |     4       |   |
       *      |_____________|___|
       *      |        6        |
       *      |_________________|
       *
       *      Possible solutions:
       *       1. subtract area1 from area 3, i.e. let the negative value of area1 count
       *
       *       Resulting state (close up):
       *       __v______________
       *      | \      |        |
       *      |   \    |   3-1  |
       *      |     \ _|        |
       *      |______|x|        |
       *      |     2|          |
       *      |      |          |
       *      |      |          |
       *      |______|__________|
       *
       *        -> area0 is a square again (area marked as x belongs to the square, i.e. it is half -1 and half +1)
       *        -> areax would belong to the negative part and is not considered to be 0!
       *        -> solution1 better than original version, but not correct!
       *
       *
       *      2. check if area1 would be negative
       *      if yes, define another area (area1*) instead, which is added to area2 (i.e. will be positive)
       *
       *       Resulting state (close up):
       *       __v______________
       *      | \    | |        |
       *      |   \  |3* 3      |
       *      |_____\|_|        |
       *      |_1*___|x|        |
       *      |     2|          |
       *      |      |          |
       *      |      |          |
       *      |______|__________|
       *
       *        -> area0 is a square again; the addition area 1* is compensated with area 3* (included in area3)
       *        -> areax belongs to area3
       *        -> solution should be approprite now!
       *
       *
       *
       *      (area x: width x height)
       *      area 1: (MinMax[1] - beta*) x (alpha+delta-MinMax[0])
       *        with  beta* = alpha+delta - MinMax[0]+beta for elements on alpha = -beta
       *              beta* = beta                         for elements not on alpha = -beta
       *      area 2: (MinMax[1] - beta) x (MinMax[0] - MinMax[2])
       *      area 3: (MinMax[3] - MinMax[1]) x (alpha+delta - MinMax[2])
       *      area 4: (MinMax[3] - beta) x (MinMax[2] - MinMax[4])
       *      area 5: (beta+delta - MinMax[3]) x (alpha+delta - MinMax[4])
       *      area 6: (beta+delta - beta) x (MinMax[4] - alpha)
       *
       *    element from (alpha,beta) to (alpha+delta,beta+delta)
       *    example input min max list (starting with min):
       *      min: beta+3*delta/4
       *      max: alpha+delta/2
       *      min: beta+delta/2
       *      min: beta+delta/4
       *      max: alpha+delta/4
       *
       *   MinMax[0] = min  MinMax[2] = min
       *       ____v________v___
       *      | \  |   |    |   |
       *      |   \|   |    |   |
       *      |    |\  |    |   |
       *      |____|__\|____|___| MinMax[1] = max
       *      |    |   |\   |   |
       *      |____|___|__\_|___| MinMax[3] = max
       *      |    |   |    \   |
       *      |____|___|____| \_|
       *              MinMax[4] = min
       *
       *    resulting areas:
       *
       *       quadratic part which is still split by diagonal -> only for elements on alpha = -beta
       *       (note that upper and lower part cancel out in sum, so that we can omit this part)
       *       _v_______________
       *      | \  | 2 |    |   |
       *      |___\|   | 4  |   |
       *      |   1|   |    |6  |
       *      |____|___|    |   |
       *      |     3  |    |   |
       *      |________|____|   |
       *      |        5    |   |
       *      |_____________|___|
       *
       *      (area x: width x height)
       *      area 1: (MinMax[0] - beta) x (alpha* - MinMax[1])
       *        with  alpha* = alpha+delta - MinMax[0] + beta for elements on alpha = -beta
       *              alpha* = alpha+delta_                   for elements not on alpha = -beta
       *      area 2: (MinMax[2] - MinMax[0]) x (alpha+delta - MinMax[1])
       *      area 3: (MinMax[2] - beta) x (MinMax[1] - MinMax[3])
       *      area 4: (MinMax[4] - MinMax[2]) x (alpha+delta - MinMax[3])
       *      area 5: (MinMax[4] - beta) x (MinMax[3] - alpha)
       *      area 6: (beta+delta - MinMax[4]) x (alpha+delta - alpha)
       *
       *  Rotation operators (-> list must not be alternating sequence of min and max)
       *    -> check two consecutive entries of list for their extremum type
       *    -> if two minima follow each other, assume maximum in between with value alpha
       *
       *      sequence of minima -> for formula set maxima in between to alpha
       *       _v_______________
       *      | 1  | 2 |  3 | 4 |
       *      |    |   |    |   |
       *      |    |   |    |   |
       *      |    |   |    |   |
       *      |    |   |    |   |
       *      |    |   |    |   |
       *      |    |   |    |   |
       *      |____|___|____|___|
       *
       *    -> if two maxima follow each other, assume minimum in between with value beta+delta_
       *
       *
       */
      for(listIt = list.begin(); listIt != list.end(); listIt++){

          entry = listIt->getVal();
          entryVec = listIt->getVec();
          isMinEntry = listIt->isMin();

          if(isMinEntry && (cnt == 0)){
            /*
             * First entry in list is a minimum; two possibilities:
             * 1. element is not wiped out yet -> minimum lies between beta and beta+delta
             * 2. element is wiped out by minimum or got initialized or overwritten with minimum -> minimum has value beta
             *
             * for case 2, area 1 will be 0 as left = right
             *
             * Note that the area 2 is already a std area, so instead of writing it once again into this case, we just
             * execute it additionally to this one here.
             */
            /*
             * area 1
             */
            left = beta;
            right = entry;
            top = alpha+delta_ - entry + beta;

            if(listIt != listEnd){
              listIt++;
              if(listIt->isMin() == isMinEntry){
                /*
                 * element of same type detected -> ok for rotList, not ok for MinMaxList
                 */
                if(isRot){
                  /*
                   * set bot to alpha
                   */
                  bot = alpha;
                } else {
                  EXCEPTION("MinMaxList has to be alternating!")
                }
              } else {
                bot = listIt->getVal();
              }
              listIt--;
            } else {
              bot = alpha;
            }

            /*
             * calculate area
             * Note: even if right = left; or top = bot, the value update may be != 0 due to numerics
             * dividng by values smaller delta with delta < 1 will increase this value, leading probably to detectable changes in the result
             * -> check if update is larger tol before dividng by delta*2
             */
            update = ((right-left)*(top-bot));
            if(update > tol_){
              update = update/(delta_*delta_);
            } else {
              update = 0.0;
            }
            /*
             * helper area is positive, so +
             */
            val += update;
            /*
             * resulting vector for the element will be the weighted sum of the stored, averaged rotation vectors
             * the weights will be the areas corresponding to those rotation vectors
             * -> only relevant for rotList
             *
             * NOTE: always add up; negative sign comes via switching operators!
             */
            vec += entryVec*update;

          } else if(!isMinEntry && (cnt == 0)) {
            /*
             * First entry in list is a maximum; two possibilities:
             * 1. element is not wiped out yet -> maximum lies between alpha and alpha+delta
             * 2. element is wiped out by maximum or got initialized or overwritten with maximum -> maximum has value alpha+delta
             *
             * for case 2, area 1 will be 0 as top = bot
             *
             * Note that the area 2 is already a std area, so instead of writing it once again into this case, we just
             * execute it additionally to this one here.
             */
            /*
             * area 1
             */
            if(listIt != listEnd){
              listIt++;
              if(listIt->isMin() == isMinEntry){
                /*
                 * element of same type detected -> ok for rotList, not ok for MinMaxList
                 */
                if(isRot){
                  /*
                   * set right to beta+delta_
                   */
                  right = beta+delta_;
                } else {
                  EXCEPTION("MinMaxList has to be alternating!")
                }
              } else {
                right = listIt->getVal();
              }
              listIt--;
            } else {
              right = beta+delta_;
            }

            left = alpha+delta_ - entry + beta;
            top = alpha+delta_;
            bot = entry;

            /*
             * calculate area
             * Note: even if right = left; or top = bot, the value update may be != 0 due to numerics
             * dividng by values smaller delta with delta < 1 will increase this value, leading probably to detectable changes in the result
             * -> check if update is larger tol before dividng by delta*2
             */
            update = ((right-left)*(top-bot));
            if(update > tol_){
              update = update/(delta_*delta_);
            } else {
              update = 0.0;
            }
            /*
             * helper area is negative, so -
             */
            val -= update;
            /*
             * NOTE: always add up; negative sign comes via switching operators!
             */
            vec += entryVec*update;
          }

          /*
           * The cases below are called independently of the counter cnt.
           * During the first iteration (cnt==0), it will compute area 2.
           */
          if(isMinEntry){
            /*
             * Current entry is negative
             */
            right = beta+delta_;
            bool found = false;

            if(listIt != listEnd){
              listIt++;
              if(listIt->isMin() == isMinEntry){
                /*
                 * element of same type detected -> ok for rotList, not ok for MinMaxList
                 */
                if(isRot){
                  /*
                   * set bottom to alpha
                   */
                  bot = alpha;
                  /*
                   * here we can directly set the next entry, too
                   * right is determined by the consecutive minimum, which just got found here
                   */
                  found = true;
                  right = listIt->getVal();
                } else {
                  EXCEPTION("MinMaxList has to be alternating!")
                }
              } else {
                bot = listIt->getVal();
              }

              if(listIt != listEnd){
                listIt++;
                if(listIt->isMin() != isMinEntry){
                  EXCEPTION("MinMaxList has to be alternating!")
                  /*
                   * here we would expect the same entry type the the one of the current iteration
                   * e.g. a minimum after a maximum which followd a minimum
                   */
                } else {
                  if(isRot&&found){
                    /*
                     * do nothing here -> we need the second min/max here but this value was already found in the previous entry
                     *
                     * in case of alternating rotList, we retrieved just the value of bot in the last case, so we need to find right
                     */
                  } else {
                    right = listIt->getVal();
                  }
                }
                listIt--;
              }
              listIt--;
            } else {
              bot = alpha;
            }

            top = alpha+delta_;
            left = entry;

            /*
             * calculate area
             * Note: even if right = left; or top = bot, the value update may be != 0 due to numerics
             * dividng by values smaller delta with delta < 1 will increase this value, leading probably to detectable changes in the result
             * -> check if update is larger tol before dividng by delta*2
             */
            update = ((right-left)*(top-bot));
            if(update > tol_){
              update = update/(delta_*delta_);
            } else {
              update = 0.0;
            }
            val -= update;
            /*
             * NOTE: always add up; negative sign comes via switching operators!
             */
            vec += entryVec*update;

          } else {
            /*
             * Current entry is positive
             */
            bot = alpha;
            bool found = false;

            if(listIt != listEnd){
              listIt++;
              if(listIt->isMin() == isMinEntry){
                /*
                 * element of same type detected -> ok for rotList, not ok for MinMaxList
                 */
                if(isRot){
                  /*
                   * set right to beta+delta_
                   */
                  right = beta+delta_;
                  /*
                   * here we can directly set the next entry, too
                   * right is determined by the consecutive minimum, which just got found here
                   */
                  found = true;
                  bot = listIt->getVal();
                } else {
                  EXCEPTION("MinMaxList has to be alternating!")
                }
              } else {
                right = listIt->getVal();
              }

              if(listIt != listEnd){
                listIt++;
                if(listIt->isMin() != isMinEntry){
                  EXCEPTION("MinMaxList has to be alternating!")
                  /*
                   * here we would expect the same entry type the the one of the current iteration
                   * e.g. a minimum after a maximum which followd a minimum
                   */
                } else {
                  if(isRot && found){
                    /*
                     * do nothing here -> we need the second min/max here but this value was already found in the previous entry
                     *
                     * in case of alternating rotList, we retrieved just the value of right in the last case, so we need to find bot
                     */
                  } else {
                    bot = listIt->getVal();
                  }
                }
                listIt--;
              }
              listIt--;
            } else {
              right =  beta+delta_;
            }

            top = entry;
            left = beta;

            /*
             * calculate area
             * Note: even if right = left; or top = bot, the value update may be != 0 due to numerics
             * dividng by values smaller delta with delta < 1 will increase this value, leading probably to detectable changes in the result
             * -> check if update is larger tol before dividng by delta*2
             */
            update = ((right-left)*(top-bot));
            if(update > tol_){
              update = update/(delta_*delta_);
            } else {
              update = 0.0;
            }
            val += update;
            /*
             * NOTE: always add up; negative sign comes via switching operators!
             */
            vec += entryVec*update;
          }
          cnt++;
      }

    } // end square
    else {
      /*
       * triangle due to diagonal alpha = beta
       *
       */
      /*
       * This case is more complicated as the one of the square, as we now have to cut the areas at the
       * diagonal alpha = beta.
       *
       * Proposal:
       *  1. apply same proceducre as for square elements but adapt areas such that at max the lower right corner
       *      lies below the diagonal alpha = beta i.e.
       *       - if the lower left corner would lie below the diagonal (area 4 and 6), restrict bottom bound to
       *            left bound (due to diagonal, alpha = beta)
       *       - if the upper right corner would lie below the diagonal (area 3 and 5), restrict the right bound
       *            to the top bound (again alpha = beta, so coordinates are symmetric w.r.t diagonal)
       *  2. subtract triangular areas which are beyond the diagonal alpha = beta
       *       this area will have the value: (right-bottom)^2/2
       *
       *      area0
       *       _v_______________
       *      | \  | 2   |  |6/ |
       *      |___\|     |4 /___|
       *      |   1|     |/_|   |
       *      |____|____/|      |
       *      |     3 /|        |
       *      |_____/__|        | <- apply std procedure, but extend areas at max till diagonal, then subtract
       *      |5  /|            |               triangles below diagonal
       *      |_/__|____________|
       *
       *
       * Discussion:
       *  Element lying on alpha = beta and alpha = -beta -> only one element, only for an odd number of rows
       *  of the Preisach plane. This element is split by two diagonals, thus requiring some extra thoughts.
       *  1. square around the split part (area0) will consist of a positive and
       *     a negative triangle of same size; these triangles will cancel out if summed, so we can omit their calculation.
       *     If these triangles overlap the diagonal alpha = beta, the overlapping parts would have to be subtracted.
       *     However, these overlapping parts cancel out, too. -> no extra treatment
       *
       *
       *     In the new evaluation version we insert a pair of min and max for all elements with alpha = -beta
       *     and on all elements on alpha = beta we insert a helper maximum which completely set the value to +1 or -1
       *     -> The points below are not true anymore!
       *
       *  2. initial value of elements on alpha = beta will be 0, but due the triangular splitting, they will have a
       *     helper max/min. This is true, also for the element alpha = beta = -beta and is different from other elements
       *     on alpha = -beta (those elements have no helper entry). The value of this helper min/max is chosen such,
       *     that elements on alpha = beta != -beta have an initial value of 0, when evaluated. This value does not hold
       *     for the case alpha = beta = -beta as we omit one part of the element (area0) which would belong completely
       *     to the positive or negative part resulting from the splitting by the helper min/max.
       *
       *      area0 -> would belong to -1 part here
       *       _v_______________
       *      | \    |        / |
       *      |   \  | area1/   |
       *      |_____\|____/_____|< helper max != alpha+delta -> area0 and area1 != 0
       *      |         /       |
       *      | area2 /         |  Note: area0 may only be omitted for alpha = beta = -beta!
       *      |     /           |        for alpha = beta != -beta, area1 has to be extended to include area0
       *      |   /             |
       *      |_/_______________|
       *
       *
       *     Solution approach: change helper min/max to beta+delta/2 / alpha+delta/2
       *      -> active area of triangle (the area above alpha = beta) is split into a square and two triangles of same
       *         size. The square will be neglected during the evaluation, as it is the area split by the diagonal
       *         alpha = -beta. The two triangles will have the same size but different signs, so they cancel out.
       *  3. as the helper min/max for alpha = beta != -beta will not be beta / alpha+delta, the algorithm of the square
       *     element case would omit one area0 although this is not wanted here.
       *      -> change beta*, alpha* from the square element algorithm
       *          alpha* = alpha+delta - MinMax[0] + beta -> alpha* = alpha + delta
       *          beta* = alpha+delta - MinMax[0] + beta -> beta* = beta
       *
       */

      for(listIt = list.begin(); listIt != list.end(); listIt++){

        entry = listIt->getVal();
        entryVec = listIt->getVec();
        isMinEntry = listIt->isMin();

        if(isMinEntry && (cnt == 0)){
          /*
           * First entry in list is a minimum; two possibilities:
           * 1. element is not wiped out yet -> minimum lies between beta and beta+delta
           * 2. element is wiped out by minimum or got initialized or overwritten with minimum -> minimum has value beta
           *
           * for case 2, area 1 will be 0 as left = right
           *
           * Note that the area 2 is already a std area, so instead of writing it once again into this case, we just
           * execute it additionally to this one here.
           */

          /*
           * area 1
           */
          left = beta;
          right = entry;

          // TODO:
          /*
           * still needed?
           *
           */
          if( (idAlpha+idBeta+1 == numRows_)&&(wipedOut_[idElem][idBeta] == false) ){
            /*
             * In contrast to the square element case, we have to distinguish between alpha = -beta elements
             * which are not wiped out and all other elements as alrea0 will not automatically get 0
             */
            top = alpha+delta_ - entry + beta;
          } else {
            top = alpha+delta_;
          }

          if(listIt != listEnd){
            listIt++;
            if(listIt->isMin() == isMinEntry){
              /*
               * element of same type detected -> ok for rotList, not ok for MinMaxList
               */
              if(isRot){
                /*
                 * set bot to alpha
                 */
                bot = alpha;
              } else {
                EXCEPTION("MinMaxList has to be alternating!")
              }
            } else {
              bot = listIt->getVal();
            }
            listIt--;
          } else {
            bot = alpha;
          }

          /*
           * restrict bound such, that only the bottom right corner lies below the diagonal
           */
          if(right-beta > top-alpha){
            /*
             * upper right corner below diagonal
             * -> set right bound to top
             */
            right = top;
          }
          if(bot-alpha < left-beta){
            /*
             * lower left corner below diagonal
             * -> set bottom bound to left
             */
            bot = left;
          }

          /*
           * calculate area
           * Note: even if right = left; or top = bot, the value update may be != 0 due to numerics
           * dividng by values smaller delta with delta < 1 will increase this value, leading probably to detectable changes in the result
           * -> check if update is larger tol before dividng by delta*2
           */
          update = ((right-left)*(top-bot));
          if(update > tol_){
            update = update/(delta_*delta_);
          } else {
            update = 0.0;
          }

          /*
           * subtract area below diagonal (if any)
           * -> check if lower right corner lies below diagonal
           */
          if(bot-alpha < right-beta){
            /*
             * subtract triangle with side length (right - bot)/delta
             */
            Double triang = 0.5*((right-bot)*(right-bot));
            if(triang > tol_){
              triang = triang/(delta_*delta_);
            } else {
              triang = 0.0;
            }
            update -= triang;
          }
          /*
           * helper area is positive, so +
           */
          val += update;

          /*
           * resulting vector for the element will be the weighted sum of the stored, averaged rotation vectors
           * the weights will be the areas corresponding to those rotation vectors
           * -> only relevant for rotList
           *
           * NOTE: always add up; negative sign comes via switching operators!
           */
          vec += entryVec*update;

        } else if(!isMinEntry && (cnt == 0)) {
          /*
           * First entry in list is a maximum; two possibilities:
           * 1. element is not wiped out yet -> maximum lies between alpha and alpha+delta
           * 2. element is wiped out by maximum or got initialized or overwritten with maximum -> maximum has value alpha+delta
           *
           * for case 2, area 1 will be 0 as top = bot
           *
           * Note that the area 2 is already a std area, so instead of writing it once again into this case, we just
           * execute it additionally to this one here.
           */
          /*
           * area 1
           */
          if(listIt != listEnd){
            listIt++;
            if(listIt->isMin() == isMinEntry){
              /*
               * element of same type detected -> ok for rotList, not ok for MinMaxList
               */
              if(isRot){
                /*
                 * set bottom to beta+delta_
                 */
                right = beta+delta_;
              } else {
                EXCEPTION("MinMaxList has to be alternating!")
              }
            } else {
              right = listIt->getVal();
            }
            listIt--;
          } else {
            right = beta+delta_;
          }

          if( (idAlpha+idBeta+1 == numRows_)&&(wipedOut_[idElem][idBeta] == false) ){
            /*
             * In contrast to the square element case, we have to distinguish between alpha = -beta elements
             * which are not wiped out and all other elements as area0 will not automatically get 0
             */
            left = alpha+delta_ - entry + beta;
          } else {
            left = beta;
          }

          top = alpha+delta_;
          bot = entry;

          /*
           * restrict bound such, that only the bottom right corner lies below the diagonal
           */
          if(right-beta > top-alpha){
            /*
             * upper right corner below diagonal
             * -> set right bound to top
             */
            right = top;
          }
          if(bot-alpha < left-beta){
            /*
             * lower left corner below diagonal
             * -> set bottom bound to left
             */
            bot = left;
          }

          /*
           * calculate area
           * Note: even if right = left; or top = bot, the value update may be != 0 due to numerics
           * dividng by values smaller delta with delta < 1 will increase this value, leading probably to detectable changes in the result
           * -> check if update is larger tol before dividng by delta*2
           */
          update = ((right-left)*(top-bot));
          if(update > tol_){
            update = update/(delta_*delta_);
          } else {
            update = 0.0;
          }

          /*
           * subtract area below diagonal (if any)
           * -> check if lower right corner lies below diagonal
           */
          if(bot-alpha < right-beta){
            /*
             * subtract triangle with side length (right - bot)/delta
             */
            Double triang = 0.5*((right-bot)*(right-bot));
            if(triang > tol_){
              triang = triang/(delta_*delta_);
            } else {
              triang = 0.0;
            }
            update -= triang;
          }

          /*
           * helper area is negative, so -
           */
          val -= update;
          /*
           * NOTE: always add up; negative sign comes via switching operators!
           */
          vec += entryVec*update;
        }

        /*
         * The cases below are called independently of the counter cnt.
         * During the first iteration (cnt==0), it will compute area 2.
         */
        if(isMinEntry){
          /*
           * Current entry is negative
           */
          right = beta+delta_;
          bool found = false;

          if(listIt != listEnd){
            listIt++;
            if(listIt->isMin() == isMinEntry){
              /*
               * element of same type detected -> ok for rotList, not ok for MinMaxList
               */
              if(isRot){
                /*
                 * set bottom to alpha
                 */
                bot = alpha;
                /*
                 * here we can directly set the next entry, too
                 * right is determined by the consecutive minimum, which just got found here
                 */
                found = true;
                right = listIt->getVal();
              } else {
                EXCEPTION("MinMaxList has to be alternating!")
              }
            } else {
              bot = listIt->getVal();
            }

            if(listIt != listEnd){
              listIt++;
              if(listIt->isMin() != isMinEntry){
                EXCEPTION("MinMaxList has to be alternating!")
                /*
                 * here we would expect the same entry type the the one of the current iteration
                 * e.g. a minimum after a maximum which followd a minimum
                 */
              } else {
                if(isRot&&found){
                  /*
                   * do nothing here -> we need the second min/max here but this value was already found in the previous entry
                   *
                   * in case of alternating rotList, we retrieved just the value of bot in the last case, so we need to find right
                   */
                } else {
                  right = listIt->getVal();
                }
              }
              listIt--;
            }
            listIt--;
          } else {
            bot = alpha;
          }

          top = alpha+delta_;
          left = entry;

          /*
           * restrict bound such, that only the bottom right corner lies below the diagonal
           */
          if(right-beta > top-alpha){
            /*
             * upper right corner below diagonal
             * -> set right bound to top
             */
            right = top;
          }
          if(bot-alpha < left-beta){
            /*
             * lower left corner below diagonal
             * -> set bottom bound to left
             */
            bot = left;
          }

          /*
           * calculate area
           * Note: even if right = left; or top = bot, the value update may be != 0 due to numerics
           * dividng by values smaller delta with delta < 1 will increase this value, leading probably to detectable changes in the result
           * -> check if update is larger tol before dividng by delta*2
           */
          update = ((right-left)*(top-bot));
          if(update > tol_){
            update = update/(delta_*delta_);
          } else {
            update = 0.0;
          }

          /*
           * subtract area below diagonal (if any)
           * -> check if lower right corner lies below diagonal
           */
          if(bot-alpha < right-beta){
            /*
             * subtract triangle with side length (right - bot)/delta
             */
            Double triang = 0.5*((right-bot)*(right-bot));
            if(triang > tol_){
              triang = triang/(delta_*delta_);
            } else {
              triang = 0.0;
            }
            update -= triang;
          }

          val -= update;
          /*
           * NOTE: always add up; negative sign comes via switching operators!
           */
          vec += entryVec*update;

        } else {
          /*
           * Current entry is positive
           */
          bot = alpha;
          bool found = false;

          if(listIt != listEnd){
            listIt++;
            if(listIt->isMin() == isMinEntry){
              /*
               * element of same type detected -> ok for rotList, not ok for MinMaxList
               */
              if(isRot){
                /*
                 * set right to beta+delta_
                 */
                right = beta+delta_;
                /*
                 * here we can directly set the next entry, too
                 * right is determined by the consecutive minimum, which just got found here
                 */
                found = true;
                bot = listIt->getVal();
              } else {
                EXCEPTION("MinMaxList has to be alternating!")
              }
            } else {
              right = listIt->getVal();
            }

            if(listIt != listEnd){
              listIt++;
              if(listIt->isMin() != isMinEntry){
                EXCEPTION("MinMaxList has to be alternating!")
                /*
                 * here we would expect the same entry type the the one of the current iteration
                 * e.g. a minimum after a maximum which followd a minimum
                 */
              } else {
                if(isRot && found){
                  /*
                   * do nothing here -> we need the second min/max here but this value was already found in the previous entry
                   *
                   * in case of alternating rotList, we retrieved just the value of right in the last case, so we need to find bot
                   */
                } else {
                  bot = listIt->getVal();
                }
              }
              listIt--;
            }
            listIt--;
          } else {
            right =  beta+delta_;
          }

          top = entry;
          left = beta;


          /*
           * restrict bound such, that only the bottom right corner lies below the diagonal
           */
          if(right-beta > top-alpha){
            /*
             * upper right corner below diagonal
             * -> set right bound to top
             */
            right = top;
          }
          if(bot-alpha < left-beta){
            /*
             * lower left corner below diagonal
             * -> set bottom bound to left
             */
            bot = left;
          }

          /*
           * calculate area
           * Note: even if right = left; or top = bot, the value update may be != 0 due to numerics
           * dividng by values smaller delta with delta < 1 will increase this value, leading probably to detectable changes in the result
           * -> check if update is larger tol before dividng by delta*2
           */
          update = ((right-left)*(top-bot));
          if(update > tol_){
            update = update/(delta_*delta_);
          } else {
            update = 0.0;
          }

          /*
           * subtract area below diagonal (if any)
           * -> check if lower right corner lies below diagonal
           */
          if(bot-alpha < right-beta){
            /*
             * subtract triangle with side length (right - bot)/delta
             */
            Double triang = 0.5*((right-bot)*(right-bot));
            if(triang > tol_){
              triang = triang/(delta_*delta_);
            } else {
              triang = 0.0;
            }
            update -= triang;
          }

          val += update;
          /*
           * NOTE: always add up; negative sign comes via switching operators!
           */
          vec += entryVec*update;
        }
        cnt++;
    }

    if(isRot){
      vec = vec*2.0;
    }

  } // end triangle

  retVal = val;
  retVec = vec;

  /*
   * After evaluating the element, we set preRotChanged_ and preMinMaxChanged_ to false.
   * Exception: if we use evalVersion_ == 4, this function is called to precompute the rotation state
   * for the evaluation of xPar. In a later step EvaluatePreisachElement_v4 is called check if the rotation
   * state changed. However, in that case the flag was already set to false here, so that we miss an evaluation.
   * Setting the flag to true after this function, does not work, as it would set this flag also for elements
   * were the rotation stayed the same.
   */
  if(isRot){
    if((evalVersion_!=4)&&(evalVersion_!=5)){
      preRotChanged_[idElem][idAlpha][idBeta] = false;
    }
  } else {
    preMinMaxChanged_[idElem][idAlpha][idBeta] = false;
  }

  return;
}


  void VectorPreisach::EvaluatePreisachElement_v4(std::list<ListEntryOld>& rotList, std::list<ListEntryOld>& mmList,
                                                  Double& retVal, Vector<Double>& retVec, Double alpha, Double beta,
                                                  UInt idElem, UInt idAlpha, UInt idBeta){
    /*
     * NEW:
     * computes the final vector value of the element alpha,beta by overlapping the already updated rotation state of the
     * element with the yet to evaluate switching state of the element
     *
     * Idea:
     * go overr all subareas which can be constructed from rotation list; for each subarea go over minmaxlist and get
     * the resulting subareas from this list; check if the minmax area overlaps with the rotation area;
     * if yes: add or subtract from final value
     * if no: try next area
     */

    Double val = 0;
    Vector<Double> vec = Vector<Double>(dim_);
    vec.Init();

    if((preMinMaxChanged_[idElem][idAlpha][idBeta] == false)&&(preRotChanged_[idElem][idAlpha][idBeta] == false)){
      /*
       * only if both, the rotation state as well as the minmax list did not change for this element, we return
       * the old value of the element
       */
      //retVal = switchingStates_[idElem][idAlpha][idBeta];-> comment this out; in case that element is not reevaluated, we get 0 back
      //-> thus we can see in the switching states, which element was reevaluated
      retVec = evaluatedState_[idElem][idAlpha][idBeta];
      return;
    }

    if(rotList.empty()||mmList.empty()){
      EXCEPTION("Lists should not be empty at this point!");
    }

    /*
     * flag for checking if we are on the diagonal alpha = beta
     * in that case we have to do some extra steps
     */
    bool triangle = false;
    if(idAlpha == idBeta){
      triangle = true;
    }
    Double update = 0.0;
    Double summedUpdate = 0.0;
    UInt rot_cnt = 0;
    UInt mm_cnt = 0;

    /*
     * we need both, the rectangular bounds of the rotation area as well as the one of the switching area
     * to calculate their overlap; for the overlap, we need bounds, too
     */
    Double rotTop = 0.0;
    Double rotBot = 0.0;
    Double rotLeft = 0.0;
    Double rotRight = 0.0;

    Double mmTop = 0.0;
    Double mmBot = 0.0;
    Double mmLeft = 0.0;
    Double mmRight = 0.0;

    Double ovTop = 0.0;
    Double ovBot = 0.0;
    Double ovLeft = 0.0;
    Double ovRight = 0.0;

    //bool isMinEntry_rot;
    bool isMinEntry_mm;
    /*
     * from minMaxList we get the scalar entry specifying the length
     * from rotList we get the rotation vector specifying the direction
     */
    //Double mmEntry;
    Vector<Double> rotEntryVec;

    std::list<ListEntryOld>::iterator rotListIt;
    std::list<ListEntryOld>::iterator mmListIt;
    /*
     * lists are not empty, so we can get iterator to last element
     * note that .end() is pointing to the element past the last element
     */
    std::list<ListEntryOld>::iterator rotListEnd = --(rotList.end());
    std::list<ListEntryOld>::iterator rotListStart = rotList.begin();

    std::list<ListEntryOld>::iterator mmListEnd = --(mmList.end());
    std::list<ListEntryOld>::iterator mmListStart = mmList.begin();

    /*
     * Outer loop -> go over rotation areas
     */
    for(rotListIt = rotList.begin(); rotListIt != rotList.end();){

      rotEntryVec = rotListIt->getVec();
      //isMinEntry_rot = rotListIt->isMin();

      getRectangleBounds(rotList,rotListIt,rotListEnd,alpha, beta, idElem, idAlpha, idBeta,
                         rot_cnt,triangle, true,rotTop,rotBot,rotLeft,rotRight);

      /*
       * Inner loop -> go over switching areas
       */
      summedUpdate = 0.0;
      update = 0.0;
      mm_cnt = 0;
      for(mmListIt = mmList.begin(); mmListIt != mmList.end(); ){

        //mmEntry = mmListIt->getVal();
        isMinEntry_mm = mmListIt->isMin();

        getRectangleBounds(mmList,mmListIt,mmListEnd,alpha, beta, idElem, idAlpha, idBeta,
                           mm_cnt,triangle, false,mmTop,mmBot,mmLeft,mmRight);

        /*
         * calculate overlapping of current switching area and current rotation area
         */
        bool success;
        success = clipRectangles(rotTop, rotBot, rotLeft, rotRight,
                       mmTop, mmBot, mmLeft, mmRight,
                       ovTop, ovBot, ovLeft, ovRight);

        if(success){
          /*
           * Areas overlap -> calculate size of area; clip if necessary; update output
           */
          update = ((ovRight-ovLeft)*(ovTop-ovBot));

          /*
           * Note: even if right = left; or top = bot, the value update may be != 0 due to numerics
           * dividng by values smaller delta with delta < 1 will increase this value, leading probably to detectable changes in the result
           * -> check if update is larger tol before dividng by delta*2
           */
  // Do this weightning to the final sum only

//          if(update > tol_){
//            update = update/(delta_*delta_);
//          } else {
//            update = 0.0;
//          }

          /*
           * subtract area below diagonal (if any)
           * -> check if lower right corner lies below diagonal
           */
          if(triangle == true){
            if(ovBot-alpha < ovRight-beta){
              /*
               * subtract triangle with side length (right - bot)/delta
               */
              Double triang = 0.5*((ovRight-ovBot)*(ovRight-ovBot));
              // Do this weightning to the final sum only

//              if(triang > tol_){
//                triang = triang/(delta_*delta_);
//              } else {
//                triang = 0.0;
//              }

              update -= triang;
            }
          }
          if(isMinEntry_mm && (mm_cnt == 0)){
            /*
             * helper area is positive, so +
             */
            summedUpdate += update;
          } else if(!isMinEntry_mm && (mm_cnt == 0)) {
            /*
             * helper area is positive, so +
             */
            summedUpdate -= update;
          } else if(isMinEntry_mm && (mm_cnt > 0)){
            /*
             * "regular" area is negative, so -
             */
            summedUpdate -= update;
          } else {
            /*
             * "regular" area is positive, so +
             */
            summedUpdate += update;
          }
        }

        if(mm_cnt > 0){
          /*
           * NOTE: area1 and area2 are both calculated using the first list entry, therefore we do not increase the iterator after
           * the first iteration
           */
           mmListIt++;
        }
        mm_cnt++;

      } // end inner loop

      /*
       * now we know the total weighting of this particular rotation area;
       * -> multiply and add up
       * NOTE: triangles do not have to be cut agailn; this is done already to the weights;
       * in the other forms (v1 - v3) the vector weights are cut, but later on multiplied with 2 again
       */
//      if((idAlpha == 0)&&(idBeta == 0)){
//        std::cout << "Element " << idAlpha << "," << idBeta << std::endl;
//        std::cout << "Outer Iteration " << rot_cnt << std::endl;
//        std::cout << "SummedUpdate: " << summedUpdate << std::endl;
//        std::cout << "rotEntryVec: " << std::endl;
//        std::cout << rotEntryVec.ToString() << std::endl;
//      }

      /*
       * the weighting with delta^2 can be problematic, if the summed update only results from
       * numerical pollution; in that case the weighting would increase its value;
       * this effect is actually visible in the final results
       */
//      if(summedUpdate > tol_){ //it does not work if we compare to tol here; summedUpdate seems to be very small ...; results fit however
//        summedUpdate = summedUpdate/(delta_*delta_);
//      } else {
//        summedUpdate = 0.0;
//      }
      summedUpdate = summedUpdate/(delta_*delta_);

      val += summedUpdate; // if done correct, val should contain the average switching state of the element
      // as it would return from version 1-3
      vec += rotEntryVec*summedUpdate;

      if(rot_cnt > 0){
        /*
         * NOTE: area1 and area2 are both calculated using the first list entry, therefore we do not increase the iterator after
         * the first iteration
         */
        rotListIt++;
      }
      rot_cnt++;

    } // end outer loop

    retVal = val;
    retVec = vec;

    /*
     * After evaluating the element, we set both preRotChanged_ and preMinMaxChanged_ to false.
     */
    preRotChanged_[idElem][idAlpha][idBeta] = false;
    preMinMaxChanged_[idElem][idAlpha][idBeta] = false;

    return;
  }

void VectorPreisach::EvaluatePreisachElement_v3(std::list<ListEntryOld>& list, Double& retVal, Vector<Double>& retVec, Double alpha, Double beta, UInt idElem, UInt idAlpha, UInt idBeta, bool isRot){

  Double val = 0;
  Vector<Double> vec = Vector<Double>(dim_);
  vec.Init();

  if((isRot)&&(preRotChanged_[idElem][idAlpha][idBeta] == false)){
    /*
     * no change this last time -> return previous value
     */
    retVal = val;
    retVec = rotationStates_[idElem][idAlpha][idBeta];
    return;
  } else if((!isRot)&&(preMinMaxChanged_[idElem][idAlpha][idBeta] == false)){
    /*
     * no change this last time -> return previous value
     */
    retVal = switchingStates_[idElem][idAlpha][idBeta];
    retVec = vec;
    return;
  }

  if(list.empty()){
    EXCEPTION("List should not be empty at this point!");
  }

  /*
   * flag for checking if we are on the diagonal alpha = beta
   * in that case we have to do some extra steps
   */
  bool triangle = false;
  if(idAlpha == idBeta){
    triangle = true;
  }
  Double update = 0.0;
  UInt cnt = 0;
  Double top = 0.0;
  Double bot = 0.0;
  Double left = 0.0;
  Double right = 0.0;
  Double factor = 1.0;

  bool isMinEntry;
  //Double entry;
  Vector<Double> entryVec;

  std::list<ListEntryOld>::iterator listIt;
  /*
   * list is not empty, so we can get iterator to last element
   * note that .end() is pointing to the element past the last element
   */
  std::list<ListEntryOld>::iterator listEnd = --(list.end());
  std::list<ListEntryOld>::iterator listStart = list.begin();

//  if(((idAlpha == 7)&&(idBeta == 7))||((idAlpha == 8)&&(idBeta == 7))||((idAlpha == 8)&&(idBeta == 8))){
//    std::cout << "idAlpha,idBeta: " << idAlpha << "," << idBeta << std::endl;
//    std::cout << "alpha: " << alpha << std::endl;
//    std::cout << "beta: " << beta << std::endl;
//    std::cout << "delta: " << delta_ << std::endl;
//    std::cout << "starting val: " << val << std::endl;
//  }

  for(listIt = list.begin(); listIt != list.end();){
    //entry = listIt->getVal();

    entryVec = listIt->getVec();
    isMinEntry = listIt->isMin();

    factor = getRectangleBounds(list,listIt,listEnd,alpha, beta, idElem, idAlpha, idBeta, cnt,triangle, isRot,top,bot,left,right);

    /*
     * calculate area
     * Note: even if right = left; or top = bot, the value update may be != 0 due to numerics
     * dividng by values smaller delta with delta < 1 will increase this value, leading probably to detectable changes in the result
     * -> check if update is larger tol before dividng by delta*2
     */
    update = ((right-left)*(top-bot));
    if(update > tol_){
      update = update/(delta_*delta_);
    } else {
      update = 0.0;
    }

//    if(((idAlpha == 7)&&(idBeta == 7))||((idAlpha == 8)&&(idBeta == 7))||((idAlpha == 8)&&(idBeta == 8))){
//      std::cout << "List entry: " << listIt->getVal() << std::endl;
//      std::cout << "left,right: " << left << "," << right << std::endl;
//      std::cout << "top,bot: " << top << "," << bot << std::endl;
//      std::cout <<  std::setprecision(12) << "update: " << update << std::endl;
//    }


    /*
     * subtract area below diagonal (if any)
     * -> check if lower right corner lies below diagonal
     */
    if(triangle == true){
      if(bot-alpha < right-beta){
        /*
         * subtract triangle with side length (right - bot)/delta
         */
        Double triang = 0.5*((right-bot)*(right-bot));
        if(triang > tol_){
          triang = triang/(delta_*delta_);
        } else {
          triang = 0.0;
        }
//        if(((idAlpha == 7)&&(idBeta == 7))||((idAlpha == 8)&&(idBeta == 7))||((idAlpha == 8)&&(idBeta == 8))){
//          std::cout <<  std::setprecision(12) << "triang: " << triang << std::endl;
//        }
        update -= triang;
      }
    }
    // factor will be positive in all cases, except if area1 would overlap area3
    // in that case, the values of top,bot,left and right are adapted by the function getRectangularBounds
    // the resulting area needs to be added to the opposite area, so factor = -1
    // do this multiplication only after the triangular area was subtracted
    update = update*factor;

//    if(((idAlpha == 7)&&(idBeta == 7))||((idAlpha == 8)&&(idBeta == 7))||((idAlpha == 8)&&(idBeta == 8))){
//      std::cout <<  std::setprecision(12) << "total update: " << update << std::endl;
//    }

    if(isMinEntry && (cnt == 0)){
      /*
       * helper area is positive, so +
       */
      val += update;
    } else if(!isMinEntry && (cnt == 0)) {
      /*
       * helper area is positive, so +
       */
      val -= update;
    } else if(isMinEntry && (cnt > 0)){
      /*
       * "regular" area is negative, so -
       */
      val -= update;
    } else {
      /*
       * "regular" area is positive, so +
       */
      val += update;
    }

    /*
     * resulting vector for the element will be the weighted sum of the stored, averaged rotation vectors
     * the weights will be the areas corresponding to those rotation vectors
     * -> only relevant for rotList
     *
     * NOTE: always add up; negative sign comes via switching operators!
     */
    vec += entryVec*update;

    if(cnt > 0){
      /*
       * NOTE: area1 and area2 are both calculated using the first list entry, therefore we do not increase the iterator after
       * the first iteration
       */
      listIt++;
    }
    cnt++;
  }

  if(isRot && (triangle == true)){
    vec = vec*2.0;
  }

//  if(((idAlpha == 7)&&(idBeta == 7))||((idAlpha == 8)&&(idBeta == 7))||((idAlpha == 8)&&(idBeta == 8))){
//    std::cout << "final value: "<<  std::setprecision(12) << val << std::endl;
//  }

  retVal = val;
  retVec = vec;

  /*
   * After evaluating the element, we set preRotChanged_ and preMinMaxChanged_ to false.
   * Exception: if we use evalVersion_ == 4, this function is called to precompute the rotation state
   * for the evaluation of xPar. In a later step EvaluatePreisachElement_v4 is called check if the rotation
   * state changed. However, in that case the flag was already set to false here, so that we miss an evaluation.
   * Setting the flag to true after this function, does not work, as it would set this flag also for elements
   * were the rotation stayed the same.
   */
  if(isRot){
    if((evalVersion_!=4)&&(evalVersion_!=5)){
      preRotChanged_[idElem][idAlpha][idBeta] = false;
    }
  } else {
    preMinMaxChanged_[idElem][idAlpha][idBeta] = false;
  }

  return;
}

  void VectorPreisach::EvaluatePreisachElement_v2(std::list<ListEntryOld>& list, Double& retVal, Vector<Double>& retVec, Double alpha, Double beta, UInt idElem, UInt idAlpha, UInt idBeta, bool isRot){

    Double val = 0;
    Vector<Double> vec = Vector<Double>(dim_);
    vec.Init();

    if((isRot)&&(preRotChanged_[idElem][idAlpha][idBeta] == false)){
      /*
       * no change this last time -> return previous value
       */
      retVal = val;
      retVec = rotationStates_[idElem][idAlpha][idBeta];
      return;
    } else if((!isRot)&&(preMinMaxChanged_[idElem][idAlpha][idBeta] == false)){
      /*
       * no change this last time -> return previous value
       */
      retVal = switchingStates_[idElem][idAlpha][idBeta];
      retVec = vec;
      return;
    }

    if(list.empty()){
      EXCEPTION("List should not be empty at this point!");
    }

    /*
     * flag for checking if we are on the diagonal alpha = beta
     * in that case we have to do some extra steps
     */
    bool triangle = false;
    if(idAlpha == idBeta){
      triangle = true;
    }
    Double update = 0.0;
    UInt cnt = 0;
    Double top = 0.0;
    Double bot = 0.0;
    Double left = 0.0;
    Double right = 0.0;

    bool isMinEntry;
    Double entry;
    Vector<Double> entryVec;

    std::list<ListEntryOld>::iterator listIt;
    /*
     * list is not empty, so we can get iterator to last element
     * note that .end() is pointing to the element past the last element
     */
    std::list<ListEntryOld>::iterator listEnd = --(list.end());

    for(listIt = list.begin(); listIt != list.end(); listIt++){

      entry = listIt->getVal();
      entryVec = listIt->getVec();
      isMinEntry = listIt->isMin();

      if(isMinEntry && (cnt == 0)){
        /*
         * First entry in list is a minimum; two possibilities:
         * 1. element is not wiped out yet -> minimum lies between beta and beta+delta
         * 2. element is wiped out by minimum or got initialized or overwritten with minimum -> minimum has value beta
         *
         * for case 2, area 1 will be 0 as left = right
         *
         * Note that the area 2 is already a std area, so instead of writing it once again into this case, we just
         * execute it additionally to this one here.
         */

        /*
         * area 1
         */
        left = beta;
        right = entry;

        if(triangle == true){
          if( (idAlpha+idBeta+1 == numRows_)&&(wipedOut_[idElem][idBeta] == false) ){
            /*
             * In contrast to the square element case, we have to distinguish between alpha = -beta elements
             * which are not wiped out and all other elements as area0 will not automatically get 0
             */
            top = alpha+delta_ - entry + beta;
          } else {
            top = alpha+delta_;
          }
        } else {
          top = alpha+delta_ - entry + beta;
        }

        if(listIt != listEnd){
          listIt++;
          if(listIt->isMin() == isMinEntry){
            /*
             * element of same type detected -> ok for rotList, not ok for MinMaxList
             */
            if(isRot){
              /*
               * set bot to alpha
               */
              bot = alpha;
            } else {
              EXCEPTION("MinMaxList has to be alternating!")
            }
          } else {
            bot = listIt->getVal();
          }
          listIt--;
        } else {
          bot = alpha;
        }

        if(triangle == true){
          /*
           * restrict bound such, that only the bottom right corner lies below the diagonal
           */
          if(right-beta > top-alpha){
            /*
             * upper right corner below diagonal
             * -> set right bound to top
             */
            right = top;
          }
          if(bot-alpha < left-beta){
            /*
             * lower left corner below diagonal
             * -> set bottom bound to left
             */
            bot = left;
          }
        }
        /*
         * calculate area
         */
        /*
         * calculate area
         * Note: even if right = left; or top = bot, the value update may be != 0 due to numerics
         * dividng by values smaller delta with delta < 1 will increase this value, leading probably to detectable changes in the result
         * -> check if update is larger tol before dividng by delta*2
         */
        update = ((right-left)*(top-bot));
        if(update > tol_){
          update = update/(delta_*delta_);
        } else {
          update = 0.0;
        }

        /*
         * subtract area below diagonal (if any)
         * -> check if lower right corner lies below diagonal
         */
        if(triangle == true){
          if(bot-alpha < right-beta){
            /*
             * subtract triangle with side length (right - bot)/delta
             */
            Double triang = 0.5*((right-bot)*(right-bot));
            if(triang > tol_){
              triang = triang/(delta_*delta_);
            } else {
              triang = 0.0;
            }
            update -= triang;
          }
        }
        /*
         * helper area is positive, so +
         */
        val += update;

        /*
         * resulting vector for the element will be the weighted sum of the stored, averaged rotation vectors
         * the weights will be the areas corresponding to those rotation vectors
         * -> only relevant for rotList
         *
         * NOTE: always add up; negative sign comes via switching operators!
         */
        vec += entryVec*update;

      } else if(!isMinEntry && (cnt == 0)) {
        /*
         * First entry in list is a maximum; two possibilities:
         * 1. element is not wiped out yet -> maximum lies between alpha and alpha+delta
         * 2. element is wiped out by maximum or got initialized or overwritten with maximum -> maximum has value alpha+delta
         *
         * for case 2, area 1 will be 0 as top = bot
         *
         * Note that the area 2 is already a std area, so instead of writing it once again into this case, we just
         * execute it additionally to this one here.
         */
        /*
         * area 1
         */
        if(listIt != listEnd){
          listIt++;
          if(listIt->isMin() == isMinEntry){
            /*
             * element of same type detected -> ok for rotList, not ok for MinMaxList
             */
            if(isRot){
              /*
               * set bottom to beta+delta_
               */
              right = beta+delta_;
            } else {
              EXCEPTION("MinMaxList has to be alternating!")
            }
          } else {
            right = listIt->getVal();
          }
          listIt--;
        } else {
          right = beta+delta_;
        }

        if(triangle == true){
          if( (idAlpha+idBeta+1 == numRows_)&&(wipedOut_[idElem][idBeta] == false) ){
            /*
             * In contrast to the square element case, we have to distinguish between alpha = -beta elements
             * which are not wiped out and all other elements as area0 will not automatically get 0
             */
            left = alpha+delta_ - entry + beta;
          } else {
            left = beta;
          }
        } else {
          left = alpha+delta_ - entry + beta;
        }

        top = alpha+delta_;
        bot = entry;

        if(triangle == true){
          /*
           * restrict bound such, that only the bottom right corner lies below the diagonal
           */
          if(right-beta > top-alpha){
            /*
             * upper right corner below diagonal
             * -> set right bound to top
             */
            right = top;
          }
          if(bot-alpha < left-beta){
            /*
             * lower left corner below diagonal
             * -> set bottom bound to left
             */
            bot = left;
          }
        }

        /*
         * calculate area
         * Note: even if right = left; or top = bot, the value update may be != 0 due to numerics
         * dividng by values smaller delta with delta < 1 will increase this value, leading probably to detectable changes in the result
         * -> check if update is larger tol before dividng by delta*2
         */
        update = ((right-left)*(top-bot));
        if(update > tol_){
          update = update/(delta_*delta_);
        } else {
          update = 0.0;
        }

        /*
         * subtract area below diagonal (if any)
         * -> check if lower right corner lies below diagonal
         */
        if(triangle == true){
          if(bot-alpha < right-beta){
            /*
             * subtract triangle with side length (right - bot)/delta
             */
            Double triang = 0.5*((right-bot)*(right-bot));
            if(triang > tol_){
              triang = triang/(delta_*delta_);
            } else {
              triang = 0.0;
            }
            update -= triang;
          }
        }

        /*
         * helper area is negative, so -
         */
        val -= update;
        /*
         * NOTE: always add up; negative sign comes via switching operators!
         */
        vec += entryVec*update;
      }

      /*
       * The cases below are called independently of the counter cnt.
       * During the first iteration (cnt==0), it will compute area 2.
       */
      if(isMinEntry){
        /*
         * Current entry is negative
         */
        right = beta+delta_;
        bool found = false;

        if(listIt != listEnd){
          listIt++;
          if(listIt->isMin() == isMinEntry){
            /*
             * element of same type detected -> ok for rotList, not ok for MinMaxList
             */
            if(isRot){
              /*
               * set bottom to alpha
               */
              bot = alpha;
              /*
               * here we can directly set the next entry, too
               * right is determined by the consecutive minimum, which just got found here
               */
              found = true;
              right = listIt->getVal();
            } else {
              EXCEPTION("MinMaxList has to be alternating!")
            }
          } else {
            bot = listIt->getVal();
          }

          if(listIt != listEnd){
            listIt++;
            if(listIt->isMin() != isMinEntry){
              EXCEPTION("MinMaxList has to be alternating!")
              /*
               * here we would expect the same entry type the the one of the current iteration
               * e.g. a minimum after a maximum which followd a minimum
               */
            } else {
              if(isRot&&found){
                /*
                 * do nothing here -> we need the second min/max here but this value was already found in the previous entry
                 *
                 * in case of alternating rotList, we retrieved just the value of bot in the last case, so we need to find right
                 */
              } else {
                right = listIt->getVal();
              }
            }
            listIt--;
          }
          listIt--;
        } else {
          bot = alpha;
        }

        top = alpha+delta_;
        left = entry;

        if(triangle == true){
          /*
           * restrict bound such, that only the bottom right corner lies below the diagonal
           */
          if(right-beta > top-alpha){
            /*
             * upper right corner below diagonal
             * -> set right bound to top
             */
            right = top;
          }
          if(bot-alpha < left-beta){
            /*
             * lower left corner below diagonal
             * -> set bottom bound to left
             */
            bot = left;
          }
        }

        /*
         * calculate area
         * Note: even if right = left; or top = bot, the value update may be != 0 due to numerics
         * dividng by values smaller delta with delta < 1 will increase this value, leading probably to detectable changes in the result
         * -> check if update is larger tol before dividng by delta*2
         */
        update = ((right-left)*(top-bot));
        if(update > tol_){
          update = update/(delta_*delta_);
        } else {
          update = 0.0;
        }

        /*
         * subtract area below diagonal (if any)
         * -> check if lower right corner lies below diagonal
         */
        if(triangle == true){
          if(bot-alpha < right-beta){
            /*
             * subtract triangle with side length (right - bot)/delta
             */
            Double triang = 0.5*((right-bot)*(right-bot));
            if(triang > tol_){
              triang = triang/(delta_*delta_);
            } else {
              triang = 0.0;
            }
            update -= triang;
          }
        }

        val -= update;
        /*
         * NOTE: always add up; negative sign comes via switching operators!
         */
        vec += entryVec*update;

      } else {
        /*
         * Current entry is positive
         */
        bot = alpha;
        bool found = false;

        if(listIt != listEnd){
          listIt++;
          if(listIt->isMin() == isMinEntry){
            /*
             * element of same type detected -> ok for rotList, not ok for MinMaxList
             */
            if(isRot){
              /*
               * set right to beta+delta_
               */
              right = beta+delta_;
              /*
               * here we can directly set the next entry, too
               * right is determined by the consecutive minimum, which just got found here
               */
              found = true;
              bot = listIt->getVal();
            } else {
              EXCEPTION("MinMaxList has to be alternating!")
            }
          } else {
            right = listIt->getVal();
          }

          if(listIt != listEnd){
            listIt++;
            if(listIt->isMin() != isMinEntry){
              EXCEPTION("MinMaxList has to be alternating!")
              /*
               * here we would expect the same entry type the the one of the current iteration
               * e.g. a minimum after a maximum which followd a minimum
               */
            } else {
              if(isRot && found){
                /*
                 * do nothing here -> we need the second min/max here but this value was already found in the previous entry
                 *
                 * in case of alternating rotList, we retrieved just the value of right in the last case, so we need to find bot
                 */
              } else {
                bot = listIt->getVal();
              }
            }
            listIt--;
          }
          listIt--;
        } else {
          right =  beta+delta_;
        }

        top = entry;
        left = beta;

        if(triangle == true){
          /*
           * restrict bound such, that only the bottom right corner lies below the diagonal
           */
          if(right-beta > top-alpha){
            /*
             * upper right corner below diagonal
             * -> set right bound to top
             */
            right = top;
          }
          if(bot-alpha < left-beta){
            /*
             * lower left corner below diagonal
             * -> set bottom bound to left
             */
            bot = left;
          }
        }

        /*
         * calculate area
         * Note: even if right = left; or top = bot, the value update may be != 0 due to numerics
         * dividng by values smaller delta with delta < 1 will increase this value, leading probably to detectable changes in the result
         * -> check if update is larger tol before dividng by delta*2
         */
        update = ((right-left)*(top-bot));
        if(update > tol_){
          update = update/(delta_*delta_);
        } else {
          update = 0.0;
        }
        /*
         * subtract area below diagonal (if any)
         * -> check if lower right corner lies below diagonal
         */
        if(triangle == true){
          if(bot-alpha < right-beta){
            /*
             * subtract triangle with side length (right - bot)/delta
             */
            Double triang = 0.5*((right-bot)*(right-bot));
            if(triang > tol_){
              triang = triang/(delta_*delta_);
            } else {
              triang = 0.0;
            }
            update -= triang;
          }
        }

        val += update;
        /*
         * NOTE: always add up; negative sign comes via switching operators!
         */
        vec += entryVec*update;
      }
      cnt++;
  }

  if(isRot && (triangle == true)){
    vec = vec*2.0;
  }

  retVal = val;
  retVec = vec;

  /*
   * After evaluating the element, we set preRotChanged_ and preMinMaxChanged_ to false.
   * Exception: if we use evalVersion_ == 4, this function is called to precompute the rotation state
   * for the evaluation of xPar. In a later step EvaluatePreisachElement_v4 is called check if the rotation
   * state changed. However, in that case the flag was already set to false here, so that we miss an evaluation.
   * Setting the flag to true after this function, does not work, as it would set this flag also for elements
   * were the rotation stayed the same.
   */
  if(isRot){
    if((evalVersion_!=4)&&(evalVersion_!=5)){
      preRotChanged_[idElem][idAlpha][idBeta] = false;
    }
  } else {
    preMinMaxChanged_[idElem][idAlpha][idBeta] = false;
  }

  return;
}

  Vector<Double> VectorPreisach::computeValue_vec(Vector<Double>& xVal, Integer idElem){

    /*
     * Determine the current rotational threshold
     */
    Double X_thres = std::pow((xVal.NormL2()/Xsaturated_),rotationalResistance_);

    /*
     * Get current direction
     */
    Vector<Double> e_u = Vector<Double>(xVal.GetSize());

    if(xVal.NormL2() > tol_){
      e_u = xVal/xVal.NormL2();
    } else {
      e_u.Init(0.0);
    }

    /*
     * Storage for return values
     */
    Double retVal;
    Vector<Double> retVec = Vector<Double>(dim_);
    retVec.Init();

    /*
     * Area of one element
     * -> Note that triangular elements already are scaled by 0.5 via the switchingState
     */
    Double area = delta_*delta_;

    /*
     * Storage for element value
     */
    Vector<Double> Yout = Vector<Double>(dim_);
    Yout.Init(0.0);

    /*
     * Optimized version of the evaluation (model-based optimization, not optimizations due to more efficient algorithms!)
     *
     * Model based facts:
     *  1. the lower triangle (alpha <= 0) and the upper triangle (beta >= 0) ALWAYS have e_u as rotation state
     *     as the corresponding switch xThresh >= 0 for all real input values and all values of rotationalResistance_
     *  2. the switching state of each element is determined only by the parallel component xPar = inner(u_in,rot_state),
     *     where rot_state is the rotational state of the element, resulting from the current and past input
     *  3. the switching state is +1 if xPar >= alpha and -1 if xPar <= beta
     *
     * Consequences:
     *  Combining 1. and 2. we can conclude, that the value xPar for all elements in the lower and upper triangle
     *  always will be equal to the absolute value of u_in. Therewith it will always be >= 0.
     *  In combination with 3., the switching state in the lower triangle will always have be +1 in each element.
     *
     * Optimizations:
     *  1. As the switching states in the lower triangle are +1 for each element and for each point in time and all
     *  rotation states are the same inside this triangle (see 1.) we can simplify the evaluation of the total
     *  rotation state of the triangle, by just multiplying the summed up weights (all belong to a switching state of
     *  +1) with the current rotation input e_u.
     *
     *  2. The rotation vector in the upper triangle is the same for each element. Therewith, we do not have to evaluated
     *  the weighted sum of the rotation states (i.e. go over each element and multiply the switching state with the
     *  rotation state; in version 4 even with clipping of the overlapping rectangles etc) but can directly evaluate
     *  the overall switching state (switches * weights) and multiply the result with e_u.
     *  This is similar to the first optimization, except that we have to reevaluate the sum during each step as the
     *  switching can change inside this upper triangle.
     *
     *  3. The remaining square part of the Preisach plane is evaluated with the evaluation function v4 (or maybe v3, if
     *  the little improvement in the resulting states is not needed in the end).l
     *
     */
    if((evalVersion_ == 5)||(evalVersion_ == 6)){

      /*
       * even number of rows
       *  _ _ _ _
       * |   |  /      lower triangle -> for i = 0; i < numRows/2; for j = 0; j <= i
       * |_ _|/        square         -> for i = numRows/2; i < numRows; for j = 0; j < numRows/2
       * |  /          upper triangle -> for i = numRows/2; i < numRows; for j = numRows/2; j <= i
       * |/
       *
       *
       * odd number of rows
       *  _ _ _ _ _
       * |     |  /                                         lower triangle -> for i = 0; i < roundDown(numRows/2)
       * |     |/                                                             for j = 0; j <= i
       * |_ _ /|< alpha = 0 in middle of this elements      square -> for i = roundDown(numRows/2); i < numRows
       * |  / ^                                                       for j = 0; j < roundUp(numRows/2)
       * |/   beta = 0 in middle of this elements           upper triangle -> for i = roundDown(numRows/2); i < numRows
       *                                                                      for j = roundUp(numRows/2); j <= i
       *
       */

      Vector<Double> rotState = Vector<Double>(dim_);
      if(isTesting_ == false){
        rotState = e_u;
      } else {
        // in the testing case, the rotation is fixed everywhere
        // instead of taking the new rotation state e_u we load
        // this fixed state
        rotState = rotStateTesting_;
      }

      /*
       * Calculate value of lower triangle
       * -> lowerTriangleValue_ was already evaluated in constructor
       */
      Yout = rotState*lowerTriangleValue_;

      /*
       * Calculate value of upper triangle
       * -> rotation state is always e_u
       */
      Double xPar = xVal.Inner(rotState);
      xPar /= Xsaturated_;

      Double upperTriangleValue = 0;
      Double alpha;
      Double beta;

      for (UInt i = (UInt) (numRows_/2.0); i < numRows_; i++){
        alpha = -1 + i*delta_;

        for (UInt j = (UInt)(numRows_/2.0) + numRows_%2; j <= i; j++){
          beta = -1 + j*delta_;
          /*
           * Here we calculate the switching state of the element as usual in two steps
           * 1. update list with current xPar
           * 2. evaluate element via EvaluatePreisachElement_v<4
           */
          std::list<ListEntryOld>& switchList = preMinMaxValues_[idElem][i][j];

          // note that the third parameter (here e_u) is not used inside UpdateList, if the last flag is false
          UpdateList(switchList, xPar, e_u, alpha, beta, idElem, i, j, false);

          EvaluatePreisachElement_v3(switchList, retVal, retVec, alpha, beta, idElem, i, j, false);

          switchingStates_[idElem][i][j] = retVal; // needed here as the stored switching state is returned by
          // EvaluaePreisachElement_v<4 if xPar did not change

          upperTriangleValue += retVal*preisachWeights_[i][j];

          rotationStates_[idElem][i][j] = rotState;
          evaluatedState_[idElem][i][j] = rotationStates_[idElem][i][j]*retVal;
        }
      }
      //std::cout << "upperTriangleValue: " << upperTriangleValue*(delta_*delta_) << std::endl;

      Yout += rotState*upperTriangleValue;

      /*
       * Calculate value of square
       */
      if(evalVersion_ == 5){
        /*
         * calculate rotation state with eval version 3 first, then
         * use EvaluatePreisachElement_v4 to overlap rotation state and switching state
         */
        for (UInt i = (UInt) (numRows_/2.0); i < numRows_; i++){
          alpha = -1 + i*delta_;

          for (UInt j = 0; j < (UInt)(numRows_/2.0) + numRows_%2; j++){
            beta = -1 + j*delta_;

           /*
            * -> std procedure like in version 4:
            * 1. if testing = false -> update rotation list
            * 2. evaluate rotation list with evaluatePreisachElement_v<4
            * 3. calculate xPar
            * 4. update switching list
            * 5. evaluate switching list and rotation list together with evaluatePreisachElement_v4
            * 6. sum up
            */

            std::list<ListEntryOld>& rotList = preRotValues_[idElem][i][j];

            if(isTesting_ == false){
              UpdateList(rotList, X_thres, e_u, alpha, beta, idElem, i, j, true);
            }

            EvaluatePreisachElement_v3(rotList, retVal, retVec, alpha, beta, idElem, i, j, true);

            rotationStates_[idElem][i][j] = retVec; // needed as EvaluatePreisachElement_v<4 with list = rotList returns this value
            // if e_u did not change since last time

            // new: use normalized vector -> explanation below; version 3/4
            Vector<Double> rotNormalized = Vector<Double>(dim_);

            if(rotationStates_[idElem][i][j].NormL2() > tol_){
              rotNormalized = rotationStates_[idElem][i][j]/rotationStates_[idElem][i][j].NormL2();
            } else {
              rotNormalized.Init(0.0);
            }


            //xPar = xVal.Inner(rotationStates_[idElem][i][j]);
            xPar = xVal.Inner(rotNormalized);
            xPar /= Xsaturated_;

            std::list<ListEntryOld>& switchList = preMinMaxValues_[idElem][i][j];

            UpdateList(switchList, xPar, e_u, alpha, beta, idElem, i, j, false);

            EvaluatePreisachElement_v4(rotList,switchList, retVal, retVec, alpha, beta, idElem, i, j);

            evaluatedState_[idElem][i][j] = retVec; // needed as EvaluatePreisachElement_v4 returns
            // this value, if neither e_u nor xPar changed during last iteration

            Yout += evaluatedState_[idElem][i][j]*preisachWeights_[i][j];
          }
        }

      } else {
        /*
         * caluclate rotation state with eval version 3 first, then
         * evaluate switching state and multiply both values
         */
        for (UInt i = (UInt) (numRows_/2.0); i < numRows_; i++){
          alpha = -1 + i*delta_;

          for (UInt j = 0; j < (UInt)(numRows_/2.0) + numRows_%2; j++){
            beta = -1 + j*delta_;

           /*
            * -> std procedure like in version < 4:
            * 1. if testing = false -> update rotation list
            * 2. evaluate rotation list with evaluatePreisachElement_v<4
            * 3. calculate xPar
            * 4. update switching list
            * 5. evaluate switching list with evaluatePreisachElement_v<4
            * 6. sum up
            */

            std::list<ListEntryOld>& rotList = preRotValues_[idElem][i][j];

            if(isTesting_ == false){
              UpdateList(rotList, X_thres, e_u, alpha, beta, idElem, i, j, true);
            }

            EvaluatePreisachElement_v3(rotList, retVal, retVec, alpha, beta, idElem, i, j, true);

            rotationStates_[idElem][i][j] = retVec; // needed as EvaluatePreisachElement_v<4 with list = rotList returns this value
            // if e_u did not change since last time

            // new: use normalized vector -> explanation below; version 3/4
            Vector<Double> rotNormalized = Vector<Double>(dim_);

            if(rotationStates_[idElem][i][j].NormL2() > tol_){
              rotNormalized = rotationStates_[idElem][i][j]/rotationStates_[idElem][i][j].NormL2();
            } else {
              rotNormalized.Init(0.0);
            }

            //xPar = xVal.Inner(rotationStates_[idElem][i][j]);
            xPar = xVal.Inner(rotNormalized);
            xPar /= Xsaturated_;

            std::list<ListEntryOld>& switchList = preMinMaxValues_[idElem][i][j];

            UpdateList(switchList, xPar, e_u, alpha, beta, idElem, i, j, false);

            EvaluatePreisachElement_v3(switchList, retVal, retVec, alpha, beta, idElem, i, j, false);

            switchingStates_[idElem][i][j] = retVal; // important as EvaluatePreisachElement_v<4 with list = switchList returns this
            // value if xPar did not change since last time

            evaluatedState_[idElem][i][j] = rotationStates_[idElem][i][j]*retVal; // needed as EvaluatePreisachElement_v4 returns
            // this value, if neither e_u nor xPar changed during last iteration

            Yout += evaluatedState_[idElem][i][j]*preisachWeights_[i][j];
          }
        }
      }
    } // end optimized version
    else {

      /*
       * Iterate over Preisach plane
       */
      Double alpha;
      Double beta;
      for (UInt i = 0; i < numRows_; i++){
        alpha = -1 + i*delta_;

        for (UInt j = 0; j <= i; j++){
          beta = -1 + j*delta_;

  //        std::cout << "Start" << std::endl;
  //        std::cout << "idAlpha, idBeta: " <<i << "," << j <<std::endl;
  //        std::cout << "preMinMaxChanged_[idElem][idAlpha][idBeta]: " <<preMinMaxChanged_[idElem][i][j] <<std::endl;
  //        std::cout << "preRotChanged_[idElem][idAlpha][idBeta]: " <<preRotChanged_[idElem][i][j] <<std::endl;

          /*
           * ROTATION PART
           */
          /*
           * Get current list storing the history of the elements rotation
           */
          std::list<ListEntryOld>& rotList = preRotValues_[idElem][i][j];

          /*
           * Overview of eval versions:
           *
           * 0. Classical approach:
           *  Data structures:
           *  - one matrix storing the switching states (std elements only +1,-1; +0.5, -0.5 on alpha = beta)
           *  - one matrix storing the rotation vector
           *  Evaluation:
           *  - determine if rotation state changes due to xThres
           *  - calculate xPar based on rotation state of the element
           *  - determine if element switched due to value of xPar
           *  -> element value = switching state * rotation vector * preisach weight
           *
           * 1.-3. List approaches:
           *  Data structures:
           *  - one matrix storing the switching states; each element has linked list which stores
           *     local minima and maxima of xPar; element value can be any double between -1 and +1
           *  - one matrix storing the rotation vector; each element has linked list which stores
           *     local minima and maxima of xThres and corresponding vector e_u;
           *     rotation state of element is reconstructed from list
           *  Evaluation:
           *  - update list of rotation states with input pair xThres, e_u
           *  - calculate rotation state and store it in matrix
           *  - calculate xPar based on the AVERAGE rotation state of the element
           *  - update list of switching states with input xPar
           *  - calculate switching state of element
           *  -> element value = average switching state * average rotation vector * preisach weight
           *
           *  Differences between version 1,2 and 3 -> evaluation function gets more compact
           *
           * 4. Stacked list approach:
           *  Data structures:
           *  - one matrix storing the switching states; each element has linked list which stores
           *    local minima and maxima of xPar; element value can be any double between -1 and +1
           *  - one matrix storing the rotation vector; each element has linked list which stores
           *    local minima and maxima of xThres and corresponding vector e_u;
           *    rotation state of element is reconstructed from list
           *  Evaluation:
           *  - update list of rotation states with input pair xThres, e_u
           *  - calculate rotation state and store it in matrix
           *  - calculate xPar based on the AVERAGE rotation state of the element
           *  - update list of switching states with input xPar
           *  - iterate over list or rotation states; for each subelement defined by a pair xThres, e_u
           *    calculate overlapping area between this subelement and the subelement resulting from a similar
           *    decomposition based on value xPar stored in switching list; return weighted average of the
           *    rotation state, where the weights are the switching states
           *  -> element value = weighted averaged rotation state * preisach weight
           *
           */
          if(evalVersion_ == 0){
            /*
             * Classical approach -> check if threshold is larger than alpha or smaller than -beta
             * -> if true, set RotationState to e_u
             */
            if((X_thres > alpha)||(-X_thres < beta)){
              if(isTesting_ == false){
                rotationStates_[idElem][i][j] = e_u;
              }
            }
          } else {
            /*
             * Weighted approach -> get rotation value by a weighted sum of old and current e_u
             */
            /*
             * In the testing case, we initialize the rotation lists with a specific value which shall not change
             * -> prohibit update of list
             */
            if(isTesting_ == false){
              /*
               * Update list of rotational operators according to new input
               */
              UpdateList(rotList, X_thres, e_u, alpha, beta, idElem, i, j, true);
            }
  //
  //          std::cout << "After update rotlist" << std::endl;
  //          std::cout << "preMinMaxChanged_[idElem][idAlpha][idBeta]: " <<preMinMaxChanged_[idElem][i][j] <<std::endl;
  //          std::cout << "preRotChanged_[idElem][idAlpha][idBeta]: " <<preRotChanged_[idElem][i][j] <<std::endl;

            /*
             * Evaluate rotation state of current element
             * -> if list did not change, the function returns rotationStates_[idElem][i][j]
             */
            if(evalVersion_ == 1){
              EvaluatePreisachElement(rotList, retVal, retVec, alpha, beta, idElem, i, j, true);
            } else if(evalVersion_ == 2) {
              EvaluatePreisachElement_v2(rotList, retVal, retVec, alpha, beta, idElem, i, j, true);
            } else if(evalVersion_ == 3) {
              EvaluatePreisachElement_v3(rotList, retVal, retVec, alpha, beta, idElem, i, j, true);
            } else if(evalVersion_ == 4){
               // we have to evaluate the switching state of the element in this case, too
               // otherwise, we cannot compute xPar which is needed to update the list for the switching
               // operators! Note that v4 is not possible here, as it does the linked evaluation of both lists
               EvaluatePreisachElement_v3(rotList, retVal, retVec, alpha, beta, idElem, i, j, true);
            }
  //          } else if(evalVersion_ == 5){
  //            if((abs(evaluatedState_[idElem][i][j][0])<tol_)&&(abs(evaluatedState_[idElem][i][j][1])<tol_)){
  //              // we have to evaluate the switching state of the element in this case, too
  //              // otherwise, we cannot compute xPar which is needed to update the list for the switching
  //              // operators! Note that v4 is not possible here, as it does the linked evaluation of both lists
  //              EvaluatePreisachElement_v3(rotList, retVal, retVec, alpha, beta, idElem, i, j, true);
  //            } else {
  //              /*
  //               * Test: instead of taking the rotation state, take the complete evaluation state of the element
  //               *
  //               */
  //              retVec = evaluatedState_[idElem][i][j];
  //            }
  //          }
            /*
             * Store current rotation state
             */

            rotationStates_[idElem][i][j] = retVec;
          }
  //        std::cout << "After evalPreisach" << std::endl;
  //        std::cout << "preMinMaxChanged_[idElem][idAlpha][idBeta]: " <<preMinMaxChanged_[idElem][i][j] <<std::endl;
  //        std::cout << "preRotChanged_[idElem][idAlpha][idBeta]: " <<preRotChanged_[idElem][i][j] <<std::endl;
  //

          //std::cout << "Resulting rotation state value: " << retVec.ToString() << std::endl;

          /*
           * SWITCHING PART
           */
          /*
           * Calculate parallel component of input
           *
           * This is the crucial part where the problem arises.
           * Although we are able to calculate the element value without any averaging (version 4)
           * the input for the minMaxList always results from such an averaging process.
           * In that way, simulations with finer Preisach planes (assuming constant weight of 0.5) return
           * other results than simulation with coarser planes. This is not the case, if the rotation is fixed
           * during the simulation as in that case, each element just has one rotation state which is equal to
           * the average state.
           * To solve this averaging problem, we would have to compute xPar for each rotation substate of the element.
           * This would be possible, but would require some sort of xPar-memory for each substate. From this memory
           * one could calculate for each substate the exact weighting, resulting in the supposingly correct result.
           * Problem: What happens to the memory if one substate gets overwritten? In the original version, assuming a very
           * fine discretization of the Preisach plane, the xPar states/switching states would be kept as they are indep-
           * endent of the rotation state. In the case of xPar-memory for each substate, the states are not stored
           * inpependetly of each other anymore. A deletion is such not easily done.
           *
           */

          Vector<Double> rotNormalized = Vector<Double>(dim_);

          if(rotationStates_[idElem][i][j].NormL2() > tol_){
            rotNormalized = rotationStates_[idElem][i][j]/rotationStates_[idElem][i][j].NormL2();
          } else {
            rotNormalized.Init(0.0);
          }

          //Double xParstd = xVal.Inner(rotationStates_[idElem][i][j]);

          /*
           * NEW 11.5.2016
           * Use normalized rotation state to determine xPar
           * Motivation:
           *  If an element is only partially switched (e.g. in the beginning, where only a smaller part of the
           *  element has rotstate e_u and the rest has rotstate 0), the final value will be decreased twice
           *  1. the average rotstate is only a fraction of the e_u
           *  2. the switching state which is determined by xPar will be lower as xPar is decreased
           *
           *  As long as xPar is still large enough to set the corresponding rotation area to +1/-1 this is no problem
           *  as the double reduction gets rescaled by the multiplication with the area. One can show that a doubling of
           *  numRows leads to
           *  1. roughly half xPar, roughly half rotstate
           *  2. four time area
           *  -> in total the result matches independently from the discretization, but only in the aformentioned case that
           *  the rotstate is e_u to some fraction and 0 otherwise. In cases that the other part has a different state, this
           *  does not hold in general. Furthermore, it only holds if the filling state is already very low. If this is not
           *  the case, the factor two does not work as it obeys (2^(n+1)-1)/(4^n) (e.g. 1>3/4>7/16>...)
           *  In tha case that we have some elements which are more densly filled (n -> 0), tests showed, that the results
           *  fit much better to each other (especially comparing an odd and an even number or rows) if we use a normalized rotation
           *  state instead of the original one, but only to determine xPar
           *
           *  To the testcase:
           *  for small input and odd number of rows, xthres sets only half of the diagonal element with center 0,0
           *  therewith, rotation state is e_u/2 and xPar is half the value it could reach
           *  -> using the normalized rotation state will just cancel out the factor 2
           *  using even number of rows, this problem does not occur as 0,0 does not lie in the center of element, but
           *  at a corner
           *
           *  Conclusion:
           *  normalization does work, but it is not clear if this is usable in general
           *  it is not clear yet if normalization is allowed at all for the calculation of xPar
           *
           *  TODO:
           *  check whats going on with version 4 -> works great for odd and even number of rows!
           *  Howerver, it just works if the normalized rotation state is used!!!! otherwise it shows the
           *  same drop in amplitude as version 3
           *
           *  UPDATE 12.5.2016
           *  Apparently the good results were more or less the consequence of the rotation vector being quite fix (only changing sign)
           *  which came from a very simple setup; already the next complexer setup (using four elements instead of one with hysteresis)
           *  shows huge differences between odd and even and between coarse and fine discretization! however, this differences
           *  only occur for decreasing hystersis state
           *  -> reason: at the beginning, when the rotstate is still 0, adding a rotstate to an element and then normalizing
           *  its length to 1 will perfectly work with the stageless switching. by the normalizing process, the whole element will
           *  act as thought it has the same rotstate (as the rest is 0) and such the overlap with the stageless switching states
           *  will result in correct weighting areas for the rotation; note that in this upswitching phase even the setting of
           *  an rotation element to e_u would have the same effect.
           *  However, during the downswitching process, the remainder of the elements is not 0 anymore! Therefore changes in the
           *  rotation state are only then recognized by the program if it is above a certain value. This value is reached if
           *  the calculation of xPar results in a value large enough to set the switching state. Then the change in the rotation state
           *  has some influence to the result (even thought it not has to be the correct one)
           *  -> CONCLUSION: during upswitching, v3/6 and v4/5 work perfectly; during downswitching, they work good (and much better than
           *  the original version v0) but are bounded to the discretization (even though this effect might barely be noticable in practice)
           *
           *  UPDATE 12.5.2016 - 2
           *  Finally, it seems like a global solution without any averaging is possible. Until now the problem was to relate the switching
           *  state to the rotation state.
           *  Suggested concept:
           *    1. Starting phase/first input
           *    1.1 start with one global rotation list; update this list using the xthresh as key and e_u as value
           *    1.2 for the first entry of the global rotation list, create an empty switching list and add xPar = x_in * e_u to it
           *    1.3 overlap rotation state and switching state like it is done in EvaluatePreisachElement_v4
           *    1.4 cut rotation state against matrix of size of the preisach weights, then evaluate the state like usual
           *
           *    2. Following input
           *    2.1 calculate xthresh and check if new rotation state overlaps with previous entries of the list
           *    2.1.1 for each entry which is completely overlapped
           *          -> keep entry including switching list (not allowed to lose history) but change the stored vector to the new one
           *          for each entry which is partially overlapped (can only be one at a time)
           *          -> add new entry to the list; this entry inherits the switching list from the overlapped area
           *    2.1.2 try to merge entries
           *          -> two entries can be merged, if rotationstate and the switching list are identical
           *    2.2 for each entry in rotlist, update the switching list with a new value of xPar (this one is not averaged but
           *          really the parallel component!)
           *    2.3 in the same loop, overlap rotation state and switching state with a function like EvaluatePreisachElement_v4
           *    2.4 cut rotation state against matrix of size of the preisach weights, then evaluate the state like usual
           *
           *  Idea for overlapping function of switching state with preisach weights:
           *    -> as currently, the overlapping algorithm of rotation state and switching state will be based on a
           *    decomposition of each of this states into plain rectangles (note that each rotation state will be decomposed
           *    into at max 2 rectangles in the new version!); the inclusion of the preisach weights could be done in the
           *    following way: for each calculated overlapping area, use the known upper, lower, right and left boundaries
           *    to determine a subpart of the weighting matrix; add up all weights of this subpart (elements which are only
           *    partially overlapped will only be added partially); this sum is then multiplied with the switching value (+1 or -1)
           *    and then applied as weighting to the rotation state; add this weighted rotation state to a summed up value;
           *    after all rotation states where handeled, this summed-up-value will store the result of the hysteresis operator
           *
           *  Idea to save runtime:
           *    for each rotation area include a flag isUpdated and a value representing the completely evaluated state, i.e.
           *    the value resulting from the rotation vector multiplied with switching areas and preisach weights
           *    if isUpdate is false, simple resue the old value
           *    isUpdate gets set to true, if
           *    1. rotation state got merged with another one
           *    2. if rotation state got overlapped by another one
           *    3. if the rotation state changed
           *    4. if the switching list got updated
           *
           */

          Double xPar = xVal.Inner(rotNormalized);
          //Double xPar = xVal.Inner(rotationStates_[idElem][i][j]);

//          if(idElem == 0){
//            // diagonal element of relevance
//            if((i == numRows_/2)&&(j == numRows_/2)){
//              std::cout << "idAlpha, idBeta:" << i << "," << j << std::endl;
//              std::cout << "alpha,xthres: " << alpha << "," << X_thres << std::endl;
//              std::cout << "beta,-xthres: " << beta << "," << -X_thres << std::endl;
//              std::cout << "e_u: " << e_u.ToString() << std::endl;
//              std::cout << "rotState: " << rotationStates_[idElem][i][j].ToString() << std::endl;
//              std::cout << "xPar (normally): " << xParstd << std::endl;
//              std::cout << "xPar (form normalized state): " << xPar << std::endl;
//            }
//          }
          /*
           * Normalize xPar by Xsaturated_
           */
          xPar /= Xsaturated_;

          /*
           * Get current list storing the history of the elements switching state
           */
          std::list<ListEntryOld>& switchList = preMinMaxValues_[idElem][i][j];

          /*
           * Evaluate switching state of current element
           * -> if list did not change, the function returns switchingStates_[idElem][i][j]
           */
          retVal = -10.0;
          if(evalVersion_ == 0){
            /*
             * classical approach -> set element state directly; no lists; only discrete values +1/-1 (except on alpha = beta)
             */
            if(xPar > alpha){
              retVal = 1.0;
              if(i == j){
                retVal = retVal/2.0;
              }
            } else if (xPar < beta){
              retVal = -1.0;
              if(i == j){
                retVal = retVal/2.0;
              }
            } else {
              retVal = switchingStates_[idElem][i][j];
            }

          } else {
            /*
             * Update list of switching states according to xPar
             * -> passed vector is not used, simply pass e_u
             */
            UpdateList(switchList, xPar, e_u, alpha, beta, idElem, i, j, false);

            if(evalVersion_ == 1){
              EvaluatePreisachElement(switchList, retVal, retVec, alpha, beta, idElem, i, j, false);
            } else if(evalVersion_ == 2){
              EvaluatePreisachElement_v2(switchList, retVal, retVec, alpha, beta, idElem, i, j, false);
            } else if(evalVersion_ == 3){
              EvaluatePreisachElement_v3(switchList, retVal, retVec, alpha, beta, idElem, i, j, false);
            } else if(evalVersion_ == 4){
              EvaluatePreisachElement_v4(rotList,switchList, retVal, retVec, alpha, beta, idElem, i, j);
            }
          }

          /*
           * Store current switching state
           */
          switchingStates_[idElem][i][j] = retVal;

          if(evalVersion_ == 4){
            /*
             * here we can directly use retVec as it contains the weighted rotation state
             */
            evaluatedState_[idElem][i][j] = retVec;
          } else {
            /*
             * DO NOT USE the old retVec directly as it was overwritten during the evaluation of the switching operators!
             */
            evaluatedState_[idElem][i][j] = rotationStates_[idElem][i][j]*retVal;
          }

          /*
           * EVALUATION PART
           */
          Double weight = preisachWeights_[i][j];

          Yout += evaluatedState_[idElem][i][j]*weight;

//          if(idElem == 0){
//            // diagonal element of relevance
//            if((i == numRows_/2)&&(j == numRows_/2)){
//              std::cout << "switchState: " << switchingStates_[idElem][i][j] << std::endl;
//              std::cout << "evalState: " << evaluatedState_[idElem][i][j].ToString() << std::endl;
//              std::cout << "area: " << area << std::endl;
//              Vector<Double> tmp = evaluatedState_[idElem][i][j]*area;
//              std::cout << "weighted state: " << tmp.ToString() << std::endl;
//            }
//          }
        }
      }
    }

    if(idElem == 0){

//    std::cout << "Switching States: " << std::endl;
//
//    for (UInt i = 0; i < numRows_; i++){
//
//      for(UInt j = 0; j <= i; j++){
//        std::cout << std::scientific << std::setprecision(12) << switchingStates_[idElem][i][j] << " ";
//      }
//      std::cout << std::endl;
//    }


//    std::cout << "Rotation States: " << std::endl;
//    std::cout << "x: " << std::endl;
//    for (UInt i = 0; i < numRows_; i++){
//
//      for (UInt j = 0; j <= i; j++){
//        std::cout << std::scientific << rotationStates_[idElem][i][j][0] << " ";
//      }
//      std::cout << std::endl;
//    }
//
//    std::cout << "y: " << std::endl;
//    for (UInt i = 0; i < numRows_; i++){
//
//      for (UInt j = 0; j <= i; j++){
//        std::cout << std::scientific << rotationStates_[idElem][i][j][1] << " ";
//      }
//      std::cout << std::endl;
//    }
//


//    std::cout << "Evaluated States: " << std::endl;
//    std::cout << "x: " << std::endl;
//    for (UInt i = 0; i < numRows_; i++){
//
//      for (UInt j = 0; j <= i; j++){
//        std::cout << std::scientific << std::setprecision(12) << evaluatedState_[idElem][i][j][0] << " ";
//      }
//      std::cout << std::endl;
//    }
//
//    std::cout << "y: " << std::endl;
//    for (UInt i = 0; i < numRows_; i++){
//
//      for (UInt j = 0; j <= i; j++){
//        std::cout << std::scientific << std::setprecision(12) << evaluatedState_[idElem][i][j][1] << " ";
//      }
//      std::cout << std::endl;
//    }
//    std::cout << "Yout: " << std::endl;
//    std::cout << std::scientific << std::setprecision(16) << Yout[0]*(area) << std::endl;
//    std::cout << std::scientific << std::setprecision(16) << Yout[1]*(area) << std::endl;
    }

    preisachSum_[idElem] = Yout*(area*YSaturated_);

    return preisachSum_[idElem];
  }

  Vector<Double> VectorPreisach::getValue_vec( Integer idElem )
  {
    return ( preisachSum_[idElem] );
  }

  Double VectorPreisach::getRectangleBounds(std::list<ListEntryOld>& list,std::list<ListEntryOld>::iterator startIt, std::list<ListEntryOld>::iterator endIt,
                                          Double alpha, Double beta, UInt idElem, UInt idAlpha, UInt idBeta, UInt idArea,
                                          bool isTriangle, bool isRot,
                                          Double& tRet, Double& bRet, Double& lRet, Double& rRet){
    /*
     * Calculates the coordinates of a rectangle lying inside element alpha,beta
     *
     * Input:
     *  list -> minmaxlist or rotlist of the current element
     *  startIt -> iterator pointing to the current entry
     *  endIt -> iterator pointing to the last entry
     *  alpha,beta -> bottom left corner coordinates of element
     *  idElem -> index of FE element
     *  idAlpha,idBeta -> indices of the element
     *  idArea -> index of area to be computed - 1 (i.e. area1 -> idArea = 0)
     *  isTriangle -> flag checking if element is on diagonal alpha = beta
     *  isRot -> rotation list or switching list?
     *
     *  tRet,bRet -> y-coordinates of the top and bottom edge of the rectangle
     *  lRet,rRet -> x-coordinates of the left and right edge of the rectangle
     */

    Double top = 0.0;
    Double bot = 0.0;
    Double left = 0.0;
    Double right = 0.0;

    Double entry = startIt->getVal();
    Double returnVal = 1.0; //normally all areas returned shall be positive; +/- will be applied later
    //exception: for overlapping areas (area1 overlapping area3, see notes below) we need a flipped sign
    //-> returnVal = -1.0

    //Vector<Double> entryVec = startIt->getVec();
    bool isMinEntry = startIt->isMin();

    if(isMinEntry && (idArea == 0)){
      /*
       * First entry in list is a minimum; two possibilities:
       * 1. element is not wiped out yet -> minimum lies between beta and beta+delta
       * 2. element is wiped out by minimum or got initialized or overwritten with minimum -> minimum has value beta
       *
       * for case 2, area 1 will be 0 as left = right
       *
       * Note that the area 2 is already a std area, so instead of writing it once again into this case, we just
       * execute it additionally to this one here.
       */

      /*
       * area 1
       */
      left = beta;
      right = entry;

      /*
       * Note to the commented out part below:
       *
       * The idea behind this code segment was the following.
       * Consider a min/max staircase pattern of the form:
       *
      *   MinMax[0] = min  MinMax[2] = min
      *       ____v________v___
      *      | \  |   |    |   |
      *      |   \|   |    |   |
      *      |    |\  |    |   |
      *      |____|__\|____|___| MinMax[1] = max
      *      |    |   |\   |   |
      *      |____|___|__\_|___| MinMax[3] = max
      *      |    |   |    \   |
      *      |____|___|____| \_|
      *              MinMax[4] = min
       *
      *
      *       quadratic part which is still split by diagonal -> only for elements on alpha = -beta
      *       (note that upper and lower part cancel out in sum, so that we can omit this part)
      *       _v_______________
      *      | \  | 2 |    |   |
      *      |___\|   | 4  |   |
      *      |   1|   |    |6  |
      *      |____|___|    |   |
      *      |     3  |    |   |
      *      |________|____|   |
      *      |        5    |   |
      *      |_____________|___|
      *
      *      (area x: width x height)
      *      area 1: (MinMax[0] - beta) x (alpha* - MinMax[1])
      *        with  alpha* = alpha+delta - MinMax[0] + beta for elements on alpha = -beta
      *              alpha* = alpha+delta_                   for elements not on alpha = -beta
      *      area 2: (MinMax[2] - MinMax[0]) x (alpha+delta - MinMax[1])
      *      area 3: (MinMax[2] - beta) x (MinMax[1] - MinMax[3])
      *      area 4: (MinMax[4] - MinMax[2]) x (alpha+delta - MinMax[3])
      *      area 5: (MinMax[4] - beta) x (MinMax[3] - alpha)
      *      area 6: (beta+delta - MinMax[4]) x (alpha+delta - alpha)
      *
      *      For elements on alpha = -beta, the upper part (area without number) is split by diagonal.
      *      -> its value is half +1 and half -1 and thus can be neglected in sum
      *      -> ONLY IF WEIGHTS ARE SYMMETRIC TO alpha = -beta!
      *      --> not relevant for these older version, as we evaluate element based, i.e. the weight
      *      has to be symmetric over the element, which is the case as it is constant in each element
      *      of the Preisach plane!
      *
      *      To omit the computation of that area, we use alpha* = alpha+delta - MinMax[0] + beta as
      *      upper boundary of area 1.
      *      Once the element got wiped out (i.e. min of value beta or max of value alpha+beta), the
      *      upper area is no longer split and thus needs to be considered. However, assuming that
      *      MinMax[0] is beta in the upper case, area 1 and the split area will be zero as
      *      MinMax[0] - beta = 0. The computation will then start with area 2, which not only encom-
      *      passes the formerly split part but furthermore has the right sign (-1) in aboves case.
      *      -> Conclusion: for elements on alpha = -beta we do NOT have to change the calculation of
      *      alpha* after a wipe out, as its value will be irrelevant then (however, changing the
      *      computation will have no bad effect either).
      *
      *      For elements not on alpha = -beta, the upper part was never split. For these elements, we
      *      are not allowed to leave out the upper part, thus alpha* should be alpha+delta. Luckily,
      *      all elements which are not on alpha = -beta already got initialized with a min of value beta
      *      or a max of value alpha+delta_. Due to this initial values, area 1 will be 0 either way and
      *      thus the value of alpha* is irrelevant again.
      *      -> Conclusion: for elements not on alpha = -beta, we do not need alpha* at all.
      *
      *      Why do we distinguish between Triangle and Square elements when checking for alpha*?
      *      -> no idea anymore. It should follow the same rules as for the square.
      */

//      if(isTriangle == true){
//        if( (idAlpha+idBeta+1 == numRows_)&&(wipedOut_[idElem][idBeta] == false) ){
//          /*
//           * In contrast to the square element case, we have to distinguish between alpha = -beta elements
//           * which are not wiped out and all other elements as area0 will not automatically get 0
//           */
//          top = alpha+delta_ - entry + beta;
//        } else {
//          std::cout << "HIT1" << std::endl;
//          top = alpha+delta_; //TODO: needed? -> see text above
//        }
//      } else {
//        top = alpha+delta_ - entry + beta;
//      }

      top = alpha+delta_ - entry + beta;

      if(startIt != endIt){
        startIt++;
        if(startIt->isMin() == isMinEntry){
          /*
           * element of same type detected -> ok for rotList, not ok for MinMaxList
           */
          if(isRot){
            /*
             * set bot to alpha
             */
            bot = alpha;
          } else {
            EXCEPTION("MinMaxList has to be alternating!")
          }
        } else {
          bot = startIt->getVal();
        }
        startIt--;
      } else {
        bot = alpha;
      }
      
//      std::cout << "top,bot: " << top << "," << bot << std::endl;
//      std::cout << "left,right: " << left << "," << right << std::endl;

      if(top < bot){
    //    std::cout << "top < bot" << std::endl;
        /*
         * NEW 10.5.2016
         *
         * check if top is smaller than bot!
         * if yes, we have the case, that the helper area overlaps with the next area (see comment under TODO, OBSERVATION)
         *
         * Note:
         * this case can only occur, if minmaxlist[1] (i.e. the second entry in the list) is smaller than the value of top
         * this can in general only happen on line alpha = -beta (other elements will always have a helper area of 0)
         * In our case: beta = -(alpha+delta) as we use the coordinates of the bottom left corner
         * -> top = -entry
         *
         * Original area:
         *  left = beta
         *  right = entry
         *  top = -entry
         *  bot = nextentry
         *
         * Correct area:
         *  left = -nextentry
         *  right = entry
         *  top = alpha+delta
         *  bot = nextentry
         *
         * Note2:
         * the resulting area must have an opposite sign compared to the normal sign as it is not added to the negative part
         * but to the positive one!
         * -> idea: switch sings of bot and top, such that the area becomes negative after this calculation
         * -> do not forget to allow negative updates!
         *
         * -> does not work as we later cut down parts which are below alpha = beta and this will trigger
         * if we multiply with -1 already at this points
         *
         * -> other idea:
         * let function return an additional factor of +/-1
         *
         */
        returnVal = -1;
        left = -bot; // bot already holds the value of nextentry
        right = entry;

        top = (alpha+delta_);
        //bot = bot; // //-1 for creating a negative area!
//        std::cout << "top,bot: " << top << "," << bot << std::endl;
//        std::cout << "left,right: " << left << "," << right << std::endl;
      }

    } else if(!isMinEntry && (idArea == 0)) {
      /*
       * First entry in list is a maximum; two possibilities:
       * 1. element is not wiped out yet -> maximum lies between alpha and alpha+delta
       * 2. element is wiped out by maximum or got initialized or overwritten with maximum -> maximum has value alpha+delta
       *
       * for case 2, area 1 will be 0 as top = bot
       *
       * Note that the area 2 is already a std area, so instead of writing it once again into this case, we just
       * execute it additionally to this one here.
       */
      /*
       * area 1
       */
      if(startIt != endIt){
        startIt++;
        if(startIt->isMin() == isMinEntry){
          /*
           * element of same type detected -> ok for rotList, not ok for MinMaxList
           */
          if(isRot){
            /*
             * set bottom to beta+delta_
             */
            right = beta+delta_;
          } else {
            EXCEPTION("MinMaxList has to be alternating!")
          }
        } else {
          right = startIt->getVal();
        }
        startIt--;
      } else {
        right = beta+delta_;
      }

      // This is a relict from older version; due to the new calculation scheme, the helper areas will vanish
      // in case of triangular elements, too
      // TODO: does wipedOut still has any effect?

//      if(isTriangle == true){
//        if( (idAlpha+idBeta+1 == numRows_)&&(wipedOut_[idElem][idBeta] == false) ){
//          /*
//           * In contrast to the square element case, we have to distinguish between alpha = -beta elements
//           * which are not wiped out and all other elements as area0 will not automatically get 0
//           */
//          left = alpha+delta_ - entry + beta;
//        } else {
//          std::cout << "HIT2" << std::endl;
//          left = beta;
//        }
//      } else {
//        left = alpha+delta_ - entry + beta;
//      }

      left = alpha+delta_ - entry + beta;

      top = alpha+delta_;
      bot = entry;

//      std::cout << "top,bot: " << top << "," << bot << std::endl;
//      std::cout << "left,right: " << left << "," << right << std::endl;

      if(right < left){
  //      std::cout << "right < left" << std::endl;
        /*
         * NEW 10.5.2016
         *
         * check if right is smaller than left!
         * if yes, we have the case, that the helper area overlaps with the next area (see comment under TODO, OBSERVATION)
         *
         * Note:
         * this case can only occur, if minmaxlist[1] (i.e. the second entry in the list) is smaller than the value of left
         * this can in general only happen on line alpha = -beta (other elements will always have a helper area of 0)
         * In our case: beta = -(alpha+delta) as we use the coordinates of the bottom left corner
         * -> left = -entry
         *
         * Original area:
         *  left = -entry
         *  right = nextentry
         *  top = alpha+delta
         *  bot = entry
         *
         * Correct area:
         *  left = beta
         *  right = nextentry
         *  top = -nextentry (mirrored at alpha = -beta
         *  bot = entry
         *
         * Note2:
         * the resulting area must have an opposite sign compared to the normal sign as it is not added to the negative part
         * but to the positive one!
         * -> idea: switch sings of left and right, such that the area becomes negative after this calculation
         * -> do not forget to allow negative updates!
         *
         * -> does not work as we later cut down parts which are below alpha = beta and this will trigger
         * if we multiply with -1 already at this points
         *
         * -> other idea:
         * let function return an additional factor of +/-1
         *
         */
        returnVal = -1;
        top = -right; // right already holds the value of nextentry
        bot = entry;

        left = beta; //-1 for creating a negative area!
        //right = right; //-1 for creating a negative area!
//        std::cout << "top,bot: " << top << "," << bot << std::endl;
//        std::cout << "left,right: " << left << "," << right << std::endl;
      }

    } else if(isMinEntry && (idArea > 0)) {
      /*
       * Current entry is negative
       */
      right = beta+delta_;
      bool found = false;

      if(startIt != endIt){
        startIt++;
        if(startIt->isMin() == isMinEntry){
          /*
           * element of same type detected -> ok for rotList, not ok for MinMaxList
           */
          if(isRot){
            /*
             * set bottom to alpha
             */
            bot = alpha;
            /*
             * here we can directly set the next entry, too
             * right is determined by the consecutive minimum, which just got found here
             */
            found = true;
            right = startIt->getVal();
          } else {
            EXCEPTION("MinMaxList has to be alternating!")
          }
        } else {
          bot = startIt->getVal();
        }

        if(startIt != endIt){
          startIt++;
          if(startIt->isMin() != isMinEntry){
            EXCEPTION("MinMaxList has to be alternating!")
            /*
             * here we would expect the same entry type the the one of the current iteration
             * e.g. a minimum after a maximum which followd a minimum
             */
          } else {
            if(isRot&&found){
              /*
               * do nothing here -> we need the second min/max here but this value was already found in the previous entry
               *
               * in case of alternating rotList, we retrieved just the value of bot in the last case, so we need to find right
               */
            } else {
              right = startIt->getVal();
            }
          }
          startIt--;
        }
        startIt--;
      } else {
        bot = alpha;
      }

      top = alpha+delta_;
      left = entry;

    } else {
      /*
       * Current entry is positive
       */
      bot = alpha;
      bool found = false;

      if(startIt != endIt){
        startIt++;
        if(startIt->isMin() == isMinEntry){
          /*
           * element of same type detected -> ok for rotList, not ok for MinMaxList
           */
          if(isRot){
            /*
             * set right to beta+delta_
             */
            right = beta+delta_;
            /*
             * here we can directly set the next entry, too
             * right is determined by the consecutive minimum, which just got found here
             */
            found = true;
            bot = startIt->getVal();
          } else {
            EXCEPTION("MinMaxList has to be alternating!")
          }
        } else {
          right = startIt->getVal();
        }

        if(startIt != endIt){
          startIt++;
          if(startIt->isMin() != isMinEntry){
            EXCEPTION("MinMaxList has to be alternating!")
            /*
             * here we would expect the same entry type the the one of the current iteration
             * e.g. a minimum after a maximum which followd a minimum
             */
          } else {
            if(isRot && found){
              /*
               * do nothing here -> we need the second min/max here but this value was already found in the previous entry
               *
               * in case of alternating rotList, we retrieved just the value of right in the last case, so we need to find bot
               */
            } else {
              bot = startIt->getVal();
            }
          }
          startIt--;
        }
        startIt--;
      } else {
        right =  beta+delta_;
      }

      top = entry;
      left = beta;
    }

    if(isTriangle == true){
      /*
       * restrict bound such, that only the bottom right corner lies below the diagonal
       */
      if(right-beta > top-alpha){
      /*
       * upper right corner below diagonal
       * -> set right bound to top
       */
      right = top;
      }
      if(bot-alpha < left-beta){
      /*
       * lower left corner below diagonal
       * -> set bottom bound to left
       */
      bot = left;
      }
    }

    tRet = top;
    bRet = bot;
    lRet = left;
    rRet = right;

    return returnVal;
  }

  bool VectorPreisach::clipRectangles(Double t1, Double b1, Double l1, Double r1, Double t2, Double b2, Double l2, Double r2, Double& tRet, Double& bRet, Double& lRet, Double& rRet){
    /*
     * Calculates the coordinates of the overlapping area
     * resulting from the overlap of two rectangles defined by their top, bottom, left and right edges
     *
     * Input:
     *        t1,b1 -> y-coordinates of the top and bottom edge of the first rectangle
     *        l1,r1 -> x-coordinates of the left and right edge of the first rectangle
     *        t2,b2 -> y-coordinates of the top and bottom edge of the second rectangle
     *        l2,r2 -> x-coordinates of the left and right edge of the second rectangle
     *
     *        tRet,bRet -> y-coordinates of the top and bottom edge of the overlapping area
     *        lRet,rRet -> x-coordinates of the left and right edge of the overlapping area
     *
     * Return:
     *        true if overlap exists
     *        false if overlap does not exist or is equals 0
     *
     */

    /*
     * check if overlap does exist
     */
    if( (t1 <= b2)||(t2 <= b1) ) {
      /*
       * top edge of one rectangle is smaller lies below or on top of the bottom edge of the other rectangle
       * -> no overlap area possible
       */
      return false;
    }
    if( (r1 <= l2)||(r2 <= l1) ){
      /*
       * right edge lies left or on top of the left edge of the other element
       * -> no overlap area possible
       */
      return false;
    }

    tRet = std::min(t1,t2);
    bRet = std::max(b1,b2);
    lRet = std::max(l1,l2);
    rRet = std::min(r1,r2);

    return true;

  }

  void VectorPreisach::switchingStateToBmp(UInt numPixel, std::string filename, UInt idElem, bool overLayWithRotState)
  {
    /*
     * directly call output function of matrix
     */
    UInt upscaling = floor(numPixel/numRows_);
    UInt realNumPixel = upscaling * numRows_;
    if(numPixel != realNumPixel){
      WARN("Image scaled to "<<realNumPixel<<" x "<<realNumPixel);
    }

    if(realNumPixel < 2){
      WARN("Image should have more than 2 x 2 pixel");
      return;
    }

    if(realNumPixel%2 != 0){
      WARN("Rounded number of pixel ("<<realNumPixel<<") to a multiple of 2 ("<<realNumPixel+1<<")");
      realNumPixel = realNumPixel + 1;
    }

    if(overLayWithRotState == true){
      /*
       * we have to transfer the rotation state from a vector*** to a matrix* first
       */
      Matrix<Double>* rotStorage = new Matrix<Double>(numRows_,numRows_);
      rotStorage->Init();

      for(UInt comp = 0; comp < dim_; comp++){

        /*
         * upper part
         */
        for ( UInt y=numRows_/2; y<numRows_; y++ )
        {
          for ( UInt x=0; x<=y; x++ )
          {
            (*rotStorage)[y][x] = rotationStates_[idElem][y][x][comp];
          }
        }
        /*
         * lower triangle -> these values are not set during calculation and therefore have
         * not the current rotation state
         * -> get rotation state of upper triangle
         */
        for ( UInt y=0; y<numRows_/2; y++ )
        {
          for ( UInt x=0; x<=y; x++ )
          {
            (*rotStorage)[y][x] = rotationStates_[idElem][numRows_-1][numRows_-1][comp];
          }
        }

        std::stringstream stream;
        if(comp == 0){
          stream << "x-" << filename;
        } else if(comp == 1){
          stream << "y-" << filename;
        } else {
          stream << "z-" << filename;
        }

        std::string filename_new = stream.str();

        switchingStates_[idElem].matrix2Bmp(upscaling,filename_new,rotStorage);
      }

      delete rotStorage;

    } else {
      switchingStates_[idElem].matrix2Bmp(upscaling,filename,NULL);
    }

  }

}//end namespace
