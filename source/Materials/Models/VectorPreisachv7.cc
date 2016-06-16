#include "VectorPreisachv7.hh"

#include <iostream>
#include <fstream>

/* See VectorPreisachv7.hh for detailed description */

namespace CoupledField
{ 

  VectorPreisachv7::VectorPreisachv7(Integer numElem, Double xSat, Double ySat,
     Matrix<Double>& preisachWeight, Double rotationalResistance , UInt dim, bool isVirgin, bool isTesting, UInt evalVersion)
     : Hysteresis(numElem)
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

      tol_ = 1e-16;

      /*
       * IMPORTANT NOTE:
       * In the new version, we (currently) require the preisachWeights_ to be symmetric to alpha = -beta.
       * The reason for this lies in the computation of the upper square region. We assume an initial state
       * of the switching plane which is +1 below alpha=-beta and -1 above. Therewith, we do not have to
       * evaluate parts of the square area which are still split by the diagonal alpha = -beta (see comments
       * to function getRectangleBounds). If the underlying weights are not symmetric, however, we have
       * to evaluate that split area, too.
       * -> TODO: implement that
       * Idea: add helper area above area1. That area is split into an upper and a lower triangle wrt alpha = -beta.
       * Clip both triangles to the Preisach plane. Use flag to indicate the special treatment or write new
       * function. Make std mapping (i.e. find out which elements are partially and completely overlapped) with
       * the full square area first. For the lower triangle make the loop such that only elements below the
       * diagonal get added (excluding the elements on the diagonal) and for the upper triangle add only
       * these elements which are above the diagonal (again excluding the diagonal). Both triangles will have
       * a stair cased hypotenuse, Between both hypotenuses we have one line of Preisach elements left which
       * are not added. This is completely fine, as these elements are symmetric (as the weighting has to be
       * constant over each element) and thus positive and negative half of these elements cancel out.
       *
       * Until this is done, check symmetry and throw exception.
       *
       */
      isSymmetric_ = preisachWeights_.IsSymmetric(tol_);
      //isSymmetric_ = false;
      // done but not tested!
//      if(!isSymmetric_){
//        EXCEPTION("Currently, version7 requires the Preisach weights to be symmetric with respect to alpha = -beta")
//      }

      // get size of preisachWeights_
      UInt M = preisachWeights_.GetNumRows();
      UInt N = preisachWeights_.GetNumCols();

      if(M != N){
        EXCEPTION("Matrix preisachWeight has dim " << M << " x " << N << " and thus is not symmetric!");
      }

      numRows_ = M;

      // resolution of Preisach plane
      delta_ = 2.0/Double(numRows_);

      rotationalResistance_ = rotationalResistance;

      lowerTriangleValue_ = 0.0;
      Evaluate_LowerTriangle();

  //    NOT IMPLEMENTED YET (AND NOT NEEDED)
  //      /*
  //       * Debugging quantities
  //       */
  //      isTesting_ = isTesting; //!true;
  //
  //      if(isTesting_){
  //        std::cout << "Testcase is for eval version >5 may not match the one of the other versions." << std::endl;
  //        std::cout << "Reason: The testcase allows the lower triangle part of the rotation plane to oppose the input vector."
  //                  << "Therewith, the switching states in this lower triangle are not longer fix at +1, which is assumed in version >5."
  //                  << std::endl;
  //      }
  //      rotStateTesting_ = Vector<Double>(dim_);
  //      rotStateTesting_.Init();
  //      rotStateTesting_[1] = 1;

      evalVersion_ = 7;

      /*
       * Local quantities, i.e. arrays storing different values for each FE element
       */
      preisachSum_ = new Vector<Double>[numElem_];
      for(UInt k = 0; k < numElem_; k++){
        preisachSum_[k] = Vector<Double>(dim_);
        preisachSum_[k].Init();
      }

      upperTriangleValue_ = Vector<Double>(numElem_);
      lastXpar_ = Vector<Double>(numElem_);

      lastEu_ = new Vector<Double>[numElem_];
      for(UInt k = 0; k < numElem_; k++){
        lastEu_[k] = Vector<Double>(dim_);
        lastEu_[k].Init();
      }

      globSwitchList_ = new std::list<ListEntry>[numElem_];
      for(UInt k = 0; k < numElem_; k++){
        globSwitchList_[k] = std::list<ListEntry>();
        Initialize_GlobalSwitchingList(k);
      }

      globRotList_ = new std::list<RotListEntry>[numElem_];
      for(UInt k = 0; k < numElem_; k++){
        globRotList_[k] = std::list<RotListEntry>();
        Initialize_GlobalRotationList(k);
      }
    }

  VectorPreisachv7::~VectorPreisachv7(){
    delete[] preisachSum_;
    delete[] globSwitchList_;
    delete[] globRotList_;
    delete[] lastEu_;
  }

  void VectorPreisachv7::Initialize_GlobalSwitchingList(UInt idElem){
    /*
     * Make sure that list is empty
     */
    globSwitchList_[idElem].clear();

    /*
     * the initial state of the upper triangle is -1
     * this can be achieved by initializing the list with
     * the smallest possible minimum in that region (i.e. beta = 0)
     */
    globSwitchList_[idElem].push_back(ListEntry(0,true,false));

    lastXpar_[idElem] = 0.0;

    Evaluate_GlobalSwitchingList(idElem);
  }

  void VectorPreisachv7::Initialize_GlobalRotationList(UInt idElem){
    /*
     * Make sure that list is empty
     */
    globRotList_[idElem].clear();

    /*
     * Initialize list with a maximum of value 0 and 0-vector
     * (Note that in previous versions, we used to store a pair of min/max on each element
     * which has its center at the line alpha = -beta; this approach is no longer needed
     * as we do not have to treat single elements anymore)
     *
     * Furthermore, we have to add a list for the switching states; this list
     * gets initialized with two entries (a min and a max) so that both an starting
     * minimum and a starting maximum can be captured; to get a correct updating behavior,
     * we have to set these starting min/max to be dummy elements (i.e. set flag isDummy)
     */
    Vector<Double> zeroVec = Vector<Double>(dim_);
    zeroVec.Init(0);

    /*
     * Create new entry for rotlist; the contained switching list is empty
     */
    std::list<ListEntry> switchingList = std::list<ListEntry>();
    RotListEntry initEntry = RotListEntry(0.0, 0.0, zeroVec, switchingList, 0.0, false, false);

    /*
     * add min/max pair to list (starting with max) and set dummy flag to true
     * (both times use 0 as value as this is the value of the smallest possible maximum
     * in region 1 and also the value of the largest possible minimum)
     */
    initEntry.getListReference().push_back(ListEntry(0,false,true));
    initEntry.getListReference().push_back(ListEntry(0,true,true));

    globRotList_[idElem].push_back(initEntry);
  }

  void VectorPreisachv7::Update_GlobalRotationList(Double xThres, Vector<Double> e_u, Vector<Double> u_in, UInt idElem){
    /*
     * function for updating the global rotation list with an entry pair xThres,e_u
     * furthermore, the original input vector u_in is passed to update the switching lists
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

    if(globRotList_[idElem].empty()){
      EXCEPTION("List may not be empty at this point");
    }

    /*
     * Step 1. Update rotation list by inserting new entries (if needed!)
     */
    std::list<RotListEntry>::iterator listIt;
    /*
     *  Note that the function list::insert adds the new entry before position listIt
     * (i.e. if listIt = list.end() it will be inserted before the end)
     */
    std::list<RotListEntry>::iterator insertPos = globRotList_[idElem].end();
    std::list<RotListEntry>::iterator listStart = globRotList_[idElem].begin();
    std::list<RotListEntry>::iterator listEnd = --(globRotList_[idElem].end());

    Double curVal;
    /*
     * lowerBound is important for RotListEntries
     * the smallest possible entry (the one which is take if an entry is appended to the end of the list)
     * is 0
     */
    Double lowerBound = 0.0;
    Vector<Double> prevState;
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
    if(xThres <= 0.0){
      xThres = 0.0;
      needsInsert = false;
    }

    for(listIt = globRotList_[idElem].begin(); listIt != globRotList_[idElem].end(); listIt++){

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

            Double tol = 1e-15;

            needsInsert = false;
            /*
             * check if previous state has the same (or nearly) the same rotation state
             * in that case, we do not have to insert this entry, as it would be merged with
             * the previous one later one
             */
            for(UInt i = 0; i< e_u.GetSize();i++){
              if(abs(e_u[i]-prevState[i])>tol){
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
            /*
             * the new rotation area will be between xThresh and curVal
             */
            lowerBound = curVal;
          }
        }
      } // xThres >= curVal
      if(posFound == true){
        /*
         * we already have found a suitable position, i.e. all following entries would be
         * overwritten by the new entry; as we want to keep the switching list, we simply overwrite
         * the rotation state
         * (if the rotation state is already e_u, the flag rotHasChanged_ is not set to true!)
         */
        listIt->setVec(e_u);
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

      Double tol = 1e-15;

      needsInsert = false;

      for(UInt i = 0; i< e_u.GetSize();i++){
        if(abs(e_u[i]-prevState[i])>tol){
          needsInsert = true;
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
       */

      std::list<ListEntry> newList = std::list<ListEntry>();
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
        /*
         * add min/max pair to list (starting with max) and set dummy flag to true
         * (both times use 0 as value as this is the value of the smallest possible maximum
         * in region 1 and also the value of the largest possible minimum)
         */
        newList.push_back(ListEntry(0,false,true));
        newList.push_back(ListEntry(0,true,true));
        wasWipedOut = false;
        startCnt = 0;
      }

      globRotList_[idElem].insert(insertPos,RotListEntry(xThres,lowerBound,e_u,newList,lastXpar,false,false,wasWipedOut,startCnt));
    }

    /*
     * Step 2. For each entry in rotation list, update the switching list with the value xPar
     */
    Double xPar;
    Double lowerVal;
    UInt updated;
    Vector<Double> rotState;

    for(listIt = globRotList_[idElem].begin(); listIt != globRotList_[idElem].end(); listIt++){
      rotState = listIt->getVecReference();
      lowerVal = listIt->getLowerVal();

      xPar = rotState.Inner(u_in);
      /*
       * Normalize to Xsaturated
       */
      xPar /= Xsaturated_;

      /*
       * check if value should be inserted into list or not
       * Issue:
       *  If the input signal decreases gradually and changes its rotation state
       *  (e.g. linear decreasing sawtooth function)
       *  both the rotation list and the switching list will be elongated each time.
       *  Merging is not possible, as the rotation list will have alternating rotation states.
       *  Reuse of old values is not possible either, as the switching list changes and such
       *  a reevaluation is triggered.
       *  -> Computational cost increases drastically!
       * Suggestion:
       *  Check if update of switching list is needed or not. An update is not necessary if the
       *  value will not lead to a change of the switching state. This is the case if xPar is
       *  smaller than the lowerVaö of the rotation element (for an increasing xPar) or
       *  if xPar is larger than -lowerVal (for decreasing input).
       *  In these cases, the resulting switching area will have no overlap with the rotation area
       *  and such cannot change the switching state. By avoiding the addition of these entries,
       *  the switching lists will only elongate for rotation entries near the end of the list. The
       *  other list entries stay unchanged and can be reused.
       * Note:
       *  As we no longer insert the same entries to all rotation areas with the same rotation state
       *  we reduce the chance for a merge. However, these cases could not merged before either as
       *  the rotation states were not the same at neighboring elements (otherwise the areas with
       *  lower upper boundaries would not have been inserted at all). The main purpose for merges
       *  is still fulfilled which is to merge rotation states together which were overwritten by an
       *  increasing input. In that case the rotation states will be equal from a certain point on and
       *  if xPar is large enough, the switching lists will be wiped out so that a merge is valid.
       */

      /*
       * we need to pass lastXpar by
       */
      updated = Update_SwitchingList(listIt->getListReference(),xPar,listIt->getLastLocalXpar(),false, listIt->wasListWipedOut(),lowerVal);

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
         * (for exmplanation of startCnt see class RotListEntry)
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

    /*
     * Step 3. Merge adjacent rotList entries
     */
    Simplify_GlobalRotationList(idElem);

    /*
     * Step 4. Simplify local switching lsits
     */
    Simplify_LocalSwitchingLists(idElem);

//    std::cout << "GlobalRotationList after simplify" << std::endl;
//
//    for(listIt = globRotList_[idElem].begin(); listIt != globRotList_[idElem].end(); listIt++){
//    std::cout << listIt->ToString() << std::endl;
//    }
//    std::cout << "##########################" << std::endl;
  }

  UInt VectorPreisachv7::Update_SwitchingList(std::list<ListEntry>& list, Double newEntry, Double lastXpar, bool isGlobal, bool wasWipedOut, Double lowerVal){
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
     */

    Double alpha,beta,delta;
    if(isGlobal == true){
      /*
       * upper triangle
       */
      alpha = 0.0;
      beta = 0.0;
    } else {
      alpha = 0.0;
      beta = -1.0;
    }
    delta = 1.0;

    std::list<ListEntry>::iterator listIt;

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
     * Check if lists gets completely overwritten (i.e. the saturated input value exceeds the bounds of the region
     */
    if(state > 0){
      /*
       * Element is to be set by maximum -> restrict to [alpha,alpha+delta]
       */
      if(newEntry >= alpha+delta){
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
        newEntry = alpha+delta;
        /*
         * add maximum to list
         */
        list.push_back(ListEntry(newEntry,false));
        /*
         * elements on alpha = -beta also loose their special state
         */
        if(isGlobal){
          /*
           * Global list has nothing to wipe out
           */
          return 1;
        } else {
          /*
           * return 2 -> list was updated and wiped
           */
          return 2;
        }
      } else if(newEntry <= alpha){
        /*
         * new entry will have no effect
         */
        return 0;
      } else if(appendDirectly){
        list.push_back(ListEntry(newEntry,false));
        /*
         * return code 3: list was not wiped, but first entry was set/replaced
         */
        return 3;
      } else if(!isGlobal){
        /*
         * for upper square, check value of lowerVal
         * if a maximum to be added has a value smaller than lowerVal, the resulting area will have
         * no effect at all onto the switching state
         * -> do not add to list
         * (similar to the check above with newEntry <= alpha)
         */
        if(newEntry <= lowerVal){
          return 0;
        }
      }
      /*
       * else -> see below
       */
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
        list.push_back(ListEntry(newEntry,true));
        /*
         * elements on alpha = -beta also loose their special state
         */
        if(isGlobal){
          /*
           * Global list has nothing to wipe out
           */
          return 1;
        } else {
          /*
           * return 2 -> list was updated and wiped
           */
          return 2;
        }
      } else if(newEntry >= beta+delta){
        /*
         * new entry will have no effect
         */
        return 0;
      } else if(appendDirectly){
        list.push_back(ListEntry(newEntry,true));
        /*
         * return code 3: list was not wiped, but first entry was set/replaced
         */
        return 3;
      } else if(!isGlobal){
        /*
         * for upper square, check value of lowerVal
         * if a minimum to be added has a value larger than -lowerVal, the resulting area will have
         * no effect at all onto the switching state
         * -> do not add to list
         * (similar to the check above with newEntry >= beta+delta)
         * NOTE: we check for -lowerVal as lowerVal is the lower boundary of the rotation area (always >= 0)
         * and -lowerVal is the right boundary of the rotation area
         */
        if(newEntry >= -lowerVal){
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
    std::list<ListEntry>::iterator lastEntry = --(list.end());

    /*
     * if we came to this point, we have to iterate through the list and check if value shall be included
     */
    bool canBeInserted = false;
    bool canBeDeleted = false;
    bool insertHelper = false;

    if( (isGlobal == false)&&(wasWipedOut == false) ){
      /*
       * for non-wiped out elements on the diagonal alpha = -beta we have to be little more careful
       *
       *  Case 1: xPar does not replace the first entry (x0) of the list but any other value later on
       *          -> std procedure
       *  Case 2: xPar replaces the first entry (x0) of the list; list is otherwise empty
       *          -> std procedure
       *  Case 3: xPar replaces the first entry, list is not empty (second element in list shall be x1)
       *          3a) xPar is min, xPar-beta <= alpha+delta-x1
       *              -> xPar replaces x0 and is so far to the left that it wipes out rest of the list
       *              -> std procedure
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

    bool listMin;
    bool firstEntrySet = false;
    UInt cnt = 0;
    Double listVal;
    ListEntry helperEntry = ListEntry(0,false);

    for(listIt = list.begin(); listIt != list.end(); listIt++){
      /*
       * NEW: For elements on alpha = -beta (isGlobal = false), we start with a pair of min and max such that the resulting area = 0
       * However, if only the minimum is overwritten, the evaluation algorithm will return 0 due to the still remaining
       * initial maximum of value = alpha. The easiest solution is to remove the maximum, too.
       * To detect this case, check for the isDummy flag, which is only set for the initial entries.
       * If case detected -> clear list completely and start anew
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

              if(newEntry-alpha >= beta+delta-nextEntry){
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
        } else {
          /*
           * This case here should never happen, as state 3 only exists, if list starts with max
           * -> i.e. when we enter state 3 for the first iteration, we either need the special treatment
           * and set canBeDeleted to true (this breaks the loop further below)
           * or we set state to 2
           */
          EXCEPTION("This point should never be reached");
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

              if(newEntry-beta <= alpha+delta-nextEntry){
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
        } else {
          EXCEPTION("This point should never be reached");
        }
      }

      if(canBeDeleted){
        /*
         * delete all entries of the list starting with the current one
         */
        list.erase(listIt,list.end());

        /*
         * set flag denoting that the first entry of the list got deleted, too.
         * -> after this fuunction, a new first entry is to be set
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
      if(insertHelper == true){
        list.push_back(helperEntry);
      }
      if(state > 0){
        /*
         * insert maximum
         */
        list.push_back(ListEntry(newEntry,false));
      } else if(state < 0){
        /*
         * insert minimum
         */
        list.push_back(ListEntry(newEntry,true));
      }
      /*
       * mark element as updated
       */
      if((firstEntrySet == true)&&(isGlobal == false)){
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

  Vector<Double> VectorPreisachv7::computeValue_vec(Vector<Double>& xVal, Integer idElem){
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
     * set value of lastEu_ (only needed for output via rotationStateToBmp
     */
    lastEu_[idElem] = e_u;

    /*
     * Storage for return values
     */
    Vector<Double> retVec = Vector<Double>(dim_);
    retVec.Init();

    /*
     * Storage for element value
     */
    Vector<Double> Yout = Vector<Double>(dim_);
    Yout.Init(0.0);

    /*
     * Compute overall value in three steps:
     * 1. lower triangular part (alpha = -1 to 0; beta = -1 to alpha
     * 2. upper triangular part (alpha = 0 to 1; beta = 0 to alpha
     * 3. upper square part (alpha = 0 to 1; beta = -1 to 0
     */

    /*
     * Part 1 - lower triangle
     * -> lowerTriangleValue_ was calculated in constructor
     */
    Yout += e_u * lowerTriangleValue_;

    /*
     * Part 2 - upper triangle
     */
    /*
     * rotation state of upper triangle is always e_u
     */
    Double xPar = xVal.Inner(e_u);
    xPar /= Xsaturated_;

    UInt updated = Update_SwitchingList(globSwitchList_[idElem], xPar, lastXpar_[idElem], true, false);

//    std::cout << "GlobalSwitchingList after update (update == " << updated << ")" << std::endl;
//    std::list<ListEntry>::iterator listIt;
//    for(listIt = globSwitchList_[idElem].begin(); listIt != globSwitchList_[idElem].end(); listIt++){
//    std::cout << listIt->ToString() << std::endl;
//    }

    if(updated > 0){
      /*
       * update was performed, evaluate list
       * -> this will set upperTriangleValue_[idElem]
       */
      Evaluate_GlobalSwitchingList(idElem);
      /*
       * lastXpar marks the last value which was used to update the switching list
       */
      lastXpar_[idElem] = xPar;
    }

    Yout += e_u * upperTriangleValue_[idElem];

    /*
     * Part 3 - upper square
     */
    Update_GlobalRotationList(X_thres, e_u, xVal, idElem);

    /*
     * Evaluate_GlobalRotationList checks for each element if it was changed or not and
     * reevaluates only the ones that did
     */
    Evaluate_GlobalRotationList(idElem, retVec);

    Yout += retVec;

//    std::cout << "###YOUT###################" << std::endl;
//    std::cout << std::setprecision(12) << std::scientific << Yout.ToString() << std::endl;
//    std::cout << "##########################" << std::endl;
    /*
     * in previous versions we scaled with delta_^2, too
     * this is no longer needed as we perform this step already in mapRectangleToPreisachWeights which is called
     * for the evaluation of all three parts
     */
    preisachSum_[idElem] = Yout*YSaturated_;

    return preisachSum_[idElem];
  }

  // calculate overlapping area of two rectangles; returns true if overlap exists
  bool VectorPreisachv7::clipRectangles(Double t1, Double b1, Double l1, Double r1, Double t2, Double b2, Double l2, Double r2, Double& tRet, Double& bRet, Double& lRet, Double& rRet){
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

  Double VectorPreisachv7::clipToElement(Double rectT, Double rectB, Double rectL, Double rectR, UInt idAlpha, UInt idBeta, Double delta){
    /*
     * Calculates the overlapping area of a rectangle (rectT,rectB,rectL,rectR) with the element defined by alphaId, betaId
     * This calculation does the following steps:
     *  1. calculate the overlapping rectangle via clipRectangles
     *  2. if success, check if lower left corner or upper right corner lie below diagonal
     *    -> if yes, cut down so that only the right corner is still below the diagonal
     *    -> subtract triangular area
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
    success = clipRectangles(rectT, rectB, rectL, rectR, elemBot+delta,elemBot, elemLeft, elemLeft+delta, ovT, ovB, ovL, ovR);

    if(!success){
      return 0.0;
    }

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

  void VectorPreisachv7::Evaluate_GlobalRotationList(UInt idElem, Vector<Double>& retVec){
    /*
     * Calculates the weighted rotation state of the upper square region
     * by evaluating the global rotation list (weighted rotation state =
     * rotation state * switching state * preisach weight)
     *
     * Evaluation is done in the following way:
     *
     * for rotEntry in globalRotList:
     *   calculate two rectangular regions corresponding to the rotation area (forming a L)
     *
     *   for entry in localSwitchList:
     *     get rectangular bounds from entry in switch list
     *     clip rectangle with first rotation area
     *     clip result to preisach plane
     *
     *     clip rectangle with second rotation area
     *     clip result to preisach plane
     *
     *     set entry of switching list to not changed
     *     sum up both parts with appropriate signs
     *
     *   set entry of rotation list to not changed
     *   multiply sum with rotation state and add product to eval state
     *
     * save eval state in retVec
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
     * iterators for local swtiching list
     */
    std::list<ListEntry>::iterator swListIt;
    std::list<ListEntry>::iterator swListEnd; // = --(globSwitchList_[idElem].end());

    /*
     * iterators for global rotation list
     */
    std::list<RotListEntry>::iterator rotListIt;
    std::list<RotListEntry>::iterator rotListEnd = --(globRotList_[idElem].end());

    /*
     * in former version we were not allowed to directly increase the rotListIt as we needed the first entry twice
     * due to new treatment, this is no longer needed; however it is still needed for the local switching lists
     */

    bool twoAreas = true;
    Double upperBound;
    Double lowerBound;

    Double sum = 0.0;
    Double sw_t,sw_b,sw_l,sw_r;
    Double rot_t1,rot_b1,rot_l1,rot_r1;
    Double rot_t2,rot_b2,rot_l2,rot_r2;
    Double ov_t,ov_b,ov_l,ov_r;
    Double tmp = 0.0;
    Double factor;
    UInt cnt = 0;

    Vector<Double> rotState;

    for(rotListIt = globRotList_[idElem].begin(); rotListIt != globRotList_[idElem].end(); rotListIt++){

      twoAreas = true;
      rotState = rotListIt->getVecReference();

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
        retVec += rotState*rotListIt->getLastEvalState();
      } else {

        /*
         * extract upper and lower bound defining the dimensions of the L-shape
         *
         *     one possible rotation state:
         *
         *                                   alpha = 1, beta = 0
         *         __ __ __ __ __ __ __ __ __ __v  _ upper bound 0
         *        |        |     |        |     |
         *        |  0     | 1   |  2     | 3   |  _ lower bound 0; upper bound 1
         *        |________|     |        |     |
         *        |              |        |     |
         *        |______________|        |     |  _ lower bound 1; upper bound 2
         *        |                       |     |
         *        |                       |     |
         *        |_______________________|     |  _ lower bound 2; upper bound 3
         *        |                             |
         *        |__ __ __ __ __ __ __ __ __ __|  _ lower bound 3
         *        ^
         *     alpha = 0, beta = -1
         *
         *
         *    area 0 -> no L-shape -> clip only against one rectangle -> twoAreas = false
         *                              left = -1, right = - lower bound i, top = 1, bot = lower bound
         *                              (note: the bounds are positive values; on beta you have to use - that value!)
         *    area i>0 -> L-shape -> decompose into two rectangles
         *                              1. left = -1, right = - lower bound i, top = upper bound i, bot = lower bound i
         *                              2. left = - upper bound i, right = - lower bound i, top = 1, bot = upper bound i
         *
         */

        upperBound = rotListIt->getVal();
        lowerBound = rotListIt->getLowerVal();

        if(upperBound >= 1){
          twoAreas = false;
        }

        /*
         * set coordinates of first and second area (even if the second one has zero size)
         * if upperBound == lowerBound the later functions will return 0 as update
         */
        rot_t1 = upperBound;
        rot_b1 = lowerBound;
        rot_l1 = -1.0;
        rot_r1 = -lowerBound;

        rot_t2 = 1.0;
        rot_b2 = upperBound;
        rot_l2 = -upperBound;
        rot_r2 = -lowerBound;

        swListEnd = --(rotListIt->getListReference().end());

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
         *  anymore (as this is the special helper area). So we startCnt in rotListEntry to mark that we
         *  deleted the first entry which originally corresponded to area 0. As all areas with id > 0 are
         *  calculated by the same scheme, the absolute value of cnt is not of interest and thus it is
         *  either 0 or 1.
         */
        cnt = rotListIt->getStartCnt();
        sum = 0.0;
        tmp = 0.0;
        bool success = false;

        /*
         * here we do not increase the iterator directly as we need to evaluate the first entry twice
         *
         * Issue:
         *  For input signals which gradually decrease and change their rotation state, the switching lists
         *  can be very very long. Due to a first optimization, we only insert values which might have an influence,
         *  i.e. entries of xPar which are > lowerVal (max) or < -lowerVal (min). Other entries would just
         *  span switching areas which cannot have an overlap with the rotation area.
         *  However, the switching list can still be very large as large lists get inherited from previous states
         *  and then get extended. These lists will have entries which might be completely unused as they
         *  encompass values which are above the upper val of the rotation area and thus cannot produce an overlap
         *  either.
         * Suggestion:
         *  Add function simplify_switchinglist. This function shall iterate over a switching list and
         *  eliminate all entries which produce no overlap. This reduces the chance of merging so that
         *  we should call this functions after we checked for a possible merge of rotation entries.
         *  The merging is still usable as its main purpose still can be trigger which is to merge
         *  adjacent rotation states after an increasing input completely overwrote the states (i.e
         *  the rotation state is the same and the switching list got overwritten).
         */
        for(swListIt = rotListIt->getListReference().begin(); swListIt != rotListIt->getListReference().end(); ){

          /*
           * get rectangular bounds corresponding to entry of switching list (compare to Evaluate_GlobalSwitchingList)
           * (do not forget to set flag isGlobal to false!)
           */
          factor = getRectangleBounds(rotListIt->getListReference(),swListIt,swListEnd,cnt,false,sw_t,sw_b,sw_l,sw_r);

          /*
           * clip resulting area against rotation area1
           */
          success = clipRectangles(rot_t1, rot_b1, rot_l1, rot_r1, sw_t, sw_b, sw_l, sw_r, ov_t, ov_b, ov_l, ov_r);

          if(success){
            /*
             *  clip overlap to PreisachPlane and sum up
             */
            tmp = mapRectangleToPreisachWeights(ov_t,ov_b,ov_l,ov_r);

            sum += factor*tmp;
          }

          if(twoAreas){
            /*
             * repeat the clipping steps for the second area
             */
            success = clipRectangles(rot_t2, rot_b2, rot_l2, rot_r2, sw_t, sw_b, sw_l, sw_r, ov_t, ov_b, ov_l, ov_r);

            if(success){
             /*
              *  clip overlap to PreisachPlane and sum up
              */
             tmp = mapRectangleToPreisachWeights(ov_t,ov_b,ov_l,ov_r);

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
              //std::cout << "Doing extra work" << std::endl;
              /*
               * set flag upperSplitSquare to true -> get area 0 instead of area 1
               */
              factor = getRectangleBounds(rotListIt->getListReference(),swListIt,swListEnd,cnt,false,sw_t,sw_b,sw_l,sw_r,true);

              if(factor != 0){
                EXCEPTION("Something got wrong here! Check function getRectangleBounds!");
              }

              /*
               * clip resulting area against rotation area1
               */
              success = clipRectangles(rot_t1, rot_b1, rot_l1, rot_r1, sw_t, sw_b, sw_l, sw_r, ov_t, ov_b, ov_l, ov_r);

              if(success){
                /*
                 *  clip overlap to PreisachPlane and sum up
                 *  -> NOTE: here we set the flag skipUpperDiagonal to true!
                 */
                tmp = mapRectangleToPreisachWeights(ov_t,ov_b,ov_l,ov_r,true);
                // factor is not needed here; mapRectangleToPreisachWeights with skipUpperDiagonal == true
                // will not only skip the diagonal entries, but also consider the right signs!
                sum += tmp;
              }

              if(twoAreas){
                /*
                 * repeat the clipping steps for the second area
                 */
                success = clipRectangles(rot_t2, rot_b2, rot_l2, rot_r2, sw_t, sw_b, sw_l, sw_r, ov_t, ov_b, ov_l, ov_r);

                if(success){
                  /*
                   *  clip overlap to PreisachPlane and sum up
                   *  -> NOTE: here we set the flag skipUpperDiagonal to true!
                   */
                  tmp = mapRectangleToPreisachWeights(ov_t,ov_b,ov_l,ov_r,true);
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
        retVec += rotState*sum;
      }
    } // rot list
  }

  void VectorPreisachv7::Evaluate_GlobalSwitchingList(UInt idElem){
    /*
     * This function calculates the overall switching state of the upper triangle;
     * This value has to be multiplied with current rotation state e_u in each iteration to get
     * the overall contribution of the upper triangle
     * -> this function has to be called in each step, if the global switching list changed
     * -> after a successful evaluation, set the value upperTriangularVal_ and set the flag
     * hasChanged_ in each entry to false
     */

    Double sum = 0.0;
    Double t,b,l,r;
    Double tmp = 0.0;
    Double factor = 0.0;
    UInt cnt = 0;

// NOT NEEDED as we call this function only after checking if update was performed or not
//      if(GlobalSwitchingChanged(idElem) == false){
//        return;
//      }

    /*
     * iterate over list of minima and maxima
     * -> from these values we can reproduce the stair case function describing the switching state
     * -> cut stair-case into rectangles and overlap these rectangles with the Preisach plane
     * -> depeding on the entry to be a min or a max, the value will be added or subtracted
     *    (the required sign will now be returend by getRectangleBounds which is called to determine
     *      the rectangular parts)
     */
    std::list<ListEntry>::iterator listIt;
    std::list<ListEntry>::iterator listEnd = --(globSwitchList_[idElem].end());

    for(listIt = globSwitchList_[idElem].begin(); listIt != globSwitchList_[idElem].end();){

      factor = getRectangleBounds(globSwitchList_[idElem],listIt,listEnd,cnt,true,t,b,l,r);

      /*
       * calculate overlapping of current switching area with Preisach plane
       * (returns 0.0 if no overlap exists)
       */
      tmp = mapRectangleToPreisachWeights(t,b,l,r);

      /*
       * In older version, where we did not clip directly to the Preisach plane, we subtracted
       * areas below the diagonal -> this is done already in mapRectangleToPreisachWeights
       */
      sum += factor*tmp;

      if(cnt > 0){
        /*
         * NOTE: area1 and area2 are both calculated using the first list entry, therefore we do not increase the iterator after
         * the first iteration
         */
         listIt++;
      }
      cnt++;

      /*
       * set flag isChanged of the current entry to false
       * (so that we later can check if list changed or not)
       */
      listIt->setHasChanged(false);

    }

    upperTriangleValue_[idElem] = sum;
  }

  void VectorPreisachv7::Evaluate_LowerTriangle(){
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
    lowerTriangleValue_ = mapRectangleToPreisachWeights(0.0,-1.0,-1.0,0.0);
  }

  Double VectorPreisachv7::mapRectangleToPreisachWeights(Double t, Double b, Double l, Double r, bool skipUpperDiagonal){
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
    UInt rowMin = (UInt) floor(b/delta_) + floor(1.0/delta_);
    UInt rowMax = (UInt) ceil(t/delta_) + floor(1.0/delta_)-1;

    UInt colMin = (UInt) floor(l/delta_) + floor(1.0/delta_);
    UInt colMax = (UInt) ceil(r/delta_) + floor(1.0/delta_)-1;

    /*
     * 2. get indices of completely overlapped elements
     */
    UInt rowMinFull = (UInt) ceil(b/delta_) + floor(1.0/delta_);
    UInt rowMaxFull = (UInt) floor(t/delta_) + floor(1.0/delta_)-1;

    UInt colMinFull = (UInt) ceil(l/delta_) + floor(1.0/delta_);
    UInt colMaxFull = (UInt) floor(r/delta_) + floor(1.0/delta_)-1;

//    std::cout << std::setprecision(12) << std::scientific << b << " " << t << " " << l << " " << r << std::endl;
//    std::cout << std::setprecision(12) << std::scientific << floor(b) << " " << ceil(t) << " " << floor(l) << " " << ceil(r) << std::endl;
//    std::cout << std::setprecision(12) << std::scientific << floor(b/delta_) << " " << ceil(t/delta_) << " " << floor(l) << " " << ceil(r) << std::endl;
 //   std::cout << rowMin << " " << rowMax << " " << colMin << " " << colMax << std::endl;
 //   std::cout << rowMinFull << " " << rowMaxFull << " " << colMinFull << " " << colMaxFull << std::endl;

    /*
     * Iterate over all elements
     */
    Double tmp = 0.0;

    if(skipUpperDiagonal == true){
      /*
       * special treatment:
       * 1. skip entries on alpha = -beta
       * 2. consider the right signs, i.e. +1 for entries below alpha = -beta
       *    and -1 for entries above!
       */
      Double sign;
      for(UInt i = rowMin; i <= rowMax; i++){
        /*
         * ensure, that j <= i -> other elements do not contribute
         */
        for(UInt j = colMin; j <= std::min(i,colMax); j++){

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
          if((i >= rowMinFull)&&(i <= rowMaxFull)){
            if((j >= colMinFull)&&(j <= colMaxFull)){
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
              tmp = tmp * clipToElement(t, b, l, r, i, j);
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
            tmp = tmp * clipToElement(t, b, l, r, i, j);
          }

          sum += sign*tmp;
        }
      }

    } else {
      /*
       * std treatment
       */
      for(UInt i = rowMin; i <= rowMax; i++){
        /*
         * ensure, that j <= i -> other elements do not contribute
         */
        for(UInt j = colMin; j <= std::min(i,colMax); j++){

          tmp = preisachWeights_[i][j];

          /*
           * Check for an inner element
           */
          if((i >= rowMinFull)&&(i <= rowMaxFull)){
            if((j >= colMinFull)&&(j <= colMaxFull)){
              /*
               * check for diagonal element on alpha = beta
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
              tmp = tmp * clipToElement(t, b, l, r, i, j);
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
            tmp = tmp * clipToElement(t, b, l, r, i, j);
          }

          sum += tmp;
        }
      }
    }

    return sum;
  }


  Double VectorPreisachv7::getRectangleBounds(std::list<ListEntry>& list,std::list<ListEntry>::iterator startIt, std::list<ListEntry>::iterator endIt,
                            UInt idArea, bool isGlobal, Double& tRet, Double& bRet, Double& lRet, Double& rRet, bool upperSplitSquare){
    /*
     * This functions is used to decompose a stair-case switching state into rectangular areas which are
     * easier to compute.
     * This functions will only be used for the switching lists (local and global)
     *
     * upperSplitSquare:
     *  If true, the coordinates of the diagonally split upper square (area 0) will be returned)
     *  (cnt will be set to 0, to trigger the right calculations)
     *
     */
    /*
     * Reuse large parts of the original getRectangularBounds function (which was used to determine the rectangles
     * inside each element). The algorithm can be reused in that sense, that we can consider the upper square (isGlobal = false)
     * or the upper triangle (isGlobal = true) to be two large elements.
     */

    if(upperSplitSquare == true){
      idArea = 0;
    }

    Double alpha, beta, delta;
    if(isGlobal == true){
      alpha = 0.0;
      beta = 0.0;
    } else {
      alpha = 0.0;
      beta = -1.0;
    }

    /*
     * use delta=1.0 in the following (compared to the original delta_)
     */
    delta = 1.0;

    /*
    * Calculates the coordinates of a rectangle lying inside element alpha,beta
    *
    * Input:
    *  list -> switching list (either global or a local one from a rotListEntry)
    *  startIt -> iterator pointing to the CURRENT entry
    *  endIt -> iterator pointing to the last entry
    *  alpha,beta -> bottom left corner coordinates of element
    *  idArea -> index of area to be computed - 1 (i.e. area1 -> idArea = 0)
    *  isGlobal -> flag checking if element is on diagonal alpha = beta (true)  or on the diagonal alpha = -beta (false)
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

    /*
    * Notes to the switching list and the decompositions of the stair-case shape
    */
    /*
    * SQUARE (i.e. isGlobal == false):
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
    *
    *       area 0
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
    *        -> solution should be appropriate now!
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
    *      NOTE:
    *       The split area can only be neglected, if the underlying weights are symmetric to alpha = -beta!
    *       Otherwise, both halves do not cancel out!
    */
    /*
    * TRIANGLE (i.e. isGlobal == true):
    *
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
    *
    *  -> not needed in new version anymore, as the element is either split along the diagonal alpha = -beta
    *  or along alpha = beta
    *
    */
    /*
    * NEW: There is no reason, why the sign could not be determined by this function directly
    * -> add according steps to end
    */

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

     top = alpha+delta - entry + beta;

     if(startIt != endIt){
       startIt++;
       if(startIt->isMin() == isMinEntry){
         EXCEPTION("MinMaxList has to be alternating!")
       } else {
         bot = startIt->getVal();
       }
       startIt--;
     } else {
       bot = alpha;
     }

     /*
      * min-helper area is positive, so +
      */
     returnVal = 1.0;

     /*
      * Special treatment for area 0
      */
     if(upperSplitSquare == true){
       /*
        * do not take area 1 but instead take the neighboring area 0 (above of area 1)
        * -> top side of area 1 = bot side of area 0
        * -> top of area 0 = alpha + delta
        * -> left and right stay the same
        */
       bot = top;
       top = alpha + delta;
       returnVal = 0.0; // not needed, but by this we can check if value was hit
     }

     if(top < bot){
       /*
        * NEW 10.5.2016
        *
        * check if top is smaller than bot!
        * if yes, we have the case, that the helper area overlaps with the next area (see comment above)
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
        * the resulting area must have an opposite sign compared to the normal sign as it is not added to the positive part
        * but to the negative one!
        */
       returnVal = -1.0;
       left = -bot; // bot already holds the value of nextentry
       right = entry;

       top = (alpha+delta);
       //bot = bot; // //-1 for creating a negative area!

       /*
        * Special treatment for area 0
        */
       if(upperSplitSquare == true){
         /*
          * do not take area 1* but instead take the neighboring area 0 (left of area 1*)
          * -> left side of area 1* = right side of area 0
          * -> left side of area 0 = beta
          * -> top and bot stay the same
          */
         right = left;
         left = beta;
         returnVal = 0.0; // not needed, but by this we can check if value was hit
       }
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
         EXCEPTION("MinMaxList has to be alternating!")
       } else {
         right = startIt->getVal();
       }
       startIt--;
     } else {
       right = beta+delta;
     }

     left = alpha+delta - entry + beta;

     top = alpha+delta;
     bot = entry;

     /*
      * max-helper area counts to the negative part, so -1
      */
     returnVal = -1.0;

     /*
      * Special treatment for area 0
      */
     if(upperSplitSquare == true){

       /*
        * do not take area 1 but instead take the neighboring area 0 (left besides area 1)
        * -> left of area 1 = right of area 0
        * -> left of area 0 = beta
        * -> top and bot stay the same
        */
       right = left;
       left = beta;
       returnVal = 0.0; // not needed, but by this we can check if value was hit
     }

     if(right < left){
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
        *
        */
       returnVal = 1.0;
       top = -right; // right already holds the value of nextentry
       bot = entry;

       left = beta; //-1 for creating a negative area!
       //right = right; //-1 for creating a negative area!

       /*
        * Special treatment for area 0
        */
       if(upperSplitSquare == true){

         /*
          * do not take area 1* but instead take the neighboring area 0 (above area 1*)
          * -> top of area 1* = bot of area 0
          * -> top of area 0 = alpha + delta
          * -> left and right stay the same
          */
         bot = top;
         top = alpha + delta;
         returnVal = 0.0; // not needed, but by this we can check if value was hit
       }
     }

    } else if(isMinEntry && (idArea > 0)) {
     /*
      * Current entry is negative
      */
     /*
      * std-min area counts to the negative part, so -1
      */
     returnVal = -1.0;

     right = beta+delta;

     if(startIt != endIt){
       startIt++;

       if(startIt->isMin() == isMinEntry){
         EXCEPTION("MinMaxList has to be alternating!")
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
           right = startIt->getVal();
         }
         startIt--;
       }
       startIt--;
     } else {

       bot = alpha;
     }

     top = alpha+delta;
     left = entry;

    } else {
     /*
      * Current entry is positive
      */
     /*
      * std-max area counts to the negative part, so -1
      */
     returnVal = 1.0;

     bot = alpha;

     if(startIt != endIt){
       startIt++;
       if(startIt->isMin() == isMinEntry){
         EXCEPTION("MinMaxList has to be alternating!")
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
           bot = startIt->getVal();
         }
         startIt--;
       }
       startIt--;
     } else {
       right =  beta+delta;
     }

     top = entry;
     left = beta;
    }

    /*
    * in previous version we checked flag isTriangle; in the new version, isTriangle would only be true for the
    * upper triangular entry which is set via isGlobal
    *
    * REMARK:
    * This cutting step might be unecessary due to the following reason:
    *   In the regions, where this would be relevant (upper triangle), the resulting rect will be mapped
    *   onto the Preisach plane directly (i.e. without further overlapping); the mapping function then
    *   divides the area into parts which overlap only a single element each; thereby, all area parts
    *   which lie below the diagonal will be omitted, so that it makes no difference if they have a
    *   triangular shape or a trapezoidal shape
    */
    //#######################
    if(isGlobal == true){
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
     * Do not cut away the areas below the diagonal yet!
     * -> it is better to clip the returning rectangle to other rectangles first and later cut away
     *    the area below diagonal
     * -> IMPORTANT: the second area to be clipped against should already be cut down like this area, i.e.
     *    only the bot-right corner should lie below the diagonal
     */
    //######################

    tRet = top;
    bRet = bot;
    lRet = left;
    rRet = right;

    return returnVal;
  }

  void VectorPreisachv7::Simplify_LocalSwitchingLists(UInt idElem){
    /*
     * This function iterates over the globalRotation list and checks each entry in the corresponding
     * local switching list for overlap with the two possible rotation areas.
     * This search is performed until at least one overlap is found. All switching entries BEFORE the
     * one which lead to the first overlap are removed from the list. If at least 1 entry got removed,
     * startCnt will be set to 1 (as the entry of the switching list which corresponds to the helper area
     * (id = 0) no longer is in the list).
     * Leave out all rotListEntries for which the flag isUpdated is false.
     */

    std::list<ListEntry>::iterator swListIt;
    std::list<ListEntry>::iterator firstToKeep;
    std::list<ListEntry>::iterator swListEnd; // = --(globSwitchList_[idElem].end());

    /*
     * iterators for global rotation list
     */
    std::list<RotListEntry>::iterator rotListIt;
    std::list<RotListEntry>::iterator rotListEnd = --(globRotList_[idElem].end());

    bool twoAreas;
    Double upperBound, lowerBound;
    Double sw_t,sw_b,sw_l,sw_r;
    Double rot_t1,rot_b1,rot_l1,rot_r1;
    Double rot_t2,rot_b2,rot_l2,rot_r2;
    Double ov_t,ov_b,ov_l,ov_r;

    UInt cnt;

    for(rotListIt = globRotList_[idElem].begin(); rotListIt != globRotList_[idElem].end(); rotListIt++){

      twoAreas = true;
      /*
       * maybe the list was already simplified once or was inherited from a previous rotListState
       * in that case we have to make sure to not test for the wrong value
       */
      cnt = rotListIt->getStartCnt();

      /*
       * check if rotListEntry changed at all
       * if not, the switching list is the same as in the last time and thus
       * can be considered as already simplified
       */
      if(rotListIt->hasChanged() == false){
        continue;
      } else {
        /*
         * Similar procedure as in evalulate_globalrotationlist
         */

        upperBound = rotListIt->getVal();
        lowerBound = rotListIt->getLowerVal();

        if(upperBound >= 1){
          twoAreas = false;
        }

        /*
         * set coordinates of first and second area (even if the second one has zero size)
         * if upperBound == lowerBound the later functions will return 0 as update
         */
        rot_t1 = upperBound;
        rot_b1 = lowerBound;
        rot_l1 = -1.0;
        rot_r1 = -lowerBound;

        rot_t2 = 1.0;
        rot_b2 = upperBound;
        rot_l2 = -upperBound;
        rot_r2 = -lowerBound;

        firstToKeep = rotListIt->getListReference().begin();
        swListEnd = --(rotListIt->getListReference().end());

        bool success1 = false;
        bool success2 = false;
        bool tmpSuccess = false;

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
           * (do not forget to set flag isGlobal to false!)
           */

          getRectangleBounds(rotListIt->getListReference(),swListIt,swListEnd,cnt,false,sw_t,sw_b,sw_l,sw_r);

          /*
           * clip resulting area against rotation area1
           */
          success1 = clipRectangles(rot_t1, rot_b1, rot_l1, rot_r1, sw_t, sw_b, sw_l, sw_r, ov_t, ov_b, ov_l, ov_r);


          if(twoAreas){
            /*
             * repeat the clipping steps for the second area
             */
            success2 = clipRectangles(rot_t2, rot_b2, rot_l2, rot_r2, sw_t, sw_b, sw_l, sw_r, ov_t, ov_b, ov_l, ov_r);
          }

          if((success1 == true)||(success2==true)){
            /*
             * the current switching area has an overlap with at least one of the two rotation areas
             * -> this entry has to be kept
             * -> end loop (as we suppose that all later entries have an overlap, too
             * This assumption is ok as only those entries are appended to the list which actually have
             * and effect, i.e. we only have to cut at the start of the list!
             */
            break;
          }

          if(cnt > 0){

            /*
             * entry has no overlap, check next one
             * (start increasing only if cnt > 0 as list entry corrsponding to cnt == 0
             * has to be checked for area 0 and 1.
             */
            swListIt++;
            firstToKeep++;
            if((tmpSuccess == true)||(success1 == true)||(success2 == true)){
              /*
               * we either have found an overlap for cnt = 0 (tmpSuccess = true) or for
               * the current
               */
            }
          }
          cnt++;
        } // sw list

        if(firstToKeep != rotListIt->getListReference().begin()){
          /*
           * At least one entry has to be removed from the list
           */
          rotListIt->setStartCnt(1);
          rotListIt->getListReference().erase(rotListIt->getListReference().begin(),firstToKeep);
        } //else
        /*
         * whole list has to be kept
         * -> continue with next rotListElement
         */
      } // rot list has changed
    } // rot list
  }


  void VectorPreisachv7::Simplify_GlobalRotationList(UInt idElem){
    /*
     * iterate over rotation list and check if two adjacent entries have the same rotation state and the
     * same list of switching entries
     * This case can happen e.g. if an input overrides multiple older rotation areas at once and the
     * resulting value of xPar is large enough to overwrite the contained switching lists
     */

    if(globRotList_[idElem].size() < 2){
      /*
       * list has not enough entries to be merged together
       */
      return;
    }

    std::list<RotListEntry>::iterator listIt;
    std::list<RotListEntry>::iterator nextListIt;
    std::list<RotListEntry>::iterator listEnd = --(globRotList_[idElem].end());

    for(listIt = globRotList_[idElem].begin(); listIt != globRotList_[idElem].end(); ){

      /*
       *  get next entry in list (if any)
       *  if not -> nothing more to do here
       */
      if(listIt != listEnd){
        nextListIt = listIt;
        nextListIt++;

        /*
         * check if entries can be merged together
         */
        if(listIt->canBeMergedWith(*nextListIt)){
          /*
           * get lower bound of second element and set lower bound of first element to that value
           * (-> i.e. we extend the area of the first element)
           */
          Double low = nextListIt->getLowerVal();
          listIt->setLowerVal(low);
          /*
           * now we can delete the second entry
           */
          globRotList_[idElem].erase(nextListIt);

          /*
           * reset iterator to last entry
           * (normally the iterators should remain valid, but I do not trust them ...
           */
          listEnd = --(globRotList_[idElem].end());

          /*
           * do not increase iterator
           * -> it might be, that the by now extended first entry matches also the next one
           */
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
  }

  void VectorPreisachv7::mapRectangleToHelperMatrix(Matrix<Double>& helper, Double t, Double b, Double l, Double r, Double factor, bool skipUpperDiagonal){
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
     */

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
    UInt rowMin = (UInt) floor(b/delta) + floor(1.0/delta);
    UInt rowMax = (UInt) ceil(t/delta) + floor(1.0/delta)-1;

    UInt colMin = (UInt) floor(l/delta) + floor(1.0/delta);
    UInt colMax = (UInt) ceil(r/delta) + floor(1.0/delta)-1;

    /*
     * 2. get indices of completely overlapped elements
     */
    UInt rowMinFull = (UInt) ceil(b/delta) + floor(1.0/delta);
    UInt rowMaxFull = (UInt) floor(t/delta) + floor(1.0/delta)-1;

    UInt colMinFull = (UInt) ceil(l/delta) + floor(1.0/delta);
    UInt colMaxFull = (UInt) floor(r/delta) + floor(1.0/delta)-1;

    if(skipUpperDiagonal == true){
      /*
       * special treatment:
       * 1. skip entries on alpha = -beta
       * 2. consider the right signs, i.e. +1 for entries below alpha = -beta
       *    and -1 for entries above!
       */
      Double sign;
      for(UInt i = rowMin; i <= rowMax; i++){
        /*
         * ensure, that j <= i -> other elements do not contribute
         */
        for(UInt j = colMin; j <= std::min(i,colMax); j++){

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
          if((i >= rowMinFull)&&(i <= rowMaxFull)){
            if((j >= colMinFull)&&(j <= colMaxFull)){
              /*
               * check for diagonal element on alpha = beta
               * (should not happen here)
               */
              if(i == j){
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
              helper[i][j] += sign*clipToElement(t, b, l, r, i, j, delta)/(delta*delta);
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
            helper[i][j] += sign*clipToElement(t, b, l, r, i, j, delta)/(delta*delta);
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
      for(UInt i = rowMin; i <= rowMax; i++){
        /*
         * ensure, that j <= i -> other elements do not contribute
         */
        for(UInt j = colMin; j <= std::min(i,colMax); j++){
          /*
           * Check for an inner element
           */
          if((i >= rowMinFull)&&(i <= rowMaxFull)){
            if((j >= colMinFull)&&(j <= colMaxFull)){
              /*
               * check for diagonal element
               */
              if(i == j){
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
              helper[i][j] += factor*clipToElement(t, b, l, r, i, j, delta)/(delta*delta);
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
            helper[i][j] += factor*clipToElement(t, b, l, r, i, j, delta)/(delta*delta);
          }
        }
      }
    }
  }

  void VectorPreisachv7::switchingStateToBmp(UInt numPixel, std::string filename, UInt idElem, bool overLayWithRotState)
  {
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
     * Fill matrix
     * 1. lower triangle part -> fill with +1 (0.5 on diagonal)
     * 2. upper triangle part -> evaluate global switching list
     * 3. square part -> evaluate global rotation list
     */
    /*
     * Part 1
     */
    for(UInt i = 0; i < numPixel/2; i++){
      for(UInt j = 0; j <= i; j++){
        if(i == j) helperMatrix[i][j] = 0.5;
        else helperMatrix[i][j] = 1.0;
      }
    }
    /*
     * Part 2
     */

    bool twoAreas = true;
    Double upperBound;
    Double lowerBound;

    Double sw_t,sw_b,sw_l,sw_r;
    Double rot_t1,rot_b1,rot_l1,rot_r1;
    Double rot_t2,rot_b2,rot_l2,rot_r2;
    Double ov_t,ov_b,ov_l,ov_r;
    Double factor;
    UInt cnt = 0;
    Double t,b,l,r;
    std::list<ListEntry>::iterator listIt;
    std::list<ListEntry>::iterator listEnd = --(globSwitchList_[idElem].end());

    for(listIt = globSwitchList_[idElem].begin(); listIt != globSwitchList_[idElem].end();){

      factor = getRectangleBounds(globSwitchList_[idElem],listIt,listEnd,cnt,true,t,b,l,r);

      mapRectangleToHelperMatrix(helperMatrix,t,b,l,r,factor);

      if(cnt > 0){
        /*
         * NOTE: area1 and area2 are both calculated using the first list entry, therefore we do not increase the iterator after
         * the first iteration
         */
         listIt++;
      }
      cnt++;
    }

    /*
     * Part 3
     */

    /*
     * iterators for local swtiching list
     */
    std::list<ListEntry>::iterator swListIt;
    std::list<ListEntry>::iterator swListEnd; // = --(globSwitchList_[idElem].end());

    /*
     * iterators for global rotation list
     */
    std::list<RotListEntry>::iterator rotListIt;
    std::list<RotListEntry>::iterator rotListEnd = --(globRotList_[idElem].end());

    bool area0 = true;
    /*
     * NEW Part 3* -> upper square part which has no rotation state yet (-> does not contribute to result,
     * however, it is diagonally split by into +1 and -1 from start and we want to visualize that, too)
     */
    /*
     * Rotation areas
     *
     *                     alpha = 1, beta = 0
     *         __ __ __ __ __ __ __ __ __ __v  _ upper bound 0
     *        |        |     |        |     |
     *        |  0     | 1   |  2     | 3   |  _ lower bound 0; upper bound 1
     *        |________|     |        |     |
     *        |              |        |     |
     *        |______________|        |     |  _ lower bound 1; upper bound 2
     *        |                       |     |
     *        |                       |     |
     *        |_______________________|     |  _ lower bound 2; upper bound 3
     *        |                             |
     *        |__ __ __ __ __ __ __ __ __ __|  _ lower bound 3
     *
     */
    /*
     * overlap area0
     * 1. with switching list (it is possible that a switching area for the first rotation state 1
     *      reaches out into rotation area 0)
     * 2. with not initialized switching state
     *
     * i.e. we check for overlap of rotation area 0 with switching areas 0,1,2,...
     *
     * Switching areas (for first min/max = max)
     *      area 0
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
     *
     * SPECIAL NOTE:
     *  Function simplify_localswitchinglists is designed such, that unneeded switching states are
     *  deleted from the list. Unneeded means, that these areas do not overlap with the rotation areas.
     *  This is fine for the computation as the upper area (even in the unsymmetric case) will not contribute
     *  if the rotation state of rotation area 0 is not yet set.
     *  However, for the output, the information of the actual switching state of rotation area 0 is lost
     *  after the simplification.
     *  Example: assume that rotation area 0 would overlap with switching area 0,1,2 and partially with 3 and 4.
     *  The function simplify_localswitchinglists will delete the minima and maxima which are needed to
     *  define area 0,1,2 because these areas do not overlap with the rotation area 1. The result of the
     *  computation will not change by this, as the actual switching state would be multiplied with a rotation
     *  state of 0. The output, however, will see an unitialized upper part consisting of switching area 0,1,2
     *  although only area 0 is diagonally split.
     *
     */

    /*
     * in former version we were not allowed to directly increase the rotListIt as we needed the first entry twice
     * due to new treatment, this is no longer needed; however it is still needed for the local switching lists
     */
    for(rotListIt = globRotList_[idElem].begin(); rotListIt != globRotList_[idElem].end(); rotListIt++){

      twoAreas = true;

      /*
       * extract upper and lower bound defining the dimensions of the L-shape
       *
       *     one possible rotation state:
       *
       *                                   alpha = 1, beta = 0
       *         __ __ __ __ __ __ __ __ __ __v  _ upper bound 0
       *        |        |     |        |     |
       *        |  0     | 1   |  2     | 3   |  _ lower bound 0; upper bound 1
       *        |________|     |        |     |
       *        |              |        |     |
       *        |______________|        |     |  _ lower bound 1; upper bound 2
       *        |                       |     |
       *        |                       |     |
       *        |_______________________|     |  _ lower bound 2; upper bound 3
       *        |                             |
       *        |__ __ __ __ __ __ __ __ __ __|  _ lower bound 3
       *        ^
       *     alpha = 0, beta = -1
       *
       *
       *    area 0 -> no L-shape -> clip only against one rectangle -> twoAreas = false
       *                              left = -1, right = - lower bound i, top = 1, bot = lower bound
       *                              (note: the bounds are positive values; on beta you have to use - that value!)
       *    area i>0 -> L-shape -> decompose into two rectangles
       *                              1. left = -1, right = - lower bound i, top = upper bound i, bot = lower bound i
       *                              2. left = - upper bound i, right = - lower bound i, top = 1, bot = upper bound i
       *
       */
      upperBound = rotListIt->getVal();
      lowerBound = rotListIt->getLowerVal();

      if(upperBound >= 1){
        twoAreas = false;
      }

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

        /*
         * area 0
         */
        rot_t1 = 1.0;
        rot_b1 = upperBound;
        rot_l1 = -1.0;
        rot_r1 = -upperBound;

        /*
         * NOTE: 16.06.2016
         * Formerly, we simply repeated all overlappings steps for an additional area 0.
         * This is wrong, as area0 takes the switching state of area1.
         * As area0 has no rotation state, all entries of xPar are 0 and thus we would only
         * have an empty list.
         * -> apply function mapRectangleToHelperMatrix to whole area0 and set flag
         *     skipUpperDiagonal to true
         * -> no overlaps with any switching states!
         */
        mapRectangleToHelperMatrix(helperMatrix,rot_t1,rot_b1,rot_l1,rot_r1,0,true);

        area0 = false;
      }

      /*
       * NEW 16.06.2016
       * Continue as usual with area 1 and 2
       */

      /*
       * set coordinates of first and second area (even if the second one has zero size)
       * if upperBound == lowerBound the later functions will return 0 as update
       */
      rot_t1 = upperBound;
      rot_b1 = lowerBound;
      rot_l1 = -1.0;
      rot_r1 = -lowerBound;

      rot_t2 = 1.0;
      rot_b2 = upperBound;
      rot_l2 = -upperBound;
      rot_r2 = -lowerBound;

      swListEnd = --(rotListIt->getListReference().end());

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
       *  anymore (as this is the special helper area). So we startCnt in rotListEntry to mark that we
       *  deleted the first entry which originally corresponded to area 0. As all areas with id > 0 are
       *  calculated by the same scheme, the absolute value of cnt is not of interest and thus it is
       *  either 0 or 1.
       */
      cnt = rotListIt->getStartCnt();
      bool success = false;

      /*
       * here we do not increase the iterator directly as we need to evaluate the first entry twice
       *
       * Issue:
       *  For input signals which gradually decrease and change their rotation state, the switching lists
       *  can be very very long. Due to a first optimization, we only insert values which might have an influence,
       *  i.e. entries of xPar which are > lowerVal (max) or < -lowerVal (min). Other entries would just
       *  span switching areas which cannot have an overlap with the rotation area.
       *  However, the switching list can still be very large as large lists get inherited from previous states
       *  and then get extended. These lists will have entries which might be completely unused as they
       *  encompass values which are above the upper val of the rotation area and thus cannot produce an overlap
       *  either.
       * Suggestion:
       *  Add function simplify_switchinglist. This function shall iterate over a switching list and
       *  eliminate all entries which produce no overlap. This reduces the chance of merging so that
       *  we should call this functions after we checked for a possible merge of rotation entries.
       *  The merging is still usable as its main purpose still can be trigger which is to merge
       *  adjacent rotation states after an increasing input completely overwrote the states (i.e
       *  the rotation state is the same and the switching list got overwritten).
       */

      for(swListIt = rotListIt->getListReference().begin(); swListIt != rotListIt->getListReference().end(); ){

        factor = getRectangleBounds(rotListIt->getListReference(),swListIt,swListEnd,cnt,false,sw_t,sw_b,sw_l,sw_r);

        /*
         * clip resulting area against rotation area1
         */
        success = clipRectangles(rot_t1, rot_b1, rot_l1, rot_r1, sw_t, sw_b, sw_l, sw_r, ov_t, ov_b, ov_l, ov_r);

        if(success){
          /*
           *  clip overlap to helper matrix
           */
          mapRectangleToHelperMatrix(helperMatrix,ov_t,ov_b,ov_l,ov_r,factor);
        }

        if(twoAreas){
          /*
           * repeat the clipping steps for the second area
           */
          success = clipRectangles(rot_t2, rot_b2, rot_l2, rot_r2, sw_t, sw_b, sw_l, sw_r, ov_t, ov_b, ov_l, ov_r);

          if(success){
           /*
            *  clip overlap to helper matrix
            */
            mapRectangleToHelperMatrix(helperMatrix,ov_t,ov_b,ov_l,ov_r,factor);
          }
        }

        /*
         * extra treatment for diagonally split area (only needed if
         * weights are not symmetric to alpha = -beta or for output visualization)
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
           * here we save the state even if the weights are symmetric
           */
          if((rotListIt->wasListWipedOut() == false)){
            /*
             * set flag upperSplitSquare to true -> get area 0 instead of area 1
             */
            factor = getRectangleBounds(rotListIt->getListReference(),swListIt,swListEnd,cnt,false,sw_t,sw_b,sw_l,sw_r,true);

            if(factor != 0){
              EXCEPTION("Something got wrong here! Check function getRectangleBounds!");
            }

            /*
             * clip resulting area against rotation area1
             */
            success = clipRectangles(rot_t1, rot_b1, rot_l1, rot_r1, sw_t, sw_b, sw_l, sw_r, ov_t, ov_b, ov_l, ov_r);

            if(success){
              /*
               *  clip overlap to PreisachPlane and sum up
               *  -> NOTE: here we set the flag skipUpperDiagonal to true!
               */
              mapRectangleToHelperMatrix(helperMatrix,ov_t,ov_b,ov_l,ov_r,factor,true);
            }

            if(twoAreas){
              /*
               * repeat the clipping steps for the second area
               */
              success = clipRectangles(rot_t2, rot_b2, rot_l2, rot_r2, sw_t, sw_b, sw_l, sw_r, ov_t, ov_b, ov_l, ov_r);

              if(success){
                /*
                 *  clip overlap to PreisachPlane and sum up
                 *  -> NOTE: here we set the flag skipUpperDiagonal to true!
                 */
                mapRectangleToHelperMatrix(helperMatrix,ov_t,ov_b,ov_l,ov_r,factor,true);
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

    } // rot list

    /*
     * now call output function of matrix
     */
    UInt upscaling = 1;

    /*
     * get information about the rotation state and overlap this information with the switching state
     * -> rotation state defines the green channel in the following way:
     *    rotation value > 0 -> green = 255
     *    rotation value <= 0 -> green = 0
     */

    if(overLayWithRotState == true){
      Matrix<Double>* rotStorage = new Matrix<Double>(numPixel,numPixel);

      for(UInt comp = 0; comp < dim_; comp++){

        rotationStateToBmp(numPixel, filename, idElem, comp, false, rotStorage);

        std::stringstream stream;
        if(comp == 0){
          stream << "x-" << filename;
        } else if(comp == 1){
          stream << "y-" << filename;
        } else {
          stream << "z-" << filename;
        }

        std::string filename_new = stream.str();

        helperMatrix.matrix2Bmp(upscaling,filename_new,rotStorage);
      }

      delete rotStorage;
    } else {
      helperMatrix.matrix2Bmp(upscaling,filename,NULL);
    }
  }

  void VectorPreisachv7::rotationStateToBmp(UInt numPixel, std::string filename, UInt idElem, UInt comp, bool writeOutput, Matrix<Double>* helperMatrix)
  {
    if(numPixel < 2){
      WARN("Image should have more than 2 x 2 pixel");
      return;
    }

    if(numPixel%2 != 0){
      WARN("Rounded number of pixel ("<<numPixel<<") to a multiple of 2 ("<<numPixel+1<<")");
      numPixel = numPixel + 1;
    }

    if(comp > dim_){
      WARN("Unallowed component requested.");
      return;
    }

    if(helperMatrix != NULL){
      Double numRowsStorage = helperMatrix->GetNumRows();
      Double numColsStorage = helperMatrix->GetNumCols();

      if((numRowsStorage != numPixel)||(numColsStorage != numPixel)){
        helperMatrix->Resize(numPixel);
      }

    } else {
      helperMatrix = new Matrix<Double>(numPixel,numPixel);
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

    /*
     * reset matrix
     */
    helperMatrix->Init();

    /*
     * Fill matrix
     * 1. lower triangle part -> fill with e_u
     * 2. upper triangle part -> fill with e_u
     * 3. square part -> evaluate global rotation list (only rotation states!)
     */

    /*
     * NOTE: we do not set a value of 0.5*rotState on the diagonal alpha = beta
     */
    /*
     * Part 1
     */
    for(UInt i = 0; i < numPixel/2; i++){
      for(UInt j = 0; j <= i; j++){
        (*helperMatrix)[i][j] = lastEu_[idElem][comp];
      }
    }

    /*
     * Part 2
     */
    for(UInt i = numPixel/2; i < numPixel; i++){
      for(UInt j = numPixel/2; j <= i; j++){
        (*helperMatrix)[i][j] = lastEu_[idElem][comp];
      }
    }

    /*
     * Part 3
     */
    bool twoAreas = true;
    Double upperBound;
    Double lowerBound;

    Double rot_t1,rot_b1,rot_l1,rot_r1;
    Double rot_t2,rot_b2,rot_l2,rot_r2;

    /*
     * iterators for global rotation list
     */
    std::list<RotListEntry>::iterator rotListIt;
    std::list<RotListEntry>::iterator rotListEnd = --(globRotList_[idElem].end());

    for(rotListIt = globRotList_[idElem].begin(); rotListIt != globRotList_[idElem].end(); rotListIt++){

      twoAreas = true;

      /*
       * extract upper and lower bound defining the dimensions of the L-shape
       *
       *     one possible rotation state:
       *
       *                                   alpha = 1, beta = 0
       *         __ __ __ __ __ __ __ __ __ __v  _ upper bound 0
       *        |        |     |        |     |
       *        |  0     | 1   |  2     | 3   |  _ lower bound 0; upper bound 1
       *        |________|     |        |     |
       *        |              |        |     |
       *        |______________|        |     |  _ lower bound 1; upper bound 2
       *        |                       |     |
       *        |                       |     |
       *        |_______________________|     |  _ lower bound 2; upper bound 3
       *        |                             |
       *        |__ __ __ __ __ __ __ __ __ __|  _ lower bound 3
       *        ^
       *     alpha = 0, beta = -1
       *
       *
       *    area 0 -> no L-shape -> clip only against one rectangle -> twoAreas = false
       *                              left = -1, right = - lower bound i, top = 1, bot = lower bound
       *                              (note: the bounds are positive values; on beta you have to use - that value!)
       *    area i>0 -> L-shape -> decompose into two rectangles
       *                              1. left = -1, right = - lower bound i, top = upper bound i, bot = lower bound i
       *                              2. left = - upper bound i, right = - lower bound i, top = 1, bot = upper bound i
       *
       *  NOTE: compared to function switchingStateToBmp, we do not have to consider area 9 separately as it
       *  would have a rotation state of 0.
       *
       */
      upperBound = rotListIt->getVal();
      lowerBound = rotListIt->getLowerVal();

      if(upperBound >= 1){
        twoAreas = false;
      }

      /*
       * set coordinates of first and second area (even if the second one has zero size)
       * if upperBound == lowerBound the later functions will return 0 as update
       */
      rot_t1 = upperBound;
      rot_b1 = lowerBound;
      rot_l1 = -1.0;
      rot_r1 = -lowerBound;

      rot_t2 = 1.0;
      rot_b2 = upperBound;
      rot_l2 = -upperBound;
      rot_r2 = -lowerBound;

      Vector<Double> curRotState = rotListIt->getVecReference();
      Double currentEntry = curRotState[comp];

      /*
       * we need not to clip against any switching state here; just map rotation areas directly to
       * helper matrix;
       * However, we do not want to set the switching state (i.e. a value between -1 and 1 but instead
       * a scaled entry of the rotation vector
       * -> use slot for parameter factor to do this (as factor will be multiplied to each entry
       * which gets evaluated in mapRectangleToHelperMatrix (except for the case skipUpperDiagonal, but
       * this is not needed here)
       *
       * NOTE: mapRectangleToHelperMatrix set a factor of 0.5 for elements on alpha = beta
       * This is wanted for switching states, but not for the rotation state (otherwise we would get
       * the factor twice). However, as numPixel is an even number all elements with alpha = beta
       * lie in one of the triangular areas. Therefore, we do not have to add a special case.
       *
       */
      mapRectangleToHelperMatrix(*helperMatrix,rot_t1,rot_b1,rot_l1,rot_r1,currentEntry,false);

      if(twoAreas){
        mapRectangleToHelperMatrix(*helperMatrix,rot_t2,rot_b2,rot_l2,rot_r2,currentEntry,false);
      }

    } // rot list

    /*
     * now call output function of matrix
     */
    UInt upscaling = 1;

    if(writeOutput){
      helperMatrix->matrix2Bmp(upscaling,filename_new);
    }
  }

}//end namespace
