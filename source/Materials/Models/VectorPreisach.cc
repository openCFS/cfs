#include "VectorPreisach.hh"

#include <iostream>
#include <fstream>

/* See VectorPreisach.hh for detailed description */

namespace CoupledField
{ 

  VectorPreisach::VectorPreisach(Integer numElem, Double xSat, Double ySat,
		       Matrix<Double>& preisachWeight, Double rotationalResistance, UInt dim,  bool isVirgin)
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

    dim_ = dim;

    preisachWeights_ = preisachWeight;

    // get size of preisachWeights_
    UInt M = preisachWeights_.GetNumRows();
    UInt N = preisachWeights_.GetNumCols();

    if(M != N){
      EXCEPTION("Matrix preisachWeight has dim " << M << " x " << N << " and thus is not symmetric!");
    }

    numRows_ = M;

    if(numRows_%2 == 1){
      /*
       * odd number of rows -> currently the element with center 0,0 is treated wrong as it is not treated as
       * initally split along alpha = -beta
       */
      std::cout << "Preisach plane has an odd number of rows (M = " << numRows_ << ")." << std::endl
      << "In that case the Preisach element with center (0,0) is not treated correctly right now." << std::endl
      << "Consider to use an even number of rows instead." << std::endl;
    }

    // resolution of Preisach plane
    delta_ = 2.0/Double(numRows_);

    rotationalResistance_ = rotationalResistance;

    /*
     * Local quantities, i.e. arrays storing different values for each FE element
     */

    preisachSum_ = new Vector<Double>[numElem];

    rotationalWeightsX_ = new Matrix<Double>[numElem];
    rotationalWeightsY_ = new Matrix<Double>[numElem];
    if(dim_ == 3){
      rotationalWeightsZ_ = new Matrix<Double>[numElem];
    }
    rotationalWeightsUpdated_ = new bool[numElem];

    switchingStates_ = new Matrix<Double>[numElem];
    switchingStatesUpdated_ = new bool[numElem];

    oldXpar_ = new Matrix<Double>[numElem];

    wipedOut_ = new bool*[numElem];

    firstMin_ = new bool**[numElem];

    preMinMaxChanged_ = new bool**[numElem];

    preMinMaxValues_ = new std::list<Double>**[numElem];

    isTesting_ = !true;

    /*
     * Initialize arrays/vectors/matrices for each element
     */

    for(UInt k = 0; k < (UInt) numElem; k++){

      preisachSum_[k] = Vector<Double>(dim);
      preisachSum_[k].Init(0.0);

      /*
       * for testing: initialize rotational weights in +y diretion and deactivate updateing below
       * -> if input alternates sign in y diretion, xPar will be alternating sign, too
       * -> weights below alpha = 0 should be getting set
       *
       * otherwise initialize all rotational weights with 0
       */
      rotationalWeightsX_[k] = Matrix<Double>(numRows_,numRows_);
      rotationalWeightsX_[k].InitValue(0.0);

      rotationalWeightsY_[k] = Matrix<Double>(numRows_,numRows_);
      if(isTesting_){
        rotationalWeightsY_[k].InitValue(1.0);
      } else {
        rotationalWeightsY_[k].InitValue(0.0);
      }

      if(dim_ == 3){
        rotationalWeightsZ_[k] = Matrix<Double>(numRows_,numRows_);
        rotationalWeightsZ_[k].InitValue(0.0);
      }
      rotationalWeightsUpdated_[k] = false;

      // initialize switchingStates_
      switchingStates_[k] = Matrix<Double>(numRows_,numRows_);

      /* current version of Preisach plane consists of small squares
       * i.e. we have no diagonal line but rather a stair-case
       *
       * Initial state
       * numCols_ = numRows_ = 4            numCols_ = numRows_ = 3
       *  _ _ _ _            _ _ _
       * |0 - - -|          |0 - -|
       * |+ 0 -  |          |+ 0  |
       * |+ +    |          |+_ _ |
       * |+_ _ _ |
       *
       * note that only for an odd number of rows, we have an element with value 0 on
       * the diagonal alpha = beta!
       */

      // initialize switchingStates along staircaseline alpha = beta
      UInt half;
      if((numRows_%2) == 0){
        half = (UInt) (numRows_/2.0);
      } else {
        half = (UInt) (numRows_/2.0); // round down
      }

      // lower triangle 0,0 - 0,half - half,half
      for(UInt i = 0; i < half; i++)
      {
        for(UInt j = 0; j <= i; j++ )
        {
          switchingStates_[k][i][j] = +1;
        }
      }

      // square 0,half - numRows_-half,half - 0,numRows_ - numRows_-half,numRows_
      for(UInt i = half; i < numRows_; i++)
      {
        for(UInt j = 0; j < numRows_-half; j++)
        {
          if(i < numRows_-j-1){
            switchingStates_[k][i][j] = +1;
          }
          else if(numRows_-j-1 < i){
            switchingStates_[k][i][j] = -1;
          }
          else{
            switchingStates_[k][i][j] = 0;
          }
        }
      }

      // triangle numRows_-half,numRows_-half - numRows_,numRows_-half - numRows_,numRows_
      for(UInt i = numRows_-half; i < numRows_; i++)
      {
        for(UInt j = numRows_-half; j <= i; j++)
        {
          switchingStates_[k][i][j] = -1;
        }
      }

      // will be set later; at first keep +1/-1 or 0 -> values needed later on
//      // scale elements on diagonal to 0.5
//      for(UInt i = 0; i < numRows_; i++)
//      {
//        for(UInt j = 0; j <= i; j++)
//        {
//          if(i == j){
//            switchingStates_[k][i][j] *= 0.5;
//          }
//        }
//      }

      switchingStatesUpdated_[k] = false;

      oldXpar_[k] = Matrix<Double>(numRows_,numRows_);
      oldXpar_[k].Init();

      /*
       * we need one entry for each element with center points lying on alpha = -beta
       * -> numDiag = half (numRows_ = even)
       * -> numDiag = half + 1 (numRows_ = odd)
       *
       * Remark: use the alpha coordinate to address this boolean (e.g. alpha = -1 -> idAlpha = 0 -> wipedOut_[k][0]
       */

      UInt numDiag;
      if(numRows_%2 == 1){
        /*
         * odd case
         */
        numDiag = half+1;
      } else {
        numDiag = half;
      }

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
        wipedOut_[k][numDiag-1] = true;
      }

      firstMin_[k] = new bool*[numRows_];
      for (UInt i = 0; i < numRows_; i++){
        firstMin_[k][i] = new bool[numRows_];

        for (UInt j = 0; j <= i; j++){
          firstMin_[k][i][j] = false;
        }
      }

      preMinMaxChanged_[k] = new bool*[numRows_];
      for (UInt i = 0; i < numRows_; i++){
        preMinMaxChanged_[k][i] = new bool[numRows_];

        for (UInt j = 0; j <= i; j++){
          preMinMaxChanged_[k][i][j] = false;
        }
      }

      preMinMaxValues_[k] = new std::list<Double>*[numRows_];
      for (UInt i = 0; i < numRows_; i++){

        preMinMaxValues_[k][i] = new std::list<Double>[numRows_];
        for (UInt j = 0; j <= i; j++){

          preMinMaxValues_[k][i][j] = std::list<Double>();
        }
      }
    }

  } // end constructor

  VectorPreisach::~VectorPreisach()
  {
    /*
     * TODO: delete everything
     */
  }

  void VectorPreisach::updateRotationalWeights(Vector<Double>& Xin, Integer idx){

    // determine computational threshold
    Double X_thres = std::pow((Xin.NormL2()/Xsaturated_),rotationalResistance_);

    // get direction of Xin
    Vector<Double> e_u = Vector<Double>(Xin.GetSize());

    if(Xin.NormL2() > 0.0){
      e_u = Xin/Xin.NormL2();
    } else {
      e_u.Init(0.0);
    }

    Double alpha = -1;

    for (UInt i = 0; i < numRows_; i++){

      Double beta = -1;
      for (UInt j = 0; j <= i; j++){

        /*
         * Switch rotational operator if middle point of element (alpha+delta_/2, beta+delta_/2) is passed
         * from bottom to top (i.e. at X_thres >= alpha+delta_/2) and from right to left (i.e. at X_thres <= beta+delta_/2).
         *
         * TODO: TO BE DICUSSED
         */
        if((X_thres >= alpha+delta_/2) || (-X_thres <= beta + delta_/2)){
          rotationalWeightsX_[idx][i][j] = e_u[0];
          rotationalWeightsY_[idx][i][j] = e_u[1];
          if(dim_ == 3){
            rotationalWeightsZ_[idx][i][j] = e_u[2];
          }
        } else {
          // keep old state
        }

        beta += delta_;
      }
      alpha += delta_;
    }

    rotationalWeightsUpdated_[idx] = true;
  }

  void VectorPreisach::updateSwitchingOperators(Vector<Double>& Xin, Integer idx){

    if(rotationalWeightsUpdated_[idx] == false){
      // update the rotationalWeights first, except in testing case
      if(isTesting_ == false){
        updateRotationalWeights(Xin, idx);
      }
    }

    // go over switchingStates_ and update entries
    Double alpha = -1;
    Double xPar;
    Double value;

    for (UInt i = 0; i < numRows_; i++){

      Double beta = -1;
      for (UInt j = 0; j <= i; j++){

        // at each combination of alpha and beta determine the component of X pointing along the current rotational vector
        xPar = 0.0;
        xPar += Xin[0]*rotationalWeightsX_[idx][i][j];
        xPar += Xin[1]*rotationalWeightsY_[idx][i][j];
        if(dim_ == 3){
          xPar += Xin[2]*rotationalWeightsZ_[idx][i][j];
        }

        //normalize Xpar to Xsaturatred
        xPar /= Xsaturated_;

        /*
         * update min max list of the current element with Xpar
         * Note: this can be costly if xPar is switching between local minima and maxima
         *        however, not all elements will really be changed.
         *        -> xPar gets cut down to the range of the element alpha,beta to alpha+delta_,beta+delta_ as it
         *           makes no difference if the single element is saturated with alpha+delta or with 1.
         *           If this cut down value is the same as the previously obtained cut down value, the list
         *           corresponding to this element gets not even touched as it will keep its value.
         */
        updateMinMax(xPar, alpha, beta, idx);

        /*
         * now we evaluate the list of min and max to calculate the value of the switching operator
         */
        value = evaluatePreisachElementValue(alpha,beta,idx);

        /*
         * set value
         */
        switchingStates_[idx][i][j] = value;

        beta = beta + delta_;
      }

      alpha = alpha + delta_;
    }

    switchingStatesUpdated_[idx] = true;

  }

  Double VectorPreisach::evaluatePreisachElementValue(Double alpha,Double beta,UInt idx)
  {
    /*
     * calculate indices of element
     */
    UInt idAlpha = (UInt) ((alpha+1.0)/delta_);
    UInt idBeta = (UInt) ((beta+1.0)/delta_);

    if(preMinMaxChanged_[idx][idAlpha][idBeta] == false){
      /*
       * min/max list has not been changed, i.e. also its value cannot have changed
       * -> return old value
       */
      return switchingStates_[idx][idAlpha][idBeta];
    }

    if(preMinMaxValues_[idx][idAlpha][idBeta].empty() ){
      /*
       * currently we have no min/max list for this element
       * this can happen if the input value xPar to the function updateMinMax was equal to the
       * intial value of the element; in that case updateMinMax returns without creating a list
       * don't worry, we can just return the initial state / current state in that case
       */
      return switchingStates_[idx][idAlpha][idBeta];

    } else {
      /*
       * iterate over list of min and max
       */
      Double value = 0.0;
      UInt cnt = 0;
      Double xCur = 0.0;
      Double xOld = 0.0;

      bool minCur = firstMin_[idx][idAlpha][idBeta];

      std::list<Double>::iterator listIt;

      /*
       * check if we have a full (i.e. square) element or if we are on the diagonal alpha = beta
       */
      if(idAlpha != idBeta){
        /*
         * square
         */

        /*
         * check if we are on the diagonal alpha = -beta (idAlpha + idBeta + 1 = numRows_)
         * -> if element has not been wiped out yet, we need a special treatment!
         */
        if( (idAlpha+idBeta+1 == numRows_)&&(wipedOut_[idx][idAlpha] == false) ){
          /*
           * special treatment
           */
          /*
           * This element is initially split along the line alpha = -beta in a positve lower triangle and a negative
           * upper triangle.
           * In contrast to the std case below, we have no min/max inserted to describe this special state (as it is not
           * possible to describe the diagonal by a vertical or horizontal line). Therefore, we just assume this splitting
           * to be present and treat the min/max values as updates to this state
           */
          for(listIt = preMinMaxValues_[idx][idAlpha][idBeta].begin(); listIt != preMinMaxValues_[idx][idAlpha][idBeta].end(); listIt++){

            xCur = *listIt;

            if(minCur && (cnt == 0)){
              /*
               * start with min; here we just flip the triangle below the diagonal and right from xPar to -1
               * Note: factor 2 and 0.5 were cancelled out
               */
              value = value - ((beta+delta_ - xCur)*(beta+delta_ - xCur)) / (delta_ * delta_);
            } else if (!minCur && (cnt == 0)){
              /*
               * start with max; here we flip the triangle above the diagonal up to xPar to +1
               * Note: factor 2 and 0.5 were cancelled out
               */
              value = value + ((xCur - alpha)*(xCur - alpha)) / (delta_ * delta_);
            } else if (minCur && (cnt == 1)){
              /*
               * we continue with a min; note that due to the former switching of a triangle, we have now already a horizontal
               * line through the element; if we stay on the right hand side of the crossing point of this line with the
               * diagonal alpha = -beta, we just flip a rectangular shape;
               * if we are on the left hand side of the crossing point, we also flip a triangle shape
               *
               *
               *
               *       state after first max      xPar-beta >= (alpha+delta_-x0)
               *               _________           ____v____
               *    diagonal >|\        |         |\   !    |
               *              |  \______|<x0      |  \_!____|
               *              |         |      -> |    !xxxx|
               *              |         |         |    !xxxx|< only this part is flipped
               *              |_________|         |____!xxxx|
               *
               *       state after first max      xPar-beta < (alpha+delta_-x0)
               *               _________           _v_______
               *    diagonal >|\        |         |\!       |
               *              |  \______|<x0      | !a\_____|< additional triangle has to be set
               *              |         |      -> | !xxxxxxx|
               *              |         |         | !xxxxxxx|
               *              |_________|         |_!xxxxxxx|
               *
               *
               */
              if(xCur-beta < alpha+delta_-xOld){
                /*
                 * subtract triangle; 0.5 and 2 have cancelled out
                 */
                value = value - ( (xOld - xCur)*(xOld - xCur) ) / (delta_ * delta_);
              }
              /*
               * subtract rectangle in both cases; don't forget the 2 as we have to switch from +1 to -1
               */
              value = value - 2*( (beta+delta_ - xCur)*(xOld - alpha) ) / (delta_ * delta_);
            } else if (!minCur && (cnt == 1)){
              /*
               * we continue with a max; we have now already a vertical line through the element
               * if we stay on the below the crossing point of this line with the
               * diagonal alpha = -beta, we just flip a rectangular shape;
               * if we are above the crossing point, we also flip a triangle shape
               *
               *
               *
               *       state after first min      xPar-alpha <= (beta+delta_-x0)
               *               _________           _________
               *    diagonal >|\        |         |\        |
               *              |  \      |         |  \      |
               *              |   |     |      -> |   |_____|
               *              |   |     |         |   |xxxxx|< only this part is flipped
               *              |___|_____|         |___|xxxxx|
               *                  ^x0
               *
               *       state after first min      xPar-alpha > (beta+delta_-x0)
               *               _________           _________
               *    diagonal >|\        |         |\ _______|
               *              |  \      |         | \axxxxxx|< additional triangle has to be set
               *              |   |     |      -> |   |xxxxx|
               *              |   |     |         |   |xxxxx|
               *              |___|_____|         |___|xxxxx|
               *                  ^x0
               *
               */
               if(xCur - alpha > (beta+delta_ - xOld)){
                 /*
                  * add triangle; 0.5 and 2 have cancelled out
                  */
                 value = value + ( (xCur - xOld)*(xCur - xOld) ) / (delta_*delta_);
               }
               /*
                * add rectangle for both cases; don't forget the 2 as we have to switch from -1 to +1
                */
               value = value + 2*( (xCur - alpha)*(beta+delta_ - xOld) ) / (delta_ * delta_);
            } else if (minCur && (cnt > 1)){
              /*
               * all further areas areas are simply switching rectangles
               */
              value = value - 2*( (beta+delta_ - xCur)*(xOld - alpha) ) / (delta_ * delta_);
            } else if (!minCur && (cnt > 1)){
              value = value + 2*( (xCur - alpha)*(beta+delta_ - xOld) ) / (delta_ * delta_);
            }

            cnt++;
            xOld = xCur;
            minCur = !minCur;
          } // end loop over list

        } else {
          /*
           * std treatment
           */
          for(listIt = preMinMaxValues_[idx][idAlpha][idBeta].begin(); listIt != preMinMaxValues_[idx][idAlpha][idBeta].end(); listIt++){

            xCur = *listIt;
            Double update = 0.0;

            if(minCur && (cnt == 0)){
              /*
               * start with min
               */
              //value = value + ( (xCur - beta) - (beta+delta_ - xCur))/delta_;
              update = ( (xCur - beta) - (beta+delta_ - xCur))/delta_;
            } else if (!minCur && (cnt == 0)){
              /*
               * start with max
               */
              //value = value + ( (xCur - alpha) - (alpha+delta_ - xCur) )/delta_;
              update = ( (xCur - alpha) - (alpha+delta_ - xCur) )/delta_;
            } else if(cnt != 0){
              if(minCur){
                update = - 2 * ((xOld - alpha) * (beta + delta_ - xCur)) / (delta_ * delta_);
                //value = value - 2 * ((xOld - alpha) * (beta + delta_ - xCur)) / (delta_ * delta_);
              } else {
                update = 2 * ((xCur - alpha) * (beta + delta_ - xOld)) / (delta_ * delta_);
                //value = value + 2 * ((xCur - alpha) * (beta + delta_ - xOld)) / (delta_ * delta_);
              }
            }
            cnt++;
            xOld = xCur;
            minCur = !minCur;
            value = value + update;

          } // end loop over list

        } // end of std treatment

      } else {
        /*
         * triangle with center on diagonal alpha = beta
         */

        /*
         * check if we are also on the diagonal alpha = -beta (idAlpha + idBeta + 1 = numRows_)
         * -> if element has not been wiped out yet, we need a special treatment!
         */
        if( (idAlpha+idBeta+1 == numRows_)&&(wipedOut_[idx][idAlpha] == false) ){
          /*
           * special treatment
           * TODO: implement me!
           */
        } else {
          /*
           * std treatment
           */

          for(listIt = preMinMaxValues_[idx][idAlpha][idBeta].begin(); listIt != preMinMaxValues_[idx][idAlpha][idBeta].end(); listIt++){

            xCur = *listIt;
            /*
             * TODO: here is a lot of potential for optimization, e.g. alpha = beta, reuse values, ...
             *        the current form was chosen for easier debugging
             */

            if(minCur && (cnt == 0)){
              /*
               * start with min
               */
              value = value - 0.5 * ( (beta + delta_ - xCur)/delta_ ) * ( (alpha + delta_ - xCur)/delta_ ); // upper triangle
              value = value + ( (xCur - beta)/delta_ ) * ( (alpha + delta_ - xCur)/delta_); // upper left rectangle
              value = value + 0.5 * ( (xCur - beta)/delta_ ) * ( (xCur - alpha)/delta_); // lower triangle
            } else if (!minCur && (cnt == 0)){
              /*
               * start with max
               */
              value = value - 0.5 * ( (beta + delta_ - xCur)/delta_ ) * ( (alpha + delta_ - xCur)/delta_ ); // upper triangle
              value = value - ( (xCur - beta)/delta_ ) * ( (alpha + delta_ - xCur)/delta_); // upper left rectangle
              value = value + 0.5 * ( (xCur - beta)/delta_ ) * ( (xCur - alpha)/delta_); // lower triangle
            } else if(cnt != 0){

              if(minCur){
                value = value - ((xOld - xCur)/delta_) * ((xOld - xCur)/delta_);
              } else {
                value = value + ((xCur - xOld)/delta_) * ((xCur - xOld)/delta_);
              }
            }
            cnt++;
            xOld = xCur;
            minCur = !minCur;
          }
        }
      }

      return value;

    } // end list not empty
  }

 Vector<Double> VectorPreisach::computeValue_vec(Vector<Double>& Xin, Integer idxElem)
  {

    if(switchingStatesUpdated_[idxElem] == false)
    {
      updateSwitchingOperators(Xin,idxElem);
    }

    Vector<Double> Yout = Vector<Double>(dim_);
    Yout.Init(0.0);

    Double weight = 0;
    Double area = delta_*delta_; // treat all areas as small squares

    // here we have to iterate over the whole Preisach plane and evaluate the integral for Yout
    for (UInt i = 0; i < numRows_; i++){

      for (UInt j = 0; j <= i; j++){
         weight = preisachWeights_[i][j]*switchingStates_[idxElem][i][j];

         // gets considered already during the calculation of the switchingStates_
//         if(i == j){
//           // entries on diagonal have only triangular areas, so scale them with an additional 0.5
//           weight *= 2;
//         }

         Yout[0] += weight*rotationalWeightsX_[idxElem][i][j];
         Yout[1] += weight*rotationalWeightsY_[idxElem][i][j];
         if(dim_ == 3){
           Yout[2] += weight*rotationalWeightsZ_[idxElem][i][j];
         }
      }
    }

    rotationalWeightsUpdated_[idxElem] = false;
    switchingStatesUpdated_[idxElem] = false;

    /*
     * reset flag for next step
     */
    for (UInt i = 0; i < numRows_; i++){
      for (UInt j = 0; j <= i; j++){
        preMinMaxChanged_[idxElem][i][j] = false;
      }
    }

    preisachSum_[idxElem] = Yout*(area*YSaturated_);

    return preisachSum_[idxElem];
  }

  void VectorPreisach::updateMinMax(Double xParUncut, Double alpha, Double beta, UInt idxElem)
  {
    /*
     * For Microswitching (i.e. input variation which stay inside a Preisach element), we need to store
     * the history separately in a map of minima and maxima, as the cell is completely overwritten each time
     * we update it (memory of old state is lost); this function updates lists according to the current state
     */
    /*
     * Get index of the element
     */
    UInt idAlpha = (UInt) ((alpha+1.0)/delta_);
    UInt idBeta = (UInt) ((beta+1.0)/delta_);
    bool isMinCur; // is current input smaller than old one?
    bool isMinHelp; // is inserted helper value a min?
    Double xPar; // the cut version of xParUncut

    /*
     * cut down xPar to the elements boundaries
     * Note: xPar can be much larger than needed to fill the element. however the exact value is not desired
     * as we need only to know up to which point the element would get filled by this value
     * if xPar defines a maximum:
     *    cut xPar to [alpha, alpha+delta]
     * if xPar defines a minimum:
     *    cut xPar to [beta, beta+delta]
     */

    /*
     * note that oldXpar_ is NOT cut down to the current element bounds
     * -> otherwise the comparison from old to new values can result in cases like:
     *  example (element alpha,beta = -1,-1, delta = 0.25):
     *    1.
     *    oldXpar = 0 (initial state)
     *    xPar = -0.5
     *    -> xPar < oldXpar -> isMin = true
     *    -> xPar > -1+delta -> xPar = -1+delta = -0.75 (as min)
     *    2.
     *    oldXpar = -0.75
     *    xPar = -0.5 (as before)
     *    -> xPar > oldXpar -> isMin = false
     *    -> xPar > -1+delta -> xPar = -0.75 (but now as max!)
     */

    if(xParUncut == oldXpar_[idxElem][idAlpha][idBeta]){
      /*
       * if input did not change since last time, we need not do anything
       */
      return;
    } else if (xParUncut > oldXpar_[idxElem][idAlpha][idBeta]){
      isMinCur = false;
    } else {
      isMinCur = true;
    }

    if(isMinCur){
      if(xParUncut < beta){
        xPar = beta;
      } else if (xParUncut > beta+delta_){
        xPar = beta+delta_;
      } else {
        xPar = xParUncut;
      }
    } else {
      if(xParUncut < alpha){
        xPar = alpha;
      } else if (xParUncut > alpha+delta_){
        xPar = alpha+delta_;
      } else {
        xPar = xParUncut;
      }
    }

    /*
     * Update list of minima and maxima
     */

    if(preMinMaxValues_[idxElem][idAlpha][idBeta].empty() ){
      /*
       * this is the case at the very beginning when no input or output has been specified yet;
       * here we set the starting values of the list
       */

      // To capture initial filling states, extract according input xParOld from filling state
      Double oldState = switchingStates_[idxElem][idAlpha][idBeta];
      Double helpXpar;

      /*
       * Distinguish between three states:
       * 1. element on diagonal alpha = -beta (excluding element (0,0))
       *      Initial state = 0.0
       * 2. element (0,0)
       *      -> note that (0,0) is the center here
       *      -> element exists only for an odd number of rows of the Preisach plane
       *      Initial state = 0.0; special treatment necessary
       * 3. other elements
       *      Initial state +/-1;
       */

      if(oldState == 0){
        /*
         * elements center is on the line alpha = -beta
         */
        /*
         * check for element with center (0,0)
         * TODO: this has to be changed, when the additional splitting along alpha = -beta shall be considered
         *       do we even have to initialize the element in that case?
         */
        UInt zeroIdx = (UInt) (numRows_/2); // rounded down

        if( ((idAlpha == idBeta)&&(idAlpha == zeroIdx)) && (numRows_%2 == 1)){
          /*
           * special treatment -> see VectorPreisach.hh
           */
          if(isMinCur){
            /*
             * input to this function specifies a minimum, so helper value has to be a maximum
             */
            helpXpar = alpha + 0.70710678118*delta_; // last term approx sqrt(2)/2*delta_
            isMinHelp = false;
          } else {
            /*
             * input to this function specifies a maximum, so helper value has to be a minimum
             */
            helpXpar = beta + 0.29289321881*delta_; // last term approx (1-sqrt(2)/2)*delta_
            isMinHelp = true;
          }
        } else {
// not needed anymore due to new treatment of alpha = -beta elements
//
//         /*
//          * std element with initial value of 0
//          */
//          if(isMinCur){
//            /*
//             * input to this function specifies a minimum, i.e. a reduction coming from the right hand side
//             * -> initialize element by maximum which fills it up to its half along a horizontal line
//             */
//            helpXpar = alpha + 0.5*delta_;
//            isMinHelp = false;
//          } else {
//            /*
//             * input to this function specifies a maximum, i.e. an increase coming from the bottom side
//             * -> initialize element by minimum which fills its right half with -1
//             */
//            helpXpar = beta + 0.5*delta_;
//            isMinHelp = true;
//          }

          /*
           * in the new treatment of the diagonal elements, we do not insert a helper min or max
           * (as it is not possible to correctly describe the diagonal splitting with it)
           * instead we directly insert the current value
           */

          /*
           * set helpXpar to any value outside the allowed range -1 to 1 and use it as flag below
           * due to the new treatment, a helpervalue is not allowed to be inserted
           */
          helpXpar = -123.0;

          /*
           * now we check if we need to insert a value at all
           */
          if( (isMinCur  == true) && (xPar == beta+delta_) ){
            /*
             * we could add this value to the list, but it will have no influence in the later computation
             * as it leaves the element untouched
             * -> do nothing here
             */
          }  else if( (isMinCur  == false) && (xPar == alpha) ){
            /*
             * same as above
             */
          } else {
            /*
             * here we add the value directly
             */
            preMinMaxValues_[idxElem][idAlpha][idBeta].push_back(xPar);
            firstMin_[idxElem][idAlpha][idBeta] = isMinCur;

            /*
             * check if input xPar is already such large, that element is completely filled with +/- 1
             * in that case its special diagonal treatment is wiped out
             *
             */
            if( (isMinCur  == true) && (xPar == beta) ){
              /*
               * completely -1
               */
              wipedOut_[idxElem][idAlpha] = true;
            } else if( (isMinCur  == false) && (xPar == alpha+delta_) ){
              /*
               * completely +1
               */
              wipedOut_[idxElem][idAlpha] = true;
            }
          }
        }
      } else if (oldState == 1.0){
        /*
         * element is completely filled with +1
         * -> insert maximum with value alpha + delta_
         */
        isMinHelp = false;
        helpXpar = alpha + delta_;
      } else if (oldState == -1.0){
        /*
         * element is completely filled with -1
         * -> insert minimum with value alpha
         */
        isMinHelp = true;
        helpXpar = beta;
      } else {
        EXCEPTION("Initial state can only be +/-1 or 0");
      }

      /*
       * add values to list and set the firstMin_ flag accordingly
       */
      if(abs(helpXpar) <= 1){
        /*
         * this condition is normally fulfilled
         * exception: elements with center on diagonal alpha = -beta will set helpXpar to a value outside the allowed
         * range as they do not need the following insertion of elements
         */

        firstMin_[idxElem][idAlpha][idBeta] = isMinHelp;
        preMinMaxValues_[idxElem][idAlpha][idBeta].push_back(helpXpar);

        if(isMinHelp == isMinCur){
          /*
           * if current input is of the same extremum type as helper entry, the larger max / the smaller min
           * is added to list; this case can only occur for completely filled elements (see above); in that case
           * the helper entry will already have the minimal min or maximal max value, so we can skip this case
           */
        } else {
          /*
           * here the input is of different type as the helper entry;
           * this is the std. case for elements with initial value 0, but can also occur if a +1 element is entered with
           * a minimum or a -1 element is entered by a maximum
           */
          preMinMaxValues_[idxElem][idAlpha][idBeta].push_back(xPar);
        }
      }

    } else {
      /*
       * list is not empty; i.e. it has had some value added before by the upper case
       */
      /*
       * iterate over current list of minima and maxima and find the right position for the current input value
       * all values after the new one get deleted from list -> wiping out property
       * Note:
       *  minima cannot wipe out maxima directly, but only if they replace another minimum. in that case, the
       *  maxima which followed the replaced minimum get removed from list (for maxima deleting minima it is similar)
       *  this property fits in almost all cases except for the case that a previously negative element gets completely
       *  swapped to positive. in that particular case, the first value of the list (the one which set the element negative)
       *  gets deleted, too (this is not necessary, as the computation would give the right values. deleting the entry
       *  will however prohibit alternating lists of +max -min +max -min and so on)
       *
       *  Special case: elements on line alpha = -beta; see below for details
       *
       */

      /*
       * check for the special case first:
       * if input = min and value = beta -> element completely set to -1; delete list, add new entry, set firstMin_ to true
       * if input = max and value = alpha+delta -> element completely set to +1; delete list, add new entry, set firstMin_ to false
       */
      if( (isMinCur  == true) && (xPar == beta) ){
        /*
         * element completely switched to -1
         */
        preMinMaxValues_[idxElem][idAlpha][idBeta].erase(preMinMaxValues_[idxElem][idAlpha][idBeta].begin(),preMinMaxValues_[idxElem][idAlpha][idBeta].end());
        preMinMaxValues_[idxElem][idAlpha][idBeta].push_back(xPar);
        firstMin_[idxElem][idAlpha][idBeta] = true;

        if(idAlpha+idBeta+1 == numRows_){
          /*
           * element on diagonal alpha = -beta -> now no special treatment needed anymore
           */
          wipedOut_[idxElem][idAlpha] = true;
        }
      } else if( (isMinCur == false) && (xPar == alpha+delta_) ){
        /*
         * element completely switched to +1
         */
        preMinMaxValues_[idxElem][idAlpha][idBeta].erase(preMinMaxValues_[idxElem][idAlpha][idBeta].begin(),preMinMaxValues_[idxElem][idAlpha][idBeta].end());
        preMinMaxValues_[idxElem][idAlpha][idBeta].push_back(xPar);
        firstMin_[idxElem][idAlpha][idBeta] = false;

        if(idAlpha+idBeta+1 == numRows_){
          /*
           * element on diagonal alpha = -beta -> now no special treatment needed anymore
           */
          wipedOut_[idxElem][idAlpha] = true;
        }
      } else if( (isMinCur  == true) && (xPar == beta+delta_) ){
        /*
         * we could add this value to the list, but it will have no influence in the later computation
         * -> do nothing here
         */
      }  else if( (isMinCur  == false) && (xPar == alpha) ){
        /*
         * same as above
         */
      } else {
        /*
         * now check list
         */
        bool loopMin_ = firstMin_[idxElem][idAlpha][idBeta];
        bool canBeDeleted = false;
        bool canBeInserted = false;
        std::list<Double>::iterator listIt;
        UInt cnt = 0;
        for(listIt = preMinMaxValues_[idxElem][idAlpha][idBeta].begin(); listIt != preMinMaxValues_[idxElem][idAlpha][idBeta].end(); ++listIt ){

          /*
           * if input and current value in the map are both minima or maxima, compare their values
           */
          if(isMinCur == loopMin_){

            if(isMinCur == true){
              // our current input is a minimum, so it deletes previous minima if its value is smaller
              if(xPar <= *listIt ){
                // our current input is a minimum and smaller than a previous minimum
                // -> delete all entries from now on, then insert xPar,xOut
                canBeDeleted = true;
                canBeInserted = true;
              }
            } else {
              // our current input is a maximum, so it deletes previous maxima if its value is larger
              if(xPar >= *listIt){
                // our current input is a maximum and larger than the previous maximum
                // -> delete all entries from now on, then insert xPar,xOut
                canBeDeleted = true;
                canBeInserted = true;
              }
            }
          }

          /*
           * list is alternating list of min and max, so we toggle the value after each iteration
           */
          loopMin_ = !loopMin_;

          if(canBeDeleted){
            /*
             * delete all values from the current one to the end of the list
             */
            if( (idAlpha+idBeta+1 == numRows_)&&(wipedOut_[idxElem][idAlpha]==false)){
              /*
               * for non-wiped out elements on the diagonal alpha = -beta we have to be little more careful
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

              if(cnt != 0){
                /*
                 * Case 1 -> we are not deleting the first entry; std procedure
                 */
                preMinMaxValues_[idxElem][idAlpha][idBeta].erase(listIt,preMinMaxValues_[idxElem][idAlpha][idBeta].end());
                break;

              } else {
                if(listIt == preMinMaxValues_[idxElem][idAlpha][idBeta].end()){
                  /*
                   * Case 2 -> we replace first entry, but list is otherwise empty -> std procedure
                   */
                  preMinMaxValues_[idxElem][idAlpha][idBeta].erase(listIt,preMinMaxValues_[idxElem][idAlpha][idBeta].end());
                  break;
                } else {
                  /*
                   * Case 3
                   */
                  Double x1 = *(listIt++); // get value of second entry
                  if(isMinCur){
                    if(xPar-beta <= alpha+delta_-x1){
                      /*
                       * Case 3a -> std
                       */
                      listIt--; // go back to beginning of list
                      preMinMaxValues_[idxElem][idAlpha][idBeta].erase(listIt,preMinMaxValues_[idxElem][idAlpha][idBeta].end());
                      break;
                    } else {
                      /*
                       * Case 3b
                       */
                      listIt--; // go back to beginning of list
                      preMinMaxValues_[idxElem][idAlpha][idBeta].erase(listIt,preMinMaxValues_[idxElem][idAlpha][idBeta].end());

                      // insert x1 as max
                      preMinMaxValues_[idxElem][idAlpha][idBeta].push_back(x1);
                      firstMin_[idxElem][idAlpha][idBeta] = false;
                      /*
                       * xPar will be inserted below
                       */
                      break;
                    }
                  } else {
                    if(xPar-alpha >= beta+delta_-x1){
                      /*
                       * Case 3c -> std
                       */
                      listIt--; // go back to beginning of list
                      preMinMaxValues_[idxElem][idAlpha][idBeta].erase(listIt,preMinMaxValues_[idxElem][idAlpha][idBeta].end());
                      break;
                    } else {
                      /*
                       * Case 3d
                       */
                      listIt--; // go back to beginning of list
                      preMinMaxValues_[idxElem][idAlpha][idBeta].erase(listIt,preMinMaxValues_[idxElem][idAlpha][idBeta].end());

                      // insert x1 as max
                      preMinMaxValues_[idxElem][idAlpha][idBeta].push_back(x1);
                      firstMin_[idxElem][idAlpha][idBeta] = true;
                      /*
                       * xPar will be inserted below
                       */
                      break;
                    }
                  }
                }
              }
            } else {
              /*
               * std procedure for elements which ly not on alpha = -beta
               */

              preMinMaxValues_[idxElem][idAlpha][idBeta].erase(listIt,preMinMaxValues_[idxElem][idAlpha][idBeta].end());
              break;
            }
          }

          cnt++;
        }

        /*
         * append current pair to end of the list
         * Note: if smaller/larger local max/min were found in the list, these values and all following once are deleted in the
         * loop above; if the current value is the largest min or the smallest max so far, nothing is deleted
         * However, in that second case we can not simply add the element as we are not allowed to put a smaller min directly
         * behind a larger min (a max has to be in between). If we would insert the min either way, it would be treated as a max later
         * as the list is evaluated as an alternating list of min and max!
         */
        if(canBeInserted == false){
          /*
           * so far we did not find a larger min or a smaller max in the list, so that we have a place for the current value.
           * however, there is another possibility for xPar to enter the list -> as its last value
           * as said above, this is only possible, though, if the current last value is of the oppositve extremum type
           */
          if(loopMin_ == isMinCur){
            /*
             * if the last element of the list is a maximum, loopMin_ will say true as it toggles at the end of the loop
             * if the current value is a min now (isMinCur == true) it can be append to the list
             */
            canBeInserted = true;
          }
        }

        if(canBeInserted){
          preMinMaxValues_[idxElem][idAlpha][idBeta].push_back(xPar);
        }
      }
    }

    // mark the current Preisach element as updated
    preMinMaxChanged_[idxElem][idAlpha][idBeta] = true;

    // finally overwrite oldXpar with current UNCUT value
    // do this in each case (except tha oldXpar_ already is xPar)
    oldXpar_[idxElem][idAlpha][idBeta] = xParUncut;
  }

  Vector<Double> VectorPreisach::getValue_vec( Integer idx )
  {
    return ( preisachSum_[idx] );
  }

}//end namespace
