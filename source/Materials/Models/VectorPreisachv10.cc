#include "VectorPreisachv10.hh"

#include <iostream>
#include <fstream>

/* See VectorPreisachv10_ListApproach.hh for detailed description */

namespace CoupledField
{ 
/*
 * BASE CLASS FUNCTIONS
 */
  VectorPreisachv10::VectorPreisachv10(Integer numElem, Double xSat, Double ySat,
         Matrix<Double>& preisachWeight, Double rotationalResistance , UInt dim, bool isVirgin,
         bool classical, Double angularDistance)
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

    tol_ = 1e-15;

    isSymmetric_ = preisachWeights_.IsSymmetric(tol_);

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

    angularDistance_ = angularDistance;

    classical_ = classical;

    if(classical_){
      std::cout << "VectorPreisach: Using classical vector model (Sutor2012)" << std::endl;
    } else {
      std::cout << "VectorPreisach: Using revised vector model (Sutor2015)" << std::endl;
    }

    /*
     * Local quantities, i.e. arrays storing different values for each FE element
     */
    preisachSum_ = new Vector<Double>[numElem_];
    for(UInt k = 0; k < numElem_; k++){
      preisachSum_[k] = Vector<Double>(dim_);
      preisachSum_[k].Init();
    }

    /*
     * Needed in context of the linesearch algorithm:
     * there we have to evaluate the hysteresis operator
     * for multiple value of eta; due to the memory property of the hysteresis operator
     * these test-steps will have a permanent effect on the output; to avoid this
     * we can use a flag which denotes if we want to overwrite or not;
     * in case of overwrite = false, we modify only a copy of the data structures and
     * store the final result in preisachSumTmp_l
     */
    preisachSumTmp_ = new Vector<Double>[numElem_];
    for(UInt k = 0; k < numElem_; k++){
      preisachSumTmp_[k] = Vector<Double>(dim_);
      preisachSumTmp_[k].Init();
    }
  }

  VectorPreisachv10::~VectorPreisachv10(){
    delete[] preisachSum_;
    delete[] preisachSumTmp_;
  }

  Vector<Double> VectorPreisachv10::evaluateNewRotationDirection(Vector<Double>& e_u_new, Vector<Double>& e_u_old, Double xVal){
    /*
     * calculates the new rotation direction for overwritten rotation states according to the
     * revised vector model from 2015
     *
     * Idea: the old rotation state e_u_old is not simply overwritten, but instead it rotates in direction of e_u_new
     *       it will however, only rotate up to an angle deltaPhi which depends on the magnitude of the current input
     *
     * Further addition:
     *       an easy axis might be given via input file; in that case we allow the material to rotate easier into that diretion
     *       an harder into the perpendicular direction -> theory not yet created
     */
    if(classical_){
      return e_u_new;
    } else {
      /*
       * Evaluate the new rotation direction in the following steps:
       * 1. check for the angular distance delta_phi which shall be kept between the new input direction e_u_new and the
       *    rotation direction to be taken e_phi
       * 2. check angle alpha between e_u_new and the currently stored direction e_u_old
       * 3. if delta_phi = 0
       *      return e_u_new
       *    else if delta_phi < alpha
       *      return e_u_old
       *    else
       *      rotate e_u_old towards e_u_new but keep an angular distance of deltaphi
       *        (or rotate e_u_new by deltaphi towards e_u_old)
       */
      Double delta_phi, alpha;

      /*
       * 0. Check old rotation state: if it is zero (i.e. unset yet), rotate completely to new direction
       */
      if(e_u_old.NormL2() < tol_){
        return e_u_new;
      }
      /*
       * 1. calculate delta_phi based on (10) Sutor2015
       * Note: original formula uses 2*abs(xVal) due to normalization to range [-0.5,0.5] instead of [-1,1];
       * xVal not normalized yet
       *
       * angularDistance_ in deg
       */
      xVal /= Xsaturated_;
      delta_phi = angularDistance_ * (1 - abs(xVal));

      /*
       * 2. get angle between e_u_new and e_u_old
       * Note: e_u_new and e_u_old already have length 1
       * Note2: delta_phi is in degree, so alpha has to be in degree, too
       */
      Double tmp = e_u_old.Inner(e_u_new);
      /*
       * due to rounding errors, the scalar product of two vectors of lenght 1
       * can be larger than 1 -> cut down to avoid NaN
       */
      if(tmp > 1.0){
        tmp = 1.0;
      } else if (tmp < -1.0){
        tmp = -1.0;
      }
      alpha = std::acos(tmp)*180/M_PI;

      /*
       * NEW: 15.2.2017
       * If alpha (the angle between new and old state is  close to zero,
       * we do not have to rotate at all as the rotation states are nearly the same)
       */
      // take much larger tolerance here (should we really go up to 1e-15°?!)
      Double angleTol = 1e-4;
      if(abs(alpha) < angleTol){
        /*
         * take new state
         */
        return e_u_new;
      }

      /*
       * 3. calculate new rotation direction depending on delta_phi and alpha
       */
      if(delta_phi < angleTol){
        /*
         * no resistance to rotation
         */
        return e_u_new;
      } else if(delta_phi >= alpha) {
        /*
         * e_u_old is already closer than the resistance angle delta_phi that should remain
         */
        return e_u_old;
      } else {
        /*
         * construct new rotation vector e_phi
         * Note: e_u_new is rotated by -delta_phi towards e_u_old
         */
        Vector<Double> e_phi = Vector<Double>(dim_);
        Matrix<Double> rotMat = Matrix<Double>(dim_,dim_);
        Double c,s;

        if(dim_ == 2){
          /*
           * in 2d we need
           *  1. direction of rotation axis (+z or -z)
           *  2. a 2x2 rotation matrix describing a rotation of -delta_phi degree around rotation axis
           */

          /*
           * in 2d rotation axis is either +z or -z; find sign via cross product of e_u_old and e_u_new but take only 3rd entry
           * Note: we build the right-hand-system e_u_old, e_u_new, e_u_old x e_u_new
           */
          Double n3;
          n3 = e_u_old[0]*e_u_new[1] - e_u_old[1]*e_u_new[0];
          n3 = n3/abs(n3);

          c = std::cos(-delta_phi/180*M_PI);
          s = std::sin(-delta_phi/180*M_PI);

          rotMat[0][0] = c;
          rotMat[0][1] = -n3*s;
          rotMat[1][0] = n3*s;
          rotMat[1][1] = c;

          e_phi = rotMat*e_u_new;

        } else {
          /*
           * in 3d we need
           *  1. rotation axis given by e_u_old x e_u_new
           *  2. a 3x3 rotation matrix describing a rotation of -delta_phi degree around rotation axis
           */
          Vector<Double> normal = Vector<Double>(dim_);
          e_u_old.CrossProduct(e_u_new,normal);
          normal = normal/normal.NormL2();

          c = std::cos(-delta_phi/180*M_PI);
          s = std::sin(-delta_phi/180*M_PI);

          /*
           * rotation matrix for a rotation around 'normal';
           * Note: e_u_old, e_u_new and normal have to have the same origin
           * (which is fulfilled here as all are vectors starting at (0,0,0))
           */
          rotMat[0][0] = normal[1]*normal[1]*(1-c) + c;
          rotMat[0][1] = normal[1]*normal[2]*(1-c) - normal[3]*s;
          rotMat[0][2] = normal[1]*normal[3]*(1-c) + normal[2]*s;

          rotMat[1][0] = normal[2]*normal[1]*(1-c) + normal[3]*s;
          rotMat[1][1] = normal[2]*normal[2]*(1-c) + c;
          rotMat[1][2] = normal[2]*normal[3]*(1-c) - normal[1]*s;

          rotMat[2][0] = normal[3]*normal[1]*(1-c) - normal[2]*s;
          rotMat[2][1] = normal[3]*normal[2]*(1-c) + normal[1]*s;
          rotMat[2][2] = normal[3]*normal[3]*(1-c) + c;

          e_phi = rotMat*e_u_new;
        }
        return e_phi;
      }
    }
  }


/*
 * MATRIX BASED IMPLEMENTATIONl
 */
  VectorPreisachv10_MatrixApproach::VectorPreisachv10_MatrixApproach(Integer numElem, Double xSat, Double ySat,
         Matrix<Double>& preisachWeight, Double rotationalResistance , UInt dim, bool isVirgin,
         bool classical, Double angularDistance)
    : VectorPreisachv10(numElem, xSat, ySat,
                        preisachWeight, rotationalResistance , dim, isVirgin,
                        classical, angularDistance)
  {

    std::cout << "Using Matrix-based approach" << std::endl;

    /*
     * Get storage for switchingStates and rotatationStates
     */
    switchingStates_ = new Matrix<Double>[numElem];
    rotationStateX_ = new Matrix<Double>[numElem];
    rotationStateY_ = new Matrix<Double>[numElem];
    rotationStateZ_ = new Matrix<Double>[numElem];

    /*
     * Initialize arrays/vectors/matrices for each element
     */
    for(UInt k = 0; k < (UInt) numElem; k++){
      switchingStates_[k] = Matrix<Double>(numRows_,numRows_);

      rotationStateX_[k] = Matrix<Double>(numRows_,numRows_);
      rotationStateX_[k].Init();

      rotationStateY_[k] = Matrix<Double>(numRows_,numRows_);
      rotationStateY_[k].Init();

      if(dim_ == 3){
        rotationStateZ_[k] = Matrix<Double>(numRows_,numRows_);
        rotationStateZ_[k].Init();
      }
    }

    /*
     * Initialize switchingStates
     */
    for(UInt k = 0; k < numElem_; k++){
      InitializeSwitchingState(k);
    }

  }

  VectorPreisachv10_MatrixApproach::~VectorPreisachv10_MatrixApproach(){
    delete[] switchingStates_;
    delete[] rotationStateX_;
    delete[] rotationStateY_;
    delete[] rotationStateZ_;
  }

  void VectorPreisachv10_MatrixApproach::InitializeSwitchingState(UInt idElem){

    /*
     * split Preisach plane along diagonal alpha=-beta in a +1 and -1 part
     */
    for(UInt i = 0; i < numRows_; i++){

      for(UInt j = 0; j <= i; j++){

        if(i+j+1 == numRows_){
          // on diagonal alpha = -beta -> leave at 0
        } else if(i+j+1 < numRows_){
          // below diagonal alpha = -beta -> +1
          switchingStates_[idElem][i][j] = +1;
        } else {
          // above diagonal alpha = -beta -> -1
          switchingStates_[idElem][i][j] = -1;
        }
        if (i == j){
          // on alpha = beta divide value by 2
          switchingStates_[idElem][i][j] = switchingStates_[idElem][i][j]/2.0;;
        }
      }
    }
  }

  void VectorPreisachv10_MatrixApproach::UpdateSwitchingStates(Vector<Double>& u_in, UInt idElem){
    /*
     * Update switching states from current input u_in; update only for FE element with index idElem
     * -> has to be called AFTER UpdateRotationStates
     */
    Vector<Double> curState = Vector<Double>(dim_);
    Double alpha,betaNext,xPar;

    /*
     * iterate over switching states
     */
    for(UInt i = 0; i < numRows_; i++){
     alpha = -1 + i*delta_;

     for(UInt j = 0; j <= i; j++){
       betaNext = -1 + (j+1)*delta_;

       /*
        * get corresponding rotation state
        */
       curState[0] = rotationStateX_[idElem][i][j];
       curState[1] = rotationStateY_[idElem][i][j];
       if(dim_ == 3){
         curState[2] = rotationStateZ_[idElem][i][j];
       }

       /*
        * calcualte xPar, which is responsible for the setting process
        */
       xPar = u_in.Inner(curState);

       xPar /= Xsaturated_;

       /*
        * check if update is needed
        */
       Double fac = 1.0;
       if(i == j){
         fac = 0.5;
       }

       if(xPar > alpha){
         switchingStates_[idElem][i][j] = fac;
       } else if (xPar < betaNext){
         switchingStates_[idElem][i][j] = -fac;
       }
     }
    }
  }

  void VectorPreisachv10_MatrixApproach::UpdateRotationStates(Double XThres, Double xVal, Vector<Double>& e_u_new, UInt idElem){
    /*
     * Update rotation states from current input u_in = xVal * e_u_new; change of rotation state indicated by xThres
     * Update only for FE element with index idElem
     */

    Vector<Double> curState = Vector<Double>(dim_);
    Vector<Double> newState;
    Double alpha,betaNext;

    /*
     * iterate over rotation states
     */
    for(UInt i = 0; i < numRows_; i++){
      alpha = -1 + i*delta_;

      for(UInt j = 0; j <= i; j++){
        betaNext = -1 + (j+1)*delta_;

        /*
         * check if update is necessary
         */
        if(classical_){
          /*
           * rotation state changes, if xThres > alpha or -xThres < beta+delta_, with alpha and beta being the
           * bottom left coordinates of each element
           * -> skip setting of rotation state, if XThres < alpha AND -xThres > beta+delta
           */
          if((XThres < alpha)&&(-XThres > betaNext)){
            continue;
          }
        } else {
          /*
           * rotation state changes, if xThres > alpha AND -xThres < beta+delta_
           * -> skip setting of rotation state, if xThres < alpha OR -xThres > beta+delta_
           */
          if((XThres < alpha)||(-XThres > betaNext)){
            continue;
          }
        }

        curState.Init();

        /*
         * extract current rotation state
         */
        curState[0] = rotationStateX_[idElem][i][j];
        curState[1] = rotationStateY_[idElem][i][j];
        if(dim_ == 3){
          curState[2] = rotationStateZ_[idElem][i][j];
        }

        /*
         * compute new rotation direction
         */
        newState = evaluateNewRotationDirection(e_u_new, curState, xVal);

        /*
         * store new state
         */
        rotationStateX_[idElem][i][j] = newState[0];
        rotationStateY_[idElem][i][j] = newState[1];
        if(dim_ == 3){
          rotationStateZ_[idElem][i][j] = newState[2];
        }
      }
    }
  }

  Vector<Double> VectorPreisachv10_MatrixApproach::computeValue_vec(Vector<Double>& u_in, Integer idElem, bool overwrite, bool overwriteDirection){

    Matrix<Double> switchingStatesSingleElemBAK;
    Matrix<Double> rotationStateXSingleElemBAK;
    Matrix<Double> rotationStateYSingleElemBAK;
    Matrix<Double> rotationStateZSingleElemBAK;

    if(overwrite == false){
      /*
       * get copies of the data structures
       */
      switchingStatesSingleElemBAK = switchingStates_[idElem];
      rotationStateXSingleElemBAK = rotationStateX_[idElem];
      rotationStateYSingleElemBAK = rotationStateY_[idElem];
      rotationStateZSingleElemBAK = rotationStateZ_[idElem];
    }

    /*
     * Determine the current rotational threshold
     */
    Double X_thres;
    if(classical_){
      X_thres = std::pow((u_in.NormL2()/Xsaturated_),rotationalResistance_);
    } else {
      X_thres = u_in.NormL2()*(rotationalResistance_/Xsaturated_);
    }

    /*
     * Get current direction
     */
    Vector<Double> e_u = Vector<Double>(u_in.GetSize());
    Double xVal = u_in.NormL2();

    if(xVal > tol_){
      e_u = u_in/xVal;
    } else {
      e_u.Init(0.0);
    }

    /*
     * Update rotation states
     */
    if(overwriteDirection){
      UpdateRotationStates(X_thres, xVal, e_u, idElem);
    }

    /*
     * Update switching states
     */
    UpdateSwitchingStates(u_in, idElem);

    /*
     * Storage for return values
     */
    Vector<Double> retVec = Vector<Double>(dim_);
    retVec.Init();

    /*
     * Storage for element value
     */
    Double x = 0;
    Double y = 0;
    Double z = 0;
    Double s = 0;

    /*
     * iterate over matrices and multiply entries together
     */
    for(UInt i = 0; i < numRows_; i++){

      for(UInt j = 0; j <= i; j++){
         s = switchingStates_[idElem][i][j];
         s *= preisachWeights_[i][j];
         x += s*rotationStateX_[idElem][i][j];
         y += s*rotationStateY_[idElem][i][j];
      }
    }

    retVec[0] = x;
    retVec[1] = y;

    if(dim_ == 3){
      for(UInt i = 0; i < numRows_; i++){

        for(UInt j = 0; j <= i; j++){
           s = switchingStates_[idElem][i][j];
           s *= preisachWeights_[i][j];
           z += s*rotationStateZ_[idElem][i][j];
        }
      }
      retVec[2] = z;
    }

    if(overwrite == false){
      /*
       * store to tmp array
       */
      preisachSumTmp_[idElem] = retVec*(YSaturated_*delta_*delta_);

      /*
       * reload backup
       */
      switchingStates_[idElem] = switchingStatesSingleElemBAK;
      rotationStateX_[idElem] = rotationStateXSingleElemBAK;
      rotationStateY_[idElem] = rotationStateYSingleElemBAK;
      rotationStateZ_[idElem] = rotationStateZSingleElemBAK;

      return preisachSumTmp_[idElem];
    } else {
      preisachSum_[idElem] = retVec*(YSaturated_*delta_*delta_);

      return preisachSum_[idElem];
    }
 }

  void VectorPreisachv10_MatrixApproach::switchingStateToBmp(UInt numPixel, std::string filename, UInt idElem, bool overLayWithRotState)
  {
    /*
     * NEW: rotation state is evaluated along with the switching states if overLayWithRotState is true
     * in the old versions, the rotation state was evaluated separately although we have to iterate over the
     * rotation list to evaluate the switching lists
     */

    if(numPixel < 2){
      WARN("Image should have more than 2 x 2 pixel");
      return;
    }

    if(numPixel%2 != 0){
      WARN("Rounded number of pixel ("<<numPixel<<") to a multiple of 2 ("<<numPixel+1<<")");
      numPixel = numPixel + 1;
    }

    /*
     * Calculate upscaling factor
     */
    UInt upscaling = (UInt) floor(numPixel/numRows_);

    /*
     * now call output function of matrix
     */

    if(overLayWithRotState == true){
      UInt version = 2;

      if(version == 1){

        for(UInt comp = 0; comp < dim_; comp++){

          std::stringstream stream;
          std::string filename_new;
          if(comp == 0){
            stream << "x-" << filename;
            filename_new = stream.str();
            switchingStates_[idElem].matrix2Bmp(upscaling,filename_new,&rotationStateX_[idElem]);
          } else if(comp == 1){
            stream << "y-" << filename;
            filename_new = stream.str();
            switchingStates_[idElem].matrix2Bmp(upscaling,filename_new,&rotationStateY_[idElem]);
          } else {
            stream << "z-" << filename;
            filename_new = stream.str();
            switchingStates_[idElem].matrix2Bmp(upscaling,filename_new,&rotationStateZ_[idElem]);
          }
        }

      } else if(version == 2) {
       /*
        * New way of outputting matrix:
        *   encode rotation state in color: 0° = red; 120° = blue; 240° = green; angles in between colored as a mix
        *   encode switching state as sign and amplitude: negative switching -> angle + 180°; absvalue < 1 -> scale final colorcombination
        *
        *   -> only for 2D rotstates, as the z component is not considered
        */

        std::stringstream stream;
        stream << "xy-" << filename;
        std::string filename_new = stream.str();

        switchingStates_[idElem].matrix2Bmp_v2(upscaling,filename_new,&rotationStateX_[idElem],&rotationStateY_[idElem]);
      }
    } else {
      switchingStates_[idElem].matrix2Bmp(upscaling,filename,NULL);
    }
  }

  /*
   * LIST BASED IMPLEMENTATIONl
   */
  VectorPreisachv10_ListApproach::VectorPreisachv10_ListApproach(Integer numElem, Double xSat, Double ySat,
     Matrix<Double>& preisachWeight, Double rotationalResistance , UInt dim, bool isVirgin,
     bool classical, Double angularDistance)
      : VectorPreisachv10(numElem, xSat, ySat,
                      preisachWeight, rotationalResistance , dim, isVirgin,
                      classical, angularDistance)
    {

      std::cout << "Using List-based approach" << std::endl;

      globRotList_ = new std::list<RotListEntryv10>[numElem_];
      for(UInt k = 0; k < numElem_; k++){
        globRotList_[k] = std::list<RotListEntryv10>();
        Initialize_GlobalRotationList(globRotList_[k]);
      }

      /*
       * lowerTriangleValue_ and lastEu_ only needed for classical_ model
       */
      lowerTriangleValue_ = 0.0;

      if(classical_){
        Evaluate_LowerTriangle();
      }

      lastEu_ = new Vector<Double>[numElem_];
      for(UInt k = 0; k < numElem_; k++){
        lastEu_[k] = Vector<Double>(dim_);
        lastEu_[k].Init();
      }
    }

  VectorPreisachv10_ListApproach::~VectorPreisachv10_ListApproach(){
    delete[] preisachSum_;
    delete[] globRotList_;
    delete[] lastEu_;
  }

  void VectorPreisachv10_ListApproach::Initialize_GlobalRotationList(std::list<RotListEntryv10>& usedList){
    /*
     * Make sure that list is empty
     */
    usedList.clear();

    /*
     * Initialize list with a maximum of value 0 and 0-vector
     */
    Vector<Double> zeroVec = Vector<Double>(dim_);
    zeroVec.Init(0);

    /*
     * Create new entry for rotlist; the contained switching list is empty
     */
    std::list<ListEntryv10> switchingList = std::list<ListEntryv10>();
    RotListEntryv10 initEntry = RotListEntryv10(0.0, 0.0, zeroVec, switchingList, 0.0, false, false);

    /*
     * For classical model, we insert a minimum of value 0
     * -> evaluation will lead to an upper triangle which is completely -1 and a split upper square
     *    together with the +1 lower triangle we get a total state of 0 (assuming symmetric weights!
     */
    if(classical_ == true){
      initEntry.getListReference().push_back(ListEntryv10(0,true,false));
    }

    usedList.push_back(initEntry);
    //lastXpar_[idElem] = 0.0;
  }

  void VectorPreisachv10_ListApproach::Update_GlobalRotationList(Double xThres, Double xVal, Vector<Double> e_u, std::list<RotListEntryv10>& usedList,bool overwriteDirection){
    /*
     * function for updating the global rotation list with an entry pair xThres,e_u
     * furthermore, the magnitude xVal of the original input vector u_in is passed to update the switching lists
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

    if(usedList.empty()){
      EXCEPTION("List may not be empty at this point");
    }

    /*
     * Step 1. Update rotation list by inserting new entries (if needed!)
     */
    std::list<RotListEntryv10>::iterator listIt;
    /*
     *  Note that the function list::insert adds the new entry before position listIt
     * (i.e. if listIt = list.end() it will be inserted before the end)
     */
    std::list<RotListEntryv10>::iterator insertPos = usedList.end();
    std::list<RotListEntryv10>::iterator listStart = usedList.begin();
    std::list<RotListEntryv10>::iterator listEnd = --(usedList.end());

    Double curVal;
    /*
     * lowerBound is important for RotListEntries
     * the smallest possible entry (the one which is take if an entry is appended to the end of the list)
     * is
     */
    Double lowerBound = 0;
    /*
     * Prevstate: rotation state of previous list entry
     * e_u_old: rotation state of current list entry but from last timestep
     */
    Vector<Double> prevState;
    Vector<Double> e_phi = Vector<Double>(dim_);
    Vector<Double> e_u_old;
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

    // for classical_ we have to add entries for xThres = 0.0 as this entry would influence
    // the result of the upper triangle!
    // -> already done during Initialize_GlobalRotationList -> no further insert needed
    if(xThres <= 0.0){
      xThres = 0.0;
      needsInsert = false;
    }

//    std::cout << "##########################" << std::endl;
//    std::cout << "GlobalRotationList pre insert" << std::endl;
//
//    for(listIt = globRotList_[idElem].begin(); listIt != globRotList_[idElem].end(); listIt++){
//    std::cout << listIt->ToString() << std::endl;
//    }

    int cntInner = 0;

    if(overwriteDirection){
      /*
       * Update rotation states
       */

      for(listIt = usedList.begin(); listIt != usedList.end(); listIt++){
        cntInner++;
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

              // ???
  //            if(insertPos == globRotList_[idElem].end()){
  //
  //            }

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

          if(classical_){
            /*
             * classical model knows no angular distance -> full rotation is performed
             */
            listIt->setVec(e_u);
          } else {
            /*
             * here we do not set the state to e_u directly but rotate it towards e_u;
             * Note that each previously stored rotation state gets rotated to a different e_phi!
             */
            e_u_old = listIt->getVecReference();
            e_phi = evaluateNewRotationDirection(e_u,e_u_old,xVal);

            listIt->setVec(e_phi);
          }
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
            /*
             * insert pos will already be pointing to the end of the list in this case (see definition above)
             */
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
         * 3. for Sutor2015: rotation state of new entry has to be created from rotating the direction of the (partially)
         *    overlaid entry towards the current direction
         */

        std::list<ListEntryv10> newList = std::list<ListEntryv10>();
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
           * e_u_old is the rotation state which will rotate towards e_u_new
           */
          e_u_old = insertPos->getVecReference();

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
          if(classical_){
            /*
             * add minimum of value 0 -> the same as in initialize_globalRotationList
             * No dummy anymore!
             */
            newList.push_back(ListEntryv10(0,true,false));
            wasWipedOut = false;
            startCnt = 0;
          }
          /*
           * new previous rotation state available; e_u_old = 0-vector
           */
          e_u_old = Vector<Double>(dim_);
        }

        /*
         * Important: do not forget to insert rotated from evaluateNewRotationDirection instead of e_u!
         * As we insert at the end of the list, we overwrite (at least partially) the last rotentry
         * -> needed rotation direction of partially overlapped rotation state -> see above
         */
        e_phi = evaluateNewRotationDirection(e_u,e_u_old,xVal);

        usedList.insert(insertPos,RotListEntryv10(xThres,lowerBound,e_phi,newList,lastXpar,false,false,wasWipedOut,startCnt));
      }
    } // if overwriteDirection = true

    /*
     * Step 2. For each entry in rotation list, update the switching list with the value xPar
     * -> this has to be done even if overwriteDirection = false
     */
    Double xPar;
    UInt updated;
    Vector<Double> rotState;
    /*
     * reset list iterator to the end; list may be extended!
     */
    listEnd = --(usedList.end());

    for(listIt = usedList.begin(); listIt != usedList.end(); listIt++){

      rotState = listIt->getVecReference();

      xPar = rotState.Inner(e_u)*xVal;// = rotState.Inner(u_in);
      /*
       * Normalize to Xsaturated
       */
      xPar /= Xsaturated_;

      if( abs(xPar) < 1e-15 ){
        xPar = 0.0;
      }

      /*
       * we need to pass lastXpar to list
       */
      bool isLastRotEntry = false;
      if(listIt == listEnd){
        /*
         * last entry of rotation list -> switching state for classical evaluation can extend into upper triangle area!
         *
         */
        isLastRotEntry = true;
      }

      Rectangle bbox = Rectangle(0,0,0,0);
      getBoundingBoxFromRotEntry(listIt, bbox, isLastRotEntry);

      updated = Update_SwitchingList(listIt->getListReference(),xPar,listIt->getLastLocalXpar(), bbox, listIt->wasListWipedOut(),isLastRotEntry);

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
         * (for exmplanation of startCnt see class RotListEntryv10)
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

    //std::cout << "##########################" << std::endl;
//    std::cout << "GlobalRotationList pre merge and simplify" << std::endl;
//
//    for(listIt = globRotList_[idElem].begin(); listIt != globRotList_[idElem].end(); listIt++){
//    std::cout << listIt->ToString() << std::endl;
//    }

    /*
     * Step 3. Merge adjacent rotList entries
     */
    Simplify_GlobalRotationList(usedList);

    /*
     * Step 4. Simplify local switching lists
     */
    Simplify_LocalSwitchingLists(usedList);

//    static int cnt = 0;
//    std::cout << "GlobalRotationList after merge and simplify -- step " << ++cnt << std::endl;
//
//    for(listIt = globRotList_[idElem].begin(); listIt != globRotList_[idElem].end(); listIt++){
//    std::cout << listIt->ToString() << std::endl;
//    }
//    std::cout << "##########################" << std::endl;

  }

  UInt VectorPreisachv10_ListApproach::Update_SwitchingList(std::list<ListEntryv10>& list, Double newEntry, Double lastXpar, Rectangle boundingBox, bool wasWipedOut, bool lastRotEntry){
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
     *  NEW VERSION (7.8.2016):
     *    Current (v7) wiping out rule:
     *      1. min override all larger minima and all maxima which followed the deleted minima
     *      2. max override all smaller maxima and all minima which followed the deleted maxima
     *      3. minima cannot delete maxima and maxima cannot delete minima DIRECTLY
     *
     *    Resulting state of min/max list:
     *      1. list will be alternating list of minima and maxima
     *      2. maxima are descending; minima are ascending
     *      3. overall amplitude does not have to be decreasing, i.e. list could
     *          be 1,-0.4,0.8,-0.2 but not 0.8,-0.4,1,-0.2
     *
     *    Issues with this rule:
     *      Asuume that a list begins with minimum/maximum and a new maximum/minimum with larger absolute value
     *      is added --- like (-0.6,0.4) add (0.8) -> (-0.6,0.8) ---. If we now replace the first entry of the list
     *      according to aboves wiping out rule --- e.g. (-0.6,0.8) add (-0.7) -> (-0.7) ---, we might
     *      loose important information.
     *
     *                  -0.7                -0.7
     *               ____v____           ____v____
     *              |\   !  ¦ |         |\   !    |
     *              |__\_!__¦_|0.8      |  \_!    |-> note the step! it cannot recreated without 0.8
     *              |    \  ¦ |      -> |   ^!    |
     *              |    ! \¦ |         |    !    |
     *              |____!__¦\|         |____!____|
     *                     -0.6
     *
     *      Currently (v<8) we checked for cases as the one above and kept important maxima/minima if needed.
     *      E.g. we would have (-0.7,0.8) then. This is, however, not the only problem. During evaluation,
     *      a list like (-0.7,0.8) will not work with the used evaluation scheme, as that scheme assumes
     *      the first entry of the list to have the largest absolute value.
     *      Take the following example:
     *
     *           MinMax[1] = min
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
     *      Here, the first entry abs(MinMax[0]) is larger than abs(MinMax[1]) and thus allows
     *      to decompose the square into subareas like
     *
     *       area 0 (square!)
     *       _v_______________
     *      | \  | 1 |    |   |
     *      |___\|___| 3  |   |
     *      |   2    |    |5  |
     *      |________|____|   |
     *      |     4       |   | < std decomposition if starting with maximum
     *      |_____________|___|
     *      |        6        |
     *      |_________________|
     *
     *      If abs(MinMax[1]) is larger than abs(MinMax[0]), however, we get the following issues
     *
     *        here is MinMax[1] = min
     *       __v______________
     *      | \|1|        |   |
     *area0>|__|\|     3  |   |
     *no    |2 |          |5  |
     *square|__|__________|   | < std decomposition does not work
     *      |     4       |   |
     *      |_____________|___|
     *      |        6        |
     *      |_________________|
     *
     *      The decomposition above cannot be treated by the used evaluation anymore as
     *        - area0 is not square anymore (which is important as we than can leave it out assuming
     *          symmetric weights)
     *        - area1 is not completely -1
     *        - area2 is too small
     *        - area3 is no rectangle
     *      This issue can be solved by checking for the special case and than tweaking the decomposition.
     *
     *      However, there seems to be a much easier way.
     *
     *    New, simpler (original?) approach:
     *      Extend wiping out rule by a forth point:
     *        4. if abs(input) is larger than abs(first list entry) -> wipe out list completely and set input
     *           as first value
     *           -> a minimum may override a maximum and vice versa, but only in that special case
     *
     *       With respect to aboves list-example, we get:
     *        (-0.6,0.4) add (0.8) -> (0.8)
     *        (0.8) add (-0.7) -> (0.8,-0.7)
     *            *
     *                  -0.7                -0.7
     *               ____v____           ____v____
     *              |\   !  ¦ |         |\   !    |
     *              |__\_!__¦_|0.8      |__\_!____|0.8
     *              |    \  ¦ |      -> |   ^!    |       < state can now be recreated
     *              |    ! \¦ |         |    !    |
     *              |____!__¦\|         |____!____|
     *                     -0.6
     *
     *      And regarding the evaluation scheme:
     *
     *         original MinMax[0]      the former MinMax[1] is now MinMax[0]
     *        ____v____________          __v______________
     *       | \|1|        |   |        |_\|          |   |
     *       |__|\|     3  |   |        |  |     2    |   |
     *       |2 |          |5  |        |1 |          |4  |
     *       |__|__________|   |  ->    |__|__________|   | < std scheme for a list starting
     *       |     4       |   |        |     3       |   |   with a minimum!
     *       |_____________|___|        |_____________|___|
     *       |        6        |        |        5        |
     *       |_________________|        |_________________|
     *
     *      Note that the former starting value MinMax[0] is not needed for the transformed state.
     *
     *    IMPORTANT REMARKS:
     *      1. The new update/wiping-out rule is only relevant for the FIRST ENTRY of the list!
     *      2. The new treatment is only for the upper square!
     *          Reason: the upper triangle is not symmetric to alpha = -beta! Here it does not work
     *                and makes problems instead!
     *          Note: in the global sense it would be ok, but as we split along beta = 0 we have to
     *                distinguish.
     *
     *
     */

    // not needed anymore as we clip against bounding box
//    Double alpha,beta,delta;
//
//    /*
//     * New treatment: switching state is always over the whole Preisach plane; it will later get clipped to
//     * fitting rotation states
//     * in that case, the min-max list can have values from -1 to +1 for alpha and beta
//     */
//    alpha = -1.0;
//    beta = -1.0;
//    delta = 2.0;

    std::list<ListEntryv10>::iterator listIt;

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
     * get outer boundings of the corresponding rotation state
     */
    Double l,r,t,b;
    boundingBox.getBounds(l,r,t,b);

//    std::cout << "xPar: " << newEntry << std::endl;
//    std::cout << "lastXpar: " << lastXpar << std::endl;
//    std::cout << "left/right: " << l << " / " << r << std::endl;
//    std::cout << "isLastRotEntry: " << lastRotEntry << std::endl;

    /*
     * Check if lists gets completely overwritten (i.e. the saturated input value exceeds the bounds of the region
     */
    if(state > 0){
      /*
       * Element is to be set by maximum
       *
       * NEW: restrict to bounding box of rotation entry!
       * Background: Assume a rotation entry inside the bounding box l,r,t,b
       *  __ __ __ _t __ __ __ __
       * |               |       |
       * |               |       |
       * |               |       |
       * |               |       |
       * l               |       r
       * |               |       |
       * |__ __ __ __ __ |       |
       * |       rotstate        |
       * |                       |
       * |__ __ __ _b __ __ __ __|
       *
       * Stairs of the switching state contribute only, if the lines pass through the bounding box
       * -> if a line lies above t, it makes no difference if its value is >t or exactly t; in both cases the whole switching states
       *      is set to +1
       * -> in a similar way it behaves for r;
       * -> for l and b, it makes no difference, if input is smaller or exactly that value
       *
       *
       *  -> restrict to [b,t]
       */
      if(newEntry >= t){
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
        newEntry = t;
        /*
         * add maximum to list
         */
        list.push_back(ListEntryv10(newEntry,false));
        /*
         * elements on alpha = -beta also loose their special state
         */
        return 2;

      } else if(newEntry <= b){
        /*
         * new entry will have no effect
         */
        return 0;
      } else if(appendDirectly){
        list.push_back(ListEntryv10(newEntry,false));
        /*
         * return code 3: list was not wiped, but first entry was set/replaced
         */
        return 3;
      } else if(classical_){
        /*
         * Works only on the L-shapes -> only for classical_
         */
        /*
         * for upper square, check value of lowerVal
         * if a maximum to be added has a value smaller than lowerVal (bottom of the rotation area), the resulting area will have
         * no effect at all onto the switching state
         * -> do not add to list
         * (similar to the check above with newEntry <= alpha)
         */
        if(newEntry <= b){
          return 0;
        }
      }
      /*
       * else -> see below
       */
    } else {
      /*
       * Element is to be set by a minimum -> restrict to [l,r]
       */
      if(newEntry < l){
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
        newEntry = l;
        /*
         * add minimum to list
         */
        list.push_back(ListEntryv10(newEntry,true));
        /*
         * elements on alpha = -beta also loose their special state
         */
        /*
         * return 2 -> list was updated and wiped
         */
        return 2;

      } else if(newEntry >= r){
        /*
         * new entry will have no effect
         */
        return 0;
      } else if(appendDirectly){
        list.push_back(ListEntryv10(newEntry,true));
        /*
         * return code 3: list was not wiped, but first entry was set/replaced
         */
        return 3;
      } else if(classical_){
        /*
         * Works only on the L-shapes -> only for classical_
         * -> make sure when calling this function, that it does set lower value to
         */
        /*
         * for upper square, check value of lowerVal
         * if a minimum to be added has a value larger than -lowerVal (right boundary), the resulting area will have
         * no effect at all onto the switching state
         * -> do not add to list
         * (similar to the check above with newEntry >= beta+delta)
         */
        if((newEntry >= r)&&(lastRotEntry == false)){
          /*
           * Note: Unlike the maxima (which are terminated by the lower bound or the beta axis),
           * the minima can NOW reach into the upper triangle region for the classical approach.
           * This is only the case for the last entry of the rot list. In that case we can allow minima
           * which do not stop at the lower value (which would be 0) but which can extend up to beta = +1
           */
          //NOTE: can only be equal r here, as the corresponding bounding box for the last entry goes up to +1
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
    std::list<ListEntryv10>::iterator lastEntry = --(list.end());

    /*
     * if we came to this point, we have to iterate through the list and check if value shall be included
     */
    bool canBeInserted = false;
    bool canBeDeleted = false;

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
    ListEntryv10 helperEntry = ListEntryv10(0,false);

    for(listIt = list.begin(); listIt != list.end(); listIt++){

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
        } else {
          /*
           * New in version 8
           * -> first entry of list can be overwritten by an extremum of different type if abs(newValue) >= abs(oldValue)!
           * Only valid for first entry!
           */
          if((cnt == 0)&&( (abs(newEntry)-abs(listVal)) >= 0)){
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
        } else {
          /*
           * New in version 8
           * -> first entry of list can be overwritten by an extremum of different type if abs(newValue) >= abs(oldValue)!
           * Only valid for first entry!
           */
          if((cnt == 0)&&( (abs(newEntry)-abs(listVal)) >= 0)){
            canBeInserted = true;
            canBeDeleted = true;
          }
        }
      }

      if(canBeDeleted){
        /*
         * delete all entries of the list starting with the current one
         */
        list.erase(listIt,list.end());

        /*
         * set flag denoting that the first entry of the list got deleted, too.
         * -> after this function, a new first entry is to be set
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

      if(state > 0){
        /*
         * insert maximum
         */
        list.push_back(ListEntryv10(newEntry,false));
      } else if(state < 0){
        /*
         * insert minimum
         */
        list.push_back(ListEntryv10(newEntry,true));
      }
      /*
       * mark element as updated
       */
      if(firstEntrySet == true){
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

  Vector<Double> VectorPreisachv10_ListApproach::computeValue_vec(Vector<Double>& u_in, Integer idElem, bool overwrite,bool overwriteDirection){
    /*
     * Determine the current rotational threshold
     */
    Double X_thres;
    if(classical_){
      X_thres = std::pow((u_in.NormL2()/Xsaturated_),rotationalResistance_);
    } else {
      X_thres = u_in.NormL2()*(rotationalResistance_/Xsaturated_);
    }

    /*
     * Get current direction
     */
    Vector<Double> e_u = Vector<Double>(u_in.GetSize());
    Double xVal = u_in.NormL2();

    if(xVal > tol_){
      e_u = u_in/xVal;
    } else {
      e_u.Init(0.0);
    }

    /*
     * set value of lastEu_ (only needed for classical_ model to get the rotation information for the lowerTriangle_)
     */
    if(overwriteDirection){
      lastEu_[idElem] = e_u;
    }
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
     * Update and evaluate global rotation list
     */

//    /*
//     * check if copy works
//     */
//    std::cout << "GlobalRotationList pre updating " << std::endl;
//
//    std::list<RotListEntryv10>::iterator listIt;
//    for(listIt = globRotList_[idElem].begin(); listIt != globRotList_[idElem].end(); listIt++){
//    std::cout << listIt->ToString() << std::endl;
//    }
//    std::cout << "##########################" << std::endl;
//

    if(overwrite == true){
      /*
       * work on std data structure
       */
      //std::cout << "Work on permanent storage" << std::endl;
      Update_GlobalRotationList(X_thres, xVal, e_u, globRotList_[idElem],overwriteDirection);
      //Update_GlobalRotationList(X_thres, xVal, e_u, globRotList_[idElem],true);

      /*
       * Evaluate_GlobalRotationList checks for each element if it was changed or not and
       * reevaluates only the ones that did
       */
      Evaluate_GlobalRotationList(globRotList_[idElem], retVec);
    } else {
      /*
       * get copy of globRotList_[idElem]
       */
      //std::cout << "Working on temporal copy" << std::endl;
      std::list<RotListEntryv10> tmpList = globRotList_[idElem];

//      std::cout << "tmpList pre updating " << std::endl;
//
//      for(listIt = tmpList.begin(); listIt != tmpList.end(); listIt++){
//      std::cout << listIt->ToString() << std::endl;
//      }
//      std::cout << "##########################" << std::endl;

      /*
       * then perform all updates on that temporal list only
       */
      std::cout << "Working only on temporal storage! " << std::endl;
      //Update_GlobalRotationList(X_thres, xVal, e_u, tmpList,true);
      Update_GlobalRotationList(X_thres, xVal, e_u, tmpList,overwriteDirection);

      Evaluate_GlobalRotationList(tmpList, retVec);

//      std::cout << "tmpList post updating " << std::endl;
//
//      for(listIt = tmpList.begin(); listIt != tmpList.end(); listIt++){
//      std::cout << listIt->ToString() << std::endl;
//      }
//      std::cout << "##########################" << std::endl;
    }

//    /*
//     * check if copy works
//     */
//    std::cout << "GlobalRotationList after updating and evaluation " << std::endl;
//
//    for(listIt = globRotList_[idElem].begin(); listIt != globRotList_[idElem].end(); listIt++){
//    std::cout << listIt->ToString() << std::endl;
//    }
//    std::cout << "##########################" << std::endl;

    Yout += retVec;

    if(classical_){
      /*
       * Add value of lower triangle
       */
      Yout += e_u * lowerTriangleValue_;
    }

  //  std::cout << "###YOUT###################" << std::endl;
  //  std::cout << std::setprecision(16) << std::scientific << Yout.ToString() << std::endl;
  //  std::cout << "##########################" << std::endl;

    /*
     * in previous versions we scaled with delta_^2, too
     * this is no longer needed as we perform this step already in mapRectangleToPreisachWeights which is called
     * for the evaluation of all three parts
     */
    if(overwrite == true){
      preisachSum_[idElem] = Yout*YSaturated_;
      return preisachSum_[idElem];
    } else {
      preisachSumTmp_[idElem] = Yout*YSaturated_;
      return preisachSumTmp_[idElem];
    }
  }

  Double VectorPreisachv10_ListApproach::clipRectangleToElement(Rectangle& source, UInt idAlpha, UInt idBeta, Double delta, bool isRotState){
    /*
     * Calculates the overlapping area of a rectangle (rectT,rectB,rectL,rectR) with the element defined by alphaId, betaId
     * This calculation does the following steps:
     *  1. calculate the overlapping rectangle via clipRectangles
     *  2. if success, check if lower left corner or upper right corner lie below diagonal
     *    -> if yes, cut down so that only the right corner is still below the diagonal
     *    -> subtract triangular area
     *
     * new flag: isRotState
     * -> if the rotation state is mapped to the helper matrix (in mapRectangleToHelperMatrix), we do not want to scale its value by 0.5 if its on the diagonal alpha = beta
     *    if we would do so, we would get the factor 0.5 twice in the final images which are created from the matrix (1. from switching state,
     *    2. from rotation state)
     * -> if isRotState == True, scale full entries on alpha = beta times 2 and for partially filled elements, skip
     *    the cutting of triangle parts below the diagonal
     *
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

    Rectangle elemRect = Rectangle(elemLeft,elemLeft+delta, elemBot+delta, elemBot);
    Rectangle target = Rectangle(0,0,0,0);

    success = elemRect.clipRectangles(source,target);

    if(!success){
      return 0.0;
    }

    target.getBounds(ovL,ovR,ovT,ovB);

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

    if(isRotState == true){
      /*
       * do not subtract the part below the diagonal -> treat element as full element
       */
      triang = 0.0;
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

  void VectorPreisachv10_ListApproach::getBoundingBoxFromRotEntry(std::list<RotListEntryv10>::iterator rotListIt, Rectangle& rect, bool lastRotListEntryv10){
    /*
     * helper function returning a bounding box including all rotation areas which belong to a given rotListEntryv10
     *
     * Idea:
     *  when new values are added to the local switching list, these values can be clipped to this bounding box
     *
     */

    Double l,r,t,b;
    Double upperBound = rotListIt->getVal(); //exi
    Double lowerBound = rotListIt->getLowerVal(); //exi+1

    if(classical_ == true){
      /*
       * Rotation states form flipped L-shapes in S_U; last entry of list extends into T_U
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
       * Bounding boxes:
       *  if Entry is not last entry:
       *    l = -1; r = -lowerVal; t = +1; b = lowerVal
       *    (Bounding box includes all areas belong to this and to former rotListEntries
       *  if Entry is last entry:
       *    l = -1; r = +1; t = +1; b = 0
       *    (Bounding box includes all areas)
       *
       *
       */
      if(lastRotListEntryv10 == true){
        /*
         * last entry in list was found -> extend right up to +1
         */
        l = -1.0;
        r = 1.0;
        b = 0.0;
        t = 1.0;

      } else {
        /*
         * inner entry (also includes i=0)
         */
        l = -1.0;
        r = -lowerBound;
        b = lowerBound;
        t = 1.0;
      }
    } else {
      /*
       * Rotation states form flipped L-shapes in the whole Preisach plane
       *
       *            S_U           alpha           T_U
       *      __ __ __ __ __ __ __ |_  __ __ __ __ __ __
       *     |                     |                  /
       *     |A02                  |                /
       *     |_ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ /_ ex1
       *     |     |   A12         |            /
       *     |A01  |_ _  _ _ _ _ _ |_ _ _ _ _ /_ ex2
       *     |     |   |    A22    |        /
       *     |     |A11|           |      /
       *     |     |   |_ _ _ _ _ _| _ _/_ _ exN
       *     |     |   |      |    |  /
       *     |_____|___|_A21__|_AN_|/____________________ beta
       *     |     |   |      |   /
       *     |     |   |      | /
       *     |     |   |      /
       *     |     |   |    /
       *     |     |   |  /
       *     |     |   |/
       *     |     |  /
       *     |     |/
       *     |    /
       *     |_ /
       *           T_L
       *
       *
       * Bounding boxes:
       *  l = -upperVal, r = upperVal, t = upperVal, b = -upperVal
       *
       */
      l = -upperBound;
      r = upperBound;
      b = -upperBound;
      t = upperBound;
    }

    rect.setBounds(l,r,t,b);
  }


  bool VectorPreisachv10_ListApproach::getRectanglesFromRotEntry(std::list<RotListEntryv10>::iterator rotListIt, Rectangle& rect1, Rectangle& rect2, bool lastRotListEntryv10){
    /*
     * -encapsulate the determination of rectangular rotation areas from a rot list entry
     * -return true if two non-zero rectangles are created; false otherwise
     * -lastRotListEntryv10 needed to check for last entry
     */

    bool twoAreas = true;
    Double l,r,t,b;
    Double l2,r2,t2,b2;
    Double upperBound = rotListIt->getVal(); //exi
    Double lowerBound = rotListIt->getLowerVal(); //exi+1

    if(classical_ == true){
      /*
       * Rotation states form flipped L-shapes in S_U; last entry of list extends into T_U
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
       *     A0:
       *      L = -1, R = -ex1, B = ex1, T = 1
       *
       *     Ai1:
       *      L = -1, R = -exi+1, B = exi+1, T = exi
       *     Ai2:
       *      L = -exi, R = -exi+1, B = exi, T = 1
       *
       *     AN1:
       *      L = -1, R = exN, B = 0, T = exN
       *     AN2:
       *      L = -exN, R = 1, B = exN, T = 1
       *
       */
      if(upperBound >= 1){
        /*
         * A0 = 0; A11 = square; A12 = 0
         */
        twoAreas = false;
      }
      if(lastRotListEntryv10 == true){
        /*
         * AN -> last entry in list was found
         */
        /*
         * AN1
         */
        l = -1.0;
        r = upperBound;
        b = 0.0;
        t = upperBound;

        /*
         * AN2
         */
        l2 = -upperBound;
        r2 = 1;
        b2 = upperBound;
        t2 = 1;
      } else {
        /*
         * Ai (also includes i=0)
         */
        l = -1.0;
        r = -lowerBound;
        b = lowerBound;
        t = upperBound;

        l2 = -upperBound;
        r2 = -lowerBound;
        b2 = upperBound;
        t2 = 1.0;
      }
    } else {
      /*
       * Rotation states form flipped L-shapes in the whole Preisach plane
       *
       *            S_U           alpha           T_U
       *      __ __ __ __ __ __ __ |_  __ __ __ __ __ __
       *     |                     |                  /
       *     |A02                  |                /
       *     |_ _ _ _ _ _ _ _ _ _ _|_ _ _ _ _ _ _ /_ ex1
       *     |     |   A12         |            /
       *     |A01  |_ _  _ _ _ _ _ |_ _ _ _ _ /_ ex2
       *     |     |   |    A22    |        /
       *     |     |A11|           |      /
       *     |     |   |_ _ _ _ _ _| _ _/_ _ exN
       *     |     |   |      |    |  /
       *     |_____|___|_A21__|_AN_|/____________________ beta
       *     |     |   |      |   /
       *     |     |   |      | /
       *     |     |   |      /
       *     |     |   |    /
       *     |     |   |  /
       *     |     |   |/
       *     |     |  /
       *     |     |/
       *     |    /
       *     |_ /
       *           T_L
       *
       *
       *     A01:
       *      L = -1, R = -ex1, B = -1, T = ex1
       *     A02:
       *      L = -1, R = 1, B = ex1, T = 1
       *
       *     Ai1:
       *      L = -exi, R = -exi+1, B = -exi, T = exi+1
       *     Ai2:
       *      L = -exi, R = exi, B = exi+1, T = exi
       *
       *     AN:
       *      L = -exN, R = exN, B = -exN, T = exN
       *
       */
      if(lastRotListEntryv10 == true){
        /*
         * AN -> last entry in list was found
         */
        twoAreas = false;
        l = -upperBound;
        r = upperBound;
        b = -upperBound;
        t = upperBound;

        l2 = 0;
        r2 = 0;
        b2 = 0;
        t2 = 0;
      } else {
        /*
         * Ai (also includes i=0)
         */
        l = -upperBound;
        r = -lowerBound;
        b = -upperBound;
        t = lowerBound;

        l2 = -upperBound;
        r2 = upperBound;
        b2 = lowerBound;
        t2 = upperBound;
      }
    }

    rect1.setBounds(l,r,t,b);
    rect2.setBounds(l2,r2,t2,b2);

    return twoAreas;
  }

  void VectorPreisachv10_ListApproach::Evaluate_GlobalRotationList(std::list<RotListEntryv10>& usedList, Vector<Double>& retVec){
    /*
     * Evaluates the weighted rotation state of
     *  a) upper square S_U and upper triangle T_U (classic)
     *  b) the whole Preisach plane (revised)
     *
     * Evaluation is done in the following way:
     *
     * for rotEntry in globalRotList:
     *   calculate two rectangular regions corresponding to the rotation area (forming a L)
     *
     *   for entry in localSwitchList:
     *     get rectangular bounds from entry in switch list
     *     clip rectangle with first rotation area
     *     clip result to Preisach plane
     *
     *     clip rectangle with second rotation area
     *     clip result to Preisach plane
     *
     *     set entry of switching list to not changed
     *     sum up both parts with appropriate signs
     *
     *   set entry of rotation list to not changed
     *   multiply sum with rotation state and add product to eval state
     *
     * save evaluated state in retVec
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
     * iterators for local switching list
     */
    std::list<ListEntryv10>::iterator swListIt;
    std::list<ListEntryv10>::iterator swListStart;
    std::list<ListEntryv10>::iterator swListEnd;

    /*
     * iterators for global rotation list
     */
    std::list<RotListEntryv10>::iterator rotListIt;
    std::list<RotListEntryv10>::iterator rotListEnd = --(usedList.end());

    bool twoAreas = true;
    bool lastRotListEntryv10;

    Rectangle rotRect1 = Rectangle(0,0,0,0);
    Rectangle rotRect2 = Rectangle(0,0,0,0);
    Rectangle swRect = Rectangle(0,0,0,0);
    Rectangle overlapRect = Rectangle(0,0,0,0);

    Double sum = 0.0;
    Double tmp = 0.0;
    Double factor;
    UInt cnt = 0;
    UInt outercnt = 1;

    Vector<Double> rotState;

    for(rotListIt = usedList.begin(); rotListIt != usedList.end(); rotListIt++){

      if(rotListIt == rotListEnd){
        lastRotListEntryv10 = true;
      } else {
        lastRotListEntryv10 = false;
      }
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
         * get rotation area(s)
         */

        swListStart = rotListIt->getListReference().begin();
        swListEnd = --(rotListIt->getListReference().end());

        twoAreas = getRectanglesFromRotEntry(rotListIt, rotRect1, rotRect2,lastRotListEntryv10);

//        std::cout << "---------------------------" << std::endl;
//        std::cout << "Rotlistentry " << outercnt << std::endl;
//        Double ltmp,rtmp,ttmp,btmp;
//        rotRect1.getBounds(ltmp,rtmp,ttmp,btmp);
//        std::cout << "Rect1: " << ltmp << ", " << rtmp << ", " << ttmp << ", " << btmp << std::endl;
//        rotRect2.getBounds(ltmp,rtmp,ttmp,btmp);
//        std::cout << "Rect2: " << ltmp << ", " << rtmp << ", " << ttmp << ", " << btmp << std::endl;
        outercnt++;

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
         *  anymore (as this is the special helper area). So we startCnt in rotListEntryv10 to mark that we
         *  deleted the first entry which originally corresponded to area 0. As all areas with id > 0 are
         *  calculated by the same scheme, the absolute value of cnt is not of interest and thus it is
         *  either 0 or 1.
         */
        cnt = rotListIt->getStartCnt();
        sum = 0.0;
        tmp = 0.0;
        bool success = false;

//        std::cout << "RotRect1" << std::endl;
//        std::cout << rotRect1.ToString() << std::endl;
//
//        std::cout << "RotRect2" << std::endl;
//        std::cout << rotRect2.ToString() << std::endl;

        for(swListIt = rotListIt->getListReference().begin(); swListIt != rotListIt->getListReference().end(); ){

          /*
           * get rectangular bounds corresponding to entry of switching list (compare to Evaluate_GlobalSwitchingList)
           */
          factor = getRectangleFromSwitchingList(rotListIt->getListReference(),swListStart, swListIt, swListEnd,cnt, swRect);

          /*
           * clip resulting area against rotation area1
           */
          success = rotRect1.clipRectangles(swRect,overlapRect);

//          std::cout << "Switching" << std::endl;
//          std::cout << swRect.ToString() << std::endl;

          if(success){
            /*
             *  clip overlap to PreisachPlane and sum up
             */
            tmp = mapRectangleToPreisachWeights(overlapRect);

//            std::cout << "Overlap1" << std::endl;
//            std::cout << overlapRect.ToString() << std::endl;
//            std::cout << "Factor: " << factor << std::endl;
//            std::cout << "Value: " << tmp << std::endl;

            sum += factor*tmp;
          }

          if(twoAreas){
            /*
             * repeat the clipping steps for the second area
             */
            success = rotRect2.clipRectangles(swRect,overlapRect);

            if(success){
             /*
              *  clip overlap to PreisachPlane and sum up
              */
             tmp = mapRectangleToPreisachWeights(overlapRect);

//             std::cout << "Overlap2" << std::endl;
//             std::cout << overlapRect.ToString() << std::endl;
//             std::cout << "Factor: " << factor << std::endl;
//             std::cout << "Value: " << tmp << std::endl;

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
              /*
               * set flag upperSplitSquare to true -> get area 0 instead of area 1
               */
              factor = getRectangleFromSwitchingList(rotListIt->getListReference(),swListStart, swListIt, swListEnd,cnt, swRect,true);

              if(factor != 0){
                EXCEPTION("Something got wrong here! Check function getRectangleBounds!");
              }

              /*
               * clip resulting area against rotation area1
               */
              success = rotRect1.clipRectangles(swRect,overlapRect);

              if(success){
                /*
                 *  clip overlap to PreisachPlane and sum up
                 *  -> NOTE: here we set the flag skipUpperDiagonal to true!
                 */
                tmp = mapRectangleToPreisachWeights(overlapRect,true);
                // factor is not needed here; mapRectangleToPreisachWeights with skipUpperDiagonal == true
                // will not only skip the diagonal entries, but also consider the right signs!
                sum += tmp;
              }

              if(twoAreas){
                /*
                 * repeat the clipping steps for the second area
                 */
                success = rotRect2.clipRectangles(swRect,overlapRect);

                if(success){
                  /*
                   *  clip overlap to PreisachPlane and sum up
                   *  -> NOTE: here we set the flag skipUpperDiagonal to true!
                   */
                  tmp = mapRectangleToPreisachWeights(overlapRect,true);
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
//        std::cout << "Sum: " << sum << std::endl;
//        std::cout << "rotState: " << rotState.ToString() << std::endl;
        retVec += rotState*sum;
//        std::cout << "retVec: " << retVec.ToString() << std::endl;
      }
    } // rot list
  }

  void VectorPreisachv10_ListApproach::Evaluate_LowerTriangle(){
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
    Rectangle rect = Rectangle(-1.0,0.0,0.0,-1.0);
    lowerTriangleValue_ = mapRectangleToPreisachWeights(rect);
  }

  Double VectorPreisachv10_ListApproach::mapRectangleToPreisachWeights(Rectangle& rect, bool skipUpperDiagonal){
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
    Double t,b,l,r;
    rect.getBounds(l,r,t,b);

 //   std::cout << "L,R,T,B: " << l << "," << r << "," << t << "," << b << "," << std::endl;

    /*
     * restrict input value to valid range
     * (Preisach plane just goes from -1 to +1 in alpha and beta)
     */
    t = std::min(t,1.0);
    b = std::max(b,-1.0);
    l = std::max(l,-1.0);
    r = std::min(r,1.0);

 //   std::cout << "L,R,T,B: " << l << "," << r << "," << t << "," << b << "," << std::endl;

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
    int rowMin =  floor(b/delta_) + floor(1.0/delta_);
    int rowMax = ceil(t/delta_) + floor(1.0/delta_)-1;

    int colMin =  floor(l/delta_) + floor(1.0/delta_);
    int colMax =  ceil(r/delta_) + floor(1.0/delta_)-1;

    /*
     * NEW: use integers instead of unsigned integer
     * REASON: colMaxFull became negative (which simply would indicate no fully overlapped elemets), but
     * unsigning took it to a very large positive value and thus -> wrong result
     */

    /*
     * 2. get indices of completely overlapped elements
     */
    int rowMinFull =  ceil(b/delta_) + floor(1.0/delta_);
    int rowMaxFull = floor(t/delta_) + floor(1.0/delta_)-1;

    int colMinFull =  ceil(l/delta_) + floor(1.0/delta_);
    int colMaxFull = floor(r/delta_) + floor(1.0/delta_)-1;

//    std::cout << std::setprecision(12) << std::scientific << b << " " << t << " " << l << " " << r << std::endl;
//    std::cout << std::setprecision(12) << std::scientific << floor(b) << " " << ceil(t) << " " << floor(l) << " " << ceil(r) << std::endl;
//    std::cout << std::setprecision(12) << std::scientific << floor(b/delta_) << " " << ceil(t/delta_) << " " << floor(l) << " " << ceil(r) << std::endl;
//    std::cout << rowMin << " " << rowMax << " " << colMin << " " << colMax << std::endl;
//    std::cout << rowMinFull << " " << rowMaxFull << " " << colMinFull << " " << colMaxFull << std::endl;

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
      for(int ii = rowMin; ii <= rowMax; ii++){
        /*
         * ensure, that j <= i -> other elements do not contribute
         */
        UInt i = (UInt) ii;

        for(int jj = colMin; jj <= std::min(ii,colMax); jj++){

          UInt j = (UInt) jj;

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
          if((ii >= rowMinFull)&&(ii <= rowMaxFull)){
            if((jj >= colMinFull)&&(jj <= colMaxFull)){
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
              tmp = tmp * clipRectangleToElement(rect, i, j);
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
            tmp = tmp * clipRectangleToElement(rect, i, j);
          }

          sum += sign*tmp;
        }
      }

    } else {
      /*
       * std treatment
       */
      for(int ii = rowMin; ii <= rowMax; ii++){
        /*
         * ensure, that j <= i -> other elements do not contribute
         */
        UInt i = (UInt) ii;

        for(int jj = colMin; jj <= std::min(ii,colMax); jj++){
          UInt j = (UInt) jj;

          tmp = preisachWeights_[i][j];

          /*
           * Check for an inner element
           */
          if((ii >= rowMinFull)&&(ii <= rowMaxFull)){
            if((jj >= colMinFull)&&(jj <= colMaxFull)){
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
              tmp = tmp * clipRectangleToElement(rect, i, j);
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
            tmp = tmp * clipRectangleToElement(rect, i, j);
          }

          sum += tmp;
        }
      }
    }

    return sum;
  }

  Double VectorPreisachv10_ListApproach::getRectangleFromSwitchingList(std::list<ListEntryv10>& list,
                            std::list<ListEntryv10>::iterator startIt, std::list<ListEntryv10>::iterator curIt, std::list<ListEntryv10>::iterator endIt,
                            UInt idArea, Rectangle& rect, bool upperSplitSquare){
    /*
     * This functions is used to decompose a stair-case switching state into rectangular areas which are
     * easier to compute.
     * This functions will only be used for the switching lists
     *
     * upperSplitSquare:
     *  If true, the coordinates of the diagonally split upper square (area 0) will be returned)
     *  (cnt will be set to 0, to trigger the right calculations)
     *
     */
    /*
     * In older version we had to distinguish between upper square and upper triangle region when
     * calling this function. Depending on the region, different outer edges were set.
     * Actually, this was only needed for the upper triangle, as the resulting rectangular
     * bounds were used directly in the later evaluation. For the upper square, the resulting
     * bounds were clipped against the rotation area and thus possible extensions into other
     * regions would have been cut away.
     * In the new version, we clip always against the rotation state and thus we can use the
     * full Preisach plane as outer range. However, we now have to use the former TRIANGLE treatment in
     * each case, i.e. we have to cut the switching areas which go over the diagonal alpha = beta
     * down to triangle, so that the total shape of areas overlapping alpha = beta become squares.
     */

    Double alpha, beta, delta;
    alpha = -1.0;
    beta = -1.0;
    delta = 2.0;

    /*
    * Calculates the coordinates of a rectangle lying inside the area
    * (alpha, alpha+delta) x (beta, beta+delta)
    *
    * Input:
    *  list -> switching list
    *  startIt -> iterator pointing to the FIRST entry
    *  curIt -> iterator pointing to the CURRENT entry
    *  endIt -> iterator pointing to the LAST entry
    *  alpha,beta -> bottom left corner coordinates of outer area
    *  idArea -> index of area to be computed - 1 (i.e. area1 -> idArea = 0)
    *
    *  tRet,bRet -> y-coordinates of the top and bottom edge of the rectangle
    *  lRet,rRet -> x-coordinates of the left and right edge of the rectangle
    */

    /*
    * Splitting into splitting areas for e1 = max
    *
    *                         alpha
    *      __ __ __ __ __ __ __ |_  __ __ __ __ __ __
    *     |\    |          |       |               /
    *     |  \  |    A0    | A2    | A4          /
    *     |_ _ \|_ _ _ _ _ |       |           /___ e1
    *     |A1              |       |         /
    *     |_ _ _ _ _  _ _ _|_ _ __ |       /___ e3
    *     |A3                      |     /
    *     |_ _ _ _ _ _ _ _ _ _ _ _ |_ _/___ e5
    *     |A5                        /
    *     |                        /
    *     |                      /_|__________________ beta
    *     |                    /   |
    *     |                  /     |
    *     |                /       e4
    *     |              / |
    *     |            /   |
    *     |          /     e2
    *     |        /
    *     |     |/
    *     |    /|
    *     |_ /  |
    *          -e1
    *
    *     set e0 = -e1
    *
    *     negative areas (i = even)
    *       L = ei;             R = min(ei+2,1) > if ei+2 does not exist: ei+2 = 1
    *       B = max(ei+1,ei);   T = 1
    *             > limit extend of areas such, that the area below the diagonal is a triangle
    *             > if ei+1 does not exist: ei+1 = -1
    *
    *     positive areas (i = odd)
    *       L = -1;             R = min(ei+1,ei) > ensure that area below diagonal is triangular; if ei+1 does not exist: ei+1 = 1
    *       B = max(ei+2,-1);   T = ei
    *             > if ei+2 does not exist: ei+2 = -1
    *
    *     upper split area:
    *       L = -1; R = -e1; B = e1; T = 1
    *
    *
    * Splitting into splitting areas for e1 = min
    *
    *                         alpha
    *      __ __ __ __ __ __ __ |_  __ __ __ __ __ __
    *     |\    |     |           |                /
    *     |  \  | A1  |   A3      | A5           /
    *     |_ _ \|     |           |           _/___ -e1
    *     |A0   |     |           |          /
    *     |_ _ _|_ _  |           |        /___ e2
    *     |A2         |           |      /
    *     |_ _ ___ _ _|_ _  _ _ _ |   _/___ e4
    *     |A4                     |  /
    *     |                       |/
    *     |                      /|__________________ beta
    *     |                    /  |
    *     |                  /    |
    *     |                /      e5
    *     |              /
    *     |           |/
    *     |          /|
    *     |        /  |
    *     |     |/    e3
    *     |    /|
    *     |_ /  |
    *           e1
    *
    *     set e0 = -e1
    *
    *     negative areas (i = odd)
    *       L = ei;             R = min(ei+2,1) > if ei+2 does not exist: ei+2 = 1
    *       B = max(ei+1,ei);   T = 1
    *             > limit extend of areas such, that the area below the diagonal is a triangle
    *             (take A5 as example: as no e6 exists, A5 would reach down to the bottom of the Preisach plane (-1)
    *              by this, the area below the diagonal would be a trapezoid standing on its side
    *              by restricting to e5, we have only a triangle below the diagonal)
    *             > if ei+1 does not exist: ei+1 = -1
    *
    *     positive areas (i = even)
    *       L = -1;             R = min(ei+1,ei) > ensure that area below diagonal is triangular; if ei+1 does not exist: ei+1 = 1
    *       B = max(ei+2,-1);   T = ei
    *             > if ei+2 does not exist: ei+2 = -1
    *
    *     upper split area:
    *       L = -1; R = e1; B = -e1; T = 1
    *
    *
    * Note regarding different starting case (e1 = min or max):
    *   The rules for calculating positive and negative areas are independent of the type of e1.
    *   The only change is w.r.t. the indices. For e1 = min, positive areas have even indices and for
    *   e1 = max, they have odd indices. However, this is true for the type of ei, too. If list starts
    *   with minimum, all odd entries will be minima and if list starts with maximum, all odd entries will
    *   be maxima. Together with the calculation rules, we see that area Ai is positive, if ei is a maximum
    *   no matter what type e1 is. The calculation can thus be done without checking for indices or the
    *   type of e1. We only have to check, whether ei is a max or a min.
    *   Exception: Area0 is obtained from e1, too. To get this area, idArea has to be set to 0.
    *
    *
    * Note regarding initial values:
    *   When the classical model is used (version 2012), the lower triangle T_L gets not evaluated
    *   via switching lists as its switching state is always +1; to compensate for the +1 triangle, we
    *   should have an initial value of -1 for the upper triangle; this can easily be achieved by
    *   initializing the list with a max(or min) of value 0.
    *
    *
    */

    Double l,r,t,b;
    Double nextVal, nextnextVal;
    Double firstVal = startIt->getVal();
    Double curVal = curIt->getVal();
    bool firstMin = startIt->isMin();
    bool curMin = curIt->isMin();

    if(upperSplitSquare == true){

      l = -1;
      if(firstMin){
        r = firstVal;
        b = -firstVal;
      } else {
        r = -firstVal;
        b = firstVal;
      }
      t = 1;

      rect.setBounds(l,r,t,b);
      return 0.0; // value 0.0 not needed, but by this we can check if value was hit
    }

    if(idArea == 0){
      /*
       * use startIt regardless of the state of curIt!
       */
      if(firstMin){
        /*
         * List starts with minimum; get area 0 according to
         *
         *       L = -1;             R = e1
         *       B = max(e2,-1);     T = -e1
         *
         *       if e2 does not exist: e2 = -1
         *
         */
        if(startIt != endIt){
          startIt++;
          if(startIt->isMin() == firstMin){
            EXCEPTION("MinMaxList has to be alternating!")
          } else {
            nextVal = startIt->getVal();
          }
          startIt--;
        } else {
          nextVal = alpha;
        }

        l = beta;
        r = firstVal;
        b = nextVal; // will automatically be std::max(nextVal,alpha);
        t = -firstVal;

        /*
         * min-helper area is positive, so +1.0
         */
        rect.setBounds(l,r,t,b);
        return 1.0;

      } else {
        /*
         * List starts with maximum; get area 0 according to
         *
         *
         *       L = -e1;       R = min(e2,1)
         *       B = e1;        T = 1
         *
         *       if e2 does not exist: e2 = 1
         *
         */
        if(startIt != endIt){
          startIt++;
          if(startIt->isMin() == firstMin){
            EXCEPTION("MinMaxList has to be alternating!")
          } else {
            nextVal = startIt->getVal();
          }
          startIt--;
        } else {
          nextVal = beta+delta;
        }

        l = -firstVal;
        r = nextVal; // will automatically be std::min(nextVal,beta+delta);
        b = firstVal;
        t = alpha+delta;

        /*
         * max-helper area is negative, so -1.0
         */
        rect.setBounds(l,r,t,b);
        return -1.0;
      }
    } else {
      /*
       * area i with i > 0 -> use curIt
       */
      if(curMin){
        /*
         * entry is minimum; get area i according to
         *
         *       L = ei;             R = min(ei+2,1) > if ei+2 does not exist: ei+2 = 1
         *       B = max(ei+1,ei);   T = 1
         *             > limit extend of areas such, that the area below the diagonal is a triangle
         *             > if ei+1 does not exist: ei+1 = -1
         *
         */
        if(curIt != endIt){
          curIt++;

          if(curIt->isMin() == curMin){
            EXCEPTION("MinMaxList has to be alternating!")
          } else {
            nextVal = curIt->getVal();
          }

          if(curIt != endIt){
            curIt++;

            if(curIt->isMin() != curMin){
              EXCEPTION("MinMaxList has to be alternating!")
              /*
               * here we would expect the same entry type the the one of the current iteration
               * e.g. a minimum after a maximum which followed a minimum
               */
            } else {
              nextnextVal = curIt->getVal();
            }
            curIt--;
          } else {
            nextnextVal = beta+delta;
          }
          curIt--;
        } else {
          nextVal = alpha;
          nextnextVal = beta+delta;
        }

        l = curVal;
        r = nextnextVal; // will automatically be std::min(nextnextVal,beta+delta);
        b = std::max(curVal,nextVal);
        t = alpha+delta;

        /*
         * min area is negative, so -1.0
         */
        rect.setBounds(l,r,t,b);
        return -1.0;

      } else {
        /*
         * entry is maximum; get area i according to
         *
         *       L = -1;             R = min(ei+1,ei) > ensure that area below diagonal is triangular; if ei+1 does not exist: ei+1 = 1
         *       B = max(ei+2,-1);   T = ei
         *             > if ei+2 does not exist: ei+2 = -1
         *
         */

        if(curIt != endIt){
          curIt++;

          if(curIt->isMin() == curMin){
            EXCEPTION("MinMaxList has to be alternating!")
          } else {
            nextVal = curIt->getVal();
          }

          if(curIt != endIt){
            curIt++;

            if(curIt->isMin() != curMin){
              EXCEPTION("MinMaxList has to be alternating!")
              /*
               * here we would expect the same entry type the the one of the current iteration
               * e.g. a minimum after a maximum which followed a minimum
               */
            } else {
              nextnextVal = curIt->getVal();
            }
            curIt--;
          } else {
            nextnextVal = alpha;
          }
          curIt--;
        } else {
          nextVal = beta+delta;
          nextnextVal = alpha;
        }

        l = beta;
        r = std::min(curVal,nextVal);
        b = nextnextVal; // will automatically be std::max(nextnextVal,alpha);
        t = curVal;

        /*
         * max area is positive, so 1.0
         */
        rect.setBounds(l,r,t,b);
        return 1.0;
      }
    }
  }

  void VectorPreisachv10_ListApproach::Simplify_LocalSwitchingLists(std::list<RotListEntryv10>& usedList){
    /*
     * This function iterates over the globalRotation list and checks each entry in the corresponding
     * local switching list for overlap with the two possible rotation areas.
     * This search is performed until at least one overlap is found. All switching entries BEFORE the
     * one which lead to the first overlap are removed from the list. If at least 1 entry got removed,
     * startCnt will be set to 1 (as the entry of the switching list which corresponds to the helper area
     * (id = 0) no longer is in the list).
     * Leave out all rotListEntries for which the flag isUpdated is false.
     */

    std::list<ListEntryv10>::iterator swListIt;
    std::list<ListEntryv10>::iterator firstToKeep;
    std::list<ListEntryv10>::iterator swListEnd; // = --(globSwitchList_[idElem].end());
    std::list<ListEntryv10>::iterator swListStart;

    /*
     * iterators for global rotation list
     */
    std::list<RotListEntryv10>::iterator rotListIt;
    std::list<RotListEntryv10>::iterator rotListEnd = --(usedList.end());

    bool twoAreas,lastRotListEntry;

    Rectangle rotRect1 = Rectangle(0,0,0,0);
    Rectangle rotRect2 = Rectangle(0,0,0,0);
    Rectangle swRect = Rectangle(0,0,0,0);
    Rectangle overlapRect = Rectangle(0,0,0,0);

    UInt cnt;

    for(rotListIt = usedList.begin(); rotListIt != usedList.end(); rotListIt++){

      if(rotListIt == rotListEnd){
        lastRotListEntry = true;
      } else {
        lastRotListEntry = false;
      }
      /*
       * maybe the list was already simplified once or was inherited from a previous rotListState
       * in that case we have to make sure to not test for the wrong value
       */
      cnt = rotListIt->getStartCnt();

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
        continue;
      } else {
        /*
         * get rotation area(s)
         */

        swListStart = rotListIt->getListReference().begin();
        swListEnd = --(rotListIt->getListReference().end());

        twoAreas = getRectanglesFromRotEntry(rotListIt, rotRect1, rotRect2,lastRotListEntry);

        firstToKeep = rotListIt->getListReference().begin();
        swListEnd = --(rotListIt->getListReference().end());

        bool success1 = false;
        bool success2 = false;
        bool gotSuccess = false;

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
           */
          getRectangleFromSwitchingList(rotListIt->getListReference(),swListStart, swListIt, swListEnd,cnt, swRect);

          /*
           * clip resulting area against rotation area1
           */
          success1 = rotRect1.clipRectangles(swRect,overlapRect);

          if(twoAreas){
            /*
             * repeat the clipping steps for the second area
             */
            success2 = rotRect2.clipRectangles(swRect,overlapRect);
          }

          if((success1 == true)||(success2==true)){
            /*
             * the current switching area has an overlap with at least one of the two rotation areas
             * -> this entry has to be kept
             * -> end loop (as we suppose that all later entries have an overlap, too
             * This assumption is ok as only those entries are appended to the list which actually have
             * and effect, i.e. we only have to cut at the start of the list!
             */
            gotSuccess = true;
            break;
          }

          if(cnt > 0){

            /*
             * entry has no overlap, check next one
             * (start increasing only if cnt > 0 as list entry corrsponding to cnt == 0
             * has to be checked for area 0 and 1)
             */
            swListIt++;
            firstToKeep++;
          }
          cnt++;
        } // sw list

        if(firstToKeep != rotListIt->getListReference().begin()){
          /*
           * At least one entry has to be removed from the list
           */
          if(!gotSuccess){
            /*
             * NEW (30.6.2016)
             * check if we had an success at all!
             * for rotResistance < 1, it is possible, that there are rotation states
             * which have no overlap at all with switching areas (as xPar is too small).
             * In that case, this function will erase the whole switching list as there is
             * no overlap (ok so far). By setting the start cnt to value > 0 we do not
             * check the rotation state against the initially half split upper triangle
             * part (will only be checked for if start cnt = 0 -> see Evaluate_GlobalRotationList).
             * However, in the case of unsymmetric weights, we have to do just that!
             * -> if we had no success at all, keep starting counter at 0, so that the
             * later functions will overlap these rotation areas with the initial half
             * split upper square area.
             *
             */
            rotListIt->setStartCnt(0);
            //std::cout << "No overlap between switching list and rotation list. Will delete whole switching list! Overlap only with initially split area!" << std::endl;
            /*
             * furthermore, we are not allowed to delete the full list!
             * Reason: we need the first entry of the list, to determine the measure of the
             * initally split upper square part
             * -> leave first entry of switching list!
             */
            if(firstToKeep != ++(rotListIt->getListReference().begin())){
              rotListIt->getListReference().erase(++(rotListIt->getListReference().begin()),firstToKeep);
            }
          } else {
            rotListIt->setStartCnt(1);
            rotListIt->getListReference().erase(rotListIt->getListReference().begin(),firstToKeep);
          }
        } //else
        /*
         * whole list has to be kept
         * -> continue with next rotListElement
         */
      } // rot list has changed
    } // rot list
  }

  void VectorPreisachv10_ListApproach::Simplify_GlobalRotationList(std::list<RotListEntryv10>& usedList){
    /*
     * New merging rule (applicable for all versions)
     *
     * Idea: due to new insertion rule for switching lists (entries get clipped to bounding box
     *       created from rotlist entries), adjacent rotlist entries very seldom have the same switching
     *       list, even if the new entries xPar would override both switching lists
     *
     *       Example (a > b):
     *        rotentry1: bounding box1 left = -a, right = 1, bot = -1, top = a; inner list: max = (a+b)/2; min = b
     *        rotentry2: bounding box2 left = -b, right = 1, bot = -1, top = b; inner list: max = b/2
     *
     *        new entry xPar = b:
     *          rotentry1: bounding box1 left = -a, right = 1, bot = -1, top = a; inner list: max = (a+b)/2;
     *          rotentry2: bounding box2 left = -b, right = 1, bot = -1, top = b; inner list: max = b
     *          -> cannot be merged
     *
     *        new entry xPar = (a+b)/2:
     *          rotentry1: bounding box1 left = -a, right = 1, bot = -1, top = a; inner list: max = (a+b)/2
     *          rotentry2: bounding box2 left = -b, right = 1, bot = -1, top = b; inner list: max = b -> because of clipping
     *                      (without clipping: inner list: max = (a+b)/2 -> could be merged)
     *          -> cannot be merged due to old merging rule in combination with new clipping rule
     *
     * New merging rule:
     *        entry_i_last shall be the last entry of the switching list switchinglist_i belonging to rotlist i
     *        entry_i+1_first shall be the first entry of the switching list switchinglist_i+1 belonging to rotlist i+1
     *
     *        if length(switchinglist_i+1) == 1    <- always the case if switchinglist_i+1 was completely overwritten
     *          if type(entry_i_last) == type(entry_i+1_first):
     *            if(entry_i_last == entry_i+1_first)   <- new list starts where old list ended
     *              OR
     *              if(type(entry_i_last == max)
     *                if(entry_i_last > entry_i+1_first) AND entry_i+1_first == top value of bounding box
     *                  -> merge
     *              else
     *                if(entry_i_last < entry_i+1_first) AND entry_i+1_first == left value of bounding box
     *                  -> merge
     *
     * Idea behind new rule:
     *       - merging in old version was more or less only possible if list of entry i+1 was completely overwritten;
     *        otherwise it is very seldom the case that both switching list can match
     *       - for the old rule, merging could only be done if the lists were the same for all entries;
     *        however, the first list can contain entries which do not affect the later rotation state
     *        (at least due to the new shape of the rotation entries coming from version 10); these entries
     *        will never be inserted into the switching list of the next following rotation state due to the
     *        new clipping and entries which was inherited will be sorted out in the simplify function for the
     *        switching list -> merging becomes nearly impossible for version 10
     *       - the new merging rule solves this problem as we already know that we can only merge if
     *        the switching list of entry i+1 is overwritten and that this will be the case no matter
     *        what entries the list of entry i has, as long as the last entry is a smaller or equal min
     *        /a larger of equal max
     *
     *       Example (a > b):
     *        rotentry1: bounding box1 left = -a, right = 1, bot = -1, top = a; inner list: max = (a+b)/2; min = b
     *        rotentry2: bounding box2 left = -b, right = 1, bot = -1, top = b; inner list: max = b/2
     *
     *        new entry xPar = b:
     *          rotentry1: bounding box1 left = -a, right = 1, bot = -1, top = a; inner list: max = (a+b)/2
     *          rotentry2: bounding box2 left = -b, right = 1, bot = -1, top = b; inner list: max = b
     *          -> can be merged now, as (a+b)/2 > b and b = top
     *
     *        new entry xPar = (a+b)/2:
     *          rotentry1: bounding box1 left = -a, right = 1, bot = -1, top = a; inner list: max = (a+b)/2
     *          rotentry2: bounding box2 left = -b, right = 1, bot = -1, top = b; inner list: max = b -> because of clipping
     *                      (without clipping: inner list: max = (a+b)/2 -> could be merged)
     *          -> can be merged now, as (a+b)/2 > b and b = top
     *
     */


    /*
     * iterate over rotation list and check if two adjacent entries have the same rotation state and the
     * same list of switching entries
     * This case can happen e.g. if an input overrides multiple older rotation areas at once and the
     * resulting value of xPar is large enough to overwrite the contained switching lists
     */

    if(usedList.size() < 2){
      /*
       * list has not enough entries to be merged together
       */
      return;
    }

    std::list<RotListEntryv10>::iterator listIt;
    std::list<RotListEntryv10>::iterator nextListIt;
    std::list<RotListEntryv10>::iterator listEnd = --(usedList.end());

    for(listIt = usedList.begin(); listIt != usedList.end(); ){

      /*
       *  get next entry in list (if any)
       *  if not -> nothing more to do here
       */
      if(listIt != listEnd){
        nextListIt = listIt;
        nextListIt++;

        /*
         * check if entries can be merged together
         * (v2 = new merging rule; without v2 = classical merging)
         */
        if(listIt->canBeMergedWith_v2(*nextListIt,classical_)){
          /*
           * get lower bound of second element and set lower bound of first element to that value
           * (-> i.e. we extend the area of the first element)
           */
          Double low = nextListIt->getLowerVal();
          listIt->setLowerVal(low);
          /*
           * now we can delete the second entry
           */
          usedList.erase(nextListIt);

          /*
           * reset iterator to last entry
           * (normally the iterators should remain valid, but I do not trust them ...
           */
          listEnd = --(usedList.end());

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

  void VectorPreisachv10_ListApproach::mapRectangleToHelperMatrix(Matrix<Double>& helper, Rectangle rect, Double factor, bool skipUpperDiagonal, bool isRotState){
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
     * new flag: isRotState
     * -> if the rotation state is mapped to the helper matrix, we do not want to scale its value by 0.5 if its on the diagonal alpha = beta
     *    if we would do so, we would get the factor 0.5 twice in the final images which are created from the matrix (1. from switching state,
     *    2. from rotation state)
     * -> if isRotState == True, scale full entries on alpha = beta times 2 and for partially filled elements, skip
     *    the cutting of triangle parts below the diagonal
     *
     */

    Double l,r,t,b;
    rect.getBounds(l,r,t,b);

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
    int rowMin =  floor(b/delta) + floor(1.0/delta);
    int rowMax = ceil(t/delta) + floor(1.0/delta)-1;

    int colMin =  floor(l/delta) + floor(1.0/delta);
    int colMax =  ceil(r/delta) + floor(1.0/delta)-1;

    /*
     * NEW: use integers instead of unsigned integer
     * REASON: colMaxFull became negative (which simply would indicate no fully overlapped elemets), but
     * unsigning took it to a very large positive value and thus -> wrong result
     */

    /*
     * 2. get indices of completely overlapped elements
     */
    int rowMinFull =  ceil(b/delta) + floor(1.0/delta);
    int rowMaxFull = floor(t/delta) + floor(1.0/delta)-1;

    int colMinFull =  ceil(l/delta) + floor(1.0/delta);
    int colMaxFull = floor(r/delta) + floor(1.0/delta)-1;


    if(skipUpperDiagonal == true){
      /*
       * special treatment:
       * 1. skip entries on alpha = -beta
       * 2. consider the right signs, i.e. +1 for entries below alpha = -beta
       *    and -1 for entries above!
       */
      Double sign;
      for(int ii = rowMin; ii <= rowMax; ii++){
        /*
         * ensure, that j <= i -> other elements do not contribute
         */
        UInt i = (UInt) ii;

        for(int jj = colMin; jj <= std::min(ii,colMax); jj++){

          UInt j = (UInt) jj;

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
          if((ii >= rowMinFull)&&(ii <= rowMaxFull)){
            if((jj >= colMinFull)&&(jj <= colMaxFull)){
              /*
               * check for diagonal element on alpha = beta
               * (should not happen here)
               */
              if((i == j)&&(isRotState==false)){
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
              helper[i][j] += sign*clipRectangleToElement(rect, i, j, delta,isRotState)/(delta*delta);
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
            helper[i][j] += sign*clipRectangleToElement(rect, i, j, delta,isRotState)/(delta*delta);
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
      for(int ii = rowMin; ii <= rowMax; ii++){
        /*
         * ensure, that j <= i -> other elements do not contribute
         */
        UInt i = (UInt) ii;

        for(int jj = colMin; jj <= std::min(ii,colMax); jj++){
          /*
           * Check for an inner element
           */
          UInt j = (UInt) jj;

          if((ii >= rowMinFull)&&(ii <= rowMaxFull)){
            if((jj >= colMinFull)&&(jj <= colMaxFull)){
              /*
               * check for diagonal element
               */
              if((i == j)&&(isRotState==false)){
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
              helper[i][j] += factor*clipRectangleToElement(rect, i, j, delta,isRotState)/(delta*delta);
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
            helper[i][j] += factor*clipRectangleToElement(rect, i, j, delta,isRotState)/(delta*delta);
          }
        }
      }
    }
  }


  void VectorPreisachv10_ListApproach::switchingStateToBmp(UInt numPixel, std::string filename, UInt idElem, bool overLayWithRotState)
  {
    /*
     * NEW: rotation state is evaluated along with the switching states if overLayWithRotState is true
     * in the old versions, the rotation state was evaluated separately although we have to iterate over the
     * rotation list to evaluate the switching lists
     */

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
     * create matrices to store the rotation state
     */
    Matrix<Double> rotX, rotY, rotZ;
    if(overLayWithRotState){
      rotX = Matrix<Double>(numPixel,numPixel);
      rotY = Matrix<Double>(numPixel,numPixel);
      rotX.Init();
      rotY.Init();
      if(dim_ == 3){
        rotZ = Matrix<Double>(numPixel,numPixel);
        rotZ.Init();
      }
    }

    /*
     * Fill matrix / matrices
     * 1. if(classical_)
     *  i. clip lowerTriangle to helperMatrix
     *
     * 2. iterate over rotation list
     *  i. store rotation direction if needed
     *  ii. iterate over switching list
     *    a. determine overlapping areas and clip them against helperMatrix
     *    b. for unsymmetric weights include triagonal areas
     */

    /*
     * Part 1
     */
    if(classical_){
      for(UInt i = 0; i < numPixel/2; i++){
        for(UInt j = 0; j <= i; j++){
          if(i == j) helperMatrix[i][j] = 0.5;
          else helperMatrix[i][j] = 1.0;
        }
      }

      if(overLayWithRotState){
        for(UInt i = 0; i < numPixel/2; i++){
          for(UInt j = 0; j <= i; j++){
            rotX[i][j] = lastEu_[idElem][0];
          }
        }
        for(UInt i = 0; i < numPixel/2; i++){
          for(UInt j = 0; j <= i; j++){
            rotY[i][j] = lastEu_[idElem][1];
          }
        }
        if(dim_ == 3){
          for(UInt i = 0; i < numPixel/2; i++){
            for(UInt j = 0; j <= i; j++){
              rotZ[i][j] = lastEu_[idElem][2];
            }
          }
        }
      }
    }

    /*
     * Part 2
     */

    /*
     * iterators for local switching list
     */
    std::list<ListEntryv10>::iterator swListIt;
    std::list<ListEntryv10>::iterator swListStart;
    std::list<ListEntryv10>::iterator swListEnd;

    /*
     * iterators for global rotation list
     */
    std::list<RotListEntryv10>::iterator rotListIt;
    std::list<RotListEntryv10>::iterator rotListEnd = --(globRotList_[idElem].end());

    bool area0 = true;
    bool twoAreas = true;
    bool lastRotListEntryv10;
    Double upperBound;

    Rectangle rotRect1 = Rectangle(0,0,0,0);
    Rectangle rotRect2 = Rectangle(0,0,0,0);
    Rectangle swRect = Rectangle(0,0,0,0);
    Rectangle overlapRect = Rectangle(0,0,0,0);

    Vector<Double> curRotState;
    Double factor;
    UInt cnt = 0;

    for(rotListIt = globRotList_[idElem].begin(); rotListIt != globRotList_[idElem].end(); rotListIt++){

      if(rotListIt == rotListEnd){
        lastRotListEntryv10 = true;
      } else {
        lastRotListEntryv10 = false;
      }

      curRotState = rotListIt->getVecReference();
      upperBound = rotListIt->getVal(); //exi

      swListStart = rotListIt->getListReference().begin();
      swListEnd = --(rotListIt->getListReference().end());

      twoAreas = getRectanglesFromRotEntry(rotListIt, rotRect1, rotRect2,lastRotListEntryv10);

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

        if(classical_){
          /*
           * area 0 consists only of one square region
           */
          Rectangle area0 = Rectangle(-1.0,-upperBound,1.0,upperBound);

          /*
           * map rectangle to HelperMatrix and set flag skipUpperDiagonal to true
           * -> by setting the flag to true, the function will automatically assume that the area shall
           * be filled with +1 up to the splitting diagonal and -1 above it
           */
          mapRectangleToHelperMatrix(helperMatrix,area0,0,true);
        } else {
          /*
           * in the revised version, area 0 has the same L-kind shape as the other rotation states;
           * however, to get it split properly, we divide this area into three parts
           * a) square part which is split into +1 and -1 -> same as in classical model
           * b) left, vertical rectangle -> completely set to +1
           * c) upper, horizontal rectangle -> completely set to -1
           *
           */
          Rectangle area0_square = Rectangle(-1.0,-upperBound,1.0,upperBound);
          mapRectangleToHelperMatrix(helperMatrix,area0_square,0,true);

          Rectangle area0_left = Rectangle(-1.0,-upperBound,upperBound,-1.0);
          mapRectangleToHelperMatrix(helperMatrix,area0_left,1.0,false);

          Rectangle area0_top = Rectangle(-upperBound,1.0,1.0,upperBound);
          mapRectangleToHelperMatrix(helperMatrix,area0_top,-1.0,false);
        }

        area0 = false;

        /*
         * Note: we do not have to write to the matrix for the rotation states as we have no rotation state here
         */
      }

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
       *  anymore (as this is the special helper area). So we startCnt in rotListEntryv10 to mark that we
       *  deleted the first entry which originally corresponded to area 0. As all areas with id > 0 are
       *  calculated by the same scheme, the absolute value of cnt is not of interest and thus it is
       *  either 0 or 1.
       */
      cnt = rotListIt->getStartCnt();

      bool success = false;

      for(swListIt = rotListIt->getListReference().begin(); swListIt != rotListIt->getListReference().end(); ){

        /*
         * get rectangular bounds corresponding to entry of switching list (compare to Evaluate_GlobalSwitchingList)
         */
        factor = getRectangleFromSwitchingList(rotListIt->getListReference(),swListStart, swListIt, swListEnd,cnt, swRect);

        /*
         * clip resulting area against rotation area1
         */
        success = rotRect1.clipRectangles(swRect,overlapRect);

        if(success){
          /*
           *  clip overlap to HelperMatrix (factor holds the value +1 or -1 and indicates how the matrix shall be filled)
           */
          mapRectangleToHelperMatrix(helperMatrix,overlapRect,factor);
        }

        if(twoAreas){
          /*
           * repeat the clipping steps for the second area
           */
          success = rotRect2.clipRectangles(swRect,overlapRect);

          if(success){
            /*
             *  clip overlap to HelperMatrix (factor holds the value +1 or -1 and indicates how the matrix shall be filled)
             */
            mapRectangleToHelperMatrix(helperMatrix,overlapRect,factor);
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

          /*
           * only difference to evaluation algorithm: do this also, if list was already wiped
           */
          if((isSymmetric_ == false)){
            /*
             * set flag upperSplitSquare to true -> get area 0 instead of area 1
             */
            factor = getRectangleFromSwitchingList(rotListIt->getListReference(),swListStart, swListIt, swListEnd,cnt, swRect,true);

            if(factor != 0){
              EXCEPTION("Something got wrong here! Check function getRectangleBounds!");
            }

            /*
             * clip resulting area against rotation area1
             */
            success = rotRect1.clipRectangles(swRect,overlapRect);

            if(success){
              /*
               *  clip overlap to HelperMatrix (factor holds the value +1 or -1 and indicates how the matrix shall be filled)
               * -> NOTE: here we set the flag skipUpperDiagonal to true!
               */
              mapRectangleToHelperMatrix(helperMatrix,overlapRect,factor,true);
            }

            if(twoAreas){
              /*
               * repeat the clipping steps for the second area
               */
              success = rotRect2.clipRectangles(swRect,overlapRect);

              if(success){
                /*
                 *  clip overlap to HelperMatrix (factor holds the value +1 or -1 and indicates how the matrix shall be filled)
                 * -> NOTE: here we set the flag skipUpperDiagonal to true!
                 */
                mapRectangleToHelperMatrix(helperMatrix,overlapRect,factor,true);
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

      if(overLayWithRotState){
        /*
         * write information about rotation state into matrices
         */
        Double currentEntry_x = curRotState[0];

        mapRectangleToHelperMatrix(rotX,rotRect1,currentEntry_x,false,true);

        if(twoAreas){
          mapRectangleToHelperMatrix(rotX,rotRect2,currentEntry_x,false,true);
        }

        Double currentEntry_y = curRotState[1];

        mapRectangleToHelperMatrix(rotY,rotRect1,currentEntry_y,false,true);

        if(twoAreas){
          mapRectangleToHelperMatrix(rotY,rotRect2,currentEntry_y,false,true);
        }

        if(dim_ == 3){
          Double currentEntry_z = curRotState[2];

          mapRectangleToHelperMatrix(rotZ,rotRect1,currentEntry_z,false,true);

          if(twoAreas){
            mapRectangleToHelperMatrix(rotZ,rotRect2,currentEntry_z,false,true);
          }
        }
      }
    } // rot list

    /*
     * now call output function of matrix
     */
    UInt upscaling = 1;

    if(overLayWithRotState == true){
      UInt version = 2;

      if(version == 1){

        for(UInt comp = 0; comp < dim_; comp++){

          std::stringstream stream;
          std::string filename_new;
          if(comp == 0){
            stream << "x-" << filename;
            filename_new = stream.str();
            helperMatrix.matrix2Bmp(upscaling,filename_new,&rotX);
          } else if(comp == 1){
            stream << "y-" << filename;
            filename_new = stream.str();
            helperMatrix.matrix2Bmp(upscaling,filename_new,&rotY);
          } else {
            stream << "z-" << filename;
            filename_new = stream.str();
            helperMatrix.matrix2Bmp(upscaling,filename_new,&rotZ);
          }
        }

      } else if(version == 2) {
       /*
        * New way of outputting matrix:
        *   encode rotation state in color: 0° = red; 120° = blue; 240° = green; angles in between colored as a mix
        *   encode switching state as sign and amplitude: negative switching -> angle + 180°; absvalue < 1 -> scale final colorcombination
        *
        *   -> only for 2D rotstates, as the z component is not considered
        */

        std::stringstream stream;
        stream << "xy-" << filename;
        std::string filename_new = stream.str();

        helperMatrix.matrix2Bmp_v2(upscaling,filename_new,&rotX,&rotY);
//
//        std::cout << "switching: \n" << helperMatrix.ToString() << std::endl;
//        std::cout << "rotx: \n " << rotX.ToString() << std::endl;
//        std::cout << "roty: \n" << rotY.ToString() << std::endl;
      }
    } else {
      helperMatrix.matrix2Bmp(upscaling,filename,NULL);
    }
  }

}//end namespace
