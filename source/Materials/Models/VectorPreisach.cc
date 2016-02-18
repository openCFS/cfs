#include "VectorPreisach.hh"

#include <iostream>
#include <fstream>

/* Based on vector hyteresis model by Sutor et al.
 * "An efficient vector Preisach hysteresis model based on a novel rotational operator"
 * doi: http://dx.doi.org/10.1063/1.3672069
 */

namespace CoupledField
{ 

  VectorPreisach::VectorPreisach(Integer numElem, Double xSat, Double ySat,
		       Matrix<Double>& preisachWeight, Double rotationalResistance, UInt dim,  bool isVirgin)
            : Hysteresis(numElem)
          //: Preisach(numElem, xSat, ySat, preisachWeight, isVirgin)
  {

    eps_ = 1e-10;

    if (xSat > 0 ) {
      Xsaturated_  = xSat;
    }
    else {
      Xsaturated_  = 1.0;
    }

    YSaturated_  = ySat;
    isVirgin_    = isVirgin;

    preisachWeights_ = preisachWeight;

    // preisachWeights_ and rotationalWeights_ have to have the same size
    UInt M = preisachWeights_.GetNumRows();
    UInt N = preisachWeights_.GetNumCols();

    if(M != N){
      EXCEPTION("Matrix preisachWeight has dim " << M << " x " << N << " and thus is not symmetric!");
    }

    preisachSum_ = new Vector<Double>[numElem];

    switchingStates_ = new Matrix<Integer>[numElem];
    switchingStatesUpdated_ = new bool[numElem];

    dim_ = dim;

    rotationalWeightsX_ = new Matrix<Double>[numElem];
    rotationalWeightsY_ = new Matrix<Double>[numElem];
    if(dim_ == 3){
      rotationalWeightsZ_ = new Matrix<Double>[numElem];
    }
    rotationalWeightsUpdated_ = new bool[numElem];

    for(UInt k = 0; k < (UInt) numElem; k++){

      // initialize switchingStates_
      switchingStates_[k] = Matrix<Integer>(M,N);

      /* current version of Preisach plane consists of small squares
       * i.e. we have no diagonal line but rather a stair-case
       *
       * Initial state
       * N = M = 4            N = M = 3
       *  _ _ _ _            _ _ _
       * |0 - - -|          |0 - -|
       * |+ 0 -  |          |+ 0  |
       * |+ +    |          |+_ _ |
       * |+_ _ _ |
       *
       */

      // initialize switchingStates along staircaseline alpha = beta
      UInt half;
      if((M%2) == 0){
        half = M/2;
      } else {
        half = M/2; // round down
      }

      // lower triangle 0,0 - 0,half - half,half
      for(UInt i = 0; i < half; i++)
      {
        for(UInt j = 0; j <= i; j++ )
        {
          switchingStates_[k][i][j] = +1;
        }
      }

      // square 0,half - M-half,half - 0,M - M-half,M
      for(UInt i = half; i < M; i++)
      {
        for(UInt j = 0; j < M-half; j++)
        {
          if(i < j) switchingStates_[k][i][j] = +1;
          else if(j > i) switchingStates_[k][i][j] = -1;
          else switchingStates_[k][i][j] = 0;
        }
      }

      // triangle M-half,M-half - M,M-half - M,M
      for(UInt i = M-half; i < M; i++)
      {
        for(UInt j = M-half; j <= i; j++)
        {
          switchingStates_[k][i][j] = -1;
        }
      }

      switchingStatesUpdated_[k] = false;

      // initialize all rotationa vector components with 0
      rotationalWeightsX_[k] = Matrix<Double>(M,N);
      rotationalWeightsX_[k].Init();

      rotationalWeightsY_[k] = Matrix<Double>(M,N);
      rotationalWeightsY_[k].Init();

      if(dim_ == 3){
        rotationalWeightsZ_[k] = Matrix<Double>(M,N);
        rotationalWeightsZ_[k].Init();
      }
      rotationalWeightsUpdated_[k] = false;

      preisachSum_[k] = Vector<Double>(dim);
      preisachSum_[k].Init(0.0);

    }

    rotationalResistance_ = rotationalResistance;
   }

  VectorPreisach::~VectorPreisach()
  {
  }

  void VectorPreisach::updateRotationalWeights(Vector<Double>& Xin, Integer idx){

    //std::cout << "Update Rot Weights for idx: " << idx << std::endl;
    // determine computational threshold
    Double X_thres = std::pow((Xin.NormL2()/Xsaturated_),rotationalResistance_);

    //std::cout << "X_thres: " << X_thres << std::endl;

    // get direction of Xin
    Vector<Double> e_u = Vector<Double>(Xin.GetSize());
    //std::cout << "Xin.NormL2(): " << Xin.NormL2() << std::endl;
    if(Xin.NormL2() > eps_){
      e_u = Xin/Xin.NormL2();
    } else {
      e_u.Init(0.0);
      // for small values of input, keep old state
      //rotationalWeightsUpdated_[idx] = true;
      //return;
    }

    //std::cout << "e_u: " << e_u << std::endl;

    // go over rotationalWeights and update entries
    UInt M = preisachWeights_.GetNumRows();
    Double delta = 2.0 / ( (Double) M );
    // take value of alpha and beta in the cell centers
    Double alpha = -1 + delta/2.0;

    for (UInt i = 0; i < M; i++){

      Double beta = -1 + delta/2.0;
      for (UInt j = 0; j <= i; j++){

        if((X_thres > alpha) || (-X_thres < beta)){
          rotationalWeightsX_[idx][i][j] = e_u[0];
          rotationalWeightsY_[idx][i][j] = e_u[1];
          if(dim_ == 3){
            rotationalWeightsZ_[idx][i][j] = e_u[2];
          }
        } else {
          // keep old state
        }

        beta += delta;
      }
      alpha += delta;
    }

    rotationalWeightsUpdated_[idx] = true;

  }

  void VectorPreisach::updateSwitchingOperators(Vector<Double>& Xin, Integer idx){
    //std::cout << "Update switches for idx: " << idx << std::endl;

    if(rotationalWeightsUpdated_[idx] == false){
      // update the rotationalWeights first
      updateRotationalWeights(Xin, idx);
    }

    // go over switchingStates_ and update entries
    UInt M = preisachWeights_.GetNumRows();
    Double delta = 2.0 / ( (Double) M );
    // take value of alpha and beta in the cell centers
    Double alpha = -1 + delta/2.0;

    Double Xpar;

    for (UInt i = 0; i < M; i++){

      Double beta = -1 + delta/2.0;
      for (UInt j = 0; j <= i; j++){

        // at each combination of alpha and beta determine the component of X pointing along the current rotational vector
        Xpar = 0;
        Xpar += Xin[0]*rotationalWeightsX_[idx][i][j];
        Xpar += Xin[1]*rotationalWeightsY_[idx][i][j];
        if(dim_ == 3){
          Xpar += Xin[2]*rotationalWeightsZ_[idx][i][j];
        }
        //normalize Xpar to Xsaturatred
        Xpar /= Xsaturated_;

        if(Xpar > alpha){
          switchingStates_[idx][i][j] = 1;
        } else if(Xpar < beta){
          switchingStates_[idx][i][j] = -1;
        } else {
          // keep old state
        }
        //std::cout << "alpha / beta = " << alpha << " / " << beta << std::endl;
        beta += delta;
      }
      alpha += delta;
    }

    switchingStatesUpdated_[idx] = true;

  }

  Vector<Double> VectorPreisach::computeValue_vec(Vector<Double>&Xin, Integer idxElem)
  {
    //std::cout << "Compute value for idx: " << idxElem << std::endl;

    if(switchingStatesUpdated_[idxElem] == false)
    {
      updateSwitchingOperators(Xin,idxElem);
    }

    Vector<Double> Yout = Vector<Double>(dim_);
    Yout.Init(0.0);

    Double weight = 0;
    UInt M = preisachWeights_.GetNumRows();
    Double delta = 2.0 / ( (Double) M );
    Double area = delta*delta; // for the beginning treat all areas as small squares
    //Double sum_weights = 0;

    // here we have to iterate over the whole Preisach plane and evaluete the integral for You
    for (UInt i = 0; i < M; i++){

      for (UInt j = 0; j <= i; j++){
         //sum_weights += preisachWeights_[i][j];
         weight = preisachWeights_[i][j]*switchingStates_[idxElem][i][j];
         if(i == j){
           // entries on diagonal have only triangular areas, so scale them with an additional 0.5
           weight *= 0.5;
           //sum_weights -= 0.5*preisachWeights_[i][j];
         }
         Yout[0] += weight*rotationalWeightsX_[idxElem][i][j];
         Yout[1] += weight*rotationalWeightsY_[idxElem][i][j];
         if(dim_ == 3){
           Yout[2] += weight*rotationalWeightsZ_[idxElem][i][j];
         }
      }
    }
    //sum_weights *= area;
    //std::cout << "Sum of all weights: " << sum_weights << std::endl;

    for (UInt k = 0; k < Yout.GetSize();k++)
      preisachSum_[idxElem][k] = area*YSaturated_*Yout[k];

    //std::cout << "area: " << area << std::endl;
    // reset update function
    rotationalWeightsUpdated_[idxElem] = false;
    switchingStatesUpdated_[idxElem] = false;

    return preisachSum_[idxElem];
  }

  Vector<Double> VectorPreisach::getValue_vec( Integer idx )
  {
    return ( preisachSum_[idx] );
  }

}
