#include <cmath>

#include "Domain/elem.hh"
#include "Domain/entityList.hh"
#include "Elements/basefe.hh"
#include "Forms/fluidMechInt.hh"
#include "General/environment.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/tools.hh"
#include "fluidMechMassInt.hh"

namespace CoupledField
{

/***********************************************************
 * constructors
 **********************************************************/
FluidMechPlaneMassNewtonInt_UV::FluidMechPlaneMassNewtonInt_UV(Double density, \
  Double kinematicViscosity, Matrix<Double> & stabilParams, \
  bool movingMesh, std::string stabilType)
: FluidMechInt(density,kinematicViscosity,movingMesh, stabilType), \
  stabilParams_(stabilParams)
{
  isSolDependent_ = true;
  isaxi_=false;
  name_ = "FluidMechPlaneMassNewtonInt_UV";
  baseType_ = MASS;
}
FluidMechPlaneMassNewtonInt_UQ::FluidMechPlaneMassNewtonInt_UQ(Double density, \
  Double kinematicViscosity, Matrix<Double> & stabilParams, \
  bool movingMesh, std::string stabilType)
: FluidMechInt(density,kinematicViscosity,movingMesh, stabilType), \
  stabilParams_(stabilParams)
{
  isSolDependent_ = true;
  isaxi_=false;
  name_ = "FluidMechPlaneMassNewtonInt_UQ";
  baseType_ = MASS;
}

/***********************************************************
 * destructors
 **********************************************************/
FluidMechPlaneMassNewtonInt_UV::~FluidMechPlaneMassNewtonInt_UV()
{
}
FluidMechPlaneMassNewtonInt_UQ::~FluidMechPlaneMassNewtonInt_UQ()
{
}

/***********************************************************
 * calc matrix
 **********************************************************/
void FluidMechPlaneMassNewtonInt_UV::CalcElementMatrix( Matrix<Double>& elemMat, \
    EntityIterator& ent1, \
    EntityIterator& ent2 )
{
  // Extract pointer to reference element and get coordinates
  ExtractElemInfo( ent1 );
  UInt elemNumber = ent1.GetElem()->elemNum;

  sol_->GetElemSolutionAsMatrix( elemResult_, ent1 );

  Matrix<Double> gridElemResult;
  if (movingMesh_) {
    gridSol_->GetElemSolutionAsMatrix( gridElemResult, ent1);
    elemResult_ -= gridElemResult;
  }

  bool computeTaus = false;
  if (stabilParams_[elemNumber-1][3] < -0.001)
  {
    computeTaus = true;
  }

  ptelem->SetAnsatzFct( ansatzFct1_ );
  const UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
  const UInt nrIntPts = ptelem->GetNumIntPoints();

  const Vector<Double>& intWeights = ptelem->GetIntWeights();

  const UInt N = spaceDim_;  // DOFs per Node (Ux,Uy)
  Double jacDet, jacDet_intWeight;
  Double lambda_k;
  Vector<Double> Vx, Vy;
  Double VxAtIP, VyAtIP, VL2, VL2AtIp, VMax;
  Double VxAtIP_multAux, VyAtIP_multAux;
  // need to be > 1e-15 and < 1e3
  Double& tau_mp  = stabilParams_[elemNumber-1][0];
  Double& tau_mu = stabilParams_[elemNumber-1][1];
  Double& tau_c  = stabilParams_[elemNumber-1][2];
  // Double tau_mp, tau_mu, tau_c;
  Double A_elem, h_k;
  if (!allocMem)
  {
    doMemAlloc(nrFncs);
  }

  Vx.Replace(nrFncs, elemResult_[0], false);
  Vy.Replace(nrFncs, elemResult_[1], false);

  if (computeTaus) {
    A_elem = ptelem->CalcVolume(ptCoord_, isaxi_);
    //A_elem /= 4.0; //wg quadratische Element

    VL2 = sqrt( (Vx.NormL2()*Vx.NormL2()) + (Vy.NormL2()*Vy.NormL2()) );
    h_k = sqrt(A_elem);

    VMax=abs(Vx[0]);
    if (abs(Vy[0])>VMax)
      VMax=abs(Vy[0]);

    for (UInt i=1; i<nrFncs; i++) {
      if (abs(Vx[i]) > abs(Vy[i]) && abs(Vx[i]) > VMax )
        VMax=abs(Vx[i]);
      else if (abs(Vx[i]) < abs(Vy[i]) && abs(Vy[i]) > VMax )
        VMax=abs(Vy[i]);
    }

    VL2AtIp = -1.0;
    ComputeTaus(tau_mp,tau_mu,tau_c,VL2,VL2AtIp,VMax,h_k,kinematicViscosity_,lambda_k,nrFncs);
    stabilParams_[elemNumber-1][3]=1.0;
  }
  // set matrix to desired size and set all elements to zero
  locElemMat_UV_.Init();

  //#pragma omp parallel for private(auxJ_,VxAtIP,VyAtIP,VzAtIP,jacDet,xiDxDy_,xiDxxDyyDxy,xiDxx,xiDyy,xiDzz,xiDx,xiDy,xiDz,xi_,A_1_)
  for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
  {
    jacDet = 0.0;

    ptelem->GetShFncAtIp(xi_, actIntPt, ent1.GetElem());
    ptelem->GetGlobDerivShFncAtIp(xiDxDy_, actIntPt, ptCoord_, jacDet, ent1.GetElem());

    jacDet_intWeight = intWeights[actIntPt-1] * jacDet;

    for (UInt i=0; i<nrFncs; ++i)
    {
      DxXi_[i] = xiDxDy_[i][0];
      DyXi_[i] = xiDxDy_[i][1];
    }

    Vx.Inner(xi_, VxAtIP);
    VxAtIP_multAux = VxAtIP * tau_mu * jacDet_intWeight;
    Vy.Inner(xi_, VyAtIP);
    VyAtIP_multAux = VyAtIP * tau_mu * jacDet_intWeight;

    // Mass Integrator
    auxI_ = xi_ * jacDet_intWeight;
    DyadicMult(A_1_, xi_, auxI_);

    auxI_ = xi_ * VxAtIP_multAux;
    DyadicMult(A_tmp_, DxXi_, auxI_);
    A_1_ += A_tmp_;

    auxI_ = xi_ * VyAtIP_multAux;
    DyadicMult(A_tmp_, DyXi_, auxI_);
    A_1_ += A_tmp_;

    // assemble element matrix
    //#pragma omp critical
    {
      AddSubMatrix(locElemMat_UV_, A_1_, 0, 0);
      AddSubMatrix(locElemMat_UV_, A_1_, nrFncs, nrFncs);
    }
  }
  ResortElementMatrix(elemMat, locElemMat_UV_, nrFncs, N);
}

void FluidMechPlaneMassNewtonInt_UQ::CalcElementMatrix( Matrix<Double>& elemMat, \
    EntityIterator& ent1, \
    EntityIterator& ent2 )
{
  // Extract pointer to reference element and get coordinates
  ExtractElemInfo( ent1 );
  UInt elemNumber = ent1.GetElem()->elemNum;

  sol_->GetElemSolutionAsMatrix( elemResult_, ent1 );

  Matrix<Double> gridElemResult;
  if (movingMesh_) {
    gridSol_->GetElemSolutionAsMatrix( gridElemResult, ent1);
    elemResult_ -= gridElemResult;
  }

  ptelem->SetAnsatzFct( ansatzFct1_ );
  const UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
  const UInt nrIntPts = ptelem->GetNumIntPoints();

  const Vector<Double>& intWeights = ptelem->GetIntWeights();

  Double jacDet;

  // derivation of shape functions after global coordinates
  Vector<Double> Vx, Vy;
  Double VL2, VL2AtIp, VMax;
  Double& tau_mp = stabilParams_[elemNumber-1][0];
  Double& tau_mu = stabilParams_[elemNumber-1][1];
  Double& tau_c  = stabilParams_[elemNumber-1][2];
  Double A_elem, h_k;
  Double lambda_k;
  bool computeTaus;
  const UInt N = spaceDim_;
  if (!allocMem)
  {
    doMemAlloc(nrFncs);
  }

  Vx.Replace(nrFncs, elemResult_[0], false);
  Vy.Replace(nrFncs, elemResult_[1], false);

  if (stabilParams_[elemNumber-1][3] < -0.001)
  {
    computeTaus = true;
  } else{
    computeTaus = false;
  }

  if (computeTaus)
  {
    A_elem = ptelem->CalcVolume(ptCoord_, isaxi_);
    //A_elem /= 4.0; //wg quadratische Element

      VL2 = sqrt( (Vx.NormL2()*Vx.NormL2()) + (Vy.NormL2()*Vy.NormL2()) );
      h_k = sqrt(A_elem);

    VMax=abs(Vx[0]);
    if (abs(Vy[0])>VMax)
      VMax=abs(Vy[0]);

    for (UInt i=1; i<nrFncs; i++)
    {
      if (abs(Vx[i])>abs(Vy[i]) && abs(Vx[i])>VMax )
        VMax=abs(Vx[i]);
      else if (abs(Vx[i])<abs(Vy[i]) && abs(Vy[i])>VMax )
        VMax=abs(Vy[i]);
    }

    VL2AtIp = -1.0;
    ComputeTaus(tau_mp,tau_mu,tau_c,VL2,VL2AtIp,VMax,h_k,kinematicViscosity_,lambda_k,nrFncs);
    stabilParams_[elemNumber-1][3]=1.0;
  }
  // set matrix to desired size and set all elements to zero
  //elemMat.Resize(nrFncs * N); elemMat.Init();
  locElemMat_UQ_.Init();

  //#pragma omp parallel for private(auxJ_,VxAtIP,VyAtIP,VzAtIP,jacDet,xiDxDy_,xiDxxDyyDxy,xiDxx,xiDyy,xiDzz,xiDx,xiDy,xiDz,xi_,A_1_)
  for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
  {
    ptelem->GetShFncAtIp(xi_, actIntPt, ent1.GetElem());
    jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_, ent1.GetElem() );

    ptelem->GetGlobDerivShFncAtIp(xiDxDy_, actIntPt, ptCoord_, jacDet, ent1.GetElem());

    const Double jacDet_intWeight_tau = intWeights[actIntPt-1] * jacDet * tau_mp;

    for (UInt i=0; i<nrFncs; i++)
    {
      DxXi_[i] = xiDxDy_[i][0];
      DyXi_[i] = xiDxDy_[i][1];
    }
    // Mass Integrator
    auxI_  = xi_ * jacDet_intWeight_tau;

    DyadicMult(A_1_, DxXi_, auxI_);
    DyadicMult(A_2_, DyXi_, auxI_);
    // assemble element matrix
    //#pragma omp critical
    {
      AddSubMatrix(locElemMat_UQ_, A_1_, 0, 0);
      AddSubMatrix(locElemMat_UQ_, A_2_, 0, nrFncs);
    }
  }
  ColResortElementMatrix(elemMat, locElemMat_UQ_, nrFncs, nrFncs, 1, N);
}

} // end namespace CoupledField
