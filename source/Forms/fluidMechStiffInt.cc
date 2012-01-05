#include <cmath>
#include <map>
#include <utility>

#include "Domain/elem.hh"
#include "Domain/entityList.hh"
#include "Elements/basefe.hh"
#include "Forms/fluidMechInt.hh"
#include "General/environment.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "PDE/timestepping.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/tools.hh"
#include "fluidMechStiffInt.hh"


namespace CoupledField {

/***********************************************************
 * constructors
 **********************************************************/
FluidMechPlaneStiffNewtonInt_UV::FluidMechPlaneStiffNewtonInt_UV(Double density, \
  Double kinematicViscosity, Matrix<Double> & stabilParams, \
  bool movingMesh, std::string stabilType)
: FluidMechInt(density, kinematicViscosity, movingMesh, stabilType), \
  stabilParams_(stabilParams)
{
  isSolDependent_ = true;
  isaxi_ = false;
  name_ = "FluidMechPlaneStiffNewtonInt_UV";
}

FluidMechPlaneStiffNewtonInt_PV::FluidMechPlaneStiffNewtonInt_PV(Double density, \
  Double kinematicViscosity, Matrix<Double> & stabilParams, \
  bool movingMesh, std::string stabilType)
: FluidMechInt(density, kinematicViscosity,movingMesh, stabilType), \
  stabilParams_(stabilParams)
{
  isSolDependent_ = true;
  isaxi_ = false;
  name_ = "FluidMechPlaneStiffNewtonInt_PV";
}

FluidMechPlaneStiffNewtonInt_UQ::FluidMechPlaneStiffNewtonInt_UQ(Double density, \
    Double kinematicViscosity, Matrix<Double> & stabilParams, \
    bool movingMesh, std::string stabilType)
: FluidMechInt(density,kinematicViscosity,movingMesh, stabilType), \
  stabilParams_(stabilParams)
{
  isSolDependent_ = true;
  isaxi_ = false;
  name_ = "FluidMechPlaneStiffNewtonInt_UQ";
}

FluidMechPlaneStiffNewtonInt_PQ::FluidMechPlaneStiffNewtonInt_PQ(Double density, \
    Double kinematicViscosity, Matrix<Double> & stabilParams, \
    bool movingMesh, std::string stabilType)
: FluidMechInt(density,kinematicViscosity,movingMesh, stabilType), \
  stabilParams_(stabilParams)
{
  isSolDependent_ = true;
  isaxi_ = false;
  name_ = "FluidMechPlaneStiffNewtonInt_PQ";
}

FluidMechPlaneIntNewton_RhsUV::FluidMechPlaneIntNewton_RhsUV(Double density, \
    Double kinematicViscosity, Matrix<Double> & stabilParams, \
    bool movingMesh, std::string stabilType)
: FluidMechInt(density, kinematicViscosity,movingMesh, stabilType), \
  stabilParams_(stabilParams)
{
  isSolDependent_ = true;
  isaxi_ = false;
  name_ = "FluidMechPlaneIntNewton_RhsUV";
  baseType_ = AUXILIARY;
}

FluidMechPlaneIntNewton_RhsUQ::FluidMechPlaneIntNewton_RhsUQ(Double density, \
    Double kinematicViscosity, Matrix<Double> & stabilParams, \
    bool movingMesh, std::string stabilType)
: FluidMechInt(density, kinematicViscosity,movingMesh, stabilType), \
  stabilParams_(stabilParams)
{
  isSolDependent_ = true;
  isaxi_ = false;
  name_ = "FluidMechPlaneIntNewton_RhsUQ";
  baseType_ = AUXILIARY;
}

/***********************************************************
 * destructors
 **********************************************************/
FluidMechPlaneStiffNewtonInt_UV::~FluidMechPlaneStiffNewtonInt_UV()
{
}

FluidMechPlaneStiffNewtonInt_PV::~FluidMechPlaneStiffNewtonInt_PV()
{
}

FluidMechPlaneStiffNewtonInt_UQ::~FluidMechPlaneStiffNewtonInt_UQ()
{
}

FluidMechPlaneStiffNewtonInt_PQ::~FluidMechPlaneStiffNewtonInt_PQ()
{
}

FluidMechPlaneIntNewton_RhsUV::~FluidMechPlaneIntNewton_RhsUV()
{
}

FluidMechPlaneIntNewton_RhsUQ::~FluidMechPlaneIntNewton_RhsUQ()
{
}

/***********************************************************
 * calc matrix
 **********************************************************/
void FluidMechPlaneStiffNewtonInt_UV::CalcElementMatrix( Matrix<Double>& elemMat, \
    EntityIterator& ent1, EntityIterator& ent2 )
{
  // Extract pointer to reference element and get coordinates
  ExtractElemInfo( ent1 );
  UInt elemNumber = ent1.GetElem()->elemNum;

  sol_->GetElemSolutionAsMatrix( elemResult_, ent1 );

  //adjust viscosity
  if ( isSpongeLayer_) {
    //get center coordinate
    Point center;
    ptelem->CalcBarycenter(ptCoord_, center);
    //take here hard coded the x-coordinate
    Double pos = (double)center[0];
    AdjustViscosity( pos );
  }

  Matrix<Double> gridElemResult;
  if (movingMesh_) {
    gridSol_->GetElemSolutionAsMatrix( gridElemResult, ent1);
    elemResult_ -= gridElemResult;
  }

  ptelem->SetAnsatzFct( ansatzFct1_ );
  const UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
  const UInt nrIntPts = ptelem->GetNumIntPoints();

  bool computeTaus = false;
  if (stabilParams_[elemNumber-1][3] < -0.001)
  {
    computeTaus = true;
  }

  Vector<Double> Vx, Vy;
  Double VxAtIP, VyAtIP, VL2, VL2AtIp, VMax;
  Double VxDxAtIP, VxDyAtIP, VxDxxAtIP, VxDyyAtIP;
  Double VyDxAtIP, VyDyAtIP, VyDxxAtIP, VyDyyAtIP;
  Double preFactor;
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

  const Vector<Double> & intWeights = ptelem->GetIntWeights();
  Double jacDet, jacDet_intWeight;

  Double lambda_k;

  if (computeTaus)
  {
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

  // time stepping variables
  const TimeStepping* const ts_alg = getTimeStepping();
  const Double ts_coefficients  = ts_alg->GetCoefficients().find(TIMESTEP_0)->second;

  //#pragma omp parallel for private(jacDet, jacDet_intWeight, preFactor, \
      VxAtIP, VxAtIP_multAux, VxDxAtIP, VxDyAtIP, VxDxxAtIP, VxDyyAtIP, \
      VyAtIP, VyAtIP_multAux, VyDxAtIP, VyDyAtIP, VyDxxAtIP, VyDyyAtIP, \
      auxI_, xiDxDy_, xiDxxDyyDxy_, \
      A_11_, A_12_, A_21_, A_22_, A_tmp_) \
      firstprivate(DxxXi_, DyyXi_, DxXi_, DyXi_, xi_, Vx, Vy) \
      shared(tau_mu, locElemMat, ts_coefficients, ent1)
  for (UInt actIntPt=1; actIntPt <= nrIntPts; ++actIntPt)
  {
    jacDet = 0.0;

//#pragma omp critical (uv_get_elem_vars)
    {
      ptelem->GetShFncAtIp(xi_, actIntPt, ent1.GetElem());
      ptelem->GetGlobDerivShFncAtIp(xiDxDy_, actIntPt, ptCoord_, jacDet, ent1.GetElem());
      ptelem->GetGlob2ndDerivShFncAtIp(xiDxxDyyDxy_, xiDxDy_, actIntPt, ptCoord_, ent1.GetElem());
    }

    jacDet_intWeight = intWeights[actIntPt-1] * jacDet;
    const Double jacDet_tauMu = jacDet_intWeight * tau_mu;

    for (UInt i=0; i<nrFncs; ++i)
    {
      DxXi_[i] = xiDxDy_[i][0];
      DyXi_[i] = xiDxDy_[i][1];

      DxxXi_[i] = xiDxxDyyDxy_[i][0];
      DyyXi_[i] = xiDxxDyyDxy_[i][1];
    }

    Vx.Inner(xi_, VxAtIP);
    VxAtIP_multAux = VxAtIP * tau_mu * jacDet_intWeight;
    Vy.Inner(xi_, VyAtIP);
    VyAtIP_multAux = VyAtIP * tau_mu * jacDet_intWeight;

    Vx.Inner(DxXi_,  VxDxAtIP);
    Vx.Inner(DyXi_,  VxDyAtIP);
    Vx.Inner(DyyXi_, VxDxxAtIP);
    Vx.Inner(DyyXi_, VxDyyAtIP);
    Vy.Inner(DxXi_,  VyDxAtIP);
    Vy.Inner(DyXi_,  VyDyAtIP);
    Vy.Inner(DyyXi_, VyDxxAtIP);
    Vy.Inner(DyyXi_, VyDyyAtIP);

    /*** A_11_ part ***/
    /* All bilinear forms with (. , DxXi_) */
    //(stab viscous term + stab diffusion + pres) * xi_i
    preFactor  = VxAtIP * ts_coefficients;
    preFactor += VxAtIP * VxDxAtIP;
    preFactor *= jacDet_tauMu;

    auxI_  = xi_ * preFactor;
    auxI_ += DxXi_ * (VxAtIP_multAux * VxAtIP);
    auxI_ += DyXi_ * (VxAtIP_multAux * VyAtIP);
    auxI_ -= DxxXi_ * (VxAtIP_multAux * kinematicViscosity_);
    auxI_ -= DyyXi_ * (VxAtIP_multAux * kinematicViscosity_);
    preFactor  = kinematicViscosity_ + tau_c;
    preFactor *= jacDet_intWeight;
    auxI_ += DxXi_ * preFactor;

    DyadicMult(A_11_, DxXi_, auxI_, nrFncs, nrFncs);

    /* All bilinear forms with (. , DyXi_) */
    auxI_  = xi_ * (VyAtIP_multAux * ts_coefficients);
    auxI_ += xi_ * (VyAtIP_multAux * VxDxAtIP);
    auxI_ += DxXi_ * (VxAtIP_multAux * VyAtIP);
    auxI_ += DyXi_ * (VyAtIP_multAux * VyAtIP);
    auxI_ -= DxxXi_ * (VyAtIP_multAux * kinematicViscosity_);
    auxI_ -= DyyXi_ * (VyAtIP_multAux * kinematicViscosity_);
    auxI_ += DyXi_ * (kinematicViscosity_ * jacDet_intWeight);

    DyadicMult(A_tmp_, DyXi_, auxI_, nrFncs, nrFncs);
    A_11_ += A_tmp_;

    /* Bilinear from (., xi_) -> timestep + viscous form*/
    preFactor  = ts_coefficients + VxDxAtIP;
    preFactor *= jacDet_intWeight;
    auxI_  = xi_ * preFactor;
    auxI_ += DxXi_ * (VxAtIP * jacDet_intWeight);
    auxI_ += DyXi_ * (VyAtIP * jacDet_intWeight);

    DyadicMult(A_tmp_, xi_, auxI_, nrFncs, nrFncs);
    A_11_ += A_tmp_;

    /*** A_12_ part ***/
    /* All bilinear forms with (. , DyXi_) */
    preFactor  = VyAtIP * VxDyAtIP;
    preFactor *= jacDet_tauMu;
    auxI_ = xi_ * preFactor;
    DyadicMult(A_12_, DyXi_, auxI_, nrFncs, nrFncs);

    /* the remaining bilinear forms */
    auxI_ = xi_ * (VxDyAtIP * jacDet_intWeight);
    DyadicMult(A_tmp_, xi_, auxI_, nrFncs, nrFncs);
    A_12_ += A_tmp_;

    auxI_  = xi_ * (VxAtIP_multAux * VxDyAtIP);
    auxI_ += DyXi_ * (tau_c * jacDet_intWeight);
    DyadicMult(A_tmp_, DxXi_, auxI_, nrFncs, nrFncs);
    A_12_ += A_tmp_;

    /*** A_21_ part ***/
    /* All bilinear forms with (. , DxXi_) */
    preFactor  = VxAtIP * VyDxAtIP;
    preFactor *= jacDet_tauMu;
    auxI_ = xi_ * preFactor;
    DyadicMult(A_21_, DxXi_, auxI_, nrFncs, nrFncs);

    /* the remaining bilinear forms */
    auxI_ = xi_ * (VyDxAtIP * jacDet_intWeight);
    DyadicMult(A_tmp_, xi_, auxI_, nrFncs, nrFncs);
    A_21_ += A_tmp_;

    auxI_  = xi_ * (VyAtIP_multAux * VyDxAtIP);
    auxI_ += DxXi_ * (tau_c * jacDet_intWeight);
    DyadicMult(A_tmp_, DyXi_, auxI_, nrFncs, nrFncs);
    A_21_ += A_tmp_;

    /*** A_22_ part ***/
    /* All bilinear forms with (. * xi_ , DyXi_) */
    //(stab viscous term + stab diffusion + pres) * xi_i
    preFactor  = VyAtIP * ts_coefficients;
    preFactor += VyAtIP * VyDyAtIP;
    preFactor *= jacDet_tauMu;

    auxI_  = xi_ * preFactor;
    auxI_ += DyXi_ * (VyAtIP_multAux * VyAtIP);
    auxI_ += DxXi_ * (VyAtIP_multAux * VxAtIP);
    auxI_ -= DxxXi_ * (VyAtIP_multAux * kinematicViscosity_);
    auxI_ -= DyyXi_ * (VyAtIP_multAux * kinematicViscosity_);
    preFactor  = kinematicViscosity_ + tau_c;
    preFactor *= jacDet_intWeight;
    auxI_ += DyXi_ * preFactor;

    DyadicMult(A_22_, DyXi_, auxI_, nrFncs, nrFncs);

    /* All bilinear forms with (. , DxXi_) */
    auxI_  = xi_ * (VxAtIP_multAux * ts_coefficients);
    auxI_ += xi_ * (VxAtIP_multAux * VyDyAtIP);
    auxI_ += DyXi_ * (VxAtIP_multAux * VyAtIP);
    auxI_ += DxXi_ * (VxAtIP_multAux * VxAtIP);
    auxI_ -= DxxXi_ * (VxAtIP_multAux * kinematicViscosity_);
    auxI_ -= DyyXi_ * (VxAtIP_multAux * kinematicViscosity_);
    auxI_ += DxXi_ * (kinematicViscosity_ * jacDet_intWeight);

    DyadicMult(A_tmp_, DxXi_, auxI_, nrFncs, nrFncs);
    A_22_ += A_tmp_;

    /* Bilinear from (., xi_) -> timestep + viscous form*/
    preFactor  = ts_coefficients + VyDyAtIP;
    preFactor *= jacDet_intWeight;
    auxI_  = xi_ * preFactor;
    auxI_ += DxXi_ * (VxAtIP * jacDet_intWeight);
    auxI_ += DyXi_ * (VyAtIP * jacDet_intWeight);

    DyadicMult(A_tmp_, xi_, auxI_, nrFncs, nrFncs);
    A_22_ += A_tmp_;

    //#pragma omp critical (uv_form_sum)
    {
      AddSubMatrix(locElemMat_UV_, A_11_, 0, 0);
      AddSubMatrix(locElemMat_UV_, A_12_, 0, nrFncs);
      AddSubMatrix(locElemMat_UV_, A_21_, nrFncs, 0);
      AddSubMatrix(locElemMat_UV_, A_22_, nrFncs, nrFncs);
    }
  }
  ResortElementMatrix(elemMat, locElemMat_UV_, nrFncs, spaceDim_);
}

void FluidMechPlaneStiffNewtonInt_PV::CalcElementMatrix( Matrix<Double>& elemMat, \
    EntityIterator& ent1, EntityIterator& ent2 )
{
  // Extract pointer to reference element and get coordinates
  ExtractElemInfo( ent1 );
  UInt elemNumber = ent1.GetElem()->elemNum;

  sol_->GetElemSolutionAsMatrix( elemResult_, ent1 );

  //adjust viscosity
  if ( isSpongeLayer_) {
    //get center coordinate
    Point center;
    ptelem->CalcBarycenter(ptCoord_, center);
    //take here hard coded the x-coordinate
    Double pos = (double)center[0];
    AdjustViscosity( pos );
  }

  if (movingMesh_)
  {
    Matrix<Double> gridElemResult_;
    gridSol_->GetElemSolutionAsMatrix( gridElemResult_, ent1);
    elemResult_-=gridElemResult_;
  }

  ptelem->SetAnsatzFct( ansatzFct1_ );
  const UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
  const UInt nrIntPts= ptelem->GetNumIntPoints();

  Vector<Double>  Vx, Vy;
  Double  VxAtIP_multAux, VyAtIP_multAux, VL2, VL2AtIp, VMax;
  Double& tau_mp  = stabilParams_[elemNumber-1][0];
  Double& tau_mu = stabilParams_[elemNumber-1][1];
  Double& tau_c  = stabilParams_[elemNumber-1][2];
  Double A_elem, h_k;

  bool computeTaus;
  if (stabilParams_[elemNumber-1][3] < -0.001)
  {
    computeTaus = true;
  } else{
    computeTaus = false;
  }

  Vx.Replace(nrFncs, elemResult_[0], false);
  Vy.Replace(nrFncs, elemResult_[1], false);

  const Vector<Double> & intWeights = ptelem->GetIntWeights();
  Double jacDet, jacDet_intWeight;
  const UInt N = spaceDim_;  // DOFs per Node (Ux, Uy)

  Double lambda_k;

  if (!allocMem)
  {
    doMemAlloc(nrFncs);
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
  //elemMat.Resize(nrFncs*N,nrFncs); elemMat.Init();
  locElemMat_PV_.Init();

  //#pragma omp parallel for private(auxI_,auxJ_,VxAtIP,VyAtIP,VzAtIP,jacDet,jacDet_intWeight,multAux,xiDxDy_,xiDxxDyyDxy_,DxxXi_,DyyXi_,xiDzz,DxXi_,DyXi_,xiDz,xi_,B1,B2,B3,B_aa,B_ab,B_ac,B_bc1,B_bc2,B_bc3)
  for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
  {
    jacDet = 0;

    ptelem->GetShFncAtIp(xi_, actIntPt, ent1.GetElem());
    ptelem->GetGlobDerivShFncAtIp(xiDxDy_, actIntPt, ptCoord_, jacDet, ent1.GetElem());
    ptelem->GetGlob2ndDerivShFncAtIp(xiDxxDyyDxy_, xiDxDy_, actIntPt, ptCoord_, ent1.GetElem());

    jacDet_intWeight = intWeights[actIntPt-1] * jacDet;

    for (UInt i=0; i<nrFncs; i++)
    {
      DxXi_[i] = xiDxDy_[i][0];
      DyXi_[i] = xiDxDy_[i][1];
    }

    Vx.Inner(xi_, VxAtIP_multAux);
    Vy.Inner(xi_, VyAtIP_multAux);
    VxAtIP_multAux *= tau_mu * jacDet_intWeight;
    VyAtIP_multAux *= tau_mu * jacDet_intWeight;

    // Pressure Integrator
#if 0 // DEBUG_SZOERNER
    auxI_ = xi_ * jacDet_intWeight * (-1.0);
    DyadicMult(A_1_, DxXi_, auxI_, nrFncs, nrFncs);

    auxI_ = xi_ * jacDet_intWeight * (-1.0);
    DyadicMult(A_2_, DyXi_, auxI_, nrFncs, nrFncs);
#else
    auxI_ = DxXi_ * jacDet_intWeight;
    DyadicMult(A_1_, xi_, auxI_, nrFncs, nrFncs);

    auxI_ = DyXi_ * jacDet_intWeight;
    DyadicMult(A_2_, xi_, auxI_, nrFncs, nrFncs);
#endif

    //Stabilisation Integrators
    auxI_ = DxXi_ * VxAtIP_multAux;
    DyadicMult(A_tmp_, DxXi_, auxI_, nrFncs, nrFncs);
    A_1_ += A_tmp_;
    auxI_ = DxXi_ * VyAtIP_multAux;
    DyadicMult(A_tmp_, DyXi_, auxI_, nrFncs, nrFncs);
    A_1_ += A_tmp_;

    auxI_ = DyXi_ * VxAtIP_multAux;
    DyadicMult(A_tmp_, DxXi_, auxI_, nrFncs, nrFncs);
    A_2_ += A_tmp_;
    auxI_ = DyXi_ * VyAtIP_multAux;
    DyadicMult(A_tmp_, DyXi_, auxI_, nrFncs, nrFncs);
    A_2_ += A_tmp_;

    //#pragma omp critical
    {
      AddSubMatrix(locElemMat_PV_, A_1_, 0, 0);
      AddSubMatrix(locElemMat_PV_, A_2_, nrFncs, 0);
    }
  } // loop over integration points
  RowResortElementMatrix(elemMat, locElemMat_PV_, nrFncs, nrFncs, N, 1);
}

void FluidMechPlaneStiffNewtonInt_UQ::CalcElementMatrix( Matrix<Double>& elemMat, \
    EntityIterator& ent1, EntityIterator& ent2 )
{
  // Extract pointer to reference element and get coordinates
  ExtractElemInfo( ent1 );
  UInt elemNumber = ent1.GetElem()->elemNum;

  sol_->GetElemSolutionAsMatrix( elemResult_, ent1);

  //adjust viscosity
  if ( isSpongeLayer_) {
    //get center coordinate
    Point center;
    ptelem->CalcBarycenter(ptCoord_, center);
    //take here hard coded the x-coordinate
    Double pos = (double)center[0];
    AdjustViscosity( pos );
  }

  if (movingMesh_)
  {
    Matrix<Double> gridElemResult_;
    gridSol_->GetElemSolutionAsMatrix( gridElemResult_, ent1 );
    elemResult_ -= gridElemResult_;
  }

  ptelem->SetAnsatzFct( ansatzFct1_ );
  const UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
  const UInt nrIntPts= ptelem->GetNumIntPoints();

  Vector<Double> Vx, Vy;
  Double VxAtIP, VyAtIP, VL2, VL2AtIp, VMax;
  Double VxDxAtIP, VxDyAtIP;
  Double VyDxAtIP, VyDyAtIP;
  Double preFactor;
  Double& tau_mp  = stabilParams_[elemNumber-1][0];
  Double& tau_mu = stabilParams_[elemNumber-1][1];
  Double& tau_c  = stabilParams_[elemNumber-1][2];
  Double A_elem, h_k;
  bool computeTaus;

  if (!allocMem)
  {
    doMemAlloc(nrFncs);
  }

  if (stabilParams_[elemNumber-1][3] < -0.001)
  {
    computeTaus = true;
  } else{
    computeTaus = false;
  }

  Vx.Replace(nrFncs, elemResult_[0], false);
  Vy.Replace(nrFncs, elemResult_[1], false);

  const Vector<Double> & intWeights = ptelem->GetIntWeights();
  Double jacDet, jacDet_intWeight;
  const UInt N = spaceDim_;  // DOFs per Node (Ux, Uy)

  Double lambda_k;

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
  locElemMat_UQ_.Init();

  // time stepping variables
  const TimeStepping* const ts_alg = getTimeStepping();
  const Double& ts_coefficients  = ts_alg->GetCoefficients().find(TIMESTEP_0)->second;

  //#pragma omp parallel for private(auxI_,VxAtIP,VyAtIP,VzAtIP,jacDet,jacDet_intWeight,multAux,xiDxDy_,xiDxxDyyDxy_,DxxXi_,DyyXi_,xiDzz,DxXi_,DyXi_,xiDz,xi_,D1,D2,D3,D_aa,D_ab,D_ac,D_bc1,D_bc2,D_bc3)
  for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
  {
    jacDet = 0;

    ptelem->GetShFncAtIp(xi_, actIntPt, ent1.GetElem());
    ptelem->GetGlobDerivShFncAtIp(xiDxDy_, actIntPt, ptCoord_, jacDet, ent1.GetElem());
    if (abs(presStabSign_)>1e-5)
      ptelem->GetGlob2ndDerivShFncAtIp(xiDxxDyyDxy_, xiDxDy_, actIntPt, ptCoord_, ent1.GetElem());

    jacDet_intWeight = intWeights[actIntPt-1] * jacDet;
    const Double jacDet_tauMq = jacDet_intWeight * tau_mp;

    for (UInt i=0; i<nrFncs; i++)
    {
      DxXi_[i] = xiDxDy_[i][0];
      DyXi_[i] = xiDxDy_[i][1];

      DxxXi_[i] = xiDxxDyyDxy_[i][0];
      DyyXi_[i] = xiDxxDyyDxy_[i][1];
    }

    Vx.Inner(xi_, VxAtIP);
    Vy.Inner(xi_, VyAtIP);
    Vx.Inner(DxXi_, VxDxAtIP);
    Vx.Inner(DyXi_, VxDyAtIP);
    Vy.Inner(DxXi_, VyDxAtIP);
    Vy.Inner(DyXi_, VyDyAtIP);

    /** A_1_ part **/
    /* All bilinear forms with (. , DxXi_) */
    preFactor = ts_coefficients + VxDxAtIP;
    auxI_  = xi_ * preFactor;
    auxI_ += DxXi_ * VxAtIP;
    auxI_ += DyXi_ * VyAtIP;
    auxI_ -= DxxXi_ * kinematicViscosity_;
    auxI_ -= DyyXi_ * kinematicViscosity_;
    auxI_ *= jacDet_tauMq;
    DyadicMult(A_1_, DxXi_, auxI_, nrFncs, nrFncs);

    /* All bilinear forms with (. , DyXi_) */
    auxI_ = xi_ * (VyDxAtIP * jacDet_tauMq);
    DyadicMult(A_tmp_, DyXi_, auxI_, nrFncs, nrFncs);
    A_1_ += A_tmp_;

    /* All bilinear forms with (. , xi_) */
    auxI_ = DxXi_ * jacDet_intWeight;
    DyadicMult(A_tmp_, xi_, auxI_, nrFncs, nrFncs);
    A_1_ += A_tmp_;

    /** A_2_ part **/
    /* All bilinear forms with (. , DyXi_) */
    preFactor = ts_coefficients + VyDyAtIP;
    auxI_  = xi_ * preFactor;
    auxI_ += DxXi_ * VxAtIP;
    auxI_ += DyXi_ * VyAtIP;
    auxI_ -= DxxXi_ * kinematicViscosity_;
    auxI_ -= DyyXi_ * kinematicViscosity_;
    auxI_ *= jacDet_tauMq;
    DyadicMult(A_2_, DyXi_, auxI_, nrFncs, nrFncs);

    /* All bilinear forms with (. , DxXi_) */
    auxI_ = xi_ * (VxDyAtIP * jacDet_tauMq);
    DyadicMult(A_tmp_, DxXi_, auxI_, nrFncs, nrFncs);
    A_2_ += A_tmp_;

    /* All bilinear forms with (. , xi_) */
    auxI_ = DyXi_ * jacDet_intWeight;
    DyadicMult(A_tmp_, xi_, auxI_, nrFncs, nrFncs);
    A_2_ += A_tmp_;

    //#pragma omp critical
    {
      AddSubMatrix(locElemMat_UQ_, A_1_, 0, 0);
      AddSubMatrix(locElemMat_UQ_, A_2_, 0, nrFncs);
    }
  } // loop integration points
  ColResortElementMatrix(elemMat, locElemMat_UQ_, nrFncs, nrFncs, 1, N);
}

void FluidMechPlaneStiffNewtonInt_PQ::CalcElementMatrix( Matrix<Double>& elemMat, \
    EntityIterator& ent1, EntityIterator& ent2 )
{
  // Extract pointer to reference element and get coordinates
  ExtractElemInfo( ent1 );
  UInt elemNumber = ent1.GetElem()->elemNum;

  sol_->GetElemSolutionAsMatrix( elemResult_, ent1);

  //adjust viscosity
  if ( isSpongeLayer_) {
    //get center coordinate
    Point center;
    ptelem->CalcBarycenter(ptCoord_, center);
    //take here hard coded the x-coordinate
    Double pos = (double)center[0];
    AdjustViscosity( pos );
  }

  if (movingMesh_) {
    Matrix<Double> gridElemResult_;
    gridSol_->GetElemSolutionAsMatrix( gridElemResult_, ent1);
    elemResult_ -= gridElemResult_;
  }
  ptelem->SetAnsatzFct( ansatzFct1_ );
  const UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
  const UInt nrIntPts= ptelem->GetNumIntPoints();

  Vector<Double>  Vx, Vy, Vz;
  Double VL2, VL2AtIp, VMax;
  Double& tau_mp  = stabilParams_[elemNumber-1][0];
  Double& tau_mu = stabilParams_[elemNumber-1][1];
  Double& tau_c  = stabilParams_[elemNumber-1][2];
  Double A_elem, h_k;

  bool computeTaus;
  if (stabilParams_[elemNumber-1][3] < -0.001) {
    computeTaus = true;
  } else {
    computeTaus = false;
  }

  Vx.Replace(nrFncs, elemResult_[0], false);
  Vy.Replace(nrFncs, elemResult_[1], false);

  const Vector<Double> & intWeights = ptelem->GetIntWeights();
  Double jacDet;
  const UInt N = 1;  // pressure DOFs per Node

  Double lambda_k;
  if (!allocMem)
  {
    doMemAlloc(nrFncs);
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
  elemMat.Resize(nrFncs * N);
  elemMat.Init();

  for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
  {
    jacDet = 0;

    ptelem->GetGlobDerivShFncAtIp(xiDxDy_, actIntPt, ptCoord_, jacDet, ent1.GetElem());

    const Double multAux = tau_mp * intWeights[actIntPt-1] * jacDet;

    for (UInt i=0; i<nrFncs; i++)
    {
      DxXi_[i] = xiDxDy_[i][0];
      DyXi_[i] = xiDxDy_[i][1];
    }

    DyadicMult(A_1_, DxXi_, DxXi_, nrFncs, nrFncs);
    DyadicMult(A_tmp_, DyXi_, DyXi_, nrFncs, nrFncs);
    A_1_ += A_tmp_;
    A_1_ *= multAux;

    {
      AddSubMatrix(elemMat, A_1_, 0, 0);
    }
  } // loop integration points
}

void FluidMechPlaneIntNewton_RhsUV::CalcElementMatrix( Matrix<Double>& elemMat, \
    EntityIterator& ent1, EntityIterator& ent2 )
{
  // Extract pointer to reference element and get coordinates
  ExtractElemInfo( ent1 );
  UInt elemNumber = ent1.GetElem()->elemNum;

  sol_->GetElemSolutionAsMatrix( elemResult_, ent1 );

  //adjust viscosity
  if ( isSpongeLayer_) {
    //get center coordinate
    Point center;
    ptelem->CalcBarycenter(ptCoord_, center);
    //take here hard coded the x-coordinate
    Double pos = (double)center[0];
    AdjustViscosity( pos );
  }

  Matrix<Double> gridElemResult;
  if (movingMesh_) {
    gridSol_->GetElemSolutionAsMatrix( gridElemResult, ent1);
    elemResult_ -= gridElemResult;
  }

  ptelem->SetAnsatzFct( ansatzFct1_ );
  const UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
  const UInt nrIntPts = ptelem->GetNumIntPoints();

  bool computeTaus = false;
  if (stabilParams_[elemNumber-1][3] < -0.001)
  {
    computeTaus = true;
  }

  Vector<Double> Vx, Vy;
  Double VxAtIP, VyAtIP, VL2, VL2AtIp, VMax;
  Double VxAtIP_multAux, VyAtIP_multAux;
  Double VxDxAtIP, VxDyAtIP, VxDxxAtIP, VxDyyAtIP;
  Double VyDxAtIP, VyDyAtIP, VyDxxAtIP, VyDyyAtIP;
  Double preFactor;
  // need to be > 1e-15 and < 1e3
  Double& tau_mp  = stabilParams_[elemNumber-1][0];
  Double& tau_mu = stabilParams_[elemNumber-1][1];
  Double& tau_c  = stabilParams_[elemNumber-1][2];
  // Double tau_mp, tau_mu, tau_c;
  Double A_elem, h_k;

  Vx.Replace(nrFncs, elemResult_[0], false);
  Vy.Replace(nrFncs, elemResult_[1], false);

  const Vector<Double> & intWeights = ptelem->GetIntWeights();
  Double jacDet, jacDet_intWeight;

  Double lambda_k;

  if (!allocMem)
  {
    doMemAlloc(nrFncs);
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

  //#pragma omp parallel for private(auxJ_,auxI_,VxAtIP,VyAtIP,VzAtIP,jacDet,jacDet_intWeight,multAux,xiDxDy_,xiDxxDyyDxy_,DxxXi_,DyyXi_,xiDzz,DxXi_,DyXi_,xiDz,xi_,A1,A2,A3,A_a,A_ba,A_ba1,A_ba2,A_ba3,A_cdef,A_ga,A_gb,A_gc,A_gd,A_ge,A_gf,A_gg,A_gh,A_gj)
  for (UInt actIntPt=1; actIntPt <= nrIntPts; ++actIntPt)
  {
    jacDet = 0.0;

    ptelem->GetShFncAtIp(xi_, actIntPt, ent1.GetElem());
    ptelem->GetGlobDerivShFncAtIp(xiDxDy_, actIntPt, ptCoord_, jacDet, ent1.GetElem());
    ptelem->GetGlob2ndDerivShFncAtIp(xiDxxDyyDxy_, xiDxDy_, actIntPt, ptCoord_, ent1.GetElem());

    jacDet_intWeight = intWeights[actIntPt-1] * jacDet;

    for (UInt i=0; i<nrFncs; ++i)
    {
      DxXi_[i] = xiDxDy_[i][0];
      DyXi_[i] = xiDxDy_[i][1];

      DxxXi_[i] = xiDxxDyyDxy_[i][0];
      DyyXi_[i] = xiDxxDyyDxy_[i][1];
    }

    Vx.Inner(xi_, VxAtIP);
    VxAtIP_multAux = VxAtIP * tau_mu * jacDet_intWeight;
    Vy.Inner(xi_, VyAtIP);
    VyAtIP_multAux = VyAtIP * tau_mu * jacDet_intWeight;

    Vx.Inner(DxXi_,  VxDxAtIP);
    Vx.Inner(DyXi_,  VxDyAtIP);
    Vx.Inner(DyyXi_, VxDxxAtIP);
    Vx.Inner(DyyXi_, VxDyyAtIP);
    Vy.Inner(DxXi_,  VyDxAtIP);
    Vy.Inner(DyXi_,  VyDyAtIP);
    Vy.Inner(DyyXi_, VyDxxAtIP);
    Vy.Inner(DyyXi_, VyDyyAtIP);

    /*** A_11_ part ***/
    /* All bilinear forms with (. , DxXi_) */
    //(stab viscous term + stab diffusion + pres) * xi_i
    preFactor  = VxAtIP_multAux * VxDxAtIP;
    auxI_ = xi_ * preFactor;
    DyadicMult(A_11_, DxXi_, auxI_, nrFncs, nrFncs);

    /* All bilinear forms with (. , DyXi_) */
    auxI_  = xi_ * VyAtIP_multAux * VxDxAtIP;
    DyadicMult(A_tmp_, DyXi_, auxI_, nrFncs, nrFncs);
    A_11_ += A_tmp_;

    /* Bilinear from (., xi_) viscous form*/
    preFactor = VxDxAtIP * jacDet_intWeight;
    auxI_ = xi_ * preFactor;
    DyadicMult(A_tmp_, xi_, auxI_, nrFncs, nrFncs);
    A_11_ += A_tmp_;

    /*** A_12_ part ***/
    /* All bilinear forms with (. , DyXi_) */

    /* the remaining bilinear forms */
    auxI_ = xi_ * (VxDyAtIP * jacDet_intWeight);
    DyadicMult(A_12_, xi_, auxI_, nrFncs, nrFncs);

    auxI_ = xi_ * (VxAtIP_multAux * VxDyAtIP);
    DyadicMult(A_tmp_, DxXi_, auxI_, nrFncs, nrFncs);
    A_12_ += A_tmp_;

    auxI_ = xi_ * (VyAtIP_multAux * VxDyAtIP);
    DyadicMult(A_tmp_, DyXi_, auxI_, nrFncs, nrFncs);
    A_12_ += A_tmp_;

    /*** A_21_ part ***/
    /* All bilinear forms with (. , DxXi_) */

    /* the remaining bilinear forms */
    auxI_ = xi_ * (VyDxAtIP * jacDet_intWeight);
    DyadicMult(A_21_, xi_, auxI_, nrFncs, nrFncs);

    auxI_ = xi_ * (VxAtIP_multAux * VyDxAtIP);
    DyadicMult(A_tmp_, DxXi_, auxI_, nrFncs, nrFncs);
    A_21_ += A_tmp_;

    auxI_ = xi_ * (VyAtIP_multAux * VyDxAtIP);
    DyadicMult(A_tmp_, DyXi_, auxI_, nrFncs, nrFncs);
    A_21_ += A_tmp_;

    /*** A_22_ part ***/
    /* All bilinear forms with (. * xi_ , DyXi_) */
    //(stab viscous term + stab diffusion + pres) * xi_i
    preFactor  = VyAtIP_multAux * VyDyAtIP;
    auxI_ = xi_ * preFactor;
    DyadicMult(A_22_, DyXi_, auxI_, nrFncs, nrFncs);

    /* All bilinear forms with (. , DxXi_) */
    auxI_  = xi_ * (VxAtIP_multAux * VyDyAtIP);
    DyadicMult(A_tmp_, DxXi_, auxI_, nrFncs, nrFncs);
    A_22_ += A_tmp_;

    /* Bilinear from (., xi_) viscous form*/
    preFactor = VyDyAtIP * jacDet_intWeight;
    auxI_ = xi_ * preFactor;
    DyadicMult(A_tmp_, xi_, auxI_, nrFncs, nrFncs);
    A_22_ += A_tmp_;

    //#pragma omp critical
    {
      AddSubMatrix(locElemMat_UV_, A_11_, 0, 0);
      AddSubMatrix(locElemMat_UV_, A_12_, 0, nrFncs);
      AddSubMatrix(locElemMat_UV_, A_21_, nrFncs, 0);
      AddSubMatrix(locElemMat_UV_, A_22_, nrFncs, nrFncs);
    }
  }
  ResortElementMatrix(elemMat, locElemMat_UV_, nrFncs, spaceDim_);
}

void FluidMechPlaneIntNewton_RhsUQ::CalcElementMatrix( Matrix<Double>& elemMat, \
    EntityIterator& ent1, EntityIterator& ent2 )
{
  // Extract pointer to reference element and get coordinates
  ExtractElemInfo( ent1 );
  UInt elemNumber = ent1.GetElem()->elemNum;

  sol_->GetElemSolutionAsMatrix( elemResult_, ent1);

  //adjust viscosity
  if ( isSpongeLayer_) {
    //get center coordinate
    Point center;
    ptelem->CalcBarycenter(ptCoord_, center);
    //take here hard coded the x-coordinate
    Double pos = (double)center[0];
    AdjustViscosity( pos );
  }

  if (movingMesh_)
  {
    Matrix<Double> gridElemResult_;
    gridSol_->GetElemSolutionAsMatrix( gridElemResult_, ent1 );
    elemResult_ -= gridElemResult_;
  }

  ptelem->SetAnsatzFct( ansatzFct1_ );
  const UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
  const UInt nrIntPts= ptelem->GetNumIntPoints();

  Vector<Double> Vx, Vy;
  Double VL2, VL2AtIp, VMax;
  Double VxDxAtIP, VxDyAtIP;
  Double VyDxAtIP, VyDyAtIP;
  Double preFactor;
  Double& tau_mp  = stabilParams_[elemNumber-1][0];
  Double& tau_mu = stabilParams_[elemNumber-1][1];
  Double& tau_c  = stabilParams_[elemNumber-1][2];
  Double A_elem, h_k;
  bool computeTaus;

  if (!allocMem)
  {
    doMemAlloc(nrFncs);
  }

  if (stabilParams_[elemNumber-1][3] < -0.001)
  {
    computeTaus = true;
  } else{
    computeTaus = false;
  }

  Vx.Replace(nrFncs, elemResult_[0], false);
  Vy.Replace(nrFncs, elemResult_[1], false);

  const Vector<Double> & intWeights = ptelem->GetIntWeights();
  Double jacDet, jacDet_intWeight;
  const UInt N = spaceDim_;  // DOFs per Node (Ux, Uy)

  Double lambda_k;

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
  locElemMat_UQ_.Init();

  //#pragma omp parallel for private(auxI_,VxAtIP,VyAtIP,VzAtIP,jacDet,jacDet_intWeight,multAux,xiDxDy_,xiDxxDyyDxy_,DxxXi_,DyyXi_,xiDzz,DxXi_,DyXi_,xiDz,xi_,D1,D2,D3,D_aa,D_ab,D_ac,D_bc1,D_bc2,D_bc3)
  for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
  {
    jacDet = 0.0;

    ptelem->GetShFncAtIp(xi_, actIntPt, ent1.GetElem());
    ptelem->GetGlobDerivShFncAtIp(xiDxDy_, actIntPt, ptCoord_, jacDet, ent1.GetElem());

    jacDet_intWeight = intWeights[actIntPt-1] * jacDet;
    const Double jacDet_tauMq = jacDet_intWeight * tau_mp;

    for (UInt i=0; i<nrFncs; i++)
    {
      DxXi_[i] = xiDxDy_[i][0];
      DyXi_[i] = xiDxDy_[i][1];
    }

    Vx.Inner(DxXi_, VxDxAtIP);
    Vx.Inner(DyXi_, VxDyAtIP);
    Vy.Inner(DxXi_, VyDxAtIP);
    Vy.Inner(DyXi_, VyDyAtIP);

    /** A_1_ part **/
    /* All bilinear forms with (. , DxXi_) */
    preFactor = VxDxAtIP * jacDet_tauMq;
    auxI_ = xi_ * preFactor;
    DyadicMult(A_1_, DxXi_, auxI_, nrFncs, nrFncs);

    /* All bilinear forms with (. , DyXi_) */
    preFactor = VyDxAtIP * jacDet_tauMq;
    auxI_ = xi_ * preFactor;
    DyadicMult(A_tmp_, DyXi_, auxI_, nrFncs, nrFncs);
    A_1_ += A_tmp_;

    /** A_2_ part **/
    /* All bilinear forms with (. , DyXi_) */
    preFactor = VyDyAtIP * jacDet_tauMq;
    auxI_ = xi_ * preFactor;
    DyadicMult(A_2_, DyXi_, auxI_, nrFncs, nrFncs);

    /* All bilinear forms with (. , DxXi_) */
    preFactor = VxDyAtIP * jacDet_tauMq;
    auxI_ = xi_ * preFactor;
    DyadicMult(A_tmp_, DxXi_, auxI_, nrFncs, nrFncs);
    A_2_ += A_tmp_;

    //#pragma omp critical
    {
      AddSubMatrix(locElemMat_UQ_, A_1_, 0, 0);
      AddSubMatrix(locElemMat_UQ_, A_2_, 0, nrFncs);
    }
  } // loop integration points
  ColResortElementMatrix(elemMat, locElemMat_UQ_, nrFncs, nrFncs, 1, N);
}

} // end namespace CoupledField
