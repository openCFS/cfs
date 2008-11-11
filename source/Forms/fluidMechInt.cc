#include <iostream>
#include <fstream>

#include "DataInOut/ParamHandling/ParamNode.hh"

#include "fluidMechInt.hh"

namespace CoupledField
{

  FluidMechInt::FluidMechInt(Double aVal, Double bVal  )
    : BaseForm(NULL), density_ (aVal), kinematicViscosity_ (bVal)
  {
    name_ = "fluidMechInt";

    eigenValue_  = false;
    resortedWay_ = true;
    multAddOpt_  = true;
    coordUpdate_ = true;
    movingMesh_  = false;

//     params->Get( "firstDt",dt_,"transient" );
    dt_=1e-2;
//     params->Get( "type", movingMeshStr_, "fluidMech", "meshUpdate" );

//     if( movingMeshStr_=="movingMesh"){
//       movingMesh_ = true;
//       Info->PrintF( "MESH-UPDATE" , "strategy = movingMesh=true\n\n");
//     }
//     else if( movingMeshStr_=="none"){
//       movingMesh_ = false;
//       Info->PrintF( "MESH-UPDATE" , "strategy = movingMesh=false\n\n");
//     }
//     else Error("Unknown mesh updating!",__FILE__,__LINE__);
    
//     params->Get( "type", stabParamEsti_, "fluidMech", "stabParamEsti" );

//     if( stabParamEsti_=="FrancaEigVal" || stabParamEsti_=="BurdaEigVal" ){
//       eigenValue_ = true;
//     }
    stabParamEsti_="Wall";
    stabilType_="SUPG";
//     if (  params->ContainElem( keyVec,attrVec,valVec ) ) {
//       params->Get( "vMin", velNormMin_, "fluidMech", "velNormBnds" );
//       params->Get( "vMax", velNormMax_, "fluidMech", "velNormBnds" );
//       if(velNormMin_>velNormMax_)
//         Error( "The specified minimal velocity norm is larger than the maximal one! ",__FILE__,__LINE__);
//     } 
//     else {
//       Warning( "Using default upper and lower bound for the velocity norm! ",__FILE__,__LINE__);
//       velNormMin_ = 0.01;
//       velNormMax_ = 10.0;
//       std::cerr << "velNormMin_="<< velNormMin_ << "; velNormMax_="<< velNormMax_ << std::endl;
//     }
    velNormMin_ = 0.01;
    velNormMax_ = 10.0;
    Info->PrintF( "Velocity-Bounds", " VelNormMin = %e, VelNormMax = %e\n\n",velNormMin_, velNormMax_ );
      
    if(stabilType_=="GLS"){
      convStabSign_=1.0;
      viscStabSign_=-1.0;   // =0 for SUPG
      presStabSign_=-1.0;
      contiStabSign_=1.0;
      
      penaltyStabSign_=-1.0;

      presTermSign_=-1.0;
      contiTermSign_=-1.0;
      penaltyTermSign_=1.0;
    }
    else if(stabilType_=="SUPG"){
      convStabSign_=1.0;
      viscStabSign_=0.0;     //=-1.0 for GLS
      presStabSign_=-1.0;
      contiStabSign_=1.0;
      
      penaltyStabSign_=-1.0;

      presTermSign_=-1.0;
      contiTermSign_=-1.0;
      penaltyTermSign_=1.0;
    }
    else if(stabilType_=="SUPG2"){
      convStabSign_=1.0;
      viscStabSign_=0.0;     //=-1.0 for GLS
      presStabSign_=-1.0;
      contiStabSign_=1.0;

      penaltyStabSign_=-1.0;

      presTermSign_=1.0;
      contiTermSign_=1.0;
      penaltyTermSign_=1.0;
    }
    else if(stabilType_=="USFEM"){
      convStabSign_=1.0;
      viscStabSign_=1.0;
      presStabSign_=-1.0;
      contiStabSign_=1.0;

      penaltyStabSign_=-1.0;

      presTermSign_=-1.0;
      contiTermSign_=-1.0;
      penaltyTermSign_=1.0;
    }

    else{
      Error("stabilType not known",__FILE__,__LINE__);
    }

//     keyVec  = "fluidMech", "penalty";
//     attrVec = "", "";  
//     valVec  = "", "";  
//     if (  params->ContainElem( keyVec,attrVec,valVec ) ) {
//       params->Get( "value", penalty_, "fluidMech", "penalty" );
//     }
//     else 
//       penalty_=1e-3;

    penalty_=1e-3;
  }
  
  FluidMechInt::~FluidMechInt()
  {
  }


// The element matrix must be resorted because during matrix setup
// the local unknown vector is ordered as: 
//[U1 U2 U3 U4 V1 V2 V3 V4 P1 P2 P3 P4]^T
// for Node 1...4
// the needed order is
//[U1 V1 P1 U2 V2 P2 U3 V3 P3 U4 V4 P4 ]^T
  void FluidMechInt::ResortElementMatrix(Matrix<Double> & sortedElemMat, 
                                         const Matrix<Double> & elemMat, 
                                         const UInt & nrOfNodes, 
                                         const UInt & dofPerNode) {

    Integer i;
    UInt j, k, dofPerElem;
    Matrix<Double> auxMat;

    dofPerElem = nrOfNodes * dofPerNode;
    sortedElemMat.Resize(dofPerElem); 
    sortedElemMat.Init();
    auxMat.Resize(dofPerElem); 
    auxMat.Init();

    //==============================================================//
    if(nrOfNodes == 8 && dofPerNode == 3){
      Error("LookUpTable is deactivated",__FILE__,__LINE__);
      
      //check size of lookuptable and Matrix
      if((dofPerElem*dofPerElem) != lookuptable[0].size())
        {
          std::cerr << "The size of the matrix is not consistent with the lookuptable" << std::endl;
          return;
        }
#pragma omp parallel for
      for(i=0;i<lookuptable[0].size();i++)
        {
          sortedElemMat[lookuptable[2][i]][lookuptable[3][i]] = elemMat[lookuptable[0][i]][lookuptable[1][i]];
        }
    }
    else{
      // resort columns
      for (j=0; j<dofPerElem; j+=dofPerNode) {
        for (k=0; k<dofPerNode; k++) {
          for (i=0; i<dofPerElem; i++) {
            auxMat[i][j+k] = elemMat[i][j/dofPerNode + k*nrOfNodes];
          }
        }
      }
      //resort rows
      for (i=0; i<dofPerElem; i+=dofPerNode) {
        for (k=0; k<dofPerNode; k++) {
          for (j=0; j<dofPerElem; j++) {
            sortedElemMat[i+k][j] = auxMat[i/dofPerNode + k*nrOfNodes][j];
          }
        }
      }
    }
  }
// The element matrix must be resorted because during matrix setup
// the local unknown vector is ordered as: 
//[U1 U2 U3 U4 V1 V2 V3 V4 P1 P2 P3 P4]^T
// for Node 1...4
// the needed order is
//[U1 V1 P1 U2 V2 P2 U3 V3 P3 U4 V4 P4 ]^T
  void FluidMechInt::ColResortElementMatrix(Matrix<Double> & sortedElemMat, 
                                            const Matrix<Double> & elemMat, 
                                            const UInt & nrOfFctRow, 
                                            const UInt & nrOfFctCol, 
                                            const UInt & dofPerNodeRow,
                                            const UInt & dofPerNodeCol) {

    Integer i;
    UInt j, k, dofPerElemRow, dofPerElemCol;
    
    dofPerElemRow = nrOfFctRow * dofPerNodeRow;
    dofPerElemCol = nrOfFctCol * dofPerNodeCol;
    sortedElemMat.Resize(dofPerElemRow,dofPerElemCol); 
    sortedElemMat.Init();

  //==============================================================//
    // resort columns
    for (j=0; j<dofPerElemCol; j+=dofPerNodeCol) {
      for (k=0; k<dofPerNodeCol; k++) {
        for (i=0; i<dofPerElemRow; i++) {
          sortedElemMat[i][j+k] = elemMat[i][j/dofPerNodeCol + k*nrOfFctCol];
        }
      }
    }
  }
// The element matrix must be resorted because during matrix setup
// the local unknown vector is ordered as: 
//[U1 U2 U3 U4 V1 V2 V3 V4 P1 P2 P3 P4]^T
// for Node 1...4
// the needed order is
//[U1 V1 P1 U2 V2 P2 U3 V3 P3 U4 V4 P4 ]^T
  void FluidMechInt::RowResortElementMatrix(Matrix<Double> & sortedElemMat, 
                                            const Matrix<Double> & elemMat, 
                                            const UInt & nrOfFctRow, 
                                            const UInt & nrOfFctCol, 
                                            const UInt & dofPerNodeRow,
                                            const UInt & dofPerNodeCol) {

    Integer i;
    UInt j, k, dofPerElemRow, dofPerElemCol;
    
    dofPerElemRow = nrOfFctRow * dofPerNodeRow;
    dofPerElemCol = nrOfFctCol * dofPerNodeCol;
    sortedElemMat.Resize(dofPerElemRow,dofPerElemCol); 
    sortedElemMat.Init();

    //==============================================================//
    //resort rows
    for (i=0; i<dofPerElemRow; i+=dofPerNodeRow) {
      for (k=0; k<dofPerNodeRow; k++) {
        for (j=0; j<dofPerElemCol; j++) {
          sortedElemMat[i+k][j] = elemMat[i/dofPerNodeRow + k*nrOfFctRow][j];
        }
      }
    }
  }

  void FluidMechInt::ComputeTaus(Double & tau_mp,
                                 Double & tau_mu,
                                 Double & tau_c,
                                 const Double & VL2, 
                                 const Double & VL2AtIp, 
                                 const Double & VMax, 
                                 const Double & h_k, 
                                 const Double & kinematicViscosity, 
                                 const Double & lambda_k,
                                 const UInt & nrFncs) {
    Double VelNorm, usedNorm;
    Double Re_k, zeta, lambda, m_k;
    Double freq1, freq2, freto;
    Double C1, C2, C3, C4, C5;
    Double myh_k=0.1;

    usedNorm = VL2;

    if (usedNorm < velNormMin_)
      VelNorm = velNormMin_;
    else if (usedNorm > velNormMax_)
      VelNorm = velNormMax_;
    else
      VelNorm = usedNorm;
    
    //VelNorm = 1.0;
    //std::cout << "V=" << VelNorm << std::endl; 

    if( stabParamEsti_=="Wall" ){

      lambda = 1.0;

      //linear elements P1:
      if(nrFncs==3)
        m_k = 1.0/3.0;

      //linear elements Q1:
      else if(nrFncs==4)
        m_k = 1.0/3.0;

      //quadratic elements P2:
      else if(nrFncs==6)
        m_k = 1.0/12.0;

      //quadratic elements Q2:
      else if(nrFncs==8)
        m_k = 1.0/12.0;
      else 
        Error("unknown element!!!",__FILE__,__LINE__);
      
      //Double C_k = 1.0/(lambda_k*h_k*h_k);
      //Double h_kF = 100.0*h_k;
      Re_k =  (m_k * VelNorm * h_k) / (2 * kinematicViscosity);

      if (Re_k >= 1e-12 && Re_k < 1.0)
        zeta = Re_k;
//       else if (Re_k < 1e-12)
//         zeta = 1e-12;
      else if (Re_k >= 1.0)
        zeta = 1.0;
      else{
	std::cerr << "m_k=" << m_k << "; VelNorm=" << VelNorm << "; h_k=" << h_k << "; kinematicViscosity=" << kinematicViscosity <<std::endl;
        std::cerr << "VL2=" << VL2 << "; VL2AtIp=" << VL2AtIp << "; VMax=" << VMax << std::endl;
        Error( "Re_k value is smaller than 1e-12;\n\n\n !!! INCREAS VMIN !!!\n", __FILE__, __LINE__);
      }

      tau_mp  = ( h_k * zeta ) / ( 2.0 * VelNorm );
      tau_c =  ( lambda * VelNorm * h_k * zeta ) / (2.0);

      if ( tau_mp > ( h_k )/(2.0*VelNorm) ) {
        //std::cout << "A\n";
        tau_mp = ( h_k )/(2.0*VelNorm);
      }
      //if ( tau_mp > dt_ ) {
      //std::cout << "!B\n";
      //  tau_mp = dt_;
      //}
      tau_mu  = tau_mp;
     }

    else if( stabParamEsti_=="Foerster" ){

      //linear elements P1:
      if(nrFncs==3)
        m_k = 1.0/3.0;

      //linear elements Q1:
      else if(nrFncs==4)
        m_k = 1.0/3.0;

      //quadratic elements P2:
      else if(nrFncs==6)
        m_k = 1.0/12.0;

      //quadratic elements Q2:
      else if(nrFncs==8)
        m_k = 1.0/12.0;
      else 
        Error("unknown element!!!",__FILE__,__LINE__);

      Double zeta1, zeta2;
      Double delta=0.5*dt_;

      Double Pe_k = (4.0*delta*kinematicViscosity)/(m_k * h_k*h_k);
      Re_k =  (m_k * VelNorm * h_k) / (2.0 * kinematicViscosity);


      if (Re_k >= 0.0 && Re_k < 1.0)
        zeta2 = Re_k;
      else if (Re_k >= 1.0)
        zeta2 = 1.0;
      else
        Error( "Re_k value is smaller than 0!", __FILE__, __LINE__);

      if (Pe_k >= 0.0 && Pe_k < 1.0)
        zeta1 = Pe_k;
      else if (Pe_k >= 1.0)
        zeta1 = 1.0;
      else
        Error( "Pe_k value is smaller than 0!", __FILE__, __LINE__);

          
      tau_mp  = ( h_k * h_k ) / ( (h_k*h_k*zeta1/delta)+(4*kinematicViscosity*zeta2/m_k) );

      tau_c =  ( delta * VelNorm * h_k * zeta2 ) / (2.0);

      tau_mu  = tau_mp;
     }

    else if ( stabParamEsti_=="Codina" ){


//     VelNorm = 1.0;

      C1=1.0;
      C2=1.0;
      C4=C1;
      C5=C2;

      freq1 = C1*4.0*kinematicViscosity/(myh_k*myh_k); //Stokes term
      freq2 = C2*2.0*VelNorm/myh_k;     //Convective term
 
      tau_mp = 0.0;
      tau_c = 0.0;
      freto = freq1+freq2;
      //freto = (freq1*freq1)+(freq2*freq2);
      if(freto >= 1e-20){
        tau_mp = 1.0/freto;
        tau_mu = tau_mp;
        freto = 1.0/(freto*myh_k*myh_k);
          //  freto /= (myh_k*myh_k);
        tau_c = 1.0/freto;
      }
    }
    else if ( stabParamEsti_=="FrancaEigVal" ){

      Re_k =  (VelNorm) / (4.0 * sqrt(lambda_k) * kinematicViscosity);

      if (Re_k >= 0.0 && Re_k < 1.0)
        zeta = Re_k;
      else if (Re_k >= 1.0)
        zeta = 1.0;
      else
        Error( "Re_k value is smaller than 1e-9!", __FILE__, __LINE__);

      tau_mp = zeta / (sqrt(lambda_k)*VelNorm );
      tau_mu = tau_mp;
      tau_c = (VelNorm*zeta) / sqrt(lambda_k);
    }
    else if ( stabParamEsti_=="BurdaEigVal" ){

      C1 = 1.0 / (4.0 * sqrt(lambda_k)*kinematicViscosity);
      C2 = 1.0 / (4.0*lambda_k*kinematicViscosity);
      C3 = 1.0 / sqrt(lambda_k);

      Re_k =  C1*VelNorm;
      
      if (Re_k >= 0.0 && Re_k < 1.0)
        tau_mp = C2;
      else if (Re_k >= 1.0){
        tau_mp = C3/VelNorm;
      }
      else
        Error( "Re_k value is smaller than 1e-9!", __FILE__, __LINE__);

      tau_mu = tau_mp;
      tau_c = 0.0;
    }
    else if ( stabParamEsti_=="Franca" ){
      lambda=1.0;
      m_k = 1.0/12.0;
      
      Double C_k = 1.0/(lambda_k*h_k*h_k);

      if ( m_k > 2.0*C_k) {
        Warning( "Using C_k !", __FILE__, __LINE__);
        m_k = 2.0*C_k;
      }
      
      Re_k =  (m_k * VelNorm * h_k) / (4.0 * kinematicViscosity);

      if (Re_k >= 0.0 && Re_k < 1.0)
        zeta = Re_k;
      else if (Re_k >= 1.0)
        zeta = 1.0;
      else
        Error( "Re_k value is smaller than 1e-9!", __FILE__, __LINE__);
      
      tau_mp  = ( h_k * zeta ) / ( 2.0 * VelNorm );
      tau_c =  ( lambda * VelNorm * h_k * zeta );

      if ( tau_mp > (( m_k*h_k*h_k )/(8.0*kinematicViscosity)) ) {
        Warning( "upper bound reached !", __FILE__, __LINE__);
        tau_mp = (( m_k*h_k*h_k )/(8.0*kinematicViscosity));
      }

      tau_mu  = tau_mp;
      
    }
//    else if ( stabParamEsti_=="specified" ){
//
//      tau_mp = tau_mpFac_;
//      tau_mu = tau_muFac_;
//      tau_c  = tau_cFac_;
//      
//    }
    else
      Error( "So what: Wall, Codina, Franca or generalized Eigenvalue???", __FILE__, __LINE__);

    //Info->PrintF( stabParamEsti_, "tau_mp = %e   tau_c = %e   h_k = %e   VelNorm = %e   \n", tau_mp, tau_c, h_k, VelNorm);
  }

//   void FluidMechInt::PrepareFluidMechIntegrators(const UInt nrFncs, 
//                                                  const UInt nrIntPts) {

//     // Extract pointer to reference element and get coordinates
//     ExtractElemInfo( ent1 );
//     elemNumber_ = ent1.GetElem()->elemNum;
    
//     sol_->GetElemSolutionAsMatrix( elemResult_, ent1 , 0);

//     if(movingMesh_){
//       gridSol_->GetElemSolutionAsMatrix( gridElemResult_, ent1 , 0);
//       elemResult_-=gridElemResult_;
//     }

//     ptelem->SetAnsatzFct( ansatzFct1_ );
    
//     if(stabilParams_[elemNumber_-1][3] < -0.001)
//       computeTaus_=true;
//     else
//       computeTaus_=false;
    

//     Vx_.Resize(nrFncs);//Vx_.Init(0.0);
//     Vy_.Resize(nrFncs);//Vy_.Init(0.0);
    
//     for (UInt i=0; i<nrFncs; i++) {
//       Vx_[i]=elemResult_[0][i];
//       Vy_[i]=elemResult_[1][i];
//     }

//     if( computeTaus_ && eigenValue_ )
//       lambda_k_ = GetLambdaMaxPlaneLapl((*ptelem), ptCoord_, ent1);
    
//     if(computeTaus_){
//       A_elem_ = ptelem->CalcVolume(ptCoord_, isaxi_);
//       //A_elem_ /= 4.0; //wg quadratische Element
//       h_k_ = sqrt(A_elem  / (PI) );

//       VL2_ = sqrt( (Vx.NormL2()*Vx.NormL2()) + (Vy.NormL2()*Vy.NormL2()) );

//       VMax_=abs(Vx_[0]);
//       if (abs(Vy_[0])>VMax_)
// 	VMax_=abs(Vy_[0]);
//       for (UInt i=1; i<nrFncs; i++) {
// 	if (abs(Vx_[i])>abs(Vy_[i]) && abs(Vx_[i])>VMax_ )
// 	  VMax_=abs(Vx_[i]);
// 	else if (abs(Vx_[i])<abs(Vy_[i]) && abs(Vy_[i])>VMax_ )
// 	  VMax_=abs(Vy_[i]);
//       }
      
//       VL2AtIp_ = -1.0;
//       ComputeTaus(tau_m_,tau_mu_,tau_c_,VL2_,VL2AtIp_,VMax_,h_k_,kinematicViscosity_,lambda_k_);
//       stabilParams_[elemNumber_-1][0]=tau_m_;
//       stabilParams_[elemNumber_-1][1]=tau_mu_;
//       stabilParams_[elemNumber_-1][2]=tau_c_;
//       stabilParams_[elemNumber_-1][3]=1.0;
//     }
//     else{
//       tau_m_  = stabilParams_[elemNumber_-1][0];
//       tau_mu_ = stabilParams_[elemNumber_-1][1];
//       tau_c_  = stabilParams_[elemNumber_-1][2];
//     }
//   }
} // end namespace CoupledField
