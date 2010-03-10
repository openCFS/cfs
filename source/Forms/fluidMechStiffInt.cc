#include <iostream>
#include <fstream>

#include "fluidMechStiffInt.hh"
#include "Utils/nodestoresol.hh"
#include "DataInOut/Logging/cfslog.hh"

namespace CoupledField {

// declare logging stream
//  DECLARE_LOG(fluidMechStiff)
//  DEFINE_LOG(fluidMechStiff, "fluidMechStiff")

//************************************************************************************************
//constructors 
//************************************************************************************************

//*******************************************************************************************************************************
//****************** stabilized FEM *********************************************************************************************
//*******************************************************************************************************************************
FluidMechPlaneStiffInt_UV::FluidMechPlaneStiffInt_UV(Double density,Double kinematicViscosity,
		Matrix<Double> & stabilParams, bool movingMesh, std::string stabilType)
: FluidMechInt(density, kinematicViscosity,movingMesh, stabilType), stabilParams_(stabilParams){
	isSolDependent_ = true;
	isaxi_=false;
	name_ = "FluidMechPlaneStiffInt_UV";
}
FluidMechPlaneStiffInt_PQ::FluidMechPlaneStiffInt_PQ(Double density,Double kinematicViscosity,
		Matrix<Double> & stabilParams, bool movingMesh, std::string stabilType)
: FluidMechInt(density,kinematicViscosity,movingMesh, stabilType), stabilParams_(stabilParams){
	isSolDependent_ = true;
	isaxi_=false;
	name_ = "FluidMechPlaneStiffInt_PQ";
}
FluidMechPlaneStiffInt_UQ::FluidMechPlaneStiffInt_UQ(Double density,Double kinematicViscosity,
		Matrix<Double> & stabilParams, bool movingMesh, std::string stabilType)
: FluidMechInt(density,kinematicViscosity,movingMesh, stabilType), stabilParams_(stabilParams) {
	isSolDependent_ = true;
	isaxi_=false;
	name_ = "FluidMechPlaneStiffInt_UQ";
}
FluidMechPlaneStiffInt_PV::FluidMechPlaneStiffInt_PV(Double density,Double kinematicViscosity,
		Matrix<Double> & stabilParams, bool movingMesh, std::string stabilType)
: FluidMechInt(density, kinematicViscosity,movingMesh, stabilType), stabilParams_(stabilParams) {
	isSolDependent_ = true;
	isaxi_=false;
	name_ = "FluidMechPlaneStiffInt_PV";
}
//*******************************************************************************************************************************
//****************** mixed FEM **************************************************************************************************
//*******************************************************************************************************************************
FluidMechPlaneMixedStiffInt_UV::FluidMechPlaneMixedStiffInt_UV(Double density,Double kinematicViscosity,
		Matrix<Double> & stabilParams, bool movingMesh, std::string stabilType)
: FluidMechInt(density, kinematicViscosity,movingMesh, stabilType), stabilParams_(stabilParams){
	isSolDependent_ = true;
	isaxi_=false;
	name_ = "FluidMechPlaneMixedStiffInt_UV";
}
FluidMechPlaneMixedStiffInt_UQ::FluidMechPlaneMixedStiffInt_UQ(Double density,Double kinematicViscosity,
		Matrix<Double> & stabilParams, bool movingMesh, std::string stabilType)
: FluidMechInt(density,kinematicViscosity,movingMesh, stabilType), stabilParams_(stabilParams) {
	isSolDependent_ = true;
	isaxi_=false;
	name_ = "FluidMechPlaneMixedStiffInt_UQ";
}
FluidMechPlaneMixedStiffInt_PV::FluidMechPlaneMixedStiffInt_PV(Double density,Double kinematicViscosity,
		Matrix<Double> & stabilParams, bool movingMesh, std::string stabilType)
: FluidMechInt(density, kinematicViscosity,movingMesh, stabilType), stabilParams_(stabilParams) {
	isSolDependent_ = true;
	isaxi_=false;
	name_ = "FluidMechPlaneMixedStiffInt_PV";
}


//*******************************************************************************************************************************
//************************************************************************************************
//destructors 
//************************************************************************************************
//*******************************************************************************************************************************

//*******************************************************************************************************************************
//****************** stabilized FEM *********************************************************************************************
//*******************************************************************************************************************************

FluidMechPlaneStiffInt_UV::~FluidMechPlaneStiffInt_UV() {
}
FluidMechPlaneStiffInt_PQ::~FluidMechPlaneStiffInt_PQ() {
}
FluidMechPlaneStiffInt_UQ::~FluidMechPlaneStiffInt_UQ() {
}
FluidMechPlaneStiffInt_PV::~FluidMechPlaneStiffInt_PV() {
}

//*******************************************************************************************************************************
//****************** mixed FEM **************************************************************************************************
//*******************************************************************************************************************************
FluidMechPlaneMixedStiffInt_UV::~FluidMechPlaneMixedStiffInt_UV() {
}
FluidMechPlaneMixedStiffInt_UQ::~FluidMechPlaneMixedStiffInt_UQ() {
}
FluidMechPlaneMixedStiffInt_PV::~FluidMechPlaneMixedStiffInt_PV() {
}


//************************************************************************************************
//methods 
//************************************************************************************************

//**************************************************************************************************************************
//***************** stabilized FEM *****************************************************************************************
//**************************************************************************************************************************
void FluidMechPlaneStiffInt_UV::CalcElementMatrix( Matrix<Double>& elemMat,
		EntityIterator& ent1, EntityIterator& ent2 ) {

	// Extract pointer to reference element and get coordinates
	ExtractElemInfo( ent1 );
	UInt elemNumber = ent1.GetElem()->elemNum;

	Matrix<Double> elemResult_;
	sol_->GetElemSolutionAsMatrix( elemResult_, ent1 );

	Matrix<Double> gridElemResult_;
	if(movingMesh_){
		gridSol_->GetElemSolutionAsMatrix( gridElemResult_, ent1);
		elemResult_-=gridElemResult_;
	}

	ptelem->SetAnsatzFct( ansatzFct1_ );
	const UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
	const UInt nrIntPts= ptelem->GetNumIntPoints();

	const UInt spaceDim = ptelem->GetDim();  

	bool computeTaus;
	if(stabilParams_[elemNumber-1][3] < -0.001){
		computeTaus=true;
		//std::cerr << "compT\n";
	}
	else{
		computeTaus=false;
		//std::cerr << "useT\n";
	}

	Vector<Double>  Vx, Vy, Vz, VGx, VGy;  
	Double  VxAtIP, VyAtIP, VzAtIP, VL2, VL2AtIp, VMax;
	Double tau_m, tau_mu, tau_c;
	Double A_elem, h_k;

	Vx.Resize(nrFncs);//Vx.Init();
	Vy.Resize(nrFncs);//Vy.Init();
	//VGx.Resize(nrFncs);//VGx.Init();
	//VGy.Resize(nrFncs);//VGy.Init();
	if(spaceDim==3){
		Vz.Resize(nrFncs);//Vz.Init();
	}

	for (UInt i=0; i<nrFncs; i++) {
		Vx[i]=elemResult_[0][i];
		Vy[i]=elemResult_[1][i];

		//if(movingMesh_){
		//  VGx[i]=gridElemResult_[0][i];
		//  VGy[i]=gridElemResult_[1][i];
		//}
		if(spaceDim==3)
			Vz[i]=elemResult_[2][i];
	}

	const Vector<Double> & intWeights = ptelem->GetIntWeights();  
	Double jacDet, jacDet_intWeight, multAux;
	//UInt N=spaceDim;  // DOFs per Node (Ux, Uy)

	Double lambda_k;

	Matrix<Double> locElemMat;

	Matrix<Double> A1, A2, A3;
	Matrix<Double> A_a; 
	Matrix<Double> A_ba, A_ba1, A_ba2, A_ba3; 
	Matrix<Double> A_cdef; 
	Matrix<Double> A_ga, A_gb, A_gc, A_gd, A_ge, A_gf, A_gg, A_gh, A_gj;
	//Matrix<Double> A_k; //A_kAll;

	Matrix<Double> xiDxDy, xiDxxDyyDxy;
	Vector<Double> xiDx, xiDy, xiDz, xi;
	Vector<Double> xiDxx, xiDyy, xiDzz;
	Vector<Double> auxJ, auxI;

	if(computeTaus){
		A_elem = ptelem->CalcVolume(ptCoord_, isaxi_);
		//A_elem /= 4.0; //wg quadratische Element

		if(spaceDim==2){
			VL2 = sqrt( (Vx.NormL2()*Vx.NormL2()) + (Vy.NormL2()*Vy.NormL2()) );
			h_k = sqrt(A_elem  / (PI) );
		}
		else if(spaceDim==3){
			VL2 = sqrt( (Vx.NormL2()*Vx.NormL2()) + (Vy.NormL2()*Vy.NormL2()) + (Vz.NormL2()*Vz.NormL2()) );
			h_k = std::pow( (0.75*A_elem/(PI)) ,(1.0/3.0) );
		}
		else
			EXCEPTION("Wrong spaceDim!");

		VMax=abs(Vx[0]);
		if (abs(Vy[0])>VMax)
			VMax=abs(Vy[0]);
		if(spaceDim==3)
			if (abs(Vz[0])>VMax)
				VMax=abs(Vz[0]);

		if(spaceDim==2)
			for (UInt i=1; i<nrFncs; i++) {
				if (abs(Vx[i])>abs(Vy[i]) && abs(Vx[i])>VMax )
					VMax=abs(Vx[i]);
				else if (abs(Vx[i])<abs(Vy[i]) && abs(Vy[i])>VMax )
					VMax=abs(Vy[i]);
			}

		else if(spaceDim==3)
			for (UInt i=1; i<nrFncs; i++) {
				if (abs(Vx[i]) > abs(Vy[i]) && abs(Vx[i]) > abs(Vz[i]) && abs(Vx[i])>VMax )
					VMax=abs(Vx[i]);
				else if (abs(Vy[i])>abs(Vx[i]) && abs(Vy[i])>abs(Vz[i]) && abs(Vy[i])>VMax )
					VMax=abs(Vy[i]);
				else if (abs(Vz[i])>abs(Vx[i]) && abs(Vz[i])>abs(Vy[i]) && abs(Vz[i])>VMax )
					VMax=abs(Vz[i]);
			}
		else
			EXCEPTION("Wrong spaceDim!");

		VL2AtIp = -1.0;
		ComputeTaus(tau_m,tau_mu,tau_c,VL2,VL2AtIp,VMax,h_k,kinematicViscosity_,lambda_k,nrFncs);
		stabilParams_[elemNumber-1][0]=tau_m;
		stabilParams_[elemNumber-1][1]=tau_mu;
		stabilParams_[elemNumber-1][2]=tau_c;
		stabilParams_[elemNumber-1][3]=1.0;
	}
	else{
		tau_m  = stabilParams_[elemNumber-1][0];
		tau_mu = stabilParams_[elemNumber-1][1];
		tau_c  = stabilParams_[elemNumber-1][2];
		//      if (tau_c<1e-15 || tau_m<1e-15 || tau_mu<1e-15) {
		//        WARN("one stabilization parameter is smaller than 1e-15");
		//        std::cerr << "tau_c="  << tau_c << "\n" 
		//                  << "tau_m="  << tau_m << "\n" 
		//                  << "tau_mu=" << tau_mu << "\n";
		//      }
		//      if (tau_c>1e3 || tau_m>1e3 || tau_mu>1e3) {
		//        WARN("one stabilization parameter is larger than 1e3");
		//        std::cerr << "tau_c="  << tau_c << "\n" 
		//                  << "tau_m="  << tau_m << "\n" 
		//                  << "tau_mu=" << tau_mu << "\n";
		//      }
	}


	// set matrix to desired size and set all elements to zero
	locElemMat.Resize(nrFncs*spaceDim); locElemMat.Init();
	//A_kAll.Resize(nrFncs); A_kAll.Init();

//#pragma omp parallel for private(auxJ,auxI,VxAtIP,VyAtIP,VzAtIP,jacDet,jacDet_intWeight,multAux,xiDxDy,xiDxxDyyDxy,xiDxx,xiDyy,xiDzz,xiDx,xiDy,xiDz,xi,A1,A2,A3,A_a,A_ba,A_ba1,A_ba2,A_ba3,A_cdef,A_ga,A_gb,A_gc,A_gd,A_ge,A_gf,A_gg,A_gh,A_gj)
	for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
	{
		jacDet = 0;

		ptelem->GetShFncAtIp(xi, actIntPt, ent1.GetElem());
		ptelem->GetGlobDerivShFncAtIp(xiDxDy, actIntPt, ptCoord_, jacDet, ent1.GetElem());
		if( abs(convStabSign_)>1e-5 || abs(viscStabSign_)>1e-5)
			ptelem->GetGlob2ndDerivShFncAtIp(xiDxxDyyDxy, xiDxDy, actIntPt, ptCoord_, ent1.GetElem());

		jacDet_intWeight = intWeights[actIntPt-1] * jacDet;

		xiDx.Resize(nrFncs); //xiDx.Init();
		xiDy.Resize(nrFncs); //xiDy.Init();
		if(spaceDim==3)
			xiDz.Resize(nrFncs); //xiDz.Init();

		if( abs(convStabSign_)>1e-5 || abs(viscStabSign_)>1e-5){
			xiDxx.Resize(nrFncs); //xiDxx.Init();
			xiDyy.Resize(nrFncs); //xiDyy.Init();
			if(spaceDim==3)
				xiDzz.Resize(nrFncs); //xiDzz.Init();
		}

		for (UInt i=0; i<nrFncs; i++) {
			xiDx[i] = xiDxDy[i][0];
			xiDy[i] = xiDxDy[i][1];
			if(spaceDim==3)
				xiDz[i] = xiDxDy[i][2];

			if( abs(convStabSign_)>1e-5 || abs(viscStabSign_)>1e-5){
				xiDxx[i] = xiDxxDyyDxy[i][0];
				xiDyy[i] = xiDxxDyyDxy[i][1];
				if(spaceDim==3)
					xiDzz[i] = xiDxxDyyDxy[i][2];
			}
		}

		Vx.Inner(xi,VxAtIP);
		Vy.Inner(xi,VyAtIP);
		if(spaceDim==3)
			Vz.Inner(xi,VzAtIP);
		//VGx.Inner(xiDx,VGxDxAtIP);
		//VGy.Inner(xiDy,VGyDyAtIP);

		//************************************************
		//       Double VGxAtIP, VGyAtIP;
		//       Matrix<Double> GLC;
		//       VGx.Inner(xi,VGxAtIP);
		//       VGy.Inner(xi,VGyAtIP);
		//       auxI =(xiDx*VGxAtIP);
		//       auxI+=(xiDy*VGyAtIP);
		//       auxJ =(xi*jacDet_intWeight);
		//       GLC.DyadicMult(auxJ,auxI);
		//       std::cout << "\nGLC\n" << GLC << std::endl;
		//************************************************

		auxI =(xiDx*VxAtIP);
		auxI+=(xiDy*VyAtIP);
		if(spaceDim==3)
			auxI+=(xiDz*VzAtIP);

		auxJ =(xi*jacDet_intWeight);
		A_a.DyadicMult(auxJ,auxI);
		A1 = A_a; 

		// LaplViscInt
		auxI = (xiDx*jacDet_intWeight*kinematicViscosity_);
		A_ba1.DyadicMult(xiDx,auxI);
		auxI = (xiDy*jacDet_intWeight*kinematicViscosity_);
		A_ba2.DyadicMult(xiDy,auxI);
		if(spaceDim==3){
			auxI = (xiDz*jacDet_intWeight*kinematicViscosity_);
			A_ba3.DyadicMult(xiDz,auxI);
		}
		A_ba = A_ba1 + A_ba2;
		if(spaceDim==3)
			A_ba+=A_ba3;
		A1 += A_ba;

		//Stabilisation Integrator A_cdef
		if( abs(convStabSign_)>1e-5 || abs(viscStabSign_)>1e-5){
			auxI =(xiDx*VxAtIP);
			auxI+=(xiDy*VyAtIP);
			auxI+=(xiDxx*kinematicViscosity_*(-1.0));
			auxI+=(xiDyy*kinematicViscosity_*(-1.0));
			if(spaceDim==3){
				auxI+=(xiDz*VzAtIP);
				auxI+=(xiDzz*kinematicViscosity_*(-1.0));
			}
			auxI*=(tau_mu*jacDet_intWeight);
			auxJ.Resize(nrFncs); auxJ.Init();
			if( abs(convStabSign_)>1e-5){
				auxJ+=(xiDx*VxAtIP*convStabSign_);
				auxJ+=(xiDy*VyAtIP*convStabSign_);
				if(spaceDim==3)
					auxJ+=(xiDz*VzAtIP*convStabSign_);

			}
			if(abs(viscStabSign_)>1e-5){
				auxJ+=(xiDxx*kinematicViscosity_*viscStabSign_);
				auxJ+=(xiDyy*kinematicViscosity_*viscStabSign_);
				if(spaceDim==3)
					auxJ+=(xiDzz*kinematicViscosity_*viscStabSign_);
			}
			A_cdef.DyadicMult(auxJ,auxI);
			A1 += A_cdef;
		}

		//auxI =(xi*(VGxDxAtIP+VGyDyAtIP)*jacDet_intWeight);
		//A_k.DyadicMult(xi,auxI);
		//A1-=A_k;
		//A_kAll+=A_k;

		if(abs(contiStabSign_)>1e-5){
			multAux=tau_c*contiStabSign_*jacDet_intWeight;

			if(spaceDim==2){
				auxI=(xiDy*multAux);
				A_gb.DyadicMult(xiDy,auxI);
				A2 = A1+A_gb;

				auxI=(xiDx*multAux);
				A_ga.DyadicMult(xiDx,auxI);
				A1 += A_ga;

				auxI=xiDx*multAux;
				A_gc.DyadicMult(xiDy,auxI);

				auxI=xiDy*multAux;
				A_gd.DyadicMult(xiDx,auxI);
			}
			else if(spaceDim==3){

				auxI=(xiDy*multAux);
				A_gb.DyadicMult(xiDy,auxI);
				A2 = A1+A_gb;

				auxI=(xiDz*multAux);
				A_gc.DyadicMult(xiDz,auxI);
				A3 = A1+A_gc;

				auxI=(xiDx*multAux);
				A_ga.DyadicMult(xiDx,auxI);
				A1 += A_ga;

				auxI=xiDx*multAux;
				A_gd.DyadicMult(xiDy,auxI);

				auxI=xiDx*multAux;
				A_ge.DyadicMult(xiDz,auxI);

				auxI=xiDy*multAux;
				A_gf.DyadicMult(xiDz,auxI);

				auxI=xiDy*multAux;
				A_gg.DyadicMult(xiDx,auxI);

				auxI=xiDz*multAux;
				A_gh.DyadicMult(xiDx,auxI);

				auxI=xiDz*multAux;
				A_gj.DyadicMult(xiDy,auxI);
			}
			else
				EXCEPTION("Wrong spaceDim!");

		}
//#pragma omp critical
{
	locElemMat.AddSubMatrix(A1,  0,          0);
	locElemMat.AddSubMatrix(A2,  nrFncs,    nrFncs);
	if(spaceDim==3)
		locElemMat.AddSubMatrix(A3,  2*nrFncs,2*nrFncs);
	if(abs(contiStabSign_)>1e-5){
		if(spaceDim==2){
			locElemMat.AddSubMatrix(A_gc,  nrFncs,    0);
			locElemMat.AddSubMatrix(A_gd,  0      ,    nrFncs);
		}
		else if(spaceDim==3){
			locElemMat.AddSubMatrix(A_gd,  nrFncs,    0);
			locElemMat.AddSubMatrix(A_ge,  2*nrFncs,  0);
			locElemMat.AddSubMatrix(A_gf,  2*nrFncs,  nrFncs);
			locElemMat.AddSubMatrix(A_gg,  0,          nrFncs);
			locElemMat.AddSubMatrix(A_gh,  0,          2*nrFncs);
			locElemMat.AddSubMatrix(A_gj,  nrFncs,    2*nrFncs);
		}
		else 
			EXCEPTION("EXCEPTION!!");
	}
}
	}
	ResortElementMatrix(elemMat, locElemMat, nrFncs, spaceDim);
}

void FluidMechPlaneStiffInt_PQ::CalcElementMatrix( Matrix<Double>& elemMat,
		EntityIterator& ent1, 
		EntityIterator& ent2 ) {

	// Extract pointer to reference element and get coordinates
	ExtractElemInfo( ent1 );
	UInt elemNumber = ent1.GetElem()->elemNum;

	Matrix<Double> elemResult_;
	sol_->GetElemSolutionAsMatrix( elemResult_, ent1);

	if(movingMesh_){
		Matrix<Double> gridElemResult_;
		gridSol_->GetElemSolutionAsMatrix( gridElemResult_, ent1);
		elemResult_-=gridElemResult_;
	}
	ptelem->SetAnsatzFct( ansatzFct1_ );
	const UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
	const UInt nrIntPts= ptelem->GetNumIntPoints();

	const UInt spaceDim = ptelem->GetDim();  

	Vector<Double>  Vx, Vy, Vz;
	Double VL2, VL2AtIp, VMax;
	Double tau_m, tau_mu, tau_c;
	Double A_elem, h_k;

	bool computeTaus;
	if(stabilParams_[elemNumber-1][3] < -0.001){
		computeTaus=true;
		//if(elemNumber==25)
			//std::cerr << "compute";
	}
	else{
		computeTaus=false;
		//if(elemNumber==25)
		//std::cerr << "saved";
	}

	Vx.Resize(nrFncs);//Vx.Init();
	Vy.Resize(nrFncs);//Vy.Init();
	if(spaceDim==3)
		Vz.Resize(nrFncs);//Vz.Init();
	for (UInt i=0; i<nrFncs; i++) {
		Vx[i]=elemResult_[0][i];
		Vy[i]=elemResult_[1][i];
		if(spaceDim==3)
			Vz[i]=elemResult_[2][i];
	}

	const Vector<Double> & intWeights = ptelem->GetIntWeights();  
	Double jacDet, jacDet_intWeight;
	UInt N=1;  // pressure DOFs per Node

	Double lambda_k;

	Matrix<Double> C_a, C_a1, C_a2, C_a3;

	Matrix<Double> xiDxDy;
	Vector<Double> xiDx, xiDy, xiDz, xi;
	Vector<Double> auxI;

	if(computeTaus){
		A_elem = ptelem->CalcVolume(ptCoord_, isaxi_);
		//A_elem /= 4.0; //wg quadratische Element

		if(spaceDim==2){
			VL2 = sqrt( (Vx.NormL2()*Vx.NormL2()) + (Vy.NormL2()*Vy.NormL2()) );
			h_k = sqrt(A_elem  / (PI) );
		}
		else if(spaceDim==3){
			VL2 = sqrt( (Vx.NormL2()*Vx.NormL2()) + (Vy.NormL2()*Vy.NormL2()) + (Vz.NormL2()*Vz.NormL2()) );
			h_k = std::pow( (0.75*A_elem/(PI)) ,(1.0/3.0) );
		}
		else
			EXCEPTION("Wrong spaceDim!");

		VMax=abs(Vx[0]);
		if (abs(Vy[0])>VMax)
			VMax=abs(Vy[0]);
		if(spaceDim==3)
			if (abs(Vz[0])>VMax)
				VMax=abs(Vz[0]);

		if(spaceDim==2)
			for (UInt i=1; i<nrFncs; i++) {
				if (abs(Vx[i])>abs(Vy[i]) && abs(Vx[i])>VMax )
					VMax=abs(Vx[i]);
				else if (abs(Vx[i])<abs(Vy[i]) && abs(Vy[i])>VMax )
					VMax=abs(Vy[i]);
			}

		else if(spaceDim==3)
			for (UInt i=1; i<nrFncs; i++) {
				if (abs(Vx[i]) > abs(Vy[i]) && abs(Vx[i]) > abs(Vz[i]) && abs(Vx[i])>VMax )
					VMax=abs(Vx[i]);
				else if (abs(Vy[i])>abs(Vx[i]) && abs(Vy[i])>abs(Vz[i]) && abs(Vy[i])>VMax )
					VMax=abs(Vy[i]);
				else if (abs(Vz[i])>abs(Vx[i]) && abs(Vz[i])>abs(Vy[i]) && abs(Vz[i])>VMax )
					VMax=abs(Vz[i]);
			}
		else
			EXCEPTION("Wrong spaceDim!");

		VL2AtIp = -1.0;
		ComputeTaus(tau_m,tau_mu,tau_c,VL2,VL2AtIp,VMax,h_k,kinematicViscosity_,lambda_k,nrFncs);
		stabilParams_[elemNumber-1][0]=tau_m;
		stabilParams_[elemNumber-1][1]=tau_mu;
		stabilParams_[elemNumber-1][2]=tau_c;
		stabilParams_[elemNumber-1][3]=1.0;
	}
	else{
		tau_m  = stabilParams_[elemNumber-1][0];
		tau_mu = stabilParams_[elemNumber-1][1];
		tau_c  = stabilParams_[elemNumber-1][2];
		//      if (tau_c<1e-15 || tau_m<1e-15 || tau_mu<1e-15) {
		//        WARN("one stabilization parameter is smaller than 1e-15");
		//        std::cerr << "tau_c="  << tau_c << "\n" 
		//                  << "tau_m="  << tau_m << "\n" 
		//                  << "tau_mu=" << tau_mu << "\n";
		//      }
		//      if (tau_c>1e3 || tau_m>1e3 || tau_mu>1e3) {
		//        WARN("one stabilization parameter is larger than 1e3");
		//        std::cerr << "tau_c="  << tau_c << "\n" 
		//                  << "tau_m="  << tau_m << "\n" 
		//                  << "tau_mu=" << tau_mu << "\n";
		//      }
	}

	// set matrix to desired size and set all elements to zero
	elemMat.Resize(nrFncs*N); elemMat.Init();
	if(abs(presStabSign_)>1e-5){
		for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++){
			jacDet = 0;

			ptelem->GetShFncAtIp(xi, actIntPt, ent1.GetElem());
			ptelem->GetGlobDerivShFncAtIp(xiDxDy, actIntPt, ptCoord_, jacDet, ent1.GetElem());

			jacDet_intWeight = intWeights[actIntPt-1] * jacDet;

			xiDx.Resize(nrFncs); //xiDx.Init();
			xiDy.Resize(nrFncs); //xiDy.Init();
			if(spaceDim==3)
				xiDz.Resize(nrFncs); //xiDz.Init();

			for (UInt i=0; i<nrFncs; i++) {
				xiDx[i] = xiDxDy[i][0];
				xiDy[i] = xiDxDy[i][1];
				if(spaceDim==3)
					xiDz[i] = xiDxDy[i][2];
			}

			//********************************************
			// Stabilisation Integrators
			auxI=xiDx*tau_m*presStabSign_*jacDet_intWeight;
			C_a1.DyadicMult(xiDx,auxI);

			auxI=xiDy*tau_m*presStabSign_*jacDet_intWeight;
			C_a2.DyadicMult(xiDy,auxI);

			if(spaceDim==3){
				auxI=xiDz*tau_m*presStabSign_*jacDet_intWeight;
				C_a3.DyadicMult(xiDz,auxI);
			}

			if(spaceDim==2)
				C_a=C_a1+C_a2;
			else if(spaceDim==3)
				C_a=C_a1+C_a2+C_a3;
			else
				EXCEPTION("EXCEPTION!!");

			{
				elemMat.AddSubMatrix(C_a,0,0);
			}
		}
	}
}

void FluidMechPlaneStiffInt_UQ::CalcElementMatrix( Matrix<Double>& elemMat,
		EntityIterator& ent1, 
		EntityIterator& ent2 ) {
	//bool is3D=false;

	// Extract pointer to reference element and get coordinates
	ExtractElemInfo( ent1 );
	UInt elemNumber = ent1.GetElem()->elemNum;

	Matrix<Double> elemResult_;
	sol_->GetElemSolutionAsMatrix( elemResult_, ent1);

	if(movingMesh_){
		Matrix<Double> gridElemResult_;
		gridSol_->GetElemSolutionAsMatrix( gridElemResult_, ent1 );
		elemResult_-=gridElemResult_;
	}

	ptelem->SetAnsatzFct( ansatzFct1_ );
	const UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
	const UInt nrIntPts= ptelem->GetNumIntPoints();

	const UInt spaceDim = ptelem->GetDim();  

	Vector<Double>  Vx, Vy, Vz;  
	Double  VxAtIP, VyAtIP, VzAtIP, VL2, VL2AtIp, VMax;
	Double tau_m, tau_mu, tau_c;
	Double A_elem, h_k;
	bool computeTaus;

	if(stabilParams_[elemNumber-1][3] < -0.001){
		computeTaus=true;
	}
	else{
		computeTaus=false;
	}

	Vx.Resize(nrFncs);//Vx.Init();
	Vy.Resize(nrFncs);//Vy.Init();
	if(spaceDim==3) 
		Vz.Resize(nrFncs);//Vz.Init();

	for (UInt i=0; i<nrFncs; i++) {
		Vx[i]=elemResult_[0][i];
		Vy[i]=elemResult_[1][i];
		if(spaceDim==3)
			Vz[i]=elemResult_[2][i];
	}

	const Vector<Double> & intWeights = ptelem->GetIntWeights();  
	Double jacDet, jacDet_intWeight, multAux;
	UInt N=spaceDim;  // DOFs per Node (Ux, Uy)

	Double lambda_k;

	Matrix<Double> locElemMat;

	Matrix<Double> D1, D2, D3;
	Matrix<Double> D_aa, D_ab, D_ac;
	Matrix<Double> D_bc1, D_bc2, D_bc3;

	Matrix<Double> xiDxDy, xiDxxDyyDxy;
	Vector<Double> xiDx, xiDy, xiDz, xi;
	Vector<Double> xiDxx, xiDyy, xiDzz;
	Vector<Double> auxI;

	if(computeTaus){
		A_elem = ptelem->CalcVolume(ptCoord_, isaxi_);
		//A_elem /= 4.0; //wg quadratische Element

		if(spaceDim==2){
			VL2 = sqrt( (Vx.NormL2()*Vx.NormL2()) + (Vy.NormL2()*Vy.NormL2()) );
			h_k = sqrt(A_elem  / (PI) );
		}
		else if(spaceDim==3) {
			VL2 = sqrt( (Vx.NormL2()*Vx.NormL2()) + (Vy.NormL2()*Vy.NormL2()) + (Vz.NormL2()*Vz.NormL2()) );
			h_k = std::pow( (0.75*A_elem/(PI)) ,(1.0/3.0) );
		}
		else
			EXCEPTION("Wrong spaceDim!");

		VMax=abs(Vx[0]);
		if (abs(Vy[0])>VMax)
			VMax=abs(Vy[0]);
		if(spaceDim==3)
			if (abs(Vz[0])>VMax)
				VMax=abs(Vz[0]);

		if(spaceDim==2)
			for (UInt i=1; i<nrFncs; i++) {
				if (abs(Vx[i])>abs(Vy[i]) && abs(Vx[i])>VMax )
					VMax=abs(Vx[i]);
				else if (abs(Vx[i])<abs(Vy[i]) && abs(Vy[i])>VMax )
					VMax=abs(Vy[i]);
			}

		else if(spaceDim==3)
			for (UInt i=1; i<nrFncs; i++) {
				if (abs(Vx[i]) > abs(Vy[i]) && abs(Vx[i]) > abs(Vz[i]) && abs(Vx[i])>VMax )
					VMax=abs(Vx[i]);
				else if (abs(Vy[i])>abs(Vx[i]) && abs(Vy[i])>abs(Vz[i]) && abs(Vy[i])>VMax )
					VMax=abs(Vy[i]);
				else if (abs(Vz[i])>abs(Vx[i]) && abs(Vz[i])>abs(Vy[i]) && abs(Vz[i])>VMax )
					VMax=abs(Vz[i]);
			}
		else
			EXCEPTION("Wrong spaceDim!");

		VL2AtIp = -1.0;
		ComputeTaus(tau_m,tau_mu,tau_c,VL2,VL2AtIp,VMax,h_k,kinematicViscosity_,lambda_k,nrFncs);
		stabilParams_[elemNumber-1][0]=tau_m;
		stabilParams_[elemNumber-1][1]=tau_mu;
		stabilParams_[elemNumber-1][2]=tau_c;
		stabilParams_[elemNumber-1][3]=1.0;
	}
	else{
		tau_m  = stabilParams_[elemNumber-1][0];
		tau_mu = stabilParams_[elemNumber-1][1];
		tau_c  = stabilParams_[elemNumber-1][2];
		//      if (tau_c<1e-15 || tau_m<1e-15 || tau_mu<1e-15) {
		//        WARN("one stabilization parameter is smaller than 1e-15");
		//        std::cerr << "tau_c="  << tau_c << "\n" 
		//                  << "tau_m="  << tau_m << "\n" 
		//                  << "tau_mu=" << tau_mu << "\n";
		//      }
		//      if (tau_c>1e3 || tau_m>1e3 || tau_mu>1e3) {
		//        WARN("one stabilization parameter is larger than 1e3");
		//        std::cerr << "tau_c="  << tau_c << "\n" 
		//                  << "tau_m="  << tau_m << "\n" 
		//                  << "tau_mu=" << tau_mu << "\n";
		//      }
	}

	// set matrix to desired size and set all elements to zero
	locElemMat.Resize(nrFncs,nrFncs*N); locElemMat.Init();

//#pragma omp parallel for private(auxI,VxAtIP,VyAtIP,VzAtIP,jacDet,jacDet_intWeight,multAux,xiDxDy,xiDxxDyyDxy,xiDxx,xiDyy,xiDzz,xiDx,xiDy,xiDz,xi,D1,D2,D3,D_aa,D_ab,D_ac,D_bc1,D_bc2,D_bc3)
	for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
	{
		jacDet = 0;

		ptelem->GetShFncAtIp(xi, actIntPt, ent1.GetElem());
		ptelem->GetGlobDerivShFncAtIp(xiDxDy, actIntPt, ptCoord_, jacDet, ent1.GetElem());
		if(abs(presStabSign_)>1e-5)
			ptelem->GetGlob2ndDerivShFncAtIp(xiDxxDyyDxy, xiDxDy, actIntPt, ptCoord_, ent1.GetElem());

		jacDet_intWeight = intWeights[actIntPt-1] * jacDet;

		xiDx.Resize(nrFncs); //xiDx.Init();
		xiDy.Resize(nrFncs); //xiDy.Init();
		if(spaceDim==3)
			xiDz.Resize(nrFncs); //xiDz.Init();


		if(abs(presStabSign_)>1e-5){
			xiDxx.Resize(nrFncs); //xiDxx.Init();
			xiDyy.Resize(nrFncs); //xiDyy.Init();
			if(spaceDim==3)
				xiDzz.Resize(nrFncs); //xiDzz.Init();
		}

		for (UInt i=0; i<nrFncs; i++) {
			xiDx[i] = xiDxDy[i][0];
			xiDy[i] = xiDxDy[i][1];
			if(spaceDim==3)
				xiDz[i] = xiDxDy[i][2];

			if(abs(presStabSign_)>1e-5){
				xiDxx[i] = xiDxxDyyDxy[i][0];
				xiDyy[i] = xiDxxDyyDxy[i][1];
				if(spaceDim==3)
					xiDzz[i] = xiDxxDyyDxy[i][2];
			}
		}

		Vx.Inner(xi,VxAtIP);
		Vy.Inner(xi,VyAtIP);
		if(spaceDim==3)
			Vz.Inner(xi,VzAtIP);

		// Continuity Integrator
		multAux=contiTermSign_*jacDet_intWeight;
		auxI=xiDx*multAux;
		D_aa.DyadicMult(xi,auxI);
		D1 = D_aa;

		auxI=xiDy*multAux;
		D_ab.DyadicMult(xi,auxI);
		D2 = D_ab;

		if(spaceDim==3){
			auxI=xiDz*multAux;
			D_ac.DyadicMult(xi,auxI);
			D3 = D_ac;
		}

		//Stabilization Integrators
		if(abs(presStabSign_)>1e-5){
			auxI =(xiDx*VxAtIP);
			auxI+=(xiDy*VyAtIP);
			if(spaceDim==3)
				auxI+=(xiDz*VzAtIP);
			auxI+=(xiDxx*(-1.0)*kinematicViscosity_);
			auxI+=(xiDyy*(-1.0)*kinematicViscosity_);
			if(spaceDim==3)
				auxI+=(xiDzz*(-1.0)*kinematicViscosity_);
			auxI*=(presStabSign_*tau_m*jacDet_intWeight);
			D_bc1.DyadicMult(xiDx,auxI);
			D_bc2.DyadicMult(xiDy,auxI);
			if(spaceDim==3)
				D_bc3.DyadicMult(xiDz,auxI);
			D1 += D_bc1;
			D2 += D_bc2;
			if(spaceDim==3)
				D3 += D_bc3;
		}

//#pragma omp critical
		{
			locElemMat.AddSubMatrix(D1,  0,0);
			locElemMat.AddSubMatrix(D2,  0,nrFncs);
			if(spaceDim==3)
				locElemMat.AddSubMatrix(D3,  0,2*nrFncs);
		}
	}
	ColResortElementMatrix(elemMat, locElemMat, nrFncs, nrFncs, 1,N);
}

void FluidMechPlaneStiffInt_PV::CalcElementMatrix( Matrix<Double>& elemMat,
		EntityIterator& ent1, 
		EntityIterator& ent2 ) {

	// Extract pointer to reference element and get coordinates
	ExtractElemInfo( ent1 );
	UInt elemNumber = ent1.GetElem()->elemNum;

	Matrix<Double> elemResult_;
	sol_->GetElemSolutionAsMatrix( elemResult_, ent1 );

	if(movingMesh_){
		Matrix<Double> gridElemResult_;
		gridSol_->GetElemSolutionAsMatrix( gridElemResult_, ent1);
		elemResult_-=gridElemResult_;
	}

	ptelem->SetAnsatzFct( ansatzFct1_ );
	const UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
	const UInt nrIntPts= ptelem->GetNumIntPoints();

	const UInt spaceDim = ptelem->GetDim();  

	Vector<Double>  Vx, Vy, Vz;  
	Double  VxAtIP, VyAtIP, VzAtIP, VL2, VL2AtIp, VMax;
	Double tau_m, tau_mu, tau_c;
	Double A_elem, h_k;

	bool computeTaus;
	if(stabilParams_[elemNumber-1][3] < -0.001){
		computeTaus=true;
	}
	else{
		computeTaus=false;
	}

	Vx.Resize(nrFncs);//Vx.Init();
	Vy.Resize(nrFncs);//Vy.Init();
	if(spaceDim==3)
		Vz.Resize(nrFncs);//Vz.Init();
	for (UInt i=0; i<nrFncs; i++) {
		Vx[i]=elemResult_[0][i];
		Vy[i]=elemResult_[1][i];
		if(spaceDim==3)
			Vz[i]=elemResult_[2][i];
	}

	const Vector<Double> & intWeights = ptelem->GetIntWeights();  
	Double jacDet, jacDet_intWeight, multAux;
	UInt N=spaceDim;  // DOFs per Node (Ux, Uy)

	Double lambda_k;

	Matrix<Double> locElemMat;

	Matrix<Double> B1, B2, B3;
	Matrix<Double> B_aa, B_ab, B_ac;
	Matrix<Double> B_bc1, B_bc2, B_bc3;

	Matrix<Double> xiDxDy, xiDxxDyyDxy;
	Vector<Double> xiDx, xiDy, xiDz, xi;
	Vector<Double> xiDxx, xiDyy, xiDzz;
	Vector<Double> auxI, auxJ;

	if(computeTaus){
		A_elem = ptelem->CalcVolume(ptCoord_, isaxi_);
		//A_elem /= 4.0; //wg quadratische Element

		if(spaceDim==2){
			VL2 = sqrt( (Vx.NormL2()*Vx.NormL2()) + (Vy.NormL2()*Vy.NormL2()) );
			h_k = sqrt(A_elem  / (PI) );
		}
		else if(spaceDim==3){
			VL2 = sqrt( (Vx.NormL2()*Vx.NormL2()) + (Vy.NormL2()*Vy.NormL2()) + (Vz.NormL2()*Vz.NormL2()) );
			h_k = std::pow( (0.75*A_elem/(PI)) ,(1.0/3.0) );
		}
		else
			EXCEPTION("Wrong spaceDim!");

		VMax=abs(Vx[0]);
		if (abs(Vy[0])>VMax)
			VMax=abs(Vy[0]);
		if(spaceDim==3)
			if (abs(Vz[0])>VMax)
				VMax=abs(Vz[0]);

		if(spaceDim==2)
			for (UInt i=1; i<nrFncs; i++) {
				if (abs(Vx[i])>abs(Vy[i]) && abs(Vx[i])>VMax )
					VMax=abs(Vx[i]);
				else if (abs(Vx[i])<abs(Vy[i]) && abs(Vy[i])>VMax )
					VMax=abs(Vy[i]);
			}

		else if(spaceDim==3)
			for (UInt i=1; i<nrFncs; i++) {
				if (abs(Vx[i]) > abs(Vy[i]) && abs(Vx[i]) > abs(Vz[i]) && abs(Vx[i])>VMax )
					VMax=abs(Vx[i]);
				else if (abs(Vy[i])>abs(Vx[i]) && abs(Vy[i])>abs(Vz[i]) && abs(Vy[i])>VMax )
					VMax=abs(Vy[i]);
				else if (abs(Vz[i])>abs(Vx[i]) && abs(Vz[i])>abs(Vy[i]) && abs(Vz[i])>VMax )
					VMax=abs(Vz[i]);
			}
		else
			EXCEPTION("Wrong spaceDim!");

		VL2AtIp = -1.0;
		ComputeTaus(tau_m,tau_mu,tau_c,VL2,VL2AtIp,VMax,h_k,kinematicViscosity_,lambda_k,nrFncs);
		stabilParams_[elemNumber-1][0]=tau_m;
		stabilParams_[elemNumber-1][1]=tau_mu;
		stabilParams_[elemNumber-1][2]=tau_c;
		stabilParams_[elemNumber-1][3]=1.0;
	}
	else{
		tau_m  = stabilParams_[elemNumber-1][0];
		tau_mu = stabilParams_[elemNumber-1][1];
		tau_c  = stabilParams_[elemNumber-1][2];
		//      if (tau_c<1e-15 || tau_m<1e-15 || tau_mu<1e-15) {
		//        WARN("one stabilization parameter is smaller than 1e-15");
		//        std::cerr << "tau_c="  << tau_c << "\n" 
		//                  << "tau_m="  << tau_m << "\n" 
		//                  << "tau_mu=" << tau_mu << "\n";
		//      }
		//      if (tau_c>1e3 || tau_m>1e3 || tau_mu>1e3) {
		//        WARN("one stabilization parameter is larger than 1e3");
		//        std::cerr << "tau_c="  << tau_c << "\n" 
		//                  << "tau_m="  << tau_m << "\n" 
		//                  << "tau_mu=" << tau_mu << "\n";
		//      }
	}


	// set matrix to desired size and set all elements to zero
	//elemMat.Resize(nrFncs*N,nrFncs); elemMat.Init();
	locElemMat.Resize(nrFncs*N,nrFncs); locElemMat.Init();

//#pragma omp parallel for private(auxI,auxJ,VxAtIP,VyAtIP,VzAtIP,jacDet,jacDet_intWeight,multAux,xiDxDy,xiDxxDyyDxy,xiDxx,xiDyy,xiDzz,xiDx,xiDy,xiDz,xi,B1,B2,B3,B_aa,B_ab,B_ac,B_bc1,B_bc2,B_bc3)
	for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
	{
		jacDet = 0;

		ptelem->GetShFncAtIp(xi, actIntPt, ent1.GetElem());
		ptelem->GetGlobDerivShFncAtIp(xiDxDy, actIntPt, ptCoord_, jacDet, ent1.GetElem());
		ptelem->GetGlob2ndDerivShFncAtIp(xiDxxDyyDxy, xiDxDy, actIntPt, ptCoord_, ent1.GetElem());

		jacDet_intWeight = intWeights[actIntPt-1] * jacDet;

		xiDx.Resize(nrFncs); //xiDx.Init();
		xiDy.Resize(nrFncs); //xiDy.Init();
		if(spaceDim==3)
			xiDz.Resize(nrFncs); //xiDz.Init();

		xiDxx.Resize(nrFncs); //xiDxx.Init();
		xiDyy.Resize(nrFncs); //xiDyy.Init();
		if(spaceDim==3)
			xiDzz.Resize(nrFncs); //xiDzz.Init();

		for (UInt i=0; i<nrFncs; i++) {
			xiDx[i] = xiDxDy[i][0];
			xiDy[i] = xiDxDy[i][1];
			if(spaceDim==3)
				xiDz[i] = xiDxDy[i][2];

			xiDxx[i] = xiDxxDyyDxy[i][0];
			xiDyy[i] = xiDxxDyyDxy[i][1];
			if(spaceDim==3)
				xiDzz[i] = xiDxxDyyDxy[i][2];
		}

		Vx.Inner(xi,VxAtIP);
		Vy.Inner(xi,VyAtIP);
		if(spaceDim==3)
			Vz.Inner(xi,VzAtIP);


		// Pressure Integrator
		multAux=presTermSign_*jacDet_intWeight;
		auxI=(xi*multAux);
		B_aa.DyadicMult(xiDx,auxI);
		B1 = B_aa;

		B_ab.DyadicMult(xiDy,auxI);
		B2 = B_ab;

		if(spaceDim==3){
			B_ac.DyadicMult(xiDz,auxI);
			B3 = B_ac;
		}

		//Stabilisation Integrators
		if(abs(convStabSign_)>1e-5 || abs(viscStabSign_)>1e-5){
			auxJ.Resize(nrFncs); auxJ.Init();
			if(abs(convStabSign_)>1e-5){
				auxJ =(xiDx*VxAtIP*convStabSign_*tau_mu);
				auxJ+=(xiDy*VyAtIP*convStabSign_*tau_mu);
				if(spaceDim==3)
					auxJ+=(xiDz*VzAtIP*convStabSign_*tau_mu);
			}
			if(abs(viscStabSign_)>1e-5){
				auxJ+=(xiDxx*kinematicViscosity_*viscStabSign_*tau_m);
				auxJ+=(xiDyy*kinematicViscosity_*viscStabSign_*tau_m);
				if(spaceDim==3)
					auxJ+=(xiDzz*kinematicViscosity_*viscStabSign_*tau_m);
			}
			auxJ*=(jacDet_intWeight);

			B_bc1.DyadicMult(auxJ,xiDx);
			B_bc2.DyadicMult(auxJ,xiDy);
			if(spaceDim==3)
				B_bc3.DyadicMult(auxJ,xiDz);

			B1 += B_bc1;
			B2 += B_bc2;
			if(spaceDim==3)
				B3 += B_bc3;

		}

//#pragma omp critical
		{
			locElemMat.AddSubMatrix(B1,  0,     0);
			locElemMat.AddSubMatrix(B2,  nrFncs,0);
			if(spaceDim==3)
				locElemMat.AddSubMatrix(B3,2*nrFncs,0);
		}
	}
	RowResortElementMatrix(elemMat, locElemMat, nrFncs, nrFncs, N,1);
}

//**************************************************************************************************************************
//***************** mixed FEM **********************************************************************************************
//**************************************************************************************************************************
void FluidMechPlaneMixedStiffInt_UV::CalcElementMatrix( Matrix<Double>& elemMat,
		EntityIterator& ent1, EntityIterator& ent2 ) {

	// Extract pointer to reference element and get coordinates
	ExtractElemInfo( ent1 );

	Matrix<Double> elemResult_;
	sol_->GetElemSolutionAsMatrix( elemResult_, ent1 );

	Matrix<Double> gridElemResult_;
	if(movingMesh_){
		gridSol_->GetElemSolutionAsMatrix( gridElemResult_, ent1);
		elemResult_-=gridElemResult_;
	}

	ptelem->SetAnsatzFct( ansatzFct1_ );
	const UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
	const UInt nrIntPts= ptelem->GetNumIntPoints();

	const UInt spaceDim = ptelem->GetDim();  

	Vector<Double>  Vx, Vy, Vz, VGx, VGy;  
	Double  VxAtIP, VyAtIP, VzAtIP;

	Vx.Resize(nrFncs);//Vx.Init();
	Vy.Resize(nrFncs);//Vy.Init();
	//VGx.Resize(nrFncs);//VGx.Init();
	//VGy.Resize(nrFncs);//VGy.Init();
	if(spaceDim==3){
		Vz.Resize(nrFncs);//Vz.Init();
	}

	for (UInt i=0; i<nrFncs; i++) {
		Vx[i]=elemResult_[0][i];
		Vy[i]=elemResult_[1][i];

		//if(movingMesh_){
		//  VGx[i]=gridElemResult_[0][i];
		//  VGy[i]=gridElemResult_[1][i];
		//}
		if(spaceDim==3)
			Vz[i]=elemResult_[2][i];
	}

	const Vector<Double> & intWeights = ptelem->GetIntWeights();  
	Double jacDet, jacDet_intWeight;

	Matrix<Double> locElemMat;

	Matrix<Double> A1, A2, A3;
	Matrix<Double> A_a; 
	Matrix<Double> A_ba, A_ba1, A_ba2, A_ba3; 

	Matrix<Double> xiDxDy;
	Vector<Double> xiDx, xiDy, xiDz, xi;
	Vector<Double> auxJ, auxI;

	// set matrix to desired size and set all elements to zero
	locElemMat.Resize(nrFncs*spaceDim); locElemMat.Init();
	//A_kAll.Resize(nrFncs); A_kAll.Init();

	for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
	{
		jacDet = 0;

		ptelem->GetShFncAtIp(xi, actIntPt, ent1.GetElem());
		ptelem->GetGlobDerivShFncAtIp(xiDxDy, actIntPt, ptCoord_, jacDet, ent1.GetElem());

		jacDet_intWeight = intWeights[actIntPt-1] * jacDet;

		xiDx.Resize(nrFncs); //xiDx.Init();
		xiDy.Resize(nrFncs); //xiDy.Init();
		if(spaceDim==3)
			xiDz.Resize(nrFncs); //xiDz.Init();

		for (UInt i=0; i<nrFncs; i++) {
			xiDx[i] = xiDxDy[i][0];
			xiDy[i] = xiDxDy[i][1];
			if(spaceDim==3)
				xiDz[i] = xiDxDy[i][2];
		}

		Vx.Inner(xi,VxAtIP);
		Vy.Inner(xi,VyAtIP);
		if(spaceDim==3)
			Vz.Inner(xi,VzAtIP);

		auxI =(xiDx*VxAtIP);
		auxI+=(xiDy*VyAtIP);
		if(spaceDim==3)
			auxI+=(xiDz*VzAtIP);

		auxJ =(xi*jacDet_intWeight);
		A_a.DyadicMult(auxJ,auxI);
		A1 = A_a; 

		// LaplViscInt
		auxI = (xiDx*jacDet_intWeight*kinematicViscosity_);
		A_ba1.DyadicMult(xiDx,auxI);
		auxI = (xiDy*jacDet_intWeight*kinematicViscosity_);
		A_ba2.DyadicMult(xiDy,auxI);
		if(spaceDim==3){
			auxI = (xiDz*jacDet_intWeight*kinematicViscosity_);
			A_ba3.DyadicMult(xiDz,auxI);
		}
		A_ba = A_ba1 + A_ba2;
		if(spaceDim==3)
			A_ba+=A_ba3;
		A1 += A_ba;

		locElemMat.AddSubMatrix(A1,  0,          0);
		locElemMat.AddSubMatrix(A1,  nrFncs,    nrFncs);
		if(spaceDim==3)
			locElemMat.AddSubMatrix(A1,  2*nrFncs,2*nrFncs);
	}
	ResortElementMatrix(elemMat, locElemMat, nrFncs, spaceDim);
}

void FluidMechPlaneMixedStiffInt_UQ::CalcElementMatrix( Matrix<Double>& elemMat,
		EntityIterator& ent1,EntityIterator& ent2 ) {
	//bool is3D=false;

	// Extract pointer to reference element and get coordinates
	ExtractElemInfo( ent1 );

	ptelem->SetAnsatzFct( ansatzFct2_ );
	const UInt nrFncsV  = ptelem->GetNumFncs( ansatzFct2_ );
	const UInt nrIntPts= ptelem->GetNumIntPoints();

	ptelem->SetAnsatzFct( ansatzFct1_,false );
	const UInt nrFncsP  = ptelem->GetNumFncs( ansatzFct1_ );

	
	const UInt spaceDim = ptelem->GetDim();  

	const Vector<Double> & intWeights = ptelem->GetIntWeights();  
	Double jacDet, jacDet_intWeight, multAux;

	Matrix<Double> locElemMat;

	Matrix<Double> D1, D2, D3;
	Matrix<Double> D_aa, D_ab, D_ac;

	Matrix<Double> xiDxDy;
	Vector<Double> xiDx, xiDy, xiDz, xi;
	Vector<Double> auxI;

	// set matrix to desired size and set all elements to zero
	locElemMat.Resize(nrFncsP,nrFncsV*spaceDim); locElemMat.Init();

	for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
	{
		jacDet = 0;

		ptelem->SetAnsatzFct( ansatzFct2_ ,false);
		ptelem->GetGlobDerivShFncAtIp(xiDxDy, actIntPt, ptCoord_, jacDet, ent1.GetElem());

		ptelem->SetAnsatzFct( ansatzFct1_, false );
		ptelem->GetShFncAtIp(xi, actIntPt, ent1.GetElem());

		jacDet_intWeight = intWeights[actIntPt-1] * jacDet;

		xiDx.Resize(nrFncsV); //xiDx.Init();
		xiDy.Resize(nrFncsV); //xiDy.Init();
		if(spaceDim==3)
			xiDz.Resize(nrFncsV); //xiDz.Init();


		for (UInt i=0; i<nrFncsV; i++) {
			xiDx[i] = xiDxDy[i][0];
			xiDy[i] = xiDxDy[i][1];
			if(spaceDim==3)
				xiDz[i] = xiDxDy[i][2];
		}

		// Continuity Integrator
		multAux=contiTermSign_*jacDet_intWeight;
		auxI=xiDx*multAux;
		D_aa.DyadicMult(xi,auxI);
		D1 = D_aa;

		auxI=xiDy*multAux;
		D_ab.DyadicMult(xi,auxI);
		D2 = D_ab;

		if(spaceDim==3){
			auxI=xiDz*multAux;
			D_ac.DyadicMult(xi,auxI);
			D3 = D_ac;
		}

		locElemMat.AddSubMatrix(D1,  0,0);
		locElemMat.AddSubMatrix(D2,  0,nrFncsV);
		if(spaceDim==3)
			locElemMat.AddSubMatrix(D3,  0,2*nrFncsV);
	}
	
	ColResortElementMatrix(elemMat, locElemMat, nrFncsP, nrFncsV, 1,spaceDim);
}

void FluidMechPlaneMixedStiffInt_PV::CalcElementMatrix( Matrix<Double>& elemMat,
		EntityIterator& ent1,EntityIterator& ent2 ) {

	// Extract pointer to reference element and get coordinates
	ExtractElemInfo( ent1 );

	ptelem->SetAnsatzFct( ansatzFct1_ );
	const UInt nrFncsV = ptelem->GetNumFncs( ansatzFct1_ );
	const UInt nrIntPts= ptelem->GetNumIntPoints();

	ptelem->SetAnsatzFct( ansatzFct2_,false );
	const UInt nrFncsP = ptelem->GetNumFncs( ansatzFct2_ );

	const UInt spaceDim = ptelem->GetDim();  

	const Vector<Double> & intWeights = ptelem->GetIntWeights();  
	Double jacDet, jacDet_intWeight, multAux;

	Matrix<Double> locElemMat;

	Matrix<Double> B1, B2, B3;
	Matrix<Double> B_aa, B_ab, B_ac;

	Matrix<Double> xiDxDy;
	Vector<Double> xiDx, xiDy, xiDz, xi;
	Vector<Double> auxI, auxJ;

	// set matrix to desired size and set all elements to zero
	//elemMat.Resize(nrFncs*N,nrFncs); elemMat.Init();
	locElemMat.Resize(nrFncsV*spaceDim,nrFncsP); locElemMat.Init();

	for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
		jacDet = 0;

		ptelem->SetAnsatzFct( ansatzFct1_,false );
		ptelem->GetGlobDerivShFncAtIp(xiDxDy, actIntPt, ptCoord_, jacDet, ent1.GetElem());

		ptelem->SetAnsatzFct( ansatzFct2_,false );
		ptelem->GetShFncAtIp(xi, actIntPt, ent1.GetElem());

		
		jacDet_intWeight = intWeights[actIntPt-1] * jacDet;

		xiDx.Resize(nrFncsV); //xiDx.Init();
		xiDy.Resize(nrFncsV); //xiDy.Init();
		if(spaceDim==3)
			xiDz.Resize(nrFncsV); //xiDz.Init();


		for (UInt i=0; i<nrFncsV; i++) {
			xiDx[i] = xiDxDy[i][0];
			xiDy[i] = xiDxDy[i][1];
			if(spaceDim==3)
				xiDz[i] = xiDxDy[i][2];
		}


		// Pressure Integrator
		multAux=presTermSign_*jacDet_intWeight;
		auxI=(xi*multAux);
		B_aa.DyadicMult(xiDx,auxI);
		B1 = B_aa;

		B_ab.DyadicMult(xiDy,auxI);
		B2 = B_ab;

		if(spaceDim==3){
			B_ac.DyadicMult(xiDz,auxI);
			B3 = B_ac;
		}

		locElemMat.AddSubMatrix(B1,  0,     0);
		locElemMat.AddSubMatrix(B2,  nrFncsV,0);

		if(spaceDim==3)
			locElemMat.AddSubMatrix(B3,2*nrFncsV,0);
	}
	RowResortElementMatrix(elemMat, locElemMat, nrFncsV, nrFncsP,spaceDim,1);
}

} // end namespace CoupledField
