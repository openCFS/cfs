#include <iostream>
#include <fstream>

#include "fluidMechMassInt.hh"
#include "Utils/nodestoresol.hh"
#include "DataInOut/Logging/cfslog.hh"

namespace CoupledField
{
// declare logging stream
//  DECLARE_LOG(fluidMechMass)
//  DEFINE_LOG(fluidMechMass, "fluidMechMass")

//************************************************************************************************
//constructors 
//************************************************************************************************
FluidMechPlaneMassInt_UV::FluidMechPlaneMassInt_UV(Double density,Double kinematicViscosity,
		Matrix<Double> & stabilParams, bool movingMesh, std::string stabilType)
: FluidMechInt(density,kinematicViscosity,movingMesh, stabilType), stabilParams_(stabilParams) {
	isSolDependent_ = true;
	isaxi_=false;
	name_ = "FluidMechPlaneMassInt_UV";
	baseType_ = MASS;
}
FluidMechPlaneMassInt_UQ::FluidMechPlaneMassInt_UQ(Double density,Double kinematicViscosity,
		Matrix<Double> & stabilParams, bool movingMesh, std::string stabilType)
: FluidMechInt(density,kinematicViscosity, movingMesh, stabilType), stabilParams_(stabilParams) {
	isSolDependent_ = true;
	isaxi_=false;
	name_ = "FluidMechPlaneMassInt_UQ";
	baseType_ = MASS;
}
FluidMechPlaneMixedMassInt_UV::FluidMechPlaneMixedMassInt_UV(Double density,Double kinematicViscosity,
		Matrix<Double> & stabilParams, bool movingMesh, std::string stabilType)
: FluidMechInt(density,kinematicViscosity,movingMesh, stabilType), stabilParams_(stabilParams) {
	isSolDependent_ = true;
	isaxi_=false;
	name_ = "FluidMechPlaneMixedMassInt_UV";
	baseType_ = MASS;
}

//************************************************************************************************
//constructors 
//************************************************************************************************
FluidMechPlaneMassInt_UV::~FluidMechPlaneMassInt_UV() {
}
FluidMechPlaneMassInt_UQ::~FluidMechPlaneMassInt_UQ() {
}
FluidMechPlaneMixedMassInt_UV::~FluidMechPlaneMixedMassInt_UV() {
}


//************************************************************************************************
//methods 
//************************************************************************************************


//*********************************************************************************************
//********* stabilized FEM ********************************************************************
//*********************************************************************************************
void FluidMechPlaneMassInt_UV::CalcElementMatrix( Matrix<Double>& elemMat,
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

	bool computeTaus;
	if(stabilParams_[elemNumber-1][3] < -0.001){
		computeTaus=true;
	}
	else{
		computeTaus=false;
	}

	Vector<Double>  Vx, Vy, Vz;  
	Double  VxAtIP, VyAtIP, VzAtIP, VL2, VL2AtIp, VMax;
	Double tau_m, tau_mu, tau_c;
	Double lambda_k, A_elem, h_k;

	Vx.Resize(nrFncs);//Vx.Init(0.0);
	Vy.Resize(nrFncs);//Vy.Init(0.0);
	if(spaceDim==3)
		Vz.Resize(nrFncs);//Vz.Init(0.0);

	for (UInt i=0; i<nrFncs; i++) {
		Vx[i]=elemResult_[0][i];
		Vy[i]=elemResult_[1][i];
		if(spaceDim==3)
			Vz[i]=elemResult_[2][i];
	}

	const Vector<Double> & intWeights = ptelem->GetIntWeights();  

	Double jacDet;
	UInt N=spaceDim;  // DOFs per Node (Ux,Uy)

	// derivation of shape functions after global coordinates 
	Matrix<Double> locElemMat;

	Matrix<Double> A1;
	Matrix<Double> A_a; 

	Matrix<Double> xiDxDy;
	Matrix<Double> xiDxxDyyDxy;

	Vector<Double>  xiDxx, xiDyy, xiDzz;
	Vector<Double>  xiDx, xiDy, xiDz;
	Vector<Double>  xi;
	Vector<Double>  auxJ;

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
			Error("Wrong spaceDim!",__FILE__,__LINE__);

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
			Error("Wrong spaceDim!",__FILE__,__LINE__);



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
		//        Warning("one stabilization parameter is smaller than 1e-15",__FILE__,__LINE__);
		//        std::cerr << "tau_c="  << tau_c << "\n" 
		//                  << "tau_m="  << tau_m << "\n" 
		//                  << "tau_mu=" << tau_mu << "\n";
		//      }
		//      if (tau_c>1e3 || tau_m>1e3 || tau_mu>1e3) {
		//        Warning("one stabilization parameter is larger than 1e3",__FILE__,__LINE__);
		//        std::cerr << "tau_c="  << tau_c << "\n" 
		//                  << "tau_m="  << tau_m << "\n" 
		//                  << "tau_mu=" << tau_mu << "\n";
		//      }
	}

	// set matrix to desired size and set all elements to zero
	//elemMat.Resize(nrFncs*N); elemMat.Init(0.0);
	locElemMat.Resize(nrFncs*N); locElemMat.Init(0.0);

//#pragma omp parallel for private(auxJ,VxAtIP,VyAtIP,VzAtIP,jacDet,xiDxDy,xiDxxDyyDxy,xiDxx,xiDyy,xiDzz,xiDx,xiDy,xiDz,xi,A1)
	for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
	{
		jacDet = 0;

		ptelem->GetShFncAtIp(xi, actIntPt, ent1.GetElem());
		ptelem->GetGlobDerivShFncAtIp(xiDxDy, actIntPt, ptCoord_, jacDet, ent1.GetElem());
		if(abs(viscStabSign_) > 1e-5)
			ptelem->GetGlob2ndDerivShFncAtIp(xiDxxDyyDxy, xiDxDy, actIntPt, ptCoord_, ent1.GetElem());

		xiDx.Resize(nrFncs); //xiDx.Init(0.0);
		xiDy.Resize(nrFncs); //xiDy.Init(0.0);
		if(spaceDim==3)
			xiDz.Resize(nrFncs); //xiDz.Init(0.0);

		if(abs(viscStabSign_) > 1e-5){
			xiDxx.Resize(nrFncs); //xiDxx.Init(0.0);
			xiDyy.Resize(nrFncs); //xiDyy.Init(0.0);
			if(spaceDim==3)
				xiDzz.Resize(nrFncs); //xiDzz.Init(0.0);
		}

		for (UInt i=0; i<nrFncs; i++) {
			xiDx[i] = xiDxDy[i][0];
			xiDy[i] = xiDxDy[i][1];
			if(spaceDim==3)
				xiDz[i] = xiDxDy[i][2];

			if(abs(viscStabSign_) > 1e-5){
				xiDxx[i] = xiDxxDyyDxy[i][0];
				xiDyy[i] = xiDxxDyyDxy[i][1];
				if(spaceDim==3)
					xiDzz[i] = xiDxxDyyDxy[i][3];
			}
		}

		Vx.Inner(xi,VxAtIP);
		Vy.Inner(xi,VyAtIP);
		if(spaceDim==3)
			Vz.Inner(xi,VzAtIP);


		// Mass + Stab Mass Integrator
		auxJ=(xi);
		if(abs(convStabSign_) > 1e-5){
			auxJ+=(xiDx*VxAtIP*tau_mu*convStabSign_);
			auxJ+=(xiDy*VyAtIP*tau_mu*convStabSign_); 
			if(spaceDim==3)
				auxJ+=(xiDz*VzAtIP*tau_mu*convStabSign_); 
		}
		if(abs(viscStabSign_) > 1e-5){
			auxJ+=(xiDxx*tau_m*kinematicViscosity_*viscStabSign_); 
			auxJ+=(xiDyy*tau_m*kinematicViscosity_*viscStabSign_); 
			if(spaceDim==3)
				auxJ+=(xiDzz*tau_m*kinematicViscosity_*viscStabSign_); 
		}
		auxJ*=(intWeights[actIntPt-1] * jacDet);

		A1.DyadicMult(auxJ,xi);
		// assemble element matrix
//#pragma omp critical
		{
			locElemMat.AddSubMatrix(A1,0,      0);
			locElemMat.AddSubMatrix(A1,nrFncs,nrFncs);
			if(spaceDim==3)
				locElemMat.AddSubMatrix(A1,2*nrFncs,2*nrFncs);
		}
	}
	ResortElementMatrix(elemMat, locElemMat, nrFncs, N);
}

void FluidMechPlaneMassInt_UQ::CalcElementMatrix( Matrix<Double>& elemMat,
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

	bool computeTaus;
	if(stabilParams_[elemNumber-1][3] < -0.001){
		computeTaus=true;
	}
	else{
		computeTaus=false;
	}

	Vector<Double>  Vx, Vy, Vz;  
	Double  VL2, VL2AtIp, VMax;
	Double tau_m, tau_mu, tau_c;
	Double lambda_k, A_elem, h_k;

	Vx.Resize(nrFncs);//Vx.Init(0.0);
	Vy.Resize(nrFncs);//Vy.Init(0.0);
	if(spaceDim==3)
		Vz.Resize(nrFncs);//Vz.Init(0.0);

	for (UInt i=0; i<nrFncs; i++) {
		Vx[i]=elemResult_[0][i];
		Vy[i]=elemResult_[1][i];
		if(spaceDim==3)
			Vz[i]=elemResult_[2][i];
	}

	const Vector<Double> & intWeights = ptelem->GetIntWeights();  

	Double jacDet, multAux;
	UInt N=spaceDim;  // DOFs per Node (Ux,Uy)

	// derivation of shape functions after global coordinates 
	Matrix<Double> locElemMat;

	Matrix<Double> D_ea, D_eb, D_ec;

	Matrix<Double> xiDxDy;

	Vector<Double>  xiDx, xiDy, xiDz;
	Vector<Double>  xi;
	Vector<Double>  auxJ;

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
			Error("Wrong spaceDim!",__FILE__,__LINE__);

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
			Error("Wrong spaceDim!",__FILE__,__LINE__);


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
		//        Warning("one stabilization parameter is smaller than 1e-15",__FILE__,__LINE__);
		//        std::cerr << "tau_c="  << tau_c << "\n" 
		//                  << "tau_m="  << tau_m << "\n" 
		//                  << "tau_mu=" << tau_mu << "\n";
		//      }
		//      if (tau_c>1e3 || tau_m>1e3 || tau_mu>1e3) {
		//        Warning("one stabilization parameter is larger than 1e3",__FILE__,__LINE__);
		//        std::cerr << "tau_c="  << tau_c << "\n" 
		//                  << "tau_m="  << tau_m << "\n" 
		//                  << "tau_mu=" << tau_mu << "\n";
		//      }
	}

	// set matrix to desired size and set all elements to zero
	//elemMat.Resize(nrFncs,nrFncs*N); elemMat.Init(0.0);

	locElemMat.Resize(nrFncs,nrFncs*N); locElemMat.Init(0.0);
	if(abs(presStabSign_)>1e-5){
//#pragma omp parallel for private(auxJ,jacDet,multAux,xiDxDy,xiDx,xiDy,xiDz,xi,D_ea,D_eb,D_ec)
		for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++){
			jacDet = 0;

			ptelem->GetShFncAtIp(xi, actIntPt, ent1.GetElem());
			ptelem->GetGlobDerivShFncAtIp(xiDxDy, actIntPt, ptCoord_, jacDet, ent1.GetElem());

			xiDx.Resize(nrFncs); //xiDx.Init(0.0);
			xiDy.Resize(nrFncs); //xiDy.Init(0.0);
			if(spaceDim==3)
				xiDz.Resize(nrFncs); //xiDz.Init(0.0);

			for (UInt i=0; i<nrFncs; i++) {
				xiDx[i] = xiDxDy[i][0];
				xiDy[i] = xiDxDy[i][1];
				if(spaceDim==3)
					xiDz[i] = xiDxDy[i][2];
			}

			multAux=tau_m * presStabSign_ * intWeights[actIntPt-1] * jacDet;
			auxJ=xiDx*multAux;
			D_ea.DyadicMult(auxJ,xi);

			auxJ=xiDy*multAux;
			D_eb.DyadicMult(auxJ,xi);

			if(spaceDim==3){
				auxJ=xiDz*multAux;
				D_ec.DyadicMult(auxJ,xi);
			}
//#pragma omp critical
			{
				locElemMat.AddSubMatrix(D_ea,0,0);
				locElemMat.AddSubMatrix(D_eb,0,nrFncs);
				if(spaceDim==3)
					locElemMat.AddSubMatrix(D_ec,0,2*nrFncs);
			}
		}
		ColResortElementMatrix(elemMat, locElemMat, nrFncs, nrFncs, 1,N);
	}
	else
		elemMat.Resize(nrFncs,nrFncs*N); elemMat.Init(0.0);

}
//*********************************************************************************************
//********* mixed FEM *************************************************************************
//*********************************************************************************************
void FluidMechPlaneMixedMassInt_UV::CalcElementMatrix( Matrix<Double>& elemMat,
		EntityIterator& ent1,EntityIterator& ent2 ) {

	// Extract pointer to reference element and get coordinates
	ExtractElemInfo( ent1 );

	ptelem->SetAnsatzFct( ansatzFct1_ );
	const UInt nrFncs = ptelem->GetNumFncs( ansatzFct1_ );
	const UInt nrIntPts= ptelem->GetNumIntPoints();

	const UInt spaceDim = ptelem->GetDim();  

	const Vector<Double> & intWeights = ptelem->GetIntWeights();  

	Double jacDet;

	// derivation of shape functions after global coordinates 
	Matrix<Double> locElemMat;

	Matrix<Double> A1;
	Vector<Double> xi;
	Vector<Double> auxJ;

	// set matrix to desired size and set all elements to zero
	//elemMat.Resize(nrFncs*N); elemMat.Init(0.0);
	locElemMat.Resize(nrFncs*spaceDim); locElemMat.Init(0.0);

	for (Integer actIntPt=1; actIntPt <= nrIntPts; actIntPt++)
	{
		jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ptCoord_,ent1.GetElem() );
		ptelem->GetShFncAtIp(xi, actIntPt, ent1.GetElem());

		// Mass Integrator
		auxJ=(xi);
		auxJ*=(intWeights[actIntPt-1] * jacDet);

		A1.DyadicMult(auxJ,xi);
		// assemble element matrix
		{
			locElemMat.AddSubMatrix(A1,0,      0);
			locElemMat.AddSubMatrix(A1,nrFncs,nrFncs);
			if(spaceDim==3)
				locElemMat.AddSubMatrix(A1,2*nrFncs,2*nrFncs);
		}
	}
	ResortElementMatrix(elemMat, locElemMat, nrFncs, spaceDim);
}
} // end namespace CoupledField
