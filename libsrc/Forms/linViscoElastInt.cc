#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "linViscoElastInt.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

namespace CoupledField
{


  LinViscoElastInt::LinViscoElastInt(BaseMaterial* matData, SubTensorType type, 
				     std::string matrixType, Double timeStep) 
    : linElastInt(matData, type)
  {
    ENTER_FCN( "LinViscoElastInt::LinViscoElastInt" );

    name_ = "LinViscoElastInt";

    SetDimensions(type);

    matrixType_ = matrixType;
    timeStep_ = timeStep;

    ptMaterial->GetScalar(dampAlpha_,RAYLEIGH_ALPHA,REAL);
    ptMaterial->GetScalar(dampBeta_,RAYLEIGH_BETA,REAL);

    StdVector<Double> fracDerivList_;
    params->GetList( "fracDeriv", fracDerivList_, "mechanic", "damping" );
    fracDeriv_ = fracDerivList_[0];
    timeStepPowerFracDeriv_ = std::pow(timeStep_,-fracDeriv_);

    StdVector<Double> elastModuleList_;
    params->GetList( "ElastModule", elastModuleList_, "mechanic", "damping" );
    elastModule_ = elastModuleList_[0];
  }
 


  LinViscoElastInt::~LinViscoElastInt()
  {
    ENTER_FCN( "LinViscoElastInt::~LinViscoElastInt" );

  }

 void LinViscoElastInt::calcDMat(Matrix<Double> & dMat)
  {
    ENTER_FCN( "LinViscoElastInt::calcDMat" );

    double E = elastModule_;
    double val = 0.0;
    
    Matrix<Double> aMat;
    Matrix<Double>  cMat;

    // compute matrix inverse(A) - 
    aMat.Resize(getDimD());
    aMat.Init();
    
    // set entries on the diagonal
    val = timeStepPowerFracDeriv_ * dampAlpha_;
    
    // set the entries on the diagonal matrix, to get the inverse, 
    // set the values on the diagonal are 1/val
    for(UInt i=0;i<getDimD();i++)
      {
        // aMat.SetEntry(i,i,1.0/(val));
        aMat.SetEntry(i,i,1.0/(val+1.0)); // The entries are so implemented that the A mtrix is already invers
      }     
  
    // compute material tensor
    calcDMat ( cMat );
	
    // differentiate, whether the modified stiffness martrix 
    // or the damping matrix for the right hand side is needed
   if(matrixType_ == "modifiedStiffness")
     { 
       val = (dampBeta_/E) * timeStepPowerFracDeriv_;
       cMat  = ((cMat *  val) + cMat);

       // compute D-matrix (inverse(A)*C)
       dMat = aMat * cMat;
     }
   else if(matrixType_ == "MatDepRHSMatrix")
     {
       val = (dampBeta_/E);
       cMat  = (cMat *  val);
       dMat = aMat * cMat;
     }      
    else
	{
	  Error("wrong matrixType_ specified", __FILE__, __LINE__);  
	}          
  }


void LinViscoElastInt::GetModifiedMaterialMat(Matrix<Double>& cMat)
   { 
     double val = 0.0;

     std::cerr << "LINVISCOELASTINT::GETMODIFIEDMATERIALMAT" << std::endl;

     StdVector<std::string> keyVec;
     std::string typeDampBeta;
     keyVec = "pdeList", "mechanic","region","damping","typeDampBeta";    
     params->Get(keyVec,typeDampBeta);    

     // compute material tensor
     calcDMat ( cMat );
     
     if(typeDampBeta == "CMat")
      {
	val = (dampBeta_/elastModule_)*timeStepPowerFracDeriv_;
	cMat  = (cMat *  val) + cMat;	   
      }
     else
      {
	val = dampBeta_ * timeStepPowerFracDeriv_;
	Matrix<Double> betaMat;
	betaMat.Resize(getDimD());
        betaMat.Init();
	//val = val/100;
	for(UInt i=0; i< getDimD();i++)
	  {
	    betaMat[i][i]  = val;
	  }	   
	cMat = cMat  + betaMat;
      }
   } 



  void LinViscoElastInt::CalcAInvMat(Matrix<Double>& aInvMat)
  {

    std::cerr << "LINVISCOELASTINT::calcAInvMat" << std::endl;
    
    aInvMat.Resize(getDimD(),true);
    aInvMat.Init();
    Double val;
    val = timeStepPowerFracDeriv_ * dampAlpha_;
	
     // set the emtries on the diagonal matrix, to get the inverse, 
     // the value on the diagonal are 1/val

     for(UInt i=0;i<getDimD();i++)
       {
	 aInvMat.SetEntry(i,i,1/(val+1));
	 //	 aInvMat.SetEntry(i,i,1/(val));
       }
  } 


  void LinViscoElastInt::SetDimensions(SubTensorType type) {

    ENTER_FCN( "linElastInt::SetDimensions" );

    if ( type == FULL ) {
      dimD_   = 6;
      nrDofs_ = 3;
    }
    else if ( type == PLANE_STRAIN ) {
      dimD_   = 3;
      nrDofs_ = 2;
    }
    else if ( type == AXI ) {
      dimD_   = 4;
      nrDofs_ = 2;
      isaxi_  = true;
    }
    else {
      Error("wrongh Type(axi,planeStrein,3d) specified", __FILE__, __LINE__);  
    }
  }

} // end namespace CoupledField


