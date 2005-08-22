#include <iostream>
#include <fstream>
#include <string>
#include <cmath>

#include "linViscoElastInt.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

namespace CoupledField
{

  LinViscoElastInt::LinViscoElastInt(BaseFE * aptelem, MaterialData & matData, std::string geomType, std::string matrixType, Double timeStep) 
    : linElastInt(aptelem, matData)
  {
    ENTER_FCN( "LinViscoElastInt::LinViscoElastInt" );

    aptelem_ = aptelem;
    geomType_ = geomType;
    matrixType_ = matrixType;
    timeStep_ = timeStep;

    dampAlpha_ = ptMaterial->GetDampingAlfa();
    dampBeta_  = ptMaterial->GetDampingBeta();

    StdVector<Double> fracDerivList_;
    params->GetList( "fracDeriv", fracDerivList_, "mechanic", "damping" );
    fracDeriv_ = fracDerivList_[0];
    //std::cerr <<"fracDeriv_=" << fracDeriv_ <<std::endl;
    timeStepPowerFracDeriv_ = std::pow(timeStep_,-fracDeriv_);

    StdVector<Double> elastModuleList_;
    params->GetList( "ElastModule", elastModuleList_, "mechanic", "damping" );
    elastModule_ = elastModuleList_[0];

  }


  LinViscoElastInt::LinViscoElastInt(MaterialData & matData, std::string geomType, std::string matrixType, Double timeStep) 
    : linElastInt(matData)
  {
    ENTER_FCN( "LinViscoElastInt::LinViscoElastInt" );

    geomType_ = geomType;
    matrixType_ = matrixType;
    timeStep_ = timeStep;

    dampAlpha_  = ptMaterial->GetDampingAlfa();
    dampBeta_  = ptMaterial->GetDampingBeta();
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
    
    // set entries on the diagonal
    val = timeStepPowerFracDeriv_ * dampAlpha_;
    
    // set the entries on the diagonal matrix, to get the inverse, 
    // set the values on the diagonal are 1/val
    for(UInt i=0;i<getDimD();i++)
      {
        // aMat.SetEntry(i,i,1.0/(val));
        aMat.SetEntry(i,i,1.0/(val+1.0)); // The entries are so implemented that the A mtrix is already invers
      }     
  
    // get the material matrix in the right dimension, depending on the geometric type
    if(geomType_ == "axi")
      {
        CalcAxiMaterialMat( cMat,  actOrientation);	
      }
    else if(geomType_ == "planeStrain")
      {
        //   std::cerr << "TEST-FRACTIONAL 1 - orientation " << actOrientation <<  std::endl;
        CalcPlaneStrainMaterialMat( cMat,  actOrientation);	 
      }
    else if(geomType_ == "3d")
      {
        Calc3DMaterialMat(cMat);	
      }      
	
    // differentiate, whether the modified stiffness martrix 
    // or the damping matrix for the right hand side is needed
   if(matrixType_ == "modifiedStiffness")
     { 
       val = (dampBeta_/E) * timeStepPowerFracDeriv_;
       cMat  = ((cMat *  val) + cMat);

       // compute D-matrix (inverse(A)*C)
       dMat.Resize(getDimD());
       dMat = aMat * cMat;
     }
   else if(matrixType_ == "MatDepRHSMatrix")
     {
       val = (dampBeta_/E);
       cMat  = (cMat *  val);
       dMat.Resize(getDimD());
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


     
     if(geomType_ == "axi")
       {
	 CalcAxiMaterialMat( cMat,  actOrientation);	
       }
     else if(geomType_ == "planeStrain")
	  {
	    //   std::cerr << "TEST-FRACTIONAL 1 - orientation " << actOrientation <<  std::endl;
	    CalcPlaneStrainMaterialMat( cMat,  actOrientation);	 	    
	  }
     else if(geomType_ == "3d")
       {
	 Calc3DMaterialMat( cMat);	
       }      
     
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
    
    aInvMat.Resize(getDimD());
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



  UInt LinViscoElastInt::getDimD()
  {
    ENTER_FCN( "LinViscoElastInt::getDimD" );
    
    if(geomType_ == "axi")
      {
    	return 4;
      }
    else if(geomType_ == "planeStrain")
      {
    	return 3;   	
      }
    else if(geomType_ == "3d")
      {
    	return 6;   	
      }
    else
      {
	Error("wrongh geomType(axi,planeStrein,3d) specified", __FILE__, __LINE__);  
      }
    return 0;
  }

  /// returns nr. of degrees of freedom
  UInt LinViscoElastInt::getNrDofs()
  {
    ENTER_FCN( "LinViscoElastInt::getNrDofs" );
   
    if(geomType_ == "axi" || geomType_ == "planeStrain")
      {
    	return 2;
      }
    else if(geomType_ == "3d")
      {
    	return 3;   	
      }
    else 
      {
	Error("wrongh geomType(axi,planeStrein,3d) specified", __FILE__, __LINE__);  
      }    
    return 0;
  }
 

} // end namespace CoupledField


