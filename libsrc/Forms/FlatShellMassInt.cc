#include <iostream>
#include <fstream>

#include "FlatShellMassInt.hh"

namespace CoupledField {

  //Constructor for normal material
  FlatShellMassInt::FlatShellMassInt(BaseMaterial * matData)
    : FlatShellInt(matData) 
  {
    
    ENTER_FCN( "FlatShellMassInt::FlatShellMassInt" );
    
    name_ = "FlatShellMassInt";
    matData->GetScalar(density_,DENSITY,REAL);
    baseType_ = MASS;
    
  }

  //Constructor for compostie material 
  FlatShellMassInt::FlatShellMassInt( Composite * composite )
    : FlatShellInt(composite) 
  {
    ENTER_FCN( "FlatShellMassInt::FlatShellMassInt" );

    name_ = "FlatShellMassInt";
    baseType_ = MASS;
  }
  
  //Destructors
  
  FlatShellMassInt::~FlatShellMassInt()
  {
    ENTER_FCN( "FlatShellMassInt::~FlatShellMassInt" );
  }

  void FlatShellMassInt::CalcElementMatrix(Matrix<Double>& elemMat,
					   EntityIterator& ent1, 
					   EntityIterator& ent2 )
  {
    ENTER_FCN( "FlatShellMassInt::CalcElemMatrix" );
  
    // Extract pointer to reference element and get coordinates
    ExtractElemInfo( ent1 );

    const UInt nrIntPts= ptelem->GetNumIntPoints();
    const UInt nrNodes = ptelem->GetNumNodes();
    const UInt nrDofs = getNrDofs();
    const Vector<Double> & intWeights = ptelem->GetIntWeights();  
    Double jacDet=0.0, t=thickness_;
    
    Double thickness_all=0.0, fac_alpha=0.0, fac_beta=0.0, fac_gamma=0.0;
        
    Vector<Double> shapeFncAtIp;
    Matrix<Double> partElemMat;
    Matrix<Double> TransMat;
    Matrix<Double> ShellCoord; 

    // set matrix to desired size and set all elements to zero
    //    partElemMat.Resize(nrNodes);
    //See Finite Element Modeling of smart cantilever plate page 167
    elemMat.Resize(nrNodes * nrDofs);
    elemMat.Init();
    partElemMat.Resize(nrNodes * nrDofs);
    partElemMat.Init();
    
    //Transformation from 3d to 2d
    FlatShellInt::CoordTrans(ptCoord_, TransMat, ShellCoord );
    
    if ( isComposite_ == TRUE ) {
    
      //Initialisation of the composite structure variables
      const UInt nrLayers  = composite_->thickness.GetSize();
        
      z_.Resize(nrLayers + 1);
      z_.Init();
         
      for (UInt k=0; k < nrLayers; k++) {
        thickness_all += composite_->thickness[k];
      }
    
      z_[0]= - (thickness_all/2.0);
        
      for(UInt i=1; i <= nrLayers; i++) {
    	z_[i] = z_[i-1] + composite_->thickness[i-1];
      }
    
      //std::cout << "z coordinates :"<< z_ << std::endl;
    
      //Loop over the layers to calculate the additional factors for the mass matrix            
      for (UInt k=1 ; k <= nrLayers ; k++) {
        
        composite_->materials[k-1]->GetScalar(density_,DENSITY,REAL);

        fac_alpha += density_*(z_[k] - z_[k-1]);
        fac_beta  += density_*((z_[k])*(z_[k]) - (z_[k-1])*(z_[k-1]));
        fac_gamma += density_*(z_[k]*z_[k]*z_[k] - z_[k-1]*z_[k-1]*z_[k-1]);
      }
      //std::cout << "fac_alpha is :"<< fac_alpha << std::endl;
      //std::cout << "fac_beta is :"<< fac_beta << std::endl;
      //std::cout << "fac_gamma is :"<< fac_gamma << std::endl;
        
      for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
           
        jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ShellCoord);
        
        ptelem-> GetShFncAtIp(shapeFncAtIp, actIntPt);
        for (UInt i = 0; i < nrNodes ; i++ ) {
	
	  for (UInt j=0; j < nrNodes ; j++ ) {
                 
	    partElemMat[nrDofs*i][nrDofs*j]         = fac_alpha*shapeFncAtIp[i]*shapeFncAtIp[j];
	    partElemMat[nrDofs*i + 1][nrDofs*j + 1] = fac_alpha*shapeFncAtIp[i]*shapeFncAtIp[j];
	    partElemMat[nrDofs*i + 2][nrDofs*j + 2] = fac_alpha*shapeFncAtIp[i]*shapeFncAtIp[j];
               
	    partElemMat[nrDofs*i][nrDofs*j + 4]     = -0.5*(fac_beta*shapeFncAtIp[i]*shapeFncAtIp[j]);
	    partElemMat[nrDofs*i + 1][nrDofs*j + 3] = -0.5*(fac_beta*shapeFncAtIp[i]*shapeFncAtIp[j]);
               
	    partElemMat[nrDofs*i + 3][nrDofs*j + 1] = -0.5*(fac_beta*shapeFncAtIp[i]*shapeFncAtIp[j]);
	    partElemMat[nrDofs*i + 4][nrDofs*j]     = -0.5*(fac_beta*shapeFncAtIp[i]*shapeFncAtIp[j]);
       
	    partElemMat[nrDofs*i + 3][nrDofs*j + 3] = 1.0/3*(fac_gamma*shapeFncAtIp[i]*shapeFncAtIp[j]);
	    partElemMat[nrDofs*i + 4][nrDofs*j + 4] = 1.0/3*(fac_gamma*shapeFncAtIp[i]*shapeFncAtIp[j]);
                        
	  }
	}
       
	partElemMat *= intWeights[actIntPt-1] * jacDet;
        
	elemMat += partElemMat;
       
      }
      //std::cout << "The element Matrix is :"<< elemMat << std::endl;      
    
    } else {
    
      Double z  =  t/2.0;
      fac_alpha = density_* t ;
      fac_beta  = density_*((z)*(z) - (z)*(z));
      fac_gamma = density_*((z)*(z)*(z) + (z)*(z)*(z));
    
      //std::cout << "fac_alpha is :"<< fac_alpha << std::endl;
      //std::cout << "fac_beta is :"<< fac_beta << std::endl;
      //std::cout << "fac_gamma is :"<< fac_gamma << std::endl;
         
    
      for (UInt actIntPt=1; actIntPt <= nrIntPts; actIntPt++) {
           
	jacDet = ptelem->CalcJacobianDetAtIp(actIntPt, ShellCoord);
        
	ptelem-> GetShFncAtIp(shapeFncAtIp, actIntPt);
      
	for (UInt i = 0; i < nrNodes ; i++ ) {
	
	  for (UInt j=0; j < nrNodes ; j++ ) {
                 
	    partElemMat[nrDofs*i][nrDofs*j]         = fac_alpha*shapeFncAtIp[i]*shapeFncAtIp[j];
	    partElemMat[nrDofs*i + 1][nrDofs*j + 1] = fac_alpha*shapeFncAtIp[i]*shapeFncAtIp[j];
	    partElemMat[nrDofs*i + 2][nrDofs*j + 2] = fac_alpha*shapeFncAtIp[i]*shapeFncAtIp[j];

	    partElemMat[nrDofs*i][nrDofs*j + 4]     = -0.5*(fac_beta*shapeFncAtIp[i]*shapeFncAtIp[j]);
	    partElemMat[nrDofs*i + 1][nrDofs*j + 3] = -0.5*(fac_beta*shapeFncAtIp[i]*shapeFncAtIp[j]);

	    partElemMat[nrDofs*i + 3][nrDofs*j + 1] = -0.5*(fac_beta*shapeFncAtIp[i]*shapeFncAtIp[j]);
	    partElemMat[nrDofs*i + 4][nrDofs*j]     = -0.5*(fac_beta*shapeFncAtIp[i]*shapeFncAtIp[j]);

	    partElemMat[nrDofs*i + 3][nrDofs*j + 3] = 1.0/3*(fac_gamma*shapeFncAtIp[i]*shapeFncAtIp[j]);
	    partElemMat[nrDofs*i + 4][nrDofs*j + 4] = 1.0/3*(fac_gamma*shapeFncAtIp[i]*shapeFncAtIp[j]);
                        
	  }
	}
      
	partElemMat *= intWeights[actIntPt-1] *  jacDet;
        
	elemMat += partElemMat;
      }
      
      //std::cout << "The element Matrix is :"<< elemMat << std::endl;
    }
  
    FlatShellInt::LocaltoGlob(elemMat,TransMat);
  }
  
  
  
} // end namespace CoupledField



