#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <limits.h>
#include <string>

#include "MHMaterialData.hh"
#include "baseMaterial.hh"
#include "../DataInOut/MHMaterialDataFile.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

namespace CoupledField
{

  MHMaterialData::MHMaterialData(StdVector<BaseMaterial*>& ptMaterial)
  {
    ENTER_FCN("MHMaterialData::MHMaterialData");
  
    ptMaterial_.Resize(ptMaterial.GetSize());

    UInt maxP=3;
    

    parameterCoeff_.Resize(10,maxP);
    parameter_.Resize(10);

    parameterCoeff_[0][0]=0.5; // c_11, a_0=1, andere a_i=0, somit konstant
    parameterCoeff_[1][0]=1.0; // c_33, a_0=1
    parameterCoeff_[1][1]=0.5; // c_33, a_1=0.5
    // parameterCoeff_[1][2]=0.25; // c_33, a_2=0.25 Polynom vom Grad 3
    parameterCoeff_[2][0]=1.0; // c_12, a_0=1, andere a_i=0, somit konstant
    parameterCoeff_[3][0]=1.0; // c_13, a_0=1, andere a_i=0, somit konstant
    parameterCoeff_[4][0]=0.5; // c_44, a_0=1, andere a_i=0, somit konstant
    parameterCoeff_[5][0]=1.0; // e_15, a_0=1, andere a_i=0, somit konstant
    parameterCoeff_[6][0]=1.0; // e_13, a_0=1, andere a_i=0, somit konstant
    parameterCoeff_[7][0]=1.0; // e_33
    parameterCoeff_[7][1]=-0.5; // e_33 Polynom vom Grad 2
    parameterCoeff_[8][0]=1.0; // eps_11, a_0=1, andere a_i=0, somit konstant
    parameterCoeff_[9][0]=3.0; // eps_33, 
    parameterCoeff_[9][1]=1.0; // eps_33, Polynom vom Grad 
    // parameterCoeff_[9][2]=-0.5; // eps_33, Polynom vom Grad 
    //    parameterCoeff_[9][3]=0.1; // eps_33, Polynom vom Grad 4

    ptMHMaterialDataFile = new MHMaterialDataFile();
    ptMHMaterialDataFile->readData();

    parameterR_ = ptMHMaterialDataFile->parameter;
    parameterC_ = ptMHMaterialDataFile->parameterC;

    for (UInt i=0;i<parameterR_.GetSize();i++)
      parameter_[i] = Complex(parameterR_[i],parameterC_[i]);

  }

 
  MHMaterialData::~MHMaterialData()
  {
    
    ENTER_FCN("MHMaterialData::~MHMaterialData");
    //    delete ptMHMaterialDataFile;
    delete ptMHMaterialDataFile;
  }


    void MHMaterialData::updateMaterialData(Vector<Complex> & parameter, 
					    StdVector<BaseMaterial*> ptMaterial){
    ENTER_FCN("MHMaterialData::updateMaterialData");    

    Error("Do not work with new material class!",__FILE__,__LINE__);

//     for(UInt i=0;i<9;i++)
//       for(UInt j=0;j<9;j++)
//         ptMaterial->SetPiezoMatrixData(i,j,0.0);

//     ptMaterial[0].SetPiezoMatrixData(0,0, parameter[0].real());
//     ptMaterial[0].SetPiezoMatrixData(1,1, parameter[0].real());
//     ptMaterial[0].SetPiezoMatrixData(2,2, parameter[1].real());
//     ptMaterial[0].SetPiezoMatrixData(0,1, parameter[2].real());
//     ptMaterial[0].SetPiezoMatrixData(1,0, parameter[2].real());
//     ptMaterial[0].SetPiezoMatrixData(0,2, parameter[3].real());
//     ptMaterial[0].SetPiezoMatrixData(2,0, parameter[3].real());
//     ptMaterial[0].SetPiezoMatrixData(1,2, parameter[3].real());
//     ptMaterial[0].SetPiezoMatrixData(2,1, parameter[3].real());
//     ptMaterial[0].SetPiezoMatrixData(3,3, parameter[4].real());
//     ptMaterial[0].SetPiezoMatrixData(4,4, parameter[4].real());
//     // std::cout<<"updateMaterialData Set Data 44"<<std::endl;
//     ptMaterial[0].SetPiezoMatrixData(5,5, 0.5*(parameter[0].real()-parameter[2].real()));
//     ptMaterial[0].SetPiezoMatrixData(6,4, parameter[5].real());
//     ptMaterial[0].SetPiezoMatrixData(7,3, parameter[5].real());
//     ptMaterial[0].SetPiezoMatrixData(4,6, parameter[5].real());
//     ptMaterial[0].SetPiezoMatrixData(3,7, parameter[5].real());
//     ptMaterial[0].SetPiezoMatrixData(8,0, parameter[6].real());
//     ptMaterial[0].SetPiezoMatrixData(8,1, parameter[6].real());
//     ptMaterial[0].SetPiezoMatrixData(0,8, parameter[6].real());
//     ptMaterial[0].SetPiezoMatrixData(1,8, parameter[6].real());
//     ptMaterial[0].SetPiezoMatrixData(8,2, parameter[7].real());
//     ptMaterial[0].SetPiezoMatrixData(2,8, parameter[7].real());
//     ptMaterial[0].SetPiezoMatrixData(6,6, parameter[8].real());
//     ptMaterial[0].SetPiezoMatrixData(7,7, parameter[8].real());
//     ptMaterial[0].SetPiezoMatrixData(8,8, parameter[9].real());

   
    // Consider poling of piezoelectric body
    Double a1, a2, a3;
    a1=a2=a3=0;

    if( params->HasValue( "x", "1", "piezo", "polingDirection" ) )
      a1=1;
    
    if( params->HasValue( "y", "1", "piezo", "polingDirection" ) )
      a2=1;
    
    if( params->HasValue( "z", "1", "piezo", "polingDirection" ) )
      a3=1;
    
    if (a1==0&&a2==0&&a3==0)
      a3=1.0;    // if no poling direction is specified, the z-direction is chosen by default 


    Error("Do not work with new material class!",__FILE__,__LINE__);
   
    //    ptMaterial->RotateMaterialMatrix(a1,a2,a3);
    //    std::cout<<" updated Material successfully " <<std::endl;        
      
   
  } 
//   } // end updateMaterialData

  void MHMaterialData::calcParameterCurveAtElement(Vector<Complex> & parameter, 
						   Matrix<Double> & parameterCoeff_, 
						   UInt element,UInt N,
						   Integer delta, UInt PP){
  std::cout<<"Enter calcParameterCurveAtElement " <<std::endl;
  UInt j;
  Complex prod;
  prod=Complex(1.0,0.0);
  UInt p;
  Double pFac;
  Matrix<Integer> exponent;
  Matrix<Integer> count;
  Double binomCoeffNenner;
  Complex innerSum;
  innerSum=Complex(0,0);
  MHMaterialDataFile *  ptMHFiles_;
  ptMHFiles_ = new MHMaterialDataFile;
//   ptMHFiles_->getExponentArray(exponent,count, N, PP, delta);
//   std::cout<<"exponent:"<<std::endl;
//   std::cout<<exponent<<std::endl;

  Complex elementSolution;
  elementSolution=Complex(2.0,0.5);

  for (UInt i=0;i<parameter.GetSize();i++){
    //    std::cout<<" Parameter Nr = " << i <<std::endl;
    parameter[i]=parameter[i]*parameterCoeff_[i][0];
    j=0;
    p=1;
    for (UInt jj=1; jj<=PP;jj++)    //Summe über Parameter
      if(parameterCoeff_[i][jj]!=0){    //Summe über Parameter
        p++;
        ptMHFiles_->getExponentArray(exponent,count, PP, p , delta);
        binomCoeffNenner=1;
          
        //        std::cout<<"exponent"<<std::endl;
        //std::cout<<exponent<<std::endl;

        for (UInt nrP=0;nrP<exponent.GetSizeRow();nrP++){ // loop over all entries in I(p,delta)
          for (UInt e=0;e<exponent.GetSizeCol();e++){   // loop over all p_n
            Integer pFacNenner=1;
             
             for(Integer ii=1;ii<=exponent[nrP][e];ii++)
               pFacNenner*=ii;    // p_n!
             
             binomCoeffNenner*=pFacNenner;  // p1!*p2!* ... *p_n!
             
             pFac=1;
             for (UInt pInd=1; pInd<=p; pInd++)
               pFac*=pInd; // p!
             
             prod=Complex(1.0,0);
             for (UInt n=0;n<2*N;n++){
               prod=prod*std::pow(elementSolution,exponent[nrP][n]);
             }
             
             prod*=pFac/binomCoeffNenner;
             //             std::cout<<"Product (p/p) * P_-2N,...2N, P= "<<prod<<std::endl;
             binomCoeffNenner=1;
             pFac=0;
           }
           innerSum+=prod;
           //           std::cout<<"inner_sum =  "<<innerSum<<", parameterCoeff["<<i<<"]["<<jj<<"] = " << parameterCoeff_[i][jj]<<std::endl;
           prod=Complex(1,0);
         }
         parameter[i]=parameter[i]+parameterCoeff_[i][jj]*innerSum;                 
         innerSum=Complex(0,0);
      }
    }
//     std::cout<<"parameter:"<<std::endl;
//     std::cout<<parameter<<std::endl;


  }

} // end of class
