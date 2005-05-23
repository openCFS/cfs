#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <stdio.h>

//#include <limits.h>

#include "MaterialData.hh"
#include "LoadMaterialData.hh"
#include "WriteInfo.hh"
#include "General/environment.hh"

namespace CoupledField
{

} // end namespace CoupledField




  /*
    void LoadMaterialData :: EulerAnglesRotate(MaterialData * material, const Vector<Double> & eulerAngles)
    {
    #ifdef TRACE
    if (trace) (*trace) << "Entering  LoadMaterialData::EulerAnglesRotate" << std::endl;
    #endif    const Double epsIndex = 6;
    const Double pi = acos(-1);
    const Double gradToRad = pi / 180;
    const Double epsilon = 1e-10;

    // indizes in for loops
    Integer i, j, k, l, m, n, o;

    Double trans[3][3];        // transformation matrix, 
    //see "Huette: Die Grundlagen der Ingeneurwissenschaften, E3"

    Double co[3][3][3][3];     // stiffness tensor
    Double cr[3][3][3][3];     // rotated stiffness tensor

    Double eo[3][3][3];        // tensor of piezoelectric moduli
    Double er[3][3][3];        // rotated tensor of piezoelectric moduli

    Double epso[3][3];         // permittivity tensor
    Double epsr[3][3];         // rotated permittivity tensor


    Double help2[3][3];        // 2D help matrix
    Double help3[3][3][3];     // 3D help matrix
    Double help4[3][3][3][3];  // 4D help matrix


    Integer abbrevSubscript[3][3] = 
    {{1 ,6 ,5},
    {6, 2 ,4},
    {5, 4, 3}};


    Matrix<Double> * fullOrigMat = material->GetMatrix();


    for ( i=0; i<3; i++)
    for ( j=0; j<3; j++)
    {
    epso[i][j] = (*fullOrigMat)((Integer)(i + 1 + epsIndex),
    (Integer)( j + 1+ epsIndex)); 

    for ( k=0; k<3; k++)
    {
    eo[i][j][k] = (*fullOrigMat)((Integer)(i + 1 + epsIndex), 
    (Integer)(abbrevSubscript[j][k]));

    for ( l=0; l<3; l++)
    co[i][j][k][l] = (*fullOrigMat)((Integer)(abbrevSubscript[i][j]), 
    (Integer)(abbrevSubscript[k][l]));
    }
    }

    //     save euler angles.
  
    Double  lambda = eulerAngles[1];
    Double  mu     = eulerAngles[2];
    Double  theta  = eulerAngles[3];

    Double  sinl = sin(lambda * gradToRad);
    Double  cosl = cos(lambda * gradToRad);
    Double  sinm = sin(mu * gradToRad);
    Double  cosm = cos(mu * gradToRad);
    Double  sint = sin(theta * gradToRad);
    Double  cost = cos(theta * gradToRad);

    // set trigonometric values to zero if they are smaller than eps

    sinl = fabs(sinl) < epsilon ? 0 : sinl ;
    cosl = fabs(cosl) < epsilon ? 0 : cosl ;
    sinm = fabs(sinm) < epsilon ? 0 : sinm ;
    cosm = fabs(cosm) < epsilon ? 0 : cosm ;
    sint = fabs(sint) < epsilon ? 0 : sint ;
    cost = fabs(cost) < epsilon ? 0 : cost ;

    //     calculate transformation matrix trans

    trans[0][0] = cosl * cost - sinl * cosm * sint;
    trans[0][1] = sinl * cost + cosl * cosm * sint;
    trans[0][2] = sinm * sint;
    trans[1][0] = -cosl * sint - sinl * cosm * cost;
    trans[1][1] = -sinl* sint + cosl * cosm * cost;
    trans[1][2] = sinm * cost;
    trans[2][0] = sinl * sinm;
    trans[2][1] = -cosl * sinm;
    trans[2][2] = cosm;


    std::cout << "================== Euler Angles =======================  " << std::endl
    << std::endl;
    << "lambda = " << lambda << std::endl;
    << "mu     = " << mu << std::endl;
    << "theta  = " << theta << std::endl << std::endl;
    

    std::cout << std::endl << "Rotation matrix = " << std::endl;
    for(l=0; l<3; l++)
    {
    for(o=0; o<3; o++)
    std::cout << std::setw(12) << trans[l][o] << "\t";
    std::cout<< std::endl;
    }
    std::cout << std::endl << std::endl;
  
  

    //     ========== calculate tensors for rotated system ==================

    //     first coordinate

    for(l=0; l<3; l++)
    for(o=0; o<3; o++)
    {
    help2[o][l] = trans[l][0] * epso[o][0]+
    trans[l][1] * epso[o][1]+
    trans[l][2] * epso[o][2];
    for(n=0; n<3; n++)
    {
    er[n][o][l] = trans[l][0] * eo[n][o][0]+
    trans[l][1] * eo[n][o][1]+
    trans[l][2] * eo[n][o][2];

    for(m=0; m<3; m++)
    help4[m][n][o][l] = trans[l][0] * co[m][n][o][0]+
    trans[l][1] * co[m][n][o][1]+
    trans[l][2] * co[m][n][o][2];
    }
    }
    
  
    //     second coordinate
  
    for(l=0; l<3; l++)
    for(k=0; k<3; k++)
    {
    epsr[k][l] = trans[k][0] * help2[0][l]+
    trans[k][1] * help2[1][l]+
    trans[k][2] * help2[2][l];

    for(n=0; n<3; n++)
    {
    help3[n][k][l] = trans[k][0] * er[n][0][l]+
    trans[k][1] * er[n][1][l]+
    trans[k][2] * er[n][2][l];

    for(m=0; m<3; m++)
    cr[m][n][k][l] = trans[k][0] * help4[m][n][0][l]+
    trans[k][1] * help4[m][n][1][l]+
    trans[k][2] * help4[m][n][2][l];
    }
    }
    

    //     third coordinate

    for(l=0; l<3; l++)
    for(k=0; k<3; k++)
    for(j=0; j<3; j++)
    {
    er[j][k][l] = trans[j][0] * help3[0][k][l] +
    trans[j][1] * help3[1][k][l]+
    trans[j][2] * help3[2][k][l];

    for(m=0; m<3; m++)
    {
    help4[m][j][k][l] = trans[j][0] * cr[m][0][k][l]+
    trans[j][1] * cr[m][1][k][l]+
    trans[j][2] * cr[m][2][k][l];
    }
    }

    //     fourth coordinate

    for(l=0; l<3; l++)
    for(k=0; k<3; k++)
    for(j=0; j<3; j++)
    for(i=0; i<3; i++)
    cr[i][j][k][l] = trans[i][0] * help4[0][j][k][l]
    + trans[i][1] * help4[1][j][k][l]
    + trans[i][2] * help4[2][j][k][l];

    // translate back into abbreviated subscripts
    for ( i=0; i<3; i++)
    for ( j=0; j<3; j++)
    {
    (*fullOrigMat)( (Integer)( i + 1 + epsIndex), 
    (Integer)(j + 1 + epsIndex)) =  epsr[i][j];

    for ( k=0; k<3; k++)
    {
    (*fullOrigMat)((Integer)(i + 1 + epsIndex), 
    abbrevSubscript[j][k]) =  er[i][j][k];
    (*fullOrigMat)(abbrevSubscript[j][k], 
    (Integer)(i + 1 + epsIndex)) = er[i][j][k];

    for ( l=0; l<3; l++)
    (*fullOrigMat)(abbrevSubscript[i][j], abbrevSubscript[k][l]) =  cr[i][j][k][l];
    }
    }
    }
  */



// ============================= nonlin: not yet used =========================================
/*
  void LoadMaterialData :: ReadMagNonLin(std::ifstream & fin, MaterialData * material)
  {
  Integer i, k;
  Integer materialNumber, nonlin;
  Integer nrItemsRead;
  Double eModule,nue;
  Double conductivity;
  Integer intervalls;
  char buffer[bufLength];  
  char materialType[bufLength];
  char materialName[bufLength];
  std::istringstream * strPtr;
  Vector <Double> * coeffPtr;
  Vector <NonlinSpline*> * splinePtr;


  ReadLine(fin,buffer);
  nrItemsRead = sscanf(buffer,"%d%s%s%d", &materialNumber, materialType, materialName, &nonlin);
  if (nrItemsRead != 4)
  std::cout << "   *** The materialfile is corrupt! *** Material: " << materialName << std::endl;
  material -> SetName(materialName);

  ReadLine(fin,buffer);
  nrItemsRead = sscanf(buffer,"%lf%lf", &eModule, &nue);
  if (nrItemsRead != 2)
  std::cout << "   *** The materialfile is corrupt! ***  Material: " << materialName  << std::endl;
  material -> SetEModule(eModule);
  material -> SetNu(nue);

  ReadLine(fin,buffer); 
  nrItemsRead = sscanf(buffer,"%lf", &conductivity);
  if (nrItemsRead != 1)
  std::cout << "   *** The materialfile is corrupt! *** Material: " << materialName << std::endl;
  material -> SetConductivity(conductivity);

  ReadLine(fin,buffer); 
  nrItemsRead = sscanf(buffer,"%i", &intervalls);
  if (nrItemsRead != 1)
  std::cout << "   *** The materialfile is corrupt! ***  Material: " << materialName << std::endl;

  splinePtr = new Vector <NonlinSpline*>;       
   
  for (i=1;i<=intervalls;i++)
  {
  Double lowerLimit,upperLimit;
  coeffPtr = new Vector <Double>;
  Integer degree;
  ReadLine(fin,buffer);
  strPtr = new std::istringstream(buffer);
  *strPtr >> lowerLimit >> upperLimit >> degree;
  if (strPtr->fail())
  std::cout << "*** The materialfile is corrupt! ***  Material: " << materialName << std::endl;
  strPtr = new std::istringstream(buffer);
  *strPtr >> lowerLimit >> upperLimit >> degree;
  if (strPtr->fail())
  std::cout << "*** The materialfile is corrupt! ***  Material: " << materialName << std::endl; 
  coeffPtr->SetSize(degree+1);
  for (k=1;k<=degree+1;k++)
  *strPtr >> coeffPtr -> Elem(k);
  if (strPtr->fail())
  std::cout << "*** The materialfile is corrupt! ***  Material: " << materialName << std::endl;
  splinePtr -> Append(new NonlinSpline(lowerLimit,upperLimit,*coeffPtr));
  delete strPtr;          
  }     
  material->SetNonLinSpline(splinePtr);
  //  fin.getline(buffer,bufLength,'\n');

  // setting the nonlin variable must be done after the spline is read, because in the linear
  // case the first parameter of the spline is used!
  material->DefLin(nonlin);
  }





  // const Vector<Double> & eulerAnglesTmp = flags.GetNumListFlag ("eulerAngles");    
  // scaleMatDat = flags.GetDefineFlag ("scaleMatDat");

    
  if (eulerAnglesTmp.size() % 3)
  {
  std::cout << std::endl << "LoadMaterialData: There must be three Euler Angles per material! " << std::endl;
  exit(EXIT_FAILURE);
  }
      
  eulerAngles.SetSize(eulerAnglesTmp.size() / 3);

  for (Integer matNr = 1; matNr<= eulerAnglesTmp.size() / 3; matNr++)
  {
  eulerAngles.Elem(matNr) = new Vector<Double>;
  eulerAngles.Elem(matNr) -> SetSize(3);
  for (Integer i = 1; i<=3; i++)
  eulerAngles.Elem(matNr) -> Elem(i) = eulerAnglesTmp[i + (matNr-1)*3];
  }


*/





