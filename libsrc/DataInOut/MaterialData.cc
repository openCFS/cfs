#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <limits.h>
#include <string>

#include "MaterialData.hh"


namespace CoupledField
{

MaterialData::MaterialData():scaledMatDat(0)
{
  ENTER_FCN("MaterialData::MaterialData");
  const int stringLength = 100;
  name = new char[stringLength];  
  piezoMatrix  = NULL;
  permeaMatrix = new Matrix<Double>(3,3);
  conducMatrix = new Matrix<Double>(3,3);
}


MaterialData::MaterialData(const MaterialData & mat)
{
  ENTER_FCN("MaterialData::MaterialData(MaterialData)")
  density         = mat.density;
  compressibility = mat.compressibility ;
  damp_alfa       = mat.damp_alfa;
  damp_beta       = mat.damp_beta;
  BoverA          = mat.BoverA;
  eModule         = mat.eModule;
  nu           = mat.nu;
//  permeability = mat.permeability;
//  conductivity = mat.conductivity;
  permMx       = mat.permMx;
  permMy       = mat.permMy;
  permMz       = mat.permMz;
  scaledMatDat = mat.scaledMatDat;
  matNr        = mat.matNr;
  nonlin       = mat.nonlin;

  const int stringLength = 100;
  name = new char[stringLength];  
  SetName(mat.name);

  piezoMatrix  = new Matrix<Double>( *mat.piezoMatrix);
  piezoMatrixC = new Matrix<Double> (*mat.piezoMatrixC);
  permeaMatrix = new Matrix<Double>( *mat.permeaMatrix);
  conducMatrix = new Matrix<Double>( *mat.conducMatrix);
}

void MaterialData::SetConductivity(const double& Conductivity) 
{
  ENTER_FCN( "MaterialData::SetConductivity" );
  SetConductivity(0,0,Conductivity);
  SetConductivity(1,1,Conductivity);
  SetConductivity(2,2,Conductivity);
//  conductivity = Conductivity;
}


void MaterialData::SetConductivity(const Integer& i, const Integer& j, const Double& value)
{
  ENTER_FCN("SetConductivity(int,int,double)");
  (*conducMatrix)(i,j) = value;
}

void MaterialData::GetConductivity(const Integer& i, const Integer& j, Double &value)
{
  ENTER_FCN("GetConductivity(int,int,double)");
  value = (*conducMatrix)(i,j);
}

MaterialData::~MaterialData()
{
  if (piezoMatrix)
   delete piezoMatrix;
  if (permeaMatrix)
   delete permeaMatrix;
  if (conducMatrix)
   delete conducMatrix;
  delete[] name;
}

void MaterialData::SetName(const char* Name)
{
  ENTER_FCN("MaterialData::SetName");
  strcpy(name,Name);
}

char * MaterialData::GetMaterialName()
{
  ENTER_FCN("MaterialData::GetMaterialName");
  return name;
}

void MaterialData::SetPermeability(const Integer& i, const Integer& j, const Double& value)
{
  ENTER_IFCN("MaterialData::SetPermeability(int,int,double");
    (*permeaMatrix)(i,j) = value;
}


void MaterialData::GetPermeability(const Integer& i, const Integer& j, Double &value)
{
  ENTER_IFCN("GetPermeability(int,int,double)");
  value = (*permeaMatrix)(i,j);
}

void MaterialData::RotateMaterialMatrix(const Double& a1, const Double& a2, const Double& a3){
  ENTER_FCN("MaterialData::RotateMaterialMatrix");
  Matrix<Double> R;
  R.Resize(3,3);   // Rotation Matrix
  Matrix<Double> Q;
  Q.Resize(6,6);  // Composed Rotation Matrix
  Matrix<Double> c;
  c.Resize(6,6);
  Matrix<Double> e; 
  e.Resize(6,3);
  Matrix<Double> eps;
  eps.Resize(3,3);
  Matrix<Double> QT;
  QT.Resize(6,6);
  Matrix<Double> RT;
  RT.Resize(3,3);

  // X-dIRECTION
  if(a1==0)
    R[0][1]=R[1][2]=R[2][0]=0.0;
  else if(a1==1)
    R[0][1]=R[1][2]=R[2][0]=1.0;
  else
    R[0][1]=R[1][2]=R[2][0]=std::cos(a1);

  //Y-DIRECTION  
  if(a2==0)
    R[0][2]=R[1][0]=R[2][1]=0.0;  
  else if(a2==1)
   R[0][2]=R[1][0]=R[2][1]=1.0;
  else
    R[0][2]=R[1][0]=R[2][1]=std::cos(a2);

  // Z-DIRECTION
  if(a3==0)
      for (Integer i=0;i<3;i++)
	R[i][i]=0.0;
  else if(a3==1)
      for (Integer i=0;i<3;i++)
	R[i][i]=1.0;
  else
      for (Integer i=0;i<3;i++)
	R[i][i]=std::cos(a3);

  for (Integer i=0;i<3;i++)
    for (Integer j=0;j<3;j++){
      Q[i][j]=R[i][j];
      Q[i+3][j+3]=R[i][j];
      Q[i][j+3]=0;
      Q[i+3][j]=0;
    }
  // Rotate real matrix ...
  for (Integer i=0;i<3;i++)
    for (Integer j=0;j<3;j++){
      c[i][j]=(*piezoMatrix)[i][j];
      c[i+3][j+3]=(*piezoMatrix)[i+3][j+3];
      e[i][j]=(*piezoMatrix)[i][6+j];
      e[i+3][j]=(*piezoMatrix)[i+3][6+j];
      eps[i][j]=(*piezoMatrix)[i+6][j+6];
    }
  Q.Transpose(QT);
  R.Transpose(RT);

  c=Q*c*QT;
  e=Q*e*RT;
  eps=R*eps*RT;

  for (Integer i=0;i<3;i++)
    for (Integer j=0;j<3;j++){
      (*piezoMatrix)[i][j]=c[i][j];
      (*piezoMatrix)[i+3][j+3]=c[i+3][j+3];
      (*piezoMatrix)[i][6+j]=e[i][j];
      (*piezoMatrix)[i+3][6+j]=e[i+3][j];
      (*piezoMatrix)[i+6][j]=e[j][i];
      (*piezoMatrix)[i+6][3+j]=e[j+3][i];
      (*piezoMatrix)[i+6][j+6]=eps[i][j];
    }

  // rotate complex Matrix

  for (Integer i=0;i<3;i++)
    for (Integer j=0;j<3;j++){
      c[i][j]=(*piezoMatrixC)[i][j];
      c[i+3][j+3]=(*piezoMatrixC)[i+3][j+3];
      e[i][j]=(*piezoMatrixC)[i][6+j];
      e[i+3][j]=(*piezoMatrixC)[i+3][6+j];
      eps[i][j]=(*piezoMatrixC)[i+6][j+6];
    }
  Q.Transpose(QT);
  R.Transpose(RT);

  c=Q*c*QT;
  e=Q*e*RT;
  eps=R*eps*RT;

  for (Integer i=0;i<3;i++)
    for (Integer j=0;j<3;j++){
      (*piezoMatrixC)[i][j]=c[i][j];
      (*piezoMatrixC)[i+3][j+3]=c[i+3][j+3];
      (*piezoMatrixC)[i][6+j]=e[i][j];
      (*piezoMatrixC)[i+3][6+j]=e[i+3][j];
      (*piezoMatrixC)[i+6][j]=e[j][i];
      (*piezoMatrixC)[i+6][3+j]=e[j+3][i];
      (*piezoMatrixC)[i+6][j+6]=eps[i][j];
    }

} // end RotateMaterialMatrix

/*
 
void MaterialData::DefLin(const int& notLin)
{

  if (notLin)
    nonlin = 1;
  else
    if (magneticSpline->Size()) // exists magneticSpline already ? 
      constPerm = magneticSpline->Get(1) -> ConstPerm(); 
    else
      cout << "You have to define the magneticSpline before the nonlin variable (in the nonlin case)!" << endl;
}


double MaterialData::GetPermiability(const double& MagB) const
{
  int i;

  i = 1;

    while (magneticSpline->Get(i)->UpperLimit() < MagB)
    {
    i++;
    if (i>magneticSpline->Size())
    {
    cout << endl << " Induction " << MagB << " too high " << endl;
    exit(EXIT_FAILURE);
    }
    }
    
    if (i == 1)
    return magneticSpline->Get(i)->ReducedValue(MagB); //at the first spline no offset is tolerated!
    else
    return magneticSpline->Get(i)->Value(MagB);
}


double MaterialData::GetErrorFunctional(const double& MagB) const
{
  int i;
  double errorfunctional=0;

  i = 1;
  while (magneticSpline->Get(i)->UpperLimit() < MagB)
    {
      if (i==magneticSpline->Size())
	{
	  cout << endl << " Induction too high ";
	  exit(EXIT_FAILURE);
	}
      errorfunctional += magneticSpline->Get(i) -> WholeIntegral();
      i++;
    }
  return errorfunctional += magneticSpline->Get(i)->IntValue(MagB);  

}
 */
 
}
