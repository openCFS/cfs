#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <math.h>

// Temporary hack for avoiding namespace conflict for old C-functions and old SGI compiler
#ifdef __sgi
#include <stdio.h>
#define SSCANF sscanf
#else
#include <cstdio>
#define SSCANF std::sscanf
#endif

//#include <limits.h>

#include "MaterialData.hh"
#include "LoadMaterialData.hh"
#include "LoadMaterialDataFile.hh"
#include "WriteInfo.hh"
#include "General/environment.hh"

namespace CoupledField
{

  LoadMaterialDataFile::LoadMaterialDataFile (const char * aFilename)
    :filename(aFilename),scaleMatDat(0)
  {
    ENTER_FCN("LoadMaterialDataFile::LoadMaterialDataFile");
  }


  // load material information from file "filename"
  void LoadMaterialDataFile::GetMaterial( MaterialData& material, const std::string matName, const std::string matType)
  {
    ENTER_FCN("LoadMaterialDataFile::GetMaterial");
    char buffer[bufLength];  
    
    // 
    char * charFileName = c_string(filename);
    char * charMatName = c_string(matName);
    char * charMatType = c_string(matType);    
    
    std::ifstream fin(charFileName);

    if (!fin.good())
      {
	std::cerr << "File " << filename << " does not exist!" << std::endl;
	exit(1);
      }

#ifdef DEBUG
    (*debug) << std::endl << "*** Load material data of " << matName 
	     << " from type " << matType << " from file " << filename << std::endl;
#endif

    
    FindMat(fin, charMatName, buffer, charMatType);

    // first line of material record: matNr. matType matName (nonLin)
    SSCANF(buffer,"%*d%s", charMatType);

	/*
	  if (strcmp(charMatType,"magnonlin") == 0)
	  {
	  material = new MaterialData;
	  ReadMagNonLin(fin, material);
	  }
	*/

    
	if (strcmp(charMatType,"piezo") == 0 )
	  {
	    ReadPiezo(fin, &material);	 
	
	    /*
	      if (eulerAngles.size())
	      {
	      EulerAnglesRotate(&material, eulerAngles[i]);
	    
	      std::cout << "after EULER ROTATION : " << std::endl
	      << "LoadMaterialDataFile::LoadMaterial: gesamte Piezo-Datenmatrix von " << matName
	      << ":" << std::endl << *material.GetMatrix() << std::endl << std::endl
	      << "density = " << material.GetDensity() << std::endl
	      << "damping coefficient alfa = " << material.GetDampingAlfa() << std::endl
	      << "damping coefficient beta = " << material.GetDampingBeta() << std::endl;

	      if (scaleMatDat)
	      std::cout << std::endl << "!!!!!! SCALING with Diag(1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5, "
	      << "1e5, 1e5, 1e5) IS ON !!!!! " << std::endl << std::endl;	   
	      }
	    */
	  }
    
	else if (strcmp(charMatType,"fluid") == 0)
	  ReadFluid(fin, &material);
    
	else if (strcmp(charMatType,"magnetic") == 0 )
	  ReadMagnetic(fin, &material);
    
	else 
	  {
	    std::cerr << "Warning: materialtype " << charMatType << " in File " << filename << " unknown!" << std::endl;
	    exit(EXIT_FAILURE);
	  }
    
    fin.close();
  }


  // read next line - ignor lines with char # in it
  void LoadMaterialDataFile::ReadLine(std::ifstream & fin, char* buffer)
  {
    ENTER_FCN("LoadMaterialDataFile::ReadLine");
    Integer found = 0;

    while (!found && !fin.eof())
      {
	fin.getline(buffer,bufLength,'\n');
	if ( strchr(buffer,'#') == NULL )
	  found = 1;
	if (fin.eof())
	  {
	    buffer = NULL;
	    //std::cout << std::endl << " ReadLine: Unexpected end of file! " << std::endl;
	    //exit(EXIT_FAILURE);
	  }
      }
  }



  void LoadMaterialDataFile :: FindMat(std::ifstream & fin, const char* matName, char* buffer, char* matType)
  {
    ENTER_FCN("LoadMaterialDataFile::FindMat");

    std::string errMsg;

    Integer found = 0;
    Integer pos;
    char tempMatName[bufLength];
    char tempMatType[bufLength];

    // set reading position to the beginning of the file
    fin.seekg(0, std::ios::beg);

    while (!found && !fin.eof() )
      {
	pos = fin.tellg();
	ReadLine(fin, buffer);
	SSCANF(buffer,"%*d%s",tempMatType);
	SSCANF(buffer,"%*d%*s%s",tempMatName);

	if ( strcmp(matName,tempMatName) == 0 &&  
	     strcmp(matType, tempMatType) == 0)
	  {
	    found = 1;
	    fin.seekg(pos, std::ios::beg);
	  }
      }

    if (!found)
      {	
	errMsg  = "FindMat: Material '";
	errMsg += matName;
	errMsg += "' not found in file '";
	errMsg += filename;
	errMsg += "'!";
	Error(errMsg.c_str(), __FILE__, __LINE__);
      }       
  }


  void LoadMaterialDataFile :: ReadPiezo(std::ifstream & fin, MaterialData * material)
  {
    ENTER_FCN("LoadMaterialDataFile::ReadPiezo");
    Integer i,j;
    Double helpval;
    Double alfa,beta;
    Double density;
    char buffer[bufLength];
    char materialName[bufLength];
    std::istringstream * strPtr;

    material->DefFull3dMatrix(); // declare matrix with 9x9 entries

    ReadLine(fin,buffer);
    SSCANF(buffer,"%*d%*s%s", materialName);  

    material -> SetName(materialName);

    // read stiffness matrix
    for (i=1; i<=6; i++)
      {
	ReadLine(fin,buffer);
	strPtr = new std::istringstream(buffer);
      
	for (j=1 ; j<=6; j++)
	  {
	    *strPtr >> helpval;
	    if (strPtr->fail())
	      std::cout << "*** The materialfile is corrupt! ***  Material: " << materialName << std::endl;

	    //	    material -> SetPiezoMatrixData(i,j, helpval);
	    material -> SetPiezoMatrixData(i-1, j-1, helpval);
	  }
	delete strPtr;	
      }
    
    

    // read piezoelectric coupling terms
    for (i=1; i<=3; i++)
      {
	ReadLine(fin,buffer);
	strPtr = new std::istringstream(buffer);
      
	for (j=1 ; j<=6; j++)
	  {
	    *strPtr >> helpval;
	    if (strPtr->fail())
	      std::cout << "*** The materialfile is corrupt! ***  Material: " << materialName << std::endl;

	    material -> SetPiezoMatrixData(i+6-1,j-1, helpval);
	    material -> SetPiezoMatrixData(j-1,i+6-1, helpval);     // writes transposed coupling terms
	    //	    material -> SetPiezoMatrixData(i+6,j, helpval);
	    //	    material -> SetPiezoMatrixData(j,i+6, helpval);     // writes transposed coupling terms
	  }
	delete strPtr;	
      }


    // read dielectric terms
    for (i=1; i<=3; i++)
      {
	ReadLine(fin,buffer);
	strPtr = new std::istringstream(buffer);
      
	for (j=1 ; j<=3; j++)
	  {
	    *strPtr >> helpval;
	    if (strPtr->fail())
	      std::cout << "*** The materialfile is corrupt! ***  Material: " << materialName << std::endl;

	    //	    material -> SetPiezoMatrixData(i+6-1,j+6-1, helpval);
	    // indizes of Matrix<Double> start at 0 !!!!!!!!!!!!!!!!!
	    material -> SetPiezoMatrixData(i+6-1,j+6-1, helpval);
	  }
	delete strPtr;	
      }

    ReadLine(fin,buffer);
    strPtr = new std::istringstream(buffer);
    *strPtr >>  density >> alfa >> beta;
    if (strPtr->fail())
      std::cout << "*** The materialfile is corrupt! ***  Material: " << materialName << std::endl;

    delete strPtr;	

    // ================ SCALING OF PIEZOELECTRIC MATRICES =====================================
    /*
    // bring order of all data elems to 1 ==>
    //            | 1e-5  0 |   | c  e^T |   | 1e-5  0 |
    // dMatNew =  |         | * |        | * |         |
    //            | 0   1e5 |   | e  eps |   | 0   1e5 |
    // additional scaling of density is necessarry !!!!
    if (scaleMatDat)
    {
    Matrix<Double> * tmpMatDat = material->GetMatrix();
    Integer nr3d = material->GetNrElems3d();

    Matrix<Double>   scaleMat( nr3d, nr3d);
    Matrix<Double>   helpMat(  nr3d, nr3d);
    scaleMat.Init();

    for(Integer i=1; i <= 6; i++)
    scaleMat(i,i) = 1e-5;

    for(Integer i=7; i <= nr3d; i++)
    scaleMat(i,i) = 1e5;

    Mult (scaleMat, *tmpMatDat, helpMat);
    Mult(helpMat, scaleMat, *tmpMatDat);

    density *= 1e-10;

    material -> SetScaledFlag();
    }
    */

    Double lambda, mu;
    //lambda = c_12 
    material->GetPiezoMatrixData(0,1,lambda);

    //mu = (c_11-c_12)/2 
    material->GetPiezoMatrixData(0,0,mu);   
    mu = 0.5*(mu-lambda);
    
    material->SetLambda(lambda);
    material->SetMu(mu);
    

    material -> SetDensity(density);
    material -> SetDampingCoeffs(alfa,beta);
    
    Info->PrintPiezoMat(*material);
    

    if (scaleMatDat)
      Info->PrintF("","\n!!!!!! SCALING with Diag(1e-5, 1e-5, 1e-5, 1e-5, 1e-5, 1e-5,1e5, 1e5, 1e5) IS ON !!!!!\n\n ");
    
  }


  void LoadMaterialDataFile :: ReadFluid(std::ifstream & fin, MaterialData * material)
  {
    ENTER_FCN("LoadMaterialDataFile::ReadFluid");
    Double alfa,beta;
    Double density, compress;   
    Double BoverA;
    std::istringstream * strPtr;
    char buffer[bufLength];
    char materialName[bufLength];

    ReadLine(fin,buffer);
    SSCANF(buffer,"%*d%*s%s", materialName);  

    material -> SetName(materialName);

    ReadLine(fin,buffer);
    strPtr = new std::istringstream(buffer);
      
    *strPtr >> density >> compress >> alfa >> beta >> BoverA;
    if (strPtr->fail())
      std::cout << "*** The materialfile is corrupt! ***  Material: " 
				<< materialName << std::endl
				<< "Please specify density, compression module, alfa, beta and BoverA."
				<< std::endl;
	
    delete strPtr;
    
    material->SetCompressibility(compress);
    material->SetDensity(density);
    material->SetDampingCoeffs(alfa,beta);
    material->SetBoverA(BoverA);

    Info->PrintFluidMat(*material);
  }


  void LoadMaterialDataFile :: ReadMagnetic(std::ifstream & fin, MaterialData * material)
  {
    ENTER_FCN("LoadMaterialDataFile::ReadMagnetic");
    Double mX, mY, mZ;
    Double conductivity, permeability;
    std::string errMsg;

    std::istringstream * strPtr;
    char buffer[bufLength];
    char materialName[bufLength];
    char aux1[bufLength];
    char aux2[bufLength];

    ReadLine(fin,buffer);
    SSCANF( buffer, "%*d%*s%s%s%s", materialName, aux1, aux2 );
    if ( strcmp(aux1,"bhapprox:") == 0 ) {
      material->SetBHCurveFileName( aux2 );
    }

    material -> SetName(materialName);

    ReadLine(fin,buffer);
    strPtr = new std::istringstream(buffer);
      
    *strPtr >> conductivity >> permeability >> mX >> mY >> mZ;
    if (strPtr->fail()) {
      errMsg  = "LoadMaterialDataFile::ReadMagnetic:  The materialfile is corrupt!\n";
      errMsg += "Material: ";
      errMsg +=  materialName;
      Error(errMsg.c_str(), __FILE__, __LINE__);
    }
    
    delete strPtr;
   
    for (int i=0; i<3; i++)
    { 
      material->SetPermeability(i,i,permeability);
      material->SetConductivity(i,i,conductivity);
    }
    material->SetPermMag(mX, mY, mZ);

    Info->PrintMagMat(*material);
  }


} // end namespace CoupledField


















  /*
    void LoadMaterialDataFile :: EulerAnglesRotate(MaterialData * material, const Vector<Double> & eulerAngles)
    {
    #ifdef TRACE
    if (trace) (*trace) << "Entering  LoadMaterialDataFile::EulerAnglesRotate" << std::endl;
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
  void LoadMaterialDataFile :: ReadMagNonLin(std::ifstream & fin, MaterialData * material)
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
  std::cout << std::endl << "LoadMaterialDataFile: There must be three Euler Angles per material! " << std::endl;
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


#undef SSCANF
