#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <math.h>

#include <cstdio>
#define SSCANF std::sscanf

//#include <limits.h>

#include "MaterialData.hh"
#include "LoadMaterialData.hh"
#include "LoadMaterialDataFile.hh"
#include "WriteInfo.hh"
#include "General/environment.hh"
#include "DataInOut/ParamHandling/BaseParamHandler.hh"

namespace CoupledField
{


  // ***************
  //   Constructor
  // ***************
  LoadMaterialDataFile::LoadMaterialDataFile( const char * aFilename )
    : filename(aFilename), scaleMatDat(0) {
    ENTER_FCN( "LoadMaterialDataFile::LoadMaterialDataFile" );
  }


  // ************************
  //   LoadMaterialDataFile
  // ************************

  // method for loading  material information from file "filename"

  void LoadMaterialDataFile::GetMaterial( MaterialData& material,
                                          const std::string matName,
                                          const std::string matType ) {

    ENTER_FCN( "LoadMaterialDataFile::GetMaterial" );

    char buffer[bufLength];

    // Open data file and check for errors
    std::ifstream fin( filename );
    Boolean matC = FALSE;

    if ( !fin.good() ) {
      std::cerr << "File " << filename << " does not exist!" << std::endl;
      exit(1);
    }

#ifdef DEBUG
    (*debug) << std::endl << "*** Load material data of " << matName 
             << " from type " << matType << " from file " << filename
             << std::endl;
#endif

    FindMat( fin, matName.c_str(), buffer, matType.c_str() );

    // first line of material record: matNr. matType matName (nonLin)
    // SSCANF(buffer,"%*d%s", charMatType );
    SSCANF(buffer,"%*d%*s" );
   
    if ( matType == "piezo" ) {

      // read real parts of piezo Material
      ReadPiezo( fin, &material, matC );

      if( params->HasValue( "type", "imagMaterialParameter", "materialDataType" ) ) {

        matC = TRUE;
        std::string matNameImag = matName + "-imag";

        // std::strcat( charMatName, "-imag" );

        FindMat( fin, matNameImag.c_str(), buffer, matType.c_str() );

        // SSCANF( buffer, "%*d%s", charMatType );
        SSCANF( buffer, "%*d%*s" );

         // Read imaginary parts of piezo Material
        ReadPiezo( fin, &material, matC );
      }

    }
    else if ( matType == "fluid" ) {
      ReadFluid(fin, &material);
    }
    else if ( matType == "magnetic" ) {
      ReadMagnetic(fin, &material);
    }
    else if ( matType == "thermic" ) {
      ReadThermic(fin, &material);
    }
    else {
      (*error) << "Warning: material type " << matType << " in File "
               << filename << " unknown!";
      Error( __FILE__, __LINE__ );
    }

    // Close data file
    fin.close();

  }


  // ************
  //   ReadLine
  // ************

  // read next line - ignore lines with char # in it
  void LoadMaterialDataFile::ReadLine( std::ifstream & fin, char* buffer ) {

    ENTER_FCN( "LoadMaterialDataFile::ReadLine" );

    Boolean found = FALSE;

    while ( !found && !fin.eof() ) {
      fin.getline( buffer, bufLength, '\n' );
      if ( strchr(buffer,'#') == NULL ) {
        found = TRUE;
        if ( fin.eof() ) {
          buffer = NULL;
          //std::cout << std::endl << " ReadLine: Unexpected end of file! "
          // << std::endl;
          //exit(EXIT_FAILURE);
        }
      }
    }
  }


  // ***********
  //   FindMat
  // ***********
  void LoadMaterialDataFile::FindMat( std::ifstream &fin, const char* matName,
                                      char* buffer, const char* matType ) {

    ENTER_FCN( "LoadMaterialDataFile::FindMat" );

    std::string errMsg;

    Boolean found = FALSE;
    Integer pos;
    char tempMatName[bufLength];
    char tempMatType[bufLength];

    // set reading position to the beginning of the file
    fin.seekg(0, std::ios::beg);

    while (!found && !fin.eof() ) {

      pos = fin.tellg();
      ReadLine(fin, buffer);
      SSCANF(buffer,"%*d%s",tempMatType);
      SSCANF(buffer,"%*d%*s%s",tempMatName);

      if ( strcmp(matName,tempMatName) == 0 &&
           strcmp(matType,tempMatType) == 0 ) {
        found = TRUE;
        fin.seekg(pos, std::ios::beg);
      }
    }

    if (!found) {       
      errMsg  = "FindMat: Material '";
      errMsg += matName;
      errMsg += "' not found in file '";
      errMsg += filename;
      errMsg += "'!";
      Error(errMsg.c_str(), __FILE__, __LINE__);
    }       
  }


  void LoadMaterialDataFile::ReadPiezo(std::ifstream &fin,
                                       MaterialData *material,
                                       Boolean &matC ) {

    ENTER_FCN( "LoadMaterialDataFile::ReadPiezo" );

    UInt i,j;
    Double helpval;
    Double alfa,beta;
    Double density;
    char buffer[bufLength];
    char materialName[bufLength];
    std::istringstream * strPtr;

    if (matC==FALSE)
      material->DefFull3dMatrix(); // declare matrix with 9x9 entries

    char nonLin[bufLength] = "no";
    std::string hystType = "no"; //char hystType[bufLength] = "no";
    char htype[bufLength];
    ReadLine(fin,buffer);
    if ( SSCANF(buffer,"%*d%*s%s%s%s", materialName, nonLin, htype) == 3 ) {
      hystType = htype;
    }

    material -> SetName(materialName);

    // read stiffness matrix
    for (i=1; i<=6; i++) {
      ReadLine(fin,buffer);
      strPtr = new std::istringstream(buffer);
      
      for (j=1 ; j<=6; j++) {
        *strPtr >> helpval;
        if (strPtr->fail())
          std::cout << "*** The materialfile is corrupt! ***  Material: "
                    << materialName << std::endl;

        // material -> SetPiezoMatrixData(i,j, helpval);

        if (matC==FALSE)
          material -> SetPiezoMatrixData(i-1, j-1, helpval);
        else if (matC==TRUE)
          material -> SetPiezoMatrixDataC(i-1, j-1, helpval);
      }
      delete strPtr;    
    }

    // read piezoelectric coupling terms
    for ( i = 1; i <= 3; i++ ) {
      ReadLine(fin,buffer);
      strPtr = new std::istringstream(buffer);
      
      for ( j = 1; j <= 6; j++ ) {
        *strPtr >> helpval;
        if (strPtr->fail())
          std::cout << "*** The materialfile is corrupt! ***  Material: "
                    << materialName << std::endl;

        if ( matC == FALSE ) {
          material -> SetPiezoMatrixData(i+6-1,j-1, helpval);
          // write transposed coupling terms
          material -> SetPiezoMatrixData(j-1,i+6-1, helpval);
          // material -> SetPiezoMatrixData(i+6,j, helpval);
          // write transposed coupling terms
          // material -> SetPiezoMatrixData(j,i+6, helpval);
        }
        else if (matC==TRUE){
          material -> SetPiezoMatrixDataC(i+6-1,j-1, helpval);
          // write transposed coupling terms
          material -> SetPiezoMatrixDataC(j-1,i+6-1, helpval);
        }
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
              std::cout << "*** The materialfile is corrupt! ***  Material: "
                        << materialName << std::endl;

            if (matC==FALSE){
              //      material -> SetPiezoMatrixData(i+6-1,j+6-1, helpval);
              // indizes of Matrix<Double> start at 0 !!!!!!!!!!!!!!!!!
              material -> SetPiezoMatrixData(i+6-1,j+6-1, helpval);
            }
            else if (matC==TRUE){
              material -> SetPiezoMatrixDataC(i+6-1,j+6-1, helpval);
            }
          }
        delete strPtr;  
      }

    // Rotation of the MaterialMatrices corresponding to the polarisation
    //  of the piezoelectric body

    Double a1, a2, a3;
    a1=a2=a3=0; 

    if( params->HasValue( "x", "1", "piezo", "polingDirection" ) )
      a1=1;
 
    if( params->HasValue( "y", "1", "piezo", "polingDirection" ) ){
      if (params->HasValue("subtype", "axi", "piezo") ){
        std::cout<< "\n Be aware, that you are treating an axisymmetric "
                 << "piezoelectric body, which does not have any y-direction."
                 << "\n Please check your xml-file. \n Press Ctrl+C to stop "
                 << "calculation, any other key to continue."
                 << std::endl;
        getchar();
      }
      a2 = 1;
    }
 
    if( params->HasValue( "z", "1", "piezo", "polingDirection" ) ) {
      a3 = 1;
    }
 
    // if no poling direction is specified,
    // the z-direction is chosen as default
    if ( a1 == 0 && a2 == 0 && a3 == 0 ) {
      a3 = 1.0;
    }

    material->RotateMaterialMatrix(a1,a2,a3);

    ReadLine(fin,buffer);
    strPtr = new std::istringstream(buffer);
    *strPtr >>  density >> alfa >> beta;
    if (strPtr->fail()) {
      std::cout << "*** The materialfile is corrupt! ***  Material: "
                << materialName << std::endl;
    }
    delete strPtr;

    if ( strcmp(nonLin,"hysteresis:") == 0 ) {
      Double Esat, Psat;
      Double aJiles, alphaJiles, kJiles, cJiles;
      UInt dirPol;
      ReadLine(fin,buffer);
      strPtr = new std::istringstream(buffer);
      *strPtr >>  Esat >> Psat >> dirPol;

      if (strPtr->fail()) {
        (*error) << "*** The materialfile is corrupt (hysteresis)! "
                 << "***  Material: " << materialName;
        Error( __FILE__, __LINE__ );
      }

      delete strPtr;
      ReadLine(fin,buffer);
      strPtr = new std::istringstream(buffer);

      if (hystType == "jiles") {
        *strPtr >> aJiles >> alphaJiles >> kJiles >> cJiles;
      }


      //check for correct direction
      std::string probGeo;
      params->Get( "type", probGeo, "geometry" );
      if ( (probGeo == "axi") || ( probGeo == "plane") ) {
        if (dirPol == 3) {
          dirPol = 2;
        }
      }

      if (hystType == "preisach") {
        material -> SetEsat(Esat);
      }
      else if ( hystType == "jiles" ) {
        material -> SetJiles_a(aJiles);
        material -> SetJiles_alpha(alphaJiles);
        material -> SetJiles_k(kJiles);
        material -> SetJiles_c(cJiles);
      }

      material -> SetPsat(Psat);
      material -> SetDirPol(dirPol);
      material -> SetHysteresisType(hystType);

      delete strPtr;    
    }



    // ================ SCALING OF PIEZOELECTRIC MATRICES =====================
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
    

    if (scaleMatDat) {
      Info->PrintF( "", "\n!!!!!! SCALING with Diag(1e-5, 1e-5, 1e-5, 1e-5,"
                    "1e-5, 1e-5,1e5, 1e5, 1e5) IS ON !!!!!\n\n ");
    }
  }


  void LoadMaterialDataFile::ReadFluid( std::ifstream &fin,
                                        MaterialData *material ) {
    ENTER_FCN("LoadMaterialDataFile::ReadFluid");

    Double density=0.0, compress=0.0, alfa=0.0, beta=0.0, BoverA=0.0;   
    std::istringstream * strPtr;
    char buffer[bufLength];
    char materialName[bufLength];


    ReadLine(fin,buffer);
    SSCANF(buffer,"%*d%*s%s", materialName);  

    material -> SetName(materialName);

    ReadLine(fin,buffer);
    strPtr = new std::istringstream(buffer);

    *strPtr >> density >> compress >> alfa >> beta >> BoverA;
    if (strPtr->fail()) {
      std::cout << "*** The materialfile is corrupt! ***  Material: " 
                << materialName << std::endl
                << "Please specify density, compression module, alfa, beta "
                << "and BoverA."
                << std::endl;
    }

    delete strPtr;
    
    material->SetCompressibility(compress);
    material->SetDensity(density);
    material->SetDampingCoeffs(alfa,beta);
    material->SetBoverA(BoverA);

    Info->PrintFluidMat(*material);
  }

  void LoadMaterialDataFile::ReadThermic( std::ifstream &fin,
                                          MaterialData *material ) {
    ENTER_FCN("LoadMaterialDataFile::ReadThermic");

    Double density=0.0, heatCapacity=0.0, thermalConductivity=0.0;
    std::istringstream * strPtr;
    char buffer[bufLength];
    char materialName[bufLength];


    ReadLine(fin,buffer);
    SSCANF(buffer,"%*d%*s%s", materialName);  

    material -> SetName(materialName);

    ReadLine(fin,buffer);
    strPtr = new std::istringstream(buffer);

    *strPtr >> density >> heatCapacity >> thermalConductivity;
    if (strPtr->fail()) {
      std::cout << "*** The materialfile is corrupt! ***  Material: " 
                << materialName << std::endl;
    }

    delete strPtr;
    
    material->SetDensity(density);
    material->SetThermic(heatCapacity,thermalConductivity);

    Info->PrintThermicMat(*material);
  }


  void LoadMaterialDataFile::ReadMagnetic( std::ifstream &fin,
                                           MaterialData *material ) {

    ENTER_FCN( "LoadMaterialDataFile::ReadMagnetic" );

    Double mX, mY, mZ;
    Double conductivity, permeability;

    std::istringstream * strPtr;
    char buffer[bufLength];
    char materialName[bufLength];
    char aux1[bufLength];
    char aux2[bufLength];

    ReadLine( fin, buffer );
    Integer numRead = SSCANF( buffer, "%*d%*s%s%s%s", materialName, aux1,aux2);
    if ( numRead > 1 ) {
      if ( strcmp(aux1,"bhapprox:") == 0 ) {
        material->SetBHCurveFileName( aux2 );
      }
    }

    material -> SetName(materialName);

    ReadLine(fin,buffer);
    strPtr = new std::istringstream(buffer);
      
    *strPtr >> conductivity >> permeability >> mX >> mY >> mZ;
    if (strPtr->fail()) {
      (*error) << "LoadMaterialDataFile::ReadMagnetic: The materialfile is "
               << "corrupt!\n";
      Error( __FILE__, __LINE__ );
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
