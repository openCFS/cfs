#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <fstream>
#include <math.h>

#include <cstdio>
#define SSCANF std::sscanf

#include "PlainMaterialHandler.hh"
#include "WriteInfo.hh"
#include "General/environment.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

// material header
#include "Materials/electroMagneticMaterial.hh"
#include "Materials/electrostaticMaterial.hh"
#include "Materials/heatMaterial.hh"
#include "Materials/acousticMaterial.hh"
#include "Materials/mechanicMaterial.hh"
#include "Materials/piezoMaterial.hh"
#include "Materials/flowMaterial.hh"

namespace CoupledField
{


  // ***************
  //   Constructor
  // ***************
  PlainMaterialHandler::PlainMaterialHandler( const std::string & fileName )
    : MaterialHandler( fileName ) {
    ENTER_FCN( "PlainMaterialHandler::PlainMaterialHandler" );

    scaleMatDat = 0;
  }


  // ************************
  //   PlainMaterialHandler
  // ************************

  // method for loading  material information from file "filename"

  BaseMaterial * PlainMaterialHandler::
  LoadMaterial( const std::string matName,
                const MaterialClass matType ) {

    ENTER_FCN( "PlainMaterialHandler::GetMaterial" );
    // Open data file and check for errors
    std::ifstream fin( fileName_.c_str() );
    bool matC = false;
    bool hasImagData_ = true;
    

    if ( !fin.good() ) {
      std::cerr << "File " << fileName_ << " does not exist!" << std::endl;
      exit(1);
    }

#ifdef DEBUG
    (*debug) << std::endl << "*** Load material data of " << matName 
             << " from type " << matType << " from file " << fileName_
             << std::endl;
#endif

    //to be conform to old material file
    BaseMaterial * material = NULL;
    std::string matTypeOld;
    if ( matType == PIEZO ) {
      matTypeOld = "piezo";
      material = new PiezoMaterial();
    }
    else if ( matType == MECHANIC ) {
      matTypeOld = "piezo";
      material = new MechanicMaterial();
    }    
    else if ( matType == FLUID ) {
      matTypeOld = "fluid";
      material = new AcousticMaterial();
    }
    else if ( matType == ELECTROMAGNETIC ) {
      matTypeOld = "magnetic";
      material = new ElectroMagneticMaterial();
    }
    else if ( matType == ELECTROSTATIC ) {
      matTypeOld = "piezo";
      material = new ElectroStaticMaterial();
    }
    else if ( matType == THERMIC ) {
      matTypeOld = "thermic";
      material = new HeatMaterial();
    }
    else if ( matType == FLOW ) {
      matTypeOld = "fluid";
      material = new FlowMaterial();
    }

    char buffer[bufLength];
    if( !FindMat( fin, matName.c_str(), buffer, matTypeOld.c_str() ) ) {
      EXCEPTION( "Could not find material '" << matName <<"' of type "
                 << matTypeOld );
    }

    // first line of material record: matNr. matType matName (nonLin)
    // SSCANF(buffer,"%*d%s", charMatType );
    SSCANF(buffer,"%*d%*s" );
   

    if ( matType == PIEZO ) {

      // read real parts of piezo Material
      ReadPiezo( fin, material, matC );

      //read imaginary part
      matC = true;
      std::string matNameImag = matName + "-imag";
      fin.seekg(0,std::ios::end);
      hasImagData_ = FindMat( fin, matNameImag.c_str(), buffer, matTypeOld.c_str() );
      if( hasImagData_ ) {
        // SSCANF( buffer, "%*d%s", charMatType );
        SSCANF( buffer, "%*d%*s" );
        // Read imaginary parts of piezo Material
        ReadPiezo( fin, material, matC );
      }

    }
    else if ( matType == MECHANIC ) {
      ReadMechanic(fin, material, matC);

      //read imaginary part
      matC = true;
      std::string matNameImag = matName + "-imag";
      fin.seekg(0,std::ios::end);
      hasImagData_ =  FindMat( fin, matNameImag.c_str(), buffer, matTypeOld.c_str() );
      if( hasImagData_ ) {
        // SSCANF( buffer, "%*d%s", charMatType );
        SSCANF( buffer, "%*d%*s" );
        // Read imaginary parts of piezo Material
        ReadMechanic( fin, material, matC );
      }
    }    
    else if ( matType == FLUID ) {
      ReadAcoustic(fin, material);
    }
    else if ( matType == ELECTROMAGNETIC ) {
      ReadMagnetic(fin, material);
    }
    else if ( matType == ELECTROSTATIC ) {
      ReadElectrostatic(fin, material, matC);

      //read imaginary part
      matC = true;
      std::string matNameImag = matName + "-imag";
      fin.seekg(0,std::ios::end);
      hasImagData_ = FindMat( fin, matNameImag.c_str(), buffer, matTypeOld.c_str() );
      if( hasImagData_ ) {
        // SSCANF( buffer, "%*d%s", charMatType );
        SSCANF( buffer, "%*d%*s" );
        // Read imaginary parts of piezo Material
        ReadElectrostatic( fin, material, matC );
      }
    }
    else if ( matType == THERMIC ) {
      ReadThermic(fin, material);
    }
    else if ( matType == FLOW ) {
      ReadFlow(fin, material);
    }
    else {
      (*error) << "Warning: material type " << matType << " in File "
               << fileName_ << " unknown!";
      Error( __FILE__, __LINE__ );
    }

    // Close data file
    fin.close();

    // Return material data
    return material;
  }


  // ************
  //   ReadLine
  // ************

  // read next line - ignore lines with char # in it
  void PlainMaterialHandler::ReadLine( std::ifstream & fin, char* buffer ) {

    ENTER_FCN( "PlainMaterialHandler::ReadLine" );

    bool found = false;

    while ( !found && !fin.eof() ) {
      fin.getline( buffer, bufLength, '\n' );
      if ( strchr(buffer,'#') == NULL ) {
        found = true;
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
  bool PlainMaterialHandler::FindMat( std::ifstream &fin, const char* matName,
                                      char* buffer, const char* matType ) {
    
    ENTER_FCN( "PlainMaterialHandler::FindMat" );
    
    std::string errMsg;
    
    bool found = false;
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
        found = true;
        fin.seekg(pos, std::ios::beg);
      }
    }

    return found;
  }


  void PlainMaterialHandler::ReadPiezo(std::ifstream &fin,
                                       BaseMaterial *material,
                                       bool &matC ) {

    ENTER_FCN( "PlainMaterialHandler::ReadPiezo" );

    UInt i,j;
    Double helpval;
    char buffer[bufLength];
    char materialName[bufLength];
    std::istringstream * strPtr;

    char nonLin[bufLength] = "no";
    std::string hystType = "no"; //char hystType[bufLength] = "no";
    char htype[bufLength];
    ReadLine(fin,buffer);

    if ( SSCANF(buffer,"%*d%*s%s%s%s", materialName, nonLin, htype) == 3 ) {
      hystType = htype;
    }

    // read stiffness matrix not necessary for direct coupled PDE
    Matrix<Double> realC;
    realC.Resize(6,6);
    realC.Init();

    for (i=1; i<=6; i++) {
      ReadLine(fin,buffer);
      strPtr = new std::istringstream(buffer);
      delete strPtr;    
    }

    // read piezoelectric coupling terms
    Matrix<Double> realP;
    realP.Resize(3,6);
    realP.Init();

    for ( i = 1; i <= 3; i++ ) {
      ReadLine(fin,buffer);
      strPtr = new std::istringstream(buffer);
      
      for ( j = 1; j <= 6; j++ ) {
        *strPtr >> helpval;
        if (strPtr->fail())
          std::cout << "*** The materialfile is corrupt! ***  Material: "
                    << materialName << std::endl;

        realP[i-1][j-1] = helpval;
      }
      delete strPtr;            
    }  

    if ( matC ) {
      material->SetTensor(realP,PIEZO_TENSOR,IMAG);
    }
    else {
      material->SetTensor(realP,PIEZO_TENSOR,REAL);
    }

    //currently not supported
    // Rotation of the MaterialMatrices corresponding to the polarisation
    //  of the piezoelectric body
    //     Double a1, a2, a3;
    //     a1=a2=a3=0; 


    //     if( params->HasValue( "x", "1", "piezo", "polingDirection" ) )
    //       a1=1;
 
    //     if( params->HasValue( "y", "1", "piezo", "polingDirection" ) ){
    //       if (params->HasValue("subtype", "axi", "piezo") ){
    //         std::cout<< "\n Be aware, that you are treating an axisymmetric "
    //                  << "piezoelectric body, which does not have any y-direction."
    //                  << "\n Please check your xml-file. \n Press Ctrl+C to stop "
    //                  << "calculation, any other key to continue."
    //                  << std::endl;
    //         getchar();
    //       }
    //       a2 = 1;
    //     }
 
    //     if( params->HasValue( "z", "1", "piezo", "polingDirection" ) ) {
    //       a3 = 1;
    //     }
 
    //     // if no poling direction is specified,
    //     // the z-direction is chosen as default
    //     if ( a1 == 0 && a2 == 0 && a3 == 0 ) {
    //       a3 = 1.0;
    //     }

    //     material->RotateMaterialMatrix(a1,a2,a3);

    // not necessary for direct coupled piezo PDE
    //     ReadLine(fin,buffer);
    //     strPtr = new std::istringstream(buffer);
    //     *strPtr >>  density >> alfa >> beta;
    //     if (strPtr->fail()) {
    //       std::cout << "*** The materialfile is corrupt! ***  Material: "
    //                 << materialName << std::endl;
    //     }
    //     delete strPtr;


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
      param->Get("domain")->Get("geometryType", probGeo);
      if ( (probGeo == "axi") || ( probGeo == "plane") ) {
        if (dirPol == 3) {
          dirPol = 2;
        }
      }

      //       if (hystType == "preisach") {
      //         material -> SetEsat(Esat);
      //       }
      //       else if ( hystType == "jiles" ) {

      if ( hystType == "jiles" ) {
        //         material -> SetJiles_a(aJiles);
        //         material -> SetJiles_alpha(alphaJiles);
        //         material -> SetJiles_k(kJiles);
        //         material -> SetJiles_c(cJiles);

	material->SetScalar(aJiles,A_JILES,REAL);
	material->SetScalar(alphaJiles,ALPHA_JILES,REAL);
	material->SetScalar(kJiles,K_JILES,REAL);
	material->SetScalar(cJiles,C_JILES,REAL);
      }

      //       material -> SetPsat(Psat);
      //       material -> SetDirPol(dirPol);
      //       material -> SetHysteresisType(hystType);

      material->SetScalar(Psat,P_SATURATION,REAL);
      material->SetScalar(Esat,E_SATURATION,REAL);
      //      material->SetScalar((Integer&)dirPol,P_DIRECTION);
      material->SetScalar(hystType,HYST_MODEL);

      delete strPtr;    
    }


    //     Double lambda, mu;
    //     //lambda = c_12 
    //     material->GetPiezoMatrixData(0,1,lambda);

    //     //mu = (c_11-c_12)/2 
    //     material->GetPiezoMatrixData(0,0,mu);   
    //     mu = 0.5*(mu-lambda);
    
    //     material->SetLambda(lambda);
    //     material->SetMu(mu);
    

    //     material -> SetDensity(density);
    //     material -> SetDampingCoeffs(alfa,beta);
    

    Info->PrintMaterial(material);

    //just for test
//     StdVector<Double> phi(3);
//     phi[0] = 45;
//     phi[1] = 0;
//     phi[2] = 45;
//     material->RotateTensorByRotationAngles( phi, PIEZO_TENSOR);

    Info->PrintMaterial(material);


    if (scaleMatDat) {
      Info->PrintF( "", "\n!!!!!! SCALING with Diag(1e-5, 1e-5, 1e-5, 1e-5,"
                    "1e-5, 1e-5,1e5, 1e5, 1e5) IS ON !!!!!\n\n "); 
    }

  }


  void PlainMaterialHandler::ReadMechanic(std::ifstream &fin,
                                          BaseMaterial* material,
                                          bool &matC ) {

    ENTER_FCN( "PlainMaterialHandler::ReadMechanic" );

    UInt i,j;
    Double helpval;
    Double alfa,beta;
    Double density;
    char buffer[bufLength];
    char materialName[bufLength];
    std::istringstream * strPtr;

    char nonLin[bufLength] = "no";
    char htype[bufLength];
    ReadLine(fin,buffer);
    SSCANF(buffer,"%*d%*s%s%s%s", materialName, nonLin, htype);

    material->SetName(materialName);

    // read stiffness matrix not necessary for direct coupled PDE
    Matrix<Double> realC;
    realC.Resize(6,6);
    realC.Init();

    for (i=1; i<=6; i++) {
      ReadLine(fin,buffer);
      strPtr = new std::istringstream(buffer);
      
      for (j=1 ; j<=6; j++) {
        *strPtr >> helpval;
        if (strPtr->fail())
          std::cout << "*** The materialfile is corrupt! ***  Material: "
                    << materialName << std::endl;

        realC[i-1][j-1] = helpval;
      }
      delete strPtr;    
    }

    if ( matC ) {
      material->SetTensor(realC,MECH_STIFFNESS_TENSOR,IMAG);
    }
    else {
      material->SetTensor(realC,MECH_STIFFNESS_TENSOR,REAL);
    }

    if ( !matC ) {
      // piezoelectric coupling terms
      for ( i = 1; i <= 3; i++ ) {
	ReadLine(fin,buffer);
      }    
      
      // read dielectric terms
      for (i=1; i<=3; i++) {
	ReadLine(fin,buffer);
      }
      
      ReadLine(fin,buffer);
      strPtr = new std::istringstream(buffer);
      *strPtr >>  density >> alfa >> beta;
      if (strPtr->fail()) {
	std::cout << "*** The materialfile is corrupt! ***  Material: "
		  << materialName << std::endl;
      }
      delete strPtr;
      
      material->SetScalar(density,DENSITY,REAL);
      material->SetScalar(alfa,RAYLEIGH_ALPHA,REAL);
      material->SetScalar(beta,RAYLEIGH_BETA,REAL); 
    }

    Info->PrintMaterial(material);
  }

  
  void PlainMaterialHandler::ReadAcoustic( std::ifstream &fin,
                                           BaseMaterial* material ) {
    ENTER_FCN("PlainMaterialHandler::ReadAcoustic");

    Double density=0.0, compress=0.0, alfa=0.0, beta=0.0, BoverA=0.0;   
    Double acousticAlpha=0.0, fracExp=0.0;
    std::istringstream * strPtr;
    char buffer[bufLength];
    char materialName[bufLength];

    ReadLine(fin,buffer);

    SSCANF(buffer,"%*d%*s%s", materialName);  

    material->SetName(materialName);

    ReadLine(fin,buffer);
    strPtr = new std::istringstream(buffer);

    *strPtr >> density >> compress >> alfa >> beta >> BoverA;
    if (strPtr->fail()) {
      std::cout << "*** Material data for acoustics not correct specified ***" 
                << materialName << std::endl
                << std::endl;
    }
    delete strPtr;

    ReadLine(fin,buffer);
    strPtr = new std::istringstream(buffer);

    *strPtr >> acousticAlpha >> fracExp;
    if (strPtr->fail()) {
      std::cout << "*** Material data for acoustics not correct specified ***" 
                << materialName << std::endl
                << std::endl;
    }
    delete strPtr;
    
    material->SetScalar(density, DENSITY, REAL);
    material->SetScalar(compress, ACOU_BULK_MODULUS, REAL);
    material->SetScalar(alfa, RAYLEIGH_ALPHA, REAL);
    material->SetScalar(beta, RAYLEIGH_BETA, REAL);
    material->SetScalar(BoverA, BOVERA, REAL);

    material->SetScalar(acousticAlpha, ACOU_ALPHA, REAL);    
    material->SetScalar(fracExp, FRACTIONAL_EXPONENT, REAL);

    Info->PrintMaterial(material);
  }

  void PlainMaterialHandler::ReadFlow( std::ifstream &fin,
                                       BaseMaterial *material ) {
    ENTER_FCN("PlainMaterialHandler::ReadFlow");

    Double density=0.0, dynViscosity=0.0;
    std::istringstream * strPtr;
    char buffer[bufLength];
    char materialName[bufLength];

    ReadLine(fin,buffer);
    SSCANF(buffer,"%*d%*s%s", materialName);  

    material->SetName(materialName);

    ReadLine(fin,buffer);
    strPtr = new std::istringstream(buffer);

    *strPtr >> density >> dynViscosity;
    if (strPtr->fail()) {
      std::cout << "*** The materialfile is corrupt! ***  Material: " 
                << materialName << std::endl
                << "Please specify density, compression module, alfa, beta "
                << "and BoverA."
                << std::endl;
    }

    delete strPtr;
    
    material->SetScalar(density, DENSITY, REAL);
    material->SetScalar(dynViscosity, DYNAMIC_VISCOSITY, REAL);

    Info->PrintMaterial(material);
  }


  void PlainMaterialHandler::ReadThermic( std::ifstream &fin,
                                          BaseMaterial *material ) {
    ENTER_FCN("PlainMaterialHandler::ReadThermic");

    Double density=0.0, heatCapacity=0.0, thermalConductivity=0.0;
    std::istringstream * strPtr;
    char buffer[bufLength];
    char materialName[bufLength];

    ReadLine(fin,buffer);
    SSCANF(buffer,"%*d%*s%s", materialName);  

    material->SetName(materialName);

    ReadLine(fin,buffer);
    strPtr = new std::istringstream(buffer);

    *strPtr >> density >> heatCapacity >> thermalConductivity;
    if (strPtr->fail()) {
      std::cout << "*** The materialfile is corrupt! ***  Material: " 
                << materialName << std::endl;
    }

    delete strPtr;
    
    //     material->SetDensity(density);
    //     material->SetThermic(heatCapacity,thermalConductivity);

    material->SetScalar(density, DENSITY, REAL);
    material->SetScalar(heatCapacity, HEAT_CAPACITY, REAL);
    material->SetScalar(thermalConductivity, HEAT_CONDUCTIVITY, REAL);
 
    Info->PrintMaterial(material);
  }


  void PlainMaterialHandler::ReadMagnetic( std::ifstream &fin,
                                           BaseMaterial *material ) {

    ENTER_FCN( "PlainMaterialHandler::ReadMagnetic" );

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
	//        material->SetBHCurveFileName( aux2 );
	material->SetNonlinFileName( aux2 );
      }
    }

    material->SetName(materialName);

    ReadLine(fin,buffer);
    strPtr = new std::istringstream(buffer);
      
    *strPtr >> conductivity >> permeability >> mX >> mY >> mZ;
    if (strPtr->fail()) {
      (*error) << "PlainMaterialHandler::ReadMagnetic: The materialfile is "
               << "corrupt!\n";
      Error( __FILE__, __LINE__ );
    }
    
    delete strPtr;
   
    //     for (int i=0; i<3; i++)
    //       { 
    //         material->SetPermeability(i,i,permeability);
    //         material->SetConductivity(i,i,conductivity);
    //       }
    //     material->SetPermMag(mX, mY, mZ);

    material->SetScalar(permeability,MAG_PERMEABILITY, REAL);
    material->SetScalar(conductivity,MAG_CONDUCTIVITY, REAL);

    Info->PrintMaterial(material);
  }


  void PlainMaterialHandler::ReadElectrostatic(std::ifstream &fin,
                                               BaseMaterial *material,
                                               bool &matC ) {

    ENTER_FCN( "PlainMaterialHandler::ReadElectrostatic" );

    UInt i,j;
    Double helpval;
    char buffer[bufLength];
    char materialName[bufLength];
    std::istringstream * strPtr;

    char nonLin[bufLength] = "no";
    std::string hystType = "no"; //char hystType[bufLength] = "no";
    char htype[bufLength];
    ReadLine(fin,buffer);
    if ( SSCANF(buffer,"%*d%*s%s%s%s", materialName, nonLin, htype) == 3 ) {
      hystType = htype;
    }

    material->SetName(materialName);

    // mechanical stiffness tensor
    for ( i = 1; i <= 6; i++ ) {
      ReadLine(fin,buffer);
    }    

    // piezoelectric coupling terms
    for ( i = 1; i <= 3; i++ ) {
      ReadLine(fin,buffer);
    }    

    Matrix<Double> realE(3,3);
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
	    realE[i-1][j-1] = helpval;
          }
        delete strPtr;  
      }

    if ( matC ) {
      material->SetTensor(realE,ELEC_PERMITTIVITY,IMAG);
    }
    else {
      material->SetTensor(realE,ELEC_PERMITTIVITY,REAL);
    }

    if ( !matC ) {
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
        param->Get("geometry")->Get("type", probGeo);
	if ( (probGeo == "axi") || ( probGeo == "plane") ) {
	  if (dirPol == 3) {
	    dirPol = 2;
	  }
	}
	
	if ( hystType == "jiles" ) {
	  material->SetScalar(aJiles,A_JILES,REAL);
	  material->SetScalar(alphaJiles,ALPHA_JILES,REAL);
	  material->SetScalar(kJiles,K_JILES,REAL);
	  material->SetScalar(cJiles,C_JILES,REAL);
	}
	
	material->SetScalar(Psat,P_SATURATION,REAL);
	material->SetScalar(Esat,E_SATURATION,REAL);
	//	material->SetScalar((Integer&)dirPol,P_DIRECTION);
	material->SetScalar(hystType,HYST_MODEL);
	
	
	delete strPtr;    
      }
    }

    Info->PrintMaterial(material);

    if (scaleMatDat) {
      Info->PrintF( "", "\n!!!!!! SCALING with Diag(1e-5, 1e-5, 1e-5, 1e-5,"
                    "1e-5, 1e-5,1e5, 1e5, 1e5) IS ON !!!!!\n\n ");
    }
  }


} // end namespace CoupledField


#undef SSCANF
