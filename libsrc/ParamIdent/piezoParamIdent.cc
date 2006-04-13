#include <iomanip>
#include <PDE/SinglePDE.hh>
#include "Domain/domain.hh"
#include "piezoParamIdent.hh"



namespace CoupledField
{

  // ========================================================================
  // ========================= piezoParamIdent - Part ===========================
  // ========================================================================

  //constructor
  // opens datafiles: measuredData.dat for input, imedCurve.dat and piezoLog.dat for output

  piezoParamIdent :: piezoParamIdent(Domain * adomain,
                                     Integer stepOffset,
                                     Double timeOffset,
                                     std::string driverTag,
                                     Boolean isPartOfSequence)
    :SingleDriver(adomain, stepOffset, timeOffset, 
                  driverTag, isPartOfSequence){

    ENTER_FCN( "piezoParamIdent::piezoParamIdent" );

    Error("Not Working",__FILE__,__LINE__);

  } // end of constructor

  // destructor
  piezoParamIdent :: ~piezoParamIdent()
  {
    ENTER_FCN( "piezoParamIdent::~piezoParamIdent" );


#ifdef USE_LAPACK
    if(lp_af77)
      DeleteArray(lp_af77);
    std::cout<<"delete wf77 in destructor ..."<<std::endl;
    if(&lp_wf77)
      DeleteArray(lp_wf77);
    std::cout<<"delete workf77"<<std::endl;
    if(lp_workf77)
      DeleteArray(lp_workf77);
    std::cout<<"delete rwork77"<<std::endl;
    if(lp_rworkf77)
      DeleteArray(lp_rworkf77);
#endif
  }

  void piezoParamIdent :: SolveProblem() {
    ENTER_FCN( "piezoParamIdent::SolveProblem" );

    
  }// End solveProblem



  // XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx - now some methods are following ...

  void piezoParamIdent::calcAbsImped(Complex & charge, Double & freq, UInt & fstep, 
				     Boolean typeOut){
  } // end norm
  
  

  Double piezoParamIdent::calcEuclidianMatrixNorm(Matrix<Complex> & mat){
    ENTER_FCN("piezoParamIdent::calcEuclidianMatrixNorm");


  } // end calcEuclidianMatrixNorm

  void piezoParamIdent::maxAndEuclNorm(Vector<Complex> & vec, Double & maxNorm, Double & euclNorm){
    ENTER_FCN("piezoParamIdent::maxAndEuclNorm");
    Double maxNormTemp = 0.0;
    maxNorm=0.0;
    euclNorm=0.0;

    for (UInt i=0;i<vec.GetSize();i++){
      maxNormTemp=std::abs(vec[i]);
      euclNorm += maxNormTemp*maxNormTemp;
      if (maxNormTemp>maxNorm)
        maxNorm=maxNormTemp;
    }
    //    euclNorm=std::sqrt(euclNorm);

  } // end maxAndEuclNorm

  void piezoParamIdent::logNorm(Vector<Complex> & vec, Double & logNorm){
    ENTER_FCN("piezoParamIdent::logNorm");
    logNorm=0.0;
    for (UInt i=0;i<vec.GetSize();i++){
      logNorm = logNorm + std::abs(std::log(vec[i]*vec[i]));
    }
    //    euclNorm=std::sqrt(euclNorm);
  } // end logNorm



  void piezoParamIdent::maxAndWeightedResNorm(Vector<Complex> & vec, Double & maxNorm, 
					      Double & wNorm, Vector<Complex> & q_meas){
    ENTER_FCN("piezoParamIdent::maxAndWeightedResNorm");
    Double maxNormTemp = 0.0;
    maxNorm=0.0;
    wNorm=0.0;
    Double Denominator=0.0;

  } // end maxAndWeightedNorm


  void piezoParamIdent::calcNorm2Resid(Vector<Complex> &res, Double & anorm, UInt nrMeasuredData){
    ENTER_FCN("piezoParamIdent::calcNorm2Resid");
    anorm=0.0;
    for (UInt i=0;i<2*nrMeasuredData;i++){
      anorm+=res[i].real()*res[i].real()+ res[i].imag()*res[i].imag();
      anorm=sqrt(anorm);
    }
  } // end calcNorm2Resid

  Double piezoParamIdent::norm2Real(Vector<Complex> &vec){
    ENTER_FCN("piezoParamIdent::realA2norm");
    Double result=0.0; 
    //      Double real_result;
    for(UInt i=0;i<vec.GetSize();i++)
      result+=vec[i].real()*vec[i].real();
    result=sqrt(result);
    return result;
  }

  Double piezoParamIdent::realA2norm(Vector<Complex> &vec){
    ENTER_FCN("piezoParamIdent::realA2norm");
    Double result=0.0; 
    Complex resultC = Complex(0.0,0.0);
    //      Double real_result;
    for(UInt i=0;i<vec.GetSize();i++)
      resultC+=vec[i]*vec[i];
    result=resultC.real();
    //    std::cout <<" \n Result in realA2norm = " << result << std::endl;
    //    result=sqrt(result);
    return result;
  }

  Double piezoParamIdent::a2norm(Vector<Complex> &vec){
    ENTER_FCN("piezoParamIdent::a2norm");
    Double result=0.0; //Complex(0.0,0.0);
    //      Double real_result;
    for(UInt i=0;i<vec.GetSize();i++)
      result+=std::abs(vec[i])*std::abs(vec[i]);
    result=sqrt(result);
    return result;
  }

  Double piezoParamIdent::a2norm(Vector<Double> &vec){
    ENTER_FCN("piezoParamIdent::a2norm");
    Double result=0.0; //Complex(0.0,0.0);
    //      Double real_result;
    for(UInt i=0;i<vec.GetSize();i++)
      result+=std::abs(vec[i])*std::abs(vec[i]);
    result=sqrt(result);
    return result;
  }


  void piezoParamIdent::typeOutSolutionOnScreen(Vector<Complex> & solElecPot,
                                                Vector<Complex> & solMechDispl){
    ENTER_FCN("piezoParamIdent::typeOutSolutionOnScreen");
    Double sol_real, sol_imag;
    //    std::cout<<"\nElecPot: Amplitude & Phase:"<<std::endl;
    for(UInt i=0;i<solElecPot.GetSize();i++){
      //      sol_real=solElecPot[i].real();
      //      sol_imag=solElecPot[i].imag();
      //   std::cout << "solElecPot("<< i<< ")=" << sol_real << " + " << sol_imag <<" i " <<std::endl;
      std::cout<<"ElecPot: Amplitude ("<< i <<") = "<< std::abs(solElecPot[i])
               << ";  Phase ("<< i <<") = "<< std::arg(solElecPot[i])*180/PI<<std::endl;
    }

    for(UInt i=0;i<solMechDispl.GetSize();i++){
      sol_real=solMechDispl[i].real();
      sol_imag=solMechDispl[i].imag();
      std::cout<<"\nMechDispl: Real & Imag :"<<std::endl;
      std::cout << "solMechDispl( " << i<< ")=" << sol_real << " + " << sol_imag <<" i " <<std::endl;
    }
  }// end typeOutSolutionOnSreen

//   void piezoParamIdent::calcInitialResidual(Vector<Complex> & res, Vector<Complex>
  //  & y_hat, Vector<Complex> & PHI_p, UInt fstep, Vector<Complex> & solElecPot, Double & meanValueMechDeformation){
//     StdVector<UInt> bcs_list;
//     ptdomain_->GetGrid()->GetNodesByName(bcs_list,"ep-top"); // for cube3dharmonic; zero, because level=0
//     //bcs_list=ptBCs->GetNodesLevel("pot", 0); // for cubexi, zero, because level=0
//     PHI_p[fstep]=solElecPot[bcs_list[0]];
//     res[fstep]=y_hat[fstep]-PHI_p[fstep];
//     res[y_hat.GetSize()+fstep]=meanValueMechDeformation;
//     std::cout << "residual ( " << fstep << ")=" << res[fstep].real() << " + " << res[fstep].imag()<<" i " <<std::endl;
//   }


  void piezoParamIdent::createMaterialTensorMatrices(Vector<Double> & parameter, 
                                                     Matrix<Double> & couplingMatrix, 
                                                     Matrix<Double> & dielectricMatrix, 
                                                     UInt spaceDim){
    ENTER_FCN("piezoParamIdent::createMaterialTensorMatrices");
    if (spaceDim==2){ // the rotational symmetric case;  couplingMatrix = e
      couplingMatrix[1][0] = couplingMatrix[1][3] = parameter[6]; //e_31
      couplingMatrix[1][1] = parameter[7];      // e_33
      couplingMatrix[0][2] = parameter[5];      //e_15
      dielectricMatrix[0][0] = parameter[8]; // \eps_11
      dielectricMatrix[1][1] = parameter[9]; // \eps_33
    }
    else if (spaceDim==3){
      couplingMatrix[2][0] = couplingMatrix[2][1] = parameter[6]; // e_31
      couplingMatrix[2][2] = parameter[7]; // e_33
      couplingMatrix[1][4] = couplingMatrix[0][5] = parameter[5]; // e_15
      dielectricMatrix[0][0] = dielectricMatrix[1][1] = parameter[8]; // \eps_11
      dielectricMatrix[2][2] = parameter[9]; // \eps_33
    }
  } // end createMaterialTensorMatrix


  void piezoParamIdent::readMeasuredData(Vector<Double> & freqs, 
                                         Vector<Double> & real,
                                         Vector<Double> & imag ,
                                         Vector<Double> & parameter, 
                                         Double & voltage,
                                         UInt & nrMeasuredData, 
                                         Double & thickness, Double & radius,
                                         Double & delta){
    ENTER_FCN( "piezoParamIdent::readMeasuredData" );
  } // end read MeasuredData


  //! Updates material data & updates system matrices!!
  void piezoParamIdent::updateMaterialData(Vector<Double> & parameter, 
					   MaterialMap& ptMaterial) {
    ENTER_FCN("piezoParamIdent::updateMaterialData");    

  } // end updateMaterialData

  void piezoParamIdent::updateComplexMaterialData(Vector<Double> & parameterC, 
						  MaterialMap& ptMaterial) {

    ENTER_FCN("piezoParamIdent::updateComplexMaterialData");    
    //    std::cout<<"updateComplexMaterialData"<<std::endl;


  } // end updateMaterialData

  void piezoParamIdent::setNewParameterSet(Vector<Double> & par,
                                           Vector<Double> &  par_new,
                                           Vector<Double> & scaling,
                                           Double & theta,
					   Vector<Double> & uStep,
                                           Vector<UInt> & whichParameterToUpdate) {

  } // end setNewParameterSet

 
} // end namespace CoupledField

 
