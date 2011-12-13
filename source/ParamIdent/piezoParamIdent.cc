// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <cmath>
#include <iostream>
#include <set>

#include "CoupledPDE/BasePairCoupling.hh"
#include "CoupledPDE/DirectCoupledPDE.hh"
#include "CoupledPDE/PiezoCoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/resultHandler.hh"
#include "Domain/domain.hh"
#include "Domain/entityList.hh"
#include "Domain/grid.hh"
#include "Domain/resultInfo.hh"
#include "General/Enum.hh"
#include "MatVec/exprt/xpr2.hh"
#include "Materials/baseMaterial.hh"
#include "PDE/SinglePDE.hh"
#include "PDE/basePDE.hh"
#include "Utils/tools.hh"
#include "piezoParamIdent.hh"

namespace CoupledField {

// ========================================================================
// ========================= piezoParamIdent - Part ===========================
// ========================================================================



class AdjointParameters;

piezoParamIdent::piezoParamIdent(UInt sequenceStep, bool isPartOfSequence) :
  SingleDriver(sequenceStep, isPartOfSequence ) {

  // Set analysistype
  analysis_ = BasePDE::HARMONIC;
  ptMyPDE_ = NULL;
  piezoCpl_ = NULL;


  residuumParIdent_=1.0e+10;
  resonanceFrequency_=0;
  antiResonanceFrequency_=0;
  globalIterationNr_=0;

  whichNormCriteria_="notSpecified";
  isPartOfSequence_ = isPartOfSequence;
  computeImpedanceCurveAfterStep_=10;
  voltage_=1.0;
  maxNumberInnerLoops_=45;

  driverNode = driverNode->Get("piezoParamIdent");
  driverNode->Get("sequenceStep")->SetValue(sequenceStep);
}

void piezoParamIdent::Init() {

  directCoupling_ = true;
  writeResults_=false;

  std::cout << "++ Opening in and output files" << std::endl;

  std::string filenameMeasuredData="measuredData.dat";
  allMeasuredData_ = new std::ifstream(filenameMeasuredData.c_str(),
      std::basic_ios<char>::in);
  if (!allMeasuredData_) {
    std::cerr << "\n File measuredData.dat does not exist!"<< std::endl;
  }

  std::string filename= "imped.dat";
  impedCurve_ = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);
  if (!impedCurve_) {
    std::cerr << "\n ImpedanceCurve.dat could not be initialized"<< std::endl;
  }

  std::string filenameSynMess= "synMess.dat";
  synMess_ = new std::ofstream(filenameSynMess.c_str(),std::basic_ios<char>::out);
  if (!synMess_) {
    std::cerr << "\n syMess.dat could not be initialized"<< std::endl;
  }

  std::string filenameMechDispl= "mechDispl.dat";
  mechDispl_ = new std::ofstream(filenameMechDispl.c_str(),std::basic_ios<char>::out);
  if (!mechDispl_) {
    std::cerr << "\n mechDispl.dat could not be initialized"<< std::endl;
  }

  std::string filenameParLog= "parLog.dat";
  parLog_ = new std::ofstream(filenameParLog.c_str(),std::basic_ios<char>::out);
  if (!parLog_) {
    std::cerr << "\n parLog.dat could not be initialized"<< std::endl;
  }

  std::string filenameParFinal= "parFinal.dat";
  parFinal_ = new std::ofstream(filenameParFinal.c_str(),std::basic_ios<char>::out);
  if (!parFinal_) {
    std::cerr << "\n parFinal.dat could not be initialized"<< std::endl;
  }

  std::string filenameAllTensors = "allTensors.dat";
  allTensors_ = new std::ofstream(filenameAllTensors.c_str(),std::basic_ios<char>::out);
  if (!allTensors_) {
    std::cerr << "\n allTensors.dat could not be initialized"<< std::endl;
  }

  std::string filenameConfInterval = "confInterval.dat";
  confInterval_ = new std::ofstream(filenameConfInterval.c_str(),std::basic_ios<char>::out);
  if (!confInterval_) {
    std::cerr << "\n confInterval.dat could not be initialized"<< std::endl;
  }

  // get parameter node
    std::string name = "sequenceStep";
    std::string idx = "index";
    std::string one = "1";

    myParam_ = param->GetByVal(name, idx, one )
    ->Get("analysis")->Get("paramIdent");
    // check, if "imagMaterialData" is used
    imagMaterialParam_ = true;
    WARN("at the moment we assume always imaginary material parameters ...");
    PtrParamNode matNode =param->Get("sequenceStep")->Get("couplingList")->Get("direct")
    ->Get("piezoDirect")->Get("materialDataType", ParamNode::PASS);

    if ( matNode ) {
      imagMaterialParam_ =matNode->Get("type")->As<std::string>() == "imagMaterialParameter";
    }

  // query parameters from "paramIdent" node (myParam_)
  myParam_->GetValue("computeCurveAfterIterationStep",
      computeImpedanceCurveAfterStep_);
  myParam_->GetValue("startFreq", startfreq_ );
  startFreq_ = startfreq_;
  myParam_->GetValue("stopFreq", stopfreq_ );
  stopFreq_ = stopfreq_;

  std::string sampling = "linear";
  myParam_->GetValue("sampling",sampling, ParamNode::PASS);
  std::cout << "\n Sampling: " << sampling << std::endl;
  String2Enum( sampling, samplingType_ );

  myParam_->GetValue("numFreq", nrfreq_ );
  // myParam_->Get("numMechMeasurements", numMechMeasurements_ );
  numMechMeasurements_=3;

  // should we calculate the impedance curve?
  myParam_->GetValue("calcImpedanceCurve", CalcImpedanceCurve_, ParamNode::PASS);
  // should we calculate the mechanical displacement curve at selected node?
  myParam_->GetValue("calcMechDisplacementCurve", CalcMechDisplCurve_, ParamNode::PASS);
  // we should always set a maximal number of iterations
  myParam_->GetValue("maxNrIterations", maxNumberNewtonLoops_, ParamNode::PASS );
  // is either the amount of data noise which we assume in our measurements
  // or is the amount of noise we add when we compute synthetically our input data
  myParam_->GetValue("artDataNoise", delta_);

  // discrepancy principle. We stop the iterations when the residual falls
  // below the value of stopRes_. It is intended to equal the quantity tau*delta
  // (see theory of inverse and ill-posed problems)
  if (myParam_->Has("stopRes"))
    myParam_->GetValue("stopRes", stopRes_ );

  // voltage which drives the ceramic
  myParam_->GetValue("excitationVoltage", voltage_);

  PtrParamNode myParamMethod;
  myParamMethod = myParam_->Get("fittingMethod");
  whichMethod_ = myParamMethod->Get("type")->As<std::string>();

  myParamMethod = myParam_->Get("fittingQuantity");
  whichNormCriteria_= myParamMethod->Get("type")->As<std::string>();

  std::cout<< "++ Starting "<< whichMethod_ << " method ..."<< std::endl;

} // end of constructor

// destructor
piezoParamIdent::~piezoParamIdent() {

  if (allMeasuredData_)
    allMeasuredData_->close();
  if (impedCurve_)
    impedCurve_->close();
    if (parLog_)
    parLog_->close();
  if (synMess_)
    synMess_->close();
    if (mechDispl_)
    mechDispl_->close();
  if (confInterval_)
    confInterval_->close();
    if (inMess_)
    inMess_.close();
  if (inMechMess_)
    inMechMess_.close();
}

void piezoParamIdent::SolveProblem(bool write_results, PtrParamNode given_analysis_id, AdjointParameters* adjointParams) {

  ResultHandler * resHandler = domain->GetResultHandler();
  InitializePDEs();

  DirectCoupledPDE* ptCoupledPDE = domain->GetDirectCoupledPDE();

  // Section added by ahauck July 4th
  // ================================

  // obtain pointer to piezo coupling object
  StdVector<BasePairCoupling*>* cpl = ptCoupledPDE->GetCouplingsObject();
  for( UInt i = 0; i < cpl->GetSize(); i++ ) {
    if( (*cpl)[i]->GetName() == "piezoDirect" ) {
      piezoCpl_ = dynamic_cast<PiezoCoupling*>((*cpl)[i]);
      break;
    }
  }
  // Check for calculation of charges
  PtrParamNode chargeNode;
  if( myParam_->Has( "calcCharge" ) ) {
    chargeNode = myParam_->Get( "calcCharge" );

    // check, if piezo coupling element is present
    if( piezoCpl_ == NULL )EXCEPTION( "No piezo coupling object found");

    // read surface elements for charge computation and neighbor region
    std::string surfName = chargeNode->Get("name")->As<std::string>();
    std::string neighborName = chargeNode->Get("neighborRegion")->As<std::string>();

    // obtain entitylist from grid
    shared_ptr<EntityList> chargeSurf = domain->GetGrid()
    ->GetEntityList( EntityList::SURF_ELEM_LIST, surfName,
                     EntityList::REGION );

    chargeNeighborRegion_ =
      domain->GetGrid()->GetRegion().Parse(neighborName);

    // obtain result info from piezo-coupling object

    shared_ptr<ResultInfo> chargeInfo;
    BasePairCoupling::ResultSet infos = piezoCpl_->GetAvailResults();
    BasePairCoupling::ResultSet::const_iterator it;
    for( it = infos.begin(); it != infos.end(); it++ ) {
      if( (*it)->resultType == ELEC_CHARGE ) {
        chargeInfo = *it;
        break;
      }
    }

    // create Result<Complex> object and store it for later use
    charges_ =
      shared_ptr<Result<Complex> >( new Result<Complex>() );
    charges_->SetResultInfo( chargeInfo );
    charges_->SetEntityList( chargeSurf );
  }

  // ===== end of section for initializing charge computation ====

  // get PDEs
  ptPDE1_=domain-> GetSinglePDE("mechanic");
  ptPDE2_=domain-> GetSinglePDE("electrostatic");

  subdomsMech_ = ptPDE1_->getPDE_subdoms();
  subdomsElec_ = ptPDE1_->getPDE_subdoms();

  // number of parameters involved.
  // pure piezo case, 1 region
  if (subdomsMech_.GetSize()==1) {
    nrParameter_ = 10;
    actNrParameter_ = 10;
    actNrParameterC_ = 10;
  }
  // piezo case with one adaption layer
  // we have additionally the Lame parameters of a
  // second material (just mechanics)
    else if (subdomsMech_.GetSize()==2) {
    nrParameter_ = 12;
    actNrParameter_ = 12;
    actNrParameterC_ = 12;
  }
    // piezo case with two adaption layers
    // we have additionally the Lame parameters of a
    // third material (just mechanics)
  else if (subdomsMech_.GetSize()==3) {
     nrParameter_ = 14;
     actNrParameter_ = 14;
     actNrParameterC_ = 14;
   }

  parameter_.Resize(nrParameter_);
  parameterC_.Resize(nrParameter_);
  parameter_.Init();
  parameterC_.Init();

  whichParameterToUpdate_.Resize(nrParameter_);
  whichParameterToUpdateC_.Resize(nrParameter_);
  whichParameterToUpdate_.Init();
  whichParameterToUpdateC_.Init();

  UInt highestAssumableNrOfMeasData=100;
  freqs_.Resize(highestAssumableNrOfMeasData);
  freqs_.Init();

  //will be overwritten in readMeasuredData ...

  // the following passage reads Data from file measuredData.dat
  // The rows are containing the values of the given frequencies, such as phase and amplitude!
  readMeasuredData();

  Vector<Double> freqsTemp;
  freqsTemp = freqs_;
  freqs_.Resize(nrMeasuredData_);
  freqs_.Init();

  for (UInt i=0; i<nrMeasuredData_; i++)
    freqs_[i]=freqsTemp[i];

  std::cout<<"\nSampling points:"<<std::endl;
  std::cout<<freqs_<<std::endl;

  y_hat_.Resize(nrMeasuredData_);

  s_.Resize(actNrParameter_+actNrParameterC_);
  scaling_.Resize(nrParameter_);
  scalingC_.Resize(nrParameter_);
  F_hat_.Resize(nrMeasuredData_);
  parameter_new_.Resize(nrParameter_);

  actNrParameter_=0;
  actNrParameterC_=0;

  y_hat_.Init();
  F_hat_.Init();
  s_.Init();
  scaling_.Init();
  scalingC_.Init();
  parameter_new_.Init();

  ptPDE_->WriteGeneralPDEdefines();

  // how many parameters are there actually to set?
  for (UInt i=0; i<whichParameterToUpdate_.GetSize(); i++)
    if (whichParameterToUpdate_[i]==1)
      actNrParameter_++;

  for (UInt i=0; i<whichParameterToUpdateC_.GetSize(); i++)
    if (whichParameterToUpdateC_[i]==1)
      actNrParameterC_++;

  if (actNrParameter_!=0) {
    whichParToUpInd_.Resize(actNrParameter_);
    whichParToUpInd_.Init();
  }

  if (actNrParameterC_!=0) {
    whichParToUpIndC_.Resize(actNrParameterC_);
    whichParToUpIndC_.Init();
  }

  UInt intTemp=0;

  for (UInt i=0; i<whichParameterToUpdate_.GetSize(); i++)
    if (whichParameterToUpdate_[i]==1) {
      whichParToUpInd_[intTemp]=i;
      intTemp++;
    }

  intTemp=0;
  for (UInt i=0; i<whichParameterToUpdateC_.GetSize(); i++)
    if (whichParameterToUpdateC_[i]==1) {
      whichParToUpIndC_[intTemp]=i;
      intTemp++;
    }


  StdVector<BasePairCoupling*>* ptCoupling = ptCoupledPDE->GetCouplingsObject();

  ptMaterialMech_ = ptPDE1_[0].getPDEMaterialData(); // Pointer to mech. MaterialData
  ptMaterialElec_ = ptPDE2_[0].getPDEMaterialData(); // Pointer to elec. MaterialData
  ptMaterialPiezo_ = (*ptCoupling)[0]->getPDEMaterialData(); // Pointer to piezo MaterialData

  if (imagMaterialParam_ ) {
    updateComplexMaterialData(parameterC_);
  }

  updateMaterialData(parameter_);

  parameterIncrement_=parameter_;

  // <<<<<<<<<<<<<< calc impedance curve <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  if (CalcImpedanceCurve_ == true|| CalcMechDisplCurve_ == true) {
    Vector<Double> freqsTemp = freqs_;
    freqs_.Resize(nrfreq_);
    Double startFreqTemp;
    startFreqTemp=startfreq_;
    Double freqincr=(stopfreq_-startfreq_)/nrfreq_;
    for (UInt i=0; i<nrfreq_; i++) {
      startFreqTemp+=freqincr;
      freqs_[i]=startFreqTemp;
    }
    calcImpedanceCurve();
    freqs_ = freqsTemp;

    if (CalcMechDisplCurve_ == true)
      std::cout
          <<"->Mechanical displacement curve (before fitting) is written into file 'mechDispl.dat'. "
          <<std::endl;
    if (CalcImpedanceCurve_ == true)
      std::cout
          <<"->Impedance curve (before fitting) is written into file 'imped.dat'. "
          <<std::endl;
      }

  // xxxxxxxxxxxxxxxxxxxxxxx Choose different regularizing solvers here xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx //

  computeScaling();

  if (whichMethod_ == "Landweber"|| whichMethod_ == "steepestDescent"
      || whichMethod_ == "minimalError") {

    if (whichMethod_ == "Landweber")
      std::cout<<"++ Nonlinear Landweber's iteration ... "<<std::endl;
    else
      std::cout<<"++ Modified Landweber's iteration with "<< whichMethod_
          << " relaxation parameter choice "<<std::endl;

    Vector<Double> newFreqs;
    readInMeasurement(newFreqs);
    freqs_ = newFreqs;
    nonlinLandweber();
  }

  else if (whichMethod_ == "Newton") {
    if (imagMaterialParam_ ) {
      std::cerr
          <<" This version of nuMethods only works for real valued paraeters "
          <<std::endl;
      std::cerr<<" Choose either nuMethods=8 in your input file or "<<std::endl;
      std::cerr<<" use only real - valued parameters "<<std::endl;

    }
    std::cout<<"++ Newton - nu Methods "<<std::endl;
    UInt nrNuMethods=0;
    newtonCounter_=0;
    nu_=1.5;

    Vector<Double> newFreqs;
    Vector<Double> newFreqsMech;
    readInMeasurement(newFreqs);
    evaluateMeasuredData(freqs_, real_, imag_, y_hat_);

    while (nrNuMethods<maxNumberNewtonLoops_) {

      if (residuumParIdent_ <= stopRes_) {
        std::cout
            <<"\nTermination of Newton iteration, since the norm of the residual "
            << residuumParIdent_ <<" is smaller than "<< stopRes_ <<std::endl;
        break;
      }

      nuMethods();
      nrNuMethods++;

      if (nrNuMethods%computeImpedanceCurveAfterStep_==0) {
        UInt nrfreqTemp=nrfreq_;
        Vector<Double> freqsTemp = freqs_;
        freqs_.Resize(nrfreqTemp);
        Double startFreqTemp;
        startFreqTemp=startfreq_;
        Double freqincr=(stopfreq_-startfreq_)/nrfreqTemp;
        for (UInt i=0; i<nrfreqTemp; i++) {
          startFreqTemp+=freqincr;
          freqs_[i]=startFreqTemp;
        }
        if (impedCurve_)
          impedCurve_->close();
        if (mechDispl_)
          mechDispl_->close();

        std::string filename= "imped";
        filename.append(".dat");
        impedCurve_ = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);
        filename="mechDispl.dat";
        mechDispl_ = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);

        calcImpedanceCurve();
        freqs_ = freqsTemp;
        nrfreq_=nrfreqTemp;

      }

      for (UInt i=0; i<parameter_.GetSize(); i++)
        *parFinal_<<parameter_[i]<<", ";
      *parFinal_<<"/"<<std::endl;
      newtonCounter_++;
    }

  }

  else if (whichMethod_=="NewtonComplex") {
    std::cout<<"++ Newton - nu MethodsC "<<std::endl;

    Vector<Double> newFreqs;
    Vector<Double> newFreqsMech;
    readInMeasurement(newFreqs);
    evaluateMeasuredData(freqs_, real_, imag_, y_hat_); // out of new measurements

    nu_=1.2;

    UInt nrNuMethodsC=0;
    newtonCounter_=0;
    while (nrNuMethodsC<maxNumberNewtonLoops_) {

      if (residuumParIdent_ <= stopRes_) {
        std::cout
            <<"Termination of Newton iteration, since the norm of the residual "
            << residuumParIdent_ <<" is smaller than "<< stopRes_ <<std::endl;
        break;
      }

      nuMethodsC2();

      if ((nrNuMethodsC+1)%computeImpedanceCurveAfterStep_==0) {
        UInt nrfreqTemp=nrfreq_;
        Vector<Double> freqsTemp = freqs_;
        freqs_.Resize(nrfreqTemp);
        Double startFreqTemp;
        startFreqTemp=startfreq_;
        Double freqincr=(stopfreq_-startfreq_)/nrfreqTemp;
        for (UInt i=0; i<nrfreqTemp; i++) {
          startFreqTemp+=freqincr;
          freqs_[i]=startFreqTemp;
        }
        if (impedCurve_)
          impedCurve_->close();
        if (mechDispl_)
          mechDispl_->close();

        std::string filename= "imped";
        filename.append(".dat");
        impedCurve_ = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);
        filename="mechDispl.dat";
        mechDispl_ = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);

        calcImpedanceCurve();
        freqs_ = freqsTemp;
        nrfreq_=nrfreqTemp;
      }

      for (UInt i=0; i<parameter_.GetSize(); i++)
        *parFinal_<<parameter_[i]<<", ";
      *parFinal_<<"/"<<std::endl;
      for (UInt i=0; i<parameterC_.GetSize(); i++)
        *parFinal_<<parameterC_[i]<<", ";
      *parFinal_<<"/"<<std::endl;

      nrNuMethodsC++;
      newtonCounter_++;

    }
  }

  else if (whichMethod_=="confidenceIntervals") {

    std::cout<<"++ Computing bounds of confidence intervals ..."<<std::endl;
    std::cout.precision(8);

    Vector<Complex> jacobi;
    Vector<Complex> jacobiH;

    Matrix<Complex>covTempWithOutWeight(actNrParameter_+actNrParameterC_,
        actNrParameter_+actNrParameterC_);
    Matrix<Complex>covWithOutWeight(actNrParameter_+actNrParameterC_,
        actNrParameter_+actNrParameterC_);

    covWithOutWeight.Init();
    covTempWithOutWeight.Init();

    for (UInt actFreq=0; actFreq<freqs_.GetSize(); actFreq++) {

      // std::cout<<freqs_<<std::endl;
      computeVecOfSingleDeriv(jacobi, freqs_[actFreq]);
      jacobiH.Resize(actNrParameter_+actNrParameterC_);
      jacobiH.Init();

      for (UInt i=0; i<jacobi.GetSize(); i++)
        jacobiH[i]=Complex(jacobi[i].real(), -jacobi[i].imag());

      globalJacobi_=jacobi;
      globalJacobiH_=jacobiH;

      for (UInt i=0; i<actNrParameter_+actNrParameterC_; i++)
        for (UInt j=0; j<actNrParameter_+actNrParameterC_; j++) {

          covTempWithOutWeight[i][j]= globalJacobiH_[i]*globalJacobi_[j]; ///
        }

      covWithOutWeight+=covTempWithOutWeight;

    } // end for actFreq ..


    globalCov_=covWithOutWeight;

    invert(covWithOutWeight);

    globalCov_=covWithOutWeight;

    Double chi=0.0;

    // Values for the chi-squared probability distribution
    if (actNrParameter_+actNrParameterC_==1)
      chi=0.00016;
    else if (actNrParameter_+actNrParameterC_==2)
      chi=0.02;
    else if (actNrParameter_+actNrParameterC_==3)
      chi=0.115;
    else if (actNrParameter_+actNrParameterC_==4)
      chi=0.3;
    else if (actNrParameter_+actNrParameterC_==5)
      chi=0.55;
    else if (actNrParameter_+actNrParameterC_==6)
      chi=0.87;
    else if (actNrParameter_+actNrParameterC_==7)
      chi=1.24;
    else if (actNrParameter_+actNrParameterC_==8)
      chi=1.65;
    else if (actNrParameter_+actNrParameterC_==9)
      chi=2.09;
    else if (actNrParameter_+actNrParameterC_==10)
      chi=2.5;
    else if (actNrParameter_+actNrParameterC_==11)
      chi=3.1;
    else if (actNrParameter_+actNrParameterC_==12)
      chi=3.6;
    else if (actNrParameter_+actNrParameterC_==13)
      chi=4.1;
    else if (actNrParameter_+actNrParameterC_==14)
      chi=4.7;
    else if (actNrParameter_+actNrParameterC_==15)
      chi=5.2;
    else if (actNrParameter_+actNrParameterC_==16)
      chi=5.8;
    else if (actNrParameter_+actNrParameterC_==17)
      chi=6.4;
    else if (actNrParameter_+actNrParameterC_==18)
      chi=7.0;
    else if (actNrParameter_+actNrParameterC_==19)
      chi=7.6;
    else if (actNrParameter_+actNrParameterC_==20)
      chi=8.3;
    else {
      std::cout<<"There are too many parameters in the game ... program has no information about chi-squared probability distribution" << std::endl;
      std::cout<<"Try again with a smaller amount of parameters" <<std::endl;
      getchar();
    }

    std::cout<<"\nUpper bounds of confidence intervals computed with"<<std::endl;
    std::cout<<"the 0.99 quantile of the chi-squared probability distribution:\n"<<std::endl;

    for (UInt ii=0; ii<actNrParameter_+actNrParameterC_; ii++){
      std::cout<<"|p("<<ii<<")_exact - p("<<ii<<")_computed| < " << sqrt(globalCov_[ii][ii].real()*chi*0.9801)<<std::endl;
      *confInterval_<<sqrt(globalCov_[ii][ii].real()*chi*0.9801)<<std::endl;
    }

    std::cout.precision(6);

  }


  else if (whichMethod_=="leastSquares") {

    std::cout<<"++ Least squares fitting started with maximal "
        << maxNumberNewtonLoops_ << " descent steps ... "<<std::endl;
    Vector<Double> newFreqs;
    Vector<Double> newFreqsMech;
    readInMeasurement(newFreqs);
    evaluateMeasuredData(freqs_, real_, imag_, y_hat_); // out of new measurements
    std::cout<<"++ Fitting will be performed with the following "
        << nrMeasuredData_ << " frequencies:"<<std::endl;
    leastSquare();
  }

  else
    std::cout
        <<"\n There was no valid fitting method specified - see your xml -file "
        <<std::endl;

  // xxxxxxxxxxxxxxxxxxxxxxx End of choice xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx //



  if (whichMethod_!="confidenceIntervals"){

    if (maxNumberNewtonLoops_!=0) {
      std::cout<<"\n\n *** FINALLY CALCULATED PARAMETERS *** \n"<<std::endl;

      *parFinal_<<"r) "<<std::endl;
      for (UInt i=0; i<parameter_.GetSize(); i++) {
        std::cout<<"par["<< i<<"]="<< parameter_[i]<<" + "<< parameterC_[i]<<"i"
          <<std::endl;
        *parFinal_<<parameter_[i]<<", ";
      }
      *parFinal_<<"/"<<std::endl;
      std::cout<<"\n===================================== \n "<<std::endl;
    }
  }

  writeTensorsInFile();

  // <<<<<<<<<<<<<< for a hopefully nicely fitted curve after identification !! <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  if ((CalcImpedanceCurve_ == true || CalcMechDisplCurve_ == true)
      && maxNumberNewtonLoops_!=0) {
    writeResults_=true;
    Vector<Double> freqsTemp(freqs_);
    freqs_.Clear(); // this is necessary to prevent a memory fault in debug mode! ??!
    freqs_.Resize(nrfreq_);
    Double freqincr=(stopfreq_-startfreq_)/nrfreq_;
    for (UInt i=0; i<nrfreq_; i++) {
      startfreq_+=freqincr;
      freqs_[i]=startfreq_;
    }

    if (impedCurve_)
      impedCurve_->close();
    if (mechDispl_)
      mechDispl_->close();

    std::string filename= "imped";
    filename.append(".dat");
    impedCurve_ = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);
    filename="mechDispl.dat";
    mechDispl_ = new std::ofstream(filename.c_str(),std::basic_ios<char>::out);

    calcImpedanceCurve();

    writeResults_=false;
    freqs_ = freqsTemp;
  }

  std::cout<<"\n-> New parameters are written into file 'parFinal.dat'. "
      <<std::endl;
  std::cout<<"-> New tensors are written into file 'allTensors.dat'. "
      <<std::endl;
  std::cout
      <<"-> Norm of residual and all iterates are written into file 'parLog.dat'. "
      <<std::endl;

  if (CalcMechDisplCurve_ == true)
    std::cout
        <<"-> Mechanical displacement curve (after fitting) is written into file 'mechDispl.dat'. "
        <<std::endl;
  if (CalcImpedanceCurve_ == true)
    std::cout
        <<"-> Impedance curve (after fitting) is written into file 'imped.dat'. "
        <<std::endl;
  if (whichMethod_=="confidenceIntervals")
    std::cout<<"-> Upper bounds of confidence interval sizes are written into file 'confInterval.dat'. " <<std::endl;

  // notify resultHandler about finishing of current sequence step
  if (CalcImpedanceCurve_ == true|| CalcMechDisplCurve_ == true) {
    resHandler->FinishMultiSequenceStep();
    resHandler->Finalize();
  }

}// End solveProblem


// XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXxx

//! Updates material data & updates system matrices!!
void piezoParamIdent::updateMaterialData(Vector<Double> & parameter_) {
  //    std::cout<<"++ updateMaterialData " <<std::endl;

  UInt dim=ptPDE1_->getPDE_spaceDim();

  //std::cout<<"spacedimesion "<< dim <<std::endl;
  Matrix<Double> stiffTensor;
  Matrix<Double> piezoTensor;
  Matrix<Double> permTensor;

  stiffTensor.Resize(6, 6);
  piezoTensor.Resize(3, 6);
  permTensor.Resize(3, 3);

  stiffTensor.Init();
  piezoTensor.Init();
  permTensor.Init();

  subdomsMech_ = ptPDE1_->getPDE_subdoms();

  if (parameter_.GetSize()==10)
      if (subdomsMech_.GetSize()!=1)
        std::cout<<"Number of parameters is not consistent with number of regions" <<std::endl;

  if (parameter_.GetSize()==12)
      if (subdomsMech_.GetSize()!=2)
        std::cout<<"Number of parameters is not consistent with number of regions" <<std::endl;

  if (parameter_.GetSize()==14)
    if (subdomsMech_.GetSize()!=3)
      std::cout<<"Number of parameters is not consistent with number of regions" <<std::endl;

  // Materialparameteridentifizierung fuer nur eine piezoelektrische keramik:
  if (subdomsMech_.GetSize()==1) {
    RegionIdType actId = subdomsMech_[0];

    stiffTensor[0][0] = parameter_[0]; //c_11
    stiffTensor[1][1] = parameter_[0]; //c_11
    stiffTensor[2][2] = parameter_[1]; //c_33
    stiffTensor[0][1] = parameter_[2]; //c_12
    stiffTensor[1][0] = parameter_[2]; //c_12
    stiffTensor[0][2] = parameter_[3]; //c_13
    stiffTensor[2][0] = parameter_[3]; //c_13
    stiffTensor[1][2] = parameter_[3]; //c_13
    stiffTensor[2][1] = parameter_[3]; //c_13
    stiffTensor[3][3] = parameter_[4]; //c_44
    stiffTensor[4][4] = parameter_[4]; //c_44
    stiffTensor[5][5] = 0.5*(parameter_[0]-parameter_[2]); //c_66

    if (dim==2) { // also axi ...
      stiffTensor[0][0] = parameter_[0]; //c_11
      stiffTensor[2][2] = parameter_[0]; //c_11
      stiffTensor[1][1] = parameter_[1]; //c_33
      stiffTensor[0][2] = parameter_[2]; //c_12
      stiffTensor[2][0] = parameter_[2]; //c_12
      stiffTensor[0][1] = parameter_[3]; //c_13
      stiffTensor[1][0] = parameter_[3]; //c_13
      stiffTensor[1][2] = parameter_[3]; //c_13
      stiffTensor[2][1] = parameter_[3]; //c_13
      stiffTensor[3][3] = parameter_[4]; //c_44
      stiffTensor[5][5] = parameter_[4]; //c_44
      stiffTensor[4][4] = 0.5*(parameter_[0]-parameter_[2]); //c_66
    }

    piezoTensor[1][3]=parameter_[5]; //e_15
    piezoTensor[0][4]=parameter_[5]; //e_15
    piezoTensor[2][0]=parameter_[6]; //e_31
    piezoTensor[2][1]=parameter_[6]; //e_31
    piezoTensor[2][2]=parameter_[7]; // e_33

    if (dim==2) {
      piezoTensor[2][3]=parameter_[5]; //e_15
      piezoTensor[0][5]=parameter_[5]; //e_15
      piezoTensor[1][0]=parameter_[6]; //e_31
      piezoTensor[1][2]=parameter_[6]; //e_31
      piezoTensor[1][1]=parameter_[7]; // e_33
    }

    permTensor[0][0] = parameter_[8]; //eps_11
    permTensor[1][1] = parameter_[8]; //eps_11
    permTensor[2][2] = parameter_[9]; //eps_33

    if (dim==2) {
      permTensor[0][0] = parameter_[8]; //eps_11
      permTensor[1][1] = parameter_[9]; //eps_33
    }

    ptMaterialPiezo_[actId]->SetTensor(piezoTensor, PIEZO_TENSOR, Global::REAL);
    ptMaterialMech_[actId]->SetTensor(stiffTensor, MECH_STIFFNESS_TENSOR, Global::REAL);
    ptMaterialElec_[actId]->SetTensor(permTensor, ELEC_PERMITTIVITY, Global::REAL);
  }

  // Materialparameteridentifizierung fuer eine piezoelektrische Keramik
  // und eine Anpassschicht
  if (subdomsMech_.GetSize()==2) {
    RegionIdType actId = subdomsMech_[0];
    RegionIdType actId1 = subdomsMech_[1];

    stiffTensor[0][0] = parameter_[0]; //c_11
    stiffTensor[1][1] = parameter_[0]; //c_11
    stiffTensor[2][2] = parameter_[1]; //c_33
    stiffTensor[0][1] = parameter_[2]; //c_12
    stiffTensor[1][0] = parameter_[2]; //c_12
    stiffTensor[0][2] = parameter_[3]; //c_13
    stiffTensor[2][0] = parameter_[3]; //c_13
    stiffTensor[1][2] = parameter_[3]; //c_13
    stiffTensor[2][1] = parameter_[3]; //c_13
    stiffTensor[3][3] = parameter_[4]; //c_44
    stiffTensor[4][4] = parameter_[4]; //c_44
    stiffTensor[5][5] = 0.5*(parameter_[0]-parameter_[2]); //c_66


    // Set Lame parameters for adaption layer
    Double lambda2nu=parameter_[11]+2*parameter_[10];
    Matrix<Double> stiffTensorAdaption(6,6);

    stiffTensorAdaption[0][0] = lambda2nu; //c_11
    stiffTensorAdaption[1][1] = lambda2nu; //c_11
    stiffTensorAdaption[2][2] = lambda2nu; //c_33
    stiffTensorAdaption[0][1] = parameter_[11]; //c_12
    stiffTensorAdaption[1][0] = parameter_[11]; //c_12
    stiffTensorAdaption[0][2] = parameter_[11]; //c_13
    stiffTensorAdaption[2][0] = parameter_[11]; //c_13
    stiffTensorAdaption[1][2] = parameter_[11]; //c_13
    stiffTensorAdaption[2][1] = parameter_[11]; //c_13
    stiffTensorAdaption[3][3] = parameter_[10]; //c_44
    stiffTensorAdaption[4][4] = parameter_[10]; //c_44
    stiffTensorAdaption[5][5] = parameter_[10]; //c_66


    if (dim==2) { // also axi ...
      stiffTensor[0][0] = parameter_[0]; //c_11
      stiffTensor[2][2] = parameter_[0]; //c_11
      stiffTensor[1][1] = parameter_[1]; //c_33
      stiffTensor[0][2] = parameter_[2]; //c_12
      stiffTensor[2][0] = parameter_[2]; //c_12
      stiffTensor[0][1] = parameter_[3]; //c_13
      stiffTensor[1][0] = parameter_[3]; //c_13
      stiffTensor[1][2] = parameter_[3]; //c_13
      stiffTensor[2][1] = parameter_[3]; //c_13
      stiffTensor[3][3] = parameter_[4]; //c_44
      stiffTensor[5][5] = parameter_[4]; //c_44
      stiffTensor[4][4] = 0.5*(parameter_[0]-parameter_[2]); //c_66

      // Set Lame parameters for adaption layer
      Double lambda2nu=parameter_[11]+2*parameter_[10];
      Matrix<Double> stiffTensorAdaption(6,6);

      stiffTensorAdaption[0][0] = lambda2nu; //c_11
      stiffTensorAdaption[2][2] = lambda2nu; //c_11
      stiffTensorAdaption[1][1] = lambda2nu; //c_33
      stiffTensorAdaption[0][2] = parameter_[11]; //c_12
      stiffTensorAdaption[2][0] = parameter_[11]; //c_12
      stiffTensorAdaption[0][1] = parameter_[11]; //c_13
      stiffTensorAdaption[1][0] = parameter_[11]; //c_13
      stiffTensorAdaption[1][2] = parameter_[11]; //c_13
      stiffTensorAdaption[2][1] = parameter_[11]; //c_13
      stiffTensorAdaption[3][3] = parameter_[10]; //c_44
      stiffTensorAdaption[5][5] = parameter_[10]; //c_44
      stiffTensorAdaption[4][4] = parameter_[10]; //c_66

    }

    piezoTensor[1][3]=parameter_[5]; //e_15
    piezoTensor[0][4]=parameter_[5]; //e_15
    piezoTensor[2][0]=parameter_[6]; //e_31
    piezoTensor[2][1]=parameter_[6]; //e_31
    piezoTensor[2][2]=parameter_[7]; // e_33

    if (dim==2) {
      piezoTensor[2][3]=parameter_[5]; //e_15
      piezoTensor[0][5]=parameter_[5]; //e_15
      piezoTensor[1][0]=parameter_[6]; //e_31
      piezoTensor[1][2]=parameter_[6]; //e_31
      piezoTensor[1][1]=parameter_[7]; // e_33
    }

    permTensor[0][0] = parameter_[8]; //eps_11
    permTensor[1][1] = parameter_[8]; //eps_11
    permTensor[2][2] = parameter_[9]; //eps_33

    if (dim==2) {
      permTensor[0][0] = parameter_[8]; //eps_11
      permTensor[1][1] = parameter_[9]; //eps_33
    }

    ptMaterialPiezo_[actId]->SetTensor(piezoTensor, PIEZO_TENSOR, Global::REAL);
    ptMaterialMech_[actId]->SetTensor(stiffTensor, MECH_STIFFNESS_TENSOR, Global::REAL);
    ptMaterialElec_[actId]->SetTensor(permTensor, ELEC_PERMITTIVITY, Global::REAL);

    ptMaterialMech_[actId1]->SetTensor(stiffTensorAdaption, MECH_STIFFNESS_TENSOR, Global::REAL);
  }

  // Materialparameteridentifizierung fuer eine piezoelektrische Keramik,
  //  Anpassschicht und backing
  if (subdomsMech_.GetSize()==3) {
    RegionIdType actId = subdomsMech_[0];
    RegionIdType actId1 = subdomsMech_[1];
    RegionIdType actId2 = subdomsMech_[2];

    stiffTensor[0][0] = parameter_[0]; //c_11
    stiffTensor[1][1] = parameter_[0]; //c_11
    stiffTensor[2][2] = parameter_[1]; //c_33
    stiffTensor[0][1] = parameter_[2]; //c_12
    stiffTensor[1][0] = parameter_[2]; //c_12
    stiffTensor[0][2] = parameter_[3]; //c_13
    stiffTensor[2][0] = parameter_[3]; //c_13
    stiffTensor[1][2] = parameter_[3]; //c_13
    stiffTensor[2][1] = parameter_[3]; //c_13
    stiffTensor[3][3] = parameter_[4]; //c_44
    stiffTensor[4][4] = parameter_[4]; //c_44
    stiffTensor[5][5] = 0.5*(parameter_[0]-parameter_[2]); //c_66


    // Set Lame parameters for adaption layer
    Double lambda2nu=parameter_[11]+2*parameter_[10];
    Matrix<Double> stiffTensorAdaption(6,6);
    Matrix<Double> stiffTensorBacking(6,6);

    stiffTensorAdaption[0][0] = lambda2nu;
    stiffTensorAdaption[1][1] = lambda2nu;
    stiffTensorAdaption[2][2] = lambda2nu;
    stiffTensorAdaption[0][1] = parameter_[11];
    stiffTensorAdaption[1][0] = parameter_[11];
    stiffTensorAdaption[0][2] = parameter_[11];
    stiffTensorAdaption[2][0] = parameter_[11];
    stiffTensorAdaption[1][2] = parameter_[11];
    stiffTensorAdaption[2][1] = parameter_[11];
    stiffTensorAdaption[3][3] = parameter_[10];
    stiffTensorAdaption[4][4] = parameter_[10];
    stiffTensorAdaption[5][5] = parameter_[10];

    lambda2nu=parameter_[13]+2*parameter_[12];

    stiffTensorBacking[0][0] = lambda2nu;
    stiffTensorBacking[1][1] = lambda2nu;
    stiffTensorBacking[2][2] = lambda2nu;
    stiffTensorBacking[0][1] = parameter_[13];
    stiffTensorBacking[1][0] = parameter_[13];
    stiffTensorBacking[0][2] = parameter_[13];
    stiffTensorBacking[2][0] = parameter_[13];
    stiffTensorBacking[1][2] = parameter_[13];
    stiffTensorBacking[2][1] = parameter_[13];
    stiffTensorBacking[3][3] = parameter_[12];
    stiffTensorBacking[4][4] = parameter_[12];
    stiffTensorBacking[5][5] = parameter_[12];


    if (dim==2) { // also axi ...
      stiffTensor[0][0] = parameter_[0]; //c_11
      stiffTensor[2][2] = parameter_[0]; //c_11
      stiffTensor[1][1] = parameter_[1]; //c_33
      stiffTensor[0][2] = parameter_[2]; //c_12
      stiffTensor[2][0] = parameter_[2]; //c_12
      stiffTensor[0][1] = parameter_[3]; //c_13
      stiffTensor[1][0] = parameter_[3]; //c_13
      stiffTensor[1][2] = parameter_[3]; //c_13
      stiffTensor[2][1] = parameter_[3]; //c_13
      stiffTensor[3][3] = parameter_[4]; //c_44
      stiffTensor[5][5] = parameter_[4]; //c_44
      stiffTensor[4][4] = 0.5*(parameter_[0]-parameter_[2]); //c_66
    }

    piezoTensor[1][3]=parameter_[5]; //e_15
    piezoTensor[0][4]=parameter_[5]; //e_15
    piezoTensor[2][0]=parameter_[6]; //e_31
    piezoTensor[2][1]=parameter_[6]; //e_31
    piezoTensor[2][2]=parameter_[7]; // e_33

    if (dim==2) {
      piezoTensor[2][3]=parameter_[5]; //e_15
      piezoTensor[0][5]=parameter_[5]; //e_15
      piezoTensor[1][0]=parameter_[6]; //e_31
      piezoTensor[1][2]=parameter_[6]; //e_31
      piezoTensor[1][1]=parameter_[7]; // e_33
    }

    permTensor[0][0] = parameter_[8]; //eps_11
    permTensor[1][1] = parameter_[8]; //eps_11
    permTensor[2][2] = parameter_[9]; //eps_33

    if (dim==2) {
      permTensor[0][0] = parameter_[8]; //eps_11
      permTensor[1][1] = parameter_[9]; //eps_33
    }


    ptMaterialPiezo_[actId]->SetTensor(piezoTensor, PIEZO_TENSOR, Global::REAL);
    ptMaterialMech_[actId]->SetTensor(stiffTensor, MECH_STIFFNESS_TENSOR, Global::REAL);
    ptMaterialElec_[actId]->SetTensor(permTensor, ELEC_PERMITTIVITY, Global::REAL);

    ptMaterialMech_[actId1]->SetTensor(stiffTensorAdaption, MECH_STIFFNESS_TENSOR, Global::REAL);
    ptMaterialMech_[actId2]->SetTensor(stiffTensorBacking, MECH_STIFFNESS_TENSOR, Global::REAL);
  }

} // end updateMaterialData

void piezoParamIdent::updateComplexMaterialData(Vector<Double> & parameterC_) {
  //    std::cout<<"++ updateComplexMaterialData " <<std::endl;

  Matrix<Double> stiffTensorC;
  Matrix<Double> piezoTensorC;
  Matrix<Double> permTensorC;

  stiffTensorC.Resize(6, 6);
  piezoTensorC.Resize(3, 6);
  permTensorC.Resize(3, 3);

  stiffTensorC.Init();
  piezoTensorC.Init();
  permTensorC.Init();

  UInt dim=ptPDE1_->getPDE_spaceDim();
  subdomsMech_ = ptPDE1_->getPDE_subdoms();

  if (subdomsMech_.GetSize()==1) {

    stiffTensorC[0][0] = parameterC_[0]; //c_11
    stiffTensorC[1][1] = parameterC_[0]; //c_11
    stiffTensorC[2][2] = parameterC_[1]; //c_33
    stiffTensorC[0][1] = parameterC_[2]; //c_12
    stiffTensorC[1][0] = parameterC_[2]; //c_12
    stiffTensorC[0][2] = parameterC_[3]; //c_13
    stiffTensorC[2][0] = parameterC_[3]; //c_13
    stiffTensorC[1][2] = parameterC_[3]; //c_13
    stiffTensorC[2][1] = parameterC_[3]; //c_13
    stiffTensorC[3][3] = parameterC_[4]; //c_44
    stiffTensorC[4][4] = parameterC_[4]; //c_44
    stiffTensorC[5][5] = 0.5*(parameterC_[0]-parameterC_[2]); //c_66

    if (dim==2) {
      stiffTensorC[0][0] = parameterC_[0]; //c_11
      stiffTensorC[2][2] = parameterC_[0]; //c_11
      stiffTensorC[1][1] = parameterC_[1]; //c_33
      stiffTensorC[0][2] = parameterC_[2]; //c_12
      stiffTensorC[2][0] = parameterC_[2]; //c_12
      stiffTensorC[0][1] = parameterC_[3]; //c_13
      stiffTensorC[1][0] = parameterC_[3]; //c_13
      stiffTensorC[1][2] = parameterC_[3]; //c_13
      stiffTensorC[2][1] = parameterC_[3]; //c_13
      stiffTensorC[3][3] = parameterC_[4]; //c_44
      stiffTensorC[5][5] = parameterC_[4]; //c_44
      stiffTensorC[4][4] = 0.5*(parameterC_[0]-parameterC_[2]); //c_66
    }

    piezoTensorC[1][3]=parameterC_[5]; //e_15
    piezoTensorC[0][4]=parameterC_[5]; //e_15
    piezoTensorC[2][0]=parameterC_[6]; //e_31
    piezoTensorC[2][1]=parameterC_[6]; //e_31
    piezoTensorC[2][2]=parameterC_[7]; // e_33

    if (dim==2) {
      piezoTensorC[2][3]=parameterC_[5]; //e_15
      piezoTensorC[0][5]=parameterC_[5]; //e_15
      piezoTensorC[1][0]=parameterC_[6]; //e_31
      piezoTensorC[1][2]=parameterC_[6]; //e_31
      piezoTensorC[1][1]=parameterC_[7]; // e_33
    }

    permTensorC[0][0] = parameterC_[8]; //eps_11
    permTensorC[1][1] = parameterC_[8]; //eps_11
    permTensorC[2][2] = parameterC_[9]; //eps_33

    if (dim==2) {
      permTensorC[0][0] = parameterC_[8]; //eps_11
      permTensorC[2][2] = parameterC_[8]; //eps_11
      permTensorC[1][1] = parameterC_[9]; //eps_33
    }

    // Get first regionId of mechanic pde
    RegionIdType actId = subdomsMech_[0];

    ptMaterialPiezo_[actId]->SetTensor(piezoTensorC, PIEZO_TENSOR, Global::IMAG);
    ptMaterialMech_[actId]->SetTensor(stiffTensorC, MECH_STIFFNESS_TENSOR, Global::IMAG);
    ptMaterialElec_[actId]->SetTensor(permTensorC, ELEC_PERMITTIVITY, Global::IMAG);

  }

  // Materialparameteridentifizierung fuer eine piezoelektrische Keramik
  // und eine Anpassschicht
  Matrix<Double> stiffTensorAdaptionC(6,6);

  if (subdomsMech_.GetSize()==2) {
    RegionIdType actId = subdomsMech_[0];
    RegionIdType actId1 = subdomsMech_[1];

    stiffTensorC[0][0] = parameterC_[0]; //c_11
    stiffTensorC[1][1] = parameterC_[0]; //c_11
    stiffTensorC[2][2] = parameterC_[1]; //c_33
    stiffTensorC[0][1] = parameterC_[2]; //c_12
    stiffTensorC[1][0] = parameterC_[2]; //c_12
    stiffTensorC[0][2] = parameterC_[3]; //c_13
    stiffTensorC[2][0] = parameterC_[3]; //c_13
    stiffTensorC[1][2] = parameterC_[3]; //c_13
    stiffTensorC[2][1] = parameterC_[3]; //c_13
    stiffTensorC[3][3] = parameterC_[4]; //c_44
    stiffTensorC[4][4] = parameterC_[4]; //c_44
    stiffTensorC[5][5] = 0.5*(parameterC_[0]-parameterC_[2]); //c_66


    // Set Lame parameters for adaption layer
    Double lambda2nu=parameterC_[11]+2*parameterC_[10];

    stiffTensorAdaptionC[0][0] = lambda2nu; //c_11
    stiffTensorAdaptionC[1][1] = lambda2nu; //c_11
    stiffTensorAdaptionC[2][2] = lambda2nu; //c_33
    stiffTensorAdaptionC[0][1] = parameterC_[11]; //c_12
    stiffTensorAdaptionC[1][0] = parameterC_[11]; //c_12
    stiffTensorAdaptionC[0][2] = parameterC_[11]; //c_13
    stiffTensorAdaptionC[2][0] = parameterC_[11]; //c_13
    stiffTensorAdaptionC[1][2] = parameterC_[11]; //c_13
    stiffTensorAdaptionC[2][1] = parameterC_[11]; //c_13
    stiffTensorAdaptionC[3][3] = parameterC_[10]; //c_44
    stiffTensorAdaptionC[4][4] = parameterC_[10]; //c_44
    stiffTensorAdaptionC[5][5] = parameterC_[10]; //c_66


    if (dim==2) { // also axi ...
      stiffTensorC[0][0] = parameterC_[0]; //c_11
      stiffTensorC[2][2] = parameterC_[0]; //c_11
      stiffTensorC[1][1] = parameterC_[1]; //c_33
      stiffTensorC[0][2] = parameterC_[2]; //c_12
      stiffTensorC[2][0] = parameterC_[2]; //c_12
      stiffTensorC[0][1] = parameterC_[3]; //c_13
      stiffTensorC[1][0] = parameterC_[3]; //c_13
      stiffTensorC[1][2] = parameterC_[3]; //c_13
      stiffTensorC[2][1] = parameterC_[3]; //c_13
      stiffTensorC[3][3] = parameterC_[4]; //c_44
      stiffTensorC[5][5] = parameterC_[4]; //c_44
      stiffTensorC[4][4] = 0.5*(parameterC_[0]-parameterC_[2]); //c_66

      // Set Lame parameters for adaption layer
      Double lambda2nu=parameterC_[11]+2*parameterC_[10];

      stiffTensorAdaptionC[0][0] = lambda2nu; //c_11
      stiffTensorAdaptionC[2][2] = lambda2nu; //c_11
      stiffTensorAdaptionC[1][1] = lambda2nu; //c_33
      stiffTensorAdaptionC[0][2] = parameterC_[11]; //c_12
      stiffTensorAdaptionC[2][0] = parameterC_[11]; //c_12
      stiffTensorAdaptionC[0][1] = parameterC_[11]; //c_13
      stiffTensorAdaptionC[1][0] = parameterC_[11]; //c_13
      stiffTensorAdaptionC[1][2] = parameterC_[11]; //c_13
      stiffTensorAdaptionC[2][1] = parameterC_[11]; //c_13
      stiffTensorAdaptionC[3][3] = parameterC_[10]; //c_44
      stiffTensorAdaptionC[5][5] = parameterC_[10]; //c_44
      stiffTensorAdaptionC[4][4] = parameterC_[10]; //c_66

    }

    piezoTensorC[1][3]=parameterC_[5]; //e_15
    piezoTensorC[0][4]=parameterC_[5]; //e_15
    piezoTensorC[2][0]=parameterC_[6]; //e_31
    piezoTensorC[2][1]=parameterC_[6]; //e_31
    piezoTensorC[2][2]=parameterC_[7]; // e_33

    if (dim==2) {
      piezoTensorC[2][3]=parameterC_[5]; //e_15
      piezoTensorC[0][5]=parameterC_[5]; //e_15
      piezoTensorC[1][0]=parameterC_[6]; //e_31
      piezoTensorC[1][2]=parameterC_[6]; //e_31
      piezoTensorC[1][1]=parameterC_[7]; // e_33
    }

    permTensorC[0][0] = parameterC_[8]; //eps_11
    permTensorC[1][1] = parameterC_[8]; //eps_11
    permTensorC[2][2] = parameterC_[9]; //eps_33

    if (dim==2) {
      permTensorC[0][0] = parameterC_[8]; //eps_11
      permTensorC[1][1] = parameterC_[9]; //eps_33
    }

    ptMaterialPiezo_[actId]->SetTensor(piezoTensorC, PIEZO_TENSOR, Global::IMAG);
    ptMaterialMech_[actId]->SetTensor(stiffTensorC, MECH_STIFFNESS_TENSOR, Global::IMAG);
    ptMaterialElec_[actId]->SetTensor(permTensorC, ELEC_PERMITTIVITY, Global::IMAG);

    ptMaterialMech_[actId1]->SetTensor(stiffTensorAdaptionC, MECH_STIFFNESS_TENSOR, Global::IMAG);
  }

  // Materialparameteridentifizierung fuer eine piezoelektrische Keramik,
  //  Anpassschicht und backing
  if (subdomsMech_.GetSize()==3) {
    RegionIdType actId = subdomsMech_[0];
    RegionIdType actId1 = subdomsMech_[1];
    RegionIdType actId2 = subdomsMech_[2];

    stiffTensorC[0][0] = parameterC_[0]; //c_11
    stiffTensorC[1][1] = parameterC_[0]; //c_11
    stiffTensorC[2][2] = parameterC_[1]; //c_33
    stiffTensorC[0][1] = parameterC_[2]; //c_12
    stiffTensorC[1][0] = parameterC_[2]; //c_12
    stiffTensorC[0][2] = parameterC_[3]; //c_13
    stiffTensorC[2][0] = parameterC_[3]; //c_13
    stiffTensorC[1][2] = parameterC_[3]; //c_13
    stiffTensorC[2][1] = parameterC_[3]; //c_13
    stiffTensorC[3][3] = parameterC_[4]; //c_44
    stiffTensorC[4][4] = parameterC_[4]; //c_44
    stiffTensorC[5][5] = 0.5*(parameterC_[0]-parameterC_[2]); //c_66


    // Set Lame parameters for adaption layer
    Double lambda2nu=parameterC_[11]+2*parameterC_[10];
    Matrix<Double> stiffTensorAdaptionC(6,6);
    Matrix<Double> stiffTensorBackingC(6,6);

    stiffTensorAdaptionC[0][0] = lambda2nu;
    stiffTensorAdaptionC[1][1] = lambda2nu;
    stiffTensorAdaptionC[2][2] = lambda2nu;
    stiffTensorAdaptionC[0][1] = parameterC_[11];
    stiffTensorAdaptionC[1][0] = parameterC_[11];
    stiffTensorAdaptionC[0][2] = parameterC_[11];
    stiffTensorAdaptionC[2][0] = parameterC_[11];
    stiffTensorAdaptionC[1][2] = parameterC_[11];
    stiffTensorAdaptionC[2][1] = parameterC_[11];
    stiffTensorAdaptionC[3][3] = parameterC_[10];
    stiffTensorAdaptionC[4][4] = parameterC_[10];
    stiffTensorAdaptionC[5][5] = parameterC_[10];

    lambda2nu=parameterC_[13]+2*parameterC_[12];

    stiffTensorBackingC[0][0] = lambda2nu;
    stiffTensorBackingC[1][1] = lambda2nu;
    stiffTensorBackingC[2][2] = lambda2nu;
    stiffTensorBackingC[0][1] = parameterC_[13];
    stiffTensorBackingC[1][0] = parameterC_[13];
    stiffTensorBackingC[0][2] = parameterC_[13];
    stiffTensorBackingC[2][0] = parameterC_[13];
    stiffTensorBackingC[1][2] = parameterC_[13];
    stiffTensorBackingC[2][1] = parameterC_[13];
    stiffTensorBackingC[3][3] = parameterC_[12];
    stiffTensorBackingC[4][4] = parameterC_[12];
    stiffTensorBackingC[5][5] = parameterC_[12];


    if (dim==2) { // also axi ...
      stiffTensorC[0][0] = parameterC_[0]; //c_11
      stiffTensorC[2][2] = parameterC_[0]; //c_11
      stiffTensorC[1][1] = parameterC_[1]; //c_33
      stiffTensorC[0][2] = parameterC_[2]; //c_12
      stiffTensorC[2][0] = parameterC_[2]; //c_12
      stiffTensorC[0][1] = parameterC_[3]; //c_13
      stiffTensorC[1][0] = parameterC_[3]; //c_13
      stiffTensorC[1][2] = parameterC_[3]; //c_13
      stiffTensorC[2][1] = parameterC_[3]; //c_13
      stiffTensorC[3][3] = parameterC_[4]; //c_44
      stiffTensorC[5][5] = parameterC_[4]; //c_44
      stiffTensorC[4][4] = 0.5*(parameterC_[0]-parameterC_[2]); //c_66
    }

    piezoTensorC[1][3]=parameterC_[5]; //e_15
    piezoTensorC[0][4]=parameterC_[5]; //e_15
    piezoTensorC[2][0]=parameterC_[6]; //e_31
    piezoTensorC[2][1]=parameterC_[6]; //e_31
    piezoTensorC[2][2]=parameterC_[7]; // e_33

    if (dim==2) {
      piezoTensorC[2][3]=parameterC_[5]; //e_15
      piezoTensorC[0][5]=parameterC_[5]; //e_15
      piezoTensorC[1][0]=parameterC_[6]; //e_31
      piezoTensorC[1][2]=parameterC_[6]; //e_31
      piezoTensorC[1][1]=parameterC_[7]; // e_33
    }

    permTensorC[0][0] = parameterC_[8]; //eps_11
    permTensorC[1][1] = parameterC_[8]; //eps_11
    permTensorC[2][2] = parameterC_[9]; //eps_33

    if (dim==2) {
      permTensorC[0][0] = parameterC_[8]; //eps_11
      permTensorC[1][1] = parameterC_[9]; //eps_33
    }



    ptMaterialPiezo_[actId]->SetTensor(piezoTensorC, PIEZO_TENSOR, Global::IMAG);
    ptMaterialMech_[actId]->SetTensor(stiffTensorC, MECH_STIFFNESS_TENSOR, Global::IMAG);
    ptMaterialElec_[actId]->SetTensor(permTensorC, ELEC_PERMITTIVITY, Global::IMAG);

    ptMaterialMech_[actId1]->SetTensor(stiffTensorAdaptionC, MECH_STIFFNESS_TENSOR, Global::IMAG);
    ptMaterialMech_[actId2]->SetTensor(stiffTensorBackingC, MECH_STIFFNESS_TENSOR, Global::IMAG);
  }

} // end updateComplexMaterialData

void piezoParamIdent::setNewParameterSet(Vector<Double> & par,
    Vector<Double> & par_new, Vector<Double> & scaling_, Double & theta,
    Vector<Double> & uStep, Vector<UInt> & whichParameterToUpdate_) {
  UInt helpInd=0;
  for (UInt i=0; i<nrParameter_; i++) {
    if (whichParameterToUpdate_[i]==1) {
      par_new[i]=par[i]+(1.0/scaling_[i])*theta*uStep[helpInd];
      helpInd++;
    }
  }

} // end setNewParameterSet


void piezoParamIdent::computeScaling() {

  // for most of the methods we scale the sought-for quantities
  for (UInt i=0; i<parameter_.GetSize(); i++)
    if (parameter_[i]!=0.0)
      scaling_[i]= std::abs(1.0/parameter_[i]);
    else
      scaling_[i]=1.0;

  for (UInt i=0; i<parameterC_.GetSize(); i++)
    if (parameterC_[i]!=0.0)
      scalingC_[i]= std::abs(1.0/parameterC_[i]);
    else
      scalingC_[i]=1.0;

}//end compute scaling


void piezoParamIdent::writeTensorsInFile() {

  Matrix<Complex> piezoMatC, stiffMatC, permMatC;

  // Get first regionId of mechanic pde
  RegionIdType actId = subdomsMech_[0];

  if(ptPDE1_->getPDE_spaceDim()==3){
    ptMaterialPiezo_[actId]->GetTensor(piezoMatC, PIEZO_TENSOR, Global::COMPLEX, FULL);
    ptMaterialMech_[actId]->GetTensor(stiffMatC, MECH_STIFFNESS_TENSOR, Global::COMPLEX, FULL);
    ptMaterialElec_[actId]->GetTensor(permMatC, ELEC_PERMITTIVITY, Global::COMPLEX, FULL);
  }
  else {
    ptMaterialPiezo_[actId]->GetTensor(piezoMatC, PIEZO_TENSOR, Global::COMPLEX, AXI);
    ptMaterialMech_[actId]->GetTensor(stiffMatC, MECH_STIFFNESS_TENSOR, Global::COMPLEX, AXI);
    ptMaterialElec_[actId]->GetTensor(permMatC, ELEC_PERMITTIVITY, Global::COMPLEX, AXI);
  }

  Matrix<Complex> sE;
  if(ptPDE1_->getPDE_spaceDim()==3)
    sE.Resize(6,6);
  else
    sE.Resize(4,4);

  sE.Init();
  for (UInt i=0; i<sE.GetNumRows(); i++)
    sE[i][i]=Complex(1.0, 0.0);

#ifdef USE_LAPACK

  lapackSysMatType LAPACK_SYS_MAT_TYPE = ZGESV;
  stiffMatC.solveWithLapack(sE,LAPACK_SYS_MAT_TYPE);

#endif

  *allTensors_<<"Mechanical Modulus Tensor cE = "<<std::endl;
  *allTensors_<<stiffMatC<<std::endl;

  *allTensors_<<"Mechanical Elasticity Tensor sE = "<<std::endl;
  *allTensors_<<sE<<std::endl;

  *allTensors_<<"Piezoelectric Coupling Tensor e = "<<std::endl;
  *allTensors_<<piezoMatC<<std::endl;

  Matrix<Complex> d(piezoMatC.GetNumRows(), piezoMatC.GetNumCols());
  d.Init();
  d=piezoMatC*sE;
  *allTensors_<<"Piezoelectric Coupling Tensor d = "<<std::endl;
  *allTensors_<<d<<std::endl;

  *allTensors_<<"Permittivity Tensor epsS = "<<std::endl;
  *allTensors_<<permMatC<<std::endl;

  Matrix<Complex> epsT(permMatC.GetNumRows(),permMatC.GetNumRows() );
  Matrix<Complex> epsTrans;
  piezoMatC.Transpose(epsTrans);

  epsT.Init();

  epsT = permMatC + d*epsTrans;

  *allTensors_<<"Permittivity Tensor epsT = "<<std::endl;
  *allTensors_<<epsT<<std::endl;

  Complex influenceConstant = Complex(1.0/8.8542e-12, 0.0);
  permMatC *= influenceConstant;
  epsT *= influenceConstant;

  *allTensors_<<"Relative Permittivity Tensor epsS = "<<std::endl;
  *allTensors_<<permMatC<<std::endl;

  *allTensors_<<"Relative Permittivity Tensor epsT = "<<std::endl;
  *allTensors_<<epsT<<std::endl;

}


  // ************************
  //   ComputeNextFrequency
  // ************************
  Double piezoParamIdent::ComputeNextFrequency( UInt freqIndex ) const {


    Double retFreq = 0.0;

    // Check for single step
    if ( nrfreq_ == 1 ) {
      retFreq = startFreq_;
    }
    else {

      switch( samplingType_ ) {

        // Linear sampling
      case LINEAR_SAMPLING:
        retFreq = (freqIndex - 1) * (stopFreq_ - startFreq_) /
          (Double)( nrfreq_ - 1 ) + startFreq_;
        break;

        // Logarithmic sampling
      case LOG_SAMPLING:
        {
          Double fac = stopFreq_ / startFreq_;
          fac = std::pow( fac, (Double)(freqIndex - 1) / (nrfreq_ - 1) );
          retFreq = startFreq_ * fac;
        }
        break;

        // Reverse logarithmic sampling
      case REVERSE_LOG_SAMPLING:
        {
          Double fac = stopFreq_ / startFreq_;
          fac = std::pow( fac, (Double)(nrfreq_ - freqIndex)
                          / (nrfreq_ - 1));
          retFreq = stopFreq_ + startFreq_ * ( 1.0 - fac );
        }
        break;

        // Something's wrong
      default:
        std::string damp;
        Enum2String( samplingType_, damp );
        EXCEPTION( "HarmonicDriver::ComputeNextFrequency: '"
                   << damp
                   << "' is not supported as sampling type" );
      }
    }

    return retFreq;
  }


} // end namespace CoupledField


