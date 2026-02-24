// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     BaseSNGR.cc
 *       \brief    <Description>
 *
 *       \date     Nov 01, 2017
 *       \author   r.krusche, Stefan Schoder, Michael Weitz
 */
//================================================================================================


#include <Filters/SynteticSources/BaseSNGR.hh>
#include <random>
#include <chrono>

namespace CFSDat{

SNGRFilter::SNGRFilter(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan)
        :BaseFilter(numWorkers,config,resMan){
  this->filtStreamType_ = FIFO_FILTER;

  globalTemp_ = 1.0;
  globalDensity_ = 1.0;

  //INPUT from CFD
  this->inTKE_ = config->Get("TKE/resultName")->As<std::string>();
  this->inTEF_ = config->Get("TEF/resultName")->As<std::string>();
  this->inVelocity_ = config->Get("meanVelocity/resultName")->As<std::string>();
  this->inDensity_ = "const1234";
  this->inTemp_ = "const1234";


  //DECLARE Input
  upResNames.insert(inTKE_);
  upResNames.insert(inVelocity_);

  if(config->Has("localDensity/resultName")){
    this->inDensity_ = config->Get("localDensity/resultName")->As<std::string>();
    upResNames.insert(inDensity_);
  } else if (config->Has("localDensity/value")) {
    globalDensity_ = config->Get("localDensity/value")->As<Double>();
  } else {
    EXCEPTION("Specify the resultName or a constant value of the density.");
  }

  if(config->Has("localTemp/resultName")){
    this->inTemp_ = config->Get("localTemp/resultName")->As<std::string>();
    upResNames.insert(inTemp_);
  } else if (config->Has("localTemp/value")) {
    globalTemp_ = config->Get("localTemp/value")->As<Double>();
  } else {
    EXCEPTION("Specify the resultName or a constant value of the temperature.");
  }

  if (this->inDensity_ == "const1234" && this->inTemp_ != "const1234"){
    EXCEPTION("Both Density and Temperature must be constant or not.");
  } else if (this->inDensity_ != "const1234" && this->inTemp_ == "const1234"){
    EXCEPTION("Both Density and Temperature must be constant or not.");
  }

  upResNames.insert(inTEF_);

  //OUTPUT mandatory
  this->outName_ = config->Get("output/resultName")->As<std::string>();
  filtResNames.insert(outName_);
  //OUTPUT optional
  if (config->Has("intermediateResult")) {
    this->interName_ = config->Get("intermediateResult/resultName")->As<std::string>();
    filtResNames.insert(interName_);
  }

  // using the TKE-criterion, the source region in the CFD-domain used for reconstructing turbulent velocity, is controlled
  // typical value in literature is 10%
  // only nodes with TKE exceeding TKEcrit_*TKEmax are used for reconstruction (typically TKEcrit_ = 0.1)
  this->TKEcrit_ = config->Get("tkeCriterion")->As<Double>();

  // if specific length of output signal is needed, define it by "signalLength"-tag
  this->sigLength_ = config->Get("signalLength")->As<Double>();
  this->fL_ = config->Get("lengthScaleFactor")->As<Double>();
  this->ft_ = config->Get("timeScaleFactor")->As<Double>();
  this->fa_ = config->Get("angularFreqFactor")->As<Double>();
  // number of Modes used for reconstruction
  this->numModes_ = config->Get("numOfModes")->As<UInt>();
  // equidstant or logarithmmic distribution of Modes in energy spectrum
  this->incrModes_ = config->Get("incrementModes")->As<std::string>();
  this->method_ = config->Get("method")->As<std::string>();

  // frequency, the user wants to dissolve acoustically
  this->maxFreq_ = config->Get("frequencyBounds/maxFreq")->As<UInt>(); //TODO is this important
  this->minFreq_ = config->Get("frequencyBounds/minFreq")->As<UInt>(); //TODO is this important

  // for maxWN and minWN percentage of the peak wave number are expected
  // typical values in literature: minWN = 0.1*peakWN; maxWN = 10*peakWN
  this->maxWN_ = config->Get("waveNumBounds/maxWN")->As<Double>();
  this->minWN_ = config->Get("waveNumBounds/minWN")->As<Double>();

  // enable ensemble-average at microphone location
  this->ensemble_ = config->Get("ensembles")->As<UInt>();


}

void SNGRFilter::PrepareCalculation(){

  //input grid
  inGrid_ = resultManager_->GetExtInfo(upResIds[0])->ptGrid;

  // set global parameters, read RANS data, etc.
  this->PrepareMethod();

  // chose method
  if(method_ == "BaillySNGR"){
    // go with Bailly, 1999
    std::cout << "Start reconstruction of synthetic turbulent velocity field using SNGR by Bailly & Juve, 1999." << std::endl;
    std::cout << "According to publication: A Stochastic Approach To Compute Subsonic Noise Using Linearized Euler's Equations." << std::endl;
    std::cout << "Doing " << ensemble_ << " independent reconstruction(s)." << std::endl;

    // TODO remove debug output
    std::cout << "readInResultName, TKE: " << inTKE_ << std::endl;

    // Initialize SNGR method
    this->InitArraysBailly();
    std::cout << "enter UpdateResultBailly" << std::endl;
    this->InitResultMethodBailly();
    std::cout << "enter SetTimelineMethodBailly" << std::endl;

  }
  else if(method_ == "BillsonSNGR"){
    // go with Billson, 2003
    std::cout << "Start reconstruction of synthetic turbulent velocity field using SNGR by Billson, Eriksson & Davidson, 2003." << std::endl;
    std::cout << "According to publication: Jet Noise Prediction Using Stochastic Turbulence Modeling." << std::endl;
    std::cout << "Doing " << ensemble_ << " independent reconstruction(s)." << std::endl;

    // Initialize SNGR method
  }
  else if(method_ == "Lafitte"){
    // go with Lafitte 2014
    std::cout << "Start reconstruction of synthetic turbulent velocity field using SNGR by Lafitte, Le Garrec, Bailly & Laurendeau, 2014." << std::endl;
    std::cout << "According to publication: Turbulence Generation from a Sweeping-Based Stochastic Model." << std::endl;
    std::cout << "Doing " << ensemble_ << " independent reconstruction(s)." << std::endl;

    this->InitArraysLafitte();
  }
  else{
    EXCEPTION("This is not a valid method.")
  }
}


bool SNGRFilter::UpdateResults(std::set<uuids::uuid>& upResults){
  resultManager_->ActivateResult(tkeId_);
  resultManager_->ActivateResult(tefId_);
  resultManager_->ActivateResult(velocityId_);
  if (this->inDensity_ != "const1234"){
    resultManager_->ActivateResult(densityId_);
  }
  if (this->inTemp_ != "const1234"){
    resultManager_->ActivateResult(temperatureId_);
  }
  resultManager_->DeactivateResult(outId_);
  if(params_->Has("intermediateResult")) resultManager_->DeactivateResult(interId_);

  //now we call for upstream data in each source
  CF::StdVector< shared_ptr<BaseFilter> >::iterator srcIter =  sources_.Begin();
  for(; srcIter != sources_.End() ; srcIter++){
    // should we check here anything for success?
    (*srcIter)->Run();
  }

  // Now we compute the downstream results for the current value
  std::set<uuids::uuid>::iterator oIter = upResults.begin();
  for(; oIter != upResults.end(); ++oIter){
    std::string resName = resultManager_->GetExtInfo(*oIter)->resultName;
//    std::cout<<"oIter: "<<resName<<std::endl;

    if (resName == resultManager_->GetResultName(outId_)) {

      const UInt actStepIndex = resultManager_->GetStepIndex(*oIter);
//      std::cout << "actStepIndex " << actStepIndex << std::endl;
      this->SetTimelineMethodBailly(actStepIndex);

      // calc current time step
//      std::cout<<"outID: "<<resultManager_->GetResultName(outId_)<<std::endl;

      //equation numbers for this result //VECTOR result
      //here, those equation numbers are equal for all in/out results
      Vector<Double>& returnVec = GetOwnResultVector<Double>(*oIter);
//      UInt size = returnVec.GetSize()/3;

      for(UInt i=0;i<idsNodesToProcess_.GetSize();++i){
          UInt j = idsNodesToProcess_[i];
          // For input values use j-index, for output quantities i-index
          returnVec[3*j] = turbReconstVelocity_[3*i];
          returnVec[3*j + 1] = turbReconstVelocity_[3*i + 1];
          returnVec[3*j + 2] = turbReconstVelocity_[3*i + 2];
      }

      resultManager_->SetResultVecUpToDate(*oIter,true);

    } else if (resName == resultManager_->GetResultName(interId_)) {
      Vector<Double>& returnVecInter = GetOwnResultVector<Double>(*oIter);
      std::cout<<"InterID: "<<resultManager_->GetResultName(interId_)<<std::endl;
      //RETURN OTHER OUTPUT intermediate result // Dummy SCALAR output
      for(UInt i=0;i<idsNodesToProcess_.GetSize();++i){
          UInt j = idsNodesToProcess_[i];

          returnVecInter[j] = reconstTKE_[i];
      }
      resultManager_->SetResultVecUpToDate(*oIter,true);

    }

  }

  return true;
}


ResultIdList SNGRFilter::SetUpstreamResults(){
  //RESULT DEFINITION
  ResultIdList generated = SetDefaultUpstreamResults();
  tkeId_ = upResNameIds[inTKE_];
  velocityId_ = upResNameIds[inVelocity_];
  tefId_ = upResNameIds[inTEF_]; //TODO is this and do we ahve to do this?
  if (this->inDensity_ != "const1234"){
    densityId_ = upResNameIds[inDensity_];
  }
  if (this->inTemp_ != "const1234"){
    temperatureId_ = upResNameIds[inTemp_];
  }

  return generated;
}

void SNGRFilter::AdaptFilterResults(){
  //We should check some validity...
  for(UInt aRes = 0; aRes < filterResIds.GetSize(); aRes++){
    //we can almost copy everything from time input (also scalar, etc.)
    if(resultManager_->GetResultName(filterResIds[aRes]) == this->outName_ ) {
      outId_ = filterResIds[aRes];
//      std::cout<<"RES1:"<<resultManager_->GetResultName(filterResIds[aRes])<<std::endl;
      resultManager_->CopyResultData(velocityId_,outId_);
      resultManager_->SetValid(outId_);
    } else if (resultManager_->GetResultName(filterResIds[aRes]) == this->interName_ ) {
    //REGISTER INTERMEDIATE result //SCALAR copied
      interId_  = filterResIds[aRes];
//      std::cout<<"RES0:"<<resultManager_->GetResultName(filterResIds[aRes])<<std::endl;
      resultManager_->CopyResultData(tkeId_,interId_);
      resultManager_->SetValid(interId_);
    }
  }
}

// ===================================================================== //
//                            Utility functions                          //
// ===================================================================== //
void SNGRFilter::GetTkeThreshold(){
  // set TKE threshold
  Double maxTKE=0.0;
  for(UInt m=0; m<numNodes_; ++m){
    if(TKE_[m]>=maxTKE) maxTKE=TKE_[m];
  }
  minTKE_ = TKEcrit_*maxTKE;
  // get node ids, that pass the tke criterion

  for(UInt i=0; i<numNodes_; ++i){
    if(TKE_[i]>=minTKE_){
			idsNodesToProcess_.Push_back(i);
    }
    else{
      idsNodesToProcessOnlyMeanVelocity_.Push_back(i); // nodes where only mean velocity is considered
    }
  }
}

void SNGRFilter::SetRandVectors(UInt k, UInt j,Double rmsOfVelFluct,Vector<Double> kn){
  // Set up random number generators
  UInt seed_uniform = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937 generator_uniform (seed_uniform);
  std::uniform_real_distribution<double> uniform01(0.0, 1.0);

  Double theta = 0.5 * (1 - cos(M_PI * uniform01(generator_uniform)));
  Double phi = 2 * M_PI * uniform01(generator_uniform);


  // phase and direction of mode.
  //srand(time(NULL));      // initialize random seed.
  //Double theta = (M_PI - 0.0)*((Double)rand()/(Double)RAND_MAX); // + 0.0;        // polar angle of wave vector, random value between 0 ... pi
  //srand(time(NULL));
  //Double phi = (2.0*M_PI - 0.0)*((Double)rand()/(Double)RAND_MAX); // + 0.0;        // azimuthal angle of wave vector, random value between 0 ... 2pi
  srand (time(NULL));
  Double alpha = (2*M_PI - 0.0)*((Double)rand()/(Double)RAND_MAX); // + 0.0;      // angle for directional vector of current mode, random value between 0 ... 2pi
  srand (time(NULL)+33.0);
  phase_[numModes_*k+j] = (2.0*M_PI - 0.0)*((Double)rand()/(Double)RAND_MAX); // + 0.0;      // phase of current mode, random value between 0 ... 2pi


  // construct a trivial random generator engine from a time-based seed:
  UInt seed_normal = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator_normal (seed_normal);
  // set angular frequency for j-th mode
  Double omegaExpectValue = rmsOfVelFluct * kn[j];
  std::normal_distribution<double> distribution(omegaExpectValue,omegaExpectValue);
  omega_[numModes_*k+j] = distribution(generator_normal);

  //TODO revisit - set wave vector
  // magnitude of wave vector is wave number of current mode
  UInt g = 3*numModes_*k+3*j; // index in array for wave and direction vector
  waveVec_[g] = kn[j]*cos(theta)*cos(phi);
  waveVec_[g+1] = kn[j]*sin(phi);
  waveVec_[g+2] = -kn[j]*sin(theta)*cos(phi);
  // set direction vector, perpendicular to wave vector, ensures incompressible velocity field
  dirVec_[g] = -sin(phi)*cos(alpha)*cos(theta)+sin(alpha)*sin(theta);
  dirVec_[g+1] = cos(phi)*cos(alpha);
  dirVec_[g+2] = sin(phi)*cos(alpha)*sin(theta)+sin(alpha)*cos(theta);

  // check if sufficiently perpendicular
  Double perp = waveVec_[g]*dirVec_[g]+waveVec_[g+1]*dirVec_[g+1]+waveVec_[g+2]*dirVec_[g+2];
  if(abs(perp) > 1.0e-9){
     perpFAIL_++;
//     std::cout << "WARNING! In node " << i <<" wave number vector not sufficiently perpendicular to direction vector!" << std::endl;
//     std::cout << "result of dot product was " << perp << std::endl;
  }
}

void SNGRFilter::SetRandVectorsBillson(){
  // TODO check if this could be included in another function
}

void SNGRFilter::SetVelocityAmplitude(UInt k, UInt j, Double peakWN,
    Double deltaWN, Double kolmogorovWN, Double rmsOfVelFluct,Vector<Double> kn){
  // compute Karman-Pao spectrum for isotropic turbulence
  Double E_kn = A_*(pow(rmsOfVelFluct,2.0)/peakWN)*pow(kn[j]/peakWN,4.0)*exp(-2.0*pow(kn[j]/kolmogorovWN,2.0))*1.0/(pow((1.0+pow(kn[j]/peakWN,2.0)),(17.0/6.0)));
  // amplitude for current mode
//  if(incrModes_ == "logarithmic") deltaWN_=dWN_[j];
  velAmplitude_[numModes_*k+j] = sqrt(E_kn*deltaWN);
  reconstTKE_[k] += E_kn*deltaWN;

//      // confirm that reconstructed velocity field is consistent with TKE from RANS analysis.
//      // literature gives reason to assume that this test fails often
//      Double diffTKE = abs(TKE_[i] - reconstTKE_[i]);
//      if(diffTKE > 1e-6){
//        std::cout << "WARNING! Deviation of reconstructed and input TKE in node " << i <<" too big!" << std::endl;
//        tkeFAIL_++;
//      }
}

void SNGRFilter::PrepareMethod(){

  PtrParamNode stepNode = params_->GetParent()->Get("stepValueDefinition");
  if(stepNode->Has("startStop")){
    PtrParamNode stNode = stepNode->Get("startStop");
    numSteps_ = stNode->Get("numSteps")->Get("value")->As<UInt>();
    deltaT_    = stNode->Get("delta")->Get("value")->As<Double>();
  } else {
    EXCEPTION("DEFINE START/STOP in PIPELINE!")
  }
  // compute time step, factor 30 for resolving a smooth curve
  maxFreq_ = 1.0/(2*deltaT_);
  minFreq_ = 1.0/((numSteps_-1)*deltaT_);
  std::cout << "frequency resolution as follows, maxFreq: " << maxFreq_ << ", minFreq: " << minFreq_ << std::endl;
//  deltaT_ = 1.0/(30*maxFreq_);
  // set length of output signal, factor 1.5 for safety
//  if (sigLength_ == 0.0){
//    sigLength_ = 1.5/(minFreq_);
//  }

  // time values
//  numSteps_ = ceil(sigLength_/deltaT_);
//  globalStepValueMap_.Resize(numSteps_);
//  for (UInt idx=0; idx<numSteps_; ++idx){
//    globalStepValueMap_[idx] = idx*deltaT_;
//  }

  // TODO remove debug output
  std::cout << "numSteps_: " << numSteps_ << std::endl;
  std::cout << "sigLength_: " << sigLength_ << std::endl;
  std::cout << "deltaT_: " << deltaT_ << std::endl;
//  std::cout << "length globalStepValueMap_: " << globalStepValueMap_.GetSize() << std::endl;
  std::cout << "stepValues(0) " << globalStepValueMap_[0] << std::endl;
  std::cout << "stepValues(1) " << globalStepValueMap_[1] << std::endl;
  std::cout << "stepValues(end-1) " << globalStepValueMap_[numSteps_-2] << std::endl;
  std::cout << "stepValues(end) " << globalStepValueMap_[numSteps_-1] << std::endl;

  // get grid
  StdVector<int> volRegionIds;
  inGrid_->GetVolRegionIds(volRegionIds);
  numNodes_ = inGrid_->GetNumNodes(volRegionIds);
  // get result for result id // TODO evtl. Vector<Double>& als Datentyp vor den Namen wieder einfügen.
  CF::StdVector<UInt> numEqation;

  TKE_ = GetUpstreamResultVector<Double>(tkeId_,numEqation);
  meanVelocity_ = GetUpstreamResultVector<Double>(velocityId_,numEqation);
  if (this->inDensity_ != "const1234"){
    localDensity_ = GetUpstreamResultVector<Double>(densityId_,numEqation);
  }
  if (this->inTemp_ != "const1234"){
    localTemp_ = GetUpstreamResultVector<Double>(temperatureId_,numEqation);
  }
  TEF_ = GetUpstreamResultVector<Double>(tefId_,numEqation);

  // get tke threshold, actually gets the ids of all nodes that pass the tke criterion
  this->GetTkeThreshold();
}

// ===================================================================== //
//                         Method Block functions                        //
// ============================ SNGR Bailly ============================ //
// ===================================================================== //
void SNGRFilter::InitArraysBailly(){
  // initialize
  UInt numNodesToProcess = idsNodesToProcess_.GetSize();
//  Vector<Double> turbLengthScale_(numModes_*nodesToProcess);
//  Vector<Double> peakWN(numNodes_);
  velAmplitude_.Resize(numModes_*numNodesToProcess);
  reconstTKE_.Resize(numNodesToProcess);
  omega_.Resize(numModes_*numNodesToProcess);
  phase_.Resize(numModes_*numNodesToProcess);
  waveVec_.Resize(numModes_*numNodesToProcess*3);
  dirVec_.Resize(numModes_*numNodesToProcess*3);
}

void SNGRFilter::InitResultMethodBailly(){
//#pragma omp parallel for
  for(UInt k=0; k<idsNodesToProcess_.GetSize(); ++k){
      UInt i = idsNodesToProcess_[k];
      // dissipation rate and viscosity (Sutherland)
      Double TDR;
      if(inTEF_ == "Turbulence_Eddy_Fr9"){
        TDR = c_mu_*TKE_[i]*TEF_[i];
      }
      else{
        TDR = TEF_[i];
      }
      Double kinVisco = (18.132941775e-6*(291.15+120.0)/(globalTemp_+120.0)*pow(globalTemp_/291.15,(3.0/2.0)))/globalDensity_;
      if (this->inDensity_ != "const1234" && this->inTemp_ != "const1234"){
        kinVisco = (18.132941775e-6*(291.15+120.0)/(localTemp_[i]+120.0)*pow(localTemp_[i]/291.15,(3.0/2.0)))/localDensity_[i];
      }

      // turbulent length scale
      Double rmsOfVelFluct = sqrt(2.0/3.0*TKE_[i]);
      //turbLengthScale_[i] = fL*pow(rmsOfVelFluct,3.0)/TDR; //
      Double turbLengthScale = fL_*pow(c_mu_,0.75)*pow(TKE_[i],1.5)/TDR;
      // Kolmogorov wave number
      Double kolmogorovWN = pow(TDR,(1.0/4.0))*1.0/(pow(kinVisco,(3.0/4.0)));
      // essential locations in engery spectrum
      Double peakWN = 9.0*M_PI/55.0*A_/turbLengthScale;
      Double kmax = maxWN_ * peakWN;
      Double kmin = minWN_ * peakWN;

      // compute wave number increments
      Vector<Double> kn;  // array of all wave number values
      kn.Resize(numModes_);
//      if(incrModes_ == "logarithmic"){
//        deltaWN_ = (log(kmax)-log(kmin))/(numModes_-1);
//        for(UInt p=0;p<numModes_;++p){
//          kn[p] = exp(log(kmin)+p*deltaWN_);
//          if(p==0) dWN_[0] = kmin;
//          else dWN_[p] = kn[p]-kn[p-1];
//        }
//      } else {
          Double deltaWN = (kmax-kmin)/(numModes_-1);
          for(UInt p=0;p<numModes_;++p){
            kn[p] = kmin + p*deltaWN;
          }
//      }
      // mode loop
      for(UInt j=0; j<numModes_; ++j){
        // set engergy and velocity amplitude for j-th mode
        this->SetVelocityAmplitude(k,j,peakWN,deltaWN,kolmogorovWN,rmsOfVelFluct,kn);

         // set wave and direction vector for j-th mode
        this->SetRandVectors(k,j,rmsOfVelFluct,kn);
      }
  }
}

void SNGRFilter::SetTimelineMethodBailly(UInt k){
  // time loop
//  for(UInt k=0; k<numSteps_; ++k){
    // initialize output
    turbReconstVelocity_.Resize(idsNodesToProcess_.GetSize()*3);
    turbReconstVelocity_.Init();
//#pragma omp parallel for
    for(UInt p=0; p<idsNodesToProcess_.GetSize(); ++p){
      UInt i = idsNodesToProcess_[p];
      // how to get a coordinate
      CF::Vector<Double> pCoord;
      inGrid_->GetNodeCoordinate3D(pCoord,i+1);
      // mode loop
      for(UInt j=0; j<numModes_; ++j){
        UInt g = numModes_*p*3+j*3;  // index in waveVec and dirVec arrays
        Double kx = waveVec_[g]*(pCoord[0]-meanVelocity_[i*3]*globalStepValueMap_[k]) + waveVec_[g+1]*(pCoord[1]-meanVelocity_[i*3+1]*globalStepValueMap_[k]) + waveVec_[g+2]*(pCoord[2]-meanVelocity_[i*3+2]*globalStepValueMap_[k]);
        Double sumElem = 2.0*velAmplitude_[numModes_*p+j]*cos(kx + phase_[numModes_*p+j] + fa_*omega_[numModes_*p+j]*globalStepValueMap_[k]);
        turbReconstVelocity_[p*3] += sumElem*dirVec_[g];
        turbReconstVelocity_[p*3+1] += sumElem*dirVec_[g+1];
        turbReconstVelocity_[p*3+2] += sumElem*dirVec_[g+2];
      }
    }

    // write node results to file
    std::cout << "Writing reconstructed velocity field to file." << std::endl;

/*  std::cout << "The reconstruction was carried out for " << flg_ << " of " << numNodes_ << " nodes, " << flg_/numNodes_*100.0 << "%." << std::endl;
    std::cout << "The remaining " << numNodes_-flg_ << "(" << (numNodes_-flg_)/numNodes_*100.0 << "%) nodes did not meet the TKE-Criterion." << std::endl;
    if(tkeFAIL_!=0){
      std::cout << "There has been a problem with the reconstructed TKE in " << tkeFAIL_ << " nodes." << std::endl;
    }
    if(perpFAIL_!=0){
      std::cout << "For " << perpFAIL_ << " nodes the wave vector and direction vector where not sufficiently perpendicular." << std::endl;
    }
*/

//  }
  std::cout << "Reconstruction of turbulent velocity field finished" << std::endl;
}

// ===================================================================== //
//                         Method Block functions                        //
// ============================ Billson SNGR =========================== //
// ===================================================================== //
void SNGRFilter::InitArraysBillson(){
	// initialize
	UInt numNodesToProcess = idsNodesToProcess_.GetSize();
  velAmplitude_.Resize(numModes_*numNodesToProcess);
  waveNumIncrements_.Resize(numModes_*numNodesToProcess);
  reconstTKE_.Resize(numNodesToProcess);
}


void SNGRFilter::UpdateResultBillson(){
	// node loop TODO pragma omp parallel for
  UInt i = 0;
  for(UInt k=0; k<idsNodesToProcess_.GetSize(); ++k){
      i = idsNodesToProcess_[k];
      // dissipation rate and viscosity (Sutherland)
      Double TDR;
      if(inTEF_ == "Turbulence_Eddy_Fr9"){
        TDR = c_mu_*TKE_[i]*TEF_[i];
      }
      else{
        TDR = TEF_[i];
      }
      Double kinVisco = (18.132941775e-6*(291.15+120.0)/(globalTemp_+120.0)*pow(globalTemp_/291.15,(3.0/2.0)))/globalDensity_;
      if (this->inDensity_ != "const1234" && this->inTemp_ != "const1234"){
        kinVisco = (18.132941775e-6*(291.15+120.0)/(localTemp_[i]+120.0)*pow(localTemp_[i]/291.15,(3.0/2.0)))/localDensity_[i];
      }
      // turbulent length scale
      Double rmsOfVelFluct = sqrt(2.0/3.0*TKE_[i]);
      //turbLengthScale_[i] = fL*pow(rmsOfVelFluct,3.0)/TDR; //
      Double turbLengthScale = fL_*pow(c_mu_,0.75)*pow(TKE_[i],1.5)/TDR;
      // Kolmogorov wave number
      Double kolmogorovWN = pow(TDR,(1.0/4.0))*1.0/(pow(kinVisco,(3.0/4.0)));
      // essential locations in engery spectrum
      Double peakWN = 9.0*M_PI/55.0*A_/turbLengthScale;
      Double kmax = maxWN_ * peakWN;
      Double kmin = minWN_ * peakWN;

      // compute wave number increments
      Vector<Double> kn;
      kn.Resize(numModes_);
//      if(incrModes_ == "logarithmic"){
//        deltaWN_ = (log(kmax)-log(kmin))/(numModes_-1);
//        for(UInt p=0;p<numModes_;++p){
//          // in constrast to Bailly's method, one has to save the wave number sampling points
//          waveNumIncrements_[numModes_*k+p] = exp(log(kmin)+p*deltaWN_);
//          kn[p] = waveNumIncrements_[numModes_*k+p];
//          if(p==0) dWN_[0] = kmin;
//          else dWN_[p] = waveNumIncrements_[numModes_*k+p]-waveNumIncrements_[numModes_*k+p-1];
//        }
//      } else {
          Double deltaWN = (kmax-kmin)/(numModes_-1);
          for(UInt p=0;p<numModes_;++p){
            waveNumIncrements_[numModes_*k+p] = kmin + p*deltaWN;
            kn[p] = waveNumIncrements_[numModes_*k+p];
          }
//      }
      // mode loop
      for(UInt j=0; j<numModes_; ++j){
        // set engergy and velocity amplitude for j-th mode
        this->SetVelocityAmplitude(k,j,peakWN,deltaWN,kolmogorovWN,rmsOfVelFluct,kn);

         // set wave and direction vector for j-th mode
        this->SetRandVectors(k,j,rmsOfVelFluct,kn);
      }
  }
}

void SNGRFilter::InitTimelineBillsonSNGR(UInt i){
// initially for every time step the random values for one mode are generated new,
// this way for the initial velocity field all time steps are independent from one another.

  for(UInt k=0;k<numSteps_;++k){
  	initVelocity_.Resize(idsNodesToProcess_.GetSize()*3); //TODO, this field needs to be stored somewhere
    for(UInt p=0; p<idsNodesToProcess_.GetSize(); ++p){
      UInt i = idsNodesToProcess_[p];
      // how to get a coordinate
      CF::Vector<Double> pCoord;
      inGrid_->GetNodeCoordinate3D(pCoord,i+1); //TODO warum i+1
      // mode loop
    	for(UInt j=0;j<numModes_;++j){


      	// set wave and direction vector for j-th mode, redo this for every timestep in the initial velocity field, the initial field should be locally white noise
      	//this->SetRandVectorsBillson(j);

      	Double velA_coskx = 2*velAmplitude_[j]*cos(waveVec_[0]*pCoord[0] + waveVec_[1]*pCoord[1] + waveVec_[2]*pCoord[2] + phase_[j]);
      	initVelocity_[k*3] += velA_coskx*dirVec_[0];
      	initVelocity_[k*3+1] += velA_coskx*dirVec_[1];
      	initVelocity_[k*3+2] += velA_coskx*dirVec_[2];
    	}
  	}

	}
}

//TODO void SNGRFilter::FinalTimelineMethodBillsonSNGR(){
// solve convection equation

// apply time filter -> final velocity field
//}


// ===================================================================== //
//                         Method Block functions                        //
// ==================== Lafitte Sweeping-Based SNGR ==================== //
// ===================================================================== //


void SNGRFilter::UpdateResultLafitte(){
  UInt i = 0;
  Double sumTDR = 0.0;
  UInt numNodesToProcess = idsNodesToProcess_.GetSize();
  velAmplitude_.Resize(numModes_*numNodesToProcess);
  reconstTKE_.Resize(numNodesToProcess);
  peakWNLafitte_.Resize(numNodesToProcess);
  cutOffWNLafitte_.Resize(numNodesToProcess);
  numLargeScaleModes_.Resize(numNodesToProcess);
  timeScale_.Resize(numNodesToProcess);

  Vector<Double> kn;
  Double deltaWN = 0;
  // node loop TODO pragma omp parallel for
  for(UInt k=0; k<numNodesToProcess; ++k){
      i = idsNodesToProcess_[k];
      // dissipation rate
      Double TDR;
      if(inTEF_ == "Turbulence_Eddy_Fr9"){
        TDR = c_mu_*TKE_[i]*TEF_[i];
      }
      else{
        TDR = TEF_[i];
      }
      sumTDR += TDR;
      // turbulent length scale
      Double turbLengthScale = fL_*pow(c_mu_,0.75)*pow(TKE_[i],1.5)/TDR;
      // essential locations in engery spectrum
      peakWNLafitte_[i] = 9.0*M_PI/55.0*A_/turbLengthScale;
      cutOffWNLafitte_[i] = 1.8*peakWNLafitte_[i];
      Double kmax = maxWN_ * peakWNLafitte_[i];
      Double kmin = minWN_ * peakWNLafitte_[i];
      // turbulent time scale
      timeScale_[i] = ft_*TKE_[i]/TDR;

      // compute wave number increments
      kn.Resize(numModes_);
//      if(incrModes_ == "logarithmic"){
//        deltaWN_ = (log(kmax)-log(kmin))/(numModes_-1);
//        for(UInt p=0;p<numModes_;++p){
//          kn[p] = exp(log(kmin)+p*deltaWN_);
//          if(p==0) dWN_[0] = kmin;
//          else dWN_[p] = kn[p]-kn[p-1];
//        }
//      } else {
          deltaWN = (kmax-kmin)/(numModes_-1);
          for(UInt p=0;p<numModes_;++p){
            kn[p] = kmin + p*deltaWN;
            if ((kn[p]>cutOffWNLafitte_[k]) && (kn[p-1]<cutOffWNLafitte_[k])){ // && p!=0
              numLargeScaleModes_[k] = p-1; // if the user chooses minWN_ > ke (e.g. minWN_>=1.8*ke) this line will fail.
            } else if((p=0) && (kn[p]>cutOffWNLafitte_[k])){
                numLargeScaleModes_[k] = 0;
            }
//          }
      }
  // TODO include kinVisco, rms, kolmogorvWN here >> call to SNGRFilter::SetVelocityAmplitude also here.
  }

  // mean value of dissipation rate in most energetic grid points
  aveTDR_ = sumTDR/numNodesToProcess;
  // TODO Achtung !!!!!!!!
  // TODO Achtung !!!!!!!!
  // TODO Achtung !!!!!!!! -- we have to check, wheather there are local changes in the number of numLargeScaleModes_
  // TODO Achtung !!!!!!!!
  // TODO Achtung !!!!!!!!
  omega_.Resize(numLargeScaleModes_.GetSize()); // TODO omega_ fixed for all nodes but still dependent on mode; only set for numLargeScaleModes_, the rest does not get an omega_

  for(UInt j=0; j<numLargeScaleModes_.GetSize(); ++j){
    omega_[j] = pow(C_k_,0.5)*pow(aveTDR_,1.0/3.0)*pow(kn[j],2.0/3.0); // TODO kn is not saved over every node, but has to be in this setup.
  }



// TODO this loop is not done yet
  i = 0;
  for(UInt k=0; k<numNodesToProcess; ++k){
      i = idsNodesToProcess_[k];
      // kinematic viscosity nu=f(localTemp_) (Sutherland)
      Double kinVisco = (18.132941775e-6*(291.15+120.0)/(localTemp_[i]+120.0)*pow(localTemp_[i]/291.15,(3.0/2.0)))/localDensity_[i];
      // turbulent length scale
      Double rmsOfVelFluct = sqrt(2.0/3.0*TKE_[i]);
      // Kolmogorov wave number
      Double TDR;
      if(inTEF_ == "Turbulence_Eddy_Fr9"){
        TDR = c_mu_*TKE_[i]*TEF_[i];
      }
      else{
        TDR = TEF_[i];
      }
      Double kolmogorovWN = pow(TDR,(1.0/4.0))*1.0/(pow(kinVisco,(3.0/4.0)));

      // mode loop
      for(UInt j=0; j<numModes_; ++j){
        // set engergy and velocity amplitude for j-th mode
        this->SetVelocityAmplitude(k,j,peakWNLafitte_[i],deltaWN,kolmogorovWN,rmsOfVelFluct,kn);

         // set wave and direction vector for j-th mode
        this->SetRandVectors(k,j,rmsOfVelFluct,kn);
      }
  }
}

// TODO numModes_ has to be replaced with something usefull.
void SNGRFilter::InitArraysLafitte(){
  // initialize
  UInt numNodesToProcess = idsNodesToProcess_.GetSize();
//  turbLengthScale_(numModes_*nodesToProcess);
//  peakWN(numNodes_);

  phase_.Resize(numModes_*numNodesToProcess);
  waveVec_.Resize(numModes_*numNodesToProcess*3);
  dirVec_.Resize(numModes_*numNodesToProcess*3);
}

}
