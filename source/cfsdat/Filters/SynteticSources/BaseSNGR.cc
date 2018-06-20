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
 *       \author   r.krusche
 */
//================================================================================================


#include <Filters/SynteticSources/BaseSNGR.hh>
#include <random>
#include <chrono>

namespace CFSDat{

SNGRFilter::SNGRFilter(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
        :BaseFilter(numWorkers,config,resMan){
  this->filtStreamType_ = FIFO_FILTER;
  //INPUT from CFD
  this->inTKE_ = config->Get("TKE/resultName")->As<std::string>();
  this->inTEF_ = config->Get("TEF/resultName")->As<std::string>();
  this->inVelocity_ = config->Get("meanVelocity/resultName")->As<std::string>();
  this->inDensity_ = config->Get("localDensity/resultName")->As<std::string>();
  this->inTemp_ = config->Get("localTemp/resultName")->As<std::string>();
  //DECLARE Input
  upResNames.insert(inTKE_);
  upResNames.insert(inVelocity_);
  upResNames.insert(inDensity_);
  upResNames.insert(inTemp_);
  upResNames.insert(inTEF_);

  //OUTPUT
  this->outName_ = config->Get("output/resultName")->As<std::string>();
  this->interName_ = config->Get("intermediateResult/resultName")->As<std::string>();
  //DECLARE Output
  filtResNames.insert(interName_);
  filtResNames.insert(outName_);

  // using the TKE-criterion, the source region in the CFD-domain used for reconstructing turbulent velocity, is controlled
  // typical value in literature is 10%
  // only nodes, where a threshold of 10% of the max TKE in the whole domain is exceeded, are processed in SNGR-Algorithm
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
  this->maxFreq_ = config->Get("frequencyBounds/maxFreq")->As<UInt>();
  this->minFreq_ = config->Get("frequencyBounds/minFreq")->As<UInt>();

  // for maxWN and minWN percentage of the peak wave number are expected
  // typical values in literature: minWN = 0.1*peakWN; maxWN = 10*peakWN
  this->maxWN_ = config->Get("waveNumBounds/maxWN")->As<Double>();
  this->minWN_ = config->Get("waveNumBounds/minWN")->As<Double>();

  // enable ensemble-average at microphone location
  this->ensemble_ = config->Get("ensembles")->As<UInt>();
  
}

void SNGRFilter::FinishInit(){

  //input Grid
  inGrid_ = resultManager_->GetExtInfo(upResIds[0])->ptGrid;

  // chose method
  if(method_ == "BaillySNGR"){
    // go with Bailly, 1999
    std::cout << "Start reconstruction of synthetic turbulent velocity field using SNGR by Bailly & Juve, 1999." << std::endl;
    std::cout << "Acording to publication: A Stochastic Approach To Compute Subsonic Noise Using Linearized Euler's Equations." << std::endl;
    std::cout << "Doing " << ensemble_ << " indipendent reconstruction(s)." << std::endl;

    // Initialize SNGR method
    this->PrepareMethodBailly();
    this->UpdateResultMethodBailly();
    this->SetTimelineMethodBailly();
  }
  else if(method_ == "BillsonSNGR"){
    // go with Billson, 2003
    std::cout << "Start reconstruction of synthetic turbulent velocity field using SNGR by Billson, Eriksson & Davidson, 2003." << std::endl;
    std::cout << "Acording to publication: Jet Noise Prediction Using Stochastic Turbulence Modeling." << std::endl;
    std::cout << "Doing " << ensemble_ << " indipendent reconstruction(s)." << std::endl;

    // Initialize SNGR method
    this->PrepareMethodBailly();
    // Node loop
//    this->FinishB();

  }
//  else if(method_ == "next"){
//    // go with
//    std::cout << "Start reconstruction of synthetic turbulent velocity field using ." << std::endl;
//    std::cout << "Acording to publication: ." << std::endl;
//    std::cout << "Doing " << ensemble_ << " indipendent reconstruction(s)." << std::endl;
//
//  }
  else{
    EXCEPTION("This is not a valid method.")
  }
}


bool SNGRFilter::UpdateResults(std::set<uuids::uuid>& upResults){

// TODO matthias checkin
//  resultManager_->SetTimeValue(tkeId,aTF);
//  resultManager_->SetTimeValue(tefId,aTF);
//  resultManager_->SetTimeValue(velocityId,aTF);
//  resultManager_->SetTimeValue(densityId,aTF);
//  resultManager_->SetTimeValue(temperatureId,aTF);
  resultManager_->ActivateResult(tkeId_);
  resultManager_->ActivateResult(tefId_);
  resultManager_->ActivateResult(velocityId_);
  resultManager_->ActivateResult(densityId_);
  resultManager_->ActivateResult(temperatureId_);
  resultManager_->DeactivateResult(outId_);
  resultManager_->DeactivateResult(interId_);


  //now we call for upstream data in each source
  CF::StdVector< str1::shared_ptr<BaseFilter> >::iterator srcIter =  sources_.Begin();
  for(; srcIter != sources_.End() ; srcIter++){
    // should we check here anything for success?
    (*srcIter)->Run();
  }

  //equation numbers for this result
  //here, those equation numbers are equal for all in/out results
  Vector<Double>& returnVec = GetOwnResultVector<Double>(outId_);

#pragma omp parallel for
  for(UInt i=0;i<meanVelocity_.GetSize();++i){
    returnVec[i] = meanVelocity_[i];
    //if(i==idsNodesToProcess_[k]){   TODO stimmt so noch nicht. Ähnlich für returnVecInter
    //  returnVec[i] = turbReconstVelocity_[k];
    //  k++;
    //}
  }

  resultManager_->SetResultVecUpToDate(outId_,true);

  Vector<Double>& returnVecInter = GetOwnResultVector<Double>(interId_);

  //RETURN OTHER OUTPUT intermediate result
#pragma omp parallel for
  for(UInt i=0;i<TKE_.GetSize();++i){
    returnVecInter[i] = TKE_[i];
  }
  resultManager_->SetResultVecUpToDate(interId_,true);

  return true;
}


ResultIdList SNGRFilter::SetUpstreamResults(){
  //RESULT DEFINITION
  ResultIdList generated = SetDefaultUpstreamResults();
  tkeId_ = upResNameIds[inTKE_];
  velocityId_ = upResNameIds[inVelocity_];
  tefId_ = upResNameIds[inTEF_];
  densityId_ = upResNameIds[inDensity_];
  temperatureId_ = upResNameIds[inTemp_];

  return generated;
}

void SNGRFilter::AdaptFilterResults(){
  //We should check some validity...
  //we can almost copy everything from time input (also scalar, etc.)
  outId_ = filterResIds[0];
  resultManager_->CopyResultData(velocityId_,outId_);
  resultManager_->SetValid(outId_);

  //REGISTER INTERMEDIATE result
  interId_  = filterResIds[1];
  resultManager_->CopyResultData(tkeId_,interId_);
  resultManager_->SetValid(interId_);
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
  UInt k=0;
  for(UInt i=0; i<numNodes_; ++i){
    if(TKE_[i]>=minTKE_){
      idsNodesToProcess_[k] = i;
      k++;
    }
  }
}

void SNGRFilter::SetRandVectors(UInt k, UInt j){
  // Set up random number generators
  UInt seed_uniform = std::chrono::system_clock::now().time_since_epoch().count();
  std::mt19937 generator_uniform (seed_uniform);
  std::uniform_real_distribution<double> uniform01(0.0, 1.0);

  Double theta = 2 * M_PI * uniform01(generator_uniform);
  Double phi = acos(1 - 2 * uniform01(generator_uniform));
//    retVec[i*2] = 2 * M_PI * uniform01(generator);
//    retVec[i*2+1] = acos(1 - 2 * uniform01(generator));

  // random values for wave vector
  //Vector<Double>& randAngles_ = SNGRUtils::GetRandValues(numModes_);

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
  Double omegaExpectValue = rmsOfVelFluct_ * kn_[j];
  std::normal_distribution<double> distribution(omegaExpectValue,omegaExpectValue);
  omega_[numModes_*k+j] = distribution(generator_normal);

  //TODO revisit - set wave vector
  // magnitude of wave vector is wave number of current mode
  UInt g = 3*numModes_*k+3*j; // index in array for wave and direction vector
  waveVec_[g] = kn_[j]*cos(theta)*cos(phi);
  waveVec_[g+1] = kn_[j]*sin(phi);
  waveVec_[g+2] = -kn_[j]*sin(theta)*cos(phi);
  // set direction vector, perpendicular to wave vector, ensures incompressible velocity field
  dirVec_[g] = -sin(phi)*cos(alpha)*cos(theta)+sin(alpha)*sin(theta);
  dirVec_[g+1] = cos(phi)*cos(alpha);
  dirVec_[g+2] = sin(phi)*cos(alpha)*sin(theta)+sin(alpha)*cos(theta);

  // check if sufficiently perperndicular
  Double perp = waveVec_[g]*dirVec_[g]+waveVec_[g+1]*dirVec_[g+1]+waveVec_[g+2]*dirVec_[g+2];
  if(abs(perp) > 1.0e-9){
     perpFAIL_++;
//     std::cout << "WARNING! In node " << i <<" wave number vector not sufficiently perpendicular to direction vector!" << std::endl;
//     std::cout << "result of dot product was " << perp << std::endl;
  }
}

void SNGRFilter::SetVelocityAmplitude(UInt k, UInt j){
  // compute Karman-Pao spectrum for isotropic turbulence
  Double E_kn = A_*(pow(rmsOfVelFluct_,2.0)/peakWN_)*pow(kn_[j]/peakWN_,4.0)*exp(-2.0*pow(kn_[j]/kolmogorovWN_,2.0))*1.0/(pow((1.0+pow(kn_[j]/peakWN_,2.0)),(17.0/6.0)));
  // amplitude for current mode
  if(incrModes_ == "logarithmic") deltaWN_=dWN_[j];
  velAmplitude_[numModes_*k+j] = sqrt(E_kn*deltaWN_);
  reconstTKE_[k] += E_kn*deltaWN_;

//      // confirm that reconstructed velocity field is consistent with TKE from RANS analysis.
//      // literature gives reason to assume that this test fails often
//      Double diffTKE = abs(TKE_[i] - reconstTKE_[i]);
//      if(diffTKE > 1e-6){
//        std::cout << "WARNING! Deviation of reconstructed and input TKE in node " << i <<" too big!" << std::endl;
//        tkeFAIL_++;
//      }
}

// ===================================================================== //
//                         Method Block functions                        //
// ============================== SNGR ================================= //
// ===================================================================== //
void SNGRFilter::PrepareMethodBailly(){
  // compute time step, factor 30 for resolving a smooth curve
  deltaT_ = 1.0/(30*maxFreq_);
  // set length of output signal, factor 1.5 for safety
  if (sigLength_ == 0.0){
    sigLength_ = 1.5/(minFreq_);
  }

  // time values
  numSteps_ = ceil(sigLength_/deltaT_);
  Vector<Double> stepValues_(numSteps_);
  for (UInt idx=0; idx<numSteps_; idx++){
    stepValues_[idx] = idx*deltaT_;
  }

  StdVector<int> volRegionIds;
  inGrid_->GetVolRegionIds(volRegionIds);
  numNodes_ = inGrid_->GetNumNodes(volRegionIds);

  // get result for result id TODO evtl. Vector<Double>& als Datentyp vor den Namen wieder einfügen.
  CF::StdVector<UInt> numEqation;
  TKE_ = GetUpstreamResultVector<Double>(tkeId_,numEqation);
  meanVelocity_ = GetUpstreamResultVector<Double>(velocityId_,numEqation);
  localDensity_ = GetUpstreamResultVector<Double>(densityId_,numEqation);
  localTemp_ = GetUpstreamResultVector<Double>(temperatureId_,numEqation);
  TEF_ = GetUpstreamResultVector<Double>(tefId_,numEqation);

  // get tke threshold, actually gets the ids of all nodes that pass the tke criterion
  this->GetTkeThreshold();
  UInt numNodesToProcess = idsNodesToProcess_.GetSize();

  // initialize
//  Vector<Double> turbLengthScale_(numModes_*nodesToProcess);
//  Vector<Double> peakWN_(numNodes_);
  Vector<Double> velAmplitude_(numModes_*numNodesToProcess);
  Vector<Double> reconstTKE_(numNodesToProcess);
  Vector<Double> omega_(numModes_*numNodesToProcess);
  Vector<Double> phase_(numModes_*numNodesToProcess);
  Vector<Double> waveVec_(numModes_*numNodesToProcess*3);
  Vector<Double> dirVec_(numModes_*numNodesToProcess*3);
}

void SNGRFilter::UpdateResultMethodBailly(){
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
      Double kinVisco = (18.132941775e-6*(291.15+120.0)/(localTemp_[i]+120.0)*pow(localTemp_[i]/291.15,(3.0/2.0)))/localDensity_[i];
      // turbulent length scale
      rmsOfVelFluct_ = sqrt(2.0/3.0*TKE_[i]);
      //turbLengthScale_[i] = fL*pow(rmsOfVelFluct,3.0)/TDR; //
      Double turbLengthScale = fL_*pow(c_mu_,0.75)*pow(TKE_[i],1.5)/TDR;
      // Kolmogorov wave number
      kolmogorovWN_ = pow(TDR,(1.0/4.0))*1.0/(pow(kinVisco,(3.0/4.0)));
      // essential locations in engery spectrum
      Double peakWN_ = 9.0*M_PI/55.0*A_/turbLengthScale;
      Double kmax = maxWN_ * peakWN_;
      Double kmin = minWN_ * peakWN_;
      // compute wave number increments
      Vector<Double> kn_(numModes_);
      if(incrModes_ == "logarithmic"){
        deltaWN_ = (log(kmax)-log(kmin))/(numModes_-1);
        for(UInt p=0;p<numModes_;p++){
          kn_[p] = exp(log(kmin)+p*deltaWN_);
          if(p==0) dWN_[0] = kmin;
          else dWN_[p] = kn_[p]-kn_[p-1];
        }
      } else {
          deltaWN_ = (kmax-kmin)/(numModes_-1);
          for(UInt p=0;p<numModes_;p++){
            kn_[p] = kmin + p*deltaWN_;
          }
      }
      // mode loop
      for(UInt j=0; j<numModes_; ++j){
        // set engergy and velocity amplitude for j-th mode
        this->SetVelocityAmplitude(k,j);

         // set wave and direction vector for j-th mode
        this->SetRandVectors(k,j);
      }
  }
}

void SNGRFilter::SetTimelineMethodBailly(){
  // time loop
  for(UInt k=0; k<numSteps_; ++k){
    // initialize output
    Vector<Double> turbReconstVelocity_(idsNodesToProcess_.GetSize()*3);
    // node loop TODO pragma omp parallel for
    for(UInt p=0; p<idsNodesToProcess_.GetSize(); ++p){
      UInt i = idsNodesToProcess_[p];
      // how to get a coordinate
      CF::Vector<Double> pCoord;
      inGrid_->GetNodeCoordinate3D(pCoord,i+1); //TODO warum i+1
      // mode loop
      for(UInt j=0; j<numModes_; ++j){
        UInt g = numModes_*p*3+j*3;  // index in waveVec and dirVec arrays
        Double kx = waveVec_[g]*(pCoord[0]-meanVelocity_[i*3]*stepValues_[k]) + waveVec_[g+1]*(pCoord[1]-meanVelocity_[i*3+1]*stepValues_[k]) + waveVec_[g+2]*(pCoord[2]-meanVelocity_[i*3+2]*stepValues_[k]);
        Double sumElem = 2.0*velAmplitude_[numModes_*p+j]*cos(kx + phase_[numModes_*p+j] + fa_*omega_[numModes_*p+j]*stepValues_[k]);
        turbReconstVelocity_[p*3] += sumElem*dirVec_[g];
        turbReconstVelocity_[p*3+1] += sumElem*dirVec_[g+1];
        turbReconstVelocity_[p*3+2] += sumElem*dirVec_[g+2];
      }
    }

    // write node results to file
    //TODO muss hier der Inhalt aus der Fkt UpdateResults angefgügt werden?
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
  }
  std::cout << "Reconstruction of turbulent velocity field finished" << std::endl;
}

// ===================================================================== //
//                         Method Block functions                        //
// ============================ Billson SNGR =========================== //
// ===================================================================== //
/*
void SNGRFilter::InitTimelineBillsonSNGR(UInt i){
// initially for every time step the random values for one mode are generated new,
// this way for the initial velocity field all time steps are independent from one another.

  // re-initialise intermidiate velocity field for every node
  Vector<Double> initVelocity_(numSteps_*3); //TODO wo soll dieses initial field zwischen gespeichert werden?

  for(UInt k=0;k<numSteps_;k++){
    for(UInt j=0;j<numModes_;j++){
      // set engergy and velocity amplitude for j-th mode
      this->SetVelocityAmplitude(i,j);

      // set wave and direction vector for j-th mode
      this->SetRandVectors(j);

      Double velA_coskx = 2*velAmplitude_[j]*cos(waveVec_[0]*pCoord_[0] + waveVec_[1]*pCoord_[1] + waveVec_[2]*pCoord_[2] + phase_[j]);
      initVelocity_[k*3] += velA_coskx*dirVec_[0];
      initVelocity_[k*3+1] += velA_coskx*dirVec_[1];
      initVelocity_[k*3+2] += velA_coskx*dirVec_[2];
    }
  }
  // call to final time line
  //this->FinalTimelineMethodBillsonSNGR
}

//TODO void SNGRFilter::FinalTimelineMethodBillsonSNGR(){
// solve convection equation

// apply time filter -> final velocity field
//}
*/

}
