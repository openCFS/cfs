// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     sngrBase.cc
 *       \brief    <Description>
 *
 *       \date     Nov 01, 2017
 *       \author   r.krusche
 */
//================================================================================================


#include <Filters/SynteticSources/BaseSNGR.hh>
#include <random>

namespace CFSDat{

SNGRFilter::SNGRFilter(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
        :BaseFilter(numWorkers,config,resMan){
  this->filtStreamType_ = FIFO_FILTER;
  //INPUT
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


  // mit dem TKE-Kriterium wird das Quellgebiet reduziert. Gänger Wert in der Literatur sind zB 10%; 
  // es werden dann der SNGR-Algorithmus nur für Knoten ausgeführt, die mindestens 10% der maximal auftretenden TKE erreichen.
  this->TKEcrit_ = config->Get("tkeCriterion")->As<Double>();

  // falls der Nutzer eine bestimmte Länge des Signals braucht kann er diese vorgeben
  this->sigLength_ = config->Get("signalLength")->As<Double>();
  this->fL_ = config->Get("lengthScaleFactor")->As<Double>();
  this->ft_ = config->Get("timeScaleFactor")->As<Double>();
  this->fa_ = config->Get("angularFreqFactor")->As<Double>();
  this->numModes_ = config->Get("numOfModes")->As<UInt>();
  this->incrModes_ = config->Get("incrementModes")->As<std::string>();

  // Frequenzen, die der Nutzer akustisch auflösen möchte. 
  // theoretisch müssen man diese in Abstimmung mit der nächsten Abfrage wählen
  this->maxFreq_ = config->Get("frequencyBounds/maxFreq")->As<UInt>();
  this->minFreq_ = config->Get("frequencyBounds/minFreq")->As<UInt>();

  // mit maxWN u. minWN werden Prozent von der Peak-Wellenzahl abgefragt. 
  // in Literatur zB üblich: minWN = 0.1*peakWN; maxWN = 10*peakWN
  this->maxWN_ = config->Get("waveNumBounds/maxWN")->As<Double>();
  this->minWN_ = config->Get("waveNumBounds/minWN")->As<Double>();

  // mit der Angabe von ensembles soll im Nachgang eine Ensembles-Mittelung an einem Mikrofonpunkt möglich sein.
  this->ensemble_ = config->Get("ensembles")->As<UInt>();
  
}

void SNGRFilter::FinishInit(){

  //input Grid
  inGrid_   = resultManager_->GetExtInfo(upResIds[0])->ptGrid;

  //Initialize SNGR method
  this->PrepareMethodBailly();

}


bool SNGRFilter::UpdateResults(std::set<uuids::uuid>& upResults){
//   --------------------------------------
//  ============== SNGR CORE =================
//  ----------------------------------------
  
  std::cout << "Start reconstruction of synthetic turbulent velocity field using SNGR by Bailly & Juve." << std::endl;
  std::cout << "Doing " << ensemble_ << " indipendent reconstruction(s)." << std::endl;
  
  // compute time step, factor 30 for resolving a smooth curve
  deltaT_ = 1.0/(30*maxFreq_);
  // set length of output signal, factor 2 for safety
  if (sigLength_ == 0.0){
    sigLength_ = 2.0/(minFreq_);
  }
  UInt numSteps = ceil(sigLength_/deltaT_);
  // time values
  Vector<Double> stepValues(numSteps);
  for (UInt idx=0; idx<numSteps; idx++){
    stepValues[idx] = idx*deltaT_;
  }

  StdVector<int> volRegionIds;
  inGrid_->GetVolRegionIds(volRegionIds);
  UInt numNodes = inGrid_->GetNumNodes(volRegionIds);

  // initialize
  UInt flg=0; // counter for number of nodes the TKE-Criterion is met
  Vector<Double> turbLengthScale(numNodes);
  Vector<Double> peakWN(numNodes);
  Vector<Double> velAmplitude(numNodes*numModes_);
  Vector<Double> reconstTKE(numNodes); 
  UInt tkeFAIL = 0;     // counter for nodes for which the deviation of reconstructed and read TKE is not small enough.
  UInt perpFAIL = 0;    // counter for nodes for which waveVec and dirVec are not sufficiently perpendicular to one another
  // constants and parameters 
  Double A = 1.452762;
  //Double fL = 2.5;
  //Double ft_ = 0.0002;
  Double c_mu = 0.09; // constand of k-epsilon model

  // get result for result id 
  CF::StdVector<UInt> numEqation;
  Vector<Double>& TKE = GetUpstreamResultVector<Double>(tkeId_,numEqation);

  //TODO von TDR nach TEF umrechnen
  Vector<Double>& TEF = GetUpstreamResultVector<Double>(tefId_,numEqation);


  Vector<Double>& meanVelocity = GetUpstreamResultVector<Double>(velocityId_,numEqation);
  Vector<Double>& localDensity = GetUpstreamResultVector<Double>(densityId_,numEqation);
  Vector<Double>& localTemp = GetUpstreamResultVector<Double>(temperatureId_,numEqation);


  // set TKE threshold
  Double maxTKE=0.0;
  for(UInt m=0; m<numNodes; m++){
    if(TKE[m]>=maxTKE) maxTKE=TKE[m];
  }
  Double minTKE = TKEcrit_*maxTKE;

  // node loop
  for(UInt i=0; i<numNodes; i++){

    if(TKE[i]>=minTKE){
      flg++;
      // initialize output
      Vector<Double> turbReconstVelocity_(numSteps*3);
      reconstTKE[i] = 0.0;
      Vector<Double> waveVec(3);
      Vector<Double> dirVec(3);

      // dissipation rate and viscosity (Sutherland)
      Double epsilon = c_mu*TKE[i]*TEF[i];
      Double kinVisco = (18.132941775e-6*(291.15+120.0)/(localTemp[i]+120.0)*pow(localTemp[i]/291.15,(3.0/2.0)))/localDensity[i];
      // turbulent length scale
      Double rmsOfVelFluct = sqrt(2.0/3.0*TKE[i]);
      //turbLengthScale[i] = fL*pow(rmsOfVelFluct,3.0)/epsilon; //
      turbLengthScale[i] = fL_*pow(c_mu,0.75)*pow(TKE[i],1.5)/epsilon;
      // Kolmogorov wave number
      Double kolmogorovWN = pow(epsilon,(1.0/4.0))*1.0/(pow(kinVisco,(3.0/4.0)));
      // essential locations in engery spectrum
      peakWN[i] = 9.0*M_PI/55.0*A/turbLengthScale[i];
      Double kmax = maxWN_ * peakWN[i];
      Double kmin = minWN_ * peakWN[i];
      // compute wave number increments
      Double deltaWN;
      Vector<Double> kn(numModes_);
      if(incrModes_ == "logarithmic"){
        deltaWN = (log(kmax)-log(kmin))/(numModes_-1);
        for(UInt p=0;p<numModes_;p++){
          kn[p] = exp(log(kmin)+p*deltaWN);
        }
      } else {
          deltaWN = (kmax-kmin)/(numModes_-1);
          for(UInt p=0;p<numModes_;p++){
            kn[p] = kmin + p*deltaWN;
          }
      }

      // how to get a coordinate
      CF::Vector<Double> pCoord;
      inGrid_->GetNodeCoordinate3D(pCoord,i+1);

      // loop over modes
      for(UInt j=0; j<numModes_; j++){
        // compute Karman-Pao spectrum for isotropic turbulence
        Double E_kn = A*(pow(rmsOfVelFluct,2.0)/peakWN[i])*pow(kn[j]/peakWN[i],4.0)*exp(-2.0*pow(kn[j]/kolmogorovWN,2.0))*1.0/(pow((1.0+pow(kn[j]/peakWN[i],2.0)),(17.0/6.0)));
        // amplitude for current mode
        velAmplitude[i*j] = sqrt(E_kn*deltaWN); // achtung hier muss delta WN entlogarithmiert werden falls incrModes==logarithmic
        reconstTKE[i] += E_kn*deltaWN;
        // select random parameters for wave vector, phase and direction of mode.
        srand(time(NULL));      // initialize random seed.
        Double theta = (M_PI - 0.0)*((Double)rand()/(Double)RAND_MAX); // + 0.0;        // polar angle of wave vector, random value between 0 ... pi
        srand(time(NULL)+33.0);
        Double phi = (2.0*M_PI - 0.0)*((Double)rand()/(Double)RAND_MAX); // + 0.0;        // azimuthal angle of wave vector, random value between 0 ... 2pi
        srand (time(NULL)+149.0);
        Double alpha = (2*M_PI - 0.0)*((Double)rand()/(Double)RAND_MAX); // + 0.0;      // angle for directional vector of current mode, random value between 0 ... 2pi
        srand (time(NULL)+273.0);
        Double phase = (2.0*M_PI - 0.0)*((Double)rand()/(Double)RAND_MAX); // + 0.0;      // phase of current mode, random value between 0 ... 2pi


        // set wave vector
        // einfache Kugelkoordinaten, der Betrag des Vektors muss der Wellenzahl der aktuellen Mode entsprechen.
        waveVec[0] = kn[j]*cos(theta)*cos(phi);
        waveVec[1] = kn[j]*sin(phi);
        waveVec[2] = -kn[j]*sin(theta)*cos(phi);
        // set direction vector, perpendicular to wave vector, ensures incompressible velocity field
        // die Einträge des Richtungsvektors habe ich händisch so berechnet, dass das Skalarprodukt beider 0 sein müsste
        // vlt. gibt es aber noch eine elegantere Lösung
        dirVec[0] = -sin(phi)*cos(alpha)*cos(theta)+sin(alpha)*sin(theta);
        dirVec[1] = cos(phi)*cos(alpha);
        dirVec[2] = sin(phi)*cos(alpha)*sin(theta)+sin(alpha)*cos(theta);

        // set angular frequency for j-th mode
        //TODO Does it match here WICHTIGH ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //TODO which Does it match here WICHTIGH ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        //TODO which Does it match here WICHTIGH ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        Double omegaExpectValue = rmsOfVelFluct * kn[j];


        // check if sufficiently perperndicular
        Double perp = waveVec[0]*dirVec[0]+waveVec[1]*dirVec[1]+waveVec[2]*dirVec[2];
        if(abs(perp) > 1.0e-9){
          perpFAIL++;
          std::cout << "WARNING! In node " << i <<" wave number vector not sufficiently perpendicular to direction vector!" << std::endl;
          std::cout << "result of dot product was " << perp << std::endl;
        }

        // compute time signal for current mode, loop over time step values
        for(UInt k=0; k<numSteps; k++){
          Double kx = waveVec[0]*(pCoord[0]-meanVelocity[i*3]*stepValues[k]) + waveVec[1]*(pCoord[1]-meanVelocity[i*3+1]*stepValues[k]) + waveVec[2]*(pCoord[2]-meanVelocity[i*3+2]*stepValues[k]);
          //TODO which omega of the mode? WICHTIGH ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          //TODO which omega WICHTIGH ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          //TODO which omega WICHTIGH ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          //TODO which omega WICHTIGH ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          Double omega = omegaExpectValue;
          std::normal_distribution<double> distribution(5.0,2.0);
          //TODO which omega WICHTIGH ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          //TODO which omega WICHTIGH ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          //TODO which omega WICHTIGH ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
          //TODO which omega WICHTIGH ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

          turbReconstVelocity_[k*3] += 2*velAmplitude[i*j]*cos(kx + phase + fa_*omega*stepValues[k])*dirVec[0]; //TODO 2* ??
          turbReconstVelocity_[k*3+1] += 2*velAmplitude[i*j]*cos(kx + phase + fa_*omega*stepValues[k])*dirVec[1];
          turbReconstVelocity_[k*3+2] += 2*velAmplitude[i*j]*cos(kx + phase + fa_*omega*stepValues[k])*dirVec[2];
        }

      }

//      // confirm that reconstructed velocity field is consistent with TKE from RANS analysis.
//      // Literatur lässt vermuten, dass das dieser Test sehr häufig negativ ausfällt.
//      Double diffTKE = abs(TKE[i] - reconstTKE[i]);
//      if(diffTKE > 1e-6){
//        std::cout << "WARNING! Deviation of reconstructed and input TKE in node " << i <<" too big!" << std::endl;
//        tkeFAIL++;
//      }

//      // final turbulent velocity component containing all modes.
//      for(UInt p=0; p<numSteps; p++){
//        turbReconstVelocity_[p*3] *= 2;
//        turbReconstVelocity_[p*3+1] *= 2;
//        turbReconstVelocity_[p*3+2] *= 2;
//      }
    } // end of tkeCriterion
  } // end of node loop

  std::cout << "Writing reconstructed velocity field to file." << std::endl;

  std::cout << "The reconstruction was carried out for " << flg << " of " << numNodes << " nodes, " << flg/numNodes*100.0 << "%." << std::endl;
  std::cout << "The remaining " << numNodes-flg << "(" << (numNodes-flg)/numNodes*100.0 << "%) nodes did not meet the TKE-Criterion." << std::endl;
  if(tkeFAIL!=0){
    std::cout << "There has been a problem with the reconstructed TKE in " << tkeFAIL << " nodes." << std::endl;
  }
  if(perpFAIL!=0){
    std::cout << "For " << perpFAIL << " nodes the wave vector and direction vector where not sufficiently perpendicular." << std::endl;
  }

  std::cout << "Reconstruction of turbulent velocity field finished" << std::endl;

//    --------------------------------------
//  ============== SNGR =================
//  ----------------------------------------

  // now we calculate all results, independently if requested or not for one time step, 
  // so set equal time steps for all results
  UInt stepIndex = resultManager_->GetStepIndex(*upResults.begin());
  for (UInt i = 0; i < filterResIds.GetSize(); i++) {
    resultManager_->SetStepIndex(filterResIds[i], stepIndex);
  }

  //equation numbers for this result
  //here, those equation numbers are equal for all in/out results
  if (calcOutId_) {
    Vector<Double>& returnVec = GetOwnResultVector<Double>(outId_);
#pragma omp parallel for
    for(UInt i=0;i<meanVelocity.GetSize();++i){
      returnVec[i] = meanVelocity[i];
    }
    resultManager_->SetResultVecUpToDate(outId_,true);
  }
  
  if (calcInterId_) {
    Vector<Double>& returnVecInter = GetOwnResultVector<Double>(interId_);
    //RETURN OTHER OUTPUT intermediate result
#pragma omp parallel for
    for(UInt i=0;i<TKE.GetSize();++i){
      returnVecInter[i] = TKE[i];
    }
    resultManager_->SetResultVecUpToDate(interId_,true);
  }
  
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
  
  if (filterResNameIds.find(outName_) != filterResNameIds.end()) {
    outId_ = filterResNameIds[outName_];
    resultManager_->CopyResultData(velocityId_,outId_);
    resultManager_->SetValid(outId_);
    calcOutId_ = true;
  } else {
    calcOutId_ = false;
  }
  
  if (filterResNameIds.find(interName_) != filterResNameIds.end()) {
    interId_ = filterResNameIds[interName_];
    resultManager_->CopyResultData(tkeId_,interId_);
    resultManager_->SetValid(interId_);
    calcInterId_ = true;
  } else {
    calcInterId_ = false;
  }
}

void SNGRFilter::PrepareMethodBailly(){

}

void SNGRFilter::UpdateResultMethodBailly(){

}


}
