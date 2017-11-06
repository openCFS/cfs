// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     Lighthill.cc
 *       \brief    <Description>
 *
 *       \date     Oct 6, 2016
 *       \author   kroppert
 */
//================================================================================================


#include "Lighthill.hh"
//#include "FeBasis/H1/H1Elems.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include <algorithm>
#include <vector>

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/basic.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#include <boost/iterator/zip_iterator.hpp>
#include <utility>


namespace CFSDat{

Lighthill::Lighthill(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
:AeroacousticBase(numWorkers,config,resMan){

#ifndef USE_CGAL
  EXCEPTION("CoefFunctionScatteredData needs to be compiled with USE_CGAL=ON!");
#endif

  this->filtStreamType_ = FIFO_FILTER;


  std::cout<<("============================================================ \n"
      "Make sure you set the ''density='' in the xml-scheme correctly, \n"
      "otherwise it is chosen to be 1.0 !! \n"
      "============================================================")<<std::endl;

  density_ = config->Get("scheme")->Get("density")->As<Double>();

  // ( LambVectorLighthill || LighthillSourceVector || LighthillSourceTerm )
  Form_ = config->Get("type")->As<std::string>();

  numNeighbors_ = (trgGrid_->GetDim() == 2) ? 4 : 8;

  p_ = config->Get("scheme")->Get("interpolationExponent")->As<UInt>();

  if(config->Get("ResultList")->Get("vorticity")->Has("resultName") == false){
    std::cout<<"Vorticity is computed by CFSDat...Curl(u)"<<std::endl;
    externVorticity_ = false;
  }else{
    externVorticity_ = true;
  }

  // add output to result manager
  std::string outRes = params_->Get("ResultList")->Get("outputQuantity")->Get("resultName")->As<std::string>();
  filtResNames.insert(outRes);

  //add velocity input result to manager
  this->res1Name = params_->Get("ResultList")->Get("velocity")->Get("resultName")->As<std::string>();
  upResNames.insert(res1Name);

  // if an external vorticity input exists, add it to manager
  if (externVorticity_ == true){
    this->res2Name = params_->Get("ResultList")->Get("vorticity")->Get("resultName")->As<std::string>();
    upResNames.insert(res2Name);
  }

  checkSum_ = false;
  if(config->Has("sourceSum")){
    checkSum_ = config->Get("sourceSum")->As<bool>();
  }
}


Lighthill::~Lighthill(){

}


bool Lighthill::Run(){
  // we deactivate every result, except for our own
  std::set<uuids::uuid> activeResults = resultManager_->GetActiveResults();
  std::set<uuids::uuid>::iterator aIter = activeResults.begin();


  for(; aIter != activeResults.end(); ++aIter){
    if(filterResIds.Find(*aIter) == -1){
      WARN(" There are still active results when reaching the derivative filter. This indicates an unexpected use of the pipeline.")
    }
    resultManager_->DeactivateResult(*aIter);
  }

  Double aTF = resultManager_->GetStepValue(filterResIds[0]);


  resultManager_->SetTimeValue(upResIds[0],aTF);
  resultManager_->ActivateResult(upResIds[0]);

  // Activate extern vorticity-input, if provided
  if(externVorticity_ == true){
    resultManager_->SetTimeValue(upResIds[1],aTF);
    resultManager_->ActivateResult(upResIds[1]);
  }

  //now we call for upstream data in each source
  CF::StdVector< str1::shared_ptr<BaseFilter> >::iterator srcIter =  sources_.Begin();
  for(; srcIter != sources_.End() ; srcIter++){
    // should we check here anything for success?
    (*srcIter)->Run();
  }


  CF::StdVector<UInt> eqnNums;
  /// this is the vector, which will be filled with the wanted quantity
  Vector<Double>& returnVec = resultManager_->GetResultVector<Double>(filterResIds[0],eqnNums);
  returnVec.Init();


  Vector<Double> tempRetVec;
  if(Form_ == "AeroacousticSource_LambVector"){LambVector(tempRetVec);}
  if(Form_ == "AeroacousticSource_LighthillSourceVector"){LighthillSourceVector(tempRetVec);}
  if(Form_ == "AeroacousticSource_LighthillSourceTerm"){LighthillSourceTerm(tempRetVec);}

  returnVec = tempRetVec * density_;

  // Check filter mesh and output values
  if(checkSum_ == 1 && Form_ == "AeroacousticSource_LighthillSourceTerm"){
    // Check output values
    // d'Alembertoperator(p)= d(Source_i)/dx_i
    // (a) INTEGRATION(Source_i) over the domain must be approximately zero
    Double intSource = returnVec.Sum();
    std::cout<<"Sum over all sources (not integrated) = "<<intSource<<std::endl;
    // (b) INTEGRATION(d(Source_i)/dx_i) over the domain must be approximately zero
    //                                   if Source_i(Wall) = 0 and Source_i(In-/Outlet)=0

    // Check mesh
    // Mesh size must be smaller than << Ma*c/(20*f) in the source region
    // to determine acoustically relevant flow structures
  }

  resultManager_->ActivateResult(filterResIds[0]);

  //now deactivate own upstream results
  for(UInt aRes=0;aRes<upResIds.GetSize();aRes++){
    resultManager_->DeactivateResult(upResIds[aRes]);
  }

  return true;
}


void Lighthill::LambVector(Vector<Double>& tempRetVec){
  CF::StdVector<UInt> eqnNums;
  Vector<Double>& inVecVel = resultManager_->GetResultVector<Double>(upResIds[0],eqnNums);

  Matrix& matrix = matrices_[matrixIndex_];
  const UInt maxNumTrgEntities = matrix.numTargets;

  Vector<Double> LambVec;
  if(externVorticity_ == true && Form_ == "AeroacousticSource_LambVector"){
    /************* case of externally provided vorticity ***************/
    // result is defined on nodes !!!
    Vector<Double>& inVecOmega = resultManager_->GetResultVector<Double>(upResIds[1],eqnNums);
    OmegaVectorProductU(inVecVel, inVecOmega, LambVec, numEquPerEnt_);
  }else{
    /************* we compute the vorticity internally ***************/
    // result is defined on elements !!!
    StdVector<CF::UInt>& targetSourceIndexNtE = matrix.targetSourceIndexNtE;
    StdVector<CF::UInt>& targetSourceNtE = matrix.targetSourceNtE;
    StdVector< CF::Matrix<CF::Double> >& targetSourceFactorCurl = matrix.targetSourceFactorCurl;
    UInt gridDim = Grid_->GetDim();
    Vector<Double> Omegatemp;
    Vector<Double> Omega;
    // Curl(u)
    CalcCurl(Omegatemp, inVecVel, numEquPerEnt_, targetSourceNtE, targetSourceIndexNtE, numNeighbors_, targetSourceFactorCurl, maxNumTrgEntities, gridDim);
    if(numEquPerEnt_ == 2){
      Omega = TwoDToScalar(Omegatemp);
    }else{
      Omega = Omegatemp;
    }

    // interpolate input-velocity to cell-centers, because the calculated curl is now defined on elements
    Vector<Double> velE;
    Node2Cell(velE, inVecVel, numEquPerEnt_, targetSourceNtE, targetSourceIndexNtE, numNeighbors_, maxNumTrgEntities, gridDim);

    OmegaVectorProductU(velE, Omega, LambVec, numEquPerEnt_);
  }
  tempRetVec = LambVec;
}




void Lighthill::LighthillSourceVector(Vector<Double>& tempRetVec){
  CF::StdVector<UInt> eqnNums;
  Vector<Double>& inVecVel = resultManager_->GetResultVector<Double>(upResIds[0],eqnNums);


  Matrix& matrix = matrices_[matrixIndex_];
  StdVector<CF::UInt>& targetSourceIndexNtE = matrix.targetSourceIndexNtE;
  StdVector<CF::UInt>& targetSourceNtE = matrix.targetSourceNtE;
  StdVector< CF::Matrix<CF::Double> >& targetSourceFactorGrad = matrix.targetSourceFactorGrad;
  UInt gridDim = Grid_->GetDim();
  Vector<Double> Omegatemp;
  Vector<Double> Omega;
  const UInt maxNumTrgEntities = matrix.numTargets;

  // interpolate input-velocity to cell-centers
  Vector<Double> velE;
  Node2Cell(velE, inVecVel, numEquPerEnt_, targetSourceNtE, targetSourceIndexNtE, numNeighbors_, maxNumTrgEntities, gridDim);

  /**************** Term 1: LambVector *****************************/
  Vector<Double> LambVec;
  this->LambVector(LambVec);


  /**************** Term 2: Grad(1/2 u*u) *****************************/
  Vector<Double> uuN;
  ScalarProduct(uuN, inVecVel, inVecVel, numEquPerEnt_, 0.5);

  Vector<Double> GraduuN;
  CalcGradient(GraduuN, uuN, 1, targetSourceNtE, targetSourceIndexNtE, numNeighbors_, targetSourceFactorGrad, maxNumTrgEntities, numEquPerEnt_);

  /**************** Term 1 + Term 2 *****************************/
  tempRetVec = GraduuN + LambVec;
}



void Lighthill::LighthillSourceTerm(Vector<Double>& tempRetVec){
  CF::StdVector<UInt> eqnNums;

  Matrix& matrix = matrices_[matrixIndex_];
  StdVector<CF::UInt>& targetSourceIndexNtE = matrix.targetSourceIndexNtE;
  StdVector<CF::UInt>& targetSourceIndexEtN = matrix.targetSourceIndexEtN;
  StdVector<CF::UInt>& targetSourceNtE = matrix.targetSourceNtE;
  StdVector<CF::UInt>& targetSourceEtN = matrix.targetSourceEtN;
  StdVector<CF::Double>& targetSourceNNFactor = matrix.targetSourceNNFactor;
  StdVector< CF::Matrix<CF::Double> >& targetSourceFactorDiv = matrix.targetSourceFactorDiv;
  Vector<Double> Omegatemp;
  Vector<Double> Omega;
  const UInt maxNumTrgEntities = matrix.numTargets;
  const UInt maxNumSrcEntities = matrix.numSources;


  Vector<Double> LightVec;
  this->LighthillSourceVector(LightVec);

  // interpolate from element-centroid to nodes
  Vector<Double> LightVecN;
  NearestNeighbourInterpolation(LightVecN, LightVec, numEquPerEnt_, targetSourceEtN, targetSourceIndexEtN, numNeighbors_, targetSourceNNFactor, maxNumSrcEntities);

  /**************** Div(Term 1 + Term 2) *****************************/
  CalcDivergence(tempRetVec, LightVecN, numEquPerEnt_, targetSourceNtE, targetSourceIndexNtE, numNeighbors_, targetSourceFactorDiv, maxNumTrgEntities);
}



// cgal copy and paste stuff for association of points to CF::UInt from
// http://doc.cgal.org/Manual/4.2/doc_html/cgal_manual/Spatial_searching/Chapter_main.html#Subsection_61.3.1
//
// later on this should be refactored into KNNSearch.hh and KNNSearch.cc
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef boost::tuple<Point_3,CF::UInt> Point_and_int;

//definition of the property map
struct My_point_property_map{
  typedef Point_3 value_type;
  typedef const value_type& reference;
  typedef const Point_and_int& key_type;
  typedef boost::readable_property_map_tag category;
};

//get function for the property map
inline My_point_property_map::reference
get(My_point_property_map,My_point_property_map::key_type p)
{return boost::get<0>(p);}

typedef CGAL::Random_points_in_cube_3<Point_3>                                          Random_points_iterator;
typedef CGAL::Search_traits_3<Kernel>                                                   Traits_base;
typedef CGAL::Search_traits_adapter<Point_and_int,My_point_property_map,Traits_base>    Traits;


typedef CGAL::Orthogonal_k_neighbor_search<Traits>                      K_neighbor_search;
typedef K_neighbor_search::Tree                                         Tree;
typedef K_neighbor_search::Distance                                     Distance;


CF::StdVector<Lighthill*> Lighthill::differentiators_;
CF::StdVector<Lighthill::Matrix> Lighthill::matrices_;




void Lighthill::PrepareCalculation(){
  std::cout << "\t ---> Lighthill preparing for calculation" << std::endl;
  //if(externVorticity_ == false || Form_ != "AeroacousticSource_LambVector"){
  std::cout << "\t\t 1/4 Obtaining source entities (nodes) " << std::endl;
  uuids::uuid upRes = upResIds[0];
  Grid_ = resultManager_->GetExtInfo(upRes)->ptGrid;

  scrMap_ = resultManager_->GetResultAdapter(upRes)->mapping;
  trgMap_ = resultManager_->GetResultAdapter(filterResIds[0])->mapping;
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  numEquPerEnt_ = scrMap_->GetNumEqnPerEnt();
  if(numEquPerEnt_ != trgGrid_->GetDim()) EXCEPTION("Grid and Values don't have the same dimension!!");
  if(numEquPerEnt_ == 1) EXCEPTION("Input (velocity) is a scalar, it should be a vector!!");


  // input MUST be defined on nodes !!
  //bool inElems = inInfo->definedOn == ExtendedResultInfo::ELEMENT;

  const CF::UInt maxNumSrcNodeEntities = scrMap_->GetNumEntities(); //velocity...upResIds[0]
  StdVector<CF::UInt> globSrcNodeEntity;
  GetUsedMappedEntities(scrMap_, globSrcNodeEntity, srcRegions_, Grid_);
  const CF::UInt numSrcNodeEntities = CountUsedEntities(globSrcNodeEntity);
  if (numSrcNodeEntities < numNeighbors_) {
    numNeighbors_ = numSrcNodeEntities;
  }

  differentiators_.Push_back(this);
  matrixIndex_ = matrices_.GetSize();
  matrices_.Resize(matrixIndex_ + 1);
  Matrix& matrix = matrices_[matrixIndex_];


  std::cout << "\t\t 2/4 Obtaining target entities (elements)" << std::endl;

  const CF::UInt maxNumTrgEntities = trgMap_->GetNumEntities();
  StdVector<CF::UInt> globTrgEntity;
  StdVector<CF::UInt> globTrgNodeEntity;
  GetUsedMappedEntities(trgMap_, globTrgEntity, trgRegions_, trgGrid_);
  const CF::UInt numTrgEntities = CountUsedEntities(globTrgEntity);

  std::cout << "\t\t\t Differentiator is dealing with " << numSrcNodeEntities <<
      " source " << ("nodes") << " and "<< numTrgEntities << " target elements" << std::endl;

  std::cout << "\t\t 3/4 Creating search tree " << std::endl;
  std::vector<Point_3> pointsNtE; //NtE...nodes to elements
  std::vector<CF::UInt> indicesNtE;
  StdVector< CF::Vector<Double> > sCoordNtE;
  sCoordNtE.Resize(maxNumSrcNodeEntities);
  for(CF::UInt srcEnt = 0; srcEnt < maxNumSrcNodeEntities; srcEnt++) {
    CF::UInt globEntityNumber = globSrcNodeEntity[srcEnt];
    if (globEntityNumber != UnusedEntityNumber) {
      CF::Vector<Double> pCoord;
      Grid_->GetNodeCoordinate3D(pCoord, globEntityNumber);
      Point_3 point(pCoord[0],pCoord[1],pCoord[2]);
      pointsNtE.push_back(point);
      indicesNtE.push_back(srcEnt);
      sCoordNtE[srcEnt] = pCoord;
    }
  }

  std::vector<Point_3> pointsEtN; //EtN...elements to nodes
  std::vector<CF::UInt> indicesEtN;
  StdVector< CF::Vector<Double> > sCoordEtN;
  sCoordEtN.Resize(maxNumTrgEntities);
  for(CF::UInt srcEnt = 0; srcEnt < maxNumTrgEntities; srcEnt++) {
    CF::UInt globEntityNumber = globTrgEntity[srcEnt];
    if (globEntityNumber != UnusedEntityNumber) {
      CF::Vector<Double> pCoord;
      Grid_->GetElemCentroid(pCoord, globEntityNumber, true);
      Point_3 point(pCoord[0],pCoord[1],pCoord[2]);
      pointsEtN.push_back(point);
      indicesEtN.push_back(srcEnt);
      sCoordEtN[srcEnt] = pCoord;
    }
  }


  Tree treeNtE(boost::make_zip_iterator(boost::make_tuple( pointsNtE.begin(),indicesNtE.begin() )),
      boost::make_zip_iterator(boost::make_tuple( pointsNtE.end(),indicesNtE.end() ) ) );

  Tree treeEtN(boost::make_zip_iterator(boost::make_tuple( pointsEtN.begin(),indicesEtN.begin() )),
      boost::make_zip_iterator(boost::make_tuple( pointsEtN.end(),indicesEtN.end() ) ) );


  std::cout << "\t\t 4/4 Creating differentiation and interpolation matrices ... can take a while " << std::endl;
  matrix.numTargets = maxNumTrgEntities;
  matrix.numSources = maxNumSrcNodeEntities;
  StdVector<CF::UInt>& targetSourceIndexNtE = matrix.targetSourceIndexNtE;
  StdVector<CF::UInt>& targetSourceIndexEtN = matrix.targetSourceIndexEtN;
  StdVector<CF::UInt>& targetSourceNtE = matrix.targetSourceNtE;
  StdVector<CF::UInt>& targetSourceEtN = matrix.targetSourceEtN;
  StdVector< CF::Matrix<CF::Double> >& targetSourceFactorDiv = matrix.targetSourceFactorDiv;
  StdVector< CF::Matrix<CF::Double> >& targetSourceFactorGrad = matrix.targetSourceFactorGrad;
  StdVector< CF::Matrix<CF::Double> >& targetSourceFactorCurl = matrix.targetSourceFactorCurl;
  StdVector<CF::Double>& targetSourceNNFactor = matrix.targetSourceNNFactor;




  // creating target -> source indices
  targetSourceIndexNtE.Resize(maxNumTrgEntities + 1);
  CF::UInt indexNtE = 0;
  for(CF::UInt trgEnt = 0; trgEnt < maxNumTrgEntities; trgEnt++) {
    targetSourceIndexNtE[trgEnt] = indexNtE;
    if (globTrgEntity[trgEnt] != UnusedEntityNumber) {
      indexNtE += numNeighbors_;
    }
  }

  targetSourceIndexEtN.Resize(maxNumSrcNodeEntities + 1);
  CF::UInt indexEtN = 0;
  for(CF::UInt trgEnt = 0; trgEnt < maxNumSrcNodeEntities; trgEnt++) {
    targetSourceIndexEtN[trgEnt] = indexEtN;
    if (globSrcNodeEntity[trgEnt] != UnusedEntityNumber) {
      indexEtN += numNeighbors_;
    }
  }


  // creating target -> source factors and filling source indices
  targetSourceNtE.Resize(indexNtE);
  targetSourceEtN.Resize(indexEtN);
  targetSourceNNFactor.Resize(indexEtN);
  targetSourceNNFactor.Init();
  targetSourceIndexNtE[maxNumTrgEntities] = indexNtE;
  targetSourceIndexEtN[maxNumSrcNodeEntities] = indexEtN;
  targetSourceFactorDiv.Resize(maxNumTrgEntities);
  targetSourceFactorDiv.Init();
  targetSourceFactorGrad.Resize(maxNumTrgEntities);
  targetSourceFactorGrad.Init();
  targetSourceFactorCurl.Resize(maxNumTrgEntities);
  targetSourceFactorCurl.Init();


//#pragma omp parallel for num_threads(NUM_CFS_THREADS)
  for(CF::UInt trgEnt = 0; trgEnt < maxNumTrgEntities; trgEnt++) {
    StdVector<CF::Double> srcDist;
    srcDist.Resize(numNeighbors_);
    srcDist.Init();
    CF::UInt globEntityNumber = globTrgEntity[trgEnt];
    if (globEntityNumber != UnusedEntityNumber) {
      CF::Vector<Double> pCoord;
      // Result is defined on elements!
      trgGrid_->GetElemCentroid(pCoord, globEntityNumber,true);
      Point_3 query(pCoord[0],pCoord[1],pCoord[2]);
      Distance tr_dist;


      //TODO do not know why but it calls the get(..) method from NearstNeighbourInterpolator, altough there is one defined here...
      K_neighbor_search searchNtE(treeNtE, query, numNeighbors_);

      const CF::UInt ITrgSrcIndex = targetSourceIndexNtE[trgEnt];
      CF::UInt* srcIndices = &targetSourceNtE[ITrgSrcIndex];


      StdVector< CF::Vector<CF::Double> > neighbourCoords;
      neighbourCoords.Resize(numNeighbors_);
      CF::UInt i = 0;
      CF::Double dmax = 0.0;
      for(K_neighbor_search::iterator it = searchNtE.begin(); it != searchNtE.end(); it++, i++) {
        srcIndices[i] = boost::get<1>(it->first);
        CF::Double distance = tr_dist.inverse_of_transformed_distance(it->second);
        srcDist[i] = distance;
        if (distance > dmax) {
          dmax = distance;
        }
        neighbourCoords[i] = sCoordNtE[srcIndices[i]];
      }

      Double alpha = dmax * 1.5e+03;
      CalcLocDivergence(targetSourceFactorDiv[trgEnt], pCoord, srcDist, neighbourCoords, alpha, numNeighbors_, numEquPerEnt_, Grid_);
      CalcLocCurl(targetSourceFactorCurl[trgEnt], pCoord, srcDist, neighbourCoords, alpha, numNeighbors_, numEquPerEnt_, Grid_);
      CalcLocGradient(targetSourceFactorGrad[trgEnt], pCoord, srcDist, neighbourCoords, alpha, numNeighbors_, 1, Grid_);
    }
  }


  // for the NN-Interpolation we need to interpolate from elements(trg) to nodes(src)
  // this means we have to "switch" target- and source-meaning
  // loop over nodes
#pragma omp parallel for num_threads(NUM_CFS_THREADS)
  for(CF::UInt trgEnt = 0; trgEnt < maxNumSrcNodeEntities; trgEnt++) {
    StdVector<CF::Double> srcDistNN;
    srcDistNN.Resize(numNeighbors_);
    srcDistNN.Init();
    CF::UInt globEntityNumber = globSrcNodeEntity[trgEnt];
    CF::Vector<Double> pCoord;
    // Result is defined on nodes!
    trgGrid_->GetNodeCoordinate3D(pCoord, globEntityNumber,true);
    Point_3 query(pCoord[0],pCoord[1],pCoord[2]);
    Distance tr_dist;

    //TODO do not know why but it calls the get(..) method from NearstNeighbourInterpolator, altough there is one defined here...
    K_neighbor_search searchEtN(treeEtN, query, numNeighbors_);
    const CF::UInt ITrgSrcIndexNN = targetSourceIndexEtN[trgEnt];
    CF::UInt* srcIndicesNN = &targetSourceEtN[ITrgSrcIndexNN];
    CF::Double* srcNNFactors = &targetSourceNNFactor[ITrgSrcIndexNN];

    StdVector< CF::Vector<CF::Double> > neighbourCoordsNN;
    neighbourCoordsNN.Resize(numNeighbors_);
    CF::UInt i = 0;
    CF::Double dmax = 0.0;
    for(K_neighbor_search::iterator it = searchEtN.begin(); it != searchEtN.end(); it++, i++) {
      srcIndicesNN[i] = boost::get<1>(it->first);
      CF::Double distance = tr_dist.inverse_of_transformed_distance(it->second);
      srcNNFactors[i] = distance;
      if (distance > dmax) {
        dmax = distance;
      }
    }

    // Apply Shepard interpolation cf. Numerical Recipes 3rd ed. p. 143ff.
    // or http://www.ems-i.com/gmshelp/Interpolation/Interpolation_Schemes \
    // /Inverse_Distance_Weighted/Shepards_Method.htm
    const CF::Double R = 1.01 * dmax;
    CF::Double weights = 0.0;
    // The point which is farthest away, should at least have a non-zero
    // weight of 0.01. If we would choose R = dmax, it would not contribute
    // at all.
    for (i = 0; i < numNeighbors_; i++) {
      Double w;
      if(srcNNFactors[i] == 0){
        w = 1.0;
      } else {
        w = std::pow((R-srcNNFactors[i])/(R*srcNNFactors[i]), p_);
      }
      srcNNFactors[i] = w;
      weights += w;
    }
    for (i = 0; i < numNeighbors_; i++) {
      srcNNFactors[i] /= weights;
    }

  }

  std::cout << "\t\t Lighthill prepared!" << std::endl;
}

ResultIdList Lighthill::SetUpstreamResults(){
  ResultIdList generated;
  //we should only have one filter Result
  CF::StdVector<uuids::uuid>::iterator aIt = filterResIds.Begin();
  std::string filterResName = resultManager_->GetExtInfo(*aIt)->resultName;

  //add velocity input result to manager
  this->res1Name = params_->Get("ResultList")->Get("velocity")->Get("resultName")->As<std::string>();
  upResNames.insert(res1Name);
  res1Id = resultManager_->AddResult(res1Name,this->filterTag_);
  //set the timeline of upstream data
  resultManager_->SetTimeLine(res1Id,(*resultManager_->GetExtInfo(*aIt)->timeLine.get()));
  generated.Push_back(res1Id);

  // if an external vorticity input exists, add it to manager
  if (externVorticity_ == true){
    this->res2Name = params_->Get("ResultList")->Get("vorticity")->Get("resultName")->As<std::string>();
    upResNames.insert(res2Name);
    res2Id = resultManager_->AddResult(res2Name,this->filterTag_);
    resultManager_->SetTimeLine(res2Id,(*resultManager_->GetExtInfo(*aIt)->timeLine.get()));
    generated.Push_back(res2Id);
  }
  return generated;
}

void Lighthill::AdaptFilterResults(){
  //some checks
  ResultManager::ConstInfoPtr inInfo1 = resultManager_->GetExtInfo(upResIds[0]);
  if(!inInfo1->isValid){
    EXCEPTION("Could not validate required input result \"" << inInfo1->resultName << "\" from upstream filters.");
  }

  //we require mesh result input
  if(!inInfo1->isMeshResult){
    EXCEPTION("Lighthill/Lamb filter requires input to be defined on mesh");
  }

  //input must be node-based
  if (resultManager_->GetDefOn(upResIds[0]) == ExtendedResultInfo::ELEMENT){
    EXCEPTION("Lighthill/Lamb filter requires velocity to be defined on nodes!");
  }


  resultManager_->CopyResultData(upResIds[0],filterResIds[0]);

  //but now, we need to overwrite some things
  resultManager_->SetRegionNames(filterResIds[0],this->trgRegions_);

  CF::StdVector<std::string> dofnames;
  if( Form_ == "AeroacousticSource_LighthillSourceTerm"){
    resultManager_->SetEntryType(filterResIds[0],ExtendedResultInfo::SCALAR);
    dofnames.Push_back("x");
  }else{
    resultManager_->SetEntryType(filterResIds[0],ExtendedResultInfo::VECTOR);
    if(trgGrid_->GetDim() == 2){
      dofnames.Push_back("x");
      dofnames.Push_back("y");
    }else{
      dofnames.Push_back("x");
      dofnames.Push_back("y");
      dofnames.Push_back("z");
    }
  }
  resultManager_->SetDofNames(filterResIds[0],dofnames);




  // choose if we have element or node-results
  if(externVorticity_ == true && Form_ == "AeroacousticSource_LambVector"){
    resultManager_->SetDefOn(filterResIds[0],ExtendedResultInfo::NODE);
  }else{
    resultManager_->SetDefOn(filterResIds[0],ExtendedResultInfo::ELEMENT);
  }

  resultManager_->SetGrid(filterResIds[0],this->trgGrid_);
  resultManager_->SetMeshResult(filterResIds[0],true);

  //validate own result
  resultManager_->SetValid(filterResIds[0]);
}


}

