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
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include <algorithm>
#include <vector>
#include <tuple>
#include <unordered_map>

#ifdef USE_CGAL
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/basic.h>
#include <CGAL/Search_traits_3.h>
#include <CGAL/Search_traits_adapter.h>
#include <CGAL/point_generators_3.h>
#include <CGAL/Orthogonal_k_neighbor_search.h>
#endif
#include <boost/iterator/zip_iterator.hpp>
#include <utility>


namespace CFSDat{

Lighthill::Lighthill(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
:AeroacousticBase(numWorkers,config,resMan){

#ifndef USE_CGAL
  EXCEPTION("Lighthill needs to be compiled with USE_CGAL=ON!");
#endif

  this->filtStreamType_ = FIFO_FILTER;





  if(config->Get("ResultList")->Get("density")->Has("resultName") == false){
    density_ = config->Get("ResultList/density/uniformValue")->As<Double>();
    std::cout<<("============================================================ \n"
        "Make sure you set density ''uniformValue='' in the xml-scheme correctly, \n"
        "otherwise it is chosen to be 1.0 !! \n"
        "============================================================")<<std::endl;
    externDensity_ = false;
  }else{
    externDensity_ = true;
  }

  // ( LambVectorLighthill || LighthillSourceVector || LighthillSourceTerm )
  Form_ = config->Get("type")->As<std::string>();

  numNeighbors_ = (trgGrid_->GetDim() == 2) ? 4 : 8;


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
  this->resVelocityName = params_->Get("ResultList")->Get("velocity")->Get("resultName")->As<std::string>();
  upResNames.insert(resVelocityName);

  // if an external vorticity input exists, add it to manager
  if (externVorticity_ == true){
    this->resVorticityName = params_->Get("ResultList")->Get("vorticity")->Get("resultName")->As<std::string>();
    upResNames.insert(resVorticityName);
  }

  // if an external density input exists, add it to manager
  // density is second or third
  if (externDensity_ == true){
    this->resDensityName = params_->Get("ResultList")->Get("density")->Get("resultName")->As<std::string>();
    upResNames.insert(resDensityName);
  }

  checkSum_ = false;
  if(config->Has("sourceSum")){
    checkSum_ = config->Get("sourceSum")->As<bool>();
  }

  epsScal_ = params_->Get("RBF_Settings")->Get("epsilonScaling")->As<Double>();
  if( params_->Get("RBF_Settings")->Has("betaScaling") ){
	  betaScal_ = params_->Get("RBF_Settings")->Get("betaScaling")->As<Double>();
  }else{
	  betaScal_ = 0.0;
  }

  if( params_->Get("RBF_Settings")->Has("betaScaling") ){
	  kScal_ = params_->Get("RBF_Settings")->Get("kScaling")->As<Double>();
  }else{
	  kScal_ = 0.0;
  }
  logEps_ = false; //no logging for this class
  //logEps_ = params_->Get("RBF_Settings")->Get("logEps")->As<bool>();

}


Lighthill::~Lighthill(){

}

#ifdef USE_CGAL
bool Lighthill::UpdateResults(std::set<uuids::uuid>& upResults) {
  /// this is the vector, which will be filled with the result
  Vector<Double>& returnVec = GetOwnResultVector<Double>(filterResIds[0]);
  stepIndex_ = resultManager_->GetStepIndex(filterResIds[0]);

  Vector<Double> tempRetVec;
  if(Form_ == "AeroacousticSource_LambVector"){LambVector(tempRetVec);}
  if(Form_ == "AeroacousticSource_LighthillSourceVector"){LighthillSourceVector(tempRetVec);}
  if(Form_ == "AeroacousticSource_LighthillSourceTerm"){LighthillSourceTerm(tempRetVec, false);}
  if(Form_ == "AeroacousticSource_LighthillSourceTensor"){LighthillSourceTerm(tempRetVec, true);}

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
  
  return true;
}


void Lighthill::LambVector(Vector<Double>& tempRetVec){
  Vector<Double>& inVecVel = GetUpstreamResultVector<Double>(resVelocityId, stepIndex_);

  Matrix& matrix = matrices_[matrixIndex_];
  const UInt maxNumTrgEntities = matrix.numTargets;

  Vector<Double> LambVec;
  if(externVorticity_ == true && Form_ == "AeroacousticSource_LambVector"){
    /************* case of externally provided vorticity ***************/
    // result is defined on nodes !!!
    Vector<Double>& inVecOmega = GetUpstreamResultVector<Double>(resVorticityId, stepIndex_);
    OmegaVectorProductU(inVecVel, inVecOmega, LambVec, numEquPerEnt_);
  }else{
    /************* we compute the vorticity internally ***************/
    // result is defined on elements !!!
    StdVector< StdVector<CF::UInt> >& sourceMCurl = matrix.targetSourceIndexCurl;
    StdVector< CF::Matrix<CF::Double> >& targetSourceFactorCurl = matrix.targetSourceFactorCurl;
    UInt gridDim = inGrid_->GetDim();
    Vector<Double> Omegatemp;
    Vector<Double> Omega;

    // Curl(u)
    CalcCurl(Omegatemp, inVecVel, numEquPerEnt_, sourceMCurl, targetSourceFactorCurl, maxNumTrgEntities, gridDim);
    if(numEquPerEnt_ == 2){
      Omega = TwoDToScalar(Omegatemp);
    }else{
      Omega = Omegatemp;
    }
    StdVector< StdVector<CF::UInt> >& sourceMN2C = matrix.targetSourceIndexNtE;

    // interpolate input-velocity to cell-centers, because the calculated curl is now defined on elements
    Vector<Double> velE;
    Node2Cell(velE, inVecVel, numEquPerEnt_, sourceMN2C, maxNumTrgEntities, gridDim);
    OmegaVectorProductU(velE, Omega, LambVec, numEquPerEnt_);
  }
  tempRetVec = LambVec;
}

void Lighthill::LighthillTensor(Vector<Double>& tempRetVec){
  Vector<Double>& inVecVel = GetUpstreamResultVector<Double>(resVelocityId, stepIndex_);
  Vector<Double> LHTensor;
  if(externDensity_ == true){
    Vector<Double>& inVecDensity = GetUpstreamResultVector<Double>(resDensityId, stepIndex_);
    TensorProduct(LHTensor, inVecVel, inVecVel, numEquPerEnt_, 1.0, inVecDensity);
  } else {
    TensorProduct(LHTensor, inVecVel, inVecVel, numEquPerEnt_, 1.0);
  }




  Matrix& matrix = matrices_[matrixIndex_];
  StdVector< StdVector<CF::UInt> >& sourceMDiv = matrix.targetSourceIndexDiv;
  StdVector< CF::Matrix<CF::Double> >& targetSourceFactorDiv = matrix.targetSourceFactorDiv;
  const UInt maxNumTrgEntities = matrix.numTargets;

  /**************** Div(renoldsStress) *****************************/
  Vector<Double> firstDivergence;
  CalcTensorDivergence(firstDivergence, LHTensor, numEquPerEnt_, sourceMDiv, targetSourceFactorDiv, maxNumTrgEntities);

  tempRetVec = firstDivergence;
}



void Lighthill::LighthillSourceVector(Vector<Double>& tempRetVec){
  Vector<Double>& inVecVel = GetUpstreamResultVector<Double>(resVelocityId, stepIndex_);


  Matrix& matrix = matrices_[matrixIndex_];
  StdVector< StdVector<CF::UInt> >& sourceMN2C = matrix.targetSourceIndexNtE;
  StdVector< CF::Matrix<CF::Double> >& targetSourceFactorGrad = matrix.targetSourceFactorGrad;
  UInt gridDim = inGrid_->GetDim();
  Vector<Double> Omegatemp;
  Vector<Double> Omega;
  const UInt maxNumTrgEntities = matrix.numTargets;

  // interpolate input-velocity to cell-centers
  Vector<Double> velE;
  Node2Cell(velE, inVecVel, numEquPerEnt_, sourceMN2C, maxNumTrgEntities, gridDim);

  /**************** Term 1: LambVector *****************************/
  Vector<Double> LambVec;
  this->LambVector(LambVec);


  /**************** Term 2: Grad(1/2 u*u) *****************************/
  Vector<Double> uuN;
  ScalarProduct(uuN, inVecVel, inVecVel, numEquPerEnt_, 0.5);


  StdVector< StdVector<CF::UInt> >& sourceM = matrix.targetSourceIndexGrad;

  Vector<Double> GraduuN;
  CalcGradient(GraduuN, uuN, 1, sourceM, targetSourceFactorGrad, maxNumTrgEntities, numEquPerEnt_);

  /**************** Term 1 + Term 2 *****************************/
  tempRetVec = GraduuN + LambVec;
}



void Lighthill::LighthillSourceTerm(Vector<Double>& tempRetVec, bool isTensorForm){
  CF::StdVector<UInt> eqnNums;

  Matrix& matrix = matrices_[matrixIndex_];
  StdVector< StdVector<CF::UInt> >& sourceMDiv = matrix.targetSourceIndexDiv;
  StdVector< StdVector<CF::UInt> >& targetSourceIndexEtN = matrix.targetSourceIndexEtN;
  //StdVector<CF::UInt>& targetSourceEtN = matrix.targetSourceEtN;
  StdVector< CF::Matrix<CF::Double> >& targetSourceNNFactor = matrix.targetSourceNNFactor;
  StdVector< CF::Matrix<CF::Double> >& targetSourceFactorDiv = matrix.targetSourceFactorDiv;
  Vector<Double> Omegatemp;
  Vector<Double> Omega;
  const UInt maxNumTrgEntities = matrix.numTargets;
  const UInt maxNumSrcEntities = matrix.numSources;


  Vector<Double> LightVec;
  if (isTensorForm) {
    this->LighthillTensor(LightVec);
  } else {
    this->LighthillSourceVector(LightVec);
  }

  // interpolate from element-centroid to nodes
  Vector<Double> LightVecN;
  NearestNeighbourLight(LightVecN, LightVec, numEquPerEnt_, targetSourceIndexEtN, targetSourceNNFactor, maxNumSrcEntities);


  /**************** Div(Term 1 + Term 2) *****************************/
  CalcDivergence(tempRetVec, LightVecN, numEquPerEnt_, sourceMDiv, targetSourceFactorDiv, maxNumTrgEntities);
}


// cgal copy and paste stuff for association of points to CF::UInt from
// http://doc.cgal.org/Manual/4.2/doc_html/cgal_manual/Spatial_searching/Chapter_main.html#Subsection_61.3.1
//
// later on this should be refactored into KNNSearch.hh and KNNSearch.cc
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::Point_3 Point_3;
typedef std::tuple<Point_3,CF::UInt> Point_and_int;

//definition of the property map
struct My_point_property_map{
  typedef Point_3 value_type;
  typedef const value_type& reference;
  typedef const Point_and_int& key_type;
  typedef boost::lvalue_property_map_tag category;
};

//get function for the property map
inline My_point_property_map::reference
get(My_point_property_map,My_point_property_map::key_type p)
{return std::get<0>(p);}

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
  uuids::uuid upRes = resVelocityId;
  inGrid_ = resultManager_->GetExtInfo(upRes)->ptGrid;

  scrMap_ = resultManager_->GetEqnMap(upRes);
  trgMap_ = resultManager_->GetEqnMap(filterResIds[0]);
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upRes);
//TODO wrong number of equations per entity for VELOCITY if VORTICITY is read-in by EnSight as CFSVarName="fluidMechVorticity"
  numEquPerEnt_ = scrMap_->GetNumEqnPerEnt();
  if(numEquPerEnt_ != inGrid_->GetDim()) EXCEPTION("Grid (dim="<<inGrid_->GetDim()<<") and Values of "<< resultManager_->GetResultName(upRes) <<" (dim="<<numEquPerEnt_<<") don't have the same dimension!!");
  if(numEquPerEnt_ == 1) EXCEPTION("Input (velocity) is a scalar, it should be a vector!!");


  // input MUST be defined on nodes !!
  //bool inElems = inInfo->definedOn == ExtendedResultInfo::ELEMENT;

  const CF::UInt maxNumSrcNodeEntities = scrMap_->GetNumEntities(); //velocity...upResIds[0]
  StdVector<CF::UInt> globSrcNodeEntity;
  GetUsedMappedEntities(scrMap_, globSrcNodeEntity, srcRegions_, inGrid_);
  const CF::UInt numSrcNodeEntities = CountUsedEntities(globSrcNodeEntity);
  // the following map is needed because the inVec in Run() doesn't know
  // which nodeNumber belongs to which entry...we also can't hardcode it
  // because we have dynamical storage of the neighbours, means we can not
  // predict the size
  std::unordered_map<UInt, UInt> sEnt;
  for(UInt i = 0; i < globSrcNodeEntity.GetSize(); ++i){
    sEnt[globSrcNodeEntity[i]] = i + 1;
  }

  differentiators_.Push_back(this);
  matrixIndex_ = matrices_.GetSize();
  matrices_.Resize(matrixIndex_ + 1);
  Matrix& matrix = matrices_[matrixIndex_];


  std::cout << "\t\t 2/4 Obtaining target entities (elements)" << std::endl;

  const CF::UInt maxNumTrgEntities = trgMap_->GetNumEntities();
  StdVector<CF::UInt> globTrgEntity;
  GetUsedMappedEntities(trgMap_, globTrgEntity, trgRegions_, trgGrid_);
  const CF::UInt numTrgEntities = CountUsedEntities(globTrgEntity);

  std::cout << "\t\t\t Differentiator is dealing with " << numSrcNodeEntities <<
      " source " << ("nodes") << " and "<< numTrgEntities << " target elements" << std::endl;


  std::cout << "\t\t 3/4 Creating interpolation matrices for derivatives " << std::endl;
  matrix.numTargets = maxNumTrgEntities;
  matrix.numSources = maxNumSrcNodeEntities;
  StdVector< StdVector<CF::UInt> >& sourceMDiv = matrix.targetSourceIndexDiv;
  StdVector< StdVector<CF::UInt> >& sourceMGrad = matrix.targetSourceIndexGrad;
  StdVector< StdVector<CF::UInt> >& sourceMCurl = matrix.targetSourceIndexCurl;
  StdVector< StdVector<CF::UInt> >& targetSourceIndexEtN = matrix.targetSourceIndexEtN;
  StdVector< StdVector<CF::UInt> >& sourceMN2C = matrix.targetSourceIndexNtE;
  StdVector<CF::UInt>& targetSourceEtN=matrix.targetSourceEtN;
  StdVector<CF::UInt>& targetSourceNtE=matrix.targetSourceNtE;


  StdVector< CF::Matrix<CF::Double> >& targetSourceFactorDiv = matrix.targetSourceFactorDiv;
  StdVector< CF::Matrix<CF::Double> >& targetSourceFactorCurl = matrix.targetSourceFactorCurl;
  StdVector< CF::Matrix<CF::Double> >& targetSourceFactorGrad = matrix.targetSourceFactorGrad;
  StdVector< CF::Matrix<CF::Double> >& targetSourceNNFactor = matrix.targetSourceNNFactor;


  sourceMN2C.Resize(numTrgEntities);
  targetSourceFactorDiv.Resize(numTrgEntities);
  sourceMDiv.Resize(numTrgEntities);
  targetSourceFactorCurl.Resize(numTrgEntities);
  sourceMCurl.Resize(numTrgEntities);
  targetSourceFactorGrad.Resize(numTrgEntities);
  sourceMGrad.Resize(numTrgEntities);
  targetSourceNNFactor.Resize(numSrcNodeEntities);
  targetSourceIndexEtN.Resize(numSrcNodeEntities);



  StdVector<RegionIdType> rId;
  rId.Init(0);
  StdVector<UInt> numNodesTrgRegs;
  numNodesTrgRegs.Init(0);
  std::set<std::string>::const_iterator destRegIt = this->trgRegions_.begin();
  for(; destRegIt != this->trgRegions_.end(); ++destRegIt ) {
    RegionIdType r = trgGrid_->GetRegion().Parse(*destRegIt);
    rId.Push_back(r);
    UInt cache = trgGrid_->GetNumNodes(r);
    numNodesTrgRegs.Push_back(cache);
  }

  StdVector<RegionIdType> src_rId;
  src_rId.Init(0);
  StdVector<UInt> numNodesSrcRegs;
  numNodesSrcRegs.Init(0);
  UInt noRegions = 0;
  std::set<std::string>::const_iterator srcRegIt = this->srcRegions_.begin();
  for(; srcRegIt != this->srcRegions_.end(); ++srcRegIt ) {
    RegionIdType r = inGrid_->GetRegion().Parse(*srcRegIt);
    src_rId.Push_back(r);
    UInt cache = inGrid_->GetNumNodes(r);
    numNodesSrcRegs.Push_back(cache);
    noRegions += 1;
  }

  bool inElems = inInfo->definedOn == ExtendedResultInfo::ELEMENT;
  for(UInt i=0; i < noRegions; i++){
	    if(numNodesSrcRegs[i]!=numNodesTrgRegs[i]){
	      std::cout << "\t\t\t Differentiator is dealing with " << numNodesSrcRegs[i] <<
	                   " source " << (inElems ? "elements" : "nodes") << " and "<< numNodesTrgRegs[i] << " target nodes" << std::endl;
	      EXCEPTION("Source and target mesh of involved regions must be consistent. Use a subsequent interpolation filter, if required otherwise.");
	    }
  }


  /********************************************************************
   ********** Prepare Differentiators and Node2Cell *******************
   ********************************************************************/

  for(CF::UInt trgEnt = 0; trgEnt < maxNumTrgEntities; trgEnt++) {
    StdVector<CF::UInt> sM, nL, nodeList;
    CF::UInt globEntityNumber = globTrgEntity[trgEnt];
    if (globEntityNumber != UnusedEntityNumber) {
        targetSourceNtE.Push_back(globEntityNumber);
        CF::Vector<Double> trgCoord;
        // Result is defined on elements!
        inGrid_->GetElemCentroid(trgCoord, globEntityNumber,true);


        inGrid_->GetElemNodes(nodeList, globEntityNumber);
        nL = nodeList;

        StdVector< CF::Vector<CF::Double> > neighbourCoords;
        StdVector<CF::Double> srcDist;
        CF::Vector<CF::Double> tmpCoords;
        UInt numSrcPoints;
        Double maxd;
        if(Form_ != "AeroacousticSource_LambVector"){
        /****************** 1) Divergence ****************************************/
        sM.Clear(false);
        maxd = 0.0;
        for(UInt i = 0; i < nodeList.GetSize(); ++i){
          if(!sM.Contains(sEnt[nodeList[i]])){
            sM.Push_back(sEnt[nodeList[i]]);
            inGrid_->GetNodeCoordinate(tmpCoords, nodeList[i], false);
            neighbourCoords.Push_back(tmpCoords);
            if(tmpCoords.GetSize() == 2) tmpCoords.Push_back(0.0);
            CF::Double d = trgCoord.NormL2(tmpCoords);
            srcDist.Push_back(d);
            if(maxd < d) maxd = d;
          }
        }
        numSrcPoints = srcDist.GetSize();
        CF::Matrix<CF::Double> tsFDiv;
        while( !CalcLocDivergence(tsFDiv, trgCoord, maxd, srcDist, neighbourCoords,
            numSrcPoints, numEquPerEnt_, inGrid_, epsScal_, betaScal_, kScal_, logEps_)){
          // find furthest point
          Double d = 0.0;
          UInt maxId = 0;
          for(UInt i = 0; i < srcDist.GetSize(); ++i){
            if(d < srcDist[i]){
              maxId = i;
            }
          }
          sM.Erase(maxId);
          srcDist.Erase(maxId);
          numSrcPoints -= 1;
          neighbourCoords.Erase(maxId);

          if( sM.GetSize() < 3) EXCEPTION("Patch-Problem, modify epsilon Divergence!");

        }// while local deriv is false
        targetSourceFactorDiv[trgEnt] = tsFDiv;
        sourceMDiv[trgEnt] = sM;

        /****************** 2) Gradient ****************************************/
        sM.Clear(false);
        neighbourCoords.Clear(false);
        srcDist.Clear(false);
        maxd = 0.0;
        for(UInt i = 0; i < nodeList.GetSize(); ++i){
          if(!sM.Contains(sEnt[nodeList[i]])){
            sM.Push_back(sEnt[nodeList[i]]);
            inGrid_->GetNodeCoordinate(tmpCoords, nodeList[i], false);
            neighbourCoords.Push_back(tmpCoords);
            if(tmpCoords.GetSize() == 2) tmpCoords.Push_back(0.0);
            CF::Double d = trgCoord.NormL2(tmpCoords);
            srcDist.Push_back(d);
            if(maxd < d) maxd = d;
          }
        }

        numSrcPoints = srcDist.GetSize();
        CF::Matrix<CF::Double> tsFGrad;
        while( !CalcLocGradient(tsFGrad, trgCoord, maxd, srcDist, neighbourCoords,
            numSrcPoints, 1, inGrid_, epsScal_, betaScal_, kScal_, logEps_)){
          // find furthest point
          Double d = 0.0;
          UInt maxId = 0;
          for(UInt i = 0; i < srcDist.GetSize(); ++i){
            if(d < srcDist[i]){
              maxId = i;
            }
          }
          sM.Erase(maxId);
          srcDist.Erase(maxId);
          numSrcPoints -= 1;
          neighbourCoords.Erase(maxId);

          if( sM.GetSize() < 3) EXCEPTION("Patch-Problem, modify epsilon Gradient!");

        }// while local deriv is false
        targetSourceFactorGrad[trgEnt] = tsFGrad;
        sourceMGrad[trgEnt] = sM;

        }
        if(externVorticity_ == false){
        /***************** 3) Curl ******************************************/
        sM.Clear(false);
        neighbourCoords.Clear(false);
        srcDist.Clear(false);
        maxd = 0.0;
        for(UInt i = 0; i < nodeList.GetSize(); ++i){
          if(!sM.Contains(sEnt[nodeList[i]])){
            sM.Push_back(sEnt[nodeList[i]]);
            inGrid_->GetNodeCoordinate(tmpCoords, nodeList[i], false);
            neighbourCoords.Push_back(tmpCoords);
            if(tmpCoords.GetSize() == 2) tmpCoords.Push_back(0.0);
            CF::Double d = trgCoord.NormL2(tmpCoords);
            srcDist.Push_back(d);
            if(maxd < d) maxd = d;
          }
        }
        numSrcPoints = srcDist.GetSize();
        CF::Matrix<CF::Double> tsFCurl;
        while( !CalcLocCurl(tsFCurl, trgCoord, maxd, srcDist, neighbourCoords,
            numSrcPoints, numEquPerEnt_, inGrid_, epsScal_, betaScal_, kScal_, logEps_)){
          // find furthest point
          Double d = 0.0;
          UInt maxId = 0;
          for(UInt i = 0; i < srcDist.GetSize(); ++i){
            if(d < srcDist[i]){
              maxId = i;
            }
          }
          sM.Erase(maxId);
          srcDist.Erase(maxId);
          numSrcPoints -= 1;
          neighbourCoords.Erase(maxId);

          if( sM.GetSize() < 3) EXCEPTION("Patch-Problem, modify epsilon Gradient!");

        }// while local deriv is false
        targetSourceFactorCurl[trgEnt] = tsFCurl;
        sourceMCurl[trgEnt] = sM;


        /****************** 4) Node2Cell ***********************************/
        sM.Clear(false);
        for(UInt i = 0; i < nL.GetSize(); ++i){
          sM.Push_back(sEnt[nodeList[i]]);
        }
        sourceMN2C[trgEnt] = sM;

        }
    }
  }




  /********************************************************************
   ********** Prepare NearestNeighbour ********************************
   ********************************************************************/
  // unfortunately the Grid-class is not able to return all elements, which
  // share a common node, therefore we use a nn-search
  if(externVorticity_ == false && Form_ != "AeroacousticSource_LambVector"){
  std::cout << "\t\t 4/4 Creating interpolation matrices for interpolators " << std::endl;
//#pragma omp parallel num_threads(CFS_NUM_THREADS)
//  {
  std::vector<Point_3> pointsEtN; //EtN...elements to nodes
  std::vector<CF::UInt> indicesEtN;
  for(CF::UInt srcEnt = 0; srcEnt < maxNumTrgEntities; srcEnt++) {
    CF::UInt globEntityNumber = globTrgEntity[srcEnt];
    if (globEntityNumber != UnusedEntityNumber) {
      CF::Vector<Double> pCoord;
//#pragma omp critical
//      {
      inGrid_->GetElemCentroid(pCoord, globEntityNumber, true);
//      }
      Point_3 point(pCoord[0],pCoord[1],pCoord[2]);
      pointsEtN.push_back(point);
      indicesEtN.push_back(srcEnt);
    }
  }

  StdVector<Point_and_int> searchPoints;
  searchPoints.Reserve(pointsEtN.size());
  size_t i = 0;
  for(const auto& pt : pointsEtN)
    searchPoints.Push_back(std::make_tuple(pt, indicesEtN[i++]));

  Tree treeEtN(searchPoints.begin(), searchPoints.end());


  for(CF::UInt trgEnt = 0; trgEnt < maxNumSrcNodeEntities; trgEnt++) {
    CF::UInt globEntityNumber = globSrcNodeEntity[trgEnt];
    if (globEntityNumber != UnusedEntityNumber) {
      targetSourceEtN.Push_back(globEntityNumber);


      CF::Vector<Double> trgCoord;
      // Result is now defined on elements, since it's an interpolation
      // from elements to nodes
//#pragma omp critical
//      {
      inGrid_->GetNodeCoordinate(trgCoord, globEntityNumber, true);
//      }
      if(trgCoord.GetSize() == 2) trgCoord.Push_back(0.0);


      Point_3 query(trgCoord[0],trgCoord[1],trgCoord[2]);
      Distance tr_dist;

      UInt numNeighbors = (numEquPerEnt_ == 2)? 4 : 8;
      //TODO do not know why but it calls the get(..) method from
      // NearstNeighbourInterpolator, altough there is one defined here...
      K_neighbor_search searchEtN(treeEtN, query, numNeighbors);


      StdVector< CF::Vector<CF::Double> > neighbourCoords;
      StdVector<CF::Double> srcDist;
      CF::Vector<CF::Double> tmpCoords;

      StdVector<CF::UInt> sM;
      CF::Double dmax = 0.0;
      for(K_neighbor_search::iterator it = searchEtN.begin(); it != searchEtN.end(); it++) {
        // the following +1 is necessary because then it's the entity-number and not index
        sM.Push_back(std::get<1>(it->first) + 1 );
        CF::Double distance = tr_dist.inverse_of_transformed_distance(it->second);
        srcDist.Push_back(distance );
        if (distance > dmax) {
          dmax = distance;
        }

      }

      CF::Matrix<CF::Double> tsNN;
      tsNN.Resize(numNeighbors, 1);
      // Apply Shepard interpolation cf. Numerical Recipes 3rd ed. p. 143ff.
      // or http://www.ems-i.com/gmshelp/Interpolation/Interpolation_Schemes \
      // /Inverse_Distance_Weighted/Shepards_Method.htm
      const CF::Double R = 1.01 * dmax;
      CF::Double weights = 0.0;
      // The point which is farthest away, should at least have a non-zero
      // weight of 0.01. If we would choose R = dmax, it would not contribute
      // at all.

      for (UInt i = 0; i < numNeighbors; i++) {
        Double w;
        if(srcDist[i] == 0){
          w = 1.0;
        } else {
          w = std::pow((R-srcDist[i])/(R*srcDist[i]), 4);
        }
        tsNN[i][0] = w;
        weights += w;
      }
      for (UInt i = 0; i < numNeighbors; i++) {
        tsNN[i][0] /= weights;
      }

      targetSourceNNFactor[trgEnt] = tsNN;
      targetSourceIndexEtN[trgEnt] = sM;
    }
  }
//  }//omp parallel


  }
  std::cout << "\t\t Lighthill prepared!" << std::endl;
}

ResultIdList Lighthill::SetUpstreamResults(){
  ResultIdList generated = SetDefaultUpstreamResults();
  resVelocityId = upResNameIds[resVelocityName];
  resVorticityId = upResNameIds[resVorticityName];
  resDensityId = upResNameIds[resDensityName];
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
    EXCEPTION("Lighthill/Lamb filter requires velocity to be defined on nodes!\n" <<
        "  Use a Cell2Node Filter to transform the filter input.");
  }


  resultManager_->CopyResultData(upResIds[0],filterResIds[0]);

  //but now, we need to overwrite some things
  resultManager_->SetRegionNames(filterResIds[0],this->trgRegions_);

  CF::StdVector<std::string> dofnames;
  if( Form_ == "AeroacousticSource_LighthillSourceTerm" ||
      Form_ == "AeroacousticSource_LighthillSourceTensor"){
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



std::cout<<"AdaptFilterResults Form_ = "<<Form_<<std::endl;
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

#else

ResultIdList Lighthill::SetUpstreamResults(){
  return SetDefaultUpstreamResults();
}
bool Lighthill::UpdateResults(std::set<uuids::uuid>& upResults) {
	return true;
}
void Lighthill::PrepareCalculation(){}
void Lighthill::AdaptFilterResults(){}
#endif


}

