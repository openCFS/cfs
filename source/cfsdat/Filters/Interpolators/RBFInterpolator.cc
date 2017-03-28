// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     RBFInterpolator.hh
 *       \brief    <Description>
 *
 *       \date     Sep 8, 2016
 *       \author   kroppert
 */
//================================================================================================


#include "RBFInterpolator.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
//#include "cfsdat/Utils/KNNSearch.hh"
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

RBFInterpolator::RBFInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
:MeshFilter(numWorkers,config,resMan){

  this->filtStreamType_ = FIFO_FILTER;
  inDim_ = 0;
  p_ = config->Get("scheme")->Get("interpolationExponent")->As<UInt>();

  numNN_ = 18;
  numNW_ = 13;
  noSlip_ = false;
  if(config->Has("noSlipWall")){noSlip_ = true;}

}

RBFInterpolator::~RBFInterpolator(){

}

bool RBFInterpolator::Run(){
  // we deactivate every result, except for our own
  std::set<uuids::uuid> activeResults = resultManager_->GetActiveResults();
  std::set<uuids::uuid>::iterator aIter = activeResults.begin();

  for(; aIter != activeResults.end(); ++aIter){
    if(filterResIds.Find(*aIter) == -1){
      WARN(" There are still active results when reaching the interpolation filter. This indicates an unexpected use of the pipeline.")
    }
    resultManager_->DeactivateResult(*aIter);
  }
  Double aTF = resultManager_->GetStepValue(filterResIds[0]);
  resultManager_->SetTimeValue(upResIds[0],aTF);
  // now we deactivate our own result and activate the others
  resultManager_->ActivateResult(upResIds[0]);

  //now we call for upstream data in each source
  CF::StdVector< str1::shared_ptr<BaseFilter> >::iterator srcIter =  sources_.Begin();
  for(; srcIter != sources_.End() ; srcIter++){
    // should we check here anything for success?
    (*srcIter)->Run();
  }


  CF::StdVector<UInt> eqnNums;

  /// this is the vector, which will be filled with the interpolation result
  Vector<Double>& returnVec = resultManager_->GetResultVector<Double>(filterResIds[0],eqnNums);
  returnVec.Init();

  // vector, containing the source data values
  Vector<Double>& inVec = resultManager_->GetResultVector<Double>(upResIds[0],eqnNums);


  Matrix& matrix = matrices_[matrixIndex_];
  const UInt maxNumTrgEntities = matrix.numTargets;
  StdVector<CF::UInt>& targetSourceIndex = matrix.targetSourceIndex;
  StdVector<CF::UInt>& targetSource = matrix.targetSource;
  StdVector<CF::Double>& targetSourceFactor = matrix.targetSourceFactor;
  StdVector<CF::Double>& targetSourceFactor2 = matrix.targetSourceFactor2;
  StdVector< CF::Matrix<Double> >& targetRBFInv = matrix.targetRBFInvMat;



  RBFInterpolation(returnVec, inVec, numEquPerEnt_, targetSource, targetSourceIndex, numNN_, targetRBFInv, targetSourceFactor, targetSourceFactor2, maxNumTrgEntities);


  resultManager_->ActivateResult(filterResIds[0]);

  //now deactivate own upstream results
  for(UInt aRes=0;aRes<upResIds.GetSize();aRes++){
    resultManager_->DeactivateResult(upResIds[aRes]);
  }

  return true;
}

/*
CF::UInt RBFInterpolator::CountUsedEntities(const StdVector<CF::UInt>& entities) {
  const CF::UInt size = entities.GetSize();
  CF::UInt numEntities = 0;
#pragma omp parallel for reduction(+:numEntities) num_threads(NUM_CFS_THREADS)
  for(CF::UInt inEnt = 0; inEnt < size; inEnt++) {
    if (entities[inEnt] != UnusedEntityNumber) {
      numEntities++;
    }
  }
  return numEntities;
}

void RBFInterpolator::GetUsedMappedEntities(const str1::shared_ptr<EqnMapSimple>& map,
                                                  StdVector<CF::UInt>& entities,
                                                  const std::set<std::string>& regions,
                                                  Grid* grid) {
  bool useElems = map->GetMapType() == ExtendedResultInfo::ELEMENT;

  const CF::UInt maxNumEntities = map->GetNumEntities();
  entities.Clear();
  entities.Resize(maxNumEntities, UnusedEntityNumber);

  std::set<std::string>::iterator sRegIter = regions.begin();
  for(;sRegIter != regions.end();++sRegIter){
    StdVector<CF::UInt> regEntities;
    if (useElems) {
      grid->GetElemNumsByName(regEntities,*sRegIter);
    } else {
      grid->GetNodesByName(regEntities, *sRegIter);
    }
    const UInt size = regEntities.GetSize();
#pragma omp parallel for num_threads(NUM_CFS_THREADS)
    for (UInt eIter = 0; eIter < size; ++eIter) {
      CF::UInt entityNumber = regEntities[eIter];
      entities[map->GetEntityIndex(entityNumber)] = entityNumber;
    }
  }
}
*/


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





CF::StdVector<RBFInterpolator*> RBFInterpolator::interpolators_;
CF::StdVector<RBFInterpolator::Matrix> RBFInterpolator::matrices_;



void RBFInterpolator::PrepareCalculation(){
  std::cout << "\t ---> RBFInterpolator preparing for interpolation" << std::endl;

  std::cout << "\t\t 1/4 Obtaining source entities " << std::endl;
  uuids::uuid upRes = upResIds[0];
  inGrid_ = resultManager_->GetExtInfo(upRes)->ptGrid;
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  scrMap_ = resultManager_->GetResultAdapter(upRes)->mapping;
  numEquPerEnt_ = scrMap_->GetNumEqnPerEnt();
  bool inElems = inInfo->definedOn == ExtendedResultInfo::ELEMENT;

  const CF::UInt maxNumSrcEntities = scrMap_->GetNumEntities();
  StdVector<CF::UInt> globSrcEntity;
  GetUsedMappedEntities(scrMap_, globSrcEntity, srcRegions_, inGrid_);
  const CF::UInt numSrcEntities = CountUsedEntities(globSrcEntity);
  if (numSrcEntities < numNN_) {
    numNN_ = numSrcEntities;
  }

  trgMap_ = resultManager_->GetResultAdapter(filterResIds[0])->mapping;


  interpolators_.Push_back(this);
  matrixIndex_ = matrices_.GetSize();
  matrices_.Resize(matrixIndex_ + 1);
  Matrix& matrix = matrices_[matrixIndex_];

  std::cout << "\t\t 2/4 Obtaining target entities " << std::endl;

  const CF::UInt maxNumTrgEntities = trgMap_->GetNumEntities();
  StdVector<CF::UInt> globTrgEntity;
  GetUsedMappedEntities(trgMap_, globTrgEntity, trgRegions_, trgGrid_);
  const CF::UInt numTrgEntities = CountUsedEntities(globTrgEntity);

  std::cout << "\t\t\t Interpolator is dealing with " << numSrcEntities <<
               " source " << (inElems ? "elements" : "nodes") << " and "<< numTrgEntities << " target nodes" << std::endl;

  std::cout << "\t\t 3/4 Creating search tree " << std::endl;
  std::vector<Point_3> points;
  std::vector<CF::UInt> indices;
  StdVector< CF::Vector<Double> > sCoord;
  sCoord.Resize(maxNumSrcEntities);
  for(CF::UInt srcEnt = 0; srcEnt < maxNumSrcEntities; srcEnt++) {
    CF::UInt globEntityNumber = globSrcEntity[srcEnt];
    if (globEntityNumber != UnusedEntityNumber) {
      CF::Vector<Double> pCoord;
      if (inElems) {
        inGrid_->GetElemCentroid(pCoord,globEntityNumber,true);
      } else {
        inGrid_->GetNodeCoordinate3D(pCoord, globEntityNumber);
      }
      Point_3 point(pCoord[0],pCoord[1],pCoord[2]);
      points.push_back(point);
      indices.push_back(srcEnt);
      sCoord[srcEnt] = pCoord;
    }
  }
  Tree tree(boost::make_zip_iterator(boost::make_tuple( points.begin(),indices.begin() )),
            boost::make_zip_iterator(boost::make_tuple( points.end(),indices.end() ) ) );


  std::cout << "\t\t 4/4 Creating interpolation matrix " << std::endl;
  matrix.numTargets = maxNumTrgEntities;
  StdVector<CF::UInt>& targetSourceIndex = matrix.targetSourceIndex;
  StdVector<CF::UInt>& targetSource = matrix.targetSource;
  StdVector<CF::Double>& targetSourceFactor = matrix.targetSourceFactor;
  StdVector<CF::Double>& targetSourceFactor2 = matrix.targetSourceFactor2;
  StdVector< CF::Matrix<Double> >& targetInvMat = matrix.targetRBFInvMat;

  // creating target -> source indices
  targetSourceIndex.Resize(maxNumTrgEntities + 1);
  CF::UInt index = 0;
  for(CF::UInt trgEnt = 0; trgEnt < maxNumTrgEntities; trgEnt++) {
    targetSourceIndex[trgEnt] = index;
    if (globTrgEntity[trgEnt] != UnusedEntityNumber) {
      index += numNN_;
    }
  }
  targetSourceIndex[maxNumTrgEntities] = index;
  targetInvMat.Resize(maxNumTrgEntities);

  // creating target -> source factors and filling source indices
  if (numNN_ > 1) {
    targetSource.Resize(index);
    targetSourceFactor.Resize(index);
    targetSourceFactor.Init();
    targetSourceFactor2.Resize(index);
    targetSourceFactor2.Init();
  }


#pragma omp parallel for num_threads(NUM_CFS_THREADS)
  for(CF::UInt trgEnt = 0; trgEnt < maxNumTrgEntities; trgEnt++) {
    StdVector<CF::Double> srcDist;
    srcDist.Resize(numNN_);
    srcDist.Init();
    CF::UInt globEntityNumber = globTrgEntity[trgEnt];
    if (globEntityNumber != UnusedEntityNumber) {
      CF::Vector<Double> pCoord;
      trgGrid_->GetNodeCoordinate3D(pCoord, globEntityNumber);
      Point_3 query(pCoord[0],pCoord[1],pCoord[2]);
      Distance tr_dist;
//TODO do not know why but it calls the get(..) method from NearstNeighbourInterpolator, altough there is one defined here...
      K_neighbor_search search(tree, query, numNN_);
      if (numNN_ > 1) {
        const CF::UInt ITrgSrcIndex = targetSourceIndex[trgEnt];
        CF::UInt* srcIndices = &targetSource[ITrgSrcIndex];
        CF::Double* srcFactors = &targetSourceFactor[ITrgSrcIndex];
        CF::Double* srcFactors2 = &targetSourceFactor2[ITrgSrcIndex];

        StdVector< CF::Vector<CF::Double> > neighbourCoords;
        neighbourCoords.Resize(numNN_);
        CF::UInt i = 0;
        CF::Double dmax = 0.0;
        for(K_neighbor_search::iterator it = search.begin(); it != search.end(); it++, i++) {
          srcIndices[i] = boost::get<1>(it->first);
          CF::Double distance = tr_dist.inverse_of_transformed_distance(it->second);
          srcDist[i] = distance;
          if (distance > dmax) {
            dmax = distance;
          }
          neighbourCoords[i] = sCoord[srcIndices[i]];
        }

        Double alpha = dmax * 1.5e+03;

        CF::Matrix<Double> InvLoc;
        CalcLocRBFInv(InvLoc, neighbourCoords, alpha, numNN_, trgGrid_);
        targetInvMat[trgEnt] = InvLoc;


        Double r_k;
        Double temp;
        Double r_Wk = srcDist[numNW_]; //former dmax influence radius of the weight function
        Double W_denom = 0.0;
        //r_Wk = targetSourceDist[srcIndices[numNN_]]; //influence radius of the weight function
        //now we can perform the computation of the weight function
        for (UInt srcIter = 0; srcIter < numNN_; ++srcIter){

          W_denom = 0.0;
          for (UInt k = 0; k < numNN_; ++k){
            r_k = srcDist[k];
            temp = (r_Wk - r_k);
            if (temp < 0.0){temp = 0.0;}
            if (r_k == 0.0){
              W_denom += 0.0;
            }else{
              W_denom += pow(temp / (r_Wk * r_k), Double(p_));
            }
          }
          r_k = srcDist[srcIter];
          temp = (r_Wk - r_k);
          if (temp < 0.0){temp = 0.0;}
          if (r_k == 0.0){
            srcFactors[srcIter] = 1.0 / W_denom; //W_k_bar
          }else{
            srcFactors[srcIter] = pow(temp / (r_Wk * r_k), Double(p_)) / W_denom;//W_k_bar
          }
          srcFactors2[srcIter] = pow(1.0 - srcDist[srcIter]/alpha, 2.0);
        }


      } else {
        targetSourceIndex[trgEnt] = boost::get<1>(search.begin()->first);
      }
    } else {
      if (numNN_ == 1) {
        targetSourceIndex[trgEnt] = UnusedEntityNumber;
      }
    }
  }


  std::cout << "\t\t Interpolation prepared!" << std::endl;
}

ResultIdList RBFInterpolator::SetUpstreamResults(){
  ResultIdList generated;
  //we should only have one filter Result
  CF::StdVector<uuids::uuid>::iterator aIt = filterResIds.Begin();
  std::string filterResName = resultManager_->GetExtInfo(*aIt)->resultName;

  //add input result to manager
  std::string inRes = params_->Get("singleResult")->Get("inputQuantity")->Get("resultName")->As<std::string>();
  uuids::uuid newId = resultManager_->AddResult(inRes,this->filterTag_);

  //set the timeline of upstream data if already set
  resultManager_->SetTimeLine(newId,(*resultManager_->GetExtInfo(*aIt)->timeLine.get()));
  generated.Push_back(newId);

  return generated;

}

void RBFInterpolator::AdaptFilterResults(){
  //some checks
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  if(!inInfo->isValid){
    EXCEPTION("Could not validate required input result \"" << inInfo->resultName << "\" from upstream filters.");
  }
  //we require mesh result input
  if(!inInfo->isMeshResult){
    EXCEPTION("RBF interpolation requires input to be defined on mesh");
  }

  //got the upstream result validated?
  if(!inInfo->isValid){
    EXCEPTION("Problem in filter pipeline detected. Interpolator input result \"" <<  inInfo->resultName << "\" could not be provided.")
  }

  resultManager_->CopyResultData(upResIds[0],filterResIds[0]);
  //but now, we need to overwrite some things
  resultManager_->SetRegionNames(filterResIds[0],this->trgRegions_);

  //after this filter we have nodal values on different regions
  //on a different grid
  resultManager_->SetDefOn(filterResIds[0],ExtendedResultInfo::NODE);
  resultManager_->SetGrid(filterResIds[0],this->trgGrid_);
  resultManager_->SetMeshResult(filterResIds[0],true);

  //validate own result
  resultManager_->SetValid(filterResIds[0]);
}


}
