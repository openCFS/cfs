// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CurlDifferentiator.cc
 *       \brief    <Description>
 *
 *       \date     Oct 6, 2016
 *       \author   kroppert
 */
//================================================================================================


#include <Filters/Derivatives/CurlDifferentiator.hh>
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

CurlDifferentiator::CurlDifferentiator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
:MeshFilter(numWorkers,config,resMan){

#ifndef USE_CGAL
    EXCEPTION("CoefFunctionScatteredData needs to be compiled with USE_CGAL=ON!");
#endif


  this->filtStreamType_ = FIFO_FILTER;

  numNeighbors_ = (trgGrid_->GetDim() == 2) ? 4 : 8;
}

CurlDifferentiator::~CurlDifferentiator(){

}

bool CurlDifferentiator::Run(){
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
  // now we deactivate our own result and activate the others
  resultManager_->ActivateResult(upResIds[0]);

  //now we call for upstream data in each source
  CF::StdVector< str1::shared_ptr<BaseFilter> >::iterator srcIter =  sources_.Begin();
  for(; srcIter != sources_.End() ; srcIter++){
    // should we check here anything for success?
    (*srcIter)->Run();
  }


  CF::StdVector<UInt> eqnNums;

  /// this is the vector, which will be filled with the derivative result
  Vector<Double>& returnVec = resultManager_->GetResultVector<Double>(filterResIds[0],eqnNums);
  returnVec.Init();

  // vector, containing the source data values
  Vector<Double>& inVec = resultManager_->GetResultVector<Double>(upResIds[0],eqnNums);



  Matrix& matrix = matrices_[matrixIndex_];
  const UInt maxNumTrgEntities = matrix.numTargets;
  StdVector<CF::UInt>& targetSourceIndex = matrix.targetSourceIndex;
  StdVector<CF::UInt>& targetSource = matrix.targetSource;
  StdVector< CF::Matrix<CF::Double> >& targetSourceFactor = matrix.targetSourceFactor;

  UInt gridDim = inGrid_->GetDim();

  CalcCurl(returnVec, inVec, numEquPerEnt_, targetSource, targetSourceIndex, numNeighbors_, targetSourceFactor, maxNumTrgEntities, gridDim);


  resultManager_->ActivateResult(filterResIds[0]);

  //now deactivate own upstream results
  for(UInt aRes=0;aRes<upResIds.GetSize();aRes++){
    resultManager_->DeactivateResult(upResIds[aRes]);
  }

  return true;
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


CF::StdVector<CurlDifferentiator*> CurlDifferentiator::differentiators_;
CF::StdVector<CurlDifferentiator::Matrix> CurlDifferentiator::matrices_;


void CurlDifferentiator::PrepareCalculation(){
  std::cout << "\t ---> CurlDifferentiator preparing for interpolation" << std::endl;


  std::cout << "\t\t 1/4 Obtaining source entities " << std::endl;
  uuids::uuid upRes = upResIds[0];
  inGrid_ = resultManager_->GetExtInfo(upRes)->ptGrid;
  scrMap_ = resultManager_->GetResultAdapter(upRes)->mapping;
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  numEquPerEnt_ = scrMap_->GetNumEqnPerEnt();
  if(numEquPerEnt_ == 1){
    EXCEPTION("Curl of a scalar is not defined in this context!");
  }

  bool inElems = inInfo->definedOn == ExtendedResultInfo::ELEMENT;
  if(inElems){
    // in this case, the nn search also finds the target-point...distance 0,
    // so add one more neighbour
    numNeighbors_ = numNeighbors_ + 1;
  }

  const CF::UInt maxNumSrcEntities = scrMap_->GetNumEntities();
  StdVector<CF::UInt> globSrcEntity;
  GetUsedMappedEntities(scrMap_, globSrcEntity, srcRegions_, inGrid_);
  const CF::UInt numSrcEntities = CountUsedEntities(globSrcEntity);
  if (numSrcEntities < numNeighbors_) {
    numNeighbors_ = numSrcEntities;
  }

  trgMap_ = resultManager_->GetResultAdapter(filterResIds[0])->mapping;


  differentiators_.Push_back(this);
  matrixIndex_ = matrices_.GetSize();
  matrices_.Resize(matrixIndex_ + 1);
  Matrix& matrix = matrices_[matrixIndex_];

  std::cout << "\t\t 2/4 Obtaining target entities " << std::endl;

  const CF::UInt maxNumTrgEntities = trgMap_->GetNumEntities();
  StdVector<CF::UInt> globTrgEntity;
  GetUsedMappedEntities(trgMap_, globTrgEntity, trgRegions_, trgGrid_);
  const CF::UInt numTrgEntities = CountUsedEntities(globTrgEntity);

  std::cout << "\t\t\t Differentiator is dealing with " << numSrcEntities <<
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
  StdVector< CF::Matrix<CF::Double> >& targetSourceFactor = matrix.targetSourceFactor;

  // creating target -> source indices
  targetSourceIndex.Resize(maxNumTrgEntities + 1);
  CF::UInt index = 0;
  for(CF::UInt trgEnt = 0; trgEnt < maxNumTrgEntities; trgEnt++) {
    targetSourceIndex[trgEnt] = index;
    if (globTrgEntity[trgEnt] != UnusedEntityNumber) {
      index += numNeighbors_;
    }
  }
  targetSourceIndex[maxNumTrgEntities] = index;
  targetSourceFactor.Resize(maxNumTrgEntities);


  // creating target -> source factors and filling source indices
  if (numNeighbors_ > 1) {
    targetSource.Resize(index);
  }


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
      K_neighbor_search search(tree, query, numNeighbors_);
      if (numNeighbors_ > 1) {
        const CF::UInt ITrgSrcIndex = targetSourceIndex[trgEnt];
        CF::UInt* srcIndices = &targetSource[ITrgSrcIndex];
        //CF::Matrix<CF::Double>* srcFactors = &targetSourceFactor[ITrgSrcIndex];

        StdVector< CF::Vector<CF::Double> > neighbourCoords;
        neighbourCoords.Resize(numNeighbors_);
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
        if(alpha == 0){
          std::cout<<dmax*1.5e+03<<std::endl;
          EXCEPTION("CurlDifferentiator::PrepareCalculation  alpha = 0!!!");
        }
        CalcLocCurl(targetSourceFactor[trgEnt], pCoord, srcDist, neighbourCoords, alpha, numNeighbors_, numEquPerEnt_, inGrid_);
      } else {
        targetSourceIndex[trgEnt] = boost::get<1>(search.begin()->first);
      }
    } else {
      if (numNeighbors_ == 1) {
        targetSourceIndex[trgEnt] = UnusedEntityNumber;
      }
    }
  }

  std::cout << "\t\t Differentiation prepared!" << std::endl;
}

ResultIdList CurlDifferentiator::SetUpstreamResults(){
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

void CurlDifferentiator::AdaptFilterResults(){
  //some checks
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  if(!inInfo->isValid){
    EXCEPTION("Could not validate required input result \"" << inInfo->resultName << "\" from upstream filters.");
  }
  //we require mesh result input
  if(!inInfo->isMeshResult){
    EXCEPTION("CurlDifferentiator requires input to be defined on mesh");
  }

  //got the upstream result validated?
  if(!inInfo->isValid){
    EXCEPTION("Problem in filter pipeline detected. Differentiator input result \"" <<  inInfo->resultName << "\" could not be provided.")
  }

  //input must be node-based
  if (resultManager_->GetDefOn(upResIds[0]) == ExtendedResultInfo::ELEMENT){
    //case of centroid-values
    std::cout<<("============================================================ \n"
                "Input of CurlDifferentiator-filter is defined on elements!! \n"
                "This works but it's not very accurate, especially at the boundary \n"
                "You better interpolate the element-values to nodes (e.g. Cell2Node) \n"
                "and differentiate afterwards\n"
                "============================================================")<<std::endl;
  }


  resultManager_->CopyResultData(upResIds[0],filterResIds[0]);
  //but now, we need to overwrite some things
  resultManager_->SetRegionNames(filterResIds[0],this->trgRegions_);

  resultManager_->SetEntryType(filterResIds[0],ExtendedResultInfo::VECTOR);

  // if we have a 2D-mesh, the output will be a scalar value (per definition, then all
  // differentiation- and aeroacoustic-filters are conistent)
  CF::StdVector<std::string> dofnames;
  if(trgGrid_->GetDim() == 2){
    dofnames.Push_back("x");
    dofnames.Push_back("y");
  }else{
    dofnames.Push_back("x");
    dofnames.Push_back("y");
    dofnames.Push_back("z");
  }
  resultManager_->SetDofNames(filterResIds[0],dofnames);

  //after this filter we have element values on different regions
  //on a different grid
  resultManager_->SetDefOn(filterResIds[0],ExtendedResultInfo::ELEMENT);
  resultManager_->SetGrid(filterResIds[0],this->trgGrid_);
  resultManager_->SetMeshResult(filterResIds[0],true);

  //validate own result
  resultManager_->SetValid(filterResIds[0]);
}


}

