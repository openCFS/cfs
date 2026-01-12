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
 *       \author   kroppert, sschoder
 */
//================================================================================================


#include "RBFInterpolator.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
//#include "cfsdat/DatUtils/KNNSearch.hh"
#include <algorithm>
#include <boost/tuple/tuple.hpp>
#include <vector>
#include <tuple>

#ifdef USE_CGAL
#include <def_use_cgal.hh>
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

RBFInterpolator::RBFInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
:MeshFilter(numWorkers,config,resMan){

  this->filtStreamType_ = FIFO_FILTER;


  if(config->Has("IntSchemeRBF") == false){
	  EXCEPTION("Scheme tag required for RBF interpolation")
  }

  globalFactor_ = config->Get("IntSchemeRBF")->Get("globalFactor")->As<Double>();
  inDim_ = 0;
  p_ = config->Get("IntSchemeRBF")->Get("interpolationExponent")->As<UInt>();
  numNN_ = config->Get("IntSchemeRBF")->Get("numNeighbours")->As<UInt>();
  numNW_ = config->Get("IntSchemeRBF")->Get("numNeighbours_weight")->As<UInt>();

  //No Slip Boundary handling on target Grid
  noSlip_ = false;
  if(config->Has("noSlipWall")){noSlip_ = true;}
  useElemAsTarget_ = false;
  if(params_->Has("useElemAsTarget")){useElemAsTarget_ = params_->Get("useElemAsTarget")->As<bool>();}

  //Target Entity Nodes/Element centrooids
  useElemAsTarget_ = false;
  if(params_->Has("useElemAsTarget")){useElemAsTarget_ = params_->Get("useElemAsTarget")->As<bool>();}

}

RBFInterpolator::~RBFInterpolator(){

}

#ifdef USE_CGAL

bool RBFInterpolator::UpdateResults(std::set<uuids::uuid>& upResults) {
  /// this is the vector, which will be filled with the result
  Vector<Double>& returnVec = GetOwnResultVector<Double>(filterResIds[0]);
  Integer stepIndex = resultManager_->GetStepIndex(filterResIds[0]);

  // vector, containing the source data values
  Vector<Double>& inVec = GetUpstreamResultVector<Double>(upResIds[0], stepIndex);


  // RBF setup
  Matrix& matrix = matrices_[matrixIndex_];
  const UInt maxNumTrgEntities = matrix.numTargets;
  StdVector<CF::UInt>& targetSourceIndex = matrix.targetSourceIndex;
  StdVector<CF::UInt>& targetSource = matrix.targetSource;
  StdVector<CF::Double>& targetSourceFactor = matrix.targetSourceFactor;
  StdVector<CF::Double>& targetSourceFactor2 = matrix.targetSourceFactor2;
  StdVector< CF::Matrix<Double> >& targetRBFInv = matrix.targetRBFInvMat;

  // RBF Interpolation
  RBFInterpolation(returnVec, inVec, numEquPerEnt_, targetSource, targetSourceIndex, targetRBFInv, targetSourceFactor, targetSourceFactor2, maxNumTrgEntities);
  //Final scaling
  returnVec.ScalarMult(globalFactor_);

  return true;
}

/*
CF::UInt RBFInterpolator::CountUsedEntities(const StdVector<CF::UInt>& entities) {
  const CF::UInt size = entities.GetSize();
  CF::UInt numEntities = 0;
#pragma omp parallel for reduction(+:numEntities) num_threads(CFS_NUM_THREADS)
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
#pragma omp parallel for num_threads(CFS_NUM_THREADS)
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





CF::StdVector<RBFInterpolator*> RBFInterpolator::interpolators_;
CF::StdVector<RBFInterpolator::Matrix> RBFInterpolator::matrices_;


//TODO Merge PrepareCGAL() and PreparePATCH() ... a lot of c&p
void RBFInterpolator::PrepareCalculation(){

  if(params_->Get("IntSchemeRBF/useCGAL4RBF")->As<std::string>()=="true"){
    PrepareCGAL();
  }else{
    PreparePATCH();
  }
  std::cout << "\n\t\t Interpolation prepared!" << std::endl;
}

//TODO prepare entities
//void RBFInterpolator::PrepareEntities(){
//
//}


void RBFInterpolator::PreparePATCH(){
  std::cout << "\t ---> RBFInterpolator preparing for interpolation" << std::endl;



  std::cout << "\t\t 1/4 Obtaining source entities " << std::endl;
  uuids::uuid upRes = upResIds[0];
  inGrid_ = resultManager_->GetExtInfo(upRes)->ptGrid;
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  scrMap_ = resultManager_->GetEqnMap(upRes);
  numEquPerEnt_ = scrMap_->GetNumEqnPerEnt();
  bool inElems = inInfo->definedOn == ExtendedResultInfo::ELEMENT;

  if ( inElems ){
	  EXCEPTION("Element to Node interpolation using RBF currently not possible");
  }



  const CF::UInt maxNumSrcEntities = scrMap_->GetNumEntities();
  StdVector<CF::UInt> globSrcEntity;
  GetUsedMappedEntities(scrMap_, globSrcEntity, srcRegions_, inGrid_);
  const CF::UInt numSrcEntities = CountUsedEntities(globSrcEntity);
  if (numSrcEntities < numNN_) {
    numNN_ = numSrcEntities;
  }

  trgMap_ = resultManager_->GetEqnMap(filterResIds[0]);


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
               " source " << (inElems ? "elements" : "nodes") << " and "<< numTrgEntities << " target "
               << (useElemAsTarget_ ? "elements" : "nodes") << std::endl;

  std::cout << "\t\t 3/4 Creating search tree " << std::endl;
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
      indices.push_back(srcEnt);
      sCoord[srcEnt] = pCoord;
    }
  }

  std::cout<< "\t\t 4/5 Boundary handling if activated " << std::endl;

  if(noSlip_){
    std::string regionName = params_->Get("noSlipWall/name")->As<std::string>();
    trgGrid_->GetNodesByRegion(boundary_,trgGrid_->GetRegionId(regionName));//TODO this is the region name...
    //TODO loop over targets...
    std::cout << "\t\t\t Interpolator considers a no slip boundary at "
        << regionName << std::endl;
  }


  std::cout << "\t\t 5/5 Creating interpolation matrix " << std::endl;
  matrix.numTargets = maxNumTrgEntities;
  StdVector<CF::UInt>& targetSourceIndex = matrix.targetSourceIndex;
  StdVector<CF::UInt>& targetSource = matrix.targetSource;
  StdVector<CF::Double>& targetSourceFactor = matrix.targetSourceFactor;
  StdVector<CF::Double>& targetSourceFactor2 = matrix.targetSourceFactor2;
  StdVector< CF::Matrix<Double> >& targetInvMat = matrix.targetRBFInvMat;

  targetSourceIndex.Resize(maxNumTrgEntities+1);
  CF::UInt index = 0;

  targetInvMat.Resize(maxNumTrgEntities);

  //TODO Fix that dirty stuff!!
  targetSource.Resize(100 * maxNumTrgEntities);
  targetSourceFactor.Resize(100 * maxNumTrgEntities);
  targetSourceFactor.Init();
  targetSourceFactor2.Resize(100 * maxNumTrgEntities);
  targetSourceFactor2.Init();


  std::cout<< "\t\t  [0% ------------ 100%]" << std::endl;
  std::cout<< "\t\t  ["<< std::flush;

  StdVector<RegionIdType>  volRegions;
  inGrid_->GetVolRegionIds(volRegions);

  // Prepare Nodes -> Element Map for the Input Grid
  // This is used for the source neighbor search
  inGrid_->SetNodesToElemsMap();

  StdVector<Vector<Double> > globCoords(maxNumTrgEntities);
  StdVector<LocPoint> lps;
  StdVector<const Elem*> elems;

//TODO this only works for one region !!
  std::string regionName = params_->Get("regions/sourceRegions/region/name")->As<std::string>();
  shared_ptr<EntityList> actSDList =  inGrid_->GetEntityList(EntityList::ELEM_LIST, regionName);
  StdVector<shared_ptr<EntityList> > inEntities(1);
  inEntities[0] = actSDList;

  for(CF::UInt trgEnt = 0; trgEnt < maxNumTrgEntities; trgEnt++) { // Loop over the Target points
    CF::UInt globEntityNumber = globTrgEntity[trgEnt];
    if (globEntityNumber != UnusedEntityNumber) {
      if(useElemAsTarget_){
        trgGrid_->GetElemCentroid(globCoords[trgEnt], globEntityNumber,true);
      } else {
        trgGrid_->GetNodeCoordinate3D(globCoords[trgEnt], globEntityNumber);
//        //TODO if belongs to entity then this is zero Inv matrix
//        if(noSlip_){//TODO
//          std::string regionName = params_->Get("noSlipWall/name")->As<std::string>();
//          trgGrid_->GetNodesByRegion(boundary,trgGrid_->GetRegionId(regionName));//TODO this is the region name...
//        }
      }

    }
  }

  // Get the Elements on the input Grid, where the Target points are located at
  inGrid_->GetElemsAtGlobalCoords( globCoords, lps, elems, inEntities);


  for(CF::UInt trgEnt = 0; trgEnt < maxNumTrgEntities; trgEnt++) { // Loop over all target points

    if((trgEnt)%(int(maxNumTrgEntities/20)) == 0){
    std::cout<< "#"<< std::flush;
    }

    targetSourceIndex[trgEnt] = index;

    StdVector<CF::Double> srcDist;
    CF::UInt globEntityNumber = globTrgEntity[trgEnt];

    if (globEntityNumber != UnusedEntityNumber) {
        CF::Vector<Double> pCoord = globCoords[trgEnt];
        StdVector<UInt> listNt ;

        // Element of the input grid that belongs to the target point
        const Elem* curE = elems[trgEnt];
        StdVector<UInt> nodeList ;
        nodeList.Resize(1);
        StdVector<const Elem*> srcElements;
        StdVector<UInt> listN ;
        srcElements.Insert(0,curE);
        inGrid_->GetNodesOfElemList(listN,srcElements,false);
        StdVector<const Elem*> gotElements;
	StdVector<const Elem*> tmp;
        for(UInt bNode =0;bNode < listN.GetSize(); ++bNode){
          inGrid_->GetElemsNextToNode(tmp,listN[bNode]);
          //StdVector<const Elem*> const & tmp = inGrid_->GetElemsByNode(listN[bNode]);
          for(UInt bElemN =0;bElemN < tmp.GetSize(); ++bElemN){
            //if(!gotElements.Contains(inGrid_->GetElem(tmp[bElemN]->elemNum)))  gotElements.Push_back(inGrid_->GetElem(tmp[bElemN]->elemNum));
            if(!gotElements.Contains(tmp[bElemN]))  gotElements.Push_back(tmp[bElemN]);
          }
        }

        if(inInfo->definedOn == ExtendedResultInfo::ELEMENT){ // ELEMENT BASED SOURCE
          numNN_ = gotElements.GetSize();
          index += numNN_;
          listNt.Resize(numNN_);
          listNt.Init();
          srcDist.Resize(numNN_);
          srcDist.Init();
          for(UInt bNode =0;bNode < gotElements.GetSize(); ++bNode){
            listNt[bNode] = gotElements[bNode]->elemNum;
          }

        }else{ // NODE BASED SOURCE

          inGrid_->GetNodesOfElemList(listNt,gotElements,false);

          numNN_ = listNt.GetSize();
          index += numNN_;
          srcDist.Resize(numNN_);
          srcDist.Init();
        }

        if (numNN_ > 1) {
          const CF::UInt ITrgSrcIndex = targetSourceIndex[trgEnt];
          CF::UInt* srcIndices = &targetSource[ITrgSrcIndex];
          CF::Double* srcFactors = &targetSourceFactor[ITrgSrcIndex];
          CF::Double* srcFactors2 = &targetSourceFactor2[ITrgSrcIndex];

          StdVector< CF::Vector<CF::Double> > neighbourCoords;
          neighbourCoords.Resize(numNN_);
          CF::Vector<Double> srcCoord;
          CF::Double dmax = 0.0;

          for(UInt i =0;i < numNN_; ++i){
            srcIndices[i] = listNt[i]-1;
            CF::Double distance = DistanceEUCLID(pCoord,sCoord[srcIndices[i]]);

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
          //Check if this is a boundary node
          if(noSlip_&&boundary_.Find(globEntityNumber)>=0){
            targetInvMat[trgEnt] *= 0.0;
          }


          Double r_k;
          Double temp;
          Double r_Wk = dmax;
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

        }
      } else {
        if (numNN_ == 1) {
          targetSourceIndex[trgEnt] = UnusedEntityNumber;
        }
      }

  }
  targetSourceIndex[maxNumTrgEntities] = index;
}


void RBFInterpolator::PrepareCGAL(){

  std::cout << "\t ---> RBFInterpolator preparing for interpolation" << std::endl;

  std::cout << "\t\t 1/4 Obtaining source entities " << std::endl;
  uuids::uuid upRes = upResIds[0];
  inGrid_ = resultManager_->GetExtInfo(upRes)->ptGrid;
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  scrMap_ = resultManager_->GetEqnMap(upRes);
  numEquPerEnt_ = scrMap_->GetNumEqnPerEnt();
  bool inElems = inInfo->definedOn == ExtendedResultInfo::ELEMENT;

  const CF::UInt maxNumSrcEntities = scrMap_->GetNumEntities();
  StdVector<CF::UInt> globSrcEntity;
  GetUsedMappedEntities(scrMap_, globSrcEntity, srcRegions_, inGrid_);
  const CF::UInt numSrcEntities = CountUsedEntities(globSrcEntity);
  if (numSrcEntities < numNN_) {
    numNN_ = numSrcEntities;
  }

  trgMap_ = resultManager_->GetEqnMap(filterResIds[0]);


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
               " source " << (inElems ? "elements" : "nodes") << " and "<< numTrgEntities << " target "
               << (useElemAsTarget_ ? "elements" : "nodes") << std::endl;

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

  StdVector<Point_and_int> searchPoints;
  searchPoints.Reserve(points.size());
  size_t i = 0;
  for(const auto& pt : points)
    searchPoints.Push_back(std::make_tuple(pt, indices[i++]));
  Tree tree(searchPoints.begin(), searchPoints.end());

  std::cout<< "\t\t 4/5 Boundary handling if activated " << std::endl;

  if(noSlip_){
    std::string regionName = params_->Get("noSlipWall/name")->As<std::string>();
    trgGrid_->GetNodesByRegion(boundary_,trgGrid_->GetRegionId(regionName));//TODO this is the region name...
    //TODO loop over targets...
    std::cout << "\t\t\t Interpolator considers a no slip boundary at "
        << regionName << std::endl;
  }

  std::cout << "\t\t 5/5 Creating interpolation matrix " << std::endl;
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


#pragma omp parallel for num_threads(CFS_NUM_THREADS)
  for(CF::UInt trgEnt = 0; trgEnt < maxNumTrgEntities; trgEnt++) {
    StdVector<CF::Double> srcDist;
    srcDist.Resize(numNN_);
    srcDist.Init();
    CF::UInt globEntityNumber = globTrgEntity[trgEnt];
    if (globEntityNumber != UnusedEntityNumber) {
      CF::Vector<Double> pCoord;
      if(useElemAsTarget_){
        trgGrid_->GetElemCentroid(pCoord, globEntityNumber,true);
      } else {
        trgGrid_->GetNodeCoordinate3D(pCoord, globEntityNumber);
      }
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
          srcIndices[i] = std::get<1>(it->first);
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
        //Check if this is a boundary node
        if(noSlip_&&boundary_.Find(globEntityNumber)>=0){
          targetInvMat[trgEnt] *= 0.0;
        }


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
        targetSourceIndex[trgEnt] = std::get<1>(search.begin()->first);
      }
    } else {
      if (numNN_ == 1) {
        targetSourceIndex[trgEnt] = UnusedEntityNumber;
      }
    }
  }


}

Double RBFInterpolator::DistanceEUCLID(CF::Vector<Double> p1, CF::Vector<Double> p2){
 Double dist = 0.0;
 if(p1.GetSize()!=p1.GetSize()){
   EXCEPTION("No distance can be computed in RBF interpolator!")
 }
 for( UInt i = 0; i<p1.GetSize(); ++i){
   dist += (p1[i]-p2[i])*(p1[i]-p2[i]);
 }
 return sqrt(dist);
}

ResultIdList RBFInterpolator::SetUpstreamResults(){
  return SetDefaultUpstreamResults();
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
  if(useElemAsTarget_){
    resultManager_->SetDefOn(filterResIds[0],ExtendedResultInfo::ELEMENT);
  } else {
    resultManager_->SetDefOn(filterResIds[0],ExtendedResultInfo::NODE);
  }


  resultManager_->SetGrid(filterResIds[0],this->trgGrid_);
  resultManager_->SetMeshResult(filterResIds[0],true);

  //validate own result
  resultManager_->SetValid(filterResIds[0]);
}

#else

ResultIdList RBFInterpolator::SetUpstreamResults(){
  return SetDefaultUpstreamResults();
}
bool RBFInterpolator::UpdateResults(std::set<uuids::uuid>& upResults) {
	return true;
}
void RBFInterpolator::PrepareCalculation(){}
void RBFInterpolator::AdaptFilterResults(){}
#endif

}
