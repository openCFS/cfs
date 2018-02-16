// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     NearestNeighbourInterpolator.cc
 *       \brief    <Description>
 *
 *       \date     Aug 8, 2016
 *       \author   kroppert
 */
//================================================================================================



#include "NearestNeighbourInterpolator.hh"
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



// check if Intel MKL is available
#include <def_use_blas.hh>
#include <boost/date_time/posix_time/posix_time_types.hpp>
#ifdef USE_MKL
  #include <mkl.h>
#endif

namespace CFSDat{

NearestNeighbourInterpolator::NearestNeighbourInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
                     :MeshFilter(numWorkers,config,resMan){

#ifndef USE_CGAL
    EXCEPTION("CoefFunctionScatteredData needs to be compiled with USE_CGAL=ON!");
#endif

  this->filtStreamType_ = FIFO_FILTER;

  p_ = config->Get("scheme")->Get("interpolationExponent")->As<UInt>();
  numNeighbors_ = config->Get("scheme")->Get("numNeighbours")->As<UInt>();


  mCheck_ = false;
}

NearestNeighbourInterpolator::~NearestNeighbourInterpolator(){

}

bool NearestNeighbourInterpolator::UpdateResults(std::set<uuids::uuid>& upResults){
  /// this is the vector, which will be filled with the derivative result
  Vector<Double>& returnVec = GetOwnResultVector<Double>(filterResIds[0]);
  Double aTF = resultManager_->GetStepValue(filterResIds[0]);

  // vector, containing the source data values
  Vector<Double>& inVec = GetUpstreamResultVector<Double>(upResIds[0], aTF);
  
  // TODO matrix interpolation here
  Matrix& matrix = matrices_[matrixIndex_];
  CF::UInt& maxNumTrgEntities = matrix.numTargets;
  StdVector<CF::UInt>& targetSourceIndex = matrix.targetSourceIndex;
  StdVector<CF::UInt>& targetSource = matrix.targetSource;
  StdVector<CF::Double>& targetSourceFactor = matrix.targetSourceFactor;
  
  NearestNeighbourInterpolation(returnVec, inVec, numEquPerEnt_, targetSource, targetSourceIndex, numNeighbors_, targetSourceFactor, maxNumTrgEntities);

  return true;
}







void NearestNeighbourInterpolator::CheckFilterResults(Vector<Double>& origVec,
                                                      Vector<Double>& newVec){
/*
#ifndef USE_MKL
  EXCEPTION("Compile with Intel MKL for mesh-check!!");
#endif
using namespace boost::posix_time;
using namespace boost::gregorian;

std::cout<<"origVec.GetSize()"<<origVec.GetSize()<<std::endl;
std::cout<<"newVec.GetSize()"<<newVec.GetSize()<<std::endl;


ptime time_startMINE(microsec_clock::local_time());

Double xi;
Double yi;
Double sumXiYi = 0.0;
Double sumXi = 0.0;
Double sumYi = 0.0;
Double sumXi_squared = 0.0;
Double sumYi_squared = 0.0;
//#pragma omp parallel shared(origVec, newVec) num_threads(NUM_CFS_THREADS)
//{
//#pragma omp for
for (UInt i = 0; i < origVec.GetSize(); ++i){
  xi = origVec[i];
  yi = newVec[i];
  sumXiYi += xi * yi;

  sumXi += xi;
  sumYi += yi;

  sumXi_squared += xi*xi;
  sumYi_squared += yi*yi;
}
//}
sumXiYi = sumXiYi * origVec.GetSize();

Double nom = sumXiYi - sumXi * sumYi;

Double denom = std::sqrt(origVec.GetSize()*sumXi_squared - sumXi*sumXi) * std::sqrt(origVec.GetSize()*sumYi_squared - sumYi*sumYi);

Double result = nom / denom;

ptime time_endMINE(microsec_clock::local_time());
time_duration durationMINE(time_endMINE - time_startMINE);
std::cout<<"MINE" << durationMINE <<std::endl;
std::cout<<"correlation=" << result <<std::endl;



ptime time_startMKL(microsec_clock::local_time());

MKL_INT xshape; //declare some variables
xshape = origVec.GetSize();
int i,j;


// now the summary statistics part
VSLSSTaskPtr task2;
MKL_INT dim = 2;
MKL_INT n = xshape;
MKL_INT X_storage;
MKL_INT cov_storage;
MKL_INT cor_storage;
double X[xshape][dim];
double cov[dim][dim], cor[dim][dim];
double mean[dim];



X_storage   = VSL_SS_MATRIX_STORAGE_COLS;
cov_storage = VSL_SS_MATRIX_STORAGE_FULL;
cor_storage = VSL_SS_MATRIX_STORAGE_FULL;


for(i = 0; i < xshape; i++){
    X[i][0] = origVec.GetDoubleEntry(i);
    X[i][1] = newVec.GetDoubleEntry(i);
}

for(i = 0; i < dim; i++)
{
    mean[i] = 0.0;
    for(j = 0; j < dim; j++)
    {
        cov[i][j] = 0;
        cor[i][j] = 0;
    }
}


vsldSSNewTask( &task2, &dim, &n, &X_storage, *X, 0, 0 );



vsldSSEditCovCor( task2, mean, *cov, &cov_storage, *cor, &cor_storage );



vsldSSCompute( task2, VSL_SS_COV|VSL_SS_COR, VSL_SS_METHOD_FAST );


vslSSDeleteTask( &task2 );

MKL_Free_Buffers();

ptime time_endMKL(microsec_clock::local_time());
time_duration durationMKL(time_endMKL - time_startMKL);
std::cout<<"MKL" << durationMKL <<std::endl;
std::cout<<"correlation=" << cor[1][0] <<std::endl;
*/


//std::cout<<"correlation="<<result<<std::endl;

}


/*
bool NearestNeighbourInterpolator::CreatesEqualMatrix(NearestNeighbourInterpolator* otherInterpolator) {
  return (otherInterpolator->numNeighbors_ == numNeighbors_) && (otherInterpolator->p_ == p_)
    && (otherInterpolator->inGrid_ == inGrid_) && (otherInterpolator->trgGrid_ == trgGrid_)
    && (otherInterpolator->srcRegions_ == srcRegions_) && (otherInterpolator->trgRegions_ == trgRegions_)
    && (otherInterpolator->scrMap_->Equals(*scrMap_)) && (otherInterpolator->trgMap_->Equals(*trgMap_));
}


CF::UInt NearestNeighbourInterpolator::CountUsedEntities(const StdVector<CF::UInt>& entities) {
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

void NearestNeighbourInterpolator::GetUsedMappedEntities(const str1::shared_ptr<EqnMapSimple>& map,
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





CF::StdVector<NearestNeighbourInterpolator*> NearestNeighbourInterpolator::interpolators_;
CF::StdVector<NearestNeighbourInterpolator::Matrix> NearestNeighbourInterpolator::matrices_;


void NearestNeighbourInterpolator::PrepareCalculation(){
  std::cout << "\t ---> NearesNeighbourInterpolator preparing for interpolation" << std::endl;
  
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
  if (numSrcEntities < numNeighbors_) {
    numNeighbors_ = numSrcEntities;
  }
  
  trgMap_ = resultManager_->GetEqnMap(filterResIds[0]);
  /*
  // checking, if interpolation matrix need to be created
  for (UInt i = 0; i < interpolators_.GetSize(); i++) {
    
    if (CreatesEqualMatrix(interpolators_[i])) {
      matrixIndex_ = i;
      std::cout << "\t\t Interpolation matrix already prepared! Skipping further steps!" << std::endl;
      return;
    }
  }
  */
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
    }
  }
  Tree tree(boost::make_zip_iterator(boost::make_tuple( points.begin(),indices.begin() )),
            boost::make_zip_iterator(boost::make_tuple( points.end(),indices.end() ) )
  );
  
  std::cout << "\t\t 4/4 Creating interpolation matrix " << std::endl;
  matrix.numTargets = maxNumTrgEntities;
  StdVector<CF::UInt>& targetSourceIndex = matrix.targetSourceIndex;
  StdVector<CF::UInt>& targetSource = matrix.targetSource;
  StdVector<CF::Double>& targetSourceFactor = matrix.targetSourceFactor;
  
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
  
  // creating target -> source factors and filling source indices
  if (numNeighbors_ > 1) {
    targetSource.Resize(index);
    targetSourceFactor.Resize(index);
  }
  #pragma omp parallel for num_threads(NUM_CFS_THREADS)
  for(CF::UInt trgEnt = 0; trgEnt < maxNumTrgEntities; trgEnt++) {
    CF::UInt globEntityNumber = globTrgEntity[trgEnt];
    if (globEntityNumber != UnusedEntityNumber) {
      CF::Vector<Double> pCoord;
      trgGrid_->GetNodeCoordinate3D(pCoord, globEntityNumber);
      Point_3 query(pCoord[0],pCoord[1],pCoord[2]);
      Distance tr_dist;
      K_neighbor_search search(tree, query, numNeighbors_);
      if (numNeighbors_ > 1) {
        const CF::UInt ITrgSrcIndex = targetSourceIndex[trgEnt];
        CF::UInt* srcIndices = &targetSource[ITrgSrcIndex];
        CF::Double* srcFactors = &targetSourceFactor[ITrgSrcIndex];
        CF::UInt i = 0;
        CF::Double dmax = 0.0;
        for(K_neighbor_search::iterator it = search.begin(); it != search.end(); it++, i++) {
          srcIndices[i] = boost::get<1>(it->first);
          CF::Double distance = tr_dist.inverse_of_transformed_distance(it->second);
          srcFactors[i] = distance;
          if (distance > dmax) {
            dmax = distance;
          }
        }
//for(UInt l = 0; l<10;++l){
//std::cout<<targetSourceFactor[l]<<std::endl;
//}
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
          if(srcFactors[i] == 0){
            w = 1.0;
          } else {
            w = std::pow((R-srcFactors[i])/(R*srcFactors[i]), p_);
          }
          srcFactors[i] = w;
          weights += w;
        }
        for (i = 0; i < numNeighbors_; i++) {
          srcFactors[i] /= weights;
        }
      } else {
        targetSourceIndex[trgEnt] = boost::get<1>(search.begin()->first);
      }
    } else {
      if (numNeighbors_ == 1) {
        targetSourceIndex[trgEnt] = UnusedEntityNumber;
      }
    }
  }
  
  std::cout << "\t\t Interpolation prepared!" << std::endl;
}


ResultIdList NearestNeighbourInterpolator::SetUpstreamResults(){
  return SetDefaultUpstreamResults();
}

void NearestNeighbourInterpolator::AdaptFilterResults(){
  //some checks
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  if(!inInfo->isValid){
    EXCEPTION("Could not validate required input result \"" << inInfo->resultName << "\" from upstream filters.");
  }
  //we require mesh result input
  if(!inInfo->isMeshResult){
    EXCEPTION("NearesNeighbour interpolation required input to be defined on mesh");
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
