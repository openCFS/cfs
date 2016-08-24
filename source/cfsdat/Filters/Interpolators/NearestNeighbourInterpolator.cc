// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     NearesNeighbourInterpolator.hh
 *       \brief    <Description>
 *
 *       \date     Aug 8, 2016
 *       \author   kroppert
 */
//================================================================================================


#include "NearestNeighbourInterpolator.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"

#include <algorithm>
#include <vector>

namespace CFSDat{

NearestNeighbourInterpolator::NearestNeighbourInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
                     :MeshBasedInterpolator(numWorkers,config,resMan){

  this->filtSteamType_ = FIFO_FILTER;
  inDim_ = 0;

  //check if the nearest neighbour search scheme has been defined in the xml-file
  if(config->Has("scheme") == false){
    EXCEPTION("============================================================ \n"
              "No <sheme> for nearest neighbour search was set in xml-file! \n"
              "============================================================");
  }
  //I'm sure, there is a much more elegant way to get the "name" of an element in the xml-scheme
  //but I couldn't manage to do it, so the following is kind of a workaround
  std::string alg = config->Get("scheme")->Get("searchAlgorithm")->As<std::string>();
  if(alg == "CGAL"){
    knnLib_ = 0;
  }else if(alg == "FLANN"){
    knnLib_ = 1;
  }else{
    EXCEPTION("no search algorithm selected!");
  }
  p_ = config->Get("scheme")->Get("interpolationExponent")->As<UInt>();
  searchRadius_ = config->Get("scheme")->Get("searchRadius")->As<Double>();
  numNeighbors_ = config->Get("scheme")->Get("numNeighbours")->As<UInt>();
}

NearestNeighbourInterpolator::~NearestNeighbourInterpolator(){

}

bool NearestNeighbourInterpolator::Run(){
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

  // this is the vector, containing the source data
  Vector<Double>& inVec = resultManager_->GetResultVector<Double>(upResIds[0],eqnNums);

  // vector containing the source values (scattered data values)
  CF::StdVector< CF::Vector<Double> >  scatteredData;
  scatteredData.Resize(sourceCoords_.GetSize());

  // Checking if input is scalar- or vector-type. This needs to be done
  // because of the input of CGAL and FLANN
  // Maybe there is a more efficient way to deal with this issue?!
  if(inVec.GetSize() == sourceCoords_.GetSize()){
    //this is the case of scalar scattered values
    inDim_=0;
    for(UInt i = 0; i < sourceCoords_.GetSize(); ++i){
        scatteredData[i].Resize(1);
        scatteredData[i][0] = inVec[i];
      }
  }
  if(inVec.GetSize() == 2 * sourceCoords_.GetSize()){
    //case of two dimensional vector
    inDim_=1;
    for(UInt i = 0; i < sourceCoords_.GetSize() ; ++i){
    scatteredData[i].Resize(2);
    scatteredData[i][0] = inVec[2 * i]; // x-component
    scatteredData[i][1] = inVec[2 * i + 1]; // y-component
    }
  }
  if(inVec.GetSize() == 3 * sourceCoords_.GetSize()){
    // case of three dimensional vector
    inDim_=2;
    for(UInt i = 0; i < sourceCoords_.GetSize(); ++i){
    scatteredData[i].Resize(3);
    scatteredData[i][0] = inVec[3 * i]; // x-component
    scatteredData[i][1] = inVec[3 * i + 1]; // y-component
    scatteredData[i][2] = inVec[3 * i + 2]; // z-component
    }
  }

  //now we can bring the scattered coordinates and data values into
  //the correct form for the nearest neighbour-search witch CGAL or FLANN
  ReadScatteredData(sourceCoords_, scatteredData);

  CF::StdVector< Vector<Double> > neighbors;
  CF::StdVector< Double > l2dists;
  CF::StdVector< Vector<Double> > vectors;
  Vector<Double> vec;

  //definition of target shape functions and corresponding node-equations
  CF::Vector<Double> shFnc;
  CF::StdVector<UInt> eqns;
  CF::shared_ptr<ElemShapeMap> eShape;
  str1::shared_ptr<EqnMapSimple> downMap = resultManager_->GetResultAdpter(filterResIds[0])->mapping;

  CF::StdVector<UInt> nodeCheck;
  nodeCheck.Resize(trgGrid_->GetNumNodes(ALL_REGIONS));
  nodeCheck.Init();
  UInt nodeIter = 0;
  bool nodeMatch;
  CF::StdVector<UInt> nodeList;
  CF::StdVector<CF::Elem*> elems;
  trgGrid_->GetVolElems(elems, ALL_REGIONS);
  trgGrid_->GetNodesOfElemList(nodeList, elems);

  // loop over all elements and over every node of each element
  for(UInt i=0;i < trgGrid_->GetNumVolElems();++i){
    InpolationStruct& aStru = interpolData_[i];
    const Elem* curE = trgGrid_->GetElem(elems[i]->elemNum);
    eShape = trgGrid_->GetElemShapeMap(curE,true);
    //eConn gives us the node numbers of element curE
    const CF::StdVector<UInt>& eConn = curE->connect;

    //get the global coordinates of the nodes of element curE
    CF::Matrix<Double> globalCoords;
    trgGrid_->GetElemNodesCoord(globalCoords, eConn, true);

    FeH1 * myElem = dynamic_cast<FeH1*>(eShape->GetBaseFE());
    shFnc.Resize(eConn.GetSize());
    shFnc.Init();

    myElem->GetShFnc(shFnc,aStru.localCoords,curE);

    //loop over every node of element curE and perform interpolation
    //BUT only if the interpolation for this certain node has not been
    //carried out before. Therefore we use a std::search in which we
    //search, if the current node has been used before
    for(UInt aNode =0;aNode < eConn.GetSize(); ++aNode){
      //extract the global point of node aNode
      CF::Vector<Double> globPoint;
      globalCoords.GetCol(globPoint, aNode);

      //that is the mentioned search, to find out if the node has already been used
      nodeMatch = nodeCheck.Contains(eConn[aNode]);
      if(nodeMatch == false){
        //add aNode to the "already-computed-list"
        nodeCheck[nodeIter]=eConn[aNode];

          // at that point we can start obtaining the nearest neighbours for every aNode
          // using CGAL or FLANN
          switch(knnLib_)
          {
          case 0: //CGAL
        #ifdef USE_CGAL
            KNNSearch_CGAL(globPoint, neighbors, l2dists, vectors);
        #else
            EXCEPTION("CoefFunctionScatteredData needs to be compiled with USE_CGAL=ON!");
        #endif
            break;

          case 1: //FLANN
        #ifdef USE_FLANN
            KNNSearch_FLANN(globPoint, neighbors, l2dists, vectors, scatteredData);

        #else
            EXCEPTION("CoefFunctionScatteredData needs to be compiled with USE_FLANN=ON!");
        #endif
            break;
          default:
            break;
          }
          vec.Resize(inDim_+1);
          Double dmin = l2dists[0];
          Double dmax = dmin;
          StdVector< Double >::iterator it, end;
          it = l2dists.Begin();
          end = l2dists.End();
          for( ; it != end; ++it) {
            Double dist = (*it);
            dmin = dmin < dist ? dmin : dist;
            dmax = dmax > dist ? dmax : dist;
          }

            // Apply Shepard interpolation cf. Numerical Recipes 3rd ed. p. 143ff.
            // or http://www.ems-i.com/gmshelp/Interpolation/Interpolation_Schemes \
            // /Inverse_Distance_Weighted/Shepards_Method.htm
            Vector<Double> sum(inDim_+1);
            Double weights = 0.0;
            // The point which is farthest away, should at least have a non-zero
            // weight of 0.01. If we would choose R = dmax, it would not contribute
            // at all.
            Double R = 1.01 * dmax;

            // report the N nearest neighbors and their distance
            // This should sort all N points by increasing distance from origin
            it = l2dists.Begin();
            for(UInt j=0; it != end; ++it, j++) {
              Double d = *it;
              Double w;
              if(d == 0){
               w = 1.0;
              }else{
               w = std::pow((R-d)/(R*d), p_);
              }
              weights += w;

              for(UInt dof=0; dof < inDim_ + 1; dof++)
              {
                sum[dof] += vectors[j][dof] * w;
              }
            }

            for(UInt dof=0; dof < inDim_ + 1; dof++)
            {
              vec[dof] = sum[dof] / weights;
            }

            //get the equation map for the nodes in eConn, in order to insert the
            //interpolation in the correct position in the result vector
            downMap->GetEquation(eqns,eConn[aNode],ExtendedResultInfo::NODE);

            //if scalar input values of scattered data->inDim_=0
            //if it is a two-dimensional vector->inDim_=1 and inDim_=2 for 3d-vector
            for(UInt aDof = 0; aDof < inDim_+1; aDof++){
              returnVec[eqns[aDof]] = vec[aDof];
            }
            //iterator for the already computed nodes
            nodeIter++;

        }// if nodeMatch == false
    }// for aNode
  }// for element

  resultManager_->ActivateResult(filterResIds[0]);

  //now deactivate own upstream results
  for(UInt aRes=0;aRes<upResIds.GetSize();aRes++){
    resultManager_->DeactivateResult(upResIds[aRes]);
  }

  return true;
}


void NearestNeighbourInterpolator::ReadScatteredData(CF::StdVector< CF::Vector<Double> > sourceCoords, CF::StdVector< CF::Vector<Double> > scatteredData)
{

  UInt n = sourceCoords_.GetSize();
  switch(knnLib_)
  {
  case 0: //CGAL
#ifdef USE_CGAL
    {
    std::vector<CGAL::Point> points;
    points.resize(n);
    for(UInt i=0; i<n; i++)
    {
      if(trgGrid_->GetDim() == 2){
        if(inDim_ == 0){
          points[i] = (CGAL::Point(sourceCoords_[i][0],
                                   sourceCoords_[i][1],
                                   0.0,
                                   scatteredData[i][0],
                                   0.0,
                                   0.0));
          }else{
            points[i] = (CGAL::Point(sourceCoords_[i][0],
                                     sourceCoords_[i][1],
                                     0.0,
                                     scatteredData[i][0],
                                     scatteredData[i][1],
                                     0.0));
            }
      }
      if(trgGrid_->GetDim() == 3){
        if(inDim_ == 0){
          points[i] = (CGAL::Point(sourceCoords_[i][0],
                                   sourceCoords_[i][1],
                                   sourceCoords_[i][2],
                                   scatteredData[i][0],
                                   0.0,
                                   0.0));
        }else{
          points[i] = (CGAL::Point(sourceCoords_[i][0],
                                   sourceCoords_[i][1],
                                   sourceCoords_[i][2],
                                   scatteredData[i][0],
                                   scatteredData[i][1],
                                   scatteredData[i][2]));
      }
    }
    searchTree_.reset(new Tree(points.begin(), points.end()));
    }
  }
#else
    EXCEPTION("CGAL not supported! Compile with USE_CGAL=ON.");
#endif
    break;

  case 1: //FLANN
#ifdef USE_FLANN
    {
        dataset_.reset(new flann::Matrix<Double>(new Double[n*3], n, 3));
        Double *dPtr = dataset_->ptr();
        if(trgGrid_->GetDim() == 3){
        for(UInt i=0; i<n; i++)
        {
          UInt idx = i*3;
          dPtr[idx+0] = sourceCoords_[i][0];
          dPtr[idx+1] = sourceCoords_[i][1];
          dPtr[idx+2] = sourceCoords_[i][2];
        }
      }else{
        for(UInt i=0; i<n; i++)
        {
          UInt idx = i*3;
          dPtr[idx+0] = sourceCoords_[i][0];
          dPtr[idx+1] = sourceCoords_[i][1];
          dPtr[idx+2] = 0.0;
        }
      }

    // construct an randomized kd-tree index using a single kd-tree
    index_.reset(new flann::Index<flann::L2<Double> >(*dataset_.get(),
                                                      flann::KDTreeSingleIndexParams(12)));
    index_->buildIndex();
    }
#else
    EXCEPTION("FLANN not supported! Compile with USE_FLANN=ON.")
#endif

    break;

  default:
    EXCEPTION("Unknown k-nearest neighbors library '" << knnLib_)
    break;
  }

  if(n < numNeighbors_)
  {
    numNeighbors_ = n;
  }
}




void NearestNeighbourInterpolator::PrepareInterpolation(){
  //1. Get get the source coordinates and the values, defined on those coordinates (Source...Src)
  //2. Get the target coordinates (trg)
  //3. Store for each trg element local Coordinates, elem number, volume, ...
  //4. Small clean-up


  std::cout << "\t ---> NearesNeighbourInterpolator preparing for interpolation" << std::endl;

  //in this filter we only have one upstream result
  uuids::uuid upRes = upResIds[0];

  Grid* inGrid   = resultManager_->GetExtInfo(upRes)->ptGrid;

  //lets declare some variables and estimate the memory
  std::vector<UInt> allSrcElems;

  //loop over source regions and add element numbers to vector
  std::set<std::string>::iterator sRegIter = srcRegions_.begin();
  for(;sRegIter != srcRegions_.end();++sRegIter){
    StdVector<UInt> curElems;
    inGrid->GetElemNumsByName(curElems,*sRegIter);
    allSrcElems.insert(allSrcElems.end(),curElems.Begin(),curElems.End());
  }

  std::cout << "\t\t\t Interpolator is dealing with " << allSrcElems.size() <<
               " source values" << std::endl;

  std::cout << "\t\t 1/5 Obtaining source coordinates " << std::endl;
  StdVector<shared_ptr<EntityList> > lists;
  std::set<std::string>::const_iterator destRegIt = this->trgRegions_.begin();
  for(; destRegIt != this->trgRegions_.end(); ++destRegIt ) {
    RegionIdType aReg = trgGrid_->GetRegion().Parse(*destRegIt);
    shared_ptr<ElemList> newList(new ElemList(trgGrid_));
    newList->SetRegion(aReg);
    lists.Push_back(newList);
  }

  //here we get some input information, because we want to find out, if the input-values
  //are defined on element-centroids (e.g. from CFD) or on nodes
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  if(inInfo->definedOn == ExtendedResultInfo::ELEMENT){
    //case of centroid-values
    sourceCoords_.Resize(allSrcElems.size());
    for(UInt i=0;i<allSrcElems.size();++i){
      CF::Vector<Double> cCoord;
      inGrid->GetElemCentroid(cCoord,allSrcElems[i],true);
      if(trgGrid_->GetDim() == 2){
        sourceCoords_[i].Resize(2);
        sourceCoords_[i][0] = cCoord[0];
        sourceCoords_[i][1] = cCoord[1];
      }else{
        sourceCoords_[i].Resize(3);
        sourceCoords_[i][0] = cCoord[0];
        sourceCoords_[i][1] = cCoord[1];
        sourceCoords_[i][2] = cCoord[2];
      }
    }
  }else{
    //that's the case, if the source values are defined on src-nodes
    sourceCoords_.Resize(inGrid->GetNumNodes(ALL_REGIONS));
    CF::StdVector<UInt> srcNodes;
    CF::StdVector<CF::Elem*> srcElems;
    inGrid->GetElems(srcElems, ALL_REGIONS);

    inGrid->GetNodesOfElemList(srcNodes, srcElems, false);
    for(UInt nIter=0; nIter < sourceCoords_.GetSize(); ++nIter){
        CF::Vector<Double> nCoord;
        inGrid->GetNodeCoordinate3D(nCoord, srcNodes[nIter]);
        if(trgGrid_->GetDim() == 2){
          sourceCoords_[nIter].Resize(2);
          sourceCoords_[nIter][0] = nCoord[0];
          sourceCoords_[nIter][1] = nCoord[1];
        }else{
          sourceCoords_[nIter].Resize(3);
          sourceCoords_[nIter][0] = nCoord[0];
          sourceCoords_[nIter][1] = nCoord[1];
          sourceCoords_[nIter][2] = nCoord[2];
        }

      }

  }

  // Obtaining target coordinates
  std::cout << "\t\t 2/5 Obtaining target coordinates" << std::endl;
  CF::StdVector<const CF::Elem*> trgElements;
  CF::StdVector< LocPoint > locPoints;
  CF::StdVector<UInt> trgNodes;
  CF::StdVector<CF::Elem*> trgElems;
  trgGrid_->GetVolElems(trgElems, ALL_REGIONS);
  trgGrid_->GetNodesOfElemList(trgNodes, trgElems);
  targetCoords_.Resize(trgNodes.GetSize());

  for(UInt nIter=0; nIter < trgNodes.GetSize(); ++nIter){
      CF::Vector<Double> SCoord;
      trgGrid_->GetNodeCoordinate3D(SCoord, trgNodes[nIter]);
      if(trgGrid_->GetDim() == 2){
        targetCoords_[nIter].Resize(2);
        targetCoords_[nIter][0] = SCoord[0];
        targetCoords_[nIter][1] = SCoord[1];
      }else{
        targetCoords_[nIter].Resize(3);
        targetCoords_[nIter][0] = SCoord[0];
        targetCoords_[nIter][1] = SCoord[1];
        targetCoords_[nIter][2] = SCoord[2];
      }
    }

  trgGrid_->GetElemsAtGlobalCoords(targetCoords_,locPoints,trgElements, lists, 1e-6, 1e-3);


  std::cout << "\t\t 3/5 Generating interpolation info ..." << std::endl;
  UInt trgSize = trgGrid_->GetNumVolElems();
  interpolData_.reserve(trgSize);
  for(UInt i = 0; i < trgSize; ++i){
      InpolationStruct newStruct;
      shared_ptr<ElemShapeMap> eShape = trgGrid_->GetElemShapeMap(trgElements[i],true);
      newStruct.localCoords = locPoints[i].coord;
      newStruct.volume = eShape->CalcVolume();
      newStruct.tENum = trgElements[i]->elemNum;
      interpolData_.push_back(newStruct);
  }

  std::cout << "\t\t 4/5 Clear generated temporary data storage ..." << std::endl;
  trgElements.Clear(false);
  allSrcElems.clear();

  std::cout << "\t\t 5/5 Remap data to equation numbers ..." << std::endl;
  str1::shared_ptr<EqnMapSimple> upMap = resultManager_->GetResultAdpter(upRes)->mapping;
  CF::StdVector<UInt> sEqn;

  std::cout << "\t\t Interpolation prepared!" << std::endl;
}

ResultIdList NearestNeighbourInterpolator::SetUpstreamResults(){
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


#ifdef USE_CGAL
void NearestNeighbourInterpolator::KNNSearch_CGAL(const Vector<Double> globPoint,
  StdVector< Vector<Double> >& neighbors,
  StdVector< Double >& l2Distances,
  StdVector< Vector<Double> >& vectors)
{
CGAL::Point query(0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
if(globPoint.GetSize() == 2)
{
  CGAL::Point query2(globPoint[0], globPoint[1], 0.0, 0.0, 0.0, 0.0);
  query = query2;
}
else
{
  CGAL::Point query3(globPoint[0], globPoint[1], globPoint[2], 0.0, 0.0, 0.0);
  query = query3;
}

K_neighbor_search search(*searchTree_.get(), query, numNeighbors_);

K_neighbor_search::iterator it = search.begin();

if(it == search.end())
{
  EXCEPTION("Could not find a nearest neighbor for " << globPoint << "!");
}

UInt nn = std::distance(search.begin(), search.end());
neighbors.Resize(nn);
l2Distances.Resize(nn);
vectors.Resize(nn);

for(UInt i=0 ; it != search.end(); ++it, i++) {
  l2Distances[i] = std::sqrt(it->second);
  neighbors[i].Resize(globPoint.GetSize());
  vectors[i].Resize(globPoint.GetSize());

    if(l2Distances[i]<searchRadius_){
      if(globPoint.GetSize() == 2)
      {
        it->first.vx(vectors[i][0]);
        it->first.vy(vectors[i][1]);
        neighbors[i][0] = it->first.x();
        neighbors[i][1] = it->first.y();
      }
      else
      {
        it->first.vx(vectors[i][0]);
        it->first.vy(vectors[i][1]);
        it->first.vz(vectors[i][2]);

        neighbors[i][0] = it->first.x();
        neighbors[i][1] = it->first.y();
        neighbors[i][2] = it->first.z();
      }
    }else{
      if(globPoint.GetSize() == 2)
      {
        vectors[i][0] = 0.0;
        vectors[i][1] = 0.0;
        neighbors[i][0] = it->first.x();
        neighbors[i][1] = it->first.y();
      }
      else
      {
        vectors[i][0] = 0.0;
        vectors[i][1] = 0.0;
        vectors[i][2] = 0.0;

        neighbors[i][0] = it->first.x();
        neighbors[i][1] = it->first.y();
        neighbors[i][2] = it->first.z();
      }
    }


}
}
#endif

#ifdef USE_FLANN
void NearestNeighbourInterpolator::KNNSearch_FLANN(const Vector<Double> globPoint,
      StdVector< Vector<Double> >& neighbors,
      StdVector< Double >& l2Distances,
      StdVector< Vector<Double> >& vectors,
      StdVector< Vector<Double> > scatteredData)
  {
  Double q[3];

  if(globPoint.GetSize()==2)
  {
    q[0] = globPoint[0];
    q[1] = globPoint[1];
    q[2] = 0.0;
  }
  else
  {
    q[0] = globPoint[0];
    q[1] = globPoint[1];
    q[2] = globPoint[2];
  }

  flann::Matrix<Double> query(q, 1,3);

  flann::Matrix<int> indices(new int[query.rows*numNeighbors_], query.rows, numNeighbors_);
  flann::Matrix<Double> dists(new Double[query.rows*numNeighbors_], query.rows, numNeighbors_);
  // do a knn search, using 3 checks
  index_->knnSearch(query, indices, dists, numNeighbors_, flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));

  neighbors.Resize(numNeighbors_);
  l2Distances.Resize(numNeighbors_);
  vectors.Resize(numNeighbors_);
  for(UInt i=0; i<indices.rows; i++)
  {
    for(UInt j=0; j<numNeighbors_; j++)
    {
      l2Distances[j] = std::sqrt(dists[i][j]);


        if(l2Distances[j]<searchRadius_){
          UInt idx = indices[i][j];

          neighbors[j].Resize(inDim_+1);
          vectors[j].Resize(inDim_+1);

          for(UInt d=0; d<inDim_+1; d++)
          {
            vectors[j][d] = scatteredData[idx][d];
            neighbors[j][d] = sourceCoords_[idx][d];
          }
        }else{
          UInt idx = indices[i][j];

          neighbors[j].Resize(inDim_+1);
          vectors[j].Resize(inDim_+1);

          for(UInt d=0; d<inDim_+1; d++)
          {
            vectors[j][d] = scatteredData[idx][d] * 0.0;
            neighbors[j][d] = sourceCoords_[idx][d];
          }
        }


    }
  }

  delete[] indices.ptr();
  delete[] dists.ptr();
  }
#endif


}

