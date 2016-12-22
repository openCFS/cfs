// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     DivergenceDifferentiator.cc
 *       \brief    <Description>
 *
 *       \date     Oct 6, 2016
 *       \author   kroppert
 */
//================================================================================================


#include "DivergenceDifferentiator.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include "cfsdat/Utils/KNNSearch.hh"
#include <algorithm>
#include <vector>
#include <math.h>


namespace CFSDat{

DivergenceDifferentiator::DivergenceDifferentiator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
                     :MeshBasedDerivative(numWorkers,config,resMan){


  this->filtSteamType_ = FIFO_FILTER;
  inDim_ = 0;

}

DivergenceDifferentiator::~DivergenceDifferentiator(){

}

bool DivergenceDifferentiator::Run(){
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

  // vector which will be filled with the source values (scattered data values), in order to
  // perform a nearest neighbour search
  CF::StdVector< CF::Vector<Double> >  scatteredData;
  scatteredData.Resize(sourceCoords_.GetSize());

  // nr. of nearest neighbours for the nodal RBF basis functions
  UInt numNN = (trgGrid_->GetDim() == 2) ? 4 : 8 ;


  // coordinate list of nearest neighbour points
  CF::StdVector< Vector<Double> > neighbors;

  // distances according to every nearest neighbour point
  CF::StdVector< Double > l2dists;

  // vector containing the values of each nearest neighbour point
  CF::StdVector< Vector<Double> > vectors;

  // actually this is the scalar value of the divergence, but Matrix class has to be used, for matrix-vector product
  CF::Matrix<Double> vec;
  vec.Resize(1);
  vec[0][0] = 0.0;

  str1::shared_ptr<EqnMapSimple> downMap = resultManager_->GetResultAdpter(filterResIds[0])->mapping;

  // Checking if input is scalar- or vector-type. This needs to be done
  // because of the input of CGAL
  // Maybe there is a more efficient way to deal with this issue?!
  if(inVec.GetSize() == sourceCoords_.GetSize()){
      //this is the case of scalar scattered values
      EXCEPTION("Divergence of a scalar field!");
    }else{
          if(inVec.GetSize() == 2 * sourceCoords_.GetSize()){
            //case of two dimensional vector
            inDim_=1;
#pragma omp parallel for shared(inVec)
            for(UInt i = 0; i < sourceCoords_.GetSize() ; ++i){
            scatteredData[i].Resize(2);
            scatteredData[i][0] = inVec[2 * i]; // x-component
            scatteredData[i][1] = inVec[2 * i + 1]; // y-component
            }
          }else{
              if(inVec.GetSize() == 3 * sourceCoords_.GetSize()){
                // case of three dimensional vector
                inDim_=2;
#pragma omp parallel for shared(inVec)
                for(UInt i = 0; i < sourceCoords_.GetSize(); ++i){
                scatteredData[i].Resize(3);
                scatteredData[i][0] = inVec[3 * i]; // x-component
                scatteredData[i][1] = inVec[3 * i + 1]; // y-component
                scatteredData[i][2] = inVec[3 * i + 2]; // z-component
                }
              }else{
                EXCEPTION("Incorrect Input Data!")
              }
          }
    }

  //now we bring the scattered coordinates and data values into
  //the correct form for the nearest neighbour-search witch CGAL

  // Object for nearest neighbor-searches and bringing the data into the correct form for CGAL search
  KNNSearch Tree;
  Tree.ReadScatteredData_Div(sourceCoords_, inDim_, trgGrid_, scatteredData);

  CF::StdVector<CF::Elem*> elems;
  trgGrid_->GetVolElems(elems, ALL_REGIONS);


  // loop over all elements
for(UInt i = 0; i < derivData_.size();++i){

    DifferentiationStruct& aStru = derivData_[i];
    const Elem* curE = trgGrid_->GetElem(aStru.tENum);

    //get the global coordinates of element centroid
    CF::Vector<Double> globalCoord;
    trgGrid_->GetElemCentroid(globalCoord, curE->elemNum , false);

          // at that point we can start obtaining the nearest neighbours for every aNode
          Tree.KNN_CGAL_Differentiation(globalCoord, neighbors, l2dists, vectors, numNN);

          // concerning radius factor alpha: we want to scale min(l2dists)/alpha~=1.5e-03
          // this means alpha = min(l2dists) * 1.5e+03
          // the first entry of l2dists is zero, if target and source node coincide, then we
          // simply take the next entry, which must be != 0
          Double alpha;
          if( l2dists[0] == 0.0 ){
            alpha = l2dists[1] * 1.5e+03;
          }else{
            alpha = l2dists[0] * 1.5e+03;
          }


          // now we have to calculate the local RBF interpolation matrix with the nn nodes
          // and also do the inversion
          CalcLocRBFDerivativeCoefs(vec, globalCoord, neighbors, l2dists, vectors, numNN, alpha);

          CF::StdVector<UInt> eqns;
            //get the equation map for the nodes in eConn, in order to insert the
            //interpolation in the correct position in the result vector
            downMap->GetEquation(eqns,curE->elemNum,ExtendedResultInfo::ELEMENT);

            //if scalar input values of scattered data->inDim_=0
            //if it is a two-dimensional vector->inDim_=1 and inDim_=2 for 3d-vector
            for(UInt aDof = 0; aDof < eqns.GetSize(); aDof++){
              returnVec[eqns[aDof]] = vec[0][0];
            }

  }// for element
  resultManager_->ActivateResult(filterResIds[0]);

  //now deactivate own upstream results
  for(UInt aRes=0;aRes<upResIds.GetSize();aRes++){
    resultManager_->DeactivateResult(upResIds[aRes]);
  }

  return true;
}



void DivergenceDifferentiator::CalcLocRBFDerivativeCoefs(CF::Matrix<Double>& vec,
                                      CF::Vector<Double>& globPoint,
                                      CF::StdVector< Vector<Double> >& neighbors,
                                      CF::StdVector< Double >& l2Distances,
                                      CF::StdVector< Vector<Double> >& vectors,
                                      UInt numNN,
                                      Double alpha){
  CF::Matrix<Double> derivCoefVec;
  CF::Matrix<Double> ALoc;
  ALoc.Resize(numNN,numNN);
  CF::Matrix<Double> vals;
  CF::Matrix<Double> derivVec; //Vector of RBF derivatives evaluated at srcPoints

  Double rNN; //distance between two src points
  for (UInt i = 0; i < numNN; ++i){
    for (UInt j = 0; j < numNN; ++j){
      if (trgGrid_->GetDim() == 3){
        rNN = sqrt(pow(neighbors[i][0]-neighbors[j][0],2.0) + pow(neighbors[i][1]-neighbors[j][1],2.0) + pow(neighbors[i][2]-neighbors[j][2],2.0));
      }else{
        rNN = sqrt(pow(neighbors[i][0]-neighbors[j][0],2.0) + pow(neighbors[i][1]-neighbors[j][1],2.0));
      }
      ALoc[i][j] = pow(1.0 - rNN/alpha, 2.0);
    }
    switch( inDim_ ){
    case 0:
      //should already be caught
      EXCEPTION("Divergence of a scalar field!");
      break;
    case 1:
      vals.Resize(numNN,2);
      vals[i][0] = vectors[i][0];
      vals[i][1] = vectors[i][1];
      if (trgGrid_->GetDim() == 2){
      derivVec.Resize(numNN,2);
      if (l2Distances[i] == 0) {
        derivVec[i][0] = 0.0;
        derivVec[i][1] = 0.0;
      }else{
        derivVec[i][0] = 2.0 * (1.0 - l2Distances[i]/alpha) * (neighbors[i][0] - globPoint[0])/(l2Distances[i]*alpha);
        derivVec[i][1] = 2.0 * (1.0 - l2Distances[i]/alpha) * (neighbors[i][1] - globPoint[1])/(l2Distances[i]*alpha);
      }
      }else{
        EXCEPTION("3D values and 2D mesh!");
      }
      break;
    case 2:
      vals.Resize(numNN,3);
      vals[i][0] = vectors[i][0];
      vals[i][1] = vectors[i][1];
      vals[i][2] = vectors[i][2];
      if (trgGrid_->GetDim() == 3){
      derivVec.Resize(numNN,3);
      if (l2Distances[i] == 0) {
        derivVec[i][0] = 0.0;
        derivVec[i][1] = 0.0;
        derivVec[i][2] = 0.0;
      }else{
        derivVec[i][0] = 2.0 * (1.0 - l2Distances[i]/alpha) * (neighbors[i][0] - globPoint[0])/(l2Distances[i]*alpha);
        derivVec[i][1] = 2.0 * (1.0 - l2Distances[i]/alpha) * (neighbors[i][1] - globPoint[1])/(l2Distances[i]*alpha);
        derivVec[i][2] = 2.0 * (1.0 - l2Distances[i]/alpha) * (neighbors[i][2] - globPoint[2])/(l2Distances[i]*alpha);
      }
      }else{
        //case of 2D mesh and 3D values
        vals.Resize(numNN,3);
        vals[i][0] = vectors[i][0];
        vals[i][1] = vectors[i][1];
        derivVec.Resize(numNN,3);
        if (l2Distances[i] == 0) {
          derivVec[i][0] = 0.0;
          derivVec[i][1] = 0.0;
          derivVec[i][2] = 0.0;
        }else{
          derivVec[i][0] = 2.0 * (1.0 - l2Distances[i]/alpha) * (neighbors[i][0] - globPoint[0])/(l2Distances[i]*alpha);
          derivVec[i][1] = 2.0 * (1.0 - l2Distances[i]/alpha) * (neighbors[i][1] - globPoint[1])/(l2Distances[i]*alpha);
          derivVec[i][2] = 0.0;
        }
      }
      break;
    }
 }

  // now we have to invert Aloc and multiply it with the according value-coloumn
  ALoc.Invert_Lapack();

  CF::Matrix<Double> temp;

  // coefficient matrix (coloumn nr. corresponding to the spatial dimension)
  derivCoefVec = ALoc * derivVec;
  vals.Transpose(temp);

  CF::Matrix<Double> tempvec;
  tempvec.Resize(1,inDim_+1);
  tempvec = temp * derivCoefVec;
  //in order to obtain the divergence of the vectorfield, we have to add tempvec[:,0]+tempvec[:,1]+tempvec[:,2]
  vec[0][0] = tempvec.Trace();

}


void DivergenceDifferentiator::PrepareDifferentiation(){
  //1. Get get the source coordinates and the values, defined on those coordinates (Source...Src)
  //2. Get the target coordinates (trg)
  //3. Store for each trg element local Coordinates, elem number, volume, ...
  //4. Small clean-up


  std::cout << "\t ---> DivergenceDifferentiator preparing for interpolation" << std::endl;

  //in this filter we only have one upstream result
  uuids::uuid upRes = upResIds[0];

  Grid* inGrid   = resultManager_->GetExtInfo(upRes)->ptGrid;


  StdVector<CF::UInt> allSrcElems;
  StdVector<CF::UInt> allSrcNodes;
  //loop over source regions and add element numbers and nodes to vector
  std::set<std::string>::iterator sRegIter = srcRegions_.begin();
  for(;sRegIter != srcRegions_.end();++sRegIter){
    StdVector<CF::UInt> curElems;
    StdVector<CF::UInt> nodeList;
    inGrid->GetElemNumsByName(curElems,*sRegIter);
    inGrid->GetNodesByName(nodeList, *sRegIter);
    //insert the element numbers into the allSrcElems list
    //unfortunately we have to loop over all entries, because in StdVector,
    //only single entries can be "pushed back"
    for (UInt eIter = 0; eIter < curElems.GetSize(); ++eIter){
      allSrcElems.Push_back(curElems[eIter]);
    }
    //do the same for the nodes
    for (UInt nIter = 0; nIter < nodeList.GetSize(); ++nIter){
      allSrcNodes.Push_back(nodeList[nIter]);
    }
  }


  //loop over target regions and add element numbers to vector and also nodes
  StdVector<const CF::Elem*> allTrgElems;
  StdVector<CF::UInt> allTrgElemNums;
  StdVector<CF::UInt> allTrgNodes;
  StdVector<shared_ptr<EntityList> > lists;
  std::set<std::string>::iterator tRegIter = trgRegions_.begin();
  for(;tRegIter != trgRegions_.end();++tRegIter){
    StdVector<UInt> curElemNums;
    StdVector<UInt> nodeList;

    trgGrid_->GetElemNumsByName(curElemNums,*tRegIter);
    trgGrid_->GetNodesByName(nodeList, *tRegIter);
    //insert the element numbers and elements into the allTrgElems list
    //unfortunately we have to loop over all entries, because in StdVector,
    //only single entries can be "pushed back"
    for (UInt eIter = 0; eIter < curElemNums.GetSize(); ++eIter){
      allTrgElemNums.Push_back(curElemNums[eIter]);

      const CF::Elem* curElem;
      curElem = trgGrid_->GetElem(curElemNums[eIter]);
      allTrgElems.Push_back(curElem);
    }
    //do the same for the nodes
    for (UInt nIter = 0; nIter < nodeList.GetSize(); ++nIter){
      allTrgNodes.Push_back(nodeList[nIter]);
    }
    RegionIdType aReg = trgGrid_->GetRegion().Parse(*tRegIter);
    shared_ptr<ElemList> newList(new ElemList(trgGrid_));
    newList->SetRegion(aReg);
    lists.Push_back(newList);
  }
  std::cout << "\t\t\t Differentiator is dealing with " << allSrcElems.GetSize() <<
               " source elements and "<< allTrgElems.GetSize() << " target elements" << std::endl;


  std::cout << "\t\t 1/5 Obtaining source coordinates " << std::endl;
  //here we also find out, if the input is defined on elements (-centroids) or nodes
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  if(inInfo->definedOn == ExtendedResultInfo::ELEMENT){
    //case of centroid-values
    EXCEPTION("============================================================ \n"
              "Input of DivergenceDifferentiator-filter must be defined on nodes!! \n"
              "Just use an interpolator to transform to node-results \n"
              "============================================================")
  }else{
    //that's the case, if the source values are defined on src-nodes
    sourceCoords_.Resize(allSrcNodes.GetSize());
    for(UInt nIter=0; nIter < allSrcNodes.GetSize(); ++nIter){
        CF::Vector<Double> nCoord;
        //inGrid->GetNodeCoordinate3D(nCoord, srcNodes[nIter]);
        inGrid->GetNodeCoordinate3D(nCoord, allSrcNodes[nIter]);
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
  // target (output) is solely defined on nodes
  targetCoords_.Resize(allTrgNodes.GetSize());
  for(UInt nIter=0; nIter < allTrgNodes.GetSize(); ++nIter){
      CF::Vector<Double> SCoord;
      trgGrid_->GetNodeCoordinate3D(SCoord, allTrgNodes[nIter]);
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

  CF::StdVector< LocPoint > locPoints;
  //tempElems are just dummy vectors, get deleted right after we get the local coordinates
  StdVector<const CF::Elem*> tempElems;
  //mapping of global point targetCoords_ to local locPoints
  trgGrid_->GetElemsAtGlobalCoords(targetCoords_,locPoints, tempElems, lists, 1e-6, 1e-3);
  //tempElems.Clear();

  std::cout << "\t\t 3/5 Generating differentiation info ..." << std::endl;
  derivData_.reserve(tempElems.GetSize());
  for(UInt aMatch = 0;aMatch < tempElems.GetSize();++aMatch){
    if(tempElems[aMatch]!= NULL){
      DifferentiationStruct newStruct;
      newStruct.localCoords = locPoints[aMatch].coord;
      newStruct.tENum = tempElems[aMatch]->elemNum;
      derivData_.push_back(newStruct);
    }
  }


  std::cout << "\t\t 4/5 Clear generated temporary data storage ..." << std::endl;
  allTrgElems.Clear(false);
  allSrcElems.Clear(false);
  allSrcNodes.Clear(false);
  allTrgElemNums.Clear(false);
  allTrgNodes.Clear(false);

  std::cout << "\t\t 5/5 Remap data to equation numbers ..." << std::endl;
  str1::shared_ptr<EqnMapSimple> upMap = resultManager_->GetResultAdpter(upRes)->mapping;

  std::cout << "\t\t Differentiation prepared!" << std::endl;
}

ResultIdList DivergenceDifferentiator::SetUpstreamResults(){
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

void DivergenceDifferentiator::AdaptFilterResults(){
  //some checks
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  if(!inInfo->isValid){
    EXCEPTION("Could not validate required input result \"" << inInfo->resultName << "\" from upstream filters.");
  }
  //we require mesh result input
  if(!inInfo->isMeshResult){
    EXCEPTION("DivergenceDifferentiator requires input to be defined on mesh");
  }

  //got the upstream result validated?
  if(!inInfo->isValid){
    EXCEPTION("Problem in filter pipeline detected. Differentiator input result \"" <<  inInfo->resultName << "\" could not be provided.")
  }


  resultManager_->CopyResultData(upResIds[0],filterResIds[0]);
  //but now, we need to overwrite some things
  resultManager_->SetRegionNames(filterResIds[0],this->trgRegions_);

  resultManager_->SetEntryType(filterResIds[0],ExtendedResultInfo::SCALAR);

  CF::StdVector<std::string> dofnames;
  dofnames.Push_back("x");
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

