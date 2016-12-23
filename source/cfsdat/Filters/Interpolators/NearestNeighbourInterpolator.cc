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
#include "cfsdat/Utils/KNNSearch.hh"
#include <algorithm>
#include <vector>

namespace CFSDat{

NearestNeighbourInterpolator::NearestNeighbourInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
                     :MeshBasedInterpolator(numWorkers,config,resMan){

  this->filtStreamType_ = FIFO_FILTER;
  inDim_ = 0;

  p_ = config->Get("scheme")->Get("interpolationExponent")->As<UInt>();
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


  // vector, containing the source data values
  Vector<Double>& inVec = resultManager_->GetResultVector<Double>(upResIds[0],eqnNums);

  // vector which will be filled with the source values (scattered data values), in order to
  // perform a nearest neighbour search
  CF::StdVector< CF::Vector<Double> >  scatteredData;
  scatteredData.Resize(sourceCoords_.GetSize());


  // vector containing the interpolated result for every trg node
  // can be 1D for scalar values or 2D/3D for vector values
  CF::Vector<Double> vec;

  str1::shared_ptr<EqnMapSimple> downMap = resultManager_->GetResultAdapter(filterResIds[0])->mapping;

  // Checking if input is scalar- or vector-type. This needs to be done
  // because of the input of CGAL and FLANN
  // Maybe there is a more efficient way to deal with this issue?!

  if(inVec.GetSize() == sourceCoords_.GetSize()){
    //this is the case of scalar scattered values
    inDim_=0;
    vec.Resize(1);
#pragma omp parallel for
    for(UInt i = 0; i < sourceCoords_.GetSize(); ++i){
      scatteredData[i].Resize(1);
      scatteredData[i][0] = inVec[i];
     }
    }else{
          if(inVec.GetSize() == 2 * sourceCoords_.GetSize()){
            //case of two dimensional vector
            inDim_=1;
            vec.Resize(2);
#pragma omp parallel for
            for(UInt i = 0; i < sourceCoords_.GetSize() ; ++i){
            scatteredData[i].Resize(2);
            scatteredData[i][0] = inVec[2 * i]; // x-component
            scatteredData[i][1] = inVec[2 * i + 1]; // y-component
            }
          }else{
              if(inVec.GetSize() == 3 * sourceCoords_.GetSize()){
                // case of three dimensional vector
                inDim_=2;
                vec.Resize(3);
#pragma omp parallel for
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

  // Object for nearest neighbor-searches and bringing the data into the correct form for CGAL search
  KNNSearch Tree;
  Tree.ReadScatteredData_Interpolation(sourceCoords_, inDim_, trgGrid_, scatteredData);

  CF::StdVector<bool> nodeCheck;
  nodeCheck.Resize(targetCoords_.GetSize());
  nodeCheck.Init(false);

  // loop over all elements and over every node of each element
  for(UInt i=0;i < interpolData_.size();++i){
    UInt nodeIter = 0;
    InpolationStruct& aStru = interpolData_[i];
    const Elem* curE = trgGrid_->GetElem(aStru.tENum);
    CF::shared_ptr<ElemShapeMap> eShape = trgGrid_->GetElemShapeMap(curE,true);
    //eConn gives us the node numbers of current element curE
    const CF::StdVector<UInt>& eConn = curE->connect;


    //get the global coordinates of the nodes of element curE
    CF::Matrix<Double> globalCoords;
    trgGrid_->GetElemNodesCoord(globalCoords, eConn, true);
    FeH1 * myElem = dynamic_cast<FeH1*>(eShape->GetBaseFE());
    CF::Vector<Double> shFnc;
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
      //eConn[aNode]-1 because nodeNumbers of eConn start with 1 and not with 0 !!
      if(nodeCheck[eConn[aNode]-1] == false){
        //add aNode to the "already-computed-list"
        nodeCheck[eConn[aNode]-1] = true;

        // coordinate list of nearest neighbour points
        CF::StdVector< Vector<Double> > neighbors;
        // distances according to every nearest neighbour point
        CF::StdVector< Double > l2dists;
        // vector containing the values of each nearest neighbour point
        CF::StdVector< Vector<Double> > vectors;

        for (UInt dof = 0; dof < inDim_ + 1; dof++){
          vec[dof] = 0.0;
        }

        // at that point we can start obtaining the nearest neighbours for every aNode
        Tree.KNN_CGAL_Interpolation(globPoint, neighbors, l2dists, vectors, numNeighbors_);

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
            CF::StdVector<UInt> eqns;
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
std::cout<<returnVec.GetSize()<<std::endl;
  return true;
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

  std::cout << "\t\t\t Interpolator is dealing with " << allSrcElems.GetSize() <<
               " source elements and "<< allTrgElems.GetSize() << " target elements" << std::endl;





  std::cout << "\t\t 1/5 Obtaining source coordinates " << std::endl;
  //here we also find out, if the input is defined on elements (-centroids) or nodes
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  if(inInfo->definedOn == ExtendedResultInfo::ELEMENT){
    //case of centroid-values
    sourceCoords_.Resize(allSrcElems.GetSize());
    for(UInt i=0;i<allSrcElems.GetSize();++i){
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
  //  tempElems.Clear();

  std::cout << "\t\t 3/5 Generating interpolation info ..." << std::endl;

    interpolData_.reserve(allTrgElems.GetSize());
    for(UInt aMatch = 0;aMatch < tempElems.GetSize();++aMatch){
      if(tempElems[aMatch]!= NULL){
        InpolationStruct newStruct;
        newStruct.localCoords = locPoints[aMatch].coord;
        newStruct.tENum = tempElems[aMatch]->elemNum;
        interpolData_.push_back(newStruct);
      }
    }


  std::cout << "\t\t 4/5 Clear generated temporary data storage ..." << std::endl;
  allTrgElems.Clear(false);
  allSrcElems.Clear(false);
  allSrcNodes.Clear(false);
  allTrgElemNums.Clear(false);
  allTrgNodes.Clear(false);

  std::cout << "\t\t 5/5 Remap data to equation numbers ..." << std::endl;
  str1::shared_ptr<EqnMapSimple> upMap = resultManager_->GetResultAdapter(upRes)->mapping;


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

}

