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

#include <algorithm>
#include <vector>
#include <math.h>


namespace CFSDat{

RBFInterpolator::RBFInterpolator(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
                     :MeshBasedInterpolator(numWorkers,config,resMan){

  this->filtSteamType_ = FIFO_FILTER;
  inDim_ = 0;

  //check if the nearest neighbour search scheme has been defined in the xml-file
  if(config->Has("scheme") == false){
    EXCEPTION("============================================================ \n"
              "No <sheme> for nearest neighbour search was set in xml-file! \n"
              "============================================================");
  }

  std::string alg = config->Get("scheme")->Get("searchAlgorithm")->As<std::string>();
  if(alg == "CGAL"){
    knnLib_ = 0;
  }else if(alg == "FLANN"){
    knnLib_ = 1;
  }else{
    EXCEPTION("no search algorithm selected!");
  }
  p_ = config->Get("scheme")->Get("interpolationExponent")->As<UInt>();
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

  // vector which will be filled with the source values (scattered data values), in order to
  // perform a nearest neighbour search
  CF::StdVector< CF::Vector<Double> >  scatteredData;
  scatteredData.Resize(sourceCoords_.GetSize());

  // find out the dimension of the mesh (2D, 3D) and assign the appropriate number of nearest neighbours
  // nr. of nearest neighbours for the nodal RBF basis functions
  UInt numNN = (trgGrid_->GetDim() == 2) ? 13 : 18 ;
  //nr. of nearest neighbours needed for the weight-calculation (must be smaller than numNN !!)
  UInt numNW = (trgGrid_->GetDim() == 3) ? 10 : 7 ;

  // vector containing the interpolated result for every trg node
  // can be 1D for scalar values or 2D/3D for vector values
  CF::Vector<Double> vec;
  CF::Vector<Double> R_k;

  str1::shared_ptr<EqnMapSimple> downMap = resultManager_->GetResultAdpter(filterResIds[0])->mapping;

  // Checking if input is scalar- or vector-type. This needs to be done
  // because of the input of CGAL and FLANN
  // Maybe there is a more efficient way to deal with this issue?!

  if(inVec.GetSize() == sourceCoords_.GetSize()){
    //this is the case of scalar scattered values
    inDim_=0;
    vec.Resize(1);
    R_k.Resize(1);
    for(UInt i = 0; i < sourceCoords_.GetSize(); ++i){
      scatteredData[i].Resize(1);
      scatteredData[i][0] = inVec[i];
     }
    }else{
          if(inVec.GetSize() == 2 * sourceCoords_.GetSize()){
            //case of two dimensional vector
            inDim_=1;
            vec.Resize(2);
            R_k.Resize(2);
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
                R_k.Resize(3);
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




  //now we can bring the scattered coordinates and data values into
  //the correct form for the nearest neighbour-search witch CGAL or FLANN
ReadScatteredData(sourceCoords_, scatteredData);



  CF::StdVector<UInt> nodeCheck;
  nodeCheck.Resize(0);
  nodeCheck.Init();
  bool nodeMatch;


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

    for(UInt aNode = 0;aNode < eConn.GetSize(); ++aNode){
      //extract the global point of node aNode
      CF::Vector<Double> globPoint;
      globalCoords.GetCol(globPoint, aNode);

      //that is the mentioned search, to find out if the node has already been used
      nodeMatch = nodeCheck.Contains(eConn[aNode]);
      if(nodeMatch == false){
        //add aNode to the "already-computed-list"
        nodeCheck.Push_back(eConn[aNode]);

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
          // using CGAL or FLANN
          switch(knnLib_)
          {
          case 0: //CGAL
        #ifdef USE_CGAL
            KNNSearch_CGAL(globPoint, neighbors, l2dists, vectors, numNN);
        #else
            EXCEPTION("CoefFunctionScatteredData needs to be compiled with USE_CGAL=ON!");
        #endif
            break;

          case 1: //FLANN
        #ifdef USE_FLANN
            KNNSearch_FLANN(globPoint, neighbors, l2dists, vectors, scatteredData, numNN);

        #else
            EXCEPTION("CoefFunctionScatteredData needs to be compiled with USE_FLANN=ON!");
        #endif
            break;
          default:
            break;
          }

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


          // local RBF interpolation coefficients vector (matrix only because
          // of the following ALoc*coefVec multiplication)
          CF::Matrix<Double> coefVec;
          CalcLocRBFCoefs(coefVec, globPoint, neighbors, l2dists, vectors, numNN, alpha);

          Double r_k;
          Double W_k_bar;
          Double temp;
          Double r_Wk = l2dists[numNW]; //influence radius of the weight function
          Double W_denom = 0.0;
      //now we can perform the computation of the weight function
      for (UInt srcIter = 0; srcIter < numNN; ++srcIter){
            r_Wk = l2dists[numNW]; //influence radius of the weight function
            W_denom = 0.0;
            for (UInt k = 0; k < numNN; ++k){
              r_k = l2dists[k];
              temp = (r_Wk - r_k);
              if (temp < 0.0){temp = 0.0;}
              if (r_k == 0.0){
                W_denom += 0.0;
              }else{
                W_denom += pow(temp / (r_Wk * r_k), Double(p_));
               }
              }
          r_k = l2dists[srcIter];
          temp = (r_Wk - r_k);
          if (temp < 0.0){temp = 0.0;}
          if (r_k == 0.0){
            W_k_bar = 1.0 / W_denom;
          }else{
            W_k_bar = pow(temp / (r_Wk * r_k), Double(p_)) / W_denom;
          }

          //RBF evaluation
          //just a clear-up of of R_k
          for (UInt dof = 0; dof < inDim_ + 1; dof++){
          R_k[dof] = 0.0;
          }

          for(UInt dof=0; dof < inDim_ + 1; dof++){
            for (UInt m = 0; m < numNN; ++m){
              R_k[dof] += coefVec[m][dof] * pow(1.0 - l2dists[m]/alpha, 2.0);
            }
          }

          for(UInt dof=0; dof < inDim_ + 1; dof++){
            vec[dof] += R_k[dof] * W_k_bar;
          }

      }
      CF::StdVector<UInt> eqns;
            //get the equation map for the nodes in eConn, in order to insert the
            //interpolation in the correct position in the result vector
            downMap->GetEquation(eqns,eConn[aNode],ExtendedResultInfo::NODE);

            //if scalar input values of scattered data->inDim_=0
            //if it is a two-dimensional vector->inDim_=1 and inDim_=2 for 3d-vector
            for(UInt aDof = 0; aDof < eqns.GetSize(); aDof++){
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



void RBFInterpolator::CalcLocRBFCoefs(CF::Matrix<Double>& coefVec,
                                      CF::Vector<Double>& globPoint,
                                      CF::StdVector< Vector<Double> >& neighbors,
                                      CF::StdVector< Double >& l2Distances,
                                      CF::StdVector< Vector<Double> >& vectors,
                                      UInt numNN,
                                      Double alpha){
  coefVec.Resize(numNN,1);
  CF::Matrix<Double> ALoc;
  ALoc.Resize(numNN,numNN);
  CF::Matrix<Double> vals;



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
      vals.Resize(numNN,1);
      vals[i][0] = vectors[i][0];
      break;
    case 1:
      vals.Resize(numNN,2);
      vals[i][0] = vectors[i][0];
      vals[i][1] = vectors[i][1];
      break;
    case 2:
      vals.Resize(numNN,3);
      vals[i][0] = vectors[i][0];
      vals[i][1] = vectors[i][1];
      vals[i][2] = vectors[i][2];
      break;
    }
 }

  // now we have to invert Aloc and multiply it with the according value-coloumn
  ALoc.Invert_Lapack();

  // coefficient matrix (coloumn nr. corresponding to the spatial dimension)
  coefVec = ALoc * vals;
}



void RBFInterpolator::ReadScatteredData(CF::StdVector< CF::Vector<Double> > sourceCoords,
                                        CF::StdVector< CF::Vector<Double> > scatteredData)
{

  UInt n = sourceCoords_.GetSize();
  switch(knnLib_)
  {
  case 0: //CGAL
#ifdef USE_CGAL
    {
UInt i=0;

    std::vector<CGAL::Point> points;
    points.resize(n);

    for( i=0; i<n; i++)
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
    }
    searchTree_.reset(new Tree(points.begin(), points.end()));

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

}

void RBFInterpolator::PrepareInterpolation(){
  //1. Get get the source coordinates and the values, defined on those coordinates (Source...Src)
  //2. Get the target coordinates (trg)
  //3. Store for each trg element local Coordinates, elem number, volume, ...
  //4. Small clean-up


  std::cout << "\t ---> RBFInterpolator preparing for interpolation" << std::endl;

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
  tempElems.Clear();

  std::cout << "\t\t 3/5 Generating interpolation info ..." << std::endl;

  interpolData_.reserve(allTrgElems.GetSize());
  for(UInt aMatch = 0;aMatch < allTrgElems.GetSize();++aMatch){
    if(allTrgElems[aMatch]!= NULL){
      InpolationStruct newStruct;
      newStruct.localCoords = locPoints[aMatch].coord;
      newStruct.tENum = allTrgElems[aMatch]->elemNum;
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
  str1::shared_ptr<EqnMapSimple> upMap = resultManager_->GetResultAdpter(upRes)->mapping;


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


#ifdef USE_CGAL
void RBFInterpolator::KNNSearch_CGAL(const Vector<Double> globPoint,
  StdVector< Vector<Double> >& neighbors,
  StdVector< Double >& l2Distances,
  StdVector< Vector<Double> >& vectors,
  UInt numNN)
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

K_neighbor_search search(*searchTree_.get(), query, numNN);

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

}
}
#endif

#ifdef USE_FLANN
void RBFInterpolator::KNNSearch_FLANN(const Vector<Double> globPoint,
      StdVector< Vector<Double> >& neighbors,
      StdVector< Double >& l2Distances,
      StdVector< Vector<Double> >& vectors,
      StdVector< Vector<Double> > scatteredData,
      UInt numNN)
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

  flann::Matrix<int> indices(new int[query.rows*numNN], query.rows, numNN);
  flann::Matrix<Double> dists(new Double[query.rows*numNN], query.rows, numNN);
  // do a knn search, using 3 checks
  index_->knnSearch(query, indices, dists, numNN, flann::SearchParams(flann::FLANN_CHECKS_UNLIMITED));

  neighbors.Resize(numNN);
  l2Distances.Resize(numNN);
  vectors.Resize(numNN);
  for(UInt i=0; i<indices.rows; i++)
  {
    for(UInt j=0; j<numNN; j++)
    {
      l2Distances[j] = std::sqrt(dists[i][j]);

          UInt idx = indices[i][j];
          neighbors[j].Resize(globPoint.GetSize());
          vectors[j].Resize(inDim_+1);

          for(UInt d=0; d<globPoint.GetSize(); d++)
          {
            neighbors[j][d] = sourceCoords_[idx][d];
          }
          for(UInt d=0; d<inDim_+1; d++)
          {
            vectors[j][d] = scatteredData[idx][d];
          }


    }
  }

  delete[] indices.ptr();
  delete[] dists.ptr();
  }
#endif


}
