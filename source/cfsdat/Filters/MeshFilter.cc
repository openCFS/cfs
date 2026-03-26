// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     MeshFilter.cc
 *       \brief    <Description>
 *
 *       \date     Jan 11, 2017
 *       \author   kroppert
 */
//================================================================================================

#include <Filters/MeshFilter.hh>
#include "DataInOut/DefineInOutFiles.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include "FeBasis/H1/H1Elems.hh"


namespace CFSDat{

MeshFilter::MeshFilter(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan)
                      :BaseMeshFilterType(numWorkers,config,resMan){

  // Aeroacoustic-Filters choose their input and output results independently because
  // there we might have two inputs (i.e. velocity and vorticity)
  if (params_-> Get("type")->As<std::string>() != "AeroacousticSource_LambVector" &&
      params_-> Get("type")->As<std::string>() != "AeroacousticSource_LighthillSourceVector" &&
      params_-> Get("type")->As<std::string>() != "AeroacousticSource_LighthillSourceTerm" &&
      params_-> Get("type")->As<std::string>() != "AeroacousticSource_LighthillSourceTensor"){
        std::string inRes = params_->Get("singleResult")->Get("inputQuantity")->Get("resultName")->As<std::string>();
        std::string outRes = params_->Get("singleResult")->Get("outputQuantity")->Get("resultName")->As<std::string>();
        filtResNames.insert(outRes);
        upResNames.insert(inRes);
  }

  ParamNodeList srcList =  params_->Get("regions")->Get("sourceRegions")->GetList("region");
  for(UInt iP = 0; iP < srcList.GetSize(); ++iP){
    srcRegions_.insert(srcList[iP]->Get("name")->As<std::string>());
  }

  ParamNodeList trgList =  params_->Get("regions")->Get("targetRegions")->GetList("region");
  for(UInt iP = 0; iP < trgList.GetSize(); ++iP){
    trgRegions_.insert(trgList[iP]->Get("name")->As<std::string>());
  }

  //Now we read the input grid. Another possibility would be to give the user the opportunity
  //to perform interpolation/differentiation onto a grid which already has results
  //a grid only input filter is not directly possible due to the concept of a result driven pipeline.

  //create grid
  CreateDummyCfsParamNode();
  PtrParamNode infoNode;
  std::string filename = params_->Get("targetMesh")->GetChild()->Get("fileName")->As<std::string>();
  trgMeshInp_ = CoupledField::DefineInOutFiles::CreateSingleInputFileObject(filename,filterId_,params_->Get("targetMesh")->GetChild(),infoNode);
  trgMeshInp_->InitModule();
  UInt maxDim = trgMeshInp_->GetDim();

  trgGrid_ = new CF::GridCFS(maxDim,dummyXMLNode_,infoNode,filterId_);


  trgMeshInp_->ReadMesh(trgGrid_);
  //it would be nice not to finish the grid here
  //in order to let other filters add some entities
  //unfortunately this is not possible as we can not access anything
  //without it another question, how can two inputs share a common grid?
  trgGrid_->FinishInit();

}

//TODO for now copy paste code from inputfilter... not very nice
void MeshFilter::CreateDummyCfsParamNode(){

  PtrParamNode meshInputNode = params_->Get("targetMesh")->GetChild();
  dummyXMLNode_.reset(new ParamNode( ParamNode::PASS, ParamNode::ELEMENT));
  CoupledField::PtrParamNode iNode = dummyXMLNode_->Get("fileFormats",ParamNode::APPEND)->Get("input",ParamNode::INSERT);
  iNode->AddChildNode(meshInputNode);

  //create domain node
  //UInt dim = inFile_->GetDim();
  CoupledField::PtrParamNode dNode = dummyXMLNode_->Get("domain",ParamNode::APPEND);
  //dNode->Get("geometryType",ParamNode::APPEND)->SetValue(itoa(dim));

}


CF::UInt MeshFilter::CountUsedEntities(const StdVector<CF::UInt>& entities) {
  const CF::UInt size = entities.GetSize();
  CF::UInt numEntities = 0;
//#pragma omp parallel for reduction(+:numEntities) num_threads(CFS_NUM_THREADS)
  for(CF::UInt inEnt = 0; inEnt < size; inEnt++) {
    if (entities[inEnt] != UnusedEntityNumber) {
      numEntities++;
    }
  }
  return numEntities;
}

void MeshFilter::GetUsedMappedEntities(const shared_ptr<EqnMapSimple>& map,
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
//#pragma omp parallel for num_threads(CFS_NUM_THREADS)
    for (UInt eIter = 0; eIter < size; ++eIter) {
      CF::UInt entityNumber = regEntities[eIter];
      entities[map->GetEntityIndex(entityNumber)] = entityNumber;
    }
  }
}

template<typename T>
void MeshFilter::Node2Cell(Vector<T>& returnVec,
                                      const boost::uuids::uuid& resId,
                                      const Vector<T>& inVec,
                                      const std::vector<QuantityStruct>& interpolData){

  returnVec.Init();
  UInt curE;
  CF::StdVector<UInt> eqns;

  shared_ptr<EqnMapSimple> downMap = resultManager_->GetEqnMap(resId);


  // for every element in the target mesh
  for(UInt i=0;i < interpolData.size();++i){
    QuantityStruct aStru = interpolData[i];
    curE = aStru.trgElemNum;

    //Sum up the contributions from all nodes of element curE
    StdVector<T> curval;
    StdVector<UInt> sEqn;
    downMap->GetEquation(eqns,curE,ExtendedResultInfo::ELEMENT);
    curval.Resize(eqns.GetSize(), 0.0);

    //UInt dim = eqns.GetSize(); //dimension of input data (scalar, 2vector, 3vector)
    for(UInt aNode =0;aNode < aStru.tNNum.GetSize(); ++aNode){
      for(UInt aDof=0;aDof < eqns.GetSize(); aDof++){
    	  Double fac = 1.0/aStru.tNNum.GetSize();
        returnVec[eqns[aDof]] += inVec[aStru.srcEqn[eqns.GetSize()*aNode + aDof]] * fac;
      }
    }
  }
}


void MeshFilter::Cell2Node(Vector<Double>& returnVec,
                          const boost::uuids::uuid& resId,
                          const Vector<Double>& inVec,
                          const std::vector<QuantityStruct>& interpolData,
                          const StdVector<UInt>& nodeNeighbours){

  //perform interpolation
  returnVec.Init();
  CF::Vector<Double> shFnc;
  CF::StdVector<UInt> eqns;
  CF::shared_ptr<ElemShapeMap> eShape;
  shared_ptr<EqnMapSimple> downMap = resultManager_->GetEqnMap(resId);

  for(UInt i=0;i < interpolData.size();++i){
    QuantityStruct aStru = interpolData[i];

    const Elem* curE = trgGrid_->GetElem(aStru.trgElemNum);
    eShape = trgGrid_->GetElemShapeMap(curE,true);

    const CF::StdVector<UInt>& eConn = curE->connect;

    FeH1 * myElem = dynamic_cast<FeH1*>(eShape->GetBaseFE());
    //we assume scalar shape functions
    shFnc.Resize(eConn.GetSize());
    shFnc.Init();
    myElem->GetShFnc(shFnc,aStru.localCoords,curE);

    Double curval = 0.0;
    for(UInt aNode =0;aNode < eConn.GetSize(); ++aNode){
      downMap->GetEquation(eqns,eConn[aNode],ExtendedResultInfo::NODE);
      curval  = eConn.GetSize() * shFnc[aNode]/nodeNeighbours[eConn[aNode]];//* aStru.volume; // We just add up the values

      for(UInt aDof=0;aDof < eqns.GetSize(); aDof++){
        returnVec[eqns[aDof]] += curval * inVec[aStru.srcEqnSingle+aDof];
      }
    }
  }
}

void MeshFilter::NearestNeighbourLight(Vector<Double>& returnVec,
                                      const Vector<Double>& inVec,
                                      const UInt& numEquPerEnt,
                                      const StdVector< StdVector<CF::UInt> >& sourceM,
                                      const StdVector< CF::Matrix<CF::Double> >& targetSourceFactor,
                                      const UInt& maxNumTrgEntities){

  UInt tDim = numEquPerEnt;
  returnVec.Resize(maxNumTrgEntities * tDim);
  returnVec.Init();


//#pragma omp parallel for num_threads(CFS_NUM_THREADS)
  for (UInt i = 0; i < maxNumTrgEntities; i++) {
    CF::UInt targetIndex = i * tDim;
    for (UInt k = 0; k < tDim; k++) {
      returnVec[targetIndex + k] = 0.0;
    }

    StdVector<CF::UInt> sM = sourceM[i];
    UInt patchSize = sM.GetSize();
    CF::Matrix<Double> vals;
    vals.Resize(tDim, patchSize);
    vals.Init();
    UInt k = 0;
    for (CF::UInt j = 0; j < patchSize; j++) {
      CF::UInt sourceIndex = sM[j];
      for(UInt d = 0; d < tDim; ++d){
        vals[d][k] = inVec[tDim * (sourceIndex-1) + d];
      }
      k += 1;
     }
    CF::Matrix<Double> sol;
    sol = vals * targetSourceFactor[i];
    for(UInt l = 0; l < tDim; ++l){
      returnVec[targetIndex + l] = sol[0][l];
    }
  }

}



void MeshFilter::RBFInterpolation(Vector<Double>& returnVec,
                                  const Vector<Double>& inVec,
                                  const UInt& numEquPerEnt,
                                  const StdVector<CF::UInt>& targetSource,
                                  const StdVector<CF::UInt>& targetSourceIndex,
                                  const StdVector< CF::Matrix<Double> >& targetRBFInv,
                                  const StdVector<CF::Double>& targetSourceFactor,
                                  const StdVector<CF::Double>& targetSourceFactor2,
                                  const UInt& maxNumTrgEntities){

  returnVec.Resize(maxNumTrgEntities * numEquPerEnt);
  returnVec.Init();


  if (numEquPerEnt == 1) {
    /***************** SCALAR DATA ****************/
#pragma omp parallel for num_threads(CFS_NUM_THREADS)
    for (UInt i = 0; i < maxNumTrgEntities; i++) {
      CF::Matrix<Double> m = targetRBFInv[i];
      UInt numNN = m.GetNumRows();
      CF::Double R_k = 0.0;
      const CF::UInt jEnd = targetSourceIndex[i + 1];
//      std::cout << targetSourceIndex[i + 1] << std::endl;
      //multiply inverse of local RBF matrix with values
      CF::Matrix<Double> vals;
      vals.Resize(numNN,1);
//      std::cout << numNN << std::endl;
      UInt t = 0;
      for (CF::UInt j = targetSourceIndex[i]; j < jEnd; j++) {
//        std::cout << j << std::endl;
         vals[t][0] = inVec[targetSource[j]];
         t += 1;
      }
      t = 0;
      CF::Matrix<Double> coefVec;


      coefVec = targetRBFInv[i] * vals;
      CF::Double factor2;
      for (CF::UInt j = targetSourceIndex[i]; j < jEnd; j++) {
         factor2 = targetSourceFactor2[j];
         R_k += coefVec[t][0] * factor2;
         t += 1;
      }

      CF::Double sum = 0.0;
      for (CF::UInt j = targetSourceIndex[i]; j < jEnd; j++) {
        sum += R_k * targetSourceFactor[j];
      }

      returnVec[i] = sum;
    }
  } else {
    /***************** 2/3D VECTOR DATA ****************/
#pragma omp parallel for num_threads(CFS_NUM_THREADS)
    for (UInt i = 0; i < maxNumTrgEntities; i++) {
      CF::Vector<Double> R_k;
      R_k.Resize(numEquPerEnt);
      R_k.Init();

      CF::Matrix<Double> m = targetRBFInv[i];
      UInt numNN = m.GetNumRows();

      CF::UInt targetIndex = i * numEquPerEnt;
      for (UInt k = 0; k < numEquPerEnt; k++) {
          returnVec[targetIndex + k] = 0.0;
      }

      CF::Matrix<Double> vals;
      vals.Resize(numNN, numEquPerEnt);
      UInt t = 0;

      const CF::UInt jEnd = targetSourceIndex[i + 1];
      for (CF::UInt j = targetSourceIndex[i]; j < jEnd; j++) {
        CF::UInt sourceIndex = targetSource[j] * numEquPerEnt;
        for(CF::UInt k = 0; k < numEquPerEnt; ++k){
           vals[t][k] = inVec[sourceIndex + k];
        }
        t += 1;
      }

      CF::Matrix<Double> coefVec;
      coefVec = targetRBFInv[i] * vals;
      t = 0;
      for (CF::UInt j = targetSourceIndex[i]; j < jEnd; j++) {
         for(UInt k = 0; k < numEquPerEnt; ++k){
             R_k[k] += coefVec[t][k] * targetSourceFactor2[j];
         }
         t += 1;
      }
      CF::Vector<Double> sum;
      sum.Resize(numEquPerEnt);
      sum.Init();
      for (CF::UInt j = targetSourceIndex[i]; j < jEnd; j++) {
        for(UInt k = 0; k < numEquPerEnt; ++k){
          sum[k] += R_k[k] * targetSourceFactor[j];
        }
      }

//      if(i == 372){
//        std::cout << numNN << std::endl;
//        std::cout << vals << std::endl;
//        std::cout << sum << std::endl;
//      }

      for(UInt k = 0; k < numEquPerEnt; ++k){
        returnVec[targetIndex + k] = sum[k];
      }
    }
  }
}



void MeshFilter::CalcLocRBFInv(CF::Matrix<Double>& matr,
                                  const StdVector< CF::Vector<CF::Double> >& neighbors,
                                  const Double& alpha,
                                  const UInt numNN,
                                  Grid* grid){

  matr.Resize(numNN,numNN);
  Double rNN; //distance between two src points
  for (UInt i = 0; i < numNN; ++i){
    for (UInt j = 0; j < numNN; ++j){
      if (grid->GetDim() == 3){
        rNN = sqrt(pow(neighbors[i][0]-neighbors[j][0],2.0) + pow(neighbors[i][1]-neighbors[j][1],2.0) + pow(neighbors[i][2]-neighbors[j][2],2.0));
      }else{
        rNN = sqrt(pow(neighbors[i][0]-neighbors[j][0],2.0) + pow(neighbors[i][1]-neighbors[j][1],2.0));
      }
      matr[i][j] = pow(1.0 - rNN/alpha, 2.0);
    }
  }
  // now we have to invert Aloc and multiply it with the according value-coloumn
  matr.Invert_Lapack();
}





bool MeshFilter::CalcLocCurl(CF::Matrix<Double>& derivCoefVec,
                            const CF::Vector<CF::Double>& globPoint,
                            const CF::Double& maxDist,
                            const StdVector<CF::Double>& l2Distances,
                            const StdVector< CF::Vector<CF::Double> >& neighbors,
                            const UInt& numNeighbors, //TODO WARUM?
                            const UInt& numEquPerEnt,
                            Grid* grid,
                            const Double epsScal,
							const Double betaScal,
							const Double kScal,
                            const bool logEps){
  bool c = false;

  derivCoefVec.Resize(numNeighbors,1);
  CF::Matrix<Double> ALoc;
  ALoc.Resize(numNeighbors,numNeighbors);
  CF::Matrix<Double> derivVec; //Vector of RBF derivatives evaluated at srcPoints

  Double rNNSquared = 0.0; //distance between two src points
  Double eps = epsScal / maxDist;
  UInt iter = 0;
  Double upperBound = 1E-14;
  Double lowerBound = 1E-17;
  while( !c){
      for (UInt i = 0; i < numNeighbors; ++i){
        for (UInt j = 0; j < numNeighbors; ++j){
          if (grid->GetDim() == 3){
            rNNSquared = pow(neighbors[i][0]-neighbors[j][0],2.0) + pow(neighbors[i][1]-neighbors[j][1],2.0) + pow(neighbors[i][2]-neighbors[j][2],2.0);
          }else{

            rNNSquared = pow(neighbors[i][0]-neighbors[j][0],2.0) + pow(neighbors[i][1]-neighbors[j][1],2.0);
          }

          ALoc[i][j] = exp(-(eps*eps * rNNSquared))+betaScal*rNNSquared+kScal;
        }
        switch( numEquPerEnt ){
        case 1:
        case 2:
          if (grid->GetDim() == 2){
            derivVec.Resize(numNeighbors,2);
            if (l2Distances[i] == 0) {
              derivVec[i][0] = 0.0;
              derivVec[i][1] = 0.0;
            }else{
                derivVec[i][0] =  exp(-(eps*eps * l2Distances[i] * l2Distances[i])) * (-2.0 * eps*eps * (globPoint[0] - neighbors[i][0])) + betaScal*(globPoint[0] - neighbors[i][0])/l2Distances[i];
                derivVec[i][1] =  exp(-(eps*eps * l2Distances[i] * l2Distances[i])) * (-2.0 * eps*eps * (globPoint[1] - neighbors[i][1])) + betaScal*(globPoint[1] - neighbors[i][1])/l2Distances[i];
            }
          }else{
            EXCEPTION("2D mesh and 3D-values!")
          }
          break;
        case 3:
          if (grid->GetDim() == 3){
            derivVec.Resize(numNeighbors,3);
            if (l2Distances[i] == 0) {
              derivVec[i][0] = 0.0;
              derivVec[i][1] = 0.0;
              derivVec[i][2] = 0.0;
            }else{
              derivVec[i][0] =  exp(-(eps*eps * l2Distances[i] * l2Distances[i])) * (-2.0 * eps*eps * (globPoint[0] - neighbors[i][0])) + betaScal*(globPoint[0] - neighbors[i][0])/l2Distances[i];
              derivVec[i][1] =  exp(-(eps*eps * l2Distances[i] * l2Distances[i])) * (-2.0 * eps*eps * (globPoint[1] - neighbors[i][1])) + betaScal*(globPoint[1] - neighbors[i][1])/l2Distances[i];
              derivVec[i][2] =  exp(-(eps*eps * l2Distances[i] * l2Distances[i])) * (-2.0 * eps*eps * (globPoint[2] - neighbors[i][2])) + betaScal*(globPoint[2] - neighbors[i][2])/l2Distances[i];
            }
          }else{
            EXCEPTION("3D mesh and 2D-values!")
          }      break;
        }
      }

      Double k; //inverse of condition number
      int inf;
      ALoc.Invert_Lapack(k, inf);
      if(k < upperBound && k > lowerBound && inf==0) c = true;
      else if(k > upperBound && inf==0) eps = eps / 2.0;
      else if(k < lowerBound || inf!=0) eps = eps * 2.0;
      if(iter > 6) upperBound = upperBound * 2.0;
      ++iter;
    }

    //log min, max distance and epsilon
    if(logEps){
      //find max and min distance
      double min = l2Distances[0], max = l2Distances[0];
      for(UInt i=0; i < l2Distances.GetSize(); ++i){
        if(min>l2Distances[i]) min = l2Distances[i];
        if(max<l2Distances[i]) max = l2Distances[i];
      }
      std::cout<<min<<"\t"<<max<<"\t"<<eps<<std::endl;
    }

    // coefficient matrix (coloumn nr. corresponding to the spatial dimension)
    derivCoefVec = ALoc * derivVec;

      return true;

}




bool MeshFilter::CalcLocGradient(CF::Matrix<Double>& derivCoefVec,
                                const CF::Vector<CF::Double>& globPoint,
                                const CF::Double& maxDist,
                                const StdVector<CF::Double>& l2Distances,
                                const StdVector< CF::Vector<CF::Double> >& neighbors,
                                const UInt& numNeighbors,
                                const UInt& numEquPerEnt,
                                Grid* grid,
                                const Double epsScal,
								const Double betaScale,
								const Double kScale,
                                const bool logEps){

  derivCoefVec.Resize(numNeighbors,1);
  CF::Matrix<Double> ALoc;
  ALoc.Resize(numNeighbors,numNeighbors);
  CF::Matrix<Double> derivVec; //Vector of RBF derivatives evaluated at srcPoints

  bool c = false;

  Double rNNSquared = 0.0; //distance between two src points
  Double eps = epsScal / maxDist;
  UInt iter = 0;
  Double upperBound = 1E-14;
  Double lowerBound = 1E-17;
  while( !c){
    for (UInt i = 0; i < numNeighbors; ++i){
      for (UInt j = 0; j < numNeighbors; ++j){
        if (grid->GetDim() == 3){
          rNNSquared = pow(neighbors[i][0]-neighbors[j][0],2.0) + pow(neighbors[i][1]-neighbors[j][1],2.0) + pow(neighbors[i][2]-neighbors[j][2],2.0);
        }else{

          rNNSquared = pow(neighbors[i][0]-neighbors[j][0],2.0) + pow(neighbors[i][1]-neighbors[j][1],2.0);
        }
        //ALoc[i][j] = exp(-(eps*eps * rNNSquared));
        ALoc[i][j] = exp(-(eps*eps * rNNSquared))+betaScale*rNNSquared+kScale;
      }
      switch( numEquPerEnt ){
      case 1:
        if (grid->GetDim() == 2){
          derivVec.Resize(numNeighbors,2);
          if (l2Distances[i] == 0) {
            derivVec[i][0] = 0.0;
            derivVec[i][1] = 0.0;
          }else{
            derivVec[i][0] =  exp(-(eps*eps * l2Distances[i] * l2Distances[i])) * (-2.0 * eps*eps * (globPoint[0] - neighbors[i][0])) + betaScale*(globPoint[0] - neighbors[i][0])/l2Distances[i];
            derivVec[i][1] =  exp(-(eps*eps * l2Distances[i] * l2Distances[i])) * (-2.0 * eps*eps * (globPoint[1] - neighbors[i][1])) + betaScale*(globPoint[1] - neighbors[i][1])/l2Distances[i];
          }
        }else{
          derivVec.Resize(numNeighbors,3);
          if (l2Distances[i] == 0) {
            derivVec[i][0] = 0.0;
            derivVec[i][1] = 0.0;
            derivVec[i][2] = 0.0;
          }else{
            derivVec[i][0] =  exp(-(eps*eps * l2Distances[i] * l2Distances[i])) * (-2.0 * eps*eps * (globPoint[0] - neighbors[i][0])) + betaScale*(globPoint[0] - neighbors[i][0])/l2Distances[i];
            derivVec[i][1] =  exp(-(eps*eps * l2Distances[i] * l2Distances[i])) * (-2.0 * eps*eps * (globPoint[1] - neighbors[i][1])) + betaScale*(globPoint[1] - neighbors[i][1])/l2Distances[i];
            derivVec[i][2] =  exp(-(eps*eps * l2Distances[i] * l2Distances[i])) * (-2.0 * eps*eps * (globPoint[2] - neighbors[i][2])) + betaScale*(globPoint[2] - neighbors[i][2])/l2Distances[i];
          }
        }
        break;
      case 2:
        EXCEPTION("Gradient of Vector not defined in this context!");
        break;
      case 3:
        EXCEPTION("Gradient of Vector not defined in this context!");
        break;
      }
    }


    Double k; //inverse of condition number
    int inf;
    ALoc.Invert_Lapack(k, inf);
    if(inf != 0){
    	EXCEPTION ("Error during lapack inversion, Error integer = "<<inf);
    }
    if(k < upperBound && k > lowerBound && inf==0) c = true;
    else if(k > upperBound && inf==0) eps = eps / 2.0;
    else if(k < lowerBound || inf!=0) eps = eps * 2.0;
    if(iter > 6) upperBound = upperBound * 2.0;
    ++iter;
  }


  //log min, max distance and epsilon
  if(logEps){
    //find max and min distance
    double min = l2Distances[0], max = l2Distances[0];
    for(UInt i=0; i < l2Distances.GetSize(); ++i){
      if(min>l2Distances[i]) min = l2Distances[i];
      if(max<l2Distances[i]) max = l2Distances[i];
    }
    std::cout<<min<<"\t"<<max<<"\t"<<eps<<std::endl;
  }

  // coefficient matrix (coloumn nr. corresponding to the spatial dimension)
  derivCoefVec = ALoc * derivVec;


  return true;
}



bool MeshFilter::CalcLocDivergence(CF::Matrix<Double>& derivCoefVec,
                                  const CF::Vector<CF::Double>& globPoint,
                                  const CF::Double& maxDist,
                                  const StdVector<CF::Double>& l2Distances,
                                  const StdVector< CF::Vector<CF::Double> >& neighbors,
                                  const UInt& numNeighbors,
                                  const UInt& numEquPerEnt,
                                  Grid* grid,
                                  const Double epsScal,
								  const Double betaScal,
								  const Double kScal,
                                  const bool logEps){

  Vector<Double> eigenVals;
  bool c = false;

  CF::Matrix<Double> ALoc;
  ALoc.Resize(numNeighbors,numNeighbors);
  CF::Matrix<Double> vals;
  CF::Matrix<Double> derivVec; //Vector of RBF derivatives evaluated at srcPoints

  Double rNNSquared = 0.0; //distance between two src points
  Double eps = epsScal / maxDist;
  UInt iter = 0;
  Double upperBound = 1E-14;
  Double lowerBound = 1E-17;
  while( !c){
    for (UInt i = 0; i < numNeighbors; ++i){
      for (UInt j = 0; j < numNeighbors; ++j){
        if (grid->GetDim() == 3){
          rNNSquared = pow(neighbors[i][0]-neighbors[j][0],2.0) + pow(neighbors[i][1]-neighbors[j][1],2.0) + pow(neighbors[i][2]-neighbors[j][2],2.0);
        }else{

          rNNSquared = pow(neighbors[i][0]-neighbors[j][0],2.0) + pow(neighbors[i][1]-neighbors[j][1],2.0);
        }
        ALoc[i][j] = exp(-(eps*eps * rNNSquared))+betaScal*rNNSquared+kScal;
      }
      switch( numEquPerEnt ){
      case 1:
        //should already be caught
        EXCEPTION("Divergence of a scalar field!");
        break;
      case 2:
        if (grid->GetDim() == 2){
        derivVec.Resize(numNeighbors,2);
        if (l2Distances[i] == 0) {
          derivVec[i][0] = 0.0;
          derivVec[i][1] = 0.0;
        }else{
            derivVec[i][0] =  exp(-(eps*eps * l2Distances[i] * l2Distances[i])) * (-2.0 * eps*eps * (globPoint[0] - neighbors[i][0])) + betaScal*(globPoint[0] - neighbors[i][0])/l2Distances[i];
            derivVec[i][1] =  exp(-(eps*eps * l2Distances[i] * l2Distances[i])) * (-2.0 * eps*eps * (globPoint[1] - neighbors[i][1])) + betaScal*(globPoint[1] - neighbors[i][1])/l2Distances[i];

        }
        }else{
          EXCEPTION("2D values and 3D mesh!");
        }
        break;
      case 3:
        if (grid->GetDim() == 3){
        derivVec.Resize(numNeighbors,3);
        if (l2Distances[i] == 0) {
          derivVec[i][0] = 0.0;
          derivVec[i][1] = 0.0;
          derivVec[i][2] = 0.0;
        }else{
          derivVec[i][0] =  exp(-(eps*eps * l2Distances[i] * l2Distances[i])) * (-2.0 * eps*eps * (globPoint[0] - neighbors[i][0])) + betaScal*(globPoint[0] - neighbors[i][0])/l2Distances[i];
          derivVec[i][1] =  exp(-(eps*eps * l2Distances[i] * l2Distances[i])) * (-2.0 * eps*eps * (globPoint[1] - neighbors[i][1])) + betaScal*(globPoint[1] - neighbors[i][1])/l2Distances[i];
          derivVec[i][2] =  exp(-(eps*eps * l2Distances[i] * l2Distances[i])) * (-2.0 * eps*eps * (globPoint[2] - neighbors[i][2])) + betaScal*(globPoint[2] - neighbors[i][2])/l2Distances[i];
        }
        }else{
           WARN("DivergenceDifferentiator.cc : Treat 3D values as 2D values, due to a 2D mesh!")
           derivVec.Resize(numNeighbors,2);
           if (l2Distances[i] == 0) {
             derivVec[i][0] = 0.0;
             derivVec[i][1] = 0.0;
           }else{
             derivVec[i][0] =  exp(-(eps*eps * l2Distances[i] * l2Distances[i])) * (-2.0 * eps*eps * (globPoint[0] - neighbors[i][0])) + betaScal*(globPoint[0] - neighbors[i][0])/l2Distances[i];
             derivVec[i][1] =  exp(-(eps*eps * l2Distances[i] * l2Distances[i])) * (-2.0 * eps*eps * (globPoint[1] - neighbors[i][1])) + betaScal*(globPoint[1] - neighbors[i][1])/l2Distances[i];
           }
        }
        break;
      }
   }

    Double k; //inverse of condition number
    int inf;
    ALoc.Invert_Lapack(k, inf);
    if(k < upperBound && k > lowerBound && inf==0) c = true;
    else if(k > upperBound && inf==0) eps = eps / 2.0;
    else if(k < lowerBound || inf!=0) eps = eps * 2.0;
    if(iter > 6) upperBound = upperBound * 2.0;
    ++iter;
  }


  //log min, max distance and epsilon
  if(logEps){
    //find max and min distance
    double min = l2Distances[0], max = l2Distances[0];
    for(UInt i=0; i < l2Distances.GetSize(); ++i){
      if(min>l2Distances[i]) min = l2Distances[i];
      if(max<l2Distances[i]) max = l2Distances[i];
    }
    std::cout<<min<<"\t"<<max<<"\t"<<eps<<std::endl;
  }

  // coefficient matrix (coloumn nr. corresponding to the spatial dimension)
  derivCoefVec = ALoc * derivVec;

  /*
      Double s1 = 0.0;
      Double s2 = 0.0;
      for(UInt i = 0; i < derivCoefVec.GetNumRows(); ++i ){
        s1 += derivCoefVec[i][0];
        s2 += derivCoefVec[i][1];
      }

      Double t = eps / 0.0;
      if( (fabs(s1)> t) || (fabs(s2) > t) ){
      ret = false;
      }else{ ret = true;}
      return ret;
  */
      return true;

}



void MeshFilter::CalcCurl(Vector<Double>& returnVec,
                                  const Vector<Double>& inVec,
                                  const UInt& numEquPerEnt,
                                  const StdVector< StdVector<CF::UInt> >& sourceM,
                                  const StdVector< CF::Matrix<CF::Double> >& targetSourceFactor,
                                  const UInt& maxNumTrgEntities,
                                  const UInt& gridDim){

  returnVec.Resize(maxNumTrgEntities * numEquPerEnt);
  returnVec.Init();

//#pragma omp parallel for num_threads(CFS_NUM_THREADS)
  for (UInt i = 0; i < maxNumTrgEntities; i++) {
    CF::UInt targetIndex = i * numEquPerEnt;
    for (UInt k = 0; k < numEquPerEnt; k++) {
      returnVec[targetIndex + k] = 0.0;
    }

    StdVector<CF::UInt> sM = sourceM[i];
    UInt patchSize = sM.GetSize();

    CF::Matrix<Double> vals;
    vals.Resize(numEquPerEnt, patchSize);
    vals.Init();
    UInt k = 0;
    for (CF::UInt j =0; j < patchSize; j++) {
      CF::UInt sourceIndex = sM[j];
      for(UInt d = 0; d < numEquPerEnt; ++d){
        vals[d][k] = inVec[numEquPerEnt * (sourceIndex-1) + d];
      }
      k += 1;
     }
    CF::Matrix<Double> tempMat;
    tempMat = vals * targetSourceFactor[i];

    StdVector<CF::Double> vec;
    vec.Resize(3);
    vec.Init(0.0);
    //now we have to combine the tempvec-entries in order to obtain the curl
    // TODO add here the variants
    if(gridDim == 2){
      // in 2D the vector actually points in z-direction
      // but here we define it to point in x !!!!
      vec[0] = tempMat[1][0] - tempMat[0][1];
    }else{
      vec[0] = tempMat[2][1] - tempMat[1][2];
      vec[1] = tempMat[0][2] - tempMat[2][0];
      vec[2] = tempMat[1][0] - tempMat[0][1];
    }
    for (UInt k = 0; k < numEquPerEnt; k++) {
      returnVec[targetIndex + k] = vec[k];
    }

  }

}


void MeshFilter::CalcDivergence(Vector<Double>& returnVec,
                                      const Vector<Double>& inVec,
                                      const UInt& numEquPerEnt,
                                      const StdVector< StdVector<CF::UInt> >& sourceM,
                                      const StdVector< CF::Matrix<CF::Double> >& targetSourceFactor,
                                      const UInt& maxNumTrgEntities){

  returnVec.Resize(maxNumTrgEntities);
  returnVec.Init();

//#pragma omp parallel for num_threads(CFS_NUM_THREADS)
  for (UInt i = 0; i < maxNumTrgEntities; i++) {
    CF::UInt targetIndex = i;
    returnVec[targetIndex] = 0.0;

    StdVector<CF::UInt> sM = sourceM[targetIndex];

    UInt patchSize = sM.GetSize();
    CF::Matrix<Double> vals;
    vals.Resize(numEquPerEnt, patchSize);
    vals.InitValue(0.0);

    UInt k = 0;
    for (CF::UInt j = 0; j < patchSize; j++) {
      CF::UInt sourceIndex = sM[j];
      for(UInt d = 0; d < numEquPerEnt; ++d){
        vals[d][k] = inVec[numEquPerEnt * (sourceIndex-1) + d];
      }
      k += 1;
     }
    CF::Matrix<Double> tempMat;
    tempMat = vals * targetSourceFactor[targetIndex];
    returnVec[targetIndex] = tempMat.Trace();

  }

}

int MeshFilter::Index2Voigt(          const UInt& dx1,
                                      const UInt& dx2,
                                      const UInt& dim) {
  if(dx1>=dim || dx2>=dim)
    EXCEPTION("Index does not exist!")
  if (dim == 3) {
    int indexT[3][3] = {{0,5,4},{5,1,3},{4,3,2}};
    return indexT[dx1][dx2];
  } else if (dim == 2){
    int indexT[2][2] = {{0,2},{2,1}} ;
    return indexT[dx1][dx2];
  }
  return dx1;

}


void MeshFilter::CalcTensorDivergence(Vector<Double>& returnVec,
                                      const Vector<Double>& inVec,
                                      const UInt& numEquPerEnt,
                                      const StdVector< StdVector<CF::UInt> >& sourceM,
                                      const StdVector< CF::Matrix<CF::Double> >& targetSourceFactor,
                                      const UInt& maxNumTrgEntities){

  returnVec.Resize(maxNumTrgEntities * numEquPerEnt);
  returnVec.Init();

  UInt numbTensorEntries = numEquPerEnt + 1;
  if (numEquPerEnt==3)numbTensorEntries = numEquPerEnt + 3;

//#pragma omp parallel for num_threads(CFS_NUM_THREADS)
  for (UInt targetIndex = 0; targetIndex < maxNumTrgEntities; targetIndex++) {

    StdVector<CF::UInt> sM = sourceM[targetIndex];

    UInt patchSize = sM.GetSize();
    CF::Matrix<Double> vals;
    vals.Resize(numEquPerEnt, patchSize);
    vals.InitValue(0.0);

    for(UInt tensorRow = 0; tensorRow < numEquPerEnt; ++tensorRow){

      for (CF::UInt j = 0; j < patchSize; j++) {
        CF::UInt sourceIndex = sM[j];
        for(UInt d = 0; d < numEquPerEnt; ++d){
          vals[d][j] = inVec[numbTensorEntries * (sourceIndex-1) + Index2Voigt(tensorRow,d,numEquPerEnt)];

        }
       }
      CF::Matrix<Double> tempMat;
      tempMat = vals * targetSourceFactor[targetIndex];
      returnVec[numEquPerEnt * targetIndex + tensorRow] = tempMat.Trace();
    }
  }

}

void MeshFilter::CalcGradient(Vector<Double>& returnVec,
                                  const Vector<Double>& inVec,
                                  const UInt& numEquPerEnt,
                                  const StdVector< StdVector<CF::UInt> >& sourceM,
                                  const StdVector< CF::Matrix<CF::Double> >& targetSourceFactor,
                                  const UInt& maxNumTrgEntities,
                                  const UInt& tDim){


  returnVec.Resize(maxNumTrgEntities * tDim);
  returnVec.Init();


//#pragma  omp parallel for num_threads(CFS_NUM_THREADS)
  for (UInt i = 0; i < maxNumTrgEntities; i++) {
    CF::UInt targetIndex = i * tDim;
    for (UInt k = 0; k < tDim; k++) {
      returnVec[targetIndex + k] = 0.0;
    }

    StdVector<CF::UInt> sM = sourceM[i];
    UInt patchSize = sM.GetSize();

    CF::Matrix<Double> vals;
    vals.Resize(1, patchSize);
    vals.Init();
    UInt k = 0;
    for (CF::UInt j = 0; j < patchSize; j++) {
      CF::UInt sourceIndex = sM[j];
      vals[0][k] = inVec[sourceIndex - 1];
      k += 1;
     }

    CF::Matrix<Double> sol;
    sol = vals * targetSourceFactor[i];

    for(UInt l = 0; l < tDim; ++l){
      returnVec[targetIndex + l] = sol[0][l];
    }
  }

}

template void MeshFilter::Node2Cell(Vector<Double>& returnVec,const boost::uuids::uuid& resId,const Vector<Double>& inVec,const std::vector<QuantityStruct>& interpolData);
template void MeshFilter::Node2Cell(Vector<Complex>& returnVec,const boost::uuids::uuid& resId,const Vector<Complex>& inVec,const std::vector<QuantityStruct>& interpolData);

}


