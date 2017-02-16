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
#include "cfsdat/Utils/KNNSearch.hh"
#include "FeBasis/H1/H1Elems.hh"


namespace CFSDat{

MeshFilter::MeshFilter(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
                      :BaseMeshFilterType(numWorkers,config,resMan){

  // Aeroacoustic-Filters choose their input and output results independently because
  // there we might have two inputs (i.e. velocity and vorticity)
  if (params_-> Get("type")->As<std::string>() != "AeroacousticSource_LambVector" &&
      params_-> Get("type")->As<std::string>() != "AeroacousticSource_LighthillSourceVector" &&
      params_-> Get("type")->As<std::string>() != "AeroacousticSource_LighthillSourceTerm"){
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

void MeshFilter::FinishInit(){
  //first go up
  CF::StdVector< str1::shared_ptr<BaseFilter> >::iterator srcIter =  sources_.Begin();
  for(; srcIter != sources_.End() ; srcIter++){
    // should we check here anything for success?
    (*srcIter)->FinishInit();
  }
  PrepareCalculation();
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



void MeshFilter::Node2Cell(Vector<Double>& returnVec,
                                      const boost::uuids::uuid& resId,
                                      const Vector<Double>& inVec,
                                      const std::vector<QuantityStruct>& interpolData){

  returnVec.Init();
  //current element
  UInt curE;
  CF::StdVector<UInt> eqns;

  str1::shared_ptr<EqnMapSimple> downMap = resultManager_->GetResultAdapter(resId)->mapping;


  // for every element in the target mesh
  for(UInt i=0;i < interpolData.size();++i){
    QuantityStruct aStru = interpolData[i];
    curE = aStru.trgElemNum;

    //Sum up the contributions from all nodes of element curE
    StdVector<Double> curval;
    StdVector<UInt> sEqn;
    downMap->GetEquation(eqns,curE,ExtendedResultInfo::ELEMENT);
    curval.Resize(eqns.GetSize(), 0.0);

    //UInt dim = eqns.GetSize(); //dimension of input data (scalar, 2vector, 3vector)
    for(UInt aNode =0;aNode < aStru.tNNum.GetSize(); ++aNode){
      for(UInt aDof=0;aDof < eqns.GetSize(); aDof++){
        returnVec[eqns[aDof]] += inVec[aStru.srcEqn[eqns.GetSize()*aNode + aDof]] / aStru.tNNum.GetSize();
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
  str1::shared_ptr<EqnMapSimple> downMap = resultManager_->GetResultAdapter(resId)->mapping;

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



void MeshFilter::NearestNeighbourInterpolation(Vector<Double>& returnVec,
                                        const CF::StdVector< CF::Vector<Double> >&  scatteredData,
                                        CF::Vector<Double>& vec,
                                        const str1::shared_ptr<EqnMapSimple>& downMap,
                                        const CF::StdVector< CF::Vector<double> >& srcCoords,
                                        const CF::StdVector< CF::Vector<double> >& trgCoords,
                                        const std::vector<QuantityStruct>& interpolationData,
                                        Grid* grid,
                                        const UInt& dim,
                                        const UInt& numNeighbors,
                                        const UInt& p){

    // Object for nearest neighbor-searches and bringing the data into the correct form for CGAL search
    KNNSearch Tree;
    Tree.ReadScatteredData_Interpolation(srcCoords, dim, grid, scatteredData);

    CF::StdVector<bool> nodeCheck;
    nodeCheck.Resize(grid->GetNumNodes(ALL_REGIONS));
    nodeCheck.Init(false);

    // loop over all elements and over every node of each element
    for(UInt i=0;i < interpolationData.size();++i){
      UInt nodeIter = 0;
      QuantityStruct aStru = interpolationData[i];
      const Elem* curE = grid->GetElem(aStru.trgElemNum);
      CF::shared_ptr<ElemShapeMap> eShape = grid->GetElemShapeMap(curE,true);
      //eConn gives us the node numbers of current element curE
      const CF::StdVector<UInt>& eConn = curE->connect;
      //std::cout<<"curE"<<curE<<std::endl;
      //std::cout<<"eConn"<<eConn<<std::endl;

      //get the global coordinates of the nodes of element curE
      CF::Matrix<Double> globalCoords;
      grid->GetElemNodesCoord(globalCoords, eConn, true);
      //std::cout<<"globalCoords"<<globalCoords<<std::endl;
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

          for (UInt dof = 0; dof < dim + 1; dof++){
            vec[dof] = 0.0;
          }

          // at that point we can start obtaining the nearest neighbours for every aNode
          Tree.KNN_CGAL_Interpolation(globPoint, neighbors, l2dists, vectors, numNeighbors);

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
              Vector<Double> sum(dim+1);
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
                 w = std::pow((R-d)/(R*d), p);
                }
                weights += w;

                for(UInt dof=0; dof < dim + 1; dof++)
                {
                  sum[dof] += vectors[j][dof] * w;
                }
              }

              for(UInt dof=0; dof < dim + 1; dof++)
              {
                vec[dof] = sum[dof] / weights;
              }
              CF::StdVector<UInt> eqns;
              eqns.Init(0);
              //get the equation map for the nodes in eConn, in order to insert the
              //interpolation in the correct position in the result vector
              //downMap->GetEquation(eqns,eConn[aNode],ExtendedResultInfo::NODE);
              downMap->GetEquation(eqns,eConn[aNode],ExtendedResultInfo::NODE);
              //if scalar input values of scattered data->inDim_=0
              //if it is a two-dimensional vector->inDim_=1 and inDim_=2 for 3d-vector
              //std::cout<<"aNode"<<aNode<<std::endl;
              //std::cout<<"eConn[aNode]"<<eConn[aNode]<<std::endl;
              //std::cout<<"eConn"<<eConn<<std::endl;
              //std::cout<<"eqns"<<eqns<<std::endl;
              //std::cout<<"returnVec.GetSize()"<<returnVec.GetSize()<<std::endl;
              //std::cout<<"vec"<<vec<<std::endl;
              for(UInt aDof = 0; aDof < dim+1; aDof++){
                returnVec[eqns[aDof]] = vec[aDof];
              }
              //iterator for the already computed nodes
              nodeIter++;

          }// if nodeMatch == false
      }// for aNode
    }// for element

}


void MeshFilter::RBFInterpolation(Vector<Double>& returnVec,
                                  const CF::StdVector< CF::Vector<Double> >&  scatteredData,
                                  CF::Vector<Double>& vec,
                                  const str1::shared_ptr<EqnMapSimple>& downMap,
                                  const CF::StdVector< CF::Vector<double> >& srcCoords,
                                  const CF::StdVector< CF::Vector<double> >& trgCoords,
                                  const std::vector<QuantityStruct>& interpolationData,
                                  Grid* grid,
                                  const UInt& dim,
                                  const UInt& p,
                                  const bool noSlip){

  // find out the dimension of the mesh (2D, 3D) and assign the appropriate number of nearest neighbours
  // nr. of nearest neighbours for the nodal RBF basis functions
  UInt numNN = (grid->GetDim() == 2) ? 13 : 18 ;
  //nr. of nearest neighbours needed for the weight-calculation (must be smaller than numNN !!)
  UInt numNW = (grid->GetDim() == 3) ? 10 : 7 ;

  //Object for nearest neighbor-searches and bringing the data into the correct form for CGAL search
  KNNSearch Tree;
  Tree.ReadScatteredData_Interpolation(srcCoords, dim, grid, scatteredData);



  CF::StdVector<bool> nodeCheck;
  nodeCheck.Resize(grid->GetNumNodes(ALL_REGIONS));
  nodeCheck.Init(false);

  /* We want to automatically set velocity at a certain region named "wall"
   * to be zero (no-slip-condition). If tag <noSlipWall/> is set in xml-scheme,
   * we  skip the computation of wall-nodes so the value is zero because
   * returnVec was initialized with 0.0
   */
  if(noSlip == true){
    CF::StdVector<UInt> wallNodes;
    grid->GetNodesByName(wallNodes, "wall");

    CoupledField::StdVector<unsigned int>::iterator nIter = wallNodes.Begin();
    //std::set<uint>::iterator nIter = wallNodes.Begin();
    for(;nIter != wallNodes.End(); ++nIter){
      // -1 because nodeNumbers start with 1
      nodeCheck[*nIter - 1] = true;
    }
  }

  CF::Vector<Double> R_k;
  R_k.Resize(dim+1);

  // loop over all elements and over every node of each element
  for(UInt i=0;i < interpolationData.size();++i){

    QuantityStruct aStru = interpolationData[i];
    const Elem* curE = grid->GetElem(aStru.trgElemNum);
    CF::shared_ptr<ElemShapeMap> eShape = grid->GetElemShapeMap(curE,true);
    //eConn gives us the node numbers of current element curE
    const CF::StdVector<UInt>& eConn = curE->connect;

    //get the global coordinates of the nodes of element curE
    CF::Matrix<Double> globalCoords;
    grid->GetElemNodesCoord(globalCoords, eConn, true);
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

        for (UInt dof = 0; dof < dim + 1; dof++){
          vec[dof] = 0.0;
        }

        // at that point we can start obtaining the nearest neighbours for every aNode
        Tree.KNN_CGAL_Interpolation(globPoint, neighbors, l2dists, vectors, numNN);

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
        CalcLocRBFCoefs(coefVec, globPoint, neighbors, l2dists, vectors, numNN, alpha, dim);

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
              W_denom += pow(temp / (r_Wk * r_k), Double(p));
            }
          }
          r_k = l2dists[srcIter];
          temp = (r_Wk - r_k);
          if (temp < 0.0){temp = 0.0;}
          if (r_k == 0.0){
            W_k_bar = 1.0 / W_denom;
          }else{
            W_k_bar = pow(temp / (r_Wk * r_k), Double(p)) / W_denom;
          }

          //RBF evaluation
          //just a clear-up of of R_k
          for (UInt dof = 0; dof < dim + 1; dof++){
            R_k[dof] = 0.0;
          }

          for(UInt dof=0; dof < dim + 1; dof++){
            for (UInt m = 0; m < numNN; ++m){
              R_k[dof] += coefVec[m][dof] * pow(1.0 - l2dists[m]/alpha, 2.0);
            }
          }

          for(UInt dof=0; dof < dim + 1; dof++){
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
      }// if nodeMatch == false


    }// for aNode

  }// for element
  //}// pragma parallel


}





void MeshFilter::CalcLocRBFCoefs(CF::Matrix<Double>& coefVec,
                                  const CF::Vector<Double>& globPoint,
                                  const CF::StdVector< Vector<Double> >& neighbors,
                                  const CF::StdVector< Double >& l2Distances,
                                  const CF::StdVector< Vector<Double> >& vectors,
                                  const UInt& numNN,
                                  const Double& alpha,
                                  const UInt& inDim){

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
    switch( inDim ){
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

void MeshFilter::CalcLocCurl(CF::Matrix<Double>& vec,
                                                  const CF::Vector<Double>& globPoint,
                                                  const CF::StdVector< Vector<Double> >& neighbors,
                                                  const  CF::StdVector< Double >& l2Distances,
                                                  const  CF::StdVector< Vector<Double> >& vectors,
                                                  const  UInt numNN,
                                                  const  Double alpha,
                                                  const UInt dim,
                                                  Grid* grid){
  CF::Matrix<Double> derivCoefVec;
  CF::Matrix<Double> ALoc;
  ALoc.Resize(numNN,numNN);
  CF::Matrix<Double> vals;
  CF::Matrix<Double> derivVec; //Vector of RBF derivatives evaluated at srcPoints

  Double rNN; //distance between two src points
  for (UInt i = 0; i < numNN; ++i){
    for (UInt j = 0; j < numNN; ++j){
      if (grid->GetDim() == 3){
        rNN = sqrt(pow(neighbors[i][0]-neighbors[j][0],2.0) + pow(neighbors[i][1]-neighbors[j][1],2.0) + pow(neighbors[i][2]-neighbors[j][2],2.0));
      }else{
        rNN = sqrt(pow(neighbors[i][0]-neighbors[j][0],2.0) + pow(neighbors[i][1]-neighbors[j][1],2.0));
      }
      ALoc[i][j] = pow(1.0 - rNN/alpha, 2.0);
    }
    switch( dim ){
    case 0:
      EXCEPTION("Curl of a scalar field!");
      break;
    case 1:
      vals.Resize(numNN,2);
      vals[i][0] = vectors[i][0];
      vals[i][1] = vectors[i][1];
      if (grid->GetDim() == 2){
        derivVec.Resize(numNN,2);
        if (l2Distances[i] == 0) {
          derivVec[i][0] = 0.0;
          derivVec[i][1] = 0.0;
        }else{
          derivVec[i][0] = 2.0 * (1.0 - l2Distances[i]/alpha) * (neighbors[i][0] - globPoint[0])/(l2Distances[i]*alpha);
          derivVec[i][1] = 2.0 * (1.0 - l2Distances[i]/alpha) * (neighbors[i][1] - globPoint[1])/(l2Distances[i]*alpha);
        }
      }else{
        EXCEPTION("2D mesh and 3D-values!")
      }
      break;
    case 2:
      vals.Resize(numNN,3);
      vals[i][0] = vectors[i][0];
      vals[i][1] = vectors[i][1];
      vals[i][2] = vectors[i][2];
      if (grid->GetDim() == 3){
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
        EXCEPTION("3D mesh and 2D-values!")
      }      break;
    }
  }

  // now we have to invert Aloc and multiply it with the according value-coloumn
  ALoc.Invert_Lapack();

  CF::Matrix<Double> temp;

  // coefficient matrix (coloumn nr. corresponding to the spatial dimension)
  derivCoefVec = ALoc * derivVec;
  vals.Transpose(temp);

  CF::Matrix<Double> tempvec;
  tempvec.Resize(1,dim+1);
  tempvec = temp * derivCoefVec;

  //now we have to combine the tempvec-entries in order to obtain the curl
  if(grid->GetDim() == 2){
    vec[0][0] = tempvec[1][0] - tempvec[0][1];
  }else{
    vec[0][0] = tempvec[2][1] - tempvec[1][2];
    vec[0][1] = tempvec[0][2] - tempvec[2][0];
    vec[0][2] = tempvec[1][0] - tempvec[0][1];
  }
}



void MeshFilter::CalcLocGradient(CF::Matrix<Double>& vec,
                                  const CF::Vector<Double>& globPoint,
                                  const CF::StdVector< Vector<Double> >& neighbors,
                                  const CF::StdVector< Double >& l2Distances,
                                  const CF::StdVector< Vector<Double> >& vectors,
                                  const UInt numNN,
                                  const Double alpha,
                                  const UInt dim,
                                  Grid* grid){

  Matrix<Double> derivCoefVec;
  derivCoefVec.Resize(numNN,1);
  CF::Matrix<Double> ALoc;
  ALoc.Resize(numNN,numNN);
  CF::Matrix<Double> vals;
  CF::Matrix<Double> derivVec; //Vector of RBF derivatives evaluated at srcPoints

  Double rNN; //distance between two src points
  for (UInt i = 0; i < numNN; ++i){
    for (UInt j = 0; j < numNN; ++j){
      if (grid->GetDim() == 3){
        rNN = sqrt(pow(neighbors[i][0]-neighbors[j][0],2.0) + pow(neighbors[i][1]-neighbors[j][1],2.0) + pow(neighbors[i][2]-neighbors[j][2],2.0));
      }else{
        rNN = sqrt(pow(neighbors[i][0]-neighbors[j][0],2.0) + pow(neighbors[i][1]-neighbors[j][1],2.0));
      }
      ALoc[i][j] = pow(1.0 - rNN/alpha, 2.0);
    }
    switch( dim ){
    case 0:
      vals.Resize(numNN,1);
      vals[i][0] = vectors[i][0];
      if (grid->GetDim() == 2){
      derivVec.Resize(numNN,2);
      if (l2Distances[i] == 0) {
        derivVec[i][0] = 0.0;
        derivVec[i][1] = 0.0;
      }else{
        derivVec[i][0] = 2.0 * (1.0 - l2Distances[i]/alpha) * (neighbors[i][0] - globPoint[0])/(l2Distances[i]*alpha);
        derivVec[i][1] = 2.0 * (1.0 - l2Distances[i]/alpha) * (neighbors[i][1] - globPoint[1])/(l2Distances[i]*alpha);
      }
      }else{
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
      }
      break;
    case 1:
      EXCEPTION("Gradient of Vector not defined in this context!");
      break;
    case 2:
      EXCEPTION("Gradient of Vector not defined in this context!");
      break;
    }
 }



  // now we have to invert Aloc and multiply it with the according value-coloumn
  ALoc.Invert_Lapack();

  CF::Matrix<Double> temp;
  // coefficient matrix (coloumn nr. corresponding to the spatial dimension)
  derivCoefVec = ALoc * derivVec;
  vals.Transpose(temp);

  vec = temp * derivCoefVec;
}

void MeshFilter::CalcLocDivergence(CF::Matrix<Double>& vec,
                                  const CF::Vector<Double>& globPoint,
                                  const CF::StdVector< Vector<Double> >& neighbors,
                                  const CF::StdVector< Double >& l2Distances,
                                  const CF::StdVector< Vector<Double> >& vectors,
                                  const UInt numNN,
                                  const Double alpha,
                                  const UInt dim,
                                  Grid* grid){
  CF::Matrix<Double> derivCoefVec;
    CF::Matrix<Double> ALoc;
    ALoc.Resize(numNN,numNN);
    CF::Matrix<Double> vals;
    CF::Matrix<Double> derivVec; //Vector of RBF derivatives evaluated at srcPoints

    Double rNN; //distance between two src points
    for (UInt i = 0; i < numNN; ++i){
      for (UInt j = 0; j < numNN; ++j){
        if (grid->GetDim() == 3){
          rNN = sqrt(pow(neighbors[i][0]-neighbors[j][0],2.0) + pow(neighbors[i][1]-neighbors[j][1],2.0) + pow(neighbors[i][2]-neighbors[j][2],2.0));
        }else{
          rNN = sqrt(pow(neighbors[i][0]-neighbors[j][0],2.0) + pow(neighbors[i][1]-neighbors[j][1],2.0));
        }
        ALoc[i][j] = pow(1.0 - rNN/alpha, 2.0);
      }
      switch( dim ){
      case 0:
        //should already be caught
        EXCEPTION("Divergence of a scalar field!");
        break;
      case 1:
        vals.Resize(numNN,2);
        vals[i][0] = vectors[i][0];
        vals[i][1] = vectors[i][1];
        if (grid->GetDim() == 2){
        derivVec.Resize(numNN,2);
        if (l2Distances[i] == 0) {
          derivVec[i][0] = 0.0;
          derivVec[i][1] = 0.0;
        }else{
          derivVec[i][0] = 2.0 * (1.0 - l2Distances[i]/alpha) * (neighbors[i][0] - globPoint[0])/(l2Distances[i]*alpha);
          derivVec[i][1] = 2.0 * (1.0 - l2Distances[i]/alpha) * (neighbors[i][1] - globPoint[1])/(l2Distances[i]*alpha);
        }
        }else{
          EXCEPTION("2D values and 3D mesh!");
        }
        break;
      case 2:
        vals.Resize(numNN,3);
        vals[i][0] = vectors[i][0];
        vals[i][1] = vectors[i][1];
        vals[i][2] = vectors[i][2];
        if (grid->GetDim() == 3){
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
           WARN("DivergenceDifferentiator.cc : Treat 3D values as 2D values, due to a 2D mesh!")
           derivVec.Resize(numNN,2);
           if (l2Distances[i] == 0) {
             derivVec[i][0] = 0.0;
             derivVec[i][1] = 0.0;
           }else{
             derivVec[i][0] = 2.0 * (1.0 - l2Distances[i]/alpha) * (neighbors[i][0] - globPoint[0])/(l2Distances[i]*alpha);
            derivVec[i][1] = 2.0 * (1.0 - l2Distances[i]/alpha) * (neighbors[i][1] - globPoint[1])/(l2Distances[i]*alpha);
           }
          //EXCEPTION("3D-values and 2D-mesh!");

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
    tempvec.Resize(1,dim+1);
    tempvec = temp * derivCoefVec;
    //in order to obtain the divergence of the vectorfield, we have to add tempvec[:,0]+tempvec[:,1]+tempvec[:,2]
    vec[0][0] = tempvec.Trace();
}

}


