// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     AeroacousticBase.cc
 *       \brief    <Description>
 *
 *       \date     Jan 7, 2017
 *       \author   kroppert
 */
//================================================================================================

#include "AeroacousticBase.hh"
#include "cfsdat/Utils/KNNSearch.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include "FeBasis/H1/H1Elems.hh"


namespace CFSDat{

AeroacousticBase::AeroacousticBase(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
:MeshFilter(numWorkers,config,resMan){


}

AeroacousticBase::~AeroacousticBase(){

}



void AeroacousticBase::OmegaVectorProductU_2D(const PhysicalEntity& data,
                                              const std::vector<QuantityStruct>& derivData,
                                              Vector<Double>& retVec){
  Vector<Double> vec;
  //if we have an extern velocity, we have node-results, otherwise element-results
  (data.externVorticity)? vec = data.U : vec = data.Uelem;

  // Check: in 2D, vorticity is scalar (pyhysically it points in z direction) and velocity 2-dimensional
  if( 2 * data.Om.GetSize() != vec.GetSize()){
    std::cout<<"size of vorticity-vector:"<<data.Om.GetSize()<<std::endl;
    std::cout<<"size of velocity-vector:"<<vec.GetSize()<<std::endl;
    EXCEPTION("AeroacousticBase.cc: Vorticity and velocity are inconsistent");
  }
  if( data.Om.GetSize() != derivData.size()){
    std::cout<<"size of derivative-data-vector:"<<derivData.size()<<std::endl;
    std::cout<<"size of vorticity-vector:"<<data.Om.GetSize()<<std::endl;
    EXCEPTION("AeroacousticBase.cc: Vorticity and velocity are inconsistent");
  }

  retVec.Resize(data.dim * derivData.size());
  retVec.Init();


  Vector<Double> om, u, t;
  om.Resize(3,0.0);
  u.Resize(3,0.0);
  t.Resize(3,0.0);

  //2D case
  for(UInt i = 0; i < derivData.size(); ++i){
    // extract vectors from Omega and U to perform the vector product
     om[2] = data.Om[i];
     u[0] = vec[2 * i];
     u[1] = vec[2 * i + 1];
     //compute the vector product
     om.CrossProduct(u,t);
     //2D case:
     retVec[2 * i] = t[0];
     retVec[2 * i + 1] = t[1];
  }
}


void AeroacousticBase::OmegaVectorProductU_3D(const PhysicalEntity& data,
                                              const std::vector<QuantityStruct>& derivData,
                                              Vector<Double>& retVec){


  // Check if velocity and vorticity have the same size
  if(data.Om.GetSize() != data.U.GetSize()){
    std::cout<<"size of vorticity-vector:"<<data.Om.GetSize()<<std::endl;
    std::cout<<"size of velocity-vector:"<<data.U.GetSize()<<std::endl;
    EXCEPTION("AeroacousticBase.cc: Vorticity and velocity are inconsistent")
  }

  retVec.Resize(3 * derivData.size());
  retVec.Init();

  Vector<Double> om, u, t;
  om.Resize(3);
  u.Resize(3);
  t.Resize(3);


  //3D case
  for(UInt i = 0; i < derivData.size(); ++i){
    // extract vectors from Omega and U to perform the vector product
    for(UInt d = 0; d < 3; ++d){ // vectors must have dimension 3, for the CrossProduct()
        om[d] = data.Om[(data.dim) * i + d];
        u[d] = data.U[(data.dim) * i + d];
      }
    //insert the vector product in retVec
    om.CrossProduct(u,t);
      for(UInt d = 0; d < data.dim; ++d){
        retVec[(data.dim) * i + d] = t[d];
      }
    }
}




void AeroacousticBase::CalcDivergence2D(Vector<Double>& returnVec,
                                    const Vector<Double>& inVec,
                                    const CF::StdVector< CF::Vector<double> >& coords,
                                    const std::vector<QuantityStruct>& derivData){

  Vector<Double> tempRetVec;
  tempRetVec.Resize(2 * derivData.size());

  // vector which will be filled with the source values (scattered data values), in order to
  // perform a nearest neighbour search
  CF::StdVector< CF::Vector<Double> >  scatteredData;
  scatteredData.Resize(coords.GetSize());

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

  str1::shared_ptr<EqnMapSimple> downMap = resultManager_->GetResultAdapter(filterResIds[0])->mapping;

  UInt dim;
  FillScatteredDataVec(scatteredData, vec, inVec, coords, dim);


  //now we bring the scattered coordinates and data values into
  //the correct form for the nearest neighbour-search witch CGAL

  // Object for nearest neighbor-searches and bringing the data into the correct form for CGAL search
  KNNSearch Tree;
  Tree.ReadScatteredData_Div(coords, dim, trgGrid_, scatteredData);


  // loop over all elements
for(UInt i = 0; i < derivData.size();++i){

    QuantityStruct aStru = derivData[i];
    const Elem* curE = trgGrid_->GetElem(aStru.trgElemNum);

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


          CalcLocDivergence(vec, globalCoord, neighbors, l2dists, vectors, numNN, alpha, dim, trgGrid_);

          // this eqns is now two-dimensional, that is the reason why we have to remap the
          // tempRetVec to a scalar-valued vector to fulfill the scalar-valued divergence
          CF::StdVector<UInt> eqns;
            //get the equation map for the nodes in eConn, in order to insert the
            //interpolation in the correct position in the result vector
            downMap->GetEquation(eqns,curE->elemNum,ExtendedResultInfo::ELEMENT);

            //if scalar input values of scattered data->inDim_=0
            //if it is a two-dimensional vector->inDim_=1 and inDim_=2 for 3d-vector
            for(UInt aDof = 0; aDof < eqns.GetSize(); aDof++){
              tempRetVec[eqns[aDof]] = vec[0][0];
            }

  }// this tempRetVec is now two-dimensional, this is the reason why we have to remap the
    // tempRetVec to a scalar-valued vector to fulfill the scalar-valued divergence
  returnVec = TwoDToScalar(tempRetVec);

}


void AeroacousticBase::CalcCurlU_2D(Vector<Double>& retVec,
                                 const CF::StdVector< Vector<Double> >&  scatteredData,
                                 const std::vector<QuantityStruct>& derivData,
                                 const CF::StdVector< CF::Vector<double> >& coords,
                                 const boost::uuids::uuid& resId,
                                 const UInt& dim){

  //resize retVec to number of nodes of target
  retVec.Resize(derivData.size());
  retVec.Init();

  Vector<Double> tempRetVec;
  tempRetVec.Resize(2*derivData.size());

  // actually this is the curl vector of the curl-operator, but Matrix class has to be used, for matrix-vector product
  CF::Matrix<Double> vec;
  vec.Resize(1);


  // Object for nearest neighbor-searches and bringing the data into the correct form for CGAL search
  KNNSearch Tree;
  Tree.ReadScatteredData_Lighthill(coords, dim - 1, trgGrid_, scatteredData);

  // nr. of nearest neighbours for the nodal RBF basis functions
  UInt numNN = (trgGrid_->GetDim() == 2) ? 4 : 8 ;

  // coordinate list of nearest neighbour points
  CF::StdVector< Vector<Double> > neighbors;

  // distances according to every nearest neighbour point
  CF::StdVector< Double > l2dists;

  // vector containing the values of each nearest neighbour point
  CF::StdVector< Vector<Double> > vectors;

  str1::shared_ptr<EqnMapSimple> downMap = resultManager_->GetResultAdapter(resId)->mapping;


  // loop over all elements and over every node of each element
  for(UInt i = 0; i < derivData.size();++i){
    QuantityStruct aStru = derivData[i];
    const Elem* curE = trgGrid_->GetElem(aStru.trgElemNum);

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


    CalcLocCurl(vec, globalCoord, neighbors, l2dists, vectors, numNN, alpha, dim - 1, trgGrid_);

    CF::StdVector<UInt> eqns;
    //get the equation map for the nodes in eConn, in order to insert the
    //interpolation in the correct position in the result vector
    downMap->GetEquation(eqns,curE->elemNum,ExtendedResultInfo::ELEMENT);

    tempRetVec[eqns[0]] = vec[0][0];

  }//for all elements

  // tempRetVec is 2D but in our convention it's 1D, so we remap it
  retVec = TwoDToScalar(tempRetVec);

}


void AeroacousticBase::CalcCurlU_3D(Vector<Double>& retVec,
                                 const CF::StdVector< Vector<Double> >&  scatteredData,
                                 const std::vector<QuantityStruct>& derivData,
                                 const CF::StdVector< CF::Vector<double> >& coords,
                                 const boost::uuids::uuid& resId,
                                 const UInt& dim){

  //resize retVec to number of nodes of target
  (dim == 2) ? retVec.Resize(2 * derivData.size()) : retVec.Resize(3 * derivData.size());
  retVec.Init();


  // actually this is the curl vector of the curl-operator, but Matrix class has to be used, for matrix-vector product
  CF::Matrix<Double> vec;
  vec.Resize(dim);


  // Object for nearest neighbor-searches and bringing the data into the correct form for CGAL search
  KNNSearch Tree;
  Tree.ReadScatteredData_Lighthill(coords, dim - 1, trgGrid_, scatteredData);

  // nr. of nearest neighbours for the nodal RBF basis functions
  UInt numNN = (trgGrid_->GetDim() == 2) ? 4 : 8 ;

  // coordinate list of nearest neighbour points
  CF::StdVector< Vector<Double> > neighbors;

  // distances according to every nearest neighbour point
  CF::StdVector< Double > l2dists;

  // vector containing the values of each nearest neighbour point
  CF::StdVector< Vector<Double> > vectors;

  str1::shared_ptr<EqnMapSimple> downMap = resultManager_->GetResultAdapter(resId)->mapping;


  // loop over all elements and over every node of each element
  for(UInt i = 0; i < derivData.size();++i){
    QuantityStruct aStru = derivData[i];
    const Elem* curE = trgGrid_->GetElem(aStru.trgElemNum);

    //get the global coordinates of element centroid
    CF::Vector<Double> globalCoord;
    trgGrid_->GetElemCentroid(globalCoord, curE->elemNum , false);

    for (UInt dof = 0; dof < dim; dof++){
      vec[0][dof] = 0.0;
    }

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


    CalcLocCurl(vec, globalCoord, neighbors, l2dists, vectors, numNN, alpha, dim - 1, trgGrid_);

    CF::StdVector<UInt> eqns;
    //get the equation map for the nodes in eConn, in order to insert the
    //interpolation in the correct position in the result vector
    downMap->GetEquation(eqns,curE->elemNum,ExtendedResultInfo::ELEMENT);


    for(UInt aDof = 0; aDof < eqns.GetSize(); aDof++){
      if(dim == 2){
        retVec[eqns[0]] = vec[0][0];
      }else{
      retVec[eqns[aDof]] = vec[0][aDof];
      }
    }
  }//for all elements


}


void AeroacousticBase::CalcGradUScalarU(Vector<Double>& retVec,
                                      const PhysicalEntity& data,
                                      const std::vector<QuantityStruct>& derivDataElem,
                                      const CF::StdVector< CF::Vector<double> >& sourceCoords,
                                      const boost::uuids::uuid& resId){

  // Check if velocity vector and source-coordinates have the same size
  if( data.U.GetSize() != data.dim * sourceCoords.GetSize()){
    std::cout<<"size of velocity-vector:"<<data.U.GetSize()<<std::endl;
    std::cout<<"size of source-coordinates-vector:"<<sourceCoords.GetSize()<<std::endl;
    EXCEPTION("AeroacousticBase.cc: Velocity vector and source coordinates are inconsistent")
  }


  UInt inDim = data.dim;

  Vector<Double> inVecVelocity = data.U;

  //resize retVec to number of element-centroids of target
  retVec.Resize(inDim * derivDataElem.size());
  retVec.Init();

  // actually this is the gradient vector, but Matrix class has to be used, for matrix-vector product
  CF::Matrix<Double> vec;
  vec.Resize(inDim);
  vec.Init();

  // Now we want to get the gradient of u*u
  // but first of all we need the scalar value u*u
  // therefore we don't use data.scatteredData because we compute the
  // scalar product right here
  CF::StdVector< CF::Vector<Double> >  scatteredData;
  scatteredData.Resize(sourceCoords.GetSize());
  if(inVecVelocity.GetSize() == sourceCoords.GetSize()){
    //this is the case of scalar scattered values
    EXCEPTION("AeroacousticBase.cc: Velocity is a scalar value -> can't be correct!")
  }else{
    if(inVecVelocity.GetSize() == 2 * sourceCoords.GetSize()){
      //case of two dimensional vector
      for(UInt i = 0; i < sourceCoords.GetSize() ; ++i){
        scatteredData[i].Resize(1);
        scatteredData[i][0] = inVecVelocity[2 * i]*inVecVelocity[2 * i] + inVecVelocity[2 * i + 1]*inVecVelocity[2 * i + 1];
      }
    }else{
      if(inVecVelocity.GetSize() == 3 * sourceCoords.GetSize()){
        // case of three dimensional vector
        for(UInt i = 0; i < sourceCoords.GetSize(); ++i){
          scatteredData[i].Resize(1);
          scatteredData[i][0] = inVecVelocity[3 * i]*inVecVelocity[3 * i] + inVecVelocity[3 * i + 1]*inVecVelocity[3 * i + 1] + inVecVelocity[3 * i + 2]*inVecVelocity[3 * i + 2];
        }
      }else{
        EXCEPTION("AeroacousticBase.cc: Inconsistency between input data and source coordinates!")
      }
    }
  }

  // now our values at the nodes are scalar quantities which means dimension = 0
  UInt dim = 0;

  // Object for nearest neighbor-searches and bringing the data into the correct form for CGAL search
  KNNSearch Tree;
  Tree.ReadScatteredData_Lighthill(sourceCoords, dim, trgGrid_, scatteredData);

  // nr. of nearest neighbours for the nodal RBF basis functions
  UInt numNN = (trgGrid_->GetDim() == 2) ? 4 : 8 ;

  // coordinate list of nearest neighbour points
  CF::StdVector< Vector<Double> > neighbors;

  // distances according to every nearest neighbour point
  CF::StdVector< Double > l2dists;

  // vector containing the values of each nearest neighbour point
  CF::StdVector< Vector<Double> > vectors;

  str1::shared_ptr<EqnMapSimple> downMap = resultManager_->GetResultAdapter(resId)->mapping;

  // loop over all elements
  for(UInt i = 0; i < derivDataElem.size();++i){

    QuantityStruct aStru = derivDataElem[i];
    const Elem* curE = trgGrid_->GetElem(aStru.trgElemNum);

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

    // build the local RBF interpolation matrix, based on the nearest neighbour source points
    // build the local interpolation-value vector and solve the system ALoc*c=vector for
    // the local RBF coefficients c
    CalcLocGradient(vec, globalCoord, neighbors, l2dists, vectors, numNN, alpha, dim, trgGrid_);

    CF::StdVector<UInt> eqns;
    //get the equation map for the nodes in eConn, in order to insert the
    //interpolation in the correct position in the result vector
    downMap->GetEquation(eqns,curE->elemNum,ExtendedResultInfo::ELEMENT);

    // dimension of the equations should be correct, 2 for 2D and 3 for 3D
    for(UInt aDof = 0; aDof < eqns.GetSize(); aDof++){
      retVec[eqns[aDof]] = vec[0][aDof];
    }
  }//for all elements
}

Vector<Double> AeroacousticBase::ScalarToTwoD(const Vector<Double>& inVec){
  Vector<Double> r;
  r.Resize(2 * inVec.GetSize());
  for(UInt i = 0; i < inVec.GetSize(); ++i){
    r[2 * i] = inVec[i];
    r[2 * i + 1] = 0.0;
  }
  return r;
}

Vector<Double> AeroacousticBase::TwoDToScalar(const Vector<Double>& inVec){
  Vector<Double> r;
  r.Resize(inVec.GetSize() / 2);
  for(UInt i = 0; i < r.GetSize(); ++i){
    r[i] = inVec[2 * i];
  }
  return r;
}


}

