// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     Lighthill.cc
 *       \brief    <Description>
 *
 *       \date     Oct 6, 2016
 *       \author   kroppert
 */
//================================================================================================


#include "Lighthill.hh"
#include "FeBasis/H1/H1Elems.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include "cfsdat/Utils/KNNSearch.hh"
#include <algorithm>
#include <vector>
#include <math.h>


namespace CFSDat{

Lighthill::Lighthill(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan)
:MeshFilter(numWorkers,config,resMan){

  this->filtStreamType_ = FIFO_FILTER;
  inDim_ = 0;

  //check if scheme has been defined in the xml-file
  if(config->Has("scheme") == false){
    EXCEPTION("============================================================ \n"
        "No <sheme> for Lighthill computation was set in xml-file! \n"
        "so default values are chosen \n"
        "============================================================");
  }


    std::cout<<("============================================================ \n"
        "Make sure you set the ''density='' in the xml-scheme correctly, \n"
        "otherwise it is chosen to be 1.0 !! \n"
        "============================================================")<<std::endl;


    density_ = config->Get("scheme")->Get("density")->As<Double>();

    Form_ = config->Get("type")->As<std::string>();

}

Lighthill::~Lighthill(){

}

bool Lighthill::Run(){
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


  /**********************************************************
   * Lighthill-source terms to be computed:
   * Div(Grad(1/2 u*u)) + rho * Div(Curl(u) x u)
   * in our case we collect the terms the following way:
   * Div{Grad(1/2 u*u) + rho * (Curl(u) x u)}
   * in this filter we compute everything but the divergence
   * therefore the elemResult output must be interpolated to nodes
   * and then we can apply the divergence filter to obtain the final
   * source-terms
   ***********************************************************/

  /**********************************************************
   * Check if the whole Lighthill source or only the Lamb-vector
   * has to be computed
   ***********************************************************/
  Vector<Double> GradOfScalar;
  if(Form_ != "AeroacousticSource_LambVector"){
  //Grad(1/2 u*u)
  CalcGradUScalarU(GradOfScalar, inVec);
  GradOfScalar.ScalarMult(0.5);
  }

  //Curl(u)
  Vector<Double> CurlU;
  CalcCurlU(CurlU, inVec);

  //u in inVec is defined on nodes but when multiplying it
  //to differentiated results (which are elemResults) we
  //need to restrict them to elem centroids
  Vector<Double> interpolateU;
  InterpolationNodeToCenter(interpolateU, inVec);

  //Curl(u) x u
  Vector<Double> OmegaCrossU;
  OmegaVectorProductU(CurlU, interpolateU, OmegaCrossU);


  if(Form_ != "AeroacousticSource_LambVector"){
  //adding all up
  returnVec = GradOfScalar+OmegaCrossU*density_;
  }else{
    returnVec = OmegaCrossU * density_;
  }


  // Check filter mesh and output values
  if(true){
	  // Check output values
	  // d'Alembertoperator(p)= d(Source_i)/dx_i
	  // (a) INTEGRATION(Source_i) over the domain must be approximately zero

	  // (b) INTEGRATION(d(Source_i)/dx_i) over the domain must be approximately zero
	  //                                   if Source_i(Wall) = 0 and Source_i(In-/Outlet)=0

	  // Check mesh
	  // Mesh size must be smaller than << Ma*c/(20*f) in the source region
	  // to determine acoustically relevant flow structures
  }

  resultManager_->ActivateResult(filterResIds[0]);

  //now deactivate own upstream results
  for(UInt aRes=0;aRes<upResIds.GetSize();aRes++){
    resultManager_->DeactivateResult(upResIds[aRes]);
  }

  return true;
}



void Lighthill::OmegaVectorProductU(const Vector<Double> Omega,
                            const Vector<Double> U,
                            Vector<Double>& retVec){

  CF::StdVector< CF::Vector<Double> >  OmegaData;
  OmegaData.Resize(derivData_.size());

  CF::StdVector< CF::Vector<Double> >  UData;
  UData.Resize(derivData_.size());

  CF::StdVector< CF::Vector<Double> >  OmegaCrossU;
  OmegaCrossU.Resize(derivData_.size());

  retVec.Resize(3*derivData_.size());
  retVec.Init();

  // case of three dimensional vector
  for(UInt i = 0; i < derivData_.size(); ++i){
    OmegaData[i].Resize(3);
    OmegaData[i][0] = Omega[3 * i]; // x-component
    OmegaData[i][1] = Omega[3 * i + 1]; // y-component
    OmegaData[i][2] = Omega[3 * i + 2]; // z-component
    UData[i].Resize(3);
    UData[i][0] = U[3 * i]; // x-component
    UData[i][1] = U[3 * i + 1]; // y-component
    UData[i][2] = U[3 * i + 2]; // z-component
  }

  //now we compute the cross product
  for(UInt i = 0; i < derivData_.size(); ++i){
    OmegaCrossU[i].Resize(3);
    OmegaCrossU[i][0] = OmegaData[i][1] * UData[i][2] - OmegaData[i][2] * UData[i][1];
    OmegaCrossU[i][1] = OmegaData[i][2] * UData[i][0] - OmegaData[i][0] * UData[i][2];
    OmegaCrossU[i][2] = OmegaData[i][0] * UData[i][1] - OmegaData[i][1] * UData[i][0];
  }

  //bring in such a returnVec-form (x1,y1,z1,x2,y2,z2,...)
  for (UInt i = 0; i < derivData_.size(); ++i){
    retVec[3 * i] = OmegaCrossU[i][0];
    retVec[3 * i + 1] = OmegaCrossU[i][1];
    retVec[3 * i + 2] = OmegaCrossU[i][2];
  }

}



void Lighthill::CalcGradUScalarU(Vector<Double>& retVec,
                                  const Vector<Double> inVec){


  //resize retVec to number of element-centroids of target
  retVec.Resize(3.0 * derivData_.size());
  retVec.Init();

  // actually this is the gradient vector, but Matrix class has to be used, for matrix-vector product
  CF::Matrix<Double> vec;
  vec.Resize(3);
  vec.Init();

  // Now we want to get the gradient of u*u
  // but first of all we need the scalar value u*u

  // vector which we fill with the source values (scattered data values), in order to
  // perform a nearest neighbour search
  CF::StdVector< CF::Vector<Double> >  scatteredData;
  scatteredData.Resize(sourceCoords_.GetSize());

  if(inVec.GetSize() == sourceCoords_.GetSize()){
    //this is the case of scalar scattered values
    EXCEPTION("not a vector value")
  }else{
    if(inVec.GetSize() == 2 * sourceCoords_.GetSize()){
      //case of two dimensional vector
      inDim_=0;
      for(UInt i = 0; i < sourceCoords_.GetSize() ; ++i){
        scatteredData[i].Resize(1);
        scatteredData[i][0] = inVec[2 * i]*inVec[2 * i] + inVec[2 * i + 1]*inVec[2 * i + 1];
      }
    }else{
      if(inVec.GetSize() == 3 * sourceCoords_.GetSize()){
        // case of three dimensional vector
        inDim_=0;
        for(UInt i = 0; i < sourceCoords_.GetSize(); ++i){
          scatteredData[i].Resize(1);
          scatteredData[i][0] = inVec[3 * i]*inVec[3 * i] + inVec[3 * i + 1]*inVec[3 * i + 1] + inVec[3 * i + 2]*inVec[3 * i + 2];
        }
      }else{
        EXCEPTION("Inconsistency between input data and source coordinates!")
      }
    }
  }

  // Object for nearest neighbor-searches and bringing the data into the correct form for CGAL search
  KNNSearch Tree;
  Tree.ReadScatteredData_Lighthill(sourceCoords_, inDim_, trgGrid_, scatteredData);

  // nr. of nearest neighbours for the nodal RBF basis functions
  UInt numNN = (trgGrid_->GetDim() == 2) ? 4 : 8 ;

  // coordinate list of nearest neighbour points
  CF::StdVector< Vector<Double> > neighbors;

  // distances according to every nearest neighbour point
  CF::StdVector< Double > l2dists;

  // vector containing the values of each nearest neighbour point
  CF::StdVector< Vector<Double> > vectors;

  str1::shared_ptr<EqnMapSimple> downMap = resultManager_->GetResultAdapter(filterResIds[0])->mapping;

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

    // build the local RBF interpolation matrix, based on the nearest neighbour source points
    // build the local interpolation-value vector and solve the system ALoc*c=vector for
    // the local RBF coefficients c
    CalcLocGradDerivativeCoefs(vec, globalCoord, neighbors, l2dists, vectors, numNN, alpha);

    CF::StdVector<UInt> eqns;
    //get the equation map for the nodes in eConn, in order to insert the
    //interpolation in the correct position in the result vector
    downMap->GetEquation(eqns,curE->elemNum,ExtendedResultInfo::ELEMENT);

    //if scalar input values of scattered data->inDim_=0
    //if it is a two-dimensional vector->inDim_=1 and inDim_=2 for 3d-vector
    for(UInt aDof = 0; aDof < eqns.GetSize(); aDof++){
      retVec[eqns[aDof]] = vec[0][aDof];
    }
  }//for all elements
}



void Lighthill::CalcCurlU(Vector<Double>& retVec,
                         const Vector<Double> inVec){

  //resize retVec to number of nodes of target
  retVec.Resize(3.0 * derivData_.size());
  retVec.Init();

  // vector which will be filled with the source values (scattered data values), in order to
  // perform a nearest neighbour search
  CF::StdVector< CF::Vector<Double> >  scatteredData;
  scatteredData.Resize(sourceCoords_.GetSize());

  // actually this is the curl vector of the curl-operator, but Matrix class has to be used, for matrix-vector product
  CF::Matrix<Double> vec;
  vec.Resize(3);


  // Checking if input is scalar- or vector-type. This needs to be done
  // because of the input of CGAL and FLANN
  // Maybe there is a more efficient way to deal with this issue?!
  if(inVec.GetSize() == sourceCoords_.GetSize()){
    //this is the case of scalar scattered values
    EXCEPTION("Divergence of a scalar field!");
  }else{
    if(inVec.GetSize() == 2 * sourceCoords_.GetSize()){
      //case of two dimensional vector
      inDim_=1;
      for(UInt i = 0; i < sourceCoords_.GetSize() ; ++i){
        scatteredData[i].Resize(2);
        scatteredData[i][0] = inVec[2 * i]; // x-component
        scatteredData[i][1] = inVec[2 * i + 1]; // y-component
        scatteredData[i][1] = 0.0;   // empty z-component
      }
    }else{
      if(inVec.GetSize() == 3 * sourceCoords_.GetSize()){
        // case of three dimensional vector
        inDim_=2;
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
  Tree.ReadScatteredData_Lighthill(sourceCoords_, inDim_, trgGrid_, scatteredData);

  // nr. of nearest neighbours for the nodal RBF basis functions
  UInt numNN = (trgGrid_->GetDim() == 2) ? 4 : 8 ;

  // coordinate list of nearest neighbour points
  CF::StdVector< Vector<Double> > neighbors;

  // distances according to every nearest neighbour point
  CF::StdVector< Double > l2dists;

  // vector containing the values of each nearest neighbour point
  CF::StdVector< Vector<Double> > vectors;

  str1::shared_ptr<EqnMapSimple> downMap = resultManager_->GetResultAdapter(filterResIds[0])->mapping;


  // loop over all elements and over every node of each element
  for(UInt i = 0; i < derivData_.size();++i){
    DifferentiationStruct& aStru = derivData_[i];
    const Elem* curE = trgGrid_->GetElem(aStru.tENum);

    //get the global coordinates of element centroid
    CF::Vector<Double> globalCoord;
    trgGrid_->GetElemCentroid(globalCoord, curE->elemNum , false);

    for (UInt dof = 0; dof < inDim_ + 1; dof++){
      vec[dof][0] = 0.0;
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

    // build the local RBF interpolation matrix, based on the nearest neighbour source points
    // build the local interpolation-value vector and solve the system ALoc*c=vector for
    // the local RBF coefficients c
    CalcLocCurlDerivativeCoefs(vec, globalCoord, neighbors, l2dists, vectors, numNN, alpha);

    CF::StdVector<UInt> eqns;
    //get the equation map for the nodes in eConn, in order to insert the
    //interpolation in the correct position in the result vector
    downMap->GetEquation(eqns,curE->elemNum,ExtendedResultInfo::ELEMENT);

    //if scalar input values of scattered data->inDim_=0
    //if it is a two-dimensional vector->inDim_=1 and inDim_=2 for 3d-vector
    for(UInt aDof = 0; aDof < eqns.GetSize(); aDof++){
      retVec[eqns[aDof]] = vec[0][aDof];
    }
  }//for all elements
}


void Lighthill::InterpolationNodeToCenter(Vector<Double>& retVec, const Vector<Double> inVec){

  //resize retVec to number of nodes of target
  retVec.Resize(3.0 * derivData_.size());
  retVec.Init();

  CF::StdVector<UInt> eqnNums;

  CF::StdVector< CF::Vector<Double> >  scatteredData;
  scatteredData.Resize(sourceCoords_.GetSize());

  UInt numNN = 15;
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
    for(UInt i = 0; i < sourceCoords_.GetSize(); ++i){
      scatteredData[i].Resize(1);
      scatteredData[i][0] = inVec[i];
    }
  }else{
    if(inVec.GetSize() == 2 * sourceCoords_.GetSize()){
      //case of two dimensional vector
      inDim_=1;
      vec.Resize(2);
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
  Tree.ReadScatteredData_Lighthill(sourceCoords_, inDim_, trgGrid_, scatteredData);

  // coordinate list of nearest neighbour points
  CF::StdVector< Vector<Double> > neighbors;
  // distances according to every nearest neighbour point
  CF::StdVector< Double > l2dists;
  // vector containing the values of each nearest neighbour point
  CF::StdVector< Vector<Double> > vectors;

  // loop over all elements and over every node of each element
  for(UInt i = 0; i < derivData_.size();++i){

    DifferentiationStruct& aStru = derivData_[i];
    const Elem* curE = trgGrid_->GetElem(aStru.tENum);

    //get the global coordinates of element centroid
    CF::Vector<Double> globalCoord;
    trgGrid_->GetElemCentroid(globalCoord, curE->elemNum , false);

    for (UInt dof = 0; dof < inDim_ + 1; dof++){
      vec[dof] = 0.0;
    }

    // at that point we can start obtaining the nearest neighbours for every aNode
    Tree.KNN_CGAL_Differentiation(globalCoord, neighbors, l2dists, vectors, numNN);


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
    //interpolation exponent
    Double p=4.0;
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
    downMap->GetEquation(eqns,curE->elemNum,ExtendedResultInfo::ELEMENT);

    //if scalar input values of scattered data->inDim_=0
    //if it is a two-dimensional vector->inDim_=1 and inDim_=2 for 3d-vector
    for(UInt aDof = 0; aDof < inDim_+1; aDof++){
      retVec[eqns[aDof]] = vec[aDof];
    }
  }// for element
}


void Lighthill::CalcLocGradDerivativeCoefs(CF::Matrix<Double>& vec,
                                            const CF::Vector<Double>& globPoint,
                                            const CF::StdVector< Vector<Double> >& neighbors,
                                            const CF::StdVector< Double >& l2Distances,
                                            const CF::StdVector< Vector<Double> >& vectors,
                                            const UInt numNN,
                                            const Double alpha){

  CF::Matrix<Double> derivCoefVec;
  derivCoefVec.Resize(numNN,1);
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
      vals.Resize(numNN,1);
      vals[i][0] = vectors[i][0];
      if (trgGrid_->GetDim() == 2){
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
  //ALoc.Invert(invAloc);

  // coefficient matrix (coloumn nr. corresponding to the spatial dimension)
  derivCoefVec = ALoc * derivVec;
  vals.Transpose(temp);

  vec = temp * derivCoefVec;

}


void Lighthill::CalcLocCurlDerivativeCoefs(CF::Matrix<Double>& vec,
                                            const CF::Vector<Double>& globPoint,
                                            const CF::StdVector< Vector<Double> >& neighbors,
                                            const CF::StdVector< Double >& l2Distances,
                                            const CF::StdVector< Vector<Double> >& vectors,
                                            const UInt numNN,
                                            const Double alpha){

  CF::Matrix<Double> derivCoefVec;
  derivCoefVec.Resize(numNN,1);
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
        EXCEPTION("2D mesh and 3D-values!")
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
  tempvec.Resize(1,inDim_+1);
  tempvec = temp * derivCoefVec;

  //now we have to combine the tempvec-entries in order to obtain the rotor
  if(trgGrid_->GetDim() == 2){
    vec[0][0] = 0.0;
    vec[0][1] = 0.0;
    vec[0][2] = tempvec[1][0] - tempvec[0][1];
  }else{
    vec[0][0] = tempvec[2][1] - tempvec[1][2];
    vec[0][1] = tempvec[0][2] - tempvec[2][0];
    vec[0][2] = tempvec[1][0] - tempvec[0][1];
  }

}


void Lighthill::PrepareCalculation(){
  //1. Get get the source coordinates and the values, defined on those coordinates (Source...Src)
  //2. Get the target coordinates (trg)
  //3. Store for each trg element local Coordinates, elem number, volume, ...
  //4. Small clean-up


  std::cout << "\t ---> Lighthill preparing for interpolation" << std::endl;

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
              "Input of Lighthill- of LambVector-filter must be defined on nodes!! \n"
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
  StdVector<const CF::Elem*> tempElems;
  //mapping of global point targetCoords_ to local locPoints
  trgGrid_->GetElemsAtGlobalCoords(targetCoords_,locPoints, tempElems, lists, 1e-6, 1e-3);
  tempElems.Clear();

  std::cout << "\t\t 3/5 Generating differentiation info ..." << std::endl;
  derivData_.reserve(allTrgElems.GetSize());
  for(UInt aMatch = 0;aMatch < allTrgElems.GetSize();++aMatch){
    if(allTrgElems[aMatch]!= NULL){
      DifferentiationStruct newStruct;
      //newStruct.localCoords = locPoints[aMatch].coord;
      newStruct.tENum = allTrgElems[aMatch]->elemNum;
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
  str1::shared_ptr<EqnMapSimple> upMap = resultManager_->GetResultAdapter(upRes)->mapping;

  std::cout << "\t\t Differentiation prepared!" << std::endl;
}

ResultIdList Lighthill::SetUpstreamResults(){
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

void Lighthill::AdaptFilterResults(){
  //some checks
  ResultManager::ConstInfoPtr inInfo = resultManager_->GetExtInfo(upResIds[0]);
  if(!inInfo->isValid){
    EXCEPTION("Could not validate required input result \"" << inInfo->resultName << "\" from upstream filters.");
  }
  //we require mesh result input
  if(!inInfo->isMeshResult){
    EXCEPTION("Lighthill requires input to be defined on mesh");
  }

  //got the upstream result validated?
  if(!inInfo->isValid){
    EXCEPTION("Problem in filter pipeline detected. Differentiator input result \"" <<  inInfo->resultName << "\" could not be provided.")
  }


  resultManager_->CopyResultData(upResIds[0],filterResIds[0]);
  //but now, we need to overwrite some things
  resultManager_->SetRegionNames(filterResIds[0],this->trgRegions_);

  resultManager_->SetEntryType(filterResIds[0],ExtendedResultInfo::VECTOR);
  CF::StdVector<std::string> dofnames;
  dofnames.Push_back("x");
  dofnames.Push_back("y");
  dofnames.Push_back("z");

  resultManager_->SetDofNames(filterResIds[0],dofnames);

  resultManager_->SetDefOn(filterResIds[0],ExtendedResultInfo::ELEMENT);
  resultManager_->SetGrid(filterResIds[0],this->trgGrid_);
  resultManager_->SetMeshResult(filterResIds[0],true);

  //validate own result
  resultManager_->SetValid(filterResIds[0]);
}


}

