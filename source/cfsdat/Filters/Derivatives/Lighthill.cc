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
:AeroacousticBase(numWorkers,config,resMan){

#ifndef USE_CGAL
    EXCEPTION("CoefFunctionScatteredData needs to be compiled with USE_CGAL=ON!");
#endif

  this->filtStreamType_ = FIFO_FILTER;


    std::cout<<("============================================================ \n"
                "Make sure you set the ''density='' in the xml-scheme correctly, \n"
                "otherwise it is chosen to be 1.0 !! \n"
                "============================================================")<<std::endl;

    density_ = config->Get("scheme")->Get("density")->As<Double>();

    // ( LambVectorLighthill || LighthillSourceVector || LighthillSourceTerm )
    Form_ = config->Get("type")->As<std::string>();

    data_.dim = trgGrid_->GetDim();

    if(data_.dim == 3) EXCEPTION("SourceTerm-implementation for 3D not yet finished  !!!")


    if(config->Get("ResultList")->Get("vorticity")->Has("resultName") == false){
       std::cout<<"Vorticity is computed by CFSDat...Curl(u)"<<std::endl;
       externVorticity_ = false;
       data_.externVorticity = false;
     }else{
       externVorticity_ = true;
       data_.externVorticity = true;
     }

    std::string outRes = params_->Get("ResultList")->Get("outputQuantity")->Get("resultName")->As<std::string>();
    filtResNames.insert(outRes);

    //add velocity input result to manager
    this->res1Name = params_->Get("ResultList")->Get("velocity")->Get("resultName")->As<std::string>();
    upResNames.insert(res1Name);

    // if an external vorticity input exists, add it to manager
    if (externVorticity_ == true){
      this->res2Name = params_->Get("ResultList")->Get("vorticity")->Get("resultName")->As<std::string>();
      upResNames.insert(res2Name);
    }



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
  resultManager_->ActivateResult(upResIds[0]);

  // Activate extern vorticity-input, if provided
  if(externVorticity_ == true){
    resultManager_->SetTimeValue(upResIds[1],aTF);
    resultManager_->ActivateResult(upResIds[1]);
  }

  //now we call for upstream data in each source
  CF::StdVector< str1::shared_ptr<BaseFilter> >::iterator srcIter =  sources_.Begin();
  for(; srcIter != sources_.End() ; srcIter++){
    // should we check here anything for success?
    (*srcIter)->Run();
  }


  //**************************** Setup-Phase **********************************


  // upMap stores the equation-information of the INPUT-DATA on the entities defined in
  // AdaptFilterResults(). Since our convention is to provide the filter with node results, the equations are
  // defined on nodes and the dimension is 2 for 2D and 3 for 3D
  str1::shared_ptr<EqnMapSimple> upMap = resultManager_->GetResultAdapter(upResIds[0])->mapping;

  // downMap stores the equation-information of the OUTPUT-DATA on the entities defined in
  // AdaptFilterResults(), this means:
  //! -) AeroacousticSource_LambVector: output is vector-valued and has the dimension of the grid (2 for 2D, 3 for 3D)
  //!    -) externVorticity == true: NODE-RESULTS
  //!    -) externVorticity == false: ELEMENT-RESULTS
  //! -) AeroacousticSource_LighthillSourceVector: output is vector-valued and has the dimension of the grid (2 for 2D, 3 for 3D)
  //!    -) ELEMENT-RESULTS no matter if externVorticity is provided or not
  //! -) AeroacousticSource_LighthillSourceTerm (only 2D yet): physically it's a scalar but due to the
  //!                                                          result-managing it is vector valued with the scalar
  //!                                                          quantity in x-direction. to extract it and get a real
  //!                                                          scalar value, use the filter PostLighthillSourceTerm
  //!    -) ELEMENT-RESULTS no matter if externVorticity is provided or not
  str1::shared_ptr<EqnMapSimple> downMap = resultManager_->GetResultAdapter(filterResIds[0])->mapping;


  CF::StdVector<UInt> eqnNums;
  /// this is the vector, which will be filled with the wanted quantity
  Vector<Double>& returnVec = resultManager_->GetResultVector<Double>(filterResIds[0],eqnNums);
  returnVec.Init();

  // vector, containing the source data values
  // We can have two inputs... velocity and vorticity BUT we assume that both are defined on the same mesh !!
  data_.U = resultManager_->GetResultVector<Double>(upResIds[0],eqnNums);


  UInt dim = data_.dim;
  if( dim == 2 && data_.U.GetSize() != dim  * sourceCoords_.GetSize()){
    // 2D case
    std::cout<<"Length of velocity-input vector:"<<data_.U.GetSize()<<"\t Length of sourceCoords:"<<sourceCoords_.GetSize()<<std::endl;
    EXCEPTION("Velocity input result vector is probably not 2D, at least it another length than number of source coordinates times dimension");
  }
  if( dim == 2 && data_.U.GetSize() != dim * sourceCoords_.GetSize()){
    //3D case
    std::cout<<"Length of velocity-input vector:"<<data_.U.GetSize()<<"\t Length of sourceCoords:"<<sourceCoords_.GetSize()<<std::endl;
    EXCEPTION("Velocity input result vector is probably not 3D, at least it another length than number of source coordinates times dimension");
  }


  //**************************** Specialization-Phase **********************************
  Vector<Double> tempRetVec;
  if(Form_ == "AeroacousticSource_LambVector"){LambVector(tempRetVec);}
  if(Form_ == "AeroacousticSource_LighthillSourceVector"){LighthillSourceVector(tempRetVec);}
  if(Form_ == "AeroacousticSource_LighthillSourceTerm"){LighthillSourceTerm(tempRetVec, upMap);}


  // Check if the temporary result vector tempResVec has the same/correct size as the
  // provided returnVec
  if(tempRetVec.GetSize() != returnVec.GetSize()){
    std::cout<<"size of computed result-vector:"<<tempRetVec.GetSize()<<std::endl;
    std::cout<<"size of provided result-vector:"<<returnVec.GetSize()<<std::endl;
    EXCEPTION("Size of computed result-vector does not coincide with the provided result-vector ... this should not happen ... indicates an implementation-error!");
  }else{
    returnVec = tempRetVec;
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



void Lighthill::LambVector(Vector<Double>& tempRetVec){
  // Remember the output is defined on:
  // -) NODE if an extern vorticity is provided
  // -) ELEMENT if no extern vorticity is provided
  // This means every downMap with filterResIds[0] gives you equations for
  // nodes/elements and VECTOR-values!!!

  // fill scatteredData vector
  CF::StdVector< Vector<Double> >  scatteredData;
  Vector<Double> vec; //neglect vec, not used here
  UInt t; //neglect dimension output, we have our own
  scatteredData.Resize(derivDataNode_.size());
  FillScatteredDataVec(scatteredData, vec, data_.U, sourceCoords_, t);
  data_.scatteredData = scatteredData;


  /**********************************************************
   * To be computed:
   * rho * (Curl(u) x u)
   ***********************************************************/

  // Term 2: rho * (Curl(u) x u)
  Vector<Double> OmegaCrossU;
  if( externVorticity_ == false){
    // Calculate the curl internally
    // here we have to use the derivative data of the element centroids derivDataElem_
    tempRetVec.Resize(data_.dim * derivDataElem_.size());
    Vector<Double> CurlU;
    if(data_.dim == 2){
      CalcCurlU_2D(CurlU, scatteredData, derivDataElem_, sourceCoords_, filterResIds[0], data_.dim);
    }else{
      CalcCurlU_3D(CurlU, scatteredData, derivDataElem_, sourceCoords_, filterResIds[0], data_.dim);
    }

    // since the curlU is now defined on element centroids, we also have
    // to restrict the velocity to element centroids (data_.Uelem)
    data_.Uelem.Resize(data_.dim * derivDataElem_.size());
    Node2Cell(data_.Uelem, filterResIds[0], data_.U, derivDataElem_);

    data_.Om = CurlU;

    OmegaVectorProductU_2D(data_, derivDataElem_, OmegaCrossU);
  }else{
    // use the provided vorticity
    tempRetVec.Resize(data_.dim * derivDataNode_.size());
    CF::StdVector<UInt> eqnNums;
    data_.Om = resultManager_->GetResultVector<Double>(upResIds[1],eqnNums);
    OmegaVectorProductU_2D(data_, derivDataNode_, OmegaCrossU);
  }

  tempRetVec = OmegaCrossU * density_;
}




void Lighthill::LighthillSourceVector(Vector<Double>& tempRetVec){
  // Remember the output is defined on ELEMENT and is VECTOR-valued!!
  // This means every downMap with filterResIds[0] gives you equations for
  // elements and VECTOR-values!!!

  tempRetVec.Resize(data_.dim * derivDataElem_.size()); //defined on elements


  // fill scatteredData vector
  CF::StdVector< Vector<Double> >  scatteredData;
  Vector<Double> vec; //neglect vec, not used here
  UInt t; //neglect dimension output, we have our own
  scatteredData.Resize(derivDataNode_.size());
  FillScatteredDataVec(scatteredData, vec, data_.U, sourceCoords_, t);
  data_.scatteredData = scatteredData;


  /**********************************************************
   * To be computed:
   * Grad(1/2 u*u)  + rho * (Curl(u) x u)
   * ### Term 1 ### + ##### Term 2  #####
   ***********************************************************/

  // Term 1: Grad(1/2 u*u)
  Vector<Double> GradOfScalar;
  CalcGradUScalarU(GradOfScalar, data_, derivDataElem_, sourceCoords_, filterResIds[0]);
  GradOfScalar.ScalarMult(0.5);

  // Term 2: rho * (Curl(u) x u)
  Vector<Double> OmegaCrossU;
  if( externVorticity_ == false){
    // Calculate the curl internally
    // here we have to use the derivative data of the element centroids derivDataElem_
    Vector<Double> CurlU;
    if(data_.dim == 2){
      CalcCurlU_2D(CurlU, scatteredData, derivDataElem_, sourceCoords_, filterResIds[0], data_.dim);
    }else{
      CalcCurlU_3D(CurlU, scatteredData, derivDataElem_, sourceCoords_, filterResIds[0], data_.dim);
    }

    // since the curlU is now defined on element centroids, we also have
    // to restrict the velocity to element centroids (data_.Uelem)
    data_.Uelem.Resize(data_.dim * derivDataElem_.size());
    Node2Cell(data_.Uelem, filterResIds[0], data_.U, derivDataElem_);

    data_.Om = CurlU;

    OmegaVectorProductU_2D(data_, derivDataElem_, OmegaCrossU);
  }else{
    // use the provided vorticity
    CF::StdVector<UInt> eqnNums;
    data_.Om = resultManager_->GetResultVector<Double>(upResIds[1],eqnNums);
    OmegaVectorProductU_2D(data_, derivDataNode_, OmegaCrossU);
  }


  // Term 1 + Term 2
  if(externVorticity_ == true){
    // OmegaCrossU is defined on nodes but GradOfScalar in defined on elements
    // so we restrict OmegaCrossU also to elements
    Vector<Double> OmegaCrossUNodes;
    OmegaCrossUNodes.Resize(data_.dim * derivDataElem_.size());
    Node2Cell(OmegaCrossUNodes, filterResIds[0], OmegaCrossU, derivDataElem_);
    tempRetVec = GradOfScalar + OmegaCrossUNodes * density_;
  }else{
    tempRetVec = GradOfScalar + OmegaCrossU * density_;
  }
}



void Lighthill::LighthillSourceTerm(Vector<Double>& tempRetVec, const str1::shared_ptr<EqnMapSimple>& upMap){
  // Remember the output is defined on:
  // -) ELEMENT no matter if externVorticity is provided
  // This means every downMap with filterResIds[0] gives you equations for
  // elements and VECTOR-values!!! To clarify...actually it's a scalar-quantity but
  // we have to write it as a vector with the scalar quantity in x-direction and extract
  // it afterwards with another filter PostLighthillSourceTerm

  tempRetVec.Resize(data_.dim * derivDataElem_.size()); //defined on elements


  // fill scatteredData vector
  CF::StdVector< Vector<Double> >  scatteredData;
  Vector<Double> vec; //neglect vec, not used here
  UInt t; //neglect dimension output, we have our own
  scatteredData.Resize(derivDataNode_.size());
  FillScatteredDataVec(scatteredData, vec, data_.U, sourceCoords_, t);
  data_.scatteredData = scatteredData;


  /**********************************************************
   * To be computed:
   * Div( Grad(1/2 u*u)  + rho * (Curl(u) x u) )
   * Div( ### Term 1 ### + ##### Term 2  ##### )
   ***********************************************************/

  // Term 1: Grad(1/2 u*u)
  Vector<Double> GradOfScalar;
  CalcGradUScalarU(GradOfScalar, data_, derivDataElem_, sourceCoords_, filterResIds[0]);
  GradOfScalar.ScalarMult(0.5);

  // Term 2: rho * (Curl(u) x u)
  Vector<Double> OmegaCrossU;
  if( externVorticity_ == false){
    // Calculate the curl internally
    // here we have to use the derivative data of the element centroids derivDataElem_
    Vector<Double> CurlU;
    if(data_.dim == 2){
      CalcCurlU_2D(CurlU, scatteredData, derivDataElem_, sourceCoords_, filterResIds[0], data_.dim);
    }else{
      CalcCurlU_3D(CurlU, scatteredData, derivDataElem_, sourceCoords_, filterResIds[0], data_.dim);
    }

    // since the curlU is now defined on element centroids, we also have
    // to restrict the velocity to element centroids (data_.Uelem)
    data_.Uelem.Resize(data_.dim * derivDataElem_.size());
    Node2Cell(data_.Uelem, filterResIds[0], data_.U, derivDataElem_);
    data_.Om = CurlU;

    OmegaVectorProductU_2D(data_, derivDataElem_, OmegaCrossU);
  }else{
    // use the provided vorticity
    CF::StdVector<UInt> eqnNums;
    data_.Om = resultManager_->GetResultVector<Double>(upResIds[1],eqnNums);
    OmegaVectorProductU_2D(data_, derivDataNode_, OmegaCrossU);
  }


  // Term 1 + Term 2
  Vector<Double> T;
  if(externVorticity_ == true){
    // OmegaCrossU is defined on nodes but GradOfScalar in defined on elements
    // so we restrict OmegaCrossU also to elements
    Vector<Double> OmegaCrossUNodes;
    OmegaCrossUNodes.Resize(data_.dim * derivDataElem_.size());
    Node2Cell(OmegaCrossUNodes, filterResIds[0], OmegaCrossU, derivDataElem_);
    T = GradOfScalar + OmegaCrossUNodes * density_;
  }else{
    T = GradOfScalar + OmegaCrossU * density_;
  }


  // T is element result, to calculate the divergence we need NODE-results,
  // this means we have to interpolate from elem -> node
  //TODO Incorporate also other elem -> node interpolators (e.g. Cell2Node or RBF, ...)

  //Prepare Interpolation
  Vector<Double> tempV;
  // fill scatteredData vector with T (LighthillSourceVector)
  scatteredData.Clear(false);
  vec.Clear(false);
  scatteredData.Resize(derivDataElem_.size());

  FillScatteredDataVec(scatteredData, vec, T, targetCoords_, t);


  // Interpolation
  tempV.Resize(data_.dim * derivDataNode_.size());
  NearestNeighbourInterpolation(tempV, scatteredData, vec, upMap, targetCoords_, sourceCoords_, derivDataElem_, trgGrid_, t, 10, 5);

  // now divergence
  Vector<Double> test;
  CalcDivergence2D(test, tempV, sourceCoords_, derivDataElem_);


  tempRetVec = ScalarToTwoD(test);

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


  std::cout << "\t\t 1/4 Obtaining source coordinates " << std::endl;
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
  std::cout << "\t\t 2/4 Obtaining target coordinates" << std::endl;
  // if Lamb-vector and extern-vorticity - then output is defined on nodes
  // else output is defined on elements

  if(Form_ == "AeroacousticSource_LambVector" && externVorticity_ == true){
    // target (output) is solely defined on nodes
    targetCoords_.Resize(allTrgNodes.GetSize());
    //for(UInt nIter=0; nIter < allTrgNums.GetSize(); ++nIter){
    for(UInt nIter=0; nIter < allTrgNodes.GetSize(); ++nIter){
      CF::Vector<Double> SCoord;
      //trgGrid_->GetElemCentroid(SCoord, allTrgElemNums[nIter], true);
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
  }else{
    // target (output) is solely defined on elements
    targetCoords_.Resize(allTrgElemNums.GetSize());
    //for(UInt nIter=0; nIter < allTrgNums.GetSize(); ++nIter){
    for(UInt nIter=0; nIter < allTrgElemNums.GetSize(); ++nIter){
      CF::Vector<Double> SCoord;
      trgGrid_->GetElemCentroid(SCoord, allTrgElemNums[nIter], true);
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
  }


  CF::StdVector< LocPoint > locPoints;
  StdVector<const CF::Elem*> tempElems;
  //mapping of global point targetCoords_ to local locPoints
  trgGrid_->GetElemsAtGlobalCoords(targetCoords_,locPoints, tempElems, lists, 1e-6, 1e-3);
  //tempElems.Clear();

  std::cout << "\t\t 3/4 Generating differentiation info ..." << std::endl;
  derivDataNode_.reserve(allTrgNodes.GetSize());
  for(UInt aMatch = 0;aMatch < allTrgNodes.GetSize();++aMatch){
    //if(tempElems[aMatch]!= NULL){
      QuantityStruct newStruct;
      //newStruct.localCoords = locPoints[aMatch].coord;
      //newStruct.trgElemNum = tempElems[aMatch]->elemNum;
      derivDataNode_.push_back(newStruct);
    //}
  }


  StdVector<UInt> tempNodeNums;
  StdVector<UInt> sEqn;
  str1::shared_ptr<EqnMapSimple> upMap = resultManager_->GetResultAdapter(upRes)->mapping;
  derivDataElem_.reserve(tempElems.GetSize());
  for(UInt aMatch = 0;aMatch < tempElems.GetSize();++aMatch){
    if(tempElems[aMatch]!= NULL){
      QuantityStruct newStruct;
      newStruct.localCoords = locPoints[aMatch].coord;
      newStruct.trgElemNum = tempElems[aMatch]->elemNum;
      trgGrid_->GetElemNodes(tempNodeNums, tempElems[aMatch]->elemNum );
      newStruct.tNNum.Resize(tempNodeNums.GetSize());
      newStruct.tNNum = tempNodeNums;
      for(UInt n = 0; n < tempNodeNums.GetSize(); ++n){
        upMap->GetEquation(sEqn, tempNodeNums[n], ExtendedResultInfo::NODE);
        for(UInt d = 0; d < sEqn.GetSize(); ++d){
          newStruct.srcEqn.Push_back(sEqn[d]);
        }
      }
      derivDataElem_.push_back(newStruct);
    }
  }


  std::cout << "\t\t 4/4 Clear generated temporary data storage ..." << std::endl;
  allTrgElems.Clear(false);
  allSrcElems.Clear(false);
  allSrcNodes.Clear(false);
  allTrgElemNums.Clear(false);
  allTrgNodes.Clear(false);


  std::cout << "\t\t Differentiation prepared!" << std::endl;
}

ResultIdList Lighthill::SetUpstreamResults(){
  ResultIdList generated;
  //we should only have one filter Result
  CF::StdVector<uuids::uuid>::iterator aIt = filterResIds.Begin();
  std::string filterResName = resultManager_->GetExtInfo(*aIt)->resultName;

  //add velocity input result to manager
  this->res1Name = params_->Get("ResultList")->Get("velocity")->Get("resultName")->As<std::string>();
  upResNames.insert(res1Name);
  res1Id = resultManager_->AddResult(res1Name,this->filterTag_);
  //set the timeline of upstream data
  resultManager_->SetTimeLine(res1Id,(*resultManager_->GetExtInfo(*aIt)->timeLine.get()));
  generated.Push_back(res1Id);

  // if an external vorticity input exists, add it to manager
  if (externVorticity_ == true){
    this->res2Name = params_->Get("ResultList")->Get("vorticity")->Get("resultName")->As<std::string>();
    upResNames.insert(res2Name);
    res2Id = resultManager_->AddResult(res2Name,this->filterTag_);
    resultManager_->SetTimeLine(res2Id,(*resultManager_->GetExtInfo(*aIt)->timeLine.get()));
    generated.Push_back(res2Id);
  }
  return generated;
}

void Lighthill::AdaptFilterResults(){
  //some checks
  ResultManager::ConstInfoPtr inInfo1 = resultManager_->GetExtInfo(upResIds[0]);
  if(!inInfo1->isValid){
    EXCEPTION("Could not validate required input result \"" << inInfo1->resultName << "\" from upstream filters.");
  }

  //we require mesh result input
  if(!inInfo1->isMeshResult){
    EXCEPTION("Lighthill/Lamb filter requires input to be defined on mesh");
  }

  //input must be node-based
  if (resultManager_->GetDefOn(upResIds[0]) == ExtendedResultInfo::ELEMENT){
    EXCEPTION("Lighthill/Lamb filter requires velocity to be defined on nodes!");
  }

  //This does not work because also if the upstream result is a scalar, we get x,y,z as dofnames
  //std::cout<<"DOFNames"<<resultManager_->GetExtInfo(upResIds[0])->dofNames<<std::endl;

  if (externVorticity_ == true){
      ResultManager::ConstInfoPtr inInfo2 = resultManager_->GetExtInfo(upResIds[1]);
      if(!inInfo2->isValid){
        EXCEPTION("Could not validate required input result \"" << inInfo2->resultName << "\" from upstream filters.");
      }
      //we require mesh result input
      if(!inInfo2->isMeshResult){
        EXCEPTION("Lighthill requires input to be defined on mesh");
      }
      //input must be node-based
      if (resultManager_->GetDefOn(upResIds[1]) == ExtendedResultInfo::ELEMENT){
        EXCEPTION("Lighthill/Lamb filter requires vorticity to be defined on nodes!");
      }
  }


  resultManager_->CopyResultData(upResIds[0],filterResIds[0]);

  //but now, we need to overwrite some things
  resultManager_->SetRegionNames(filterResIds[0],this->trgRegions_);


  resultManager_->SetEntryType(filterResIds[0],ExtendedResultInfo::VECTOR);
  CF::StdVector<std::string> dofnames;
  if(trgGrid_->GetDim() == 2){
    dofnames.Push_back("x");
    dofnames.Push_back("y");
  }else{
    dofnames.Push_back("x");
    dofnames.Push_back("y");
    dofnames.Push_back("z");
  }
  resultManager_->SetDofNames(filterResIds[0],dofnames);

  // choose if we have element or node-results
  if(externVorticity_ == true && Form_ == "AeroacousticSource_LambVector"){
    resultManager_->SetDefOn(filterResIds[0],ExtendedResultInfo::NODE);
  }else{
    resultManager_->SetDefOn(filterResIds[0],ExtendedResultInfo::ELEMENT);
  }

  resultManager_->SetGrid(filterResIds[0],this->trgGrid_);
  resultManager_->SetMeshResult(filterResIds[0],true);

  //validate own result
  resultManager_->SetValid(filterResIds[0]);
}


}

