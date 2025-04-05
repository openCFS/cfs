// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionCurvilinearPML.hh
 *       \brief    Implements a custom CoefFunction for a curvilinear PML that is supposed to 
 *                 function in arbitrarily-shaped, convex domains. 
 *                 The CoefFunctionCurvilinearPML requires geometric properties (e.g, normal 
 *                 vectors, principal-direction vectors and principal curvatures) that must be 
 *                 provided or can be computed via the layerGeneration property.
 *       \date     Mar 2023
 *       \author   pheidegger
 */
//================================================================================================

#include "CoefFunctionCurvilinearPML.hh"

namespace CoupledField{
  //================================================================================================
  // Curvilinear PML
  //================================================================================================

  template<typename T>
  CoefFunctionCurvilinearPML<T>::CoefFunctionCurvilinearPML(PtrParamNode pmlDef, PtrCoefFct speedOfSound, shared_ptr<EntityList> EntList,
                        StdVector<RegionIdType> pdeDomains, OutputType outputType) : CoefFunctionPMLBase<T>(pmlDef, speedOfSound, EntList, pdeDomains) {
    // set name and type
    this->name_ = "CoefFunctionCurvilinearPML";
    this->formulationType_ = CURVILINEAR;
    // get the grid pointer
    grid_ =  this->entities_[0]->GetGrid();
    // set the type of output
    outputType_ = outputType;
    // assign the dimType_ of the coefFunction according to the output type
    SetDimType();
    // read from PML node
    ReadDataPML(pmlDef);
    // check for valid declaration in XML
    CheckForInvalidParams(pmlDef);
    // check if PML layer was generated automatically
    CheckForLayerGenerationNode(pmlDef);
    // trigger computation of geometry data and get pointer to it (is only computed once)
    grid_->GetGeometryOnRegionNodes(ptrNodeGeom_, layerGenNode_, ifSurfRegion_);
    // get connected iso surfaces of the volume region
    grid_->GetConnectedSurfaceRegions(surfRegions_, volRegion_);
    // get the lowest and highest node numbers of the automatically generated iso surfaces
    // these will be needed for the index searching in GetIdxByNodeId()
    minIds_.Resize(numLayers_+1);
    maxIds_.Resize(numLayers_+1);
    StdVector<UInt> nodeList;
    for (UInt iLayers = 0; iLayers <= numLayers_; ++iLayers) {
      grid_->GetNodesByRegion(nodeList, surfRegions_[iLayers]);
      minIds_[iLayers] = nodeList[0];
      maxIds_[iLayers] = nodeList[numSurfNodes_-1];
      nodeList.Clear();
      //std::cout << "Layer " << ii << ":\n" << "min=" << nodeList[0] << "\n" << "max=" << nodeList[numSurfNodes_-1] << "\n";
    }
    // get nodal distances to the PML interface and store in a StdVector
    GetThicknessOnNodes();
    // set identity operators for mapping nodal values to lpm
    CreateMappingOperators();
  }


  template<typename T>
  CoefFunctionCurvilinearPML<T>::~CoefFunctionCurvilinearPML() { };


  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetTensor(Matrix<Complex>& tensor, const LocPointMapped& lpm ) {
    switch(outputType_) {
      case OutputType::TENSOR: {
        // get a pointer to the considered element from the lpm
        const Elem* ptrElem = NULL;
        ptrElem = lpm.ptEl;
        // the global element ID
        UInt elemId = ptrElem->GetElemNum();
        // get the nodes connected to the element
        StdVector<UInt> nodeIds;
        grid_->GetElemNodes(nodeIds, elemId);

        // index of a node id
        UInt nodeIdx = 0;
        // index of the corresponding node on the interface
        UInt nodeIdxIface = 0;
        // index for vectorial quantities on nodes stored in Vector
        UInt vecIdx;
        // number of nodes in element
        UInt numElemNodes = nodeIds.GetSize();

        // get the ElemShapeMap of the current element from the grid (attention! updated geometry is not set!)
        shared_ptr<ElemShapeMap> ptrEsm = this->grid_->GetElemShapeMap(ptrElem);

        // get Base Fe, which provides the shape functions for interpolating with the identity operator 
        BaseFE* ptrFe = ptrEsm->GetBaseFE();

        if( this->dim_ == 2 )
        {
          // Initialize stacked geometry tensors
          Vector<Double> normMats(numElemNodes * pow(this->dim_,2));
          Vector<Double> tmaxMats(numElemNodes * pow(this->dim_,2));

          // Scalar quantities
          Vector<Double> maxPrincCurv(numElemNodes);
          Vector<Double> dist(numElemNodes);

          // Interpolated tensor storage
          Vector<Double> N(this->dim_);
          Vector<Double> Tmax(this->dim_);

          // Scalars for interpolation
          Vector<Double> kmax(1);
          Vector<Double> d(1);

          // Loop over element nodes and compute interpolated geometry on lpm
          for (UInt iNodes = 0; iNodes < numElemNodes; ++iNodes) 
          {
              nodeIdx = GetIdxByNodeId(nodeIds[iNodes], (nodeIdx < numElemNodes) ? 0 : nodeIdx - numElemNodes);
              nodeIdxIface = nodeIdx % numSurfNodes_;

              // Compute dyadic products for normal and tangential vectors
              for (UInt iRow = 0; iRow < this->dim_; ++iRow) 
              {
                  for (UInt iCol = 0; iCol < this->dim_; ++iCol) 
                  {
                      vecIdx = iCol + this->dim_ * iRow + pow(this->dim_, 2) * iNodes;
                      normMats[vecIdx] = ptrNodeGeom_->normalVectors_[nodeIdxIface][iRow] * ptrNodeGeom_->normalVectors_[nodeIdxIface][iCol];
                      tmaxMats[vecIdx] = ptrNodeGeom_->maxPrincipalVectors_[nodeIdxIface][iRow] * ptrNodeGeom_->maxPrincipalVectors_[nodeIdxIface][iCol];
                  }
              }

              // Get scalar values
              maxPrincCurv[iNodes] = ptrNodeGeom_->maxPrincipalCurvatures_[nodeIdxIface];
              dist[iNodes] = thicknessOnNodes_[nodeIdx];
          }

          // Perform interpolation
          this->tensorMappingOperator_->ApplyOp(N, lpm, ptrFe, normMats);
          this->tensorMappingOperator_->ApplyOp(Tmax, lpm, ptrFe, tmaxMats);
          this->scalarMappingOperator_->ApplyOp(kmax, lpm, ptrFe, maxPrincCurv);
          this->scalarMappingOperator_->ApplyOp(d, lpm, ptrFe, dist);

          // Compute damping function
          Double dampFunc = this->dampFunction_->ComputeFactor(d[0], layerThickness_);
          Double intDampFunc = this->dampFunction_->ComputeIntegralFactor(d[0], layerThickness_);

          // Get speed of sound at current location
          Double sos;
          this->speedOfSound_->GetScalar(sos, lpm);

          // Compute the wave number
          Double K = this->omega_ / sos;

          // Compute helper functions
          Vector<Double> h(this->dim_);
          h[0] = dampFunc;
          h[1] = intDampFunc * kmax[0] / (1.0 + d[0] * kmax[0]);

          // Compute inverse of the metric tensor
          Vector<Complex> s(this->dim_);
          for (UInt iDim = 0; iDim < this->dim_; ++iDim) {
              s[iDim] = Complex(1.0, 0.0) / Complex(1.0, -h[iDim] / K);
          }

          // Compute the tensor and unstack
          tensor.Resize(this->dim_, this->dim_);
          for (UInt iRow = 0; iRow < this->dim_; ++iRow) {
              for (UInt iCol = 0; iCol < this->dim_; ++iCol) {
                  vecIdx = iCol + this->dim_ * iRow;
                  tensor[iRow][iCol] = s[0] * N[vecIdx] + s[1] * Tmax[vecIdx];
              }
          }
        }
        else if( this->dim_ == 3 )
        {

          // initialize stacked geometry tensors (needs to be stacked as interpolation requires it)
          // vector entries: [N1_T11, N1_T12,..., N1_T33, N2_T11, ...], 
          // where Nx is the node and Txx are entries in respective tensor
          Vector<Double> normMats(numElemNodes * pow(this->dim_,2));
          Vector<Double> tminMats(numElemNodes * pow(this->dim_,2));
          Vector<Double> tmaxMats(numElemNodes * pow(this->dim_,2));
          // vectors for scalar quantities
          Vector<Double> minPrincCurv(numElemNodes);
          Vector<Double> maxPrincCurv(numElemNodes);
          Vector<Double> dist(numElemNodes);

          // variable to store interpolated tensors
          Vector<Double> N(this->dim_);
          Vector<Double> Tmin(this->dim_);
          Vector<Double> Tmax(this->dim_);
          // create 1 element vectors for scalars because ApplyOp only takes Vectors
          Vector<Double> kmin(1);
          Vector<Double> kmax(1);
          Vector<Double> d(1);

          // loop over element nodes and compute interpolated geometry on lpm
          for (UInt iNodes = 0; iNodes < numElemNodes; ++iNodes) {
            // get node idx and provide an initial guess for faster search
            if (nodeIdx < numElemNodes) {
              nodeIdx = GetIdxByNodeId(nodeIds[iNodes], 0);
            } else {
              nodeIdx = GetIdxByNodeId(nodeIds[iNodes], nodeIdx-numElemNodes);
            }
            nodeIdxIface = nodeIdx % numSurfNodes_;
            //std::cout << nodeIdx << "\n" << nodeIdxIface;
            
            // store dyadic products of current node in a stacked vector
            for (UInt iRow = 0; iRow < this->dim_; ++iRow) {
              for (UInt iCol = 0; iCol < this->dim_; ++iCol) {
                vecIdx = iCol + this->dim_ * iRow + pow(this->dim_,2) * iNodes;
                normMats[vecIdx] = ptrNodeGeom_->normalVectors_[nodeIdxIface][iRow] * ptrNodeGeom_->normalVectors_[nodeIdxIface][iCol];
                tminMats[vecIdx] = ptrNodeGeom_->minPrincipalVectors_[nodeIdxIface][iRow] * ptrNodeGeom_->minPrincipalVectors_[nodeIdxIface][iCol];
                tmaxMats[vecIdx] = ptrNodeGeom_->maxPrincipalVectors_[nodeIdxIface][iRow] * ptrNodeGeom_->maxPrincipalVectors_[nodeIdxIface][iCol];
              }
            }
            // get scalar values
            minPrincCurv[iNodes] = ptrNodeGeom_->minPrincipalCurvatures_[nodeIdxIface];
            maxPrincCurv[iNodes] = ptrNodeGeom_->maxPrincipalCurvatures_[nodeIdxIface];
            dist[iNodes] = thicknessOnNodes_[nodeIdx];
          }

          // perform interpolation
          this->tensorMappingOperator_->ApplyOp(N,lpm,ptrFe,normMats);
          this->tensorMappingOperator_->ApplyOp(Tmin,lpm,ptrFe,tminMats);
          this->tensorMappingOperator_->ApplyOp(Tmax,lpm,ptrFe,tmaxMats);
          this->scalarMappingOperator_->ApplyOp(kmin,lpm,ptrFe,minPrincCurv);
          this->scalarMappingOperator_->ApplyOp(kmax,lpm,ptrFe,maxPrincCurv);
          this->scalarMappingOperator_->ApplyOp(d,lpm,ptrFe,dist);

          // get the damping function and its integral at the mapped distance
          Double dampFunc = this->dampFunction_->ComputeFactor(d[0], layerThickness_);
          Double intDampFunc = this->dampFunction_->ComputeIntegralFactor(d[0], layerThickness_);

          // get speed of sound at current lpm
          Double sos;
          this->speedOfSound_->GetScalar(sos,lpm);

          // compute the current wave number
          Double K = this->omega_ / sos;

          // compute helper functions
          Vector<Double> h = Vector<Double>(this->dim_);
          h[0] = dampFunc;
          h[1] = intDampFunc * kmin[0] / (1.0 + d[0] * kmin[0]);
          h[2] = intDampFunc * kmax[0] / (1.0 + d[0] * kmax[0]);
          // get the entries of (I+jD)^-1 (is named s but resembles 1/s_i)
          Vector<Complex> s = Vector<Complex>(3);
          for (UInt iDim = 0; iDim < 3; ++iDim) {
            s[iDim] = Complex(1.0,0.0) / Complex(1.0, -h[iDim]/K);
          }
          // compute the tensor and unstack
          tensor.Resize(this->dim_, this->dim_);
          for (UInt iRow = 0; iRow < this->dim_; ++iRow) {
            for (UInt iCol = 0; iCol < this->dim_; ++iCol) {
              vecIdx = iCol + this->dim_ * iRow;
              tensor[iRow][iCol] = s[0]*N[vecIdx] + s[1]*Tmin[vecIdx] + s[2]*Tmax[vecIdx];
            }
          }
        }
        break;
      }
      default:
        EXCEPTION("CoefFunctionCurvilinearPML::GetTensor(Complex...) is used for OutputType: " <<
                  "TENSOR only.");
    }
  };


  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetScalar(Complex& val, const LocPointMapped& lpm ) {
    switch(outputType_) {
      case OutputType::DETERMINANT: {
        // here we only need curvarute and distances...
        // get a pointer to the considered element from the lpm
        const Elem* ptrElem = NULL;
        ptrElem = lpm.ptEl;
        // the global element ID
        UInt elemId = ptrElem->GetElemNum();
        // get the nodes connected to the element
        StdVector<UInt> nodeIds;
        grid_->GetElemNodes(nodeIds, elemId);

        // index of a node id
        UInt nodeIdx;
        // index of the corresponding node on the interface
        UInt nodeIdxIface = 0;
        // number of nodes in element
        UInt numElemNodes = nodeIds.GetSize();

        // get the ElemShapeMap of the current element from the grid (attention! updated geometry is not set!)
        shared_ptr<ElemShapeMap> ptrEsm = this->grid_->GetElemShapeMap(ptrElem);

        // get Base Fe, which provides the shape functions for interpolating with the identity operator 
        BaseFE* ptrFe = ptrEsm->GetBaseFE();

        // in 2D we only use the tmin_ as tangential vector and kmin_ as curvature
        if (this->dim_ == 2)
        {
          // Vectors to store geometry data at nodes
          Vector<Double> maxPrincCurv(numElemNodes);
          Vector<Double> dist(numElemNodes);

          // Loop over element nodes and get necessary quantities
          for (UInt iNodes = 0; iNodes < numElemNodes; ++iNodes) {
            nodeIdx = GetIdxByNodeId(nodeIds[iNodes], 0);
            nodeIdxIface = nodeIdx % numSurfNodes_;
            
            maxPrincCurv[iNodes] = ptrNodeGeom_->maxPrincipalCurvatures_[nodeIdxIface];
            dist[iNodes] = thicknessOnNodes_[nodeIdx];
          }

          // Create single-element vectors for interpolation since ApplyOp only takes vectors
          Vector<Double> kmax(1);
          Vector<Double> d(1);

          // Interpolate curvature and distance to the local point
          this->scalarMappingOperator_->ApplyOp(kmax, lpm, ptrFe, maxPrincCurv);
          this->scalarMappingOperator_->ApplyOp(d, lpm, ptrFe, dist);

          // Compute the damping function and its integral
          Double dampFunc = this->dampFunction_->ComputeFactor(d[0], layerThickness_);
          Double intDampFunc = this->dampFunction_->ComputeIntegralFactor(d[0], layerThickness_);

          // Get speed of sound at current lpm
          Double sos;
          this->speedOfSound_->GetScalar(sos, lpm);

          // Compute wave number
          Double K = this->omega_ / sos;

          // Compute helper functions
          Vector<Double> h(2);
          h[0] = dampFunc;
          h[1] = intDampFunc * kmax[0] / (1.0 + d[0] * kmax[0]);

          // Compute determinant in 2D:
          // The eigenvalues in 2D are of the form s_i = (1 + 1i/K * h[i])
          val = Complex(1.0, -h[0] / K) * Complex(1.0, -h[1] / K);
        } // dim_ == 2
        else if (this->dim_ == 3) {
          // vectors with extracted geometry data
          Vector<Double> minPrincCurv(numElemNodes);
          Vector<Double> maxPrincCurv(numElemNodes);
          // and on node-to-interface distances
          Vector<Double> dist(numElemNodes);

          // loop over element nodes and get quantities
          for (UInt iNodes = 0; iNodes < numElemNodes; ++iNodes) {
            nodeIdx = GetIdxByNodeId(nodeIds[iNodes],0);
            nodeIdxIface = nodeIdx % numSurfNodes_;
            minPrincCurv[iNodes] = ptrNodeGeom_->minPrincipalCurvatures_[nodeIdxIface];
            maxPrincCurv[iNodes] = ptrNodeGeom_->maxPrincipalCurvatures_[nodeIdxIface];
            dist[iNodes] = thicknessOnNodes_[nodeIdx];
          }

          // create 1 element vectors for scalars because ApplyOp only takes Vectors
          Vector<Double> kmin(1);
          Vector<Double> kmax(1);
          Vector<Double> d(1);

          // interpolate from nodes to the local point mapped
          this->scalarMappingOperator_->ApplyOp(kmin,lpm,ptrFe,minPrincCurv);
          this->scalarMappingOperator_->ApplyOp(kmax,lpm,ptrFe,maxPrincCurv);
          this->scalarMappingOperator_->ApplyOp(d,lpm,ptrFe,dist);

          // get the damping function and its integral at the mapped distance
          Double dampFunc = this->dampFunction_->ComputeFactor(d[0], layerThickness_);
          Double intDampFunc = this->dampFunction_->ComputeIntegralFactor(d[0], layerThickness_);

          // get speed of sound at current lpm
          Double sos;
          this->speedOfSound_->GetScalar(sos,lpm);

          // compute the current wave number
          Double K = this->omega_ / sos;

          // compute helper funcions
          Vector<Double> h = Vector<Double>(3);
          h[0] = dampFunc;
          h[1] = intDampFunc * kmin[0] / (1.0 + d[0] * kmin[0]);
          h[2] = intDampFunc * kmax[0] / (1.0 + d[0] * kmax[0]);

          // compute the determinant...
          // the eigenvalues are of the form s_i = (1 + 1i/K * h[i])
          // where h[i] contains the damping function and the metric correction for the curvature
          val = Complex(1.0, -h[0]/K) * Complex(1.0, -h[1]/K) * Complex(1.0, -h[2]/K);
        } // dim_ == 3
        break;
      }
      case OutputType::DISTANCE: {
        // actually here it would be good to only call the GetScalar(Double...) function.
        // however, CFS calls the GetScalar(Complex...) functions even though that the 
        // template is set to Double.
        Double temp;
        GetScalar(temp, lpm);
        val = Complex(temp, 0.0);
        break;
      }
      default:
        EXCEPTION("CoefFunctionCurvilinearPML::GetScalar(Complex...) is used for OutputType: " <<
                  "DETERMINANT only.");
    }
  };


  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetScalar(Double& val, const LocPointMapped& lpm ) {
    switch(outputType_) {
      case OutputType::DISTANCE: {
        // here we need to interpolate only the nodal distances to the interface...
        // get a pointer to the considered element from the lpm
        const Elem* ptrElem = NULL;
        ptrElem = lpm.ptEl;
        // the global element ID
        UInt elemId = ptrElem->GetElemNum();
        // get the nodes connected to the element
        StdVector<UInt> nodeIds;
        grid_->GetElemNodes(nodeIds, elemId);

        // index of a node id
        UInt nodeIdx;
        // number of nodes in element
        UInt numElemNodes = nodeIds.GetSize();

        // get the ElemShapeMap of the current element from the grid (attention! updated geometry is not set!)
        shared_ptr<ElemShapeMap> ptrEsm = this->grid_->GetElemShapeMap(ptrElem);

        // get Base Fe, which provides the shape functions for interpolating with the identity operator 
        BaseFE* ptrFe = ptrEsm->GetBaseFE();

        Vector<Double> dist(numElemNodes);
        // loop over element nodes and get quantities
        for (UInt iNodes = 0; iNodes < numElemNodes; ++iNodes) {
          nodeIdx = GetIdxByNodeId(nodeIds[iNodes],0);
          dist[iNodes] = thicknessOnNodes_[nodeIdx];
          //std::cout << dist[iNodes];
        }
        // create 1 element vector for scalar because ApplyOp only takes Vectors
        Vector<Double> d(1);
        // interpolate from nodes to the local point mapped
        this->scalarMappingOperator_->ApplyOp(d,lpm,ptrFe,dist);
        // extract from vector
        val = d[0];
        break;
      }
      default:
        EXCEPTION("CoefFunctionCurvilinearPML::GetScalar(Double...) is used for OutputType: " <<
                  "DISTANCE only.");
    }
  };


  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetVector(Vector<Complex>& val, const LocPointMapped& lpm ) {
    switch(outputType_) {
      case OutputType::DAMP_FACTOR: {
        // here we don't need to interpolate the geometry vectors...
        // get a pointer to the considered element from the lpm
        const Elem* ptrElem = NULL;
        ptrElem = lpm.ptEl;
        // the global element ID
        UInt elemId = ptrElem->GetElemNum();
        // get the nodes connected to the element
        StdVector<UInt> nodeIds;
        grid_->GetElemNodes(nodeIds, elemId);

        // index of a node id
        UInt nodeIdx;
        // index of corresponding node on the interface
        UInt nodeIdxIface;
        // number of nodes in element
        UInt numElemNodes = nodeIds.GetSize();

        // get the ElemShapeMap of the current element from the grid (attention! updated geometry is not set!)
        shared_ptr<ElemShapeMap> ptrEsm = this->grid_->GetElemShapeMap(ptrElem);

        // get Base Fe, which provides the shape functions for interpolating with the identity operator 
        BaseFE* ptrFe = ptrEsm->GetBaseFE();
        
        if (this->dim_ == 2) 
        {
          // Vectors to store geometry data at nodes
          Vector<Double> maxPrincCurv(numElemNodes);
          Vector<Double> dist(numElemNodes);

          // Loop over element nodes and get necessary quantities
          for (UInt iNodes = 0; iNodes < numElemNodes; ++iNodes) 
          {
            nodeIdx = GetIdxByNodeId(nodeIds[iNodes], 0);
            nodeIdxIface = nodeIdx % numSurfNodes_;
            
            maxPrincCurv[iNodes] = ptrNodeGeom_->maxPrincipalCurvatures_[nodeIdxIface];
            dist[iNodes] = thicknessOnNodes_[nodeIdx];
          }

          // Create single-element vectors for interpolation since ApplyOp only takes vectors
          Vector<Double> kmax(1);
          Vector<Double> d(1);

          // Interpolate curvature and distance to the local point
          this->scalarMappingOperator_->ApplyOp(kmax, lpm, ptrFe, maxPrincCurv);
          this->scalarMappingOperator_->ApplyOp(d, lpm, ptrFe, dist);

          // Compute the damping function and its integral
          Double dampFunc = this->dampFunction_->ComputeFactor(d[0], layerThickness_);
          Double intDampFunc = this->dampFunction_->ComputeIntegralFactor(d[0], layerThickness_);

          // Get speed of sound at current lpm
          Double sos;
          this->speedOfSound_->GetScalar(sos, lpm);

          // Compute wave number
          Double K = this->omega_ / sos;

          // Compute helper functions
          Vector<Double> h(2);
          h[0] = dampFunc;
          h[1] = intDampFunc * kmax[0] / (1.0 + d[0] * kmax[0]);

          // the eigenvalues are of the form s_i = (1 + 1i/K * h[i])
          //TODO CO isnt there a minus sign missing in the imaginary part?
          val.Resize(this->dim_);
          for (UInt iDim = 0; iDim < this->dim_; ++iDim) {
            val[iDim] = Complex(1, h[iDim] / K);
          }
        } // dim_ == 2
        else if (this->dim_ == 3) {
          // vectors with extracted geometry data
          Vector<Double> minPrincCurv(numElemNodes);
          Vector<Double> maxPrincCurv(numElemNodes);
          // and on node-to-interface distances
          Vector<Double> dist(numElemNodes);

          // loop over element nodes and get quantities
          for (UInt iNodes = 0; iNodes < numElemNodes; ++iNodes) {
            nodeIdx = GetIdxByNodeId(nodeIds[iNodes],0);
            nodeIdxIface = nodeIdx % numSurfNodes_;
            minPrincCurv[iNodes] = ptrNodeGeom_->minPrincipalCurvatures_[nodeIdxIface];
            maxPrincCurv[iNodes] = ptrNodeGeom_->maxPrincipalCurvatures_[nodeIdxIface];
            dist[iNodes] = thicknessOnNodes_[nodeIdx];
          }

          // create 1-element vectors for scalars because ApplyOp only takes Vectors
          Vector<Double> kmin(1);
          Vector<Double> kmax(1);
          Vector<Double> d(1);

          // interpolate from nodes to the local point mapped
          this->scalarMappingOperator_->ApplyOp(kmin,lpm,ptrFe,minPrincCurv);
          this->scalarMappingOperator_->ApplyOp(kmax,lpm,ptrFe,maxPrincCurv);
          this->scalarMappingOperator_->ApplyOp(d,lpm,ptrFe,dist);

          // get the damping function and its integral at the mapped distance
          Double dampFunc = this->dampFunction_->ComputeFactor(d[0], layerThickness_);
          Double intDampFunc = this->dampFunction_->ComputeIntegralFactor(d[0], layerThickness_);

          // get speed of sound at current lpm
          Double sos;
          this->speedOfSound_->GetScalar(sos,lpm);

          // compute the current wave number
          Double K = this->omega_ / sos;

          // resize output vector
          val.Resize(this->dim_);

          // compute helper functions
          Vector<Double> h = Vector<Double>(3);
          h[0] = dampFunc;
          h[1] = intDampFunc * kmin[0] / (1.0 + d[0] * kmin[0]);
          h[2] = intDampFunc * kmax[0] / (1.0 + d[0] * kmax[0]);

          // the eigenvalues are of the form s_i = (1 + 1i/K * h[i])
          val.Resize(this->dim_);
          for (UInt iDim = 0; iDim < this->dim_; ++iDim) {
            val[iDim] = Complex(1, h[iDim] / K);
          }
        } // dim_ == 3
        break;
      }
      default:
        EXCEPTION("CoefFunctionCurvilinearPML::GetVector(Vector<Complex>...) is used for OutputType: " <<
                  "DAMP_FACTOR only.");
    }
  };


  template<typename T>
  UInt CoefFunctionCurvilinearPML<T>::GetIdxByNodeId(const UInt& nodeId, unsigned int guess) const {
    // here we assume that the interface surface nodes are the numSurfNodes_ lowest entries. 
    // The geometry is computed only there, and we assume the automatically generated layer
    // has the same number of nodes in every iso surface. Thus, we use the modulo to access 
    // the respective interface node geometry for every node in the layer.
    // The initial guess can be provided to speed up the computation, as nodes in elements will be 
    // located at similar indices in the ptrNodeGeom_
    
    UInt idx = 0;
    // in case the node is at the interface, use the find algorithm as nodes may not be ordered here
    if (nodeId < minIds_[1]) {
      int tmpIdx = ptrNodeGeom_->nodeIds_.Find(nodeId, guess, false);
      idx = tmpIdx;
      //std::cout << "idx is:" << idx << "\n" << "id is:" << nodeId << "\n";
      return idx;
    } else {
      // if the nodeId lies within the region of the generated layers, we can speed up the search...
      for (UInt iLayers = 1; iLayers <= numLayers_; ++iLayers) {
        if (nodeId >= minIds_[iLayers] && nodeId <= maxIds_[iLayers]) {
          idx = (nodeId - minIds_[iLayers]) + (numSurfNodes_ * iLayers);
          //std::cout << "idx is:" << idx << "\n" << "id is:" << nodeId << "\n";
          return idx;
        }
      }
    }
    // we should not arrive here...
    EXCEPTION("Node Id could not be determined to lie on an iso surface of the auto generated layer");
  };


  template<typename T>
  void CoefFunctionCurvilinearPML<T>::GetThicknessOnNodes() {
    if (layerGenNode_) {
      // every iso surface in the automatically generated layer has the 
      // same number of nodes. Hence, finding the thickness at nodes is easy...
      // iterate over all surface layers (generated and the interface) and assign
      for (UInt iLayers = 0; iLayers <= numLayers_; ++iLayers) {
        for (UInt iNodes = 0; iNodes < numSurfNodes_; ++iNodes) {
          thicknessOnNodes_.Push_back(elemHeight_ * iLayers);
        }
      }
    } else {
      // room for future implementation that e.g. reads data from a file
      EXCEPTION("CoefFunctionCurvilinearPML::GetThicknessOnNodes currently only implemented with auto layer generation specified.");
    }
  };


  template<typename T>
  void CoefFunctionCurvilinearPML<T>::CreateMappingOperators(){
    // create identity operators for 2D or 3D problems
    if (this->dim_ == 2) {
      this->scalarMappingOperator_.reset(new IdentityOperator<FeH1,2,1,Double>());
      // in the vector case we also need to interpolate vector quantites
      if (this->dimType_ == CoefFunction::VECTOR)
        this->vectorMappingOperator_.reset(new IdentityOperator<FeH1,2,2,Double>());
      // in the tensor case we also need to interpolate tensorial quantites
      if (this->dimType_ == CoefFunction::TENSOR)
        this->tensorMappingOperator_.reset(new IdentityOperator<FeH1,2,4,Complex>());
    }
    else if (this->dim_ == 3) {
      this->scalarMappingOperator_.reset(new IdentityOperator<FeH1,3,1,Double>());
      // in thevector case we also need to interpolate vector quantites
      if (this->dimType_ == CoefFunction::VECTOR)
        this->vectorMappingOperator_.reset(new IdentityOperator<FeH1,3,3,Double>());
      // in the tensor case we also need to interpolate tensorial quantites
      if (this->dimType_ == CoefFunction::TENSOR)
        this->tensorMappingOperator_.reset(new IdentityOperator<FeH1,3,9,Double>());
    }
  };


  template<typename T>
  void CoefFunctionCurvilinearPML<T>::SetDimType() {
    switch (outputType_) {
      case OutputType::TENSOR:
        this->dimType_ = CoefFunction::TENSOR;
        break;
      case OutputType::DETERMINANT:
        this->dimType_ = CoefFunction::SCALAR;
        break;
      case OutputType::DAMP_FACTOR:
        this->dimType_ = CoefFunction::VECTOR;
        break;
      case OutputType::DISTANCE:
        this->dimType_ = CoefFunction::SCALAR;
        break;
      default:
        EXCEPTION("In CoefFunctionCurvilinearPML: outputType '" << outputType_ <<
                        "' is invalid.");
    }
  };


  template<typename T>
  void CoefFunctionCurvilinearPML<T>::ReadDataPML(PtrParamNode pmlDef){
    // read and set coordinate system 
    std::string cSysId;
    pmlDef->GetValue("coordSysId",cSysId,ParamNode::PASS);
    // currently this type wil only be implemented to work in Cartesian coordinates
    // to work with other coordinate systems we need to compute different transformation matrices
    if(cSysId == "" || cSysId == "default")
      cSysId = "default";
    else
      EXCEPTION("Curvilinear PML currently only works in the default Cartesian coordinate system!");
    this->SetCoordinateSystem(domain->GetCoordSystem(cSysId));

    // type of PML damping function
    std::string typeOfPml;
    pmlDef->GetValue("type",typeOfPml);
    this->pmlType_ = DampFunction::DampingTypeEnum.Parse(typeOfPml);
    // set the DampingFunction object to the corresponding type 
    CreateDampFunction();

    // check for damping factor
    Double dampFactor;
    pmlDef->GetValue("dampFactor",dampFactor);
    this->dampFunction_->DampFactor = dampFactor;

    // extract the ID of the PML volume region
    volRegion_ = this->entities_[0]->GetRegion();
  }


  template<typename T>
  void CoefFunctionCurvilinearPML<T>::CheckForInvalidParams(PtrParamNode pmlDef) {
    // check for propRegion, scaling or frequency shift coeff in the xml
    // if these are set for curvilinear PML, ignore and warn
    PtrParamNode propRegionNode = pmlDef->Get("propRegion", ParamNode::PASS);
    if (propRegionNode) {
      std::string propRegionNodeFormul; 
      propRegionNode->GetParent()->GetValue("formulation", propRegionNodeFormul, ParamNode::PASS);
      if (propRegionNodeFormul == "curvilinear")
        WARN("propRegion is currently not implemented for curvilinear PML and will be ignored!");
    }
    PtrParamNode scalingNode = pmlDef->Get("scalingCoef", ParamNode::PASS);
    if (scalingNode) {
      std::string scalingNodeFormul; 
      scalingNode->GetParent()->GetValue("formulation", scalingNodeFormul, ParamNode::PASS);
      if (scalingNodeFormul == "curvilinear")
        WARN("scalingCoef is currently not implemented for curvilinear PML and will be ignored!");
    }
    PtrParamNode shiftNode = pmlDef->Get("frqShiftCoef", ParamNode::PASS);
    if (shiftNode) {
      std::string shiftNodeFormul; 
      shiftNode->GetParent()->GetValue("formulation", shiftNodeFormul, ParamNode::PASS);
      if (shiftNodeFormul == "curvilinear")
        WARN("frqShiftCoef is currently not implemented for curvilinear PML and will be ignored!");
    }
  };


  template<typename T>
  void CoefFunctionCurvilinearPML<T>::CheckForLayerGenerationNode(PtrParamNode pmlDef) {
    // check wether the PML layer has been generated using autoLayerGenerationList and store the parameter node
    PtrParamNode layerGenNode = pmlDef->GetParent()->GetParent()->GetParent()->
                                        GetParent()->GetParent()->Get("domain")->
                                        Get("layerGenerationList", ParamNode::PASS);
    // check if the layerGenerationList has actually generated the PML region
    bool nodeFound = false;
    if (layerGenNode) {
      std::string generatedLayerName;
      RegionIdType generatedLayerId;
      PtrParamNode actLayerGenNode;
      // iterate over the autoLayerGenerationList...
      UInt numGenLayers = layerGenNode->GetChildren().size();
      for (UInt iLayers = 0; iLayers < numGenLayers; iLayers++) {
        PtrParamNode actLayerGenNode = layerGenNode->GetChildren()[iLayers];
        // get name of generated region
        actLayerGenNode->GetValue("name", generatedLayerName);
        // get its ID
        generatedLayerId = grid_->GetRegionId(generatedLayerName);
        // check if generated ID fits to the PML region
        if (volRegion_ == generatedLayerId) {
          layerGenNode_ = actLayerGenNode;

          // assign layer generation parameters for later use
          std::string surfRegionName;
          layerGenNode_->Get("extrusionParameters")->GetValue("numLayers", numLayers_);
          layerGenNode_->Get("extrusionParameters")->GetValue("elemHeight", elemHeight_);
          layerGenNode_->Get("sourceSurfRegion")->GetValue("name", surfRegionName);
          ifSurfRegion_ = grid_->GetRegionId(surfRegionName);

          // set the total layer thickness
          layerThickness_ = numLayers_ * elemHeight_;
          // obtain the number of nodes on the interface region
          numSurfNodes_ = grid_->GetNumNodes(ifSurfRegion_);
          // set bool to succeed search
          nodeFound = true;
          break;
        }
      }
    }
    else if (nodeFound == false || !layerGenNode) {
      EXCEPTION("Specified PML layer not found in layerGenerationList");
    }
  };
  template class CoefFunctionCurvilinearPML<Double>;
  template class CoefFunctionCurvilinearPML<Complex>;
}