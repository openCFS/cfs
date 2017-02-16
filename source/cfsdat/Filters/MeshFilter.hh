// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     MeshFilter.hh
 *       \brief    <Description>
 *
 *       \date     Jan 11, 2017
 *       \author   kroppert
 */
//================================================================================================

#pragma once

#include <Filters/BaseMeshFilterType.hh>

namespace CFSDat{


//! Base class for Grid based interpolation- or differentiation-schemes

//! Some procedures require mesh results for
//! input and output. This class serves as a base for algorithms
//! requiring this kind of result representation
class MeshFilter : public BaseMeshFilterType{

public:

  //! Struct for storing various quantities, needed for interpolation
  //! and differentiation. TODO Maybe in the next cleanup-sweep, the PrepareCalculations()
  //! can be somehow unified, because they are quite messy
  struct QuantityStruct{
    CF::StdVector<CF::UInt> tNNum; //node numbers of specific element
    CF::StdVector<CF::UInt> srcEqn; //all coordinates
    CF::Vector<Double> localCoords;
    UInt srcEqnSingle; //only the first coordinate
    UInt trgElemNum;
    UInt srcElemNum;
    Double volume;

    QuantityStruct() : srcEqnSingle(0), trgElemNum(0), srcElemNum(0), volume(.0){
      tNNum.Resize(1);
      srcEqn.Resize(0);
      localCoords.Resize(3);
    }

    bool operator < (const QuantityStruct& str) const
    {
        return (srcEqnSingle < str.srcEqnSingle);
    }
  };

  //! Constructor, which reads the input- and output-result names from the xml-file and stores it
  //! in a string for further checking, if the "connectivity" of filters is correct
  //! Also the source- and target-regions are read from the xml-file
  //! Last step of this constructor is to read the target-mesh and create
  //! a new GridCFS-object
  MeshFilter(UInt numWorkers, CF::PtrParamNode config, str1::shared_ptr<ResultManager> resMan);

  virtual ~MeshFilter(){
    delete trgGrid_;
  }

  virtual bool Run() = 0;

  virtual void FinishInit();



protected:

  virtual void PrepareCalculation() = 0;

  virtual ResultIdList SetUpstreamResults()=0;

  virtual void AdaptFilterResults()=0;

  //TODO for now copy paste code from inputfilter... not very nice
  void CreateDummyCfsParamNode();

  /// XML node for creating a CFS input object
  PtrParamNode dummyXMLNode_;

  ///Set of input/src regions
  std::set<std::string> srcRegions_;

  ///Set of target regions
  std::set<std::string> trgRegions_;

  //Siminput pointer to target mesh
  str1::shared_ptr<CoupledField::SimInput> trgMeshInp_;

  ///Pointer to new grid object
  Grid* trgGrid_;


  //***************************************************************************
  //********************** Interpolation Section ******************************
  //***************************************************************************

  //! Interpolation from nodes to element results by averaging over element-nodes

  //! \param returnVec (out) vector containing the interpolated results
  //! \param resId (in) result-id of the current filter results
  //! \param inVec (in) Vector containing the node-results
  //! \param interpolData (in) information of elements of target-regions
  void Node2Cell(Vector<Double>& returnVec,
                const boost::uuids::uuid& resId,
                const Vector<Double>& inVec,
                const std::vector<QuantityStruct>& interpolData);



  //! Interpolation from cell-center (element-results) to node-results by
  //! averaging across elements

  //! \param returnVec (out) vector containing the interpolated results
  //! \param resId (in) result-id of the current filter results
  //! \param inVec (in) Vector containing the node-results
  //! \param interpolData (in) information of elements of target-regions
  //! \param nodeNeighbours (in) for every source point (element-result)
  //!                            the number of nodes is stored (e.g. for tet 4
  //!                            or hex 8)
  void Cell2Node(Vector<Double>& returnVec,
                  const boost::uuids::uuid& resId,
                  const Vector<Double>& inVec,
                  const std::vector<QuantityStruct>& interpolData,
                  const StdVector<UInt>& nodeNeighbours);

  //! Nearest neighbour interpolation

  //! \param returnVec (out) vector containing the interpolated results
  //! \param vec (in) vector containing the interpolated result for every
  //!                 target-node, can be 1D for scalar values or 2D/3D
  //!                 for vector values, no return-parameter, only needed
  //!                 for intern-calculation
  //! \param downMap (in) pointer to equation numbers of target-mesh
  //! \param srcCoords (in) coordinate matrix of source-mesh
  //! \param trgCoords (in) coordinate matrix of target-mesh
  //! \param interpolationData (in) information of elements of target-regions
  //! \param grid (in) pointer to the target-grid
  //! \param numNeighbors (in) number of points, which shalle be considered for
  //!                          the nearest neighbour search
  //! \param p (in) Interpolation exponent
  void NearestNeighbourInterpolation(Vector<Double>& returnVec,
                                    const CF::StdVector< CF::Vector<Double> >&  scatteredData,
                                    CF::Vector<Double>& vec,
                                    const str1::shared_ptr<EqnMapSimple>& downMap,
                                    const CF::StdVector< CF::Vector<double> >& srcCoords,
                                    const CF::StdVector< CF::Vector<double> >& trgCoords,
                                    const std::vector<QuantityStruct>& interpolationData,
                                    Grid* grid,
                                    const UInt& dim,
                                    const UInt& numNeighbors,
                                    const UInt& p);


  //! Perform interpolation, based on a local radial basis functions (RBF)-approach
  //! Performs good, if the source mesh is not significantly bigger than the target
  //! mesh, otherwise, oscillations can occur in the interpolated results,
  //! if this happens, the FEBasedInterpolator is better suited

  //! \param returnVec (out) vector containing the interpolated results
  //! \param vec (in) vector containing the interpolated result for every
  //!                 target-node, can be 1D for scalar values or 2D/3D
  //!                 for vector values, no return-parameter, only needed
  //!                 for intern-calculation
  //! \param downMap (in) pointer to equation numbers of target-mesh
  //! \param srcCoords (in) coordinate matrix of source-mesh
  //! \param trgCoords (in) coordinate matrix of target-mesh
  //! \param interpolationData (in) information of elements of target-regions
  //! \param grid (in) pointer to the target-grid
  //! \param dim (in) Dimension of the data values
  //! \param p (in) Interpolation exponent, also controls the "locality" of the
  //!               RBF Ansatz-functions
  void RBFInterpolation(Vector<Double>& returnVec,
                        const CF::StdVector< CF::Vector<Double> >&  scatteredData,
                        CF::Vector<Double>& vec,
                        const str1::shared_ptr<EqnMapSimple>& downMap,
                        const CF::StdVector< CF::Vector<double> >& srcCoords,
                        const CF::StdVector< CF::Vector<double> >& trgCoords,
                        const std::vector<QuantityStruct>& interpolationData,
                        Grid* grid,
                        const UInt& dim,
                        const UInt& p,
                        const bool noSlip);


  //! Build the local RBF interpolation matrix, based on the nearest neighbour source points,
  //! build the local interpolation-value vector and solve the system ALoc*c=vector for
  //! the local RBF coefficients

  //! \param coefVec (out) vector containing the RBF-coefficients for the specific point
  //! \param globPoint (in) coordinate of the point, where the interpolation is performed
  //! \param neighbors (in) coordinates of the nearest neighbors of globPoint
  //! \param l2Distances (in) distances of the nearest neighbors of globPoint
  //! \param vectors (in) containing the values of the nearest neighbors
  //! \param numNN (in) number of nearest neighbor points
  //! \alpha (in) shape parameter for local radial basis functions, controls
  //! the locality, the bigger alpha is, the more "global" the basis function becomes
  //! \param inDim (in) Dimension of input-data
  void CalcLocRBFCoefs(CF::Matrix<Double>& coefVec,
                       const CF::Vector<Double>& globPoint,
                       const CF::StdVector< Vector<Double> >& neighbors,
                       const CF::StdVector< Double >& l2Distances,
                       const CF::StdVector< Vector<Double> >& vectors,
                       const UInt& numNN,
                       const Double& alpha,
                       const UInt& inDim);




  //***************************************************************************
  //***************** Spatial Derivative Section ******************************
  //***************************************************************************
  // Here, only the local methods are defined, because the spatial-
  // derivative filters have specialized methods to combine the local results
  // into global ones (e.g. different handling of physical (AeroacousticBase) and
  // generic (spatial) values

  //! Calculate the local spatial-derivatives and combine it in such
  //! a way that we obtain the local curl

  //! \param vec (out) vector containing the result for the specific point
  //!                  will be inserted into the result-vector
  //! \param globPoint (in) coordinate of the point, where the curl shall be calculated
  //! \param neighbors (in) coordinates of the nearest neighbors of globPoint
  //! \param l2Distances (in) distances of the nearest neighbors of globPoint
  //! \param vectors (in) containing the values of the nearest neighbors
  //! \param numNN (in) number of nearest neighbor points
  //! \alpha (in) shape parameter for local radial basis functions, controls
  //! the locality, the bigger alpha is, the more "global" the basis function becomes
  //! \param dim (in) dimension of the data (scalar=0, 2d=1, 3d=2)
  //! \param grid (in) pointer to the target-grid
  void CalcLocCurl(CF::Matrix<Double>& vec,
                                const CF::Vector<Double>& globPoint,
                                const CF::StdVector< Vector<Double> >& neighbors,
                                const CF::StdVector< Double >& l2Distances,
                                const CF::StdVector< Vector<Double> >& vectors,
                                const UInt numNN,
                                const Double alpha,
                                const UInt dim,
                                Grid* grid);

  //! Calculate the local spatial-derivatives and combine it in such
  //! a way that we obtain the local gradient

  //! \param vec (out) vector containing the result for the specific point
  //!                  will be inserted into the result-vector
  //! \param globPoint (in) coordinate of the point, where the curl shall be calculated
  //! \param neighbors (in) coordinates of the nearest neighbors of globPoint
  //! \param l2Distances (in) distances of the nearest neighbors of globPoint
  //! \param vectors (in) containing the values of the nearest neighbors
  //! \param numNN (in) number of nearest neighbor points
  //! \alpha (in) shape parameter for local radial basis functions, controls
  //! the locality, the bigger alpha is, the more "global" the basis function becomes
  //! \param dim (in) dimension of the data (scalar=0, 2d=1, 3d=2)
  //! \param grid (in) pointer to the target-grid
  void CalcLocGradient(CF::Matrix<Double>& vec,
                                const CF::Vector<Double>& globPoint,
                                const CF::StdVector< Vector<Double> >& neighbors,
                                const CF::StdVector< Double >& l2Distances,
                                const CF::StdVector< Vector<Double> >& vectors,
                                const UInt numNN,
                                const Double alpha,
                                const UInt dim,
                                Grid* grid);


//! Calculate the local spatial-derivatives and combine it in such
//! a way that we obtain the local divergence

//! \param vec (out) vector containing the result for the specific point
//!                  will be inserted into the result-vector
//! \param globPoint (in) coordinate of the point, where the curl shall be calculated
//! \param neighbors (in) coordinates of the nearest neighbors of globPoint
//! \param l2Distances (in) distances of the nearest neighbors of globPoint
//! \param vectors (in) containing the values of the nearest neighbors
//! \param numNN (in) number of nearest neighbor points
//! \alpha (in) shape parameter for local radial basis functions, controls
//! the locality, the bigger alpha is, the more "global" the basis function becomes
//! \param dim (in) dimension of the data (scalar=0, 2d=1, 3d=2)
//! \param grid (in) pointer to the target-grid
void CalcLocDivergence(CF::Matrix<Double>& vec,
                              const CF::Vector<Double>& globPoint,
                              const CF::StdVector< Vector<Double> >& neighbors,
                              const CF::StdVector< Double >& l2Distances,
                              const CF::StdVector< Vector<Double> >& vectors,
                              const UInt numNN,
                              const Double alpha,
                              const UInt dim,
                              Grid* grid);


//! Bring the input-results and the input-coordinate into the
//! appropriate form a nearest neighbour search with CGAL.
//! The problem is that the input-results, provided by result-manager, is
//! only a vector containing the results in the form (u1,v1,w1,u2,v2,w2,...)
//! and for CGAL we need something like ([u1,v1,w1],[u2,v2,w2],...)

//! \param scatteredData (out) Vector containing the correct ordered input-data
//! \param vec (in) vector containing the result for every
//!                 target-node, can be 1D for scalar values or 2D/3D
//!                 for vector values, no return-parameter, only needed
//!                 for intern-calculation
//! \param inVec (in) vector containing the input-results, provided by the
//!                   result-manager
//! \coords (in) locations, where the values of inVec are defined
template<class T>
void FillScatteredDataVec(CF::StdVector< CF::Vector<Double> >& scatteredData,
                          T& vec,
                          const Vector<Double>& inVec,
                          const CF::StdVector< CF::Vector<double> >& coords,
                          UInt& inDim){


  if(inVec.GetSize() == coords.GetSize()){
    //this is the case of scalar scattered values
    inDim=0;
    vec.Resize(1);
    for(UInt i = 0; i < coords.GetSize(); ++i){
      scatteredData[i].Resize(1);
      scatteredData[i][0] = inVec[i];
     }
    }else{
          if(inVec.GetSize() == 2 * coords.GetSize()){
            //case of two dimensional vector
            inDim=1;
            vec.Resize(2);
            for(UInt i = 0; i < coords.GetSize() ; ++i){
            scatteredData[i].Resize(2);
            scatteredData[i][0] = inVec[2 * i]; // x-component
            scatteredData[i][1] = inVec[2 * i + 1]; // y-component
            }
          }else{
              if(inVec.GetSize() == 3 * coords.GetSize()){
                // case of three dimensional vector
                inDim=2;
                vec.Resize(3);
                for(UInt i = 0; i < coords.GetSize(); ++i){
                scatteredData[i].Resize(3);
                scatteredData[i][0] = inVec[3 * i]; // x-component
                scatteredData[i][1] = inVec[3 * i + 1]; // y-component
                scatteredData[i][2] = inVec[3 * i + 2]; // z-component
                }
              }else{
                std::cout<<"\t\t Size of input-vector:"<<inVec.GetSize()<<std::endl;
                std::cout<<"\t\t Number of source-points:"<<coords.GetSize()<<std::endl;
                EXCEPTION("Incorrect Input Data!")
              }
          }
    }
}


};
}
