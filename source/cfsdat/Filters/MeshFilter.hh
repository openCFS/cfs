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

#include <limits>

namespace CFSDat{

//! Base class for Grid based interpolation- or differentiation-schemes

//! Some procedures require mesh results for
//! input and output. This class serves as a base for algorithms
//! requiring this kind of result representation
class MeshFilter : public BaseMeshFilterType{

public:


  //! Entity number for unused indices in in- and output results
  static const CF::UInt UnusedEntityNumber = UINT_MAX;

  // TODO this struct is obsolete, it's only used by Cell2Node and Node2Cell
  // we should change those two to the new conventions
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
  MeshFilter(UInt numWorkers, CF::PtrParamNode config, shared_ptr<ResultManager> resMan);

  virtual ~MeshFilter(){
    delete trgGrid_;
  }


protected:

  //! Count the number of entities, which are used for the interpolation
  CF::UInt CountUsedEntities(const StdVector<CF::UInt>& entities);


  //! Collects global entity numbers (nodeNums or elemNums), according to the specified regions of the grid.
  //! If an equation-number is not used for the calculation, the according entry in entities is UInt_MAX
   void GetUsedMappedEntities(const shared_ptr<EqnMapSimple>& map,
                             StdVector<CF::UInt>& entities,
                             const std::set<std::string>& regions,
                             Grid* grid);


  //TODO for now copy paste code from inputfilter... not very nice
  void CreateDummyCfsParamNode();

  /// XML node for creating a CFS input object
  PtrParamNode dummyXMLNode_;

  ///Set of input/src regions
  std::set<std::string> srcRegions_;

  ///Set of target regions
  std::set<std::string> trgRegions_;

  //Siminput pointer to target mesh
  shared_ptr<CoupledField::SimInput> trgMeshInp_;

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
  template<typename T>
  void Node2Cell(Vector<T>& returnVec,
                const boost::uuids::uuid& resId,
                const Vector<T>& inVec,
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

  //! Perform interpolation, based on a local radial basis functions (RBF)-approach
  //! Performs good, if the source mesh is not significantly bigger than the target
  //! mesh, otherwise, oscillations can occur in the interpolated results,
  //! if this happens, the FEBasedInterpolator is better suited
  void RBFInterpolation(Vector<Double>& returnVec,
                        const Vector<Double>& inVec,
                        const UInt& numEquPerEnt,
                        const StdVector<CF::UInt>& targetSource,
                        const StdVector<CF::UInt>& targetSourceIndex,
                        const StdVector< CF::Matrix<Double> >& targetRBFInv,
                        const StdVector<CF::Double>& targetSourceFactor,
                        const StdVector<CF::Double>& targetSourceFactor2,
                        const UInt& maxNumTrgEntities);


  void CalcLocRBFInv(CF::Matrix<Double>& matr,
                    const StdVector< CF::Vector<CF::Double> >& neighbors,
                    const Double& alpha,
                    const UInt numNN,
                    Grid* grid);

  // Light- NearestNeighbour interpolation for the use in AeroSourceTerms
  void NearestNeighbourLight(Vector<Double>& returnVec,
                      const Vector<Double>& inVec,
                      const UInt& numEquPerEnt,
                      const StdVector< StdVector<CF::UInt> >& sourceM,
                      const StdVector< CF::Matrix<CF::Double> >& targetSourceFactor,
                      const UInt& maxNumTrgEntities);


  //***************************************************************************
  //***************** Spatial Derivative Section ******************************
  //***************************************************************************

  bool CalcLocCurl(CF::Matrix<Double>& derivCoefVec,
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
                      const bool logEps);


  bool CalcLocGradient(CF::Matrix<Double>& derivCoefVec,
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
                      const bool logEps);


  bool CalcLocDivergence(CF::Matrix<Double>& derivCoefVec,
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
                      const bool logEps);


  void CalcCurl(Vector<Double>& returnVec,
                const Vector<Double>& inVec,
                const UInt& numEquPerEnt,
                const StdVector< StdVector<CF::UInt> >& sourceM,
                const StdVector< CF::Matrix<CF::Double> >& targetSourceFactor,
                const UInt& maxNumTrgEntities,
                const UInt& gridDim);

  void CalcDivergence(Vector<Double>& returnVec,
                      const Vector<Double>& inVec,
                      const UInt& numEquPerEnt,
                      const StdVector< StdVector<CF::UInt> >& sourceM,
                      const StdVector< CF::Matrix<CF::Double> >& targetSourceFactor,
                      const UInt& maxNumTrgEntities);

  void CalcTensorDivergence(Vector<Double>& returnVec,
                      const Vector<Double>& inVec,
                      const UInt& numEquPerEnt,
                      const StdVector< StdVector<CF::UInt> >& sourceM,
                      const StdVector< CF::Matrix<CF::Double> >& targetSourceFactor,
                      const UInt& maxNumTrgEntities);

  int Index2Voigt(const UInt& dx1,
                  const UInt& dx2,
                  const UInt& dim);

  void CalcGradient(Vector<Double>& returnVec,
                      const Vector<Double>& inVec,
                      const UInt& numEquPerEnt,
                      const StdVector< StdVector<CF::UInt> >& sourceM,
                      const StdVector< CF::Matrix<CF::Double> >& targetSourceFactor,
                      const UInt& maxNumTrgEntities,
                      const UInt& tDim);
};
}
