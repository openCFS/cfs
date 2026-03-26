// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionGridNodalInterp.hh
 *       \brief    The coefficient function for external grids which need interpolation
 *
 *       \date     Jan. 20, 2013
 *       \author   Andreas Hueppe
 */
//================================================================================================

#ifndef COEFFUNCTIONGRIDNODALINTERP_HH
#define COEFFUNCTIONGRIDNODALINTERP_HH

#include "CoefFunctionGridNodal.hh"
#include <boost/random/random_device.hpp>
#include <boost/random/mersenne_twister.hpp>
#include <boost/random/uniform_int_distribution.hpp>

#include <boost/smart_ptr/enable_shared_from_this.hpp>
#include <set>
#include <map>

namespace CoupledField{
/*! \class CoefFunctionGridNodalInterp
 *  \brief This class represents the specialization of interpolation functionality if the
 *         source date is on an external grid and needs interpolation in space
 *  \tparam DATA_TYPE Can be Complex or Double
 *  @author A. Hueppe
 *  @date 01/2013
 *
 *         This class implements the interpolation of external data to the default grid.
 *         It can operate in two modes STANDARD or CONSERVATIVE
 *         The conservative case is used e.g. for the interpolation of aeroacoustic sources
 *         it directly sets the FeFunction vector based upon the supplied data
 *         The standard case makes use of another feFunction defined on the default grid.
 *         This function is just used for interpolation and makes it possible to interpolate using
 *         higher order even though the calculation is done in a lower order setting. This makes
 *         sense because it can always happen that we integrate higher order and the evaluation of
 *         values at the integration points is crucial wrt accuracy.
 *         TODO:
 *          - Tested with 2D data which works fine 3D tests are still to come
 *          - Parametrization within the xml file
 *
 */
template<typename DATA_TYPE>
class CoefFunctionGridNodalInterp : public CoefFunctionGridNodal<DATA_TYPE>
                                   ,public enable_shared_from_this<CoefFunctionGridNodalInterp<DATA_TYPE> >{
public:

  CoefFunctionGridNodalInterp(Domain* ptDomain,
                              PtrParamNode configNode,PtrParamNode curInfo, shared_ptr<RegionList> regions);

  virtual ~CoefFunctionGridNodalInterp() { };

  virtual string GetName() const { return "CoefFunctionGridNodalInterp"; }

  // ========================
  //  ACCESS METHODS
  // ========================
  //@{ \name Access Methods

  ///Getting a Tensor Value at the requested lpm, INvalid call for Conservative interpolation
  virtual void GetTensor(Matrix<DATA_TYPE>& CoefMat,
                         const LocPointMapped& lpm );

  ///Getting a vector Value at the requested lpm, Invalid call for Conservative interpolation
  virtual void GetVector(Vector<DATA_TYPE>& CoefMat,
                         const LocPointMapped& lpm );

  ///Getting a scalar value at the requested lpm, INvalid call for Conservative interpolation
  virtual void GetScalar(DATA_TYPE& CoefMat,
                         const LocPointMapped& lpm );
  //@}

  //! Prepares the interpolation, former done in GetVector method but there we had problems with OMP
  virtual void PrepareInterpolation();


  //@{ \name Shortcuts

  //! Map Conservative to FeFunction Vector
  virtual void MapConservative( shared_ptr<FeSpace> targetSpace,
                                    Vector<DATA_TYPE>& feFncVec);


  //! Determine if coefFunction has conservative mapping
  virtual bool IsConservative(){
    if(this->curInterpType_ == CoefFunctionGrid::CONSERVATIVE)
      return true;
    else
      return false;
  }

  //! Give Values at global coordinate locations
  virtual void GetVectorValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                        StdVector< Vector<DATA_TYPE> >& values, 
                                        Grid* ptGrid,
                                        const StdVector<shared_ptr<EntityList> >& srcEntities =
                                        StdVector<shared_ptr<EntityList> >(),
                                        bool updatedGeo = false);

  //! Give Values at global coordinate locations
  virtual void GetScalarValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                        StdVector< DATA_TYPE >& values, 
                                        Grid* ptGrid,
                                        const StdVector<shared_ptr<EntityList> >& srcEntities =
                                        StdVector<shared_ptr<EntityList> >()  );

  //@}

  //! \copydoc CoefFunction::SetConservative
  virtual void SetConservative(bool value){
    CoefFunctionGrid::SetConservative(value);
    globalTol_ = 0.0;
    this->myConfigNode_->GetValue("globalTolerance", globalTol_, ParamNode::PASS);
    localTol_ = 1e-3;
    this->myConfigNode_->GetValue("localTolerance", localTol_, ParamNode::PASS);
  }

protected:
  
  //! sets the destination regions
  virtual void SetDestRegions(shared_ptr<RegionList> regions);
  
private:

  // =====================================
  // Functions and Variables for conservative interpolation
  // =====================================
  ///Compute datastructures for conservative interpolation
  void MapElemNodesConservative();

  //! print unmaped nodes to CSV file for paraview read in
  virtual void PrintNodesToCSV(const StdVector<const Elem*>& foundElements,
                                  const StdVector< Vector<Double> >& nodeGlobCoords);

  //@{ \name Conservative interpolation datastructures
  ///flag to determine if everything is ready to go
  bool consInterpReady_;

  ///save the association targetElement to source nodes
  std::map< UInt , std::vector <UInt> > elemNodeAssoc_;

  ///save the local coordinates of mapped points. We will use this only temporarily
  std::map<UInt, LocPoint > localCoordsNodeAssoc_;

  ///Matrix storing the weight for conservative interpolation
  shared_ptr<BaseMatrix> consInterpMat_;
  
  //! global tolerance (size in meters of bounding box of nodes)
  Double globalTol_;
  
  //! local tolerance (interval of local coordinates outside of element)
  Double localTol_;
  
  //! For 3D => 2D interpolation use x-y-plane at this z-coordinate
  Double xyPlaneAtZ_;
  
  //! Tolerance for z-coordinate of x-y-plane
  Double zTol_;
  //@}

  // =====================================
  // Functions and Variables for standard interpolation
  // =====================================
  //@{ \name Standard interpolation functions and datastructures
  ///Compute (interpolated) solutions at point do not distiguish between structures
  void ComputeSolutionAtPoint(Vector<DATA_TYPE>& solAtPoint,
                              const LocPointMapped& lpm){

  }

  ///we create the Interpolation space and set its polynomial oder
  void PrepareForStdInterp(PtrParamNode motherNode);

  ///flag indicating if the interpolation space is finalized
  bool stdInterpReady_;
  //@}

  ///pointer to FeFunction defined on the target grid
  shared_ptr<FeFunction<DATA_TYPE> > interpolFunction_;

  ///pointer to own polyIdNode used for configuration of feFunct_, allows higher order interpolation
  PtrParamNode polyNode_;

  ///IntegId node just for compatibility and we set it always to lowest order
  PtrParamNode integNode_;

  ///Set storing the destination region IDs
  std::set<RegionIdType> destRegions_;

  //================================================================================
  // Interpolation Information
  //================================================================================
private:
  void ReadXMLNode(PtrParamNode configNode);

  //! names of all destination regions separated by commas
  std::string destRegionNames_;

  //! Pointer to destination grid for interpolation (usually "default" grid)
  Grid * destGrid_;

  //! cache global coords in case of standard interpolation of transient data
  StdVector<LocPoint> localCoords_;

  // cache found elements in case of standard interpolation of transient data
  StdVector< const Elem* > foundElements_;
  
  // file name to store not interpoleted nodes
  std::string notInterpolatedNodesFileName_;

  //! count the given index in the solvec for verbose sum
  std::vector<bool> countSolvecIndex_;
  
  //! add the given index to this dimension
  StdVector<UInt> solVecIndexDim_;
  
};



}
#endif
