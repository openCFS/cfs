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
 *          - User supplied tolerances are ignored at the moment
 *
 */
template<typename DATA_TYPE>
class CoefFunctionGridNodalInterp : public CoefFunctionGridNodal<DATA_TYPE>
                                   ,public boost::enable_shared_from_this<CoefFunctionGridNodalInterp<DATA_TYPE> >{
public:

  CoefFunctionGridNodalInterp(Domain* ptDomain,
                              PtrParamNode configNode,PtrParamNode curInfo);

  virtual ~CoefFunctionGridNodalInterp(){
  };

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

  //! Dump coefficient function to string
  std::string ToString() const{

    return "";
  };

  virtual void AddEntityList(shared_ptr<EntityList> ents);

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
                                        Grid* ptGrid );

  //! Give Values at global coordinate locations
  virtual void GetScalarValuesAtCoords( const StdVector<Vector<Double> >& globCoord,
                                        StdVector< DATA_TYPE >& values, 
                                        Grid* ptGrid );

  //@}

  //! \copydoc CoefFunction::SetConservative
  virtual void SetConservative(bool value){
    CoefFunctionGrid::SetConservative(value);
    globalTol_ = 0.0;
    this->myConfigNode_->GetValue("globalTolerance", globalTol_, ParamNode::PASS);
    localTol_ = 1e-3;
    this->myConfigNode_->GetValue("localTolerance", localTol_, ParamNode::PASS);
  }


private:

  // =====================================
  // Functions and Variables for conservative interpolation
  // =====================================
  ///Compute datastructures for conservative interpolation
  void MapElemNodesConservative();

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



};



}
#endif
