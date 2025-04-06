// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionGridElemDefault.hh 
 *       \brief    CoefFunction specialized for the case that the source Grid is the default grid
 *
 *       \date     Sep. 29, 2024
 *       \author   Dominik Mayrhofer
 */
//================================================================================================


#ifndef COEFFUNCTIONGRIDELEMDEFAULT_HH
#define COEFFUNCTIONGRIDELEMDEFAULT_HH

#include "CoefFunctionGridElem.hh"

namespace CoupledField{

/*! \class CoefFunctionGridNodalDefault
 *  \brief This class represents the specialization of interpolation functionality if the 
 *         source date is on the default grid or a topological equivalent
 *  \tparam DATA_TYPE Can be Complex or Double
 *  @author Dominik Mayrhofer
 *     @date 09/2024
 * 
 *  The class is completely based upon two assumptions:
 *   \li The source data is stored on the default Grid or an topological equivalent
 *   \li The results to be mapped are stored corresponding to elements
 */

template<class DATA_TYPE>
class CoefFunctionGridElemDefault : public CoefFunctionGridElem<DATA_TYPE>{
public:
  
  CoefFunctionGridElemDefault(Domain* ptDomain,
                               PtrParamNode configNode,PtrParamNode curInfo, shared_ptr<RegionList> regions,ResultInfo::EntryType type);

  virtual ~CoefFunctionGridElemDefault(){
  };

  virtual string GetName() const { return "CoefFunctionGridElemDefault"; }


  // ========================
  //  ACCESS METHODS
  // ========================
  //@{ \name Access Methods

  ///Getting a Tensor Value at the requested lpm, Invalid call for Conservative interpolation
  virtual void GetTensor(Matrix<DATA_TYPE>& CoefMat,
                         const LocPointMapped& lpm );

  ///Getting a vector Value at the requested lpm, Invalid call for Conservative interpolation
  virtual void GetVector(Vector<DATA_TYPE>& CoefMat,
                         const LocPointMapped& lpm );

  ///Getting a scalar value at the requested lpm, INvalid call for Conservative interpolation
  virtual void GetScalar(DATA_TYPE& CoefMat,
                         const LocPointMapped& lpm );
  //@}

  // COLLECTION ACCESS
  virtual void GetScalarValuesAtPoints( const StdVector<Vector<Double> >  & points,
                                             StdVector<DATA_TYPE >  & vals);

  virtual void GetVectorValuesAtPoints( const StdVector<Vector<Double> >  & points,
                                             StdVector<Vector<DATA_TYPE> >  & vals);

  virtual void GetTensorValuesAtPoints( const StdVector<Vector<Double> >  & points,
                                            StdVector<Matrix<DATA_TYPE> >  & vals);

  virtual void MapConservative( shared_ptr<FeSpace> targetSpace,
                                    Vector<DATA_TYPE>& feFncVec);

  //! Determine if coefFunction has conservative mapping
  virtual bool IsConservative(){
    return (this->curInterpType_==CoefFunctionGrid::CONSERVATIVE);
  }
  
protected:
  
  //! sets the source region and destination regions
  virtual void SetRegions(shared_ptr<RegionList> regions);

private:
  ///method searching for elements for given global points
  void GetElemsForPoints(const StdVector<Vector<Double> >  & points, StdVector< const Elem* > & elements, StdVector<LocPoint> & locals);

  //! build associative arrays for conservative mapping
  void BuildNodeIdxAssoc(shared_ptr<FeSpace> targetSpace);
  
  //! stores if a value from the solution vector has to be copied to target space vector
  std::vector<bool> copyValueIndex_;
  
  //! stores if to which index of the target space vector a value from 
  // the solution vector has to copied to
  StdVector<UInt> valueTargetIndex_;

  //!Flag indicating if Conservative is ready yet
  bool conservativeReady_;

  ///flag indicating if target region is surface
  bool srcIsSurface_;
};
   

}
#endif
