// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionGridElementDefault.hh
 *       \brief    CoefFunction specialized for the case that the source Grid is the default grid
 *
 *                 NOTE: Still quite a dirty hack, since it's nearly an identical copy
 *                 of the nodal-version.
 *
 *       \date     Apr. 10, 2018
 *       \author   Klaus Roppert
 */
//================================================================================================


#ifndef COEFFUNCTIONGRIDELEMENTDEFAULT_HH
#define COEFFUNCTIONGRIDELEMENTDEFAULT_HH

#include "CoefFunctionGridElement.hh"


namespace CoupledField{

template<typename DATA_TYPE>
class CoefFunctionGridElementDefault : public CoefFunctionGridElement<DATA_TYPE>{
public:
  
  CoefFunctionGridElementDefault(Domain* ptDomain,
                               PtrParamNode configNode,PtrParamNode curInfo, shared_ptr<RegionList> regions);

  virtual ~CoefFunctionGridElementDefault(){
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
  virtual std::string ToString() const{
    return "";
  };

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
