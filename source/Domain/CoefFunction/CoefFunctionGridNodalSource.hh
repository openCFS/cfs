// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionGridNodalDefault.hh 
 *       \brief    CoefFunction specialized for the case that the source Grid is the default grid
 *
 *       \date     May 2015
 *       \author   Manfred Kaltenbbacher
 */
//================================================================================================


#ifndef COEFFUNCTIONGRIDNODALSOURCE_HH
#define COEFFUNCTIONGRIDNODALSOURCE_HH

#include "CoefFunctionGridNodal.hh"


namespace CoupledField{

/*! \class CoefFunctionGridNodalDefault
 *  \brief This class represents the specialization of interpolation functionality if the 
 *         source date is on the default grid or a topological equivalent
 *  \tparam DATA_TYPE Can be Complex or Double
 *  @author A. Hueppe
 *  @date 01/2013
 * 
 *  The class is completely based upon two assumptions:
 *   \li The source data is stored on the default Grid or an topological equivalent
 *   \li The results to be mapped are stored corresponding to Nodes
 *  
 *  NOTE: Extensions needed for Results stored by elements
 */

template<typename DATA_TYPE>
class CoefFunctionGridNodalSource : public CoefFunctionGridNodal<DATA_TYPE>{
public:
  
  CoefFunctionGridNodalSource(Domain* ptDomain,
                               PtrParamNode configNode,PtrParamNode curInfo);

  virtual ~CoefFunctionGridNodalSource(){
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

  //! adds entities to the coeffunction
  virtual void AddEntityList(shared_ptr<EntityList> ent);
  

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

  //! computes the optimality condition
  virtual void ComputeOptCondition(Double& optAmp, Double& optPhase);

  //! update the source values (amplitude and phase)
  virtual void UpdateSource(Double& stepLength, bool lineSearch);

  //! computes the L2 norm of error
  virtual void ComputeTikh(Double& funcVal, Double& resSquared,
                           bool adjustAlpha, bool adjustBeta);

  //! computes the L2 norm of error
  virtual void ComputeDiff2Meas( Double& error );

  //! compute square of L2-norm of measured pressure at mic-positions
  virtual void ComputeMeasL2squared( Double& vaL2 );

  //! computes the L2 norm of error
  virtual void SetInverseParam( Double& alpha, Double& beta, Double& qExp,
		                        Double& freq, std::string fileNameMeasdata);


protected:


private:
  ///method searching for elements for given global points
  void GetElemsForPoints(const StdVector<Vector<Double> >  & points, StdVector< const Elem* > & elements, StdVector<LocPoint> & locals);

  //! build associative array for conservative
  void BuildNodeIdxAssoc(shared_ptr<FeSpace> targetSpace);

  //! computes the source values
  void ComputeSourceData(Vector<DATA_TYPE>& actPDEsol);

  //! for conservative mapping, we store just the pair of own
  //! solution vector to coeffunction vector
  StdVector< std::pair<UInt,UInt> > fctSolAssoc_;

  //! nodes of microphone arrays (measured values)
  StdVector<UInt> measNodes_;

  //!Flag indicating if Conservative is ready yet
  bool conservativeReady_;

  ///flag indicating if target region is surface
  bool srcIsSurface_;

  //! conatins all  nodes of the source region
  shared_ptr<EntityList> nodeListSource_;

  //!contains the amplitude of the sources at the nodes
  Vector<Double> sourceAmp_;

  //!contains the incremental amplitude of the sources at the nodes
   Vector<Double> sourceAmpDelta_;

   //!saves amplitude of the sources at the nodes (due to line search)
   Vector<Double> sourceAmpSave_;

  //!contains the phase of the sources at the nodes
  Vector<Double> sourcePhi_;

  //!contains the incremental phase of the sources at the nodes
  Vector<Double> sourcePhiDelta_;

  //!saves phase of the sources at the nodes (due to line search)
  Vector<Double> sourcePhiSave_;

  //! regularization parameter 1
  Double alpha_;

  //! regularization parameter 2
  Double beta_;

  //! regularization parameter 3
  Double qExp_;

  //! defines, if measured data is read at micro-positions!
  bool isDataReadFromFile_;

  //! measurement data read from file
  Vector<DATA_TYPE> readMeasVec_;

  //! stores the measured data
  Vector<DATA_TYPE> measVec_;

  StdVector<bool> isMeasuredNode_;

  std::string typeCoeff_;

  Double meanElemVol_;

  Double freq_;

  Double scalingHesse_;

 };
   

}
#endif
