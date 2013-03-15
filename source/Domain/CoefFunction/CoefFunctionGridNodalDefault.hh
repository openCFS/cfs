// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionGridNodalDefault.hh 
 *       \brief    CoefFunction specialized for the case that the source Grid is the default grid
 *
 *       \date     Jan. 20, 2013
 *       \author   Andreas Hueppe
 */
//================================================================================================


#ifndef COEFFUNCTIONGRIDNODALDEFAULT_HH
#define COEFFUNCTIONGRIDNODALDEFAULT_HH

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
class CoefFunctionGridNodalDefault : public CoefFunctionGridNodal<DATA_TYPE>{
public:
  
  CoefFunctionGridNodalDefault(PtrParamNode configNode);

  virtual ~CoefFunctionGridNodalDefault(){
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

protected:

private:
  ///method searching for elements for given global points
  void GetElemsForPoints(const StdVector<Vector<Double> >  & points, StdVector< const Elem* > & elements, StdVector<LocPoint> & locals);

};
   

}
#endif
