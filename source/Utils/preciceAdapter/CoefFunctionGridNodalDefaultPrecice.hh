// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*- 
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     CoefFunctionGridNodalDefaultPrecice.hh 
 *       \brief    CoefFunction specialized for the case that the source Grid is the default grid
 *                 and when using precice coupling adapter
 *
 *       \date     Feb. 18, 2025
 *       \author   Klaus Roppert
 */
//================================================================================================


#ifndef COEFFUNCTIONGRIDNODEALDEFAULTPRECICE_HH
#define COEFFUNCTIONGRIDNODEALDEFAULTPRECICE_HH

#include "def_use_precice.hh"

#include "Domain/CoefFunction/CoefFunctionGridNodal.hh"

namespace CoupledField{

  class IPreciceAdapter;

template<class DATA_TYPE>
class CoefFunctionGridNodalDefaultPrecice : public CoefFunctionGridNodal<DATA_TYPE>{
public:
  
  CoefFunctionGridNodalDefaultPrecice(Domain* ptDomain,
                                      PtrParamNode configNode,
                                      PtrParamNode curInfo,
                                      shared_ptr<RegionList> regions,
                                      ResultInfo::EntryType type);

  virtual ~CoefFunctionGridNodalDefaultPrecice(){
  };

  virtual string GetName() const { return "CoefFunctionGridNodalDefaultPrecice"; }


  // ========================
  //  ACCESS METHODS
  // ========================
  //@{ \name Access Methods

  ///Getting a vector Value at the requested lpm, Invalid call for Conservative interpolation
  virtual void GetVector(Vector<DATA_TYPE>& CoefMat,
                         const LocPointMapped& lpm );

  ///Getting a scalar value at the requested lpm, INvalid call for Conservative interpolation
  virtual void GetScalar(DATA_TYPE& CoefMat,
                         const LocPointMapped& lpm );
  //@}

  
protected:
  
  virtual void DetermineResult(std::string inputID,UInt seqStep) override;

  //! sets the source region and destination regions
  virtual void SetRegions(shared_ptr<RegionList> regions);

private:
  IPreciceAdapter* preciceAdapter_;
};
   

}
#endif
