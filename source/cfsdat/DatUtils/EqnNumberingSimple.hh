// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     EqnNumberingSimple.hh
 *       \brief    Declaration of class to provide a fast, simple equation mapping
 *
 *       \date     Aug 24, 2015
 *       \author   ahueppe
 */
//================================================================================================

#ifndef EQNNUMBERINGSIMPLE_HH_
#define EQNNUMBERINGSIMPLE_HH_

#include "General/defs.hh"
#include "General/Environment.hh"
#include "Utils/ThreadLocalStorage.hh"
#include <memory>
#include "Domain/Results/BaseResults.hh"
#include "cfsdat/DatUtils/Defines.hh"


namespace CFSDat{

//!  This class implements a simple equation mapping for CFS grids.
/*!
  The Equation mapping is based on the assumption that
  almost every entity in the input grid (elements or nodes)
  has an equation number. Therefore we have only minimal overhead
  if we resize the vector of equations to the maximum number of
  entities of that type e.g. numNodes in the grid.

  Linear access time.

  Higher order elements are not supported.
*/
class EqnMapSimple{

public:
   EqnMapSimple(CoupledField::ResultInfo::EntityUnknownType type,
                     CoupledField::Grid* myGrid,
                     CoupledField::UInt eqnPerEnt = 1, bool isOneBased = true);

  virtual ~EqnMapSimple(){
  }

  void AddRegion(CoupledField::RegionIdType region){
    regions_.Push_back(region);
  }

  virtual void Finalize();

  virtual CoupledField::UInt GetNumEquations() const {
    return numEqns_;
  }
  
  virtual CoupledField::UInt GetNumEntities() const {
    return numEqns_ / eqnPerEnt_;
  }
  
  virtual CoupledField::UInt GetNumEqnPerEnt() const {
    return eqnPerEnt_;
  }
  
  CF::ResultInfo::EntityUnknownType GetMapType() const {
    return mapType_;
  }
  
  bool Equals(EqnMapSimple& other) {
    if ((ptGrid_ != other.ptGrid_) || (mapType_ != other.mapType_)
      || (regions_ != other.regions_) || (zeroOne_/eqnPerEnt_) != (other.zeroOne_/other.eqnPerEnt_)) {
      return false;
    }
    return true;
  }

  virtual void GetEquation(CF::StdVector<UInt> & eqns,
                           const UInt globalEntNum,
                           CF::ResultInfo::EntityUnknownType type) const;

  virtual CoupledField::UInt GetEntityIndex(const UInt globalEntNum) const;
  
  virtual bool IsEntityUsed(const UInt globalEntNum) const;
                           
  void GetRegionEquations(CF::StdVector<UInt> & eqns, CF::RegionIdType region) const;

  void GetSubsetEquations(CF::StdVector<UInt> & eqns, CF::StdVector<UInt> & globalEntNumbers) const;
  
  void GetReverseEntityMap(CF::StdVector<UInt>& revMap) const;
  
  
protected:
  //! equation indices
  CF::StdVector<CF::UInt> entityEquations_;
  
  //! stores if entity is used/mapped by this map
  std::vector<bool> usedEntity_;

  //! Array of regions managed by the eqnMap
  CF::StdVector<CF::RegionIdType> regions_;

  //! pointer to grid associated to equations
  CF::Grid* ptGrid_;

  //! stores entry type of eqnMap
  CF::ResultInfo::EntityUnknownType mapType_;

  //! Equations per entity
  UInt eqnPerEnt_;

  //! number of equations
  UInt numEqns_;

  //! maximum number of equations
  UInt maxNumEqns_;

  //! offset switch for zero and one based entity numbers
  UInt zeroOne_;

  //! flag indicating if equation numbering is finalized
  bool isFinalized_;

};

}

#endif /* EQNNUMBERINGSIMPLE_HH_ */
