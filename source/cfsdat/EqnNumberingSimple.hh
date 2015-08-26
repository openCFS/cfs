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
#include "boost/shared_ptr.hpp"

namespace CFSDat{

namespace CF = CoupledField;

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
   EqnMapSimple(CF::ResultInfo::EntityUnknownType type,
                     boost::shared_ptr<CF::Grid> myGrid,
                     CF::UInt eqnPerEnt = 1, bool isOneBased = true);

  virtual ~EqnMapSimple(){
    delete entityEquations_;
  }

  void AddRegion(CF::RegionIdType region){
    regions_.Push_back(region);
  }

  virtual void Finalize();

  virtual void GetEquation(CF::StdVector<CF::UInt> & eqns,
                           const CF::UInt globalEntNum,
                           CF::ResultInfo::EntityUnknownType type) const;

  void GetRegionEquations(CF::StdVector<CF::UInt> & eqns,
                          CF::RegionIdType region) const;

protected:
  //! equation indices, 0 indicates "not mapped"
  CF::UInt* entityEquations_;

  //! Array of regions managed by the eqnMap
  CF::StdVector<CF::RegionIdType> regions_;

  //! pointer to grid associated to equations
  boost::shared_ptr<CF::Grid> ptGrid_;

  //! stores entry type of eqnMap
  CF::ResultInfo::EntityUnknownType mapType_;

  //! Equations per entity
  CF::UInt eqnPerEnt_;

  //! number of equations
  CF::UInt numEqns_;

  //! maximum number of equations
  CF::UInt maxNumEqns_;

  //! offset switch for zero and one based entity numbers
  CF::UInt zeroOne_;

  //! flag indicating if equation numbering is finalized
  bool isFinalized_;

};

}

#endif /* EQNNUMBERINGSIMPLE_HH_ */
