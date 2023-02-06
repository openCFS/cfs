// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     LinFlowHeatCoupling.hh
 *       \brief    Coupling between linearized flow equations and heat equation
 *       \date     Oct 24, 2018
 *       \author   Manfred Kaltenbacher
 */
//================================================================================================

#ifndef FILE_LINFLOWHEATCOUPLING_HH
#define FILE_LINFLOWHEATCOUPLING_HH

#include "BasePairCoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/defs.hh"
#include "General/Environment.hh"

namespace CoupledField {
class BaseResult;
class EntityIterator;
template <class TYPE> class Matrix;
template <class TYPE> class Vector;

// Forward declarations
class BaseMaterial;
class SinglePDE;
class BiLinearForm;

  //! Implements the definition of pairwise linearized flow - heat coupling
  
  //! This class implements the linearized flow - heat coupling
  //! Within this object, pde1_ refers to the linear flow PDE,
  //! whereas  pde2_ refers to the heatConduction PDE.
  //! This is fixed in Domian.cc / CreateDirectCoupledPDEs

  class LinFlowHeatCoupling : public BasePairCoupling
  {
  public:
    //! Constructor
    //! \param pde1 pointer to linearized flow PDE
    //! \param pde2 pointer to heat PDE
    //! \param paramNode pointer to "couplinglist/direct/piezoDirect" element
	LinFlowHeatCoupling( SinglePDE *pde1, SinglePDE *pde2, PtrParamNode paramNode,
                      PtrParamNode infoNode,
                      shared_ptr<SimState> simState, Domain* domain );

    //! Destructor
    virtual ~LinFlowHeatCoupling();

  protected:

    //! Definition of the (bi)linear forms
    void DefineIntegrators();

    //! Subtype of related mechanical PDE
    std::string subType_;
    
  private:

    //! Whether to use a symmetric formulation or not
    bool isCouplingFormulationSymmetric_;

  };


} // end of namespace

#endif
