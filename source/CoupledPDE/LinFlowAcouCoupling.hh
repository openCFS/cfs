// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     LinFlowAcouCoupling.hh
 *       \brief    Coupling between linearized flow equations and wave equation
 *       \date     Jan 11, 2019
 *       \author   Manfred Kaltenbacher
 */
//================================================================================================

#ifndef FILE_LINFLOWACOUCOUPLING_HH
#define FILE_LINFLOWACOUCOUPLING_HH

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
  
  //! This class implements the linearized flow - acoustic coupling
  //! Within this object, pde1_ refers to the linear flow PDE,
  //! whereas  pde2_ refers to the acoustic PDE.
  //! This is fixed in Domian.cc / CreateDirectCoupledPDEs

  class LinFlowAcouCoupling : public BasePairCoupling
  {
  public:
    //! Constructor
    //! \param pde1 pointer to linearized flow PDE
    //! \param pde2 pointer to heat PDE
    //! \param paramNode pointer to "couplinglist/direct/piezoDirect" element
	LinFlowAcouCoupling( SinglePDE *pde1, SinglePDE *pde2, PtrParamNode paramNode,
                      PtrParamNode infoNode,
                      shared_ptr<SimState> simState, Domain* domain );

    //! Destructor
    virtual ~LinFlowAcouCoupling();

  protected:

    //! Definition of the (bi)linear forms
    void DefineIntegrators();

    //! Subtype of related mechanical PDE
    std::string subType_;
    
  private:

  };


} // end of namespace

#endif
