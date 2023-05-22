// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
// ================================================================================================
/*!
 *       \file     LinFlowMechCoupling.hh
 *       \brief    Coupling between linearized flow equations and mechanical equation
 *       \date     Feb 08, 2019
 *       \author   Manfred Kaltenbacher
 */
//================================================================================================

#ifndef FILE_LINFLOWMECHCOUPLING_HH
#define FILE_LINFLOWMECHCOUPLING_HH

#include "BasePairCoupling.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
#include "Forms/BiLinForms/BiLinearForm.hh"

namespace CoupledField {
class BaseResult;
class EntityIterator;
template <class TYPE> class Matrix;
template <class TYPE> class Vector;

// Forward declarations
class BaseMaterial;
class SinglePDE;
class BiLinearForm;
class MortarInterface;

  //! Implements the definition of pairwise piezo-coupling
  
  //! This class implements the piezoelectric coupling based on the 
  //! volumetric material law.
  //! Within this object, pde1_ refers to the mechanical PDE,
  //! whereas  pde2_ refers to the electric PDE.
  class LinFlowMechCoupling : public BasePairCoupling
  {
  public:
    //! Constructor
    //! \param pde1 pointer to first coupling PDE
    //! \param pde2 pointer to second coupling PDE
    //! \param paramNode pointer to "couplinglist/direct/piezoDirect" element
    LinFlowMechCoupling( SinglePDE *pde1, SinglePDE *pde2, PtrParamNode paramNode,
                       PtrParamNode infoNode,
                       shared_ptr<SimState> simState, Domain* domain );

    //! Destructor
    virtual ~LinFlowMechCoupling();

  protected:

    //! Definition of the (bi)linear forms
    void DefineIntegrators();

    //! Define available results
    void DefineAvailResults();

    //! Define available results
    void DefinePrimaryResults();

    //! Create FeSpaces according to formulation
     virtual void CreateFeSpaces( const std::string&  type,
                                  PtrParamNode infoNode,
                                  std::map<SolutionType, shared_ptr<FeSpace> >& crSpaces);

    //! define all (bilinearform) integrators needed for this pdewith template
    //! for the space dimension
    void DefineDampingIntegrators(const std::string& name,
                               shared_ptr<BaseFeFunction>& dispFct,
                               shared_ptr<BaseFeFunction>& lmFct,
                               shared_ptr<SurfElemList>& actSDList,
                               const std::map< RegionIdType, PtrCoefFct >& muOverDensityFuncs,
                               const std::set< RegionIdType >& flowRegions);
    //! define all (bilinearform) integrators needed for this pdewith template
    //! for the space dimension
    void DefineStiffnessIntegrators(const std::string& name,
                                    shared_ptr<BaseFeFunction>& dispFct,
                                    shared_ptr<BaseFeFunction>& velFct,
                                    shared_ptr<BaseFeFunction>& lmFct,
                                    shared_ptr<SurfElemList>& actSDList,
                                    const std::map< RegionIdType, PtrCoefFct >& densityFuncs,
                                    const std::map< RegionIdType, PtrCoefFct >& muFuncs,
                                    const std::map< RegionIdType, PtrCoefFct >& oneFuncs,
                                    const std::set< RegionIdType >& flowRegions);

    //! Define Mortar integrator (coupling across NC interface)
    void DefineMortarIntNC(const std::string& name,
                           shared_ptr<BaseFeFunction>& fct1,
                           shared_ptr<BaseFeFunction>& fct2,
                           MortarInterface* mortarIf,
                           Double factor,
                           PtrCoefFct scalCoef,
                           FEMatrixType matType,
                           BiLinearForm::CouplingDirection cplDir);

    //! Define Mortar integrator (FeFunctions on same surface)
    void DefineMortarIntNCSecondary(const std::string& name,
                              shared_ptr<BaseFeFunction>& fct1,
                              shared_ptr<BaseFeFunction>& fct2,
                              shared_ptr<SurfElemList>& actSDList,
                              Double factor,
                              PtrCoefFct scalCoef,
                              FEMatrixType matType);

   //! Subtype of related mechanical PDE
    std::string subType_;
    
  private:
    //! solution type
    SolutionType formulation_;

    //! results
    ResultInfoList results_;

    //! Lagrange multiplier order is same as velocity?
    bool lmOrderSameAsVel_;
    bool IsLagrangeMultiplierMethod_;
    bool hasMortarIface_;
    bool IncludeSymmetrizationTerm_;
    
	//! polyIds for the LinFlow side
    std::string presPolyId_;
    std::string velPolyId_;
    std::string lagrangeMultPolyId_;
    
    std::string presIntegId_;
    std::string velIntegId_;
    std::string lagrangeMultIntegId_;
  };


} // end of namespace

#endif
