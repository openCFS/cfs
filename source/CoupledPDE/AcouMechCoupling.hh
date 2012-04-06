// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_ACOUMECHCOUPLING_HH
#define FILE_ACOUMECHCOUPLING_HH

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

  //! Implements the definition of pairwise piezo-coupling
  
  //! This class implements the piezoelectric coupling based on the 
  //! volumetric material law.
  //! Within this object, pde1_ refers to the mechanical PDE,
  //! whereas  pde2_ refers to the electric PDE.
  class AcouMechCoupling : public BasePairCoupling
  {
  public:
    //! Constructor
    //! \param pde1 pointer to first coupling PDE
    //! \param pde2 pointer to second coupling PDE
    //! \param paramNode pointer to "couplinglist/direct/piezoDirect" element
    AcouMechCoupling( SinglePDE *pde1, SinglePDE *pde2, PtrParamNode paramNode );

    //! Destructor
    virtual ~AcouMechCoupling();

    //! Trigger calculation of postprocessing results
    void CalcResults( shared_ptr<BaseResult> result );

  protected:

    //! Definition of the (bi)linear forms
    void DefineIntegrators();

    //! Define available results
    void DefineAvailResults();

    //! define all (bilinearform) integrators needed for this pdewith template
    //! for the space dimension
    void DefinePresOrPotIntegrators(const std::string& name,
                                    FEMatrixType matType,
                                    shared_ptr<BaseFeFunction>& dispFct,
                                    shared_ptr<BaseFeFunction>& acouFct,
                                    shared_ptr<SurfElemList>& actSDList,
                                    const std::map< RegionIdType, PtrCoefFct >& coefFuncs,
                                    const std::set< RegionIdType >& acouRegions);
    
    //! Returns a stiffness integrator appropriate to the actual problem (e.g. 3D)
    BiLinearForm * GetStiffIntegrator( BaseMaterial* actSDMat,
                                       RegionIdType regionId,
                                       bool isComplex );


    //! Subtype of related mechanical PDE
    std::string subType_;
    
  private:

  };


} // end of namespace

#endif
