// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_PIEZOCOUPLING_HH
#define FILE_PIEZOCOUPLING_HH

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
  class PiezoCoupling : public BasePairCoupling
  {
  public:
    //! Constructor
    //! \param pde1 pointer to first coupling PDE
    //! \param pde2 pointer to second coupling PDE
    //! \param paramNode pointer to "couplinglist/direct/piezoDirect" element
    PiezoCoupling( SinglePDE *pde1, SinglePDE *pde2, PtrParamNode paramNode );

    //! Destructor
    virtual ~PiezoCoupling();

    //! Trigger calculation of postprocessing results
    void CalcResults( shared_ptr<BaseResult> result );

  protected:

    //! Definition of the (bi)linear forms
    void DefineIntegrators();

    //! Define available results
    void DefineAvailResults();
    
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
