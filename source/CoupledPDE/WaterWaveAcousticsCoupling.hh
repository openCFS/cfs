// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_WATERWAVESACOUSTICSCOUPLING_HH
#define FILE_WATERWAVESACOUSTICSCOUPLING_HH

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
  
  //! This class implements the coupling between water waves 
  //! and acoustic waves in air via a Lagrange multiplier
  //! Within this object, pde1_ refers to the water wave PDE,
  //! whereas  pde2_ refers to the acoustic PDE.
  class WaterWaveAcousticCoupling : public BasePairCoupling
  {
  public:
    //! Constructor
    //! \param pde1 pointer to first coupling PDE
    //! \param pde2 pointer to second coupling PDE
    //! \param paramNode pointer to "couplinglist/direct/piezoDirect" element
    WaterWaveAcousticCoupling( SinglePDE *pde1, SinglePDE *pde2, PtrParamNode paramNode,
                       PtrParamNode infoNode,
                       shared_ptr<SimState> simState, Domain* domain );

    //! Destructor
    virtual ~WaterWaveAcousticCoupling();

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

    //! Create a particular bilinear form integrator
    void DefCouplInt( const std::string& name,
                      bool assembleTransposed,
                      Double factor,
                      FEMatrixType matType,
                      shared_ptr<BaseFeFunction>& fnc1,
                      shared_ptr<BaseFeFunction>& fnc2,
                      shared_ptr<SurfElemList>& actSDList,
                      const std::map< RegionIdType, PtrCoefFct >& coefFuncs,
                      const std::set< RegionIdType >& acouRegions) ;
    
  private:
    //! solution type
    SolutionType formulation_;

    //! results
    ResultInfoList results_;

  };


} // end of namespace

#endif
