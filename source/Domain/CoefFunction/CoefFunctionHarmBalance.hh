// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// vim: set ts=2 sw=2 et nu ai ft=cpp cindent !:
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef COEF_FUNCTION_HARMBALANCE_HH
#define COEF_FUNCTION_HARMBALANCE_HH

#define _USE_MATH_DEFINES
#include <cmath>

#include <string>
#include <map>
#include <boost/tr1/type_traits.hpp>


#include "CoefFunction.hh"
#include "FeBasis/BaseFE.hh"
#include "Materials/BaseMaterial.hh"
#include "FeBasis/FeFunctions.hh"
#include "Utils/mathParser/mathParser.hh"

// To switch between time- and frequency domain
#include "Domain/Results/MHTimeFreqResult.hh"


namespace CoupledField  {

// ==========================================================================
//  COEFFICIENT FUNCTION HARMONIC BALANCE
// ==========================================================================

//!   PLEASE DESCRIBE ME BRIEFLY
//
//!   OI! DESCRIBE ME IN DETAIL

template <class T>
class CoefFunctionHarmBalance : public CoefFunction {

public:
  //! Constructor
  CoefFunctionHarmBalance(shared_ptr<BaseFeFunction> feFct,
                          shared_ptr<FeSpace> feSpc,
                          const StdVector<RegionIdType>& regions,
                          std::map<RegionIdType, BaseMaterial*>& materials,
                          Grid* ptGrid,
                          PtrCoefFct magFluxCoef,
                          const UInt& N,
                          const UInt& M,
                          const Double& baseFreq,
                          const UInt& nFFT);

  //! Destructor
  virtual ~CoefFunctionHarmBalance();


  //! \copydoc CoefFunction::GetScalar
  virtual void GetScalar(T& coefScal, const LocPointMapped& lpm );

  //! \copydoc CoefFunction::GetVector
  virtual void GetVector(Vector<T>& coefVec, const LocPointMapped& lpm);

  //! \copydoc CoefFunction::GetTensor
  virtual void GetTensor(Matrix<T>& coefMat, const LocPointMapped& lpm);

  //! \copydoc CoefFunction::ToString
  virtual std::string ToString() const;

  //! \copydoc CoefFunction::GetVecSize
  virtual UInt GetVecSize() const;

  PtrCoefFct GenerateMatCoefFnc(const UInt& iRegion,
                                const std::string& name,
                                const bool nonLin,
                                shared_ptr<ElemList> actSDList);

  virtual CoefDimType GetDimType() const{
    return dimType_;
  }

protected:

  void UpdateHarm();

  void UpdateSolution();

  //! DESCRIBE ME
  //! ===========================================
  struct HBRegionHelper
  {
    //! Store the elements of each region (need information
    //! about the integration points later on)
    shared_ptr<ElemList> elemListPerRegion;

    //! Attribute handling info on material data. From StdPDE
    BaseMaterial* material = NULL;

    //! CoefFunction for nonlinear reluctivity evaluation
    PtrCoefFct nonLinNuCoefMap;

    //! Region, the PDE is defined on. From StdPDE
    RegionIdType region;
  };


  //! Method is called with a specific subdomain in the PDE
  //! and registers this region
  //! Method is called in the integrator-definition section of the
  //! specific PDE
  void RegisterElemsInRegion(shared_ptr<ElemList> actSDList,
                             const UInt& iRegion,
                             HBRegionHelper& regStruc);


  //! FeFunction containing the coefficients
  shared_ptr<FeFunction<T> > feFct_;

  //! Pointer to math parser instance
  MathParser* mp_;

  //! Handle for the harmonic callback-mechanism
  MathParser::HandleType harmHandle_;

  //! Handle for the solution cache callback-mechanism
  MathParser::HandleType solHandle_;

  //! CoefFunction for magnetic flux density
  PtrCoefFct magFluxCoef_ = NULL;

  //! CoefFunction for nonlinear reluctivity evaluation
  //std::map<UInt, PtrCoefFct> nonLinNuCoefMap_;

  //! Vector containing all regions the PDE is defined on. From StdPDE
//  StdVector<RegionIdType> regions_;

  //! Attribute handling info on material data
  //! Maps regions and (simple) materials. From StdPDE
  //std::map<RegionIdType, BaseMaterial*> materials_;

  //! Pointer to grid object
  Grid * ptGrid_;

  //! Result-vector in frequency domain. Outer vector is frequency (harmonic),
  //! inner one contains the spatial dof's
  //StdVector< Vector<Complex> > freqResult_;
  Matrix<Complex> freqResult_;

  //! Number of regions
  UInt numRegions_;

  //! Pointer to the FeSpace for quantity MAG_FLUX_DENSITY,
  //! in order to retrieve information about integration points
  shared_ptr<FeSpace> fluxFESpace_;

  //! Safety flag to check if the elements of all regions were registered
  bool regionRegistration_;

  //! Store the elements of each region (need information
  //! about the integration points later on)
  //StdVector< shared_ptr<ElemList> > elemListPerRegion_;
  //StdVector< Integer > regionList_;

  //! Total number of elements
  UInt numElems_;

  //! Store how many times, the solution vector was read
  UInt updateIter_;

  unsigned int maxInt_;

  MHTimeFreqResult freqTimeRes_;

  //! For every region we create one HBRegionHelper struct.
  //! The iRegion variable is the key (comes from PDE)
  StdVector<HBRegionHelper> hbRegion_;

};
} //end of namespace

#endif
