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
  CoefFunctionHarmBalance();

  //! Destructor
  virtual ~CoefFunctionHarmBalance();

  virtual void Init(shared_ptr<BaseFeFunction> feFct,
          shared_ptr<FeSpace> feSpc,
          const StdVector<RegionIdType>& regions,
          std::map<RegionIdType, BaseMaterial*>& materials,
          Grid* ptGrid,
          PtrCoefFct magFluxCoef,
          const UInt& N,
          const UInt& M,
          const Double& baseFreq,
          const UInt& nFFT,
          std::string modelName,
          PtrCoefFct matModelCoef);



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

  void FinishCash();

  void CashTimeResult();

  //! ===========================================
  //! ===========================================
  //! DESCRIBE ME
  struct HBRegionHelper
  {
    //! Store the elements of each region (need information
    //! about the integration points later on)
    shared_ptr<ElemList> elemListPerRegion;

    //! Attribute handling info on material data. From StdPDE
    BaseMaterial* material = NULL;

    //! CoefFunction for nonlinear reluctivity evaluation
    PtrCoefFct nonLinNuCoefMap;

    //! CoefFunction for linear reluctivity evaluation
    PtrCoefFct linNuCoefMap;

    //! Region, the PDE is defined on. From StdPDE
    RegionIdType region;

    bool isNonLin;
  };
  //! DESCRIBE ME
  //! ===========================================
  //! ===========================================


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

  //! Handle for the cashing callback-mechanism
  MathParser::HandleType cashHandle_;

  //! Handle for the solution cache callback-mechanism
  MathParser::HandleType solHandle_;

  //! Handle for the harmonic callback-mechanism
  MathParser::HandleType harmonicHandle_;

  //! Pointer to grid object
  Grid * ptGrid_;

  //! CoefFunction for magnetic flux density
  PtrCoefFct BField = NULL;

  PtrCoefFct magFluxCoef_ = NULL;

  //! Number of regions
  UInt numRegions_;

  UInt numRegisteredRegions_;

  //! Pointer to the FeSpace for quantity MAG_FLUX_DENSITY,
  //! in order to retrieve information about integration points
  shared_ptr<FeSpace> fluxFESpace_;

  //! Safety flag to check if the elements of all regions were registered
  bool regionRegistration_ = false;

  //! Total number of elements
  UInt numElems_;

  //! Store how many times, the solution vector was read
  UInt updateIter_;

  int maxInt_;

  //! Number of harmonics
  UInt N_;

  UInt M_;

  MHTimeFreqResult freqTimeRes_;

  //! Boolean to check if we even have a multiharmonic analysis
  bool isMH_;

  //! For every region we create one HBRegionHelper struct.
  //! The iRegion variable is the key (comes from PDE)
  StdVector<HBRegionHelper> hbRegion_;

  //! Temporary vector nuFreqTmp_[i] contains the nu(timestep) of element nuFreqTmpELEM_[i]
  Vector<Double> nuFreqTmp_;

  //! Map for getting the position of a certain elemNum in the freqTimeRes_ construct
  boost::unordered_map<UInt, UInt> positionOfElem_;

  std::string modelName_;

  PtrCoefFct matModelCoef_;
};



template <class T>
class CoefFunctionHarmBalanceEval : public CoefFunction {

public:
  //! Constructor
  CoefFunctionHarmBalanceEval(PtrCoefFct magFluxCoef,
                              Grid* ptGrid);

  //! Destructor
  virtual ~CoefFunctionHarmBalanceEval();


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

private:

  //! Pointer to math parser instance
  MathParser* mp_;

  //! Pointer to grid object
  Grid * ptGrid_;

  //! CoefFunction for magnetic flux density
  PtrCoefFct magFluxCoef_ = NULL;

};
} //end of namespace

#endif
