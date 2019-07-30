// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_MAG_BASE_PDE_HH
#define FILE_MAG_BASE_PDE_HH

/*
 * MagBasePDE
 * > Idea taken from NACS
 * > Collect all functions that are common to MagneticPDE and MagEdgePDE (eventually also interesting for MagEdgeMxiedAVPDE)
 *
 * > Currently only used for hysteresis related functions (other parts commented out)
 */

#include "SinglePDE.hh"

namespace CoupledField
{

  class MagBasePDE: public SinglePDE {

  public:

    //! Constructor
    MagBasePDE( Grid *aGrid, PtrParamNode paramNode,
      PtrParamNode infoNode,
      shared_ptr<SimState> simState, Domain* domain );

    //!  Destructor
    virtual ~MagBasePDE();

    /*
     * Common functions for NODAL and EDGE formulation
     */
    //! define the SoltionStep-Driver
    void DefineSolveStep();

    virtual void InitNonLin();

    //! Init the time stepping
    void InitTimeStepping();

    //! read special boundary conditions (coils, magnets)
    void ReadSpecialBCs();

    void DefineMagFluxDensity();

    // =======================================================================
    //   COILS
    // =======================================================================
    //@{ \name Attributes related to coils
    //! Map CoilID to coil definition
    std::map<Coil::IdType, shared_ptr<Coil> > coils_;

    //! Map regionIds to coil definitions
    typedef std::map<RegionIdType, shared_ptr<Coil> > CoilRegionMap;
    CoilRegionMap coilRegions_;

    //! Coefficients holding the current density for each coil
    std::map<RegionIdType, PtrCoefFct> coilCurrentDens_;

    //! Storage for CoefFunctions of external current density as source
    std::map<shared_ptr<Coil::Part>, PtrCoefFct> coilPartsExtJ_;

    //! Tells if there are coils excited by voltage
    bool hasVoltCoils_;
    //@}

    /*
     * Functions related to coils
     * Mostly copied from NACS (thanks and credits to Simetris!)
     */
    //! Get mehtod for specific coils. Needed e.g. by the SinglePDE for
    //! specifying coil results.
    shared_ptr<Coil> GetCoilById(const Coil::IdType& id);

    //! Query parameter object for information on coils
    void ReadCoils();

    //! define special integrators for coils
    //! addition compared to NACS-version: Scaling factor
    //! (needed for magstrict coupling where we have to multiply with -1.0)
    void DefineCoilIntegrators(Double scaling);

    //! This method returns the RHS-integrator for a given current density,
    //! which is basically a BU-Integrator with a given identity-operator
    //! depending on the formulation.
    virtual LinearForm* GetCurrentDensityInt( Double factor,
                                              PtrCoefFct coef ) = 0;

    // =======================================================================
    //   HYSTERESIS
    // =======================================================================
    /*
     * Objects and variables related to or required for hysteresis
     */
    //! CoefFunctions for postprocessing results
    shared_ptr<CoefFunctionMulti> hysteresisPolarization_;
    shared_ptr<CoefFunctionMulti> hysteresisMagnetization_;
    shared_ptr<CoefFunctionMulti> hysteresisFieldIntensity_;

    //! In case of hysteresis, we set the region approximation in InitHystCoefs which has to be called
    //! prior to DefineIntegrators (where the regionApproximation usually is set)
    bool regionApproxSet_;
    std::map<RegionIdType,shared_ptr<ElemList> > SDLists_;

    /*
     * Functions related to hysteresis
     */
    //! Check all regions for hysteretic material behavior;
    //!for each found hysteretic region create a CoefFunctionHyst
    //!and add it to the map hysteresisCoefs_
    void InitHystCoefs();

    PtrCoefFct GenerateHystereticCoefFunctions(RegionIdType actRegion){

      PtrCoefFct hystPol = hysteresisCoefs_->GetRegionCoef(actRegion);
      // shared_ptr<CoefFunction >
      PtrCoefFct curCoef = hystPol->GenerateMatCoefFnc("Reluctivity");

//      PtrCoefFct hystOutput = hystPol->GenerateOutputCoefFnc("MagPolarization");
//      hysteresisPolarization_->AddRegion( actRegion, hystOutput);
//
//      PtrCoefFct hystOutput2 = hystPol->GenerateOutputCoefFnc("MagMagnetization");
//      hysteresisMagnetization_->AddRegion( actRegion, hystOutput2);
//
//      PtrCoefFct hystOutput3 = hystPol->GenerateOutputCoefFnc("MagFieldIntensityHyst");
//      hysteresisFieldIntensity_->AddRegion( actRegion, hystOutput3);

      return curCoef;
    }

    virtual LinearForm* GetRHSHystInt( Double factor, PtrCoefFct rhsMag, bool fullEvaluation ) = 0;

    virtual BaseBDBInt* GeHystStiffInt( Double factor, PtrCoefFct tensorReluctivity ) = 0;

    // =======================================================================
    //   MISC
    // =======================================================================
    /*
     * Common variables for NODAL and EDGE formulation
     */
    //! Formulation type of magnetic PDE
    typedef enum { BASE, NODAL, EDGE } MagFormulationType;

    //!true, if 3d computation
    bool is3d_;

    //! formulation
    MagFormulationType formulation_;

    bool fluxDensityDefined_;

    //! In case of a transient / harmonic 3D simulation using NODAL formulation, we need
    //! an additional scalar potential to ensure HCurl compatibility (A-V-formulation)
    bool isMixed_;

    //! in mixed case, we have to check if there are any regions with conductivity at all
    bool anyRegionHasConductivity_;

    //! Coefficient function, containing the overall reluctivity
    shared_ptr<CoefFunctionMulti> relucTensor_;

    //! Coefficient function, containing the overall reluctivity
    shared_ptr<CoefFunctionMulti> reluc_;

    //! Coefficient function, containing the overall conductivity
    shared_ptr<CoefFunctionMulti> conduc_;

  };
} // end of namespace
#endif

//
//#include <map>
//
//#include "SinglePDE.hh"
//#include "Utils/Coil.hh"
//
//#include "Domain/CoefFunction/CoefFunctionMulti.hh"
//#include "Forms/Operators/GradientOperator.hh"
//#include "Domain/Results/ResultFunctor.hh"
//
//namespace NACS {
//
//
//  // Forward declaration of classes
//  class MagMechB0Coupling;
//  class BaseFeFunction;
//
//  //! Abstract class for magnetic field problems without specific formulation
//
//  //! This abstract base class for magnetic field problems encapsulates common
//  //! methods for coil handling and provides interfaces for calculating
//  //! quantities, independent of the underlying formulation (nodal, edge).
//  class MagBasePDE: public SinglePDE {
//
//  public:
//
//    //! Make MagMechB0 Coupling friend class
//    friend class MagMechB0Coupling;
//
//    //! Formulation type of magnetic PDE
//    typedef enum { NODAL, EDGE } MagFormulationType;
//
//    //! Constructor
//    MagBasePDE( Grid *aGrid, PtrParamNode paramNode,
//                PtrParamNode infoNode,
//                shared_ptr<SimState> simState, Domain* domain );
//
//    //!  Destructor
//    virtual ~MagBasePDE();
//
//    //! Get method for specific coils. Needed e.g. by the SinglePDE for
//    //! specifying coil results.
//    shared_ptr<Coil> GetCoilById(const Coil::IdType& id);
//
//  protected:
//
//    //! Initialize NonLinearities
//    void InitNonLin();
//
//    //! define special integrators for coils
//    void DefineCoilIntegrators();
//
//    //! Return RHS-integrator for current density excitation
//
//    //! This method returns the RHS-integrator for a given current density,
//    //! which is basically a BU-Integrator with a given identity-operator
//    //! depending on the formulation.
//    virtual LinearForm* GetCurrentDensityInt( Double factor,
//                                              PtrCoefFct coef ) = 0;
//
//    // =======================================================================
//    //  Initialization
//    // =======================================================================
//
//     //! Define fepsace for voltage driven coils
//    void DefineCoilFeSpace(std::map<SolutionType, shared_ptr<FeSpace> >& crSpaces,
//                          PtrParamNode infoNode);
//
//    //! Define primary coil results (in case of voltage-loaded coil)
//    void DefinePrimaryCoilResults();
//
//    //! Define available results for coils
//    void DefinePostProcCoilResults();
//
//    //! Finalize the postprocessing results (e.g. coil results)
//    void FinalizePostProcCoilResults();
//
//    //! Query parameter object for information on coils
//    void ReadCoils();
//
//    //! Define time-stepping for coils
//    void InitCoilTimeStepping();
//
//
//    // =======================================================================
//    //   COILS
//    // =======================================================================
//
//    //@{ \name Attributes related to coils
//    //! Map CoilID to coil definition
//    std::map<Coil::IdType, shared_ptr<Coil> > coils_;
//
//    //! Map regionIds to coil definitions
//    typedef std::map<RegionIdType, shared_ptr<Coil> > CoilRegionMap;
//    CoilRegionMap coilRegions_;
//
//    //! Coefficients holding the current density for each coil
//    std::map<RegionIdType, PtrCoefFct> coilCurrentDens_;
//
//    //! Flag if voltage-driven coils are present
//    bool hasVoltCoils_;
//    //@}
//
//    //!true, if 3d computation
//    bool is3d_;
//
//    //! formulation
//    MagFormulationType formulation_;
//
//    //! Coefficient function, containing the overall reluctivity (full tensor)
//    shared_ptr<CoefFunctionMulti> reluc_;
//
//    //! Coefficient function, containing the overall conductivity
//    shared_ptr<CoefFunctionMulti> conduc_;
//
//    //! Map with regions on which a current density is prescribed
//    std::map<RegionIdType, PtrCoefFct> explicitCurrentRegions_;
//
//    // Coefficient function containing the magnetization of permanent magnets
//    shared_ptr<CoefFunctionMulti> magnetization_;
//
//    //! Pointer to info node for coils
//    PtrParamNode coilInfo_;
//  };
//
//
//
//  // ========================================================================
//  //  Class for calculating coil results
//  // ========================================================================
//
//  //! ResultFunctor class for calculating various coil results
//  template<class TYPE>
//  class CoilResultFunctor : public ResultFunctor {
//  public:
//
//    //! Constructor
//    CoilResultFunctor( shared_ptr<MathParser> mp,
//                       Grid * ptGrid,
//                       shared_ptr<BaseFeFunction> feFct,
//                       shared_ptr<ResultInfo> info,
//                       const std::map<Coil::IdType, shared_ptr<Coil> >& coilMap,
//                       PtrCoefFct potential,
//                       PtrCoefFct potentialD1,
//                       shared_ptr<BaseFeFunction> currentFct,
//                       shared_ptr<CoefFunctionMulti> conductivity );
//
//    //! Destructor
//    virtual ~CoilResultFunctor();
//
//    //! Evaluate result for complete entity list
//    virtual void EvalResult( shared_ptr<BaseResult> res );
//
//
//  protected:
//
//    //! Math parser instance
//    shared_ptr<MathParser> mp_;
//
//    //! Math parser handle for expression "t"
//    MathParser::HandleType mphT_;
//
//    //! Time of the current step
//    Double curTime_;
//
//    //! Time of the previous step
//    Double prevTime_;
//
//    //! Grid pointer
//    Grid* ptGrid_;
//
//    //! Pointer to feSpace of vector potential
//    shared_ptr<BaseFeFunction> feFct_;
//
//    //! Map CoilID to coil definition
//    std::map<Coil::IdType, shared_ptr<Coil> > coils_;
//
//    //! Stores the current per coil for the current time step
//    std::map<shared_ptr<Coil>, TYPE> coilCurrent_;
//
//    //! Stores the current per coil for the previous time step
//    std::map<shared_ptr<Coil>, TYPE> prevCoilCurrent_;
//
//    //! CoefFunction for vector potential
//    PtrCoefFct potential_;
//
//    //! CoefFunction for vector potential derivative
//    PtrCoefFct potentialD1_;
//
//    //! CoefFunction for (A * e_J)
//    shared_ptr<CoefFunctionMulti> aTimesJ_;
//
//    //! CoefFunction for (dA/dt * e_J)
//    shared_ptr<CoefFunctionMulti> aD1TimesJ_;
//
//    //! Coefficient function, containing the overall conductivity
//    shared_ptr<CoefFunctionMulti> conduc_;
//
//    //! Pointer to feFcuntion for coil currents
//    shared_ptr<BaseFeFunction> currentFunc_;
//
//    //! Result functor for integrating the flux
//    shared_ptr<ResultFunctorIntegrate<TYPE> > fluxInt_;
//
//    //! Result functor for integartion the first derivative of the flux
//    shared_ptr<ResultFunctorIntegrate<TYPE> > fluxD1Int_;
//
//
//    //! Calculate the current for a given coil
//    TYPE CalcSingleCoilCurrent( EntityIterator& coilIt ) ;
//
//    //! Calculate the voltage for a given coil
//    TYPE CalcSingleCoilVoltage( EntityIterator& coil );
//
//    //! Calculate the inductance for a given coil
//    TYPE CalcSingleCoilInductance( EntityIterator& coilIt );
//
//    //! Calculate the flux / derivative of a given coil
//    TYPE CalcSingleCoilFlux( EntityIterator& coilIt,
//                             bool timeDeriv );
//  };
//
//  // ========================================================================
//  //  CoefFunction for current density of voltage driven coils
//  // ========================================================================
//
//  //! Coefficient function for current density of voltage driven coil
//  template<class TYPE>
//  class CoefFunctionJVoltCoil :
//    public CoefFunction,
//    public boost::enable_shared_from_this<CoefFunctionJVoltCoil<TYPE> > {
//  public:
//    //! Constructor
//    CoefFunctionJVoltCoil( Grid* ptGrid,
//                           shared_ptr<Coil> coil,
//                           shared_ptr<BaseFeFunction> currentFnc );
//
//    //! Destructor
//    virtual ~CoefFunctionJVoltCoil();
//
//    //! \copydoc CoefFunction::GetVector
//    virtual void GetVector( Vector<TYPE>& vec,
//                            const LocPointMapped& lpm );
//
//    //! \copydoc CoefFunction::GetVecSize
//    virtual UInt GetVecSize() const { return dim_;}
//
//    //! \copydoc CoefFunction::IsZero
//    bool IsZero() const { return false;}
//
//    //! \copydoc CoefFunction::ToString
//    std::string ToString() const;
//
//  private:
//
//    //! Dimension
//    UInt dim_;
//
//    //! Point to coil
//    shared_ptr<Coil> coil_;
//
//    //! Coil list
//    shared_ptr<CoilList> coilList_;
//
//    //! Pointer to feFcuntion for coil currents
//    shared_ptr<BaseFeFunction> currentFunc_;
//  };
//
//  // ========================================================================
//  // CoefFunction for magnetic field intensity
//  // ========================================================================
//
//  //! CoefFunction for magnetic field intensity with correction for permanent magnets
//  template<class TYPE>
//  class CoefFunctionMagFieldIntensity :
//      public CoefFunction,
//      public boost::enable_shared_from_this< CoefFunctionMagFieldIntensity<TYPE> >
//  {
//    public:
//
//      //! Constructor
//      CoefFunctionMagFieldIntensity(PtrCoefFct hField,
//                                    PtrCoefFct magnetization,
//                                    PtrCoefFct reluctivity);
//
//      virtual ~CoefFunctionMagFieldIntensity() {};
//
//      //! \copydoc CoefFunction::GetVector
//      virtual void GetVector( Vector<TYPE>& vec, const LocPointMapped& lpm );
//
//      //! \copydoc CoefFunction::GetVecSize
//      virtual UInt GetVecSize() const;
//
//    private:
//
//      //! Pointer to magnetic field intensity function
//      PtrCoefFct hField_;
//
//      //! Pointer to magnetization function for permanent magnets
//      PtrCoefFct magnetization_;
//
//      //! Pointer to reluctivity function
//      PtrCoefFct reluc_;
//  };
//
//  //! This class represents coefficient functions, which are defined just on a
//  //! surface and computes the force defined by Maxwell's stress tensor
//  //! It#s derived from CoefFunctionSurf
//  template<class FE>
//  class CoefFunctionMaxwell : public CoefFunction {
//  public:
//
//    //! Constructor
//
//    //! Constructor for the class
//    //! \param mapNormal If true, only the normal component w.r.t. to the
//    //!                  surface element is taken into account.By default,
//    //!                  the normal direction points OUT of the related volumes.
//    //! \param material parameter
//    //! \param factor Additional scaling factor
//    //! \param surfInfo Result info object for surface result
//
//    CoefFunctionMaxwell(PtrCoefFct Bfield,
//                        std::map<SolutionType, shared_ptr<CoefFunctionMulti> > matCoefs,
//                        Grid* ptGrid, MagBasePDE::MagFormulationType formulation,
//                        shared_ptr<FeSpace> feSpace);
//
//    //! Destructor
//    virtual ~CoefFunctionMaxwell() {}
//
//    //! \copydoc CoefFunction::GetTensor
//    void GetTensor(Matrix<Double>& coefMat, const LocPointMapped& lpm ) {
//      EXCEPTION("CoefFunctionSurfMaxwell:GetTensor not implemented");
//    }
//
//    //! \copydoc CoefFunction::GetTensor
//    void GetTensor(Matrix<Complex>& coefMat, const LocPointMapped& lpm ) {
//      EXCEPTION("CoefFunctionSurfMaxwell:GetTensor not implemented");
//    }
//
//    //! \copydoc CoefFunction::GetVector
//    void GetVector(Vector<Double>& coefVec, const LocPointMapped& lpm);
//
//    //! \copydoc CoefFunction::GetVector
//    void GetVector(Vector<Complex>& coefVec, const LocPointMapped& lpm );
//
//    //! \copydoc CoefFunction::GetScalar
//    void GetScalar(Double& coefScalar, const LocPointMapped& lpm ) {
//        EXCEPTION("CoefFunctionSurfMaxwell:GetScalar not implemented");
//    }
//
//    //! \copydoc CoefFunction::GetScalar
//    void GetScalar(Complex& coefScalar, const LocPointMapped& lpm ) {
//        EXCEPTION("CoefFunctionSurfMaxwell:GetScalar not implemented");
//    };
//
//    //! Return row and columns size of tensor if coefficient function is a tensor
//    UInt GetVecSize() const {
//      return ptGrid_->GetDim();
//    }
//
//  private:
//
//    //! Pointer to B field
//    PtrCoefFct bField_;
//
//    //! coef-function as defined in PDE
//    std::map<SolutionType, shared_ptr<CoefFunctionMulti> > matCoef_;
//
//    //! pointer to the grid
//    Grid* ptGrid_;
//
//    //! Pointer to FeSpace of PDE
//    shared_ptr<FeSpace> feSpace_;
//  };
//
//} // end of namespace
//#endif
