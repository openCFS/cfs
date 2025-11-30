// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NEWBASESMOOTHPDE
#define FILE_NEWBASESMOOTHPDE

#include <map>

#include "Domain/CoefFunction/CoefFunction.hh"
#include "SinglePDE.hh"
#include "Forms/BiLinForms/BiLinearForm.hh"

namespace CoupledField
{

  class LinearFormContext;
  class set;

  //! Class for smooth equation (no adaptivity)
  class SmoothPDE: public SinglePDE
  {

  public:

    /** Constructor. here we read integration parameters
     * @param aGrid pointer to grid */
    SmoothPDE( Grid *aGrid, PtrParamNode paramNode,
             PtrParamNode infoNode,
             shared_ptr<SimState> simState, Domain* domain );

    //!  Deconstructor
    virtual ~SmoothPDE();


    /** @param see SinglePDE::GetSubTensorType() */
    SubTensorType GetSubTensorType() const { return tensorType_; }


    /** return the von Mises matrix (stress^T * M * stress = von Mises Stress)
     * @param dim desired dimension. axis is ignored currently :( */
    const Matrix<double>& GetVonMisesMatrix(int dim);

protected:

    //! Initialize NonLinearities
    void InitNonLin();

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators();
    
    //! Defines the integrators needed for ncInterfaces
    void DefineNcIntegrators();

    //! define surface integrators needed for this pde (currently only ABC)
    void DefineSurfaceIntegrators( );

    //! Define all RHS linearforms for load / excitation 
    void DefineRhsLoadIntegrators(PtrParamNode input = PtrParamNode());
    
    //! define the SoltionStep-Driver
    void DefineSolveStep();

    //! Read special results definition
    void ReadSpecialResults();

    //! \copydoc SinglePDE::CreateFeSpaces
    virtual std::map<SolutionType, shared_ptr<FeSpace> >  CreateFeSpaces(const std::string&  formulation, PtrParamNode infoNode );

    /** Returns a stiffness integrator appropriate to the actual problem (e.g. 3D)
     * //! Definition of convective integrators (Pierce Operator)
     * @param actSDMat  region material
     * @param regionID  current region node
     * @param feFct     fe Function (primary)
     * @param isComplex indicates if the problem is real- or complex valued */
    BaseBDBInt* GetStiffIntegrator(BaseMaterial* actSDMat, RegionIdType regionId, shared_ptr<BaseFeFunction> feFct, bool isComplex);

    //! Return flux integrator used for Nitsche coupling
    template<typename DATA_TYPE>
    BiLinearForm* GetFluxIntegrator(PtrCoefFct scalCoefFucn, PtrCoefFct coefFuncPMLVec, Double factor,
                                    BiLinearForm::CouplingDirection cplDir, bool fluxOpA, bool icModes, bool preStress = false);

    //! Return penalty integrator used for Nitsche coupling
    template<typename DATA_TYPE>
    BiLinearForm* GetPenaltyIntegrator(PtrCoefFct scalCoefFunc, Double factor, BiLinearForm::CouplingDirection cplDir);

    //! Return strain operator 
    BaseBOperator* GetStrainOperator( bool isComplex, bool icModes);

    /** @see virtual SinglePDE::GetNativeSolutionType() */
    SolutionType GetNativeSolutionType() const { return SMOOTH_DISPLACEMENT; }
    // ========================
    // set solution information
    // ========================

    //! Define contact between two surfaces
    void ReadContact(StdVector<std::string>& surfList1,
                          StdVector<std::string>& surfList2,
                          StdVector<std::string>& volumeList,
                          StdVector<std::string>& contactLawList,
                          StdVector<bool>& useSurfaceMidpointsList);

  private:

    //! Initialize time stepping method
    void InitTimeStepping();

    //! Define available primary result types
    void DefinePrimaryResults();
    
    //! Define available postprocessing results
    void DefinePostProcResults();
    


    //! Stores the linear stiffness for each region
    std::map<RegionIdType, PtrCoefFct > regionStiffness_;

    //! Dimension of stresses
    UInt stressDim_;
    
    //! Tensor type
    SubTensorType tensorType_;

    StdVector<std::string> dofNames_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class
  // =========================================================================

  //! \class SmoothPDE
  //!
  //! \purpose
  //! This class defines the smooth field PDE and the according
  //! postprocessing methods.

#endif

} // end of namespace
#endif

