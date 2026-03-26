// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NEWBASEMECHPDE
#define FILE_NEWBASEMECHPDE

#include <map>

#include "Domain/CoefFunction/CoefFunction.hh"
#include "SinglePDE.hh"
#include "Forms/BiLinForms/BiLinearForm.hh"

namespace CoupledField
{

  class LinearFormContext;
  class set;

  //! Class for mechanic equation (no adaptivity)
  class MechPDE: public SinglePDE
  {

  public:

    /** Constructor. here we read integration parameters
     * @param aGrid pointer to grid */
    MechPDE( Grid *aGrid, PtrParamNode paramNode,
             PtrParamNode infoNode,
             shared_ptr<SimState> simState, Domain* domain );

    //!  Deconstructor
    virtual ~MechPDE();

    /** constants for test-strains, used for homogenization. We depend on the int values! */
    typedef enum { NOT_SET=-1, X=0, Y=1, Z=2, YZ=3, XZ=4, XY=5 } TestStrain;

    /** @param see SinglePDE::GetSubTensorType() */
    SubTensorType GetSubTensorType() const { return tensorType_; }


    /** return the von Mises matrix (stress^T * M * stress = von Mises Stress)
     * @param dim desired dimension. axis is ignored currently :( */
    const Matrix<double>& GetVonMisesMatrix(int dim);

    /** return the Hill Mandel matrix (stress^T * M * stress = Hill Mandel norm)
     * @param dim desired dimension. */
    const Matrix<double>& GetHillMandelMatrix(int dim);

    /** Add the integrators for the test strains for homogenization to the linear forms, similar as in multiple load case;
     * called from Excitation::ReadLoads or Excitation::SetHomogenizationTestStrains() (optimization)
     * @param test is an enum
     * @param linForms set to append linear Forms to, if NULL use assemble_ */
    void DefineTestStrainIntegrator(const TestStrain test, StdVector<LinearFormContext*>* linForms = NULL);

    static Enum<TestStrain> testStrain;

protected:

    //! read in damping information, see SinglePDE.cc  and SinglePDE.hh
    void ReadDampingInformation();

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
     * @param isComplex either from complex material or bloch mode */
    BaseBDBInt* GetStiffIntegrator(BaseMaterial* actSDMat, RegionIdType regionId, bool isComplex);

    /** Returns a stiffness integrator appropriate to the actual problem (e.g. 3D) with the material tensor scaled by a given factor
     * @param isComplex either from complex material or bloch mode
     * @param scalingFactor is a factor the material tensor to be multiplied by */
    BaseBDBInt* GetStiffIntegrator(BaseMaterial* actSDMat, RegionIdType regionId, bool isComplex, PtrCoefFct scalingFactor);

    /** Returns a Kelvin-Voigt damping integrator appropriate to the actual problem (e.g. 3D)
     * @param isComplex either from complex material or bloch mode */
    BaseBDBInt *GetKelvinVoigtDampingIntegrator(BaseMaterial *actSDMat, RegionIdType regionId, bool isComplex);

    /** Returns an integrator for prestressing
     * @param preStressFct is a tensor coefficient function for prestressing
     * @param isComplex either from complex material or bloch mode */
    BaseBDBInt* GetPreStressIntegrator(PtrCoefFct preStressFct, RegionIdType regionId, bool isComplex, Double factor = 1.0);

    /** Returns an integrator for prestressing
     * @param preStressFct is a tensor coefficient function for prestressing
     * @param isComplex either from complex material or bloch mode
     * @param scalingFactor is a factor the prestressing tensor to be multiplied by */
    BaseBDBInt* GetPreStressIntegrator(PtrCoefFct preStressFct, RegionIdType regionId, bool isComplex, PtrCoefFct scalingFactor);

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
    SolutionType GetNativeSolutionType() const { return MECH_DISPLACEMENT; }
    // ========================
    // set solution information
    // ========================

  private:

    //! Initialize time stepping method
    void InitTimeStepping();
    
    //! read in softening types
    void ReadSoftening();

    //! Define available primary result types
    void DefinePrimaryResults();
    
    //! Define available postprocessing results
    void DefinePostProcResults();
    
    //! Define concentrated mechanical elements
    //! Models a node attached to a spring. The given value is the spring constant in N/m which is added directly to the appropriate stiffness matrix entry. example: testsuite Optimization/Static/output
    void DefineConcentratedElems();

    //! Create CoefFunction for preStressing
    PtrCoefFct CreatePreStressFct( bool isComplex, PtrParamNode stressNode);

    //! Returns a MechStress CoefFunction from a preceding sequence step
    PtrCoefFct GetStressCoefFromSeqStep(UInt seqStep);

    //! Stores softening for each region
    std::map<RegionIdType, std::string> regionSoftening_;
    
    //! Stores Rayleigh damping definition for each region
    std::map<RegionIdType, RaylDampingData > regionRaylDamping_;

    //! Stores Kelvin-Voigt damping definition (viscous tensor as coef function) for each region
    std::map<RegionIdType, PtrCoefFct > regionKelvinVoigtDamping_;

    //! Stores the linear stiffness for each region
    std::map<RegionIdType, PtrCoefFct > regionStiffness_;

    //! Stores the prestressing for each region
    std::map<RegionIdType, PtrCoefFct > regionPreStress_;

    //! Dimension of stresses
    UInt stressDim_;
    
    //! Tensor type
    SubTensorType tensorType_;

    //! Returns SubTensorType based on string
    SubTensorType GetSubTensorType();

    //! coefFunctzion for thermal strain
    shared_ptr<CoefFunctionMulti> thermalStrain_;

    //! coefFunctzion for thermal stress
    shared_ptr<CoefFunctionMulti> thermalStress_;

    /** static matrix which allows the calculation of the von Mises Stress from the stresses in 2D.
     * See Kocvara and Stingl; 2007 */
    Matrix<double> vonMisesMatrix_2d_;
    Matrix<double> vonMisesMatrix_3d_;

    /** static matrix which allows the calculation of the Hill Mandel norm of a vector. */
    Matrix<double> hillMandelMatrix_2d_;
    Matrix<double> hillMandelMatrix_3d_;

    StdVector<std::string> dofNames_;

    // ========================================================================
    //  Class for calculating the 2nd Piola-Kirchhoff stress tensor
    // ========================================================================

    //! CoefficientFunction for 2nd Piola-Kirchhoff stress tensor
    class CoefFunction2ndPiolaTensor : public CoefFunction {
        
      public:
        
        //! Constructor
        CoefFunction2ndPiolaTensor(SubTensorType &subType,
                                   PtrCoefFct stiffness,
                                   shared_ptr<BaseFeFunction> displ);
        
        //! Destructor
        ~CoefFunction2ndPiolaTensor();
        
        virtual string GetName() const { return "CoefFunction2ndPiolaTensor"; }

        //! Return tensor in Voigt vector notation
        virtual void GetVector(Vector<Double>& vec, 
                               const LocPointMapped& lpm );
        
        //! Return full tensor at integration point
        virtual void GetTensor(Matrix<Double>& tensor, 
                               const LocPointMapped& lpm );
        
        //! Return row and columns size of tensor if coefficient function is a tensor
        virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const;
        
      protected:
        
        //! Tensor type
        SubTensorType tensorType_;
        
        //! Stiffness CoefFunction
        PtrCoefFct stiffCoef_;
        
        //! Real-valued displacement CoefFunction
        shared_ptr< FeFunction<Double> > dispCoefReal_;
        
        //! Complex-values displacement CoefFunction
        shared_ptr< FeFunction<Complex> > dispCoefComplex_;
        
        //! Linear strain operator
        BaseBOperator *linOp_;
        
        //! Nonlinear strain operator
        BaseBOperator *nonLinOp_;
    };
    
  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class
  // =========================================================================

  //! \class MechPDE
  //!
  //! \purpose
  //! This class defines the mechanical field PDE and the according
  //! postprocessing methods.

#endif

} // end of namespace
#endif

