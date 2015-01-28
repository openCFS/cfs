// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NEWBASEMECHPDE
#define FILE_NEWBASEMECHPDE

#include <map>

#include "Domain/CoefFunction/CoefFunction.hh"
#include "SinglePDE.hh"

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

    /** small helper which translates the test strain code on a vector of size 3/6 with one entry 1.0, the other zero
     * @param reduced in 2d case the result size is 3 otherwise 6 as it is always in 6 */
    Vector<Double> CalcTestStrainVector(TestStrain ts, bool reduced = false);

    /** return the von Mises matrix (stress^T * M * stress = von Mises Stress)
     * @param dim desired dimension. axis is ignored currently :( */
    const Matrix<double>& GetVonMisesMatrix(int dim);

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

    //! define surface integrators needed for this pde
    void DefineSurfaceIntegrators( ){};

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
    void DefineConcentratedElems();

    //! Create CoefFunction for preStressing
    PtrCoefFct CreatePreStressFct( bool isComplex, PtrParamNode stressNode);

    //! Returns a MechStress CoefFunction from a preceeding sequence step
    PtrCoefFct GetStressCoefFromSeqStep(UInt seqStep);

    //! Stores softening for each region
    std::map<RegionIdType, std::string> regionSoftening_;
    
    //! Stores Rayleigh damping definition for each region
    std::map<RegionIdType, RaylDampingData > regionRaylDamping_;
    
    //! Stores the linear stiffness for each region
    std::map<RegionIdType, PtrCoefFct > regionStiffness_;
    
    //! Dimension of stresses
    UInt stressDim_;
    
    //! Tensor type
    SubTensorType tensorType_;

    //! coefFunctzion for thermal strain
    shared_ptr<CoefFunctionMulti> thermalStrain_;

    //! coefFunctzion for thermal stress
    shared_ptr<CoefFunctionMulti> thermalStress_;

    StdVector<std::string> dofNames_;

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

