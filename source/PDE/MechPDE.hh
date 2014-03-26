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

  class BaseForm;
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
    void DefineRhsLoadIntegrators();
    
    //! define the SoltionStep-Driver
    void DefineSolveStep();

    //! Read special results definition
    void ReadSpecialResults();

    //! \copydoc SinglePDE::CreateFeSpaces
    virtual std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string&  formulation,
                    PtrParamNode infoNode );
    //! Returns a stiffness integrator appropriate to the actual problem (e.g. 3D)
    BaseBDBInt * GetStiffIntegrator( BaseMaterial* actSDMat,
                                     RegionIdType regionId,
                                     bool isComplex );
    
    //! Return strain operator 
    BaseBOperator * GetStrainOperator( bool isComplex, bool icModes );

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
  //!
  //! \collab
  //!
  //! \implement
  //!
  //! \status In use
  //!
  //! \unused
  //!
  //! \improve
  //!

#endif

} // end of namespace
#endif

