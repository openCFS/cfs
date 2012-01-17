// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NEWBASEFLUIDMECHPDE
#define FILE_NEWBASEFLUIDMECHPDE

#include <map>
#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/defs.hh"
#include "General/environment.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "SinglePDE.hh"
#include "Domain/bcs.hh"
#include "Utils/StdVector.hh"
#include "Utils/nodestoresol.hh"

namespace CoupledField {
class BaseFE;
class BaseMaterial;
class BaseNodeStoreSol;
class BaseResult;
class EntityIterator;
class EntityList;
class Grid;
class PDECoupling;
struct Elem;
}  // namespace CoupledField

namespace CoupledField
{

  class BaseForm;

  //! Class for fluid mechanic equation (no adaptivity)
  class FluidMechPDE: public SinglePDE
  {

  public:

    //!  Constructor. here we read integration parameters
    /*!
      \param aGrid pointer to grid
    */
    FluidMechPDE( Grid *aGrid, PtrParamNode paramNode );

    //!  Deconstructor
    virtual ~FluidMechPDE();

    //! Read special boundary conditions
    void ReadSpecialBCs();

    //! Read special results definition
    void ReadSpecialResults();

    //! read in damping information, see SinglePDE.cc  and SinglePDE.hh
    void ReadDampingInformation();

    //! Initialize NonLinearities
    void InitNonLin();

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators( );

    //! define the SoltionStep-Driver
    void DefineSolveStep();

    //! Initialize staqbilization parameters
    void InitStabParams();
    void PrintStabParams();

    /// returns a stiffness integrator appropriate to the actual problem (e.g.3D)
    BaseForm * GetStiffIntegrator(BaseMaterial* actSDMat,
                                  RegionIdType regionId,
                                  bool reducedInt=false);

    /// returns stiffness integrator appropriate to the actual problem (e.g.3D)
    BaseForm * GetStiffIntegrator_UV(Double density,Double kinematicViscosity);
    BaseForm * GetStiffIntegrator_PQ(Double density,Double kinematicViscosity);
    BaseForm * GetStiffIntegrator_UQ(Double density,Double kinematicViscosity);
    BaseForm * GetStiffIntegrator_PV(Double density,Double kinematicViscosity);

    /// returns mass integrator appropriate to the actual problem (e.g.3D)
    BaseForm * GetMassIntegrator_UV(Double density,Double kinematicViscosity);
    BaseForm * GetMassIntegrator_UQ(Double density,Double kinematicViscosity);
    // these two are necessary for the newton method
    BaseForm * GetAuxIntegrator_UV(Double density,Double kinematicViscosity);
    BaseForm * GetAuxIntegrator_UQ(Double density,Double kinematicViscosity);

    //! Calculate the acoustic sources due to the flow
    void AcouSourceCalc();

    // ======================================================
    // COUPLING SECTION
    // ======================================================

    //! initalize PDE coupling
    void InitCoupling(PDECoupling * Coupling);

    //! Fill in input coupling terms
    void CalcInputCoupling();

    //! calculate coupling terms
    void CalcOutputCoupling();


    //! returns if PDE can compute the quantity
    bool HasOutput(SolutionType output);



    // ======================================================
    // POSTPROC SECTION
    // ======================================================

    //! Calculate result for given result class
    void CalcResults( shared_ptr<BaseResult> results );

  protected:

    // ======================================================
    // POSTPROCESSING METHODS
    // ======================================================

    //computes mechanical deformation energy TODO
    template <class TYPE>
    void CalcEnergy( shared_ptr<BaseResult> vals );

    //computes mechanical stresses TODO
    template <class TYPE>
    void CalcStresses(  shared_ptr<BaseResult> vals );

    //computes fluid mechanical strain rates
    template <class TYPE>
    void CalcStrainRates(  shared_ptr<BaseResult> vals );

    //! compute volume above a deformed surface TODO
    template <class TYPE>
    void ComputeVolDefSurf( shared_ptr<BaseResult> vals );

    //! computes averaged volume of an deformed element
    template<class TYPE>
    TYPE ComputeVolElem(BaseFE * ptSurfEl, Matrix<Double>& SurfCoord,
                        Vector<TYPE> disp);


    //! Nodestoresol for RHS
    BaseNodeStoreSol * rhs_;



    // ======================================================
    // SCRIPTING SECTION
    // ======================================================
    //@{
    //! \name Scripting Methods

    //! Register scriptable functions
    void RegisterFunctions();
    //@}

    //! Obtain information on desired output quantities from parameter file

    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    //! The output quantities currently supported by the mechanics PDE are
    //! given in the following table. Here 'Keyword' and 'Result Type' refer
    //! to the XML parameter file, while 'Class Attribute' refers to the
    //! internal attribute of the FluidMechPDE class that is set, if the keyword
    //! is specified.\n\n
    //! <table border="1">
    //!   <tr>
    //!     <td><b>Keyword</b></td>
    //!     <td><b>Result Type</b></td>
    //!     <td><b>Class Attribute</b></td>
    //!   </tr>
    //!   <tr>
    //!     <td>displacement</td>
    //!     <td>nodeResults</td>
    //!     <td>savesol_</td>
    //!   </tr>
    //!   <tr>
    //!     <td>velocity</td>
    //!     <td>nodeResults</td>
    //!     <td>savederiv_</td>
    //!   </tr>
    //!   <tr>
    //!     <td>acceleration</td>
    //!     <td>nodeResults</td>
    //!     <td>savederiv2_</td>
    //!   </tr>
    //! </table>
    void ReadStoreResults();

    //! Init the time stepping
    void InitTimeStepping();

    //! Contains grid velocity
    NodeStoreSol<Double> * gridSol_;
    NodeStoreSol<Double> acou_src_;

    Double forceFac_;

  private:

    /// returns the vector of the fluid mechanical pressure solution
    /// belonging to all nodes of the actual element
    void GetPresSolVecOfElement(Vector<Double>& sol, StdVector<UInt>& connect_PDE);

    /// returns the vector of the fluid mechanical velocities solution
    /// belonging to all nodes of the actual element
    template<class TYPE>
    void GetVeloSolMatOfElement(Matrix<TYPE>& sol, StdVector<UInt>& connect_PDE);

    /// returns the vector of the time derivative of fluid mechanical velocities solution
    /// belonging to all nodes of the actual element
    void GetVeloDeriv1SolVecOfElement(Vector<Double>& sol, const EntityIterator& it);


    //! calculate the vector of coupling forces to the mechanical PDE
    void CalcMechCouplingRHS( StdVector<Elem*> * couplingElems,
                              StdVector<UInt> & couplingNodes,
                              Vector<Double>& elemCouplingSols,
                              UInt couplingdof );

    //! calculate the vector of coupling RHS for the acoustic PDE
    //! due to the moving fluid-structure inerface
    void CalcAcouSurfSourceCouplingRHS( StdVector<Elem*> * couplingElems,
                                        StdVector<UInt> & couplingNodes,
                                        Vector<Double>& elemCouplingSols,
                                        UInt couplingdof );


    //! read pressure loads
    void ReadPressureLoads();


    //! define available result types
    void DefineAvailResults();

    /// stores an algsys_ vector into a StdVector and returns that L2-norm
    void StoreAlgsysToVec(Vector<Double>& vec, Double * pt);

    /// returns that L2-norm of an algsys vector
    Double AlgsysL2Norm(Double * pt);

    /// returns the solution matrix belonging to all nodes of the actual element
    void GetSolOfElement( Matrix<Double>& elDisp, StdVector<UInt>& connect_PDE);


    //! Number of dimension for stresses
    UInt stressDim_;

    //! surface of pressure loads
    StdVector<shared_ptr<EntityList> > pressSurf_;

    //! values of the pressure loads
    StdVector<std::string>  pressVals_;

    //! phase of the pressure loads
    StdVector<std::string>  pressPhase_;

    //! dynamics of pressure loads
    StdVector<std::string> pressFnc_;

    //! stabilization parameters
    Matrix<Double> stabilParams_;

    //! list of prestressing types
    std::map<RegionIdType,std::string> preStressList_;

    //! prestress-values: numSubdoms x 3
    std::map< RegionIdType, Vector<Double> > preStressVal_;

     //@{ \name Attributes related to post-processing

     //! Contains mechanic velocity
    NodeStoreSol<Double> solDeriv1_;

    //! Contains mechanic acceleration
    NodeStoreSol<Double> solDeriv2_;

    //! surface elements with absorbing boundary conditions
    StdVector<shared_ptr<EntityList> > absBCs_; 

    //! Contains the directions for which the deformed volume is computed
    std::map<shared_ptr<EntityList>,std::string> volAboveDefSurfDir_;

    std::string approxType_;           // string to specify approximation type (Lagrange/Legendre/Taylor-Hood/)
    std::string stabilizationType_;    // stabilization type (none/SUPG)
    // defines the non linear method. Needs to be set here to define forms.
    std::string nonLinMethod_;

    bool movingMesh_; // flag to indicate treatment of moving meshes (ALE)

    // ========================
    // coupling
    // ========================
    //! assigns each coupling element node the according Coupling Node number
    StdVector<StdVector<StdVector<UInt> > > elemNodeToCouplingNode_;


    bool saveAcouSrc_; //Calculate the acoustic sourceterms
    bool saveAcouSrcInFile_; //Shall the source_terms be written in a save file?
    bool coordUp_;

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class
  // =========================================================================

  //! \class FluidMechPDE
  //!
  //! \purpose
  //! This class defines the fluid mechanical field PDE (Navier-Stokes equations)
  //! and the according
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

