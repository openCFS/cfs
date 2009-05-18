// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_NEWBASEMECHPDE
#define FILE_NEWBASEMECHPDE

#include <map>

#include "SinglePDE.hh"

namespace CoupledField
{

  class BaseForm;
  class LinearFormContext;

  //! Class for mechanic equation (no adaptivity)
  class MechPDE: public SinglePDE
  {

  public:

    //!  Constructor. here we read integration parameters
    /*!
      \param aGrid pointer to grid
    */
    MechPDE( Grid *aGrid, ParamNode* paramNode );

    //!  Deconstructor
    virtual ~MechPDE();

    //! read in damping information, see SinglePDE.cc  and SinglePDE.hh
    void ReadDampingInformation();

    //! Initialize NonLinearities
    void InitNonLin();

    //! define all (bilinearform) integrators needed for this pde
    void DefineIntegrators( );

    //! define the SoltionStep-Driver
    void DefineSolveStep();

    //! Read special boundary conditions
    void ReadSpecialBCs();

    //! Read special results definition
    void ReadSpecialResults();

    /// returns a stiffness integrator appropriate to the actual problem (e.g.3D)
    BaseForm * GetStiffIntegrator(BaseMaterial* actSDMat,
                                  RegionIdType regionId,
                                  bool reducedInt=false);

    /** Returns the protected damping list - don't alter! */
    const std::map<RegionIdType,DampingType>& GetDampingList() { return dampingList_; }

    // ======================================================
    // COUPLING SECTION
    // ======================================================

    //! initalize PDE coupling
    void InitCoupling(PDECoupling * Coupling);

    //! calculate coupling terms
    void CalcOutputCoupling();

    //! returns if PDE can compute the quantity
    bool HasOutput(SolutionType output);

    /// setup source term
    void SetupRHS( );

    //! Turn on thermo-elastic coupling
    void SetHeatCoupling() {
      isHeatCoupled_ = true;
    }

    // ======================================================
    // POSTPROC SECTION
    // ======================================================

    //! Calculate result for given result class
    void CalcResults( shared_ptr<BaseResult> results );
    
    // ======================================================
    // INTERFACE USED BY OPTIMIZATION, MULTI-LOAD
    // ======================================================
    
    /** export of methods to generate Linear Forms used in multiload-cases by optimization
     * @param pressSurf as returned from ReadPressureLoadsFromXML
     * @param pressVals as returned from ReadPressureLoadsFromXML
     * @param pressPhase as returned from ReadPressureLoadsFromXML
     * @param linForms set to append linear Forms to, if NULL use assemble_ */
    void DefinePressureIntegrators(StdVector<shared_ptr<EntityList> >& pressSurf, StdVector<std::string>& pressVals, StdVector<std::string>& pressPhase, std::set<LinearFormContext*>* linForms = NULL);
    
    /** export of methods to generate Linear Forms used in multiload-cases by optimization 
     * @param regionLoads as returned from ReadRegionLoadsFromXML
     * @param linForms set to append linear Forms to, if NULL use assemble_ */
    void DefineRegionLoadIntegrators(std::map<RegionIdType, RegionLoad>& regionLoads, std::set<LinearFormContext*>* linForms = NULL);

    /** Does the actual reading of pressure loads, also called from optimization 
     * @param bcNode paramnode that has "pressure" nodes as children 
     * @param pressSurf StdVector containing the information
     * @param pressVals StdVector containing the information
     * @param pressPhase StdVector containing the information */
    void ReadPressureLoadsFromXML(ParamNode* bcNode, StdVector<shared_ptr<EntityList> >& pressSurf, StdVector<std::string>& pressVals, StdVector<std::string>& pressPhase);

  protected:

    // ======================================================
    // POSTPROCESSING METHODS
    // ======================================================

    //computes mechanical deformation energy
    template <class TYPE>
    void CalcEnergy( shared_ptr<BaseResult> vals );

    //computes mechanical stresses
    template <class TYPE>
    void CalcStresses(  shared_ptr<BaseResult> vals );

    //computes mechanical strains
    template <class TYPE>
    void CalcStrains(  shared_ptr<BaseResult> vals );

    //! compute volume above a deformed surface
    template <class TYPE>
    void ComputeVolDefSurf( shared_ptr<BaseResult> vals );

    //! computes averaged volume of an deformed element
    template<class TYPE>
    TYPE ComputeVolElem(BaseFE * ptSurfEl, Matrix<Double>& SurfCoord,
                        Vector<TYPE> disp);
    void ExtractNodeOffset(shared_ptr<BaseResult> result);


    //! Nodestoresol for RHS
    BaseNodeStoreSol * rhs_;


    // ========================
    // set solution information
    // ========================
	//! flag for thermo-elastic coupling
    bool isHeatCoupled_;

    //! vector containing regionIds of non-conforming interfaces
    StdVector<RegionIdType> ncIFaces_;


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
    //! internal attribute of the MechPDE class that is set, if the keyword
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

    Double displFac_;
    bool useAitken_, FSI_;
    Double aitkenOmega_, fixedOmega_;
    Double aitkenMu_;
    Vector<Double> actDelta_, oldDelta_;
    bool firstTime_;
    Vector<Double> gSolOld_;

  private:

    // calc rhs coupling to acoustic pde
    // void CalcAcousticCouplingRHS(StdVector<Elem*> * couplingElems,
    // Vector<Double>& forceOnElem);

    /// calc rhs coupling to acoustic pde
    void CalcAcousticCouplingRHS( StdVector<Elem*> * couplingElems,
                                  StdVector<BaseMaterial*> & materials,
                                  StdVector<UInt>& couplingNodes,
                                  Vector<Double> & forceOnElem,
                                  UInt couplingdof );



    /// does a line search and returns the optimal residual norm
    Double LineSearch(Vector<Double>& solIncrement, Vector<Double>& actSol,
                      Double& etaLineSearch, bool trans=false);


    /// Write nonlin iteration norms to the cla-file
    void WriteClaNlNorms( const UInt iterationCounter,
                          const Double residualL2Norm,
                          const Double extForcesL2Norm, const Double residualErr,
                          const Double solIncrL2Norm, const Double actSolL2Norm,
                          const Double incrementalErr );


    //! read in the domains with prestressing
    void ReadPreStressing();

    //! read in surface stress boundary conditions
    void ReadSurfStress();

    //! read pressure loads
    void ReadPressureLoads();
    
    //! read in softening types
    void ReadSoftening();

    //! define available result types
    void DefineAvailResults();

    struct SurfStress {

      //! Name of surface elements
      RegionIdType surface;

      //! Name of neighbouring region
      RegionIdType region;

      //! Vector of load
      Vector<Double> stress;

      //! Phase value
      std::string phase;
    };


    //! List of surface stresses
    std::map<RegionIdType, SurfStress> surfStresses_;

    //! Read information about springs and define related integrators
    void DefineSprings();

    /// stores an algsys_ vector into a StdVector and returns that L2-norm
    void StoreAlgsysToVec(Vector<Double>& vec, Double * pt);

    /// returns that L2-norm of an algsys vector
    Double AlgsysL2Norm(Double * pt);

    /// flag for reduced Integration for each subdomain
    StdVector<std::string> reducedIntegration_;

    /// returns the solution matrix belonging to all nodes of the actual element
    void GetSolOfElement( Matrix<Double>& elDisp, StdVector<UInt>& connect_PDE);

    //! Number of dimension for stresses
    UInt stressDim_;

    //! Flag indicating the use of fractional damping
    //bool fracDamping_;

    //! surface of pressure loads
    StdVector<shared_ptr<EntityList> > pressSurf_;

    //! values of the pressure loads
    StdVector<std::string>  pressVals_;

    //! phase of the pressure loads
    StdVector<std::string>  pressPhase_;

    //! list of prestressing types
    std::map<RegionIdType,std::string> preStressList_;

    //! prestress-values: numSubdoms x 3
    std::map< RegionIdType, Vector<Double> > preStressVal_;

     //@{ \name Attributes related to post-processing

     //! Contains mechanic velocity
    NodeStoreSol<Double> solDeriv1_;

    //! Contains mechanic acceleration
    NodeStoreSol<Double> solDeriv2_;

    //! Contains mechanic displacement of last time step
    NodeStoreSol<Double> sol_tn_1_;

    //! Contains the directions for which the deformed volume is computed
    std::map<shared_ptr<EntityList>,std::string> volAboveDefSurfDir_;

    //! Stores softening for each region
    std::map<RegionIdType, std::string> regionSoftening_;

    //! Flag indicating use of penalty dof for plate formulation
    bool usePlatePenaltyDof_;

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

