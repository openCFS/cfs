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
  class set;

  //! Class for mechanic equation (no adaptivity)
  class MechPDE: public SinglePDE
  {

  public:

    /** constants for test-strains, used for homogenization. We depend on the int values! */
    typedef enum { NOT_SET=-1, X=0, Y=1, Z=2, YZ=3, XZ=4, XY=5 } TestStrain;

    static Enum<TestStrain> testStrain;

    /** Constructor. here we read integration parameters
     * @param aGrid pointer to grid */
    MechPDE( Grid *aGrid, PtrParamNode paramNode );

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
    void DefinePressureIntegrators(StdVector<shared_ptr<EntityList> >& pressSurf, StdVector<std::string>& pressVals, StdVector<std::string>& pressPhase, StdVector<LinearFormContext*>* linForms = NULL);
    
    /** Add the integrators for the test strains for homogenization to the linear forms, similar as in multiple load case;
     * called from Excitation::ReadLoads or Excitation::SetHomogenizationTestStrains() (optimization)
     * @param test is an enum
     * @param linForms set to append linear Forms to, if NULL use assemble_ */
    void DefineTestStrainIntegrator(const TestStrain test, StdVector<LinearFormContext*>* linForms = NULL);

    /** small helper which translates the test strain code on a vector of size 3/6 with one entry 1.0, the other zero
     * @param reduced in 2d case the result size is 3 otherwise 6 as it is always in 6 */
    Vector<Double> CalcTestStrainVector(TestStrain ts, bool reduced = false);
    
    /** export of methods to generate Linear Forms used in multiload-cases by optimization 
     * @param regionLoads as returned from ReadRegionLoadsFromXML
     * @param linForms set to append linear Forms to, if NULL use assemble_ */
    void DefineRegionLoadIntegrators(std::map<RegionIdType, RegionLoad>& regionLoads, StdVector<LinearFormContext*>* linForms = NULL);

    /** Does the actual reading of pressure loads, also called from optimization 
     * @param bcNode paramnode that has "pressure" nodes as children 
     * @param pressSurf StdVector containing the information
     * @param pressVals StdVector containing the information
     * @param pressPhase StdVector containing the information */
    void ReadPressureLoadsFromXML(PtrParamNode bcNode, StdVector<shared_ptr<EntityList> >& pressSurf, StdVector<std::string>& pressVals, StdVector<std::string>& pressPhase);
    
    /** add the integrators for the polarization matrix to the linear forms, similar as in multiple load case;
      * called from Excitation::SetPolarizationMatrixRHS
      * @param vals contains the values from the xml test strains
      * @param linForms set to append linear Forms to, if NULL use assemble_ */
    void DefinePolarizationMatrixIntegrators(const Vector<Double> &vals,
        StdVector<LinearFormContext*> *linForms, const int num);

    /** @see virtual SinglePDE::GetNativeSolutionType() */
    SolutionType GetNativeSolutionType() const { return MECH_DISPLACEMENT; }

    /** @see virtual SinglePDE::GetNativeDOF() */
    virtual UInt GetNativeDOF() const { return dim_; }

    /** return the von Mises matrix (stress^T * M * stress = von Mises Stress)
     * @param dim desired dimension. axis is ignored currently :( */
    const Matrix<double>& GetVonMisesMatrix(int dim);

    /** Scalar version of the strains/stresses. Called from Optimization/ErsatzMaterial.
     * @param ss wither MECH_STRESS or MECH_STRAIN */
    template <class TYPE>
    void CalcVonMises(shared_ptr<BaseResult> res, SolutionType ss);

  protected:

    // ======================================================
    // POSTPROCESSING METHODS
    // ======================================================

    /** map nodal results to element results */
    template <class TYPE>
    void CalcLumpedDisplacement(shared_ptr<BaseResult> result);

    void CalcHomogenizedTensor(shared_ptr<BaseResult> base_result);

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
    Double aitkenOmega_, aitkenOmegaPrevIter_, fixedOmega_;
    Double aitkenMu_;
    Vector<Double> actDelta_, oldDelta_;
    bool firstTime_;
    Vector<Double> gSolOld_;

    Vector<Double> deltaTildeDisplPrevIter_;

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

    /** checks and processes the xml file for test strain, to do kind of manual (mathematical) homogenization */
    void ReadTestStrains();

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
    
    //! Stores Rayleigh damping definition for each rgion
    std::map<RegionIdType, RaylDampingData > regionRaylDamping_;

    //! Flag indicating use of penalty dof for plate formulation
    bool usePlatePenaltyDof_;

    /** static matrix which allows the calculation of the von Mises Stress from the stresses in 2D.
     * See Kocvara and Stingl; 2007 */
    Matrix<double> vonMisesMatrix_2d_;
    Matrix<double> vonMisesMatrix_3d_;

    void calcAitkenOmega(const Vector<Double>& displTilde, const Vector<Double>& displPrevIter);

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

