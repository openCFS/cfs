#ifndef FILE_SINGLEPDE
#define FILE_SINGLEPDE

#include "PDE/StdPDE.hh"

#include <list>
#include <map>

#include "Domain/Results/ResultInfo.hh"
#include "Domain/BCs.hh"
#include "Domain/Results/BaseResults.hh"

namespace CoupledField
{
  // forward class declaration
  class SpaceErrorEstimator;
  class DirectCoupledPDE;
  class Assemble;
  class BiLinearForm;
  class BaseBDBInt;
  class ResultFunctor;
  class BaseFieldFunctor;
  class LinearFormContext;
  class CoefFunction;
  class CoefFunctionSurf;
  class CoefFunctionMulti;
  class CoefFunctionFormBased;
  class IterCoupledPDE;
  class MathParser;
  
  //! Base class for all kinds of single field problems.

  class SinglePDE : public StdPDE {
  
  public:

    // friend declaration
    friend class DirectCoupledPDE;
    friend class BasePairCoupling;
    friend class IterCoupledPDE;

    
    typedef StdVector<shared_ptr<BaseResult> > ResultList;
    typedef std::map<shared_ptr<ResultInfo> , ResultList > ResultMap;

    //! Initialize PDEs (1st stage)
    
    //! In this phase of initialization, the following things are performed:
    //! - definition of regions, materials, damping
    //! - creation of Assemble class
    //! - definition of primary fefunctions, fespaces and time stepping
    //! - definition of post-processing results
    //! @param base pointer to InfoNode of this PDE */
    virtual void Init_Stage1( UInt sequenceStep, PtrParamNode base = PtrParamNode() ); 
    
    //! Initialize PDEs (2st stage)
    
    //! The second phase of initialization depends on all other SinglePDEs to
    //! have finished stage1. Within stage 2, the following things are performed
    //! - definition of boundary conditions
    //! - definition of integrators
    //! - definition of loads and RHS integrators 
    virtual void Init_Stage2();

    //! Initialize PDEs (3rd stage)
    
    //! In the third phase of initialization, the following tasks are performed:
    //! - Finalization of spaces and registration at OLAS
    //! - Initialization of time schemes
    //! - Finalization of postprocessing results
    //! - Incorporation of initial conditions
    void Init_Stage3( );

    // ---------------------- ***** --------------------------------

    //! Destructor
    virtual ~SinglePDE();
  
  
    // ======================================================
    // GET /SET  METHODS
    // ======================================================

    //! Return if PDE is formulated in updated Lagrangian coordinates
    bool IsUpdatedGeo() {
      return updatedGeo_;
    }
    
    //! check, if PDE has complex material parameters;
    //! has to be overwritten by specific PDE
    virtual bool IsMaterialComplex()  {
    	return isMaterialComplex_;
    }

    //! Set Direct coupling information
    virtual void SetDirectCoupling();

    //! Set iterative coupled PDE 
    
    //! Set the iterative coupled PDE this SinglePDE belongs to. This is
    //! needed in order to retrieve coupling quantities.
    void SetIterCoupledPDE( IterCoupledPDE* iterCplPde );

    IterCoupledPDE* GetIterCoupledPDE() {
      return iterCplPde_;
    }

    //! set boundary condition
    void SetBCs();

    //! set special PDE dependent boundary conditions
    virtual void SetSpecialBCs(){ return; }

    //! Update PDE due to updated step in multistep solution strategy
    virtual void UpdateToSolStrategy();
    
    /** Write general defines (BCs, loads, etc.) to info.xml.
     * Note, that only the current state is (over) written! */
    void WriteGeneralPDEdefines();

    ResultMap& GetResults() { return resultLists_; }
    
    /** Return the native solution type, MECH_DISPLACEMENT, ... */
    virtual SolutionType GetNativeSolutionType() const { EXCEPTION("not implemented"); }

    /**<p>This is part of ReadStoreResults(). If candidate is defined in the xml file
     * it is added to resultLists_.</p>
     * <p>This method is to be called by ReadStoreResults() for every element in
     * availResults_. Additionally an Optimization instance calls when there a
     * result element defines one of the solution types optResult_*in more detail
     * @param candidate normally an element of the (mathematical) set availResults_
     * @return true if in xml and added */
    bool CheckStoreResult(shared_ptr<ResultInfo> candidate);

    //! Define a field result
    void DefineFieldResult( PtrCoefFct coef, shared_ptr<ResultInfo> res );

    /** Shortcut for DefineFieldResult() */
    void DefineFieldResult(SolutionType solType, ResultInfo::EntryType entryType, ResultInfo::EntityUnknownType definedOn, const std::string& dofNames, bool fromOptimization);

    //! Manually set an entry to fieldCoefs_. This can be helpful when a result is defined
    //! in the xml file (via postproc results), where DefineFieldResults() is skipped.
    //! In such cases, the coefFunction is not available, which can be done later with this
    //! function
    void SetFieldCoef( PtrCoefFct coef, SolutionType resultType);

    //! Obtain coefficient function of given type
    PtrCoefFct GetCoefFct( SolutionType solType );
    
    /** return sub type. The string is stored internally any we need to convert. :(
     * @return if StdPDE::subType_ is not set we return NO_TENSOR  */
    virtual SubTensorType GetSubTensorType() const;

    //! Read single RHS excitation

    //! This method reads an xml element for a general RHS excitation and
    //! returns the entityList and CoefFunction.
    //! \param elemName Name of ParamNode within <bcsAndLoads> to be read
    //! \param compNames Names of the components (vector, tensor)
    //! \param type Type of CoefFunction to be read in (scalar, vector, tensor)
    //! \param isComplex Denotes  if a complex valued coef-function is to be
    //!                  generated
    //! \param entities Vector of entityLists of the boundary condition
    //! \param coef Vector of coefficients function for the values
    //! \param updateGeo Flag indicating, if coefficient function is defined
    //!                  on an updated geometry (e.g. due to iterative coupling).
    //! \param input default nullptr (will search in bcsAndLoad)
    //!                  if specified it will search elemName in input-Nodeset
    void ReadRhsExcitation( const std::string& elemName,
                            const StdVector<std::string>& compNames,
                            ResultInfo::EntryType type,
                            bool isComplex,
                            StdVector<shared_ptr<EntityList> >& entities,
                            StdVector<PtrCoefFct>& coef,
                            bool& updateGeo,
                            PtrParamNode input = PtrParamNode());

    void ReadRhsExcitation( const std::string& elemName,
                                const StdVector<std::string>& compNames,
                                ResultInfo::EntryType type,
                                StdVector<shared_ptr<EntityList> >& entities,
                                StdVector<PtrCoefFct>& coef,
                                bool& updateGeo,
                                PtrParamNode input = PtrParamNode());

    void ReadRhsExcitation( const std::string& elemName,
                                const StdVector<std::string>& compNames,
                                ResultInfo::EntryType type,
                                bool isComplex,
                                StdVector<shared_ptr<EntityList> >& entities,
                                StdVector<PtrCoefFct>& coef,
                                bool& updateGeo,
                                StdVector<std::string>& volumeRegions);
    void ReadRhsExcitation( const std::string& elemName,
                                    const StdVector<std::string>& compNames,
                                    ResultInfo::EntryType type,
                                    StdVector<shared_ptr<EntityList> >& entities,
                                    StdVector<PtrCoefFct>& coef,
                                    bool& updateGeo,
                                    StdVector<std::string>& volumeRegions);


    //! Read general external field information from given xml node
    //! The node has to contain either a values tag, a number of comp tags or
    //! a grid node
    //! \param[in] list EntityList the Field should be applied to
    //! \param[in] valueNode The xml node of the user parameters
    //! \param[in] compNames Names of the components (vector, tensor)
    //! \param[in] type Type of CoefFunction to be read in (scalar, vector, tensor)
    //! \param[in] isComplex Indicates if we need to account for Complex results
    //! \param[out] coef The generated coefficient function
    //! \param[out] definedDofs Set containing all defined dofs in case of a
    //!             vector-valued quantity.
    //! \param[out] updateGeo Flag indicating, if coefficient function is defined
    //!                  on an updated geometry (e.g. due to iterative coupling).
    //! \param[out] harm Only for MultiharmonicCase: specifies in which harmonic what value
    void ReadUserFieldValues( shared_ptr<EntityList> list,
                              PtrParamNode valueNode,
                              const StdVector<std::string>& compNames,
                              ResultInfo::EntryType type,
                              bool isComplex,
                              PtrCoefFct & coef,
                              std::set<UInt>& definedDofs,
                              bool& updateGeo,
                              PtrCoefFct & harm);

    //! as ReadUserFieldValues but determine isComplex from xml-input
    void ReadUserFieldValues( shared_ptr<EntityList> list,
                                  PtrParamNode valueNode,
                                  const StdVector<std::string>& compNames,
                                  ResultInfo::EntryType type,
                                  PtrCoefFct & coef,
                                  std::set<UInt>& definedDofs,
                                  bool& updateGeo);

    void ReadUserFieldValues( shared_ptr<EntityList> list,
                              PtrParamNode valueNode,
                              const StdVector<std::string>& compNames,
                              ResultInfo::EntryType type,
                              bool isComplex,
                              PtrCoefFct & coef,
                              std::set<UInt>& definedDofs,
                              bool& updateGeo);
    //! Read history result
    template<typename T>
    void ReadUserHistValues(  PtrParamNode valueNode,
        ResultInfo::EntryType type,
        Vector<T>& res,
        std::string regionName);


    /** Define all RHS linearforms for load / excitation
     * @param input for multiple load optimization we point to the multipleExcitation excitiation definition. Default is from bscAndLoads() */
    virtual void DefineRhsLoadIntegrators(PtrParamNode input) { }
    virtual void DefineRhsLoadIntegrators() { DefineRhsLoadIntegrators(PtrParamNode()); } // Only where we do optimization we use the parameter

    /** identify this pde for logging debug purpose */
    std::string ToString() const;

  protected:

    //! Constructor
    /*!
      \param aptgrid pointer to grid
    */
    SinglePDE( Grid *aptgrid, PtrParamNode, PtrParamNode infoNode,
               shared_ptr<SimState> simState,
               Domain* domain);

    //! private copy constructor
    SinglePDE & operator= (const StdPDE & myPDE) {
      EXCEPTION( "Not implemented" );

      // For compiler
      return *this;
      }

    //! Helper function for ReadRhsLoad ...
    void ReadEntities( const std::string& elemName,
                                         const StdVector<std::string>& compNames,
                                         ResultInfo::EntryType type,
                                         StdVector<shared_ptr<EntityList> >& entities,
                                         StdVector<PtrParamNode >& xmls,
                                         StdVector<PtrCoefFct >& coef,
                                         bool& updateGeo,
                                         PtrParamNode input);

    //! Helper function for ReadRhsLoad ...
    void ReadVolumeRegions( const std::string& elemName, StdVector<std::string>& volumeRegions);

    // ======================================================
    // INITIALIZATION METHODS
    // ======================================================

    //! Define primary results
    
    //! Initialize the primary results, i.e. the results corresponding to the
    //! primary variables and their unknowns.
    virtual void DefinePrimaryResults() { };
    
    //! Define post-processing results
    
    //! This method defines the post-processing results, which are computed
    //! mostly using the bilinearform of the main problem.
    //! \note This method gets called very early in the initialization 
    //! process, so that other PDEs can access already CoefFunctions using
    //! the GetCoefFct() method. The finalization is performed in the 
    //! method FinalizePostProcResults()
    virtual void DefinePostProcResults() {};
    
    
    //! Finalize post-processing results
     
    //! This method finalizes the setup of postprocessing results, e.g. 
    //! assigning bdb-forms to the calculation routines etc. This method 
    //! gets called at a very late stage.
    virtual void FinalizePostProcResults();

    //! Obtain information on desired output quantities from parameter file
    //! This method is used to query the parameter handling object for the
    //! desired output quantities and translate their literal description into
    //! the internal format by setting the corresponding class attributes.
    void ReadStoreResults();

    //! This function updates coefFunctions from previous definitions
    //! We might have defined some generic results which depend on the feFunction.
    //! Since they get initialized in stage 2, we had to pass the feFunction without
    //! knowing the actual coefFunction. We collected these BCs and update their
    //! coefFunction know accordingly
    void UpdateCoefFuncsForPostProcResults();

    //! define all (bilinearform) integrators needed for this pde
    virtual void DefineIntegrators( )=0;
    
    //! define surface integrators needed for this pde
    virtual void DefineSurfaceIntegrators( )=0;


    //! Read material depenecy information

    //! This method reads an xml element for a general material dependency and
    //! returns the CoefFunction.
    //! \param elemName Name of ParamNode within <matDependencyList> to be read
    //! \param compNames Names of the components (vector, tensor)
    //! \param type Type of CoefFunction to be read in (scalar, vector, tensor)
    //! \param isComplex Denotes  if a complex valued coef-function is to be
    //!                  generated
    //! \param entities Vector of entityLists
    //! \param coef Vector of coefficients function for the values
    //! \param updateGeo Flag indicating, if coefficient function is defined
    //!                  on an updated geometry (e.g. due to iterative coupling).
    void ReadMaterialDependency( const std::string& elemName,
                            const StdVector<std::string>& compNames,
                            ResultInfo::EntryType type,
                            bool isComplex,
                            shared_ptr<EntityList>& entity,
                            PtrCoefFct& coef,
                            bool& updateGeo );

    //! Read results information for interpolation of continuous fields
    virtual void ReadSensorArrayResults();
    
    // =======================================================================
    //   INTERPOLATION OF FIELD VARIABLES
    // =======================================================================
    
    //! Helper struct for interpolating field variables at arbitrary points
    struct FieldAtPoints {

      //! Physical Quantity
      shared_ptr<ResultInfo> resultInfo;
     
      //! Filename where points get written to
      std::string fileName;

      //! Format output file as CSV (comma separated values)
      bool csv;

      //! Delimiter for CSV fields
      char delim;

      //! Pointer to coordinate system
      CoordSystem * coordSys;
      
      //! Vector with elements 
      StdVector<const Elem*> elems;
      
      //! Vector with local points
      StdVector<LocPoint> locPoints;

      //! Flux values for each point
      SingleVector* field;
    };
    
    //! List of fields to be interpolated
    StdVector<FieldAtPoints> sensors_;
    
    //! read damping information
    virtual void ReadDampingInformation( ){
    };
    
    //! read material data
    virtual void ReadMaterialData();

    //! Read initial conditions / values 
    void ReadInitialConditions();
    
    //! read from config-file info about BCs
    void ReadBCs();

    //! Read periodic BC
    void ReadPeriodicBC(PtrParamNode prNode);
    
    //! overloaded version of ReadBCs for special
    //! boundary conditions in derived classes
    virtual void ReadSpecialBCs(){}
    
    //! write results in file
    void WriteResultsInFile( const UInt kstep, 
                             const Double actTimeFreq );
    
    //! Initialize NonLinearities
    virtual void InitNonLin();

    //! Initialize material dependencies
    virtual void InitMaterialDependencies();

    //! Define the time FeFunctions for this PDE according to the
    //! definition in the XML file
    virtual void DefineFeFunctions();

    //! Create FeSpaces according to formulation
    virtual std::map<SolutionType, shared_ptr<FeSpace> > 
    CreateFeSpaces( const std::string&  formulation,
                    PtrParamNode infoNode ) = 0;

    //@{
    
  public:
    //! Struct that collects the Rayleigh alpha and beta coefficients
    struct RaylDampingData {
      //! Damping parameters used for MASS and STIFFNESS integrator, respectively.
      //! The parameters are stored as strings that contain mathParser expressions.
      std::string alpha, beta;
    };

  protected:
    
     
    //@}

    // ======================================================
    // DATA SECTION
    // ======================================================

    // -----------------------------------------------------------------------
    // Storing information
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes connected to storing information
    
    
    //! Define result based on the time derivative of the main results
    void DefineTimeDerivResult( SolutionType derivSolType,
                                UInt timeDerivOrder,
                                SolutionType primSolType );
    
    //! Map containing the result types and the results
    ResultMap resultLists_;

    //@}
    
    // -----------------------------------------------------------------------
    // Miscellaneous parameters
    // -----------------------------------------------------------------------

    //@{
    //! \name Miscellaneous parameters

    //! flag for direct coupling
    bool isDirectCoupled_;

    //! flag indicating if Init() was already called
    bool isInitialized_;

    //! pointer to iterative coupled PDE
    IterCoupledPDE* iterCplPde_;
    
    //! Flag, if PDE used updated geometry (updated Lagrangian formulation)
    bool updatedGeo_;
    
    //! flag indicating that material parametters are complex
    bool isMaterialComplex_;

    //! Map for storing the primary BDB integrators of the problem
    
    //! This map stores the primary BDB integrators, which can be used for 
    //! calculating spatial derivatives, fluxes and energy.
    //! Here we use a multimap in order to enable storing multiple bdbInts (see e.g. LinFlowPDE)
    std::multimap<RegionIdType, BaseBDBInt*> bdbInts_;
    
    //! Map for storing the auxiliary primary BDB integrators of the problem

    //! This map stores the auxiliary primary BDB integrators, which can be used for
    //! calculating spatial derivatives, fluxes and energy.
    //! This is necessary because we can have a SinglePDE with more than
    //! one unknown (FeFunction) and also need the primary BOperator for that
    //! secondary FeFunction
    std::map<RegionIdType, BaseBDBInt*> bdbIntsAux1_;

    //! Map for storing the primary mass integrator of the problem
    
    //! This map stores the primary MASS integrators, which can be used for 
    //! calculating spatial derivatives, fluxes and energy.
    std::map<RegionIdType, BaseBDBInt*> massInts_;

    //! true, if analysistype is multiharmonic
    bool isMultHarm_;

    // -----------------------------------------------------------------------
    //  Result Handling
    // -----------------------------------------------------------------------
    //@{ \name Data for Result Handling
    
    //! Store field coefficient functions

    //! This map contains all coefficient functions, which calculate a "field"
    //! result, i.e. a spatially varying result (in contrast to an integrated 
    //! result like e.g. energy) 
    std::map<SolutionType, PtrCoefFct > fieldCoefs_;

    //! Vector of coefFunctions which need to be updated since they have been
    //! initialized with a dummy function
    StdVector<PtrCoefFct> coefsToUpdate_;
    
    //! Map for storing "material" parameters as coefficient functions
    
    //! This map holds the coefficient functions for the different material
    //! parameters and PML parameters. As material parameters are typically
    //! defined per region, we store them in a CoefFunctionMulti object
    std::map<SolutionType, shared_ptr<CoefFunctionMulti> > matCoefs_;
    
    //! Map storing functors for calculating general results
    
    //! This map stores the result functors for non-spatially varying results,
    //! e.g. energy, total force etc.
    std::map<SolutionType, shared_ptr<ResultFunctor> > resultFunctors_;
    
    //! stores the functors to field average results
    std::map<SolutionType, shared_ptr<ResultFunctor> > fieldAverageFunctors_;

    //! Store bilinarform-based coefficient function for stiffness integrator
    
    //! In this set we store all coefficient functions, which compute by the
    //! help of the primary stiffness matrix (e.g. flux values, energy values).
    //! In the method SinglePde::FinalizePostProcResults() every CoefFunction in
    //! this map gets related to the stiffness integrator on each region.
    std::set<shared_ptr<CoefFunctionFormBased> > stiffFormCoefs_;
    
    std::set<shared_ptr<CoefFunctionFormBased> > stiffFormCoefsAux1_;


    //! Store bilinarform-based coefficient function for mass integrator
    
    //! In this set we store all coefficient functions, which compute by the
    //! help of the mass matrix (e.g. kinetic energy, eddy current density).
    //! In the method SinglePde::FinalizePostProcResults() every CoefFunction in
    //! this map gets related to the mass integrator on each region.
    std::set<shared_ptr<CoefFunctionFormBased> > massFormCoefs_;
    
    //! Store result functors related to stiffness integrator
    std::set<shared_ptr<ResultFunctor> > stiffFormFunctors_;

    std::set<shared_ptr<ResultFunctor> > stiffFormFunctorsAux1_;

    //! Store bilinarform-based coefficient function for mass integrator
    std::set<shared_ptr<ResultFunctor> > massFormFunctors_;
    
    //! Store volume coefficient functions for surface CoefFunctions
    
    //! This map stores fore every coefficient function on a surface (e.g. surface
    //! charge density, normal velocity) the related volume CoefFunction,
    //! which has to be set for each region the surface is neighboring to.
    //! This is performed in the method  SinglePde::FinalizePostProcResults().
    std::map<shared_ptr<CoefFunctionSurf>, PtrCoefFct > surfCoefFcts_;

    std::map<shared_ptr<CoefFunctionSurf>, PtrCoefFct > surfCoefFctsAux1_;

    
    //! Map containing the input states and the related domains
    std::map<shared_ptr<SimState>, Domain* > inputs_;
    //@}

    // -----------------------------------------------------------------------
    //  Non-conforming interfaces
    // -----------------------------------------------------------------------
    //@{ \name Functions and data for non-conforming interfaces
    
  protected:
    
    //! Reads ncInterfaces defined in the XML file
    virtual void ReadNcInterfaces();
    
    //! Creates FeSpaces for additional unknowns used by ncInterfaces
    
    //! Default behavior in SinglePDE is to create the same FeSpace as for
    //! the first primary result. So this function must be overridden in PDEs
    //! with more than one primary result.
    virtual void CreateNcFeSpaces(std::map<SolutionType, shared_ptr<FeSpace> >
                                      &spaces,
                                  const std::string &formulation,
                                  PtrParamNode infoNode);
    
    //! Defines auxiliary results for non-conforming interfaces
    
    //! Default behavior in SinglePDE is to add the Lagrange multiplier as
    //! a primary result. This can be overridden in a derived class.
    virtual void DefineNcAuxResults();
    
    //! Defines the integrators needed for ncInterfaces
    
    //! If a derived class does not support ncInterfaces, it should implement
    //! an exception and not simply ignore this.
    virtual void DefineNcIntegrators() = 0;
    
    //! Defines integrators for Mortar coupling of an unknown on one specific
    //! interface.
    template<UInt DIM, UInt D_DOF>
    void DefineMortarCoupling(SolutionType solType, NcInterfaceInfo &iface);

    //! Defines integrators for Nitsche coupling of an unknown on one specific
    //! interface.
    template<UInt DIM, UInt D_DOF>
    void DefineNitscheCoupling(SolutionType solType, NcInterfaceInfo &iface,
                                shared_ptr<CoefFunctionMulti> additionalCoef = NULL);
                               
    //! When the mortarInterface (Nitsche or Mortar) has ALESystem enabled, we assemble additional convective terms
    //! that counter-rotate the field against the grid motion, resulting in a stationary material behavior.
    //! The system is defined similar to the convective terms in the AcouPDE, which is based on 
    //! "Kaltenbacher, Hüppe, 2018: Advanced Finite Element Formulation for the Convective Wave Equation"
    template<UInt DIM, UInt D_DOF>
    void DefineEulerianSystem(SolutionType solType, NcInterfaceInfo &iface);

    //! Vector containing all ncInterfaces for this PDE
    StdVector< NcInterfaceInfo > ncInterfaces_;
    
    //@}
  private:

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class SinglePDE
  //! 
  //! \purpose 
  //! This class serves as base class for all single field problems, 
  //! like electrostatic,  acoustic, mechanic and others.
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
