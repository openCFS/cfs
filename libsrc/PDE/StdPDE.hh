#ifndef FILE_STDPDE
#define FILE_STDPDE

#include "PDE/basePDE.hh"

#include "DataInOut/writeresults.hh"
#include "Driver/assemble.hh"
#include "PDE/timestepping.hh"
#include "Driver/stdSolveStep.hh"
#include "nodeEQN.hh"

namespace CoupledField {


  // forward class declarations
  class PDECoupling;
  
  //! Base class for all single-field and direct-coupled problems

  //! This class serves as base class for all single-field and 
  //! direct coupled problems. The idea is, that an iterative coupled
  //! pde can have one or more objects of StdPDE, which are then
  //! solved iteratively.
  
  class StdPDE : public BasePDE {
  
  public:

    // friend class declarations
    friend class StdSolveStep;
    friend class PDEMemento;
    friend class PDECoupling;

    //! Virtual destructor
    virtual ~StdPDE() {
      ENTER_FCN( "StdPDE::~StdPDE" );
    }

    // ------------------- ASK TOM TO CLEAN UP THIS SECTION -------------------

    //@{
    //! \todo The following Methods are used only durig parameter
    //!       identificationprocess! Ask Tom to clean this section up!
    MaterialData * getPDEMaterialData()
    {return materialData_;};
  
    BaseNodeStoreSol * getPDESolution()
    {return sol_;};
  
    BCs * getPDE_BCs()
    {return ptBCs_;};
  
    BaseSystem * getPDE_algsys()
    {return algsys_;};
  
    Integer getPDE_numElems()
    {return numElems_;};

    Integer getPDE_dofspernode()
    {return dofspernode_;};
  
    Integer getPDE_numPDENodes(){return numPDENodes_;};
  
    Integer getPDE_spaceDim(){return dim_;};
  
    NodeEQN * getPDE_eqnData(){return eqnData_;}; 
  
    Grid * getPDE_grid(){return ptgrid_;};  
  
    Assemble * getPDE_assemble(){return assemble_;}
  
    StdVector<std::string> getPDE_subdoms(){return subdoms_;}
  
    Vector<Complex> complexValuedCharge_;
  
    Vector<Complex> getPDE_complexValuedCharge(){return complexValuedCharge_;};

    void setBCs_id_phase_(Integer i, Double & phase){
      ENTER_FCN("SinlgePDE::setBCs_id_phase");
      if (bcs_id_phase_.GetSize() <= i)
        Error("no such index in Vector bcs_id_phase_",__FILE__,__LINE__);
      bcs_id_phase_[i]=phase;
    };
  
    void setPDE_actFrequency(Double & freq){
      ENTER_FCN("SinglePDE::setPDE_actFrequency");
      actFrequency_ = freq;
    };
    void setPDE_actFreqStep(Integer & fstep){
      ENTER_FCN("SinglePDE::setPDE_actFreqStep");
      actFreqStep_ = fstep;
    };

    piezoMaterialType piezoMaterialType_; 
  
    void setPDE_piezoMaterialType(piezoMaterialType & pMatType){
      piezoMaterialType_ = pMatType;};
  
    piezoMaterialType getPDE_piezoMaterialType(){return piezoMaterialType_;}
    //@}

    // ======================================================
    // ALGSYS SECTION (SOLVER, ...)
    // ======================================================
  
    //! define algebraic system 
    virtual void SetAlgSys();
  
    //! Create the matrices and Solver as well as Preconditioner
    virtual void CreateMatrices_Solver();
  
    //! deletes the algebraic system
    void DeleteAlgSys(int as_id)
    {if (algsys_)
      assemble_->DeleteAlgSys();
    };

    // ======================================================
    // COUPLING SECTION
    // ======================================================
  
    //! initalize PDE coupling (only done once)
    virtual void InitCoupling(PDECoupling * Coupling) = 0;

    //! reset coupling counters and data (done after each timestep)
    virtual void ResetCoupling() = 0;
  
    //! Fill in input coupling terms
    virtual void CalcInputCoupling() = 0;
  
  
    //! calculate coupling terms
    virtual void CalcOutputCoupling() = 0;
  


  // ======================================================
  // POSTPROCESSING METHODS
  // ======================================================
  
  //! compute volume above a deformed surface
  virtual void ComputeVolDefSurf(StdVector<std::string> surfRegions, 
				 StdVector<std::string> strDir);

  //!
  virtual Double ComputeVolElem(BaseFE * ptSurfEl, Matrix<Double>& SurfCoord, 
				Vector<Double> disp);
  
    // ======================================================
    // GETTER METHODS
    // ======================================================
    //! get the encapsulated state of the PDE
  
    //! returns the current state of the PDE (solution, derivative,
    //! coupling-objects) in an encapsulated object. This is needed to
    //! enable full MultiSequence simulation, where from one step to 
    //! another the solution, the derivative and perhaps coupling 
    //! values like geometry update have to be passed. 
    //! The PDEMemento object encapsulates this information. 
    //! Later on the information can be given back to the PDE
    //! with the method SetMemento();
    //! \param memento (output) Object where the current state gets saved
    virtual void GetMemento(PDEMemento & memento);
  
    //! set the encapsulated state of the PDE
  
    //! set the current state of this PDE (solution, derivative,
    //! coupling-objects) from an encapsulated object. This is needed to
    //! enable full MultiSequence simulation, where from one step to 
    //! another the solution, the derivative and perhaps coupling 
    //! values like geometry update have to be passed. 
    //! The PDEMemento object encapsulates this information. 
    //! With this method the previous stored information can be set
    //! to the current PDE.
    //! \param memento (input) Previously saved state of the PDE
    //! \param transFromTo (input) : standard or complexToReal 
    virtual void SetMemento(PDEMemento & memento, std::string transFromTo);

    //! Get number of time step
    virtual Integer GetTimeStepCounter()
    { return laststepcalc_; }

    //! return pointer to vector with subdomains, on which we calculate the PDE
    virtual StdVector<std::string> * getSDsPDE()
    { return &subdoms_;}
  
    //! returns if PDE can compute the quantity
    virtual Boolean HasOutput(SolutionType output)
    {
      Error("not implemented",__FILE__,__LINE__);
      return FALSE;
    }
  
    //! return pointer to vector with first derivative of solution
    virtual const Vector<Double> & getS1() const;
  
    //! return pointer to vector with first derivative of solution
    virtual const Vector<Double> & getS2() const;
  
    //! return size of solution
    virtual Integer getSize() const 
    {
      Error("Function getSize is not overloaded in this class",__FILE__,
            __LINE__);
      return 0;

    }

    //! Get coefficient for damping matrix in fractional damping model
    virtual Double GetFracDampMatrixCoeff(Integer actSD);


    //! computes the coordinates of an element including the delta
    //! \param connect (input) global node numbers of element
    //! \param ptCoord (output) coordinates of the element nodes
    //!                (nrNodes \f$\times\f$ spaceDim);
    //! \param level (input) index for multilevel hierarchy
    virtual void GetElemCoords(const StdVector< Integer > connect,
                               Matrix< Double > &coordMat,
                               const Integer level);
  
    //!
    virtual void ComputeRHS(const Double atime) {;};
  
    //! set boundary condition
    //! \param level             level of grid
    //! \param atimestep         time step of claculation
    virtual  void SetBCs(const Integer level, 
                         const Double atimestep) = 0;
  
  
    //! set time step
    //! \params dt Current time step
    virtual void SetTimeStep(const Double dt)
    {TS_alg_->Init(matrix_factor_, dt);}
  
  protected:
  
    //! Constructor
    /*!
      \param aptgrid pointer to grid
      \param aptBCs pointer to boundary condition object
      \param aInFile pointer to class FileType. input data.
      \param aOutFile  pointer to class WriteResults. output data.
      \param aTimeFunc pointer to class TimeFunc
    */
    StdPDE(Grid *aptgrid, BCs *aptBCs, FileType *aInFile,
           WriteResults *aOutFile, TimeFunc *aTimeFunc);
  
    //! private copy constructor
    StdPDE & operator= (const StdPDE & myPDE) {
    };
  
  
    //! return index of dof defined by keyword (e.g. 'ux')
    virtual Integer GetBCDof(const std::string keyword);

    //! stores an algsys_ vector into a StdVector and returns that L2-norm
    void StoreAlgsysToVec(Vector<Double>& vec, Double * pt);

    //! calculates L2-norm of RHS regarding entries due to penalty formulation
    Double RhsL2Norm(Vector<Double>& stdVec);
  
    /// returns the time derivative of the solution belonging to all nodes of the actual element
    void GetDerivSolOfElement(Matrix<Double>& sol, StdVector<Integer>& connect_PDE);

    /// returns the vector of the solution belonging to all nodes of the actual element
    void GetSolVecOfElement(Vector<Double>& sol, StdVector<Integer>& connect_PDE);

    /// returns the vector of time derivative of the solution belonging to all nodes of the actual element
    void GetDerivSolVecOfElement(Vector<Double>& sol, StdVector<Integer>& connect_PDE);

    /// returns the vector of 2nd time derivative of the solution belonging to all nodes 
    /// of the actual element
    void GetDeriv2SolVecOfElement(Vector<Double>& sol, StdVector<Integer>& connect_PDE);
  

    //!
    void sortStresses(Vector<Double>& unsorted, Vector<Double>& sorted);
    void sortStresses(Vector<Complex>& unsorted, Vector<Complex>& sorted);

    // ======================================================
    // DATA SECTION
    // ======================================================
  
    // -----------------------------------------------------------------------
    // Geometry & node numbering
    // -----------------------------------------------------------------------
    //@{
    //! \name Attributes related to geometry and node numbering

    Grid * ptgrid_;        //!< pointer to grid object

    //! subdomain-levels belonging to PDE
    StdVector<std::string> subdoms_;
  
    //! surface-domain-levels belongig to PDE
    StdVector<std::string> surfdoms_;
  
    StdVector<std::string> pressSurf_;  //!< surface of pressure loads
    StdVector<Double>      pressVals_;  //!< values of the pressure loads
    StdVector<std::string> pressFnc_;   //!< function names of pressure loads

    //! pointer to boundary condition object
    BCs *ptBCs_;
    //@}

  
    // -----------------------------------------------------------------------
    // Geometry & node numbering
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes related to geometry and node numbering
    EQNType eqnType_;      //!< type of equation numbering used
    Integer dofspernode_;  //!< number of unknowns per node
    Integer dofsperedge_;  //!< number of unknowns per edge
    Integer numPDENodes_;  //!< number of nodes in subdomains
    Integer numElems_;     //!< number of elements in subdomains
  
    //! defines subtype of mechanic PDE: plainStrain, 3d, ...
    std::string subType_;
    //@}

    // -----------------------------------------------------------------------
    // Boundary conditions
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes connected to the handling of boundary conditions
  
    //! names of interfaces with homogeneous Dirichlet boundary conditions
    StdVector<std::string> bcs_hd_;
  
    //! names of interfaces with inhomogeneous Dirichlet boundary conditions
    StdVector<std::string> bcs_id_;
  
    //! names of surfaces with homogeneous von Neumann boundary conditions
    StdVector<std::string> bcs_nh_;
  
    //! names of surfaces with inhomogeneous von Neumann boundary conditions
    StdVector<std::string> bcs_ni_;
  
    //! names of surfaces with inhomogeneous Robin boundary conditions
    StdVector<std::string> bcs_rh_;
  
    //! names of surfaces with inhomogeneous Robin boundary conditions
    StdVector<std::string> bcs_ri_;
  
    //! values of solution/magnitude at inhomogeneous Dirichlet boundaries
    StdVector<Double> val_id_;
  
    //! values of phases at inhomogeneous Dirichlet boundaries
    StdVector<Double> bcs_id_phase_;
  
    //! names of the time functions for inhomogeneous Dirichlet BCs
    StdVector<std::string> fncnames_id_;
  
    //! degrees of freedom (e.g. ux) of homogenous Dirichlet BC
  
    //! This is a vector of strings, which describe the degrees of freedom
    //! (e.g. ux) that are fixed at a certain interface by the homogenous
    //! Dirichlet boundary conditions.
    StdVector<std::string> homDirichDof_;
  
    //! This is a vector of strings, which describe the degrees of freedom
    //! (e.g. ux) that are fixed at a certain interface by the inhomogenous
    //! Dirichlet boundary conditions.
    StdVector<std::string> inhomDirichDof_;

    //! nnumber of dirichlet boundary conditions
    Integer numDirichletBCs_;

    //@}

    // -----------------------------------------------------------------------
    // Non-linearity
    // -----------------------------------------------------------------------

    //@{
    //! \name Attributes connected to nonlinearity
    Boolean nonLin_;           //!< flag for nonlinear calculations
    Boolean geoUpdate_;        //!< flag for geometric update
    Double incStopCrit_;       //!< stopping criterion for incremental error
    Double residualStopCrit_;  //!< stopping criterion for residual error
    Integer nonLinMaxIter_;    //!< maximal number of NL-iterations
    std::string nonLinMethod_; //!< method for handling the non-linearity
    Boolean nonLinLogging_;    //!< log progress of non-linear iterations

    std::string lineSearch_;   //!< switch for lineSearch
    StdVector<NonLinPDE> nonLinPDEName_;//!< some PDEs carry a name (->acoustics!)
    //@}

    // -----------------------------------------------------------------------
    // Material data
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes handling info on material data
    // Char * charMaterialFileName_;    //!< name of material file
    // LoadMaterialData *loadMaterial_; //!< material reader
    MaterialData *materialData_;     //!< material data structure
    std::string pdematerialclass_;    //!< material class
  
    //@}


    // -----------------------------------------------------------------------
    // PDE coupling
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes connected to handling PDE coupling
    Boolean pdeIsCoupled_;        //!< PDE couples with others
    Matrix<Double> deltCoords_;    //!< offset to grid coordinates
    Vector<Double> matParam_;      //!< change to material parameter
    Boolean updateCouplingBCs_ ;  //!< flag if coupling BC were already set
    Integer couplingBCsCounter_;  //!< counter for number of coupling BCs
    PDECoupling *ptCoupling_;     //!< pointer to coupling object
  
    //! nodes at which coupling terms are calculated
    std::list<Integer> couplingNodes;    
  
    //! elements at which coupling terms are calculated
    StdVector<Elem*> couplingElements;
  
    //! iteration counter for coupled PDE solution process
    Integer iterCoupledCounter_;
    //@}


    // -----------------------------------------------------------------------
    // Time stepping
    // -----------------------------------------------------------------------
  
    //@{
    //! \name Attributes connected to time stepping
    TimeStepping * TS_alg_;       //!< handles the time stepping
    Boolean effectiveMass_;       //!< use effective mass formulation for transient analysis
    Boolean firstTimeStepStatic_; //!< needed for coupled, iterative methods
    Double lasttimecalc_;    //!< Last time on which we have calculated solution
  
    //! Number of last timestep on which we have calculated our solution
    Integer laststepcalc_;
  
    Double  actFrequency_; //!< current frequency for harmonic analysis
    Integer actFreqStep_;  //!< current frequency step for harmonic analysis
    //@}


    FileType * inFile_;       //!< pointer to input file
    WriteResults * outFile_;  //!< pointer to output file


    // -----------------------------------------------------------------------
    // Miscellaneous
    // -----------------------------------------------------------------------

    //@{
    //! \name Miscellaneous attributes
    AnalysisType analysistype_; //!< analysis type
    Boolean isAlwaysStatic_;    //!< flag for static PDEs (like electrostatic)
    ShortInt dim_;              //!< space dimension of pde
    Boolean isaxi_;             //!< TRUE: axisymmetric problem
    Boolean recalc_;            //!< flag indicating reassembling of system matrix
    Boolean isComplex_;         //!< true, if some part of PDE is complex (Material, solution)
    NodeEQN * eqnData_;         //!< equation handling

    BaseNodeStoreSol * sol_;    //!< solution

    //! specifies the type of damping model (see environment.hh)
    DampingType dampingType_;
    Boolean needsDampingMatrix_;
  
    Vector<Double> solIncr_;      //! needed in iterative coupled computation 
    Vector<Double> actSol_;       //! needed in iterative coupled computation 
    Boolean isIncrFormulation_;   //! checks, if we have for the coupling a incremental solution
    

    //! factors for constructing effective mass / stiffness matrix
    Double matrix_factor_[4];

    //! actual level (actual = current?)
    Integer actlevel_;

    //! Pointer to object of analysis (Static, Trans, Harm or Eig)
    Assemble * assemble_;
    BaseSystem * algsys_;     //!< pointer to algebraic system
    TimeFunc * ptTimeFunc_;   //!< pointer to time functions

  
    OLAS_Params * olasParams_; //! pointer to paramter object of OLAS
    OLAS_Report * olasReport_; //! pointer to report object of OLAS
    //@}

}; // class StdPDE

} // end of namespace

#endif
