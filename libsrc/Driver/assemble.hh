#ifndef FILE_ANALYSIS
#define FILE_ANALYSIS

/**************************************************************************/
/* File:   analysis.hh                                                    */
/* Author: Fred Hofer                                                     */
/* Date:   24. Nov. 2003                                                  */
/*                                                                        */
/* Handles the assembling of the various integrators to the appropriate   */
/* matrices and initiates the basic analysis "stepping"                   */
/**************************************************************************/


#include "Utils/StdVector.hh"
#include "Forms/baseForm.hh" 
#include "DataInOut/timefunc.hh"
#include "Utils/nodestoresol.hh"
#include "Utils/elemstoresol.hh"
#include "Domain/bcs.hh"
#include "PDE/nodeEQN.hh"

namespace CoupledField
{


  //! additional information for every integrator
  class BaseIntDescriptor {

  public:

    /// constructor
    BaseIntDescriptor();
      
    /// constructor
    BaseIntDescriptor(BaseForm * aIntegrator, const Boolean aNonLin=FALSE);
  
    /// destructor
    virtual ~BaseIntDescriptor();
  
    /// is the integrator nonlinear?
    Boolean IsNonLin() {return nonLin;};

    /// is the integrator set to reduced integration?
    Boolean IsReducedInt() {return reducedIntegration_;};

    /// set the integrator to reduced integration
    void SetReducedInt() {reducedIntegration_ = TRUE;};
  
    /// returns the integrator
    BaseForm * GetIntegrator(){return integrator;};

  protected:

    /// pointer to integrator
    BaseForm * integrator;
  
    /// is the integrator a nonlinear one?
    Boolean nonLin;

    //! reduced integration flag
    Boolean reducedIntegration_;
  };
  

  /// additional information for every integrator
  class IntegratorDescriptor : public BaseIntDescriptor {

  public:

    /// constructor
    IntegratorDescriptor();

#ifdef USE_OLAS      
    /// constructor
    IntegratorDescriptor(BaseForm * aIntegrator, FEMatrixType aDestMat,
			 const Boolean aNonLin=FALSE);
#else
    /// constructor
    IntegratorDescriptor(BaseForm * aIntegrator,  enum MatrixType aDestMat,
			 const Boolean aNonLin=FALSE);
#endif
      
    /// destructor
    virtual ~IntegratorDescriptor();
      
#ifdef USE_OLAS
    /// returns the destination matrix
    FEMatrixType DestMat() {return destinationMatrix;}; 
#else
    /// returns the destination matrix
    MatrixType DestMat() {return destinationMatrix;};
#endif


#ifdef USE_OLAS
    /// sets the destination matrix
    void SetDestMat(FEMatrixType destMat)
    {destinationMatrix = destMat;};
#else
    /// sets the destination matrix
    void SetDestMat(enum MatrixType destMat)
    {destinationMatrix = destMat;};
#endif
 
#ifdef USE_OLAS
    /// defines a secondary destination for the calculated element marices of an integrator      
    void SetSecondaryMat(FEMatrixType aSecMat, Double aSecMatFac,
			 AnalysisType analysisType);

#else
    /// defines a secondary destination for the calculated element marices of an integrator      
    void SetSecondaryMat(enum MatrixType aSecMat, Double aSecMatFac,
			 AnalysisType analysisType);
#endif

#ifdef USE_OLAS
    // returns matrix type of the secondary matrix (if there is any, otherwise NOTYPE=0)
    FEMatrixType GetSecondaryMat() const {return secondaryMatrix;} 
#else
    /// returns matrix type of the secondary matrix (if there is any, otherwise NOTYPE=0)
    MatrixType GetSecondaryMat() const {return secondaryMatrix;} 
#endif

    /// returns matrix type of the secondary matrix (if there is any, otherwise NOTYPE=0)
    Double GetSecMatFac() const {return secMatFac;} 
      
#ifdef USE_OLAS  
    void SetOrigMatrixType(FEMatrixType matType)
    {origMatrixType_ = matType;};

    FEMatrixType GetOrigMatrixType()
    {return origMatrixType_;};

    void SetOrigSecMatrixType(FEMatrixType matType)
    {origSecondMatrixType_ = matType;};

    FEMatrixType GetOrigSecMatrixType()
    {return origSecondMatrixType_;};

#else
    void SetOrigMatrixType(MatrixType matType)
    {origMatrixType_ = matType;};

    MatrixType GetOrigMatrixType()
    {return origMatrixType_;};

    void SetOrigSecMatrixType(MatrixType matType)
    {origSecondMatrixType_ = matType;};

    MatrixType GetOrigSecMatrixType()
    {return origSecondMatrixType_;};
#endif
      
  private:

#ifdef USE_OLAS
    /// holds the destination matrix
    FEMatrixType destinationMatrix;

    /// holds the secondary destination matrix
    FEMatrixType secondaryMatrix;

    //! hold the original matrix types (just used in harmonic analysis!!)
    FEMatrixType origMatrixType_;

    //! hold the original secondary matrix types (just used in harmonic analysis!!)
    FEMatrixType origSecondMatrixType_;
#else
    /// holds the destination matrix
    enum MatrixType destinationMatrix;

    /// holds the secondary destination matrix
    enum MatrixType secondaryMatrix;

    //! hold the original matrix types (just used in harmonic analysis!!)
    MatrixType origMatrixType_;

    //! hold the original secondary matrix types (just used in harmonic analysis!!)
    MatrixType origSecondMatrixType_;
#endif

    /// holds the matrix factor for secondaryMatrix
    Double secMatFac;
  };
  

  //! Class for anaylsis handling
  class Assemble {
    
  public:
    
    //!  Constructor
    Assemble(BaseSystem * algsys, Grid * aptgrid);
    
    //!  Deconstructor
    virtual ~Assemble();
    
#ifdef USE_OLAS
    /// adds integrators to the pde
    virtual void AddIntegrator(BaseForm * integrator,
			       const std::string & subdomain,
			       const FEMatrixType destinationMatrix,
			       const Integer nonLin)=0;
#else
    /// adds integrators to the pde
    virtual void AddIntegrator(BaseForm * integrator,
			       const std::string & subdomain,
			       const enum MatrixType destinationMatrix,
			       const Integer nonLin)=0;
#endif

    /// adds integrators to the pde
    virtual void AddIntegrator(IntegratorDescriptor * intDescr,
			       const std::string & subdomain)=0;

#ifdef USE_OLAS
    /// adds surface integrators to the pde
    virtual void AddSurfIntegrator(BaseForm * integrator,
				   const std::string & subdomain,
				   const FEMatrixType destinationMatrix,
				   const Integer nonLin)=0;
#else
    /// adds surface integrators to the pde
    virtual void AddSurfIntegrator(BaseForm * integrator, const std::string & subdomain,
				   const enum MatrixType destinationMatrix, const Integer nonLin)=0;
#endif

    //! specify type of system matrix for AlgebraicSystem
    /*! \param level (input) level of Grid     */
    virtual void AssembleMatrices(const Integer level);
    

    /// setup source term
    void AssembleSrcRHS(const Integer level, const Double time=0);
    

    /// assemble integral sources
    void AssembleRHSIntegralSources(const Integer level, const Double time=0);


    /// assembling nodal sources
    void AssembleRHSNodalSources(const Integer level, const Double time=0);
    

    ///  assemble a nonlinear RHS part
    void AssembleNLRHS(const Integer level, const Double time=0);


    //! computes the coordinates of an element including the delta
    /*!
      \param connect (input) global node numbers of element
      \param ptCoord (output) coordinates of the element nodes (nrNodes \f$\times\f$ spaceDim);
      \param level (input) index for multilevel hierarchy
    */
    void GetElemCoords(const StdVector<Integer> connect, 
		       Matrix<Double> &coordMat, const Integer level); 



    //! set boundary condition
    /*!
      \param level level of grid
      \param update indicator: do we update boundary condition in algebraic system or set new
      \param atimestep time step of calculation
    */
    virtual void SetBCs(const Integer level, const Integer update, 
			const Double atimestep){};



    
    //!
    //     void SetPtrDeltaCoordinates(Array<Double>* deltCoords)
    void SetPtrDeltaCoordinates(Matrix<Double> * deltCoords)
    {deltaCoords_ = deltCoords;};

    Matrix<Double> * GetPtrDeltaCoordinates(){return deltaCoords_;}

    void SetNonlinGeo()
    { nonLinGeo = TRUE;};


    /// sets the number of nodes in the pde    
    void SetNumberPDENodes(Integer aNumPDENodes)
    { numPDENodes_=aNumPDENodes; };


    //! compute rhs
    /*! \param atime time of calculation    */
    virtual void ComputeRHS(const Double atime) {;};
    

    /// parameters set by PDE
    void SetGeneralParams(const std::string & pdename, 
			  const Integer dofsPerNode,
			  const Integer numPDENodes, 
			  const StdVector<std::string> subdoms,
			  const StdVector<std::string> surfdoms,
			  const std::string bcSequenceId);
    
    

    void SetNumDirichlet(Integer numDirichlet)
    {numDirichletBCs_ = numDirichlet;};

    /// define RHS integrators
    void AddRhsIntegrator(BaseForm * integrator, const std::string & subDomName, 
			  const Integer nonLin=FALSE);

    /// define RHS integrators (static and transient case)
    void AddRhsSrcIntegrator(BaseForm * integrator, const std::string & subDomName, 			  
			     const std::string fncname="---not-defined--",
			     const Integer nonLin=FALSE);

    /// define RHS surface integrators (static and transient case)
    void AddRhsSrcSurfIntegrator(BaseForm * integrator, const std::string & subDomName,
			     const std::string fncname="---not-defined--",
			     const Integer nonLin=FALSE);

    /// define RHS integrators (harmonic case)
    void AddRhsSrcIntegrator(BaseForm * integrator, const std::string & subDomName,
			     const Double phaseval, const Integer nonLin=FALSE);

     /// define RHS surface integrators (harmonic case)
    void AddRhsSrcSurfIntegrator(BaseForm * integrator, const std::string & subDomName,
			     const Double phaseval, const Integer nonLin=FALSE);  

    /// set ptr to time function
    void SetPtr2TimeFnc(TimeFunc * aPtTimeFunc)
    {ptTimeFunc_ = aPtTimeFunc;};

    // set ptr to equation data
    void SetPtr2EQNData(NodeEQN * aPtNodeEQN);
  

    // ======================================================
    // STUFF BELONGING TO ALGSYS (matrices, graphs, ...)
    // ======================================================
  

    //! define algebraic system 
    void SetAlgSys();
    

    //! deletes the algebraic system 
    void DeleteAlgSys(){if (algsys_) delete algsys_;};

    
    /// Initialize all necessary matrices 
    void InitMatrices();


    /// Initialize all matrices with nonlinear behavior
    void InitNonLinMatrices();


    /// establish matrices
    void CreateMatrices();

   Matrix<Double> GetElemMat(){return elemmat;};
    


    //! define discrete PDE
    virtual void MatrixSettings() = 0;

#ifdef USE_OLAS
    //! define entry type of matrices (DOUBLE, COMPLEX)
    virtual void SetMatrixEntryType(MatrixEntryType etype)
    {entryType_ = etype;};
    
    //! return entry type of the matrix
    MatrixEntryType GetMatrixEntryType()
    {return entryType_;};
    
    //! define storage type of matrices (SPARSE_SYM, SPARSE_NONSYM, ...)
    virtual void SetMatrixStorageType(MatrixStorageType stype)
    {storageType_ = stype;};

    //! return storage type of the matrix
    MatrixStorageType GetMatrixStorageType()
    {return storageType_;};
    
#else
    //! define matrix type
    virtual void SetMatrixType(Integer matType)
    {matrixType_=matType;};


    //! return matrix type
    Integer GetMatrixType(){return matrixType_;};
#endif    



#ifdef USE_OLAS
    //! Sets the type of the matrix graph
    void SetGraphType(GraphType aGraphType)
    {graphType_ = aGraphType;};
#else
    //! Sets the type of the matrix graph
    void SetGraphType(enum GraphType aGraphType)
    {graphType_ = aGraphType;};
#endif
    

    //! constructes the matrix graph by providing to the algebraic system the element connectivities
    //    void SetupMatrixGraph(Integer numeq, Integer graphtype);
    void SetupMatrixGraph(Integer numeq);


    //! set information for algebraic system about PDE. set matrix factors
    virtual void SetMatrixFactors()=0;


    /// set actual calculation level
    void SetLevel(Integer alevel){actlevel_ = alevel;};
    

    /// sets the pointer to the BCs
    void SetPtrBCs(BCs* aptBCs) { ptBCs_ = aptBCs;}      


    /// require damping matrix
    void NeedDampingMatrix(){dampingMatrix_ = TRUE;};
    

    /// set solution 
    void SetPtr2Sol(BaseNodeStoreSol * aSol){sol_ = aSol;};
    

    /// return index to dof
    Integer GetBCDof(const std::string dofString);
    

    //! sets the actual frequency (just needed for harmonic analysis)
    virtual void SetFrequency(Double actFreq) {;};

	//! transform element matrix to account for spezial RHS during parameter Identification process
    virtual void TransformMatrix2HarmonicRHS_for_paramIdent(Vector<Double>& harmMat,
					  Matrix<Double> origMat) {;};

    //sets all finite elements to reduced integration
    void SetFE2ReducedInt();


    //sets all finite elements back to standard integration
    void SetFE2StandardInt();

    // ====================================================
    // DATA SECTION 
    // ====================================================
    
    BaseSystem * algsys_;                //!< pointer to algebraic system  
    Grid * ptgrid_;                      //!< pointer to Grid
    NodeEQN * ptEQN_;                    //!< pointer to equation data
    Vector<Double> harmonicRHSVec;	 //! special right Hand Side Vector needed for calc
    Matrix<Double> elemmat;
#ifdef USE_OLAS
    OLAS_Params * olasParams_;               //!< pointer to parameter object of OLAS
    OLAS_Report * olasReport_;               //!< pointer ro report object of OLAS
#endif

    //! paramters for discrete PDE
#ifdef USE_OLAS
    MatrixStructureType structuretype_;  //!< type of Matrix (SuperBlockMarix=SBM, Standard)
    MatrixEntryType entryType_;          //!< type of matrix entries (double, complex)
    MatrixStorageType storageType_;      //!< storage type of matrix (sparse, symmetric,..)
    GraphType graphType_;                //!< type of graph (nodal, edge,..)
#else
    Boolean matrixType_;                 //!< type of matrix (real, complex, etc.)
    Boolean graphType_;                  //!< type of graph (nodal, edge,..)
#endif

    Boolean systemMatrix_;               //!< need system matrix (TRUE/FALSE)
    Boolean stiffnessMatrix_;            //!< need stiffness matrix (TRUE/FALSE)
    Boolean dampingMatrix_;              //!< need damping matrix (TRUE/FALSE)
    Boolean massMatrix_;                 //!< need mass matrix (TRUE/FALSE)
    Boolean convectionMatrix_;           //!< need convective matrix (TRUE/FALSE)

    Integer dofsPerNode_;                //!< number of unknowns per node
    Integer numDirichletBCs_;            //!< number of dirichlet boundary conditions
    Integer numPDENodes_;                //!< number of nodes in pde

    std::string pdename_;                //!< name of calling pde
    std::string bcSequenceTag_;          //!< name of tag for loads/boundary condition
  
    StdVector<std::string> subdoms_;  //!< subdomain-levels belongig to PDE
    StdVector<std::string> surfdoms_; //!< surface-domain-levels belongig to PDE

    StdVector<std::string> loadDom_;  //!< load subdomains
    StdVector<std::string> loadDof_;  //!< dofs of loads
    StdVector<Double>      loadVals_; //!< values of the load condition
    StdVector<std::string> fncname_loads_; //!< function names of the loads
    StdVector<std::string> fncname_rhs_; //!< function names for RHS integrators
    StdVector<std::string> fncname_rhsSurf_; //!< function names for RHS surface integrators

    StdVector<Double> rhsSrcPhase_;      //!< contains the phase values in harmonic case;
    StdVector<Double> rhsSrcSurfPhase_;  //!< contains the phase values in harmonic case;
         
    TimeFunc * ptTimeFunc_;             //!< ptr to time function
    
    BCs *ptBCs_;                       //!< pointer to Boundary Condition  Object

    Integer actlevel_;                 //!< actual level of calculation
    
    
    /// vector of all needed integrators (every subdomain needs one "list of integrators")
    StdVector< StdVector<IntegratorDescriptor *>* > integrators_;

    /// vector of all needed surface integrators (every surface needs one "list of surfaceintegrators")
    StdVector< StdVector<IntegratorDescriptor *>* > surfintegrators_;

    /// vector of all needed integrators (every subdomain needs one "list of integrators")
    StdVector< StdVector<BaseIntDescriptor *>* > rhsIntegrators_;

    /// vector of all needed RHS src-intergators (not every subdomain needs a "list of rhs_source_integrators")
    StdVector< StdVector<BaseIntDescriptor *>* > rhsSrcIntegrators_;

    /// vector of all needed RHS src-surfaceintergators (not every subdomain needs a "list of rhs_source_integrators")
    StdVector< StdVector<BaseIntDescriptor *>* > rhsSrcSurfIntegrators_;

    /// ptr to solution

    BaseNodeStoreSol * sol_;
    Matrix<Double> * deltaCoords_;

    //! nonlinear parameters;
    Boolean firstTime_;
    Boolean oneIntIsNonlin_;
    Integer nrMatrices_;
    StdVector<Boolean> reassembleMat_;
    Boolean nonLinGeo;

    Double actFreq_; //!< contains the frequency multiplied by 2*pi

    void SetReassemble()
    {firstTime_ = TRUE;};

    void SetAnalysisType(AnalysisType analysis)
    {analysisType_ = analysis;};

    // ==============================================
    // AUXILIARY METHODS
    // ==============================================



  protected:
    //! calculates the index of the subdoman with name "subDomName" in the subdomain-list
    Integer SubDomIndex(const std::string & subDomName);

    //! calculates the index of the surfdoman with name "surfDomName" in the surface-domain-list
    Integer SurfDomIndex(const std::string & surfDomName);

#ifdef USE_OLAS
    //! transform element matrix to account for harmonic analysis
    virtual void TransformMatrix2Harmonic(Vector<Double>& harmMat,
					  Matrix<Double> origMat,
					  const FEMatrixType matrixType) {;};
#else
    //! transform element matrix to account for harmonic analysis
    virtual void TransformMatrix2Harmonic(Vector<Double>& harmMat,
					  Matrix<Double> origMat,
					  const MatrixType matrixType) {;};
#endif

    //! transform element vector to account for harmonic analysis
    virtual void TransformVector2Harmonic(Vector<Double>& harmMat,
					  Vector<Double> origVec,
					  const Double valPhase) {;};




  private:

    //! returns the index of the named dof
    Integer GetNrBCDof (const std::string & dofStartString);

    //! set analysis type
    AnalysisType analysisType_;
  };
    
      



  class StaticAssemble : public Assemble
  {
  public:
    StaticAssemble(BaseSystem * algsys, Grid * agrid);
    
    virtual ~StaticAssemble(){};

    //! define discrete PDE
    virtual void MatrixSettings(){};

    //! set information for algebraic system about PDE. set matrix factors
    virtual void SetMatrixFactors(){};

#ifdef USE_OLAS
    virtual void AddIntegrator(BaseForm * integrator, const std::string & subdomain,
			       const FEMatrixType destinationMatrix, const Integer nonLin);
#else
    virtual void AddIntegrator(BaseForm * integrator, const std::string & subdomain,
			       const enum MatrixType destinationMatrix, const Integer nonLin);
#endif

    /// adds integrators to the pde
    virtual void AddIntegrator(IntegratorDescriptor * intDescr, const std::string & subdomain);

#ifdef USE_OLAS
    /// adds surface integrators to the pde
    virtual void AddSurfIntegrator(BaseForm * integrator, const std::string & subdomain,
				   const FEMatrixType destinationMatrix, const Integer nonLin);
#else
    /// adds surface integrators to the pde
    virtual void AddSurfIntegrator(BaseForm * integrator, const std::string & subdomain,
				   const enum MatrixType destinationMatrix, const Integer nonLin);
#endif


  };




  class TransientAssemble : public Assemble
  {
    
  public:
    TransientAssemble(BaseSystem * algsys, Grid * agrid);
    
    virtual ~TransientAssemble(){};

    //! define discrete PDE
    virtual void MatrixSettings(){};

    //! set information for algebraic system about PDE. set matrix factors
    virtual void SetMatrixFactors(){};  
    
#ifdef USE_OLAS
    /// adds integrators to the pde
    virtual void AddIntegrator(BaseForm * integrator, const std::string & subdomain,
			       const FEMatrixType destinationMatrix, const Integer nonLin);
#else
    /// adds integrators to the pde
    virtual void AddIntegrator(BaseForm * integrator, const std::string & subdomain,
			       const enum MatrixType destinationMatrix, const Integer nonLin);
#endif

    /// adds integrators to the pde
    virtual void AddIntegrator(IntegratorDescriptor * intDescr, const std::string & subdomain);

#ifdef USE_OLAS
    /// adds surface integrators to the pde
    virtual void AddSurfIntegrator(BaseForm * integrator, const std::string & subdomain,
				   const FEMatrixType destinationMatrix, const Integer nonLin);
#else
    /// adds surface integrators to the pde
    virtual void AddSurfIntegrator(BaseForm * integrator, const std::string & subdomain,
				   const enum MatrixType destinationMatrix, const Integer nonLin);
#endif

  };


  class HarmonicAssemble : public Assemble
  {
  public:
    HarmonicAssemble(BaseSystem * algsys, Grid * agrid);
    
    virtual ~HarmonicAssemble(){};

    //! define discrete PDE
    virtual void MatrixSettings(){};

    //! set information for algebraic system about PDE. set matrix factors
    virtual void SetMatrixFactors(){};

    //!
    virtual void SetFrequency(Double actFreq);

#ifdef USE_OLAS
    //! transform element matrix to account for harmonic analysis
    virtual void TransformMatrix2Harmonic(Vector<Double>& harmMat,
					  Matrix<Double> origMat,
					  const FEMatrixType matrixType);
#else
    //! transform element matrix to account for harmonic analysis
    virtual void TransformMatrix2Harmonic(Vector<Double>& harmMat,
					  Matrix<Double> origMat,
					  const MatrixType matrixType);
#endif


    //! transform element matrix to account for harmonic analysis
    virtual void TransformMatrix2HarmonicRHS_for_paramIdent(Vector<Double>& harmMat,
					  Matrix<Double> origMat);


    //! transform element vector to account for harmonic analysis
    virtual void TransformVector2Harmonic(Vector<Double>& harmMat, 
					  Vector<Double> origVec,
					  const Double valPhase);

#ifdef USE_OLAS
    virtual void AddIntegrator(BaseForm * integrator, const std::string & subdomain,
			       const FEMatrixType destinationMatrix, const Integer nonLin);
#else
    virtual void AddIntegrator(BaseForm * integrator, const std::string & subdomain,
			       const enum MatrixType destinationMatrix, const Integer nonLin);
#endif

    /// adds integrators to the pde
    virtual void AddIntegrator(IntegratorDescriptor * intDescr, const std::string & subdomain);

#ifdef USE_OLAS
    /// adds surface integrators to the pde
    virtual void AddSurfIntegrator(BaseForm * integrator, const std::string & subdomain,
				   const FEMatrixType destinationMatrix, const Integer nonLin);
#else
    /// adds surface integrators to the pde
    virtual void AddSurfIntegrator(BaseForm * integrator, const std::string & subdomain,
				   const enum MatrixType destinationMatrix, const Integer nonLin);
#endif


  };


} // end of namespace
#endif
