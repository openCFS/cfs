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

#include <Forms/baseForm.hh>
 
namespace CoupledField
{


  /// additional information for every integrator
class BaseIntDescriptor
{
public:
  /// constructor
  BaseIntDescriptor();
      
  /// constructor
  BaseIntDescriptor(BaseForm * aIntegrator, const Boolean aNonLin=FALSE);
  
  /// destructor
  virtual ~BaseIntDescriptor();
  
  /// is the integrator nonlinear?
  Boolean IsNonLin() {return nonLin;};
  
  /// returns the integrator
  BaseForm * GetIntegrator(){return integrator;};
  
  
protected:
  /// pointer to integrator
  BaseForm * integrator;
  
  /// is the integrator a nonlinear one?
  Boolean nonLin;
};
  
  


  /// additional information for every integrator
class IntegratorDescriptor : public BaseIntDescriptor
    {
    public:
      /// constructor
      IntegratorDescriptor();
      
      /// constructor
      IntegratorDescriptor(BaseForm * aIntegrator,  enum MatrixType aDestMat, const Boolean aNonLin=FALSE);

      /// destructor
      virtual ~IntegratorDescriptor();
      

      /// returns the destination matrix
      MatrixType DestMat() {return destinationMatrix;};

      /// sets the destination matrix
      void SetDestMat(enum MatrixType destMat)
      {destinationMatrix = destMat;};
      

      /// defines a secondary destination for the calculated element marices of an integrator      
      void SetSecondaryMat(enum MatrixType aSecMat, Double aSecMatFac)
      {
	secondaryMatrix = aSecMat;
	secMatFac = aSecMatFac;
      };

      /// returns matrix type of the secondary matrix (if there is any, otherwise NOTYPE=0)
      MatrixType GetSecondaryMat() const {return secondaryMatrix;} 

      /// returns matrix type of the secondary matrix (if there is any, otherwise NOTYPE=0)
      Double GetSecMatFac() const {return secMatFac;} 
      
      
      
      
    private:
      /// holds the destination matrix
      enum MatrixType destinationMatrix;

      /// holds the secondary desination matrix
      enum MatrixType secondaryMatrix;

      /// holds the matrix factor for secondaryMatrix
      Double secMatFac;
    };
  




  //! Class for anaylsis handling
  class Assemble
  {
    
  public:
    
    //!  Constructor
    Assemble(BaseSystem * algsys, Grid * aptgrid);
    
    
    //!  Deconstructor
    virtual ~Assemble();
    
    
    /// adds integrators to the pde
    virtual void AddIntegrator(BaseForm * integrator, const std::string & subdomain,
			       const enum MatrixType destinationMatrix, const Integer nonLin)=0;

    /// adds integrators to the pde
    virtual void AddIntegrator(IntegratorDescriptor * intDescr, const std::string & subdomain)=0;


    /// adds surface integrators to the pde
    virtual void AddSurfIntegrator(BaseForm * integrator, const std::string & subdomain,
			       const enum MatrixType destinationMatrix, const Integer nonLin)=0;

    //! specify type of system matrix for AlgebraicSystem
    /*! \param level (input) level of Grid     */
    virtual void AssembleMatrices(const Integer level);
    

    
    
    /// setup source term
    void AssembleRHS(const Integer level, const Double time=0);
    

    /// assemble integral sources
    void AssembleRHSIntegralSources(const Integer level, const Double time=0);


    /// assembling nodal sources
    void AssembleRHSNodalSources(const Integer level, const Double time=0);
    
  

    //! computes the coordinates of an element including the delta
    /*!
      \param connect (input) global node numbers of element
      \param ptCoord (output) coordinates of the element nodes (nrNodes \f$\times\f$ spaceDim);
      \param level (input) index for multilevel hierarchy
    */
    void GetElemCoords(const Vector<Integer> connect, 
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
    void SetPtrDeltaCoordinates(Array<Double>* deltCoords)
    {deltaCoords_ = deltCoords;};

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
			  const std::vector<std::string> subdoms,
			  const std::vector<std::string> surfdoms);
    
    

    void SetNumDirichlet(Integer numDirichlet)
    {numDirichletBCs_ = numDirichlet;};

    /// define RHS integrators
    void AddRhsIntegrator(BaseForm * integrator, const std::string & subDomName, 
			  const Integer nonLin=FALSE);

     /// define RHS integrators
    void AddRhsSrcIntegrator(BaseForm * integrator, const std::string & subDomName, 			  
			     const std::string fncname="---not-defined--",
			     const Integer nonLin=FALSE);
   
    /// set ptr to time function
    void SetPtr2TimeFnc(TimeFunc * aPtTimeFunc)
    {ptTimeFunc_ = aPtTimeFunc;};
  

    // ======================================================
    // STUFF BELONGING TO ALGSYS (matrices, graphs, ...)
    // ======================================================
  

    //! define algebraic system 
    void SetAlgSys();
    

    //! deletes the algebraic system 
    void DeleteAlgSys(){if (algsys_) delete algsys_;};

    
    /// Initialize all necessary matrices 
    void InitMatrices();


    /// establish matrices
    void CreateMatrices();
    


    //! define discrete PDE
    virtual void MatrixSettings() = 0;


    //! define matrix type
    virtual void SetMatrixType(Integer matType)
    {matrixType_=matType;};


    //! return matrix type
    Integer GetMatrixType(){return matrixType_;};
    

    //! Sets the type of the matrix graph
    void SetGraphType(enum GraphType aGraphType)
    {graphType_ = aGraphType;};
    

    //! constructes the matrix graph by providing to the algebraic system the element connectivities
    //    void SetupMatrixGraph(Integer numeq, Integer graphtype);
    void SetupMatrixGraph(Integer numeq);


     //! set information for algebraic system about PDE. set matrix factors
     virtual void SetMatrixFactors()=0;



    //! returns the local PDE number of an array of nodes
    /*!
      \param PDENodes (output) Vector of PDE node numbers
      \param MeshNodes (input) Vector of mesh node numbers
      \param Mesh2PDENode (input) Vector assigning PDE to mesh node numbers    
    */
    void Mesh2PDENode(Vector<Integer> & PDENodes, 
		      const Vector<Integer> & MeshNodes,
		      const std::vector<Integer> & Mesh2PDENode);
    


    void SetMesh2PDENode(std::vector<Integer> * aMesh2PDENode)
    {mesh2PDENode_ = aMesh2PDENode;};


    /// set actual calculation level
    void SetLevel(Integer alevel){actlevel_ = alevel;};
    

    /// sets the pointer to the BCs
    void SetPtrBCs(BCs* aptBCs) { ptBCs_ = aptBCs;}      


    /// require damping matrix
    void NeedDampingMatrix(){dampingMatrix_ = TRUE;};
    

    /// set solution 
    void SetPtr2Sol(Array<Double> * aSol){sol_ = aSol;};
    

    /// extracts solution belonging to nodes in connect_PDE
    void GetSolOfElement(Matrix<Double>& elDisp, Vector<Integer>& connect_PDE);

    /// return index to dof
    Integer GetBCDof(const std::string dofString);
    

    // ====================================================
    // DATA SECTION 
    // ====================================================
    
    BaseSystem * algsys_;         //!< pointer to algebraic system  
    Grid * ptgrid_;               //!< pointer to Grid

    //! paramters for discrete PDE
    Boolean matrixType_;          //!< type of matrix (real, complex, etc.)
    Boolean graphType_;           //!< type of graph (nodal, edge,..)
    Boolean systemMatrix_;        //!< need system matrix (TRUE/FALSE)
    Boolean stiffnessMatrix_;     //!< need stiffness matrix (TRUE/FALSE)
    Boolean dampingMatrix_;       //!< need damping matrix (TRUE/FALSE)
    Boolean massMatrix_;          //!< need mass matrix (TRUE/FALSE)
    Boolean convectionMatrix_;    //!< need convective matrix (TRUE/FALSE)

    Integer dofsPerNode_;         //!< number of unknowns per node
    Integer numDirichletBCs_;     //!< number of dirichlet boundary conditions
    Integer numPDENodes_;         //!< number of nodes in pde

    std::string pdename_;         //!< name of calling pde
    std::vector<Integer> * mesh2PDENode_; //!< array containing PDE (=local) node numbers

    std::vector<std::string> subdoms_;  //!< subdomain-levels belongig to PDE
    std::vector<std::string> surfdoms_; //!< surface-domain-levels belongig to PDE

    std::vector<std::string> loadDom_;  //!< load subdomains
    std::vector<std::string> loadDof_;  //!< dofs of loads
    std::vector<Double>      loadVals_; //!< values of the load condition
    std::vector<std::string> fncname_loads_; //!< function names of the loads
    std::vector<std::string> fncname_rhs_; //!< function names for RHS integrators

    TimeFunc * ptTimeFunc_;             //!< ptr to time function
    
    BCs *ptBCs_;                       //!< pointer to Boundary Condition  Object

    Integer actlevel_;                 //!< actual level of calculation
    
    
    /// vector of all needed integrators (every subdomain needs one "list of integrators")
    std::vector< std::vector<IntegratorDescriptor *>* > integrators_;

    /// vector of all needed surface integrators (every surface needs one "list of surfaceintegrators")
    std::vector< std::vector<IntegratorDescriptor *>* > surfintegrators_;

    /// vector of all needed integrators (every subdomain needs one "list of integrators")
    std::vector< std::vector<BaseIntDescriptor *>* > rhsIntegrators_;

    /// vector of all needed RHS src-intergators (not every subdomain needs a "list of rhs_source_integrators")
    std::vector< std::vector<BaseIntDescriptor *>* > rhsSrcIntegrators_;

    /// ptr to solution
    Array<Double> * sol_;
    Array<Double> * deltaCoords_;

    //! nonlinear parameters;
    Boolean firstTime_;
    Boolean oneIntIsNonlin_;
    Integer nrMatrices_;  
    Vector<Boolean> reassembleMat_;
    Boolean nonLinGeo;

    void SetReassemble()
    {firstTime_ = TRUE;};

    // ==============================================
    // AUXILIARY METHODS
    // ==============================================



  protected:
    //! calculates the index of the subdoman with name "subDomName" in the subdomain-list
    Integer SubDomIndex(const std::string & subDomName);

    //! calculates the index of the surfdoman with name "surfDomName" in the surface-domain-list
    Integer SurfDomIndex(const std::string & surfDomName);

  private:

    //! returns the index of the named dof
    Integer GetNrBCDof (const std::string & dofStartString);


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

    virtual void AddIntegrator(BaseForm * integrator, const std::string & subdomain,
			const enum MatrixType destinationMatrix, const Integer nonLin);

    /// adds integrators to the pde
    virtual void AddIntegrator(IntegratorDescriptor * intDescr, const std::string & subdomain);

    /// adds surface integrators to the pde
    virtual void AddSurfIntegrator(BaseForm * integrator, const std::string & subdomain,
			       const enum MatrixType destinationMatrix, const Integer nonLin);

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

    /// adds integrators to the pde
     virtual void AddIntegrator(BaseForm * integrator, const std::string & subdomain,
				const enum MatrixType destinationMatrix, const Integer nonLin);

    /// adds integrators to the pde
    virtual void AddIntegrator(IntegratorDescriptor * intDescr, const std::string & subdomain);

    /// adds surface integrators to the pde
    virtual void AddSurfIntegrator(BaseForm * integrator, const std::string & subdomain,
				   const enum MatrixType destinationMatrix, const Integer nonLin);

  };

} // end of namespace
#endif
