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

  //! Class for anaylsis handling

  class Analysis
  {
    
  public:
    
    //!  Constructor
    Analysis(Grid * aptgrid);
    
    
    //!  Deconstructor
    virtual ~Analysis();
    
    
    /// adds integrators to the pde
    void AddIntegrator(BaseForm * integrator, const std::string & subdomain,
		       const enum MatrixType destinationMatrix, const Integer nonLin);



     //! specify type of system matrix for AlgebraicSystem
     /*! \param level (input) level of Grid     */
     virtual void SetupMatrices(const Integer level);



    /// setup source term
     void SetupRHS(const Integer level);



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
			  const std::vector<std::string> subdoms);
    
    


    // ======================================================
    // STUFF BELONGING TO ALGSYS (solvers, matrices, graphs, ...)
    // ======================================================
  

    //! define algebraic system 
    void SetAlgSys();
    

    //! deletes the algebraic system 
    void DeleteAlgSys();

    
    /// Initialize all necessary matrices 
    void InitMatrices();


    //! define discrete PDE
    virtual void MatrixSettings() = 0;

    
    //! Create the matrices and Solver as well as Preconditioner
    void CreateMatrices_Solver();

    //! Sets the type of the matrix graph
    void SetGraphType(enum GraphType aGraphType)
    {graphType_ = aGraphType;};
    
    //! constructes the matrix graph by providing to the algebraic system the element connectivities
    void SetupMatrixGraph(Integer numeq, Integer graphtype);


     //! set information for algebraic system about PDE. set matrix factors
     virtual void SetMatrixFactors()=0;


    //! define algebraic system solver Parameters
    virtual void SetSolverParameters();


    //! returns the local PDE number of an array of nodes
    /*!
      \param PDENodes (output) Vector of PDE node numbers
      \param MeshNodes (input) Vector of mesh node numbers
      \param Mesh2PDENode (input) Vector assigning PDE to mesh node numbers    
    */
    void Mesh2PDENode(Vector<Integer> & PDENodes, 
		      const Vector<Integer> & MeshNodes,
		      const std::vector<Integer> & Mesh2PDENode);
    

    /// allocates matrix space
    void NeedMatrix(const enum MatrixType matType)
    {algsys_->InitMatrix(matType);};



    void SetMesh2PDENode(std::vector<Integer> * aMesh2PDENode)
    {mesh2PDENode_ = aMesh2PDENode;};


    // set actual calculation level
    void SetLevel(Integer alevel)
    {actlevel_ = alevel;};
    

    // "time stepping" for solver
    virtual void SolveStep()=0;
    



    // ====================================================
    // DATA SECTION 
    // ====================================================


    /// additional information for every integrator
    struct integratorDescriptor
    {
      BaseForm * integrator;
      enum MatrixType destinationMatrix;
      Integer nonLin;
    };
    
    
    /// vector of all needed integrators (every subdomain needs one "list of integrators")
    //    Vector< Vector<integratorDescriptor *>* > integrators_;
    std::vector< std::vector<integratorDescriptor *>* > integrators_;

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
    Integer actlevel_;             //!< actual level of calculation
    




    // ==============================================
    // AUXILIARY METHODS
    // ==============================================

  private:
    //! calculates the index of the subdoman with name "subDomName" in the subdomain-list
    Integer SubDomIndex(const std::string & subDomName);
  };
    
      



  class StaticAnalysis : public Analysis
  {
    
  public:
    StaticAnalysis(Grid * agrid):Analysis(agrid){};
    
    virtual ~StaticAnalysis(){};

    //! define discrete PDE
    virtual void MatrixSettings();

     //! set information for algebraic system about PDE. set matrix factors
     virtual void SetMatrixFactors(){};

    /// "time stepping" for solver
    void SolveStep();

  };


} // end of namespace
#endif
