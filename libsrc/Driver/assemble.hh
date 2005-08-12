#ifndef FILE_ANALYSIS
#define FILE_ANALYSIS

/**************************************************************************/
/* File:   analysis.hh                                                    */
/* Author: Fred Hofei                                                     */
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
#include "PDE/nodeEQN.hh"

#include "olas.hh"

namespace CoupledField
{

  // Forward declaration of StdPDE
  class StdPDE;
  class SinglePDE;

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

    //!
    void SetPDEIds(SinglePDE * aPDE1, SinglePDE * aPDE2) {
      myPDE1_ = aPDE1;
      myPDE2_ = aPDE2;
    }

    //! returns pointer to PDE1
    SinglePDE * GetPDE1 () 
    { return myPDE1_; };

    //! returns pointer to PDE2
    SinglePDE * GetPDE2 () 
    { return myPDE2_; };

    //! set function for SetCounterPart
    void SetCounterPart(bool setCounterPart) 
    { setCounterPart_ = setCounterPart; };

    //! check, if element matrix has to be assembled to upper and lower part of global matrix
    bool IsSetCounterPart() {return setCounterPart_; };

  protected:

    /// pointer to integrator
    BaseForm * integrator;
  
    /// is the integrator a nonlinear one?
    Boolean nonLin;

    //! reduced integration flag
    Boolean reducedIntegration_;

    SinglePDE * myPDE1_;                    //!< pointer to PDE 1
    SinglePDE * myPDE2_;                    //!< pointer to PDE 2

    //! TRUE: add to upper and lower part of global matrix
    bool setCounterPart_;
  };
  

  /// additional information for every integrator
  class IntegratorDescriptor : public BaseIntDescriptor {

  public:

    /// constructor
    IntegratorDescriptor();

    /// constructor
    IntegratorDescriptor(BaseForm * aIntegrator, FEMatrixType aDestMat,
                         const Boolean aNonLin=FALSE);
      
    /// destructor
    virtual ~IntegratorDescriptor();
      
    /// returns the destination matrix
    FEMatrixType DestMat() {return destinationMatrix;}; 

    /// sets the destination matrix
    void SetDestMat(FEMatrixType destMat)
    {destinationMatrix = destMat;};

    /// defines a secondary destination for the calculated element marices of an integrator      
    void SetSecondaryMat(FEMatrixType aSecMat, Double aSecMatFac,
                         AnalysisType analysisType);

    // returns matrix type of the secondary matrix (if there is any, otherwise NOTYPE=0)
    FEMatrixType GetSecondaryMat() const {return secondaryMatrix;} 

    /// returns matrix type of the secondary matrix (if there is any, otherwise NOTYPE=0)
    Double GetSecMatFac() const {return secMatFac;} 

    /// piezoMaterialType_ contains information wheather real
    /// or complex valued material parameters will be considered
    piezoMaterialType piezoMaterialType_;

    piezoMaterialType GetPiezoMaterialType(){return piezoMaterialType_;};

    void SetPiezoMaterialType(piezoMaterialType &pMatType){
      piezoMaterialType_ = pMatType;};

    void SetOrigMatrixType(FEMatrixType matType)
    {origMatrixType_ = matType;};

    FEMatrixType GetOrigMatrixType()
    {return origMatrixType_;};

    void SetOrigSecMatrixType(FEMatrixType matType)
    {origSecondMatrixType_ = matType;};

    FEMatrixType GetOrigSecMatrixType()
    {return origSecondMatrixType_;};
      
  private:

    /// holds the destination matrix
    FEMatrixType destinationMatrix;

    /// holds the secondary destination matrix
    FEMatrixType secondaryMatrix;

    //! hold the original matrix types (just used in harmonic analysis!!)
    FEMatrixType origMatrixType_;

    //! hold the original secondary matrix types (just used in harmonic analysis!!)
    FEMatrixType origSecondMatrixType_;

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
    
    /// adds integrators to the pde
    virtual void AddIntegrator(IntegratorDescriptor * intDescr,
                               const RegionIdType subdomain)=0;

    /// adds surface integrators to the pde
    virtual void AddSurfIntegrator(IntegratorDescriptor * intDescr,
                               const RegionIdType subdomain)=0;
    
    //! specify type of system matrix for AlgebraicSystem
    virtual void AssembleMatrices();
    

    /// setup source term
    void AssembleSrcRHS(const Double time=0.0);
    

    /// assemble integral sources
    void AssembleRHSIntegralSources(const Double time = 0.0);


    /// assembling nodal sources
    void AssembleRHSNodalSources(const Double time = 0.0);
    

    ///  assemble a nonlinear RHS part
    void AssembleNLRHS(const Double time = 0.0);

    ///  assemble a spring into the system matrix
    void AssembleSprings(const Double time = 0.0);

    //! computes the coordinates of an element including the delta
    /*!
      \param connect (input) global node numbers of element
      \param ptCoord (output) coordinates of the element nodes (nrNodes \f$\times\f$ spaceDim);
    */
    void GetElemCoords(const StdVector<UInt> connect, 
                       Matrix<Double> &coordMat); 



    //! set boundary condition
    /*!
      \param update indicator: do we update boundary condition in algebraic system or set new
      \param atimestep time step of calculation
    */
    virtual void SetBCs(const Boolean update, const Double atimestep){};



    
    //!
    //     void SetPtrDeltaCoordinates(Array<Double>* deltCoords)
    void SetPtrDeltaCoordinates(Matrix<Double> * deltCoords)
    {deltaCoords_ = deltCoords;};

    Matrix<Double> * GetPtrDeltaCoordinates(){return deltaCoords_;}

    void SetNonlinGeo()
    { nonLinGeo = TRUE;};


    /// sets the number of nodes in the pde    
    void SetNumberPDENodes(UInt aNumPDENodes)
    { numPDENodes_=aNumPDENodes; };


    //! compute rhs
    /*! \param atime time of calculation    */
    virtual void ComputeRHS(const Double atime) {;};
    

    /// parameters set by PDE
    void SetGeneralParams(const std::string & pdename, 
                          const UInt dofsPerNode,
                          const StdVector<RegionIdType> & subdoms,
                          const StdVector<RegionIdType> & surfdoms,
                          const std::string bcSequenceId);
    
    

    /// define RHS integrators
    void AddRhsIntegrator(BaseForm * integrator, const RegionIdType regionId, 
                          const Boolean nonLin=FALSE);

    /// define RHS integrators (static and transient case)
    void AddRhsSrcIntegrator(BaseForm * integrator, const RegionIdType regionId,                       
                             const std::string fncname="---not-defined--",
                             const Boolean nonLin=FALSE);

    /// define RHS surface integrators (static and transient case)
    void AddRhsSrcSurfIntegrator(BaseForm * integrator, const RegionIdType regionId,
                                 const std::string fncname="---not-defined--",
                                 const Boolean nonLin=FALSE);

    /// define RHS integrators (harmonic case)
    void AddRhsSrcIntegrator(BaseForm * integrator, const RegionIdType regionId,
                             const Double phaseval, const Boolean nonLin=FALSE);

    /// define RHS surface integrators (harmonic case)
    void AddRhsSrcSurfIntegrator(BaseForm * integrator, const RegionIdType regionId,
                                 const Double phaseval, const Boolean nonLin=FALSE);  

    /// set ptr to time function
    void SetPtr2TimeFnc(TimeFunc * aPtTimeFunc)
    {ptTimeFunc_ = aPtTimeFunc;};

    // set ptr to equation data
    void SetPtr2EQNData(NodeEQN * aPtNodeEQN1,
                        NodeEQN * aPtNodeEQN2 = NULL );

    //! set the identification tag of the pde
    //! \param id1 identification tag for PDE 1
    //! \param id2 identification tag for PDE 2 (only in coupled case)
    void SetPDEId( const PdeIdType id1,
                   const PdeIdType id2 = NO_PDE_ID );

  

    // ======================================================
    // STUFF BELONGING TO ALGSYS (matrices, graphs, ...)
    // ======================================================
  

    // TODO: Delete the following mathod
    //! deletes the algebraic system 
    //void DeleteAlgSys(){if (algsys_) delete algsys_;};

    
    // TODO: Delete the following mathod
    /// Initialize all necessary matrices 
    //void InitMatrices();


    /// Initialize all matrices with nonlinear behavior
    void InitNonLinMatrices();


    // TODO: Delete the following mathod
    /// establish matrices
    //void CreateMatrices();


    //Matrix<Double> GetElemMat(){return elemmat;};
    


    //! define discrete PDE
    virtual void MatrixSettings() = 0;

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

    //! Generates the matrix graph

    //! The method generates the matrix graph by passing the element
    //! connectivities to the algebraic system.
    //! \param insertCounterPart by default we always insert the complete
    //!                          connectivity of an for a PDE coupling, i.e.
    //!                          for forward and backward coupling; setting
    //!                          this flag to false will only insert the
    //!                          forward coupling
    void SetupMatrixGraph( bool insertCounterPart = true );


    //! set information for algebraic system about PDE. set matrix factors
    virtual void SetMatrixFactors()=0;

    void SetMaterialPointer(MaterialData * ptMat){ptMaterial_ = ptMat;};    

    void SetAlternatingMaterial(Boolean boolVar){alternateMaterialData_=boolVar;};

    
    /// set solution 
    void SetPtr2Sol(BaseNodeStoreSol * aSol){sol_ = aSol;};
    

    //! sets the actual frequency (just needed for harmonic analysis)
    virtual void SetFrequency(Double actFreq) {;};

    //! transform element matrix to account for spezial RHS during parameter Identification process
    virtual void TransformMatrix2HarmonicRHS_for_paramIdent(Vector<Double>& harmMat,
                                                            Matrix<Double> origMat) {;};

    //sets all finite elements to reduced integration
    void SetFE2ReducedInt();


    //sets all finite elements back to standard integration
    void SetFE2StandardInt();

    //! set the PDE pointer
    void SetPDEPointer(StdPDE * aptPDE)
    {ptPDE_ = aptPDE;};

    //! returns private variable matArray_
    Matrix<Double>* GetPDEMatArray (){
      return   matArray_;
    }

    // ====================================================
    // DATA SECTION 
    // ====================================================

        
    StdPDE * ptPDE_;                    //!< pointer to class StdPDE
    
    BaseSystem * algsys_;                //!< pointer to algebraic system  
    Grid * ptgrid_;                      //!< pointer to Grid
    NodeEQN * ptEQN1_;                    //!< pointer to equation data of pde1
    NodeEQN * ptEQN2_;                    //!< pointer to equation data of pde2
    Vector<Double> harmonicRHSVec;       //! special right Hand Side Vector needed for calc
    Matrix<Double> elemmat;
    OLAS_Params * olasParams_;               //!< pointer to parameter object of OLAS
    OLAS_Report * olasReport_;               //!< pointer ro report object of OLAS

    //! paramters for discrete PDE
    MatrixStructureType structuretype_;  //!< type of Matrix (SuperBlockMarix=SBM, Standard)
    MatrixEntryType entryType_;          //!< type of matrix entries (double, complex)
    MatrixStorageType storageType_;      //!< storage type of matrix (sparse, symmetric,..)

    //! set defining which type of matrices (stiffness, mass,...) is used
    std::set<FEMatrixType> matrixTypes_;

    UInt dofsPerNode_;                //!< number of unknowns per node
    UInt numPDENodes_;                //!< number of nodes in pde

    std::string pdename_;                //!< name of calling pde
    std::string bcSequenceTag_;          //!< name of tag for loads/boundary condition
  
    StdVector<RegionIdType> subdoms_;  //!< subdomain-levels belongig to PDE
    StdVector<RegionIdType> surfdoms_; //!< surface-domain-levels belongig to PDE

    StdVector<std::string> loadDom_;  //!< load subdomains
    StdVector<std::string> loadDof_;  //!< dofs of loads
    StdVector<Double>      loadVals_; //!< values of the load condition
    StdVector<std::string> fncname_loads_; //!< function names of the loads
    StdVector<std::string> fncname_rhs_; //!< function names for RHS integrators
    StdVector<std::string> fncname_rhsSurf_; //!< function names for RHS surface integrators

    StdVector<std::string> springDom_;  //!< spring subdomains
    StdVector<std::string> springDof_;  //!< dofs of springs
    StdVector<Double>      springMassVals_; //!< values of the spring mass
    StdVector<Double>      springDampVals_; //!< values of the spring damping
    StdVector<Double>      springStiffVals_; //!< values of the spring stiffness
    StdVector<std::string> fncname_springs_; //!< function names of the loads

    StdVector<Double> rhsSrcPhase_;      //!< contains the phase values in harmonic case;
    StdVector<Double> rhsSrcSurfPhase_;  //!< contains the phase values in harmonic case;
         
    TimeFunc * ptTimeFunc_;             //!< ptr to time function
    
    MaterialData * ptMaterial_;              //!< pointer to material
    
    
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
    UInt nrMatrices_;
    StdVector<Boolean> reassembleMat_;
    Boolean nonLinGeo;

    //! Default = FALSE, if content of ptMaterial changes, set TRUE
    Boolean alternateMaterialData_;

    Double actFreq_; //!< contains the frequency multiplied by 2*pi
    Double startFreq_; //!< contains the starting frequency multiplied by 2*pi within a harmonic analysis

    //!
    void SetReassemble()
    {firstTime_ = TRUE;};

    //!
    void SetAnalysisType(AnalysisType analysis)
    {analysisType_ = analysis;};

    //!
    void SetStartFrequency(Double freq)
    {startFreq_ = 2*PI*freq;};

    // ==============================================
    // AUXILIARY METHODS
    // ==============================================

    //! set material parameter array used in forms
    void SetMaterialArray(Matrix<Double>* matOfElements) 
    { matArray_ = matOfElements;};

  protected:
    //! calculates the index of the subdoman 
    //! with name "subDomName" in the subdomain-list
    UInt SubDomIndex(const RegionIdType subDomName);

    //! calculates the index of the surfdoman 
    //! with name "surfDomName" in the surface-domain-list
    UInt SurfDomIndex(const RegionIdType surfDomName);

    //! transform element matrix to account for harmonic analysis
    virtual 
    void TransformMatrix2Harmonic(Vector<Double>& harmMat,
                                  Matrix<Double> origMat,
                                  const FEMatrixType matrixType,
                                  const piezoMaterialType piezoMatType)
    {;};

    //! transform element vector to account for harmonic analysis
    virtual void TransformVector2Harmonic(Vector<Double>& harmMat,
                                          Vector<Double> origVec,
                                          const Double valPhase) {;};

    //! identifier for the first PDE

    //! needed to identify the PDE uniquely in the algebraic system
    PdeIdType pdeId1_;

    //! identifier for the second PDE

    //! needed to identify the PDE uniquely in the algebraic system
    PdeIdType pdeId2_;

    //!
    Matrix<Double>* matArray_;


  private:

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

    /// adds integrators to the pde
    virtual void AddIntegrator(IntegratorDescriptor * intDescr, 
                               const RegionIdType subdomain);

    /// adds surface integrators to the pde
    virtual void AddSurfIntegrator(IntegratorDescriptor * intDescr, 
				   const RegionIdType subdomain);

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
    virtual void AddIntegrator(IntegratorDescriptor * intDescr, 
                               const RegionIdType subdomain);

    /// adds surface integrators to the pde
    virtual void AddSurfIntegrator(IntegratorDescriptor * intDescr, 
                               const RegionIdType subdomain);
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

    //! transform element matrix to account for harmonic analysis
    virtual void TransformMatrix2Harmonic(Vector<Double>& harmMat,
                                          Matrix<Double> origMat,
                                          const FEMatrixType matrixType,
                                          const piezoMaterialType 
                                          piezoMatType);

    //! transform element vector to account for harmonic analysis
    virtual void TransformVector2Harmonic(Vector<Double>& harmMat, 
                                          Vector<Double> origVec,
                                          const Double valPhase);

    /// adds integrators to the pde
    virtual void AddIntegrator(IntegratorDescriptor * intDescr, 
                               const RegionIdType subdomain);

    /// adds surface integrators to the pde
    virtual void AddSurfIntegrator(IntegratorDescriptor * intDescr, 
				   const RegionIdType subdomain);

  };


} // end of namespace

#endif
