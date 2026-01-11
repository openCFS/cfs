#ifndef FILE_CFS_FEFUNCTION_HH
#define FILE_CFS_FEFUNCTION_HH

#include "General/Environment.hh"

#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/BCs.hh"
#include "Domain/Results/ResultInfo.hh"
#include "Domain/CoefFunction/CoefFunction.hh"
#include "Domain/Mesh/Grid.hh"

#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "Utils/ThreadLocalStorage.hh"

namespace CoupledField {


  // forward class declarations
  class FeSpace;
  class SinglePDE;
  class MathParser;
  class BaseResult;
  class AlgebraicSys;
  class BaseBOperator;
  class LinearForm;
  class BiLinearForm;
  class SBM_Vector;
  class BaseTimeScheme;
  class SimState;
  class MathParser;

//!  Base class for a function approximated by Finite Elements 
/*!
  The FeFunction class represents an function, approximated by in a finite
  element space (replaces the old Node-/ElemStoreSol class)
  Therefore it holds the following information:
    - Pointer to related Finite Element Space (which knows about the equation 
      numbers)
    - Description of the function (i.e. ResultInfo)
    - List of entities it is defined on (region, element list)
    - Boundary conditions (Dirichlet, Constraints, etc.)
    - Coefficient vector
  
  The class has methods to
  - Aquire the solution of an element
  - Add Result
*/


class BaseFeFunction : public CoefFunction, 
                       public boost::enable_shared_from_this<BaseFeFunction> {
public:
  
  //! Make simulation state class friend 
  friend class SimState;

  //! Constructor
  BaseFeFunction(MathParser* mp);
  
  virtual string GetName() const { return "BaseFeFunction"; }

  //! Destructor
  virtual ~BaseFeFunction();
  
  //! Finalize initialization
  virtual void Finalize() = 0;
  
  //! Cleanup date. To b called before destruction. 
  virtual void CleanUp() = 0;
  
  //! Query for complex-valued results
  virtual bool IsComplex() const = 0;
  
  /** query the constraints if there are periodic boundary conditions */
  bool HasPeriodicBC() const;

  // ========================================================================
  //  Function Meta Information
  // ========================================================================
  //@{ \name Function Meta Information 
    
  //! Set result information
  virtual void SetResultInfo( shared_ptr<ResultInfo> info ) = 0;
  
  //! Get ResultInfo
  shared_ptr<ResultInfo> GetResultInfo();
  
  //! Set FeSpace
  void SetFeSpace( shared_ptr<FeSpace> space );
  
  //! Get FeSpace
  shared_ptr<FeSpace> GetFeSpace();
  
  //! Add EntityList
  virtual void AddEntityList( shared_ptr<EntityList> list );
  
  //! Get EntityList
  virtual StdVector< shared_ptr<EntityList> > GetEntityList();
  
  //! Get regions the FeFunction is defined on
  const std::set<RegionIdType>& GetRegions() const;

  //! Set the function Id, as assigned by OLAS
  void SetFctId(const FeFctIdType id );
  
  //! Return the function I
  FeFctIdType GetFctId() const {return fctId_;}
  
  //! Set the PDE pointer of the function
  void SetPDE(SinglePDE* pde);

  //! Get the PDE Pointer
  SinglePDE* GetPDE();

  //! Set the PDE pointer of the function
  void SetGrid(Grid*  grid);

  //! Get the PDE Pointer
  Grid*  GetGrid();

  //! Set the algebraic System
  void SetSystem( AlgebraicSys* sys );

  //! Get the algebraic System
  AlgebraicSys * GetSystem();

  shared_ptr<BaseTimeScheme> GetTimeScheme();

  void SetTimeScheme(shared_ptr<BaseTimeScheme> scheme);
  
  //@}
  
  
  // ========================================================================
  //  Boundary Conditions
  // ========================================================================
  //@{ \name Boundary Conditions
  
  //! Add hom. Dirichlet boundary condition
  void AddHomDirichletBc( shared_ptr<HomDirichletBc> bc );

  //! Add inhom. Dirichlet boundary condition
  void AddInhomDirichletBc( shared_ptr<InhomDirichletBc> bc );
  
  //! Add constraint boundary condition
  void AddConstraint( shared_ptr<Constraint> bc );

  //@{
  //! Add a Coefficient function to fill the FeFunction
  void AddLoadCoefFunction( PtrCoefFct coef,
                            const StdVector<shared_ptr<EntityList> >& list);
  
  void AddLoadCoefFunction( PtrCoefFct coef,
                            shared_ptr<EntityList >& list);
  //@}

  //@{
  //! Add an external data function to fill the FE function
  void AddExternalDataSource( PtrCoefFct coef,
                              const StdVector<shared_ptr<EntityList> >& list);
  
  void AddExternalDataSource( PtrCoefFct coef,
                              shared_ptr<EntityList> list);
  //@}

  //! Remove external data sources
  void RemoveExternalDataSource();

  /** Get Homogeneous Boundary Conditions
  * Note that we return a reference to allow the modification of the list for optimization purpose */
  HdBcList& GetHomDirichletBCs(){
    return hdBcs_;
  }

  /** Get Inhomogenious Dirichlet Boundary Conditions.
   * Note that we return a reference to allow the modification of the list for optimization purpose */
  IdBcList& GetInHomDirichletBCs(){
    return idBcs_;
  }

  
  //! Get Constraint Boundary Conditions
  const ConstraintList GetConstraints(){
    return constraints_;
  }

  /** searches for constraint */
  bool HasConstraint(std::string& name, unsigned int dof) const;

  //! Get Load CoefFunctions
  const LoadCoefList GetLoadCoefFunctions(){
    return loadCoefs_;
  }

  //! Get external Data CoefFunctions
  const ExternalDataCoefList GetExternalDataSource(){
    return externalDataCoefs_;
  }
  //@}
  
  // ========================================================================
  //  Coefficient Function Interface
  // ========================================================================
  //@{ \name Access Methods of CoefFunction interface

  //! \copydoc CoefFunction::GetVecSize
  virtual UInt GetVecSize() const;

  //! \copydoc CoefFunction::GetTensorSize
  virtual void GetTensorSize( UInt& numRows, UInt& numCols ) const;
  
  //! \copydoc CoefFunction::ToString 
  virtual std::string ToString() const;
  
  //@}
  
  // ========================================================================
  //  Function Values
  // ========================================================================
  //@{ \name Function Values
  
  //! Get solution for specific entity
  virtual void GetEntitySolution( SingleVector& elemSol, 
                        const EntityIterator& it ) = 0;

  /** shortcut for GetEntitySolution() with EntityIterator */
  void GetEntitySolution(SingleVector& elemSol, const Elem* elem);


  //! Get solution as matrix for specific entity
  virtual void GetEntitySolutionAsMatrix( DenseMatrix& elemSol,
                                  const EntityIterator& it ) = 0;
  
  //! Copy values into result object
  virtual void ExtractResult( shared_ptr<BaseResult> res )=0;
  
  //! Return vector containing the function coefficients
  virtual SingleVector* GetSingleVector() = 0;
  
  //! Compute the BC values and hand them over to the Algebraic System
  virtual void ApplyBC() = 0;

  //! Incorporate load conditions, the characteristic here is that the values will be
  //! added, not set as in ApplyBCs
  virtual void ApplyLoads() = 0;
   
  //! Set the feFunction to values obtained by external data sources
  virtual void ApplyExternalData() = 0;
  //@}
  
  //! generates an interpolation operator by determining the space used
  virtual BaseBOperator* GenerateInterpolationOperator(UInt dim, UInt dofDim)=0;

  //! Generate interpolation bilinear form
  virtual BiLinearForm* GenerateInterpolBilinForm( UInt spaceDim, UInt dofDim ,
                                           bool updatedGeo )=0;
  
  //! Generate interpolation linear form
  virtual LinearForm* GenerateInterpolLinForm( UInt spaceDim, UInt dofDim, PtrCoefFct ,
                                           bool updatedGeo )=0;
  
  //! Copy the values from another FeFunction
  
  //! This method tries to copy the coefficient values from another FeFunction.
  //! If both associated spaces have the same approximation type, this is merely
  //! a simple equation mapping and copy operation. If the spaces have different
  //! approximation types, the general mapping mechanism will be utilized.
  virtual void InitFromFeFunction( shared_ptr<BaseFeFunction> feFct ) = 0;

  struct EqNodeGeom{
    UInt indexNum; //index number in OLAS
    UInt nodeNum; //one nodeNum can contain several eqNum's, e.g. 2D, 3D vector values
    Vector<Double> coord; //coordinate of the node with nodeNum
    //node-coordinates of both nodes of one edge, must be stored in positive direction
    StdVector< Vector<Double> > eCoords;
    StdVector<Integer> eNodes; // edge nodes

  };

  //! Create a connection between OLAS-equations and geometry

  //! Checks if nodal- (Lagrangian-) FE or edge- (Nédéléc-) is used
  //! and fills the appropriate info into the algebraic system
  virtual void ApplyGeomInfo() = 0;

protected:

  //! Identifier for the function
  FeFctIdType fctId_;
  
  //! Pointer to finite element function space
  shared_ptr<FeSpace> feSpace_;
  
  //! Set with all regions the function is defined on
  std::set<RegionIdType> regions_;
  
  //! Support of the CoefFunction. Only needed for grid/solution results
  StdVector<shared_ptr<EntityList> > entities_;
  
  /** to speed up AddEntityList() for periodic B.C. we use this structure.
   * This is meant for pairs of two node numbers. The first is the key for the map, the second
   * for the list. The list will be short and we insert only if not already contained.
   * A list has less overhead than a set and does not copy as an array. A C++11 forward_list would be more efficient*/
  std::map<unsigned int, std::list<unsigned int> > two_node_entries_cache_;

  //! Homogeneous Dirichlet BCs
  HdBcList hdBcs_;
  
  //! Inhomogeneous Dirichlet BCs
  IdBcList idBcs_;
  
  //! Constraints
  ConstraintList constraints_;
  
  //! LoadCoefFunctions
  LoadCoefList loadCoefs_;

  //! coefFunction providing external Data
  ExternalDataCoefList externalDataCoefs_;

  //! The ResultInfo for the Function
  shared_ptr<ResultInfo> result_;

  //! Pointer to the creating PDE
  //shared_ptr<StdPDE> pde_;
  SinglePDE* pde_;

  //! Pointer to the grid 
  Grid* grid_;
  
  //! Handle for MathParser object
  unsigned int mHandle_;

  //! pointer to algebraic system
  AlgebraicSys* algsys_;
  
  //! Pointer to MathParser
  MathParser * mp_;

  //! pointer to timestepping scheme
  shared_ptr<BaseTimeScheme> timeScheme_;
  
};


///////////////////////////////////////////////////////////////////
// Templatized Version
///////////////////////////////////////////////////////////////////

//! Templatized FeFunction (real/complex)
template<typename T>
class FeFunction : public BaseFeFunction {
public:
  
  //! Make simulation state class friend 
  friend class SimState;
  
  FeFunction(MathParser* mp);

  virtual ~FeFunction();

  //! Explicitly set the time derivative for the FeFunction (only complex)
  void SetTimeDerivOrder( UInt i, shared_ptr<FeFunction<T> > feFct  );
  
  //! \see BaseFeFunction::Finalize()
  void Finalize();
  
  //! \see BaseFeFunction::CleanUp 
  void CleanUp();
  
  
  virtual bool IsComplex() const {
    return std::is_same<T,Complex>::value;
    
  }
  //! \copydoc BaseFeFunction::SetResultInfo 
  void SetResultInfo( shared_ptr<ResultInfo> info );
  
  /////////////////////////////////////////////////////////////////
  // Solution Access functions 
  /////////////////////////////////////////////////////////////////
  //
  
  //! Copy values into result object
  void ExtractResult( shared_ptr<BaseResult> res ); 
  
  //! Get solution for specific entity
  void GetEntitySolution( SingleVector& elemSol, 
                          const EntityIterator& it );
  
   //! Get solution as matrix for specific entity
  void GetEntitySolutionAsMatrix( DenseMatrix& elemSol,
                                   const EntityIterator& it );
  
  
  //! Get solution for specific element 
  void GetElemSolution( Vector<T>& elemSol,
                        const Elem* elem );
  
  
  //! Return vector containing the function coefficients
  SingleVector* GetSingleVector() {return coeffs_;}
    
  //! Return complete coefficient vector
  Vector<T>& GetVector() {return *coeffs_;}
  
  //! Return complete coefficient vector (const)
  const Vector<T>& GetVector() const {return *coeffs_;}

  //! Incorporate boundary conditions
  virtual void ApplyBC();
  
  //! Incorporate load conditions, the characteristic here is that the values will be
  //! added, not set as in ApplyBCs
  virtual void ApplyLoads();

  //! incorporate loads defined in values
  virtual void ApplyLoads(PtrCoefFct& values);

  //! Set the feFunction to values obtained by external data sources
  virtual void ApplyExternalData();


  //! Create a connection between OLAS-equations and geometry

  //! Checks if nodal- (Lagrangian-) FE or edge- (Nédéléc-) is used
  //! and fills the appropriate info into the algebraic system
  virtual void ApplyGeomInfo();

  // ========================================================================
   //  Coefficient Function Interface
   // ========================================================================
   //@{ \name Access Methods of CoefFunction interface

  //! Return complex-valued vector at integration point
  virtual void GetVector(Vector<T>& vec, 
                         const LocPointMapped& lpm );

  //! Return real-valued element averaged value
  virtual void GetAvgElemValue(T & vec, 
                         const Elem* elem); 

  //! Return complex-valued scalar at integration point
  virtual void GetScalar(T& scalar, 
                         const LocPointMapped& lpm );
   
   //@}
  
  //! generates an interpolation operator by determining the space used
  BaseBOperator* GenerateInterpolationOperator(UInt dim, UInt dofDim);

  //! Generate interpolation bilinear form
  BiLinearForm* GenerateInterpolBilinForm( UInt spaceDim, UInt dofDim,
                                           bool updatedGeo );
  
  //! Generate interpolation linear form
  LinearForm* GenerateInterpolLinForm( UInt spaceDim, UInt dofDim, PtrCoefFct,
                                       bool updatedGeo );
  
  //! \copydoc BaseFeFunction::InitFromFeFunction
  void InitFromFeFunction( shared_ptr<BaseFeFunction> feFct );
  
protected:

  //! Update factor for time derivative (only complex case)
  void UpdateTimeDeriv();
  
  //! Coefficient vector
  Vector<T> * coeffs_;
  
  //! Interpolation operator
  BaseBOperator* idOp_;

  //! Time derivative order (only complex valued case)
  UInt timeDerivOrder_;
  
  //! Factor for time derivative
  T factor_;

  //! Thread-local work buffer for GetElemSolution to avoid per-call allocations
  mutable CfsTLS<StdVector<Integer>> work_eqns_;

  //! Thread-local cache for element solution to avoid redundant recomputation
  //! When the same element is queried multiple times (e.g., from different
  //! coefficient functions in UpdateXpr), the cached result is returned.
  mutable CfsTLS<UInt> lastElemNum_;
  mutable CfsTLS<Vector<T>> lastElemSol_;
};


}  // namespace CoupledField
#endif
