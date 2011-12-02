#ifndef FILE_CFS_FEFUNCTION_HH
#define FILE_CFS_FEFUNCTION_HH

#include "General/Environment.hh"

#include "Domain/ElemMapping/EntityLists.hh"
#include "Domain/BCs.hh"
#include "Domain/Results/ResultInfo.hh"
#include "Domain/Mesh/Grid.hh"

#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"

#include "Utils/mathParser/mathParser.hh"

namespace CoupledField {


  // forward class declarations
  class FeSpace;
  class SinglePDE;
  class MathParser;
  class BaseResult;
  class AlgebraicSys;

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


class BaseFeFunction {
public:
  

  //! Constructor
  BaseFeFunction();
  
  //! Destructor
  virtual ~BaseFeFunction();
  
  //! Finalize initialization
  virtual void Finalize() = 0;
  
  // ========================================================================
  //  Function Meta Information
  // ========================================================================
  //@{ \name Function Meta Information 
    
  //! Set result information
  void SetResultInfo( shared_ptr<ResultInfo> info );
  
  //! Get ResultInfo
  shared_ptr<ResultInfo> GetResultInfo();
  
  //! Set FeSpace
  void SetFeSpace( shared_ptr<FeSpace> space );
  
  //! Get FeSpace
  shared_ptr<FeSpace> GetFeSpace();
  
  //! Add EntityList
  void AddEntityList( shared_ptr<EntityList> list );
  
  //! Get EntityList
  StdVector< shared_ptr<EntityList> > GetEntityList();

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
  void AddInhomNeumannBC( shared_ptr<InhomNeumannBc> bc );

  //! Add constraint boundary condition
  void AddConstraint( shared_ptr<Constraint> bc );

  //! Get Homogenious Boundary Conditions
  const HdBcList GetHomDirichletBCs(){
    return hdBcs_;
  }

  //! Get Inhomogenious Dirichlet Boundary Conditions
  const IdBcList GetInHomDirichletBCs(){
    return idBcs_;
  }

  //! Get Inhomogenious Neumann Boundary Conditions
  const InBcList AddInhomNeumannBC(){
    return inBcs_;
  }
  
  //! Get Constraint Boundary Conditions
  const ConstraintList GetConstraints(){
    return constraints_;
  }

  //@}
  
  // ========================================================================
  //  Function Values
  // ========================================================================
  //@{ \name Function Values
  
  //! Get solution for specific entity
  virtual void GetEntitySolution( SingleVector& elemSol, 
                        const EntityIterator& it ) = 0;
                        
  //! Get solution as matrix for specific entity
  virtual void GetEntitySolutionAsMatrix( DenseMatrix& elemSol,
                                  const EntityIterator& it ) = 0;
  
  //! Copy values into result object
  virtual void ExtractResult( shared_ptr<BaseResult> res ) = 0;
  
  //! Return vector containing the function coefficients
  virtual SingleVector* GetSingleVector() = 0;
  
  //! Compute the BC values and hand them over to the Algebraic System
   virtual void ApplyBC() = 0;
  //@}
  

protected:

  //! Identifier for the function
  FeFctIdType fctId_;
  
  //! Pointer to finite element function space
  shared_ptr<FeSpace> feSpace_;
  
  //! Entitylists (elements, nodes, etc.) the function is defined on
  StdVector<shared_ptr<EntityList> > entities_;
  
  //! Homogeneous Dirichlet BCs
  HdBcList hdBcs_;
  
  //! Inhomogeneous Dirichlet BCs
  IdBcList idBcs_;
  
  //! Inhomogeneous Neumann BCs
  InBcList inBcs_;

  //! Constraints
  ConstraintList constraints_;
  
  //! The ResultInfo for the Function
  shared_ptr<ResultInfo> result_;

  //! Pointer to the creating PDE
  //shared_ptr<StdPDE> pde_;
  SinglePDE* pde_;

  //! Pointer to the grid 
  Grid* grid_;
  
  //! Handle for MathParser object
  MathParser::HandleType mHandle_;

  //! pointer to algebraic system
  AlgebraicSys* algsys_;      

};


///////////////////////////////////////////////////////////////////
// Templatized Version
///////////////////////////////////////////////////////////////////

//! Templatized FeFunction (real/complex)
//TODO: Concept for Applying higher order Boundary conditions 
template<typename T>
class FeFunction : public BaseFeFunction {
public:
  FeFunction();

  virtual ~FeFunction();

  //! \see BaseFeFunction::Finalize()
  void Finalize();
  
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
  SingleVector* GetSingleVector() {return &coeffs_;}
    
  //! Return complete coefficient vector
  Vector<T>& GetVector() {return coeffs_;}
  
  //! Return complete coefficient vector (const)
  const Vector<T>& GetVector() const {return coeffs_;}

  virtual void ApplyBC();
  
  
protected:

  //! Coefficient vector
  Vector<T> coeffs_;
};


}  // namespace CoupledField
#endif
