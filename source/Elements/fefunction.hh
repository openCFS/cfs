// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;
#ifndef FILE_CFS_FEFUNCTION_HH
#define FILE_CFS_FEFUNCTION_HH

#include "General/environment.hh"
#include "MatVec/vector.hh"
#include "Domain/entityList.hh"
#include "Domain/bcs.hh"
#include "Domain/resultInfo.hh"
#include "MatVec/matrix.hh"
//#include "PDE/StdPDE.hh"
#include "Domain/grid.hh"
#include "Utils/mathParser/mathParser.hh"

namespace CoupledField {


  // forward class declarations
  class FeSpace;
//  class StdPDE;
  class SinglePDE;
  class MathParser;
//  class BaseSystem;

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
  ~BaseFeFunction();
  
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

  //! Set the PDE pointer of the function
  void SetPDE(shared_ptr<SinglePDE> pde);

  //! Get the PDE Pointer
  shared_ptr<SinglePDE> GetPDE();

  //! Set the PDE pointer of the function
  void SetGrid(shared_ptr<Grid>  grid);

  //! Get the PDE Pointer
  shared_ptr<Grid>  GetGrid();

  //! Set the algebraic System
  void SetSystem( shared_ptr<BaseSystem> sys );

  //! Get the algebraic System
  shared_ptr<BaseSystem> GetSystem();

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
                        const EntityIterator& it );
                        
  //! Get solution as matrix for specific entity
  virtual void GetEntitySolutionAsMatrix( DenseMatrix& elemSol,
                                  const EntityIterator& it );
  
  //! Compute the BC values and hand them over to the Algebraic System
  virtual void ApplyBC(){};
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
  shared_ptr<SinglePDE> pde_;

  //! Pointer to the grid 
  shared_ptr<Grid> grid_;
  
  //! Handle for MathParser object
  MathParser::HandleType mHandle_;

  //! pointer to algebraic system
  shared_ptr<BaseSystem> algsys_;      

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

  ~FeFunction();

  /////////////////////////////////////////////////////////////////
  // Solution Access functions 
  /////////////////////////////////////////////////////////////////
  //
  //! Get solution for specific entity
  void GetEntitySolution( SingleVector& elemSol, 
                        const EntityIterator& it );
                        
  //! Get solution as matrix for specific entity
  void GetEntitySolutionAsMatrix( DenseMatrix& elemSol,
                                  const EntityIterator& it );

  virtual void ApplyBC(){};
protected:

  //! Coefficient vector
  Vector<T> coeffs_;
};

/////////////////////////////////////////////////////////////////
// Specialized version for  different treatment of Boundary Conditions
/////////////////////////////////////////////////////////////////
template<>
class FeFunction<Double> : public BaseFeFunction {
public:
  FeFunction(){};

  ~FeFunction(){};
    virtual void ApplyBC();
};

template<>
class FeFunction<Complex> : public BaseFeFunction {
public:
  FeFunction(){};

  ~FeFunction(){};
    virtual void ApplyBC();
};


}  // namespace CoupledField
#endif
