#ifndef FILE_ABSTRACTAS_2001
#define FILE_ABSTRACTAS_2001
 
namespace CoupledField
{

class AbstractAlgSys
{
public:
  //! Constructor with parameter - pointer to Grid and tolerancy for alg. system
 AbstractAlgSys(Grid<Point2D> * aptgrid, const Double aeps) { ptGrid=aptgrid; eps=aeps;}

 //! Restore solution after applying penalty method for zero boundary condition
 virtual void Restore()=0 ; 

  //! Assemble system matrix from element matrix and boundary condition
virtual void AssembleSysMatrix(const Double CoefLaplace, const Double CoefMass)=0;

   //! Assemble right hand side
  virtual void SetRHS(const Vector<Double> & CoefM=Vector<Double>(), const Vector<Double> & CoefS=Vector<Double>(), const Vector<Double> & f=Vector<Double>())=0;

  //! Set Dirichlet condition to right hand side, penalty method
  virtual void SetDirichletBoundaryCondRHS_PenaltyMethod(const Double val_tf)=0;

   //! Set Dirichlet condition to system matrix, penalty method
 virtual void SetDirichletBoundaryCondSysMat_PenaltyMethod()=0;

    //! Set zero Dirichlet boundary condition, cut row and colomn
  virtual void SetDirichletBoundaryCondZero_Cut()=0;

  //! Define array of nodes at which we have Dirichlet boundary condition
  virtual void SetNodesBoundaryCondition(FileType * aptFileType)=0;

    //! Conjugate Gradient for positive symmetric matrix
  virtual Boolean CG(const Integer maxIter, enum precond typePrecond)=0;

  //! Generalized minimum residual method
  virtual Boolean GMRes_m(const Integer maxIter, enum precond typePrecond, const Integer m)=0;

  //! BiCGSTAB method, can be used for all types of matrix
//  virtual Boolean BiCGSTAB(const Integer maxIter, enum precond typePrecond)=0;

  //! Return size of system matrix
   virtual Integer getSize()=0;

  // Print system matrix and right hand side
  virtual void printAb(std::ostream * out, const std::string & title) const =0;

  //Print solution
   virtual void printx(std::ostream * out, const Double time) const = 0;

  // Return solution
   virtual  Vector<Double> getSolution()=0;

protected:
  
  Grid<Point2D> * ptGrid;
  Double eps;
} ;

} // end of namespace

#endif
