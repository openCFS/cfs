#ifndef FILE_STORESOLUTION_2004
#define FILE_STORESOLUTION_2004

#include "basestoresol.hh"
#include <Matrix/matrix.hh>
#include <Utils/vector.hh>

namespace CoupledField{

//! This class is the new interface for handling the 
//! (multidimensional) solutions of PDEs (instead of 
//! old Array-class).
  
template<class TYPE>
class StoreSol : public BaseStoreSol{
public:

  //! Constructor for empty class
  StoreSol();

  //! Constructor
  StoreSol(Integer numNodes, 
	   std::vector<SolutionType> solutionTypes, 
	   std::vector<Integer> solutionDofs);
  
  //! Copy Constructor
  StoreSol(const StoreSol & x);

  //! Destructor
  virtual ~StoreSol();

  //! initializes the object with given value
  void Init(const TYPE val);

  //!
  void SetNumSolutions(const Integer nSols);
  
  //!
  void SetNumNodes(const Integer nNodes);

  //!
  void SetSolutionType(const SolutionType solTypes, const Integer numSolution = 0);

  //!
  void SetNumDofs(const Integer dof, const SolutionType sol = NO_SOLUTION_TYPE);

  //!
  TYPE& operator()(Integer node, Integer dof);

  //!
  TYPE operator()(Integer node, Integer dof) const;

  //! 
  void GetSolVector(const SolutionType type, CFSVector & val) const;

  //!
  void SetSolVector(const SolutionType type, const CFSVector & val) ;

  //!
  void GetSolution(const SolutionType type, BaseStoreSol & val) const;

  //!
  void SetNodalResult(const Integer nodeNr, const CFSVector &val);
 
  //!
  void GetNodalResult(const Integer nodeNr, CFSVector & val) const;


  //! Get vector of a single solution of one special degree of freedom
  void GetSolVectorSingleDof(const SolutionType type, const Integer dof, CFSVector & val) const;
  
  //! Get vector of one special degree of freedom
  void GetSolVectorSingleDof(const Integer dof, CFSVector & val) const;

  //!
  void Get(const Integer nodeNr, const Integer dof, TYPE & ret) const;

  //!
  void Get(const SolutionType type, const Integer nodeNr, const Integer dof, TYPE & ret) const;

  //!
  void Set(const SolutionType type, const Integer nodeNr, const Integer dof, const TYPE val);

  //!
  void Add(const SolutionType type, const Integer nodeNr, const Integer dof, const TYPE val) const;

  //!
  void SetCompleteVector(const CFSVector & val);

  //!
  void GetCompleteVector(CFSVector & val) const;

  //!
  void GetVectorPointer(CFSVector* &ptrToVec);

  //!
  void GetVectorPointer(Vector<TYPE>* &ptrToVec);

  //!
  void CopyFromDataPointer(TYPE * ptr);
  
  //!
  void SetDataPointer(TYPE * ptr);
  
  //! 
  void GetDataPointer(TYPE * &ptr);
  
  //!
  Double * GetDoublePointer();

  //!
  Vector<TYPE>&  GetCompleteVector()
  {return data_;}

  /////////////////////////////////////////
  // Transformation routines for mapping //
  // pde to mesh solution and vice versa //
  /////////////////////////////////////////

  //!
  void GetElemSolutionAsMatrix(CFSMatrix & elemSol, Vector<Integer> & connect) const;

  //!
  void TransformNodeSolution(BaseStoreSol & transformedSolution,
			     const std::vector<Integer> & mapping,
			     Grid * ptGrid,
			     const Integer level) const;
  
  //!
  void TransformElemSolution(BaseStoreSol & transformedSolution,
			     const std::vector<Elem*> & elems,
			     Grid * ptGrid,
			     const Integer level) const;
  
  //! maps the local element solution to the global mesh solution
  /*!
    \param MeshSol (output) Solution vector referring to mesh node numbers    
    \param PDESol (input) Solution vector
    \param Elems (input) Vector of subdomains to which PDESol belongs to
     */
  void TransformElemSolution(BaseStoreSol & MeshSol, 
			     const std::vector<std::string> & SD,
			     Grid * ptGrid,
			     const Integer level) const;
  
  
  //! maps the local node solution to the coupling nodes
  void NodeSolutionToCoupling(BaseStoreSol & couplingSol,
			      const std::vector<Integer>& nodeNumbers,
			      const std::vector<Integer> & mapping) const;
  
  //! maps the local element solution to the coupling nodes
  void ElemSolutionToCoupling(BaseStoreSol & couplingSol,
			      const std::vector<Elem*>& elements,
			      const CFSVector & elemSol) const;
  

  ///////////////
  // Operators //
  ///////////////
  
  //! assignment operator
  StoreSol & operator= (const StoreSol & x);
  
  //! assignment operator for Base-clas
  BaseStoreSol & operator= (const BaseStoreSol & x);
  

protected:
  
  //! contains the solution itself
  Vector<TYPE> data_;
    
};

 
// Explicit template instantiation 
// needed for linking 
template class StoreSol<Integer>;
template class StoreSol<Double>;
template class StoreSol<Complex>;
 


} //end of namespace

#endif
