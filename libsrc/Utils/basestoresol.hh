#ifndef FILE_BASESTORESOL_2004
#define FILE_BASESTORESOL_2004

#include <Utils/tools.hh>
#include <vector>
#include <map>

// Not really needed
//#include <typeinfo>
//#include <General/defs.hh>



namespace CoupledField{

class Elem;
class Grid;
class CFSVector;
class CFSMatrix;
template<class TYPE> class Matrix;
template<class TYPE> class Vector;
template<class TYPE> class Array;


//////////////////
// Enumerations //
//////////////////
typedef enum{NOSOLUTIONTYPE, COUPLING_TYPE,\
	     MECH_DISPLACEMENT, MECH_ACCELERATION, MECH_VELOCITY,\
	     ELEC_POTENTIAL, ELEC_FIELD,\
	     SMOOTH_DISPLACEMENT,\
	     ACOU_POTENTIAL, ACOU_VELOCITY, ACOU_PRESSURE,\
	     MAG_POTENTIAL, MAG_FIELD, MAG_EDDY_CURRENT} SolutionType;
	       
//! This class is the new interface for handling the 
//! (multidimensional) solutions of PDEs (instead of 
//! old Array-class).
//! In principle this class is only a wrapper around 
//! a CFSVector class. Additionally it contains
//! informations about the information stored in it and
//! it can selectively be accessed, set and read.
  
class BaseStoreSol
{
public:

   
  //! Initializes the object AFTER the other layout
  //! functions have been called
  //! \note Only after calling Init, the object can
  //! store information
  virtual void Init(const Double val)
  {Error("BaseStoreSol::Init() not implemented here",__FILE__,__LINE__);}
  
  //! Set the number of different solutiontypes
  //! \note All entries of this object are deleted
  virtual void SetNumSolutions(const Integer nSols) = 0;
  
  //! Set the number of solution nodes
  //! \note All entries of this object are deleted
  virtual void SetNumNodes(const Integer nNodes) = 0;

  //! Set the solution types of the different solutions
  virtual void SetSolutionType(const SolutionType solType, Integer numSolution = 0) = 0;
  
  //!
  //! \note All entries of this object are deleted
  virtual void SetDof(const Integer dof, const SolutionType sol = NOSOLUTIONTYPE) = 0;

  //! Return the accumulated number of Dofs for
  //! all solution types
  inline Integer GetTotalNumDofs()
  {
    return totalDofs_;
  }

  //!
  inline Integer GetDof(const SolutionType sol = NOSOLUTIONTYPE) const
  {
    if (numSolutions_ == 1 && sol == NOSOLUTIONTYPE)
      return (*solDofs_.begin()).second;
    else
      return (*solDofs_.find(sol)).second;
  }
  
  //! 
  inline Integer GetNumNodes() const 
  {return numNodes_;}

  //!
  inline Integer GetSize() const 
  {return length_;}
  
  //!
  inline Integer GetNumSolutions() const 
  {return numSolutions_;}

  //! Return solution vector for all nodes of one type as a vector
  virtual void GetSolVector(const SolutionType type, CFSVector & val) const = 0;
   
  //! Set complete solution vector for all nodes of one solutiontype
  virtual void SetSolVector(const SolutionType type, const CFSVector & val) = 0;
 
  //! Set point velus
  virtual void SetNodalResult(const Integer nodeNr, const CFSVector &val) = 0;

  //! 
  virtual void GetNodalResult(const Integer nodeNr, CFSVector & val) const = 0;

  //! Get vector of a single solution of one special degree of freedom
  virtual void GetSolVectorSingleDof(const SolutionType type, const Integer dof, CFSVector & val) const = 0;
  
  //! Get vector of one special degree of freedom
  virtual void GetSolVectorSingleDof(const Integer dof, CFSVector & val) const = 0;

  //! Get single type of result in a new BaseStoreSol-class
  virtual void GetSolution(const SolutionType type, BaseStoreSol & val) const = 0;

  //! Get a single entry
  //! \NOTE: This method must only be used with solutions with only one entrytype
  virtual void Get(const Integer nodeNr, const Integer dof, Double & ret) const
  {Error("BaseStoreSol::Get() not implemented here",__FILE__,__LINE__);} 

  //! Get a single entry of a given solutiontype and a given dof
  virtual void Get(const SolutionType type, const Integer nodeNr, const Integer dof, Double & ret) const
  {Error("BaseStoreSol::Get() not implemented here",__FILE__,__LINE__);}

  //! Set a single entry of a given solutiontype and a given dof
  virtual void Set(const SolutionType type, const Integer nodeNr, const Integer dof, const Double val)
  {Error("BaseStoreSol::Set() not implemented here",__FILE__,__LINE__);}
 
  //! Ass a value tp  a single entry of a given solutiontype and a given dof
  virtual void Add(const SolutionType type, const Integer nodeNr, const Integer dof, const Double val) const
  {Error("BaseStoreSol::Add() not implemented here",__FILE__,__LINE__);} 
  
  //!
  virtual void SetCompleteVector(const CFSVector & val) = 0;

  //!
  virtual void GetCompleteVector(CFSVector & val) const = 0;
  
  //! Get the pointer to the basevector inside this class
  virtual void GetVectorPointer(CFSVector* &ptrToVec) = 0;

  //! Set the pointer to the data inside this class
  //! to an extern array of numbers.
  //! \NOTE: This method is very intrusive and should only be used
  //! when one can ensure, that the internal layout of the solution
  //! matches to the one of the given array. This is the case e.g. for
  //! the solution of the algebraic system.
  virtual void SetDataPointer(Double * ptr) 
  {Error("BaseStoreSol::SetDataPointer() not implemented here",__FILE__,__LINE__);}
  
  //!
  virtual void CopyFromDataPointer(Double * ptr)
  {Error("BaseStoreSol::CopyFromDataPointer() not implemented here",__FILE__,__LINE__);}  
  
  //!
  virtual void GetDataPointer(Double* &ptr)
  {Error("BaseStoreSol::GetDataPointer() not implemented here",__FILE__,__LINE__);}
  
  //! This function is only preliminary, until there is a REAL
  //! interfacew between CFS++ and OLAS
  virtual Double* GetDoublePointer() = 0;
  
  /////////////////////////////////////////
  // Transformation routines for mapping //
  // pde to mesh solution and vice versa //
  /////////////////////////////////////////
  //!
  virtual void GetElemSolutionAsMatrix(CFSMatrix & elemSol, Vector<Integer> & connect) const = 0;

  //!
  virtual void TransformNodeSolution(BaseStoreSol & transformedSolution,
				     const std::vector<Integer> & mapping,
				     Grid * ptGrid,
				     const Integer level) const = 0;
  
  //!
  virtual void TransformElemSolution(BaseStoreSol & transformedSolution,
				     const std::vector<Elem*> & elems,
				     Grid * ptGrid,
				     const Integer level) const = 0;
  
  //! maps the local element solution to the global mesh solution
  /*!
    \param MeshSol (output) Solution vector referring to mesh node numbers    
    \param PDESol (input) Solution vector
    \param Elems (input) Vector of subdomains to which PDESol belongs to
     */
  virtual void TransformElemSolution(BaseStoreSol & MeshSol, 
				     const std::vector<std::string> & SD,
				     Grid * ptGrid,
				     const Integer level) const = 0;

  
  //! maps the local node solution to the coupling nodes
  virtual void NodeSolutionToCoupling(BaseStoreSol & couplingSol,
				      const std::vector<Integer>& nodeNumbers,
				      const std::vector<Integer> & mapping) const  = 0;
  
  //! maps the local element solution to the coupling nodes
  virtual void ElemSolutionToCoupling(BaseStoreSol & couplingSol,
				      const std::vector<Elem*>& elements,
				      const CFSVector & elemSol) const = 0;

  ////////////////////////////////////////
  // Operators for dealing with vectors //
  // and matrices                       //
  ////////////////////////////////////////

  
#define DEFINE_BASESTORESOL_FCT(TYPE)							\
  virtual void Init(const TYPE val)							\
  {Error("BaseStoreSol::Init() not implemented here",__FILE__,__LINE__);}		\
  virtual void Get(const Integer nodeNr, TYPE & ret) const                              \
  {Error("BaseStoreSol::Get() not implemented here",__FILE__,__LINE__);}                \
  virtual void Get(const SolutionType type, const Integer nodeNr, const Integer dof, TYPE & ret) const	\
  {Error("BaseStoreSol::Get not implemented here",__FILE__,__LINE__);}			\
  virtual void Set(const SolutionType type, const Integer nodeNr, const Integer dof, const TYPE val) const	\
  {Error("BaseStoreSol::Set not implemented here",__FILE__,__LINE__);}                  \
  virtual void Add(const SolutionType type, const Integer nodeNr, const Integer dof, const TYPE val) const         \
  {Error("BaseStoreSol::Add() not implemented here",__FILE__,__LINE__);}                \
   virtual void CopyFromDataPointer(TYPE * ptr)                                         \
  {Error("BaseStoreSol::CopyFromDataPointer() not implemented here",__FILE__,__LINE__);}\
  virtual void SetDataPointer(TYPE * ptr)                                               \
  {Error("BaseStoreSol::SetDataPointer() not implemented here",__FILE__,__LINE__);}     \
   virtual void GetDataPointer(TYPE* &ptr)                                              \
  {Error("BaseStoreSol::GetDataPointer() not implemented here",__FILE__,__LINE__);}   
  
DEFINE_BASESTORESOL_FCT(Complex)

protected:

  //! number of nodes
  Integer numNodes_;

  //! number of different solutions
  Integer numSolutions_;
  
  //! total number of entries
  //! (= numNodes * sumOf(solDofs_[i]))
  Integer length_;
  
  //! contains the types of the solution quantities stored 
  //! in the object (ref. enum SolutionType)
  //std::vector<SolutionType> solTypes_;
  std::map<SolutionType,Integer> solTypes_;


  //! stores the relative position of one
  //! result w.r.t. to the beginning
  //! (depends on the number of dofs
  //! of all previous results)
  std::map<SolutionType,Integer> solOffset_;

  //! contains the number of degrees of freedom for
  //! each quantity
  std::map<SolutionType,Integer> solDofs_;
  //std::vector<Integer> solDofs_;
  
  //! total number of dofs.
  //! The total size of the vector calculates as
  //! totalDofs * numNodes;
  Integer totalDofs_;
  

  //!
  void ResizeData(){};
};

} //end of namespace



#endif
