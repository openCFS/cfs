#ifndef FILE_BASECOUPLEDPDE_2003
#define FILE_BASECOUPLEDPDE_2003


#include <list>
#include <Domain/bcs.hh>
#include <DataInOut/filetype.hh>
#include <DataInOut/writeresults.hh>
#include <PDE/basepde.hh>

namespace CoupledField
{

 
class TimeErrorEstimator;
class PDECoupling;

template<Integer dim>
class SpaceErrorEstimator;

 //! Base class for coupling between different PDEs
  /*! Class BaseCoupledPDE is base class for coupled PDE problems
   */
class BaseCoupledPDE
{
public:

  //! Constructor
  /*!    
    \param PDEs vector with pointers to pdes
    \param aptBCs pointer to boundary condition object
    \param aInFile pointer to class FileType. input data.
    \param aOutFile  pointer to class WriteResults. output data.
    \param aTimeFunc pointer to class TimeFunc
  */
  BaseCoupledPDE(std::vector<BasePDE*> & PDEs,
		 std::vector<PDECoupling*> & Couplings,
		 Grid *aptgrid, 
		 BCs *aptBCs, 
		 FileType *aInFile, 
		 WriteResults *aOutFile); 

  //! Deconstructor
  virtual ~BaseCoupledPDE();

 
  //! define algebraic system solver Parameters
  virtual void SetSolverParameters();
  
  //! calculates coupling interfaces
  virtual void InitCoupling(Integer level)=0;
  
  //! Solve static step
  virtual void SolveStepStatic(const Integer level)=0;
  
  //! solve transient step
  virtual void SolveStepTrans(const Integer kstep, const Double asteptime, const Integer level, 
			      const Boolean updatesysmat)=0;
  
  //! write results in file
  virtual void WriteResultsInFile()=0;

  
  //------------------------ get functions--------------------------------------------------------

  //! return pointer to vector with subdomains, on which we calculate the PDE
  //virtual std::vector<std::string> * getSDsPDE()
  //{ 
  //  return &subdoms_;
  // }

  
  //! returns the local PDE number of an array of nodes
  /*!
    \param PDENodes (output) Vector of PDE node numbers
    \param MeshNodes (input) Vector of mesh node numbers
    \param Mesh2PDENode (input) Vector assigning PDE to mesh node numbers
  */
  virtual void Mesh2PDENode(Vector<Integer> & PDENodes, 
			    const Vector<Integer> & MeshNodes,
			    const std::vector<Integer> & Mesh2PDENode);
  
  //! returns the local global Mesh node numbers of an array of nodes
  /*!
    \param MeshNodes (output) Vector of mesh node numbers    
    \param PDENodes (input) Vector of PDE node numbers
    \param PDE2MeshNode (input) Vector assigning mesh to PDE  node numbers
  */
  virtual void PDE2MeshNode(Vector<Integer> & MeshNodes, 
			    const Vector<Integer> & PDENodes,
			    const std::vector<Integer> & PDE2MeshNode);

protected:
  
  //! maps the local node solution to the global solution
  /*!
    \param MeshSol (output) Solution vector referring to mesh node numbers    
    \param PDESol (input) Solution vector referring to PDE node numbers
    \param nrDofs (input) Dimension of solution values
  */
  virtual void TransformNodeSolution(Vector<Double> & MeshSol, 
				     const Vector<Double> & PDESol,
				     const std::vector<Integer> & PDE2MeshNode,
				     const Integer nrDofs);

  //! maps the local element solution to the global mesh solution
  /*!
    \param MeshSol (output) Solution vector referring to mesh node numbers    
    \param PDESol (input) Solution vector referring to PDE node numbers
    \param Elems (input) Vector of Elements to which the PDESol belongs to (same ordering!)
    \param nrDofs (input) Dimension of solution values
  */
  virtual void TransformElemSolution(Vector<Double> & MeshSol, 
				     const Vector<Double> & PDESol, 
				     const std::vector<Elem*> & Elems,
				     const Integer nrDofs); 

  //! maps the local element solution to the global mesh solution
  /*!
    \param MeshSol (output) Solution vector referring to mesh node numbers    
    \param PDESol (input) Solution vector
    \param Elems (input) Vector of subdomains to which PDESol belongs to
    \param nrDofs (input) Dimension of solution values
  */
  virtual void TransformElemSolution(Vector<Double> & MeshSol, 
				     const Vector<Double> & PDESol, 
				     const std::vector<std::string> & SD,
				     const Integer nrDofs);


  //! assign local PDE node numbers to subdomains
  /*!
    \param Mesh2PDENode (output) Vector assigning mesh to PDE node numbers
    \param PDE2MeshNode (output) Vector assigning PDE to mesh node numbers
    \param subdoms (input) Vector of subdomains which are to be mapped
  */
  virtual void AssignPDENodeNumbers(std::vector<Integer> & Mesh2PDENode,
			    std::vector<Integer> & PDE2MeshNode,
			    const std::vector<std::string> & subdoms);

  //! assign local PDE node numbers to a list of elements (e.g. interface) 
  /*!
    \param Mesh2PDENode (output) Vector assigning mesh to PDE node numbers
    \param PDE2MeshNode (output) Vector assigning PDE to mesh node numbers
    \param subdoms (input) Vector of elements which are to be mapped
  */
  virtual void AssignPDENodeNumbers(std::vector<Integer> & Mesh2PDENode,
			    std::vector<Integer> & PDE2MeshNode,
			    const std::vector<Elem*> &Elements );
  

  // general PDE parameters
  AnalysisType analysistype_;           //!< type of analysis
  std::string pdename_;                 //!< type of PDE (set in the derived classes)
  std::vector<BasePDE *> PDEs_;         //!< list of belonging PDEs
  std::vector<PDECoupling*> Couplings_;   //!< vector of coupling objects
  Integer NumPDEs_;                     //!< number of PDEs 
  std::vector<std::string> *subdoms_;   //!< subdomains belonging to CoupledPDE
  Integer actlevel_;                    //!< current level (for multigrid)

  // coupling elements / nodes
  std::vector<Elem*> * CoupleElems_;    //!< vector of coupling elements 
  std::list<Integer> * CoupleNodes_;  //!< vector of coupling nodes

  Integer NumPDENodes_;   //!< number of nodes in subdomains
  Integer NumElems_;      //!< number of elements in subdomains 

  // pointers to objects
  Grid * ptgrid_;           //!< pointer to Grid
  BCs *ptBCs_;              //!< pointer to Boundary Condition  Object
  FileType * InFile_;       //!< pointer tio input file
  WriteResults * OutFile_;  //!< pointer to output file

  // Assignment MeshNodeNumers <-> PDENodeNumbers
  std::vector<Integer> Mesh2PDENode_;  //!< array containing PDE (=local) node numbers
  std::vector<Integer> PDE2MeshNode_;  //!< array containing Mesh (=global) node numbers 
 
  //!solver parameters
  Integer maxnumiter_;    //!< maximum of iterations (for iterative solver)
  Double  eps_;           //!< accuracy

  
};

} // end of namespace

#endif
