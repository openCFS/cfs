#ifndef FILE_SCALARNODEEQN_2004
#define FILE_SCALARNODEEQN_2004

#include "nodeEQN.hh"

namespace CoupledField {

  //! This class is used to perform the handling of equation numbers in the
  //! case of scalar PDEs, i.e. PDEs which only have one degree of freedom per
  //! node. In this case every unknown receives a unique equation number.
  //! Unknowns whose value is fixed by homogeneous Dirichlet boundary
  //! conditions are eliminated. Unknowns given by inhomogeneous Dirichlet
  //! boundary conditions are treated with the penalty approach. The matrix
  //! representing the linear system associated with the discretised PDE will
  //! have scalar double or complex entries.
  class ScalarNodeEQN : public NodeEQN {

  public:
  
    //! Constructor
    ScalarNodeEQN( Grid * aptGrid, 
		   BCs * aptBCs,
		   StdVector<std::string>& asubdoms, 
		   Integer actlevel, 
		   Integer dofsPerNode );
  
    //! Destructor
    virtual ~ScalarNodeEQN();
  
    //! Calculate the mapping after Dirichlet and
    //! constraint nodes were set
    void CalcMapping();

    //! Calculate only the mapping from global to
    //! local node numbers and back.
    //! This is a reduced form of "CalcMapping()"
    void CalcMpcciMapping();
  
    //! Print the mapping nodes <-> EQNs
    void Print(std::ostream & out) const;
  
    // -----------------------------------------------------------------------
    // Functions for mapping node numbers <-> EQN numbers
    // -----------------------------------------------------------------------
  
    //! Map vector of equation numbers to 
    //! positions in global solution vector
    void EQN2SolVectorPos(const StdVector<Integer> &eqnNr, 
			  StdVector<Integer> &pos) const;

    //! Map node number and dof to according equation number
    void Node2EQN(const Integer nodeNr, 
		  const Integer dof,
		  Integer & eqnNr,
		  Integer & eqnDof) const;

    //! Map node number to according equation number(s)
    void Node2EQN(const Integer nodeNr, StdVector<Integer> &eqns) const;
  
    //! Map vector of node numbers to according
    //! vector of equiation numbers
    void Node2EQN(const StdVector<Integer> &nodeNr,
		  StdVector<Integer> &eqnNr) const;

  private:

    //! Default constructor is disallowed
    ScalarNodeEQN() {};

    //! Mapping global node number->EQN
    StdVector<Integer> pdeNode2EQN_;

  };

}  // end of namespace

#endif // FILE_SCALARNODEEQN
