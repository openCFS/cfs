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
                   StdVector<RegionIdType>& asubdoms, 
                   UInt dofsPerNode );
  
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
  
    //! Map node number and dof to according equation number
    void Node2EQN(const UInt nodeNr, 
                  const UInt dof,
                  Integer & eqnNr,
                  UInt & eqnDof) const;

    //! Map node number to according equation number(s)
    void Node2EQN(const UInt nodeNr, StdVector<Integer> &eqns) const;
  
    //! Map vector of node numbers to according
    //! vector of equiation numbers
    void Node2EQN(const StdVector<UInt> &nodeNr,
                  StdVector<Integer> &eqnNr) const;

    //! Maps the equation numbers according to the reordering

    //! If the algebraic system has re-ordered the equations, e.g. to improve
    //! the performance of a sparse direct solver, this method must be used to
    //! incorporate the resulting change into the dof to eqn map.
    //! \param order permutation vector containing the re-ordering, i.e.
    //!              order[k-1] is the new equation number for the dof that
    //!              was previously associated to the k-th equation number.
    //! \note It is safe to pass a NULL pointer to the method. In this case
    //!       the method assumes that no re-ordering was performed by the
    //!       algebraic system.
    void ReorderMapping(Integer *order);

  private:

    //! Default constructor is disallowed
    ScalarNodeEQN() {};

    //! Mapping global node number->EQN
    StdVector<Integer> pdeNode2EQN_;

  };

}  // end of namespace

#endif // FILE_SCALARNODEEQN
