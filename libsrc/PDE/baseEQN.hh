#ifndef FILE_BASEEQN_2004
#define FILE_BASEEQN_2004

#include "General/environment.hh"
#include "Domain/grid.hh"
#include "Utils/StdVector.hh"

namespace CoupledField
{

  //! a base class for equation data handling

  class BaseEQN
  {
  public:

    //! Constructor
    BaseEQN(Grid * aptGrid,
            StdVector<RegionIdType>& asubdoms, 
            UInt dofsPerNode);

    //! Destructor
    virtual ~BaseEQN();
  
    //! Return true if eqns are assigned
    inline Boolean IsInitialized() {return isInitialized_;}

    //! Return the number of equations
    //! (= number of nodes * total number Dofs <br>
    //! - hom. Dirichlet nodes <br>
    //! - constraint nodes)
    inline UInt GetNumEQNs() const {return numEqns_;}

    //! Return the number of real equations
    //! (= number of nodes * total number Dofs
    //! - hom. Dirichlet nodes
    //! - inhom. Dirichlet nodes
    //! - constraint nodes)
    inline UInt GetNumRealEQNs() const {
      return numRealEqns_;
    }
  
    //! Return number of degree of freedoms per node
    inline UInt GetNumDofs() const {return dofsPerNode_;}

    //! Return number of dofs per equation
    inline UInt GetNumDofsPerEQN() const {return dofsPerEQN_;}
  
    //! Return number of build in Dirichlet values

    //! Returns the number of Dirichlet values that were assigned
    //! a 0 and therefore have not to set explicitly to 
    //! the algebraic system
    inline UInt GetNumBuildInDirichletEQNs() const
    {return numBuildInDirichletEQNs_;}

    //! Returns true, if nodal mapping is applied
    inline Boolean IsNodalMapped() const {return isNodalMapped_;}

    //! Returns true, if block mapping is applied
    inline Boolean IsBlockMapped() const {return isBlockMapped_;}

    //! Set the node numbers and dofs 
    //! with homogeneous Dirchlet BC
    void SetHomoDirichletBCs(const StdVector<std::string> & nodeNrs,
                             const StdVector<std::string> & dofs);

    //! Set the node numbers and dof for inhomogeneous Dirchlet BC
    void SetInhomDirichletBCs( const StdVector<std::string> & nodeNrs,
                               const StdVector<std::string> & dofs );

    //! Set the node numbers and dofs which are
    //! slave nodes w.r.t. to constraints
    void SetConstraints(const StdVector<UInt> & slaveNodeNrs,
                        const StdVector<UInt> & masterodeNrs,
                        const StdVector<std::string> & dofs);
                              
    //! Calculate the mapping after Dirichlet and
    //! constraint nodes were set
    virtual void CalcMapping() = 0;

    //! Calculate only the mapping from global to
    //! local node numbers and back.
    //! This is a reduced form of "CalcMapping()"
    virtual void CalcMpcciMapping(){};

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
    virtual void ReorderMapping(Integer *order) = 0;

    //! Print the mapping nodes <->EQNs
    virtual void Print(std::ostream & out) const = 0;
  
  protected:

    //! Status flag; TRUE indicates state after computation of mappings
    Boolean isInitialized_;
  
    //! Flag for indicating blockwise numbering
    Boolean isBlockMapped_;

    //! Flag for indicating nodal numbering vs. edge numbering
    Boolean isNodalMapped_;

    //! Flag for turning on ordering of equation numbers

    //! If this flag is TRUE, the equation numbers are ordered in the following
    //! fashion. The set of highest equation numbers belongs to "unknowns" whose
    //! values are fixed by inhomogeneous Dirichlet boundary conditions.
    Boolean sortEQNs_;

    //! Pointer to Grid
    Grid * ptGrid_;  

    //! Number of nodes in PDE
    UInt numPDENodes_;

    //! Number of dofs per node
    UInt dofsPerNode_;
 
    //! Number of dofs per eqn
    UInt dofsPerEQN_;

    //! Number of equations in PDE
    UInt numEqns_;

    //! Number of "real" equations in PDE

    //! This attribute stores the number of "real" equations in the linear
    //! system associated with the underlying PDE. "Real" here means that the
    //! value of the unknown that corresponds to an equation is not fixed by
    //! a Dirichlet boundary condition.
    //! \note This value is only computed, when sortEqns_ is TRUE, otherwise
    //!       we store a -1.
    UInt numRealEqns_;
  
    //! Number of Dirichlet values
    //! which have been eliminated
    //! by equation class
    UInt numBuildInDirichletEQNs_;

    //! Vector containing subdomain names
    //! of according PDE
    StdVector<RegionIdType> subdoms_;  

    //! Vector with indices of nodes with hom. Dirichlet boundary conditions
    StdVector<UInt> homoDirichletNodes_;

    //! Vector containing dofs associated
    //! with hom. Dirichlet nodes
    StdVector<UInt> homoDirichletDofs_;

    //! Vector with indices of nodes with inhom. Dirichlet boundary conditions
    StdVector<UInt> inhomDirichletNodes_;

    //! Vector containing dofs associated with inhom. Dirichlet nodes
    StdVector<UInt> inhomDirichletDofs_;
  
    //! Vector containing constraint master nodes
    StdVector<UInt> constraintMasterNodes_;

    //! Vector containing according slave nodes
    StdVector<UInt> constraintSlaveNodes_;
  
    //! Vector containing dofs associated
    //! with slave nodes
    StdVector<UInt> constraintDofs_;

    //! Mapping of string-dof to integer (1-based)
    UInt GetBCDof(const std::string dofString) const;

  protected:

    //! Default constructor is disallowed
    BaseEQN() {}; 


  };

}

#endif // FILE_BASEEQN
