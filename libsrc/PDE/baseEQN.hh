#ifndef FILE_BASEEQN_2004
#define FILE_BASEEQN_2004

#include "General/environment.hh"
#include "Domain/grid.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {

  //! a base class for equation data handling
  class BaseEQN {

  public:

    //! Constructor
    BaseEQN( Grid *aptGrid, StdVector<RegionIdType>& asubdoms,
             UInt dofsPerNode, bool sortEQNs );

    //! Destructor
    virtual ~BaseEQN();
  
    //! Return true if eqns are assigned
    inline Boolean IsInitialized() {return isInitialized_;}

    //! Return the total number of equations

    //! This method can be used to query the total number of equations. This
    //! value coincides with the highest equation number that was assigned
    //! to a degree of freedom. It represents the sum of free degrees of
    //! freedom and those fixed by inhomogeneous Dirichlet boundary
    //! conditions. Dofs fixed by homogeneous Dirichlet boundary conditions
    //! or constraints do not enter this value.
    inline UInt GetNumEQNs() const {
      return numEqns_;
    }

    //! Return the equation number of the last unfixed degree of freedom

    //! This method returns the equation number of the last free degree of
    //! freedom. We consider a degree of freedom to be free, if it is not
    //! fixed by either
    //! - homogeneous Dirichlet boundary conditions
    //! - inhomogeneous Dirichlet boundary conditions
    //! - constraints (master - slave dof relations)
    //! The return value actually depends on the status of the sortEQNs_ flag.
    //! If the latter is true, then we return numRealEqns_. If it is false,
    //! then inhomogeneous Dirichlet boundary conditions are treated by the
    //! penalty method and must be considered free dofs as well. Thus, in this
    //! case we return numEqns_.
    inline UInt GetNumLastFreeDof() {
      ENTER_FCN( "BaseEQN::GetNumLastFreeDof" );
      return ( ( sortEQNs_ == TRUE ) ? numRealEqns_ : numEqns_ );
    }
  
    //! Return number of degree of freedoms per node
    inline UInt GetNumDofs() const {return dofsPerNode_;}

    //! Return number of dofs per equation
    inline UInt GetNumDofsPerEQN() const {return dofsPerEQN_;}
  
    //! Return number of dropped degrees of freedom

    //! This method returns the number of degrees of freedom that were
    //! dropped, i.e. they were either assigned a zero equation number,
    //! since they represent homogeneous Dirichlet boundary conditions,
    //! or they were multiply defined, or their nodes do not belong to
    //! the associated PDE.
    inline UInt GetNumDroppedDofs() const {
      return numDroppedDofs_;
    }

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

    //! If this flag is TRUE, the equation numbers are ordered such that
    //! the set of highest equation numbers belongs to "unknowns" whose
    //! values are fixed by inhomogeneous Dirichlet boundary conditions.
    Boolean sortEQNs_;

    //! Pointer to Grid
    Grid *ptGrid_;  

    //! Number of nodes in PDE
    UInt numPDENodes_;

    //! Number of dofs per node
    UInt dofsPerNode_;
 
    //! Number of dofs per eqn
    UInt dofsPerEQN_;

    //! Number of equations in PDE
    UInt numEqns_;

    //! Number of "real" equations in PDE

    //! This attribute stores the number of 'real' equation numbers. By the
    //! latter we understand equation numbers for degrees of freedom that
    //! are not fixed by either
    //! - homogeneous Dirichlet boundary conditions
    //! - inhomogeneous Dirichlet boundary conditions
    //! - constraints (master - slave dof relations)
    UInt numRealEqns_;
  
    //! Number of Dirichlet values
    //! which have been eliminated
    //! by equation class
    UInt numDroppedDofs_;

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

  protected:

    //! Default constructor is disallowed
    BaseEQN() {}; 

  };

#ifdef DOXYGEN_DETAILED_DOC

  // =========================================================================
  //     Detailed description of the class 
  // =========================================================================

  //! \class BaseEQN
  //! 
  //! \purpose
  //!
  //! This class is the base class which defines the interface of all
  //! objects that deal with the assignment of equation numbers to degrees
  //! of freedom. The following gives an overview on the meaning of the
  //! different equation numbers
  //! <center><img src="../AddDoc/splitting.png"></center>
  //!
  //! \collab 
  //!
  //! \implement 
  //!
  //! \status In use
  //!
  //! \unused 
  //!
  //! \improve
  //! 

#endif

}

#endif // FILE_BASEEQN
