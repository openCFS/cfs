// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef IDBC_HANDLER_PENALTY_HH
#define IDBC_HANDLER_PENALTY_HH

#include <iostream>
#include <set>
#include <utility>

#include "General/defs.hh"
#include "General/environment.hh"
#include "MatVec/basematrix.hh"
#include "baseidbchandler.hh"
#include "def_expl_templ_inst.hh"

namespace CoupledField {


  // forward declarations
class BaseVector;
  struct BaseEntryManipulator;


  //! The %IDBC_HandlerPenalty deals with degrees of freedom that are fixed
  //! by <b>I</b>nhomogeneous <b>D</b>irichlet <b>B</b>oundary
  //! <b>C</b>onditions by the so called penalty approach. The latter works
  //! in the following fashion.
  //! In order to deal with the Dirichlet boundary conditions in the penalty
  //! approach we first set up the linear system \f$A^Nx=b^N\f$ in the fashion
  //! we would use if the nodes with Dirichlet boundary conditions were in
  //! fact nodes with homogeneous Neumann boundary conditions. We then
  //! transform this system into a new system \f$Ax=b\f$. We demonstrate this
  //! for a \f$3\times 3\f$ system, where the node corresponding to \f$x_1\f$
  //! is a Dirichlet node.
  //! \f[
  //! \underbrace{
  //! \left(\begin{array}{ccc}
  //! \ast & \ast & \ast \\
  //! \ast & \ast & \ast \\
  //! \ast & \ast & \ast
  //! \end{array}\right)
  //! }_{A^N}
  //! \enspace,\quad
  //! \underbrace{
  //! \left(\begin{array}{c} \ast \\ \ast \\ \ast \end{array}\right)
  //! }_{b^N}
  //! \qquad\longrightarrow\qquad
  //! \underbrace{
  //! \left(\begin{array}{ccc}
  //! \diamondsuit & \ast & \ast \\
  //! \ast & \ast & \ast \\
  //! \ast & \ast & \ast
  //! \end{array}\right)
  //! }_{A}
  //! \enspace,\quad
  //! \underbrace{
  //! \left(\begin{array}{c} \otimes \\ \ast \\ \ast \end{array}\right)
  //! }_{b}
  //! \f]
  //! where \f$\ast\f$ represents an arbitrary value, \f$\diamondsuit\f$
  //! represents a sufficiently large positive value (determined in %OLAS as
  //! product of the largest magnitude of a diagonal matrix entry and the
  //! value of the #PENALTY macro) and
  //! \f[
  //! \otimes = d_1 \cdot\,\diamondsuit\enspace,
  //! \f]
  //! with \f$d_1\f$ being the specified Dirichlet boundary value at the node.
  template<typename T>
  class IDBC_HandlerPenalty : public BaseIDBC_Handler {

  public:

    //! Adapt system matrix

    //! This method can be used by those approaches that require a
    //! modification of the system matrix in order to incorporate
    //! inhomogeneous Dirichlet boundary conditions into the linear system,
    //! like e.g. the penalty approach.
    void AdaptSystemMatrix( BaseMatrix &sysMat );

    //! Incorporate inhomogeneous Dirichlet BCs into right hand side

    //! Call this method in order to update an exsisting right hand side
    //! such that it incorporates the modification resulting from
    //! treating the dofs for inhomogeneous Dirichlet boundary conditions
    //! by the penalty approach.
    //! \note The method over-writes existing entries in the positions of
    //!       the Dirichlet dofs by the modified boundary values. It is the
    //!       caller's responsibility to ensure that this does not destroy
    //!       any valuable information.
    //! \param rhs vector with right-hand side entries
    inline void AddIDBCToRHS( BaseVector *rhs );

    //! Set value for a Dirichlet boundary condition

    //! This method can be used to set the value of a degree of freedom that
    //! is fixed by a homogeneous Dirichlet boundary condition.
    //! \param pdeID identifier of PDE; this is only used in the case of an
    //!              SBM_System in order to identify the sub-vector in which
    //!              the value must be stored
    //! \param eqnNo equation number for the degree of freedom whose value
    //!              should be set
    //! \param val   inhomogeneous Dirichlet value
    void SetIDBC( PdeIdType pdeID, UInt eqnNo, const T &val );

    //! Set fixed dofs to specified Dirichlet boundary values

    //! This method replaces the values of all fixed degrees of freedom in the
    //! specified input vector by new values. These new values are taken to
    //! be the values specified via the inhomogeneous Dirichlet boundary
    //! condition that fixes the respective degrre of freedom.
    void SetDofsToIDBC( BaseVector *vec );


    // =======================================================================
    // CONSTRUCTION, DESTRUCTION and CLEARING
    // =======================================================================

    //! \name Constructors, destructor and initialisation methods

    //@{

    //! Admissable constructor for the case of StdMatrices

    //! This constructor is intended to by used in the case of a linear
    //! system specified with the help of StdMatrices.
    //! \param numIDBC   number of degrees of freedom fixed by inhomogeneous
    //!                  Dirichlet boundary conditions
    //! \param blockSize block size of matrix entries. This is > 1 in the
    //!                  case that dof blocking is employed and e.g. the three
    //!                  displacements at a single node in a mechanics problem
    //!                  are treated as one entity in the linear system
    //!                  (Note only required for setting up the assembler_
    //!                  object!)
    IDBC_HandlerPenalty( UInt numIDBC, UInt blockSize );

    //! Admissable constructor for the case of SBM_Matrices

    //! This constructor is intended to by used in the case of a linear
    //! system specified with the help of SBM_Matrices.
    //! \param numIDBC   number of degrees of freedom fixed by inhomogeneous
    //!                  Dirichlet boundary conditions
    //! \param numPDEs   number of PDEs in the SBM_System
    //! \param bcOffsets array for initialising the bcOffsets_ attribute
    IDBC_HandlerPenalty( UInt numIDBC, UInt numPDEs, UInt *bcOffsets );

    //! Destructor
    ~IDBC_HandlerPenalty();

    //! Re-set vector of Dirichlet values

    //! Calling this method deletes all information stored by the object
    //! internally on the values and degrees of freedom that are fixed
    //! by inhomogeneous Dirichlet boundary conditions.
    void InitDirichletValues();

    //@}


    // =======================================================================
    // STUBS
    // =======================================================================

    //! \name Stubs
    //! The following method are not required for this derivation of the
    //! BaseIDBC_Handler class. For simplicity we decided to make the methods
    //! stubs, i.e. if they are called, they will simply return without
    //! performing any real action.

    //@{

    //! Combine different FE matrices into a single system matrix
    void BuiltSystemMatrix( const std::map<FEMatrixType, Double> &factors ) {
    }

    //! Remove inhomogeneous Dirichlet BCs from right hand side
    inline void RemoveIDBCFromRHS( BaseVector *rhs ) {
    }

    //! Add weight of coupling between a fixed and a free dof into matrix
    void AddWeightFixedToFree( FEMatrixType matID,
                               PdeIdType pdeID1,
                               PdeIdType pdeID2,
                               UInt rowInd,
                               UInt colInd,
                               const T& value ) {
    };

    //! Set weight of coupling between a fixed and a free dof into matrix
    void SetWeightFixedToFree( FEMatrixType matID,
                               PdeIdType pdeID1,
                               PdeIdType pdeID2,
                               UInt rowInd,
                               UInt colInd,
                               const T& value ) {
    };

    //! Get weight of coupling between a fixed and a free dof from matrix
    void GetWeightFixedToFree( FEMatrixType matID, PdeIdType pdeID1,
                               PdeIdType pdeID2, UInt rowInd, UInt colInd,
                               T & value ) const {
    };
    
    //! Set the value of all coupling weights of a free dof to its fixed ones
    void SetRowWeights( FEMatrixType matID, PdeIdType pdeID, UInt rowInd,
                        Double realPart, Double imagPart = 0.0 ) {
    }
    
    
    //! Set the value of all coupling weights of a fixed dof to its free ones
    void SetColWeights( FEMatrixType matID, PdeIdType pdeID,UInt colInd,
                        Double realPart, Double imagPart = 0.0 ) {
    };

    

    //! Re-set specified internal matrix to zero
    void InitMatrix( FEMatrixType matrixID ) {
    };

    /** Our penalty term */
    Double GetPenalty() const {
      return penaltyTerm_;
    }
    
    void Dump()
    {
      std::map<PdeIdType, std::map<Integer, UInt> >::iterator i;
      for(i = bcIndices_.begin(); i != bcIndices_.end(); i++)
      {
        std::cout << " pde_type = " << i->first << std::endl;
        std::map<Integer, UInt>::iterator j;
        for(j = i->second.begin(); j != i->second.end(); j++)
          std::cout << " " << j->second;
        std::cout << std::endl;
      }
    }

    /** Get for a equation number the penalty dirichlet value.
     * @param equation if invalid we return false
     * @param dirichlet_value because it can be double/complex as paramter. Unset if return false
     * @return true if equation is valid, then dirichlet_value ist set. */
    bool GetIDBC(PdeIdType pde_type, int equation, T& dirichlet_value) 
    {
      // search for the equation and check 
      std::map<Integer, UInt>::iterator iter = bcIndices_[pde_type].find(equation);
      if(iter == bcIndices_[pde_type].end()) return false;

      // set dirichlet value
      dirichletValue_->GetEntry(iter->second, dirichlet_value);

      return true;
    }
    
    
    //@}


  private:

    //! Flag indicating whether we are faced with SBM_Matrices/Vectors
    bool sbmCase_;

    //! Our private assembling object
    BaseEntryManipulator *assembler_;

    //! Number of inhomogeneous Dirichlet boundary conditions

    //! This attribute stores the number of degrees of freedom that are
    //! fixed by inhomogeneous Dirichlet boundary conditions. Note that in
    //! the case of dof-blocking, i.e. when dofs at a node a treated as one
    //! algebraic entity, this number may be smaller than the number of
    //! CFS++ equation numbers for inhomogeneous Dirichlet boundary
    //! conditions. This number is the total number, i.e. the sum over all
    //! PDEs.
    UInt numIDBC_;

    //! Array used for storing equation numbers for fixed dofs

    //! This array is used to store for each degree of freedom that is
    //! fixed by an inhomogeneous Dirichlet boundary condition the related
    //! equation number assigned to this dof (scalar case) or the block of
    //! dofs the dof belongs to (block case).
    UInt *dirichletEQN_;

    //! Next index for dirichlet condition set
    UInt nextIndex_;

    //! Map for mapping eqn number to position in dirichletValue_ vector
    std::map<PdeIdType, std::map<Integer, UInt> > bcIndices_;

    //! Array for storing inhomogeneous Dirichlet values
    BaseVector *dirichletValue_;

    //! Value of the penalty term in the penalty approach
    Double penaltyTerm_;

    //! Attribute storing run-time information on class template
    BaseMatrix::EntryType eType_;

    //! Array storing index pointers into internal arrays

    //! This array stores for each PDE the index (minus 1) at which the
    //! information for this PDE starts in the dirichletComponent_ and
    //! dirichletEQN_ arrays. Since storage is contiguous we can also use this
    //! information to compute the number of IDBCs for an individual PDE.
    //! The array is enlarged by one entry for this purpose.
    //! \note This array is only used in the SBM case. In the case of an
    //!       IDBC handler for a standard system, this will be a NULL
    //!       pointer.
    UInt *offset_;

    //! Method for formatted output of internal arrays

    //! This method will print the information on the inhomogeneous
    //! Dirichlet boundary conditions (equation numbers, components,
    //! values) stored in the internal arrays to the specified stream.
    //! It is primarily intended for debugging purposes.
    void PrintInternalState( std::ostream &os );

  };

}

#ifndef EXPLICIT_TEMPLATE_INSTANTIATION
//#include "idbchandlerpenalty.cc"
#endif

#endif
