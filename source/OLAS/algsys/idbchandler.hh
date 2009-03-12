// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef IDBC_HANDLER_HH
#define IDBC_HANDLER_HH

#include <set>

#include <def_expl_templ_inst.hh>

#include "OLAS/algsys/baseidbchandler.hh"


namespace CoupledField {


  // forward declarations
  class BaseGraphManager;


  // ========================================================================
  // IDBC_HANDLER CLASS
  // ========================================================================


  //! The IDBCHandler class is responsible for storing the weights with
  //! which real degrees of freedom couple to those that are fixed by
  //! <b>I</b>nhomogeneous <b>D</b>irichlet <b>B</b>oundary <b>C</b>onditions
  //! and for updating the right hand side of the linear system using these
  //! weights and the prescribed Dirichlet values.
  //! \note The CFS++/OLAS interface arranges equation numbers for degrees of
  //! freedom in such a manner that those degrees of freedom that are fixed
  //! by inhomogeneous Dirichlet boundary conditions start at (limit+1)\n\n
  //! <center><img src="../AddDoc/splitting.png"></center>\n
  //! In order to store the coupling weights between these fixed and free
  //! degrees of freedom as well as the boundary values themselves, we need
  //! one-based indices for the dofs. Thus, we subtract a corresponding offset
  //! (=limit)  from the equation numbers.
  template<typename T>
  class IDBC_Handler : public BaseIDBC_Handler {

  public:


    // =======================================================================
    // CONSTRUCTION, DESTRUCTION and CLEARING
    // =======================================================================

    //! \name Constructors, destructor and initialisation methods

    //@{

    //! Admissable constructor for both the %StdMatrix and %SBM_Matrix case

    //! This is currently the only allowed constructor for the IDBC_Handler
    //! class.
    //! \param usedFEMatrices set containing specifiers for those FE matrices
    //!                       (DAMPING, STIFFNESS, ...) that are required to
    //!                       generate the final SYSTEM matrix
    //! \param graphManager   pointer to a graph manager object, the latter
    //!                       is used to obtain size and sparsity pattern
    //!                       information for the internal auxilliary matrices
    //! \param numPDEs        number of PDEs for which we administrate
    //!                       inhomogeneous Dirichlet boundary conditions
    //!                       (ignored in the %StandardSystem case)
    //! \param sbmCase        flag indicating whether we are confronted with
    //!                       a Standard- or an %SBM_System
    IDBC_Handler( const std::set<FEMatrixType> &usedFEMatrices,
                  BaseGraphManager *graphManager, UInt numPDEs,
                  bool sbmCase );

    //! Destructor
    ~IDBC_Handler();

    //! Re-set specified internal matrix to zero

    //! Calling this method re-sets the specified internal matrix of the
    //! object to its initial state, i.e. Init() is called on that matrix
    //! object.
    void InitMatrix( FEMatrixType matrixID ) {


      // Check for validity of matrix identifier
      mySetIterator it = myMatrices_.find( matrixID );
      if ( it == myMatrices_.end() ) {
        EXCEPTION( "IDBC_Handler::InitMatrix: Received invalid matrix "
                 << "identifier '" << matrixID << "'");
      }
      else {
        auxMat_[ *it ]->Init();
      }

      // Adapt internal status flags
      addIDBCPossible_ = false;
      remIDBCPossible_ = false;
    }

    //! Re-set vector of Dirichlet values

    //! Calling this method deletes all information stored by the object
    //! internally on the values and degrees of freedom that are fixed
    //! by inhomogeneous Dirichlet boundary conditions.
    void InitDirichletValues() {
      vecIDBC_->Init();
    };

    //@}

    // =======================================================================

    //! Combine different FE matrices into a single system matrix

    //! This method must be called to combine the different internal FE
    //! matrices into a single system matrix. The system matrix is constructed
    //! by computing a weighted sum of the internal FE matrices using the
    //! weights supplied in the factors map input parameter.
    void BuiltSystemMatrix( const std::map<FEMatrixType, Double> &factors );

    //! Incorporate inhomogeneous Dirichlet BCs into right hand side

    //! Call this method in order to update an exsisting right hand side
    //! such that it incorporates the modification resulting from eliminating
    //! the dofs for inhomogeneous Dirichlet boundary conditions from the
    //! linear system.
    //! \note It is the caller's responsibility to ensure that the initial
    //!       right-hand side passed to the method does not contain Dirichlet
    //!       values from previous calls, since the method will not eliminate
    //!       them.
    //! \param rhs vector with right-hand side entries
    inline void AddIDBCToRHS( BaseVector *rhs ) {
      if ( addIDBCPossible_ == false ) {
        EXCEPTION( "IDBCHandler::AddIDBCToRHS: Internal error! Refusing to "
                 << "add Dirichlet BCs, since addIDBCPossible_ = false! "
                 << "Did you call BuiltSystemMatrix()?");
      }
      auxMat_[SYSTEM]->MultSub( *vecIDBC_, *rhs );
      remIDBCPossible_ = true;
    }

    //! Remove inhomogeneous Dirichlet BCs from right hand side

    //! Call this method in order to update an exisiting right hand side
    //! such that it incorporates the modification resulting from eliminating
    //! the dofs for inhomogeneous Dirichlet boundary conditions from the
    //! linear system.
    //! \note It is the caller's responsibility to ensure that the initial
    //!       right-hand side passed to the method does not contain Dirichlet
    //!       values from previous calls, since the method cannot eliminate
    //!       them.
    //! \param rhs vector with right-hand side entries
    inline void RemoveIDBCFromRHS( BaseVector *rhs ) {
      if ( remIDBCPossible_ == false ) {
    	EXCEPTION( "IDBCHandler::RemoveIDBCFromRHS: Internal error! "
                 << "Refusing to remove Dirichlet BCs, since "
                 << "remIDBCPossible_ = false! "
                 << "Either FE matrices or Dirichlet values have changed "
                 << "since last call to AddIDBCToRHS()!");
      }
      auxMat_[SYSTEM]->MultAdd( *vecIDBC_, *rhs );
    }

    //! Set value for a Dirichlet boundary condition

    //! This method can be used to set the value of a degree of freedom that
    //! is fixed by a homogeneous Dirichlet boundary condition.
    //! \param pdeID identifier of PDE; this is only used in the case of an
    //!              SBM_System in order to identify the sub-vector in which
    //!              the value must be stored
    //! \param eqnNo equation number for the degree of freedom whose value
    //!              should be set
    //! \param val   inhomogeneous Dirichlet value
    //! \note The bcNum parameter is only supported for compatibility with
    //!       the %SetIDBC() method of the IDBC_HandlerPenalty class. It is
    //!       ignored by the current implementation.
    void SetIDBC( PdeIdType pdeID, UInt eqnNo, const T &val );

    //! Add weight of coupling between a fixed and a free dof into matrix

    //! This method provides an interface to add the weight of the coupling
    //! between a free degree of freedom and a degree of freedom fixed by
    //! an inhomogeneous Dirichlet boundary condition into the desired
    //! Finite Element matrix.
    //! \param matID    specifies the FE matrix type of the auxilliary matrix
    //! \param pdeID1   row index of sub-matrix in the SBM_Matrix case
    //! \param pdeID2   column index of sub-matrix in the SBM_Matrix case
    //! \param rowInd   row index of entry to be set, i.e. the equation
    //!                 number of the free degree of freedom
    //! \param colInd   column index of entry to be set, i.e. the equation
    //!                 number of the fixed degree of freedom; must be the
    //!                 real equation number, %IDBC_Handler transforms this
    //!                 to a one-based index itself.
    //! \param val      value of the weight of the coupling
    void AddWeightFixedToFree( FEMatrixType matID, PdeIdType pdeID1,
                               PdeIdType pdeID2, UInt rowInd, UInt colInd,
                               const T& val );
    
    //! Set weight of coupling between a fixed and a free dof into matrix

    //! This method provides an interface to set the weight of the coupling
    //! between a free degree of freedom and a degree of freedom fixed by
    //! an inhomogeneous Dirichlet boundary condition into the desired
    //! Finite Element matrix.
    //! \param matID    specifies the FE matrix type of the auxilliary matrix
    //! \param pdeID1   row index of sub-matrix in the SBM_Matrix case
    //! \param pdeID2   column index of sub-matrix in the SBM_Matrix case
    //! \param rowInd   row index of entry to be set, i.e. the equation
    //!                 number of the free degree of freedom
    //! \param colInd   column index of entry to be set, i.e. the equation
    //!                 number of the fixed degree of freedom; must be the
    //!                 real equation number, %IDBC_Handler transforms this
    //!                 to a one-based index itself.
    //! \param val      value of the weight of the coupling
    void SetWeightFixedToFree( FEMatrixType matID, PdeIdType pdeID1,
                               PdeIdType pdeID2, UInt rowInd, UInt colInd,
                               const T& val );

    //! Get weight of coupling between a fixed and a free dof from matrix

    //! This method provides an interface to get the weight of the coupling
    //! between a free degree of freedom and a degree of freedom fixed by
    //! an inhomogeneous Dirichlet boundary condition from the desired
    //! Finite Element matrix.
    //! \param matID    specifies the FE matrix type of the auxiliary matrix
    //! \param pdeID1   row index of sub-matrix in the SBM_Matrix case
    //! \param pdeID2   column index of sub-matrix in the SBM_Matrix case
    //! \param rowInd   row index of entry to get, i.e. the equation
    //!                 number of the free degree of freedom
    //! \param colInd   column index of entry to get, i.e. the equation
    //!                 number of the fixed degree of freedom; must be the
    //!                 real equation number, %IDBC_Handler transforms this
    //!                 to a one-based index itself.
    //! \param val      value of the weight of the coupling
    void GetWeightFixedToFree( FEMatrixType matID, PdeIdType pdeID1,
                               PdeIdType pdeID2, UInt rowInd, UInt colInd,
                               T& val ) const;
    
    //! Set the value of all coupling weights of a free dof to its fixed ones

    //! This method provides an interface to set the weights of the coupling
    //! between a free degree of freedom and all its degrees of freedom fixed
    //! by an inhomogeneous Dirichlet boundary condition.
    //! \param matID    specifies the FE matrix type of the auxiliary matrix
    //! \param pdeID    row index of sub-matrix in the SBM_Matrix case
    //! \param rowInd   row index of entry to get, i.e. the equation
    //!                 number of the free degree of freedom
    //! \param val      value of the weight of the coupling
    void SetRowWeights( FEMatrixType matID, PdeIdType pdeID, UInt rowInd,
                        const T& val );

    
    //! Set the value of all coupling weights of a fixed dof to its free ones
    
    //! This method provides an interface to set the weights of the coupling
    //! between a degree of freedom fixed by an inhomogeneous Dirichlet 
    //! boundary  condition and all free dofs which depend on it.
    //! \param matID    specifies the FE matrix type of the auxiliary matrix
    //! \param pdeID    column index of sub-matrix in the SBM_Matrix case
    //! \param colInd   column index of entry to get, i.e. the equation
    //!                 number of the fixed degree of freedom
    //! \param realPart real valued part of the weight of the coupling
    void SetColWeights( FEMatrixType matID, PdeIdType pdeID,UInt colInd,
                        const T& val );

    //! Set fixed dofs to specified Dirichlet boundary values

    //! This method will insert the values specified for fixed degrees of
    //! freedom into the given vector.
    //! \note This is a temporary solution for inserting the Dirichlet
    //!       values into an enlarged solution vector until the StoreSolution
    //!       classes in CFS++ handle these values themselves
    void SetDofsToIDBC( BaseVector *vec );


    // =======================================================================
    // STUBS
    // =======================================================================

    //! \name Stubs
    //! The following method are not required for this derivation of the
    //! BaseIDBC_Handler class. For simplicity we decided to make the methods
    //! stubs, i.e. if they are called, they will simply return without
    //! performing any real action.

    //@{

    //! Adapt system matrix
    void AdaptSystemMatrix( BaseMatrix &sysMat ) {
    }

    //@}


    // =======================================================================
    // METHODS FOR FACTORY CONCEPT
    // =======================================================================

    //! \name Methods required for factory concept
    //! The following methods are required for our implementation of the
    //! factory concept and our idea of compiler-independent template class
    //! instantiation.

  private:

    //! Private shortcut for set iterator
    typedef std::set<FEMatrixType>::iterator mySetIterator;

    //! Default constructor is dis-allowed
    IDBC_Handler() {
      EXCEPTION( "Call to forbidden default constructor of "
               << "IDBC_Handler class" );
    }

    //! Copy constructor is dis-allowed
    IDBC_Handler( const IDBC_Handler &idbcHandler ) {
      EXCEPTION( "Call to forbidden copy constructor of IDBC_Handler class" );
    }

    //! Set containing markers for the potentially required matrices
    std::set<FEMatrixType> myMatrices_;

    //! Array containing links to all internal auxilliary matrices
    BaseMatrix **auxMat_;

    //! %Vector for storing the inhomogeneous Dirichlet values
    BaseVector *vecIDBC_;

    // Entry type of auxilliary matrices

    // We use this attribute to store the entry type of the auxilliary
    // matrices. This actually is then also an enumeration describing
    // the value of the template parameter of the instance of this class.
    // MatrixEntryType eType_;

    //! Array storing the number of nodes fixed by inhom. Dirichlet BCs

    //! This array contains for every sub-system the number of degrees of
    //! of freedom that are fixed by inhomogeneous Dirichlet boundary
    //! conditions. In the StdMatrix case there is only a single system,
    //! so the length of the array is one. In the SBM_Matrix case the length
    //! corresponds to the number of diagonal matrix blocks.
    UInt *numFixedDofs_;

    //! Array containing for equation numbers of last free dofs

    //! This array contains for every sub-system the equation number for
    //! the last free degree of freedom. In the StdMatrix case there is
    //! only a single system, so the length of the array is one. In the
    //! SBM_Matrix case the length corresponds to the number of diagonal
    //! matrix blocks.
    UInt *numFreeDofs_;

    //! Flag signalling whether built-in of Dirichlet BCs is possible

    //! This attribute is used as a flag for the internal state of the
    //! %IDBC_Handler object. If addIDBCPossible_ = true, then it is safe
    //! to call AddIDBCToRHS(). The flag monitors the status of the internal
    //! system matrixis. Thus, it is set to false, whenever one of the other
    //! FE matrices is changed and set to true by calling
    //! BuiltSystemMatrix().
    bool addIDBCPossible_;

    //! Flag signalling whether removal of Dirichlet BCs is possible

    //! This attribute is used as a flag for the internal state of the
    //! %IDBC_Handler object. If remIDBCPossible_ = true, then it is safe
    //! to call RemoveIDBCFromRHS(). The flag monitors the status of the
    //! internal system matrix and vecIDBC_, i.e. the vector containing the
    //! Dirichlet values. The flag is set to true, whenever AddIDBCToRHS()
    //! is called and re-set to false, whenever one of the FE matrices or
    //! an entry of vecIDBC_ is changed.
    bool remIDBCPossible_;

    //! Flag indicating whether we are faced with SBM_Matrices/Vectors
    bool sbmCase_;

  };
}

#ifndef EXPLICIT_TEMPLATE_INSTANTIATION
//#include "idbchandler.cc"
#endif

#endif
