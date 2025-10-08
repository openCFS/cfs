#ifndef IDBC_HANDLER_HH
#define IDBC_HANDLER_HH

#include <set>

#include "BaseIDBC_Handler.hh"


#include "MatVec/SBM_Matrix.hh"



namespace CoupledField {


  // forward declarations
  class GraphManager;

  // ========================================================================
  // IDBC_HANDLER CLASS
  // ========================================================================


  //! The IDBCHandler class is responsible for storing the weights with
  //! which real degrees of freedom couple to those that are fixed by
  //! <b>I</b>nhomogeneous <b>D</b>irichlet <b>B</b>oundary <b>C</b>onditions
  //! and for updating the right hand side of the linear system using these
  //! weights and the prescribed Dirichlet values.
  //! \note The openCFS interface arranges equation numbers for degrees of
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
    //! \param numBlocks      number of sbmBlocks for which we administrate
    //!                       inhomogeneous Dirichlet boundary conditions
    IDBC_Handler( const std::set<FEMatrixType> &usedFEMatrices,
                  GraphManager *graphManager, UInt numBlocks);

    //MH-Constructor
    IDBC_Handler( const std::set<FEMatrixType> &usedFEMatrices,
                  GraphManager *graphManager, UInt numBlocks, UInt sbmRow, UInt sbmCol);
    //! Destructor
    ~IDBC_Handler();

    //! @copydoc BaseIDBC_Handler::InitMatrix()
    void InitMatrix( FEMatrixType matrixType ) {


      // Check for validity of matrix identifier
      mySetIterator it = myMatrices_.find( matrixType );
      if ( it == myMatrices_.end() ) {
        EXCEPTION( "IDBC_Handler::InitMatrix: Received invalid matrix "
                 << "identifier '" << matrixType << "'");
      }
      else {
        auxMat_[ *it ]->Init();
      }

      // Adapt internal status flags
      addIDBCPossible_ = false;
      remIDBCPossible_ = false;
    }

    //! @copydoc BaseIDBC_Handler::InitDirichletValues()
    void InitDirichletValues() {
      vecIDBC_->Init();
      vecOldIDBC_->Init();
    };

    void SetOldDirichletValues();
    void ToString();

    //@}

    // =======================================================================

    //! @copydoc BaseIDBC_Handler::BuiltSystemMatrix()
    void BuildSystemMatrix( const std::map<FEMatrixType, Double> &factors,
                            std::map<UInt, std::set<UInt> >& indicesPerBlock );
    //! OVERLOAD: BuildSystemMatrix with complex factors
    void BuildSystemMatrix( const std::map<FEMatrixType, Complex> &factors,
                            std::map<UInt, std::set<UInt> >& indicesPerBlock );
    //! @copydoc BaseIDBC_Handler::AddIDBCToRHS()
    void AddIDBCToRHS( SBM_Vector *rhs, bool deltaIDBC = false  );

    //! @copydoc BaseIDBC_Handler::RemoveIDBCFromRHS
    void RemoveIDBCFromRHS( SBM_Vector *rhs, bool deltaIDBC = false );

    //! @copydoc BaseIDBC_Handler::SetIDBC()
    void SetIDBC( UInt blockNum, UInt index, const T &val );
    
    //! @copydoc BaseIDBC_Handler::GetIDBC()
    void GetIDBC( UInt blockNum, UInt index, T &val, bool deltaIDBC=false );

    // Prints the IDBC-Vector
    void PrintIDBCvec();

    //! @copydoc BaseIDBC_Handler::AddWeightFixedToFree()
    void AddWeightFixedToFree( FEMatrixType matID,
                               UInt rowBlock,
                               UInt colBlock,
                               UInt rowInd,
                               UInt colInd,
                               const T& val );
    
    //! @copydoc BaseIDBC_Handler::SetWeightFixedToFree()
    void SetWeightFixedToFree( FEMatrixType matID,
                               UInt rowBlock,
                               UInt colBlock,
                               UInt rowInd,
                               UInt colInd,
                               const T& val );

    //! @copydoc BaseIDBC_Handler::GetWeightFixedToFree()
    void GetWeightFixedToFree( FEMatrixType matID,
                               UInt rowBlock,
                               UInt colBlock,
                               UInt rowInd,
                               UInt colInd,
                               T& val );
    
    //! @copydoc BaseIDBC_Handler::AddFixedToFreeRHS()
    void AddFixedToFreeRHS( FEMatrixType matID, UInt colBlock,
                                        UInt colInd, SBM_Vector *rhs, const T& val );

    //! @copydoc BaseIDBC_Handler::SetDofsToIDBC()
    void SetDofsToIDBC( SBM_Vector *vec, bool deltaIDBC = false );


    // =======================================================================
    // STUBS
    // =======================================================================

    //! \name Stubs
    //! The following method are not required for this derivation of the
    //! BaseIDBC_Handler class. For simplicity we decided to make the methods
    //! stubs, i.e. if they are called, they will simply return without
    //! performing any real action.

    //@{

    //! @copydoc BaseIDBC_Handler::AdaptSystemMatrix
    void AdaptSystemMatrix( SBM_Matrix &sysMat ) {
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

    //! Array containing links to all internal auxiliary matrices
    std::map<FEMatrixType, SBM_Matrix *> auxMat_;
    
    //! Array containing links to all internal auxiliary vectors
    
    //! These vectors are used for the multiplication of the auxiliary
    //! matrix with a given IDBC vector
    SBM_Vector  *auxVec_;

    //! Vector for storing the inhomogeneous Dirichlet values
    SBM_Vector *vecIDBC_;

    //! Vector for storing the inhomogeneous Dirichlet values from a previous
    //! time step or iteration -> needed to compute delta IDBC values, i.e.
    //! vecIDBC-vecOldIDBC
    SBM_Vector *vecOldIDBC_;

    // Entry type of auxiliary matrices

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
    StdVector<UInt> numFixedDofs_;

    //! Array containing for equation numbers of last free dofs

    //! This array contains for every sub-system the equation number for
    //! the last free degree of freedom. In the StdMatrix case there is
    //! only a single system, so the length of the array is one. In the
    //! SBM_Matrix case the length corresponds to the number of diagonal
    //! matrix blocks.
    StdVector<UInt> numFreeDofs_;

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
  };
}

#endif
