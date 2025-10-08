#ifndef IDBC_HANDLER_PENALTY_HH
#define IDBC_HANDLER_PENALTY_HH

#include <set>


#include "BaseIDBC_Handler.hh"

namespace CoupledField {


  // forward declarations
  class BaseGraphManager;


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

     //! @copydoc BaseIDBC_Handler::AdaptSystemMatrix()
    void AdaptSystemMatrix( SBM_Matrix &sysMat );

    //! @copydoc BaseIDBC_Handler::AddIDBCToRHS()
    void AddIDBCToRHS( SBM_Vector *rhs, bool deltaIDBC = false );

    //! @copydoc BaseIDBC_Handler::SetIDBC()
    void SetIDBC( UInt rowBlock, UInt rowNum, const T &val );
    
    //! @copydoc BaseIDBC_Handler::GetIDBC()
    void GetIDBC( UInt rowBlock, UInt rowNum, T &val, bool deltaIDBC=false );

    //! @copydoc BaseIDBC_Handler::SetDofsToIDBC()
    void SetDofsToIDBC( SBM_Vector *vec, bool deltaIDBC = false );

    // =======================================================================
    // CONSTRUCTION, DESTRUCTION and CLEARING
    // =======================================================================

    //! \name Constructors, destructor and initialisation methods

    //@{

    //! Admissible constructor for the case of SBM_Matrices

    //! This constructor is intended to by used in the case of a linear
    //! system specified with the help of SBM_Matrices.
    //! \param numIDBC   vector containing for each SBM row the number of
    //!                  inhomogeneous Dirichlet values (length of vector
    //!                  = number of sbmRows)
    IDBC_HandlerPenalty( StdVector<UInt> numIDBC);

    //! Destructor
    ~IDBC_HandlerPenalty();

    //! @copydoc BaseIDBC_Handler::InitDirichletValues()
    void InitDirichletValues();

    void SetOldDirichletValues();

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

    ///! @copydoc BaseIDBC_Handler::BuiltSystemMatrix()
    void BuildSystemMatrix( const std::map<FEMatrixType, Double> &factors,
                            std::map<UInt, std::set<UInt> >& indicesPerBlock ) {
    }
    void BuildSystemMatrix( const std::map<FEMatrixType, Complex> &factors, std::map<UInt, std::set<UInt> >& indicesPerBlock ) {};

    //! @copydoc BaseIDBC_Handler::RemoveIDBCFromRHS()
    inline void RemoveIDBCFromRHS( SBM_Vector *rhs, bool deltaIDBC = false ) {
    }

    //! @copydoc BaseIDBC_Handler::AddWeightFixedToFree()
    void AddWeightFixedToFree( FEMatrixType matID,
                               UInt rowBlock,
                               UInt colBlock,
                               UInt rowInd,
                               UInt colInd,
                               const T& val ) {
    };

    //! @copydoc BaseIDBC_Handler::SetWeightFixedToFree()
    void SetWeightFixedToFree( FEMatrixType matID,
                               UInt rowBlock,
                               UInt colBlock,
                               UInt rowInd,
                               UInt colInd,
                               const T& val ) {
    };

    //! @copydoc BaseIDBC_Handler::GetWeightFixedToFree()
    void GetWeightFixedToFree( FEMatrixType matID,
                               UInt rowBlock,
                               UInt colBlock,
                               UInt rowInd,
                               UInt colInd,
                               T & val ) {
    };
    
    //! @copydoc BaseIDBC_Handler::InitMatrix()
    void InitMatrix( FEMatrixType matrixType) {
    };

    //! Return penalty term
    Double GetPenalty() const {
      return penaltyTerm_;
    }
    
    //! Dump information to std output
    void Dump()
    {
      std::map<UInt, std::map<Integer, UInt> >::iterator i;
      for(i = bcIndices_.begin(); i != bcIndices_.end(); i++)
      {
        std::cout << " sbmBlock = " << i->first << std::endl;
        std::map<Integer, UInt>::iterator j;
        for(j = i->second.begin(); j != i->second.end(); j++)
          std::cout << " " << dirichletEqns_[i->first][j->second];
        std::cout << std::endl;
      }
    }
    //@}

    void PrintIDBCvec(){
      EXCEPTION("Not implemented here");
    }


  private:
    //! Keep for every SBM-block the eqn numbers
    std::map<UInt, StdVector<UInt> > dirichletEqns_;
    
    //! Map for mapping eqn number to position in dirichletValue_ vector
    std::map<UInt, std::map<Integer, UInt> > bcIndices_;

    //! Array for storing inhomogeneous Dirichlet values
    SBM_Vector *dirichletValue_;
    SBM_Vector *oldDirichletValue_;

    //! Value of the penalty term in the penalty approach
    Double penaltyTerm_;

    //! Attribute storing run-time information on class template
    BaseMatrix::EntryType eType_;

    //! Method for formatted output of internal arrays

    //! This method will print the information on the inhomogeneous
    //! Dirichlet boundary conditions (equation numbers, components,
    //! values) stored in the internal arrays to the specified stream.
    //! It is primarily intended for debugging purposes.
    void PrintInternalState( std::ostream &os );

  };

}

#endif
