#ifndef OLAS_ILU0_PRECOND_HH
#define OLAS_ILU0_PRECOND_HH

#include <def_expl_templ_inst.hh>

#include "BasePrecond.hh"
#include "BNPrecond.hh"

namespace CoupledField {

  template<typename> class CRS_Matrix;

  //! ILU(0) Preconditioner
  
  //! This class implements a simple ILU(0) preconditioner. In an ILU
  //! preconditioner an approximate factorisation of the system matrix \f$A\f$
  //! is computed as
  //! \f[
  //! A = L\cdot U + R
  //! \f]
  //! where \f$L\f$ is a lower triangular, \f$U\f$ an upper triangular and
  //! \f$R\f$ a remainder matrix. The preconditioner is then given by
  //! \f$M^{-1}= (LU)^{-1}\f$.
  //! The notion ILU(0) indicates that no fill-in
  //! is allowed in the factorisation. Thus the matrices \f$A\f$ and \f$L+U\f$
  //! will have the same sparsity pattern.
  //! Application of the ILU preconditioner means computation of
  //! \f$M^{-1}r=z\f$. This is achieved by solving \f$Mz=r\f$ for \f$z\f$.
  //! The latter can simply be performed by a pair of backward/forward
  //! substitutions using the approximate factors \f$L\f$ and \f$U\f$.
  //! \note Currently this class only supports preconditioning for matrices in
  //! the CRS format with scalar complex of real entries.
  template <typename T>
  class ILU0Precond : public BNPrecond<ILU0Precond<T>,CRS_Matrix<T>,T> {

  public:

    using BNPrecond<ILU0Precond<T>,CRS_Matrix<T>,T>::Apply;
    using BNPrecond<ILU0Precond<T>,CRS_Matrix<T>,T>::Setup;

    //!
    typedef typename AssocType<T>::T_Mtype T_Mtype;
    //!
    typedef typename AssocType<T>::T_Vtype T_Vtype;
    //!
    typedef typename AssocType<T>::T_Stype T_Stype;

    //! Constructor (for use in GenerateStdPrecondObject)
    ILU0Precond( const StdMatrix &mat, PtrParamNode solverNode,
		 PtrParamNode olasInfo );

    //! Deep Destructor
    ~ILU0Precond();

    //! Apply ILU preconditioner
    
    //! Application of the ILU preconditioner means computation of
    //! \f$M^{-1}r=z\f$. This is achieved by solving \f$Mz=r\f$ for \f$z\f$.
    //! The latter can simply be done by performing a pair of backward/forward
    //! substitution using the approximate factors \f$L\f$ and \f$U\f$.
    void Apply( const CRS_Matrix<T> &sysmat, const Vector<T> &r,
		Vector<T> &z );

    //! Triggers setup of the ILU Preconditioner

    //! In the setup phase the matrix is factored into the product of a lower
    //! triangular matrix L and an upper triangular matrix U. In order to keep
    //! the factors sparse, the factorisation is only incomplete, i.e.
    //! A = LU + R with a remainder matrix R. LU preserves the sparsity pattern
    //! of A, so this is an ILU(0) factorisation.
    void Setup( CRS_Matrix<T> & sysmat, PtrParamNode analysis_id );

    //! Query type of preconditioner object

    //! When called this method returns the type of the preconditioner object.
    //! In the case of an object of this class the return value is ILU0.
    BasePrecond::PrecondType GetPrecondType() const {
      return BasePrecond::ILU0;
    };
    
    //! Export ILU factorisation to a file in MatrixMarket format

    //! The method will export the factorisation matrix to an ascii file
    //! according to the MatrixMarket specifications. By factorisation
    //! matrix we understand the matrix \f$F=L+U-I\f$.
    //! For details of the specification see http://math.nist.gov/MatrixMarket
    //! \param fname name of output file
    void ExportILUFactorisation( const std::string& fileName );


  private:
 
    /*! \name CRS format for the incomplete factors L and U
     */
    //@{
	
    //! dimension of system
    UInt dim_;
    
    //! nonzero entries
    T_Mtype *ilu_data_;
    
    //! row pointer
    UInt *ilu_rptr_;
    
    //! column indices
    UInt *ilu_cidx_;
    
    //!positions of the diagonal elements in ILU structure
    UInt *diagPos_;

    //! stores the number of unknowns
    UInt size_;
    
    //! flag for logging information
    bool logging_;
    
    //@}

  };

}//namespace

#ifndef EXPLICIT_TEMPLATE_INSTANTIATION
//#include "ilu0precond.cc"
#endif

#endif

