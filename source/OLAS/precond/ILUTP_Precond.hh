#ifndef OLAS_ILUTP_PRECOND_HH
#define OLAS_ILUTP_PRECOND_HH



#include "OLAS/utils/math/CroutLU.hh"

#include "BasePrecond.hh"
#include "BNPrecond.hh"

namespace CoupledField {

  //! This class implements an \f$\mbox{ILU}(\tau,p)\f$ preconditioner

  //! This class implements an \f$\mbox{ILU}(\tau,p)\f$ preconditioner. This
  //! kind of preconditioner uses an incomplete factorisation of the problem
  //! matrix in the form
  //! \f[
  //! A = LU - R \enspace,
  //! \f]
  //! where L and U are the factors of the incomplete factorisation and R is
  //! a remainder matrix that is not explicitely computed. The factorisation
  //! is incomplete, since some of the fill-in generated during the
  //! factorisation is not allowed, but discarded. The \f$\mbox{ILU}(\tau,p)\f$
  //! preconditioner uses the Crout variant of Gaussian elimination and employs
  //! the following following <em>dropping strategy</em> for discarding
  //! fill-in:\n
  //! -# Compute a new row \f$u_k\f$ of U
  //! -# Determine threshold value \f$\tau_k := \tau \cdot \|u_k\|_p\f$
  //! -# Drop all entries from \f$u_k\f$ with magnitude smaller than
  //!    \f$\tau_k\f$
  //! -# From the remaining entries keep the diagonal entry plus the fill ones
  //!    largest in absolute magnitude
  //!
  //! The dropping strategy for the columns of the second factor L works in
  //! precisely the same fashion. In the current implementation the norm used
  //! for step 2 of the above scheme is \f$p=1\f$, i.e. we use the mean of the
  //! absolute magnitude of all entries.
  //!
  //! An ILUTP_Precond object reads the following parameters from its myParams_
  //! object:
  //! \n\n
  //! <center>
  //!   <table border="1" width="80%" cellpadding="10">
  //!     <tr>
  //!       <td colspan="4" align="center">
  //!         <b>Parameters for %ILUTP_Precond object</b>
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td align="center"><b>ID string in OLAS_Params</b></td>
  //!       <td align="center"><b>range</b></td>
  //!       <td align="center"><b>default value</b></td>
  //!       <td align="center"><b>description</b></td>
  //!     </tr>
  //!     <tr>
  //!       <td>ILUTP_tau</td>
  //!       <td align="center">(0,1]</td>
  //!       <td align="center">0.01</td>
  //!       <td>This parameter determines the threshold \f$\tau\f$ used for
  //!           computing the row / column threshold \f$\tau_k\f$ in step 2
  //!           of the dropping strategy above.</td>
  //!     </tr>
  //!     <tr>
  //!       <td>ILUTP_fillVal</td>
  //!       <td align="center">whole number</td>
  //!       <td align="center">-2</td>
  //!       <td>This integer value is used to determned the upper bound on the
  //!           number of entries in the factors L and U. Each row / column of
  //!           U / L will have at most maxFill_ + 1 entries, see step 4 above.
  //!           maxFill_ is determined from fillVal by the following formula
  //!           \f[
  //!           \mbox{maxFill\_} = \left\{ \begin{array}{ll}
  //!           \displaystyle\mbox{fillVal} &
  //!           \displaystyle\mbox{ for } \mbox{fillVal} \geq 0\\
  //!           \displaystyle\mbox{fillVal} \cdot \left(
  //!           \frac{\mbox{nnz}(A)}{n(A)} - 1\right) &
  //!           \displaystyle\mbox{ for } \mbox{fillVal} < 0
  //!           \end{array}\right.
  //!           \f]
  //!       Here \f$\mbox{nnz}(A)\f$ denotes the number of non-zero entries
  //!       in the system matrix and \f$n(A)\f$ its dimension.
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>CROUT_saveFacToFile</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">false</td>
  //!       <td>If the value of this boolean parameter is true, then
  //!           factorisation matrix \f$F=L+U-I\f$ computed in the setup phase
  //!           of the solver, will be exported to a file in MatrixMarket
  //!           format. This IO operation will take place immediately after
  //!           the end of the setup phase.</td>
  //!     </tr>
  //!     <tr>
  //!       <td>CROUT_facFileName</td>
  //!       <td align="center">string</td>
  //!       <td align="center">-</td>
  //!       <td>File name for exporting the factorisation matrix. This value
  //!           is only required, when CROUT_saveFacToFile is 'true'.
  //!           </td>
  //!     </tr>
  //!   </table>
  //! </center>
  //!
  //! \note
  //! - Pivoting is currently not supported during the factorisation.
  //! - The class only works together with CRS_Matrices with scalar entries.
  template <typename T>
  class ILUTP_Precond : public BNPrecond<ILUTP_Precond<T>, CRS_Matrix<T>, T >,
			public CroutLU<T> {

  public:

    using BNPrecond<ILUTP_Precond<T>, CRS_Matrix<T>, T >::Apply;
    using BNPrecond<ILUTP_Precond<T>, CRS_Matrix<T>, T >::Setup;

    //! Constructor
    ILUTP_Precond( const StdMatrix &stdMat, PtrParamNode precondNode,
		   PtrParamNode olasInfo );

    //! Destructor
    ~ILUTP_Precond();

    //! Query type of preconditioner

    //! When called, this method returns the type of the preconditioner
    //! object. In the case of an object of this class the return
    //! value is ILUTP.
    BasePrecond::PrecondType GetPrecondType() const {
      return BasePrecond::ILUTP;
    };

    //! Setup function inherited from class BasePrecond
    void Setup( CRS_Matrix<T>& sysmatrix );

    //! Applies the preconditioner by solving \f$ LU z = r\f$ for \f$z\f$

    //! This method applies the preconditioner. In the case of the
    //! %ILUTP_Precond this mean that the equation \f$ LU z = r\f$ with the
    //! incomplete factors \f$L\f$ and \f$U\f$ is solved for \f$z\f$ using
    //! a pair of forward/backward substitutions.
    //! \param sysMat problem matrix (ignored)
    //! \param res    residual vector for current iteration step (\f$r\f$)
    //! \param sol    output vector computed by the preconditioner (\f$z\f$)
    void Apply( const CRS_Matrix<T> &sysMat, const Vector<T> &res,
		Vector<T> &sol );

  private:

    //! Dropping strategy for fill-in in the incomplete decomposition

    //! The CroutLU class requires that a derived class implements this method,
    //! which is used to drop fill-in and compute an incomplete LU
    //! decomposition. In the case of the ILUTP_Precond class we implement
    //! the strategy described above.
    //! \param k        index of current row / column
    //! \param vecZ     dense vector storing current row of factor U
    //! \param vecZFill vector containing indices of non-zero entries in vecZ
    //! \param vecW     dense vector storing current column of factor L
    //! \param vecWFill vector containing indices of non-zero entries in vecW
    void DropEntries( UInt k, std::vector<T> &vecZ,
		      std::vector<UInt> &vecZFill,
		      std::vector<T> &vecW, std::vector<UInt> &vecWFill );

    //! Default constructor

    //! The default constructor is not allowed, since we need size information
    //! and pointers to communication objects for corrected initialisation.
    ILUTP_Precond() {
      EXCEPTION( "Default constructor of ILUTP_Precond should never be called!" );
    };

    //! Keep track on status of factorisation

    //! The Setup routine will set this boolean to true, once it has computed
    //! a factorisation. This gives us a confortable way of keeping track of
    //! the status of the solver and the internal memory allocation.
    bool amFactorised_;

    //! \f$\tau\f$ parameter for ILU decomposition

    //! This attribute stores the threshold parameter \f$\tau\f$ used for
    //! computing the row / column threshold \f$\tau_k\f$ in step 2 of the
    //! dropping strategy in the ILU decomposition, see algoritm above.
    //! The value is obtained from the parameter object myParams_ each time
    //! Setup() is called.
    Double tau_;

    //! maxFill parameter for ILU decomposition

    //! This attribute stores the value of the maxFill parameter used for
    //! determining an upper bound on the number of non-zero entries in the
    //! factors L and U of the ILU decomposition, see step 4 of the algoritm
    //! above.
    //! The value is obtained from the parameter object myParams_ each time
    //! Setup() is called.
    UInt maxFill_;

    //! Auxilliary index vector
    std::vector<UInt> indexVec_;

    //! This class provides the predicate for sorting

    //! This class provides the predicate for sorting the auxilliary index
    //! vector used in determining the maxFill_ entries of the k-th row /
    //! column of the factor U resp. L with largest modulus after dropping all
    //! entries below the threshold \f$\tau\f$.
    class FindMaxEntries {

    public:

      //! Constructor takes reference to dense vector
      FindMaxEntries( std::vector<T> &denseVec ) : dVec_(denseVec) {};

      //! The predicate method for sorting

      //! This is the predicate method for sorting. The method will return
      //! true for the test \f$i < j\f$, if
      //! \f$|\mbox{dVec}_i| > |\mbox{dVec}_j|\f$ or
      //! \f$|\mbox{dVec}_i| = |\mbox{dVec}_j|\f$ and \f$i < j\f$. Thus, in the
      //! case of a tie in the magnitude of the dense vector entries, we give
      //! preference to those entries that are closer to the diagonal
      bool operator()( const UInt &i, const UInt &j ) const {
	bool retVal = false;

	// A problem is here with the abs (g++ using int abs(int)!!!)
	if ( abs(dVec_[i]) > abs(dVec_[j]) ) {
	  retVal = true;
	}
	else if ( abs(dVec_[i]) == abs(dVec_[j]) ) {
	  if ( i < j ) {
	    retVal = true;
    	  }
	}
	return retVal;
      };

      //! Reference to dense vector

      //! This class stores a reference to the dense vector, whose entries
      //! are used for the sorting.
      std::vector<T> &dVec_;
    };

  private:

    //! Output contents of an STL container to screen (for debugging purposes)
    template<typename conT>
    void PrintContainer( const conT &con, const char *conName ) const {
#ifdef DEBUG_ILUPRECOND
      (*debug) << " " << conName << ": ";
      typename conT::iterator it;
      for ( it = con.begin(); it != con.end(); it++ ) {
        (*debug) << *it << " ";
      }
      (*debug) << std::endl;
#else
      std::cerr << "DEBUG_ILUPRECOND must be defined to enjoy "
                   "method ILUTP_Precond::PrintContainer\n";
#endif
    }

  };

} // namespace CoupledField

#endif
