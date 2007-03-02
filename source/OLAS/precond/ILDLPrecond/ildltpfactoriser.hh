#ifndef ILDLTP_FACTORISER_HH
#define ILDLTP_FACTORISER_HH

#include "utils/utils.hh"
#include "matvec/matvec.hh"
#include "precond/ILDLPrecond/baseildlfactoriser.hh"


namespace OLAS {

  //! This class implements an incomplete LDL(tau,p) factorisation

  //! This class implements an incomplete LDL factorisation. In this special
  //! variant the dropping of entries is based on their relative size,
  //! ILDL(tau,p) and so forth
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
  //!       <td>ILDLPRECOND_tau</td>
  //!       <td align="center">(0,1]</td>
  //!       <td align="center">0.01</td>
  //!       <td>This parameter determines the threshold \f$\tau\f$ used for
  //!           computing the row / column threshold \f$\tau_k\f$ in step 2
  //!           of the dropping strategy above.</td>
  //!     </tr>
  //!     <tr>
  //!       <td>ILDLPRECOND_fillVal</td>
  //!       <td align="center">whole number</td>
  //!       <td align="center">-2</td>
  //!       <td>This integer value is used to determine the upper bound on the
  //!           number of entries in the factor L. Each row of L will have at
  //!           most maxFill + 1 entries, see step 4 above.
  //!           maxFill is determined from fillVal by the following formula
  //!           \f[
  //!           \mbox{maxFill\_} = \left\{ \begin{array}{ll}
  //!           \displaystyle\mbox{fillVal} &
  //!           \displaystyle\mbox{ for } \mbox{fillVal} \geq 0\\
  //!           \displaystyle\mbox{fillVal} \cdot \left(
  //!           \frac{\mbox{nnz}(A) - n(A)}{2\,n(A)}\right) &
  //!           \displaystyle\mbox{ for } \mbox{fillVal} < 0
  //!           \end{array}\right.
  //!           \f]
  //!       Here \f$\mbox{nnz}(A)\f$ denotes the number of non-zero entries
  //!       in the system matrix and \f$n(A)\f$ its dimension.
  //!       </td>
  //!     </tr>
  //!   </table>
  //! </center>
  template <class T>
  class ILDLTPFactoriser : public BaseILDLFactoriser<T> {

  public:

    //! Default constructor
    ILDLTPFactoriser();

    //! Standard constructor
    ILDLTPFactoriser( OLAS_Params *myParams, OLAS_Report *myReport = NULL );

    //! Default destructor
    ~ILDLTPFactoriser();

    //! Incomplete \f$LDL^T\f$ factorisation of a square matrix

    //! This is the central method of the class. It computes an
    //! incomplete \f$LDL^T\f$ factorisation of a square matrix using the
    //! ILDL(tau,p) variant as special dropping strategy.
    //! \param sysMat     the matrix to be factorised
    //! \param dataD
    //! \param rptrU
    //! \param cidxU
    //! \param dataU
    //! \param newPattern for the ILDL(tau,p) factoriser it does not matter,
    //!                   whether the matrix pattern has changed or not; a new
    //!                   numerical factorisation must be computed in any case
    //!                   thus, the parameter is ignored
    void Factorise( SCRS_Matrix<T> &sysMat, std::vector<T> &dataD,
                    std::vector<UInt> &rptrU, std::vector<UInt> &cidxU,
                    std::vector<T> &dataU, bool newPattern );

    //! This class provides the predicate for sorting

    //! This class provides the predicate for sorting the auxilliary index
    //! vector used in determining the maxFill_ entries of the k-th row /
    //! column of the factor U resp. L with largest modulus after dropping all
    //! entries below the threshold \f$\tau\f$.
    class FindMaxEntries {

    public:

      //! Constructor takes reference to dense vector
      FindMaxEntries( T *numValue ) : numValue_(numValue) {};

      //! The predicate method for sorting

      //! This is the predicate method for sorting. The comparison is based on
      //! the magnitude of the row entries associated with the entries (column
      //! indices) in the STL vector that is sorted. In the case of a tie it
      //! gives preference to those entries that are closer to the diagonal
      //! \return
      //! \f[
      //! (i < j) = \mbox{true, if} \left\{ \begin{array}{lll}
      //! v_i > v_j & &\\[2ex]
      //! v_i = v_j & \mbox{and} & i < j
      //! \end{array}\right.
      //! \f]
      bool operator()( const UInt &i, const UInt &j ) const {
	bool retVal = false;

	// A problem is here with the abs (g++ using int abs(int)!!!)
	if ( abs(this->dVec_[i]) > abs(this->dVec_[j]) ) {
	  retVal = true;
	}
	else if ( abs(this->dVec_[i]) == abs(this->dVec_[j]) ) {
	  if ( i < j ) {
	    retVal = true;
    	  }
	}
	return retVal;
      };

      //! Reference to array with numerical values of row entries

      //! This class stores a reference to a vector, in which the numerical
      //! values for all row entries are stored. This is the value array of
      //! for the linked list, resp. to the STL vector containing their
      //! column indices
      T* numValue_;
    };

  private:

    //! List with indices of rows to be scanned

    //! This list contains the indices of those rows that must be scanned
    //! when updating the auxilliary index vectors after each elimination step.
    //! The entries in the list are in ascending row index order.
    // UInt *scanList_;

    //! List with indices of rows contributing in current elimination step

    //! List with indices of rows contributing in current elimination step.
    //! The entries in the list are in ascending row index order and
    //! activeList_ is a subset of scanList_.
    // UInt *activeList_;

    //! Drop entries from factor row depending on their modulus

    //! This method turns the full LDL factorisation into an incomplete one by
    //! dropping entries from the factor rows.
    void DropEntries( UInt *listIDX, T *listVAL, Double tau, UInt maxFill );

  };

}

#endif
