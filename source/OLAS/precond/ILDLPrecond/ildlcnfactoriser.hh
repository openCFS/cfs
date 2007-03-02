#ifndef ILDLCN_FACTORISER_HH
#define ILDLCN_FACTORISER_HH

#include "utils/utils.hh"
#include "matvec/matvec.hh"
#include "precond/ILDLPrecond/baseildlfactoriser.hh"


namespace OLAS {

  //! This class implements a inverse-based incomplete LDL factorisation

  //! This class implements an incomplete LDL factorisation. In this special
  //! variant the dropping of entries is based on an estimation of their
  //! influence on the error in the inverse of the factor L. An entry
  //! \f$l_{i,k}\f$ in the factor \f$L\f$ of the incomplete factorisation
  //! \f$A=LDL^T - R\f$ will be discarded, if
  //! \f[
  //! |l_{i,k}| \cdot \| e_k^T L_k^{-1} \|_\infty \leq \tau
  //! \enspace.
  //! \f]
  //! Here \f$e_k\f$ denotes the k-th canonical base vector, i.e.
  //! \f$e_k^T L_k^{-1}\f$ represents the k-th row in the inverse of \f$L\f$.
  //! The inverse of \f$L\f$ is of course not directly available an algorithm
  //! borrowed from the field of condition number estimation is employed to
  //! approximate \f$\| e_k^T L_k^{-1} \|_\infty\f$.
  //!
  //! An ILUCN_Precond object reads the following parameters from its myParams_
  //! object:
  //! \n\n
  //! <center>
  //!   <table border="1" width="80%" cellpadding="10">
  //!     <tr>
  //!       <td colspan="4" align="center">
  //!         <b>Parameters for %ILUCN_Precond object</b>
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
  //!   </table>
  //! </center>
  template <class T>
  class ILDLCNFactoriser : public BaseILDLFactoriser<T> {

  public:

    //! Default constructor
    ILDLCNFactoriser();

    //! Standard constructor
    ILDLCNFactoriser( OLAS_Params *myParams, OLAS_Report *myReport = NULL );

    //! Default destructor
    ~ILDLCNFactoriser();

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

    //! Drop entries from factor row depending on influence on \f$L^{-1}\f$

    //! This method turns the full LDL factorisation into an incomplete one by
    //! dropping entries from the factor rows.
    void DropEntries( UInt *listIDX, T *listVAL, Double tau, Double norm );

  };

}

#endif
