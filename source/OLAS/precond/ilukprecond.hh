// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef ILUK_PRECOND_HH
#define ILUK_PRECOND_HH

#include <vector>

#include <def_expl_templ_inst.hh>

#include "OLAS/utils/math/croutlu.hh"

#include "baseprecond.hh"
#include "bnprecond.hh"

namespace CoupledField {

  //! This class implements an ILU(k) preconditioner

  //! This class implements an ILU(k) preconditioner. This kind of
  //! preconditioner uses an incomplete factorisation of the problem
  //! matrix in the form
  //! \f[
  //! A = LU - R \enspace,
  //! \f]
  //! where L and U are the factors of the incomplete factorisation and R is
  //! a remainder matrix that is not explicitely computed. The factorisation
  //! is incomplete, since some of the fill-in generated during the
  //! factorisation is not allowed, but discarded. This ILU(k) preconditioner
  //! uses the Crout variant of Gaussian elimination and employs the following
  //! <em>dropping strategy</em> for discarding fill-in. Each matrix entry
  //! initially receives a fill-value according to the following rule:
  //! \f[
  //! \mbox{lev}_{ij} = \left\{ \begin{array}{ll}
  //!                   0 & \mbox{ if } a_{ij} \neq 0 \mbox{ or } i = j \\
  //!                   \infty & \mbox{ otherwise}
  //!                   \end{array}\right.
  //! \f]
  //! The fill-value of an entry is re-computed when it is touched during
  //! the factorisation based on the fill-value of the parent entries
  //! responsible for the operation according to the following formula
  //! \f[
  //! \mbox{lev}_{ij} = \min\{ \mbox{lev}_{ij}, \mbox{lev}_{ik} +
  //! \mbox{lev}_{kj} + 1 \} \enspace.
  //! \f]
  //! Thus, the fill-value of an entry cannot grow and the original entries
  //! in the matrix pattern will always remain at \f$\mbox{lev}_{ij} = 0 \f$.
  //! Fill-in is discarded, if its final fill-value, i.e. the fill-value once
  //! the entry has completely been computed, is larger than the threshold
  //! value k.
  //!
  //! An %ILUK_Precond object reads the following parameters from its myParams_
  //! object:
  //! \n\n
  //! <center>
  //!   <table border="1" width="80%" cellpadding="10">
  //!     <tr>
  //!       <td colspan="4" align="center">
  //!         <b>Parameters for %ILUK_Precond object</b>
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td align="center"><b>ID string in OLAS_Params</b></td>
  //!       <td align="center"><b>range</b></td>
  //!       <td align="center"><b>default value</b></td>
  //!       <td align="center"><b>description</b></td>
  //!     </tr>
  //!     <tr>
  //!       <td>ILUK_level</td>
  //!       <td align="center">non-negative integer</td>
  //!       <td align="center">1</td>
  //!       <td>This integer value determines the threshold for the level
  //!           \f$\mbox{lev}_{ij}\f$ up to which fill-in is allowed. If
  //!           \f$\mbox{lev}_{ij}>\mbox{ILUK\_level}\f$ then the fill-in
  //!           at position (i,j) is discarded.
  //!       </td>
  //!     </tr>
  //!     <tr>
  //!       <td>ILUK_logging</td>
  //!       <td align="center">true or false</td>
  //!       <td align="center">true</td>
  //!       <td>If the parameter is true, then the object will log some
  //!           information on the factorisation to the standard %OLAS
  //!           logstream.
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
  class ILUK_Precond : public BNPrecond<ILUK_Precond<T>, CRS_Matrix<T>, T >,
                       public CroutLU<T> {

  public:

    using BNPrecond<ILUK_Precond<T>, CRS_Matrix<T>, T >::Apply;
    using BNPrecond<ILUK_Precond<T>, CRS_Matrix<T>, T >::Setup;

    //! Constructor
    ILUK_Precond( const StdMatrix &stdMat, PtrParamNode solverNode,
                   PtrParamNode olasInfo );

    //! Destructor
    ~ILUK_Precond();

    //! Query type of preconditioner

    //! When called, this method returns the type of the preconditioner
    //! object. In the case of an object of this class the return
    //! value is ILUK.
    BasePrecond::PrecondType GetPrecondType() const {
      return BasePrecond::ILUK;
    };

    //! Setup function inherited from class BasePrecond
    void Setup( CRS_Matrix<T>& sysmatrix );

    //! Applies the preconditioner by solving \f$ LU z = r\f$ for \f$z\f$

    //! This method applies the preconditioner. In the case of the
    //! %ILUK_Precond this mean that the equation \f$ LU z = r\f$ with the
    //! incomplete factors \f$L\f$ and \f$U\f$ is solved for \f$z\f$ using
    //! a pair of forward/backward substitutions.
    //! \param sysMat problem matrix (ignored)
    //! \param res    residual vector for current iteration step (\f$r\f$)
    //! \param sol    output vector computed by the preconditioner (\f$z\f$)
    void Apply( const CRS_Matrix<T> &sysMat, const Vector<T> &res,
                Vector<T> &sol ) const;

  private:

    //! Dropping strategy for fill-in in the incomplete decomposition

    //! The CroutLU class requires that a derived class implements this method,
    //! which is used to drop fill-in and compute an incomplete LU
    //! decomposition. In the case of the ILUK_Precond class we implement
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
    ILUK_Precond() {
      EXCEPTION( "Default constructor of ILUK_Precond should never be called!" );
    };

    //! Dimension of system matrix
    UInt problemDim_;

    //! Class providing predicate for discarding fill-in in ILUK_Precond

    //! This class provides the predicate for discarding all fill-entries
    //! whose fill-level is larger than the user specified threshold in the
    //! ILUK_Precond class.
    class TestFillLevel {

    public:

      //! Constructor takes threshold, fill-levels and values of vector entries
      TestFillLevel( UInt threshold, std::vector<UInt> &levels,
		     std::vector<T> &entries, UInt maxVal ) {
	thresh_ = threshold;
	level_  = levels;
	entry_  = entries;
	maxVal_ = maxVal;
      }

      //! The predicate method for testing fill-in

      //! This is the predicate method for testing the fill-in. It receives
      //! an entry of the index set of non-zero vector entries and compares
      //! the fill-level of the corresponding entry to the threshold. If the
      //! test succeeds, the method returns true, zeros the vector entry and
      //! re-sets its fill-level to a maximal unattainable value.
      //! \return
      //! - true if level > maxLevel_
      //! - false if level <= maxLevel_
      bool operator()( const UInt &index ) {
	bool retVal = false;

#ifdef DEBUG_ILUKPRECOND
	(*debug) << "Entry " << index << ": lvl = " << level_[index];
#endif

	// Level larger than threshold
	if ( level_[index] > thresh_ ) {

#ifdef DEBUG_ILUKPRECOND
	  (*debug) << "> " << thresh_ << "------- DROPPED" << std::endl;
#endif

	  // test failed
	  retVal = true;

	  // zero vector entry
	  entry_[index] = 0.0;

	  // set level to maximal value
	  level_[index] = maxVal_;
	}
	else {

#ifdef DEBUG_ILUKPRECOND
	  (*debug) << std::endl;
#endif

	}
	return retVal;

      };

      //! Maximal fill-level
      UInt thresh_;

      //! Fill-levels of vector entries
      std::vector<UInt> level_;

      //! Values of vector entries
      std::vector<T> entry_;

      //! Penalty value
      UInt maxVal_;

    };

    //! Maximum allowed fill-level

    //! This attribute stores the threshold for the maximum allowed fill-level
    //! supplied by the user. It is (re-)set each time Setup() is called and
    //! is only required to be passed to the constructor of the TestFillLevel
    //! predicate class
    UInt maxLevel_;

    //! Flag, if logging is performed
    bool logging_;
  };

}

#ifndef EXPLICIT_TEMPLATE_INSTANTIATION
//#include "ilukprecond.cc"
#endif

#endif
