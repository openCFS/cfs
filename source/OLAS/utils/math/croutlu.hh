// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef CROUT_LU_HH
#define CROUT_LU_HH

#include <vector>

#include <def_expl_templ_inst.hh>

#include "General/defs.hh"

namespace CoupledField {

  template<typename> class CRS_Matrix;
  template<typename> class Vector;
  

  //! (I)LU factorisation of a matrix by Crout variant of Gaussian Elimination

  //! This class implements an LU decomposition \f$A = LU\f$ of a matrix and
  //! the associated solution of a linear system \f$LUx=y\f$ by a pair of
  //! forward-backward substitutions. The LU decomposition can either be done
  //! in a complete fashion for use as sparse direct solver, or in an
  //! incomplete fashion, so that it can be employed as preconditioner.
  //! This implementation of the (I)LU decomposition uses the Crout variant of
  //! Gaussian elimination. For details on this approach and the connection
  //! to ILU see e.g.\n
  //! <b>Youssef Saad, <em>Iterative Methods for Sparse Linear Systems</em>,
  //! SIAM</b>.\n\n
  //! \note
  //! - This class is <b>only</b> a base class and no objects of this class
  //!   can be instantiated, because the method DropEntries() is purely
  //!   virtual.
  //! - While the interface currently requires to supply a matrix in
  //!   compressed row storage (CRS) format, the implementation of the Crout
  //!   method does not depend on this. Knowledge of the special data layout
  //!   of the system matrix is completely encapsulated in the two private
  //!   methods GetRowOfA() and GetColOfA().
  //! - The implementation does currently not make use of pivoting for
  //!   numerical stability.
  template <typename T>
  class CroutLU {

  public:

    // Default constructor
    CroutLU();

    // Default denstructor
    virtual ~CroutLU();

    //! Export ILU factorisation to a file in MatrixMarket format

    //! The method will export the factorisation matrix to an ascii file
    //! according to the MatrixMarket specifications. By factorisation
    //! matrix we understand the matrix \f$F=L+U-I\f$.
    //! For details of the specification see http://math.nist.gov/MatrixMarket
    //! \param fname name of output file
    void ExportILUFactorisation( const char *fname );

  protected:

    //! Dimension of the square system matrix
    UInt sysMatDim_;

    //! Number of non-zero entries of the system matrix
    UInt sysMatNNZ_;

    //! Compute an (I)LU factorisation of the input matrix
    virtual void Factorise( CRS_Matrix<T> &sysMat );

    //! Solve linear system with LU matrix
    void Solve( const Vector<T> &rhs, Vector<T> &sol ) const;

    //! Clear the object

    //! This method can be used to reset all internal data structures to
    //! their initial state and delete the factors L and U computed by a
    //! previous factorisation, thus preparing the object to compute a
    //! new (I)LU decomposition.
    void Reset();

    //! Dropping strategy for fill-in in an incomplete decomposition

    //! This method is called each time a new row / column of the
    //! decomposition was computed by the Crout scheme. A derived class can
    //! implement this method according to its liking. One can do nothing,
    //! as would be appropriate for a sparse direct solver, or apply some
    //! dropping strategy in order to limit fill-in and obtain an incomplete
    //! decomposition to be used as preconditioner.
    //! \param k        index of current row / column
    //! \param vecZ     dense vector storing current row of factor U
    //! \param vecZFill vector containing indices of non-zero entries in vecZ
    //! \param vecW     dense vector storing current column of factor L
    //! \param vecWFill vector containing indices of non-zero entries in vecW
    virtual void DropEntries( UInt k, std::vector<T> &vecZ,
			      std::vector<UInt> &vecZFill,
			      std::vector<T> &vecW,
			      std::vector<UInt> &vecWFill ) = 0;

    //! An estimate for the growth of the number of non-zero entries

    //! This parameter is used to predict the growth of the number of
    //! non-zero entries during the elimination. This prediction is used for
    //! the initial reservation of memory in dynamic data structures. The
    //! number of non-zero entries in the factor U is e.g. estimated as
    //! \code (sysMatNNZ_ + sysMatDim_ ) / 2 * memGrowthEstimate_ \endcode
    UInt memGrowthEstimate_;


    // =======================================================================
    // AUXILLIARY VECTORS FOR ILUK DECOMPOSITION
    // =======================================================================

    //@{
    //! \name Auxilliary vectors and attributes for ILU(k) decomposition
    //! These vectors are used to store the fill-levels of the entries that
    //! are computed for the current row of U and column of L in the Crout
    //! factorisation and of the values in the final factors. Memory is only
    //! allocated and computations performed, if compFillLevels_ is true.

    //! Array containing fill-levels of current row vector of U
    std::vector<UInt> vecWLevel_;

    //! Array containing fill-levels of current column vector of L
    std::vector<UInt> vecZLevel_;

    //! Array containing fill-levels of non-zero entries of U
    std::vector<UInt> levelU_;

    //! Array containing fill-levels of non-zero entries of L
    std::vector<UInt> levelL_;

    //! Control setting of fill-pattern in GetColOfA() and GetRowOfA()
    bool compFillLevels_;
    //@}

  private:

    //! This attribute is used to keep track of the factorisation status
    bool storingFactors_;


    // =======================================================================
    // CCS DATA STRUCTURE FOR FACTOR L
    // =======================================================================

    //@{
    //! \name Data structure for storing L factor
    //! The lower triangular factor L of the (I)LU decomposition is stored
    //! in compressed column storage format. Additionally the entries in each
    //! column are stored in an ascending lexicographic ordering with respect
    //! to their row indices.
    //! \note The ones on the main diagonal of L are <b>not</b> stored. This
    //! must be taken into account e.g. in the forward substitution step.

    //! Array containing starting indices for the columns in the CCS format 
    std::vector<UInt> cptrL_;

    //! Array containing row indices of non-zero entries of L
    std::vector<UInt> ridxL_;

    //! Array containing values of non-zero entries of L
    std::vector<T> entryL_;

    //@}


    // =======================================================================
    // CRS DATA STRUCTURE FOR FACTOR U
    // =======================================================================

    //@{
    //! \name Data structure for storing U factor
    //! The lower upper factor U of the (I)LU decomposition is stored
    //! in compressed row storage format. Additionally the entries in each
    //! row are stored in an ascending lexicographic ordering with respect
    //! to their column indices.

    //! Array containing starting indices for the rows in the CRS format 
    std::vector<UInt> rptrU_;

    //! Array containing column indices of non-zero entries of U
    std::vector<UInt> cidxU_;

    //! Array containing values of non-zero entries of U
    std::vector<T> entryU_;
    //@}


    // =======================================================================
    // AUXILLIARY METHODS FOR CROUT FACTORISATION
    // =======================================================================

    //@{ \name Methods that depend on the storage format of the system matrix

    //! Copy k-th row of A into a dense vector

    //! This method copies the right part of the k-th row of the system matrix
    //! \f$A\f$, i.e. \f$\{a_{ki}\neq 0\,|\, k\leq i \leq n\}\f$ into the
    //! dense vector vecZ. The non-zero entries of the row are inserted into
    //! vecZ at the positions given by their respective column indices
    //! \f$i\f$.
    //! \param sysMat   (input)  system matrix in CRS storage format. Passing
    //!                          this parameter is currently not necessary,
    //!                          since we store the pointers to the CRS
    //!                          structure as attributes anyhow. However, for
    //!                          separating the format of the system matrix
    //!                          from the Crout this would be necessary.
    //! \param vecZ     (output) dense vector containing the copied non-zero
    //!                          entries
    //! \param k        (input)  index of the row to be copied
    //! \note We make the following assumptions:
    //! - the output vector vecZ contains only zeros on input
    //! - the system matrix is stored in OLAS' default CRS format, i.e. the
    //!   diagonal entry of a row is stored first in that row, followed by the
    //!   remaining off-diagonal entries in lexicographic ordering
    //! \return column index of right most non-zero entry in row
    UInt GetRowOfA( const CRS_Matrix<T> &sysMat, std::vector<T> &vecZ,
                    UInt k );

    //! Copy k-th column of A into a dense vector

    //! This method copies the lower part of the k-th column of the system
    //! matrix \f$A\f$, i.e. \f$\{a_{ik}\neq 0\,|\, k < i\leq n\}\f$ into the
    //! dense vector vecW. The non-zero entries of the column are inserted
    //! into vecW at the positions given by their respective row indices
    //! \f$i\f$.
    //! \param sysMat   (input)  system matrix in CRS storage format. Passing
    //!                          this parameter is currently not necessary,
    //!                          since we store the pointers to the CRS
    //!                          structure as attributes anyhow. However, for
    //!                          separating the format of the system matrix
    //!                          from the Crout this would be necessary.
    //! \param vecW     (output) dense vector containing the copied non-zero
    //!                          entries
    //! \param k        (input)  index of the column to be copied
    //! \note We make the following assumptions:
    //! - the method is called successively for k = 1,2,3, ... as is the case
    //!   in the Crout method
    //! - the output vector vecW contains only zeros on input
    //! - the system matrix is stored in OLAS' default CRS format, i.e. the
    //!   diagonal entry of a row is stored first in that row, followed by the
    //!   remaining off-diagonal entries in lexicographic ordering
    //! \return largest row index of non-zero entry in row
    UInt GetColOfA( const CRS_Matrix<T> &sysMat, std::vector<T> &vecW,
                    UInt k );

    //@}

    // Free all dynamically allocated memory
    void DeleteDynamicDataStructures();

    #ifdef DEBUG_CROUTLU
    //! Output content of an STL vector to screen (for debugging purposes)
    template<typename vecT>
    void PrintVector( const vecT &vec, const char *vecName ) const {
      (*debug) << " " << vecName << ": ";
      for ( UInt i = 0; i < vec.size(); i++ ) {
	(*debug) << vec[i] << " ";
      }
      (*debug) << std::endl;
    };
    #endif

  };

}

#ifndef EXPLICIT_TEMPLATE_INSTANTIATION
//#include "croutlu.cc"
#endif

#endif
