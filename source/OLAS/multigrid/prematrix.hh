// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

/* $Id$ */

#ifndef _PREMATRIX_HH_
#define _PREMATRIX_HH_

/**********************************************************
 * The size of the buffer used for each row in a prematrix grows
 * dynamically with the number of entries in a row. It grows in steps
 * of size PREMATRIX_ROW_GRANULARITY (see prematrix.cc). To determine,
 * when the row buffer must be enlarged, there exist two possibilities
 * in this class. Either by storing the current buffer size in an
 * additional array or by testing the current row size with a modulo
 * operation against PREMATRIX_ROW_GRANULARITY, which saves some memory,
 * but might be slightly slower in heavy usage of setting new entries.
 * To switch between these two strategies, define the preprocessor flag
 * PREMATRIX_USE_MODULO in order to use the modulo strategy, or leave
 * it undefined to stick to the first version.
 **********************************************************/
// #define PREMATRIX_USE_MODULO
/**********************************************************/

#include "matvec/matvec.hh"

namespace CoupledField {
/**********************************************************/

//! matrix class with dynamic memory allocation

/*! This class can be used to construct CRS or other matrices without
 *  previously constructing the connectivity graph. It is simply an
 *  array of dynamically (in chunks) growing rows. These rows are
 *  organized like the rows in a CRS matrix, so that they can be copied
 *  direktly into the CRSMatrix object for conversion.
 */
template <typename T>
class PreMatrix
{
    public:

        //! create an empty matrix
        PreMatrix();
        //! create a new matrix with num_rows rows
        PreMatrix( const Integer num_rows );
        
        //! destroy the matrix
        ~PreMatrix();

        //! returns the number of set (i.e. nonzero) entries
        inline Integer GetNumNonzeros() const;        
        //! returns the number of rows
        inline Integer GetNumRows() const;        
        //! returns the number of non-zero entries in the i-th row
        inline Integer GetRowSize( const Integer row ) const;
        /*! \brief returns a pointer to the data of the i-th row
         *  (one-based indexed).
         */
        inline const T* GetRowData( const Integer row ) const;
        /*! \brief returns a pointer to the columns of the i-th row
         *  (one-based indexed)
         */
        inline const Integer* GetRowCols( const Integer row ) const;
        //! returns the entry (i,j) of the matrix
        inline T GetEntry( const Integer i,
                                 const Integer j ) const;

        /*! \brief sets the number of rows. Independend on previously
         *  contained data the result of a successfull call will be
         *  an initialized empty matrix. Therefore you do not need to
         *  invoke PreMatrix::Reset before.
         */
        inline bool SetNumRows( const Integer nrows ) {
            return CreateArrays( nrows );
        }

        //! sets the matrix entry (row,col) to value "value"
        void SetEntry( const Integer  row, const Integer col,
                       const T& value );
        //! adds "value" to entry (row,col)
        void AddToEntry( const Integer  row, const Integer col,
                         const T& value );

        //! places the diagonal entries at leading position
        void SortDiagonal();
        //! sort the rows of the pre-matrix (leading diagonal entry)
        void Sort();
        //! resets the whole class
        void Reset();
        
        /*! \name methods concerning the sorting of matrix entries */
        //@{
        //! find out if the matrix conforms to our CRS layout
        
        /*! returns true if all rows are sorted ascendingly and
            the diagonal element of each row is stored in the first
            position
        */
        bool GetLayoutFlag() const { return hasDiagFirst_ && hasRowsSorted_; }
        //! returns true if all rows are sorted ascendingly
        bool HasRowsSorted() const { return hasRowsSorted_; }
        //! returns true if the diagonal elements in each row
        //! are stored in the first position
        bool HasDiagFirst() const { return hasDiagFirst_; }
        //! sets the diagonal flag directly (use this method carefully!)
        void SetDiagFlag( const bool flag = true ) { hasDiagFirst_ = flag; }
        //! sets the sorted flag directly (use this method carefully!)
        void SetSortedFlag( const bool flag = true ) { hasRowsSorted_ = flag; }
        //@}

#ifdef PREMATRIX_ENABLE_IMPORT
        //! imports a matrix from a matrix market ASCII format

        /*! This function imports a matrix from a file in matrix market
         *  ASCII format, like for example exported by CRS_Matrix::Export.
         *  The standard instantiation of this methods works for matrices
         *  with entry type \c Double. To create further specializations,
         *  e.g. for \c Complex, specialize the auxiliary method
         *  PreMatrix::ImportEntryASCII.
         *  \param filename name of the matrix file
         *  \return \c true, if the import was successfull, otherwise
         *          \c false
         */
        bool ImportASCII( const char* const filename );
#endif

        //! prints out the matrix
        std::ostream& Print( std::ostream& out ) const;

        //! method to force instantiation (should never be invoked)
        void InstantiatePublicMethods();

    protected:

        //! creates and initializes the needed memory

        /*! Creates and initializes the memory needed for a matrix
         *  with num_rows rows. If the PreMatrix object already
         *  contains some entries, the are removed. Present memory
         *  arrays are reused, if possible.
         */
        bool CreateArrays( const Integer num_rows );
        //! enlarges the row buffer for the passed row, if necessary
        inline void CheckBufferSize( const Integer row );
        //! quick sort for sorting the rows
        static void QuickSort(       Integer *const cols,
                                     T *const data,
                               const Integer        length );

#ifdef PREMATRIX_ENABLE_IMPORT
        //! auxiliary function for ImportASCII, which imports one single entry

        /*! This auxiliary function is exclusively used by method
         *  PreMatrix::ImportASCII and imports one single entry from
         *  a line containing row index, column index, value. The standard
         *  implementation expects the entry to be one \c Double value.
         *  To adapt ImportASCII to different entry types, specialize
         *  method ImportEntryASCII. ImportEntryASCII should write the
         *  row and column indices and the matrix entry itself into the
         *  passed parameters.
         *  \param line buffer containing the matrix entry line
         *  \param row row index
         *  \param col column index
         *  \param matix entry
         *  \return \c true, if the matrix entry could be read successfully,
         *          \c false otherwise
         */
        bool ImportEntryASCII( const char* const line,
                               Integer &row, Integer &col, T &val );
#endif
        

        T  **rows_;     //!< array of pointers to the rows
        Integer  **cols_;     //!< array of pointers to the column indices
        Integer   *rowSize_;  //!< array containing current row sizes
#ifndef PREMATRIX_USE_MODULO
        Integer   *buffSize_; //!< array containing the current buffer length
#endif
        Integer    numRows_;  //!< number of rows
        bool       hasRowsSorted_; //!< if true the rows are sorted
        bool       hasDiagFirst_; //!< is the diag elem first in each row?

};    

/**********************************************************/
} // namespace CoupledField

#endif // _PREMATRIX_HH_
