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

#include "OLAS/multigrid/transfer.hh"


namespace CoupledField {
/**********************************************************/

//! matrix class with dynamic memory allocation

/*! This class can be used to construct CRS or other matrices without
 *  previously constructing the connectivity graph. It is simply an
 *  array of dynamically (in chunks) growing rows. These rows are
 *  organized like the rows in a CRS matrix, so that they can be copied
 *  directly into the CRSMatrix object for conversion.
 */
template <typename T>
class PreMatrix
{
    public:

        //! create an empty matrix
        PreMatrix();
        //! create a new matrix with num_rows rows
        PreMatrix( const UInt num_rows );
        
        //! destroy the matrix
        ~PreMatrix();

        //! returns the number of set (i.e. nonzero) entries
        UInt GetNumNonzeros() const;
        //! returns the number of rows
        UInt GetNumRows() const;
        //! returns the number of non-zero entries in the i-th row
        UInt GetRowSize( const UInt row ) const;
        /*! \brief returns a pointer to the data of the i-th row
         *  (one-based indexed).
         */
        const T* GetRowData( const UInt row ) const;
        /*! \brief returns a pointer to the columns of the i-th row
         *  (one-based indexed)
         */
        const UInt* GetRowCols( const UInt row ) const;
        //! returns the entry (i,j) of the matrix
        T GetEntry( const UInt i,
                                 const UInt j ) const;

        /*! \brief sets the number of rows. Independent on previously
         *  contained data the result of a successful call will be
         *  an initialized empty matrix. Therefore you do not need to
         *  invoke PreMatrix::Reset before.
         */
        inline bool SetNumRows( const UInt nrows ) {
            return CreateArrays( nrows );
        }

        //! sets the matrix entry (row,col) to value "value"
        void SetEntry( const UInt  row, const UInt col,
                       const T& value );
        //! adds "value" to entry (row,col)
        void AddToEntry( const UInt  row, const UInt col,
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
        bool CreateArrays( const UInt num_rows );
        //! enlarges the row buffer for the passed row, if necessary
        inline void CheckBufferSize( const UInt row );
        //! quick sort for sorting the rows
        static void QuickSort(       UInt *const cols,
                                     T *const data,
                               const UInt        length );


        T  **rows_;     //!< array of pointers to the rows
        UInt  **cols_;     //!< array of pointers to the column indices
        UInt   *rowSize_;  //!< array containing current row sizes
#ifndef PREMATRIX_USE_MODULO
        UInt   *buffSize_; //!< array containing the current buffer length
#endif
        UInt    numRows_;  //!< number of rows
        bool       hasRowsSorted_; //!< if true the rows are sorted
        bool       hasDiagFirst_; //!< is the diag elem first in each row?

};    

/**********************************************************/
} // namespace CoupledField

#endif // _PREMATRIX_HH_
