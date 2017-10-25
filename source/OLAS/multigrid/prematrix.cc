// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "OLAS/multigrid/prematrix.hh"


#include "OLAS/multigrid/hierarchylevel.hh"
#include "OLAS/multigrid/gaussseidel.hh"
#include "OLAS/multigrid/jacobi.hh"

#include "OLAS/solver/generatesolver.hh"
#include "OLAS/solver/BaseSolver.hh"
#include "OLAS/precond/IdPrecondStd.hh"



/**********************************************************/
#ifdef DEBUG_TO_CERR
#ifndef DEBUG_PREMATRIX
#define DEBUG_PREMATRIX
#endif // DEBUG_PREMATRIX
#define  debug  &std::cerr
#endif // DEBUG_PREMATRIX

/**********************************************************/
#define  PREMATRIX_ROW_GRANULARITY  32
/**********************************************************/

namespace CoupledField {
/**********************************************************/

template <typename T>
PreMatrix<T>::PreMatrix()
            : rows_( NULL ),
              cols_( NULL ),
              rowSize_( NULL ),
#ifndef PREMATRIX_USE_MODULO
              buffSize_( NULL ),
#endif
              numRows_( 0 ),
              hasRowsSorted_(false),
              hasDiagFirst_(false)
{
}


template <typename T>
PreMatrix<T>::PreMatrix( const UInt num_rows)
            : rows_( NULL ),
              cols_( NULL ),
              rowSize_( NULL ),
#ifndef PREMATRIX_USE_MODULO
              buffSize_( NULL ),
#endif
              numRows_( num_rows ),
              hasRowsSorted_(false),
              hasDiagFirst_(false)
{
 SetNumRows(num_rows);

}


template <typename T>
PreMatrix<T>::~PreMatrix()
{
    
    Reset();
}

/**********************************************************/

template <typename T>
UInt PreMatrix<T>::GetNumNonzeros() const
{
    UInt nnz = 0;
    for( UInt i = 0; i < numRows_; i++ )  nnz += rowSize_[i];
    return nnz;
}

/**********************************************************/

template <typename T>
UInt PreMatrix<T>::GetRowSize( const UInt row ) const
{

    return rowSize_[row];
}

/**********************************************************/

template <typename T>
const T* PreMatrix<T>::GetRowData( const UInt row ) const
{

    return rows_[row];
}

/**********************************************************/

template <typename T>
const UInt* PreMatrix<T>::
GetRowCols( const UInt row ) const
{

    return cols_[row];
}

template <typename T>
inline UInt PreMatrix<T>::
GetNumRows() const
{
    return numRows_;
}

/**********************************************************/

template <typename T>
T PreMatrix<T>::GetEntry( const UInt i,
                        const UInt j ) const
{
    const UInt* const cols = GetRowCols(i);
    for( UInt ij = 0; ij < GetRowSize(i); ij++ ) {
        if( cols[ij] == j ) {
            return GetRowData(i)[ij];
        }
    }
    return (T)0;
}

/**********************************************************/

template <typename T>
void PreMatrix<T>::SetEntry( const UInt  row, const UInt col,
                             const T& value )
{

    // find eventually allready present entry
    for( UInt ic = 0; ic < rowSize_[row]; ic++ ) {
        if( cols_[row][ic] == col ) {
            rows_[row][ic] = value;
            return;
        }
    }

    // entry not yet present
    CheckBufferSize( row );
    // add new entry
    rows_[row][++rowSize_[row]] = value;
    cols_[row][  rowSize_[row]] = col;
}

/**********************************************************/

template <typename T>
void PreMatrix<T>::AddToEntry( const UInt  row,
                               const UInt  col,
                               const T& value )
{
    if( row < 0 || row > GetNumRows() ){
      EXCEPTION("PreMatrix::AddToEntry recieved a wrong row index !!")
    }
    // find eventually allready present entry
    for( UInt ic = 0; ic < rowSize_[row] ; ic++ ) {
        if( cols_[row][ic] == col ) {
            rows_[row][ic] += value;
            return;
        }
    }

    // entry not yet present
    CheckBufferSize( row );
    // add new entry
    rows_[row][++rowSize_[row]] = value;
    cols_[row][rowSize_[row]] = col;

}

/**********************************************************/

template <typename T>
void PreMatrix<T>::SortDiagonal()
{

    T tdat;
    for( UInt i = 0; i < numRows_; i++ ) {
        // search the diagonal entry
        for( UInt ij = 0; ij < rowSize_[i]; ij++ ) {
            // place the diagonal entry at leading position
            if( cols_[i][ij] == i ) {
                // swap the column indices
                cols_[i][ij] = cols_[i][0];
                cols_[i][0]  = i;
                // swap the data
                tdat         = rows_[i][ij];
                rows_[i][ij] = rows_[i][0];
                rows_[i][0]  = tdat;
                // jump to next row
                break;
            }

        }
    }

    hasDiagFirst_ = true;
}

/**********************************************************/

template <typename T>
void PreMatrix<T>::Sort()
{
EXCEPTION("Attention! Sort() not yet ported from 1 to 0-based indexing !!!")
    hasDiagFirst_ = true;

    T tdat;
    for( UInt i = 0; i < numRows_; i++ ) {
//std::cerr
//<< "sorting row "<<i<<" / "<<numRows_<<": size = "<<rowSize_[i]<<std::endl;
        // search the diagonal entry
        UInt ij;
        for( ij = 0; ij < rowSize_[i]; ij++ ) {
            // place the diagonal entry at leading position
            if( cols_[i][ij] == i ) {
                // swap the column indices
                cols_[i][ij] = cols_[i][1];
                cols_[i][1]  = i;
                // swap the data
                tdat         = rows_[i][ij];
                rows_[i][ij] = rows_[i][1];
                rows_[i][1]  = tdat;
                // sort the remaining row using quick sort
                QuickSort( cols_[i]+2, rows_[i]+2, rowSize_[i]-1 );
                // next row
                break;
            }
        }
        // if we could not find a diagonal entry in this row,
        // simply sort the whole row
        if( ij > rowSize_[i] ) {
            QuickSort( cols_[i]+1, rows_[i]+1, rowSize_[i] );
            hasDiagFirst_ = false;
        }
    }

    hasRowsSorted_ = true;
}

/**********************************************************/

template <typename T>
void PreMatrix<T>::Reset()
{

    // delete array with row data
    if( rows_ != NULL ) {
        for( UInt i = 0; i < numRows_; i++ ) {
          delete [] ( rows_[i] );  rows_[i]  = NULL;
				}
        delete [] ( rows_ );  rows_  = NULL;
    }
    // delete array with column data
    if( cols_ != NULL ) {
        for( UInt i = 0; i < numRows_; i++ )
				{
          delete [] ( cols_[i] );  cols_[i]  = NULL;
				}
        delete [] ( cols_ );  cols_  = NULL;
    }
    // delete array with row sizes
    delete [] ( rowSize_ );  rowSize_  = NULL;
#ifndef PREMATRIX_USE_MODULO
    delete [] ( buffSize_ );  buffSize_  = NULL;
#endif
    // reset number of rows
    numRows_ = 0;
    
    hasDiagFirst_  = false;
    hasRowsSorted_ = false;
}



template <typename T>
std::ostream& PreMatrix<T>::Print( std::ostream& out ) const
{

    out << "\033[1mPreMatrix:\033[0m" << std::endl
        << "  "<<GetNumRows()<<" rows " << std::endl
        << "  "<<GetNumNonzeros()<<" entries: " << std::endl;
    
    for( UInt r = 0; r < GetNumRows(); r++ ) {
        for( UInt ic = 0; ic < GetRowSize(r); ic++ ) {
            out << "("<<r<<","<<GetRowCols(r)[ic]<<") = "
                << GetRowData(r)[ic] << std::endl;
        }
    }
    
    return out;
}

/**********************************************************/

template <typename T>
bool PreMatrix<T>::CreateArrays( const UInt num_rows )
{

    if( !num_rows ) { Reset();  return true; }

    // delete rows and the rows array, if it cannot be reused
    if( rows_ != NULL ) {
        for( UInt i = 0; i < numRows_; i++ ) {
            delete [] ( rows_[i] );  rows_[i]  = NULL;
				}
        if( num_rows != numRows_ ) {
            delete [] ( rows_ );  rows_  = NULL;
				}
    }
    // new rows array
    if( rows_ == NULL ) {
        NEWARRAY( rows_, T*, num_rows );
        for( UInt i = 0; i < num_rows; i++ )  rows_[i] = NULL;
    }

    // delete column indices and columns array, if it cannot be reused
    if( cols_ != NULL ) {
        for( UInt i = 0; i < numRows_; i++ ) {
            delete [] ( cols_[i] );  cols_[i]  = NULL;
				}
        if( num_rows != numRows_ ) {
            delete [] ( cols_ );  cols_  = NULL;
				}
    }
    // new column array
    if( cols_ == NULL ) {
        NEWARRAY( cols_, UInt*, num_rows );
        for( UInt i = 0; i < num_rows; i++ )  cols_[i] = NULL;
    }

    // new array for row sizes
    if( num_rows != numRows_ || !rowSize_) {
        delete [] ( rowSize_ );  rowSize_  = NULL;
        NEWARRAY( rowSize_, UInt, num_rows );
#ifndef PREMATRIX_USE_MODULO
        delete [] ( buffSize_ );  buffSize_  = NULL;
        NEWARRAY( buffSize_, UInt, num_rows );
#endif
    }
    // must be initialized or reset in any case
    for( UInt i = 0; i < num_rows; i++ )  rowSize_[i] = 0;
#ifndef PREMATRIX_USE_MODULO
    for( UInt i = 0; i < num_rows; i++ )  buffSize_[i] = 0;
#endif
    numRows_ = num_rows;

    return true;
}

/**********************************************************/

template <typename T>
void PreMatrix<T>::InstantiatePublicMethods()
{
    PreMatrix<T> dummyPreMatrix( 1 );
    GetNumNonzeros();
    GetNumRows();
    GetRowSize( 1 );
    GetRowData( 1 );
    GetRowCols( 1 );
    GetEntry( 1, 1 );
    SetNumRows( 1 );
    SetEntry( 1, 1, (T)1 );
    AddToEntry( 1, 1, (T)1 );
    SortDiagonal();
    Sort();
    Reset();
    Print( std::cout );
}

/**********************************************************/

template <typename T>
inline void PreMatrix<T>::CheckBufferSize( const UInt row )
{

#ifdef PREMATRIX_USE_MODULO
    // enlarge row, if necessary
    if( rowSize_[row] % PREMATRIX_ROW_GRANULARITY == 0 ) {
#else
    if( rowSize_[row] == buffSize_[row] ) {
#endif
        // enlarge data array
        T *newRow = NULL;
#ifdef PREMATRIX_USE_MODULO
        NEWARRAY( newRow, T,
                  rowSize_[row] + PREMATRIX_ROW_GRANULARITY );
#else
        NEWARRAY( newRow, T,
                  buffSize_[row] += PREMATRIX_ROW_GRANULARITY );
#endif
        for( UInt i = 0; i < rowSize_[row]; i++ ) {
            newRow[i] = rows_[row][i]; }
        delete [] ( rows_[row] );  rows_[row]  = NULL;
        rows_[row] = newRow;
        // enlarge column array
        UInt *newCols = NULL;
#ifdef PREMATRIX_USE_MODULO
        NEWARRAY( newCols, UInt,
                  rowSize_[row] + PREMATRIX_ROW_GRANULARITY );
#else
        NEWARRAY( newCols, UInt, buffSize_[row] );
#endif
        for( UInt i = 0; i < rowSize_[row]; i++ ) {
            newCols[i] = cols_[row][i]; }
        delete [] ( cols_[row] );  cols_[row]  = NULL;
        cols_[row] = newCols;
    }
}

/**********************************************************/

template <typename T>
void PreMatrix<T>::QuickSort(       UInt *const cols,
                                    T *const data,
                              const UInt        length )
{

    if( length <= 1 )  return;

    const UInt last         = length-1,
                  splitter     = cols[last];
    const T splitterData = data[last];

    UInt tcol, i = 0;
    T tdat;

    // splitt array
    for( UInt j = 0; j < last; j++ ) {
        if( cols[j] < splitter ) {
            // swap the column index
            tcol    = cols[j];
            cols[j] = cols[i];
            cols[i] = tcol;
            // swap the data
            tdat    = data[j];
            data[j] = data[i];
            data[i] = tdat;
            i++;
        }
    }
    // swap the splitter
    cols[last] = cols[i];
    cols[i]    = splitter;
    // and the corresponding data
    data[last] = data[i];
    data[i]    = splitterData;

    // quicksort the two parts
    QuickSort( cols,   data,   i++        );
    QuickSort( cols+i, data+i, length - i );
}


#undef  PREMATRIX_ROW_GRANULARITY
/**********************************************************/
} // namespace CoupledField

/**********************************************************
 * print-out operator
 **********************************************************/

namespace std {
template <typename T> std::ostream& operator <<
( std::ostream& out, const CoupledField::PreMatrix<T>& A ) {
    return A.Print( out ); }

}

/**********************************************************/
