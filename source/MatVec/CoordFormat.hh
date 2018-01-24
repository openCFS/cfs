#ifndef COORD_FORMAT_HH
#define COORD_FORMAT_HH

#include <vector>
#include <iterator>
#include <algorithm>
#include <cstdio>

#include "General/defs.hh"
#include "General/Exception.hh"

namespace CoupledField {

  //! Class for storing a sparse matrix in coordinate format

  //! This class allows storing a matrix in coordinate format, i.e. we store
  //! for every non-zero matrix entry its row and column index and value.
  //! This matrix class is intended as a pure container and offers no
  //! arithmetic methods. It is therefore not derived from the BaseMatrix
  //! class.\n\n
  //! In order to be memory efficient, it is possible to exploit the symmetry
  //! of a matrix and to store only its upper triangular part. In this case
  //! the constructor should be passed symmStorage=true.
  template<typename T>
  class CoordFormat {

  public:

    //! Enumeration data type describing storage layout of the matrix entries

    //! This enumeration data type describes the storage layout of the matrix
    //! entries in the CoordFormat. Currently there are the following values
    //! - NOLAYOUT
    //! - PURELEX
    typedef enum { NOLAYOUT, PURELEX } Layout;

    //! Struct describing a matrix entry by its two indices and its value
    typedef struct {
      UInt rowIndex;
      UInt colIndex;
      T value;
    } triple;

    //! Allowed Constructor
    CoordFormat( UInt nrows, UInt ncols, UInt nnz, bool symmStorage,
                 Layout myLayout = NOLAYOUT ) {


      nrows_ = nrows;
      ncols_ = ncols;
      nnz_   = nnz;

      assemblyDone_ = false;
      myLayout_     = myLayout;

      // Set flag for symmetry and if symmetric storage is requested
      // make sure that matrix is square
      symmStorage_  = symmStorage;
      if ( symmStorage == true && (ncols_ != nrows_) ) {
        EXCEPTION("You requested symmetric storage for a non-square "
                 << "matrix\n"
                 << "Dimensions are " << nrows_ << " x " << ncols_);
      }

      // Since we know the number of entries to be stored in advance,
      // we can "allocate" enough memory for the internal vector
      if ( symmStorage_ == true ) {
        entries_.reserve( (nnz_+nrows)/2 );
      }
      else {
        entries_.reserve( nnz_ );
      }

    }

    //! Constructor which imports the matrix based on a matrix market format file
    CoordFormat( std::string fname,
                 Layout myLayout = NOLAYOUT ){

      //ENTER_FCN( "ImportRealMM" );

         char singleChar;
         char myString[200];
         UInt nrows;
         UInt ncols;
         UInt nentries;
         UInt ridx;
         UInt cidx;
         Double value;
         symmStorage_ = false;

         assemblyDone_ = false;
         myLayout_     = myLayout;

         // Open input file and check for errors

         // Open input file and check for errors
         FILE *fp = fopen( fname.c_str(), "r" );
         if ( fp == NULL ) {
           EXCEPTION( "Cannot open file '" << fname << "' for reading");
         }

         // *******************************************************
         //   Read Matrix Market header and validate it is one of
         //
         //   %%MatrixMarket matrix coordinate real general
         //   %%MatrixMarket matrix coordinate real symmetric
         // *******************************************************
         fscanf( fp, "%c", &singleChar );
         if ( singleChar != '%' ) {
           EXCEPTION(  "Input file '" << fname << "' seems faulty.\n" //
              << " Problems reading its header!");
         }
         fscanf( fp, "%c", &singleChar );
         if ( singleChar != '%' ) {
           EXCEPTION( "Input file '" << fname << "' seems faulty.\n" //
              << " Problems reading its header!");
         }
         fscanf( fp, "%s", myString );
         if ( strcmp( myString, "MatrixMarket" ) != 0 ) {
           EXCEPTION( "Input file '" << fname << "' seems faulty.\n"
              << " Problems reading its header!");
         }
         fscanf( fp, "%s", myString );
         if ( strcmp( myString, "matrix" ) != 0 ) {
           EXCEPTION( "Input file '" << fname << "' seems faulty.\n"
              << "Contains a '" << myString << "' and not a 'matrix'");
         }
         fscanf( fp, "%s", myString );
         if ( strcmp( myString, "coordinate" ) != 0 ) {
           EXCEPTION("Input file '" << fname << "' seems faulty.\n"
              << "Format is '" << myString << "' and not 'coordinate'");
         }
         fscanf( fp, "%s", myString );
         if ( strcmp( myString, "real" ) != 0 ) {
           EXCEPTION("Input file '" << fname << "' seems faulty.\n"
              << "Entries are of type '" << myString << "' and not 'real'");
         }
         fscanf( fp, "%s", myString );
         if ( strcmp( myString, "general" ) == 0 ) {
           symmStorage_ = false;
         }
         else if ( strcmp( myString, "symmetric") == 0 ) {
           symmStorage_ = true;
         }
         else {
           EXCEPTION("Input file '" << fname << "' seems faulty.\n" //
              << "Format is '" << myString << "' and not 'general' " //
              << "or 'symmetric'");
         }

         // Finish this line (by reading characters until the '\n' appears
         char auxChar;
         do {
           fscanf( fp, "%c", &auxChar );
         } while ( auxChar != '\n' );

         // Read lines, until we have skipped comments
         do {
           fscanf( fp, "%c", &auxChar );
           if ( auxChar == '%' ) {
       do {
         fscanf( fp, "%c", &auxChar );
       } while ( auxChar != '\n' );
           }
           else {
       ungetc( auxChar, fp );
       break;
           }
         } while ( true );

         // This line (hopefully) contains the number of rows, columns and non-zeros
         fscanf( fp, "%u%u%u", &nrows, &ncols, &nentries );

         // Minimal consistency test
         if ( nentries > (nrows * ncols) ) {
           EXCEPTION("Input file '" << fname << "' seems faulty!\n" //
              << "#rows       = " << nrows << ",\n" //
              << "#cols       = " << ncols << ", but\n" //
              << "#numEntries = " << nentries);
         }
         else {
           std::cout << " Input file '" << fname << "' contains a matrix with"
            << "\n  + rows       = " << nrows
            << "\n  + cols       = " << ncols
            << "\n  + numEntries = " << nentries << std::endl;
         }


         UInt nnz = 0;
         if ( symmStorage_ == true ) {
           nnz = ( nentries - nrows ) * 2 + nrows;
           std::cout << "  + format     = symmetric\n";
         }
         else {
           nnz = nentries;
           std::cout << "  + format     = general\n";
         }
         std::cout << "  + nnz        = " << nnz << std::endl;
         this->nnz_ = nnz;
         this->nrows_ = nrows;
         this->ncols_ = ncols;
         entries_.reserve( nnz_ );

         // Read nonzero matrix entries and their indices
         for ( unsigned int k = 0; k < nentries; k++ ) {
           fscanf( fp, "%u%u%lf", &ridx, &cidx, &value );
           this->AddEntry( ridx, cidx, value );
         }

         // The matrix is assembled now
         this->FinaliseAssembly();

         // *********************
         //   Close output file
         // *********************
         if ( fclose( fp ) == EOF ) {
           WARN( "Could not close file %s after reading!" << fname );
         }
    }

    //! Default destructor
    ~CoordFormat() {
    }



    //! Add a matrix entry

    //! This adds the entry (i,k) to the matrix and assigns the value val to
    //! it. In the case of a symmetric matrix, if the symmStorage_ flag is
    //! true, all AddEntry requests for matrix entries in the lower triangular
    //! part will silently be ignored.
    //! \note
    //! - The entry is simply pushed onto the stack of matrix entries
    //!   that have already been added. No checking is performed, if it
    //!   already exists.
    //! - If tumble is set to true, then the method will swap the row and
    //!   column indices of the entry in the case of symmetric storage.
    void AddEntry( UInt i, UInt k, T val, bool tumble = false ) {

      UInt swap;

      if ( assemblyDone_ == true ) {
        EXCEPTION("Attempt to AddEntry to CoordFormat object after "
                 << "matrix assembly was finalised by a call to "
                 << "FinaliseAssembly!");
      }

      if ( symmStorage_ == true && i > k ) {
        if ( tumble == false ) {
          return;
        }
        else {
          swap = i;
          i = k;
          k = swap;
        }
      }

      triple newEntry;
      newEntry.rowIndex = i;
      newEntry.colIndex = k;
      newEntry.value = val;
      entries_.push_back( newEntry );

    }

    //! Finalises the assembly of the matrix

    //! Call this method after all matrix entries have been set. The object
    //! will then perform some simple consistency checks.
    //! \note
    //! - It is an error to call AddEntry after FinaliseAssembly was called.
    //! - It is an error to call any one of the access methods for matrix
    //!   entries before FinaliseAssembly was called.
    void FinaliseAssembly() {

      // Determine that number of stored entries is correct
      UInt numEntries = 0;
      if ( symmStorage_ == true ) {
        numEntries = (nnz_ + ncols_) / 2;
      }
      else {
        numEntries = nnz_;
      }
      if ( numEntries > entries_.size() ) {
        EXCEPTION("Matrix is stored incompletely in CoordFormat object!\n"
                 << " Expected " << numEntries << " entries, but got only "
                 << entries_.size());
      }
      else if  ( numEntries < entries_.size() ) {
    	EXCEPTION( "Matrix is stored incorrectly in CoordFormat object!\n"
                 << " Expected " << numEntries << " entries, but got "
                 << entries_.size());
      }

      // g++ reports syntax error before ;" in the for line
      // Check that all row and column indices are inside bounds
      // for ( std::vector<CoordFormat::triple>::iterator it=entries_.begin();
      //      it != entries_.end(); it++ ) {
      //  if ( it.ridx > nrows_ || it.ridx < 1 ) {
      //    EXCEPTION("Row index of entry out of bounds:\n"
      //           << "row index    = " << it.ridx
      //           << " should be in [1," << nrows_ << "]\n"
      //           << "column index = " << it.cidx
      //             << " should be in [1," << ncols_ << "]\n"
      //             << " value       = " << it.value);
      //  }
      //  if ( it.cidx > ncols_ || it.cidx < 1 ) {
      //    EXCEPTION("Column index of entry out of bounds:\n"
      //             << "row index    = " << it.ridx
      //             << " should be in [1," << nrows_ << "]\n"
      //             << "column index = " << it.cidx
      //             << " should be in [1," << ncols_ << "]\n"
      //             << " value       = " << it.value);
      //  }
      //}

      // Object passed all tests
      assemblyDone_ = true;
    }

    //! Query number of matrix rows
    UInt GetNumRows() const {
      return nrows_;
    }

    //! Query number of matrix columns
    UInt GetNumCols() const {
      return ncols_;
    }

    //! Query number of non-zero matrix entries
    UInt GetNNZ() const {
      return nnz_; 
    }

    //! Query number of stored matrix entries
    UInt GetNumEntries() const {
      return entries_.size();
    }

    //! Query status of flag for symmetric storage
    bool GetSymmStorageFlag() const {
      return symmStorage_;
    }

    //! Read only access to row index vector
    UInt ridx( UInt i ) const {
      if ( assemblyDone_ == false ) {
        EXCEPTION( "Attempt to access matrix entry, before assembly was "
                 << "finalised by a call to FinaliseAssembly!" );
      }
      return entries_[i].rowIndex;
    }

    //! Read only access to column index vector
    UInt cidx( UInt i ) const {
      if ( assemblyDone_ == false ) {
        EXCEPTION( "Attempt to access matrix entry, before assembly was "
                 << "finalised by a call to FinaliseAssembly!" );
      }
      return entries_[i].colIndex;
    }

    //! Read only access to matrix entry vector
    T val( UInt i ) const {
      if ( assemblyDone_ == false ) {
        EXCEPTION( "Attempt to access matrix entry, before assembly was "
                 << "finalised by a call to FinaliseAssembly!" );
      }
      return entries_[i].value;
    }

    //! Re-sort entries to conform to a special layout

    //! This method can be used to re-order the storage layout of the matrix
    //! entries, such that it conforms to certain requirements. If the matrix
    //! is already stored in that form, calling the method does not incur any
    //! operations.
    //! \param newLayout Speficies the desired layout
    void Sort( Layout newLayout ) {

      // Test, if we need to re-sort
      if ( myLayout_ == newLayout ) {
        return;
      }

      // We need to re-sort
      switch( newLayout ) {

        // Nothing to be done
      case NOLAYOUT:
        break;

        // Pure lexicographic layout
      case PURELEX:
        // std::sort( entries_.begin(), entries_.end(), PureLex );
        std::sort( entries_.begin(), entries_.end(), CoordFormat::PureLex );
        break;

        // Should not occur
      default:
        EXCEPTION( "Request for unsupported layout type" );
      }
    }

  private:

    //! Storage layout of matrix entries

    //! This attribute keeps track of the layout in which the matrix entries
    //! are stored
    Layout myLayout_;

    //! Number of matrix rows
    unsigned int nrows_;

    //! Number of matrix columns
    unsigned int ncols_;

    //! Number of non-zero matrix entries
    unsigned int nnz_;

    //! Vector containing the matrix entries
    std::vector<triple> entries_;

    //! Flag indicating whether the matrix is symmetric or not

    //! This flag signals that the matrix is symmetric and that this symmetry
    //! is exploited by storing only its upper triangle. In this case the
    //! of number of stored matrix entries is given by entries and not by nnz.
    //! If the flag is not set, then the matrix might still be symmetric,
    //! however, all entries are stored.
    bool symmStorage_;

    //! Flag indicating whether assembly of matrix was finalised
    bool assemblyDone_;

    //! Use of the default constructor is forbidden.
    CoordFormat() {};

    //! Comparison operator of for a pure lexicographic ordering

    //! This comparison operator can be used to sort the matrix entries, such
    //! that they are stored in a pure lexicographic ordering. If (i,j) is the
    //! index pair of matrix entry A and (l,m) that of matrix entry B, then
    //! PureLex(A,B) will give true, if either (i<l) or if (i=l) and
    //! (j<m).
    static bool PureLex( const triple& a, const triple& b ) {

      bool retVal = false;

      if ( a.rowIndex < b.rowIndex ) {
        retVal = true;
      }
      else if ( a.rowIndex > b.rowIndex ) {
        retVal = false;
      }
      else {
        if ( a.colIndex < b.colIndex ) {
          retVal = true;
        }
        else {
          retVal = false;
        }
      }
      return retVal;
    }
  };

}

#endif
