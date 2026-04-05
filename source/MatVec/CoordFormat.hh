#ifndef COORD_FORMAT_HH
#define COORD_FORMAT_HH

#include <vector>
#include <iterator>
#include <algorithm>
#include <iostream>

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
    CoordFormat( UInt nrows, UInt ncols, UInt nnz, bool symmStorage, Layout myLayout = NOLAYOUT );

    //! Constructor which imports the matrix based on a matrix market format file
    CoordFormat( const std::string& fname, Layout myLayout = NOLAYOUT );
    
    //! Default destructor
    ~CoordFormat() { }

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
    inline void AddEntry( UInt i, UInt k, T val, bool tumble = false ) {

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
    void FinaliseAssembly();

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
    inline void Sort( Layout newLayout ) {

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
    CoordFormat() {}

    //! Comparison operator of for a pure lexicographic ordering

    //! This comparison operator can be used to sort the matrix entries, such
    //! that they are stored in a pure lexicographic ordering. If (i,j) is the
    //! index pair of matrix entry A and (l,m) that of matrix entry B, then
    //! PureLex(A,B) will give true, if either (i<l) or if (i=l) and
    //! (j<m).
    static bool PureLex( const triple& a, const triple& b );
}; // class CoordFormat

} // namespace CoupledField 

#endif
