// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef OLAS_VECTOR_HH
#define OLAS_VECTOR_HH

#include <string>
#include <iostream>
#include <fstream>
#include <cmath>

#include "utils/utils.hh"
#include "matvec/typedefs.hh"
#include "matvec/opdefs.hh"
#include "matvec/SparseVector.hh"

namespace OLAS {

  //! Class for dense array-based vector
  //! dof can be specified by using T = tvmet::Vector<T, dof>

  template <typename T>
  class Vector : public SparseVector {

  public:

    //! Typename of matrix entries compatible with T
    typedef typename assocType<T>::T_Mtype T_Mtype;

    //! Typename of the vector entries (=T)
    typedef typename assocType<T>::T_Vtype T_Vtype;

    //! Scalar of the same primitive data type as matrix
    typedef typename assocType<T>::T_Stype T_Stype;


    // =======================================================================
    // CONSTRUCTION, DESTRUCTION AND REPLACEMENT
    // =======================================================================

    //! \name Methods for construction, destruction and replacement
    //! The vector class can also be used as a form of container. We
    //! provide methods that allow to add to a simple data array the
    //! functionality of this vector class, especially its arithmetic and
    //! access operations.

    //@{

    //! Default Constructor
    Vector() : SparseVector(), data_(NULL) {
      dof_ = BlockSize<T>::size;
      size_ = 0;
      dataSize_ = 0;

      // We set this to true be default, to be able to do a re-size
      // afterwards. If a Replace() is done, it will adapt it accordingly.
      memBelongsToMe_ = true;
    }

    //! Copy Constructor

    //! This is a deep copy constructor. It will allocate its own data_
    //! array and generate an actual copy of the entries of the original
    //! vector. Consequently it sets the memBelongsToMe_ attribute to true.
    //! This constructor is mainly provided to allow for certain sparse
    //! arithmetic operations of the SBM_Vector class.
    Vector( const Vector<T> &origVec ) {


      dof_ = BlockSize<T>::size;

      // Obtain size info and allocate memory
      size_     = origVec.size_;
      dataSize_ = size_;
      NewArray( data_, T, size_ );
      memBelongsToMe_ = true;

      // Copy vector entries
      for ( UInt i = 1; i <= size_; i++ ) {
        data_[i] = origVec.data_[i];
      }
    }

    //! Constructor for vector of specified length
    
    //! This constructor generates a vector of the given length. If the input
    //! parameter mem is true, then memory will be allocated and the vector
    //! initialised with zeros.
    Vector( Integer length );

    //! Destructor

    //! The default destructor must be deep, i.e. it must free all dynamically
    //! allocated memory.
    ~Vector();
    
    //! Return the Entry type of the vector
      
    //! The method returns the entry type of the vector (i.e. Double, Complex,
    //! or whatever it is). This is encoded as a value of the enumeration data
    //! type MatrixEntryType.
    MatrixEntryType GetEntryType() const {
      return  EntryType<T>::M_EntryType;
    }

    //! Re-size the vector

    //! This method can be used to change the length of the vector.
    //! \param newSize the new length of the vector
    //! \param init    if this boolean is yes, then the vector entries
    //!                will be set to zero by a call to Init().
    //! \note
    //! - A re-allocation of memory will only be triggered by Resize(),
    //!   if the new length of the vector exceeds the length of the
    //!   internal data array as given by dataSize_.
    //! - If init = false, then the vector entries will be undefined
    //!   if re-allocation was performed. Otherwise they will retain
    //!   their old values.
    //! - Re-size will currently refuse to perform a re-size operation,
    //!   if it is not responsible for the memory management of the data_
    //!   array.
    void Resize( UInt newSize, bool init = false );

    //! Add functionality of vector class to a data array

    //! This method allows to add the functionality of the Vector class,
    //! especially its arithmetic and access methods, to a plain data array.
    //! Calling this method will replace the internal data_ array by the
    //! entries vector and re-set the internal attributes, like e.g. size_.
    //! Note that, since the vector is one-based, the pointer to the new
    //! memory array must have an appropriate off-set and that length may
    //! not be larger than the actual length of the allocated memory block.
    //! Calling Replace is memory safe in the sense that the old data_ array
    //! will be de-allocated, if this is the responsibility of the
    //! corresponding vector object. Responsibility for de-allocating the
    //! new data_ array is transferred to the vector object by setting the
    //! transferMem parameter to true.
    //! \param length      the length of the new vector
    //! \param entries     pointer (one-based) to the data array containing
    //!                    the new vector entries
    //! \param transferMem flag signalling transfer of responsibility for
    //!                    memory management
    void Replace( UInt length, T_Vtype* entries, bool transferMem );

    //! Withdraw responsibility for memory management from vector object

    //! Calling this method will relieve the vector object from the task of
    //! managing the memory block containing the entries of the vector.
    //! It will also return a (one-based) pointer to that memory block.
    //! \note Responsibility for de-allocating the memory should belong to
    //!       a unique object. Therefore, DecoupleMem will issue a warning
    //!       when an attempt is made to obtain that responsibilty from the
    //!       vector object and the latter does not own it.
    T_Vtype* DecoupleMem() {
      if ( memBelongsToMe_ == false ) {
	(*warning) << "DecoupleMem was called on a vector object not "
		   << "responsible for managing its memory block! Memory "
		   << "problems may arise!";
	Warning( __FILE__, __LINE__ );
      }
      memBelongsToMe_ = false;
      return data_;
    }

    //! Clear the vector

    //! Calling this method will clear the vector, i.e. the internal data_
    //! array will be de-allocated, if it belongs to the vector object, and
    //! the internal attributes will be re-set to the state we also obtain
    //! from the default constructor.
    void Clear() {

      if ( memBelongsToMe_ == true ) {
        DeleteArray( data_ );
      }

      dof_      = BlockSize<T>::size;
      size_     = 0;
      dataSize_ = 0;
      data_     = NULL;

      // We set this to true be default, to be able to do a re-size
      // afterwards. If a Replace() is done, it will adapt it accordingly.
      memBelongsToMe_ = true;
    }

    //@}

    // =======================================================================
    // ARITHMETIC OPERATIONS
    // =======================================================================

    //! \name Arithmetic operations

    //@{
    //! Add vec to this vector object
    void Add(const SparseVector &vec);

    //! Add a scaled version of a vector to this vector object

    //! The method takes this vector object \f$x\f$ and replaces it with
    //! the result of \f$x + \alpha v\f$.
    void Add(T_Stype a, const SparseVector &vec);   

    //! Replace this vector object by the sum of two scaled vectors

    //! This method replaces this vector object by the sum
    //! \f$\alpha x +\beta y\f$.
    void Add(T_Stype a, const SparseVector &vec1,  
             T_Stype b, const SparseVector &vec2);         
                                                
    //! Compute inner product

    //! The method computes the value of the inner product between this vector
    //! and the input vector vec. The value is returned in sum.
    void Inner(const SparseVector& vec, T_Stype& sum) const;

    //! Compute Euclidean norm of this vector object
    Double NormEuclid() const;

    //! Same as the BLAS functions of the same name

    //! The method assumes that this vector is x and performs the
    //! classical BLAS function AXPY, i.e. it scales the vector x by the
    //! factor alpha and adds the vector y to it. The result will over-write
    //! the vector x.
    void Axpy( const T_Stype alpha, const SparseVector &y );

    //! Overload assignment operator
    SparseVector &operator= ( const SparseVector &stdvec );

    //! Divide each vector entry by specified real-valued scalar
    void ScalarDiv( const Double factor );

    //! Multiply each vector entry by specified real-valued scalar
    void ScalarMult( const Double factor );

    //! Divide each vector entry by specified complex-valued scalar
    void ScalarDiv( const Complex factor );

    //! Multiply each vector entry by specified complex-valued scalar
    void ScalarMult( const Complex factor );

    //@}

    // =======================================================================
    // I/O OPERATIONS
    // =======================================================================

    //! \name Methods for I/O operations

    //@{
    /** Prints the content for Logging
     * @return comma seperated list from index 1 to <= size! */
    std::string ToString() const;
    
    /** Print vector to output stream os
     * @depreciated to be replaces by ToString and logging */
    void Print(std::ostream &os) const;

    //! Export vector to file

    //! This method can be used to export the vector to an ascii file. The
    //! format is extremely simple. If the vector is of dimension
    //! \f$n\times 1\f$, then the output file will contain \f$n+1\f$ rows.
    //! The first row contains the dimension \f$n\f$, while the remaining
    //! rows contain the vector's entries, so row (k+1) contains entry
    //! \f$a_k\f$.
    virtual void Export( const Char *fname ) const;

    //@}

    // =======================================================================
    // OBTAIN / MANIPULATE MATRIX ENTRIES
    // =======================================================================

    //! \name Obtain/manipulate matrix entries
    //@{

    //! Set all vector entries to zero

    //! All entries of the vector are set/initialised to zero (on a scalar
    //! level).
    void Init();

    //! Return a reference to i-th entry
    inline T_Vtype & operator[](Integer i) {
      return data_[i];
    };

    //! Return a reference to i-th entry
    inline const T_Vtype& operator[](Integer i) const {
      return data_[i];
    };

    //! Set the value of a vector entry

    //! This method sets the entry of the vector at position i to the
    //! specified value.
    void SetVectorEntry( Integer i, const T_Vtype &val );

    //! Get a reference to a vector entry

    //! This method sets val to the value of vector entry i
    void GetVectorEntry( Integer i, T_Vtype &val ) const;

    //! Add val to the value of a vector entry

    //! This method adds to the entry of the vector at position i the
    //! specified value val.
    void AddToVectorEntry( Integer i, const T_Vtype &val );

    //@}

    // =======================================================================
    // DIRECT ACCESS TO INTERNAL DATA STRUCTURE
    // =======================================================================

    //! \name Direct access to internal data structure
    //! The following methods allow a direct manipulation (read or write) of
    //! the internal data structure. Though this does violate the concept of
    //! data encapsulation, these methods are provided for the sake of
    //! efficiency.

    //@{

    //! Return pointer to internal data array data_
    T_Vtype* GetPointer(){
      return data_;
    }

    //! Return pointer to internal data array data_ (const)
    const T_Vtype* GetPointer() const{
      return data_;
    }

    //@}

    // =======================================================================
    // MISCELLANEOUS METHODS
    // =======================================================================

    //! \name Miscellaneous methods
    //@{

    //! Query size of internal data array
    UInt GetDataSize() const {
      return dataSize_;
    }

    //! Query tag for matrix entries compatible with vector
    std::string GetTagM() {
      return assocType<T>::tagM;
    }
 
    //! Query tag for vector entries compatible with vector
    std::string GetTagV() {
      return assocType<T>::tagV;
    }

    //! Query tag for scalar entries compatible with vector
    std::string GetTagS() {
      return assocType<T>::tagS;
    }

    //! Method to force instantiation of members in templated vectors

    //! This method is used by the GenerateVectorObject function to force the
    //! instantiation of the public member functions of this templated vector
    //! class.
    //! \todo Complete the implementation of this method.
    void InstantiatePublicMethods();

    //@}

  protected:

    //! 1-based array storing the vector entries.
    T_Vtype *data_;

    //! Actual size of the internal data array.

    //! This attribute is used to keep track of the size of the internal
    //! data_ array. Note that dataSize_ will initially be identical to
    //! size_, but might start to differ from it due to re-size operations.
    UInt dataSize_;

    //! Flag signaling whether management of data array is done by this object

    //! This attribute is used to keep track on the fact whether the object
    //! is responsible for managing the memory of the data_ array, especially
    //! its deallocation.
    bool memBelongsToMe_;

  };


}

#endif // OLAS_VECTOR_HH
