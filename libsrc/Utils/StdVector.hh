#ifndef FILE_STDVECTOR_2004
#define FILE_STDVECTOR_2004

#include "General/environment.hh"
#include "Utils/tools.hh"
#include "DataInOut/WriteInfo.hh"

namespace CoupledField {


// forward class declaration
template <class TYPE> class StdVectorListInitializer;

//! General template class for storing elements of same type
//! in sequential order. In contrast to the CFS-Vector, this class
//! requires no arithmetic functionality for its elements.
//! This class replaces the use of std::vector<>-class and therefore
//! provides a similar interface
template<class TYPE>
class StdVector{
public:

  //! Constructor
  StdVector();

  //! Constructor with inital size.
  //! All entries are filled with zeroes
  StdVector(Integer size);

  //! Copy constructor
  StdVector(const StdVector<TYPE> & vec);

  //! Copy constructor with std::vector
  StdVector(const std::vector<TYPE> & vec);

  //! Destructor
  ~StdVector();  

  //! Deletes all entries in the vector
  //! and frees its memory
  void Clear();

  //! Initalizes the vector with a given entry
  /*!
    \param entry (input) Entry vector gets inialized with
  */
  //! \note this method does not change the size of the vector!
  void Init(const TYPE entry = TYPE());

  //! True, if vector is empty
  inline Boolean IsEmpty() const {return (size_? FALSE : TRUE);}
  
  //! Returns capacity of the vector
  inline Integer Capacity() const {return capacity_;}

  //! 
  void Reserve(Integer size);

  //! Get the length of the vector
  inline Integer GetSize() const {return size_;}

  //! Set the lenght of the vector
  /*!
    \param size (input) Lengh of vector
  */
  //! \note the entries are set to zero afterwards!
  void Resize(const Integer size);
 
  //! Overloading of operation =
  StdVector	&operator=	(const StdVector &);

  //! Overloading for operator=

  //! This method is needed to intialize a StdVector like this
  //! \verb StdVector<Integer> A;
  //! \verb A = 1,2,5,10
  inline StdVectorListInitializer<TYPE> operator=(const TYPE x);

  //! Build vector from std::vector
  StdVector & operator= (const std::vector<TYPE> & vec);
  
  //! General access operator
  inline TYPE &operator[] (const Integer i);

  //! Read-Only access operator
  inline TYPE operator[] (const Integer i) const;

  //! Return pointer p to array 
  inline TYPE*  GetPointer() const;

  //! Add element of the same type at position pos
  void Insert(const TYPE & y, Integer pos);

  //! Add element of the same type at position pos
  //! numCopies times
  void Insert(const TYPE & y, Integer pos, Integer numCopies);

  //! Add element of the same type at the end of the vector
  void Push_back(const TYPE & y);

  //! Delete element from vector on position pos
  void Erase (const Integer pos);

  //! Delete elements from position pos1 to pos2, on pos1, pos2 too
  void Erase (const Integer pos1, const Integer pos2);


  //! Return the position number of element x in the vector

  //! Finds the element x in the vector and returns the 
  //! position. If no element was found, it returns -1.
  Integer Find(const TYPE &x) const;

  //! Overloading of operation equal for Vector
  Boolean operator== (const StdVector &) const;
  
  //! Overloading of operation not equal for Vector
  Boolean operator!= (const StdVector &) const;
  

  // ******************************************************
  // MISCELANEOUS FUNCTIONS
  // ******************************************************
  
  //! Sort of vector: v - vec.p, n - vec.size
  template <class S> void Sort(S* v, Integer n);
  
  //! Swap 2 elements in vector Ex Swap(v[i],v[j])
  template<class T2> void Swap(T2& a, T2 & b);

 protected:

  //! Length of the vector
  Integer size_;

  //! Data of the vector
  TYPE* data_;
  
  //! Capacity of the vector
  Integer capacity_;

};

  // ******************************************************
  // HELPER CLASS FOR INITALIZING StdVector
  // (ref. 'Techniques for Scientific C++' 
  // by Todd Veldhuizen, page. 43ff
  // ******************************************************
  template <class TYPE>
  class StdVectorListInitializer
  {
  public:

    //! Constructor
    StdVectorListInitializer(StdVector<TYPE> * vec)
      :vec_(vec) {};
    
    //! Overloading of comma operator
    StdVectorListInitializer<TYPE> operator, (const TYPE x)
    {
      vec_->Push_back(x);
      return StdVectorListInitializer<TYPE>(vec_);
    }

  private:
    //! pointer to vector
    StdVector<TYPE> * vec_;
  };


  // ******************************************************
  // INLINE FUNCTION DECLARATION
  // ******************************************************
  
  //! Element can be referred to as v[i]
  template<class TYPE>
  TYPE & StdVector<TYPE>::StdVector::operator[] (const Integer i)
  {	
#ifdef CHECK_INDEX
    if (i >= size_){
      std::string errorMsg;
      errorMsg =  "Vector: invalid access to element ";
      errorMsg += Info->GenStr(i+1);
      errorMsg += "! \n Length of vector: ";
      errorMsg += Info->GenStr(size_);
      Error(errorMsg.c_str(),__FILE__,__LINE__);
    }
#endif
    return  data_[i];
  }

  //! Element can be referred to as v[i]
  template<class TYPE>
  TYPE StdVector<TYPE>::operator[] (const Integer i) const
  {	
#ifdef CHECK_INDEX
    if (i >= size_){
      std::string errorMsg;
      errorMsg =  "Vector: invalid access to element ";
      errorMsg += Info->GenStr(i+1);
      errorMsg += "! \n Length of vector: ";
      errorMsg += Info->GenStr(size_);
      Error(errorMsg.c_str(),__FILE__,__LINE__);
    }
#endif
    return  data_[i];
  }
  

  //! 
  template<class TYPE>
  StdVectorListInitializer<TYPE> StdVector<TYPE>::operator= (const TYPE x)
  {
    Clear();
    Push_back(x);
    return StdVectorListInitializer<TYPE>(this);
  }

  //! Return pointer p to array 
  template<class TYPE>
  TYPE* StdVector<TYPE>::GetPointer() const
  {
#ifdef CHECK_MEMORY
      if (!data_)
	Error("Vector: undefined Vector",__FILE__,__LINE__);
#endif
      return data_;
  }


// ******************************************************
// * Additional functions related with handling vectors *
// ******************************************************
  
//! Overloading << for class vector
template<class TYPE>  std::ostream& operator << ( std::ostream & , const StdVector<TYPE> &);

} // end of namespace

// Inclusion of Implementation
#include "StdVector.cc"

#endif
