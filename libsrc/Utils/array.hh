#ifndef _FILE_Array_2003
#define _FILE_Array_2003

#include <vector>
#include <General/environment.hh>
#include <Utils/vector.hh>
#include <iostream>

namespace CoupledField
{

  template<class TYPE> class Matrix;

  template<class TYPE> class Array;
  
  template<class TYPE> class Vector;

  template<class TYPE> std::ostream &operator << (std::ostream &, const Array<TYPE> &);

  //! This class stores values in the format row(dimension) x column(values)
template<class TYPE> class Array{

public:
  
  friend class Matrix<TYPE>;
  friend class Vector<TYPE>;

  //! constructor
  Array();

  //! constructor with shape
  Array(ShortInt dim, Integer size);

  //! copy constructor
  Array(const Array& X);

  //! constructor with std-vector
  Array(std::vector<TYPE> & v);

  //! constructor with Vector
  Array(Vector<TYPE> & v);

  //! destructor
  ~Array();

  //! erases all values
  void clear();

  //! initalizes array with given value
  void init();

  //! changes dimension of array
  void redim(Integer dim);

  //! changes size of array
  void resize(Integer size);

  //! changes dimension and size of array
  void reshape(ShortInt Dim, Integer size);

  //! sets row-values of array (= values belonging to one point)
  void setValuesRow(Vector<TYPE> & v, Integer pos);

  //! set column-values of array (= values belonging to one dimension)
  void setValuesColumn(Vector<TYPE> & v, Integer pos);

  //! adds an array at the end
  void push_back(Array & arr);

  //! adds an vector at the and
  void push_back(Vector<TYPE> & v);

  //! returns the size
  Integer size() const {return size_;}

  //! return the dimension
  Integer dim() const {return dim_;}

  //! access operator
  inline Vector<TYPE>& operator[] (Integer dim) const 
  {
    if (dim >= dim_)
      Error("Array: index out of bounds",__FILE__,__LINE__);

    return sol_[dim];
  }

  

//   //!
//   inline Vector<TYPE> operator[] (Integer dim) const
//   {
//     if (dim >= dim_)
//       Error("Array: index out of bounds",__FILE__,__LINE__);
    
//     return sol_[dim];
//   }

  //! assignment operator with array
  Array& operator= (const Array & x);

  //! assignment operator with std-vector
  Array& operator= (const std::vector<TYPE> x);

  //! assignent operator with Vector
  Array& operator= (const Vector<TYPE> x);


  // Overloading of +, +=
  //!
  Array  operator+ () const;
  //! 
  Array operator+ (const Array & x) const ;
  //!
  Array& operator+= (const Array & x);


  // Overloading of -, -=
  //!
  Array operator- () const;
  //! 
  Array operator- (const Array & x) const;
  //! 
  Array& operator-= (const Array & x);


  //! scalar product
  TYPE operator* (const Array & x) const;


  // Overloading of *, *=
  //! 
  Array operator* (TYPE x) const;
  //!
  Array& operator*= (TYPE x);
  

  // Overloading of /, /=
  //!  
  Array operator/ (TYPE x) const;

  //! 
  Array& operator/= (TYPE x);  
  
  //! Calculation of L2 Norm at position i
  Double normL2 (Integer pos);

  //! Calcultion of summed L2 Norm
  Double normL2 ();


  // -------- conversion functions ------------

  //! converts array to matrix
  void toMatrix(Matrix<TYPE> & m);

  //! converts array to std-vector
  void toStdVector(std::vector<TYPE> & v, Integer dim = 0);

  //! converts array to Vector
  void toVector(Vector<TYPE> & v, Integer dim = 0);

 
private:

  ShortInt dim_;                          //!< dimension of array

  Integer size_;                          //!< size of array

  Vector<TYPE> * sol_;                    //!< values


};

#ifdef __GNUC__
template class Array<Double>;
template class Array<Integer>;
#endif

} // end of namespace

#endif
