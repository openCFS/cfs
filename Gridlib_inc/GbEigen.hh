/*----------------------------------------------------------------------
|
| $Id$
|
+---------------------------------------------------------------------*/

#ifndef  GBEIGEN_HH
#define  GBEIGEN_HH

#include "GbDefines.hh"

/*----------------------------------------------------------------------
|       includes
+---------------------------------------------------------------------*/

#include <typeinfo>

#include "GbTypes.hh"
#include "GbMath.hh"

/*----------------------------------------------------------------------
|       declaration
+---------------------------------------------------------------------*/

template <class T>
class GbEigen
{
public:
  GbEigen(int iSize);
  ~GbEigen();

  //! Set the matrix for eigensolving
  INLINE T& matrix (int iRow, int iCol);
  INLINE void setMatrix (T** aafMat);

  //! Get the results of eigensolving (eigenvectors are columns of matrix)
  INLINE int getSize() const;
  INLINE T getEigenvalue (int i) const;
  INLINE T getEigenvector (int iRow, int iCol) const;
  INLINE T* getEigenvalue ();
  INLINE T** getEigenvector ();

  //! Solve eigensystem
  INLINE void eigenSolv();

  //! Solve eigensystem, use decreasing sort on eigenvalues
  INLINE void decrSortEigenSolv();

  //! Solve eigensystem, use increasing sort on eigenvalues
  INLINE void incrSortEigenSolv();

  //! This operator displays the values
  friend std::ostream& operator<<(std::ostream&, const GbEigen<T>&);

private:
  int size_;
  T** matrix_;
  T* diagonal_;
  T* subdiagonal_;

  //! Householder reduction to tridiagonal form
  void tridiagonal2 ();
  void tridiagonal3 ();
  void tridiagonal4 ();
  void tridiagonalN ();

  //! QL algorithm with implicit shifting, applies to tridiagonal matrices
  GbBool QLAlgorithm ();

  //! sort eigenvalues from largest to smallest
  void decreasingSort ();

  //! sort eigenvalues from smallest to largest
  void increasingSort ();
};

template<class T>
std::ostream&
operator<<(std::ostream& s, const GbEigen<T>& v)
{
  s<<"eigenvalues  ="<<std::endl<<"\t";
  for (int i=0; i<v.getSize(); ++i)
    s<<v.getEigenvalue(i)<<" ";
  s<<std::endl<<"eigenvectors ="<<std::endl;
  for (int i=0; i<v.getSize(); ++i) {
    s<<"\t";
    for (int j=0; j<v.getSize(); ++j)
      s<<v.getEigenvector(i,j)<<" ";
    s<<std::endl;
  }
  return s;
}


#ifdef __GNUG__

#include "GbEigen.in"
#include "GbEigen.T"

#else

#pragma instantiate GbEigen<float>
#pragma instantiate GbEigen<double>
#pragma instantiate std::ostream& operator<<(std::ostream&, const GbEigen<float>&)
#pragma instantiate std::ostream& operator<<(std::ostream&, const GbEigen<double>&)

#ifndef OUTLINE
#include "GbEigen.in"
#endif

#endif  // g++

#endif  // GBEIGEN_HH
/*----------------------------------------------------------------------
|
| $Log$
| Revision 1.1  2002/02/22 14:47:56  elena
| new: dir Gridlib_inc
|
| Revision 1.1  2001/01/02 15:02:11  prkipfer
| introduced new classes
|
|
+---------------------------------------------------------------------*/
