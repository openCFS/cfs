#ifndef FILE_ABSMATRIX_2001
#define FILE_ABSMATRIX_2001
 
namespace CoupledField
{      

  //! Abstract class for working with matrixes
  /*!  It is a base class for different types of matrices. Principles of organization
    of storage matrix is similar for all types of matrices. The data of a matrix are allocated as a one contigues memory block and there are pointers to the first element of each row.
  */

template<class T_leaftype, class TYPE>
class AbsMatrix
{
public:
  	
  T_leaftype & asLeaf()
  { return static_cast<T_leaftype&>(*this); }

  //! All elements are zero
  void Init();   

  //! Return size of matrix
  Integer getSize();

  //! Is symmetric or not
  Boolean IsSymmetric();

};

template<class T_leaftype, class TYPE>
inline void AbsMatrix<T_leaftype, TYPE>::Init()
{
 asLeaf.Init() ;
}

template<class T_leaftype, class TYPE>
inline Boolean AbsMatrix<T_leaftype, TYPE>::IsSymmetric()
{
 return FALSE ;
}

template<class T_leaftype, class TYPE>
inline Integer AbsMatrix<T_leaftype, TYPE>::getSize()
{
 return asLeaf.getSize() ;
}

} //end of namespace
#endif	// FILE_MATRIX

