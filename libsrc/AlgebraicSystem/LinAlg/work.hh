#ifndef FILE_WORK_2001
#define FILE_WORK_2001

#include "linsystem.hh"
#include "assemble.hh" 

namespace CoupledField
{

  //! General class for merging methods of LinSystem and Assemble classes
  /*! This class is derived from class LinSystem and Assemble. It is created only with aim to have for classes LinSystem and Assemble the same base class SystemMatrix
   */ 

template<Integer Dim, class T_Matrix>
class WorkWithSysMat: virtual public Assemble<Dim,T_Matrix>, virtual public LinSystem<Double, T_Matrix> 
{
public:
  //! constructor (does nothing)
  WorkWithSysMat(Grid * aptgrid, const Integer level, const Double eps);

  //! Returns solution
  Vector<Double> getSolution(){ return x;}
};

template<Integer Dim, class T_Matrix>
inline WorkWithSysMat<Dim, T_Matrix>::WorkWithSysMat(Grid * aptgrid, const Integer level, const Double eps)
 : Assemble<Dim, T_Matrix>(aptgrid, level), LinSystem<Double, T_Matrix>(eps) 
{
#ifdef TRACE
 (*trace) << " entering WorkWithSysMat<Dim>::WorkWithSysMat " << std::endl;
#endif
 ;
}

} // end of namespace

#endif
