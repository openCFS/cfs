#ifndef FILE_WORK_2001
#define FILE_WORK_2001
 
namespace CoupledField
{

  //! General class for merging methods of LinSystem and Assemble classes
  /*! This class is derived from class LinSystem and Assemble. It is created only with aim to have for classes LinSystem and Assemble the same base class SystemMatrix
   */ 

template< class Dim, class T_Matrix>
class WorkWithSysMat: virtual public Assemble<Dim,T_Matrix>, virtual public LinSystem<Double, T_Matrix> 
{
public:
  //! constructor (does nothing)
  WorkWithSysMat(Grid<Dim> * aptgrid, const Double eps);

  //! Returns solution
  Vector<Double> getSolution(){ return x;}
};

template<class Dim, class T_Matrix>
inline WorkWithSysMat<Dim, T_Matrix>::WorkWithSysMat(Grid<Dim> * aptgrid, const Double eps)
 : Assemble<Dim, T_Matrix>(aptgrid), LinSystem<Double, T_Matrix>(eps) 
{
#ifdef TRACE
 (*trace) << " entering WorkWithSysMat<Dim>::WorkWithSysMat " << endl;
#endif
 ;
}

} // end of namespace

#endif
