#ifndef FILE_DAMPINT
#define FILE_DAMPINT

#include "baseform.hh" 

namespace CoupledField
{
//! class for calculation of boundary damping matrix
template<Integer dim>
class DampInt : public BaseForm<dim>
{
public:
  //! Constructor with pointer to BaseElem
	/*!
	\param aptelem pointer to Elems
	\param ndofs number of dofs
	*/
  DampInt(BaseFE * aptelem, const ShortInt ndofs);

  //! Deconstructor
  virtual ~DampInt();

  //! Function for calculation damping matrix \f $\int_surface N_i N_j d\Gamma$ 
  /*!
	\param ptCoord pointer to coordinates of element's nodes
	\param elemmat (output) element boundary damping matrix
*/
  virtual void CalcElemMatrix(Point<dim> * ptCoord, Matrix<Double> & elemmat);

protected:

private:
  //! dofs per node
  ShortInt DofsPerNode;
};

#ifdef __GNUC__
template class DampInt<2>;
template class DampInt<3>;
#endif


}

#endif // FILE_DAMPINT
