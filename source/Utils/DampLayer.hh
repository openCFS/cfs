#ifndef FILE_DAMPLAYER
#define FILE_DAMPLAYER

#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "FeBasis/BaseFE.hh"
#include "Domain/ElemMapping/EntityLists.hh"

namespace CoupledField
{

  /// Class for calculation of damping factor depending on position

class DampLayer
{
public:

  //! Constructor
  DampLayer(std::string type);

  /// 
  virtual ~DampLayer();

  //! Calculation of stiffmess matrix
  void CalcDampFactor( Double& factor,
		       EntityIterator& ent);

  //! set center and starting points for damping
  void SetDampingParams(Vector<Double> point, 
			Double& dampFactor, 
			Double& dampFactorMax, 
			Double& startRadius, 
			Double& endRadius);


private:

  //! evalate the damping function
  Double EvalDampFnc(Vector<Double>& pos);

  //! type of bilinear form
  std::string dampFncType_;

  //! damping factor
  Double dampingFactor_;

  //! damping factor at saturation
  Double dampingFactorMax_;

  //!center point
  Vector<Double> midPoint_;

  //!start radius
  Double startRadius_;

  //! outer radius
  Double endRadius_;

  //! thickness of layer
  Double layerThickness_;

  //! matrix, holding the coordinates
  Matrix<Double> ptCoord_;

};


}

#endif // FILE_DAMPLAYER
