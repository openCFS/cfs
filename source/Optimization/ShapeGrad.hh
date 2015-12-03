#ifndef SHAPEGRAD_HH_
#define SHAPEGRAD_HH_

#include <iostream>
#include <map>
#include <string>

#include "General/Environment.hh"
#include "Optimization/ErsatzMaterial.hh"
#include "Optimization/Optimization.hh"
#include "Utils/StdVector.hh"

namespace CoupledField {
class SingleVector;
template <class TYPE> class Vector;
}  // namespace CoupledField

namespace CoupledField
{
class MechMat;
class linElastInt;

/** Optimization via ShapeGradient and Level-Set method, not by Parametrization */
class ShapeGrad : public ErsatzMaterial
{
public:
  /** Up to now w/o parameters */
  explicit ShapeGrad();
  virtual ~ShapeGrad() {}

  /** Calculates the Lame material parameters from poisson ratio and elasticity module
  * defined in xml */
  void GetMaterialParameters(double &lambda, double &mu) const;

  /** Calculates the strains for the forward and adjoint solution on every element
  @param forward for the strains of the forward solution
  @param adjoint for the strains of the adjoint solution
  */
  void GetElementSolution(Vector<double> &vecforward, Vector<double> &vecadjoint, 
                          const unsigned int e,
                          const SubTensorType type = PLANE_STRAIN,
                          App::Type app = App::MECH);

  StdVector<SingleVector*>& getSolutionVectors(const bool forward_solution = true)
  {
    if(forward_solution)
      return forward.Get(0)->elem[App::MECH];
    else // adjoint
      return adjoint.Get(0)->elem[App::MECH];
  }

  /** called in LevelSet::CalcShapeGradientOnAllElements() */
  linElastInt* getBDBForm();

  int getMaxVolumeToRemove() const { return max_volume_to_remove_; }
  
  void PrepareExteriorPiezoProblem() { std::cout << "ShapeGrad, prepare the problem!!" << std::endl; }

  //virtual std::string LogFileHeader() { return ""; }
  virtual void LogFileLine(std::ofstream* out) {}

protected:

  // kind of second phase constructor
  void PostInit();

private:
  /** \var int max_volume_to_remove_
   * \brief How much volume should be removed in total (from xml) */
  int max_volume_to_remove_;
};

} // namespace


#endif /*SHAPEGRAD_HH_*/
