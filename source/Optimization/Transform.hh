#ifndef OPTIMIZATION_TRANSFORM_HH_
#define OPTIMIZATION_TRANSFORM_HH_

#include <cstddef>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Optimization/Design/DesignElement.hh"
#include "General/Enum.hh"
#include "Utils/Point.hh"

namespace CoupledField
{

class DesignSpace;

/** A transformation builds a physical value after filtering and application of the TransferFunction.
 * A transformation has nothing to to with the TransferFunction which is just a poor name for an interpolation/ penalization of the design.
 *
 *  <transform>
 *     <rotate center="center" angle="0.5">
 *   </transform>
 * */

class Transform
{
public:
  /** for StdVector only */
  Transform() { space_ = NULL; type_ = ROTATION; index = -1; };

  /** reads a transform child, e.g. rotate but not the list transform! */
  Transform(PtrParamNode pn, DesignSpace* space);

  /** finds the source of the transformation of element de by doing the inverse transformation
   * It is clearly faster to generate a standardized regular mesh in 2D and 3D!
   * @param tf to get the physical design */
  DesignElement* FindSource(const DesignElement* de) const;

  void ToInfo(PtrParamNode in);

  /** @param level: -1 is debug, 0 is very brief (only angle), 1 adds also the index */
  std::string ToString(int level = 0);

  typedef enum { ROTATION } Type;

  static Enum<Type> type;

  /** at time of constructor call there are no multiple excitations yet */
  std::string excitation_str;

  /** this is the index of the transformation, shall coincide with the excitation_str_ representation */
  int index;

private:

  /** helper for FindSource. Traverses the regular design space in one axis direction via vicinity elements
   * @param dir e.g. X_P or X_N
   * @param val for dir is X_P or X_N the component 0 of the target element
   * @return NULL if last element within space was too far (using h_)*/
  DesignElement* SearchDesignSpace(const DesignElement* start, VicinityElement::Neighbour dir, double val) const;

  /** counterclockwise in rad */
  double angle_;


  std::string center_str_;

  Point center_;

  Type type_;

  DesignSpace* space_;

  /** the largest edge size for the center element */
  double h_;
};


} // end of namespace

#endif /* OPTIMIZATION_TRANSFORM_HH_ */
