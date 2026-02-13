#ifndef TOPGRAD_HH_
#define TOPGRAD_HH_

#include <functional>
#include <vector>

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/Environment.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Optimization.hh"
#include "boost/shared_ptr.hpp"

namespace CoupledField {
class ShapeGrad;
}  // namespace CoupledField

namespace CoupledField
{
class LevelSet;

/** \class ElementValues
 *  \brief This class links element numbers to gradient values; it is used for shape gradients and topological gradients
 *  @param elem_num the element number
 *  @param gradvalue value of gradient for this element
 * */
class ElementValues
{
public:
  explicit ElementValues(const unsigned int number = 0, const double val = 0.0) : 
    elem_num_(number),
    value_(val)
    { }
  ~ElementValues() {}
  inline double GetValue() const { return value_; }
  inline unsigned int GetElemNum() const { return elem_num_; }
  inline void SetValue(const double val) { value_ = val; }

private:
  unsigned int elem_num_;
  double value_;
};

/** \class CompareElementValues 
 *  \brief Function object for sorting ElementValues comparing their values */
struct CompareElementValues
{
  bool operator()(const ElementValues &lhs, const ElementValues &rhs) const
  {
    return lhs.GetValue() > rhs.GetValue();
  }
};


struct FindElementValuesByNumber
{
  bool operator()(const ElementValues &lhs, const unsigned int elemNum) const
  {
    return lhs.GetElemNum() == elemNum;
  }
};

/** \class TopGrad
 * \brief Implements topology optimization using the topological gradient.
 *
 * This is an implementation of an optimization that uses the topological
 * derivative as means of finding the optimal distribution of material under
 * a given volume constraint for different cost functions. */
class TopGrad
{
public:
  /** \fn TopGrad(Optimization* optimization, PtrParamNode pn)
   * \brief Constructor of class TopGrad.
   * @param optimization the problem we optimize
   * @param pn here we must have options
   * @param ls use topgrad together with levelset */
  explicit TopGrad(Optimization* opt, PtrParamNode pn, const bool ls = false);

  /** \fn ~TopGrad()
   * \brief Destructor of class TopGrad.*/
  ~TopGrad() {}

  /** \fn void SolveProblem()
   * \brief Solves the optimization problem using the topological derivative.
   *
   * Everything (including evaluations of the state problem) is done
   * within this method.
   *
   * The algorithm in short:
   * - solve state and adjoint problem
   * - calculate topological derivative (real number) elementwise
   * - find the lowest values in whole domain
   * - either:
   *      + set pseudo-density of these elements to eps > 0
   * - or:
   *      + create a hole in the levelset
   * - check volume constraint ( if(volume < V) finished; )
   * - goto 1.
   */
  void SolveProblem(const unsigned int iter);
  /** if we use the levelset method, we also need a pointer to the levelset */
  void SolveProblem(const unsigned int iter, boost::shared_ptr<LevelSet> lsptr);
  /** contains the computation that are needed for both of the above methods */
  void SolveProblemCommon(const unsigned int iter);
  
  /** for testing the topology gradient with polarization matrix */
  void SolvePoissonProblem(const unsigned int curr_iter);

  /** With given strain vectors, calculate the topgrad for element nr. e */
  double CalcTopGradOnElement() const;

  double CalcPoissonTopGradOnElement(const unsigned int e) const;
  
  /** helper function; walks over all elements, calcs the topgrad and adds to topGrads vector */
  void CalcTopGrads(SubTensorType sub = PLANE_STRAIN, App::Type app = App::MECH);
  
  inline double getMaxVolumeToRemove() const { return max_volume_to_remove_; }

private:
  /** forbid instance from standard ctor */
  explicit TopGrad();
	/** do not copy the topgrad! */
	TopGrad(const TopGrad &other);
	TopGrad& operator=(const TopGrad &other);
  
  /** \var double max_volume_to_remove_
   * \brief How much volume should be removed in total (from xml) */
  double max_volume_to_remove_;

  /** \var unsigned int elems_removed_per_iteration_
   * \brief The number of elements to be removed in one topGrad iteration (from xml) */
  unsigned int elems_removed_per_iteration_;

  /** \var double pivot_percent_
   * \brief How much volume should be removed in first iteration (from xml).
   *
   * This also determines the number of elements removed in subsequent iterations.
   * After calculating and sorting the vector with the topgradvalues we choose a
   * value according to this percentage in the vector. */
  double pivot_percent_;

  /** \var bool levelset_
   * \brief do we use topgrad in combination with levelset (from xml) */
  bool levelset_;

  /** \var bool adaptive_
   * \brief signals to adaptively change the number of removed elements per iteration (from xml) */
  bool adaptive_;

  /** \var bool more_output_
   * \brief show more output on the console (from xml) */
  bool more_output_;

  /** \var unsigned int elements_
  * \brief number of design elements, taken from designspace */
  unsigned int numelems_;

  /** \var double lambda_
  * \brief Lame material parameter from xml file, read during construction */
  double lambda_;

  /** \var double mu_
   * \brief Lame material parameter from xml file, read during construction */
  double mu_;
  
  double c2;
  double c3;

  /** \var double small_topgrad_value_
   * \brief in the first iteration we determine a topgradvalue so that a given portion of
   * material is removed and use this value in the upcoming iterations as a pivot value */
  double small_topgrad_value_;

  /** \var SubTensorType subtype_
   * \brief determines the type of PDE we are handling */
  SubTensorType subtype_;

  /** A list of all topological gradient values on each element;
  * we remove elements from this list in each iteration */
  std::vector<ElementValues> topGrads;
  
  std::vector<TopGradElement> elements;
  
  /** for convenience */
  Vector<double> elem_forw;
  Vector<double> elem_adj;
  
  Matrix<double> pol;

  /** \var ShapeGrad* optimization
   * \brief Reference to our optimization problem */
  ShapeGrad* optimization;
};

} // end of namespace

#endif /*TOPGRAD_HH_*/
