#ifndef SHAPEGRAD_HH_
#define SHAPEGRAD_HH_

#include "Optimization/ErsatzMaterial.hh"

namespace CoupledField
{
class MechPDE;

/** \class TopGradElementValues
 * \brief This class links element numbers to topGradValues
 * @param elem_num the element number
 * @param topgradvalue value of topological gradient for this element
 */
class TopGradElementValues {
public:
  TopGradElementValues() { elem_num = 0; topgradvalue = 0.0; }
  TopGradElementValues(unsigned int number, double value) { elem_num = number; topgradvalue = value; }
  ~TopGradElementValues() {}
  void UpdateTopGrad(double value) { topgradvalue = value; }
  bool operator<(TopGradElementValues& a) { return (topgradvalue < a.topgradvalue); }
  bool operator>(TopGradElementValues& a) { return (topgradvalue > a.topgradvalue); }

  unsigned int elem_num;
  double topgradvalue;
};



/** Optimization via ShapeGradient and Level-Set method, not by Parametrization */
class ShapeGrad : public ErsatzMaterial
{
public:
  /** Up to now w/o parameters */
  ShapeGrad();
  
  virtual ~ShapeGrad()
  {
  }


  /** ... */
  void CalcObjectiveGradient(double* grad_out);
  
  double CalcTopGrad(unsigned int iter, unsigned int elems_removed_per_iteration, 
      std::list<TopGradElementValues>& tgvalues, std::list<TopGradElementValues>& remtgvalues);
  
  // Bastian Schmidt: MultiplyForwardAndAdjointSolution is always called with mechStiffness
  double MultiplyForwardAndAdjointSolution(Solution* sol_forward,
      Solution* sol_adjoint, unsigned int iter,
      unsigned int elems_removed_per_iteration);

  /** A list of all topological gradient values on each element;
   * we remove elements from this list in each iteration
   */
  std::list<TopGradElementValues> topGrads;

  /** Here we remember all the elements we removed from topGrads */
  std::list<TopGradElementValues> removedTopGrads;

  int getMaxVolumeToRemove() { return max_volume_to_remove_; }

  virtual std::string LogFileHeader();

  virtual void LogFileLine(std::ofstream* out);
  
  Matrix<double> mechStiffness;
  
private:
  
  /** stores in DesignElement */
  void CalcTopGrad();
  
  /** evaluate */
  bool topgrad;
  
  /** \var int max_volume_to_remove_
   * \brief How much volume should be removed in total (from xml) */
  int max_volume_to_remove_;

  /** \var int elems_removed_per_iteration_
   * \brief The number of elements to be removed in one topGrad iteration (from xml) */
  int elems_removed_per_iteration_;

  /** \var int remove_volume_fraction_
   * \brief How much volume should be removed in one iteration (from xml) */
  int remove_volume_fraction_;

  /** \var bool volume_not_elems_
   * \brief signals to use remove_volume_fraction instead of directly specifying the number of
   * elements to be removed in one iteration */
  bool volume_not_elems_;

  /** \var bool enable_damping_
   * \brief signals to use a damping factor to reduce the number of removed elements
   * the more iterations we do (from xml)
   */
  bool enable_damping_;
};


} // namespace


#endif /*SHAPEGRAD_HH_*/
