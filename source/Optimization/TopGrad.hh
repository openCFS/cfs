#ifndef TOPGRAD_HH_
#define TOPGRAD_HH_

#include "Optimization/BaseOptimizer.hh"
#include "Optimization/SIMP.hh"
#include "Utils/vector.hh"

namespace CoupledField
{
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


/** \class TopGrad 
 * \brief Implements topology optimization using the topological gradient.
 * 
 * This is an implementation of an optimization that uses the topological
 * derivative as means of finding the optimal distribution of material under
 * a given volume contstraint for different cost functions. */
class TopGrad : public BaseOptimizer
{
public:
  /** \fn TopGrad(Optimization* optimization, ParamNode* pn)
   * \brief Constructor of class TopGrad.
   * @param optimization the problem we optimize
   * @param pn here we can have options - might be NULL! */
  TopGrad(Optimization* optimization, ParamNode* pn);

  /** \fn ~TopGrad()
   * \brief Destructor of class TopGrad. 
   * Deletes topGrads */
  ~TopGrad();

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
   * - set pseudo-density of these elements to eps > 0
   * - check volume contraint ( if(volume < V) finished; )
   * - goto 1.
   */ 
  void SolveProblem();
    
  /** A list of all topological gradient values on each element; 
   * we remove elements from this list in each iteration
   */
  std::list<TopGradElementValues> topGrads;
  
  /** Here we remember all the elements we removed from topGrads */
  std::list<TopGradElementValues> removedTopGrads;
  
  int getMaxVolumeToRemove() { return max_volume_to_remove_; }

  virtual std::string LogFileHeader();

  virtual void LogFileLine(std::ofstream* out);

private: 
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
  
  /** \var SIMP* simp
   * \brief Reference to our optimization problem */
  SIMP* simp;
};

} // end of namespace

#endif /*TOPGRAD_HH_*/
