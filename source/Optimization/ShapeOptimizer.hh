#ifndef SHAPEOPTIMIZER_HH_
#define SHAPEOPTIMIZER_HH_

#include "Optimization/BaseOptimizer.hh"

#include "boost/date_time/posix_time/posix_time.hpp"
#include <boost/shared_ptr.hpp>
#include <string>

namespace CoupledField
{
class Optimization;
class ParamNode;
class TopGrad;
class LevelSet;

/** Optimizer for shape- and topology optimization! */
class ShapeOptimizer : public BaseOptimizer
{
public:
  explicit ShapeOptimizer(Optimization* optimization, ParamNode* pn);
  virtual ~ShapeOptimizer();

  /** Implements virtual function from BaseOptimizer */
  void SolveProblem();
  
  /** pointer to topgrad */
  boost::shared_ptr<TopGrad> ptrTG_;
  /** pointer to the level set */
  boost::shared_ptr<LevelSet> ptrLS_;
  
  /** static function for use in topgrad and levelset classes; converts a time interval
   *  to a formatted string suitable for command line output */
  static const std::string GetTimeString(const boost::posix_time::time_duration period);
  
private:
  ShapeOptimizer(); // forbid empty standard constructor
  
  /** for timing the shape optimizer */
  boost::posix_time::ptime start_time;
  
  ParamNode* shoptpn;
  
  /** evaluate TopGrad? Read from XML */
  bool topgrad_;
  /** use level set? mandatory for shape optimization! Read from XML */
  bool levelset_;
  /** do shape optimization? Read from XML */
  bool shapeopt_;
};


} // end of namespace

#endif /*SHAPEOPTIMIZER_HH_*/
