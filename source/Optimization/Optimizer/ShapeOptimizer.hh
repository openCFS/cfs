#ifndef SHAPEOPTIMIZER_HH_
#define SHAPEOPTIMIZER_HH_

#include <string>
#include "Optimization/Optimizer/BaseOptimizer.hh"
#include "boost/date_time/posix_time/posix_time.hpp"
#include "boost/shared_ptr.hpp"

namespace CoupledField
{
class LevelSet;
class Optimization;
class TopGrad;

/** Optimizer for shape- and topology optimization! */
class ShapeOptimizer : public BaseOptimizer
{
public:
  explicit ShapeOptimizer(Optimization* optimization, PtrParamNode pn);
  virtual ~ShapeOptimizer();

  
  /** pointer to topgrad */
  boost::shared_ptr<TopGrad> ptrTG_;
  /** pointer to the level set */
  boost::shared_ptr<LevelSet> ptrLS_;

protected:

  /** Implements virtual function from BaseOptimizer */
  void SolveProblem();
  
private:
  ShapeOptimizer(); // forbid empty standard constructor
  
  /** for timing the shape optimizer */
  boost::posix_time::ptime start_time;
  
  PtrParamNode shoptpn;
  
  /** evaluate TopGrad? Read from XML */
  bool topgrad_;
  /** use level set? mandatory for shape optimization! Read from XML */
  bool levelset_;
  /** do shape optimization? Read from XML */
  bool shapeopt_;
};


} // end of namespace

#endif /*SHAPEOPTIMIZER_HH_*/
