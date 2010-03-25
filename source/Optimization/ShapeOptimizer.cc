  #include "Optimization/ShapeOptimizer.hh"
#include "Optimization/Optimization.hh"
#include "Domain/domain.hh"
#include "Optimization/TopGrad.hh"
#include "Optimization/LevelSet.hh"
#include "Utils/Timer.hh"

namespace CoupledField
{
using std::cout;
using std::endl;
using std::string;

using boost::posix_time::second_clock;
using boost::posix_time::time_duration;

ShapeOptimizer::ShapeOptimizer(Optimization* optimization, PtrParamNode pn) :
  BaseOptimizer(optimization, pn),
  topgrad_(false),
  levelset_(false),
  shapeopt_(false)
{
  start_time = second_clock::local_time();
  /** ParamNode for reading info from XML-file */
  shoptpn = pn->Get(Optimization::optimizer.ToString(Optimization::SHAPE_SOLVER), ParamNode::PASS);

  if(shoptpn != NULL)
  {
    if(shoptpn->Has("enableTopgrad"))
    {
      topgrad_ = shoptpn->Get("enableTopgrad")->As<bool>();
    }
    if(shoptpn->Has("useLevelSet"))
    {
      levelset_ = shoptpn->Get("useLevelSet")->As<bool>();
    }
    if(shoptpn->Has("useShapeOptimization"))
    {
      shapeopt_ = shoptpn->Get("useShapeOptimization")->As<bool>();
      // shape optimization only works with a levelset
      if(shapeopt_) levelset_ = true;
    }
  }
  else
  {
    throw Exception("No configuration for ShapeOptimizer found!");
  }

  cout << "Use shapeoptimizer with the following properties:" << endl;
  cout << "\t topgrad: " << (topgrad_  ? "yes" : "no") << "\t";
  cout << "\t levelset: " << (levelset_ ? "yes" : "no") << "\t";
  cout << "\t shapeopt: " << (shapeopt_ ? "yes" : "no") << endl;

  PostInit(1.0, true);
}

ShapeOptimizer::~ShapeOptimizer()
{
  cout << "*** shapeoptimizer finished after " << Timer::GetTimeString(second_clock::local_time() - start_time) 
       <<  " ***" << endl;
}

void ShapeOptimizer::SolveProblem()
{
  const int max_iter(optimization->GetMaxIterations());
  int curr_iter(optimization->GetCurrentIteration());

  if(topgrad_)
  {
    // make a new TopGrad object
    boost::shared_ptr<TopGrad> ptrTGtmp(new TopGrad(optimization, shoptpn, levelset_));
    std::swap(ptrTG_, ptrTGtmp);
    assert(ptrTG_ != NULL);
  }

  if(levelset_)
  {
    // build LevelSet (mandatorily needed for shape optimization)
    boost::shared_ptr<LevelSet> ptrLStmp(new LevelSet(optimization, shoptpn));
    std::swap(ptrLS_, ptrLStmp);
    assert(ptrLS_ != NULL);
  }

  optimization->SolveStateProblem();

  while(curr_iter <= max_iter && !optimization->DoStopOptimization())
  {
    // in every iteration we need to solve the state problem again
    if(curr_iter > 0) optimization->SolveStateProblem();
    
    if(dynamic_cast<ErsatzMaterial*>(optimization)->ToPDE(Optimization::MECH, false) != NULL)
    {
      if(topgrad_ && curr_iter > 0)
      {
        if(levelset_)
        {
          // when using topgrad with levelset we also need a pointer to the levelset_
          // so that we can put holes in it
          ptrTG_->SolveProblem(curr_iter, ptrLS_);
        }
        else
        {
          // now we call the TopGrad optimization
          ptrTG_->SolveProblem(curr_iter);
        }
      }
      if(shapeopt_ && curr_iter > 0)
      {
        // now we call the shape optimization
        ptrLS_->SolveProblem(curr_iter);
      }
    }
    else
    {
      // for now this is only the poisson equation
      ptrTG_->SolvePoissonProblem(curr_iter);
    }
    optimization->CalcObjective();
    optimization->CommitIteration(false);
    curr_iter = optimization->GetCurrentIteration();
  }
}

} // end namespace
