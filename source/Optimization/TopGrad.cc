#include "Optimization/TopGrad.hh"
#include "Optimization/LevelSet.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/Condition.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "MatVec/vector.hh"
#include "Utils/myclock.hh"

#include "boost/date_time/posix_time/posix_time.hpp"
#include <string>
#include <functional>

namespace CoupledField
{

DECLARE_LOG(tg)
DEFINE_LOG(tg, "topGrad")

using std::abs;
using std::endl;
using std::cout;
using std::vector;
using std::string;

using boost::posix_time::ptime;
using boost::posix_time::second_clock;
using boost::posix_time::microsec_clock;

TopGrad::TopGrad(Optimization* opt, ParamNode* pn, const bool ls) :
  max_volume_to_remove_(0.5), elems_removed_per_iteration_(2),
  pivot_percent_(0.3), levelset_(ls), adaptive_(false),
  more_output_(false), elements_(0), lambda_(0.0),
  mu_(0.0), pi(std::acos(-1.0)), c2(0.0), c3(0.0), small_topgrad_value_(0.2)
{
  optimization = dynamic_cast<ShapeGrad*>(opt);
  assert(optimization != NULL);

  // get volume constraint from xml
  Condition& condition = optimization->GetConstraint(Condition::VOLUME);
  max_volume_to_remove_ = condition.value;

  // cache number of elements
  elements_ = optimization->GetDesign()->GetNumberOfElements();
  topGrads.reserve(elements_);

  // get Lame material parameters
  optimization->GetMaterialParameters(lambda_, mu_);
  c2 = pi * (lambda_ + 2 * mu_) / (2 * (lambda_ + mu_) * mu_);
  c3 = pi * (lambda_ + 2 * mu_) / ((9 * lambda_ + 14 * mu_) * mu_);

  // transform the subType of the pde (FULL, PLANE_STRAIN etc. -> defined in General/environment.hh)
  optimization->GetSubTensorType(subtype_);

  // reduce to our actual ParamNode
  pn = pn->Get("topGrad", false);
  if(pn != NULL)
  {
    // read necessary parameters from XML
    // note: the lower bound for the pseudo-density is set in the options for the design!
    if(pn->Has("tg_elemsremove")) elems_removed_per_iteration_ = pn->Get("tg_elemsremove")->AsInt();
    if(pn->Has("tg_initpercent")) pivot_percent_ = pn->Get("tg_initpercent")->AsDouble();
    if(pn->Has("tg_enableadaptive")) adaptive_ = pn->Get("tg_enableadaptive")->AsBool();
    if(pn->Has("tg_moreoutput")) more_output_ = pn->Get("tg_moreoutput")->AsBool();
  }
  else
    cout << "No ParamNode found, using default values" << endl;

  if(more_output_)
  {
    cout << "TopGrad-Optimization has the following parameters:" << endl;
    cout << "  max_volume_to_remove_ = " << max_volume_to_remove_
        << ", elemsrem = " << elems_removed_per_iteration_
        << ", percent = " << pivot_percent_
        << ", adaptive = " << (adaptive_ ? "yes" : "no")
        << ", levelset = " << (levelset_ ? "yes" : "no")
        << endl;
  }
}

void TopGrad::SolveProblemCommon(const unsigned int iter)
{
  // will contain the strains on one element
  Vector<double> forward_strain;
  Vector<double> adjoint_strain;

  // first we have to calculate the topgrad values
  // we save them in our own vector of ElementValues
  // in the first iteration we calculate the values on all elements_
  // in higher iterations we only update the remaining elements
  if(iter == 1)
  {
    for(unsigned int e = 0; e < elements_; ++e)
    {
      // call helper to get strains
      optimization->GetStrainsOnElement(forward_strain, adjoint_strain, e, subtype_);
      // do the actual calculation of the topgrad
      const double tg(CalcTopGradOnElement(forward_strain, adjoint_strain));
      // save topgradvalues and element numbers
      topGrads.push_back(ElementValues(e, tg));
    }

    // sort the topgrads
    std::nth_element(topGrads.begin(), topGrads.end() - elems_removed_per_iteration_,
                      topGrads.end(), CompareElementValues());

    if(adaptive_)
    {
      // we only remove a maximum of elems_removed_per_iteration_ of elements in one iteration
      elems_removed_per_iteration_ = elements_ * pivot_percent_;
    }

    // value for comparison
    small_topgrad_value_ = (*(topGrads.end() - elems_removed_per_iteration_)).GetValue();

    if(more_output_)
    {
      cout << "we have " << elements_ << " elements and remove " << elems_removed_per_iteration_ <<
          " of them per iteration" << endl;
      cout << "setting pivot topgradvalue to " << small_topgrad_value_ << endl;
    }
  }
  else //iter > 1
  {
    // if we are not in the first iteration any more, we can speed up the topgrad calculation
    // by only looking at the elements that have not already been removed
    for(vector<ElementValues>::iterator it = topGrads.begin(); it != topGrads.end(); ++it)
    {
      const unsigned int e((*it).GetElemNum());
      optimization->GetStrainsOnElement(forward_strain, adjoint_strain, e, subtype_);
      const double tg(CalcTopGradOnElement(forward_strain, adjoint_strain));
      (*it).SetValue(tg);
    }

  }

  // we are now ready to remove values according to the topgradvalue
  // for sorting we use a function object defined in ElementValues.hh 
  std::nth_element(topGrads.begin(), topGrads.end() - elems_removed_per_iteration_, 
                    topGrads.end(), CompareElementValues());
}

void TopGrad::SolveProblem(const unsigned int iter, boost::shared_ptr<LevelSet> lsptr)
{
  assert(lsptr != NULL);
  ptime before_time = second_clock::local_time();
  
  SolveProblemCommon(iter);

  // walk over all the elements and update the topgradvalue for postprocessing
  for(vector<ElementValues>::iterator it = topGrads.begin(); it != topGrads.end(); ++it)
  {
    // get the corresponding de
    DesignElement* de = &optimization->GetDesign()->data[(*it).GetElemNum()];
    assert(de != NULL); // because we picked an existing number
    assert(de->lse_ != NULL); // because we are in the levelset case
    de->lse_->topGradValue = (*it).GetValue();
  }
  
  vector<ElementValues>::iterator max_it = (topGrads.end() - elems_removed_per_iteration_);
  
  unsigned int c(0);
  //for(unsigned int i = 0; i < elems_removed_per_iteration_; ++i)
  while(max_it != topGrads.end())
  {
    DesignElement* de = &optimization->GetDesign()->data[topGrads.back().GetElemNum()];
    const unsigned int n_corn(de->lse_->nodes_.size());

    for(unsigned int k = 0; k < n_corn; ++k)
    {
      unsigned int nodenr(de->lse_->nodes_[k]->GetNumber());
      LOG_DBG(tg) << "trying to make a hole at node number " << nodenr;
      // now set the hole in the levelset using the lsptr
      lsptr->MakeTrivialHole(nodenr);
    }
    // remove element from list of all elements
    topGrads.pop_back();
    ++c;
  }
  
  cout << "removed " << c << " elements (in " << 
    ShapeOptimizer::GetTimeString(second_clock::local_time() - before_time) << ") - ";
  // reset counter
  before_time = microsec_clock::local_time();
  
  // now, recompute the signed distance function
  lsptr->BuildSignedDistanceFunction();
  // then, update the simp-densities
  lsptr->UpdateDesignForAllLevelSetElements();
  // clock.GetTime(wtime, ctime);

  cout << " (in " << ShapeOptimizer::GetTimeString(microsec_clock::local_time() - before_time) << ")" << endl;
  return;
}

void TopGrad::SolveProblem(const unsigned int iter)
{
  ptime before_time = microsec_clock::local_time();
  
  SolveProblemCommon(iter);

  // now remove material (= set density to lower bound) where topgradvalue is smallest
  // at max we remove elems_removed_per_iteration_ elements
  vector<ElementValues>::iterator max_it = (topGrads.end() - elems_removed_per_iteration_);
  
  unsigned int c(0);
	//if(iter == 150) elems_removed_per_iteration_ = 100;
  //FSFSFS: have to build in the max_vol_to_remove constraint!!
  
  // for(unsigned int i = 0; i < elems_removed_per_iteration_; ++i) // for sort
  while(max_it != topGrads.end()) // for nth_element
  {
    DesignElement* de = &optimization->GetDesign()->data[topGrads.back().GetElemNum()];
    // use the simp densities and set design to lower bound
    de->SetDesign(de->GetLowerBound());
    // remove element from list of all elements
    topGrads.pop_back();
    ++c;
  }
  
  cout << ShapeOptimizer::GetTimeString(microsec_clock::local_time() - before_time) << "; removed " << c << " elements - ";
  // reduce number of removed elements in each iteration
  //elems_removed_per_iteration_ -= 40;
}

double TopGrad::CalcTopGradOnElement(const Vector<double> &forward_strain, const Vector<double> &adjoint_strain) const
{
  double topG(0.0);
  double tr_eu2(0.0);
  double mu2(40.0 * mu_ * mu_);

  switch(subtype_)
  {
    case PLANE_STRAIN: // top-grad 2d *neumann* allaire!
      /*
       * formular is
       * c_2 = \pi * (lambda_ + 2mu_)/(2mu_(lambda_ + mu_))
       * tg = c_2(4mu_ * Ae(u)*e(u) + (lambda_ - mu_)tr(Ae(u))tr(e(u)))
       * entries of e(u) are contained in forward_strain[0 - 2]      */
      assert(forward_strain.GetSize() == 3);
      assert(adjoint_strain.GetSize() == 3);
   
      tr_eu2 = forward_strain[0] + forward_strain[1];
      tr_eu2 *= adjoint_strain[0] + adjoint_strain[1];
      tr_eu2 *= lambda_ - mu_;
      tr_eu2 *= 2.0*lambda_ + 2.0*mu_;

      topG = (lambda_ + 2.0*mu_) * (forward_strain[0]*adjoint_strain[0] + forward_strain[1]*adjoint_strain[1]);
      topG += lambda_ * (forward_strain[1]*adjoint_strain[0] + forward_strain[0]*adjoint_strain[1]);
      topG += 4.0 * mu_ * (forward_strain[2]*adjoint_strain[2]);
      topG *= 4.0 * mu_;
      
      topG += tr_eu2;
      
      topG *= c2;
      break;
     

      /*
      tr_eu2 = forward_strain[0] + forward_strain[1];
      tr_eu2 *= tr_eu2;

      topG  =  lambda_ * lambda_ * tr_eu2;
      topG +=  5 * mu_ * lambda_ * tr_eu2;
      topG -=  2 * mu_ * mu_ * tr_eu2;

      topG +=  mu1 * (forward_strain[0] * forward_strain[0]);
      topG +=  mu1 * (forward_strain[1] * forward_strain[1]);
      topG += 2 * mu1 * (forward_strain[2] * forward_strain[2]);

      topG *= c2;
      break;
        */    

    case FULL: // top-grad 3d *neumann*
      assert(forward_strain.GetSize() == 6);
      tr_eu2 = forward_strain[0] + forward_strain[1] + forward_strain[2];
      tr_eu2 *= tr_eu2;

      topG =  20.0 * mu_ * lambda_;
      topG += 9.0 * lambda_ * lambda_;
      topG -= 4.0 * mu_ * mu_;
      topG *= tr_eu2;

      topG += mu2 * (forward_strain[0] * forward_strain[0]);
      topG += mu2 * (forward_strain[1] * forward_strain[1]);
      topG += mu2 * (forward_strain[2] * forward_strain[2]);

      topG += 2.0 * mu2 * (forward_strain[3] * forward_strain[3]);
      topG += 2.0 * mu2 * (forward_strain[4] * forward_strain[4]);
      topG += 2.0 * mu2 * (forward_strain[5] * forward_strain[5]);

      topG *= c3;
      break;

    case AXI:
      EXCEPTION("SubType AXI not implemented!");
    case PLANE_STRESS:
      EXCEPTION("SubType PLANE_STRESS not implemented!");
    default:
      EXCEPTION("SubType not implemented!");
  }
  return topG;
}

} // end of namespace
