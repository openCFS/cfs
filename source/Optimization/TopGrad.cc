#include "Optimization/TopGrad.hh"
#include "Optimization/SIMP.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/DesignElement.hh"
#include "General/exception.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"

#include <string>

using namespace CoupledField;
//using std::abs;

// DECLARE_LOG(tg)
// DEFINE_LOG(tg, "topGrad")

TopGrad::TopGrad(Optimization* optimization, ParamNode* pn) : BaseOptimizer(optimization, pn)
{
  // die volumenbeschraenkung sollte man vielleicht besser von GetConstraint holen
  // Condition& condition = optimization->GetConstraint(Condition::VOLUME);
  
  std::cout << "...using topology gradient optimization\n";
  simp = dynamic_cast<SIMP*>(optimization);

  // this defines how many elements should be removed per topGrad-Iteration
  // value is overwritten by definition in XML-file
  this->elems_removed_per_iteration_ = 1; 
  // alternatively we remove a fixed fraction of the volume in each iteration
  this->remove_volume_fraction_ = 10;
  this->volume_not_elems_ = true;
  // do not enable damping by default
  this->enable_damping_ = false;
  // by default we remove 50 percent of the volume
  this->max_volume_to_remove_ = 50;

  // reduce to our actual ParamNode
  pn = pn->Get(Optimization::optimizer.ToString(Optimization::TOPOLOGY_GRADIENT), false);
  if(pn != NULL)
  {
    // read necessary parameters from XML
    // note: the lower bound for the pseudo-density is set 
    // in the options for the design!
    this->max_volume_to_remove_ = pn->Get("max_vol")->AsInt();
    this->elems_removed_per_iteration_ = pn->Get("remove_elems_perit")->AsInt();
    this->remove_volume_fraction_ = pn->Get("remove_volfrac")->AsInt();
    this->volume_not_elems_ = pn->Get("vol_not_elems")->AsBool();
    this->enable_damping_ = pn->Get("use_damping")->AsBool();
  }
  
  if(this->volume_not_elems_)
  {
    std::cout << optimization->GetDesign()->GetNumberOfElements() << " elements; removing " << this->remove_volume_fraction_ << "% of volume per iteration";
    this->elems_removed_per_iteration_ = this->remove_volume_fraction_ * optimization->GetDesign()->GetNumberOfElements() / 100;
    std::cout << " = " << this->elems_removed_per_iteration_ << " elements";
  }
  else std::cout << "Removing " << this->elems_removed_per_iteration_ << " elements per iteration";
  
  std::cout << (this->enable_damping_ ? " with damping" : "");
  std::cout << std::endl;
}

TopGrad::~TopGrad()
{
  while(topGrads.empty() == false) topGrads.pop_back();
  while(removedTopGrads.empty() == false) removedTopGrads.pop_back();
}

void TopGrad::SolveProblem()
{
  unsigned int iter = 1;
  unsigned int max_iter = optimization->GetMaxIterations();
  double vol = 0.0, volprev = 0.0;
  int numberofelems = 0, numberofelemsfirstit = 0;
  
  while(iter <= max_iter && (int) vol < this->max_volume_to_remove_)
  {
    /* In each iteration for the topology optimization
     * calculate the topology gradient, save using DesignElement::SetTopGradient
     * and set pseudo-density of elements with lowest value to lower bound
     */
    volprev = vol;
    vol = optimization->CalcTopGrad(iter, elems_removed_per_iteration_, topGrads, removedTopGrads);

    // CommitIteration produces the ouput for gid/gmv and the logfiles
    optimization->CommitIteration();

    if(this->enable_damping_)
    {
      // in each iteration we reduce the number of removed elements by a certain factor
      if(iter == 1) 
      {
        numberofelemsfirstit = this->elems_removed_per_iteration_;
        numberofelems = this->elems_removed_per_iteration_;
      }

      if(iter ==  3) numberofelems = 2*numberofelemsfirstit/3;
      if(iter ==  6) numberofelems = numberofelemsfirstit/2;
      if(iter ==  9) numberofelems = numberofelemsfirstit/3;
      if(iter == 18) numberofelems = numberofelemsfirstit/4;

      this->elems_removed_per_iteration_ = numberofelems;
    }

    // if we do not change anything we might as well break
    if(vol - volprev == 0.0) break;
    iter++;
  }
  
  if((int) vol >= this->max_volume_to_remove_) std::cout << "topGrad: max volume removed (" << max_volume_to_remove_ << "%)" << std::endl;
  else if(iter >= max_iter) std::cout << "topGrad: max iterations reached (" << max_iter << ")" << std::endl;
  else if(vol - volprev == 0.0) std::cout << "topGrad: no more elements were removed" << std::endl;
}

