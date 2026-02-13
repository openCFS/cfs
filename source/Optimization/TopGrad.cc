#include <assert.h>
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>
#include <string>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/Mesh/Grid.hh"
#include "FeBasis/BaseFE.hh"
#include "General/defs.hh"
#include "General/Exception.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Design/DesignElement.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Function.hh"
#include "Optimization/LevelSet.hh"
#include "Optimization/ShapeGrad.hh"
#include "Optimization/TopGrad.hh"
#include "PDE/SinglePDE.hh"
#include "Utils/StdVector.hh"
#include "Utils/Timer.hh"
#include "Utils/tools.hh"
#include <chrono>

namespace CoupledField
{

DEFINE_LOG(tg, "topGrad")

using std::abs;
using std::endl;
using std::cout;
using std::vector;
using std::string;


TopGrad::TopGrad(Optimization* opt, PtrParamNode pn, const bool ls) :
  max_volume_to_remove_(0.5),
  elems_removed_per_iteration_(2),
  pivot_percent_(0.3),
  levelset_(ls),
  adaptive_(false),
  more_output_(false),
  numelems_(0),
  lambda_(0.0),
  mu_(0.0),
  c2(0.0),
  c3(0.0),
  small_topgrad_value_(0.2),
  optimization(dynamic_cast<ShapeGrad*>(opt))
{
  assert(optimization != NULL);

  // get volume constraint from xml
  Condition* condition = optimization->constraints.Get(Condition::VOLUME);
  max_volume_to_remove_ = condition->GetBoundValue();

  // cache number of elements
  numelems_ = optimization->GetDesign()->GetNumberOfElements();
  topGrads.reserve(numelems_);
  elements.reserve(numelems_);
  
  // walk over all design elements and attach topgrad elements
  for(unsigned int e = 0; e < numelems_; ++e)
  {
    elements.push_back(TopGradElement(1.0));
    // link the new element to the design element
    (optimization->GetDesign()->data[e]).tge = &elements.back();
  }
  
  // get Lame material parameters
  optimization->GetMaterialParameters(lambda_, mu_);
  c2 = M_PI * (lambda_ + 2.0 * mu_) / (2.0 * (lambda_ + mu_) * mu_);
  c3 = M_PI * (lambda_ + 2.0 * mu_) / ((9.0 * lambda_ + 14.0 * mu_) * mu_);

  // transform the subType of the pde (FULL, PLANE_STRAIN etc. -> defined in General/Environment.hh)
  subtype_ = optimization->context->stt;
  
  // FIXME resize member vectors and matrix according to subtype

  // reduce to our actual ParamNode
  pn = pn->Get("topGrad", ParamNode::PASS);
  if(pn != NULL)
  {
    // read necessary parameters from XML
    // note: the lower bound for the pseudo-density is set in the options for the design!
    if(pn->Has("tg_elemsremove")) elems_removed_per_iteration_ = pn->Get("tg_elemsremove")->As<Integer>();
    if(pn->Has("tg_initpercent")) pivot_percent_ = pn->Get("tg_initpercent")->As<Double>();
    if(pn->Has("tg_enableadaptive")) adaptive_ = pn->Get("tg_enableadaptive")->As<bool>();
    if(pn->Has("tg_moreoutput")) more_output_ = pn->Get("tg_moreoutput")->As<bool>();
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
  // first we have to calculate the topgrad values
  // we save them in our own vector of ElementValues
  // in the first iteration we calculate the values on all numelems_
  // in higher iterations we only update the remaining elements
  if(iter == 1)
  {
    CalcTopGrads(subtype_, App::MECH);

    // sort the topgrads
    std::sort(topGrads.begin(), topGrads.end(), CompareElementValues());
    // std::nth_element(topGrads.begin(), topGrads.end() - elems_removed_per_iteration_,
    //                  topGrads.end(), CompareElementValues());

    if(adaptive_)
    {
      // we only remove a maximum of elems_removed_per_iteration_ of elements in one iteration
      elems_removed_per_iteration_ = static_cast<unsigned int>(static_cast<double>(numelems_) * pivot_percent_);
    }

    // value for comparison
    // small_topgrad_value_ = (*(topGrads.end() - elems_removed_per_iteration_)).GetValue();
    small_topgrad_value_ = (topGrads[topGrads.size() - elems_removed_per_iteration_]).GetValue();

    if(more_output_)
    {
      cout << "we have " << numelems_ << " elements and remove " << elems_removed_per_iteration_ <<
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
      optimization->GetElementSolution(elem_forw, elem_adj, e, subtype_);
      const double tg(CalcTopGradOnElement());
      (*it).SetValue(tg);
    }

  }

  // we are now ready to remove values according to the topgradvalue
  // for sorting we use a function object defined in ElementValues.hh 
  // std::nth_element(topGrads.begin(), topGrads.end() - elems_removed_per_iteration_, 
  //                  topGrads.end(), CompareElementValues());
  std::sort(topGrads.begin(), topGrads.end(), [] (const ElementValues &lhs, const ElementValues &rhs) {return lhs.GetValue() > rhs.GetValue();});

}

void TopGrad::SolveProblem(const unsigned int iter, boost::shared_ptr<LevelSet> lsptr)
{
  assert(lsptr != NULL);
  auto before_time = std::chrono::steady_clock::now();
  
  SolveProblemCommon(iter);

  // walk over all the elements and update the topgradvalue for postprocessing
  for(vector<ElementValues>::iterator it = topGrads.begin(); it != topGrads.end(); ++it)
  {
    // get the corresponding de
    DesignElement* de = &optimization->GetDesign()->data[(*it).GetElemNum()];
    assert(de != NULL); // because we picked an existing number
    assert(de->tge != NULL); // because we are calculating the topgrad
    assert(de->lse_ != NULL); // because we are in the levelset case
    de->tge->value = (*it).GetValue();
  }
  
  // vector<ElementValues>::iterator max_it = (topGrads.end() - elems_removed_per_iteration_);
  
  unsigned int c(0);
  
  // while(max_it != topGrads.end())
  for( ; c < elems_removed_per_iteration_; ++c)
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
  }
  
  cout << "removed " << c << " elements (in " << 
    Timer::Duration(std::chrono::steady_clock::now() - before_time) << ") - ";
  // reset counter
  before_time = std::chrono::steady_clock::now();
  
  // now, recompute the signed distance function
  lsptr->BuildSignedDistanceFunction();
  // then, update the simp-densities
  lsptr->UpdateDesignForAllLevelSetElements();
  // clock.GetTime(wtime, ctime);

  cout << " (in " << Timer::Duration(std::chrono::steady_clock::now() - before_time) << ")" << endl;
  return;
}

void TopGrad::SolveProblem(const unsigned int iter)
{  
  auto before_time = std::chrono::steady_clock::now();
  
  SolveProblemCommon(iter);

  // now remove material (= set density to lower bound) where topgradvalue is smallest
  // at max we remove elems_removed_per_iteration_ elements
  // vector<ElementValues>::iterator max_it = (topGrads.end() - elems_removed_per_iteration_);
  
  unsigned int c(0);
	//if(iter == 150) elems_removed_per_iteration_ = 100;
  //FSFSFS: have to build in the max_vol_to_remove constraint!!
  
  // while(max_it != topGrads.end()) // for nth_element
  for( ; c < elems_removed_per_iteration_; ++c) // for sort
  {
    DesignElement* de = &optimization->GetDesign()->data[topGrads.back().GetElemNum()];
    // use the simp densities and set design to lower bound
    de->SetDesign(de->GetLowerBound());
    // remove element from list of all elements
    topGrads.pop_back();
  }
  
  cout << Timer::Duration(std::chrono::steady_clock::now() - before_time) << "; removed " << c << " elements - ";
  // reduce number of removed elements in each iteration
  //elems_removed_per_iteration_ -= 40;
}

void TopGrad::SolvePoissonProblem(const unsigned int curr_iter)
{
  Matrix<double> A1(2, 2); // material 1 = identity
  A1[0][0] = 1.0;
  A1[1][1] = 1.0;
  Matrix<double> A0(2, 2); // material 0
  A0[0][0] = 0.01;
  A0[1][1] = 100.0;
  
  Matrix<double> D(A1);
  D -= A0;
  
  Matrix<double> tmp(A1);
  tmp -= D;
  Matrix<double> inv;
  tmp.Invert(inv);
  
  tmp = D * inv;
  tmp *= D;
  
  pol.Resize(2, 2); // polarization matrix
  pol = tmp;
  pol += D;
  pol *= -1;
  
  cout << "Pol = " << pol.ToString() << endl;
  
  CalcTopGrads(PLANE_STRAIN, App::HEAT);
}

double TopGrad::CalcTopGradOnElement() const
{
  double topG(0.0);
  double tr_eu2(0.0);
  const double mu2(40.0 * mu_ * mu_);

  switch(subtype_)
  {
    case PLANE_STRAIN: // top-grad 2d *neumann* allaire!
      /*
       * formula is
       * c_2 = \pi * (lambda_ + 2mu_)/(2mu_(lambda_ + mu_))
       * tg = c_2(4mu_ * Ae(u)*e(u) + (lambda_ - mu_)tr(Ae(u))tr(e(u)))
       * entries of e(u) are contained in forward_strain[0 - 2]
       * */
      assert(elem_forw.GetSize() == 3);
      assert(elem_adj.GetSize() == 3);
   
      tr_eu2 = elem_forw[0] + elem_forw[1];
      tr_eu2 *= elem_adj[0] + elem_adj[1];
      tr_eu2 *= lambda_ - mu_;
      tr_eu2 *= 2.0*lambda_ + 2.0*mu_;

      topG =  (lambda_ + 2.0*mu_); 
      topG *= (elem_forw[0]*elem_adj[0] + elem_forw[1]*elem_adj[1]);
      topG += lambda_ * (elem_forw[1]*elem_adj[0] + elem_forw[0]*elem_adj[1]);
      topG += 4.0 * mu_ * (elem_forw[2]*elem_adj[2]);
      topG *= 4.0 * mu_;
      
      topG += tr_eu2;
      
      topG *= c2;
      break;

    case FULL: // top-grad 3d *neumann*
      // FIXME not correct if forw != adjoint
      assert(elem_forw.GetSize() == 6);
      tr_eu2 = elem_forw[0] + elem_forw[1] + elem_forw[2];
      tr_eu2 *= tr_eu2;

      topG =  20.0 * mu_ * lambda_;
      topG += 9.0 * lambda_ * lambda_;
      topG -= 4.0 * mu_ * mu_;
      topG *= tr_eu2;

      topG += mu2 * (elem_forw[0] * elem_forw[0]);
      topG += mu2 * (elem_forw[1] * elem_forw[1]);
      topG += mu2 * (elem_forw[2] * elem_forw[2]);

      topG += 2.0 * mu2 * (elem_forw[3] * elem_forw[3]);
      topG += 2.0 * mu2 * (elem_forw[4] * elem_forw[4]);
      topG += 2.0 * mu2 * (elem_forw[5] * elem_forw[5]);

      topG *= c3;
      break;

    case AXI:
    case PLANE_STRESS:
    default:
      EXCEPTION("SubType not implemented!");
  }
  return topG;
}

double TopGrad::CalcPoissonTopGradOnElement(const unsigned int e) const
{
  // FIXME xiDx is constant for all elements!!!
  // FIXME move out of the loop
  double partElemMat(0.0);
  Matrix<double> xiDx;
  Matrix<double> xiDxT;
  Matrix<double> ptCoord;

  assert(false);
  /* FIXME
  
  Elem* de = optimization->GetDesign()->data[e].elem; // FIXME loop over all elems!!
  
  const Vector<double> &intWeights = de->ptElem->GetIntWeights();
  double jacobi(0.0);
  
  domain->GetGrid()->GetElemNodesCoord(ptCoord, de->connect); //, coordUpdate_);

  for(unsigned int ap = 1, nr = de->ptElem->GetNumIntPoints(); ap <= nr; ++ap)
  {
    // xiDx is 4x2 matrix
    de->ptElem->GetGlobDerivShFncAtIp(xiDx, ap, ptCoord, jacobi, de);
    xiDx.Transpose(xiDxT);
    // v is something like grad(u)
    Vector<double> v(2);
    xiDxT.Mult(elem_forw, v);
    
    partElemMat  = v[0] * (v[0] * pol[0][0] + v[1] * pol[1][0]);
    partElemMat += v[1] * (v[0] * pol[0][1] + v[1] * pol[1][1]);
    
    partElemMat *= intWeights[ap-1] * jacobi;
  }
  */
  return -partElemMat; // FIXME sign??
}

void TopGrad::CalcTopGrads(SubTensorType sub, App::Type app)
{
  for(unsigned int e = 0; e < numelems_; ++e)
  {
    // in poisson case, the subtype PLANE_STRAIN is ignored
    optimization->GetElementSolution(elem_forw, elem_adj, e, sub, app);
    // do the actual calculation of the topgrad
    double tg(0.0);
    switch(app)
    {
    case App::MECH:
      tg = CalcTopGradOnElement();
      break;
    case App::HEAT:
      tg = CalcPoissonTopGradOnElement(e);
      break;
    default:
      assert(false);
    }
    
    // save topgradvalues and element numbers
    topGrads.push_back(ElementValues(e, tg));
    elements[e].value = tg;
  }
}

} // end of namespace
