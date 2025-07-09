#include <cmath>
#include <complex>
#include <cstdlib>
#include <iostream>

#include "DataInOut/Logging/LogConfigurator.hh"
#include "Domain/Domain.hh"
#include "Domain/ElemMapping/Elem.hh"
#include "Domain/Mesh/Grid.hh"
#include "FeBasis/BaseFE.hh"
#include "General/defs.hh"
#include "General/Environment.hh"
#include "General/Exception.hh"
#include "MatVec/SingleVector.hh"
#include "MatVec/exprt/xpr1.hh"
#include "MatVec/exprt/xpr2.hh"
#include "MatVec/Matrix.hh"
#include "Optimization/Condition.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Function.hh"
#include "Optimization/LevelSet.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/ShapeGrad.hh"
#include "Utils/Point.hh"
#include "Utils/tools.hh"

namespace CoupledField
{
using std::abs;
using std::endl;
using std::cout;
using std::flush;
using std::min;
using std::max;
using std::vector;
using std::pair;
using std::make_pair;

DEFINE_LOG(ls, "levelSet")

/** Set in Optimization::SetEnums() */
Enum<LevelSet::Action::Type>  LevelSet::Action::type;

bool IsClose(const double d1, const double d2)
{ 
  return abs(d1-d2) < 1e-8;
}

bool HasIntersection(const LevelSetNode *node1, const LevelSetNode *node2)
{
  assert(node1 != NULL && node2 != NULL); // this should only be called with existing nodes
  
  if((node1->value == 0.0) && (node2->value == 0.0)) return true;
  if((node1->value == 0.0) || (node2->value == 0.0)) return false;
  // no zero left
  return (node1->value > 0.0) != (node2->value > 0.0);
}

bool HasIntersection(const lsnodepair &p)
{
  return HasIntersection(p.first, p.second);
}

/******** LEVELSETNODE ***************/
LevelSetNode::LevelSetNode(const double val, const unsigned int cfs_number) :
  state(LS_NONE), value(val), phi_temp(0.0), intersection_length(0.0),
  f_ext(0.0), shapegrad(0.0), global_node_number(cfs_number)
{
  static const unsigned int dim(domain->GetDim());
  neighbours_.resize((dim == 2 ? 4 : 6), NULL);
}

bool LevelSetNode::HasAcceptedNeighbour(const LevelSetNode::State state_accepted) const
{
  // check if a neighbour node is located at the boundary -> needed for fast marching
  for(unsigned int i = 0, nei_size = neighbours_.size(); i < nei_size; ++i)
  {
    if(neighbours_[i] == NULL) continue;
		if(neighbours_[i]->state == state_accepted)	return true;
  }
  return false;
}

bool LevelSetNode::IsBoundary() const
{
	// node is at the boundary if it has a value of 0.0
  // or the 0-levelset intersects the edge to a neighbour
	if(this->value == 0.0) return true;
  for(unsigned int i = 0, nei_size = neighbours_.size(); i < nei_size; ++i)
  {
    if(neighbours_[i] == NULL) continue;
    if(HasIntersection(this, neighbours_[i])) return true;
  }
  return false;
}
/******** LEVELSETNODE ***************/


/******** LEVELSETELEMENT ***************/
LevelSetElement::LevelSetElement(DesignElement* de, const double val) :
  de_(de),
  shapeGradValue(0.0)
{
  // it is important here to call resize, not reserve, because we assume the existence
  // of this object at a later point where we acquire a reference to it
  const unsigned int dim(domain->GetDim());
  nodes_.resize((dim == 2 ? 4 : 8), 0);
  objects.reserve((dim == 2) ? 2 : 4);
}

void LevelSetElement::UpdateDesign()
{
  // walk over all attached nodes and sum up
  // 1. absvalue: |val| (for denom)
  // 2. negvalue: all negative levelset values from the nodes
  // with 1 and 2: we get the correct percentage of the coverage of the element
  // by computing negvalue/absvalue;
  double absvalue(0.0), negvalue(0.0);
  for(unsigned int i = 0, nodes_size = nodes_.size(); i < nodes_size; ++i)
  {
    const double val(nodes_[i]->value);
    absvalue += abs(val);
    if(val < 0.0) negvalue -= val;
  }

  // now update the design
  const double low(de_->GetLowerBound());
  if(IsClose(absvalue, 0.0))
  {   
    // this should not happen very often
    LOG_DBG(ls) << "element " << de_->elem->elemNum << " had only nodes with tiny levelset values";
    // if negvalue is 0.0, then we have only positive levelset values in all nodes
    // which corresponds to material
    negvalue = (negvalue == 0.0 ? de_->GetUpperBound() : low );
  }
  else
  {
    negvalue /= absvalue;
  }
  de_->SetDesign(low > negvalue ? low : negvalue);
}

bool LevelSetElement::ContainsFront() const
{
  assert(nodes_[0] != NULL && nodes_[1] != NULL && nodes_[2] != NULL && nodes_[3] != NULL);
  // test the edges of the element
  // front face in 3D, whole element in 2D
  for(unsigned int i = 0; i < 4; ++i)
    if(HasIntersection(GetEdgeNodes(i))) return true;

  if(domain->GetDim() == 3)
  {
    assert(nodes_[4] != NULL && nodes_[5] != NULL && nodes_[6] != NULL && nodes_[7] != NULL);
    // back face and middle edges
    for(unsigned int i = 4; i < 12; ++i)
      if(HasIntersection(GetEdgeNodes(i))) return true;
  }
  
  return false; // no sign change within our nodes
}

const lsnodepair LevelSetElement::GetEdgeNodes(const unsigned int idx) const
{
  // check the index, idx >= 0 trivially fulfilled because of unsigned
  assert((domain->GetDim() == 2) ? (idx < 4) : (idx < 12));

  switch(idx)
  {
  case 0:
    return make_pair(nodes_[0], nodes_[1]);
  case 1:
    return make_pair(nodes_[1], nodes_[2]);
  case 2:
    return make_pair(nodes_[2], nodes_[3]);
  case 3:
    return make_pair(nodes_[3], nodes_[0]);
  case 4:
    return make_pair(nodes_[4], nodes_[5]);
  case 5:
    return make_pair(nodes_[5], nodes_[6]);
  case 6:
    return make_pair(nodes_[6], nodes_[7]);
  case 7:
    return make_pair(nodes_[7], nodes_[4]);
  case 8:
    return make_pair(nodes_[0], nodes_[4]);
  case 9:
    return make_pair(nodes_[1], nodes_[5]);
  case 10:
    return make_pair(nodes_[2], nodes_[6]);
  case 11:
    return make_pair(nodes_[3], nodes_[7]);
  default:
    EXCEPTION("wrong index " << idx);
  }
}

void LevelSetElement::CalcShapeGrad(const Vector<double> &elem_sol, linElastInt* bdb_form)
{
  // forget former intersection objects and reset shapegrad values
  objects.clear();

  // elem_sol = u1_x u1_y u2_x u2_y  ... u4_x u4_y
  // if we are not at the boundary we do nothing
  if(!ContainsFront()) return;

  CalcIntersectionObjects(elem_sol);
  
  double length(0.0); // of intersection line
  for(unsigned int k = 0; k < objects.size(); ++k)
  {
    // do this for both nodes at the end of the line!
    length += IntegrateIntersectionObject(objects[k], bdb_form);
    LOG_DBG3(ls) << "integrating intersection object on element " << objects[k].elemNum;
    LOG_DBG3(ls) << "obtained value = " << objects[k].integralValue;
  }
  
  // average the shape gradient values in the intersection objects and
  // save the result in the member variable
  const unsigned int numObj(objects.size());
  for(unsigned int i = 0; i < numObj; ++i)
    shapeGradValue += objects[i].integralValue;
  
  shapeGradValue /= static_cast<double>(numObj);
  
  // FIXME brutally write this value into the attached levelset nodes
  // if they do not already have a value = first come first serve
  for(unsigned int i = 0, s = nodes_.size(); i < s; ++i)
  {
    if(nodes_[i]->shapegrad == 0.0)
    {
      nodes_[i]->shapegrad = shapeGradValue;
      nodes_[i]->intersection_length = length;
    }
  }
}

void LevelSetElement::CalcIntersectionObjects(const Vector<double> &elem_sol)
{
  static const unsigned int dim(domain->GetDim());
  assert(dim == 2); // FIXME
  
  // calc the intersection objects
  vector<Vector<double> > displacements;
  displacements.reserve(4);
  vector<unsigned int> indices;
  indices.reserve(4);
  vector<Vector<double> > intersection_points;
  intersection_points.reserve(4);

  for(unsigned int outer = 0; outer < 4; ++outer)
  {
    const lsnodepair p(GetEdgeNodes(outer));
    if(HasIntersection(p))
    {
      // check if both node values are 0.0, this is a special case that almost always
      // happens only directly after constructing the levelset elements!
      // we have to add two intersection points so that we may get the line
      // connecting the 0.0-nodes
      if(p.first->value == 0.0 && p.second->value == 0.0)
      {
        // we have to remember both end points of the line
        Vector<double> point1;
        Vector<double> point2;
        domain->GetGrid()->GetNodeCoordinate(point1,  p.first->GetNumber());
        domain->GetGrid()->GetNodeCoordinate(point2, p.second->GetNumber());

        Vector<double> disp1;
        for(unsigned int j = 0; j < dim; ++j)
        {
          disp1.Push_back(elem_sol[2*outer + j]);
        }

        unsigned int idx(2*outer + dim);
        if(!(idx < elem_sol.GetSize())) idx = 0;

        Vector<double> disp2;
        for(unsigned int j = 0; j < dim; ++j)
        {
          disp2.Push_back(elem_sol[idx + j]);
        }

        intersection_points.push_back(point1);
        intersection_points.push_back(point2);
        displacements.push_back(disp1);
        displacements.push_back(disp2);
        indices.push_back(2*outer);
        indices.push_back(idx);
      }
      else
      {
        // we get one intersection point at max for each edge
        Vector<double> point;
        double weight = CalcBarycenterNode(point, p);
        intersection_points.push_back(point);
        Vector<double> tmpdisp;
        CalcBarycenterDisplacement(tmpdisp, weight, elem_sol, 2*outer);
        displacements.push_back(tmpdisp);
        indices.push_back(outer*2);
      }
    }
  } // end loop over edges


  // we now have to build the intersection objects from 
  // the nodes contained in intersection_points
  // we do this by taking two consecutive nodes and building an intersection object with them
  if(intersection_points.size() != 2 && intersection_points.size() != 4)
  {
    LOG_DBG(ls) << "element " << this->de_->elem->elemNum << " has " << intersection_points.size() 
                << " intersection points";
    return;
  }
  assert(intersection_points.size() == 2 || intersection_points.size() == 4);
  assert(intersection_points.size() == indices.size());
  assert(displacements.size() == indices.size());

  // from values obtained above, build intersection object
  for(unsigned int l = 0; l < displacements.size()/2; ++l)
  {
    IntersectionObject tmpObj;
    const unsigned int idx(2*l);

    tmpObj.elemNum = de_->elem->elemNum;
    tmpObj.integralValue = 0.0;
    tmpObj.displacements.assign(displacements.begin() + idx, displacements.begin() + idx + 2);
    tmpObj.points.assign(intersection_points.begin() + idx, intersection_points.begin() + idx + 2);
    tmpObj.indices.assign(indices.begin() + idx, indices.begin() + idx + 2);

    this->objects.push_back(tmpObj);
  }
}

double LevelSetElement::CalcBarycenterNode(Vector<double> &outpoint, const lsnodepair &p) const
{
  // find points
  Vector<double> c1;
  Vector<double> c2;

  domain->GetGrid()->GetNodeCoordinate(c1, p.first->GetNumber());
  domain->GetGrid()->GetNodeCoordinate(c2, p.second->GetNumber());
  outpoint.Resize(c1.GetSize());
  
  assert(c1.GetSize() == c2.GetSize());

  // find intersection weight
  double weight(0.0);
  if(p.first->value != 0.0 || p.second->value != 0.0) // FIXME could be wrong if node2->GetValue() == 0.0
  {
    weight = abs(p.first->value);
    weight /= (abs(p.first->value) + abs(p.second->value));
  }

  // outpoint is a barycentrical combination of c1 and c2
  for(unsigned int i = 0, n = c1.GetSize(); i < n; ++i)
  {
    outpoint[i]  = c1[i];
    outpoint[i] -= c2[i];
    outpoint[i] *= weight;
    outpoint[i] += c2[i];
  }
  
  // return weight, because we need it for the displacement
  return weight;
}

void LevelSetElement::CalcBarycenterDisplacement(Vector<double> &out_u, const double w1,
                                                 const Vector<double> &elem_sol, const int idx) const
{
  static const unsigned int dim(domain->GetDim());
  out_u.Resize(dim);
  for(unsigned int i = 0; i < dim; ++i)
  {
    out_u[i] = w1 * elem_sol[idx + i] + (1.0-w1) * elem_sol[idx + i];
  }
}

double LevelSetElement::IntegrateIntersectionObject(IntersectionObject &o, linElastInt *bdb_form)
{
  assert(false);
  /* FIXME
  static const unsigned int dim(domain->GetDim());
  assert(dim == 2); // FIXME

  Matrix<double> B; // 2D: 3 * 8
  Matrix<double> C; // = [c] or E_ijkl or dMat in BDBInt or A in shape grad papers
  Vector<double> b_u;
  Vector<double> c_b_u;
  
  for(unsigned int num = 0; num < 2; ++num)
  {

    // set B
    bdb_form->CalcBMatOnly(B, o.points[num], de_->elem);
    // LOG_DBG3(ls) << "B = " << B.ToString();

    // set c(\rho) !!!
    bdb_form->calcDMat(C, de_->elem);

    // in the full element u = 2*4 in lin 2D we have all 0 but u at idx
    Vector<double> full_u(8, 0.0);
    // set correct indices != 0
    for(unsigned int d = 0; d < dim; ++d)
      full_u[o.indices[num] + d] = (o.displacements[num])[d];

    b_u.Resize(C.GetNumCols());
    b_u = B*full_u;
    // must multiply with [c] here
    c_b_u.Resize(C.GetNumCols());
    c_b_u = C * b_u;
    double integrand = b_u * c_b_u;

    o.integralValue += integrand;
  }
  o.integralValue *= 0.5;
  
  // return volume
  Vector<double> difference(o.points[0]);
  difference -=  o.points[1];
  return difference.NormL2();
  */
  return -1.0;
}

const std::string ToString(const LevelSetElement &elem)
{
  std::stringstream ss;
  ss << "levelset element[" << elem.de_->elem->elemNum << "]";
  const unsigned int nodes_size(elem.nodes_.size());
  ss << " consists of " << nodes_size << " nodes:" << endl;
  for(unsigned int i = 0; i < nodes_size; ++i)
  {
    const LevelSetNode *node(elem.nodes_[i]);
    ss << "\t" << i+1 << ". number: " << node->GetNumber() << ": value = " << node->value;
    ss << "; neighbour node numbers = [";
    ss << "r ";
    if(node->GetNeighbour(VicinityElement::X_P) != NULL)
      ss << node->GetNeighbour(VicinityElement::X_P)->GetNumber();
    else ss << "NA";
    ss << ", l ";
    if(node->GetNeighbour(VicinityElement::X_N) != NULL)
      ss << node->GetNeighbour(VicinityElement::X_N)->GetNumber();
    else ss << "NA";
    ss << ", u ";
    if(node->GetNeighbour(VicinityElement::Y_P) != NULL)
      ss << node->GetNeighbour(VicinityElement::Y_P)->GetNumber();
    else ss << "NA";
    ss << ", d ";
    if(node->GetNeighbour(VicinityElement::Y_N) != NULL)
      ss << node->GetNeighbour(VicinityElement::Y_N)->GetNumber();
    else
      ss << "NA";
    if(domain->GetDim() == 3)
    {
      ss << ", f ";
      if(node->GetNeighbour(VicinityElement::Z_P) != NULL)
        ss << node->GetNeighbour(VicinityElement::Z_P)->GetNumber();
      else ss << "NA";
      ss << ", b ";
      if(node->GetNeighbour(VicinityElement::Z_N) != NULL)
        ss << node->GetNeighbour(VicinityElement::Z_N)->GetNumber();
      else ss << "NA";
    }

    ss <<"]" << endl;
  }
  return ss.str();
}
/******** END OF LEVELSETELEMENT ***************/



LevelSet::LevelSet(Optimization* opt, PtrParamNode pn) : 
  last_node_index_(0), dump_fast_marching_(false), dump_lselement_(0)
{
  cout << "Levelset: " << flush;
  optimization = dynamic_cast<ShapeGrad*>(opt);
  assert(optimization != NULL);
  design_ = &optimization->GetDesign()->data;
  assert(design_ != NULL);

  // cache the element widths, assumes a uniform grid!!
  domain->GetGrid()->GetElemShapeMap((*design_)[0].elem, false)->GetEdgeLength(edge_length_);
  assert(edge_length_.GetSize() == (domain->GetDim() == 2 ? 2 : 3));
  LOG_DBG(ls) << "edge_length_: x = " << edge_length_[0] << ", y = " << edge_length_[1];
  if(edge_length_.GetSize() == 3)
  {
    LOG_DBG(ls) << "edge_length_: z = " << edge_length_[2];
  }

  if(pn != NULL)
  {
    // reduce to our actual ParamNode
    pn = pn->Get("levelSet");
    dump_fast_marching_ = pn->Get("ls_dump_fast_marching")->As<bool>();
    dump_lselement_ = pn->Get("ls_element")->As<Integer>();
    time_step_ = pn->Get("ls_time_step")->As<Double>();
    bool make_trivial_holes(false);
    if(pn->Has("actionList"))
    {
      ParamNodeList list = pn->Get("actionList")->GetChildren();
      const unsigned int list_size(list.GetSize());
      LOG_DBG(ls) << "read " << list_size << " actions from xml";
      // FSFSFS: maybe give some defaults if list_size == 0
      action_.reserve(list_size);
      for(unsigned int i = 0; i < list_size; ++i)
      {
        action_.push_back(Action(list[i]));
        if(action_.back().GetType() == Action::TRIVIAL_HOLE) make_trivial_holes = true;
      }
    }
    
    if(make_trivial_holes)
    {
      if(pn->Has("nodeList"))
      {
        ParamNodeList list = pn->Get("nodeList")->GetChildren();
        const unsigned int list_size(list.GetSize());
        LOG_DBG(ls) << "trivial hole: read " << list_size << " node numbers from xml";
        trivial_holes_.reserve(list_size);
        for(unsigned int i = 0; i < list_size; ++i)
        {
          trivial_holes_.push_back(list[i]->Get("value")->As<Integer>());
          LOG_DBG(ls) << " node " << i+1 << ": number " << trivial_holes_[i];
        }
        assert(trivial_holes_.size() == list_size);
      }
    }
  }
  else
  {
    // FSFSFS: maybe substitute this with a useful default
    throw Exception("the element <levelSet> is mandatory for level-set optimization");
  }

  // This is just outsourced constructor stuff
  SetupSpace();
}

void LevelSet::SolveProblem(const unsigned int iter)
{
  if(iter == 1)
  {
    cout << "Start Level-Set optimization ..." << endl;
  }

  Action* action = NULL;
  if((action = IsTriggered(Action::SIGNED_DISTANCE_FIELD, iter)) != NULL)
  {
    //BuildSignedDistanceFunction(dump_fast_marching_);
		for(unsigned int p = 0, max = action->GetPerform(); p < max; ++p)
			BuildSignedDistanceFunction(false);
  }

  if((action = IsTriggered(Action::TRIVIAL_HOLE, iter)) != NULL)
  {
    cout << "  creating some trivial holes" << endl;

    //std::for_each(trivial_holes_.begin(), trivial_holes_.end(), MakeTrivialHole());
    for(unsigned int j = 0, nr_holes = trivial_holes_.size(); j < nr_holes; ++j)
      MakeTrivialHole(trivial_holes_[j]);

    // after creating the hole, we have to rebuild the signed distance
    LOG_DBG(ls) << " BuildSignedDistanceFunction() after MakeTrivialHole()";
    BuildSignedDistanceFunction(dump_fast_marching_);
    UpdateDesignForAllLevelSetElements();
  }

  if((action = IsTriggered(Action::DO_SHAPE_STEP, iter)) != NULL)
  {
    cout << "  transporting levelset" << endl;
    
    // since this function is only called when doing shape optimization, 
    // we calculate the shape gradients here
    CalcShapeGradientOnAllElements();

    // now we can transport the levelset
    for(unsigned int p = 0, max = action->GetPerform(); p < max; ++p)
      TransportLevelSet(time_step_);
    
    BuildSignedDistanceFunction(false);
    UpdateDesignForAllLevelSetElements();
  }
}

void LevelSet::BuildSignedDistanceFunction(const bool really_dump)
{
  cout << "  build signed distance" << endl;

  // fill the sets with pointers to the nodes
  bool one_accepted(false);
  one_accepted = SetSignedDistanceNodeStates();
  
  if(!one_accepted)
  {
    LOG_DBG(ls) << "no boundary nodes could be found";
    
    // set all levelset values to a reasonable default 
    const double new_val((nodes_[0].value < 0.0) ? -1.0 : 1.0);
    // in this case there is not one accepted node = boundary node left
    for(unsigned int n = 0, ss = nodes_.size(); n < ss; ++n)
    {
      nodes_[n].value = new_val;
    }
    
    return;
  }

  bool negative(true);
  ExecuteFastMarching( negative, really_dump);
  ExecuteFastMarching(!negative, really_dump);
}

void LevelSet::MakeTrivialHole(const unsigned int node_nr)
{
  // node_nr is the cfs node number!
  LOG_DBG(ls) << "calling MakeTrivialHole() with node number " << node_nr;
  std::map<unsigned int, LevelSetNode*>::iterator it(lsnode_map_.find(node_nr));

  if(it == lsnode_map_.end())
  {
    LOG_DBG(ls) << "node number " << node_nr << " could not be found in the node_map!";
    return;
  }
  
  LevelSetNode *canditate(it->second);
  assert(canditate != NULL);

  // if we are outside the domain we cannot make a hole
  if(canditate->value > 0.0)
    return;

  double minimum(min(edge_length_[0], edge_length_[1]));
  assert(domain->GetDim() == edge_length_.GetSize());
  if(domain->GetDim() == 3)
  { 
    minimum = min(minimum, edge_length_[2]);
  }    
  // this node will be in accepted so we set the correct value here
  canditate->value = 0.5 * minimum;
  for(unsigned int k = 0, size = canditate->neighbours_.size(); k < size; ++k)
  {
    if(canditate->neighbours_[k] == NULL) continue;

		if(canditate->neighbours_[k]->value < 0.0)
			canditate->neighbours_[k]->value = -0.5 * edge_length_[GetAxisFromIndex(k)];
		else
			canditate->neighbours_[k]->value =  0.5 * edge_length_[GetAxisFromIndex(k)];
  }

  // immediately commit this iteration
  if(dump_fast_marching_)
    optimization->CommitIteration();
}

LevelSetNode* LevelSet::GetNodePointer(const unsigned int nodeNumber) const
{
  // this function needs the cfs node number!!!
  std::map<unsigned int, LevelSetNode*>::const_iterator ittmp = lsnode_map_.find(nodeNumber);
  assert(ittmp != lsnode_map_.end());
  
  return ittmp->second;
}

double LevelSet::GetGradientAtNode(const unsigned int node_nr, const unsigned int dir) const
{
  // 0 = X_P, 1 = X_N, 2 = Y_P, 3 = Y_N, 4 = Z_P, 5 = Z_N
  const LevelSetNode *node(&nodes_[node_nr]);
  assert(node != NULL);
  // assume that we have at least one neighbour in each direction
  assert((node->neighbours_[0] != NULL) || (node->neighbours_[1] != NULL));
  assert((node->neighbours_[2] != NULL) || (node->neighbours_[3] != NULL));
  assert(domain->GetDim() != 3 ? true : ((node->neighbours_[4] != NULL) || (node->neighbours_[5] != NULL)));
  
  if(node->neighbours_[dir] != NULL)
    return (node->neighbours_[dir]->value - node->value) / edge_length_[GetAxisFromIndex(dir)];
  
  unsigned int other_index(dir + 1);
  if(dir == 1 || dir == 3 || dir == 5) 
    other_index -= 2;
  
  return (node->value - node->neighbours_[other_index]->value) / edge_length_[GetAxisFromIndex(dir)];
}

double LevelSet::GetNormGradientAtNode(const unsigned int node_nr) const
{
  double gx(GetGradientAtNode(node_nr, 0));
  double gy(GetGradientAtNode(node_nr, 2));
  double gz(0.0);
  if(domain->GetDim() == 3)
    gz = GetGradientAtNode(node_nr, 4);
  gz = std::sqrt(gx*gx + gy*gy + gz*gz);
  LOG_DBG3(ls) << "gradient at node " << node_nr << ": " << gz;
  return gz;
}

/** Debug output for the complete space */
const std::string ToString(const vector<LevelSetElement> &space)
{
  std::stringstream sstr;
  for(unsigned int i = 0, ss = space.size(); i < ss; ++i)
    sstr << "#" << i << "=" << ToString(space[i]);
  return sstr.str();
}

void LevelSet::SetupSpace()
{
  // we need the vicinity in the design elements set (if not set yet)

  // Assume no periodic B.C. Otherwise organize DesignStructure from ErsatzMaterial or create
  // own instance!
  VicinityElement::Init(optimization->GetDesign(), NULL);
  assert((*design_)[0].vicinity != NULL);
  
  cout << "Levelset space... " << flush;

  // next we create a levelset node for every gridnode in the correct region
  StdVector<RegionIdType> reg(1);
  reg[0] = optimization->GetDesign()->GetRegionId();
  const unsigned int region_size(domain->GetGrid()->GetNumNodes(reg));
  // if we do not reserve enough, we risk to generate a copy of this vector
  // which would invalidate all pointers! 
  nodes_.reserve(region_size+1);
  
  // create our level set element space.
  // The number of level set elements equals the number of design nodes.
  // e.g.  9 design elements = 9 finite elements
  //    = 16 nodes for the finite elements
  //    = 16 level set nodes = 9 design elements
  //    =  9 level set elements
  const unsigned int design_size(design_->GetSize());
  space_.reserve(design_size);

  unsigned int factor(1);
  for(unsigned int i = 0; i < design_size; ++i)
  {
    if(5*i == factor*design_size) //draw a dot every 1/5 elements to cout
    {
      cout << "." << flush;
      ++factor;
    }
    // in this algorithm we define a design node to be the upper right node of the
    // to be created element. If we have a +x or +y ghost neighbour we create
    // and additional element if both we also create a corner element. (see example above!)
    DesignElement *de = &(*design_)[i];
    assert(de != NULL);

    // add a levelset element for the current designelement
    LevelSetElement lse(de);

    // compute the correct order of the levelset nodes and attach them to the lse
    // creates the levelset nodes and a map of global node number with pointer to the lsn
    AddOrderedLevelsetNodesToLevelsetElement(lse);

    space_.push_back(lse);
    // link the levelset element to the designelement
    de->lse_ = &space_.back();
  }

  // we now need the neighbours of the levelset nodes as well
  // walk over all elements again
  // build neighbourhood according to information from the element vicinity
  // and based on the sorting done above
  // at the same time, we can update the levelset values in the levelset nodes
  // and initialize the set accepted
  factor = 1;
  for(unsigned int i = 0; i < design_size; ++i)
  {
    if(5*i == factor*design_size) //draw a dot every 1/5 elements to cout
    {
      cout << "." << flush;
      ++factor;
    }
    BuildNeighbourhoodOfLevelsetNodes(&(*design_)[i]);
  }
  cout << "done." << endl;

  // debug: for one element show all the neighbours
  assert(dump_lselement_ < space_.size()); // not very nice
  LOG_DBG(ls) << ToString(space_[dump_lselement_]);
}

void LevelSet::ExecuteFastMarching(const bool negative, const bool really_dump)
{
  LOG_DBG(ls) << "fast marching mode is " << (negative ? "NEGATIVE" : "POSITIVE");
  
  // make some tests more readable
  LevelSetNode::State state_accepted = (negative ? LevelSetNode::LS_ACCEPTED : LevelSetNode::LS_ACCEPTED_POS);
  LevelSetNode::State state_close = (negative ? LevelSetNode::LS_CLOSE : LevelSetNode::LS_CLOSE_POS);
  LevelSetNode::State state_far = (negative ? LevelSetNode::LS_FAR : LevelSetNode::LS_FAR_POS);
  
  // take a node from close
  while(true)
  {
    if(really_dump)
    {
      optimization->CommitIteration();
      LOG_DBG3(ls) << "iteration = " << optimization->GetCurrentIteration();
    }
    
    /* steps 1 through 4 are taken from Adalsteinsson/Sethian:
    * the fast construction of extension velocities in level set methods */
    // 1. let trial be the element in close with the largest value (or smallest if we are positive)
    LevelSetNode *trial(NULL);
    unsigned int start(0), number_nodes(nodes_.size());
    // find the first close node
    for(start = 0; start < number_nodes; ++start)
    {
      if(nodes_[start].state == state_close)
			{
				trial = &nodes_[start];
				break;
			}
    }
    if(trial == NULL) break; // no more close points left

    for(unsigned int n = start; n < number_nodes; ++n)
    {
      if(!(nodes_[n].state == state_close))
        continue;
 
			if(negative)
			{
				// in the negative case we need the largest value
				if(trial->value < nodes_[n].value)
					trial = &nodes_[n];
			}
			else
			{
				if(trial->value > nodes_[n].value)
					trial = &nodes_[n];
			}
    }
    assert(trial != NULL);
    assert(trial->state == state_close);
    LOG_DBG3(ls) << "trial is (" << trial->GetNumber() << ", " << trial->value << ")";
    
    static const unsigned int ss(trial->neighbours_.size()); // number of neighbours never changes, so static
    
    // 2. move all neighbours of trial that are in far into close
    for(unsigned int k = 0; k < ss; ++k)
    {
			LevelSetNode *curr_nei(trial->neighbours_[k]);
      if(curr_nei == NULL)
        continue;
      
      if(curr_nei->state == state_far)
      {
        curr_nei->state = state_close;
        LOG_DBG3(ls) << "   -> found a neighbour of trial in far (" << curr_nei->GetNumber() << ", " 
                     << curr_nei->value << ") -> moving to close (state = " << curr_nei->state << ")";
      }
    }
    
    // 3. recompute the value of the levelset function at each neighbour of trial in close
    // reminder: ss = trial->neighbours_.size());
    for(unsigned int outer = 0; outer < ss; ++outer)
    {
			assert(trial->state == state_close);
      if(trial->neighbours_[outer] == NULL)
        continue;
      
      if(trial->neighbours_[outer]->state == state_close)
      {
        // check if we need to update the value at all
        bool update_value(false);
        // new_value = "infty" because we need to redefine it
        // FIXME this value set here should not be used!!!
        double new_value((negative ? -1.0 : 1.0) * 100);
        // trial has neighbour in close
        // recompute levelset value of this node
        // find neighbour node in accepted
        for(unsigned int inner = 0; inner < ss; ++inner)
        {
					LevelSetNode *ptrNei(trial->neighbours_[outer]->neighbours_[inner]);
          if(ptrNei == NULL) continue;

          if(ptrNei->state == state_accepted)
          {
            update_value = true;
            // the edge length can differ in x and y (and z) direction
            if(negative)
              new_value = max(new_value, ptrNei->value - edge_length_[GetAxisFromIndex(inner)]);
            else
              new_value = min(new_value, ptrNei->value + edge_length_[GetAxisFromIndex(inner)]);
          }
				} // end of inner loop
        
        if(!update_value)
        {
          // here, trial->neighbours_[outer] is 
          // a point with no accepted neighbours!
          // calculate value from trial
          // if we don't do this, then there exist cases where the algorithm fails
          update_value = true;
          // the edge length can differ in x and y (and z) direction
          if(negative)
            new_value = max(new_value, trial->value - edge_length_[GetAxisFromIndex(outer)]);
          else
            new_value = min(new_value, trial->value + edge_length_[GetAxisFromIndex(outer)]);
        }

        // if(update_value) // update_value is always true here
        // {
        trial->neighbours_[outer]->value = new_value;
        LOG_DBG3(ls) << "   -> new levelset value for (" << trial->neighbours_[outer]->GetNumber() 
        << ", " << trial->neighbours_[outer]->value << ")";
        // }
      }
    } // end of outer loop

    // 4. move trial into accepted and continue from the top of the loop
    assert(trial->state == state_close);
    trial->state = state_accepted;
    LOG_DBG3(ls) << "  -> accepted node (" << trial->GetNumber() << ", " << trial->value << ")";
  } // end of while(!close.empty())
}

bool LevelSet::SetSignedDistanceNodeStates()
{
  const unsigned int nodes_size(nodes_.size());
  
  // clear the old states first
  for(LevelSetNode& node : nodes_)
    node.state = LevelSetNode::LS_NONE;
  
  // find out if we have at least one node at the boundary
  bool one_accepted(false);
  // first, find all accepted nodes
  for(unsigned int n = 0; n < nodes_size; ++n)
  {
    LevelSetNode &node = nodes_[n];
    assert(node.state == LevelSetNode::LS_NONE);
    // first find all the nodes for accepted
    if(node.IsBoundary())
    {
      node.state = ((node.value <= 0.0) ? LevelSetNode::LS_ACCEPTED : LevelSetNode::LS_ACCEPTED_POS);
      one_accepted = true;
    }
  }
  
  if(!one_accepted) return one_accepted;

  // third, walk again over all nodes and put them into close or far
  // for the close vectors we compute on the fly the correct levelset value
  // using the levelset values from nodes in accepted
  for(unsigned int n = 0; n < nodes_size; ++n)
  {
    LevelSetNode &node = nodes_[n];
    const bool neg(node.value <= 0.0);
    if(node.state == (neg ? LevelSetNode::LS_ACCEPTED : LevelSetNode::LS_ACCEPTED_POS))
      continue;

    // at this point we are either close or far
    // if there is a neighbour node that is at the boundary we are close
    if(node.HasAcceptedNeighbour(neg ? LevelSetNode::LS_ACCEPTED : LevelSetNode::LS_ACCEPTED_POS))
    { 
      node.state = (neg ? LevelSetNode::LS_CLOSE : LevelSetNode::LS_CLOSE_POS);
      FindNewNodeValueFromNeighbours(&node);
    }
    else
    {
      // node has no accepted neighbour, so we are at far
      // far must be the smallest or largest value, or else the algorithm will not work!
      // FIXME better would be to get the grid dimension and set it to the largest/smallest possible value
      node.state = (neg ? LevelSetNode::LS_FAR : LevelSetNode::LS_FAR_POS);
      node.value = (neg ? -20.0 : 20.0);
    }
  }
  
#ifndef NDEBUG
  for(unsigned int n = 0; n < nodes_size; ++n)
  {
    if(nodes_[n].state == LevelSetNode::LS_NONE)
      EXCEPTION("some nodes have not been assigned a state!");
  }
#endif
  
  return one_accepted;
}

void LevelSet::FindNewNodeValueFromNeighbours(LevelSetNode *nodeptr)
{
  assert(nodeptr != NULL);

  vector<double> dist;
  dist.reserve(6);
  vector<int> orientation;
  orientation.reserve(6);
  for(unsigned int i = 0, nei_size = nodeptr->neighbours_.size(); i < nei_size; ++i)
  {
    LevelSetNode *curr_nei = nodeptr->neighbours_[i];
    if(curr_nei == NULL) continue;
    // check if neighbour is accepted
    if(curr_nei->state == ((curr_nei->value <= 0.0) ? LevelSetNode::LS_ACCEPTED : LevelSetNode::LS_ACCEPTED_POS))
    {
      dist.push_back(abs(curr_nei->value));
      orientation.push_back(GetAxisFromIndex(i));
    }
  }

  assert(dist.size() == orientation.size() && dist.size() > 0);
  double result(100000000.0);
  
  /*
  for(unsigned int k = 0, size = dist.size(); k < size; ++k)
  {
    result = min(result, dist[k] + edge_length_[orientation[k]]);
  }
  */

  // this is taken from
  // adalsteinsson/sethian: fast construction of extension velocities, p. 11
  
  double s(dist[0] + edge_length_[orientation[0]]);
  double t(0.0);
  switch(dist.size())
  {
  case 1:
    result = s;
    break;
  case 2:
    t = dist[1] + edge_length_[orientation[1]];
    if(orientation[0] == orientation[1])
      result = std::min(s, t);
    else
      result = s * t / std::sqrt(s*s + t*t);
    break;
  case 3:
    if(orientation[0] == orientation[1])
    {
      s = std::min(s, dist[1] + edge_length_[orientation[1]]);
      t = dist[2] + edge_length_[orientation[2]];
    }
    if(orientation[0] == orientation[2])
    {
      s = std::min(s, dist[2] + edge_length_[orientation[2]]);
      t = dist[1] + edge_length_[orientation[1]];
    }
    if(orientation[1] == orientation[2])
    {
      t = std::min(dist[2] + edge_length_[orientation[2]], 
          dist[1] + edge_length_[orientation[1]]);
    }
    result = s * t / std::sqrt(s*s + t*t);
    break;
  case 4:
    assert((orientation[0] == orientation[1]) && (orientation[2] == orientation[3]));
    s = std::min(s, dist[1] + edge_length_[orientation[1]]);
    t = std::min(dist[2] + edge_length_[orientation[2]], 
        dist[3] + edge_length_[orientation[3]]);
    result = s * t / std::sqrt(s*s + t*t);
    break;
  default:
    assert(false);
  }
  
  result *= ((nodeptr->value < 0.0) ? -1.0 : 1.0);
  nodeptr->value = result;

}

void LevelSet::AddOrderedLevelsetNodesToLevelsetElement(LevelSetElement &lse)
{
  // we get the node numbers of the corresponding de from cfs
  StdVector<unsigned int> cfs_node_numbers;
  domain->GetGrid()->GetElemNodes(cfs_node_numbers, lse.de_->elem->elemNum);

  // add pointers to the nodes according to the node numbers
  const unsigned int nn(cfs_node_numbers.GetSize());
  assert(nn == (domain->GetDim() == 2 ? 4 : 8));
  vector<Point> point_coords;
  point_coords.reserve(nn);
  
  Point tmpPoint; // Points are always 3D!
  assert(false);
  // FIXME domain->GetGrid()->GetNodeCoordinate(tmpPoint, cfs_node_numbers[0], false);
  double mins[3] = { tmpPoint.data[0], tmpPoint.data[1], tmpPoint.data[2] };
  point_coords.push_back(tmpPoint);
  
  for(unsigned int s = 1; s < nn; ++s)
  {
    // get coordinates for all the points sorted as in cfs_node_numbers
    assert(false);
    // FIXME domain->GetGrid()->GetNodeCoordinate(tmpPoint, cfs_node_numbers[s], false);
    // remember the minimal coordinates for sorting
    for(unsigned int i = 0; i < 3; ++i)
    {
      if(tmpPoint.data[i] < mins[i]) mins[i] = tmpPoint.data[i];
    }
    point_coords.push_back(tmpPoint);
  }
  assert(point_coords.size() == nn);

  vector<unsigned int> node_numbers_sorted;
  node_numbers_sorted.resize(nn, 0);
  // now find out the correct order of the points using the minimal values obtained above
  // we use the lower left corner as first entry, then go counterclockwise
  // we use IsClose() here to compare the coordinates
  for(unsigned int s = 0; s < nn; ++s)
  {
    unsigned int position(0);
    // we are either in 2D, where minz == 0.0 for all points
    // or in 3D in the "front" face of the cube
    if(!IsClose(point_coords[s].data[2], mins[2])) position += 4;

    // if we are at the upper edge, we have to add another offset to the
    // position to get the correct counterclockwise order of the points
    if(!IsClose(point_coords[s].data[1], mins[1]))
    {
      position += 2;
      if(IsClose(point_coords[s].data[0], mins[0])) position += 1;
    }
    else
    {
      if(!IsClose(point_coords[s].data[0], mins[0])) position += 1;
    }
    assert(position < nn); // position >= 0 trivial because of unsigned

    node_numbers_sorted[position] = cfs_node_numbers[s];
  }

  for(unsigned int s = 0; s < nn; ++s)
  {
    // finally we are able to attach the ordered node number to the lse
    const unsigned int node_number(node_numbers_sorted[s]);
    LevelSetNode *lsntmp(NULL);
    std::map<unsigned int, LevelSetNode*>::const_iterator ittmp = lsnode_map_.find(node_number);
    // check if we already have created a levelset node for the given node number
    if(ittmp == lsnode_map_.end())
    {
      // create a new levelset node
      // be careful not to cause a reallocation of nodes_
      LevelSetNode tmplsnode(-1.0, node_number);
      if(nodes_.size() < nodes_.capacity())
      {
        nodes_.push_back(tmplsnode);
      }
      else
      {
        EXCEPTION("We have a problem! Increase the reserved size for nodes_ from LevelSet!! All pointers are now invalid!");
      }

      lsntmp = &nodes_.back();
      // important: insert this new node into the map
      lsnode_map_.insert(std::make_pair(node_number, lsntmp));
    }
    else
    {
      // node already exists, so we use the pointer from the map
      lsntmp = ittmp->second;
    }
    // now add the pointer to this node to the levelset element
    lse.nodes_[s] = lsntmp;
  }
}

void LevelSet::BuildNeighbourhoodOfLevelsetNodes(const DesignElement* de)
{
  assert(de != NULL);
  const bool threeDAndZPIsNULL((domain->GetDim() == 3) && 
                               (de->vicinity->GetNeighbour(VicinityElement::Z_P) == NULL));
  
  // we always have to set the neighbours of this node
  // low left node has index 0 in the front row, index 4 in the back row
  SetNeighboursOfLowLeftNode(0, de);
  if(threeDAndZPIsNULL)
  {
    // Z_P == 0: set nbh of node 4
    SetNeighboursOfLowLeftNode(4, de);
  }

  // check for neighbours of the designelement
  // there are 3 possible cases (in 2D)
  // 1. X_P == 0:
  // also set neighbourhood of lower right node
  if(de->vicinity->GetNeighbour(VicinityElement::X_P) == NULL)
  {
    // low right node has index 1 in the front row, index 5 in the back row
    SetNeighboursOfLowRightNode(1, de);
    if(threeDAndZPIsNULL)
      SetNeighboursOfLowRightNode(5, de);
  }
  
  // 2. Y_P == 0
  // also set neighbourhood of upper left node
  if(de->vicinity->GetNeighbour(VicinityElement::Y_P) == NULL)
  {
    // up left node has index 3 in the front row, index 7 in the back row
    SetNeighboursOfUpLeftNode(3, de);
    if(threeDAndZPIsNULL)
      SetNeighboursOfUpLeftNode(7, de);
  }
  
  // 3. X_P == 0 && Y_P == 0
  // also set neighbourhood of lower right (done), upper left (done) and upper right node
  if((de->vicinity->GetNeighbour(VicinityElement::X_P) == NULL) &&
     (de->vicinity->GetNeighbour(VicinityElement::Y_P) == NULL))
  {
    // up right node has index 2 in the front row, index 6 in the back row
    SetNeighboursOfUpRightNode(2, de);
    if(threeDAndZPIsNULL)
      SetNeighboursOfUpRightNode(6, de);
  }
}

void LevelSet::SetNeighboursOfLowLeftNode(const int startIdx, const DesignElement *de)
{
  assert(startIdx == 0 || startIdx == 4);
  LevelSetNode *ptrnodes = de->lse_->nodes_[startIdx];
  assert(ptrnodes != NULL);

  // set neighbourhood of lower left levelset node which has index 0
  // if we have at least one neighbour missing we know that we are at the boundary
  vector<LevelSetNode*> &low_left = ptrnodes->neighbours_;
  assert(low_left.size() == ((domain->GetDim() == 2) ? 4 : 6));

  // right is always in the same element
  low_left[VicinityElement::X_P] = de->lse_->nodes_[startIdx+1];
  // left neighbour is the lower left node in neighbour element to the left
  DesignElement *tmp_de = de->vicinity->GetNeighbour(VicinityElement::X_N);
  if(tmp_de != NULL)
  {
    low_left[VicinityElement::X_N] = tmp_de->lse_->nodes_[startIdx];
  }
  else
  {
    low_left[VicinityElement::X_N] = NULL;
    ptrnodes->value = 0.0;
  }

  // upper is always in the same element
  low_left[VicinityElement::Y_P] = de->lse_->nodes_[startIdx+3];
  // lower neighbour is the lower left node in neighbour element downwards
  tmp_de = de->vicinity->GetNeighbour(VicinityElement::Y_N);
  if(tmp_de != NULL)
  {
    low_left[VicinityElement::Y_N] = tmp_de->lse_->nodes_[startIdx];
  }
  else
  {
    low_left[VicinityElement::Y_N] = NULL;
    ptrnodes->value = 0.0;
  }

  if((startIdx == 0) && (domain->GetDim() == 3))
  {
    // also add z-neighbours
    assert(low_left.size() == 6);
    low_left[VicinityElement::Z_P] = de->lse_->nodes_[4];
    tmp_de = de->vicinity->GetNeighbour(VicinityElement::Z_N);
    if(tmp_de != NULL)
    {
      low_left[VicinityElement::Z_N] = tmp_de->lse_->nodes_[0];
    }
    else
    {
      low_left[VicinityElement::Z_N] = NULL;
      ptrnodes->value = 0.0;
    }
    return;
  }

  if(startIdx == 4) // only possible in 3D
  {
    assert(low_left.size() == 6);
    low_left[VicinityElement::Z_P] = NULL; // is only called when this is true
    low_left[VicinityElement::Z_N] = de->lse_->nodes_[0];
    ptrnodes->value = 0.0;
  }
}

void LevelSet::SetNeighboursOfLowRightNode(const int startIdx, const DesignElement *de)
{
  // X_P == 0
  assert(startIdx == 1 || startIdx == 5);
  LevelSetNode *ptrnodes = de->lse_->nodes_[startIdx]; // upper left node has index 7
  assert(ptrnodes != NULL);
  ptrnodes->value = 0.0;
  
  vector<LevelSetNode*> & low_right = ptrnodes->neighbours_;
  assert(low_right.size() == ((domain->GetDim() == 2) ? 4 : 6));
   
  low_right[VicinityElement::X_P] = NULL;
  low_right[VicinityElement::X_N] = de->lse_->nodes_[startIdx-1];
  low_right[VicinityElement::Y_P] = de->lse_->nodes_[startIdx+1];
  DesignElement *tmp_de = de->vicinity->GetNeighbour(VicinityElement::Y_N);
  if(tmp_de != NULL)
  {
    low_right[VicinityElement::Y_N] = tmp_de->lse_->nodes_[startIdx];
  }
  else
  {
    low_right[VicinityElement::Y_N] = NULL;
  }
  
  if((startIdx == 1) && (domain->GetDim() == 3))
  {
    assert(low_right.size() == 6);
    low_right[VicinityElement::Z_P] = de->lse_->nodes_[5];
    tmp_de = de->vicinity->GetNeighbour(VicinityElement::Z_N);
    if(tmp_de != NULL)
    {
      low_right[VicinityElement::Z_N] = tmp_de->lse_->nodes_[1];
    }
    else
    {
      low_right[VicinityElement::Z_N] = NULL;
    }
    return;
  }

  if(startIdx == 5) // only possible in 3D
  {
    assert(low_right.size() == 6);
    low_right[VicinityElement::Z_P] = NULL;
    low_right[VicinityElement::Z_N] = de->lse_->nodes_[1];
  }
}

void LevelSet::SetNeighboursOfUpLeftNode(const int startIdx, const DesignElement *de)
{
  // Y_P == 0
  assert(startIdx == 3 || startIdx == 7);
  LevelSetNode *ptrnodes = de->lse_->nodes_[startIdx]; // upper left node has index 7
  assert(ptrnodes != NULL);
  ptrnodes->value = 0.0;
  
  vector<LevelSetNode*> & up_left = ptrnodes->neighbours_;
  assert(up_left.size() == ((domain->GetDim() == 2) ? 4 : 6));
  
  up_left[VicinityElement::X_P] = de->lse_->nodes_[startIdx-1];

  DesignElement *tmp_de = de->vicinity->GetNeighbour(VicinityElement::X_N);
  if(tmp_de != NULL)
  {
    up_left[VicinityElement::X_N] = tmp_de->lse_->nodes_[startIdx];
  }
  else
  {
    up_left[VicinityElement::X_N] = NULL;
  }

  up_left[VicinityElement::Y_P] = NULL;
  up_left[VicinityElement::Y_N] = de->lse_->nodes_[startIdx-3];

  if((startIdx == 3) && (domain->GetDim() == 3))
  {
    assert(up_left.size() == 6);
    up_left[VicinityElement::Z_P] = de->lse_->nodes_[7];
    tmp_de = de->vicinity->GetNeighbour(VicinityElement::Z_N);
    if(tmp_de != NULL)
    {
      up_left[VicinityElement::Z_N] = tmp_de->lse_->nodes_[3];
    }
    else
    {
      up_left[VicinityElement::Z_N] = NULL;
    }
    return;
  }

  if(startIdx == 7) // only possible in 3D
  {
    assert(up_left.size() == 6);
    up_left[VicinityElement::Z_P] = NULL;
    up_left[VicinityElement::Z_N] = de->lse_->nodes_[3];
  }
}

void LevelSet::SetNeighboursOfUpRightNode(const int startIdx, const DesignElement *de)
{
  // X_P == 0 && Y_P == 0
  assert(startIdx == 2 || startIdx == 6); // handles only those two indexes
  LevelSetNode *ptrnodes = de->lse_->nodes_[startIdx]; // upper right node
  assert(ptrnodes != NULL);
  ptrnodes->value = 0.0;
  
  vector<LevelSetNode*> &up_right = ptrnodes->neighbours_;
  assert(up_right.size() == ((domain->GetDim() == 2) ? 4 : 6));
  
  up_right[VicinityElement::X_P] = NULL; // only called when this is true
  up_right[VicinityElement::X_N] = de->lse_->nodes_[startIdx+1];
  up_right[VicinityElement::Y_P] = NULL; // only called when this is true
  up_right[VicinityElement::Y_N] = de->lse_->nodes_[startIdx-1];
  
  if((startIdx == 2) && (domain->GetDim() == 3))
  {
    assert(up_right.size() == 6);
    up_right[VicinityElement::Z_P] = de->lse_->nodes_[6];
    const DesignElement *tmp_de = de->vicinity->GetNeighbour(VicinityElement::Z_N);
    if(tmp_de != NULL)
    {
      up_right[VicinityElement::Z_N] = tmp_de->lse_->nodes_[2];
    }
    else
    {
      up_right[VicinityElement::Z_N] = NULL;
    }
    return;
  }
  
  if(startIdx == 6) // only possible in 3D
  {
    assert(up_right.size() == 6);
    up_right[VicinityElement::Z_P] = NULL;
    up_right[VicinityElement::Z_N] = de->lse_->nodes_[2];
  }
}

void LevelSet::CalcShapeGradientOnAllElements()
{
  LOG_DBG(ls) << "calculating shape gradients";
  
  // first, reset shape gradient values in the nodes
  for(unsigned int n = 0; n < nodes_.size(); ++n)
  {
    // FIXME maybe also the other node values have to be reset...
    nodes_[n].shapegrad = 0.0;
    nodes_[n].intersection_length = 0.0;
  }
  
  // disable the transfer functions to get the correct material tensors
  // for elements with a density smaller than 1.0 we need the tensors to be treated
  // as if the density was 1.0
  // this is realized via GetErsatzMaterialFactor from DesignSpace
  optimization->GetDesign()->DisableTransferFunctions();

  // first, get all element solution vectors
  StdVector<SingleVector*> &solution_vectors = optimization->getSolutionVectors();
  // we assume base = 0 -> see CalcU1KU2 in ErsatzMaterial
  const unsigned int elements(optimization->GetDesign()->GetNumberOfElements());
  assert(elements == space_.size());
  linElastInt *form = optimization->getBDBForm();
  for(unsigned int i = 0; i < elements; ++i)
  {
    // FIXME do we need the number of curr_elem at this point?
    // solution_vectors could be ordered differently than space_...
    // const unsigned int curr_elem(space_[i].de_->elem->elemNum - 1);
    // old Vector<double>& u1_vec = dynamic_cast<Vector<double>& >(*solution_vectors[i]);
    Vector<double>& u1_vec = dynamic_cast<Vector<double>& >(*solution_vectors[i]);
    // Vector<double>& u1_vec = dynamic_cast<Vector<double>& >(*solution_vectors[curr_elem]);
    space_[i].CalcShapeGrad(u1_vec, form);
  }

  // re-enable the transfer functions
  optimization->GetDesign()->EnableTransferFunctions();
}

void LevelSet::TransportLevelSet(const double dt)
{
  static const unsigned int ns(nodes_[0].neighbours_.size());
  // \phi_{ij}^{n+1} = \phi_{ij}^n - F * dt * |\grad(\phi_{ij}^n)|
  // FIXME F = shapegrad here, but should be f_ext!!
  // get the volume constraint from the optimization
  const double penalty(optimization->constraints.Get(Condition::VOLUME)->penalty);
  const double vol(optimization->constraints.Get(Condition::VOLUME)->GetBoundValue());
  LOG_DBG3(ls) << "using volume constraint of " << vol;
  
  double constraint(0.0);
  for(unsigned int n = 0; n < nodes_.size(); ++n)
  {
    constraint += nodes_[n].intersection_length;
  }
  constraint -= vol;
  constraint *= penalty;
  
  for(unsigned int n = 0; n < nodes_.size(); ++n)
  {
    LevelSetNode &node = nodes_[n];
    LOG_DBG3(ls) << "node " << n << " (cfs-node " << node.GetNumber() << ")";
    LOG_DBG3(ls) << "old levelset value = " << node.value << ", shapegrad = " << node.shapegrad;
    double update(GetNormGradientAtNode(n));
    update *= dt;
    update *= node.shapegrad + constraint;
    
    // check for nodes at the mesh boundary and truncate negative values to 0.0
    bool boundary(false);
    for(unsigned int k = 0; k < ns && !boundary; ++k)
    {
      if(node.neighbours_[k] == NULL) boundary = true;
    }
    
    if(boundary && (node.value - update) < 0.0)
      node.value = 0.0;
    else
      node.value -= update;
    
    LOG_DBG3(ls) << "new levelset value = " << node.value << ", update = " << update;
  }
}

/**************** ACTION ************************/
LevelSet::Action::Action(PtrParamNode const pn) : modulus_(1), perform_(1), first_(true)
{
  type_     = type.Parse(pn->Get("type")->As<std::string>());
  modulus_  = pn->Get("modulus")->As<Integer>();
  perform_  = pn->Get("perform")->As<Integer>();
  first_    = pn->Get("first")->As<bool>();
}

LevelSet::Action* LevelSet::IsTriggered(const Action::Type type, const int iteration)
{
  for(unsigned int i = 0, ss = action_.size(); i < ss; ++i)
  {
     Action& action = action_[i];
     if(action.GetType() == type)
     {
       if(iteration == 1)
			 {
				 if(action.GetFirst()) return &action;
				 else return NULL;
			 }
       // we trigger the action according to the modulus
       if(iteration % action.GetModulus() == 0) return &action;
     }
  }
  // not found
  return NULL;
}

} // end of namespace
