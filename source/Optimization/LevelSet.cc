#include "Optimization/LevelSet.hh"
#include "Optimization/Optimization.hh"
#include "Optimization/DesignSpace.hh"
#include "Optimization/DesignElement.hh"
#include "Domain/domain.hh"
#include "Domain/elem.hh"
#include "Domain/grid.hh"
#include "DataInOut/Logging/cfslog.hh"

using namespace std;

using namespace CoupledField;

DECLARE_LOG(ls)
DEFINE_LOG(ls, "levelSet")

/** Set in Optimization::SetEnums() */
Enum<LevelSet::Action::Type>  LevelSet::Action::type;

StdVector<LevelSetElement::DistanceHelper> LevelSetElement::distance_tmp_; 

LevelSet::LevelSet(Optimization* optimization, ParamNode* pn)
 : BaseOptimizer(optimization, pn)
{
  last_node_index_ = 0;
  design = &optimization->GetDesign()->data;

  // no (auto) scaling
  PostInit(1.0, true);  
  
  // reduce to our actual ParamNode
  pn = pn->Get(Optimization::optimizer.ToString(Optimization::LEVEL_SET));
  
  // in constrast to other optimizers we require this element for the mandatory actionList
  // there is not meaningfull default 
  if(pn == NULL) throw Exception("the element <levelSet> is mandatory for level-set optimization");
  StdVector<ParamNode*> list = pn->Get("actionList")->GetChildren();
  for(unsigned int i = 0; i < list.GetSize(); i++)
    action_.Push_back(Action(list[i]));
  
  // This is just outsourced constructor stuff
  SetupDesign();
  SetupSpace(pn->Get("subLevelSpacing")->AsDouble());
}

void LevelSet::SetupDesign()
{
  // we need the vicinity in the design elements set (if not set yet)
  // mark! assues only one design variable!!
  DesignElement& de = (*design)[0];
  if(de.vicinity == NULL) VicinityElement::Init(optimization->GetDesign());
  assert(de.vicinity != NULL);

  // create the level set nodes firs. 
  for(unsigned int i = 0; i < design->GetSize(); i++)
    (*design)[i].lsn = new LevelSetNode(&(*design)[i], -1);   // create level set element -> interior default initialization = -1

  // set the ghost nodes
  int ghosts = 0; // count the ghosts for resizing space
  for(unsigned int i = 0; i < design->GetSize(); i++)
    ghosts += (*design)[i].lsn->CreateGhostElements(1); // Set ghost element(s) -> boundary elements are 1. is stored in design element vicinity!

  assert(ghosts >= 1);
  LOG_DBG(ls) << "ghosts: " << ghosts;  
  LOG_DBG3(ls) << "LevelSet::nodes with ghosts w/o vicinity " << optimization->GetDesign()->ToString();

  // in another run, link the ghost elements -> set their vicinity information, 
  // the non-ghost elements have this already
  for(unsigned int i = 0; i < design->GetSize(); i++)
    (*design)[i].lsn->SetGhostVicinity();

  LOG_DBG3(ls) << "LevelSet::nodes with ghosts w/  vicinity " << optimization->GetDesign()->ToString();
  
  // now set set the "outer design array", not of that much use
  for(unsigned int n = 0; n < design->GetSize(); n++)  
  {
    LevelSetNode* lsn = (*design)[n].lsn;
    nodes.Push_back(lsn);
    for(unsigned int g = 0; g < lsn->GetGhosts().GetSize(); g++)
      nodes.Push_back(lsn->GetGhosts()[g]->lsn);
  }
}

void LevelSet::SetupSpace(double subLevelSpace)
{
  // create our level set element space. 
  // The number of level set elements equals the number of desing nodes.
  // e.g. 9 design elements = 9 finite elemens
  //    = 16 nodes for the finite lements
  //    = 25 level set nodes = 9 design elements + ghost elements containing a lsn
  //    = 16 level set elements
  for(unsigned int i = 0; i < design->GetSize(); i++)
  {
    // in this algorith we define a design node to be the upper right node of the
    // to be created element. If we have a +x or +y ghost neighbour we create
    // and additional element if both we also create a corner element. (see example above!)
    DesignElement* de = &(*design)[i];

    AddLevelSetElement(de, subLevelSpace);

    // is our right neighbour a ghost, then create additional element
    if(de->vicinity->GetNeighbour(VicinityElement::X_P)->lsn->IsGhost())
      AddLevelSetElement(de->vicinity->GetNeighbour(VicinityElement::X_P), subLevelSpace);

    if(de->vicinity->GetNeighbour(VicinityElement::Y_P)->lsn->IsGhost())
      AddLevelSetElement(de->vicinity->GetNeighbour(VicinityElement::Y_P), subLevelSpace);

    // in 2D there is one corner left!
    if(de->vicinity->GetNeighbour(VicinityElement::X_P)->lsn->IsGhost()
        && de->vicinity->GetNeighbour(VicinityElement::Y_P)->lsn->IsGhost()) 
      AddLevelSetElement(de->vicinity->GetNeighbour(VicinityElement::X_P, VicinityElement::Y_P ), subLevelSpace);
  }

  LOG_DBG3(ls) << "LevelSet::space " << LevelSetElement::ToString(space);  
}

void LevelSet::AddLevelSetElement(DesignElement* upper_right, double subLevelSpace)
{
  space.Push_back(LevelSetElement(subLevelSpace));
  LevelSetElement& lse = space.Last();
  assert(lse.nodes.GetSize() == 0);
  // "this node" is always part of the element (as first node)
  lse.nodes.Push_back(upper_right->lsn); // start upper right then counter clock wiese
  lse.nodes.Push_back(upper_right->vicinity->GetNeighbour(VicinityElement::X_N)->lsn);
  lse.nodes.Push_back(upper_right->vicinity->GetNeighbour(VicinityElement::X_N, VicinityElement::Y_N)->lsn);    
  lse.nodes.Push_back(upper_right->vicinity->GetNeighbour(VicinityElement::Y_N)->lsn);
  
  // now set the element node reference, we take the node coordinates of the first non-ghost node
  // and the the one with the smallest distance to the center location.
  // note, that the center location of ghost nodes is an extrapolation
  // 
  // the idea to take the common node of our lsn doesn't work straight forward as ghost elements
  // hava an extrapolated center location but no mesh nodes
  
  // find the coordinate of this element by averaging
  Point center;
  for(unsigned int i = 0; i < lse.nodes.GetSize(); i++)
    center += *(lse.nodes[i]->de_->GetLocation());
  center /= lse.nodes.GetSize();

  // find a non-ghost element - one of the mesh nodes shall coincide with the lse center
  DesignElement* de = NULL;
  for(unsigned int i = 0; i < lse.nodes.GetSize(); i++)
    if(lse.nodes[i]->IsGhost()) continue;
                      else de = lse.nodes[i]->de_; 
  
  Matrix<double>  coords;
  domain->GetGrid()->GetElemNodesCoord(coords, de->elem->connect, false );
  
  assert(domain->GetGrid()->GetDim() == 2);
  lse.nodeNumber = 0;
  for(unsigned int c = 0; lse.nodeNumber == 0 && c < coords.GetNumCols(); c++)
  {
    //LOG_DBG3(ls) << "LevelSet::AddLevelSetElement compare (" << coords[0][c] << ", " <<  coords[1][c] 
    //             << ") against (" << center.data[0] << ", " << center.data[1] << ")"; 
    if(abs(coords[0][c] - center.data[0]) > 1e-10) continue;
    if(abs(coords[1][c] - center.data[1]) > 1e-10) continue;
    // we have our candiate
    lse.nodeNumber = de->elem->connect[c];
  }
  if(lse.nodeNumber == 0) EXCEPTION("could not find mesh node for level set element based on de " << upper_right->ToString());
  
  LOG_DBG3(ls) << "LevelSet::AddLevelSetElement: center: " << center.ToString() << " node_number: " << lse.nodeNumber << " lse: " << lse.ToString();
}

LevelSetElement* LevelSet::Find(unsigned int node)
{
  // go further
  for(unsigned int i = last_node_index_; i < space.GetSize(); i++)
  {
    if(space[i].nodeNumber == node)
    {
      last_node_index_ = i;
      return &space[i];
    }
  }
  
  // go beyond
  for(unsigned int i = 0; i < last_node_index_; i++)
  {
    if(space[i].nodeNumber == node)
    {
      last_node_index_ = i;
      return &space[i];
    }
  }
  
  EXCEPTION("none of the " << space.GetSize() << " level set elements is assigned to node " << node);
}


void LevelSet::SolveProblem()
{
  // up to now we simply solve the forward problem once, like in evaluateInitialDesign
  std::cout << "Start Level-Set optimization ..." << std::endl;

  // initialize to start values, the level set node values were initialized on the construction of the elements
  // this are cheap operations, they might be done again in the first operations
  UpdateLevelSetNodes(); // gradients, normal, curvature
  
  while(!optimization->IsMinimumReached() && optimization->GetCurrentIteration() < optimization->GetMaxIterations())
  {
    int iter = optimization->GetCurrentIteration();
    optimization->SolveStateProblem();
    
    optimization->CalcObjectiveGradient(NULL);
    if(optimization->constraints.GetSize() > 0)
      optimization->CalcConstraintGradient(NULL);

    Action* action = NULL;
    if(IsTriggered(Action::SIGNED_DISTANCE_FIELD, iter) != NULL)
      CalcSignedDistanceField();

    if((action = IsTriggered(Action::FIRST_ORDER_FD, iter)) != NULL)
      for(int i = 0; i < action->perform; i++)
        EvolveFirstOrderFiniteDifferences();
    
    if((action = IsTriggered(Action::TRIVIAL_HOLE, iter)) != NULL)
      for(int i = 0; i < action->perform; i++)
        MakeTrivialHole();
    
    // effects of action
    UpdateLevelSetNodes(); // gradients, normal, curvature, design   
    
    optimization->CalcObjective();
    optimization->CommitIteration();
  }
}

void LevelSet::UpdateLevelSetNodes()
{
  for(unsigned int i = 0; i < design->GetSize(); i++)
    (*design)[i].lsn->CalcGradients();
  
  for(unsigned int i = 0; i < design->GetSize(); i++)
  {
    LevelSetNode* lsn = (*design)[i].lsn; 
    lsn->UpdateNormal();
    lsn->UpdateCurvature();
    lsn->UpdateDesign();   // "export" to the design element
  }
}


void LevelSet::CalcSignedDistanceField()
{
  // we assume the level set values to be upd to date
  // > quadratic effort! :(
  
  // We search for the distance of level-set nodes (including ghosts) to 
  // the 0-level within level set elements.
  for(unsigned outer = 0; outer < nodes.GetSize(); outer++)
  {
    LevelSetNode* lsn = nodes[outer];
    double min_dist = 1e15; 
    
    for(unsigned int inner = 0; inner < space.GetSize(); inner++)
    {
      LevelSetElement& lse = space[inner];
      if(!lse.ContainsFront()) continue;
      double distance = lse.CalcDistance(lsn);
      min_dist = min(min_dist, distance);
    }
    // the signed distance does not change its sign!
    double signed_distance = (lsn->GetValue() < 0 ? -1.0 : 1.0) * min_dist;
    lsn->SetValue(signed_distance, false); // we also update ghost elements!
  }
}

double LevelSet::FindMaxTimeStep()
{
  double max_velocity = 0.0; // we make it always positive!
  
  for(unsigned int i = 0; i < design->GetSize(); i++)
    max_velocity = max(max_velocity, abs((*design)[i].GetObjectiveGradient(DesignElement::PLAIN)));
  assert(max_velocity > 0);
  
  // we dont' call this method too often, so we retrieve the values here instead of spoiling the class attributes
  // assume a uniform grid
  double min, max;
  Matrix<double>  coords; 
  domain->GetGrid()->GetElemNodesCoord(coords, (*design)[0].elem->connect, false );
  (*design)[0].elem->ptElem->GetMaxMinEdgeLength(coords, max, min);
  
  // the CFL condition is, that we do not cross a element. Assuming the worst case direction
  double delta_t = min / max_velocity;
  LOG_DBG(ls) << "LevelSet::FindMaxTimeStep max |velocity|=" << max_velocity << " min h=" << min << " -> delta_t = " << delta_t;
  return delta_t;
}

void LevelSet::EvolveFirstOrderFiniteDifferences()
{
  double delta_t = FindMaxTimeStep();

  for(unsigned int i = 0; i < design->GetSize(); i++)
    (*design)[i].lsn->CalcNextFirstOrderHamiltonJacobiStep(delta_t);
}

void LevelSet::MakeTrivialHole()
{
  // check all nodes, and make a hole to the first element where the full neighbourhood is solid
  StdVector<DesignElement*> neighbourhood;
  for(unsigned int e = 0; e < design->GetSize(); e++)
  {
    DesignElement& de = (*design)[e]; 
    de.vicinity->GetFullNeighbourhood(neighbourhood);
    bool canditate = true;
    for(unsigned int n = 0; canditate && n < neighbourhood.GetSize(); n++)
    {
      assert(neighbourhood[n] != NULL); // there are always ghosts!
      LevelSetNode* lsn = neighbourhood[n]->lsn;
      if(!lsn->IsSolid()) canditate = false;
    }
    
    if(canditate)
    {
      de.lsn->SetHole();
      return;
    }
  }
  
  throw Exception("No material left to create another (trivial chosen) hole");
}


LevelSet::Action* LevelSet::IsTriggered(Action::Type type, int iteration)
{
  for(unsigned int i = 0; i < action_.GetSize(); i++)
  {
     Action& action = action_[i];
     if(action.type_ == type)
     {
       if(iteration % action.modulus == 0) return &action;
     }
  }
  
  // not found
  return NULL;
}

LevelSet::Action::Action(ParamNode* pn)
{
  type_    = type.Parse(pn->Get("type"));
  modulus  = pn->Get("modulus")->AsInt();
  perform  = pn->Get("perform")->AsInt();
}

LevelSetElement::LevelSetElement(double subLevelSpacing)
{
  subLevelSpacing_ = subLevelSpacing;
  nodes.Reserve(domain->GetGrid()->GetDim() == 2 ? 4 : 8);
  nodes.Init(NULL);
}

double LevelSetElement::GetCenterValue()
{
  double result = 0;
  for(unsigned int i = 0; i < nodes.GetSize(); i++)
    result += nodes[i]->GetValue();
  
  return result / nodes.GetSize();
}


bool LevelSetElement::ContainsFront()
{
  int pos = nodes[0]->GetValue() >= 0;
  
  for(unsigned int i = 1; i < nodes.GetSize(); i++)
    if((nodes[i]->GetValue() >= 0) != pos)
      return true;
  
  return false; // no sign change within our nodes
}

double LevelSetElement::CalcDistance(LevelSetNode* lsn)
{
  assert(ContainsFront());
  assert(domain->GetGrid()->GetDim() == 2);
  
  const double h = subLevelSpacing_; // this is our subdivision
  const double step = 1.0/h;
  distance_tmp_.Resize(static_cast<UInt>((h+1) * (h+1)));
  
  Point* other = lsn->de_->GetLocation();
  
  // opposite edge corners
  Point* upper_right = nodes[0]->de_->GetLocation();
  Point* lower_left  = nodes[2]->de_->GetLocation(); 
  
  double d_x = upper_right->data[0] - lower_left->data[0];
  double d_y = upper_right->data[1] - lower_left->data[1];
  
  // subdivide this element in a h*h grid
  // we could save effor by not storing the points with the same value sign as lsn.
  //
  // Concurrently we search for the closest distance.
  // Our assumption is, that we take the closest point which has another sign than the test node
  // An alternative implementation would be to take distance to the sub-node with the smallest value
  // of different sign - but the first idea seems to be more robust.
  double closest = 1e30;
  for(unsigned int x = 0; x <= h; x++)
  {
    for(unsigned int y = 0; y <= h; y++)
    {
      DistanceHelper& helper = distance_tmp_[static_cast<UInt>(x * h + y)]; 
      helper.sub_location.data[0] = lower_left->data[0] + x * step * d_x; 
      helper.sub_location.data[1] = lower_left->data[1] + y * step * d_y;
      helper.sub_location.data[2] = 0.0;
      helper.value = nodes[0]->GetValue() * (x * step) * (y * step)              // upper_right
                  +  nodes[1]->GetValue() * (1.0 - x * step) * (y * step)        // upper_left
                  +  nodes[2]->GetValue() * (1.0 - x * step) * (1.0 - y * step)  // lower_left
                  +  nodes[3]->GetValue() * (x * step) * (1.0 - y * step);       // lower_right

      LOG_DBG3(ls) << "LevelSetElement::CalcDistance sum 0: " << nodes[0]->GetValue() << "*" << (x * step) * (y * step)
                                                    << " 1: " << nodes[1]->GetValue() << "*" << (1.0 - x * step) * (y * step)
                                                    << " 2: " << nodes[2]->GetValue() << "*" << (1.0 - x * step) * (1.0 - y * step)
                                                    << " 3: " << nodes[3]->GetValue() << "*" << (x * step) * (1.0 - y * step);                                                     
      
      helper.distance = other->dist(helper.sub_location);
      
      if(!SameSign(helper.value, lsn->GetValue()))
        closest = min(closest, helper.distance);
      
      LOG_DBG3(ls) << "LevelSetElement::CalcDistance x=" << x << " y=" << y << " point=" 
                         << helper.sub_location.ToString() << " value=" << helper.value 
                         << " distance=" << helper.distance << " closest=" << closest;
    }
  }

  
  LOG_DBG2(ls) << "LevelSetElement::CalcDistance elem=" << ToString() << " other=" << lsn->de_->ToString() 
               << " level set value=" << lsn->GetValue() << " closest=" << closest; 
  assert(closest < 1e30);
  return closest;
}


/** Debug output for a single element */
std::string LevelSetElement::ToString()
{
  stringstream ss;
  ss << "node[";
  for(unsigned int i = 0; i < nodes.GetSize(); i++)
    ss << " value=" << nodes[i]->GetValue() 
       << " num=" << (nodes[i]->IsGhost() ? 0 : nodes[i]->de_->elem->elemNum);
  ss << "]";
  return ss.str();
}

/** Debug output for the complete space */
std::string LevelSetElement::ToString(StdVector<LevelSetElement>& space)
{
  stringstream ss;
  for(unsigned int i = 0; i < space.GetSize(); i++)
    ss << "#" << i << "=" << space[i].ToString();
  return ss.str();
}
