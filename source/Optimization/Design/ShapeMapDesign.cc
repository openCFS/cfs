/*
 * ShapeMapDesign.cc
 *
 *  Created on: Mar 11, 2016
 *      Author: fwein
 */
#include "Optimization/Design/ShapeMapDesign.hh"
#include "Optimization/Design/DensityFile.hh"
#include "Optimization/Excitation.hh"
#include "DataInOut/ProgramOptions.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Utils/tools.hh"
#include "Utils/mathParser/mathParser.hh"

using std::string;
using std::to_string;
using std::exp;
using std::log;

namespace CoupledField {

DEFINE_LOG(SMD, "shapeMapDesign")

StdVector<Vector<double> > ShapeMapDesign::newtonCotes = GetNewtonCotes();

ShapeMapDesign::ShapeMapDesign(StdVector<RegionIdType>& regionIds, PtrParamNode pn, ErsatzMaterial::Method method)
: FeaturedDesign(regionIds, pn, method)
{
  setup_timer_->Start();

  intStrategy_.Add(TAILORED, "tailored");

  this->combine_ = combine.Parse(pn->Get("shapeMap/overlap")->As<string>());
  this->boundary_ = boundary.Parse(pn->Get("shapeMap/shape")->As<string>());
  if(boundary_ == TANH || combine_ == TANH_SUM) {
    if(!pn->Has("shapeMap/beta"))
      throw Exception("'shapeMap' attribute 'beta' mandatory for 'shape'='tanh' or 'overlap'='tanh_sum'");
    this->beta_ = pn->Get("shapeMap/beta")->As<double>();
  }
  if(boundary_ == LINEAR && combine_ == TANH_SUM)
    info_->Get(ParamNode::HEADER)->SetWarning("'overlap'='tanh_sum' with 'shape'='linear' is questionable!");
  this->enforce_bounds_ = pn->Get("shapeMap/enforce_bounds")->As<bool>();
  this->relative_node_bound_ = pn->Get("shapeMap/relative_node_bound")->As<double>();
  this->relative_profile_bound_ = pn->Get("shapeMap/relative_profile_bound")->As<double>();
  this->export_leveset_ = pn->Has("shapeMap/export") ? pn->Get("shapeMap/export/enable")->As<bool>() : false;

  if(pn->Get("shapeMap/gradplot")->As<bool>())
    gradplot_.open((progOpts->GetSimName() + ".grad.dat").c_str()); // the auto destructor does the job.

  // set shape_, shape_param_ and map_, does not apply the mapping yet
  SetupDesign(pn->Get("shapeMap"));

  // numInt had to wait for n_, note that we give the this pointer within the constructor
  this->numInt_.Init(this, pn->Get("shapeMap"), info_->Get("shapeMap/numInt"));

  if(IsProfileFixed() && relative_profile_bound_ >= 0.0) {
    info_->Get(ParamNode::HEADER)->SetWarning("reset 'relative_profile_bound' as the profile is fixed");
    this->relative_profile_bound_ = -1.0;
  }

  // give a warning that relative bounds are set even when we have no given non-trivial initial design.
  // the the relative_bounds my be a mistake we need to warn. It must not be an error as for continuation without
  // initial design we need that setting
  if((relative_node_bound_ >= 0.0 || relative_profile_bound_ >= 0.0) && !DensityFile::NeedLoadErsatzMaterial())
    info_->Get(ParamNode::HEADER)->SetWarning("relative shape map bounds overwrite design bounds");

  // validate result settings
  if((dim_ == 3 && GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::SHAPE_MAP_GRAD, DesignElement::SM_NODE) != -1)   ||
      (dim_ == 2 && GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::SHAPE_MAP_GRAD, DesignElement::SM_NODE_A) != -1) ||
      (dim_ == 2 && GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::SHAPE_MAP_GRAD, DesignElement::SM_NODE_B) != -1))
    info_->Get(ParamNode::HEADER)->SetWarning("'result' attribute 'detail' with 'node' for 2D and 'node_a' and 'node_b' for 3D");

  // from the shapes create the rho equivalent ShapeMapVariable variables
  SetupShapeParam();

  // copy the non-fixed stuff to opt_shape_param_ and reflect symmetry
  SetupOptParam();

  // now with possibly induced shapes we may map to the design to be ready for initial evaluation
  MapFeatureToDensity();

  setup_timer_->Stop();

  LOG_DBG(SMD) << "regions: " << regionIds.ToString();
}


ShapeMapDesign::ShapeParam* ShapeMapDesign::InduceSymmetryNodeHelper(ShapeParam& ref_node)
{
  // naked as new born
  // make sure our push back does not
  shape_.Push_back(ShapeParam(), false);

  ShapeParam* ind = &(shape_.Last());
  // add to ref_node
  ref_node.induced.Push_back(ind);
  ind->induce.master = &ref_node;

  // set properties
  ind->idx = shape_.GetSize()-1;
  ind->CopyProperties(&ref_node);

  // don't forget to set further ind.induce flags
  return ind;
}

void ShapeMapDesign::InduceSymmetryNodeHelper(ShapeParam& first, ShapeParam& second)
{
  // set induced and copies data
  ShapeParam* ind_first  = InduceSymmetryNodeHelper(first);
  ShapeParam* ind_second = InduceSymmetryNodeHelper(second);

  assert(ind_first == first.induced.Last());
  assert(ind_first->dof == first.dof);
  assert(ind_first->orientation == first.orientation);
  assert(ind_first->orientation != ShapeParamElement::NOT_SET);
  assert(ind_first->induce.master == &first);
  assert(ind_first->induce.reciprocal == false); // set later for only one node

  ind_first->other_center = ind_second;
  ind_second->other_center = ind_first;
}


void ShapeMapDesign::Induce2DSymmetryNodes(ShapeParam& ref_node)
{
  LOG_DBG(SMD) << "I2SN ref=" << ref_node.ToString() << " SIMS=" << ref_node.ShallInduceMirrorSymmetry()
                   << " SICS=" << ref_node.ShallInduceCloneSymmetry() << " SIDS=" << ref_node.ShallInduceDiagonalSymmetry();

  assert(ref_node.Is2DShape());

  if(ref_node.ShallInduceMirrorSymmetry())
  {
    ShapeParam* ind = InduceSymmetryNodeHelper(ref_node);
    assert(ref_node.induced.GetSize() >= 1);
    assert(ind == ref_node.induced.Last());
    ind->induce.reciprocal = true;
    assert(ind->type == NODE);
    assert(ind->ShallInduceMirrorSymmetry());
    LOG_DBG(SMD)<< "ISM -> mirror sym " << ind->ToString();
  }

  if(ref_node.ShallInduceDiagonalSymmetry())
  {
    ShapeParam* ind = InduceSymmetryNodeHelper(ref_node);
    // not reciprocal but copied values
    ind->FlipOrientation();
    assert(ind->type == NODE);
    assert(ind->dof != ref_node.dof);
    assert(ind->orientation != ref_node.orientation);
    LOG_DBG(SMD) << "ISM -> diag sym " << ind->ToString();
  }
  // third is diagonal and then orthogonal (to the side of mapped of original)
  if(ref_node.ShallInduceDiagonalSymmetry() && ref_node.ShallInduceMirrorSymmetry())
  {
    ShapeParam* ind = InduceSymmetryNodeHelper(ref_node);
    ind->FlipOrientation();
    ind->induce.reciprocal = true; // from the mirroring part
    assert(ind->type == NODE); // only nodes are reciprocal
    assert(ind->dof != ref_node.dof);
    assert(ind->orientation != ref_node.orientation);
    LOG_DBG(SMD) << "SSD diag mirror sym " << ind->ToString();
  }
}



void ShapeMapDesign::InduceCenterSymmetryNodes(ShapeParam& first, ShapeParam& second)
{
  LOG_DBG(SMD) << "ICNSN first=" << first.ToString() << " SIMS=" << first.ShallInduceMirrorSymmetry()
                   << " SICS=" << first.ShallInduceCloneSymmetry() << " SIDS=" << first.ShallInduceDiagonalSymmetry();
  LOG_DBG(SMD) << "ICNSN second=" << second.ToString() << " SIMS=" << second.ShallInduceMirrorSymmetry()
                   << " SICS=" << second.ShallInduceCloneSymmetry() << " SIDS=" << second.ShallInduceDiagonalSymmetry();

  assert(first.IsFirstCenterNode() && second.IsSecondCenterNode());
  assert(first.GetSecondCenterNode() == &second);
  assert(first.orientation == second.orientation && first.orientation != ShapeParamElement::NOT_SET);

  // we handle square symmetry as required for 3D Bloch.
  // we simplify to either all directions sym_mirroring and diag or either aribtrary sym_mirror and diag
  // extend if there is really a use case for arbitrary combinations
  bool square_sym = first.diag && (first.x_sym || first.y_sym || first.z_sym);
  if(square_sym && !(first.x_sym && first.y_sym && first.z_sym))
    throw Exception("currently 'diagonal_sym' is only implemented with all other symmetries or none");

  // square symmetry is the following procedure
  // 1) induce the original center shape three times (mirror/clone, clone/mirror, mirror/mirror)
  // 2) first diagonal induced shape
  // 3) mirror induce the first diagonal shape three times
  // 4) second diagonal induced shape
  // 5) mirror induce the second diagonal shape three times

  // without diag we have the combinations clone/mirror, mirror/clone and mirror/mirror.
  // for this two necessary ShallInduceMirrorSymmetry() need to be set. The third is mapping
  if(first.ShallInduceMirrorSymmetry() && second.ShallInduceCloneSymmetry())
  {
    InduceSymmetryNodeHelper(first, second);
    // now the mirror node becomes reciprocal
    assert(first.ShallInduceMirrorSymmetry());
    assert(first.induced.GetSize() >= 1);
    first.induced.Last()->induce.reciprocal = true; // the other of false by default

    LOG_DBG(SMD) << "ICNSN ind_first=" << first.induced.Last()->ToString();
    LOG_DBG(SMD) << "ICNSN ind_second=" << second.induced.Last()->ToString();
  }

  if(first.ShallInduceCloneSymmetry() && second.ShallInduceMirrorSymmetry())
  {
    InduceSymmetryNodeHelper(first, second);
    second.induced.Last()->induce.reciprocal = true; // the other of false by default

    LOG_DBG(SMD) << "ICNSN ind_first=" << first.induced.Last()->ToString() << " ind_second=" << second.induced.Last()->ToString();
  }

  // now the mirror/mirror variant
  if(first.ShallInduceMirrorSymmetry() && second.ShallInduceMirrorSymmetry())
  {
    InduceSymmetryNodeHelper(first, second);
    first.induced.Last()->induce.reciprocal = true;
    second.induced.Last()->induce.reciprocal = true;

    LOG_DBG(SMD) << "ICNSN ind_first=" << first.induced.Last()->ToString() << " ind_second=" << second.induced.Last()->ToString();
  }

  // we known only a common 3D center nodes diagonal mirroring which always flips twice.
  // out of a shape a cross is induced
  assert((first.ShallInduceDiagonalSymmetry() &&  second.ShallInduceDiagonalSymmetry())
      || (!first.ShallInduceDiagonalSymmetry() && !second.ShallInduceDiagonalSymmetry()));
  if(first.ShallInduceDiagonalSymmetry())
  {
    // we induce two shapes, where one of the dofs is replaced by orientation
    // assume a horizontal shape with first.dof=y and second.dof=z and orientation=x
    // we first induce a vertial shape with first.dof=x, second.dof=z and orientation=y
    InduceSymmetryNodeHelper(first, second);
    first.induced.Last()->FlipOrientation(0); // flips also second, do NOT flip second manually, it will flip second and first back!!

    assert(second.induced.Last()->orientation == first.induced.Last()->orientation);
    assert(first.induced.Last()->orientation != first.orientation);
    LOG_DBG(SMD) << "ICNSN diag1 ind_first=" << first.induced.Last()->ToString();
    LOG_DBG(SMD) << "ICNSN diag1 ind_second=" << second.induced.Last()->ToString();

    if(square_sym)
    {
      assert(first.type == NODE);
      InduceSymmetryNodeHelper(first, second);
      first.induced.Last()->FlipOrientation(0);
      first.induced.Last()->induce.reciprocal = true;
      LOG_DBG(SMD) << "ICNSN diag1 square1 ind_first=" << first.induced.Last()->ToString();
      LOG_DBG(SMD) << "ICNSN diag1 square1 ind_second=" << second.induced.Last()->ToString();

      InduceSymmetryNodeHelper(first, second);
      LOG_DBG(SMD) << "ICNSN induced=" << first.induced.Last()->ToString();
      first.induced.Last()->FlipOrientation(0);
      second.induced.Last()->induce.reciprocal = true;
      LOG_DBG(SMD) << "ICNSN diag1 square2 ind_first=" << first.induced.Last()->ToString();
      LOG_DBG(SMD) << "ICNSN diag1 square2 ind_second=" << second.induced.Last()->ToString();

      InduceSymmetryNodeHelper(first, second);
      first.induced.Last()->FlipOrientation(0);
      first.induced.Last()->induce.reciprocal = true;
      second.induced.Last()->induce.reciprocal = true;
      LOG_DBG(SMD) << "ICNSN diag1 square3 ind_first=" << first.induced.Last()->ToString();
      LOG_DBG(SMD) << "ICNSN diag1 square3 ind_second=" << second.induced.Last()->ToString();
    }

    // now the lateral induced shape with first.dof=y, second.dof=x and orientation=z
    InduceSymmetryNodeHelper(first, second);
    first.induced.Last()->FlipOrientation(1); // again, flip only first, it will flip second, too

    LOG_DBG(SMD) << "ICNSN diag2 ind_first=" << first.induced.Last()->ToString();
    LOG_DBG(SMD) << "ICNSN diag2 ind_second=" << second.induced.Last()->ToString();

    if(square_sym)
    {
      InduceSymmetryNodeHelper(first, second);
      first.induced.Last()->FlipOrientation(1);
      first.induced.Last()->induce.reciprocal = true;
      LOG_DBG(SMD) << "ICNSN diag2 square1 ind_first=" << first.induced.Last()->ToString();
      LOG_DBG(SMD) << "ICNSN diag2 square1 ind_second=" << second.induced.Last()->ToString();

      InduceSymmetryNodeHelper(first, second);
      first.induced.Last()->FlipOrientation(1);
      second.induced.Last()->induce.reciprocal = true;
      LOG_DBG(SMD) << "ICNSN diag2 square2 ind_first=" << first.induced.Last()->ToString();
      LOG_DBG(SMD) << "ICNSN diag2 square2 ind_second=" << second.induced.Last()->ToString();

      InduceSymmetryNodeHelper(first, second);
      first.induced.Last()->FlipOrientation(1);
      first.induced.Last()->induce.reciprocal = true;
      second.induced.Last()->induce.reciprocal = true;
      LOG_DBG(SMD) << "ICNSN diag2 square3 ind_first=" << first.induced.Last()->ToString();
      LOG_DBG(SMD) << "ICNSN diag2 square3 ind_second=" << second.induced.Last()->ToString();
    }
  }
}

/* <shapeMap beta="2">
    <center>
      <node lower=".499" initial="0.5" upper=".501" dof="x" />
      <node lower=".499" initial="0.5" upper=".501" dof="y"/>
    </center>
    <node dof="x" lower="0" upper=".5" initial=".25"/>
    <profile lower=".01" upper=".1" initial=".1"/>
  </shapeMap> */
void ShapeMapDesign::SetupDesign(PtrParamNode pn)
{
  // check rhos which should be already be set in DesignSpace::data
  assert(GetRegionIds().GetSize() == 1); // more is not implemented yet
  StdVector<int> elem_to_idx;
  StdVector<int> idx_to_elem;
  // n_ = domain->GetGrid()->GetBoundaries(GetRegionIds().First());
  n_ = SetupLexicographicMesh(domain->GetGrid(), GetRegionIds().First(), elem_to_idx, idx_to_elem);
  nx_ = n_[0];
  ny_ = n_[1];
  nz_ = n_[2];
  assert(data.GetSize() == nx_ * ny_ * nz_); // DesignSpace::data has an element for each FEM-Cell
  assert(!(dim_ == 2 && nz_ != 1));

  // read shapeParam to shape_
  ParamNodeList nodes = pn->GetList("node"); // there must be at least one for 2D
  ParamNodeList centers = pn->GetList("center"); // only 3D
  PtrParamNode  profile = pn->Get("profile"); // there must be one

  if(dim_ == 2 && nodes.IsEmpty())
    throw Exception("for 2D at least one 'node' entry for 'shapeMap' is mandatory");
  if(dim_ == 3 && nodes.IsEmpty() && centers.IsEmpty())
    throw Exception("provide at least a 'center' or 'node' entry for 'shapeMap'");

  // first we add the shapes, as we might induce additional shapes during ParseAndInit we add the profiles later. We must not resize during push_back, therefore check the initial guess
  unsigned int estimate = 0;
  if(dim_ == 2)
    estimate = 2 * nodes.GetSize() * 4; // 2 for profile and 4 for symmetry
  else
    estimate = (2 * nodes.GetSize() + 3 * centers.GetSize()) * 12; // 2 for profile, 3 for two nodes and one profile and 12(!) for square symmetry
  shape_.Reserve(estimate);

  // 3d only: centers exist only in xml and have two node childs
  for(unsigned int ci = 0; ci < centers.GetSize(); ci++)
  {
    // note that base is never stored itself but only its sub nodes to be found by FindCenters()
    ShapeParam base;
    base.idx = -1;
    base.ParseAndInit(centers[ci], NULL);
    assert(base.type == CENTER);

    ParamNodeList cn = centers[ci]->GetList("node");
    assert(cn.GetSize() == 2);

    shape_.Push_back(ShapeParam());
    ShapeParam& first = shape_.Last();
    first.idx = shape_.GetSize()-1;
    // we cannot set everything yet - e.g. the orientation, other_center, ... can only be determined after second
    first.ParseAndInit(cn[0], &base);

    assert(first.dof != ShapeParamElement::NOT_SET);
    assert(first.orientation == ShapeParamElement::NOT_SET);
    assert(first.other_center == NULL);
    assert(first.partner == NULL);

    shape_.Push_back(ShapeParam());
    ShapeParam& second = shape_.Last();
    second.idx = shape_.GetSize()-1;
    second.ParseAndInit(cn[1], &base);

    if(first.dof == second.dof)
      throw Exception("the 'nodes' of a 'center' 'ShapeMap' must not have the same dof");
    if(first.slave && second.slave)
      throw Exception("only one of the 'center' nodes can be slave (shares the same variable");
    // copy slave data like bounds, ... but don't overwrite slave and dof
    if(first.slave)
      first.CopyProperties(&second, true);
    if(second.slave)
      second.CopyProperties(&first, true);

    first.other_center = &second;
    second.other_center = &first;

    // only now we can now about the orientation of a center node!
    base.orientation = Flip(first.dof, second.dof);
    first.orientation = base.orientation;
    second.orientation = base.orientation;

    // only now, we can induce the nodes
    InduceCenterSymmetryNodes(first, second);
  } // end of looping center param nodes

  // nodes are mandatory in 2D and in 3D surfaces beside the center rods
  for(unsigned int i = 0; i < nodes.GetSize(); i++)
  {
    shape_.Push_back(ShapeParam());
    ShapeParam& item = shape_.Last(); // the item to be processed. Capacity needs to be large enough such that this reference is not fucked up
    item.idx = shape_.GetSize()-1;
    item.ParseAndInit(nodes[i], NULL); // size might be != i when we induce
    LOG_DBG(SMD) << "SSD " << item.ToString() << " : i=" << i  << " osi=" << item.ShallInduceMirrorSymmetry()
                      << " dsi=" << item.ShallInduceDiagonalSymmetry() << " ind=" << item.IsInduced();

    // first added is orthogonal
    Induce2DSymmetryNodes(item);
  }
  assert(shape_.GetSize() <= estimate);
  num_node_shapes_ = shape_.GetSize();

  // now add the profiles where 3D center nodes share on profile!
  for(unsigned int i = 0; i < (unsigned int) num_node_shapes_; i++)
  {
    ShapeParam* nodal = &shape_[i];
    if(nodal->IsSecondCenterNode())
    {
      // we are second center node and the first center node added the profile before!
      assert(nodal->GetFirstCenterNode() != NULL);
      assert(nodal->GetFirstCenterNode() != nodal);
      assert(nodal->GetFirstCenterNode()->partner != NULL && nodal->GetFirstCenterNode()->partner->type == PROFILE);
      assert(nodal == nodal->GetSecondCenterNode());
      // link the profile added before by the first center node to the second center node as partner
      nodal->partner = nodal->GetFirstCenterNode()->partner;
    }
    else
    {
      // first center node case or any other 2D or 3D case with a 1:1 node to profile link
      shape_.Push_back(ShapeParam());
      ShapeParam* prof = &shape_.Last(); // don't shadow the profile ParamNode
      prof->idx = shape_.GetSize()-1;
      nodal->partner = prof;
      prof->partner = nodal;
      prof->ParseAndInit(profile, nodal); // one profile for each node shape, give reverence to node to copy sym and orientation
      prof->dof = ShapeParamElement::NOT_SET;
      assert(prof->other_center == NULL); // we don't have this in profile

      // handle the induced stuff
      prof->induce.reciprocal = false; // only nodes can be reciprocal
      assert(prof->induce.master == NULL); // not set yet
      if(nodal->IsInduced())
      {
        assert(nodal->induce.master != NULL);
        assert(nodal->induce.master->partner != NULL);
        assert(nodal->induce.master->partner->type == PROFILE);
        assert(nodal->induce.master->partner->idx >= 1); // already set
        prof->induce.master = nodal->induce.master->partner;
        prof->induce.master->induced.Push_back(prof);
        assert(prof->induced.IsEmpty()); // induced shapes have this vector empty
      }
    }
  }
  assert(shape_.GetSize() <= estimate);
  assert((int) shape_.GetSize() <= 2 * num_node_shapes_);

  for(unsigned int i = 0; i < shape_.GetSize(); i++)
    LOG_DBG(SMD) << "SSD: " << shape_[i].ToString();
}

void ShapeMapDesign::SetupShapeParam()
{
  // we must not Resize() shape_param_!
  // for 3D with surfaces instead of rods we estimate to low
  unsigned int max_size = n_.Max();
  LOG_DBG(SMD) << "SSP: max_size=" << max_size << " from " << n_.ToString();
  assert(max_size >= nx_ && max_size >= ny_ && max_size >= nz_);
  assert(max_size == nx_ || max_size == ny_ || max_size == nz_);

  shape_param_.Reserve(shape_.GetSize() * (max_size+1));

  // take the shapes in the order they are stored in shape_ as read from xml plus induced shapes
  for(int s = 0; s < num_node_shapes_; s++)
  {
    ShapeParam& node = shape_[s];
    node.start_param = shape_param_.GetSize();
    assert(node.idx == s);

    // up to now we do not support surfaces like dof=x in 3D which means (ny_+1)*(nz_+1) variables in the yz plane
    assert(node.dof >= ShapeParamElement::X && node.dof <= ShapeParamElement::Z);
    assert(node.orientation >= ShapeParamElement::X && node.orientation <= ShapeParamElement::Z);
    assert(node.orientation != node.dof);
    assert(ShapeParamElement::X == 0 && ShapeParamElement::Z == 2); // we assume this to index n_[]

    // when in 2D dof=x then we traverse the ny_+1, with 3D center pairs x and y with traverse nz_+1 for both, with 3D surface x it will be (ny_+1)*(nz_+1)
    unsigned int end = (int) n_[node.orientation] + 1;
    for (unsigned int e = 0; e < end; e++)
      CreateShapeVariable(&node, node.orientation, e, e == 0 || e == (end - 1)); // makes a push_back to shape_param_
    node.end_param = shape_param_.GetSize();
    LOG_DBG(SMD) << "SSP: node s=" << s << " end = " << end << " shape=" << node.ToString();
  }
  num_node_shape_params_ = shape_param_.GetSize();

  // add the profiles
  for(int n = 0; n < num_node_shapes_; n++)
  {
    ShapeParam& node = shape_[n];
    // two 3D center nodes share one profile. Skip the second
    if(node.IsSecondCenterNode())
      continue;
    ShapeParam& prof = *(node.partner);
    assert(dim_ == 3 || prof.idx == num_node_shapes_ + n);
    assert(node.type == NODE && prof.type == PROFILE);
    prof.start_param = shape_param_.GetSize();
    for(int e = 0; e < node.end_param - node.start_param; e++)
      CreateShapeVariable(&prof, prof.orientation, e, e == 0 || e == (node.end_param - node.start_param - 1));
    prof.end_param = shape_param_.GetSize();
    LOG_DBG(SMD) << "SSP: prof n=" << n << " end = " << (node.end_param - node.start_param) << " nsp=" << num_node_shape_params_ << " prof=" << prof.ToString();
    assert(prof.end_param - prof.start_param == node.end_param - node.start_param); // same size
    assert(prof.start_param >= num_node_shape_params_);
    // note that the profile 0 connects to node 0 but profile 1 connects to node 2 as the nodes order is a0,b0,a1,b1,...
  }
  assert(dim_ == 3 || (int ) shape_param_.GetSize() == 2 * num_node_shape_params_); // doubles variables for 2D. In 3D different for center nodes

  // set map_ to map from shape_param to DesignSpace::data, fill with shape_param_
  map_.Resize(data.GetSize());
  StdVector<Elem*> designElems;
  domain->GetGrid()->GetElems(designElems, GetRegionIds().First()); // FIXME assumes elements in designElems are ordered!
  assert(designElems.GetSize() <= nx_ * ny_ * nz_);
  assert(map_.GetSize() == designElems.GetSize());
  assert(num_node_shapes_ > 0);

  // num_nodes has no 3D second center nodes and is then smaller num_node_shapes_
  int num_nodes = num_node_shapes_ - FindCenters().GetSize();
  LOG_DBG(SMD) << "SSP nns=" << num_node_shapes_ << " fc=" << FindCenters().GetSize() << " -> nn=" << num_nodes;

  for(unsigned int i = 0, n = map_.GetSize(); i < n; i++)
  {
    map_[i].elemval = &(data[Find(designElems[i]->elemNum)]); // is very fast and gives a layer for arbitrary element ordering in the mesh
    // each design node connects to two density elements, also for 3D center node rods
    // this comes from the bilinear interpolation.
    map_[i].nodes.Resize(num_nodes);
    map_[i].min_corner_value.Resize(num_nodes);
    map_[i].max_corner_value.Resize(num_nodes);
    for(int s = 0; s < num_nodes; s++)
      map_[i].nodes[s].Reserve(dim_ == 2 ? 2 : 4);
  }

  // setup  coord_*_ stuff
  Matrix<double>    coords;   // within the element coordinates we perform the integration
  domain->GetGrid()->GetElemNodesCoord(coords, map_.First().elemval->elem->connect, false); // no deformed mesh
  coords.GetColMin(coord_min_);
  LOG_DBG(SMD) << "SSP data=" << " min coords=" << coords.ToString();
  LOG_DBG(SMD) << "SSP data=" << " min=" << coord_min_.ToString();
  domain->GetGrid()->GetElemNodesCoord(coords, map_.Last().elemval->elem->connect, false); // no deformed mesh
  coords.GetColMax(coord_max_);
  LOG_DBG(SMD) << "SSP data=" << " max coords=" << coords.ToString();
  LOG_DBG(SMD) << "SSP data=" << " max=" << coord_max_.ToString();

  assert              (coord_min_[0] <= coords[0][0] && coord_min_[0] <= coords[0][1] && coord_min_[0] <= coords[0][2]  && coord_min_[0] <= coords[0][3]);
  assert              (coord_min_[1] <= coords[1][0] && coord_min_[1] <= coords[1][1] && coord_min_[1] <= coords[1][2]  && coord_min_[1] <= coords[1][3]);
  assert(dim_ == 2 || (coord_min_[2] <= coords[2][0] && coord_min_[2] <= coords[2][1] && coord_min_[2] <= coords[2][2]  && coord_min_[2] <= coords[2][3]));
  assert(dim_ == 2 || (coord_min_[0] <= coords[0][4] && coord_min_[0] <= coords[0][5] && coord_min_[0] <= coords[0][6]  && coord_min_[0] <= coords[0][7]));
  assert(dim_ == 2 || (coord_min_[1] <= coords[1][4] && coord_min_[1] <= coords[1][5] && coord_min_[1] <= coords[1][6]  && coord_min_[1] <= coords[1][7]));
  assert(dim_ == 2 || (coord_min_[2] <= coords[2][4] && coord_min_[2] <= coords[2][5] && coord_min_[2] <= coords[2][6]  && coord_min_[2] <= coords[2][7]));

  assert(              coord_max_[0] >= coords[0][0] && coord_max_[0] >= coords[0][1] && coord_max_[0] >= coords[0][2]  && coord_max_[0] >= coords[0][3]);
  assert(              coord_max_[1] >= coords[1][0] && coord_max_[1] >= coords[1][1] && coord_max_[1] >= coords[1][2]  && coord_max_[1] >= coords[1][3]);
  assert(dim_ == 2 || (coord_max_[2] >= coords[2][0] && coord_max_[2] >= coords[2][1] && coord_max_[2] >= coords[2][2]  && coord_max_[2] >= coords[2][3]));
  assert(dim_ == 2 || (coord_max_[0] >= coords[0][4] && coord_max_[0] >= coords[0][5] && coord_max_[0] >= coords[0][6]  && coord_max_[0] >= coords[0][7]));
  assert(dim_ == 2 || (coord_max_[1] >= coords[1][4] && coord_max_[1] >= coords[1][5] && coord_max_[1] >= coords[1][6]  && coord_max_[1] >= coords[1][7]));
  assert(dim_ == 2 || (coord_max_[2] >= coords[2][4] && coord_max_[2] >= coords[2][5] && coord_max_[2] >= coords[2][6]  && coord_max_[2] >= coords[2][7]));

  coord_step_.Resize(dim_);
  for(unsigned int i = 0; i < dim_; i++)
    coord_step_[i] = (coord_max_[i] - coord_min_[i]) / (double) n_[i];

  LOG_DBG(SMD) << "SSP data=" << data.GetSize() << " map=" << map_.GetSize() << " n_=" << n_.ToString();

  // now set the nodes idx in map_::nodes. For every rho we store here the two (2D) or 4 (3D surface) ShapeMapVariables which have an bilinear interpolation
  // within the node per shape

  if(dim_ == 2)
  {
    // Assume 2D with dof=x. Then for every y there is a x-variable
    // The mapping is now such, that
    // x0 variable effects all rho from the bottom horizontal row.
    // x1 affects the bottom row (idx=0) and the one above (idx=1).
    // x2 affects the rows with idx=1 and idx=2,
    // x3 affects the upper row only.
    // -----------------
    // |    x3         |
    // |     x2        |
    // |      x1       |
    // |     x0        |
    // -----------------

    for(int i = 0; i < num_node_shape_params_; i++)
    {
      ShapeMapVariable* spe = dynamic_cast<ShapeMapVariable*>(shape_param_[i]);
      ShapeParam* shape = FindShape(spe); // GetShape() is not initialized yet
      assert(shape->IsPart(spe));
      assert(shape->idx >= 0);

      LOG_DBG2(SMD) << "SSP: spe=" << spe->ToString() << " idx=" << spe->idx.ToString();

      switch(spe->dof_)
      {
      case ShapeParamElement::X:
      {
        int y = spe->idx[1];
        assert(y >= 0);

        // exclude the case where the y node is at the upper boundary (one node more than elements!)
        for(unsigned int x = 0; x < nx_; x++){
          LOG_DBG2(SMD) << "SSP: x=" << x << ", y=" << y << " -> " << DensityIdx(x, y);
          if(y < (int) ny_) // are we the topmost node ontop of the last element
            map_[DensityIdx(x, y)].nodes[shape->idx].Push_back(spe);
          if(y-1 >= (int) 0) // node 2 participates to the lower element (2-1) and the upper element (2)
            map_[DensityIdx(x, y-1)].nodes[shape->idx].Push_back(spe);
        }
        break;
      }
      case ShapeParamElement::Y:
      {
        int x = spe->idx[0];
        for(unsigned y = 0; y < ny_; y++)  {
          LOG_DBG2(SMD) << "SSP: x=" << x << (x < (int) nx_ ? "+" : "!") << ", y=" << y << " -> " << DensityIdx(x, y);
          if(x < (int) nx_)
            map_[DensityIdx(x, y)].nodes[shape->idx].Push_back(spe);
          if(x-1 >= 0)
            map_[DensityIdx(x-1, y)].nodes[shape->idx].Push_back(spe);
        }
        break;
      }
      default:
        assert(false);
      } // end switch
    } // end num_node_shape_params_ loop
  } // end dim == 2 case

  if(dim_ == 3)
  {
    for(int i = 0; i < num_node_shape_params_; i++) // traverse the ShapeMapVariables an set it to the proper map_::nodes
    {
      ShapeMapVariable* spe1 = dynamic_cast<ShapeMapVariable*>(shape_param_[i]); // will be first node param
      ShapeParam* shape = FindShape(spe1); // GetShape() not yet ready
      assert(spe1->GetType() == Convert(shape->type));

      // the structure idx is different from the shape idx and counts the structures within Item::nodes.
      // In 3D with center nodes we store only the fist center node shapes 0 and 2 for but the struct_idx is 0 and 1
      assert((int) FindCenters().GetSize() * 2 == num_node_shapes_);
      unsigned int struct_idx = shape->idx/2; // change when we have not only center nodes, e.g. to attribute

      // the 3D surface case is not yet implemented
      assert(shape->IsFirstCenterNode() || shape->IsSecondCenterNode());

      LOG_DBG2(SMD) << "SSP: spe1=" << spe1->ToString() << " idx=" << spe1->idx.ToString() << " s=" << shape->ToString();

      // we operate only on the first center node and search for the second
      if(shape->IsFirstCenterNode())
      {
        assert((int) shape->orientation < (int) spe1->idx.GetSize());
        assert((int) shape->orientation >= 0);
        assert((int) shape->orientation < (int) dim_);
        assert(shape->orientation != shape->dof && shape->orientation != shape->GetSecondCenterNode()->dof);
        assert(shape->GetFirstCenterNode()->orientation == shape->GetSecondCenterNode()->orientation);

        ShapeMapVariable* spe2 = GetSecondCenterNodeParam(shape, spe1);

        // in 2D a dof=x node connects to 1 or 2 y lines.
        // for a 3D rod with we have always dof pairs and they hold for each of the common orientation
        switch(shape->orientation)
        {
        case ShapeParamElement::X:
        {
          assert((spe1->dof_ == ShapeParamElement::Y || spe1->dof_ == ShapeParamElement::Z));
          assert((spe2->dof_ == ShapeParamElement::Y || spe2->dof_ == ShapeParamElement::Z));
          assert(spe1->dof_ != spe2->dof_);

          int x = spe1->idx[0];
          assert(spe1->idx[0] == spe2->idx[0]);

          for(unsigned int z = 0; z < nz_; z++)
          {
            for(unsigned int y = 0; y < ny_; y++)
            {
              LOG_DBG2(SMD) << "SSP: x=" << x << (x < (int) nx_ ? "+" : "!") << ", y=" << y << " z=" << z << " -> " << DensityIdx(x, y);
              if(x < (int) nx_)  // not yet topmost plane
                map_[DensityIdx(x, y, z)].nodes[struct_idx].Push_back(spe1, spe2);
              if(x-1 >= 0)  // also not the lowest plane
                map_[DensityIdx(x-1, y, z)].nodes[struct_idx].Push_back(spe1, spe2);
            }
          }
          break;
        }

        case ShapeParamElement::Y:
        {
          int y = spe1->idx[1];

          for(unsigned int z = 0; z < nz_; z++)
          {
            for(unsigned int x = 0; x < nx_; x++)
            {
              if(y < (int) ny_)
                map_[DensityIdx(x, y, z)].nodes[struct_idx].Push_back(spe1, spe2);
              if(y-1 >= 0)
                map_[DensityIdx(x, y-1, z)].nodes[struct_idx].Push_back(spe1, spe2);
            }
          }
          break;
        }

        case ShapeParamElement::Z:
        {
          int z = spe1->idx[2];

          for(unsigned int y = 0; y < ny_; y++)
          {
            for(unsigned int x = 0; x < nx_; x++)
            {
              if(z < (int) nz_)
                map_[DensityIdx(x, y, z)].nodes[struct_idx].Push_back(spe1, spe2);
              if(z-1 >= 0)
                map_[DensityIdx(x, y, z-1)].nodes[struct_idx].Push_back(spe1, spe2);
            }
          }
          break;
        } // end of case
        case ShapeParamElement::NOT_SET:
          assert(false);
          break;
        } // end of switch
      } // first center node case
    } // end node loop
  } // end dim == 3 case

  // might make it a little smaller
  shape_param_map_.Resize(shape_param_.GetSize());
  for(unsigned int s = 0; s < shape_.GetSize(); s++)
    for(int p = shape_[s].start_param; p < shape_[s].end_param; p++)
      shape_param_map_[p] = &shape_[s];

  assert(!shape_param_map_.Contains(NULL));
}


void ShapeMapDesign::SetupOptParam()
{
  // copy the non-fixed stuff opt_shape_param_
  // TODO one could do a better assumption using symmetries
  opt_shape_param_.Reserve(IsProfileFixed() ? num_node_shape_params_ : shape_param_.GetSize());
  for(unsigned int s = 0; s < shape_.GetSize(); s++)
  {
    ShapeParam& shape = shape_[s];
    LOG_DBG(SMD) << "SOSP opt " << shape.ToString() << " start_param=" << " center=" << shape.IsCenterShape()
                     << " first=" << shape.IsFirstCenterNode() << " second=" << shape.IsSecondCenterNode();
    LOG_DBG(SMD) << "SOSP start_param=" << shape.start_param << " end_param=" << shape.end_param
        << " start_opt=" << shape.start_opt << " end_opt=" << shape.end_opt << " ind=" << shape.IsInduced();
    // fixed and symmetry induced shapes don't belong to opt_shape_param_ by definition
    if(!shape.fixed && !shape.IsInduced() && !shape.slave)
    {
      // we copy the element when we have symmetry but not the above case (0->6, 1->5, 2->4, ...)
      bool sym_map = shape.ShallMapHalfShape();

      // of the center nodes shares its variable with the other to achieve cubical symmetry
      bool slave = shape.ShallBeMasterOfSlave(); // we are not slave, but is the other center node one?
      assert(!slave || shape.other_center != NULL);

      // see SetupShapeDesign()
      LOG_DBG(SMD)<< "SOSP sym_map=" << sym_map << " slave=" << slave;
      // add one in the mirror case with odd design elements for not mirrored center element - note 0-based counting!
      bool odd_elements = (shape.end_param - shape.start_param) % 2 == 0 ? false : true;
      shape.start_opt = opt_shape_param_.GetSize();
      unsigned int end = sym_map ? shape.start_param + (shape.end_param - shape.start_param) / 2 + (odd_elements ? 1 : 0) : shape.end_param;
      LOG_DBG(SMD)<< "SOSP opt cand odd=" << odd_elements << " sym_map=" << sym_map << " start_opt=" << shape.start_opt << " end_opt=" << shape.end_opt;

      // traverse all elements of the master shape only.
      // For the mapping case only half and then the center element is mapped or not (depending on odd or even shape elements)
      for(unsigned int p = shape.start_param; p < end; p++)
      {
        ShapeMapVariable* opt = dynamic_cast<ShapeMapVariable*>(shape_param_[p]);
        opt_shape_param_.Push_back(opt); // when we have fixed nodes we shall handle the index!!
        opt->sym = new ShapeMapDesign::ElementSymmetry(opt);
        assert(opt->sym->hidden.IsEmpty());
        ElementSymmetry* sym = opt->sym; // and this corresponding structure has all additional virtual elements
        unsigned int idx = p - shape.start_param;
        bool odd_center_elem = odd_elements && p == end - 1; // is this for an odd number the center element?

        // prior to the induced shape stuff
        assert(!(slave && shape.other_center == NULL));
        assert(!(slave && !shape.other_center->slave));
        if(slave)
          sym->AddSymmetryReference(shape_param_[shape.other_center->start_param + idx], shape.other_center, false); // not reciprocal

        // if we simply map we go up the proper center element
        if(sym_map && !odd_center_elem)
        {
          // don't add symmetry for the center element of an odd row
          // take the element from the end
          sym->AddSymmetryReference(shape_param_[shape.end_param - 1 - idx], &shape, false); // not reciprocal
          if(slave)
            sym->AddSymmetryReference(shape_param_[shape.other_center->end_param - 1 - idx], shape.other_center, false); // not reciprocal
        }

        // now we loop all induced shapes. In 2D this is for the square symmetry case parallel, diagonal and diagonal+parallel
        // we go from start to end and then in the mapping case also for the mapped element
        for(unsigned int is = 0; is < shape.induced.GetSize(); is++)
        {
          ShapeParam* induced = shape.induced[is];
          assert(!slave || induced->other_center != NULL);
          assert(!(slave && induced->other_center->induce.master == NULL));
          // only parallel (included diagonal+parallel) node shapes are reciprocal
          assert(!(induced->induce.reciprocal && induced->type != NODE)); // only mirror nodes can be reciprocal
          sym->AddSymmetryReference(shape_param_[induced->start_param + idx], induced, induced->induce.reciprocal);
          if(slave)
            sym->AddSymmetryReference(shape_param_[induced->other_center->start_param + idx], induced->other_center, induced->other_center->induce.reciprocal);

          // for symmetry mapping we come from the back but for odd elements (e.g. 7 = 0...6) we don't mat the center element (3) to itself
          if(sym_map && !odd_center_elem)
          {
            sym->AddSymmetryReference(shape_param_[induced->end_param - 1 - idx], induced, induced->induce.reciprocal);
            if(slave)
              sym->AddSymmetryReference(shape_param_[induced->other_center->end_param - 1 - idx], induced->other_center, induced->other_center->induce.reciprocal);
          }
        }

        // apply the value to the symmetry stuff, especially for the induced shape but also slaves!
        sym->ApplyDesign();
        LOG_DBG(SMD)<< "SOSP opt shape=" << shape.idx << " type=" << shape.type << " el=" << opt->GetIndex()
                 << " #oel=" << opt_shape_param_.GetSize() << sym->ToString();
      }
      shape.end_opt = opt_shape_param_.GetSize(); // luckily the complex counting of symmetry references is not necessary

      LOG_DBG(SMD)<< "SOSP shape=" << shape.idx << " type=" << shape.type << " start=" << shape.start_param << " end=" << shape.end_param << " this end=" << end
          << " start_opt=" << shape.start_opt << " end_opt=" << shape.end_opt << " odd=" << odd_elements << " sym_map=" << sym_map;
    }
  }
  // set for the design subject to optimization the opt_index_ as BaseDesignElement::index_ may not reflect the design space seen
  // by external optimizers when we have symmetry -> necessary for the sparse Jacobian.
  for(unsigned int i = 0, n = opt_shape_param_.GetSize(); i < n; i++)
    opt_shape_param_[i]->SetOptIndex(i);

  LOG_DBG(SMD)<< "SOSP osp=" << opt_shape_param_.GetSize() << " sp=" << shape_param_.GetSize() << " data=" << data.GetSize();
}

void ShapeMapDesign::SetupVirtualShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& vem, Function::Local::Locality locality)
{
  assert(f != NULL);
  assert(f->IsLocal(f->GetType()));
  // we shall be called by Local::PostInit() therefore local shall exist
  assert(f->GetLocal() != NULL);

  // we assume fixed only for profile
  if(f->GetDesignType() == BaseDesignElement::PROFILE && IsProfileFixed())
    throw Exception("Configuration error: cannot have local constraint of shape map design 'profile' when this design is fixed.");

  assert(locality == Function::Local::NEXT || locality == Function::Local::PREV_NEXT_AND_REVERSE || locality == Function::Local::PREV_NEXT || locality == Function::Local::NEXT_AND_REVERSE || locality == Function::Local::PREV_NEXT_AND_REVERSE);

  // a lot copy&paste from Function::SetupVirtualElementMap()
  bool prev = locality == Function::Local::PREV_NEXT_AND_REVERSE || locality == Function::Local::PREV_NEXT;
  // next is always true!
  bool two_signs = locality == Function::Local::NEXT_AND_REVERSE || locality == Function::Local::PREV_NEXT_AND_REVERSE;

  int sign_1 = two_signs ? 1 : Function::Local::Identifier::NO_SIGN;
  int sign_2 = -1;

  // we don't set Function::Local::element_dimension_, it would be 2 (dim==1 * two signs)
  // we wont't use the full space as the individual shape_ are not connected
  // FIXME consider symmetry via opt_shape_param_ ?!
  vem.Reserve(num_node_shape_params_ * (two_signs ? 2 : 1)); // separately for node and profile but both of same size

  // traverse shape_ to have proper start and end. check for node or profile
  for(unsigned int s = GetFirstShapeIdx(f), n = GetEndShapeIdx(f); s < n; s++)
  {
    const ShapeParam& shape = shape_[s];

    if(!shape.fixed && !shape.IsInduced() && !shape.slave)
    {
      assert(shape.start_opt >= 0); // - 1 if not applicable
      assert(shape.end_opt > 0);

      bool periodic = shape.ShallMapHalfShape() ? false : f->GetLocal()->periodic; // symmetric stuff has special periodicity handling below

      LOG_DBG(SMD) << "SVSEM f=" << f->ToString() << " s=" << s << " ts=" << two_signs << " prev=" << prev << " per=" << periodic << " sp=" << shape.start_param << " ep=" << shape.end_param << " so=" << shape.start_opt << " eo=" << shape.end_opt;

      // skip the last element as we want only 'full' elements with next when we are not periodic
      for(int e = shape.start_opt + (prev ? (periodic ? 0 : 1) : 0), n = shape.end_opt - (periodic ? 0 : +1); e < n; e++)
      {
        BaseDesignElement* bde = opt_shape_param_[e];
        assert(f->GetDesignType() == bde->GetType());
        assert((!periodic && e < shape.end_opt-1) || (periodic && e < shape.end_opt));

        // note that opt_shape_param can be a fragmented variant of shape_param which is an issue with the index which needs to consecutive in opt_shape_param
        BaseDesignElement* prev_de = prev ? opt_shape_param_[e == shape.start_opt ? shape.end_opt-1 : e-1] : NULL; // if not prev take last
        BaseDesignElement* next_de =        opt_shape_param_[e == shape.end_opt-1 ? shape.start_opt : e+1]; // we next cannot be next we take first (only if periodic)

        LOG_DBG(SMD) << "SVSEM s=" << s << " n=" << n << " po=" << (prev_de != NULL ? (int) prev_de->GetOptIndex() : -1) << " eo=" << e << " no=" << (next_de != NULL ? (int) next_de->GetOptIndex() : -1)
                                                          << " p="  << (prev_de != NULL ? (int) prev_de->GetIndex() : -1)    << " e=" << bde->GetIndex() << " n=" << (next_de != NULL ? (int) next_de->GetIndex() : -1);


        vem.Push_back(Function::Local::Identifier(bde, prev_de, next_de, sign_1));
        if(two_signs)
          vem.Push_back(Function::Local::Identifier(bde, prev_de, next_de, sign_2));
      }

      // special case for curvature on a symmetric shape where we opt only half (e.g. 6) then a curvature 5 - 6 - 5 is not possible
      // as snopt complains. Therefore we do a slope constraint 5 - 6 with the bound for the curvature
      if(prev && shape.ShallMapHalfShape() && (shape.end_opt - shape.start_opt) > 2)
      {
        // curvature over the inner element
        // note that the slope should actually have half the curvature bound!
        BaseDesignElement* bde  = opt_shape_param_[shape.end_opt-2];
        BaseDesignElement* next = opt_shape_param_[shape.end_opt-1];

        LOG_DBG(SMD) << "SVSEM s=" << s << " curvature slope: eo=" << bde->GetOptIndex() << " no=" << next->GetOptIndex() << " e=" << bde->GetIndex() << " n=" << next->GetIndex();
        vem.Push_back(Function::Local::Identifier(bde, NULL, next, sign_1));
        if(two_signs)
          vem.Push_back(Function::Local::Identifier(bde, NULL, next, sign_2));

        // curvature over the outer element last - 0 - 1 which would be 1 - 0 - 1 which is slope 0 - 1 (should be actually half bound)
        bde  = opt_shape_param_[shape.start_opt];
        next = opt_shape_param_[shape.start_opt+1];

        LOG_DBG(SMD) << "SVSEM s=" << s << " curvature slope: eo=" << bde->GetOptIndex() << " no=" << next->GetOptIndex() << " e=" << bde->GetIndex() << " n=" << next->GetIndex();
        vem.Push_back(Function::Local::Identifier(bde, NULL, next, sign_1));
        if(two_signs)
          vem.Push_back(Function::Local::Identifier(bde, NULL, next, sign_2));
      }
    }
  }

  LOG_DBG(SMD) << "SVSEM final f=" << f->ToString() << " loc=" << locality << " ts=" << two_signs << " prev=" << prev << " -> vem=" << vem.GetSize();
}

void ShapeMapDesign::SetupVirtualMultiShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& vem, Function::Local::Locality locality)
{
  assert(f != NULL);
  assert(f->IsLocal(f->GetType()));
  // we shall be called by Local::PostInit() therefore local shall exist
  assert(f->GetLocal() != NULL);
  assert(!f->GetLocal()->periodic);
  // we assume fixed only for profile
  if(IsProfileFixed())
    throw Exception("Configuration error: cannot have local constraint of 'shape_map' when 'profile' is fixed.");

  assert(locality == Function::Local::MULT_DESIGNS_PREV_NEXT_AND_REVERSE || locality == Function::Local::MULT_DESIGNS_PREV_NEXT || locality == Function::Local::MULT_DESIGNS_NEXT_AND_REVERSE || locality == Function::Local::MULT_DESIGNS_NEXT);

  // a lot copy&paste from SetupVirtualShapeElementMap()
  bool prev = locality == Function::Local::MULT_DESIGNS_PREV_NEXT_AND_REVERSE || locality == Function::Local::MULT_DESIGNS_PREV_NEXT;
  // next is always true!
  bool two_signs = locality == Function::Local::MULT_DESIGNS_NEXT_AND_REVERSE || locality == Function::Local::MULT_DESIGNS_NEXT_AND_REVERSE;

  int sign_1 = two_signs ? 1 : Function::Local::Identifier::NO_SIGN;
  int sign_2 = -1;

  // in principle other functions are also possible but these two are for either vertical or horizontal structures.
  assert(f->GetType() == Function::OVERHANG_HOR || f->GetType() == Function::OVERHANG_VERT);
  assert(f->GetDesignType() == DesignElement::SHAPE_MAP);
  ShapeParamElement::Dof dof = f->GetType() == Function::OVERHANG_HOR ? ShapeParamElement::Y : ShapeParamElement::X;
  StdVector<ShapeParam*> shapes = FindShape(NODE, dof);
  assert(num_node_shapes_ >= (int) shapes.GetSize() / 2);

  if(shapes.IsEmpty())
    throw Exception("There are no shape variables for function '" + f->ToString() + "'");

  vem.Reserve(shapes.GetSize() * 2 * (two_signs ? 2 : 1)); // common for node and profile

  StdVector<BaseDesignElement*> buddies; // to be reused temporary vector

  // traverse nodes only and the the corresponding profiles implicitly
  // do
  for(unsigned int si = 0; si < shapes.GetSize(); si++)
  {
    unsigned int s = shapes[si]->idx;
    const ShapeParam& node = shape_[s];
    const ShapeParam& prof  = *node.partner;
    // don't deal with the complicated stuff!
    assert(!node.fixed);
    assert(!prof.fixed);
    assert(!node.IsInduced());
    assert(!node.ShallMapHalfShape());
    assert((node.dof == 0 && f->GetType() == Function::OVERHANG_VERT) || (node.dof == 1 && f->GetType() == Function::OVERHANG_HOR));
    // LOG_DBG(SMD) << "SVSEM f=" << f->ToString() << " s=" << s << " ts=" << two_signs << " prev=" << prev << " per=" << periodic << " sp=" << shape.start_param << " ep=" << shape.end_param << " so=" << shape.start_opt << " eo=" << shape.end_opt;

    // skip the last element as we want only 'full' elements with next when we are not periodic
    for(int e = node.start_opt + (prev ? 1 : 0), n = node.end_opt - 1; e < n; e++) // we always assume 'next'
    {
      BaseDesignElement* bde = opt_shape_param_[e];
      assert(bde->GetType() == DesignElement::NODE);

      // note that opt_shape_param can be a fragmented variant of shape_param which is an issue with the index which needs to consecutive in opt_shape_param
      BaseDesignElement* prev_de = prev ? opt_shape_param_[e == node.start_opt ? node.end_opt-1 : e-1] : NULL; // if not prev take last
      BaseDesignElement* next_de =        opt_shape_param_[e == node.end_opt-1 ? node.start_opt : e+1]; // we next cannot be next we take first (only if periodic)

      // the profile for the nodes
      BaseDesignElement* bde_pr = opt_shape_param_[prof.start_opt + (bde->GetOptIndex() - node.start_opt)];
      BaseDesignElement* prev_pr = prev_de != NULL ? opt_shape_param_[prof.start_opt + (prev_de->GetOptIndex() - node.start_opt)] : NULL;
      BaseDesignElement* next_pr = next_de != NULL ? opt_shape_param_[prof.start_opt + (next_de->GetOptIndex() - node.start_opt)] : NULL;

      buddies.Clear(true); // this(node)->elem is implicit, then this(profile), then prev(node) and prev(profile) if exist, then next(node) and next(profile)
      buddies.Push_back(bde_pr);
      if(prev_de && prev_pr) {
        buddies.Push_back(prev_de);
        buddies.Push_back(prev_pr);
      }
      buddies.Push_back(next_de);
      buddies.Push_back(next_pr);

      vem.Push_back(Function::Local::Identifier(bde, buddies, sign_1));
      if(two_signs)
        vem.Push_back(Function::Local::Identifier(bde, buddies, sign_2));


      LOG_DBG(SMD) << "SVMSEM s=" << s << " n=" << n << " prev_opt_node=" << (prev_de != NULL ? (int) prev_de->GetOptIndex() : -1) << " node =" << e
          << " next_opt_node=" << (next_de != NULL ? (int) next_de->GetOptIndex() : -1)
          << " prev_opt_prof=" << (prev_pr != NULL ? (int) prev_pr->GetOptIndex() : -1) << " prof =" << bde_pr->GetOptIndex()
          << " next_opt_prof=" << (next_pr != NULL ? (int) next_pr->GetOptIndex() : -1);

    }
  }

  LOG_DBG(SMD) << "SVMSEM final f=" << f->ToString() << " loc=" << locality << " ts=" << two_signs << " prev=" << prev << " -> vem=" << vem.GetSize();
}



void ShapeMapDesign::SetupCyclicVirtualShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& vem, Function::Local::Locality locality)
{
  // we assume fixed only for profile
  if(f->GetDesignType() == BaseDesignElement::PROFILE && IsProfileFixed())
    throw Exception("cannot have local constraint of shape map design 'profile' when this design is fixed.");

  assert(locality == Function::Local::CYCLIC);

  vem.Reserve(shape_.GetSize());
  assert(vem.GetCapacity() > 0);

  for(unsigned int s = GetFirstShapeIdx(f), n = GetEndShapeIdx(f); s < n; s++)
  {
    const ShapeParam& shape = shape_[s];

    assert(!(f->GetDesignType() == BaseDesignElement::NODE && shape.type == PROFILE));
    assert(!(f->GetDesignType() == BaseDesignElement::PROFILE && shape.type == NODE));
    if(!shape.fixed && !shape.IsInduced())
    {
      assert(f->GetDesignType() == Convert(shape.type)); // NODE or PROFILE

      BaseDesignElement* left  = opt_shape_param_[shape.start_opt];
      BaseDesignElement* right = opt_shape_param_[shape.end_opt-1]; // now take the last which is end-1

      vem.Push_back(Function::Local::Identifier(left, NULL, right)); // makes a neighborhood of size 1
    }
  }

  LOG_DBG(SMD) << "SCVSEM f=" << f->ToString() << " loc=" << locality << " -> vem=" << vem.GetSize();
}


int ShapeMapDesign::ReadDesignFromExtern(const double* space_in, bool setAndWriteCurrent)
{
  // TODO: This is almost in FeaturedDesign()


  assert(!std::isnan(scaling_));
  int old_design = design_id;

  // write aux design variables (slack and alpha if any) last
  assert(export_fe_design_ == false); // we do shape map
  assert(DesignSpace::GetNumberOfVariables() > 0); // we need this variables but they are hidden!

  bool new_design = false;

  for(unsigned int i = 0, n = opt_shape_param_.GetSize(); i < n; i++)
  {
    double v = space_in[i] * scaling_;
    assert(!std::isnan(v));
    if(!new_design && v != opt_shape_param_[i]->GetPlainDesignValue())
      new_design = true;

    opt_shape_param_[i]->SetDesign(v);

    if(dynamic_cast<ShapeMapVariable*>(opt_shape_param_[i])->sym->HasSymmetry())
      dynamic_cast<ShapeMapVariable*>(opt_shape_param_[i])->sym->ApplyDesign();

    LOG_DBG3(SMD) << "RDFE: i=" << i << "=" << opt_shape_param_[i]->ToString() << " -> " << v;
  }

  // append aux design, might also change design_id
  AuxDesign::ReadDesignFromExtern(space_in, setAndWriteCurrent); // note the asserts above!

  if(new_design && design_id <= old_design)
    design_id++; // if new design and not already changed by AuxDesign

  // the design are shape parameters, map them to rho
  if(mapped_design_ != design_id)
    MapFeatureToDensity();

  LOG_DBG(SMD) << "RDFE: di -> " << design_id;
  return design_id;
}


int ShapeMapDesign::WriteDesignToExtern(double* space_out, bool scaling) const
{
  FeaturedDesign::WriteDesignToExtern(space_out, scaling);

  if(export_leveset_)
    ExportLevelSet();

  return design_id;
}

void ShapeMapDesign::WriteGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Function* f, bool scaling)
{
  // TODO: Almost same as FeaturedDesign::WriteGradientToExtern()


  LOG_DBG2(SMD) << "WGTE: f=" << f->ToString() << " ad=" << aux_design_.GetSize() << " osp=" << opt_shape_param_.GetSize() << " out=" << out.GetSize() << " owst=" << out.window.GetStart() << " owsz=" << out.window.GetSize();
  assert(out.window.GetStart() + out.window.GetSize() <= out.GetSize());

  // we cannot cache easily for mapped_constr_gradient_ as we would need it for each function.
  // MapShapeGradient would be good to perform it for all functions concurrently, however this is not possible as it is not the case that first all
  // simp function gradients are called and then all exported. This would need rewriting some stuff in cfs!
  if(f->IsObjective() || !Function::IsLocal(f->GetType())) // don't map local functions
    MapFeatureGradient(f); // see comment above for what is necessary to cache the stuff

  assert(f != NULL);
  assert(opt_->objectives.data.GetSize() == 1); // implement multi objective and be careful!
  assert(!(vs == DesignElement::COST_GRADIENT && !f->IsObjective()));
  assert(vs == DesignElement::COST_GRADIENT || vs == DesignElement::CONSTRAINT_GRADIENT);
  assert(!opt_->GetMultipleExcitation()->DoMetaExcitation(f->ctxt)); // robustness and transformation don't make sense for shape map. In DesignSpace we do f->GetExcitation()->Apply()
  assert(design.GetSize() == 1); // only pseudo density, nothing else implemented yet
  assert(regions.GetSize() == 1); // needs to be as design shall be 1

  assert(out.window.Initialized());
  unsigned int base = out.window.GetStart();
  // out contains the Jacobian for a possibly many functions. Snopt combines cost function (first) and then the constraints derivatives
  // Here we assume that out has a window set where to write to for the given function.
  if(f->HasDenseJacobian())
  {
    LOG_DBG(SMD) << "WGTE: fvi=" << GetFirstOptVarIdx(f) << " evi=" << GetEndOptVarIdx(f) << " ad=" << aux_design_.GetSize() << " w=" << out.window.GetSize();
    assert(GetEndOptVarIdx(f) - GetFirstOptVarIdx(f) + aux_design_.GetSize() == out.window.GetSize());
    for(unsigned int s = GetFirstOptVarIdx(f), n = GetEndOptVarIdx(f); s < n ; s++) // for node (from 0) and profile (later) or default for both
    {
      assert(out.InWindow(base + s));

      double opt = opt_shape_param_[s]->GetPlainGradient(f);
      double sym = dynamic_cast<ShapeMapVariable*>(opt_shape_param_[s])->sym->GetPlainSymGradient(f);

      LOG_DBG3(SMD) << "WGTE oe=" << opt_shape_param_[s]->ToString();
      assert(!std::isnan(opt));
      assert(!std::isnan(sym));

      out[base + s] = (opt + sym) * scaling;

      LOG_DBG3(SMD) << "WGTE f=" << f->ToString() << " ws=" << out.window.GetStart() << " s=" << s << " opt=" << opt << " sym=" << sym << " -> " << out[base + s];
    }
    // add slack stuff. No need to cheat window size
    AuxDesign::WriteGradientToExtern(out, vs, access, f, scaling);
  }
  else
  {
    assert(f->GetDesignType() == BaseDesignElement::NODE || f->GetDesignType() == BaseDesignElement::PROFILE || f->GetDesignType() == BaseDesignElement::SHAPE_MAP);
    // uses opt_index_!
    StdVector<unsigned int>& sparsity = f->GetSparsityPattern();
    LOG_DBG2(SMD) << "WGTE f=" << f->ToString() << " sparsity=" << sparsity.ToString();
    assert(out.window.GetSize() == sparsity.GetSize());
    for(unsigned int i = 0; i < sparsity.GetSize(); i++)
    {
      unsigned int s = sparsity[i];
      assert(s < opt_shape_param_.GetSize());
      LOG_DBG3(SMD) << "WGTE i=" << i << " s=" << s << " base=" << base << " b+s=" << (base+s);

      assert(out.InWindow(base + i));
      double scale = scaling ? scaling_ : 1.0;
      assert(vs == BaseDesignElement::CONSTRAINT_GRADIENT);

      double opt = opt_shape_param_[s]->GetPlainGradient(f);
      double sym = dynamic_cast<ShapeMapVariable*>(opt_shape_param_[s])->sym->GetPlainSymGradient(f);

      out[base + i] = (opt + sym) * scale;
    }
  }
}

void ShapeMapDesign::Reset(DesignElement::ValueSpecifier vs, DesignElement::Type design)
{
  assert(design == BaseDesignElement::DEFAULT || design == BaseDesignElement::NODE); // extend to check for profile

  FeaturedDesign::Reset(vs, design);

  for(ShapeParamElement* spe : opt_shape_param_) {
    dynamic_cast<ShapeMapVariable*>(spe)->sym->Reset(vs);
    assert(!spe->costGradient.IsEmpty());
  }
}


int ShapeMapDesign::FindDesign(DesignElement::Type dt, bool throw_exception) const
{
  // check for DENSITY, ...
  int idx = DesignSpace::FindDesign(dt, false);
  if(idx >= 0)
    return idx;

  assert(dt == DesignElement::NODE || dt == DesignElement::PROFILE || dt == DesignElement::SHAPE_MAP);

  if(dt == DesignElement::SHAPE_MAP)
    dt = DesignElement::NODE; // return the node index

  for(unsigned int i = 0; i < shape_.GetSize(); i++)
    if(Convert(shape_[i].type) == dt)
      return i;

  if(throw_exception)
    EXCEPTION("Design " << DesignElement::type.ToString(dt) << " no FEM based and no shape mapping design.");
  return -1;
}

void ShapeMapDesign::ReadDensityXml(PtrParamNode set, double& lower_violation, double& upper_violation)
{
  ParamNodeList list = set->GetList("shapeParamElement");
  if(list.IsEmpty())
    throw Exception("no 'shapeParamElement' found in provided density.xml");

  assert(!shape_param_.IsEmpty());
  if(list.GetSize() != shape_param_.GetSize() && list.GetSize() != shape_param_.GetSize()/2)
    EXCEPTION("incompatible shape variables in density.xml (" << list.GetSize() << ") with " << shape_param_.GetSize() << " variables expected");

  // needs to be checked while reading
  bool read_profile = list.GetSize() != shape_param_.GetSize() / 2; // FIXME! No issue known but dangerous
  info_->Get("ersatzMaterialFile")->Get("load")->SetValue(read_profile ? "node_and_profile" : "only_node");
  for(unsigned int i = 0; i < list.GetSize(); i++)
  {
    const PtrParamNode pn = list[i];
    ShapeMapVariable* spe = dynamic_cast<ShapeMapVariable*>(shape_param_[i]);
    ShapeParam* shape = GetShape(spe);
    assert(spe->GetIndex() == i);
    assert(pn->Get("nr")->As<unsigned int>() == i);

    double val = pn->Get("design")->As<double>();
    BaseDesignElement::Type dt = BaseDesignElement::type.Parse(pn->Get("type")->As<string>());
    LOG_DBG2(SMD) << "RDX i=" << i << " dt=" << dt << " dof=" << (pn->Has("dof") ? pn->Get("dof")->As<string>() : "-") << " val=" << val << " spe=" << spe->ToString();

    // a profile has no dof
    if(dt != spe->GetType())
      EXCEPTION("shapeParamElement nr=" << i << " has type " << pn->Get("type")->As<string>()
          << " but we expect type " << BaseDesignElement::type.ToString(spe->GetType()));

    if(!read_profile && dt == BaseDesignElement::PROFILE)
      EXCEPTION("no shapeParamElement of type profile expected in density.xml");

    if(dt == BaseDesignElement::NODE && pn->Get("dof")->As<string>() != ShapeParamElement::dof.ToString(spe->dof_))
      EXCEPTION("shapeParamElement nr " << i << " has not the expected dof " << ShapeParamElement::dof.ToString(spe->dof_));

    if(!shape->fixed && !shape->IsInduced()) // with fixed we don't read bounds
    {
      lower_violation = std::max(lower_violation, spe->GetLowerBound() - val);
      upper_violation = std::max(upper_violation, val - spe->GetUpperBound());
      LOG_DBG2(SMD) << "RDX: e=" << spe->GetIndex() << " v=" << val << " lb=" << spe->GetLowerBound() << " ub=" << spe->GetUpperBound() << " lv=" << lower_violation << " uv=" << upper_violation;

      if(enforce_bounds_) {
        spe->SetDesign(std::max(spe->GetLowerBound(), std::min(spe->GetUpperBound(), val)));
        val = spe->GetPlainDesignValue();
      }
    }
    spe->SetDesign(val); // possibly corrected val

    // Get value of the relative bound for current design variable. If value not set, db = -1.
    double rb = spe->GetType() == DesignElement::NODE ? relative_node_bound_ : relative_profile_bound_;
    assert(!(shape->fixed && rb >= 0.0));

    // consider ShapeParam::clamped which overwrites relative_*_bound for the first and last element
    if(shape->clamp >= 0.0 && ((int) i == shape->start_param || (int) i == shape->end_param-1)) {
      LOG_DBG(SMD) << "RDX: set clamped idx=" << spe->GetIndex() << " shape=" << shape->ToString() << " rb=" << rb << " -> " << shape->clamp;
      rb = shape->clamp;
    }

    if(rb >= 0.0)
    {
      LOG_DBG2(SMD) << "RDX: before i=" << i << " v=" << val << " rb=" << rb << " lb = " << spe->GetLowerBound() << " ub=" << spe->GetUpperBound();
      // if a relative_bound is set in the xml file, upper and lower bound are overwritten
      // assume that the initial value is out of original bound (e.g. too small), this needs to be catched
      spe->SetUpperBound(std::min(spe->GetUpperBound(), std::max(val + rb, spe->GetLowerBound())));
      spe->SetLowerBound(std::max(spe->GetLowerBound(), std::min(val - rb, spe->GetUpperBound())));
    }
    LOG_DBG2(SMD) << "RDX: e=" << i << spe->ToString() << " v=" << val << " rb=" << rb << " lb=" << spe->GetLowerBound() << " ub=" << spe->GetUpperBound();
    assert(spe->GetLowerBound() <= spe->GetUpperBound());
  }

  // now apply symmetries. The density.xml doesn't make a difference between opt and sym variables
  for(ShapeParamElement* var : opt_shape_param_)
  {
    ShapeMapVariable* opt = dynamic_cast<ShapeMapVariable*>(var);
    if(opt->sym->HasSymmetry())
      opt->sym->ApplyDesign();
  }

  MapFeatureToDensity();
}

inline void ShapeMapDesign::MapIdxToCoords(const Vector<unsigned int>& idx, Vector<double>& out) const
{
  out[0] = coord_min_[0] + idx[0] * coord_step_[0];
  out[1] = coord_min_[1] + idx[1] * coord_step_[1];
  if(dim_ == 3)
    out[2] = coord_min_[2] + idx[2] * coord_step_[2];

  assert(out[0] <= coord_max_[0] + 1e-10);
  assert(out[1] <= coord_max_[1] + 1e-10);
  assert(dim_ == 2 || out[2] <= coord_max_[2] + 1e-10);
}

/** This is the inverse to DensityIdx */
inline void ShapeMapDesign::DensityIdx(unsigned int i, Vector<unsigned int>& idx) const
{
  assert(idx.GetSize() >= dim_);
  if(dim_ == 2) // i = y * nx_ + x;
  {
    idx[1] = i / nx_;
    idx[0] = i - nx_ * idx[1];

    assert(DensityIdx(idx[0], idx[1]) == i);

    if(idx.GetSize() == 3)
      idx[2] = 0;
  }
  else // i = z * nx_*ny_ + y * nx_ + x
  {
    idx[2] = i / (nx_*ny_);
    unsigned int r1 = i - (nx_*ny_) * idx[2];
    idx[1] = r1 / nx_;
    idx[0] = r1 - (idx[1] * nx_);

    assert(DensityIdx(idx[0], idx[1], idx[2]) == i);
  }
}


BaseDesignElement::Type ShapeMapDesign::Convert(Type type)
{
  switch(type)
  {
  case NODE:
    return BaseDesignElement::NODE;
  case PROFILE:
    return BaseDesignElement::PROFILE;
  case CENTER:
    assert(false); // shall not happen as CENTER is only the virtual combination of two NDOE shapes
    break;
  }
  assert(false);
  return BaseDesignElement::NO_TYPE;
}

inline bool ShapeMapDesign::IsProfileFixed() const
{
  // assume all nodes are concurrently fixed?!
  assert(shape_.Last().type == PROFILE);
  return shape_.Last().fixed;
}

bool ShapeMapDesign::IsAllNodeFixed() const
{
  for(int i = 0; i < num_node_shapes_; i++)
    if(!shape_[i].fixed)
      return false;

  return true;
}


unsigned int ShapeMapDesign::GetFirstOptVarIdx(const Function* f) const
{
  const ShapeParam* shape = FindShape(f, true);
  assert(shape != NULL);
  return shape->start_opt;
}

/** small helper which gives the  index *after* the element based on type (node or profile) so*/
unsigned int ShapeMapDesign::GetEndOptVarIdx(const Function* f) const
{
  assert(f->GetDesignType() == BaseDesignElement::DEFAULT || f->GetDesignType() == BaseDesignElement::DENSITY || f->GetDesignType() == BaseDesignElement::PROFILE || f->GetDesignType() == BaseDesignElement::NODE);
  // either opt_shape_param_.GetSize() if not only NODE or end of NODE
  if(f->GetDesignType() != BaseDesignElement::NODE || IsProfileFixed())
    return opt_shape_param_.GetSize();

  // end of node var is the first profile var which is shape_[num_node_shapes_].start_opt
  assert(f->GetDesignType() != BaseDesignElement::PROFILE);
  assert(num_node_shapes_ < (int) shape_.GetSize());
  assert(shape_[num_node_shapes_].type == PROFILE);
  assert(shape_[num_node_shapes_].start_opt >= 0);

  return shape_[num_node_shapes_].start_opt;
}

unsigned int ShapeMapDesign::GetFirstShapeIdx(const Function* f) const
{
  const ShapeParam* shape = FindShape(f, true);
  return shape->idx;
}

/** small helper which gives the  index *after* the element based on type (node or profile) so*/
unsigned int ShapeMapDesign::GetEndShapeIdx(const Function* f) const
{
  // NODE, PROFILE or SHAPE_MAP!
  BaseDesignElement::Type dt = f->GetDesignType();
  assert(dt == BaseDesignElement::DEFAULT || dt == BaseDesignElement::DENSITY || BaseDesignElement::IsShapeMapType(dt));
  assert(2 * num_node_shapes_ >= (int) shape_.GetSize()); // end of profile
  if(dt == BaseDesignElement::NODE || IsProfileFixed())
    return num_node_shapes_;
  assert(num_node_shapes_ < (int) shape_.GetSize());
  assert(shape_[num_node_shapes_].type == PROFILE);
  assert(!shape_[num_node_shapes_].fixed);
  return shape_.GetSize();
}


void ShapeMapDesign::ToInfo(ErsatzMaterial* em)
{
  AuxDesign::ToInfo(em);

  PtrParamNode sm = info_->Get("shapeMap");
  sm->Get("overlap")->SetValue(combine.ToString(combine_));
  PtrParamNode msh = sm->Get("mesh");
  msh->Get("n")->SetValue(n_.ToString());
  msh->Get("min")->SetValue(coord_min_.ToString());
  msh->Get("max")->SetValue(coord_max_.ToString());
  msh->Get("step")->SetValue(coord_step_.ToString());
  numInt_.ToInfo(sm->Get("numInt"));
  PtrParamNode base = info_->Get("designVariables");
  for(unsigned int i = 0; i < shape_.GetSize(); i++)
    shape_[i].ToInfo(base->Get("shapeParam", ParamNode::APPEND));
}


// also used by LatticeBoltzmannPDE !
Vector<unsigned int> ShapeMapDesign::SetupLexicographicMesh(Grid* grid, const RegionIdType design_reg, StdVector<int>& elem_to_idx, StdVector<int>& idx_to_elem)
{
  // for 2D n_z_=1

  // TODO use FeaturedDesign:SetupMeshStructure()

  StdVector<unsigned int> t = grid->GetRegularDiscretization(design_reg);
  LOG_DBG(SMD) << "SLM r=" << design_reg << " -> " << t.ToString();
  Vector<unsigned int> n(3);
  Copy(t, n);

  unsigned int n_elems = n[0] * n[1] * n[2];

  // TODO here we assume that the whole mesh is for LBM and the mesh is of lexicographic ordering.
  // To be good this needs to handled by element neighbors!
  // this shall better be assuered technically!

  idx_to_elem.Resize(n_elems);
  for(unsigned int i = 0; i < n_elems; i++)
    idx_to_elem[i] = i+1; // one based

  elem_to_idx.Resize(n_elems + 1); // one-based elem_nr
  for(unsigned int i = 0, n = elem_to_idx.GetSize(); i < n; i++)
    elem_to_idx[i] = i-1; // -1 for not appropriate idx

  return n;
}


int ShapeMapDesign::Item::GetOrder(Vector<int>& order, const ShapeMapDesign::NumInt& ni) const
{
  assert(order.GetSize() == min_corner_value.GetSize());
  assert(order.GetSize() == max_corner_value.GetSize());
  int max = 0;
  if(ni.strategy_ == FeaturedDesign::TAILORED)
  {
    for(unsigned int s = 0; s < order.GetSize(); s++)
    {
      double min_val = min_corner_value[s];
      double max_val = max_corner_value[s];

      assert(min_val >= 0);
      assert(max_val >= min_val);
      assert(max_val <= 1);

      if(max_val < ni.sensitivity_)
        order[s] = 0;
      else if(min_val > 1-ni.sensitivity_)
        order[s] = 1;
      else
        order[s] = ni.GetTailoredOrder(max_val - min_val);

      max = std::max(max, order[s]);
    }
  }
  else
  {
    max = FeaturedDesign::Item::GetOrder(order, ni);
  }

  return max;
}


inline double ShapeMapDesign::Item::SetIPGetWeight(const ShapeMapDesign* smd, StdVector<double>& ip, int ip_x, int ip_y, int ip_z, int order)
{
  assert(ip.GetSize() == dim_);
  assert(ip_x >= 0 && ip_x < order);
  assert(ip_y >= 0 && ip_y < order);
  assert((dim_ == 2 && ip_z == 0) || (dim_ == 3 && ip_z >= 0 && ip_z < order));

  ip[0] = (double) ip_x / (double) (order-1);
  ip[1] = (double) ip_y / (double) (order-1);
  if(dim_ == 3)
    ip[2] = (double) ip_z / (double) (order-1);

  assert(ip[0] >= 0 && ip[0] <= 1);
  assert(ip[1] >= 0 && ip[1] <= 1);
  assert(dim_ == 2 || (ip[2] >= 0 && ip[2] <= 1));

  if(smd->GetBoundary() == TANH)
  {
    assert(order >= 2  && order <= (int) ShapeMapDesign::newtonCotes.GetSize());
    const Vector<double>& w = ShapeMapDesign::newtonCotes[order-1];
    assert((int) w.GetSize() == order);

    return w[ip_x] * w[ip_y] * (dim_ == 3 ? w[ip_z] : 1.0); // might be negative
  }
  else
  {
    assert(smd->GetBoundary() == LINEAR);
    assert(order > 1);

    double weight = 1.0;
    weight *= 1./(order-1) * (ip_x == 0 || ip_x == order-1 ? 0.5 : 1.0);
    weight *= 1./(order-1) * (ip_y == 0 || ip_y == order-1 ? 0.5 : 1.0);
    if(dim_ == 3)
      weight *= 1./(order-1) * (ip_z == 0 || ip_z == order-1 ? 0.5 : 1.0);

    return weight;
  }
}

void ShapeMapDesign::NumInt::Init(ShapeMapDesign* smd, PtrParamNode pn, PtrParamNode info)
{
  this->sf_  = smd->GetBoundary();
  this->beta_ = smd->GetBeta();
  Vector<unsigned int> n = smd->GetDiscretization();

  FeaturedDesign::NumInt::Init(smd, pn, info);

  if(max_order_ < 2)
    info->SetWarning("minimal value for 'max_order' is '2'");
  if(smd->GetBoundary() == TANH && max_order_ > (int) ShapeMapDesign::newtonCotes.Last().GetSize())
    info->SetWarning("maximal value for 'max_order' is " + to_string(ShapeMapDesign::newtonCotes.Last().GetSize()));

  assert(n[n.GetSize()-1] != 0);

  unsigned int res = std::min(n[0], domain->GetDim() == 2 ? n[1] : std::min(n[1], n[2]));
  h = 1.0/res;
  assert(res >= 2);

  if(strategy_ == TAILORED)
  {
    if(sf_ == TANH)
      SetTailoredTanh(info);
    else
      SetLinearIntOrder(info);
  }
}

/** searches tailored_bounds and returns the appropriate tailored_order content */
inline int ShapeMapDesign::NumInt::GetTailoredOrder(double max_min) const
{
  assert(!(sf_ == LINEAR && linear_int_order_ < 2));
  if(sf_ == LINEAR)
    return linear_int_order_;

  for(int oi = tailored_bounds_.GetSize()-1; oi > 0; oi--) // 10,9, ...,1 which corresponds to 1,0.1, 0.01, 0.001, ...
    if(max_min > tailored_bounds_[oi])
      return tailored_order_[oi];

  return tailored_order_[0];
}


void ShapeMapDesign::NumInt::SetTailoredTanh(PtrParamNode info)
{
  // see svn+ssh://eamc080/home/svn_repo/repository/publications/geometry_projection/plots/tanh.py

  tailored_bounds_ = LogspaceBase(-10, 0, 11);
  tailored_order_.Resize(11, 2);

  // the bounds for pos from 0.5 to the next element by 10 steps
  for(double pos = 0.5; pos <= 0.5 + h; pos += h/10)
  {
    // sample the elements
    for(double x1 = 0; x1 < 1; x1+=h)
    {
      double x2 = x1+h;
      double max_min = Func(x1, pos) - Func(x2,pos);
      assert(max_min >= 0);
      assert(max_min <= 1);

      int order = FindOrder(x1, x2, pos, sensitivity_);
      // search where we are in the logscale range
      for(int oi = tailored_bounds_.GetSize()-1; oi >= 0; oi--) // 10,9, ...,1,0 which corresponds to 1,0.1, 0.01, 0.001, ...
      {
        if(max_min > tailored_bounds_[oi])
        {
          // e.g. oi = 9, max_min = .34 which is > 0.1
          tailored_order_[oi] = std::max(tailored_order_[oi], order);
          break;
        }
      }
      LOG_DBG2(SMD) << "NI:ST p=" << pos << " x1=" << x1 << " x2=" << x2
          << " mm=" << max_min << " o=" << order << " -> " << tailored_order_.ToString();
    } // x1 loop
  } // pos loop

  // check for too high ordering and assure monotone ordering!
  int max_possible_order = ShapeMapDesign::newtonCotes.Last().GetSize();
  LOG_DBG(SMD) << "NI:ST mpo=" << max_possible_order << " mo" << max_order_ << " max=" << tailored_order_.Max()
           << " tb=" << tailored_bounds_.ToString() << " org to=" << tailored_order_.ToString();

  if(tailored_order_.Max() > max_possible_order)
    info->SetWarning("resolution and beta combination allows no sufficient numerical integration for sensitvity " + lexical_cast<string>(sensitivity_), true);

  if(tailored_order_.Max() > max_order_ && max_order_ < max_possible_order)
    info->SetWarning("configuration 'integration_order=" + to_string(max_order_) + "' where " + to_string(std::min(max_possible_order, tailored_order_.Max())) + " is required", true);
  assert(tailored_order_[0] >= 2);
  for(unsigned int i = 1; i < tailored_order_.GetSize(); i++) {
    tailored_order_[i-1] = std::min(tailored_order_[i-1], std::min(max_order_, max_possible_order));
    tailored_order_[i] = std::max(tailored_order_[i], tailored_order_[i-1]);
  }

  LOG_DBG(SMD) << "NI:ST clean to=" << tailored_order_.ToString();
}


void ShapeMapDesign::NumInt::SetLinearIntOrder(PtrParamNode info)
{
  // consider a 1D piecewise linear shape function. It has slope h = 1/n.
  // We do numerical integration to handle the arbitrary complex aggregation of shapes.
  // As the gray region is only 1/n there is not much to integrate with TAILORED.
  //
  // We check the error for the gradient which is constant 0 or h.
  // The integration points are at an h/o spacing with o the order.
  // When the jump is within the element the grad value is h on one side of the interval and 0 on the other side.
  // numerical integration of first order gives an integral of 1/2 * h * h/o.
  // The extreme error is when the jump is at the left or right side. The error in this case is .5*h^2/o.
  // Setting error to the sensitivity, the order results in o = .5 * h^2 / sensitivity

  assert(sensitivity_ > 1e-14);

  linear_int_order_ = .5 * (h*h) / sensitivity_;
  linear_int_order_ = std::max(linear_int_order_, 2);

  if(linear_int_order_ > max_order_) {
    std::stringstream ss;
    ss << "tailored linear integration order would require " << linear_int_order_ << " with sensitivity="
        << sensitivity_ << " for h=" << h << " -> cut to " << max_order_;
    info->SetWarning(ss.str());
    linear_int_order_ = max_order_;
  }
}

void ShapeMapDesign::NumInt::ToInfo(PtrParamNode info) const
{
  FeaturedDesign::NumInt::ToInfo(info);

  if(strategy_ == TAILORED)
  {
    if(sf_ == TANH)
    {
      // we print only downwards up to the first order 2
      std::stringstream ss;
      bool done = false;
      for(int oi = tailored_order_.GetSize()-1; oi >= 0 && !done; oi--) {
        ss << tailored_bounds_[oi] << "->" << tailored_order_[oi] << " ";
        if(tailored_order_[oi] <= 2)
          done = true;
      }
      info->Get("tailored_order")->SetValue(ss.str());
    }
    else
    {
      info->Get("h")->SetValue(h);
      info->Get("order")->SetValue(linear_int_order_);
    }
  }
}

double ShapeMapDesign::NumInt::Func(double x, double pos) const
{
  return 1/(exp(beta_*(x-pos)) + 1);
}

double ShapeMapDesign::NumInt::GradFunc(double x, double pos) const
{
  double e = exp(beta_*(x-pos));
  return -(beta_ * e)/((e+1)*(e+1));
}

double ShapeMapDesign::NumInt::IntGradError(double x1, double x2, double pos, int order) const
{
  assert(x2 > x1);
  // the integral of GradFunc is clearly Func itself :)
  double corr = (Func(x2, pos) - Func(x1, pos)) / (x2 - x1);

  LOG_DBG3(SMD) << "NI:IGE x1=" << x1 << " x2=" << x2 << " p=" << pos << " b=" << beta_ <<  " F2=" << Func(x2, pos) << " F1=" << Func(x1, pos);

  assert(order >= 2  && order <= (int) ShapeMapDesign::newtonCotes.GetSize());
  const Vector<double>& w = ShapeMapDesign::newtonCotes[order-1];
  assert((int) w.GetSize() == order);
  double sum = 0;
  for(int i = 0; i < order; i++)
    sum += w[i] * GradFunc(x1 + i * (x2-x1)/(order-1), pos);

  LOG_DBG3(SMD) << "NI:IGE x1=" << x1 << " x2=" << x2 << " p=" << pos
      << " o=" << order <<  " corr=" << corr << " sum=" << sum << "->" << std::abs(corr-sum);

  return std::abs(corr-sum);
}


int ShapeMapDesign::NumInt::FindOrder(double x1, double x2, double pos, double accuracy) const
{
  // our order is the number of int points which is one more than the official newton cotes order
  assert(ShapeMapDesign::newtonCotes[0].GetSize() == 0);
  assert(ShapeMapDesign::newtonCotes[1].GetSize() == 2);
  int limit = ShapeMapDesign::newtonCotes.Last().GetSize();
  for(int o = 2; o <= limit; o++) // note the the limit itself is a valid number!
    if(IntGradError(x1,x2,pos,o) < accuracy)
      return o;
  // not found. Error needs to be handled!
  LOG_DBG(SMD) << "NI:FO order exceeded, return " << limit+1;
  return limit+1;
}


void ShapeMapDesign::MapFeatureToDensity()
{
  assert(data.GetSize() == map_.GetSize());
  mapping_timer_->Start();

  LOG_DBG(SMD) << "MSTD: di=" << design_id;


  int res_idx_r = GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::SHAPE_MAP_ORDER);
  int res_idx_c = GetSpecialResultIndex(DesignElement::DEFAULT, DesignElement::SHAPE_MAP_CORNER);

  // prepare all corner values
  EvalAllCornerValues();

  // statistics - be sure to handle them correctly when running parallel.
  // Set them to NumInt::int_cells_* after the parallel loop, omp just does not allow the class attributes
  int cells_cnt = 0;
  int cells_order_sum = 0;

  #pragma omp parallel num_threads(CFS_NUM_THREADS)
  {
    // these are thread local objects reused over the for loop iterations
    Vector<unsigned int> idx(dim_);
    StdVector<double>    ip(dim_); // the current ip within the element
    // num_node_shapes_ can be larger Item::nodes as we are not interested in 3D center second shapes.
    // for each shape of the Item this is order obtained form the cornver_val. 0 = void, 1 = solid, >=2 need integration
    Vector<int>          order(map_[0].nodes.GetSize());
    // this helps us evaluation
    EvalAtIp eval(this);

     // the integration effort is not evenly distributed
     // dynamic scheduling fails for gcc9 and is 4 times slower!!!
     #pragma omp for schedule(static) reduction(+:cells_cnt,cells_order_sum)
     for(int r = 0; r < (int) map_.GetSize(); r++)
     {
       Item& item = map_[r];
       DesignElement* de = item.elemval;

      DensityIdx(r, idx);

      // check what we need to integrate. The order is the maximal order of all relevant shapes
      // when we integrate we integrate all shapes at the same integration points.
      // If this could be relaxed we might be able to save quite some time!!
      int max_order = item.GetOrder(order, numInt_); // sets order. 0, 1 or >= 2

      if(res_idx_r >= 0)
        de->specialResult[res_idx_r] = max_order;

      if(res_idx_c >= 0)
        de->specialResult[res_idx_c] = item.MaxDiffCornerValue();

      double rho = -1.0; // indicator value!
      if(max_order == 0)
        rho = 0.0;
      if(combine_ == MAX && max_order ==1)
        rho = 1.0;
      // in the max sum case we need no integration if there is a 1.
      if(combine_ == MAX && max_order >= 2 && order.Contains(1))
        rho = 1.0;

      LOG_DBG2(SMD) << "MSTD: de=" << de->elem->elemNum << " mo=" << max_order << " o=" << order.ToString() << " rho=" << rho;

      // the number of integration points per dimension is actually max_order
      // but we need to exclude the coded max_order=1 when overlapping != MAX
      int num_ip = std::max(2, max_order);

      // we really need to integrate when we found no special case yet
      if(rho == -1.0)
      {
        cells_cnt++;
        cells_order_sum += num_ip;

        rho = 0.0; // such that we can sum the ip to it
        // it makes sense to traverse first the ip and then the variables
        for(int ip_x = 0; ip_x < num_ip; ip_x++)
        {
          for(int ip_y = 0; ip_y < num_ip; ip_y++)
          {
            for(int ip_z = 0; ip_z < (dim_ == 2 ? 1 : num_ip); ip_z++)
            {
              double ip_rho = 0.0; // the final value for one integration point. shall be not much larger one
              double weight = Item::SetIPGetWeight(this, ip, ip_x, ip_y, ip_z, num_ip);
              switch(combine_)
              {
              case MAX:
                assert(order.Max() >= 2);
                for(unsigned int si = 0; si < item.nodes.GetSize(); si++)
                {
                  assert(order[si] == 0 || order[si] >= 2);
                  if(order[si] >= 2) // otherwise it is 0 as 1 is checked above
                  {
                    eval.Setup(item.nodes[si], idx, ip, beta_);
                    double t = eval.ShapeFunc();
                    if(t >= ip_rho)  // >= is important! > may result in ip_eval == -1
                      ip_rho = t;
                    LOG_DBG3(SMD) << "MSTD: de=" << de->elem->elemNum << " ip=" << ip.ToString() << " si=" << si << " t=" << t << " max=" << ip_rho;
                  }
                }
                break;
              case TANH_SUM:
                // the original sum but with half beta
                for(unsigned int si = 0; si < item.nodes.GetSize(); si++)
                {
                  if(order[si] == 1)
                    ip_rho += 1.0;
                  if(order[si] >= 2){
                    eval.Setup(item.nodes[si], idx, ip, 0.5 * beta_); // half beta as it is applied to tanh_sum_.map()
                    ip_rho += eval.ShapeFunc();
                  }
                }
                // correct ip_rho by mapping <= 1. This is not exact, it might be slightly larger 1. See TanhSum()
                ip_rho = tanh_sum_.map(ip_rho);
                break;
              default:
                break;
              } // end of switch(combine_)
              assert(ip_rho >= 0 && ip_rho <= 1.02); // allow small overlap for TANH_SUM
              rho += weight * ip_rho; // apply weight only here at the end such that tanh_sum_.map() can be applied
              LOG_DBG3(SMD) << "MSTD: de=" << de->elem->elemNum << " ip=" << ip.ToString() << " mo=" << max_order << " ni=" << num_ip << " w=" << weight << " ip_rho=" << ip_rho << " -> " << rho;
            } // end ip_z
          } // end ip_y
        } // end ip_x
      } // end real integration
      assert(rho >= 0 && rho <= 1.02);
      de->SetDesign(de->GetLowerBound() + (de->GetUpperBound() - de->GetLowerBound()) * rho); // we assume 0 <= v <= 1
      LOG_DBG2(SMD) << "MS2D: -> el=" << de->elem->elemNum << " -> avg=" << de->GetPlainDesignValue()
                         << " delta=" << (de->GetPlainDesignValue() - de->GetLowerBound());
      assert(!(combine_ == MAX && de->GetPlainDesignValue() >= de->GetUpperBound() + 1e-10));
      assert(!(combine_ == TANH_SUM && de->GetPlainDesignValue() >= de->GetUpperBound() + 0.01)); // allow a slight overshot at overlap
      assert(de->GetPlainDesignValue() >= de->GetLowerBound() - 1e-10);
    } // end loop over density elements
  } // end of omp parallel section

  // omp just does not allow to reduce numInt.int_cells_cnt, ...
  numInt_.int_cells_cnt_ = cells_cnt;
  numInt_.int_cells_order_sum_ = cells_order_sum;
  numInt_.ToInfo(info_->Get("shapeMap/numInt"));
  mapped_design_ = design_id;
  mapping_timer_->Stop();
}

void ShapeMapDesign::MapFeatureGradient(const Function* f)
{
  assert(design_id == mapped_design_); // we need the Item setting from MapShapeDesign for the current design!
  assert(!(!f->IsObjective() && dynamic_cast<const Condition*>(f)->IsLocalCondition())); // it makes no sense for a local condition!!

  gradient_timer_->Start();

  // fixme! We do the job of dtanh_da for each function! However if we do it common
  // Optimization::EvalObjectiveConstraints() triggers MapShapeGradient() via WriteGradientToExtern() but rho::constraintGrad might not be set yet
  // the mapping goes shape->rho(->state), e.g. to compliance
  // from DesignSpace and ErsatzMaterial we have d_function/d_rho.
  // By chain rule we get d_function/d_shape by summing up for each shape_var the corresponding d_function/d_rho.
  // Note that we do this based on each integration point! It counts for each ip the shape_var which has larges rho(ip).
  // This was set in MapShapeToDensity() to Item::ip_param_id
  // ip_param_id is -1 if rho indicates zero gradient d_mapping/d_shape
  //
  // !! Also all simp gradients are not computed before the first WriteGradientToExtern() which triggers this MapShapeGradient

  // shall we store max_grad for output? -1 if not

  int res_idx_da = GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::SHAPE_MAP_GRAD, dim_ == 2 ? DesignElement::SM_NODE : DesignElement::SM_NODE_A);
  int res_idx_db = GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::SHAPE_MAP_GRAD, DesignElement::SM_NODE_B);
  int res_idx_dw = GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::SHAPE_MAP_GRAD, DesignElement::SM_PROFILE);

  // in the tanh_sum case the inner sum is constructed by 1/2* normal beta
  double beta = combine_ == TANH_SUM ? 0.5 * beta_ : beta_;


  bool node_grad = !IsAllNodeFixed();
  bool profile_grad = !IsProfileFixed();
  assert(node_grad || profile_grad);

  // to speed up performance and allow parallelization we have to not call BaseDesignElement::AddGradient()
  // within each integration point but have a flat vector which is added after the map loop
  // The variables a0,a1, b0, ... might be non-opt variables in the symmetry case as the symmetry mapping is done later
  Vector<double> shape_f_grad(shape_param_.GetSize());
  shape_f_grad.Init(0.0);

  #pragma omp parallel num_threads(CFS_NUM_THREADS)
  {
    // this are thread private constructs to be reused over the for loop iterations

    // as we have no OpenMP 4.5 but 3.1 we cannot use array reduction or reduction functions :(
    Vector<double> private_shape_f_grad(shape_f_grad); // the local array is initialized with 0
    assert(private_shape_f_grad.Max() == 0.0);

    // within the element coordinates we perform the integration
    Vector<unsigned int> idx(3);
    StdVector<double>    ip(dim_); // the current ip within the element

    Vector<int>          order(map_[0].nodes.GetSize()); // see MapShapeToDensity()
    StdVector<EvalAtIp>  eval(map_[0].nodes.GetSize());
    for(unsigned int i = 0; i < eval.GetSize(); i++)
      eval[i].Init(this);


     #pragma omp for schedule(static) // dynamic fails for gcc9 and is much slower!
     for(Integer r = 0; r < (Integer) map_.GetSize(); r++) // traverse all rho design elements
     {
       Item& item = map_[r];
       DesignElement* de = item.elemval;
       double log_da = 0.0;
       double log_db = 0.0;
       double log_dw = 0.0;

      int max_order = item.GetOrder(order, numInt_); // sets order. 0, 1 or >= 2

      LOG_DBG2(SMD) << "MSG: f=" << f->ToString() << " de=" << de->elem->elemNum << " rho=" << de->GetPlainDesignValue() << " mo=" << max_order << " o=" << order.ToString();

      DensityIdx(r, idx);

      // max_order = 0=void there is no gradient to add.
      // if max_order=1=solid and we have strict solid and the gradient is also zero.
      if((combine_ != MAX && max_order >= 1) || (combine_ == MAX && max_order >= 2 && !order.Contains(1)))
      {
        // case were we have a non-zero gradient!
        double de_plain_f_grad = de->GetPlainGradient(f);
        assert(!std::isnan(de_plain_f_grad));

        assert(max_order >= 2 || (max_order == 1 && combine_ == TANH_SUM));
        assert(!(combine_ == MAX && order.Contains(1)));

        // @see MapShapeToDensity()
        int num_ip = std::max(max_order, 2);

        for(int ip_x = 0; ip_x < num_ip; ip_x++)
        {
          for(int ip_y = 0; ip_y < num_ip; ip_y++)
          {
            for(int ip_z = 0; ip_z < (dim_ == 2 ? 1 : num_ip); ip_z++)
            {
              double weight = Item::SetIPGetWeight(this, ip, ip_x, ip_y, ip_z, num_ip);

              // find for MAX the idx of the shape with the maximal value
              // find for TANH_SUM the sum of all shape evals
              // we don't save this from MapShapeToDensity() as the number can be extremely large
              int max_idx = 0; // such it works also for TANH_SUM
              double eval_sum = 0;
              switch(combine_)
              {
              case MAX:
              {
                double max = -1;
                for(unsigned int s = 0; s < item.nodes.GetSize(); s++)
                {
                  eval[s].Setup(item.nodes[s], idx, ip, beta);
                  double t = eval[s].SmartShapeFunc(order[s]);
                  if(t > max) {
                    max = t;
                    max_idx = s;
                  }
                }
                break;
              }
              case TANH_SUM:
                for(unsigned int s = 0; s < item.nodes.GetSize(); s++)
                {
                  eval[s].Setup(item.nodes[s], idx, ip, beta);
                  eval_sum += eval[s].SmartShapeFunc(order[s]);
                }
                break;
              default:
                break;
              }

              double da = 0.0;
              double db = 0.0;
              double dw = 0.0;

              // in the MAX case we operate exactly on one shape, for TANH_SUM we loop over all
              assert(!(combine_ == TANH_SUM && max_idx != 0));
              for(unsigned int si = max_idx; si < (combine_ == MAX ? max_idx+1 : item.nodes.GetSize()); si++)
              {
                double t = eval[si].Setup(item.nodes[si], idx, ip, beta);
                assert(t >= 0 && t <= 1);
                assert(!(ip_x == 0 && ip_y == 0 && ip_z == 0 && t != 0));

                if(node_grad)
                {
                  da = eval[si].SmartGradShapeFunc(order[si], true, false, false); // dtanh_da
                  assert(!std::isnan(da) && !std::isinf(da));
                  LOG_DBG3(SMD) << "MSG: da d=" << de->GetIndex() << " p=" << de->GetLocation()->ToString() << " ip=" << ip.ToString() << " da=" << da;
                  if(dim_ == 3) // FIXME assumes center nodes!
                    db = eval[si].SmartGradShapeFunc(order[si], false, true, false); // dtanh_db

                  if(combine_ == TANH_SUM) { // tanh_sum shapes the sum approx to 0...1. sum is ip_eval
                    da = tanh_sum_.d_map(eval_sum, da);
                    if(dim_ == 3)
                      db = tanh_sum_.d_map(eval_sum, db);
                  }
                }
                if(profile_grad)
                {
                  dw = eval[si].SmartGradShapeFunc(order[si], false, false, true); // dtanh_dw
                  if(combine_ == TANH_SUM)
                    dw = tanh_sum_.d_map(eval_sum, dw);
                }

                // even if we have no node_grad, we need a0 for profile_grad
                ShapeMapVariable* a0 = item.nodes[si][0]; // due to symmetrie, this element is not necessarily a direct optimization variable
                ShapeMapVariable* a1 = item.nodes[si][dim_ == 3 ? 2 : 1];
                ShapeMapVariable* b0 = dim_ == 3 ? item.nodes[si][1] : NULL;
                ShapeMapVariable* b1 = dim_ == 3 ? item.nodes[si][3] : NULL;

                if(node_grad)
                {
                  double da_norm = (de->GetUpperBound() - de->GetLowerBound()) * da * weight;
                  double db_norm = (de->GetUpperBound() - de->GetLowerBound()) * db * weight;
                  log_da += da_norm;
                  log_db += db_norm;

                  private_shape_f_grad[a0->GetIndex()] += de_plain_f_grad * (1-t) * da_norm; // a0->AddGradient(f, de_plain_f_grad * (1-t) * da_norm);
                  private_shape_f_grad[a1->GetIndex()] += de_plain_f_grad * t * da_norm;
                  if(dim_ == 3) {
                    private_shape_f_grad[b0->GetIndex()] += de_plain_f_grad * (1-t) * db_norm;
                    private_shape_f_grad[b1->GetIndex()] += de_plain_f_grad * t * db_norm;
                  }
                }
                if(profile_grad)
                {
                  // a and b share common w
                  double dw_norm = (de->GetUpperBound() - de->GetLowerBound()) * dw * weight;

                  log_dw += dw_norm;
                  private_shape_f_grad[GetProfile(a0)->GetIndex()] += de_plain_f_grad * (1-t) * dw_norm;
                  private_shape_f_grad[GetProfile(a1)->GetIndex()] += de_plain_f_grad * t * dw_norm;
                  LOG_DBG3(SMD) << "MSG: prof de=" << de->elem->elemNum
                      << " p0=" << GetProfile(item.nodes[si][0])->GetIndex()
                      << " p1=" << GetProfile(item.nodes[si][1])->GetIndex()
                      << " s0=" << item.nodes[si][0]->GetIndex() << " s1=" << item.nodes[si][1]->GetIndex()
                      << " ip=" << ip.ToString() << " t=" << t << " dwn=" << dw_norm;
                }

                LOG_DBG3(SMD) << "MSG: el=" << de->elem->elemNum << " ip=" << ip.ToString() << " da=" << da << " dw=" << dw
                    << " fg=" << de_plain_f_grad << " da0=" << shape_f_grad[a0->GetIndex()] << " da1=" << shape_f_grad[a1->GetIndex()]
                                                                                                                       << " dw0=" << shape_f_grad[GetProfile(a0)->GetIndex()] << " dw1=" << shape_f_grad[GetProfile(a1)->GetIndex()];
              } // shape loop
            } // end ip_z
          } // end ip_y
        } // end ip_x
      } // normalize by integration points.
      LOG_DBG2(SMD) << "MSG: el=" << de->elem->elemNum << " rho=" << de->GetPlainDesignValue() << " sum da=" << log_da << " sum dw=" << log_dw;
      if(res_idx_da >= 0)
        de->specialResult[res_idx_da] = log_da;
      if(res_idx_db >= 0)
        de->specialResult[res_idx_db] = log_db;
      if(res_idx_dw >= 0)
        de->specialResult[res_idx_dw] = log_dw;
    } // end loop over density elements

    // now gather private_shape_f_grad to shape_f_grad as we may not use array reduction yet
    // https://stackoverflow.com/questions/20413995/reducing-on-array-in-openmp
    #pragma omp critical
    {
      // we are in the omp parallel section with n threads. The critical only of the n threads at one time
      // but this is executed for all.
      shape_f_grad.Add(private_shape_f_grad);
    }
  } // end of omp parallel

  // write back shape_f_grad
  LOG_DBG3(SMD) << "MSG: f=" << f->ToString() << " sfg=" << shape_f_grad.ToString();

  assert(shape_f_grad.GetSize() == shape_param_.GetSize());
  for(unsigned int i = 0; i < shape_f_grad.GetSize(); i++)
    if(shape_f_grad[i] != 0) // if it is not an opt variable DesignElement::*Gradient has size zero. Note negative gradients!
      shape_param_[i]->AddGradient(f, shape_f_grad[i]);

  gradient_timer_->Stop();
}

void ShapeMapDesign::WriteGradientFile()
{
  // plot the stuff like this:
  // plot "shape_map_mech.grad.dat" u 1:6 every ::0::40 w lp, "shape_map_mech.grad.dat" u ($1-41):6 every ::41::82 w lp

  if(!gradplot_.is_open())
    return; // obviously the option was not set

  // for shape map only the state of a single iteration makes sense, hence we overwrite for every iteration
  gradplot_.seekp(0);

  gradplot_.precision(8);
  gradplot_.flags(std::ios::scientific);

  assert(opt_->objectives.data.GetSize() == 1);
  gradplot_ << "# iteration: " << opt_->GetCurrentIteration() << std::endl;
  gradplot_ << "# gnuplot: plot \"" << progOpts->GetSimName() << ".grad.dat\"  u 1:6 every ::0::" << (shape_[0].end_opt-1) << " w lp";
  gradplot_ << ", \"" + progOpts->GetSimName() << ".grad.dat\"  u ($1-" << shape_[0].end_opt << "):6 every ::"
      << shape_[0].end_opt << "::" << (2*shape_[0].end_opt) << " w lp" << std::endl;
  gradplot_ << "#(1) el \t(2) var \t(3) shape \t(4) dof \t(5) val \t(6) " << opt_->objectives.data[0]->ToString();

  int cnt = 6; // will be preincremented
  for(unsigned int g = 0; g < opt_->constraints.all.GetSize(); g++)
    if(opt_->constraints.all[g]->HasDenseJacobian())
      gradplot_ << " ("  << lexical_cast<string>(++cnt) << ") " + ToValidXML(opt_->constraints.all[g]->ToString()) + "\t";
  gradplot_ << std::endl;

  Function* c = opt_->objectives.data[0];
  for(unsigned int e = 0; e < opt_shape_param_.GetSize(); e++)
  {
    ShapeParamElement* spe = opt_shape_param_[e];
    ShapeParam* shape = GetShape(spe);
    gradplot_ << e << " \t" << spe->type.ToString(spe->GetType()) << " \t" << shape->idx << " \t";
    gradplot_ << spe->dof_ << " \t" << spe->GetPlainDesignValue() << " \t" << (spe->GetPlainGradient(c) + dynamic_cast<ShapeMapVariable*>(opt_shape_param_[e])->sym->GetPlainSymGradient(c)) << " \t";

    for(unsigned int g = 0; g < spe->constraintGradient.GetSize(); g++)
      if(opt_->constraints.all[g]->HasDenseJacobian())
        gradplot_ << spe->constraintGradient[g] << " \t";
    gradplot_ << std::endl;
  }
  // trust on the destructor to close the file
}


void ShapeMapDesign::ExportLevelSet() const
{
  // this is still a dump implementation
  assert(export_leveset_);

  std::ofstream out;
  string name = progOpts->GetSimName() + ".leveset.vtk";
  out.open(name.c_str());
  out.precision(8);
  out.flags(std::ios::scientific);

  out << "# vtk DataFile Version 3.0\n";
  out << "Para0\n";
  out << "ASCII\n";
  out << "DATASET RECTILINEAR_GRID\n";
  assert(dim_ == 2);
  out << "DIMENSIONS " << (nx_+1) << " " << (ny_+1) << " 1\n";
  out << "X_COORDINATES " << (nx_+1) << " int\n";
  for(unsigned int i = 0; i <= nx_; i++)
    out << i << " ";
  out << "\nY_COORDINATES " << (ny_+1) << " int\n";
  for(unsigned int i = 0; i <= ny_; i++)
    out << i << " ";
  out << "\nZ_COORDINATES 1 int\n0\n\n";
  out << "POINT_DATA " << ((nx_+1) * (ny_+1) * 1) << "\n";
  out << "SCALARS lsf_0 float 1\n";
  out << "LOOKUP_TABLE default\n";
  // first shot simply takes the density, cooler options are nodes and Gaussian points
  // write the upper right corner with special handling for -1
  // lowest y-line first
  out << DensityToLevelSet(0,0) << std::endl;
  for(int x = 0; x < (int) nx_; x++)
    out << DensityToLevelSet(x,0) << std::endl;
  for(int y = 0; y < (int) ny_; y++)
  {
    out << DensityToLevelSet(0,y) << std::endl;
    for(int x = 0; x < (int) nx_; x++)
      out << DensityToLevelSet(x,y) << std::endl;
  }
  out.close();
}

inline double ShapeMapDesign::DensityToLevelSet(int x, int y) const
{
  return 2.0 * map_[DensityIdx(x,y)].elemval->GetPlainDesignValue() - 1.0;
}

void ShapeMapDesign::EvalAtIp::Init(ShapeMapDesign* smd)
{
  smd_ = smd;
  dim = domain->GetDim();
  coord_.Resize(dim, -1.0);

  // calc h which we need only for shapeFunc == LINEAR
  if(smd->GetBoundary() == LINEAR)
  {
    Matrix<double> box = domain->GetGrid()->CalcGridBoundingBox(NULL, true); // force 3d (0 size for z)
    Vector<unsigned int> n = smd->GetDiscretization();
    assert(n.GetSize() == box.GetNumRows());
    assert(box.GetNumCols() == 2); // min and max for every dim
    assert(n.GetSize() == 3);
    assert(n.Min() >= 1);

    Vector<double> spacing(3);
    for(unsigned int i = 0; i < 3; i++)
      spacing[i] = (box[i][1] - box[i][0]) / n[i];

    h = spacing.Max();
    assert(h > 0);
  }
}

inline double ShapeMapDesign::EvalAtIp::Setup(const StdVector<ShapeMapVariable*>& nodes, const Vector<unsigned int>& idx, const StdVector<double>& ip, double beta)
{
  assert(smd_ != NULL);

  if(dim == 2)
    return Setup2d(nodes, idx, ip, beta);
  else
    return Setup3d(nodes, idx, ip, beta);

}

inline double ShapeMapDesign::EvalAtIp::Setup2d(const StdVector<ShapeMapVariable*>& nodes, const Vector<unsigned int>& idx, const StdVector<double>& ip, double beta)
{
  assert(dim == 2); // there is also a 3D version of Eval
  assert(nodes.GetSize() == 2);
  // this stuff is 3D stuff
  const ShapeMapVariable* s1 = nodes[0];
  const ShapeMapVariable* s2 = nodes[1];
  assert(s1->dof_ == s2->dof_);
  assert(s1->GetType() == BaseDesignElement::NODE && s2->GetType() == BaseDesignElement::NODE);

  ShapeParamElement::Dof dof = s1->dof_;
  ShapeParamElement::Dof dir = Flip(dof);
  assert(dof != dir);

  smd_->MapIdxToCoords(idx, coord_);
  double start = coord_[s1->dof_];

  // the parameters
  double a1 = s1->GetPlainDesignValue();
  double a2 = s2->GetPlainDesignValue();
  // the profiles
  assert(smd_->GetProfile(s1)->GetType() == BaseDesignElement::PROFILE);
  double w1 = smd_->GetProfile(s1)->GetPlainDesignValue();
  double w2 = smd_->GetProfile(s2)->GetPlainDesignValue();

  // when dof = X the tanh has x as parameter, when we interpolate a1 and a2 these are x values applied at different y positions
  assert(ip.GetSize() == 2);
  assert(ip[0] >= 0 && ip[0] <= 1 && ip[1] >= 0 && ip[1] <= 1);

  x = start + ip[dof] * smd_->coord_step_[s1->dof_]; // for dof=0 (x) xy is x and might be far away from a

  // a = (1-t)*a1 + t*a2 = a1 - t*a1 + t*a2 = a1+t*(a2-a1): t=0 -> a1, t=1 -> a2
  t = ip[dir];
  a  = a1 + t * (a2-a1);         // for dof=1 (c) a is y and with a1=a2 we have the same value for a
  w  = w1 + t * (w2-w1);

  if(smd_->boundary_ == TANH) {
    assert(beta > 0);
    this->beta = beta;
    exapw = exp(beta*(x-a+w));
    examw = exp(beta*(x-a-w));
  }
  assert(b == -1);
  assert(y == -1);
  assert(r == -1);
  assert(erw == -1);

  return t;
}


inline double ShapeMapDesign::EvalAtIp::Setup3d(const StdVector<ShapeMapVariable*>& nodes, const Vector<unsigned int>& idx, const StdVector<double>& ip, double beta)
{
  assert(smd_ != NULL);
  assert(dim == 3);
  assert(nodes.GetSize() == 4);
  assert(ip.GetSize() == 3);
  assert(ip[0] >= 0 && ip[0] <= 1 && ip[1] >= 0 && ip[1] <= 1 && ip[2] >= 0 && ip[2] <= 1);

  const ShapeMapVariable* sa1 = nodes[0];
  const ShapeMapVariable* sb1 = nodes[1];
  const ShapeMapVariable* sa2 = nodes[2];
  const ShapeMapVariable* sb2 = nodes[3];

  LOG_DBG3(SMD) << "EAI:S3 sa1=" << sa1->ToString();

  assert(sa1->dof_ == sa2->dof_ && sb1->dof_ == sb2->dof_ && sa1->dof_ != sb1->dof_);
  assert(sa1->GetType() == sa2->GetType() && sa1->GetType() == sb1->GetType() && sa1->GetType() == sb2->GetType() && sa1->GetType() == BaseDesignElement::NODE);
  assert(smd_->GetShape(sa1)->IsCenterShape() && smd_->GetShape(sa1)->type == NODE);

  assert(smd_->GetProfile(sa1) == smd_->GetProfile(sb1)); // two center nodes share a profile

  // the parameters
  double a1 = sa1->GetPlainDesignValue();
  double a2 = sa2->GetPlainDesignValue();
  double b1 = sb1->GetPlainDesignValue();
  double b2 = sb2->GetPlainDesignValue();
  double w1 = smd_->GetProfile(sa1)->GetPlainDesignValue(); // a and b share profile
  double w2 = smd_->GetProfile(sa2)->GetPlainDesignValue();
  assert(!std::isnan(a1));
  assert(!std::isnan(smd_->GetProfile(sa1)->GetPlainDesignValue()));
  assert(smd_->GetProfile(sa1)->GetPlainDesignValue() == smd_->GetProfile(sb1)->GetPlainDesignValue());
  assert(smd_->GetProfile(sa2)->GetPlainDesignValue() == smd_->GetProfile(sb2)->GetPlainDesignValue());

  // sa and sb have different dof and define an ab-plane.
  // the complementary dof dir is where we interpolate within start and end
  smd_->MapIdxToCoords(idx, coord_);
  double start_x = coord_[sa1->dof_];
  double start_y = coord_[sb1->dof_];

  // we are in the ab-plane and call it xy-plane. We test of the point (x,y) on the xy-plane.
  int dir = Flip(sa1->dof_, sb1->dof_); // orthogonal to the xy-plane. E.g. the z-axis

  x = start_x + ip[sa1->dof_] * smd_->coord_step_[sa1->dof_];;
  y = start_y + ip[sb1->dof_] * smd_->coord_step_[sb1->dof_];

  t = ip[dir]; // Setup2d()
  assert(t >= 0 && t <= 1);

  a = a1 + t * (a2 - a1);
  b = b1 + t * (b2 - b1);
  w = w1 + t * (w2 - w1);

  r = sqrt((a-x)*(a-x)+(b-y)*(b-y)); // note that r can be 0!! Important for gradients!!

  if(smd_->GetBoundary() == TANH) // for 3D LINEAR we just need r
  {
    this->beta = beta;
    erw = exp(beta * (r-w));
  }

  assert(examw == -1);
  assert(exapw == -1);

  return t;
}

inline double ShapeMapDesign::EvalAtIp::ShapeFunc() const
{
  if(smd_->GetBoundary() == TANH)
    return dim == 2 ? EvalTanh2d() : EvalTanh3d();
  else
    return dim == 2 ? EvalLinear2d() : EvalLinear3d();
}

inline double ShapeMapDesign::EvalAtIp::GradShapeFunc(bool grad_a, bool grad_b, bool grad_w) const
{
  if(smd_->GetBoundary() == TANH)
    return dim == 2 ? EvalTanhGrad2d(grad_a, grad_w) : EvalTanhGrad3d(grad_a, grad_b, grad_w);
  else
    return dim == 2 ? EvalLinearGrad2d(grad_a, grad_w) : EvalLinearGrad3d(grad_a, grad_b, grad_w);
}

inline double ShapeMapDesign::EvalAtIp::SmartShapeFunc(int order) const
{
  assert(order >= 0);

  if(order == 0)
    return 0.0;
  if(order == 1)
    return 1.0;
  return ShapeFunc();
}


inline double ShapeMapDesign::EvalAtIp::EvalTanh2d() const
{
  // set xrange[0:1]; a = 0.5; w=0.1; beta=30
  // plot 1-1/(exp(beta*(x-a+w)) + 1), 1/(exp(beta*(x-a-w)) + 1)
  //
  // ta(x)=1-1/(exp(beta*(x-a+w)) + 1)
  if(x <= a)
    return 1.0 - 1/(exapw+1);
  else
    return 1/(examw+1);
}

inline double ShapeMapDesign::EvalAtIp::SmartGradShapeFunc(int order, bool grad_a, bool grad_b, bool grad_w) const
{
  assert(order >= 0);

  if(order == 0 || order == 1)
    return 0.0;
  return GradShapeFunc(grad_a, grad_b, grad_w);
}

double ShapeMapDesign::EvalAtIp::EvalTanhGrad2d(bool grad_a, bool grad_w) const
{
  // this is d_tanh_da
  // set xrange[0:1]; a = 0.5; beta=30; w=0.1
  // plot -1* (exp(beta*(x-a+w)) + 1)**-2 * beta*exp(beta*(x-a+w)), (exp(beta*(x-a-w))+1)**-2 beta*exp(beta*(x-a-w))
  //
  // ta(x)=1-1/(exp(beta*(x-a+w)) + 1)
  //
  // plot
  // matlab:
  //if x <= a
  //  f = -1* (exp(beta*(x-a+d)) + 1)^-2 *beta*exp(beta*(x-a+d));
  //else
  ///  f = (exp(beta*(x-a-d))+1)^-2 *2*beta*exp(beta*(x-a-d));
  // end

  // the difference between da and dw is only that da for x <= a has the factor -1. all other cases are without -1
  if(x <= a)
    return (grad_a ? -1.0 : 1.0) * 1.0/((exapw+1) * (exapw+1)) * beta*exapw;
  else
    return 1/((examw+1) * (examw+1)) * beta*examw;
}

inline double ShapeMapDesign::EvalAtIp::EvalTanh3d() const
{
  // a:(1-t)*a1+t*a2;
  // b:(1-t)*b1+t*b2;
  // r: sqrt((a-x)^2+(b-y)^2);
  // t:1/(exp(beta*(r-w))+1);
  return 1/(erw +1);
}

inline double ShapeMapDesign::EvalAtIp::EvalTanhGrad3d(bool grad_a, bool grad_b, bool grad_w) const
{
  assert(!std::isnan(erw));
  assert(!std::isinf(erw));
  assert(erw != -1);
  assert(erw > 0);
  assert(beta != -1);
  // r: sqrt((a-x)^2+(b-y)^2); might be zero, see below!
  // e = exp(beta * (r-w));
  // t:1/(exp(beta*(r-w))+1);

  // beta=100; w=.2; set xrange[0:1]; set yrange[0:1]
  // splot 1/(exp(beta*(sqrt((.5-x)**2+(.5-y)**2)-w))+1

  // da = -1 * beta*(a-x) * e / (r * (e+1)*(e+1))
  // db = -1 * beta*(b-x) * e / (r * (e+1)*(e+1))
  // dw = beta * e / ((e+1)*(e+1))

  // note! da computes to 1/0 for r = 0 -> x=a and y=b!!
  // from gnuplut we see, that da(r=0) = 0, this needs to be reflected in implementation!

  // wolframalpha
  // differentiate 1/(exp(beta*(sqrt((a-x)^2+(b-y)^2)-w))+1) wrt a

  // gnuplot
  // a = 0.5; b = 0.5; beta=100; w=.2; set xrange[0:1]; set yrange[0:1]; set isosample 40
  // r(x,y)=sqrt((a-x)**2+(b-y)**2)
  // e(x,y)=exp(beta * (r(x,y)-w))
  // f(x,y)=1/(exp(beta*(r(x,y)-w))+1)
  // da(x,y)=-1 * beta*(a-x) * e(x,y) / (r(x,y) * (e(x,y)+1)**2)
  // db(x,y)=-1 * beta*(b-y) * e(x,y) / (r(x,y) * (e(x,y)+1)**2)
  // dw(x,y)=beta * e(x,y) / (e(x,y)+1)**2
  // splot da(x,y)

  if(grad_w)
    return beta * erw / ((erw+1)*(erw+1));

  // see comment above about 1/0 case!
  if(r < 1e-13)
    return 0;

  assert(grad_a || grad_b);
  double grad = grad_a ? -1 * beta*(a-x) * erw / (r * (erw+1)*(erw+1)) :
      -1 * beta*(b-y) * erw / (r * (erw+1)*(erw+1));

  LOG_DBG3(SMD) << "ETG3: da=" << grad_a << " a=" << a << " b=" << b << " x=" << x << " y=" << y << " r=" << r << " erw=" << erw << " -> " << grad;
  return grad;
}

inline double ShapeMapDesign::EvalAtIp::EvalLinear2d() const
{
  //         __________
  //        /          \
  //       /            \
  //  ____/              \_____
  //     a0  a1       b1 b0

  // a1 -> b1 = 2*w + h
  // a0 -> a1 = h

  double a0 = a - w - h/2;
  double b1 = a + w - h/2;

  if(x < a0)
    return 0;
  if(x < a0 + h)
    return (x-a0)/h;
  if(x < b1)
    return 1;
  if(x < b1+h)
    return 1-(x-b1)/h;
  return 0;
}

inline double ShapeMapDesign::EvalAtIp::EvalLinearGrad2d(bool grad_a, bool grad_w) const
{
  double a0 = a - w - h/2;
  double b1 = a + w - h/2;

  if(x < a0)
    return 0;
  if(x < a0 + h)
    return grad_a ? -1/h : 1/h;; // spacing h means for dx=h dy=h. Same value for grad_a and grad_w
  if(x < b1)
    return 0;
  if(x < b1+h)
    return 1/h; // grad_a and grad_w

  return 0;
}

inline double ShapeMapDesign::EvalAtIp::EvalLinear3d() const
{
  // similar to 2d but with radius r
  // |_____
  // |     \
  // |      \
  // |       \_____
  // |    r0  r1

  double r0 = w - h/2;

  if(r < r0)
    return 1;
  if(r < r0 + h)
    return 1-(r-r0)/h;
  return 0; // x >= r1
}

inline double ShapeMapDesign::EvalAtIp::EvalLinearGrad3d(bool grad_a, bool grad_b, bool grad_w) const
{
  // r = sqrt((a-x)*(a-x)+(b-y)*(b-y))
  // r0= w - h/2
  // f = 1-(r-r0)/h
  // f = 1-(sqrt((a-x)^2+(b-y)^2)-(w - h/2))/h
  // wolframalpha: diff 1-(sqrt((a-x)^2+(b-y)^2)-(w - h/2))/h by a

  double r0 = w - h/2;
  // outside the gray region 0
  if(r < r0 || r >= r0 + h)
    return 0;

  if(grad_w)
    return 1/h;

  // see EvalTanhGrad3d()
  if(r < 1e-13)
    return 0;

  if(grad_a)
    return -(a-x) / (h*r);
  if(grad_b)
    return -(b-y) / (h*r);

  assert(false);
  return -1;
}


void ShapeMapDesign::EvalAllCornerValues()
{
  // for a nx_*ny_*nz_ mesh we evaluate (nx_+1)*(ny_+1)*(nz_+1) points for each shape and set each value for all adjacent Item. Hence the value is repeated almost
  // 4 times in 2D and 8 times in 3D
  assert(dim_ == 3 || nz_ == 1);

  StdVector<double> ip(dim_); // integration location 0 ... 1
  Vector<unsigned int> idx(3);

  int num_nodes = map_[0].nodes.GetSize();
  assert((unsigned int) num_nodes == num_node_shapes_ - FindCenters().GetSize());

  // temporary array. Could be used permanent instead of replicating it 4 to 8 times in Item::corner_vals
  StdVector<Vector<double> > glob(num_nodes);
  for(int s = 0; s < num_nodes; s++)
    glob[s].Resize((nx_+1) * (ny_+1) * (dim_ == 3 ? nz_+1 : 1), -1.0);

  // we index by z * (nx_+1)*(ny_+1) + y * (nx_+1) + x
  int zb = (dim_ == 3 ? (nx_+1)*(ny_+1) : 0);
  int yb = nx_ + 1;
  // z * zb + y * yb + x for 3D and 2D

  EvalAtIp eval(this);

  for(unsigned int z = 0; z < nz_+1; z++)
  {
    for(unsigned int y = 0; y < ny_+1; y++)
    {
      for(unsigned int x = 0; x < nx_+1; x++)
      {
        // normaly we evaluate the smallest corner. For the last elements we use the outer corner of the element before
        idx[0] = x < nx_ ? x : x-1;
        idx[1] = y < ny_ ? y : y-1;
        idx[2] = z < nz_ ? z : z-1;

        ip[0] =  x < nx_ ? 0 : 1.0;
        ip[1] =  y < ny_ ? 0 : 1.0;
        if(dim_ == 3)
          ip[2] =  z < nz_ ? 0 : 1.0;

        Item& item = map_[DensityIdx(idx[0], idx[1], idx[2])]; // ignores z in 2D

        for(int s = 0; s < num_nodes; s++)
        {
          assert(item.nodes[s].GetSize() >= 2);

          eval.Setup(item.nodes[s], idx, ip, 2 * beta_);
          glob[s][z * zb + y * yb + x] = eval.ShapeFunc();
        } // end shape loop
      } // x
    } // y
  } // z


  //std::cout << glob[0].ToString() << std::endl;
  assert(glob[0].Min() >= 0.0);

  // we remap from glob to vals
  Vector<double> vals(dim_ == 2 ? 4 : 8);

  for(unsigned int z = 0; z < nz_; z++)
  {
    for(unsigned int y = 0; y < ny_; y++)
    {
      for(unsigned int x = 0; x < nx_; x++)
      {
        Item& item = map_[DensityIdx(x, y, z)]; // ignores z in 2D

        // NOTE: the 3D mapping is different form the coords definition of CFS as used in Eval() !
        // 2D -------------      : 3D --------- coords is 8 columns with 3 rows
        //
        //  a=0, b=1, c=2, d=3   :  a=0  b=1  c=2  d=3  e=4  f=5  g=6  h=7
        //
        //                       :   d----------c
        //                       :  /|         /|
        //  d ------- c          : h ------- g  |
        //  |         |          : | a--------|-b
        //  |         |          : |/         |/
        //  a ------- b          : e ------- f
        //  x=a->b, y=a->d       : x=a->b, y=a->d, z=a->e, zero = a

        for(int s = 0; s < num_nodes; s++)
        {
          assert(item.nodes[s].GetSize() >= 2);

          vals[0] = glob[s][z * zb + y * yb + x];          // a
          vals[1] = glob[s][z * zb + y * yb + (x+1)];      // b
          vals[2] = glob[s][z * zb + (y+1) * yb + (x+1)];  // c
          vals[3] = glob[s][z * zb + (y+1) * yb + x];     // d
          if(dim_ == 3) {
            vals[4] = glob[s][(z+1) * zb + y * yb + x];          // e
            vals[5] = glob[s][(z+1) * zb + y * yb + (x+1)];      // f
            vals[6] = glob[s][(z+1) * zb + (y+1) * yb + (x+1)];  // g
            vals[7] = glob[s][(z+1) * zb + (y+1) * yb + x];     // h
          }
          item.min_corner_value[s] = vals.Min();
          item.max_corner_value[s] = vals.Max();
          LOG_DBG3(SMD) << "EACV: x=" << x << " y=" << " z=" << z << " s=" << s
              << " min=" << item.min_corner_value[s] << " max=" << item.max_corner_value[s] << " -> " << vals.ToString();
          assert(vals.Min() >= 0.0);
        } //s
      } // x
    } // y
  } // z

  // DumpMap();
}



inline ShapeParamElement::Dof ShapeMapDesign::Flip(ShapeParamElement::Dof dof)
{
  assert(dof != ShapeParamElement::NOT_SET);
  assert(dim_ == 2);
  assert(dof == ShapeParamElement::X || dof == ShapeParamElement::Y);
  return dof == ShapeParamElement::X ? ShapeParamElement::Y : ShapeParamElement::X;
}


inline ShapeParamElement::Dof ShapeMapDesign::Flip(ShapeParamElement::Dof first, ShapeParamElement::Dof second)
{
  assert(first == ShapeParamElement::X || first == ShapeParamElement::Y || first == ShapeParamElement::Z);
  assert(second == ShapeParamElement::X || second == ShapeParamElement::Y || second == ShapeParamElement::Z);
  assert(first != second);
  switch(first) {
  case ShapeParamElement::X:
    return second == ShapeParamElement::Y ? ShapeParamElement::Z : ShapeParamElement::Y;
  case ShapeParamElement::Y:
    return second == ShapeParamElement::Z ? ShapeParamElement::X : ShapeParamElement::Z;
  case ShapeParamElement::Z:
    return second == ShapeParamElement::X ? ShapeParamElement::Y : ShapeParamElement::X;
  default:
    assert(false);
  }
  return ShapeParamElement::NOT_SET; // must not happen
}

void ShapeMapDesign::DumpMap()
{
  for(unsigned int z = 0; z < nz_; z++)
  {
    for(unsigned int y = 0; y < ny_; y++)
    {
      for(unsigned int x = 0; x < nx_; x++)
      {
        Item i = map_[DensityIdx(x, y, z)];
        std::cout << "z=" << z << " y=" << y << " x=" << x << " elidx=" << DensityIdx(x, y, z) << " rho=" << i.elemval->GetPlainDesignValue()<< std::endl;
        // " min=" << i.min_corner_value << " max=" << i.max_corner_value << std::endl;
        for(unsigned int n = 0; n < i.nodes.GetSize(); n++)
          for(unsigned int q = 0; q < i.nodes[n].GetSize(); q++)
            std::cout << "-> n=" << n << " q=" << q << " si=" << i.nodes[n][q]->GetIndex() << " dof=" << i.nodes[n][q]->dof_ << " p=" << GetProfile(i.nodes[n][q])->GetIndex() << std::endl;
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }
  }
}

void ShapeMapDesign::CreateShapeVariable(const ShapeParam* param,  ShapeParamElement::Dof free_dof, int free_idx, bool start_end)
{
  // note that free corresponds to the node counter, not element counter as max free is ny_ and not ny_-1

  // this is one example of a shape param with dof=x, value= 4 and free 0...5 (y)
  // We have 5 elements E0...E4 and 6 nodes 0..5
  //  y (node)
  //  5              M              E4
  //  4              L           E3 E4
  //  3              K        E2 E3
  //  2              Z     E1 E2
  //  1              Y  E0 E1
  //  0              X  E0
  //  x  0  1  2  3  4  5  6  7  8  9 (node)
  //
  // E0(X,Y), E1(Y,Z), E2(Z,K), E3(K,L), E4(L,M)
  //
  // X has free=0 and is for the 9 elements 0..8 with y=0 and y=1 (E0)
  // Y has free=1 and is for the 9 elements 0..8 with y=1 and y=2 (E0 and E1)
  // ...
  // L has free=4 and is for the 9 elements 0..8 with y=3 and y=4 (E3 and E4)
  // M has free=4 and is for the 9 elements 0..8 with y=4 and y=5 (E4)


  assert((int) free_dof >= 0 && (int) free_dof <= 2);
  assert(free_dof == ShapeParamElement::X || free_dof == ShapeParamElement::Y || free_dof == ShapeParamElement::Z); // equivalent to the assert above
  assert(free_dof >= 0);
  assert(free_dof <= (int) n_[(int) free_dof]);

  // 3D cubical symmetry has one center node as slave. Then we need for fixed, initial, lower, upper, clamp the data from the master shape
  // but we still need the original param dof at the end
  assert(!param->slave || param->other_center != NULL);
  const ShapeParam* data = param->slave ? param->other_center : param;
  assert(param->type == data->type);

  // add element to shape_param_
  assert(shape_param_.GetCapacity() > shape_param_.GetSize());
  ShapeMapVariable* temp = new ShapeMapVariable(Convert(data->type), shape_param_.GetSize());
  shape_param_.Push_back(temp);
  ShapeMapVariable* spe = dynamic_cast<ShapeMapVariable*>(shape_param_.Last());

  MathParser* mp = domain->GetMathParser();
  unsigned int handle = mp->GetNewHandle();

  // the value might be a formula like "0.5*xi/nx", to evaluate this we set xi to the free variable. Note: in 2D for dof=x the free_idx variable are the y-nodes
  // TODO this has limitations in 3D with surfacs
  std::string var = ShapeParamElement::dof.ToString(free_dof) + "i";
  assert(var == "xi" || var == "yi" || var == "zi");
  // set the variable which might be in the formula
  mp->SetValue(handle, var, free_idx);

  // this might contain a formula or simply a value
  mp->SetExpr(handle, data->value); // TODO: make faster by doing this outside the look

  double value = mp->Eval(handle);

  double lower = -1;
  double upper = -1;

  if(!data->fixed)
  {
    mp->SetExpr(handle, data->lower);
    lower = mp->Eval(handle);

    mp->SetExpr(handle, data->upper);
    upper = mp->Eval(handle);
  }

  spe->SetDesign(value);
  if(data->clamp >= 0.0 && start_end)
  {
    spe->SetLowerBound(value - data->clamp/2);
    spe->SetUpperBound(value + data->clamp/2);
    LOG_DBG(SMD) << "CSV el=" << (shape_param_.GetSize() - 1) << " shape=" << param->ToString() << " data=" << data->ToString() << " clamped! lb=" << spe->GetLowerBound() << " ub=" << spe->GetUpperBound();
  }
  else
  {
    double rb = spe->GetType() == DesignElement::NODE ? relative_node_bound_ : relative_profile_bound_;
    if(rb >= 0 && !DensityFile::NeedLoadErsatzMaterial()) // don't set the bounds relative to initial when we later load an external design
    {
      spe->SetUpperBound(std::min(upper, std::max(value + rb, lower)));
      spe->SetLowerBound(std::max(lower, std::min(value - rb, upper)));
    }
    else
    {
      spe->SetLowerBound(lower);
      spe->SetUpperBound(upper);
    }
  }

  spe->dof_ = param->dof; // here read param, not data!
  spe->coord[(int) free_dof] = (double) free_idx / n_[(int) free_dof];
  spe->idx[(int) free_dof] = free_idx;

  mp->ReleaseHandle(handle);
  // PostInit() sets arrays for objective and constraint gradients

  LOG_DBG2(SMD)<< "CSV el=" << (shape_param_.GetSize() - 1) << " dof=" << spe->dof_ << " free_idx=" << free_idx << " d=" << spe->GetPlainDesignValue() << " coord=" << spe->coord.ToString();
}

void ShapeMapDesign::PostInit(int objectives, int constraints)
{
  // add aux_design to full_data
  FeaturedDesign::PostInit(objectives, constraints);

  for(ShapeParamElement* opt : opt_shape_param_)
    dynamic_cast<ShapeMapVariable*>(opt)->sym->PostInit(objectives, constraints);
}

StdVector<ShapeMapDesign::ShapeParam*> ShapeMapDesign::FindShape(Type type, ShapeParamElement::Dof dof)
{
  StdVector<ShapeParam*> res;
  for(unsigned int i = 0; i < shape_.GetSize(); i++) // could be smarter bur more complex if we consider the type
    if(shape_[i].type == type && shape_[i].dof == dof)
      res.Push_back(&shape_[i]);
  return res;
}

ShapeMapDesign::ShapeMapVariable* ShapeMapDesign::GetSecondCenterNodeParam(ShapeParam* shape, ShapeMapVariable* test)
{
  assert(shape->IsCenterShape());
  assert(shape->type == NODE);
  unsigned int test_idx = test->GetIndex();
  assert(shape->IsPart(test) || shape->other_center->IsPart(test));

  // we assume that test is a first center node
  if(shape->IsFirstCenterNode() && shape->IsPart(test))
    return dynamic_cast<ShapeMapVariable*>(shape_param_[shape->GetSecondCenterNode()->start_param + test_idx - shape->start_param]);

  // obviously test is already a second center node
  assert(shape->GetSecondCenterNode()->IsPart(test));
  return test;
}


inline const ShapeMapDesign::ShapeMapVariable* ShapeMapDesign::GetProfile(const ShapeMapVariable* node) const
{
  const ShapeParam* shape = GetShape(node);
  assert(node->GetType() == ShapeParamElement::NODE && shape->type == NODE);
  assert(shape->IsPart(node));
  assert(node->GetIndex() - shape->start_opt >= 0);
  assert(shape->partner != NULL);
  assert(shape->partner->start_param >= shape->start_param);

  return dynamic_cast<ShapeMapVariable*>(shape_param_[shape->partner->start_param + node->GetIndex() - shape->start_param]);
}

inline ShapeMapDesign::ShapeMapVariable* ShapeMapDesign::GetProfile(const ShapeMapVariable* node)
{
  const ShapeParam* shape = GetShape(node);
  return dynamic_cast<ShapeMapVariable*>(shape_param_[shape->partner->start_param + node->GetIndex() - shape->start_param]);
}


ShapeMapDesign::ShapeParam* ShapeMapDesign::FindShape(const ShapeMapVariable* spe)
{
  for(unsigned int s = 0; s < shape_.GetSize(); s++) {
    assert(shape_[s].end_param > 0);
    if((int) spe->GetIndex() < shape_[s].end_param)
      return &shape_[s];
  }

  return NULL;
}

const ShapeMapDesign::ShapeParam* ShapeMapDesign::FindShape(const Function* f, bool opt) const
{
  assert(num_node_shape_params_ > 0);
  assert(!(f->GetDesignType() == BaseDesignElement::PROFILE && IsProfileFixed()));

  for(unsigned int i = 0; i < shape_.GetSize(); i++)
  {
    const ShapeParam& shape = shape_[i];

    if(f->GetDesignType() == BaseDesignElement::PROFILE && shape.type != PROFILE)
      continue;

    if(opt && shape.start_opt == -1)
      continue;

    return &shape;
  }
  assert(false);
  return NULL;
}




StdVector<std::pair<ShapeMapDesign::ShapeParam*, ShapeMapDesign::ShapeParam*> > ShapeMapDesign::FindCenters()
{
  StdVector<std::pair<ShapeMapDesign::ShapeParam*, ShapeMapDesign::ShapeParam*> > result;

  for(unsigned int i = 0; i < shape_.GetSize(); i++)
  {
    ShapeParam& cand = shape_[i];
    assert(!(cand.other_center != NULL && cand.other_center->other_center == NULL));
    assert(!(cand.other_center != NULL && cand.other_center->other_center != &cand));
    if(cand.other_center != NULL && cand.idx < cand.other_center->idx) // pairs start with the lower idx such that we we have no doubles
      result.Push_back(std::make_pair(&cand, cand.other_center));
  }

  return result;
}

void ShapeMapDesign::ShapeParam::ParseSymmetry(ShapeParam* target, const PtrParamNode& pn)
{
  target->x_sym = pn->Get("left_right_sym")->As<string>() != "none";
  target->y_sym = pn->Get("bottom_up_sym")->As<string>() != "none";
  target->z_sym = pn->Has("front_back_sym") ? pn->Get("front_back_sym")->As<string>() != "none" : false;
  target->diag  = pn->Get("diagonal_sym")->As<string>() != "none";
}

void ShapeMapDesign::ShapeParam::ParseBounds(ShapeParam* target, const PtrParamNode& pn)
{
  // this is only for 3D center nodes when doing square symmetry, the 2D node doesn't know it
  if(pn->Has("slave"))
  {
    target->slave = pn->Get("slave")->As<bool>();
    if(target->slave && (pn->Has("initial") || pn->Has("lower") || pn->Has("upper") || pn->Has("fixed")))
      throw Exception("a slave node for cubical symmetry must not have initial/lower/upper/fixed");
  }

  if(!target->slave && pn->Has("initial")) // in the slave case we dont't have the properties but we don't want to return early
  {
    if(pn->Has("fixed"))
      throw Exception("shapeParam cannot have 'initial' and 'fixed' concurrently.");

    target->value = pn->Get("initial")->As<string>();
    if(!pn->Has("lower") || !pn->Has("upper"))
      throw Exception("shapeParam which is not fixed needs 'lower' and 'upper'");

    target->lower = pn->Get("lower")->As<string>();
    target->upper = pn->Get("upper")->As<string>();
    target->fixed = false;
  }

  if(!target->slave && pn->Has("fixed"))
  {
    if(pn->Has("initial"))
      throw Exception("shapeParam cannot have 'initial' and 'fixed' concurrently.");

    target->value = pn->Get("fixed")->As<string>();
    target->fixed = true;
  }

  target->clamp = pn->Get("clamp")->As<double>();

  if(!target->slave && !pn->Has("initial") && !pn->Has("fixed"))
    throw Exception("shapeParam needs to have either 'initial' or 'fixed'");

  if(target->fixed && target->clamp > 0)
    throw Exception("don't use 'clamp' for shapeParam together with 'fixed'.");
}

void ShapeMapDesign::ShapeParam::CopyBaseCenterProperties(ShapeParam* base)
{
  assert(base != NULL);
  if(base->orientation != ShapeParamElement::NOT_SET)
    orientation = base->orientation; // prevent overriding for induced center nodes

  x_sym = base->x_sym;
  y_sym = base->y_sym;
  z_sym = base->z_sym;
  diag = base->diag;
}

void ShapeMapDesign::ShapeParam::FlipOrientation(int center_node)
{
  assert(type == NODE);

  if(Is2DShape())
  {
    assert(domain->GetDim());
    assert(center_node == -1);
    std::swap(x_sym, y_sym);

    // the Flip() hase some more checks than just std::swap()
    dof = ShapeMapDesign::Flip(dof);
    orientation = ShapeMapDesign::Flip(dof);
  }
  else
  {
    assert(center_node == 0 || center_node == 1);
    assert(IsCenterShape());
    assert(GetFirstCenterNode() == this || GetSecondCenterNode() == this);
    assert(GetFirstCenterNode()->orientation == GetSecondCenterNode()->orientation);
    assert(GetFirstCenterNode()->dof != GetSecondCenterNode()->dof);

    ShapeParamElement::Dof new_o = center_node == 0 ? GetFirstCenterNode()->dof : GetSecondCenterNode()->dof;
    assert(new_o != ShapeParamElement::NOT_SET);

    // we spaw with one shape dof and orientation and set the new orientation also for the other shape
    if(center_node == 0)
      GetFirstCenterNode()->dof = GetFirstCenterNode()->orientation;
    else
      GetSecondCenterNode()->dof = GetFirstCenterNode()->orientation;

    GetFirstCenterNode()->orientation = new_o;
    GetSecondCenterNode()->orientation = new_o;
  }
}


void ShapeMapDesign::ShapeParam::ParseAndInit(PtrParamNode pn, ShapeParam* base)
{
  type = ShapeMapDesign::type.Parse(pn->GetName()); // can be CENTER, NODE and PROFILE!

  assert(!(type == CENTER && base != NULL));

  // is this pn part of a center?
  bool part_of_center = base != NULL && base->type == CENTER;

  // not for CENTER and PROFILE
  if(type == NODE)
    dof = ShapeParamElement::dof.Parse(pn->Get("dof")->As<std::string>());

  // not for part of center nodes
  //if((type == CENTER || type == NODE) && !part_of_center)
  if((type == CENTER || type == NODE) && !part_of_center)
    ParseSymmetry(this, pn);

  // a center child node
  if((type == NODE && part_of_center) || type == PROFILE)
    CopyBaseCenterProperties(base);

  // a real node
  if(type == NODE && !part_of_center)
    orientation = ShapeMapDesign::Flip(dof);

  // reads initial, upper, lower, ... for real node, center node and profile
  if(type != CENTER)
    ParseBounds(this, pn);

  LOG_DBG(SMD) << "SP:PAI idx=" << idx << " type=" << type
      << " base=" << (base == NULL ? -1 : base->type)
      << " dof=" << dof << " orientation=" << orientation
      << " x_sym=" << x_sym << " y_sym=" << y_sym << " z_sym=" << z_sym
      << " diag=" << diag
      << " induce_mirror=" << ShallInduceMirrorSymmetry()
      << " induce_clone=" << ShallInduceCloneSymmetry()
      << " induce_diag=" << ShallInduceDiagonalSymmetry();

  LOG_DBG(SMD) << "SP:PAI idx=" << idx << " type=" << type << " base=" << (base != NULL ? base->idx : -1) << " x_sym=" << x_sym
      << " y_sym=" << y_sym << " z_sym=" << z_sym << " orientation=" << orientation;
}

void ShapeMapDesign::ShapeParam::CopyProperties(const ShapeParam* ref, bool copy_master_data)
{
  assert(ref->type == CENTER || ref->type == NODE || ref->type == PROFILE);
  assert(type == NODE); // this shall be the default from the "constructor"
  assert(!(copy_master_data && !slave)); // makes only sense when we are slave
  if(ref->type != CENTER)
    type = ref->type;
  if(!copy_master_data)
    dof = ref->dof;
  if(type == PROFILE)
    dof = ShapeParamElement::NOT_SET;
  orientation = ref->orientation;
  lower = ref->lower;
  upper = ref->upper;
  value = ref->value;
  clamp = ref->clamp;
  max   = ref->max;
  fixed = ref->fixed;
  if(!copy_master_data) // don't make the slave a master
    slave = ref->slave;
  x_sym = ref->x_sym;
  y_sym = ref->y_sym;
  z_sym = ref->z_sym;
  diag  = ref->diag;
}

bool ShapeMapDesign::ShapeParam::ShallInduceMirrorSymmetry() const
{
  if(Is2DShape())
  {
    if(x_sym && orientation == ShapeParamElement::Y)
      return true;

    if(y_sym && orientation == ShapeParamElement::X)
      return true;
  }
  if(IsCenterShape())
  {
    // standing shape (dof=x,z, orientation=y) and a left_right_sym/x_sym:
    //  for dof=x there is mirror symmetry (x -> 1-x)
    //  for dof=z there is a clone symmetry (z -> z)
    // standing shape (dof=x,z, orientation=y) and a bottom_up_sym/y_sym:  Only mapping!
    if(x_sym && dof == ShapeParamElement::X)
      return true;

    if(y_sym && dof == ShapeParamElement::Y)
      return true;

    if(z_sym && dof == ShapeParamElement::Z)
      return true;
  }
  assert(!IsSurfaceShape());

  return false;
}

bool ShapeMapDesign::ShapeParam::ShallInduceCloneSymmetry() const
{
  if(Is2DShape())
    return false;

  assert(!IsSurfaceShape());

  // see comments in ShallInduceMirrorSymmetry
  if(x_sym && dof != ShapeParamElement::X && orientation != ShapeParamElement::X)
    return true;

  if(y_sym && dof != ShapeParamElement::Y && orientation != ShapeParamElement::Y)
    return true;

  if(z_sym && dof != ShapeParamElement::Z && orientation != ShapeParamElement::Z)
    return true;

  return false;
}


bool ShapeMapDesign::ShapeParam::ShallInduceDiagonalSymmetry() const
{
  assert(!IsSurfaceShape());

  // for 2D always true, for 3D we do not differentite about the angle
  return diag;
}


bool ShapeMapDesign::ShapeParam::ShallMapHalfShape() const
{
  assert(!IsSurfaceShape());

  if(x_sym && orientation == ShapeParamElement::X)
    return true;

  if(y_sym && orientation == ShapeParamElement::Y)
    return true;

  if(z_sym && orientation == ShapeParamElement::Z)
    return true;

  return false;
}


bool ShapeMapDesign::ShapeParam::ShallBeMasterOfSlave() const
{
  return other_center != NULL ? other_center->slave : false;
}

inline bool ShapeMapDesign::ShapeParam::IsFirstCenterNode() const
{
  return other_center != NULL && idx < other_center->idx;
}

inline bool ShapeMapDesign::ShapeParam::Is2DShape() const
{
  assert(!IsSurfaceShape());
  // IsCenterShape() might not be worked for not fully read center nodes
  return domain->GetDim() == 2;
}

/** are we first or second 3d center node or its profile? */
inline bool ShapeMapDesign::ShapeParam::IsCenterShape() const
{
  if(type == NODE)
    return other_center != NULL; // 2d and 3d
  else
    return partner != NULL && partner->other_center != NULL;
}


inline ShapeMapDesign::ShapeParam* ShapeMapDesign::ShapeParam::GetSecondCenterNode()
{
  assert(!(other_center != NULL && other_center->other_center != this));
  assert(!(other_center != NULL && type == PROFILE)); // due to unsymetry only nodes have the link
  assert(!(other_center != NULL && idx == other_center->idx));
  assert(idx >= 0 && (other_center == NULL || other_center->idx >= 0));
  assert(!(type == PROFILE && partner != NULL));
  if(other_center != NULL)
    return idx < other_center->idx ? other_center : this;
  if(type == PROFILE && partner->other_center != NULL)
    return partner->idx < partner->other_center->idx ? partner->other_center : partner;
  return NULL;
}

inline bool ShapeMapDesign::ShapeParam::IsSecondCenterNode() const
{
  return other_center != NULL && idx > other_center->idx;
}

inline bool ShapeMapDesign::ShapeParam::IsInduced() const
{
  return induce.master != NULL;
}

std::string ShapeMapDesign::ShapeParam::ToString() const
{
  std::stringstream ss;
  ss << ShapeMapDesign::type.ToString(this->type) << " idx=" << idx << " dof=" << dof << " o=" << orientation
      << " oc=" << (other_center == NULL ? -1 : other_center->idx)
      << " p=" << (partner == NULL ? -1 : partner->idx)
      << " sp=" << start_param << " ep=" << end_param
      << " so=" << start_opt << " eo=" << end_opt << " sl=" << slave
      << " ind=" << (IsInduced() ? (induce.reciprocal ? "reciprocal" : "yes") : "no");
  return ss.str();
}

void ShapeMapDesign::ShapeParam::ToInfo(PtrParamNode in)
{
  in->Get("idx")->SetValue(idx);
  in->Get("ref")->SetValue(GetReferenceId());
  in->Get("type")->SetValue(ShapeMapDesign::type.ToString(type));
  if(type == NODE)
  {
    in->Get("dof")->SetValue(ShapeParamElement::dof.ToString(dof));
    in->Get("orientation")->SetValue(ShapeParamElement::dof.ToString(orientation));
  }
  if(slave)
  {
    in->Get("slave")->SetValue(true);
  }
  else
  {
    if(fixed)
      in->Get("fixed")->SetValue(value);
    else
    {
      in->Get("lower")->SetValue(lower);
      in->Get("upper")->SetValue(upper);
      in->Get("clamp")->SetValue(clamp >= 0 ? std::to_string(clamp) : "no");
    }
    in->Get("variables")->SetValue(end_param - start_param);
    in->Get("design")->SetValue(end_opt - start_opt);
  }
  in->Get("induced")->SetValue(IsInduced() ? (induce.reciprocal ? "reciprocal" : "yes") : "no"); // we have this shape only as a result of a shape inducing symmetry
}

ShapeMapDesign::ShapeMapVariable::ShapeMapVariable(Type type, unsigned int index)
 : ShapeParamElement(type, index)
{
  coord.Resize(domain->GetDim(), -1.0);
  idx.Resize(domain->GetDim(), -1);
}


ShapeMapDesign::TanhSum::TanhSum()
{
  // set xrange[0:1]; a = 0.5; b=0.8; w=0.1; beta=30
  // ta(x) = 1-1/(exp(beta*(x-a+w)) + 1)
  // tb(x) = 1-1/(exp(beta*(x-b+w)) + 1)
  // l(x) = 1-1/(exp(11*(x-.5)) + 1)
  // plot ta(x), l(x), l(ta(x))
  // plot ta(x)+tb(x), l(x), l(ta(x)+tb(x))

  assert(beta > 0);
  // see Filter::SetNonLinCorrection()

  // it shall match the overlap of two shapes (2)
  scale = 1/(tanh(1.0)-tanh(0.0));
  offset = -scale * tanh(0.0);

  LOG_DBG(SMD) << "TS:TS scale=" << scale << " offset=" << offset << " tanh(0)=" << tanh(0) << " map(0)=" << map(0)
                   << " tanh(1)=" << tanh(1) << " map(1)=" << map(1) << " tanh(2)=" << tanh(2) << " map(2)=" << map(2);

  assert(close(map(0.0), 0.0, 1e-10));
  // assert(map(0.0) >= 0.0); TODO: fails for test suite -> check why!
  assert(close(map(1.0), 1.0, 1e-10));
  assert(map(2.0) <= 1.01);
}

inline double ShapeMapDesign::TanhSum::tanh(double x)
{
  return 1.0 - 1/(exp(beta*(x-0.5))+1);
}

inline double ShapeMapDesign::TanhSum::map(double x)
{
  // set xrange[0:1]; a = 0.5; b=0.8; w=0.1; beta=30
  // ta(x,beta) = 1-1/(exp(beta*(x-a+w)) + 1)
  // da(x,beta) = -1* (exp(beta*(x-a+w)) + 1)**-2 *beta*exp(beta*(x-a+w))
  // l(x) = 1-1/(exp(11*(x-.5)) + 1)
  // plot ta(x,beta), l(x), l(ta(x,beta/2))
  // plot ta(x)+tb(x), l(x), l(ta(x)+tb(x))

  return offset + scale * tanh(x);
}

inline double ShapeMapDesign::TanhSum::d_map(double x, double dx)
{
  // set xrange[0:1]; a = 0.5; b=0.8; w=0.1; beta=30
  // ta(x,beta) = 1-1/(exp(beta*(x-a+w)) + 1)
  // da(x,beta) = -1* (exp(beta*(x-a+w)) + 1)**-2 *beta*exp(beta*(x-a+w))
  // l(x) = 1-1/(exp(11*(x-.5)) + 1)
  // plot ta(x,beta), l(x), l(ta(x,beta/2)), da(x,beta), 11*exp(11*(ta(x,beta/2)-0.5))*(exp(11*(ta(x,beta/2)-0.5))+1)**-2 * da(x,beta/2)
  //
  // maxima:
  // f(x) := 1-1/(exp(11*(g(x)-0.5))+1);
  // diff(f(x),x);
  double e = exp(beta*(x-0.5));
  return scale * e * beta * dx / ((e+1) * (e+1));
  // return scale * std::pow(exp(beta*(x-0.5))+1, -2) * exp(beta*(x-0.5)) * beta * dx;
}

ShapeMapDesign::ElementSymmetry::ElementSymmetry(ShapeMapVariable* base)
{
  this->base = base;
}

void ShapeMapDesign::ElementSymmetry::AddSymmetryReference(ShapeMapVariable* elem, ShapeParam* shape, bool reciprocal)
{
  Virtual virt;
  virt.elem = elem;
  virt.shape = shape;
  virt.reciprocal = reciprocal;

  assert(!(reciprocal && shape->type != NODE));
  assert(elem->GetType() == Convert(shape->type));

  hidden.Push_back(virt);
}

inline void ShapeMapDesign::ElementSymmetry::ApplyDesign()
{
  for(unsigned int i = 0; i < hidden.GetSize(); i++)
  {
    Virtual& vir = hidden[i];
    assert(!(vir.reciprocal && vir.elem->GetType() != DesignElement::NODE));
    double val = vir.reciprocal ? (vir.shape->max - base->GetPlainDesignValue()) : base->GetPlainDesignValue();
    vir.elem->SetDesign(val);
    LOG_DBG2(SMD) << "ES:AD base=" << base->GetIndex() << "/" << base->dof_ << " vir=" << vir.elem->GetIndex() << "/" << vir.elem->dof_
        << " r=" << vir.reciprocal << " t=" << vir.elem->GetType() << " o=" << base->GetPlainDesignValue() << " -> " << val;
  }
}


double ShapeMapDesign::ElementSymmetry::GetPlainSymGradient(const Function* f) const
{
  double res = 0.0;

  for(unsigned int i = 0; i < hidden.GetSize(); i++)
  {
    const Virtual& vir = hidden[i];
    assert(!(vir.reciprocal && vir.elem->GetType() != DesignElement::NODE));
    res += (vir.reciprocal ? -1.0 : 1.0) * vir.elem->GetPlainGradient(f);
  }
  return res;
}

void ShapeMapDesign::ElementSymmetry::Reset(DesignElement::ValueSpecifier vs)
{
  for(unsigned int i = 0; i < hidden.GetSize(); i++)
    hidden[i].elem->Reset(vs);
}

void ShapeMapDesign::ElementSymmetry::PostInit(int objectives, int constraints)
{
  for(unsigned int i = 0; i < hidden.GetSize(); i++)
    hidden[i].elem->PostInit(objectives, constraints);
}


std::string ShapeMapDesign::ElementSymmetry::ToString(bool grad) const
{
  std::stringstream ss;
  ss << " ES b=" << base->GetIndex();
  for(unsigned int i = 0; i < hidden.GetSize(); i++)
  {
    ss << " [" << i << " ve=" << hidden[i].elem->GetIndex() << " s=" << hidden[i].shape->idx << " r=" << hidden[i].reciprocal;
    if(grad)
      ss << " dJ=" << hidden[i].elem->costGradient[0];
    ss << "] ";
  }
  return ss.str();
}

} // end of namespace


