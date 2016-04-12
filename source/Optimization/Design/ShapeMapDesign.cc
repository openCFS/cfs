/*
 * ShapeMapDesign.cc
 *
 *  Created on: Mar 11, 2016
 *      Author: fwein
 */
#include "Optimization/Design/ShapeMapDesign.hh"
#include "Optimization/Excitation.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"

namespace CoupledField {

DECLARE_LOG(SMD)
DEFINE_LOG(SMD, "shapeMapDesign")


ShapeMapDesign::ShapeMapDesign(StdVector<RegionIdType>& regionIds, PtrParamNode pn, ErsatzMaterial::Method method)
: AuxDesign(regionIds, pn, method), order_(10) // order 2 is technical minimum but too poor for practial use
{
  this->beta_ = pn->Get("shapeMap/beta")->As<double>();
  this->dim_ = domain->GetGrid()->GetDim();
  this->exoprt_fe_design_ = false; // we use the original design but don't communicate it via ReadDesignFromExtern(), ...
  this->tailing_aux_design_ = true; // we want shape_param_ or better opt_shape_param_ to take the role of DesignSpace::data

  assert(order_ >= 2); // too poor for technical use. 10 is nice
  if(order_ <= 3)
    info_->Get(ParamNode::HEADER)->Get(ParamNode::WARNING)->SetValue("low integration order for shape map");

  // set shape_, shape_param_ and map_
  SetupShapeDesign(pn->Get("shapeMap"));
  // copy the non-fixed stuff opt_shape_param_
  opt_shape_param_.Reserve(IsProfileFixed() ? num_node_shape_params_ : shape_param_.GetSize());
  for(unsigned int s = 0; s < shape_.GetSize(); s++)
    for(unsigned int p = shape_[s].start_param, n = shape_[s].end_param; !shape_[s].fixed && p < n; p++)
      opt_shape_param_.Push_back(&shape_param_[p]); // when we have fixed nodes we shall handle the index!!

  LOG_DBG(SMD) << "SMP rho_desig=" << data.GetSize();
  LOG_DBG(SMD) << "regions: " << regionIds.ToString();
}

int ShapeMapDesign::ReadDesignFromExtern(const double* space_in)
{
  int old_design = design_id;

  // write aux design variables (slack and alpha if any) last
  assert(exoprt_fe_design_ == false); // we do shape map
  assert(DesignSpace::GetNumberOfVariables() > 0); // we need this variables but they are hidden!

  bool new_design = false;

  for(unsigned int i = 0, n = opt_shape_param_.GetSize(); i < n; i++)
  {
    double v = space_in[i] * scaling_;
    if(!new_design && v != opt_shape_param_[i]->GetDesign(BaseDesignElement::PLAIN))
      new_design = true;

    opt_shape_param_[i]->SetDesign(v);
    LOG_DBG(SMD) << "RDFE: i=" << i << "-> " << v;
  }

  // append aux design, might also change design_id
  AuxDesign::ReadDesignFromExtern(space_in); // note the asserts above!

  if(new_design && design_id <= old_design)
    design_id++; // if new design and not already changed by AuxDesign

  // the design are shape parameters, map them to rho
  if(mapped_design_ != design_id)
    MapShapeToDensity();

  return design_id;
}

bool ShapeMapDesign::CompareDesign(const double* space_in)
{
  for(unsigned int i=0; i < opt_shape_param_.GetSize(); i++)
  {
    double v = space_in[i] * scaling_;
    if(v != opt_shape_param_[i]->GetDesign(BaseDesignElement::PLAIN))
      return false;
  }

  // no change here, let aux_design decide
  return AuxDesign::CompareDesign(space_in);
}

int ShapeMapDesign::WriteDesignToExtern(double* space_out, bool scale) const
{
  double rscaling = scale ? 1.0 / scaling_ : 1.0;

  for(unsigned int i=0; i < opt_shape_param_.GetSize(); i++)
  {
    space_out[i] = opt_shape_param_[i]->GetDesign(BaseDesignElement::PLAIN) * rscaling;
    LOG_DBG(SMD) << "WDTE: out[" << i << "]=" << space_out[i];
  }

  AuxDesign::WriteDesignToExtern(space_out, scale);

  return design_id;
}

void ShapeMapDesign::WriteGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Function* f, bool scaling)
{
  LOG_DBG(SMD) << "WGTE: ad=" << aux_design_.GetSize() << " osp=" << opt_shape_param_.GetSize() << " owst=" << out.window.GetStart() << " owsz=" << out.window.GetSize();
  // we cannot cache easily for mapped_constr_gradient_ as we would need it for each function.
  // MapShapeGradient would be good to perform it for all functions concurrently, however this is not possible as it is not the case that first all
  // simp function gradients are called and then all exported. This would need rewriting some stuff in cfs!
  if(f->IsObjective() || !Function::IsLocal(f->GetType())) // don't map local functions
    MapShapeGradient(f); // see comment above for what is necessary to cache the stuff

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
    assert(GetEndVarIdx(f,true) - GetFirstVarIdx(f,true) + aux_design_.GetSize() == out.window.GetSize());
    for(unsigned int s = GetFirstVarIdx(f,true), n = GetEndVarIdx(f,true); s < n ; s++) // for node (from 0) and profile (later) or default for both
    {
      assert(out.InWindow(base + s));
      out[base + s] = opt_shape_param_[s]->GetPlainGradient(f) * scaling;
      LOG_DBG2(SMD) << "WGTE f=" << f->ToString() << " ws=" << out.window.GetStart() << " s=" << s << " -> " << out[base + s];
    }
    // add slack stuff. No need to cheat window size
    AuxDesign::WriteGradientToExtern(out, vs, access, f, scaling);
  }
  else
  {
    assert(f->GetDesignType() == BaseDesignElement::NODE || f->GetDesignType() == BaseDesignElement::PROFILE);
    StdVector<unsigned int>& sparsity = f->GetSparsityPattern();
    assert(out.window.GetSize() == sparsity.GetSize());
    for(unsigned int i = 0; i < sparsity.GetSize(); i++)
    {
      unsigned int s = sparsity[i];
      assert(s < opt_shape_param_.GetSize());
      assert(out.InWindow(base + i));
      double scale = scaling ? scaling_ : 1.0;
      assert(vs == BaseDesignElement::CONSTRAINT_GRADIENT);
      out[base + i] = opt_shape_param_[s]->GetPlainGradient(f) * scale;
    }
  }
}

void ShapeMapDesign::Reset(DesignElement::ValueSpecifier vs, DesignElement::Type design)
{
  assert(design == BaseDesignElement::DEFAULT || design == BaseDesignElement::NODE); // extend to check for profile
  for(unsigned int i=0; i < opt_shape_param_.GetSize(); i++)
    opt_shape_param_[i]->Reset(vs);

  AuxDesign::Reset(vs, design);
}

void ShapeMapDesign::WriteBoundsToExtern(double* x_l, double* x_u) const
{
  for(unsigned int i=0; i < opt_shape_param_.GetSize(); i++)
  {
    x_l[i] = opt_shape_param_[i]->GetLowerBound() / scaling_;
    x_u[i] = opt_shape_param_[i]->GetUpperBound() / scaling_;
    LOG_DBG3(SMD) << "WBTE: l[" << i << "]=" << x_l[i] << " u[" << i << "]=" << x_u[i];
  }

  AuxDesign::WriteBoundsToExtern(x_l, x_u);
}

inline unsigned int ShapeMapDesign::GetNumberOfVariables() const
{
  return aux_design_.GetSize() + opt_shape_param_.GetSize();
}


inline BaseDesignElement* ShapeMapDesign::GetDesignElement(unsigned int idx)
{
  if(idx < opt_shape_param_.GetSize())
    return opt_shape_param_[idx];
  else
    return AuxDesign::GetDesignElement(idx); // handles its offset properly
}

int ShapeMapDesign::FindDesign(DesignElement::Type dt, bool throw_exception) const
{
  // check for DENSITY, ...
  int idx = DesignSpace::FindDesign(dt, false);
  if(idx >= 0)
    return idx;

  assert(dt == DesignElement::NODE || dt == DesignElement::PROFILE);

  for(unsigned int i = 0; i < shape_.GetSize(); i++)
    if(Convert(shape_[i].type) == dt)
      return i;

  if(throw_exception)
    EXCEPTION("Design " << DesignElement::type.ToString(dt) << " no FEM based and no shape mapping design.");
  return -1;
}

inline BaseDesignElement::Type ShapeMapDesign::Convert(Type type) const
{
  switch(type)
  {
  case NODE:
    return BaseDesignElement::NODE;
  case PROFILE:
    return BaseDesignElement::PROFILE;
  }
  assert(false);
  return BaseDesignElement::NO_TYPE;
}

inline bool ShapeMapDesign::IsProfileFixed() const
{
  // assume all nodes are concurrently fixed?!
  return shape_[num_node_shapes_].fixed;
}

unsigned int ShapeMapDesign::GetFirstVarIdx(const Function* f, bool opt) const
{
  assert(num_node_shape_params_ > 0);
  if(f->GetDesignType() == BaseDesignElement::DEFAULT || f->GetDesignType() == BaseDesignElement::NODE)
    return 0;
  assert(f->GetDesignType() == BaseDesignElement::PROFILE);
  assert(!IsProfileFixed()); // don't call if fixed
  return num_node_shape_params_; // assume no fixed node
}

/** small helper which gives the  index *after* the element based on type (node or profile) so*/
unsigned int ShapeMapDesign::GetEndVarIdx(const Function* f, bool opt) const
{
  assert(2 * num_node_shape_params_ == (int) shape_param_.GetSize()); // end of profile
  if(f->GetDesignType() == BaseDesignElement::NODE)
    return num_node_shape_params_; // assume no fixed node
  assert(f->GetDesignType() == BaseDesignElement::DEFAULT || f->GetDesignType() == BaseDesignElement::PROFILE);
  return opt ? opt_shape_param_.GetSize() : shape_param_.GetSize();
}

unsigned int ShapeMapDesign::GetFirstShapeIdx(const Function* f, bool opt) const
{
  assert(num_node_shapes_ == (int) shape_.GetSize() / 2);
  if(f->GetDesignType() == BaseDesignElement::DEFAULT || f->GetDesignType() == BaseDesignElement::NODE)
    return 0;
  assert(f->GetDesignType() == BaseDesignElement::PROFILE);
  return num_node_shapes_;
}

/** small helper which gives the  index *after* the element based on type (node or profile) so*/
unsigned int ShapeMapDesign::GetEndShapeIdx(const Function* f, bool opt) const
{
  assert(2 * num_node_shapes_ == (int) shape_.GetSize()); // end of profile
  if(f->GetDesignType() == BaseDesignElement::NODE)
    return num_node_shapes_;
  assert(f->GetDesignType() == BaseDesignElement::DEFAULT || f->GetDesignType() == BaseDesignElement::PROFILE);
  assert(shape_[shape_.GetSize() / 2].type == PROFILE);
  assert(!shape_[shape_.GetSize() / 2].fixed);
  return opt && IsProfileFixed() ? shape_.GetSize() / 2 : shape_.GetSize();
}

void ShapeMapDesign::SetupVirtualShapeElementMap(Function* f, StdVector<Function::Local::Identifier>& vem, Function::Local::Locality locality, Function::Local::Phase ph)
{
  assert(f != NULL);
  assert(f->IsLocal(f->GetType()));
  // note that we are called by the Function::Local() constructor, therefore Function::GetLocal() cannot work!
  assert(f->GetLocal() == NULL);

  // we assume fixed only for profile
  if(f->GetDesignType() == BaseDesignElement::PROFILE && IsProfileFixed())
    throw Exception("cannot have local constraint of shape map design 'profile' when this design is fixed.");

  // a lot copy&paste from Function::SetupVirtualElementMap()
  bool prev = locality == Function::Local::PREV_NEXT_AND_REVERSE || locality == Function::Local::PREV_NEXT;
  // next is always true!
  bool two_signs = locality == Function::Local::NEXT_AND_REVERSE || locality == Function::Local::PREV_NEXT_AND_REVERSE;

  int sign_1 = ph != Function::Local::BOTH ? (int) ph : two_signs ? 1 : Function::Local::Identifier::NO_SIGN;
  int sign_2 = ph != Function::Local::BOTH ? (int) ph : -1;

  // we don't set Function::Local::element_dimension_, it would be 2 (dim==1 * two signs)
  // we wont't use the full space as the individual shape_ are not connected
  vem.Reserve(num_node_shape_params_ * (two_signs ? 2 : 1)); // separately for node and profile but both of same size

  // traverse shape_ to have proper start and end. check for node or profile
  for(unsigned int s = GetFirstShapeIdx(f,true), n = GetEndShapeIdx(f,true); s < n; s++)
  {
    const ShapeParam& param = shape_[s];
    assert(f->GetDesignType() == Convert(param.type)); // NODE or PROFILE
    assert(!param.fixed);

    // skip the last element as we want only 'full' elements with next
    for(int e = param.start_param + (prev ? 1 : 0); e < param.end_param -1; e++)
    {
      BaseDesignElement& bde = shape_param_[e];
      assert(f->GetDesignType() == bde.GetType());

      BaseDesignElement* prev_de = prev ? &shape_param_[e-1] : NULL;
      BaseDesignElement* next_de = &shape_param_[e+1];

      vem.Push_back(Function::Local::Identifier(&bde, prev_de, next_de, sign_1));
      if(two_signs)
        vem.Push_back(Function::Local::Identifier(&bde, prev_de, next_de, sign_2));
    }
  }

  LOG_DBG(SMD) << "SVSEM f=" << f->ToString() << " loc=" << locality << " ts=" << two_signs << " prev=" << prev << " -> vem=" << vem.GetSize();

}

void ShapeMapDesign::ToInfo(PtrParamNode in, ErsatzMaterial* em)
{
  AuxDesign::ToInfo(in, em);
  PtrParamNode base_ = in->Get("designVariables");

  for(unsigned int i = 0; i < shape_.GetSize(); i++)
    shape_[i].ToInfo(base_->Get("shapeParam", ParamNode::APPEND));
}


StdVector<unsigned int> ShapeMapDesign::SetupLexicographicMesh(Grid* grid, const RegionIdType design_reg, StdVector<int>& elem_to_idx, StdVector<int>& idx_to_elem)
{
  // for 2D n_z_=1
  StdVector<unsigned int> n = grid->GetBoundaries(design_reg);
  unsigned int n_elems = n[0] * n[1] * n[2];

  // TODO here we assume that the whole mesh is for LBM and the mesh is of lexicographic ordering.
  // To be good this needs to handled by element neighbors!
  if(grid->GetNumElems() != n_elems)
    EXCEPTION("the current implementation assumes the whole mesh to be used for LBM and lexicographic ordering. Mesh has " << grid->GetNumElems() << " but we assume " << n_elems);

  idx_to_elem.Resize(n_elems);
  for(unsigned int i = 0; i < n_elems; i++)
    idx_to_elem[i] = i+1; // one based

  elem_to_idx.Resize(n_elems + 1); // one-based elem_nr
  for(unsigned int i = 0, n = elem_to_idx.GetSize(); i < n; i++)
    elem_to_idx[i] = i-1; // -1 for not appropriate idx

  return n;
}

/* <shapeMap beta="2">
    <node dof="x" lower="0" upper=".5" initial=".25"/>
    <profile lower=".01" upper=".1" initial=".1"/>
  </shapeMap> */
 void ShapeMapDesign::SetupShapeDesign(PtrParamNode pn)
{
   // check rhos which should be already be set in DesignSpace::data
   assert(GetRegionIds().GetSize() == 1); // more is not implemented yet
   StdVector<int> elem_to_idx;
   StdVector<int> idx_to_elem;
   StdVector<unsigned int> n = SetupLexicographicMesh(domain->GetGrid(), GetRegionIds().First(), elem_to_idx, idx_to_elem);
   nx_ = n[0];
   ny_ = n[1];
   nz_ = n[2];
   assert(data.GetSize() == nx_ * ny_ * nz_); // DesignSpace::data has an element for each FEM-Cell
   assert(n[0] == n[1]); // up to now
   assert(n[2] == 1); // 2D

   // read shapeParam to shape_
   assert(beta_ >= 0.0); // shall be already set
   ParamNodeList nodes = pn->GetList("node"); // there must be at least one
   PtrParamNode  profile = pn->Get("profile"); // there must be one
   LOG_DBG(SMD) << "SSD: childs of shapeParam: " << nodes.GetSize();
   assert(!nodes.IsEmpty());
   shape_.Resize(nodes.GetSize() * 2); // list is for nodes which are doubled by profile
   for(unsigned int i = 0, n = nodes.GetSize(); i < n; i++) {
     shape_[i].Init(nodes[i], i); // empty constructor due to StdVector/(
     shape_[n+i].Init(profile, n+i);
   }
   num_node_shapes_ = nodes.GetSize();

   // setup nodes within shape_param_
   StdVector<ShapeParam*> shape_x = FindShape(NODE, 0);
   StdVector<ShapeParam*> shape_y = FindShape(NODE, 1);

   num_node_shape_params_ = (nx_+1) * shape_x.GetSize() + (ny_+1) * shape_y.GetSize(); // one node more than elements
   shape_param_.Reserve(2 * num_node_shape_params_); // all doubled for profile
   assert(shape_param_.Capacity() > 0);
   for(unsigned int s = 0; s < shape_x.GetSize(); s++)
   {
     shape_x[s]->start_param = shape_param_.GetSize();
     for(unsigned int y = 0; y < ny_+1; y++)
       CreateShapeVariable(shape_x[s], y);
     shape_x[s]->end_param = shape_param_.GetSize();
   }
   for(unsigned int s = 0; s < shape_y.GetSize(); s++)
   {
     shape_y[s]->start_param = shape_param_.GetSize();
     for(unsigned int x = 0; x < nx_+1; x++)
       CreateShapeVariable(shape_y[s], x);
     shape_y[s]->end_param = shape_param_.GetSize();
   }
   // add the profiles
   for(int n = 0; n < num_node_shapes_; n++)
   {
     ShapeParam& node = shape_[n];
     ShapeParam& prof = shape_[num_node_shapes_ + n];
     assert(node.type == NODE && prof.type == PROFILE);
     prof.start_param = shape_param_.GetSize();
     for(int e = 0; e < node.end_param - node.start_param; e++)
       CreateShapeVariable(&prof, e);
     prof.end_param = shape_param_.GetSize();
     assert(prof.end_param - prof.start_param == node.end_param - node.start_param);
     assert(prof.start_param - num_node_shape_params_ == node.start_param);
   }

   // set map_ to map from shape_param to DesignSpace::data, fill with shape_param_
   map_.Resize(data.GetSize());
   for(unsigned int i = 0, n = map_.GetSize(); i < n; i++)
   {
     map_[i].rho = &(data[Find(idx_to_elem[i])]); // is very fast and gives a layer for arbitrary element ordering in the mesh
     map_[i].nodes.Reserve(2 * num_node_shapes_); // each design node connects to two density elements
     map_[i].ip_param_idx.Resize(order_ * order_);
     map_[i].ip_param_idx.Init(-1);
   }

   for(int i = 0; i < num_node_shape_params_; i++)
   {
     ShapeParamElement& spe = shape_param_[i];
     if(spe.dof == 0) // x
     {
       int y = spe.idx[1];
       assert(y >= 0);

       // exclude the case where the y node is at the upper boundary (one node more than elements!)
       for(unsigned int x = 0; x < nx_; x++){
         if(y < (int) ny_) // are we the topmost node ontop of the last element
           map_[DensityIdx(x, y)].nodes.Push_back(&spe);
         if(y-1 >= (int) 0) // node 2 participates to the lower element (2-1) and the upper element (2)
           map_[DensityIdx(x, y-1)].nodes.Push_back(&spe);
       }
     }

     if(spe.dof == 1) // y
     {
       int x = spe.idx[0];
       for(unsigned y = 0; y < ny_; y++)  {
         if(x < (int) nx_)
           map_[DensityIdx(x, y)].nodes.Push_back(&spe);
         if(x-1 >= 0)
           map_[DensityIdx(x-1, y)].nodes.Push_back(&spe);
       }
     }
   }
   //DumpMap();
   MapShapeToDensity();
}

 void ShapeMapDesign::MapShapeToDensity()
 {
   LOG_DBG(SMD) << "MSTD: di=" << design_id;
   Grid* grid = domain->GetGrid();
   // within the element coordinates we perform the integration
   Matrix<double> coords;
   // this are shapes we need to integrated. If too far away we dont't need them if all then the 0....n/2 with n size of item.para.GetSize()
   StdVector<unsigned int> shapes;

   assert(data.GetSize() == map_.GetSize());
   for(unsigned int r = 0, n = map_.GetSize(); r < n; r++)
   {
     Item& item = map_[r];
     DesignElement* de = item.rho;
     grid->GetElemNodesCoord(coords, de->elem->connect, false); // no deformed mesh
     // the index we store for the gradient. We gather the information when computing rho and save computing all the tanh again
     item.ip_param_idx.Init(-1);

     // reduce integration by identifying the structures close enough to zero. One could save more tanh if also the close to one are skipped.
     shapes.Clear(true); // keep capacity
     assert(item.nodes.GetSize() % 2 == 0); // we expect to have pairs
     for(unsigned int s = 0; s < item.nodes.GetSize(); s+=2)
       if(ApproxMaxRho(item.nodes[s], item.nodes[s+1], coords) > 1e-10)
         shapes.Push_back(s);

     LOG_DBG3(SMD) << "MS2D: -> el=" << de->elem->elemNum << " shape=" << shapes.ToString();

     double rho = 0.0; // sum up over all ips (if shapes is large enough :))

     // it makes sense to traverse first the ip and then the variables
     for(unsigned int ip_x = 0; ip_x < order_; ip_x++)
     {
       for(unsigned int ip_y = 0; ip_y < order_; ip_y++)
       {
         double max_ip_rho = 0.0; // we need 0.0 for a valid final rho if t is too small everywhere but we need to overwrite the default idx -1!
         for(unsigned int si = 0; si < shapes.GetSize(); si++)
         {
           double t = Eval(item.nodes[shapes[si]], item.nodes[shapes[si]+1], coords, ip_x, ip_y, false, false); // no derivative
           if(t >= max_ip_rho) { // equal is important to overwrite index -1 with a 0.0 rho and not to keep it
             max_ip_rho = t;
             item.ip_param_idx[ip_y*order_+ip_x] = shapes[si]; // x fastest, easy to extend to 3D
           }
         }
         rho += max_ip_rho;
       } // end ip_y
     } // end ip_x

     de->SetDesign(de->GetLowerBound() + (de->GetUpperBound() - de->GetLowerBound()) * (rho / (order_ * order_))); // we assume 0 <= v <= 1
     assert(de->GetPlainDesignValue() <= de->GetUpperBound() + 1e-10);
     assert(de->GetPlainDesignValue() >= de->GetLowerBound() - 1e-10);
     LOG_DBG3(SMD) << "MS2D: -> el=" << de->elem->elemNum << " -> avg=" << de->GetPlainDesignValue(); // << " pi=" << item.ip_param_idx.ToString();
     // Matrix<double> m;
     // m.Assign(v, order_, order_, true);
     // LOG_DBG3(SMD) << m.ToString();
   } // end loop over density elements

   mapped_design_ = design_id;
 }

 void ShapeMapDesign::MapShapeGradient(const Function* f)
 {
   assert(design_id == mapped_design_); // we need the Item setting from MapShapeDesign for the current design!
   assert(!(!f->IsObjective() && dynamic_cast<const Condition*>(f)->IsLocalCondition())); // it makes no sense for a local condition!!

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
   int res_idx_da = GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::SHAPE_MAP_GRAD, DesignElement::SM_NODE);
   int res_idx_dw = GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::SHAPE_MAP_GRAD, DesignElement::SM_PROFILE);

   Grid* grid = domain->GetGrid();
   // within the element coordinates we perform the integration
   Matrix<double> coords;
   for(unsigned int r = 0, n = map_.GetSize(); r < n; r++) // traverse all rho design elements
   {
     Item& item = map_[r];
     DesignElement* de = item.rho;
     grid->GetElemNodesCoord(coords, de->elem->connect, false); // no deformed mesh
     double log_da = 0.0;
     double log_dw = 0.0;

     // either all ip_param_idx are -1 or none
     if(item.ip_param_idx[0] > -1)
     {
       for(unsigned int ip_x = 0; ip_x < order_; ip_x++)
       {
         for(unsigned int ip_y = 0; ip_y < order_; ip_y++)
         {
           unsigned int ip_idx    = ip_y*order_+ip_x;
           int shape_idx = item.ip_param_idx[ip_idx];
           assert(shape_idx >= 0); // all -1 or none
           assert(shape_idx % 2 == 0); // shall be even
           // select the dominant shape parameter which lead to max rho for this ip
           ShapeParamElement* s1 = item.nodes[shape_idx];
           ShapeParamElement* s2 = item.nodes[shape_idx+1];

           double da = Eval(s1, s2, coords, ip_x, ip_y, true, false); // dtanh_da
           // normalize FIXME! For a reason I don't understand the 0.5 makes the snopt gradient check work?!
           double da_norm = (de->GetUpperBound() - de->GetLowerBound()) * 0.5 * da / (order_ * order_);
           log_da += da_norm;

           assert(opt_ != NULL);
           assert(opt_->objectives.data.GetSize() == s1->costGradient.GetSize());
           s1->AddGradient(f, de->GetPlainGradient(f) * da_norm);
           s2->AddGradient(f, de->GetPlainGradient(f) * da_norm);

           if(!IsProfileFixed())
           {
             ShapeParamElement* w1 = GetProfile(s1);
             ShapeParamElement* w2 = GetProfile(s2);
             double dw = Eval(s1, s2, coords, ip_x, ip_y, false, true); // dtanh_dw
             // the first 0.5 is not understood as above and the second 0.5 is because we apply 0.5*profile to tanh
             double dw_norm = (de->GetUpperBound() - de->GetLowerBound()) * 0.5 * 0.5 * dw / (order_ * order_);
             log_dw += dw_norm;

             w1->AddGradient(f, de->GetPlainGradient(f) * dw_norm);
             w2->AddGradient(f, de->GetPlainGradient(f) * dw_norm);
           }
           //LOG_DBG3(SMD) << "MSG: -> el=" << de->elem->elemNum << " rho=" << de->GetPlainDesignValue() << " ip_x=" << ip_x
           //    << " ip_y=" << ip_y << " ip_idx=" << ip_idx << " si=" << shape_idx << " s1=" << s1->GetIndex()
           //    << " s2=" << s2->GetIndex() << " da=" << da << " da_norm=" << da_norm << " dw=" << dw << " dw_norm=" << dw_norm;

         } // end ip_y
       } // end ip_x
       // normalize by integration points.
     }
     else
     {
       LOG_DBG2(SMD) << "MSG: f=" << f->ToString() << " el=" << de->elem->elemNum << " item.ip_param_idx[0]=" << item.ip_param_idx[0] << " rho=" << de->GetPlainDesignValue() << " sum da=" << log_da;
       // when one ip_param_idx is -1, all are -1.
       // the reason is, that all are -1 if ApproxMaxRho() to small for all shapes. Otherwise at least one shape has the ips
       assert(item.ip_param_idx.Sum() == -1 * (int) (order_*order_));
     }
     LOG_DBG2(SMD) << "MSG: el=" << de->elem->elemNum << " rho=" << de->GetPlainDesignValue() << " sum da=" << log_da << " sum dw=" << log_dw;
     if(res_idx_da >= 0)
       de->specialResult[res_idx_da] = log_da;
     if(res_idx_dw >= 0)
       de->specialResult[res_idx_dw] = log_dw;
   } // end loop over density elements
 }

 inline double ShapeMapDesign::Eval(const ShapeParamElement* s1, const ShapeParamElement* s2, const Matrix<double>& coords, unsigned int ip_x, unsigned int ip_y, bool grad_a, bool grad_w) const
 {
   assert(dim_ == 2);
   assert(s1->dof == s2->dof);
   assert(s1->GetType() == BaseDesignElement::NODE && s2->GetType() == BaseDesignElement::NODE);
   assert(!(grad_a == true && grad_w == true)); // the other three combinations are allowed
   // we assume the elements to be oriented ccw starting lower left
   assert(coords.GetNumRows() == 2); // dim == 2
   assert(coords.GetNumCols() == 4); // dim == 4 nodes
   assert(close(coords[0][0], coords[0][3])); // x-pos: lower left and upper left
   assert(close(coords[0][1], coords[0][2])); // x-pos: lower right and upper right
   assert(coords[0][0] < coords[0][1]);       // x-pos: lower left is smaller than lower right
   assert(close(coords[1][0], coords[1][1])); // y-pos: lower left and lower right
   assert(close(coords[1][2], coords[1][3])); // y-pos: upper right and upper left
   assert(coords[1][1] < coords[1][2]);       // y-pos: lower right smaller upper right

   int dof = s1->dof;
   // element position
   double start = dof == 0 ? coords[0][0] : coords[1][0];
   double end   = dof == 0 ? coords[0][1] : coords[1][3];
   // the parameters
   double a1 = s1->GetDesign(BaseDesignElement::PLAIN);
   double a2 = s2->GetDesign(BaseDesignElement::PLAIN);
   // the profiles
   assert(GetProfile(s1)->GetType() == BaseDesignElement::PROFILE);
   double w1 = 0.5 * GetProfile(s1)->GetDesign(BaseDesignElement::PLAIN); // half profile as tanh wants the half width
   double w2 = 0.5 * GetProfile(s2)->GetDesign(BaseDesignElement::PLAIN);

   double xy = start + (dof == 0 ? ip_x : ip_y) /(order_-1.) * (end-start); // for dof=0 (x) xy is x and might be far away from a
   double a  = a1 + (dof == 0 ? ip_y : ip_x)/(order_-1.) * (a2-a1);         // for dof=1 (c) a is y and with a1=a2 we have the same value for a
   double w  = w1 + (dof == 0 ? ip_y : ip_x)/(order_-1.) * (w2-w1);
   double val = -1.0;
   if(grad_a)
     val = d_tanh_da(xy, a, w);
   if(grad_w)
     val = d_tanh_dw(xy, a, w);
   if(!grad_a && !grad_w)
     val = tanh(xy, a, w);
   LOG_DBG3(SMD) << "E: s1=" << s1->GetIndex() << " s2=" << s2->GetIndex() << " dof=" << dof << " start=" << start << " end=" << end  << " a=" << a1 << "..." << a2
                 << " ip_x=" << ip_x << " ip_y=" << ip_y << " xy=" << xy << " a=" << a << " da:" << grad_a << " dw=" << grad_w << " -> " << val;
   return val;
 }

 inline double ShapeMapDesign::ApproxMaxRho(const ShapeParamElement* s1, const ShapeParamElement* s2, const Matrix<double>& coords) const
 {
   assert(dim_ == 2);
   double max        = Eval(s1, s2, coords, 0,               0, false, false); // false=no derivative
   max = std::max(max, Eval(s1, s2, coords, 0,        order_-1, false, false));
   max = std::max(max, Eval(s1, s2, coords, order_-1,        0, false, false));
   max = std::max(max, Eval(s1, s2, coords, order_-1, order_-1, false, false));
   return max;
 }


 inline double ShapeMapDesign::tanh(double x, double a, double w) const
 {
   // plot 1-1/(exp(2*beta*(x-a+w)) + 1), 1/(exp(2*beta*(x-a-w)) + 1)
   if(x <= a)
     return 1.0 - 1/(std::exp(2*beta_*(x-a+w))+1);
  else
    return 1/(std::exp(2*beta_*(x-a-w))+1);
 }

 inline double ShapeMapDesign::d_tanh_da(double x, double a, double w) const
 {
   // plot -1* (exp(2*beta*(x-a+w)) + 1)**-2 *2*beta*exp(2*beta*(x-a+w)), (exp(2*beta*(x-a-w))+1)**-2 *2*beta*exp(2*beta*(x-a-w))
   //
   // matlab:
   //if x <= a
   //  f = -1* (exp(2*beta*(x-a+d)) + 1)^-2 *2*beta*exp(2*beta*(x-a+d));
   //else
   ///  f = (exp(2*beta*(x-a-d))+1)^-2 *2*beta*exp(2*beta*(x-a-d));
   // end

   if(x <= a)
     return -1./((std::exp(2*beta_*(x-a+w))+1) * (std::exp(2*beta_*(x-a+w))+1)) * 2*beta_*std::exp(2*beta_*(x-a+w));
  else
    return 1/((std::exp(2*beta_*(x-a-w))+1) * (std::exp(2*beta_*(x-a-w))+1)) * 2*beta_*std::exp(2*beta_*(x-a-w));
 }

 inline double ShapeMapDesign::d_tanh_dw(double x, double a, double w) const
 {
   if(x <= a)
     return 1./((std::exp(2*beta_*(x-a+w))+1) * (std::exp(2*beta_*(x-a+w))+1)) * 2*beta_*std::exp(2*beta_*(x-a+w));
  else
    return 1/((std::exp(2*beta_*(x-a-w))+1) * (std::exp(2*beta_*(x-a-w))+1)) * 2*beta_*std::exp(2*beta_*(x-a-w));
 }


 void ShapeMapDesign::DumpMap()
 {
    for(unsigned int y = 0; y < ny_; y++)
    {
      for(unsigned int x = 0; x < nx_; x++)
      {
        Item i = map_[DensityIdx(x, y)];
        std::cout << "y=" << y << " x=" << x << " elidx=" << DensityIdx(x, y);
        for(unsigned int s = 0; s < i.nodes.GetSize(); s++)
           std::cout << " [nr=" << i.nodes[s]->GetIndex() << " dof=" << i.nodes[s]->dof << " free=" << i.nodes[s]->idx[1-i.nodes[s]->dof] << "]"; // works only for 2D
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }
 }

void ShapeMapDesign::CreateShapeVariable(const ShapeParam* param, int free)
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


  // add element to shape_param_
  shape_param_.Push_back(ShapeParamElement(Convert(param->type), shape_param_.GetSize()));
  ShapeParamElement& spe = shape_param_.Last();
  spe.SetLowerBound(param->lower);
  spe.SetUpperBound(param->upper);
  spe.SetDesign(param->value);
  spe.dof = param->dof;
  if(spe.dof == 0) {
    spe.coord[1] = (double) free / ny_; // free is max ny_
    spe.idx[1] = free;
  }
  if(spe.dof == 1) {
    spe.coord[0] = (double) free / nx_;
    spe.idx[0] = free;
  }

  // PostInit() sets arrays for objective and constraint gradients

  LOG_DBG2(SMD) << "CSV el=" << (shape_param_.GetSize() - 1) << " dof=" << spe.dof << " free=" << free << " d=" << spe.GetDesign(BaseDesignElement::PLAIN) << " coord=" << spe.coord.ToString();
}

void ShapeMapDesign::PostInit(int objectives, int constraints)
{
  // full_data is only used for internal optimizers to have a consecutive design array
  full_data.Clear(false); // release memory
  full_data.Reserve(GetNumberOfVariables()); // shape param + aux
  full_data.Resize(opt_shape_param_.GetSize()); // was set in DesignSpace and will be released
  for(unsigned int i = 0; i < opt_shape_param_.GetSize(); i++)
    full_data[i] = opt_shape_param_[i];

  // add aux_design to full_data
  AuxDesign::PostInit(objectives, constraints);

  for(unsigned int i = 0, n = opt_shape_param_.GetSize(); i < n; i++)
    opt_shape_param_[i]->PostInit(objectives, constraints);

  if(domain->GetOptimization() != NULL)
    opt_ = domain->GetOptimization();
}

StdVector<ShapeMapDesign::ShapeParam*> ShapeMapDesign::FindShape(Type type, int dof)
{
  StdVector<ShapeParam*> res;
  for(unsigned int i = 0; i < shape_.GetSize(); i++) // could be smarter bur more complex if we consider the type
    if(shape_[i].type == type && shape_[i].dof == dof)
      res.Push_back(&shape_[i]);
  return res;
 }

 void ShapeMapDesign::ShapeParam::Init(PtrParamNode pn, unsigned int idx_)
 {
   type = ShapeMapDesign::type.Parse(pn->GetName());
   idx = (int) idx_;
   if(type == NODE)
   {
     if(!pn->Has("dof"))
       throw Exception("shapeParam of type 'node' requires 'dof'");
     std::string std_dof = pn->Get("dof")->As<std::string>();
     if(std_dof == "x")
       dof = 0;
     if(std_dof == "y")
       dof = 1;
     if(dof == -1) // default, no z yet
       throw Exception("shapeParam of type 'node' has invalid dof '" + std_dof + "'");
   }
   else
     if(pn->Has("dof"))
       throw Exception("shapeParam knows 'dof' only for type 'node'");
   if(pn->Has("initial"))
   {
     if(pn->Has("fixed"))
       throw Exception("shapeParam cannot have 'initial' and 'fixed' concurrently.");
     value = pn->Get("initial")->As<double>();
     if(!pn->Has("lower") || !pn->Has("upper"))
       throw Exception("shapeParam which is not fixed needs 'lower' and 'upper'");
     lower = pn->Get("lower")->As<double>();
     upper = pn->Get("upper")->As<double>();
     fixed = false;
   }
   if(pn->Has("fixed"))
   {
     if(pn->Has("initial"))
       throw Exception("shapeParam cannot have 'initial' and 'fixed' concurrently.");
     value = pn->Get("fixed")->As<double>();
     fixed = true;
   }
   if(!pn->Has("initial") && !pn->Has("fixed"))
     throw Exception("shapeParam needs to have either 'initial' or 'fixed'");
}

void ShapeMapDesign::ShapeParam::ToInfo(PtrParamNode in)
{
  in->Get("idx")->SetValue(idx);
  in->Get("type")->SetValue(ShapeMapDesign::type.ToString(type));
  in->Get("dof")->SetValue(dof == 0 ? "x" : "y");
  in->Get("lower")->SetValue(lower);
  in->Get("upper")->SetValue(upper);
  assert(start_param >= 0 && end_param > 0);
  in->Get("variables")->SetValue(end_param - start_param);
}


} // end of namespace



