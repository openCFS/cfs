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
  this->alsomatopt_ = false; // we use the original design but don't communicate it via ReadDesignFromExtern(), ...

  assert(order_ >= 2); // too poor for technical use. 10 is nice
  if(order_ <= 3)
    info_->Get(ParamNode::HEADER)->Get(ParamNode::WARNING)->SetValue("low integration order for shape map");

  SetupShapeDesign(pn->Get("shapeMap"));

  LOG_DBG(SMD) << "SMP rho_desig=" << data.GetSize();
  LOG_DBG(SMD) << "regions: " << regionIds.ToString();
}

int ShapeMapDesign::ReadDesignFromExtern(const double* space_in)
{
  int old_design = design_id;

  // write aux design variables (slack and alpha if any) first
  assert(alsomatopt_ == false); // we do shape map
  assert(DesignSpace::GetNumberOfVariables() > 0); // we need this variables but they are hidden!
  AuxDesign::ReadDesignFromExtern(space_in); // note the asserts above!

  // design_id might be changed above in AuxDesign::ReadDesignFromExtern()
  bool new_design = old_design != design_id;

  unsigned int offset = aux_design_.GetSize();

  for(unsigned int i = 0, n = shape_param_.GetSize(); i < n; i++)
  {
    double v = space_in[offset + i] * scaling_;
    if(!new_design && v != shape_param_[i].GetDesign())
      new_design = true;

    shape_param_[i].SetDesign(v);
    LOG_DBG(SMD) << "RDFE: i=" << i << "-> " << v;
  }
  if(new_design && design_id <= old_design)
    design_id++; // if new design and not already changed by AuxDesign

  // the design are shape parameters, map them not to rho
  MapShapeToDensity();

  return design_id;
}

bool ShapeMapDesign::CompareDesign(const double* space_in)
{
  if(!AuxDesign::CompareDesign(space_in))
    return false;

  unsigned int offset = aux_design_.GetSize();

  for(unsigned int i=0; i < shape_param_.GetSize(); i++)
  {
    double v = space_in[offset + i] * scaling_;
    if(v != shape_param_[i].GetDesign())
      return false;
  }

  return true;
}

int ShapeMapDesign::WriteDesignToExtern(double* space_out, bool scale) const
{
  AuxDesign::WriteDesignToExtern(space_out, scale);

  double rscaling = scale ? 1.0 / scaling_ : 1.0;
  unsigned int offset = aux_design_.GetSize();

  for(unsigned int i=0; i < shape_param_.GetSize(); i++)
  {
    space_out[offset + i] = shape_param_[i].GetDesign() * rscaling;
    LOG_DBG(SMD) << "WDTE: out[" << i << "]=" << space_out[offset +i];
  }
  return design_id;
}

void ShapeMapDesign::WriteGradientToExtern(StdVector<double>& out, DesignElement::ValueSpecifier vs, DesignElement::Access access, Function* f, bool scaling)
{
  LOG_DBG(SMD) << "WGTE: ad=" << aux_design_.GetSize() << " sp=" << shape_param_.GetSize() << " owst=" << out.window.GetStart() << " owsz=" << out.window.GetSize();

  assert(design_id == mapped_design_); // our mapping is the current one and we can use the information for the gradient
  if(f->IsObjective() && design_id != mapped_obj_gradient_)
    MapShapeGradient(true); // needs to run for the first function and processes for all functions
  if(!f->IsObjective() && design_id != mapped_constr_gradient_)
    MapShapeGradient(false); // shall be done smarter !

  assert((f->IsObjective() && design_id == mapped_obj_gradient_) || (!f->IsObjective() && design_id == mapped_constr_gradient_));

  assert(f != NULL);
  assert(opt_->objectives.data.GetSize() == 1); // implement multi objective and be careful!
  assert(!(vs == DesignElement::COST_GRADIENT && !f->IsObjective()));
  assert(vs == DesignElement::COST_GRADIENT || vs == DesignElement::CONSTRAINT_GRADIENT);
  assert(!opt_->GetMultipleExcitation()->DoMetaExcitation(f->ctxt)); // robustness and transformation don't make sense for shape map. In DesignSpace we do f->GetExcitation()->Apply()
  assert(design.GetSize() == 1); // only pseudo density, nothing else implemented yet
  assert(regions.GetSize() == 1); // needs to be as design shall be 1

  assert(out.window.Initialized());
  // out contains the Jacobian for a possibly many functions. Snopt combies cost function (first) and then the constraints derivatives
  // Here we assume that out has a window set where to write to for the given function.

  // cheat window to write the aux variables (slack and alpha)
  StdVector<double>::Window org_window = out.window;
  out.window.Set(out.window.GetStart(), aux_design_.GetSize());
  AuxDesign::WriteGradientToExtern(out, vs, access, f, scaling);
  out.window = org_window;

  if(f->HasDenseJacobian())
  {
    unsigned int n0 = out.window.GetStart() + aux_design_.GetSize(); // to grow up to the total number of design variables

    for(unsigned int s = 0, n = shape_param_.GetSize(); s < n ; s++)
    {
      assert(out.InWindow(n0 + s));
      out[n0 + s] = shape_param_[s].GetPlainGradient(f) * scaling;
    }
  }
  else
  {
   assert(false); // not yet implemented
  }
}

void ShapeMapDesign::Reset(DesignElement::ValueSpecifier vs, DesignElement::Type design)
{
  AuxDesign::Reset(vs, design);

  for(unsigned int i=0; i < shape_param_.GetSize(); i++)
    shape_param_[i].Reset(vs);
}

void ShapeMapDesign::WriteBoundsToExtern(double* x_l, double* x_u) const
{
  AuxDesign::WriteBoundsToExtern(x_l, x_u);

  unsigned int offset = aux_design_.GetSize();

  for(unsigned int i=0; i < shape_param_.GetSize(); i++)
  {
    x_l[offset + i] = shape_param_[i].GetLowerBound() / scaling_;
    x_u[offset + i] = shape_param_[i].GetUpperBound() / scaling_;
    LOG_DBG3(SMD) << "WBTE: l[" << (offset + i) << "]=" << x_l[offset + i] << " u[" << (offset + i) << "]=" << x_u[offset + i];
  }
}

inline unsigned int ShapeMapDesign::GetNumberOfVariables() const
{
  return aux_design_.GetSize() + shape_param_.GetSize();
}


inline BaseDesignElement* ShapeMapDesign::GetDesignElement(unsigned int idx)
{
  if(idx < aux_design_.GetSize())
    return AuxDesign::GetDesignElement(idx);
  else
    return &shape_param_[idx - aux_design_.GetSize()];
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
    <shapeParam type="node" dof="x" lower="0" upper=".5" initial=".25"/>

    <shapeParam type="profile" fixed=".1" />
  </shapeMap> */
 void ShapeMapDesign::SetupShapeDesign(PtrParamNode pn)
{
   // read shapeParam to shape_
   assert(beta_ >= 0.0); // shall be already set
   ParamNodeList list = pn->GetList("shapeParam");
   LOG_DBG(SMD) << "SSD: childs of shapeParam: " << list.GetSize();
   assert(!list.IsEmpty());
   shape_.Resize(list.GetSize());
   for(unsigned int i = 0; i < list.GetSize(); i++)
     shape_[i].Init(list[i], i); // empty constructor due to StdVector/(

   // check rhos which should be already be set in DesignSpace::data
   assert(GetRegionIds().GetSize() == 1); // more is not implemented yet
   StdVector<int> elem_to_idx;
   StdVector<int> idx_to_elem;
   StdVector<unsigned int> n = SetupLexicographicMesh(domain->GetGrid(), GetRegionIds().First(), elem_to_idx, idx_to_elem);
   nx_ = n[0];
   ny_ = n[1];
   nz_ = n[2];
   assert(data.GetSize() == nx_ * ny_ * nz_); // DesignSpace::data has an element for each FEM-Cell

   // set map_ to map from shape_param to DesignSpace::data, fill with shape_param_
   map_.Resize(data.GetSize());
   for(unsigned int i = 0; i < map_.GetSize(); i++)
     map_[i].param.Reserve(shape_.GetSize());

   // setup shape_param_
   assert(n[0] == n[1]); // up to now
   assert(n[2] == 1); // 2D
   StdVector<ShapeParam*> shape_x = FindShape(NODE, 0);
   StdVector<ShapeParam*> shape_y = FindShape(NODE, 1);

   shape_param_.Reserve((nx_+1) * shape_x.GetSize()+ (ny_+1) * shape_y.GetSize()); // one node more than elements
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

   // setup map_
   map_.Resize(data.GetSize());
   for(unsigned int i = 0, n = map_.GetSize(); i < n; i++)
   {
     map_[i].rho = &(data[Find(idx_to_elem[i])]); // is very fast and gives a layer for arbitrary element ordering in the mesh
     map_[i].param.Reserve(2 * shape_.GetSize()); // each design node connects to two density elements
     map_[i].ip_param_idx.Resize(order_ * order_);
     map_[i].ip_param_idx.Init(-1);
   }

   for(unsigned int i = 0, n = shape_param_.GetSize(); i < n; i++)
   {
     ShapeParamElement& spe = shape_param_[i];
     if(spe.dof == 0) // x
     {
       int y = spe.idx[1];
       assert(y >= 0);

       // exclude the case where the y node is at the upper boundary (one node more than elements!)
       for(unsigned int x = 0; x < nx_; x++){
         if(y < (int) ny_) // are we the topmost node ontop of the last element
           map_[DensityIdx(x, y)].param.Push_back(&spe);
         if(y-1 >= (int) 0) // node 2 participates to the lower element (2-1) and the upper element (2)
           map_[DensityIdx(x, y-1)].param.Push_back(&spe);
       }
     }

     if(spe.dof == 1) // y
     {
       int x = spe.idx[0];
       for(unsigned y = 0; y < ny_; y++)  {
         if(x < (int) nx_)
           map_[DensityIdx(x, y)].param.Push_back(&spe);
         if(x-1 >= 0)
           map_[DensityIdx(x-1, y)].param.Push_back(&spe);
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
     assert(item.param.GetSize() % 2 == 0); // we expect to have pairs
     for(unsigned int s = 0; s < item.param.GetSize(); s+=2)
       if(ApproxMaxRho(item.param[s], item.param[s+1], coords) > 1e-10)
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
           double t = Eval(item.param[shapes[si]], item.param[shapes[si]+1], coords, ip_x, ip_y, false);
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

 void ShapeMapDesign::MapShapeGradient(bool obj)
 {
   assert(mapped_design_ == design_id); // the map_::Item stuff needs to be up to date
   LOG_DBG(SMD) << "MSG: obj=" << obj << " di=" << design_id << " md=" << mapped_design_ << " mog=" << mapped_obj_gradient_ << " mcg=" << mapped_constr_gradient_;
   // fixme! We do the double job of dtanh_da for obj and for constraints. However if we do it common
   // Optimization::EvalObjectiveConstraints() triggers MapShapeGradient() via WriteGradientToExtern() but rho::constraintGrad might not be set yet
   // the mapping goes shape->rho(->state), e.g. to compliance
   // from DesignSpace and ErsatzMaterial we have d_function/d_rho.
   // By chain rule we get d_function/d_shape by summing up for each shape_var the corresponding d_function/d_rho.
   // Note that we do this based on each integration point! It counts for each ip the shape_var which has larges rho(ip).
   // This was set in MapShapeToDensity() to Item::ip_param_id
   // ip_param_id is -1 if rho indicates zero gradient d_mapping/d_shape

   // shall we store max_grad for output? -1 if not
   int res_idx_da = GetSpecialResultIndex(DesignElement::DENSITY, DesignElement::SHAPE_MAP_GRAD, DesignElement::SM_NODE);

   Grid* grid = domain->GetGrid();
   // within the element coordinates we perform the integration
   Matrix<double> coords;
   for(unsigned int r = 0, n = map_.GetSize(); r < n; r++) // traverse all rho design elements
   {
     Item& item = map_[r];
     DesignElement* de = item.rho;
     grid->GetElemNodesCoord(coords, de->elem->connect, false); // no deformed mesh
     double log_da = 0.0;

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
           ShapeParamElement* s1 = item.param[shape_idx];
           ShapeParamElement* s2 = item.param[shape_idx+1];

           double da = Eval(s1, s2, coords, ip_x, ip_y, true); // dtanh_da
           // normalize FIXME! For a reason I don't understand the 0.5 makes the snopt gradient check work?!
           double da_norm = 0.5 * da / (order_ * order_);
           log_da += da_norm;

           LOG_DBG3(SMD) << "MSG: -> el=" << de->elem->elemNum << " rho=" << de->GetPlainDesignValue() << " ip_x=" << ip_x
                         << " ip_y=" << ip_y << " ip_idx=" << ip_idx << " si=" << shape_idx << " s1=" << s1->GetIndex()
                         << " s2=" << s2->GetIndex() << " da=" << da << " da_norm=" << da_norm;

           assert(opt_ != NULL);
           assert(opt_->objectives.data.GetSize() == s1->costGradient.GetSize());
           if(obj)
           {
             for(unsigned int c = 0; c < opt_->objectives.data.GetSize(); c++)
             {
               const Objective* f = opt_->objectives.data[c];
               LOG_DBG3(SMD) << "MSG: c=" << f->GetName() << " s1=" << s1->GetIndex() << " rho_grad=" << de->costGradient[c]
                             << " old =" << s1->costGradient[c] << " -> " << s1->costGradient[c] * (1.0 + f->GetPenalty() * de->costGradient[c] * da_norm);
               LOG_DBG3(SMD) << "MSG: c=" << f->GetName() << " s2=" << s2->GetIndex() << " rho_grad=" << de->costGradient[c]
                             << " old =" << s2->costGradient[c] << " -> " << s2->costGradient[c] * (1.0 + f->GetPenalty() * de->costGradient[c] * da_norm);
               s1->AddGradient(opt_->objectives.data[c], de->costGradient[c] * da_norm);
               s2->AddGradient(opt_->objectives.data[c], de->costGradient[c] * da_norm);
             }
           }
           else
           {
             assert(opt_->constraints.active.GetSize() == s1->constraintGradient.GetSize());
             for(unsigned int g = 0; g < opt_->constraints.active.GetSize(); g++)
             {
               assert(!opt_->constraints.active[g]->IsLocalCondition()); // FIXME local functions!!
               const Condition* gf = opt_->constraints.active[g];
               LOG_DBG3(SMD) << "MSG: g=" << gf->ToString() << " s1=" << s1->GetIndex() << " rho_grad=" << de->constraintGradient[g]
                             << " old=" << s1->constraintGradient[g] << " -> " << s1->constraintGradient[g] * (1.0 + de->constraintGradient[g] * da_norm);
               LOG_DBG3(SMD) << "MSG: g=" << gf->ToString() << " s2=" << s2->GetIndex() << " rho_grad=" << de->constraintGradient[g]
                             << " old=" << s2->constraintGradient[g] << " -> " << s2->constraintGradient[g] * (1.0 + de->constraintGradient[g] * da_norm);
               s1->AddGradient(opt_->constraints.active[g], de->constraintGradient[g] * da_norm);
               s2->AddGradient(opt_->constraints.active[g], de->constraintGradient[g] * da_norm);
             }
           }
         } // end ip_y
       } // end ip_x
       // normalize by integration points.
     }
     else
     {
       LOG_DBG2(SMD) << "MSG: obj=" << obj << " el=" << de->elem->elemNum << " item.ip_param_idx[0]=" << item.ip_param_idx[0] << " rho=" << de->GetPlainDesignValue() << " sum da=" << log_da;
       // when one ip_param_idx is -1, all are -1.
       // the reason is, that all are -1 if ApproxMaxRho() to small for all shapes. Otherwise at least one shape has the ips
       assert(item.ip_param_idx.Sum() == -1 * (int) (order_*order_));
     }
     LOG_DBG2(SMD) << "MSG: el=" << de->elem->elemNum << " rho=" << de->GetPlainDesignValue() << " sum da=" << log_da;
     if(res_idx_da >= 0)
       de->specialResult[res_idx_da] = log_da;

   } // end loop over density elements

   if(obj)
     mapped_obj_gradient_ = design_id;
   else
     mapped_constr_gradient_ = design_id;
 }

 inline double ShapeMapDesign::Eval(const ShapeParamElement* s1, const ShapeParamElement* s2, const Matrix<double>& coords, unsigned int ip_x, unsigned int ip_y, bool grad) const
 {
   assert(dim_ == 2);
   assert(s1->dof == s2->dof);
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
   double a1 = s1->GetDesign();
   double a2 = s2->GetDesign();

   double xy = start + (dof == 0 ? ip_x : ip_y) /(order_-1.) * (end-start); // for dof=0 (x) xy is x and might be far away from a
   double a  = a1 + (dof == 0 ? ip_y : ip_x)/(order_-1.) * (a2-a1);         // for dof=1 (c) a is y and with a1=a2 we have the same value for a
   double val = grad ? d_tanh_da(xy, a, 0.05) : tanh(xy, a, 0.05);
   LOG_DBG3(SMD) << "E: s1=" << s1->GetIndex() << " s2=" << s2->GetIndex() << " dof=" << dof << " start=" << start << " end=" << end  << " a=" << a1 << "..." << a2
                 << " ip_x=" << ip_x << " ip_y=" << ip_y << " xy=" << xy << " a=" << a << " grad:" << grad << " -> " << val;
   return val;
 }

 inline double ShapeMapDesign::ApproxMaxRho(const ShapeParamElement* s1, const ShapeParamElement* s2, const Matrix<double>& coords) const
 {
   assert(dim_ == 2);
   double max        = Eval(s1, s2, coords, 0,               0, false); // false=no derivative
   max = std::max(max, Eval(s1, s2, coords, 0,        order_-1, false));
   max = std::max(max, Eval(s1, s2, coords, order_-1,        0, false));
   max = std::max(max, Eval(s1, s2, coords, order_-1, order_-1, false));
   return max;
 }


 inline double ShapeMapDesign::tanh(double x, double a, double w) const
 {
   if(x <= a)
     return 1.0 - 1/(std::exp(2*beta_*(x-a+w))+1);
  else
    return 1/(std::exp(2*beta_*(x-a-w))+1);
 }

 inline double ShapeMapDesign::d_tanh_da(double x, double a, double w) const
 {
   //if x <= a
   //  f = -1* (exp(2*beta*(x-a+d)) + 1)^-2 *2*beta*exp(2*beta*(x-a+d));
   //else
   ///  f = (exp(2*beta*(x-a-d))+1)^-2 *2*beta*exp(2*beta*(x-a-d));
   // end

   if(x <= a)
     return -1./((std::exp(2*beta_*(x-a+w))+1) * (std::exp(2*beta_*(x-a+w))+1)) * 2*beta_*std::exp(2*beta_*(x-a+w));
  else
    return 1/((std::exp(2*beta_*(x-a-w))+1) * (std::exp(2*beta_*(x-a-w))+1)) * 2*beta_*std::exp(2*beta_*(x-a-w)) ;
 }

 void ShapeMapDesign::DumpMap()
 {
    for(unsigned int y = 0; y < ny_; y++)
    {
      for(unsigned int x = 0; x < nx_; x++)
      {
        Item i = map_[DensityIdx(x, y)];
        std::cout << "y=" << y << " x=" << x << " elidx=" << DensityIdx(x, y);
        for(unsigned int s = 0; s < i.param.GetSize(); s++)
           std::cout << " [nr=" << i.param[s]->GetIndex() << " dof=" << i.param[s]->dof << " free=" << i.param[s]->idx[1-i.param[s]->dof] << "]"; // works only for 2D
        std::cout << std::endl;
      }
      std::cout << std::endl;
    }
 }

void ShapeMapDesign::CreateShapeVariable(const ShapeParam* param, int free)
{
  assert(param->type == NODE);
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
  shape_param_.Push_back(ShapeParamElement(BaseDesignElement::NODE, shape_param_.GetSize()));
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

  LOG_DBG2(SMD) << "CSV el=" << (shape_param_.GetSize() - 1) << " dof=" << spe.dof << " free=" << free << " d=" << spe.GetDesign() << " coord=" << spe.coord.ToString();
}

void ShapeMapDesign::PostInit(int objectives, int constraints)
{
  AuxDesign::PostInit(objectives, constraints);
  for(unsigned int i = 0, n = shape_param_.GetSize(); i < n; i++)
    shape_param_[i].PostInit(objectives, constraints);

  if(domain->GetOptimization() != NULL)
    opt_ = domain->GetOptimization();
}

StdVector<ShapeMapDesign::ShapeParam*> ShapeMapDesign::FindShape(Type type, int dof)
{
   StdVector<ShapeParam*> res;
   for(unsigned int i = 0; i < shape_.GetSize(); i++)
     if(shape_[i].type == type && shape_[i].dof == dof)
       res.Push_back(&shape_[i]);
   return res;
 }

 void ShapeMapDesign::ShapeParam::Init(PtrParamNode pn, unsigned int idx_)
 {
   type = ShapeMapDesign::type.Parse(pn->Get("type")->As<std::string>());
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



