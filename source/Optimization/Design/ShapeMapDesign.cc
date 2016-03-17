/*
 * ShapeMapDesign.cc
 *
 *  Created on: Mar 11, 2016
 *      Author: fwein
 */
#include "Optimization/Design/ShapeMapDesign.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"

namespace CoupledField {

DECLARE_LOG(SMD)
DEFINE_LOG(SMD, "shapeMapDesign")


ShapeMapDesign::ShapeMapDesign(StdVector<RegionIdType>& regionIds, PtrParamNode pn, ErsatzMaterial::Method method)
: AuxDesign(regionIds, pn, method), order_(10)
{
  this->beta_ = pn->Get("shapeMap/beta")->As<double>();
  this->dim_ = domain->GetGrid()->GetDim();

  if(order_ <= 3)
    info_->Get(ParamNode::HEADER)->Get(ParamNode::WARNING)->SetValue("low integration order for shape map");

  SetupShapeDesign(pn->Get("shapeMap"));

  LOG_DBG(SMD) << "SMP rho_desig=" << data.GetSize();
  LOG_DBG(SMD) << "regions: " << regionIds.ToString();
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
     shape_[i].Init(list[i]); // empty constructor due to StdVector/(

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
   StdVector<ShapeParam> shape_x = FindShape(NODE, 0);
   StdVector<ShapeParam> shape_y = FindShape(NODE, 1);

   shape_param_.Reserve((nx_+1) * shape_x.GetSize()+ (ny_+1) * shape_y.GetSize()); // one node more than elements
   assert(shape_param_.Capacity() > 0);
   for(unsigned int s = 0; s < shape_x.GetSize(); s++)
     for(unsigned int y = 0; y < ny_+1; y++)
       CreateShapeVariable(shape_x[s], y);
   for(unsigned int s = 0; s < shape_y.GetSize(); s++)
     for(unsigned int x = 0; x < nx_+1; x++)
       CreateShapeVariable(shape_y[s], x);

   // setup map_
   map_.Resize(data.GetSize());
   for(unsigned int i = 0, n = map_.GetSize(); i < n; i++)  {
     map_[i].rho = &(data[Find(idx_to_elem[i])]); // is very fast and gives a layer for arbitrary element ordering in the mesh
     map_[i].param.Reserve(2 * shape_.GetSize()); // each design node connects to two density elements
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
   Grid* grid = domain->GetGrid();
   // here we store the value for each integration point.
   Vector<double> v(order_ * order_);
   Matrix<double> coords;

   for(unsigned int r = 0, n = data.GetSize(); r < n; r++)
   {
     v.Init(0.0); // reset such we can do max
     const Item& item = map_[r];
     DesignElement* de = item.rho;
     grid->GetElemNodesCoord(coords, de->elem->connect, false); // no
     //LOG_DBG3(SMD) << "MS2D: r=" << r << " coord=" << coords.ToString(2, false);
     assert(item.param.GetSize() % 2 == 0); // we expect to have pairs
     for(unsigned int s = 0; s < item.param.GetSize(); s+=2)
     {
       ShapeParamElement* s1 = item.param[s];
       ShapeParamElement* s2 = item.param[s+1];
       assert(dim_ == 2);
       assert(s1->dof == s2->dof);
       // element position
       double start = s1->dof == 0 ? coords[0][0] : coords[1][0];
       double end   = s1->dof == 0 ? coords[0][1] : coords[1][3];
       // the parameters
       double a1 = s1->GetDesign();
       double a2 = s2->GetDesign();

       LOG_DBG3(SMD) << "MS2D: r=" << r << " s=" << s << " a=" << a1 << "..." << a2 << " dof=" << s1->dof << " var=" << (1-s1->dof) << " start=" << start << " end=" << end;

       for(unsigned int c = 0; c < order_; c++)
       {
         for(unsigned int p = 0; p < order_; p++)
         {
           double xy = start +  c/(order_-1.) * (end-start); // for dof=0 (x) xy is x and might be far away from a
           double a  = a1 + p/(order_-1.) * (a2-a1);         // for dof=1 (c) a is y and with a1=a2 we have the same value for a
           double val = tanh(xy, a, 0.05);
           // v represents a matrix but the view depends on the dof. With dof==0, the fast variable shall be x for col_major
           unsigned int idx = s1->dof == 0 ? p * order_  + c : c * order_ + p;
           v[idx] = std::max(v[idx], val);
           LOG_DBG3(SMD) << "MS2D: -> c=" << c << " p=" << p << " idx=" << idx << " xy=" << xy << " a=" << a << " v=" << val << " -> " << v[idx];
         }
       }

       assert(v.Avg() >= 0 && v.Avg() <= 1);
       de->SetDesign(de->GetLowerBound() + (de->GetUpperBound() - de->GetLowerBound()) * v.Avg()); // we assume 0 <= v <= 1
       LOG_DBG3(SMD) << "MS2D: -> el=" << de->elem->elemNum << " -> avg=" << de->GetPlainDesignValue();
       // Matrix<double> m;
       // m.Assign(v, order_, order_, true);
       // LOG_DBG3(SMD) << m.ToString();
      }
   }
 }


 inline double ShapeMapDesign::tanh(double x, double a, double w) const
 {
   if(x <= a)
     return 1.0 - 1/(std::exp(2*beta_*(x-a+w))+1);
  else
    return 1/(std::exp(2*beta_*(x-a-w))+1);
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

void ShapeMapDesign::CreateShapeVariable(ShapeParam& param, int free)
{
  assert(param.type == NODE);
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
  spe.SetLowerBound(param.lower);
  spe.SetUpperBound(param.upper);
  spe.SetDesign(param.value);
  spe.dof = param.dof;
  if(spe.dof == 0) {
    spe.coord[1] = (double) free / ny_; // free is max ny_
    spe.idx[1] = free;
  }
  if(spe.dof == 1) {
    spe.coord[0] = (double) free / nx_;
    spe.idx[0] = free;
  }

  LOG_DBG2(SMD) << "CSV el=" << (shape_param_.GetSize() - 1) << " dof=" << spe.dof << " free=" << free << " d=" << spe.GetDesign() << " coord=" << spe.coord.ToString();
}

StdVector<ShapeMapDesign::ShapeParam> ShapeMapDesign::FindShape(Type type, int dof)
{
   StdVector<ShapeParam> res;
   for(unsigned int i = 0; i < shape_.GetSize(); i++)
     if(shape_[i].type == type && shape_[i].dof == dof)
       res.Push_back(shape_[i]);
   return res;
 }

 void ShapeMapDesign::ShapeParam::Init(PtrParamNode pn)
 {
   type = ShapeMapDesign::type.Parse(pn->Get("type")->As<std::string>());
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


} // end of namespace



