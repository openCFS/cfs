#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION

#include "Optimization/Design/FeatureMappingDesign.hh"
#include "DataInOut/Logging/LogConfigurator.hh"

using std::string;
using std::pair;
using std::make_pair;
using std::to_string;

namespace CoupledField {

DEFINE_LOG(fm, "featureMapping")

FeatureMappingDesign::FeatureMappingDesign(StdVector<RegionIdType>& regionIds, PtrParamNode empn, ErsatzMaterial::Method method)
: FeaturedDesign(regionIds, empn, method)
{
  setup_timer_->Start();

  assert(method == ErsatzMaterial::FEATURE_MAPPING);

  PtrParamNode pn = empn->Get(ErsatzMaterial::method.ToString(method)); 

  combine_ = combine.Parse(pn->Get("combine")->As<string>());
  if(combine_ != MAX)
  {
    if(!pn->Has("p"))
      throw Exception("featureMapping attribute 'p' is mandatory for combine " + combine.ToString(combine_));
    cmb_param = pn->Get("p")->As<double>();
  } 
  boundary_ = POLY;

  //if (method == ErsatzMaterial::SPAGHETTI_PARAM_MAT)
  //  orientation_ = orientation.Parse(pn->Get("orientation")->As<string>());

  order = pn->Get("order")->As<unsigned int>();
  transition = pn->Get("transition")->As<double>(); // make optional

  // map_
  SetupMapping();

  // setup pills with pills and shape_param_ and opt_shape_param_
  SetupDesign(pn);

  // the framework for efficient distance field determination - we need to know the number of features 
  SetupFixedIntegrationPoints();

  // get rhomin from first density element and check transfer function
  assert(!data.IsEmpty());
  const DesignElement& de = data.First();
  rhomin = de.GetLowerBound(); // at least in python use for all
  rhomax = de.GetUpperBound();

  bnd_fnc = std::make_unique<Poly>(transition, rhomin);
  switch(combine_)
  {
  case MAX:
    cmb_fnc = std::make_unique<Max>();
    if(features_.GetSize() > 1)
      info_->Get("features/combine")->SetWarning("Using non differentiable combine=max with multiple features is not recommended, use p-norm, KS or softmax");
    break;
  case P_NORM:
    cmb_fnc = std::make_unique<P_Norm>(cmb_param);
    break;
  case KS:
    cmb_fnc = std::make_unique<KreisselmeierSteinhauser>(cmb_param);
    break;
  case SOFTMAX:
    cmb_fnc = std::make_unique<SoftMax>(cmb_param);
    break;
  default:
    assert(false);
    break;
  }
  assert(mapped_design_ != design_id);

  MapFeatureToDensity(); // only so late because of python -> calls PythonUpdateSpaghetti()

  LOG_DBG(fm) << "FMD data=" << data.GetSize() << " aux=" << aux_design_.GetSize() << " shape=" << opt_shape_param_.GetSize() << " -> N=" << GetNumberOfVariables();
  assert(GetNumberOfVariables() > 0);

  setup_timer_->Stop(); // python and map are subtimers, so count them here
}

void FeatureMappingDesign::ToInfo(ErsatzMaterial* em)
{
  FeaturedDesign::ToInfo(em);

  PtrParamNode in = info_->Get("features");
  assert(combine_ != NO_COMBINE);
  if(combine_ != MAX)
    in->Get("combine/param")->SetValue(cmb_param);
  in->Get("boundary/transition")->SetValue(transition);
}


void FeatureMappingDesign::SetupMapping()
{
  FeaturedDesign::SetupMapping();
  assert(map_.GetSize() == nx_ * ny_); 

  for(Item& item : map_)
  {
    ItemIP* item_ip = new ItemIP(); // we use ItemIP as extension
    if(order > 1)
      item_ip->corner.Reserve(dim_ == 2 ? 4 : 8); // max 4 corners for 2D, 8 for 3D
    item.extension = item_ip; // ~Item() will delete this
  } 
}

inline FeatureMappingDesign::ItemIP* FeatureMappingDesign::GetItemIP(unsigned int item_idx)
{
  assert(map_[item_idx].extension != nullptr);
  assert(dynamic_cast<ItemIP*>(map_[item_idx].extension) != nullptr);
  return dynamic_cast<ItemIP*>(map_[item_idx].extension);
}

FeatureMappingDesign::IntegrationPoint* FeatureMappingDesign::SetupDesignCreateAddIP(FeatureMappingDesign::ItemIP::Storage storage, const Point& ref, FeatureMappingDesign::ItemIP* item_ip, double dx, double dy, int num)
{
  IntegrationPoint* ip = nullptr;
  if(storage == ItemIP::CORNER)
  {
    corners.Push_back(IntegrationPoint());
    ip = &corners.Last();
    assert(item_ip->corner.GetCapacity() > item_ip->corner.GetSize()); 
    item_ip->corner.Push_back(ip); // the item shares the corner point
  }
  else
  {
    ip = &(inners.emplace_front()); // inners has unpredicted size for order > 2, hence no resize allowed and we use a list
    assert(item_ip->inner.GetCapacity() > item_ip->inner.GetSize()); 
    item_ip->inner.Push_back(ip);
  }
  ip->dist.Resize(num);
  ip->part.Resize(num, FeatureVariable::NO_TIP);
  ip->loc.Set(ref[0]+dx, ref[1]+dy); 

  item_ip->rho.Resize(num);
  
  // LOG_DBG3(fm) << "SDCAI c=" << center.ToString() << " dx=" << dx << " dy=" << dy << " ext=" << item_ip->ToString();
  return ip;
}

void FeatureMappingDesign::SetupDesignAddIP(FeatureMappingDesign::ItemIP::Storage storage, FeatureMappingDesign::Item& item, FeatureMappingDesign::IntegrationPoint* ip)
{
  ItemIP* item_ip = dynamic_cast<ItemIP*>(item.extension);
  assert(item_ip != nullptr); 
  StdVector<IntegrationPoint*>& vec = (storage == ItemIP::CORNER) ? item_ip->corner : item_ip->inner;
  assert(vec.GetCapacity() > vec.GetSize()); 
  assert(ip != nullptr);
  vec.Push_back(ip); 
  LOG_DBG3(fm) << "SDCAI item=" << item.lexicographic_pos << " s=" << (storage == ItemIP::CORNER ? "corner" : "inner") << " ip=" << ip->loc.ToString() 
               << " item_ip=" << item_ip->ToString() << " vec.size=" << vec.GetSize();
}

void FeatureMappingDesign::SetupDesign(PtrParamNode pn)
{
  FeaturedDesign::SetupDesign(pn);
 
  // the order is nodes, profile, [alpha|
  for(Feature& pill : pills)
  {
    for(auto& p : pill.points) // P, [inner], Q
      for(auto& v : p)
        AddVariable(&v);

    AddVariable(&(pill.p));

    // if(pill.HasAlpha())
    //   AddVariable(&(pill.alpha));
  }

  // post process mapping
  for(Feature& pill : pills)
  {
    for(FeatureVariable* var : pill.GetAllVariables()) 
    {
      if(var->key != "")
        if(CountKey(var->key) > 1)
          throw Exception("feature key '" + var->key + "' occurs more than once. Use unique keys for each variable."); 
      if(var->map != "")
      {
        var->map_to = FindKey(var->map); // important to traverse a non-resizable vector
        if(var->map_to == nullptr)
          throw Exception("There is no feature key '" + var->map + "' to map to");
        var->SetDesign(var->map_to->GetPlainDesignValue()); // set the design value to the mapped variable
      }
    } 
  }

  // we did exact space "estimation"
  LOG_DBG(fm) << "SD: shape_param=" << shape_param_.GetSize() << " opt_shape_param=" << opt_shape_param_.GetSize() 
            << " opt_indices=" << opt_indices.GetSize() << " pills=" << pills.GetSize();
  assert(shape_param_.GetCapacity() == shape_param_.GetSize());
  assert(opt_shape_param_.GetCapacity() == opt_shape_param_.GetSize());
}

void FeatureMappingDesign::SetupFixedIntegrationPoints()  
{
  assert(dim_ == 2); 

  int num = features_.GetSize();
  assert(num >= 1);
  Grid* grid = domain->GetGrid();
  Matrix<double> coords;
  Point center;

  assert(map_.GetSize() == nx_ * ny_); 

  if(order == 1)
  {
    for(unsigned int y = 0; y < ny_; y++)
    {
      for(unsigned int x = 0; x < nx_; x++)
      {
        Item& item = map_[y * nx_ + x];
        DesignElement* de = item.elemval;
        assert(de != nullptr);
        grid->GetElemNodesCoord(coords, de->elem->connect, false); // obtain element nodal coords, no update. )
        LagrangeElemShapeMap::CalcBarycenter(center, coords, domain); // get the barycenter of the element
        LOG_DBG3(fm) << "SD y=" << y << " x=" << x << " e=" << de->elem->elemNum << " c=" << center.ToString();// << " coords=" << coords.ToString();
      
        ItemIP* ext = dynamic_cast<ItemIP*>(item.extension);
        assert(ext->inner.GetCapacity() == 0); 
        assert(ext->corner.GetCapacity() == 0); 
        ext->inner.Reserve(1);
        ext->center = center;
        SetupDesignCreateAddIP(ItemIP::INNER, center, ext, 0.0, 0.0, num); // barycenter
      }
    }
  }
  if(order >= 2) // for order > 2 we have corners fixed and inner dynamically handled in MapFeatureToDensity()
  {
    corners.Reserve((nx_+1) * (ny_+1)); // corners are the element nodes, hence larger. We don't care about the ordering

    // we traverse the map elements and set the corners - usually the right lower one 
    // Attention!!! We assume lexicographic ordering here!! 
    for(unsigned int y = 0; y < ny_; y++)
    {
      for(unsigned int x = 0; x < nx_; x++)
      {
        Item& item = map_[y * nx_ + x];
        DesignElement* de = item.elemval;
        assert(de != nullptr);
        // we are not sure about node ordering for pathological cases and therefore recompute everything from the element barycenter
        // the same time we assume that y/x ordering is lexicographic!!                                    
        grid->GetElemNodesCoord(coords, de->elem->connect, false); // obtain element nodal coords, no update. )
        LagrangeElemShapeMap::CalcBarycenter(center, coords, domain); // get the barycenter of the element
        LOG_DBG3(fm) << "SD y=" << y << " x=" << x << " e=" << de->elem->elemNum << " c=" << center.ToString();// << " coords=" << coords.ToString();
        
        ItemIP* iip = dynamic_cast<ItemIP*>(item.extension);
        iip->center = center;

        // o ----- A
        // |       |
        // |   c   | c == center
        // |       |
        // B ----- X
  
        // the "home" ip for this element is X - max three other elements share this ip
        IntegrationPoint* ip = SetupDesignCreateAddIP(ItemIP::CORNER, center, iip, dx_/2, -dx_/2, num); // right lower corner X
        if(y > 0)
          SetupDesignAddIP(ItemIP::CORNER, map_[(y-1) * nx_ + x], ip); // the item below shares the corner point X
        if(y > 0 && x < nx_-1)
          SetupDesignAddIP(ItemIP::CORNER, map_[(y-1) * nx_ + (x+1)], ip); // diagonal below/right
        if(x < nx_-1)
          SetupDesignAddIP(ItemIP::CORNER, map_[y * nx_ + (x+1)], ip); // right
        
        // special handling for upper most row  
        if(y == ny_-1) // right side but above
        {  
          // add IP A (right but above)
          ip = SetupDesignCreateAddIP(ItemIP::CORNER, center, iip, dx_/2, dx_/2, num); // right upper corner A
          if(x < nx_-1)
            SetupDesignAddIP(ItemIP::CORNER, map_[y * nx_ + (x+1)], ip); // right neighbor
        }
        // special handling vor first column
        if(x == 0) // left side 
        {
          // add IP V (below but left)
          ip = SetupDesignCreateAddIP(ItemIP::CORNER, center, iip, -dx_/2, -dx_/2, num); // left lower corner B
          if(y > 0)
            SetupDesignAddIP(ItemIP::CORNER, map_[(y-1) * nx_ + x], ip); 
        }
        if(x == 0 && y == ny_-1) // left upper corner o
        {
          // this node o is not shared by anyone
          SetupDesignCreateAddIP(ItemIP::CORNER, center, iip, -dx_/2, +dx_/2, num); // left upper corner o
        }
      } 
    } // end of y,x loop
  } // order >= 2 - corner

  // optional debug logging
  for(IntegrationPoint& ip : corners)
    LOG_DBG3(fm) << "SD corners IP loc=" << ip.loc.ToString();

 for(IntegrationPoint& ip : inners)
    LOG_DBG3(fm) << "SD inners IP loc=" << ip.loc.ToString();

  for(Item& item : map_)
    LOG_DBG3(fm) << "SD map " << item.lexicographic_pos << " ext: " << item.extension->ToString();
}

std::string FeatureMappingDesign::ItemIP::ToString()
{
  std::stringstream ss;
  
  ss << "c=" << center.ToString();
  
  ss << " corner #" << corner.GetSize() << "/" << corner.GetCapacity() << "->[";
  for(IntegrationPoint* ip : corner)
    ss << ip->loc.ToString() << " d=" << ip->dist.ToString() << ", ";
  ss << "]";
  return ss.str();
}

void FeatureMappingDesign::ItemIP::FillDist(Vector<double>& data, unsigned int num) const
{
  data.Resize(corner.GetSize() + inner.GetSize()); // no data copied
  for(unsigned int i = 0; i < corner.GetSize(); i++)
    data[i]=corner[i]->dist[num];
  for(unsigned int i = 0; i < inner.GetSize(); i++)
    data[corner.GetSize()+i]=inner[i]->dist[num];
}

void FeatureMappingDesign::ItemIP::GetAllIP(StdVector<IntegrationPoint*>& out)
{
  out.ResizeNoCopy(corner.GetSize() + inner.GetSize());
  for(unsigned int i = 0; i < corner.GetSize(); i++)
    out[i]=corner[i];
  for(unsigned int i = 0; i < inner.GetSize(); i++)
    out[corner.GetSize()+i]=inner[i];
}

StdVector<FeatureMappingDesign::IntegrationPoint*> FeatureMappingDesign::ItemIP::GetAllIP()
{
  StdVector<IntegrationPoint*> out;
  GetAllIP(out);
  return out;
}

void FeatureMappingDesign::MapFeatureToDensity()
{
  mapping_timer_->Start();
  LOG_DBG(fm) << "MFTD called";

  assert(!map_.IsEmpty());

  // the DesignElement might be updated, update the Feature internal representation
  for(Feature* f : features_)
    f->Update();
  
  // for order=1 we have a fixed ip in the barycenter for each Item in the constant inners list
  // for order=2 we have fixed nodal ip in inners and refer to Item - constant
  // for order>2 we clear here inners and keep corners constant


  // We first make inners/corners to determine the distances for each Item
  // when we next traverse map_, we know where we need higher integration
  // finally, when we have all IP with distances, we can compute the density
  if(order == 1)
  {
    for(IntegrationPoint& ip : inners)
    {
      for(unsigned int i = 0; i < pills.GetSize(); i++)
        ip.dist[i] = pills[i].Distance(ip.loc, &(ip.part[i]));
      LOG_DBG3(fm) << "MFTD barycenter IP: " << ip.loc.ToString() << " d=" << ip.dist.ToString();
    } 
  }
  else
  {
    for(IntegrationPoint& ip : corners) // for order >= 2
    {
      for(unsigned int i = 0; i < pills.GetSize(); i++)
        ip.dist[i] = pills[i].Distance(ip.loc, &(ip.part[i]));
      LOG_DBG3(fm) << "MFTD corner IP: " << ip.loc.ToString() << " d=" << ip.dist.ToString();
    } 

    if(order > 2)
    {
      // we do not know yet the shape to rho mapping. 
      // we store an InterpolationPoint ip in the inners forward_list and point to them from the ItemIP::inner vector
      // this is necessary as we don't the final number of elements and must not resize a vector.
      // 1) we check for an item (mesh element) if it will become gray. We do this by using the above calculated corner distances
      // 2) we add inner integration points for this element. The difficulty is, that each ip is shared by up to 4 other items.
      //    The difficulty is to handle that when we create the ip, we add it to the other involved items.
      // 3) we calculate the distance via the inners list such that the order * order ip of the element (of which are 4 corners) are present
      // 4) we can now apply the boundary function (polynomial) to the distances (by feature) and integrate the rho for each feature
      // 6) finally, we compute the density by combining the rho values for each feature
      
      // we clear the inner points, as we will add them dynamically
      // we reserve the space for all inner points and use the capcity as flag to fill them dynamically
      // we do this in advance to not have to call IsConstant() when we check the two edges below
      for(Item& item : map_)
      {
        ItemIP* ext = dynamic_cast<ItemIP*>(item.extension);
        ext->inner.Clear(false); // clear and do not keep capacity
        assert(ext->inner.GetCapacity() == 0); 
        if(!bnd_fnc->IsConstant(ext))
          ext->inner.Reserve(order * order - 4); // we reserve the space for all inner points, but not the corners

        LOG_DBG3(fm) << "MFTD reserve inner order=" << order << " item=" << item.lexicographic_pos 
                  << " bnd_fct=" << (bnd_fnc->IsConstant(ext) ? "constant" : "non-constant") 
                  << " inner capacity=" << ext->inner.GetCapacity();
      }
      inners.clear(); // we cleared all inner and nothing points to the outdated list any more

      // we traverse the corners and create inner points for each Item where it is necessary
      int num = features_.GetSize();
      Point base((int) dim_); // the left lower corner of the element
      assert(map_.GetSize() == nx_ * ny_);
      for(unsigned int y = 0; y < ny_; y++)
      {
        for(unsigned int x = 0; x < nx_; x++)
        {
          Item& item = map_[y * nx_ + x];
          ItemIP* ext = dynamic_cast<ItemIP*>(item.extension);
        
          if(ext->inner.GetCapacity() > 0) // this is our marker
          {
            // G --H-- I
            // |       |
            // D   E   F 
            // |       |
            // A --B-- C
  
            base.Set(ext->center[0]-dx_/2, ext->center[1]-dx_/2); // A left lower corner of the element
            // we think from the left lower corner and for inner elements set all IP which are not left and lower edge
            // hence, when we are inner, not constant, and our neighbors are also not constant, we set E,F,H (I is corner) only
            // because the lower horizontal edge B is shared from the element below  (A,C are corner) and D from the element left (A,G are corner)
            // but we share only, if the lower or left element is not constant, i.e. gets inner points
            for(unsigned int iy = 0; iy < order; iy++) 
            {
              for(unsigned int ix = 0; ix < order; ix++) 
              {
                // skip corner integration points (A,C,G,I) which are aready set and constant
                if((ix == 0 && iy == 0) || (ix == order-1 && iy == order-1) || (ix == 0 && iy == order-1) || (ix == order-1 && iy == 0))
                  continue; // skip corners
              
                // when the element below is not constant, we get B from it and continue - no corners left (no sharing across corners)
                if(iy == 0 && y > 0 && GetItemIP((y-1) * nx_ + x)->inner.GetCapacity() > 0) {
                  LOG_DBG3(fm) << "MFTD item=" << item.lexicographic_pos << " ix=" << ix << " iy="  << iy << " -> (" << base[0]+ix*(dx_/(order-1)) << "," << base[1]+iy*(dx_/(order-1)) << ") "
                               << " element below item=" << map_[(y-1) * nx_ + x].lexicographic_pos << " capacity=" << GetItemIP((y-1) * nx_ + x)->inner.GetCapacity()
                               << " will be integrated and we get B from it our inner=" << ext->inner.GetSize();
                  continue; // skip B, as it is shared with the element below
                }
                  
                if(ix == 0 && x > 0 && GetItemIP(y * nx_ + (x-1))->inner.GetCapacity() > 0) {
                  LOG_DBG3(fm) << "MFTD item=" << item.lexicographic_pos << " ix=" << ix << " iy="  << iy << " -> (" << base[0]+ix*(dx_/(order-1)) << "," << base[1]+iy*(dx_/(order-1)) << ") "
                               << " element left item=" << map_[y * nx_ + (x-1)].lexicographic_pos << " capacity=" << GetItemIP(y * nx_ + (x-1))->inner.GetCapacity()
                               << " will be integrated and we get D from it our inner=" << ext->inner.GetSize();
                  continue; // skip D, as it is shared with the element left
                }

                // actually integrate integration point
                IntegrationPoint* ip = SetupDesignCreateAddIP(ItemIP::INNER, base, ext, ix*(dx_/(order-1)), iy*(dx_/(order-1)), num);
                LOG_DBG3(fm) << "MFTD item=" << item.lexicographic_pos << " ix=" << ix << " iy=" << iy << " ip=" << ip->loc.ToString() << " created and added to inner=" << ext->inner.GetSize()
                             << " capacity=" << ext->inner.GetCapacity();

                // shall we share the upper edge (G),H,I with the upper element? Because of this not parallelizable
                if(iy == order-1 && y < ny_-1)
                {
                  Item& above = map_[(y+1) * nx_ + x]; 
                  if(dynamic_cast<ItemIP*>(above.extension)->inner.GetCapacity() > 0) { // upper edge. The capacity is our marker for non-constant elements
                    LOG_DBG3(fm) << "MFTD item=" << item.lexicographic_pos << " ix=" << ix << " iy=" << iy << " share ip " << ip->loc.ToString() << " from " << item.lexicographic_pos 
                                  << " to above " << above.lexicographic_pos << " old inner=" << GetItemIP((y+1) * nx_ + x)->inner.GetSize();
                    SetupDesignAddIP(ItemIP::INNER, above, ip); // share with upper element
                  }
                }

                // shall we share the right edge C,F,I with the right element?
                if(ix == order-1 && x < nx_-1)
                {
                  Item& right = map_[y * nx_ + (x+1)];
                  if(dynamic_cast<ItemIP*>(right.extension)->inner.GetCapacity() > 0) { // right edge
                    LOG_DBG3(fm) << "MFTD tem=" << item.lexicographic_pos << " ix=" << ix << " iy=" << iy << " share ip " << ip->loc.ToString() << " from " << item.lexicographic_pos 
                                 << " with right " << right.lexicographic_pos << " old inner=" << GetItemIP(y * nx_ + (x+1))->inner.GetSize();;
                    SetupDesignAddIP(ItemIP::INNER, right, ip); // share with right element
                  }
                }
              }
            } // end iy,ix order loop   
          } // end filling current element
        } // end loop of item
      } // end of y,x loop

      // the dynamic inner points are set, we can now finally compute the distances - parallelizable  
      for(IntegrationPoint& ip : inners)
        for(unsigned int i = 0; i < pills.GetSize(); i++)
          ip.dist[i] = pills[i].Distance(ip.loc, &(ip.part[i]));

      LOG_DBG(fm) << "MFTD order=" << order << " inners=" << std::distance(inners.begin(), inners.end());
      for(IntegrationPoint& ip : inners) 
        LOG_DBG3(fm) << "MFTD inners loc=" << ip.loc.ToString();
    } // end of order > 2
  } // end of order >= 2

  Vector<double> all_dist; // all dist from ItemIP::corner and inner for a given feature
  Vector<double> bnd;      // all_dist mapped by boundary function

  LOG_DBG(fm) << "MFTD test boundary func: " << bnd_fnc->DebugLog();

  for(Item& item : map_)
  {
    ItemIP* ext = dynamic_cast<ItemIP*>(item.extension);
    LOG_DBG3(fm) << "MFTD: item=" << item.lexicographic_pos << " center=" << ext->center.ToString() << " corner=" << ext->corner.GetSize() << " inner=" << ext->inner.GetSize() << " -> " << IntegrationPoint::ToString(ext->inner);

  }

  LOG_DBG3(fm) << "MFTD: start calculating distances for all features";
  // all integration points sets and their distances calculated, we can now integrate the density
  for(Item& item : map_)
  {
    ItemIP* ext = dynamic_cast<ItemIP*>(item.extension);
    LOG_DBG3(fm) << "MFTD: item=" << item.lexicographic_pos << " center=" << ext->center.ToString() << " corner=" << ext->corner.GetSize() << " inner=" << ext->inner.GetSize() << " -> " << IntegrationPoint::ToString(ext->inner);
    assert(order <= 2 || (ext->GetAllIP().GetSize() == 4 || ext->GetAllIP().GetSize() ==std::pow(order,dim_))); // we have at least corners, but for order > 2 we have also inner points
    for(unsigned int i = 0; i < pills.GetSize(); i++)
    {
      // collect all computed distances from all element integration points by feature number
      ext->FillDist(all_dist, i); 
      LOG_DBG3(fm) << "MFTD: item=" << item.lexicographic_pos << " f=" << i << " d=" << all_dist.ToString();
      // project these distances by the boundary function (rhomin ... 1)
      bnd_fnc->Eval(all_dist, bnd);
      // integrate these projections to a single rho by feature -> before combine
      ext->rho[i] = bnd_fnc->Integrate(bnd);
      LOG_DBG3(fm) << "MFTD: item=" << item.lexicographic_pos << " c=" << ext->center.ToString() << " f=" << i << " prj:#" << bnd.GetSize() << "->" << bnd.ToString() << " -> " << ext->rho[i];
    }
    // combine for the element the rho of all features to the remaining element rho
    double rho = cmb_fnc->Eval(ext->rho);
    LOG_DBG3(fm) << "MFTD: item=" << item.lexicographic_pos << " c=" << ext->center.ToString() << " proj=" << ext->rho.ToString() << " -> rho=" << rho;
    item.elemval->SetDesign(rho);
  }

  LOG_DBG(fm)  << "MFTD: data=" << data.GetSize() << " first=" << data[0].ToString() << " simp=" << data[0].simp;
  LOG_DBG(fm)  << "MFTD: old mid=" << mapped_design_ << " current did=" << design_id;

  mapped_design_ = design_id;
  mapping_timer_->Stop();
}

void FeatureMappingDesign::MapFeatureGradient(const Function* func)
{
  // for many density based functions we could cache stuff, as all is same but d_J/d_rho_e where J is func

  assert(mapped_design_ == design_id);
  gradient_timer_->Start();

  LOG_DBG(fm) << "MFG md= " << mapped_design_;

  int N_ip = std::pow(order,dim_); // number of integration points 

  // d_J/d_s = sum_e d_J/d_rho_e * d_rho_e/d_s 
  // with d_rho_e/d_s = 1/N_ip sum_i d_H/d_d * d_d/d_s
  // the expensive part, the distance calculation is already done, therefore we parallelize only by feature, not by variable

  // this is the case for max_rho_e via p-norm: rho_e = (sum_f rho_f^p)^(1/p)
  // d_max_rho_r/d_s = (sum_f rho_f^p)^(1/p -1) * rho_f(s)^(p-1) * d_rho_f(s)/d_s 
  // d_J/d_s = sum_e d_J/d_rho_e * (sum_f rho_f^p)^(1/p -1) * rho_f(s)^(p-1) * 1/N_ip sum_i d_H/d_d * d_d/d_s

  // we contribute to d_J/d_s and in the end move it to the opt variables (e.g. considering linking)
  // fixed is 0

  // assume all features are equal - this is not opt variables. We skip later fixed optionally link later
  unsigned int num_var_by_feature = features_.First()->GetAllVariables().GetSize(); 
  Vector<double> dJ_ds(features_.GetSize() * num_var_by_feature, 0.0);

  StdVector<IntegrationPoint*> elem_ip; // elem inner + corner
  StdVector<FeatureVariable*> f_si;     // all variables of a feature
  StdVector<double> ddist_ds(num_var_by_feature); // we calculate all gradients at once and store it then in the FeatureVariables
  // for each element
  for(Item& item : map_)
  {
    ItemIP* ext = dynamic_cast<ItemIP*>(item.extension);
    // SIMP gradient: function by rho_e
    double dJ_drho_e = item.elemval->GetPlainGradient(func);

    // for each feature
    for(unsigned int f = 0; f < pills.GetSize(); f++)
    {
      // the rho_e of this feature 
      double rho_e_f = ext->rho[f];
      // skip when we have constant feature rho_e (gradient is zero)
      if(rho_e_f == 1 || rho_e_f == rhomin)
        continue;

      Feature& pill = pills[f];

      // derivative of rho_f aggregation of all features (1.0 single feature)
      double dmax_drho_f = cmb_fnc->Grad(ext->rho, f);

      // obtain all (5) variables of the feature to use it for each IP. 
      // we handle fixed below and handle connections later
      pill.GetAllVariables(f_si);
      assert(f_si.GetSize() == num_var_by_feature); // all constant

      // for all integration points of the element
      ext->GetAllIP(elem_ip);
      LOG_DBG3(fm) << "MFG: item=" << item.lexicographic_pos << " f=" << f << " inner=" << ext->inner.GetSize() << " elem_ip=" << elem_ip.GetSize() 
                   << " dJ_drho_e=" << dJ_drho_e << " dmax_drho_f=" << dmax_drho_f 
                   << " rho_e_f=" << rho_e_f << " num_var_by_feature=" << num_var_by_feature;
      assert(elem_ip.GetSize() ==4 || (int) elem_ip.GetSize() == N_ip); 
      for(IntegrationPoint* ip : elem_ip)
      {
        // now we have the actual point X and distance for each IP know which part is closest

        // derivative of boundary function (poly) by distance
        double dH_ddist = bnd_fnc->Grad(ip->dist[f]);

        // for all variables of the specific feature (P_x, P_y, Q_x, Q_y and p at once)
        // this is the most tricky part in the derivatives
        pill.GradDistance(ip->loc, ip->part[f], ddist_ds);

        for(unsigned int s = 0; s < num_var_by_feature; s++)
        {
          FeatureVariable* var = f_si[s];
          if(var->fixed) // skip fixed but process mapped variables
            continue;
          
          dJ_ds[f * num_var_by_feature + s] += dJ_drho_e * dmax_drho_f * 1/N_ip * dH_ddist * ddist_ds[s];
        }
      }
    }
  }

  // copy from dJ_ds to our feature variables - handle mapping
  for(Feature& pill : pills)
  {
    pill.GetAllVariables(f_si);
    assert(f_si.GetSize() == num_var_by_feature);
    for(unsigned int i = 0; i < f_si.GetSize(); i++)
    {
      FeatureVariable* target = f_si[i]->map_to != nullptr ? f_si[i]->map_to : f_si[i]; // map to target variable if available
      target->AddGradient(func, dJ_ds[pill.idx * num_var_by_feature + i]);
    }
  }

  assert(dJ_ds.GetSize() == features_.GetSize() * num_var_by_feature);
  assert(dJ_ds.GetSize() == shape_param_.GetSize()); 
  LOG_DBG(fm) << "MFG: func=" << func->ToString() << " dJ_ds=" << dJ_ds.GetSize() << " sp=" << shape_param_.GetSize()  << " osp=" << opt_shape_param_.GetSize();
  for(unsigned int i = 0; i < dJ_ds.GetSize(); i++)
    LOG_DBG2(fm) << "MFG: func=" << func->ToString() << " dJ_ds[" << i << "]=" << dJ_ds[i] << " var=" << dynamic_cast<FeatureVariable*>(shape_param_[i])->ToString(); 
  for(unsigned int i = 0; i < opt_shape_param_.GetSize(); i++)
    LOG_DBG2(fm) << "MFG: func=" << func->ToString() << " opt[" << i << "]=" << opt_shape_param_[i]->GetPlainGradient(func) << " var=" << dynamic_cast<FeatureVariable*>(shape_param_[i])->ToString(); 

  gradient_timer_->Stop();
}


int FeatureMappingDesign::ReadDesignFromExtern(const double* space_in, bool setAndWriteCurrent)
{
  // practically updates the optimization Variable objects within features. Set's design_id
  FeaturedDesign::ReadDesignFromExtern(space_in, setAndWriteCurrent); 

  // all variables are set, we can now perform the optional mapping to mapped variables which are not within opt_shape_param_
  for(ShapeParamElement* spe : shape_param_)
  {
    FeatureVariable* var = dynamic_cast<FeatureVariable*>(spe);
    assert(var != nullptr); // all shape_param_ are FeatureVariables  

    FeatureVariable* org = var->map_to;
    if(org) // yes, we map
    {
      assert(var->map != "" && org->key != "" && var->map == org->key);
      var->SetDesign(org->GetPlainDesignValue()); // here we could apply advanced mapping like -key, 1-key, key+5, etc.
    }
  }

  LOG_DBG(fm) << "RDFE mid=" << mapped_design_ << " did=" << design_id;

  if(mapped_design_ != design_id)
    MapFeatureToDensity();

  return design_id;
}

void FeatureMappingDesign::ReadDensityXml(PtrParamNode set, double& lower_violation, double& upper_violation)
{
  // we expect all nodes (start x, start y, end x, end y, ...) then all profiles then all normals (if any).
  // by the shape attribute we could alll collect from random order but for easy we expect this ordering
  assert(features_.GetSize() > 0);

  ParamNodeList list = set->GetList("shapeParamElement");
  if(list.GetSize() != shape_param_.GetSize())
    EXCEPTION("expect " << shape_param_.GetSize() << " 'shapeParamElement' in density.xml but got " << list.GetSize());

  for(unsigned int i = 0; i < shape_param_.GetSize(); i++)
  {
    FeatureVariable*    var = dynamic_cast<FeatureVariable*>(shape_param_[i]);
    PtrParamNode pn  = list[i];

    // do a lot of validation
    const std::string pnnr = pn->Get("nr")->As<string>();
    if(pn->Get("nr")->As<unsigned int>() != i)
      EXCEPTION("shapeParamElement nr=" << pnnr << " in density.xml has not expected 'nr' value " << i);

    if(pn->Get("shape")->As<int>() != var->feature)
      EXCEPTION("shapeParamElement nr=" << pnnr << " in density.xml has not expected 'shape' value " << var->feature);

    if(DesignElement::type.Parse(pn->Get("type")->As<string>()) != var->GetType())
      EXCEPTION("shapeParamElement nr=" << pnnr << " in density.xml has not expected 'type' value " << DesignElement::type.ToString(var->GetType()));

    if(var->GetType() == FeatureVariable::NODE)
    {
      if(ShapeParamElement::dof.Parse(pn->Get("dof")->As<string>()) != var->dof_)
        EXCEPTION("shapeParamElement nr=" << pnnr << " in density.xml has not expected 'dof' value " << ShapeParamElement::dof.ToString(var->dof_));
      if(FeatureVariable::tip_enum.Parse(pn->Get("tip")->As<string>()) != var->tip)
        EXCEPTION("shapeParamElement nr=" << pnnr << " in density.xml has not expected 'tip' value " << FeatureVariable::tip_enum.ToString(var->tip));
    }

    var->SetDesign(pn->Get("design")->As<double>());
  }

  // it shall be save to update the design_id, because we changed stuff
  design_id++;

  MapFeatureToDensity();
}

FeatureVariable* FeatureMappingDesign::FindKey(const std::string& key) const
{
  for(const Feature* f : features_)
    for(FeatureVariable* var : f->GetAllVariables()) // important that this variable is not subject to a resizable vector
      if(var->key == key)
        return var;
  return nullptr; // not found
}

int FeatureMappingDesign::CountKey(const std::string& key) const
{
  int count = 0;
  for(const Feature* f : features_)
    for(FeatureVariable* var : f->GetAllVariables())
      if(var->key == key)
        count++;
  return count;
}

void FeatureMappingDesign::SetupParsedFeatures(PtrParamNode base)
{
  StdVector<PtrParamNode> pnl = base->GetList("pill");
  pills.Resize(pnl.GetSize());
  features_.Resize(pnl.GetSize());
  for(unsigned int i = 0; i < pnl.GetSize(); i++)
  {
    pills[i].Parse(pnl[i], i);
    features_[i] = &pills[i];
  }
}

 std::string FeatureMappingDesign::IntegrationPoint::ToString(StdVector<IntegrationPoint*>& ips)
{
  std::stringstream ss;
  for(const IntegrationPoint* ip : ips)
    ss << ip->loc.ToString() + " ";
  return ss.str();
}

void FeatureMappingDesign::Pill::Update()
{
  assert(domain->GetDim() == 2);
  assert(points.GetSize() >= 2);
  assert(P.data.GetSize() >= 2); // make it 2D
  P.Set(points.First()[0].GetPlainDesignValue(), points.First()[1].GetPlainDesignValue());
  Q.Set(points.Last()[0].GetPlainDesignValue(), points.Last()[1].GetPlainDesignValue());

  length = P.Dist(Q); // || Q - P ||

  U = (Q-P) / length;
  V.Set(U[1], -U[0]);

  dx = Q[0] - P[0]; // x-component of the vector from start to end
  dy = Q[1] - P[1]; // y-component
  dp_norm = Q[0] * P[1] - Q[1] * P[0]; // end_x*start_y - end_y*start_x

  LOG_DBG(fm) << "P:U: P=" << P.ToString() << " Q=" << Q.ToString() << " l=" << length << " U=" << U.ToString() << " V=" << V.ToString() 
              << " dx= " << dx << " dy=" << dy << " dp_nom=" << dp_norm;
}

double FeatureMappingDesign::Pill::Distance(const Point& X, FeatureVariable::Tip* part) const
{
  assert(points.GetSize() >= 2);
  assert(X.data.GetSize() >= 2);
  assert(P.data.GetSize() >= 2);
  assert(Q.data.GetSize() >= 2);
  assert(dim_ == 2);

  double pval = this->p.GetPlainDesignValue();
  double dist = -1;

  // test if we are in the region of the straight bar
  // if np.dot((X - self.Q), self.U0) <= 0 and np.dot((X - self.P), self.U0) >= 0: 
  //    projected_distance = abs(np.dot((X - self.Q), self.V0)) 
  // Point XQ = X-Q;
  // LOG_DBG(fm) << "P:CD U.Dot(X-Q)=" << U.Dot(X-Q) << " U.Dot(X-P)=" << U.Dot(X-P) << " X-Q=" << XQ.ToString();
  if(length > 0 && U.Dot(X-Q) <= 0 && U.Dot(X-P) >= 0)
  {
    // distance to line segment https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line
    dist = std::abs(V.Dot(X-Q)) - pval;
    if(part) 
      *part = FeatureVariable::INNER;
    // double dl = (std::abs(dy*X[0] - dx*X[0] + dp_norm) / length) - pval;
    LOG_DBG(fm) << "P:CD X=" << X.ToString() << " p=" << pval << " P=" << P.ToString() << " Q=" << P.ToString() << " -> line=" << dist;
  }
  else
  {
    double dp =  X.Dist(P)-pval; // distance to start point
    double dq =  X.Dist(Q)-pval; 
  
    dist = dp < dq ? dp : dq;
    if(part)
      *part = dp < dq ? FeatureVariable::START : FeatureVariable::END;

    LOG_DBG(fm) << "P:CD from " << X.ToString() << " p=" << pval << " P=" << P.ToString() << "->"  << dp 
              << " Q=" << P.ToString() << "->" << dq << " tip min=" << dist << " part=" << (part ? *part : -1);
  }  
  return dist;
}

void FeatureMappingDesign::Pill::GradDistance(const Point& X, FeatureVariable::Tip part, StdVector<double>& out) const
{
  assert(part != FeatureVariable::NO_TIP);
  assert(points.GetSize() == 2);
  assert(out.GetSize() == 5);

  switch(part)
  {
  case FeatureVariable::START:
    out[0] = (P[0]-X[0])/X.Dist(P); 
    out[1] = (P[1]-X[1])/X.Dist(P);
    out[2] = 0.0;
    out[3] = 0.0;
    out[4] = -1.0; // p
    break;
  case FeatureVariable::END:
    out[0] = 0.0; 
    out[1] = 0.0;
    out[2] = (Q[0]-X[0])/X.Dist(Q);
    out[3] = (Q[1]-X[1])/X.Dist(Q);
    out[4] = -1.0; // p
    break;
  case FeatureVariable::INNER:
  {
    // wolframalpha
    // diff (V0*(X0-Q0)+V1*(X1-Q1)) by Q1
    // V0 = ((Q1-P1) / sqrt((Q0-P0)**2 + (Q1-P1)**2))
    // V1 = ((P0-Q0) / sqrt((Q0-P0)**2 + (Q1-P1)**2))
    // diff (((Q1-P1) / sqrt((Q0-P0)**2 + (Q1-P1)**2))*(X0-Q0)+((P0-Q0) / sqrt((Q0-P0)**2 + (Q1-P1)**2))*(X1-Q1)) by Q1    
    // 
    // d/dQ1(((Q1 - P1) (X0 - Q0))/sqrt((Q0 - P0)^2 + (Q1 - P1)^2) + ((P0 - Q0) (X1 - Q1))/sqrt((Q0 - P0)^2 + (Q1 - P1)^2)) 
    // = ((Q0 - P0) (P0^2 - P0 (Q0 + X0) + P1^2 - P1 (Q1 + X1) + Q0 X0 + Q1 X1))/((P0 - Q0)^2 + (P1 - Q1)^2)^(3/2)
    
    // d/dP0 =  ((Q1 - P1) (P0 (Q0 - X0) + (P1 - Q1) (Q1 - X1) - Q0^2 + Q0 X0)) / ||P-Q||^3
    // d/dP1 = -((Q0 - P0) (P0 (Q0 - X0) + (P1 - Q1) (Q1 - X1) - Q0^2 + Q0 X0)) / ||P-Q||^3
    // d/dQ0 = -((Q1 - P1) (P0^2 - P0 (Q0 + X0) + P1^2 - P1 (Q1 + X1) + Q0 X0 + Q1 X1)) / ||P-Q||^3
    // d/dQ1 =  ((Q0 - P0) (P0^2 - P0 (Q0 + X0) + P1^2 - P1 (Q1 + X1) + Q0 X0 + Q1 X1)) / ||P-Q||^3

    Point QmP = Q - P;
    Point QmX = Q - X;
    Point QpX = Q + X;

    double denom = length*length*length;
    double sg = V.Dot(X-Q) > 0 ? 1.0 : -1.0; // to differentiate the abs from std::abs(V.Dot(X-Q)) 

    double dP_common = (P[0]*QmX[0]-QmP[1]*QmX[1]-Q[0]*Q[0]+Q[0]*X[0]);
    double dQ_common = (P[0]*P[0]-P[0]*QpX[0]+P[1]*P[1]-P[1]*QpX[1]+Q.Dot(X));

    // we assume the compiler sees that 
    out[0] = sg *  (QmP[1] * dP_common) / denom; // dP_x
    out[1] = sg * -(QmP[0] * dP_common) / denom; // dP_y
    out[2] = sg * -(QmP[1] * dQ_common) / denom; // dQ_x
    out[3] = sg *  (QmP[0] * dQ_common) / denom; // dQ_y
    out[4] = -1.0; // p
    break;
  } 
  case FeatureVariable::NO_TIP:
    assert(false);
    break;
  }
}


inline void FeatureMappingDesign::BoundaryFunction::Eval(const Vector<double>& dist, Vector<double>& eval) const
{
  eval.Resize(dist.GetSize());
  for(unsigned int i = 0; i < dist.GetSize(); i++)
    eval[i] = Eval(dist[i]);
}

inline bool FeatureMappingDesign::BoundaryFunction::IsConstant(const ItemIP* item_ip) const
{
  for(const IntegrationPoint* ip : item_ip->corner)
  {
    for(double d : ip->dist)
      if(!IsConstant(d))
        return false;
  }
  return true;
}

std::string FeatureMappingDesign::BoundaryFunction::DebugLog()
{
  std::stringstream ss;
  ss << "BF:DL t=" << 2*h << " h=" << h << " rhomin=" << rhomin << ": ";

  for(double d : {2*h, h, h/2, h/4, 0.0, -h/2, -h})
    ss << " f(" << d << ")=" << Eval(d) << ", ";
  return ss.str();
}

FeatureMappingDesign::Poly::Poly(double transition, double rhomin)
{
  this->h = transition/2; // half transition zone!
  this->rhomin = rhomin;
  this->fac = 3.0*(1.0-rhomin)/4.0;
  this->bias = ((1+rhomin)/2);
  this->h3 = h*h*h;
}


inline double FeatureMappingDesign::Poly::Eval(double d) const
{
  // note that the Feature-Mapping Review paper uses phi(x) which is -dist(x). We use here distance!
  if(IsOutside(d))
    return rhomin;
  if(IsInside(d))  
    return 1;

  return fac * (((d*d*d)/(3*h3))-(d/h))+bias; // d = -phi!
}

inline double FeatureMappingDesign::Poly::Grad(double d) const
{
  if(IsConstant(d))
    return 0.0;

  return fac * (d*d - h*h)/h3; 
}

inline double FeatureMappingDesign::Max::Eval(const Vector<double>& projected) const
{
  return projected.Max();
}

inline double FeatureMappingDesign::Max::Grad(const Vector<double>& projected, unsigned int num) const
{
  // Max combination is not really differentiable as only one feature can contribute, the others are 0
  // with more features 1, it is almost random. std::max_element() find another max than our Vector::Max() 
  // but as we are not differentiable anyway, we don't care
  unsigned int idx = std::distance(projected.GetPointer(), std::max_element(projected.GetPointer(), projected.GetPointer()+projected.GetSize()));
  LOG_DBG3(fm) << "Max:G prj=" << projected.ToString() << " by " << num << " idx=" << idx << " -> " << (num == idx ? 1.0 : 0.0);
  return num == idx ? 1.0 : 0.0; // when we have 1 feature only, this is the plain rho 
}

double FeatureMappingDesign::P_Norm::Eval(const Vector<double>& projected) const
{
   double sum = 0.0;
   for(double v : projected)
     sum += std::pow(v, p);
   double res = std::pow(sum, 1.0/p);
   LOG_DBG2(fm) << "PN:E " << projected.ToString() << " p=" << p << " sum=" << sum << " res=" << res;
   return res;
}

double FeatureMappingDesign::P_Norm::Grad(const Vector<double>& projected, unsigned int num)  const
{
  assert(projected.GetSize() > num);
  double sum = 0.0;
  for(double v : projected)
    sum += std::pow(v, p);
  double res = std::pow(sum, 1.0/p - 1.0) * std::pow(projected[num],p-1); // to be multiplied with d_rho_s/d_s
  LOG_DBG3(fm) << "PN:G by " << num << " v=" << projected.ToString() << " p=" << p << " -> " << res; 
  return res;
}    

double FeatureMappingDesign::SoftMax::Eval(const Vector<double>& projected) const
{
  // (sum_f rho_f*exp(beta*rho_f)) / (sum_f exp(beta*rho_f))
  double num = 0.0;
  double denum = 0.0;

  for(double x : projected)
  {
      double exp = std::exp(beta * x);
      num += x * exp;
      denum += exp;
  }
  return num/denum;
}

double FeatureMappingDesign::SoftMax::Grad(const Vector<double>& projected, unsigned int num)  const
{
  double sum_x_exp = 0.0;
  double sum_exp = 0.0;

  for(double x : projected)
  {
      double exp = std::exp(beta * x);
      sum_x_exp += x * exp;
      sum_exp += exp;
  }
  // u = sum_x_exp (sum x_i*exp(beta_x_i))
  // v = sum_exp (sum exp(beta*x_i))
  // u' = exp(beta*x_j)+x_j*beta*exp(beta*x_j)=(1+beta*x_j)*exp(beta*x_j)
  // v' = beta*exp(beta*x_j)
  assert(projected.GetSize() > num);
  double x_j = projected[num];
  return std::exp(beta * x_j)*((1+beta*x_j)*sum_exp - beta*sum_x_exp)/(sum_exp*sum_exp);
} 

} // end of namespace

