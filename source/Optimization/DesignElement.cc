#include "General/exception.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/Condition.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"

using namespace CoupledField;

// the static enum
Enum DesignElement::filter;
Enum DesignElement::type;
Enum DesignElement::valueSpecifier;
Enum DesignElement::access;
Enum DesignElement::detail;

DesignElement::DesignElement()
{
    filter_         = false;
    design          = 0.0;
    weight          = 1.0;
    cost_gradient   = -1.0;
}


double DesignElement::GetLowerBound() const
{
  switch(type_)
  {
    // this is the default in literature - playing with it
    // might be worse - especially when doing penalization!
    case DENSITY: return 1.0e-3;
    case POLARIZATION: return -1.0;
    default: throw Exception("type not implemented");
  }
} 

double DesignElement::GetUpperBound() const
{
  switch(type_)
  {
    case DENSITY: return 1.0;
    case POLARIZATION: return 1.0;
    default: throw Exception("type not implemented");
  }
} 

void DesignElement::SetConstraintGradient(const Condition* condition, double value) 
{ 
  this->constraintGradient[condition->GetIndex()] = value; 
}

double DesignElement::GetConstraintGradient(const Condition* condition) 
{ 
  return this->constraintGradient[condition->GetIndex()]; 
}


double DesignElement::GetValue(ResultDescription* rd)
{
  // check for special result
  if(rd->value == OBJECTIVE || (rd->value == COST_GRADIENT && rd->detail != NONE))
  {
    switch(rd->solutionType)
    {
      case OPT_RESULT_1: return specialResult[0];
      case OPT_RESULT_2: return specialResult[1];
      case OPT_RESULT_3: return specialResult[2];
      default: throw Exception("solution type not handled");
    }
  }

  return GetValue(rd->value, rd->access); 
}

double DesignElement::GetValue(ValueSpecifier vs, Access access)
{
  if(access == PLAIN || !filter_) 
     return GetValue(vs);
  else 
     return GetFilteredValue(vs);
}

double DesignElement::GetValue(ValueSpecifier sp)
{
   switch(sp)
   {
       case DESIGN:               return design;
       case DESIGN_COST_GRADIENT: return design * cost_gradient;
       case COST_GRADIENT:        return cost_gradient;
       case WEIGHT:               return weight;
       case CONSTRAINT_GRADIENT:  
            throw Exception("for constraint gradient we need an index!");
            
       default: throw Exception("value specifier not implemented"); 
   }
}

double DesignElement::GetFilteredValue(ValueSpecifier sp)
{
    // We filter over this element and the neighbours. 
     
    double factor_sum = 0.0;
    double weight_sum = 0.0;

    for(UInt i = 0; i < neighbourhood.GetSize(); i++) 
    {
        const NeighbourElement& ne = neighbourhood[i];
        double factor = ne.weight * ne.neighbour->GetValue(sp);
        factor_sum += factor;
        weight_sum += ne.weight;
    }
    
    // add our part :)
    factor_sum += GetValue(sp);
    weight_sum += 1.0;
    
    // in the 99-lines paper this is inverse of this_value
    double result = factor_sum / weight_sum;

    // std::cout << "return " << result << " instead of " << GetValue(sp) << std::endl;
    return result;      
}


double DesignElement::GetDesign(Access access)
{
    if(access == PLAIN || !filter_) 
       return design;
    else 
       return GetFilteredValue(DESIGN);
}

double DesignElement::GetObjectiveGradient(Access access)
{
    if(access == PLAIN || !filter_) 
       return GetValue(COST_GRADIENT);
    else
       // see the filter definition for gradients in the
       // 99-lines paper why we use here "design TIMES cost_gradient" 
       return GetFilteredValue(DESIGN_COST_GRADIENT);
}

void DesignElement::Dump()
{
   double weight_sum = weight;
   double distance_avg = 0.0;
   for(unsigned int i = 0; i < neighbourhood.GetSize(); i++)
   {
      weight_sum   += neighbourhood[i].weight;
      distance_avg += neighbourhood[i].distance;
   }
   distance_avg /= neighbourhood.GetSize();

   std::cout << "\nelement: " << elem->elemNum << " location " << location_.ToString() 
             << " weight sum " << weight_sum << " this weight " << weight <<" distance avg " << distance_avg << std::endl;
   for(unsigned int i = 0; i < neighbourhood.GetSize(); i++)
       std::cout << "  n[" << i << "]: elem " << neighbourhood[i].neighbour->elem->elemNum << " location " << neighbourhood[i].neighbour->location_.ToString() << " dist=" << neighbourhood[i].distance << " w=" << neighbourhood[i].weight << std::endl; 
}


void DesignElement::SetEnums()
{
   filter.SetName("DesignElement::Filter");
   filter.Add(RADIUS, "radius");
   filter.Add(NEIGHBOURS, "neighbours");
   
   type.SetName("DesignElement::Type");
   type.Add(DEFAULT, "default");
   type.Add(DENSITY, "density");
   type.Add(POLARIZATION, "polarization");

   access.SetName("DesignElement::Access");
   access.Add(PLAIN, "plain");
   access.Add(SMART, "smart");
   
   valueSpecifier.SetName("DesignElement::ValueSpecifier");
   valueSpecifier.Add(DESIGN, "design");
   valueSpecifier.Add(DESIGN_COST_GRADIENT, "designTimesCostGradient");
   valueSpecifier.Add(COST_GRADIENT, "costGradient");
   valueSpecifier.Add(CONSTRAINT_GRADIENT, "constraintGradient");
   valueSpecifier.Add(WEIGHT, "weight");
   valueSpecifier.Add(OBJECTIVE, "objective");
   
   detail.SetName("DesignElement::Detail");
   detail.Add(NONE, "none");
   detail.Add(UKU, "uKu");
   detail.Add(UKP, "uKp");
   detail.Add(PKU, "pKu");
   detail.Add(PKP, "pKp");      
}
  

void DesignElement::InitFilter(StdVector<DesignElement>& data, ParamNode* pn)
{
    Type   type   = (Type) DesignElement::type.Parse(pn->Get("design")->AsString());
    Filter filter = (Filter) DesignElement::filter.Parse(pn->Get("type")->AsString());
    double value  = pn->Get("value")->AsDouble();
    if(filter == NEIGHBOURS && ((int) value) == 0) 
      throw Exception("Filtering with 0 neighbours not possible");
 
    // set the positions
    Matrix<double>  coords;
    Grid* grid = domain->GetGrid();
    
    for(UInt i = 0; i < data.GetSize(); i++)
    {
      // do only the filtering for our design type - density might be different from polarization
      if(data[i].type_ != type) continue;
      
      StdVector<UInt>& connect = data[i].elem->connect;
      // do not use updated coordinates up to now!!
      grid->GetElemNodesCoord(coords, connect, false );
       
      // a barycenter is simply the average of all coordinates -> location_ gets the Point out
      BaseFE::CalcBarycenter(coords, data[i].location_);
    }
    
    // look for the neighbours. This is n^2 but Manfred said, that there
    // will be a smarter structure soon in CFS to do it faster.
    // for the same reason GridCFS<DIM>::GetElemsNextToNodes() is not used.
    for(UInt element = 0; element < data.GetSize(); element++)
    {
        DesignElement* de = &data[element];  

        de->filter_ = true;
        Point& loc = de->location_;
        
        
        for(UInt n = 0; n < data.GetSize(); n++)
        {
            // an element is not a neighbour of itself!
            if(de->elem->elemNum == data[n].elem->elemNum) continue;
              
            double distance = dist(loc, data[n].location_);
            if(filter == RADIUS && distance < value)
            {
                // value is here a double radius
                // this is the implementation from Bendsoe/ Sigmund
                NeighbourElement ne;
                ne.neighbour = &data[n];
                ne.weight    = value - distance;
                ne.distance  = distance;
                de->neighbourhood.Push_back(ne); // let the default copy constructor work!
                // std::cout << "element " << data[element].elem->elemNum << " " 
                //          << data[element].location_.ToString() << " has neighbour "
                //          << data[n].elem->elemNum << " " <<  data[n].location_.ToString()
                //          << " distance= " << distance << " weight=" << ne.weight << std::endl; 
            }
            // with neighbours we don't know the weight a priori as we don't know the furthest distance
            if(filter == NEIGHBOURS)
            {
                unsigned int size = (int) value; // value is here a int size
                // if neighbourhood smaller size we are in
                if(de->neighbourhood.GetSize() < size) 
                {
                    NeighbourElement ne;
                    ne.neighbour = &data[n];
                    ne.distance  = distance;
                    de->neighbourhood.Push_back(ne);
                } 
                else 
                {
                   // check if we can replace any of the neighbours. search biggest
                   NeighbourElement* biggest = &de->neighbourhood[0];
                   // first search smallest
                   for(unsigned int e = 1; e < de->neighbourhood.GetSize(); e++) 
                   {
                      if(de->neighbourhood[e].distance > biggest->distance) biggest = &de->neighbourhood[e]; 
                   }
                   // compare our value with smallest and replace - dont't set the weight yet
                   if(distance < biggest->distance) 
                   {
                      // replace the attributes
                      biggest->neighbour = &data[n];
                      biggest->distance  = distance;
                   }
                }
            } // end filter == NEIGHBOURS
        } // here the element is calculated
        
        // in the NEIGHBOURS case the weights were not set. We cannot use the furthest distance
        if(filter == NEIGHBOURS)
        {
           // search largest distance
           double max_distance = 0.0;
           for(unsigned int e = 0; e < de->neighbourhood.GetSize(); e++)
              max_distance = std::max(de->neighbourhood[e].distance, max_distance);  

           // to what we did with radius but the 0-line is 1.2 * max_distance
           for(unsigned int e = 0; e < de->neighbourhood.GetSize(); e++)
           { 
              de->neighbourhood[e].weight = (max_distance * 1.2) - de->neighbourhood[e].distance;
              if(de->neighbourhood[e].weight < 0.0)
                throw Exception("we have a negative weight!");
           }  
        }
        
        // now normalize the weights. The weight of this element is by definition 1.0
        double sum = 1.0;  
        for(unsigned int e = 0; e < de->neighbourhood.GetSize(); e++) 
           sum += de->neighbourhood[e].weight;
           
        // now normalize all
        de->weight /= sum;
        for(unsigned int e = 0; e < de->neighbourhood.GetSize(); e++) 
           de->neighbourhood[e].weight /= sum;  
                 
        // de->Dump();
    }
    
    // for direct debug output: determin neighbourhood statistics
    UInt max_neighbours = 0;
    UInt min_neighbours = 0xffffff;
    UInt avg_neighbours = 0;
    for(UInt i = 0; i < data.GetSize(); i++) {
       UInt s = data[i].neighbourhood.GetSize();
       max_neighbours = std::max(max_neighbours, s);
       min_neighbours = std::min(min_neighbours, s);
       avg_neighbours += s;

    }
    avg_neighbours /= data.GetSize();
    
    std::cout << "Neighbourhoods: min = " << min_neighbours << " max = " << max_neighbours << " avg = " << avg_neighbours << std::endl;  
}


ResultDescription::ResultDescription()
{
  // does not set all!!
  access = DesignElement::PLAIN;
  value  = DesignElement::DESIGN;
  design = DesignElement::DEFAULT;
}

ResultDescription::ResultDescription(ParamNode* pn)
{
   // killme: do this nice when there is no more environment.hh
   SolutionType st; 
   String2Enum(pn->Get("id")->AsString(), st);
   solutionType = st;

   design = DesignElement::DEFAULT;
   if(pn->Has("design"))
     design = (DesignElement::Type) DesignElement::type.Parse(pn->Get("design"));
 
   access = (DesignElement::Access) DesignElement::access.Parse(pn->Get("access"));
   
   value = (DesignElement::ValueSpecifier) DesignElement::valueSpecifier.Parse(pn->Get("value"));
   
   detail = (DesignElement::Detail) DesignElement::detail.Parse(pn->Get("detail")); 
}

