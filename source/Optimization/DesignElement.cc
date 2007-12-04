#include "General/exception.hh"
#include "Optimization/DesignElement.hh"
#include "Optimization/Condition.hh"
#include "Domain/domain.hh"
#include "Domain/grid.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/Logging/cfslog.hh"

using namespace CoupledField;

DECLARE_LOG(filterLog)
DEFINE_LOG(filterLog, "filter")

// the static enum
Enum<DesignElement::Filter>         DesignElement::filter;
Enum<DesignElement::Type>           DesignElement::type;
Enum<DesignElement::ValueSpecifier> DesignElement::valueSpecifier;
Enum<DesignElement::Access>         DesignElement::access;
Enum<DesignElement::Detail>         DesignElement::detail;

/** The default constructor is only for the StdVector */
DesignElement::DesignElement()
{
  // dont' use this empty constructor!
}

DesignElement::DesignElement(ParamNode* pn)
{
  filter_         = false;
  design          = 0.0;
  weight          = 1.0;
  cost_gradient   = -1.0;

  // it is a little slow to perform this code for every DesignElement but the
  // implementations are rater fast and it should be not measurable in the end
  type_ = type.Parse(pn->Get("name"));  

  upper_ = 1.0;
  // eventually overwrite
  pn->Get("upper", upper_, false);
  
  lower_ = type_ == POLARIZATION ? -1.0 : 0.001;
  pn->Get("lower", lower_, false);
}


void DesignElement::SetConstraintGradient(const Condition* condition, double value) 
{ 
  this->constraintGradient[condition->GetIndex()] = value; 
}

double DesignElement::GetConstraintGradient(const Condition* condition) const
{ 
  return this->constraintGradient[condition->GetIndex()]; 
}


double DesignElement::GetValue(ResultDescription* rd) const
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

double DesignElement::GetValue(ValueSpecifier vs, Access access) const
{
  double val = 0.0;
  if(access == PLAIN || !filter_)
    val = GetValue(vs);
  else 
    val = GetFilteredValue(vs, false);
  
  LOG_DBG2(filterLog) << elem->elemNum << " GetValue(" << valueSpecifier.ToString(vs) << ", " 
                      << access << " -> " << val;
  return val;
}

double DesignElement::GetObjectiveGradient(Access access) const
{
  double val = 0;
  if(access == PLAIN || !filter_) 
    val = GetValue(COST_GRADIENT);
  else
    // see the filter definition for gradients in the
    // 99-lines paper why we use here "design TIMES cost_gradient" 
    val = GetFilteredValue(COST_GRADIENT, true);
  
  LOG_DBG2(filterLog) << elem->elemNum << " GetObjectiveGradient() -> " << val;
  return val;
}

double DesignElement::GetDesign(Access access) const
{
    if(access == PLAIN || !filter_) 
       return design;
    else 
       return GetFilteredValue(DESIGN, false);
}



double DesignElement::GetValue(ValueSpecifier sp) const
{
   switch(sp)
   {
       case DESIGN:               return design;
       case DESIGN_COST_GRADIENT: return design * cost_gradient;
       case COST_GRADIENT:        return cost_gradient;
       case WEIGHT:               return weight;
       case NUM_NEIGHBOURS:       return neighbourhood.GetSize();
       case CONSTRAINT_GRADIENT:  
            throw Exception("for constraint gradient we need an index!");
            
       default: throw Exception(valueSpecifier.ToString(sp) + " not stored in design element"); 
   }
}

double DesignElement::GetFilteredValue(ValueSpecifier sp, bool design_weighted) const
{
    // We filter over this element and the neighbours. 
     
    double factor_sum = 0.0;
    double weight_sum = 0.0;

    for(UInt i = 0; i < neighbourhood.GetSize(); i++) 
    {
        const NeighbourElement& ne = neighbourhood[i];
        double des_weight = design_weighted ?  ne.neighbour->design : 1.0;
        double factor = ne.weight * ne.neighbour->GetValue(sp) * des_weight;
        factor_sum += factor;
        weight_sum += ne.weight;
        LOG_DBG3(filterLog) << this->elem->elemNum << ": neighbour " << ne.neighbour->elem->elemNum 
                            << " weight = " << ne.weight
                            << " des_weight " << des_weight
                            << " value = " << ne.neighbour->GetValue(sp);
    }
    // add our part :)
    double des_weight = design_weighted ? this->design : 1.0;

    LOG_DBG2(filterLog) << this->elem->elemNum << ": factor_sum = " << factor_sum << " + " << GetValue(sp)
                        << " weight_sum = " << weight_sum << " + " << weight << " -> " 
                        << ((factor_sum+(GetValue(sp)*des_weight))/((weight_sum+weight)*des_weight)) << " instead of " << GetValue(sp);
    
    factor_sum += GetValue(sp) * des_weight;
    weight_sum += this->weight;
    
    // in the 99-lines paper this is inverse of this_value
    double result = factor_sum / (weight_sum * des_weight);

    return result;      
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
   filter.Add(VOLUME_RADIUS, "volumeRadius");
   filter.Add(MAX_EDGE, "maxEdge");
   
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
   valueSpecifier.Add(NUM_NEIGHBOURS, "neighbours");
   
   detail.SetName("DesignElement::Detail");
   detail.Add(NONE, "none");
   detail.Add(UKU, "uKu");
   detail.Add(UKP, "uKp");
   detail.Add(PKU, "pKu");
   detail.Add(PKP, "pKp");
   detail.Add(SYMMETRY, "symmetry");
}
  

void DesignElement::InitFilter(StdVector<DesignElement>& data, ParamNode* pn)
{
    Type   type   = DesignElement::type.Parse(pn->Get("design"));
    Filter filter = DesignElement::filter.Parse(pn->Get("type"));
    double value  = pn->Get("value")->AsDouble();

    double avg_radius = 0;
    double avg_neighbours = 0;
    
    // set the positions
    Matrix<double>  coords;
    Grid* grid = domain->GetGrid();
    double dim = (double) grid->GetDim();
    
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
        
        double radius = value;
        if(filter == VOLUME_RADIUS)
        {
          // TODO really check for axis symmetry off
          double volume = de->elem->ptElem->CalcVolume(coords, false);
          // The radius is <value> times square/cube edge length where the
          // square/cube has the volume of the element
          radius = value * std::pow(volume, 1.0/dim);
          LOG_DBG3(filterLog) << "from volume " << volume << " to radius " << radius;
        }
        if(filter == MAX_EDGE)
        {
          double min;
          de->elem->ptElem->GetMaxMinEdgeLength(coords, radius, min);
          LOG_DBG3(filterLog) << "Element " << de->elem->elemNum << " edge max=" << radius << " min=" << min; 
        }
        
        for(UInt n = 0; n < data.GetSize(); n++)
        {
            // an element is not a neighbour of itself!
            if(de->elem->elemNum == data[n].elem->elemNum) continue;
              
            double distance = Point::dist(loc, data[n].location_);
            if(distance < radius)
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
        } // here the element is calculated
        
        // now normalize the weights. The weight of this element is by definition 1.0
        double sum = 1.0;  
        for(unsigned int e = 0; e < de->neighbourhood.GetSize(); e++) 
           sum += de->neighbourhood[e].weight;
           
        // now normalize all
        de->weight /= sum;
        for(unsigned int e = 0; e < de->neighbourhood.GetSize(); e++) 
           de->neighbourhood[e].weight /= sum;  
                 
        avg_radius += radius;
        avg_neighbours += de->neighbourhood.GetSize();
    }
    
    // for direct debug output: determin neighbourhood statistics
    std::cout << "Filter: avg radius: " << (avg_radius / data.GetSize()) 
              << " avg neighbourhood: " << (avg_neighbours / data.GetSize()) << std::endl;
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
     design = DesignElement::type.Parse(pn->Get("design"));
 
   access = DesignElement::access.Parse(pn->Get("access"));
   
   value = DesignElement::valueSpecifier.Parse(pn->Get("value"));
   
   detail = DesignElement::detail.Parse(pn->Get("detail")); 
}

std::string ResultDescription::ToString()
{
   // killme: do this nice when there is no more environment.hh
   std::string string;
   Enum2String((SolutionType) solutionType, string);
   return string;
}
