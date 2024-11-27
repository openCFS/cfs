#include "Domain/Domain.hh"
#include "Domain/Mesh/Grid.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "Utils/mathParser/mathParser.hh"
#include "Optimization/Transform.hh"
#include "Optimization/Design/DesignSpace.hh"
#include "Optimization/Design/DesignElement.hh"

using namespace CoupledField;

Enum<Transform::Type> Transform::type;

DEFINE_LOG(transform, "transform")

Transform::Transform(PtrParamNode pn, DesignSpace* space)
{
  space_ = space;

  type_ = type.Parse(pn->GetName());

  index = -1; // to be set properly from outside

  assert(type_ == ROTATION); // nothing else implemented yet

  // the angle can be described as a formula
  MathParser* parser = domain->GetMathParser();
  unsigned int handle = parser->GetNewHandle();
  parser->SetExpr(handle, pn->Get("angle")->As<std::string>());
  angle_ = parser->Eval(handle);  
  parser->ReleaseHandle(handle);
  
  center_str_     = pn->Get("center")->As<std::string>();

  // during constructor call we have no multiple excitations yet set up.
  excitation_str = pn->Get("excitation")->As<std::string>();

  // we don't case excitation_str to the index as there is no easy way to check it here w/o any context

  StdVector<Elem*> elems;
  domain->GetGrid()->GetElemsByName(elems, center_str_);

  if(elems.GetSize() != 1)
    throw Exception("transform: the named element '" + center_str_ + "' for center of rotation shall be only 1 node!");

  // most probably the barycenters are not calculated yet
  // domain->GetGrid()->GetElemShapeMap(elems[0], false)->CalcBarycenter(center_);
  // elems[0]->ptrShapeMap->CalcBarycenter(center_);
  elems[0]->GetElemShapeMap(domain->GetGrid(), false)->CalcBarycenter(center_);

  double tmp;
  // domain->GetGrid()->GetElemShapeMap(elems[0], false)->GetMaxMinEdgeLength(h_, tmp); // h_ is the largest edge size
  // elems[0]->ptrShapeMap->GetMaxMinEdgeLength(h_, tmp);
  elems[0]->GetElemShapeMap(domain->GetGrid(), false)->GetMaxMinEdgeLength(h_, tmp);

  LOG_DBG(transform) << "T a='" << pn->Get("angle")->As<std::string>() << " -> " << angle_ << " c=" << elems[0]->elemNum;

}

DesignElement* Transform::FindSource(const DesignElement* de) const
{
  // we assume that the barycenters are calculated!
  assert(domain->GetGrid()->regionData[space_->regions[0][0].regionId].barycenters == true);

  assert(domain->GetGrid()->GetDim() == 2); // nothing else implemented yet

  // calculate source coordinates
  double x = de->elem->extended->barycenter[0];
  double y = de->elem->extended->barycenter[1];

  double cx = center_[0];
  double cy = center_[1];

  double x2 = std::sin(angle_) * (y - cy) + std::cos(angle_) * (x - cx) + cx;
  double y2 = std::cos(angle_) * (y - cy) - std::sin(angle_) * (x - cx) + cx;

  LOG_DBG2(transform) << "FS de=" << de->ToString(true) << " -> (" << x2 << ", " << y2 << ")";

  // relay make it faster by mapping the design to a defined oriented regular grid!

  // find the element behind x2, y2 by assuming a regular grid and traversing the coordinate.
  // in case we have a circular design region the order (first x, then y) or vice versa is important

  // first try x then y
  DesignElement* end = NULL;
  DesignElement* tmp = SearchDesignSpace(de, x2 > x ? VicinityElement::X_P : VicinityElement::X_N, x2);
  // go only in second direction if the first search was good
  if(tmp != NULL)
    end = SearchDesignSpace(tmp, y2 > y ? VicinityElement::Y_P : VicinityElement::Y_N, y2);
  else {
    // search the other way
    tmp= SearchDesignSpace(de, y2 > y ? VicinityElement::Y_P : VicinityElement::Y_N, y2);
    if(tmp != NULL)
      end = SearchDesignSpace(tmp, x2 > x ? VicinityElement::X_P : VicinityElement::X_N, x2);
  }


  LOG_DBG2(transform) << "FS org=" << de->ToString(true) << " tmp=" << DesignElement::ToString(tmp,true) << " end=" << DesignElement::ToString(end,true);

  return end;
}


DesignElement* Transform::SearchDesignSpace(const DesignElement* start, VicinityElement::Neighbour dir, double val) const
{
  assert(start != NULL);

  // the index within a Point
  int axis = VicinityElement::ToMainAxis(dir);

  const DesignElement* cand = start; // the start might be what we look for (e.g. if we rotate the center around the center)
  double closest = std::abs(val - start->elem->extended->barycenter[axis]); // let's hope negative coordinates make no problems :(

  bool stop = false;
  while(!(stop || closest < 1e-15))
  {
    DesignElement* test = cand->vicinity->GetNeighbour(dir);

    if(test == NULL) {
      stop = true; // let's check before returning if the last element of the domain was what we were looking for
      continue;
    }

    double dist = std::abs(test->elem->extended->barycenter[axis] - val);

    LOG_DBG3(transform) << "SDS start=" << start->ToString(true) << " dir=" << dir << " << val=" << val << " cand=" << cand->ToString(true) << " test=" << test->ToString(true) << " c=" << closest << " d=" << dist;

    // we assume that we get linearly closer and then grow again. The growing again indicates the before was closest
    if(dist <= closest) { // still not growing again
      closest = dist;
      cand = test;
    }
    else { // we went one step too far and are growing
      stop = true;
      continue;
    }
  }

  if(closest <= 1.1 * h_)
    return const_cast<DesignElement*>(cand);
  else
    return NULL;
}

std::string Transform::ToString(int level)
{
  std::stringstream ss;
  ss.precision(3);
  if(level < 0)
    ss << "angle=" << angle_ << " center=" << center_.ToString() << " ex=" << excitation_str;
  if(level == 0)
    ss << angle_;
  else
    ss << excitation_str << "_" << angle_;
  return ss.str();
}

void Transform::ToInfo(PtrParamNode in)
{
  PtrParamNode in_ = in->Get(type.ToString(type_), ParamNode::APPEND);
  in_->Get("angle")->SetValue(angle_);
  in_->Get("center")->SetValue(center_str_);
  in_->Get("coord")->SetValue(center_.ToString());
  in_->Get("excitation")->SetValue(excitation_str);
}




