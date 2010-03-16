#ifndef ParamNode2_HH_
#define ParamNode2_HH_

#include "DataInOut/ParamHandling/ParamNode.hh"
#include <complex>
#include <boost/any.hpp>
#include <map>

namespace CoupledField
{

class InfoNode;
/** Definitions of pointers, using boost::shared_ptr */
typedef boost::shared_ptr<InfoNode> PtrParamNode;
typedef StdVector<boost::shared_ptr<InfoNode> > InfoNodeList;

/** global parameter instance */
extern PtrParamNode info;

/** Add class description */
class InfoNode : public ParamNode {
public:

  /** Default constructor */
  InfoNode( ActionType defaultAction = INSERT, NodeType type = UNDEF );

  /** Destructor */
  virtual ~InfoNode();
};

} // end of namespace


#endif /*ParamNode2_HH_*/
