#include "DataInOut/ParamHandling/ParamNode.hh"
namespace CoupledField 
{
 

  /** This is the global instance of the global pointer */
  PtrParamNode info;
  
  InfoNode::InfoNode( ActionType defaultAction, NodeType type ) 
    : ParamNode(defaultAction, type )
  {
  }
  
  InfoNode::~InfoNode()
  {
  }
}
