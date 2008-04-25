#include "DataInOut/ParamHandling/InfoNode.cc"

using namespace CoupledField;

InfoNode::InfoNode()
{
  parent_      = NULL;
  log_         = NOT_SET;
  part_        = NO_PART;
  cardinality_ = USE_EXISTING; // default
  type_        = UNKNOWN; // to be set with SetValue()
}


InfoNode* InfoNode::Create(const std::string& name, const std::string& caption, 
                 Cardinality card, Part part, LogLevel log)
{
  if(name == "" && caption == "") EXCEPTION("need at least name or caption");
  // use label to search for existing nodes
  std::string label = name == "" ? ToValidLabel(caption) : name;
  
  // decide on cardinality
  // search for an existing InfoNode
  InfoNode* node = NULL;
  if(Has(label)) node = Get(label);
  
  switch(card)
  {
  case USE_EXISTING:
    if(node != NULL) return old;
    break;
    
  case REPLACE:
    if(node != NULL) 
    {
      children_.remove(old);
      delete node;
      node = NULL;
    }
    break;
    
  case APPEND:
    if(node != NULL && node->cardinality_ == UNIQUE)
      EXCEPTION("You wanted to add an InfoNode '" << label 
                << "' but there is already an UNIQUE one");
    break; // we don't care if there is an old one
    
  case UNIQUE:
    if(node != NULL) 
      EXCEPTION("InfoNode '" << label << "' already exists but you requested unique");
    break;
  }
  
  // create and add a new one if old == NULL
  if(node == NULL)
  {
    node = new InfoNode();
    children_.push_back(node);
  }
  
  node->parent_ = this;
  
  node->caption_ = caption == "" ? name : caption;
  node->name_    = label; // see above!
  
  node->cardinality_ = card;
  node->part_        = part;
  node->log_         = log;
  
  return node;
}
