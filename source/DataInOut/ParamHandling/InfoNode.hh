#ifndef INFONODE_HH_
#define INFONODE_HH_

#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField
{

class InfoNode : public ParamNode
{
public:
  
  
  /** Default constructor for StdVector */
  InfoNode();

  /** ParamNode handles the destructor stuff */
  ~InfoNode() {};
  
  typedef enum LogLevel { NOT_SET, DISABLED, BRIEF, DETAILED, VERBOSE, DEBUG, ALL};
  
  /** For main elements (PDE, Optimization, ...) there is a seperation into
   * three parts. XSLT can make use of it! */ 
  typedef enum Part { STARTUP, MAIN, FINALIZE, NO_PART };
  
  typedef enum Cardinality { USE_EXISTING, REPLACE, APPEND, UNIQUE};
  
  typedef enum ValueType { BOOL, INTEGER, REAL, COMPLEX, STRING, MATRIX, VECTOR, UNKNOWN };
  
  /** Gives back The InfoNode, depending on the cardinality this is a old one.
   * One has to give either name or caption or both. If name is not given, it
   * is created out of caption via ToValidLabel()
   * @param name either a valid XML name (no spaces, brackets, ...) or "", then caption is used
   * @param caption the label for the output ('Number of iterations'). If "" then name is used */
  InfoNode* Create(const std::string& name, const std::string& caption, 
                   Cardinality card = USE_EXISTING, Part part = NO_PART
                   LogLevel log = NOT_SET); 
  
  /** The ParamNode::SetString() shall not be used, This SetValue() methods also set
   * valueType (explicit and implicit)
   * todo: implement CDATA stuff :)
   * @param string automatically checks if it has to become CDATA if conating &lt; , ... */
  void SetValue(const std::string& string)
  {
    value_ = string;
    valueType_ = STRING;
  }
  
private:
  
  
  
  
  /** Makes a valid XML name, e.g. removes spaces. Used to create
   * automatically an ParamNode::name_ value out of caption_ 
   * @param in might contain spaces, e.g. "Number of iterations"
   * @return for this example 'NumberOfIterations' */
  std::string ToValidLabel(const std::string& in) const;
  
  /** This is the label of the entry.
   * If a ParamNode::name_ is not explicit give, it correlates via ToValidLabel()
   * with capation_ */
  std::string caption_;
  
  /** We have to know about out parent for parameter inheritance!
   *  NULL only for ROOT */
  InfoNode* parent_;

  /** if NOT_SET the log level inherits from the first set parent */
  LogLevel log_;

  /** which part/ sequence are we? */
  Part part_;
  
  Cardinality cardinality_;
  
  ValueType type_;
};


} // end of namespace


#endif /*INFONODE_HH_*/
