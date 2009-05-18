#ifndef INFONODE_HH_
#define INFONODE_HH_

#include "DataInOut/ParamHandling/ParamNode.hh"
#include <complex>
#include <map>

namespace CoupledField
{

/** Declare the global pointer info created in cfs.cc */
class InfoNode;
extern InfoNode* info;

class InfoNode : public ParamNode
{
public:
  /** Default constructor for StdVector.
   * @param filename if given one can do a easy (repeated) ToFile() w/p parameters.
   * @param preamble is written at the beginning of the file */
  InfoNode(const std::string& filename = "", const std::string& preamble = "");

  /** ParamNode handles the destructor stuff */
  ~InfoNode() {};

  //typedef enum LogLevel { NOT_SET, DISABLED, BRIEF, DETAILED, VERBOSE, DEBUG, ALL};
  enum LogLevel { NOT_SET, DISABLED, BRIEF, DETAILED, VERBOSE};

  enum Cardinality { USE_EXISTING, REPLACE, APPEND, UNIQUE};

  enum ValueType { BOOL, INTEGER, REAL, COMPLEX, STRING, MATRIX, VECTOR, UNKNOWN };

  /** Gives back The InfoNode, depending on the cardinality this is a old one.
   * @param name a valid xml name
   * @param comment if given a child element is created, one can doe by hand
   * @param name either a valid XML name (no spaces, brackets, ...) or "", then caption is used
   *             Multiple levels might be spilt via slashes like in xpath.
   * @param card applies only on the last level if there are multiple via slashes. */
  InfoNode* Get(const std::string& name, const std::string& comment, Cardinality card = USE_EXISTING)
  {
    return Get(name, comment, "", "", card);
  }

  /** Works like the Has(parent, child, value) in ParamNode but handles the '/' delemiter! */
  InfoNode* Get(const std::string& name, const std::string& child,
                const std::string& value, Cardinality card = USE_EXISTING)
  {
    return Get(name, "", child, value, card);    
  }

  
  InfoNode* Get(const std::string& name, Cardinality card = USE_EXISTING)
  {
    return Get(name, "", "", "", card);
  }


  /** The ParamNode::SetString() shall not be used, This SetValue() methods also set
   * valueType (explicit and implicit)
   * todo: implement CDATA stuff :)
   * @param string automatically checks if it has to become CDATA if conating &lt; , ... */
  void SetValue(const std::string& string);

  /** @param if c_str is null  an exception is thrown conditionally */
  void SetValue(const char* c_str, bool throw_exception = true) 
  {
    if(c_str == NULL)
    {
      if(throw_exception) 
        EXCEPTION("invalid attempt to set NULL string als value")
      else SetValue("");  
    }
    else    
      SetValue(std::string(c_str));
  }

  void SetValue(bool param);

  void SetValue(int param);

  void SetValue(unsigned int param);

  void SetValue(double param);

  void SetValue(const std::complex<double>& param);

  /** copies the node with all children. Overwrites the current name */
  void SetValue(ParamNode* node);
  
  /** Creates a sub-node with the content */
  void SetComment(const std::string& string);

  /** Prints this as xml element to the stream. Builds a tree. Shall no be directly
   * called for an attribute
   * @param depth number of ident steps, will be interpreted as one space */
  void ToXML(std::ostream& os, int depth = 0) const;

  /** Write the current state to the file. Uses existing filename and preabmle if given
   * @param filename if not set the must be one given in the constructor,
   *        if given overwrites an existing value
   * @param preamble is written at the beginning of the file. Logic as for filename */
  void ToFile(const std::string& filename = "", const std::string& preamble = "");

  /** This string constant shall mark the header of a part - e.g. defaults and assumptions */
  const static std::string HEADER;

  /** This string constant shall mark the logging part which contains "dynamic" data, e.g. current iteration */
  const static std::string PROCESS;
  /** this string constann shall mark summary information, e.g. total number of iterations */
  const static std::string SUMMARY;

  const static std::string COMMENT;
  const static std::string WARNING;
  const static std::string ERROR;

private:

  InfoNode* Get(const std::string& name, const std::string& comment,
                const std::string& child, const std::string& value,
                Cardinality card = USE_EXISTING);

  /** Init COMMENT, ... */
  static void InitStatic();

  /** limits the siginificant digits of real and complex */
  std::string GetFormatedValue() const;

  /** Makes a valid XML name, e.g. removes spaces. Used to create
   * automatically an ParamNode::name_ value out of caption_
   * @param in might contain spaces, e.g. "Number of iterations"
   * @return for this example 'NumberOfIterations' */
  std::string ToValidLabel(const std::string& in) const;

  /** due to a strange compiler error, try to remove */
  bool IsAttribute() const { return attribute_; }

  /** Is the string value good for a attribute (length, quotation mark).
   * To decide if comment and value are XML-elements or attributes */
  bool IsGoodAttribute(const std::string& string) const;

  /** this static counter is for the number of ToFile() calls by filename. */
  static std::map<const std::string, unsigned int> writeCounter_;
  
  /** This is the label of the entry.
   * If a ParamNode::name_ is not explicit give, it correlates via ToValidLabel()
   * with capation_ */
  std::string caption_;

  /** We have to know about out parent for parameter inheritance!
   *  NULL only for ROOT */
  InfoNode* parent_;

  /** if NOT_SET the log level inherits from the first set parent */
  LogLevel log_;

  Cardinality cardinality_;

  ValueType type_;

  // e.g. the root element can be given a filename for (intermedia) saves
  std::string filename_;
  // is written as preamble to files by WriteToFile();
  std::string preamble_;


};


} // end of namespace


#endif /*INFONODE_HH_*/
