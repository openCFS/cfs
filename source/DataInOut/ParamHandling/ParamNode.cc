#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/ColoredConsole.hh"
#include "General/defs.hh"
#include "MatVec/Matrix.hh"
#include "MatVec/Vector.hh"
#include "Utils/Timer.hh"
#include "Utils/mathParser/mathParser.hh"
#include "DataInOut/ProgramOptions.hh"
#include "Domain/Domain.hh"
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/device/file_descriptor.hpp>
#include <fstream>
#include <string>
#include <utility>

using namespace std;
using namespace boost;

namespace CoupledField
{

/** set the global names for fields */
const string ParamNode::HEADER = "header";
const string ParamNode::PROCESS = "process";
const string ParamNode::SUMMARY = "summary";

const string ParamNode::WARNING = "warning";
const string ParamNode::FAIL = "error";

/** This is our global pointer of the root ParamNode holding the XML file.
 *  Filed in cfs.cc. The corresponding
 * "extern PtrParamNode param;" is in ParamNode.hh */


ParamNode::ParamNode(ActionType defaultAction, NodeType type) :
  precision_(5), 
  name_("DD"), 
  type_(type),
  defaultAction_(defaultAction),
  lastresultidx_(-1)

{ }

ParamNode::~ParamNode()
{
  if(rootNode) {
    delete rootNode;
    rootNode = NULL;
  }

}

PtrParamNode ParamNode::GenerateWriteNode(const string& root, const string& filename, ActionType defaultAction, bool lazy_write, bool add_counter)
{
  assert(defaultAction != DEFAULT);

  PtrParamNode node = PtrParamNode(new ParamNode(defaultAction, ParamNode::ELEMENT));
  node->SetName("cfsInfo");
  node->rootNode = new RootNode();
  node->rootNode->filename = filename;
  #ifndef NDEBUG
    node->rootNode->lazy_write = lazy_write;
  #else
    node->rootNode->lazy_write = false;
  #endif
  node->rootNode->add_counters = add_counter;
  return node;
}

/*****************static PtrParamNode GenerateWriteNode(const string& root, const string& filename, bool lazy_write, bool add_counter);*******************************************************
 * S E T    M E T H O D S
 *************************************************************************/

void ParamNode::SetValue(const boost::any& value, bool cerr_warning)
{
  //std::cout << "ParamNode::SetValue(const boost::any&)" << std::endl;
  this->value_ = value;

  // note that we may not check for valid strings here, as e.g. < or > are valid inputs for attribute values for expressions
  if(this->name_ == WARNING && cerr_warning)
  {
    if(progOpts == nullptr || !progOpts->IsQuiet())
      std::cerr  << std::endl;
    std::cerr  << fg_red << "WARNING: " << boost::any_cast<std::string>(value_)<< fg_reset << std::endl;
  }
}
void ParamNode::SetValue(const char* value)
{
  std::string str(value);
  SetValue(str);
}

void ParamNode::SetValue(const double value, const int precision)
{
  this->precision_ = precision;
  SetValue(value);
}

void ParamNode::SetValue(PtrParamNode node, bool overwrite_name, bool cerr_warning)
{
  //std::cout << "ParamNode::SetValue(PtrParamNode) name" << node->GetName() << std::endl;
  if(overwrite_name)
    SetName(node->GetName());

  assert(name_ != "");

  // set the value
  SetValue(node->value_, cerr_warning);

  ParamNodeList& children = node->GetChildren();
  // reserve own children size
  children_.Resize(children.GetSize());

  // run recursively through all children
  for(UInt i = 0; i < children.GetSize(); i++)
  {
    // add new element
    PtrParamNode other = children[i];
    PtrParamNode new_node = SetNewChild(other->GetName(), i); // faster for large arrays
    new_node->SetValue(other, false, cerr_warning);
  }
}

/** Set ParamaNodes as child nodes. Works recursively */
void ParamNode::SetValue(ParamNodeList nodes)
{
  for(auto add : nodes) {
    PtrParamNode pn = Get(add->GetName(), APPEND);
    //std::cout << "ParamNode::SetValue(ParamNodeList) call SetValue" << std::endl;
    pn->SetValue(add, false, false);
  }
}

void ParamNode::SetComment(const std::string& comment)
{
  // check if we already have the commend
  if(HasByVal("comment", comment))
    return;

  // add new child node with type COMMENT
  PtrParamNode newChild(new ParamNode(defaultAction_, COMMENT));
  newChild->SetName("comment");
  newChild->SetValue(comment);
  newChild->parent_ = weak_from_this();
  newChild->defaultAction_ = defaultAction_;
  children_.Push_back(newChild);
}

void ParamNode::SetWarning(const std::string& msg, bool append)
{
  PtrParamNode pn = Get(ParamNode::WARNING, append ? APPEND : DEFAULT);
  pn->SetType(ELEMENT); // attributes would not be able to append
  pn->SetValue(msg);
}

void ParamNode::AddChildNode(PtrParamNode child)
{
  child->parent_ = weak_from_this();
  if (child->defaultAction_ == DEFAULT)
  {
    child->defaultAction_ = defaultAction_;
  }
  children_.Push_back(child);
}

PtrParamNode ParamNode::SetNewChild(const std::string& name, unsigned int index)
{

  PtrParamNode node(new ParamNode());
  node->SetName(name);
  node->parent_ = weak_from_this();
  node->defaultAction_ = defaultAction_;
  children_[index] = node;
  return node;
}

PtrParamNode ParamNode::ReplaceChild(PtrParamNode node, unsigned int index)
{
  NodeType type = children_[index]->type_;
  node->parent_ = weak_from_this();
  node->SetType(type);
  children_[index] = node;
  return node;
}

void ParamNode::ClearChildren()
{
  children_.Clear();
}

void ParamNode::ClearChildren(const string& child_name)
{
  ParamNodeList tmp;
  tmp.Reserve(children_.GetSize()); // avoid expensive Push_back
  for(unsigned int i = 0; i < children_.GetSize(); i++)
    if(children_[i]->GetName() != child_name)
      tmp.Push_back(children_[i]);

  if(tmp.GetSize() == children_.GetSize())
    return; // nothing matched, nothing to do

  children_.Clear(true); // keep capacity

  // copy back if there is some important stuff we wanted to keep
  for(unsigned int i = 0; i < tmp.GetSize(); i++)
    children_.Push_back(tmp[i]);
}


/************************************************************************
 * N O D E   A C C E S S     M E T H O D S
 ************************************************************************/
PtrParamNode ParamNode::GetChild()
{
  // check that only one child is present
  if (children_.GetSize() != 1)
  {
    EXCEPTION( "Element '" << name_ << "' must only have 1 child");
  }
  return children_[0];
}

PtrParamNode ParamNode::Get(const string& name_raw, ActionType action)
{
  PtrParamNode result;
  // make sure we have valid element and attribute names.
  string myName = ToValidLabel(name_raw);

  if (action == DEFAULT)
    action = defaultAction_;

  // check if name contains special "/" signs to
  if (ContainsTokens(myName))
  {
    //result = TokenizedHasAndGet(name, std::string(""), false ); // has some overhead
    StdVector<string> tokens = SplitIntoTokens(myName);
    PtrParamNode actNode = shared_from_this();
    for (unsigned int i = 0; i < tokens.GetSize(); i++)
    {
      myName = tokens[i];

      // here we have to be careful: if we have the actiontype APPEND,
      // we want to perform this just on the last element and re-use
      // the previous ones
      if (action == APPEND && i != tokens.GetSize() - 1)
      {
        action = INSERT;
      }
      // check, if we still have a result
      if (actNode)
      {
        result = actNode->Get(myName, action);
        actNode = result;
      }
    }
    return result;
    // set name to last element
    //myName = tokens[tokens.GetSize()-1];
    // setting the 'ptr' is done below!

  }
  else
  {
    const int chsize(children_.GetSize());

    // check if the last result is close to the currently requested node
    // This is a hard-coded optimization coming from F. Wein
    if (lastresultidx_ > 0 && lastresultidx_ < chsize)
    {
      if (children_[lastresultidx_]->name_ == myName)
        result = children_[lastresultidx_];
      else if (lastresultidx_ + 1 < chsize)
      {
        if (children_[lastresultidx_ + 1]->name_ == myName)
          result = children_[++lastresultidx_];
      }
    }

    for (int i = 0; i < chsize && result == NULL; i++)
    {
      if (children_[i]->name_ == myName)
      {
        lastresultidx_ = i;
        result = children_[i];
      }
    }
  }

  // perform action, in case no node was found
  if (result == NULL || action == APPEND) {
    // depending on ActionType:
    switch (action) {
    case DEFAULT:
      EXCEPTION("Some action has to be performed");
      break;
    case EX: {
      std::string matName;
      GetValue("name", matName, ParamNode::EX);
      EXCEPTION("None of the " << children_.GetSize() << " childs of element '" << this->name_ << "' with name '" << matName << "' has a child '" << myName << "' ");
      break;
    }
    case PASS:
      break;
    case INSERT:
    case APPEND:
      // create new node containing "default"
      PtrParamNode newChild(new ParamNode(defaultAction_));
      newChild->SetName(myName);
      // ATTENTION: Do NOT set an empty string as value to this node, as
      // std::string("").empty() != boost::any(st::string("")).empty()
      newChild->parent_ = weak_from_this();
      newChild->defaultAction_ = defaultAction_;
      children_.Push_back(newChild);
      result = newChild;
      break;
    }
  }

  return result;
}

PtrParamNode ParamNode::GetRoot() {
  PtrParamNode parent_locked = parent_.lock();
  if( parent_locked ) {
    return parent_locked->GetRoot();
  } else {
    return shared_from_this();
  }
}


string ParamNode::GetLocation() const
{
  string loc;
  const ParamNode* pn = this;
  while(pn != nullptr)
  {
    loc = "/" + pn->name_ + loc;
    pn = pn->parent_.lock().get();
  }
  return loc;
}


template<typename TYPE>
PtrParamNode ParamNode::GetByVal(const string& parent_raw, const string& child_raw,
    const TYPE& value, ActionType action)
{
  string parent = ToValidLabel(parent_raw);
  string child  = ToValidLabel(child_raw);

  // Strategy depending on action:
  // EX / PASS / INSERT: Search for node. If not found
  //   EX / PASS : throw exception / pass
  //   INSERT: INSERT content with current one
  // APPEND: Directly insert node, regardless of existing ones
  if (action == DEFAULT)
    action = defaultAction_;

  bool insertNew = false;
  if (action != APPEND)
  {
    StdVector<PtrParamNode> result = GetListByVal(parent, child, value);

    if (result.GetSize() == 1)
      return result[0];

    if (result.GetSize() > 1)
    {
      if (action == PASS || action == INSERT)
        return result[0];
      else
      EXCEPTION("element '" << parent << "' has more than one (" << result.GetSize()
          << ") childs '" << child << "' with value '" << value << "'");
    }

    // result size = 0!
    if (action == INSERT)
      insertNew = true;
    if (action == EX)
      EXCEPTION("element '" << parent << "' has no child '" << child
          << "' with value '" << value << "'");
  }

  // Insert new node with given value
  // action is either APPEND or INSERT, i.e. the following unconstrained
  // Get()-calls will always create a child elements
  if(insertNew || action == APPEND)
  {
    //changed 2nd parameter to APPEND
    PtrParamNode ret = this->Get(parent, ParamNode::APPEND);
    ret->Get(child, ParamNode::APPEND)->SetValue(value);
    return ret;
  }

  return PtrParamNode();
}

PtrParamNode ParamNode::GetByVal(const std::string& parent,
    const std::string& child, const char* value, const ActionType action)
{
  return this->GetByVal(parent, child, std::string(value), action);
}


PtrParamNode ParamNode::GetByVal(const string& parent_raw, const string& child1_raw,  const string& value1,
                                                           const string& child2_raw,  const string& value2, ActionType action)
{
  string parent = ToValidLabel(parent_raw);
  string child1  = ToValidLabel(child1_raw);
  string child2  = ToValidLabel(child2_raw);

  if (action == DEFAULT)
    action = defaultAction_;

  ParamNodeList l = GetListByVal(parent_raw, child1, value1);

  if(action == APPEND || (l.IsEmpty() && action == INSERT))
  {
    //changed 2nd parameter to APPEND
    PtrParamNode ret = this->Get(parent, ParamNode::APPEND);
    ret->Get(child1, ParamNode::APPEND)->SetValue(value1);
    ret->Get(child2, ParamNode::APPEND)->SetValue(value2);
    return ret;
  }
  if(l.IsEmpty())
    EXCEPTION("parent " << parent_raw << " has no child " << child1 << " with value " << value1);

  for(unsigned int i = 0; i < l.GetSize(); i++)
  {
    // assert(l[i]->HasByVal(child1, value1));
    if(l[i]->HasByVal(child2, value2))
      return l[i]; // stupid, as we ignore more matches
  }

  if(action == INSERT)
  {
    PtrParamNode ret = this->Get(parent, ParamNode::APPEND);
    ret->Get(child1, ParamNode::APPEND)->SetValue(value1);
    ret->Get(child2, ParamNode::APPEND)->SetValue(value2);
    return ret;
  }
  else
    EXCEPTION("parent " << parent_raw << " has  child " << child1 << " with value " << value1 <<
            " but no child " << child2 << " with value " << value2);
}

ParamNodeList ParamNode::GetList(const string& name)
{
  string vl = ToValidLabel(name);
  StdVector<PtrParamNode> result;
  result.Reserve(children_.GetSize());

  for (unsigned int i = 0; i < children_.GetSize(); ++i)
    if (children_[i]->name_ == vl)
      result.Push_back(children_[i]);

  return result; // copy-constructor magic stuff!
}


ParamNodeList ParamNode::GetListByChild(const ParamNodeList& base, const std::string&  name)
{
  string vl = ToValidLabel(name);
  StdVector<PtrParamNode> result;
  for(unsigned int b = 0; b < base.GetSize(); b++)
  {
    if(base[b]->name_ == vl)
      result.Push_back(base[b]);
   }
  return result;
}


ParamNodeList ParamNode::GetListByGrandChild(const ParamNodeList& base, const std::string&  name)
{
  string vl = ToValidLabel(name);
  StdVector<PtrParamNode> result;
  for(unsigned int b = 0; b < base.GetSize(); b++)
  {
    PtrParamNode pn = base[b];
    for(unsigned int c = 0; c < pn->children_.GetSize(); c++)
    {
      if(pn->children_[c]->name_ == vl)
        result.Push_back(pn->children_[c]);
    }
  }
  return result;
}

template<typename TYPE>
ParamNodeList ParamNode::GetListByVal(const string& parent_raw,
    const string& child_raw, const TYPE& value)
{
  StdVector<PtrParamNode> result;

  string parent = ToValidLabel(parent_raw);
  string child  = ToValidLabel(child_raw);

  for (unsigned int p = 0, chsize = children_.GetSize(); p < chsize; ++p)
  {
    PtrParamNode ch = children_[p];
    // do we have parent name?
    if (ch->GetName() == parent)
    {
      // children_[i] has zero, one or more child matching child elements,
      // we loop by hand to handle the "more" case
      const unsigned int chchsize(ch->children_.GetSize());
      for (unsigned int c = 0; c < chchsize; ++c)
      {
        const PtrParamNode chch = ch->children_[c];
        if (chch->GetName() == child && chch->As<TYPE> () == value)
        {
          // parent, child and value matches. Take it!
          result.Push_back(ch);
        }
      }
    }
  }

  return result; // copy-constructor magic stuff!
}

ParamNodeList ParamNode::GetListByVal(const string& parent,
    const string& child, const char * value)
{
  return GetListByVal(parent, child, std::string(value));
}

/************************************************************************
 * D A T A    G E T   M E T H O D S
 ************************************************************************/

StdVector<std::string>& ParamNode::GetFastBulkBlock()
{
  if(value_.empty())
  {
    StdVector<std::string> tmp;
    SetValue(tmp);
    return boost::any_cast<StdVector<std::string>&>(value_); // As() is const!
  }
  else
  {
    if(value_.type() == typeid(StdVector<std::string>))
      return boost::any_cast<StdVector<std::string>&>(value_);

    EXCEPTION("cannot provide fast bulk block, value already set to '" << ToString(0) << "'");
  }
}

template<typename TYPE>
TYPE ParamNode::As() const
{
  return boost::any_cast<TYPE>(value_);

}
#define AS_INTEGRAL( TYPE   )\
  template<>\
  TYPE ParamNode::As<TYPE>() const {\
    if( value_.empty()) return TYPE();\
    if( value_.type() == typeid(TYPE) ) {\
      return boost::any_cast<TYPE>(value_);\
    } else {\
      try {\
        if( value_.type() == typeid(int) ) {\
          int val = boost::any_cast<int>(value_);\
          return boost::lexical_cast<TYPE>(val);\
        }\
        if( value_.type() == typeid(unsigned int) ) {\
          unsigned int val = boost::any_cast<unsigned int>(value_);\
          return boost::lexical_cast<TYPE>(val);\
        }\
        if( value_.type() == typeid(double) ) {\
          double val = boost::any_cast<double>(value_);\
          return boost::lexical_cast<TYPE>(val);\
        }\
        if( value_.type() == typeid(std::string) ) {\
            std::string val = boost::any_cast<std::string>(value_);\
            return boost::lexical_cast<TYPE>(val);\
          } else {\
            EXCEPTION("cannot cast value in XML node '" << name_\
                              << "' to desired type.");\
          }\
      } catch (boost::bad_lexical_cast & ) {\
        EXCEPTION("cannot cast value in XML node '" << name_\
                  << "' to desired type.");\
      }\
    }\
    TYPE dummy = TYPE();\
     return dummy;\
  }

AS_INTEGRAL(Integer)
AS_INTEGRAL(UInt)
AS_INTEGRAL(Double)
AS_INTEGRAL(std::string)

// special implementation for bool
template<>
bool ParamNode::As<bool>() const
{
  if(value_.type() == typeid(bool))
    return boost::any_cast<bool>(value_);

  if(value_.type() == typeid(std::string))
  {
    std::string str = boost::any_cast<std::string>(value_);
    if(str == "yes" || str == "true" || str == "on" || str == "enable" || str == "1")
      return true;
    if(str == "no" || str == "false" || str == "off" || str == "disable" || str == "0" || str == "none")
      return false;

   EXCEPTION("Cannot convert node '" << name_ << "' with value '" << str << "' to boolean");
  }
  EXCEPTION("Cannot convert node '" << name_ << "' to boolean, it's neither string nor bool");
}

boost::shared_ptr<Timer> ParamNode::AsTimer()
{
  if(value_.empty())
    value_ = boost::shared_ptr<Timer>(new Timer());

  if(value_.type() != typeid(boost::shared_ptr<Timer>))
    EXCEPTION("param node " << name_ << " cannot be returned as timer");

  return boost::any_cast<boost::shared_ptr<Timer> >(value_);
}

template<typename TYPE>
const TYPE& ParamNode::AsConst() const
{
  // first, try to directly use any_cast<Integer>
  try
  {
    // first check, if type is same as returned value
    return boost::any_cast<const TYPE&>(value_);

  } catch (boost::bad_any_cast &)
  {
    EXCEPTION("Cannot determine the data type of xml node '" << name_ << "'.");

  }
}

template<typename TYPE>
TYPE ParamNode::MathParse() const
{
  if(value_.empty())
    return TYPE();

  // obtain handle
  MathParser* parser = domain->GetMathParser();
  unsigned int handle = parser->GetNewHandle(false);

  std::string expr = "";
  if(value_.type() == typeid(std::string))
    expr = boost::any_cast<std::string>(value_);
  else
    EXCEPTION("XML node '" << name_ << "' can not be parsed, as it is no string type.");

  // Set expression and evaluate
  parser->SetExpr(handle, expr);
  TYPE ret = TYPE();
  ret = static_cast<TYPE>(parser->Eval(handle));
  
  // release handle
  parser->ReleaseHandle(handle);
  
  return ret;
}

template<typename TYPE>
void ParamNode::GetValue(const std::string& name, TYPE& ret, ActionType action)
{

  if (action == DEFAULT)
    action = defaultAction_;

  if (action != APPEND)
  {
    // First, check if node exists
    PtrParamNode node = this->Get(name, PASS);

    if (node)
    {
      // => node exists: convert value and return
      ret = node->As<TYPE> ();
    }
    else
    {
      switch (action)
      {
      case PASS:
        return;
        break;
      case EX:
        EXCEPTION("element '" << name_ << "' has no child '" << name << "'")
        ;
        break;
      case INSERT:
      default:
        break;
      } // switch
    } // if node exists
  } // if APPEND

  // At this point we know, that we have to create a new entry, either
  // setting it uniquely (action == INSERT) or appending it (action==APPEND)
  PtrParamNode myNode = this->Get(name, action);
  myNode->SetValue(ret);
}

/************************************************************************
 *  Q U E R Y     M E T H O D S
 ************************************************************************/
bool ParamNode::Has(const string& name) const
{
  // check in a fast way if we have tokens for trivial xpath
  if(ContainsTokens(name)) // check for path, if so, we cannot compare for name
    return TokenizedHasAndGet(name, string(""), false) == NULL ? false : true ;

  return GetIndex(name) >= 0;
}

int ParamNode::GetIndex(const string& name) const
{
  assert(!ContainsTokens(name));

  for(unsigned int i = 0, chsize = children_.GetSize(); i < chsize; i++)
  {
    if(children_[i]->name_ == name)
      return i;
  }
  return -1;
}


template<typename TYPE>
bool ParamNode::HasByVal(const string& parent, const string& child, const TYPE& value) const
{
  if(ContainsTokens(parent))
    EXCEPTION("HasByVal(parent, child, value) does not allow for tokenized search strings!");

  // see GetList() for comments
  for (unsigned int p = 0, chsize = children_.GetSize(); p < chsize; ++p)
  {
    PtrParamNode ch = children_[p];
    if (ch->name_ == parent)
    {
      const unsigned int chchsize(ch->children_.GetSize());
      for (unsigned int c = 0; c < chchsize; ++c)
      {
        PtrParamNode chch(ch->children_[c]);
        if (chch->GetName() == child && chch->As<TYPE> () == value)
          return true;
      }
    }
  }

  // nothing found
  return false;
}

bool ParamNode::HasByVal(const string& parent, const string& child,
    const char* value) const
{
  return HasByVal(parent, child, std::string(value));
}

template<typename TYPE>
bool ParamNode::HasByVal(const string& name, const TYPE& value) const
{
  if (ContainsTokens(name))
  {
    return TokenizedHasAndGet(name, value, true) == NULL ? false : true ;
  }
  else
  {
    for (unsigned int i = 0, chsize = children_.GetSize(); i < chsize; i++)
      if (children_[i]->name_ == name && children_[i]->As<TYPE> () == value)
        return true;

    return false;
  }
}

bool ParamNode::HasByVal(const string& name, const char* value) const
{
  return HasByVal(name, std::string(name));
}

bool ParamNode::HasChildren() const
{
  return (bool) children_.GetSize();
}

unsigned int ParamNode::Count(const string& name) const
{
  unsigned int count = 0;

  for (unsigned int i = 0, chsize = children_.GetSize(); i < chsize; ++i)
    if (children_[i]->name_ == name)
      ++count;

  return count;
}
/************************************************************************
 * M I S C    M E T H O D S
 ************************************************************************/

template<typename TYPE>
PtrParamNode ParamNode::TokenizedHasAndGet(const string& name,
    const TYPE& value, bool restrictToVal) const
{
  StdVector<string> tokens = SplitIntoTokens(name);
  PtrParamNode ptr = boost::const_pointer_cast<ParamNode>(shared_from_this());
  for (unsigned int i = 0; i < tokens.GetSize(); i++)
  {

    // consider the value only for the last token!
    if (i == tokens.GetSize() - 1)
    {
      // only look for specified value, if value is not empty
      // otherwise, we have already the searched for value
      if (restrictToVal)
      {
        if (!ptr->HasByVal(tokens[i], value))
          return PtrParamNode();
      }
      else
      {
        if (!ptr->Has(tokens[i]))
          return PtrParamNode();
      }
    }
    else
    {
      if (!ptr->Has(tokens[i]))
        return PtrParamNode();
    }
    ptr = ptr->Get(tokens[i]);
  }

  return ptr;
}

std::string ParamNode::ToString(int depth) const
{
  std::string tmp;
  ToString(tmp, depth);
  return tmp;
}


void ParamNode::ToString(std::string& ret, int depth) const
{
  // This method is currently very hard coded.
  // In the future we rely on a general formatter class,
  // which can convert the special types
  // default case for not value (but children)
  if (value_.type() == typeid(std::string))
  {
    ret = boost::any_cast<std::string>(value_);
    // https://www.w3.org/TR/2006/REC-xml11-20060816/
    std::replace(ret.begin(), ret.end(), '"', '\''); // " -> '
    boost::replace_all(ret, "&", "&amp;");
    boost::replace_all(ret, "<", "&lt;");
    boost::replace_all(ret, ">", "&gt;");
    return;
  }
  if (value_.type() == typeid(double))
  {
    double val = boost::any_cast<double>(value_);
    std::stringstream stream;
    stream.precision(precision_);
    stream << val;
    ret = stream.str();
    return;
  }
  if (value_.type() == typeid(Complex))
  {
    Complex val = boost::any_cast<Complex >(value_);
    std::stringstream stream;
    stream.precision(precision_);
    stream << val;
    ret = stream.str();
    return;
  }
  if (value_.type() == typeid(UInt))
  {
    UInt val = boost::any_cast<UInt>(value_);
    ret = boost::lexical_cast<std::string>(val);
    return;
  }
  if (value_.type() == typeid(Integer))
  {
    Integer val = boost::any_cast<Integer>(value_);
    ret = boost::lexical_cast<std::string>(val);
    return;
  }
  if (value_.type() == typeid(bool))
  {
    bool val = boost::any_cast<bool>(value_);
    ret = val ? "yes" : "no";
    return;
  }
  // special conversion types for vector
  if (value_.type() == typeid(Vector<Double> ))
  {
    Vector<Double> vec = boost::any_cast<Vector<Double> >(value_);
    ret = vec.ToString();
    return;
  }
  if (value_.type() == typeid(Vector<Complex> ))
  {
    Vector<Complex> vec = boost::any_cast<Vector<Complex> >(value_);
    ret = vec.ToString();
    return;
  }
  // special conversion types for matrix
  if (value_.type() == typeid(Matrix<Double> ))
  {
    Matrix<Double> mat = boost::any_cast<Matrix<Double> >(value_);
    ret = mat.ToXMLFormat(name_, depth);
    return;
  }
  if (value_.type() == typeid(Matrix<Complex> ))
  {
    Matrix<Complex> mat = boost::any_cast<Matrix<Complex> >(value_);
    ret = mat.ToXMLFormat(name_, depth);
    return;
  }

  // special conversion types for vector
  if (value_.type() == typeid(Vector<Double>*))
  {
    Vector<Double> & vec = *(boost::any_cast<Vector<Double> *>(value_));
    ret = vec.ToString();
    return;
  }
  if (value_.type() == typeid(Vector<Complex>*))
  {
    Vector<Complex> & vec = *(boost::any_cast<Vector<Complex>*>(value_));
    ret = vec.ToString();
    return;
  }
  // special conversion types for matrix
  if (value_.type() == typeid(Matrix<Double>*))
  {
    Matrix<Double> & mat = *(boost::any_cast<Matrix<Double>*>(value_));
    ret = mat.ToXMLFormat(name_, depth);
    return;
  }
  if (value_.type() == typeid(Matrix<Complex>*))
  {
    Matrix<Complex> & mat = *(boost::any_cast<Matrix<Complex>*>(value_));
    ret = mat.ToXMLFormat(name_, depth);
    return;
  }

  // special conversion for timer
  if(value_.type() == typeid(boost::shared_ptr<Timer>))
  {
    boost::shared_ptr<Timer> timer = boost::any_cast<boost::shared_ptr<Timer> >(value_);
    // Note, that we are a SELF_XML type
    ret = timer->ToXMLFormat(name_);
    return;
  }

  if(value_.type() == typeid(StdVector<std::string>))
  {
    ret = "error in fast bulk block writing"; // this should not be printed
    return;
  }
}

void ParamNode::ToXML(std::ostream& os, int depth, bool adjust_element_type)
{
  // normally this does not need to be called!
  if(adjust_element_type)
    AdjustElementType(); // is recursively!

  assert(name_ != "");
  // note, that this is an recursive method!

  std::string strValue;
  ToString(strValue, depth);

  if (type_ == ATTRIBUTE)
  {
    // makes only sense in an recursive call
    os << " " << name_ << "=\"" << strValue << "\"";
    return;
  }

  // we are NO Attribute

  // The sorting-thing should be performed by an external instance
  // in the future
  //Sort(); // (HEADER before) PROCESS before SUMMARY

  // if we start a new element or are part of an element is same/same
  os << std::endl << string(std::max(depth, 0), ' ');

  if (type_ != COMMENT && type_ != SELF_XML)
  {
    os << "<" << name_;
  }

  // go through all children an check if we close with "../>" or </name_>
  // note, that the own XML-element value can also be a attribute
  bool only_attributes = true;

  const UInt chsize(children_.GetSize());
  for (UInt i = 0; i < chsize; i++)
  {
    PtrParamNode actNode = children_[i];
    if (actNode->type_ != ATTRIBUTE)
      only_attributes = false;
    break;
  }

  // first loop, do the attributes
  for (UInt i = 0; i < chsize; i++)
  {
    PtrParamNode actNode = children_[i];
    if (actNode->type_ == ATTRIBUTE)
    {
      children_[i]->ToXML(os);
    }
    else
    {
      only_attributes = false; // no, there is a non-attribute element
    }
  }

  // SELF_XML types have their own element opening and closing with this name as element name
  if (type_ == SELF_XML)
  {
    os << strValue;
    return; // done, everything is printed!
  }

  // the special fast bulk block mode
  if(type_ == BULK)
  {
    os << ">" << std::endl;
    StdVector<std::string>& block = GetFastBulkBlock();
    assert(!block.IsEmpty());
    for(UInt i = 0, n = block.GetSize(); i < n; i++)
    {
      os << string(std::max(depth + 2, 0), ' ');
      os << block[i] << std::endl;
    }
    os << string(std::max(depth, 0), ' ');
    os << "</" << name_ << ">";
    return;
  }



  if (only_attributes)
  {
    if (type_ == COMMENT)
    {
      os << "<!-- " << strValue << "-->";
    }
    else
    {
      // check, if value is set
      if (!value_.empty())
      {
        assert(true); // this should not be possible and be checked in AdjustType()
        os << "> " << strValue << "</" << name_ << ">" << std::endl;
      }
      else
      {
        os << "/>"; // "quick close"
      }
    }
  }
  else
  {
    os << ">";

    // second loop, child elements
    for (unsigned int i = 0; i < chsize; i++)
    {
      // attributes are already done
      if (children_[i]->type_ == ATTRIBUTE)
        continue;
      children_[i]->ToXML(os, depth + 2);
    }
    bool endl_written = false;

    // own element
    //      if(type_ == COMMENT) {
    //        os << "<!-- " << ToString() << "-->" << std::endl;
    //      } else
    //      if (type_ == ELEMENT ) {
    //        // check, if value is set
    //        if (!value_.empty() ) {
    //          os << string(depth + 2, ' ') << strValue;
    //        }
    //        os << std::endl;
    //        endl_written = true;
    //      }

    // do we close in the same line?
    if (chsize != 0 && !endl_written)
      os << std::endl;
    os << string(std::max(depth, 0), ' ');
    if (type_ != COMMENT)
    {
      os  << "</" << name_ << ">";
    }
  }
}

void ParamNode::ToFile(const std::string& filename, bool force)
{
  assert(rootNode != NULL || !filename.empty());

  // silent mode, meant for dummy
  if(rootNode != NULL && rootNode->filename.empty())
    return;
  
  std::string myFileName = !filename.empty() ? filename : rootNode->filename;

  if(rootNode != NULL)
  {
    // only really write the file if at least a certain amount of time has passed since last write
    // or if forced. Write always in the debug mode
    if(rootNode->lazy_write && !force && rootNode->write_timer->IsRunning() && rootNode->write_timer->GetWallTime() < 2.0)
    {
      rootNode->reject_counter++;
      return;
    }

    // we will write in some microseconds
    rootNode->write_timer->ResetStart();
    rootNode->write_counter++;
  }

  bool compress = ends_with(myFileName, ".gz");

  // a little complicated by transparent for encoding
  ofstream file(myFileName.c_str(), ios_base::out);
  boost::iostreams::filtering_ostream out;
  if(compress)
    out.push(boost::iostreams::gzip_compressor()); // gzip might not compress as good as bzip2 but is ways faster!
  out.push(file);

  // write preamble
  out << "<?xml version=\"1.0\"?>" << std::endl;

  // store how often we are written -> if the number is too high one should cancel some ToFile() calls
  //    if(writeCounter_.count(filename_) == 0) writeCounter_[filename_] = 0;
  //    Get("writeCounter")->SetValue(++writeCounter_[filename_]);
  // just for info.xml put write counter to info.xml
  if(rootNode != NULL && rootNode->add_counters) {
    Get("infoWriteCounter")->SetValue(rootNode->write_counter);
    Get("infoRejectCounter")->SetValue(rootNode->reject_counter);
  }

  // First we determine the type of nodes
  AdjustElementType();

  // Then we print the tree
  ToXML(out);

  out << std::endl;
}

void ParamNode::Dump(int level) const
{
  for (int i = 0; i < level; i++)
    cout << "   ";

  std::string strValue;
  ToString(strValue, level);
  // Format element name
  switch (type_)
  {
  case UNDEF:
    cout << name_ << ": " << strValue << std::endl;
    break;
  case ELEMENT:
    cout << "<" << name_ << "> " << strValue << std::endl;
    break;
  case SELF_XML:
    cout << strValue << std::endl;
    break;
  case ATTRIBUTE:
    cout << "@" << name_ << ": " << strValue << std::endl;
    break;
  case COMMENT:
    cout << "<!-- " << strValue << " -->" << std::endl;
    break;
  case BULK:
    {
      const StdVector<std::string>& block = boost::any_cast<const StdVector<std::string>&>(value_);
      for(UInt i = 0, n = block.GetSize(); i < n; i++)
        cout << " bulk: " << block[i] << std::endl;
      break;
    }
  }
  for (unsigned int i = 0, chsize = children_.GetSize(); i < chsize; i++)
    children_[i]->Dump(level + 1);
}


StdVector<std::pair<std::string, std::string>> ParamNode::ToStringList(int max_level) const
{
  StdVector<std::pair<std::string, std::string>> res;
  ToStringList(res, max_level, 0); // let the recursive game begin
  return res;
}


void ParamNode::ToStringList(StdVector<std::pair<std::string, std::string> >& list, int max_level, int level) const
{
  if(level > max_level)
    return;

  // level 0 name is not printed, level 1 is printed and only from level 2 we have a parent name chain
  string parentname;
  for(int add = level; add >= 2; add--) // e.g. we have level 2
  {
    PtrParamNode base = parent_.lock(); // this is then level 1
    for(int search = 1; search < add-2; search++)
    {
      base = base->parent_.lock();
      assert(base);
    }
    parentname += base->name_ + ":";
  }


  string val = As<string>();
  if(val != "")
    list.Push_back(std::make_pair(parentname + name_, val));


  for(auto child : children_)
      child->ToStringList(list, max_level, level +1);
}

StdVector<string> ParamNode::SplitIntoTokens(const string& input) const
{
  StdVector<string> out;
  char_separator<char> sep("/");
  tokenizer<char_separator<char> > tokens(input, sep);

  for (tokenizer<char_separator<char> >::iterator it = tokens.begin(); it != tokens.end(); ++it)
    out.Push_back(*it);
  return out;
}

void ParamNode::AdjustElementType()
{

  // loop over all children and call recursivley method
  // to determine type of nodes
  ParamNodeList::iterator it = children_.Begin();
  for (; it != children_.End(); it++)
  {
    (*it)->AdjustElementType();
  }

  // Adjust element type
  if (type_ == UNDEF)
  {

    // if children are present -> ELEMENT
    if (children_.GetSize())
    {
      type_ = ELEMENT;
    }
    else
    {
      // if no children are present
      //   value_ not set -> ELEMENT
      //   value_ set     -> ATTRIBUTE
      if (value_.empty())
      {
        type_ = ELEMENT;
      }
      else
      {
        type_ = ATTRIBUTE;
      }
    }
  }
  // Check if node is ELEMENT, has children and a value
  // -> not possible in an XML tree
  if (type_ == ELEMENT && children_.GetSize() && !value_.empty() && value_.type() != typeid(StdVector<std::string>))
  {
    string value;
    ToString(value, 0);
    if(value != "")
    {
      WARN("Node '" << name_ << "' has children AND a non-empty value '" << value << "'. This is not possible in an xml tree!");
      assert(false); // find the stuff. Maybe attribute and element with the same name!
    }
  }

  // Check if node is ATTRIBUTE and has children
  // -> not possible in XML tree
  if (type_ == ATTRIBUTE && children_.GetSize())
  {
    WARN("Node '" << name_ << "' is a ATTRIBUTE and has child nodes. "
        << "This is not possible in a xml tree!");
    assert(false); // find the stuff
  }

  // there are special elements (Matrix, Timer(, Vector if someone implements it)
  // which do their own xml formatting by ToXMLFormat(), they are of type SELF_XML
  if(   value_.type() == typeid(Matrix<Double>*) || value_.type() == typeid(Matrix<Complex>*)
     || value_.type() == typeid(Matrix<Double>)  || value_.type() == typeid(Matrix<Complex>)
     || value_.type() == typeid(boost::shared_ptr<Timer>))
  {
    type_ = SELF_XML;
  }

  if(value_.type() == typeid(StdVector<std::string>))
    type_ = BULK;

  // ToDO: What is currently missing is the check for valid labels etc.
  // Maybe Fabian is willing to assist here ...
}

//  bool ParamNode::IsGoodAttribute(const string& attr) const
//  {
//    //if(attr == DATA_IS_IN_MATRIX_VALUE) return false;
//    //if(attr == DATA_IS_IN_TIMER_VALUE)  return false;
//
//    if(attr.size() > 30) return false;
//    if(attr.find('\"') != string::npos) return false;
//    return true;
//  }

inline std::string ParamNode::ToValidLabel(std::string out)
{
  boost::trim(out);

  // remove spaces and don't make upper case yet! :(
  boost::erase_all(out, " ");
  boost::erase_all(out, "\"");
  boost::erase_all(out, "(");
  boost::erase_all(out, ")");
  // don't touch the slash '/', it is a 'xpath' element
  return out;
}

ParamNode::RootNode::RootNode()
{
  write_timer = new Timer();
}


ParamNode::RootNode::~RootNode()
{
  if(write_timer != NULL) {
    delete write_timer;
    write_timer = NULL;
  }
}

// ========================================================================
//  Instantiate public methods
// ========================================================================
// A) INSTANTIATION OF ALL TYPES
// -----------------------------
// Just ParamNode::As<> and ParamNode::GetData<> can be used for
// general types, as we
// we rely on operator==() and operator<<() for the other types
// (the integral types can be casted using the boost::lexical_cast<>
// statement
#define INSTANTIATE_METHOD_AS(TYPE)\
  template\
  TYPE ParamNode::As<TYPE>() const;
INSTANTIATE_METHOD_AS(Vector<Double>)
INSTANTIATE_METHOD_AS(Vector<Complex>)
INSTANTIATE_METHOD_AS(Matrix<Double>)
INSTANTIATE_METHOD_AS(Matrix<Complex>)
INSTANTIATE_METHOD_AS(Vector<Double>*)
INSTANTIATE_METHOD_AS(Vector<Complex>*)
INSTANTIATE_METHOD_AS(Matrix<Double>*)
INSTANTIATE_METHOD_AS(Matrix<Complex>*)
INSTANTIATE_METHOD_AS(boost::shared_ptr<Timer>)

#define INSTANTIATE_METHOD_MATH_PARSE(TYPE)\
  template\
  TYPE ParamNode::MathParse<TYPE>() const;
INSTANTIATE_METHOD_MATH_PARSE(Double)
INSTANTIATE_METHOD_MATH_PARSE(UInt)
INSTANTIATE_METHOD_MATH_PARSE(Integer)


#define INSTANTIATE_METHOD_GETVALUE(TYPE)\
  template void ParamNode::GetValue<TYPE>\
  (const std::string& name, TYPE& ret, const ActionType);

INSTANTIATE_METHOD_GETVALUE(Integer)
INSTANTIATE_METHOD_GETVALUE(Double)
INSTANTIATE_METHOD_GETVALUE(UInt)
INSTANTIATE_METHOD_GETVALUE(bool)
INSTANTIATE_METHOD_GETVALUE(std::string)
INSTANTIATE_METHOD_GETVALUE(Vector<Double>)
INSTANTIATE_METHOD_GETVALUE(Vector<Complex>)
INSTANTIATE_METHOD_GETVALUE(Matrix<Double>)
INSTANTIATE_METHOD_GETVALUE(Matrix<Complex>)
INSTANTIATE_METHOD_GETVALUE(Vector<Double>*)
INSTANTIATE_METHOD_GETVALUE(Vector<Complex>*)
INSTANTIATE_METHOD_GETVALUE(Matrix<Double>*)
INSTANTIATE_METHOD_GETVALUE(Matrix<Complex>*)
INSTANTIATE_METHOD_GETVALUE(boost::shared_ptr<Timer>)
// B)Now just the integral types
// -----------------------------
#define INSTANTIATE_METHODS_INT(TYPE)\
    template\
     PtrParamNode ParamNode\
    ::GetByVal<TYPE>(const string& parent,\
               const string& child,\
               const TYPE& value,\
               const ActionType action );\
    template\
    ParamNodeList ParamNode::GetListByVal<TYPE>\
    (const string& parent, const string& child, const TYPE& value);\
  template  bool ParamNode::HasByVal<TYPE>\
  ( const std::string& name, const std::string& child,const TYPE& value ) const;\
  template  bool ParamNode::HasByVal<TYPE>\
  (const std::string& name, const TYPE& value) const;\
  template PtrParamNode ParamNode::TokenizedHasAndGet<TYPE>\
      (const std::string& name, const TYPE& value, bool ) const;
INSTANTIATE_METHODS_INT(Integer)
INSTANTIATE_METHODS_INT(Double)
INSTANTIATE_METHODS_INT(UInt)
INSTANTIATE_METHODS_INT(bool)
INSTANTIATE_METHODS_INT(std::string)

} // end of namespace


