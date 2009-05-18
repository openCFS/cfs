#include "DataInOut/ParamHandling/InfoNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/coloredConsole.hh"

#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>
#include <complex>

using std::string;

namespace CoupledField 
{
  
  // declare class specific logging stream
  DECLARE_LOG(infoNode)
  DEFINE_LOG(infoNode, "infoNode")


  /** This is the global instance of the global pointer */
  InfoNode* info = NULL;
  
  /** set the global names for fields */
  const string InfoNode::HEADER = "header";
  const string InfoNode::PROCESS = "process";
  const string InfoNode::SUMMARY = "summary";
  
  const string InfoNode::COMMENT = "comment";
  const string InfoNode::WARNING = "warning";
  const string InfoNode::ERROR   = "error";
  
  std::map<const string, unsigned int> InfoNode::writeCounter_ = std::map<const string, unsigned int>();

  InfoNode::InfoNode(const string& filename, const string& preamble) :
    parent_(NULL),
    log_(NOT_SET),
    cardinality_(USE_EXISTING), // default
    type_(UNKNOWN), // to be set with SetValue()
    filename_(filename),
    preamble_(preamble)
  {
  }
  
  InfoNode* InfoNode::Get(const string& org_name, const string& comment,
                          const string& child, const string& value, Cardinality card)
  {
    InfoNode* ptr = this;
    string&   name = const_cast<string&>(org_name);
    
    if(ContainsTokens(org_name))
    {
      // call us recursive!
      if(child != "" || value != "")
        EXCEPTION("Name '" << org_name << "' containts delemiter '/', which is not allowed when child and value are given");
      
      StdVector<string> tokens = SplitIntoTokens(org_name);
      
      for(unsigned int i = 0; i < tokens.GetSize()-1; i++)
      {
        name = tokens[i];
        ptr = ptr->Get(name, USE_EXISTING);
      }
      // set name to last element
      name = tokens[tokens.GetSize()-1];
      // setting the 'ptr' is done below!
    }
    
    if(name == "") EXCEPTION("Invalid name");
    
    // just in case
    string label = ToValidLabel(name); // might contain slashes
    
    // decide on cardinality
    // search for an existing InfoNode
    InfoNode* node = NULL;
    if(child == "" || value == "")  
    {
      if(ptr->Has(label)) node = dynamic_cast<InfoNode*>(ptr->ParamNode::Get(label, true)); // don't call this Get recursive!
    } 
    else 
    {
      if(Has(label, child, value)) node = dynamic_cast<InfoNode*>(ptr->ParamNode::Get(label, child, value, true));
    }
    
    switch(card)
    {
    case USE_EXISTING:
      if(node != NULL) return node;
      break;
      
    case REPLACE:
      if(node != NULL)
      {
        ptr->children_.Erase(ptr->children_.Find(node));
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
    if(node == NULL || card == APPEND)
    {
      node = new InfoNode();
      ptr->children_.Push_back(node);
      ptr->attribute_ = false;
    }
    
    node->parent_ = this;
    
    node->name_    = label; // see above!
    
    node->cardinality_ = card;
    
    if(child != "" && value != "")
      node->Get(child)->SetValue(value);
    
    // set comment if we have on
    if(comment != "")
    {
      InfoNode* in = node->Get(COMMENT);
      in->SetValue(comment);
      in->attribute_ = IsGoodAttribute(comment);
    }
    
    return node;
  }
  
  void InfoNode::ToXML(std::ostream& os, int depth) const
  {
    // note, that this is an recursive method!
    if(attribute_)
    {
      // makes only sense in an recursive call
      os << " " << name_ << "=\"" << GetFormatedValue() << "\"";
      return;
    }
    
    // we are NO Attribute
    // if we start a new element or are part of an element is same/same
    os << std::endl << string(depth, ' ') << "<" << name_;
    // go through all children an check if we close with "../>" or </name_>
    // note, that the own XML-element value can also be a attribute
    bool only_attributes = value_ != "" && !IsGoodAttribute(value_) ? false : true;
    
    for(unsigned int i = 0; i < children_.GetSize(); i++)
      if(!dynamic_cast<InfoNode*>(children_[i])->IsAttribute())
        only_attributes = false;
    
    // first loop, do the attributes
    for(unsigned int i = 0; i < children_.GetSize(); i++)
    {
      InfoNode* in = dynamic_cast<InfoNode*>(children_[i]);
      // process the attributs
      if(in->IsAttribute())
        in->ToXML(os);
      else
        only_attributes = false; // no, there is a non-attribute element
    }
    // the own value might be an attribute
    if(value_ != "" && IsGoodAttribute(value_))
      os << " value=\"" << GetFormatedValue() << "\"";
    
    
    if(only_attributes)
      os << "/>"; // "quick close"
    else
    {
      os << ">";
      // second loop, child elements
      for(unsigned int i = 0; i < children_.GetSize(); i++)
      {
        InfoNode* in = dynamic_cast<InfoNode*>(children_[i]);
        if(in->IsAttribute()) continue; // attributes are already done
        in->ToXML(os, depth + 1);
      }
      // own element
      if(value_ != "" && !IsGoodAttribute(value_))
      {
        if(name_ == COMMENT) os << value_;
        else os << std::endl << string(depth + 1, ' ') << "<value>" << GetFormatedValue() << "</value>" << std::endl;
      }
      // do we close in the same line?
      if(children_.GetSize() != 0) os << std::endl << string(depth, ' ');
      os << "</" << name_ << ">";
    }
  }
 
  void InfoNode::ToFile(const string& filename, const string& preamble)
  {
    if(filename != "") filename_ = filename;
    if(preamble != "") preamble_ = preamble;
    
    if(filename_ == "") EXCEPTION("Cannot write info.xml file, a filename was never given");
    
    std::ofstream info_file(filename_.c_str());
    info_file << preamble_;
    // store how often we are written -> if the number is too high one should cancel some ToFile() calls
    if(writeCounter_.count(filename_) == 0) writeCounter_[filename_] = 0;
    Get("writeCounter")->SetValue(++writeCounter_[filename_]);
    ToXML(info_file);
    info_file.close();
  }
  
  
  bool InfoNode::IsGoodAttribute(const string& attr) const
  {
    if(attr.size() > 30) return false;
    if(attr.find('\"') != string::npos) return false;
    return true;
  }
  
  string InfoNode::ToValidLabel(const string& in) const
  {
    string out = in;
    
    
    boost::trim(out);
    
    // remove spaces and don't make upper case yet! :(
    boost::erase_all(out, " ");
    boost::erase_all(out, "\"");
    boost::erase_all(out, "(");
    boost::erase_all(out, ")");
    // don't touch the slash '/', it is a 'xpath' element
    return out;
  }
  
  string InfoNode::GetFormatedValue() const
  {
    if(type_ != REAL && type_ != COMPLEX)
      return value_;
    
    
    std::stringstream ss;
    ss.precision(6);
    
    try
    {
      if(type_ == REAL)
        ss << boost::lexical_cast<double>(value_);
      else
        ss << boost::lexical_cast<std::complex<double> >(value_);
    }
    catch(bad_lexical_cast&)
    {
      // there is something like inf, ... give back the original string then
      LOG_DBG(infoNode) << "InfoNode::GetFormatedValue(): boost::lexical_cast failed to convert " << value_;
      return value_;
    }
    
    return ss.str();
  }
  
  void InfoNode::SetValue(ParamNode* node) 
  { 
    SetName(node->GetName()); 
    // set the value  
    SetValue(node->AsString()); 

    StdVector<ParamNode*>& children = node->GetChildren(); 
    attribute_ = children.GetSize() == 0; 

    // run recusively through all children 
    for(unsigned int i = 0; i < children.GetSize(); i++) 
    { 
      // add new element 
      ParamNode* other = children[i]; 
      InfoNode* new_node = Get(other->GetName(), APPEND); 
      new_node->SetValue(other); 
    } 
  } 

  void InfoNode::SetValue(const string& value)
  {
    value_ = value;
    type_ = STRING;
    attribute_ = IsGoodAttribute(value);
    
    // handle warnings also als console print
    if(name_ == WARNING)
    {
      std::cerr << fg_magenta << "Warning: /" << fg_reset;
      
      InfoNode* in = parent_; // beyond warning
      StdVector<string> tree;
      while(in != NULL)
      {
        tree.Push_back(in->name_);
        in = in->parent_;
      }
      for(unsigned int i = tree.GetSize() -1; i > 0; --i)
      {
        std::cerr << tree[i];
        if(i > 1) std::cerr << "/";
      }
      
      std::cerr << ": " << fg_magenta << "\"" << fg_reset << value << fg_magenta << "\"" << fg_reset << std::endl;
    }
  }
  
  void InfoNode::SetValue(bool param)
  {
    value_ = param ? "true" : "false";
    type_ = BOOL;
    attribute_ = true;
  }
  
  
  void InfoNode::SetValue(unsigned int param)
  {
    SetValue((int) param);
  }
  
  void InfoNode::SetValue(int param)
  {
    value_ = boost::lexical_cast<string>(param);
    type_ = INTEGER;
    attribute_ = true;
  }
  
  void InfoNode::SetValue(double param)
  {
    value_ = boost::lexical_cast<string>(param);
    type_ = REAL;
    attribute_ = true;
  }
  
  void InfoNode::SetValue(const std::complex<double>& param)
  {
    value_ = boost::lexical_cast<string>(param);
    type_ = COMPLEX;
    attribute_ = true;
  }
  
  void InfoNode::SetComment(const string& value)
  {
    // are we alrey a comment ?
    if(name_ == COMMENT) throw Exception("You cannot create a comment info node for a comment");
    
    InfoNode* in = Get(COMMENT);
    in->SetValue(value);
    in->attribute_ = IsGoodAttribute(value);
  }
  
}
