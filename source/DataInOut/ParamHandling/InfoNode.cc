#include "DataInOut/ParamHandling/InfoNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/coloredConsole.hh"

#include <boost/algorithm/string.hpp>
#include <iostream>
#include<fstream>
#include <complex>



namespace CoupledField 
{
  
  // declare class specific logging stream
  DECLARE_LOG(infoNode)
  DEFINE_LOG(infoNode, "infoNode")


  /** This is the global instance of the global pointer */
  InfoNode* info = NULL;
  
  /** set the global names for fields */
  const std::string InfoNode::HEADER = "header";
  const std::string InfoNode::PROCESS = "process";
  const std::string InfoNode::SUMMARY = "summary";
  
  const std::string InfoNode::COMMENT = "comment";
  const std::string InfoNode::WARNING = "warning";
  const std::string InfoNode::ERROR   = "error";
  


  InfoNode::InfoNode(const std::string& filename, const std::string& preamble) :
    parent_(NULL),
    log_(NOT_SET),
    cardinality_(USE_EXISTING), // default
    type_(UNKNOWN), // to be set with SetValue()
    filename_(filename),
    preamble_(preamble)
  {
  }


  InfoNode* InfoNode::Get(const std::string& name, const std::string& comment,
                          const std::string& child, const std::string& value, Cardinality card)
  {
    if(name == "") EXCEPTION("Invalid name");
    
    // just in case
    std::string label = ToValidLabel(name);
    
    // decide on cardinality
    // search for an existing InfoNode
    InfoNode* node = NULL;
    if(child == "" || value == "")  {
      if(Has(label)) node = dynamic_cast<InfoNode*>(ParamNode::Get(label, true)); // don't call this Get recursive!
    } else {
      if(Has(label, child, value)) node = dynamic_cast<InfoNode*>(ParamNode::Get(label, child, value, true));
    }
    
    switch(card)
    {
    case USE_EXISTING:
      if(node != NULL) return node;
      break;
      
    case REPLACE:
      if(node != NULL)
      {
        children_.Erase(children_.Find(node));
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
      children_.Push_back(node);
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
    os << std::endl << std::string(depth, ' ') << "<" << name_;
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
        else os << std::endl << std::string(depth + 1, ' ') << "<value>" << GetFormatedValue() << "</value>" << std::endl;
      }
      // do we close in the same line?
      if(children_.GetSize() != 0) os << std::endl << std::string(depth, ' ');
      os << "</" << name_ << ">";
    }
  }
  
  void InfoNode::ToFile(const std::string& filename, const std::string& preamble)
  {
    if(filename != "") filename_ = filename;
    if(preamble != "") preamble_ = preamble;
    
    if(filename_ == "") EXCEPTION("Cannot write info.xml file, a filename was never given");
    
    std::ofstream info_file(filename_.c_str());
    info_file << preamble_;
    info->ToXML(info_file);
    info_file.close();
  }
  
  
  bool InfoNode::IsGoodAttribute(const std::string& attr) const
  {
    if(attr.size() > 30) return false;
    if(attr.find('\"') != std::string::npos) return false;
    return true;
  }
  
  std::string InfoNode::ToValidLabel(const std::string& in) const
  {
    std::string out = in;
    
    
    boost::trim(out);
    
    // remove spaces and don't make upper case yet! :(
    boost::erase_all(out, " ");
    boost::erase_all(out, "\"");
    boost::erase_all(out, "(");
    boost::erase_all(out, ")");
    
    return out;
  }
  
  std::string InfoNode::GetFormatedValue() const
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
  
  void InfoNode::SetValue(const std::string& string)
  {
    value_ = string;
    type_ = STRING;
    attribute_ = IsGoodAttribute(string);
    
    // handle warnings also als console print
    if(name_ == WARNING)
    {
      std::cerr << fg_magenta << "Warning: /" << fg_reset;
      
      InfoNode* in = parent_; // beyond warning
      StdVector<std::string> tree;
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
      
      std::cerr << ": " << fg_magenta << "\"" << fg_reset << string << fg_magenta << "\"" << fg_reset << std::endl;
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
    value_ = boost::lexical_cast<std::string>(param);
    type_ = INTEGER;
    attribute_ = true;
  }
  
  void InfoNode::SetValue(double param)
  {
    value_ = boost::lexical_cast<std::string>(param);
    type_ = REAL;
    attribute_ = true;
  }
  
  void InfoNode::SetValue(const std::complex<double>& param)
  {
    value_ = boost::lexical_cast<std::string>(param);
    type_ = COMPLEX;
    attribute_ = true;
  }
  
  void InfoNode::SetComment(const std::string& string)
  {
    // are we alrey a comment ?
    if(name_ == COMMENT) throw Exception("You cannot create a comment info node for a comment");
    
    InfoNode* in = Get(COMMENT);
    in->SetValue(string);
    in->attribute_ = IsGoodAttribute(string);
  }
  
}
