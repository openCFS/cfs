#include "DataInOut/ParamHandling/InfoNode.hh"
#include "DataInOut/Logging/cfslog.hh"
#include "DataInOut/coloredConsole.hh"
#include "MatVec/denseMatrix.hh"
#include "Utils/Timer.hh"
#include "MatVec/matrix.hh"

#include <boost/algorithm/string.hpp>
#include <iostream>
#include <fstream>
#include <complex>

using std::string;
using std::complex;

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
  
  const string DATA_IS_IN_MATRIX_VALUE = "data is in matrix_value_!";
  const string DATA_IS_IN_TIMER_VALUE  = "data is in timer_value_!";

  std::map<const string, unsigned int> InfoNode::writeCounter_ = std::map<const string, unsigned int>();

  InfoNode::InfoNode(const string& filename, const string& preamble) :
    parent_(NULL),
    log_(NOT_SET),
    cardinality_(USE_EXISTING), // default
    type_(UNKNOWN), // to be set with SetValue()
    matrix_value_(NULL),
    timer_value_(NULL),
    filename_(filename),
    preamble_(preamble)
  {
  }
  
  InfoNode::~InfoNode()
  {
    delete matrix_value_; matrix_value_ = NULL;
    delete timer_value_;  timer_value_ = NULL;
  }

  InfoNode* InfoNode::Get(const string& org_name, const string& comment,
                          const string& child, const string& value, Cardinality card)
  {
    InfoNode* ptr = this;
    string&   name = const_cast<string&>(org_name);
    
    if(ContainsTokens(org_name))
    {
      // call us recursive!
      if(!child.empty() || !value.empty())
        EXCEPTION("Name '" << org_name << "' contains delimiter '/', which is not allowed when child and value are given");
      
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
    
    if(name.empty()) EXCEPTION("Invalid name");
    
    // just in case
    string label = ToValidLabel(name); // might contain slashes
    
    // decide on cardinality
    // search for an existing InfoNode
    InfoNode* node = NULL;
    if(card != APPEND) // for APPEND we need a new node anyway
    {
      if(child.empty() || value.empty())  
        node = dynamic_cast<InfoNode*>(ptr->ParamNode::Get(label, false)); // don't call this Get recursive!
      else 
        node = dynamic_cast<InfoNode*>(ptr->ParamNode::Get(label, child, value, false));
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
    
    if(!child.empty() && !value.empty())
      node->Get(child)->SetValue(value);
    
    // set comment if we have on
    if(!comment.empty())
    {
      InfoNode* in = node->Get(COMMENT);
      in->SetValue(comment);
      in->attribute_ = IsGoodAttribute(comment);
    }
    
    return node;
  }
  
  InfoNode* InfoNode::SetNewChild(const string& name, unsigned int index)
  {
    InfoNode* in = new InfoNode();
    children_[index] = in;
    in->parent_ = this;
    in->name_   = name;

    return in;
  }


  void InfoNode::ToXML(std::ostream& os, int depth)
  {
    // note, that this is an recursive method!
    if(attribute_)
    {
      // makes only sense in an recursive call
      os << " " << name_ << "=\"" << GetFormatedValue(depth) << "\"";
      return;
    }
    
    // we are NO Attribute

    Sort(); // (HEADER before) PROCESS before SUMMARY

    // if we start a new element or are part of an element is same/same
    os << std::endl << string(depth, ' ') << "<" << name_;
    // go through all children an check if we close with "../>" or </name_>
    // note, that the own XML-element value can also be a attribute
    bool only_attributes = (!value_.empty() && !IsGoodAttribute(value_) ? false : true);
    
    const unsigned int chsize(children_.GetSize());
    for(unsigned int i = 0; i < chsize; i++)
      if(!dynamic_cast<InfoNode*>(children_[i])->IsAttribute())
        only_attributes = false;
    
    // first loop, do the attributes
    for(unsigned int i = 0; i < chsize; i++)
    {
      InfoNode* in = dynamic_cast<InfoNode*>(children_[i]);
      // process the attributes
      if(in->IsAttribute())
        in->ToXML(os);
      else
        only_attributes = false; // no, there is a non-attribute element
    }
    // the own value might be an attribute
    if(!value_.empty() && IsGoodAttribute(value_))
      os << " value=\"" << GetFormatedValue(depth) << "\"";
    
    
    if(only_attributes)
      os << "/>"; // "quick close"
    else
    {
      os << ">";
      // second loop, child elements
      for(unsigned int i = 0; i < chsize; i++)
      {
        InfoNode* in = dynamic_cast<InfoNode*>(children_[i]);
        if(in->IsAttribute()) continue; // attributes are already done
        in->ToXML(os, depth + 2);
      }
      bool endl_written = false;
      // own element
      if(!value_.empty() && !IsGoodAttribute(value_))
      {
        if(name_ == COMMENT)
        {
          os << value_;
        }
        else
        {
          os << std::endl << string(depth + 2, ' ') << GetFormatedValue(depth+2) << std::endl;
          endl_written = true;
        }
      }
      // do we close in the same line?
      if(chsize != 0 && !endl_written) os << std::endl;
      os << string(depth, ' ') << "</" << name_ << ">";
    }
  }
 
  void InfoNode::ToFile(const string& filename, const string& preamble)
  {
    if(!filename.empty()) filename_ = filename;
    if(!preamble.empty()) preamble_ = preamble;
    
    if(filename_.empty()) EXCEPTION("Cannot write info.xml file, a filename was never given");
    
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
    if(attr == DATA_IS_IN_MATRIX_VALUE) return false;
    if(attr == DATA_IS_IN_TIMER_VALUE)  return false;

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
  

  string InfoNode::GetFormatedValue(unsigned int depth) const
  {
    switch(type_)
    {
      case REAL:
      case COMPLEX:
      {
        std::stringstream stream;
        stream.precision(6);

        try
        {
          if(type_ == REAL)
            stream << boost::lexical_cast<double>(value_);
          else
            stream << boost::lexical_cast<complex<double> >(value_);
        }
        catch(bad_lexical_cast&)
        {
          // there is something like inf, ... give back the original string then
          LOG_DBG(infoNode) << "InfoNode::GetFormatedValue(): boost::lexical_cast failed to convert " << value_;
          stream << value_;
        }

        return stream.str();
      }

      case MATRIX:
        return matrix_value_->ToXML(depth);

      case TIMER:
        return timer_value_->ToXML();

      default:
        return value_;
    }
    
  }

  bool InfoNode::SortCheck(const int idx, int& this_idx, int& other_idx)
  {
    if(other_idx != -1 && other_idx < idx)
    {
      children_.Swap(idx, other_idx);
      this_idx  = other_idx;
      other_idx = idx;
      return true;
    }
    return false;
  }

  /** First the order of HEADER, PROCESS and SUMMARY is ensured
   * then HEADER is set first and SUMMARY last.
   * If there is a performance issue this could be done more complex ?! */
  void InfoNode::Sort()
  {
    int head_idx = -1;
    int proc_idx = -1;
    int summ_idx = -1;

    for(unsigned int i = 0; i < children_.GetSize(); i++)
    {
      ParamNode* cand = children_[i];
      if(cand->GetName() == InfoNode::HEADER)
      {
        if(SortCheck(i, head_idx, proc_idx))
        {
          i--;
          continue;
        }

        if(SortCheck(i, head_idx, summ_idx))
        {
          i--;
          continue;
        }
      }

      if(cand->GetName() == InfoNode::PROCESS)
      {
        if(SortCheck(i, proc_idx, summ_idx))
        {
          i--;
          continue;
        }
      }

      // all sorting done
      if(cand->GetName() == InfoNode::SUMMARY)
        summ_idx = i;
    }

    // not the order is correct. Assure HEADER is first and Summary LAST
    if(head_idx != -1 && head_idx != 0)
      children_.Swap(0, head_idx);

    if(summ_idx != -1 && summ_idx != (int) children_.GetSize() - 1)
      children_.Swap(children_.GetSize() - 1, summ_idx);
  }

  void InfoNode::SetValue(ParamNode* node)
  { 
    SetName(node->GetName()); 
    // set the value  
    SetValue(node->AsString()); 

    StdVector<ParamNode*>& children = node->GetChildren();
    const unsigned int chsize(children_.GetSize());
    attribute_ = chsize == 0; 

    // run recusively through all children 
    for(unsigned int i = 0; i < chsize; i++) 
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
  
  void InfoNode::SetValue(long unsigned int param)
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
  
  void InfoNode::SetValue(const complex<double>& param)
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

  DenseMatrix* InfoNode::SetValue(DenseMatrix* matrix)
  {
    assert(matrix != NULL);

    value_        = DATA_IS_IN_MATRIX_VALUE;
    matrix_value_ = matrix;
    type_         = MATRIX;
    attribute_    = false;

    return matrix_value_;
  }

  Timer* InfoNode::SetValue(Timer* timer)
  {
    assert(timer != NULL);

    value_       = DATA_IS_IN_TIMER_VALUE;
    timer_value_ = timer;
    type_        = TIMER;
    attribute_   = false;

    return timer_value_;
  }

  DenseMatrix* InfoNode::AsMatrix()
  {
    assert(type_ == MATRIX);
    return matrix_value_;
  }

  Timer* InfoNode::AsTimer()
  {
    assert(type_ == TIMER);
    return timer_value_;
  }
  
}
