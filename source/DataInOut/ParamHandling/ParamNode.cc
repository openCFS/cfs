// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/exception.hh"
#include "boost/lexical_cast.hpp"

namespace CoupledField
{
    
  /** This is our global pointer of the root ParamNode holding the XML file. 
   *  Filed in cfs.cc. The correpsonding
   * "extern ParamNode* param;" is in ParamNode.hh */
  ParamNode* param;     
    
  ParamNode::ParamNode(bool attribute)
  {
    this->attribute_ = attribute;
  } 
    
    
  ParamNode::~ParamNode()
  {
    // values are no problem, but we have to delete our childs recursively
    for(unsigned int i = 0; i < children_.GetSize(); i++)
      {
        if(children_[i] != NULL) 
          { 
            delete children_[i]; 
            children_[i] = NULL; 
          }
      }
  }
  
  int ParamNode::AsInt() const
  {
    try
    {
      return boost::lexical_cast<int>(value_);
    }
    catch(boost::bad_lexical_cast &)
    {
      EXCEPTION("cannot cast value '" << value_ << "' in XML node " << name_
                << " to an integer.");
    }
  }
  
  unsigned int ParamNode::AsUInt() const
  {
    try
    {
      return boost::lexical_cast<unsigned int>(value_);
    }
    catch(boost::bad_lexical_cast &)
    {
      EXCEPTION("cannot cast value '" << value_ << "' in XML node " << name_
                << " to an unsigned integer.");
    }
  }

  double ParamNode::AsDouble() const
  {
    try
    {
      return boost::lexical_cast<double>(value_);
    }
    catch(boost::bad_lexical_cast &)
    {
      EXCEPTION("cannot cast value '" << value_ << "' in XML node " << name_
                << " to a double.");
    }
  }
       
  bool ParamNode::AsBool() const
  {
    if(value_ == "yes" || value_ == "true"  || value_ == "on")  return true;
    if(value_ == "no"  || value_ == "false" || value_ == "off") return false;
    
    EXCEPTION("cannot interpret " << value_ << " as boolean");
  }
     
  ParamNode* ParamNode::Get(const std::string& name, const bool throwException)
  {
    for(unsigned int i = 0; i < children_.GetSize(); i++)
      {
        if(children_[i]->name_ == name) return children_[i];
      }

    if( throwException) 
      {
        EXCEPTION("None of the " << children_.GetSize() << " childs of element '"
                  << this->name_ << "'" << " has a child '" << name << "'");   
      } else {
      return NULL;
    }
  }

  void ParamNode::Get(const std::string&  name, std::string& ret, const bool throwException)
  {
    ParamNode * actNode = Get(name,  throwException);
    if( actNode ) {
      ret = actNode->AsString();
    }
  }
   

  void ParamNode::Get(const std::string&  name, int& ret, const bool throwException)
  {
    ParamNode * actNode = Get(name, throwException);
    if( actNode ) {
      ret = actNode->AsInt();
    }
    
  }

  void ParamNode::Get(const std::string&  name, unsigned int& ret, const bool throwException)
  {
    ParamNode * actNode = Get(name, throwException);
    if( actNode ) {
      ret = actNode->AsUInt();
    }

  }

  void ParamNode::Get(const std::string&  name, double& ret, const bool throwException)
  {
    ParamNode * actNode = Get(name, throwException);
    if( actNode ) {
      ret = actNode->AsDouble();
    }
  }

  void ParamNode::Get(const std::string&  name, bool& ret, const bool throwException)
  {

    ParamNode * actNode = Get(name, throwException);
    if( actNode ) {
      ret = actNode->AsBool();
    }
  }

  bool ParamNode::Has(const std::string& name) const
  {
    for(unsigned int i = 0; i < children_.GetSize(); i++)
      if(children_[i]->name_ == name) return true;

    return false;
  }

     
  StdVector<ParamNode*> ParamNode::GetList(const std::string& name)
  {
    StdVector<ParamNode*> result;

    for(unsigned int i = 0; i < children_.GetSize(); i++)
      if(children_[i]->name_ == name) result.Push_back(children_[i]);
   
    return result; // copy-constructor magic stuff!
  }

  unsigned int ParamNode::Count(const std::string& name) const
  {
    unsigned int count = 0;

    for(unsigned int i = 0; i < children_.GetSize(); i++)
      if(children_[i]->name_ == name) count++;
   
    return count;
  }


  ParamNode* ParamNode::GetChild()
  {

    // check that only one child is present
    if( children_.GetSize() != 1 )
      {
        EXCEPTION( "Element '" << name_ << "' must only have 1 child");
      }
    return children_[0];
  }


  StdVector<ParamNode*> ParamNode::GetList(const std::string& parent, const std::string& child, const std::string& value)
  {
    StdVector<ParamNode*> result;

    for(unsigned int p = 0; p < children_.GetSize(); p++)
      {
        // do we have parent name? 
        if(children_[p]->name_ == parent) 
          {
            // children_[i] has zero, one or more child matching child elements, 
            // we loop by hand to handle the "more" case
            for(unsigned int c = 0; c < children_[p]->children_.GetSize(); c++)
              {
                if(children_[p]->children_[c]->GetName() == child && children_[p]->children_[c]->AsString() == value)
                  { 
                    // parent, child and value matches. Take it!
                    result.Push_back(children_[p]);
                  }
              }
          }
      }
       
   
    return result; // copy-constructor magic stuff!
  }

  ParamNode* ParamNode::Get(const std::string& parent, const std::string& child, 
                            const std::string& value, const bool throwException)
  {
    ParamNode * result = NULL;
    
    for(unsigned int p = 0; p < children_.GetSize(); p++)
      {
        // do we have parent name? 
        if(children_[p]->name_ == parent) 
          {
            // children_[i] has zero, one or more child matching child elements, 
            // we loop by hand to handle the "more" case
            for(unsigned int c = 0; c < children_[p]->children_.GetSize(); c++)
              {
                if(children_[p]->children_[c]->GetName() == child && children_[p]->children_[c]->AsString() == value)
                  { 
                    // parent, child and value matches. Now check, if value was already set or
                    // throw exception
                    if( result != NULL && throwException == true) {
                      EXCEPTION( "More than one match found" );
                    }
                    result = children_[p];
                    break;
                  }
              }
          }
      }
       
    if( result == NULL && throwException == true) {
      EXCEPTION( "No match found" );
    }
    return result;
  }

  bool ParamNode::Has(const std::string& parent, const std::string& child, 
                      const std::string& value ) const
  {
    for(unsigned int p = 0; p < children_.GetSize(); p++)
      {
        // do we have parent name? 
        if(children_[p]->name_ == parent) 
          {
            // children_[i] has zero, one or more child matching child elements, 
            // we loop by hand to handle the "more" case
            for(unsigned int c = 0; c < children_[p]->children_.GetSize(); c++)
              {
                if(children_[p]->children_[c]->GetName() == child && children_[p]->children_[c]->AsString() == value)
                  { 
                    return true;
                  }
              }
          }
      }
    
    return false;
  }
  

  unsigned int ParamNode::Count(const std::string& parent, const std::string& child, const std::string& value)  
  {
    // better slow than copy-paste
    StdVector<ParamNode*> result = GetList(parent, child, value);
    return result.GetSize();
  }

  void ParamNode::ToXML(std::ostream& os) const
  {
    // note, that this is an recursive method!
    if(attribute_)
    {
      // makes only sense in an recursive call
      os << " " << name_ << "=\"" << value_ << "\"";
    }
    else
    {
      // if we start a new element or are part of an element is same/same
      os << "<" << name_;
      // go through all children an check if we close with "../>" or </name_>
      bool only_attributes = true;
      for(unsigned int i = 0; i < children_.GetSize(); i++)
      {
        if(!children_[i]->attribute_)
        {
           only_attributes = false;
           os << "/>";
        }
        children_[i]->ToXML(os); // I love recursive calls :)
      }
      
      // we are no attribute, but we might have a own value "<posDef>no</posDef>"
      if(value_ != "")
      {
        // see if we already close the beginning element
        if(only_attributes)
        {
          only_attributes = false;
          os << "/>";
        }
        os << value_;
      }
      
      // how to close
      if(only_attributes) os << "/>" << std::endl;
                     else os << "</" << name_ << ">" << std::endl;
    }
  }


  std::string ParamNode::ToString() const
  {
    std::ostringstream os;
    os << name_ << " = '" << value_ << "'" << " attribute: " << (attribute_ ? "true" : "false");
    return os.str();
  }
       
  void ParamNode::Dump(int level) const
  {
    for(int i = 0; i < level; i++) std::cout << "   ";
    std::cout << ToString() << std::endl;
    
    for(unsigned int i = 0; i < children_.GetSize(); i++) 
      children_[i]->Dump(level + 1);
  }


} // end of namespace


