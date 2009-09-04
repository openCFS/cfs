// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include <boost/tokenizer.hpp>
#include "DataInOut/ParamHandling/ParamNode.hh"

using namespace std;
using namespace boost;

namespace CoupledField
{
    
  /** This is our global pointer of the root ParamNode holding the XML file. 
   *  Filed in cfs.cc. The corresponding
   * "extern ParamNode* param;" is in ParamNode.hh */
  ParamNode* param = NULL;
    
  ParamNode::ParamNode(bool attribute) :
    attribute_(attribute),
    lastresultidx_(-1)
  {} 
    
    
  ParamNode::~ParamNode()
  {
    // values are no problem, but we have to delete our childs recursively
    for(unsigned int i = 0, chsize = children_.GetSize(); i < chsize; ++i)
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
     
  ParamNode* ParamNode::Get(const string& name, const bool throwException)
  {
    ParamNode* result = NULL;
    
    if(ContainsTokens(name))
    {
       result = TokenizedHasAndGet(name, "", false); // has some overhead   
    }
    else
    {
      const int chsize(children_.GetSize());
      
      // check if the last result is close to the currently requested node
      if(lastresultidx_ > 0 && lastresultidx_ < chsize)
      {
        if(children_[lastresultidx_]->name_ == name)
          result = children_[lastresultidx_];
        else
          if(lastresultidx_ + 1 < chsize)
          {
            if(children_[lastresultidx_ + 1]->name_ == name)
              result = children_[++lastresultidx_];
          }   
      }
      
      for(int i = 0; i < chsize && result == NULL; i++)
      {
        if(children_[i]->name_ == name)
        {
          lastresultidx_ = i;
          result = children_[i];
        }
      }
    }
    
    if(result == NULL && throwException) 
      EXCEPTION("None of the " << children_.GetSize() << " childs of element '"
                << this->name_ << "'" << " has a child '" << name << "'");
    
    return result;
  }

  void ParamNode::Get(const string&  name, string& ret, const bool throwException)
  {
    ParamNode * actNode = Get(name,  throwException);
    if(actNode)
      ret = actNode->AsString();
  }
   

  void ParamNode::Get(const string&  name, int& ret, const bool throwException)
  {
    ParamNode * actNode = Get(name, throwException);
    if( actNode )
      ret = actNode->AsInt();
  }

  void ParamNode::Get(const string&  name, unsigned int& ret, const bool throwException)
  {
    ParamNode * actNode = Get(name, throwException);
    if( actNode )
      ret = actNode->AsUInt();
  }

  void ParamNode::Get(const string&  name, double& ret, const bool throwException)
  {
    ParamNode * actNode = Get(name, throwException);
    if( actNode )
      ret = actNode->AsDouble();
  }

  void ParamNode::Get(const string&  name, bool& ret, const bool throwException)
  {
    ParamNode * actNode = Get(name, throwException);
    if( actNode )
      ret = actNode->AsBool();
  }

  
  ParamNode* ParamNode::TokenizedHasAndGet(const string& name, const string& value, bool has_bool_value, bool bool_value) const
  {
    StdVector<string> tokens = SplitIntoTokens(name);
    
    ParamNode* ptr = const_cast<ParamNode*>(this);
    for(unsigned int i = 0; i < tokens.GetSize(); i++)
    {
      // ptr->Dump();
      // cout << "token = " << tokens[i] << endl;
      
      // consider the value only for the last token!
      if(i == tokens.GetSize()-1 && (!value.empty() || has_bool_value))
      {
        if(has_bool_value)
        {
          if(!ptr->Has(tokens[i], bool_value)) return NULL;
        }
        else
        {
          if(!ptr->Has(tokens[i], value)) return NULL;
        }
      }
      else
      {
        if(!ptr->Has(tokens[i])) return NULL;
      }
      ptr = ptr->Get(tokens[i]);
    }
    
    return ptr;
  }
  
  
  bool ParamNode::Has(const string& name) const
  {
    // check in a fast way if we have tokens for trivial xpath
    if(ContainsTokens(name)) 
    {
      return TokenizedHasAndGet(name, string(""), false, false) != NULL;
    }
    else
    {
      for(unsigned int i = 0, chsize = children_.GetSize(); i < chsize; i++)
        if(children_[i]->name_ == name) return true;

      return false;
    }
  }

  bool ParamNode::Has(const string& name, const string& value) const
  {
    if(ContainsTokens(name))
    {
      return TokenizedHasAndGet(name, value, false) != NULL;
    }
    else
    {
      for(unsigned int i = 0, chsize = children_.GetSize(); i < chsize; i++)
        if(children_[i]->name_ == name && children_[i]->value_ == value) return true;

      return false;
    }
  }

  bool ParamNode::Has(const string& name, bool value) const
  {
    if(ContainsTokens(name))
    {
      return TokenizedHasAndGet(name, "", true, value) != NULL;
    }
    else
    {
      for(unsigned int i = 0, chsize = children_.GetSize(); i < chsize; i++)
        if(children_[i]->name_ == name && children_[i]->AsBool() == value) return true;

      return false;
    }
  }

  
  StdVector<ParamNode*> ParamNode::GetList(const string& name)
  {
    const unsigned int chsize(children_.GetSize());
    StdVector<ParamNode*> result;
    result.Reserve(chsize);

    for(unsigned int i = 0; i < chsize; ++i)
      if(children_[i]->name_ == name)
        result.Push_back(children_[i]);
   
    return result; // copy-constructor magic stuff!
  }

  unsigned int ParamNode::Count(const string& name) const
  {
    unsigned int count = 0;

    for(unsigned int i = 0, chsize = children_.GetSize(); i < chsize; ++i)
      if(children_[i]->name_ == name) ++count;
   
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
  
  
  StdVector<ParamNode*> ParamNode::GetList(const string& parent, const string& child, const string& value)
  {
    StdVector<ParamNode*> result;
    
    for(unsigned int p = 0, chsize = children_.GetSize(); p < chsize; ++p)
    {
      ParamNode *ch = children_[p];
      // do we have parent name? 
      if(ch->GetName() == parent)
      {
        // children_[i] has zero, one or more child matching child elements,
        // we loop by hand to handle the "more" case
        const unsigned int chchsize(ch->children_.GetSize());
        for(unsigned int c = 0; c < chchsize; ++c)
        {
          const ParamNode *chch = ch->children_[c];
          if(chch->GetName() == child && chch->AsString() == value)
          {
            // parent, child and value matches. Take it!
            result.Push_back(ch);
          }
        }
      }
    }

    return result; // copy-constructor magic stuff!
  }

  ParamNode* ParamNode::Get(const string& parent, const string& child, 
                            const string& value, const bool throwException)
  {
    StdVector<ParamNode*> result = GetList(parent, child, value);
    
    if(result.GetSize() == 1)
      return result[0];
    
    if(result.GetSize() > 1)
    {
      if(!throwException)
        return result[0];
      else
        EXCEPTION("element '" << parent << "' has more than one (" << result.GetSize() 
                   << ") childs '" << child << "' with value '" << value << "'");
    }
    
    // reslt size = 0!
    if(throwException)
      EXCEPTION("element '" << parent << "' has no child '" << child
                << "' with value '" << value << "'");
    return NULL;
  }
  

  bool ParamNode::Has(const string& parent, const string& child, 
                      const string& value ) const
  {
    // see GetList() for comments
    for(unsigned int p = 0, chsize = children_.GetSize(); p < chsize; ++p)
    {
      ParamNode *ch = children_[p];
      if(ch->name_ == parent)
      {
        const unsigned int chchsize(ch->children_.GetSize());
        for(unsigned int c = 0; c < chchsize; ++c)
        {
          ParamNode *chch(ch->children_[c]);
          if(chch->GetName() == child && chch->AsString() == value)
            return true;
        }
      }
    }
    
    // nothing found
    return false;
  }

  string ParamNode::ToString() const
  {
    ostringstream os;
    os << name_ << " = '" << value_ << "'" << " attribute: " << (attribute_ ? "true" : "false");
    return os.str();
  }
       
  void ParamNode::Dump(int level) const
  {
    for(int i = 0; i < level; i++) cout << "   ";
    cout << ToString() << endl;
    
    for(unsigned int i = 0, chsize = children_.GetSize(); i < chsize; i++) 
      children_[i]->Dump(level + 1);
  }

  StdVector<string> ParamNode::SplitIntoTokens(const string& input) const
  {
    StdVector<string> out;
    char_separator<char> sep("/");
    tokenizer<char_separator<char> > tokens(input, sep);
    
    for(tokenizer<char_separator<char> >::iterator it = tokens.begin(); it != tokens.end(); ++it)
      out.Push_back(*it);
    
    return out;
  }
} // end of namespace


