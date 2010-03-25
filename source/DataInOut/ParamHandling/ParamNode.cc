// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#include "DataInOut/ParamHandling/ParamNode.hh"
#include "DataInOut/programOptions.hh"
#include "General/defs.hh"
#include "MatVec/matrix.hh"
#include "MatVec/vector.hh"
#include "Utils/Timer.hh"

#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include <fstream>

using namespace std;
using namespace boost;

namespace CoupledField
{

/** set the global names for fields */
  const string ParamNode::HEADER = "header";
  const string ParamNode::PROCESS = "process";
  const string ParamNode::SUMMARY = "summary";
  
  const string ParamNode::WARNING = "warnings";
  const string ParamNode::ERROR   = "error";
  
  /** This is our global pointer of the root ParamNode holding the XML file. 
   *  Filed in cfs.cc. The corresponding
   * "extern PtrParamNode param;" is in ParamNode.hh */
  PtrParamNode param;
  PtrParamNode info;
      
  ParamNode::ParamNode(ActionType defaultAction, NodeType type ) :
      name_(""),
      type_(type),
      defaultAction_(defaultAction),
      lastresultidx_(-1)
  {}    
    
  ParamNode::~ParamNode()
  {
   // explicit delete is not needed anymore
  }
  /************************************************************************
   * S E T    M E T H O D S
   *************************************************************************/
  
  void ParamNode::SetValue(const boost::any& value ) {
    this->value_ = value;
  }
  
  void ParamNode::SetValue(const char* value ) {
    std::string toSet(value);
    this->value_ = toSet;
  }    

  void ParamNode::SetComment(const std::string& comment) {
    // add new child node with type COMMENT
    PtrParamNode newChild (new ParamNode(defaultAction_, COMMENT));
    newChild->SetName("comment");
    newChild->SetValue(comment);
    newChild->parent_ = shared_from_this();
    newChild->defaultAction_ = defaultAction_;
    children_.Push_back(newChild);
  }

   
  void ParamNode::AddChildNode( PtrParamNode child) {
    child->parent_ = shared_from_this();
    if( child->defaultAction_ == DEFAULT ) {
      child->defaultAction_ = defaultAction_;
    }
    children_.Push_back(child);
  }
  
  PtrParamNode ParamNode::SetNewChild(const std::string& name, 
                                      unsigned int index ) {
    
    PtrParamNode node (new ParamNode());
    node->SetName(name);
    node->parent_ = shared_from_this();
    node->defaultAction_ = defaultAction_;
    children_[index] = node;
    return node;
  }
  
  /************************************************************************
  * N O D E   A C C E S S     M E T H O D S
  ************************************************************************/
  PtrParamNode ParamNode::GetChild()
    {
      // check that only one child is present
      if( children_.GetSize() != 1 )
      {
        EXCEPTION( "Element '" << name_ << "' must only have 1 child");
      }
      return children_[0];
    }
  
  PtrParamNode ParamNode::Get(const string& name, 
                              ActionType action)
  {
    PtrParamNode result;
    std::string myName = name;
    if( action == DEFAULT ) 
      action = defaultAction_;

    // check if name contains special "/" signs to
    if(ContainsTokens(name))
    {
      
      
      //result = TokenizedHasAndGet(name, std::string(""), false ); // has some overhead
      StdVector<string> tokens = SplitIntoTokens(name);
      PtrParamNode actNode = shared_from_this();
      for(unsigned int i = 0; i < tokens.GetSize(); i++)
      {
        myName = tokens[i];

        // here we have to be careful: if we have the actiontype APPEND,
        // we want to perform this just on the last element and re-use
        // the previous ones
        if( action == APPEND && i !=  tokens.GetSize()-1) {
          action = INSERT;
        }
        // check, if we still have a result
        if( actNode ) {
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
      // This is a harc-coded optimization coming from F. Wein
      if(lastresultidx_ > 0 && lastresultidx_ < chsize)
      {
        if(children_[lastresultidx_]->name_ == myName)
          result = children_[lastresultidx_];
        else
          if(lastresultidx_ + 1 < chsize)
          {
            if(children_[lastresultidx_ + 1]->name_ == myName)
              result = children_[++lastresultidx_];
          }   
      }

      for(int i = 0; i < chsize && result == NULL; i++)
      {
        if(children_[i]->name_ == myName)
        {
          lastresultidx_ = i;
          result = children_[i];
        }
      }
    }
    
    // perform action, in case no node was found
    if(result == NULL || action == APPEND) 
    {
      // depending on ActionType:
      switch(action)
      {
        case DEFAULT:
          EXCEPTION("Some action has to be performed");
          break;
        case EX:
          EXCEPTION("None of the " << children_.GetSize() << " childs of element '"
                          << this->name_ << "'" << " has a child '" << myName << "'");
          break;
        case PASS:
          break;
        case INSERT:
        case APPEND:
          // create new node containing "default" 
          PtrParamNode newChild (new ParamNode(defaultAction_));
          newChild->SetName(myName);
          // ATTENTION: Do NOT set an empty string as value to this node, as 
          // std::string("").empty() != boost::any(st::string("")).empty()
          newChild->parent_ = shared_from_this();
          newChild->defaultAction_ = defaultAction_;
          children_.Push_back(newChild);
          result = newChild;
          break;
      }
    }

    return result;
  }
  
  template<typename TYPE>
  PtrParamNode ParamNode::GetByVal(const string& parent, 
                                   const string& child, 
                                   const TYPE& value, 
                                   ActionType action )
  {

    // Strategy depending on action:
    // EX / PASS / INSERT: Search for node. If not found 
    //   EX / PASS : throw exception / pass
    //   INSERT: INSERT content with current one
    // APPEND: Directly insert node, regardless of existing ones
    if( action == DEFAULT ) action = defaultAction_;
    
    bool insertNew = false;
    if( action != APPEND) {
      StdVector<PtrParamNode> result = GetListByVal(parent, child, value);

      if(result.GetSize() == 1)
        return result[0];

      if(result.GetSize() > 1)
      {
        if( action == PASS || action == INSERT )
          return result[0];
        else
          EXCEPTION("element '" << parent << "' has more than one (" << result.GetSize() 
                    << ") childs '" << child << "' with value '" << value << "'");
      }

      // result size = 0!
      if (action == INSERT )
        insertNew = true;
      if( action == EX )
        EXCEPTION("element '" << parent << "' has no child '" << child
                  << "' with value '" << value << "'");
    }
    
    // Insert new node with given value
    // action is either APPEND or INSERT, i.e. the following unconstrained
    // Get()-calls will always create a child elements
    if( insertNew || action == APPEND ) {
      
      PtrParamNode ret =  this->Get(parent, action);
      ret->Get(child,action)->SetValue(value);
      return ret;
    }
    
    
    
    return PtrParamNode();
  }
  PtrParamNode ParamNode::GetByVal( const std::string& parent, 
                                   const std::string& child, 
                                   const char* value,
                                   const ActionType action) {
    return this->GetByVal( parent, child, std::string(value), action);
  }
    
  StdVector<PtrParamNode> ParamNode::GetList(const string& name)
  {
    const unsigned int chsize(children_.GetSize());
    StdVector<PtrParamNode> result;
    result.Reserve(chsize);

    for(unsigned int i = 0; i < chsize; ++i)
      if(children_[i]->name_ == name)
        result.Push_back(children_[i]);

    return result; // copy-constructor magic stuff!
  }
  
  template<typename TYPE>
  StdVector<PtrParamNode> ParamNode::GetListByVal(const string& parent, 
                                                  const string& child, 
                                                  const TYPE& value)
    {
      StdVector<PtrParamNode> result;
      
      for(unsigned int p = 0, chsize = children_.GetSize(); p < chsize; ++p)
      {
        PtrParamNode ch = children_[p];
        // do we have parent name? 
        if(ch->GetName() == parent)
        {
          // children_[i] has zero, one or more child matching child elements,
          // we loop by hand to handle the "more" case
          const unsigned int chchsize(ch->children_.GetSize());
          for(unsigned int c = 0; c < chchsize; ++c)
          {
            const PtrParamNode chch = ch->children_[c];
            if(chch->GetName() == child && chch->As<TYPE>() == value)
            {
              // parent, child and value matches. Take it!
              result.Push_back(ch);
            }
          }
        }
      }

      return result; // copy-constructor magic stuff!
    }

  StdVector<PtrParamNode> ParamNode::GetListByVal(const string& parent, 
                                                  const string& child, 
                                                  const char * value) {
    return GetListByVal( parent, child, std::string(value));
  }

  /************************************************************************
  * D A T A    G E T   M E T H O D S
  ************************************************************************/
  
  
  template<typename TYPE>
    TYPE ParamNode::As() const {
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
template<>\
  bool ParamNode::As<bool>() const {
  bool retVal = false;
  if( value_.type() == typeid(bool) ) {
       return boost::any_cast<bool>(value_);
     } else {
       if( value_.type() == typeid(std::string) ) {
         std::string str = boost::any_cast<std::string>(value_);
         if(str == "yes" || str == "true"  || str == "on")  retVal = true;
         if(str == "no"  || str == "false" || str == "off") retVal = false;
       } else {
         EXCEPTION("Could not convert node '" << name_ << "' to bool value");
       }
     }
  return retVal;
}

  template<typename TYPE>
  const TYPE& ParamNode::AsConst() const
  {

    // first, try to directly use any_cast<Integer>
    try
    {
      // first check, if type is same as returned value
      return boost::any_cast<const TYPE&>(value_);

    } catch (boost::bad_any_cast & ) {
      EXCEPTION("Cannot determine the data type of xml node '" << name_ << "'.");

    }
  }

  template<typename TYPE>
  void ParamNode::GetValue(const std::string& name, TYPE& ret, 
                           ActionType action) {

    if (action == DEFAULT )
      action = defaultAction_;

    if (action != APPEND ) {
      // First, check if node exists
      PtrParamNode node = this->Get(name, PASS );

      if(node) {
        // => node exists: convert value and return
        ret = node->As<TYPE>();
      } else {
        switch(action) {
          case PASS:
            return;
            break;
          case EX:
            EXCEPTION("element '" << name_ << "' has no child '" << name << "'");
            break;
          case INSERT:
          default:
            break;
        } // switch
      } // if node exists
    } // if APPEND

    // At this point we know, that we have to create a new entry, either
    // setting it uniquely (action == INSERT) or appending it (action==APPEND)
    PtrParamNode myNode =  this->Get(name,action);
    myNode->SetValue(ret);
  }
         
  

  /************************************************************************
  *  Q U E R Y     M E T H O D S
  ************************************************************************/
  bool ParamNode::Has(const string& name) const
  {
    // check in a fast way if we have tokens for trivial xpath
    if(ContainsTokens(name)) 
    {
      return TokenizedHasAndGet(name, string(""), false );
    }
    else
    {
      for(unsigned int i = 0, chsize = children_.GetSize(); i < chsize; i++)
        if(children_[i]->name_ == name) return true;

      return false;
    }
  }
  
  template<typename TYPE>
  bool ParamNode::HasByVal(const string& parent, 
                            const string& child, 
                            const TYPE& value ) const
    {
    if(ContainsTokens(parent)) {
      EXCEPTION("HasByVal(parent, child, value) does not allow for "
                 "tokenized search strings! ");
    }
    // see GetList() for comments
    for(unsigned int p = 0, chsize = children_.GetSize(); p < chsize; ++p)
    {
      PtrParamNode ch = children_[p];
      if(ch->name_ == parent)
      {
        const unsigned int chchsize(ch->children_.GetSize());
        for(unsigned int c = 0; c < chchsize; ++c)
        {
          PtrParamNode chch(ch->children_[c]);
          if(chch->GetName() == child && chch->As<TYPE>() == value)
            return true;
        }
      }
    }

    // nothing found
    return false;
    }
  
  bool ParamNode::HasByVal(const string& parent, 
                            const string& child, 
                            const char* value ) const {
    return HasByVal( parent, child, std::string(value));
  }
  
  template<typename TYPE>
    bool ParamNode::HasByVal(const string& name, const TYPE& value) const
    {
      if(ContainsTokens(name))
      {
        return TokenizedHasAndGet(name, value, true);
      }
      else
      {
        for(unsigned int i = 0, chsize = children_.GetSize(); i < chsize; i++)
          if(children_[i]->name_ == name && children_[i]->As<TYPE>() == value) 
            return true;

        return false;
      }
    }

  bool ParamNode::HasByVal(const string& name, const char* value) const
  {
    return HasByVal(name, std::string(name));
  }
  
  bool ParamNode::HasChildren() const {
    return (bool)children_.GetSize();
  }
  
  unsigned int ParamNode::Count(const string& name) const
  {
     unsigned int count = 0;

     for(unsigned int i = 0, chsize = children_.GetSize(); i < chsize; ++i)
       if(children_[i]->name_ == name) ++count;
    
     return count;
   }
/************************************************************************
* M I S C    M E T H O D S
************************************************************************/
  

template<typename TYPE>
PtrParamNode ParamNode::TokenizedHasAndGet( const string& name, 
                                            const TYPE& value,
                                            bool restrictToVal ) const
{
  StdVector<string> tokens = SplitIntoTokens(name);
  PtrParamNode ptr = boost::const_pointer_cast<ParamNode>(shared_from_this());
  for(unsigned int i = 0; i < tokens.GetSize(); i++)
  {

    // consider the value only for the last token!
    if(i == tokens.GetSize()-1)
    {
      // only look for specified value, if value is not empty
      // otherwise, we have already the searched for value
      if( restrictToVal ) { 
        if(!ptr->HasByVal(tokens[i], value)) 
          return PtrParamNode();
      } else {
          if(!ptr->Has(tokens[i])) return PtrParamNode();
      }
    }
    else
    {
      if(!ptr->Has(tokens[i])) return PtrParamNode();
    }
    ptr = ptr->Get(tokens[i]);
  }

  return ptr;
}


  void ParamNode::ToString(std::string& ret) const
  {
    // This method is currently very hardcoded.
    // In the future we rely on a general formatter class,
    // which can convert the special types 
    if( value_.type() == typeid(std::string) ) {
      ret =  boost::any_cast<std::string>(value_);
      return;
    }    
    if( value_.type() == typeid(double) ) {
         double val = boost::any_cast<double>(value_);
         std::stringstream stream;
         stream.precision(6);
         stream << val;
         ret =  stream.str();
         return;
    }
    if( value_.type() == typeid(UInt) ) {
      UInt val = boost::any_cast<UInt>(value_);
      ret =  boost::lexical_cast<std::string>(val);
      return;
    }
    if( value_.type() == typeid(Integer) ) {
      Integer val = boost::any_cast<Integer>(value_);
      ret =  boost::lexical_cast<std::string>(val);
      return;
    }
    if( value_.type() == typeid(bool) ) {
      bool val = boost::any_cast<bool>(value_);
      ret = val ? "yes" : "no";
      return;
    }
  
    // special conversion tpyes for vector
    if( value_.type() == typeid(Vector<Double>) ) {
      Vector<Double>  vec  = boost::any_cast<Vector<Double> >(value_);
      ret = vec.ToString();
      return;
    }
    if( value_.type() == typeid(Vector<Complex>) ) {
      Vector<Complex>  vec  = boost::any_cast<Vector<Complex> >(value_);
      ret = vec.ToString();
      return;
    }
    // special conversion types for matrix
    if( value_.type() == typeid(Matrix<Double>) ) {
     Matrix<Double>  mat  = boost::any_cast<Matrix<Double> >(value_);
      ret = mat.ToString();
      return;
    }
    if( value_.type() == typeid(Matrix<Complex>) ) {
      Matrix<Complex> mat  = boost::any_cast<Matrix<Complex> >(value_);
      ret = mat.ToString();
      return;
    }
    
    
    // special conversion types for vector
    if( value_.type() == typeid(Vector<Double>*) ) {
      Vector<Double> & vec  = *(boost::any_cast<Vector<Double> *>(value_));
      ret = vec.ToString();
      return;
    }
    if( value_.type() == typeid(Vector<Complex>*) ) {
      Vector<Complex> & vec  = *(boost::any_cast<Vector<Complex>* >(value_));
      ret = vec.ToString();
      return;
    }
    // special conversion types for matrix
    if( value_.type() == typeid(Matrix<Double>*) ) {
      Matrix<Double> & mat  = *(boost::any_cast<Matrix<Double>* >(value_));
      ret = mat.ToString();
      return;
    }
    if( value_.type() == typeid(Matrix<Complex>*) ) {
      Matrix<Complex> & mat  = *(boost::any_cast<Matrix<Complex>* >(value_));
      ret = mat.ToString();
      return;
    }
    
    // special conversion for timer
    if( value_.type() == typeid(Timer*) ) {
      Timer* timer = boost::any_cast<Timer* >(value_);
     ret = timer->ToXML();
     return;
    }

  }
  
  void ParamNode::ToXML( std::ostream& os , int depth ) {
    // note, that this is an recursive method!
    
    std::string strValue;
    ToString(strValue);
    
    if(type_ == ATTRIBUTE )
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
    os << std::endl << string(depth, ' ');
    if( type_ != COMMENT ) {
      //os << "<" << ToValidLabel(name_);
      os << "<" << name_;
    }

    // go through all children an check if we close with "../>" or </name_>
    // note, that the own XML-element value can also be a attribute
    bool only_attributes = true;  

    const UInt chsize(children_.GetSize());
    for(UInt i = 0; i < chsize; i++) {
      PtrParamNode actNode = children_[i];
      if (actNode->type_ != ATTRIBUTE)
        only_attributes = false;
      break;
    }

    // first loop, do the attributes
    for(UInt i = 0; i < chsize; i++) {
      PtrParamNode actNode = children_[i];
      if(actNode->type_ == ATTRIBUTE ) {
        children_[i]->ToXML(os);
      }
      else {
        only_attributes = false; // no, there is a non-attribute element
      }
    }

    // The following code should not be needed ..
    // the own value might be an attribute
    //      if(!value_.empty() && IsGoodAttribute(ToString()))
    //        os << " value=\"" << ToString() << "\"";


    if(only_attributes) {
      if(type_ == COMMENT) {
        os << "<!-- " << strValue << "-->" << std::endl;
      } else {
        // check, if value is set
        if (!value_.empty() ) {
          os <<"> " <<  strValue << "</" << name_ << ">" << std::endl;
        } else {
          os << "/>"; // "quick close"
        }
      }
    }
    else {
      os << ">";
      
      // second loop, child elements
      for(unsigned int i = 0; i < chsize; i++) {
        // attributes are already done
        if(children_[i]->type_ == ATTRIBUTE ) continue; 
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
      if(chsize != 0 && !endl_written) os << std::endl;
      os << string(depth, ' ');
      if( type_ != COMMENT) {
        //os  << "</" << ToValidLabel(name_) << ">";
        os  << "</" << name_ << ">";
      }
    }
  }
  
  void ParamNode::ToFile( const std::string& filename) {
    
    // determine correct filename (if not given)
    std::string myFileName;
    if(!filename.empty()) {
      myFileName = filename;
    } else {
      myFileName = progOpts->GetSimName()+".info.xml";
    }
    // write preamble

    std::ofstream info_file(myFileName.c_str());
    info_file << "<?xml version=\"1.0\"?>" << std::endl;
    // store how often we are written -> if the number is too high one should cancel some ToFile() calls
//    if(writeCounter_.count(filename_) == 0) writeCounter_[filename_] = 0;
//    Get("writeCounter")->SetValue(++writeCounter_[filename_]);
    
    // First we determine the type of nodes
    AdjustElementType();
    
    // Then we print the tree
    ToXML(info_file);
    info_file.close();
  }
       
  void ParamNode::Dump(int level) const
  {
    for(int i = 0; i < level; i++) cout << "   ";
    
    std::string strValue;
    ToString(strValue);
    // Format element name
    switch( type_ ) {
      case UNDEF:
        cout << name_ << ": " << strValue << std::endl;
        break;
      case ELEMENT:
        cout << "<" << name_ << "> "<<  strValue << std::endl;
        break;
      case ATTRIBUTE:
        cout << "@" << name_ << ": " <<  strValue << std::endl;
        break;
      case COMMENT:
        cout << "<!-- " << strValue << " -->" << std::endl;
      break;
    }
    for(unsigned int i = 0, chsize = children_.GetSize(); i < chsize; i++) 
      children_[i]->Dump(level + 1);
  }

  StdVector<string> ParamNode::SplitIntoTokens(const string& input) const
  {
    StdVector<string> out;
    char_separator<char> sep("/");
    tokenizer<char_separator<char> > tokens(input, sep);
    
    for(tokenizer<char_separator<char> >::iterator it = tokens.begin(); 
        it != tokens.end(); ++it)
      out.Push_back(*it);
    return out;
  }
  
  void ParamNode::AdjustElementType() {
    
    // loop over all children and call recursivley method
    // to determine type of nodes
    ParamNodeList::iterator it = children_.Begin();
    for( ; it != children_.End(); it++ ) {
      (*it)->AdjustElementType();
    }

    // Adjust element type
    if( type_ == UNDEF ) {
     

      // if children are present -> ELEMENT
      if( children_.GetSize() ) {
        type_ = ELEMENT;
      } else {
        // if no children are present
        //   value_ not set -> ELEMENT
        //   value_ set     -> ATTRIBUTE
        if( value_.empty() ) { 
          type_ = ELEMENT;
        } else {
          type_ = ATTRIBUTE;
        }
      }
    }
    // Check if node is ELEMENT, has children and a value
    // -> not possible in a XML tree
    if( type_ == ELEMENT && children_.GetSize() && !value_.empty() ) {
      WARN("Node '" << name_ << "' has children AND a non-empty value. " 
           << "This is not possible in a xml tree!");
    }

    // Check if node is ATTRIBUTE and has children
    // -> not possible in XML tree
    if( type_ == ATTRIBUTE && children_.GetSize() ) {
      WARN("Node '" << name_ << "' is a ATTRIBUTE and has child nodes. "
           << "This is not possible in a xml tree!");
    }
    
    
    // HARD-CODED Intermediate check for Timer element:
    // As the method "ToXML()" of the timer currently returns hard-coded 
    // the string representation of a xml element, we have to ensure that
    // the type of ParamNodes containing a Timer-object always ELEMENT is
    if( value_.type() == typeid(Timer*) ) {
      type_ = ELEMENT;
    }
    
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

 void ParamNode::ToValidLabel(const string& in, std::string& out) const
  {
    out = in;
    boost::trim(out);

    // remove spaces and don't make upper case yet! :(
    boost::erase_all(out, " ");
    boost::erase_all(out, "\"");
    boost::erase_all(out, "(");
    boost::erase_all(out, ")");
    // don't touch the slash '/', it is a 'xpath' element
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
  INSTANTIATE_METHOD_AS(Timer*)

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
  INSTANTIATE_METHOD_GETVALUE(Timer*)
  
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
    StdVector<PtrParamNode> ParamNode::GetListByVal<TYPE>\
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


