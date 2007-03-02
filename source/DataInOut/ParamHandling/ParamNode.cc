#include "DataInOut/ParamHandling/ParamNode.hh"
#include "General/exception.hh"
#include "Utils/tools.hh"

namespace CoupledField
{
    
  /** This is our global pointer of the root ParamNode holding the XML file. 
   *  Filed in cfs.cc. The correpsonding
   * "extern ParamNode* param;" is in ParamNode.hh */
  ParamNode* param;     
    
    
  ParamNode::~ParamNode()
  {
    // values are no problem, but we have to delete our childs recursively
    for(UInt i = 0; i < children_.GetSize(); i++)
      {
        if(children_[i] != NULL) 
          { 
            delete children_[i]; 
            children_[i] = NULL; 
          }
      }
  }
  
  Integer ParamNode::AsInt() const
  {
    return String2Int(value_);
  }
  
  UInt ParamNode::AsUInt() const
  {
    return String2UInt(value_);    
  }

  Double ParamNode::AsDouble() const
  {
    return String2Double(value_);    
  }
       
  bool ParamNode::AsBool() const
  {
    if(value_ == "yes" || value_ == "true"  || value_ == "on")  return true;
    if(value_ == "no"  || value_ == "false" || value_ == "off") return false;
    
    EXCEPTION("cannot interpret " << value_ << " as boolean");
  }
     
  ParamNode* ParamNode::Get(const std::string& name, const bool throwException)
  {
    for(UInt i = 0; i < children_.GetSize(); i++)
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
   

  void ParamNode::Get(const std::string&  name, Integer& ret, const bool throwException)
  {
    ParamNode * actNode = Get(name, throwException);
    if( actNode ) {
      ret = actNode->AsInt();
    }
    
  }

  void ParamNode::Get(const std::string&  name, UInt& ret, const bool throwException)
  {
    ParamNode * actNode = Get(name, throwException);
    if( actNode ) {
      ret = actNode->AsUInt();
    }

  }

  void ParamNode::Get(const std::string&  name, Double& ret, const bool throwException)
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
 
  void ParamNode::GetDim1xDim2Tensor( const std::string& name, 
                                      const unsigned int &dim1,
                                      const unsigned int &dim2, 
                                      Matrix<Integer>& ret,
                                      const bool throwException )
  {
    // Fetch entry as string
    ParamNode * actNode = Get(name, throwException);

    if( actNode ) {
      StdVector<std::string> strVec;
      SplitStringList( actNode->AsString(), strVec, ' ' );
      
      Matrix<Integer> helpMat;
      helpMat.Resize( dim1, dim2 );
      helpMat.Init();
      
      if (strVec.GetSize() == dim1*dim2) {
        for ( UInt i = 0; i < dim1; i++ ) {
          for ( UInt j = 0; j < dim2; j++ ) {
            helpMat[i][j]=( String2Int( strVec[i*dim2+j] ) );
          }
        }
        ret = helpMat;
      } else{
        if( throwException)
          EXCEPTION("Wrong size of matrix '" << name 
                    << "'!. It contains of " << strVec.GetSize() 
                    << " entries which cannot be converted into matrix of size") ;
      }
    }
  }
  
  void ParamNode::GetDim1xDim2Tensor( const std::string& name, 
                                      const unsigned int &dim1,
                                      const unsigned int &dim2, 
                                      Matrix<UInt>& ret,
                                      const bool throwException )
  {
    EXCEPTION(" Not imeplemented" );
  }
  
  void ParamNode::GetDim1xDim2Tensor( const std::string& name, 
                                      const unsigned int &dim1,
                                      const unsigned int &dim2, 
                                      Matrix<Double>& ret,
                                      const bool throwException ) 
  {
    EXCEPTION(" Not imeplemented" );
  }
  

  bool ParamNode::Has(const std::string& name) const
  {
    for(UInt i = 0; i < children_.GetSize(); i++)
      if(children_[i]->name_ == name) return true;

    return false;
  }

     
  StdVector<ParamNode*> ParamNode::GetList(const std::string& name)
  {
    StdVector<ParamNode*> result;

    for(UInt i = 0; i < children_.GetSize(); i++)
      if(children_[i]->name_ == name) result.Push_back(children_[i]);
   
    return result; // copy-constructor magic stuff!
  }

  UInt ParamNode::Count(const std::string& name) const
  {
    UInt count = 0;

    for(UInt i = 0; i < children_.GetSize(); i++)
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

    for(UInt p = 0; p < children_.GetSize(); p++)
      {
        // do we have parent name? 
        if(children_[p]->name_ == parent) 
          {
            // children_[i] has zero, one or more child matching child elements, 
            // we loop by hand to handle the "more" case
            for(UInt c = 0; c < children_[p]->children_.GetSize(); c++)
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
    
    for(UInt p = 0; p < children_.GetSize(); p++)
      {
        // do we have parent name? 
        if(children_[p]->name_ == parent) 
          {
            // children_[i] has zero, one or more child matching child elements, 
            // we loop by hand to handle the "more" case
            for(UInt c = 0; c < children_[p]->children_.GetSize(); c++)
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
    for(UInt p = 0; p < children_.GetSize(); p++)
      {
        // do we have parent name? 
        if(children_[p]->name_ == parent) 
          {
            // children_[i] has zero, one or more child matching child elements, 
            // we loop by hand to handle the "more" case
            for(UInt c = 0; c < children_[p]->children_.GetSize(); c++)
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
  

  UInt ParamNode::Count(const std::string& parent, const std::string& child, const std::string& value)  
  {
    // better slow than copy-paste
    StdVector<ParamNode*> result = GetList(parent, child, value);
    return result.GetSize();
  }


  std::string ParamNode::ToString() const
  {
    std::ostringstream os;
    os << name_ << " = '" << value_ << "'";
    return os.str();
  }
       
  void ParamNode::Dump(int level) const
  {
    for(int i = 0; i < level; i++) std::cout << "   ";
    std::cout << ToString() << std::endl;
    
    for(UInt i = 0; i < children_.GetSize(); i++) 
      children_[i]->Dump(level + 1);
  }


} // end of namespace


