// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef PARAMNODE_HH_
#define PARAMNODE_HH_

#include <string>
#include "Utils/StdVector.hh"
#include "boost/lexical_cast.hpp"
#include "General/exception.hh"

namespace CoupledField
{
  /** The global pointer of the ParamNode (tree) holding the XML file */
  class ParamNode;
  template <class TYPE> class Matrix;  
  extern ParamNode* param;

  /** This class realizes the following concept of param handling, mainly the representation
   * of the XML file.
   * <ul>
   *   <li>The whole XML structure is contained in a tree of ParamNode elements</li>
   *   <li>Schema validation has to be done externaly</li>
   *   <li>The construction, e.g. from a DOM tree has to be done externaly</li>
   *   <li>An attribute is equal to a simple XML element content (<element type="1"/>
   *       is equal to <element><type>1</type></element>)</li>
   *   <li>We internally store everything as a string but offer conversion methods. 
   *       Only leaf nodes have a value set</li>
   *   <li>There are leaf nodes which are defined as having no children</li>
   *   <li>The children can be unsorted leaf nodes or "more complex" nodes</li>
   *   <li>You can alter the content (modify value or modify children)</li>
   * </ul> */  
  class ParamNode
  {
  public:
    /** The default constructor, name and value are to be set via the setter methods */
    ParamNode(bool attribute = false); 

    /** Recursively delete the child nodes */
    ~ParamNode();

    /** @return the name of the attribute or XML element */
    const std::string& GetName() const { return name_;} 
       
    /** Only when we fill the global root element from xerces we need this method */
    void SetName(const std::string& name) { this->name_ = name; }
       
    /** @return a string, a int/Double as string or empty if this is not a leaf node (we have children) */
    const std::string& AsString() const { return value_;} 

    /** set the value */
    void SetValue(const std::string& value) { this->value_ = value; }

    /** @return the integer if this is convertable
     * @throws an exception if the value is not set or not converible */
    int  AsInt() const;

    /** @return the UInteger if this is convertable
     * @throws an exception if the value is not set or not convertible or negative */
    unsigned int AsUInt() const;

    /** @return the Double if this is convertable
     * @throws an exception if the value is not set or not converible */
    double AsDouble() const;
       
    /** @return checks the value for "yes", "true", "on" respectively "no", "false", "off".
     * @throws exception if no value or none of the above */
    bool AsBool() const;
     
    /** returns the only direct child which has the name.<br>
     * Is only valid, if the corresponding GetList() would return one value. Exception if 0 or greater 1.<br>
     * In case check before with Has() and Count().
     * example: "optimization" is a complex element which is a direct child of the root: param.Get("optimization") 
     * @throws exception if there is not such a direct child, e.g. if this is a leaf node OR if there more than only
     * one of such elements (e.g. simple xml elements). */     
    ParamNode* Get(const std::string&  name, const bool throwException = true ); 
    

    //@{
    /** gets the only direct child which has the name and stores its value in ret <br>
     * Is only valid, if the corresponding GetList() would return one value. Exception if 0 or greater 1
     * or if the found value can not be converted into type of return value.<br>
     * If no matching element is found an exception will be thrown if throwException is set to true. 
     * Otherwise the function will return silently and the original value in ret will be kept.
     * example: "optimization" is a complex element which is a direct child of the root: param.Get("optimization") 
     * @throws exception if there is not such a direct child, e.g. if this is a leaf node OR if there more than only
     * one of such elements (e.g. simple xml elements). */
    void Get(const std::string& name, std::string& ret, const bool throwException = true);

    void Get(const std::string& name, int& ret, const bool throwException = true);

    void Get(const std::string& name, unsigned int& ret, const bool throwException = true);

    void Get(const std::string& name, double& ret, const bool throwException = true);

    void Get(const std::string& name, bool& ret, const bool throwException = true);

    
    /** Checks if there is at least one direct child with the given name.<br>
     * Note, that even when not specified in the XML file, the value might come from the default value in the
     * XML schema definition.<br>
     * Does not differentiate between one or more than one occurence (then Get() will throw an exception). Note,
     * that there might be only one XML attribute but multiple XML simple elements and we do not differentiate 
     * between this two types. 
     * @return true if there is at least one direct child (leaf or "complex") with the given name */
    bool Has(const std::string& name) const;  

    /** Checks if a direct child (e.g. attribute) exists and has a special value
     * @see Has(const std::string&) const */
    bool Has(const std::string& name, const std::string& value) const;

    /** Checks if a direct child (e.g. attribute) exists and is set to "true"/"false", "yes"/"no", "on"/"off"
    * @see Has(const std::string&) const */
    bool Has(const std::string& name, bool value) const;
    
    /** Checks if there is at least one direct child with the given name and an attribute with the given value.<br>
     * Note, that even when not specified in the XML file, the value might come from the default value in the
     * XML schema definition.<br>
     * Does not differentiate between one or more than one occurence (then Get() will throw an exception). Note,
     * that there might be only one XML attribute but multiple XML simple elements and we do not differentiate 
     * between this two types. 
     * @return true if there is at least one direct child (leaf or "complex") with the given name */
    bool Has(const std::string& name,
             const std::string& child, 
             const std::string& value ) const;  
        
    /** Get all direct childs of a name
     * example: param.Get("pdeList").Get("mechanic").Get("bcsAndLoads").GetList("dirichletInHom") */
    StdVector<ParamNode*> GetList(const std::string&  name);

    /** Returns the number of entries, the corresponding GetList() would return.
     * @See GetList(const std::string&) */
    unsigned int Count(const std::string&  name) const; 
       
    /** Get all direct childs which have an attribute with a given value or in other words:
     *  Get all direct childs where the grandchildren are es specified.<br>
     *  example: param.Get("pdeList").Get("mechanic").Get("bcsAndLoads").GetList("dirichletInHom", "name", "fixed")
     *  @return the direct childs, not the grandchildren */
    StdVector<ParamNode*> GetList(const std::string& parent, const std::string& child, const std::string& value);

    /** Get the direct childs which has an attribute with a given value or in other words:
     *  Get the direct child where the grandchildren are as specified.<br>
     *  example: param.Get("pdeList").Get("mechanic").Get("bcsAndLoads").Get("dirichletInHom", "name", "fixed")
     *  @return the direct child, not the grandchildren */
    ParamNode* Get( const std::string& parent, 
                    const std::string& child, 
                    const std::string& value,
                    const bool throwException = true );


    /** This Get() version overwrites with throwException=false a preset value and does nothing if the
     *  value does not exist. So you do not need to check with Has() first.<br>
     *  Get the direct child which has an attribute with a given value and store it in the
     *  return variable or in other words:
     *  Get the direct child where the grandchildren are as specified.<br>
     *  example: 
     *  <pre>
        // instead of:
        double manual_scaling = node->Has("option", "name", "obj_scaling_factor") ?
               node->Get("option", "name", "obj_scaling_factor")->Get("value")->AsDouble() : 1.0;
        // you can do
        double manual_scaling = 1.0;
        node->Get<double>("option", "name", "obj_scaling_factor", manual_scaling, "value", false);
        </pre> 
     */
    template <class TYPE>
    void Get(const std::string& parent, const std::string& child, const std::string& value, 
             TYPE& ret, const std::string& value_node, bool throwException = true);

    /** Returns the number of entries, the corresponding GetList() would return.
     * @See GetList(const std::string&, const std::string&, const std::string&) */
    unsigned int Count( const std::string& name, const std::string& child, const std::string& value ) const; 

    /** Returns all children which are attributes and simple xml elements (cannot be differentiated) or in
     * other words leaf nodes - and without any sorting complex ParamNodes which habe children by themself. 
     * If this element is already a leaf node, the list is empty.<br>
     * Be careful when editing this list! */
    StdVector<ParamNode*>& GetChildren() { return children_;}
    
    /** Returns the only child of an element which might be an attribute or simple xml element or in
     * other words a leaf node - and without any sorting a complex ParamNode which has children by itself. 
     * If the element has more than one child node, an exception is thrown.
     * If this element is already a leaf node, the return value is NULL. */
    ParamNode* GetChild();
    
    /** returns name and value, ans child summary information */
    std::string ToString() const;
         
    /** Prints this as xml element to the stream. Builds a tree. Shall no be directly
     * called for an attribute */
    void ToXML(std::ostream& os) const;      
         
    /** This is a recursive Dump of the tree to std::cout
     * @param level start with 0, is used for ident */
    void Dump(int level = 0) const;


  private:

    /** The real content (attribute or simple type content */
    std::string value_;
    
    /** The name of the xml element/ attribute */
    std::string name_;
       
    /** This are the children of this element, either simple types (simple element
     * or attributes) or complex elements. If this element is a simple type the 
     * vector is empty. */
    StdVector<ParamNode*> children_;
    
    /** Was this originally an attribute or no ? */
    bool attribute_;  
  }; 

  
  template <class TYPE>
  void ParamNode::Get(const std::string& parent, const std::string& child, const std::string& value, 
           TYPE& ret, const std::string& value_node, bool throwException)
  {
    ParamNode * base = Get(parent, child, value, throwException);
    if(base != NULL) 
    {
      try
      {
        ret = boost::lexical_cast<TYPE>(base->Get(value_node)->AsString());
      }
      catch(boost::bad_lexical_cast &)
      {
        EXCEPTION("cannot cast value '" << value_ << "' in XML node " << parent << "/" 
                  << child << "=" << value);
      }
    }// else do nothing
  }
  
  
} // end of namespace


#endif /*PARAMNODE_HH_*/
