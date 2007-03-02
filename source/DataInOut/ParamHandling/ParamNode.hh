#ifndef PARAMNODE_HH_
#define PARAMNODE_HH_

#include <string>
#include "Utils/StdVector.hh"

namespace CoupledField
{
  /** The global pointer of the ParamNode (tree) holding the XML file */
  class ParamNode;
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
    ParamNode() {};

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
    Integer               AsInt() const;

    /** @return the UInteger if this is convertable
     * @throws an exception if the value is not set or not convertible or negative */
    UInt       AsUInt() const;

    /** @return the Double if this is convertable
     * @throws an exception if the value is not set or not converible */
    Double             AsDouble() const;
       
    /** @return checks the value for "yes", "true", "on" respectively "no", "false", "off".
     * @throws exception if no value or none of the above */
    bool               AsBool() const;
     
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

    void Get(const std::string& name, Integer& ret, const bool throwException = true);

    void Get(const std::string& name, UInt& ret, const bool throwException = true);

    void Get(const std::string& name, Double& ret, const bool throwException = true);

    void Get(const std::string& name, bool& ret, const bool throwException = true);
    //@}
    
    //@{
    /** gets the only direct child which has the name and stores its value in ret as Matrix<br>
     * Is only valid, if the corresponding GetList() would return one value. Exception if 0 or greater 1
     * ,if the found value can not be converted into type of return value or if the dimensions
     * of of the matrix provided do not match.<br>
     * If no matching element is found an exception will be thrown if throwException is set to true. 
     * Otherwise the function will return silently and the original value in ret will be kept.
     * example: "optimization" is a complex element which is a direct child of the root: param.Get("optimization") 
     * @throws exception if there is not such a direct child, e.g. if this is a leaf node OR if there more than only
     * one of such elements (e.g. simple xml elements). */     
    void GetDim1xDim2Tensor( const std::string& name, const unsigned int &dim1,
                             const unsigned int &dim2, Matrix<Integer>& ret,
                             const bool throwException = true );

    void GetDim1xDim2Tensor( const std::string& name, const unsigned int &dim1,
                             const unsigned int &dim2, Matrix<UInt>& ret,
                             const bool throwException = true );

    void GetDim1xDim2Tensor( const std::string& name, const unsigned int &dim1,
                             const unsigned int &dim2, Matrix<Double>& ret,
                             const bool throwException = true );
    //@}
    
    /** Checks if there is at least one direct child with the given name.<br>
     * Note, that even when not specified in the XML file, the value might come from the default value in the
     * XML schema definition.<br>
     * Does not differentiate between one or more than one occurence (then Get() will throw an exception). Note,
     * that there might be only one XML attribute but multiple XML simple elements and we do not differentiate 
     * between this two types. 
     * @return true if there is at least one direct child (leaf or "complex") with the given name */
    bool Has(const std::string& name) const;  

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
    UInt Count(const std::string&  name) const; 
       
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
                    const bool throwException = false );

    //@{
    /** Get the direct child which has an attribute with a given value and store it in the
     *  retrn variable or in other words:
     *  Get the direct child where the grandchildren are as specified.<br>
     *  example: param.Get("pdeList").Get("mechanic").Get("bcsAndLoads").Get("dirichletInHom", "name", "fixed") */
    void Get( const std::string& parent, const std::string& child, const std::string& value,
              std::string& ret,
              const bool throwException = false );

    void Get( const std::string& parent, const std::string& child, const std::string& value,
              Integer& ret,
              const bool throwException = false );

    void Get( const std::string& parent, const std::string& child, const std::string& value,
              UInt& ret,
              const bool throwException = false );

    void Get( const std::string& parent, const std::string& child, const std::string& value,
              Double& ret,
              const bool throwException = false );

    void Get( const std::string& parent, const std::string& child, const std::string& value,
              bool& ret,
              const bool throwException = false );
    //@}
    
    
    /** Returns the number of entries, the corresponding GetList() would return.
     * @See GetList(const std::string&, const std::string&, const std::string&) */
    UInt Count( const std::string& name, const std::string& child, const std::string& value ); 

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
  }; 

} // end of namespace


#endif /*PARAMNODE_HH_*/
