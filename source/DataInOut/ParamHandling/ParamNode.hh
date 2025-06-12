#ifndef ParamNode_HH_
#define ParamNode_HH_

#include <string>
#include <boost/enable_shared_from_this.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/tokenizer.hpp>
#include <boost/any.hpp>

#include "Utils/StdVector.hh"
#include "General/Exception.hh"

using std::string;

namespace CoupledField
{
  
  class Timer;
  class ParamNode;
  
  /** Definitions of pointers, using boost::shared_ptr */
  typedef boost::shared_ptr<ParamNode> PtrParamNode;
  typedef StdVector<boost::shared_ptr<ParamNode> > ParamNodeList;
  
  /** This class realizes the following concept of param handling, mainly the representation
   * of the XML file.
   * <ul>
   *   <li>The whole XML structure is contained in a tree of ParamNode elements</li>
   *   <li>Schema validation has to be done externally</li>
   *   <li>The construction, e.g. from a DOM tree has to be done externally</li>
   *   <li>An attribute is equal to a simple XML element content (<element type="1"/>
   *       is equal to <element><type>1</type></element>)</li>
   *   <li>We internally store everything as a string but offer conversion methods. 
   *       Only leaf nodes have a value set</li>
   *   <li>There are leaf nodes which are defined as having no children</li>
   *   <li>The children can be unsorted leaf nodes or "more complex" nodes</li>
   *   <li>You can alter the content (modify value or modify children)</li>
   * </ul> 
   * 
   * The ParamNode internally disinguishes the types of nodes (ATTRIBUTE, ELEMENT, COMMENT), 
   * this is just needed when performing i/o in xml format.
   * Things to note:
   *  - node type normally not of ineteres
   *  - const char* implicitly converted to std::string
   *  - 
   *
   * ToDO:
   * - Implement special case for template<const char*> to be converted into normal strings
   * - Merge missing methods from ParamNode
   *
   *
   * Some open questions remain:
   *  - Should we use shared_ptr<> instead of plain *? This would ease the use withing the
   *    program and prevent memory leaks, but it would introduce some addtional runtime 
   *    overhead
   *  - Should we clearly distinguish between set / get methods? At the moment, most get-
   *    methods implicitly set the value passed, if INSERT / ADD is given
   *  - Should all methods direclty return the PtrParamNode they are working on as return value?
   *    This would allow for chains in the form:
   *      myNode->SetValue("aValue")->Get("myChildNode")->SetValue("test");
   *  - Is it necessary to have a separate ParamNode? One possibility is to re-implement
   *    the methods of the ParamNode, but instead of the dfeault argument EX they have a default
   *    argument APPEND / REPLACE, i.e. by default we just add new elements, if they are not present
   *  - How to treat comments? Should we have a method SetComment(), wich just creates a comment node
   *    as child?
   *  
   * */  
  class ParamNode : public boost::enable_shared_from_this<ParamNode>
  {
  public:
    
    /** This string constant shall mark the logging part which contains "dynamic" data, e.g. current iteration */
    const static std::string HEADER;
    const static std::string PROCESS;
    /** this string constants shall mark summary information, e.g. total number of iterations */
    const static std::string SUMMARY;
    const static std::string WARNING;
    const static std::string FAIL;
    
    /** Define behavior for get() methods in case the element does not exist:
    * DEFAULT: Apply default action of current node (gets inherited)
    * EX:   Throw exception, if element does not exist
    * PASS: If element does not exist, just leave method and leave all variables
    *       in their original state. Otherwise return element.
    * INSERT: Create new element with default value if the element does not exist.
    *          If it exists, return value.
    * APPEND: Always create new node with default value.
    */
    
    typedef enum { DEFAULT, EX, PASS, INSERT, APPEND } ActionType;
    
    /** Possible node types
     * UNDEF: In most cases we do not care about the type of the paramnode,
     *        so UNDEF is a valid state. Only when converting a node to an XML
     *        structure, we have to convert it to another type
     * ELEMENT: This node represents a XML element
     * ATTRIBUTE: This node represents a XML attribute
     * COMMENT: This node represents just a comment without syntactical meaning.
     *          Note: Comment-nodes are not allowed to have children!
     * SELF_XML: Complex elements like Timer or Matrix which create an
     *           xml string by themselves.
     * BULK:    A pre-formated block for fast large files (density.xml in optimization)
     */
    typedef enum { UNDEF, ELEMENT, ATTRIBUTE, COMMENT, SELF_XML, BULK } NodeType;
    
    /** The default constructor.
     * Can also be used to generate a root node. For a root node of a to be written structure (info.xml) see GenerateWriteNode())
     * @param type Type of the node (defaults to UNDEF)
     * @return */
    ParamNode(ActionType defaultAction = EX, NodeType type = UNDEF);

    /** Recursively delete the child nodes */
    virtual ~ParamNode();

    /** Generates the root node for a tree to be written, e.g. info.xml
     * @param root the name of the root element
     * @param filename if empty any ToFile() will return silently
     * @param lazy_write will in non debug only write when the last ToFile() was sufficiently long ago
     * @param add_counter adds write and reject counter to the root element */
    static PtrParamNode GenerateWriteNode(const string& root, const string& filename, ActionType defaultAction = INSERT, bool lazy_write = false, bool add_counter = false);

    /************************************************************************
    * S E T    M E T H O D S
    *************************************************************************/
    
    /** Set name of the current node */
    void SetName(const std::string& name) { this->name_ = name; }

    /** Set the type of the current node */
    void SetType(const NodeType type) { this->type_ = type; }
    
    /** Set a valid value (native type, Timer* and Matrix and Vector instances and pointers.
     * Any other type which has no special SetValue() implementation is lost!! */
    void SetValue(const boost::any& value, bool cerr_warning = true);

    /** a string is a string :) */
    void SetValue(const char* value);

    /** Special version which handles the precision of the value. See implementation note! */
    void SetValue(const double, const int precision);

    /** Set a ParamNode to an expandable InfoNode. Works recursively
     * @param overwrite_name if the the name of the parent node should be used
     * @param cerr_warning shall a warning be written to cerr */
    void SetValue(PtrParamNode node, bool overwrite_name = false, bool cerr_warning = true);

    /** Set ParamaNodes as child nodes. Works recursively */
    void SetValue(ParamNodeList nodes);

    /** Creates a sub-node with the content */
    void SetComment(const std::string& string);

    /** Creates a ParamNode::WARNING node and sets the msg as content. */
    void SetWarning(const std::string& msg, bool append = false);

    /** Returns the root node. */
    PtrParamNode GetRoot();

    /** Returns father element. If this element is the root node, NULL is returned */
    PtrParamNode  GetParent() { return parent_;}
    
    /** @return the name of the attribute or XML element */
    const std::string& GetName() const { return name_;} 

    /** return the path up to the current element */
    std::string GetLocation() const;

    /** Add child parameter nod*/
    void AddChildNode( PtrParamNode child);
    
    /** Returns all children which are attributes and simple xml elements (cannot be differentiated) or in
    * other words leaf nodes - and without any sorting complex ParamNodes which have children themselves. 
    * If this element is already a leaf node, the list is empty.<br>
    * Be careful when editing this list! */
    ParamNodeList& GetChildren() { return children_;}
    
    /** If you really know what you do, you can set a child manually.
      * An existing value is overwritten and not deleted!
      * @param name an element with this name will exist afterwards
      * @param index the children list must be large enough.
      * @return the newly created object. */
    PtrParamNode SetNewChild(const std::string& name, unsigned int index);
    
    PtrParamNode ReplaceChild(PtrParamNode node, unsigned int index);

    /** Convenience function to clear all children. Here smart pointers are really nice */
    void ClearChildren();

    /** clear all children by given name. This sorts all entries if something is cleared */
    void ClearChildren(const string& name);

    /** Returns the only child of an element which might be an attribute or simple xml element or in
    * other words a leaf node - and without any sorting a complex ParamNode which has children by itself. 
    * If the element has more than one child node, an exception is thrown.
    * If this element is already a leaf node, the return value is NULL. */
    PtrParamNode GetChild();

    /** returns the only direct child which has the name.<br>
    * Is only valid, if the corresponding GetList() would return one value. Exception if 0 or greater 1.<br>
    * In case check before with Has() and Count().
    * example: "optimization" is a complex element which is a direct child of the root: param.Get("optimization")
    * @param name might contain several levels by the '/' separator. 
    *             Get("optimization/ersatzMaterial") is equivalent to Get("optimization")->Get("ersatzMaterial").
    *             Check with Has() first! 
    * @throws exception if there is not such a direct child, e.g. if this is a leaf node OR if there more than only
    * one of such elements (e.g. simple xml elements). */     
    PtrParamNode Get(const std::string&  name, ActionType = DEFAULT ); 
    
    //@{
    /** Get the direct childs which has an attribute with a given value or in other words:
    *  Get the direct child where the grandchildren are as specified.<br>
    *  example: param->Get("pdeList/mechanic/bcsAndLoads/dirichletInHom", "name", "fixed")
    *  will return that dirichletInHom where @name == fixed.
    *  @return the direct child, not the grandchildren */
    template<typename TYPE >
    PtrParamNode GetByVal( const std::string& parent, 
                          const std::string& child, 
                          const TYPE& value,
                          ActionType action = DEFAULT );
    PtrParamNode GetByVal( const std::string& parent, 
                          const std::string& child, 
                          const char* value,
                          ActionType action = DEFAULT );
    //@}

    /** Get an element identified by two childs with prescribed value. Extend to the full standard GetByVal() features if you need them :) */
    PtrParamNode GetByVal(const std::string& parent_raw, const std::string& child1,  const std::string& value1,
                                                         const std::string& child2,  const std::string& value2, ActionType action = DEFAULT);

    /** Get all direct childs of a name
    * example: param.Get("pdeList").Get("mechanic").Get("bcsAndLoads").GetList("dirichletInHom") */
    ParamNodeList GetList(const std::string&  name);
    
    /** search for the given token in the names of children */
    static ParamNodeList GetListByChild(const ParamNodeList& base, const std::string&  name);
    /** search for the given token in the names of children of base-children (= base grandchilds) */
    static ParamNodeList GetListByGrandChild(const ParamNodeList& base, const std::string&  name);


    //@{
    /** Get all direct childs which have an attribute with a given value or in other words:
    *  Get all direct childs where the grandchildren are specified.<br>
    *  example: param.Get("pdeList").Get("mechanic").Get("bcsAndLoads").GetList("dirichletInHom", "name", "fixed")
    *  @return the direct childs, not the grandchildren */
    template<typename TYPE>
    ParamNodeList GetListByVal(const std::string& parent, const std::string& child, const TYPE & value);

    ParamNodeList GetListByVal(const std::string& parent, const std::string& child, const char* value);
    //@}

    /************************************************************************
    * D A T A    G E T   M E T H O D S
    ************************************************************************/

    /** by "fast bulk block" or better "dirty bulk block" we denote optional content as xml formated strings.
     * For high load applications as grid export (-G), streaming (grid and results and optimization design export
     * this can save a significant amount of time.
     * Logically the strings are written (ToXML() after param node childs but not considered in Getxx(). So it is
     * rather dirty thand fast ! :(.
     * The fast bulk block, a StdVector of std::string is stored as a boost::any type.
     * @return create the block if it does not exist already. Exception if boost::any has already another type. */
    StdVector<std::string>& GetFastBulkBlock();

    /** @return converted value
    * @throws an exception if the value is not set or not convertible
    * version for all types instead of integral ones */
    template<typename TYPE>
    TYPE As() const;

    /** This is a special case where the timer is created if it does not already exist.
     * Therefore no constness */
    boost::shared_ptr<Timer> AsTimer();

    /** @return the integer if this is convertible
    * @throws an exception if the value is not set or not convertible */
    template<typename TYPE>
    const TYPE& AsConst() const;
    
    /** @return the value of the node, parsed as number by mathParser **/
    template<typename TYPE>
    TYPE MathParse() const;

    /** Directly access the value of the current node as a given type.
     *
     * The most prominent practical usage is GetValue("hans", value_, ParamNode::PASS).
     * If "hans" is given in xml, the value is read and written to value_. If "hans" is not
     * present, value_ is not touched (and the C++ default remains). Makes only sense,
     * if there is no schema default. For a schema default you can use
     * value_ = Get("hans")->As<TYPE>();
     *
     *  In case the node does not exist the following action is taken:
     *  PASS: Method returns without warning/error if the element is not contained and does not touch reg
     *  EX: An exception is generated
     *  INSERT: A new node is generated with the value of the ret variable as new value
     *  APPEND: A new node is generated anyway and the value of the ret variable  is set as new value */
    template<typename TYPE>
    void GetValue(const std::string& name, TYPE& ret, ActionType = EX );        
        
    /** This Get() version overwrites with throwException=false a preset value and does nothing if the
    *  value does not exist. So you do not need to check with Has() first.<br>
    *  Get the direct child which has an attribute with a given value and store it in the
    *  return variable or in other words:
    *  Get the direct child where the grandchildren are as specified.<br>
    *  example: 
    *  <pre>
    // instead of:
    double manual_scaling = node->Has("option", "name", "obj_scaling_factor") ?
    node->Get("option", "name", "obj_scaling_factor")->Get("value")->As<Double>() : 1.0;
    // you can do
    double manual_scaling = 1.0;
    node->Get<double>("option", "name", "obj_scaling_factor", manual_scaling, "value", false);
    </pre> 
//    */
//    template <typename TYPE, typename TYPE2>
//    void Get( const std::string& parent, 
//              const std::string& child, 
//              const std::string& value, 
//              TYPE& ret, 
//              const TYPE2& value_node, 
//              const ActionType = EX );
    
    /************************************************************************
    *  Q U E R Y     M E T H O D S
    ************************************************************************/

    /** Checks if there is at least one direct child with the given name.<br>
    * Note, that even when not specified in the XML file, the value might come from the default value in the
    * XML schema definition.<br>
    * Does not differentiate between one or more than one occurrence (then Get() will throw an exception). Note,
    * that there might be only one XML attribute but multiple XML simple elements and we do not differentiate 
    * between this two types. 
    * @param name might contains several levels by the '/' token.
    * @return index of the first direct child (leaf or "complex") with the given name or -1 if none */
    int GetIndex(const std::string& name) const;

    /* @seee GetIndex()
    * @return true if there is at least one direct child (leaf or "complex") with the given name */
    bool Has(const std::string& name) const;


    //@{
    /** Checks if there is at least one direct child with the given name and an attribute with the given value.<br>
    * Note, that even when not specified in the XML file, the value might come from the default value in the
    * XML schema definition.<br>
    * Does not differentiate between one or more than one occurrence (then Get() will throw an exception). Note,
    * that there might be only one XML attribute but multiple XML simple elements and we do not differentiate 
    * between this two types. 
    * @return true if there is at least one direct child (leaf or "complex") with the given name */
    template<typename TYPE>
    bool HasByVal( const std::string& name,
                   const std::string& child, 
                   const TYPE& value ) const;  
    bool HasByVal( const std::string& name,
                   const std::string& child, 
                   const char* value ) const;
    //@}

    //@{
    /** Checks if a direct child (e.g. attribute) exists and has a special value
    * @see Has(const std::string&) const */
    template<typename TYPE>
    bool HasByVal(const std::string& name, const TYPE& value) const;
    bool HasByVal(const std::string& name, const char* value) const;
    //@}
    
    /** Checks if current node has children */
    bool HasChildren() const;

    /** Returns the number of entries, the corresponding GetList() would return.
    * @See GetList(const std::string&) */
    unsigned int Count(const std::string&  name) const; 

    /************************************************************************
    * M I S C    M E T H O D S
    ************************************************************************/

    /** returns name and value, and child summary information
     * @param depth if the elemen is matrix type the depth is mandatory for Matrix::ToXMLFormat(name, depth) */
    void ToString(std::string& ret, int depth ) const;
    
    /** variant of other ToString() with more copy operations but also more convenient */
    std::string ToString(int depth) const;

    /** Prints this as xml element to the stream. Builds a tree. Shall not be directly
     * called for an attribute.
     * Might change the order as Sort() is called to ensure HEADER < PROCESS < SUMMARY
     * @param depth number of ident steps, will be interpreted as one space
    * @param adjust_element_type this is done in ToFile(), here only for the stream output writer */
    void ToXML(std::ostream& os, int depth = 0, bool adjust_element_type = false);
    
    /** Print node tree to file. If the name ends with .gz it will be automatically compressed!
     * Normally the tree will be created by GenerateWriteNode() where a name will be given. But you can write any element
     * @param filename if not given an created by GenerateWriteNode() we will use the name from there (no name=silent no write)
     * @param force if not debug and lacy write this to File might ignored
     * Note: This is just a temporary solution, until we move the serialization of the 
     * ParamNode to the WriteInfo class  */
    void ToFile(const std::string& filename = std::string(), bool force = false);

    void ToFile(bool force) { ToFile(std::string(), force); }

    /** This is a recursive Dump of the tree to std::cout
    * @param level start with 0, is used for ident */
    void Dump(int level = 0) const;

    /** Similar to recursive dump, this adds all childs to the flat list.
     * colons (:) separates (sub) child names. */
    StdVector<std::pair<std::string, std::string>> ToStringList(int max_level = 9999) const;

    /** The ParamNode::Dump() shows all sub-content. This shows the parent path */
    void DumpParentPath();
    
    /** Helper method that splits a string by slashes '/'. Is trivial xpath stuff
    * @param input might be empty and containt slashes with strings inbetween
    * @return an empty list only for an empty parameter */
    StdVector<std::string> SplitIntoTokens(const std::string& input) const;

    /** Speed up the tokenizer stuff from SplitIntoTokens() by asking if it contains tokens. */
    bool ContainsTokens(const std::string& input) const
    {
      return input.find('/') != std::string::npos;
    }

  private:

    /** Helper that implements a Get() as same as Has(). Call this if ContainsTokens() returns true. */
    template<typename TYPE>
    PtrParamNode TokenizedHasAndGet(const std::string& name, const TYPE& value, bool restrictToVal ) const;
    
    /** Makes a valid XML name, e.g. removes spaces. Used to create
    * automatically an ParamNode::name_ value out of caption_
    * @param in might contain spaces, e.g. "Number of iterations"
    * @return for this example 'NumberOfIterations' */
    static std::string ToValidLabel(std::string in);

    /** Determine recursively the suitable type for nodes. This method
     * iterates recursively over all nodes and determines for all nodes, which 
     * have a type UNDEF a suitable type (ELEMENT, ATTRIBUTE) */
    void AdjustElementType();

    /** recursive helper */
    void ToStringList(StdVector<std::pair<std::string, std::string> >& list, int max_level, int level) const;

    /** The real content (attribute or simple type content) */
    boost::any value_;

    /** the precision for numerical values, to be optionally set! */
    int precision_;

    /** The name of the xml element/ attribute */
    std::string name_;
    
    /** Type of the node (ELEMENT, ATTRIBUTE, COMMENT)*/
    NodeType type_;
       
    /** This are the children of this element, either simple types (simple element
     * or attributes) or complex elements. If this element is a simple type the 
     * vector is empty. */
    ParamNodeList children_;
    
    /** pointer to father node */
    PtrParamNode parent_;
    
    /** default action for non-existing nodes */
    ActionType defaultAction_;
    
    /** cache variable for last Get() result */
    int lastresultidx_;

    /** Some attributes are only necessary for the root node of a writing ParamNode.
     * Collect them here to save footprint for all nodes. Easier than class hierarchy */
    struct RootNode
    {
      RootNode();
      ~RootNode();

      /** if ToString() provides another name this one is not overwritten but just temporarily ignored.
       * Nothing is written if filename is empty. */
      std::string filename;

      /** lazy write means that we only write for ToFile() w/o name and w/o force when the last write
       * was too long ago and when we are not in debug */
      bool lazy_write = false;

      /** shall we add the write_counter and recejt_counter information to the root element as attributes */
      bool add_counters = true;

      /** write_timer restricts the number of times the info.xml file is written
       *  if not enough time has passed, the file is not written
       *  this only affects ParamNode::ToFile()
       *  writing can be force via new parameter of ToFile */
      Timer* write_timer = NULL;

      /** how often the file is actually written */
      unsigned int write_counter = 0;
      /** how often we rejected writing the info.xml-file */
      unsigned int reject_counter = 0;
    };
    
    /** only a root node which is meant to be written with ToFile() needs it.
     * @see GenerateWriteNode() */
    RootNode* rootNode = NULL;
  }; 


} // end of namespace


#endif /*ParamNode_HH_*/
