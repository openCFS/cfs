#ifndef ENUM_HH_
#define ENUM_HH_

#include <map>
#include <string>

#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField 
{

/** There is a case where for a key we have multiple names. Then the
 * map key is not an enum but an large negative number which cannot be
 * retrieved. Therfore we store the real enum key in the value part.
 * As the Enum is Templated we cannot define this struct inside Enum<T> :( */
struct EnumTupel
{
  std::string string;
  int key; // this is always an enum key!
};

/** <p>This class is a convenient wrapper for enums, it is somehow oriented to the 
   * Java Enum type. It can also be used for general text string handling.</p>
   * <p>The basic idea for an <b>class internal</b> enum is as following:</p>
   * <pre>
   * class Foo 
   * {
   *   public:
   *     Foo() {
   *       SetEnums();
   *       filter_ = NONE; 
   *     }
   *   private:
   *     void foo(); 
   *     typedef enum { NONE, RADIUS, VOLUME_RADIUS } Filter;
   *     Enum<Filter> filter;
   *     void SetEnums() {
   *        filter.SetName("filter");
   *        filter.Add(NONE, "none");
   *        filter.Add(RADIUS, "radius");
   *        filter.Add(VOLUME_RADIUS);
   *     }
   *     Filter filter_; 
   * };
   * 
   * void Foo::foo() {
   *    filter_ = filter.Parse("radius");
   *    filter_ = filter.Parse("param->Get("radius"));
   *    bool ok = filter.isValid("none");
   *    bool ko = filter.isValid((Radius) 5);
   *    EXCEPTION("not handled: " << filter.ToString());
   * }
   * </pre>
   * <p>In case you want to use the enum <b>outside the class</b> such that you give
   * the enum a namespace with its class you have to modify it like follows<p>
   * <pre>
   * class Foo 
   * {
   *   public:
   *     typedef enum { NONE, RADIUS, VOLUME_RADIUS } Filter;
   *     static Enum<Filter> filter;
   *     static void SetEnums();
   * };
   *
   * Enum<Foo::Filter> Foo::filter;
   * 
   * class Boo
   * {
   *    void boo() {
   *       // make somehow sure that Foo::SetEnums() was called exactly one!
   *       Foo::Filter filter = Foo::filter.Parse("radius");
   *    }
   * };
   * </pre> */
  template<typename T>
  class Enum 
  {
     public:
        Enum()
        {
          this->name_ = "";
        }
        
        ~Enum()
        {
        }

        /** optional constructor sets the name of the type */
        Enum(const std::string& name)
        {
          this->name_ = name;
        }

        
        std::map<int, EnumTupel> map;
        
        /** checks if the string has an representation */
        bool IsValid(const std::string& string)
        {
          std::map<int,EnumTupel>::iterator it;
          for(it = map.begin(); it != map.end(); it++)
          {
            if(it->second.string == string) return true;
          } 

          return false;
        }
        
        /** checks if the integer has an representation */
        bool IsValid(T key)
        {
          return map.find((int) key) != map.end();
        }
        
        /** converts the string representation to the type.
         * You can use isValid() first to avoid the exceptzion  
         * @throw an exception if the string is not found. */
        T Parse(const std::string& value)
        {
          std::map<int,EnumTupel>::iterator it;
          for(it = map.begin(); it != map.end(); it++)
          {
            // Check comments on Tupel!!
            if(it->second.string == value) return (T) it->second.key;
          }

          throw Exception("There is no enum key '" + value + "' for '" + name_ + "'"); 
        }

        /** Takes the ->ToString() value */
        T Parse(ParamNode* pn)
        {
          return (T) Parse(pn->AsString());
        }

        /** converts the enumeration to the string.
         * @param type as enums are not typesafe actually an arbitrary int
         * @throw an exception if the type is invalid */
        const std::string& ToString(T type)
        {
          std::map<int,EnumTupel>::iterator it = map.find((int) type);
          if(it != map.end()) return it->second.string;

          EXCEPTION("Invalid key " << type << " for '" + name_ + "'");
        }  
        
        /** Inserts a key/value pair which has to be by default unique. If one
         * allows additional values the map get's a generated/hidden key internaly */
        void Add(T key, const std::string& value, bool force_uniqueness = true)
        {
          // for further values this is a generated key
          int working_key = (int) key;

          if(working_key <= -3092006)
            EXCEPTION("enum values <= -3092006 are reserved");
          
          // check if we already had the key
          if(IsValid(key))
          {
            if(force_uniqueness)
              EXCEPTION("You want to set " << key << ":'" << value << "' in enum "
                        << name_ << " but key is already used");
            working_key = -3092006 - map.size();
          }
          
          // check if we already has the value
          if(IsValid(value))
            EXCEPTION("You want to set " << key << ":'" << value << "' in enum "
                      << name_ << " but value is already used");
          
          map[working_key].key    = (int) key;
          map[working_key].string = value;
        }

        
        void SetName(const std::string& name) 
        {
            this->name_ = name;
        };   
         
     private:
        std::string name_;   
  };

} // end of namespace  


#endif /*ENUM_HH_*/

