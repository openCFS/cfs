#ifndef ENUM_HH_
#define ENUM_HH_

#include <map>
#include <string>

#include "General/defs.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField 
{

/** There is a case where for a key we have multiple names. 
 * Something about ParamIdent <-> Harmonic analysis types. Therefore we
 * use the C++ multimap type. 
 * As the Enum is Templated we cannot define this struct inside Enum<T> :( */
typedef std::pair<int, std::string> EnumTuple;
typedef std::multimap<int, std::string> EnumMap;


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
          name_ = "";
        }
        
        ~Enum()
        {
        }

        /** optional constructor sets the name of the type */
        Enum(const std::string& name)
        {
          name_ = name;
        }

        Enum( const std::string& name,
               UInt numTuples,
               EnumTuple * tuples,
               bool forceUniqueness = true)
        {        
          name_ = name;
          for(UInt i=0; i<numTuples; i++) {
        	Add(T(tuples[i].first), tuples[i].second, forceUniqueness);  
          }
        }
        
        EnumMap map;
        
        /** checks if the string has an representation */
        bool IsValid(const std::string& string)
        {
          EnumMap::iterator it, end;
          it = map.begin();
          end = map.end();
          
          for( ; it != end; it++)
          {
            if(it->second == string) return true;
          } 

          return false;
        }
        
        /** checks if the integer has an representation */
        bool IsValid(T key)
        {
          return map.find(static_cast<Integer>(key)) != map.end();
        }
        
        /** converts the string representation to the type.
         * You can use isValid() first to avoid the exceptzion  
         * @throw an exception if the string is not found. */
        T Parse(const std::string& value)
        {
          EnumMap::iterator it, end;
          it = map.begin();
          end = map.end();
            
          for( ; it != end; it++)
          {
            // Check comments on Tupel!!
            if(it->second == value) return static_cast<T>(it->first);
          }

          throw Exception("There is no enum key '" + value + "' for '" + name_ + "'"); 
        }

        /** Takes the ->ToString() value */
        T Parse(ParamNode* pn)
        {
          return static_cast<T>(Parse(pn->AsString()));
        }

        /** converts the enumeration to the string.
         * @param type as enums are not typesafe actually an arbitrary int
         * @throw an exception if the type is invalid */
        const std::string& ToString(T type)
        {
          EnumMap::iterator it = map.find(static_cast<Integer>(type));
          if(it != map.end()) return it->second;

          EXCEPTION("Invalid key " << type << " for '" + name_ + "'");
        }  
        
        /** Inserts a key/value pair which has to be by default unique. */
        void Add(T key, const std::string& value, bool force_uniqueness = true)
        {
          // check if we already had the key
          if(IsValid(key))
          {
            if(force_uniqueness)
              EXCEPTION("You want to set " << key << ":'" << value << "' in enum "
                        << name_ << " but key is already used");
          }
          
          // check if we already have the value
          if(IsValid(value))
            EXCEPTION("You want to set " << key << ":'" << value << "' in enum "
                      << name_ << " but value is already used");
          
          map.insert( EnumTuple(key,value) );
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

