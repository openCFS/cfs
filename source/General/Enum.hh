#ifndef ENUM_HH_
#define ENUM_HH_

#include <map>
#include <string>
#include <cassert>

#include "General/defs.hh"
#include "General/Exception.hh"
#include "Utils/StdVector.hh"

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
   *    filter_ = filter.Parse("paramNode->Get("radius"));
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
        bool IsValid(const std::string& string) const
        {
          EnumMap::const_iterator it, end;
          it = map.begin();
          end = map.end();
          
          for( ; it != end; it++)
          {
            if(it->second == string) return true;
          } 

          return false;
        }
        
        /** checks if the integer has an representation */
        bool IsValid(T key) const
        {
          return map.find(static_cast<Integer>(key)) != map.end();
        }
        
        /** Unfortunately, the following code gives link errors :(
         *  Uses the string represenation of the ParamNode value and calls
         * @see Parse(const std::string&) */
        // T Parse(const PtrParamNode& pn) const
        // {
        //  return Parse(pn->AsConst<std::string>());
        // }

        /** converts the string representation to the type.
         * You can use isValid() first to avoid the exception
         * @throw an exception if the string is not found. */
        T Parse(const std::string& value) const
        {
          EnumMap::const_iterator it, end;
          it = map.begin();
          end = map.end();
            
          for( ; it != end; it++)
          {
            // Check comments on Tupel!!
            if(it->second == value) return static_cast<T>(it->first);
          }

          throw Exception("There is no enum key '" + value + "' for '" + name_ + "'"); 
        }

        /** A batch processing version of Parse.
         * @param values all names are processed here
         * @param keys_out will be resized to values.GetSize() */
        void Parse(const StdVector<std::string>& values, StdVector<T>& keys_out) const
        {
          keys_out.Resize(values.GetSize());

          for(unsigned int i = 0, vs = values.GetSize(); i < vs; i++)
            keys_out[i] = Parse(values[i]);
        }

        /** converts the string representation to the type if possible
         * if not returns the given default value
         * checking isvalid before is more expensive
         */
        T Parse(const std::string& value, T def) const
        {
          EnumMap::const_iterator it, end;
          it = map.begin();
          end = map.end();
            
          for( ; it != end; it++)
          {
            // Check comments on Tupel!!
            if(it->second == value) return static_cast<T>(it->first);
          }
          
          return(def);
        }

        /** similar to standard parse but does not throw an error */
        T TryParse(const std::string& value, T backup) const
        {
          EnumMap::const_iterator it, end;
          it = map.begin();
          end = map.end();
            
          for( ; it != end; it++)
          {
            // Check comments on Tupel!!
            if(it->second == value) return static_cast<T>(it->first);
          }

          return backup;
        }

        // Commented out, as this introduces another dependency
        // on the ParamNode, which leads to problems in the h5tool.
        // (ahauck, 2010-03-15)
//        /** Takes the ->ToString() value */
//        T Parse(PtrParamNode pn) const
//        {
//          return static_cast<T>(Parse(pn->As<std::string>()));
//        }

        /** converts the enumeration to the string.
         * @param type as enums are not typesafe actually an arbitrary int
         * @throw an exception if the type is invalid */
        const std::string& ToString(T type) const
        {
          EnumMap::const_iterator it = map.find(static_cast<Integer>(type));

          if(it != map.end())
            return it->second;

          EXCEPTION("Invalid key " << type << " for '" + name_ + "'");
        }  
        
        /** Batch processing of ToString().
         * See the Parse() batch processing version */
        void ToString(const StdVector<T>& keys_in, StdVector<std::string>& values_out) const
        {
          values_out.Resize(keys_in.GetSize());

          for(unsigned int i = 0, ks = keys_in.GetSize(); i < ks; i++)
            values_out[i] = ToString(keys_in[i]);
        }

        /** Inserts a key/value pair which has to be by default unique. */
        void Add(T key, const std::string& value, bool force_uniqueness = true)
        {
          // check if we already had the key
          if(IsValid(key))
          {
            if(force_uniqueness)
              EXCEPTION("You want to set " << key << ":'" << value << "' in enum "
                        << name_ << " but key is already used with value " << ToString(key));
          }
          
          // check if we already have the value
          if(force_uniqueness && IsValid(value))
            EXCEPTION("You want to set " << key << ": '" << value << "' in enum "
                      << name_ << " but value is already used");
          
          map.insert( EnumTuple(key,value) );
        }

        
        void SetName(const std::string& name) 
        {
            this->name_ = name;
        };   

        /** for debug purpose */
        void Dump() const
        {
          for(EnumMap::const_iterator it = map.begin() ; it != map.end(); ++it)
            std::cout << "enum " << name_ << " key=" << it->first << " value=" << it->second << std::endl;
        }

     private:
        std::string name_;   
  };

} // end of namespace  


#endif /*ENUM_HH_*/

