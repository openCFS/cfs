#ifndef TENUM_HH_
#define TENUM_HH_

#include <map>
#include <string>
#include <iostream>

#include "utils/Exception.hh"

namespace OLAS 
{
  template<typename TYPE>
  class TEnum 
  {
     public:
        TEnum()
        {
            this->name = "";
        };


        /** optional constructor sets the name of the type */
        TEnum(const std::string& name)
        {
            this->name = name;
        };

        /** checks if the string has an representation */
        bool IsValid(const std::string& string)
        {
            //std::map<TYPE,std::string>::iterator it;
            //for(it = map.begin(); it != map.end(); it++)
            {
              //  if(it->second == string) return true;
            } 
        
            return false;
        };
        
        
        /** checks if the integer has an representation */
        bool IsValid(TYPE key)
        {
            return map.find(key) != map.end();
        };


        

        /** converts the string representation to the type.
         * You can use isValid() first to avoid the exceptzion  
         * @throw an exception if the string is not found. */
        TYPE Parse(const std::string& string)
        {
            //std::map<TYPE,std::string>::iterator it;
            //for(it = map.begin(); it != map.end(); it++)
            {
//                if(it->second == string) return it->first;
            }
            
            throw Exception("There is no enum key for value", string, name, __LINE__, __PRETTY_FUNCTION__); 
        };

        

        /** converts the enumeration to the string.
         * @param type as enums are not typesafe actually an arbitrary int
         * @throw an exception if the type is invalid */
        const std::string& ToString(TYPE type)
        {
            //std::map<TYPE,std::string>::iterator it = map.find(type);
            //if(it != map.end()) return it->second;
            
            throw Exception("There is no enum for the type", name, type, __LINE__, __PRETTY_FUNCTION__);
        };

        
        
        /** a simple shortcut to map.insert */
        void Add(TYPE key, const std::string& value)
        {
           // check if we already had the value
           if(map.find(key) != map.end()) throw Exception("enum key already set", value, key, __LINE__, __PRETTY_FUNCTION__, Exception::AUTO); 
           map[key] = value;
        };   

        std::map<TYPE, std::string> map;
        
        void SetName(const std::string& name) 
        {
            this->name = name;
        };  
         
     private:
        /** an optional name */
        std::string name;   
  };

} // end of namespace  



#endif /*TENUM_HH_*/
