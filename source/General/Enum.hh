#ifndef ENUM_HH_
#define ENUM_HH_

#include <map>
#include <string>
#include <iostream>

namespace CoupledField 
{
  class ParamNode;
  
  class Enum 
  {
     public:
        Enum();
        
        ~Enum();

        /** optional constructor sets the name of the type */
        Enum(const std::string& name); 

        /** checks if the string has an representation */
        bool IsValid(const std::string& string);

        /** checks if the integer has an representation */
        bool IsValid(int key);

        /** converts the string representation to the type.
         * You can use isValid() first to avoid the exceptzion  
         * @throw an exception if the string is not found. */
        int Parse(const std::string& value);

        /** Takes the ->ToString() value */
        int Parse(ParamNode* pn);

        /** converts the enumeration to the string.
         * @param type as enums are not typesafe actually an arbitrary int
         * @throw an exception if the type is invalid */
        const std::string& ToString(int type);
        
        /** a simple shortcut to map.insert */
        void Add(int key, const std::string& value); 

        /** This gives the larges (positive!) key value. Useful when the enum is used as
         * an index. Note, that his is NOT the number of keys!
         * @return the largest positive key
         * @exception if empty */
        int LargestKey();  

        std::map<int, std::string> map;
        
        void SetName(const std::string& name) 
        {
            this->name = name;
        };   
         
     private:
        /** an optional name */
        std::string name;   
  };

} // end of namespace  


#endif /*ENUM_HH_*/

