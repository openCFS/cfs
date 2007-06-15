#include "General/Enum.hh"
#include "General/exception.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

#include <iostream>

namespace CoupledField 
{

Enum::Enum()
{
    this->name = "";
}

Enum::~Enum()
{
    // nothing dynamic
}

Enum::Enum(const std::string& name)
{
    this->name = name;
}

bool Enum::IsValid(const std::string& string)
{
    std::map<int,std::string>::iterator it;
    for(it = map.begin(); it != map.end(); it++)
    {
        if(it->second == string) return true;
    } 

    return false;
}

bool Enum::IsValid(int key)
{
    return map.find(key) != map.end();
}

int Enum::Parse(ParamNode* pn)
{
  return Parse(pn->AsString());
}


int Enum::Parse(const std::string& string)
{
    std::map<int,std::string>::iterator it;
    for(it = map.begin(); it != map.end(); it++)
    {
        if(it->second == string) return it->first;
    }
    
    EXCEPTION("There is no enum key " << string << " for enum " << name); 
}

const std::string& Enum::ToString(int type)
{
    std::map<int,std::string>::iterator it = map.find(type);
    if(it != map.end()) return it->second;
    
    EXCEPTION("There is no enum value for key " << type << " in enum " << name)
}

void Enum::Add(int key, const std::string& value)
{
    // check if we already had the value
    if(map.find(key) != map.end()) 
       EXCEPTION("enum key " << key << " with value " << value << " in " 
                  << name << " already set"); 
    map[key] = value;
}

int Enum::LargestKey()
{
  if(map.size() == 0)
    throw Exception("no largest key in empty enum " + name);
  
  std::map<int,std::string>::iterator it;
  int largest = 0; // only for positive elements!
  
  for(it = map.begin(); it != map.end(); it++)
    largest = it->first > largest ? it->first : largest;
    
  return largest;  
}

    
} // end of namespace
