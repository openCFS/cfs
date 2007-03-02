#include "utils/Enum.hh"
#include "utils/Exception.hh"

#include <iostream>

namespace OLAS 
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


int Enum::Parse(const std::string& string)
{
    std::map<int,std::string>::iterator it;
    for(it = map.begin(); it != map.end(); it++)
    {
        if(it->second == string) return it->first;
    }
    
    throw Exception("There is no enum key for value", string, name, __LINE__, __PRETTY_FUNCTION__); 
}

const std::string& Enum::ToString(int type)
{
    std::map<int,std::string>::iterator it = map.find(type);
    if(it != map.end()) return it->second;
    
    throw Exception("There is no enum for the type", name, type, __LINE__, __PRETTY_FUNCTION__);
}

void Enum::Add(int key, const std::string& value)
{
    // check if we already had the value
    if(map.find(key) != map.end()) throw Exception("enum key already set", value, key, __LINE__, __PRETTY_FUNCTION__, Exception::AUTO); 
    map[key] = value;
}

    
} // end of namespace
