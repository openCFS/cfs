#include <fstream>
#include <sstream>
#include <iostream>
#include <cstring>
#include <boost/filesystem/operations.hpp>

#include "General/exception.hh"
#include "RapidReader.hh"

using std::cout;
using std::endl;
using std::vector;
using std::string;
using namespace rapidxml;
namespace fs = boost::filesystem;

namespace CoupledField
{

RapidReader::RapidReader(const std::string& filename) :
    filename_(filename)
{
  // do some checks here first if file exists at all
  // create canonical path from native-representation of the
  // file and the schema path
  fs::path filePath = fs::path( filename_, fs::native );
  if(!fs::exists(filePath))
      EXCEPTION("xml file " << filename_ << " doesn't exist");

  // Read file into temporary vector<char>
  std::ifstream densfile(filename_.c_str());
  vector<char> buffer((std::istreambuf_iterator<char>(densfile)), std::istreambuf_iterator<char>( ));
  // essential: do not forget to append terminating character!!
  buffer.push_back('\0');
  // swap the temporary vector with the member variable
  buffer_.swap(buffer);
  
  // parse the density file that is now in buffer_
  doc_.parse<0>(&buffer_[0]);
}


PtrParamNode RapidReader::CreateHeaderParamNode()
{
  // we get the header which is hardcoded here
  xml_node<> *header = doc_.first_node()->first_node();
  if(strcmp("header", header->name()) != 0)
    EXCEPTION("rapidxml: header of density file not found.")
  
  // now we need to traverse the header and construct the ParamNode for all children etc.
  // we do a slightly simpler version for the density file here, since the structure 
  // is very predictable
  PtrParamNode out (new ParamNode(ParamNode::EX, ParamNode::ELEMENT ) );
  out->SetName("header");
  FillChildren(header->first_node(), out);
  // std::cout << std::endl; out->Dump();
  return out;
}

void RapidReader::FillChildren(xml_node<> * node, PtrParamNode parent)
{
  while(node != NULL)
  {
    PtrParamNode pn = PtrParamNode(new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
    pn->SetName(string(node->name()));
    pn->SetValue(node->value());
    parent->AddChildNode(pn);

    // add the attributes of the element
    xml_attribute<> *attr;
    for(attr = node->first_attribute(); attr; attr = attr->next_attribute())
    {
      PtrParamNode newNode = PtrParamNode(new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
      newNode->SetName(string(attr->name()));
      newNode->SetValue(string(attr->value()));
      pn->AddChildNode(newNode);
    }
    
    // check for children
    if(node->first_node() != NULL)
      FillChildren(node->first_node(), pn);
    
    // do the next sibling
    node = node->next_sibling();
  }
}



xml_node<>* RapidReader::FindSet(const std::string& setid)
{
  // we get the first set which is hardcoded here
  xml_node<> *set = doc_.first_node()->first_node()->next_sibling();

  // check if this a set at all 
  // this should be an EXCEPTION!
  if(set == NULL)
    cout << "no sets in density file!" << endl;

  // if the setid is "first", we are done
  bool done("first" == setid);

  // check all other sets for the setid we are searching
  xml_attribute<> *attr;
  while(!done && set->next_sibling() != NULL)
  {
    for(attr = set->first_attribute(); attr; attr = attr->next_attribute())
      if((strcmp("id", attr->name()) == 0) && (setid == string(attr->value())))
        done = true;
    if(!done) // necessary for first set!
      set = set->next_sibling();
  }

  // we either return the set we found or, if we did not find it, the last set
  return set;
}


}
