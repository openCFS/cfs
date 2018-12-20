/*
 * LibXml2.cc
 *
 *  Created on: 28.05.2016
 *      Author: fwein
 */

#include <sstream>
#include <stdexcept>
#include <string>
#include <iostream>
#include <boost/algorithm/string/trim.hpp>

#include "DataInOut/ParamHandling/LibXml2.hh"
#include "General/Exception.hh"
#include "DataInOut/Logging/LogConfigurator.hh"
#include "DataInOut/Logging/log.hpp"

DECLARE_LOG(libxml2)
DEFINE_LOG(libxml2, "libxml2")

using std::string;

namespace CoupledField
{

PtrParamNode LibXml2::ParseFile(const std::string& file, const std::string& schema)
{
  if(schema != "")
  {
    xmlDoc* doc = xmlParseFile(file.c_str());
    if(doc == NULL)
      throw Exception("cannot open and parse '" + file + "'");
    return ParseAndValidated(doc, schema); // frees doc
  }
  else
  {
    xmlTextReader* reader = xmlNewTextReaderFilename(file.c_str());
    if (reader == NULL)
      throw Exception("cannot open and parse '" + file + "'");
    return ParseReader(reader); // frees the reader
  }
}

PtrParamNode LibXml2::ParseString(const std::string& data, const std::string& schema)
{
  if(schema != "")
  {
    xmlDoc* doc = xmlParseMemory(data.c_str(), data.size());
    if(doc == NULL)
      throw Exception("cannot initialize libxml2 dom and parse memory data");
    return ParseAndValidated(doc, schema);
  }
  else
  {
    xmlTextReader* reader = xmlReaderForMemory(data.c_str(), data.size(), "", NULL, 0);
    if (reader == NULL)
      throw Exception("cannot initialize libxml2 reader and parse memory data");
    return ParseReader(reader); // frees the reader
  }
}


void HandleValidationError(void *ctx, const char *format, ...)
{
    const int buffer_size = 512;
    char *errMsg = new char[512];
    va_list args;
    va_start(args, format);
    vsnprintf(errMsg, buffer_size, format, args);
    va_end(args);
    string msg(errMsg);
    free(errMsg);

    // remove any '{http://www.cfs++.org}' it doesn't improve readability
    // C++ has such a poor lib :(
    for (string::size_type i = msg.find("{http://www.cfs++.org}"); i != string::npos; i = msg.find("{http://www.cfs++.org}"))
       msg.erase(i, strlen("{http://www.cfs++.org}"));
    throw Exception("validation error " + msg);
}

PtrParamNode LibXml2::ParseAndValidated(xmlDoc* doc, const std::string& schemafile)
{
  // see http://stackoverflow.com/questions/6284827/why-does-this-xml-validation-via-xsd-fail-in-libxml2-but-succeed-in-xmllint-an

  xmlSchemaParserCtxt* parserCtxt = NULL;
  xmlSchema* schema = NULL;
  xmlSchemaValidCtxt* validCtxt = NULL;

  xmlNode *root_element = NULL;

  // do we catch an exception? Store it here to rethrow after cleanup, there is no finally in C++
  const Exception* caught = NULL; // Note that we may not create an instance here as this will segfault intenionally with cfs -f

  PtrParamNode pn(new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
  try
  {
    assert(doc != NULL);

    parserCtxt = xmlSchemaNewParserCtxt(schemafile.c_str());
    if(parserCtxt == NULL)
      throw Exception("failed to load schema file " + schemafile);

    schema = xmlSchemaParse(parserCtxt);
    if(schema == NULL)
      throw Exception("failed to parse XSD schema");

    validCtxt = xmlSchemaNewValidCtxt(schema);
    if (!validCtxt)
      throw Exception("Could not create XSD schema validation context");

    xmlSetStructuredErrorFunc(NULL, NULL);
    xmlSetGenericErrorFunc(NULL, HandleValidationError);
    xmlThrDefSetStructuredErrorFunc(NULL, NULL);
    xmlThrDefSetGenericErrorFunc(NULL, HandleValidationError);

    // !!!! set default values !!!!!!
    xmlSchemaSetValidOptions(validCtxt, XML_SCHEMA_VAL_VC_I_CREATE);

    int result = xmlSchemaValidateDoc(validCtxt, doc);

    if(result == 0)
    {
      // Get the root element node
      root_element = xmlDocGetRootElement(doc);
      // parse the stuff to our own structure
      RecursiveFill(root_element, pn, 0);
    }
  }
  catch(const Exception& e)
  {
    // cheat missing finally
    caught = &e;
  }

  // free up the resulting stuff, hopefully in the right order
  if(parserCtxt)
    xmlSchemaFreeParserCtxt(parserCtxt);

  if(schema)
    xmlSchemaFree(schema);

  if(validCtxt != NULL)
    xmlSchemaFreeValidCtxt(validCtxt);

  if(doc != NULL)
    xmlFreeDoc(doc);

  xmlCleanupParser(); // take care for multithreading. Would work without!

  // was there an error on the way, suppressed to clean-up?
  if(caught != NULL)
    throw *caught;

  return pn;
}

void LibXml2::RecursiveFill(xmlNode* node, PtrParamNode parent, int level)
{
  // std::cout << "Fill parent=" << parent->name << " level=" << level << std::endl;
  // in general we are a set of syblings on the same level, a special
  // case is if we are the root elememnt, then we have no sybling
  assert(!(level == 0 && node->next != NULL));

  for(xmlNode* cur = node; cur != NULL; cur = cur->next)
  {
    // std::cout << "ct=" << cur->type << " cn=" << (cur->next == NULL ? "NULL" : "not-null") << std::endl;

    PtrParamNode pn;

    if(cur->type == XML_ELEMENT_NODE)
    {
      if(level > 0)
      {
        PtrParamNode tmp(new ParamNode(ParamNode::EX, ParamNode::ELEMENT));
        parent->AddChildNode(tmp);
        pn = tmp;
      }
      else
      {
        pn = parent;
        assert(cur->next == NULL);
      }
      pn->SetName(string((char*) cur->name));

      //if(level > 0)
      //  std::cout  << "set " << parent->children.back().name << " level=" << parent->children.back().level << std::endl;

      xmlAttr* att = cur->properties;
      while(att != NULL)
      {
        PtrParamNode tmp(new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
        tmp->SetName(string((char*) att->name));
        tmp->SetValue(string((char*) att->children->content));
        pn->AddChildNode(tmp);
        att = att->next;
      }
    }
    if(node->type == XML_TEXT_NODE)
    {
      char* c_ptr = (char*) node->content;
      // one might search the data for numbers, letters and symbols
      // to be faster than convert to string and trim and then throw aways
      std::string cnt(const_cast<const char*>(c_ptr));
      boost::trim(cnt);
      if(cnt.size() > 0)
        parent->SetValue(cnt);
    }
    RecursiveFill(cur->children, pn ? pn : parent, level +1);
  }
}

void LibXml2::Dump(xmlNode* node, int level)
{
    for (xmlNode* cur_node = node; cur_node != NULL; cur_node = cur_node->next)
    {
      if(cur_node->type == XML_ELEMENT_NODE)
      {
        std::cout << "level=" << level << " elem " << cur_node->name << " att " << (cur_node->properties != NULL ? "true" : "false") << std::endl;
        // if(cur_node->content != NULL)
         // std::cout << "  ==> " << cur_node->content << std::endl;
        // if node.lsCountNode() < 2 and len(node.content) > 0:
        //  print "  => " + node.content
      }
      else if(cur_node->type == XML_TEXT_NODE)
      {
        char* c_ptr = (char*) cur_node->content;
        if((int) c_ptr[0] != 10)
          std::cout << "level=" << level << " text " << cur_node->name << " first=" <<  (int) c_ptr[0]  << std::endl;
      }
      Dump(cur_node->children, level + 1);
    }
}


PtrParamNode LibXml2::ParseReader(xmlTextReader* reader)
{
  // http://xmlsoft.org/xmlreader.html
  // for the die hards: SAX shall even be faster: http://xmlsoft.org/xmlmem.html
  LOG_DBG(libxml2) << "PR: ParseReader";
  PtrParamNode root(new ParamNode(ParamNode::EX, ParamNode::ELEMENT));

  PtrParamNode current = root;

  int ret = xmlTextReaderRead(reader);
  int last_depth = -1; // depth of current, shall be 0 but fake here for special case first element in while below
  while (ret == 1)
  {
    // process the current element, if this is an xml element create and add it to current
    // if this is text it is the value of current
    int depth = xmlTextReaderDepth(reader); // defines who's pn child this element becomes
    int type  = xmlTextReaderNodeType(reader); //  1 for start element, 15 for end of element, 2 for attributes, 3 for text nodes, ...
    assert(depth <= last_depth+1); // +1 means child of current, <1 means child of current parent (or grand-parent, ...)

    if(type == 1) // # Element
    {
      //std::cout << "PR: create new element d=" << depth << " t=" << type  << " c=" << current->GetName() << " ld=" << last_depth << std::endl;
      //PtrParamNode element = last_depth > 0 ? PtrParamNode(new ParamNode(ParamNode::EX, ParamNode::ELEMENT)) : root; // there is only one root element
      PtrParamNode element;
      if(root->GetChildren().IsEmpty())
        element = root;
      else
        element = PtrParamNode(new ParamNode(ParamNode::EX, ParamNode::ELEMENT));

      xmlChar* name = xmlTextReaderName(reader);
      element->SetName(string((char*) name));
      xmlFree(name);

      // depth==last_depth+1 -> element shall be child of current, nothing to do, we can simply add element to current
      // depth==last_depth   -> element shall be sibling of current, make current it's parent such that we can add element to new current
      // depth < last_depth  -> element shall be child of n-th grand-parent of current
      // std::cout << "PR: d=" << depth << " t=" << type << " n=" << element->GetName() << " c=" << current->GetName() << " ld=" << last_depth << std::endl;
      while(last_depth+1 > depth && current != NULL)
      {
        current = current->GetParent(); // to handle special case of first element with depth=0 we started with last_depth=-1
        last_depth--; // last_depth is the level of current
        // std::cout << " ld=" << last_depth << " c=" << current->GetName() << "#" << std::endl;
      }
      assert(current != NULL);
      assert(last_depth+1 == depth);
      if(last_depth != -1)
        current->AddChildNode(element);

      assert(!(current->GetChildren().IsEmpty() && last_depth != -1));

      // we assume no value but if, then a successive text element

      // add possible attributes
      while(xmlTextReaderMoveToNextAttribute(reader))
      {
        xmlChar* nme = xmlTextReaderName(reader);
        xmlChar* val = xmlTextReaderValue(reader);

        PtrParamNode tmp(new ParamNode(ParamNode::EX, ParamNode::ATTRIBUTE));
        tmp->SetName(string((char*) nme));
        tmp->SetValue(string((char*) val));
        element->AddChildNode(tmp);

        xmlFree(nme);
        xmlFree(val);
      }
      current = element;
      last_depth = depth; // this is the information where current is located,
    }
    if(type == 3) // # text node means value of current element
    {
      // most next nodes are simply new lines or empty lines :(
      xmlChar* value = xmlTextReaderValue(reader);
      // one might search the data for numbers, letters and symbols
      // to be faster than convert to string and trim and then throw aways
      std::string msg((char*) value);
      boost::trim(msg);
      if(msg.size() > 0)
        current->SetValue(msg);
    }

    ret = xmlTextReaderRead(reader);
  }
  xmlFreeTextReader(reader);
  xmlCleanupParser(); // take care for multithreading. Would work without!

  if(ret != 0)
    throw Exception("failed to parse libxml2 reader object");

  return root;
}


} // end of namespace




