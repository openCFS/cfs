/*
 * XmlReader.cc
 *
 *  Created on: 28.05.2016
 *      Author: fwein
 */
#include <def_use_xerces.hh>
#include <def_use_libxml2.hh>

#include <fstream>
#include <iostream>
#include "DataInOut/ParamHandling/XmlReader.hh"
#include <boost/algorithm/string.hpp>
#include <boost/iostreams/filtering_stream.hpp>

#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

#ifdef USE_XERCES
#include "DataInOut/ParamHandling/Xerces.hh"
#endif

#ifdef USE_LIBXML2
#include "DataInOut/ParamHandling/LibXml2.hh"
#endif

namespace CoupledField
{

PtrParamNode XmlReader::ParseFile(const std::string& filename,
                                  const std::string& schema,
                                  const std::string &schemaUrl)
{
  bool compress = boost::algorithm::ends_with(filename, ".gz"); // see ParamNode::ToFile()

  if(compress)
  {
    // when we are compress we decompress to a stringstream an parse the memeory data
    std::stringstream ss;
    std::ifstream file(filename.c_str(), std::ios_base::in | std::ios_base::binary);
    boost::iostreams::filtering_istream in;
    in.push(boost::iostreams::gzip_decompressor());
    in.push(file);
    boost::iostreams::copy(in, ss);
    return ParseString(ss.str(), schema);
  }


#ifdef USE_XERCES
  Xerces xerces(schema, schemaUrl);
  xerces.SetFile(filename);
  return xerces.CreateParamNodeInstance();
#endif

#ifdef USE_LIBXML2
  return LibXml2::ParseFile(filename, schema);
#endif
  assert(false);
  return PtrParamNode();
}

/** same as parse file but from memory */
PtrParamNode XmlReader::ParseString(const std::string& str,
                                    const std::string& schema,
                                    const std::string &schemaUrl)
{
#ifdef USE_XERCES
  Xerces xerces(schema, schemaUrl);
  xerces.SetString(str);
  return xerces.CreateParamNodeInstance();
#endif

#ifdef USE_LIBXML2
  return LibXml2::ParseString(str, schema);
#endif
  assert(false);
  return PtrParamNode();
}

} // end of namespace

