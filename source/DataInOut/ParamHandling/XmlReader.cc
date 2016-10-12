/*
 * XmlReader.cc
 *
 *  Created on: 28.05.2016
 *      Author: fwein
 */
#include <def_use_xerces.hh>
#include <def_use_libxml2.hh>

#include "DataInOut/ParamHandling/XmlReader.hh"

#ifdef USE_XERCES
#include "DataInOut/ParamHandling/Xerces.hh"
#endif

#ifdef USE_LIBXML2
#include "DataInOut/ParamHandling/LibXMl2.hh"
#endif

namespace CoupledField
{

PtrParamNode XmlReader::ParseFile(const std::string& file, const std::string& schema)
{
#ifdef USE_XERCES
  Xerces xerces(schema);
  xerces.SetFile(file);
  return xerces.CreateParamNodeInstance();
#endif

#ifdef USE_LIBXML2
  return LibXml2::ParseFile(file, schema);
#endif

}

/** same as parse file but from memory */
PtrParamNode XmlReader::ParseString(const std::string& str, const std::string& schema)
{
#ifdef USE_XERCES
  Xerces xerces(schema);
  xerces.SetString(str);
  return xerces.CreateParamNodeInstance();
#endif

#ifdef USE_LIBXML2
  return LibXml2::ParseString(str, schema);
#endif

}

} // end of namespace

