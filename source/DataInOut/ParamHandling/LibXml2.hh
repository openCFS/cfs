/*
 * LibXml2.hh
 *
 *  Created on: 28.05.2016
 *      Author: fwein
 */

#ifndef DATAINOUT_PARAMHANDLING_LIBXML2_HH_
#define DATAINOUT_PARAMHANDLING_LIBXML2_HH_

#include <libxml/parser.h>
#include <libxml/valid.h>
#include <libxml/xmlschemas.h>
#include <libxml/xmlreader.h>

#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {


/** this class is a wrapper for libxml2 for xml parsing, alternatively to Xerces.
 * Xerces fails to work on OS-X with some strange xsd schema issues */
class LibXml2
{
public:
  /** parse a file with DOM when the schema is given or faster when no schema is given */
  static PtrParamNode ParseFile(const std::string& file, const std::string& schema = "");

  /** parse memory data with DOM when the schema is given or faster when no schema is given */
  static PtrParamNode ParseString(const std::string& data, const std::string& schema = "");

private:
  /** use the schema to validate and fill defaults, the DOM based approach */
  static PtrParamNode ParseAndValidated(xmlDoc* doc, const std::string& schema);

  /** recursive helper for ParseAndValidated() */
  static void RecursiveFill(xmlNode* node, PtrParamNode parent, int level);

  /** for debug purpose an recursice dump of an libxml2 dom node */
  static void Dump(xmlNode* node, int level);

  /** this is the libxml2 reader interface, similar to SAX for fast parsing without schema
   * and validation. Meant for density.xml, log.xml and that stuff. Files might be really huge!! */
  static PtrParamNode ParseReader(xmlTextReader* reader);

};

} // end of namespace

#endif /* DATAINOUT_PARAMHANDLING_LIBXML2_HH_ */
