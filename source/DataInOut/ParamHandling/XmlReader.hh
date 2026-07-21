/*
 * XmlReader.hh
 *
 *  Created on: 28.05.2016
 *      Author: fwein
 */

#ifndef DATAINOUT_PARAMHANDLING_XMLREADER_HH_
#define DATAINOUT_PARAMHANDLING_XMLREADER_HH_

#include "General/defs.hh"

namespace CoupledField {

  /** This class holds service methods to read xml data to a ParamNode tree.
   * Depending on compile time configutation xerces or libxml2 are used. */
  class XmlReader
  {
  public:

    /** open, parse and close the file to a ParamNode tree. Without schema a fast parser mode
     * is used.
     * @param name if ends with .gz, the file is automatically decompressed
     * @param schema xsd schema if present */
    static PtrParamNode ParseFile(const std::string &file,
                                  const std::string &schema = "",
                                  const std::string &schemaUrl = "");

    /** same as parse file but from memory */
    static PtrParamNode ParseString(const std::string& str,
                                    const std::string& schema = "",
                                    const std::string &schemaUrl = "");

  private:

  }; // end of class


} // end of namespace





#endif /* DATAINOUT_PARAMHANDLING_XMLREADER_HH_ */
