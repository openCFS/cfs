/*
 * XmlReader.hh
 *
 *  Created on: 28.05.2016
 *      Author: fwein
 */

#ifndef DATAINOUT_PARAMHANDLING_XMLREADER_HH_
#define DATAINOUT_PARAMHANDLING_XMLREADER_HH_

#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {

  /** This class holds service methods to read xml data to a ParamNode tree.
   * Depending on compile time configutation xerces or libxml2 are used. */
  class XmlReader
  {
  public:

    /** open, parse and close the file to a ParamNode tree. Without schema a fast parser mode
     * is used.
     * @schema xsd schema if present */
    static PtrParamNode ParseFile(const std::string& file, const std::string& schema = "");

    /** same as parse file but from memory */
    static PtrParamNode ParseString(const std::string& str, const std::string& schema = "");

  private:

  }; // end of class


} // end of namespace





#endif /* DATAINOUT_PARAMHANDLING_XMLREADER_HH_ */
