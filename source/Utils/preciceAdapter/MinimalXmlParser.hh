#ifndef MINIMAL_XML_PARSER_HH
#define MINIMAL_XML_PARSER_HH

#include <string>
#include <vector>
#include <map>
#include <stdexcept>

namespace CoupledField {

  struct XmlNode {
    std::string name;                           // Tag name (e.g. "participant")
    std::map<std::string, std::string> attrs;   // Attributes (e.g. name="mechanicSolver")
    std::string content;                        // Text content (if any)
    std::vector<XmlNode> children;              // Child nodes
  };

  /**
   * @brief Minimal, dependency-free XML parser for the preCICE configuration file.
   *
   * Parses the precice-config.xml into a tree of XmlNode (tags, attributes, text)
   * without pulling a full XML library into the adapter. Only the small XML
   * subset used by preCICE configuration files is supported (no namespaces,
   * no CDATA, no DTD).
   */
  class MinimalXmlParser {
  public:
    MinimalXmlParser(const std::string &filename);
    XmlNode parse();

  private:
    std::string xmlContent;
    size_t pos;

    void skipWhitespace();
    void skipComment();
    std::string parseTagName();
    std::map<std::string, std::string> parseAttributes();
    XmlNode parseElement();
  };
}
#endif // MINIMAL_XML_PARSER_HH
