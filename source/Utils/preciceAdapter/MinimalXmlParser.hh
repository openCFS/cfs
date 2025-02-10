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
