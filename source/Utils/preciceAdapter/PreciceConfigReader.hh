#ifndef FILE_CFS_PRECICECONFIGREADER_MINIMAL_HH
#define FILE_CFS_PRECICECONFIGREADER_MINIMAL_HH

#include "MinimalXmlParser.hh"
#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>

namespace CoupledField {


  struct ParticipantConfig {
    std::string name;
    std::string provideMesh;
    // Later we can add receiveMesh, readData, writeData, etc.
    // but for now, I am happy to get the basics correct
  };

  class MyPreciceConfigReader {
  public:
    MyPreciceConfigReader(const std::string &filename) {
      MinimalXmlParser parser(filename);
      root = parser.parse();
      parseParticipants();
    }

    const std::vector<ParticipantConfig>& getParticipants() const {
      return participants;
    }

  private:
    XmlNode root;
    std::vector<ParticipantConfig> participants;

    // iterates over the root children for nodes with the given name.
    std::vector<XmlNode> findChildren(const XmlNode &node, const std::string &name) {
      std::vector<XmlNode> results;
      for (const auto &child : node.children) {
        if (child.name == name)
          results.push_back(child);
      }
      return results;
    }

    // Parse the <participant> elements under the root.
    void parseParticipants() {
      if (root.name != "precice-configuration") {
        EXCEPTION("Expected precice-configuration as root node");
      }
      auto parts = findChildren(root, "participant");
      for (const auto &p : parts) {
        ParticipantConfig pc;
        // Assume the participant node has an attribute "name"
        if (p.attrs.find("name") != p.attrs.end())
          pc.name = p.attrs.at("name");
        else
          EXCEPTION("Participant node missing 'name' attribute");
        // Look for a child element <provide-mesh>
        auto provideNodes = findChildren(p, "provide-mesh");
        if (!provideNodes.empty()) {
          // Assume the first <provide-mesh> has an attribute "name"
          if (provideNodes[0].attrs.find("name") != provideNodes[0].attrs.end())
            pc.provideMesh = provideNodes[0].attrs.at("name");
        }
        std::cout << "Found participant: " << pc.name 
                  << " provides mesh: " << pc.provideMesh << std::endl;
        participants.push_back(pc);
      }
    }
  };

} // end namespace CoupledField
#endif // FILE_CFS_PRECICECONFIGREADER_MINIMAL_HH
