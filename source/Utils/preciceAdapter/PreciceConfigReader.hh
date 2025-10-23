#ifndef FILE_CFS_PRECICECONFIGREADER_MINIMAL_HH
#define FILE_CFS_PRECICECONFIGREADER_MINIMAL_HH

#include "MinimalXmlParser.hh"
#include <vector>
#include <string>
#include <stdexcept>
#include <iostream>

namespace CoupledField {


  // A simple structure for mesh-related entries.
  struct MeshEntry {
    std::string name;  // The mesh name (e.g., "mechanicSolver-Mesh")
    std::string from;  // For receive-mesh only, the "from" attribute (empty if not applicable)
  };

  // A structure for data exchange entries.
  struct DataEntry {
    std::string name;  // The data quantity name (e.g., "Temperature")
    std::string mesh;  // The mesh on which this data is defined
  };

  struct ParticipantConfig {
    std::string name;
    std::vector<MeshEntry> provideMeshes;   // All provide-mesh elements
    std::vector<MeshEntry> receiveMeshes;   // All receive-mesh elements
    std::vector<DataEntry> readData;        // All read-data elements
    std::vector<DataEntry> writeData;       // All write-data elements
  };


  class PreciceConfigReader {
  public:
    PreciceConfigReader(const std::string &filename) {
      MinimalXmlParser parser(filename);
      root = parser.parse();
      parseParticipants();
    }

    const std::vector<ParticipantConfig>& getParticipants() const {
      std::cout << "Number of participants: " << participants.size() << std::endl;
      for (const auto &p : participants){
        std::cout << "Participant: " << p.name << std::endl;
      }
      return participants;
    }

  private:
    XmlNode root;
    std::vector<ParticipantConfig> participants;

    // Helper: find all child nodes of a given node that have the specified name.
    std::vector<XmlNode> findChildren(const XmlNode &node, const std::string &childName) {
      std::vector<XmlNode> result;
      for (const auto &child : node.children) {
        // For simplicity, we compare names directly.
        if (child.name == childName) {
          result.push_back(child);
        }
      }
      return result;
    }



    // Parse the <participant> elements under the root.
    void parseParticipants() {
      if (root.name != "precice-configuration") {
        throw std::runtime_error("Expected 'precice-configuration' as root node");
      }

      // Get all participant nodes.
      auto participantNodes = findChildren(root, "participant");
      for (const auto &pNode : participantNodes) {
        ParticipantConfig pc;
        // Participant must have a "name" attribute.
        auto it = pNode.attrs.find("name");
        if (it == pNode.attrs.end()) {
          throw std::runtime_error("Participant node missing 'name' attribute");
        }
        pc.name = it->second;

        // Process provide-mesh elements.
        auto provideNodes = findChildren(pNode, "provide-mesh");
        for (const auto &node : provideNodes) {
          MeshEntry me;
          auto nameIt = node.attrs.find("name");
          if (nameIt != node.attrs.end()) {
            me.name = nameIt->second;
          }
          // For provide-mesh, the "from" field is not used.
          pc.provideMeshes.push_back(me);
        }

        // Process receive-mesh elements.
        auto receiveNodes = findChildren(pNode, "receive-mesh");
        for (const auto &node : receiveNodes) {
          MeshEntry me;
          auto nameIt = node.attrs.find("name");
          if (nameIt != node.attrs.end()) {
            me.name = nameIt->second;
          }
          auto fromIt = node.attrs.find("from");
          if (fromIt != node.attrs.end()) {
            me.from = fromIt->second;
          }
          pc.receiveMeshes.push_back(me);
        }

        // Process read-data elements.
        auto readNodes = findChildren(pNode, "read-data");
        for (const auto &node : readNodes) {
          DataEntry de;
          auto nameIt = node.attrs.find("name");
          if (nameIt != node.attrs.end()) {
            de.name = nameIt->second;
          }
          auto meshIt = node.attrs.find("mesh");
          if (meshIt != node.attrs.end()) {
            de.mesh = meshIt->second;
          }
          pc.readData.push_back(de);
        }

        // Process write-data elements.
        auto writeNodes = findChildren(pNode, "write-data");
        for (const auto &node : writeNodes) {
          DataEntry de;
          auto nameIt = node.attrs.find("name");
          if (nameIt != node.attrs.end()) {
            de.name = nameIt->second;
          }
          auto meshIt = node.attrs.find("mesh");
          if (meshIt != node.attrs.end()) {
            de.mesh = meshIt->second;
          }
          pc.writeData.push_back(de);
        }

        // Add the participant configuration to the list.
        participants.push_back(pc);
      }
    }
  };

} // end namespace CoupledField
#endif // FILE_CFS_PRECICECONFIGREADER_MINIMAL_HH
