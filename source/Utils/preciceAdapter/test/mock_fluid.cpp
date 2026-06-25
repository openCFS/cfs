// ---------------------------------------------------------------------------
// mock_fluid - a tiny, generic preCICE "partner" for openCFS coupling tests.
//
// It replaces a heavyweight external solver (e.g. OpenFOAM) in CI: it provides
// the partner meshes and exchanges analytic, fully deterministic data, so the
// openCFS result can be compared against a committed reference (.h5ref).
//
// Everything it does is driven by a small text "spec" file (see make_mesh.py),
// so the same binary can serve many coupled tests:
//
//     participant <name>
//     dimensions  <d>
//     loadVector  <v0> <v1> [<v2>]      # base value written on every "write" mesh
//     rampTime    <T>                   # linear ramp factor min(t/T, 1)
//     mesh <name> write <dataName>
//     <nVerts>
//     <x y [z]> ...
//     mesh <name> read  <dataName>
//     <nVerts>
//     <x y [z]> ...
//     end
//
// Usage:  mock_fluid <precice-config.xml> <spec-file>
// ---------------------------------------------------------------------------
#include <precice/precice.hpp>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct MeshSpec {
  std::string name;
  std::string dataName;
  bool write = false;                       // true: write data, false: read data
  std::vector<double> coords;               // flattened, size = nVerts * dim
  std::vector<precice::VertexID> ids;       // filled by setMeshVertices
};

struct Spec {
  std::string participant;
  int dim = 2;
  std::vector<double> loadVector;
  double rampTime = 1.0;
  std::vector<MeshSpec> meshes;
};

Spec parseSpec(const std::string &path) {
  std::ifstream in(path);
  if (!in) {
    throw std::runtime_error("mock_fluid: cannot open spec file '" + path + "'");
  }
  Spec spec;
  std::string line;
  while (std::getline(in, line)) {
    std::istringstream ls(line);
    std::string key;
    if (!(ls >> key) || key.empty() || key[0] == '#') {
      continue;
    }
    if (key == "participant") {
      ls >> spec.participant;
    } else if (key == "dimensions") {
      ls >> spec.dim;
    } else if (key == "rampTime") {
      ls >> spec.rampTime;
    } else if (key == "loadVector") {
      double v;
      while (ls >> v) spec.loadVector.push_back(v);
    } else if (key == "end") {
      break;
    } else if (key == "mesh") {
      MeshSpec m;
      std::string role;
      ls >> m.name >> role >> m.dataName;
      m.write = (role == "write");
      // next line: vertex count, then that many coordinate lines
      int nv = 0;
      std::string countLine;
      std::getline(in, countLine);
      std::istringstream(countLine) >> nv;
      m.coords.reserve(static_cast<size_t>(nv) * spec.dim);
      for (int i = 0; i < nv; ++i) {
        std::getline(in, line);
        std::istringstream cs(line);
        double c;
        for (int k = 0; k < spec.dim; ++k) {
          cs >> c;
          m.coords.push_back(c);
        }
      }
      m.ids.resize(static_cast<size_t>(nv));
      spec.meshes.push_back(std::move(m));
    }
  }
  if (spec.participant.empty()) {
    throw std::runtime_error("mock_fluid: spec is missing a 'participant' entry");
  }
  if (static_cast<int>(spec.loadVector.size()) < spec.dim) {
    spec.loadVector.resize(spec.dim, 0.0);
  }
  return spec;
}

} // namespace

int main(int argc, char **argv) {
  if (argc != 3) {
    std::cerr << "usage: " << argv[0] << " <precice-config.xml> <spec-file>\n";
    return 1;
  }
  try {
    const std::string configFile = argv[1];
    const Spec spec = parseSpec(argv[2]);

    precice::Participant participant(spec.participant, configFile, 0, 1);

    for (auto &m : const_cast<Spec &>(spec).meshes) {
      participant.setMeshVertices(m.name, m.coords, m.ids);
    }

    const int dim = spec.dim;
    auto rampFactor = [&](double t) {
      return spec.rampTime > 0.0 ? std::min(t / spec.rampTime, 1.0) : 1.0;
    };

    // helper: fill a "write" mesh's buffer with loadVector * factor
    auto fillWrite = [&](const MeshSpec &m, double factor, std::vector<double> &buf) {
      const size_t nv = m.ids.size();
      buf.assign(nv * dim, 0.0);
      for (size_t i = 0; i < nv; ++i)
        for (int k = 0; k < dim; ++k)
          buf[i * dim + k] = spec.loadVector[k] * factor;
    };

    std::vector<double> buf;

    if (participant.requiresInitialData()) {
      for (const auto &m : spec.meshes) {
        if (m.write) {
          fillWrite(m, rampFactor(0.0), buf);
          participant.writeData(m.name, m.dataName, m.ids, buf);
        }
      }
    }

    participant.initialize();

    double t = 0.0;
    while (participant.isCouplingOngoing()) {
      if (participant.requiresWritingCheckpoint()) {
        // serial-explicit: nothing to store; data is purely analytic
      }

      const double dt = participant.getMaxTimeStepSize();
      const double factor = rampFactor(t + dt);

      for (const auto &m : spec.meshes) {
        if (m.write) {
          fillWrite(m, factor, buf);
          participant.writeData(m.name, m.dataName, m.ids, buf);
        } else {
          buf.assign(m.ids.size() * dim, 0.0);
          participant.readData(m.name, m.dataName, m.ids, dt, buf);
        }
      }

      participant.advance(dt);
      t += dt;

      if (participant.requiresReadingCheckpoint()) {
        // serial-explicit: never requested
      }
    }

    participant.finalize();
    return 0;
  } catch (const std::exception &e) {
    std::cerr << "mock_fluid ERROR: " << e.what() << "\n";
    return 1;
  }
}
