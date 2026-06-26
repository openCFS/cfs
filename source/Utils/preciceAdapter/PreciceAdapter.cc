
#include "IPreciceAdapter.hh"
#include "PreciceAdapter.hh"
#include "MinimalXmlParser.hh"
#include "NodeResult.hh"
#include "ElementResult.hh"

#include "DataInOut/Logging/LogConfigurator.hh"

#include "Domain/Domain.hh"
#include "Domain/Results/ResultInfo.hh"
#include "Domain/Mesh/GridCFS/GridCFS.hh"
#include "TransientDriverPrecice.hh"
#include "PDE/SinglePDE.hh"
#include "Utils/StdVector.hh"
#include "MatVec/Vector.hh"
#include <algorithm>
#include <numeric>
#include <set>
#include <unordered_map>

namespace CoupledField
{

    namespace {

    static ResultType mapDefinedOnToResultType(ResultInfo::EntityUnknownType definedOn)
    {
        switch (definedOn) {
        case ResultInfo::NODE:
            return ResultType::NODE;
        case ResultInfo::EDGE:
        case ResultInfo::FACE:
        case ResultInfo::ELEMENT:
        case ResultInfo::SURF_ELEM:
            // The adapter stores all non-nodal fields in element-style containers.
            return ResultType::ELEMENT;
        default:
            {
                std::string defOn;
                ResultInfo::Enum2String(definedOn, defOn);
                EXCEPTION("PreciceAdapter: Unsupported definedOn type '" << defOn
                          << "' for data exchange. Only nodal and element/surface-element fields are supported.");
            }
        }
    }

    static bool tryGetDefinedOnFromRegisteredResult(
        Domain* domain,
        const std::string& cfsResultName,
        ResultInfo::EntityUnknownType& definedOn)
    {
        if (!domain || !domain->GetResultHandler()) {
            return false;
        }

        auto* resultContextsPtr = domain->GetResultHandler()->GetResultContexts();
        if (!resultContextsPtr) {
            return false;
        }

        for (const auto& entry : *resultContextsPtr) {
            const shared_ptr<BaseResult>& baseResult = entry.first;
            if (!baseResult || !baseResult->GetResultInfo()) {
                continue;
            }
            if (baseResult->GetResultInfo()->resultName == cfsResultName) {
                definedOn = baseResult->GetResultInfo()->definedOn;
                return true;
            }
        }
        return false;
    }

    } // namespace

    // Define / declare logging stream
    DEFINE_LOG(preciceAdapter, "preciceAdapter")

    PreciceAdapter::PreciceAdapter(boost::shared_ptr<ParamNode> paramNode)
   
         : isInit_(false),
        paramNode_(paramNode),
#ifdef USE_PRECICE
      participant_(nullptr),
#endif
      configFileName_(""),
      participantName_(""),
      participantMeshName_(""),
      participantElemMeshName_(""),
      sequenceStep_(),
      nodeNumCoordMap_(),
      cfsNodeNumsVec_(),
      flatCoords_(),
      preciceNodeNumsVec_(),
      elemNumCoordMap_(),
      cfsElemNumsVec_(),
      flatElemCoords_(),
      preciceElemNumsVec_(),
      needElemMesh_(false),
      activeParticipantConfig_(),
      runtimeReadResults_(),
      runtimeWriteResults_(),
      rank_(0),
      size_(1),
      domain_(nullptr),
      pdes_(),
      externalResultsRegistered_(false)
    {
    }

    PreciceAdapter::~PreciceAdapter()
    {
        finalize();
    }

    const PreciceAdapter::CouplingMeshData& PreciceAdapter::getMeshData(const std::string& meshName) const
    {
        auto it = meshDataByName_.find(meshName);
        if (it == meshDataByName_.end()) {
            EXCEPTION("PreciceAdapter: No coupling mesh data available for mesh '" << meshName
                      << "'. Please define entries in fileFormats/preciceCoupling/participantMeshList.");
        }
        return it->second;
    }

    PreciceAdapter::CouplingMeshData& PreciceAdapter::getMeshData(const std::string& meshName)
    {
        auto it = meshDataByName_.find(meshName);
        if (it == meshDataByName_.end()) {
            EXCEPTION("PreciceAdapter: No coupling mesh data available for mesh '" << meshName
                      << "'. Please define entries in fileFormats/preciceCoupling/participantMeshList.");
        }
        return it->second;
    }


    // Helper function: add a unique node to the node container and map
    static void addUniqueNode(GridCFS* gridCFS, int nodenum,
                            std::vector<int>& nodeVec,
                            std::map<unsigned int, Vector<double>>& nodeMap) {
        if (nodeMap.find(nodenum) == nodeMap.end()) {
            nodeVec.push_back(nodenum);
            Vector<double> coord;
            gridCFS->GetNodeCoordinate(coord, nodenum, true);
            nodeMap[nodenum] = coord;
        }
    }

    // Helper function: add a unique element to the element container and map
    static void addUniqueElem(GridCFS* gridCFS, int elemNum,
                            std::vector<int>& elemVec,
                            std::map<unsigned int, Vector<double>>& elemMap) {
        if (elemMap.find(elemNum) == elemMap.end()) {
            elemVec.push_back(elemNum);
            Vector<double> center;
            gridCFS->GetElemCentroid(center, elemNum, true);
            elemMap[elemNum] = center;
        }
    }


    void PreciceAdapter::readConfigurationParameters() {
        // Retrieve parameters from paramNode_ (participant name, mesh list, config file, sequence step).
        paramNode_->Get("fileFormats/preciceCoupling/participantName")->GetValue("name", participantName_);
        paramNode_->Get("fileFormats/preciceCoupling/precice_configFile")->GetValue("file", configFileName_);
        sequenceStep_ = paramNode_->Get("fileFormats/preciceCoupling")->Get("sequenceStep")->As<int>();

        meshDataByName_.clear();

        if (!paramNode_->Get("fileFormats/preciceCoupling")->Has("participantMeshList")) {
            EXCEPTION("PreciceAdapter: Missing required fileFormats/preciceCoupling/participantMeshList.");
        }

        // Expected structure:
        // fileFormats/preciceCoupling/participantMeshList/participantMeshName name="<preCICE mesh>" region="<region>"
        ParamNodeList meshNodes = paramNode_->Get("fileFormats/preciceCoupling")
                                    ->Get("participantMeshList")->GetChildren();
        if (meshNodes.GetSize() == 0) {
            EXCEPTION("PreciceAdapter: participantMeshList must contain at least one participantMeshName entry.");
        }

        participantMeshName_.clear();
        for (int iMesh = 0; iMesh < static_cast<int>(meshNodes.GetSize()); ++iMesh) {
            std::string meshName = meshNodes[iMesh]->Get("name")->As<std::string>();
            std::string regionName = meshNodes[iMesh]->Get("region")->As<std::string>();

            if (participantMeshName_.empty()) {
                // Legacy anchor mesh used by backward-compatible code paths.
                participantMeshName_ = meshName;
            }

            CouplingMeshData& meshData = meshDataByName_[meshName];
            if (std::find(meshData.regionNames.begin(), meshData.regionNames.end(), regionName) == meshData.regionNames.end()) {
                meshData.regionNames.push_back(regionName);
            }
        }

        std::cout << "PreciceAdapter: Configuration parameters:" << "\n";
        std::cout << "  Participant Name: " << participantName_ << "\n";
        std::cout << "  Participant Mesh Entries: " << meshNodes.GetSize() << "\n";
        std::cout << "  Precice Config File: " << configFileName_ << "\n";
        std::cout << "  Sequence Step: " << sequenceStep_ << "\n";
    }

    void PreciceAdapter::readPreciceConfiguration() {
        std::cout << "PreciceAdapter: reading precice configuration from file: " << configFileName_ << "\n";
        // Use PreciceConfigReader to read additional configuration.
        try {
            PreciceConfigReader configReader(configFileName_);
            const auto &participants = configReader.getParticipants();
            bool found = false;
            for (const auto &pc : participants) {
                if (pc.name == participantName_) {
                    activeParticipantConfig_ = pc; // Store the relevant configuration
                    // Mesh usage is resolved later from ResultInfo::MapSolTypeToDefinedOn.
                    needElemMesh_ = false;
                    participantElemMeshName_.clear();
                    LOG_DBG(preciceAdapter) << "Loaded configuration for participant: " << pc.name;
                    found = true;
                    break;
                }
            }
            if (!found) {
                EXCEPTION("Participant " << participantName_ << " not found in precice-config.xml");
            }
        } catch (const std::exception &e) {
            EXCEPTION("Error reading precice-config.xml: " << e.what());
        }
    }


    std::tuple<std::string, SolutionType> PreciceAdapter::convertResultNamesToCFS(const std::string& precicename, bool isWrite) {
        // Characteristic (impedance-matched) coupling: the openCFS quantity is determined by
        // ROLE, not by the (owner-tagged) data name. Every openCFS participant writes its own
        // outgoing invariant w_out = p' + rho0 c0 u_n' (acouCharacteristic, a post-processing
        // SURF_ELEM result) and reads the partner's outgoing invariant as its incoming w_in
        // (acouCharacteristicCoupling, consumed by the characteristicCouplingBC). Handling it by
        // role lets two openCFS acoustic participants couple symmetrically, while the single
        // openCFS<->OpenFOAM case is unchanged (Acoustic writes ...FromAcou -> acouCharacteristic,
        // reads ...FromFluid -> acouCharacteristicCoupling). The data names stay owner-tagged in
        // the preCICE config so the physical meaning never flips with wave direction.
        if (precicename.rfind("AcousticCharacteristic", 0) == 0) {
            if (isWrite)
                return {"acouCharacteristic", SolutionType::ACOU_CHARACTERISTIC};
            else
                return {"acouCharacteristicCoupling", SolutionType::ACOU_CHARACTERISTIC_COUPLING};
        }

        static const std::unordered_map<std::string, std::tuple<std::string, SolutionType>> conversionMap = {
            {"Temperature", {"heatTemperature", SolutionType::HEAT_TEMPERATURE}},
            {"Displacement", {"mechDisplacement", SolutionType::MECH_DISPLACEMENT}},
            // Whole-domain mesh deformation read into the SmoothPDE (prescribed-displacement
            // mode): an externally computed displacement field (e.g. OpenFOAM's volume
            // pointDisplacement, written with 'locations volumeNodes') routed into
            // SMOOTH_DISPLACEMENT. The "Displacement" prefix makes the OpenFOAM FSI adapter
            // reuse its existing Displacement writer (matchingStrings is a prefix match), while
            // the distinct full name avoids clashing with the FSI "Displacement" exchange.
            {"DisplacementVolume", {"smoothDisplacement", SolutionType::SMOOTH_DISPLACEMENT}},
            {"Force", {"mechForce", SolutionType::MECH_FORCE}},
            {"Stress", {"mechNormalStress", SolutionType::MECH_NORMAL_STRESS}},
            {"MagneticFluxDensity", {"magFluxDensity", SolutionType::MAG_FLUX_DENSITY}},
            {"MagneticFieldIntensity", {"magFieldIntensity", SolutionType::MAG_FIELD_INTENSITY}},
            {"JouleLossDensity", {"magJouleLossPowerDensity", SolutionType::MAG_JOULE_LOSS_POWER_DENSITY}},
            // "Pressure" = full thermodynamic pressure (p' + p0), as produced by
            // fluidMechAbsolutePressure. This is what external solvers such as
            // OpenFOAM rhoPimpleFoam expect as a fixedValue BC on p.
            {"Pressure", {"fluidMechAbsolutePressure", SolutionType::FLUIDMECH_ABSOLUTE_PRESSURE}},
            // "PerturbationPressure" = acoustic perturbation p' only (fluidMechPressure).
            // Use this for openCFS-to-openCFS LinFlow coupling where both sides
            // work with perturbation quantities.
            {"PerturbationPressure", {"fluidMechPressure", SolutionType::FLUIDMECH_PRESSURE}},
            {"Velocity", {"fluidMechVelocity", SolutionType::FLUIDMECH_VELOCITY}},
            {"PressureTemporalDerivative", {"acouRhsLoad", SolutionType::ACOU_RHS_LOAD}}
            // NOTE: the characteristic invariants (data names starting with
            // "AcousticCharacteristic") are handled role-based above, not via this map.
        };

        auto it = conversionMap.find(precicename);
        if (it != conversionMap.end())
            return it->second;
        else
            EXCEPTION("Invalid quantity: " << precicename
                    << ". Currently the adapter supports: \"Temperature\", \"Displacement\", "
                    << "\"Force\", \"Stress\", \"MagneticFluxDensity\", \"MagneticFieldIntensity\", "
                    << "\"JouleLossDensity\", \"Pressure\", \"PerturbationPressure\", \"Velocity\", "
                    << "\"PressureTemporalDerivative\", \"AcousticCharacteristicFromAcou\", "
                    << "\"AcousticCharacteristicFromFluid\"");
    }


    void PreciceAdapter::createPreciceParticipant() {
        std::cout << "PreciceAdapter: creating PreCICE participant..." << "\n";
        if (domain_->GetGridMap().size() > 1)
        {
            EXCEPTION("PreciceAdapter: works only with one grid per simulation - until now");
        }

        // add a check that we dont have a multidriver and that it's transient
        if( !dynamic_cast<TransientDriverPrecice*>(domain_->GetSingleDriver()) )
        {
            EXCEPTION("PreciceAdapter: this only works when using a transient simulation")
        }

        // Create the PreCICE participant instance.
        participant_ = std::make_unique<precice::Participant>(participantName_, configFileName_, rank_, size_);

        std::cout << "PreciceAdapter: created PreCICE participant '" << participantName_
                  << "' with config file '" << configFileName_ << "'" << "\n";

        // Create runtime data structures for all configured read/write entries.
        runtimeWriteResults_.clear();
        runtimeReadResults_.clear();

        for (const auto &entry : activeParticipantConfig_.writeData) {
            ResultConfig config;
            config.precicename = entry.name;
            auto [tempName, tempSolutionType] = convertResultNamesToCFS(config.precicename, /*isWrite=*/true);
            config.cfsname = tempName;
            config.solutiontype = tempSolutionType;
            config.meshName = entry.mesh;
            config.quantitydim = participant_->getDataDimensions(config.meshName, config.precicename);

            ResultInfo::EntityUnknownType mappedDefinedOn = ResultInfo::MapSolTypeToDefinedOn(config.solutiontype);
            ResultType mappedType = mapDefinedOnToResultType(mappedDefinedOn);

            // Sanity check: if a result with this name is already registered in openCFS,
            // compare its definedOn against the canonical SolutionType mapping.
            ResultInfo::EntityUnknownType declaredDefinedOn;
            if (tryGetDefinedOnFromRegisteredResult(domain_, config.cfsname, declaredDefinedOn)) {
                ResultType declaredType = mapDefinedOnToResultType(declaredDefinedOn);
                if (declaredType != mappedType) {
                    std::string mappedDefOn;
                    std::string declaredDefOn;
                    ResultInfo::Enum2String(mappedDefinedOn, mappedDefOn);
                    ResultInfo::Enum2String(declaredDefinedOn, declaredDefOn);
                    EXCEPTION("PreciceAdapter: Inconsistent definition for result '" << config.cfsname
                              << "'. SolutionType maps to definedOn='" << mappedDefOn
                              << "' but registered result uses definedOn='" << declaredDefOn << "'.");
                }
            }

            std::cout << "PreciceAdapter: Preparing to register write-data entry: " << config.precicename
                      << " on mesh: " << config.meshName
                      << " with dimension: " << config.quantitydim << "\n";

            std::cout << "CFSName is: " << config.cfsname << "\n";

            if (mappedType == ResultType::ELEMENT) {
                runtimeWriteResults_.push_back(std::make_unique<ElementResult>(config));
                runtimeWriteResults_.back()->allocateData(getMeshData(config.meshName).cfsElemNums.size());
            } else {
                runtimeWriteResults_.push_back(std::make_unique<NodeResult>(config));
                runtimeWriteResults_.back()->allocateData(getMeshData(config.meshName).cfsNodeNums.size());
            }
        }

        std::cout << "PreciceAdapter: Registered " << runtimeWriteResults_.size() << " write-data entries." << "\n";

        for (const auto &entry : activeParticipantConfig_.readData) {
            ResultConfig config;
            config.precicename = entry.name;
            auto [tempName, tempSolutionType] = convertResultNamesToCFS(config.precicename, /*isWrite=*/false);
            config.cfsname = tempName;
            config.solutiontype = tempSolutionType;
            config.meshName = entry.mesh;
            config.quantitydim = participant_->getDataDimensions(config.meshName, config.precicename);

            ResultInfo::EntityUnknownType mappedDefinedOn = ResultInfo::MapSolTypeToDefinedOn(config.solutiontype);
            ResultType mappedType = mapDefinedOnToResultType(mappedDefinedOn);

            // Sanity check against existing openCFS result definitions, if available.
            ResultInfo::EntityUnknownType declaredDefinedOn;
            if (tryGetDefinedOnFromRegisteredResult(domain_, config.cfsname, declaredDefinedOn)) {
                ResultType declaredType = mapDefinedOnToResultType(declaredDefinedOn);
                if (declaredType != mappedType) {
                    std::string mappedDefOn;
                    std::string declaredDefOn;
                    ResultInfo::Enum2String(mappedDefinedOn, mappedDefOn);
                    ResultInfo::Enum2String(declaredDefinedOn, declaredDefOn);
                    EXCEPTION("PreciceAdapter: Inconsistent definition for result '" << config.cfsname
                              << "'. SolutionType maps to definedOn='" << mappedDefOn
                              << "' but registered result uses definedOn='" << declaredDefOn << "'.");
                }
            }

            std::cout << "PreciceAdapter: Preparing to register read-data entry: " << config.precicename
                      << " on mesh: " << config.meshName
                      << " with dimension: " << config.quantitydim << "\n";

            std::cout << "CFSName is: " << config.cfsname << "\n";

            if (mappedType == ResultType::ELEMENT) {
                runtimeReadResults_.push_back(std::make_unique<ElementResult>(config));
                runtimeReadResults_.back()->allocateData(getMeshData(config.meshName).cfsElemNums.size());
            } else {
                runtimeReadResults_.push_back(std::make_unique<NodeResult>(config));
                runtimeReadResults_.back()->allocateData(getMeshData(config.meshName).cfsNodeNums.size());
                std::cout << "PreciceAdapter: Registered " << runtimeReadResults_.size() << " read-data entries." << "\n";
            }
        }

        for (auto& meshEntry : meshDataByName_) {
            const std::string& meshName = meshEntry.first;
            CouplingMeshData& md = meshEntry.second;

            if (md.needsNodeData && md.needsElemData) {
                EXCEPTION("PreciceAdapter: mesh '" << meshName
                          << "' is used for both nodal and element-based exchange. "
                          << "Please split into two meshes (e.g. '<name>-Nodes' and '<name>-Faces').");
            }

            if (md.needsNodeData) {
                md.preciceNodeNums = md.cfsNodeNums;
                participant_->setMeshVertices(meshName, md.flatNodeCoords, md.preciceNodeNums);
                std::cout << "PreciceAdapter: Registered nodal mesh '" << meshName
                          << "' with " << md.cfsNodeNums.size() << " nodes.\n";
            } else if (md.needsElemData) {
                md.preciceElemNums = md.cfsElemNums;
                participant_->setMeshVertices(meshName, md.flatElemCoords, md.preciceElemNums);
                std::cout << "PreciceAdapter: Registered element mesh '" << meshName
                          << "' with " << md.cfsElemNums.size() << " elements.\n";
            }
        }

        // Backward compatibility with legacy members.
        if (meshDataByName_.find(participantMeshName_) != meshDataByName_.end()) {
            CouplingMeshData& md = meshDataByName_[participantMeshName_];
            preciceNodeNumsVec_ = md.preciceNodeNums;
            cfsNodeNumsVec_ = md.cfsNodeNums;
            flatCoords_ = md.flatNodeCoords;
            cfsElemNumsVec_ = md.cfsElemNums;
            flatElemCoords_ = md.flatElemCoords;
        }

    }

    void PreciceAdapter::setupMeshAndCoordinates() {
        std::cout << "PreciceAdapter: setting up mesh and coordinates..." << "\n";
        GridCFS* gridCFS = dynamic_cast<GridCFS*>(domain_->GetGrid());
        if (!gridCFS) {
            EXCEPTION("PreciceAdapter: Domain's grid is not of type GridCFS as expected.");
        }

        // Determine mesh usage (node-based vs element-based) directly from active preCICE data definitions.
        std::set<std::string> usedMeshNames;
        auto collectUsage = [&](const std::vector<DataEntry>& entries, bool isWrite) {
            for (const auto& entry : entries) {
                usedMeshNames.insert(entry.mesh);
                auto mapped = convertResultNamesToCFS(entry.name, isWrite);
                ResultType rt = mapDefinedOnToResultType(ResultInfo::MapSolTypeToDefinedOn(std::get<1>(mapped)));
                CouplingMeshData& md = meshDataByName_[entry.mesh];
                if (rt == ResultType::NODE) {
                    md.needsNodeData = true;
                } else {
                    md.needsElemData = true;
                }
            }
        };
        collectUsage(activeParticipantConfig_.readData, /*isWrite=*/false);
        collectUsage(activeParticipantConfig_.writeData, /*isWrite=*/true);

        if (usedMeshNames.empty()) {
            EXCEPTION("PreciceAdapter: no read-data/write-data mesh usage found for participant '" << participantName_ << "'.");
        }

        int cfsGridDim = gridCFS->GetDim();

        for (const auto& meshName : usedMeshNames) {
            CouplingMeshData& md = meshDataByName_[meshName];

            md.cfsNodeNums.clear();
            md.flatNodeCoords.clear();
            md.preciceNodeNums.clear();
            md.cfsElemNums.clear();
            md.flatElemCoords.clear();
            md.preciceElemNums.clear();

            std::map<unsigned int, Vector<double>> nodeCoordMap;
            std::map<unsigned int, Vector<double>> elemCoordMap;

            const std::vector<std::string>& regions = md.regionNames;
            if (regions.empty()) {
                EXCEPTION("PreciceAdapter: mesh '" << meshName << "' has no assigned regions. "
                          << "Define fileFormats/preciceCoupling/participantMeshList entries with region attribute.");
            }

            for (const auto& regionName : regions) {
                shared_ptr<EntityList> elemlist = gridCFS->GetEntityList(EntityList::ListType::ELEM_LIST, regionName);
                shared_ptr<EntityList> nodelist = gridCFS->GetEntityList(EntityList::ListType::NODE_LIST, regionName);

                if (md.needsNodeData && nodelist) {
                    EntityIterator it = nodelist->GetIterator();
                    for (it.Begin(); !it.IsEnd(); it++) {
                        addUniqueNode(gridCFS, it.GetNode(), md.cfsNodeNums, nodeCoordMap);
                    }
                }

                if ((md.needsElemData || md.needsNodeData) && elemlist) {
                    EntityIterator it = elemlist->GetIterator();
                    for (it.Begin(); !it.IsEnd(); it++) {
                        const Elem* el = it.GetElem();
                        if (md.needsElemData) {
                            addUniqueElem(gridCFS, el->GetElemNum(), md.cfsElemNums, elemCoordMap);
                        }
                        if (md.needsNodeData) {
                            for (int iNode = 0; iNode < static_cast<int>(el->connect.GetSize()); ++iNode) {
                                addUniqueNode(gridCFS, static_cast<int>(el->connect[iNode]), md.cfsNodeNums, nodeCoordMap);
                            }
                        }
                    }
                }
            }

            if (md.needsNodeData && md.cfsNodeNums.empty()) {
                EXCEPTION("PreciceAdapter: mesh '" << meshName << "' requires nodal coupling data but has no nodes.");
            }
            if (md.needsElemData && md.cfsElemNums.empty()) {
                EXCEPTION("PreciceAdapter: mesh '" << meshName << "' requires element coupling data but has no elements.");
            }

            md.flatNodeCoords.resize(md.cfsNodeNums.size() * cfsGridDim);
            for (size_t i = 0; i < md.cfsNodeNums.size(); ++i) {
                const Vector<double>& point = nodeCoordMap[md.cfsNodeNums[i]];
                for (int j = 0; j < cfsGridDim; ++j) {
                    md.flatNodeCoords[i * cfsGridDim + j] = point[j];
                }
            }

            md.flatElemCoords.resize(md.cfsElemNums.size() * cfsGridDim);
            for (size_t i = 0; i < md.cfsElemNums.size(); ++i) {
                const Vector<double>& center = elemCoordMap[md.cfsElemNums[i]];
                for (int j = 0; j < cfsGridDim; ++j) {
                    md.flatElemCoords[i * cfsGridDim + j] = center[j];
                }
            }

            std::cout << "PreciceAdapter: prepared mesh '" << meshName
                      << "' nodes=" << md.cfsNodeNums.size()
                      << " elems=" << md.cfsElemNums.size() << "\n";
        }

        // Backward compatibility for code paths still relying on legacy vectors.
        if (meshDataByName_.find(participantMeshName_) != meshDataByName_.end()) {
            const CouplingMeshData& md = meshDataByName_[participantMeshName_];
            cfsNodeNumsVec_ = md.cfsNodeNums;
            flatCoords_ = md.flatNodeCoords;
            cfsElemNumsVec_ = md.cfsElemNums;
            flatElemCoords_ = md.flatElemCoords;
        }

    }


    void PreciceAdapter::initialize(Domain *domain, SinglePDE* pde) {
#ifdef USE_PRECICE
    domain_ = domain;

    // Participant-level setup runs exactly once. initialize() is invoked once per
    // coupled PDE; re-running readConfigurationParameters()/setupMeshAndCoordinates()
    // would clear meshDataByName_ (including the preCICE vertex IDs populated by
    // createPreciceParticipant) and corrupt subsequent data exchange.
    if(!isInit_){
        // Step 1: Retrieve configuration parameters from the openCFS ParamNode.
        readConfigurationParameters();

        // Step 2: Parse the precice-config.xml file
        readPreciceConfiguration();

        // Bail out if this driver is not (yet) in the precice sequence step.
        TransientDriverPrecice* tp = dynamic_cast<TransientDriverPrecice*>(domain_->GetSingleDriver());
        if(!tp || !(tp->GetActSequenceStep() == sequenceStep_)){
            std::cout << "PreciceAdapter: Not the correct sequence step, skipping precice initialization" << "\n";
            return;
        }

        // Step 3: Extract and flatten node/element coordinates.
        setupMeshAndCoordinates();

        // Step 4: Create the PreCICE participant.
        std::cout << "Not yet initialized" << "\n";

        createPreciceParticipant();

        std::cout << "PreciceAdapter: Precice participant created." << "\n";
        
        // This is such a DIRTY HACK!!!! This should actually be set by precice!!!!!!!!!
        // This ONLY WORKS for initializing the "Temperature"
        if(participant_->requiresInitialData()){
            std::cout << "PreciceAdapter: Participant requires initial data - initializing 'Temperature' to 293.15K" << "\n";
            for (auto &result : runtimeWriteResults_) {
                if(result->getConfig().precicename == "Temperature"){
                    if (result->getResultType() == ResultType::NODE) {
                        NodeResult* nr = dynamic_cast<NodeResult*>(result.get());
                        if (nr) {
                            const CouplingMeshData& md = getMeshData(nr->getConfig().meshName);
                            std::vector<double> initVec(nr->getFlatData().size());
                            std::fill(initVec.begin(), initVec.end(), 293.15);
                            participant_->writeData(nr->getConfig().meshName, nr->getConfig().precicename,
                                                    md.preciceNodeNums, initVec);
                        }
                    } else { // ELEMENT-based result.
                        ElementResult* er = dynamic_cast<ElementResult*>(result.get());
                        if (er) {
                            const CouplingMeshData& md = getMeshData(er->getConfig().meshName);
                            std::vector<double> initVec(er->getFlatData().size());
                            std::fill(initVec.begin(), initVec.end(), 293.15);
                            participant_->writeData(er->getConfig().meshName, er->getConfig().precicename,
                                                    md.preciceElemNums, initVec);
                        }
                    }
                }
            }
        }

        participant_->initialize();
        isInit_ = true;
    } else {
        // Already initialized: only register further coupled PDEs that belong to
        // the precice sequence step.
        TransientDriverPrecice* tp = dynamic_cast<TransientDriverPrecice*>(domain_->GetSingleDriver());
        if(!tp || !(tp->GetActSequenceStep() == sequenceStep_)){
            return;
        }
    }

    RegisterSinglePDE(pde);

    RegisterExternalResults();
    
#endif
    }



    void PreciceAdapter::RegisterSinglePDE(SinglePDE* pde){
        // initialize() is called once per coupled PDE; collect them all instead
        // of keeping only the last one, so a participant can drive several
        // (iteratively coupled) PDEs.
        if (std::find(pdes_.begin(), pdes_.end(), pde) == pdes_.end()) {
            pdes_.push_back(pde);
            std::cout << "PreciceAdapter: Registered SinglePDE (" << pdes_.size() << " total)." << "\n";
        }
    }


    void PreciceAdapter::RegisterTimeStepReadData(){
        // initial state case: wait until all coupled PDEs are initialized before
        // pushing read data into their result contexts.
        if(pdes_.empty()){
            return;
        }
        for(auto* pde : pdes_){
            if(!pde->IsInitialized()){
                return;
            }
        }
        // could happen if we are not in the right sequencestep
        TransientDriverPrecice* tp = dynamic_cast<TransientDriverPrecice*>(domain_->GetSingleDriver());
        if(!tp){
            return;
        }else{
            if(!(tp->GetActSequenceStep() == sequenceStep_)){
                return;
            }
        }

       
        //TODO handle nodal and element cases. Maybe handle it based on the particular result
        // like expectResultType nodes or elements
        // loop over the required results that we need - defined in the precice config
        for (auto &result : runtimeReadResults_) {
            LOG_DBG(preciceAdapter) << "Participant " << activeParticipantConfig_.name
                                    << " reads data: " << result->getConfig().cfsname
                                    << " on mesh: " << result->getConfig().meshName;

            const CouplingMeshData& md = getMeshData(result->getConfig().meshName);

            // Depending on whether the result is node- or element-based, call participant_->readData
            if (result->getResultType() == ResultType::ELEMENT ||
                ResultInfo::MapSolTypeToDefinedOn(result->getConfig().solutiontype) == ResultInfo::SURF_ELEM) {
                result->allocateData(md.cfsElemNums.size());
                participant_->readData(result->getConfig().meshName,
                                        result->getConfig().precicename,
                                        md.preciceElemNums,
                                        dynamic_cast<TransientDriverPrecice*>(domain_->GetSingleDriver())->GetDeltaT(),
                                        const_cast<std::vector<double>&>(result->getFlatData()));
                // Let the result object convert the flat data into its internal map.
                result->readData(md.cfsElemNums,
                                dynamic_cast<TransientDriverPrecice*>(domain_->GetSingleDriver())->GetDeltaT());
            }
            else { // Node-based result
                result->allocateData(md.cfsNodeNums.size());
                participant_->readData(result->getConfig().meshName,
                                        result->getConfig().precicename,
                                        md.preciceNodeNums,
                                        dynamic_cast<TransientDriverPrecice*>(domain_->GetSingleDriver())->GetDeltaT(),
                                        const_cast<std::vector<double>&>(result->getFlatData()));
                result->readData(md.cfsNodeNums,
                                dynamic_cast<TransientDriverPrecice*>(domain_->GetSingleDriver())->GetDeltaT());
            }
        }


        // ================================== START: Update OpenCFS result contexts

        ResultHandler* resHandler = domain_->GetResultHandler();
        auto* resultContextsPtr = resHandler->GetResultContexts();

        // Iterate over each result context from OpenCFS
        for (auto &entry : *resultContextsPtr) {
            shared_ptr<BaseResult> baseResult = entry.first;
            // Compare the result name with the configuration of our result objects.
            std::string cfsResultName = baseResult->GetResultInfo()->resultName;
            
            bool initializedVec = false;
            int quantityDim = 0;

            for (const auto &result : runtimeReadResults_) {
                if (result->getConfig().cfsname != cfsResultName) {
                    continue;
                }

                Result<Double>& actSol = static_cast<Result<Double>&>(*baseResult);
                EntityIterator it = actSol.GetEntityList()->GetIterator();
                Vector<Double>& vec = actSol.GetVector();
                if (!initializedVec) {
                    quantityDim = result->getConfig().quantitydim;
                    vec.Resize(it.GetSize() * quantityDim);
                    vec.Init();
                    initializedVec = true;
                }

                if (result->getResultType() == ResultType::ELEMENT ||
                    ResultInfo::MapSolTypeToDefinedOn(result->getConfig().solutiontype) == ResultInfo::SURF_ELEM) {
                    ElementResult* elemResult = dynamic_cast<ElementResult*>(result.get());
                    if (!elemResult) {
                        EXCEPTION("Expected an ElementResult for element-based data.");
                    }
                    const auto& elemMap = elemResult->getElementResultMap();
                    for (it.Begin(); !it.IsEnd(); it++) {
                        const Elem* el = it.GetElem();
                        auto search = elemMap.find(el->GetElemNum());
                        if (search != elemMap.end()) {
                            const Vector<double>& tempField = search->second;
                            for (int iDim = 0; iDim < quantityDim; iDim++) {
                                vec[it.GetPos() * quantityDim + iDim] = tempField[iDim];
                            }
                        }
                    }
                } else if (result->getResultType() == ResultType::NODE) {
                    NodeResult* nodeResult = dynamic_cast<NodeResult*>(result.get());
                    if (!nodeResult) {
                        EXCEPTION("Expected a NodeResult for node-based data.");
                    }
                    const auto& nodeMap = nodeResult->getNodeResultMap();
                    for (it.Begin(); !it.IsEnd(); it++) {
                        int nodeNum = it.GetNode();
                        auto search = nodeMap.find(nodeNum);
                        if (search != nodeMap.end()) {
                            const Vector<double>& tempField = search->second;
                            for (int iDim = 0; iDim < quantityDim; iDim++) {
                                vec[it.GetPos() * quantityDim + iDim] = tempField[iDim];
                            }
                        }
                    }
                }
            }
        }


    }



    void PreciceAdapter::RegisterTimeStepWriteData(){
        TransientDriverPrecice* tp = dynamic_cast<TransientDriverPrecice*>(domain_->GetSingleDriver());
        if(!tp){
            return;
        }else{
            if(!(tp->GetActSequenceStep() == sequenceStep_)){
                return;
            }
        }

        // Get the result contexts from the OpenCFS result handler.
        ResultHandler* resHandler = domain_->GetResultHandler();
        auto* resultContextsPtr = resHandler->GetResultContexts();

        // Build a vector of contexts using boost::shared_ptr.
        std::vector<std::pair<boost::shared_ptr<BaseResult>, boost::shared_ptr<ResultHandler::ResultContext>>> contexts;
        for (auto &entry : *resultContextsPtr) {
            contexts.push_back(std::make_pair(entry.first, entry.second));
        }

        // Process each write result.
        for (auto &result : runtimeWriteResults_) {
            const CouplingMeshData& md = getMeshData(result->getConfig().meshName);
            if (result->getResultType() == ResultType::NODE) {
                NodeResult* nr = dynamic_cast<NodeResult*>(result.get());
                std::cout << "WE ACTUALLY HAVE NODE\n"; 
                if (nr) {
                    nr->updateFromOpenCFS(contexts, md.cfsNodeNums);
                    nr->writeData(md.cfsNodeNums);
                    participant_->writeData(nr->getConfig().meshName, nr->getConfig().precicename,
                                            md.preciceNodeNums, nr->getFlatData());
                }
            } else { // ELEMENT-based result.
                ElementResult* er = dynamic_cast<ElementResult*>(result.get());
                std::cout << "WE ACTUALLY HAVE ELEMENT\n"; 
                if (er) {
                    er->updateFromOpenCFS(contexts, md.cfsElemNums);
                    er->writeData(md.cfsElemNums);
                    participant_->writeData(er->getConfig().meshName, er->getConfig().precicename,
                                            md.preciceElemNums, er->getFlatData());
                }
            }
        }

    }





    void PreciceAdapter::MarkReadResultsUpdated()
    {
#ifdef USE_PRECICE
        if (runtimeReadResults_.empty()) return;

        ResultHandler* resHandler = domain_->GetResultHandler();
        auto* resultContextsPtr = resHandler->GetResultContexts();

        for (auto &entry : *resultContextsPtr) {
            shared_ptr<BaseResult> baseResult = entry.first;
            std::string cfsResultName = baseResult->GetResultInfo()->resultName;

            for (const auto &result : runtimeReadResults_) {
                if (result->getConfig().cfsname == cfsResultName) {
                    std::cout << "PreciceAdapter::MarkReadResultsUpdated: marking '"
                              << cfsResultName << "' as updated\n";
                    resHandler->UpdateResult(baseResult);
                    break;
                }
            }
        }
#endif
    }

    void PreciceAdapter::RegisterExternalResults()
    {
    #ifdef USE_PRECICE
    // The read results are participant-level and fully known once the participant
    // has been created. initialize() runs once per coupled PDE, so register only
    // on the first call and skip on subsequent ones.
    if (externalResultsRegistered_) return;

    std::cout << "PreciceAdapter: Registering external results with OpenCFS result handler..." << "\n";

    ResultHandler* resHandler = domain_->GetResultHandler();
    if(!resHandler)
        EXCEPTION("PreciceAdapter::RegisterExternalResults: No result handler available.");

    std::cout << "runtimeReadResults.size() = " << runtimeReadResults_.size() << "\n";

    if(runtimeReadResults_.size() == 0) { externalResultsRegistered_ = true; return; }

    GridCFS* gridCFS = dynamic_cast<GridCFS*>(domain_->GetGrid());

    // Each read result is bound to the openCFS region(s) of its preCICE mesh, as
    // declared in fileFormats/preciceCoupling/participantMeshList. The coupling
    // region (not any single PDE) defines where the field lives, so registration
    // is independent of which / how many PDEs consume it.
    std::set<std::string> registeredResultKeys;
    for(auto &result : runtimeReadResults_) {
        const CouplingMeshData& md = getMeshData(result->getConfig().meshName);
        if (md.regionNames.empty()) {
            EXCEPTION("PreciceAdapter::RegisterExternalResults: mesh '"
                      << result->getConfig().meshName << "' has no assigned regions.");
        }

        // Select entity list type and definedOn based on node- vs element-based result.
        EntityList::ListType listType;
        ResultInfo::EntityUnknownType definedOnType;
        if (result->getResultType() == ResultType::NODE) {
            listType = EntityList::ListType::NODE_LIST;
            definedOnType = ResultInfo::NODE;
        } else {
            listType = EntityList::ListType::ELEM_LIST;
            definedOnType = ResultInfo::ELEMENT;
        }

        for (const auto& regionName : md.regionNames) {
            std::string regKey = result->getConfig().cfsname + "#"
                + std::to_string(static_cast<int>(result->getResultType())) + "#" + regionName;
            if (registeredResultKeys.find(regKey) != registeredResultKeys.end()) {
                continue;
            }
            registeredResultKeys.insert(regKey);

            shared_ptr<EntityList> entityList = gridCFS->GetEntityList(listType, regionName);
            if (!entityList) {
                EXCEPTION("PreciceAdapter::RegisterExternalResults: Could not get entity list for result "
                          << result->getConfig().cfsname << " on region " << regionName);
            }

            shared_ptr<BaseResult> sol(new Result<Double>());
            sol->SetEntityList(entityList);

            shared_ptr<ResultInfo> ri(new ResultInfo());
            ri->resultName = result->getConfig().cfsname;
            ri->resultType = result->getConfig().solutiontype;
            ri->definedOn  = definedOnType;
            if (result->getConfig().quantitydim <= 1) {
                ri->entryType = ResultInfo::SCALAR;
                ri->dofNames = "";
            } else {
                ri->entryType = ResultInfo::VECTOR;
                if (result->getConfig().quantitydim == 2 || result->getConfig().quantitydim == 3) {
                    // preCICE coupling is Cartesian in current adapter usage
                    ri->SetVectorDOFs(result->getConfig().quantitydim, false);
                } else {
                    // Generic fallback for uncommon vector dimensions.
                    ri->dofNames.Clear();
                    for (int i = 0; i < result->getConfig().quantitydim; ++i) {
                        ri->dofNames.Push_back("comp" + std::to_string(i));
                    }
                }
            }
            sol->SetResultInfo(ri);

            StdVector<std::string> outDest;
            outDest.Push_back("");

            resHandler->RegisterResult(sol, /*functor*/ nullptr, sequenceStep_, 0, 1,
                                       domain_->GetSingleDriver()->GetNumSteps(),
                                       outDest, "", true, false);
        }
    }

    externalResultsRegistered_ = true;
    #endif
    }

    Vector<Double> PreciceAdapter::GetElemResult(SolutionType solType, int elemNum)
    {
        std::string solTypeName = SolutionTypeEnum.ToString(solType);
        bool hasMatchingSolutionType = false;
        bool hasElementTypedMatch = false;
        size_t totalCoupledElements = 0;
        for (const auto &result : runtimeReadResults_) {
            if (result->getConfig().solutiontype == solType) {
                hasMatchingSolutionType = true;
                if (result->getResultType() == ResultType::ELEMENT ||
                    ResultInfo::MapSolTypeToDefinedOn(result->getConfig().solutiontype) == ResultInfo::SURF_ELEM) {
                    hasElementTypedMatch = true;
                    ElementResult* er = dynamic_cast<ElementResult*>(result.get());
                    if (er) {
                        const auto& mapRef = er->getElementResultMap();
                        totalCoupledElements += mapRef.size();
                        auto it = mapRef.find(elemNum);
                        if (it != mapRef.end()) {
                            return it->second;
                        }
                        continue;
                    }
                    EXCEPTION("PreciceAdapter: internal error while accessing element result for solution type '"
                              << solTypeName << "' (" << solType << ")"
                              << " (dynamic cast to ElementResult failed).");
                }
            }
        }
        if (hasElementTypedMatch) {
            EXCEPTION("PreciceAdapter: solution type '" << solTypeName << "' (" << solType << ")"
                      << " is available, but element " << elemNum << " was not found in any coupled surface mesh. "
                      << "Received values for " << totalCoupledElements << " coupled elements across all matching meshes.");
        }
        if (hasMatchingSolutionType) {
            EXCEPTION("PreciceAdapter: solution type '" << solTypeName << "' (" << solType << ")"
                      << " is registered as non-element result. "
                      << "GetElemResult() was called with incompatible result type.");
        }
        EXCEPTION("PreciceAdapter: no result with solution type '" << solTypeName
                  << "' (" << solType << ") found");
    }

    Vector<Double> PreciceAdapter::GetNodeResult(SolutionType solType, int nodeNum)
    {
        std::string solTypeName = SolutionTypeEnum.ToString(solType);
        bool hasMatchingSolutionType = false;
        bool hasNodeTypedMatch = false;
        size_t totalCoupledNodes = 0;
        for (const auto &result : runtimeReadResults_) {
            if (result->getConfig().solutiontype == solType) {
                hasMatchingSolutionType = true;
                if (result->getResultType() == ResultType::NODE) {
                    hasNodeTypedMatch = true;
                    NodeResult* er = dynamic_cast<NodeResult*>(result.get());
                    if (er) {
                        const auto& mapRef = er->getNodeResultMap();
                        totalCoupledNodes += mapRef.size();
                        auto it = mapRef.find(nodeNum);
                        if (it != mapRef.end()) {
                            return it->second;
                        }
                    }
                }
            }
        }
        if (hasNodeTypedMatch) {
            EXCEPTION("PreciceAdapter: solution type '" << solTypeName << "' (" << solType << ")"
                      << " is available, but node " << nodeNum << " was not found in any coupled surface mesh. "
                      << "Received values for " << totalCoupledNodes << " coupled nodes across all matching meshes.");
        }
        if (hasMatchingSolutionType) {
            EXCEPTION("PreciceAdapter: solution type '" << solTypeName << "' (" << solType << ")"
                      << " is registered as non-node result. "
                      << "GetNodeResult() was called with incompatible result type.");
        }
        EXCEPTION("PreciceAdapter: no result with solution type '" << solTypeName
                  << "' (" << solType << ") found");
    }

    void PreciceAdapter::finalize()
    {
#ifdef USE_PRECICE
        if (participant_)
        {
            participant_->finalize();
            participant_.reset();
            LOG_DBG(preciceAdapter) << "PreCICE participant finalized.";
        }
#endif
    }



} // end of namespace
