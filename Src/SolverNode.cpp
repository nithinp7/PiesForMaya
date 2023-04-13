#include "SolverNode.h"

#include <Pies/Solver.h>
#include <maya/MDataHandle.h>
#include <maya/MFnArrayAttrsData.h>
#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnStringData.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MGlobal.h>
#include <maya/MIOStream.h>
#include <maya/MIntArray.h>
#include <maya/MNodeCacheDisablingInfoHelper.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MStatus.h>
#include <maya/MMatrix.h>

#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#define McheckErr(stat, msg)                                                   \
  if (!stat) {                                                                 \
    MGlobal::displayInfo(msg);                                                 \
    stat.perror(msg);                                                          \
    return stat;                                                               \
  }

/*static*/
MString SolverNode::pluginDir;

/*static*/
MTypeId SolverNode::id(0x00000001);

/*static*/
MObject SolverNode::time;

/*static*/
MObject SolverNode::prevTime;

/*static*/
MObject SolverNode::stepSize;

/*static*/
MObject SolverNode::iterations;

/*static*/
MObject SolverNode::simulationEnabled;

/*static*/
MObject SolverNode::simulationStartTime;

/*static*/
MObject SolverNode::meshArray;

/*static*/
MObject SolverNode::ouptutPositions;

static void setupTestScene(Pies::Solver& solver) {
  solver.createTetBox(
      glm::vec3(0.0f, 30.0f, 0.0f),
      1.5f,
      glm::vec3(0.0f),
      1000.0f,
      1.0f,
      true);

  solver.createTetBox(
      glm::vec3(0.0f, 40.0f, 0.0f),
      2.0f,
      glm::vec3(0.0f),
      1000.0f,
      1.0f,
      false);
  solver.createTetBox(
      glm::vec3(10.0f, 40.0f, 0.0f),
      2.0f,
      glm::vec3(0.0f),
      1000.0f,
      1.0f,
      false);
  solver.createTetBox(
      glm::vec3(-10.0f, 40.0f, 0.0f),
      2.0f,
      glm::vec3(0.0f),
      1000.0f,
      1.0f,
      false);

  solver.createTetBox(
      glm::vec3(0.0f, 40.0f, 10.0f),
      2.0f,
      glm::vec3(0.0f),
      1000.0f,
      1.0f,
      false);
  solver.createTetBox(
      glm::vec3(10.0f, 40.0f, 10.0f),
      2.0f,
      glm::vec3(0.0f),
      1000.0f,
      1.0f,
      false);
  solver.createTetBox(
      glm::vec3(-10.0f, 40.0f, 10.0f),
      2.0f,
      glm::vec3(0.0f),
      1000.0f,
      1.0f,
      false);

  solver.createSheet(glm::vec3(-20.0f, 20.0f, -20.0f), 1.0f, 1.0f, 1000000.0f);
}

SolverNode::SolverNode() {}

MStatus SolverNode::compute(const MPlug& plug, MDataBlock& data) {
  MStatus status = MStatus::kSuccess;

  if (plug == ouptutPositions) {
    MDataHandle outputHandle = data.outputValue(ouptutPositions, &status);
    McheckErr(status, "ERROR getting output positions handle.");

    MDataHandle outputPrevTimeHandle = data.outputValue(prevTime, &status);
    McheckErr(status, "ERROR getting output prevTime handle.");

    MDataHandle stepSizeHandle = data.inputValue(stepSize, &status);
    McheckErr(status, "ERROR getting step size handle.");

    MDataHandle itersHandle = data.inputValue(iterations, &status);
    McheckErr(status, "ERROR getting iterations handle.");

    MDataHandle timeHandle = data.inputValue(time, &status);
    McheckErr(status, "ERROR getting time handle.");

    MDataHandle prevTimeHandle = data.inputValue(prevTime, &status);
    McheckErr(status, "ERROR getting prevTime handle.");

    MDataHandle simulationEnabledHandle =
        data.inputValue(simulationEnabled, &status);
    McheckErr(status, "ERROR getting simulationEnabledHandle handle.");

    MDataHandle simulationStartTimeHandle =
        data.inputValue(simulationStartTime, &status);
    McheckErr(status, "ERROR getting simulationStartTimeHandle handle.");

    MArrayDataHandle meshArrayHandle = data.inputArrayValue(meshArray, &status);
    McheckErr(status, "ERROR getting meshArray handle.");

    MTime currentTime = timeHandle.asTime();
    const MTime previousTime = prevTimeHandle.asTime();

    double timeSeconds = currentTime.asUnits(MTime::Unit::kSeconds);

    bool simulationEnabledValue = simulationEnabledHandle.asBool();
    float currentStepSize = stepSizeHandle.asFloat();
    int currentIters = itersHandle.asInt();

    // The time check should technically not be necessary
    if (this->_resetSimulation) {
      this->_simulationTime = 0.0f;
      this->_resetSimulation = false;

      Pies::SolverOptions options{};
      options.gridSpacing = 1.0f;
      options.floorHeight = 0.0f;
      options.fixedTimestepSize = 1.0f / 60.0f;
      options.timeSubsteps = 1;

      this->_pSolver = std::make_unique<Pies::Solver>(options);
      setupTestScene(*this->_pSolver);

      std::vector<glm::vec3> newNodes;
      MPointArray pointArr;
      for (uint32_t meshIndex = 0; meshIndex < meshArrayHandle.elementCount();
           ++meshIndex) {
        MFnMesh mesh = meshArrayHandle.inputValue().asMesh();
        const MMatrix& transform = mesh.transformationMatrix();
        mesh.getPoints(pointArr);
        newNodes.resize(pointArr.length());

        for (uint32_t pointIndex = 0; pointIndex < pointArr.length();
             ++pointIndex) {
          MPoint point = transform * pointArr[pointIndex];
          newNodes[pointIndex] = glm::vec3(
              static_cast<float>(point.x),
              static_cast<float>(point.y),
              static_cast<float>(point.z));
        }

        this->_pSolver->addNodes(newNodes);

        pointArr.clear();
        newNodes.clear();
        meshArrayHandle.next();
      }
    }

    if (simulationEnabledValue && this->_simulationTime < timeSeconds) {
      this->_pSolver->tick(0.012f);
      this->_simulationTime = timeSeconds;

      MFnArrayAttrsData particlesAAD;
      MObject particlesObj = particlesAAD.create();
      MVectorArray particlesPosArr = particlesAAD.vectorArray("position");
      MVectorArray particlesScaleArr = particlesAAD.vectorArray("scale");
      MDoubleArray particlesIdArr = particlesAAD.doubleArray("id");

      const std::vector<Pies::Solver::Vertex> piesVerts =
          this->_pSolver->getVertices();

      particlesPosArr.setLength(piesVerts.size());
      particlesScaleArr.setLength(piesVerts.size());
      particlesIdArr.setLength(piesVerts.size());

      for (size_t i = 0; i < piesVerts.size(); ++i) {
        const Pies::Solver::Vertex& vertex = piesVerts[i];

        particlesIdArr[i] = static_cast<double>(i);
        particlesPosArr[i] =
            MVector(vertex.position.x, vertex.position.y, vertex.position.z);
        particlesScaleArr[i] =
            MVector(vertex.radius, vertex.radius, vertex.radius);
      }

      outputHandle.setMObject(particlesObj);
    }

    data.setClean(ouptutPositions);

    outputPrevTimeHandle.setMTime(currentTime);
    data.setClean(prevTime);

    data.setClean(plug);

    return MStatus::kSuccess;
  }

  return MStatus::kUnknownParameter;
}

void SolverNode::getCacheSetup(
    const MEvaluationNode& evaluationNode,
    MNodeCacheDisablingInfo& disablingInfo,
    MNodeCacheSetupInfo& setupInfo,
    MObjectArray& monitoredAttributes) const {
  // Check if the attribute is animated.  If it is animated, requestSimulation
  // will be true so that the behavior is the same as if the scene was played
  // back from the beginning of the time range.  That being said, it makes very
  // little sense to animate whether the simulation is active or not.
  const bool requestSimulation =
      MNodeCacheDisablingInfoHelper::testBooleanAttribute(
          nullptr,
          monitoredAttributes,
          evaluationNode,
          simulationEnabled,
          false);

  if (requestSimulation) {
    // Request simulation support
    setupInfo.setRequirement(MNodeCacheSetupInfo::kSimulationSupport, true);
    setupInfo.setPreference(MNodeCacheSetupInfo::kWantToCacheByDefault, true);
  }
}

void SolverNode::configCache(
    const MEvaluationNode& evalNode,
    MCacheSchema& schema) const {
  if (evalNode.dirtyPlugExists(ouptutPositions)) {
    schema.add(ouptutPositions);
  }
}

MTimeRange SolverNode::transformInvalidationRange(
    const MPlug& source,
    const MTimeRange& input) const {
  // Never call the parent class method, it's not meant to be called from the
  // base class, the implementation is only used to detect whether there is an
  // override or not.

  // When setting infinite values, do not use max()/min() directly, could
  // overflow easily.
  static constexpr MTime::MTick kMaximumTimeTick =
      std::numeric_limits<MTime::MTick>::max() / 2;
  static constexpr MTime::MTick kMinimumTimeTick =
      std::numeric_limits<MTime::MTick>::min() / 2 + 1;
  static const MTime kMaximumTime{
      kMaximumTimeTick / static_cast<double>(MTime::ticksPerSecond()),
      MTime::kSeconds};
  static const MTime kMinimumTime{
      kMinimumTimeTick / static_cast<double>(MTime::ticksPerSecond()),
      MTime::kSeconds};

  // Get the start time and whether simulation is enabled., but it should NOT be
  // animated, so it should be clean.
  MDataBlock data = const_cast<SolverNode*>(this)->forceCache();
  if (!data.isClean(simulationStartTime) || !data.isClean(simulationEnabled) ||
      !data.isClean(meshArray)) {
    this->_resetSimulation = true;
    return MTimeRange{kMinimumTime, kMaximumTime};
  }

  MStatus returnStatus;
  // Whether the simulation is enabled or not is clean, get its value.
  MDataHandle simulationEnabledData =
      data.inputValue(simulationEnabled, &returnStatus);
  const bool simulationEnabled =
      returnStatus ? simulationEnabledData.asBool() : true;
  if (!simulationEnabled) {
    return input;
  }

  // The start time is clean, get its value.
  MDataHandle simulationStartTimeData =
      data.inputValue(simulationStartTime, &returnStatus);
  const MTime simulationStartTime =
      returnStatus ? simulationStartTimeData.asTime() : kMinimumTime;
  const MTime simulationEndTime = kMaximumTime;

  if (input.intersects(simulationStartTime, simulationEndTime)) {
    // We invalidate the whole simulation time range, because something in the
    // range was invalidated, therefore we have to recompute the whole
    // simulation.
    this->_resetSimulation = true;
    return input | MTimeRange{simulationStartTime, simulationEndTime};
  } else {
    // Since the invalidation range does not touch the simulation,
    // We don't need to invalidate anything.
    return MTimeRange{};
  }
}

/*static*/
void* SolverNode::creator() { return new SolverNode(); }

/*static*/
MStatus SolverNode::initialize() {
  MStatus status = MStatus::kSuccess;

  MFnUnitAttribute timeAttr;
  SolverNode::time =
      timeAttr.create("time", "t", MFnUnitAttribute::kTime, 0.0, &status);
  McheckErr(status, "ERROR creating attribute SolverNode::time.");

  timeAttr.setWritable(true);
  timeAttr.setStorable(true);
  timeAttr.setReadable(true);
  timeAttr.setKeyable(true);
  status = addAttribute(SolverNode::time);
  McheckErr(status, "ERROR adding attribute SolverNode::time.");

  MFnUnitAttribute prevTimeAttr;
  SolverNode::prevTime =
      prevTimeAttr
          .create("prevTime", "prevT", MFnUnitAttribute::kTime, 0.0, &status);
  McheckErr(status, "ERROR creating attribute SolverNode::prevTime.");

  prevTimeAttr.setWritable(false);
  prevTimeAttr.setStorable(false);
  prevTimeAttr.setHidden(true);
  status = addAttribute(SolverNode::prevTime);
  McheckErr(status, "ERROR adding attribute SolverNode::prevTime.");

  MFnNumericAttribute stepSizeAttr;
  SolverNode::stepSize =
      stepSizeAttr
          .create("stepSize", "step", MFnNumericData::kFloat, 0.0, &status);
  McheckErr(status, "ERROR creating attribute SolverNode::stepSize.");

  stepSizeAttr.setWritable(true);
  status = addAttribute(SolverNode::stepSize);
  McheckErr(status, "ERROR adding attribute SolverNode::stepSize.");

  MFnNumericAttribute itersAttr;
  SolverNode::iterations =
      itersAttr.create("iterations", "iters", MFnNumericData::kInt, 0, &status);
  McheckErr(status, "ERROR creating attribute SolverNode::iterations.");

  itersAttr.setWritable(true);
  status = addAttribute(SolverNode::iterations);
  McheckErr(status, "ERROR adding attribute SolverNode::iterations.");

  MFnNumericAttribute simulationEnabledAttr;
  SolverNode::simulationEnabled = simulationEnabledAttr.create(
      "simulationEnabled",
      "enable",
      MFnNumericData::kBoolean,
      true,
      &status);
  McheckErr(status, "ERROR creating attribute SolverNode::simulationEnabled.");

  simulationEnabledAttr.setWritable(true);
  simulationEnabledAttr.setKeyable(false);
  status = addAttribute(SolverNode::simulationEnabled);
  McheckErr(status, "ERROR adding attribute SolverNode::simulationEnabled.");

  MFnNumericAttribute simulationStartTimeAttr;
  SolverNode::simulationStartTime = simulationStartTimeAttr.create(
      "simulationStartTime",
      "start",
      MFnNumericData::kFloat,
      0.0,
      &status);
  McheckErr(
      status,
      "ERROR creating attribute SolverNode::simulationStartTime.");

  simulationStartTimeAttr.setWritable(true);
  status = addAttribute(SolverNode::simulationStartTime);
  McheckErr(status, "ERROR adding attribute SolverNode::simulationStartTime.");

  MFnTypedAttribute meshArrayAttr;
  SolverNode::meshArray = meshArrayAttr.create(
      "meshArray",
      "meshes",
      MFnData::kMesh,
      MObject::kNullObj,
      &status);
  McheckErr(status, "ERROR creating attribute SolverNode::meshArray.");

  meshArrayAttr.setWritable(true);
  meshArrayAttr.setArray(true);
  meshArrayAttr.setUsesArrayDataBuilder(true); // ??
  meshArrayAttr.setIndexMatters(false);
  status = addAttribute(SolverNode::meshArray);
  McheckErr(status, "ERROR adding attribute SolverNode::meshArray.");

  MFnTypedAttribute outputPositionsAttr;
  SolverNode::ouptutPositions = outputPositionsAttr.create(
      "outputPositions",
      "out",
      MFnData::kDynArrayAttrs,
      MObject::kNullObj,
      &status);
  McheckErr(status, "ERROR creating attribute SolverNode::outputPositions.");

  outputPositionsAttr.setReadable(true);
  outputPositionsAttr.setWritable(false);
  outputPositionsAttr.setStorable(true);
  outputPositionsAttr.setCached(true);
  status = addAttribute(SolverNode::ouptutPositions);
  McheckErr(status, "ERROR adding attribute SolverNode::ouptutPositions.");

  // Setup input-output dependencies
  status = attributeAffects(time, ouptutPositions);
  McheckErr(status, "ERROR attributeAffects(time, ouptutPositions).");

  status = attributeAffects(stepSize, ouptutPositions);
  McheckErr(status, "ERROR attributeAffects(stepSize, ouptutPositions).");

  status = attributeAffects(iterations, ouptutPositions);
  McheckErr(status, "ERROR attributeAffects(iterations, ouptutPositions).");

  status = attributeAffects(simulationEnabled, ouptutPositions);
  McheckErr(
      status,
      "ERROR attributeAffects(simulationEnabled, ouptutPositions).");

  status = attributeAffects(simulationStartTime, ouptutPositions);
  McheckErr(
      status,
      "ERROR attributeAffects(simulationStartTime, ouptutPositions).");

  status = attributeAffects(meshArray, ouptutPositions);
  McheckErr(
      status,
      "ERROR attributeAffects(simulationStartTime, ouptutPositions).");

  return status;
}