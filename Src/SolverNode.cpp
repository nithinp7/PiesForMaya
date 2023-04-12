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
MObject SolverNode::ouptutPositions;

SolverNode::SolverNode() {
  Pies::SolverOptions options{};
  options.gridSpacing = 1.0f;
  options.floorHeight = -8.0f;

  this->_pSolver = std::make_unique<Pies::Solver>(options);
  this->_pSolver->createTetBox(
      glm::vec3(0.0f),
      1.0f,
      glm::vec3(0.0f),
      1000.0f,
      1.0f,
      true);
}

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

    MTime currentTime = timeHandle.asTime();
    const MTime previousTime = prevTimeHandle.asTime();

    // if (currentTime == previousTime) {
    //   // TODO: ??
    //   data.setClean(ouptutPositions);
    //   return MStatus::kSuccess;
    // }

    double timeSeconds = currentTime.asUnits(MTime::Unit::kSeconds);

    bool simulationEnabledValue = simulationEnabledHandle.asBool();
    float currentStepSize = stepSizeHandle.asFloat();
    int currentIters = itersHandle.asInt();

    if (simulationEnabledValue && this->_simulationTime < timeSeconds) {
      // TODO: Is this correct??
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
  if (!data.isClean(simulationStartTime) || !data.isClean(simulationEnabled)) {
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

  return status;
}