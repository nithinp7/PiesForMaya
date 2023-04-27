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
#include <maya/MMatrix.h>
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
MObject SolverNode::collisionIterations;

/*static*/
MObject SolverNode::collisionDistance;

/*static*/
MObject SolverNode::collisionThickness;

/*static*/
MObject SolverNode::gridSpacing;

/*static*/
MObject SolverNode::gravity;

/*static*/
MObject SolverNode::damping;

/*static*/
MObject SolverNode::friction;

/*static*/
MObject SolverNode::floorHeight;

/*static*/
MObject SolverNode::threadCount;

/*static*/
MObject SolverNode::simulationEnabled;

/*static*/
MObject SolverNode::simulationStartTime;

/*static*/
MObject SolverNode::meshArray;

/*static*/
MObject SolverNode::outputPositions;

/*static*/
MObject SolverNode::outputMesh;

static void setupTestScene(Pies::Solver& solver) {
  // solver.createTetBox(
  //     glm::vec3(0.0f, 30.0f, 0.0f),
  //     1.5f,
  //     glm::vec3(0.0f),
  //     1000.0f,
  //     1.0f,
  //     true);

  // solver.createTetBox(
  //     glm::vec3(0.0f, 40.0f, 0.0f),
  //     2.0f,
  //     glm::vec3(0.0f),
  //     1000.0f,
  //     1.0f,
  //     false);
  // solver.createTetBox(
  //     glm::vec3(10.0f, 40.0f, 0.0f),
  //     2.0f,
  //     glm::vec3(0.0f),
  //     1000.0f,
  //     1.0f,
  //     false);
  // solver.createTetBox(
  //     glm::vec3(-10.0f, 40.0f, 0.0f),
  //     2.0f,
  //     glm::vec3(0.0f),
  //     1000.0f,
  //     1.0f,
  //     false);

  // solver.createTetBox(
  //     glm::vec3(0.0f, 40.0f, 10.0f),
  //     2.0f,
  //     glm::vec3(0.0f),
  //     1000.0f,
  //     1.0f,
  //     false);
  // solver.createTetBox(
  //     glm::vec3(10.0f, 40.0f, 10.0f),
  //     2.0f,
  //     glm::vec3(0.0f),
  //     1000.0f,
  //     1.0f,
  //     false);
  // solver.createTetBox(
  //     glm::vec3(-10.0f, 40.0f, 10.0f),
  //     2.0f,
  //     glm::vec3(0.0f),
  //     1000.0f,
  //     1.0f,
  //     false);

  solver.createSheet(glm::vec3(0.0f, 20.0f, 0.0f), 2.0f, 1.0f, 800.0f);
}

SolverNode::SolverNode() {}

MStatus SolverNode::compute(const MPlug& plug, MDataBlock& data) {
  MStatus status = MStatus::kSuccess;

  if (plug == outputPositions || plug == outputMesh) {
    MDataHandle outputHandle = data.outputValue(outputPositions, &status);
    McheckErr(status, "ERROR getting output positions handle.");

    MDataHandle outMeshHandle = data.outputValue(outputMesh, &status);
    McheckErr(status, "ERROR getting output mesh handle.");

    MDataHandle outputPrevTimeHandle = data.outputValue(prevTime, &status);
    McheckErr(status, "ERROR getting output prevTime handle.");

    MDataHandle stepSizeHandle = data.inputValue(stepSize, &status);
    McheckErr(status, "ERROR getting step size handle.");

    MDataHandle itersHandle = data.inputValue(iterations, &status);
    McheckErr(status, "ERROR getting iterations handle.");

    MDataHandle colItersHandle = data.inputValue(collisionIterations, &status);
    McheckErr(status, "ERROR getting collision iterations handle.");

    MDataHandle colDistHandle = data.inputValue(collisionDistance, &status);
    McheckErr(status, "ERROR getting collision distance handle.");

    MDataHandle colThicknessHandle = data.inputValue(collisionThickness, &status);
    McheckErr(status, "ERROR getting collision thickness handle.");

    MDataHandle gridHandle = data.inputValue(gridSpacing, &status);
    McheckErr(status, "ERROR getting grid spacing handle.");

    MDataHandle gravityHandle = data.inputValue(gravity, &status);
    McheckErr(status, "ERROR getting gravity handle.");
    
    MDataHandle dampingHandle = data.inputValue(damping, &status);
    McheckErr(status, "ERROR getting damping handle.");

    MDataHandle frictionHandle = data.inputValue(friction, &status);
    McheckErr(status, "ERROR getting friction handle.");

    MDataHandle floorHeightHandle = data.inputValue(floorHeight, &status);
    McheckErr(status, "ERROR getting floorHeight handle.");

    MDataHandle threadCountHandle = data.inputValue(threadCount, &status);
    McheckErr(status, "ERROR getting threadCount handle.");

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
      options.fixedTimestepSize = stepSizeHandle.asFloat();
      options.iterations = itersHandle.asInt();
      options.collisionStabilizationIterations = colItersHandle.asInt();
      options.collisionThresholdDistance = colDistHandle.asFloat();
      options.collisionThickness = colThicknessHandle.asFloat();
      options.gridSpacing = gridHandle.asFloat();
      options.gravity = gravityHandle.asFloat();
      options.damping = dampingHandle.asFloat();
      options.friction = frictionHandle.asFloat();
      options.floorHeight = floorHeightHandle.asFloat();
      options.threadCount = threadCountHandle.asInt();

      this->_pSolver = std::make_unique<Pies::Solver>(options);
      setupTestScene(*this->_pSolver);

      std::vector<glm::vec3> vertices;
      std::vector<uint32_t> indices;
      MPointArray pointArr;
      MIntArray triCountPerPoly;
      MIntArray triIndices;
      for (uint32_t meshIndex = 0; meshIndex < meshArrayHandle.elementCount();
           ++meshIndex) {
        MDataHandle meshHandle = meshArrayHandle.inputValue();
        const MMatrix& transform = meshHandle.geometryTransformMatrix();
        MFnMesh mesh = meshHandle.asMesh();
        mesh.getPoints(pointArr, MSpace::kWorld);
        mesh.getTriangles(triCountPerPoly, triIndices);
        vertices.resize(pointArr.length());
        indices.resize(triIndices.length());

        for (uint32_t pointIndex = 0; pointIndex < pointArr.length();
             ++pointIndex) {
          MPoint point = pointArr[pointIndex];
          vertices[pointIndex] = glm::vec3(
              static_cast<float>(point.x),
              static_cast<float>(point.y),
              static_cast<float>(point.z));
        }

        for (uint32_t i = 0; i < triIndices.length(); i += 3) {
          indices[i] = static_cast<uint32_t>(triIndices[i]);
          indices[i + 1] = static_cast<uint32_t>(triIndices[i + 1]);
          indices[i + 2] = static_cast<uint32_t>(triIndices[i + 2]);
        }

        // this->_pSolver->addNodes(vertices);
        // this->_pSolver->addTriMeshVolume(vertices, indices, 1000.0f);
        this->_pSolver->addClothMesh(vertices, indices, 5000.0f);

        pointArr.clear();
        triCountPerPoly.clear();
        triIndices.clear();
        vertices.clear();
        indices.clear();
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
      const std::vector<Pies::Triangle> piesTris =
          this->_pSolver->getTriangles();

      particlesPosArr.setLength(piesVerts.size());
      particlesScaleArr.setLength(piesVerts.size());
      particlesIdArr.setLength(piesVerts.size());

      MFnMeshData meshCreator;
      MObject meshData = meshCreator.create(&status);
      McheckErr(status, "ERROR creating output mesh.");

      MPointArray points;
      MIntArray faceCounts;
      MIntArray faceConnects;

      points.setLength(piesVerts.size());
      faceCounts.setLength(piesTris.size());
      faceConnects.setLength(3 * piesTris.size());

      for (size_t i = 0; i < piesVerts.size(); ++i) {
        const Pies::Solver::Vertex& vertex = piesVerts[i];

        particlesIdArr[i] = static_cast<double>(i);
        particlesPosArr[i] =
            MVector(vertex.position.x, vertex.position.y, vertex.position.z);
        particlesScaleArr[i] =
            MVector(vertex.radius, vertex.radius, vertex.radius);

        points[i] =
            MPoint(vertex.position.x, vertex.position.y, vertex.position.z);
      }

      for (size_t i = 0; i < piesTris.size(); ++i) {
        const Pies::Triangle& triangle = piesTris[i];

        faceCounts[i] = 3;
        faceConnects[3 * i] = triangle.nodeIds[0];
        faceConnects[3 * i + 1] = triangle.nodeIds[1];
        faceConnects[3 * i + 2] = triangle.nodeIds[2];
      }

      MFnMesh mesh;
      mesh.create(
          points.length(),
          faceCounts.length(),
          points,
          faceCounts,
          faceConnects,
          meshData,
          &status);
		  McheckErr(status, "ERROR reconstructing mesh.");

      outMeshHandle.set(meshData);
      outputHandle.setMObject(particlesObj);
    }

    data.setClean(outputPositions);
    data.setClean(outputMesh);

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
  if (evalNode.dirtyPlugExists(outputMesh)) {
    schema.add(outputMesh);
  }

  if (evalNode.dirtyPlugExists(outputPositions)) {
    schema.add(outputPositions);
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
  if (!data.isClean(simulationStartTime) || 
      !data.isClean(simulationEnabled) ||
      !data.isClean(meshArray) ||
      !data.isClean(stepSize) ||
      !data.isClean(iterations) ||
      !data.isClean(collisionIterations) ||
      !data.isClean(collisionDistance) ||
      !data.isClean(collisionThickness) ||
      !data.isClean(gridSpacing) ||
      !data.isClean(gravity) ||
      !data.isClean(damping) ||
      !data.isClean(friction) ||
      !data.isClean(floorHeight) ||
      !data.isClean(threadCount)) {
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
          .create("stepSize", "step", MFnNumericData::kFloat, 0.012, &status);
  McheckErr(status, "ERROR creating attribute SolverNode::stepSize.");

  stepSizeAttr.setWritable(true);
  status = addAttribute(SolverNode::stepSize);
  McheckErr(status, "ERROR adding attribute SolverNode::stepSize.");

  MFnNumericAttribute itersAttr;
  SolverNode::iterations =
      itersAttr.create("iterations", "iters", MFnNumericData::kInt, 4, &status);
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

  
  MFnNumericAttribute collisionItersAttr;
  SolverNode::collisionIterations =
      collisionItersAttr.create("collisionIterations", "cIters", MFnNumericData::kInt, 4, &status);
  McheckErr(status, "ERROR creating attribute SolverNode::collisionIterations.");

  collisionItersAttr.setWritable(true);
  status = addAttribute(SolverNode::collisionIterations);
  McheckErr(status, "ERROR adding attribute SolverNode::collisionIterations.");

  MFnNumericAttribute collisionDistanceAttr;
  SolverNode::collisionDistance =
      collisionDistanceAttr
          .create("collisionDistance", "cDist", MFnNumericData::kFloat, 0.1, &status);
  McheckErr(status, "ERROR creating attribute SolverNode::collisionDistance.");

  collisionDistanceAttr.setWritable(true);
  status = addAttribute(SolverNode::collisionDistance);
  McheckErr(status, "ERROR adding attribute SolverNode::collisionDistance.");

  MFnNumericAttribute collisionThicknessAttr;
  SolverNode::collisionThickness =
      collisionThicknessAttr
          .create("collisionThickness", "cThick", MFnNumericData::kFloat, 0.05, &status);
  McheckErr(status, "ERROR creating attribute SolverNode::collisionThickness.");

  collisionThicknessAttr.setWritable(true);
  status = addAttribute(SolverNode::collisionThickness);
  McheckErr(status, "ERROR adding attribute SolverNode::collisionThickness.");

  MFnNumericAttribute gridSpacingAttr;
  SolverNode::gridSpacing =
      gridSpacingAttr
          .create("gridSpacing", "grid", MFnNumericData::kFloat, 2.0, &status);
  McheckErr(status, "ERROR creating attribute SolverNode::gridSpacing.");

  gridSpacingAttr.setWritable(true);
  status = addAttribute(SolverNode::gridSpacing);
  McheckErr(status, "ERROR adding attribute SolverNode::gridSpacing.");
  
  MFnNumericAttribute gravityAttr;
  SolverNode::gravity =
      gravityAttr
          .create("gravity", "g", MFnNumericData::kFloat, 10.0, &status);
  McheckErr(status, "ERROR creating attribute SolverNode::gravity.");

  gravityAttr.setWritable(true);
  status = addAttribute(SolverNode::gravity);
  McheckErr(status, "ERROR adding attribute SolverNode::gravity.");

  MFnNumericAttribute dampingAttr;
  SolverNode::damping =
      dampingAttr
          .create("damping", "d", MFnNumericData::kFloat, 0.006, &status);
  McheckErr(status, "ERROR creating attribute SolverNode::damping.");

  dampingAttr.setWritable(true);
  status = addAttribute(SolverNode::damping);
  McheckErr(status, "ERROR adding attribute SolverNode::damping.");

  MFnNumericAttribute frictionAttr;
  SolverNode::friction =
      frictionAttr
          .create("friction", "f", MFnNumericData::kFloat, 0.01, &status);
  McheckErr(status, "ERROR creating attribute SolverNode::friction.");

  frictionAttr.setWritable(true);
  status = addAttribute(SolverNode::friction);
  McheckErr(status, "ERROR adding attribute SolverNode::friction.");

  MFnNumericAttribute floorAttr;
  SolverNode::floorHeight =
      floorAttr
          .create("floorHeight", "floor", MFnNumericData::kFloat, 0.0, &status);
  McheckErr(status, "ERROR creating attribute SolverNode::floorHeight.");

  floorAttr.setWritable(true);
  status = addAttribute(SolverNode::floorHeight);
  McheckErr(status, "ERROR adding attribute SolverNode::floorHeight.");

  MFnNumericAttribute threadCountAttr;
  SolverNode::threadCount = 
      threadCountAttr.create("threadCount", "threads", MFnNumericData::kInt, 8, &status);
  McheckErr(status, "ERROR creating attribute SolverNode::threadCount.");

  threadCountAttr.setWritable(true);
  status = addAttribute(SolverNode::threadCount);
  McheckErr(status, "ERROR adding attribute SolverNode::threadCount.");

  MFnTypedAttribute outputPositionsAttr;
  SolverNode::outputPositions = outputPositionsAttr.create(
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
  status = addAttribute(SolverNode::outputPositions);
  McheckErr(status, "ERROR adding attribute SolverNode::outputPositions.");

  MFnTypedAttribute outputMeshAttr;
  SolverNode::outputMesh = outputMeshAttr.create(
      "outputMesh",
      "outMesh",
      MFnData::kMesh,
      MObject::kNullObj,
      &status);
  McheckErr(status, "ERROR creating attribute SolverNode::outputMesh.");

  outputMeshAttr.setReadable(true);
  outputMeshAttr.setWritable(false);
  outputMeshAttr.setStorable(true);
  outputMeshAttr.setCached(true);
  status = addAttribute(SolverNode::outputMesh);
  McheckErr(status, "ERROR adding attribute SolverNode::outputMesh.");

  // Setup input-output dependencies
  status = attributeAffects(time, outputPositions);
  McheckErr(status, "ERROR attributeAffects(time, outputPositions).");

  status = attributeAffects(stepSize, outputPositions);
  McheckErr(status, "ERROR attributeAffects(stepSize, outputPositions).");

  status = attributeAffects(iterations, outputPositions);
  McheckErr(status, "ERROR attributeAffects(iterations, outputPositions).");

  status = attributeAffects(simulationEnabled, outputPositions);
  McheckErr(
      status,
      "ERROR attributeAffects(simulationEnabled, outputPositions).");

  status = attributeAffects(simulationStartTime, outputPositions);
  McheckErr(
      status,
      "ERROR attributeAffects(simulationStartTime, outputPositions).");

  status = attributeAffects(meshArray, outputPositions);
  McheckErr(
      status,
      "ERROR attributeAffects(simulationStartTime, outputPositions).");

  status = attributeAffects(collisionIterations, outputPositions);
  McheckErr(status, "ERROR attributeAffects(collisionIterations, outputPositions).");

  status = attributeAffects(collisionDistance, outputPositions);
  McheckErr(status, "ERROR attributeAffects(collisionDistance, outputPositions).");

  status = attributeAffects(collisionThickness, outputPositions);
  McheckErr(status, "ERROR attributeAffects(collisionThickness, outputPositions).");

  status = attributeAffects(gridSpacing, outputPositions);
  McheckErr(
      status,
      "ERROR attributeAffects(gridSpacing, outputPositions).");

  status = attributeAffects(gravity, outputPositions);
  McheckErr(
      status,
      "ERROR attributeAffects(gravity, outputPositions).");

  status = attributeAffects(damping, outputPositions);
  McheckErr(
      status,
      "ERROR attributeAffects(damping, outputPositions).");

  status = attributeAffects(friction, outputPositions);
  McheckErr(
      status,
      "ERROR attributeAffects(friction, outputPositions).");

  status = attributeAffects(floorHeight, outputPositions);
  McheckErr(
      status,
      "ERROR attributeAffects(floorHeight, outputPositions).");
      
  status = attributeAffects(threadCount, outputPositions);
  McheckErr(
      status,
      "ERROR attributeAffects(threadCount, outputPositions).");

  // Setup input-output dependencies
  status = attributeAffects(time, outputMesh);
  McheckErr(status, "ERROR attributeAffects(time, outputMesh).");

  status = attributeAffects(stepSize, outputMesh);
  McheckErr(status, "ERROR attributeAffects(stepSize, outputMesh).");

  status = attributeAffects(iterations, outputMesh);
  McheckErr(status, "ERROR attributeAffects(iterations, outputMesh).");

  status = attributeAffects(simulationEnabled, outputMesh);
  McheckErr(status, "ERROR attributeAffects(simulationEnabled, outputMesh).");

  status = attributeAffects(simulationStartTime, outputMesh);
  McheckErr(status, "ERROR attributeAffects(simulationStartTime, outputMesh).");

  status = attributeAffects(meshArray, outputMesh);
  McheckErr(status, "ERROR attributeAffects(simulationStartTime, outputMesh).");
  
  status = attributeAffects(collisionIterations, outputMesh);
  McheckErr(status, "ERROR attributeAffects(collisionIterations, outputMesh).");

  status = attributeAffects(collisionDistance, outputMesh);
  McheckErr(status, "ERROR attributeAffects(collisionDistance, outputMesh).");

  status = attributeAffects(collisionThickness, outputMesh);
  McheckErr(status, "ERROR attributeAffects(collisionThickness, outputMesh).");

  status = attributeAffects(gridSpacing, outputMesh);
  McheckErr(
      status,
      "ERROR attributeAffects(gridSpacing, outputMesh).");

  status = attributeAffects(gravity, outputMesh);
  McheckErr(
      status,
      "ERROR attributeAffects(gravity, outputMesh).");

  status = attributeAffects(damping, outputMesh);
  McheckErr(
      status,
      "ERROR attributeAffects(damping, outputMesh).");

  status = attributeAffects(friction, outputMesh);
  McheckErr(
      status,
      "ERROR attributeAffects(friction, outputMesh).");

  status = attributeAffects(floorHeight, outputMesh);
  McheckErr(
      status,
      "ERROR attributeAffects(floorHeight, outputMesh).");
      
  status = attributeAffects(threadCount, outputMesh);
  McheckErr(
      status,
      "ERROR attributeAffects(threadCount, outputMesh).");
  return status;
}