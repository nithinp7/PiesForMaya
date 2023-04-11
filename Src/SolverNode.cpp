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
MObject SolverNode::stepSize;

/*static*/
MObject SolverNode::iterations;

/*static*/
MObject SolverNode::ouptutPositions;

/*static*/
std::unique_ptr<Pies::Solver> SolverNode::_pSolver = nullptr;

SolverNode::SolverNode() {
  if (!this->_pSolver) {
    Pies::SolverOptions options{};
    options.gridSpacing = 1.0f;
    options.floorHeight = -8.0f;

    SolverNode::_pSolver = std::make_unique<Pies::Solver>(options);
    SolverNode::_pSolver->createTetBox(
        glm::vec3(0.0f),
        1.0f,
        glm::vec3(0.0f),
        1000.0f,
        1.0f,
        true);
  }
}

MStatus SolverNode::compute(const MPlug& plug, MDataBlock& data) {
  MStatus status = MStatus::kSuccess;

  // TODO: this is a hack
  SolverNode::_pSolver->tick(0.012f);

  if (plug == ouptutPositions) {
    MDataHandle outputHandle = data.outputValue(ouptutPositions, &status);
    McheckErr(status, "ERROR getting output positions handle.");

    MDataHandle stepSizeHandle = data.inputValue(stepSize, &status);
    McheckErr(status, "ERROR getting step size handle.");

    MDataHandle itersHandle = data.inputValue(iterations, &status);
    McheckErr(status, "ERROR getting iterations handle.");

    MDataHandle timeHandle = data.inputValue(time, &status);
    McheckErr(status, "ERROR getting time handle.");

    double timeSeconds = timeHandle.asTime().asUnits(MTime::Unit::kSeconds);

    float currentStepSize = stepSizeHandle.asFloat();
    int currentIters = itersHandle.asInt();

    MFnArrayAttrsData particlesAAD;
    MObject particlesObj = particlesAAD.create();
    MVectorArray particlesPosArr = particlesAAD.vectorArray("position");
    MVectorArray particlesScaleArr = particlesAAD.vectorArray("scale");
    MDoubleArray particlesIdArr = particlesAAD.doubleArray("id");

    const std::vector<Pies::Solver::Vertex> piesVerts =
        SolverNode::_pSolver->getVertices();

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

    data.setClean(plug);

    return MStatus::kSuccess;
  }

  return MStatus::kUnknownParameter;
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
  status = addAttribute(SolverNode::time);
  McheckErr(status, "ERROR adding attribute SolverNode::time.");

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
  status = addAttribute(SolverNode::ouptutPositions);
  McheckErr(status, "ERROR adding attribute SolverNode::ouptutPositions.");

  // Setup input-output dependencies
  status = attributeAffects(time, ouptutPositions);
  McheckErr(status, "ERROR attributeAffects(time, ouptutPositions).");

  status = attributeAffects(stepSize, ouptutPositions);
  McheckErr(status, "ERROR attributeAffects(stepSize, ouptutPositions).");

  status = attributeAffects(iterations, ouptutPositions);
  McheckErr(status, "ERROR attributeAffects(iterations, ouptutPositions).");

  return status;
}