#include "SolverNode.h"

#include <maya/MStatus.h>
#include <maya/MIOStream.h>
#include <maya/MDataHandle.h>
#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MIntArray.h>
#include <maya/MFnStringData.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MGlobal.h>

#include <Pies/Solver.h>

#include <vector>
#include <string>
#include <iostream>

#define McheckErr(stat, msg)   \
  if (!stat)                   \
  {                            \
    MGlobal::displayInfo(msg); \
    stat.perror(msg);          \
    return stat;               \
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
MObject SolverNode::angle;

/*static*/
MObject SolverNode::iterations;

/*static*/
MObject SolverNode::grammar;

/*static*/
MObject SolverNode::outputMesh;

MStatus SolverNode::compute(const MPlug &plug, MDataBlock &data)
{
  MStatus status = MStatus::kSuccess;

  Pies::SolverOptions solverOptions{};
  Pies::Solver solver(solverOptions);
  solver.createTetBox(glm::vec3(0.0f), 1.0f, glm::vec3(0.0f), 1000.0f);

  if (plug == outputMesh)
  {
    MDataHandle outputHandle = data.outputValue(outputMesh, &status);
    McheckErr(status, "ERROR getting output mesh handle.");

    MDataHandle stepSizeHandle = data.inputValue(stepSize, &status);
    McheckErr(status, "ERROR getting step size handle.");

    MDataHandle angleHandle = data.inputValue(angle, &status);
    McheckErr(status, "ERROR getting angle handle.");

    MDataHandle itersHandle = data.inputValue(iterations, &status);
    McheckErr(status, "ERROR getting iterations handle.");

    MDataHandle grammarHandle = data.inputValue(grammar, &status);
    McheckErr(status, "ERROR getting grammar handle.");

    MDataHandle timeHandle = data.inputValue(time, &status);
    McheckErr(status, "ERROR getting time handle.");

    double timeSeconds = timeHandle.asTime().asUnits(MTime::Unit::kSeconds);

    float currentStepSize = stepSizeHandle.asFloat();
    double currentAngle = angleHandle.asDouble();
    int currentIters = itersHandle.asInt();
    MString currentGrammar = grammarHandle.asString();

    MFnMeshData meshCreator;
    MObject meshData = meshCreator.create(&status);
    McheckErr(status, "ERROR creating output mesh.");

    MPointArray points;
    MIntArray faceCounts;
    MIntArray faceConnects;

    MFnMesh mesh;
    mesh.create(points.length(), faceCounts.length(), points, faceCounts, faceConnects, meshData, &status);
    McheckErr(status, "ERROR reconstructing mesh.");

    outputHandle.set(meshData);
    data.setClean(plug);

    return MStatus::kSuccess;
  }

  return MStatus::kUnknownParameter;
}

/*static*/
void *SolverNode::creator()
{
  return new SolverNode();
}

/*static*/
MStatus SolverNode::initialize()
{
  MStatus status = MStatus::kSuccess;

  MFnUnitAttribute timeAttr;
  SolverNode::time = timeAttr.create("time", "t", MFnUnitAttribute::kTime, 0.0, &status);
  McheckErr(status, "ERROR creating attribute SolverNode::time.");

  timeAttr.setWritable(true);
  status = addAttribute(SolverNode::time);
  McheckErr(status, "ERROR adding attribute SolverNode::time.");

  MFnTypedAttribute grammarAttr;
  SolverNode::grammar = grammarAttr.create("grammar", "grm", MFnStringData::kString, MObject::kNullObj, &status);
  McheckErr(status, "ERROR creating attribute SolverNode::grammar.");

  grammarAttr.setWritable(true);
  status = addAttribute(SolverNode::grammar);
  McheckErr(status, "ERROR adding attribute SolverNode::grammar.");

  MFnNumericAttribute stepSizeAttr;
  SolverNode::stepSize = stepSizeAttr.create("stepSize", "step", MFnNumericData::kFloat, 0.0, &status);
  McheckErr(status, "ERROR creating attribute SolverNode::stepSize.");

  stepSizeAttr.setWritable(true);
  status = addAttribute(SolverNode::stepSize);
  McheckErr(status, "ERROR adding attribute SolverNode::stepSize.");

  MFnNumericAttribute angleAttr;
  SolverNode::angle = angleAttr.create("angle", "ang", MFnNumericData::kDouble, 0.0, &status);
  McheckErr(status, "ERROR creating attribute SolverNode::angle.");

  angleAttr.setWritable(true);
  status = addAttribute(SolverNode::angle);
  McheckErr(status, "ERROR adding attribute SolverNode::angle.");

  MFnNumericAttribute itersAttr;
  SolverNode::iterations = itersAttr.create("iterations", "iters", MFnNumericData::kInt, 0, &status);
  McheckErr(status, "ERROR creating attribute SolverNode::iterations.");

  itersAttr.setWritable(true);
  status = addAttribute(SolverNode::iterations);
  McheckErr(status, "ERROR adding attribute SolverNode::iterations.");

  MFnTypedAttribute outputMeshAttr;
  SolverNode::outputMesh = outputMeshAttr.create("outputMesh", "out", MFnMeshData::kMesh, MObject::kNullObj, &status);
  McheckErr(status, "ERROR creating attribute SolverNode::outputMesh.");

  outputMeshAttr.setReadable(true);
  outputMeshAttr.setWritable(false);
  outputMeshAttr.setStorable(true);
  status = addAttribute(SolverNode::outputMesh);
  McheckErr(status, "ERROR adding attribute SolverNode::outputMesh.");

  // Setup input-output dependencies
  status = attributeAffects(time, outputMesh);
  McheckErr(status, "ERROR attributeAffects(time, outputMesh).");

  status = attributeAffects(stepSize, outputMesh);
  McheckErr(status, "ERROR attributeAffects(stepSize, outputMesh).");

  status = attributeAffects(angle, outputMesh);
  McheckErr(status, "ERROR attributeAffects(angle, outputMesh).");

  status = attributeAffects(grammar, outputMesh);
  McheckErr(status, "ERROR attributeAffects(grammar, outputMesh).");

  status = attributeAffects(iterations, outputMesh);
  McheckErr(status, "ERROR attributeAffects(iterations, outputMesh).");

  return status;
}