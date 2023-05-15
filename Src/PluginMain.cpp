#include "PiesSoftBodyNode.h"
#include "SolverNode.h"

#include <cuda.h>

#include <maya/MArgList.h>
#include <maya/MDGModifier.h>
#include <maya/MDoubleArray.h>
#include <maya/MFnDependencyNode.h>
#include <maya/MFnNurbsCurve.h>
#include <maya/MFnPlugin.h>
#include <maya/MGlobal.h>
#include <maya/MIOStream.h>
#include <maya/MPlugArray.h>
#include <maya/MPoint.h>
#include <maya/MPointArray.h>
#include <maya/MPxCommand.h>
#include <maya/MSimple.h>
#include <maya/MStreamUtils.h>
#include <maya/MString.h>
#include <maya/MStringArray.h>
#include <maya/MVector.h>

#include <iostream>
#include <list>

static CUcontext PiesCudaContext;

MStatus initializePlugin(MObject obj) {
  MStatus status = MStatus::kSuccess;
  MFnPlugin plugin(obj, "PiesForMayaPlugin", "1.0", "Any");

  SolverNode::pluginDir = plugin.loadPath(&status);
  if (!status) {
    status.perror("loadPath");
    return status;
  }

  // TODO: Eventually implement persistent threaded cuda context...
  // See:
  // https://www.nvidia.com.tw/content/apacevents/siggraph-asia-2012/developing-an-optimized-maya-plugin-using-cuda-and-opengl-WBraithwaite.pdf
  CUdevice dev;
  // TODO: Double check this is picking the discrete GPU
  cuDeviceGet(&dev, 0);
  cuCtxCreate(&PiesCudaContext, 0, dev);
  cuCtxPopCurrent(0);

  SolverNode::PiesCudaContext = PiesCudaContext;

  std::string pluginPath(SolverNode::pluginDir.asChar());

  std::string scriptPath = pluginPath + "/../../Scripts/Pies.mel";
  std::string addScriptPathCmd = "source \"" + scriptPath + "\";";

  MGlobal::executeCommand(MString(addScriptPathCmd.c_str()));

  status = plugin.registerNode(
      "SolverNode",
      SolverNode::id,
      SolverNode::creator,
      SolverNode::initialize);
  if (!status) {
    status.perror("registerNode");
    return status;
  }

  status = plugin.registerNode(
      "PiesSoftBodyNode",
      PiesSoftBodyNode::id,
      PiesSoftBodyNode::creator,
      PiesSoftBodyNode::initialize);
  if (!status) {
    status.perror("registerNode");
    return status;
  }

  return status;
}

MStatus uninitializePlugin(MObject obj) {
  MStatus status = MStatus::kSuccess;
  MFnPlugin plugin(obj);

  cuCtxPushCurrent(PiesCudaContext);
  cuCtxDestroy(PiesCudaContext);

  status = plugin.deregisterNode(SolverNode::id);
  if (!status) {
    status.perror("deregisterNode");
    return status;
  }

  status = plugin.deregisterNode(PiesSoftBodyNode::id);
  if (!status) {
    status.perror("deregisterNode");
    return status;
  }

  return status;
}
