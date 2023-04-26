#include "SolverNode.h"

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

MStatus initializePlugin(MObject obj) {
  MStatus status = MStatus::kSuccess;
  MFnPlugin plugin(obj, "PiesForMayaPlugin", "1.0", "Any");

  SolverNode::pluginDir = plugin.loadPath(&status);
  if (!status) {
    status.perror("loadPath");
    return status;
  }

  std::string pluginPath(SolverNode::pluginDir.asChar());
  // TODO: Add scripts from folder

  std::string scriptPath = pluginPath + "/../PiesMenu.mel";
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

  return status;
}

MStatus uninitializePlugin(MObject obj) {
  MStatus status = MStatus::kSuccess;
  MFnPlugin plugin(obj);

  status = plugin.deregisterNode(SolverNode::id);
  if (!status) {
    status.perror("deregisterNode");
    return status;
  }

  return status;
}
