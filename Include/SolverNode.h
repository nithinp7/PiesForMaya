#pragma once

#include <Pies/Solver.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MObject.h>
#include <maya/MPxNode.h>
#include <maya/MStatus.h>
#include <maya/MString.h>
#include <maya/MTypeId.h>

#include <memory>
#include <string>

// TODO: Rename to PiesSolvernode
class SolverNode : public MPxNode {
public:
  SolverNode();
  MStatus compute(const MPlug& plug, MDataBlock& data) override;
  void getCacheSetup(
      const MEvaluationNode& evaluationNode,
      MNodeCacheDisablingInfo& disablingInfo,
      MNodeCacheSetupInfo& setupInfo,
      MObjectArray& monitoredAttributes) const override;
  void configCache(const MEvaluationNode& evalNode, MCacheSchema& schema) const;
  MTimeRange transformInvalidationRange(
      const MPlug& source,
      const MTimeRange& input) const override;

  static void* creator();
  static MStatus initialize();

  static MString pluginDir;
  static MTypeId id;
  static MObject time;
  static MObject prevTime;
  static MObject stepSize;
  static MObject iterations;
  static MObject simulationEnabled;
  static MObject simulationStartTime;

  static MObject meshArray;

  static MObject ouptutPositions;

private:
  std::unique_ptr<Pies::Solver> _pSolver;
  mutable bool _resetSimulation = true;
  float _simulationTime = 0.0f;
};
