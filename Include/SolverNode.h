#pragma once

#include <Pies/Solver.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MObject.h>
#include <maya/MPxNode.h>
#include <maya/MStatus.h>
#include <maya/MString.h>
#include <maya/MTypeId.h>

#include <cuda.h>

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

  static CUcontext PiesCudaContext;

  static MString pluginDir;
  static MTypeId id;
  static MObject time;
  static MObject prevTime;

  static MObject stepSize;
  static MObject substeps;
  static MObject iterations;

  // Soft body compound attribute
  static MObject strainStiffness;
  static MObject strainRange;
  static MObject volStiffness;
  static MObject volRange;
  static MObject softBodyDensity;
  static MObject softBodyVelocity;
  static MObject softBodyMesh;
  static MObject softBodyArray;
  
  static MObject collisionIterations;
  static MObject collisionDistance;
  static MObject collisionThickness;
  static MObject gridSpacing;
  static MObject gravity;
  static MObject damping;
  static MObject friction;
  static MObject floorHeight;
  static MObject threadCount;
  
  static MObject simulationEnabled;
  static MObject simulationStartTime;

  static MObject fixedRegionsArray;
  static MObject linkedRegionsArray;

  static MObject outputPositions;
  static MObject outputMesh;

private:
  std::unique_ptr<Pies::Solver> _pSolver;
  mutable bool _resetSimulation = true;
  float _simulationTime = 0.0f;
};
