#pragma once

#include <Pies/Solver.h>

#include <maya/MStatus.h>
#include <maya/MPxNode.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MObject.h>
#include <maya/MTypeId.h>
#include <maya/MString.h>

#include <string>
#include <memory>

// TODO: Rename to PiesSolvernode
class SolverNode : public MPxNode {
public:
	SolverNode();
	virtual MStatus compute(const MPlug& plug, MDataBlock& data);
	static void* creator();
	static MStatus initialize();

	static MString pluginDir;
	static MTypeId id;
	static MObject time;
	static MObject stepSize;
	static MObject iterations;
	static MObject ouptutPositions;

private:
  static std::unique_ptr<Pies::Solver> _pSolver;
};
