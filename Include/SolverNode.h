#pragma once

#include <maya/MStatus.h>
#include <maya/MPxNode.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MObject.h>
#include <maya/MTypeId.h>
#include <maya/MString.h>

#include <string>

class SolverNode : public MPxNode {
public:
	SolverNode() = default;
	virtual MStatus compute(const MPlug& plug, MDataBlock& data);
	static void* creator();
	static MStatus initialize();

	static MString pluginDir;
	static MTypeId id;
	static MObject time;
	static MObject stepSize;
	static MObject angle;
	static MObject iterations;
	static MObject grammar;
	static MObject outputMesh;
};
