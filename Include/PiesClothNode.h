#pragma once

#include <maya/MPxNode.h>
#include <maya/MFnTypedAttribute.h>
#include <maya/MFnUnitAttribute.h>
#include <maya/MObject.h>
#include <maya/MPxNode.h>
#include <maya/MStatus.h>
#include <maya/MString.h>
#include <maya/MTypeId.h>

#include <memory>
#include <string>

class PiesClothNode : public MPxNode {
public:
  PiesClothNode() = default;
  MStatus compute(const MPlug& plug, MDataBlock& data) override;
  static void* creator();
  static MStatus initialize();

  static MTypeId id;

  static MObject inMesh;

  static MObject stretchStiffness; // 

  static MObject bendStiffness;
  // removed

  static MObject outMesh;

  static MObject outStretchStiffness;

  static MObject outBendStiffness;

  static MObject outCompound;
};
