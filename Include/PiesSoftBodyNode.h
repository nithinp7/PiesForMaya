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

class PiesSoftBodyNode : public MPxNode {
public:
  PiesSoftBodyNode() = default;
  MStatus compute(const MPlug& plug, MDataBlock& data) override;
  static void* creator();
  static MStatus initialize();

  static MTypeId id;

  static MObject inMesh;

  static MObject strainStiffness;
  static MObject minStrain;
  static MObject maxStrain;

  static MObject volStiffness;
  static MObject volMultiplier;

  static MObject outMesh;

  static MObject outStrainStiffness;
  static MObject outMinStrain;
  static MObject outMaxStrain;

  static MObject outVolStiffness;
  static MObject outVolMultiplier;

  static MObject outCompound;
};