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
  static MObject strainRange;

  static MObject volStiffness;
  static MObject volRange;

  static MObject density;
  static MObject velocity;
  
  static MObject outMesh;

  static MObject outStrainStiffness;
  static MObject outStrainRange;

  static MObject outVolStiffness;
  static MObject outVolRange;

  static MObject outDensity;
  static MObject outVelocity;

  static MObject outCompound;
};