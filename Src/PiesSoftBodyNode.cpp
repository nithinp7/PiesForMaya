#include "PiesSoftBodyNode.h"

#include <maya/MFloatArray.h>
#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MGlobal.h>
#include <maya/MIOStream.h>
#include <maya/MStringArray.h>
#include <maya/MFnCompoundAttribute.h>

#define McheckErr(stat, msg)                                                   \
  if (!stat) {                                                                 \
    MGlobal::displayInfo(msg);                                                 \
    stat.perror(msg);                                                          \
    return stat;                                                               \
  }

/*static*/
MTypeId PiesSoftBodyNode::id(0x00000003);

/*static*/
MObject PiesSoftBodyNode::inMesh;

/*static*/
MObject PiesSoftBodyNode::strainStiffness;

/*static*/
MObject PiesSoftBodyNode::minStrain;

/*static*/
MObject PiesSoftBodyNode::maxStrain;

/*static*/
MObject PiesSoftBodyNode::volStiffness;

/*static*/
MObject PiesSoftBodyNode::volMultiplier;

/*static*/
MObject PiesSoftBodyNode::outMesh;

/*static*/
MObject PiesSoftBodyNode::outStrainStiffness;

/*static*/
MObject PiesSoftBodyNode::outMinStrain;

/*static*/
MObject PiesSoftBodyNode::outMaxStrain;

/*static*/
MObject PiesSoftBodyNode::outVolStiffness;

/*static*/
MObject PiesSoftBodyNode::outVolMultiplier;

/*static*/
MObject PiesSoftBodyNode::outCompound;

/*static*/
void* PiesSoftBodyNode::creator() {
  return new PiesSoftBodyNode();
}

MStatus PiesSoftBodyNode::compute(const MPlug& plug, MDataBlock& data) {
  MStatus status = MStatus::kSuccess;

  if (plug == outCompound) {
    float strainStiffnessValue =
        data.inputValue(strainStiffness, &status).asFloat();
    float minStrainValue = data.inputValue(minStrain, &status).asFloat();
    float maxStrainValue = data.inputValue(maxStrain, &status).asFloat();

    float volStiffValue = data.inputValue(volStiffness, &status).asFloat();
    float volMultValue = data.inputValue(volMultiplier, &status).asFloat();

    MDataHandle inMeshHandle = data.inputValue(inMesh, &status);
    McheckErr(status, "ERROR getting inMesh handle.");

    MDataHandle compoundHandle = data.outputValue(outCompound, &status);
    McheckErr(status, "ERROR getting outCompound handle.");
    
    MDataHandle outMeshHandle = compoundHandle.child(outMesh);

    outMeshHandle.copy(inMeshHandle);
    MFnMesh meshOut(outMeshHandle.asMesh());

    compoundHandle.child(outStrainStiffness).setFloat(strainStiffnessValue);
    compoundHandle.child(outMinStrain).setFloat(minStrainValue);
    compoundHandle.child(outMaxStrain).setFloat(maxStrainValue);
    compoundHandle.child(outVolStiffness).setFloat(volStiffValue);
    compoundHandle.child(outVolMultiplier).setFloat(volMultValue);

    data.setClean(outCompound);
  } else {
    return MStatus::kUnknownParameter;
  }

  return status;
}

/*static*/
MStatus PiesSoftBodyNode::initialize() {
  MStatus status = MStatus::kSuccess;

  // Input attributes

  MFnNumericAttribute strainAttr;
  PiesSoftBodyNode::strainStiffness = strainAttr.create(
      "strainStiffness",
      "strainStiff",
      MFnNumericData::kFloat,
      1000.0,
      &status);
  McheckErr(status, "ERROR creating attribute PiesSoftBodyNode::strainStiffness.");

  strainAttr.setWritable(true);
  strainAttr.setKeyable(false);
  status = addAttribute(PiesSoftBodyNode::strainStiffness);
  McheckErr(status, "ERROR adding attribute PiesSoftBodyNode::strainStiffness.");

  MFnNumericAttribute minStrainAttr;
  PiesSoftBodyNode::minStrain =
      minStrainAttr
          .create("minStrain", "minS", MFnNumericData::kFloat, 0.8, &status);
  McheckErr(status, "ERROR creating attribute PiesSoftBodyNode::minStrain.");

  minStrainAttr.setWritable(true);
  minStrainAttr.setKeyable(false);
  status = addAttribute(PiesSoftBodyNode::minStrain);
  McheckErr(status, "ERROR adding attribute PiesSoftBodyNode::minStrain.");

  MFnNumericAttribute maxStrainAttr;
  PiesSoftBodyNode::maxStrain =
      maxStrainAttr
          .create("maxStrain", "maxS", MFnNumericData::kFloat, 1.0, &status);
  McheckErr(status, "ERROR creating attribute PiesSoftBodyNode::maxStrain.");

  maxStrainAttr.setWritable(true);
  maxStrainAttr.setKeyable(false);
  status = addAttribute(PiesSoftBodyNode::maxStrain);
  McheckErr(status, "ERROR adding attribute PiesSoftBodyNode::maxStrain.");

  MFnNumericAttribute volStiffnessAttr;
  PiesSoftBodyNode::volStiffness = volStiffnessAttr.create(
      "volumeStiffness",
      "volStiff",
      MFnNumericData::kFloat,
      1000.0,
      &status);
  McheckErr(status, "ERROR creating attribute PiesSoftBodyNode::volStiffness.");

  volStiffnessAttr.setWritable(true);
  volStiffnessAttr.setKeyable(false);
  status = addAttribute(PiesSoftBodyNode::volStiffness);
  McheckErr(status, "ERROR adding attribute PiesSoftBodyNode::volStiffness.");

  MFnNumericAttribute volMultiplierAttr;
  PiesSoftBodyNode::volMultiplier = volMultiplierAttr.create(
      "volumeMultiplier",
      "volMult",
      MFnNumericData::kFloat,
      1.0,
      &status);
  McheckErr(status, "ERROR creating attribute PiesSoftBodyNode::volMultiplier.");

  volMultiplierAttr.setWritable(true);
  volMultiplierAttr.setKeyable(false);
  status = addAttribute(PiesSoftBodyNode::volMultiplier);
  McheckErr(status, "ERROR adding attribute PiesSoftBodyNode::volMultiplier.");

  MFnTypedAttribute inMeshAttr;
  PiesSoftBodyNode::inMesh = inMeshAttr.create(
      "inMesh",
      "inMesh",
      MFnData::kMesh,
      MObject::kNullObj,
      &status);
  McheckErr(status, "ERROR creating attribute PiesSoftBodyNode::inMesh.");

  inMeshAttr.setWritable(true);
  inMeshAttr.setKeyable(false);
  status = addAttribute(PiesSoftBodyNode::inMesh);
  McheckErr(status, "ERROR adding attribute PiesSoftBodyNode::inMesh.");

  // Elements of the compound output attribute

  MFnNumericAttribute outStrainAttr;
  PiesSoftBodyNode::outStrainStiffness = outStrainAttr.create(
      "outStrainStiffness",
      "ss",
      MFnNumericData::kFloat,
      1000.0,
      &status);
  outStrainAttr.setWritable(false);
  outStrainAttr.setStorable(true);
  McheckErr(status, "ERROR creating attribute PiesSoftBodyNode::outStrainStiffness.");

  MFnNumericAttribute outMinStrainAttr;
  PiesSoftBodyNode::outMinStrain =
      outMinStrainAttr
          .create("outMinStrain", "oMinS", MFnNumericData::kFloat, 0.8, &status);
  outMinStrainAttr.setWritable(false);
  outMinStrainAttr.setStorable(true);
  McheckErr(status, "ERROR creating attribute PiesSoftBodyNode::outMinStrain.");

  MFnNumericAttribute outMaxStrainAttr;
  PiesSoftBodyNode::outMaxStrain =
      outMaxStrainAttr
          .create("outMaxStrain", "oMaxS", MFnNumericData::kFloat, 1.0, &status);
  outMaxStrainAttr.setWritable(false);
  outMaxStrainAttr.setStorable(true);
  McheckErr(status, "ERROR creating attribute PiesSoftBodyNode::outMaxStrain.");

  MFnNumericAttribute outVolStiffnessAttr;
  PiesSoftBodyNode::outVolStiffness = outVolStiffnessAttr.create(
      "outVolumeStiffness",
      "vs",
      MFnNumericData::kFloat,
      1000.0,
      &status);
  outVolStiffnessAttr.setWritable(false);
  outVolStiffnessAttr.setStorable(true);
  McheckErr(status, "ERROR creating attribute PiesSoftBodyNode::outVolStiffness.");

  MFnNumericAttribute outVolMultiplierAttr;
  PiesSoftBodyNode::outVolMultiplier = outVolMultiplierAttr.create(
      "outVolumeMultiplier",
      "vm",
      MFnNumericData::kFloat,
      1.0,
      &status);
  outVolMultiplierAttr.setWritable(false);
  outVolMultiplierAttr.setStorable(true);
  McheckErr(status, "ERROR creating attribute PiesSoftBodyNode::outVolMultiplier.");

  MFnTypedAttribute outMeshAttr;
  PiesSoftBodyNode::outMesh = outMeshAttr.create(
      "outMesh",
      "outMesh",
      MFnData::kMesh,
      MObject::kNullObj,
      &status);
  outMeshAttr.setWritable(false);
  outMeshAttr.setStorable(true);
  outMeshAttr.setCached(true);
  McheckErr(status, "ERROR creating attribute PiesSoftBodyNode::outMesh.");

  MFnCompoundAttribute compoundAttr;
  PiesSoftBodyNode::outCompound = compoundAttr.create(
      "output",
      "out",
      &status);
  McheckErr(status, "ERROR creating attribute PiesSoftBodyNode::outCompound.");
  status = compoundAttr.addChild(PiesSoftBodyNode::outStrainStiffness);
  McheckErr(status, "ERROR adding compound children.");
  status = compoundAttr.addChild(PiesSoftBodyNode::outMinStrain);
  McheckErr(status, "ERROR adding compound children.");
  status = compoundAttr.addChild(PiesSoftBodyNode::outMaxStrain);
  McheckErr(status, "ERROR adding compound children.");
  status = compoundAttr.addChild(PiesSoftBodyNode::outVolStiffness);
  McheckErr(status, "ERROR adding compound children.");
  status = compoundAttr.addChild(PiesSoftBodyNode::outVolMultiplier);
  McheckErr(status, "ERROR adding compound children.");
  status = compoundAttr.addChild(PiesSoftBodyNode::outMesh);
  McheckErr(status, "ERROR adding compound children.");

  compoundAttr.setWritable(false);
  compoundAttr.setStorable(true);
  compoundAttr.setCached(true);
  status = addAttribute(PiesSoftBodyNode::outCompound);
  McheckErr(status, "ERROR adding compound attribute.");
  
  // Setup input-output dependencies
  status = attributeAffects(inMesh, outCompound);
  McheckErr(status, "ERROR attributeAffects(inMesh, outCompound).");

  status = attributeAffects(strainStiffness, outCompound);
  McheckErr(status, "ERROR attributeAffects(strainStiffness, outCompound).");

  status = attributeAffects(minStrain, outCompound);
  McheckErr(status, "ERROR attributeAffects(minStrain, outCompound).");

  status = attributeAffects(maxStrain, outCompound);
  McheckErr(status, "ERROR attributeAffects(maxStrain, outCompound).");

  status = attributeAffects(volStiffness, outCompound);
  McheckErr(status, "ERROR attributeAffects(volStiffness, outCompound).");

  status = attributeAffects(volMultiplier, outCompound);
  McheckErr(status, "ERROR attributeAffects(volMultiplier, outCompound).");

  return status;
}