#include "PiesSoftBodyNode.h"

#include <maya/MFloatArray.h>
#include <maya/MFnCompoundAttribute.h>
#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MGlobal.h>
#include <maya/MIOStream.h>
#include <maya/MStringArray.h>

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
MObject PiesSoftBodyNode::strainRange;

/*static*/
MObject PiesSoftBodyNode::volStiffness;

/*static*/
MObject PiesSoftBodyNode::volRange;

/*static*/
MObject PiesSoftBodyNode::density;

/*static*/
MObject PiesSoftBodyNode::velocity;

/*static*/
MObject PiesSoftBodyNode::outMesh;

/*static*/
MObject PiesSoftBodyNode::outStrainStiffness;

/*static*/
MObject PiesSoftBodyNode::outStrainRange;

/*static*/
MObject PiesSoftBodyNode::outVolStiffness;

/*static*/
MObject PiesSoftBodyNode::outVolRange;

/*static*/
MObject PiesSoftBodyNode::outDensity;

/*static*/
MObject PiesSoftBodyNode::outVelocity;

/*static*/
MObject PiesSoftBodyNode::outCompound;

/*static*/
void* PiesSoftBodyNode::creator() { return new PiesSoftBodyNode(); }

MStatus PiesSoftBodyNode::compute(const MPlug& plug, MDataBlock& data) {
  MStatus status = MStatus::kSuccess;

  if (plug == outCompound) {
    float strainStiffnessValue =
        data.inputValue(strainStiffness, &status).asFloat();
    const float2& strainRangeValue =
        data.inputValue(strainRange, &status).asFloat2();

    float volStiffValue = data.inputValue(volStiffness, &status).asFloat();
    const float2& volRangeValue = data.inputValue(volRange, &status).asFloat2();

    float densityValue = data.inputValue(density, &status).asFloat();

    const float3& velocityValue = data.inputValue(velocity, &status).asFloat3();

    MDataHandle inMeshHandle = data.inputValue(inMesh, &status);
    McheckErr(status, "ERROR getting inMesh handle.");

    MDataHandle compoundHandle = data.outputValue(outCompound, &status);
    McheckErr(status, "ERROR getting outCompound handle.");

    MDataHandle outMeshHandle = compoundHandle.child(outMesh);

    outMeshHandle.copy(inMeshHandle);
    MFnMesh meshOut(outMeshHandle.asMesh());

    compoundHandle.child(outStrainStiffness).setFloat(strainStiffnessValue);
    compoundHandle.child(outStrainRange)
        .set2Float(strainRangeValue[0], strainRangeValue[1]);
    compoundHandle.child(outVolStiffness).setFloat(volStiffValue);
    compoundHandle.child(outVolRange)
        .set2Float(volRangeValue[0], volRangeValue[1]);
    compoundHandle.child(outDensity).setFloat(densityValue);
    compoundHandle.child(outVelocity)
        .set3Float(velocityValue[0], velocityValue[1], velocityValue[2]);

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
  McheckErr(
      status,
      "ERROR creating attribute PiesSoftBodyNode::strainStiffness.");

  strainAttr.setWritable(true);
  strainAttr.setKeyable(false);
  status = addAttribute(PiesSoftBodyNode::strainStiffness);
  McheckErr(
      status,
      "ERROR adding attribute PiesSoftBodyNode::strainStiffness.");

  MFnNumericAttribute strainRangeAttr;
  PiesSoftBodyNode::strainRange =
      strainRangeAttr.create("strainRange", "strainR", MFnNumericData::k2Float);

  strainRangeAttr.setDefault(1.0f, 1.0f);
  strainRangeAttr.setWritable(true);
  strainRangeAttr.setKeyable(false);
  status = addAttribute(PiesSoftBodyNode::strainRange);
  McheckErr(status, "ERROR adding attribute PiesSoftBodyNode::strainRange.");

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

  MFnNumericAttribute volRangeAttr;
  PiesSoftBodyNode::volRange =
      volRangeAttr.create("volRange", "volRange", MFnNumericData::k2Float);
  volRangeAttr.setDefault(0.5f, 1.5f);
  volRangeAttr.setWritable(true);
  volRangeAttr.setKeyable(false);
  status = addAttribute(PiesSoftBodyNode::volRange);
  McheckErr(status, "ERROR adding attribute PiesSoftBodyNode::volRange.");

  MFnNumericAttribute densityAttr;
  PiesSoftBodyNode::density =
      densityAttr.create("density", "p", MFnNumericData::kFloat, 1.0, &status);
  McheckErr(status, "ERROR creating attribute PiesSoftBodyNode::density.");

  densityAttr.setWritable(true);
  densityAttr.setKeyable(false);
  status = addAttribute(PiesSoftBodyNode::density);
  McheckErr(status, "ERROR adding attribute PiesSoftBodyNode::density.");

  MFnNumericAttribute velocityAttr;
  PiesSoftBodyNode::velocity =
      velocityAttr.create("velocity", "vel", MFnNumericData::k3Float);
  velocityAttr.setDefault(0.0f, 0.0f, 0.0f);
  velocityAttr.setWritable(true);
  velocityAttr.setKeyable(false);
  status = addAttribute(PiesSoftBodyNode::velocity);
  McheckErr(status, "ERROR adding attribute PiesSoftBodyNode::velocity.");

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
  McheckErr(
      status,
      "ERROR creating attribute PiesSoftBodyNode::outStrainStiffness.");

  MFnNumericAttribute outStrainRangeAttr;
  PiesSoftBodyNode::outStrainRange = outStrainRangeAttr.create(
      "outStrainRange",
      "oStrainR",
      MFnNumericData::k2Float);
  outStrainRangeAttr.setDefault(0.8f, 1.0f);
  outStrainRangeAttr.setWritable(false);
  outStrainRangeAttr.setStorable(true);
  McheckErr(status, "ERROR creating attribute PiesSoftBodyNode::outMinStrain.");

  MFnNumericAttribute outVolStiffnessAttr;
  PiesSoftBodyNode::outVolStiffness = outVolStiffnessAttr.create(
      "outVolumeStiffness",
      "vs",
      MFnNumericData::kFloat,
      1000.0,
      &status);
  outVolStiffnessAttr.setWritable(false);
  outVolStiffnessAttr.setStorable(true);
  McheckErr(
      status,
      "ERROR creating attribute PiesSoftBodyNode::outVolStiffness.");

  MFnNumericAttribute outVolRangeAttr;
  PiesSoftBodyNode::outVolRange = outVolRangeAttr.create(
      "outVolRange",
      "oVolRange",
      MFnNumericData::k2Float);
  outVolRangeAttr.setDefault(1.0f, 1.0f);
  outVolRangeAttr.setWritable(false);
  outVolRangeAttr.setStorable(true);
  McheckErr(status, "ERROR creating attribute PiesSoftBodyNode::outVolRange.");

  MFnNumericAttribute outDensityAttr;
  PiesSoftBodyNode::outDensity =
      outDensityAttr
          .create("outDensity", "outP", MFnNumericData::kFloat, 1.0, &status);
  outDensityAttr.setWritable(false);
  outDensityAttr.setStorable(true);
  McheckErr(status, "ERROR creating attribute PiesSoftBodyNode::outDensity.");

  MFnNumericAttribute outVelocityAttr;
  PiesSoftBodyNode::outVelocity =
      outVelocityAttr.create("outVelocity", "outVel", MFnNumericData::k3Float);
  outVelocityAttr.setDefault(0.0f, 0.0f, 0.0f);
  outVelocityAttr.setWritable(false);
  outVelocityAttr.setStorable(true);
  McheckErr(status, "ERROR creating attribute PiesSoftBodyNode::outVelocity.");

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
  PiesSoftBodyNode::outCompound = compoundAttr.create("output", "out", &status);
  McheckErr(status, "ERROR creating attribute PiesSoftBodyNode::outCompound.");
  status = compoundAttr.addChild(PiesSoftBodyNode::outStrainStiffness);
  McheckErr(status, "ERROR adding compound children.");
  status = compoundAttr.addChild(PiesSoftBodyNode::outStrainRange);
  McheckErr(status, "ERROR adding compound children.");
  status = compoundAttr.addChild(PiesSoftBodyNode::outVolStiffness);
  McheckErr(status, "ERROR adding compound children.");
  status = compoundAttr.addChild(PiesSoftBodyNode::outVolRange);
  McheckErr(status, "ERROR adding compound children.");
  status = compoundAttr.addChild(PiesSoftBodyNode::outDensity);
  McheckErr(status, "ERROR adding compound children.");
  status = compoundAttr.addChild(PiesSoftBodyNode::outVelocity);
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

  status = attributeAffects(strainRange, outCompound);
  McheckErr(status, "ERROR attributeAffects(strainRange, outCompound).");

  status = attributeAffects(volStiffness, outCompound);
  McheckErr(status, "ERROR attributeAffects(volStiffness, outCompound).");

  status = attributeAffects(volRange, outCompound);
  McheckErr(status, "ERROR attributeAffects(volRange, outCompound).");

  status = attributeAffects(density, outCompound);
  McheckErr(status, "ERROR attributeAffects(density, outCompound).");

  status = attributeAffects(velocity, outCompound);
  McheckErr(status, "ERROR attributeAffects(velocity, outCompound).");

  return status;
}