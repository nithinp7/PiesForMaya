#include "PiesClothNode.h"

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
MTypeId PiesClothNode::id(0x00000004);

/*static*/
MObject PiesClothNode::inMesh;

/*static*/
MObject PiesClothNode::stretchStiffness;

/*static*/
MObject PiesClothNode::bendStiffness;

/*static*/
MObject PiesClothNode::outMesh;

/*static*/
MObject PiesClothNode::outStretchStiffness;

/*static*/
MObject PiesClothNode::outBendStiffness;

/*static*/
MObject PiesClothNode::outCompound;

/*static*/
void* PiesClothNode::creator() {
  return new PiesClothNode();
}

MStatus PiesClothNode::compute(const MPlug& plug, MDataBlock& data) {
  MStatus status = MStatus::kSuccess;

  if (plug == outCompound) {
    float stretchStiffnessValue =
        data.inputValue(stretchStiffness, &status).asFloat();

    float bendStiffnessValue = data.inputValue(bendStiffness, &status).asFloat();

    MDataHandle inMeshHandle = data.inputValue(inMesh, &status);
    McheckErr(status, "ERROR getting inMesh handle.");

    MDataHandle compoundHandle = data.outputValue(outCompound, &status);
    McheckErr(status, "ERROR getting outCompound handle.");
    
    MDataHandle outMeshHandle = compoundHandle.child(outMesh);

    outMeshHandle.copy(inMeshHandle);
    MFnMesh meshOut(outMeshHandle.asMesh());

    compoundHandle.child(outStretchStiffness).setFloat(stretchStiffnessValue);
    compoundHandle.child(outBendStiffness).setFloat(bendStiffnessValue);

    data.setClean(outCompound);
  } 
  else {
    return MStatus::kUnknownParameter;
  }

  return status;
}

/*static*/
MStatus PiesClothNode::initialize() {
  MStatus status = MStatus::kSuccess;

  // Input attributes

  MFnNumericAttribute stretchAttr;
  PiesClothNode::stretchStiffness = stretchAttr.create(
      "stretchStiffness",
      "stretchStiff",
      MFnNumericData::kFloat,
      1000.0,
      &status);
  McheckErr(status, "ERROR creating attribute PiesClothNode::stretchStiffness.");

  stretchAttr.setWritable(true);
  stretchAttr.setKeyable(false);
  status = addAttribute(PiesClothNode::stretchStiffness);
  McheckErr(status, "ERROR adding attribute PiesClothNode::stretchStiffness.");

  MFnNumericAttribute bendStiffnessAttr;
  PiesClothNode::bendStiffness = bendStiffnessAttr.create(
      "bendStiffness",
      "bendStiff",
      MFnNumericData::kFloat,
      1000.0,
      &status);
  McheckErr(status, "ERROR creating attribute PiesClothNode::bendStiffness.");

  bendStiffnessAttr.setWritable(true);
  bendStiffnessAttr.setKeyable(false);
  status = addAttribute(PiesClothNode::bendStiffness);
  McheckErr(status, "ERROR adding attribute PiesClothNode::bendStiffness.");

  MFnTypedAttribute inMeshAttr;
  PiesClothNode::inMesh = inMeshAttr.create(
      "inMesh",
      "inMesh",
      MFnData::kMesh,
      MObject::kNullObj,
      &status);
  McheckErr(status, "ERROR creating attribute PiesClothNode::inMesh.");

  inMeshAttr.setWritable(true);
  inMeshAttr.setKeyable(false);
  status = addAttribute(PiesClothNode::inMesh);
  McheckErr(status, "ERROR adding attribute PiesClothNode::inMesh.");

  // Elements of the compound output attribute

  MFnNumericAttribute outStretchAttr;
  PiesClothNode::outStretchStiffness = outStretchAttr.create(
      "outStretchStiffness",
      "ss",
      MFnNumericData::kFloat,
      1000.0,
      &status);
  outStretchAttr.setWritable(false);
  outStretchAttr.setStorable(true);
  McheckErr(status, "ERROR creating attribute PiesClothNode::outStretchStiffness.");

  MFnNumericAttribute outBendStiffnessAttr;
  PiesClothNode::outBendStiffness = outBendStiffnessAttr.create(
      "outBendStiffness",
      "bs",
      MFnNumericData::kFloat,
      1000.0,
      &status);
  outBendStiffnessAttr.setWritable(false);
  outBendStiffnessAttr.setStorable(true);
  McheckErr(status, "ERROR creating attribute PiesClothNode::outBendStiffness.");

  MFnTypedAttribute outMeshAttr;
  PiesClothNode::outMesh = outMeshAttr.create(
      "outMesh",
      "outMesh",
      MFnData::kMesh,
      MObject::kNullObj,
      &status);
  outMeshAttr.setWritable(false);
  outMeshAttr.setStorable(true);
  outMeshAttr.setCached(true);
  McheckErr(status, "ERROR creating attribute PiesClothNode::outMesh.");

  MFnCompoundAttribute compoundAttr;
  PiesClothNode::outCompound = compoundAttr.create(
      "output",
      "out",
      &status);
  McheckErr(status, "ERROR creating attribute PiesClothNode::outCompound.");
  status = compoundAttr.addChild(PiesClothNode::outStretchStiffness);
  McheckErr(status, "ERROR adding compound children.");
  status = compoundAttr.addChild(PiesClothNode::outBendStiffness);
  McheckErr(status, "ERROR adding compound children.");
  status = compoundAttr.addChild(PiesClothNode::outMesh);
  McheckErr(status, "ERROR adding compound children.");

  compoundAttr.setWritable(false);
  compoundAttr.setStorable(true);
  compoundAttr.setCached(true);
  status = addAttribute(PiesClothNode::outCompound);
  McheckErr(status, "ERROR adding compound attribute.");
  
  // Setup input-output dependencies
  status = attributeAffects(inMesh, outCompound);
  McheckErr(status, "ERROR attributeAffects(inMesh, outCompound).");

  status = attributeAffects(stretchStiffness, outCompound);
  McheckErr(status, "ERROR attributeAffects(stretchStiffness, outCompound).");

  status = attributeAffects(bendStiffness, outCompound);
  McheckErr(status, "ERROR attributeAffects(bendStiffness, outCompound).");

  return status;
}
