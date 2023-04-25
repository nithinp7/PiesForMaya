#include "PiesVolumeNode.h"

#include <maya/MFloatArray.h>
#include <maya/MFnMesh.h>
#include <maya/MFnMeshData.h>
#include <maya/MFnNumericAttribute.h>
#include <maya/MGlobal.h>
#include <maya/MIOStream.h>
#include <maya/MStringArray.h>
#include <maya/adskDataAssociations.h>
#include <maya/adskDataStream.h>
#include <maya/adskDataStructure.h>

#define McheckErr(stat, msg)                                                   \
  if (!stat) {                                                                 \
    MGlobal::displayInfo(msg);                                                 \
    stat.perror(msg);                                                          \
    return stat;                                                               \
  }

/*static*/
MTypeId PiesVolumeNode::id(0x00000002);

/*static*/
MObject PiesVolumeNode::inMesh;

/*static*/
MObject PiesVolumeNode::strainStiffness;

/*static*/
MObject PiesVolumeNode::minStrain;

/*static*/
MObject PiesVolumeNode::maxStrain;

/*static*/
MObject PiesVolumeNode::preserveVolume;

/*static*/
MObject PiesVolumeNode::volStiffness;

/*static*/
MObject PiesVolumeNode::volMultiplier;

/*static*/
MObject PiesVolumeNode::outMesh;

// MStatus PiesVolumeNode::compute(const MPlug& plug, MDataBlock& data) {
  
MStatus	deform(MDataBlock& block,
          MItGeometry& iter,
          const MMatrix& mat,
          unsigned int	multiIndex) {
            
  MStatus status = MStatus::kSuccess;

  if (plug == outMesh) {
    float strainStiffnessValue =
        data.inputValue(strainStiffness, &status).asFloat();
    float minStrainValue = data.inputValue(minStrain, &status).asFloat();
    float maxStrainValue = data.inputValue(maxStrain, &status).asFloat();

    float presVolValue = data.inputValue(preserveVolume, &status).asFloat();
    float volStiffValue = data.inputValue(volStiffness, &status).asFloat();
    float volMultValue = data.inputValue(volMultiplier, &status).asFloat();

    MDataHandle inMeshHandle = data.inputValue(inMesh, &status);
    McheckErr(status, "ERROR getting inMesh handle.");

    MDataHandle outMeshHandle = data.inputValue(outMesh, &status);
    McheckErr(status, "ERROR getting outMesh handle.");

    outMeshHandle.copy(inMeshHandle);

    MFnMesh meshOut(outMeshHandle.asMesh());

    MStringArray longNames;
    longNames.append("strainStiffness");
    longNames.append("minStrain");
    longNames.append("maxStrain");

    longNames.append("preserveVolume");
    longNames.append("volStiffness");
    longNames.append("volMultiplier");

    MStringArray shortNames;
    shortNames.append("strainStiff");
    shortNames.append("minStrain");
    shortNames.append("maxStrain");

    shortNames.append("preserveVolume");
    shortNames.append("volStiffness");
    shortNames.append("volMultiplier");

    MStringArray formatNames;
    formatNames.append("float");
    formatNames.append("float");
    formatNames.append("float");

    formatNames.append("float");
    formatNames.append("float");
    formatNames.append("float");

    meshOut.createBlindDataType(1, longNames, shortNames, formatNames);

    // Define struct encapsulating soft body data
    adsk::Data::Structure* PiesSoftBodyStruct = adsk::Data::Structure::create();
    PiesSoftBodyStruct->setName("PiesSoftBody");
    PiesSoftBodyStruct->addMember(
        adsk::Data::Member::eDataType::kFloat,
        1,
        "strainStiffness");
    PiesSoftBodyStruct->addMember(
        adsk::Data::Member::eDataType::kFloat,
        1,
        "minStrain");
    PiesSoftBodyStruct->addMember(
        adsk::Data::Member::eDataType::kFloat,
        1,
        "maxStrain");
    PiesSoftBodyStruct->addMember(
        adsk::Data::Member::eDataType::kBoolean,
        1,
        "preserveVolume");
    PiesSoftBodyStruct->addMember(
        adsk::Data::Member::eDataType::kFloat,
        1,
        "volStiffness");
    PiesSoftBodyStruct->addMember(
        adsk::Data::Member::eDataType::kFloat,
        1,
        "volMultiplier");


    adsk::Data::Handle softBodyInfo(*PiesSoftBodyStruct);
    softBodyInfo.setPositionByMemberIndex(0);
    *softBodyInfo.asFloat() = strainStiffnessValue;
    softBodyInfo.setPositionByMemberIndex(1);
    *softBodyInfo.asFloat() = minStrainValue;
    softBodyInfo.setPositionByMemberIndex(2);
    *softBodyInfo.asFloat() = maxStrainValue;
    softBodyInfo.setPositionByMemberIndex(3);
    *softBodyInfo.asBoolean() = presVolValue;
    softBodyInfo.setPositionByMemberIndex(4);
    *softBodyInfo.asFloat() = volStiffValue;
    softBodyInfo.setPositionByMemberIndex(5);
    *softBodyInfo.asFloat() = volMultValue;


    // Enumerate unique stream names?
    adsk::Data::Stream piesData(*PiesSoftBodyStruct, "piesData");

    adsk::Data::Channel channel("PiesChannel");
    channel.setDataStream(piesData);
    
    adsk::Data::Associations* assoc = adsk::Data::Associations::create();
    assoc->setChannel(channel);

    // meshOut.
    //  adsk::Data::Associations metadata;
    //  adsk::Data::Channel channel("strainStiffness");
    //  adsk::Data::Structure piesStruct();

    // meshOut.setMetadata()

    // MFloatArray strainStiffnessArr, minStrainArr;
    // strainStiffnessArr.setLength(meshOut.numVertices());

    McheckErr(status, "ERROR copying inMesh to outMesh.");
  }

  return status;
}

/*static*/
MStatus PiesVolumeNode::initialize() {
  MStatus status = MStatus::kSuccess;

  MFnNumericAttribute strainAttr;
  PiesVolumeNode::strainStiffness = strainAttr.create(
      "strainStiffness",
      "strainStiff",
      MFnNumericData::kFloat,
      1000.0,
      &status);
  McheckErr(status, "ERROR creating attribute PiesVolumeNode::stepSize.");

  strainAttr.setWritable(true);
  strainAttr.setKeyable(false);
  status = addAttribute(PiesVolumeNode::strainStiffness);
  McheckErr(status, "ERROR adding attribute PiesVolumeNode::strainStiffness.");

  MFnNumericAttribute minStrainAttr;
  PiesVolumeNode::minStrain =
      minStrainAttr
          .create("minStrain", "minS", MFnNumericData::kFloat, 0.8, &status);
  McheckErr(status, "ERROR creating attribute PiesVolumeNode::minStrain.");

  minStrainAttr.setWritable(true);
  minStrainAttr.setKeyable(false);
  status = addAttribute(PiesVolumeNode::minStrain);
  McheckErr(status, "ERROR adding attribute PiesVolumeNode::minStrain.");

  MFnNumericAttribute maxStrainAttr;
  PiesVolumeNode::maxStrain =
      maxStrainAttr
          .create("maxStrain", "maxS", MFnNumericData::kFloat, 1.0, &status);
  McheckErr(status, "ERROR creating attribute PiesVolumeNode::maxStrain.");

  maxStrainAttr.setWritable(true);
  maxStrainAttr.setKeyable(false);
  status = addAttribute(PiesVolumeNode::maxStrain);
  McheckErr(status, "ERROR adding attribute PiesVolumeNode::maxStrain.");

  MFnNumericAttribute preserveVolumeAttr;
  PiesVolumeNode::preserveVolume = preserveVolumeAttr.create(
      "preserveVolume",
      "presVol",
      MFnNumericData::kBoolean,
      true,
      &status);
  McheckErr(status, "ERROR creating attribute PiesVolumeNode::preserveVolume.");

  preserveVolumeAttr.setWritable(true);
  preserveVolumeAttr.setKeyable(false);
  status = addAttribute(PiesVolumeNode::preserveVolume);
  McheckErr(status, "ERROR adding attribute PiesVolumeNode::preserveVolume.");

  MFnNumericAttribute volStiffnessAttr;
  PiesVolumeNode::volStiffness = volStiffnessAttr.create(
      "volumeStiffness",
      "volStiff",
      MFnNumericData::kFloat,
      1000.0,
      &status);
  McheckErr(status, "ERROR creating attribute PiesVolumeNode::volStiffness.");

  volStiffnessAttr.setWritable(true);
  volStiffnessAttr.setKeyable(false);
  status = addAttribute(PiesVolumeNode::volStiffness);
  McheckErr(status, "ERROR adding attribute PiesVolumeNode::volStiffness.");

  MFnNumericAttribute volMultiplierAttr;
  PiesVolumeNode::volMultiplier = volMultiplierAttr.create(
      "volumeMultiplier",
      "volMult",
      MFnNumericData::kFloat,
      1.0,
      &status);
  McheckErr(status, "ERROR creating attribute PiesVolumeNode::volMultiplier.");

  volMultiplierAttr.setWritable(true);
  volMultiplierAttr.setKeyable(false);
  status = addAttribute(PiesVolumeNode::volMultiplier);
  McheckErr(status, "ERROR adding attribute PiesVolumeNode::volMultiplier.");

  MFnTypedAttribute inMeshAttr;
  PiesVolumeNode::inMesh = inMeshAttr.create(
      "inMesh",
      "inMesh",
      MFnData::kMesh,
      MObject::kNullObj,
      &status);
  McheckErr(status, "ERROR creating attribute PiesVolumeNode::inMesh.");

  inMeshAttr.setWritable(true);
  inMeshAttr.setKeyable(false);
  status = addAttribute(PiesVolumeNode::inMesh);
  McheckErr(status, "ERROR adding attribute PiesVolumeNode::inMesh.");

  MFnTypedAttribute outMeshAttr;
  PiesVolumeNode::outMesh = outMeshAttr.create(
      "outMesh",
      "outMesh",
      MFnData::kMesh,
      MObject::kNullObj,
      &status);
  McheckErr(status, "ERROR creating attribute PiesVolumeNode::outMesh.");

  outMeshAttr.setReadable(true);
  outMeshAttr.setWritable(false);
  outMeshAttr.setStorable(true);
  outMeshAttr.setCached(true);
  status = addAttribute(PiesVolumeNode::outMesh);
  McheckErr(status, "ERROR adding attribute PiesVolumeNode::outMesh.");

  // Setup input-output dependencies
  status = attributeAffects(inMesh, outMesh);
  McheckErr(status, "ERROR attributeAffects(inMesh, outMesh).");

  status = attributeAffects(strainStiffness, outMesh);
  McheckErr(status, "ERROR attributeAffects(strainStiffness, outMesh).");

  status = attributeAffects(minStrain, outMesh);
  McheckErr(status, "ERROR attributeAffects(minStrain, outMesh).");

  status = attributeAffects(maxStrain, outMesh);
  McheckErr(status, "ERROR attributeAffects(maxStrain, outMesh).");

  status = attributeAffects(preserveVolume, outMesh);
  McheckErr(status, "ERROR attributeAffects(preserveVolume, outMesh).");

  status = attributeAffects(volStiffness, outMesh);
  McheckErr(status, "ERROR attributeAffects(volStiffness, outMesh).");

  status = attributeAffects(volMultiplier, outMesh);
  McheckErr(status, "ERROR attributeAffects(volMultiplier, outMesh).");

  return status;
}