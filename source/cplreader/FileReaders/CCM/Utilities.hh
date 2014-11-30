/* 
 * File:   Utilities.h
 * Author: Mace
 *
 * Created on 1. Juni 2013, 03:04
 */

#ifndef UTILITIES_H
#define	UTILITIES_H

#include <vector>
#include <string>

#include "ccmio.h"
#include "ccmiocore.h"
#include "ccmiotypes.h"

#include "General/environment.hh"


namespace CCM {

  bool IsError(CCMIOError& err);

  void CheckFileError(CCMIOError &err, std::string operation, bool mayExit);

  void CheckFileError(CCMIOError &err, std::string operation);
  
  void CollectChildEntities(CCMIOID& parent, std::vector<CCMIOID>& children, CCMIOEntity entityType);

  void OpenLabel(CCMIOID& id, std::string& label);
  
  typedef enum { DimNullData = 0, ScalarData = 1, VectorData, TensorData } DataDimensionality;
  typedef enum { VertexDataLoc = 0, CellDataLoc, FaceDataLoc } DataLocation;
  typedef enum { FloatData = 0, DoubleData, IntData, OtherData} DataType;

  struct Data {
    int id;
    std::string label;
  };
  
  struct Map : public Data {
    uint size;
    uint maxValue;
    int* value;
  };
  
  struct MappedData : public Data {
    int mapID;
    std::string mapLabel;
  };

  struct FieldData : public MappedData {
    DataDimensionality dimension;
    DataLocation location;
    DataType type;
    uint size;
    double* doubleData;
    float* floatData;
    int* intData;
    
    int fieldSetID;
    std::string fieldSetLabel;
    int fieldPhaseID;
    std::string fieldPhaseLabel;
    int fieldID;
    std::string fieldLabel;
    std::string fieldLabelShort;
  };
  
  struct FaceData : public MappedData {
    bool isBoundary;
    uint size;
    int* faceCells;
    uint faceVerticesSize;
    int* faceVertices;
    
    int topolgyID;
    std::string topologyLabel;
  };
  
  struct VertexData : public MappedData {
    uint size;
    int dim;
    float scale;
    double* coordinates;
  };
  
  struct CellData : public MappedData {
    uint size;
    
    int topolgyID;
    std::string topologyLabel;
  };
  
  
}


#endif	/* UTILITIES_H */

