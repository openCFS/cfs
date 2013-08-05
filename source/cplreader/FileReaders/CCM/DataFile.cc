/* 
 * File:   DataFile.cpp
 * Author: tz
 * 
 * Created on 31. Mai 2013, 10:29
 */

#include <vector>
#include <string>
#include <iostream>

#include "DataFile.hh"

#include "ccmio.h"
#include "ccmiocore.h"
#include "ccmiotypes.h"

#include "Utilities.hh"
#include "FieldDataAcceptor.hh"
#include "NameFieldDataAcceptor.hh"

namespace CCM {
  
  DataFile::DataFile(std::string fileName):
        isOpen_(false),err_(kCCMIONoErr),fileName_(fileName),verbose_(true) {
    Open();
  }
  
  DataFile::DataFile(std::string fileName, bool verbose):
        isOpen_(false),err_(kCCMIONoErr),fileName_(fileName),verbose_(verbose) {
    Open();
  }
  
  DataFile::DataFile(const DataFile& orig) {
      std::cerr << "DataFile: Copy constructor not implemented" << std::endl;
  }
  
  void DataFile::SetVerbose(bool verbose) {
    verbose_ = verbose;
  }
  
  bool DataFile::GetVerbose() {
    return verbose_;
  }
  
  bool DataFile::IsOpen() {
    return isOpen_;
  }
  
  void DataFile::ReadMaps(std::vector<Map>* maps) {
    std::cout << "Reading Maps:" << std::endl;
    Map map;
    std::vector<CCMIOID> mapIDs;
    CollectChildEntities(root_, mapIDs, kCCMIOMap);
    for (uint i=0; i < mapIDs.size(); i++) {
      CCMIOID ccmid = mapIDs[i];
      map.id = ccmid.id;
      OpenLabel(ccmid, map.label);
      CCMIOSize_t maxValue;
      CCMIOSize_t size;
      CCMIOEntitySize(&err_, ccmid, &size, &maxValue);
      map.maxValue = TOINT64(maxValue);
      map.size = TOINT64(size);
      std::cout << "  " << map.label << " (size: "  << map.size << ", maxvalue: " << map.maxValue << ")" << std::endl;
      map.value = new int[map.size];
      CCMIOReadMap(NULL, ccmid, map.value, CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));
      maps->push_back(map);
    }
  }
  
  void DataFile::ReadFields(std::vector<FieldData>* data, std::string name) {
    NameFieldDataAcceptor acceptor(name);
    ReadFields(data, &acceptor);
  }
  
  void DataFile::ReadFields(std::vector<FieldData>* data, FieldDataAcceptor* acceptor) {
    if (verbose_) {
      std::cout << "Loading Scalar Fields:" << std::endl;
    }
    
    FieldData scalarData;
    
    std::vector<CCMIOID> fieldSetIDs;
    CollectChildEntities(root_, fieldSetIDs, kCCMIOFieldSet);
    for (uint iSet = 0; iSet < fieldSetIDs.size(); iSet++) {
      CCMIOID fieldSetID = fieldSetIDs[iSet];
      scalarData.fieldSetID = fieldSetID.id;
      OpenLabel(fieldSetID, scalarData.fieldSetLabel);
      if (verbose_) {
        std::cout << "  Field Set: " << scalarData.fieldSetLabel << " (id: " << scalarData.fieldSetID << ")" << std::endl;
      }
      
      std::vector<CCMIOID> fieldPhaseIDs;
      CollectChildEntities(fieldSetID, fieldPhaseIDs, kCCMIOFieldPhase);
      for (uint iPhase = 0; iPhase < fieldPhaseIDs.size(); iPhase++) {
        CCMIOID fieldPhaseID = fieldPhaseIDs[iPhase];
        scalarData.fieldPhaseID = fieldPhaseID.id;
        OpenLabel(fieldPhaseID, scalarData.fieldPhaseLabel);
        if (verbose_) {
          std::cout << "    Field Phase: " << scalarData.fieldPhaseLabel << " (id: " << scalarData.fieldPhaseID << ")" << std::endl;
        }
        
        std::vector<CCMIOID> fieldIDs;
        CollectChildEntities(fieldPhaseID, fieldIDs, kCCMIOField);
        for (uint iField = 0; iField < fieldIDs.size(); iField++) {
          CCMIOID fieldID = fieldIDs[iField];
          scalarData.fieldID = fieldID.id;
          //OpenLabel(fieldID, scalarData.fieldLabel);
          
          // getting dimension, datatype and shortname
          CCMIODimensionality dims;
          CCMIODataType dataType;
          char shortName[kCCMIOProstarShortNameLength+1];
          char name[100];
          CCMIOReadField(&err_, fieldID, name, shortName, &dims, &dataType);
          scalarData.fieldLabel = std::string(name);
          if (dims == kCCMIOScalar) {
            scalarData.dimension = ScalarData;
            if (verbose_) {
              std::cout << "dim: scalar, ";
            }
          } else if (dims == kCCMIOVector) {
            scalarData.dimension = VectorData;
            if (verbose_) {
              std::cout << "dim: vector, ";
            }
          } else if (dims == kCCMIOTensor) {
            scalarData.dimension = TensorData;
            if (verbose_) {
              std::cout << "dim: tensor, ";
            }
          } else if (dims == kCCMIODimNull) {
            scalarData.dimension = DimNullData;
            if (verbose_) {
              std::cout << "dim: null, ";
            }
          }
          if (dataType == kCCMIOFloat64) {
            scalarData.type = DoubleData;
            if (verbose_) {
              std::cout << "type: double, ";
            }
          } else if (dataType == kCCMIOFloat32) {
            scalarData.type = FloatData;
            if (verbose_) {
              std::cout << "type: float, ";
            }
          } else if (dataType == kCCMIOInt32) {
            scalarData.type = IntData;
            if (verbose_) {
              std::cout << "type: int, ";
            }
          } else {
            scalarData.type = OtherData;
            if (verbose_) {
              std::cout << "type: other, ";
            }
          }
          scalarData.fieldLabelShort = std::string(shortName);
          if (verbose_) {
            std::cout << "      Field: " << scalarData.fieldLabel << " (id: " << scalarData.fieldID << ", ";
            std::cout << "short: " << scalarData.fieldLabelShort << ")" << std::endl;
          }
          if (dims == kCCMIOScalar) {
            scalarData.dimension = ScalarData;
            if (verbose_) {
              std::cout << "dim: scalar, ";
            }
          } else if (dims == kCCMIOVector) {
            scalarData.dimension = VectorData;
            if (verbose_) {
              std::cout << "dim: vector, ";
            }
          } else if (dims == kCCMIOTensor) {
            scalarData.dimension = TensorData;
            if (verbose_) {
              std::cout << "dim: tensor, ";
            }
          } else if (dims == kCCMIODimNull) {
            scalarData.dimension = DimNullData;
            if (verbose_) {
              std::cout << "dim: null, ";
            }
          }
          if (dataType == kCCMIOFloat64) {
            scalarData.type = DoubleData;
            if (verbose_) {
              std::cout << "type: double, ";
            }
          } else if (dataType == kCCMIOFloat32) {
            scalarData.type = FloatData;
            if (verbose_) {
              std::cout << "type: float, ";
            }
          } else if (dataType == kCCMIOInt32) {
            scalarData.type = IntData;
            if (verbose_) {
              std::cout << "type: int, ";
            }
          } else {
            scalarData.type = OtherData;
            if (verbose_) {
              std::cout << "type: other, ";
            }
          }
          
          std::vector<CCMIOID> fieldDataIDs;
          CollectChildEntities(fieldID, fieldDataIDs, kCCMIOFieldData);
          for (uint iFieldData = 0; iFieldData < fieldDataIDs.size(); iFieldData++) {
            CCMIOID fieldDataID = fieldDataIDs[iFieldData];
            scalarData.id = fieldDataID.id;
            OpenLabel(fieldDataID, scalarData.label);
            CCMIOSize_t size;
	    CCMIOEntitySize(&err_, fieldDataID, &size, NULL);
	    scalarData.size = TOINT64(size);
    	      
            CCMIOID mapID;
            CCMIODataLocation loc;
            if (scalarData.type == DoubleData) {
              CCMIOReadFieldDatad(&err_, fieldDataID, &mapID, &loc, NULL, CCMIOINDEXC(0), CCMIOINDEXC(0));
            } else if (scalarData.type == FloatData) {
              CCMIOReadFieldDataf(&err_, fieldDataID, &mapID, &loc, NULL, CCMIOINDEXC(0), CCMIOINDEXC(0));
            } else if (scalarData.type == IntData) {
              CCMIOReadFieldDatai(&err_, fieldDataID, &mapID, &loc, NULL, CCMIOINDEXC(0), CCMIOINDEXC(0));
            }
            if (verbose_) {
              std::cout << "        Field Data: " << scalarData.label << " (id: " << scalarData.id << ", size: " << scalarData.size << ", ";
            }
            if (loc == kCCMIOCell) {
              scalarData.location = CellDataLoc;
              if (verbose_) {
                std::cout << "location: cell)" << std::endl;
              }
            } else if (loc == kCCMIOFace) {
              scalarData.location = FaceDataLoc;
              if (verbose_) {
                std::cout << "location: face)" << std::endl;
              }
            } else if (loc == kCCMIOVertex) {
              scalarData.location = VertexDataLoc;
              if (verbose_) {
                std::cout << "location: vertex)" << std::endl;
              }
            }
            scalarData.mapID = mapID.id;
            OpenLabel(mapID, scalarData.mapLabel);
            if (verbose_) {
              std::cout << "          Map: " << scalarData.mapLabel << " (id: " << scalarData.mapID << ")" << std::endl;
            }
            
            bool alreadyAllocated = false;
            if (acceptor->Accept(scalarData, alreadyAllocated)) {
              std::cout << "  Reading Field Data: " << scalarData.label << std::endl;
              if (!alreadyAllocated) {
                if (scalarData.type == DoubleData) {
                  scalarData.doubleData = new double[scalarData.size];
                } else if (scalarData.type == FloatData) {
                  scalarData.floatData = new float[scalarData.size];
                } else if (scalarData.type == IntData) {
                  scalarData.intData = new int[scalarData.size];
                }
              }
              if (scalarData.type == DoubleData) {
                CCMIOReadFieldDatad(&err_, fieldDataID, NULL, NULL, scalarData.doubleData, CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));
              } else if (scalarData.type == FloatData) {
                CCMIOReadFieldDataf(&err_, fieldDataID, NULL, NULL, scalarData.floatData, CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));
              } else if (scalarData.type == IntData) {
                CCMIOReadFieldDatai(&err_, fieldDataID, NULL, NULL, scalarData.intData, CCMIOINDEXC(kCCMIOStart), CCMIOINDEXC(kCCMIOEnd));
              }
              if (data != NULL) {
                data->push_back(scalarData);
              }
            } else {
              if (verbose_) {
                std::cout << "  Skipping Field Data: " << scalarData.label << std::endl;              
	      }
            }
          }
        }
      }
    }
  }
  
  void DataFile::Open() {
    err_ = CCMIOOpenFile(NULL, fileName_.c_str(), kCCMIORead, &root_);
    std::cout << "Opening File: " << fileName_ << std::endl;
    CheckError("Opening file ");
    isOpen_ = true;
  }

  bool DataFile::Close() {
    if (!isOpen_) {
      return false;;
    }
    CCMIOCloseFile(&err_, root_);
    std::cout << "Closing File " << fileName_ << std::endl << std::endl;
    CheckError("Closing file ");
    isOpen_ = false;
    return true;
  }
  
  void DataFile::CheckError(std::string operation) {
    CheckFileError(err_, operation, false);
  }

  uint DataFile::GetMapSize(CCMIOID mapID) {
    CCMIOSize_t size;
    CCMIOEntitySize(&err_, mapID, (CCMIOSize_t*) &size, NULL);
    return TOINT64(size);
  }
  
  DataFile::~DataFile() {
    Close();
  }
  
  void DataFile::ResetError() {
    err_ = kCCMIONoErr;
  }
  
}

