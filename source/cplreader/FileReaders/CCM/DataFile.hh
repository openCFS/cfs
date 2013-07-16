/* 
 * File:   CCMDataFile.h
 * Author: tz
 *
 * Created on 31. Mai 2013, 10:29
 */

#ifndef DATAFILE_H
#define	DATAFILE_H

#include <vector>
#include <string>

#include "LibCCMIO/ccmio.h"
#include "LibCCMIO/ccmiocore.h"
#include "LibCCMIO/ccmiotypes.h"

#include "Utilities.hh"
#include "FieldDataAcceptor.hh"

namespace CCM {
  
  class DataFile {
    public:
      DataFile(std::string fileName);
      DataFile(std::string fileName, bool verbose);
      DataFile(const DataFile& orig);
      virtual ~DataFile();
      
      void ReadMaps(std::vector<Map>* maps);
      void ReadFields(std::vector<FieldData>* data, FieldDataAcceptor* acceptor);
      void ReadFields(std::vector<FieldData>* data, std::string fieldLabel);
      
      bool Close();
      bool IsOpen();
      
      void SetVerbose(bool verbose);
      
      
    protected:
      void Open();
      void CheckError(std::string job);
      void ResetError();
      
      uint GetMapSize(CCMIOID mapID);
      
      bool isOpen_;
      CCMIOError err_;
      CCMIOID root_;
      std::string fileName_;
      bool verbose_;
  };

}

#endif	/* DATAFILE_H */

