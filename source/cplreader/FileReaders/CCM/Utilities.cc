
#include <iostream>
#include <stdlib.h>

#include "Utilities.hh"

#include "General/environment.hh"


namespace CCM {
  
  static const char * CCMIOErrorNames[] = {
    "kCCMIONoErr",
    "kCCMIONoFileErr",
    "kCCMIOPermissionErr",
    "kCCMIOCorruptFileErr",
    "kCCMIOBadLinkErr",
    "kCCMIONoNodeErr",
    "kCCMIODuplicateNodeErr",
    "kCCMIOWrongDataTypeErr",
    "kCCMIONoDataErr",
    "kCCMIOWrongParentErr",
    "kCCMIOBadParameterErr",
    "kCCMIONoMemoryErr",
    "kCCMIOIOErr",
    "kCCMIOTooManyFacesErr",
    "kCCMIOVersionErr",
    "kCCMIOArrayDimensionToLargeErr",
    "kCCMIOInternalErr"
  };

  bool IsError(CCMIOError& err) {
    return err != kCCMIONoErr;
  }

  void CheckFileError(CCMIOError &err, std::string operation, bool mayExit) {
    if (!IsError(err)) {
      return;
    }
    std::cerr << "Error while reading ccm file: " << operation << " (error " << err << "  " << CCMIOErrorNames[err] << ")" << std::endl;
    if (mayExit) {
      exit(1);
    }
  }

  void CheckFileError(CCMIOError &err, std::string operation) {
    CheckFileError(err, operation, true);
  }
  
  void CollectChildEntities(CCMIOID& parent, std::vector<CCMIOID>& children, CCMIOEntity entityType) {
    children.clear();
    CCMIOSize_t i = CCMIOSIZEC(0);
    CCMIOID id;
    CCMIOError err = kCCMIONoErr;
    while (!IsError(err)) {
      err = CCMIONextEntity(&err, parent, entityType, (CCMIOSize_t*) &i, &id);
      if (!IsError(err)) {
        children.push_back(id);
      }
    }
  }
  
  void OpenLabel(CCMIOID& id, std::string& label) {
    CCMIOError err = kCCMIONoErr;
    CCMIOSize_t length;
    CCMIOEntityLabel(&err, id, (CCMIOSize_t*) &length, NULL);
    char * labelchars = new char[TOINT64(length) + 1];
    CCMIOEntityLabel(&err, id, NULL, labelchars);
    label = std::string(labelchars, TOINT64(length));
  }

}

