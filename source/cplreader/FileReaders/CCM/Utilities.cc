
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

  void CheckError(CCMIOError &err, std::string operation, bool mayExit) {
    if (!IsError(err)) {
      return;
    }
    std::cout << operation << " (error " << err << "  " << CCMIOErrorNames[err] << ")" << std::endl;
    if (mayExit) {
      exit(1);
    }
  }

  void CheckError(CCMIOError &err, std::string operation) {
    CheckError(err, operation, true);
  }
  
  void CollectChildEntities(CCMIOID& parent, std::vector<CCMIOID>& children, CCMIOEntity entityType) {
    children.clear();
    int i = 0;
    CCMIOID id;
    CCMIOError err = kCCMIONoErr;
    while (!IsError(err)) {
      err = CCMIONextEntity(&err, parent, entityType, &i, &id);
      if (!IsError(err)) {
        children.push_back(id);
      }
    }
  }
  
  void OpenLabel(CCMIOID& id, std::string& label) {
    CCMIOError err = kCCMIONoErr;
    uint length;
    CCMIOEntityLabel(&err, id, &length, NULL);
    char * labelchars = new char[length + 1];
    CCMIOEntityLabel(&err, id, NULL, labelchars);
    label = std::string(labelchars, length);
  }

}

