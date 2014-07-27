#ifndef DYNAMICOBJECT_H
#define DYNAMICOBJECT_H

namespace CoupledField {

  class DynamicObject
  {
  private:
    // Callback function that should be called to delete dynamic object
    void (*_deleteObject)(void*);
  public:
    // The constructor sets the callback function to use
    DynamicObject(void (*delObj)(void*));
    
    // The destructor
    virtual ~DynamicObject(void);
    
    // Sends "this" to the callback destructor function.
    void deleteSelf(void);
  };
  
  DynamicObject::DynamicObject(void (*delObj)(void*))
    : _deleteObject(delObj)
  {
  }
  
  DynamicObject::~DynamicObject(void)
  {
  }
  
  void
  DynamicObject::deleteSelf(void)
  {
    (*_deleteObject)(reinterpret_cast<void*>(this));
  }
  
}

#endif
