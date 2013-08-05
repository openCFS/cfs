/* 
 * File:   FieldDataAcceptor.h
 * Author: Mace
 *
 * Created on 1. Juni 2013, 23:55
 */

#ifndef FIELDDATAACCEPTOR_H
#define	FIELDDATAACCEPTOR_H

#include "Utilities.hh"

namespace CCM {
  
  class FieldDataAcceptor {
  public:
    FieldDataAcceptor();
    FieldDataAcceptor(const FieldDataAcceptor& orig);
    virtual ~FieldDataAcceptor();
    
    virtual bool Accept(FieldData& data, bool& alreadyAllocated);
  private:
  
  };
    
}

#endif	/* FIELDDATAACCEPTOR_H */

