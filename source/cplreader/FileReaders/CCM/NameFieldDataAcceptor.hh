/* 
 * File:   NameFieldDataAcceptor.h
 * Author: Mace
 *
 * Created on 2. Juni 2013, 00:36
 */

#ifndef NAMEFIELDDATAACCEPTOR_H
#define	NAMEFIELDDATAACCEPTOR_H

#include <string>

#include "FieldDataAcceptor.hh"

namespace CCM {
  
  class NameFieldDataAcceptor : public FieldDataAcceptor {
  public:
    NameFieldDataAcceptor(std::string name);
    NameFieldDataAcceptor(const NameFieldDataAcceptor& orig);
    virtual ~NameFieldDataAcceptor();
    bool Accept(FieldData& data, bool& alreadyAllocated);
  private:
    std::string name_;
  };
  
}

#endif	/* NAMEFIELDDATAACCEPTOR_H */

