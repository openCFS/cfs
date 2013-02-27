/*
 * MortarInterface.hh
 *
 *  Created on: 23.02.2013
 *      Author: jens
 */

#ifndef _MORTARINTERFACE_HH_
#define _MORTARINTERFACE_HH_

#include "BaseNcInterface.hh"
#include "DataInOut/ParamHandling/ParamNode.hh"

namespace CoupledField {

class MortarInterface : public BaseNcInterface {
    
  public:
    
    MortarInterface(Grid* grid, PtrParamNode nciNode);

    virtual ~MortarInterface();
    
    void UpdateInterface();
    
    bool IsPlanar() const { return coplanar_; }

  protected:
    
    bool coplanar_;

};

} /* namespace CoupledField */
#endif /* _MORTARINTERFACE_HH_ */
