#include <fstream>
#include <iostream>
#include <string>

#include "basepde.hh"
 
namespace CoupledField
{
BasePDE::BasePDE(FileType * ptFileType, Material * aptMaterial)
{
#ifdef TRACE
  (*trace) << "entering BasePDE::BasePDE" << std::endl;
#endif

  ptMaterial=aptMaterial;
  //  ptFileType->ReadIntegrationParam(alpha, beta, gamma_hyperbolic, gamma_parabolic);
 
}

void BasePDE :: CalcIntegrationParam(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering PDE::CalcIntegrationParam" << std::endl;
#endif

  cinteg[1]  = alpha;
  cinteg[2]  = 1.0 - alpha;
  cinteg[3]  = cinteg[2]*dt;
  cinteg[4]  = 0.5*cinteg[3]*(1.0 - 2*beta)*dt;
  cinteg[5]  = cinteg[3]*(1.0 - gamma_hyperbolic);
  cinteg[6]  = cinteg[3]*(1.0 - gamma_parabolic);
  cinteg[7]  = cinteg[3]*beta*dt;
  cinteg[8]  = cinteg[3]*gamma_hyperbolic;
  cinteg[9]  = cinteg[3]*gamma_parabolic;
  cinteg[10] = 1.0/cinteg[2];
  cinteg[11] = 1.0/cinteg[7];
  cinteg[12] = gamma_hyperbolic/(cinteg[3]*beta);
}

} // end of namespace
