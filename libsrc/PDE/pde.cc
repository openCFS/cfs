#include <fstream>
#include <iostream>
#include <string>

#include <general_head.hh> 
#include <utils_head.hh>
#include <datainout_head.hh>
#include <domain_head.hh>
#include <linalg_head.hh>
#include "pde.hh"
 
namespace CoupledField
{
PDE::PDE(FileType * ptFileType, Material * aptMaterial)
{
#ifdef TRACE
  (*trace) << "entering PDE::PDE" << std::endl;
#endif

  ptMaterial=aptMaterial;
  ptFileType->ReadIntegrationParam(alpha, beta, gamma);
 
}

void PDE :: CalcParamForNewmarkMethod(const Double dt)
{
#ifdef TRACE
  (*trace) << "entering PDE::CalcParamForNewmarkMethod" << std::endl;
#endif
 
 a2=1.0/(beta*dt);
 a0=a2*(1/dt);
 a1=gamma*a2;
 a3=0.5/(beta)-1.0;
 a4=gamma/beta-1.0;
 a5=0.5*dt*(a4-1.0);
 a6=dt*(1-gamma);
 a7=gamma*dt;
}

} // end of namespace
