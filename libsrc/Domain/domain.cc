#include <iostream>
#include <fstream>
#include <string>

#include "domain.hh"

namespace CoupledField
{

template <class Dim>
Domain<Dim> :: Domain(Integer anumsubdomain, FileType * const aptFileType)
{
#ifdef TRACE
  (*trace) << "entering Domain::Domain" << std::endl;
#endif
  
 numsubdomain = anumsubdomain;
//  grid=new GridInterfaceCFS<Dim>(aptFileType);
 
#ifdef DEBUG
// grid->PrintCoordinate(0,debug);
#endif
}

template <class Dim>
Domain<Dim> :: ~Domain()
{
#ifdef TRACE
  (*trace) << "entering Domain::~Domain" << std::endl;
#endif

  ;
}

template<class Dim>
void Domain<Dim> :: SetSubdomains()
{
#ifdef TRACE
  (*trace) << "entering Domain::SetSubdomains" << std::endl;
#endif

}

template<class Dim>
void Domain<Dim> :: PrintDomain()
{
#ifdef TRACE
  (*trace) << "entering Domain::PrintDomain" << std::endl;
#endif
}

}
