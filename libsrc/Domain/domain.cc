#include <stdlib.h>
#include <iostream.h>
#include <fstream.h>
#include <string>

#include <general_head.hh>
#include <utils_head.hh>
#include "grid_cfs.hh"
#include "grid.hh" 
#include "grid_interface_cfs.hh"
#include "domain.hh"

namespace CoupledField
{

template <class Dim>
Domain<Dim> :: Domain(Integer anumsubdomain, FileType * const aptFileType)
{
#ifdef TRACE
  (*trace) << "entering Domain::Domain" << endl;
#endif
  
 numsubdomain = anumsubdomain;
//  grid=new GridInterfaceCFS<Dim>(aptFileType);
 
#ifdef DEBUG
 grid->PrintCoordinate(0,debug);
#endif
}

template <class Dim>
Domain<Dim> :: ~Domain()
{
#ifdef TRACE
  (*trace) << "entering Domain::~Domain" << endl;
#endif

  ;
}

template<class Dim>
void Domain<Dim> :: SetSubdomains()
{
#ifdef TRACE
  (*trace) << "entering Domain::SetSubdomains" << endl;
#endif

}

template<class Dim>
void Domain<Dim> :: PrintDomain()
{
#ifdef TRACE
  (*trace) << "entering Domain::PrintDomain" << endl;
#endif
}

}
