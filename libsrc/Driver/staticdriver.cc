#include <fstream>
#include <iostream>
#include <string>

#include "staticdriver.hh"

namespace CoupledField
{

template<class Dim>
StaticDriver<Dim> :: StaticDriver(Domain<Dim> * adomain)
:BaseDriver<Dim>(adomain)
{
#ifdef TRACE
  (*trace) << "entering StaticDriver::StaticDriver" << std::endl;
#endif

}

template<class Dim>
StaticDriver<Dim> :: ~StaticDriver()
{
#ifdef TRACE
  (*trace) << "entering StaticDriver::~StaticDriver" << std::endl;
#endif

}

template<class Dim>
void StaticDriver<Dim> :: SetupMatricesPDE(const Integer pdenumber, const Integer type)
{
#ifdef TRACE
  (*trace) << "entering StaticDriver::SetUpMatricesPDE" << std::endl;
#endif

  ptdomain_->GetPDE(pdenumber)->SetupMatrices(type);
}

template<class Dim>
void StaticDriver<Dim> :: SolveProblem()
{
#ifdef TRACE
  (*trace) << "entering StaticDriver::SolveProblem" << std::endl;
#endif
  Integer level=0;
  Integer pdenumber  = 0;
  Integer matrixtype = 0;
  Integer nsys = 0;

  ptdomain_->GetPDE(pdenumber)->SolveStepStatic(ptdomain_->GetBCs(), level);

  ptdomain_->GetPDE(pdenumber)->WriteResultsInFile();
}

}
