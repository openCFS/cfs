#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <math.h>
#include <limits.h>
#include <string>

#include "MaterialData.hh"


namespace CoupledField
{

MaterialData::MaterialData():scaledMatDat(0)
{
#ifdef TRACE
  if (trace) (*trace) << "Entering MaterialData::MaterialData" << std::endl;
#endif
  const int stringLength = 100;
  name = new char[stringLength];  
}


void MaterialData::SetConductivity(const double& Conductivity) 
{
#ifdef TRACE
  if (trace) (*trace) << "Entering MaterialData::SetConductivity" << std::endl;
#endif
  conductivity = Conductivity;
}

void MaterialData::SetName(const char* Name)
{
#ifdef TRACE
  if (trace) (*trace) << "Entering MaterialData::SetName" << std::endl;
#endif

  strcpy(name,Name);
}

char * MaterialData::GetMaterialName()
{
#ifdef TRACE
  if (trace) (*trace) << "Entering MaterialData::GetMaterialName" << std::endl;
#endif
  return name;
}


/*
 
void MaterialData::DefLin(const int& notLin)
{

  if (notLin)
    nonlin = 1;
  else
    if (magneticSpline->Size()) // exists magneticSpline already ? 
      constPerm = magneticSpline->Get(1) -> ConstPerm(); 
    else
      cout << "You have to define the magneticSpline before the nonlin variable (in the nonlin case)!" << endl;
}


double MaterialData::GetPermiability(const double& MagB) const
{
  int i;

  i = 1;

    while (magneticSpline->Get(i)->UpperLimit() < MagB)
    {
    i++;
    if (i>magneticSpline->Size())
    {
    cout << endl << " Induction " << MagB << " too high " << endl;
    exit(EXIT_FAILURE);
    }
    }
    
    if (i == 1)
    return magneticSpline->Get(i)->ReducedValue(MagB); //at the first spline no offset is tolerated!
    else
    return magneticSpline->Get(i)->Value(MagB);
}


double MaterialData::GetErrorFunctional(const double& MagB) const
{
  int i;
  double errorfunctional=0;

  i = 1;
  while (magneticSpline->Get(i)->UpperLimit() < MagB)
    {
      if (i==magneticSpline->Size())
	{
	  cout << endl << " Induction too high ";
	  exit(EXIT_FAILURE);
	}
      errorfunctional += magneticSpline->Get(i) -> WholeIntegral();
      i++;
    }
  return errorfunctional += magneticSpline->Get(i)->IntValue(MagB);  

}
 */
 
}
