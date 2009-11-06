/*
 * SnOptBase.cc
 *
 *  Created on: Sep 17, 2009
 *      Author: fschury
 */
#include "Optimization/SnOptInterface.hh"
#include "Optimization/SnOptBase.hh"
#include "DataInOut/Logging/cfslog.hh"

namespace CoupledField
{
// declare class specific logging stream
DECLARE_LOG(snoptbase)
DEFINE_LOG(snoptbase, "snoptbase")

/** global pointer to the snoptbase base class */
SnOptBase* static_snopt_base = NULL;

SnOptBase::SnOptBase() :
  f_evals(0),       // number of function evaluations
  g_evals(0),       // number of gradient evaluations
  nxname(1),
  nFname(1),
  iSumm(6),         // variable for snoptbase output
  iPrint(1),       // variable for file output
  iSpecs(4),        // variable for file output
  DerOpt(1),
  minrw(500),
  miniw(500),
  mincw(500),
  lenrw(500),
  leniw(500),
  lencw(500),
  n(0),              // number of variables
  nF(0),
  nG(0),
  lenG(1),
  nA(0),
  lenA(1),
  nS(0),
  nInf(0),
  sInf(0.0),
  Start(0),          // start with a cold start
  stop(false),
  INFO(0),
  ObjRow(1),      // objective row, must be 1!!
  ObjAdd(0.0)
{
  static_snopt_base = this;
}

SnOptBase::~SnOptBase()
{ }

void SnOptBase::SetIntegerValue(const std::string& key, int value)
{
  // FIXME check INFO, if after setting an option INFO > 0 we have had an error!
  
  string option;
  
  if(key == "major_iterations_limit")
  {
    LOG_TRACE(snoptbase) << "adjusted major iterations limit";
    // limit for inner iterations
    option = "Major iteration limit";
  }
  else if(key == "minor_iterations_limit")
  {
    LOG_TRACE(snoptbase) << "adjusted minor iterations limit";
    // limit for inner iterations
    option = "Minor iterations limit";
  }
  else if(key == "superbasics_limit")
  {
    LOG_TRACE(snoptbase) << "adjusted superbasics limit";
    option = "Superbasics limit";
  }
  else if(key == "timing_level")
  {
    LOG_TRACE(snoptbase) << "adjusted timing level";
    // prints cpu times
    option = "Timing level";
  }
  else if(key == "verify_level")
  {
    LOG_TRACE(snoptbase) << "adjusted verify level";
    // checks gradient information using finite differences
    // verifylvl = -1,...,3
    option = "Verify level";
  }
  else if(key == "major_print_level")
  {
    LOG_TRACE(snoptbase) << "adjusted major print level";
    option = "Major print level";
  }
  else if(key == "minor_print_level")
  {
    LOG_TRACE(snoptbase) << "adjusted minor print level";
    option = "Minor print level";
  }
  else if(key == "maximize")
  {
    LOG_TRACE(snoptbase) << "adjusted maximization";
    // maximize or minimize?
    option = "Maximize";
  }
  
  if(!option.empty())
  {
    // set the option
    snseti_(option.c_str(), &value, &iPrint, &iSumm, &INFO,
        &cw[0], &lencw, &iw[0], &leniw, &rw[0], &lenrw, option.size(), 8*lencw);
  }
  else // nothing found
    EXCEPTION("trying to set an snopt-option that is unknown!");
}

void SnOptBase::SetNumericValue(const std::string& key, double value)
{
  // FIXME check INFO, if after setting an option INFO > 0 we have had an error!
  
  string option;
  
  if(key == "major_optimality_tolerance")
  {
    LOG_TRACE(snoptbase) << "adjusted major optimality tolerance";
    option = "Major optimality tolerance";
  }
  else if(key == "major_feasibility_tolerance")
  {
    LOG_TRACE(snoptbase) << "adjusted major feasibility tolerance";
    option = "Major feasibility tolerance";
  }
  else if(key == "minor_feasibility_tolerance")
  {
    LOG_TRACE(snoptbase) << "adjusted minor feasibility tolerance";
    option = "Minor feasibility tolerance";
  }
  else if(key == "linesearch_tolerance")
  {
    LOG_TRACE(snoptbase) << "adjusted linesearch tolerance";
    option = "Linesearch tolerance";
  }
     
  if(!option.empty())
  {
    // set the option
    snsetr_(option.c_str(), &value, &iPrint, &iSumm, &INFO,
        &cw[0], &lencw, &iw[0], &leniw, &rw[0], &lenrw, option.size(), 8*lencw );
  }
  else // nothing found
    EXCEPTION("trying to set an snopt-option that is unknown!");
}
} // namespace

