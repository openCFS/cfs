#ifndef FILE_PDES_HEADER_2003
#define FILE_PDES_HEADER_2003


#ifdef NEWBASEPDE

#include "newElecPDE.hh"
#include "newMechPDE.hh"
#include "newAcousticPDE.hh"
#else // NEWBASEPDE

#include "acousticPDE.hh"
#include "acouFlowNoise.hh"
#include "mechPDE.hh"
#include "magEdgePDE.hh"
#include "smoothPDE.hh"
#include "smoothlaplacePDE.hh"

#include "elecPDE.hh"

#endif // NEWBASEPDE


#endif
