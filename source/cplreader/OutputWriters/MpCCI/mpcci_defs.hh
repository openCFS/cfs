<<<<<<< .working
#ifndef FILE_MPCCIDEFS_2006
#define FILE_MPCCIDEFS_2006

// #undef MpCCI
#ifdef MpCCI
#include "mpcci.h"
#else
 /**
  *   *  Element types:
  *     */
#define CCI_ELEM_TRIANGLE            -1500
#define CCI_ELEM_TRIANGLE6           -1502
#define CCI_ELEM_SPH_TRIANGLE        -1504
#define CCI_ELEM_SPH_QUAD            -1506
#define CCI_ELEM_QUAD                -1510
#define CCI_ELEM_QUAD8               -1512
#define CCI_ELEM_QUAD9               -1514
#define CCI_ELEM_PENTAGON            -1520
#define CCI_ELEM_HEXAGON             -1530
#define CCI_ELEM_LINE                -1540
#define CCI_ELEM_TETRAHEDRON         -1550
#define CCI_ELEM_HEXAHEDRON          -1560
#define CCI_ELEM_PYRAMID             -1570
#define CCI_ELEM_PRISM               -1580

#define CCI_CONTINUE                     1
#endif

#endif
=======
// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef FILE_MPCCIDEFS_2006
#define FILE_MPCCIDEFS_2006

// #undef MpCCI
#ifdef MpCCI
#include "mpcci.h"
#else
 /**
  *   *  Element types:
  *     */
#define CCI_ELEM_TRIANGLE            -1500
#define CCI_ELEM_TRIANGLE6           -1502
#define CCI_ELEM_SPH_TRIANGLE        -1504
#define CCI_ELEM_SPH_QUAD            -1506
#define CCI_ELEM_QUAD                -1510
#define CCI_ELEM_QUAD8               -1512
#define CCI_ELEM_QUAD9               -1514
#define CCI_ELEM_PENTAGON            -1520
#define CCI_ELEM_HEXAGON             -1530
#define CCI_ELEM_LINE                -1540
#define CCI_ELEM_TETRAHEDRON         -1550
#define CCI_ELEM_HEXAHEDRON          -1560
#define CCI_ELEM_PYRAMID             -1570
#define CCI_ELEM_PRISM               -1580

#define CCI_CONTINUE                     1
#endif

#endif
>>>>>>> .merge-right.r8654
