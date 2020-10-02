// -*- mode: c++; coding: utf-8; indent-tabs-mode: nil; -*-
// kate: space-indent on; indent-width 2; encoding utf-8;
// kate: auto-brackets on; mixedindent off; indent-mode cstyle;

#ifndef _GMVREADH_
#define _GMVREADH_


// This code is also based on the gmvread library from the official GMV
// website at: http://www-xdiv.lanl.gov/XCM/gmv/GMVHome.html

#include <vector>

#include "DataInOut/Logging/LogConfigurator.hh"

namespace CoupledField {

#ifndef RDATA_INIT
#define EXTERN extern
#else
#define EXTERN /**/
#endif

/*  Keyword types.  */
#define RAYS        1
#define RAYIDS      2
#define RAYEND     51

#define NRAYVARS 20

struct gmvray
{
  int    npts;
  double *x;
  double *y;
  double *z;
  double *field[NRAYVARS];
};

struct gmvray_data_struct 
{
  int     nrays;    /*  Number of rays in the file.  */
  int     nvars;    /*  Number of ray variable fields.  */
  char    *varnames;  /*  Ray variable names.  */
  short   vartype[NRAYVARS];   /*  Variable types.  */
  /*  0 = nodes, 1 = line segments.  */
  int     *rayids;
  struct gmvray *gmvrays;
} EXTERN gmvray_data;

/*  C, C++ prototypes.  */

int gmvrayread_open(char *filnam);

void gmvrayread_close(void);

void gmvrayread_data(void);


  /*  Keyword types.  */
#define NODES       1
#define CELLS       2
#define FACES       3
#define VFACES      4
#define XFACES      5
#define MATERIAL    6
#define VELOCITY    7
#define VARIABLE    8
#define FLAGS       9
#define POLYGONS   10
#define TRACERS    11
#define PROBTIME   12
#define CYCLENO    13
#define NODEIDS    14
#define CELLIDS    15
#define SURFACE    16
#define SURFMATS   17
#define SURFVEL    18
#define SURFVARS   19
#define SURFFLAG   20
#define UNITS      21
#define VINFO      22
#define TRACEIDS   23
#define GROUPS     24
#define FACEIDS    25
#define SURFIDS    26
#define CELLPES    27
#define SUBVARS    28
#define GHOSTS     29
#define CODENAME   48
#define CODEVER    49
#define VECTOR     50
#define SIMDATE    51
#define GMVEND     52
#define INVALIDKEYWORD 53
#define GMVERROR 54


  /*  Data types for Nodes:  */
#define UNSTRUCT  100
#define STRUCT 101
#define LOGICALLY_STRUCT 102
#define AMR 103
#define VFACES2D 104
#define VFACES3D 105
#define NODE_V   106

  /*  Data types for Cells:  */
#define GENERAL 110
#define REGULAR 111
#define VFACE2D 112
#define VFACE3D 113

  /*  Data types for vectors, variables, materials, flags, tracers, groups: */
#define NODE       200
#define CELL       201
#define FACE       202
#define SURF       203
#define XYZ        204
#define TRACERDATA 205
#define VEL        206
#define ENDKEYWORD 207
#define FROMFILE   208


  enum CellType
  {
    TRI = 1,
    QUAD = 2,
    TET = 3,
    HEX = 4,
    PRISM = 5,
    PYRAMID = 6,
    LINE = 7,
    PHEX8 = 8,
    PHEX20 = 9,
    PPYRMD5 = 10,
    PPYRMD13 = 11,
    PPRISM6 = 12,
    PPRISM15 = 13,
    PTET4 = 14,
    PTET10 = 15,
    TRI6 = 16,
    QUAD8 = 17,
    LINE3 = 18,
    PHEX27 = 19
  };

  extern std::vector<CellType> gmv_cell_types;
  extern std::string gmv_base_dir;
  extern std::vector< std::string > gmv_vector_components;

  class GMVReadException
  {
  public:
    GMVReadException() {}
    ~GMVReadException() {}
  };

  struct gmv_data_struct
  {
    int     keyword;    /*  See above for definitions.  */
    int     datatype;   /*  See above for definitions.  */
    char    name1[33];  /*  hex, tri, etc, flag name, field name.  */
    long    num;        /*  nnodes, ncells, nsurf, ntracers.  */
    long    num2;       /*  no. of faces, number of vertices.  */

    long    ndoubledata1;
    double  *doubledata1;
    long    ndoubledata2;
    double  *doubledata2;
    long    ndoubledata3;
    double  *doubledata3;

    long    nlongdata1;
    long    *longdata1;
    long    nlongdata2;
    long    *longdata2;

    int     nchardata1;   /*  Number of 33 character string.  */
    char    *chardata1;   /*  Array of 33 character strings.  */
    int     nchardata2;   /*  Number of 33 character string.  */
    char    *chardata2;   /*  Array of 33 character strings.  */
  } EXTERN gmv_data;


  struct gmv_meshdata_struct
  {
    long    nnodes; 
    long    ncells;
    long    nfaces;
    long    totfaces;
    long    totverts;
    int     intype;  /* CELLS, FACES, STRUCT, LOGICALLY_STRUCT, AMR. */
    int     nxv;  /*  nxv, nyv, nzv for STRUCT,  */
    int     nyv;  /*  LOGICALLY_STRUC and AMR.   */
    int     nzv;

    double  *x;  /*  Node x,y,zs, nnodes long.  */
    double  *y;
    double  *z;

    long    *celltoface;   /*  Cell to face pointer, ncells+1 long. */
    long    *cellfaces;    /*  Faces in cells, totfaces+1 long.  */
    long    *facetoverts;  /*  Face to verts pointer, nfaces+1 long.*/
    long    *faceverts;    /*  Verts per face, totverts long.  */
    long    *facecell1;    /*  First cell face attaches to.  */
    long    *facecell2;    /*  Second cell, nfaces long.  */
    long    *vfacepe;      /*  Vface pe no.  */
    long    *vfaceoppface;  /*  Vface opposite face no.  */
    long    *vfaceoppfacepe;  /*  Vface opposite face pe no.  */
    long    *cellnnode;    /*  No. of nodes per cell (regular cells).  */
    long    *cellnodes;    /*  Node list per cell (regular cells).  */
  } EXTERN gmv_meshdata;

  /*  C, C++ prototypes.  */

  int gmvread_checkfile(const char *filnam);

  int gmvread_open(const char *filnam);

  int gmvread_open_fromfileskip(const char *filnam);

  void gmvread_close(void);

  void gmvread_data(void);

  void gmvread_mesh(void);

  void gmvread_printon();

  void gmvread_printoff();

  void struct2face(void);

  void struct2vface(void);


  // helper functions

  int chk_gmvend(FILE *fin);
  void fromfilecheck(int fkeyword);
  void rdxfaces();
  void checkfromfile();
  void rdfaces();
  void rdcells(int nodetype_in);
  void fillmeshdata(long nc);
  void rdvfaces(long nc);
  void vfacecell(long icell, long nc);
  void regcell(long icell, long nc);
  void gencell(long icell, long nc);
  void fillcellinfo(long nc, long *facecell1, long *facecell2);
  int chk_rayend(FILE *fin);

}

#endif
