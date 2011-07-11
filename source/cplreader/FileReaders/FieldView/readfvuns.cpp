#include <iostream>
#include <stdlib.h>
#include <string.h>
#include <vector>
#include <assert.h>
#include <stdio.h>

#include "fv_reader_tags.h"
#include "readfvuns.h"

#define READ(c,n) fread(c, sizeof(*c), n, outfp); if (endian_conv) endian_convert((char*) c, n*sizeof(*c));

// defines for decoding of element headers
#define MAX_NUM_ELEM_FACES  6
#define BITS_PER_WALL       3
#define ELEM_TYPE_BIT_SHIFT (MAX_NUM_ELEM_FACES*BITS_PER_WALL)


#define VERBOSE 0


class Part
{
  public:
    Part()
    {
      num_vars = 0;
      num_bvars = 0;
      nnodes = 0;
      
      xyz = 0;
      
      vars = 0;
      
      num_bndpatches_tri_quad = 0; 
      num_bndpatches_poly = 0;    

      ntet = 0;
      nhex = 0;
      npri = 0;
      npyr = 0;
      tet_ids = 0;
      hex_ids = 0;
      pri_ids = 0;
      pyr_ids = 0;
    }
    
    ~Part()
    {
      if (xyz)
        delete [] xyz;

      if (vars)
        delete [] vars;

      if (face_ids.size())
        for (unsigned int i = 0; i < face_ids.size(); ++i)
          if (face_ids[i])
            delete [] face_ids[i];
      
      if (boundary_vars.size())
        for (unsigned int i = 0; i < boundary_vars.size(); ++i)
          if (boundary_vars[i])
            delete [] boundary_vars[i];
      
      if (tet_ids)
        delete [] tet_ids;
        
      if (hex_ids)
        delete [] hex_ids;
        
      if (pri_ids)
        delete [] pri_ids;
        
      if (pyr_ids)
        delete [] pyr_ids;
    }
    
    int num_vars;  // number of region variables
    int num_bvars; // number of boundary variables
    int nnodes;    // number of grid nodes
    
    float *xyz; // vertex coordinates for this part
    
    float *vars; // variable data for this part
    
    std::vector<int> num_faces; // number of boundary faces of a specific boundary type
    std::vector<int> boundary_types; // boundary type the faces belong to
    std::vector<int*> face_ids; // vertex ids constituting boundary faces
    std::vector<float*> boundary_vars; // variable data on boundaries
    
    int num_bndpatches_tri_quad; // number of boundaries having tri/tet face elements
    int num_bndpatches_poly;     // number of boundaries having polygonal face elements
    
    int ntet; // number of tet cells in this part
    int nhex; // number of hex cells in this part
    int npri; // number of prism cells in this part
    int npyr; // number of pyramid cells in this part
    int *tet_ids; // vertex ids constituting tet cells
    int *hex_ids; // vertex ids constituting hex cells
    int *pri_ids; // vertex ids constituting prism cells
    int *pyr_ids; // vertex ids constituting pyramid cells
};


template <class T>
void smartalloc(T *&ptr, int n)
{
  if (n)
  {
    ptr = new T[n];
    assert(ptr);
  }
  else
    ptr = 0;
}


void endian_convert(char *cs, int n)
{
  for (int i = 0; i < n; ++i)
  {
    char tmp = cs[i];
    cs[i] = cs[n-1-i];
    cs[n-1-i] = tmp;
  }
}

const char *file_name;
int num_grids;
Part *parts;
int * results_flag;
int * normals_flag;

void Initreadfvuns();
void cleanup();
void Initreadfvuns()
{
  //char *file_name = "pendulum.fvuns";
  //= "/home/lse24/lsestud/plorenz/starccm/FFS2D_V5_10_/FS/FFS2D_V5_10_0_exampleFS.fvuns";
  FILE *outfp;
  int i;
  int endian_conv = 0;
  float time;

  char cbuf[80];
  int ibuf[10];
  float fbuf[4];

  std::vector <int> bf,nb,apf;

  if ((outfp = fopen(file_name, "rb")) == NULL)
  {
    std::cout << "Cannot open input file\n";
    perror(NULL);
    exit(-1);
  }
  fread(ibuf, sizeof(int), 1, outfp); // FV_MAGIC
  int fv_magic = ibuf[0];
  if (fv_magic != FV_MAGIC)
  {
    endian_convert((char*) &fv_magic, sizeof(int));
    if (fv_magic != FV_MAGIC)
    {
      std::cout << "Wrong value for FV_MAGIC\n";
      std::cout << "FV_MAGIC=" << ibuf[0] << std::endl;
      exit(-1);
    }
    else
    {
      endian_conv = 1;
      if (VERBOSE)
        std::cout << "Using endian conversion.\n";
    }
  }

  READ(cbuf, 80); // 80 character string "FIELDVIEW"
  if (VERBOSE)
    std::cout << "FV     =" << cbuf << std::endl;

  READ(ibuf, 2); // Version, should be 3.0
  if (ibuf[0] != 3 || ibuf[1] != 0)
  {
    std::cout << "Wrong file version\n";
    std::cout << "Version=" << ibuf[0] << "." << ibuf[1] << std::endl;
    exit(-1);
  }

  READ(ibuf, 1); // File type
  int file_type = ibuf[0];
  if (file_type != FV_GRIDS_FILE && file_type != FV_RESULTS_FILE && file_type != FV_COMBINED_FILE)
  {
    std::cout << "Wrong file type\n";
    std::cout << "File type=" << ibuf[0] << std::endl;
    exit(-1);
  }

  READ(ibuf, 1); // Reserved field, should be zero
  if (ibuf[0] != 0)
  {
    std::cout << "Unexpected value for reserved field\n";
    std::cout << "Reserved field=" << ibuf[0] << std::endl;
    exit(-1);
  }

  if (file_type != FV_GRIDS_FILE)
  {
    READ(fbuf, 4); // constants for time, fsmach, alpha and re
    time = fbuf[0];
    if (VERBOSE)
    {
      std::cout << "time =" << fbuf[0] << std::endl;   
      std::cout << "fsmach =" << fbuf[1] << std::endl;   
      std::cout << "alpha =" << fbuf[2] << std::endl;           
      std::cout << "re =" << fbuf[3] << std::endl;
    }
  }

  READ(ibuf, 1); // number of grids
  /*int*/ num_grids=ibuf[0];
  if (VERBOSE)
    std::cout << "Number of meshes=" << num_grids << std::endl;

  /*int * */results_flag = 0;
  /*int * */normals_flag = 0;
  int num_face_types = 0;
  if (file_type != FV_RESULTS_FILE)
  {
    READ(ibuf, 1); // number of boundary ids
    num_face_types=ibuf[0];
    if (VERBOSE)
      std::cout << "Number of face types=" << num_face_types << std::endl;

    results_flag = new int[(int)num_face_types];   //num_face_types als int gecastet!!
    normals_flag = new int[(int)num_face_types];  //num_face_types als int gecastet!!
    assert(results_flag && normals_flag);
    
    for (i = 0; i < num_face_types; i++)
    {    
      READ(ibuf, 2); // results_flag and normals_flag

      if (VERBOSE)
      {
        std::cout << "results_flag=" << ibuf[0] << std::endl;
        std::cout << "normals_flag=" << ibuf[1] << std::endl;
      }

      results_flag[i] = ibuf[0];
      normals_flag[i] = ibuf[1];

      READ(cbuf, 80); // boundary name
      if (VERBOSE)
        std::cout << "Boundary name     =" << cbuf << std::endl;
    }
  }

  int num_vars = 0;
  int num_bvars = 0;
  if (file_type != FV_GRIDS_FILE)
  {
    READ(ibuf, 1); // number of variables
    num_vars=ibuf[0];
    if (VERBOSE)
      std::cout << "Number of variables=" << num_vars << std::endl;    

    for (i = 0; i < num_vars; i++)
    {
      READ(cbuf, 80); // variable name
      if (VERBOSE)
        std::cout << "Variable name =" << cbuf << std::endl;      
    }

    READ(ibuf, 1); // number of boundary variables
    num_bvars=ibuf[0];
    if (VERBOSE)
      std::cout << "Number of boundary variables=" << num_bvars << std::endl;    

    for (i = 0; i < num_bvars; i++)
    {
      READ(cbuf, 80); // boundary variable name
      if (VERBOSE)
        std::cout << "boundary variable names =" << cbuf << std::endl;      
    }
  }

  //read the grid data
  /*Part*/ parts = new Part[num_grids];
  assert(parts);
  for (int grid = 0; grid < num_grids; grid++)
  {
    READ(ibuf, 2); // read flag (should be FV_NODES) and number of nodes
    if (ibuf[0] != FV_NODES)
    {
      std::cerr << "Wrong grid flag\n";
      std::cerr << "Grid flag=" << ibuf[0] << std::endl;
      exit(-1);
    }
    int nnodes=ibuf[1];
    if (VERBOSE)
      std::cout << "Number of grid nodes=" << nnodes << std::endl;
    
    parts[grid].num_vars = num_vars;
    parts[grid].num_bvars = num_bvars;
    parts[grid].nnodes = nnodes;
    

    if (file_type != FV_RESULTS_FILE)
    {
      float *xyz= new float[nnodes*3];
      assert(xyz);

      READ(xyz, nnodes*3); // read x, y and z coordinates of all nodes
    
      parts[grid].xyz = xyz;
      //std::cout<<"x,y,z_ "<<parts[grid].xyz[0]<<" "<<parts[grid].xyz[1]<<" "<<parts[grid].xyz[2]<<std::endl;
      if (VERBOSE)
        for(i=0;i<nnodes;i++){
          
          if(i<5)
          std::cout << "vertices =" << i << " " << 
            xyz[i] << " " << xyz[i+nnodes] << " " << xyz[i+2*nnodes] << std::endl;
        }
    
      bool read_boundaryfaces_done = false;
      int num_bndpatches_tri_quad = 0;
      int num_bndpatches_poly = 0;
      while (!read_boundaryfaces_done)
      {
        int boundary_type;
        int num_faces;
        std::vector<int> ivec;
        int *ids;
        READ(ibuf, 1); // read boundary face flag 
        switch (ibuf[0])
        {
          case FV_FACES: // this is a boundary patch of tri/quad faces
            READ(ibuf, 2); // boundary type (plus 1) and number of faces of that type
            boundary_type = ibuf[0] - 1;
            num_faces = ibuf[1];
            num_bndpatches_tri_quad += 1;
            parts[grid].boundary_types.push_back(boundary_type);
            parts[grid].num_faces.push_back(num_faces);
            if (VERBOSE)
            {
              //if(i<5)
              std::cout << "boundary type =" << boundary_type << std::endl;
              std::cout << "number of tri/quad faces =" << num_faces << std::endl;
            }
            
            ids = new int[num_faces*4];
            assert(ids);
            
            READ(ids, num_faces*4);
            for (i = 0; i < num_faces*4; ++i)
              ids[i] -= 1; // convert to zero-based indexing

            if (VERBOSE)
              for (i = 0; i < num_faces; ++i)
              {  
                 //if(i<5)
                std::cout << "face= " << i << " boundary vert =" << 
                  ids[i*4] << " " << ids[i*4+1] << " " << ids[i*4+2] << " " << ids[i*4+3] << std::endl;
              }
            
            parts[grid].face_ids.push_back(ids);
            
            break;
            
          case FV_ARB_POLY_FACES: // this is a boundary patch of poly faces
            READ(ibuf, 2); // boundary type (plus 1) and number of faces of that type
            boundary_type = ibuf[0] - 1;
            num_faces = ibuf[1];
            num_bndpatches_poly += 1;
            parts[grid].boundary_types.push_back(boundary_type);
            parts[grid].num_faces.push_back(num_faces);
            if (VERBOSE)
            {
            //  if(i<5)
              std::cout << "boundary type =" << boundary_type << std::endl;
              std::cout << "number of poly faces =" << num_faces << std::endl;
            }
            
            // format of ids is:
            // for face=1,num_faces; int number of vertex ids; list of int vertex ids; endfor
            for (int face = 0; face < num_faces; ++face)
            {              
              READ(ibuf, 1);
              int nvert = ibuf[0];
              ivec.push_back(nvert);
              for (int ivert = 0; ivert < nvert; ++ivert)
              {
                READ(ibuf, 1);
                ivec.push_back(ibuf[0] - 1);
              }
            }
            
            ids = new int[ivec.size()];
            assert(ids);
            for (i = 0; i < (int) ivec.size(); ++i)
              ids[i] = ivec[i];
            
            if (VERBOSE)
            {
              int idx = 0;
              for (int face = 0; face < num_faces; ++face)
              {
                //if(i<5)        
                std::cout << "face= " << face << " boundary vert =";
                int nvert = ids[idx++];
                for (i = 0; i < nvert; ++i)
                  std::cout << ids[idx++] << " ";
                std::cout << std::endl;
              }
            }
            
            parts[grid].face_ids.push_back(ids);
            ivec.clear();
            
            break;
            
          case FV_ELEMENTS:
            read_boundaryfaces_done = true;
            break;
          
          default:
            std::cerr << "Unexpected tag\n";
            std::cerr << "tag=" << ibuf[0] << std::endl;
            exit(-1);
        }
      }
      parts[grid].num_bndpatches_tri_quad = num_bndpatches_tri_quad;
      parts[grid].num_bndpatches_poly = num_bndpatches_poly;

      READ(ibuf, 4);
      int ntet  = ibuf[0];
      int nhex  = ibuf[1];
      int npri  = ibuf[2];
      int npyr  = ibuf[3];
      parts[grid].ntet = ntet;
      parts[grid].nhex = nhex;
      parts[grid].npri = npri;
      parts[grid].npyr = npyr;
      if (VERBOSE)
      {
        std::cout << "n tet =" << ntet << std::endl;
        std::cout << "n hex =" << nhex << std::endl;
        std::cout << "n pri =" << npri << std::endl;
        std::cout << "n pyr =" << npyr << std::endl;
      }
      int ncells = ntet+nhex+npri+npyr;
      
      int *tetcells, *hexcells, *pricells, *pyrcells;
      smartalloc<int>(tetcells, ntet*4);
      smartalloc<int>(hexcells, nhex*8);
      smartalloc<int>(pricells, npri*6);
      smartalloc<int>(pyrcells, npyr*5);
      
      parts[grid].tet_ids = tetcells;
      parts[grid].hex_ids = hexcells;
      parts[grid].pri_ids = pricells;
      parts[grid].pyr_ids = pyrcells;
      
      ntet = 0;
      nhex = 0;
      npri = 0;
      npyr = 0;

      for (int icell = 0; icell < ncells; icell++)
      {   
        READ(ibuf, 1); // element header
        
        // determine element type from element header
        // Note: wall face info (Bits 1 to ELEM_TYPE_BIT_SHIFT) is discarded!!!
        unsigned int elem_header = (unsigned int) ibuf[0];
        int elem_type = elem_header >> ELEM_TYPE_BIT_SHIFT;
        
        switch (elem_type)
        {
          case 1: // tet cell
            READ(ibuf, 4); // vertex ids of this cell
            
            for (i = 0; i < 4; ++i)
              tetcells[4*ntet+i] = ibuf[i] - 1;
            
            if (VERBOSE)
            {
              std::cout << "cell " << icell << " vertex=";
              for (i=0;i<4;i++){
                std::cout << tetcells[4*ntet+i] << " ";
              }
              std::cout << std::endl;
            }
            
            ntet += 1;
            
            break;
            
          case 2: // pyramid cell
            READ(ibuf, 5); // vertex ids of this cell
            
            for (i = 0; i < 5; ++i)
              pyrcells[5*npyr+i] = ibuf[i] - 1;
            
            if (VERBOSE)
            {
              std::cout << "cell " << icell << " vertex=";
              for (i=0;i<5;i++){
                std::cout << pyrcells[5*npyr+i] << " ";
              }
              std::cout << std::endl;
            }
            
            npyr += 1;
            
            break;
            
          case 3: // prism cell
            READ(ibuf, 6); // vertex ids of this cell
            
            for (i = 0; i < 6; ++i)
              pricells[6*npri+i] = ibuf[i] - 1;
            
            if (VERBOSE)
            {
              std::cout << "cell " << icell << " vertex=";
              for (i=0;i<6;i++){
                std::cout << pricells[6*npri+i] << " ";
              }
              std::cout << std::endl;
            }
            
            npri += 1;
            
            break;
            
          case 4: // hex cell
            READ(ibuf, 8); // vertex ids of this cell
            
            for (i = 0; i < 8; ++i)
              hexcells[8*nhex+i] = ibuf[i] - 1;
            
            if (VERBOSE)
            {
              //if(icell<5)
              std::cout << "cell " << icell << " vertex=";
              for (i=0;i<8;i++){
              
                //              if(icell<5)
                std::cout << hexcells[8*nhex+i] << " ";
              }
              //if(icell<5)
              std::cout << std::endl;
            }
            
            nhex += 1;
            
            break;
            
          default:
            std::cerr << "Unknown element type\n";
            std::cerr << "type=" << elem_type << std::endl;
            exit(-1);
        }
      }
      
      // Note: no support for polyhedral cells at the moment
    }

    if (file_type != FV_GRIDS_FILE)
    {
      READ(ibuf, 1); // read next tag
      if (ibuf[0] == FV_ARB_POLY_ELEMENTS)
      {
        std::cerr << "No support for polyhedral cells\n";
        exit(-1);
      }

      if (ibuf[0] != FV_VARIABLES)
      {
        std::cerr << "Unexpected tag\n";
        std::cerr << "tag=" << ibuf[0] << std::endl;
        exit(-1);
      }
      
      float *vars= new float[num_vars*nnodes*sizeof(float)];

      READ(vars, num_vars*nnodes);
      
      parts[grid].vars = vars;

      if (VERBOSE)
        for (int ivar = 0; ivar < num_vars; ++ivar)
        {
          std::cout << "Variable " << ivar << std::endl;
          for(i = 0; i < nnodes; i++)
          {
          
                  //        if(i<5)
            std::cout << i << " = " << vars[i+ivar*nnodes] << std::endl;
          }
        }

      READ(ibuf, 1); // FV_BNDRY_VARS
      if (ibuf[0] != FV_BNDRY_VARS)
      {
        std::cerr << "Unexpected tag\n";
        std::cerr << "tag=" << ibuf[0] << std::endl;
        exit(-1);
      }
      
      // get number of all boundary tri/quad faces
      int num_faces_total = 0;
      for (int ipatch = 0; ipatch < parts[grid].num_bndpatches_tri_quad; 
           ++ipatch)
        num_faces_total += parts[grid].num_faces[ipatch];
      
      // read all variables on all boundary tri/quad faces
      if (num_bvars && num_faces_total)
      {
        float *bvars= new float[num_bvars*num_faces_total*sizeof(float)];

        READ(bvars, num_bvars*num_faces_total);

        parts[grid].boundary_vars.push_back(bvars);

        if (VERBOSE)
          for (int ivar = 0; ivar < num_bvars; ++ivar)
          {
            std::cout << "Boundary Variable " << ivar << std::endl;
            int patchstart = 0;
            for (int ipatch = 0; ipatch < parts[grid].num_bndpatches_tri_quad; 
                 ++ipatch)
            {
              int num_faces = parts[grid].num_faces[ipatch];

              std::cout << "Patch " << ipatch << std::endl;
              for(i = 0; i < num_faces; i++)
              {
              //                if(i<5)
                std::cout << i << " = " 
                          << bvars[i+patchstart+ivar*num_faces_total] << std::endl;
              }patchstart += num_faces;
            }
          }
      }

      READ(ibuf, 1); // FV_ARB_POLY_BNDRY_VARS
      if (ibuf[0] != FV_ARB_POLY_BNDRY_VARS)
      {
        std::cout << "Unexpected tag\n";
        std::cout << "tag=" << ibuf[0] << std::endl;
        exit(-1);
      }

      // get number of all boundary poly faces
      num_faces_total = 0;
      for (unsigned int ipatch = parts[grid].num_bndpatches_tri_quad; 
           ipatch < parts[grid].num_faces.size(); ++ipatch)
        num_faces_total += parts[grid].num_faces[ipatch];
      
      // read all variables on all boundary poly faces
      if (num_bvars && num_faces_total)
      {
        float *bvars= new float[num_bvars*num_faces_total*sizeof(float)];

        READ(bvars, num_bvars*num_faces_total);

        parts[grid].boundary_vars.push_back(bvars);

        if (VERBOSE)
          for (int ivar = 0; ivar < num_bvars; ++ivar)
          {
            std::cout << "Boundary Variable " << ivar << std::endl;
            int patchstart = 0;
            for (unsigned int ipatch = parts[grid].num_bndpatches_tri_quad; 
                 ipatch < parts[grid].num_faces.size(); ++ipatch)
            {
              int num_faces = parts[grid].num_faces[ipatch];

              std::cout << "Patch " << ipatch << std::endl;
              for(i = 0; i < num_faces; i++)
              {
              
//                              if(i<5)
                std::cout << i << " = " 
                          << bvars[i+patchstart+ivar*num_faces_total] << std::endl;
              }patchstart += num_faces;
              
            }
          }
      }
    }
  }
  
  fclose(outfp);
} // void Initreadfvuns()

void cleanup(){
  // cleanup
  delete [] parts;
  
  if (results_flag)
    delete [] results_flag;
  
  if (normals_flag)
    delete [] normals_flag;
  
}
  



