/***************************************************************************
    File        : MeshReader.cpp
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Tue Jan 22 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#include "MeshReader.hh"


namespace grd {


//===========================================================================
//
//  Unit Square - a square with coords in [0,1]x[0,1]
//
//===========================================================================


void
MeshReader::createTriangleMesh()
{
  double x[8][3] = {{0.,0.,0.},{0.5,0.,0.},
                    {0.,0.5,0.},{0.5,0.5,0.},{1.,0.5,0.},
                    {0.,1.,0.},{0.5,1.,0.},{1.,1.,0.}};

  Vertex* v0 = new Vertex(x[0]);
  Vertex* v1 = new Vertex(x[1]);
  Vertex* v2 = new Vertex(x[2]);
  Vertex* v3 = new Vertex(x[3]);
  Vertex* v4 = new Vertex(x[4]);
  Vertex* v5 = new Vertex(x[5]);
  Vertex* v6 = new Vertex(x[6]);
  Vertex* v7 = new Vertex(x[7]);

  // Set Ids
  v0->setId(1);
  v1->setId(2);
  v2->setId(3);
  v3->setId(4);
  v4->setId(5);
  v5->setId(6);
  v6->setId(7);
  v7->setId(8);

  // Create Triangles
  Triangle* tri0 = new Triangle;
  Triangle* tri1 = new Triangle;
  Triangle* tri2 = new Triangle;
  Triangle* tri3 = new Triangle;
  Triangle* tri4 = new Triangle;
  Triangle* tri5 = new Triangle;

  // Set vertices
  tri0->setVertex(0,v0);
  tri0->setVertex(1,v1);
  tri0->setVertex(2,v2);
  
  tri1->setVertex(0,v2);
  tri1->setVertex(1,v1);
  tri1->setVertex(2,v3);

  tri2->setVertex(0,v2);
  tri2->setVertex(1,v3);
  tri2->setVertex(2,v5);

  tri3->setVertex(0,v5);
  tri3->setVertex(1,v3);
  tri3->setVertex(2,v6);

  tri4->setVertex(0,v3);
  tri4->setVertex(1,v4);
  tri4->setVertex(2,v6);

  tri5->setVertex(0,v6);
  tri5->setVertex(1,v4);
  tri5->setVertex(2,v7);

  // Set lists
  // Set vertex list
  vertex.resize(8);

  vertex[0] = v0;
  vertex[1] = v1;
  vertex[2] = v2;
  vertex[3] = v3;
  vertex[4] = v4;
  vertex[5] = v5;
  vertex[6] = v6;
  vertex[7] = v7;
  
  element.resize(6);
  element[0] = tri0;
  element[1] = tri1;
  element[2] = tri2;
  element[3] = tri3;
  element[4] = tri4;
  element[5] = tri5;
  
}

void
MeshReader::createUnitSquareWithTris()
{
  double x[8][3] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0}};

  // Create the eight vertices
  Vertex* v0 = new Vertex(x[0]);
  Vertex* v1 = new Vertex(x[1]);
  Vertex* v2 = new Vertex(x[2]);
  Vertex* v3 = new Vertex(x[3]);

  // Set vertex ids
  v0->setId(1);
  v1->setId(2);
  v2->setId(3);
  v3->setId(4);
  
  // Set vertex list
  vertex.resize(4);

  vertex[0] = v0;
  vertex[1] = v1;
  vertex[2] = v2;
  vertex[3] = v3;
  
  // create element
  Triangle* tri0 = new Triangle;
  Triangle* tri1 = new Triangle;
  
  // Set vertices
  tri0->setVertex(0,v0);
  tri0->setVertex(1,v1);
  tri0->setVertex(2,v2);
  
  tri1->setVertex(0,v0);
  tri1->setVertex(1,v2);
  tri1->setVertex(2,v3);
  
  // Set element list
  element.resize(2);
  element[0] = tri0;
  element[1] = tri1;
}


void 
MeshReader::createUnitSquareWithQuad()
{
    double x[8][3] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0}};

  // Create the eight vertices
  Vertex* v0 = new Vertex(x[0]);
  Vertex* v1 = new Vertex(x[1]);
  Vertex* v2 = new Vertex(x[2]);
  Vertex* v3 = new Vertex(x[3]);

  // Set vertex ids
  v0->setId(1);
  v1->setId(2);
  v2->setId(3);
  v3->setId(4);
  
  // Set vertex list
  vertex.resize(4);

  vertex[0] = v0;
  vertex[1] = v1;
  vertex[2] = v2;
  vertex[3] = v3;
  
  // create element
  Quadrangle* quad = new Quadrangle;
  
  // Set vertices
  quad->setVertex(0,v0);
  quad->setVertex(1,v1);
  quad->setVertex(2,v2);
  quad->setVertex(3,v3);
  
  // Set element list
  element.resize(1);
  element[0] = quad;  
}


//===========================================================================
//
//  Unit Cube - a cube with coord in [0,1]x[0,1]x[0,1]
//
//===========================================================================

//
//
void
MeshReader::createUnitCubeWithTets()
{
  double x[8][3] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                    {0,0,1},{1,0,1},{1,1,1},{0,1,1}};

  // Create the eight vertices
  Vertex* v0 = new Vertex(x[0]);
  Vertex* v1 = new Vertex(x[1]);
  Vertex* v2 = new Vertex(x[2]);
  Vertex* v3 = new Vertex(x[3]);
  Vertex* v4 = new Vertex(x[4]);
  Vertex* v5 = new Vertex(x[5]);
  Vertex* v6 = new Vertex(x[6]);
  Vertex* v7 = new Vertex(x[7]);

  // Set vertex id
  v0->setId(1);
  v1->setId(2);
  v2->setId(3);
  v3->setId(4);
  v4->setId(5);
  v5->setId(6);
  v6->setId(7);
  v7->setId(8);

  // Set vertex list
  vertex.resize(8);

  vertex[0] = v0;
  vertex[1] = v1;
  vertex[2] = v2;
  vertex[3] = v3;
  vertex[4] = v4;
  vertex[5] = v5;
  vertex[6] = v6;
  vertex[7] = v7;

  // Create element list
  Tetrahedron* tet0 = new Tetrahedron;
  Tetrahedron* tet1 = new Tetrahedron;
  Tetrahedron* tet2 = new Tetrahedron;
  Tetrahedron* tet3 = new Tetrahedron;
  Tetrahedron* tet4 = new Tetrahedron;

  // Set tets vertices
  tet0->setVertex(0,v0);
  tet0->setVertex(1,v7);
  tet0->setVertex(2,v4);
  tet0->setVertex(3,v5);

  tet1->setVertex(0,v0);
  tet1->setVertex(1,v1);
  tet1->setVertex(2,v2);
  tet1->setVertex(3,v5);

  tet2->setVertex(0,v2);
  tet2->setVertex(1,v6);
  tet2->setVertex(2,v7);
  tet2->setVertex(3,v5);

  tet3->setVertex(0,v0);
  tet3->setVertex(1,v2);
  tet3->setVertex(2,v3);
  tet3->setVertex(3,v7);

  tet4->setVertex(0,v0);
  tet4->setVertex(1,v2);
  tet4->setVertex(2,v7);
  tet4->setVertex(3,v5);


  // Set element list
  element.resize(5);
  element[0] = tet0;
  element[1] = tet1;
  element[2] = tet2;
  element[3] = tet3;
  element[4] = tet4;


}

// Description
//  Create a unit cube and
//  generate a Hex element
void
MeshReader::createUnitCubeWithHexs()
{
  double x[8][3] = {{0,0,0},{1,0,0},{1,1,0},{0,1,0},
                    {0,0,1},{1,0,1},{1,1,1},{0,1,1}};

  // Create the eight vertices
  Vertex* v0 = new Vertex(x[0]);
  Vertex* v1 = new Vertex(x[1]);
  Vertex* v2 = new Vertex(x[2]);
  Vertex* v3 = new Vertex(x[3]);
  Vertex* v4 = new Vertex(x[4]);
  Vertex* v5 = new Vertex(x[5]);
  Vertex* v6 = new Vertex(x[6]);
  Vertex* v7 = new Vertex(x[7]);

  // Set vertex id
  v0->setId(1);
  v1->setId(2);
  v2->setId(3);
  v3->setId(4);
  v4->setId(5);
  v5->setId(6);
  v6->setId(7);
  v7->setId(8);

  // Set vertex list
  vertex.resize(8);

  vertex[0] = v0;
  vertex[1] = v1;
  vertex[2] = v2;
  vertex[3] = v3;
  vertex[4] = v4;
  vertex[5] = v5;
  vertex[6] = v6;
  vertex[7] = v7;

  // Create element list
  Hexahedron* hex = new Hexahedron;
  element.resize(1);
  element[0] = hex;
  // Set vertices
  hex->setVertex(0,v0);
  hex->setVertex(1,v1);
  hex->setVertex(2,v2);
  hex->setVertex(3,v3);
  hex->setVertex(4,v4);
  hex->setVertex(5,v5);
  hex->setVertex(6,v6);
  hex->setVertex(7,v7);
}



//===========================================================================
//
//  CAPA - unstructured mixed meshes
//
//===========================================================================

//
//
void
MeshReader::readCAPA(std::string& filename)
{
  std::ifstream theFile;
  theFile.open(filename.c_str());

  std::vector<int>     vertexId;
  std::vector<double*> vertexValues;
  std::vector<int>  elemId;
  std::vector<int*> elemTopology;

  // Check if file could be opend
  if (theFile.fail())
  {
    Error("can't open input file: ",filename.c_str());
  }

  while (!theFile.eof())
  {
    char buffer[1024];
    theFile.getline(buffer,1024);
    std::istringstream cline(buffer);
    std::string cword;
    std::string minusEins("-1");
    std::string key781("781");
    std::string key780("780");
    cline >> cword;
    if (cword == minusEins)
      Info("Found ",cword);
    else if (cword == key781)  {
      Info("Found ",cword);
      readCAPAVertex(theFile,vertexId,vertexValues);
    }
    else if (cword == key780) {
      Info("Found ",cword);
      readCAPAElement(theFile,elemId,elemTopology);
    }
  }

  setCAPAVertexAndElement(vertexId,vertexValues,elemId,elemTopology);
}



// Description
// Read CAPA vertices
void
MeshReader::readCAPAVertex(std::ifstream& theFile,std::vector<int>& vertexId,std::vector<double*>& vertexValues)
{
  while (!theFile.eof()) {
    char buffer[1024];
    theFile.getline(buffer,1024);
    std::istringstream cline(buffer);
    std::string cword;
    std::string minusEins("-1");
    cline >> cword;
    if (cword == minusEins)
      return;
    // vertex id
    int vrtId = atoi(cword.c_str());

    // vertex coordinates
    double* vtVal = new double[3];
    theFile.getline(buffer,1024);
    std::istringstream vline(buffer);
    vline >> cword;
    vtVal[0] = atof(cword.c_str());
    vline >> cword;
    vtVal[1] = atof(cword.c_str());
    vline >> cword;
    vtVal[2] = atof(cword.c_str());
    vertexValues.push_back(vtVal);
    vertexId.push_back(vrtId);
  }
}

// Description
// Read CAPA Elements
void
MeshReader::readCAPAElement(std::ifstream& theFile,std::vector<int>& elemId,std::vector<int*>& elemTopology)
{
  int i;

  while (!theFile.eof()) {
    // read buffer
    char buffer[1024];
    // read element description line
    theFile.getline(buffer,1024);
    std::istringstream dline(buffer);
    std::string cword;
    std::string minusEins("-1");
    dline >> cword;
    if (cword == minusEins)
      return;
    int id = atoi(cword.c_str());
    int type;
    for (i = 0; i < 7; i++)
    {
       dline >> cword;
    }
    type = atoi(cword.c_str());

    // read element
    theFile.getline(buffer,1024);
    std::istringstream vline(buffer);

    int* vtr = new int[type+1];
    vtr[0] = type;
    for (i = 0; i < type; i++)
    {
      vline >> cword;
      vtr[i+1] = atoi(cword.c_str());
    }

    elemTopology.push_back(vtr);
    elemId.push_back(id);
  }
}


//
//
void
MeshReader::setCAPAVertexAndElement(std::vector<int>& vertexId,std::vector<double*>& vertexValues,
                                    std::vector<int>& elemId,std::vector<int*>& elem)
{
  size_t i;
  size_t noOfVertices = vertexId.size();
  size_t noOfElements = elemId.size();

  // Create Vertices
  vertex.resize(noOfVertices);
  for (i = 0; i < noOfVertices; i++)
  {
    Vertex* tmpVt = new Vertex(vertexValues[i]);
    tmpVt->setId(vertexId[i]);
    vertex[vertexId[i]-1] = tmpVt;
  }

  // Create Elements
  element.resize(noOfElements);
  for (i = 0; i < noOfElements; i++)
  {
    int type = elem[i][0];
    if (type == 3)
    {
      Triangle* tmpElem = new Triangle;
      tmpElem->setVertex(0,vertex[elem[i][1]-1]);
      tmpElem->setVertex(1,vertex[elem[i][2]-1]);
      tmpElem->setVertex(2,vertex[elem[i][3]-1]);
      element[elemId[i]-1] = tmpElem;
    }
    else if (type == 4)
    {
      Quadrangle* tmpElem = new Quadrangle;
      tmpElem->setVertex(0,vertex[elem[i][1]-1]);
      tmpElem->setVertex(1,vertex[elem[i][2]-1]);
      tmpElem->setVertex(2,vertex[elem[i][3]-1]);
      tmpElem->setVertex(3,vertex[elem[i][4]-1]);
      element[elemId[i]-1] = tmpElem;
    }
  }

}







//===========================================================================
//
//  GRUMMP - Tetrahedral meshes
//
//===========================================================================

// Description
void
MeshReader::readGRUMMP(const std::string& filename)
{
  int i,len;
  std::ifstream theFile;
  theFile.open(filename.c_str());
  std::vector<int*> elemTopology;
  std::vector<double*> vertexValues;

  // Check if file could be opend
  if (theFile.fail()) {
    Error("can't open input file",filename.c_str());
  }

  while (!theFile.eof()) {
    char buffer[1024];
    theFile.getline(buffer,1024);

    // Read the Element topology
    if (strstr(buffer,"volumeelements")) {
      theFile.getline(buffer,1024);
      int noOfElements = atoi(buffer);
      elemTopology.resize(noOfElements);
      for (i = 0; i < noOfElements; i++)
      {
        int* vtr = new int[4];
        len = 0;
        while (len <= 0) {
          theFile.getline(buffer,1024);
          len = strlen(buffer);
          if (len > 0) {
            char* sep = " ";
            char* ptr = buffer;
            char* theToken;
            int counter = 0;

            while ((theToken = strtok(ptr,sep)) != NULL) {
              ptr = NULL;
              if (counter > 1)
                vtr[counter-2] = atoi(theToken);
              counter++;
              if (counter > 5)
                break;
            } // while()
          } // if (len > 0)
        } // while (len )
        elemTopology[i] = vtr;
      } // for noOfElements
    } // if (volumenelements)

    // Read the Vertices coordinates
    if (strstr(buffer,"points")) {
      theFile.getline(buffer,1024);
      int noOfVertices = atoi(buffer);
      vertexValues.resize(noOfVertices);
      for (i = 0; i < noOfVertices; i++)
      {
        double* vtr = new double[3];
        len = 0;
        while(len <= 0){
          theFile.getline(buffer,1024);
          len = strlen(buffer);
          if (len > 0) {
            char* sep = " ";
            char* ptr = buffer;
            char* theToken;
            int counter = 0;
            while ((theToken = strtok(ptr,sep)) != NULL) {
              ptr = NULL;
              vtr[counter] = (double)atof(theToken);
              counter++;
              if (counter > 2 )
                break;
            } // switch()
          } // if (len > 4)
        } // while(theFile.eof())
        vertexValues[i] = vtr;
      } // for noOfVertices
    }
  }

  setVertexAndElement(elemTopology,vertexValues);
}



// Description
// Read LSE-Mesh data
void 
MeshReader::readLSEMesh(const std::string& filename)
{
  int noNodes;
  int no3DElems;
  int i,len;
  std::ifstream theFile;
  theFile.open(filename.c_str());
  std::vector<int*>    elemTopology;
  std::vector<double*> vertexValues;

  // Check if file could be opend
  if (theFile.fail()) {
    Error("can't open input file",filename.c_str());
  }

  while (!theFile.eof()) {
    char buffer[1024];
    buffer[0] = 0;
    theFile.getline(buffer,1024);

    // Check if comment line
    if (buffer[0] == '#')
      continue;
    
    // Read the Element topology
    // Variables 
    //   NumNodes
    //   Num3DElements
    //   Num2DElements
    if (strstr(buffer,"NumNodes")) 
    {
      char* sep = " ";
      char* ptr = buffer;
      char* theToken;
      
      // call this function twice      
      theToken = strtok(ptr,sep);
      ptr = NULL;
      theToken = strtok(ptr,sep);
      noNodes = atoi(theToken);
    }
    else if (strstr(buffer,"Num3DElements")) 
    {
      char* sep = " ";
      char* ptr = buffer;
      char* theToken;
      
      // call this function twice      
      theToken = strtok(ptr,sep);
      ptr = NULL;
      theToken = strtok(ptr,sep);
      no3DElems = atoi(theToken);
    }

    // Read element topology
    else if (strstr(buffer,"[3D Elements]"))
    {
      elemTopology.resize(no3DElems);
      for (i = 0; i < no3DElems; i++)
      {
        int* vtr = new int[4];
        len = 0;
        while (len <= 0) 
        {
          buffer[0] = 0;  
          theFile.getline(buffer,1024);
          while (buffer[0] == '#')
          {
            buffer[0] = 0;  
            theFile.getline(buffer,1024);
          }
          // read another line
          theFile.getline(buffer,1024);
          len = strlen(buffer);
          if (len > 0)
          {
            char* sep = " ";
            char* ptr = buffer;
            char* theToken;

            // call this function four times
            theToken = strtok(ptr,sep);
            ptr = NULL;
            vtr[0] = atoi(theToken);
            // the first coordinate
            theToken = strtok(ptr,sep);
            ptr = NULL;
            vtr[1] = atoi(theToken);
            // the second coordinate
            theToken = strtok(ptr,sep);
            ptr = NULL;
            vtr[2] = atoi(theToken);
            // the third coordinate
            theToken = strtok(ptr,sep);
            ptr = NULL;
            vtr[3] = atoi(theToken);
          } // if (len > 0)
        } // while (len )
        elemTopology[i] = vtr;
      } // for noOfElements
    } // if (volumenelements)

    // Read the Vertices coordinates
    else if (strstr(buffer,"[Nodes]")) 
    {
      vertexValues.resize(noNodes);
      for (i = 0; i < noNodes; i++)
      {
        double* vtr = new double[3];
        buffer[0] = 0;  
        theFile.getline(buffer,1024);
        while (buffer[0] == '#')
        {
          buffer[0] = 0;  
          theFile.getline(buffer,1024);
        }
        len = strlen(buffer);
        if (len > 0) 
        {
          char* sep = " ";
          char* ptr = buffer;
          char* theToken;
          // call this function four times
          theToken = strtok(ptr,sep);
          ptr = NULL;
          // the first coordinate
          theToken = strtok(ptr,sep);
          ptr = NULL;
          vtr[0] = (double)atof(theToken);
          // the second coordinate
          theToken = strtok(ptr,sep);
          ptr = NULL;
          vtr[1] = (double)atof(theToken);
          // the third coordinate
          theToken = strtok(ptr,sep);
          ptr = NULL;
          vtr[2] = (double)atof(theToken);
        } // if (len > 0)
        else
          std::cerr << "ERROR: wrong file format in LSE-mesh" << std::endl;
        vertexValues[i] = vtr;
      } // for noOfVertices
    }
  }

  setVertexAndElement(elemTopology,vertexValues);  
}

//
//
void
MeshReader::setVertexAndElement(std::vector<int*>& elem,std::vector<double*>& vtr)
{
  size_t i;
  size_t noOfVertices = vtr.size();
  size_t noOfElements = elem.size();

  // Create Vertices
  vertex.resize(noOfVertices);
  for (i = 0; i < noOfVertices; i++)
  {
    Vertex* tmpVt = new Vertex(vtr[i]);
    tmpVt->setId(i+1);
    vertex[i] = tmpVt;
  }

  // Create Elements
  element.resize(noOfElements);
  for (i = 0; i < noOfElements; i++)
  {
    Tetrahedron* tmpElem = new Tetrahedron;
    tmpElem->setVertex(0,vertex[elem[i][0]-1]);
    tmpElem->setVertex(1,vertex[elem[i][1]-1]);
    tmpElem->setVertex(2,vertex[elem[i][2]-1]);
    tmpElem->setVertex(3,vertex[elem[i][3]-1]);
    element[i] = tmpElem;
  }

}


} // namespace
