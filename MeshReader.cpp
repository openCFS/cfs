/***************************************************************************
    File        : MeshReader.cpp
    Description : 

 ---------------------------------------------------------------------------
    Begin       : Tue Jan 22 2002
    Author(s)   : Roberto Grosso
 ***************************************************************************/

#include "MeshReader.h"

namespace grd {

//===========================================================================
//
//  SRG - Hierarchical Tri-Mesehs
//
//===========================================================================


//
//
void
MeshReader::readSRG(string& filename)
{
  ifstream theFile;
  theFile.open(filename.c_str());

  vector<int>     vertexId;
  vector<double*> vertexValues;
  vector<int>  elemId;
  vector<int*> elemTopology;

  // Check if file could be opend
  if (theFile.fail()) {
    cerr << "ERROR: can't open input file: " << filename.c_str() << '\n';
    exit(0);
  }

  // read head
  char buffer[1024];
  theFile.getline(buffer,1024);
  theFile.getline(buffer,1024);
  int noOfLevel = atoi(buffer);

  // read Level 0
  readCoarseMeshSRG(theFile);
  // read next leves
  mlVertex.resize(noOfLevel-1);
  mlElement.resize(noOfLevel-1);
  for (int i = 0; i < noOfLevel-1; i++)
  {
    readMeshLevelSRG(theFile,i);
  }

  // build multilevel boundary
  buildMultilevelBoundary();
}


//
//
void
MeshReader::readCoarseMeshSRG(ifstream& theFile)
{
  char buffer[1024];
  theFile.getline(buffer,1024);
  istringstream headLine(buffer);
  string noOfVertexString;
  string noOfTrisString;
  headLine >> noOfVertexString;
  headLine >> noOfTrisString;

  int noOfVertex = atoi(noOfVertexString.c_str());
  int noOfTris   = atoi(noOfTrisString.c_str());
  // Read vertices
  vertex.resize(noOfVertex);
  for (int vt = 0; vt < noOfVertex; vt++)
  {
    theFile.getline(buffer,1024);
    istringstream cline(buffer);
    string cword;
    // get vertex position
    double pos[3];
    cline >> cword;
    pos[0] = atof(cword.c_str());
    cline >> cword;
    pos[1] = atof(cword.c_str());
    cline >> cword;
    pos[2] = atof(cword.c_str());
    Vertex* tmpVt = new Vertex(pos);
    tmpVt->setId(vt+1);

    // set vertex into list
    vertex[vt] = tmpVt;
  }

  // Read elements
  element.resize(noOfTris);
  for (int el = 0; el < noOfTris; el++)
  {
    char buffer[1024];
    theFile.getline(buffer,1024);
    istringstream cline(buffer);
    string cword;
    // get vertex position
    int vt[3];
    cline >> cword; // here is a 3
    cline >> cword;
    vt[0] = atoi(cword.c_str());
    cline >> cword;
    vt[1] = atoi(cword.c_str());
    cline >> cword;
    vt[2] = atoi(cword.c_str());
    Triangle* tmpTri = new Triangle;
    tmpTri->setVertex(0,vertex[vt[0]]);
    tmpTri->setVertex(1,vertex[vt[1]]);
    tmpTri->setVertex(2,vertex[vt[2]]);

    // set element into list
    element[el] = tmpTri;
  }

}



//
//
void
MeshReader::readMeshLevelSRG(ifstream& theFile,int level)
{
  char buffer[1024];
  theFile.getline(buffer,1024);
  istringstream headLine(buffer);
  string noOfVertexString;
  string noOfTrisString;
  headLine >> noOfVertexString;
  headLine >> noOfTrisString;

  int noOfVertex = atoi(noOfVertexString.c_str());
  int noOfTris   = atoi(noOfTrisString.c_str());

  // Read vertices
  double** pos = new double*[noOfVertex];
  for (int vt = 0; vt < noOfVertex; vt++)
  {
    theFile.getline(buffer,1024);
    istringstream cline(buffer);
    string cword;
    // get vertex position
    pos[vt] = new double[3];
    cline >> cword;
    pos[vt][0] = atof(cword.c_str());
    cline >> cword;
    pos[vt][1] = atof(cword.c_str());
    cline >> cword;
    pos[vt][2] = atof(cword.c_str());
  }

  // set stack
  mlVertex[level] = pos;

  // Read elements
  int** vtId = new int*[noOfTris];
  for (int el = 0; el < noOfTris; el++)
  {
    char buffer[1024];
    theFile.getline(buffer,1024);
    istringstream cline(buffer);
    string cword;
    // get vertex position
    vtId[el] = new int[3];
    cline >> cword; // here is a 3
    cline >> cword;
    vtId[el][0] = atoi(cword.c_str());
    cline >> cword;
    vtId[el][1] = atoi(cword.c_str());
    cline >> cword;
    vtId[el][2] = atoi(cword.c_str());
  }
  mlElement[level] = vtId;
}






//===========================================================================
//
//  CAPA - unstructured mixed meshes
//
//===========================================================================

//
//
void
MeshReader::readCAPA(string& filename)
{
  ifstream theFile;
  theFile.open(filename.c_str());

  vector<int>     vertexId;
  vector<double*> vertexValues;
  vector<int>  elemId;
  vector<int*> elemTopology;

  // Check if file could be opend
  if (theFile.fail()) {
    cerr << "ERROR: can't open input file: " << filename.c_str() << '\n';
    exit(0);
  }

  while (!theFile.eof()) {
    char buffer[1024];
    theFile.getline(buffer,1024);
    istringstream cline(buffer);
    string cword;
    string minusEins("-1");
    string key781("781");
    string key780("780");
    cline >> cword;
    if (cword == minusEins)
      cerr << "Found " << cword << '\n';
    else if (cword == key781)  {
      cerr << "Found " << cword << '\n';
      readCAPAVertex(theFile,vertexId,vertexValues);
    }
    else if (cword == key780) {
      cerr << "Found " << cword << '\n';
      readCAPAElement(theFile,elemId,elemTopology);
    }
  }

  setCAPAVertexAndElement(vertexId,vertexValues,elemId,elemTopology);
}



// Description
// Read CAPA vertices
void
MeshReader::readCAPAVertex(ifstream& theFile,vector<int>& vertexId,vector<double*>& vertexValues)
{
  while (!theFile.eof()) {
    char buffer[1024];
    theFile.getline(buffer,1024);
    istringstream cline(buffer);
    string cword;
    string minusEins("-1");
    cline >> cword;
    if (cword == minusEins)
      return;
    // vertex id
    int vrtId = atoi(cword.c_str());

    // vertex coordinates
    double* vtVal = new double[3];
    theFile.getline(buffer,1024);
    istringstream vline(buffer);
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
MeshReader::readCAPAElement(ifstream& theFile,vector<int>& elemId,vector<int*>& elemTopology)
{
  int i;

  while (!theFile.eof()) {
    // read buffer
    char buffer[1024];
    // read element description line
    theFile.getline(buffer,1024);
    istringstream dline(buffer);
    string cword;
    string minusEins("-1");
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
    istringstream vline(buffer);

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
MeshReader::setCAPAVertexAndElement(vector<int>& vertexId,vector<double*>& vertexValues,
                                    vector<int>& elemId,vector<int*>& elem)
{
  int i;
  int noOfVertices = vertexId.size();
  int noOfElements = elemId.size();

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
MeshReader::readGRUMMP(string& filename)
{
  int i,len;
  ifstream theFile;
  theFile.open(filename.c_str());
  vector<int*> elemTopology;
  vector<double*> vertexValues;

  // Check if file could be opend
  if (theFile.fail()) {
    cerr << "ERROR: can't open input file\n";
    exit(0);
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

//
//
void
MeshReader::setVertexAndElement(vector<int*>& elem,vector<double*>& vtr)
{
  int i;
  int noOfVertices = vtr.size();
  int noOfElements = elem.size();

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




//
//
//
void
MeshReader::buildMultilevelBoundary()
{
  int counter = 0;
  Element* tri;
  typedef vector<Element*>::iterator TriI;
  double** pos = mlVertex[0];
  int**    elm = mlElement[0];

  int noOfLevel = mlVertex.size();
  list<TriSkeleton*> triSk1;
  list<TriSkeleton*> triSk2;

  for (TriI p = element.begin(); p != element.end(); ++p)
  {
    tri = *p;
    TriSkeleton* skl = new TriSkeleton;
    TriSkeleton* ch0 = new TriSkeleton;
    TriSkeleton* ch1 = new TriSkeleton;
    TriSkeleton* ch2 = new TriSkeleton;
    TriSkeleton* ch3 = new TriSkeleton;

    for (int n = 0; n < 3; n++)
    {
      ch0->pos[0][n] = pos[elm[counter + 0][0]][n];
      ch0->pos[1][n] = pos[elm[counter + 0][1]][n];
      ch0->pos[2][n] = pos[elm[counter + 0][2]][n];

      ch1->pos[0][n] = pos[elm[counter + 1][0]][n];
      ch1->pos[1][n] = pos[elm[counter + 1][1]][n];
      ch1->pos[2][n] = pos[elm[counter + 1][2]][n];

      ch2->pos[0][n] = pos[elm[counter + 2][0]][n];
      ch2->pos[1][n] = pos[elm[counter + 2][1]][n];
      ch2->pos[2][n] = pos[elm[counter + 2][2]][n];

      ch3->pos[0][n] = pos[elm[counter + 3][0]][n];
      ch3->pos[1][n] = pos[elm[counter + 3][1]][n];
      ch3->pos[2][n] = pos[elm[counter + 3][2]][n];
    }

    // set position of root TriSkeleton
    for (int i = 0; i < 3; i++)
    {
      Vertex* v = tri->getVertex(i);
      skl->pos[i][0] = (*v)[0];
      skl->pos[i][1] = (*v)[1];
      skl->pos[i][2] = (*v)[2];
    }
    // set children
    skl->child[0] = ch0;
    skl->child[1] = ch1;
    skl->child[2] = ch2;
    skl->child[3] = ch3;

    ((Triangle*)tri)->setSkeleton(skl);
    counter += 4;

    triSk1.push_back(ch0);
    triSk1.push_back(ch1);
    triSk1.push_back(ch2);
    triSk1.push_back(ch3);
  }

  list<TriSkeleton*>::iterator q;
  for (int level = 1; level < noOfLevel; level++)
  {
    cerr << "Processing level: " << level << '\n';
    cerr << "List size: " << (triSk1.size()) << '\n';
    pos = mlVertex[level];
    elm = mlElement[level];
    counter = 0;
    int chIdx = 0;
    for (q = triSk1.begin(); q != triSk1.end(); ++q)
    {
      TriSkeleton* skl = *q;

      TriSkeleton* ch0 = new TriSkeleton;
      TriSkeleton* ch1 = new TriSkeleton;
      TriSkeleton* ch2 = new TriSkeleton;
      TriSkeleton* ch3 = new TriSkeleton;

      int chLoc = (chIdx%4);
      if (chLoc == 3) chLoc = 0;
      int chCh = 0;
      for (int n = 0; n < 3; n++)
      {
        ch0->pos[0][n] = pos[elm[counter][0]][n];
        ch0->pos[1][n] = pos[elm[counter][1]][n];
        ch0->pos[2][n] = pos[elm[counter][2]][n];
        chCh++;
        ch1->pos[0][n] = pos[elm[counter + 1][0]][n];
        ch1->pos[1][n] = pos[elm[counter + 1][1]][n];
        ch1->pos[2][n] = pos[elm[counter + 1][2]][n];
        chCh++;
        ch2->pos[0][n] = pos[elm[counter + 2][0]][n];
        ch2->pos[1][n] = pos[elm[counter + 2][1]][n];
        ch2->pos[2][n] = pos[elm[counter + 2][2]][n];
        chCh++;
        ch3->pos[0][n] = pos[elm[counter + 3][0]][n];
        ch3->pos[1][n] = pos[elm[counter + 3][1]][n];
        ch3->pos[2][n] = pos[elm[counter + 3][2]][n];
      }

      // set children
      skl->child[0] = ch0;
      skl->child[1] = ch1;
      skl->child[2] = ch2;
      skl->child[3] = ch3;

      triSk2.push_back(ch0);
      triSk2.push_back(ch1);
      triSk2.push_back(ch2);
      triSk2.push_back(ch3);
      counter += 4;
      chIdx++;
    }

    // copy lists
    triSk1.clear();
    for (q = triSk2.begin(); q != triSk2.end(); ++q)
    {
      triSk1.push_back(*q);
    }
    triSk2.clear();
  }
}


} // namespace